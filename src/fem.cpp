// file: fem.cpp
// ============================================================
// Finite Element Method implementation for motional induction
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#include "fem.h"
#include "MTUtils.h"
#include <stdexcept>
#include <algorithm>
#include "postProcessor.h"
#include "PwMatrixCoefficient.h"

using namespace mfem;

// ============================================================================
// Construction / Destruction
// ============================================================================
FiniteElementSolver::FiniteElementSolver(
    ParameterHandler& input_params,
    Source& current_source,
    ParMesh& mesh_instance
)
    : rank_id(0)
    , num_procs(1)
    , parameters(&input_params)
    , motional_current(&current_source)
    , parallel_mesh(&mesh_instance)
    , fe_space(nullptr)
    , stiffness_coefficient(nullptr)
    , mass_coefficient(nullptr)
    , tensor_mass_coefficient(nullptr)
    , stiffness_negative(nullptr)
    , bilinear_form(nullptr)
    , linear_form(nullptr)
    , solution_field(nullptr)
    , dual_linear_form(nullptr)
    , dual_solution_field(nullptr)
    , preconditioner(nullptr)
{
    if (parameters->amr_method == 0)
    {
        solve_dual_problem = true;   // Goal-oriented AMR
    }
    else
    {
        solve_dual_problem = false;  // Non-goal-oriented or global AMR
    }
}

FiniteElementSolver::~FiniteElementSolver()
{
    delete fe_space;
    delete bilinear_form;
    delete linear_form;
    delete solution_field;
    delete dual_linear_form;
    delete dual_solution_field;
    delete preconditioner;
    delete stiffness_coefficient;
    delete mass_coefficient;
    delete tensor_mass_coefficient;
    delete stiffness_negative;
}

// ============================================================================
// Initialization methods
// ============================================================================
void FiniteElementSolver::setup_coefficients()
{
    stiffness_coefficient = new ConstantCoefficient(1.0);
    stiffness_negative = new ConstantCoefficient(-1.0);

    int attribute_count = parallel_mesh->attributes.Size();
    int spatial_dim = parallel_mesh->SpaceDimension();
    
    tensor_mass_coefficient = new PiecewiseMatrixCoefficient(spatial_dim, spatial_dim);
    
    for (int idx = 0; idx < attribute_count; ++idx)
    {
        int attribute = parallel_mesh->attributes[idx];
        double angular_freq = motional_current->angular_frequency;
        
        DenseMatrix conductivity_tensor(spatial_dim, spatial_dim);
        parameters->get_element_conductivity_tensor(attribute, conductivity_tensor);
        
        conductivity_tensor *= (angular_freq * MT::kMuZero);
        
        auto* coefficient = new MatrixConstantCoefficient(conductivity_tensor);
        tensor_mass_coefficient->add_coefficient(attribute, coefficient);
    }
}

void FiniteElementSolver::initialize(int polynomial_order)
{
    // MPI setup
    mpi_communicator = parallel_mesh->GetComm();
    MPI_Comm_size(mpi_communicator, &num_procs);
    MPI_Comm_rank(mpi_communicator, &rank_id);
    
    // Create H(curl) finite element space
    int spatial_dim = parallel_mesh->SpaceDimension();
    fe_space = new ParFiniteElementSpace(
        parallel_mesh,
        new ND_FECollection(polynomial_order, spatial_dim)
    );
    
    setup_coefficients();

    // Bilinear form setup
    bilinear_form = new ParSesquilinearForm(
        fe_space,
        ComplexOperator::BLOCK_SYMMETRIC
    );
    bilinear_form->AddDomainIntegrator(
        new CurlCurlIntegrator(*stiffness_coefficient),
        nullptr
    );
    bilinear_form->AddDomainIntegrator(
        nullptr,
        new VectorFEMassIntegrator(*tensor_mass_coefficient)
    );
    
    // Linear forms
    linear_form = new ParComplexLinearForm(
        fe_space,
        ComplexOperator::BLOCK_SYMMETRIC
    );
    linear_form->Vector::operator=(0.0);
    
    solution_field = new ParComplexGridFunction(fe_space);
    *solution_field = 0.0;

    // Dual problem components
    dual_linear_form = new ParComplexLinearForm(
        fe_space,
        ComplexOperator::BLOCK_SYMMETRIC
    );
    dual_linear_form->Vector::operator=(0.0);
    
    dual_solution_field = new ParComplexGridFunction(fe_space);
    *dual_solution_field = 0.0;

    // Preconditioner setup
    preconditioner = new ParBilinearForm(fe_space);
    preconditioner->AddDomainIntegrator(
        new CurlCurlIntegrator(*stiffness_coefficient)
    );
    preconditioner->AddDomainIntegrator(
        new VectorFEMassIntegrator(*tensor_mass_coefficient)
    );

    // Boundary condition handling
    int boundary_attrs = parallel_mesh->bdr_attributes.Size();
    if (boundary_attrs < 1 && rank_id == 0)
    {
        throw std::invalid_argument(
            "FiniteElementSolver::initialize(): Mesh must have at least one boundary attribute!"
        );
    }
    
    int dirichlet_location = -1;
    for (int idx = 0; idx < boundary_attrs; ++idx)
    {
        if (parallel_mesh->bdr_attributes[idx] == parameters->boundary_marker)
        {
            dirichlet_location = idx;
            break;
        }
    }
    
    if (dirichlet_location == -1 && rank_id == 0)
    {
        throw std::invalid_argument(
            "FiniteElementSolver::initialize(): Cannot find Dirichlet boundary attribute!"
        );
    }
    
    dirichlet_marker.SetSize(boundary_attrs);
    dirichlet_marker = 0;
    dirichlet_marker[dirichlet_location] = 1;
    
    fe_space->GetEssentialTrueDofs(dirichlet_marker, dirichlet_dof_list);
}

// ============================================================================
// Utility methods
// ============================================================================
void FiniteElementSolver::print_elapsed_time(
    double seconds,
    const std::string& message
)
{
    std::cout << message << seconds << " (s)\n";
}

HYPRE_Int FiniteElementSolver::get_total_degrees_of_freedom()
{
    return 2 * fe_space->GlobalTrueVSize();
}

void FiniteElementSolver::print_problem_statistics()
{
    HYPRE_Int global_size = fe_space->GlobalTrueVSize();
    long global_elements = parallel_mesh->GetGlobalNE();
    
    if (rank_id == 0)
    {
        std::cout << "Number of tetrahedral elements: " << global_elements << "\n";
        std::cout << "Number of complex-valued unknowns: " << global_size << "\n";
        std::cout << "Total Degrees of Freedom (DOFs): " << 2 * global_size << "\n";
    }
}

// ============================================================================
// Assembly methods
// ============================================================================
void FiniteElementSolver::assemble_system_matrix()
{
    bilinear_form->Assemble();
    bilinear_form->Finalize();
}

void FiniteElementSolver::assemble_right_hand_side()
{
    linear_form->Vector::operator=(0.0);
    assemble_source_term();
    
    if (solve_dual_problem)
    {
        dual_linear_form->Vector::operator=(0.0);
        assemble_dual_source_term();
    }
}

void FiniteElementSolver::assemble_source_term()
{
    int active_tet_count = motional_current->element_count;
    double angular_freq = motional_current->angular_frequency;
    
    for (int elem_id = 0; elem_id < parallel_mesh->GetNE(); ++elem_id)
    {
        int global_id = parallel_mesh->GetAttribute(elem_id);
        auto iter = motional_current->element_id_to_index.find(global_id);
        
        if (iter == motional_current->element_id_to_index.end())
        {
            continue;
        }
        
        int source_idx = iter->second;
        
        double current_real[3] = {
            motional_current->jx_real[source_idx],
            motional_current->jy_real[source_idx],
            motional_current->jz_real[source_idx]
        };
        
        double current_imag[3] = {
            motional_current->jx_imag[source_idx],
            motional_current->jy_imag[source_idx],
            motional_current->jz_imag[source_idx]
        };

        const auto* element = fe_space->GetFE(elem_id);
        auto* transform = fe_space->GetElementTransformation(elem_id);
        auto* rule = &IntRules.Get(
            element->GetGeomType(),
            3 * element->GetOrder() + transform->OrderW()
        );
        
        int num_dofs = element->GetDof();
        int dimension = element->GetDim();
        
        DenseMatrix shape_matrix(num_dofs, dimension);
        Vector real_contrib(num_dofs), imag_contrib(num_dofs);
        real_contrib = 0.0;
        imag_contrib = 0.0;
        
        for (int quad_idx = 0; quad_idx < rule->GetNPoints(); ++quad_idx)
        {
            const auto& point = rule->IntPoint(quad_idx);
            transform->SetIntPoint(&point);
            
            double weight_factor = transform->Weight() * point.weight;
            element->CalcPhysVShape(*transform, shape_matrix);
            
            double scale_factor = weight_factor * angular_freq * MT::kMuZero;
            
            for (int dof_idx = 0; dof_idx < num_dofs; ++dof_idx)
            {
                double shape_product = shape_matrix(dof_idx, 0) * current_imag[0] +
                                       shape_matrix(dof_idx, 1) * current_imag[1] +
                                       shape_matrix(dof_idx, 2) * current_imag[2];
                
                real_contrib[dof_idx] += scale_factor * shape_product;
                
                shape_product = shape_matrix(dof_idx, 0) * current_real[0] +
                                shape_matrix(dof_idx, 1) * current_real[1] +
                                shape_matrix(dof_idx, 2) * current_real[2];
                
                imag_contrib[dof_idx] += -scale_factor * shape_product;
            }
        }
        
        Array<int> dof_indices;
        fe_space->GetElementVDofs(elem_id, dof_indices);
        linear_form->real().AddElementVector(dof_indices, real_contrib);
        linear_form->imag().AddElementVector(dof_indices, imag_contrib);
    }
}

void FiniteElementSolver::assemble_dual_source_term()
{
    for (int idx = 0; idx < local_station_elements.Size(); ++idx)
    {
        int elem_id = local_station_elements[idx];
        
        const auto* element = fe_space->GetFE(elem_id);
        auto* rule = &IntRules.Get(
            element->GetGeomType(),
            2 * element->GetOrder()
        );
        auto* transform = fe_space->GetElementTransformation(elem_id);

        double volume = parallel_mesh->GetElementVolume(elem_id);
        double inv_volume = 1.0 / volume;
        
        int num_dofs = element->GetDof();
        int dimension = element->GetDim();
        
        DenseMatrix shape_matrix(num_dofs, dimension);
        Array<int> dof_indices;
        fe_space->GetElementVDofs(elem_id, dof_indices);
        
        Vector element_vector(num_dofs);
        element_vector = 0.0;
        
        for (int quad_idx = 0; quad_idx < rule->GetNPoints(); ++quad_idx)
        {
            const auto& point = rule->IntPoint(quad_idx);
            transform->SetIntPoint(&point);
            
            element->CalcPhysVShape(*transform, shape_matrix);
            double weight_factor = transform->Weight() * point.weight;
            
            for (int dof_idx = 0; dof_idx < num_dofs; ++dof_idx)
            {
                element_vector[dof_idx] += weight_factor * inv_volume * (
                    shape_matrix(dof_idx, 0) +
                    shape_matrix(dof_idx, 1) +
                    shape_matrix(dof_idx, 2)
                );
            }
        }
        
        dual_linear_form->real().AddElementVector(dof_indices, element_vector);
    }
}

// ============================================================================
// Station location methods
// ============================================================================
void FiniteElementSolver::locate_station_elements()
{
    local_station_elements.DeleteAll();
    local_station_radii.DeleteAll();
    local_station_thetas.DeleteAll();
    local_station_phis.DeleteAll();

    std::ifstream site_stream(parameters->sites_file);
    if (!site_stream)
    {
        throw std::runtime_error(
            "FiniteElementSolver::locate_station_elements(): failed to open sites file: "
            + parameters->sites_file
        );
    }
    
    int total_sites;
    if (!(site_stream >> total_sites) || total_sites <= 0)
    {
        throw std::runtime_error(
            "FiniteElementSolver::locate_station_elements(): invalid station count in sites file: "
            + parameters->sites_file
        );
    }
    
    Array<double> global_radii(total_sites), global_thetas(total_sites), global_phis(total_sites);
    Array<double> global_x(total_sites), global_y(total_sites), global_z(total_sites);
    
    for (int idx = 0; idx < total_sites; ++idx)
    {
        double r_val, theta_val, phi_val;
        if (!(site_stream >> r_val >> theta_val >> phi_val))
        {
            throw std::runtime_error(
                "FiniteElementSolver::locate_station_elements(): malformed station record in "
                + parameters->sites_file
            );
        }
        
        if (r_val < 0.0 || theta_val < 0.0 || theta_val > 180.0 ||
            phi_val < 0.0 || phi_val > 360.0)
        {
            throw std::invalid_argument(
                "FiniteElementSolver::locate_station_elements(): station coordinates are out of range."
            );
        }
        
        global_radii[idx] = r_val;
        global_thetas[idx] = theta_val;
        global_phis[idx] = phi_val;
        
        double theta_rad = theta_val / 180.0 * MT::kPi;
        double phi_rad = phi_val / 180.0 * MT::kPi;
        
        global_x[idx] = r_val * sin(theta_rad) * cos(phi_rad);
        global_y[idx] = r_val * sin(theta_rad) * sin(phi_rad);
        global_z[idx] = r_val * cos(theta_rad);
    }
    
    DenseMatrix point_matrix(3, total_sites);
    Array<int> element_indices;
    Array<IntegrationPoint> integration_points;
    
    for (int idx = 0; idx < total_sites; ++idx)
    {
        point_matrix(0, idx) = global_x[idx];
        point_matrix(1, idx) = global_y[idx];
        point_matrix(2, idx) = global_z[idx];
    }
    
    parallel_mesh->FindPoints(point_matrix, element_indices, integration_points);
    
    for (int idx = 0; idx < element_indices.Size(); ++idx)
    {
        int elem_id = element_indices[idx];
        if (elem_id == -1 && rank_id == 0)
        {
            throw std::invalid_argument(
                "FiniteElementSolver::locate_station_elements(): Some stations not found!"
            );
        }
        if (elem_id != -2)
        {
            local_station_elements.Append(elem_id);
            local_station_radii.Append(global_radii[idx]);
            local_station_thetas.Append(global_thetas[idx]);
            local_station_phis.Append(global_phis[idx]);
        }
    }
}

void FiniteElementSolver::locate_station_elements_with_gslib()
{
#ifdef ENABLE_GSLIB
    local_station_elements.DeleteAll();
    local_station_radii.DeleteAll();
    local_station_thetas.DeleteAll();
    local_station_phis.DeleteAll();

    std::ifstream site_stream(parameters->sites_file);
    if (!site_stream)
    {
        throw std::runtime_error(
            "FiniteElementSolver::locate_station_elements_with_gslib(): failed to open sites file: "
            + parameters->sites_file
        );
    }
    
    int total_sites;
    if (!(site_stream >> total_sites) || total_sites <= 0)
    {
        throw std::runtime_error(
            "FiniteElementSolver::locate_station_elements_with_gslib(): invalid station count in sites file: "
            + parameters->sites_file
        );
    }
    
    std::vector<double> global_radii(total_sites), global_thetas(total_sites), global_phis(total_sites);
    std::vector<double> global_x(total_sites), global_y(total_sites), global_z(total_sites);
    
    for (int idx = 0; idx < total_sites; ++idx)
    {
        double r_val, theta_val, phi_val;
        if (!(site_stream >> r_val >> theta_val >> phi_val))
        {
            throw std::runtime_error(
                "FiniteElementSolver::locate_station_elements_with_gslib(): malformed station record in "
                + parameters->sites_file
            );
        }
        
        if (r_val < 0.0 || theta_val < 0.0 || theta_val > 180.0 ||
            phi_val < 0.0 || phi_val > 360.0)
        {
            throw std::invalid_argument(
                "FiniteElementSolver::locate_station_elements_with_gslib(): station coordinates are out of range."
            );
        }
        
        global_radii[idx] = r_val;
        global_thetas[idx] = theta_val;
        global_phis[idx] = phi_val;
        
        double theta_rad = theta_val / 180.0 * MT::kPi;
        double phi_rad = phi_val / 180.0 * MT::kPi;
        
        global_x[idx] = r_val * sin(theta_rad) * cos(phi_rad);
        global_y[idx] = r_val * sin(theta_rad) * sin(phi_rad);
        global_z[idx] = r_val * cos(theta_rad);
    }

    Vector coordinates(total_sites * 3);
    for (int idx = 0; idx < total_sites; ++idx)
    {
        coordinates(idx) = global_x[idx];
        coordinates(total_sites + idx) = global_y[idx];
        coordinates(2 * total_sites + idx) = global_z[idx];
    }

    FindPointsGSLIB point_locator(parallel_mesh->GetComm());
    point_locator.Setup(*parallel_mesh);
    point_locator.FindPoints(coordinates);
    
    Array<unsigned int> status_codes = point_locator.GetCode();
    Array<unsigned int> processor_ids = point_locator.GetProc();
    Array<unsigned int> element_ids = point_locator.GetElem();

    int missing_flag = status_codes.Find(2);
    if (missing_flag != -1)
    {
        throw std::invalid_argument(
            "FiniteElementSolver::locate_station_elements_with_gslib(): Some stations not found!"
        );
    }

    for (int idx = 0; idx < total_sites; ++idx)
    {
        if (processor_ids[idx] == static_cast<unsigned int>(rank_id))
        {
            local_station_elements.Append(element_ids[idx]);
            local_station_radii.Append(global_radii[idx]);
            local_station_thetas.Append(global_thetas[idx]);
            local_station_phis.Append(global_phis[idx]);
        }
    }

    point_locator.FreeData();
#endif
}

// ============================================================================
// Linear system setup
// ============================================================================
void FiniteElementSolver::form_linear_system()
{
    assemble_system_matrix();
    assemble_right_hand_side();
    
    double start_time, end_time;
    if (rank_id == 0)
    {
        std::cout << "\nStep 1: Forming linear system...\n";
    }
    
    Vector zero_vector(parallel_mesh->SpaceDimension());
    zero_vector = 0.0;
    VectorConstantCoefficient zero_coefficient(zero_vector);
    
    solution_field->ProjectBdrCoefficientTangent(
        zero_coefficient,
        zero_coefficient,
        dirichlet_marker
    );

    start_time = MPI_Wtime();
    bilinear_form->FormLinearSystem(
        dirichlet_dof_list,
        *solution_field,
        *linear_form,
        system_matrix,
        primal_solution_vector,
        primal_rhs_vector
    );

    if (solve_dual_problem)
    {
        dual_solution_field->ProjectBdrCoefficientTangent(
            zero_coefficient,
            zero_coefficient,
            dirichlet_marker
        );
        
        bilinear_form->FormLinearSystem(
            dirichlet_dof_list,
            *dual_solution_field,
            *dual_linear_form,
            system_matrix,
            dual_solution_vector,
            dual_rhs_vector
        );
    }

    end_time = MPI_Wtime();
    if (rank_id == 0)
    {
        print_elapsed_time(end_time - start_time, "form_linear_system: ");
    }
}

void FiniteElementSolver::setup_preconditioner(SolverParameters& solver_params)
{
    double start_time, end_time;
    if (rank_id == 0)
    {
        std::cout << "\nStep 2: Setting up preconditioner...\n";
        if (solver_params.prec_type == 0)
            std::cout << "Preconditioner: AMS\n";
        else if (solver_params.prec_type == 1)
            std::cout << "Preconditioner: PCG-AMS\n";
        else
            std::cout << "Preconditioner: Multigrid\n";
    }
    
    start_time = MPI_Wtime();
    preconditioner->Assemble();
    preconditioner->Finalize();
    end_time = MPI_Wtime();
    
    if (rank_id == 0)
    {
        print_elapsed_time(end_time - start_time, "Preconditioner assembly: ");
    }

    start_time = MPI_Wtime();
    preconditioner->FormSystemMatrix(dirichlet_dof_list, prec_matrix);
    end_time = MPI_Wtime();
    
    if (rank_id == 0)
    {
        print_elapsed_time(end_time - start_time, "FormSystemMatrix: ");
    }
}

void FiniteElementSolver::solve()
{
    SolverParameters solver_params(parameters->linear_opts_file);
    setup_preconditioner(solver_params);
    
    Array<int> block_offsets;
    block_offsets.SetSize(3);
    block_offsets[0] = 0;
    block_offsets[1] = prec_matrix.Height();
    block_offsets[2] = prec_matrix.Height();
    block_offsets.PartialSum();
    
    BlockDiagonalPreconditioner block_preconditioner(block_offsets);

    if (solver_params.prec_type == 0)
    {
        auto* ams_real = new HypreAMS(prec_matrix, fe_space);
        solver_params.apply_to_ams(*ams_real);
        auto* ams_imag = new HypreAMS(prec_matrix, fe_space);
        solver_params.apply_to_ams(*ams_imag);
        
        block_preconditioner.SetDiagonalBlock(0, ams_real);
        block_preconditioner.SetDiagonalBlock(1, ams_imag);
    }
    else if (solver_params.prec_type == 1)
    {
        auto* ams = new HypreAMS(prec_matrix, fe_space);
        auto* pcg_real = new HyprePCG(prec_matrix);
        pcg_real->SetPreconditioner(*ams);
        solver_params.apply_to_pcg(*pcg_real);
        
        auto* pcg_imag = new HyprePCG(prec_matrix);
        pcg_imag->SetPreconditioner(*ams);
        solver_params.apply_to_pcg(*pcg_imag);
        
        block_preconditioner.SetDiagonalBlock(0, pcg_real);
        block_preconditioner.SetDiagonalBlock(1, pcg_imag);
    }
    else if (solver_params.prec_type == 2)
    {
        OperatorHandle prec_operator;
        preconditioner->FormSystemMatrix(dirichlet_dof_list, prec_operator);
        
        auto* multigrid_real = new HypreAMS(
            *prec_operator.As<HypreParMatrix>(),
            fe_space
        );
        
        auto* multigrid_imag = new ScaledOperator(
            multigrid_real,
            (ComplexOperator::BLOCK_SYMMETRIC == ComplexOperator::HERMITIAN) ? -1.0 : 1.0
        );
        
        block_preconditioner.SetDiagonalBlock(0, multigrid_real);
        block_preconditioner.SetDiagonalBlock(1, multigrid_imag);
    }
    
    block_preconditioner.owns_blocks = 1;

    FlexibleGMRES fgmres_solver(mpi_communicator);
    fgmres_solver.SetOperator(*system_matrix.Ptr());
    fgmres_solver.SetPreconditioner(block_preconditioner);
    fgmres_solver.SetMaxIter(solver_params.fgmres_max_iterations);
    fgmres_solver.SetRelTol(solver_params.fgmres_primal_tolerance);
    fgmres_solver.SetKDim(solver_params.fgmres_restart_dimension);
    fgmres_solver.SetPrintLevel(solver_params.fgmres_output_level);
    
    double start_time, end_time;
    if (rank_id == 0)
    {
        std::cout << "\nStep 3: Solving linear system...\n";
    }
    
    start_time = MPI_Wtime();
    fgmres_solver.Mult(primal_rhs_vector, primal_solution_vector);
    end_time = MPI_Wtime();
    
    if (rank_id == 0)
    {
        print_elapsed_time(end_time - start_time);
    }
    
    bilinear_form->RecoverFEMSolution(
        primal_solution_vector,
        *linear_form,
        *solution_field
    );

    if (solve_dual_problem)
    {
        fgmres_solver.SetRelTol(solver_params.fgmres_dual_tolerance);
        
        if (rank_id == 0)
        {
            std::cout << "Solving dual problem...\n";
        }
        
        start_time = MPI_Wtime();
        fgmres_solver.Mult(dual_rhs_vector, dual_solution_vector);
        end_time = MPI_Wtime();
        
        if (rank_id == 0)
        {
            std::cout << "Dual solve time: " << end_time - start_time << " (s)\n";
        }
        
        bilinear_form->RecoverFEMSolution(
            dual_solution_vector,
            *dual_linear_form,
            *dual_solution_field
        );
    }
}

// ============================================================================
// Error estimation and adaptation
// ============================================================================
void FiniteElementSolver::estimate_error()
{
    FaceJumpEstimator estimator;
    
    primal_error_indicators = estimator.get_error_estimate(
        *parameters,
        *solution_field
    );
    
    int num_elements = primal_error_indicators.Size();
    final_error_indicators.SetSize(num_elements);

    if (parameters->amr_method == 0)
    {
        Vector dual_error = estimator.get_error_estimate(
            *parameters,
            *dual_solution_field
        );
        
        for (int elem_idx = 0; elem_idx < num_elements; ++elem_idx)
        {
            final_error_indicators[elem_idx] =
                primal_error_indicators[elem_idx] * dual_error[elem_idx];
        }
    }
    else
    {
        for (int elem_idx = 0; elem_idx < num_elements; ++elem_idx)
        {
            final_error_indicators[elem_idx] = primal_error_indicators[elem_idx];
        }
    }
}

void FiniteElementSolver::update()
{
    fe_space->Update();

    solution_field->Update();
    *solution_field = 0.0;
    
    dual_solution_field->Update();
    *dual_solution_field = 0.0;

    bilinear_form->Update();
    linear_form->Update();
    dual_linear_form->Update();
    preconditioner->Update();
}

// ============================================================================
// Main forward modeling loop
// ============================================================================
void FiniteElementSolver::run_forward_modeling(int source_identifier)
{
    double total_start = MPI_Wtime();
    
    int iteration = 0;
    bool max_dofs_reached = false;

    initialize();
    
    for (iteration = 0; iteration < parameters->max_refinement_iterations; ++iteration)
    {
        if (rank_id == 0)
        {
            std::cout << "-------Goal-oriented hp-adaptive refinement loop ------- #"
                      << iteration + 1 << "\n";
        }
        
        int mesh_sequence = parallel_mesh->GetSequence();
        
        if (rank_id == 0)
        {
            std::cout << "Polynomial order: 1\n";
        }

        print_problem_statistics();

        if (rank_id == 0)
        {
            std::cout << "Locating station elements...\n";
        }
        
        double start_time = MPI_Wtime();
        
#ifdef ENABLE_GSLIB
        locate_station_elements_with_gslib();
#else
        locate_station_elements();
#endif
        
        double end_time = MPI_Wtime();
        
        if (rank_id == 0)
        {
            std::cout << "Location time: " << end_time - start_time << " (s)\n";
        }

        if (rank_id == 0)
        {
            std::cout << "Forming linear system and preconditioner...\n";
        }
        
        start_time = MPI_Wtime();
        form_linear_system();
        end_time = MPI_Wtime();
        
        if (rank_id == 0)
        {
            std::cout << "Assembly time: " << end_time - start_time << " (s)\n";
        }

        solve();

        if (rank_id == 0)
        {
            std::cout << "\n############### Post-processing ###############\n";
        }
        
        double post_start = MPI_Wtime();
        int poly_order = parameters->polynomial_order;
        
        std::ostringstream solution_filename;
        solution_filename << "Forward_solutions/"
                          << "MT_source" << source_identifier
                          << "_iter" << iteration << ".sol";
        
        std::ofstream solution_file(solution_filename.str());
        solution_file.precision(8);

        PostProcessor post_processor(
            *parameters,
            *parallel_mesh,
            *fe_space,
            solution_field->real(),
            solution_field->imag(),
            motional_current->angular_frequency
        );
        
        post_processor.execute();
        post_processor.save_as_single_file(solution_file);

        estimate_error();

        if (parameters->print_vtk)
        {
            std::ostringstream vtk_filename;
            vtk_filename << "Forward_solutions/"
                         << "period_amr" << mesh_sequence << "_p1.vtk";
            
            std::ofstream vtk_file(vtk_filename.str());
            write_vtk_file(vtk_file);
        }

        if (parameters->amr_method == 0 || parameters->amr_method == 1)
        {
            double local_max = final_error_indicators.Max();
            double global_max;
            
            MPI_Allreduce(
                &local_max,
                &global_max,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                parallel_mesh->GetComm()
            );

            double threshold = parameters->refinement_threshold_ratio * global_max;
            parallel_mesh->RefineByError(final_error_indicators, threshold);
        }

        update();

        if (get_total_degrees_of_freedom() > parameters->max_dofs)
        {
            max_dofs_reached = true;
            
            if (rank_id == 0)
            {
                std::cout << "\nStopping: Maximum DOFs reached\n";
                std::cout << "Current DOFs: " << get_total_degrees_of_freedom() << "\n";
                std::cout << "Maximum DOFs: " << parameters->max_dofs << "\n";
            }
            break;
        }

        double total_end = MPI_Wtime();
        
        if (rank_id == 0)
        {
            std::cout << "Cumulative time: " << total_end - total_start << " (s)\n";
            std::cout << "\n --------------------- Next loop --------------------\n";
        }
    }
}

// ============================================================================
// VTK output
// ============================================================================
void FiniteElementSolver::write_vtk_file(std::ofstream& output_stream)
{
    if (rank_id == 0)
    {
        output_stream << "# vtk DataFile Version 3.0\n"
                      << "Generated by Liangyu Xie\n"
                      << "ASCII\n"
                      << "DATASET UNSTRUCTURED_GRID\n";
    }

    MPI_Comm communicator = parallel_mesh->GetComm();
    int comm_size;
    MPI_Comm_size(communicator, &comm_size);
    int spatial_dim = parallel_mesh->SpaceDimension();

    // ------------------------------------------------------------------------
    // Gather vertex coordinates
    // ------------------------------------------------------------------------
    int local_vertices = parallel_mesh->GetNV();
    std::vector<int> vertex_counts;
    if (rank_id == 0)
    {
        vertex_counts.resize(comm_size);
    }
    MPI_Gather(
        &local_vertices,
        1,
        MPI_INT,
        (rank_id == 0) ? vertex_counts.data() : nullptr,
        1,
        MPI_INT,
        0,
        communicator
    );

    std::vector<double> local_vertex_data(local_vertices * spatial_dim);
    int vertex_offset_local = 0;
    for (int idx = 0; idx < local_vertices; ++idx)
    {
        double* vertex = parallel_mesh->GetVertex(idx);
        for (int dim = 0; dim < spatial_dim; ++dim)
        {
            local_vertex_data[vertex_offset_local++] = vertex[dim];
        }
    }

    std::vector<int> vertex_value_counts;
    std::vector<int> vertex_value_displacements;
    std::vector<double> global_vertex_data;
    int global_vertices = 0;

    if (rank_id == 0)
    {
        vertex_value_counts.resize(comm_size);
        vertex_value_displacements.resize(comm_size);

        int displacement = 0;
        for (int proc = 0; proc < comm_size; ++proc)
        {
            vertex_value_counts[proc] = vertex_counts[proc] * spatial_dim;
            vertex_value_displacements[proc] = displacement;
            displacement += vertex_value_counts[proc];
            global_vertices += vertex_counts[proc];
        }

        global_vertex_data.resize(displacement);
    }

    MPI_Gatherv(
        local_vertex_data.data(),
        local_vertices * spatial_dim,
        MPI_DOUBLE,
        (rank_id == 0) ? global_vertex_data.data() : nullptr,
        (rank_id == 0) ? vertex_value_counts.data() : nullptr,
        (rank_id == 0) ? vertex_value_displacements.data() : nullptr,
        MPI_DOUBLE,
        0,
        communicator
    );

    // ------------------------------------------------------------------------
    // Gather element connectivity
    // ------------------------------------------------------------------------
    int local_elements = parallel_mesh->GetNE();
    int local_vertices_per_element =
        (local_elements > 0) ? parallel_mesh->GetElement(0)->GetNVertices() : 0;
    int vertices_per_element = 0;
    MPI_Allreduce(
        &local_vertices_per_element,
        &vertices_per_element,
        1,
        MPI_INT,
        MPI_MAX,
        communicator
    );

    std::vector<int> element_counts;
    if (rank_id == 0)
    {
        element_counts.resize(comm_size);
    }
    MPI_Gather(
        &local_elements,
        1,
        MPI_INT,
        (rank_id == 0) ? element_counts.data() : nullptr,
        1,
        MPI_INT,
        0,
        communicator
    );

    std::vector<int> local_element_vertices(local_elements * vertices_per_element);
    for (int idx = 0; idx < local_elements; ++idx)
    {
        Array<int> elem_vertices;
        parallel_mesh->GetElementVertices(idx, elem_vertices);
        for (int jdx = 0; jdx < vertices_per_element; ++jdx)
        {
            local_element_vertices[idx * vertices_per_element + jdx] = elem_vertices[jdx];
        }
    }

    std::vector<int> element_displacements;
    std::vector<int> element_value_counts;
    std::vector<int> element_value_displacements;
    std::vector<int> global_element_vertices;
    int global_elements = 0;

    if (rank_id == 0)
    {
        element_displacements.resize(comm_size);
        element_value_counts.resize(comm_size);
        element_value_displacements.resize(comm_size);

        int element_disp = 0;
        int value_disp = 0;
        for (int proc = 0; proc < comm_size; ++proc)
        {
            element_displacements[proc] = element_disp;
            element_value_displacements[proc] = value_disp;
            element_value_counts[proc] = element_counts[proc] * vertices_per_element;

            element_disp += element_counts[proc];
            value_disp += element_value_counts[proc];
            global_elements += element_counts[proc];
        }

        global_element_vertices.resize(value_disp);
    }

    MPI_Gatherv(
        local_element_vertices.data(),
        local_elements * vertices_per_element,
        MPI_INT,
        (rank_id == 0) ? global_element_vertices.data() : nullptr,
        (rank_id == 0) ? element_value_counts.data() : nullptr,
        (rank_id == 0) ? element_value_displacements.data() : nullptr,
        MPI_INT,
        0,
        communicator
    );

    // ------------------------------------------------------------------------
    // Gather per-element scalars and tensors
    // ------------------------------------------------------------------------
    double local_max = primal_error_indicators.Max();
    double global_max = 0.0;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, communicator);

    std::vector<double> local_error_data(local_elements);
    std::vector<int> local_marker_data(local_elements);
    std::vector<double> local_volume_data(local_elements);
    std::vector<double> local_tensor_data(local_elements * 9);

    for (int idx = 0; idx < local_elements; ++idx)
    {
        local_error_data[idx] = primal_error_indicators[idx];
        local_marker_data[idx] = parallel_mesh->GetAttribute(idx);
        local_volume_data[idx] = parallel_mesh->GetElementVolume(idx);

        DenseMatrix conductivity(3, 3);
        parameters->get_element_conductivity_tensor(
            parallel_mesh->GetAttribute(idx),
            conductivity
        );

        local_tensor_data[idx * 9 + 0] = conductivity(0, 0);
        local_tensor_data[idx * 9 + 1] = conductivity(0, 1);
        local_tensor_data[idx * 9 + 2] = conductivity(0, 2);
        local_tensor_data[idx * 9 + 3] = conductivity(1, 0);
        local_tensor_data[idx * 9 + 4] = conductivity(1, 1);
        local_tensor_data[idx * 9 + 5] = conductivity(1, 2);
        local_tensor_data[idx * 9 + 6] = conductivity(2, 0);
        local_tensor_data[idx * 9 + 7] = conductivity(2, 1);
        local_tensor_data[idx * 9 + 8] = conductivity(2, 2);
    }

    std::vector<double> global_error_data;
    std::vector<int> global_marker_data;
    std::vector<double> global_volume_data;
    std::vector<double> global_tensor_data;
    std::vector<int> tensor_value_counts;
    std::vector<int> tensor_value_displacements;

    if (rank_id == 0)
    {
        global_error_data.resize(global_elements);
        global_marker_data.resize(global_elements);
        global_volume_data.resize(global_elements);

        tensor_value_counts.resize(comm_size);
        tensor_value_displacements.resize(comm_size);

        int tensor_disp = 0;
        for (int proc = 0; proc < comm_size; ++proc)
        {
            tensor_value_counts[proc] = element_counts[proc] * 9;
            tensor_value_displacements[proc] = tensor_disp;
            tensor_disp += tensor_value_counts[proc];
        }
        global_tensor_data.resize(tensor_disp);
    }

    MPI_Gatherv(
        local_error_data.data(),
        local_elements,
        MPI_DOUBLE,
        (rank_id == 0) ? global_error_data.data() : nullptr,
        (rank_id == 0) ? element_counts.data() : nullptr,
        (rank_id == 0) ? element_displacements.data() : nullptr,
        MPI_DOUBLE,
        0,
        communicator
    );

    MPI_Gatherv(
        local_marker_data.data(),
        local_elements,
        MPI_INT,
        (rank_id == 0) ? global_marker_data.data() : nullptr,
        (rank_id == 0) ? element_counts.data() : nullptr,
        (rank_id == 0) ? element_displacements.data() : nullptr,
        MPI_INT,
        0,
        communicator
    );

    MPI_Gatherv(
        local_volume_data.data(),
        local_elements,
        MPI_DOUBLE,
        (rank_id == 0) ? global_volume_data.data() : nullptr,
        (rank_id == 0) ? element_counts.data() : nullptr,
        (rank_id == 0) ? element_displacements.data() : nullptr,
        MPI_DOUBLE,
        0,
        communicator
    );

    MPI_Gatherv(
        local_tensor_data.data(),
        local_elements * 9,
        MPI_DOUBLE,
        (rank_id == 0) ? global_tensor_data.data() : nullptr,
        (rank_id == 0) ? tensor_value_counts.data() : nullptr,
        (rank_id == 0) ? tensor_value_displacements.data() : nullptr,
        MPI_DOUBLE,
        0,
        communicator
    );

    // ------------------------------------------------------------------------
    // Root writes VTK content
    // ------------------------------------------------------------------------
    if (rank_id == 0)
    {
        output_stream << "POINTS " << global_vertices << " double\n";
        for (int idx = 0; idx < global_vertices; ++idx)
        {
            for (int dim = 0; dim < spatial_dim; ++dim)
            {
                output_stream << global_vertex_data[idx * spatial_dim + dim] << " ";
            }
            output_stream << "\n";
        }

        output_stream << "CELLS " << global_elements << " "
                      << global_elements * (vertices_per_element + 1) << "\n";

        std::vector<int> vertex_displacements(comm_size, 0);
        for (int proc = 1; proc < comm_size; ++proc)
        {
            vertex_displacements[proc] = vertex_displacements[proc - 1] + vertex_counts[proc - 1];
        }

        for (int proc = 0; proc < comm_size; ++proc)
        {
            int value_offset = element_value_displacements[proc];
            for (int elem_idx = 0; elem_idx < element_counts[proc]; ++elem_idx)
            {
                output_stream << vertices_per_element;
                for (int jdx = 0; jdx < vertices_per_element; ++jdx)
                {
                    int local_vertex_id =
                        global_element_vertices[value_offset + elem_idx * vertices_per_element + jdx];
                    output_stream << " " << vertex_displacements[proc] + local_vertex_id;
                }
                output_stream << "\n";
            }
        }

        output_stream << "CELL_TYPES " << global_elements << "\n";
        for (int idx = 0; idx < global_elements; ++idx)
        {
            output_stream << "10\n";  // VTK_TETRA = 10
        }

        output_stream << "CELL_DATA " << global_elements << "\n";
        output_stream << "SCALARS element_relative_error_indicator double 1\n"
                      << "LOOKUP_TABLE default\n";
        double error_scale = (global_max > 0.0) ? global_max : 1.0;
        for (int idx = 0; idx < global_elements; ++idx)
        {
            output_stream << global_error_data[idx] / error_scale << "\n";
        }

        output_stream << "SCALARS element_marker int 1\n"
                      << "LOOKUP_TABLE default\n";
        for (int idx = 0; idx < global_elements; ++idx)
        {
            output_stream << global_marker_data[idx] << "\n";
        }

        output_stream << "SCALARS element_process_id int 1\n"
                      << "LOOKUP_TABLE default\n";
        for (int proc = 0; proc < comm_size; ++proc)
        {
            for (int idx = 0; idx < element_counts[proc]; ++idx)
            {
                output_stream << proc << "\n";
            }
        }

        output_stream << "SCALARS element_volume double 1\n"
                      << "LOOKUP_TABLE default\n";
        for (int idx = 0; idx < global_elements; ++idx)
        {
            output_stream << global_volume_data[idx] << "\n";
        }

        output_stream << "TENSORS element_conductivity_tensor double\n";
        for (int idx = 0; idx < global_elements; ++idx)
        {
            int base = idx * 9;
            output_stream << global_tensor_data[base + 0] << " "
                          << global_tensor_data[base + 1] << " "
                          << global_tensor_data[base + 2] << "\n";
            output_stream << global_tensor_data[base + 3] << " "
                          << global_tensor_data[base + 4] << " "
                          << global_tensor_data[base + 5] << "\n";
            output_stream << global_tensor_data[base + 6] << " "
                          << global_tensor_data[base + 7] << " "
                          << global_tensor_data[base + 8] << "\n";
        }

        output_stream.close();
    }
}
