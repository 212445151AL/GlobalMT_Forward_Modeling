// file: fem.h
// ============================================================
// Finite Element Method interface for motional induction
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#ifndef FINITE_ELEMENT_SOLVER_HEADER
#define FINITE_ELEMENT_SOLVER_HEADER

#include "config.h"
#include <mpi.h>
#include "parameterHandler.h"
#include "MTCurrentSource.h"
#include "solver.h"
#include "mfem.hpp"
#include "PwMatrixCoefficient.h"
#include <fstream>
#include "errorEstimators.h"
#include <vector>

// ============================================================================
// FiniteElementSolver - Main solver class for motional induction problems
// ============================================================================
class FiniteElementSolver
{
public:
    // ------------------------------------------------------------------------
    // Construction / Destruction
    // ------------------------------------------------------------------------
    FiniteElementSolver(
        ParameterHandler& input_params,
        Source& current_source,
        mfem::ParMesh& mesh_instance
    );
    
    virtual ~FiniteElementSolver();

    // ------------------------------------------------------------------------
    // Initialization
    // ------------------------------------------------------------------------
    void initialize(int polynomial_order = 1);
    void setup_coefficients();

    // ------------------------------------------------------------------------
    // Utility methods
    // ------------------------------------------------------------------------
    void print_elapsed_time(double seconds, const std::string& message = "MPI_Wtime: ");
    HYPRE_Int get_total_degrees_of_freedom();
    void print_problem_statistics();

    // ------------------------------------------------------------------------
    // Assembly methods
    // ------------------------------------------------------------------------
    void assemble_source_term();
    void form_linear_system();
    void setup_preconditioner(SolverParameters& solver_params);
    void solve();

    // ------------------------------------------------------------------------
    // Station location methods
    // ------------------------------------------------------------------------
    void locate_station_elements();
    void locate_station_elements_with_gslib();

    // ------------------------------------------------------------------------
    // Assembly helpers
    // ------------------------------------------------------------------------
    void assemble_system_matrix();
    void assemble_right_hand_side();
    void assemble_dual_source_term();

    // ------------------------------------------------------------------------
    // Error estimation and adaptation
    // ------------------------------------------------------------------------
    void estimate_error();
    void update();
    void run_forward_modeling(int source_identifier);

    // ------------------------------------------------------------------------
    // Output methods
    // ------------------------------------------------------------------------
    void write_vtk_file(std::ofstream& output_stream);

    // ------------------------------------------------------------------------
    // Public member variables
    // ------------------------------------------------------------------------

    // MPI variables
    int rank_id;
    int num_procs;
    MPI_Comm mpi_communicator;

    // Class ParamAdministrator for handlering global MT input parameters
    ParameterHandler* parameters;

    // Motional induction current source
    Source* motional_current;

    // Parallel mesh
    mfem::ParMesh* parallel_mesh;

    // Parallel H(curl) finite-element space
    mfem::ParFiniteElementSpace* fe_space;

    // Ccefficients of governing curl-curl equation
    mfem::Coefficient* stiffness_coefficient;
    mfem::Coefficient* mass_coefficient;
    PiecewiseMatrixCoefficient* tensor_mass_coefficient;
    mfem::Coefficient* stiffness_negative;

    // Weak formulation B(U,V) = D(V)
    mfem::ParSesquilinearForm* bilinear_form;
    mfem::ParComplexLinearForm* linear_form;
    mfem::ParComplexGridFunction* solution_field;

    // Weak formulation B(W,V) = L(V) for dual problem 
    mfem::ParComplexLinearForm* dual_linear_form;
    mfem::ParComplexGridFunction* dual_solution_field;

    // Linear system KU = F
    mfem::OperatorHandle system_matrix;
    mfem::Vector primal_rhs_vector, primal_solution_vector;
    mfem::Vector dual_rhs_vector, dual_solution_vector;

    // Auxiliary-space block-diagonal preconditioner
    mfem::ParBilinearForm* preconditioner;
    mfem::HypreParMatrix prec_matrix;

    // Dirichlet boundary condition on the outer boundary
    mfem::Array<int> dirichlet_marker;
    mfem::Array<int> dirichlet_dof_list;

    // Solver control and error estimation
    bool solve_dual_problem;
    mfem::Vector final_error_indicators;
    mfem::Vector primal_error_indicators;

    mfem::Array<int> local_station_elements;
    mfem::Array<double> local_station_radii;
    mfem::Array<double> local_station_thetas;
    mfem::Array<double> local_station_phis;
};

#endif // FINITE_ELEMENT_SOLVER_HEADER
