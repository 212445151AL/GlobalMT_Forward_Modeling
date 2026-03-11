// file: solver.cpp
// ============================================================
// Solver parameters and custom FGMRES implementation
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#include "solver.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>

using namespace mfem;

// ============================================================================
// SolverParameters - Constructor / Destructor
// ============================================================================
SolverParameters::SolverParameters(const std::string& parameter_filename)
{
    load_parameters(parameter_filename);
}

SolverParameters::~SolverParameters()
{
    // Nothing to clean up
}

// ============================================================================
// Parameter loading
// ============================================================================
void SolverParameters::load_parameters(const std::string& parameter_filename)
{
    std::ifstream input_stream(parameter_filename);
    std::string placeholder;

    // FGMRES parameters
    input_stream >> placeholder >> fgmres_max_iterations;
    input_stream >> placeholder >> fgmres_primal_tolerance;
    input_stream >> placeholder >> fgmres_dual_tolerance;
    input_stream >> placeholder >> fgmres_restart_dimension;
    input_stream >> placeholder >> fgmres_output_level;

    // Preconditioner type
    input_stream >> placeholder >> prec_type;
    
    if (prec_type != 0 && prec_type != 1 && prec_type != 2)
    {
        throw std::invalid_argument(
            "SolverParameters::load_parameters(): Invalid preconditioner type!"
        );
    }

    // PCG parameters
    input_stream >> placeholder >> pcg_max_iterations;
    input_stream >> placeholder >> pcg_tolerance;
    input_stream >> placeholder >> pcg_output_level;

    // AMS parameters
    input_stream >> placeholder >> ams_cycle_type;
    input_stream >> placeholder >> ams_max_iterations;
    input_stream >> placeholder >> ams_tolerance;
    input_stream >> placeholder >> ams_output_level;
    
    // Smoothing options
    input_stream >> placeholder >> a_relaxation_type;
    input_stream >> placeholder >> a_relaxation_sweeps;
    input_stream >> placeholder >> a_relaxation_weight;
    input_stream >> placeholder >> a_relaxation_omega;

    // AMG options
    input_stream >> placeholder >> amg_coarsening_type;
    input_stream >> placeholder >> amg_aggressive_levels;
    input_stream >> placeholder >> amg_strength_threshold;
    input_stream >> placeholder >> amg_interpolation_type;
    input_stream >> placeholder >> amg_max_elements_per_row;
    input_stream >> placeholder >> amg_relaxation_type;

    input_stream.close();
}

// ============================================================================
// Parameter output
// ============================================================================
void SolverParameters::write_to_file(const std::string& output_filename)
{
    std::ofstream output_stream(output_filename);

    output_stream << fgmres_max_iterations << "\n";
    output_stream << fgmres_primal_tolerance << "\n";
    output_stream << fgmres_dual_tolerance << "\n";
    output_stream << fgmres_restart_dimension << "\n";
    output_stream << fgmres_output_level << "\n\n";

    output_stream << prec_type << "\n\n";

    output_stream << pcg_max_iterations << "\n";
    output_stream << pcg_tolerance << "\n";
    output_stream << pcg_output_level << "\n\n";

    output_stream << ams_cycle_type << "\n";
    output_stream << ams_max_iterations << "\n";
    output_stream << ams_tolerance << "\n";
    output_stream << ams_output_level << "\n\n";
    
    output_stream << a_relaxation_type << "\n";
    output_stream << a_relaxation_sweeps << "\n";
    output_stream << a_relaxation_weight << "\n";
    output_stream << a_relaxation_omega << "\n\n";

    output_stream << amg_coarsening_type << "\n";
    output_stream << amg_aggressive_levels << "\n";
    output_stream << amg_strength_threshold << "\n";
    output_stream << amg_interpolation_type << "\n";
    output_stream << amg_max_elements_per_row << "\n";
    output_stream << amg_relaxation_type << "\n";

    output_stream.close();
}

// ============================================================================
// Parameter application
// ============================================================================
void SolverParameters::apply_to_pcg(HyprePCG& pcg_solver)
{
    pcg_solver.SetMaxIter(pcg_max_iterations);
    pcg_solver.SetTol(pcg_tolerance);
    pcg_solver.SetPrintLevel(pcg_output_level);
}

void SolverParameters::apply_to_ams(HypreAMS& ams_solver)
{
    HYPRE_AMSSetCycleType(ams_solver, ams_cycle_type);
    HYPRE_AMSSetMaxIter(ams_solver, ams_max_iterations);
    HYPRE_AMSSetTol(ams_solver, ams_tolerance);
    HYPRE_AMSSetPrintLevel(ams_solver, ams_output_level);
    
    HYPRE_AMSSetSmoothingOptions(
        ams_solver,
        a_relaxation_type,
        a_relaxation_sweeps,
        a_relaxation_weight,
        a_relaxation_omega
    );
    
    HYPRE_AMSSetAlphaAMGOptions(
        ams_solver,
        amg_coarsening_type,
        amg_aggressive_levels,
        amg_relaxation_type,
        amg_strength_threshold,
        amg_interpolation_type,
        amg_max_elements_per_row
    );
    
    HYPRE_AMSSetBetaAMGOptions(
        ams_solver,
        amg_coarsening_type,
        amg_aggressive_levels,
        amg_relaxation_type,
        amg_strength_threshold,
        amg_interpolation_type,
        amg_max_elements_per_row
    );

    HYPRE_AMSSetAlphaAMGCoarseRelaxType(ams_solver, amg_relaxation_type);
    HYPRE_AMSSetBetaAMGCoarseRelaxType(ams_solver, amg_relaxation_type);
}

// ============================================================================
// FlexibleGMRES - Constructors
// ============================================================================
FlexibleGMRES::FlexibleGMRES()
    : FGMRESSolver()
{
    // Nothing else
}

FlexibleGMRES::FlexibleGMRES(MPI_Comm communicator)
    : FGMRESSolver(communicator)
{
    // Nothing else
}

// ============================================================================
// FlexibleGMRES - Modified solver with residual output
// ============================================================================
void FlexibleGMRES::Mult(const Vector& rhs, Vector& solution) const
{
    int restart_dim = this->m;
    int max_iter = this->max_iter;
    double rel_tol = this->rel_tol;
    double abs_tol = this->abs_tol;
    int print_level = this->print_level;
    bool iterative_mode = this->iterative_mode;
    const Operator* oper = this->oper;
    const Solver* prec = this->prec;
    
    DenseMatrix hessenberg(restart_dim + 1, restart_dim);
    Vector s_vector(restart_dim + 1);
    Vector cs_vector(restart_dim + 1);
    Vector sn_vector(restart_dim + 1);
    
    Vector residual_vector(rhs.Size());

    int outer_iter, inner_iter, krylov_idx;

    if (iterative_mode)
    {
        oper->Mult(solution, residual_vector);
        subtract(rhs, residual_vector, residual_vector);
    }
    else
    {
        solution = 0.0;
        residual_vector = rhs;
    }
    
    double beta_norm = Norm(residual_vector);
    MFEM_ASSERT(IsFinite(beta_norm), "beta_norm = " << beta_norm);

    double final_norm = std::max(rel_tol * beta_norm, abs_tol);

    if (beta_norm <= final_norm)
    {
        this->final_norm = beta_norm;
        this->final_iter = 0;
        this->converged = 1;
        return;
    }
    
    double norm_rhs = Norm(rhs);
    
    if (print_level == 1)
    {
        mfem::out << "   Pass: " << std::setw(2) << 1
                  << "   Iteration: " << std::setw(3) << 0
                  << "   ||r|| = " << beta_norm << "\t"
                  << "||r||/||b|| = " << beta_norm / norm_rhs << std::endl;
    }

    this->Monitor(0, beta_norm, residual_vector, solution);

    Array<Vector*> krylov_vectors(restart_dim + 1);
    Array<Vector*> preconditioned_vectors(restart_dim + 1);
    
    for (int idx = 0; idx <= restart_dim; ++idx)
    {
        krylov_vectors[idx] = nullptr;
        preconditioned_vectors[idx] = nullptr;
    }

    inner_iter = 1;
    while (inner_iter <= max_iter)
    {
        if (krylov_vectors[0] == nullptr)
        {
            krylov_vectors[0] = new Vector(rhs.Size());
        }
        
        *krylov_vectors[0] = 0.0;
        krylov_vectors[0]->Add(1.0 / beta_norm, residual_vector);
        
        s_vector = 0.0;
        s_vector(0) = beta_norm;

        for (outer_iter = 0; outer_iter < restart_dim && inner_iter <= max_iter;
             ++outer_iter, ++inner_iter)
        {
            if (preconditioned_vectors[outer_iter] == nullptr)
            {
                preconditioned_vectors[outer_iter] = new Vector(rhs.Size());
            }
            
            *preconditioned_vectors[outer_iter] = 0.0;

            if (prec)
            {
                prec->Mult(*krylov_vectors[outer_iter],
                           *preconditioned_vectors[outer_iter]);
            }
            else
            {
                *preconditioned_vectors[outer_iter] = *krylov_vectors[outer_iter];
            }
            
            oper->Mult(*preconditioned_vectors[outer_iter], residual_vector);

            for (krylov_idx = 0; krylov_idx <= outer_iter; ++krylov_idx)
            {
                hessenberg(krylov_idx, outer_iter) = Dot(
                    residual_vector,
                    *krylov_vectors[krylov_idx]
                );
                
                residual_vector.Add(
                    -hessenberg(krylov_idx, outer_iter),
                    *krylov_vectors[krylov_idx]
                );
            }

            hessenberg(outer_iter + 1, outer_iter) = Norm(residual_vector);
            
            if (krylov_vectors[outer_iter + 1] == nullptr)
            {
                krylov_vectors[outer_iter + 1] = new Vector(rhs.Size());
            }
            
            *krylov_vectors[outer_iter + 1] = 0.0;
            krylov_vectors[outer_iter + 1]->Add(
                1.0 / hessenberg(outer_iter + 1, outer_iter),
                residual_vector
            );

            for (krylov_idx = 0; krylov_idx < outer_iter; ++krylov_idx)
            {
                apply_plane_rotation(
                    hessenberg(krylov_idx, outer_iter),
                    hessenberg(krylov_idx + 1, outer_iter),
                    cs_vector(krylov_idx),
                    sn_vector(krylov_idx)
                );
            }

            generate_plane_rotation(
                hessenberg(outer_iter, outer_iter),
                hessenberg(outer_iter + 1, outer_iter),
                cs_vector(outer_iter),
                sn_vector(outer_iter)
            );
            
            apply_plane_rotation(
                hessenberg(outer_iter, outer_iter),
                hessenberg(outer_iter + 1, outer_iter),
                cs_vector(outer_iter),
                sn_vector(outer_iter)
            );
            
            apply_plane_rotation(
                s_vector(outer_iter),
                s_vector(outer_iter + 1),
                cs_vector(outer_iter),
                sn_vector(outer_iter)
            );

            double residual_norm = fabs(s_vector(outer_iter + 1));
            MFEM_ASSERT(IsFinite(residual_norm), "residual_norm = " << residual_norm);

            if (print_level == 1)
            {
                mfem::out << "   Pass: " << std::setw(2) << (inner_iter - 1) / restart_dim + 1
                          << "   Iteration: " << std::setw(3) << inner_iter
                          << "   ||r|| = " << residual_norm << "\t"
                          << "||r||/||b|| = " << residual_norm / norm_rhs << std::endl;
            }
            
            this->Monitor(inner_iter, residual_norm, residual_vector, solution,
                          residual_norm <= final_norm);

            if (residual_norm <= final_norm)
            {
                update_solution(solution, outer_iter, hessenberg, s_vector,
                               preconditioned_vectors);
                
                this->final_norm = residual_norm;
                this->final_iter = inner_iter;
                this->converged = 1;

                if (print_level == 2)
                {
                    mfem::out << "FGMRES iterations: " << this->final_iter << std::endl;
                }

                for (int idx = 0; idx <= restart_dim; ++idx)
                {
                    delete krylov_vectors[idx];
                    delete preconditioned_vectors[idx];
                }
                
                return;
            }
        }

        if (print_level == 1)
        {
            mfem::out << "Restarting..." << std::endl;
        }

        update_solution(solution, outer_iter - 1, hessenberg, s_vector,
                       preconditioned_vectors);

        oper->Mult(solution, residual_vector);
        subtract(rhs, residual_vector, residual_vector);
        
        beta_norm = Norm(residual_vector);
        MFEM_ASSERT(IsFinite(beta_norm), "beta_norm = " << beta_norm);
        
        if (beta_norm <= final_norm)
        {
            this->final_norm = beta_norm;
            this->final_iter = inner_iter;
            this->converged = 1;

            if (print_level == 2)
            {
                mfem::out << "FGMRES iterations: " << this->final_iter << std::endl;
            }

            for (int idx = 0; idx <= restart_dim; ++idx)
            {
                delete krylov_vectors[idx];
                delete preconditioned_vectors[idx];
            }
            
            return;
        }
    }

    for (int idx = 0; idx <= restart_dim; ++idx)
    {
        delete krylov_vectors[idx];
        delete preconditioned_vectors[idx];
    }
    
    this->converged = 0;

    if (print_level >= 0)
    {
        mfem::out << "FGMRES: No convergence!" << std::endl;
    }

    return;
}

// ============================================================================
// Helper functions (direct copies from MFEM)
// ============================================================================
void generate_plane_rotation(double& dx, double& dy, double& cs, double& sn)
{
    if (dy == 0.0)
    {
        cs = 1.0;
        sn = 0.0;
    }
    else if (fabs(dy) > fabs(dx))
    {
        double temp = dx / dy;
        sn = 1.0 / sqrt(1.0 + temp * temp);
        cs = temp * sn;
    }
    else
    {
        double temp = dy / dx;
        cs = 1.0 / sqrt(1.0 + temp * temp);
        sn = temp * cs;
    }
}

void apply_plane_rotation(double& dx, double& dy, double& cs, double& sn)
{
    double temp = cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = temp;
}

void update_solution(
    Vector& solution,
    int dimension,
    DenseMatrix& hessenberg,
    Vector& s_vector,
    Array<Vector*>& vectors
)
{
    Vector y_vector(s_vector);

    for (int idx = dimension; idx >= 0; --idx)
    {
        y_vector(idx) /= hessenberg(idx, idx);
        for (int jdx = idx - 1; jdx >= 0; --jdx)
        {
            y_vector(jdx) -= hessenberg(jdx, idx) * y_vector(idx);
        }
    }

    for (int jdx = 0; jdx <= dimension; ++jdx)
    {
        solution.Add(y_vector(jdx), *vectors[jdx]);
    }
}
