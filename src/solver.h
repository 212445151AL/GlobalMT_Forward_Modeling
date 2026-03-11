// file: solver.h
// ============================================================
// Solver parameters and custom FGMRES interface
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#ifndef CUSTOM_SOLVER_HEADER
#define CUSTOM_SOLVER_HEADER

#include <string>
#include <mpi.h>
#include "mfem.hpp"

// ============================================================================
// class SolverParm： Solver's parameters handler. Including MFEM's FGMRES, 
// Hypre's PCG, and Hypre's auxiliary-space Maxwell solver (AMS)
// ============================================================================
class FlexibleGMRES;

// ============================================================================
// SolverParameters - Manages solver configuration
// ============================================================================
class SolverParameters
{
public:
    // ------------------------------------------------------------------------
    // FGMRES parameters
    // ------------------------------------------------------------------------
    int fgmres_max_iterations;
    double fgmres_primal_tolerance;   // Relative tolerance for primal problem
    double fgmres_dual_tolerance;      // Relative tolerance for dual problem
    int fgmres_restart_dimension;
    int fgmres_output_level;

    // ------------------------------------------------------------------------
    // Preconditioner type
    // ------------------------------------------------------------------------
    int prec_type;   // 0: AMS, 1: PCG-AMS, 2: Multigrid

    // ------------------------------------------------------------------------
    // PCG parameters
    // ------------------------------------------------------------------------
    int pcg_max_iterations;
    double pcg_tolerance;
    int pcg_output_level;

    // ------------------------------------------------------------------------
    // AMS parameters
    // ------------------------------------------------------------------------
    int ams_cycle_type;
    int ams_max_iterations;
    double ams_tolerance;
    int ams_output_level;

    // Smoothing options
    int a_relaxation_type;
    int a_relaxation_sweeps;
    double a_relaxation_weight;
    double a_relaxation_omega;

    // AMG options (used for both B_Pi and B_G)
    int amg_coarsening_type;
    int amg_aggressive_levels;
    double amg_strength_threshold;
    int amg_interpolation_type;
    int amg_max_elements_per_row;
    int amg_relaxation_type;

    // ------------------------------------------------------------------------
    // Construction / Destruction
    // ------------------------------------------------------------------------
    explicit SolverParameters(const std::string& parameter_filename);
    ~SolverParameters();

    // ------------------------------------------------------------------------
    // Methods
    // ------------------------------------------------------------------------
    void load_parameters(const std::string& parameter_filename);
    void write_to_file(const std::string& output_filename = "input_solver_parm.log");
    
    void apply_to_pcg(mfem::HyprePCG& pcg_solver);
    void apply_to_ams(mfem::HypreAMS& ams_solver);
};

// ============================================================================
// FlexibleGMRES - Custom implementation with residual output
// 
// Based on MFEM's FGMRESSolver with modifications for relative residual output
// ============================================================================
class FlexibleGMRES : public mfem::FGMRESSolver
{
public:
    FlexibleGMRES();
    explicit FlexibleGMRES(MPI_Comm communicator);
    
    void Mult(const mfem::Vector& rhs, mfem::Vector& solution) const override;
};

// ============================================================================
// Helper functions (from MFEM)
// ============================================================================
void generate_plane_rotation(double& dx, double& dy, double& cs, double& sn);
void apply_plane_rotation(double& dx, double& dy, double& cs, double& sn);
void update_solution(
    mfem::Vector& solution,
    int dimension,
    mfem::DenseMatrix& hessenberg,
    mfem::Vector& s_vector,
    mfem::Array<mfem::Vector*>& vectors
);

#endif // CUSTOM_SOLVER_HEADER
