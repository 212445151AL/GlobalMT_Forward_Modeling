// file: errorEstimators.h
// ============================================================
// Adaptive mesh refinement error estimator implementation
// ============================================================

// Implements face-jump error estimator based on normal current density
// Adapted from Zhengyong Ren et al. (2013) for secondary field formulation
// J = sigma * E_s + (sigma - sigma_0) * E_0

// Reference:
// Zhengyong Ren, Thomas Kalscheuer, Stewart Greenhalgh, Hansruedi Maurer
// "A goal-oriented adaptive finite-element approach for plane wave 3-D
// electromagnetic modelling", Geophysical Journal International, 194, 700-718

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#ifndef ERROR_ESTIMATORS_INTERFACE
#define ERROR_ESTIMATORS_INTERFACE

#include "postProcessor.h"
#include "parameterHandler.h"
#include "mfem.hpp"

// ============================================================================
// FaceJumpEstimator - Computes error based on normal current density jumps

// field formulation: J = sigma * E_s + (sigma - sigma_0) * E_0
// ============================================================================
class FaceJumpEstimator
{
public:
    // ------------------------------------------------------------------------
    // Construction / Destruction
    // ------------------------------------------------------------------------
    FaceJumpEstimator();
    ~FaceJumpEstimator();

    // ------------------------------------------------------------------------
    // Public API methods
    // ------------------------------------------------------------------------

    /**
     * @brief Compute estimated error of solution U (E field or its dual field)
     * 
     * @param params ParameterHandler
     * @param solution_field Complex grid function
     * @return Vector Estimated error vector
     */
    mfem::Vector compute_error_estimate(
        ParameterHandler& params,
        mfem::ParComplexGridFunction& solution_field
    );

    /**
     * @brief Get non-goal-oriented estimated error
     * 
     * @param params ParameterHandler
     * @param solution_field Complex grid function
     * @return Vector Estimated error vector
     */
    mfem::Vector get_error_estimate(
        ParameterHandler& params,
        mfem::ParComplexGridFunction& solution_field
    );

    /**
     * @brief Get goal-oriented estimated error
     * 
     * @param params ParameterHandler
     * @param primal_solution Primal problem solution
     * @param dual_solution Dual problem solution
     * @return Vector Goal-oriented estimated error vector
     */
    mfem::Vector get_goal_oriented_error_estimate(
        ParameterHandler& params,
        mfem::ParComplexGridFunction& primal_solution,
        mfem::ParComplexGridFunction& dual_solution
    );

private:
    // ------------------------------------------------------------------------
    // Core implementation
    // ------------------------------------------------------------------------
    
    /// Perform the actual error computation
    mfem::Vector perform_error_computation(
        ParameterHandler& params,
        mfem::ParComplexGridFunction& solution_field
    );
    
    /// Initialize a 3x3 tensor to zero
    void initialize_zero_tensor(mfem::DenseMatrix& tensor);
};

// ------------------------------------------------------------------------
// Utility: Evaluate field on a surface/interface
// ------------------------------------------------------------------------
inline void evaluate_field_on_surface(
    mfem::ParComplexGridFunction& field,
    mfem::ParFiniteElementSpace& fe_space,
    mfem::ElementTransformation& elem_transform,
    const mfem::FiniteElement& elem,
    int elem_id,
    double real_values[],
    double imag_values[],
    bool is_mpi_neighbor = false
)
{
    int num_dofs = elem.GetDof();
    int dimension = elem.GetDim();
    
    mfem::DenseMatrix shape_matrix(num_dofs, dimension);
    elem.CalcPhysVShape(elem_transform, shape_matrix);
    
    mfem::Array<int> dof_indices;
    mfem::Vector real_dofs, imag_dofs;
    
    if (!is_mpi_neighbor)
    {
        fe_space.GetElementVDofs(elem_id, dof_indices);
        field.real().GetSubVector(dof_indices, real_dofs);
        field.imag().GetSubVector(dof_indices, imag_dofs);
    }
    else
    {
        fe_space.GetFaceNbrElementVDofs(elem_id, dof_indices);
        field.real().FaceNbrData().GetSubVector(dof_indices, real_dofs);
        field.imag().FaceNbrData().GetSubVector(dof_indices, imag_dofs);
    }
    
    // Field = sum(N_j * U_j) over all DOFs
    for (int dof_idx = 0; dof_idx < num_dofs; ++dof_idx)
    {
        real_values[0] += shape_matrix(dof_idx, 0) * real_dofs[dof_idx];
        real_values[1] += shape_matrix(dof_idx, 1) * real_dofs[dof_idx];
        real_values[2] += shape_matrix(dof_idx, 2) * real_dofs[dof_idx];
        
        imag_values[0] += shape_matrix(dof_idx, 0) * imag_dofs[dof_idx];
        imag_values[1] += shape_matrix(dof_idx, 1) * imag_dofs[dof_idx];
        imag_values[2] += shape_matrix(dof_idx, 2) * imag_dofs[dof_idx];
    }
}

#endif // ERROR_ESTIMATORS_INTERFACE
