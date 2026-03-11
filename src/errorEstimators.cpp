// file: errorEstimators.cpp
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

#include "errorEstimators.h"

using namespace mfem;

// ============================================================================
// Constructor / Destructor
// ============================================================================
FaceJumpEstimator::FaceJumpEstimator()
{
    // Nothing to initialize
}

FaceJumpEstimator::~FaceJumpEstimator()
{
    // Nothing to clean up
}

// ============================================================================
// Public interface methods
// ============================================================================
Vector FaceJumpEstimator::compute_error_estimate(
    ParameterHandler& params,
    ParComplexGridFunction& solution_field)
{
    return perform_error_computation(params, solution_field);
}

Vector FaceJumpEstimator::get_error_estimate(
    ParameterHandler& params,
    ParComplexGridFunction& solution_field)
{
    return perform_error_computation(params, solution_field);
}

Vector FaceJumpEstimator::get_goal_oriented_error_estimate(
    ParameterHandler& params,
    ParComplexGridFunction& primal_solution,
    ParComplexGridFunction& dual_solution)
{
    Vector primal_error = perform_error_computation(params, primal_solution);
    Vector dual_error = perform_error_computation(params, dual_solution);
    
    int num_elements = primal_error.Size();
    Vector goal_error(num_elements);
    
    for (int elem_idx = 0; elem_idx < num_elements; ++elem_idx)
    {
        goal_error[elem_idx] = primal_error[elem_idx] * dual_error[elem_idx];
    }
    
    return goal_error;
}

// ============================================================================
// Core error computation implementation
// ============================================================================
Vector FaceJumpEstimator::perform_error_computation(
    ParameterHandler& params,
    ParComplexGridFunction& solution_field)
{
    auto* fe_space = solution_field.ParFESpace();
    auto* parallel_mesh = fe_space->GetParMesh();
    int spatial_dim = parallel_mesh->SpaceDimension();
    int elem_count = parallel_mesh->GetNE();
    
    Vector element_errors(elem_count);
    element_errors = 0.0;

    // ------------------------------------------------------------------------
    // Process interior faces
    // ------------------------------------------------------------------------
    int total_faces = parallel_mesh->GetNFaces();
    for (int face_idx = 0; face_idx < total_faces; ++face_idx)
    {
        auto* face_transform = parallel_mesh->GetInteriorFaceTransformations(face_idx);
        if (face_transform == nullptr)
        {
            continue;
        }
        
        DenseMatrix conductivity_tensor1(spatial_dim, spatial_dim);
        DenseMatrix conductivity_tensor2(spatial_dim, spatial_dim);
        
        initialize_zero_tensor(conductivity_tensor1);
        initialize_zero_tensor(conductivity_tensor2);

        auto* face_trans = face_transform->Face;
        auto* face_rule = &IntRules.Get(
            face_trans->GetGeometryType(),
            2 * face_trans->Order() - 1
        );
        
        // Element 1 properties
        int elem1_id = face_transform->Elem1No;
        auto* elem1_fe = fe_space->GetFE(elem1_id);
        auto* elem1_trans = face_transform->Elem1;
        params.get_element_conductivity_tensor(elem1_trans->Attribute, conductivity_tensor1);
        
        // Element 2 properties
        int elem2_id = face_transform->Elem2No;
        auto* elem2_fe = fe_space->GetFE(elem2_id);
        auto* elem2_trans = face_transform->Elem2;
        params.get_element_conductivity_tensor(elem2_trans->Attribute, conductivity_tensor2);

        // Integrate over face quadrature points
        for (int quad_idx = 0; quad_idx < face_rule->GetNPoints(); ++quad_idx)
        {
            const auto& int_point = face_rule->IntPoint(quad_idx);
            face_trans->SetIntPoint(&int_point);
            
            Vector jacobian_normal(spatial_dim);
            CalcOrtho(face_trans->Jacobian(), jacobian_normal);
            
            Vector unit_normal(spatial_dim);
            unit_normal.Set(1.0 / jacobian_normal.Norml2(), jacobian_normal);
            
            double weight_factor = face_trans->Weight() * int_point.weight;

            // Evaluate fields on both sides
            IntegrationPoint elem1_ip;
            face_transform->Loc1.Transform(int_point, elem1_ip);
            elem1_trans->SetIntPoint(&elem1_ip);
            
            double field1_real[3] = {0.0, 0.0, 0.0};
            double field1_imag[3] = {0.0, 0.0, 0.0};
            evaluate_field_on_surface(
                solution_field, *fe_space, *elem1_trans, *elem1_fe,
                elem1_id, field1_real, field1_imag
            );

            IntegrationPoint elem2_ip;
            face_transform->Loc2.Transform(int_point, elem2_ip);
            elem2_trans->SetIntPoint(&elem2_ip);
            
            double field2_real[3] = {0.0, 0.0, 0.0};
            double field2_imag[3] = {0.0, 0.0, 0.0};
            evaluate_field_on_surface(
                solution_field, *fe_space, *elem2_trans, *elem2_fe,
                elem2_id, field2_real, field2_imag
            );

            // Form complex field vectors
            MT::ComplexDouble field1_complex[3], field2_complex[3];
            for (int comp = 0; comp < 3; ++comp)
            {
                field1_complex[comp] = MT::ComplexDouble(field1_real[comp], field1_imag[comp]);
                field2_complex[comp] = MT::ComplexDouble(field2_real[comp], field2_imag[comp]);
            }

            // Compute current densities using tensor conductivity
            MT::ComplexDouble current_density1[3] = {0.0, 0.0, 0.0};
            MT::ComplexDouble current_density2[3] = {0.0, 0.0, 0.0};
            
            for (int row = 0; row < 3; ++row)
            {
                for (int col = 0; col < 3; ++col)
                {
                    current_density1[row] += conductivity_tensor1(row, col) * field1_complex[col];
                    current_density2[row] += conductivity_tensor2(row, col) * field2_complex[col];
                }
            }

            // Compute normal component jump
            MT::ComplexDouble normal_jump = 0.0;
            for (int comp = 0; comp < 3; ++comp)
            {
                normal_jump += unit_normal[comp] * (current_density1[comp] - current_density2[comp]);
            }

            // Accumulate error contribution
            double integral_contrib = 0.5 * weight_factor * std::pow(std::abs(normal_jump), 2);
            element_errors[elem1_id] += integral_contrib;
            element_errors[elem2_id] += integral_contrib;
        }
    }

    // ------------------------------------------------------------------------
    // Process shared faces (MPI boundary)
    // ------------------------------------------------------------------------
    solution_field.real().ExchangeFaceNbrData();
    solution_field.imag().ExchangeFaceNbrData();
    
    int shared_face_count = parallel_mesh->GetNSharedFaces();
    for (int face_idx = 0; face_idx < shared_face_count; ++face_idx)
    {
        auto* face_transform = parallel_mesh->GetSharedFaceTransformations(face_idx);
        auto* face_trans = face_transform->Face;
        auto* face_rule = &IntRules.Get(
            face_trans->GetGeometryType(),
            2 * face_trans->Order() - 1
        );
        
        DenseMatrix conductivity_tensor1(spatial_dim, spatial_dim);
        DenseMatrix conductivity_tensor2(spatial_dim, spatial_dim);
        initialize_zero_tensor(conductivity_tensor1);
        initialize_zero_tensor(conductivity_tensor2);

        // Element 1 
        int elem1_id = face_transform->Elem1No;
        auto* elem1_fe = fe_space->GetFE(elem1_id);
        auto* elem1_trans = face_transform->Elem1;
        params.get_element_conductivity_tensor(elem1_trans->Attribute, conductivity_tensor1);

        // Element 2 
        int elem2_id = face_transform->Elem2No - parallel_mesh->GetNE();
        auto* elem2_fe = fe_space->GetFaceNbrFE(elem2_id);
        auto* elem2_trans = face_transform->Elem2;
        params.get_element_conductivity_tensor(elem2_trans->Attribute, conductivity_tensor2);

        for (int quad_idx = 0; quad_idx < face_rule->GetNPoints(); ++quad_idx)
        {
            const auto& int_point = face_rule->IntPoint(quad_idx);
            face_trans->SetIntPoint(&int_point);
            
            Vector jacobian_normal(spatial_dim);
            CalcOrtho(face_trans->Jacobian(), jacobian_normal);
            
            Vector unit_normal(spatial_dim);
            unit_normal.Set(1.0 / jacobian_normal.Norml2(), jacobian_normal);
            
            double weight_factor = face_trans->Weight() * int_point.weight;

            // Evaluate on element 1
            IntegrationPoint elem1_ip;
            face_transform->Loc1.Transform(int_point, elem1_ip);
            elem1_trans->SetIntPoint(&elem1_ip);
            
            double field1_real[3] = {0.0, 0.0, 0.0};
            double field1_imag[3] = {0.0, 0.0, 0.0};
            evaluate_field_on_surface(
                solution_field, *fe_space, *elem1_trans, *elem1_fe,
                elem1_id, field1_real, field1_imag
            );

            // Evaluate on element 2
            IntegrationPoint elem2_ip;
            face_transform->Loc2.Transform(int_point, elem2_ip);
            elem2_trans->SetIntPoint(&elem2_ip);
            
            double field2_real[3] = {0.0, 0.0, 0.0};
            double field2_imag[3] = {0.0, 0.0, 0.0};
            evaluate_field_on_surface(
                solution_field, *fe_space, *elem2_trans, *elem2_fe,
                elem2_id, field2_real, field2_imag, true
            );

            MT::ComplexDouble field1_complex[3], field2_complex[3];
            for (int comp = 0; comp < 3; ++comp)
            {
                field1_complex[comp] = MT::ComplexDouble(field1_real[comp], field1_imag[comp]);
                field2_complex[comp] = MT::ComplexDouble(field2_real[comp], field2_imag[comp]);
            }

            MT::ComplexDouble current_density1[3] = {0.0, 0.0, 0.0};
            MT::ComplexDouble current_density2[3] = {0.0, 0.0, 0.0};
            
            for (int row = 0; row < 3; ++row)
            {
                for (int col = 0; col < 3; ++col)
                {
                    current_density1[row] += conductivity_tensor1(row, col) * field1_complex[col];
                    current_density2[row] += conductivity_tensor2(row, col) * field2_complex[col];
                }
            }

            MT::ComplexDouble normal_jump = 0.0;
            for (int comp = 0; comp < 3; ++comp)
            {
                normal_jump += unit_normal[comp] * (current_density1[comp] - current_density2[comp]);
            }

            double integral_contrib = 0.5 * weight_factor * std::pow(std::abs(normal_jump), 2);
            element_errors[elem1_id] += integral_contrib;
        }
    }

    // ------------------------------------------------------------------------
    // Finalize error values (take square root)
    // ------------------------------------------------------------------------
    Vector final_errors(elem_count);
    final_errors = 0.0;
    
    for (int elem_idx = 0; elem_idx < elem_count; ++elem_idx)
    {
        final_errors[elem_idx] = std::sqrt(element_errors[elem_idx]);
    }

    return final_errors;
}

// ============================================================================
// Helper method to zero-initialize tensor
// ============================================================================
void FaceJumpEstimator::initialize_zero_tensor(DenseMatrix& tensor)
{
    for (int row = 0; row < 3; ++row)
    {
        for (int col = 0; col < 3; ++col)
        {
            tensor(row, col) = 0.0;
        }
    }
}
