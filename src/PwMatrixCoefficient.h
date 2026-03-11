// file: PwMatrixCoefficient.h
// ============================================================
// Piecewise matrix coefficient for element-wise varying properties
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#ifndef PIECEWISE_MATRIX_COEFFICIENT_HEADER
#define PIECEWISE_MATRIX_COEFFICIENT_HEADER

#include "mfem.hpp"
#include <map>

// ============================================================================
// PiecewiseMatrixCoefficient - Returns different matrix coefficients
// based on mesh element attributes
// ============================================================================
class PiecewiseMatrixCoefficient : public mfem::MatrixCoefficient
{
private:
    std::map<int, mfem::MatrixCoefficient*> attribute_to_coefficient;
    mfem::DenseMatrix default_matrix;

public:
    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    /// @param rows Number of rows
    /// @param cols Number of columns
    PiecewiseMatrixCoefficient(int rows, int cols)
        : mfem::MatrixCoefficient(rows, cols), default_matrix(rows, cols)
    {
        default_matrix = 0.0;
        for (int idx = 0; idx < rows && idx < cols; ++idx)
        {
            default_matrix(idx, idx) = 1.0;
        }
    }

    // ------------------------------------------------------------------------
    // Coefficient management
    // ------------------------------------------------------------------------
    /// @brief Add a coefficient for a specific attribute
    /// @param attribute Mesh element attribute
    /// @param coefficient Matrix coefficient to use
    void add_coefficient(int attribute, mfem::MatrixCoefficient* coefficient)
    {
        attribute_to_coefficient[attribute] = coefficient;
    }

    /// @brief Set the default matrix for unassigned attributes
    /// @param matrix Default matrix
    void set_default_matrix(const mfem::DenseMatrix& matrix)
    {
        default_matrix = matrix;
    }

    // ------------------------------------------------------------------------
    // Evaluation
    // ------------------------------------------------------------------------
    void Eval(
        mfem::DenseMatrix& output_matrix,
        mfem::ElementTransformation& transform,

        
        const mfem::IntegrationPoint& point
    ) override
    {
        int attribute = transform.Attribute;
        
        auto iterator = attribute_to_coefficient.find(attribute);
        
        if (iterator != attribute_to_coefficient.end())
        {
            iterator->second->Eval(output_matrix, transform, point);
        }
        else
        {
            output_matrix = default_matrix;
        }
    }

    // ------------------------------------------------------------------------
    // Destructor
    // ------------------------------------------------------------------------
    ~PiecewiseMatrixCoefficient() override
    {
        for (auto& pair : attribute_to_coefficient)
        {
            delete pair.second;
        }
    }
};

#endif // PIECEWISE_MATRIX_COEFFICIENT_HEADER
