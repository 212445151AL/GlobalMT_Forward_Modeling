// file: postProcessor.h
// ============================================================
// Post-processing interface for EM fields
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#ifndef POST_PROCESSOR_HEADER
#define POST_PROCESSOR_HEADER

#include <fstream>
#include "MTUtils.h"
#include "parameterHandler.h"
#include "mfem.hpp"

// ------------------------------------------------------------------------
// Point location enumeration
// ------------------------------------------------------------------------
enum class PointLocation
{
    INSIDE,
    OUTSIDE,
    ON_SURFACE
};

// ============================================================================
// PostProcessor - Computes magnetic field B from electric field E
// ============================================================================
class PostProcessor
{
public:
    // ------------------------------------------------------------------------
    // Construction / Destruction
    // ------------------------------------------------------------------------
    PostProcessor(
        ParameterHandler& input_params,
        mfem::ParMesh& mesh_instance,
        mfem::ParFiniteElementSpace& fe_space_instance,
        mfem::ParGridFunction& real_solution,
        mfem::ParGridFunction& imag_solution,
        double angular_frequency
    );
    
    ~PostProcessor();

    // ------------------------------------------------------------------------
    // Geometric utilities
    // ------------------------------------------------------------------------
    mfem::Vector cross_product(const mfem::Vector& vec1, const mfem::Vector& vec2);
    mfem::Vector vector_difference(const mfem::Vector& vec1, const mfem::Vector& vec2);
    double vector_length(const mfem::Vector& vec);
    double dot_product(const mfem::Vector& vec1, const mfem::Vector& vec2);
    
    // ------------------------------------------------------------------------
    // Point location
    // ------------------------------------------------------------------------
    PointLocation check_point_in_tetrahedron(const mfem::Vector& point, int element_id);

    // ------------------------------------------------------------------------
    // Main processing
    // ------------------------------------------------------------------------
    void execute();

    // ------------------------------------------------------------------------
    // Output methods
    // ------------------------------------------------------------------------
    void write_local_results(std::ofstream& output_stream);
    void save_as_single_file(std::ofstream& output_stream);

    // ------------------------------------------------------------------------
    // Public member variables
    // ------------------------------------------------------------------------
    ParameterHandler* parameters;
    mfem::ParMesh* parallel_mesh;
    mfem::ParFiniteElementSpace* fe_space;
    mfem::ParGridFunction& real_field;
    mfem::ParGridFunction& imag_field;
    double frequency;
    
    int station_count;
    mfem::Array<double> global_radii, global_thetas, global_phis;
    mfem::Array<double> global_x, global_y, global_z;

    mfem::Array<double> local_station_radii, local_station_thetas, local_station_phis;

    // Field components in spherical coordinates (Re and Im)
    mfem::Array<double> h_r_real, h_r_imag;
    mfem::Array<double> h_theta_real, h_theta_imag;
    mfem::Array<double> h_phi_real, h_phi_imag;
    
    mfem::Array<double> e_r_real, e_r_imag;
    mfem::Array<double> e_theta_real, e_theta_imag;
    mfem::Array<double> e_phi_real, e_phi_imag;
};

#endif // POST_PROCESSOR_HEADER
