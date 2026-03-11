// file: parameterHandler.h
// ============================================================
// Parameter management interface
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#ifndef PARAMETER_HANDLER_HEADER
#define PARAMETER_HANDLER_HEADER

#include <fstream>
#include <string>
#include <map>
#include <vector>
#include "mfem.hpp"

namespace mfem
{
    class DenseMatrix;
}

/**
 * @class ParamAdministrator
 * @brief Class for managing all input parameters and conductivity models
 * 
 * This class reads and stores all input parameters from configuration files,
 * including mesh files, source files, conductivity models, and AMR settings.
 */
class ParameterHandler
{
public:
    // ------------------------------------------------------------------------
    // Construction / Destruction
    // ------------------------------------------------------------------------
    explicit ParameterHandler(const char* parameter_filename);
    ~ParameterHandler();

    /**
     * @brief Skip comments and empty lines in input stream
     * 
     * @param inputStream Input stream to read from
     * @param parameterStringVector Vector to store non-comment lines
     * @param commentString Comment identifier (default: "#")
     */
    void skip_comments(
        std::istream& input_stream,
        std::vector<std::string>& output_lines,
        const std::string& comment_marker = "#"
    );
    
    /**
     * @brief Read model information from parameter file
     * @param parameterFile Path to parameter file
     */
    void read_parameter_file(const char* parameter_filename);

    /**
     * @brief Get conductivity tensor for a given marker (attribute)
     * @param marker Element marker/attribute
     * @param output_tensor Output conductivity tensor
     */
    void get_element_conductivity_tensor(
        int marker,
        mfem::DenseMatrix& output_tensor
    );
 
    /**
     * @brief Get tag value for a given marker (attribute)
     * @param marker Element marker/attribute
     * @return double Tag value
     */
    double get_element_tag(int marker);

    /**
     * @brief Output all parameters to file for verification
     * @param outputFileName Output file name (default: "input_parameter_list.log")
     */
    void write_to_file(const std::string& output_filename = "input_parameter_list.log");

    // ------------------------------------------------------------------------
    // Input parameters from configuration file
    // ------------------------------------------------------------------------
    std::string source_file1;
    std::string source_file2;
    std::string source_file3;
    
    int polynomial_order;                      // Polynomial order
    std::string linear_opts_file;      // FGMRES and AMS options
    std::string conductivity_model_file; // Conductivity model
    std::string mesh_file;              // Mesh file
    std::string sites_file;              // Station locations
    int boundary_marker;                 // Boundary marker for Dirichlet BC

    // AMR parameters
    int amr_method;    // 0: goal-oriented, 1: non-goal-oriented, 2: global
    int max_refinement_iterations;        // Number of h-refinement iterations
    int max_dofs;      // Maximum allowed DOFs
    double refinement_threshold_ratio;       // Marking parameter for refinement
    int print_vtk;     // Whether to output VTK files

    // ------------------------------------------------------------------------
    // Model parameters
    // ------------------------------------------------------------------------
    std::string marker_type;           // "region_marker" or "element_id"
    int num_regions;                    // Number of regions/elements
    std::vector<int> marker_list;       // List of markers
    
    std::map<int, std::vector<double>> region_conductivity;  // Conductivity tensors
    std::map<int, double> region_tag;   // Region tags for measurement points
};

#endif // PARAMETER_HANDLER_HEADER
