// file: parameterHandler.cpp
// ============================================================
// Parameter management implementation
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#include "parameterHandler.h"
#include "MTUtils.h"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
template <typename ValueType>
void parse_config_line(
    const std::vector<std::string>& parameter_lines,
    int line_index,
    ValueType& value,
    const std::string& field_name
)
{
    std::istringstream line_stream(parameter_lines[line_index]);
    if (!(line_stream >> value))
    {
        throw std::runtime_error(
            "ParameterHandler::read_parameter_file(): failed to parse field '"
            + field_name + "' at line index " + std::to_string(line_index)
        );
    }
}
}

// ============================================================================
// Construction / Destruction
// ============================================================================
ParameterHandler::ParameterHandler(const char* parameter_filename)
{
    read_parameter_file(parameter_filename);
}

ParameterHandler::~ParameterHandler()
{
    // Nothing to clean up
}

// ============================================================================
// File parsing utilities
// ============================================================================
void ParameterHandler::skip_comments(
    std::istream& input_stream,
    std::vector<std::string>& output_lines,
    const std::string& comment_marker
)
{
    std::string line;
    
    while (std::getline(input_stream, line))
    {
        for (char& character : line)
        {
            if (character == '\t' || character == ',' ||
                character == ';' || character == '\r' || character == '\n')
            {
                character = ' ';
            }
        }
        
        line.erase(0, line.find_first_not_of(" "));
        line.erase(line.find_last_not_of(" ") + 1);
        
        size_t comment_pos = line.find_first_of(comment_marker);
        if (comment_pos != std::string::npos)
        {
            line.erase(comment_pos);
        }
        
        if (line.empty())
        {
            continue;
        }
        
        output_lines.push_back(line);
    }
}

// ============================================================================
// Main parameter reading
// ============================================================================
void ParameterHandler::read_parameter_file(const char* parameter_filename)
{
    std::ifstream input_stream(parameter_filename);
    if (!input_stream)
    {
        throw std::runtime_error(
            "ParameterHandler::read_parameter_file(): failed to open parameter file: "
            + std::string(parameter_filename)
        );
    }
    
    const int expected_lines = 14;
    std::vector<std::string> parameter_lines;
    
    skip_comments(input_stream, parameter_lines);
    if (static_cast<int>(parameter_lines.size()) != expected_lines)
    {
        throw std::runtime_error(
            "ParameterHandler::read_parameter_file(): expected "
            + std::to_string(expected_lines)
            + " non-comment lines but got "
            + std::to_string(parameter_lines.size())
        );
    }
    
    input_stream.close();

    // Parse parameters
    parse_config_line(parameter_lines, 0, source_file1, "source_file1");
    parse_config_line(parameter_lines, 1, source_file2, "source_file2");
    parse_config_line(parameter_lines, 2, source_file3, "source_file3");
    parse_config_line(parameter_lines, 3, polynomial_order, "polynomial_order");
    parse_config_line(parameter_lines, 4, linear_opts_file, "linear_opts_file");
    parse_config_line(parameter_lines, 5, conductivity_model_file, "conductivity_model_file");
    parse_config_line(parameter_lines, 6, mesh_file, "mesh_file");
    parse_config_line(parameter_lines, 7, sites_file, "sites_file");
    parse_config_line(parameter_lines, 8, boundary_marker, "boundary_marker");
    parse_config_line(parameter_lines, 9, amr_method, "amr_method");
    parse_config_line(parameter_lines, 10, max_refinement_iterations, "max_refinement_iterations");
    parse_config_line(parameter_lines, 11, max_dofs, "max_dofs");
    parse_config_line(parameter_lines, 12, refinement_threshold_ratio, "refinement_threshold_ratio");
    parse_config_line(parameter_lines, 13, print_vtk, "print_vtk");

    write_to_file();

    // ------------------------------------------------------------------------
    // Read conductivity model
    // ------------------------------------------------------------------------
    std::ifstream cond_stream(conductivity_model_file.c_str());
    if (!cond_stream)
    {
        throw std::runtime_error(
            "ParameterHandler::read_parameter_file(): failed to open conductivity model file: "
            + conductivity_model_file
        );
    }
    
    if (!(cond_stream >> marker_type))
    {
        throw std::runtime_error(
            "ParameterHandler::read_parameter_file(): failed to read marker_type from "
            + conductivity_model_file
        );
    }

    std::transform(
        marker_type.begin(),
        marker_type.end(),
        marker_type.begin(),
        [](unsigned char character)
        {
            return static_cast<char>(std::tolower(character));
        }
    );
    
    if (marker_type != "region_marker" && marker_type != "element_id")
    {
        throw std::invalid_argument(
            "ParameterHandler::read_parameter_file(): invalid marker_type: " + marker_type
        );
    }
    
    if (!(cond_stream >> num_regions) || num_regions <= 0)
    {
        throw std::runtime_error(
            "ParameterHandler::read_parameter_file(): invalid region count in "
            + conductivity_model_file
        );
    }
    
    for (int idx = 0; idx < num_regions; ++idx)
    {
        int marker_value;
        std::vector<double> conductivity_tensor(9);
        double tag_value;
        
        if (!(cond_stream >> marker_value))
        {
            throw std::runtime_error(
                "ParameterHandler::read_parameter_file(): failed to read marker at region index "
                + std::to_string(idx)
            );
        }
        
        for (int comp = 0; comp < 9; ++comp)
        {
            if (!(cond_stream >> conductivity_tensor[comp]))
            {
                throw std::runtime_error(
                    "ParameterHandler::read_parameter_file(): failed to read conductivity tensor at region index "
                    + std::to_string(idx)
                );
            }
        }
        
        if (!(cond_stream >> tag_value))
        {
            throw std::runtime_error(
                "ParameterHandler::read_parameter_file(): failed to read tag value at region index "
                + std::to_string(idx)
            );
        }
        
        marker_list.push_back(marker_value);
        region_conductivity[marker_value] = conductivity_tensor;
        region_tag[marker_value] = tag_value;
    }
    
    cond_stream.close();

    // ------------------------------------------------------------------------
    // Verify input files exist
    // ------------------------------------------------------------------------
    auto check_readable = [](const std::string& file_path, const std::string& label)
    {
        std::ifstream check_stream(file_path);
        if (!check_stream)
        {
            throw std::runtime_error(
                "ParameterHandler::read_parameter_file(): cannot open " + label + ": " + file_path
            );
        }
    };

    check_readable(source_file1, "source_file1");
    check_readable(source_file2, "source_file2");
    check_readable(source_file3, "source_file3");
    check_readable(mesh_file, "mesh_file");
    check_readable(linear_opts_file, "linear_opts_file");
    check_readable(sites_file, "sites_file");
}

// ============================================================================
// Data access methods
// ============================================================================
void ParameterHandler::get_element_conductivity_tensor(
    int marker,
    mfem::DenseMatrix& output_tensor
)
{
    if (region_conductivity.empty())
    {
        throw std::logic_error(
            "ParameterHandler::get_element_conductivity_tensor(): conductivity map is empty."
        );
    }
    
    auto iterator = region_conductivity.find(marker);
    if (iterator == region_conductivity.end())
    {
        throw std::out_of_range(
            "ParameterHandler::get_element_conductivity_tensor(): missing marker "
            + std::to_string(marker)
        );
    }
    
    std::vector<double>& tensor_values = iterator->second;
    
    output_tensor(0, 0) = tensor_values[0];
    output_tensor(0, 1) = tensor_values[1];
    output_tensor(0, 2) = tensor_values[2];
    output_tensor(1, 0) = tensor_values[3];
    output_tensor(1, 1) = tensor_values[4];
    output_tensor(1, 2) = tensor_values[5];
    output_tensor(2, 0) = tensor_values[6];
    output_tensor(2, 1) = tensor_values[7];
    output_tensor(2, 2) = tensor_values[8];
}

double ParameterHandler::get_element_tag(int marker)
{
    if (region_tag.empty())
    {
        throw std::logic_error(
            "ParameterHandler::get_element_tag(): region_tag map is empty."
        );
    }
    
    auto iterator = region_tag.find(marker);
    if (iterator == region_tag.end())
    {
        throw std::out_of_range(
            "ParameterHandler::get_element_tag(): missing marker "
            + std::to_string(marker)
        );
    }
    
    return iterator->second;
}

// ============================================================================
// Output methods
// ============================================================================
void ParameterHandler::write_to_file(const std::string& output_filename)
{
    std::ofstream output_stream(output_filename);
    if (!output_stream)
    {
        throw std::runtime_error(
            "ParameterHandler::write_to_file(): failed to open output file: "
            + output_filename
        );
    }
    
    output_stream << source_file1 << "\n";
    output_stream << source_file2 << "\n";
    output_stream << source_file3 << "\n";
    output_stream << polynomial_order << "\n";
    output_stream << linear_opts_file << "\n";
    output_stream << conductivity_model_file << "\n";
    output_stream << mesh_file << "\n";
    output_stream << sites_file << "\n";
    output_stream << boundary_marker << "\n";
    output_stream << amr_method << "\n";
    output_stream << max_refinement_iterations << "\n";
    output_stream << max_dofs << "\n";
    output_stream << refinement_threshold_ratio << "\n";
    output_stream << print_vtk << "\n";
}
