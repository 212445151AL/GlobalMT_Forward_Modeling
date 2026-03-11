// file: MTCurrentSource.cpp
// ============================================================
// Motional induction current source implementation
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#include "MTCurrentSource.h"
#include "MTUtils.h"
#include <fstream>
#include <stdexcept>

// ============================================================================
// Construction / Destruction
// ============================================================================
Source::Source(const std::string& source_filename)
    : source_file(source_filename)
{
    load_source_data();
}

Source::~Source()
{
    // Nothing to clean up
}

// ============================================================================
// Data loading
// ============================================================================
void Source::load_source_data()
{
    std::ifstream input_stream(source_file);
    if (!input_stream)
    {
        throw std::runtime_error(
            "Source::load_source_data(): failed to open source file: " + source_file
        );
    }
    
    if (!(input_stream >> source_name >> period_hours) || period_hours <= 0.0)
    {
        throw std::runtime_error(
            "Source::load_source_data(): failed to read source name / period from "
            + source_file
        );
    }
    
    double period_seconds = period_hours * 3600.0;
    angular_frequency = 2.0 * MT::kPi / period_seconds;
    
    if (!(input_stream >> element_count) || element_count <= 0)
    {
        throw std::runtime_error(
            "Source::load_source_data(): invalid element count in " + source_file
        );
    }
    
    element_ids.resize(element_count);
    jx_real.resize(element_count);
    jx_imag.resize(element_count);
    jy_real.resize(element_count);
    jy_imag.resize(element_count);
    jz_real.resize(element_count);
    jz_imag.resize(element_count);
    
    for (int idx = 0; idx < element_count; ++idx)
    {
        if (!(input_stream >> element_ids[idx]
                          >> jx_real[idx] >> jx_imag[idx]
                          >> jy_real[idx] >> jy_imag[idx]
                          >> jz_real[idx] >> jz_imag[idx]))
        {
            throw std::runtime_error(
                "Source::load_source_data(): malformed source record at index "
                + std::to_string(idx) + " in " + source_file
            );
        }
        
        const std::pair<std::map<int, int>::iterator, bool> insert_result =
            element_id_to_index.insert(std::make_pair(element_ids[idx], idx));
        if (!insert_result.second)
        {
            throw std::invalid_argument(
                "Source::load_source_data(): duplicate element id "
                + std::to_string(element_ids[idx]) + " in " + source_file
            );
        }
    }
    
    input_stream.close();
}
