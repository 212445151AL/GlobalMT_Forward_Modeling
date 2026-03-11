// file: MTCurrentSource.h
// ============================================================
// Motional induction current source interface
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================
#ifndef MOTIONAL_CURRENT_SOURCE_HEADER
#define MOTIONAL_CURRENT_SOURCE_HEADER

#include <string>
#include <vector>
#include <map>

// ============================================================================
// Source - Represents the current source for MT 
// ============================================================================
class Source
{
public:
    // ------------------------------------------------------------------------
    // Construction / Destruction
    // ------------------------------------------------------------------------
    explicit Source(const std::string& source_filename);
    ~Source();

    /**
     * @brief Load motional current data from file
     * 
     * Reads the current density values for each element from the input file
     * and populates the member vectors.
     */
    void load_source_data();

    // ------------------------------------------------------------------------
    // Public member variables
    // ------------------------------------------------------------------------
    std::string source_file;
    std::string source_name;
    double period_hours;                // Period in hours
    double angular_frequency;           // Angular frequency (rad/s)
    int element_count;                  // Number of elements with current
    
    std::vector<int> element_ids;       // Global element IDs
    std::map<int, int> element_id_to_index;  // Map from element ID to index
    
    std::vector<double> jx_real, jx_imag;
    std::vector<double> jy_real, jy_imag;
    std::vector<double> jz_real, jz_imag;
};

#endif // MOTIONAL_CURRENT_SOURCE_HEADER
