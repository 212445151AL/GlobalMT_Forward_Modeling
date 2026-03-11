// file: MTUtils.h
// ============================================================
// Electromagnetic utilities and namespace definitions
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#ifndef EM_UTILITIES_HEADER
#define EM_UTILITIES_HEADER

#include <complex>
#include <cmath>

namespace MT
{
    using ComplexDouble = std::complex<double>;
    static const ComplexDouble kImagUnit = ComplexDouble(0.0, 1.0);
    static const double kPi = 3.1415926535897932384626433832795;
    static const double kMuZero = 4.0 * kPi * 1e-7;
}

// ---------------------------------------------------------------------
// Trigonometric functions with degree inputs
// ---------------------------------------------------------------------
inline double cosine_degrees(double angle_deg)
{
    return std::cos(angle_deg / 180.0 * MT::kPi);
}

inline double sine_degrees(double angle_deg)
{
    return std::sin(angle_deg / 180.0 * MT::kPi);
}

inline double tangent_degrees(double angle_deg)
{
    return std::tan(angle_deg / 180.0 * MT::kPi);
}

// ---------------------------------------------------------------------
// Inverse trigonometric functions with degree outputs
// ---------------------------------------------------------------------
inline double arccos_degrees(double rad_val)
{
    return std::acos(rad_val) / MT::kPi * 180.0;
}

inline double arcsin_degrees(double rad_val)
{
    return std::asin(rad_val) / MT::kPi * 180.0;
}

inline double arctan_degrees(double rad_val)
{
    return std::atan(rad_val) / MT::kPi * 180.0;
}

inline double arctan2_degrees(double y_val, double x_val)
{
    return std::atan2(y_val, x_val) / MT::kPi * 180.0;
}

/**
 * @brief Convert spherical coordinates (r,theta,phi) to Cartesian (x,y,z)
 * 
 * @param spherical[0] = radius in meters
 * @param spherical[1] = theta in degrees [0,180]
 * @param spherical[2] = phi in degrees [0,360]
 * @param cartesian output array of size 3
 */
inline void spherical_to_cartesian(const double spherical[], double cartesian[])
{
    double radius = spherical[0];
    double theta_deg = spherical[1];
    double phi_deg = spherical[2];
    
    cartesian[0] = radius * sine_degrees(theta_deg) * cosine_degrees(phi_deg);
    cartesian[1] = radius * sine_degrees(theta_deg) * sine_degrees(phi_deg);
    cartesian[2] = radius * cosine_degrees(theta_deg);
}

/**
 * @brief Convert Cartesian (x,y,z) to spherical (r,theta,phi) coordinates
 * 
 * @param cartesian[0] = x in meters
 * @param cartesian[1] = y in meters
 * @param cartesian[2] = z in meters
 * @param spherical output array of size 3
 */
inline void cartesian_to_spherical(const double cartesian[], double spherical[])
{
    double x_coord = cartesian[0];
    double y_coord = cartesian[1];
    double z_coord = cartesian[2];
    
    spherical[0] = std::sqrt(x_coord * x_coord + y_coord * y_coord + z_coord * z_coord);
    spherical[1] = arccos_degrees(z_coord / spherical[0]);
    spherical[2] = arctan2_degrees(y_coord, x_coord);
    
    if (spherical[2] < 0.0)
    {
        spherical[2] += 360.0;
    }
}

#endif // EM_UTILITIES_HEADER
