// file: config.h
// ============================================================
// Configuration settings for external library integration
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#ifndef CONFIG_HEADER_FILE
#define CONFIG_HEADER_FILE

// =========================================================================
// GSLIB integration for point location
// Source: https://github.com/Nek5000/gslib
//
// This library provides robust interpolation capabilities.
// Specifically, GSLIB-FindPoints is used to find the elements 
// that contain the measuring point
// =========================================================================
#define ENABLE_GSLIB

// TODO: Implement MUMPS direct solver integration

#endif // CONFIG_HEADER_FILE