// file: postProcessor.cpp
// ============================================================
// Post-processing implementation for EM fields
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#include "postProcessor.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <mpi.h>
#include <stdexcept>
#include <vector>

using namespace mfem;

// ============================================================================
// Construction / Destruction
// ============================================================================
PostProcessor::PostProcessor(
    ParameterHandler& input_params,
    ParMesh& mesh_instance,
    ParFiniteElementSpace& fe_space_instance,
    ParGridFunction& real_solution,
    ParGridFunction& imag_solution,
    double angular_frequency
)
    : parameters(&input_params)
    , parallel_mesh(&mesh_instance)
    , fe_space(&fe_space_instance)
    , real_field(real_solution)
    , imag_field(imag_solution)
    , frequency(angular_frequency)
{
    // Nothing else to initialize
}

PostProcessor::~PostProcessor()
{
    // Nothing to clean up
}

// ============================================================================
// Geometric utilities
// ============================================================================
Vector PostProcessor::cross_product(const Vector& vec1, const Vector& vec2)
{
    MFEM_ASSERT(vec1.Size() == 3, "cross_product expects 3D vectors.");
    MFEM_ASSERT(vec2.Size() == 3, "cross_product expects 3D vectors.");
    
    Vector result(3);
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = -vec1[0] * vec2[2] + vec1[2] * vec2[0];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    
    return result;
}

Vector PostProcessor::vector_difference(const Vector& vec1, const Vector& vec2)
{
    MFEM_ASSERT(vec1.Size() == 3, "vector_difference expects 3D vectors.");
    MFEM_ASSERT(vec2.Size() == 3, "vector_difference expects 3D vectors.");
    
    Vector result(3);
    result[0] = vec2[0] - vec1[0];
    result[1] = vec2[1] - vec1[1];
    result[2] = vec2[2] - vec1[2];
    
    return result;
}

double PostProcessor::vector_length(const Vector& vec)
{
    MFEM_ASSERT(vec.Size() == 3, "vector_length expects a 3D vector.");
    
    return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

double PostProcessor::dot_product(const Vector& vec1, const Vector& vec2)
{
    MFEM_ASSERT(vec1.Size() == 3, "dot_product expects 3D vectors.");
    MFEM_ASSERT(vec2.Size() == 3, "dot_product expects 3D vectors.");
    
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

// ============================================================================
// Point-in-tetrahedron test
// ============================================================================
PointLocation PostProcessor::check_point_in_tetrahedron(
    const Vector& point,
    int element_id
)
{
    Array<int> vertex_indices;
    parallel_mesh->GetElementVertices(element_id, vertex_indices);

    double volume = parallel_mesh->GetElementVolume(element_id);
    
    Array<double*> vertex_coords(4);
    for (int idx = 0; idx < vertex_indices.Size(); ++idx)
    {
        vertex_coords[idx] = parallel_mesh->GetVertex(vertex_indices[idx]);
    }

    Vector pa(3), pb(3), pc(3), pd(3);
    pa = vertex_coords[0];
    pb = vertex_coords[1];
    pc = vertex_coords[2];
    pd = vertex_coords[3];

    Vector ab = vector_difference(pa, pb);
    Vector ba = vector_difference(pb, pa);
    Vector ac = vector_difference(pa, pc);
    Vector ad = vector_difference(pa, pd);
    Vector bc = vector_difference(pb, pc);
    Vector cb = vector_difference(pc, pb);
    Vector bd = vector_difference(pb, pd);
    Vector db = vector_difference(pd, pb);
    Vector cd = vector_difference(pc, pd);
    
    Vector oa = vector_difference(pa, point);
    Vector ob = vector_difference(pb, point);

    Vector nabc = cross_product(ab, ac);
    Vector nacd = cross_product(ac, ad);
    Vector nadb = cross_product(ab, ad);
    Vector ncbd = cross_product(bc, bd);

    // Orient normals outward
    double value = dot_product(nabc, ad);
    if (value > 0.0) nabc.Neg();
    
    value = dot_product(nacd, ab);
    if (value > 0.0) nacd.Neg();
    
    value = dot_product(nadb, ac);
    if (value > 0.0) nadb.Neg();
    
    value = dot_product(ncbd, ba);
    if (value > 0.0) ncbd.Neg();

    double vabc = dot_product(nabc, oa);
    double vacd = dot_product(nacd, oa);
    double vadb = dot_product(nadb, oa);
    double vcbd = dot_product(ncbd, ob);
    
    double tolerance = 1e-8;

    if (vabc > 0.0 && vacd > 0.0 && vadb > 0.0 && vcbd > 0.0)
    {
        return PointLocation::INSIDE;
    }
    else if (vabc >= -tolerance && vacd >= -tolerance &&
             vadb >= -tolerance && vcbd >= -tolerance &&
             std::abs(vabc * vacd * vadb * vcbd) <= tolerance)
    {
        return PointLocation::ON_SURFACE;
    }
    else
    {
        return PointLocation::OUTSIDE;
    }
}

// ============================================================================
// Main post-processing routine
// ============================================================================
void PostProcessor::execute()
{
    // ------------------------------------------------------------------------
    // Load global station coordinates
    // ------------------------------------------------------------------------
    std::ifstream station_stream(parameters->sites_file);
    if (!station_stream)
    {
        throw std::runtime_error(
            "PostProcessor::execute(): failed to open sites file: "
            + parameters->sites_file
        );
    }
    
    int station_count;
    if (!(station_stream >> station_count) || station_count <= 0)
    {
        throw std::runtime_error(
            "PostProcessor::execute(): invalid station count in sites file: "
            + parameters->sites_file
        );
    }
    
    Array<double> global_radii(station_count), global_thetas(station_count), global_phis(station_count);
    Array<double> global_x(station_count), global_y(station_count), global_z(station_count);
    
    for (int idx = 0; idx < station_count; ++idx)
    {
        double r_val, theta_val, phi_val;
        if (!(station_stream >> r_val >> theta_val >> phi_val))
        {
            throw std::runtime_error(
                "PostProcessor::execute(): malformed station record in "
                + parameters->sites_file
            );
        }
        
        if (r_val < 0.0 || theta_val < 0.0 || theta_val > 180.0 ||
            phi_val < 0.0 || phi_val > 360.0)
        {
            throw std::invalid_argument(
                "PostProcessor::execute(): station coordinates are out of range."
            );
        }
        
        global_radii[idx] = r_val;
        global_thetas[idx] = theta_val;
        global_phis[idx] = phi_val;
        
        double theta_rad = theta_val / 180.0 * MT::kPi;
        double phi_rad = phi_val / 180.0 * MT::kPi;
        
        global_x[idx] = r_val * sin(theta_rad) * cos(phi_rad);
        global_y[idx] = r_val * sin(theta_rad) * sin(phi_rad);
        global_z[idx] = r_val * cos(theta_rad);
    }

    // ------------------------------------------------------------------------
    // Locate stations using GSLIB
    // ------------------------------------------------------------------------
    Vector coordinates(station_count * 3);
    
    for (int idx = 0; idx < station_count; ++idx)
    {
        coordinates(idx) = global_x[idx];
        coordinates(station_count + idx) = global_y[idx];
        coordinates(2 * station_count + idx) = global_z[idx];
    }

    int my_rank;
    MPI_Comm_rank(parallel_mesh->GetComm(), &my_rank);
    
    FindPointsGSLIB point_locator(parallel_mesh->GetComm());
    point_locator.Setup(*parallel_mesh);
    point_locator.FindPoints(coordinates);
    
    Array<int> status_codes = point_locator.GetCode();
    Array<int> processor_ids = point_locator.GetProc();
    Array<int> element_ids = point_locator.GetElem();

    Array<int> local_station_elements;
    
    int missing_flag = status_codes.Find(2);
    if (missing_flag != -1)
    {
        throw std::invalid_argument(
            "PostProcessor::execute(): Some stations not found!"
        );
    }
    
    for (int idx = 0; idx < station_count; ++idx)
    {
        if (processor_ids[idx] == my_rank)
        {
            local_station_elements.Append(element_ids[idx]);
            local_station_radii.Append(global_radii[idx]);
            local_station_thetas.Append(global_thetas[idx]);
            local_station_phis.Append(global_phis[idx]);
        }
    }
    
    point_locator.FreeData();

    // ------------------------------------------------------------------------
    // Compute magnetic induction B = -1/(i*omega) * Curl(E)
    // ------------------------------------------------------------------------
    int local_station_count = local_station_radii.Size();
    
    Array<double> hx_real(local_station_count);
    Array<double> hx_imag(local_station_count);
    Array<double> hy_real(local_station_count);
    Array<double> hy_imag(local_station_count);
    Array<double> hz_real(local_station_count);
    Array<double> hz_imag(local_station_count);
    
    Array<double> ex_real(local_station_count);
    Array<double> ex_imag(local_station_count);
    Array<double> ey_real(local_station_count);
    Array<double> ey_imag(local_station_count);
    Array<double> ez_real(local_station_count);
    Array<double> ez_imag(local_station_count);

    hx_real = 0.0; hx_imag = 0.0;
    hy_real = 0.0; hy_imag = 0.0;
    hz_real = 0.0; hz_imag = 0.0;
    ex_real = 0.0; ex_imag = 0.0;
    ey_real = 0.0; ey_imag = 0.0;
    ez_real = 0.0; ez_imag = 0.0;

    for (int idx = 0; idx < local_station_count; ++idx)
    {
        double radius = local_station_radii[idx];
        double theta_deg = local_station_thetas[idx];
        double phi_deg = local_station_phis[idx];
        
        double theta_rad = theta_deg / 180.0 * MT::kPi;
        double phi_rad = phi_deg / 180.0 * MT::kPi;

        int element_id = local_station_elements[idx];
        auto* transform = fe_space->GetElementTransformation(element_id);
        const auto* element = fe_space->GetFE(element_id);
        
        Vector point(3);
        point[0] = radius * sin(theta_rad) * cos(phi_rad);
        point[1] = radius * sin(theta_rad) * sin(phi_rad);
        point[2] = radius * cos(theta_rad);
        
        IntegrationPoint integration_point;
        transform->TransformBack(point, integration_point);
        transform->SetIntPoint(&integration_point);

        int num_dofs = element->GetDof();
        int dimension = element->GetDim();
        
        DenseMatrix curl_shape(num_dofs, dimension);
        element->CalcPhysCurlShape(*transform, curl_shape);
        
        DenseMatrix shape_matrix(num_dofs, dimension);
        element->CalcPhysVShape(*transform, shape_matrix);
        
        Array<int> dof_indices;
        fe_space->GetElementVDofs(element_id, dof_indices);
        
        Vector real_dofs, imag_dofs;
        real_field.GetSubVector(dof_indices, real_dofs);
        imag_field.GetSubVector(dof_indices, imag_dofs);

        Vector h_real(3), h_imag(3);
        Vector e_real(3), e_imag(3);
        
        h_real = 0.0; h_imag = 0.0;
        e_real = 0.0; e_imag = 0.0;

        double inv_frequency = 1.0 / frequency;
        
        for (int dof_idx = 0; dof_idx < num_dofs; ++dof_idx)
        {
            Vector shape_row;
            shape_matrix.GetRow(dof_idx, shape_row);
            
            Vector curl_row;
            curl_shape.GetRow(dof_idx, curl_row);
            
            for (int comp = 0; comp < 3; ++comp)
            {
                e_real[comp] += shape_row[comp] * real_dofs[dof_idx];
                e_imag[comp] += shape_row[comp] * imag_dofs[dof_idx];
                
                h_real[comp] += -inv_frequency * curl_row[comp] * imag_dofs[dof_idx];
                h_imag[comp] += inv_frequency * curl_row[comp] * real_dofs[dof_idx];
            }
        }

        hx_real[idx] = h_real[0];
        hx_imag[idx] = h_imag[0];
        hy_real[idx] = h_real[1];
        hy_imag[idx] = h_imag[1];
        hz_real[idx] = h_real[2];
        hz_imag[idx] = h_imag[2];
        
        ex_real[idx] = e_real[0];
        ex_imag[idx] = e_imag[0];
        ey_real[idx] = e_real[1];
        ey_imag[idx] = e_imag[1];
        ez_real[idx] = e_real[2];
        ez_imag[idx] = e_imag[2];
    }

    // ------------------------------------------------------------------------
    // Convert to spherical components
    // ------------------------------------------------------------------------
    h_r_real.SetSize(local_station_count);
    h_r_imag.SetSize(local_station_count);
    h_theta_real.SetSize(local_station_count);
    h_theta_imag.SetSize(local_station_count);
    h_phi_real.SetSize(local_station_count);
    h_phi_imag.SetSize(local_station_count);
    
    e_r_real.SetSize(local_station_count);
    e_r_imag.SetSize(local_station_count);
    e_theta_real.SetSize(local_station_count);
    e_theta_imag.SetSize(local_station_count);
    e_phi_real.SetSize(local_station_count);
    e_phi_imag.SetSize(local_station_count);

    h_r_real = 0.0; h_r_imag = 0.0;
    h_theta_real = 0.0; h_theta_imag = 0.0;
    h_phi_real = 0.0; h_phi_imag = 0.0;
    e_r_real = 0.0; e_r_imag = 0.0;
    e_theta_real = 0.0; e_theta_imag = 0.0;
    e_phi_real = 0.0; e_phi_imag = 0.0;

    for (int idx = 0; idx < local_station_radii.Size(); ++idx)
    {
        double theta_deg = local_station_thetas[idx];
        double phi_deg = local_station_phis[idx];

        h_r_real[idx] = sine_degrees(theta_deg) * cosine_degrees(phi_deg) * hx_real[idx] +
                        sine_degrees(theta_deg) * sine_degrees(phi_deg) * hy_real[idx] +
                        cosine_degrees(theta_deg) * hz_real[idx];
        
        h_r_imag[idx] = sine_degrees(theta_deg) * cosine_degrees(phi_deg) * hx_imag[idx] +
                        sine_degrees(theta_deg) * sine_degrees(phi_deg) * hy_imag[idx] +
                        cosine_degrees(theta_deg) * hz_imag[idx];
        
        h_theta_real[idx] = cosine_degrees(theta_deg) * cosine_degrees(phi_deg) * hx_real[idx] +
                             cosine_degrees(theta_deg) * sine_degrees(phi_deg) * hy_real[idx] -
                             sine_degrees(theta_deg) * hz_real[idx];
        
        h_theta_imag[idx] = cosine_degrees(theta_deg) * cosine_degrees(phi_deg) * hx_imag[idx] +
                             cosine_degrees(theta_deg) * sine_degrees(phi_deg) * hy_imag[idx] -
                             sine_degrees(theta_deg) * hz_imag[idx];
        
        h_phi_real[idx] = -sine_degrees(phi_deg) * hx_real[idx] + cosine_degrees(phi_deg) * hy_real[idx];
        h_phi_imag[idx] = -sine_degrees(phi_deg) * hx_imag[idx] + cosine_degrees(phi_deg) * hy_imag[idx];
    
        e_r_real[idx] = sine_degrees(theta_deg) * cosine_degrees(phi_deg) * ex_real[idx] +
                        sine_degrees(theta_deg) * sine_degrees(phi_deg) * ey_real[idx] +
                        cosine_degrees(theta_deg) * ez_real[idx];
        
        e_r_imag[idx] = sine_degrees(theta_deg) * cosine_degrees(phi_deg) * ex_imag[idx] +
                        sine_degrees(theta_deg) * sine_degrees(phi_deg) * ey_imag[idx] +
                        cosine_degrees(theta_deg) * ez_imag[idx];
        
        e_theta_real[idx] = cosine_degrees(theta_deg) * cosine_degrees(phi_deg) * ex_real[idx] +
                             cosine_degrees(theta_deg) * sine_degrees(phi_deg) * ey_real[idx] -
                             sine_degrees(theta_deg) * ez_real[idx];
        
        e_theta_imag[idx] = cosine_degrees(theta_deg) * cosine_degrees(phi_deg) * ex_imag[idx] +
                             cosine_degrees(theta_deg) * sine_degrees(phi_deg) * ey_imag[idx] -
                             sine_degrees(theta_deg) * ez_imag[idx];
        
        e_phi_real[idx] = -sine_degrees(phi_deg) * ex_real[idx] + cosine_degrees(phi_deg) * ey_real[idx];
        e_phi_imag[idx] = -sine_degrees(phi_deg) * ex_imag[idx] + cosine_degrees(phi_deg) * ey_imag[idx];
        
        // Convert Tesla to nanoTesla
        h_r_real[idx] *= 1e9;
        h_r_imag[idx] *= 1e9;
        h_theta_real[idx] *= 1e9;
        h_theta_imag[idx] *= 1e9;
        h_phi_real[idx] *= 1e9;
        h_phi_imag[idx] *= 1e9;
    }
}

// ============================================================================
// Output methods
// ============================================================================
void PostProcessor::write_local_results(std::ofstream& output_stream)
{
    for (int idx = 0; idx < local_station_radii.Size(); ++idx)
    {
        output_stream.setf(std::ios::scientific);
        output_stream.precision(6);
        
        output_stream << std::setw(15) << local_station_radii[idx] / 1000.0 << "\t"
                      << std::setw(15) << local_station_thetas[idx] << "\t"
                      << std::setw(15) << local_station_phis[idx] << "\t"
                      << std::setw(15) << h_r_real[idx] << "\t"
                      << std::setw(15) << h_r_imag[idx] << "\t"
                      << std::setw(15) << h_theta_real[idx] << "\t"
                      << std::setw(15) << h_theta_imag[idx] << "\t"
                      << std::setw(15) << h_phi_real[idx] << "\t"
                      << std::setw(15) << h_phi_imag[idx] << "\t"
                      << std::setw(15) << e_r_real[idx] << "\t"
                      << std::setw(15) << e_r_imag[idx] << "\t"
                      << std::setw(15) << e_theta_real[idx] << "\t"
                      << std::setw(15) << e_theta_imag[idx] << "\t"
                      << std::setw(15) << e_phi_real[idx] << "\t"
                      << std::setw(15) << e_phi_imag[idx] << "\n";
    }
}

void PostProcessor::save_as_single_file(std::ofstream& output_stream)
{
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Comm communicator = fe_space->GetComm();
    
    int my_rank, total_ranks;
    MPI_Comm_size(communicator, &total_ranks);
    MPI_Comm_rank(communicator, &my_rank);

    const int kFieldCount = 15;
    const int local_count = local_station_radii.Size();

    std::vector<double> local_buffer(local_count * kFieldCount);
    for (int idx = 0; idx < local_count; ++idx)
    {
        const int base = idx * kFieldCount;
        local_buffer[base + 0] = local_station_radii[idx];
        local_buffer[base + 1] = local_station_thetas[idx];
        local_buffer[base + 2] = local_station_phis[idx];
        local_buffer[base + 3] = h_r_real[idx];
        local_buffer[base + 4] = h_r_imag[idx];
        local_buffer[base + 5] = h_theta_real[idx];
        local_buffer[base + 6] = h_theta_imag[idx];
        local_buffer[base + 7] = h_phi_real[idx];
        local_buffer[base + 8] = h_phi_imag[idx];
        local_buffer[base + 9] = e_r_real[idx];
        local_buffer[base + 10] = e_r_imag[idx];
        local_buffer[base + 11] = e_theta_real[idx];
        local_buffer[base + 12] = e_theta_imag[idx];
        local_buffer[base + 13] = e_phi_real[idx];
        local_buffer[base + 14] = e_phi_imag[idx];
    }

    std::vector<int> station_counts;
    if (my_rank == 0)
    {
        station_counts.resize(total_ranks);
    }

    MPI_Gather(
        &local_count,
        1,
        MPI_INT,
        (my_rank == 0) ? station_counts.data() : nullptr,
        1,
        MPI_INT,
        0,
        communicator
    );

    std::vector<int> recv_counts;
    std::vector<int> displacements;
    std::vector<double> global_buffer;

    if (my_rank == 0)
    {
        recv_counts.resize(total_ranks);
        displacements.resize(total_ranks);

        int total_station_count = 0;
        int displacement = 0;
        for (int proc = 0; proc < total_ranks; ++proc)
        {
            recv_counts[proc] = station_counts[proc] * kFieldCount;
            displacements[proc] = displacement;
            displacement += recv_counts[proc];
            total_station_count += station_counts[proc];
        }

        global_buffer.resize(total_station_count * kFieldCount);
    }

    MPI_Gatherv(
        local_buffer.data(),
        local_count * kFieldCount,
        MPI_DOUBLE,
        (my_rank == 0) ? global_buffer.data() : nullptr,
        (my_rank == 0) ? recv_counts.data() : nullptr,
        (my_rank == 0) ? displacements.data() : nullptr,
        MPI_DOUBLE,
        0,
        communicator
    );

    if (my_rank == 0)
    {
        output_stream.setf(std::ios::scientific);
        output_stream.precision(6);

        for (int offset = 0; offset < static_cast<int>(global_buffer.size()); offset += kFieldCount)
        {
            output_stream << std::setw(15) << global_buffer[offset + 0] / 1000.0 << "\t"
                          << std::setw(15) << global_buffer[offset + 1] << "\t"
                          << std::setw(15) << global_buffer[offset + 2] << "\t"
                          << std::setw(15) << global_buffer[offset + 3] << "\t"
                          << std::setw(15) << global_buffer[offset + 4] << "\t"
                          << std::setw(15) << global_buffer[offset + 5] << "\t"
                          << std::setw(15) << global_buffer[offset + 6] << "\t"
                          << std::setw(15) << global_buffer[offset + 7] << "\t"
                          << std::setw(15) << global_buffer[offset + 8] << "\t"
                          << std::setw(15) << global_buffer[offset + 9] << "\t"
                          << std::setw(15) << global_buffer[offset + 10] << "\t"
                          << std::setw(15) << global_buffer[offset + 11] << "\t"
                          << std::setw(15) << global_buffer[offset + 12] << "\t"
                          << std::setw(15) << global_buffer[offset + 13] << "\t"
                          << std::setw(15) << global_buffer[offset + 14] << "\n";
        }
    }
}
