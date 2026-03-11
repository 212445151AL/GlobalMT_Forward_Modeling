// file: GlobalMT.cpp
// ============================================================
// Main driver for global motional induction forward modeling
// ============================================================

// Author  :    Liangyu Xie,Chaojian Chen,Juanyu Wang
//             
// Institute:   Central South University (CSU)
// Email   :    8211221219@csu.edu.cn
// Date    :    2026/01/05

// GitHub Page: https://github.com/212445151AL
// ============================================================

#include <mpi.h>
#include <sys/stat.h>
#include "MTUtils.h"
#include "parameterHandler.h"
#include "MTCurrentSource.h"
#include "fem.h"
#include "postProcessor.h"
#include "mfem.hpp"

using namespace mfem;

// ============================================================================
// Main program entry point
// ============================================================================
int main(int argc, char* argv[])
{
    // ------------------------------------------------------------------------
    // MPI initialization
    // ------------------------------------------------------------------------
    int num_processes, process_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

    double total_start_time = MPI_Wtime();

    // ------------------------------------------------------------------------
    // Check command line arguments
    // ------------------------------------------------------------------------
    if (argc < 2)
    {
        if (process_rank == 0)
        {
            std::cout << "Usage: " << argv[0] << " input_model_filename\n";
        }
        MPI_Finalize();
        return 1;
    }

    // ------------------------------------------------------------------------
    // Step 1: Load model parameters
    // ------------------------------------------------------------------------
    ParameterHandler model_parameters(argv[1]);

    Source mt_source1(model_parameters.source_file1);
    Source mt_source2(model_parameters.source_file2);
    Source mt_source3(model_parameters.source_file3);

    int polynomial_order = model_parameters.polynomial_order;
    
    if (process_rank == 0)
    {
        std::cout << "\n*************** Motional Induction Forward Modeling ***************\n";
        std::cout << mt_source1.source_name << "\t" << mt_source1.period_hours << " h\n";
        std::cout << mt_source2.source_name << "\t" << mt_source2.period_hours << " h\n";
        std::cout << mt_source3.source_name << "\t" << mt_source3.period_hours << " h\n";
        std::cout << "Polynomial order: " << polynomial_order << "\n\n";
        
        mkdir("Forward_solutions", S_IRWXU);
    }

    // ------------------------------------------------------------------------
    // Step 2: Read and partition mesh
    // ------------------------------------------------------------------------
    double mesh_start_time = MPI_Wtime();
    
    if (process_rank == 0)
    {
        std::cout << "Reading mesh...\n";
    }
    
    auto* serial_mesh = new Mesh(model_parameters.mesh_file.c_str(), 1, 1);

    if (serial_mesh->GetNodes() == nullptr)
    {
        serial_mesh->EnsureNodes();
    }
    
    if (model_parameters.marker_type == "element_id")
    {
        for (int idx = 0; idx < serial_mesh->GetNE(); ++idx)
        {
            serial_mesh->SetAttribute(idx, model_parameters.marker_list[idx]);
        }
        serial_mesh->SetAttributes();
    }

    // Create parallel meshes (one per source)
    auto* parallel_mesh1 = new ParMesh(MPI_COMM_WORLD, *serial_mesh);
    auto* parallel_mesh2 = new ParMesh(MPI_COMM_WORLD, *serial_mesh);
    auto* parallel_mesh3 = new ParMesh(MPI_COMM_WORLD, *serial_mesh);

    delete serial_mesh;

    // Save mesh statistics
    std::ostringstream stats_filename;
    stats_filename << "Forward_solutions/FEM_mesh.stats";
    std::ofstream stats_file(stats_filename.str());
    stats_file.precision(8);

    parallel_mesh1->PrintInfo(stats_file);
    parallel_mesh2->PrintInfo(stats_file);
    parallel_mesh3->PrintInfo(stats_file);

    double mesh_end_time = MPI_Wtime();
    
    if (process_rank == 0)
    {
        std::cout << "Mesh loading time: " << mesh_end_time - mesh_start_time << " (s)\n\n";
    }

    // ------------------------------------------------------------------------
    // Step 3: Run finite element solvers
    // ------------------------------------------------------------------------
    if (polynomial_order > 1)
    {
        parallel_mesh1->ReorientTetMesh();
        parallel_mesh2->ReorientTetMesh();
        parallel_mesh3->ReorientTetMesh();
    }

    // Source 1
    FiniteElementSolver fem_solver1(model_parameters, mt_source1, *parallel_mesh1);
    fem_solver1.initialize(polynomial_order);
    fem_solver1.run_forward_modeling(1);

    // Source 2
    FiniteElementSolver fem_solver2(model_parameters, mt_source2, *parallel_mesh2);
    fem_solver2.initialize(polynomial_order);
    fem_solver2.run_forward_modeling(2);

    // Source 3
    FiniteElementSolver fem_solver3(model_parameters, mt_source3, *parallel_mesh3);
    fem_solver3.initialize(polynomial_order);
    fem_solver3.run_forward_modeling(3);

    // ------------------------------------------------------------------------
    // Cleanup and finalize
    // ------------------------------------------------------------------------
    delete parallel_mesh1;
    delete parallel_mesh2;
    delete parallel_mesh3;

    double total_end_time = MPI_Wtime();
    
    if (process_rank == 0)
    {
        std::cout << "\nTotal execution time: " << total_end_time - total_start_time << " (s)\n";
    }

    MPI_Finalize();
    return 0;
}
