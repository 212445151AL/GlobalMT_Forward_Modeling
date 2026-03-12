# Global MT Forward Modeling(Considering anisotropic media)

A parallel finite element solver for global MT forward modeling with adaptive mesh refinement.

## Features

- **3D MT Modeling**: Solves the curl-curl equation for secondary electric fields
- **Adaptive Mesh Refinement**: Goal-oriented and non-goal-oriented h-refinement strategies
- **Parallel Computing**: MPI-based distributed memory parallelism
- **Multiple Source Periods**: Handles three independent motional current sources
- **Flexible Solver Options**: FGMRES with AMS or PCG-AMS preconditioning
- **Station Interpolation**: GSLIB-based point location for field extraction at observation points
- **VTK Output**: Visualization-ready output with error indicators and conductivity tensors

## Dependencies

- **MFEM** (≥ 4.5): Finite element library
- **Hypre** (≥ 2.24.0): Linear solvers and preconditioners
- **METIS** (≥ 5.1.0): Mesh partitioning
- **GSLIB** (optional): Enhanced point location (Nek5000/gslib)
- **MPI**: Message Passing Interface (e.g., OpenMPI, MPICH)

## Building

### CMake Build (Recommended)

```bash
git clone https://github.com/your-repo/motional-mt.git
cd motional-mt
mkdir build && cd build
cmake .. -DCMAKE_PREFIX_PATH="/path/to/mfem;/path/to/hypre"
make -j4
```

### Manual Build

```bash
mpicxx -std=c++11 -O3 \
  -I/path/to/mfem -I/path/to/hypre/include \
  -o GlobalMT \
  GlobalMT.cpp MTCurrentSource.cpp parameterHandler.cpp \
  postProcessor.cpp solver.cpp errorEstimators.cpp fem.cpp \
  -L/path/to/mfem/lib -lmfem -L/path/to/hypre/lib -lHYPRE \
  -L/path/to/metis/lib -lmetis
```

## Input Files

### 1. Main Parameter File

```
source1.txt          # Source file 1
source2.txt          # Source file 2
source3.txt          # Source file 3
1                    # Polynomial order (currently fixed to 1)
solver.opts          # Solver parameters file
conductivity.mod     # Conductivity model file
mesh.msh             # Mesh file (MFEM format)
stations.dat         # Station coordinates file
1                    # Boundary marker for Dirichlet BC
0                    # AMR method: 0=goal-oriented, 1=non-goal, 2=global
5                    # Maximum refinement iterations
5000000              # Maximum DOFs
0.3                  # Refinement threshold ratio
1                    # Print VTK files (0/1)
```

### 2. Motional Current Source File

```
# Format: source_name period_hours
Ocean_P1_S1 12.0
# element_count
3
# element_id Jx_real Jx_imag Jy_real Jy_imag Jz_real Jz_imag (You need to consider the coordinate system transformation)
1001 1.2e-3 0.5e-3 2.1e-3 0.3e-3 0.0 0.0
1002 0.8e-3 0.2e-3 1.5e-3 0.1e-3 0.0 0.0
1003 0.5e-3 0.0   0.9e-3 0.0   0.0 0.0
```

### 3. Conductivity Model File

```
# marker_type (region_marker or element_id)
region_marker
# num_regions
3
# marker sigma_xx sigma_xy sigma_xz sigma_yx sigma_yy sigma_yz sigma_zx sigma_zy sigma_zz tag
1 3.3e6 0.0 0.0 0.0 3.3e6 0.0 0.0 0.0 3.3e6 1.0
2 1.0e-2 0.0 0.0 0.0 1.0e-2 0.0 0.0 0.0 1.0e-2 0.0
3 1.0e-4 0.0 0.0 0.0 1.0e-4 0.0 0.0 0.0 1.0e-4 -1.0
```

### 4. Station Coordinates File

```
# station_count
5
# r(km) theta(deg) phi(deg)
6371.0 90.0 0.0
6371.0 90.0 90.0
6371.0 45.0 0.0
6371.0 45.0 90.0
6371.0 0.0 0.0
```

### 5. Solver Parameters File

```
# FGMRES parameters
fgmres_max_iter 500
fgmres_primal_tol 1e-6
fgmres_dual_tol 1e-4
fgmres_restart 50
fgmres_print 2

# Preconditioner type: 0=AMS, 1=PCG-AMS, 2=Multigrid
prec_type 0

# PCG parameters (used only with prec_type=1)
pcg_max_iter 100
pcg_tol 1e-2
pcg_print 0

# AMS parameters
ams_cycle 1          # 1=multiplicative, 2=additive
ams_max_iter 20
ams_tol 1e-2
ams_print 1

# Smoothing options
a_relax_type 3       # 0=none, 1=Jacobi, 2=Gauss-Seidel, 3=hybrid
a_relax_sweeps 2
a_relax_weight 1.0
a_relax_omega 1.0

# AMG options (for both B_Pi and B_G)
amg_coarsen_type 10
amg_agg_levels 0
amg_strength 0.25
amg_interp_type 6
amg_max_row_elem 3
amg_relax_type 8
```

## Running

```bash
mpirun -np 4 ./GlobalMT input.txt
```

### Output Files

All output files are created in the `Forward_solutions/` directory:

- `MT_source1_iter0.sol`, `MT_source2_iter0.sol`, `MT_source3_iter0.sol`: Field values at stations (nT for B, V/m for E)
- `period_amr0_p1.vtk`: VTK visualization files for each iteration
- `FEM_mesh.stats`: Mesh statistics
- `input_parameter_list.log`: Echo of input parameters
- `input_solver_parm.log`: Echo of solver parameters

### Station Output Format

Each line contains (tab-separated):

| Column | Field | Units |
|--------|-------|-------|
| 1 | Radius | km |
| 2 | Theta | degrees |
| 3 | Phi | degrees |
| 4-5 | B_r (real, imag) | nT |
| 6-7 | B_theta (real, imag) | nT |
| 8-9 | B_phi (real, imag) | nT |
| 10-11 | E_r (real, imag) | V/m |
| 12-13 | E_theta (real, imag) | V/m |
| 14-15 | E_phi (real, imag) | V/m |

### Error Estimation

The face-jump error estimator (Ren et al., 2013) is based on the discontinuity of the normal component of the current density:
```
η_K = ∫_{∂K} ⟦n·(σE_s)⟧² dS
```

## Developer Information

**Authors**: Liangyu Xie, Chaojian Chen, Juanyu Wang  
**Institute**: Central South University (CSU), China  
**Email**: 8211221219@csu.edu.cn  
**GitHub**: [212445151AL](https://github.com/212445151AL)
