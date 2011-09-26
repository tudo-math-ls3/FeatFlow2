Directory structure
-------------------
(Relative to the main directory)

./
--> User provided routines.
    - indat2d.f
    - parq2d.f
    - fictbdry.f
    - fictbdrypc.f

./include
--> Include files of global importance; Memory management, 
    Parametrisation, Constants

./src/core
--> The main application(s).

./src/init
--> General initialisation routines for the main program.
    Read parameters from DAT files, initialise triangulations.

./src/optimization
--> Optimisation components; not used in default project

./src/nonstationary
--> The nonstationary Navier-Stokes solver NONST2

./src/stationary
--> The stationary Navier-Stokes solver NSDEF2

./src/solvers
--> Support for solvers with extended calling convention.
    A collection of black-box solvers/smoothers  (Multigrid, BiCGStab)

./src/discr_nonlinear
--> High-level routines for the discretisation process of the
    nonlinear equations / global system.
    Assembly of nonlinear matrices, global matrix,
    initialization/generation of matrices/vectors on all levels

./src/gridadapt
--> Grid adaption support

./src/discr_linear
--> Low-level routines for the discretisation process of the linear
    equations.
    Element routines, Matrix-vector multiplication, implementation
    of boundary conditions, elementary multigrid routines.

./src/parametrization
--> Support for 
    - 2D OMEGA and FEAT parametrisations
    - general fictitious boundaries

./src/triangulation
--> Support for extended 2D triangulations, scalar evaluation
    of FEM functions, calculation of cutlines, tracers, GMV-output

./src/geometry
--> Support for fictitious boundary geometries

./src/postprocessing
--> Low level pre- and postprocessing routines
    (Reading/Writing solutions from/to  disc, particle tracing,
     Calculation of integrals and boundary forces, Calculation of 
     (a-posteriori) error)
     
./src/misc
--> General low-level routines, not connected to anything
    (Support for strings, file-output, reading from DAT-files,
     Heap-Sort,...)
    
