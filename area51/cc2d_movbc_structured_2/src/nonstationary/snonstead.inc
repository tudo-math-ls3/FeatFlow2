************************************************************************
* This file describes the solver structure for the nonlinear
* iteration for solving the time dependend Navier Stokes equation.
* As the structure is an extension to the general solver structure,
* the file "SSOLVERS.INC" has to be included before.
*
* It defines two parameter blocks for integer and double parameters,
* respectively. Both parameter blocks contain three parts:
* - An input-block which has to be initialized by the caller
* - An output-block which returns the result of the computation
* - A status- & temporary block, which is only valid
*   during the computation and which gives additional information
*   to the callback routines, like simulation time, current
*   triangulations,... Therefore, the include file STRIA.INC
*   must be included prior to this file.
* The caller can attach additional elements to the parameter blocks
* with further information for the callback-routines, e.g. matrix-
* handles, flags, tags,...
*
* The structures themselves are realized as integer/double precision
* arrays:
*   IPARAM = array [1..SZISDI] of integer
*          = TInstatDefIParams
*   DPARAM = array [1..SZISDD] of double
*          = TInstatDefDParams
*
* These structures are direct extensions to the original solver
* structure, i.e. there has no SZxxxx-constant to be added to
* the offset variable, the variables can be acessed directly
* by a construction like IPARAM(OITE).
* The structures are organized as follows. The first part of the
* structure is as usual a block with general solver parameters.
* Some of the general parameters are not used as indicated with "=0"
* below.
*
* The whole structure can be initializes with the INIISD subroutine
* in NONSTEAD2.F.
*
* TInstatDefIParams = record
*   // The first part of TInstatDefIParams is the whole TSolverIParams-
*   // structure, i.e. TInstatDefIParams is a direct extension of
*   // TSolverIParams:
*
*   // Standard-Output:
*
*   STATUS : Result of the solving process.
*            =0: Solver converged successfully
*            =1: Nonlinear solver broke down
*            =2: Linear solver broke down
*   ITE    : Number of performed iterations without repetitions
*            of time steps
*
*   // Standard-Input:
*
*   SLTAG  : Solver tag.
*            =21: Nonstationary Navier-Stokes.
*   NITMIN : = 0
*   NITMAX : maximum number of iterations = time steps
*   IFILT  : = 0
*   OPCTAG : = 0
*   INRM   : = 0
*   ITIM   : = 0
*   IASRLN : = 0
*   MSGTRM : Level of output of the solver of algorithm; standard = 2
*   MFLTRM : = 0
*   MTRMRS : = 1
*   MATI   : = 0
*   PRECI  : = 0
*
*   // Standard-Status
*
*   CRITE  : Current number of the iteration / time step
*   NEQ    : Number of equations in the solution/RHS vectors
*            (on finest level)
*
*   // The second part now describes the configuration of the
*   // nonlinear solver specific settings. It again consists of an
*   // Output-block, an Input-block and a Status-block.
*
*   // Nonstationary solver, Output:
*
*   TITTS  : Total number of time steps including all repetitions
*   TITNL  : Total number of iterations of nonlinear solver
*   TITLI  : Total number of iterations of linear solver
*
*   // Nonstationary solver, Input:
*
*   IRHS   : = 0: homogenuous RHS
*            = 1: steady inhomogenuous RHS
*            = 2: time-dependent RHS
*   IBDR   : Type of boundary conditions. Bitfield.
*            =0: Dirichlet + not moving homogenuous Neumann
*                boundary conditions
*            Bit 0: activate inhomogenuous Neumann boundary
*                   (=pressure drop)
*            Bit 1: activate time-dependent location of
*                   Dirichlet/Neumann boundary components,
*                   time-dependent fictitious boundaries
*   IDMESH : Type of dynamic mesh construction
*            =0: Fixed mesh in all time steps
*            =1: Mesh adaption in the very first time step to generate
*                basic mesh; kept fix during the whole simulation.
*            =2: Moving mesh, deformed in every time step, starting from
*                the initial mesh
*            =3: Moving mesh, deformed in every time step, starting from
*                the previous time step
*            =4: Dynamic mesh, rebuild in every time step
*                (not implemented)
*   IMSCRT  : mesh adaption criterion if IDMESH>=2
*             =0: geometric
*             =1: H1-error
*   IMSREP  : Number of repetitions of every time step if IDMESH>=2.
*             The mesh is readapted in every step depending on the
*             criterion IMSCRT
*   IMSSTP  : Number of time steps between each grid adaption:
*             The grid is adapted each IMSSTP steps.
*             IMSSTP=0 adapts the mesh only in the very first iteration,
*             even if IDMESH=2,3.
*
*   // Time-stepping input parameters:
*
*   IFRSTP : Type of time discretization to use.
*            =0: Use One-step scheme, depending on SSTHETA
*            =1: Use standard Fractional step, don't use SSTHETA
*
*   // Nonstationary solver, Status:
*
* end
*
* TInstatDefDParams = record
*   // The first part of TInstatDefDParams is the whole TSolverDParams-
*   // structure, i.e. TInstatDefDParams is a direct extension of
*   // TSolverDParams:
*
*   // Standard-Output:
*
*   DEFINI : = 0
*   DEFFIN : = 0
*   RHO    : = 0
*   RHOASM : = 0
*
*   TMTOT  : total time for nonlinear solver
*   TMFILT : = 0
*
*   // Standard-Input:
*
*   OMEGA  : = 0
*   OMGPC  : = 0
*
*   // The standard stopping criteria of the solver structure are not
*   // used! Instead, the structure defines a couple of "specialized"
*   // stopping criteria for the nonlinear solver, which cannot be
*   // acquires by the standard stopping criteria - see below.
*
*   EPSREL : = 0
*   EPSABS : = 0
*   DIVREL : = 0
*   DIVABS : = 0
*
*   VECZER : = 0
*   MATD   : = 0
*   PRECD  : = 0
*
*   // Standard-Status;
*   // nothing defined here.
*
*   // -----------------------------------------------------------------
*
*   // The second part now describes the configuration of the
*   // nonlinear solver specific settings. It again consists of an
*   // Output-block, an Input-block and a Status-block.
*
*   // Nonstationary solver, Output:
*
*   TNSTIM : array [1..SZTIMG] of double
*            TTimings-structure that collects statistical information
*            about needed time for the calculation.
*   TTGRID : Time for processing of the grid (mesh generation,...)
*
*   // Nonstationary solver, Input:
*
*   TIMEMX : Maximum simulation time; simulation is stopped upon
*            reaching time >= TIMEMX.
*
*   SSTHETA: Theta parameter for One-step (single step) Theta scheme.
*            Only used if IFRSTP=0.
*
*   DMSMTM : If mesh adaption is activated:
*            Maximum simulational time for mesh adaption. If simulational
*            time reaches DMSMTM, the mesh is kept fix.
*            =-1D0: no maximum time
*
*   // Nonstationary solver, Status:
*
*   REINST : Set to RE, based on the RE parameter of the stationary
*            solver.
*   NYINST : Set to NU=1/RE, based of NY parameter of the stationary
*            solver.
*
*   TIMEST : Initial simulation time
*
* end
************************************************************************

      INTEGER OIRHS,OIBDR,OIDMESH,OTITTS,OTITNL,OTITLI,OFRSTP
      INTEGER OIMSCRT,OIMSREP,OIMSSTP

      INTEGER SZISDI
      
      PARAMETER (OTITTS  = SZSLVI +  1)
      PARAMETER (OTITNL  = SZSLVI +  2)
      PARAMETER (OTITLI  = SZSLVI +  3)
      PARAMETER (OIRHS   = SZSLVI +  4)
      PARAMETER (OIBDR   = SZSLVI +  5)
      PARAMETER (OIDMESH = SZSLVI +  6)
      PARAMETER (OIMSCRT = SZSLVI +  7)
      PARAMETER (OIMSREP = SZSLVI +  8)
      PARAMETER (OIMSSTP = SZSLVI +  9)
      
      PARAMETER (OFRSTP  = SZSLVI + 10)
      
      PARAMETER (SZISDI  = SZSLVI + 16)

      INTEGER OTIMEST,OTNSTIM,OTTGRID
      INTEGER OTIMEMX,OSSTHETA,OREINST,ONYINST,ODMSMTM
      
      INTEGER SZISDD
      
      PARAMETER (OTIMEST = SZSLVD  +  1)
      PARAMETER (OTIMEMX = SZSLVD  +  2)
      PARAMETER (OSSTHETA= SZSLVD  +  3)
      PARAMETER (ODMSMTM = SZSLVD  +  4)
      
      PARAMETER (OTNSTIM = SZSLVD  +  5)
      PARAMETER (OTTGRID = SZSLVD  +  SZTIMG + 5)
      PARAMETER (OREINST = SZSLVD  +  SZTIMG + 6)
      PARAMETER (ONYINST = SZSLVD  +  SZTIMG + 7)
      
      PARAMETER (SZISDD  = SZSLVD  +  SZTIMG + 16)
      