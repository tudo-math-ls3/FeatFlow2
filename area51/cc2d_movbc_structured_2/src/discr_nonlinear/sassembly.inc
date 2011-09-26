************************************************************************
* This file describes a structure that collects all information
* regarding the Finite Element discretization of the linear and
* nonlinear terms and right-hand-side vectors.
*
* The structures themselves are realized as integer/double array
* combination:
*   IASMBL = array [1..SZASMI] of integer 
*          = TIntAssembly
*   DASMBL = array [1..SZASMD] of integer 
*          = TDoubleAssembly
* where the entries can be accessed with the offsets defined in this
* file. In detail the structures contain the following data:
*
* TIntAssembly = record
*
*   // Problem related information - describes what to discretize
*
*   ISTOK  : Whether to discretize Stokes or Navier Stokes
*            = 0: Discretize Navier Stokes equation
*            = 1: Discretize Stokes equation
*
*   // Discretization related information - describes how to discretize
*   
*   IELEMT : Element identifier to prescribe the element that should be
*            used for the discretization of the system matrix.
*            = 0: E030
*            = 1: E031
*            = 2: EM30 (standard)
*            = 3: EM31
*   IAPRM  : configure adaptive matrix generation. Bitfield.
*            0   =no adaptive matrix generation = standard
*            Bit2=switch to constant mat. gen., depending on DMTEP
*            Bit5(additionally to Bit2)
*                =switch mat. gen. depenting on size of neighbour 
*                 element, too
*   ICUBM  : Cubature formula for building the real mass matric
*            Standard = 4
*   ICUBA  : Cubature formula for building the Stokes (Laplace) matrix
*            Standard = 4
*   ICUBN  : Cubature formula for building the convective matrix
*            Standard = 4
*   ICUBB  : Cubature formula for building the Gradient matrix B1,B2
*            Standard = 4
*   ICUBF  : Cubature formula for building the right hand side
*            Standard = 4
*   IMASS  : = 0: use lumped mass matrix = standard
*            = 1: use real mass matrix
*   IMASSL : If IMASS=0:
*            = 0: Use standard mass lumping
*            = 1: Use diagonal lumping (add non-diagonal entries to
*                 the diagonal)
*   IUPW   : =0: Use Upwinding technique to build nonlinearity; 
*                standard
*            =1: Use streamline diffusion to build nonlinearity
*   IPRECA : Method how to set up the system matrix
*            =0,
*            =1: Precalculation of Stokes (=Laplace)-matrix.
*                Nonlinear matrix is assembled using the precalculated
*                matrix.
*            =2,
*            =3: Precalculation of Stokes matrix, stored on disc.
*                When setting up the nonlinear matrix, the linear
*                Stokes matrix is read into memory and used for
*                assembling
*            =4: Nonlinear system matrix is completely rebuild in
*                every iteration without being assembled separately
*                by a linear and nonlinear part (not implemented).
*            Standard=1
*   IPRECB : Method how to set up pressure coupling matrices.
*            =0,
*            =1: Calculation by exact quadrature rule 
*                Works only with Q1~/Q0
*            =2: (undefined)
*            =3,
*            =4: Calculation with complete evaluation 
*                (slower, but works for all FE pairs)
*            Standard=1
*
*   // Information about the discretization on the boundary
*
*   IINFCF : Type of inflow profile.
*            Used in INDAT2D.F to choose the type of the inflow.
*            Setting for standard INDAT2D:
*            =0: benchmark configuration, parabolic inflow from 
*                the left; 
*            =1: vertical column, parabolic inflow from the bottom
*            Can be user defined if INDAT2D is changed!
*
*   // Dynamic information
*   // Information that is dynamically calculated and updated during
*   // the solution process, depending on the current situation.
*   // The user has not to maintain information in these variables.
*   // During the solution process, discretization routines might
*   // access the data stored here to do the discretization correctly.
*
*   INEUM  : Indicates whether the assembling of the matrices need
*            to handle other types of boundaries than Dirichlet.
*            =0: The problem contains only Dirichlet boundary data;
*                standard
*            =1: The problem contains Dirichlet and Neumann boundaries
*            Updated in the initialization and before each time step
*
*   IALPHA : The linear part of the nonlinear iteration matrix 
*            has the form
*                    [ alpha*M + THETA*L ] u
*            IALPHA=ALPHA=0/1 decides, whether the mass matrix is 
*            added to the iteration matrix.
*            =0: mass matrix is not added (typical for stationary
*                simulation, standard)
*            =1: mass matrix is added (typical for nonstationary
*                simulation, e.g. if NSDEF is a sub-solver of a 
*                nonstationary solver)
*
*   IALE   : The convective part of the nonlinear iteration has the
*            form
*                    U * grad(U)
*            If IALE=1, this term is modified to subtract a "grid
*            velocity" in case that the computational grid is moving.
*            The convective part is then discretized in the form
*                    U * grad(U-V)
*            with 
*                    V = array [1..2,1..NVT] of double
*            being an approximation to the velocity of the corner
*            vertices of the grid (X- and Y-velocity).
*            This velocity field is interpreted by the
*            discretization routines to fit to the degrees of freedom
*            of the current discretization.
*
*   LGRVEL : If IALE=1, during the solution process this variable
*            receives a handle to an array [1..2,1..NVT] of double,
*            which specifies the grid velocity in every corner
*            vertex of the grid.
*            The memory is released when the solver is finished.
*            If IALE=0, this variable is not used.
*
* end
*
* TDoubleAssembly = record
*
*   // Problem related information - describes what to discretize
*
*   RE     : Reynolds-number 1/NU for the equation
*   NY     : Must be set to NU=1/RE; helps to skip some divisions.
*
*   // Discretization related information   
*
*   UPSAM  : Parameter to control the creation of the nonlinear terms.
*            IUPW=0: UPSAM controls upwinding.
*                    < 0: simple upwinding
*                    >=0: Smarskji-Upwinding
*            IUPW=1: UPSAM controls streamline diffusion
*   DMTEP  : treshold parameter: switch construction of coarse grid 
*            matrices from standard finite element approach to locally
*            constant interpolation for all rows belonging to elements
*            with aspect ratio >= DMTEPS.
*            0D0=all elements
*            20D0=standard
*            -1D0=infinity=no matrix modifications
*
*   // Information about the discretization on the boundary
*
*   DPUMAX : Inflow velocity.
*            Used in INDAT2D.F to set the inflow velocity.
*            Setting for standard INDAT2D:
*            =0.3 : benchmark configuration, horizontal channel
*            =10.0: vertical channel
*            Can be user defined if INDAT2D is changed!
*
*   // Dynamic information
*   // Maintained by the solver, depending on the solver.
*   // The user has not to maintain information in these variables.
*   // During the solution process, discretization routines might
*   // access the data stored here to do the discretization correctly.
*
*   TIMENS : Current simulation time.
*            Not used for stationary simulation, =0 there.
*            Updated by the nonlinear solver for the assembling process.
*   THWEIG : double. Theta-scheme weight.
*            For the creation of the system matrix, this value defines
*            the weight to multiply the Stokes/nonlinear submatrix with.
*            For stationary simulations this parameter must be set
*            to 1D0 to include the full nonlinear matrix into the
*            system matrix:
*                [ IALPHA * M  +  THWEIG * N(u) ] u + ... = ...
* end
************************************************************************

      INTEGER OIELEMT,OIAPRM ,OICUBM ,OICUBA ,OICUBN ,OICUBB ,OICUBF 
      INTEGER OIMASS ,OIMASSL,OIUPW  ,OIPRECA,OIPRECB,OISTOK, OINEUM
      INTEGER OIALPHA,OIALE  ,OLGRVEL,OIINFCF
      
      PARAMETER (OISTOK  =  1)
      PARAMETER (OIELEMT =  2)
      PARAMETER (OIAPRM  =  3)
      PARAMETER (OICUBM  =  4)
      PARAMETER (OICUBA  =  5)
      PARAMETER (OICUBN  =  6)
      PARAMETER (OICUBB  =  7)
      PARAMETER (OICUBF  =  8)
      PARAMETER (OIMASS  =  9)
      PARAMETER (OIMASSL = 10)
      PARAMETER (OIUPW   = 11)
      PARAMETER (OIPRECA = 12)
      PARAMETER (OIPRECB = 13)
      PARAMETER (OIINFCF = 14)
      
      PARAMETER (OINEUM  = 15)
      PARAMETER (OIALPHA = 16)
      PARAMETER (OIALE   = 17)
      PARAMETER (OLGRVEL = 18)
      
      INTEGER SZASMI
      
      PARAMETER (SZASMI  = 20)


      INTEGER ORE,ONY,OUPSAM,ODMTEP,OTIMENS,OTHWEIG,ODPUMAX
      
      PARAMETER (ORE    =  1)
      PARAMETER (ONY    =  2)
      
      PARAMETER (OUPSAM =  3)
      PARAMETER (ODMTEP =  4)
      
      PARAMETER (ODPUMAX = 5)
      
      PARAMETER (OTIMENS=  6)
      PARAMETER (OTHWEIG=  7)
      
      INTEGER SZASMD
      
      PARAMETER (SZASMD =  16)
      