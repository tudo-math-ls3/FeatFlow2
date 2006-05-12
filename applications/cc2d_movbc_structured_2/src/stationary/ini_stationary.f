***********************************************************************
* This file contains initialization routines for the problem of
* solving stationary Navier Stokes equations.
*
* The routines here initialise/generate structure information about
* the discretiation of the problem which is to solve.
***********************************************************************

***********************************************************************
* Initialize discretisation structures for stationary solver
*
* This routine allocates memory for the structures used in the
* stationary solver. It initializes all the structures, but does
* not compute anything.
*
* In:
*   MFILE  : handle to an file where output is written to
*   NLMIN  : minimum level where calculation is allowed
*   NLMAX  : maximum level where calculation is allowed
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV in its ummodified form, i.e.
*            without any boundary conditions implemented, probably
*            undeformed (if the user wants to start the grid
*            deformation in every time step from the original grids).
*
* Out:
*   MATDAT : array [1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialized with data
*   VECDAT : array [1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*   IMGPAR : array [1..SZ020I+2*SZSLVI] of integer
*   DMGPAR : array [1..SZ020D+2*SZSLVD] of double
*            Integer and double parameter blocks for the linear
*            sub-solver. 
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. 
*
*   LUP    : Handle to array [1..NEQ] of double
*            Iteration vector; not initialized
*   LRHS   : Handle to array [1..NEQ] of double
*            Right hand side vector
*   NEQ    : Length of solution/RHS vector
*
* INMAIN must be called prior to this routine.
* Parametrization and triangulations must be initialized prior to
* calling this routine.
***********************************************************************
      
      SUBROUTINE INSTSL (NLMIN,NLMAX,
     *                   TRIAS,
     *                   MATDAT,VECDAT,
     *                   ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,
     *                   IASMBL,DASMBL,
     *                   LUP,LRHS,NEQ)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      
      INCLUDE 'stiming.inc'
      INCLUDE 'snsdef.inc'
      
      INCLUDE 'dstrings.inc'
      
      INCLUDE 'cpostproc.inc'
      INCLUDE 'clinsol_cc2d.inc'
      
C     parameters

      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,*)
      INTEGER MATDAT(SZN2MI,*),VECDAT(SZN2VI,*),ISTPAR(*)
      INTEGER IMGPAR(*),IASMBL(*)
      INTEGER LUP,LRHS,NEQ
      DOUBLE PRECISION DSTPAR(*),DMGPAR(*),DASMBL(*)
      
C     Initialize assembly parameters

      CALL INIASM (IASMBL,DASMBL,1)

C     Initialize matrix parameter block

      CALL ININSM (NLMIN,NLMAX,TRIAS,IASMBL,DASMBL,MATDAT)

C     Initialize vector parameter block

      CALL ININSV (NLMIN,NLMAX,TRIAS,IASMBL,DASMBL,VECDAT)
      
C     Initialize parameter blocks of linear solver.
C     Maybe for multigrid or other linear solver.
C     The global variable ILINSL decides on the type of the solver.

      CALL INLSST (NLMIN,NLMAX,IMGPAR,DMGPAR,ILINSL,1)

C     Initialize NSDEF-structures

      CALL INISTS (ISTPAR,DSTPAR,1)
      
C     Initialize basic filename for autosave

      ISTPAR (OLFLAUT)  = STNEWC (.TRUE.,CFLAUT)
        
C     Get the dimension from the vector parameter block

      NEQ  = VECDAT(ONEQV,NLMAX)
      
C     Allocate memory for solution and RHS vector

      CALL ZNEW (NEQ,1,LUP,'DUP   ')
      CALL ZNEW (NEQ,1,LRHS,'RHS   ')
      
C     Remark:
C     The RHS vector is filled with 0. In particular, the pressure
C     part is filled with 0. During the whole computation, the pressure
C     part is normally not touched and will therefore always stay at 0!
      
      END
      
***********************************************************************
* Clean up discretisation structures for stationary solver
*
* This routine cleans up the structures for the nonstationary solver,
* releasing all allocated memory.
*
*   MFILE  : handle to an file where output is written to
*   NLMIN  : minimum level where calculation is allowed
*   NLMAX  : maximum level where calculation is allowed
*
*   MATDAT : array [1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialized with data
*   VECDAT : array [1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*
*   LUP    : Handle to array [1..NUVP] of double
*            Iteration vector; not initialized
*   LRHS   : Handle to array [1..NUVP] of double
*            Right hand side vector
*
* Out:
*   All handles to vectors / matrices / ... in the above structures
*   are invalid after calling this routine.
***********************************************************************
      
      SUBROUTINE DNSTSL (NLMIN,NLMAX,
     *                   MATDAT,VECDAT,
     *                   ISTPAR,DSTPAR,
     *                   LUP,LRHS)
      
      IMPLICIT NONE

      INCLUDE 'stria.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'

      INCLUDE 'stiming.inc'
      INCLUDE 'snsdef.inc'

C     parameters

      INTEGER NLMIN,NLMAX
      INTEGER MATDAT(SZN2MI,*),VECDAT(SZN2VI,*),ISTPAR(*)
      INTEGER LUP,LRHS
      DOUBLE PRECISION DSTPAR(*)
      
C     Release memory of solution/RHS vectors

      CALL ZDISP (0,LRHS,'RHS   ')
      CALL ZDISP (0,LUP,'DUP   ')
      
C     Release matrix parameter block

      CALL DONNSM (NLMIN,NLMAX,MATDAT)

C     Release vector parameter block

      CALL DONNSV (NLMIN,NLMAX,VECDAT)
      
C     Release basic filename for autosave

      IF (ISTPAR (OLFLAUT).NE.0) CALL STDIS (ISTPAR (OLFLAUT))
      
      END
      
************************************************************************
* Generate geometry and boundary information for stationary solver
*
* This routine sets up boundary releated structures for the solver.
* Calculates the shortcut nodal property array, modifies KNXPR,...
* If necessary arrays do not exist, they are allocated. Existing
* arrays are updated.
*
* As this preparation phase already needs information about the problem
* and the discretization, the routine has major access to the parameters
* of the solver, which have to be initialized prior to calling this
* routine by the appropriate initialization routine!
*
* The routine uses the given triangulation for the discretization.
* If mesh adaption should be performed, it must be done prior to
* this routine.
*
* For nonstationary simulation, this routine is called by the solver
* automatically. As time index, TIMENS from the DASMBL data
* structures is used.
*
* In:
*   MFILE  : handle to an file where output is written to
*   NLMIN  : minimum level where calculation is allowed
*   NLMAX  : maximum level where calculation is allowed
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV in its ummodified form, i.e.
*            without any boundary conditions implemented, probably
*            undeformed (if the user wants to start the grid
*            deformation in every time step from the original grids).
*
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. 
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to fictitious boundary
*            routines. Not used in this routine.
*
* Out:
*   - Generates precalculated information about the 
*     current geometries, based on the current triangulations. A link to
*     the precalculated information is saved in each triangulation
*     structure. Any previously generated precalculated information
*     is destroyed. Memory is allocated if necessary.
*   - Implements boundary information into the KXNPR-array.
*     Figures out whether there are Neumann boundary components
*     on the boundary in the current situation.
*   - Allocates (if necessary) and calculates a shortcut nodal property 
*     depending on the geometry and triangulation (see bdry2.f). Any 
*     previously existing array containing shortcut nodal properties is
*     destroyed and rebuild.
*   - Updates the assembly data structure, whether there are Neumann
*     boundary components or not (-> INEUM).
************************************************************************

      SUBROUTINE GNBDGE (NLMIN,NLMAX,
     *                   TRIAS,
     *                   IASMBL,DASMBL,
     *                   IGEOM,DGEOM)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'

      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'sassembly.inc'

C     parameters

      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,*),IGEOM(*)
      INTEGER IASMBL(*)
      DOUBLE PRECISION DASMBL(*),DGEOM(*)
      LOGICAL BDC

C     local variables

      INTEGER KXNPR,INEUM,I,KSCNPR
      
C     externals

      DOUBLE PRECISION PARX,PARY,TMAX,UE
      EXTERNAL PARX,PARY,TMAX,UE

C     -----------------------------------------------------------------
C     Step 1: Generate precalculated information from the geometries/
C             triangulations
C
C     Release old precalculated information if it exists:
      
      DO I=NLMAX,NLMIN,-1

C       If there is any precalculated information, delete it - we want
C       to recalculate

        IF (TRIAS(OTRIUD+2,I).NE.0) THEN
          CALL ZDISP(0,TRIAS(OTRIUD+2,I),'LDARR ')
        END IF

        IF (TRIAS(OTRIUD+3,I).NE.0) THEN
          CALL ZDISP(0,TRIAS(OTRIUD+3,I),'LIARR ')
        END IF

      END DO
      
C     Calculate precalculated information
      
      DO I=NLMIN,NLMAX

C       Precalculate information. This returns handles to int/double-
C       arrays. These are saved directly into the TRIA structure into
C       the first two variables of the user-defined block.
        
        CALL FBDPRC (TRIAS(1,I),TRIAS(OTRIUD+2,I),TRIAS(OTRIUD+3,I),
     *               IGEOM,DGEOM)
      
      END DO

C     The precalculated information allowes us to perform faster
C     calculations when setting up matrices/vectors in a geometry
C     containing fictitious boundaries.

C     -----------------------------------------------------------------
C     Step 2: Implement boundary condition into KXNPR.
C             Calculate shortcut nodal property array.
C
      DO I=NLMIN,NLMAX

C       Allocate memory for the shortcut nodal property, if the memory
C       is not already allocated.

        IF (TRIAS(OTRIUD+1,I).EQ.0) THEN
          CALL ZNEW (TRIAS(ONMT,I),3,TRIAS(OTRIUD+1,I),'KSCNPR')
        END IF

C       Update the boundary conditions using the update routine and
C       parameter block of the stationary solver. Use the time
C       stamp in DASMBL to identify the current simulational time.
C       The routine returns INEUM, writes information about the
C       shortcut nodal property array KSCNPR and updates KXNPR 
C       to the current situation, based on the KXNPR array in TRIAS. 
C       We provide the KXNPR array of the TRIAS structure as output
C       here, too, which results in updating the KXNPR array.
      
        KXNPR = L(TRIAS(OLXNPR,I))
        KSCNPR = L(TRIAS(OTRIUD+1,I))
      
        CALL UPDBDX (TRIAS(1,I),DASMBL(OTIMENS),INEUM,KWORK(KXNPR),
     *               TRIAS(OTRIUD,I),KWORK(KSCNPR), 
     *               IASMBL,DASMBL,IGEOM,DGEOM)

C       We must not release unused memory of the shortcut nodal property
C       array!
C       Reason: Because of fictitious boundaries, due to modification of 
C       the objects in the domain, every edge may get a Dirichlet edge 
C       and thus being added to this array. So depending on the object
C       (what we don't know), the shortcut nodal property array might
C       receive up to all edges of the whole triangulation!

CCC     CALL ZDISP(TRIAS(OTRIUD,I),TRIAS(OTRIUD+1,I),'KSCNPR')
        
      END DO
      
C     INEUM on the finest level defines whether there are
C     Neumann boundary conditions in our problem at all.

      IASMBL(OINEUM) = INEUM

      END

************************************************************************
* Release boundary/geometry information
*
* This routine releases the memory for any type of information
* which was calculated in GNBDGE.
*
* In:
*   MFILE  : handle to an file where output is written to
*   NLMIN  : minimum level where calculation is allowed
*   NLMAX  : maximum level where calculation is allowed
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV in its ummodified form, i.e.
*            without any boundary conditions implemented, probably
*            undeformed (if the user wants to start the grid
*            deformation in every time step from the original grids).
*
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. 
************************************************************************

      SUBROUTINE DNBDGE (NLMIN,NLMAX,
     *                   TRIAS,
     *                   IASMBL,DASMBL)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'

      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'sassembly.inc'

C     parameters

      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,*)
      INTEGER IASMBL(*)
      DOUBLE PRECISION DASMBL(*)

C     local variables

      INTEGER I
      
C     externals

      DOUBLE PRECISION PARX,PARY,TMAX,UE
      EXTERNAL PARX,PARY,TMAX,UE

C     Release old precalculated geometry information if it exists:
      
      DO I=NLMAX,NLMIN,-1

        IF (TRIAS(OTRIUD+2,I).NE.0) THEN
          CALL ZDISP(0,TRIAS(OTRIUD+2,I),'LDARR ')
        END IF

        IF (TRIAS(OTRIUD+2,I).NE.0) THEN
          CALL ZDISP(0,TRIAS(OTRIUD+3,I),'LIARR ')
        END IF

      END DO
      
C     Remove any existing array containing the shortcut nodal property

      DO I = NLMAX,NLMIN,-1
        IF (TRIAS(OTRIUD+1,I).NE.0) THEN
          CALL ZDISP(0,TRIAS(OTRIUD+1,I),'KSCNPR')
        END IF
      END DO
      
      END

************************************************************************
* Generate (constant) matrix information for stationary solver
*
* This routine calculates all "missing" matrix information for the call
* to NSDEF. It adapts the grid (if necessary), generates a RHS-vector,
* and generates matrix entries. Afterwards the solver can be called
* to solve the problem.
*
* As this preparation phase already needs information about the problem
* and the discretization, the routine has major access to the parameters
* of the solver, which have to be initialized prior to calling this
* routine by the appropriate initialization routine!
*
* The routine uses the given triangulation for the discretization.
* If mesh adaption should be performed, it must be done prior to
* this routine.
*
* In:
*   MFILE  : handle to an file where output is written to
*   NLMIN  : minimum level where calculation is allowed
*   NLMAX  : maximum level where calculation is allowed
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV in its ummodified form, i.e.
*            without any boundary conditions implemented, probably
*            undeformed (if the user wants to start the grid
*            deformation in every time step from the original grids).
*   MATDAT : array [1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialized with data
*
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. 
*
*   NEQU   : length of velocity vector (number of equations in u)
*   NEQP   : length of pressure vector (numer of equations in p)
*
* Out:
*   - Allocates memory and generates precalculated information about the 
*     current geometries, based on the current triangulations. A link to
*     the precalculated information is saved in each triangulation
*     structure. Any previously generated precalculated information
*     is destroyed.
*   - Implements boundary information into the KXNPR-array.
*   - Allocates and calculates a shortcut nodal property depending on
*     the geometry and triangulation (see bdry2.f). Any previously
*     existing array containing the shortcut nodal property is
*     destroyed and rebuild.
*   - Missing matrices in MATDAT are calculated (Stokes matrix, mass
*     matrix)
************************************************************************

      SUBROUTINE GNMATD (NLMIN,NLMAX,
     *                   TRIAS,
     *                   MATDAT,
     *                   ISTPAR,DSTPAR,
     *                   IASMBL,DASMBL)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'

      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'stiming.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'snsdef.inc'

C     parameters

      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,*)
      INTEGER MATDAT(SZN2MI,*)
      INTEGER IASMBL(*),ISTPAR(*)
      DOUBLE PRECISION DSTPAR(*),DASMBL(*)

C     Generate matrix entries

      CALL GENNSM (NLMIN,NLMAX,TRIAS,IASMBL,DASMBL,255,MATDAT)

      END

************************************************************************
* Generate vector information for the stationary solver
*
* This routine calculates all "missing" information for the call
* to NSDEF. It generates a RHS-vector, and implements boundary
* conditions if desired.
*
* As this preparation phase already needs information about the problem,
* and the discretization, the discretization structure must be set up
* completely before calling this routine.
*
* The routine needs also precalculated information about boundaries,
* so the appropriate matrix/boundary generation routine has to be
* executed before.
*
* In:
*   MFILE  : handle to an file where output is written to
*   NLMIN  : minimum level where calculation is allowed
*   NLMAX  : maximum level where calculation is allowed
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV in its ummodified form, i.e.
*            without any boundary conditions implemented, probably
*            undeformed (if the user wants to start the grid
*            deformation in every time step from the original grids).
*   VECDAT : array [1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. 
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*
*   IBDTYP : Type of boundary conditions to implement into RHS / solution
*            vector. Bitfield.
*            =0: Don't implement boundary conditions;
*                standard for generating the first RHS for time-
*                dependent simulation.
*            Bit 0: whether to implement pressure drop into RHS-vectors
*            Bit 1: whether to implement dirichlet boundary conditions
*                   into RHS vectors
*            Bit 2: whether to implement dirichlet boundary conditions
*                   into solution vectors
*
* Out:
*   - The RHS vector in LRHS is calculated
*   - Implements boundary conditions into RHS/solution if desired
************************************************************************

      SUBROUTINE GNSTVC (NLMIN,NLMAX,
     *                   TRIAS,
     *                   VECDAT,
     *                   ISTPAR,DSTPAR,
     *                   IASMBL,DASMBL,IGEOM,DGEOM,
     *                   IBDTYP,
     *                   LUP,LRHS)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'

      INCLUDE 'stria.inc'
      
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      

C     parameters

      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,*),IGEOM(*)
      INTEGER VECDAT(SZN2VI,*),ISTPAR(*),IASMBL(*)
      INTEGER LUP,LRHS,IBDTYP
      DOUBLE PRECISION DSTPAR(*),DASMBL(*),DGEOM(*)
                             
C     externals
      
      DOUBLE PRECISION PARX,PARY,TMAX,UE
      EXTERNAL PARX,PARY,TMAX,UE
                             
C     local variables

      INTEGER KRHS,KU,NEQU,NEQP

C     Get the vector size of velocity and pressure part:

      NEQU = VECDAT(ONU,NLMAX)
      NEQP = VECDAT(ONP,NLMAX)

C     -----------------------------------------------------------------
C     Step 1: Generate RHS vector on finest level

      KRHS = L(LRHS)

      CALL GNRHSV (TRIAS(1,NLMAX),
     *             ISTPAR,DSTPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *             VECDAT(1,NLMAX),DWORK(KRHS))

C     -----------------------------------------------------------------
C     Step 2: Implement boundary conditions into RHS vector

      IF (IAND(IBDTYP,1).NE.0) THEN
        CALL PDSETX (DWORK(KRHS),DWORK(KRHS+NEQU),0D0,1D0,
     *               DASMBL(ORE),TRIAS(1,NLMAX),IASMBL,DASMBL,
     *               IGEOM,DGEOM)
      END IF

      IF (IAND(IBDTYP,2).NE.0) THEN
        CALL BDRSTX (DWORK(KRHS),DWORK(KRHS+NEQU),
     *               KWORK(L(TRIAS(OLXNPR,NLMAX))),
     *               PARX,PARY,UE,0D0,
     *               DASMBL(ORE),
     *               TRIAS(1,NLMAX),IASMBL,DASMBL,IGEOM,DGEOM) 
      END IF
      
C     -----------------------------------------------------------------
C     Step 3: Implement boundary conditions into solution vector

      KU = L(LUP)

      IF (IAND(IBDTYP,4).NE.0) THEN
        CALL BDRSTX (DWORK(KU),DWORK(KU+NEQU),
     *               KWORK(L(TRIAS(OLXNPR,NLMAX))),
     *               PARX,PARY,UE,0D0,
     *               DASMBL(ORE),
     *               TRIAS(1,NLMAX),IASMBL,DASMBL,IGEOM,DGEOM) 
      END IF

      END 
      
