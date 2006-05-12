***********************************************************************
* This file contains initialization routines for the problem of
* solving stationary and instationary Navier Stokes equations.
*
* To run the stationary/instationary Navier Stokes solver, the
* initialization/release routines are used as follows:
*
*  CALL INMAIN          ! Initialize main structures
*  CALL INMTRI          ! Initialize triangulations
*  CALL INSTSL/INISSL   ! Initialize solver
*  CALL INITRK          ! Open output for tracking of solutions
*  CALL GNMATD          ! Generate matrix entries for const. matrices
*
*  DO 
*
*    IF (stationary) THEN
*      CALL GNBDGE          ! Generate boundary/geometry information
*    END IF                 ! for stationary solver
*
*    CALL GNSTVC/GNISVC   ! Generate RHS/solution vector data
*
*    CALL the solver
*
*    do something with the solution, e.g. clear it
*
*  WHILE (something)
*
*  IF (stationary) THEN
*    CALL DNBDGE          ! Release boundary/geometry information
*  END IF
*
*  CALL DONTRK          ! Close files for tracking solution values
*  CALL DNSTSL/DNISSL   ! Release solver structures
*  CALL DNMTRI          ! Release triangulations
*  CALL DNMAIN          ! Release main structures
*
* The GNBDGE/DNBDGE routines must only be used for the call of the
* stationary solver. the nonstationary solver calls these routines
* automatically.
***********************************************************************

***********************************************************************
* Customization - Implementational details
* ----------------------------------------
* Here we give a general overview about some implementational details
* that differ from the standard implementation of structures/algorithms.
* Some of the structures are customized to the special situation
* of CC2D with moving boundaries and geometry support.
*
* 1.) Triangulation structures
* ----------------------------
* The STRIA-structures for all levels are generally kept as described in
* STRIA.INC. The following optional arrays are activated and used
* during the nonlinear iteration:
*   TRIA.DMBDP  - boundary midpoints
*   TRIA.DCORMG - midpoint coordinates
*   TRIA.DAREA  - Area of the elements
*
* The user-defined part of the TRIA structure is maintained as follows:
*   TRIA.TRIUD(1) = NBDMT
*                 = Number of boundary vertices
*   TRIA.TRIUD(2) = LSCNPR
*                 = Handle to array [1..NMT] of integer = KSCNPR
*     Shortcut nodal property for all edges. This contains a list of
*     all boundary nodes. The first NVBD entries contain the numbers
*     of all boundary edges. Then all fictitious boundary edges
*     follow, attached to that list. NBDMT is the actual size of
*     this array, as it may vary over time if the fictitious boundary
*     changes. As the maximum number of edges in the geometry is NMT,
*     NBDMT<=NMT. The entries follow the following rule:
*      |Entry|   : Number of the edge on the boundary.
*                  Value 1..NMT (not NVT+1..NVT+NMT !!!)
*      Entry > 0 : The edge is a Dirichlet edge
*      Entry < 0 : The edge is a Neumann edge
*
*     The shortcut nodal property must be recalculated if the geometry
*     changes. An array large enough for all the egdes is allocated in
*     advance together with the triangulation structure.
*
*   TRIA.TRIUD(3) = LIARR
*                   Handle to precalculated integer information about
*                   the current geometry or 0 if not used
*   TRIA.TRIUD(4) = LDARR
*                   Handle to precalculated double information about
*                   the current geometry or 0 if not used
*
***********************************************************************

***********************************************************************
* Initialize output channels
*
* This routine initializes the global output channels. Like ZINIT,
* it must be called prior to any other routines!
*
* The routine will read information about the output from the
* output DAT file, initialize file- and terminal output and returns
* in MSHOW the output level during the initialization phase.
*
* In:
*     -
* Out:
*   MSHOW  - Level of output to the terminal
*
* File handle 79 is used for reading the initialization file.
***********************************************************************

      SUBROUTINE OINIT (MSHOW)
      
      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'

      INTEGER MSHOW
      
C     local variables
      
      INTEGER MDATA
      CHARACTER CDATA*60

C     Read parameters and initialise the output structures.
C
C     Use unit number 79 for file output and open the file as 
C     prescribed in the DAT file.

      MDATA=79
      CDATA='data/output.dat'
      
      CALL RDOUT (MDATA,CDATA,62,MSHOW)

      END

***********************************************************************
* Main initialization routine
*
* This routine initializes the main global problem structures without
* acutally computing anything.
* Initializes structures, reads paramterization, reads DAT files.
*
* The routine has to be adapted problem-specifically. Remark that
* no matrix/vector initialization is done here, and no reading of
* grids/computation of subdivision is done! Only problem-specific
* data that has nothing to do with the actual discretization has to
* be initialized here!
*
* In:
*   MSHOW  - Level of output to the terminal
* Out:
*   IGEOM  - array [1..*] of integer 
*            = TIGeometryData
*   DGEOM  - array [1..*] of double 
*            = TDGeometryData
*            Integer- and double-precision parameter blocks with
*            geometry information.
***********************************************************************

      SUBROUTINE INMAIN (MSHOW,IGEOM,DGEOM)
      
      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cparametrization.inc'
      
      INTEGER MSHOW,IGEOM(*),DGEOM(*)
      
C     local variables
      
      INTEGER MDATA
      CHARACTER CDATA*60
      CHARACTER CGADAT*60,CGEODT*60
      
      CHARACTER CSTR*255

C     Now read parameters and initialise the global variables.
C
C     Use unit number 79 for file output and open the file as 
C     prescribed in the DAT file.

      MDATA=79
      
      CDATA = 'data/paramtriang.dat'
      CALL RDPARD (MDATA,MSHOW,CDATA)
      
      CDATA = 'data/discretization.dat'
      CALL RDDISC (MDATA,MSHOW,CDATA)
      
      CDATA = 'data/linsol_cc2d.dat'
      CALL RDLSCC (MDATA,MSHOW,CDATA)

      CDATA = 'data/nonlinsol_cc2d.dat'
      CALL RDNLCC (MDATA,MSHOW,CDATA)
      
      CDATA = 'data/timediscr.dat'
      CALL RDTIDC (MDATA,MSHOW,CDATA)

      CDATA = 'data/postprocessing.dat'
      CALL RDPOST (MDATA,MSHOW,CDATA)
      
C     The variables for output are read in. Now we can start reading
C     all the various variables...

C     Read and set up geometry data

      CGEODT = 'data/geometry.dat'
      CALL RDGEOM (MDATA,MSHOW,CGEODT)
      CALL INIGEO (IGEOM,DGEOM,1)
      
C     Initialize fictitious boundary structures

      CALL FBDCRE (IGEOM,DGEOM)
      
C     Read parameters of the grid adaption DAT file
C     Use a string variable as some compilers have problems using a
C     direct string as parameter!

      CGADAT = 'data/gridadapt.dat'
      CALL RDGSMT (MDATA,MSHOW,CGADAT)
      
C     Read parameters for the parametrisation

      CALL GENPAR (.TRUE.,IMESH,CPARM)
           
      WRITE (CSTR,9000) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
9000  FORMAT(79('-'))

      END
      
***********************************************************************
* Cleanup routine for the structures allocated in INMAIN
*
* This routine cleans up the memory that was allocated in the function
* INMAIN. It cleans up the paramterization, geometry information
* and some global structures.
*
* In:
*   IGEOM  - array [1..*] of integer 
*            = TIGeometryData
*   DGEOM  - array [1..*] of double 
*            = TDGeometryData
*            Integer- and double-precision parameter blocks with
*            geometry information.
***********************************************************************

      SUBROUTINE DNMAIN (IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      
C     Clean up the parametrization

      CALL DISPAR ()
      
C     Clean up fictitious boundary structures

      CALL FBDDIS (IGEOM,DGEOM)
      
C     Delete allocated memory of geometries

      CALL DONGEO (IGEOM)
      
C     Nothing to be cleaned up here for the grid adaption.
C     Nothing to be cleaned up here for the optimizer.
C     Nothing to be cleaned up here for the geometry.
C     Nothing to be cleaned up here for the parameters of the DAT file.
C     Nothing to be cleaned up here for the output channels.

      END

***********************************************************************
* Initialize triangulations for stationary and/or nonstationary
* solver.
*
* This routine initializes a TRIAS structure with information about
* all levels where the computation should be performed. The coarse
* grid is initialized by reading in a file, the finer grids
* by reading or by regular refinement. The minimum/maximum level
* is corrected if necessary.
*
* On each level, an array of length NMT(level) is allocated to take
* the shortcut nodal property on each edge.
*
* In:
*   MSHOW   - Level of output
*   MFILE   - Handle of file where to write additional output to
*
*   NLMIN   - Minimum level where the later computation will be
*             performed; if <= 0, the minimum level is calculated
*             as NLMAX-NLMIN
*   NLMAX   - Maximum level where the later computation will be
*             performed
*   IMETH   - Method how to generate the grids.
*             =0: Read the coarse mesh from a file and refine; standard
*             =1: Read all meshes from a file
*   CFILE   - Name of the file
*
* Out:
*   NLMIN   - If NLMAX>9 or NLMIN<0, NLMIN/NLMAX is corrected to be in
*             the range 1..9; the meshes are set up appropriately.
*   NLMAX   - If NLMAX>9 or NLMIN<0, NLMIN/NLMAX is corrected to be in
*             the range 1..9; the meshes are set up appropriately.
*   TRIAS   - array [1..SZRIA,1..NNLEV] of integer
*             Triangulation structure; Level NLMIN..NLMAX is filled
*             with data.
***********************************************************************

      SUBROUTINE INMTRI (MSHOW,TRIAS,NLMIN,NLMAX,IMETH,CFILE)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'
      
C     parameters
      
      INTEGER MSHOW,TRIAS(SZTRIA,*),NLMIN,NLMAX,IMETH
      CHARACTER CFILE*(*)
      
C     local variables
      
      INTEGER I,IREF,MF
      CHARACTER CSTR*(255)
      
C     A negative NLMIN corresponds to a "relative" minimum level; correct
C     this:

      IF (NLMIN.LE.0) NLMIN = NLMAX+NLMIN
      
C     Correct if the maximum level is too large.
C     IREF specifies the number of pre-refinements of
C     the coarse grid if the mesh is not read in from a file. If NLMAX is
C     set > 9, we shift NLMIN and NLMAX by some levels so that NLMAX is 9.

      IREF = 0
      IF (NLMAX.GT.9) THEN
        IREF = NLMAX-9
        NLMIN = MAX(1,NLMIN-IREF)
        NLMAX = 9
        WRITE (MTERM,'(A)') 'Warning: NLMAX too large.'
        WRITE (MTERM,'(A)') 'Performing level shift for calculation '//
     *                      'on finer levels.'
        WRITE (MTERM,'(A,I3,A,I3)') 'Calculation will be on level ',
     *                      NLMIN+IREF,'..',NLMAX+IREF
      END IF
      
C     For IMETH=0, read the coarse mesh and refine it:

      IF (IMETH.EQ.0) THEN
      
C       Clear the TRIA-structure

        CALL LCL3(TRIAS(1,NLMIN),SZTRIA)
        
C       Set up the coarse grid, read it from the file, pre-refine it

        MF = 0
        CALL GENTRI (TRIAS(1,NLMIN), 2, 0, IREF+NLMIN-1, 1,
     *               0, TRIAS(1,NLMIN), 
     *               MF, 1, CFILE)
     
C       Create the other grids by refinement.
C       Share the DCORVG-information between all levels.

        DO I=NLMIN+1,NLMAX
          CALL LCL3(TRIAS(1,I),SZTRIA)
          CALL GENTRI (TRIAS(1,I), 1, 0, 1, 1,
     *                 1, TRIAS(1,I-1), 
     *                 MFILE, 1, CFILE)
       
        END DO
        
C       Only the LCORVG handle of the maximum level is valid at the
C       moment. Restore LCORVG of the lower levels to a valid handle
C       so that DCORVG is appropriately shared between all levels

        DO I=NLMAX-1,NLMIN,-1
          CALL TRICSH (TRIAS(1,I+1),TRIAS(1,I))
        END DO
      
      ELSE
      
C       For IMETH=1, read in all grids from one single file.
C      
C       Clear the TRIA-structure

        CALL LCL3(TRIAS(1,NLMIN),SZTRIA)
        
C       Read all grids. The first call will open the file.

        MF = 0
        
        DO I=NLMIN,NLMAX
          CALL LCL3(TRIAS(1,I),SZTRIA)
          CALL GENTRI (TRIAS(1,I), 3, 0, 1, 1,
     *                 0, TRIAS(1,I), 
     *                 MF, 1, CFILE)
        END DO
        
C       Close the file afterwards

        CLOSE (MF)
      
      END IF
      
C     print mesh statistics

      DO I = NLMIN,NLMAX
        WRITE(CSTR,'(A,5I10)') 'ILEV,NVT,NMT,NEL,NVBD: ', I+IREF,
     *        TRIAS(ONVT,I),TRIAS(ONMT,I),TRIAS(ONEL,I),TRIAS(ONVBD,I)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      END DO

      END
      
***********************************************************************
* Clean up triangulations for stationary and/or nonstationary
* solver.
*
* This routine cleans up the triangulation structures,
* releasing all allocated memory.
*
* In:
*   NLMIN   - Minimum level where the later computation will be
*             performed; if <= 0, the minimum level is calculated
*             as NLMAX-NLMIN
*   NLMAX   - Maximum level where the later computation will be
*             performed
*   TRIAS   - array [1..SZRIA,1..NNLEV] of integer
*             Triangulation structure; Level NLMIN..NLMAX is filled
*             with data.
* Out:
*   All handles in TRIAS will be invalud
***********************************************************************

      SUBROUTINE DNMTRI (NLMIN,NLMAX,TRIAS)

      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
C     parameters
      
      INTEGER TRIAS(SZTRIA,*),NLMIN,NLMAX
      
C     local variables

      INTEGER I
      
C     Loop through all essential levels from back to front...
C
C     Delete the triangulations.

      DO I=NLMAX,NLMIN,-1
        CALL TRIDEL(TRIAS(1,I))
      END DO

      END 
      
      