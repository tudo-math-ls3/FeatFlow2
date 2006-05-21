***********************************************************************
* Initialize triangulations 
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
        
C       Put the TRIA structure to the COMMON block - so the FEAT library
C       assumes that there is no triangulation (because everything is 0)!
        CALL TRIA2C(TRIAS(1,NLMIN))
        
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
      
