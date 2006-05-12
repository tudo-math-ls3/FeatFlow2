***********************************************************************
* This file contains special output routines for writing only the grid
* to a GMV file. These routines can be used in the grid deformation
* process to visualise the grid deformation in time.
***********************************************************************

************************************************************************
* Colorize macros
*
* Determine the macro in the coarse mesh belonging 
* to an element on the finest level
* Uses logical function ISINEL in file search.f
*
* In:
*  NELM   - Number of elements in macro (=coarse) grid
*  KVERTM - Vertex information on coarse grid
*  NEL    - Number of elements in current grid
*  DCORVG - Coordinate information on current grid
*  KVERT  - Vertex information on current grid
*
* Out:
*  KCOLOR - array [1..NEL] of integer
*           Each entry KCOLOR(IEL) represents the number of the 
*           element on the coarse grid that contains the element IEL.
************************************************************************

      SUBROUTINE COLMCR(NELM,KVERTM,DCORVG,NEL,KVERT,KCOLOR)
      
      IMPLICIT NONE
      INCLUDE 'cbasictria.inc'
      
      DIMENSION KVERTM(NNVE,*),KVERT(NNVE,*),KCOLOR(*),DCORVG(2,*)
      INTEGER NELM,NEL,KVERTM,KVERT,KCOLOR
      DOUBLE PRECISION DCORVG
      INTEGER IEL,JEL
      DOUBLE PRECISION DXM,DYM
      LOGICAL ISINEL
      
      DO IEL=1,NEL
        DXM=0.5D0*(DCORVG(1,KVERT(1,IEL))+DCORVG(1,KVERT(3,IEL)))
        DYM=0.5D0*(DCORVG(2,KVERT(1,IEL))+DCORVG(2,KVERT(3,IEL)))
        DO JEL=1,NELM
          IF (ISINEL(DXM,DYM,JEL,DCORVG,KVERTM)) THEN
            KCOLOR(IEL)=JEL
            GOTO 310
          END IF
        END DO
310   END DO

      END 

************************************************************************** 
* Grid adaption GMV output
*
* This routine writes out GMV-files with data about the current
* grids while performing the grid adaption.
*
* In:
*  TRIA   - Triangulation structure describing the grid to write
*  IDX    - General index, user defined, to include in the filename
*  ILV    - Number of current level, to include in the filename;
*           =-1: don't include in the filename.
*  IMST   - Number of current macro step, to include in the filename
*           =-1: don't include in the filename.
*  ISMST  - Number of smooting step, to include in the filename
*           =-1: don't include in the filename.
*
* The following handles are input handles to various information
* arrays for the grid adaption. If <> 0, the corresponding information
* is written to the GMV file. If = 0, the corresponding information
* is not written to the GMV file.
*
*  NELM   - Number of elements in macro (=coarse)-grid
*  LCOLOR - Handle to: array [1..TRIA.NEL] of integer.
*           For each cell: Array with number of macro (cell in coarse grid)
*           that contains the cell.
*  LSIZE  - Handle to array [1..NVT] of double,
*           Average area around each vertex
*  LSIZEM - Handle to array [1..NEL] of double,
*           Area of each cell
*  LRATIO - Handle to array [1..NEL] of double,
*           Area ratio of one cell to neighbour cells
*  LMON   - Handle to array [1..NVT] of double,
*           Monitor function in the vertices 
*  LF     - Handle to array [1..NVT] of double,
*           RHS of linear system in the vertices
*  LU     - Handle to array [1..NVT] of double,
*           Solution of the linear system in the vertices
*  LGX    - Handle to array [1..NVT] of double,
*           X-gradient vector.
*  LGY    - Handle to array [1..NVT] of double,
*           Y-gradient vector.
*
*  The last parameter defines the name of the output files:
*
*  CFN    - The basic name of the files
*
* Out:
*  A GMV-file.
*
* Remarks: Unit no. 69 will be used for writing the file.
************************************************************************** 

      SUBROUTINE GAGMVS (TRIA,IDX,ILV,IMST,ISMST,
     *                   NELM,LCOLOR,LSIZE,LSIZEM,LRATIO,
     *                   LMON,LF,LU,LGX,LGY,
     *                   CFN)
     
      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cmem.inc'

      INCLUDE 'stria.inc'
      
      INCLUDE 'dstrings.inc'

      INTEGER IDX,ILV,IMST,ISMST,TRIA(SZTRIA)
      INTEGER NELM,LCOLOR,LMON,LF,LU,LGX,LGY,LSIZE,LRATIO,LSIZEM
      INTEGER NCELLS,NVERTS
      CHARACTER CFN*(*)
      
C local variables

      CHARACTER CFILE*(64)
      INTEGER LSTR, MFILE
      
      EXTERNAL GMVDG1,GMVDG2

C Create the filename

      LSTR = STNEWC (.TRUE.,CFN)
      CALL STCATC (LSTR,.TRUE.,'.gmv')
      IF (IDX.GT.-1) THEN
        CALL STCATC (LSTR,.TRUE.,'.')
        CALL STCATI (LSTR,IDX,1,.FALSE.)
      END IF
      IF (ILV.GT.-1) THEN
        CALL STCATC (LSTR,.TRUE.,'.')
        CALL STCATI (LSTR,ILV,2,.TRUE.)
      END IF
      IF (IMST.GT.-1) THEN
        CALL STCATC (LSTR,.TRUE.,'.')
        CALL STCATI (LSTR,IMST,3,.TRUE.)
      END IF
      IF (ISMST.GT.-1) THEN
        CALL STCATC (LSTR,.TRUE.,'.')
        CALL STCATI (LSTR,ISMST,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, CFILE)
      CALL STDIS (LSTR)

      IF (MT.GE.2) THEN
        WRITE (MTERM,'(A,A)') 'Writing GMV-file: ',CFILE
      END IF

C     Open the GMV-file for writing

      MFILE = 69
      CALL GMVOF0 (MFILE,-1,CFILE)
      CALL GMVHEA (MFILE)
      
C     Write current triangulation

      CALL GMVTRI (MFILE,TRIA,0,NCELLS,NVERTS)
      
C     Write material number for all cells.
C     The material colors the cell, so we either use the coloring
C     array to define the materials or (if this does not exist)
C     we use the dummy material number 1.

      IF (LCOLOR.NE.0) THEN
        CALL GMVMAT (MFILE,TRIA,0,NCELLS,TRIA(ONEL),
     *               -MAX(NELM,1),'',GMVDG2,
     *               KWORK(L(LCOLOR)),0)
      ELSE
        CALL GMVMAT (MFILE,TRIA,0,NCELLS,NELM,
     *               -1,'',GMVDG1,0,0)
      END IF
      
C     Write monitor function

      IF (LMON.NE.0) THEN
        CALL GMVSCA (MFILE,TRIA,1,NVERTS,TRIA(ONVT),
     *               DWORK(L(LMON)),'Monitor_function')
      END IF
      
C     Area of each cell

      IF (LSIZEM.NE.0) THEN
        CALL GMVSCA (MFILE,TRIA,0,NCELLS,TRIA(ONEL),
     *               DWORK(L(LSIZEM)),'Area')
      END IF

C     Average area around each vertex

      IF (LSIZE.NE.0) THEN
        CALL GMVSCA (MFILE,TRIA,1,NVERTS,TRIA(ONVT),
     *               DWORK(L(LSIZE)),'Mean_area')
      END IF

C     Area ratio of one cell to the neighbour cell

      IF (LRATIO.NE.0) THEN
        CALL GMVSCA (MFILE,TRIA,0,NCELLS,TRIA(ONEL),
     *               DWORK(L(LRATIO)),'Area_Ratio')
      END IF

C     Write solution/rhs 

      IF (LU.NE.0) THEN
        CALL GMVSCA (MFILE,TRIA,1,NVERTS,TRIA(ONVT),
     *               DWORK(L(LU)),'Solution')
      END IF

      IF (LF.NE.0) THEN
        CALL GMVSCA (MFILE,TRIA,1,NVERTS,TRIA(ONVT),
     *               DWORK(L(LF)),'RHS')
      END IF
      
C     Write gradients

      IF (LGX.NE.0) THEN
        CALL GMVSCA (MFILE,TRIA,1,NVERTS,TRIA(ONVT),
     *               DWORK(L(LGX)),'X_Gradient')
      END IF
      
      IF (LGY.NE.0) THEN
        CALL GMVSCA (MFILE,TRIA,1,NVERTS,TRIA(ONVT),
     *               DWORK(L(LGY)),'Y_Gradient')
      END IF
      
C     Close the GMV-file, finish

      CALL GMVFOT (MFILE)
      CLOSE (MFILE)
      
      END
      