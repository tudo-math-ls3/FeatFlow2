**********************************************************************
* (Depending on the implementation:)   User defined monitor function
*
* This file defines the monitor functions for the grid deformation
* process. These functions are dependent of the geometry of the
* moving boundary objects and have to be defined properly in order
* to deform the grid correctly.
*
* The standard implementation of the monitor function uses the
* fictitious boundary routines to determine the interface of the
* objects inside of the domain. If a special monitor function has
* to be defined, the user has to implement it either here or 
* using the fictitious boundary routines that are called here.
*
* Legacy code that does not use the fictitious boundary "library"
* can directly define the monitor function here.
**********************************************************************

**********************************************************************
* Purpose: Defines the monitor function and so the cell distribution
*
* Returns for the given point (X,Y) a value in the interval
* [EPS0,1] for the monitor function, where EPS0>0 is a given parameter
* by the user in the .DAT-file.
*
* In:
*  X,Y    - coordinates of the point
*  INODE  - Number of the node corresponding to (X,Y). Can be 0.
*           If INODE=0, the monitor function should be evaluated
*                       for a general point, only given by (X,Y)
*           If INODE>0, the parameter TRIA must also be defined.
*                       In this case INODE is a node number in the
*                       triangulation TRIA with coordinates (X,Y)
*  EPS0   - minimum monitor function values; values below that
*           will be truncated.
*  TRIA   - array [1..SZTRIA] of integer
*           If INODE>0, TRIA defines the triangulation INODE refers to.
*           If INODE=0, this parameter can be a dummy parameter.
*  LIARR  - integer; might be 0 
*           Handle to integer array with precalculated information.
*  LDARR  - integer; might be 0 
*           Handle to double array with precalculated information.
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*  IDATA  - array [1..*] of integer 
*  DDATA  - array [1..*] of double 
*           User defined integer- and double-precision parameter 
*           blocks, which was passed to XGASTA.
*
* LIARR and LDARR are handles to integer/double arrays with
* precalculated information. If INODE=0, these variables are not used.
* If INODE>0 these variables can be 0 to indicate that there is no
* precalculated information.
*
* Out:
*  Return value = value of the monitor function.
*
**********************************************************************
      
      DOUBLE PRECISION FUNCTION FMON(X,Y,INODE,EPS0GS,TRIA,LIARR,LDARR,
     *                               IGEOM,DGEOM,IDATA,DDATA)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'

C parameters

      DOUBLE PRECISION X,Y,EPS0GS
      INTEGER INODE,LIARR,LDARR,TRIA(*)
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
C externals
      
      DOUBLE PRECISION FBDMON
      INTEGER NFBDYC
      EXTERNAL FBDMON,NFBDYC
      DOUBLE PRECISION FBMMON
      EXTERNAL FBMMON
      EXTERNAL E011

C local variables 

      DOUBLE PRECISION F,TMP1,TMP2,TMP3,G
      INTEGER I,MST

C      INCLUDE 'ciniopt.inc'

C      DOUBLE PRECISION X,Y
      
C      DOUBLE PRECISION DIST

C      DIST=DSQRT((X-DCXPOS)**2+(Y-DCYPOS)**2)
!      FMON=MIN(1D0,MAX(ABS(DIST-0.05D0)*4D0,EPS0GS))
!      FMON=MIN(1D0,MAX(ABS(DIST-DRAD2)*4D0,EPS0GS))
!      FMON=MIN(1D0,MAX(ABS(DIST-DRAD2)*8D0,EPS0GS))

C      IF (DIST.LT.DCRAD) THEN
C        FMON=MIN(1D0,MAX(ABS(DIST-DCRAD)*5D0,EPS0GS))
C      ELSE
C        FMON=MIN(1D0,MAX(ABS(DIST-DCRAD)*1D0,EPS0GS))
C      END IF

C      IF (DIST.LT.DCRAD) THEN
C        FMON=MIN(1D0,MAX(ABS(DIST-DCRAD)*3D0,EPS0GS))
C      ELSE
C        FMON=MIN(1D0,MAX(ABS(DIST-DCRAD)*1D0,EPS0GS))
C      END IF

!      FMON=X+0.01D0

C Get the EPS0-parameter grom the GRIDSMOOTH-Common-blocks.
C Go through the fictitious boundary components and search for the
C lowest value. EPS0>0 is the minimum bound.
C All returned values are restricted to the interval [EPS0,1].
C
C When using an error distribution as monitor function,
C we use IDATA and DDATA to pass necessary information
C to the monitor-function callback routine. IDATA has the
C following internal structure:
C
C IDATA(1) - =0: only use geometry details for monitor function
C            =1: only use error distribution as monitor function
C            =2: use both, geometric details and error distribution
C IDATA(2) - Minimum level in SOLTRI
C IDATA(3) - Maximum level in SOLTRI, corresponding to solution vector
C IDATA(4..*) - This saves a copy of the mesh information in SOLTRI,
C               which defines the underlying mesh for the error
C               distribution
C DDATA(1..NVT(NLMAX)) - Saves the error distribuition in all
C                        vertices of the triangulation SOLTRI.
C                        All values are assumed to be in the 
C                        range [0..1].
C
C     Should we use geometric features for the monitor function?

      F = 1D0
      
      IF ((IDATA(1).EQ.0).OR.(IDATA(1).EQ.2)) THEN

C       Handle general fictitious boundary case:

        IF ((INODE.EQ.0).OR.(LIARR.EQ.0).OR.(LDARR.EQ.0)) THEN
        
C         No precalculation objects of mesh-specific information is used:

          F = MAX(EPS0GS, MIN(F,FBDMON(X,Y,0,EPS0GS,IGEOM,DGEOM)))
          
        ELSE
        
C         Precalculation objects are used:

          F = MAX(EPS0GS,MIN(F,FBMMON(INODE,TRIA,LIARR,LDARR,0,EPS0GS,
     *                                  IGEOM,DGEOM)))
     
        END IF
        
C       Loop through the fictitious boundary components to get 
C       the lowest value

        DO I=TRIA(ONBCT)+1,TRIA(ONBCT)+NFBDYC(IGEOM,DGEOM)
        
          IF ((INODE.EQ.0).OR.(LIARR.EQ.0).OR.(LDARR.EQ.0)) THEN
          
C           No precalculation objects of mesh-specific information is used:

            F = MAX(EPS0GS, MIN(F,FBDMON(X,Y,I,EPS0GS,IGEOM,DGEOM)))
            
          ELSE
          
C           Precalculation objects are used:

            F = MAX(EPS0GS,MIN(F,FBMMON(INODE,TRIA,LIARR,LDARR,I,EPS0GS,
     *                                  IGEOM,DGEOM)))
     
          END IF

        END DO
        
      END IF

C      F = MIN(2D0*F,1D0)

C     Should we use an error distribution for the monitor function?

      IF ((IDATA(1).EQ.1).OR.(IDATA(1).EQ.2)) THEN
      
C       DDATA is the error distribution, represented Q1 elements.
C       IDATA(2..*) describes the corresponding triangulation.
C       We need the value in the vertex X,Y. For that purpose,
C       we have to know the element that contains X,Y.
C       Start the hierarchical search routine to search for that
C       point!

        CALL PSRCH5(IDATA(4),IDATA(2),IDATA(3),X,Y,0,I)
        
        IF (I.NE.0) THEN
        
C         Calculate the starting address of the triangulation
C         structure corresponding to "NLMAX" inside of IDATA:

          MST = 4-1+(IDATA(3)-1)*SZTRIA

C         Evaluate the FE-function given in DDATA in X,Y.

          CALL SCEVLQ (IDATA(MST+ONVT),DDATA,E011,.FALSE.,X,Y,I,
     *                 IDATA(MST+1),DWORK(L(IDATA(MST+OLCORVG))),
     *                 KWORK(L(IDATA(MST+OLVERT))),
     *                 KWORK(L(IDATA(MST+OLMID))),
     *                 TMP1, TMP2, TMP3)
     
C         TMP1 is the result - and will be in the range [0..1],
C         as this is the only allowed range for the whole error
C         distribution. Take this as additional monitor function
C         value!

          TMP1 = MAX(EPS0GS,TMP1)
C          TMP1 = (-TMP1+2D0*TMP1*EPS0GS-EPS0GS)/(-1D0+EPS0GS)

          F = MIN(F,TMP1)
         
        END IF

      END IF
      
C     Final check:

      IF (F.LE.0D0) THEN
      
C       Should not happen - if the user sets the EPS-parameter correctly!
C       We don't test for EPS0=0 here because we allow EPS0=0 in the 
C       case if the fictitious boundary implementation is hopefully 
C       working correctly...

        WRITE (MTERM,*) 'Error: Monitor function <= 0 for point ',X,Y
        WRITE (MTERM,*) 'Set EPS0-parameter > 0 to avoid this!'
        STOP
        
      END IF

      FMON = F

      END 
