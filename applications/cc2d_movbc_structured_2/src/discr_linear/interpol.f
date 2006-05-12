************************************************************************
* Constant <-> nonconforming interpolation
* 
* This converts between constant and nonconforming piecewise linear
* FE pressure functions.
* IPAR=0 : convert piecewise constant DPC to piecewise linear DPL
* IPAR=1 : convert piecewise linear DPL to piecewise constant DPC
*
* DPC    : array [1..NEL] of double
* DPL    : array [1..NMT] of double
************************************************************************

      SUBROUTINE C2N2DM (DPC,DPL,KMID,KADJ,NEL,NMT,NVT,IPAR)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      DOUBLE PRECISION DPC(*),DPL(*)
      INTEGER KMID(NNVE,*),KADJ(NNVE,*)
      INTEGER NEL,NMT,NVT,IPAR

C local variables

      DOUBLE PRECISION DPH
      INTEGER IEL, IVE, IADJ, IMID

C-----------------------------------------------------------------------
C
C *** constant to nonconforming = 0
      IF (IPAR.EQ.0) THEN
C
       DO 10 IEL=1,NEL
       DPH  =DPC(IEL)
C
       DO 20 IVE=1,4
       IADJ=KADJ(IVE,IEL)
       IMID=KMID(IVE,IEL)-NVT
C
       IF (IADJ.EQ.0)   DPL(IMID)=DPH
       IF (IADJ.GT.IEL) DPL(IMID)=0.5D0*(DPH+DPC(IADJ))
C
20     CONTINUE
10     CONTINUE
C
      ELSE
C
       DO 110 IEL=1,NEL
       DPC(IEL)=0.25D0*( DPL(KMID(1,IEL)-NVT)+DPL(KMID(2,IEL)-NVT)
     *                  +DPL(KMID(3,IEL)-NVT)+DPL(KMID(4,IEL)-NVT))
110    CONTINUE
C
      ENDIF      

      END
