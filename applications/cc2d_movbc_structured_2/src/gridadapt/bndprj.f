************************************************************************
* Projects boundary nodes back to the boundary (on the old bndry. component)
*
* In: 
*  DXN,DYN (dp): x/y coordinate of a point
*  DPAR (dp): (old) parameter value of the point
*  IBCT (int): index of the boundary container
*
* Out:
*  DXN,DYN: new coordinates of the point after projection to the boundary
*  DPAR: new parameter value after projection
************************************************************************

      SUBROUTINE BNDPRJ(DXN,DYN,DPAR,IBCT)

      IMPLICIT NONE

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INCLUDE 'cmem.inc'
      INCLUDE 'cparqdata.inc'
      
      DOUBLE PRECISION DXN,DYN,DPAR
      INTEGER IBCT

      INTEGER ICOMP,ICPTR,IICOMP,ITYP,IPPTR
      INTEGER KXPAR,KYPAR
      DOUBLE PRECISION DNORM1,DNORM2,DPROD
      DOUBLE PRECISION DXM,DYM,DRAD,DPHI1,DPHI2,DX1,DX2,DY1,DY2,DPHI
      
      DOUBLE PRECISION PARX,PARY
      
C *** determine old parameter value and the segment (boundary component)
      ICOMP=DPAR+1D0
C *** determine ityp
      ICPTR=KWORK(L(LICPTR)+IBCT-1)
      IICOMP=ICPTR+ICOMP-1
      ITYP=KWORK(L(LITYP)+IICOMP)
      IPPTR=KWORK(L(LIPPTR)+IICOMP)
      KXPAR=L(LXPAR)+IPPTR
      KYPAR=L(LYPAR)+IPPTR
C
C
C *** ITYP=1 LINE
      IF (ITYP.EQ.1) THEN
        DNORM1=DSQRT(DWORK(KXPAR+1)*DWORK(KXPAR+1)+
     *               DWORK(KYPAR+1)*DWORK(KYPAR+1))
        DXN=DXN-DWORK(KXPAR)
        DYN=DYN-DWORK(KYPAR)
        DPROD=DXN*DWORK(KXPAR+1)+DYN*DWORK(KYPAR+1)
        DXN=DPROD*DWORK(KXPAR+1)/(DNORM1*DNORM1)+DWORK(KXPAR)
        DYN=DPROD*DWORK(KYPAR+1)/(DNORM1*DNORM1)+DWORK(KYPAR)
        DNORM2=DSQRT((DXN-DWORK(KXPAR))*
     *               (DXN-DWORK(KXPAR))+
     *               (DYN-DWORK(KYPAR))*
     *               (DYN-DWORK(KYPAR)))
        DPAR=DNORM2/DNORM1
        IF (DPAR.GT.1D0) DPAR=0.9999d0
        DPAR=DPAR+ICOMP-1
C
C
C *** ITYP=2 PART OF A CIRCLE
      ELSE IF (ITYP.EQ.2) THEN
        DXM=DWORK(KXPAR)
        DYM=DWORK(KYPAR)
        DRAD=DWORK(KXPAR+1)
        DPHI1=DWORK(KXPAR+2)
        DPHI2=DWORK(KYPAR+2)
C
        DX1=DXN
        DX2=DYN
        DY1=DXM+DRAD*COS(DPHI1)
        DY2=DYM+DRAD*SIN(DPHI1)
C      
        DPHI=ACOS((DX1*DY1+DX2*DY2)/
     *             (DSQRT(DX1*DX1+DX2*DX2)*DSQRT(DY1*DY1+DY2*DY2)))
        IF (DY1*DX2-DY2*DX1.LT.0D0) DPHI=2*PI-DPHI
C
c        DPAR=DPHI/(DPHI2-DPHI1)
        DPAR=DPHI/(DPHI2-DPHI1)+ICOMP-1D0
C
C
      ELSE 
        WRITE (*,FMT='(A)') 'ityp=3 not supportet yet'
        STOP
      END IF !ITYP

      DXN=PARX(DPAR,IBCT)
      DYN=PARY(DPAR,IBCT)
C
      END SUBROUTINE
