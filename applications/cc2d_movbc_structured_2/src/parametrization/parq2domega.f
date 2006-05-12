************************************************************************
* This file contains the OMEGA parametrisation routines for the
* computational domain:
*   OPARX - Obtain the X-coordinate of a point from its parameter value
*   OPARX - Obtain the Y-coordinate of a point from its parameter value
*   OTMAX - Obtain the maximum parameter value on a boundary component
*
* The implementation here is used as the standard implementation in case
* of that the domain is read in by RDPARM from an OMEGA-style file.
************************************************************************
*   parxypre.f:
*     - contains the parametrization functions:
*           OPARX, OPARY, OTMAX, OTNBC
*       for the preprocessing variant ('rdparm')
*     - the COMMON-block /TDATA/ is used frequently and has to be
*       initialized by the modul 'initxx.f'
************************************************************************

************************************************************************
* Get X-coordinate of parameter value.
*
* This routine returns the X-coordinate of a parameter value.
*
* In:
*   T      - parameter value
*   IBCT   - number of boundary component
*
* Return:
*   X-coordinate of parameter value T on boundary component IBCT.
************************************************************************

      DOUBLE PRECISION FUNCTION OPARX(T,IBCT)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cparqdata.inc'

C parameters

      DOUBLE PRECISION T
      INTEGER IBCT

C local variables
      
      INTEGER NCOMP,ICPTR,ICOMP,IICOMP,ITYP,NSPLIN,NPAR,IPPTR
      INTEGER KXPAR,KYPAR
      DOUBLE PRECISION TDIFF,X1,XDIFF,XM,R,PHI,PHI1,PHI2,PHI3,PHI4
      
      INTEGER ISPLIN,KK
      DOUBLE PRECISION TTDIF,TT,C1,C2,C3,C4
      
C=======================================================================
C *** IBCT dependent quantities
      NCOMP=KWORK(L(LNCOMP)+IBCT-1)
      ICPTR=KWORK(L(LICPTR)+IBCT-1)
C
C *** determine boundary component ICOMP of part IBCT and
C     the ICOMP-dependent quantities
C
      ICOMP=T+1.d0
      tdiff=T-ICOMP+1.d0
      if(tdiff.lt.0d0 .or. tdiff.ge.1d0) then
        write(*,*)'error in OPARX (userpre): conflict with ICOMP'
        stop
      endif
      iicomp=ICPTR+ICOMP-1
      ITYP=KWORK(L(LITYP)+iicomp)
      NSPLIN=KWORK(L(LNSPLN)+iicomp)
      NPAR=KWORK(L(LNPAR)+iicomp)
      IPPTR=KWORK(L(LIPPTR)+iicomp)
C
C *** pointer for the LXPAR/LYPAR-arrays
      KXPAR=L(LXPAR)+IPPTR
      KYPAR=L(LYPAR)+IPPTR      
C=======================================================================
C   ITYP=1:  line      
C=======================================================================
      if(ITYP.eq.1) then
C      
        if(NSPLIN.ne.1 .or. NPAR.ne.2) then
          write(*,*)'error in OPARX (userpre): conflict at ITYP=1'
          stop
        endif
        x1=DWORK(KXPAR)
        xdiff=DWORK(KXPAR+1)
        OPARX=x1 + xdiff*tdiff        
C=======================================================================
C   ITYP=2:  part of a circle      
C=======================================================================        
      else if(ITYP.eq.2) then
C      
        if(NSPLIN.ne.1 .or. NPAR.ne.3) then
          write(*,*)'error in OPARX (userpre): conflict at ITYP=2'
          stop
        endif
        xm=DWORK(KXPAR)
        r=DWORK(KXPAR+1)
        phi1=DWORK(KXPAR+2)
        phi2=DWORK(KYPAR+2)
        phi=phi1 + tdiff*(phi2-phi1)
        OPARX=xm + r*COS(phi)
C=======================================================================
C   ITYP=3:  spline
C=======================================================================        
      else if(ITYP.eq.3) then
C      
C *** determine "number of the subspline -1"=ISPLIN and tt=t_tilde
        ttdif=tdiff*NSPLIN
        ISPLIN=ttdif
        tt=ttdif-ISPLIN
        if(tt.lt.0d0 .or. tt.ge.1d0 .or. ISPLIN.gt.NSPLIN-1) then
          write(*,*)'error in OPARX (userpre): conflict with ISPLIN'
          stop
        endif
C      
C *** determine cubic spline-function values at 'tt'
        phi1=(2d0*tt*tt - 3d0*tt)*tt + 1d0
        phi2=(-2d0*tt + 3d0)*tt*tt
        phi3=(tt*tt - 2d0*tt + 1d0)*tt
        phi4=(tt - 1d0)*tt*tt
C      
C *** determine parameter constants
        kk=KXPAR + 4*ISPLIN
        c1=DWORK(kk)
        c2=DWORK(kk+1)
        c3=DWORK(kk+2)
        c4=DWORK(kk+3)
C
        OPARX=c1*phi1+c2*phi2+c3*phi3+c4*phi4
C=======================================================================      
      else
        write(*,*)'error in OPARX (userpre): wrong ITYP-value'
        stop
      endif
C
      END


************************************************************************
* Get Y-coordinate of parameter value.
*
* This routine returns the Y-coordinate of a parameter value.
*
* In:
*   T      - parameter value
*   IBCT   - number of boundary component
*
* Return:
*   Y-coordinate of parameter value T on boundary component IBCT.
************************************************************************

      DOUBLE PRECISION FUNCTION OPARY(T,IBCT)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cparqdata.inc'

C parameters

      DOUBLE PRECISION T
      INTEGER IBCT

C local variables
      
      INTEGER NCOMP,ICPTR,ICOMP,IICOMP,ITYP,NSPLIN,NPAR,IPPTR
      INTEGER KXPAR,KYPAR
      DOUBLE PRECISION TDIFF,R,PHI,PHI1,PHI2,PHI3,PHI4
      
      INTEGER ISPLIN,KK
      DOUBLE PRECISION TTDIF,TT,C1,C2,C3,C4,Y1,YDIFF,YM
C=======================================================================
C *** IBCT dependent quantities
      NCOMP=KWORK(L(LNCOMP)+IBCT-1)
      ICPTR=KWORK(L(LICPTR)+IBCT-1)
C
C *** determine boundary component ICOMP of part IBCT and
C     the ICOMP-dependent quantities
C
      ICOMP=T+1.d0
      tdiff=T-ICOMP+1.d0
      if(tdiff.lt.0d0 .or. tdiff.ge.1d0) then
        write(*,*)'error in OPARY (userpre): conflict with ICOMP'
        stop
      endif
      iicomp=ICPTR+ICOMP-1
      ITYP=KWORK(L(LITYP)+iicomp)
      NSPLIN=KWORK(L(LNSPLN)+iicomp)
      NPAR=KWORK(L(LNPAR)+iicomp)
      IPPTR=KWORK(L(LIPPTR)+iicomp)
C
C *** pointer for the LXPAR/LYPAR-arrays
      KXPAR=L(LXPAR)+IPPTR
      KYPAR=L(LYPAR)+IPPTR      
C=======================================================================
C   ITYP=1:  line      
C=======================================================================
      if(ITYP.eq.1) then
C      
        if(NSPLIN.ne.1 .or. NPAR.ne.2) then
          write(*,*)'error in OPARY (userpre): conflict at ITYP=1'
          stop
        endif
        y1=DWORK(KYPAR)
        ydiff=DWORK(KYPAR+1)
        OPARY=y1 + ydiff*tdiff        
C=======================================================================
C   ITYP=2:  part of a circle      
C=======================================================================        
      else if(ITYP.eq.2) then
C      
        if(NSPLIN.ne.1 .or. NPAR.ne.3) then
          write(*,*)'error in OPARY (userpre): conflict at ITYP=2'
          stop
        endif
        ym=DWORK(KYPAR)
        r=DWORK(KXPAR+1)
        phi1=DWORK(KXPAR+2)
        phi2=DWORK(KYPAR+2)
        phi=phi1 + tdiff*(phi2-phi1)
        OPARY=ym + r*SIN(phi)
C=======================================================================
C   ITYP=3:  spline
C=======================================================================        
      else if(ITYP.eq.3) then
C      
C *** determine "number of the subspline -1"=ISPLIN and tt=t_tilde
        ttdif=tdiff*NSPLIN
        ISPLIN=ttdif
        tt=ttdif-ISPLIN
        if(tt.lt.0d0 .or. tt.ge.1d0 .or. ISPLIN.gt.NSPLIN-1) then
          write(*,*)'error in OPARY (userpre): conflict with ISPLIN'
          stop
        endif
C      
C *** determine cubic spline-function values at 'tt'
        phi1=(2d0*tt*tt - 3d0*tt)*tt + 1d0
        phi2=(-2d0*tt + 3d0)*tt*tt
        phi3=(tt*tt - 2d0*tt + 1d0)*tt
        phi4=(tt - 1d0)*tt*tt
C      
C *** determine parameter constants
        kk=KYPAR + 4*ISPLIN
        c1=DWORK(kk)
        c2=DWORK(kk+1)
        c3=DWORK(kk+2)
        c4=DWORK(kk+3)
C
        OPARY=c1*phi1+c2*phi2+c3*phi3+c4*phi4
C=======================================================================      
      else
        write(*,*)'error in OPARY (userpre): wrong ITYP-value'
        stop
      endif
C
      END


************************************************************************
* Get maximum parameter value
*
* This routine returns the maximum parameter value on a boundary
* component
*
* In:
*   IBCT   - number of boundary component
*
* Return:
*   Maximum parameter value on that boundary component.
************************************************************************

      DOUBLE PRECISION FUNCTION OTMAX(IBCT)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cparqdata.inc'

C parameters

      INTEGER IBCT
      
C local variables 

      INTEGER NCOMP

C=======================================================================

      NCOMP=KWORK(L(LNCOMP)+IBCT-1)

      OTMAX=NCOMP
      
      END

************************************************************************
* Get number of real boundary components
*
* This routine returns the number of real boundary components in the
* parametrization.
*
* In:
*   -
*
* Return:
*   Number of boundary components
************************************************************************

      INTEGER FUNCTION OTNBC ()

      IMPLICIT NONE
      
      INCLUDE 'cparqdata.inc'

      OTNBC = NBCT

      END
