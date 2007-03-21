*=======================================================================
      SUBROUTINE RDPARM (CDATA, MDATA)
*=======================================================================
*  Purpose:  - reads the parameter file 'cdata' and initializes
*              the corresponding data structure
*=======================================================================
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cparqdata.inc'
      
C parameters

      CHARACTER*(*) CDATA
      INTEGER MDATA

C local variables

      INTEGER IBCT,IIBCT,NCOMP,NNPAR,ICOMP,ITYP,NSPLIN
      INTEGER NPAR,KITYP,KNSPLN,KNPAR,KIPPTR,KXPAR,KYPAR,ICPTR
      INTEGER IPPTR,KKPAR,IPAR
      DOUBLE PRECISION XPAR,YPAR

C=======================================================================
C   open file 'cdata'
C=======================================================================
      OPEN(UNIT=MDATA,FILE=CDATA,ERR=222)
      GOTO 221
222   WRITE(*,*)'ERROR IN RDPARM: OPEN FAILED'
      STOP
221   CONTINUE
C=======================================================================
C   determine dimensions for creating arrays
C=======================================================================
      READ(MDATA,*,ERR=111)
      READ(MDATA,*,ERR=111) NBCT
C *** NBCT=0 MEANS: RDPARM IS A DUMMY PROCESS (NO PREPROCESSING)
      IF (NBCT.EQ.0) GOTO 99999
C *** LOOP OVER NBCT BOUNDARY PARTS
      NNNPAR=0
      NNCOMP=0
      DO IBCT=1,NBCT
        READ(MDATA,*,ERR=111)
        READ(MDATA,*,ERR=111) IIBCT
        IF(IBCT.NE.IIBCT) THEN
          WRITE(*,*)'ERROR IN RDPARM: CONFLICT WITH IBCT'
          STOP
        ENDIF
        READ(MDATA,*,ERR=111)
        READ(MDATA,*,ERR=111) NCOMP
        NNCOMP=NNCOMP + NCOMP
C
C ***   LOOP OVER NCOMP COMPONENTS OF PART IBCT
        NNPAR=0
        READ(MDATA,*,ERR=111)
        DO ICOMP=1,NCOMP
          READ(MDATA,*,ERR=111) ITYP, NSPLIN, NPAR
          NNPAR=NNPAR + NPAR
        ENDDO
        NNNPAR=NNNPAR + NNPAR
C
      ENDDO
C=======================================================================
C      read error
C=======================================================================
        GOTO 103
 111    WRITE(*,*)'ERROR IN RDPARM: READ ERROR --> STOP'
	STOP
 103    CONTINUE
C=======================================================================
C   allocation of all needed arrays
C=======================================================================
        CALL ZNEW(NBCT,3,LNCOMP,'NCOMP ')
        IF (IER.NE.0) GOTO 99998
        CALL ZNEW(NBCT,3,LICPTR,'ICPTR ')
        IF (IER.NE.0) GOTO 99998
        CALL ZNEW(NNCOMP,3,LITYP,'ITYP  ')
        IF (IER.NE.0) GOTO 99998
        CALL ZNEW(NNCOMP,3,LNSPLN,'NSPLN ')
        IF (IER.NE.0) GOTO 99998
        CALL ZNEW(NNCOMP,3,LNPAR,'NPAR  ')
        IF (IER.NE.0) GOTO 99998
        CALL ZNEW(NNCOMP,3,LIPPTR,'IPPTR ')
        IF (IER.NE.0) GOTO 99998
C
        CALL ZNEW(NNNPAR,1,LXPAR,'DXPAR ')
        IF (IER.NE.0) GOTO 99998
        CALL ZNEW(NNNPAR,1,LYPAR,'DYPAR ')
        IF (IER.NE.0) GOTO 99998
C=======================================================================
C   fill the arrays
C=======================================================================
        CLOSE(MDATA)
        OPEN(UNIT=MDATA,FILE=CDATA,ERR=222)
C--------------------------------------------------------------
C     determine KWORK/DWORK-pointers 
C--------------------------------------------------------------
      KITYP=L(LITYP)
      KNSPLN=L(LNSPLN)
      KNPAR=L(LNPAR)
      KIPPTR=L(LIPPTR)
      KXPAR=L(LXPAR)
      KYPAR=L(LYPAR)
C--------------------------------------------------------------
      READ(MDATA,*,ERR=111)
      READ(MDATA,*,ERR=111) 
C *** LOOP OVER NBCT BOUNDARY PARTS
      ICPTR=0
      IPPTR=0
      KKPAR=0
      DO IBCT=1,NBCT
C       SET POINTER ICPTR FOR THE BOUNDARY PART 'IBCT'
        KWORK(L(LICPTR)+IBCT-1)=ICPTR
C
        READ(MDATA,*,ERR=111)
        READ(MDATA,*,ERR=111) IIBCT
        IF(IBCT.NE.IIBCT) THEN
          WRITE(*,*)'ERROR IN RDPARM: CONFLICT WITH IBCT'
          STOP
        ENDIF
        READ(MDATA,*,ERR=111)
        READ(MDATA,*,ERR=111) NCOMP
        KWORK(L(LNCOMP)+IBCT-1)=NCOMP
C
C ***   LOOP OVER NCOMP COMPONENTS OF PART IBCT
        NNPAR=0
        READ(MDATA,*,ERR=111)
        DO ICOMP=0,NCOMP-1
          READ(MDATA,*,ERR=111) ITYP, NSPLIN, NPAR
          NNPAR=NNPAR + NPAR
          KWORK(KITYP+ICPTR)=ITYP
          KWORK(KNSPLN+ICPTR)=NSPLIN
          KWORK(KNPAR+ICPTR)=NPAR
          KWORK(KIPPTR+ICPTR)=IPPTR
          ICPTR=ICPTR+1
          IPPTR=IPPTR + NPAR
        ENDDO
C
      ENDDO
C
C ***   LOOP OVER NNNPAR PARAMETER PAIRS OF BOUNDARY PART IBCT
        READ(MDATA,*,ERR=111)
        KXPAR=L(LXPAR)
        KYPAR=L(LYPAR)
        DO IPAR=0,NNNPAR-1
          READ(MDATA,*,ERR=111) XPAR, YPAR
          DWORK(KXPAR+IPAR)=XPAR
          DWORK(KYPAR+IPAR)=YPAR
        ENDDO
C
C=======================================================================
C   error case for ZNEW
C=======================================================================
      GOTO 99999
99998 WRITE(*,*) 'ERROR IN RDPARM: ZNEW WITH IER=', IER
      STOP
C=======================================================================
99999 CONTINUE
C
      END
