************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* ORAn                                                                 *
*                                                                      *
* Purpose  Get arrays from file, previously stored by XOWA or OWAn     *
*                                                                      *
* Subroutines/functions called  ORA0                                   *
*                                                                      *
* Version from  03/21/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----    ----                                                        *
* MFILE   I*4    I/O unit - previouly opened (by OF0)                  *
* IFMT    I*4    0  Format free input                                  *
*                1  Formatted input                                    *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* DX      R*8    Double precision array - version ORA1                 *
* VX      R*4    Single precision array - version ORA2                 *
* KX      I*4    Integer array          - version ORA3                 *
* ARR     C*6    Name of array read from data file                     *
* IER     I*4    Error indicator                                       *
*                -110  Error while reading from unit  MFILE            *
*                -113  Wrong value of  ITYPE  read from  MFILE         *
*                                                                      *
************************************************************************
C
      SUBROUTINE ORA1(DX,ARR,MFILE,IFMT)
C
C *** Double precision version
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6,CFORM*15
      DIMENSION DX(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='ORA1'
      IF (ICHECK.GE.997) CALL OTRC('ORA1  ','03/21/89')
C
      IER=0
      BFMT=IFMT.EQ.1
C
      IF (BFMT) THEN
       READ (MFILE,'(2A15,2I15)',ERR=99997,END=99997)
     *       ARR,CFORM,ITYPE,ILEN
      ELSE
       READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN,ILEN8,JRECL
      ENDIF
C
      IF (ITYPE.NE.1) GOTO 99997
C
      IF (BFMT) THEN
       READ(MFILE,CFORM,ERR=99997,END=99997) (DX(I),I=1,ILEN)
      ELSE
       CALL ORA0(DX,ILEN8,MFILE,JRECL)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(7,'ORA1  ')
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORA  ')
C
99999 END
C
C
C
      SUBROUTINE ORA2(VX,ARR,MFILE,IFMT)
C
C *** Double precision version
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6,CFORM*15
      DIMENSION VX(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='ORA2'
      IF (ICHECK.GE.997) CALL OTRC('ORA2  ','03/21/89')
C
      IER=0
      BFMT=IFMT.EQ.1
C
      IF (BFMT) THEN
       READ (MFILE,'(2A15,2I15)',ERR=99997,END=99997)
     *       ARR,CFORM,ITYPE,ILEN
      ELSE
       READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN,ILEN8,JRECL
      ENDIF
C
      IF (ITYPE.NE.2) GOTO 99997
C
      IF (BFMT) THEN
       READ(MFILE,CFORM,ERR=99997,END=99997) (VX(I),I=1,ILEN)
      ELSE
       CALL ORA0(VX,ILEN8,MFILE,JRECL)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(7,'ORA1  ')
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORA  ')
C
99999 END
C
C
C
      SUBROUTINE ORA3(KX,ARR,MFILE,IFMT)
C
C *** Integer version
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6,CFORM*15
      DIMENSION KX(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='ORA3'
      IF (ICHECK.GE.997) CALL OTRC('ORA3  ','03/21/89')
C
      IER=0
      BFMT=IFMT.EQ.1
C
      IF (BFMT) THEN
       READ (MFILE,'(2A15,2I15)',ERR=99997,END=99997)
     *       ARR,CFORM,ITYPE,ILEN
      ELSE
       READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN,ILEN8,JRECL
      ENDIF
C
      IF (ITYPE.NE.3) GOTO 99997
C
      IF (BFMT) THEN
       READ(MFILE,CFORM,ERR=99997,END=99997) (KX(I),I=1,ILEN)
      ELSE
       CALL ORA0(KX,ILEN8,MFILE,JRECL)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(7,'ORA1  ')
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORA  ')
C
99999 END
