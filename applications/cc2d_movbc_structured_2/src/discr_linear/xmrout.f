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
* ORALCn                                                               *
*                                                                      *
* Purpose  Get arrays from file, previously stored by XOWA or OWAn     *
*                                                                      *
* Subroutines/functions called  ORALC0                                 *
*                                                                      *
* Version from  10/11/94                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----    ----                                                        *
* MFILE   I*4    I/O unit - previouly opened (by OF0)                  *
* IFMT    I*4    0  Format free input                                  *
*                1  Formatted input                                    *
* A1      R*8    Scaling factor: DX=DX+A1*DX(readin)                   *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* DX      R*8    Double precision array - version ORA1                 *
* VX      R*4    Single precision array - version ORA2                 *
* ARR     C*6    Name of array read from data file                     *
* IER     I*4    Error indicator                                       *
*                -110  Error while reading from unit  MFILE            *
*                -113  Wrong value of  ITYPE  read from  MFILE         *
*                                                                      *
************************************************************************
C
      SUBROUTINE ORALC1(DX,A1,ARR,MFILE,IFMT)

C *** Double precision version
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

C parameters
      
      CHARACTER ARR*6
      DOUBLE PRECISION DX(*)
      DOUBLE PRECISION A1
      INTEGER MFILE,IFMT
      
C local variables

      INTEGER ITYPE, ILEN, ILEN8, JRECL, LNR

      SUB='ORALC1'
      IF (ICHECK.GE.997) CALL OTRC('ORALC1','10/11/94')
C
      IER=0
C
      IF (IFMT.NE.0) GOTO 99997
C
      READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN,ILEN8,JRECL
C
      IF (ITYPE.NE.1) GOTO 99997
C
      CALL ORALCD(DX,A1,ILEN8,MFILE,JRECL)
      IF (IER.NE.0) GOTO 99999
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(7,'ORALC1')
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORALC')
C
99999 END
C
C
C
      SUBROUTINE ORALC2(VX,A1,ARR,MFILE,IFMT)

C *** Single precision version
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

C parameters
      
      CHARACTER ARR*6
      REAL VX(*)
      DOUBLE PRECISION A1
      INTEGER MFILE,IFMT
      
C local variables

      INTEGER ITYPE, ILEN, ILEN8, JRECL, LNR
C
      SUB='ORALC2'
      IF (ICHECK.GE.997) CALL OTRC('ORALC2','10/11/94')
C
      IER=0
C
      IF (IFMT.NE.0) GOTO 99997
C
      READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN,ILEN8,JRECL
C
      IF (ITYPE.NE.2) GOTO 99997
C
      CALL ORALCV(VX,A1,ILEN8,MFILE,JRECL)
      IF (IER.NE.0) GOTO 99999
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(7,'ORALC2')
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORALC')

99999 END


      SUBROUTINE ORALCD(DX,A1,N,MFILE,JRECL8)
      
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'

C parameters
      
      INTEGER N
      DOUBLE PRECISION DX(N)
      DOUBLE PRECISION A1
      INTEGER MFILE,JRECL8
      
C local variables

      INTEGER IREC,JREC,I,J1
      DOUBLE PRECISION DXH(512)

      IF (ICHECK.EQ.999) CALL OTRC('ORALCD','10/11/94')
C
      IREC=N/JRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*JRECL8
      READ (MFILE,ERR=99998,END=99998) (DXH(I),I=1,JRECL8)
      DO 1 I=1,JRECL8
      DX(J1+I)=DX(J1+I)+DXH(I)
1     CONTINUE
      IF (MOD(N,JRECL8).EQ.0) GOTO 99999
      J1=IREC*JRECL8
      READ (MFILE,ERR=99998,END=99998) (DXH(I),I=1,MOD(N,JRECL8))
      DO 2 I=1,MOD(N,JRECL8)
      DX(J1+I)=DX(J1+I)+DXH(I)
2     CONTINUE
      GOTO 99999
C
99998 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'ORALCD')

99999 END


      SUBROUTINE ORALCV(VX,A1,N,MFILE,JRECL8)
      
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'

C parameters
      
      INTEGER N
      REAL VX(N)
      DOUBLE PRECISION A1
      INTEGER MFILE,JRECL8
      
C local variables

      INTEGER IREC,JREC,I,J1
      REAL VXH(512)

      IF (ICHECK.EQ.999) CALL OTRC('ORALCV','10/11/94')
C
      IREC=N/JRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*JRECL8
      READ (MFILE,ERR=99998,END=99998) (VXH(I),I=1,JRECL8)
      DO 1 I=1,JRECL8
      VX(J1+I)=VX(J1+I)+VXH(I)
1     CONTINUE
      IF (MOD(N,JRECL8).EQ.0) GOTO 99999
      J1=IREC*JRECL8
      READ (MFILE,ERR=99998,END=99998) (VXH(I),I=1,MOD(N,JRECL8))
      DO 2 I=1,MOD(N,JRECL8)
      VX(J1+I)=VX(J1+I)+VXH(I)
2     CONTINUE
      GOTO 99999
C
99998 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'ORALCV')
C
99999 END
C
C
C
      SUBROUTINE OWA2V(VX,ARR,NX,MFILE,IFMT)

      IMPLICIT NONE

      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
C parameter

      INTEGER NX, MFILE, IFMT
      REAL VX(NX)
      CHARACTER ARR*6

C local variables

      CHARACTER CFORM*15
      INTEGER I,LNR
      LOGICAL BFMT

      SUB='OWA2V'
      IF (ICHECK.GE.997) CALL OTRC('OWA2V ','10/11/94')
C
      IER=0
      BFMT=IFMT.EQ.1
C
      IF (BFMT) THEN
       CFORM=FMT(1)
       WRITE (MFILE,'(2A15,2I15)') ARR,CFORM,1,NX
      ELSE
       WRITE (MFILE) ARR,2,NX,NX,IRECL8
      ENDIF
C
      IF (BFMT) THEN
       WRITE (MFILE,CFORM) (VX(I),I=1,NX)
      ELSE
       CALL OWA0V(VX,NX,MFILE,IRECL8)
      ENDIF
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(8,'XOWA2V')
C
99999 END


      SUBROUTINE OWA0V(VX,N,MFILE,IRECL8)
      
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      
C parameter
      
      INTEGER N,MFILE,IRECL8
      REAL VX(N)

C local variables

      INTEGER IREC,JREC,I,J1

      IF (ICHECK.EQ.999) CALL OTRC('OWA0V ','10/11/94')
C
      IREC=N/IRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*IRECL8
1     WRITE (MFILE) (VX(J1+I),I=1,IRECL8)
      IF (MOD(N,IRECL8).EQ.0) GOTO 99999
      J1=IREC*IRECL8
      WRITE (MFILE) (VX(J1+I),I=1,MOD(N,IRECL8))

99999 END
