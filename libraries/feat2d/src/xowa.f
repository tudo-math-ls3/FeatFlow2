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
* XOWA                                                                 *
*                                                                      *
* Purpose  Store array                                                 *
*                                                                      *
* Subroutines/functions called  OF0, OWA0, ZLEN, ZLEN8, ZTYPE          *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----    ----                                                        *
* LNR     I*4    Number of array to be stored                          *
* ARR     C*6    Name of array                                         *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*(*)  Filename                                              *
*                CFILE.EQ.'SCRATCH' means temporary file               *
* IFMT    I*4    0  Format free input                                  *
*                1  Formatted input                                    *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* MFILE   I*4    Unit number used                                      *
* CFILE   C*(*)  Filename used (length not larger than for input)      *
* IER     I*4    Error indicator                                       *
*                -106  Value of  LNR  invalid                          *
*                -110  MFILE  exceeds range  16,...,80                 *
*                -112  File  CFILE  could not be opened                *
*                                                                      *
************************************************************************
C
      SUBROUTINE XOWA(LNR,ARR,MFILE,CFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6,CFILE*(*),CFORM*15
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XOWA'
      IF (ICHECK.GE.997) CALL OTRC('XOWA  ','12/11/89')
C
      IER=0
      BWARN=.FALSE.
      BFMT=IFMT.EQ.1
C
C *** Valid array number - 0 means that array does not yet exist
C
      IF ((LNR.LT.1).OR.(LNR.GT.NNARR)) THEN
       WRITE (CPARAM,'(A6)') ARR
       CALL WERR(-106,'XOWA  ')
       GOTO 99999
      ENDIF
C
C *** Open I/O-file
C
      CALL OF0(MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99999
C
      CALL ZTYPE(LNR,ITYPE)
      CALL ZLEN(LNR,ILEN)
      CALL ZLEN8(LNR,ILEN8)
C
      IF (BFMT) THEN
       CFORM=FMT(ITYPE)
       WRITE (MFILE,'(2A15,2I15)') ARR,CFORM,ITYPE,ILEN
      ELSE
       WRITE (MFILE) ARR,ITYPE,ILEN,ILEN8,IRECL8
      ENDIF
C
      L1=L(LNR)
      IF (BFMT) THEN
       GOTO (111,222,333),ITYPE
111    WRITE (MFILE,CFORM) (DWORK(L1+I),I=0,ILEN-1)
       GOTO 555
222    WRITE (MFILE,CFORM) (VWORK(L1+I),I=0,ILEN-1)
       GOTO 555
333    WRITE (MFILE,CFORM) (KWORK(L1+I),I=0,ILEN-1)
555    CONTINUE
      ELSE
       IF (ITYPE.GT.1) L1=(L1+1)/2
       CALL OWA0(DWORK(L1),ILEN8,MFILE,IRECL8)
      ENDIF
C
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(8,'XOWA  ')
C
99999 END
C
C
C
      SUBROUTINE OWA0(DX,N,MFILE,IRECL8)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('OWA0  ','12/11/89')
C
      IREC=N/IRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*IRECL8
1     WRITE (MFILE) (DX(J1+I),I=1,IRECL8)
      IF (MOD(N,IRECL8).EQ.0) GOTO 99999
      J1=IREC*IRECL8
      WRITE (MFILE) (DX(J1+I),I=1,MOD(N,IRECL8))
C
99999 END
