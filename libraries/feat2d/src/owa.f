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
* OWAn                                                                 *
*                                                                      *
* Purpose  Write array to file                                         *
*                                                                      *
* Subroutines/functions called  OWA0                                   *
*                                                                      *
* Version from  03/21/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----    ----                                                        *
* DX      R*8    Double precision array - version OWA1                 *
* VX      R*4    Single precision array - version OWA2                 *
* KX      I*4    Integer array          - version OWA3                 *
* MFILE   I*4    I/O unit - previouly opened (by OF0)                  *
* IFMT    I*4    0  Format free output                                 *
*                1  Formatted output                                   *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
*                                                                      *
************************************************************************
C
      SUBROUTINE OWA1(DX,ARR,NX,MFILE,IFMT)
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
      SUB='OWA1'
      IF (ICHECK.GE.997) CALL OTRC('OWA1  ','03/21/89')
C
      IER=0
      BFMT=IFMT.EQ.1
C
      IF (BFMT) THEN
       CFORM=FMT(1)
       WRITE (MFILE,'(2A15,2I15)') ARR,CFORM,1,NX
      ELSE
       WRITE (MFILE) ARR,1,NX,NX,IRECL8
      ENDIF
C
      IF (BFMT) THEN
       WRITE (MFILE,CFORM) (DX(I),I=1,NX)
      ELSE
       CALL OWA0(DX,NX,MFILE,IRECL8)
      ENDIF
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(8,'XOWA  ')
C
99999 END
C
C
C
      SUBROUTINE OWA2(VX,ARR,NX,MFILE,IFMT)
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
      SUB='OWA2'
      IF (ICHECK.GE.997) CALL OTRC('OWA2  ','03/21/89')
C
      IER=0
      BFMT=IFMT.EQ.1
C
      IF (BFMT) THEN
       CFORM=FMT(2)
       WRITE (MFILE,'(2A15,2I15)') ARR,CFORM,2,NX
      ELSE
       ILEN8=(NX+1)/2
       WRITE (MFILE) ARR,2,NX,ILEN8,IRECL8
      ENDIF
C
      IF (BFMT) THEN
       WRITE (MFILE,CFORM) (VX(I),I=1,NX)
      ELSE
       CALL OWA0(VX,ILEN8,MFILE,IRECL8)
      ENDIF
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(8,'XOWA  ')
C
99999 END
C
C
C
      SUBROUTINE OWA3(KX,ARR,NX,MFILE,IFMT)
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
      SUB='OWA3'
      IF (ICHECK.GE.997) CALL OTRC('OWA3  ','03/21/89')
C
      IER=0
      BFMT=IFMT.EQ.1
C
      IF (BFMT) THEN
       CFORM=FMT(3)
       WRITE (MFILE,'(2A15,2I15)') ARR,CFORM,3,NX
      ELSE
       ILEN8=(NX+1)/2
       WRITE (MFILE) ARR,3,NX,ILEN8,IRECL8
      ENDIF
C
      IF (BFMT) THEN
       WRITE (MFILE,CFORM) (KX(I),I=1,NX)
      ELSE
       CALL OWA0(KX,ILEN8,MFILE,IRECL8)
      ENDIF
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(8,'XOWA  ')
C
99999 END
