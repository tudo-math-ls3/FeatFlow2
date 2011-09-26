************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D  (Release 1.1)               *
*                                                                      *
* Authors: J. Harig, S. Turek                                          *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XSB0X                                                                *
*                                                                      *
* Purpose  Call of the following subdivision routines                  *
*                                                                      *
* Subroutines/functions called  XSBCA, XSB0, ZDISP                     *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ICHECK   I*4    =0 Skip check of subdivision                         *
* NFINE    I*4    Number of regular refinements (SB0)                  *
* ISCAD    I*4    =1 Determine array KADJ from coarse grid             *
* ISE      I*4    =1 Determine numbers of midpoints                    *
* ISA      I*4    =0 Release array KADJ on return                      *
*                    after determination of the new subdivision        *
* ISEEL     I*4   =1 Determine numbers of elements meeting at each     *
*                    edge                                              *
* ISAEL     I*4   =1 Determine numbers of elements meeting at each     *
*                    area                                              *
* ISVEL     I*4   =1 Determine numbers of elements meeting at each     *
*                    vertex                                            *
* IDISP    I*4    =1 Release free space on all arrays after completion *
*                                                                      *
* For the description of the remaining parameters see SB0              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE XSB0X(NFINE,ISCAD,ISE,ISA,ISVEL,ISEEL,ISAEL,
     *                 ISVED,ISAED,ISVAR,ISEAR,ISEVE,ISAVE,
     *                 ISVBD,ISEBD,ISABD,IDISP,PARX,PARY,PARZ,
     *                 SEDB,SADB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      EXTERNAL PARX,PARY,PARZ,SEDB,SADB
      SAVE /ERRCTL/,/TRIAA/
C
      IF (ICHECK.GE.997) CALL OTRC('XSB0X ','12/11/89')
C
      IF (ISCAD.GT.0) THEN
       CALL XSBCA(IDISP)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL XSB0(NFINE,ISE,ISA,ISVEL,ISEEL,ISAEL,ISVED,ISAED,
     *          ISVAR,ISEAR,ISEVE,ISAVE,ISVBD,ISEBD,ISABD,IDISP,
     *          PARX,PARY,PARZ,SEDB,SADB)
      IF (IER.NE.0) GOTO 99999
C
      IF ((ISE.EQ.0).AND.(LEDGE.GT.0)) THEN
       CALL ZDISP(0,LEDGE,'KEDGE ')
       LEDGE=0
      ENDIF
C
99999 END
