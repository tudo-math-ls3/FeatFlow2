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
* NDFG                                                                 *
*                                                                      *
* Purpose  Determination of the total number of degrees of freedom     *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  07/22/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IELTYP   I*4    Type of element corresponding to  Ennn               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* NDFG     I*4    Total number of degrees of freedom                   *
* IER      I*4    Error indicator                                      *
*                 -120 Wrong value of IELTYP                           *
*                                                                      *
************************************************************************
C
      FUNCTION NDFG(IELTYP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='NDFG'
      IF (ICHECK.EQ.999) CALL OTRC('NDFG  ','07/22/91')
C
      IER=0
      GOTO  (999,999,999,999,999,999,999,999,999,999,
     *       100,110,999,999,999,999,999,999,999,999,
     *       999,999,999,999,999,999,999,999,999,999,
     *       300,310,999,999,999,999,999,999,999,999), IELTYP+1
C
999   CONTINUE
      NDFG=0
      WRITE (CPARAM,'(I15)') IELTYP
      CALL WERR(-120,'NDFG  ')
      GOTO 99999
C
100   CONTINUE
      NDFG=NEL
      GOTO 99999
C
110   CONTINUE
      NDFG=NVT
      GOTO 99999
C
130   CONTINUE
      NDFG=NVT+NET+NAT+NEL
      GOTO 99999
C
300   CONTINUE
310   CONTINUE
      NDFG=NAT
C
99999 END
