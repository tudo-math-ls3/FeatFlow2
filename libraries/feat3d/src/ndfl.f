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
* NDFL                                                                 *
*                                                                      *
* Purpose  Determination of the number of degrees of freedom           *
*          per element                                                 *
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
* NDFL     I*4    Number of degrees of freedom per element             *
* IER      I*4    Error indicator                                      *
*                 -120 Wrong value of IELTYP                           *
*                                                                      *
************************************************************************
C
      FUNCTION NDFL(IELTYP)
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
      SUB='NDFL'
      IF (ICHECK.EQ.999) CALL OTRC('NDFL  ','07/22/91')
      IER=0
C
      GO TO (999,999,999,999,999,999,999,999,999,999,
     *       100,110,999,999,999,999,999,999,999,999,
     *       999,999,999,999,999,999,999,999,999,999,
     *       300,310,999,999,999,999,999,999,999,999), IELTYP+1
C
999   CONTINUE
      NDFL=0
      WRITE (CPARAM,'(I15)') IELTYP
      CALL WERR(-120,'NDFL  ')
      GOTO 99999
C
100   CONTINUE
      NDFL=1
      GOTO 99999
C
110   CONTINUE
      NDFL=NVE
      GOTO 99999
C
130   CONTINUE
      NDFL=NVE+NEE+NAE+1
      GOTO 99999
C
300   CONTINUE
310   CONTINUE
      NDFL=NAE
C
99999 END
