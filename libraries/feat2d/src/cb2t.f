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
* CB2T                                                                 *
*                                                                      *
* Purpose  Set coordinates and weights for numerical integration       *
*          over a triangle                                             *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  03/21/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ICUB   I*4    Number of desired cubature formula                     *
*               ICUB    NCUBP    DEGREE                                *
*                 1       1        2       1x1 Gauss formula           *
*                 2       3        2       Trapezoidal rule            *
*                 3       3        3       3-point Gauss               *
*                 4       3        3       Collatz formula             *
*                 5       6        4       As 4 plus midpoints of edge *
*                 6       7        4       vertices, midpoints, center *
*                 7                        Not yet implemented         *
*                 8       7        6                                   *
*                 9      12        7                                   *
*                10                        Not yet implemented         *
*                11       4        2       Piecewise 1x1 Gauss formula *
*                12       6        2       Piecewise trapezoidal rule  *
*                13       9        3       Piecewise 3-point Gauss     *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DXI      R*8    Barycentric coordinates of cubature points           *
* DOMEGA   R*8    Weights                                              *
* NCUBP    I*4    Number of cubature points                            *
* IER      I*4    Error indicator                                      *
*                 -119 Wrong value of ICUB                             *
*                                                                      *
************************************************************************
C
      SUBROUTINE CB2T (ICUB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNCUBP=36)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/CUB/
C
      SUB='CB2T'
      IF (ICHECK.GE.997) CALL OTRC('CB2T  ','03/21/89')
C
      IER=0
C
      GOTO (10,20,30,40,50,60,999,80,90,999,110,120,130) ICUB
C
999   WRITE (CPARAM,'(I15)') ICUB
      CALL WERR(-119,'CB2T  ')
      GOTO 99999
C
C
C
10    DO 11 I=1,3
11    DXI(1,I)=.3333333333333333D0
      DOMEGA(1)=0.5D0
      NCUBP=1
      GOTO 99999
C
C
C
20    DXI(1,1)=1D0
      DXI(1,2)=0D0
      DXI(1,3)=0D0
      DXI(2,1)=0D0
      DXI(2,2)=1D0
      DXI(2,3)=0D0
      DXI(3,1)=0D0
      DXI(3,2)=0D0
      DXI(3,3)=1D0
      DO 21 I=1,3
21    DOMEGA(I)=.1666666666666667D0
      NCUBP=3
      GOTO 99999
C
C
C
30    DXI(1,1)=.5D0
      DXI(1,2)=.5D0
      DXI(1,3)= 0D0
      DXI(2,1)= 0D0
      DXI(2,2)=.5D0
      DXI(2,3)=.5D0
      DXI(3,1)=.5D0
      DXI(3,2)= 0D0
      DXI(3,3)=.5D0
      DO 31 I=1,3
31    DOMEGA(I)=.1666666666666667D0
      NCUBP=3
      GOTO 99999
C
C
C
40    DXI(1,1)=.6666666666666667D0
      DXI(1,2)=.1666666666666667D0
      DXI(1,3)=.1666666666666667D0
      DXI(2,1)=.1666666666666667D0
      DXI(2,2)=.6666666666666667D0
      DXI(2,3)=.1666666666666667D0
      DXI(3,1)=.1666666666666667D0
      DXI(3,2)=.1666666666666667D0
      DXI(3,3)=.6666666666666667D0
      DO 41 I=1,3
41    DOMEGA(I)=.1666666666666667D0
      NCUBP=3
      GOTO 99999
C
C
C
50    DXI(1,1)=.5D0
      DXI(1,2)=.5D0
      DXI(1,3)= 0D0
      DXI(2,1)= 0D0
      DXI(2,2)=.5D0
      DXI(2,3)=.5D0
      DXI(3,1)=.5D0
      DXI(3,2)= 0D0
      DXI(3,3)=.5D0
C
      DXI(4,1)=.6666666666666667D0
      DXI(4,2)=.1666666666666667D0
      DXI(4,3)=.1666666666666667D0
      DXI(5,1)=.1666666666666667D0
      DXI(5,2)=.6666666666666667D0
      DXI(5,3)=.1666666666666667D0
      DXI(6,1)=.1666666666666667D0
      DXI(6,2)=.1666666666666667D0
      DXI(6,3)=.6666666666666667D0
C
      DO 51 I=1,3
      DOMEGA(I)  =.1666666666666667D-1
51    DOMEGA(I+3)=.15D0
      NCUBP=6
      GOTO 99999
C
C
C
60    DO 61 I=1,3
61    DXI(1,I)=.3333333333333333D0
C
      DXI(2,1)=.5D0
      DXI(2,2)=.5D0
      DXI(2,3)= 0D0
      DXI(3,1)= 0D0
      DXI(3,2)=.5D0
      DXI(3,3)=.5D0
      DXI(4,1)=.5D0
      DXI(4,2)= 0D0
      DXI(4,3)=.5D0
C
      DXI(5,1)=1D0
      DXI(5,2)=0D0
      DXI(5,3)=0D0
      DXI(6,1)=0D0
      DXI(6,2)=1D0
      DXI(6,3)=0D0
      DXI(7,1)=0D0
      DXI(7,2)=0D0
      DXI(7,3)=1D0
C
      DOMEGA(1)=.225D0
      DO 62 I=2,4
      DOMEGA(I)  =.6666666666666667D-1
62    DOMEGA(I+3)=.25D-01
C
      NCUBP=7
      GOTO 99999
C
C
C
80    DO 81 I=1,3
81    DXI(1,I)=.3333333333333333D0
C
      DXI(2,1)=.5971587178976983D-1
      DXI(2,2)=.4701420641051151D0
      DXI(2,3)=.4701420641051151D0
      DXI(3,1)=.4701420641051151D0
      DXI(3,2)=.5971587178976983D-1
      DXI(3,3)=.4701420641051151D0
      DXI(4,1)=.4701420641051151D0
      DXI(4,2)=.4701420641051151D0
      DXI(4,3)=.5971587178976983D-1
C
      DXI(5,1)=.7974269853530872D0
      DXI(5,2)=.1012865073234563D0
      DXI(5,3)=.1012865073234563D0
      DXI(6,1)=.1012865073234563D0
      DXI(6,2)=.7974269853530872D0
      DXI(6,3)=.1012865073234563D0
      DXI(7,1)=.1012865073234563D0
      DXI(7,2)=.1012865073234563D0
      DXI(7,3)=.7974269853530872D0
C
      DOMEGA(1)=.1125D0
      DO 82 I=2,4
      DOMEGA(I)   =.661970763942531D-1
82    DOMEGA(I+3) =.6296959027241355D-1
C
      NCUBP=7
      GOTO 99999
C
C
C
90    DXI(1,1)=.27771616697639178D+00
      DXI(1,2)=.51584233435359177D+00
      DXI(1,3)=.20644149867001643D+00
      DXI(2,2)=.27771616697639178D+00
      DXI(2,3)=.51584233435359177D+00
      DXI(2,1)=.20644149867001643D+00
      DXI(3,3)=.27771616697639178D+00
      DXI(3,1)=.51584233435359177D+00
      DXI(3,2)=.20644149867001643D+00
C
      DXI(4,1)=.32150249385198182D+00
      DXI(4,2)=.55225456656926611D-01
      DXI(4,3)=.62327204949109156D+00
      DXI(5,2)=.32150249385198182D+00
      DXI(5,3)=.55225456656926611D-01
      DXI(5,1)=.62327204949109156D+00
      DXI(6,3)=.32150249385198182D+00
      DXI(6,1)=.55225456656926611D-01
      DXI(6,2)=.62327204949109156D+00
C
      DXI(7,1)=.30472650086816719D+00
      DXI(7,2)=.66094919618673565D+00
      DXI(7,3)=.34324302945097146D-01
      DXI(8,2)=.30472650086816719D+00
      DXI(8,3)=.66094919618673565D+00
      DXI(8,1)=.34324302945097146D-01
      DXI(9,3)=.30472650086816719D+00
      DXI(9,1)=.66094919618673565D+00
      DXI(9,2)=.34324302945097146D-01
C
      DXI(10,1)=.67517867073916085D-01
      DXI(10,2)=.62382265094402118D-01
      DXI(10,3)=.87009986783168179D+00
      DXI(11,2)=.67517867073916085D-01
      DXI(11,3)=.62382265094402118D-01
      DXI(11,1)=.87009986783168179D+00
      DXI(12,3)=.67517867073916085D-01
      DXI(12,1)=.62382265094402118D-01
      DXI(12,2)=.87009986783168179D+00
C
      DO 91 I=1,3
      DOMEGA(I)=  .67493187009802774D-01
      DOMEGA(I+3)=.43881408714446055D-01
      DOMEGA(I+6)=.28775042784981585D-01
91    DOMEGA(I+9)=.26517028157436251D-01
C
      NCUBP=12
      GOTO 99999
C
C
C
110   DO 111 I=1,3
111   DXI(1,I)=0.3333333333333333D0
C
      DXI(2,1)=.66666666666666667D0
      DXI(2,2)=.16666666666666667D0
      DXI(2,3)=.16666666666666667D0
      DXI(3,1)=.16666666666666667D0
      DXI(3,2)=.66666666666666667D0
      DXI(3,3)=.16666666666666667D0
      DXI(4,1)=.16666666666666667D0
      DXI(4,2)=.16666666666666667D0
      DXI(4,3)=.66666666666666667D0
C
      DO 113 I=1,4
113   DOMEGA(I)=0.125D0
C
      NCUBP=4
      GOTO 99999
C
C
C
120   DXI(1,1)=.5D0
      DXI(1,2)=.5D0
      DXI(1,3)= 0D0
      DXI(2,1)= 0D0
      DXI(2,2)=.5D0
      DXI(2,3)=.5D0
      DXI(3,1)=.5D0
      DXI(3,2)= 0D0
      DXI(3,3)=.5D0
C
      DXI(4,1)= 1D0
      DXI(4,2)= 0D0
      DXI(4,3)= 0D0
      DXI(5,1)=.5D0
      DXI(5,2)=.5D0
      DXI(5,3)= 0D0
      DXI(6,1)=.5D0
      DXI(6,2)= 0D0
      DXI(6,3)=.5D0
C
      DXI(7,1)= 0D0
      DXI(7,2)= 1D0
      DXI(7,3)= 0D0
      DXI(8,1)= 0D0
      DXI(8,2)=.5D0
      DXI(8,3)=.5D0
      DXI(9,1)=.5D0
      DXI(9,2)=.5D0
      DXI(9,3)= 0D0
C
      DXI(10,1)= 0D0
      DXI(10,2)= 0D0
      DXI(10,3)= 1D0
      DXI(11,1)=.5D0
      DXI(11,2)= 0D0
      DXI(11,3)=.5D0
      DXI(12,1)= 0D0
      DXI(12,2)=.5D0
      DXI(12,3)=.5D0
C
      DO 121 I=1,12
121   DOMEGA(I)=.4166666666666667D-1
      NCUBP=12
      GOTO 99999
C
C
C
130   DXI(1,1)= .5D0
      DXI(1,2)=.25D0
      DXI(1,3)=.25D0
      DXI(2,1)=.25D0
      DXI(2,2)= .5D0
      DXI(2,3)=.25D0
      DXI(3,1)=.25D0
      DXI(3,2)=.25D0
      DXI(3,3)= .5D0
C
      DXI(4,1)=.75D0
      DXI(4,2)=.25D0
      DXI(4,3)=  0D0
      DXI(5,1)= .5D0
      DXI(5,2)=.25D0
      DXI(5,3)=.25D0
      DXI(6,1)=.75D0
      DXI(6,2)=  0D0
      DXI(6,3)=.25D0
C
      DXI(7,1)=  0D0
      DXI(7,2)=.75D0
      DXI(7,3)=.25D0
      DXI(8,1)=.25D0
      DXI(8,2)= .5D0
      DXI(8,3)=.25D0
      DXI(9,1)=.25D0
      DXI(9,2)=.75D0
      DXI(9,3)= 0D0
C
      DXI(10,1)=.25D0
      DXI(10,2)=  0D0
      DXI(10,3)=.75D0
      DXI(11,1)=.25D0
      DXI(11,2)=.25D0
      DXI(11,3)=  5D0
      DXI(12,1)=  0D0
      DXI(12,2)=.25D0
      DXI(12,3)=.75D0
C
      DO 131 I=1,12
131   DOMEGA(I)=.41666666666666667D-1
      NCUBP=12
C
99999 END
