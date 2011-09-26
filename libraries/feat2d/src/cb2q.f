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
* CB2Q                                                                 *
*                                                                      *
* Purpose  Set coordinates and weights for numerical integration       *
*          over the unit square [-1,1]x[-1,1]                          *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  04/19/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ICUB     I*4    Number of cubature formula                           *
*                 ICUB    NCUB    DEGREE                               *
*                   1       1        2      1x1 Gaussian formula       *
*                   2       3        2      Trapezoidal rule           *
*                   3       4        1      Midpoints of edges         *
*                   4       4        4      2x2 Gaussian formula       *
*                   5       4        4      Not symmetric              *
*                   6       6        5      Not symmetric              *
*                   7       7        6      Not symmetric              *
*                   8       9        6      3x3 Gaussian formula       *
*                   9      12        7      Gaussian formula           *
*                  10      16        8      4x4 Gaussian formula       *
*                  11      25       10      5x5 Gaussian formula       *
*                  12       4        2      Piecewise 1x1 Gauss formula*
*                  13       9        2      Piecewise trapezoidal rule *
*                  14      16        4      Piecewise 2x2 Gauss formula*
*                  15      36        6      Piecewise 3x3 Gauss formula*
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DXI      R*8    Cartesian coordinates of cubature points             *
* DOMEGA   R*8    Weights                                              *
* NCUBP    I*4    Number of cubature points                            *
* IER      I*4    Error indicator                                      *
*                 -119 Wrong value of ICUB                             *
*                                                                      *
************************************************************************
C
      SUBROUTINE CB2Q (ICUB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNCUBP=36)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/CUB/
C
      SUB='CB2Q'
      IF (ICHECK.GE.997) CALL OTRC('CB2Q  ','04/19/89')
C
      IER=0
C
      GOTO (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150) ICUB
C
      WRITE (CPARAM,'(I15)') ICUB
      CALL WERR(-119,'CB2Q  ')
      GOTO 99999
C
C
C
10    DXI(1,1)= 0D0
      DXI(1,2)= 0D0
      DOMEGA(1)=4D0
      NCUBP=1
      GOTO 99999
C
C
C
20    DXI(1,1)=-1D0
      DXI(1,2)=-1D0
      DXI(2,1)= 1D0
      DXI(2,2)=-1D0
      DXI(3,1)= 1D0
      DXI(3,2)= 1D0
      DXI(4,1)=-1D0
      DXI(4,2)= 1D0
      DO 21 I=1,4
21    DOMEGA(I)=1D0
      NCUBP=4
      GOTO 99999
C
C
C
30    DXI(1,1)= 0D0
      DXI(1,2)=-1D0
      DXI(2,1)= 1D0
      DXI(2,2)= 0D0
      DXI(3,1)= 0D0
      DXI(3,2)= 1D0
      DXI(4,1)=-1D0
      DXI(4,2)= 0D0
      DO 31 I=1,4
31    DOMEGA(I)=1D0
      NCUBP=4
      GOTO 99999
C
C
C
40    DXI(1,1)= 0.577350269189626D0
      DXI(1,2)= 0.577350269189626D0
      DXI(2,1)=-0.577350269189626D0
      DXI(2,2)= 0.577350269189626D0
      DXI(3,1)=-0.577350269189626D0
      DXI(3,2)=-0.577350269189626D0
      DXI(4,1)= 0.577350269189626D0
      DXI(4,2)=-0.577350269189626D0
      DO 41 I=1,4
41    DOMEGA(I)=1D0
      NCUBP=4
      GOTO 99999
C
C
C
50    DXI(1,1)= 0.816496580927726D0
      DXI(1,2)= 0D0
      DXI(2,1)=-0.816496580927726D0
      DXI(2,2)= 0D0
      DXI(3,1)= 0D0
      DXI(3,2)=-0.816496580927726D0
      DXI(4,1)= 0D0
      DXI(4,2)= 0.816496580927726D0
      DO 51 I=1,4
51    DOMEGA(I)=1D0
      NCUBP=4
      GOTO 99999
C
C
C
60    DXI(1,1)= 0D0
      DXI(1,2)= 0.934172358962716D0
      DXI(2,1)= 0D0
      DXI(2,2)=-0.356822089773090D0
      DXI(3,1)= 0.774596669241483D0
      DXI(3,2)= 0.390885162530070D0
      DXI(4,1)=-0.774596669241483D0
      DXI(4,2)= 0.390885162530070D0
      DXI(5,1)= 0.774596669241483D0
      DXI(5,2)=-0.852765377881771D0
      DXI(6,1)=-0.774596669241483D0
      DXI(6,2)=-0.852765377881771D0
C
      DOMEGA(1)=0.491365692888926D0
      DOMEGA(2)=0.128641208488885D1
      DOMEGA(3)=0.761883709085613D0
      DOMEGA(4)=0.761883709085613D0
      DOMEGA(5)=0.349227402025498D0
      DOMEGA(6)=0.349227402025498D0
      NCUBP=6
      GOTO 99999
C
C
C
70    DXI(1,1)= 0D0
      DXI(1,2)= 0D0
      DXI(2,1)= 0D0
      DXI(2,2)= 0.966091783079296D0
      DXI(3,1)= 0D0
      DXI(3,2)=-0.966091783079296D0
      DXI(4,1)= 0.774596669241483D0
      DXI(4,2)= 0.774596669241483D0
      DXI(5,1)=-0.774596669241483D0
      DXI(5,2)= 0.774596669241483D0
      DXI(6,1)= 0.774596669241483D0
      DXI(6,2)=-0.774596669241483D0
      DXI(7,1)=-0.774596669241483D0
      DXI(7,2)=-0.774596669241483D0
C
      DOMEGA(1)=1.142857142857143D0
      DOMEGA(2)=0.317460317460317D0
      DOMEGA(3)=0.317460317460317D0
      DO 71 I=4,7
71    DOMEGA(I)=0.555555555555556D0
      NCUBP=7
      GOTO 99999
C
C
C
80    DXI(1,1)= 0.774596669241483D0
      DXI(1,2)= 0.774596669241483D0
      DXI(2,1)=-0.774596669241483D0
      DXI(2,2)= 0.774596669241483D0
      DXI(3,1)= 0.774596669241483D0
      DXI(3,2)=-0.774596669241483D0
      DXI(4,1)=-0.774596669241483D0
      DXI(4,2)=-0.774596669241483D0
      DXI(5,1)= 0.774596669241483D0
      DXI(5,2)= 0D0
      DXI(6,1)=-0.774596669241483D0
      DXI(6,2)= 0D0
      DXI(7,1)= 0D0
      DXI(7,2)= 0.774596669241483D0
      DXI(8,1)= 0D0
      DXI(8,2)=-0.774596669241483D0
      DXI(9,1)= 0D0
      DXI(9,2)= 0D0
C
      DO 81 I=1,4
81    DOMEGA(I)=0.308641975308642D0
      DO 82 I=5,8
82    DOMEGA(I)=0.493827160493827D0
      DOMEGA(9)=0.790123456790123D0
      NCUBP=9
      GOTO 99999
C
C
C
90    DXI(1,1)= .92582009977255146D0
      DXI(1,2)=0D0
      DXI(2,1)=-.92582009977255146D0
      DXI(2,2)=0D0
      DXI(3,1)=0D0
      DXI(3,2)= .92582009977255146D0
      DXI(4,1)=0D0
      DXI(4,2)=-.92582009977255146D0
C
      DXI(5,1)= .38055443320831566D0
      DXI(5,2)= .38055443320831566D0
      DXI(6,1)= .38055443320831566D0
      DXI(6,2)=-.38055443320831566D0
      DXI(7,1)=-.38055443320831566D0
      DXI(7,2)= .38055443320831566D0
      DXI(8,1)=-.38055443320831566D0
      DXI(8,2)=-.38055443320831566D0
C
      DXI(9,1)=  .80597978291859874D0
      DXI(9,2)=  .80597978291859874D0
      DXI(10,1)= .80597978291859874D0
      DXI(10,2)=-.80597978291859874D0
      DXI(11,1)=-.80597978291859874D0
      DXI(11,2)= .80597978291859874D0
      DXI(12,1)=-.80597978291859874D0
      DXI(12,2)=-.80597978291859874D0
C
      DO 91 I=1,4
91    DOMEGA(I)=.24197530864197531D+00
      DO 92 I=5,8
92    DOMEGA(I)=.52059291666739446D0
      DO 93 I=9,12
93    DOMEGA(I)=.23743177469063023D0
      NCUBP=12
      GOTO 99999
C
C
C
100   DXI(1,1)= 0.861136311594053003D0
      DXI(1,2)= 0.339981043584855994D0
      DXI(2,1)=-0.861136311594053003D0
      DXI(2,2)=-0.339981043584855994D0
      DXI(3,1)=-0.861136311594053003D0
      DXI(3,2)= 0.339981043584855994D0
      DXI(4,1)= 0.861136311594053003D0
      DXI(4,2)=-0.339981043584855994D0
      DXI(5,1)=-0.339981043584855994D0
      DXI(5,2)=-0.861136311594053003D0
      DXI(6,1)= 0.339981043584855994D0
      DXI(6,2)= 0.861136311594053003D0
      DXI(7,1)= 0.339981043584855994D0
      DXI(7,2)=-0.861136311594053003D0
      DXI(8,1)=-0.339981043584855994D0
      DXI(8,2)= 0.861136311594053003D0
C
      DXI(9,1)= -0.339981043584855994D0
      DXI(9,2)=  0.339981043584855994D0
      DXI(10,1)= 0.339981043584855994D0
      DXI(10,2)=-0.339981043584855994D0
      DXI(11,1)= 0.339981043584855994D0
      DXI(11,2)= 0.339981043584855994D0
      DXI(12,1)=-0.339981043584855994D0
      DXI(12,2)=-0.339981043584855994D0
C
      DXI(13,1)= 0.861136311594053003D0
      DXI(13,2)=-0.861136311594053003D0
      DXI(14,1)=-0.861136311594053003D0
      DXI(14,2)= 0.861136311594053003D0
      DXI(15,1)=-0.861136311594053003D0
      DXI(15,2)=-0.861136311594053003D0
      DXI(16,1)= 0.861136311594053003D0
      DXI(16,2)= 0.861136311594053003D0
C
      DO 101 I=1,8
101   DOMEGA(I)=0.226851851851851888D0
      DO 102 I=9,12
      DOMEGA(I)  =0.425293303010694096D0
102   DOMEGA(I+4)=0.121002993285602101D0
C
      NCUBP=16
      GOTO 99999
C
C
C
110   DXI(1,1)= 0.906179845938664005D0
      DXI(1,2)= 0.906179845938664005D0
      DXI(2,1)= 0.906179845938664005D0
      DXI(2,2)=-0.906179845938664005D0
      DXI(3,1)=-0.906179845938664005D0
      DXI(3,2)= 0.906179845938664005D0
      DXI(4,1)=-0.906179845938664005D0
      DXI(4,2)=-0.906179845938664005D0
C
      DXI(5,1)= -0.906179845938664005D0
      DXI(5,2)= -0.538469310105682997D0
      DXI(6,1)= -0.538469310105682997D0
      DXI(6,2)= -0.906179845938664005D0
      DXI(7,1)=  0.906179845938664005D0
      DXI(7,2)=  0.538469310105682997D0
      DXI(8,1)= -0.906179845938664005D0
      DXI(8,2)=  0.538469310105682997D0
      DXI(9,1)=  0.538469310105682997D0
      DXI(9,2)= -0.906179845938664005D0
      DXI(10,1)=-0.538469310105682997D0
      DXI(10,2)= 0.906179845938664005D0
      DXI(11,1)= 0.538469310105682997D0
      DXI(11,2)= 0.906179845938664005D0
      DXI(12,1)= 0.906179845938664005D0
      DXI(12,2)=-0.538469310105682997D0
C
      DXI(13,1)= 0.538469310105682997D0
      DXI(13,2)= 0D0
      DXI(14,1)= 0D0
      DXI(14,2)= 0.538469310105682997D0
      DXI(15,1)=-0.538469310105682997D0
      DXI(15,2)= 0D0
      DXI(16,1)= 0D0
      DXI(16,2)=-0.538469310105682997D0
C
      DXI(17,1)= 0.538469310105682997D0
      DXI(17,2)= 0.538469310105682997D0
      DXI(18,1)= 0.538469310105682997D0
      DXI(18,2)=-0.538469310105682997D0
      DXI(19,1)=-0.538469310105682997D0
      DXI(19,2)= 0.538469310105682997D0
      DXI(20,1)=-0.538469310105682997D0
      DXI(20,2)=-0.538469310105682997D0
C
      DXI(21,1)= 0.906179845938664005D0
      DXI(21,2)= 0D0
      DXI(22,1)= 0D0
      DXI(22,2)= 0.906179845938664005D0
      DXI(23,1)=-0.906179845938664005D0
      DXI(23,2)= 0D0
      DXI(24,1)= 0D0
      DXI(24,2)=-0.906179845938664005D0
C
      DXI(25,1)= 0D0
      DXI(25,2)= 0D0
C
      DO 111 I=1,4
      DOMEGA(I)   = 0.561343488624285944D-1
      DOMEGA(I+4) = 0.113399999999999834D0
      DOMEGA(I+8) = 0.113399999999999834D0
      DOMEGA(I+12)= 0.272286532550750485D0
      DOMEGA(I+16)= 0.229085404223990666D0
111   DOMEGA(I+20)= 0.134785072387520868D0
      DOMEGA(25)= 0.323634567901234682D0
C
      NCUBP=25
      GOTO 99999
C
C
C
120   DXI(1,1)= 0.5D0
      DXI(1,2)= 0.5D0
      DXI(2,1)= 0.5D0
      DXI(2,2)=-0.5D0
      DXI(3,1)=-0.5D0
      DXI(3,2)= 0.5D0
      DXI(4,1)=-0.5D0
      DXI(4,2)=-0.5D0
      DO 121 I=1,4
121   DOMEGA(I)=1D0
      NCUBP=4
      GOTO 99999
C
C
C
130   DXI (1,1)=-1D0
      DXI (1,2)=-1D0
      DXI (2,1)= 0D0
      DXI (2,2)=-1D0
      DXI (3,1)= 0D0
      DXI (3,2)= 0D0
      DXI (4,1)=-1D0
      DXI (4,2)= 0D0
C
      DXI (5,1)= 1D0
      DXI (5,2)=-1D0
      DXI (6,1)= 1D0
      DXI (6,2)= 0D0
      DXI (7,1)= 0D0
      DXI (7,2)= 0D0
      DXI (8,1)= 0D0
      DXI (8,2)=-1D0
C
      DXI (9,1)= 1D0
      DXI (9,2)= 1D0
      DXI(10,1)= 0D0
      DXI(10,2)= 1D0
      DXI(11,1)= 0D0
      DXI(11,2)= 0D0
      DXI(12,1)= 1D0
      DXI(12,2)= 0D0
C
      DXI(13,1)=-1D0
      DXI(13,2)= 1D0
      DXI(14,1)=-1D0
      DXI(14,2)= 0D0
      DXI(15,1)= 0D0
      DXI(15,2)= 0D0
      DXI(16,1)= 0D0
      DXI(16,2)= 1D0
C
      DO 131 I=1,16
131   DOMEGA(I)=0.25D0
      NCUBP=16
      GOTO 99999
C
C
C
140   DXI(1,1)=  0.7886751345948130D0
      DXI(1,2)=  0.7886751345948130D0
      DXI(2,1)=  0.2113248654051870D0
      DXI(2,2)=  0.7886751345948130D0
      DXI(3,1)=  0.2113248654051870D0
      DXI(3,2)=  0.2113248654051870D0
      DXI(4,1)=  0.7886751345948130D0
      DXI(4,2)=  0.2113248654051870D0
C
      DXI(5,1)= -0.7886751345948130D0
      DXI(5,2)=  0.7886751345948130D0
      DXI(6,1)= -0.2113248654051870D0
      DXI(6,2)=  0.7886751345948130D0
      DXI(7,1)= -0.2113248654051870D0
      DXI(7,2)=  0.2113248654051870D0
      DXI(8,1)= -0.7886751345948130D0
      DXI(8,2)=  0.2113248654051870D0
C
      DXI(9,1)=  0.7886751345948130D0
      DXI(9,2)= -0.7886751345948130D0
      DXI(10,1)= 0.2113248654051870D0
      DXI(10,2)=-0.7886751345948130D0
      DXI(11,1)= 0.2113248654051870D0
      DXI(11,2)=-0.2113248654051870D0
      DXI(12,1)= 0.7886751345948130D0
      DXI(12,2)=-0.2113248654051870D0
C
      DXI(13,1)=-0.7886751345948130D0
      DXI(13,2)=-0.7886751345948130D0
      DXI(14,1)=-0.2113248654051870D0
      DXI(14,2)=-0.7886751345948130D0
      DXI(15,1)=-0.2113248654051870D0
      DXI(15,2)=-0.2113248654051870D0
      DXI(16,1)=-0.7886751345948130D0
      DXI(16,2)=-0.2113248654051870D0
C
      DO 141 I=1,16
141   DOMEGA(I)=.25D0
      NCUBP=16
      GOTO 99999
C
C
C
150   DXI(1,1)=  0.8872983346207415D0
      DXI(1,2)=  0.8872983346207415D0
      DXI(2,1)=  0.1127016653792585D0
      DXI(2,2)=  0.8872983346207415D0
      DXI(3,1)=  0.8872983346207415D0
      DXI(3,2)=  0.1127016653792585D0
      DXI(4,1)=  0.1127016653792585D0
      DXI(4,2)=  0.1127016653792585D0
      DXI(5,1)=  0.8872983346207415D0
      DXI(5,2)=  0.5D0
      DXI(6,1)=  0.1127016653792585D0
      DXI(6,2)=  0.5D0
      DXI(7,1)=  0.5D0
      DXI(7,2)=  0.8872983346207415D0
      DXI(8,1)=  0.5D0
      DXI(8,2)=  0.1127016653792585D0
      DXI(9,1)=  0.5D0
      DXI(9,2)=  0.5D0
C
      DXI(10,1)=-0.8872983346207415D0
      DXI(10,2)= 0.8872983346207415D0
      DXI(11,1)=-0.1127016653792585D0
      DXI(11,2)= 0.8872983346207415D0
      DXI(12,1)=-0.8872983346207415D0
      DXI(12,2)= 0.1127016653792585D0
      DXI(13,1)=-0.1127016653792585D0
      DXI(13,2)= 0.1127016653792585D0
      DXI(14,1)=-0.8872983346207415D0
      DXI(14,2)= 0.5D0
      DXI(15,1)=-0.1127016653792585D0
      DXI(15,2)= 0.5D0
      DXI(16,1)=-0.5D0
      DXI(16,2)= 0.8872983346207415D0
      DXI(17,1)=-0.5D0
      DXI(17,2)= 0.1127016653792585D0
      DXI(18,1)=-0.5D0
      DXI(18,2)= 0.5D0
C
      DXI(19,1)= 0.8872983346207415D0
      DXI(19,2)=-0.8872983346207415D0
      DXI(20,1)= 0.1127016653792585D0
      DXI(20,2)=-0.8872983346207415D0
      DXI(21,1)= 0.8872983346207415D0
      DXI(21,2)=-0.1127016653792585D0
      DXI(22,1)= 0.1127016653792585D0
      DXI(22,2)=-0.1127016653792585D0
      DXI(23,1)= 0.8872983346207415D0
      DXI(23,2)=-0.5D0
      DXI(24,1)= 0.1127016653792585D0
      DXI(24,2)=-0.5D0
      DXI(25,1)= 0.5D0
      DXI(25,2)=-0.8872983346207415D0
      DXI(26,1)= 0.5D0
      DXI(26,2)=-0.1127016653792585D0
      DXI(27,1)= 0.5D0
      DXI(27,2)=-0.5D0
C
      DXI(28,1)=-0.8872983346207415D0
      DXI(28,2)=-0.8872983346207415D0
      DXI(29,1)=-0.1127016653792585D0
      DXI(29,2)=-0.8872983346207415D0
      DXI(30,1)=-0.8872983346207415D0
      DXI(30,2)=-0.1127016653792585D0
      DXI(31,1)=-0.1127016653792585D0
      DXI(31,2)=-0.1127016653792585D0
      DXI(32,1)=-0.8872983346207415D0
      DXI(32,2)=-0.5D0
      DXI(33,1)=-0.1127016653792585D0
      DXI(33,2)=-0.5D0
      DXI(34,1)=-0.5D0
      DXI(34,2)=-0.8872983346207415D0
      DXI(35,1)=-0.5D0
      DXI(35,2)=-0.1127016653792585D0
      DXI(36,1)=-0.5D0
      DXI(36,2)=-0.5D0
C
      DO 151 J=0,27,9
      DO 151 I=1,4
      DOMEGA(I+J)  =0.7716049382716050D-1
151   DOMEGA(I+J+4)=0.1234567901234570D0
      DO 152 I=9,36,9
152   DOMEGA(I)=0.1975308641975310D0
C
      NCUBP=36
C
99999 END
