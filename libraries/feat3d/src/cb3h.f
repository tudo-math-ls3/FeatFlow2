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
* CB3H                                                                 *
*                                                                      *
* Purpose  Set coordinates and weights for numerical integration       *
*          over the unit cube [-1,1]x[-1,1]x[-1,1]                     *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  09/09/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ICUB     I*4    Number of cubature formula                           *
*                 ICUB    NCUB    DEGREE                               *
*                   1       1        1      1x1x1 Gaussian formula     *
*                   2       6        1      Midpoints of areas         *
*                   3       8        1      Product trapezoidal rule   *
*                   4      12        1      Midpoints of edges (nyi)   *
*                   5       4        2      Stroud                     *
*                   6       6        3      Stroud                     *
*                   7       8        3      2x2x2 Gaussian formula     *
*                   8      14        5      Hammer and Stroud          *
*                   9      27        5      3x3x3 Gaussian formula     *
*                  10      34        7      Sarma and Stroud           *
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
      SUBROUTINE CB3H (ICUB)
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
      SUB='CB3H'
      IF (ICHECK.GE.997) CALL OTRC('CB3H  ','09/09/89')
C
      IER=0
C
      GOTO (10,20,30,999,50,60,70,80,90,100) ICUB
C
999   WRITE (CPARAM,'(I15)') ICUB
      CALL WERR(-119,'CB3H  ')
      GOTO 99999
C
C
C
10    DXI(1,1)= 0D0
      DXI(1,2)= 0D0
      DXI(1,3)= 0D0
      DOMEGA(1)=8D0
      NCUBP=1
      GOTO 99999
C
C
C
20    DO 21 I=1,6
      DO 21 J=1,3
21    DXI(I,J)= 0D0
      DXI(1,2)= 0D0
      DXI(1,3)=-1D0
      DXI(2,1)=-1D0
      DXI(3,2)=-1D0
      DXI(4,1)= 1D0
      DXI(5,2)= 1D0
      DXI(6,3)= 1D0
      DO 22 I=1,6
22    DOMEGA(I)=1.333333333333333D0
      NCUBP=6
      GOTO 99999
C
C
C
30    DO 31 I=1,5,4
      DXI(I  ,1)=-1D0
      DXI(I  ,2)=-1D0
      DXI(I+1,1)= 1D0
      DXI(I+1,2)=-1D0
      DXI(I+2,1)= 1D0
      DXI(I+2,2)= 1D0
      DXI(I+3,1)=-1D0
31    DXI(I+3,2)= 1D0
      DO 32 I=1,4
32    DXI(I,3)=-1D0
      DO 33 I=5,8
33    DXI(I,3)= 1D0
      DO 34 I=1,8
34    DOMEGA(I)=1D0
      NCUBP=8
      GOTO 99999
C
C
C
50    DXI(1,1)= 0.816496580927726D0
      DXI(1,2)= 0D0
      DXI(1,3)= 0.577350269189626D0
      DXI(2,1)= 0D0
      DXI(2,2)= 0.816496580927726D0
      DXI(2,3)=-0.577350269189626D0
      DXI(3,1)=-0.816496580927726D0
      DXI(3,2)= 0D0
      DXI(3,3)= 0.577350269189626D0
      DXI(4,1)= 0D0
      DXI(4,2)=-0.816496580927726D0
      DXI(4,3)=-0.577350269189626D0
      DO 51 I=1,4
51    DOMEGA(I)=2D0
      NCUBP=4
      GOTO 99999
C
C
C
60    DXI(1,1)= 0.408248290463863D0
      DXI(1,2)= 0.707106781186547D0
      DXI(1,3)=-0.577350269189626D0
      DXI(2,1)=-0.408248290463863D0
      DXI(2,2)= 0.707106781186547D0
      DXI(2,3)= 0.577350269189626D0
      DXI(3,1)=-0.816496580927726D0
      DXI(3,2)= 0D0
      DXI(3,3)=-0.577350269189626D0
      DXI(4,1)=-0.408248290463863D0
      DXI(4,2)=-0.707106781186547D0
      DXI(4,3)= 0.577350269189626D0
      DXI(5,1)= 0.408248290463863D0
      DXI(5,2)=-0.707106781186547D0
      DXI(5,3)=-0.577350269189626D0
      DXI(6,1)= 0.816496580927726D0
      DXI(6,2)= 0D0
      DXI(6,3)= 0.577350269189626D0
      DO 61 I=1,6
61    DOMEGA(I)=1.333333333333333D0
      NCUBP=6
      GOTO 99999
C
C
C
70    DO 71 I=1,5,4
      DXI(I  ,1)= 0.577350269189626D0
      DXI(I  ,2)= 0.577350269189626D0
      DXI(I+1,1)=-0.577350269189626D0
      DXI(I+1,2)= 0.577350269189626D0
      DXI(I+2,1)=-0.577350269189626D0
      DXI(I+2,2)=-0.577350269189626D0
      DXI(I+3,1)= 0.577350269189626D0
71    DXI(I+3,2)=-0.577350269189626D0
      DO 72 I=1,4
72    DXI(I,3)=-0.577350269189626D0
      DO 73 I=5,8
73    DXI(I,3)= 0.577350269189626D0
      DO 74 I=1,8
74    DOMEGA(I)=1D0
      NCUBP=8
      GOTO 99999
C
C
C
80    DO 81 I=1,6
      DO 81 J=1,3
81    DXI(I,J)= 0D0
      DXI(1,1)= 0.795822425754221D0
      DXI(2,2)= 0.795822425742215D0
      DXI(3,3)= 0.795822425754221D0
      DXI(4,1)=-0.795822425754221D0
      DXI(5,2)=-0.795822425754221D0
      DXI(6,3)=-0.795822425754221D0
      DO 82 I=7,11,4
      DXI(I  ,1)= 0.758786910639328D0
      DXI(I  ,2)= 0.758786910639328D0
      DXI(I+1,1)=-0.758786910639328D0
      DXI(I+1,2)= 0.758786910639328D0
      DXI(I+2,1)=-0.758786910639328D0
      DXI(I+2,2)=-0.758786910639328D0
      DXI(I+3,1)= 0.758786910639328D0
82    DXI(I+3,2)=-0.758786910639328D0
      DO 83 I=7,10
83    DXI(I,3)=-0.758786910639328D0
      DO 84 I=11,14
84    DXI(I,3)= 0.758786910639328D0
      DO 85 I=1,6
85    DOMEGA(I)=0.886426592797784D0
      DO 86 I=7,14
86    DOMEGA(I)=0.335180055401662D0
      NCUBP=14
      GOTO 99999
C
C
C
90    DO 91 I=1,19,9
      DXI(I  ,1)= 0.774596669241483D0
      DXI(I  ,2)= 0.774596669241483D0
      DXI(I+1,1)=-0.774596669241483D0
      DXI(I+1,2)= 0.774596669241483D0
      DXI(I+2,1)= 0.774596669241483D0
      DXI(I+2,2)=-0.774596669241483D0
      DXI(I+3,1)=-0.774596669241483D0
      DXI(I+3,2)=-0.774596669241483D0
      DXI(I+4,1)= 0.774596669241483D0
      DXI(I+4,2)= 0D0
      DXI(I+5,1)=-0.774596669241483D0
      DXI(I+5,2)= 0D0
      DXI(I+6,1)= 0D0
      DXI(I+6,2)= 0.774596669241483D0
      DXI(I+7,1)= 0D0
      DXI(I+7,2)=-0.774596669241483D0
      DXI(I+8,1)= 0D0
91    DXI(I+8,2)= 0D0
      DO 920 I=1,9
920   DXI(I,3)=-0.774596669241483D0
      DO 921 I=10,18
921   DXI(I,3)= 0D0
      DO 922 I=19,27
922   DXI(I,3)= 0.774596669241483D0
      DO 930 I=1,4
930   DOMEGA(I)= 0.17146776406036D0
      DO 931 I=5,8
931   DOMEGA(I)= 0.27434842249657D0
      DOMEGA(9)= 0.43895747599451D0
      DO 940 I=10,13
940   DOMEGA(I)= 0.27434842249657D0
      DO 941 I=14,17
941   DOMEGA(I)= 0.43895747599451D0
      DOMEGA(18)= 0.70233196159122D0
      DO 950 I=19,22
950   DOMEGA(I)= 0.17146776406036D0
      DO 951 I=23,26
951   DOMEGA(I)= 0.27434842249657D0
      DOMEGA(27)= 0.43895747599451D0
      NCUBP=27
      GOTO 99999
C
C
C
100   DO 101 I=1,6
      DO 101 J=1,3
101   DXI(I,J)= 0D0
      DXI(1,1)= 0.9317380000D0
      DXI(2,2)= 0.9317380000D0
      DXI(3,3)= 0.9317380000D0
      DXI(4,1)=-0.9317380000D0
      DXI(5,2)=-0.9317380000D0
      DXI(6,3)=-0.9317380000D0
      DO 102 I=7,12
      DO 102 J=1,3
102   DXI(I,J)= 0.9167441779D0
      DXI( 7,3)= 0D0
      DXI( 8,2)= 0D0
      DXI( 9,1)= 0D0
      DXI(10,2)=-0.9167441779D0
      DXI(10,3)= 0D0
      DXI(11,2)= 0D0
      DXI(11,3)=-0.9167441779D0
      DXI(12,1)= 0D0
      DXI(12,3)=-0.9167441779D0
      DO 103 I=13,18
      DO 103 J=1,3
103   DXI(I,J)=-0.9167441779D0
      DXI(13,3)= 0D0
      DXI(14,2)= 0D0
      DXI(15,1)= 0D0
      DXI(16,2)= 0.9167441779D0
      DXI(16,3)= 0D0
      DXI(17,2)= 0D0
      DXI(17,3)= 0.9167441779D0
      DXI(18,1)= 0D0
      DXI(18,3)= 0.9167441779D0
      DO 104 I=19,23,4
      DXI(I  ,1)= 0.4086003800D0
      DXI(I  ,2)= 0.4086003800D0
      DXI(I+1,1)=-0.4086003800D0
      DXI(I+1,2)= 0.4086003800D0
      DXI(I+2,1)=-0.4086003800D0
      DXI(I+2,2)=-0.4086003800D0
      DXI(I+3,1)= 0.4086003800D0
104   DXI(I+3,2)=-0.4086003800D0
      DO 105 I=19,22
105   DXI(I,3)=-0.4086003800D0
      DO 106 I=23,26
106   DXI(I,3)= 0.4086003800D0
      DO 107 I=27,31,4
      DXI(I  ,1)= 0.7398529500D0
      DXI(I  ,2)= 0.7398529500D0
      DXI(I+1,1)=-0.7398529500D0
      DXI(I+1,2)= 0.7398529500D0
      DXI(I+2,1)=-0.7398529500D0
      DXI(I+2,2)=-0.7398529500D0
      DXI(I+3,1)= 0.7398529500D0
107   DXI(I+3,2)=-0.7398529500D0
      DO 108 I=27,30
108   DXI(I,3)=-0.7398529500D0
      DO 109 I=31,34
109   DXI(I,3)= 0.7398529500D0
      DO 110 I=1,6
110   DOMEGA(I)=0.284654471680D0
      DO 111 I=7,18
111   DOMEGA(I)=0.998314216000D-1
      DO 112 I=19,26
112   DOMEGA(I)=0.422941839280D0
      DO 113 I=27,34
113   DOMEGA(I)=0.213820174560D0
      NCUBP=34
C
99999 END
