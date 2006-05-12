************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek, M. Koester         *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* MSRNA                                                                *
*                                                                      *
* Purpose: Calculate number of matrix entries in MSR matrix            *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* JLU      I*8(*)     Index set JLU of MSR structure                   *
* NEQ      I*4        Number of equations                              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* return   I*4        Number of elements NA in MSR matrix              *
*                                                                      *
************************************************************************

      INTEGER FUNCTION MSRNA (JLU,NEQ)
      IMPLICIT NONE
      INTEGER JLU(*),NEQ

      MSRNA = JLU(NEQ+1)-JLU(1)+NEQ

      END      

************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek, M. Koester         *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* MSR2T7                                                               *
*                                                                      *
* Purpose: Convert matrix in SPLIB's MSR structure to matrix           *
*          format 7 of FEAT                                            *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LU       I*8(*)     Pointer to MSR matrix                            *
* JLU      I*8(*)     Index set JLU of MSR structure                   *
* NEQ      I*4        Number of equations                              *
* ICPMT    I*4        =1: copy content of matrix to DA                 *
*                     =0: only build the structure of the matrix.      *
*                         DA and LU is ignored.                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DA       R*8(NA)    The resulting matrix in memory technique 7       *
* KCOL     I*4(NA)    Column description of matrix                     *
* KLD      I*4(NEQ+1) Row description of matrix                        *
*                                                                      *
* This subroutine assumes that the caller provides the correct amount  *
* of memory in DA/KCOL/KLD; it only fills this structure. Therefore    *
* for converting a matrix, the caller has to provide NA elements       *
* in either KCOL and DA and NEQ+1 elements in KLD. To calculate NA,    *
* the above subroutine MSRNA can be used.                              *
*                                                                      *
************************************************************************

      SUBROUTINE MSR2T7 (LU, JLU, DA, KCOL, KLD, NEQ, ICPMT)
      
      IMPLICIT NONE
      DOUBLE PRECISION LU(*), DA(*)
      INTEGER I, JLU(*), KCOL(*), KLD(*), NEQ
      INTEGER NLCNT,ICPMT

        KLD(1) = 1
      
        DO I=1,NEQ

C Number of elements in source row    
          NLCNT = JLU(I+1)-JLU(I)
  
C Calculate KLD entry; add 1 for the diagonal          

          KLD(I+1) = KLD(I)+NLCNT+1

C Insert number of diagonal and copy the column numbers of 
C the current row from MSR matrix structure to KCOL

          KCOL(KLD(I)) = I
          CALL LCP3(JLU(JLU(I)),KCOL(KLD(I)+1),NLCNT)
      
C If necessary, copy the elements below the diagonal to DA

          IF (ICPMT.NE.0) THEN
          
C   Invert the diagonal entry and copy it - because in SPLIB format the
C   diagonal entry is stored inverted.

            DA(KLD(I)) = 1/LU(I)
            CALL LCP1 (LU(JLU(I)), DA(KLD(I)+1),NLCNT)
          END IF
        END DO
      
      END

