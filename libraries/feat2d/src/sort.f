************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek, M.Koester          *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* VECSRT                                                               *
*                                                                      *
* Purpose: Resorting of a vector                                       *
*                                                                      *
* Resorts the entries in the vector DX corresponding to KTR1/KTR2.     *
* IPAR2 describes the type of resorting. The result is written to DD.  *
*                                                                      *
* Version from  10/15/04                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*4(NEQ)   Source vector to be resorted.                    *
* NEQ      I*4        Number of equations/entries in DX/KTRx           *
* KTR1     I*4(NEQ)   array with permitation 1                         *
* KTR2     I*4(NEQ)   array with permitation 2                         *
* IPAR     I*4        Type of resorting:                               *
*                     =1: resort vector corresponding to permutation   *
*                         KTR1 (1..NEQ)                                *
*                     =2: resort vector corresponding to permutation   *
*                         KTR2 (1..NEQ)                                *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DD       R*4(NEQ)   The resorted vector.                             *
*                                                                      *
************************************************************************

      SUBROUTINE VECSRT(DX,DD,KTR1,KTR2,NEQ,IPAR)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*),DD(*),KTR1(*),KTR2(*)
      
        IF (IPAR.EQ.1) THEN
          DO IEQ=1,NEQ
            DD(IEQ)=DX(KTR1(IEQ))
          END DO  
        ELSE
          DO IEQ=1,NEQ
            DD(IEQ)=DX(KTR2(IEQ))
          END DO
        ENDIF
        
      END

************************************************************************
*                                                                      *
* MTSRTD                                                               *
*                                                                      *
* Purpose: Resorting of a matrix                                       *
*                                                                      *
* Resorts the entries of the given matrix, corresponding to KTR1/KTR2. *
* Double precision version, storage technique 7 only.                  *
*                                                                      *
* Version from  10/15/04                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DAH      R*4(NEQ)   Source matrix to be resorted.                    *
* KCOLH    I*4(NEQ)   Column structure of source matrix                *
* KLDH     I*4(NEQ)   Row positions of source matrix                   *
* KTR1     I*4        Permutation of 1..NEQ describing how to resort.  *
* KTR2     I*4(NEQ)   Permutation of 1..NEQ describing how to sort     *
*                     the matrix back.                                 *
* IPAR     I*4        Type of resorting:                               *
*                     =1: resort vector corresponding to permutation   *
*                         KTR1 (1..NEQ)                                *
*                     =2: resort vector corresponding to permutation   *
*                         KTR2 (1..NEQ)                                *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DAH      R*4(NEQ)   Resorted matrix                                  *
* KCOLH    I*4(NEQ)   Column structure of resorted matrix              *
* KLDH     I*4(NEQ)   Row positions of resorted matrix                 *
*                                                                      *
************************************************************************

      SUBROUTINE MTSRTD(DA,DAH,KCOL,KCOLH,KLD,KLDH,KTR1,KTR2,NEQ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),DAH(*),KCOL(*),KCOLH(*),KLD(*),KLDH(*),
     *          KTR1(*),KTR2(*)
      DIMENSION KH1(100),KH2(100)

        ILD=1
        DO I=1,NEQ
          I1       =KTR1(I)
          DA(ILD)  =DAH(KLDH(I1))
          KLD(I)   =ILD
          KCOL(ILD)=I
          ILD      =ILD+1
          IH1=KLDH(I1)+1
          IH2=KLDH(I1+1)-1
          ID1=IH2-IH1+1
          DO JH=IH1,IH2
            KH1(JH-IH1+1)=KTR2(KCOLH(JH))
            KH2(JH-IH1+1)=JH
          END DO
          CALL KHSORT(KH1,KH2,ID1)
          DO J=1,ID1
            DH       =DAH(KH2(J))
            ICOL     =KH1(J)
            DA(ILD)  =DH
            KCOL(ILD)=ICOL
            ILD      =ILD+1
          END DO
        END DO
        KLD(NEQ+1)=KLDH(NEQ+1)

      
      END

************************************************************************
*                                                                      *
* MTSRTV                                                               *
*                                                                      *
* Purpose: Resorting of a matrix                                       *
*                                                                      *
* Resorts the entries of the given matrix, corresponding to KTR1/KTR2. *
* Single precision version, storage technique 7 only.                  *
*                                                                      *
* Version from  10/15/04                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DAH      R*4(NEQ)   Source matrix to be resorted.                    *
* KCOLH    I*4(NEQ)   Column structure of source matrix                *
* KLDH     I*4(NEQ)   Row positions of source matrix                   *
* KTR1     I*4        Permutation of 1..NEQ describing how to resort.  *
* KTR2     I*4(NEQ)   Permutation of 1..NEQ describing how to sort     *
*                     the matrix back.                                 *
* IPAR     I*4        Type of resorting:                               *
*                     =1: resort vector corresponding to permutation   *
*                         KTR1 (1..NEQ)                                *
*                     =2: resort vector corresponding to permutation   *
*                         KTR2 (1..NEQ)                                *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DAH      R*4(NEQ)   Resorted matrix                                  *
* KCOLH    I*4(NEQ)   Column structure of resorted matrix              *
* KLDH     I*4(NEQ)   Row positions of resorted matrix                 *
*                                                                      *
************************************************************************

      SUBROUTINE MTSRTV(VA,VAH,KCOL,KCOLH,KLD,KLDH,KTR1,KTR2,NEQ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),VAH(*),KCOL(*),KCOLH(*),KLD(*),KLDH(*),
     *          KTR1(*),KTR2(*)
      DIMENSION KH1(100),KH2(100)

        ILD=1
        DO I=1,NEQ
          I1       =KTR1(I)
          VA(ILD)  =VAH(KLDH(I1))
          KLD(I)   =ILD
          KCOL(ILD)=I
          ILD      =ILD+1
          IH1=KLDH(I1)+1
          IH2=KLDH(I1+1)-1
          ID1=IH2-IH1+1
          DO JH=IH1,IH2
            KH1(JH-IH1+1)=KTR2(KCOLH(JH))
            KH2(JH-IH1+1)=JH
          END DO
          CALL KHSORT(KH1,KH2,ID1)
          DO J=1,ID1
            VH       =VAH(KH2(J))
            ICOL     =KH1(J)
            VA(ILD)  =VH
            KCOL(ILD)=ICOL
            ILD      =ILD+1
          END DO
        END DO
        KLD(NEQ+1)=KLDH(NEQ+1)

      END

************************************************************************
*                                                                      *
* MTSRTR                                                               *
*                                                                      *
* Purpose: Resorting back the matrix structure                         *
*                                                                      *
* Resorts the entries of the given matrix, corresponding to KTR1/KTR2: *
* Sorts the structure of the given matrix back.                        *
* Storage technique 7. Only the structure of a given matrix is         *
* resorted, the routine doesn't handle the entries.                    *
*                                                                      *
* Version from  10/15/04                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KCOLH    I*4(NEQ)   Column structure of source matrix                *
* KLDH     I*4(NEQ)   Row positions of source matrix                   *
* KTR1     I*4        Permutation of 1..NEQ describing how to resort.  *
* KTR2     I*4(NEQ)   Permutation of 1..NEQ describing how to sort     *
*                     the matrix back.                                 *
* IPAR     I*4        Type of resorting:                               *
*                     =1: resort vector corresponding to permutation   *
*                         KTR1 (1..NEQ)                                *
*                     =2: resort vector corresponding to permutation   *
*                         KTR2 (1..NEQ)                                *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KCOLH    I*4(NEQ)   Column structure of resorted matrix              *
* KLDH     I*4(NEQ)   Row positions of resorted matrix                 *
*                                                                      *
************************************************************************

      SUBROUTINE MTSRTR(KCOL,KCOLH,KLD,KLDH,KTR1,KTR2,NEQ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KCOL(*),KCOLH(*),KLD(*),KLDH(*),
     *          KTR1(*),KTR2(*)
      DIMENSION KH1(100),KH2(100)

        ILD=1
        DO I=1,NEQ
          I1       =KTR2(I)
          KLD(I)   =ILD
          KCOL(ILD)=I
          ILD      =ILD+1
          IH1=KLDH(I1)+1
          IH2=KLDH(I1+1)-1
          ID1=IH2-IH1+1
          DO JH=IH1,IH2
            KH1(JH-IH1+1)=KTR1(KCOLH(JH))
            KH2(JH-IH1+1)=JH
          END DO
          CALL KHSORT(KH1,KH2,ID1)
          DO J=1,ID1
            ICOL     =KH1(J)
            KCOL(ILD)=ICOL
            ILD      =ILD+1
          END DO
        END DO
        KLD(NEQ+1)=KLDH(NEQ+1)

      END
C
C *******************************************************************
C *******************************************************************

      SUBROUTINE KHSORT(KH1,KH2,IDIM)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KH1(*),KH2(*)

        BMORE=.TRUE.
        
C WHILE loop, realized with IF...
        
5       IF (.NOT.BMORE) GOTO 99999
        BMORE=.FALSE.
        DO ICOMP=1,IDIM-1
          IF (KH1(ICOMP).GT.KH1(ICOMP+1)) THEN
            IAUX1=KH1(ICOMP)
            IAUX2=KH2(ICOMP)
            KH1(ICOMP)=KH1(ICOMP+1)
            KH2(ICOMP)=KH2(ICOMP+1)
            KH1(ICOMP+1)=IAUX1
            KH2(ICOMP+1)=IAUX2
            BMORE=.TRUE.
          ENDIF
        END DO

        GOTO 5

99999 END
