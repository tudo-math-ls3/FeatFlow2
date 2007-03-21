
! Utility routines for the FEM-FCT solver

MODULE afcutil

  USE fsystem

  IMPLICIT NONE

CONTAINS

  ! Associate a pointer with a Feat2D array
  
  FUNCTION feat_htpdouble(n,handle) RESULT(p)
    INTEGER, INTENT(IN) :: n,handle
    DOUBLE PRECISION, DIMENSION(:), POINTER :: p
    INCLUDE 'cmem.inc'
    
    p => feat_vr(n,DWORK(L(handle)))
    
  END FUNCTION

  FUNCTION feat_htpdouble2D(n,m,handle) RESULT(p)
    INTEGER, INTENT(IN) :: n,m,handle
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: p
    INCLUDE 'cmem.inc'
    
    p => feat_mr(n,m,DWORK(L(handle)))
    
  END FUNCTION

  FUNCTION feat_htpint(n,handle) RESULT(p)
    INTEGER, INTENT(IN) :: n,handle
    INTEGER, DIMENSION(:), POINTER :: p
    INCLUDE 'cmem.inc'
    
    p => feat_vi(n,KWORK(L(handle)))
    
  END FUNCTION

  FUNCTION feat_htpint2D(n,m,handle) RESULT(p)
    INTEGER, INTENT(IN) :: n,m,handle
    INTEGER, DIMENSION(:,:), POINTER :: p
    INCLUDE 'cmem.inc'
    
    p => feat_mi(n,m,KWORK(L(handle)))
    
  END FUNCTION
  

  FUNCTION feat_vr(n,arr) RESULT(p)
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, DIMENSION(n), TARGET :: arr
    DOUBLE PRECISION, DIMENSION(:), POINTER :: p

    p => arr

  END FUNCTION feat_vr

  FUNCTION feat_mr(n,m,arr) RESULT(p)
    INTEGER, INTENT(IN) :: n,m
    DOUBLE PRECISION, DIMENSION(n,m), TARGET :: arr
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: p

    p => arr

  END FUNCTION feat_mr

  FUNCTION feat_vi(n,arr) RESULT(p)
    INTEGER, INTENT(IN) :: n
    INTEGER, DIMENSION(n), TARGET :: arr
    INTEGER, DIMENSION(:), POINTER :: p

    p => arr

  END FUNCTION feat_vi

  FUNCTION feat_mi(n,m,arr) RESULT(p)
    INTEGER, INTENT(IN) :: n,m
    INTEGER, DIMENSION(n,m), TARGET :: arr
    INTEGER, DIMENSION(:,:), POINTER :: p

    p => arr

  END FUNCTION feat_mi

END MODULE 
