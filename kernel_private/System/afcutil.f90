
! Utility routines for the FEM-FCT solver

module afcutil

  use fsystem

  implicit none

contains

  ! Associate a pointer with a Feat2D array
  
  function feat_htpdouble(n,handle) result(p)
    integer, intent(IN) :: n,handle
    double precision, dimension(:), pointer :: p
    include 'cmem.inc'
    
    p => feat_vr(n,DWORK(L(handle)))
    
  end function

  function feat_htpdouble2D(n,m,handle) result(p)
    integer, intent(IN) :: n,m,handle
    double precision, dimension(:,:), pointer :: p
    include 'cmem.inc'
    
    p => feat_mr(n,m,DWORK(L(handle)))
    
  end function

  function feat_htpint(n,handle) result(p)
    integer, intent(IN) :: n,handle
    integer, dimension(:), pointer :: p
    include 'cmem.inc'
    
    p => feat_vi(n,KWORK(L(handle)))
    
  end function

  function feat_htpint2D(n,m,handle) result(p)
    integer, intent(IN) :: n,m,handle
    integer, dimension(:,:), pointer :: p
    include 'cmem.inc'
    
    p => feat_mi(n,m,KWORK(L(handle)))
    
  end function
  

  function feat_vr(n,arr) result(p)
    integer, intent(IN) :: n
    double precision, dimension(n), target :: arr
    double precision, dimension(:), pointer :: p

    p => arr

  end function feat_vr

  function feat_mr(n,m,arr) result(p)
    integer, intent(IN) :: n,m
    double precision, dimension(n,m), target :: arr
    double precision, dimension(:,:), pointer :: p

    p => arr

  end function feat_mr

  function feat_vi(n,arr) result(p)
    integer, intent(IN) :: n
    integer, dimension(n), target :: arr
    integer, dimension(:), pointer :: p

    p => arr

  end function feat_vi

  function feat_mi(n,m,arr) result(p)
    integer, intent(IN) :: n,m
    integer, dimension(n,m), target :: arr
    integer, dimension(:,:), pointer :: p

    p => arr

  end function feat_mi

end module 
