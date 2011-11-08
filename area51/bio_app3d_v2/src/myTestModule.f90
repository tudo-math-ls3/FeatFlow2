module myTestModule

  use fsystem
  use storage
  use genoutput
  use linearsolver
  use boundary
  use triangulation
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use derivatives
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use bcassembly
  use element
  implicit none
	
  ! Defining pi
  real(DP), parameter, public :: PI = 3.141592654_DP

contains

  !<subroutine>
  ! test purposes: print matrix (written in the format 9)
  subroutine testPrintMatrix9 (rmatrix)
  
    type(t_matrixScalar) , intent(IN) :: rmatrix
    
    ! some local variables
    integer :: i, j, nvt
    integer(PREC_VECIDX), dimension(:), pointer :: p_Kcol
    integer(PREC_MATIDX), dimension(:), pointer :: p_Kld , p_Kdiagonal
    real(DP), dimension(:), pointer :: p_Da_mass, p_Da
    integer , dimension(:) , allocatable :: iaux

    call lsyssc_getbase_Kld (rmatrix,p_Kld)
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
    call lsyssc_getbase_double (rmatrix,p_Da)
    nvt = rmatrix%NEQ
    
    DO i = 1,nvt
        Do j = p_Kld(i), p_Kld(i+1)-1
            print *,i,": ",p_Kcol(j),": ",p_Da( j )
        end Do
    end do
    
  end subroutine
  !</subroutine>


  !<subroutine>
  ! test purposes: print vektor
  subroutine testPrintVector (rvector)
  
    type(t_vectorScalar) , intent(IN) :: rvector
    real(DP), dimension(:), pointer :: p_vector
        
    ! some local variables
    integer :: i, j, nvt

    call lsyssc_getbase_double (rvector,p_vector)
    nvt = rvector%NEQ
    
    DO i = 1,nvt
        print *,i,": ",p_vector( i )
    end do
    
  end subroutine
  !</subroutine>


  !<subroutine>
  ! multiplication of the inverse lumped matrix by vector
  subroutine multInvLumpedVector (rmatrix, rvector)

    type(t_matrixScalar) , intent(IN) :: rmatrix
    type(t_vectorScalar) , intent(INOUT) :: rvector
    real(DP), dimension(:), pointer :: p_vector
    
    ! some local variables
    integer :: i, nvt
    integer(PREC_VECIDX), dimension(:), pointer :: p_Kcol
    integer(PREC_MATIDX), dimension(:), pointer :: p_Kld , p_Kdiagonal
    real(DP), dimension(:), pointer :: p_Da_mass, p_Da

    !matrix initialization things
    !call lsyssc_getbase_Kld (rmatrix,p_Kld)
    !call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
    call lsyssc_getbase_double (rmatrix,p_Da)
    nvt = rmatrix%NEQ
    
    !vector initialization things
    call lsyssc_getbase_double (rvector,p_vector)
    nvt = rvector%NEQ
    
    DO i = 1,nvt
        p_vector( i ) = (1.0_DP/p_Da(p_Kdiagonal(i)))*p_vector( i )
    end do
    
  end subroutine

end module
