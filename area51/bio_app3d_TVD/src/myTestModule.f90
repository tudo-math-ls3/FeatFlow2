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


  !<subroutine>
  ! test purposes: print vektor
  subroutine vector_writeToFile (rtria, rvector, fileName, fileNumber, itimestep)
  
    type(t_vectorScalar), intent(IN) :: rvector
    character (LEN=*), intent(IN) :: fileName
    integer,  intent(INOUT), optional :: fileNumber
    integer,  intent(IN), optional :: itimestep
    type(t_triangulation), intent(IN) :: rtria
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    !auxiliary
    real(DP), dimension(:), pointer :: p_vector
    integer :: i, j, nvt
    integer :: fileNumber_local
  
    ! set variables
    call lsyssc_getbase_double (rvector,p_vector)
    nvt = rvector%NEQ
    if(present(fileNumber)) then
	fileNumber_local=fileNumber
    else
	fileNumber_local=1313
    endif

    ! get coordinates
    call storage_getbase_double2d (rtria%h_DvertexCoords,p_DvertexCoords)

    !open(fileNumber,FILE=fileName//trim(sys_si0L(itimestep,5)),STATUS="NEW")
    if(present(itimestep)) then
	open(fileNumber_local,FILE=fileName//trim('_y05')//trim(sys_si0L(itimestep,5)),STATUS="REPLACE")
    else
    	open(fileNumber_local,FILE=fileName//trim('_y05'),STATUS="REPLACE")
    endif
    DO i = 1,nvt
	!print fileNumber,p_vector(i)
	!write(mfile,'(A,1X,I10)') 'nodes',rexport%nvertices
	if( p_DvertexCoords(2,i).eq.0.5D0 ) then
	    write(fileNumber_local,'(ES16.8E3,1X,ES16.8E3,1X,ES16.8E3)') &
				     p_DvertexCoords(1,i), p_DvertexCoords(2,i), p_vector(i)
        endif
    end do
    close(fileNumber_local)

    !write x=0.5 cutline
    !open(fileNumber,FILE=fileName//trim(sys_si0L(itimestep,5)),STATUS="NEW")
    if(present(itimestep)) then
	open(fileNumber_local,FILE=fileName//trim('_x05')//trim(sys_si0L(itimestep,5)),STATUS="REPLACE")
    else
    	open(fileNumber_local,FILE=fileName//trim('_x05'),STATUS="REPLACE")
    endif
    DO i = 1,nvt
	!print fileNumber,p_vector(i)
	!write(mfile,'(A,1X,I10)') 'nodes',rexport%nvertices
	if( p_DvertexCoords(1,i).eq.0.5D0 ) then
	    write(fileNumber_local,'(ES16.8E3,1X,ES16.8E3,1X,ES16.8E3)') &
				     p_DvertexCoords(1,i), p_DvertexCoords(2,i), p_vector(i)
        endif
    end do
    close(fileNumber_local)
	
  end subroutine
  !</subroutine>


end module
