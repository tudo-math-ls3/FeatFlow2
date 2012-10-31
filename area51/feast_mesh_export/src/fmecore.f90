module fmecore

  use fsystem
  use genoutput
  use paramlist
  use storage
  use triangulation

  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fme_writeCoordsChunk(iunit, rtria, ndigits)
  integer, intent(IN) :: iunit
  type(t_triangulation), intent(INOUT) :: rtria
  integer, intent(in) :: ndigits

!</subroutine>

  integer :: i, j
  real(DP), dimension(:,:), pointer :: p_Dvtx
  
    ! fetch vertex coords array
    call storage_getbase_double2d(rtria%h_DvertexCoords, p_Dvtx)
    
    ! write chunk starter
    write(iunit, '(A)') ' <coords>'
    
    ! write coordinates
    select case(rtria%ndim)
    case (1)
      do i = 1, rtria%NVT
        write(iunit, '(A)') '  ' // &
          trim(sys_sdel(p_Dvtx(1,i), ndigits))
      end do
      
    case (2)
      do i = 1, rtria%NVT
        write(iunit, '(A)') '  ' // &
          trim(sys_sdel(p_Dvtx(1,i), ndigits)) // ' ' // &
          trim(sys_sdel(p_Dvtx(2,i), ndigits))
      end do

    case (3)
      do i = 1, rtria%NVT
        write(iunit, '(A)') '  ' // &
          trim(sys_sdel(p_Dvtx(1,i), ndigits)) // ' ' // &
          trim(sys_sdel(p_Dvtx(2,i), ndigits)) // ' ' // &
          trim(sys_sdel(p_Dvtx(3,i), ndigits))
      end do

    end select
    
    ! write chunk terminator
    write(iunit, '(A)') ' </coords>'

  end subroutine
  
  subroutine fma_writeAdjacencyChunkEdge(iunit, h_Idx, n)
  integer, intent(in) :: iunit
  integer, intent(in) :: h_Idx
  integer, intent(in) :: n
  
  integer :: i,j
  integer, dimension(:,:), pointer :: p_Idx
  
    ! fetch the index array
    call storage_getbase_int2d(h_Idx, p_Idx)
    
    ! write chunk starter
    write(iunit,'(A)') ' <vert@edge>'
    
    do i = 1, n
      write(iunit,'(A)') '  ' // &
        trim(sys_sil(p_Idx(1,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(2,i)-1,16))
    end do
    
    ! write chunk terminator
    write(iunit,'(A)') ' </vert@edge>'
    
  end subroutine
  
  subroutine fma_writeAdjacencyChunkTria(iunit, h_Idx, n)
  integer, intent(in) :: iunit
  integer, intent(in) :: h_Idx
  integer, intent(in) :: n
  
  integer :: i,j
  integer, dimension(:,:), pointer :: p_Idx
  
    ! fetch the index array
    call storage_getbase_int2d(h_Idx, p_Idx)
    
    ! write chunk starter
    write(iunit,'(A)') ' <vert@tria>'
    
    do i = 1, n
      write(iunit,'(A)') '  ' // &
        trim(sys_sil(p_Idx(1,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(2,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(3,i)-1,16))
    end do
    
    ! write chunk terminator
    write(iunit,'(A)') ' </vert@tria>'
    
  end subroutine
  
  subroutine fma_writeAdjacencyChunkQuad(iunit, h_Idx, n)
  integer, intent(in) :: iunit
  integer, intent(in) :: h_Idx
  integer, intent(in) :: n
  
  integer :: i,j
  integer, dimension(:,:), pointer :: p_Idx
  
    ! fetch the index array
    call storage_getbase_int2d(h_Idx, p_Idx)
    
    ! write chunk starter
    write(iunit,'(A)') ' <vert@quad>'
    
    do i = 1, n
      write(iunit,'(A)') '  ' // &
        trim(sys_sil(p_Idx(1,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(2,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(4,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(3,i)-1,16))
    end do
    
    ! write chunk terminator
    write(iunit,'(A)') ' </vert@quad>'
    
  end subroutine
  
  subroutine fma_writeAdjacencyChunkTriaQuad(iunit, h_Idx, n)
  integer, intent(in) :: iunit
  integer, intent(in) :: h_Idx
  integer, intent(in) :: n
  
  integer :: i,j
  integer, dimension(:,:), pointer :: p_Idx
  
    ! fetch the index array
    call storage_getbase_int2d(h_Idx, p_Idx)
    
    ! write chunk starter
    write(iunit,'(A)') ' <vert@tria>'
    
    do i = 1, n
      if(p_Idx(4,i) .eq. 0) then
      write(iunit,'(A)') '  ' // &
        trim(sys_sil(p_Idx(1,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(2,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(3,i)-1,16))
      end if
    end do
    
    ! write chunk terminator
    write(iunit,'(A)') ' </vert@tria>'

        ! write chunk starter
    write(iunit,'(A)') ' <vert@quad>'
    
    do i = 1, n
      if(p_Idx(4,i) .ne. 0) then
      write(iunit,'(A)') '  ' // &
        trim(sys_sil(p_Idx(1,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(2,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(4,i)-1,16)) // ' ' // &
        trim(sys_sil(p_Idx(3,i)-1,16))
      end if
    end do
    
    ! write chunk terminator
    write(iunit,'(A)') ' </vert@quad>'

  end subroutine

end module
