!##############################################################################
!# ****************************************************************************
!# <name> lumping </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises different methods to create lumped mass matrices.
!#
!# The following routines can be found here:
!#
!# 1.) lump_extractDiagonal
!#     -> Lumps an existing mass matrix by extracting the diagonal entries
!#
!# 1.) lump_sumToDiagonal
!#     -> Lumps an existing mass matrix summing up the values to the diagonal
!# </purpose>
!##############################################################################

module lumping

  use fsystem
  use genoutput
  use storage
  use linearsystemscalar
  
  implicit none
  private
  
  public :: lump_extractDiagonal
  public :: lump_sumToDiagonal

contains

  !****************************************************************************

!<subroutine>

  subroutine lump_extractDiagonal(rmatrixSource,rmatrixDest)

!<description>  
  ! Creates a lumped (mass) matrix by extracting the diagonal of a source matrix.
  !
  ! The destination matrix shall be empty (filled by zero)!
!</description>

!<input>
  ! Source matrix.
  type(t_matrixScalar), intent(in) :: rmatrixSource
!</input>

!<inputoutput>
  ! Destination matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixDest
!</inputoutput>

!</subroutine>
    integer, dimension(:), pointer :: p_KdiagSource,p_KdiagDest
    real(DP), dimension(:), pointer :: p_Dsource, p_Ddest
    integer :: irow

    if (.not. lsyssc_hasMatrixContent(rmatrixDest)) then
      ! Auto-allocate memory.
      call lsyssc_allocEmptyMatrix (rmatrixDest,.true.)
    end if

    if (rmatrixSource%cdataType .ne. rmatrixDest%cdataType) then
      call output_line("Source and destination matrix have different data type",&
          OU_CLASS_ERROR,OU_MODE_STD,"lump_extractDiagonal")
      call sys_halt()
    end if

    if (rmatrixSource%cdataType .ne. ST_DOUBLE) then
      call output_line("Currently, only double precision data is supported",&
          OU_CLASS_ERROR,OU_MODE_STD,"lump_extractDiagonal")
      call sys_halt()
    end if

    ! What is the source matrix format?
    select case (rmatrixSource%cmatrixFormat)
    case (LSYSSC_MATRIX9)

      ! Get the structure/data
      call lsyssc_getbase_Kdiagonal (rmatrixSource,p_KdiagSource)
      call lsyssc_getbase_double (rmatrixSource,p_Dsource)
      
      call lsyssc_getbase_double (rmatrixDest,p_Ddest)
    
      ! Destination matrix format?
      select case (rmatrixDest%cmatrixFormat)
      case (LSYSSC_MATRIX9)
      
        call lsyssc_getbase_Kdiagonal (rmatrixDest,p_KdiagDest)
        
        ! Extract the entries        
        do irow = 1,rmatrixSource%NEQ
          p_Ddest(p_KdiagDest(irow)) = &
              p_Ddest(p_KdiagDest(irow)) + p_Dsource(p_KdiagSource(irow))
        end do
        
      case (LSYSSC_MATRIXD)
      
        ! Extract the entries        
        do irow = 1,rmatrixSource%NEQ
          p_Ddest(irow) = p_Ddest(irow) + p_Dsource(p_KdiagSource(irow))
        end do

      case (LSYSSC_MATRIX1)
      
        ! Extract the entries        
        do irow = 1,rmatrixSource%NEQ
          p_Ddest((irow-1)*rmatrixSource%NEQ+irow) = &
              p_Ddest((irow-1)*rmatrixSource%NEQ+irow) + p_Dsource(p_KdiagSource(irow))
        end do

      case (LSYSSC_MATRIX7)
      
        call lsyssc_getbase_Kld (rmatrixDest,p_KdiagDest)
        
        ! Extract the entries        
        do irow = 1,rmatrixSource%NEQ
          p_Ddest(p_KdiagDest(irow)) = &
              p_Ddest(p_KdiagDest(irow)) + p_Dsource(p_KdiagSource(irow))
        end do

      case default
        call output_line("Unsupported matrix format",&
            OU_CLASS_ERROR,OU_MODE_STD,"lump_extractDiagonal")
        call sys_halt()

      end select
      
    case (LSYSSC_MATRIX7)

      ! Get the structure/data
      call lsyssc_getbase_Kld (rmatrixSource,p_KdiagSource)
      call lsyssc_getbase_double (rmatrixSource,p_Dsource)

      call lsyssc_getbase_double (rmatrixDest,p_Ddest)
    
      ! Destination matrix format?
      select case (rmatrixDest%cmatrixFormat)
      case (LSYSSC_MATRIX9)
      
        call lsyssc_getbase_Kdiagonal (rmatrixDest,p_KdiagDest)
        
        ! Extract the entries        
        do irow = 1,rmatrixSource%NEQ
          p_Ddest(p_KdiagDest(irow)) = &
              p_Ddest(p_KdiagDest(irow)) + p_Dsource(p_KdiagSource(irow))
        end do
        
      case (LSYSSC_MATRIXD)
      
        ! Extract the entries        
        do irow = 1,rmatrixSource%NEQ
          p_Ddest(irow) = p_Ddest(irow) + p_Dsource(p_KdiagSource(irow))
        end do

      case (LSYSSC_MATRIX1)
      
        ! Extract the entries        
        do irow = 1,rmatrixSource%NEQ
          p_Ddest((irow-1)*rmatrixSource%NEQ+irow) = &
              p_Ddest((irow-1)*rmatrixSource%NEQ+irow) + p_Dsource(p_KdiagSource(irow))
        end do

      case (LSYSSC_MATRIX7)
      
        call lsyssc_getbase_Kld (rmatrixDest,p_KdiagDest)
        
        ! Extract the entries        
        do irow = 1,rmatrixSource%NEQ
          p_Ddest(p_KdiagDest(irow)) = &
              p_Ddest(p_KdiagDest(irow)) + p_Dsource(p_KdiagSource(irow))
        end do
        
      case default
        call output_line("Unsupported matrix format",&
            OU_CLASS_ERROR,OU_MODE_STD,"lump_extractDiagonal")
        call sys_halt()

      end select

    case default
      call output_line("Unsupported matrix format",&
          OU_CLASS_ERROR,OU_MODE_STD,"lump_extractDiagonal")
      call sys_halt()

    end select
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lump_sumToDiagonal(rmatrixSource,rmatrixDest)

!<description>  
  ! Creates a lumped (mass) matrix by summing up the entries of a row to
  ! the diagonal of the destination matrix.
  !
  ! The destination matrix shall be empty (filled by zero)!
!</description>

!<input>
  ! Source matrix.
  type(t_matrixScalar), intent(in) :: rmatrixSource
!</input>

!<inputoutput>
  ! Destination matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixDest
!</inputoutput>

!</subroutine>
    integer, dimension(:), pointer :: p_KdiagDest
    integer, dimension(:), pointer :: p_Kld
    real(DP), dimension(:), pointer :: p_Dsource, p_Ddest
    real(DP) :: dval
    integer :: irow,icol

    if (.not. lsyssc_hasMatrixContent(rmatrixDest)) then
      ! Auto-allocate memory.
      call lsyssc_allocEmptyMatrix (rmatrixDest,.true.)
    end if

    if (rmatrixSource%cdataType .ne. rmatrixDest%cdataType) then
      call output_line("Source and destination matrix have different data type",&
          OU_CLASS_ERROR,OU_MODE_STD,"lump_sumToDiagonal")
      call sys_halt()
    end if

    if (rmatrixSource%cdataType .ne. ST_DOUBLE) then
      call output_line("Currently, only double precision data is supported",&
          OU_CLASS_ERROR,OU_MODE_STD,"lump_sumToDiagonal")
      call sys_halt()
    end if

    ! What is the source matrix format?
    select case (rmatrixSource%cmatrixFormat)
    case (LSYSSC_MATRIX9,LSYSSC_MATRIX7)

      ! Get the structure/data
      call lsyssc_getbase_Kld (rmatrixSource,p_Kld)
      call lsyssc_getbase_double (rmatrixSource,p_Dsource)
      
      call lsyssc_getbase_double (rmatrixDest,p_Ddest)
    
      ! Destination matrix format?
      select case (rmatrixDest%cmatrixFormat)
      case (LSYSSC_MATRIX9)
      
        call lsyssc_getbase_Kdiagonal (rmatrixDest,p_KdiagDest)
        
        ! Add all off-diagonals to the diagonal
        do irow = 1,rmatrixSource%NEQ
          dval = 0.0_DP
          do icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            dval = dval + p_Dsource(icol)
          end do
          p_Ddest(p_KdiagDest(irow)) = p_Ddest(p_KdiagDest(irow)) + dval
        end do
        
      case (LSYSSC_MATRIXD)
      
        ! Add all off-diagonals to the diagonal
        do irow = 1,rmatrixSource%NEQ
          dval = 0.0_DP
          do icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            dval = dval + p_Dsource(icol)
          end do
          p_Ddest(irow) = p_Ddest(irow) + dval
        end do

      case (LSYSSC_MATRIX1)
      
        ! Add all off-diagonals to the diagonal
        do irow = 1,rmatrixSource%NEQ
          dval = 0.0_DP
          do icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            dval = dval + p_Dsource(icol)
          end do
          p_Ddest((irow-1)*rmatrixSource%NEQ+irow) = &
              p_Ddest((irow-1)*rmatrixSource%NEQ+irow) + dval
        end do

      case (LSYSSC_MATRIX7)
      
        call lsyssc_getbase_Kld (rmatrixDest,p_KdiagDest)
        
        ! Add all off-diagonals to the diagonal
        do irow = 1,rmatrixSource%NEQ
          dval = 0.0_DP
          do icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            dval = dval + p_Dsource(icol)
          end do
          p_Ddest(p_KdiagDest(irow)) = p_Ddest(p_KdiagDest(irow)) + dval
        end do

      case default
        call output_line("Unsupported matrix format",&
            OU_CLASS_ERROR,OU_MODE_STD,"lump_extractDiagonal")
        call sys_halt()

      end select
      
    case default
      call output_line("Unsupported matrix format",&
          OU_CLASS_ERROR,OU_MODE_STD,"lump_extractDiagonal")
      call sys_halt()

    end select

  end subroutine

end module
