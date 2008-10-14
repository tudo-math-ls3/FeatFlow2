!##############################################################################
!# ****************************************************************************
!# <name> interpolation </name>
!# ****************************************************************************
!#
!# <purpose>
!#   This module contains some interpolation routines
!# </purpose>
!##############################################################################

module interpolation

  use fsystem
  use triangulation

implicit none

contains

!<subroutine>

  subroutine Interpol_Press(DpEl,DpMid,rtriangulation,ipar)

  !<description>
    ! This subroutine enables the interpolation of the pressure
    ! from piecewise constant approximation to piecewice linear
    ! and back, depending on the ipar value
    ! ipar == 0 : convert piecewise constant DPC to piecewise linear DPL
    ! ipar == 1 : convert piecewise linear DPL to piecewise constant DPC
  !</description>

  !<input>
    type (t_triangulation), intent(IN) :: rtriangulation
    real(DP), intent(INOUT), dimension(:) :: DpEl,DpMid
    integer(I32),intent(IN) ::  ipar
  !</input>
  
  !<output>

  !</output>

!</subroutine>

  ! local variables
  integer(I32), dimension(:,:), pointer :: p_Kmid, p_Kadj

  integer(I32) :: iel,iadj,imid,ive
  real(DP) :: dph

  call storage_getbase_int2D(rtriangulation%h_IedgesAtElement,p_Kmid)
  call storage_getbase_int2D(rtriangulation%h_IneighboursAtElement,p_Kadj)

  if (ipar == 0) then

    do iel=1,rtriangulation%NEL
      dph = DpEl(iel)

      do ive=1,4
        iadj=p_Kadj(ive,iel)
        imid=p_Kmid(ive,iel)-rtriangulation%NVT

        if (iadj == 0) DpMid(imid)=dph
        if (iadj >iel) DpMid(imid)=0.5d0*(dph+DpEl(iadj))

      end do
    end do

  else

    do iel=1,rtriangulation%NEL
      DpEl(iel)=0.25d0*(DpMid(p_Kmid(1,iel)-rtriangulation%NVT)+ &
                        DpMid(p_Kmid(2,iel)-rtriangulation%NVT)+ &
                        DpMid(p_Kmid(3,iel)-rtriangulation%NVT)+ &
                        DpMid(p_Kmid(4,iel)-rtriangulation%NVT))
    end do

  endif  

  end subroutine Interpol_Press

end module

