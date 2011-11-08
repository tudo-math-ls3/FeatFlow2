!##############################################################################
!# ****************************************************************************
!# <name> disto2d_aux </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a hand full of auxiliary routines for the
!# elemdbg2d_testX debugger modules.
!#
!# The following routines can be found in this module:
!#
!# .) disto2d_distortQuadLineX
!#    -> Distorts the vertices of a QUAD mesh by X-line-wise index-based
!#       stochastical distortion.
!#
!# .) disto2d_distortQuadLineY
!#    -> Distorts the vertices of a QUAD mesh by Y-line-wise index-based
!#       stochastical distortion.
!#
!# </purpose>
!##############################################################################

module disto2d_aux

use fsystem
use genoutput
use storage
use cubature
use transformation
use triangulation

implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine disto2d_distortQuadLineX(rtria, ddist, dsecDist)

!<description>
  ! Distorts the vertices of a QUAD mesh by X-line-wise index-based
  ! stochastical distortion.
!</description>

!<input>
  ! The distortion parameters
  real(DP), intent(IN) :: ddist, dsecDist
!</input>

!<inputoutput>
  ! The CUBE mesh that is to be distorted.
  type(t_triangulation), intent(INOUT) :: rtria
!</inputoutput>

!</subroutine>

  real(DP), parameter :: dtol = 1e-8_DP

  ! local variables
  real(DP), dimension(:,:), pointer :: p_Dvtx
  real(DP) :: dhdist,dhdist2,dh2
  integer :: ivt,iX,iY,invt

    ! Calculate number of vertices in each dimension
    dh2 = sqrt(real(rtria%NVT,DP))
    
    ! Get number of vertices
    invt = int(dh2) + 10
    
    ! Calculate distortion parameters
    dhdist = ddist / (dh2 + 1.0_DP)
    dhdist2 = dsecDist / (dh2 + 1.0_DP)
    if((dhdist .eq. 0.0_DP) .and. (dhdist2 .eq. 0.0_DP)) return
    
    ! Get arrays from the triangulation
    call storage_getbase_double2d (rtria%h_DvertexCoords, p_Dvtx)
    
    ! Loop over the vertices
    do ivt = 1, rtria%NVT
    
      ! Calculate the X- and Y-position of this vertice
      iX = int(p_Dvtx(1,ivt) * dh2) + 1
      iY = int(p_Dvtx(2,ivt) * dh2) + 1
      
      ! Distort the vertice's coordiantes.
      if((p_Dvtx(1,ivt) .gt. dtol) .and. (p_Dvtx(1,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(1,ivt) = p_Dvtx(1,ivt) + real((-1)**mod(iX,17),DP)*dhdist
      if((p_Dvtx(2,ivt) .gt. dtol) .and. (p_Dvtx(2,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(2,ivt) = p_Dvtx(2,ivt) + real((-1)**mod(ivt+7,17),DP)*dhdist2
        
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine disto2d_distortQuadLineY(rtria, ddist, dsecDist)

!<description>
  ! Distorts the vertices of a QUAD mesh by Y-line-wise index-based
  ! stochastical distortion.
!</description>

!<input>
  ! The distortion parameters
  real(DP), intent(IN) :: ddist, dsecDist
!</input>

!<inputoutput>
  ! The CUBE mesh that is to be distorted.
  type(t_triangulation), intent(INOUT) :: rtria
!</inputoutput>

!</subroutine>

  real(DP), parameter :: dtol = 1e-8_DP

  ! local variables
  real(DP), dimension(:,:), pointer :: p_Dvtx
  real(DP) :: dhdist,dhdist2,dh2
  integer :: ivt,iX,iY,invt

    ! Calculate number of vertices in each dimension
    dh2 = sqrt(real(rtria%NVT,DP))
    
    ! Get number of vertices
    invt = int(dh2) + 10
    
    ! Calculate distortion parameters
    dhdist = ddist / (dh2 + 1.0_DP)
    dhdist2 = dsecDist / (dh2 + 1.0_DP)
    if((dhdist .eq. 0.0_DP) .and. (dhdist2 .eq. 0.0_DP)) return
    
    ! Get arrays from the triangulation
    call storage_getbase_double2d (rtria%h_DvertexCoords, p_Dvtx)
    
    ! Loop over the vertices
    do ivt = 1, rtria%NVT
    
      ! Calculate the X- and Y-position of this vertice
      iX = int(p_Dvtx(1,ivt) * dh2) + 1
      iY = int(p_Dvtx(2,ivt) * dh2) + 1
      
      ! Distort the vertice's coordiantes.
      if((p_Dvtx(1,ivt) .gt. dtol) .and. (p_Dvtx(1,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(1,ivt) = p_Dvtx(1,ivt) + real((-1)**mod(ivt+7,17),DP)*dhdist2
      if((p_Dvtx(2,ivt) .gt. dtol) .and. (p_Dvtx(2,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(2,ivt) = p_Dvtx(2,ivt) + real((-1)**mod(iX,17),DP)*dhdist
        
    end do

  end subroutine

end module