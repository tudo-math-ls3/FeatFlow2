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
use random

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

  ! local variables
  type(t_random) :: rrng
  real(DP), dimension(:,:), pointer :: p_Dvtx
  real(DP) :: dhdist,dhdist2,dh,dtol,dx,dy
  integer :: i,ivt,invt
  integer(i32) :: iv, idx

    ! Compute inverse mesh width
    dh = sqrt(real(rtria%NVT,DP)) - 1.0_DP
    
    ! compute tolerance for boundary vertices
    dtol = 0.5_DP / dh
    
    ! Compute number of vertices per dimension
    invt = int(dh) + 1
    
    ! Compute distortion parameters
    dhdist = ddist / dh
    dhdist2 = dsecDist / dh
    if((dhdist .eq. 0.0_DP) .and. (dhdist2 .eq. 0.0_DP)) return
    
    ! Get arrays from the triangulation
    call storage_getbase_double2d (rtria%h_DvertexCoords, p_Dvtx)
    
    ! Loop over the vertices
    do ivt = 1, rtria%NVT

      ! Calculate line index
      idx = int(p_Dvtx(1,ivt) * dh) + 1

      ! seed the rng with the line index and advance it a few steps
      call rng_init(rrng, idx)
      do i = 1, mod(idx,7)
        call rng_advance(rrng)
      end do
      
      ! compute distortion
      call rng_get_int32(rrng, iv)
      dx = real(1 - iand(ishft(iv,-1), 2), DP)*dhdist
      dy = real(1 - iand(ishft(iv,-2), 2), DP)*dhdist2
      
      ! Distort the vertice's coordiantes.
      if((p_Dvtx(1,ivt) .gt. dtol) .and. (p_Dvtx(1,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(1,ivt) = p_Dvtx(1,ivt) + dx
      if((p_Dvtx(2,ivt) .gt. dtol) .and. (p_Dvtx(2,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(2,ivt) = p_Dvtx(2,ivt) + dy

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


  ! local variables
  type(t_random) :: rrng
  real(DP), dimension(:,:), pointer :: p_Dvtx
  real(DP) :: dhdist,dhdist2,dh,dtol,dx,dy
  integer :: i,ivt,invt
  integer(i32) :: iv, idx

    ! Compute inverse mesh width
    dh = sqrt(real(rtria%NVT,DP)) - 1.0_DP
    
    ! compute tolerance for boundary vertices
    dtol = 0.5_DP / dh
    
    ! Compute number of vertices per dimension
    invt = int(dh) + 1
    
    ! Compute distortion parameters
    dhdist = ddist / dh
    dhdist2 = dsecDist / dh
    if((dhdist .eq. 0.0_DP) .and. (dhdist2 .eq. 0.0_DP)) return
    
    ! Get arrays from the triangulation
    call storage_getbase_double2d (rtria%h_DvertexCoords, p_Dvtx)
    
    ! Loop over the vertices
    do ivt = 1, rtria%NVT

      ! Calculate line index
      idx = int(p_Dvtx(2,ivt) * dh) + 1

      ! seed the rng with the line index and advance it a few steps
      call rng_init(rrng, idx)
      do i = 1, mod(idx,7)
        call rng_advance(rrng)
      end do
      
      ! compute distortion
      call rng_get_int32(rrng, iv)
      dx = real(1 - iand(ishft(iv,-1), 2), DP)*dhdist2
      dy = real(1 - iand(ishft(iv,-2), 2), DP)*dhdist
      
      ! Distort the vertice's coordiantes.
      if((p_Dvtx(1,ivt) .gt. dtol) .and. (p_Dvtx(1,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(1,ivt) = p_Dvtx(1,ivt) + dx
      if((p_Dvtx(2,ivt) .gt. dtol) .and. (p_Dvtx(2,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(2,ivt) = p_Dvtx(2,ivt) + dy

    end do

  end subroutine

end module
