!##############################################################################
!# ****************************************************************************
!# <name> interpolation </name>
!# ****************************************************************************
!#
!# <purpose>
!#   This module contains some interpolation routines
!# </purpose>
!##############################################################################

MODULE interpolation

  USE fsystem
  USE triangulation

IMPLICIT NONE

CONTAINS

!<subroutine>

  SUBROUTINE Interpol_Press(DpEl,DpMid,rtriangulation,ipar)

  !<description>
    ! This subroutine enables the interpolation of the pressure
    ! from piecewise constant approximation to piecewice linear
    ! and back, depending on the ipar value
    ! ipar == 0 : convert piecewise constant DPC to piecewise linear DPL
    ! ipar == 1 : convert piecewise linear DPL to piecewise constant DPC
  !</description>

  !<input>
    TYPE (t_triangulation), INTENT(IN) :: rtriangulation
    REAL(DP), INTENT(INOUT), DIMENSION(:) :: DpEl,DpMid
    INTEGER(I32),INTENT(IN) ::  ipar
  !</input>
  
  !<output>

  !</output>

!</subroutine>

  ! local variables
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_Kmid, p_Kadj

  INTEGER(I32) :: iel,iadj,imid,ive
  REAL(DP) :: dph

  CALL storage_getbase_int2D(rtriangulation%h_IedgesAtElement,p_Kmid)
  CALL storage_getbase_int2D(rtriangulation%h_IneighboursAtElement,p_Kadj)

  IF (ipar == 0) THEN

    DO iel=1,rtriangulation%NEL
      dph = DpEl(iel)

      DO ive=1,4
        iadj=p_Kadj(ive,iel)
        imid=p_Kmid(ive,iel)-rtriangulation%NVT

        IF (iadj == 0) DpMid(imid)=dph
        IF (iadj >iel) DpMid(imid)=0.5d0*(dph+DpEl(iadj))

      END DO
    END DO

  ELSE

    DO iel=1,rtriangulation%NEL
      DpEl(iel)=0.25d0*(DpMid(p_Kmid(1,iel)-rtriangulation%NVT)+ &
                        DpMid(p_Kmid(2,iel)-rtriangulation%NVT)+ &
                        DpMid(p_Kmid(3,iel)-rtriangulation%NVT)+ &
                        DpMid(p_Kmid(4,iel)-rtriangulation%NVT))
    END DO

  ENDIF  

  END SUBROUTINE Interpol_Press

END MODULE

