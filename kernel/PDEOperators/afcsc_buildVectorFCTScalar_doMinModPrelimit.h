#include "../feat2macros.h"
#include "../template.h"
#include "afc.h"

    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes using minmod limiter


    subroutine FEAT2_PP_TEMPLATE1(doMinModPrelimit,__AFCName__)(IedgeListIdx, IedgeList,&
        NEDGE, Dflux, DfluxPrel, Dalpha)

      ! input parameters
      __AFCData__, dimension(:), intent(in) :: Dflux,DfluxPrel
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      __AFCData__, dimension(:), intent(inout) :: Dalpha

      ! local variables
      integer :: iedge

      ! Loop over all edges
      !$omp parallel do default(shared)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Check if the magnitude of the antidiffusive flux is larger
        ! than an absolute tolerance; otherwise no prelimiting is done
        if (abs(Dflux(iedge)) .gt. AFCSTAB_PRELIMABS) then
          ! Check if the antidiffusive flux is directed down the gradient
          !   $f_ij*fp_ij < 0$
          if (Dflux(iedge)*DfluxPrel(iedge) .lt. FEAT2_PP_CONST(0.0,TemplateType_AFC)) then
            ! Then, cancel the antidiffusive flux completely
            Dalpha(iedge) = FEAT2_PP_CONST(0.0,TemplateType_AFC)
          elseif (abs(Dflux(iedge)) .gt. abs(DfluxPrel(iedge))) then
            ! Check if the magnitude of the raw antidiffusive flux
            ! exceeds the magnitude of the prelimiting flux
            !   $|f_ij| > |fp_ij|$
            ! then set the correction factor as follows
            Dalpha(iedge) = min(Dalpha(iedge),DfluxPrel(iedge)/Dflux(iedge))
          end if
        end if
      end do
      !$omp end parallel do

      end subroutine
