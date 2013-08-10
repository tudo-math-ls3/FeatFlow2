#include "kernel/feat2macros.h" 
#include "kernel/template.h" 
#include "kernel/PDEOperators/afc.h"

    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes the standard way, as
    ! suggested by Boris and Book in their first FCT algorithm

    subroutine FEAT2_PP_TEMPLATE1(doStdPrelimit,__AFCName__)(IedgeListIdx,&
        IedgeList, NEDGE, Dflux, DfluxPrel, Dalpha)

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

        ! Check if the antidiffusive flux is directed down the gradient
        !   $f_ij*(u_i-u_j) < 0$
        ! and if its magnitude is larger than an absolute tolerance
        !  $ |f_ij| > tol$
        ! In this case, cancel the flux completely.
        if ((Dflux(iedge)*DfluxPrel(iedge) .lt. FEAT2_PP_CONST(0.0,TemplateType_AFC)) .and.&
            abs(Dflux(iedge)) .gt. AFCSTAB_PRELIMABS)&
            Dalpha(iedge) = FEAT2_PP_CONST(0.0,TemplateType_AFC)
      end do
      !$omp end parallel do

    end subroutine
