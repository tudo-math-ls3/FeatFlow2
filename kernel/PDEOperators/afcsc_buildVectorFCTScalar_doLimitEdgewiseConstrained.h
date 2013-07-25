#include "../feat2macros.h" 
#include "../template.h" 
#include "afc.h"

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes

    subroutine FEAT2_PP_TEMPLATE1(doLimitEdgewiseConstrained,__AFCName__)(IedgeList, NEDGE,&
        Dflux1, Dflux2, Drp, Drm, Dalpha)

      ! input parameters
      __AFCData__, dimension(:), intent(in) :: Dflux1,Dflux2
      __AFCData__, dimension(:), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE

      ! input/output parameters
      __AFCData__, dimension(:), intent(inout) :: Dalpha

      ! local variables
      __AFCData__ :: f1_ij,f2_ij,r_ij
      integer :: iedge,i,j

      ! Loop over all edges
      !$omp parallel do default(shared) private(i,j,f1_ij,f2_ij,r_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Get precomputed raw antidiffusive fluxes
        f1_ij = Dflux1(iedge)
        f2_ij = Dflux2(iedge)

        ! Compute nodal correction factors
        if (f1_ij*f2_ij .le. FEAT2_PP_CONST(0.0,TemplateType_AFC)) then
          r_ij = FEAT2_PP_CONST(0.0,TemplateType_AFC)
        else
          if (f1_ij .ge. FEAT2_PP_CONST(0.0,TemplateType_AFC)) then
            r_ij = min(FEAT2_PP_CONST(1.0,TemplateType_AFC), f1_ij/f2_ij*min(Drp(i),Drm(j)))
          else
            r_ij = min(FEAT2_PP_CONST(1.0,TemplateType_AFC), f1_ij/f2_ij*min(Drp(j),Drm(i)))
          end if
        end if

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * r_ij
      end do
      !$omp end parallel do

    end subroutine
