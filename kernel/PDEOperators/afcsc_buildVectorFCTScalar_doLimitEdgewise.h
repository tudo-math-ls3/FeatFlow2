#include "../feat2macros.h" 
#include "../template.h" 
#include "afc.h"

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes

    subroutine FEAT2_PP_TEMPLATE1(doLimitEdgewise,__AFCName__)(IedgeList, NEDGE,&
        Dflux, Drp, Drm, Dalpha)

      ! input parameters
      __AFCData__, dimension(:), intent(in) :: Dflux
      __AFCData__, dimension(:), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factors from previous
      !           multiplicative correction steps
      ! On exit: the edge-wise correction factors resulting from
      !          the nodal correction factors Rp and Rm
      __AFCData__, dimension(:), intent(inout) :: Dalpha

      ! local variables
      __AFCData__ :: f_ij,r_ij
      integer :: iedge,i,j

      ! Loop over all edges
      !$omp parallel do default(shared) private(i,j,f_ij,r_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Get precomputed raw antidiffusive fluxes
        f_ij = Dflux(iedge)

        ! Compute nodal correction factors
        if (f_ij .gt. AFCSTAB_EPSABS) then
          r_ij = min(Drp(i),Drm(j))
        elseif (f_ij .lt. -AFCSTAB_EPSABS) then
          r_ij = min(Drp(j),Drm(i))
        else
          r_ij = FEAT2_PP_CONST(1.0,TemplateType_AFC)
        end if

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * r_ij
      end do
      !$omp end parallel do

    end subroutine
