#include "../feat2macros.h" 
#include "../template.h" 
#include "afc.h"

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes without prelimiting

    subroutine FEAT2_PP_TEMPLATE1(doADIncrements,__AFCName__)(IedgeListIdx, IedgeList,&
        NEDGE, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      __AFCData__, dimension(:), intent(in) :: Dflux
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      __AFCData__, dimension(:), intent(inout) :: Dalpha

      ! The sums of positive/negative antidiffusive increments
      __AFCData__, dimension(:), intent(out) :: Dpp,Dpm

      ! local variables
      __AFCData__ :: f_ij,fp_ij,fm_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,f_ij,fp_ij,fm_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear P`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp end sections

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Get node numbers
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Apply multiplicative correction factor
          f_ij = Dalpha(iedge) * Dflux(iedge)

          ! Separate fluxes into positive/negative contributions
          fp_ij = max(FEAT2_PP_CONST(0.0,__AFCType__),f_ij)
          fm_ij = min(FEAT2_PP_CONST(0.0,__AFCType__),f_ij)

          ! Compute the sums of antidiffusive increments
          Dpp(i) = Dpp(i) + fp_ij   ! += max(FEAT2_PP_CONST(0.0,__AFCType__), f_ij)
          Dpp(j) = Dpp(j) - fm_ij   ! += max(FEAT2_PP_CONST(0.0,__AFCType__),-f_ij)
          Dpm(i) = Dpm(i) + fm_ij   ! += min(FEAT2_PP_CONST(0.0,__AFCType__), f_ij)
          Dpm(j) = Dpm(j) - fp_ij   ! += min(FEAT2_PP_CONST(0.0,__AFCType__),-f_ij)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine
