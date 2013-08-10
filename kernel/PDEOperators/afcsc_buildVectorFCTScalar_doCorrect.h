#include "kernel/feat2macros.h" 
#include "kernel/template.h" 
#include "kernel/PDEOperators/afc.h"

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine FEAT2_PP_TEMPLATE2(doCorrect,__VectorName__,__AFCName__)(IedgeListIdx, IedgeList,&
        NEDGE, dscale, Dalpha, Dflux, Dy)

      ! input parameters
      __AFCData__, dimension(:), intent(in) :: Dalpha,Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! input/output parameters
      __VectorData__, dimension(:), intent(inout) :: Dy

      ! local variables
      __AFCData__ :: f_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,f_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

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

          ! Correct antidiffusive flux
          f_ij = dscale * Dalpha(iedge) * Dflux(iedge)

          ! Apply limited antidiffusive fluxes
          Dy(i) = Dy(i) + f_ij
          Dy(j) = Dy(j) - f_ij
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine
