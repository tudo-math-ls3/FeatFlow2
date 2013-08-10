#include "kernel/feat2macros.h" 
#include "kernel/template.h" 
#include "kernel/PDEOperators/afc.h"

    !**************************************************************
    ! Assemble the local bounds from the predicted solution

    subroutine FEAT2_PP_TEMPLATE2(doBounds,__VectorName__,__AFCName__)(IedgeListIdx, IedgeList, NEDGE, Dx, Dqp, Dqm)

      ! input parameters
      __VectorData__, dimension(:), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! The local upper/lower bounds computed from Dx
      __AFCData__, dimension(:), intent(out) :: Dqp,Dqm

      ! local variables
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Initialise Q`s by solution
      !$omp sections
      !$omp section
      call lalg_copyVector(Dx, Dqp)
      !$omp section
      call lalg_copyVector(Dx, Dqm)
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

          ! Compute local upper and lower bounds
#if TemplateType_Vector == TemplateType_AFC
          Dqp(i) = max(Dqp(i), Dx(j))
          Dqm(i) = min(Dqm(i), Dx(j))
          Dqp(j) = max(Dqp(j), Dx(i))
          Dqm(j) = min(Dqm(j), Dx(i))
#else
          Dqp(i) = max(Dqp(i), real(Dx(j),__AFCType__))
          Dqm(i) = min(Dqm(i), real(Dx(j),__AFCType__))
          Dqp(j) = max(Dqp(j), real(Dx(i),__AFCType__))
          Dqm(j) = min(Dqm(j), real(Dx(i),__AFCType__))
#endif
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine
