#include "../feat2macros.h" 
#include "../template.h" 
#include "afc.h"

    !**************************************************************
    ! Assemble solution difference used for classical prelimiting.

	subroutine FEAT2_PP_TEMPLATE2(doDifferences,__VectorName__,__AFCName__)(IedgeList, NEDGE, Dx, Dflux)

      ! input parameters
      __VectorData__, dimension(:), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE

      ! output parameters
      __AFCData__, dimension(:), intent(out) :: Dflux

      ! local variables
      integer :: iedge,i,j

      !$omp parallel do default(shared) private(i,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Determine indices
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute solution difference; in contrast to the literature,
        ! we compute the solution difference $u_i-u_j$ and check if
        ! $f_{ij}(u_i-u_j)<0$ in the prelimiting step.
        Dflux(iedge) = Dx(i)-Dx(j)
      end do
      !$omp end parallel do

    end subroutine
