#include "../feat2macros.h" 
#include "../template.h" 
#include "afc.h"

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the CSR-matrix stored in Dmatrix


	subroutine FEAT2_PP_TEMPLATE3(doFluxesByMatrix,__MatrixName__,__VectorName__,__AFCName__)(IedgeList, NEDGE,&
        Dmatrix, Dx, dscale, bclear, Dflux)

      ! input parameters
	  __MatrixData__, dimension(:), intent(in) :: Dmatrix
	  __VectorData__, dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear

      ! input/output parameters
     __AFCData__, dimension(:), intent(inout) :: Dflux

      ! local variables
      integer :: iedge,i,j,ij

      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux, NEDGE)

      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)

            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)

            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) + Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if

      elseif (dscale .eq. -1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)

            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dmatrix(ij) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)

            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) + Dmatrix(ij) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        end if

      else

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)

            ! Compute raw antidiffusive flux
            Dflux(iedge) = dscale * Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)

            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + dscale * Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if

      end if

    end subroutine
