#include "../feat2macros.h" 
#include "../template.h" 
#include "afc.h"

    !**************************************************************
    ! Compute nodal correction factors with constraints

    subroutine FEAT2_PP_TEMPLATE2(doLimitNodalConstrained,__VectorName__,__AFCName__)(NEQ, dscale,&
        ML, Dx, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      ! input parameters
      __VectorData__, dimension(:), intent(in) :: ML,Dx
      __AFCData__, dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ

      ! output parameters
      __AFCData__, dimension(:), intent(out) :: Drp,Drm

      ! local variables
      integer :: ieq

      !$omp parallel sections default(shared) private(ieq)

      !$omp section

      !$omp parallel do default(shared) private(ieq)
      do ieq = 1, NEQ
        Drp(ieq) = min(FEAT2_PP_CONST(1.0,TemplateType_AFC), (Dqp(ieq)-Dx(ieq)+AFCSTAB_EPSABS) /&
                               (Dpp(ieq)+AFCSTAB_EPSABS) * (ML(ieq)/dscale))
      end do
      !$omp end parallel do

      !$omp section

      !$omp parallel do default(shared) private(ieq)
      do ieq = 1, NEQ
        Drm(ieq) = min(FEAT2_PP_CONST(1.0,TemplateType_AFC), (Dqm(ieq)-Dx(ieq)-AFCSTAB_EPSABS) /&
                               (Dpm(ieq)-AFCSTAB_EPSABS) * (ML(ieq)/dscale))
      end do
      !$omp end parallel do

      !$omp end parallel sections

    end subroutine
