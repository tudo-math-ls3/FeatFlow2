!##############################################################################
!# ****************************************************************************
!# <Name> zpinch_meshadaptation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to perform
!# mesh adaptation for the time-dependent magnetohydrodynamic
!# equations in the one-, two- or three-dimensional domain $\Omega$.
!#
!# The following routines are available:
!#
!# 1.) zpinch_adaptTriangulation
!#     -> Performs h-adaptation for the given triangulation
!#
!#
!# Th following auxiliary routines are available:
!#
!# 1.) zpinch_hadaptCallbackScalar2d
!#     -> Performs application specific tasks in the adaptation
!#        algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 2.) zpinch_hadaptCallbackBlock2d
!#     -> Performs application specific tasks in the adaptation
!#        algorithm in 2D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module zpinch_meshadaptation

#include "../../flagship.h"
#include "hydro.h"

!$use omp_lib
  use basicgeometry
  use collection
  use flagship_basic
  use flagship_callback
  use fsystem
  use genoutput
  use hadaptaux
  use hadaptivity
  use zpinch_callback2d
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use triangulation

  implicit none

  private

  public :: zpinch_adaptTriangulation

contains

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_adaptTriangulation(rparlist, ssectionName, rhadapt,&
      rtriangulationSrc, rindicator, rcollection, rtriangulationDest)

!<description>
    ! This subroutine performs h-adaptation for the given triangulation
!</description>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! adaptation structure
    type(t_hadapt), intent(inout) :: rhadapt

    ! source triangulation structure
    type(t_triangulation), intent(inout), target :: rtriangulationSrc

    ! element-wise indicator
    type(t_vectorScalar), intent(inout) :: rindicator

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! OPTIONAL: destination triangulation structure
    ! If it is not given, the source triangulation is updated
    type(t_triangulation), intent(out), optional, target :: rtriangulationDest
!</output>
!</subroutine>


    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    integer :: isystemFormat


    ! Check if adaptation structure has been prepared
    if (rhadapt%iSpec .eq. HADAPT_HAS_PARAMETERS) then

      ! Initialise adaptation structure from triangulation
      call hadapt_initFromTriangulation(rhadapt, rtriangulationSrc)

    else

      ! Refresh adaptation structure
      call hadapt_refreshAdaptation(rhadapt, rtriangulationSrc)

    end if

    ! Get parameters from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'isystemformat', isystemFormat)

    ! What type of system format are we?
    select case (isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)

      if (rtriangulationSrc%ndim .eq. NDIM2D) then
        call zpinch_hadaptCallbackScalar2D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, zpinch_hadaptCallbackScalar2D)
      else
        call output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_adaptTriangulation')
        call sys_halt()
      end if

    case (SYSTEM_BLOCKFORMAT)

      if (rtriangulationSrc%ndim .eq. NDIM2D) then
        call zpinch_hadaptCallbackBlock2D(&
            HADAPT_OPR_INITCALLBACK, rcollection)
        call hadapt_performAdaptation(rhadapt, rindicator,&
            rcollection, zpinch_hadaptCallbackBlock2D)
      else
        call output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_adaptTriangulation')
        call sys_halt()
      end if

    case default

      call output_line('Invalid type of system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_adaptTriangulation')
      call sys_halt()
    end select


    ! Do we have a destination triangulation or should the source
    ! triangulation structure be updated?
    if (present(rtriangulationDest)) then
      p_rtriangulation => rtriangulationDest
    else
      p_rtriangulation => rtriangulationSrc
    end if

    ! Generate raw mesh from adaptation structure
    call hadapt_generateRawMesh(rhadapt, p_rtriangulation)

  end subroutine zpinch_adaptTriangulation

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_hadaptCallbackScalar2d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
    ! to be store in scalar interleave format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   IquickAccess(1):     NEQ or ivt
    !   IquickAccess(2:5):   ivt1,...,ivt5
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolutionHydro, rsolutionTransport
    real(DP), dimension(:), pointer, save :: p_DsolutionHydro, p_DsolutionTransport
    integer :: ivar


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vectors from colletion and set pointer
      rsolutionHydro     => rcollection%p_rvectorQuickAccess1
      rsolutionTransport => rcollection%p_rvectorQuickAccess2

      ! Check if solution is stored in interleave format
      if (rsolutionHydro%nblocks .ne. 1) then
        call output_line('Vector is not in interleave format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_hadaptCallbackScalar2d')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vectors
      nullify(rsolutionHydro, p_DsolutionHydro)
      nullify(rsolutionTransport, p_DsolutionTransport)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .ne. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro,&
            NVAR2D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if

      ! Resize solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .ne. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro,&
            NVAR2D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if
      do ivar = 1, NVAR2D
        p_DsolutionHydro((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = &
            0.5_DP*(p_DsolutionHydro((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)+&
                    p_DsolutionHydro((rcollection%IquickAccess(3)-1)*NVAR2D+ivar))
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(rcollection%IquickAccess(1)) =&
          0.5_DP*(p_DsolutionTransport(rcollection%IquickAccess(2))+&
                  p_DsolutionTransport(rcollection%IquickAccess(3)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)
      

    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro,&
            NVAR2D*rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if
      do ivar = 1, NVAR2D
        p_DsolutionHydro((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = &
            0.25_DP*(p_DsolutionHydro((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)+&
                     p_DsolutionHydro((rcollection%IquickAccess(3)-1)*NVAR2D+ivar)+&
                     p_DsolutionHydro((rcollection%IquickAccess(4)-1)*NVAR2D+ivar)+&
                     p_DsolutionHydro((rcollection%IquickAccess(5)-1)*NVAR2D+ivar))
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(rcollection%IquickAccess(1)) =&
          0.25_DP*(p_DsolutionTransport(rcollection%IquickAccess(2))+&
                   p_DsolutionTransport(rcollection%IquickAccess(3))+&
                   p_DsolutionTransport(rcollection%IquickAccess(4))+&
                   p_DsolutionTransport(rcollection%IquickAccess(5)))


      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution for the hydrodynamic model
      if (rcollection%IquickAccess(2) .ne. 0) then
        do ivar = 1, NVAR2D
          p_DsolutionHydro((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = &
              p_DsolutionHydro((rcollection%IquickAccess(2)-1)*NVAR2D+ivar)
        end do
      else
        do ivar = 1, NVAR2D
          p_DsolutionHydro((rcollection%IquickAccess(1)-1)*NVAR2D+ivar) = 0.0_DP
        end do
      end if

      ! Remove vertex from solution for the scalar transport model
      if (rcollection%IquickAccess(2) .ne. 0) then
        p_DsolutionTransport(rcollection%IquickAccess(1)) =&
            p_DsolutionTransport(rcollection%IquickAccess(2))
      else
        p_DsolutionTransport(rcollection%IquickAccess(1)) = 0.0_DP
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine zpinch_hadaptCallbackScalar2d

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_hadaptCallbackBlock2d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
    ! to be store in block format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   IquickAccess(1):     NEQ or ivt
    !   IquickAccess(2:5):   ivt1,...,ivt5
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolutionHydro, rsolutionTransport
    real(DP), dimension(:), pointer, save :: p_DsolutionHydro, p_DsolutionTransport
    integer :: ivar,neq


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vectors from colletion and set pointer
      rsolutionHydro     => rcollection%p_rvectorQuickAccess1
      rsolutionTransport => rcollection%p_rvectorQuickAccess2

      ! Check if solution is stored in interleave format
      if (rsolutionHydro%nblocks .ne. NVAR2D) then
        call output_line('Vector is not in block format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_hadaptCallbackBlock2d')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vectors
      nullify(rsolutionHydro, p_DsolutionHydro)
      nullify(rsolutionTransport, p_DsolutionTransport)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .ne. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro,&
            NVAR2D*rcollection%IquickAccess(1), .false., .true.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if

      ! Resize solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .ne. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro, NVAR2D&
            *rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if
      neq = rsolutionHydro%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(1)) = &
            0.5_DP*(p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(2))+&
                    p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(3)) )
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(rcollection%IquickAccess(1)) =&
          0.5_DP*(p_DsolutionTransport(rcollection%IquickAccess(2))+&
                  p_DsolutionTransport(rcollection%IquickAccess(3)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector for the hydrodynamic model
      if (rsolutionHydro%NEQ .lt. NVAR2D*rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionHydro, NVAR2D&
            *rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionHydro, p_DsolutionHydro)
      end if
      neq = rsolutionHydro%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(1)) =&
            0.25_DP*(p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(2))+&
                     p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(3))+&
                     p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(4))+&
                     p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(5)) )
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(rcollection%IquickAccess(1)) =&
          0.25_DP*(p_DsolutionTransport(rcollection%IquickAccess(2))+&
                   p_DsolutionTransport(rcollection%IquickAccess(3))+&
                   p_DsolutionTransport(rcollection%IquickAccess(4))+&
                   p_DsolutionTransport(rcollection%IquickAccess(5)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution for the hydrodynamic model
      if (rcollection%IquickAccess(2) .ne. 0) then
        neq = rsolutionHydro%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(1)) = &
              p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(2))
        end do
      else
        neq = rsolutionHydro%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_DsolutionHydro((ivar-1)*neq+rcollection%IquickAccess(1)) = 0.0_DP
        end do
      end if

      ! Remove vertex from solution for the scalar transport model
      if (rcollection%IquickAccess(2) .ne. 0) then
        p_DsolutionTransport(rcollection%IquickAccess(1)) =&
            p_DsolutionTransport(rcollection%IquickAccess(2))
      else
        p_DsolutionTransport(rcollection%IquickAccess(1)) = 0.0_DP
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine zpinch_hadaptCallbackBlock2d

end module zpinch_meshadaptation
