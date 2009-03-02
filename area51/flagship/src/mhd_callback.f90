!##############################################################################
!# ****************************************************************************
!# <name> mhd_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the simplified MHD equations in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) mhd_nlsolverCallback
!#     -> Callback routine for the nonlinear solver
!#
!# 2.) mhd_calcVelocityField
!#     -> Calculates the velocity field for the scalar transport model
!#
!# 3.) mhd_hadaptCallbackScalar2d
!#     -> Performs application specific tasks in the adaptation
!#        algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 4.) mhd_hadaptCallbackBlock2d
!#     -> Performs application specific tasks in the adaptation
!#        algorithm in 2D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module mhd_callback

  use collection
  use euler_basic
  use euler_callback
  use flagship_basic
  use flagship_callback
  use fsystem
  use genoutput
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solveraux
  use statistics
  use storage
  use timestepaux

  implicit none

  private
  public :: mhd_nlsolverCallback
  public :: mhd_calcVelocityField
  public :: mhd_hadaptCallbackScalar2d
  public :: mhd_hadaptCallbackBlock2d

contains

  !*****************************************************************************

!<subroutine>

  subroutine mhd_nlsolverCallback(rproblemLevel, rtimestep, rsolver,&
                                  rsolution, rsolutionInitial,&
                                  rrhs, rres, istep, ioperationSpec,&
                                  rcollection, istatus)

!<description>
    ! This subroutine is called by the nonlinear solver and it is responsible
    ! to assemble preconditioner, right-hand side vector, residual vector, etc.
!</description>

!<input>
    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial
    
    ! number of solver step
    integer, intent(IN) :: istep
    
    ! specifier for operations
    integer(I32), intent(IN) :: ioperationSpec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
    
    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep
    
    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution
    
    ! right-hand side vector
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres
        
    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>

!<output>
    ! status flag
    integer, intent(OUT) :: istatus
!</output>
!</subroutine>
    
    ! local variable
    type(t_vectorBlock), pointer :: p_rsolutionTransport
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DdataTransport
    real(DP), dimension(:), pointer :: p_DdataEuler
    real(DP), dimension(:), pointer :: p_DdataResidual
    real(DP), dimension(:), pointer :: p_DdataMassMatrix
    integer :: neq, nvar, lumpedMassMatrix, isystemFormat

      
    ! Do we have to calculate the preconditioner?
    ! --------------------------------------------------------------------------
    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) .or.&
        (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0)) then
      
      call euler_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
                                    rsolution, rcollection)
    end if
    
    ! Do we have to calculate the residual and the constant right-hand side
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      
      call euler_calcResidual(rproblemLevel, rtimestep, rsolver,&
                              rsolution, rsolutionInitial,&
                              rrhs, rres, istep, rcollection)
      
      ! Get solution from scalar transport model
      p_rsolutionTransport => rcollection%p_rvectorQuickAccess1
      
      ! Set pointer to global solution vectors
      call lsysbl_getbase_double(rsolution, p_DdataEuler)
      call lsysbl_getbase_double(rres, p_DdataResidual)
      call lsysbl_getbase_double(p_rsolutionTransport, p_DdataTransport)
      
      ! Get lumped mass matrix
      lumpedMassMatrix = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
      call lsyssc_getbase_double(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                 p_DdataMassMatrix)
      
      ! Set pointer to the vertex coordinates
      call storage_getbase_double2D(&
          rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)
      
      ! Set dimensions
      neq  = p_rsolutionTransport%NEQ
      nvar = euler_getNVAR(rproblemLevel)
      
      isystemFormat = collct_getvalue_int(rcollection, 'isystemformat')
      select case(isystemFormat)
      case (SYSTEM_INTERLEAVEFORMAT)
        call calcSourceTermInterleaveFormat(rtimestep%dTime, rtimestep%dStep,&
                                            neq, nvar, p_DvertexCoords,&
                                            p_DdataMassMatrix, p_DdataTransport,&
                                            p_DdataEuler, p_DdataResidual)
      case (SYSTEM_BLOCKFORMAT)
        call calcSourceTermBlockFormat(rtimestep%dTime, rtimestep%dStep,&
                                       neq, nvar, p_DvertexCoords,&
                                       p_DdataMassMatrix, p_DdataTransport,&
                                       p_DdataEuler, p_DdataResidual)
      case DEFAULT
        call output_line('Invalid system format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'mhd_nlsolverCallback')
        call sys_halt()
      end select
    end if

    
    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      
      call euler_setBoundary(rproblemLevel, rtimestep, rsolver,&
                             rsolution, rsolutionInitial, rres, rcollection)
    end if

    ! Set status flag
    istatus = 0
    
  contains
    
    ! Here, the real working routines follow
    
    !**************************************************************
    
    subroutine calcSourceTermInterleaveFormat(dtime, dstep, neq, nvar, DvertexCoords,&
                                              DdataMassMatrix, DdataTransport,&
                                              DdataEuler, DdataResidual)

      real(DP), dimension(nvar,neq), intent(IN) :: DdataEuler
      real(DP), dimension(:,:), intent(IN) :: DvertexCoords
      real(DP), dimension(:), intent(IN) :: DdataMassMatrix
      real(DP), dimension(:), intent(IN) :: DdataTransport
      real(DP), intent(IN) :: dtime, dstep
      integer, intent(IN) :: neq, nvar
      
      real(DP), dimension(nvar,neq), intent(INOUT) :: DdataResidual
      
      ! local variables
      real(DP) :: dradius, daux, dscale, v1, v2, x1, x2
      integer :: ieq
      
      
      ! Compute the scaling parameter
      dscale = -dstep * 12.0 * (1.0-dtime**4) * dtime**2
!!$      dscale = -dstep * 12.0 * dtime*dtime
      
      ! Loop over all nodal values
      do ieq = 1, neq
        
        ! Get coodrinates
        x1 = DvertexCoords(1, ieq)
        x2 = DvertexCoords(2, ieq)
          
        ! Compute distance from origin
        dradius = sqrt(x1*x1 + x2*x2)
        
        ! Compute unit vector emanating from the origin
        if (x1> 0.0_DP) then
          daux = atan(x2/x1)
          x1   = cos(daux)
          x2   = sin(daux)
        elseif (x1 < 0.0_DP) then
          daux = atan(x2/x1)
          x1   = -cos(daux)
          x2   = -sin(daux)
        else
          if (x2 > 0.0) then
            x1 = 0.0_DP
            x2 = 1.0_DP
          else
            x1 =  0.0_DP
            x2 = -1.0_DP
          end if
        end if
        
        ! Compute auxiliary quantity
        daux = dscale * DdataMassMatrix(ieq) * DdataTransport(ieq) / max(dradius, 1.0e-4_DP)
        
        ! Compute velocities
        v1 = DdataEuler(2, ieq)/DdataEuler(1, ieq)
        v2 = DdataEuler(3, ieq)/DdataEuler(1, ieq)
        
        ! Impose source values into global vector
        DdataResidual(2, ieq) = DdataResidual(2, ieq) + daux * x1
        DdataResidual(3, ieq) = DdataResidual(3, ieq) + daux * x2
        DdataResidual(4, ieq) = DdataResidual(4, ieq) + max(daux * (x1*v1 + x2*v2), 0.0_DP)
      end do
      
    end subroutine calcSourceTermInterleaveFormat

    !**************************************************************
    
    subroutine calcSourceTermBlockFormat(dtime, dstep, neq, nvar, DvertexCoords,&
                                         DdataMassMatrix, DdataTransport,&
                                         DdataEuler, DdataResidual)
      
      real(DP), dimension(neq,nvar), intent(IN) :: DdataEuler
      real(DP), dimension(:,:), intent(IN) :: DvertexCoords
      real(DP), dimension(:), intent(IN) :: DdataMassMatrix
      real(DP), dimension(:), intent(IN) :: DdataTransport
      real(DP), intent(IN) :: dtime, dstep
      integer, intent(IN) :: neq, nvar
      
      real(DP), dimension(neq,nvar), intent(INOUT) :: DdataResidual
      
      ! local variables
      real(DP) :: dradius, daux, dscale, v1, v2, x1, x2
      integer :: ieq
      
      
      ! Compute the scaling parameter
      dscale = -dstep * 12.0 * (1.0-dtime**4) * dtime**2
!!$      dscale = -dstep * 12.0 * dtime*dtime
      
      ! Loop over all nodal values
      do ieq = 1, neq
        
        ! Get coodrinates
        x1 = DvertexCoords(1, ieq)
        x2 = DvertexCoords(2, ieq)
          
        ! Compute distance from origin
        dradius = sqrt(x1*x1 + x2*x2)
        
        ! Compute unit vector emanating from the origin
        if (x1> 0.0_DP) then
          daux = atan(x2/x1)
          x1   = cos(daux)
          x2   = sin(daux)
        elseif (x1 < 0.0_DP) then
          daux = atan(x2/x1)
          x1   = -cos(daux)
          x2   = -sin(daux)
        else
          if (x2 > 0.0) then
            x1 = 0.0_DP
            x2 = 1.0_DP
          else
            x1 =  0.0_DP
            x2 = -1.0_DP
          end if
        end if
        
        ! Compute auxiliary quantity
        daux = dscale * DdataMassMatrix(ieq) * DdataTransport(ieq) / max(dradius, 1.0e-4_DP)
        
        ! Compute velocities
        v1 = DdataEuler(ieq, 2)/DdataEuler(ieq, 1)
        v2 = DdataEuler(ieq, 3)/DdataEuler(ieq, 1)
        
        ! Impose source values into global vector
        DdataResidual(ieq, 2) = DdataResidual(ieq, 2) + daux * x1
        DdataResidual(ieq, 3) = DdataResidual(ieq, 3) + daux * x2
        DdataResidual(ieq, 4) = DdataResidual(ieq, 4) + max(daux * (x1*v1 + x2*v2), 0.0_DP)
      end do
      
    end subroutine calcSourceTermBlockFormat
    
  end subroutine mhd_nlsolverCallback
  
  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcVelocityField(rproblemLevel, rsolution, rcollection, nlminOpt)

!<description>
    ! This subroutine calculates the velocity field from the solution
    ! of the compressible Euler model. The result is stored separately
    ! for each problem level.
!</description>

!<input>
    ! solution vector of compressible Euler model
    type(t_vectorBlock), intent(IN) :: rsolution

    ! OPTIONAL: minimum problem level
    integer, intent(IN), optional :: nlminOpt
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT), target :: rproblemLevel

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: velocityfield
    integer :: neq, ndim
    

    ! Get parameter from collection
    velocityfield = collct_getvalue_int(rcollection, 'velocityfield')

    ! Get number of degrees of freedom and spatial dimension
    neq  = rproblemLevel%rtriangulation%NVT
    ndim = rproblemLevel%rtriangulation%ndim

    ! Create/resize velocity vector if required
    if (rproblemLevel%RvectorBlock(velocityfield)%NEQ .eq. 0) then
      call lsysbl_createVectorBlock(&
          rproblemLevel%rvectorBlock(velocityfield), neq, ndim, .true.)
    elseif (rproblemLevel%RvectorBlock(velocityfield)%NEQ .ne. neq*ndim) then
      call lsysbl_resizeVectorBlock(&
          rproblemLevel%rvectorBlock(velocityfield), neq, .true.)
    end if

    
    ! Set x-velocity
    call euler_getVariable(rsolution, 'velocity_x',&
        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(1))

    ! Set y-velocity
    call euler_getVariable(rsolution, 'velocity_y',&
        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(2))

    ! Set update notification in problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     PROBLEV_MSPEC_UPDATE)

  end subroutine mhd_calcVelocityField

  !*****************************************************************************

!<subroutine>

  subroutine mhd_hadaptCallbackScalar2d(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
    ! to be store in scalar interleave format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(IN) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolutionEuler, rsolutionTransport
    real(DP), dimension(:), pointer, save :: p_DsolutionEuler, p_DsolutionTransport
    integer :: ivar,neq


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the solution of the
      ! Euler model is stored in the second quick access string and
      ! the name of the solution of the scalar transport model ost
      ! storen in the third quick access string.

      ! Retrieve solution vectors from colletion and set pointer
       rsolutionEuler     => collct_getvalue_vec(rcollection,&
                                           trim(rcollection%SquickAccess(2)))
      rsolutionTransport => collct_getvalue_vec(rcollection,&
                                           trim(rcollection%SquickAccess(3)))
     
      ! Check if solution is stored in interleave format
      if (rsolutionEuler%nblocks .ne. 1) then
        call output_line('Vector is not in interleave format!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'mhd_hadaptCallbackScalar2d')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)
      
      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vectors
      nullify(rsolutionEuler, p_DsolutionEuler)
      nullify(rsolutionTransport, p_DsolutionTransport)
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)
      

    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector for the Euler model
      if (rsolutionEuler%NEQ .ne. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D*Ivertices(1), .false., .true.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if

      ! Resize solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .ne. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector for the Euler model
      if (rsolutionEuler%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D*Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if
      do ivar = 1, NVAR2D
        p_DsolutionEuler((Ivertices(1)-1)*NVAR2D+ivar) = &
            0.5_DP*(p_DsolutionEuler((Ivertices(2)-1)*NVAR2D+ivar)+&
                    p_DsolutionEuler((Ivertices(3)-1)*NVAR2D+ivar))
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(Ivertices(1)) =&
          0.5_DP*(p_DsolutionTransport(Ivertices(2))+&    
                  p_DsolutionTransport(Ivertices(3)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector for the Euler model
      if (rsolutionEuler%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D*Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if
      do ivar = 1, NVAR2D
        p_DsolutionEuler((Ivertices(1)-1)*NVAR2D+ivar) = &
            0.25_DP*(p_DsolutionEuler((Ivertices(2)-1)*NVAR2D+ivar)+&
                     p_DsolutionEuler((Ivertices(3)-1)*NVAR2D+ivar)+&
                     p_DsolutionEuler((Ivertices(4)-1)*NVAR2D+ivar)+&
                     p_DsolutionEuler((Ivertices(5)-1)*NVAR2D+ivar))
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(Ivertices(1)) =&
          0.25_DP*(p_DsolutionTransport(Ivertices(2))+&
                   p_DsolutionTransport(Ivertices(3))+&
                   p_DsolutionTransport(Ivertices(4))+&
                   p_DsolutionTransport(Ivertices(5)))    
      

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution for the Euler model
      if (Ivertices(2) .ne. 0) then
        do ivar = 1, NVAR2D
          p_DsolutionEuler((Ivertices(1)-1)*NVAR2D+ivar) = &
              p_DsolutionEuler((Ivertices(2)-1)*NVAR2D+ivar)
        end do
      else
        do ivar = 1, NVAR2D
          p_DsolutionEuler((Ivertices(1)-1)*NVAR2D+ivar) = 0.0_DP
        end do
      end if

      ! Remove vertex from solution for the scalar transport model
      if (Ivertices(2) .ne. 0) then
        p_DsolutionTransport(Ivertices(1)) = p_DsolutionTransport(Ivertices(2))
      else
        p_DsolutionTransport(Ivertices(1)) = 0.0_DP
      end if
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

    end select
    
  end subroutine mhd_hadaptCallbackScalar2d

  !*****************************************************************************

!<subroutine>

  subroutine mhd_hadaptCallbackBlock2d(rcollection, iOperation, Ivertices, Ielements)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D. The solution vector is assumed
    ! to be store in block format.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(IN) :: iOperation

    ! Array of vertices involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ivertices

    ! Array of elements involved in the adaptivity step
    integer, dimension(:), intent(IN) :: Ielements
!</input>

!<inputoutput>
    ! Collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolutionEuler, rsolutionTransport
    real(DP), dimension(:), pointer, save :: p_DsolutionEuler, p_DsolutionTransport
    integer :: ivar,neq


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! This subroutine assumes that the name of the solution of the
      ! Euler model is stored in the second quick access string and
      ! the name of the solution of the scalar transport model ost
      ! storen in the third quick access string.

      ! Retrieve solution vectors from colletion and set pointer
       rsolutionEuler     => collct_getvalue_vec(rcollection,&
                                           trim(rcollection%SquickAccess(2)))
      rsolutionTransport => collct_getvalue_vec(rcollection,&
                                           trim(rcollection%SquickAccess(3)))
      
      ! Check if solution is stored in interleave format
      if (rsolutionEuler%nblocks .ne. NVAR2D) then
        call output_line('Vector is not in block format!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'mhd_hadaptCallbackBlock2d')
        call sys_halt()
      end if

      ! Set pointers
      call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)
      
      
    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vectors
      nullify(rsolutionEuler, p_DsolutionEuler)
      nullify(rsolutionTransport, p_DsolutionTransport)
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)
      

    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector for the Euler model
      if (rsolutionEuler%NEQ .ne. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D*Ivertices(1), .false., .true.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if

      ! Resize solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .ne. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector for the Euler model
      if (rsolutionEuler%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D*Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if
      neq = rsolutionEuler%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_DsolutionEuler((ivar-1)*neq+Ivertices(1)) = &
            0.5_DP*(p_DsolutionEuler((ivar-1)*neq+Ivertices(2))+&
                    p_DsolutionEuler((ivar-1)*neq+Ivertices(3)) )
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(Ivertices(1)) =&
          0.5_DP*(p_DsolutionTransport(Ivertices(2))+&    
                  p_DsolutionTransport(Ivertices(3)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

      
    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector for the Euler model
      if (rsolutionEuler%NEQ .lt. NVAR2D*Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionEuler, NVAR2D*Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionEuler, p_DsolutionEuler)
      end if
      neq = rsolutionEuler%NEQ/NVAR2D
      do ivar = 1, NVAR2D
        p_DsolutionEuler((ivar-1)*neq+Ivertices(1)) =&
            0.25_DP*(p_DsolutionEuler((ivar-1)*neq+Ivertices(2))+&
                     p_DsolutionEuler((ivar-1)*neq+Ivertices(3))+&
                     p_DsolutionEuler((ivar-1)*neq+Ivertices(4))+&
                     p_DsolutionEuler((ivar-1)*neq+Ivertices(5)) )
      end do

      ! Insert vertex into solution vector for the scalar transport model
      if (rsolutionTransport%NEQ .lt. Ivertices(1)) then
        call lsysbl_resizeVectorBlock(rsolutionTransport, Ivertices(1), .false.)
        call lsysbl_getbase_double(rsolutionTransport, p_DsolutionTransport)
      end if
      p_DsolutionTransport(Ivertices(1)) =&
          0.25_DP*(p_DsolutionTransport(Ivertices(2))+&
                   p_DsolutionTransport(Ivertices(3))+&
                   p_DsolutionTransport(Ivertices(4))+&
                   p_DsolutionTransport(Ivertices(5)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution for the Euler model
      if (Ivertices(2) .ne. 0) then
        neq = rsolutionEuler%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_DsolutionEuler((ivar-1)*neq+Ivertices(1)) = &
              p_DsolutionEuler((ivar-1)*neq+Ivertices(2))
        end do
      else
        neq = rsolutionEuler%NEQ/NVAR2D
        do ivar = 1, NVAR2D
          p_DsolutionEuler((ivar-1)*neq+Ivertices(1)) = 0.0_DP
        end do
      end if

      ! Remove vertex from solution for the scalar transport model
      if (Ivertices(2) .ne. 0) then
        p_DsolutionTransport(Ivertices(1)) = p_DsolutionTransport(Ivertices(2))
      else
        p_DsolutionTransport(Ivertices(1)) = 0.0_DP
      end if
      
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)


    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(rcollection, iOperation, Ivertices, Ielements)

    end select
    
  end subroutine mhd_hadaptCallbackBlock2d

end module mhd_callback
