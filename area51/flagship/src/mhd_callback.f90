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
!# 3.) mhd_calcSourceTerm
!#     -> Calculates the source term
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
  public :: mhd_calcSourceTerm
  public :: mhd_hadaptCallbackScalar2d
  public :: mhd_hadaptCallbackBlock2d

contains

  !*****************************************************************************

!<subroutine>

  subroutine mhd_nlsolverCallback(rproblemLevel, rtimestep, rsolver,&
                                  rsolution, rsolutionInitial,&
                                  rrhs, rres, istep, ioperationSpec,&
                                  rcollection, istatus, rb)

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
    
    ! OPTIONAL: constant right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rb
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
    type(t_vectorBlock) :: rvector
    type(t_vectorBlock), pointer :: p_rsolutionTransport
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DdataTransport, p_DdataEuler, p_Ddata, p_MC, p_ML
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    integer :: neq, nvar, isystemFormat, lumpedMassMatrix, consistentMassMatrix

    print *, "NOT USED ANY MORE"
    stop
      
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
      
      ! In the zeroth iteration the explicit source term is applied to
      ! right-hand side vector and the residual vector
      if (istep .eq. 0) then

        ! Get solution from scalar transport model
        p_rsolutionTransport => rcollection%p_rvectorQuickAccess1
        call lsysbl_getbase_double(p_rsolutionTransport, p_DdataTransport)
        
        ! Set pointer to global solution vector
        call lsysbl_getbase_double(rsolution, p_DdataEuler)
        call lsysbl_createVectorBlock(rsolution, rvector, .true., ST_DOUBLE)
        call lsysbl_getbase_double(rvector, p_Ddata)
        
        ! Get mass matrices
        lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
        consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
        call lsyssc_getbase_double(rproblemLevel%Rmatrix(lumpedMassMatrix), p_ML)
        call lsyssc_getbase_double(rproblemLevel%Rmatrix(consistentMassMatrix), p_MC)

        call lsyssc_getbase_Kld(rproblemLevel%Rmatrix(consistentMassMatrix), p_Kld)
        call lsyssc_getbase_Kcol(rproblemLevel%Rmatrix(consistentMassMatrix), p_Kcol)

        ! Set pointer to the vertex coordinates
        call storage_getbase_double2D(&
            rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)
    
        ! Set dimensions
        neq  = p_rsolutionTransport%NEQ
        nvar = euler_getNVAR(rproblemLevel)
        
        isystemFormat = collct_getvalue_int(rcollection, 'isystemformat')
        select case(isystemFormat)
        case (SYSTEM_INTERLEAVEFORMAT)
          call calcSourceTermInterleaveFormat(rtimestep%dTime, rtimestep%dStep, neq, nvar,&
                                              p_DvertexCoords, p_Kld, p_Kcol, p_MC,&
                                              p_DdataTransport, p_DdataEuler, p_Ddata)

        case DEFAULT
          call output_line('Invalid system format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'mhd_nlsolverCallback')
          call sys_halt()
        end select

        ! Apply source term to constant right-hand side and defect vector
        call lsysbl_vectorLinearComb(rvector, rrhs, 1.0_DP, 1.0_DP)
        call lsysbl_vectorLinearComb(rvector, rres, 1.0_DP, 1.0_DP)
      end if
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
                                              Kld, Kcol, MC, DdataTransport, DdataEuler, Ddata)

      real(DP), dimension(:,:), intent(IN) :: DvertexCoords
      real(DP), dimension(:), intent(IN) :: MC,DdataTransport
      real(DP), dimension(nvar,neq), intent(IN) :: DdataEuler
      real(DP), intent(IN) :: dtime, dstep
      integer, dimension(:), intent(IN) :: Kld, Kcol
      integer, intent(IN) :: neq, nvar
      
      real(DP), dimension(nvar,neq), intent(OUT) :: Ddata
      
      ! local variables
      real(DP) :: drad, dang, daux, dscale, x1, x2, p, rq
      integer :: i,j,ij

      ! Compute the scaling parameter
      dscale = -dstep * 12.0 * (1.0-dtime**4) * dtime**2
!!$      dscale = -dstep * 12.0 * dtime*dtime

      Ddata = 0.0_DP

      ! Loop over all rows
      do i = 1, neq

        ! Loop over all columns
        do ij = Kld(i), Kld(i+1)-1

          ! Get columns number
          j = Kcol(ij)
        
          ! Get coodrinates at node j
          x1 = DvertexCoords(1, j)
          x2 = DvertexCoords(2, j)
          
          ! Compute polar coordinates
          drad = sqrt(x1*x1 + x2*x2)
          dang = atan2(x2, x1)
          
          ! Compute unit vector into origin
          x1 = cos(dang)
          x2 = sin(dang)
          
          ! Compute source term 
          daux = dscale * MC(ij) * max(DdataTransport(j)*DdataEuler(1, j), 0.0_DP) / max(drad, 1.0e-4_DP)
                    
          ! Impose source values
          Ddata(2, i) = Ddata(2, i) + daux * x1
          Ddata(3, i) = Ddata(3, i) + daux * x2
          Ddata(4, i) = Ddata(4, i) + daux * ( x1 * DdataEuler(2, j) + x2 * DdataEuler(3, j)) / DdataEuler(1, j)

        end do
      end do

    end subroutine calcSourceTermInterleaveFormat
    
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

  subroutine mhd_calcSourceTerm(rproblemLevel, rtimestep,&
                                rsolutionTransport, rsolutionEuler, rcollection)

!<description>
    ! Apply the source term and update the solution from the Euler system
!</description>

!<input>
    ! solution vector for Transport model
    type(t_vectorBlock), intent(IN) :: rsolutionTransport
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
    
    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep
        
    ! solution vector for Euler model
    type(t_vectorBlock), intent(INOUT) :: rsolutionEuler
    
    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>


    ! local variable
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DdataTransport, p_DdataEuler, p_MC, p_ML
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    integer :: neq, nvar, isystemFormat, consistentMassMatrix, lumpedMassMatrix
    
    ! Get lumped and consistent mass matrix
    lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
    consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(lumpedMassMatrix), p_ML)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(consistentMassMatrix), p_MC)

    call lsyssc_getbase_Kld(rproblemLevel%Rmatrix(consistentMassMatrix), p_Kld)
    call lsyssc_getbase_Kcol(rproblemLevel%Rmatrix(consistentMassMatrix), p_Kcol)
    
    ! Set pointer to global solution vectors
    call lsysbl_getbase_double(rsolutionEuler, p_DdataEuler)
    call lsysbl_getbase_double(rsolutionTransport, p_DdataTransport)
    
    ! Set pointer to the vertex coordinates
    call storage_getbase_double2D(&
        rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)
    
    ! Set dimensions
    neq  = rsolutionTransport%NEQ
    nvar = euler_getNVAR(rproblemLevel)
    
    isystemFormat = collct_getvalue_int(rcollection, 'isystemformat')
    select case(isystemFormat)
    case (SYSTEM_INTERLEAVEFORMAT)
      call calcSourceTermInterleaveFormat(rtimestep%dTime, rtimestep%dStep, neq, nvar,&
                                          p_DvertexCoords, p_Kld, p_Kcol, p_MC, p_ML,&
                                          p_DdataTransport, p_DdataEuler)
    case DEFAULT
      call output_line('Invalid system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcSourceTerm')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the real working routines follow
    
    !**************************************************************
    
    subroutine calcSourceTermInterleaveFormat(dtime, dstep, neq, nvar, DvertexCoords,&
                                              Kld, Kcol, MC, ML, DdataTransport, DdataEuler)

      real(DP), dimension(:,:), intent(IN) :: DvertexCoords
      real(DP), dimension(:), intent(IN) :: MC,ML,DdataTransport
      real(DP), intent(IN) :: dtime, dstep
      integer, dimension(:), intent(IN) :: Kld, Kcol
      integer, intent(IN) :: neq, nvar
      
      real(DP), dimension(nvar,neq), intent(INOUT) :: DdataEuler
      
      ! local variables
      real(DP), dimension(:,:), allocatable :: DsourceTerm
      real(DP) :: drad, dang, daux, dscale, x1, x2, p, rq
      integer :: i,j,ij


      ! Compute the scaling parameter
      dscale = -dstep * 12.0 * (1.0-dtime**4) * dtime**2
!!$      dscale = -dstep * 12.0 * dtime*dtime

      allocate(DsourceTerm(2,neq)); DsourceTerm = 0.0_DP

      ! Loop over all rows
      do i = 1, neq

        ! Loop over all columns
        do ij = Kld(i), Kld(i+1)-1

          ! Get columns number
          j = Kcol(ij)
        
          ! Get coodrinates at node j
          x1 = DvertexCoords(1, j)
          x2 = DvertexCoords(2, j)
          
          ! Compute polar coordinates
          drad = sqrt(x1*x1 + x2*x2)
          dang = atan2(x2, x1)
          
          ! Compute unit vector into origin
          x1 = cos(dang)
          x2 = sin(dang)
          
          ! Compute source term
          if (DdataTransport(j) > sqrt(SYS_EPSREAL)) then
            daux = dscale * MC(ij) * DdataTransport(j) / max(drad, 1.0e-4_DP)
          else
            daux = 0.0_DP
          end if
                    
          ! Impose source values
          DsourceTerm(1, i) = DsourceTerm(1, i) + daux * x1
          DsourceTerm(2, i) = DsourceTerm(2, i) + daux * x2

        end do
      end do

      do i = 1, neq

        ! Compute kinetic energy from momentum values without source term
        rq = 0.5 * ( DdataEuler(2, i)*DdataEuler(2, i) +&
                     DdataEuler(3, i)*DdataEuler(3, i) ) / DdataEuler(1, i)

        ! Compute pressure value
        p = DdataEuler(4, i) - rq

        ! Update momentum equations
        DdataEuler(2, i) = DdataEuler(2, i) + DsourceTerm(1, i)/ML(i)
        DdataEuler(3, i) = DdataEuler(3, i) + DsourceTerm(2, i)/ML(i)

        ! Compute kinetic energy from momentum values with source term
        rq = 0.5 * ( DdataEuler(2, i)*DdataEuler(2, i) +&
                     DdataEuler(3, i)*DdataEuler(3, i) ) / DdataEuler(1, i)

        ! Update total energy equation
        DdataEuler(4, i) = p + rq

      end do
      
      deallocate(DsourceTerm)

    end subroutine calcSourceTermInterleaveFormat

  end subroutine mhd_calcSourceTerm

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
