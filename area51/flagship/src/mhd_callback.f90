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
!# </purpose>
!##############################################################################

module mhd_callback

  use collection
  use euler_basic
  use euler_callback
  use flagship_basic
  use fsystem
  use genoutput
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
    integer :: neq, nvar, lumpedMassMatrix

      
    ! Do we have to calculate the preconditioner?
    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) .or.&
        (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0)) then
      
      call euler_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
                                    rsolution, rcollection)
    end if
    
    ! Do we have to calculate the residual and the constant right-hand side
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
      
      call calcSourceTermInterleaveFormat(rtimestep%dTime, rtimestep%dStep,&
                                          neq, nvar, p_DvertexCoords,&
                                          p_DdataMassMatrix, p_DdataTransport,&
                                          p_DdataEuler, p_DdataResidual)
    end if
    
    ! Do we have to impose boundary conditions?
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
!!$      dscale = -dstep * 12.0 * (1.0-dtime**4) * dtime**2
      dscale = -dstep * 12.0 * dtime*dtime
      
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

end module mhd_callback
