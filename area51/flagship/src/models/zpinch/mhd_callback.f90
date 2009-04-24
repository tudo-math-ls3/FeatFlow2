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
!# 4.) mhd_calcLinearizedFCT
!#     -> Calculates the linearized FCT correction
!#
!# 5.) mhd_hadaptCallbackScalar2d
!#     -> Performs application specific tasks in the adaptation
!#        algorithm in 2D, whereby the vector is stored in interleave format
!#
!# 6.) mhd_hadaptCallbackBlock2d
!#     -> Performs application specific tasks in the adaptation
!#        algorithm in 2D, whereby the vector is stored in block format
!#
!# </purpose>
!##############################################################################

module mhd_callback

  use boundaryfilter
  use collection
  use euler_basic
  use euler_callback
  use euler_callback2d
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
  use transport_callback2d

  implicit none

  private
  public :: mhd_calcVelocityField
  public :: mhd_calcSourceTerm
  public :: mhd_calcTracerIndicator
  public :: mhd_calcLinearizedFCT
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

    ! Set global solution vector as external vector for the transport model
    call transp_setVariable2d(rsolution%RvectorBlock(1), 3)

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
          if (drad .gt. 1e-4) then
            x1 = cos(dang)
            x2 = sin(dang)
          else
            x1 = 0.0; x2 = 0.0
          end if
          
          ! Compute source term
          if (DdataTransport(j) > SYS_EPSREAL) then
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

  subroutine mhd_calcTracerIndicator(rvector, rerror)

!<description>
    ! This subroutine computes the error indicator based on the
    ! magnitude of the tracer quantity
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(IN) :: rvector
!</input>
      
!<inputoutput>
    ! local error distribution
    type(t_vectorScalar), intent(INOUT) :: rerror
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata, p_Derror
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisActiveElement
    integer :: iel,jel,ive,ivt,iprotectLayer
    
    ! Set pointer to the underlying triangulation
    p_rtriangulation => rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! Create empty error
    call lsyssc_createVector(rerror, p_rtriangulation%NEL, .true.)
    
    ! Set pointers
    call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2D(p_rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    call lsysbl_getbase_double(rvector, p_Ddata)
    call lsyssc_getbase_double(rerror, p_Derror)
    
    ! Compute the mean value in each elements
    do iel = 1, p_rtriangulation%NEL
      do ive = 1, tria_getNVE(p_IverticesAtElement, iel)
        ivt = p_IverticesAtElement(ive, iel)
        p_Derror(iel) = p_Derror(iel) + p_Ddata(ivt)
      end do
    end do

    ! Scale by the number of vertices
    do iel = 1, p_rtriangulation%NEL
      p_Derror(iel) = p_Derror(iel) / tria_getNVE(p_IverticesAtElement, iel)
    end do

    ! Add protection layers
    allocate(p_BisActiveElement(p_rtriangulation%NEL))
    
    do iprotectLayer = 1, 10

      p_BisActiveElement = .false.

      do iel = 1, p_rtriangulation%NEL

        if (p_BisactiveElement(iel)) cycle
        if (p_Derror(iel) .le. 1e4) cycle

        do ive = 1, tria_getNVE(p_IverticesAtElement, iel)
          jel = p_IneighboursAtElement(ive, iel)
          if (jel .eq. 0) cycle
          if (p_BisactiveElement(jel)) then
            p_Derror(jel) = max(p_Derror(jel), p_Derror(iel))
          else
            if (p_Derror(jel) .lt. 1e4) then
              p_Derror(jel) = max(p_Derror(jel), p_Derror(iel))
              p_BisactiveElement(jel) = .true.
            end if
          end if
        end do
      end do
    end do

    deallocate(p_BisActiveElement)
    
  end subroutine mhd_calcTracerIndicator

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcLinearizedFCT(rbdrCondEuler, rbdrCondTransport, rproblemLevel,&
                                   rtimestep, rsolutionEuler, rsolutionTransport, rcollection)

!<description>
    ! This subroutine calculates the linearized FCT correction for the
    ! Euler model and the scalar transport model simultaneously
!</description>

!<input>
    ! boundary condition structures
    type(t_boundaryCondition), intent(IN) :: rbdrCondEuler, rbdrCondTransport

    ! problem level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep
!</input>

!<inputoutput>
    ! solution vectors
    type(t_vectorBlock), intent(INOUT) :: rsolutionEuler, rsolutionTransport

    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_vectorScalar) :: rfluxEuler0, rfluxEuler, rfluxTransport0, rfluxTransport, ralpha
    type(t_vectorBlock) :: rdataEuler, rdataTransport
    real(DP), dimension(:), pointer :: p_MC, p_ML, p_Cx, p_Cy, p_solEuler, p_solTransport
    real(DP), dimension(:), pointer :: p_fluxEuler0, p_fluxEuler, p_fluxTransport0, p_fluxTransport, p_alpha
    real(DP), dimension(:), pointer :: p_dataEuler, p_dataTransport
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal, p_Ksep
    integer :: h_Ksep, templatematrix, lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, nedge

    templateMatrix       = collct_getvalue_int(rcollection, 'templatematrix')
    consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentMassMatrix')
    lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedMassMatrix')
    coeffMatrix_CX       = collct_getvalue_int(rcollection, 'coeffMatrix_CX')
    coeffMatrix_CY       = collct_getvalue_int(rcollection, 'coeffMatrix_CY')
    
    p_rmatrix => rproblemLevel%Rmatrix(templatematrix)

    ! Set pointers
    call lsyssc_getbase_Kld(p_rmatrix, p_Kld)
    call lsyssc_getbase_Kcol(p_rmatrix, p_Kcol)
    call lsyssc_getbase_Kdiagonal(p_rmatrix, p_Kdiagonal)
    
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(consistentMassMatrix), p_MC)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(lumpedMassMatrix), p_ML)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(coeffMatrix_CX), p_Cx)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(coeffMatrix_CY), p_Cy)

    ! Create diagonal separator
    h_Ksep = ST_NOHANDLE
    call storage_copy(p_rmatrix%h_Kld, h_Ksep)
    call storage_getbase_int(h_Ksep, p_Ksep, p_rmatrix%NEQ+1)

    ! Compute number of edges
    nedge = int(0.5*(p_rmatrix%NA-p_rmatrix%NEQ))

    ! Create auxiliary vectors
    call lsyssc_createVector(rfluxEuler0, nedge, NVAR2D, .true., ST_DOUBLE)
    call lsyssc_createVector(rfluxEuler,  nedge, NVAR2D, .true., ST_DOUBLE)
    call lsyssc_createVector(rfluxTransport0, nedge, .true., ST_DOUBLE)
    call lsyssc_createVector(rfluxTransport,  nedge, .true., ST_DOUBLE)
    call lsyssc_createVector(ralpha, nedge, .false., ST_DOUBLE)
    call lsysbl_createVectorBlock(rsolutionEuler, rdataEuler, .true.)
    call lsysbl_createVectorBlock(rsolutionTransport, rdataTransport, .true.)
    
    ! Set pointers
    call lsysbl_getbase_double(rsolutionEuler, p_solEuler)
    call lsysbl_getbase_double(rsolutionTransport, p_solTransport)
    call lsysbl_getbase_double(rdataEuler, p_dataEuler)
    call lsysbl_getbase_double(rdataTransport, p_dataTransport)
    call lsyssc_getbase_double(rfluxEuler, p_fluxEuler)
    call lsyssc_getbase_double(rfluxTransport, p_fluxTransport)
    call lsyssc_getbase_double(rfluxEuler0, p_fluxEuler0)
    call lsyssc_getbase_double(rfluxTransport0, p_fluxTransport0)
    call lsyssc_getbase_double(ralpha, p_alpha)

    ! >>> SYNCHRONIZED IMPLEMENTATION <<<

    ! Initialize alpha with ones
    p_alpha = 1.0_DP
      
    ! Build the fluxes
    call buildFluxCons2d(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                         NVAR2D, nedge, p_solEuler, rtimestep%dStep,&
                         p_MC, p_ML, p_Cx, p_Cy, p_DataEuler, p_fluxEuler, p_fluxEuler0)

    call buildFlux2d(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                     nedge, p_solTransport, rtimestep%dStep,&
                     p_MC, p_ML, p_Cx, p_Cy, p_DataTransport, p_fluxTransport, p_fluxTransport0)

    ! Build the correction factors
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                             NVAR2D, nedge, p_ML, p_fluxEuler, p_fluxEuler0, 1, p_alpha, p_solEuler)
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                             NVAR2D, nedge, p_ML, p_fluxEuler, p_fluxEuler0, 4, p_alpha, p_solEuler)
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                             1, nedge, p_ML, p_fluxTransport, p_fluxTransport0, 1, p_alpha, p_solTransport)

    ! Apply correction to low-order solutions
    call applyCorrection(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                         NVAR2D, nedge, p_ML, p_fluxEuler, p_alpha, p_DataEuler, p_solEuler)
    call applyCorrection(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                         1, nedge, p_ML, p_fluxTransport, p_alpha, p_DataTransport, p_solTransport)

    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCondEuler, rproblemLevel%rtriangulation,&
                                   rsolutionEuler, rtimestep%dTime,&
                                   rproblemLevel%p_rproblem%rboundary,&
                                   euler_calcBoundaryvalues2d)

    call bdrf_filterVectorExplicit(rbdrCondTransport, rproblemLevel%rtriangulation,&
                                   rsolutionTransport, rtimestep%dTime,&
                                   rproblemLevel%p_rproblem%rboundary)
    
    ! Release flux vectors
    call storage_free(h_Ksep)
    call lsyssc_releaseVector(rfluxEuler0)
    call lsyssc_releaseVector(rfluxTransport0)
    call lsyssc_releaseVector(rfluxEuler)
    call lsyssc_releaseVector(rfluxTransport)
    call lsysbl_releaseVector(rdataEuler)
    call lsysbl_releaseVector(rdataTransport)
    call lsyssc_releaseVector(ralpha)

  contains
    
    !***************************************************************************

    subroutine buildFluxCons2d(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, NEDGE, u,&
                               dscale, MC, ML, Cx, Cy, troc, flux, flux0)

      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: MC,ML,Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR,NEDGE
      
      integer, dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(NVAR,NEDGE), intent(INOUT) :: flux0,flux
      
      real(DP), dimension(NVAR,NEQ), intent(OUT) :: troc     
      
      ! local variables
      real(DP), dimension(NVAR) :: K_ij,K_ji,D_ij,Diff,F_ij,F_ji
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      integer :: ij,ji,i,j,iedge
      
      ! Initialize time rate of change
      call lalg_clearVector(troc)
      
      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1; iedge = iedge+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Calculate low-order flux
          call euler_calcFluxRusanov2d(u(:,i), u(:,j), C_ij, C_ji, i, j, dscale, F_ij, F_ji)
          
          ! Update the time rate of change vector
          troc(:,i) = troc(:,i) + F_ij
          troc(:,j) = troc(:,j) + F_ji

          ! Calculate diffusion coefficient
          call euler_calcMatrixRusanovDiag2d(u(:,i), u(:,j), C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)
          
          ! Compute solution difference
          Diff = u(:,j)-u(:,i)

          ! Compute the raw antidiffusive flux
          flux0(:,iedge) = -D_ij*Diff
          
        end do
      end do

      ! Scale the time rate of change by the lumped mass matrix
      do i = 1, NEQ
        troc(:,i) = troc(:,i)/ML(i)
      end do

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
          
          ! Apply mass antidiffusion
          flux(:,iedge) = flux0(:,iedge) + MC(ij)*(troc(:,i)-troc(:,j))
          
          ! Update edge counter
          iedge = iedge-1
          
        end do
      end do

    end subroutine buildFluxCons2d

    !***************************************************************************

    subroutine buildFlux2d(Kld, Kcol, Kdiagonal, Ksep, NEQ, NEDGE, u,&
                           dscale, MC, ML, Cx, Cy, troc, flux, flux0)

      real(DP), dimension(:), intent(IN) :: MC,ML,Cx,Cy,u
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE
      
      integer, dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT) :: flux0,flux
      
      real(DP), dimension(:), intent(OUT) :: troc     
      
      ! local variables
      real(DP) :: k_ii,k_ij,k_ji,d_ij,diff,f_ij,f_ji
      real(DP), dimension(NDIM2D) :: C_ii,C_ij,C_ji
      integer :: ii,ij,ji,i,j,iedge
      
      ! Initialize time rate of change
      call lalg_clearVector(troc)
      
      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficient
        C_ii(1) = Cx(ii);   C_ii(2) = Cy(ii)

        ! Compute convection coefficients
        call transp_calcMatrixPrimalConst2d(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii, d_ij)

        ! Update the time rate of change vector
        troc(i) = troc(i) + dscale*k_ii*u(i)

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1; iedge = iedge+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Compute convection coefficients
          call transp_calcMatrixPrimalConst2d(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

          ! Compute solution difference
          diff = u(j)-u(i)

          ! Compute fluxes
          f_ij = dscale * (k_ij*u(j) + d_ij*diff)
          f_ji = dscale * (k_ji*u(i) - d_ij*diff)

          ! Update the time rate of change vector
          troc(i) = troc(i) + f_ij
          troc(j) = troc(j) + f_ji

          ! Compute the raw antidiffusive flux
          flux0(iedge) = -d_ij*diff
          
        end do
      end do

      ! Scale the time rate of change by the lumped mass matrix
      do i = 1, NEQ
        troc(i) = troc(i)/ML(i)
      end do

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
          
          ! Apply mass antidiffusion
          flux(iedge) = flux0(iedge) + MC(ij)*(troc(i)-troc(j))
          
          ! Update edge counter
          iedge = iedge-1
          
        end do
      end do

    end subroutine buildFlux2d
    
    !***************************************************************************
    
    subroutine  buildCorrectionCons(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, NEDGE,&
                                    ML, flux, flux0, ivar, alpha, u)

      real(DP), dimension(NVAR,NEDGE), intent(IN) :: flux0
      real(DP), dimension(:), intent(IN) :: ML
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR,NEDGE,ivar
      
      real(DP), dimension(NVAR,NEDGE), intent(INOUT) :: flux
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: u
      real(DP), dimension(:), intent(INOUT) :: alpha
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(:), allocatable :: pp,pm,qp,qm,rp,rm
      real(DP), dimension(NVAR) :: Cdiff
      real(DP) :: f_ij,f0_ij,diff,diff_ij,diff_ji,aux,p_ij,p_ji,p0_ij,p0_ji,u_i,u_j,v_i,v_j,r_i,r_j
      integer :: ij,ji,i,j,iedge

      allocate(pp(neq), pm(neq), qp(neq), qm(neq), rp(neq), rm(neq))
      
      pp = 0.0_DP; pm = 0.0_DP
      qp = 0.0_DP; qm = 0.0_DP
      rp = 1.0_DP; rm = 1.0_DP

      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; iedge = iedge+1
          
          ! Flux correction in conservative variables
          f_ij  = flux(ivar,iedge)
          f0_ij = flux0(ivar,iedge)
          diff  = u(ivar,j)-u(ivar,i)
          
          ! MinMod prelimiting of antidiffusive fluxes
          if (f_ij > SYS_EPSREAL .and. f0_ij > SYS_EPSREAL) then
            aux = min(f_ij, 2*f0_ij)
            alpha(iedge) = min(alpha(iedge), aux/f_ij)
            f_ij = aux
          elseif (f_ij < - SYS_EPSREAL .and. f0_ij < -SYS_EPSREAL) then
            aux = max(f_ij, 2*f0_ij)
            alpha(iedge) = min(alpha(iedge), aux/f_ij)
            f_ij = aux
          else
            f_ij = 0.0; alpha(iedge) = 0.0
          end if
          
          ! Sums of raw antidiffusive fluxes
          pp(i) = pp(i) + max(0.0_DP,  f_ij)
          pp(j) = pp(j) + max(0.0_DP, -f_ij)
          pm(i) = pm(i) + min(0.0_DP,  f_ij)
          pm(j) = pm(j) + min(0.0_DP, -f_ij)

          ! Sums of admissible edge contributions
          qp(i) = max(qp(i),  diff)
          qp(j) = max(qp(j), -diff)
          qm(i) = min(qm(i),  diff)
          qm(j) = min(qm(j), -diff)
          
!!$          ! Flux correction in primitive variables
!!$          select case(ivar)
!!$          case default
!!$            
!!$            ! Conservative variables
!!$            p_ij =  flux(ivar,iedge)
!!$            p_ji = -flux(ivar,iedge)
!!$
!!$            p0_ij =  flux0(ivar,iedge)
!!$            p0_ji = -flux0(ivar,iedge)
!!$
!!$            diff_ij = u(ivar,j)-u(ivar,i)
!!$            diff_ji = u(ivar,i)-u(ivar,j)
!!$
!!$            ! MinMod prelimiting of antidiffusive fluxes
!!$            if ((p_ij >  SYS_EPSREAL .and. p0_ij >  SYS_EPSREAL) .or.&
!!$                (p_ij < -SYS_EPSREAL .and. p0_ij < -SYS_EPSREAL)) then
!!$              aux = min(p_ij, 2*p0_ij)
!!$              alpha(iedge) = min(alpha(iedge), aux/p_ij)
!!$              p_ij = aux
!!$            else
!!$              p_ij = 0.0; alpha(iedge) = 0.0
!!$            end if
!!$
!!$            if ((p_ji >  SYS_EPSREAL .and. p0_ji >  SYS_EPSREAL) .or.&
!!$                (p_ji < -SYS_EPSREAL .and. p0_ji < -SYS_EPSREAL)) then
!!$              aux = min(p_ji, 2*p0_ji)
!!$              alpha(iedge) = min(alpha(iedge), aux/p_ji)
!!$              p_ji = aux
!!$            else
!!$              p_ji = 0.0; alpha(iedge) = 0.0
!!$            end if
!!$            
!!$            ! Sums of raw antidiffusive fluxes
!!$            pp(i) = pp(i) + max(0.0_DP, p_ij)
!!$            pp(j) = pp(j) + max(0.0_DP, p_ji)
!!$            pm(i) = pm(i) + min(0.0_DP, p_ij)
!!$            pm(j) = pm(j) + min(0.0_DP, p_ji)
!!$            
!!$            ! Sums of admissible edge contributions
!!$            qp(i) = max(qp(i), diff_ij)
!!$            qp(j) = max(qp(j), diff_ji)
!!$            qm(i) = min(qm(i), diff_ij)
!!$            qm(j) = min(qm(j), diff_ji)
!!$
!!$          case default
!!$            
!!$            ! Conservative variables
!!$            f_ij  = flux(ivar,iedge)
!!$            f0_ij = flux0(ivar,iedge)
!!$            diff  = u(ivar,j)-u(ivar,i)
!!$
!!$            ! MinMod prelimiting of antidiffusive fluxes
!!$            if (f_ij > SYS_EPSREAL .and. f0_ij > SYS_EPSREAL) then
!!$              aux = min(f_ij, 2*f0_ij)
!!$              alpha(iedge) = min(alpha(iedge), aux/f_ij)
!!$              f_ij = aux
!!$            elseif (f_ij < - SYS_EPSREAL .and. f0_ij < -SYS_EPSREAL) then
!!$              aux = max(f_ij, 2*f0_ij)
!!$              alpha(iedge) = min(alpha(iedge), aux/f_ij)
!!$              f_ij = aux
!!$            else
!!$              f_ij = 0.0; alpha(iedge) = 0.0
!!$            end if
!!$            
!!$            ! Sums of raw antidiffusive fluxes
!!$            pp(i) = pp(i) + max(0.0_DP,  f_ij)
!!$            pp(j) = pp(j) + max(0.0_DP, -f_ij)
!!$            pm(i) = pm(i) + min(0.0_DP,  f_ij)
!!$            pm(j) = pm(j) + min(0.0_DP, -f_ij)
!!$            
!!$            ! Sums of admissible edge contributions
!!$            qp(i) = max(qp(i),  diff)
!!$            qp(j) = max(qp(j), -diff)
!!$            qm(i) = min(qm(i),  diff)
!!$            qm(j) = min(qm(j), -diff)
!!$            
!!$          case (4)
!!$            
!!$            ! Velocities
!!$            u_i = u(2,i)/u(1,i);   v_i = u(3,i)/u(1,i)
!!$            u_j = u(2,j)/u(1,j);   v_j = u(3,j)/u(1,j)
!!$
!!$            ! Pressure variables
!!$            p_ij = (GAMMA-1) * (flux(4, iedge) + &
!!$                         0.5 * (u_i*u_i + v_i*v_i)*flux(1, iedge) -&
!!$                                u_i*flux(2, iedge) - v_i*flux(3, iedge) )
!!$            p_ji = -(GAMMA-1) * (flux(4, iedge) + &
!!$                          0.5 * (u_j*u_j + v_j*v_j)*flux(1, iedge) -&
!!$                                 u_j*flux(2, iedge) - v_j*flux(3, iedge) )
!!$            p0_ij = (GAMMA-1) * (flux0(4, iedge) + &
!!$                          0.5 * (u_i*u_i + v_i*v_i)*flux0(1, iedge) -&
!!$                                 u_i*flux0(2, iedge) - v_i*flux0(3, iedge) )
!!$            p0_ji = -(GAMMA-1) * (flux0(4, iedge) + &
!!$                           0.5 * (u_j*u_j + v_j*v_j)*flux0(1, iedge) -&
!!$                                  u_j*flux0(2, iedge) - v_j*flux0(3, iedge) )
!!$
!!$            Cdiff = u(:,j)-u(:,i)
!!#
!!$            diff_ij = (GAMMA-1) * (Cdiff(4) + &
!!$                            0.5 * (u_i*u_i + v_i*v_i)*Cdiff(1) -&
!!$                                   u_i*Cdiff(2) - v_i*Cdiff(3) )
!!$
!!$            diff_ji = -(GAMMA-1) * (Cdiff(4) + &
!!$                             0.5 * (u_j*u_j + v_j*v_j)*Cdiff(1) -&
!!$                                    u_j*Cdiff(2) - v_j*Cdiff(3) )
!!$
!!$            ! Pressure variables
!!$            p_ij =   (GAMMA-1) * ( u(4,i)*flux(1,iedge) - u(2,i)*flux(2,iedge) &
!!$                                  -u(3,i)*flux(3,iedge) + u(1,i)*flux(4,iedge) )
!!$            p_ji =  -(GAMMA-1) * ( u(4,j)*flux(1,iedge) - u(2,j)*flux(2,iedge) &
!!$                                  -u(3,j)*flux(3,iedge) + u(1,j)*flux(4,iedge) )
!!$            p0_ij =  (GAMMA-1) * ( u(4,i)*flux0(1,iedge) - u(2,i)*flux0(2,iedge) &
!!$                                  -u(3,i)*flux0(3,iedge) + u(1,i)*flux0(4,iedge) )
!!$            p0_ji = -(GAMMA-1) * ( u(4,j)*flux0(1,iedge) - u(2,j)*flux0(2,iedge) &
!!$                                  -u(3,j)*flux0(3,iedge) + u(1,j)*flux0(4,iedge) )
!!$
!!$            ! Solution differences
!!$            Cdiff = u(:,j)-u(:,i)
!!$
!!$            diff_ij =  (GAMMA-1) * ( u(4,i)*Cdiff(1) - u(2,i)*Cdiff(2) &
!!$                                    -u(3,i)*Cdiff(3) + u(1,i)*Cdiff(4) )
!!$            diff_ji = -(GAMMA-1) * ( u(4,j)*Cdiff(1) - u(2,j)*Cdiff(2) &
!!$                                    -u(3,j)*Cdiff(3) + u(1,j)*Cdiff(4) )
!!$            
!!$            ! MinMod prelimiting of antidiffusive fluxes
!!$            if ((p_ij >  SYS_EPSREAL .and. p0_ij >  SYS_EPSREAL) .or.&
!!$                (p_ij < -SYS_EPSREAL .and. p0_ij < -SYS_EPSREAL)) then
!!$              aux = min(p_ij, 2*p0_ij)
!!$              alpha(iedge) = min(alpha(iedge), aux/p_ij)
!!$              p_ij = aux
!!$            else
!!$              p_ij = 0.0; alpha(iedge) = 0.0
!!$            end if
!!$
!!$            if ((p_ji >  SYS_EPSREAL .and. p0_ji >  SYS_EPSREAL) .or.&
!!$                (p_ji < -SYS_EPSREAL .and. p0_ji < -SYS_EPSREAL)) then
!!$              aux = min(p_ji, 2*p0_ji)
!!$              alpha(iedge) = min(alpha(iedge), aux/p_ji)
!!$              p_ji = aux
!!$            else
!!$              p_ji = 0.0; alpha(iedge) = 0.0
!!$            end if
!!$            
!!$            ! Sums of raw antidiffusive fluxes
!!$            pp(i) = pp(i) + max(0.0_DP, p_ij)
!!$            pp(j) = pp(j) + max(0.0_DP, p_ji)
!!$            pm(i) = pm(i) + min(0.0_DP, p_ij)
!!$            pm(j) = pm(j) + min(0.0_DP, p_ji)
!!$            
!!$            ! Sums of admissible edge contributions
!!$            qp(i) = max(qp(i), diff_ij)
!!$            qp(j) = max(qp(j), diff_ji)
!!$            qm(i) = min(qm(i), diff_ij)
!!$            qm(j) = min(qm(j), diff_ji)
!!$          end select
            
        end do
      end do

      ! Compute nodal correction factors
      do i = 1, NEQ
        qp(i) = qp(i)*ML(i)
        qm(i) = qm(i)*ML(i)

        if (pp(i) > qp(i) + SYS_EPSREAL) rp(i) = qp(i)/pp(i)
        if (pm(i) < qm(i) - SYS_EPSREAL) rm(i) = qm(i)/pm(i)
      end do
      
      ! Loop over all rows (backward)
      do i = NEQ, 1, -1
        
        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Flux correction in conservative variables
          f_ij = flux(ivar,iedge)

          ! Limit conservative fluxes
          if (f_ij > SYS_EPSREAL) then
            alpha(iedge) = min(alpha(iedge), rp(i), rm(j))
          elseif (f_ij < -SYS_EPSREAL) then
            alpha(iedge) = min(alpha(iedge), rm(i), rp(j))
          end if

!!$          ! Flux correction in primitive variables
!!$          select case(ivar)
!!$          case default
!!$            
!!$            ! Flux correction in conservative variables
!!$            f_ij = flux(ivar,iedge)
!!$            
!!$            ! Limit conservative fluxes
!!$            if (f_ij > SYS_EPSREAL) then
!!$              alpha(iedge) = min(alpha(iedge), rp(i), rm(j))
!!$            elseif (f_ij < -SYS_EPSREAL) then
!!$              alpha(iedge) = min(alpha(iedge), rm(i), rp(j))
!!$            end if
!!$            
!!$            p_ij =  flux(ivar,iedge)
!!$            p_ji = -flux(ivar,iedge)
!!$
!!$          case (4)
!!$            
!!$            ! Velocities
!!$            u_i = u(2,i)/u(1,i);   v_i = u(3,i)/u(1,i)
!!$            u_j = u(2,j)/u(1,j);   v_j = u(3,j)/u(1,j)
!!$
!!$            ! Pressure variables
!!$            p_ij = (GAMMA-1) * (flux(4, iedge) + &
!!$                         0.5 * (u_i*u_i + v_i*v_i)*flux(1, iedge) -&
!!$                                u_i*flux(2, iedge) - v_i*flux(3, iedge) )
!!$            p_ji = -(GAMMA-1) * (flux(4, iedge) + &
!!$                          0.5 * (u_j*u_j + v_j*v_j)*flux(1, iedge) -&
!!$                                 u_j*flux(2, iedge) - v_j*flux(3, iedge) )
!!$
!!$            ! Pressure variables
!!$            p_ij = (GAMMA-1) * ( u(4,i)*flux(1,iedge) - u(2,i)*flux(2,iedge) &
!!$                                -u(3,i)*flux(3,iedge) + u(1,i)*flux(4,iedge) )
!!$            
!!$            p_ji = -(GAMMA-1) * ( u(4,j)*flux(1,iedge) - u(2,j)*flux(2,iedge) &
!!$                                 -u(3,j)*flux(3,iedge) + u(1,j)*flux(4,iedge) )
!!$
!!$            if (p_ij > SYS_EPSREAL) then
!!$              r_i = rp(i)
!!$            elseif (p_ij < -SYS_EPSREAL) then
!!$              r_i = rm(i)
!!$            else
!!$              r_i = 1.0
!!$            end if
!!$
!!$            if (p_ji > SYS_EPSREAL) then
!!$              r_j = rp(j)
!!$            elseif (p_ji < -SYS_EPSREAL) then
!!$              r_j = rm(j)
!!$            else
!!$              r_j = 1.0
!!$            end if
!!$
!!$            alpha(iedge) = min(alpha(iedge), r_i, r_j)
!!$
!!$          end select
          
          ! Update edge counter
          iedge = iedge-1
          
        end do
      end do

    end subroutine buildCorrectionCons

    !***************************************************************************

    subroutine applyCorrection(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
                               NEDGE, ML, flux, alpha, data, u)
      
      real(DP), dimension(NVAR,NEDGE), intent(IN) :: flux
      real(DP), dimension(:), intent(IN) :: ML,alpha
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR,NEDGE
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: u,data
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(NVAR) :: f_ij
      integer :: ij,ji,i,j,iedge

      ! Initialize correction
      call lalg_clearVector(data)

      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; iedge = iedge+1

          ! Limit raw antidiffusive flux
          f_ij = alpha(iedge)*flux(:,iedge)
          
          ! Apply correction
          data(:,i) = data(:,i) + f_ij
          data(:,j) = data(:,j) - f_ij
        end do
      end do

      ! Loop over all rows
      do i = 1, NEQ
        u(:,i) = u(:,i) + data(:,i)/ML(i)
      end do

      ! Just to be sure that this routine can be called repeatedly
      Ksep = Kld
      
    end subroutine applyCorrection

    !***************************************************************************

    pure elemental function minmod(a,b)
      real(DP), intent(IN) :: a,b
      real(DP) :: minmod

      if (a > 0 .and. b > 0) then
        minmod = min(a,b)
      elseif (a < 0 .and. b < 0) then
        minmod = max(a,b)
      else
        minmod = 0
      end if
    end function minmod

  end subroutine mhd_calcLinearizedFCT

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
                         OU_CLASS_ERROR,OU_MODE_STD,'mhd_hadaptCallbackScalar2d')
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
                         OU_CLASS_ERROR,OU_MODE_STD,'mhd_hadaptCallbackBlock2d')
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
