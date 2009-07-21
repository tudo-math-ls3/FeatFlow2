!##############################################################################
!# ****************************************************************************
!# <name> zpinch_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the simplified MHD equations in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) zpich_nlsolverCallbackTransport
!#     -> Callback routine for the nonlinear solver used in the
!#        computation of the scalar transport model
!#
!# 2.) zpinch_initVelocityField
!#     -> Initializes the velocity field for the transport model
!#
!# 3.) zpinch_initDensityAveraging
!#     -> Initializes the density averaged mass matrices
!#
!# 4.) zpinch_calcLinearizedFCT
!#     -> Calculates the linearized FCT correction
!#
!# </purpose>
!##############################################################################

module zpinch_callback

  use boundaryfilter
  use collection
  use derivatives
  use euler_basic
  use euler_callback2d
  use fsystem
  use basicgeometry
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use scalarpde
  use solveraux
  use storage
  use timestepaux
  use transport_callback
  use transport_callback2d
  use trilinearformevaluation

  implicit none

  private
  public :: zpinch_nlsolverCallbackTransport
  public :: zpinch_initVelocityField
  public :: zpinch_initDensityAveraging
  public :: zpinch_calcLinearizedFCT
  
contains

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_nlsolverCallbackTransport(rproblemLevel,&
      rtimestep, rsolver, rsolution, rsolutionInitial, rrhs, rres,&
      istep, ioperationSpec, rcollection, istatus, rb)

!<description>
    ! This subroutine is called by the nonlinear solver and it is responsible
    ! to assemble preconditioner, right-hand side vector, residual vector, etc.
!</description>

!<input>
    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolutionInitial
    
    ! number of solver step
    integer, intent(in) :: istep
    
    ! specifier for operations
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorBlock), intent(in), optional :: rb
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel
    
    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep
    
    ! solver structure
    type(t_solver), intent(inout) :: rsolver
    
    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution
        
    ! right-hand side vector
    type(t_vectorBlock), intent(inout) :: rrhs

    ! residual vector
    type(t_vectorBlock), intent(inout) :: rres

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! status flag
    integer, intent(out) :: istatus
!</output>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_vectorBlock), pointer :: p_rsolutionEuler
    integer(i32) :: iSpec
    integer :: jacobianMatrix

    !###########################################################################
    ! REMARK: The order in which the operations are performed is
    ! essential. This is due to the fact that the calculation of the
    ! residual/rhs requires the discrete transport operator to be
    ! initialized which is assembled in the calculation of the
    ! preconditioner. To prevent the re-assembly of the
    ! preconditioner twice, we remove the specifier
    ! NLSOL_OPSPEC_CALCPRECOND if the residual/rhs vector is built.
    !###########################################################################

    ! Make a local copy
    iSpec = ioperationSpec
    
    ! Do we have to calculate the constant right-hand side?
    ! --------------------------------------------------------------------------
    if ((iand(iSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then

      ! Compute the preconditioner
      call transp_calcPreconditioner(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)
      
      ! Compute the right-hand side
      call transp_calcRHS(rproblemLevel, rtimestep, rsolver,&
          rsolution, rsolutionInitial, rrhs, istep, rcollection, rb)
      
      ! Remove specifier for the preconditioner (if any)
      iSpec = iand(iSpec, not(NLSOL_OPSPEC_CALCPRECOND))
    end if
    
    
    ! Do we have to calculate the residual?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      if (istep .eq. 0) then
        ! Assemble the constant right-hand side
        call transp_calcExplicitRHS(rproblemLevel, rtimestep, rsolver,&
            rsolution, rrhs, rcollection, rb)

        ! Set pointer to parameter list
        p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
        
        ! Set pointer to solution vector
        p_rsolutionEuler => rcollection%p_rvectorQuickAccess1
        
        ! Calculate density averaged mass matrices
        call zpinch_initDensityAveraging(p_rparlist,&
            rcollection%SquickAccess(1),&
            rproblemLevel, p_rsolutionEuler, rcollection)
      end if

      ! Compute the preconditioner
      call transp_calcPreconditioner(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)
      
      ! Compute the residual
      call transp_calcResidual(rproblemLevel, rtimestep, rsolver,&
          rsolution, rsolutionInitial, rrhs, rres, istep, rcollection)

      ! Remove specifier for the preconditioner (if any)
      iSpec = iand(iSpec, not(NLSOL_OPSPEC_CALCPRECOND))
    end if
    

    ! Do we have to calculate the preconditioner?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then
     
      ! Compute the preconditioner
      call transp_calcPreconditioner(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)
    end if

    
    ! Do we have to calculate the Jacobian operator?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCJACOBIAN) .ne. 0) then
      
      ! Compute the Jacobian matrix
      call transp_calcJacobian(rproblemLevel, rtimestep, rsolver,&
          rsolution, rsolutionInitial, rcollection)
    end if
    
    
    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      
      ! Impose boundary conditions
      call transp_setBoundaryConditions(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolutionInitial, rres, rcollection)
    end if
    

    ! Do we have to apply the Jacobian operator?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_APPLYJACOBIAN) .ne. 0) then

      p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'jacobianMatrix', jacobianMatrix)

      ! Apply Jacobian matrix
      call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(jacobianMatrix),&
          rsolution%RvectorBlock(1), rres%RvectorBlock(1), 1.0_DP,&
          1.0_DP)
    end if
    
    
    ! Set status flag
    istatus = 0

  end subroutine zpinch_nlsolverCallbackTransport

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_initVelocityField(rparlist, ssectionName,&
      rproblemLevel, rsolution, rcollection)

!<description>
    ! This subroutine initializes the velocity field from the solution
    ! of the compressible Euler model. The result is stored separately
    ! for each problem level.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! solution vector of compressible Euler model
    type(t_vectorBlock), intent(in) :: rsolution
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: velocityfield
    integer :: neq, ndim
    
    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'velocityfield', velocityfield)
    
    ! Get number of degrees of freedom and spatial dimension
    neq  = rproblemLevel%rtriangulation%NVT
    ndim = rproblemLevel%rtriangulation%ndim
    
    ! Create/resize velocity vector if required
    if (rproblemLevel%RvectorBlock(velocityfield)%NEQ .eq. 0) then
      call lsysbl_createVectorBlock(rproblemLevel&
          %rvectorBlock(velocityfield), neq, ndim, .true.)
    elseif (rproblemLevel%RvectorBlock(velocityfield)%NEQ .ne. neq*ndim) then
      call lsysbl_resizeVectorBlock(rproblemLevel&
          %rvectorBlock(velocityfield), neq, .true.)
    end if

    ! Set x-velocity, i.e., momentum in x-direction
    call euler_getVariable(rsolution, 'momentum_x', rproblemLevel&
        %RvectorBlock(velocityfield)%RvectorBlock(1))

    ! Set y-velocity, i.e., momentum in y-direction
    call euler_getVariable(rsolution, 'momentum_y', rproblemLevel&
        %RvectorBlock(velocityfield)%RvectorBlock(2))

    ! Set global solution vector as external vector for the transport model
    call transp_setVariable2d(rsolution%RvectorBlock(1), 3)
    
    ! Set update notification in problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     PROBLEV_MSPEC_UPDATE)

  end subroutine zpinch_initVelocityField

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_initDensityAveraging(rparlist,&
      ssectionNameTransport, rproblemlevel, rsolutionEuler,&
      rcollection)

!<description>
    ! This subroutine initializes the density averaged mass matrices
    ! for the transport model based on the solution from the Euler
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section names in parameter list
    character(LEN=*), intent(in) :: ssectionNameTransport

    ! solution vector for Euler model
    type(t_vectorBlock), intent(in) :: rsolutionEuler
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_trilinearform) :: rform
    type(t_vectorScalar) :: rvector
    integer :: lumpedMassMatrix, consistentMassMatrix

    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionNameTransport, 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionNameTransport, 'consistentmassmatrix', consistentMassMatrix)
    
    ! Get density distribution from the solution of the Euler model
    !  and create block vector which is attached to the collection
    call euler_getVariable(rsolutionEuler, 'density', rvector)
    
    ! We have variable coefficients
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff    = .true.

    ! Initialize the bilinear form
    rform%itermCount = 1
    rform%Dcoefficients(1)  = 1.0_DP
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    rform%Idescriptors(3,1) = DER_FUNC

    ! Create density averaged consistent mass matrix
    call trilf_buildMatrixScalar(rform, .true.,&
        rproblemLevel%Rmatrix(consistentMassMatrix), rvector)
    
    ! Create density averaged lumped mass matrix
    call lsyssc_duplicateMatrix(rproblemLevel&
        %Rmatrix(consistentMassMatrix), rproblemLevel &
        %Rmatrix(lumpedMassMatrix), LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
    
    call lsyssc_lumpMatrixScalar(rproblemLevel&
        %Rmatrix(lumpedMassMatrix), LSYSSC_LUMP_DIAG)
    
    ! Release temporal vector
    call lsyssc_releaseVector(rvector)
    
  end subroutine zpinch_initDensityAveraging

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcLinearizedFCT(rbdrCondEuler,&
      rbdrCondTransport, rproblemLevel, rtimestep, rsolutionEuler,&
      rsolutionTransport, rcollection)

!<description>
    ! This subroutine calculates the linearized FCT correction for the
    ! Euler model and the scalar transport model simultaneously
!</description>

!<input>
    ! boundary condition structures
    type(t_boundaryCondition), intent(in) :: rbdrCondEuler, rbdrCondTransport

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep
!</input>

!<inputoutput>
    ! solution vectors
    type(t_vectorBlock), intent(inout) :: rsolutionEuler, rsolutionTransport

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
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
    call buildFluxCons2d(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, NVAR2D, nedge, p_solEuler, rtimestep%dStep,&
        p_MC, p_ML, p_Cx, p_Cy, p_DataEuler, p_fluxEuler, p_fluxEuler0)

    call buildFlux2d(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix&
        %NEQ, nedge, p_solTransport, rtimestep%dStep, p_MC, p_ML,&
        p_Cx, p_Cy, p_DataTransport, p_fluxTransport, p_fluxTransport0)

    ! Build the correction factors
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, NVAR2D, nedge, p_ML, p_fluxEuler, p_fluxEuler0&
        , 1, p_alpha, p_solEuler)
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, NVAR2D, nedge, p_ML, p_fluxEuler, p_fluxEuler0&
        , 4, p_alpha, p_solEuler)
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, 1, nedge, p_ML, p_fluxTransport,&
        p_fluxTransport0, 1, p_alpha, p_solTransport)

    ! Apply correction to low-order solutions
    call applyCorrection(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, NVAR2D, nedge, p_ML, p_fluxEuler, p_alpha,&
        p_DataEuler, p_solEuler)
    call applyCorrection(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, 1, nedge, p_ML, p_fluxTransport, p_alpha,&
        p_DataTransport, p_solTransport)

    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCondEuler, rsolutionEuler,&
        rtimestep%dTime, euler_calcBoundaryvalues2d)

    call bdrf_filterVectorExplicit(rbdrCondTransport,&
        rsolutionTransport, rtimestep%dTime)
    
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

    subroutine buildFluxCons2d(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        NEDGE, u, dscale, MC, ML, Cx, Cy, troc, flux, flux0)

      real(DP), dimension(NVAR,NEQ), intent(in) :: u
      real(DP), dimension(:), intent(in) :: MC,ML,Cx,Cy
      real(DP), intent(in) :: dscale
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NVAR,NEDGE
      
      integer, dimension(:), intent(inout) :: Ksep
      real(DP), dimension(NVAR,NEDGE), intent(inout) :: flux0,flux
      
      real(DP), dimension(NVAR,NEQ), intent(out) :: troc     
      
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

      real(DP), dimension(:), intent(in) :: MC,ML,Cx,Cy,u
      real(DP), intent(in) :: dscale
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE
      
      integer, dimension(:), intent(inout) :: Ksep
      real(DP), dimension(:), intent(inout) :: flux0,flux
      
      real(DP), dimension(:), intent(out) :: troc     
      
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
    
    subroutine  buildCorrectionCons(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        NVAR, NEDGE, ML, flux, flux0, ivar, alpha, u)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: flux0
      real(DP), dimension(:), intent(in) :: ML
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NVAR,NEDGE,ivar
      
      real(DP), dimension(NVAR,NEDGE), intent(inout) :: flux
      real(DP), dimension(NVAR,NEQ), intent(inout) :: u
      real(DP), dimension(:), intent(inout) :: alpha
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(:), allocatable :: pp,pm,qp,qm,rp,rm
      real(DP) :: f_ij,f0_ij,diff,aux,p_ij,p_ji,p0_ij,p0_ji,u_i,u_j,v_i,v_j,r_i,r_j
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
          
!!$          ! Flux correction in primitive variables
!!$          select case(ivar)
!!$          case default

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
          
            
!!$          case (4)
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
!!$            ! Solution differences
!!$            diff = (GAMMA-1) * (u(4,j) + &
!!$                         0.5 * (u_j*u_j + v_j*v_j)*u(1,j) -&
!!$                                u_j*u(2,j) - v_j*u(3,j) )&
!!$                 - (GAMMA-1) * (u(4,i) + &
!!$                         0.5 * (u_i*u_i + v_i*v_i)*u(1,i) -&
!!$                                u_i*u(2,i) - v_i*u(3,i) )
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
!!$            diff =  (GAMMA-1) * ( u(4,j)*u(1,j) - u(2,j)*u(2,j) &
!!$                                 -u(3,j)*u(3,j) + u(1,j)*u(4,j) )&
!!$                   -(GAMMA-1) * ( u(4,i)*u(1,i) - u(2,i)*u(2,i) &
!!$                                 -u(3,i)*u(3,i) + u(1,i)*u(4,i) )
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
!!$            qp(i) = max(qp(i),  diff)
!!$            qp(j) = max(qp(j), -diff)
!!$            qm(i) = min(qm(i),  diff)
!!$            qm(j) = min(qm(j), -diff)
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

!!$          ! Flux correction in primitive variables
!!$          select case(ivar)
!!$          case default
            
            ! Flux correction in conservative variables
            f_ij = flux(ivar,iedge)
            
            ! Limit conservative fluxes
            if (f_ij > SYS_EPSREAL) then
              alpha(iedge) = min(alpha(iedge), rp(i), rm(j))
            elseif (f_ij < -SYS_EPSREAL) then
              alpha(iedge) = min(alpha(iedge), rm(i), rp(j))
            end if

!!$          case (4)
!!$            
!!$            ! Flux correction in primitive variables
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
      
      real(DP), dimension(NVAR,NEDGE), intent(in) :: flux
      real(DP), dimension(:), intent(in) :: ML,alpha
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NVAR,NEDGE
      
      real(DP), dimension(NVAR,NEQ), intent(inout) :: u,data
      integer, dimension(:), intent(inout) :: Ksep
      
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
      real(DP), intent(in) :: a,b
      real(DP) :: minmod

      if (a > 0 .and. b > 0) then
        minmod = min(a,b)
      elseif (a < 0 .and. b < 0) then
        minmod = max(a,b)
      else
        minmod = 0
      end if
    end function minmod

  end subroutine zpinch_calcLinearizedFCT

end module zpinch_callback
