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
!# 1.) zpinch_calcLinearizedFCT
!#     -> Calculates the linearized FCT correction
!#
!# </purpose>
!##############################################################################

module zpinch_callback

  use boundaryfilter
  use collection
  use euler_basic
  use euler_callback2d
  use fsystem
  use basicgeometry
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solveraux
  use storage
  use timestepaux
  use transport_callback2d

  implicit none

  private
  public :: zpinch_calcLinearizedFCT
  
contains

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
    call bdrf_filterVectorExplicit(rbdrCondEuler,&
                                   rsolutionEuler, rtimestep%dTime,&
                                   euler_calcBoundaryvalues2d)

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

    subroutine buildFluxCons2d(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, NEDGE, u,&
                               dscale, MC, ML, Cx, Cy, troc, flux, flux0)

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
    
    subroutine  buildCorrectionCons(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, NEDGE,&
                                    ML, flux, flux0, ivar, alpha, u)

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
