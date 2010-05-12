!##############################################################################
!# ****************************************************************************
!# <name> spacetimematvec </name>
!# ****************************************************************************
!#
!# <purpose>
!# Application of and with space-time matrices.
!# </purpose>
!##############################################################################

module spacetimematvec

  use fsystem
  use storage
  use genoutput
  use linearalgebra
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use collection
  use vectorio
  use derivatives
  use timediscretisation
  
  use stdoperators
  use bilinearformevaluation
  
  use spatialoperators
  use spacetimevectors
  
  use physics
  use spacetimebc
  
  implicit none

  ! Encapsules a space-time matrix
  type t_spaceTimeMatrix
  
    ! Type of the matrix.
    ! =0: standard matrix.
    ! =1: matrix of the Frechet-deerivative of the operator.
    integer :: cmatrixType
  
    ! Underlying space discretisation
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr
    
    ! Underlying time discretisation
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ! Template matrices for this space-time level
    type(t_matvecTemplates), pointer :: p_rmatVecTempl
    
    ! Physics of the problem.
    type(t_physics), pointer :: p_rphysics
    
    ! Boundary conditions
    type(t_spacetimeBC), pointer :: p_rboundaryCond
  
  end type


contains

  ! ***************************************************************************

  subroutine stmv_createMatrix (cmatrixType,rspaceDiscr,rtimeDiscr,rphysics,&
      rboundaryCond,rmatVecTempl,rmatrix)

  ! Creates a space-time matrix.
  
  ! Type of the matrix.
  ! =0: standard matrix.
  ! =1: matrix of the Frechet-deerivative of the operator.
  integer, intent(in) :: cmatrixType
  
  ! Underlying spatial discretisation structure
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr

  ! Underlying time discretisation structure
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr

  ! Physics of the problem
  type(t_physics), intent(in), target :: rphysics

  ! Boundary condition structure.
  type(t_spacetimeBC), intent(in), target :: rboundaryCond
  
  ! Template matrices in space.
  type(t_matvecTemplates), intent(in), target :: rmatVecTempl
  
  ! Space-time matrix to be created.
  type(t_spaceTimeMatrix), intent(out) :: rmatrix
    
    rmatrix%cmatrixType = cmatrixType
    rmatrix%p_rphysics => rphysics
    rmatrix%p_rspaceDiscr => rspaceDiscr
    rmatrix%p_rtimeDiscr => rtimeDiscr
    rmatrix%p_rboundaryCond => rboundaryCond
    rmatrix%p_rmatVecTempl => rmatVecTempl
    
  end subroutine

  ! ***************************************************************************

  subroutine stmv_releaseMatrix (rmatrix)

  ! Releases a space-time matrix.
  
  ! Space-time matrix to be created.
  type(t_spaceTimeMatrix), intent(inout) :: rmatrix
    
    rmatrix%cmatrixType = 0
    nullify(rmatrix%p_rphysics)
    nullify(rmatrix%p_rspaceDiscr)
    nullify(rmatrix%p_rtimeDiscr)
    nullify(rmatrix%p_rboundaryCond)
    
  end subroutine

  ! ***************************************************************************

  subroutine stmv_matvec (rmatrix, rx, ry, cx, cy)
  
  ! Matrix-vector multiplication with a space-time matrix:
  ! ry = cx*A*rx + cy*ry
  
  ! Underlying space-time matrix
  type(t_spaceTimeMatrix), intent(in) :: rmatrix
  
  ! Input vector x.
  type(t_spaceTimeVector), intent(in) :: rx
  
  ! Input and output vector y.
  type(t_spaceTimeVector), intent(inout) :: ry
  
  ! Weight for x.
  real(DP), intent(in) :: cx
  
  ! Weioght for y
  real(DP), intent(in) :: cy
  
    ! local variables
    integer :: istep
    type(t_vectorBlock) :: rtempVecY,rtempVecX1,rtempVecX2,rtempVecX3
    type(t_matrixBlock) :: rsubmatrix
    real(DP), dimension(:), pointer :: p_Dx1, p_Dx2, p_Dx3, p_Dy
    
    ! Allocate a temp vectors and matrices
    call lsysbl_createVectorBlock(rmatrix%p_rspaceDiscr,rtempVecX1)
    call lsysbl_createVectorBlock(rmatrix%p_rspaceDiscr,rtempVecX2)
    call lsysbl_createVectorBlock(rmatrix%p_rspaceDiscr,rtempVecX3)
    call lsysbl_createVectorBlock(rmatrix%p_rspaceDiscr,rtempVecY)
    
    call lsysbl_getbase_double (rtempVecX1,p_Dx1)
    call lsysbl_getbase_double (rtempVecX2,p_Dx2)
    call lsysbl_getbase_double (rtempVecX3,p_Dx3)
    call lsysbl_getbase_double (rtempVecY,p_Dy)
    
    call spop_allocateMatrix (rmatrix%p_rmatVecTempl,rmatrix%p_rphysics,rsubmatrix)
    
    ! Loop over all timesteps
    do istep=1,rx%NEQtime
    
      ! Load the RHS into the temp vector. Scale it by cy; we cannot 
      ! include cy in the space-time MV below since tghe weights would
      ! be wrong for subsequent MV's otherwise.
      call sptivec_getTimestepData(ry,istep,rtempVecY)
      call lsysbl_scaleVector (rtempVecY,cy)
      
      ! There are three MV's in each timestep...
      
      if (istep .gt. 1) then
        ! Left subdiagonal
        call sptivec_getTimestepData(rx,istep-1,rtempVecX1)
        call stmv_getSubmatrix (rmatrix, istep, istep-1, rsubmatrix)
        call lsysbl_blockMatVec (rsubmatrix,rtempVecX1,rtempVecY,cx,1.0_DP)
      end if
      
      ! Diagonal
      call sptivec_getTimestepData(rx,istep,rtempVecX2)
      call stmv_getSubmatrix (rmatrix, istep, istep, rsubmatrix)
      call lsysbl_blockMatVec (rsubmatrix,rtempVecX2,rtempVecY,cx,1.0_DP)

      if (istep .lt. rx%NEQtime) then
        ! Right subdiagonal
        call sptivec_getTimestepData(rx,istep+1,rtempVecX3)
        call stmv_getSubmatrix (rmatrix, istep, istep+1, rsubmatrix)
        call lsysbl_blockMatVec (rsubmatrix,rtempVecX3,rtempVecY,cx,1.0_DP)
      end if
      
      call sptivec_setTimestepData(ry,istep,rtempVecY)

    end do
    
    ! Release the temp vector.
    call lsysbl_releaseMatrix (rsubmatrix)
    call lsysbl_releaseVector (rtempVecX3)
    call lsysbl_releaseVector (rtempVecX2)
    call lsysbl_releaseVector (rtempVecX1)
    call lsysbl_releaseVector (rtempVecY)
  
  end subroutine

  ! ***************************************************************************

  subroutine stmv_allocSubmatrix (cmatrixType,rspaceDiscr,rphysics,&
      rdiscreteBC,rmatVecTempl,rsubmatrix)

  ! Creates a space-time matrix.
  
  ! Type of the matrix.
  ! =0: standard matrix.
  ! =1: matrix of the Frechet-deerivative of the operator.
  integer, intent(in) :: cmatrixType
  
  ! Underlying spatial discretisation structure
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr

  ! Physics of the problem
  type(t_physics), intent(in), target :: rphysics

  ! Boundary condition structure in space.
  type(t_discreteBC), intent(in), target :: rdiscreteBC
  
  ! Template matrices in space.
  type(t_matvecTemplates), intent(in), target :: rmatVecTempl
  
  ! Spatial destination matrix; memory has to be allocated in advance.
  type(t_matrixBlock), intent(out) :: rsubmatrix

    ! Create the basic matrix  
    call lsysbl_createMatBlockByDiscr (rspaceDiscr,rsubmatrix)
    rsubmatrix%p_rdiscreteBC => rdiscreteBC

    ! Create submatrices.
    !
    ! Type of equation?
    select case (rphysics%cequation)
    case (0)
      ! Heat equation
      
      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
    end select

  end subroutine

  ! ***************************************************************************

  subroutine stmv_getSubmatrix (rmatrix, irow, icol, rsubmatrix)
  
  ! Assembles a sub-blockmatrix of the global matrix at position (irow,icol)
  
  ! Underlying space-time matrix
  type(t_spaceTimeMatrix), intent(in) :: rmatrix
  
  ! Row and column
  integer, intent(in) :: irow, icol
  
  ! Spatial destination matrix; memory has to be allocated in advance.
  type(t_matrixBlock), intent(inout) :: rsubmatrix
  
    ! local variables
    real(DP) :: dtheta, dtstep, dalpha, dgamma
    real(DP) :: dcoupleDualToPrimal,dcouplePrimalToDual
  
    ! Clear the destination matrix.
    call lsysbl_clearMatrix(rsubmatrix)
  
    ! Matrix depends on the physics, on the time discretisation
    ! and on the matrix type (Frechet-derivative or direct operator)!
    select case (rmatrix%p_rtimeDiscr%ctype)
    
    case (TDISCR_ONESTEPTHETA)
    
      dtheta = rmatrix%p_rtimeDiscr%dtheta
      dtstep = rmatrix%p_rtimeDiscr%dtstep
      dalpha = rmatrix%p_rphysics%doptControlAlpha
      dgamma = rmatrix%p_rphysics%doptControlGamma
      dcoupleDualToPrimal = rmatrix%p_rphysics%dcoupleDualToPrimal
      dcouplePrimalToDual = rmatrix%p_rphysics%dcouplePrimalToDual
    
      ! Standard 1-step theta scheme.
      select case (rmatrix%p_rphysics%cequation)
      case (0)
        ! Heat equation. Equation is linear, so the Frechet-derivative
        ! matches the standard matrix!
        !
        ! We have three cases...
        if (irow .eq. 1) then
          ! ##################
          ! First timestep
          ! ##################
          
          if (irow .eq. icol) then
            ! -------------------------------
            ! PRIMAL DEFECT, DIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            ! -------------------------------
            ! DUAL DEFECT, DIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, (-1.0_DP)*dcouplePrimalToDual,&
                rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
                
          end if

          if (irow+1 .eq. icol) then
            ! -------------------------------
            ! RIGHT OFFDIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, -(1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual
            !call lsyssc_matrixLinearComb (&
            !    rmatrix%p_rmatVecTempl%rmatrixMassA11, -(1.0_DP-dtheta)*dcouplePrimalToDual,&
            !    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
            !    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
          end if

        else if (irow .eq. rmatrix%p_rtimeDiscr%nintervals+1) then

          ! ##################
          ! Last timestep
          ! ##################

          if (irow-1 .eq. icol) then
            ! -------------------------------
            ! LEFT OFFDIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, -(1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual
            !call lsyssc_matrixLinearComb (&
            !    rmatrix%p_rmatVecTempl%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP-dtheta)/dalpha,&
            !    rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
            !    rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
          end if
          
          if (irow .eq. icol) then
            ! -------------------------------
            ! PRIMAL DEFECT, DIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            ! Coupling of the dual to the primal
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP/dalpha),&
                rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)

            ! -------------------------------
            ! DUAL DEFECT, DIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, -dcouplePrimalToDual*(1.0_DP+dgamma/dtstep),&
                rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
                
          end if

        else
        
          ! ######################
          ! Intermediate timestep.
          ! ######################

          if (irow-1 .eq. icol) then
            ! -------------------------------
            ! LEFT OFFDIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, -(1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual
            !call lsyssc_matrixLinearComb (&
            !    rmatrix%p_rmatVecTempl%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP-dtheta)/dalpha,&
            !    rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
            !    rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
          end if
          
          if (irow .eq. icol) then
            ! -------------------------------
            ! PRIMAL DEFECT, DIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP/dalpha),&
                rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)

            ! -------------------------------
            ! DUAL DEFECT, DIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, dcouplePrimalToDual*(-1.0_DP),&
                rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
                
          end if

          if (irow+1 .eq. icol) then
            ! -------------------------------
            ! RIGHT OFFDIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl%rmatrixMassA11, -(1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual
            !call lsyssc_matrixLinearComb (&
            !    rmatrix%p_rmatVecTempl%rmatrixMassA11, -dcouplePrimalToDual*(1.0_DP-dtheta),&
            !    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
            !    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
          end if

        end if
        
      end select

    end select      
    
  end subroutine

end module
