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
  use matrixfilters
  use matrixio
  
  use stdoperators
  use bilinearformevaluation
  
  use spatialoperators
  use spacetimevectors
  
  use physics
  use spacetimebc
  
  use fespacehierarchy
  use spacetimehierarchy
  
  implicit none

  ! Encapsules a space-time matrix
  type t_spaceTimeMatrix
  
    ! Type of the matrix.
    ! =0: standard matrix.
    ! =1: matrix of the Frechet-deerivative of the operator.
    integer :: cmatrixType
    
    ! Level of the matrix in the space-time hierarchy
    integer :: ilevel
    
    ! Underlying space-time hierarchy
    type(t_spacetimeHierarchy), pointer :: p_rspaceTimeHierarchy
  
    ! Template matrices for this space-time level
    type(t_matvecTemplates), dimension(:), pointer :: p_rmatVecTempl
    
    ! Physics of the problem.
    type(t_physics), pointer :: p_rphysics
    
    ! Boundary conditions
    type(t_spacetimeBC), pointer :: p_rboundaryCond
  
  end type


contains

  ! ***************************************************************************

  subroutine stmv_createMatrix (cmatrixType,ilevel,rspaceTimeHierarchy,rphysics,&
      rboundaryCond,rmatVecTempl,rmatrix)

  ! Creates a space-time matrix.
  
  ! Type of the matrix.
  ! =0: standard matrix.
  ! =1: matrix of the Frechet-deerivative of the operator.
  integer, intent(in) :: cmatrixType

  ! Level of the matrix in the space-time hierarchy
  integer, intent(in) :: ilevel
  
  ! Underlying space-time hierarchy
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierarchy

  ! Physics of the problem
  type(t_physics), intent(in), target :: rphysics

  ! Boundary condition structure.
  type(t_spacetimeBC), intent(in), target :: rboundaryCond
  
  ! Array of template matrices in space for all space-levels.
  type(t_matvecTemplates), dimension(:), intent(in), target :: rmatVecTempl
  
  ! Space-time matrix to be created.
  type(t_spaceTimeMatrix), intent(out) :: rmatrix
    
    rmatrix%cmatrixType = cmatrixType
    rmatrix%p_rphysics => rphysics
    rmatrix%ilevel = ilevel
    rmatrix%p_rspaceTimeHierarchy => rspaceTimeHierarchy
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
    nullify(rmatrix%p_rspaceTimeHierarchy)
    nullify(rmatrix%p_rmatVecTempl)
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
    integer :: istep,ispaceLevel
    type(t_vectorBlock) :: rtempVecY,rtempVecX1,rtempVecX2,rtempVecX3
    type(t_matrixBlock) :: rsubmatrix
    real(DP) :: dnormy
    real(DP), dimension(:), pointer :: p_Dx1, p_Dx2, p_Dx3, p_Dy
    type(t_feSpaceLevel), pointer :: p_rfeSpaceLevel
    
    ! Level in space of the space-time matrix
    call sth_getLevel(rmatrix%p_rspaceTimeHierarchy,rmatrix%ilevel,&
        p_rfeSpaceLevel=p_rfeSpaceLevel,ispaceLevel=ispaceLevel)
    
    ! Allocate a temp vectors and matrices
    call lsysbl_createVectorBlock(p_rfeSpaceLevel%p_rdiscretisation,rtempVecX1)
    call lsysbl_createVectorBlock(p_rfeSpaceLevel%p_rdiscretisation,rtempVecX2)
    call lsysbl_createVectorBlock(p_rfeSpaceLevel%p_rdiscretisation,rtempVecX3)
    call lsysbl_createVectorBlock(p_rfeSpaceLevel%p_rdiscretisation,rtempVecY)
    
    call lsysbl_getbase_double (rtempVecX1,p_Dx1)
    call lsysbl_getbase_double (rtempVecX2,p_Dx2)
    call lsysbl_getbase_double (rtempVecX3,p_Dx3)
    call lsysbl_getbase_double (rtempVecY,p_Dy)
    
    call stmv_allocSubmatrix (rmatrix%cmatrixType,rmatrix%p_rphysics,&
        rmatrix%p_rmatVecTempl(ispaceLevel),rsubmatrix)
    
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
        call stmv_getSubmatrix (rmatrix, ispaceLevel, istep, istep-1, rsubmatrix)
        call lsysbl_blockMatVec (rsubmatrix,rtempVecX1,rtempVecY,cx,1.0_DP)
      end if
      
      ! Diagonal
      call sptivec_getTimestepData(rx,istep,rtempVecX2)
      call stmv_getSubmatrix (rmatrix, ispaceLevel, istep, istep, rsubmatrix)
      call lsysbl_blockMatVec (rsubmatrix,rtempVecX2,rtempVecY,cx,1.0_DP)
      
      !call matio_writeBlockMatrixHR (rsubmatrix, "matrix",&
      !    .true., 0, "matrix", "(E20.10)")

      if (istep .lt. rx%NEQtime) then
        ! Right subdiagonal
        call sptivec_getTimestepData(rx,istep+1,rtempVecX3)
        call stmv_getSubmatrix (rmatrix, ispaceLevel, istep, istep+1, rsubmatrix)
        call lsysbl_blockMatVec (rsubmatrix,rtempVecX3,rtempVecY,cx,1.0_DP)
      end if
      
      call sptivec_setTimestepData(ry,istep,rtempVecY)

      dnormy = lsysbl_vectorNorm(rtempVecY,LINALG_NORML2)

    end do
    
    ! Release the temp vector.
    call lsysbl_releaseMatrix (rsubmatrix)
    call lsysbl_releaseVector (rtempVecX3)
    call lsysbl_releaseVector (rtempVecX2)
    call lsysbl_releaseVector (rtempVecX1)
    call lsysbl_releaseVector (rtempVecY)
  
  end subroutine

  ! ***************************************************************************

  subroutine stmv_allocSubmatrix (cmatrixType,rphysics,&
      rmatVecTempl,rsubmatrix,rdiscreteBC)

  ! Creates a space matrix.
  
  ! Type of the matrix.
  ! =0: standard matrix.
  ! =1: matrix of the Frechet-deerivative of the operator.
  integer, intent(in) :: cmatrixType
  
  ! Physics of the problem
  type(t_physics), intent(in), target :: rphysics

  ! Template matrices in space.
  type(t_matvecTemplates), intent(in), target :: rmatVecTempl
  
  ! Spatial destination matrix; memory has to be allocated in advance.
  type(t_matrixBlock), intent(out) :: rsubmatrix

  ! OPTIONAL: Boundary condition structure in space.
  type(t_discreteBC), intent(in), target, optional :: rdiscreteBC
  
    ! Create the basic matrix  
    call lsysbl_createMatBlockByDiscr (rmatvecTempl%p_rspaceDiscr,rsubmatrix)
    if (present(rdiscreteBC)) then
      rsubmatrix%p_rdiscreteBC => rdiscreteBC
    end if

    ! Create submatrices.
    !
    ! Type of equation?
    select case (rphysics%cequation)
    case (0,2)
      ! Heat equation
      
      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
    case (1)
      ! Stokes equation
      
      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(4,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)


      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(5,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateA11,&
          rsubmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)


      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateB,&
          rsubmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateB,&
          rsubmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateD,&
          rsubmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateD,&
          rsubmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)


      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateB,&
          rsubmatrix%RmatrixBlock(4,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateB,&
          rsubmatrix%RmatrixBlock(5,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateD,&
          rsubmatrix%RmatrixBlock(6,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateD,&
          rsubmatrix%RmatrixBlock(6,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)


      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateC,&
          rsubmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatVecTempl%rmatrixTemplateC,&
          rsubmatrix%RmatrixBlock(6,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

    case default
    
      call output_line ("Equation not supported.")
      call sys_halt()

    end select

  end subroutine

  ! ***************************************************************************

  subroutine stmv_reduceDiagToPrimal (rsubmatrix)
  
  ! Removes the coupling and dual part from a spatial matrix.
  
  ! Underlying space-time matrix
  type(t_matrixBlock), intent(inout) :: rsubmatrix
  
    integer :: i,k
  
    ! Decouple.
    k = rsubmatrix%nblocksPerRow
    
    rsubmatrix%RmatrixBlock(k/2+1:,1:k/2)%dscaleFactor = 0.0_DP
    rsubmatrix%RmatrixBlock(1:k/2,k/2+1:)%dscaleFactor = 0.0_DP
    rsubmatrix%RmatrixBlock(k/2+1:,k/2+1:)%dscaleFactor = 0.0_DP
  
    ! Initialise the other matrices with identities on the diagonal.
    do i=rsubmatrix%nblocksPerRow/2+1,rsubmatrix%nblocksPerRow
      rsubmatrix%RmatrixBlock(i,i)%dscaleFactor = 1.0_DP
      call lsyssc_clearMatrix (rsubmatrix%RmatrixBlock(i,i))
      call lsyssc_initialiseIdentityMatrix (rsubmatrix%RmatrixBlock(i,i))
    end do
  
  end subroutine
  
  ! ***************************************************************************

  subroutine stmv_reduceDiagToDual (rsubmatrix)
  
  ! Removes the coupling and primal part from a spatial matrix.
  
  ! Underlying space-time matrix
  type(t_matrixBlock), intent(inout) :: rsubmatrix
  
    integer :: i,k
  
    ! Decouple.
    k = rsubmatrix%nblocksPerRow
    
    rsubmatrix%RmatrixBlock(k/2+1:,1:k/2)%dscaleFactor = 0.0_DP
    rsubmatrix%RmatrixBlock(1:k/2,k/2+1:)%dscaleFactor = 0.0_DP
    rsubmatrix%RmatrixBlock(1:k/2,1:k/2)%dscaleFactor = 0.0_DP
  
    ! Initialise the other matrices with identities on the diagonal.
    do i=1,rsubmatrix%nblocksPerRow/2
      rsubmatrix%RmatrixBlock(i,i)%dscaleFactor = 1.0_DP
      call lsyssc_clearMatrix (rsubmatrix%RmatrixBlock(i,i))
      call lsyssc_initialiseIdentityMatrix (rsubmatrix%RmatrixBlock(i,i))
    end do
  
  end subroutine

  ! ***************************************************************************

  subroutine stmv_getSubmatrix (rmatrix, ispaceLevel, irow, icol, rsubmatrix)
  
  ! Assembles a sub-blockmatrix of the global matrix at position (irow,icol)
  
  ! Underlying space-time matrix
  type(t_spaceTimeMatrix), intent(in) :: rmatrix
  
  ! Level in space, corresponding to rsubmatrix.
  ! =0: Use the default level of the space-time matrix.
  integer, intent(in) :: ispaceLevel
  
  ! Row and column
  integer, intent(in) :: irow, icol
  
  ! Spatial destination matrix; memory has to be allocated in advance.
  type(t_matrixBlock), intent(inout) :: rsubmatrix
  
    ! local variables
    real(DP) :: dtheta, dtstep, dalpha, dgamma
    real(DP) :: dcoupleDualToPrimal,dcouplePrimalToDual,dcoupleTermCond
    integer :: ithetaschemetype, ilev
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
  
    ! Get the level in space
    ilev = ispaceLevel
    if (ilev .eq. 0) then
      call sth_getLevel(rmatrix%p_rspaceTimeHierarchy,rmatrix%ilevel,ispaceLevel=ilev)
    end if
  
    ! Clear the destination matrix.
    rsubMatrix%RmatrixBlock(:,:)%dscaleFactor = 1.0_DP
    call lsysbl_clearMatrix(rsubmatrix)
    
    call sth_getLevel(rmatrix%p_rspaceTimeHierarchy,rmatrix%ilevel,p_rtimeDiscr=p_rtimeDiscr)
  
    ! Matrix depends on the physics, on the time discretisation
    ! and on the matrix type (Frechet-derivative or direct operator)!
    select case (p_rtimeDiscr%ctype)
    
    case (TDISCR_ONESTEPTHETA)
    
      dtheta = p_rtimeDiscr%dtheta
      dtstep = p_rtimeDiscr%dtstep
      dalpha = rmatrix%p_rphysics%doptControlAlpha
      dgamma = rmatrix%p_rphysics%doptControlGamma
      dcoupleDualToPrimal = rmatrix%p_rphysics%dcoupleDualToPrimal
      dcouplePrimalToDual = rmatrix%p_rphysics%dcouplePrimalToDual
      dcoupleTermCond = rmatrix%p_rphysics%dcoupleTermCond
      ithetaschemetype = p_rtimeDiscr%itag
    
      ! Standard 1-step theta scheme.
      select case (rmatrix%p_rphysics%cequation)
      case (0,2)
        ! ###############################################################################
        ! Heat equation. Equation is linear, so the Frechet-derivative
        ! matches the standard matrix!
        ! ###############################################################################
        !
!        ! We have three cases...
!        if (irow .eq. 1) then
!          ! ##################
!          ! First timestep
!          ! ##################
!          
!          if (irow .eq. icol) then
!            ! -------------------------------
!            ! PRIMAL, DIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, &
!                (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!
!            ! -------------------------------
!            ! DUAL, DIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!
!            ! Coupling of the primal to the dual
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (dtheta-1.0_DP)*dcouplePrimalToDual,&
!                rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
!                
!          end if
!
!          if (irow+1 .eq. icol) then
!            ! -------------------------------
!            ! RIGHT OFFDIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!
!            ! Coupling of the dual to the primal
!            !call lsyssc_matrixLinearComb (&
!            !    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP-dtheta)*dcouplePrimalToDual,&
!            !    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
!            !    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
!          end if
!
!        else if (irow .eq. rmatrix%p_rtimeDiscr%nintervals+1) then
!
!          ! ##################
!          ! Last timestep
!          ! ##################
!
!          if (irow-1 .eq. icol) then
!            ! -------------------------------
!            ! LEFT OFFDIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!
!            ! Coupling of the primal to the dual
!            !call lsyssc_matrixLinearComb (&
!            !    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP-dtheta)/dalpha,&
!            !    rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
!            !    rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
!          end if
!          
!          if (irow .eq. icol) then
!            ! -------------------------------
!            ! PRIMAL, DIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!
!            ! Coupling of the dual to the primal
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP/dalpha),&
!                rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
!
!            ! -------------------------------
!            ! DUAL, DIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!
!            ! Coupling of the primal to the dual
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual*(dcoupleTermCond*dtheta+dgamma/dtstep),&
!                rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
!                
!          end if
!
!        else
        
          ! ######################
          ! Intermediate timestep.
          ! ######################

          if (irow-1 .eq. icol) then
            ! -------------------------------
            ! LEFT OFFDIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            if (irow .gt. 1) then
              if (ithetaschemetype .eq. 2) then
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(0.25_DP/dalpha),&
                    rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
                    rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
              end if
            end if

            ! Coupling of the primal to the dual
            !call lsyssc_matrixLinearComb (&
            !    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP-dtheta)/dalpha,&
            !    rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
            !    rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
          end if
          
          if (irow .eq. icol) then
            ! -------------------------------
            ! PRIMAL, DIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual.
            ! No coupling in the first timestep.
            if (irow .gt. 1) then
              if (ithetaschemetype .eq. 2) then
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(0.5_DP/dalpha),&
                    rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
                    rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
              else
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP/dalpha),&
                    rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
                    rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
              end if
            end if

            ! -------------------------------
            ! DUAL, DIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual.
            ! Different weights in the last timestep(s).
            
            if (ithetaschemetype .eq. 1) then
              if (irow .eq. 1) then
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(1.0_DP-dtheta),&
                    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
              else if (irow .eq. p_rtimeDiscr%nintervals) then
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(1.0_DP+(1.0_DP-dtheta)*dgamma/dtstep),&
                    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
              else if (irow .eq. p_rtimeDiscr%nintervals+1) then
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(dcoupleTermCond*dtheta+dtheta*dgamma/dtstep),&
                    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
              else
                ! Standard weights.
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    dcouplePrimalToDual*(-1.0_DP),&
                    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
              end if
              
            else if (ithetaschemetype .eq. 3) then
              if (irow .eq. p_rtimeDiscr%nintervals+1) then
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(dcoupleTermCond*dtheta+dgamma/dtstep),&
                    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
              else
                ! Standard weights.
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    dcouplePrimalToDual*(-dtheta),&
                    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
              end if
            
            else
              if (irow .eq. p_rtimeDiscr%nintervals+1) then
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(dcoupleTermCond+dgamma/dtstep),&
                    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
              else
                ! Standard weights.
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    dcouplePrimalToDual*(-1.0_DP),&
                    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
              end if
            end if
                
          end if

          if (irow+1 .eq. icol) then
            ! -------------------------------
            ! RIGHT OFFDIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            if (irow .gt. 1) then
              if (ithetaschemetype .eq. 2) then
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(0.25_DP/dalpha),&
                    rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
                    rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
              end if
            end if
            
            if (ithetaschemetype .eq. 3) then
              ! Standard weights.
              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                  dcouplePrimalToDual*(-(1.0_DP-dtheta)),&
                  rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                  rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
            end if
            

            ! Coupling of the primal to the dual
            !call lsyssc_matrixLinearComb (&
            !    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual*(1.0_DP-dtheta),&
            !    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
            !    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
          end if

!        end if
        
        ! DEBUG!!!
        !call lsysbl_scalematrix (rsubMatrix,dtstep)

      case (1)
        ! ###############################################################################
        ! Stokes equation. Equation is linear, so the Frechet-derivative
        ! matches the standard matrix!
        ! ###############################################################################
        !
!        ! We have three cases...
!        if (irow .eq. 1) then
!          ! ##################
!          ! First timestep
!          ! ##################
!          
!          if (irow .eq. icol) then
!            ! -------------------------------
!            ! PRIMAL, DIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, &
!                (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, &
!                (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!                
!            ! B/D
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixB1, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,3),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,3),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixB2, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,3),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,3),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixD1, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(3,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(3,1),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixD2, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(3,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(3,2),.false.,.false.,.true.,.true.)
!
!            ! -------------------------------
!            ! DUAL, DIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(4,4),1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,4),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(5,5),1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,5),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(4,4),1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,4),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(5,5),1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,5),.false.,.false.,.true.,.true.)
!
!            ! Coupling of the primal to the dual
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (dtheta-1.0_DP)*dcouplePrimalToDual,&
!                rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (dtheta-1.0_DP)*dcouplePrimalToDual,&
!                rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
!                
!            ! B/D
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixB1, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,6),1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,6),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixB2, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,6),1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,6),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixD1, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(6,4),1.0_DP,&
!                rsubMatrix%RmatrixBlock(6,4),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixD2, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(6,5),1.0_DP,&
!                rsubMatrix%RmatrixBlock(6,5),.false.,.false.,.true.,.true.)
!
!          end if
!
!          if (irow+1 .eq. icol) then
!            ! -------------------------------
!            ! RIGHT OFFDIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
!                rsubMatrix%RmatrixBlock(4,4),1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,4),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
!                rsubMatrix%RmatrixBlock(5,5),1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,5),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(4,4),1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,4),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(5,5),1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,5),.false.,.false.,.true.,.true.)
!
!            ! Coupling of the primal to the dual
!            !call lsyssc_matrixLinearComb (&
!            !    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP-dtheta)*dcouplePrimalToDual,&
!            !    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
!            !    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
!          end if
!
!        else if (irow .eq. rmatrix%p_rtimeDiscr%nintervals+1) then
!
!          ! ##################
!          ! Last timestep
!          ! ##################
!
!          if (irow-1 .eq. icol) then
!            ! -------------------------------
!            ! LEFT OFFDIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!                
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!
!            ! Coupling of the primal to the dual
!            !call lsyssc_matrixLinearComb (&
!            !    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP-dtheta)/dalpha,&
!            !    rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
!            !    rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
!          
!          end if
!          
!          if (irow .eq. icol) then
!            ! -------------------------------
!            ! PRIMAL, DIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
!
!            ! Coupling of the dual to the primal
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP/dalpha),&
!                rsubMatrix%RmatrixBlock(1,4),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,4),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP/dalpha),&
!                rsubMatrix%RmatrixBlock(2,5),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,5),.false.,.false.,.true.,.true.)
!
!            ! B/D
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixB1, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,3),1.0_DP,&
!                rsubMatrix%RmatrixBlock(1,3),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixB2, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,3),1.0_DP,&
!                rsubMatrix%RmatrixBlock(2,3),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixD1, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(3,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(3,1),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixD2, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(3,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(3,2),.false.,.false.,.true.,.true.)
!
!            ! -------------------------------
!            ! DUAL, DIAGONAL
!            ! -------------------------------
!            
!            ! Laplace
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(4,4),1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,4),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
!                rsubMatrix%RmatrixBlock(5,5),1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,5),.false.,.false.,.true.,.true.)
!                
!            ! Mass
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(4,4),1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,4),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
!                rsubMatrix%RmatrixBlock(5,5),1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,5),.false.,.false.,.true.,.true.)
!
!            ! Coupling of the primal to the dual
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual*(dcoupleTermCond*dtheta+dgamma/dtstep),&
!                rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual*(dcoupleTermCond*dtheta+dgamma/dtstep),&
!                rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
!                
!            ! B/D
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixB1, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,6),1.0_DP,&
!                rsubMatrix%RmatrixBlock(4,6),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixB2, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,6),1.0_DP,&
!                rsubMatrix%RmatrixBlock(5,6),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixD1, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(6,4),1.0_DP,&
!                rsubMatrix%RmatrixBlock(6,4),.false.,.false.,.true.,.true.)
!
!            call lsyssc_matrixLinearComb (&
!                rmatrix%p_rmatVecTempl(ilev)%rmatrixD2, 1.0_DP,&
!                rsubMatrix%RmatrixBlock(6,5),1.0_DP,&
!                rsubMatrix%RmatrixBlock(6,5),.false.,.false.,.true.,.true.)
!                
!          end if
!
!        else
        
          ! ######################
          ! Intermediate timestep.
          ! ######################

          if (irow-1 .eq. icol) then
            ! -------------------------------
            ! LEFT OFFDIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual
            !call lsyssc_matrixLinearComb (&
            !    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP-dtheta)/dalpha,&
            !    rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
            !    rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
          end if
          
          if (irow .eq. icol) then
            ! -------------------------------
            ! PRIMAL, DIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)
                
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(1,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(2,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual.
            ! No coupling in the first timestep.
            
            if (irow .gt. 1) then
              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP/dalpha),&
                  rsubMatrix%RmatrixBlock(1,4),1.0_DP,&
                  rsubMatrix%RmatrixBlock(1,4),.false.,.false.,.true.,.true.)

              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP/dalpha),&
                  rsubMatrix%RmatrixBlock(2,5),1.0_DP,&
                  rsubMatrix%RmatrixBlock(2,5),.false.,.false.,.true.,.true.)
            end if

            ! B/D
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixB1, 1.0_DP,&
                rsubMatrix%RmatrixBlock(1,3),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,3),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixB2, 1.0_DP,&
                rsubMatrix%RmatrixBlock(2,3),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,3),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixD1, 1.0_DP,&
                rsubMatrix%RmatrixBlock(3,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(3,1),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixD2, 1.0_DP,&
                rsubMatrix%RmatrixBlock(3,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(3,2),.false.,.false.,.true.,.true.)

            ! -------------------------------
            ! DUAL, DIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(4,4),1.0_DP,&
                rsubMatrix%RmatrixBlock(4,4),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (rmatrix%p_rphysics%dviscosity*dtheta),&
                rsubMatrix%RmatrixBlock(5,5),1.0_DP,&
                rsubMatrix%RmatrixBlock(5,5),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(4,4),1.0_DP,&
                rsubMatrix%RmatrixBlock(4,4),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, (1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(5,5),1.0_DP,&
                rsubMatrix%RmatrixBlock(5,5),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual.
            ! Different weights in the last timestep(s).
            if (ithetaschemetype .eq. 1) then
              if (irow .eq. 1) then
                
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(1.0_DP-dtheta),&
                    rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)

                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(1.0_DP-dtheta),&
                    rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
                    rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
              
              else if (irow .eq. p_rtimeDiscr%nintervals) then
                
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(1.0_DP+(1.0_DP-dtheta)*dgamma/dtstep),&
                    rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)

                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(1.0_DP+(1.0_DP-dtheta)*dgamma/dtstep),&
                    rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
                    rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
              
              else if (irow .eq. p_rtimeDiscr%nintervals+1) then
              
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(dcoupleTermCond*dtheta+dtheta*dgamma/dtstep),&
                    rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)

                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(dcoupleTermCond*dtheta+dtheta*dgamma/dtstep),&
                    rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
                    rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
                    
              else             
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual,&
                    rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)

                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual,&
                    rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
                    rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
              end if
              
            else
              
              if (irow .eq. p_rtimeDiscr%nintervals+1) then
              
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(dcoupleTermCond+dgamma/dtstep),&
                    rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)

                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                    -dcouplePrimalToDual*(dcoupleTermCond+dgamma/dtstep),&
                    rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
                    rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
                    
              else             
              
                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual,&
                    rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
                    rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)

                call lsyssc_matrixLinearComb (&
                    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual,&
                    rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
                    rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
                    
              end if
            end if
                
            ! B/D
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixB1, 1.0_DP,&
                rsubMatrix%RmatrixBlock(4,6),1.0_DP,&
                rsubMatrix%RmatrixBlock(4,6),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixB2, 1.0_DP,&
                rsubMatrix%RmatrixBlock(5,6),1.0_DP,&
                rsubMatrix%RmatrixBlock(5,6),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixD1, 1.0_DP,&
                rsubMatrix%RmatrixBlock(6,4),1.0_DP,&
                rsubMatrix%RmatrixBlock(6,4),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixD2, 1.0_DP,&
                rsubMatrix%RmatrixBlock(6,5),1.0_DP,&
                rsubMatrix%RmatrixBlock(6,5),.false.,.false.,.true.,.true.)

          end if

          if (irow+1 .eq. icol) then
            ! -------------------------------
            ! RIGHT OFFDIAGONAL
            ! -------------------------------
            
            ! Laplace
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(4,4),1.0_DP,&
                rsubMatrix%RmatrixBlock(4,4),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, rmatrix%p_rphysics%dviscosity*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(5,5),1.0_DP,&
                rsubMatrix%RmatrixBlock(5,5),.false.,.false.,.true.,.true.)
                
            ! Mass
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(4,4),1.0_DP,&
                rsubMatrix%RmatrixBlock(4,4),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -(1.0_DP/dtstep),&
                rsubMatrix%RmatrixBlock(5,5),1.0_DP,&
                rsubMatrix%RmatrixBlock(5,5),.false.,.false.,.true.,.true.)

            ! Coupling of the primal to the dual
            !call lsyssc_matrixLinearComb (&
            !    rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual*(1.0_DP-dtheta),&
            !    rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
            !    rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
          end if

!        end if

      case default
      
        call output_line ("Equation not supported.")
        call sys_halt()
        
      end select

    end select
    
    if (irow .eq. icol) then
    
      ! The matrix is a diagonal matrix in the supermatrix.
      rsubmatrix%imatrixSpec = iand(rsubmatrix%imatrixSpec,not(LSYSBS_MSPEC_OFFDIAGSUBMATRIX))
      
    else
      
      ! The matrix is an offdiagonalmatrix in the supermatrix.
      rsubmatrix%imatrixSpec = ior(rsubmatrix%imatrixSpec,LSYSBS_MSPEC_OFFDIAGSUBMATRIX)
    
    end if
    
    ! DEBUG: Scale the matrix by dtstep.
    !call lsysbl_scaleMatrix (rsubMatrix,dtstep)
    
  end subroutine

  ! ***************************************************************************

  subroutine stmv_implementDefBCSubmatrix (rboundaryCond, ispaceLevel, rsubmatrix, irow, icol, rdiscreteBC)
  
  ! Implements boundary conditions into a submatrix of the global matrix.
  
  ! Boundary conditions.
  type(t_spacetimeBC), intent(in), target :: rboundaryCond

  ! Space-level corresponding to rsubmatrix.
  ! =0: Use the space level of the current space-time level in rboundaryCond.
  integer, intent(in) :: ispaceLevel

  ! Spatial submatrix of the global space-time matrix
  type(t_matrixBlock), intent(inout) :: rsubmatrix
  
  ! Row and column in the global space-time matrix corresponding to
  ! rsubmatrix.
  integer, intent(in) :: irow, icol
  
  ! Temporary boundary condition structure.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
  
    call bcasm_clearDiscreteBC (rdiscreteBC)
    call spop_assembleSpaceBC (rboundaryCond, ispaceLevel, irow, SPOP_DEFECT, rdiscreteBC)
    call matfil_discreteBC (rsubmatrix,rdiscreteBC)

  end subroutine
  
end module
