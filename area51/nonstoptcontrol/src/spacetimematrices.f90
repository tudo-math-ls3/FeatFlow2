!##############################################################################
!# ****************************************************************************
!# <name> spacetimematrices </name>
!# ****************************************************************************
!#
!# <purpose>
!# Spatial operators used to create space-time matrices.
!# </purpose>
!##############################################################################

module spacetimematrices

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
  use spacetimebc
  
  use physics
  
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

  subroutine stmat_createMatrix (cmatrixType,ilevel,rspaceTimeHierarchy,rphysics,&
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

  subroutine stmat_releaseMatrix (rmatrix)

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

  subroutine stmat_allocSubmatrix (cmatrixType,rphysics,&
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

  subroutine stmat_reduceDiagToPrimal (rsubmatrix)
  
  ! Removes the coupling and dual part from a spatial matrix.
  
  ! Underlying space-time matrix
  type(t_matrixBlock), intent(inout) :: rsubmatrix
  
    integer :: i,j,k
  
    ! Decouple.
    k = rsubmatrix%nblocksPerRow
    
    ! All matrix except the primal diagonal block := 0.
    ! Dual diagonal block = identity.
    do i=1,k
      do j=1,k
        if (((i .ge. k/2+1) .or. (j .ge. k/2+1)) .and. &
            lsysbl_isSubmatrixPresent (rsubmatrix,i,j)) then
          call lsyssc_clearMatrix (rsubmatrix%RmatrixBlock(i,j))
          if (i .eq. j) then
            call lsyssc_initialiseIdentityMatrix (rsubmatrix%RmatrixBlock(i,i))
          end if
        end if
      end do
    end do
    
  end subroutine
  
  ! ***************************************************************************

  subroutine stmat_reduceDiagToDual (rsubmatrix)
  
  ! Removes the coupling and primal part from a spatial matrix.
  
  ! Underlying space-time matrix
  type(t_matrixBlock), intent(inout) :: rsubmatrix
  
    integer :: i,j,k
  
    ! Decouple.
    k = rsubmatrix%nblocksPerRow
    
    ! All matrix except the dual diagonal block := 0.
    ! Primal diagonal block = identity.
    do i=1,k
      do j=1,k
        if (((i .le. k/2) .or. (j .le. k/2)) .and. &
            lsysbl_isSubmatrixPresent (rsubmatrix,i,j)) then
          call lsyssc_clearMatrix (rsubmatrix%RmatrixBlock(i,j))
          if (i .eq. j) then
            call lsyssc_initialiseIdentityMatrix (rsubmatrix%RmatrixBlock(i,i))
          end if
        end if
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

  subroutine stmat_initDualIdentity (rsubmatrix)
  
  ! Removes the coupling of the dual from the primal.
  ! Initialises the dual main part by an identity matrix.
  
  ! Underlying space-time matrix
  type(t_matrixBlock), intent(inout) :: rsubmatrix
  
    integer :: i,k
  
    ! Decouple.
    k = rsubmatrix%nblocksPerRow
    
    ! Switch off all dual matrices.
    rsubmatrix%RmatrixBlock(k/2+1:,:)%dscaleFactor = 0.0_DP
  
    ! Switch on only the diagonal matrices. Initialise with identity.
    do i=rsubmatrix%nblocksPerRow/2+1,rsubmatrix%nblocksPerRow
      rsubmatrix%RmatrixBlock(i,i)%dscaleFactor = 1.0_DP
      call lsyssc_clearMatrix (rsubmatrix%RmatrixBlock(i,i))
      call lsyssc_initialiseIdentityMatrix (rsubmatrix%RmatrixBlock(i,i))
    end do
  
  end subroutine
  
  ! ***************************************************************************

  subroutine stmat_initPrimalIdentity (rsubmatrix)
  
  ! Removes the coupling of the primal from the dual.
  ! Initialises the dual main part by an identity matrix.
  
  ! Underlying space-time matrix
  type(t_matrixBlock), intent(inout) :: rsubmatrix
  
    integer :: i,k
  
    ! Decouple.
    k = rsubmatrix%nblocksPerRow
    
    ! Switch off all dual matrices.
    rsubmatrix%RmatrixBlock(1:k/2,:)%dscaleFactor = 0.0_DP
  
    ! Switch on only the diagonal matrices. Initialise with identity.
    do i=1,rsubmatrix%nblocksPerRow/2
      rsubmatrix%RmatrixBlock(i,i)%dscaleFactor = 1.0_DP
      call lsyssc_clearMatrix (rsubmatrix%RmatrixBlock(i,i))
      call lsyssc_initialiseIdentityMatrix (rsubmatrix%RmatrixBlock(i,i))
    end do
  
  end subroutine
  
  ! ***************************************************************************

  subroutine stmat_getSubmatrix (rmatrix, ispaceLevel, irow, icol, rsubmatrix)
  
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

          if (ithetaschemetype .eq. 0) then
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*((1.0_DP-dtheta)/dalpha),&
                rsubMatrix%RmatrixBlock(1,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,2),.false.,.false.,.true.,.true.)
          end if

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
            if (ithetaschemetype .eq. 0) then
              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(dtheta/dalpha),&
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
                  dcouplePrimalToDual*(-dtheta),&
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

          if (ithetaschemetype .eq. 0) then
            ! Coupling of the primal to the dual
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(2,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,1),.false.,.false.,.true.,.true.)
          end if

        end if
      
      ! DEBUG!!!
      !call lsysbl_scalematrix (rsubMatrix,dtstep)

      case (1)
        ! ###############################################################################
        ! Stokes equation. Equation is linear, so the Frechet-derivative
        ! matches the standard matrix!
        ! ###############################################################################
        !
        
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
          if (ithetaschemetype .eq. 0) then
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP-dtheta)/dalpha,&
                rsubMatrix%RmatrixBlock(1,4),1.0_DP,&
                rsubMatrix%RmatrixBlock(1,4),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP-dtheta)/dalpha,&
                rsubMatrix%RmatrixBlock(2,5),1.0_DP,&
                rsubMatrix%RmatrixBlock(2,5),.false.,.false.,.true.,.true.)
          end if
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
            if (ithetaschemetype .eq. 0) then
            
              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(dtheta/dalpha),&
                  rsubMatrix%RmatrixBlock(1,4),1.0_DP,&
                  rsubMatrix%RmatrixBlock(1,4),.false.,.false.,.true.,.true.)

              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(dtheta/dalpha),&
                  rsubMatrix%RmatrixBlock(2,5),1.0_DP,&
                  rsubMatrix%RmatrixBlock(2,5),.false.,.false.,.true.,.true.)

            else
            
              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP/dalpha),&
                  rsubMatrix%RmatrixBlock(1,4),1.0_DP,&
                  rsubMatrix%RmatrixBlock(1,4),.false.,.false.,.true.,.true.)

              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, dcoupleDualToPrimal*(1.0_DP/dalpha),&
                  rsubMatrix%RmatrixBlock(2,5),1.0_DP,&
                  rsubMatrix%RmatrixBlock(2,5),.false.,.false.,.true.,.true.)
            
            end if
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
            
          else if (ithetaschemetype .eq. 0) then
            
            if (irow .eq. p_rtimeDiscr%nintervals+1) then
            
              ! Impose a divergence free projection here.
            
              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                  -dcouplePrimalToDual*(dgamma/dtstep),&
                  rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
                  rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)

              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, &
                  -dcouplePrimalToDual*(dgamma/dtstep),&
                  rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
                  rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
                  
              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (dgamma*rmatrix%p_rphysics%dviscosity*dtheta),&
                  rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
                  rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)

              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixLaplaceA11, (dgamma*rmatrix%p_rphysics%dviscosity*dtheta),&
                  rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
                  rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
                  
            else             
            
              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual*dtheta,&
                  rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
                  rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)

              call lsyssc_matrixLinearComb (&
                  rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual*dtheta,&
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
          if (ithetaschemetype .eq. 0) then
            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(4,1),1.0_DP,&
                rsubMatrix%RmatrixBlock(4,1),.false.,.false.,.true.,.true.)

            call lsyssc_matrixLinearComb (&
                rmatrix%p_rmatVecTempl(ilev)%rmatrixMassA11, -dcouplePrimalToDual*(1.0_DP-dtheta),&
                rsubMatrix%RmatrixBlock(5,2),1.0_DP,&
                rsubMatrix%RmatrixBlock(5,2),.false.,.false.,.true.,.true.)
          end if
        end if

!      end if

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

  subroutine stmat_implementDefBCSubmatrix (rboundaryCond, ispaceLevel, rsubmatrix, irow, icol, rdiscreteBC)
  
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
  
  ! ***************************************************************************

  subroutine stmat_applyOperator (rmatrix, irow, icol, cprimaldual, rvector, rdestVector)

  ! Applies the operator to a vector.
  
  ! Space-time matrix that defines the operators.
  type(t_spaceTimeMatrix), intent(in) :: rmatrix
  
  ! Row/column of the operator in the space-time matrix which should be applieed
  ! to rvector
  integer, intent(in) :: irow, icol
  
  ! Determines whether the primal and/or dual operator should be applied
  ! to the vector.
  ! =1: only primal, decoupled. =2: only dual, decoupled. 
  ! =3: primal + dual, both decoupled.
  integer, intent(in) :: cprimaldual
  
  ! Vector, the operator should be applied to
  type(t_vectorBlock), intent(inout) :: rvector

  ! Destination vector that receives the result. The content of the primal
  ! and/or dual part is overwritten, depending on cprimaldual
  type(t_vectorBlock), intent(inout) :: rdestVector
  
    ! local varibales
    type(t_vectorBlock) :: rvector2,rvector3
    type(t_matrixBlock) :: rsubmatrix
    integer :: ispaceLevel

    call lsysbl_createVectorBlock(rvector,rvector2,.true.)
    call lsysbl_createVectorBlock(rvector,rvector3,.true.)
    
    ! Multiply the solution with the matrix in the 1st timestep to
    ! get the correct RHS.
    call sth_getLevel(rmatrix%p_rspaceTimeHierarchy,rmatrix%ilevel,&
        ispaceLevel=ispaceLevel)
    call stmat_allocSubmatrix (rmatrix%cmatrixType,rmatrix%p_rphysics,&
      rmatrix%p_rmatVecTempl(ispaceLevel),rsubmatrix)
    call stmat_getSubmatrix (rmatrix, ispaceLevel, irow, icol, rsubmatrix)
    
    ! Primal part.
    if (iand(cprimaldual,1) .ne. 0) then
    
      ! Copy the primal part to the input vector, delete dual part.
      call lsysbl_clearVector (rvector2)
      select case (rmatrix%p_rphysics%cequation)
      case (0,2)
        ! 1D/2D Heat equation
        call lsyssc_copyVector (rvector%RvectorBlock(1),rdestVector%RvectorBlock(1))
        
      case (1)
        ! 2D Stokes equation
        call lsyssc_copyVector (rvector%RvectorBlock(1),rdestVector%RvectorBlock(1))
        call lsyssc_copyVector (rvector%RvectorBlock(2),rdestVector%RvectorBlock(2))
        call lsyssc_copyVector (rvector%RvectorBlock(3),rdestVector%RvectorBlock(3))
        
      end select
    
      call lsysbl_blockMatVec (rsubmatrix,rvector2,rvector3,1.0_DP,0.0_DP)
      
      ! Copy back.
      select case (rmatrix%p_rphysics%cequation)
      case (0,2)
        ! 1D/2D Heat equation
        call lsyssc_copyVector (rvector3%RvectorBlock(1),rvector%RvectorBlock(1))
        
      case (1)
        ! 2D Stokes equation
        call lsyssc_copyVector (rvector3%RvectorBlock(1),rvector%RvectorBlock(1))
        call lsyssc_copyVector (rvector3%RvectorBlock(2),rvector%RvectorBlock(2))
        call lsyssc_copyVector (rvector3%RvectorBlock(3),rvector%RvectorBlock(3))
        
      end select
      
    end if

    ! Dual part.
    if (iand(cprimaldual,2) .ne. 0) then
    
      ! Copy the dual part to the input vector, delete dual part.
      call lsysbl_clearVector (rvector2)
      select case (rmatrix%p_rphysics%cequation)
      case (0,2)
        ! 1D/2D Heat equation
        call lsyssc_copyVector (rvector%RvectorBlock(2),rvector2%RvectorBlock(2))
        
      case (1)
        ! 2D Stokes equation
        call lsyssc_copyVector (rvector%RvectorBlock(4),rvector2%RvectorBlock(4))
        call lsyssc_copyVector (rvector%RvectorBlock(5),rvector2%RvectorBlock(5))
        call lsyssc_copyVector (rvector%RvectorBlock(6),rvector2%RvectorBlock(6))
        
      end select
    
      call lsysbl_blockMatVec (rsubmatrix,rvector2,rvector3,1.0_DP,0.0_DP)
      
      ! Copy back.
      select case (rmatrix%p_rphysics%cequation)
      case (0,2)
        ! 1D/2D Heat equation
        call lsyssc_copyVector (rvector3%RvectorBlock(2),rdestVector%RvectorBlock(2))
        
      case (1)
        ! 2D Stokes equation
        call lsyssc_copyVector (rvector3%RvectorBlock(4),rdestVector%RvectorBlock(4))
        call lsyssc_copyVector (rvector3%RvectorBlock(5),rdestVector%RvectorBlock(5))
        call lsyssc_copyVector (rvector3%RvectorBlock(6),rdestVector%RvectorBlock(6))
        
      end select
      
    end if

    call lsysbl_releaseMatrix (rsubmatrix)
    call lsysbl_releaseVector (rvector3)
    call lsysbl_releaseVector (rvector2)

  end subroutine

end module
