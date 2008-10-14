!##############################################################################
!# ****************************************************************************
!# <name> l2projection </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to calculate $L_2$ projections of analytically
!# given functions, i.e. to calculate the Finite Element representation of
!# an analytic function.
!#
!# One can find the following subroutines here:
!#
!# 1.) l2prj_analytL2projectionByMass
!#     -> Performs defect correction with the consistent and lumped mass 
!#        matrix of a FE space to calculate the consistent $L_2$ projection
!#        of an analytically given function.
!# </purpose>
!##############################################################################

module l2projection

  use fsystem
  use linearalgebra
  use linearsystemscalar
  use linearformevaluation
  
  implicit none

!<types>

!<typeblock>

  ! Configuration block for the function anevl_L2projectionByMass which carries
  ! out the $L_2$ projection of an analytically given function by defect correction.
  type t_configL2ProjectionByMass
  
    ! Relative error criterium. Standard = $10^-{5}$.
    ! anevl_L2projectionByMass carries out the iteration until the relative as well
    ! the absolute error criterium is reached. 
    ! A value of 0.0 disables the check against the relative residuum.
    real(DP) :: depsRel = 1.0E-5_DP

    ! Absolute error criterium. Standard = $10^-{5}$.
    ! anevl_L2projectionByMass carries out the iteration until the relative as well
    ! the absolute error criterium is reached. 
    ! A value of 0.0 disables the check against the absolute residuum.
    real(DP) :: depsAbs = 1.0E-5

    ! Maximum iterations to be carried out. Standard = 10.
    ! The iteration stops prematurely if the number of iterations reaches this number.
    integer :: nmaxIterations = 10
    
    ! Type of norm to use for measuring errors. Standard is LINALG_NORML2.
    integer :: cnorm = LINALG_NORML2
    
    ! Damping parameter for the iteration. Standard = 1.0.
    real(DP) :: domega = 1.0_DP
    
    ! Output: Returns the initial residuum.
    ! This value is only set if depsRel > 0; otherwise, the relative error is
    ! not computed and left to 1.0_DP.
    real(DP) :: dinitResiduum = 0.0_DP

    ! Output: Returns the final relative error.
    ! This value is only set if depsRel > 0; otherwise, the relative error is
    ! not computed.
    real(DP) :: drelError = 0.0_DP

    ! Output: Returns the final absolute error.
    real(DP) :: dabsError = 0.0_DP

    ! Output: Returns the number of performed iterations
    integer :: iiterations = 0
  end type
  
!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine l2prj_analytL2projectionByMass (rvector,rmatrixMass,rmatrixMassLumped,&
      rvectorTemp1,rvectorTemp2,fcoeff_buildVectorSc_sim,rcollection,&
      rL2ProjectionConfig)
      
!<description>
  ! Converts an analytically given FE function fcoeff_buildVectorSc_sim
  ! to a finite element vector rvector by using mass matrices.
  ! So the resulting vector is the consistent $L_2$ projection of the
  ! analytically given function.
  ! The following iteration is performed for this purpose:
  !
  !    $$ x_{n+1}  :=  x_n  +  \omega M_l^{-1} ( f - M x_n ) $$
  !
  ! with $f$ being a RHS vector generated with fcoeff_buildVectorSc_sim,
  ! $M$ being the consistent mass matrix and $M_l$ being the
  ! lumped mass matrix. The iteration is performed until the error criteria or
  ! are fulfilled or the given number of steps is reached.
!</description>

!<input>
  ! The consistent mass matrix of the FE space of rvector.
  type(t_matrixScalar), intent(IN) :: rmatrixMass
  
  ! The lumped mass matrix of the FE space of rvector.
  type(t_matrixScalar), intent(IN) :: rmatrixMassLumped

  ! A callback routine for the function to be discretised. The callback routine
  ! has the same syntax as that for evaluating analytic functions for the 
  ! computation of RHS vectors.
  include 'intf_coefficientVectorSc.inc'

  ! OPTIONAL: A pointer to a collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(INOUT), target, optional :: rcollection

!</input>

!<inputoutput>
  ! A scalar vector that receives the $L_2$ projection of the function.
  type(t_vectorScalar), intent(INOUT) :: rvector

  ! A temporary vector of the same size as rvector.
  type(t_vectorScalar), intent(INOUT) :: rvectorTemp1

  ! A second temporary vector of the same size as rvector.
  type(t_vectorScalar), intent(INOUT) :: rvectorTemp2
  
  ! OPTIONAL: A configuration block for the iteration.
  ! If not specified, the standard settings are used.
  type(t_configL2ProjectionByMass), intent(INOUT), optional :: rL2ProjectionConfig
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_linearForm) :: rlinform
    integer :: iiteration
    type(t_configL2ProjectionByMass) :: rconfig    
    real(DP) :: depsAbs, depsRel, dresInit
    
    ! Evaluate the optional arguments as far as possible
    if (present(rL2ProjectionConfig)) rconfig = rL2ProjectionConfig
    ! otherwise use the standard initialisation of that structure!
    
    ! We want to solve the system:
    !
    !     (rvector,phi) = (f,phi)
    !
    ! <=>     M rvector = F
    !
    ! So we have to invert M. As long as the condition number of M is not 
    ! too bad, we don't need tricks like multigrid, but we can use a standard
    ! defect correction!
    !
    ! Clear the initial solution
    call lsyssc_clearVector (rvector)

    ! Assemble a RHS vector using the analytically given function in rtempVector1.
    
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    call linf_buildVectorScalar (&
              rmatrixMass%p_rspatialDiscrTest,rlinform,.true.,&
              rvectorTemp1,fcoeff_buildVectorSc_sim,rcollection)
              
    ! This is also the initial defect because x=0!
    call lsyssc_copyVector (rvectorTemp1,rvectorTemp2)
              
    depsRel = rconfig%depsRel
    depsAbs = rconfig%depsAbs
    
    ! Probably calculate the initial residuum.
    if (rconfig%depsRel .ne. 0.0_DP) then
      dresInit = lsyssc_vectorNorm (rvectorTemp2,rconfig%cnorm)
      if (dresInit .eq. 0.0_DP) dresInit = 1.0_DP
      rconfig%dinitResiduum = dresInit
    else
      ! Set dinitResiduum=1 and dEpsRel = depsAbs, so the relative
      ! and absolute error criterion is identical!
      depsRel = depsAbs
      dresInit = 1.0_DP
    end if

    ! Now, let's start the iteration.
    do iiteration = 1,rconfig%nmaxIterations
    
      ! Multiply by the inverse of the lumped mass matrix:
      ! d := M_l^-1 d
      call lsyssc_invertedDiagMatVec (rmatrixMassLumped,rvectorTemp2,1.0_DP,rvectorTemp2)
    
      ! Add to the main vector:  x = x + omega*d
      call lsyssc_vectorLinearComb (rvectorTemp2,rvector,rconfig%domega,1.0_DP)
      
      ! Set up the defect: d := b-Mx
      call lsyssc_copyVector (rvectorTemp1,rvectorTemp2)
      call lsyssc_scalarMatVec (rmatrixMass, rvector, rvectorTemp2, -1.0_DP, 1.0_DP)
      
      ! Check norms?
      if ((rconfig%depsAbs .ne. 0.0_DP) .or. (rconfig%depsRel .ne. 0.0_DP)) then
      
        rconfig%dabsError = lsyssc_vectorNorm (rvectorTemp2,rconfig%cnorm)
        rconfig%drelError = rconfig%dabsError / dresInit
      
        if (((rconfig%dabsError .le. depsAbs) .or. (depsAbs .eq. 0.0_DP)) .and. &
            (rconfig%dabsError .le. depsRel*dresInit)) then
          ! Quit the loop
          exit
        end if
        
      end if
    
    end do

    ! There is iiterations = niterations+1 if the loop is carried out completely!
    rconfig%iiterations = min(iiteration,rconfig%nmaxIterations)
    
    ! Return the configuration block
    if (present(rL2ProjectionConfig)) rL2ProjectionConfig = rconfig

  end subroutine

end module
