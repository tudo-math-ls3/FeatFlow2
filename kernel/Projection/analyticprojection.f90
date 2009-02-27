!##############################################################################
!# ****************************************************************************
!# <name> analyticprojection </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to calculate projections of analytically
!# given functions, i.e. to calculate the Finite Element representation of
!# an analytic function.
!#
!# One can find the following subroutines here:
!#
!# 1.) anprj_analytL2projectionByMass
!#     -> Performs defect correction with the consistent and lumped mass 
!#        matrix of a FE space to calculate the consistent $L_2$ projection
!#        of an analytically given function.
!#
!# 2.) anprj_discrDirect
!#     -> Evaluates an analytically given function and creates the 
!#        corresponding FE representation. 
!#
!# </purpose>
!##############################################################################

module analyticprojection

  use fsystem
  use linearalgebra
  use linearsystemscalar
  use linearformevaluation
  
  implicit none
  
  private
  
  public :: t_configL2ProjectionByMass
  public :: anprj_analytL2projectionByMass
  public :: anprj_discrDirect
  

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

    ! Maximum iterations to be carried out. Standard = 100.
    ! The iteration stops prematurely if the number of iterations reaches this number.
    integer :: nmaxIterations = 100
    
    ! Type of norm to use for measuring errors. Standard is LINALG_NORML2.
    integer :: cnorm = LINALG_NORML2
    
    ! Type of preconditioner to use.
    ! =0: Use damped Jacobi or lumped mass matrix (depending on whether
    !     rmatrixMassLumped is specified or not). 
    ! =1: Use damped Jacobi.
    ! =2: Use lumped mass matrix. rmatrixMassLumped must be specified 
    !     in the call to l2prj_analytL2projectionByMass.
    ! Standard = 0.
    integer :: cpreconditioner = 0
    
    ! Damping parameter for the iteration. 
    ! If SYS_MAXREAL is specified, the standard damping parameters are used,
    ! which are: = 1.0, if the lumped mass matrix is used for preconditioning,
    !            = 0.7, if the Jacobi preconditioner is used for preconditioning.
    real(DP) :: domega = SYS_MAXREAL
    
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

  subroutine anprj_analytL2projectionByMass (rvector,rmatrixMass,&
      fcoeff_buildVectorSc_sim,rcollection,&
      rL2ProjectionConfig,rmatrixMassLumped,rvectorTemp1,rvectorTemp2)
      
!<description>
  ! Converts an analytically given function fcoeff_buildVectorSc_sim
  ! to a finite element vector rvector by using mass matrices.
  ! So the resulting vector is the consistent $L_2$ projection of the
  ! analytically given function.
  ! The following iteration is performed for this purpose:
  !
  !    $$ x_{n+1}  :=  x_n  +  \omega C^{-1} ( f - M x_n ) $$
  !
  ! with $f$ being a RHS vector generated with fcoeff_buildVectorSc_sim,
  ! $M$ being the consistent mass matrix and $C$ being either the
  ! lumped mass matrix $M_l$ (if specified as rmatrixMassLumped)
  ! or the Jacobi preconditioner. The iteration is performed until the 
  ! error criteria are fulfilled or the given number of steps is reached.
!</description>

!<input>
  ! The consistent mass matrix of the FE space of rvector.
  type(t_matrixScalar), intent(IN) :: rmatrixMass
  
  ! A callback routine for the function to be discretised. The callback routine
  ! has the same syntax as that for evaluating analytic functions for the 
  ! computation of RHS vectors.
  include '../DOFMaintenance/intf_coefficientVectorSc.inc'

  ! OPTIONAL: A pointer to a collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(INOUT), target, optional :: rcollection

  ! OPTIONAL: The lumped mass matrix of the FE space of rvector.
  ! If not specified, the damped Jacobi method will be used.
  type(t_matrixScalar), intent(IN), optional :: rmatrixMassLumped

!</input>

!<inputoutput>
  ! A scalar vector that receives the $L_2$ projection of the function.
  type(t_vectorScalar), intent(INOUT) :: rvector

  ! OPTIONAL: A configuration block for the iteration.
  ! If not specified, the standard settings are used.
  type(t_configL2ProjectionByMass), intent(INOUT), optional :: rL2ProjectionConfig

  ! OPTIONAL: A temporary vector of the same size as rvector.
  type(t_vectorScalar), intent(INOUT), target, optional :: rvectorTemp1

  ! OPTIONAL: A second temporary vector of the same size as rvector.
  type(t_vectorScalar), intent(INOUT), target, optional :: rvectorTemp2
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_linearForm) :: rlinform
    integer :: iiteration
    type(t_configL2ProjectionByMass) :: rconfig    
    real(DP) :: depsAbs, depsRel, dresInit
    type(t_vectorScalar), pointer :: p_rvectorTemp1,p_rvectorTemp2
    real(dp) :: domega
    
    ! Evaluate the optional arguments as far as possible
    if (present(rL2ProjectionConfig)) rconfig = rL2ProjectionConfig
    ! otherwise use the standard initialisation of that structure!
    
    ! Check the parameters.
    if ((rconfig%cpreconditioner .eq. 2) .and. &
        (.not. present(rmatrixMassLumped))) then
      call output_line('Lumped mass matrix not present!',&
          OU_CLASS_ERROR,OU_MODE_STD,'l2prj_analytL2projectionByMass')
      call sys_halt()
    end if
    
    ! Create temporary vectors if necessary.
    if (present(rvectorTemp1)) then
      p_rvectorTemp1 => rvectorTemp1
    else
      allocate(p_rvectorTemp1)
      call lsyssc_duplicateVector (rvector,p_rvectorTemp1,&
          LSYSSC_DUP_COPY,LSYSSC_DUP_EMPTY)
    end if

    if (present(rvectorTemp2)) then
      p_rvectorTemp2 => rvectorTemp2
    else
      allocate(p_rvectorTemp2)
      call lsyssc_duplicateVector (rvector,p_rvectorTemp2,&
          LSYSSC_DUP_COPY,LSYSSC_DUP_EMPTY)
    end if
    
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
              p_rvectorTemp1,fcoeff_buildVectorSc_sim,rcollection)
              
    ! This is also the initial defect because x=0!
    call lsyssc_copyVector (p_rvectorTemp1,p_rvectorTemp2)
              
    ! Get some parameters.
    domega = rconfig%domega
    if (domega .eq. SYS_MAXREAL) then
      if (present(rmatrixMassLumped)) then
        domega = 1.0_DP
      else
        domega = 0.7_DP
      end if
    end if

    depsRel = rconfig%depsRel
    depsAbs = rconfig%depsAbs
    
    ! Probably calculate the initial residuum.
    if (rconfig%depsRel .ne. 0.0_DP) then
      dresInit = lsyssc_vectorNorm (p_rvectorTemp2,rconfig%cnorm)
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
    
      select case (rconfig%cpreconditioner)
      case (0)
        if (present(rmatrixMassLumped)) then
          ! Multiply by the inverse of the lumped mass matrix:
          ! d := M_l^-1 d
          call lsyssc_invertedDiagMatVec (rmatrixMassLumped,p_rvectorTemp2,&
              1.0_DP,p_rvectorTemp2)
        else
          ! Multiply by the inverse of the diagonal: d := D^-1 d
          call lsyssc_invertedDiagMatVec (rmatrixMass,p_rvectorTemp2,&
              1.0_DP,p_rvectorTemp2)
        end if

      case (1)
        ! Multiply by the inverse of the diagonal: d := D^-1 d
        call lsyssc_invertedDiagMatVec (rmatrixMass,p_rvectorTemp2,&
            1.0_DP,p_rvectorTemp2)

      case (2)
        ! Multiply by the inverse of the lumped mass matrix:
        ! d := M_l^-1 d
        call lsyssc_invertedDiagMatVec (rmatrixMassLumped,p_rvectorTemp2,&
            1.0_DP,p_rvectorTemp2)
      end select
    
      ! Add to the main vector:  x = x + omega*d
      call lsyssc_vectorLinearComb (p_rvectorTemp2,rvector,domega,1.0_DP)
      
      ! Set up the defect: d := b-Mx
      call lsyssc_copyVector (p_rvectorTemp1,p_rvectorTemp2)
      call lsyssc_scalarMatVec (rmatrixMass, rvector, p_rvectorTemp2, -1.0_DP, 1.0_DP)
      
      ! Check norms?
      if ((rconfig%depsAbs .ne. 0.0_DP) .or. (rconfig%depsRel .ne. 0.0_DP)) then
      
        rconfig%dabsError = lsyssc_vectorNorm (p_rvectorTemp2,rconfig%cnorm)
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

    ! Release temp vectors.
    if (.not. present(rvectorTemp2)) then
      call lsyssc_releaseVector(p_rvectorTemp2)
      deallocate(p_rvectorTemp2)
    end if

    if (.not. present(rvectorTemp1)) then
      call lsyssc_releaseVector(p_rvectorTemp1)
      deallocate(p_rvectorTemp1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine anprj_discrDirect (rvector,&
      ffunctionReference,rcollection,iorder)
      
!<description>
  ! Converts an analytically given function fcoeff_buildVectorSc_sim
  ! to a finite element vector rvector by direct point evaluation.
  !
  ! Note: This function may not work for all finite elements!
  ! Standard 'immediate' elements (Q0,Q1,..., Q1~) are 
  ! nevertheless supported.
  !
  ! Note: For conformal finite elements (Q0,Q1,...) the result of this
  ! routine is by FE-theory the same as the L2-projection of the function,
  ! but computed much faster.
  ! For nonconformal finite elements, the result is different and also
  ! depending on iorder!
!</description>

!<input>
  ! A callback routine for the function to be discretised. The callback routine
  ! has the same syntax as that for evaluating analytic functions for the 
  ! computation of RHS vectors.
  include '../Postprocessing/intf_refFunctionSc.inc'

  ! OPTIONAL: A pointer to a collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(INOUT), target, optional :: rcollection

  ! OPTIONAL: Order of the evaluation. Standard = 0.
  ! =0: Automatic; use highest possible evaluation (currently =2) to
  !     calculate as exact as possible.
  ! =1: Use pure point evaluation. This is exact for Q0, Q1, Q2,... as well as
  !     for midpoint based Q1~ finite elements.
  ! =2: Use point evaluation for point-value based finite elements 
  !     (Q0, Q1, Q2,...). Use exact integral mean value evaluation for integral
  !     mean value based elements like Q1~ (E030, EM30).
  ! Example: If the given vector is an integral mean value based Q1~
  ! vector and iorder=1 is specified, it will be created as "1st order
  ! approximation" to the given function -- which means by evaluating
  ! in edge midpoints.
  integer, intent(in), optional :: iorder
!</input>

!<inputoutput>
  ! A scalar vector that receives the $L_2$ projection of the function.
  type(t_vectorScalar), intent(INOUT), target :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,k,icurrentElementDistr, ICUBP
    integer :: IEL, IELmax, IELset
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    real(dp), dimension(:), pointer :: p_Ddata
    integer :: h_Dweight, ccub, idof
    real(dp), dimension(:), pointer :: p_Dweight
    integer :: iactualorder
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(CUB_MAXCUBP) :: Domega
    
    ! number of cubature points on the reference element
    integer :: ncubp
    
    ! Number of local degees of freedom for test functions
    integer :: indofTrial
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:), allocatable :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer(PREC_DOFIDX), dimension(:,:), allocatable, target :: IdofsTrial
  
    ! Type of transformation from the reference to the real element 
    integer(I32) :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag
    
    ! Evaluate optional parameters.
    iactualorder = 0
    if (present(iorder)) iactualorder = iorder
    if (iactualorder .eq. 0) iactualorder = 2
    
    ! We choose the midpoint rule for evaluation. Actually, we don't compute
    ! integrals but point values...
    Bder = .false.
    Bder(DER_FUNC) = .true.
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rvector%p_rspatialDiscr%p_rtriangulation
    p_rdiscretisation => rvector%p_rspatialDiscr
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(1000,p_rtriangulation%NEL)
    
    ! Get the data of the FE function we want to initialise.
    ! Clear the FE function in-advance.
    call lsyssc_clearVector (rvector)
    call lsyssc_getbase_double (rvector,p_Ddata)
    
    ! Allocate an array that contains the DOF weight;
    ! Initialise with zero. For 'primal' elements e.g., the entries count
    ! how often a DOF was touched. The contributions are summed up and
    ! later on divided by this value.
    call storage_new1D ('anprj_discrDirectEx31', 'Dweight', size(p_Ddata), &
        ST_DOUBLE, h_Dweight, ST_NEWBLOCK_ZERO)
    call storage_getbase_double(h_Dweight,p_Dweight)

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1,p_rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => p_rdiscretisation%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Now a big element and dimension dependent part: Evaluation points of the
      ! node functionals. The DOF's of most finite elements can be filled by the
      ! correct evaluation of the analytic function -- which is done by
      ! calculating values in different points and creating that functional.
      !
      ! For Q0, Q1, Q2,... we just have to evaluate the vertices, edges,...
      ! and we have the node values.
      !
      ! For Q1~ it depends. If it's the midpoint based version, this is done by
      ! taking the values in the midpoints of the edges. For the integral
      ! mean value based variant, we have to evaluate line integrals.
      ! Exception: If the 'order' is to low, we also take midpoint values
      ! in the integral mean value case of Q1~.
      select case (p_relementDistribution%celement)

      case (EL_E030,EL_EM30)
        if (iactualorder .eq. 1) then
          ! We evaluate in the edge midpoints.
          ! That can be archieved by getting the cubature points from the midpoint rule,
          ! these are exactly the midpoints!
          allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))
          call cub_getCubPoints(CUB_MID, ncubp, Dxi, Domega)

          ! Reformat the cubature points; they are in the wrong shape!
          do i=1,ncubp
            do k=1,ubound(p_DcubPtsRef,1)
              p_DcubPtsRef(k,i) = Dxi(i,k)
            end do
          end do
        else
          ! In the 'full order' approximation, we compute the 2-point Gauss formula
          ! on all 4 lines of the element shape. Get the coordinates of the
          ! cubature points on a line [-1,1].
          call cub_getCubPoints(CUB_G2_1D, ncubp, Dxi, Domega)
          
          ! Transfer the coordinates to the four edges of the reference element
          allocate(p_DcubPtsRef(NDIM2D,CUB_MAXCUBP))
          p_DcubPtsRef(1,1) = Dxi(1,1)
          p_DcubPtsRef(2,1) = -1.0_DP
          p_DcubPtsRef(1,2) = Dxi(2,1)
          p_DcubPtsRef(2,2) = -1.0_DP

          p_DcubPtsRef(1,3) = 1.0_DP
          p_DcubPtsRef(2,3) = Dxi(1,1)
          p_DcubPtsRef(1,4) = 1.0_DP
          p_DcubPtsRef(2,4) = Dxi(2,1)

          p_DcubPtsRef(1,5) = -Dxi(1,1)
          p_DcubPtsRef(2,5) = 1.0_DP
          p_DcubPtsRef(1,6) = -Dxi(2,1)
          p_DcubPtsRef(2,6) = 1.0_DP

          p_DcubPtsRef(1,7) = -1.0_DP
          p_DcubPtsRef(2,7) = -Dxi(1,1)
          p_DcubPtsRef(1,8) = -1.0_DP
          p_DcubPtsRef(2,8) = -Dxi(2,1)
          
          ncubp = 8
        
        end if

      case default
        ! For most of the standard finite elements based on point values,
        ! the evaluation points coincide with the cubature points of that
        ! cubature formula that lump the mass matrix in that FE space.
        ! So try to get them...
        ccub = spdiscr_getLumpCubature (p_relementDistribution%celement)
        if (ccub .eq. 0) then
          ! That FE space does not support lumped mass matrix and so
          ! we have no point coordinates available to get the DOF values :(
          call output_line ('Element not supported!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'anprj_discrDirect')
          call sys_halt()
        else
          ! Get the coordinates of that points.
          ! The weights Domega are ignored in the following...
          allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))
          call cub_getCubPoints(ccub, ncubp, Dxi, Domega)

          ! Reformat the cubature points; they are in the wrong shape!
          do i=1,ncubp
            do k=1,ubound(p_DcubPtsRef,1)
              p_DcubPtsRef(k,i) = Dxi(i,k)
            end do
          end do
        end if

      end select
      
      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the function values
      allocate(Dcoefficients(ncubp,nelementsPerBlock))
    
      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! We don't want to evaluate the element, but we need some information
      ! from the preparation routine for the callback routine.
      ! So we 'pretend' that we want to evaluate.
      !
      ! Create an element evaluation tag that computes us the element corners
      ! (we may need them for integration), the coordinates on the
      ! reference and on the real elements. Jacobian mapping, Jacobian determinants
      ! etc. are not needed since we don't evaluate...
      cevaluationTag = EL_EVLTAG_COORDS + EL_EVLTAG_REFPOINTS + EL_EVLTAG_REALPOINTS
                      
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL
    
      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, 1000
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        IELmax = min(NEL,IELset-1+1000)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
                                     
        ! Prepare the call to the evaluation routine of the analytic function.    
        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => IdofsTrial
        rintSubset%celement = p_relementDistribution%celement
    
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! It's time to call our coefficient function to calculate the
        ! function values in the cubature points:  u(x,y)
        call ffunctionReference (DER_FUNC,p_rdiscretisation, &
                    int(IELmax-IELset+1),ncubp,&
                    revalElementSet%p_DpointsReal,&
                    IdofsTrial,rintSubset,&
                    Dcoefficients(:,1:IELmax-IELset+1_I32),rcollection)

        ! Another element dependent part: evaluation of the functional.
        select case (p_relementDistribution%celement)
        
        case (EL_E030,EL_EM30)
        
          if (iactualorder .eq. 1) then
          
            ! Loop through elements in the set and for each element,
            do IEL=1,IELmax-IELset+1
            
              do icubp = 1, ncubp
              
                ! Sum up the calculated value to the existing value of the 
                ! corresponding DOF.
                p_Ddata(IdofsTrial(icubp,IEL)) = p_Ddata(IdofsTrial(icubp,IEL)) + &
                  Dcoefficients(icubp,IEL)
                  
                ! Count the entry
                p_Dweight(IdofsTrial(icubp,IEL)) = p_Dweight(IdofsTrial(icubp,IEL)) + 1.0_DP

              end do ! ICUBP 

            end do ! IEL
            
          else
          
            ! Loop through elements in the set and for each element,
            do IEL=1,IELmax-IELset+1
            
              ! Loop through the DOF's on the current element.
              do idof = 1,indofTrial 
              
                ! Calculate the DOF. For that purpose, calculate the line
                ! integral (using the cubature points calculated above).
                ! Weight by 0.5 to get the integral corresponding to an interval
                ! of length 1 instead of 2 (what the cubature formula
                ! uses: [-1,1]). That's already the integral mean
                ! value we save as DOF...
                p_Ddata(IdofsTrial(idof,IEL)) = p_Ddata(IdofsTrial(idof,IEL)) + &
                  0.5_DP * (Dcoefficients(2*(idof-1)+1,IEL) * Domega(1) + &
                            Dcoefficients(2*(idof-1)+2,IEL) * Domega(2))
              
                ! Count the entry
                p_Dweight(IdofsTrial(idof,IEL)) = p_Dweight(IdofsTrial(idof,IEL)) + 1.0_DP

              end do

            end do ! IEL

          end if

        case default
        
          ! Loop through elements in the set and for each element,
          do IEL=1,IELmax-IELset+1
          
            do icubp = 1, ncubp
            
              ! Sum up the calculated value to the existing value of the 
              ! corresponding DOF.
              p_Ddata(IdofsTrial(icubp,IEL)) = p_Ddata(IdofsTrial(icubp,IEL)) + &
                Dcoefficients(icubp,IEL)
                
              ! Count the entry
              p_Dweight(IdofsTrial(icubp,IEL)) = p_Dweight(IdofsTrial(icubp,IEL)) + 1.0_DP

            end do ! ICUBP 

          end do ! IEL
          

        end select        
        
        ! Release the temporary domain integration structure again
        call domint_doneIntegration (rintSubset)
    
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)
      
      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)

    end do ! icurrentElementDistr
    
    ! Take the mean value in all entries. We just summed up all contributions and now
    ! this divides by the number of contributions...
    ! All DOF's should be touched, so we assume that there is Dweight != 0 everywhere.
    do i=1,size(p_Ddata)
      p_Ddata(i) = p_Ddata(i) / p_Dweight(i)
    end do
    
    ! Release temp data
    call storage_free (h_Dweight)

  end subroutine

end module
