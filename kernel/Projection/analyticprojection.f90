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
!#        matrix of a FE space to calculate the consistent <tex>$ L_2 $</tex>
!#        projection of an analytically given function.
!#
!# 2.) anprj_analytRewProjection
!#     -> Performs a restricted element-wise L2-projection of an analytically
!#        given function.
!#
!# 3.) anprj_discrDirect
!#     -> Evaluates an analytically given function and creates the
!#        corresponding FE representation.
!#
!# 4.) anprj_charFctRealBdComp
!#     -> Calculate the characteristic function of a boundary region
!#
!# 5.) anprj_analytL2projectionConstr
!#     -> Performs an L2-projection using the lumped mass matrix and applies
!#        constrained mass antidiffusive to improve the resolution
!# </purpose>
!##############################################################################

module analyticprojection

!$ use omp_lib
  use afcstabbase
  use afcstabscalar
  use afcstabscalarfct
  use basicgeometry
  use boundary
  use collection
  use cubature
  use derivatives
  use dofmapping
  use domainintegration
  use element
  use elementpreprocessing
  use fsystem
  use genoutput
  use groupfemscalar
  use linearalgebra
  use linearformevaluation
  use linearsystemscalar
  use perfconfig
  use scalarpde
  use spatialdiscretisation
  use storage
  use transformation
  use triangulation

  implicit none

  private

  public :: t_configL2ProjectionByMass
  public :: anprj_analytL2projectionByMass
  public :: anprj_analytL2projectionConstr
  public :: anprj_analytRewProjection
  public :: anprj_discrDirect
  public :: anprj_charFctRealBdComp
  
!<constants>

!<constantblock description="averaging specifiers for rew-projection">

  ! use arithmetic average for rew-projection
  integer, parameter, public :: ANPRJ_REW_ARITHMETIC = 0
  
  ! use volume average for rew-projection
  integer, parameter, public :: ANPRJ_REW_VOLUME     = 1

!</constantblock>

!</constants>

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
    ! If SYS_MAXREAL_DP is specified, the standard damping parameters are used,
    ! which are: = 1.0, if the lumped mass matrix is used for preconditioning,
    !            = 0.7, if the Jacobi preconditioner is used for preconditioning.
    real(DP) :: domega = SYS_MAXREAL_DP

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

  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: anprj_perfconfig

  !************************************************************************

contains

  ! ***************************************************************************

!<subroutine>

  subroutine anprj_analytL2projectionByMass (rvector, rmatrixMass,&
      fcoeff_buildVectorSc_sim, rcollection,&
      rL2ProjectionConfig, rmatrixMassLumped, rvectorTemp1, rvectorTemp2)

!<description>
  ! Converts an analytically given function fcoeff_buildVectorSc_sim
  ! to a finite element vector rvector by using mass matrices.
  ! So the resulting vector is the consistent $L_2$ projection of the
  ! analytically given function.
  ! The following iteration is performed for this purpose:
  !
  !    <tex> $$ x_{n+1}  :=  x_n  +  \omega C^{-1} ( f - M x_n ) $$ </tex>
  !
  ! with $f$ being a RHS vector generated with fcoeff_buildVectorSc_sim,
  ! $M$ being the consistent mass matrix and $C$ being either the
  ! lumped mass matrix $M_l$ (if specified as rmatrixMassLumped)
  ! or the Jacobi preconditioner. The iteration is performed until the
  ! error criteria are fulfilled or the given number of steps is reached.
!</description>

!<input>
  ! The consistent mass matrix of the FE space of rvector.
  type(t_matrixScalar), intent(in) :: rmatrixMass

  ! A callback routine for the function to be discretised. The callback routine
  ! has the same syntax as that for evaluating analytic functions for the
  ! computation of RHS vectors.
  include '../DOFMaintenance/intf_coefficientVectorSc.inc'

  ! OPTIONAL: A pointer to a collection structure. This structure is
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! OPTIONAL: The lumped mass matrix of the FE space of rvector.
  ! If not specified, the damped Jacobi method will be used.
  type(t_matrixScalar), intent(in), optional :: rmatrixMassLumped
!</input>

!<inputoutput>
  ! A scalar vector that receives the $L_2$ projection of the function.
  type(t_vectorScalar), intent(inout) :: rvector

  ! OPTIONAL: A configuration block for the iteration.
  ! If not specified, the standard settings are used.
  type(t_configL2ProjectionByMass), intent(inout), optional :: rL2ProjectionConfig

  ! OPTIONAL: A temporary vector of the same size as rvector.
  type(t_vectorScalar), intent(inout), target, optional :: rvectorTemp1

  ! OPTIONAL: A second temporary vector of the same size as rvector.
  type(t_vectorScalar), intent(inout), target, optional :: rvectorTemp2
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_linearForm) :: rlinform
    integer :: iiteration
    type(t_configL2ProjectionByMass) :: rconfig
    real(DP) :: depsAbs, depsRel, dresInit
    type(t_vectorScalar), pointer :: p_rvectorTemp1,p_rvectorTemp2
    real(dp) :: domega
    type(t_scalarCubatureInfo) :: rcubatureInfo

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
    ! too bad, we do not need tricks like multigrid, but we can use a standard
    ! defect correction!
    !
    ! Clear the initial solution
    call lsyssc_clearVector (rvector)

    ! Assemble a RHS vector using the analytically given function in rtempVector1.
    call spdiscr_createDefCubStructure (rvector%p_rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO)

    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    call linf_buildVectorScalar (rlinform,.true.,&
        p_rvectorTemp1,rcubatureInfo,fcoeff_buildVectorSc_sim,rcollection)

    call spdiscr_releaseCubStructure (rcubatureInfo)

    ! This is also the initial defect because x=0!
    call lsyssc_copyVector (p_rvectorTemp1,p_rvectorTemp2)

    ! Get some parameters.
    domega = rconfig%domega
    if (domega .eq. SYS_MAXREAL_DP) then
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

    ! Now, let us start the iteration.
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
      call lsyssc_matVec (rmatrixMass, rvector, p_rvectorTemp2, -1.0_DP, 1.0_DP)

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

  subroutine anprj_analytRewProjection(rvector, fcoeff_buildVectorSc_sim, &
                                       rcubature, cavrgType, rcollection, rperfconfig)

!<description>
  ! Projects an analytically given function fcoeff_buildVectorSc_sim
  ! to a finite element vector rvector by using restricted element-wise
  ! L2-projection.
!</description>

!<input>
  ! A callback routine for the function to be discretised. The callback routine
  ! has the same syntax as that for evaluating analytic functions for the
  ! computation of RHS vectors.
  include '../DOFMaintenance/intf_coefficientVectorSc.inc'

  ! A cubature info structure describing the cubature rule to be used
  ! for integration.
  type(t_scalarCubatureInfo), intent(in) :: rcubature
  
  ! OPTIONAL: An averaging specifier, one of the ANPRJ_REW_XXX constants.
  ! If not given, ANPRJ_REW_VOLUME is used.
  integer, optional, intent(in) :: cavrgType
  
  ! OPTIONAL: A pointer to a collection structure. This structure is
  ! given to the callback function.
  type(t_collection), optional, intent(inout) :: rcollection

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
  ! A scalar vector that receives the rew-projection of the function.
  type(t_vectorScalar), target, intent(inout) :: rvector
!</output>

  ! local variables
  type(t_linearForm) :: rform
  type(t_vectorScalar) :: rvecWeight
  type(t_spatialDiscretisation), pointer :: p_rdisc
  type(t_elementDistribution), pointer :: p_relDist
  real(DP), dimension(:), pointer :: p_Dx, p_Dw, Domega, Dval
  real(DP), dimension(:,:), pointer :: DcubPtsRef, p_Ddetj, Dmass, Dvec
  real(DP), dimension(:,:,:), pointer :: Dcoeff
  real(DP), dimension(:,:,:,:), pointer :: Dbas
  integer, dimension(:), pointer :: p_IelemList, Ipivot
  integer, dimension(:,:), pointer :: Idofs
  type(t_domainIntSubset) :: rintSubset
  type(t_evalElementSet) :: reval
  real(DP) :: dom, daux, dvol
  integer :: ctype, ielDist, ndof, ncubp, nelTodo, nelDone, nelSize, &
      iel, icpt,i,j, info, iblockSize
  integer(I32) :: celement, ctrafo, cevalTag, ccubature
  logical, dimension(EL_MAXNDER) :: Bder
    
    ! choose averaging type
    ctype = ANPRJ_REW_VOLUME
    if(present(cavrgType)) ctype = cavrgType
    
    ! choose a block size
    if(present(rperfconfig)) then
      iblockSize = rperfconfig%NELEMSIM
    else
      iblockSize = anprj_perfconfig%NELEMSIM
    end if
    
    ! set up a linearform descriptor
    Bder = .false.
    Bder(DER_FUNC) = .true.
    rform%itermCount = 1
    rform%Idescriptors(1) = DER_FUNC
    
    ! fetch the discretisation
    p_rdisc => rvector%p_rspatialDiscr
    
    ! allocate a weight vector
    call lsyssc_createVector(rvecWeight, rvector%NEQ, .true.)
    call lsyssc_getbase_double(rvector, p_Dx)
    call lsyssc_getbase_double(rvecWeight, p_Dw)
    
    ! loop over all element distributions
    do ielDist = 1, p_rdisc%inumFESpaces
    
      ! fetch the element distribution
      p_relDist => p_rdisc%RelementDistr(ielDist)
      
      if(p_relDist%NEL .eq. 0) cycle
      nelTodo = 1
      
      ! fetch the element list
      call storage_getbase_int(p_relDist%h_IelementList, p_IelemList)
      
      ! fetch the element and trafo type
      celement = p_relDist%celement
      ccubature = rcubature%p_RinfoBlocks(ielDist)%ccubature
      ctrafo = elem_igetTrafoType(celement)
      cevalTag = ior(elem_getEvaluationTag(celement), EL_EVLTAG_REALPOINTS)
      
      ! fetch the number of dofs and cubature points
      ndof = elem_igetNDofLoc(celement)
      ncubp = cub_igetNumPts(ccubature)
      
      ! allocate two arrays for the cubature formula
      allocate(Domega(ncubp))
      allocate(DcubPtsRef(trafo_igetReferenceDimension(ctrafo),ncubp))

      ! Get the cubature formula
      call cub_getCubature(ccubature,DcubPtsRef, Domega)
      
      ! allocate basis function array
      allocate(Dbas(ndof,elem_getMaxDerivative(celement), ncubp, iblockSize))

      ! Allocate memory for the DOF`s of all the elements.
      allocate(Idofs(ndof,iblockSize))

      ! Allocate memory for the coefficients.
      allocate(Dcoeff(1,ncubp,iblockSize))
      
      allocate(Dmass(ndof, ndof))
      allocate(Dvec(ndof,iblockSize))
      allocate(Ipivot(ndof))
      allocate(Dval(ndof))
      
      ! initialise weighting values for arithmetic average
      if(ctype .eq. ANPRJ_REW_ARITHMETIC) &
        Dval = 1.0_DP
  
      ! Initialisation of the element set.
      call elprep_init(reval)
      
      ! loop over all element blocks
      do nelDone = 0, p_relDist%NEL, iblockSize
      
        ! compute the number of elements to be processed
        nelTodo = min(p_relDist%NEL, neldone+iblockSize)
        nelSize = nelTodo - nelDone
        
        ! perform dof-mapping
        call dof_locGlobMapping_mult(p_rdisc, p_IelemList(nelDone+1:nelTodo), Idofs)
        
        ! prepare element for evaluation
        call elprep_prepareSetForEvaluation (reval, cevalTag, p_rdisc%p_rtriangulation,&
            p_IelemList(nelDone+1:nelTodo), ctrafo, DcubPtsRef)
        p_Ddetj => reval%p_Ddetj
        
        ! prepare domain integration
        call domint_initIntegrationByEvalSet (reval,rintSubset)
        rintSubset%ielementStartIdx     =  nelDone+1
        rintSubset%p_Ielements          => p_IelemList(nelDone:nelTodo)
        rintSubset%p_IdofsTrial         => Idofs
        rintSubset%celement             =  celement
      
        ! evaluate callback function
        call fcoeff_buildVectorSc_sim (p_rdisc,rform, &
            nelTodo-nelDone+1,ncubp,reval%p_DpointsReal(:,:,1:nelSize), &
            Idofs,rintSubset, Dcoeff(:,:,1:nelSize), rcollection)

        ! release domain integration
        call domint_doneIntegration(rintSubset)
        
        ! Calculate the values of the basis functions.
        call elem_generic_sim2 (celement, reval, Bder, Dbas)
        
        ! loop over all elements in the current batch
        do iel = 1, nelSize
        
          ! clear local vector, mass matrix and element volume
          Dvec(:,iel) = 0.0_DP
          Dmass = 0.0_DP
          dvol = 0.0_DP
          
          ! loop over all cubature points
          do icpt = 1, ncubp
          
            ! compute integration weight
            dom = Domega(icpt) * abs(p_Ddetj(icpt, iel))
            
            ! get the function value pre-multiplied by integration weight
            daux = dom * Dcoeff(1, icpt, iel)
            
            ! update element volume
            dvol = dvol + dom
            
            ! loop over all local dofs and build the local mass matrix and rhs vector
            do i = 1, ndof
              ! build rhs vector entry
              Dvec(i,iel) = Dvec(i,iel) + daux * Dbas(i,DER_FUNC,icpt,iel)
              do j = 1, ndof
                ! build mass matrix entry
                Dmass(i,j) = Dmass(i,j) + dom * Dbas(i,DER_FUNC,icpt,iel)*Dbas(j,DER_FUNC,icpt,iel)
              end do
            end do
          end do
          
          ! solve the local system
          call DGESV(ndof,1, Dmass,ndof, Ipivot, Dvec(:,iel), ndof, info)
          
          ! choose weights
          if(ctype .eq. ANPRJ_REW_VOLUME) &
            Dval = dvol
          
          ! incorporate local solution into global vector
          do i = 1, ndof
            j = Idofs(i, iel)
            ! update solution
            p_Dx(j) = p_Dx(j) + Dval(i)*Dvec(i,iel)
            ! update weights
            p_Dw(j) = p_Dw(j) + Dval(i)
          end do
        
        end do ! iel
        
      end do
      
      deallocate(Ipivot)
      deallocate(Dvec)
      deallocate(Dmass)
      deallocate(Dcoeff)
      deallocate(Idofs)
      deallocate(Dbas)
      deallocate(DcubPtsRef)
      deallocate(Domega)
      
    end do
    
    ! scale solution by sumed weights
    do i = 1, rvector%NEQ
      p_Dx(i) = p_Dx(i) / p_Dw(i)
    end do
    
    ! release weight vector
    call lsyssc_releaseVector(rvecWeight)
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine anprj_discrDirect (rvector,&
      ffunctionReference, rcollection, iorder, ntempArrays, rperfconfig)

!<description>
  ! Converts an analytically given function ffunctionReference
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
  ! A callback routine for the function to be discretised.
  include '../Postprocessing/intf_refFunctionSc.inc'

  ! OPTIONAL: A pointer to a collection structure. This structure is
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection

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

  ! OPTIONAL: Number of temp arrays.
  ! If this is specified (and larger than zero), a number of temporary
  ! arrays is allocated and provided to the callback routine.
  ! This temporary memory is user defined and can be used during the
  ! assembly, e.g., for the temporary evaluation of FEM functions.
  !
  ! The temporary memory is available in the callback function as
  !    rdomainIntSubset%p_DtempArrays !
  integer, intent(in), optional :: ntempArrays

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! A scalar vector that receives the $L_2$ projection of the function.
  type(t_vectorScalar), intent(inout), target :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,k,ielemGroup, ICUBP
    integer :: IEL, IELmax, IELset
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    real(dp), dimension(:), pointer :: p_Ddata
    integer :: h_Dweight, idof
    integer(I32) :: ccub
    real(dp), dimension(:), pointer :: p_Dweight
    integer :: iactualorder
    integer(I32) :: celement

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

    ! Number of elements in the current element group
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:), allocatable :: Dcoefficients

    ! Number of elements in a block. Normally =NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock

    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet

    ! An allocateable array accepting the DOF`s of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial

    ! Type of transformation from the reference to the real element
    integer(I32) :: ctrafoType

    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! Pointer to temporary memory for callback functions.
    real(DP), dimension(:,:,:), pointer :: p_DtempArrays

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => anprj_perfconfig
    end if

    ! Evaluate optional parameters.
    iactualorder = 0
    if (present(iorder)) iactualorder = iorder
    if (iactualorder .eq. 0) iactualorder = 2

    ! We choose the midpoint rule for evaluation. Actually, we do not compute
    ! integrals but point values...
    Bder = .false.
    Bder(DER_FUNC) = .true.

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rvector%p_rspatialDiscr%p_rtriangulation
    p_rdiscretisation => rvector%p_rspatialDiscr

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(p_rperfconfig%NELEMSIM,p_rtriangulation%NEL)

    ! Get the data of the FE function we want to initialise.
    ! Clear the FE function in-advance.
    call lsyssc_clearVector (rvector)
    call lsyssc_getbase_double (rvector,p_Ddata)

    ! Allocate an array that contains the DOF weight;
    ! Initialise with zero. For 'primal' elements e.g., the entries count
    ! how often a DOF was touched. The contributions are summed up and
    ! later on divided by this value.
    call storage_new ('anprj_discrDirectEx31', 'Dweight', size(p_Ddata), &
        ST_DOUBLE, h_Dweight, ST_NEWBLOCK_ZERO)
    call storage_getbase_double(h_Dweight,p_Dweight)

    ! Now loop over the different element groups (=combinations
    ! of trial and test functions) in the discretisation.

    do ielemGroup = 1,spdiscr_getNelemGroups(p_rdiscretisation)

      ! Activate the current element group
      !
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call spdiscr_getElemGroupInfo (p_rdiscretisation,ielemGroup,celement,NEL,p_IelementList)

      ! Cancel if this element group is empty.
      if (NEL .eq. 0) cycle

      ! Get the number of local DOF`s for trial functions
      indofTrial = elem_igetNDofLoc(celement)

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(celement)

      ! Now a big element and dimension dependent part: Evaluation points of the
      ! node functionals. The DOF`s of most finite elements can be filled by the
      ! correct evaluation of the analytic function -- which is done by
      ! calculating values in different points and creating that functional.
      !
      ! For Q0, Q1, Q2,... we just have to evaluate the vertices, edges,...
      ! and we have the node values.
      !
      ! For Q1~ it depends. If it is the midpoint based version, this is done by
      ! taking the values in the midpoints of the edges. For the integral
      ! mean value based variant, we have to evaluate line integrals.
      ! Exception: If the 'order' is to low, we also take midpoint values
      ! in the integral mean value case of Q1~.
      select case (celement)

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

      case (EL_Q2)
        ! Get the coordinates of that points.
        ! Here, we have to follow the order of the local DOFs of Q2!
        allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

        p_DcubPtsRef(1,1) = -1.0_DP
        p_DcubPtsRef(2,1) = -1.0_DP

        p_DcubPtsRef(1,2) =  1.0_DP
        p_DcubPtsRef(2,2) = -1.0_DP

        p_DcubPtsRef(1,3) =  1.0_DP
        p_DcubPtsRef(2,3) =  1.0_DP

        p_DcubPtsRef(1,4) = -1.0_DP
        p_DcubPtsRef(2,4) =  1.0_DP

        p_DcubPtsRef(1,5) =  0.0_DP
        p_DcubPtsRef(2,5) = -1.0_DP

        p_DcubPtsRef(1,6) =  1.0_DP
        p_DcubPtsRef(2,6) =  0.0_DP

        p_DcubPtsRef(1,7) =  0.0_DP
        p_DcubPtsRef(2,7) =  1.0_DP

        p_DcubPtsRef(1,8) = -1.0_DP
        p_DcubPtsRef(2,8) =  0.0_DP

        p_DcubPtsRef(1,9) =  0.0_DP
        p_DcubPtsRef(2,9) =  0.0_DP

        ncubp = 9
        
      case (EL_P2)
        ! Get the coordinates of that points.
        ! Here, we have to follow the order of the local DOFs of P2!
        allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

        ! vertices
        p_DcubPtsRef(1,1) = 1.0_DP
        p_DcubPtsRef(2,1) = 0.0_DP
        p_DcubPtsRef(3,1) = 0.0_DP
        p_DcubPtsRef(1,2) = 0.0_DP
        p_DcubPtsRef(2,2) = 1.0_DP
        p_DcubPtsRef(3,2) = 0.0_DP
        p_DcubPtsRef(1,3) = 0.0_DP
        p_DcubPtsRef(2,3) = 0.0_DP
        p_DcubPtsRef(3,3) = 1.0_DP

        ! edge midpoints
        p_DcubPtsRef(1,4) = 0.5_DP
        p_DcubPtsRef(2,4) = 0.5_DP
        p_DcubPtsRef(3,4) = 0.0_DP
        p_DcubPtsRef(1,5) = 0.0_DP
        p_DcubPtsRef(2,5) = 0.5_DP
        p_DcubPtsRef(3,5) = 0.5_DP
        p_DcubPtsRef(1,6) = 0.5_DP
        p_DcubPtsRef(2,6) = 0.0_DP
        p_DcubPtsRef(3,6) = 0.5_DP

        ncubp = 6

      case default
        ! For most of the standard finite elements based on point values,
        ! the evaluation points coincide with the cubature points of that
        ! cubature formula that lump the mass matrix in that FE space.
        ! So try to get them...
        ccub = spdiscr_getLumpCubature (celement)
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

      ! Allocate memory for the DOF`s of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the function values
      allocate(Dcoefficients(ncubp,nelementsPerBlock))

      ! Allocate memory for user defined data.
      nullify(p_DtempArrays)
      if (present(ntempArrays)) then
        allocate(p_DtempArrays(ncubp,nelementsPerBlock,ntempArrays))
      end if

      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! We do not want to evaluate the element, but we need some information
      ! from the preparation routine for the callback routine.
      ! So we 'pretend' that we want to evaluate.
      !
      ! Create an element evaluation tag that computes us the element corners
      ! (we may need them for integration), the coordinates on the
      ! reference and on the real elements. Jacobian mapping, Jacobian determinants
      ! etc. are not needed since we do not evaluate...
      cevaluationTag = EL_EVLTAG_COORDS + EL_EVLTAG_REFPOINTS + EL_EVLTAG_REALPOINTS

      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, p_rperfconfig%NELEMSIM

        ! We always handle NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most NELEMSIM
        ! elements simultaneously.
        IELmax = min(NEL,IELset-1+p_rperfconfig%NELEMSIM)

        ! Calculate the global DOF`s into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF`s of our NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)

        ! In the next loop, we do not have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! Prepare the call to the evaluation routine of the analytic function.
        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)

        ! Associate temp memory for the callback.
        if (associated(p_DtempArrays)) then
          call domint_setTempMemory (rintSubset,p_DtempArrays)
        end if

        !rintSubset%ielemGroupibution = ielemGroup
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => IdofsTrial
        rintSubset%celement = celement

        ! It is time to call our coefficient function to calculate the
        ! function values in the cubature points:  u(x,y)
        call ffunctionReference (DER_FUNC,p_rdiscretisation, &
                    int(IELmax-IELset+1),ncubp,&
                    revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
                    IdofsTrial,rintSubset,&
                    Dcoefficients(:,1:IELmax-IELset+1),rcollection)

        ! Another element dependent part: evaluation of the functional.
        select case (celement)

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

              ! Loop through the DOF`s on the current element.
              do idof = 1,indofTrial

                ! Calculate the DOF. For that purpose, calculate the line
                ! integral (using the cubature points calculated above).
                ! Weight by 0.5 to get the integral corresponding to an interval
                ! of length 1 instead of 2 (what the cubature formula
                ! uses: [-1,1]). That is already the integral mean
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
      
      ! Release the temp memory
      if (associated(p_DtempArrays)) then
        deallocate(p_DtempArrays)
      end if

    end do ! ielemGroup

    ! Take the mean value in all entries. We just summed up all contributions and now
    ! this divides by the number of contributions...
    ! All DOF`s should be touched, so we assume that there is Dweight != 0 everywhere.
    do i=1,size(p_Ddata)
      p_Ddata(i) = p_Ddata(i) / p_Dweight(i)
    end do

    ! Release temp data
    call storage_free (h_Dweight)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine anprj_charFctRealBdComp (rboundaryRegion, rvector)

!<description>
  ! Calculates the characteristic function of a boundary region.
  ! rvector is a FE solution vector, discretised by P1,Q1,P2,Q2 or Ex3x,
  ! rboundaryRegion a boundary region on the boundary of the underlying
  ! domain. The routine will now set the DOF's of all vertices to 1 which
  ! belong to this boundary region. All other DOF's are ignored, so the
  ! routine can be called for multiple boundary regions to accumulate
  ! a characteristic vector.
!</description>

!<input>
  ! A boundary region representing a part of the boudary.
  type(t_boundaryRegion), intent(in), target :: rboundaryRegion
!</input>

!<inputoutput>
  ! A vector that should receive the characteristic function of that
  ! boundary region in the underlying FEM-space.
  ! Only point-based FEM spaces are allowed.
  ! Ontegral mean value based FEM-spaces (Q1~) will be treated as point based.
  type(t_vectorScalar), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    integer(I32) :: celement, celement2
    
    ! This task is a bit tricky and depends on the discretisation.
    ! We have to bypass DOFMapping!
    !
    ! We only support a special set of discretisations, it get`s too
    ! complicated otherwise...
    if (spdiscr_getNelemGroups(rvector%p_rspatialDiscr) .eq. 1) then

      call spdiscr_getElemGroupInfo (rvector%p_rspatialDiscr,1,celement)

      if ((celement .eq. EL_P1_2D) .or. &
          (celement .eq. EL_Q1_2D)) then
        ! P1 or Q1
        call charFctRealBdComp2d_P1Q1 (rboundaryRegion,rvector)
        return

      else if ((celement .eq. EL_P2_2D) .or. &
               (celement .eq. EL_Q2_2D)) then
        ! P2 or Q2
        call charFctRealBdComp2d_P2Q2 (rboundaryRegion,rvector)
        return

      else if ((celement .eq. EL_P3_2D) .or. &
               (celement .eq. EL_Q3_2D)) then
        ! P3 or Q3
        call charFctRealBdComp2d_P3Q3 (rboundaryRegion,rvector)
        return

      else if ((elem_getPrimaryElement(celement) .eq. EL_Q1T_2D) .or. &
               (elem_getPrimaryElement(celement) .eq. EL_Q1TB_2D)) then
        ! Q1~ or Q1~bubble
        call charFctRealBdComp2d_Ex3x (rboundaryRegion,rvector)
        return

      end if

    else if (spdiscr_getNelemGroups(rvector%p_rspatialDiscr) .eq. 2) then
    
      call spdiscr_getElemGroupInfo (rvector%p_rspatialDiscr,1,celement)
      call spdiscr_getElemGroupInfo (rvector%p_rspatialDiscr,2,celement2)

      if ((celement .eq. EL_P1_2D) .and. &
          (celement2 .eq. EL_Q1_2D)) then
        ! Mixed P1/Q1
        call charFctRealBdComp2d_P1Q1 (rboundaryRegion,rvector)
        return
      end if
    end if

    call output_line ('Discretisation not supported!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'anprj_charFctRealBdComp')
    call sys_halt()

  end subroutine

  ! ---------------------------------------------------------------

  subroutine charFctRealBdComp2d_P1Q1 (rboundaryRegion,rvector)
    ! Worker routine for P1/Q1 space.
    type(t_boundaryRegion), intent(in), target :: rboundaryRegion
    type(t_vectorScalar), intent(inout), target :: rvector

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: p_IverticesAtBoundary,p_IboundaryCpIdx
    real(dp), dimension(:), pointer :: p_DvertexParameterValue
    integer :: ibct,ivbd
    real(DP), dimension(:), pointer :: p_Ddata

    ! Get some data
    p_rtriangulation => rvector%p_rspatialDiscr%p_rtriangulation
    call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
    call storage_getbase_int(p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_double(p_rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)
    call lsyssc_getbase_double (rvector,p_Ddata)

    ! Loop through the boundary components.
    do ibct = 1,p_rtriangulation%nbct
      ! Loop through the vertices
      do ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-1
        ! Check if the vertex is in the region
        if (boundary_isInRegion (rboundaryRegion,ibct,&
            p_DvertexParameterValue(ivbd))) then
          ! Add a 1 to the vector
          p_Ddata(p_IverticesAtBoundary(ivbd)) = &
              p_Ddata(p_IverticesAtBoundary(ivbd)) + 1.0_DP
        end if
      end do
    end do

  end subroutine

  ! ---------------------------------------------------------------

  subroutine charFctRealBdComp2d_P2Q2(rboundaryRegion,rvector)
    ! Worker routine for P1/Q1 space.
    type(t_boundaryRegion), intent(in), target :: rboundaryRegion
    type(t_vectorScalar), intent(inout), target :: rvector

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: p_IverticesAtBoundary,p_IboundaryCpIdx
    real(dp), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    real(dp), dimension(:), pointer :: p_DedgeParameterValue
    integer :: ibct,ivbd,iebd
    real(DP), dimension(:), pointer :: p_Ddata

    ! Get some data
    p_rtriangulation => rvector%p_rspatialDiscr%p_rtriangulation
    call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
    call storage_getbase_int(p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_double(p_rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)
    call storage_getbase_int(p_rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    call storage_getbase_int(p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_double(p_rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)
    call lsyssc_getbase_double (rvector,p_Ddata)

    ! Loop through the boundary components.
    do ibct = 1,p_rtriangulation%nbct
      ! Loop through the vertices
      do ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-1
        ! Check if the vertex is in the region
        if (boundary_isInRegion (rboundaryRegion,ibct,&
            p_DvertexParameterValue(ivbd))) then
          ! Add a 1 to the vector
          p_Ddata(p_IverticesAtBoundary(ivbd)) = &
              p_Ddata(p_IverticesAtBoundary(ivbd)) + 1.0_DP
        end if
      end do

      ! Loop through the edges
      do iebd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-1
        ! Check if the edge is in the region
        if (boundary_isInRegion (rboundaryRegion,ibct,&
            p_DedgeParameterValue(iebd))) then
          ! Add a 1 to the vector
          p_Ddata(p_IedgesAtBoundary(iebd)+p_rtriangulation%NVT) = &
              p_Ddata(p_IedgesAtBoundary(iebd)+p_rtriangulation%NVT) + 1.0_DP
        end if
      end do
    end do

  end subroutine

  ! ---------------------------------------------------------------

  subroutine charFctRealBdComp2d_P3Q3(rboundaryRegion,rvector)
    ! Worker routine for P1/Q1 space.
    type(t_boundaryRegion), intent(in), target :: rboundaryRegion
    type(t_vectorScalar), intent(inout), target :: rvector

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: p_IverticesAtBoundary,p_IboundaryCpIdx
    real(dp), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    real(dp), dimension(:), pointer :: p_DedgeParameterValue
    integer :: ibct,ivbd,iebd
    real(DP), dimension(:), pointer :: p_Ddata

    ! Get some data
    p_rtriangulation => rvector%p_rspatialDiscr%p_rtriangulation
    call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
    call storage_getbase_int(p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_double(p_rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)
    call storage_getbase_int(p_rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    call storage_getbase_int(p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_double(p_rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)
    call lsyssc_getbase_double (rvector,p_Ddata)

    ! Loop through the boundary components.
    do ibct = 1,p_rtriangulation%nbct
      ! Loop through the vertices
      do ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-1
        ! Check if the vertex is in the region
        if (boundary_isInRegion (rboundaryRegion,ibct,&
            p_DvertexParameterValue(ivbd))) then
          ! Add a 1 to the vector
          p_Ddata(p_IverticesAtBoundary(ivbd)) = &
              p_Ddata(p_IverticesAtBoundary(ivbd)) + 1.0_DP
        end if
      end do

      ! Loop through the edges
      do iebd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-1
        ! Check if the edge is in the region
        if (boundary_isInRegion (rboundaryRegion,ibct,&
            p_DedgeParameterValue(iebd))) then
          ! Add a 1 to the vector
          p_Ddata(2*p_IedgesAtBoundary(iebd)-1+p_rtriangulation%NVT) = &
              p_Ddata(2*p_IedgesAtBoundary(iebd)-1+p_rtriangulation%NVT) + 1.0_DP

          p_Ddata(2*p_IedgesAtBoundary(iebd)  +p_rtriangulation%NVT) = &
              p_Ddata(2*p_IedgesAtBoundary(iebd)  +p_rtriangulation%NVT) + 1.0_DP
        end if
      end do
    end do

  end subroutine

  ! ---------------------------------------------------------------

  subroutine charFctRealBdComp2d_Ex3x (rboundaryRegion,rvector)
    ! Worker routine for Ex3x~-spaces.
    type(t_boundaryRegion), intent(in), target :: rboundaryRegion
    type(t_vectorScalar), intent(inout), target :: rvector

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: p_IedgesAtBoundary,p_IboundaryCpIdx
    real(dp), dimension(:), pointer :: p_DedgeParameterValue
    integer :: ibct,iebd
    real(DP), dimension(:), pointer :: p_Ddata

    ! Get some data
    p_rtriangulation => rvector%p_rspatialDiscr%p_rtriangulation
    call storage_getbase_int(p_rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    call storage_getbase_int(p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_double(p_rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)
    call lsyssc_getbase_double (rvector,p_Ddata)

    ! Loop through the boundary components.
    do ibct = 1,p_rtriangulation%nbct
      ! Loop through the vertices
      do iebd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-1
        ! Check if the edge is in the region
        if (boundary_isInRegion (rboundaryRegion,ibct,&
            p_DedgeParameterValue(iebd))) then
          ! Add a 1 to the vector
          p_Ddata(p_IedgesAtBoundary(iebd)) = &
              p_Ddata(p_IedgesAtBoundary(iebd)) + 1.0_DP
        end if
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine anprj_analytL2projectionConstr (rvector, rmatrixMass,&
      fcoeff_buildVectorSc_sim, rcollection,&
      rL2ProjectionConfig, rmatrixMassLumped, rvectorTemp1, rvectorTemp2)

!<description>
  ! Converts an analytically given function fcoeff_buildVectorSc_sim
  ! to a finite element vector rvector by using mass matrices.
  ! So the resulting vector is the lumped $L_2$ projection of the
  ! analytically given function. The error induced by mass lumping is
  ! reduced by adding some amount of mass antidiffusion a posteriori.
  ! Flux limiting techniques are employed to determine the amount of
  ! mass antidiffusion that can be applied without violating the
  ! upper and lower bounds set of by the lumped-mass L2-projection.
!</description>

!<input>
  ! The consistent mass matrix of the FE space of rvector.
  type(t_matrixScalar), intent(in) :: rmatrixMass

  ! A callback routine for the function to be discretised. The callback routine
  ! has the same syntax as that for evaluating analytic functions for the
  ! computation of RHS vectors.
  include '../DOFMaintenance/intf_coefficientVectorSc.inc'

  ! OPTIONAL: A pointer to a collection structure. This structure is
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! OPTIONAL: The lumped mass matrix of the FE space of rvector.
  ! If not specified, it is computed on the fly
  type(t_matrixScalar), intent(in), optional, target :: rmatrixMassLumped
!</input>

!<inputoutput>
  ! A scalar vector that receives the $L_2$ projection of the function.
  type(t_vectorScalar), intent(inout) :: rvector

  ! OPTIONAL: A configuration block for the iteration.
  ! If not specified, the standard settings are used.
  type(t_configL2ProjectionByMass), intent(inout), optional :: rL2ProjectionConfig

  ! OPTIONAL: A temporary vector of the same size as rvector.
  type(t_vectorScalar), intent(inout), target, optional :: rvectorTemp1

  ! OPTIONAL: A second temporary vector of the same size as rvector.
  type(t_vectorScalar), intent(inout), target, optional :: rvectorTemp2
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_configL2ProjectionByMass) :: rconfig
    type(t_vectorScalar) :: rvectorAux
    type(t_matrixScalar), pointer :: p_rmatrixMassLumped
    type(t_afcstab) :: rafcstab

    ! Retrieve lumped mass matrix or create it on the fly
    if (present(rmatrixMassLumped)) then
      p_rmatrixMassLumped => rmatrixMassLumped
    else
      allocate(p_rmatrixMassLumped)
      call lsyssc_duplicateMatrix(rmatrixMass, p_rmatrixMassLumped,&
          LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
      call lsyssc_lumpMatrix(p_rmatrixMassLumped, LSYSSC_LUMP_DIAG)
    end if

    ! We need to perform s single defect correction step
    ! preconditioned by the lumped mass matrix
    rconfig%nmaxIterations  = 1
    rconfig%cpreconditioner = 2
    rconfig%domega          = 1.0_DP

    ! Compute the mass-lumped low-order L2-projection
    call anprj_analytL2projectionByMass(rvector, p_rmatrixMassLumped,&
        fcoeff_buildVectorSc_sim, rcollection, rconfig,&
        p_rmatrixMassLumped, rvectorTemp1, rvectorTemp2)

    ! Compute the consistent L2-projection
    call lsyssc_copyVector(rvector, rvectorAux)
    call anprj_analytL2projectionByMass(rvectorAux, rmatrixMass,&
        fcoeff_buildVectorSc_sim, rcollection, rL2ProjectionConfig,&
        p_rmatrixMassLumped, rvectorTemp1, rvectorTemp2)

    ! Initialise the stabilisation structure
    rafcstab%istabilisationSpec= AFCSTAB_NONE
    rafcstab%cprelimitingType = AFCSTAB_PRELIMITING_NONE
    rafcstab%cafcstabType = AFCSTAB_LINFCT_MASS
    call afcsc_initStabilisation(rmatrixMass, rafcstab)
    call afcstab_genEdgeList(rmatrixMass, rafcstab)

    ! Compute the fluxes for the raw mass antidiffusion
    call afcsc_buildFluxFCT(rafcstab, rvectorAux,&
        0.0_DP, 0.0_DP, 1.0_DP, .true., .true.,&
        AFCSTAB_FCTFLUX_EXPLICIT, rmatrix=rmatrixMass)

    ! Apply flux correction to improve the low-order L2-projection
    call afcsc_buildVectorFCT(rafcstab,&
        p_rmatrixMassLumped, rvector, 1._DP, .false.,&
        AFCSTAB_FCTALGO_STANDARD+AFCSTAB_FCTALGO_SCALEBYMASS, rvector)

    ! Release stabilisation structure
    call afcstab_releaseStabilisation(rafcstab)

    ! Release auxiliary vector
    call lsyssc_releaseVector(rvectorAux)

    ! Release temporal lumped mass matrix (if required)
    if (.not.present(rmatrixMassLumped)) then
      call lsyssc_releaseMatrix(p_rmatrixMassLumped)
      deallocate(p_rmatrixMassLumped)
    end if

  end subroutine

end module
