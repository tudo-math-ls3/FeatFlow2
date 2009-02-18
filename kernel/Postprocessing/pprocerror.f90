!#########################################################################
!# ***********************************************************************
!# <name> pprocerror </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains various routines for calculating errors and norms
!# of finite element functions.
!#
!# The following routines can be found in this module:
!#
!# 1.) pperr_scalar
!#     -> Calculate $L_1$-error, $L_2$-error or $H_1$-error to an
!#        analytic reference function or the $L_1$-norm, $L_2$-norm
!#        or $H_1$-norm of a FE function:
!#   $$ int_\Omega u-u_h dx , \qquad int_\Omega \nabla u-\nabla u_h dx $$
!#
!# 2.) pperr_scalarBoundary2d
!#     -> On a 2D boundary segment, calculate $L_1$-error, $L_2$-error 
!#        or $H_1$-error to an analytic reference function or the 
!#        $L_1$-norm, $L_2$-norm or $H_1$-norm of a FE function.
!#   $$ int_\Gamma u-cu_h dx , \qquad int_\Gamma \nabla u-c\nabla u_h dx $$
!#
!# 3.) pperr_scalarErrorEstimate
!#     -> Calculate error to two different scalar vectors of a  FE function:
!#   $$ int_\Omega u_h-u_ref dx $$
!#        where $u_h$ denotes the FE solution vector and $u_ref$ is 
!#        some reference solution vector which is supposed to be a 
!#        better approximation of the true solution.
!#
!# 4.) pperr_blockErrorEstimate
!#     -> Calculate error to two different block vectors of a  FE function:
!#   $$ int_\Omega u_h-u_ref dx $$
!#        where $u_h$ denotes the FE solution vector and $u_ref$ is
!#        some reference solution vector which is supposed to be a 
!#        better approximation of the true solution.
!#
!# 5.) pperr_scalarStandardDeviation
!#     -> Calculate the standard deviation of a scalar vector
!#
!# 6.) pperr_blockStandardDeviation
!#     -> Calculate the standard deviation of a block vector
!#
!# 7.) pperr_scalarTargetFunc
!#     -> Calculate the target functional for a FE function:
!#   $$ int_\Omega w(x)u(x) - w(x)u_h(x) dx $$
!#        where w(x) is a weighting function
!#
!# </purpose>
!#########################################################################

module pprocerror

  use fsystem
  use storage
  use boundary
  use cubature
  use triangulation
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use scalarpde
  use spatialdiscretisation
  use domainintegration
  use elementpreprocessing
  use feevaluation
  use collection

  implicit none

!<constants>

!<constantblock description = "Identifiers for the type of error to be computed.">

  ! $L_2$-error/norm
  integer, parameter :: PPERR_L2ERROR = 1
  
  ! $H_1$-error/norm
  integer, parameter :: PPERR_H1ERROR = 2

  ! $L_1$-error/norm
  integer, parameter :: PPERR_L1ERROR = 3
  
  
!</constantblock>

!<constantblock description="Constants defining the blocking of the error calculation.">

  ! Number of elements to handle simultaneously when building vectors
  integer :: PPERR_NELEMSIM   = 1000
  
!</constantblock>

!</constants>

!<types>

!<typeblock>

  type t_errorScVec
  
    ! IN: An array of scalar coefficient vectors which represents a
    ! FE function or vector field. If p_rdiscr is null, then all vectors in
    ! the array must have the same spatial discretisation!
    type(t_vectorScalar), dimension(:), pointer :: p_RvecCoeff => null()

    ! IN, OPTIONAL: A spatial discretisation structure that is to be used
    ! for the scalar coefficient vectors. If given, the total number of DOFs
    ! of the spatial discretisation must be equal to the number of entries of
    ! each scalar vector given in p_RvecCoeff. If not given, the spatial
    ! discretisation of p_RvecCoeff is used.
    type(t_spatialDiscretisation), pointer :: p_rdiscr => null()
    
    ! OUT, OPTIONAL: If given, recieves the calculated L2-errors for each
    ! component of p_RvecCoeff.
    real(DP), dimension(:), pointer :: p_DerrorL2 => null()
    
    ! OUT, OPTIONAL: If given, recieves the calculated H1-errors for each
    ! component of p_RvecCoeff. If given, but the spatial discretisation does
    ! not offer first derivatives, then all components of p_DerrorH1 are set
    ! to SYS_INFINITY to indicate that the H1-error cannot be calculated.
    real(DP), dimension(:), pointer :: p_DerrorH1 => null()
    
    ! OUT, OPTIONAL: If given, recieves the calculated L1-errors for each
    ! component of p_RvecCoeff.
    real(DP), dimension(:), pointer :: p_DerrorL1 => null()
    
    ! OUT, OPTIONAL: An array of scalar vectors which recieve the L2-errors
    ! of each component per element.
    type(t_vectorScalar), dimension(:), pointer :: p_RvecErrorL2 => null()
    
    ! OUT, OPTIONAL: An array of scalar vectors which recieve the H1-errors
    ! of each component per element.
    type(t_vectorScalar), dimension(:), pointer :: p_RvecErrorH1 => null()
    
    ! OUT, OPTIONAL: An array of scalar vectors which recieve the L1-errors
    ! of each component per element.
    type(t_vectorScalar), dimension(:), pointer :: p_RvecErrorL1 => null()
  
  end type

!</typeblock>

!</types>

contains

  !****************************************************************************

!<subroutine>
  
  subroutine pperr_scalarVec(rerror, frefFunction, rcollection)

!<description>
  ! This routine calculates the errors of a set of given FE functions given
  ! analytical callback function frefFunction.
  ! In contrast to pperr_scalar, this routine is able to compute multiple
  ! errors in multiple components of a vector field at once.
!</description>

!<input>
  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_refFunctionScVec.inc'
  optional :: frefFunction
  
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information. 
  type(t_collection), intent(INOUT), target, optional :: rcollection
!</input>

!<inputoutput>
  ! A structure which defines what errors are to be calculated.
  type(t_errorScVec), intent(INOUT), target :: rerror
!</inputoutput>

!</subroutine>

  ! A pointer to the discretisation that is to be used
  type(t_spatialDiscretisation), pointer :: p_rdiscr
  
  ! A pointer to an element distribution
  type(t_elementDistribution), pointer :: p_relemDist
  
  ! An array holding the element list
  integer, dimension(:), pointer :: p_IelemList, p_IcurElemList
  
  ! A pointer to the triangulation
  type(t_triangulation), pointer :: p_rtria
  
  ! Which errors do we have to calculate?
  logical :: bcalcL2, bcalcH1, bcalcL1
  
  ! The total number of components
  integer :: ncomp,icomp
  
  ! Indices of first and last derivative
  integer :: ifirstDer, ilastDer
  
  ! Arrays concerning the cubature formula
  real(DP), dimension(:,:), allocatable :: DcubPts
  real(DP), dimension(:), allocatable :: Domega
  
  ! An array needed for the DOF-mapping
  integer, dimension(:,:), target, allocatable :: Idofs

  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  type(t_domainIntSubset) :: rintSubset
  type(t_evalElementSet) :: reval
  
  ! Two arrays for the function values and derivatives
  real(DP), dimension(:,:), allocatable :: DvalFunc
  real(DP), dimension(:,:,:), allocatable :: DvalDer
  
  ! Two arrays for the element evaluation
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: Dbas
  
  ! A pointer to the data array of the currently active coefficient vector
  real(DP), dimension(:), pointer :: p_Dcoeff
  
  ! Pointers to the arrays that recieve th element-wise errors
  real(DP), dimension(:), pointer :: p_DerrL2, p_DerrH1, p_DerrL1
  
  ! Some other local variables
  integer :: i,j,k,NEL,ieldist,ndofs,ncubp,ctrafo,ccubature,iel,ider
  integer :: NELtodo,NELdone,NELbatchSize
  integer(I32) :: cevalTag, celement
  real(DP) :: derrL2, derrH1, derrL1, dom, daux, daux2
  real(DP), dimension(:,:), pointer :: p_Ddetj

    ! Make sure we have the coefficient vectors
    if(.not. associated(rerror%p_RvecCoeff)) then
      call output_line('Coefficient vectors missing!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
      call sys_halt()
    end if
    
    ! Get the number of components
    ncomp = ubound(rerror%p_RvecCoeff,1)
    
    ! Do we have a separate discretisation?
    if(associated(rerror%p_rdiscr)) then
      ! Yes, so check whether the discretisation is compatible to the
      ! coefficient vectors.
      p_rdiscr => rerror%p_rdiscr
      k = dof_igetNDofGlob(p_rdiscr)
      do i = 1, ncomp
        if(k .ne. rerror%p_RvecCoeff(i)%NEQ) then
          call output_line('Discretisation and coefficient vectors incompatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
          call sys_halt()
        end if
      end do
      
    else
      ! No, so grab the discretisation of the first coefficient vector.
      p_rdiscr => rerror%p_RvecCoeff(1)%p_rspatialDiscr
      if(.not. associated(p_rdiscr)) then
        call output_line('No discretisation assigned!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
        call sys_halt()
      end if
    end if
    
    ! Get the triangulation and the number of elements
    p_rtria => p_rdiscr%p_rtriangulation
    NEL = p_rtria%NEL
    
    ! Okay, now that we have the discretisation, determine the dimension
    ! to figure out which is the first and the last derivative we need to
    ! evaluate in the case that we want to compute H1-errors.
    select case(p_rdiscr%ndimension)
    case (NDIM1D)
      ifirstDer = DER_DERIV1D_X
      ilastDer  = DER_DERIV1D_X
    case (NDIM2D)
      ifirstDer = DER_DERIV2D_X
      ilastDer  = DER_DERIV2D_Y
    case (NDIM3D)
      ifirstDer = DER_DERIV3D_X
      ilastDer  = DER_DERIV3D_Z
    case default
      call output_line('Invalid discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
      call sys_halt()
    end select
    
    ! Now determine which errors we are going to calculate
    bcalcL2 = .false.
    bcalcH1 = .false.
    bcalcL1 = .false.
    if(associated(rerror%p_DerrorL2)) then
      bcalcL2 = .true.
      if(ubound(rerror%p_DerrorL2,1) .lt. ncomp) then
        call output_line('Dimension of p_DerrorL2 array is too small!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
        call sys_halt()
      end if
      rerror%p_DerrorL2 = 0.0_DP
    end if
    if(associated(rerror%p_DerrorH1)) then
      bcalcH1 = .true.
      if(ubound(rerror%p_DerrorH1,1) .lt. ncomp) then
        call output_line('Dimension of p_DerrorH1 array is too small!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
        call sys_halt()
      end if
      rerror%p_DerrorH1 = 0.0_DP
    end if
    if(associated(rerror%p_DerrorL1)) then
      bcalcL1 = .true.
      if(ubound(rerror%p_DerrorL1,1) .lt. ncomp) then
        call output_line('Dimension of p_DerrorL1 array is too small!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
        call sys_halt()
      end if
      rerror%p_DerrorL1 = 0.0_DP
    end if
    if(associated(rerror%p_RvecErrorL2)) then
      bcalcL2 = .true.
      if(ubound(rerror%p_RvecErrorL2,1) .lt. ncomp) then
        call output_line('Dimension of p_RvecErrorL2 array is too small!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
        call sys_halt()
      end if
      do i = 1, ncomp
        if(rerror%p_RvecErrorL2(i)%NEQ .lt. NEL) then
          call output_line('Length of p_RvecErrorL2('//trim(sys_siL(i,4))//&
              ') is too small!',OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
          call sys_halt()
        end if
        call lsyssc_clearVector(rerror%p_RvecErrorL2(i))
      end do
    end if
    if(associated(rerror%p_RvecErrorH1)) then
      bcalcH1 = .true.
      if(ubound(rerror%p_RvecErrorH1,1) .lt. ncomp) then
        call output_line('Dimension of p_RvecErrorH1 array is too small!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
        call sys_halt()
      end if
      do i = 1, ncomp
        if(rerror%p_RvecErrorH1(i)%NEQ .lt. NEL) then
          call output_line('Length of p_RvecErrorH1('//trim(sys_siL(i,4))//&
              ') is too small!',OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
          call sys_halt()
        end if
        call lsyssc_clearVector(rerror%p_RvecErrorH1(i))
      end do
    end if
    if(associated(rerror%p_RvecErrorL1)) then
      bcalcL1 = .true.
      if(ubound(rerror%p_RvecErrorL1,1) .lt. ncomp) then
        call output_line('Dimension of p_RvecErrorL1 array is too small!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
        call sys_halt()
      end if
      do i = 1, ncomp
        if(rerror%p_RvecErrorL1(i)%NEQ .lt. NEL) then
          call output_line('Length of p_RvecErrorL1('//trim(sys_siL(i,4))//&
              ') is too small!',OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarVec')
          call sys_halt()
        end if
        call lsyssc_clearVector(rerror%p_RvecErrorL1(i))
      end do
    end if
    
    ! Don't we have anything to do?
    if(.not. (bcalcL2 .or. bcalcH1 .or. bcalcL1)) return
    
    ! Do we have to calculate H1 errors?
    if(bcalcH1) then
    
      ! Okay, in this case all element distributions of the spatial
      ! discretisation must support first derivatives.
      do i = 1, p_rdiscr%inumFEspaces
        
        ! Get a pointer to the element distribution
        p_relemDist => p_rdiscr%RelementDistr(i)
        
        ! If the maximum supported derivative is 1, then the element does not
        ! support first derivatives!
        if(elem_getMaxDerivative(p_relemDist%celement) .le. 1) then
          bcalcH1 = .false.
          exit
        end if
      end do
      
      ! Now if bcalcH1 is .false. now, then at least one element distribution
      ! does not support first derivatives...
      if(.not. bcalcH1) then
        ! If p_DerrorH1 is given, set its entries to SYS_INFTY to indicate that
        ! the H1-errors are not available
        if(associated(rerror%p_DerrorH1)) &
          rerror%p_DerrorH1 = SYS_INFINITY
        
      end if
      
    end if
    
    ! Set up the Bder array
    Bder = .false.
    Bder(DER_FUNC) = bcalcL2 .or. bcalcL1
    if(bcalcH1) Bder(ifirstDer:ilastDer) = .true.
    
    ! Okay, let's loop over all element distributions
    do ieldist = 1, p_rdiscr%inumFEspaces
    
      ! Get a pointer to the element distribution
      p_relemDist => p_rdiscr%RelementDistr(ieldist)
      if(p_relemDist%NEL .le. 0) cycle
      
      ! Get a pointer to the element list
      call storage_getbase_int(p_relemDist%h_IelementList, p_IelemList)
      
      ! Calculate element batch size
      NEL = p_relemDist%NEL
      NELbatchSize = min(NEL, PPERR_NELEMSIM)
      
      ! Get the element
      celement = p_relemDist%celement
      
      ! Get the number of dofs per element
      ndofs = elem_igetNDofLoc(celement)
      
      ! Get the trafo
      ctrafo = elem_igetTrafoType(celement)
      
      ! Get the evaluation tag
      cevalTag = elem_getEvaluationTag(celement)
      
      ! If a reference function is given, it will surely need real points.
      if(present(frefFunction)) &
        cevalTag = ior(cevalTag, EL_EVLTAG_REALPOINTS)
      
      ! And we definately need jacobian determinants for integration.
      cevalTag = ior(cevalTag, EL_EVLTAG_DETJ)
      
      ! Get the cubature rule
      ccubature = p_relemDist%ccubTypeEval

      ! Allocate the arrays for the cubature formula
      ncubp = cub_igetNumPts(ccubature)
      allocate(Domega(ncubp))
      allocate(DcubPts(trafo_igetReferenceDimension(ctrafo),ncubp))
      
      ! Get the cubature formula
      call cub_getCubature(ccubature,DcubPts, Domega)
      
      ! Allocate an array for the DOF-mapping
      allocate(Idofs(ndofs,NELbatchSize))
      
      ! Allocate an array that recieves the evaluated basis functions
      allocate(Dbas(ndofs,elem_getMaxDerivative(celement),ncubp,NELbatchsize))
      
      ! Allocate two arrays for the evaluation
      if(bcalcL2 .or. bcalcL1) &
        allocate(DvalFunc(ncubp,NELbatchSize))
      if(bcalcH1) &
        allocate(DvalDer(ncubp,NELbatchSize,ifirstDer:ilastDer))
      
      ! Initialise the element evaluation set
      call elprep_init(reval)
      
      ! Okay, loop over all elements
      NELdone = 0
      do while(NELdone .lt. NEL)
      
        ! How many elements do we process this time?
        NELtodo = min(NELbatchSize, NEL-NELdone)
        
        ! Get a pointer to the current element list
        p_IcurElemList => p_IelemList(NELdone+1:NELdone+NELtodo)
        
        ! First, let's perform the DOF-mapping
        call dof_locGlobMapping_mult(p_rdiscr, p_IcurElemList, Idofs)
        
        ! Prepare the element for evaluation
        call elprep_prepareSetForEvaluation (reval, cevalTag, p_rtria, &
            p_IcurElemList, ctrafo, DcubPts)
        p_Ddetj => reval%p_Ddetj
        
        ! Remove the ref-points eval tag for the next loop iteration
        cevalTag = iand(cevalTag,not(EL_EVLTAG_REFPOINTS))

        ! Prepare the domain integration structure
        call domint_initIntegrationByEvalSet(reval, rintSubset)
        rintSubset%ielementDistribution = ielDist
        rintSubset%ielementStartIdx = NELdone+1
        rintSubset%p_Ielements => p_IcurElemList
        rintSubset%p_IdofsTrial => Idofs
        rintSubset%celement = celement
        
        ! Evaluate the element
        call elem_generic_sim2(celement, reval, Bder, Dbas)
        
        ! Now loop over all vector components
        do icomp = 1, ncomp
        
          ! Get the coefficient vector's data array
          call lsyssc_getbase_double(rerror%p_RvecCoeff(icomp), p_Dcoeff)
          
          ! Get the element-wise error arrays, if given
          if(associated(rerror%p_RvecErrorL2)) then
            call lsyssc_getbase_double(rerror%p_RvecErrorL2(icomp), p_DerrL2)
          else
            nullify(p_DerrL2)
          end if
          if(associated(rerror%p_RvecErrorH1)) then
            call lsyssc_getbase_double(rerror%p_RvecErrorH1(icomp), p_DerrH1)
          else
            nullify(p_DerrH1)
          end if
          if(associated(rerror%p_RvecErrorL1)) then
            call lsyssc_getbase_double(rerror%p_RvecErrorL2(icomp), p_DerrL1)
          else
            nullify(p_DerrL1)
          end if

          ! Reset errors for this component
          derrL2 = 0.0_DP
          derrH1 = 0.0_DP
          derrL1 = 0.0_DP
          
          ! Evaluate function values?
          if(allocated(DvalFunc)) then
          
            ! Do we have a reference function? If yes, then evaluate it,
            ! otherwise simply format DvalFunc to zero.
            if(present(frefFunction)) then
              call frefFunction(icomp, DER_FUNC, p_rdiscr, NELtodo, ncubp, &
                  reval%p_DpointsReal, rintSubset, DvalFunc, rcollection)
            else
              DvalFunc = 0.0_DP
            end if
            
            ! Now subtract the function values of the FE function.
            !$omp parallel do private(i,k,daux)
            do j = 1, NELtodo
              do i = 1, ncubp
                daux = 0.0_DP
                do k = 1, ndofs
                  daux = daux + Dbas(k,DER_FUNC,i,j)*p_Dcoeff(Idofs(k,j))
                end do ! k
                DvalFunc(i,j) = DvalFunc(i,j) - daux
              end do ! i
            end do ! j
            !$omp end parallel do
          
          end if ! function values evaluation
          
          ! Evaluate derivatives?
          if(allocated(DvalDer)) then
          
            ! Do we have a reference function? If yes, then evaluate its
            ! derivatives.
            if(present(frefFunction)) then
              do ider = ifirstDer, ilastDer
                call frefFunction(icomp, ider, p_rdiscr, NELtodo, ncubp, &
                    reval%p_DpointsReal, rintSubset, DvalDer(:,:,ider), rcollection)
              end do
            else
              DvalDer = 0.0_DP
            end if
            
            ! Now subtract the derivatives of the FE function.
            !$omp parallel do private(i,ider,k,daux)
            do j = 1, NELtodo
              do i = 1, ncubp
                do ider = ifirstDer, ilastDer
                  daux = 0.0_DP
                  do k = 1, ndofs
                    daux = daux + Dbas(k,ider,i,j)*p_Dcoeff(Idofs(k,j))
                  end do ! k
                  DvalDer(i,j,ider) = DvalDer(i,j,ider) - daux
                end do ! ider
              end do ! i
            end do ! j
            !$omp end parallel do
          
          end if ! derivatives evaluation
            
          ! Do we calculate L2-errors?
          if(bcalcL2 .and. associated(p_DerrL2)) then
            !$omp parallel do private(i,iel,daux,dom) reduction(+:derrL2)
            do j = 1, NELtodo
              iel = p_IcurElemList(j)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*DvalFunc(i,j)**2
              end do ! i
              p_DerrL2(iel) = sqrt(daux)
              derrL2 = derrL2 + daux
            end do ! j
            !$omp end parallel do
          else if(bcalcL2) then
            !$omp parallel do private(i,daux,dom) reduction(+:derrL2)
            do j = 1, NELtodo
              daux = 0.0_DP
              do i = 1, ncubp
                dom = Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*DvalFunc(i,j)**2
              end do ! i
              derrL2 = derrL2 + daux
            end do ! j
            !$omp end parallel do
          end if
          
          ! Do we calculate H1-errors?
          if(bcalcH1 .and. associated(p_DerrH1)) then
            !$omp parallel do private(i,ider,iel,daux,daux2,dom) reduction(+:derrH1)
            do j = 1, NELtodo
              iel = p_IcurElemList(j)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = Domega(i) * abs(p_Ddetj(i,j))
                daux2 = 0.0_DP
                do ider = ifirstDer, ilastDer
                  daux2 = daux2 + DvalDer(i,j,ider)**2
                end do ! ider
                daux = daux + dom*daux2
              end do ! i
              p_DerrH1(iel) = sqrt(daux)
              derrH1 = derrH1 + daux
            end do ! j
            !$omp end parallel do
          else if(bcalcH1) then
            !$omp parallel do private(i,ider,daux,daux2,dom) reduction(+:derrH1)
            do j = 1, NELtodo
              daux = 0.0_DP
              do i = 1, ncubp
                dom = Domega(i) * abs(p_Ddetj(i,j))
                daux2 = 0.0_DP
                do ider = ifirstDer, ilastDer
                  daux2 = daux2 + DvalDer(i,j,ider)**2
                end do ! ider
                daux = daux + dom*daux2
              end do ! i
              derrH1 = derrH1 + daux
            end do ! j
            !$omp end parallel do
          end if

          ! Do we calculate L1-errors?
          if(bcalcL1 .and. associated(p_DerrL1)) then
            !$omp parallel do private(i,iel,daux,dom) reduction(+:derrL1)
            do j = 1, NELtodo
              iel = p_IcurElemList(j)
              daux = 0.0_DP
              do i = 1, ncubp
                dom = Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*abs(DvalFunc(i,j))
              end do ! i
              p_DerrL1(iel) = daux
              derrL2 = derrL2 + daux
            end do ! j
            !$omp end parallel do
          else if(bcalcL1) then
            !$omp parallel do private(i,daux,dom) reduction(+:derrL1)
            do j = 1, NELtodo
              daux = 0.0_DP
              do i = 1, ncubp
                dom = Domega(i) * abs(p_Ddetj(i,j))
                daux = daux + dom*abs(DvalFunc(i,j))
              end do ! i
              derrL1 = derrL1 + daux
            end do ! j
            !$omp end parallel do
          end if
          
          ! Incorporate errors
          if(bcalcL2 .and. associated(rerror%p_DerrorL2)) &
            rerror%p_DerrorL2(icomp) = rerror%p_DerrorL2(icomp) + derrL2
          if(bcalcH1 .and. associated(rerror%p_DerrorH1)) &
            rerror%p_DerrorH1(icomp) = rerror%p_DerrorH1(icomp) + derrH1
          if(bcalcL1 .and. associated(rerror%p_DerrorL1)) &
            rerror%p_DerrorL1(icomp) = rerror%p_DerrorL1(icomp) + derrL1
        
        end do ! icomp
      
        ! Release the domain integration structure
        call domint_doneIntegration (rintSubset)
        
        NELdone = NELdone + NELtodo

      end do ! while(NELdone .lt. NEL)

      ! Release the element evaluation set
      call elprep_releaseElementSet(reval)
      
      ! Deallocate all arrays
      if(allocated(DvalDer)) deallocate(DvalDer)
      if(allocated(DvalFunc)) deallocate(DvalFunc)
      deallocate(Dbas)
      deallocate(Idofs)
      deallocate(DcubPts)
      deallocate(Domega)
    
    end do ! ieldist
    
    ! Don't forget to take the square roots of the L2- and H1-errors
    if(associated(rerror%p_DerrorL2) .and. bcalcL2) then
      do icomp = 1, ncomp
        rerror%p_DerrorL2(icomp) = sqrt(rerror%p_DerrorL2(icomp))
      end do
    end if
    if(associated(rerror%p_DerrorH1) .and. bcalcH1) then
      do icomp = 1, ncomp
        rerror%p_DerrorH1(icomp) = sqrt(rerror%p_DerrorH1(icomp))
      end do
    end if
    
    ! That's it

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalar (rvectorScalar, cerrortype, derror,&
                           ffunctionReference, rcollection,&
                           rdiscretisation, relementError)

!<description>
  ! This routine calculates the error or the norm, respectively, of a given 
  ! finite element function in rvector to a given analytical 
  ! callback function ffunctionReference.
  !
  ! If ffunctionReference is specified, the routine calculates
  !   $$ ||y-z||_{L_1}, ||y-z||_{L_2}  \textrm{ or }  ||y-z||_{H_1}$$
  ! with $y$=rvectorScalar and $z$=ffunctionReference.
  !
  ! If ffunctionReference is not specified, the routine calculates
  !   $$ ||y||_{L_1}, ||y||_{L_2}  \textrm{ or }  ||y||_{H_1}.$$
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the
  ! cubature formula to use for each element distribution.
  !
  ! If the H1-error is desired and the element does not provide first
  ! derivatives, then this routine sets derror to -1 to indicate that the
  ! calculation of the H1-error is not available.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(IN), target :: rvectorScalar
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(IN) :: cerrortype
  
  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_refFunctionSc.inc'
  optional :: ffunctionReference
  
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information. 
  type(t_collection), intent(INOUT), target, optional :: rcollection
  
  ! OPTIONAL: A discretisation structure specifying how to compute the error.
  ! If not specified, the discretisation structure in the vector is used.
  ! If specified, the discretisation structure must be 'compatible' to the
  ! vector (concerning NEQ,...). pperr_scalar uses the cubature formula
  ! specifier of the linear form in rdiscretisation to compute the integrals
  ! for the error.
  type(t_spatialDiscretisation), intent(IN), target, optional :: rdiscretisation
!</input>

!<inputoutput>
  ! OPTIONAL: A scalar vector which holds the calculated error per element
  type(t_vectorScalar), intent(INOUT), optional :: relementError
!</inputoutput>

!<output>
  ! The calculated error.
  real(DP), intent(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    integer :: i, imaxder
    integer(I32) :: celement
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    
    ! Get the correct discretisation structure and check if we can use it.
    if (present(rdiscretisation)) then
      p_rdiscretisation => rdiscretisation
      call lsyssc_checkDiscretisation (rvectorScalar,p_rdiscretisation)
    else
      p_rdiscretisation => rvectorScalar%p_rspatialdiscr
    end if
    
    if (.not. associated(p_rdiscretisation)) then
      call output_line('No discretisation structure!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar')
      call sys_halt()
    end if
    
    ! The vector must be unsorted, otherwise we can't set up the vector.
    if (rvectorScalar%isortStrategy .gt. 0) then
      call output_line('Vector must be unsorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar')
      call sys_halt()
    end if
    
    ! There is one thing we need to assure: If the user wants to calculate
    ! an H1-error, the corresponding element must have first derivatives!
    if(cerrortype .eq. PPERR_H1ERROR) then
      
      ! So, let's loop through all element distributions of the spatial
      ! discretisation structure.
      do i = 1, p_rdiscretisation%inumFESpaces
        
        ! Get the element of the FE space
        celement = p_rdiscretisation%RelementDistr(i)%celement
        
        ! Get the maximum derivative of the element
        imaxder = elem_getMaxDerivative(celement)
        
        ! If the maximum derivate is 1 then the element's derivatives can
        ! not be evaluated (since they are zero). We will simply return
        ! an error of -1 to indicate that the H1-error cannot be computed.
        if(imaxder .le. 1) then
          derror = -1.0_DP
          return
        end if
      end do
      
      ! If we come out here, then the element provides first derivatives, and
      ! we are ready to calculate the H1-error.
    
    end if
  
    ! Do we have a uniform triangulation? Would simplify a lot...
    if ((p_rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) .or.&
        (p_rdiscretisation%ccomplexity .eq. SPDISC_CONFORMAL)) then 
    
      select case(rvectorScalar%cdataType)

      case (ST_DOUBLE)
        ! Do we have a 1D, 2D or 3D discretisation here?
        select case(p_rdiscretisation%ndimension)
        case (NDIM1D)
          call pperr_scalar1d_conf (rvectorScalar,cerrortype,derror,&
                                    p_rdiscretisation,ffunctionReference,&
                                    rcollection, relementError)
        case (NDIM2D)
          call pperr_scalar2d_conf (rvectorScalar,cerrortype,derror,&
                                    p_rdiscretisation,ffunctionReference,&
                                    rcollection, relementError)
        case (NDIM3D)
          call pperr_scalar3d_conf (rvectorScalar,cerrortype,derror,&
                                    p_rdiscretisation,ffunctionReference,&
                                    rcollection, relementError)
        end select

      case DEFAULT
        call output_line('Single precision vectors currently not supported!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar')
        call sys_halt()
      end select
    
    else
      call output_line('General discretisation not implemented!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar')
      call sys_halt()
    end if

  end subroutine pperr_scalar

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalar1d_conf (rvectorScalar, cerrortype, derror,&
                                  rdiscretisation, ffunctionReference,&
                                  rcollection, relementError)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 1D version for double-precision vectors.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(IN), target :: rvectorScalar
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(IN) :: cerrortype
  
  ! A discretisation structure specifying how to compute the error.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisation
  
  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_refFunctionSc.inc'
  optional :: ffunctionReference
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional :: rcollection

  ! OPTIONAL: A scalar vector which holds the calculated error per element
  type(t_vectorScalar), intent(INOUT), optional :: relementError
!</inputoutput>

!<output>
  ! Array receiving the calculated error.
  real(DP), intent(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    integer :: icurrentElementDistr, ICUBP, NVE
    integer :: IEL, IELmax, IELset, IELGlobal
    real(DP) :: OM
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(:), allocatable :: Domega
    
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

    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet
    
    ! Type of transformation from the reference to the real element 
    integer :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial
  
    ! Pointer to the element-wise error
    real(DP), dimension(:), pointer :: p_Derror


    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.

    Bder = .false.
    select case (cerrortype)
    case (PPERR_L1ERROR, PPERR_L2ERROR) 
      Bder(DER_FUNC1D) = .true.
    case (PPERR_H1ERROR) 
      Bder(DER_DERIV1D_X) = .true.
    case default
      call output_line('Unknown error type identifier!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar1d_conf')
      call sys_halt()
    end select
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPERR_NELEMSIM,p_rtriangulation%NEL)
                               
    ! Set the current error to 0 and add the error contributions of each element
    ! to that.
    derror = 0.0_DP

    ! Set pointer to element-wise error
    if (present(relementError)) then
      call lsyssc_getbase_double(relementError, p_Derror)
    end if

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1, rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Get the number of cubature points for the cubature formula
      ncubp = cub_igetNumPts(p_relementDistribution%ccubTypeEval)
      
      ! Allocate two arrays for the points and the weights
      allocate(Domega(ncubp))
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType), ncubp))
      
      ! Get the cubature formula
      call cub_getCubature(p_relementDistribution%ccubTypeEval, p_DcubPtsRef, Domega)

      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial, nelementsPerBlock))

      ! Allocate memory for the coefficients
      allocate(Dcoefficients(ncubp, nelementsPerBlock, 2))
    
      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      if (present(ffunctionReference)) then
        ! Evaluate real coordinates if necessary.
        cevaluationTag = ior(cevaluationTag, EL_EVLTAG_REALPOINTS)
      end if

      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL
    
      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, PPERR_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
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
        p_Ddetj => revalElementSet%p_Ddetj

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        select case (cerrortype)
        
        case (PPERR_L2ERROR)
          
          ! L2-error uses only the values of the function.
          
          if (present(ffunctionReference)) then
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x)
            ! The result is saved in Dcoefficients(:,:,1)
            call ffunctionReference (DER_FUNC1D, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp, &
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset, &
                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
          else
            Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP
          end if

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x).
          ! Save the result to Dcoefficients(:,:,2)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC1D,&
                  Dcoefficients(:,1:IELmax-IELset+1,2))
                  
          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (u-u_h,u-u_h) dx
          
          if (present(relementError)) then

            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
                
                IELGlobal = p_IelementList(IELset+IEL-1)
                
                p_Derror(IELGlobal) = OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2

                derror = derror + p_Derror(IELGlobal)
                
              end do ! ICUBP 
              
            end do ! IEL

          else

            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
                
                derror = derror + &
                         OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2
                
              end do ! ICUBP 
              
            end do ! IEL

          end if

        case (PPERR_L1ERROR)
          
          ! L1-error uses only the values of the function.
          
          if (present(ffunctionReference)) then
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x)
            ! The result is saved in Dcoefficients(:,:,1)
            call ffunctionReference (DER_FUNC1D, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
          else
            Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP
          end if

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x).
          ! Save the result to Dcoefficients(:,:,2)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC1D,&
                  Dcoefficients(:,1:IELmax-IELset+1,2))

          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega abs(u-u_h) dx
          
          if (present(relementError)) then

            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! L1-error is:   int_... abs(u-u_h) dx
                
                IELGlobal = p_IelementList(IELset+IEL-1)

                p_Derror(IELGlobal) = OM * abs(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))

                derror = derror + p_Derror(IELGlobal)
                
              end do ! ICUBP 
              
            end do ! IEL

          else

            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! L1-error is:   int_... abs(u-u_h) dx
                
                derror = derror + &
                         OM * abs(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))
                
              end do ! ICUBP 
              
            end do ! IEL
            
          end if

        case (PPERR_H1ERROR)

          ! H1-error uses only 1st derivative of the function.

          if (present(ffunctionReference)) then          
            ! It's time to call our coefficient function to calculate the
            ! X-derivative values in the cubature points:  u(x,y)
            ! The result is saved in Dcoefficients(:,:,1)
            call ffunctionReference (DER_DERIV1D_X, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
                        
          else
            Dcoefficients(:,1:IELmax-IELset+1,1:2) = 0.0_DP
          end if
          
          ! Calculate the X/Y-derivative of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,3)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV1D_X,&
                  Dcoefficients(:,1:IELmax-IELset+1,2))

          ! Subtraction of Dcoefficients(:,:,1..2) from Dcoefficients(:,:,3..4) gives
          ! the error "grad(u-u_h)(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
          
          if (present(relementError)) then

            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
                
                IELGlobal = p_IelementList(IELset+IEL-1)

                p_Derror(IELGlobal) = OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2

                derror = derror + p_Derror(IELGlobal)

              end do ! ICUBP 
              
            end do ! IEL

          else
            
            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
                
                derror = derror + OM * &
                         (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2
                
              end do ! ICUBP 
              
            end do ! IEL

          end if
        
        case default
          call output_line('Unknown error type identifier!',&
              OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar1d_conf')
          call sys_halt()
        end select
        
        ! Release the temporary domain integration structure again
        call domint_doneIntegration (rintSubset)
    
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)
      deallocate(Domega)

    end do ! icurrentElementDistr

    ! derror is ||error||^2, so take the square root at last.
    if ((cerrortype .eq. PPERR_L2ERROR) .or.&
        (cerrortype .eq. PPERR_H1ERROR)) then
      derror = sqrt(derror)
      if (present(relementError)) then
        do IEL = 1, size(p_Derror,1)
          p_Derror(IEL) = sqrt(p_Derror(IEL))
        end do
      end if
    end if

  end subroutine pperr_scalar1d_conf

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalar2d_conf (rvectorScalar, cerrortype, derror,&
                                  rdiscretisation, ffunctionReference,&
                                  rcollection, relementError)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(IN), target :: rvectorScalar
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(IN) :: cerrortype
  
  ! A discretisation structure specifying how to compute the error.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisation

  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_refFunctionSc.inc'
  optional :: ffunctionReference
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional 
  ! information for callback routines.
  type(t_collection), intent(INOUT), optional :: rcollection

  ! OPTIONAL: A scalar vector which holds the calculated error per element
  type(t_vectorScalar), intent(INOUT), optional :: relementError
!</inputoutput>

!<output>
  ! Array receiving the calculated error.
  real(DP), intent(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    integer :: icurrentElementDistr, ICUBP, NVE
    integer :: IEL, IELmax, IELset, IELGlobal
    real(DP) :: OM
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(:), allocatable :: Domega
    
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

    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial
  
    ! Type of transformation from the reference to the real element 
    integer :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! Pointer to the element-wise error
    real(DP), dimension(:), pointer :: p_Derror


    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.

    Bder = .false.
    select case (cerrortype)
    case (PPERR_L1ERROR, PPERR_L2ERROR) 
      Bder(DER_FUNC) = .true.
    case (PPERR_H1ERROR) 
      Bder(DER_DERIV_X) = .true.
      Bder(DER_DERIV_Y) = .true.
    case default
      call output_line('Unknown error type identifier!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar2d_conf')
      call sys_halt()
    end select
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPERR_NELEMSIM, p_rtriangulation%NEL)
    
    ! Set the current error to 0 and add the error contributions of each element
    ! to that.
    derror = 0.0_DP

    ! Set pointer to element-wise error
    if (present(relementError)) then
      call lsyssc_getbase_double(relementError, p_Derror)
    end if

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1, rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Get the number of cubature points for the cubature formula
      ncubp = cub_igetNumPts(p_relementDistribution%ccubTypeEval)
      
      ! Allocate two arrays for the points and the weights
      allocate(Domega(ncubp))
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType), ncubp))
      
      ! Get the cubature formula
      call cub_getCubature(p_relementDistribution%ccubTypeEval, p_DcubPtsRef, Domega)
      
      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial, nelementsPerBlock))

      ! Allocate memory for the coefficients
      allocate(Dcoefficients(ncubp, nelementsPerBlock, 4))
    
      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      if (present(ffunctionReference)) then
        ! Evaluate real coordinates if not necessary.
        cevaluationTag = ior(cevaluationTag, EL_EVLTAG_REALPOINTS)
      end if
                      
      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag, EL_EVLTAG_DETJ)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL
    
      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, PPERR_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
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
        p_Ddetj => revalElementSet%p_Ddetj

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        select case (cerrortype)
        
        case (PPERR_L2ERROR)
          
          ! L2-error uses only the values of the function.
          
          if (present(ffunctionReference)) then
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x,y)
            ! The result is saved in Dcoefficients(:,:,1)
            call ffunctionReference (DER_FUNC, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
          else
            Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP
          end if

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,2)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1,2))

          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (u-u_h,u-u_h) dx
          
          if (present(relementError)) then

            do IEL=1,IELmax-IELset+1
          
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! L2-error is:   int_... (u-u_h)*(u-u_h) dx

                IELGlobal = p_IelementList(IELset+IEL-1)
                
                p_Derror(IELGlobal) = OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2

                derror = derror + p_Derror(IELGlobal)
                
              end do ! ICUBP 
              
            end do ! IEL

          else

            do IEL=1,IELmax-IELset+1
          
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
                
                derror = derror + &
                         OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2
                
              end do ! ICUBP 
              
            end do ! IEL

          end if

        case (PPERR_L1ERROR)
          
          ! L1-error uses only the values of the function.
          
          if (present(ffunctionReference)) then
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x,y)
            ! The result is saved in Dcoefficients(:,:,1)
            call ffunctionReference (DER_FUNC, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
          else
            Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP
          end if

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,2)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1,2))
          
          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega abs(u-u_h) dx
          
          if (present(relementError)) then

            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                
                OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)
                
                ! L1-error is:   int_... abs(u-u_h) dx
                
                IELGlobal = p_IelementList(IELset+IEL-1)

                p_Derror(IELGlobal) = OM * abs(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))

                derror = derror + p_Derror(IELGlobal)
                
              end do ! ICUBP 
              
            end do ! IEL

          else
            
            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                
                OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)
                
                ! L1-error is:   int_... abs(u-u_h) dx
                
                derror = derror + &
                         OM * abs(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))
                
              end do ! ICUBP 
              
            end do ! IEL

          end if

        case (PPERR_H1ERROR)

          ! H1-error uses only 1st derivative of the function.

          if (present(ffunctionReference)) then          
            ! It's time to call our coefficient function to calculate the
            ! X-derivative values in the cubature points:  u(x,y)
            ! The result is saved in Dcoefficients(:,:,1)
            call ffunctionReference (DER_DERIV_X, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
                        
            ! Calculate the Y-derivative to Dcoefficients(:,:,2)

            call ffunctionReference (DER_DERIV_Y,rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1,2), rcollection)
          else
            Dcoefficients(:,1:IELmax-IELset+1,1:2) = 0.0_DP
          end if
          
          ! Calculate the X/Y-derivative of the FE function in the
          ! cubature points: u_h(x,y).
          ! Save the result to Dcoefficients(:,:,3..4)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV_X,&
                  Dcoefficients(:,1:IELmax-IELset+1,3))

          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV_Y,&
                  Dcoefficients(:,1:IELmax-IELset+1,4))

          ! Subtraction of Dcoefficients(:,:,1..2) from Dcoefficients(:,:,3..4) gives
          ! the error "grad(u-u_h)(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
          
          if (present(relementError)) then

            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
                
                IELGlobal = p_IelementList(IELset+IEL-1)

                p_Derror(IELGlobal) = OM * ((Dcoefficients(icubp,IEL,3)-Dcoefficients(icubp,IEL,1))**2 + &
                                            (Dcoefficients(icubp,IEL,4)-Dcoefficients(icubp,IEL,2))**2)

                derror = derror + p_Derror(IELGlobal)

              end do ! ICUBP 
              
            end do ! IEL

          else
            
            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
                
                derror = derror + OM * &
                         ((Dcoefficients(icubp,IEL,3)-Dcoefficients(icubp,IEL,1))**2 + &
                          (Dcoefficients(icubp,IEL,4)-Dcoefficients(icubp,IEL,2))**2)

              end do ! ICUBP 
              
            end do ! IEL

          end if
        
        case default
          call output_line('Unknown error type identifier!',&
              OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar2d_conf')
          call sys_halt()
        end select
        
        ! Release the temporary domain integration structure again
        call domint_doneIntegration (rintSubset)
    
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)
      deallocate(Domega)

    end do ! icurrentElementDistr

    ! derror is ||error||^2, so take the square root at last.
    if ((cerrortype .eq. PPERR_L2ERROR) .or.&
        (cerrortype .eq. PPERR_H1ERROR)) then
      derror = sqrt(derror)
      if (present(relementError)) then
        do IEL = 1, size(p_Derror,1)
          p_Derror(IEL) = sqrt(p_Derror(IEL))
        end do
      end if
    end if

  end subroutine pperr_scalar2d_conf

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalar3d_conf (rvectorScalar, cerrortype, derror,&
                                  rdiscretisation, ffunctionReference,&
                                  rcollection, relementError)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 3D version for double-precision vectors.
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(IN), target :: rvectorScalar
  
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(IN) :: cerrortype
  
  ! A discretisation structure specifying how to compute the error.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisation

  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_refFunctionSc.inc'
  optional :: ffunctionReference
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure to provide additional 
  ! information for callback routines.
  type(t_collection), intent(INOUT), optional :: rcollection

  ! OPTIONAL: A scalar vector which holds the calculated error per element
  type(t_vectorScalar), intent(INOUT), optional :: relementError
!</inputoutput>

!<output>
  ! Array receiving the calculated error.
  real(DP), intent(OUT) :: derror 
!</output>

!</subroutine>

    ! local variables
    integer :: icurrentElementDistr, ICUBP, NVE
    integer :: IEL, IELmax, IELset, IELGlobal
    real(DP) :: OM
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(:), allocatable :: Domega
    
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
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef

    ! Arrays for saving Jacobian determinants 
    real(DP), dimension(:,:), pointer :: p_Ddetj
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_domainIntSubset) :: rintSubset
    type(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial
  
    ! Type of transformation from the reference to the real element 
    integer :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! Pointer to the element-wise error
    real(DP), dimension(:), pointer :: p_Derror


    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDER
    ! according to these.

    Bder = .false.
    select case (cerrortype)
    case (PPERR_L1ERROR, PPERR_L2ERROR) 
      Bder(DER_FUNC3D) = .true.
    case (PPERR_H1ERROR) 
      Bder(DER_DERIV3D_X) = .true.
      Bder(DER_DERIV3D_Y) = .true.
      Bder(DER_DERIV3D_Z) = .true.
    case default
      call output_line('Unknown error type identifier!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar3d_conf')
      call sys_halt()
    end select
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPERR_NELEMSIM, p_rtriangulation%NEL)
    
    ! Set the current error to 0 and add the error contributions of each element
    ! to that.
    derror = 0.0_DP

    ! Set pointer to element-wise error
    if (present(relementError)) then
      call lsyssc_getbase_double(relementError, p_Derror)
    end if

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1, rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
    
      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Get the number of cubature points for the cubature formula
      ncubp = cub_igetNumPts(p_relementDistribution%ccubTypeEval)
      
      ! Allocate two arrays for the points and the weights
      allocate(Domega(ncubp))
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType), ncubp))
      
      ! Get the cubature formula
      call cub_getCubature(p_relementDistribution%ccubTypeEval, p_DcubPtsRef, Domega)
      
      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial, nelementsPerBlock))

      ! Allocate memory for the coefficients
      allocate(Dcoefficients(ncubp, nelementsPerBlock, 6))
    
      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      if (present(ffunctionReference)) then
        ! Evaluate real coordinates if not necessary.
        cevaluationTag = ior(cevaluationTag, EL_EVLTAG_REALPOINTS)
      end if
                      
      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag, EL_EVLTAG_DETJ)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL
    
      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, PPERR_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
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
        p_Ddetj => revalElementSet%p_Ddetj

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        select case (cerrortype)
        
        case (PPERR_L2ERROR)
          
          ! L2-error uses only the values of the function.
          
          if (present(ffunctionReference)) then
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x,y,z)
            ! The result is saved in Dcoefficients(:,:,1)
            call ffunctionReference (DER_FUNC3D, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
          else
            Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP
          end if

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x,y,z).
          ! Save the result to Dcoefficients(:,:,2)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC3D,&
                  Dcoefficients(:,1:IELmax-IELset+1,2))
          
          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (u-u_h,u-u_h) dx
          
          if (present(relementError)) then

            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
                
                IELGlobal = p_IelementList(IELset+IEL-1)

                p_Derror(IELGlobal) = OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2

                derror = derror + p_Derror(IELGlobal)

              end do ! ICUBP 
              
            end do ! IEL

          else
            
            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
                
                derror = derror + &
                         OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2
                
              end do ! ICUBP 
              
            end do ! IEL

          end if

        case (PPERR_L1ERROR)
          
          ! L1-error uses only the values of the function.
          
          if (present(ffunctionReference)) then
            ! It's time to call our coefficient function to calculate the
            ! function values in the cubature points:  u(x,y,z)
            ! The result is saved in Dcoefficients(:,:,1)
            call ffunctionReference (DER_FUNC3D, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
          else
            Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP
          end if

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x,y,z).
          ! Save the result to Dcoefficients(:,:,2)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC3D,&
                  Dcoefficients(:,1:IELmax-IELset+1,2))
          
          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2) gives
          ! the error "u-u_h(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega abs(u-u_h) dx

          if (present(relementError)) then

            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! L1-error is:   int_... abs(u-u_h) dx
                
                IELGlobal = p_IelementList(IELset+IEL-1)
              
                p_Derror(IELGlobal) = OM * abs(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))
                
                derror = derror + p_Derror(IELGlobal)

              end do ! ICUBP 
              
            end do ! IEL

          else
            
            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! L1-error is:   int_... abs(u-u_h) dx
                
                derror = derror + &
                         OM * abs(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))
                
              end do ! ICUBP 
              
            end do ! IEL
            
          end if

        case (PPERR_H1ERROR)

          ! H1-error uses only 1st derivative of the function.

          if (present(ffunctionReference)) then          
            ! It's time to call our coefficient function to calculate the
            ! X-derivative values in the cubature points:  u(x,y,z)
            ! The result is saved in Dcoefficients(:,:,1)
            call ffunctionReference (DER_DERIV3D_X, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
                        
            ! Calculate the Y-derivative to Dcoefficients(:,:,2)
            call ffunctionReference (DER_DERIV3D_Y, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset,&
                        Dcoefficients(:,1:IELmax-IELset+1,2), rcollection)

            ! Calculate the Z-derivative to Dcoefficients(:,:,3)
            call ffunctionReference (DER_DERIV3D_Z, rdiscretisation, &
                        int(IELmax-IELset+1), ncubp,&
                        revalElementSet%p_DpointsReal,&
                        IdofsTrial, rintSubset, &
                        Dcoefficients(:,1:IELmax-IELset+1,3), rcollection)
          else
            Dcoefficients(:,1:IELmax-IELset+1,1:3) = 0.0_DP
          end if
          
          ! Calculate the X/Y/Z-derivative of the FE function in the
          ! cubature points: u_h(x,y,z).
          ! Save the result to Dcoefficients(:,:,4..6)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV3D_X,&
                  Dcoefficients(:,1:IELmax-IELset+1,4))

          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV3D_Y,&
                  Dcoefficients(:,1:IELmax-IELset+1,5))

          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_DERIV3D_Z,&
                  Dcoefficients(:,1:IELmax-IELset+1,6))

          ! Subtraction of Dcoefficients(:,:,1..3) from Dcoefficients(:,:,4..6) gives
          ! the error "grad(u-u_h)(cubature pt.)"!
          !        
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx

          if (present(relementError)) then

            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
                IELGlobal = p_IelementList(IELset+IEL-1)

                p_Derror(IELGlobal) = OM * ((Dcoefficients(icubp,IEL,4)-Dcoefficients(icubp,IEL,1))**2 + &
                                            (Dcoefficients(icubp,IEL,5)-Dcoefficients(icubp,IEL,2))**2 + &
                                            (Dcoefficients(icubp,IEL,6)-Dcoefficients(icubp,IEL,3))**2)

                derror = derror + p_Derror(IELGlobal)

              end do ! ICUBP 
              
            end do ! IEL

          else
            
            do IEL=1,IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
                
                derror = derror + OM * &
                         ((Dcoefficients(icubp,IEL,4)-Dcoefficients(icubp,IEL,1))**2 + &
                          (Dcoefficients(icubp,IEL,5)-Dcoefficients(icubp,IEL,2))**2 + &
                          (Dcoefficients(icubp,IEL,6)-Dcoefficients(icubp,IEL,3))**2)
                
              end do ! ICUBP 
              
            end do ! IEL

          end if
        
        case default
          call output_line('Unknown error type identifier!',&
              OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar3d_conf')
          call sys_halt()
        end select
    
        ! Release again the temporary domain integration subset
        call domint_doneIntegration (rintSubset)
    
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)
      deallocate(Domega)

    end do ! icurrentElementDistr

    ! derror is ||error||^2, so take the square root at last.
    if ((cerrortype .eq. PPERR_L2ERROR) .or.&
        (cerrortype .eq. PPERR_H1ERROR)) then
      derror = sqrt(derror)
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarBoundary2D (cerrortype,ccubType,&
      derror,rboundaryRegion,rvectorScalar,ffunctionReference,&
      rcollection,rdiscretisation)

!<description>
  ! This routine calculates the error or the norm, respectively, of a given 
  ! finite element function in rvector to a given analytical 
  ! callback function ffunctionReference.
  !
  ! If ffunctionReference is specified, the routine calculates
  !   $$ ||y-z||_{L_2}  \textrm{ or }  ||y-z||_{L_1}  \textrm{ or }  ||y-z||_{H_1}$$
  ! with $y$=rvectorScalar and $z$=ffunctionReference.
  !
  ! If ffunctionReference is not specified, the routine calculates
  !   $$ ||y||_{L_2}  \textrm{ or }  ||y||_{L_1}  \textrm{ or }  ||y||_{H_1}.$$
  !
  ! If the vector rvectorScalar is not specified, it's assumed to be =0.
  !
  ! rboundaryRegion is a t_boundaryRegion object that allows to
  ! specify the boundary region where the error should be computed.
  ! If not specified, the error is computed over the whole boundary.
!</description>

!<input>
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(IN) :: cerrortype
  
  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer, intent(IN) :: ccubType
  
  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  type(t_boundaryRegion), intent(IN), optional :: rboundaryRegion
  
  ! OPTIONAL: The FE solution vector. Represents a scalar FE function.
  ! If omitted, the function is assumed to be constantly =1.
  type(t_vectorScalar), intent(IN), optional, target :: rvectorScalar
  
  ! OPTIONAL: A callback function that provides the analytical reference 
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_refFunctionScBoundary.inc'
  optional :: ffunctionReference
  
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information. 
  type(t_collection), intent(INOUT), target, optional :: rcollection

  ! OPTIONAL: A discretisation structure specifying how to compute the error.
  ! Must be specified if rvectorScalar is not specified as this
  ! describes the domain/triangulation/...
  type(t_spatialDiscretisation), intent(IN), target, optional :: rdiscretisation
!</input>

!<output>
  ! The calculated error.
  real(DP), intent(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    type(t_boundaryRegion) :: rboundaryReg
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    real(DP) :: dlocalError
    integer :: ibdc
    
    ! Get the correct discretisation structure and check if we can use it.
    if (present(rdiscretisation)) then
      p_rdiscretisation => rdiscretisation
      call lsyssc_checkDiscretisation (rvectorScalar,p_rdiscretisation)
    else
      p_rdiscretisation => rvectorScalar%p_rspatialdiscr
    end if
    
    if (.not. associated(p_rdiscretisation)) then
      call output_line('No discretisation structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarBoundary2D')
      call sys_halt()
    end if
    
    if (p_rdiscretisation%ndimension .ne. NDIM2D) then
      call output_line('Only 2D discretisations allowed.',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarBoundary2D')
      call sys_halt()
    end if

    ! The vector must be unsorted, otherwise we can't set up the vector.
    if (rvectorScalar%isortStrategy .gt. 0) then
      call output_line('Vector must be unsorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarBoundary2D')
      call sys_halt()
    end if
  
    ! If the boundary region is specified, call pperr_scalarBoundary2d_conf
    ! for that boundary region. Otherwise, call pperr_scalarBoundary2d_conf
    ! for all possible boundary regions and sum up the errors.
    if (present(rboundaryRegion)) then
      call pperr_scalarBoundary2d_conf (cerrortype,ccubType,derror,&
        rboundaryRegion,rvectorScalar,ffunctionReference,&
        p_rdiscretisation,rcollection)
    else
      derror = 0.0_DP
      ! Create a boundary region for each boundary component and call
      ! the calculation routine for that.
      do ibdc=1,boundary_igetNBoundComp(p_rdiscretisation%p_rboundary)
        call boundary_createRegion (p_rdiscretisation%p_rboundary, &
            ibdc, 0, rboundaryReg)
        call pperr_scalarBoundary2d_conf (cerrortype,ccubType,dlocalError,&
          rboundaryReg,rvectorScalar,ffunctionReference,&
          p_rdiscretisation,rcollection)
        derror = derror + dlocalError
      end do
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarBoundary2d_conf (cerrortype,ccubType,derror,&
      rboundaryRegion,rvectorScalar,ffunctionReference,&
      rdiscretisation,rcollection)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! Type of error to compute. Bitfield. This is a combination of the
  ! PPERR_xxxx-constants, which specifies what to compute.
  ! Example: PPERR_L2ERROR computes the $L_2$-error.
  integer, intent(IN) :: cerrortype
  
  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer, intent(IN) :: ccubType

  ! A t_boundaryRegion specifying the boundary region where
  ! to calculate. 
  type(t_boundaryRegion), intent(IN), optional :: rboundaryRegion

  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(IN), optional, target :: rvectorScalar
  
  ! A discretisation structure specifying how to compute the error.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisation
  
  ! Optional: A collection structure to provide additional 
  ! information for callback routines.
  type(t_collection), intent(INOUT), optional :: rcollection

  ! OPTIONAL: A callback function that provides a coefficient in front
  ! of the FE function. If not specified, a value of 1 is assumed.
  include 'intf_refFunctionScBoundary.inc'
  optional :: ffunctionReference
!</input>

!<output>
  ! Array receiving the calculated error.
  real(DP), intent(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    integer(I32), dimension(:), allocatable :: IelementOrientation
    integer, dimension(:), allocatable :: Ielements
    real(DP), dimension(:,:), allocatable :: DedgePosition
    
    integer :: ibdc,ibdcoffset,iedge,ilocaledge
    integer :: NEL,NELbdc,iel
    integer(I32) :: ctrafoType
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:), pointer :: p_IelementsAtBoundary
    real(DP), dimension(:), pointer :: p_DedgeParameterValue
    real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Arrays for cubature points
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    real(DP), dimension(:,:,:), allocatable :: Dxi2D,Dpoints,DpointsRef
    real(DP), dimension(CUB_MAXCUBP) :: Domega1D
    real(DP), dimension(:,:,:), allocatable :: Dvalues
    real(DP), dimension(NDIM2D,TRIA_MAXNVE) :: Dcoord
    integer :: ncubp,ipoint,ieltype
    integer(I32) :: icoordSystem
    real(DP) :: dlen,dpar1,dpar2
    
    ! Arrays for element distributions for every element
    integer(I32), dimension(:), pointer :: p_IelementDistr

    ! List of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

    ! Get some pointers and arrays for quicker access
    p_rtriangulation => rdiscretisation%p_rtriangulation
    p_RelementDistribution => rdiscretisation%RelementDistr
    
    call storage_getbase_int (p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_int (p_rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    call storage_getbase_int (p_rtriangulation%h_IelementsAtBoundary,&
        p_IelementsAtBoundary)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
        p_DvertexCoordinates)
    call storage_getbase_double (p_rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)
    call storage_getbase_double (p_rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)
        
    ! Boundary component?
    ibdc = rboundaryRegion%iboundCompIdx
    
    ! Number of elements on that boundary component?
    NELbdc = p_IboundaryCpIdx(ibdc+1)-p_IboundaryCpIdx(ibdc)
    
    ! Position of the boundary component?
    ibdcoffset = p_IboundaryCpIdx(ibdc)
        
    ! In a first step, we figure out the elements on the boundary and their
    ! orientation. Allocate arrays that are large enough to hold
    ! even all elements on the boundary if necessary.
    allocate(Ielements(NELbdc), IelementOrientation(NELbdc))
    
    ! Allocate an array saving the start- and end-parameter values
    ! of the edges on the boundary.
    allocate(DedgePosition(2,NELbdc))
    
    ! Loop through the edges on the boundary component ibdc.
    ! If the edge is inside, remember the element number and figure out
    ! the orientation of the edge.
    ! NEL counts the total number of elements in the region.
    NEL = 0
    do iedge = 1,NELbdc
      if (boundary_isInRegion(rboundaryRegion,ibdc,&
          p_DedgeParameterValue(iedge))) then
        NEL = NEL + 1
        
        ! Element number
        Ielements(NEL) = p_IelementsAtBoundary(iedge)
        
        ! Element orientation; i.e. the local number of the boundary edge 
        do ilocaledge = 1,ubound(p_IedgesAtElement,1)
          if (p_IedgesAtElement(ilocaledge,p_IelementsAtBoundary(iedge)) .eq. &
              p_IedgesAtBoundary(iedge)) exit
        end do
        IelementOrientation(NEL) = ilocaledge
        
        ! Save the start parameter value of the edge -- in length
        ! parametrisation.
        dpar1 = p_DvertexParameterValue(iedge)
        
        ! Save the end parameter value. Be careful: The last edge
        ! must be treated differently!
        if (iedge .ne. NELbdc) then
          dpar2 = p_DvertexParameterValue(iedge+1)
        else
          dpar2 = boundary_dgetMaxParVal(&
            rdiscretisation%p_rboundary,ibdc)
        end if
        
        DedgePosition(1,NEL) = &
          boundary_convertParameter(rdiscretisation%p_rboundary, &
            ibdc, dpar1, rboundaryRegion%cparType, BDR_PAR_LENGTH)
            
        DedgePosition(2,NEL) = &
          boundary_convertParameter(rdiscretisation%p_rboundary, &
            ibdc, dpar2, rboundaryRegion%cparType, BDR_PAR_LENGTH)
         
      end if
    end do
    
    ! Get the parameter values of the 1D cubature formula
    ! as well as the number of cubature points ncubp
    call cub_getCubPoints(ccubType, ncubp, Dxi1D, Domega1D)
    
    ! Map the 1D cubature points to the edges in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,NEL))
    if (rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then
      ! All elements have the same type. Get the type of the
      ! corresponding coordinate system and transform the coordinates.
      icoordSystem = elem_igetCoordSystem(&
        rdiscretisation%RelementDistr(1)%celement)
      do iel = 1,NEL
        call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(iel), &
            ncubp, Dxi1D, Dxi2D(:,:,iel))
      end do
    else
      ! The type of the coordinate system may change with every element.
      ! So we may have to switch... celements in the discretisation
      ! structure informs us about the element type.
      call storage_getbase_int (rdiscretisation%h_IelementDistr,&
          p_IelementDistr)
      do iel = 1,NEL
        ieltype = p_RelementDistribution(p_IelementDistr(Ielements(iel)))%celement
        icoordSystem = elem_igetCoordSystem(ieltype)
        call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(iel), &
            ncubp, Dxi1D, Dxi2D(:,:,iel))
      end do
    end if
    
    ! Transpose the coordinate array such that we get coordinates
    ! we can work with.
    allocate(DpointsRef(NDIM2D+1,ncubp,NEL))
    do iel=1,NEL
      DpointsRef(:,:,iel) = transpose(Dxi2D(:,:,iel))
    end do
    
    ! Dxi2D is not needed anymore.
    deallocate(Dxi2D)
    
    ! If the reference function exists, calculate the coordinates of the
    ! points on world coordinates
    if (present(ffunctionReference)) then
      
      ! We need the real coordinates of the points.
      allocate(Dpoints(ncubp,NDIM2D+1,NEL))
      
      if (rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then
        ! All elements with the same transformation
        ctrafoType = elem_igetTrafoType(&
            rdiscretisation%RelementDistr(1)%celement)
        do iel = 1,NEL

          ! Get the points forming the element
          do ipoint = 1,ubound(p_IverticesAtElement,1)
            Dcoord(1,ipoint) = &
                p_DvertexCoordinates(1,p_IverticesAtElement(ipoint,iel))
            Dcoord(2,ipoint) = &
                p_DvertexCoordinates(2,p_IverticesAtElement(ipoint,iel))
          end do

          ! Transform the cubature points
          do ipoint = 1,ncubp
            call trafo_calcRealCoords (ctrafoType,Dcoord,&
                DpointsRef(:,ipoint,iel),Dpoints(:,ipoint,iel))
          end do
        end do
      else
        ! Transformation can be different for all elements
        do iel = 1,NEL

          ! Get the points forming the element
          do ipoint = 1,ubound(p_IverticesAtElement,1)
            Dcoord(1,ipoint) = &
                p_DvertexCoordinates(1,p_IverticesAtElement(ipoint,iel))
            Dcoord(2,ipoint) = &
                p_DvertexCoordinates(2,p_IverticesAtElement(ipoint,iel))
          end do

          ! Transform the cubature points
          do ipoint = 1,ncubp
            ieltype = p_RelementDistribution(&
                p_IelementDistr(Ielements(iel)))%celement
            ctrafoType = elem_igetTrafoType(ieltype)
            call trafo_calcRealCoords (ctrafoType,Dcoord,&
                DpointsRef(:,ipoint,iel),Dpoints(:,ipoint,iel))
          end do
        end do
      end if
    
    end if    
    
    ! So Dxi2 defines the coordinates on the reference element for all
    ! elements. Generally said, we have to evaluate the elements in these
    ! points now. That can be done by using fevl_evaluate_mult.
    !
    ! Which type of integral is to calculate? H1 or L2 or L1?
    select case (cerrortype)
    case (PPERR_L2ERROR)
    
      allocate (Dvalues(ncubp,NEL,2))
      Dvalues = 0.0_DP
      
      ! If the FE function exists, evaluate it.
      if (present(rvectorScalar)) then
        do iel = 1,NEL
          ! Evaluate on the element, write results to Dvalues
          call fevl_evaluate_mult (DER_FUNC, Dvalues(:,iel,1), rvectorScalar, &
              Ielements(iel), DpointsRef(:,:,iel))
        end do
      end if

      ! If the reference function exists, evaluate it.
      if (present(ffunctionReference)) then
        
        ! Evaluate the reference function on the boundary
        call ffunctionReference (DER_FUNC,rdiscretisation, &
            DpointsRef,Dpoints, Ielements(1:NEL), Dvalues(:,:,2),rcollection)
            
      end if
      
      ! Linear combination to get the actual values in the cubature points.
      do iel = 1,NEL
        do ipoint = 1,ncubp
          Dvalues(ipoint,iel,1) = Dvalues(ipoint,iel,2)-Dvalues(ipoint,iel,1)
        end do
      end do
    
      ! Now, Dvalues1 contains in Dvalues1(:,:,1) the term
      ! "u(x,y)-u_h(x,y)" -- in every cubature point on every
      ! element. We are finally able to calculate the integral!
      ! That means, run over all the edges and sum up...
      ! (ok, if rvectorScalar is not specified, we have
      !  -u_h(x,y) in Dvalues1(:,:,1), but as we take the square,
      !  it doesn't matter if we have u_h or -u_h there!)
      
      derror = 0.0_DP
      do iel = 1,NEL
      
        ! Get the length of the edge. Let's use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit inverval [-1,1] has length 2.
        dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))
      
        do ipoint = 1,ncubp
          derror = derror + dlen * Domega1D(ipoint) * (Dvalues(ipoint,iel,1)**2)
        end do
      end do
      
      deallocate(Dvalues)

    case (PPERR_L1ERROR)
    
      allocate (Dvalues(ncubp,NEL,2))
      Dvalues = 0.0_DP
      
      ! If the FE function exists, evaluate it.
      if (present(rvectorScalar)) then
        do iel = 1,NEL
          ! Evaluate on the element, write results to Dvalues
          call fevl_evaluate_mult (DER_FUNC, Dvalues(:,iel,1), rvectorScalar, &
              Ielements(iel), DpointsRef(:,:,iel))
        end do
      end if

      ! If the reference function exists, evaluate it.
      if (present(ffunctionReference)) then
        
        ! Evaluate the reference function on the boundary
        call ffunctionReference (DER_FUNC,rdiscretisation, &
            DpointsRef,Dpoints, Ielements(1:NEL), Dvalues(:,:,2),rcollection)
            
      end if
      
      ! Linear combination to get the actual values in the cubature points.
      do iel = 1,NEL
        do ipoint = 1,ncubp
          Dvalues(ipoint,iel,1) = Dvalues(ipoint,iel,2)-Dvalues(ipoint,iel,1)
        end do
      end do
    
      ! Now, Dvalues1 contains in Dvalues1(:,:,1) the term
      ! "u(x,y)-u_h(x,y)" -- in every cubature point on every
      ! element. We are finally able to calculate the integral!
      ! That means, run over all the edges and sum up...
      ! (ok, if rvectorScalar is not specified, we have
      !  -u_h(x,y) in Dvalues1(:,:,1), but as we take the square,
      !  it doesn't matter if we have u_h or -u_h there!)
      
      derror = 0.0_DP
      do iel = 1,NEL
      
        ! Get the length of the edge. Let's use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit inverval [-1,1] has length 2.
        dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))
      
        do ipoint = 1,ncubp
          derror = derror + dlen * Domega1D(ipoint) * abs(Dvalues(ipoint,iel,1))
        end do
      end do
      
      deallocate(Dvalues)
      
    case (PPERR_H1ERROR)
    
      allocate (Dvalues(ncubp,NEL,4))
      Dvalues = 0.0_DP
      
      ! If the FE function exists, evaluate it.
      if (present(rvectorScalar)) then
        do iel = 1,NEL
          ! Evaluate on the element, write results to Dvalues.
          !
          ! X-derivative
          call fevl_evaluate_mult (DER_DERIV_X, Dvalues(:,iel,1), rvectorScalar, &
              Ielements(iel), DpointsRef(:,:,iel))
              
          ! Y-derivative
          call fevl_evaluate_mult (DER_DERIV_Y, Dvalues(:,iel,2), rvectorScalar, &
              Ielements(iel), DpointsRef(:,:,iel))
        end do
      end if

      ! If the reference function exists, evaluate it.
      if (present(ffunctionReference)) then

        ! Evaluate the reference function on the boundary
        !
        ! X-derivative
        call ffunctionReference (DER_DERIV_X,rdiscretisation, &
            DpointsRef,Dpoints, Ielements(1:NEL), Dvalues(:,:,3),rcollection)

        ! Y-derivative
        call ffunctionReference (DER_DERIV_Y,rdiscretisation, &
            DpointsRef,Dpoints, Ielements(1:NEL), Dvalues(:,:,4),rcollection)
            
      end if

      ! Linear combination to get the actual values in the cubature points.
      ! ||u-u_h||_H1 = int ( grad(u-u_h) * grad(u-u_h) )
      !              ~ sum grad_x(u-u_h)**2 + grad_y(u-u_h)
      do iel = 1,NEL
        do ipoint = 1,ncubp
          Dvalues(ipoint,iel,1) = (Dvalues(ipoint,iel,1)-Dvalues(ipoint,iel,3))**2 + &
                                  (Dvalues(ipoint,iel,2)-Dvalues(ipoint,iel,4))**2 
        end do
      end do
      
      ! Now, Dvalues1 contains in Dvalues1(:,:,1) the term
      ! "grad(u(x,y))-grad(u_h(x,y))" -- in every cubature point on every
      ! element. We are finally able to calculate the integral!
      ! That means, run over all the edges and sum up...
      
      derror = 0.0_DP
      do iel = 1,NEL
      
        ! Get the length of the edge. Let's use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit inverval [-1,1] has length 2.
        dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))
      
        do ipoint = 1,ncubp
          derror = derror + dlen * Domega1D(ipoint) * (Dvalues(ipoint,iel,1)**2)
        end do
      end do
      
      deallocate(Dvalues)  
      
    end select

    ! Release memory

    if (present(ffunctionReference)) deallocate(Dpoints)
      
    deallocate(DedgePosition)
    deallocate(Ielements, IelementOrientation)
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarErrorEstimate (rvector,rvectorRef,ctype,derror,&
                                        rdiscretisationRef,relementError)

!<description>
  ! This routine calculates the error of a given FE function in
  ! rvector and a reference vector given in rvectorRef. As an example,
  ! one can think of the consistent FE gradient and some recovered
  ! reference gradient, c.f. ZZ-technique.
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the 
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution vector
    type(t_vectorScalar), intent(IN), target :: rvector

    ! FE reference solution vector
    type(t_vectorScalar), intent(IN), target :: rvectorRef
            
    ! Type of error to compute. A PPERR_xxERROR constant.
    ! PPERR_L2ERROR computes the L2-error, PPERR_L1ERROR the L1-error.
    integer, intent(IN) :: ctype
            
    ! OPTIONAL: A discretisation structure specifying how to compute the error.
    ! If not specified, the discretisation structure in the reference gradient 
    ! is used. If specified, the discretisation structure must be 'compatible'
    ! to the two gradient vectors (concerning NEQ,...). pperr_gradient uses the
    ! cubature formula specifier of the linear form in rdiscretisation to 
    ! compute the integrals for the error.
    type(t_spatialDiscretisation), intent(IN), target, optional :: rdiscretisationRef
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated error on each element.
    type(t_vectorScalar), intent(INOUT), optional :: relementError
!</inputoutput>

!<output>
    ! The calculated error.
    real(DP), intent(OUT) :: derror 
!</output>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rvectorBlock,rvectorBlockRef
    type(t_blockDiscretisation) :: rDiscr,rDiscrRef

    ! Create block discretisations with one component
    if (associated(rvector%p_rspatialdiscr)) then
      call spdiscr_createBlockDiscrInd(rvector%p_rspatialdiscr, rDiscr)
    else
      call output_line('Vector does not provide a spatial discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarL2ErrorEstimate')
      call sys_halt()
    end if

    if (associated(rvectorRef%p_rspatialdiscr)) then
      call spdiscr_createBlockDiscrInd(rvectorRef%p_rspatialdiscr, rDiscrRef)
    else
      call output_line('Reference vector does not provide a spatial discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarL2ErrorEstimate')
      call sys_halt()
    end if

    ! Create block vectors with one block
    call lsysbl_createVecFromScalar(rvector, rvectorBlock, rDiscr)
    call lsysbl_createVecFromScalar(rvectorRef, rvectorBlockRef, rDiscrRef)

    ! Call block version
    call pperr_blockErrorEstimate(rvectorBlock, rvectorBlockRef,ctype,&
        derror, rdiscretisationRef,relementError)

    ! Release auxiliary block discretisations
    call spdiscr_releaseBlockDiscr(rDiscr)
    call spdiscr_releaseBlockDiscr(rDiscrRef)

    ! Release auxiliary block vectors
    call lsysbl_releaseVector(rvectorBlock)
    call lsysbl_releaseVector(rvectorBlockRef)

  end subroutine 

  !****************************************************************************

!<subroutine>

  subroutine pperr_blockErrorEstimate (rvector,rvectorRef,ctype,derror,&
                                       rdiscretisationRef,relementError)

!<description>
  ! This routine calculates the error of a given FE function in rvector
  ! and a reference vector given in rvectorRef. Both vectors must have the
  ! same number of blocks. As an example, one can think of the consistent
  ! FE gradient and some recovered reference gradient, c.f. ZZ-technique.
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the 
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution block vector
    type(t_vectorBlock), intent(IN), target :: rvector

    ! FE reference solution block vector
    type(t_vectorBlock), intent(IN), target :: rvectorRef
            
    ! Type of error to compute. A PPERR_xxERROR constant.
    ! PPERR_L2ERROR computes the L2-error.
    ! PPERR_L1ERROR computes the L1-error.
    ! PPERR_H1ERROR computes the H1-error.
    integer, intent(IN) :: ctype

    ! OPTIONAL: A discretisation structure specifying how to compute the error.
    ! If not specified, the discretisation structure in the reference gradient 
    ! is used. If specified, the discretisation structure must be 'compatible'
    ! to the two gradient vectors (concerning NEQ,...). pperr_gradient uses the
    ! cubature formula specifier of the linear form in rdiscretisation to 
    ! compute the integrals for the error.
    type(t_spatialDiscretisation), intent(IN), target, optional :: rdiscretisationRef
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated error on each element.
    type(t_vectorScalar), intent(INOUT), optional :: relementError
!</inputoutput>

!<output>
    ! The calculated error.
    real(DP), intent(OUT) :: derror 
!</output>
!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisationRef
    integer :: i,k,icurrentElementDistr,iblock,ICUBP,NVE
    integer :: IEL, IELmax, IELset,IELGlobal
    real(DP) :: OM,delementError

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
    integer :: indofTrial,indofTrialRef
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    
    ! Jacobian determinant in all points
    real(DP), dimension(:,:), pointer :: p_Ddetj

    ! Pointer to the element error
    real(DP), dimension(:), pointer :: p_DelementError
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution,p_relementDistributionRef
    
    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    
    ! Type of transformation from the reference to the real element 
    integer :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! Element evaluation set that collects element specific information
    ! for the evaluation on the cells.
    type(t_evalElementSet) :: revalElementSet
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial,IdofsTrialRef
    
    ! Get the correct discretisation structure for the solution vector
    p_rdiscretisation => rvector%p_rblockDiscr%RspatialDiscr(1)
    do iblock=2,rvector%nblocks
      call lsyssc_checkDiscretisation (rvector%RvectorBlock(iblock),p_rdiscretisation)
    end do

    ! Get the correct discretisation structure for the reference 
    ! vector and check if we can use it.
    if (present(rdiscretisationRef)) then
      p_rdiscretisationRef => rdiscretisationRef
      call lsyssc_checkDiscretisation (rvectorRef%RvectorBlock(1),p_rdiscretisationRef)
    else
      p_rdiscretisationRef => rvectorRef%p_rblockDiscr%RspatialDiscr(1)
    end if
    do iblock=2,rvectorRef%nblocks
      call lsyssc_checkDiscretisation (rvectorRef%RvectorBlock(iblock),p_rdiscretisationRef)
    end do
    
    if (.not. associated(p_rdiscretisation) .or.&
        .not. associated(p_rdiscretisationRef)) then
      call output_line('No discretisation structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
      call sys_halt()
    end if
    
    ! The vectors must have the same number of blocks
    if (rvector%nblocks .ne. rvectorRef%nblocks) then
      call output_line('Vectors have different number of blocks!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
      call sys_halt()
    end if
    
    ! The vector must be unsorted.
    do iblock=1,rvector%nblocks
      if (rvector%RvectorBlock(iblock)%isortStrategy    .gt. 0 .or.&
          rvectorRef%RvectorBlock(iblock)%isortStrategy .gt. 0) then
        call output_line('Vectors must be unsorted!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
        call sys_halt()
      end if
    end do

    ! We only need the function values of basis functions
    Bder = .false.
    Bder(DER_FUNC) = .true.

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscretisationRef%p_rtriangulation

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPERR_NELEMSIM,p_rtriangulation%NEL)

    ! Set the current error to 0 and add the error contributions of each element to that.
    derror = 0.0_DP

    ! Get a pointer to the element error (if required)
    if (present(relementError)) then
      call lsyssc_getbase_double(relementError,p_DelementError)
      call lalg_clearVectorDble (p_DelementError)
    end if

    ! Check that both discretisations have the same number of element distributions
    if (p_rdiscretisation%inumFESpaces .ne. &
        p_rdiscretisationRef%inumFESpaces) then
      call output_line('Number of element distributions mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
      call sys_halt()
    end if

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1,p_rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution    => p_rdiscretisation%RelementDistr(icurrentElementDistr)
      p_relementDistributionRef => p_rdiscretisationRef%RelementDistr(icurrentElementDistr)

      ! Check if element distrbutions have different number of elements
      if (p_relementDistribution%NEL .ne. &
          p_relementDistributionRef%NEL) then
        call output_line('Number of elements in distributions mismatch!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
        call sys_halt()
      end if

      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial    = elem_igetNDofLoc(p_relementDistribution%celement)
      indofTrialRef = elem_igetNDofLoc(p_relementDistributionRef%celement)

      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
     
      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      call cub_getCubPoints(p_relementDistributionRef%ccubTypeEval, ncubp, Dxi, Domega)

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistributionRef%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      do i=1,ncubp
        do k=1,ubound(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        end do
      end do
      
      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))
      allocate(IdofsTrialRef(indofTrialRef,nelementsPerBlock))

      ! Allocate memory for the coefficients, that is, two times the spatial dimension
      allocate(Dcoefficients(ncubp,nelementsPerBlock,2*rvector%nblocks))
    
      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
      cevaluationTag = ior(cevaluationTag,&
                      elem_getEvaluationTag(p_relementDistributionRef%celement))

      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)
    
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL

      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, PPERR_NELEMSIM
  
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
        call dof_locGlobMapping_mult(p_rdiscretisationRef, p_IelementList(IELset:IELmax), &
                                     IdofsTrialRef)
                                     
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj
            
        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! Depending on ctype, choose now the error to compute.
        select case (ctype)
        case (PPERR_L1ERROR)

          ! L1-error uses only the values of the function.

          ! Calculate the values of the FE solution vector and the reference solution 
          ! vector in the cubature points: u_h(x,y) and u_ref(x,y)
          ! Save the result to Dcoefficients(:,:,2*iblock-1) and 
          ! Dcoefficients(:,:,2*iblock)

          do iblock=1,rvector%nblocks
            
            ! solution vector
            call fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                    p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1,2*iblock))

            ! solution reference vector
            call fevl_evaluate_sim3 (rvectorRef%RvectorBlock(iblock), revalElementSet,&
                    p_relementDistributionRef%celement, IdofsTrialRef, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1,2*iblock-1))

          end do

          ! Subtraction of Dcoefficients(:,:,2*iblock-1) from Dcoefficients(:,:,2*iblock)
          ! and summing over all iblock=1,..,nblocks gives the error 
          ! $u_h(cubature pt.) - u_ref(cubature pt.)$
          
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (u_h-u_ref,u_h-u_ref) dx

          do IEL=1,IELmax-IELset+1

            ! Initialise element error by 0
            delementError = 0.0_DP

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp
              
              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!
              
              OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
              
              ! L1-error is:   int_... abs(u_h-u_ref) dx
              
              do iblock=1,rvector%nblocks
                delementError = delementError + &
                    OM * abs(Dcoefficients(icubp,IEL,2*iblock-1)-&
                            Dcoefficients(icubp,IEL,2*iblock))
              end do

            end do ! ICUBP 

            ! Apply to global error
            derror = derror + delementError

            ! Store in element error (if required)
            if (present(relementError)) then
              IELGlobal = p_IelementList(IELset+IEL-1)
              p_DelementError(IELGlobal) = delementError
            end if

          end do ! IEL
        
        case (PPERR_L2ERROR)
        
          ! L2-error uses only the values of the function.

          ! Calculate the values of the FE solution vector and the reference solution 
          ! vector in the cubature points: u_h(x,y) and u_ref(x,y)
          ! Save the result to Dcoefficients(:,:,2*iblock-1) and 
          ! Dcoefficients(:,:,2*iblock)

          do iblock=1,rvector%nblocks
            
            ! solution vector
            call fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                    p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1,2*iblock))

            ! solution reference vector
            call fevl_evaluate_sim3 (rvectorRef%RvectorBlock(iblock), revalElementSet,&
                    p_relementDistributionRef%celement, IdofsTrialRef, DER_FUNC,&
                    Dcoefficients(:,1:IELmax-IELset+1,2*iblock-1))

          end do

          ! Subtraction of Dcoefficients(:,:,2*iblock-1) from Dcoefficients(:,:,2*iblock)
          ! and summing over all iblock=1,..,nblocks gives the error 
          ! $u_h(cubature pt.) - u_ref(cubature pt.)$
          
          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (u_h-u_ref,u_h-u_ref) dx

          do IEL=1,IELmax-IELset+1

            ! Initialise element error by 0
            delementError = 0.0_DP

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp
              
              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!
              
              OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
              
              ! L2-error is:   int_... (u_h-u_ref)*(u_h-u_ref) dx
              
              do iblock=1,rvector%nblocks
                delementError = delementError + &
                    OM * (Dcoefficients(icubp,IEL,2*iblock-1)-&
                          Dcoefficients(icubp,IEL,2*iblock))**2
              end do

            end do ! ICUBP 
            
            ! Apply to global error
            derror = derror + delementError

            ! Store in element error (if required)
            if (present(relementError)) then
              IELGlobal = p_IelementList(IELset+IEL-1)
              p_DelementError(IELGlobal) = sqrt(delementError)
            end if

          end do ! IEL
          
        case DEFAULT
          
          call output_line('Requested error estimate not implemented!',&
              OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockErrorEstimate')
          call sys_halt()

        end select
        
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial,IdofsTrialRef)
      
    end do ! icurrentElementDistr

    if (ctype .ne. PPERR_L1ERROR) then
      ! derror is ||error||^2, so take the square root at last.
      derror = sqrt(derror)
    end if

  end subroutine 

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarStandardDeviation (rvector,ddeviation,relementDeviation)

!<description>
  ! This routine calculates the standard deviation
  ! 
  ! $$ \sigma=\sqrt{\int_\Omega r^2 u dx} $$
  !
  ! of a given FE function $u$ in rvector, whereby
  !
  ! $$ r^2=(x-\hat x)^2 + (y-\hat y)^2 + (z-\hat z)^2 $$
  !
  ! and each component is computed from the following relation
  !
  ! $$ \hat x=\int_\Omega x u dx $$
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the 
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution vector
    type(t_vectorScalar), intent(IN), target :: rvector
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated deviation on each element.
    type(t_vectorScalar), intent(INOUT), optional :: relementDeviation
!</inputoutput>

!<output>
    ! The calculated standard deviation.
    real(DP), intent(OUT) :: ddeviation
!</output>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rvectorBlock
    type(t_blockDiscretisation) :: rDiscr

    ! Create block discretisations with one component
    if (associated(rvector%p_rspatialdiscr)) then
      call spdiscr_createBlockDiscrInd(rvector%p_rspatialdiscr, rDiscr)
    else
      call output_line('Vector does not provide a spatial discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarStandardDeviation')
      call sys_halt()
    end if

    ! Create block vectors with one block
    call lsysbl_createVecFromScalar(rvector, rvectorBlock, rDiscr)

    ! Call block version
    call pperr_blockStandardDeviation(rvectorBlock, ddeviation, relementDeviation)

    ! Release auxiliary block discretisations
    call spdiscr_releaseBlockDiscr(rDiscr)

    ! Release auxiliary block vectors
    call lsysbl_releaseVector(rvectorBlock)

  end subroutine pperr_scalarStandardDeviation

  !****************************************************************************

!<subroutine>

  subroutine pperr_blockStandardDeviation (rvector,ddeviation,relementDeviation)

!<description>
  ! This routine calculates the standard deviation
  ! 
  ! $$ \sigma=\sqrt{\int_\Omega r^2 u dx} $$
  !
  ! of a given FE function $u$ in rvector, whereby
  !
  ! $$ r^2=(x-\hat x)^2 + (y-\hat y)^2 + (z-\hat z)^2 $$
  !
  ! and each component is computed from the following relation
  !
  ! $$ \hat x=\int_\Omega x u dx $$
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the 
  ! cubature formula to use for each element distribution.
!</description>

!<input>
    ! FE solution block vector
    type(t_vectorBlock), intent(IN), target :: rvector
!</input>

!<inputoutput>
    ! OPTIONAL: Scalar vector that stores the calculated deviation on each element.
    type(t_vectorScalar), intent(INOUT), optional :: relementDeviation
!</inputoutput>

!<output>
    ! The calculated deviation.
    real(DP), intent(OUT) :: ddeviation
!</output>
!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    integer :: i,k,icurrentElementDistr,iblock,ICUBP,NVE,idim
    integer :: IEL, IELmax, IELset,IELGlobal
    real(DP) :: OM,delementDeviation

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
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef

    ! An array receiving the coordinates of cubature points on
    ! the real element for all elements in a set.
    real(DP), dimension(:,:,:), pointer :: p_DcubPtsReal

    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    
    ! Pointer to the element deviation
    real(DP), dimension(:), pointer :: p_DelementDeviation
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution
    
    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients

    ! Mathematical expectation of the center of mass
    real(DP), dimension(NDIM3D) :: DmassCenter
    
    ! Type of transformation from the reference to the real element 
    integer :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! Element evaluation set that collects element specific information
    ! for the evaluation on the cells.
    type(t_evalElementSet) :: revalElementSet

    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial
    
    ! Get the correct discretisation structure for the solution vector
    p_rdiscretisation => rvector%p_rblockDiscr%RspatialDiscr(1)
    do iblock=2,rvector%nblocks
      call lsyssc_checkDiscretisation (rvector%RvectorBlock(iblock),p_rdiscretisation)
    end do

    if (.not. associated(p_rdiscretisation)) then
      call output_line('No discretisation structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockStandardDeviation')
      call sys_halt()
    end if

    ! The vector must be unsorted.
    do iblock=1,rvector%nblocks
      if (rvector%RvectorBlock(iblock)%isortStrategy .gt. 0) then
        call output_line('Vectors must be unsorted!',&
            OU_CLASS_ERROR,OU_MODE_STD,'pperr_blockStandardDeviation')
        call sys_halt()
      end if
    end do

    ! We only need the function values of basis functions
    Bder = .false.
    Bder(DER_FUNC) = .true.

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => p_rdiscretisation%p_rtriangulation

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPERR_NELEMSIM,p_rtriangulation%NEL)

    ! Set the mathematical expectation of the center of mass to 0
    DmassCenter = 0.0_DP

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1,p_rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => p_rdiscretisation%RelementDistr(icurrentElementDistr)

      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)

      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
     
      ! Initialise the cubature formula.
      ! Get cubature weights and point coordinates on the reference element
      call cub_getCubPoints(p_relementDistribution%ccubTypeEval,&
          ncubp, Dxi, Domega)

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      do i=1,ncubp
        do k=1,ubound(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        end do
      end do

      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients.
      allocate(Dcoefficients(ncubp,nelementsPerBlock,rvector%nblocks))
    
      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REALPOINTS)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL

      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, PPERR_NELEMSIM
  
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj
        p_DcubPtsReal => revalElementSet%p_DpointsReal
            
        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! Standard deviation uses only the values of the function.

        ! Calculate the values of the FE solution vector in the cubature 
        ! points u_h(x,y) and save the result to Dcoefficients(:,:,iblock)

        do iblock=1,rvector%nblocks
          
          ! solution vector
          call fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1,iblock))
                  
        end do

        ! Calculate the mathematical expectation of the center of mass
        ! $\hat x_h=\int_\Omega x u_h dx$

        ! Loop through elements in the set and for each element,
        ! loop through the DOF's and cubature points to calculate the
        ! integral: int_Omega x*u_h dx, for x,y and z
        
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp
            
            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that's normal!
            
            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
            
            ! Mathematical expectation of the center of mass is:
            ! int_... x*u_h dx
            
            do iblock=1,rvector%nblocks
              do idim=1,p_rdiscretisation%ndimension
                DmassCenter(idim) = DmassCenter(idim) + &
                    OM * p_DcubPtsReal(idim,icubp,IEL) * &
                         Dcoefficients(icubp,IEL,iblock)
              end do
            end do

          end do ! ICUBP 
        
        end do ! IEL
        
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)      
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)
      
    end do ! icurrentElementDistr

    ! Ok, we have the mathematical expectation of the center of mass.
    ! Let's compute the standard deviation.

    Ddeviation = 0.0_DP

    ! Get a pointer to the element deviation (if required)
    if (present(relementDeviation)) then
      call lsyssc_getbase_double(relementDeviation,p_DelementDeviation)
      call lalg_clearVectorDble (p_DelementDeviation)
    end if

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

    ! Make sure that we have determinants.
    cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1,p_rdiscretisation%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistribution => p_rdiscretisation%RelementDistr(icurrentElementDistr)

      ! Cancel if this element distribution is empty.
      if (p_relementDistribution%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)

      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistribution%celement)
     
      ! Initialise the cubature formula.
      ! Get cubature weights and point coordinates on the reference element
      call cub_getCubPoints(p_relementDistribution%ccubTypeEval,&
          ncubp, Dxi, Domega)

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

      ! Allocate some memory to hold the cubature points on the reference element
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      do i=1,ncubp
        do k=1,ubound(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        end do
      end do

      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))

      ! Allocate memory for the coefficients.
      allocate(Dcoefficients(ncubp,nelementsPerBlock,rvector%nblocks))
    
      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REALPOINTS)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistribution%NEL

      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, PPERR_NELEMSIM
  
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
      
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => revalElementSet%p_Ddetj
        p_DcubPtsReal => revalElementSet%p_DpointsReal

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! Standard deviation uses only the values of the function.

        ! Calculate the values of the FE solution vector in the cubature 
        ! points u_h(x,y) and save the result to Dcoefficients(:,:,iblock)

        do iblock=1,rvector%nblocks
          
          ! solution vector
          call fevl_evaluate_sim3 (rvector%RvectorBlock(iblock), revalElementSet,&
                  p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
                  Dcoefficients(:,1:IELmax-IELset+1,iblock))
                  
        end do

        ! Calculate the standard deviation
        ! $\int_\Omega ((x-\hat x_h)^2 + (y-\hat y_h)^2 + (z-\hat z_h)^2) u_h dx$

        ! Loop through elements in the set and for each element,
        ! loop through the DOF's and cubature points to calculate the
        ! integral: int_Omega (x-\hat x_h)*(x-\hat x_h)*u_h dx
        ! and sum up for x,y and z
        
        do IEL=1,IELmax-IELset+1

          ! Initialise element deviation by 0
          delementDeviation = 0.0_DP

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp
            
            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that's normal!
            
            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
            
            ! Standard deviation is: int_... (x-\hat x_h)*(x-\hat x_h)*u_h dx
            ! summed up for all x,y and z
            
            do iblock=1,rvector%nblocks
              do idim=1,p_rdiscretisation%ndimension
                delementDeviation = delementDeviation + &
                    abs(OM * Dcoefficients(icubp,IEL,iblock)) * &
                       (p_DcubPtsReal(idim,icubp,IEL)-DmassCenter(idim))**2
              end do
            end do

          end do ! ICUBP 
        
          ! Apply to global deviation
          ddeviation = ddeviation + delementDeviation

          ! Store in element deviation (if required)
          if (present(relementDeviation)) then
            IELGlobal = p_IelementList(IELset+IEL-1)
            p_DelementDeviation(IELGlobal) = sqrt(delementDeviation)
          end if

        end do ! IEL
        
      end do ! IELset
      
      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)      
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)
      
    end do ! icurrentElementDistr
    
    ! ddeviation is ||deviation||^2, so take the square root at last.
    if (ddeviation .ge. huge(ddeviation)) then
      ddeviation = 0._DP
    else
      ddeviation = sqrt(ddeviation)
    end if

  end subroutine pperr_blockStandardDeviation

  !****************************************************************************

!<subroutine>

  subroutine pperr_scalarTargetFunc (rvectorScalar, derror, ftargetFuncReference,&
                                     rcollection, rdiscretisation, rerror)

!<description>

  ! This routine calculates the local weighted mean value target
  ! functional of a given finite element function in rvectorScalar to
  ! a given analytical callback function ftargetFuncReference using
  ! the weighting function ftargetFuncWeight.
  !
  ! $$ \int_\Omega w(x)*( u(x)-u_h(x) ) dx $$
  !
  ! Note: For the evaluation of the integrals, ccubTypeEval from the
  ! element distributions in the discretisation structure specifies the
  ! cubature formula to use for each element distribution. 
!</description>

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorScalar), intent(IN), target :: rvectorScalar
  
  ! OPTIONAL: A callback function that provides the analytical
  !           reference function to which the error should be
  !           computed.  If not specified, the reference function is
  !           assumed to be zero!
  include 'intf_refTargetFunctionSc.inc'
  optional :: ftargetFuncReference
  
  ! OPTIONAL: A collection structure. This structure is given to the
  !           callback function to provide additional information.
  type(t_collection), intent(INOUT), target, optional :: rcollection
  
  ! OPTIONAL: A discretisation structure specifying how to compute the
  !           error. If not specified, the discretisation structure
  !           in the vector is used.  If specified, the discretisation
  !           structure must be 'compatible' to the vector (concerning
  !           NEQ,...). pperr_scalar uses the cubature formula
  !           specifier of the linear form in rdiscretisation to
  !           compute the integrals for the error.
  type(t_spatialDiscretisation), intent(IN), target, optional :: rdiscretisation
!</input>

!<inputoutput>
  ! OPTIONAL: A scalar vector which holds the calculated error per element
  type(t_vectorScalar), intent(INOUT), optional :: rerror
!</inputoutput>

!<output>
  ! The calculated error
  real(DP), intent(OUT) :: derror
!</output>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation

    ! Get the correct discretisation structure and check if we can use it.
    if (present(rdiscretisation)) then
      p_rdiscretisation => rdiscretisation
      call lsyssc_checkDiscretisation (rvectorScalar,p_rdiscretisation)
    else
      p_rdiscretisation => rvectorScalar%p_rspatialdiscr
    end if
    
    if (.not. associated(p_rdiscretisation)) then
      call output_line('No discretisation structure!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarTargetFunc')
      call sys_halt()
    end if
    
    ! The vector must be unsorted, otherwise we can't set up the vector.
    if (rvectorScalar%isortStrategy .gt. 0) then
      call output_line('Vector must be unsorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarTargetFunc')
      call sys_halt()
    end if

    ! Do we have a uniform triangulation? Would simplify a lot...
    if ((p_rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) .or.&
        (p_rdiscretisation%ccomplexity .eq. SPDISC_CONFORMAL)) then 

      select case(rvectorScalar%cdataType)
        
      case (ST_DOUBLE)
        call scalarTargetFunc_conf(rvectorScalar, derror,&
                                   p_rdiscretisation, ftargetFuncReference,&
                                   rcollection, rerror)

      case DEFAULT
        call output_line('Single precision vectors currently not supported!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarTargetFunc')
        call sys_halt()
      end select
    
    else
      call output_line('General discretisation not implemented!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalarTargetFunc')
      call sys_halt()
    end if

  contains

    ! Here, the real working routines follow

    !**************************************************************

    subroutine scalarTargetFunc_conf(rvectorScalar, derror,&
                                     rdiscretisation, ftargetFuncReference,&
                                     rcollection, rerror)

      type(t_vectorScalar), intent(IN), target :: rvectorScalar
      type(t_spatialDiscretisation), intent(IN), target :: rdiscretisation
      type(t_collection), intent(INOUT), optional :: rcollection
      type(t_vectorScalar), intent(INOUT), optional :: rerror
      real(DP), intent(OUT) :: derror

      include 'intf_refTargetFunctionSc.inc'
      optional :: ftargetFuncReference

      ! local variables
      integer :: icurrentElementDistr, ICUBP, NVE
      integer :: IEL, IELmax, IELset, IELGlobal
      real(DP) :: OM
      
      ! Array to tell the element which derivatives to calculate
      logical, dimension(EL_MAXNDER) :: Bder
      
      ! For every cubature point on the reference element,
      ! the corresponding cubature weight
      real(DP), dimension(:), allocatable :: Domega
      
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
      
      ! Arrays for saving Jacobian determinants and matrices
      real(DP), dimension(:,:), pointer :: p_Ddetj
      
      ! Current element distribution
      type(t_elementDistribution), pointer :: p_relementDistribution
      
      ! Number of elements in the current element distribution
      integer :: NEL
      
      ! Pointer to the values of the function that are computed by the callback routine.
      real(DP), dimension(:,:,:), allocatable :: Dcoefficients
      
      ! Number of elements in a block. Normally =BILF_NELEMSIM,
      ! except if there are less elements in the discretisation.
      integer :: nelementsPerBlock
      
      ! A t_domainIntSubset structure that is used for storing information
      ! and passing it to callback routines.
      type(t_domainIntSubset) :: rintSubset
      type(t_evalElementSet) :: revalElementSet
      
      ! Type of transformation from the reference to the real element 
      integer :: ctrafoType
      
      ! Element evaluation tag; collects some information necessary for evaluating
      ! the elements.
      integer(I32) :: cevaluationTag
      
      ! An allocateable array accepting the DOF's of a set of elements.
      integer, dimension(:,:), allocatable, target :: IdofsTrial
      
      ! Pointer to the element-wise target functional
      real(DP), dimension(:), pointer :: p_Derror

      
      ! We only need function values
      Bder = .false.
      Bder(DER_FUNC) = .true.

      ! Get a pointer to the triangulation - for easier access.
      p_rtriangulation => rdiscretisation%p_rtriangulation
      
      ! For saving some memory in smaller discretisations, we calculate
      ! the number of elements per block. For smaller triangulations,
      ! this is NEL. If there are too many elements, it's at most
      ! BILF_NELEMSIM. This is only used for allocating some arrays.
      nelementsPerBlock = min(PPERR_NELEMSIM,p_rtriangulation%NEL)
      
      ! Set the current target functional to 0 
      derror = 0.0_DP

      ! Set pointer to element-wise target functional
      if (present(rerror)) then
        call lsyssc_getbase_double(rerror, p_Derror)
      end if


      ! Now loop over the different element distributions
      ! (=combinations of trial and test functions) in the
      ! discretisation.
      
      do icurrentElementDistr = 1, rdiscretisation%inumFESpaces
        
        ! Activate the current element distribution
        p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
        
        ! Cancel if this element distribution is empty.
        if (p_relementDistribution%NEL .eq. 0) cycle
        
        ! Get the number of local DOF's for trial functions
        indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
        
        ! Get the number of corner vertices of the element
        NVE = elem_igetNVE(p_relementDistribution%celement)
        
        ! Get from the trial element space the type of coordinate system
        ! that is used there:
        ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
        
        ! Get the number of cubature points for the cubature formula
        ncubp = cub_igetNumPts(p_relementDistribution%ccubTypeEval)
        
        ! Allocate two arrays for the points and the weights
        allocate(Domega(ncubp))
        allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType), ncubp))
        
        ! Get the cubature formula
        call cub_getCubature(p_relementDistribution%ccubTypeEval, p_DcubPtsRef, Domega)
        
        ! Allocate memory for the DOF's of all the elements.
        allocate(IdofsTrial(indofTrial, nelementsPerBlock))
        
        ! Allocate memory for the coefficients
        allocate(Dcoefficients(ncubp,nelementsPerBlock, 3))
        
        ! Initialisation of the element set.
        call elprep_init(revalElementSet)
        
        ! Get the element evaluation tag of all FE spaces. We need it to evaluate
        ! the elements later. All of them can be combined with OR, what will give
        ! a combined evaluation tag. 
        cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
      
        ! Evaluate real coordinates if necessary.
        cevaluationTag = ior(cevaluationTag, EL_EVLTAG_REALPOINTS)

        ! Make sure that we have determinants.
        cevaluationTag = ior(cevaluationTag, EL_EVLTAG_DETJ)

        ! p_IelementList must point to our set of elements in the discretisation
        ! with that combination of trial functions
        call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                  p_IelementList)
        
        ! Get the number of elements there.
        NEL = p_relementDistribution%NEL
        
        ! Loop over the elements - blockwise.
        do IELset = 1, NEL, PPERR_NELEMSIM
          
          ! We always handle LINF_NELEMSIM elements simultaneously.
          ! How many elements have we actually here?
          ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
          ! elements simultaneously.
          
          IELmax = min(NEL, IELset-1+PPERR_NELEMSIM)
          
          ! Calculate the global DOF's into IdofsTrial.
          !
          ! More exactly, we call dof_locGlobMapping_mult to calculate all the
          ! global DOF's of our LINF_NELEMSIM elements simultaneously.
          call dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
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
          p_Ddetj => revalElementSet%p_Ddetj
          
          ! In the next loop, we don't have to evaluate the coordinates
          ! on the reference elements anymore.
          cevaluationTag = iand(cevaluationTag, not(EL_EVLTAG_REFPOINTS))

          ! Calculate the values of the FE function in the
          ! cubature points: u_h(x).
          ! Save the result to Dcoefficients(:,:,1)
          
          call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                                   p_relementDistribution%celement,&
                                   IdofsTrial, DER_FUNC,&
                                   Dcoefficients(:,1:IELmax-IELset+1,1))
          

          if (present(ftargetFuncReference)) then
            ! Calculate the values of the coefficient function in the
            ! cubature points: u(x).
            ! Save the result to Dcoefficients(:,:,2)
            
            call ftargetFuncReference(1, rdiscretisation, int(IELmax-IELset+1), ncubp, &
                                      revalElementSet%p_DpointsReal,&
                                      IdofsTrial, rintSubset, &
                                      Dcoefficients(:,1:IELmax-IELset+1,2), rcollection)

            ! Calculate the values of the weighting function in the
            ! cubature points: w(x).
            ! Save the result to Dcoefficients(:,:,3)

            call ftargetFuncReference(0, rdiscretisation, int(IELmax-IELset+1), ncubp, &
                                      revalElementSet%p_DpointsReal,&
                                      IdofsTrial, rintSubset, &
                                      Dcoefficients(:,1:IELmax-IELset+1,3), rcollection)
          else
            Dcoefficients(:,1:IELmax-IELset+1,2) = 0.0_DP
            Dcoefficients(:,1:IELmax-IELset+1,3) = 1.0_DP
          end if

          ! Subtract Dcoefficients(:,:,2) from Dcoefficients(:,:,1)
          ! and multiply the result by Dcoefficients(:,:,3) to obtain
          ! "w*(u_h-u) (cubature pt.)"!

          if (present(rerror)) then
            
            do IEL = 1, IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                IELGlobal = p_IelementList(IELset+IEL-1)
                
                p_Derror(IELGlobal) = OM * ( Dcoefficients(icubp,IEL,1) - &
                                             Dcoefficients(icubp,IEL,2)   ) * &
                                             Dcoefficients(icubp,IEL,3)

                derror = derror + p_Derror(IELGlobal)
                
              end do ! ICUBP 
              
            end do ! IEL
            
          else
            
            do IEL = 1, IELmax-IELset+1
              
              ! Loop over all cubature points on the current element
              do icubp = 1, ncubp
                
                ! calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Take the absolut value of the determinant of the mapping.
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that's normal!
                
                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))
                
                derror = derror + &
                         OM * ( Dcoefficients(icubp,IEL,1) - &
                                Dcoefficients(icubp,IEL,2)   ) * &
                                Dcoefficients(icubp,IEL,3)
                
              end do ! ICUBP 
              
            end do ! IEL

          end if

          ! Release the temporary domain integration structure again
          call domint_doneIntegration (rintSubset)
          
        end do ! IELset
        
        ! Release memory
        call elprep_releaseElementSet(revalElementSet)
        
        deallocate(p_DcubPtsRef)
        deallocate(Dcoefficients)
        deallocate(IdofsTrial)
        deallocate(Domega)
        
      end do ! icurrentElementDistr

    end subroutine scalarTargetFunc_conf

  end subroutine pperr_scalarTargetFunc
end module
