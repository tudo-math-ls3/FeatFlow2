!##############################################################################
!# ****************************************************************************
!# <name> trilinearformevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the discretisation of trilinear forms,
!# i.e. the creation of matrices and matrix structures. It contains the
!# following set of routines:
!#
!# 1.) trilf_buildMatrixScalar
!#     -> Assembles the entries of a matrix, which structure was build
!#        with bilf_createMatrixStructure before.
!#
!# A 'trilinear form' is in our context a bilinear form whose coefficient
!# function may depend on a Finite Element function:
!# <tex>
!#   $$ a(u,phi_i,psi_j)  =  \int c(x,y) f(u) g(\psi_j) h(\phi_i) $$
!# </tex>
!# with f,g,h derivative operators. Therefore, there is no
!# special 'create structure' routine; the structure is generated based on the
!# underlying bilinear form, which steps from the idea of treating u as
!# a variable coefficient:
!# <tex>
!#   $$ a_u(phi_i,psi_j)  =  \int c_u(x,y) g(\psi_j) h(\phi_i), $$
!#
!#   $$ c_u = c*f(u) $$
!# </tex>
!# </purpose>
!##############################################################################

module trilinearformevaluation

  use fsystem
  use genoutput
  use storage
  use linearsystemscalar
  use spatialdiscretisation
  use scalarpde
  use derivatives
  use cubature
  use domainintegration
  use element
  use elementpreprocessing
  use linearalgebra
  use transformation
  use triangulation
  use dofmapping
  use collection, only: t_collection
  
  implicit none
  
  private

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of elements to handle simultaneously when building matrices
#ifdef TRILF_NELEMSIM
  integer, parameter, public :: TRILF_NELEMSIM   = TRILF_NELEMSIM
#else
  integer, parameter, public :: TRILF_NELEMSIM   = 1000
#endif
  
!</constantblock>
!</constants>

  public :: trilf_buildMatrixScalar

contains

  !****************************************************************************

!<subroutine>

  subroutine trilf_buildMatrixScalar (rform,bclear,rmatrixScalar,rvector,&
                                      fcoeff_buildTrilMatrixSc_sim,rcollection)
  
!<description>
  ! This routine calculates the entries of a finite element matrix using
  ! a trilinear form
  ! <tex>
  !     $$ a(u,phi_i,psi_j)  =  \int c(x,y) f(u) g(\psi_j) h(\phi_i) $$
  ! </tex>
  ! The function $u$ is specified in rvector. The derivative quantifier
  ! rform(Idescriptors(1,.) specifies f(.), i.e. whether to take function 
  ! values, derivatives or what else to get the function value f(u).
  !
  ! The matrix structure must be prepared with bilf_createMatrixStructure
  ! in advance.
  ! In case the array for the matrix entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries.
  !
  ! For setting up the entries, the discretisation structure attached to
  ! the matrix is used (rmatrixScalar%p_rdiscretisation). This is
  ! normally attached to the matrix by bilf_createMatrixStructure.
  !
  ! For the routine to work properly, it is important that the discretisation
  ! structure of rvector is 'compatible' with the discretisation structure
  ! of the matrix! I.e., the element distributions must be the same 
  ! (in number and ordering of the elements) except for the element type!
  !
  ! The matrix must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
!</description>

!<input>

  ! The trilinear form specifying the underlying PDE of the discretisation.
  type(t_trilinearForm), intent(in) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! A finite element function $u$ to be used as multiplier in the trilinear form.
  ! Must be a double precision vector.
  type(t_vectorScalar), intent(in) :: rvector

  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  include 'intf_coefficientTrilMatrixSc.inc'
  optional :: fcoeff_buildTrilMatrixSc_sim
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! The matrix must be unsorted, otherwise we can not set up the matrix.
  ! Note that we cannot switch off the sorting as easy as in the case
  ! of a vector, since there is a structure behind the matrix! So the caller
  ! has to make sure, the matrix is unsorted when this routine is called.
  if (rmatrixScalar%isortStrategy .gt. 0) then
    call output_line('Matrix-structure must be unsorted!',&
        OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrixScalar')
    call sys_halt()
  end if

  if ((.not. associated(rmatrixScalar%p_rspatialDiscrTest)) .or. &
      (.not. associated(rmatrixScalar%p_rspatialDiscrTrial))) then
    call output_line('No discretisation associated!',&
        OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrixScalar')
    call sys_halt()
  end if

  ! Do we have a uniform triangulation? Would simplify a lot...
  select case (rmatrixScalar%p_rspatialDiscrTest%ccomplexity)
  case (SPDISC_UNIFORM) 
    ! Uniform discretisation; only one type of elements, e.g. P1 or Q1
    select case (rmatrixScalar%cdataType)
    case (ST_DOUBLE) 
      ! Which matrix structure do we have?
      select case (rmatrixScalar%cmatrixFormat) 
      case (LSYSSC_MATRIX9)
        !IF (PRESENT(fcoeff_buildMatrixSc_sim)) THEN
          call trilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,rvector,&  
                                          fcoeff_buildTrilMatrixSc_sim,rcollection)
      case (LSYSSC_MATRIX7)
        ! Convert structure 7 to structure 9.
        call lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX9)
        
        ! Create the matrix in structure 9
        call trilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,rvector,&  
                                        fcoeff_buildTrilMatrixSc_sim,rcollection)
                                       
        ! Convert back to structure 7
        call lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)
                                       
      case DEFAULT
        call output_line('Not supported matrix structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrixScalar')
        call sys_halt()
      end select
    case DEFAULT
      call output_line('Single precision matrices currently not supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrixScalar')
      call sys_halt()
    end select
    
  case (SPDISC_CONFORMAL) 
    
    ! Conformal discretisation; may have mixed P1/Q1 elements e.g.
    select case (rmatrixScalar%cdataType)
    case (ST_DOUBLE) 
      ! Which matrix structure do we have?
      select case (rmatrixScalar%cmatrixFormat) 
      case (LSYSSC_MATRIX9)
        !IF (PRESENT(fcoeff_buildMatrixSc_sim)) THEN
          call trilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,rvector,&  
                                          fcoeff_buildTrilMatrixSc_sim,rcollection)
        
      case (LSYSSC_MATRIX7)
        ! Convert structure 7 to structure 9
        call lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX9)
        
        ! Create the matrix in structure 9
        call trilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,rvector,&  
                                        fcoeff_buildTrilMatrixSc_sim,rcollection)
                                       
        ! Convert back to structure 7
        call lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)

      case DEFAULT
        call output_line('Not supported matrix structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrixScalar')
        call sys_halt()
      end select
    case DEFAULT
      call output_line('Single precision matrices currently not supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrixScalar')
      call sys_halt()
    end select
  case DEFAULT
    call output_line('General discretisation not implemented!',&
        OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrixScalar')
    call sys_halt()
  end select

  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine trilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,rvector,&
                                        fcoeff_buildTrilMatrixSc_sim,rcollection)
  
!<description>
  ! This routine calculates the entries of a finite element matrix.
  ! The matrix structure must be prepared with bilf_createMatrixStructure
  ! in advance. The discretisation is assumed to be conformal, i.e. the DOF`s
  ! of all finite elements must 'match'. Trial and test functions may be
  ! different.
  ! In case the array for the matrix entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries.
  !
  ! For setting up the entries, the discretisation structure attached to
  ! the matrix is used (rmatrixScalar%p_rdiscretisation). This is
  ! normally attached to the matrix by bilf_createMatrixStructure.
  !
  ! For the routine to work properly, it is important that the discretisation
  ! structure of rvector is 'compatible' with the discretisation structure
  ! of the matrix! I.e., the element distributions must be the same 
  ! (in number and ordering of the elements) except for the element type!
  !
  ! Double-precision version.
!</description>

!<input>
  ! The trilinear form specifying the underlying PDE of the discretisation.
  type(t_trilinearForm), intent(in) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! A finite element function $u$ to be used as multiplier in the trilinear form.
  ! Must be a double precision vector.
  type(t_vectorScalar), intent(in) :: rvector

  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  include 'intf_coefficientTrilMatrixSc.inc'
  optional :: fcoeff_buildTrilMatrixSc_sim
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variablesvboxdrv
  integer :: i,i1,icurrentElementDistr,JDFG, ICUBP, IALBET, IA, IB, ifunc
  logical :: bIdenticalTrialAndTest, bIdenticalFuncAndTrial, bIdenticalFuncAndTest
  integer :: IEL, IELmax, IELset, IDOFE, JDOFE
  integer :: JCOL0,JCOL,idertype
  real(DP) :: OM,AUX, DB
  
  ! Array to tell the element which derivatives to calculate
  logical, dimension(EL_MAXNDER) :: BderTrialTempl, BderTestTempl, BderFuncTempl
  logical, dimension(EL_MAXNDER) :: BderTrial, BderTest, BderFunc

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(:), allocatable :: Domega
  
  ! number of cubature points on the reference element
  integer :: ncubp
  
  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_KLD, p_KCOL
  real(DP), dimension(:), pointer :: p_DA
  
  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsTest, IdofsTrial
  integer, dimension(:,:), allocatable, target :: IdofsFunc
  integer, dimension(:,:), pointer :: p_IdofsTrial,p_IdofsFunc
  !INTEGER, DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
  !INTEGER, DIMENSION(:,:), POINTER :: p_IdofsTrial
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: DbasFunc,DbasTest,DbasTrial
  real(DP), dimension(:,:,:,:), pointer :: p_DbasFunc,p_DbasTrial
  
  ! Number of entries in the matrix - for quicker access
  integer :: NA,NVE
  integer :: NEQ
  
  ! Type of transformation from the reference to the real element 
  integer(I32) :: ctrafoType
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag
  
  ! Number of local degees of freedom for trial and test functions
  integer :: indofFunc, indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList
  
  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:), allocatable :: Dentry
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

  ! Pointer to the jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj
  
  ! Current element distribution for discretisation and function $u$.
  type(t_elementDistribution), pointer :: p_elementDistrTrial
  type(t_elementDistribution), pointer :: p_elementDistrTest
  type(t_elementDistribution), pointer :: p_elementDistrFunc
  
  ! Number of elements in the current element distribution
  integer :: NEL
  
  ! DOF-Data of the vector
  real(DP), dimension(:), pointer :: p_Ddata
  
  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock
  
  ! Some variables to support nonconstant coefficients in the matrix.
  
  ! Pointer to the coefficients that are computed by the callback routine.
  real(DP), dimension(:,:,:), allocatable :: Dcoefficients
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  type(t_domainIntSubset) :: rintSubset
  type(t_evalElementSet) :: revalElementSet
  logical :: bcubPtsInitialised
  
  ! The discretisation - for easier access
  type(t_spatialDiscretisation), pointer :: p_rdiscrTest,p_rdiscrTrial
  type(t_spatialDiscretisation), pointer :: p_rdiscrFunc
  
  if ((.not. associated(rmatrixScalar%p_rspatialDiscrTest)) .or. &
      (.not. associated(rmatrixScalar%p_rspatialDiscrTrial))) then
    call output_line('No discretisation associated!',&
        OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
    call sys_halt()
  end if

  ! Which derivatives of basis functions are needed?
  ! Check the descriptors of the bilinear form and set BDERxxxx
  ! according to these.

  BderTrialTempl = .false.
  BderTestTempl = .false.
  BderFuncTempl = .false.
  
  ! Loop through the additive terms
  do i=1,rform%itermCount
    ! The desriptor Idescriptors gives directly the derivative
    ! which is to be computed! Build templates for BDER.
    ! We do not compute the actual BDER here, as there might be some special
    ! processing if trial/test functions are identical!
    !
    ! At first build the descriptors for the trial functions
    I1=rform%Idescriptors(1,I)
    
    if ((I1 .lt.0) .or. (I1 .gt. DER_MAXNDER)) then
      call output_line('Invalid descriptor!',&
          OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
      call sys_halt()
    endif
    
    if (I1 .ne. 0) BderFuncTempl(I1)=.true.

    I1=rform%Idescriptors(2,I)
    
    if ((I1 .le.0) .or. (I1 .gt. DER_MAXNDER)) then
      call output_line('Invalid descriptor!',&
          OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
      call sys_halt()
    endif
    
    BderTrialTempl(I1)=.true.

    ! Then those of the test functions
    I1=rform%Idescriptors(3,I)
    
    if ((I1 .le.0) .or. (I1 .gt. DER_MAXNDER)) then
      call output_line('Invalid descriptor!',&
          OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
      call sys_halt()
    endif
    
    BderTestTempl(I1)=.true.
  end do
  
  ! Get information about the matrix:
  NA = rmatrixScalar%NA
  NEQ = rmatrixScalar%NEQ
  
  ! We need KCOL/KLD of our matric
  call lsyssc_getbase_Kcol (rmatrixScalar,p_KCOL)
  call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
  
  ! Check if the matrix entries exist. If not, allocate the matrix.
  if (rmatrixScalar%h_DA .eq. ST_NOHANDLE) then

    ! Clear the entries in the matrix - we need to start with zero
    ! when assembling a new matrix!
    call storage_new ('trilf_buildMatrix9d_conf2', 'DA', &
                        NA, ST_DOUBLE, rmatrixScalar%h_DA, &
                        ST_NEWBLOCK_ZERO)
    call lsyssc_getbase_double (rmatrixScalar,p_DA)

  else
  
    call lsyssc_getbase_double (rmatrixScalar,p_DA)

    ! If desired, clear the matrix before assembling.
    if (bclear) then
      call lalg_clearVectorDble (p_DA)
    end if
    
  end if
  
  ! Get the discretisation(s)
  p_rdiscrTest => rmatrixScalar%p_rspatialDiscrTest
  p_rdiscrTrial => rmatrixScalar%p_rspatialDiscrTrial
  p_rdiscrFunc => rvector%p_rspatialDiscr
  
  if (.not. associated(p_rdiscrTest)) then
    call output_line('No discretisation attached to the matrix!',&
        OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
    call sys_halt()
  end if
  
  if (.not. associated(p_rdiscrTrial)) then
    call output_line('No discretisation attached to the matrix!',&
        OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
    call sys_halt()
  end if
  
  if (.not. associated(p_rdiscrFunc)) then
    call output_line('No discretisation attached to the matrix!',&
        OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
    call sys_halt()
  end if
  
  if ((p_rdiscrTest%inumFESpaces .ne. p_rdiscrFunc%inumFESpaces) .or. &
      (p_rdiscrTrial%inumFESpaces .ne. p_rdiscrFunc%inumFESpaces) .or. &
      (p_rdiscrTrial%inumFESpaces .ne. p_rdiscrTest%inumFESpaces)) then
    call output_line('Discretisations not compatible!',&
        OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
    call sys_halt()
  end if
  
  ! Get a pointer to the triangulation - for easier access.
  p_rtriangulation => p_rdiscrTest%p_rtriangulation
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it is at most
  ! BILF_NELEMSIM. This is only used for allocating some arrays.
  nelementsPerBlock = min(TRILF_NELEMSIM,p_rtriangulation%NEL)
  
  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.

  do icurrentElementDistr = 1,p_rdiscrTest%inumFESpaces
  
    ! Activate the current element distribution(s)
    p_elementDistrTest => &
      p_rdiscrTest%RelementDistr(icurrentElementDistr)

    p_elementDistrTrial => &
      p_rdiscrTrial%RelementDistr(icurrentElementDistr)

    p_elementDistrFunc => &
      p_rdiscrFunc%RelementDistr(icurrentElementDistr)
  
    ! Cancel if this element distribution is empty.
    if (p_elementDistrTest%NEL .eq. 0) cycle

    ! Get the number of local DOF`s for trial and test functions
    indofFunc = elem_igetNDofLoc(p_elementDistrFunc%celement)
    indofTrial = elem_igetNDofLoc(p_elementDistrTrial%celement)
    indofTest = elem_igetNDofLoc(p_elementDistrTest%celement)
    
    ! Get the number of corner vertices of the element
    NVE = elem_igetNVE(p_elementDistrTrial%celement)
    if (NVE .ne. elem_igetNVE(p_elementDistrTest%celement)) then
      call output_line('Element spaces incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
      call sys_halt()
    end if
    if (NVE .ne. elem_igetNVE(p_elementDistrFunc%celement)) then
      call output_line('Element spaces incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
      call sys_halt()
    end if
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_elementDistrTest%celement)
    
    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(p_elementDistrTest%ccubTypeBilForm)
    
    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp))
    
    ! Get the cubature formula
    call cub_getCubature(p_elementDistrTest%ccubTypeBilForm,p_DcubPtsRef, Domega)
    
    ! Open-MP-Extension: Open threads here.
    ! Each thread will allocate its own local memory...
    !
    !%omp parallel private(rintSubset, revalElementSet,&
    !%omp   p_Ddetj,kentry, dentry,DbasTest,DbasTrial, &
    !%omp   IdofsTest,IdofsTrial,cevaluationTag, &
    !%omp   Dcoefficients, bIdenticalTrialandTest, p_IdofsTrial, &
    !%omp   p_DbasTrial, BderTrial, BderTest, ielmax,IEL, idofe,jdofe,jcol0, &
    !%omp   jcol,JDFG,ICUBP, IALBET,OM,ia,ib,ifunc,aux, db,IdofsFunc,DbasFunc, &
    !%omp   bIdenticalFuncAndTrial,p_IdofsFunc, &
    !%omp   p_DbasFunc,iderType,p_IelementList,p_Ddata)    
    
    ! Quickly check if one of the specified derivatives is out of the allowed range:
    !%omp single
    do IALBET = 1,rform%itermcount
      ifunc = rform%Idescriptors(1,IALBET)
      IA = rform%Idescriptors(2,IALBET)
      IB = rform%Idescriptors(3,IALBET)      

      if ((ifunc.lt.0) .or. &
          (ifunc .gt. elem_getMaxDerivative(p_elementDistrFunc%celement))) then
         call output_line('Specified function-derivative '&
             //trim(sys_siL(ifunc,10))//' not available',&
             OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
        call sys_halt()
      end if
      
      if ((IA.le.0) .or. &
          (IA .gt. elem_getMaxDerivative(p_elementDistrTrial%celement))) then
        call output_line('Specified trial-derivative '&
             //trim(sys_siL(IA,10))//' not available',&
             OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
        call sys_halt()
      end if

      if ((IB.le.0) .or. &
          (IB .gt. elem_getMaxDerivative(p_elementDistrTest%celement))) then
        call output_line('Specified test-derivative '&
             //trim(sys_siL(IB,10))//' not available',&
             OU_CLASS_ERROR,OU_MODE_STD,'trilf_buildMatrix9d_conf2')
        call sys_halt()
      end if
    end do
    !%omp end single
    
    ! Allocate an array saving the coordinates of corner vertices of elements
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    
    allocate(DbasTest(indofTest,&
        elem_getMaxDerivative(p_elementDistrTest%celement),&
        ncubp,nelementsPerBlock))
    allocate(DbasTrial(indofTrial,&
        elem_getMaxDerivative(p_elementDistrTrial%celement), &
        ncubp,nelementsPerBlock))
    allocate(DbasFunc(indofTrial,&
        elem_getMaxDerivative(p_elementDistrFunc%celement), &
        ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(IdofsTest(indofTest,nelementsPerBlock))
    allocate(IdofsTrial(indofTrial,nelementsPerBlock))
    allocate(IdofsFunc(indofFunc,nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    allocate(Kentry(indofTrial,indofTest,nelementsPerBlock))
    allocate(Dentry(indofTrial,indofTest))
    
    ! Allocate an array for the coefficients c(x,y) f(u(x,y)). Because
    ! of the finite element function u included there, we have definitely
    ! no constant coefficients.
    allocate(Dcoefficients(rform%itermCount,ncubp,nelementsPerBlock))
                    
    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! Indicate that cubature points must still be initialised in the element set.
    bcubPtsInitialised = .false.
    
    ! p_IdofsTest points either to the just created array or to the
    ! array with the DOF`s of the trial functions - when trial and
    ! test functions are identical.
    ! We do not rely on bidenticalTrialAndTest purely, as this does not
    ! indicate whether there are identical trial and test functions
    ! in one block!
    bIdenticalTrialAndTest = p_elementDistrTrial%celement .eq. &
                             p_elementDistrTest%celement
                             
    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF`s of the test
    ! space IdofTest (if both spaces are identical). 
    ! We create a pointer for the trial space and not for the test space to
    ! prevent pointer-arithmetic in the innerst loop below!
    if (bIdenticalTrialAndTest) then
      p_IdofsTrial => IdofsTest
      p_DbasTrial  => DbasTest
      ! Build the actual combination of what the element should calculate.
      ! As we evaluate only once, what the element must calculate is an
      ! OR combination of the BDER from trial and test functions.
      BderTrial = BderTrialTempl .or. BderTestTempl
      BderTest = BderTrial
    else
      p_IdofsTrial => IdofsTrial
      p_DbasTrial  => DbasTrial
      
      ! Build the actual combination of what the element should calculate.
      ! Copy BDERxxxx to BDERxxxxAct
      BderTrial = BderTrialTempl
      BderTest = BderTestTempl
    end if
    
    ! Test whether the u trial functions are identical to the trial functions
    ! of the discretisation.
    bIdenticalFuncAndTest = p_elementDistrTest%celement .eq. &
                            p_elementDistrFunc%celement
    bIdenticalFuncAndTrial = p_elementDistrTrial%celement .eq. &
                             p_elementDistrFunc%celement
                             
    ! If yes, we can use the data calculated for the trial functions.
    if (bIdenticalFuncAndTest) then
      p_IdofsFunc => IdofsTest
      p_DbasFunc  => DbasTest
      ! Build the actual combination of what the element should calculate.
      ! As we evaluate only once, what the element must calculate is an
      ! OR combination of the BDER from trial and test functions.
      BderTest = BderTest .or. BderFuncTempl
    else if (bIdenticalFuncAndTrial) then
      p_IdofsFunc => p_IdofsTrial
      p_DbasFunc  => p_DbasTrial
      BderTrial = BderTrial .or. BderFuncTempl
    else
      p_IdofsFunc => IdofsFunc
      p_DbasFunc  => DbasFunc
      BderFunc = BderFuncTempl
    end if

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_elementDistrTest%h_IelementList, &
                              p_IelementList)
                         
    ! Get the data array from the vector
    call lsyssc_getbase_double(rvector,p_Ddata)
                              
    ! Get the number of elements there.
    NEL = p_elementDistrTest%NEL

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so TRILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !%omp do schedule(static,1)
    do IELset = 1, NEL, TRILF_NELEMSIM
    
      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(NEL,IELset-1+TRILF_NELEMSIM)
    
      ! --------------------- DOF SEARCH PHASE ------------------------
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF`s on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | IEL |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element IEL, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our BILF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(p_rdiscrTest, p_IelementList(IELset:IELmax), &
                                   IdofsTest)
                                   
      ! If the DOF`s for the test functions are different, calculate them, too.
      if (.not.bIdenticalTrialAndTest) then
        call dof_locGlobMapping_mult(p_rdiscrTrial, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
      end if

      ! If the DOF`s for the coefficient function values are different, 
      ! calculate them, too.
      if ((.not. bIdenticalFuncAndTest) .and. (.not. bIdenticalFuncAndTrial)) then
        call dof_locGlobMapping_mult(p_rdiscrFunc, p_IelementList(IELset:IELmax), &
                                     IdofsFunc)
      end if
      
      ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
      
      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF`s per element and
      ! indofTest test DOF`s per element. Therefore there are
      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
      ! "active" (i.e. have common support) on our current element, each 
      ! giving an additive contribution to the system matrix.
      !
      ! We build a quadratic indofTrial*indofTest local matrix:
      ! Kentry(1..indofTrial,1..indofTest) receives the position 
      !   in the global system matrix, where the corresponding value 
      !   has to be added to.
      ! (The corresponding contrbutions can be saved separately, 
      !  but we directly add them to the global matrix in this 
      !  approach.)
      !
      ! We build local matrices for all our elements 
      ! in the set simultaneously.
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      do IEL=1,IELmax-IELset+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"`s), as these
        ! define the rows in the matrix.
        do IDOFE=1,indofTest
        
          ! Row IDOFE of the local matrix corresponds 
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"`s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(IdofsTest(IDOFE,IEL))
          
          ! Now we loop through the other DOF`s on the current element
          ! (the "O"`s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          do JDOFE=1,indofTrial
            
            ! Get the global DOF of the "X" which interacts with 
            ! our "O".
            
            JDFG=p_IdofsTrial(JDOFE,IEL)
            
            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the DO loop if we find the column.
            
            do JCOL=JCOL0,NA
              if (p_KCOL(JCOL) .eq. JDFG) exit
            end do

            ! Because columns in the global matrix are sorted 
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! JCOL0=JCOL+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry/DENTRY this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            Kentry(JDOFE,IDOFE,IEL)=JCOL
            
          end do ! IDOFE
          
        end do ! JDOFE
        
      end do ! IEL
      
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.
      !
      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_elementDistrTrial%celement)
      cevaluationTag = ior(cevaluationTag,&
                      elem_getEvaluationTag(p_elementDistrTest%celement))
      cevaluationTag = ior(cevaluationTag,&
                      elem_getEvaluationTag(p_elementDistrFunc%celement))
                      
      if (.not. rform%ballCoeffConstant) then
        ! Evaluate real coordinates if not necessary.
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REALPOINTS)
      end if
      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! Note: Why not using
      !   IF (IELset .EQ. 1) THEN
      ! here, but this strange concept with the boolean variable?
      ! Because the IF-command does not work with OpenMP! bcubPtsInitialised
      ! is a local variable and will therefore ensure that every thread
      ! is initialising its local set of cubature points!
      if (.not. bcubPtsInitialised) then
        bcubPtsInitialised = .true.
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => revalElementSet%p_Ddetj

      ! If the matrix has nonconstant coefficients c(x,y) , calculate the 
      ! coefficients now.
      if (.not. rform%ballCoeffConstant) then
        ! Prepare an evaluation set that passes some information to the
        ! callback routine about where to evaluate
        call domint_initIntegrationByEvalSet(revalElementSet,rintSubset)
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => p_IdofsTrial
        rintSubset%celement = p_elementDistrTrial%celement
        
        call fcoeff_buildTrilMatrixSc_sim (p_rdiscrTest,p_rdiscrTrial,rform, &
                  IELmax-IELset+1,ncubp,&
                  rintSubset%p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
                  p_IdofsTrial,IdofsTest,rintSubset, &
                  Dcoefficients(:,:,1:IELmax-IELset+1),rcollection)
        call domint_doneIntegration(rintSubset)
      end if
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_elementDistrTest%celement, &
          revalElementSet, BderTest, DbasTest)
      
      ! Omit the calculation of the trial function values if they
      ! are identical to the test function values.
      if (.not. bidenticalTrialAndTest) then
        call elem_generic_sim2 (p_elementDistrTrial%celement, &
            revalElementSet, BderTrial, DbasTrial)
      end if
      
      ! Omit the calculation of the coefficient function values if they
      ! are identical to the trial function values.
      if ((.not. bIdenticalFuncAndTest) .and. (.not. bIdenticalFuncAndTrial)) then
        call elem_generic_sim2 (p_elementDistrFunc%celement, &
            revalElementSet, BderFunc, DbasFunc)
      end if
      
      ! ----------------- COEFFICIENT EVALUATION PHASE ---------------------
      
      ! Now its time to form the actual coefficients for the integral:
      ! Dcoefficients = c(x,y) f(u(x,y)). So we must evaluate f(u(x,y))
      ! and multiply it with the coefficients in Dcoefficient to get the
      ! correct coefficients for later use.
      if (rform%ballCoeffConstant) then
        ! Constant coefficients. Take the coefficients from the bilinear form
        ! and multiply with the values of f(u).
        do IALBET = 1,rform%itermcount
          iderType = rform%Idescriptors(1,IALBET)
          if (iderType .ne. 0) then
            do iel=1,IELmax-IELset+1
              do ICUBP = 1,ncubp
                ! Calculate the value in the point
                DB = 0.0_DP
                do IDOFE = 1,indofTrial
                  DB = DB + &
                    p_Ddata(p_IdofsFunc(IDOFE,iel)) * p_DbasFunc(IDOFE,iderType,ICUBP,iel)
                end do
                ! Save the value in the point, multiplied with the coefficient
                Dcoefficients(IALBET,ICUBP,iel) = rform%Dcoefficients(IALBET) * DB
              end do
            end do
          else
            Dcoefficients(IALBET,1:ncubp,1:IELmax-IELset+1) = rform%Dcoefficients(IALBET)
          end if
        end do
      else
        ! Nonconstant coefficients. Take the calculated coefficients in Dcoefficients
        ! and multiply with the values of f(u).      
        do IALBET = 1,rform%itermcount
          iderType = rform%Idescriptors(1,IALBET)
          if (iderType .ne. 0) then
            do iel=1,IELmax-IELset+1
              do ICUBP = 1,ncubp
                ! Calculate the value in the point
                DB = 0.0_DP
                do IDOFE = 1,indofTrial
                  DB = DB + &
                    p_Ddata(p_IdofsFunc(IDOFE,iel)) * p_DbasFunc(IDOFE,iderType,ICUBP,iel)
                end do
                ! Save the value in the point, multiplied with the existing coefficient
                Dcoefficients(IALBET,ICUBP,iel) = Dcoefficients(IALBET,ICUBP,iel) * DB
              end do
            end do
          end if
        end do
      end if
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate! As we have never the case of constant coefficients
      ! (because of the FE function u involved), we have only the 'complex'
      ! loop here in comparison to the standard bilinear form.
      !
      ! Loop over the elements in the current set.
      do IEL=1,IELmax-IELset+1
        
        ! Clear the local matrix
        Dentry = 0.0_DP
        ! Loop over all cubature points on the current element
        do ICUBP = 1, ncubp

          ! calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!

          OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

          ! Loop over the additive factors in the bilinear form.
          do IALBET = 1,rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be added to the matrix entry:
            !
            ! a_ij  =  int_... ( psi_j )_IA  *  ( phi_i )_IB
            !
            ! -> Ix=0: function value, 
            !      =1: first derivative, ...
            !    as defined in the module 'derivative'.
            
            IA = rform%Idescriptors(2,IALBET)
            IB = rform%Idescriptors(3,IALBET)
            
            ! Multiply OM with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            AUX = OM * Dcoefficients(IALBET,ICUBP,IEL)
            
            ! Now loop through all possible combinations of DOF`s
            ! in the current cubature point. The outer loop
            ! loops through the "O" in the above picture,
            ! the test functions:

            do IDOFE=1,indofTest
              
              ! Get the value of the (test) basis function 
              ! phi_i (our "O") in the cubature point:
              DB = DbasTest(IDOFE,IB,ICUBP,IEL)
              
              ! Perform an inner loop through the other DOF`s
              ! (the "X"). 

              do JDOFE=1,indofTrial
            
                ! Get the value of the basis function 
                ! psi_j (our "X") in the cubature point. 
                ! Them multiply:
                !    DB * DBAS(..) * AUX
                ! ~= phi_i * psi_j * coefficient * cub.weight
                ! Summing this up gives the integral, so the contribution
                ! to the global matrix. 
                !
                ! Simply summing up DB * DBAS(..) * AUX would give
                ! the coefficient of the local matrix. We save this
                ! contribution in the local matrix of element IEL.

                !JCOLB = Kentry(JDOFE,IDOFE,IEL)
                !p_DA(JCOLB) = p_DA(JCOLB) + DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
                Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE)+DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
              
              end do
            
            end do ! JDOFE
            
          end do ! IALBET

        end do ! ICUBP 
        
        ! Incorporate the local matrices into the global one.
        ! Kentry gives the position of the additive contributions in Dentry.
        !%omp critical
        do IDOFE=1,indofTest
          do JDOFE=1,indofTrial
            p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + Dentry(JDOFE,IDOFE)
          end do
        end do
        !%omp end critical

      end do ! IEL

    end do ! IELset
    !%omp end do
    
    ! Release memory
    call elprep_releaseElementSet(revalElementSet)

    deallocate(Dcoefficients)
    deallocate(IdofsTrial)
    deallocate(IdofsTest)
    deallocate(IdofsFunc)
    deallocate(DbasTrial)
    deallocate(DbasTest)
    deallocate(DbasFunc)
    deallocate(Kentry)
    deallocate(Dentry)
    
    !%omp end parallel

    deallocate(p_DcubPtsRef)
    deallocate(Domega)

  end do ! icurrentElementDistr

  ! Finish
  
  end subroutine
  
end module
