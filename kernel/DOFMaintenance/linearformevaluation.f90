!##############################################################################
!# ****************************************************************************
!# <name> linearformevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the discretisation of linear forms,
!# i.e. the creation of vectors (which usually appear as RHS vectors). 
!#
!# It contains the following set of routines:
!#
!# 1.) linf_buildVectorScalar
!#     -> Assembles the entries of a vector according to a linear form
!#        defined in terms of a volume integral.
!#
!# 2.) linf_buildVectorScalarBdr2d
!#     -> Assembles the entries of a vector according to a linear form
!#        defined in terms of a boundary integral.
!#
!# 3.) linf_buildVecIntlScalarBdr2d
!#     -> Assembles the entries of a vector according to a linear form
!#        defined in terms of a boundary integral.
!#
!# 4.) linf_buildVectorBlockBdr2d
!#     -> Assembles the entries of a vector according to a linear form
!#        defined in terms of a boundary integral.
!#
!# 5.) linf_initAssembly
!#     -> Initialise a vector assembly structure for assembling a linear form
!#
!# 6.) linf_doneAssembly
!#     -> Clean up a vector assembly structure.
!#
!# 7.) linf_assembleSubmeshVector
!#     -> Assembles the vector entries for a submesh.
!#
!# 8.) linf_assembleSubmeshVectorBdr2D
!#     -> Assembles the vector entries for a submesh.
!#
!# 9.) linf_assembleSubmeshVecIntl
!#     -> Assembles the vector entries for a submesh.
!#
!# 10.) linf_assembleSubmeshVecIntlBdr2D
!#      -> Assembles the vector entries for a submesh.
!#
!# 11.) linf_assembleSubmeshVectorBlock
!#     -> Assembles the vector entries for a submesh.
!#
!# 12.) linf_assembleSubmeshVectorBlockBdr2D
!#     -> Assembles the vector entries for a submesh.
!#
!# 13.) linf_buildVectorScalar2
!#      -> Assembles the entries of a vector according to a linear form
!#         defined in terms of a volume integral. This subroutine is a 
!#         replacement of the previous version linf_buildVectorScalar.
!#
!# 12.) linf_buildVecIntlScalar2
!#      -> Assembles the entries of a vector according to a linear form
!#         defined in terms of a volume integral.
!#
!# 13.) linf_buildVectorBlock2
!#      -> Assembles the entries of a vector according to a linear form
!#         defined in terms of a volume integral.
!#
!# It contains the following set of auxiliary routines:
!#
!# 1.) linf_buildVectorDble_conf
!#     -> Assembles the entries of a vector according to a linear form
!#        defined in terms of a volume integral for conforming discretisations.
!#
!# 2.) linf_allocAssemblyData
!#     -> Allocate 'local' memory, needed for assembling vector entries.
!#
!# 3.) linf_releaseAssemblyData
!#     -> Release 'local' memory, needed for assembling vector entries.
!#
!# Frequently asked questions \\
!# -------------------------- \\
!#
!# 1.) How to assemble a rhs-vector?
!#
!#  To assemble a vector, you first have to specify a linear form, a vector and
!#  a spatial discretisation structure that defines the FE spaces to use.
!#
!#  <code>
!#    type(t_linearForm) :: rform
!#    type(t_linfVectorAssembly) :: rvector
!#    type(t_spatialDiscretisation) :: rdiscretisation
!#  </code>
!#
!#  We set up a linear form structure, e.g. for the Laplace operator:
!#
!#  <code>
!#    rform%itermCount = 1
!#    rform%Idescriptors(1) = DER_FUNC
!#  </code>
!#
!#  This e.g. initialises a linear form for evaluating a function.
!#  In the next step, use the linear form to create the vector entries:
!#
!#  <code>
!#    call linf_buildVectorScalar (rdiscretisation,rform,.true.,rvector,coeff_RHS)
!#  </code>
!#
!#  where coeff_RHS is a callback routine returning the RHS function. That is it.
!#
!# 2.) What is the 'manual matrix assembly'?
!#
!#  This is a possibility to assemble parts of the vector by specifying
!#  an element, a cubature formula and a list of elements where to assemble.
!#  The call
!#
!#  <code>
!#    call linf_buildVectorScalar2 (rform,.true.,rvector,coeff_RHS)
!#  </code>
!#
!#  assembles a vector just like the call to linf_buildVectorScalar, but using
!#  another technique which you can also use if you want to assemble parts
!#  of the vector on your own. 
!# 
!#  To 'manually' assemble parts of the matrix, you can use the
!#  linf_initAssembly / linf_doneAssembly / linf_assembleSubmeshVector
!#  subroutines in conjunction with the linf_assembleSubmeshVector structure.
!#  Assume e.g. that elements 501..750 of a mesh are discretised with Q1
!#  and the Gauss 2x2 cubature formula in 2D. We now want to assemble a
!#  RHS on these elements. We start like before, defining
!#  a linear form, some structures and create the matrix structure:
!#
!#  <code>
!#    type(t_linearForm) :: rform
!#    type(t_vectorScalar) :: rvector
!#    type(t_vectorAssembly) :: rvectorAssembly
!#    type(t_spatialDiscretisation) :: rdiscretisation
!#  </code>
!#
!#  We set up a linear form structure, e.g. for the Laplace operator:
!#
!#  <code>
!#    rform%itermCount = 1
!#    rform%Idescriptors(1) = DER_FUNC
!#  </code>
!#
!#  Create the vector:
!#
!#  <code>
!#    call lsyssc_createVecByDiscr (rdiscretisation,rx,.true.)
!#  </code>
!#
!#  Initialise the vector assembly structure for the assembly
!#
!#  <code>
!#    call linf_initAssembly(rvectorAssembly,rform,EL_Q1,CUB_G2_2D)
!#  </code>
!#
!#  and assemble only on the elements 501..750 which we specify in Ielements:
!#
!#  <code>
!#    do i=501,750
!#      Ielements(i-500) = i
!#    end do
!#    call linf_assembleSubmeshVector (rvectorAssembly,rvector,IelementList)
!#  </code>
!#
!#  Finally, we release the assembly structure.
!#
!#    call linf_doneAssembly(rvectorAssembly)
!#
!#  So the result is a vector with the RHS assembled only in
!#  the elements 501..750. That way, the user has the ability to
!#  specify element type, element numbers and cubature formula manually
!#  without having to specify everything in a discretisation structure;
!#  that is the reason why this assembly is called 'manual' assembly.
!#
!#  This method is extremely useful when one wants to assemble vectors with
!#  adaptive/summed cubature formulas. Some parts of the domain can that way
!#  be assembled with a cubature formula which is high enough to capture the
!#  behaviour of the integral for nonsmooth, nonconstant coefficients, while
!#  other parts of the domain may be assembled with the standard cubature 
!#  formula.
!# </purpose>
!##############################################################################

module linearformevaluation

  use basicgeometry
  use boundary
  use boundaryaux
  use collection, only: t_collection
  use cubature
  use derivatives
  use dofmapping
  use domainintegration
  use element
  use elementpreprocessing
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use mprimitives
  use scalarpde
  use spatialdiscretisation
  use storage
  use transformation
  use triangulation
  
  implicit none
  
  private

!<types>
!<typeblock>

  ! A vector assembly structure that saves crucial data during the vector assembly.
  type t_linfVectorAssembly
  
    ! The bilinear form specifying the underlying PDE of the discretisation.
    type(t_linearForm) :: rform

    ! Number of local DOF`s.
    integer :: indof
    
    ! Array to tell the element which derivatives to calculate.
    logical, dimension(EL_MAXNDER) :: Bder

    ! Maximum number of elements to handle simultaneously.
    integer :: nelementsPerBlock
    
    ! Number of vertices per element
    integer :: NVE
    
    ! Type of element to evaluate.
    integer(I32) :: celement
    
    ! Type of transformation
    integer(I32) :: ctrafoType
    
    ! Basic evaluation tag of the element spaces
    integer(I32) :: cevaluationTag
    
    ! Type of cubature formula to use.
    integer(I32) :: ccubType
    
    ! Number of cubature points per element
    integer :: ncubp
    
    ! The number of elements in revalElementSet whose cubature points on the
    ! reference element(s) have already been initialised.
    integer :: iinitialisedElements
    
    ! Cubature weights
    real(DP), dimension(:), pointer :: p_Domega
    
    ! An array that takes coordinates of the cubature formula on the reference element
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    
    ! Arrays for the basis function values in the cubature points
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    
    ! Arrays saving the DOF`s in the elements
    integer, dimension(:,:), pointer :: p_Idofs

    ! Element set used for evaluating elements
    type(t_evalElementSet) :: revalElementSet
    
    ! Pointer to the coefficients that are computed by the callback routine.
    real(DP), dimension(:,:,:,:), pointer :: p_Dcoefficients
  
  end type
  
!</typeblock>

  public :: t_linfVectorAssembly

!</types>

!<constants>
!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of elements to handle simultaneously when building vectors
#ifdef LINF_NELEMSIM
  integer,public :: LINF_NELEMSIM   = LINF_NELEMSIM
#else
  integer,public :: LINF_NELEMSIM   = 1000
#endif
  
!</constantblock>
!</constants>

  public :: linf_buildVectorScalar
  public :: linf_buildVectorScalar2
  public :: linf_buildVectorScalarBdr2D
  public :: linf_buildVecIntlScalar2
  public :: linf_buildVecIntlScalarBdr2D
  public :: linf_buildVectorBlock2
  public :: linf_buildVectorBlockBdr2D
  public :: linf_initAssembly
  public :: linf_doneAssembly
  public :: linf_assembleSubmeshVector
  public :: linf_assembleSubmeshVectorBdr2D
  public :: linf_assembleSubmeshVecIntl
  public :: linf_assembleSubmeshVecIntlBdr2D
  public :: linf_assembleSubmeshVectorBlock
  public :: linf_assembleSubmeshVectorBlockBdr2D
  public :: linf_allocAssemblyData
  public :: linf_releaseAssemblyData

contains

  !****************************************************************************

!<subroutine>

  subroutine linf_buildVectorScalar (rdiscretisation, rform, bclear, rvectorScalar,&
                                     fcoeff_buildVectorSc_sim, rcollection)
  
!<description>
  ! This routine assembles the entries of a vector according to a linear form
  ! (typically used for assembling RHS vectors).
  !
  ! If bclear=TRUE, the vector is cleared before the assembly and any 
  ! sorting of the entries is switched off - the vector is set up unsorted.
  !
  ! If bclear=FALSE, the vector must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
!</description>

!<input>
  ! The underlying discretisation structure which is to be used to
  ! create the vector.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  
  ! The linear form specifying the underlying PDE of the discretisation.
  type(t_linearForm), intent(in) :: rform
  
  ! Whether to clear the vector before calculating the entries.
  ! If .FALSE., the new entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! A callback routine for the function to be discretised.
  include 'intf_coefficientVectorSc.inc'
  optional :: fcoeff_buildVectorSc_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated entries are imposed to this vector.
  type(t_vectorScalar), intent(inout) :: rvectorScalar

  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection
!</inputoutput>

!</subroutine>
  
  ! If the vector is not set up as new vector, it has to be unsorted.
  ! If it is a new vector, we switch off the sorting.
  if (bclear) then
    rvectorScalar%isortStrategy = -abs(rvectorScalar%isortStrategy)
  end if
  
  ! The vector must be unsorted, otherwise we can not set up the vector.
  if (rvectorScalar%isortStrategy .gt. 0) then
    call output_line('Vector must be unsorted!',&
        OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalar')
    call sys_halt()
  end if

  ! Do we have a uniform triangulation? Would simplify a lot...
  if ((rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) .or.&
      (rdiscretisation%ccomplexity .eq. SPDISC_CONFORMAL)) then 
    
    select case(rvectorScalar%cdataType)
      
    case(ST_DOUBLE)
      call linf_buildVectorDble_conf (rdiscretisation, rform, bclear, rVectorScalar,&  
                                      fcoeff_buildVectorSc_sim, rcollection)

    case DEFAULT
      call output_line('Single precision vectors currently not supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalar')
      call sys_halt()
    end select
    
  else
    call output_line('General discretisation not implemented!',&
        OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalar')
    call sys_halt()
  end if

  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine linf_buildVectorDble_conf (rdiscretisation, rform, bclear, rVectorScalar,&
                                        fcoeff_buildVectorSc_sim, rcollection)

!<description>
  ! This routine calculates the entries of a discretised finite element vector.
  ! The discretisation is assumed to be conformal, i.e. the DOF`s
  ! of all finite elements must 'match'. 
  ! The linear form is defined by
  !        <tex>(f,$\phi_i$), i=1..*</tex>
  ! with <tex>$\phi_i$</tex> being the test functions defined in the discretisation
  ! structure.
  ! In case the array for the vector entries does not exist, the routine
  ! allocates memory in size of the vector of the heap for the vector entries
  ! and initialises all necessary variables of the vector according to the
  ! parameters (NEQ, pointer to the discretisation,...)
  !
  ! Double-precision version.
!</description>

!<input>
  ! The underlying discretisation structure which is to be used to
  ! create the vector.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  
  ! The linear form specifying the underlying PDE of the discretisation.
  type(t_linearForm), intent(in) :: rform
  
  ! Whether to clear the vector before calculating the entries.
  ! If .FALSE., the new vector entries are added to the existing entries.
  logical, intent(in) :: bclear

  
  ! A callback routine which is able to calculate the values of the
  ! function <tex>$f$</tex> which is to be discretised.
  include 'intf_coefficientVectorSc.inc'
  optional :: fcoeff_buildVectorSc_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated vector entries are added to this vector.
  type(t_vectorScalar), intent(inout) :: rvectorScalar
 
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,i1,icurrentElementDistr, ICUBP, IALBET, IA
  integer :: IEL, IELmax, IELset, IDOFE
  real(DP) :: OM,AUX
  
  ! Array to tell the element which derivatives to calculate
  logical, dimension(EL_MAXNDER) :: Bder
  
  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(:), allocatable :: Domega
  
  ! number of cubature points on the reference element
  integer :: ncubp
  
  ! Pointer to the vector entries
  real(DP), dimension(:), pointer :: p_Ddata

  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsTest
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test space.
  real(DP), dimension(:,:,:,:), allocatable, target :: DbasTest
  
  ! Number of entries in the vector - for quicker access
  integer :: NEQ
  
  ! Type of transformation from the reference to the real element 
  integer(I32) :: ctrafoType
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag

  ! Number of local degees of freedom for test functions
  integer :: indofTest
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList
  
  ! A small vector holding only the additive contributions of
  ! one element
  real(DP), dimension(EL_MAXNBAS) :: DlocalData
  
  ! An array that takes coordinates of the cubature formula on the reference element
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

  ! Pointer to the jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj
  
  ! Current element distribution
  type(t_elementDistribution), pointer :: p_elementDistribution
  
  ! Number of elements in the current element distribution
  integer :: NEL

  ! Number of elements in a block. Normally =LINF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock
  
  ! Pointer to the coefficients that are computed by the callback routine.
  real(DP), dimension(:,:,:), allocatable :: Dcoefficients
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  type(t_domainIntSubset) :: rintSubset
  type(t_evalElementSet) :: revalElementSet
  logical :: bcubPtsInitialised
  
  
  ! Which derivatives of basis functions are needed?
  ! Check the descriptors of the bilinear form and set BDER
  ! according to these.
  
  Bder = .false.
  
  ! Loop through the additive terms
  do i = 1, rform%itermCount
    ! The desriptor Idescriptors gives directly the derivative
    ! which is to be computed!
    I1 = rform%Idescriptors(i)
    
    if ((I1 .le.0) .or. (I1 .gt. DER_MAXNDER)) then
      call output_line('Invalid descriptor',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorDble_conf')
      call sys_halt()
    endif
    
    Bder(I1)=.true.
  end do
  
  if (rvectorScalar%h_Ddata .eq. ST_NOHANDLE) then
  
    ! Get the size of the vector and put it to the vector structure.
    NEQ = dof_igetNDofGlob(rdiscretisation)
    
    ! Initialise the vector parameters
    rvectorScalar%NEQ             = NEQ
    rvectorScalar%iidxFirstEntry  = 1
    rvectorScalar%p_rspatialDiscr => rdiscretisation
    rvectorScalar%cdataType       = ST_DOUBLE

    ! Clear the entries in the vector - we need to start with zero
    ! when assembling a new vector.
    call storage_new ('linf_buildVectorDble_conf', 'vector', &
                        NEQ, ST_DOUBLE, rvectorScalar%h_Ddata, &
                        ST_NEWBLOCK_ZERO)
    call lsyssc_getbase_double (rvectorScalar,p_Ddata)

  else
  
    ! Get information about the vector:
    NEQ = rvectorScalar%NEQ
  
    ! Get the data array.
    call lsyssc_getbase_double (rvectorScalar,p_Ddata)
    
    ! If desired, clear the vector before assembling.
    if (bclear) then
      call lalg_clearVectorDble (p_Ddata)
    end if
    
  end if
  
  ! Get a pointer to the triangulation - for easier access.
  p_rtriangulation => rdiscretisation%p_rtriangulation
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it is at most
  ! LINF_NELEMSIM. This is only used for allocating some arrays.
  nelementsPerBlock = min(LINF_NELEMSIM, p_rtriangulation%NEL)
  
  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.
  do icurrentElementDistr = 1,rdiscretisation%inumFESpaces
  
    ! Activate the current element distribution
    p_elementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
  
    ! Cancel if this element distribution is empty.
    if (p_elementDistribution%NEL .eq. 0) cycle

    ! Get the number of local DOF`s for trial and test functions
    indofTest = elem_igetNDofLoc(p_elementDistribution%celement)
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_elementDistribution%celement)

    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(p_elementDistribution%ccubTypeLinForm)

    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp))
    
    ! Get the cubature formula
    call cub_getCubature(p_elementDistribution%ccubTypeLinForm,p_DcubPtsRef, Domega)
    
    ! Open-MP-Extension: Open threads here.
    ! Each thread will allocate its own local memory...
    !
    !%OMP PARALLEL PRIVATE(rintSubset, revalElementSet,&
    !%OMP   p_Ddetj,DbasTest, cevaluationTag, bcubPtsInitialised,&
    !%OMP   IdofsTest,&
    !%OMP   Dcoefficients, &
    !%OMP   IELmax,IEL, idofe, &
    !%OMP   ICUBP, IALBET,OM,IA,aux)    
    
    ! Quickly check if one of the specified derivatives is out of the allowed range:
    do IALBET = 1,rform%itermcount
      IA = rform%Idescriptors(IALBET)
      if ((IA.lt.0) .or. &
          (IA .gt. elem_getMaxDerivative(p_elementDistribution%celement))) then
        call output_line('Specified test-derivative '//trim(sys_siL(IA,3))//' not available!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorDble_conf')
        call sys_halt()
      end if
    end do

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    allocate(DbasTest(indofTest,elem_getMaxDerivative(p_elementDistribution%celement),&
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(IdofsTest(indofTest,nelementsPerBlock))

    ! Allocate memory for the coefficients
    allocate(Dcoefficients(rform%itermCount,ncubp,nelementsPerBlock))
  
    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! Indicate that cubature points must still be initialised in the element set.
    bcubPtsInitialised = .false.
    
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_elementDistribution%h_IelementList, &
                              p_IelementList)
                              
    ! Get the number of elements there.
    NEL = p_elementDistribution%NEL
  
  
    ! Loop over the elements - blockwise.
    !%OMP do schedule(static,1)
    do IELset = 1, NEL, LINF_NELEMSIM
    
      ! We always handle LINF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(NEL,IELset-1+LINF_NELEMSIM)
    
      ! Calculate the global DOF`s into IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                   IdofsTest)
      
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! Ok, we found the positions of the local vector entries
      ! that we have to change.
      ! To calculate the vector contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_elementDistribution%celement)
                      
      ! Evaluate real coordinates; they are needed in the callback function.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REALPOINTS)
      
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

      ! Now it is time to call our coefficient function to calculate the
      ! function values in the cubature points:
      call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
      rintSubset%ielementDistribution = icurrentElementDistr
      rintSubset%ielementStartIdx = IELset
      rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
      rintSubset%p_IdofsTrial => IdofsTest
      rintSubset%celement = p_elementDistribution%celement
      
      if (present(fcoeff_buildVectorSc_sim)) then
        call fcoeff_buildVectorSc_sim (rdiscretisation,rform, &
            IELmax-IELset+1,ncubp,revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1), &
            IdofsTest,rintSubset, &
            Dcoefficients(:,:,1:IELmax-IELset+1),rcollection)
      else
        Dcoefficients(:,:,1:IELmax-IELset+1) = 1.0_DP
      end if
        
      ! Release the domain integration subset again
      call domint_doneIntegration(rintSubset)
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_elementDistribution%celement, &
          revalElementSet, Bder, DbasTest)
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the DOF`s and cubature points to calculate the
      ! integral:
      
      do IEL = 1,IELmax-IELset+1
        
        ! We make a 'local' approach, i.e. we calculate the values of the
        ! integral into the vector DlocalData and add them later into
        ! the large solution vector.
        
        ! Clear the output vector.
        DlocalData(1:indofTest) = 0.0_DP

        ! Loop over all cubature points on the current element
        do ICUBP = 1, ncubp
        
          ! calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!

          OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

          ! Loop over the additive factors in the linear form.
          do IALBET = 1,rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be:
            !
            ! int_...  f * ( phi_i )_IA
            !
            ! -> IA=0: function value, 
            !      =1: first derivative, 
            !      =2: 2nd derivative,...
            !    as defined in the module 'derivative'.
            
            IA = rform%Idescriptors(IALBET)
            
            ! Multiply OM with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            AUX = OM * Dcoefficients(IALBET,ICUBP,IEL)
          
            ! Now loop through all possible combinations of DOF`s
            ! in the current cubature point. 

            do IDOFE = 1,indofTest
            
              ! Get the value of the basis function 
              ! phi_o in the cubature point. 
              ! Them multiply:
              !    DBAS(..) * AUX
              ! ~= phi_i * coefficient * cub.weight
              ! Summing this up gives the integral, so the contribution
              ! to the vector. 
              !
              ! Simply summing up DBAS(..) * AUX would give
              ! the additive contribution for the vector. We save this
              ! contribution in the local array.

              DlocalData(IDOFE) = DlocalData(IDOFE)+DbasTest(IDOFE,IA,ICUBP,IEL)*AUX
            
            end do ! IDOFE
            
          end do ! IALBET

        end do ! ICUBP 

        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        do IDOFE = 1,indofTest
          p_Ddata(IdofsTest(IDOFE,IEL)) = p_Ddata(IdofsTest(IDOFE,IEL)) + DlocalData(IDOFE)
        end do

      end do ! IEL
      
    end do ! IELset
    !%OMP END DO
    
    ! Release memory
    deallocate(Dcoefficients)
    deallocate(IdofsTest)
    deallocate(DbasTest)

    call elprep_releaseElementSet(revalElementSet)
    
    !%OMP END PARALLEL

    deallocate(p_DcubPtsRef)
    deallocate(Domega)

  end do ! icurrentElementDistr
  
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine linf_buildVectorScalarBdr2d (rform, ccubType, bclear, rvectorScalar,&
                                          fcoeff_buildVectorScBdr2D_sim,&
                                          rboundaryRegion, rcollection)
  
!<description>
  ! This routine assembles the entries of a vector according to a linear form
  ! (typically used for assembling RHS vectors).
  !
  ! If bclear=TRUE, the vector is cleared before the assembly and any 
  ! sorting of the entries is switched off - the vector is set up unsorted.
  !
  ! If bclear=FALSE, the vector must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
!</description>

!<input>
  ! The linear form specifying the underlying PDE of the discretisation.
  type(t_linearForm), intent(in) :: rform

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer(I32), intent(in) :: ccubType
  
  ! Whether to clear the vector before calculating the entries.
  ! If .FALSE., the new entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion
  
  ! A callback routine for the function to be discretised.
  include 'intf_coefficientVectorScBdr2D.inc'
  optional :: fcoeff_buildVectorScBdr2D_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated entries are imposed to this vector.
  type(t_vectorScalar), intent(inout) :: rvectorScalar
  
  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_linfVectorAssembly) :: rvectorAssembly
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundaryRegion) :: rboundaryReg
    real(DP), dimension(:,:), pointer :: DedgePosition
    integer, dimension(:), pointer :: IelementList, IelementOrientation
    integer :: ibdc,ielementDistr,NELbdc
    
    ! If the vector does not exist, stop here.
    if (rvectorScalar%h_Ddata .eq. ST_NOHANDLE) then  
      call output_line('Vector not available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalarBdr2D')
    end if
    
    ! If the vector os stored in interleave format, stop here.
    if (rvectorScalar%NVAR .ne. 1) then
      call output_line('Vector must not be in interleaved format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalarBdr2D')
    end if
    
    ! The vector must be unsorted, otherwise we can not set up the vector.
    if (rvectorScalar%isortStrategy .gt. 0) then
      call output_line('Vector must be unsorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalarBdr2D')
      call sys_halt()
    end if
    
    ! Clear the vector if necessary.
    if (bclear) call lsyssc_clearVector (rvectorScalar)
    
    ! The vector must provide a discretisation structure
    if (.not. associated(rvectorScalar%p_rspatialDiscr)) then
      call output_line('No discretisation associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalarBdr2D')
      call sys_halt()
    end if
    
    ! The discretisation must provide a triangulation structure
    if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rtriangulation)) then
      call output_line('No triangulation associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalarBdr2d')
      call sys_halt()
    end if
    
    ! The discretisation must provide a boundary structure
    if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rboundary)) then
      call output_line('No boundary associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalarBdr2d')
      call sys_halt()
    end if
    
    ! Set pointers for quicker access
    p_rboundary => rvectorScalar%p_rspatialDiscr%p_rboundary
    p_rtriangulation => rvectorScalar%p_rspatialDiscr%p_rtriangulation
    
    ! Do we have a uniform triangulation? Would simplify a lot...
    if ((rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) .or.&
        (rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_CONFORMAL)) then 
      
      select case(rvectorScalar%cdataType)
        
      case(ST_DOUBLE)
        
        if (present(rboundaryRegion)) then
          
          ! Calculate total number of elements adjacent to the boundary region
          NELbdc = bdraux_getNELAtRegion(rboundaryRegion, p_rtriangulation)
          
          ! Allocate memory for element list, element orientation and
          ! the start- and end-parameter values of edges at the boundary
          allocate(IelementList(NELbdc), IelementOrientation(NELbdc))
          allocate(DedgePosition(2,NELbdc))
          
          ! Loop over the element distributions.
          do ielementDistr = 1,rvectorScalar%p_rspatialDiscr%inumFESpaces
            
            ! Calculate the list of elements adjacent to the boundary
            call bdraux_getElementsAtRegion(rboundaryRegion,&
                rvectorScalar%p_rspatialDiscr, NELbdc,&
                IelementList, IelementOrientation, DedgePosition,&
                rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%celement)
            
            ! Check if element distribution is empty
            if (NELbdc .le. 0) cycle
            
            ! Initialise a vector assembly structure for that element distribution
            call linf_initAssembly(rvectorAssembly, rform,&
                rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%celement,&
                ccubType, min(LINF_NELEMSIM, NELbdc))
            
            ! Assemble the data for all elements in this element distribution
            call linf_assembleSubmeshVectorBdr2D (rvectorAssembly, rvectorScalar,&
                rboundaryRegion, IelementList(1:NELbdc), IelementOrientation(1:NELbdc),&
                DedgePosition(:,1:NELbdc), fcoeff_buildVectorScBdr2D_sim, rcollection)
            
            ! Release the assembly structure.
            call linf_doneAssembly(rvectorAssembly)
            
          end do
          
          ! Release memory
          deallocate(IelementList, IelementOrientation, DedgePosition)
          
        else
          
          ! Loop over the element distributions.
          do ielementDistr = 1,rvectorScalar%p_rspatialDiscr%inumFESpaces
            
            ! Check if element distribution is empty
            if (rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%NEL .le. 0) cycle
            
            ! Initialise a vector assembly structure for that element distribution
            call linf_initAssembly(rvectorAssembly, rform,&
                rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%celement,&
                ccubType, LINF_NELEMSIM)
            
            ! Create a boundary region for each boundary component and call
            ! the calculation routine for that.
            do ibdc = 1,boundary_igetNBoundComp(p_rboundary)
              call boundary_createRegion (p_rboundary, ibdc, 0, rboundaryReg)
              
              ! Calculate number of elements adjacent to the boundary region
              NELbdc = bdraux_getNELAtRegion(rboundaryReg, p_rtriangulation)
              
              ! Check if element distribution is empty
              if (NELbdc .le. 0) cycle
              
              ! Allocate memory for element list, element orientation and
              ! the start- and end-parameter values of edges at the boundary
              allocate(IelementList(NELbdc), IelementOrientation(NELbdc))
              allocate(DedgePosition(2,NELbdc))
              
              ! Calculate the list of elements adjacent to the boundary
              call bdraux_getElementsAtRegion(rboundaryReg,&
                  rvectorScalar%p_rspatialDiscr, NELbdc,&
                  IelementList, IelementOrientation, DedgePosition,&
                  rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%celement)
              
              if (NELbdc .gt. 0) then
                
                ! Assemble the data for all elements in this element distribution
                call linf_assembleSubmeshVectorBdr2D (rvectorAssembly, rvectorScalar,&
                    rboundaryReg, IelementList(1:NELbdc), IelementOrientation(1:NELbdc),&
                    DedgePosition(:,1:NELbdc), fcoeff_buildVectorScBdr2D_sim, rcollection)
                
              end if
              
              ! Deallocate memory
              deallocate(IelementList, IelementOrientation, DedgePosition)
              
            end do ! ibdc
            
            ! Release the assembly structure.
            call linf_doneAssembly(rvectorAssembly)
            
          end do ! ielementDistr
          
        end if
        
      case DEFAULT
        call output_line('Single precision vectors currently not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalarBdr2D')
      end select
      
    else
      call output_line('General discretisation not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalarBdr2D')
      call sys_halt()
    end if
 
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine linf_buildVecIntlScalarBdr2d (rform, ccubType, bclear, rvectorScalar,&
                                           fcoeff_buildVectorBlBdr2D_sim,&
                                           rboundaryRegion, rcollection)
  
!<description>
  ! This routine assembles the entries of a vector according to a linear form
  ! (typically used for assembling RHS vectors).
  !
  ! If bclear=TRUE, the vector is cleared before the assembly and any 
  ! sorting of the entries is switched off - the vector is set up unsorted.
  !
  ! If bclear=FALSE, the vector must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
!</description>

!<input>
  ! The linear form specifying the underlying PDE of the discretisation.
  type(t_linearForm), intent(in) :: rform

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer(I32), intent(in) :: ccubType
  
  ! Whether to clear the vector before calculating the entries.
  ! If .FALSE., the new entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion
  
  ! A callback routine for the function to be discretised.
  include 'intf_coefficientVectorBlBdr2D.inc'
  optional :: fcoeff_buildVectorBlBdr2D_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated entries are imposed to this vector.
  type(t_vectorScalar), intent(inout) :: rvectorScalar
  
  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_linfVectorAssembly) :: rvectorAssembly
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundaryRegion) :: rboundaryReg
    real(DP), dimension(:,:), pointer :: DedgePosition
    integer, dimension(:), pointer :: IelementList, IelementOrientation
    integer :: ibdc,ielementDistr,NELbdc
    
    ! If the vector does not exist, stop here.
    if (rvectorScalar%h_Ddata .eq. ST_NOHANDLE) then  
      call output_line('Vector not available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalarBdr2D')
    end if
    
    ! The vector must be unsorted, otherwise we can not set up the vector.
    if (rvectorScalar%isortStrategy .gt. 0) then
      call output_line('Vector must be unsorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalarBdr2D')
      call sys_halt()
    end if
    
    ! Clear the vector if necessary.
    if (bclear) call lsyssc_clearVector (rvectorScalar)
    
    ! The vector must provide a discretisation structure
    if (.not. associated(rvectorScalar%p_rspatialDiscr)) then
      call output_line('No discretisation associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalarBdr2D')
      call sys_halt()
    end if
    
    ! The discretisation must provide a triangulation structure
    if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rtriangulation)) then
      call output_line('No triangulation associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalarBdr2d')
      call sys_halt()
    end if
    
    ! The discretisation must provide a boundary structure
    if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rboundary)) then
      call output_line('No boundary associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalarBdr2d')
      call sys_halt()
    end if
    
    ! Set pointers for quicker access
    p_rboundary => rvectorScalar%p_rspatialDiscr%p_rboundary
    p_rtriangulation => rvectorScalar%p_rspatialDiscr%p_rtriangulation
    
    ! Do we have a uniform triangulation? Would simplify a lot...
    if ((rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) .or.&
        (rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_CONFORMAL)) then 
      
      select case(rvectorScalar%cdataType)
        
      case(ST_DOUBLE)
        
        if (present(rboundaryRegion)) then
          
          ! Calculate total number of elements adjacent to the boundary region
          NELbdc = bdraux_getNELAtRegion(rboundaryRegion, p_rtriangulation)
          
          ! Allocate memory for element list, element orientation and
          ! the start- and end-parameter values of edges at the boundary
          allocate(IelementList(NELbdc), IelementOrientation(NELbdc))
          allocate(DedgePosition(2,NELbdc))
          
          ! Loop over the element distributions.
          do ielementDistr = 1,rvectorScalar%p_rspatialDiscr%inumFESpaces
            
            ! Calculate the list of elements adjacent to the boundary
            call bdraux_getElementsAtRegion(rboundaryRegion,&
                rvectorScalar%p_rspatialDiscr, NELbdc,&
                IelementList, IelementOrientation, DedgePosition,&
                rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%celement)
            
            ! Check if element distribution is empty
            if (NELbdc .le. 0) cycle
            
            ! Initialise a vector assembly structure for that element distribution
            call linf_initAssembly(rvectorAssembly, rform,&
                rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%celement,&
                ccubType, min(LINF_NELEMSIM, NELbdc))
            
            ! Assemble the data for all elements in this element distribution
            call linf_assembleSubmeshVecIntlBdr2D (rvectorAssembly, rvectorScalar,&
                rboundaryRegion, IelementList(1:NELbdc), IelementOrientation(1:NELbdc),&
                DedgePosition(:,1:NELbdc), fcoeff_buildVectorBlBdr2D_sim, rcollection)
            
            ! Release the assembly structure.
            call linf_doneAssembly(rvectorAssembly)
            
          end do
          
          ! Release memory
          deallocate(IelementList, IelementOrientation, DedgePosition)
          
        else
          
          ! Loop over the element distributions.
          do ielementDistr = 1,rvectorScalar%p_rspatialDiscr%inumFESpaces
            
            ! Check if element distribution is empty
            if (rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%NEL .le. 0) cycle
            
            ! Initialise a vector assembly structure for that element distribution
            call linf_initAssembly(rvectorAssembly, rform,&
                rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%celement,&
                ccubType, LINF_NELEMSIM)
            
            ! Create a boundary region for each boundary component and call
            ! the calculation routine for that.
            do ibdc = 1,boundary_igetNBoundComp(p_rboundary)
              call boundary_createRegion (p_rboundary, ibdc, 0, rboundaryReg)
              
              ! Calculate number of elements adjacent to the boundary region
              NELbdc = bdraux_getNELAtRegion(rboundaryReg, p_rtriangulation)
              
              ! Check if element distribution is empty
              if (NELbdc .le. 0) cycle
              
              ! Allocate memory for element list, element orientation and
              ! the start- and end-parameter values of edges at the boundary
              allocate(IelementList(NELbdc), IelementOrientation(NELbdc))
              allocate(DedgePosition(2,NELbdc))
              
              ! Calculate the list of elements adjacent to the boundary
              call bdraux_getElementsAtRegion(rboundaryReg,&
                  rvectorScalar%p_rspatialDiscr, NELbdc,&
                  IelementList, IelementOrientation, DedgePosition,&
                  rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%celement)
              
              if (NELbdc .gt. 0) then
                
                ! Assemble the data for all elements in this element distribution
                call linf_assembleSubmeshVecIntlBdr2D (rvectorAssembly, rvectorScalar,&
                    rboundaryReg, IelementList(1:NELbdc), IelementOrientation(1:NELbdc),&
                    DedgePosition(:,1:NELbdc), fcoeff_buildVectorBlBdr2D_sim, rcollection)
                
              end if
              
              ! Deallocate memory
              deallocate(IelementList, IelementOrientation, DedgePosition)
              
            end do ! ibdc
            
            ! Release the assembly structure.
            call linf_doneAssembly(rvectorAssembly)
            
          end do ! ielementDistr
          
        end if
        
      case DEFAULT
        call output_line('Single precision vectors currently not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalarBdr2D')
      end select
      
    else
      call output_line('General discretisation not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalarBdr2D')
      call sys_halt()
    end if
 
  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine linf_initAssembly(rvectorAssembly,rform,&
      celement,ccubType,nelementsPerBlock)

!<description>
  ! Initialise a vector assembly structure for assembling a linear form.
!</description>

!<input>
  ! The bilinear form specifying the underlying PDE of the discretisation.
  type(t_linearForm), intent(in) :: rform
  
  ! Type of element in the test space.
  integer(I32), intent(in) :: celement
  
  ! Type of cubature formula to use.
  integer(I32), intent(in) :: ccubType
  
  ! Optional: Maximum number of elements to process simultaneously.
  ! If not specified, LINF_NELEMSIM is assumed.
  integer, intent(in), optional :: nelementsPerBlock
!</input>

!<output>
  ! A vector assembly structure.
  type(t_linfVectorAssembly), intent(out) :: rvectorAssembly
!</output>

!</subroutine>
  
    ! local variables
    integer :: i,i1
  
    ! Initialise the structure.
    rvectorAssembly%rform = rform
    rvectorAssembly%ccubType = ccubType
    rvectorAssembly%nelementsPerBlock = LINF_NELEMSIM
    if (present(nelementsPerBlock)) &
        rvectorAssembly%nelementsPerBlock = nelementsPerBlock
    rvectorAssembly%celement = celement
    
    ! Get the number of local DOF`s for trial and test functions
    rvectorAssembly%indof = elem_igetNDofLoc(celement)
    
    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDERxxxx
    ! according to these.
    rvectorAssembly%Bder(:) = .false.
    
    ! Loop through the additive terms
    do i = 1,rform%itermCount
      ! The desriptor Idescriptors gives directly the derivative
      ! which is to be computed! Build templates for BDER.
      ! We do not compute the actual BDER here, as there might be some special
      ! processing if trial/test functions are identical!
      !
      ! At first build the descriptors for the trial functions
      I1=rform%Idescriptors(i)
      
      if ((I1 .le.0) .or. (I1 .gt. DER_MAXNDER)) then
        call output_line ('Invalid descriptor!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_initAssembly')
        call sys_halt()
      endif
      
      rvectorAssembly%Bder(I1)=.true.
    end do

    ! Get the number of vertices of the element, specifying the transformation
    ! form the reference to the real element.
    rvectorAssembly%NVE = elem_igetNVE(celement)
    
      ! Get from the element space the type of coordinate system
    ! that is used there:
    rvectorAssembly%ctrafoType = elem_igetTrafoType(celement)
    
    ! Get the number of cubature points for the cubature formula
    rvectorAssembly%ncubp = cub_igetNumPts(ccubType)

    ! Allocate two arrays for the points and the weights
    allocate(rvectorAssembly%p_Domega(rvectorAssembly%ncubp))
    allocate(rvectorAssembly%p_DcubPtsRef(&
        trafo_igetReferenceDimension(rvectorAssembly%ctrafoType),&
        rvectorAssembly%ncubp))
    
    ! Get the cubature formula
    call cub_getCubature(ccubType,rvectorAssembly%p_DcubPtsRef,rvectorAssembly%p_Domega)

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    rvectorAssembly%cevaluationTag = elem_getEvaluationTag(rvectorAssembly%celement)
        
  end subroutine
  
  !****************************************************************************

!<subroutine>
  
  subroutine linf_doneAssembly(rvectorAssembly)

!<description>
  ! Clean up a vector assembly structure.
!</description>  

!<inputoutput>
  ! Matrix assembly structure to clean up
  type(t_linfVectorAssembly), intent(inout) :: rvectorAssembly
!</inputoutput>  

!</subroutine>
  
    deallocate(rvectorAssembly%p_DcubPtsRef)
    deallocate(rvectorAssembly%p_Domega)
  
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine linf_allocAssemblyData(rvectorAssembly, nblocks)

!<description>
  ! Auxiliary subroutine.
  ! Allocate 'local' memory, needed for assembling vector entries.
!</description>

!<input>
  ! Optional: Number of blocks present in block or interleaved vectors
  !           If not present, nblocks=1 is assumed.
  integer, intent(in), optional :: nblocks
!</input>

!<inputoutput>
  ! A vector assembly structure.
  type(t_linfVectorAssembly), intent(inout) :: rvectorAssembly
!</inputoutput>

!</subroutine>

    ! Evaluate real coordinates.
    rvectorAssembly%cevaluationTag = &
        ior(rvectorAssembly%cevaluationTag,EL_EVLTAG_REALPOINTS)
    
    ! Allocate an array for the coefficients computed by the callback routine.
    if (present(nblocks)) then
      allocate(rvectorAssembly%p_Dcoefficients(max(1,nblocks),&
               rvectorAssembly%rform%itermCount,rvectorAssembly%ncubp,&
               rvectorAssembly%nelementsPerBlock))
    else
      allocate(rvectorAssembly%p_Dcoefficients(1,&
               rvectorAssembly%rform%itermCount,rvectorAssembly%ncubp,&
               rvectorAssembly%nelementsPerBlock))
    end if

    ! Allocate arrays for the values of the test functions.
    allocate(rvectorAssembly%p_Dbas(rvectorAssembly%indof,&
             elem_getMaxDerivative(rvectorAssembly%celement),&
             rvectorAssembly%ncubp,rvectorAssembly%nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(rvectorAssembly%p_Idofs(&
             rvectorAssembly%indof,rvectorAssembly%nelementsPerBlock))

    ! Initialisation of the element set.
    call elprep_init(rvectorAssembly%revalElementSet)
  
    ! No cubature points initalised up to now.
    rvectorAssembly%iinitialisedElements = 0
  
  end subroutine
  
  !****************************************************************************
  
!<subroutine>
  
  subroutine linf_releaseAssemblyData(rvectorAssembly)

!<description>
  ! Auxiliary subroutine.
  ! Release 'local' memory.
!</description>

!<inputoutput>
  ! Vector assembly structure to clean up
  type(t_linfVectorAssembly), intent(inout) :: rvectorAssembly
!</inputoutput>

!</subroutine>
  
    ! Release all information in the structure.
    call elprep_releaseElementSet(rvectorAssembly%revalElementSet)

    deallocate(rvectorAssembly%p_Dcoefficients)
    deallocate(rvectorAssembly%p_Idofs)
    deallocate(rvectorAssembly%p_Dbas)   
       
  end subroutine
  
  !****************************************************************************
  
!<subroutine>  
  
  subroutine linf_assembleSubmeshVector (rvectorAssembly, rvector, IelementList,&
      fcoeff_buildVectorSc_sim, rcollection)

!<description>

  ! Assembles the vector entries for a submesh by integrating over the domain.

!</description>
  
!<input>
  
  ! List of elements where to assemble the bilinear form.
  integer, dimension(:), intent(in), target :: IelementList
  
  ! A callback routine which is able to calculate the values of the
  ! function $f$ which is to be discretised.
  include 'intf_coefficientVectorSc.inc'
  optional :: fcoeff_buildVectorSc_sim
  
!</input>

!<inputoutput>
  
  ! A vector assembly structure prepared with linf_initAssembly.
  type(t_linfVectorAssembly), intent(inout), target :: rvectorAssembly
  
  ! A vector where to assemble the contributions to.
  type(t_vectorScalar), intent(inout) :: rvector
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection

!</inputoutput>

!</subroutine>
  
    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: indof,ncubp
    
    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax
    integer :: iel,icubp,ialbet,ia,idofe
    real(DP) :: domega,daux
    integer(I32) :: cevaluationTag
    type(t_linfVectorAssembly), target :: rlocalVectorAssembly
    type(t_domainIntSubset) :: rintSubset
    real(DP), dimension(:,:), pointer :: p_Ddetj
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
    integer, dimension(:),pointer :: p_Idescriptors
    integer, dimension(:,:), pointer :: p_Idofs
    type(t_evalElementSet), pointer :: p_revalElementSet

    ! A small vector holding only the additive contributions of
    ! one element
    real(DP), dimension(EL_MAXNBAS) :: DlocalData
  
    ! Get some pointers for faster access
    call lsyssc_getbase_double (rvector,p_Ddata)
    indof = rvectorAssembly%indof
    ncubp = rvectorAssembly%ncubp

    ! Copy the assembly data to the local assembly data,
    ! where we can allocate memory.
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the vector
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original assembly
    ! stucture or disturbing the data of the other processors.
    rlocalVectorAssembly = rvectorAssembly
    call linf_allocAssemblyData(rlocalVectorAssembly)
    
    ! Get some more pointers to local data.
    p_Domega => rlocalVectorAssembly%p_Domega
    p_Dbas => rlocalVectorAssembly%p_Dbas
    p_Dcoefficients => rlocalVectorAssembly%p_Dcoefficients(1,:,:,:)
    p_Idescriptors => rlocalVectorAssembly%rform%Idescriptors
    p_Idofs => rlocalVectorAssembly%p_Idofs
    p_revalElementSet => rlocalVectorAssembly%revalElementSet
        
    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !%OMP do schedule(static,1)
    do IELset = 1, size(IelementList), rlocalVectorAssembly%nelementsPerBlock
    
      ! We always handle nelementsPerBlock elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(IelementList),IELset-1+rlocalVectorAssembly%nelementsPerBlock)
    
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
      !        |     | iel |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element iel, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(rvector%p_rspatialDiscr, &
          IelementList(IELset:IELmax), p_Idofs)
                                   
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! To calculate the element contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = rlocalVectorAssembly%cevaluationTag
      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! If the cubature points are already initialised, do not do it again.
      ! We check this by taking a look to iinitialisedElements which
      ! gives the current maximum of initialised elements.
      if (IELmax .gt. rlocalVectorAssembly%iinitialisedElements) then

        ! (Re-)initialise!
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)

        ! Remember the new number of initialised elements
        rlocalVectorAssembly%iinitialisedElements = IELmax

      else
        ! No need.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (p_revalElementSet,&
          cevaluationTag, rvector%p_rspatialDiscr%p_rtriangulation, &
          IelementList(IELset:IELmax), rlocalVectorAssembly%ctrafoType, &
          rlocalVectorAssembly%p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => p_revalElementSet%p_Ddetj
      
      ! Now it is time to call our coefficient function to calculate the
      ! function values in the cubature points:
      if (present(fcoeff_buildVectorSc_sim)) then
        call domint_initIntegrationByEvalSet (p_revalElementSet,rintSubset)
        rintSubset%ielementDistribution = 0
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => p_Idofs
        rintSubset%celement = rlocalVectorAssembly%celement
        call fcoeff_buildVectorSc_sim (rvector%p_rspatialDiscr,&
            rlocalVectorAssembly%rform,  IELmax-IELset+1, ncubp,&
            p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
            p_Idofs, rintSubset,&
            p_Dcoefficients(:,:,1:IELmax-IELset+1), rcollection)
        call domint_doneIntegration (rintSubset)
      else
        p_Dcoefficients(:,:,1:IELmax-IELset+1) = 1.0_DP
      end if
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (rlocalVectorAssembly%celement, &
          p_revalElementSet, rlocalVectorAssembly%Bder, &
          rlocalVectorAssembly%p_Dbas)
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the DOF`s and cubature points to calculate the
      ! integral:

      do iel = 1,IELmax-IELset+1
        
        ! We make a 'local' approach, i.e. we calculate the values of the
        ! integral into the vector DlocalData and add them later into
        ! the large solution vector.
        
        ! Clear the output vector.
        DlocalData(1:indof) = 0.0_DP

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!

          domega = p_Domega(icubp)*abs(p_Ddetj(icubp,iel))

          ! Loop over the additive factors in the bilinear form.
          do ialbet = 1,rlocalVectorAssembly%rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be:
            !
            ! int_...  f * ( phi_i )_IA
            !
            ! -> IA=0: function value, 
            !      =1: first derivative, 
            !      =2: 2nd derivative,...
            !    as defined in the module 'derivative'.
            
            ia = p_Idescriptors(ialbet)
            
            ! Multiply domega with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            daux = domega * p_Dcoefficients(ialbet,icubp,iel)
          
            ! Now loop through all possible combinations of DOF`s
            ! in the current cubature point. 

            do idofe = 1,indof
              
              ! Get the value of the basis function 
              ! phi_o in the cubature point. 
              ! Them multiply:
              !    DBAS(..) * AUX
              ! ~= phi_i * coefficient * cub.weight
              ! Summing this up gives the integral, so the contribution
              ! to the vector. 
              !
              ! Simply summing up DBAS(..) * AUX would give
              ! the additive contribution for the vector. We save this
              ! contribution in the local array.
              
              DlocalData(idofe) = DlocalData(idofe)+p_Dbas(idofe,ia,icubp,iel)*daux
              
            end do ! jdofe
            
          end do ! ialbet

        end do ! icubp 
        
        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        do IDOFE = 1,indof
          p_Ddata(p_Idofs(idofe,iel)) = p_Ddata(p_Idofs(idofe,iel)) + DlocalData(idofe)
        end do

      end do ! iel

    end do ! IELset
    
    ! Release the local vector assembly structure
    call linf_releaseAssemblyData(rlocalVectorAssembly)
  
  end subroutine

  !****************************************************************************
  
!<subroutine>  
  
  subroutine linf_assembleSubmeshVecIntl (rvectorAssembly, rvector, IelementList,&
      fcoeff_buildVectorBl_sim, rcollection)

!<description>

  ! Assembles the vector entries for a submesh by integrating over the domain.

!</description>
  
!<input>
  
  ! List of elements where to assemble the bilinear form.
  integer, dimension(:), intent(in), target :: IelementList
  
  ! A callback routine which is able to calculate the values of the
  ! function $f$ which is to be discretised.
  include 'intf_coefficientVectorBl.inc'
  optional :: fcoeff_buildVectorBl_sim
  
!</input>

!<inputoutput>
  
  ! A vector assembly structure prepared with linf_initAssembly.
  type(t_linfVectorAssembly), intent(inout), target :: rvectorAssembly
  
  ! A vector where to assemble the contributions to.
  type(t_vectorScalar), intent(inout) :: rvector
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection

!</inputoutput>

!</subroutine>
  
    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: indof,ncubp
    
    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax
    integer :: iel,icubp,ialbet,ia,idofe,ivar
    real(DP) :: domega
    integer(I32) :: cevaluationTag
    type(t_linfVectorAssembly), target :: rlocalVectorAssembly
    type(t_domainIntSubset) :: rintSubset
    real(DP), dimension(:,:), pointer :: p_Ddetj
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP), dimension(:,:,:,:), pointer :: p_Dcoefficients
    integer, dimension(:),pointer :: p_Idescriptors
    integer, dimension(:,:), pointer :: p_Idofs
    type(t_evalElementSet), pointer :: p_revalElementSet

    ! A small vector holding only the additive contributions of one element
    real(DP), dimension(:,:), allocatable :: DlocalData
    real(DP), dimension(:), allocatable :: Daux
  
    ! Allocate temporal array
    allocate(DlocalData(rvector%NVAR,EL_MAXNBAS), Daux(rvector%NVAR))

    ! Get some pointers for faster access
    call lsyssc_getbase_double (rvector,p_Ddata)
    indof = rvectorAssembly%indof
    ncubp = rvectorAssembly%ncubp

    ! Copy the assembly data to the local assembly data,
    ! where we can allocate memory.
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the vector
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original assembly
    ! stucture or disturbing the data of the other processors.
    rlocalVectorAssembly = rvectorAssembly
    call linf_allocAssemblyData(rlocalVectorAssembly, rvector%NVAR)
    
    ! Get some more pointers to local data.
    p_Domega => rlocalVectorAssembly%p_Domega
    p_Dbas => rlocalVectorAssembly%p_Dbas
    p_Dcoefficients => rlocalVectorAssembly%p_Dcoefficients
    p_Idescriptors => rlocalVectorAssembly%rform%Idescriptors
    p_Idofs => rlocalVectorAssembly%p_Idofs
    p_revalElementSet => rlocalVectorAssembly%revalElementSet
        
    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !%OMP do schedule(static,1)
    do IELset = 1, size(IelementList), rlocalVectorAssembly%nelementsPerBlock
    
      ! We always handle nelementsPerBlock elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(IelementList),IELset-1+rlocalVectorAssembly%nelementsPerBlock)
    
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
      !        |     | iel |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element iel, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(rvector%p_rspatialDiscr, &
          IelementList(IELset:IELmax), p_Idofs)
                                   
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! To calculate the element contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = rlocalVectorAssembly%cevaluationTag
      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! If the cubature points are already initialised, do not do it again.
      ! We check this by taking a look to iinitialisedElements which
      ! gives the current maximum of initialised elements.
      if (IELmax .gt. rlocalVectorAssembly%iinitialisedElements) then

        ! (Re-)initialise!
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)

        ! Remember the new number of initialised elements
        rlocalVectorAssembly%iinitialisedElements = IELmax

      else
        ! No need.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (p_revalElementSet,&
          cevaluationTag, rvector%p_rspatialDiscr%p_rtriangulation, &
          IelementList(IELset:IELmax), rlocalVectorAssembly%ctrafoType, &
          rlocalVectorAssembly%p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => p_revalElementSet%p_Ddetj
      
      ! Now it is time to call our coefficient function to calculate the
      ! function values in the cubature points:
      if (present(fcoeff_buildVectorBl_sim)) then
        call domint_initIntegrationByEvalSet (p_revalElementSet,rintSubset)
        rintSubset%ielementDistribution = 0
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => p_Idofs
        rintSubset%celement = rlocalVectorAssembly%celement
        call fcoeff_buildVectorBl_sim (rvector%p_rspatialDiscr,&
            rlocalVectorAssembly%rform,  IELmax-IELset+1, ncubp,&
            p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
            p_Idofs, rintSubset,&
            p_Dcoefficients(:,:,:,1:IELmax-IELset+1), rcollection)
        call domint_doneIntegration (rintSubset)
      else
        p_Dcoefficients(:,:,:,1:IELmax-IELset+1) = 1.0_DP
      end if
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (rlocalVectorAssembly%celement, &
          p_revalElementSet, rlocalVectorAssembly%Bder, &
          rlocalVectorAssembly%p_Dbas)
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the DOF`s and cubature points to calculate the
      ! integral:

      do iel = 1,IELmax-IELset+1
        
        ! We make a 'local' approach, i.e. we calculate the values of the
        ! integral into the vector DlocalData and add them later into
        ! the large solution vector.
        
        ! Clear the output vector.
        DlocalData(:,1:indof) = 0.0_DP

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!

          domega = p_Domega(icubp)*abs(p_Ddetj(icubp,iel))

          ! Loop over the additive factors in the bilinear form.
          do ialbet = 1,rlocalVectorAssembly%rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be:
            !
            ! int_...  f * ( phi_i )_IA
            !
            ! -> IA=0: function value, 
            !      =1: first derivative, 
            !      =2: 2nd derivative,...
            !    as defined in the module 'derivative'.
            
            ia = p_Idescriptors(ialbet)
            
            ! Multiply domega with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            Daux = domega * p_Dcoefficients(:,ialbet,icubp,iel)
          
            ! Now loop through all possible combinations of DOF`s
            ! in the current cubature point. 

            do idofe = 1,indof
              
              ! Get the value of the basis function 
              ! phi_o in the cubature point. 
              ! Them multiply:
              !    DBAS(..) * AUX
              ! ~= phi_i * coefficient * cub.weight
              ! Summing this up gives the integral, so the contribution
              ! to the vector. 
              !
              ! Simply summing up DBAS(..) * AUX would give
              ! the additive contribution for the vector. We save this
              ! contribution in the local array.
              
              DlocalData(:,idofe) = DlocalData(:,idofe)+p_Dbas(idofe,ia,icubp,iel)*Daux
              
            end do ! jdofe
            
          end do ! ialbet

        end do ! icubp 
        
        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        do IDOFE = 1,indof
          do ivar = 1, rvector%NVAR
            p_Ddata(rvector%NVAR*(p_Idofs(idofe,iel)-1)+ivar) =&
                p_Ddata(rvector%NVAR*(p_Idofs(idofe,iel)-1)+ivar) + DlocalData(ivar,idofe)
          end do
        end do

      end do ! iel

    end do ! IELset
    
    ! Release the local vector assembly structure
    call linf_releaseAssemblyData(rlocalVectorAssembly)
  
    ! Deallocate memory
    deallocate(DlocalData, Daux)
    
  end subroutine

  !****************************************************************************
  
!<subroutine>  
  
  subroutine linf_assembleSubmeshVectorBdr2D (rvectorAssembly, rvector,&
      rboundaryRegion, IelementList, IelementOrientation, DedgePosition,&
      fcoeff_buildVectorScBdr2D_sim, rcollection)

!<description>

  ! Assembles the vector entries for a submesh by integration over the boundary region.

!</description>

!<input>
  
  ! A boundary region where to assemble the contribution
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! List of elements where to assemble the linear form.
  integer, dimension(:), intent(in), target :: IelementList
  
  ! List of element orientations where to assemble the linear form.
  integer, dimension(:), intent(in) :: IelementOrientation

  ! List of start- and end-parameter values of the edges on the boundary
  real(DP), dimension(:,:), intent(in) :: DedgePosition

  ! A callback routine which is able to calculate the values of the
  ! function $f$ which is to be discretised.
  include 'intf_coefficientVectorScBdr2D.inc'
  optional :: fcoeff_buildVectorScBdr2D_sim 
  
!</input>

!<inputoutput>
  
  ! A vector assembly structure prepared with linf_initAssembly.
  type(t_linfVectorAssembly), intent(inout), target :: rvectorAssembly
  
  ! A vector where to assemble the contributions to.
  type(t_vectorScalar), intent(inout) :: rvector  
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection

!</inputoutput>
  
!</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: indof,ncubp
    
    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax,ibdc,k
    integer :: iel,icubp,ialbet,ia,idofe
    real(DP) :: domega,daux,dlen
    integer(I32) :: cevaluationTag
    type(t_linfVectorAssembly), target :: rlocalVectorAssembly
    type(t_domainIntSubset) :: rintSubset
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    integer, dimension(:),pointer :: p_Idescriptors
    integer, dimension(:,:), pointer :: p_Idofs
    type(t_evalElementSet), pointer :: p_revalElementSet

    ! A small vector holding only the additive contributions of
    ! one element
    real(DP), dimension(EL_MAXNBAS) :: DlocalData
  
    ! Arrays for cubature points 1D->2D
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    real(DP), dimension(:,:,:), allocatable :: Dxi2D,DpointsRef
    real(DP), dimension(:,:), allocatable :: DpointsPar
    
    integer(i32) :: icoordSystem

    ! Boundary component?
    ibdc = rboundaryRegion%iboundCompIdx

    ! Get some pointers for faster access
    call lsyssc_getbase_double (rvector, p_Ddata)
    indof = rvectorAssembly%indof
    ncubp = rvectorAssembly%ncubp

    ! Copy the assembly data to the local assembly data,
    ! where we can allocate memory.
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the vector
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original assembly
    ! stucture or disturbing the data of the other processors.
    rlocalVectorAssembly = rvectorAssembly
    call linf_allocAssemblyData(rlocalVectorAssembly)
    
    ! Get some more pointers to local data.
    p_Domega => rlocalVectorAssembly%p_Domega
    p_Dbas => rlocalVectorAssembly%p_Dbas
    p_Dcoefficients => rlocalVectorAssembly%p_Dcoefficients(1,:,:,:)
    p_DcubPtsRef => rlocalVectorAssembly%p_DcubPtsRef
    p_Idescriptors => rlocalVectorAssembly%rform%Idescriptors
    p_Idofs => rlocalVectorAssembly%p_Idofs
    p_revalElementSet => rlocalVectorAssembly%revalElementSet
    
    ! Transpose the coordinate array such that we get coordinates we
    ! can work with in the mapping between 1D and 2D.
    do k = 1, ubound(p_DcubPtsRef,1)
      do icubp = 1,ncubp
        Dxi1D(icubp,k) = p_DcubPtsRef(k,icubp)
      end do
    end do

    ! Allocate memory for the cubature points in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,rlocalVectorAssembly%nelementsPerBlock))

    ! Allocate memory for the coordinates of the reference points
    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalVectorAssembly%nelementsPerBlock))

    ! Allocate memory for the parameter values of the points on the boundary
    allocate(DpointsPar(ncubp,rlocalVectorAssembly%nelementsPerBlock))

    ! Get the type of coordinate system
    icoordSystem = elem_igetCoordSystem(rlocalVectorAssembly%celement)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !%OMP do schedule(static,1)
    do IELset = 1, size(IelementList), rlocalVectorAssembly%nelementsPerBlock
    
      ! We always handle nelementsPerBlock elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(IelementList),IELset-1+rlocalVectorAssembly%nelementsPerBlock)
      
      ! Map the 1D cubature points to the edges in 2D.
      do iel = 1,IELmax-IELset+1
        call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(IELset+iel-1), &
            ncubp, Dxi1D, Dxi2D(:,:,iel))
      end do

      ! Calculate the parameter values of the points
      do iel = 1,IELmax-IELset+1
        do icubp = 1,ncubp
          ! Dxi1D is in [-1,1] while the current edge has parmeter values
          ! [DedgePosition(1),DedgePosition(2)]. So do a linear
          ! transformation to transform Dxi1D into that interval, this 
          ! gives the parameter values in length parametrisation
          call mprim_linearRescale(Dxi1D(icubp,1), -1.0_DP, 1.0_DP,&
              DedgePosition(1,IELset+iel-1), DedgePosition(2,IELset+iel-1),&
              DpointsPar(icubp,iel))
        end do
      end do
      
      ! Transpose the coordinate array such that we get coordinates we
      ! can work with.
      do iel = 1,IELmax-IELset+1
        do icubp = 1,ncubp
          do k = 1,ubound(DpointsRef,1)
            DpointsRef(k,icubp,iel) = Dxi2D(icubp,k,iel)
          end do
        end do
      end do
      
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
      !        |     | iel |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element iel, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(rvector%p_rspatialDiscr, &
          IelementList(IELset:IELmax), p_Idofs)
                                   
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! To calculate the element contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = rlocalVectorAssembly%cevaluationTag
      
      ! The cubature points are already initialised by 1D->2D mapping.
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

      ! Calculate all information that is necessary to evaluate the
      ! finite element on all cells of our subset. This includes the
      ! coordinates of the points on the cells.
      call elprep_prepareSetForEvaluation (p_revalElementSet,&
          cevaluationTag, rvector%p_rspatialDiscr%p_rtriangulation, &
          IelementList(IELset:IELmax), rlocalVectorAssembly%ctrafoType, &
          DpointsRef=DpointsRef)
      
      ! Now it is time to call our coefficient function to calculate the
      ! function values in the cubature points:
      if (present(fcoeff_buildVectorScBdr2D_sim)) then
        call domint_initIntegrationByEvalSet (p_revalElementSet, rintSubset)
        rintSubset%ielementDistribution = 0
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => p_Idofs
        rintSubset%celement = rlocalVectorAssembly%celement
        call fcoeff_buildVectorScBdr2D_sim (rvector%p_rspatialDiscr,&
            rlocalVectorAssembly%rform,  IELmax-IELset+1, ncubp,&
            p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
            ibdc, DpointsPar(:,1:IELmax-IELset+1),&
            p_Idofs, rintSubset, &
            p_Dcoefficients(:,:,1:IELmax-IELset+1), rcollection)
        call domint_doneIntegration (rintSubset)
      else
        p_Dcoefficients(:,:,1:IELmax-IELset+1) = 1.0_DP
      end if
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (rlocalVectorAssembly%celement, &
          p_revalElementSet, rlocalVectorAssembly%Bder, &
          rlocalVectorAssembly%p_Dbas)
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the DOF`s and cubature points to calculate the
      ! integral:

      do iel = 1,IELmax-IELset+1
        
        ! We make a 'local' approach, i.e. we calculate the values of the
        ! integral into the vector DlocalData and add them later into
        ! the large solution vector.
        
        ! Clear the output vector.
        DlocalData(1:indof) = 0.0_DP

        ! Get the length of the edge. Let us use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit interval [-1,1] has length 2.
        dlen = 0.5_DP*(DedgePosition(2,IELset+iel-1)-DedgePosition(1,IELset+iel-1))

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature
          ! formula in that cubature point.

          domega = dlen * p_Domega(icubp)
          
          ! Loop over the additive factors in the bilinear form.
          do ialbet = 1,rlocalVectorAssembly%rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be:
            !
            ! int_...  f * ( phi_i )_IA
            !
            ! -> IA=0: function value, 
            !      =1: first derivative, 
            !      =2: 2nd derivative,...
            !    as defined in the module 'derivative'.
            
            ia = p_Idescriptors(ialbet)
            
            ! Multiply domega with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            daux = domega * p_Dcoefficients(ialbet,icubp,iel)
          
            ! Now loop through all possible combinations of DOF`s
            ! in the current cubature point. 

            do idofe = 1,indof
              
              ! Get the value of the basis function 
              ! phi_o in the cubature point. 
              ! Them multiply:
              !    DBAS(..) * AUX
              ! ~= phi_i * coefficient * cub.weight
              ! Summing this up gives the integral, so the contribution
              ! to the vector. 
              !
              ! Simply summing up DBAS(..) * AUX would give
              ! the additive contribution for the vector. We save this
              ! contribution in the local array.
              
              DlocalData(idofe) = DlocalData(idofe)+p_Dbas(idofe,ia,icubp,iel)*daux
              
            end do ! idofe
            
          end do ! ialbet

        end do ! icubp 
        
        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        do idofe = 1,indof
          p_Ddata(p_Idofs(idofe,iel)) = p_Ddata(p_Idofs(idofe,iel)) + DlocalData(idofe)
        end do

      end do ! iel

    end do ! IELset
    
    ! Release the local vector assembly structure
    call linf_releaseAssemblyData(rlocalVectorAssembly)

    ! Deallocate memory
    deallocate(Dxi2D, DpointsRef, DpointsPar)

  end subroutine

  !****************************************************************************
  
!<subroutine>  
  
  subroutine linf_assembleSubmeshVecIntlBdr2D (rvectorAssembly, rvector,&
      rboundaryRegion, IelementList, IelementOrientation, DedgePosition,&
      fcoeff_buildVectorBlBdr2D_sim, rcollection)

!<description>

  ! Assembles the vector entries for a submesh by integration over the boundary region.

!</description>

!<input>
  
  ! A boundary region where to assemble the contribution
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! List of elements where to assemble the linear form.
  integer, dimension(:), intent(in), target :: IelementList
  
  ! List of element orientations where to assemble the linear form.
  integer, dimension(:), intent(in) :: IelementOrientation

  ! List of start- and end-parameter values of the edges on the boundary
  real(DP), dimension(:,:), intent(in) :: DedgePosition

  ! A callback routine which is able to calculate the values of the
  ! function $f$ which is to be discretised.
  include 'intf_coefficientVectorBlBdr2D.inc'
  optional :: fcoeff_buildVectorBlBdr2D_sim 
  
!</input>

!<inputoutput>
  
  ! A vector assembly structure prepared with linf_initAssembly.
  type(t_linfVectorAssembly), intent(inout), target :: rvectorAssembly
  
  ! A vector where to assemble the contributions to.
  type(t_vectorScalar), intent(inout) :: rvector  
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection

!</inputoutput>
  
!</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: indof,ncubp
    
    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax,ibdc,k
    integer :: iel,icubp,ialbet,ia,idofe,ivar
    real(DP) :: domega,dlen
    integer(I32) :: cevaluationTag
    type(t_linfVectorAssembly), target :: rlocalVectorAssembly
    type(t_domainIntSubset) :: rintSubset
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP), dimension(:,:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    integer, dimension(:),pointer :: p_Idescriptors
    integer, dimension(:,:), pointer :: p_Idofs
    type(t_evalElementSet), pointer :: p_revalElementSet

    ! A small vector holding only the additive contributions of one element
    real(DP), dimension(:,:), allocatable :: DlocalData
    real(DP), dimension(:), allocatable :: Daux
    
    ! Arrays for cubature points 1D->2D
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    real(DP), dimension(:,:,:), allocatable :: Dxi2D,DpointsRef
    real(DP), dimension(:,:), allocatable :: DpointsPar
    
    integer(i32) :: icoordSystem

    ! Boundary component?
    ibdc = rboundaryRegion%iboundCompIdx

    ! Allocate temporal array
    allocate(DlocalData(rvector%NVAR,EL_MAXNBAS), Daux(rvector%NVAR))

    ! Get some pointers for faster access
    call lsyssc_getbase_double (rvector, p_Ddata)
    indof = rvectorAssembly%indof
    ncubp = rvectorAssembly%ncubp

    ! Copy the assembly data to the local assembly data,
    ! where we can allocate memory.
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the vector
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original assembly
    ! stucture or disturbing the data of the other processors.
    rlocalVectorAssembly = rvectorAssembly
    call linf_allocAssemblyData(rlocalVectorAssembly, rvector%NVAR)
    
    ! Get some more pointers to local data.
    p_Domega => rlocalVectorAssembly%p_Domega
    p_Dbas => rlocalVectorAssembly%p_Dbas
    p_Dcoefficients => rlocalVectorAssembly%p_Dcoefficients
    p_DcubPtsRef => rlocalVectorAssembly%p_DcubPtsRef
    p_Idescriptors => rlocalVectorAssembly%rform%Idescriptors
    p_Idofs => rlocalVectorAssembly%p_Idofs
    p_revalElementSet => rlocalVectorAssembly%revalElementSet
    
    ! Transpose the coordinate array such that we get coordinates we
    ! can work with in the mapping between 1D and 2D.
    do k = 1, ubound(p_DcubPtsRef,1)
      do icubp = 1,ncubp
        Dxi1D(icubp,k) = p_DcubPtsRef(k,icubp)
      end do
    end do

    ! Allocate memory for the cubature points in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,rlocalVectorAssembly%nelementsPerBlock))

    ! Allocate memory for the coordinates of the reference points
    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalVectorAssembly%nelementsPerBlock))

    ! Allocate memory for the parameter values of the points on the boundary
    allocate(DpointsPar(ncubp,rlocalVectorAssembly%nelementsPerBlock))

    ! Get the type of coordinate system
    icoordSystem = elem_igetCoordSystem(rlocalVectorAssembly%celement)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !%OMP do schedule(static,1)
    do IELset = 1, size(IelementList), rlocalVectorAssembly%nelementsPerBlock
    
      ! We always handle nelementsPerBlock elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(IelementList),IELset-1+rlocalVectorAssembly%nelementsPerBlock)
      
      ! Map the 1D cubature points to the edges in 2D.
      do iel = 1,IELmax-IELset+1
        call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(IELset+iel-1), &
            ncubp, Dxi1D, Dxi2D(:,:,iel))
      end do

      ! Calculate the parameter values of the points
      do iel = 1,IELmax-IELset+1
        do icubp = 1,ncubp
          ! Dxi1D is in [-1,1] while the current edge has parmeter values
          ! [DedgePosition(1),DedgePosition(2)]. So do a linear
          ! transformation to transform Dxi1D into that interval, this 
          ! gives the parameter values in length parametrisation
          call mprim_linearRescale(Dxi1D(icubp,1), -1.0_DP, 1.0_DP,&
              DedgePosition(1,IELset+iel-1), DedgePosition(2,IELset+iel-1),&
              DpointsPar(icubp,iel))
        end do
      end do
      
      ! Transpose the coordinate array such that we get coordinates we
      ! can work with.
      do iel = 1,IELmax-IELset+1
        do icubp = 1,ncubp
          do k = 1,ubound(DpointsRef,1)
            DpointsRef(k,icubp,iel) = Dxi2D(icubp,k,iel)
          end do
        end do
      end do
      
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
      !        |     | iel |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element iel, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(rvector%p_rspatialDiscr, &
          IelementList(IELset:IELmax), p_Idofs)
                                   
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! To calculate the element contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = rlocalVectorAssembly%cevaluationTag
      
      ! The cubature points are already initialised by 1D->2D mapping.
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

      ! Calculate all information that is necessary to evaluate the
      ! finite element on all cells of our subset. This includes the
      ! coordinates of the points on the cells.
      call elprep_prepareSetForEvaluation (p_revalElementSet,&
          cevaluationTag, rvector%p_rspatialDiscr%p_rtriangulation, &
          IelementList(IELset:IELmax), rlocalVectorAssembly%ctrafoType, &
          DpointsRef=DpointsRef)
      
      ! Now it is time to call our coefficient function to calculate the
      ! function values in the cubature points:
      if (present(fcoeff_buildVectorBlBdr2D_sim)) then
        call domint_initIntegrationByEvalSet (p_revalElementSet, rintSubset)
        rintSubset%ielementDistribution = 0
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => p_Idofs
        rintSubset%celement = rlocalVectorAssembly%celement
        call fcoeff_buildVectorBlBdr2D_sim (rvector%p_rspatialDiscr,&
            rlocalVectorAssembly%rform,  IELmax-IELset+1, ncubp,&
            p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
            ibdc, DpointsPar(:,1:IELmax-IELset+1),&
            p_Idofs, rintSubset, &
            p_Dcoefficients(:,:,:,1:IELmax-IELset+1), rcollection)
        call domint_doneIntegration (rintSubset)
      else
        p_Dcoefficients(:,:,:,1:IELmax-IELset+1) = 1.0_DP
      end if
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (rlocalVectorAssembly%celement, &
          p_revalElementSet, rlocalVectorAssembly%Bder, &
          rlocalVectorAssembly%p_Dbas)
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the DOF`s and cubature points to calculate the
      ! integral:

      do iel = 1,IELmax-IELset+1
        
        ! We make a 'local' approach, i.e. we calculate the values of the
        ! integral into the vector DlocalData and add them later into
        ! the large solution vector.
        
        ! Clear the output vector.
        DlocalData(:,1:indof) = 0.0_DP

        ! Get the length of the edge. Let us use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit interval [-1,1] has length 2.
        dlen = 0.5_DP*(DedgePosition(2,IELset+iel-1)-DedgePosition(1,IELset+iel-1))

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature
          ! formula in that cubature point.

          domega = dlen * p_Domega(icubp)
          
          ! Loop over the additive factors in the bilinear form.
          do ialbet = 1,rlocalVectorAssembly%rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be:
            !
            ! int_...  f * ( phi_i )_IA
            !
            ! -> IA=0: function value, 
            !      =1: first derivative, 
            !      =2: 2nd derivative,...
            !    as defined in the module 'derivative'.
            
            ia = p_Idescriptors(ialbet)
            
            ! Multiply domega with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            Daux = domega * p_Dcoefficients(:,ialbet,icubp,iel)
          
            ! Now loop through all possible combinations of DOF`s
            ! in the current cubature point. 

            do idofe = 1,indof
              
              ! Get the value of the basis function 
              ! phi_o in the cubature point. 
              ! Them multiply:
              !    DBAS(..) * AUX
              ! ~= phi_i * coefficient * cub.weight
              ! Summing this up gives the integral, so the contribution
              ! to the vector. 
              !
              ! Simply summing up DBAS(..) * AUX would give
              ! the additive contribution for the vector. We save this
              ! contribution in the local array.
              
              DlocalData(:,idofe) = DlocalData(:,idofe)+p_Dbas(idofe,ia,icubp,iel)*Daux
              
            end do ! idofe
            
          end do ! ialbet

        end do ! icubp 
        
        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        do idofe = 1,indof
          do ivar = 1, rvector%NVAR
            p_Ddata(rvector%NVAR*(p_Idofs(idofe,iel)-1)+ivar) =&
                p_Ddata(rvector%NVAR*(p_Idofs(idofe,iel)-1)+ivar) + DlocalData(ivar,idofe)
          end do
        end do

      end do ! iel

    end do ! IELset
    
    ! Release the local vector assembly structure
    call linf_releaseAssemblyData(rlocalVectorAssembly)

    ! Deallocate memory
    deallocate(Dxi2D, DpointsRef, DpointsPar, DlocalData, Daux)

  end subroutine

  !****************************************************************************
  
!<subroutine>  
  
  subroutine linf_assembleSubmeshVectorBlock (rvectorAssembly, rvector,&
      IelementList, fcoeff_buildVectorBl_sim, rcollection)

!<description>

  ! Assembles the vector entries for a submesh by integrating over the domain.

!</description>
  
!<input>
  
  ! List of elements where to assemble the bilinear form.
  integer, dimension(:), intent(in), target :: IelementList
  
  ! A callback routine which is able to calculate the values of the
  ! vector-valued function $f$ which is to be discretised.
  include 'intf_coefficientVectorBl.inc'
  optional :: fcoeff_buildVectorBl_sim
  
!</input>

!<inputoutput>
  
  ! A vector assembly structure prepared with linf_initAssembly.
  type(t_linfVectorAssembly), intent(inout), target :: rvectorAssembly
  
  ! A block vector where to assemble the contributions to.
  type(t_vectorBlock), intent(inout) :: rvector
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection

!</inputoutput>

!</subroutine>
  
    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: indof,ncubp
    
    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax
    integer :: iel,icubp,ialbet,ia,idofe,iblock
    real(DP) :: domega
    integer(I32) :: cevaluationTag
    type(t_linfVectorAssembly), target :: rlocalVectorAssembly
    type(t_domainIntSubset) :: rintSubset
    real(DP), dimension(:,:), pointer :: p_Ddetj
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP), dimension(:,:,:,:), pointer :: p_Dcoefficients
    integer, dimension(:),pointer :: p_Idescriptors
    integer, dimension(:,:), pointer :: p_Idofs
    type(t_evalElementSet), pointer :: p_revalElementSet

    ! A small vector holding only the additive contributions of one element
    real(DP), dimension(:,:), allocatable :: DlocalData
    real(DP), dimension(:), allocatable :: Daux
  
    ! Allocate temporal array
    allocate(DlocalData(rvector%nblocks,EL_MAXNBAS), Daux(rvector%nblocks))

    ! Get some pointers for faster access
    call lsysbl_getbase_double (rvector,p_Ddata)
    indof = rvectorAssembly%indof
    ncubp = rvectorAssembly%ncubp

    ! Copy the assembly data to the local assembly data,
    ! where we can allocate memory.
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the vector
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original assembly
    ! stucture or disturbing the data of the other processors.
    rlocalVectorAssembly = rvectorAssembly
    call linf_allocAssemblyData(rlocalVectorAssembly, rvector%nblocks)
    
    ! Get some more pointers to local data.
    p_Domega => rlocalVectorAssembly%p_Domega
    p_Dbas => rlocalVectorAssembly%p_Dbas
    p_Dcoefficients => rlocalVectorAssembly%p_Dcoefficients
    p_Idescriptors => rlocalVectorAssembly%rform%Idescriptors
    p_Idofs => rlocalVectorAssembly%p_Idofs
    p_revalElementSet => rlocalVectorAssembly%revalElementSet

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !%OMP do schedule(static,1)
    do IELset = 1, size(IelementList), rlocalVectorAssembly%nelementsPerBlock
    
      ! We always handle nelementsPerBlock elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(IelementList),IELset-1+rlocalVectorAssembly%nelementsPerBlock)
    
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
      !        |     | iel |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element iel, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(&
          rvector%p_rblockDiscr%RspatialDiscr(1), &
          IelementList(IELset:IELmax), p_Idofs)
                                   
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! To calculate the element contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = rlocalVectorAssembly%cevaluationTag
      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! If the cubature points are already initialised, do not do it again.
      ! We check this by taking a look to iinitialisedElements which
      ! gives the current maximum of initialised elements.
      if (IELmax .gt. rlocalVectorAssembly%iinitialisedElements) then

        ! (Re-)initialise!
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)

        ! Remember the new number of initialised elements
        rlocalVectorAssembly%iinitialisedElements = IELmax

      else
        ! No need.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (p_revalElementSet,&
          cevaluationTag, rvector%p_rblockDiscr%RspatialDiscr(1)%p_rtriangulation, &
          IelementList(IELset:IELmax), rlocalVectorAssembly%ctrafoType, &
          rlocalVectorAssembly%p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => p_revalElementSet%p_Ddetj
      
      ! Now it is time to call our coefficient function to calculate the
      ! function values in the cubature points:
      if (present(fcoeff_buildVectorBl_sim)) then
        call domint_initIntegrationByEvalSet (p_revalElementSet,rintSubset)
        rintSubset%ielementDistribution = 0
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => p_Idofs
        rintSubset%celement = rlocalVectorAssembly%celement
        call fcoeff_buildVectorBl_sim (rvector%p_rblockDiscr%RspatialDiscr(1),&
            rlocalVectorAssembly%rform,  IELmax-IELset+1, ncubp,&
            p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
            p_Idofs, rintSubset,&
            p_Dcoefficients(:,:,:,1:IELmax-IELset+1), rcollection)
        call domint_doneIntegration (rintSubset)
      else
        p_Dcoefficients(:,:,:,1:IELmax-IELset+1) = 1.0_DP
      end if
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (rlocalVectorAssembly%celement, &
          p_revalElementSet, rlocalVectorAssembly%Bder, &
          rlocalVectorAssembly%p_Dbas)
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the DOF`s and cubature points to calculate the
      ! integral:

      do iel = 1,IELmax-IELset+1
        
        ! We make a 'local' approach, i.e. we calculate the values of the
        ! integral into the vector DlocalData and add them later into
        ! the large solution vector.
        
        ! Clear the output vector.
        DlocalData(:,1:indof) = 0.0_DP

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!

          domega = p_Domega(icubp)*abs(p_Ddetj(icubp,iel))

          ! Loop over the additive factors in the bilinear form.
          do ialbet = 1,rlocalVectorAssembly%rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be:
            !
            ! int_...  f * ( phi_i )_IA
            !
            ! -> IA=0: function value, 
            !      =1: first derivative, 
            !      =2: 2nd derivative,...
            !    as defined in the module 'derivative'.
            
            ia = p_Idescriptors(ialbet)
            
            ! Multiply domega with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            Daux = domega * p_Dcoefficients(:,ialbet,icubp,iel)
          
            ! Now loop through all possible combinations of DOF`s
            ! in the current cubature point. 

            do idofe = 1,indof
              
              ! Get the value of the basis function 
              ! phi_o in the cubature point. 
              ! Them multiply:
              !    DBAS(..) * AUX
              ! ~= phi_i * coefficient * cub.weight
              ! Summing this up gives the integral, so the contribution
              ! to the vector. 
              !
              ! Simply summing up DBAS(..) * AUX would give
              ! the additive contribution for the vector. We save this
              ! contribution in the local array.
              
              DlocalData(:,idofe) = DlocalData(:,idofe)+p_Dbas(idofe,ia,icubp,iel)*Daux
              
            end do ! jdofe
            
          end do ! ialbet

        end do ! icubp 
        
        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        do iblock =1,rvector%nblocks
          do IDOFE = 1,indof
            p_Ddata(rvector%RvectorBlock(iblock)%iidxFirstEntry+p_Idofs(idofe,iel)-1)=&
                p_Ddata(rvector%RvectorBlock(iblock)%iidxFirstEntry+p_Idofs(idofe,iel)-1)+&
                DlocalData(iblock,idofe)
          end do
        end do
        
      end do ! iel

    end do ! IELset
    
    ! Release the local vector assembly structure
    call linf_releaseAssemblyData(rlocalVectorAssembly)
  
    ! Deallocate memory
    deallocate(DlocalData, Daux)

  end subroutine

  !****************************************************************************
  
!<subroutine>  
  
  subroutine linf_assembleSubmeshVectorBlockBdr2D (rvectorAssembly, rvector,&
      rboundaryRegion, IelementList, IelementOrientation, DedgePosition,&
      fcoeff_buildVectorBlBdr2D_sim, rcollection)

!<description>

  ! Assembles the vector entries for a submesh by integration over the boundary region.

!</description>

!<input>
  
  ! A boundary region where to assemble the contribution
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! List of elements where to assemble the linear form.
  integer, dimension(:), intent(in), target :: IelementList
  
  ! List of element orientations where to assemble the linear form.
  integer, dimension(:), intent(in) :: IelementOrientation

  ! List of start- and end-parameter values of the edges on the boundary
  real(DP), dimension(:,:), intent(in) :: DedgePosition

  ! A callback routine which is able to calculate the values of the
  ! function $f$ which is to be discretised.
  include 'intf_coefficientVectorBlBdr2D.inc'
  optional :: fcoeff_buildVectorBlBdr2D_sim 
  
!</input>

!<inputoutput>
  
  ! A vector assembly structure prepared with linf_initAssembly.
  type(t_linfVectorAssembly), intent(inout), target :: rvectorAssembly
  
  ! A vector where to assemble the contributions to.
  type(t_vectorBlock), intent(inout) :: rvector  
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection

!</inputoutput>
  
!</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: indof,ncubp
    
    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax,ibdc,k
    integer :: iel,icubp,ialbet,ia,idofe,iblock
    real(DP) :: domega,dlen
    integer(I32) :: cevaluationTag
    type(t_linfVectorAssembly), target :: rlocalVectorAssembly
    type(t_domainIntSubset) :: rintSubset
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP), dimension(:,:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    integer, dimension(:),pointer :: p_Idescriptors
    integer, dimension(:,:), pointer :: p_Idofs
    type(t_evalElementSet), pointer :: p_revalElementSet

    ! A small vector holding only the additive contributions of one element
    real(DP), dimension(:,:), allocatable :: DlocalData
    real(DP), dimension(:), allocatable :: Daux
    
    ! Arrays for cubature points 1D->2D
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    real(DP), dimension(:,:,:), allocatable :: Dxi2D,DpointsRef
    real(DP), dimension(:,:), allocatable :: DpointsPar
    
    integer(i32) :: icoordSystem

    ! Boundary component?
    ibdc = rboundaryRegion%iboundCompIdx

    ! Allocate temporal array
    allocate(DlocalData(rvector%nblocks,EL_MAXNBAS), Daux(rvector%nblocks))

    ! Get some pointers for faster access
    call lsysbl_getbase_double (rvector, p_Ddata)
    indof = rvectorAssembly%indof
    ncubp = rvectorAssembly%ncubp

    ! Copy the assembly data to the local assembly data,
    ! where we can allocate memory.
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the vector
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original assembly
    ! stucture or disturbing the data of the other processors.
    rlocalVectorAssembly = rvectorAssembly
    call linf_allocAssemblyData(rlocalVectorAssembly, rvector%nblocks)
    
    ! Get some more pointers to local data.
    p_Domega => rlocalVectorAssembly%p_Domega
    p_Dbas => rlocalVectorAssembly%p_Dbas
    p_Dcoefficients => rlocalVectorAssembly%p_Dcoefficients
    p_DcubPtsRef => rlocalVectorAssembly%p_DcubPtsRef
    p_Idescriptors => rlocalVectorAssembly%rform%Idescriptors
    p_Idofs => rlocalVectorAssembly%p_Idofs
    p_revalElementSet => rlocalVectorAssembly%revalElementSet
    
    ! Transpose the coordinate array such that we get coordinates we
    ! can work with in the mapping between 1D and 2D.
    do k = 1, ubound(p_DcubPtsRef,1)
      do icubp = 1,ncubp
        Dxi1D(icubp,k) = p_DcubPtsRef(k,icubp)
      end do
    end do

    ! Allocate memory for the cubature points in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,rlocalVectorAssembly%nelementsPerBlock))

    ! Allocate memory for the coordinates of the reference points
    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalVectorAssembly%nelementsPerBlock))

    ! Allocate memory for the parameter values of the points on the boundary
    allocate(DpointsPar(ncubp,rlocalVectorAssembly%nelementsPerBlock))

    ! Get the type of coordinate system
    icoordSystem = elem_igetCoordSystem(rlocalVectorAssembly%celement)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !%OMP do schedule(static,1)
    do IELset = 1, size(IelementList), rlocalVectorAssembly%nelementsPerBlock
    
      ! We always handle nelementsPerBlock elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(IelementList),IELset-1+rlocalVectorAssembly%nelementsPerBlock)
      
      ! Map the 1D cubature points to the edges in 2D.
      do iel = 1,IELmax-IELset+1
        call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(IELset+iel-1), &
            ncubp, Dxi1D, Dxi2D(:,:,iel))
      end do

      ! Calculate the parameter values of the points
      do iel = 1,IELmax-IELset+1
        do icubp = 1,ncubp
          ! Dxi1D is in [-1,1] while the current edge has parmeter values
          ! [DedgePosition(1),DedgePosition(2)]. So do a linear
          ! transformation to transform Dxi1D into that interval, this 
          ! gives the parameter values in length parametrisation
          call mprim_linearRescale(Dxi1D(icubp,1), -1.0_DP, 1.0_DP,&
              DedgePosition(1,IELset+iel-1), DedgePosition(2,IELset+iel-1),&
              DpointsPar(icubp,iel))
        end do
      end do
      
      ! Transpose the coordinate array such that we get coordinates we
      ! can work with.
      do iel = 1,IELmax-IELset+1
        do icubp = 1,ncubp
          do k = 1,ubound(DpointsRef,1)
            DpointsRef(k,icubp,iel) = Dxi2D(icubp,k,iel)
          end do
        end do
      end do
      
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
      !        |     | iel |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element iel, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(&
          rvector%p_rblockDiscr%RspatialDiscr(1), &
          IelementList(IELset:IELmax), p_Idofs)
                                   
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! To calculate the element contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = rlocalVectorAssembly%cevaluationTag
      
      ! The cubature points are already initialised by 1D->2D mapping.
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

      ! Calculate all information that is necessary to evaluate the
      ! finite element on all cells of our subset. This includes the
      ! coordinates of the points on the cells.
      call elprep_prepareSetForEvaluation (p_revalElementSet,&
          cevaluationTag, rvector%p_rblockDiscr%RspatialDiscr(1)%p_rtriangulation, &
          IelementList(IELset:IELmax), rlocalVectorAssembly%ctrafoType, &
          DpointsRef=DpointsRef)
      
      ! Now it is time to call our coefficient function to calculate the
      ! function values in the cubature points:
      if (present(fcoeff_buildVectorBlBdr2D_sim)) then
        call domint_initIntegrationByEvalSet (p_revalElementSet, rintSubset)
        rintSubset%ielementDistribution = 0
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => p_Idofs
        rintSubset%celement = rlocalVectorAssembly%celement
        call fcoeff_buildVectorBlBdr2D_sim (rvector%p_rblockDiscr%RspatialDiscr(1),&
            rlocalVectorAssembly%rform,  IELmax-IELset+1, ncubp,&
            p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
            ibdc, DpointsPar(:,1:IELmax-IELset+1),&
            p_Idofs, rintSubset, &
            p_Dcoefficients(:,:,:,1:IELmax-IELset+1), rcollection)
        call domint_doneIntegration (rintSubset)
      else
        p_Dcoefficients(:,:,:,1:IELmax-IELset+1) = 1.0_DP
      end if
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (rlocalVectorAssembly%celement, &
          p_revalElementSet, rlocalVectorAssembly%Bder, &
          rlocalVectorAssembly%p_Dbas)
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the DOF`s and cubature points to calculate the
      ! integral:

      do iel = 1,IELmax-IELset+1
        
        ! We make a 'local' approach, i.e. we calculate the values of the
        ! integral into the vector DlocalData and add them later into
        ! the large solution vector.
        
        ! Clear the output vector.
        DlocalData(:,1:indof) = 0.0_DP

        ! Get the length of the edge. Let us use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit interval [-1,1] has length 2.
        dlen = 0.5_DP*(DedgePosition(2,IELset+iel-1)-DedgePosition(1,IELset+iel-1))

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature
          ! formula in that cubature point.

          domega = dlen * p_Domega(icubp)
          
          ! Loop over the additive factors in the bilinear form.
          do ialbet = 1,rlocalVectorAssembly%rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be:
            !
            ! int_...  f * ( phi_i )_IA
            !
            ! -> IA=0: function value, 
            !      =1: first derivative, 
            !      =2: 2nd derivative,...
            !    as defined in the module 'derivative'.
            
            ia = p_Idescriptors(ialbet)
            
            ! Multiply domega with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            Daux = domega * p_Dcoefficients(:,ialbet,icubp,iel)
          
            ! Now loop through all possible combinations of DOF`s
            ! in the current cubature point. 

            do idofe = 1,indof
              
              ! Get the value of the basis function 
              ! phi_o in the cubature point. 
              ! Them multiply:
              !    DBAS(..) * AUX
              ! ~= phi_i * coefficient * cub.weight
              ! Summing this up gives the integral, so the contribution
              ! to the vector. 
              !
              ! Simply summing up DBAS(..) * AUX would give
              ! the additive contribution for the vector. We save this
              ! contribution in the local array.
              
              DlocalData(:,idofe) = DlocalData(:,idofe)+p_Dbas(idofe,ia,icubp,iel)*Daux
              
            end do ! idofe
            
          end do ! ialbet

        end do ! icubp 
        
        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        do iblock =1,rvector%nblocks
          do idofe = 1,indof
            p_Ddata(rvector%RvectorBlock(iblock)%iidxFirstEntry+p_Idofs(idofe,iel)-1)=&
                p_Ddata(rvector%RvectorBlock(iblock)%iidxFirstEntry+p_Idofs(idofe,iel)-1)+&
                DlocalData(iblock,idofe)
          end do
        end do
        
      end do ! iel

    end do ! IELset
    
    ! Release the local vector assembly structure
    call linf_releaseAssemblyData(rlocalVectorAssembly)

    ! Deallocate memory
    deallocate(Dxi2D, DpointsRef, DpointsPar, DlocalData, Daux)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine linf_buildVectorScalar2 (rform, bclear, rvectorScalar,&
                                      fcoeff_buildVectorSc_sim, rcollection)
  
!<description>
  ! This routine assembles the entries of a vector according to a linear form
  ! (typically used for assembling RHS vectors).
  !
  ! If bclear=TRUE, the vector is cleared before the assembly and any 
  ! sorting of the entries is switched off - the vector is set up unsorted.
  !
  ! If bclear=FALSE, the vector must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
  !
  ! IMPLEMENTATIONAL REMARK:
  ! This is a new implementation of the vector assembly using element subsets.
  ! In contrast to linf_buildVectorScalar, this routine loops itself about
  ! the element subsets and calls linf_initAssembly/
  ! linf_assembleSubmeshVector/linf_doneAssembly to assemble vector entries of a
  ! submesh.
  ! The linf_assembleSubmeshVector interface allows to assemble parts of a
  ! vector based on an arbitrary element list which is not bound to an
  ! element distribution.
!</description>

!<input>
  ! The linear form specifying the underlying PDE of the discretisation.
  type(t_linearForm), intent(in) :: rform
  
  ! Whether to clear the vector before calculating the entries.
  ! If .FALSE., the new entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! A callback routine for the function to be discretised.
  include 'intf_coefficientVectorSc.inc'
  optional :: fcoeff_buildVectorSc_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated entries are imposed to this vector.
  ! The vector must exist before being passed to this routine.
  type(t_vectorScalar), intent(inout) :: rvectorScalar
  
  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_linfVectorAssembly) :: rvectorAssembly
    integer :: ielementDistr
    integer, dimension(:), pointer :: p_IelementList
    
    ! If the vector does not exist, stop here.
    if (rvectorScalar%h_Ddata .eq. ST_NOHANDLE) then  
      call output_line('Vector not available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalar2')
    end if
    
    ! If the vector os stored in interleave format, stop here.
    if (rvectorScalar%NVAR .ne. 1) then
      call output_line('Vector must not be in interleaved format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalar2')
    end if
    
    ! The vector must be unsorted, otherwise we can not set up the vector.
    if (rvectorScalar%isortStrategy .gt. 0) then
      call output_line('Vector must be unsorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalar2')
      call sys_halt()
    end if
    
    ! Clear the vector if necessary.
    if (bclear) call lsyssc_clearVector (rvectorScalar)
    
    ! The vector must provide a discretisation structure
    if (.not. associated(rvectorScalar%p_rspatialDiscr)) then
      call output_line('No discretisation associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalar2')
      call sys_halt()
    end if
    
    ! Do we have a uniform triangulation? Would simplify a lot...
    if ((rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) .or.&
        (rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_CONFORMAL)) then 
      
      select case(rvectorScalar%cdataType)
        
      case(ST_DOUBLE)
        ! Loop over the element distributions.
        do ielementDistr = 1,rvectorScalar%p_rspatialDiscr%inumFESpaces
          
          ! Check if element distribution is empty
          if (rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%NEL .le. 0) cycle
          
          ! Get list of elements in distribution
          call storage_getbase_int(&
              rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%h_IelementList,&
              p_IelementList)
          
          ! Initialise a vector assembly structure for that element distribution
          call linf_initAssembly(rvectorAssembly,rform,&
              rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%celement,&
              rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%ccubTypeLinForm,&
              min(LINF_NELEMSIM, size(p_IelementList)))
          
          ! Assemble the data for all elements in this element distribution
          call linf_assembleSubmeshVector (rvectorAssembly,rvectorScalar,&
              p_IelementList,fcoeff_buildVectorSc_sim,rcollection)
          
          ! Release the assembly structure.
          call linf_doneAssembly(rvectorAssembly)
        end do
        
      case DEFAULT
        call output_line('Single precision vectors currently not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalar2')
      end select
      
    else
      call output_line('General discretisation not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalar2')
      call sys_halt()
    end if

  end subroutine  

  !****************************************************************************

!<subroutine>

  subroutine linf_buildVecIntlScalar2 (rform, bclear, rvectorScalar,&
                                       fcoeff_buildVectorBl_sim, rcollection)
  
!<description>
  ! This routine assembles the entries of a vector according to a linear form
  ! (typically used for assembling RHS vectors).
  !
  ! If bclear=TRUE, the vector is cleared before the assembly and any 
  ! sorting of the entries is switched off - the vector is set up unsorted.
  !
  ! If bclear=FALSE, the vector must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
  !
  ! IMPLEMENTATIONAL REMARK:
  ! This is a new implementation of the vector assembly using element subsets.
  ! In contrast to linf_buildVectorScalar, this routine loops itself about
  ! the element subsets and calls linf_initAssembly/
  ! linf_assembleSubmeshVector/linf_doneAssembly to assemble vector entries of a
  ! submesh.
  ! The linf_assembleSubmeshVector interface allows to assemble parts of a
  ! vector based on an arbitrary element list which is not bound to an
  ! element distribution.
!</description>

!<input>
  ! The linear form specifying the underlying PDE of the discretisation.
  type(t_linearForm), intent(in) :: rform
  
  ! Whether to clear the vector before calculating the entries.
  ! If .FALSE., the new entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! A callback routine for the function to be discretised.
  include 'intf_coefficientVectorBl.inc'
  optional :: fcoeff_buildVectorBl_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated entries are imposed to this vector.
  ! The vector must exist before being passed to this routine.
  type(t_vectorScalar), intent(inout) :: rvectorScalar
  
  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_linfVectorAssembly) :: rvectorAssembly
    integer :: ielementDistr
    integer, dimension(:), pointer :: p_IelementList

    ! If the vector does not exist, stop here.
    if (rvectorScalar%h_Ddata .eq. ST_NOHANDLE) then  
      call output_line('Vector not available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalar2')
    end if
    
    ! The vector must be unsorted, otherwise we can not set up the vector.
    if (rvectorScalar%isortStrategy .gt. 0) then
      call output_line('Vector must be unsorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalar2')
      call sys_halt()
    end if
    
    ! Clear the vector if necessary.
    if (bclear) call lsyssc_clearVector (rvectorScalar)
    
    ! The vector must provide a discretisation structure
    if (.not. associated(rvectorScalar%p_rspatialDiscr)) then
      call output_line('No discretisation associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalar2')
      call sys_halt()
    end if
    
    ! Do we have a uniform triangulation? Would simplify a lot...
    if ((rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) .or.&
        (rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_CONFORMAL)) then 
      
      select case(rvectorScalar%cdataType)
        
      case(ST_DOUBLE)
        ! Loop over the element distributions.
        do ielementDistr = 1,rvectorScalar%p_rspatialDiscr%inumFESpaces
          
          ! Check if element distribution is empty
          if (rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%NEL .le. 0) cycle
          
          ! Get list of elements in distribution
          call storage_getbase_int(&
              rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%h_IelementList,&
              p_IelementList)
          
          ! Initialise a vector assembly structure for that element distribution
          call linf_initAssembly(rvectorAssembly,rform,&
              rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%celement,&
              rvectorScalar%p_rspatialDiscr%RelementDistr(ielementDistr)%ccubTypeLinForm,&
              min(LINF_NELEMSIM, size(p_IelementList)))
          
          ! Assemble the data for all elements in this element distribution
          call linf_assembleSubmeshVecIntl (rvectorAssembly,rvectorScalar,&
              p_IelementList,fcoeff_buildVectorBl_sim,rcollection)
          
          ! Release the assembly structure.
          call linf_doneAssembly(rvectorAssembly)
        end do
        
      case DEFAULT
        call output_line('Single precision vectors currently not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalar2')
      end select
      
    else
      call output_line('General discretisation not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVecIntlScalar2')
      call sys_halt()
    end if
    
  end subroutine  

  !****************************************************************************

!<subroutine>

  subroutine linf_buildVectorBlock2 (rform, bclear, rvectorBlock,&
                                     fcoeff_buildVectorBl_sim, rcollection)
  
!<description>
  ! This routine assembles the entries of a vector according to a linear form
  ! (typically used for assembling RHS vectors).
  !
  ! If bclear=TRUE, the vector is cleared before the assembly and any 
  ! sorting of the entries is switched off - the vector is set up unsorted.
  !
  ! If bclear=FALSE, the vector must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
  !
  ! IMPLEMENTATIONAL REMARK:
  ! This is a new implementation of the vector assembly using element subsets.
  ! In contrast to linf_buildVectorBlock, this routine loops itself about
  ! the element subsets and calls linf_initAssembly/
  ! linf_assembleSubmeshVector/linf_doneAssembly to assemble vector entries of a
  ! submesh.
  ! The linf_assembleSubmeshVector interface allows to assemble parts of a
  ! vector based on an arbitrary element list which is not bound to an
  ! element distribution.
!</description>

!<input>
  ! The linear form specifying the underlying PDE of the discretisation.
  type(t_linearForm), intent(in) :: rform
  
  ! Whether to clear the vector before calculating the entries.
  ! If .FALSE., the new entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! A callback routine for the function to be discretised.
  include 'intf_coefficientVectorBl.inc'
  optional :: fcoeff_buildVectorBl_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated entries are imposed to this vector.
  ! The vector must exist before being passed to this routine.
  type(t_vectorBlock), intent(inout) :: rvectorBlock
  
  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_linfVectorAssembly) :: rvectorAssembly
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    integer :: ielementDistr, iblock
    integer, dimension(:), pointer :: p_IelementList
    logical :: bcompatible
    
    ! If the vector does not exist, stop here.
    if (rvectorBlock%h_Ddata .eq. ST_NOHANDLE) then  
      call output_line('Vector not available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlock2')
    end if
    
    ! The vector must be unsorted, otherwise we can not set up the vector.
    do iblock = 1, rvectorBlock%nblocks
      if (rvectorBlock%RvectorBlock(iblock)%isortStrategy .gt. 0) then
        call output_line('Vector must be unsorted!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlock2')
        call sys_halt()
      end if
    end do
    
    ! Clear the vector if necessary.
    if (bclear) call lsysbl_clearVector (rvectorBlock)
    
    ! The vector must provide a block discretisation structure
    if (.not. associated(rvectorBlock%p_rblockDiscr)) then
      call output_line('No block discretisation associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlock2')
      call sys_halt()
    end if
    
    ! The discretisation of all scalar subvectors must be the same
    do iblock = 2, rvectorBlock%nblocks
      call spdiscr_isDiscrCompatible(rvectorBlock%p_rblockDiscr%RspatialDiscr(1),&
          rvectorBlock%p_rblockDiscr%RspatialDiscr(iblock), bcompatible)
      if (.not.bcompatible) then
        call output_line('All scalar subvectors must have the same spatial discretisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlock2')
        call sys_halt()
      end if
    end do

    ! Set pointer to first spatial discretisation which is used as template
    p_rspatialDiscr => rvectorBlock%p_rblockDiscr%RspatialDiscr(1)

    ! Do we have a uniform triangulation? Would simplify a lot...
    if ((p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) .or.&
        (p_rspatialDiscr%ccomplexity .eq. SPDISC_CONFORMAL)) then 
    
      select case(rvectorBlock%cdataType)

      case(ST_DOUBLE)
        ! Loop over the element distributions.
        do ielementDistr = 1,p_rspatialDiscr%inumFESpaces
          
          ! Check if element distribution is empty
          if (p_rspatialDiscr%RelementDistr(ielementDistr)%NEL .le. 0) cycle
          
          ! Get list of elements in distribution
          call storage_getbase_int(&
              p_rspatialDiscr%RelementDistr(ielementDistr)%h_IelementList,&
              p_IelementList)

          ! Initialise a vector assembly structure for that element distribution
          call linf_initAssembly(rvectorAssembly,rform,&
              p_rspatialDiscr%RelementDistr(ielementDistr)%celement,&
              p_rspatialDiscr%RelementDistr(ielementDistr)%ccubTypeLinForm,&
              min(LINF_NELEMSIM, size(p_IelementList)))

          ! Assemble the data for all elements in this element distribution
          call linf_assembleSubmeshVectorBlock (rvectorAssembly,rvectorBlock,&
              p_IelementList,fcoeff_buildVectorBl_sim,rcollection)
          
          ! Release the assembly structure.
          call linf_doneAssembly(rvectorAssembly)
        end do

      case DEFAULT
        call output_line('Single precision vectors currently not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlock2')
      end select
      
    else
      call output_line('General discretisation not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlock2')
      call sys_halt()
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine linf_buildVectorBlockBdr2d (rform, ccubType, bclear, rvectorBlock,&
                                         fcoeff_buildVectorBlBdr2D_sim,&
                                         rboundaryRegion, rcollection)
  
!<description>
  ! This routine assembles the entries of a vector according to a linear form
  ! (typically used for assembling RHS vectors).
  !
  ! If bclear=TRUE, the vector is cleared before the assembly and any 
  ! sorting of the entries is switched off - the vector is set up unsorted.
  !
  ! If bclear=FALSE, the vector must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
!</description>

!<input>
  ! The linear form specifying the underlying PDE of the discretisation.
  type(t_linearForm), intent(in) :: rform

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer(I32), intent(in) :: ccubType
  
  ! Whether to clear the vector before calculating the entries.
  ! If .FALSE., the new entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion
  
  ! A callback routine for the function to be discretised.
  include 'intf_coefficientVectorBlBdr2D.inc'
  optional :: fcoeff_buildVectorBlBdr2D_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated entries are imposed to this vector.
  type(t_vectorBlock), intent(inout) :: rvectorBlock
  
  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_linfVectorAssembly) :: rvectorAssembly
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundaryRegion) :: rboundaryReg
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    real(DP), dimension(:,:), pointer :: DedgePosition
    integer, dimension(:), pointer :: IelementList, IelementOrientation
    integer :: ibdc,ielementDistr,NELbdc,iblock
    logical :: bcompatible
    
    ! If the vector does not exist, stop here.
    if (rvectorBlock%h_Ddata .eq. ST_NOHANDLE) then  
      call output_line('Vector not available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlockBdr2D')
    end if
    
    ! The vector must be unsorted, otherwise we can not set up the vector.
    do iblock = 1, rvectorBlock%nblocks
      if (rvectorBlock%RvectorBlock(iblock)%isortStrategy .gt. 0) then
        call output_line('Vector must be unsorted!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlockBdr2D')
        call sys_halt()
      end if
    end do
    
    ! Clear the vector if necessary.
    if (bclear) call lsysbl_clearVector (rvectorBlock)
    
    ! The vector must provide a block discretisation structure
    if (.not. associated(rvectorBlock%p_rblockDiscr)) then
      call output_line('No block discretisation associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlockBdr2D')
      call sys_halt()
    end if
    
    ! The discretisation of all scalar subvectors must be the same
    do iblock = 2, rvectorBlock%nblocks
      call spdiscr_isDiscrCompatible(rvectorBlock%p_rblockDiscr%RspatialDiscr(1),&
          rvectorBlock%p_rblockDiscr%RspatialDiscr(iblock), bcompatible)
      if (.not.bcompatible) then
        call output_line('All scalar subvectors must have the same spatial discretisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlockBdr2D')
        call sys_halt()
      end if
    end do

    ! Set pointer to first spatial discretisation which is used as template
    p_rspatialDiscr => rvectorBlock%p_rblockDiscr%RspatialDiscr(1)
    
    ! The discretisation must provide a triangulation structure
    if (.not. associated(p_rspatialDiscr%p_rtriangulation)) then
      call output_line('No triangulation associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlockBdr2d')
      call sys_halt()
    end if
    
    ! The discretisation must provide a boundary structure
    if (.not. associated(p_rspatialDiscr%p_rboundary)) then
      call output_line('No boundary associated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlockBdr2d')
      call sys_halt()
    end if
    
    ! Set pointers for quicker access
    p_rboundary => p_rspatialDiscr%p_rboundary
    p_rtriangulation => p_rspatialDiscr%p_rtriangulation
    
    ! Do we have a uniform triangulation? Would simplify a lot...
    if ((p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) .or.&
        (p_rspatialDiscr%ccomplexity .eq. SPDISC_CONFORMAL)) then 
      
      select case(rvectorBlock%cdataType)
        
      case(ST_DOUBLE)
        
        if (present(rboundaryRegion)) then
          
          ! Calculate total number of elements adjacent to the boundary region
          NELbdc = bdraux_getNELAtRegion(rboundaryRegion, p_rtriangulation)
          
          ! Allocate memory for element list, element orientation and
          ! the start- and end-parameter values of edges at the boundary
          allocate(IelementList(NELbdc), IelementOrientation(NELbdc))
          allocate(DedgePosition(2,NELbdc))
          
          ! Loop over the element distributions.
          do ielementDistr = 1,p_rspatialDiscr%inumFESpaces
            
            ! Calculate the list of elements adjacent to the boundary
            call bdraux_getElementsAtRegion(rboundaryRegion,&
                p_rspatialDiscr, NELbdc,&
                IelementList, IelementOrientation, DedgePosition,&
                p_rspatialDiscr%RelementDistr(ielementDistr)%celement)
            
            ! Check if element distribution is empty
            if (NELbdc .le. 0) cycle
            
            ! Initialise a vector assembly structure for that element distribution
            call linf_initAssembly(rvectorAssembly, rform,&
                p_rspatialDiscr%RelementDistr(ielementDistr)%celement,&
                ccubType, min(LINF_NELEMSIM, NELbdc))
            
            ! Assemble the data for all elements in this element distribution
            call linf_assembleSubmeshVectorBlockBdr2D (rvectorAssembly, rvectorBlock,&
                rboundaryRegion, IelementList(1:NELbdc), IelementOrientation(1:NELbdc),&
                DedgePosition(:,1:NELbdc), fcoeff_buildVectorBlBdr2D_sim, rcollection)
            
            ! Release the assembly structure.
            call linf_doneAssembly(rvectorAssembly)
            
          end do
          
          ! Release memory
          deallocate(IelementList, IelementOrientation, DedgePosition)
          
        else
          
          ! Loop over the element distributions.
          do ielementDistr = 1,p_rspatialDiscr%inumFESpaces
            
            ! Check if element distribution is empty
            if (p_rspatialDiscr%RelementDistr(ielementDistr)%NEL .le. 0) cycle
            
            ! Initialise a vector assembly structure for that element distribution
            call linf_initAssembly(rvectorAssembly, rform,&
                p_rspatialDiscr%RelementDistr(ielementDistr)%celement,&
                ccubType, LINF_NELEMSIM)
            
            ! Create a boundary region for each boundary component and call
            ! the calculation routine for that.
            do ibdc = 1,boundary_igetNBoundComp(p_rboundary)
              call boundary_createRegion (p_rboundary, ibdc, 0, rboundaryReg)
              
              ! Calculate number of elements adjacent to the boundary region
              NELbdc = bdraux_getNELAtRegion(rboundaryReg, p_rtriangulation)
              
              ! Check if element distribution is empty
              if (NELbdc .le. 0) cycle
              
              ! Allocate memory for element list, element orientation and
              ! the start- and end-parameter values of edges at the boundary
              allocate(IelementList(NELbdc), IelementOrientation(NELbdc))
              allocate(DedgePosition(2,NELbdc))
              
              ! Calculate the list of elements adjacent to the boundary
              call bdraux_getElementsAtRegion(rboundaryReg,&
                  p_rspatialDiscr, NELbdc,&
                  IelementList, IelementOrientation, DedgePosition,&
                  p_rspatialDiscr%RelementDistr(ielementDistr)%celement)
              
              if (NELbdc .gt. 0) then
                
                ! Assemble the data for all elements in this element distribution
                call linf_assembleSubmeshVectorBlockBdr2D (rvectorAssembly, rvectorBlock,&
                    rboundaryReg, IelementList(1:NELbdc), IelementOrientation(1:NELbdc),&
                    DedgePosition(:,1:NELbdc), fcoeff_buildVectorBlBdr2D_sim, rcollection)
                
              end if
              
              ! Deallocate memory
              deallocate(IelementList, IelementOrientation, DedgePosition)
              
            end do ! ibdc
            
            ! Release the assembly structure.
            call linf_doneAssembly(rvectorAssembly)
            
          end do ! ielementDistr
          
        end if
        
      case DEFAULT
        call output_line('Single precision vectors currently not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlockBdr2D')
      end select
      
    else
      call output_line('General discretisation not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorBlockBdr2D')
      call sys_halt()
    end if
 
  end subroutine
  
end module linearformevaluation
