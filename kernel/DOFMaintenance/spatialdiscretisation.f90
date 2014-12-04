!##############################################################################
!# ****************************************************************************
!# <name> spatialdiscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic discretisation structures, collecting
!# information about which geometric elements are discretised by which types
!# of finite elements.
!#
!# The following routines can be found in this module:
!#
!# 1.) spdiscr_initBlockDiscr
!#     -> Initialises a block discretisation structure for the
!#        discretisation of multiple equations
!#
!# 2.) spdiscr_initDiscr_simple
!#     -> Initialise a scalar discretisation structure. One element type for
!#        all geometric elements, for test and trial functions.
!#
!# 3.) spdiscr_initDiscr_triquad
!#     -> Initialise a scalar discretisation structure for a mixed
!#        triangle/quad discretisation.
!#
!# 4.) spdiscr_deriveBlockDiscr
!#     -> Creates a block discretisation structure as a subset of
!#        another block discretisation structure
!#
!# 5.) spdiscr_deriveSimpleDiscrSc
!#     -> Based on an existing discretisation structure, derive a new
!#        discretisation structure with different trial elements
!#
!# 6.) spdiscr_deriveDiscr_triquad
!#     -> Based on an existing discretisation structure, derive a new
!#        discretisation structure for a mixed triangular/quad mesh.
!#
!# 7.) spdiscr_releaseDiscr
!#     -> Release a scalar discretisation structure.
!#
!# 8.) spdiscr_createBlockDiscrInd
!#     -> Creates a block discretisation with one component from an
!#        existing scalar discretisation.
!#
!# 9.) spdiscr_releaseBlockDiscr
!#     -> Releases a block discretisation structure from memory
!#
!# 10.) spdiscr_checkCubature
!#      -> Checks if a cubature formula is compatible to an element
!#         distribution.
!#
!# 11.) spdiscr_duplicateDiscrSc
!#      -> Copies a spatial discretisation structure to another
!#
!# 12.) spdiscr_duplicateBlockDiscr
!#      -> Copies a block discretisation structure to another
!#
!# 13.) spdiscr_getLumpCubature
!#      -> Try to get a cubature formula for an element type that leads to
!#         diagonal lumping when setting up a mass matrix with that element
!#
!# 14.) spdiscr_getStdCubature
!#      -> Try to get the typical cubature formula for an element
!#
!# 15.) spdiscr_infoBlockDiscr
!#      -> Outputs information about the block discretisation
!#         (mostly for debugging)
!#
!# 16.) spdiscr_infoDiscr
!#      -> Outputs information about the spatial discretisation
!#         (mostly for debugging)
!#
!# 17.) spdiscr_infoElementDistr
!#      -> Outputs information about the element group
!#         (mostly for debugging)
!#
!# 18.) spdiscr_igetNDofLocMax
!#      -> Calculate the maximum number of local DOF`s in a discretisation
!#
!# 19.) spdiscr_checkCubature
!#      -> Checks if the cubature formula of type icubType can be applied
!#         to the elements of the type celement
!#
!# 20.) spdiscr_releaseDofMapping
!#      -> Release precomputed DOF-mapping data
!#
!# 21.) spdscr_edgeBlocking2D
!#      -> Blocks the edges in a triangulation according to the FE spaces
!#         of the adjacent elements
!#
!# 22.) spdscr_releaseEdgeBlocking
!#      -> Releases the edge blocking structure allocated by
!#         spdscr_edgeBlocking2D
!#
!# 23.) spdiscr_concatBlockDiscr
!#     -> Concatenates two block discretisation structures to a new
!#        block discretisation structure.
!#
!# 24.) spdisc_isBlockDiscrCompatible
!#      -> Checks whether two block discretisation structures are compatible
!#
!# 25.) spdiscr_isDiscrCompatible
!#      -> Checks whether two spatial discretisation structures are compatible
!#
!# 26.) spdiscr_isElemDistrCompatible
!#      -> Checks whether two element group structures are compatible
!#
!# 25.) spdiscr_createDefCubStructure
!#      -> Create a default cubature information structure based on a
!#         discretisation structure.
!#
!# 26.) spdiscr_releaseCubStructure
!#      -> Release an cubature information structure.
!#
!# 27.) spdiscr_getStdDiscrInfo
!#      -> Obtain typical information from discretisation/assembly information
!#         structure which is used during the assembly.
!#
!# 28.) spdiscr_defineCubature
!#      -> Initialise a cubature information structure based on an automatic
!#         cubature formula.
!#
!# 29.) spdiscr_getElementCubMapping
!#      -> Calculate for every element its cubature formula
!#
!# 30.) spdiscr_initBlockDiscr_free
!#      -> Initialise a "free" block discretisation.
!#
!# 31.) spdiscr_initDiscr_free
!#      -> Initialise a "free" scalar discretisation
!#
!# 32.) spdiscr_createDefTrafoStructure
!#      -> Creates a default transformation structure based on a discretisation
!#
!# 33.) spdiscr_createDefTrafoByCubInfo
!#      -> Creates a default transformation structure based on a 
!#         cubature structure
!#
!# 34.) spdiscr_defaultTrafo
!#      -> Auxiliary routine: Default initialisation of a transformation structure
!#         based on a cubature or discretisation structure 
!#
!# 35.) spdiscr_releaseTrafoStructure
!#      -> Releases a transformation structure
!#
!# 36.) spdiscr_getElementTrafoMapping
!#      -> Obtains for every element the usesd transformation
!#
!# 37.) spdiscr_getStdTrafoInfo
!#      -> Returns basic information of a transformation structure
!#
!# 38.) spdiscr_infoCubatureInfo
!#      -> Outputs information about the cubature info structure
!#         (mostly for debugging)
!#
!# 39.) spdiscr_infoCubatureInfoBlock
!#      -> Outputs information about the block of a cubature info structure
!#         (mostly for debugging)
!#
!# 40.) spdiscr_infoTrafoInfo
!#      -> Outputs information about the trafo info structure
!#         (mostly for debugging)
!#
!# 41.) spdiscr_infoTrafoInfoBlock
!#      -> Outputs information about the block of a trafo info structure
!#         (mostly for debugging)
!#
!# 42.) spdiscr_appendBlockComponent
!#      -> Add a new spatial discretisation as component to a block discretisation
!#
!# 43.) spdiscr_commitBlockDiscr
!#      -> Commits a block discretisation
!#
!# 44.) spdiscr_getNelemGroups
!#      -> Determins the number of element groups. Each element group has an
!#         element type and an underlying transformation.
!#
!# 45.) spdiscr_getElemGroupInfo
!#      -> Determins information about an element group.
!# 
!#   The cubature information structure \\
!# -------------------------------------- \\
!#
!# A t_scalarCubatureInfo structure is associated to a discretisation
!# structure. For the assembly of a matrix or vector with a specific cubature
!# formula e.g., the code can create a default assembly structure, modify the
!# cubature formula and assemble with this:
!#
!# <code>
!#   ! Use 4-pt Gauss formula. We assume for this example: 2D QUAD mesh
!#
!#   ! Get a structure and modify the cubature formulas.
!#   call spdiscr_createDefCubStructure (rdiscretisation,rcubatureInfo,CUB_GEN_AUTO_G4)
!#
!#   ! Assemble a matrix based on this.
!#   call bilf_buildMatrixScalar2 (rbilform,.true.,rmatrix,&
!#       rcubatureInfo=rcubatureInfo)
!#
!#   ! Assemble a RHS vector based on this.
!#   call linf_buildVectorScalar (rlinform,.true.,rvector,fcallback,&
!#       rcubatureInfo=rcubatureInfo)
!#
!#   ! Release the info structure.
!#   call spdiscr_releaseCubStructure (rcubatureInfo)
!# </code>
!#
!# In this example, a default rcubatureInfo structure is created based on the
!# cubature formula CUB_GEN_AUTO_G4. This constant represents the automatic
!# 4-pt Gauss formula which automatically adapts to the underlying finite
!# element: It automatically chooses the proper 4-pt Gauss formula, independent
!# of whether the finite element is 2D, 3D, triangular, quadrilateral or whatever.
!# 
!# If the user wants to choose a specific cubature formula, the corresponding
!# constant can directly be written into the structure as follows:
!#
!# <code>
!#   ! Use 4-pt Gauss formula. We assume for this example: 2D QUAD mesh
!#
!#   ! Get a structure and modify the cubature formulas.
!#   call spdiscr_createDefCubStructure (rdiscretisation,rcubatureInfo)
!#   rcubatureInfo%p_RinfoBlocks(:)%ccubature = CUB_G4_2D
!#
!#   ! Assemble a matrix based on this.
!#   call bilf_buildMatrixScalar2 (rbilform,.true.,rmatrix,&
!#       rcubatureInfo=rcubatureInfo)
!#
!#   ! Assemble a RHS vector based on this.
!#   call linf_buildVectorScalar (rlinform,.true.,rvector,fcallback,&
!#       rcubatureInfo=rcubatureInfo)
!#
!#   ! Release the info structure.
!#   call spdiscr_releaseCubStructure (rcubatureInfo)
!# </code>
!#
!# However it has to be noted that in this example, it is not possible to
!# write "...%ccubature = CUB_GEN_AUTO_G4". The automatic cubature formula
!# only works in combination with "spdiscr_createDefCubStructure".
!# </purpose>
!##############################################################################

module spatialdiscretisation

!$ use omp_lib
  use fsystem
  use storage
  use triangulation
  use basicgeometry
  use boundary
  use transformation
  use element
  use cubature
  use genoutput
  use sort

  implicit none

  private

!<constants>

!<constantblock description="Constants defining the complexity of the discretisation">

  ! Free discretisation. No FEM space associated. No elements associated.
  ! Just a set of DOFs.
  integer, parameter, public :: SPDISC_FREE      = -1

  ! Uniform discretisation: All elements are of the same type.
  integer, parameter, public :: SPDISC_UNIFORM   = 0

  ! Conformal discretisation: Elements of different FE spaces are mixed,
  ! but the DOF`s 'fit together' (e.g. quads/tri elements with same DOF`s,
  ! isoparametric elements on the boundary).
  integer, parameter, public :: SPDISC_CONFORMAL = 1

  ! Mixed discretisation: Elements of different FE spaces, arbitrary mixed.
  integer, parameter, public :: SPDISC_MIXED     = 2

!</constantblock>

!<constantblock description="Additional constants for cubature formulas">

  ! Automatically determine cubature formula for a discretisation.
  integer(I32), parameter, public :: SPDISC_CUB_AUTOMATIC = CUB_GEN_AUTO

  ! Cubature formula stays unchanged.
  integer(I32), parameter, public :: SPDISC_CUB_NOCHANGE  = -20_I32

!</constantblock>

!<constantblock description="Operator types">

  ! Mass matrix
  integer(I32), parameter, public :: SPDISC_OPTP_MASS = 0_I32

  ! Laplace matrix
  integer(I32), parameter, public :: SPDISC_OPTP_LAPLACE = 1_I32

  ! RHS
  integer(I32), parameter, public :: SPDISC_OPTP_RHS = 2_I32

  ! Convection matrix
  integer(I32), parameter, public :: SPDISC_OPTP_CONVEC = 3_I32

!</constantblock>

!<constantblock variable="ccubType" description="Additional generic cubature formulas">

  ! DEPRECATED: Use ccubTypeBilForm from the discretisation structure.
  integer(I32), parameter, public :: CUB_GEN_DEPR_BILFORM = CUB_TP_DEPR + 10_I32

  ! DEPRECATED: Use ccubTypeLinForm from the discretisation structure.
  integer(I32), parameter, public :: CUB_GEN_DEPR_LINFORM = CUB_TP_DEPR + 11_I32

  ! DEPRECATED: Use ccubTypeEval from the discretisation structure.
  integer(I32), parameter, public :: CUB_GEN_DEPR_EVAL = CUB_TP_DEPR + 12_I32

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! element group structure. This structure collects for one type
  ! of element (e.g. <tex>$Q_1$</tex>), on which geometric element primitives it is
  ! to be used. In the t_spatialDiscretisation there is a list of these
  ! element structures for each type of element. This e.g. allows to use
  ! triangular elements with <tex>$P_1$</tex>, quad elements with <tex>$Q_1$</tex> and possible
  ! isoparametric "surface" elements to be mixed in the triangulation.
  !
  ! The structure is assigned to a triangulation by means of the 'parent'
  ! structure t_spatialDiscretisation, which contains a pointer to it.

  type t_elementDistribution

    ! Element identifier for Finite Element functions to use in this
    ! element list.
    integer(I32) :: celement        = EL_UNDEFINED

    ! DEPRECATED: Cubature formula to use for the discretisation of this element pair
    ! during the evaluation of bilinear forms (matrix generation).
    ! Note: When evaluating bilinear forms, the ccubTypeBilForm
    ! constant of the test space decides about the cubature formula
    ! to be used!
    integer(I32) :: ccubTypeBilForm      = CUB_UNDEFINED

    ! DEPRECATED: Cubature formula to use for the discretisation of this element pair
    ! during the evaluation of linear forms (RHS generation).
    integer(I32) :: ccubTypeLinForm      = CUB_UNDEFINED

    ! DEPRECATED: Cubature formula to use for the evaluation of integrals over an FE
    ! function. This is used e.g. in postprocessing routines to calculate
    ! an integral to get an error to a reference solution.
    integer(I32) :: ccubTypeEval         = CUB_UNDEFINED

    ! Type of transformation to use from the reference element to
    ! the real element. One of the TRAFO_IDxxxx constants of the module
    ! 'transformation' identifying the type of transformation.
    ! The same transformation is used for both, the trial and the test
    ! space, during the evaluation of linear as well as bilinear forms
    ! (matrix/RHS generation).
    integer(I32) :: ctrafoType      = TRAFO_ID_UNKNOWN

    ! Number of elements in the list p_IelementList.
    ! May vary from the actual length of p_IelementList!
    integer :: NEL = 0

    ! Handle to list of element numbers that are discretised with this
    ! combination of trial/test functions.
    ! If NEL=0, the element list is empty, i.e. h_IelementList = ST_NOHANDLE!
    integer :: h_IelementList       = ST_NOHANDLE

  end type

  public :: t_elementDistribution

!</typeblock>

!<typeblock>

  ! The central discretisation structure corresponding to one mesh level.
  ! Here, all information about the discretisation are collected (mesh
  ! information, trial functions, test functions,...).
  !
  ! Remark: The structure realises only the discretisation of `a scalar
  !  equation`. For multidimensional problems, there are multiple of
  !  these structures, each one for one PDE. I this case, the structure
  !  is part of the block discretisation structure below and
  !  'hung into' each scalar matrix that discretises that equation.

  type t_spatialDiscretisation

    ! Dimension of the discretisation. 0=not initialised,
    ! 1=1D discretisation, 2=2D discretisation, 3=3D discretisation
    ! -1=arbitrary (used for free discretisations)
    integer                          :: ndimension             = 0

    ! Whether the discretisation structure is a copy of another discretisation
    ! structure. If set to TRUE, the structure was derived from another one
    ! and shares the same dynamic information (element lists,...).
    ! (This prevents the release-routine from releasing memory when
    ! cleaning up the structure.)
    logical                          :: bisCopy                = .false.

    ! Pointer to the domain that is discretised
    type(t_boundary), pointer        :: p_rboundary            => null()

    ! Pointer to the underlying triangulation of the mesh (2D)
    type(t_triangulation), pointer   :: p_rtriangulation       => null()

    ! Complexity of the discretisation. One of the SPDISC_xxxx constants
    integer                          :: ccomplexity            = SPDISC_UNIFORM

    ! Handle to the element group identifier list.
    ! For every geometric element i, IelemGroupIDs(i) specifies the
    ! number of the element group that contains that element.
    ! In a uniform discretisation (ccomplexity=SPDISC_UNIFORM), this
    ! handle is ST_NOHANDLE as all elements are in the
    ! element group 1.
    integer                          :: h_IelemGroupIDs       = ST_NOHANDLE

    ! Handle to an 'element counter' array. For every element of every
    ! type, there is a unique running number given to that element in the
    ! corresponding element subset.
    !
    ! Example: Think of a mixed mesh of triangles and quads, discretised
    !  with <tex>$P_1$</tex> and <tex>$Q_1$</tex>. Then there are two disjunct element sets,
    !  one with triangles, one with quads. Every triangle gets a running
    !  number (1,2,3,...) and every quad gets a running number (1,2,3,...).
    !  These numbers are stored here, corresponding to each element.
    !
    ! Note that the handle may be =ST_NOHANDLE. Whether it is set up or not
    ! depends on the discretisation. The information is usually used in
    ! the DOFMapping-routines to compute the mapping between local and
    ! global degrees of freedom.
    integer                          :: h_IelementCounter      = ST_NOHANDLE

    ! Number of different FE spaces mixed in this discretisation.
    ! This is the number of elements occupied in RelementDistibution.
    integer                          :: inumFESpaces           = 0

    ! List of element group structures for every element type
    ! that is used in the discretisation.
    type(t_elementDistribution), dimension(:), pointer :: RelementDistr => null()

    ! Specifies whether the DOF-mapping is precomputed.
    logical                          :: bprecompiledDofMapping = .false.

    ! Number of DOF`s total. Only available if bprecompiledDofMapping=true
    ! or ccomplexity = SPDISC_FREE. A "free" discretisation just contains
    ! ndof degrees of freedom, not associated with any element or whatever.
    integer                          :: ndof = 0

    ! Specifies for every element a list of DOF`s how to map a local DOF
    ! to a global DOF. p_IelementDofIdx is a list with starting indices
    ! for every element in this list.
    ! Only available if bprecompiledDofMapping=true.
    integer :: h_IelementDofs = ST_NOHANDLE

    ! List of starting indices. p_IelementDofIdx(iel) specifies the index
    ! in p_IelementDofs where the global DOF`s of element iel start.
    ! DIMENSION(nelements+1).
    ! Only available if bprecompiledDofMapping=true.
    integer :: h_IelementDofIdx = ST_NOHANDLE

  end type

  public :: t_spatialDiscretisation

!</typeblock>

!<typeblock>

  ! The block discretisation realises the discretsation of the actual PDE,
  ! where one large soution vector consists of one or multiple solution
  ! components (e.g. $(u_x,u_y,p)^T$) on one mesh level. There is a
  ! pointer to the underlying domain, the underlying triangulation
  ! and a list of spatial discretisation structures, each responsible
  ! for one solution component.
  !
  ! Additionally, the block discretisation structure contains information
  ! about the boundary conditions, as boundary conditions affect always
  ! the full PDE system.

  type t_blockDiscretisation

    ! Dimension of the discretisation. 0=not initialised,
    ! 1=1D discretisation, 2=2D discretisation, 3=3D discretisation
    integer                          :: ndimension             = 0

    ! Complexity of the discretisation. One of the SPDISC_xxxx constants.
    ! SPDISC_UNIFORM = all elements in each discretisation
    !   substructure RspatialDiscr(:) are the same.
    ! SPDISC_CONFORMAL = Elements of different FE spaces are mixed,
    !   but the DOF`s 'fit together'. Each discretisation substructure
    !   RspatialDiscr(:) has exactly the same number of element
    !   distributions, and each element group
    !     RspatialDiscr(1)%Relementistributions(i),
    !     RspatialDiscr(2)%Relementistributions(i),
    !     RspatialDiscr(3)%Relementistributions(i),...
    !   describe exactly the same set of elements (Same size, same type,
    !   same order in the element lists,...).
    integer                          :: ccomplexity            = SPDISC_UNIFORM

    ! Pointer to the domain that is discretised
    type(t_boundary), pointer        :: p_rboundary            => null()

    ! Pointer to the underlying triangulation of the mesh (2D)
    type(t_triangulation), pointer   :: p_rtriangulation     => null()

    ! Number of solution components maintained by this structure.
    integer                          :: ncomponents = 0

    ! A list of up to ncomponents scalar spatial discretisation structures.
    ! Each structure corresponds to one solution component and defines
    ! trial-/test-functions, complexity of the discretisation etc.
    type(t_spatialDiscretisation), dimension(:), pointer :: RspatialDiscr => null()

  end type

  public :: t_blockDiscretisation

!</typeblock>

!<typeblock>

  ! Defines the blocking of a set of edges.
  type t_edgeBlocking

    ! Number of edges described by this structure.
    integer :: nedges

    ! Number of FE spaces
    integer :: nelementDistributions

    ! Number of combinations of FE spaces
    integer :: nfecombinations

    ! A list of all edges, ordered in blocks in such a way that all edges in a block
    ! have the same FE spaces at the adjacent elements.
    integer :: h_Iedges

    ! Array of length (#combinations+1). Defines the start positions in Iedges
    ! of every block of edges that have the same type of elements adjacent.
    ! #combinations = #element groups * (#element groups+1)/2
    integer :: h_IdistPositions

    ! Array of dimension(2,#combinations/2).
    ! Returns for each block described by IdistPositions the number of the
    ! element groups of the elements on each side of the edge.
    ! #combinations = #element groups * (#element groups+1)/2
    integer :: h_Idistributions
  end type

  public :: t_edgeBlocking

!</typeblock>

!<typeblock>

  ! Contains information that configures the transformation from the reference
  ! to the real element(s) in the assembly of matrices and vectors 
  ! for a set of elements.
  type t_scalarTrafoInfoBlock

    ! Transformation to use.
    integer(I32) :: ctrafoType = TRAFO_ID_UNKNOWN

    ! Id of the element set, this assembly block refers to.
    integer :: ielemGroup = 0
    
    ! Number of elements in this block.
    ! =-1: structure not initialised.
    integer :: NEL = -1

    ! If h_IelementList != ST_NOHANDLE, this specifies a handle to a list of elements 
    ! where to apply the above cubature rule. All elements shall be in the element set
    ! ielemGroup!
    ! If h_IelementList=ST_NOHANDLE, the element list for the above cubature rule 
    ! is the complete list of elements in the element group ielemGroup.
    integer :: h_IelementList = ST_NOHANDLE

    ! Ownership flag. This flag is set to .true. by the routines in this module
    ! if h_IelementList is internally created. Assures that
    ! spdiscr_releaseCubStructure does not accidentally release memory
    ! that was not allocated here.
    logical :: blocalElementList = .false.

  end type

!</typeblock>

!<typeblock>

  ! Contains information that configures the transformation from the
  ! reference to the real element(s) in the assembly
  ! of matrices and vectors.
  type t_scalarTrafoInfo

    ! Number of cubature information blocks in p_RinfoBlocks.
    integer :: ninfoBlockCount = 0

    ! A list of cubature information blocks.
    type(t_scalarTrafoInfoBlock), dimension(:), pointer :: p_RinfoBlocks => null()

  end type

!</typeblock>

!<typeblock>

  ! Contains information that configures the cubature in the assembly 
  ! of matrices and vectors for a set of elements.
  type t_scalarCubatureInfoBlock

    ! Cubature rule to use.
    integer(I32) :: ccubature = CUB_UNDEFINED

    ! Id of the element group, this assembly block refers to.
    integer :: ielemGroup = 0
    
    ! Id of the transformation block, this assembly block refers to.
    integer :: itrafoBlock = 0

    ! Number of elements in this block.
    ! =-1: structure not initialised.
    integer :: NEL = -1

    ! If h_IelementList != ST_NOHANDLE, this specifies a handle to a list of elements 
    ! where to apply the above cubature rule. All elements shall be in the element set
    ! ielemGroup!
    ! If h_IelementList=ST_NOHANDLE, the element list for the above cubature rule 
    ! is the complete list of elements in the element group ielemGroup.
    integer :: h_IelementList = ST_NOHANDLE

    ! Ownership flag. This flag is set to .true. by the routines in this module
    ! if h_IelementList is internally created. Assures that
    ! spdiscr_releaseCubStructure does not accidentally release memory
    ! that was not allocated here.
    logical :: blocalElementList = .false.

  end type

!</typeblock>

!<typeblock>

  ! Contains information that configures the cubature in the assembly
  ! of matrices and vectors.
  type t_scalarCubatureInfo

    ! Number of cubature information blocks in p_RinfoBlocks.
    integer :: ninfoBlockCount = 0

    ! A list of cubature information blocks.
    type(t_scalarCubatureInfoBlock), dimension(:), pointer :: p_RinfoBlocks => null()

  end type

!</typeblock>

!</types>

  public :: spdiscr_initBlockDiscr
  public :: spdiscr_initBlockDiscr_free
  public :: spdiscr_initDiscr_simple
  public :: spdiscr_initDiscr_free
  public :: spdiscr_initDiscr_triquad
  public :: spdiscr_deriveBlockDiscr
  public :: spdiscr_deriveSimpleDiscrSc
  public :: spdiscr_deriveDiscr_triquad
  public :: spdiscr_releaseDiscr
  public :: spdiscr_createBlockDiscrInd
  public :: spdiscr_releaseBlockDiscr
  public :: spdiscr_checkCubature
  public :: spdiscr_duplicateDiscrSc
  public :: spdiscr_duplicateBlockDiscr
  public :: spdiscr_getLumpCubature
  public :: spdiscr_getStdCubature
  public :: spdiscr_infoBlockDiscr
  public :: spdiscr_infoDiscr
  public :: spdiscr_infoElementDistr
  public :: spdiscr_igetNDofLocMax
  public :: spdiscr_releaseDofMapping
  public :: spdscr_edgeBlocking2D
  public :: spdscr_releaseEdgeBlocking
  public :: spdiscr_concatBlockDiscr
  public :: spdiscr_isBlockDiscrCompatible
  public :: spdiscr_isDiscrCompatible
  public :: spdiscr_isElemDistrCompatible
  
  public :: spdiscr_getNelemGroups
  public :: spdiscr_getElemGroupInfo
  public :: spdiscr_getElemGroupIDs

  public :: t_scalarCubatureInfoBlock
  public :: t_scalarCubatureInfo
  public :: spdiscr_createDefCubStructure
  public :: spdiscr_releaseCubStructure
  public :: spdiscr_defineCubature
  public :: spdiscr_getElementCubMapping
  public :: spdiscr_infoCubatureInfo
  public :: spdiscr_infoCubatureInfoBlock

  public :: spdiscr_getStdDiscrInfo
  
  public :: t_scalarTrafoInfoBlock
  public :: t_scalarTrafoInfo
  public :: spdiscr_createDefTrafoStructure
  public :: spdiscr_createDefTrafoByCubInfo
  public :: spdiscr_defaultTrafo
  public :: spdiscr_releaseTrafoStructure
  public :: spdiscr_getElementTrafoMapping
  public :: spdiscr_getStdTrafoInfo
  public :: spdiscr_infoTrafoInfo
  public :: spdiscr_infoTrafoInfoBlock
  
  public :: spdiscr_appendBlockComponent
  public :: spdiscr_commitBlockDiscr
  
  interface spdiscr_initDiscr_simple
    module procedure spdiscr_initDiscr_simple_old
    module procedure spdiscr_initDiscr_simple_new
  end interface

  interface spdiscr_initDiscr_triquad
    module procedure spdiscr_initDiscr_triquad_old
    module procedure spdiscr_initDiscr_triquad_new
  end interface

  interface spdiscr_deriveSimpleDiscrSc
    module procedure spdiscr_deriveSimpleDiscrSc_old
    module procedure spdiscr_deriveSimpleDiscrSc_new
  end interface

  interface spdiscr_deriveDiscr_triquad
    module procedure spdiscr_deriveDiscr_triquad_old
    module procedure spdiscr_deriveDiscr_triquad_new
  end interface
  
  interface spdiscr_initBlockDiscr
    module procedure spdiscr_initBlockDiscr_fix
    module procedure spdiscr_initBlockDiscr_open
  end interface

contains

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_checkCubature (ccubType,celement)

!<description>

  ! This routine checks if the cubature formula of type icubType can be applied
  ! to the elements of the type celement.
  ! If this is not possible, an error is thrown.

!</description>

!<input>
  ! The cubature formula to be tested
  integer(I32), intent(in)                       :: ccubType

  ! The element type the cubature formula should be checked against
  integer(I32), intent(in)                       :: celement
!</input>

!</subroutine>

  logical :: bcompatible
  integer(I32) :: ishapeEL, ishapeCUB

    ! Get the shape identifiers for both the element and the cubature formula
    ishapeEL = elem_igetShape(celement)
    ishapeCUB = cub_igetShape(ccubType)

    ! The element and cubature formula are compatible if both shapes are equal,
    ! and of course the shape must be valid.
    bcompatible = (ishapeEL .eq. ishapeCUB) .and. &
                  (ishapeEL .ne. BGEOM_SHAPE_UNKNOWN)


! 'Old' Implementation follows

!  integer :: NVE, idim
!  logical :: bcompatible
!
!  ! Get from the element group the trial space and from that
!  ! the number of vertices, the element expects.
!  NVE = elem_igetNVE(celement)
!  idim = elem_igetDimension(celement)
!
!  bcompatible = .true.
!
!  ! Now we directly access the cubature constants in cubature.f90!
!  ! This is the only point in the kernel where this is necessary.
!
!  ! 1D: Line?
!  if (ccubType .le. 99) then
!    if ((NVE .ne. 2) .or. (idim .ne. NDIM1D)) bcompatible = .false.
!  end if
!
!  ! 2D: Quad?
!  if ((ccubType .ge. 200) .and. (ccubType .le. 249)) then
!    ! Tri?
!    if ((NVE .ne. 4) .or. (idim .ne. NDIM2D)) bcompatible = .false.
!  end if
!
!  ! 2D: Tri?
!  if ((ccubType .ge. 250) .and. (ccubType .le. 299)) then
!    ! Quad?
!    if ((NVE .ne. 3) .or. (idim .ne. NDIM2D)) bcompatible = .false.
!  end if
!
!  ! 3D: Hexa?
!  if ((ccubType .ge. 300) .and. (ccubType .le. 349)) then
!    if ((NVE .ne. 8) .or. (idim .ne. NDIM3D)) bcompatible = .false.
!  end if
!
!  ! 3D: Tetra?
!  if ((ccubType .ge. 350) .and. (ccubType .le. 499)) then
!    if ((NVE .ne. 4) .or. (idim .ne. NDIM3D)) bcompatible = .false.
!  end if

  ! Q2T with bubble does not work with G1X1, Trapezoidal rule and
  ! G2X2 -- Laplace matrices would be indefinite because of the definition
  ! of the bubble function!
  if (celement .eq. EL_Q2TB) then
    if ((ccubType .ge. 201) .and. (ccubType .le. 204)) bcompatible = .false.
  end if

  if (.not. bcompatible) then
    call output_line ('Element and cubature formula not compatible!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_checkCubature')
    call sys_halt()
  end if

  end subroutine

  ! ***************************************************************************

!<function>

  integer(I32) function spdiscr_getLumpCubature (celement) result (ccubType)

!<description>
  ! This routine tries to determine a cubature formula identifier according
  ! to a given element type, such that the corresponding mass matrix will
  ! get diagonal (mass lumping). If this is not possible, 0 is returned.
!</description>

!<input>
  ! An element type identifier
  integer(I32), intent(in)                       :: celement
!</input>

!<result>
  ! A cubature formula identifier that will diagonalise the mass matrix,
  ! or 0 if such an identifier is unknown / not possible.
!</result>

!</function>

    select case (elem_igetDimension(celement))
    case (NDIM1D)

      select case (elem_getPrimaryElement(celement))
      case (EL_P0_1D)
        ! Use Gauss-1
        ccubType = CUB_G1_1D

      case (EL_P1_1D)
        ! Use trapezoidal rule
        ccubType = CUB_TRZ_1D

      case (EL_P2_1D)
        ! Use Gauss-2
        ccubType = CUB_G2_1D

      case (EL_S31_1D)
        ! Use Gauss-4
        ccubType = CUB_G4_1D

      case default
        ccubType = 0
      end select

    case (NDIM2D)

      select case (elem_getPrimaryElement(celement))
      case (EL_P0)
        ! Use Gauss 1X1
        ccubType = CUB_G1_T

      case (EL_P1)
        ! Use trapezoidal rule
        ccubType = CUB_TRZ_T

      case (EL_P1T)
        ! Use Gauss-3pt
        ccubType = CUB_Collatz

      case (EL_Q0)
        ! Use Gauss 1X1
        ccubType = CUB_G1X1

      case (EL_Q1)
        ! Use trapezoidal rule
        ccubType = CUB_TRZ

      case (EL_Q1T)
        ! Use midpoint rule
        ccubType = CUB_MID

      case (EL_Q2)
        ! Summed trapezoidal rule
        ccubType = CUB_PTRZ

      case default
        ccubType = 0
      end select

    case (NDIM3D)

      select case (elem_getPrimaryElement(celement))
      case (EL_P0_3D)
        ! Use Gauss 1X1
        ccubType = CUB_G1_3D_T

      case (EL_P1_3D)
        ! Use trapezoidal rule
        ccubType = CUB_TRZ_3D_T

!      Not implemented in 3D
!      case (EL_P1T)
!        ! Use Gauss-3pt
!        ccubType = CUB_G3_T

      case (EL_Q0_3D)
        ! Use Gauss 1X1
        ccubType = CUB_G1_3D

      case (EL_Q1_3D)
        ! Use trapezoidal rule
        ccubType = CUB_TRZ_3D

      case (EL_Q1T_3D)
        ! Use midpoint rule
        ccubType = CUB_MIDAREA_3D

      case default
        ccubType = 0
      end select

    case default
      ccubType = 0
    end select

  end function

  ! ***************************************************************************

!<function>

  elemental integer(I32) function spdiscr_getStdCubature (celement,iopertype) result (ccubType)

!<description>
  ! This routine returns a standard cubature formula for an element which
  ! can be used as default when setting up matrices/vectors.
  ! If this is not possible, 0 is returned.
!</description>

!<input>
  ! An element type identifier
  integer(I32), intent(in) :: celement

  ! OPTIONAL: Type of operator which should be assembled using this cubature
  ! formula. One of the SPDISC_OPTP_xxxx constants. If not specified,
  ! SPDISC_OPTP_MASS is the default.
  integer, intent(in), optional :: iopertype
!</input>

!<result>
  ! A standard cubature formula for the assembly of matrices/vectors
  ! with the specified element celement.
!</result>

!</function>

    integer :: ioperation

    ioperation = SPDISC_OPTP_MASS
    if (present(iopertype)) ioperation = iopertype

    if ((celement .eq. EL_QPW4P0_2D) .or.&
        (celement .eq. EL_QPW4P1_2D) .or.&
        (celement .eq. EL_QPW4P1T_2D) .or.&
        (celement .eq. EL_QPW4P1TVDF_2D) .or.&
        (celement .eq. EL_QPW4DCP1_2D)) then
      ! Piecewise linear cubature on sub-triangles on a quad
      ccubType = CUB_QPW4G3T_2D
      return
    end if

    if (celement .eq. EL_QPW4P2_2D) then
      ! Piecewise linear cubature on sub-triangles on a quad
      ccubType = CUB_QPW4G3T_2D
      return
    end if

    select case (elem_igetDimension(celement))
    case (NDIM1D)

      select case (elem_getPrimaryElement(celement))
      case (EL_P0_1D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 2-point Gauss
          ccubType = CUB_G2_1D
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 1-point Gauss
          ccubType = CUB_G1_1D
        end select

      case (EL_P1_1D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 3-point Gauss
          ccubType = CUB_G3_1D
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 2-point Gauss
          ccubType = CUB_G2_1D
        end select

      case (EL_P2_1D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 4-point Gauss
          ccubType = CUB_G4_1D
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 3-point Gauss
          ccubType = CUB_G3_1D
        end select

      case (EL_S31_1D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 5-point Gauss
          ccubType = CUB_G5_1D
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 4-point Gauss
          ccubType = CUB_G4_1D
        end select

      case default
        ccubType = 0
      end select

    case (NDIM2D)

      select case (elem_getPrimaryElement(celement))
      case (EL_P0)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! Use Gauss 3pt
          ccubType = CUB_G3_T
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! Use Gauss 1X1
          ccubType = CUB_G1_T
        end select

      case (EL_P1)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! Use VMC
          ccubType = CUB_VMC
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! Use Gauss-3pt
          ccubType = CUB_G3_T
        end select

      case (EL_P2,EL_P2E)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! Use Gauss-4pt
          ccubType = CUB_VMC
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! Use Gauss-3pt
          ccubType = CUB_G3_T
        end select

      case (EL_P3)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! Use Gauss-4pt
          ccubType = CUB_VMC
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! Use Gauss-4pt
          ccubType = CUB_VMC
        end select

      case (EL_P1T)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! Use Gauss-4pt
          ccubType = CUB_VMC
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! Use Gauss-3pt
          ccubType = CUB_G3_T
        end select

      case (EL_Q0)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 2x2 Gauss formula
          ccubType = CUB_G2X2
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 1x1 Gauss formula
          ccubType = CUB_G1X1
        end select

      case (EL_Q1,EL_Q1TB)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 3x3 Gauss formula
          ccubType = CUB_G3X3
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 2x2 Gauss formula
          ccubType = CUB_G2X2
        end select

      case (EL_Q2,EL_Q2T,EL_Q2TB)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 4x4 Gauss formula
          ccubType = CUB_G4X4
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 3x3 Gauss formula
          ccubType = CUB_G3X3
        end select

      case (EL_Q3,EL_Q3T_2D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 5x5 Gauss formula; in particular needed for EL_Q3T_2D
          ccubType = CUB_G5X5
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 4x4 Gauss formula
          ccubType = CUB_G4X4
        end select

      case (EL_Q1T)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 3x3 Gauss formula
          ccubType = CUB_G3X3
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 2x2 Gauss formula
          ccubType = CUB_G2X2
        end select

      case (EL_QP1, EL_DCQP1_2D, EL_DG_T1_2D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 3x3 Gauss formula
          ccubType = CUB_G3X3
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 2x2 Gauss formula
          ccubType = CUB_G2X2
        end select

      case (EL_DCQP2_2D,EL_DG_T2_2D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 4x4 Gauss formula
          ccubType = CUB_G4X4
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 3x3 Gauss formula
          ccubType = CUB_G3X3
        end select

      case default
        ccubType = 0
      end select

    case (NDIM3D)

      select case (elem_getPrimaryElement(celement))
      case (EL_Q0_3D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 2pt Gauss formula
          ccubType = CUB_G2_3D
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 1pt Gauss formula
          ccubType = CUB_G1_3D
        end select

      case (EL_Q1_3D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 3x3 Gauss formula
          ccubType = CUB_G3_3D
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 2x2 Gauss formula
          ccubType = CUB_G2_3D
        end select

      case (EL_QP1_3D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 2-pt Gauss formula
          ccubType = CUB_G2_3D
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 1-pt Gauss formula
          ccubType = CUB_G1_3D
        end select

      case (EL_Q1T_3D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 3-pt Gauss formula
          ccubType = CUB_G3_3D
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 2-pt Gauss formula
          ccubType = CUB_G2_3D
        end select

      case (EL_Q2T_3D)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 4-pt Gauss formula
          ccubType = CUB_G4_3D
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 3-pt Gauss formula
          ccubType = CUB_G3_3D
        end select

      case default
        ccubType = 0
      end select

    case default
      ccubType = 0
    end select

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_initBlockDiscr_fix (rblockDiscr,ncomponents,&
      rtriangulation, rboundary)

!<description>

  ! This routine initialises a block discretisation structure accept ncomponents
  ! solution components. Pointers to the triangulation, domain and boundary
  ! conditions are saved in the structure.
  !
  ! The routine performs only basic initialisation. The caller must
  ! separately initialise the specific scalar discretisation structures
  ! of each solution component (as collected in the RspatialDiscr
  ! array of the rblockDiscr structure).

!</description>

!<input>

  ! Number of solution components maintained by the block structure
  integer, intent(in)                          :: ncomponents

  ! OPTIONAL: The triangulation structure underlying to the discretisation.
  type(t_triangulation), intent(in), optional, target :: rtriangulation

  ! OPTIONAL: The underlying domain.
  type(t_boundary), intent(in), target, optional :: rboundary

!</input>

!<output>

  ! The block discretisation structure to be initialised.
  type(t_blockDiscretisation), intent(out) :: rblockDiscr

!</output>

!</subroutine>

    ! Initialise the variables of the structure for the simple discretisation
    rblockDiscr%ccomplexity      = SPDISC_UNIFORM
    
    if (present(rtriangulation)) then
      rblockDiscr%ndimension       = rtriangulation%ndim
      rblockDiscr%p_rtriangulation => rtriangulation
    else
      ! Unknown dimension
      rblockDiscr%ndimension       = -1
      nullify(rblockDiscr%p_rtriangulation)
    end if
    
    if (present(rboundary)) then
      rblockDiscr%p_rboundary    => rboundary
    else
      nullify(rblockDiscr%p_rboundary)
    end if

    rblockDiscr%ncomponents      = ncomponents
    allocate(rblockDiscr%RspatialDiscr(ncomponents))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_initBlockDiscr_free (rblockDiscr,ndof)

!<description>
  ! Initialises a "free" block discretisation with one component and ndof
  ! degrees of freedom. This is a general purpose discretisation without
  ! any triangulation/boundary/FEM space attached.
  !
  ! This routine is a replacement for creating a block discretisation
  ! with one component by spdiscr_initBlockDiscr and creating one
  ! component with a "free" scalar discretisation using spdiscr_initDiscr_free.
!</description>

!<input>
  ! Number of degrees of freedom in the space.
  integer, intent(in) :: ndof
!</input>

!<output>
  ! The block discretisation structure to be initialised.
  type(t_blockDiscretisation), intent(out) :: rblockDiscr
!</output>

!</subroutine>

    call spdiscr_initBlockDiscr (rblockDiscr,1)
    call spdiscr_initDiscr_free (rblockDiscr%RspatialDiscr(1),ndof)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_initBlockDiscr_open (rblockDiscr, rtriangulation, rboundary)

!<description>
  ! This routine initialises an "open" block discretisation structure.
  ! Pointers to the triangulation, domain and boundary conditions are 
  ! saved in the structure. The caller can add spatial discretisation
  ! structures using spdiscr_appendBlockComponent. After the discretisation
  ! is finished, spdiscr_commitBlockDiscr shall be called.
!</description>

!<input>

  ! OPTIONAL: The triangulation structure underlying to the discretisation.
  type(t_triangulation), intent(in), optional, target :: rtriangulation

  ! OPTIONAL: The underlying domain.
  type(t_boundary), intent(in), target, optional :: rboundary

!</input>

!<output>

  ! The block discretisation structure to be initialised.
  type(t_blockDiscretisation), intent(out) :: rblockDiscr

!</output>

!</subroutine>

    ! Initialise the variables of the structure for the simple discretisation
    rblockDiscr%ccomplexity      = SPDISC_UNIFORM
    
    if (present(rtriangulation)) then
      rblockDiscr%ndimension       = rtriangulation%ndim
      rblockDiscr%p_rtriangulation => rtriangulation
    else
      ! Unknown dimension
      rblockDiscr%ndimension       = -1
      nullify(rblockDiscr%p_rtriangulation)
    end if
    
    if (present(rboundary)) then
      rblockDiscr%p_rboundary    => rboundary
    else
      nullify(rblockDiscr%p_rboundary)
    end if

    ! Allocate a minimum amount of discretisation structures.
    ! Memory will be reduced later.
    rblockDiscr%ncomponents      = 0
    allocate(rblockDiscr%RspatialDiscr(16))

    ! That is it.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_appendBlockComponent (rblockDiscr, rspatialDiscr, ntensorDim, bshare)

!<description>
  ! Adds a new discretisation structure to the block discretisation.
!</description>

!<input>
  ! Source discretisation structure to be added.
  type(t_spatialDiscretisation), intent(in) :: rspatialDiscr
  
  ! OPTIONAL: Defines a tensor dimension. If specified, rspatialDiscr will
  ! be added ntensorDim times. Default =1.
  integer, intent(in), optional :: ntensorDim
  
  ! OPTIONAL: Defines whether or not the block discretisation should "share"
  ! the data with rspatialDiscr. =TRUE by default.
  logical, intent(in), optional :: bshare
!</input>

!<inputoutput>
  ! The block discretisation structure where to add rspatialDiscr.
  type(t_blockDiscretisation), intent(inout) :: rblockDiscr
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), dimension(:), pointer :: p_RspatialDiscr
    integer :: isize,ntensor,i
    
    ntensor = 1
    if (present(ntensorDim)) ntensor = ntensorDim
    
    do i=1,ntensor
      
      ! Check if there is enough space.
      if (rblockDiscr%ncomponents .ge. ubound (rblockDiscr%RspatialDiscr,1)) then

        ! Reallocate
        isize = ubound (rblockDiscr%RspatialDiscr,1)
        allocate (p_RspatialDiscr(isize+16))
        p_RspatialDiscr (1:isize) = rblockDiscr%RspatialDiscr(1:isize)
        deallocate (rblockDiscr%RspatialDiscr)
        rblockDiscr%RspatialDiscr => p_RspatialDiscr

      end if
      
      ! New block
      rblockDiscr%ncomponents = rblockDiscr%ncomponents + 1
      
      ! Duplicate the spatial discretisation into this block
      call spdiscr_duplicateDiscrSc (rspatialDiscr, &
          rblockDiscr%RspatialDiscr(rblockDiscr%ncomponents), bshare)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_commitBlockDiscr (rblockDiscr)

!<description>
  ! "Finishes" a block discretisation after being set up with 
  ! spdiscr_appendBlockComponent.
!</description>

!<inputoutput>
  ! The block discretisation structure to be committed.
  type(t_blockDiscretisation), intent(inout) :: rblockDiscr
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), dimension(:), pointer :: p_RspatialDiscr
    integer :: isize
    
    ! Reallocate to the correct size
    isize = ubound (rblockDiscr%RspatialDiscr,1)
    if (rblockDiscr%ncomponents .eq. isize) return

    ! Reallocate
    isize = rblockDiscr%ncomponents

    allocate (p_RspatialDiscr(isize))
    p_RspatialDiscr (1:isize) = rblockDiscr%RspatialDiscr(1:isize)
    deallocate (rblockDiscr%RspatialDiscr)
    rblockDiscr%RspatialDiscr => p_RspatialDiscr

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_createBlockDiscrInd (rspatialDiscr,rblockDiscr)

!<description>
  ! This routine creates a block discretisation structure with one block from
  ! a scalar discretisation structure. The scalar discretisation structure
  ! is embedded as a 'shared copy' into the first component of the
  ! block discretisation.
  ! (So, releasing the block discretisation will not destroy the original
  ! spatial discretisation.)
!</description>

!<input>

  ! Spatial discretisation structure that is embedded into the
  ! block discretisation.
  type(t_spatialDiscretisation), intent(in) :: rspatialDiscr

!</input>

!<output>

  ! The block discretisation structure to be initialised.
  type(t_blockDiscretisation), intent(out) :: rblockDiscr

!</output>

!</subroutine>

    ! Initialise a new block discretisation with one component.
    if (associated(rspatialDiscr%p_rboundary)) then
      call spdiscr_initBlockDiscr (rblockDiscr, 1,&
          rspatialDiscr%p_rtriangulation, rspatialDiscr%p_rboundary)
    else
      call spdiscr_initBlockDiscr (rblockDiscr, 1,&
          rspatialDiscr%p_rtriangulation)
    end if

    ! Put a copy of the spatial discretisation to first component
    ! of the the block discretisation. Share the data.
    call spdiscr_duplicateDiscrSc (rspatialDiscr,&
        rblockDiscr%RspatialDiscr(1), .true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_deriveBlockDiscr (rsourceDiscr, rdestDiscr, &
      ifirstBlock, ilastBlock, rtriangulation, rboundary, bshare)

!<description>
  ! This routine derives a block discretisation structure from another one.
  !
  ! rsourceDiscr is a given block discretisation structure.
  ! ifirstBlock is the number of the block in rsourceDiscr that should be
  ! used as first block in rdestDiscr.
  ! ilastBlock is the number of the block in rsourceDiscr that should be
  ! used as last block in rdestDiscr.
  !
  ! rdestDiscr will therefore contain the blocks
  ! ifirstBlock..ilastBlock of rsourceDiscr. No memory will be allocated by
  ! this procedure, rdestDiscr will simply share all handles with
  ! rsourceDiscr.
!</description>

!<input>
  ! A source discretisation structure that should be used as template
  type(t_blockDiscretisation), intent(in), target :: rsourceDiscr

  ! OPTIONAL: Number of the block in rsourceDiscr that should be
  ! used as first block in rdestDiscr. Default value is =1.
  integer, intent(in), optional :: ifirstBlock

  ! OPTIONAL: Number of the last block in rsourceDiscr that should be
  ! used as last block in rdestDiscr. Default value is the
  ! number of components in rsourceDiscr.
  integer, intent(in), optional :: ilastBlock

  ! OPTIONAL: Reference to a new triangulation, the new discretisation should use.
  ! If not specified, the triangulation in rsourceDiscr is used.
  type(t_triangulation), intent(in), target, optional :: rtriangulation

  ! OPTIONAL: Reference to a new domain, the new discretisation should use.
  ! If not specified, the domain in rsourceDiscr is used.
  type(t_boundary), intent(in), target, optional :: rboundary

  ! OPTIONAL: Whether the new discretisation structure should share its information
  ! with rsourceDiscr.
  ! =FALSE: Create a complete copy of a subdiscretisation of rsourceDiscr
  !  which is independent of rsourceDiscr.
  ! =TRUE: The new discretisation will not be a complete new structure, but a
  !  'derived' structure, i.e. it uses the same dynamic information
  !  (handles and therefore element lists) as rsourceDiscr.
  ! If not specified, TRUE is assumed.
  logical, intent(in), optional :: bshare
!</input>

!<inputoutput>
  ! The discretisation structure to be initialised.
  ! Any old existing information in rdestDiscr is released if necessary.
  type(t_blockDiscretisation), intent(inout), target :: rdestDiscr
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ifirst, ilast, ncount, i
    logical :: bshr

    ! Check that the source discretisation structure is valid.
    if (rsourceDiscr%ndimension .le. 0) then
      call output_line ('Source structure invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveBlockDiscr')
      call sys_halt()
    end if

    ! Release old information if present
    call spdiscr_releaseBlockDiscr(rdestDiscr,.true.)

    ! Evaluate the optional parameters
    ifirst = 1
    ilast  = rsourceDiscr%ncomponents

    if (present(ifirstBlock)) then
      ifirst = min(max(ifirst,ifirstBlock),ilast)
    end if

    if (present(ilastBlock)) then
      ilast = max(min(ilast,ilastBlock),ifirst)
    end if

    ncount = ilast-ifirst+1

    bshr = .true.
    if (present(bshare)) bshr = bshare

    ! Copy all information from the source discretisation structure

    rdestDiscr%ndimension       =  rsourceDiscr%ndimension
    rdestDiscr%ccomplexity      =  rsourceDiscr%ccomplexity
    rdestDiscr%p_rboundary      => rsourceDiscr%p_rboundary
    rdestDiscr%p_rtriangulation => rsourceDiscr%p_rtriangulation
    if (present(rboundary)) rdestDiscr%p_rboundary => rsourceDiscr%p_rboundary
    if (present(rtriangulation)) rdestDiscr%p_rtriangulation => rsourceDiscr%p_rtriangulation
    rdestDiscr%ncomponents      =  ncount

    ! Copy all substructures -- from ifirstBlock to ilastBlock.
    ! Use spdiscr_duplicateDiscrSc which savely copies the scalar discretisation
    ! structures. We set bshare=.TRUE. here, so the information is shared
    ! between the source and destination structure; the dynamic information
    ! 'belongs' to rdiscrSource and not to the newly created rdiscrDest!
    allocate(rdestDiscr%RspatialDiscr(ncount))
    do i = 1, ncount
      call spdiscr_duplicateDiscrSc (rsourceDiscr%RspatialDiscr(ifirst+i-1), &
                                     rdestDiscr%RspatialDiscr(i), bshr)
    end do

    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_releaseBlockDiscr (rblockDiscr, breleaseSubstruc)

!<description>
  ! This routine releases a block discretisation structure from memory.
!</description>

!<input>
  ! OPTIONAL: Release substructures.
  ! If set to TRUE, the memory of all scalar spatial discretisation structures
  !   in rblockDiscr is also released from memory. This is the standard setting.
  ! Is set to FALSE, only rblockDiscr is cleaned up, the substructures
  !   are ignored.
  logical, intent(in), optional :: breleaseSubstruc
!</input>

!<inputoutput>
  ! The block discretisation structures to be released.
  type(t_blockDiscretisation), intent(inout), target :: rblockDiscr
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  logical :: brelsub

  if (rblockDiscr%ndimension .ne. 0) then
    brelsub = .true.
    if (present(breleaseSubstruc)) brelsub = breleaseSubstruc

    ! Cut the connection to the other structures
    nullify(rblockDiscr%p_rtriangulation)
    nullify(rblockDiscr%p_rboundary)

    ! Release substructures?
    if (associated(rblockDiscr%RspatialDiscr)) then
      if (brelsub) then
        do i = 1, rblockDiscr%ncomponents
          call spdiscr_releaseDiscr (rblockDiscr%RspatialDiscr(i))
        end do
      end if
      deallocate(rblockDiscr%RspatialDiscr)
    end if
    rblockDiscr%ncomponents = 0

    ! Structure not initialised anymore
    rblockDiscr%ndimension = 0
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_initDiscr_simple_old (rspatialDiscr,celement, ccubType,&
                                           rtriangulation, rboundary)

!<description>
  ! This routine initialises a discretisation structure for a uniform
  ! discretisation with one element for all geometric element primitives,
  ! for trial as well as for test functions.
  !
  ! If rspatialDiscr is NULL(), a new structure will be created. Otherwise,
  ! the existing structure is recreated/updated.
!</description>

!<input>
  ! The element type identifier that is to be used for all elements.
  integer(I32), intent(in) :: celement

  ! Cubature formula CUB_xxxx to use for calculating integrals.
  ! Alternatively, the value CUB_GEN_AUTO means:
  ! automatically determine cubature formula.
  integer(I32), intent(in) :: ccubType

  ! The triangulation structure underlying to the discretisation.
  type(t_triangulation), intent(in), target :: rtriangulation

  ! The underlying domain.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<intputoutput>
  ! The discretisation structure to be initialised.
  type(t_spatialDiscretisation), intent(inout), target :: rspatialDiscr
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  integer, dimension(:), pointer :: p_Iarray
  type(t_elementDistribution), pointer :: p_relementDistr
  integer(I32) :: ccub

  ! Automatically determine cubature formula if necessary
  ccub = ccubType

  if (ccub .eq. CUB_GEN_AUTO) &
      ccub = spdiscr_getStdCubature(celement)

  if (ccub .eq. CUB_GEN_AUTO_LUMPMASS) &
      ccub = spdiscr_getLumpCubature(celement)

  ! Do we have a structure?
  if (rspatialDiscr%ndimension .ne. 0) then
    ! Release the old structure.
    call spdiscr_releaseDiscr(rspatialDiscr)
  end if

  ! Initialise the variables of the structure for the simple discretisation
  rspatialDiscr%ndimension       =  rtriangulation%ndim
  rspatialDiscr%p_rtriangulation => rtriangulation
  if (present(rboundary)) then
    rspatialDiscr%p_rboundary    => rboundary
  else
    nullify(rspatialDiscr%p_rboundary)
  end if
  rspatialDiscr%ccomplexity      =  SPDISC_UNIFORM

  ! All trial elements are celement:

!  CALL storage_new ('spdiscr_initDiscr_simple', 'h_ItrialElements', &
!        rtriangulation%NEL, ST_INT, rspatialDiscr%h_ItrialElements,   &
!        ST_NEWBLOCK_NOINIT)
!  CALL storage_getbase_int (rspatialDiscr%h_ItrialElements,p_Iarray)
!  DO i=1,rtriangulation%NEL
!    p_Iarray(i) = celement
!  END DO
  rspatialDiscr%h_IelemGroupIDs = ST_NOHANDLE

  ! Initialise the first element group
  rspatialDiscr%inumFESpaces = 1
  allocate(rspatialDiscr%RelementDistr(rspatialDiscr%inumFESpaces))
  p_relementDistr => rspatialDiscr%RelementDistr(1)

  ! Initialise FE space for that block
  p_relementDistr%celement        = celement
  p_relementDistr%ccubTypeBilForm = ccub
  p_relementDistr%ccubTypeLinForm = ccub
  p_relementDistr%ccubTypeEval    = ccub

  ! Get the typical transformation used with the element
  p_relementDistr%ctrafoType = elem_igetTrafoType(celement)

  ! Check the cubature formula against the element group.
  ! This stops the program if this is not fulfilled.
  call spdiscr_checkCubature(ccub,celement)

  ! Initialise an 'identity' array containing the numbers of all elements.
  ! This list defines the sequence how elements are processed, e.g. in the
  ! assembly of matrices/vectors.
  call storage_new ('spdiscr_initDiscr_simple', 'h_IelementList', &
      rtriangulation%NEL, ST_INT, p_relementDistr%h_IelementList,   &
      ST_NEWBLOCK_NOINIT)
  call storage_getbase_int (p_relementDistr%h_IelementList,p_Iarray)
  do i = 1, rtriangulation%NEL
    p_Iarray(i) = i
  end do

  ! Save the number of elements in that element list.
  p_relementDistr%NEL = rtriangulation%NEL

  ! This is a complete new structure, everything 'belongs' to this.
  rspatialDiscr%bisCopy = .false.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_initDiscr_simple_new (rspatialDiscr,celement,&
                                           rtriangulation, rboundary)

!<description>
  ! This routine initialises a discretisation structure for a uniform
  ! discretisation with one element for all geometric element primitives,
  ! for trial as well as for test functions.
  !
  ! If rspatialDiscr is NULL(), a new structure will be created. Otherwise,
  ! the existing structure is recreated/updated.
!</description>

!<input>
  ! The element type identifier that is to be used for all elements.
  integer(I32), intent(in) :: celement

  ! The triangulation structure underlying to the discretisation.
  type(t_triangulation), intent(in), target :: rtriangulation

  ! The underlying domain.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<inputoutput>
  ! The discretisation structure to be initialised.
  type(t_spatialDiscretisation), intent(inout), target :: rspatialDiscr
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    integer, dimension(:), pointer :: p_Iarray
    type(t_elementDistribution), pointer :: p_relementDistr

    ! Do we have a structure?
    if (rspatialDiscr%ndimension .ne. 0) then
      ! Release the old structure.
      call spdiscr_releaseDiscr(rspatialDiscr)
    end if

    ! Initialise the variables of the structure for the simple discretisation
    rspatialDiscr%ndimension       =  rtriangulation%ndim
    rspatialDiscr%p_rtriangulation => rtriangulation
    if (present(rboundary)) then
      rspatialDiscr%p_rboundary    => rboundary
    else
      nullify(rspatialDiscr%p_rboundary)
    end if
    rspatialDiscr%ccomplexity      =  SPDISC_UNIFORM

    ! All trial elements are celement:

  !  CALL storage_new ('spdiscr_initDiscr_simple', 'h_ItrialElements', &
  !        rtriangulation%NEL, ST_INT, rspatialDiscr%h_ItrialElements,   &
  !        ST_NEWBLOCK_NOINIT)
  !  CALL storage_getbase_int (rspatialDiscr%h_ItrialElements,p_Iarray)
  !  DO i=1,rtriangulation%NEL
  !    p_Iarray(i) = celement
  !  END DO
    rspatialDiscr%h_IelemGroupIDs = ST_NOHANDLE

    ! Initialise the first element group
    rspatialDiscr%inumFESpaces = 1
    allocate(rspatialDiscr%RelementDistr(rspatialDiscr%inumFESpaces))
    p_relementDistr => rspatialDiscr%RelementDistr(1)

    ! Initialise FE space for that block
    p_relementDistr%celement        = celement

    ! Get the typical transformation used with the element
    p_relementDistr%ctrafoType = elem_igetTrafoType(celement)

    ! Initialise an 'identity' array containing the numbers of all elements.
    ! This list defines the sequence how elements are processed, e.g. in the
    ! assembly of matrices/vectors.
    call storage_new ('spdiscr_initDiscr_simple', 'h_IelementList', &
        rtriangulation%NEL, ST_INT, p_relementDistr%h_IelementList,   &
        ST_NEWBLOCK_NOINIT)
    call storage_getbase_int (p_relementDistr%h_IelementList,p_Iarray)
    do i = 1, rtriangulation%NEL
      p_Iarray(i) = i
    end do

    ! Save the number of elements in that element list.
    p_relementDistr%NEL = rtriangulation%NEL

    ! This is a complete new structure, everything 'belongs' to this.
    rspatialDiscr%bisCopy = .false.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_initDiscr_free (rspatialDiscr,ndof)

!<description>
  ! Initialises a "free" discretisation. A "free" discretisation is a
  ! discretisation which creates a solution space with ndof degrees of
  ! freedom without any domain/triangulation/FEM space attached.
  ! This can be used as a "general purpose" discretisation.
!</description>

!<input>
  ! Number of degrees of freedom
  integer, intent(in) :: ndof
!</input>

!<inputoutput>
  ! The discretisation structure to be initialised.
  type(t_spatialDiscretisation), intent(inout), target :: rspatialDiscr
!</inputoutput>

!</subroutine>

    ! Do we have a structure?
    if (rspatialDiscr%ndimension .ne. 0) then
      ! Release the old structure.
      call spdiscr_releaseDiscr(rspatialDiscr)
    end if

    ! Initialise the variables of the structure for the simple discretisation
    rspatialDiscr%ndimension  = -1
    rspatialDiscr%ccomplexity = SPDISC_FREE
    rspatialDiscr%ndof        = ndof

    ! This is a complete new structure, everything 'belongs' to this.
    rspatialDiscr%bisCopy = .false.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_initDiscr_triquad_old (rspatialDiscr, ieltyptri, ieltypquad,&
                                        ccubTypeTri, ccubTypeQuad,&
                                        rtriangulation, rboundary)

!<description>
  ! This routine initialises a discretisation structure for a conformal
  ! discretisation, mixed triangular/quad mesh with one element type for all
  ! triangles and one element type for all quads -- for trial as well as
  ! for test functions.
  !
  ! If rspatialDiscr is NULL(), a new structure will be created. Otherwise,
  ! the existing structure is recreated/updated.
!</description>

!<input>
  ! The element type identifier that is to be used for all triangular elements.
  integer(I32), intent(in) :: ieltypTri

  ! The element type identifier that is to be used for all quadrilateral elements.
  integer(I32), intent(in) :: ieltypQuad

  ! Cubature formula CUB_xxxx to use for calculating integrals
  ! on triangular elements
  ! Alternatively, the value CUB_GEN_AUTO means:
  ! automatically determine cubature formula.
  integer(I32), intent(in) :: ccubTypeTri

  ! Cubature formula CUB_xxxx to use for calculating integrals on
  ! quadrilateral elements
  ! Alternatively, the value CUB_GEN_AUTO means:
  ! automatically determine cubature formula.
  integer(I32), intent(in) :: ccubTypeQuad

  ! The triangulation structure underlying to the discretisation.
  type(t_triangulation), intent(in), target :: rtriangulation

  ! The underlying domain.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<inputoutput>
  ! The discretisation structure to be initialised.
  type(t_spatialDiscretisation), intent(inout), target :: rspatialDiscr
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,j
  integer, dimension(2) :: IelemCount
  integer, dimension(:), pointer :: p_Iarray,p_IelementCounter
  type(t_elementDistribution), pointer :: p_relementDistrTria,p_relementDistrQuad
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer(I32) :: ccubTri,ccubQuad

  ! Automatically determine cubature formula if necessary
  ccubTri = ccubTypeTri

  if (ccubTri .eq. CUB_GEN_AUTO) &
      ccubTri = spdiscr_getStdCubature(ieltypTri)

  if (ccubTri .eq. CUB_GEN_AUTO_LUMPMASS) &
      ccubTri = spdiscr_getLumpCubature(ieltypTri)

  ccubQuad = ccubTypeQuad

  if (ccubQuad .eq. CUB_GEN_AUTO) &
      ccubQuad = spdiscr_getStdCubature(ieltypQuad)

  if (ccubQuad .eq. CUB_GEN_AUTO_LUMPMASS) &
      ccubQuad = spdiscr_getLumpCubature(ieltypQuad)

  ! Do we have a structure?
  if (rspatialDiscr%ndimension .ne. 0) then
    ! Release the old structure.
    call spdiscr_releaseDiscr(rspatialDiscr)
  end if

  ! Initialise the variables of the structure for the simple discretisation
  rspatialDiscr%ndimension       = NDIM2D
  rspatialDiscr%p_rtriangulation => rtriangulation
  if (present(rboundary)) then
    rspatialDiscr%p_rboundary    => rboundary
  else
    nullify(rspatialDiscr%p_rboundary)
  end if
  rspatialDiscr%ccomplexity      = SPDISC_CONFORMAL

  ! Allocate an array containing the element group for each element
  call storage_new ('spdiscr_initDiscr_triquad', 'h_ItrialElements', &
      rtriangulation%NEL, ST_INT, rspatialDiscr%h_IelemGroupIDs,   &
      ST_NEWBLOCK_NOINIT)
  call storage_getbase_int (rspatialDiscr%h_IelemGroupIDs,p_Iarray)

  ! Allocate an array with an element counter for every element type.
  call storage_new ('spdiscr_initDiscr_triquad', 'h_IelementCounter', &
      rtriangulation%NEL, ST_INT, rspatialDiscr%h_IelementCounter,   &
      ST_NEWBLOCK_NOINIT)
  call storage_getbase_int (rspatialDiscr%h_IelementCounter,p_IelementCounter)

  ! Create both arrays simultaneously.
  call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

  IelemCount(:) = 0
  if (ubound(p_IverticesAtElement,1) .ge. 4) then
    ! There are quads and probably triangles in the mesh
    do i = 1, rtriangulation%NEL
      if (p_IverticesAtElement (4,i) .eq. 0) then
        ! Triangular elements are in element group 1
        p_Iarray(i) = 1

        ! This is the IelemCount(1)-th triangle
        IelemCount(1) = IelemCount(1)+1
        p_IelementCounter(i) = IelemCount(1)
      else
        ! Quad elements are in element group 2
        p_Iarray(i) = 2

        ! This is the IelemCount(2)-th quad
        IelemCount(2) = IelemCount(2)+1
        p_IelementCounter(i) = IelemCount(2)
      end if
    end do
  else
    ! Pure triangular mesh
    do i = 1, rtriangulation%NEL
      ! Triangular elements are in element group 1
      p_Iarray(i) = 1

      ! This is the IelemCount(1)-th triangle
      IelemCount(1) = IelemCount(1)+1
      p_IelementCounter(i) = IelemCount(1)
    end do
  end if

  ! Initialise the first element group
  rspatialDiscr%inumFESpaces = 2
  allocate(rspatialDiscr%RelementDistr(rspatialDiscr%inumFESpaces))
  p_relementDistrTria => rspatialDiscr%RelementDistr(1)
  p_relementDistrQuad => rspatialDiscr%RelementDistr(2)

  ! Initialise test and trial space for that block
  p_relementDistrTria%celement        = ieltypTri
  p_relementDistrTria%ccubTypeBilForm = ccubTri
  p_relementDistrTria%ccubTypeLinForm = ccubTri
  p_relementDistrTria%ccubTypeEval    = ccubTri

  p_relementDistrQuad%celement        = ieltypQuad
  p_relementDistrQuad%ccubTypeBilForm = ccubQuad
  p_relementDistrQuad%ccubTypeLinForm = ccubQuad
  p_relementDistrQuad%ccubTypeEval    = ccubQuad

  ! Get the typical transformation used with the element
  p_relementDistrTria%ctrafoType = elem_igetTrafoType(ieltypTri)
  p_relementDistrQuad%ctrafoType = elem_igetTrafoType(ieltypQuad)

  ! Check the cubature formula against the element group.
  ! This stops the program if this is not fulfilled.
  call spdiscr_checkCubature(ccubTri,ieltypTri)
  call spdiscr_checkCubature(ccubQuad,ieltypQuad)

  ! Save the number of elements in the two element lists.
  p_relementDistrTria%NEL = rtriangulation%InelOfType(TRIA_NVETRI2D)
  p_relementDistrQuad%NEL = rtriangulation%InelOfType(TRIA_NVEQUAD2D)

  ! Initialise an 'identity' array containing the numbers of all elements.
  ! This list defines the sequence how elements are processed, e.g. in the
  ! assembly of matrices/vectors.

  ! We have to collect all triangles to the first and all quads to the second
  ! element group. j counts how many elements we found
  !
  ! Collect all triangles
  j = 0
  if (rtriangulation%InelOfType(TRIA_NVETRI2D) .ne. 0) then

    call storage_new ('spdiscr_initDiscr_triquad', 'h_IelementList', &
          rtriangulation%InelOfType(TRIA_NVETRI2D), &
          ST_INT, p_relementDistrTria%h_IelementList,   &
          ST_NEWBLOCK_NOINIT)

    call storage_getbase_int (p_relementDistrTria%h_IelementList,p_Iarray)

    if (ubound(p_IverticesAtElement,1) .ge. TRIA_NVEQUAD2D) then
      ! There are quads and probably triangles in the mesh
      do i = 1, rtriangulation%NEL
        if (p_IverticesAtElement(TRIA_NVEQUAD2D,i) .eq. 0) then
          j = j+1
          p_Iarray(j) = i
        end if
      end do
    else
      ! Pure triangular mesh
      do i = 1, rtriangulation%NEL
        j = j+1
        p_Iarray(j) = i
      end do
    end if
  end if

  ! Collect all quads
  j = 0
  if (rtriangulation%InelOfType(TRIA_NVEQUAD2D) .ne. 0) then

    call storage_new ('spdiscr_initDiscr_triquad', 'h_IelementList', &
          rtriangulation%InelOfType(TRIA_NVEQUAD2D), &
          ST_INT, p_relementDistrQuad%h_IelementList,   &
          ST_NEWBLOCK_NOINIT)

    call storage_getbase_int (p_relementDistrQuad%h_IelementList,p_Iarray)

    ! Because of the IF above, there are for sure quads in the mesh!
    do i = 1, rtriangulation%NEL
      if (p_IverticesAtElement(4,i) .ne. 0) then
        j = j+1
        p_Iarray(j) = i
      end if
    end do

  end if

  rspatialDiscr%bisCopy = .false.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_initDiscr_triquad_new (rspatialDiscr, ieltyptri, ieltypquad,&
                                            rtriangulation, rboundary)

!<description>
  ! This routine initialises a discretisation structure for a conformal
  ! discretisation, mixed triangular/quad mesh with one element type for all
  ! triangles and one element type for all quads -- for trial as well as
  ! for test functions.
  !
  ! If rspatialDiscr is NULL(), a new structure will be created. Otherwise,
  ! the existing structure is recreated/updated.
!</description>

!<input>
  ! The element type identifier that is to be used for all triangular elements.
  integer(I32), intent(in) :: ieltypTri

  ! The element type identifier that is to be used for all quadrilateral elements.
  integer(I32), intent(in) :: ieltypQuad

  ! The triangulation structure underlying to the discretisation.
  type(t_triangulation), intent(in), target :: rtriangulation

  ! The underlying domain.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<inpuoutput>
  ! The discretisation structure to be initialised.
  type(t_spatialDiscretisation), intent(inout), target :: rspatialDiscr
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,j
  integer, dimension(2) :: IelemCount
  integer, dimension(:), pointer :: p_Iarray,p_IelementCounter
  type(t_elementDistribution), pointer :: p_relementDistrTria,p_relementDistrQuad
  integer, dimension(:,:), pointer :: p_IverticesAtElement

  ! Do we have a structure?
  if (rspatialDiscr%ndimension .ne. 0) then
    ! Release the old structure.
    call spdiscr_releaseDiscr(rspatialDiscr)
  end if

  ! Initialise the variables of the structure for the simple discretisation
  rspatialDiscr%ndimension       = NDIM2D
  rspatialDiscr%p_rtriangulation => rtriangulation
  if (present(rboundary)) then
    rspatialDiscr%p_rboundary    => rboundary
  else
    nullify(rspatialDiscr%p_rboundary)
  end if
  rspatialDiscr%ccomplexity      = SPDISC_CONFORMAL

  ! Allocate an array containing the element group for each element
  call storage_new ('spdiscr_initDiscr_triquad', 'h_ItrialElements', &
      rtriangulation%NEL, ST_INT, rspatialDiscr%h_IelemGroupIDs,   &
      ST_NEWBLOCK_NOINIT)
  call storage_getbase_int (rspatialDiscr%h_IelemGroupIDs,p_Iarray)

  ! Allocate an array with an element counter for every element type.
  call storage_new ('spdiscr_initDiscr_triquad', 'h_IelementCounter', &
      rtriangulation%NEL, ST_INT, rspatialDiscr%h_IelementCounter,   &
      ST_NEWBLOCK_NOINIT)
  call storage_getbase_int (rspatialDiscr%h_IelementCounter,p_IelementCounter)

  ! Create both arrays simultaneously.
  call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

  IelemCount(:) = 0
  if (ubound(p_IverticesAtElement,1) .ge. 4) then
    ! There are quads and probably triangles in the mesh
    do i = 1, rtriangulation%NEL
      if (p_IverticesAtElement (4,i) .eq. 0) then
        ! Triangular elements are in element group 1
        p_Iarray(i) = 1

        ! This is the IelemCount(1)-th triangle
        IelemCount(1) = IelemCount(1)+1
        p_IelementCounter(i) = IelemCount(1)
      else
        ! Quad elements are in element group 2
        p_Iarray(i) = 2

        ! This is the IelemCount(2)-th quad
        IelemCount(2) = IelemCount(2)+1
        p_IelementCounter(i) = IelemCount(2)
      end if
    end do
  else
    ! Pure triangular mesh
    do i = 1, rtriangulation%NEL
      ! Triangular elements are in element group 1
      p_Iarray(i) = 1

      ! This is the IelemCount(1)-th triangle
      IelemCount(1) = IelemCount(1)+1
      p_IelementCounter(i) = IelemCount(1)
    end do
  end if

  ! Initialise the first element group
  rspatialDiscr%inumFESpaces = 2
  allocate(rspatialDiscr%RelementDistr(rspatialDiscr%inumFESpaces))
  p_relementDistrTria => rspatialDiscr%RelementDistr(1)
  p_relementDistrQuad => rspatialDiscr%RelementDistr(2)

  ! Initialise test and trial space for that block
  p_relementDistrTria%celement        = ieltypTri
  p_relementDistrQuad%celement        = ieltypQuad

  ! Get the typical transformation used with the element
  p_relementDistrTria%ctrafoType = elem_igetTrafoType(ieltypTri)
  p_relementDistrQuad%ctrafoType = elem_igetTrafoType(ieltypQuad)

  ! Save the number of elements in the two element lists.
  p_relementDistrTria%NEL = rtriangulation%InelOfType(TRIA_NVETRI2D)
  p_relementDistrQuad%NEL = rtriangulation%InelOfType(TRIA_NVEQUAD2D)

  ! Initialise an 'identity' array containing the numbers of all elements.
  ! This list defines the sequence how elements are processed, e.g. in the
  ! assembly of matrices/vectors.

  ! We have to collect all triangles to the first and all quads to the second
  ! element group. j counts how many elements we found
  !
  ! Collect all triangles
  j = 0
  if (rtriangulation%InelOfType(TRIA_NVETRI2D) .ne. 0) then

    call storage_new ('spdiscr_initDiscr_triquad', 'h_IelementList', &
          rtriangulation%InelOfType(TRIA_NVETRI2D), &
          ST_INT, p_relementDistrTria%h_IelementList,   &
          ST_NEWBLOCK_NOINIT)

    call storage_getbase_int (p_relementDistrTria%h_IelementList,p_Iarray)

    if (ubound(p_IverticesAtElement,1) .ge. TRIA_NVEQUAD2D) then
      ! There are quads and probably triangles in the mesh
      do i = 1, rtriangulation%NEL
        if (p_IverticesAtElement(TRIA_NVEQUAD2D,i) .eq. 0) then
          j = j+1
          p_Iarray(j) = i
        end if
      end do
    else
      ! Pure triangular mesh
      do i = 1, rtriangulation%NEL
        j = j+1
        p_Iarray(j) = i
      end do
    end if
  end if

  ! Collect all quads
  j = 0
  if (rtriangulation%InelOfType(TRIA_NVEQUAD2D) .ne. 0) then

    call storage_new ('spdiscr_initDiscr_triquad', 'h_IelementList', &
          rtriangulation%InelOfType(TRIA_NVEQUAD2D), &
          ST_INT, p_relementDistrQuad%h_IelementList,   &
          ST_NEWBLOCK_NOINIT)

    call storage_getbase_int (p_relementDistrQuad%h_IelementList,p_Iarray)

    ! Because of the IF above, there are for sure quads in the mesh!
    do i = 1, rtriangulation%NEL
      if (p_IverticesAtElement(4,i) .ne. 0) then
        j = j+1
        p_Iarray(j) = i
      end if
    end do

  end if

  rspatialDiscr%bisCopy = .false.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_deriveSimpleDiscrSc_old (rsourceDiscr, celement, ccubType, &
                                              rdestDiscr)

!<description>
  ! This routine derives a discretisation structure from another one.
  !
  ! rsourceDiscr is a given uniform discretisation structure. The structure
  ! will be copied to rdestDiscr and changed in such a way, that the element
  ! type and cubature formula are changed according to the parameters.
  !
  ! The new discretisation will also be a uniform discretisation based
  ! on the element celement. It is not a complete new structure, but a
  ! 'derived' structure, i.e. it uses the same dynamic information
  ! (element lists) as rsourceDiscr.
!</description>

!<input>
  ! A source discretisation structure that should be used as template
  type(t_spatialDiscretisation), intent(in), target :: rsourceDiscr

  ! The element type identifier that is to be used for all elements
  ! in the new discretisation structure
  integer(I32), intent(in) :: celement

  ! Cubature formula to use for calculating integrals
  ! in the new discretisation structure
  ! Alternatively, the value CUB_GEN_AUTO means:
  ! automatically determine cubature formula.
  ! A value SPDISC_CUB_NOCHANGE means:
  ! take the cubature formula from the source discretisation.
  integer(I32), intent(in) :: ccubType
!</input>

!<inputoutput>
  ! The discretisation structure to be initialised.
  ! Any old existing discretisation information in rdestDiscr
  ! is released if necessary.
  type(t_spatialDiscretisation), intent(inout), target :: rdestDiscr
!</inputoutput>

!</subroutine>

  ! local variables
  ! INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
  ! TYPE(t_elementDistribution), POINTER :: p_relementDistr
  integer(I32) :: ccub

  ! Automatically determine cubature formula if necessary
  ccub = ccubType

  if (ccub .eq. CUB_GEN_AUTO) &
      ccub = spdiscr_getStdCubature(celement)

  if (ccub .eq. CUB_GEN_AUTO_LUMPMASS) &
      ccub = spdiscr_getLumpCubature(celement)

  ! Check that the source discretisation structure is valid.
  if (rsourceDiscr%ndimension .le. 0) then
    call output_line ('Source structure invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveSimpleDiscr')
    call sys_halt()
  end if

  ! Check that the discretisation structure is really uniform.
  ! More complex situations are not supported by this routine.
  if (rsourceDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
    call output_line ('Only uniform discretisations supported!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveSimpleDiscr')
    call sys_halt()
  end if

  if (elem_igetDimension(rsourceDiscr%RelementDistr(1)%celement) .ne. &
      elem_igetDimension(celement)) then
    call output_line ('Element dimension different!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveSimpleDiscr')
    call sys_halt()
  end if

  ! Release old information if present
  call spdiscr_releaseDiscr(rdestDiscr)

  if (ccub .ne. SPDISC_CUB_NOCHANGE) then
    ! Check the cubature formula against the element group.
    ! This stops the program if this is not fulfilled.
    call spdiscr_checkCubature(ccub,celement)
  end if

  ! Copy the source structure to the destination.
  ! This copies all handles and hence all dynamic information
  rdestDiscr = rsourceDiscr

  ! Allocate a new element group and copy content from source
  allocate(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
  rdestDiscr%RelementDistr(1:rdestDiscr%inumFESpaces) = &
      rsourceDiscr%RelementDistr(1:rsourceDiscr%inumFESpaces)

  ! Change the element type of all trial functions to celement
  rdestDiscr%RelementDistr(1)%celement = celement

  ! Init the cubature rule
  if (ccub .eq. SPDISC_CUB_NOCHANGE) then
    ! Copy the old cubature formula
    rdestDiscr%RelementDistr(1)%ccubTypeBilForm = &
        rsourceDiscr%RelementDistr(1)%ccubTypeBilForm
    rdestDiscr%RelementDistr(1)%ccubTypeLinForm = &
        rsourceDiscr%RelementDistr(1)%ccubTypeLinForm
    rdestDiscr%RelementDistr(1)%ccubTypeEval    = &
        rsourceDiscr%RelementDistr(1)%ccubTypeEval
  else
    rdestDiscr%RelementDistr(1)%ccubTypeBilForm = ccub
    rdestDiscr%RelementDistr(1)%ccubTypeLinForm = ccub
    rdestDiscr%RelementDistr(1)%ccubTypeEval = ccub
  end if

  ! Get the typical transformation used with the element
  rdestDiscr%RelementDistr(1)%ctrafoType = elem_igetTrafoType(celement)

  ! Mark the new discretisation structure as 'copy', to prevent
  ! the dynamic information to be released.
  ! The dynamic information 'belongs' to rdiscrSource and not to the
  ! newly created rdiscrDest!
  rdestDiscr%bisCopy = .true.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_deriveSimpleDiscrSc_new (rsourceDiscr, celement, rdestDiscr)

!<description>
  ! This routine derives a discretisation structure from another one.
  !
  ! rsourceDiscr is a given uniform discretisation structure. The structure
  ! will be copied to rdestDiscr and changed in such a way, that the element
  ! type and cubature formula are changed according to the parameters.
  !
  ! The new discretisation will also be a uniform discretisation based
  ! on the element celement. It is not a complete new structure, but a
  ! 'derived' structure, i.e. it uses the same dynamic information
  ! (element lists) as rsourceDiscr.
!</description>

!<input>
  ! A source discretisation structure that should be used as template
  type(t_spatialDiscretisation), intent(in), target :: rsourceDiscr

  ! The element type identifier that is to be used for all elements
  ! in the new discretisation structure
  integer(I32), intent(in) :: celement
!</input>

!<inputoutput>
  ! The discretisation structure to be initialised.
  ! Any old existing discretisation information in rdestDiscr
  ! is released if necessary.
  type(t_spatialDiscretisation), intent(inout), target :: rdestDiscr
!</inputoutput>

!</subroutine>

  ! local variables
  ! INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
  ! TYPE(t_elementDistribution), POINTER :: p_relementDistr

  ! Check that the source discretisation structure is valid.
  if (rsourceDiscr%ndimension .le. 0) then
    call output_line ('Source structure invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveSimpleDiscr')
    call sys_halt()
  end if

  ! Check that the discretisation structure is really uniform.
  ! More complex situations are not supported by this routine.
  if (rsourceDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
    call output_line ('Only uniform discretisations supported!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveSimpleDiscr')
    call sys_halt()
  end if

  if (elem_igetDimension(rsourceDiscr%RelementDistr(1)%celement) .ne. &
      elem_igetDimension(celement)) then
    call output_line ('Element dimension different!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveSimpleDiscr')
    call sys_halt()
  end if

  ! Release old information if present
  call spdiscr_releaseDiscr(rdestDiscr)

  ! Copy the source structure to the destination.
  ! This copies all handles and hence all dynamic information
  rdestDiscr = rsourceDiscr

  ! Allocate a new element group and copy content from source
  allocate(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
  rdestDiscr%RelementDistr(1:rdestDiscr%inumFESpaces) = &
      rsourceDiscr%RelementDistr(1:rsourceDiscr%inumFESpaces)

  ! Change the element type of all trial functions to celement
  rdestDiscr%RelementDistr(1)%celement = celement

  ! Get the typical transformation used with the element
  rdestDiscr%RelementDistr(1)%ctrafoType = elem_igetTrafoType(celement)

  ! Mark the new discretisation structure as 'copy', to prevent
  ! the dynamic information to be released.
  ! The dynamic information 'belongs' to rdiscrSource and not to the
  ! newly created rdiscrDest!
  rdestDiscr%bisCopy = .true.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_deriveDiscr_triquad_old (rsourceDiscr, ieltypTri, ieltypQuad,&
                                              ccubTypeTri, ccubTypeQuad, rdestDiscr)

!<description>
  ! This routine derives a discretisation structure from another one.
  !
  ! rsourceDiscr is a given uniform discretisation structure. The structure
  ! will be copied to rdestDiscr and changed in such a way, that the element
  ! type and cubature formula are changed according to the parameters.
  !
  ! The new discretisation will also be a uniform discretisation based
  ! on the element celement. It is not a complete new structure, but a
  ! 'derived' structure, i.e. it uses the same dynamic information
  ! (element lists) as rsourceDiscr.
!</description>

!<input>
  ! A source discretisation structure that should be used as template
  type(t_spatialDiscretisation), intent(in), target :: rsourceDiscr

  ! The element type identifier that is to be used for all triangular
  ! elements in the new discretisation structure
  integer(I32), intent(in) :: ieltypTri

  ! The element type identifier that is to be used for all quad
  ! elements in the new discretisation structure
  integer(I32), intent(in) :: ieltypQuad

  ! Cubature formula to use for calculating integrals on triangular
  ! elements in the new discretisation structure.
  ! Alternatively, the value CUB_GEN_AUTO means:
  ! automatically determine cubature formula.
  integer(I32), intent(in) :: ccubTypeTri

  ! Cubature formula to use for calculating integrals on quad
  ! elements in the new discretisation structure.
  ! Alternatively, the value CUB_GEN_AUTO means:
  ! automatically determine cubature formula.
  integer(I32), intent(in) :: ccubTypeQuad
!</input>

!<inputoutput>
  ! The discretisation structure to be initialised.
  ! Any old existing discretisation information in rdestDiscr
  ! is released if necessary.
  type(t_spatialDiscretisation), intent(inout), target :: rdestDiscr
!</inputoutput>

!</subroutine>

    ! local variables
    ! INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
    ! TYPE(t_elementDistribution), POINTER :: p_relementDistr
    integer(I32) :: ccubTri,ccubQuad
    integer :: idistr,nve

    ! Automatically determine cubature formula if necessary
    ccubTri = ccubTypeTri

    if (ccubTri .eq. CUB_GEN_AUTO) &
        ccubTri = spdiscr_getStdCubature(ieltypTri)

    if (ccubTri .eq. CUB_GEN_AUTO_LUMPMASS) &
        ccubTri = spdiscr_getLumpCubature(ieltypTri)

    ccubQuad = ccubTypeQuad

    if (ccubQuad .eq. CUB_GEN_AUTO) &
        ccubQuad = spdiscr_getStdCubature(ieltypQuad)

    if (ccubQuad .eq. CUB_GEN_AUTO_LUMPMASS) &
        ccubQuad = spdiscr_getLumpCubature(ieltypQuad)

    ! Check that the source discretisation structure is valid.
    if (rsourceDiscr%ndimension .ne. NDIM2D) then
      call output_line ('Source structure invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdiscr_deriveSimpleDiscr_triquad')
      call sys_halt()
    end if

    ! Check that the discretisation structure is really uniform.
    ! More complex situations are not supported by this routine.
    if ((rsourceDiscr%ccomplexity .ne. SPDISC_UNIFORM) .and. &
        (rsourceDiscr%ccomplexity .ne. SPDISC_CONFORMAL)) then
      call output_line ('Only uniform discretisations supported!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdiscr_deriveSimpleDiscr_triquad')
      call sys_halt()
    end if

    ! Release old information if present
    call spdiscr_releaseDiscr(rdestDiscr)

    ! Check the cubature formula against the element group.
    ! This stops the program if this is not fulfilled.
    call spdiscr_checkCubature(ccubTri,ieltypTri)
    call spdiscr_checkCubature(ccubQuad,ieltypQuad)

    ! Copy the source structure to the destination.
    ! This copies all handles and hence all dynamic information
    rdestDiscr = rsourceDiscr

    ! Allocate a new element group
    allocate(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
    rdestDiscr%RelementDistr(1:rdestDiscr%inumFESpaces) = &
        rsourceDiscr%RelementDistr(1:rsourceDiscr%inumFESpaces)

    ! Loop through the element groups...
    do idistr = 1,rdestDiscr%inumFESpaces

      ! Check the element there. If it is a triangular element,
      ! change the element type to ielTypTri. If it is a quad
      ! element, change the element type to ielTypQuad.
      nve = elem_igetNVE(rsourceDiscr%RelementDistr(idistr)%celement)
      select case (nve)
      case (TRIA_NVETRI2D)
        rdestDiscr%RelementDistr(idistr)%celement = ieltypTri

        ! Init the cubature rule
        rdestDiscr%RelementDistr(idistr)%ccubTypeBilForm = ccubTri
        rdestDiscr%RelementDistr(idistr)%ccubTypeLinForm = ccubTri
        rdestDiscr%RelementDistr(idistr)%ccubTypeEval = ccubTri

        ! Get the typical transformation used with the element
        rdestDiscr%RelementDistr(idistr)%ctrafoType = &
            elem_igetTrafoType(ieltypTri)

      case (TRIA_NVEQUAD2D)
        rdestDiscr%RelementDistr(idistr)%celement = ieltypQuad

        ! Init the cubature rule
        rdestDiscr%RelementDistr(idistr)%ccubTypeBilForm = ccubQuad
        rdestDiscr%RelementDistr(idistr)%ccubTypeLinForm = ccubQuad
        rdestDiscr%RelementDistr(idistr)%ccubTypeEval = ccubQuad

        ! Get the typical transformation used with the element
        rdestDiscr%RelementDistr(idistr)%ctrafoType = &
            elem_igetTrafoType(ieltypQuad)
      end select

    end do

    ! Mark the new discretisation structure as 'copy', to prevent
    ! the dynamic information to be released.
    ! The dynamic information 'belongs' to rdiscrSource and not to the
    ! newly created rdiscrDest!
    rdestDiscr%bisCopy = .true.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_deriveDiscr_triquad_new (rsourceDiscr, ieltypTri, ieltypQuad,&
                                              rdestDiscr)

!<description>
  ! This routine derives a discretisation structure from another one.
  !
  ! rsourceDiscr is a given uniform discretisation structure. The structure
  ! will be copied to rdestDiscr and changed in such a way, that the element
  ! type and cubature formula are changed according to the parameters.
  !
  ! The new discretisation will also be a uniform discretisation based
  ! on the element celement. It is not a complete new structure, but a
  ! 'derived' structure, i.e. it uses the same dynamic information
  ! (element lists) as rsourceDiscr.
!</description>

!<input>
  ! A source discretisation structure that should be used as template
  type(t_spatialDiscretisation), intent(in), target :: rsourceDiscr

  ! The element type identifier that is to be used for all triangular
  ! elements in the new discretisation structure
  integer(I32), intent(in) :: ieltypTri

  ! The element type identifier that is to be used for all quad
  ! elements in the new discretisation structure
  integer(I32), intent(in) :: ieltypQuad
!</input>

!<inputoutput>
  ! The discretisation structure to be initialised.
  ! Any old existing discretisation information in rdestDiscr
  ! is released if necessary.
  type(t_spatialDiscretisation), intent(inout), target :: rdestDiscr
!</inputoutput>

!</subroutine>

    ! local variables
    ! INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
    ! TYPE(t_elementDistribution), POINTER :: p_relementDistr
    integer :: idistr,nve

    ! Check that the source discretisation structure is valid.
    if (rsourceDiscr%ndimension .ne. NDIM2D) then
      call output_line ('Source structure invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdiscr_deriveSimpleDiscr_triquad')
      call sys_halt()
    end if

    ! Check that the discretisation structure is really uniform.
    ! More complex situations are not supported by this routine.
    if ((rsourceDiscr%ccomplexity .ne. SPDISC_UNIFORM) .and. &
        (rsourceDiscr%ccomplexity .ne. SPDISC_CONFORMAL)) then
      call output_line ('Only uniform discretisations supported!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdiscr_deriveSimpleDiscr_triquad')
      call sys_halt()
    end if

    ! Release old information if present
    call spdiscr_releaseDiscr(rdestDiscr)

    ! Copy the source structure to the destination.
    ! This copies all handles and hence all dynamic information
    rdestDiscr = rsourceDiscr

    ! Allocate a new element group
    allocate(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
    rdestDiscr%RelementDistr(1:rdestDiscr%inumFESpaces) = &
        rsourceDiscr%RelementDistr(1:rsourceDiscr%inumFESpaces)

    ! Loop through the element groups...
    do idistr = 1,rdestDiscr%inumFESpaces

      ! Check the element there. If it is a triangular element,
      ! change the element type to ielTypTri. If it is a quad
      ! element, change the element type to ielTypQuad.
      nve = elem_igetNVE(rsourceDiscr%RelementDistr(idistr)%celement)
      select case (nve)
      case (TRIA_NVETRI2D)
        rdestDiscr%RelementDistr(idistr)%celement = ieltypTri

        ! Get the typical transformation used with the element
        rdestDiscr%RelementDistr(idistr)%ctrafoType = &
            elem_igetTrafoType(ieltypTri)

      case (TRIA_NVEQUAD2D)
        rdestDiscr%RelementDistr(idistr)%celement = ieltypQuad

        ! Get the typical transformation used with the element
        rdestDiscr%RelementDistr(idistr)%ctrafoType = &
            elem_igetTrafoType(ieltypQuad)
      end select

    end do

    ! Mark the new discretisation structure as 'copy', to prevent
    ! the dynamic information to be released.
    ! The dynamic information 'belongs' to rdiscrSource and not to the
    ! newly created rdiscrDest!
    rdestDiscr%bisCopy = .true.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_releaseDiscr (rspatialDiscr)

!<description>
  ! This routine releases a discretisation structure from memory.
!</description>

!<inputoutput>
  ! The discretisation structure to be released.
  type(t_spatialDiscretisation), intent(inout), target :: rspatialDiscr
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  type(t_elementDistribution), pointer :: p_relementDistr

  if (rspatialDiscr%ndimension .ne. 0) then
    ! Cut the connection to the other structures
    nullify(rspatialDiscr%p_rtriangulation)
    nullify(rspatialDiscr%p_rboundary)

    if (.not. rspatialDiscr%bisCopy) then

      ! Release element group lists.
      if (rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
        call storage_free (rspatialDiscr%h_IelemGroupIDs)
      end if

    else
      rspatialDiscr%h_IelemGroupIDs = ST_NOHANDLE
    end if

    ! Loop through all element groups
    do i = 1, rspatialDiscr%inumFESpaces

      p_relementDistr => rspatialDiscr%RelementDistr(i)

      ! If the element group is empty, skip it
      if (p_relementDistr%NEL .ne. 0) then

        ! Release the element list there.
        ! Take care: If the current structure is a copy of another one, the
        ! element list 'belongs' to another structure, and so we must not
        ! delete it from memory!
        if (.not. rspatialDiscr%bisCopy) then
          if (p_relementDistr%h_IelementList .ne. ST_NOHANDLE) &
            call storage_free (p_relementDistr%h_IelementList)
        else
          p_relementDistr%h_IelementList = ST_NOHANDLE
        end if

      end if

      p_relementDistr%celement = EL_UNDEFINED

    end do

    if (.not. rspatialDiscr%bisCopy) then
      if (rspatialDiscr%h_IelementCounter .ne. ST_NOHANDLE) &
        call storage_free (rspatialDiscr%h_IelementCounter)
    else
      rspatialDiscr%h_IelementCounter = ST_NOHANDLE
    end if

    ! No FE-spaces in here anymore...
    if (associated(rspatialDiscr%RelementDistr)) then
      deallocate(rspatialDiscr%RelementDistr)
    end if
    rspatialDiscr%inumFESpaces = 0

    ! Release the DOF-Mapping if necessary
    call spdiscr_releaseDofMapping(rspatialDiscr)

    ! Structure not initialised anymore
    rspatialDiscr%ndimension = 0
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_duplicateDiscrSc (rsourceDiscr, rdestDiscr, bshare)

!<description>
  ! This routine creates a copy of the discretisation structure rsourceDiscr.
  ! Depending on bshare, the destination structure rdestDiscr will either
  ! obtain a 'simple' copy (i.e. sharing all handles and all information
  ! with rsourceDiscr) or a separate copy (which costs memory for all the
  ! element information!).
!</description>

!<input>
  ! A source discretisation structure that should be used as template
  type(t_spatialDiscretisation), intent(in) :: rsourceDiscr

  ! OPTIONAL: Whether the new discretisation structure should share its information
  ! with rsourceDiscr.
  ! =FALSE: Create a complete copy of rsourceDiscr which is independent
  !  of rsourceDiscr.
  ! =TRUE: The new discretisation will not be a complete new structure, but a
  !  'derived' structure, i.e. it uses the same dynamic information
  !  (handles and therefore element lists) as rsourceDiscr.
  ! If not specified, TRUE is assumed.
  logical, intent(in), optional :: bshare
!</input>

!<inputoutput>
  ! The new discretisation structure. Any old existing information in rdestDiscr
  ! is released if necessary.
  type(t_spatialDiscretisation), intent(inout), target :: rdestDiscr
!</inputoutput>

!</subroutine>

    logical :: bshr
    integer :: i

    bshr = .true.
    if (present(bshare)) bshr = bshare

    ! Release old information if present
    call spdiscr_releaseDiscr(rdestDiscr)

    ! Copy all information
    rdestDiscr = rsourceDiscr

    ! Duplicate the element group structure
    allocate(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
    rdestDiscr%RelementDistr = rsourceDiscr%RelementDistr
    
    ! Do we have to copy data explicitly?
    if (bshr) then

      ! Mark the new discretisation structure as 'copy', to prevent
      ! the dynamic information to be released.
      ! The dynamic information 'belongs' to rdiscrSource and not to the
      ! newly created rdiscrDest!
      rdestDiscr%bisCopy = .true.
      
    else

      rdestDiscr%bisCopy = .false.

      ! element groups?
      if (rsourceDiscr%h_IelemGroupIDs .ne. ST_NOHANDLE) then
        rdestDiscr%h_IelemGroupIDs = ST_NOHANDLE
        call storage_copy(rsourceDiscr%h_IelemGroupIDs, rdestDiscr%h_IelemGroupIDs)
      end if

      ! Element counters?
      if (rsourceDiscr%h_IelementCounter .ne. ST_NOHANDLE) then
        rdestDiscr%h_IelementCounter = ST_NOHANDLE
        call storage_copy(rsourceDiscr%h_IelementCounter, rdestDiscr%h_IelementCounter)
      end if

      ! Element DOFs?
      if (rsourceDiscr%h_IelementDofs .ne. ST_NOHANDLE) then
        rdestDiscr%h_IelementDofs = ST_NOHANDLE
        call storage_copy(rsourceDiscr%h_IelementDofs, rdestDiscr%h_IelementDofs)
      end if

      ! Element DOFs indices?
      if (rsourceDiscr%h_IelementDofIdx .ne. ST_NOHANDLE) then
        rdestDiscr%h_IelementDofIdx = ST_NOHANDLE
        call storage_copy(rsourceDiscr%h_IelementDofIdx, rdestDiscr%h_IelementDofIdx)
      end if

      do i=1,rsourceDiscr%inumFESpaces
        if (rsourceDiscr%RelementDistr(i)%h_IelementList .ne. ST_NOHANDLE) then
          rdestDiscr%RelementDistr(i)%h_IelementList = ST_NOHANDLE
          call storage_copy(rsourceDiscr%RelementDistr(i)%h_IelementList,&
                            rdestDiscr%RelementDistr(i)%h_IelementList)
        end if
      end do

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_duplicateBlockDiscr (rsourceDiscr, rdestDiscr, bshare)

!<description>
  ! This routine creates a copy of the dblock iscretisation structure rsourceDiscr.
  ! Depending on bshare, the destination structure rdestDiscr will either
  ! obtain a 'simple' copy (i.e. sharing all handles and all information
  ! with rsourceDiscr) or a separate copy (which costs memory for all the
  ! element information in all the blocks!).
  !
  ! The routine does a similar job as
  ! spdiscr_deriveBlockDiscr(rsourceDiscr,rdestDiscr), but in contrast,
  ! discretisation specific information like boundary conditions are copied, too.
!</description>

!<input>
  ! A source discretisation structure that should be used as template
  type(t_blockDiscretisation), intent(in) :: rsourceDiscr

  ! OPTIONAL: Whether the new discretisation structure should share its information
  ! with rsourceDiscr.
  ! =FALSE: Create a complete copy of rsourceDiscr which is independent
  !  of rsourceDiscr.
  ! =TRUE: The new discretisation will not be a complete new structure, but a
  !  'derived' structure, i.e. it uses the same dynamic information
  !  (handles and therefore element lists) as rsourceDiscr.
  ! If not specified, TRUE is assumed.
  logical, intent(in), optional :: bshare
!</input>

!<inputoutput>
  ! The new discretisation structure. Any old existing information in rdestDiscr
  ! is released if necessary.
  type(t_blockDiscretisation), intent(inout), target :: rdestDiscr
!</inputoutput>

!</subroutine>

    logical :: bshr
    integer :: i

    bshr = .true.
    if (present(bshare)) bshr = bshare

    ! Release old information if present
    call spdiscr_releaseBlockDiscr(rdestDiscr,.true.)

    ! Check that the source discretisation structure is valid.
    if (rsourceDiscr%ndimension .le. 0) then
      call output_line ('Source structure invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_duplicateBlockDiscr')
      call sys_halt()
    end if

    ! At first, derive a new block discretisation strucutr

    ! Copy all information from the source discretisation structure containing
    ! all basic information.
    call spdiscr_deriveBlockDiscr(rsourceDiscr,rdestDiscr)

    ! Concerning the substructures, at the moment we share the information.
    ! If bshare = false, we have to create copies.
    if (.not. bshr) then
      do i=1,rsourceDiscr%ncomponents
        call spdiscr_duplicateDiscrSc (rsourceDiscr%RspatialDiscr(i), &
            rdestDiscr%RspatialDiscr(i), .false.)
      end do
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_infoBlockDiscr (rdiscr)

!<description>
    ! This subroutine outputs information about the block discretisation
!</description>

!<input>
    ! block discretisation
    type(t_blockDiscretisation), intent(in) :: rdiscr
!</input>
!</subroutine>

    ! local variables
    integer :: icomponent

    call output_lbrk()
    call output_line ('BlockDiscretisation:')
    call output_line ('--------------------')
    call output_line ('ndimension:  '//trim(sys_siL(rdiscr%ndimension,1)))
    call output_line ('ccomplexity: '//trim(sys_siL(rdiscr%ccomplexity,1)))
    call output_line ('ncomponents: '//trim(sys_siL(rdiscr%ncomponents,3)))

    if (associated(rdiscr%RspatialDiscr)) then
      do icomponent = 1, rdiscr%ncomponents
        call spdiscr_infoDiscr(rdiscr%RspatialDiscr(icomponent))
      end do
    end if

  end subroutine spdiscr_infoBlockDiscr

  ! ***************************************************************************

!<subroutine>

   subroutine spdiscr_infoDiscr (rspatialDiscr)

!<description>
     ! This subroutine outputs information about the spatial discretisation
!</description>

!<input>
     ! spatial discretisation
     type(t_spatialDiscretisation), intent(in) :: rspatialDiscr
!</input>
!</subroutine>

     ! local variable
     integer :: inumFESpace

     call output_lbrk()
     call output_line ('SpatialDiscretisation:')
     call output_line ('----------------------')
     call output_line ('ndimension:             '&
         //trim(sys_siL(rspatialDiscr%ndimension,1)))
     call output_line ('bisCopy:                '&
         //trim(sys_sl(rspatialDiscr%bisCopy)))
     call output_line ('ccomplexity:            '&
         //trim(sys_siL(rspatialDiscr%ccomplexity,1)))
     call output_line ('inumFESpaces:           '&
         //trim(sys_siL(rspatialDiscr%inumFESpaces,15)))
     call output_line ('h_IelemGroupIDs:        '&
         //trim(sys_siL(rspatialDiscr%h_IelemGroupIDs,15)))
     call output_line ('h_IelementCounter:      '&
         //trim(sys_siL(rspatialDiscr%h_IelementCounter,15)))

     if (associated(rspatialDiscr%RelementDistr)) then
       do inumFESpace = 1, rspatialDiscr%inumFESpaces
         call spdiscr_infoElementDistr(rspatialDiscr%RelementDistr(inumFESpace))
       end do
     end if

   end subroutine spdiscr_infoDiscr

   ! ***************************************************************************

!<subroutine>

   subroutine spdiscr_infoElementDistr (relementDistr)

!<description>
     ! This subroutine outputs information about the spatial discretisation
!</description>

!<input>
     ! element group
     type(t_elementDistribution), intent(in) :: relementDistr
!</input>
!</subroutine>

     call output_lbrk()
     call output_line ('ElementDistribution:')
     call output_line ('--------------------')
     call output_line ('ielement:        '//trim(sys_siL(int(relementDistr%celement),15))//&
                                      ' ('//trim(elem_getName(relementDistr%celement))//')')
     call output_line ('ccubTypeBilForm: '//trim(sys_siL(int(relementDistr%ccubTypeBilForm),15)))
     call output_line ('ccubTypeLinForm: '//trim(sys_siL(int(relementDistr%ccubTypeLinForm),15)))
     call output_line ('ccubTypeEval:    '//trim(sys_siL(int(relementDistr%ccubTypeEval),15)))
     call output_line ('ctrafoType:      '//trim(sys_siL(int(relementDistr%ctrafoType),15)))
     call output_line ('NEL:             '//trim(sys_siL(relementDistr%NEL,15)))
     call output_line ('h_IelementList:  '//trim(sys_siL(relementDistr%h_IelementList,15)))

   end subroutine spdiscr_infoElementDistr

  ! ***************************************************************************
!<subroutine>

  subroutine spdiscr_releaseDofMapping(rdiscretisation)

!<description>
  ! Releases precomputed DOF-mapping arrays.
!</description>

!<input>
  ! The discretisation structure that specifies the (scalar) discretisation.
  type(t_spatialDiscretisation), intent(inout) :: rdiscretisation
!</input>

!</subroutine>

    ! Release old arrays of necessary.
    if (rdiscretisation%bprecompiledDofMapping) then
      rdiscretisation%bprecompiledDofMapping = .false.
      call storage_free (rdiscretisation%h_IelementDofs)
      call storage_free (rdiscretisation%h_IelementDofIdx)
    end if

  end subroutine

    ! ***************************************************************************
!<function>

  elemental integer function spdiscr_igetNDofLocMax(rdiscretisation)

!<description>
  ! Calculates the maximum number of local DOF`s in this discretisation.
!</description>

!<input>
  ! The discretisation structure where information should be printed.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
!</input>

!</function>

    ! local variables
    integer :: i,imax

    imax = 0

    ! Loop through the element groups and calculate the maximum
    do i=1,rdiscretisation%inumFESpaces
      imax = max(imax,elem_igetNDofLoc(rdiscretisation%RelementDistr(i)%celement))
    end do

    spdiscr_igetNDofLocMax = imax

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine spdscr_edgeBlocking2D (rdiscretisation,redgeBlocking)

!<description>
  ! Blocks all edges in a triangulation in such a way that all edges in a block
  ! have the same FE spaces at the adjacent elements.
!</description>

  !<input>
    ! A discretisation structure that defines for all elements the element types.
    type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  !</input>

  !<output>
    ! Edge blocking structure that defines the blocking of the edges.
    type(t_edgeBlocking), intent(out) :: redgeBlocking
  !<output>
!</subroutine>

    ! local variables

    ! A list of all edges, ordered in blocks in such a way that all edges in a block
    ! have the same FE spaces at the adjacent elements.
    integer, dimension(:), pointer :: p_Iedges

    ! Array of length (#combinations+1). Defines the start positions in Iedges
    ! of every block of edges that have the same type of elements adjacent.
    ! #combinations = #element groups * (#element groups+1)/2
    integer, dimension(:), pointer :: p_IdistPositions

    ! Array of dimension(2,#combinations/2).
    ! Returns for each block described by IdistPositions the number of the
    ! element groups of the elements on each side of the edge.
    ! #combinations = #element groups * (#element groups+1)/2
    integer, dimension(:,:), pointer :: p_Idistributions

    integer, dimension(2) :: Isize

    ! Allocate memory in the structure
    redgeBlocking%nelementDistributions = rdiscretisation%inumFESpaces
    redgeBlocking%nedges = rdiscretisation%p_rtriangulation%NMT

    ! Number of combinations is n*(n+1)2 + number of element spaces.
    ! The forst coefficient is for the inner edges, the last for the
    ! boundary edges.
    redgeBlocking%nfecombinations = redgeBlocking%nelementDistributions * &
        (redgeBlocking%nelementDistributions+1) / 2 + &
        redgeBlocking%nelementDistributions

    call storage_new ('spdscr_edgeBlocking2D', 'Iedges', redgeBlocking%nedges,&
        ST_INT, redgeBlocking%h_Iedges, ST_NEWBLOCK_NOINIT)
    call storage_new ('spdscr_edgeBlocking2D', 'IedgePositions', &
        redgeBlocking%nfecombinations+1,&
        ST_INT, redgeBlocking%h_IdistPositions, ST_NEWBLOCK_NOINIT)
    Isize = (/2,redgeBlocking%nfecombinations/)
    call storage_new ('spdscr_edgeBlocking2D', 'Idistributions', &
        Isize, ST_INT, redgeBlocking%h_Idistributions, ST_NEWBLOCK_NOINIT)

    ! Calculate the blocking
    call storage_getbase_int(redgeBlocking%h_Iedges,p_Iedges)
    call storage_getbase_int(redgeBlocking%h_IdistPositions,p_IdistPositions)
    call storage_getbase_int2d(redgeBlocking%h_Idistributions,p_Idistributions)

    call spdscr_doEdgeBlocking2D (rdiscretisation,p_Iedges,p_IdistPositions,p_Idistributions)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdscr_releaseEdgeBlocking (redgeBlocking)

!<description>
  ! Releases memory used by an edge blocking structure.
!</description>

  !<inputoutput>
    ! Edge blocking structure to be cleaned up
    type(t_edgeBlocking), intent(inout) :: redgeBlocking
  !<inputoutput>
!</subroutine>

    call storage_free(redgeBlocking%h_Iedges)
    call storage_free(redgeBlocking%h_IdistPositions)
    call storage_free(redgeBlocking%h_Idistributions)
    redgeBlocking%nelementDistributions = 0
    redgeBlocking%nedges = 0
    redgeBlocking%nfecombinations = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdscr_doEdgeBlocking2D (rdiscretisation,Iedges,IdistPositions,Idistributions)

!<description>
  ! Blocks all edges in a triangulation in such a way that all edges in a block
  ! have the same FE spaces at the adjacent elements.
  ! Worker routine.
!</description>

  !<input>
    ! A discretisation structure that defines for all elements the element types.
    type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  !</input>

  !<output>
    ! A list of all edges, ordered in blocks in such a way that all edges in a block
    ! have the same FE spaces at the adjacent elements.
    integer, dimension(:), intent(out) :: Iedges

    ! Array of length (#combinations+1). Defines the start positions in Iedges
    ! of every block of edges that have the same type of elements adjacent.
    ! #combinations = #element groups * (#element groups+1)/2
    integer, dimension(:), intent(out) :: IdistPositions

    ! Array of dimension(2,#combinations/2).
    ! Returns for each block described by IdistPositions the number of the
    ! element groups of the elements on each side of the edge.
    ! Idistributions(2,:) is =0 for boundary edges.
    ! #combinations = #element groups * (#element groups+1)/2
    integer, dimension(:,:), intent(out) :: Idistributions
  !</output>

!</subroutine>

    ! local variables
    integer :: i,j,k
    type(t_triangulation), pointer :: p_rtria
    integer, dimension(:,:), allocatable :: IsortArray
    integer, dimension(:), pointer :: p_ielemGroup
    integer, dimension(:,:), pointer :: p_IelementsAtEdge

    p_rtria => rdiscretisation%p_rtriangulation

    ! Is that a uniform discretisation? If yes, we can easily fill the output
    ! arrays.
    if (rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then
      do i=1,p_rtria%NMT
        Iedges(i) = i
      end do
      IdistPositions(1) = 1
      IdistPositions(2:) = p_rtria%NMT+1
      Idistributions(1:2,1) = 1
      return
    end if

    ! Get some arrays
    call storage_getbase_int(rdiscretisation%h_IelemGroupIDs,p_ielemGroup)
    call storage_getbase_int2d(p_rtria%h_IelementsAtEdge,p_IelementsAtEdge)

    ! Generate an array that contains all edge numbers and the adjacent
    ! FE space identifiers from the discretisation.
    allocate(IsortArray(3,p_rtria%NMT))

    ! Put the edge number to coordinate 1.
    ! Put the id if the first FE space to coordinate 2.
    ! Put the id of the 2nd FE space to coordinate 3.
    do i=1,p_rtria%NMT
      IsortArray(1,i) = i
      j = p_ielemGroup(p_IelementsAtEdge(1,i))

      ! The edge may be a boundary edge.
      if (p_IelementsAtEdge(2,i) .eq. 0) then
        k = 0

        IsortArray(2,i) = j
        IsortArray(3,i) = k
      else
        k = p_ielemGroup(p_IelementsAtEdge(2,i))

        if (j .lt. k) then
          IsortArray(2,i) = j
          IsortArray(3,i) = k
        else
          IsortArray(2,i) = k
          IsortArray(3,i) = j
        end if
      end if

    end do

    ! Not sort the array -- first for the 3rd coordinate, then for the 2nd.
    call arraySort_sortByIndex_int (IsortArray, 3, SORT_STABLE)
    call arraySort_sortByIndex_int (IsortArray, 2, SORT_STABLE)

    ! Now count how many elements belong to each set.
    IdistPositions(1) = 1
    k = 1
    i = 1
    blockloop: do
      do j=i+1,p_rtria%NMT
        if ((IsortArray(2,i) .ne. IsortArray(2,j)) .or. &
            (IsortArray(3,i) .ne. IsortArray(3,j))) then
          k=k+1
          IdistPositions(k) = j
          i = j
          cycle blockloop
        end if
      end do

      ! We reached the end
      k = k + 1
      IdistPositions(k:) = p_rtria%NMT + 1
      exit
    end do blockloop

    ! Finally, collect the element numbers
    do i=1,k-1
      Idistributions(1,i) = IsortArray(2,IdistPositions(i))
      Idistributions(2,i) = IsortArray(3,IdistPositions(i))
      do j=IdistPositions(i),IdistPositions(i+1)-1
        Iedges(j) = IsortArray(1,j)
      end do
    end do
    Idistributions(:,k:) = 0

    ! That is it.
    deallocate(IsortArray)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_concatBlockDiscr (rsourceDiscr1, rsourceDiscr2, rdestDiscr,&
      rtriangulation, rboundary)

!<description>
  ! Concatenates two block discretisation structures to a new one
!</description>

!<input>
  ! The first block discretisation structure to concatenate.
  type(t_blockDiscretisation), intent(in), target :: rsourceDiscr1

  ! The second block discretisation structure to concatenate.
  ! Must have the same domain, triangulation and type FE structure
  ! as rsourceDiscr1. (So the spatial subdiscretisations are at least "conformal" to
  ! those in rsourceDiscr1).
  type(t_blockDiscretisation), intent(in), target :: rsourceDiscr2

  ! OPTIONAL: Reference to a new triangulation, the new discretisation should use.
  ! If not specified, the triangulation in rsourceDiscr1 is used.
  ! If specified, the new triangulation must be compatible to both triangulations,
  ! in rsourceDiscr1 and rsourceDiscr2.
  type(t_triangulation), intent(in), target, optional :: rtriangulation

  ! OPTIONAL: Reference to a new domain, the new discretisation should use.
  ! If not specified, the domain in rsourceDiscr1 is used.
  ! If specified, the new domain must be compatible to both domains,
  ! in rsourceDiscr1 and rsourceDiscr2.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<inputoutput>
  ! The discretisation structure to be initialised.
  ! Receives the concatenated block discretisation "rsourceDiscr1 + rsourceDiscr2".
  ! Any old existing information in rdestDiscr is released if necessary.
  ! The new discretisation shares its information with rsourceDiscr1 and
  ! rsourceDiscr2.
  type(t_blockDiscretisation), intent(inout), target, optional :: rdestDiscr
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Check that the source discretisation structure is valid.
    if ((rsourceDiscr1%ndimension .le. 0) .or. (rsourceDiscr2%ndimension .le. 0) .or.&
        (rsourceDiscr1%ndimension .ne. rsourceDiscr2%ndimension)) then
      call output_line ('Source structures invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_concatBlockDiscr')
      call sys_halt()
    end if

    ! Copy all information from the source discretisation structures

    rdestDiscr%ndimension       =  rsourceDiscr1%ndimension
    rdestDiscr%ccomplexity      =  SPDISC_UNIFORM
    if ((rsourceDiscr1%ccomplexity .eq. SPDISC_CONFORMAL) .or. &
        (rsourceDiscr2%ccomplexity .eq. SPDISC_CONFORMAL)) then
      rdestDiscr%ccomplexity      =  SPDISC_CONFORMAL
    end if
    rdestDiscr%p_rboundary      => rsourceDiscr1%p_rboundary
    rdestDiscr%p_rtriangulation => rsourceDiscr1%p_rtriangulation
    if (present(rboundary)) rdestDiscr%p_rboundary => rsourceDiscr1%p_rboundary
    if (present(rtriangulation)) rdestDiscr%p_rtriangulation => rsourceDiscr1%p_rtriangulation
    rdestDiscr%ncomponents      =  rsourceDiscr1%ncomponents + rsourceDiscr2%ncomponents

    ! Copy all substructures -- from ifirstBlock to ilastBlock.
    ! Use spdiscr_duplicateDiscrSc which savely copies the scalar discretisation
    ! structures. We set bshare=.TRUE. here, so the information is shared
    ! between the source and destination structure; the dynamic information
    ! 'belongs' to rdiscrSource and not to the newly created rdiscrDest!
    allocate(rdestDiscr%RspatialDiscr(rdestDiscr%ncomponents))
    do i = 1, rsourceDiscr1%ncomponents
      call spdiscr_duplicateDiscrSc (rsourceDiscr1%RspatialDiscr(i), &
          rdestDiscr%RspatialDiscr(i), .true.)
    end do

    do i = 1, rsourceDiscr2%ncomponents
      call spdiscr_duplicateDiscrSc (rsourceDiscr2%RspatialDiscr(i), &
          rdestDiscr%RspatialDiscr(rsourceDiscr1%ncomponents+i), .true.)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_isBlockDiscrCompatible (rdiscr1,rdiscr2,bcompatible)

!<description>
  ! Checks whether two block discretisations are compatible to each other.
  ! Two block discretisations are compatible if
  ! - they have the same dimension,
  ! - they have the same complexity,
  ! - they have the same number of components,
  ! - they have the same structure of spatial discretisations
!</description>

!<input>
  ! The first block discretisation
  type(t_blockDiscretisation), intent(in) :: rdiscr1

  ! The second block discretisation
  type(t_blockDiscretisation), intent(in) :: rdiscr2
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether the block discretisations are compatible or not.
  ! If not given, an error will inform the user if the two block
  ! discretisations are not compatible and the program will halt.
  logical, intent(out), optional :: bcompatible
!</output>

!</subroutine>

    ! local variables
    integer :: i

    ! We assume that we are compatible
    if (present(bcompatible)) bcompatible = .true.

    ! Block discretisation structures must have the same dimension,
    ! complexity and number of components
    if ((rdiscr1%ndimension .ne. rdiscr2%ndimension) .or.&
        (rdiscr1%ccomplexity .ne. rdiscr2%ccomplexity) .or.&
        (rdiscr1%ncomponents .ne. rdiscr2%ncomponents)) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Block discretisation structures are not compatible!', &
            OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_isBlockDiscrCompatible')
        call sys_halt()
      end if
    end if

    ! All spatial discretisations must be compatible.
    do i = 1, rdiscr1%ncomponents
      call spdiscr_isDiscrCompatible(rdiscr1%RspatialDiscr(i),&
          rdiscr2%RspatialDiscr(i), bcompatible)
    end do

  end subroutine spdiscr_isBlockDiscrCompatible

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_isDiscrCompatible (rdiscr1,rdiscr2,bcompatible)

!<description>
  ! Checks whether two spatial discretisations are compatible to each other.
  ! Two spatial discretisations are compatible if
  ! - they have the same dimension,
  ! - they have the same complexity,
  ! - they have the same number of different FE spaces,
  ! - they have both precomputed DOF-mapping or not,
  ! - they have the same structure of element groups
!</description>

!<input>
  ! The first spatial discretisation
  type(t_spatialDiscretisation), intent(in) :: rdiscr1

  ! The second spatial discretisation
  type(t_spatialDiscretisation), intent(in) :: rdiscr2
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether the spatial discretisations are compatible or not.
  ! If not given, an error will inform the user if the two spatial
  ! discretisations are not compatible and the program will halt.
  logical, intent(out), optional :: bcompatible
!</output>

!</subroutine>

    ! local variables
    integer :: i

    ! We assume that we are compatible
    if (present(bcompatible)) bcompatible = .true.

    ! Spatial discretisation structures must have the same dimension,
    ! complexity, number of components and DOF-mapping must be precomputed
    ! either for both discretisations or not.
    if ((rdiscr1%ndimension .ne. rdiscr2%ndimension) .or.&
        (rdiscr1%ccomplexity .ne. rdiscr2%ccomplexity) .or.&
        (rdiscr1%inumFESpaces .ne. rdiscr2%inumFESpaces) .or.&
        (rdiscr1%ndof .ne. rdiscr2%ndof) .or.&
        (rdiscr1%bprecompiledDofMapping .neqv. rdiscr2%bprecompiledDofMapping)) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Spatial discretisation structures are not compatible!', &
            OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_isDiscrCompatible')
        call sys_halt()
      end if
    end if

    ! All element groups must be compatible.
    do i = 1, rdiscr1%inumFESpaces
      call spdiscr_isElemDistrCompatible(rdiscr1%RelementDistr(i),&
          rdiscr2%RelementDistr(i), bcompatible)
    end do

  end subroutine spdiscr_isDiscrCompatible

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_isElemDistrCompatible (relemDistr1,relemDistr2,bcompatible)

!<description>
  ! Checks whether two element groups are compatible to each other.
  ! Twoelement groups are compatible if
  ! - they have the same type of element,
  ! - they have the same cubature formulas,
  ! - they have the same type of transformation,
  ! - they have the same number of elements
!</description>

!<input>
  ! The first element group
  type(t_elementDistribution), intent(in) :: relemDistr1

  ! The second spatial discretisation
  type(t_elementDistribution), intent(in) :: relemDistr2
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending
  ! on whether the element groups are compatible or not.
  ! If not given, an error will inform the user if the two element
  ! distributions are not compatible and the program will halt.
  logical, intent(out), optional :: bcompatible
!</output>

!</subroutine>

  ! We assume that we are compatible
    if (present(bcompatible)) bcompatible = .true.

    ! element groups structures must have the same dimension,
    ! complexity, number of components and DOF-mapping must be precomputed
    ! either for both discretisations or not.
    if ((relemDistr1%celement .ne. relemDistr2%celement) .or.&
        (relemDistr1%ctrafoType .ne. relemDistr1%ctrafoType) .or.&
        (relemDistr1%NEL .ne. relemDistr1%NEL)) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('element groups are not compatible!', &
            OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_isElemDistrCompatible')
        call sys_halt()
      end if
    end if

  end subroutine spdiscr_isElemDistrCompatible

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_createDefCubStructure (rdiscretisation, rcubatureInfo, &
      ccubType, nlevels, rtrafoInfo)

!<description>
  ! Creates a default cubature information structure based on a discretisation.
  ! All elements in an element set are assembled with the same cubature formula.
!</description>

!<input>
  ! A discretisation structure, the cubature information structure should be
  ! associated to.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! OPTIONAL: Cubature formula to use. If not specified, CUB_GEN_AUTO is used,
  ! which automatically choses the default cubature formula. 
  !
  ! DEPRECATED: Additionally, the following constants are allowed for compatibility to
  ! the old discretisation structure. 
  ! =CUB_GEN_DEPR_BILFORM: Take ccubTypeBilForm from the discretisation structure.
  ! =CUB_GEN_DEPR_LINFORM: Take ccubTypeLinForm from the discretisation structure.
  ! =CUB_GEN_DEPR_EVAL: Take ccubTypeEval from the discretisation structure
  integer(I32), intent(in), optional :: ccubType
  
  ! OPTIONAL: Number of refinements of the reference element if a summed
  ! cubature rule should be employed.
  integer, intent(in), optional :: nlevels

  ! OPTIONAL: A transformation structure that defines the transformation
  ! from the reference to the real element(s). If not specified, the default
  ! transformation is used.
  type(t_scalarTrafoInfo), intent(in), optional :: rtrafoInfo
!</input>

!<output>
  ! Cubature information structure to be created.
  type(t_scalarCubatureInfo), intent(out) :: rcubatureInfo
!</output>

!</subroutine>
    ! local variables
    integer :: i
    integer(I32) :: ccub

    if (.not. present(rtrafoInfo)) then
    
      ! We take as many blocks as element sets.
      rcubatureInfo%ninfoBlockCount = rdiscretisation%inumFESpaces
      allocate(rcubatureInfo%p_RinfoBlocks(rcubatureInfo%ninfoBlockCount))

      ! Loop through the element sets and insert the default cubature formula.
      do i = 1,rcubatureInfo%ninfoBlockCount
        rcubatureInfo%p_RinfoBlocks(i)%ielemGroup = i

        ! Handle all elements in the same way...
        rcubatureInfo%p_RinfoBlocks(i)%h_IelementList = ST_NOHANDLE

        rcubatureInfo%p_RinfoBlocks(i)%NEL = rdiscretisation%RelementDistr(i)%NEL

        ! Default transformation
        rcubatureInfo%p_RinfoBlocks(i)%itrafoBlock = 0
      end do
      
    else
    
      ! We take as many blocks as transformation sets
      rcubatureInfo%ninfoBlockCount = rtrafoInfo%ninfoBlockCount
      allocate(rcubatureInfo%p_RinfoBlocks(rcubatureInfo%ninfoBlockCount))

      ! Loop through the element sets and insert the default cubature formula.
      do i = 1,rcubatureInfo%ninfoBlockCount
        rcubatureInfo%p_RinfoBlocks(i)%ielemGroup = &
            rtrafoInfo%p_RinfoBlocks(i)%ielemGroup

        ! Handle all elements in the same way...
        rcubatureInfo%p_RinfoBlocks(i)%h_IelementList = &
            rtrafoInfo%p_RinfoBlocks(i)%h_IelementList

        rcubatureInfo%p_RinfoBlocks(i)%NEL = rtrafoInfo%p_RinfoBlocks(i)%NEL

        ! Id of the referring transformation block
        rcubatureInfo%p_RinfoBlocks(i)%itrafoBlock = i

      end do
      
    end if    

    ! Initialise the cubature formula
    ! Default: Automatic.
    ccub = CUB_GEN_AUTO

    ! Is the cubature formula specified?
    if (present(ccubType)) then
      select case (ccubType)
      case (CUB_GEN_DEPR_BILFORM,CUB_GEN_DEPR_LINFORM,CUB_GEN_DEPR_EVAL)
        ! Use automatic cubature formula.
        ! Nothing to do...
      case default
        ! Use predefined cubature formula.
        ccub = ccubType
      end select
    end if

    ! Initialise the cubature formula
    call spdiscr_defineCubature (rdiscretisation, rcubatureInfo, ccub, nlevels)

    ! Post-correction. Is the cubature formula specified?      
    if (present(ccubType)) then
      select case (ccubType)
      case (CUB_GEN_DEPR_BILFORM,CUB_GEN_DEPR_LINFORM,CUB_GEN_DEPR_EVAL)
        ! Fetch the DEPRECATED cubature formula from the
        ! discretisation structure.
        do i = 1,rcubatureInfo%ninfoBlockCount
          select case (ccubType)
          case (CUB_GEN_DEPR_BILFORM)
            if (rdiscretisation%RelementDistr(i)%ccubTypeBilForm .ne. CUB_UNDEFINED) then
              rcubatureInfo%p_RinfoBlocks(i)%ccubature = &
                  rdiscretisation%RelementDistr(i)%ccubTypeBilForm

#if WARN_DEPREC
              call output_line ("Using deprecated feature. Please update your code.", &
                  OU_CLASS_WARNING,OU_MODE_STD,"spdiscr_createDefCubStructure")
#endif
            end if

          case (CUB_GEN_DEPR_LINFORM)
            if (rdiscretisation%RelementDistr(i)%ccubTypeLinForm .ne. CUB_UNDEFINED) then
              rcubatureInfo%p_RinfoBlocks(i)%ccubature = &
                  rdiscretisation%RelementDistr(i)%ccubTypeLinForm

#if WARN_DEPREC
              call output_line ("Using deprecated feature. Please update your code.", &
                  OU_CLASS_WARNING,OU_MODE_STD,"spdiscr_createDefCubStructure")
#endif
            end if

          case (CUB_GEN_DEPR_EVAL)
            if (rdiscretisation%RelementDistr(i)%ccubTypeEval .ne. CUB_UNDEFINED) then
              rcubatureInfo%p_RinfoBlocks(i)%ccubature = &
                  rdiscretisation%RelementDistr(i)%ccubTypeEval

#if WARN_DEPREC
              call output_line ("Using deprecated feature. Please update your code.", &
                  OU_CLASS_WARNING,OU_MODE_STD,"spdiscr_createDefCubStructure")
#endif
            end if
          end select
        end do

      end select
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_defineCubature (rdiscretisation, rcubatureInfo, ccubType, nlevels)

!<description>
  ! Initialises the cubature formula in rcubatureInfo according to ccubType.
  ! ccubType should be a generic cubature formula like "CUB_GEN_AUTO_xxxx".
  ! The cubature formulas for all element spaces in rcubatureInfo
  ! are initialised with a cubature rule compatible to ccubType.
  ! If no compatible cubature rule exists, the cubature formula is
  ! initialised with CUB_UNDEFINED.
!</description>

!<input>
  ! A discretisation structure, the cubature information structure should be
  ! associated to.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! Generic cubature formula constant. One of the "CUB_GEN_AUTO(_xxxx)" constants.
  ! E.g., CUB_GEN_AUTO initialises the cubature formula with the default
  ! cubature rule for each FEM space.
  integer(I32), intent(in) :: ccubType

  ! OPTIONAL: Number of refinements of the reference element if a summed
  ! cubature rule should be employed.
  integer, intent(in), optional :: nlevels
!</input>

!<inputoutput>
  ! Cubature information structure to be modified.
  type(t_scalarCubatureInfo), intent(inout) :: rcubatureInfo
!</inputoutput>

!</subroutine>
    ! local variables
    integer :: i,ielemGroup

    ! Loop through the element sets and insert the default cubature formula.
    select case (cub_isExtended(ccubType))
    case (0)
      ! Not generic. Check every cubature formula if it is valid
      ! for the current element. If yes, initialise.
      do i = 1,rcubatureInfo%ninfoBlockCount
        ielemGroup = rcubatureInfo%p_RinfoBlocks(i)%ielemGroup
        call spdiscr_checkCubature (ccubType,&
            rdiscretisation%RelementDistr(ielemGroup)%celement)
        rcubatureInfo%p_RinfoBlocks(i)%ccubature = ccubType
      end do

    case (1)
      ! Generic, resolve using the element shape.
      do i = 1,rcubatureInfo%ninfoBlockCount
        ielemGroup = rcubatureInfo%p_RinfoBlocks(i)%ielemGroup
        rcubatureInfo%p_RinfoBlocks(i)%ccubature = &
            cub_resolveGenericCubType(&
              elem_igetShape(rdiscretisation%RelementDistr(ielemGroup)%celement),&
              ccubType)
      end do

    case (2)
      ! Special case
      select case (ccubType)

      case (CUB_GEN_AUTO)
        ! Standard cubature formula for that element set.
        do i = 1,rcubatureInfo%ninfoBlockCount
          ielemGroup = rcubatureInfo%p_RinfoBlocks(i)%ielemGroup
          rcubatureInfo%p_RinfoBlocks(i)%ccubature = &
              spdiscr_getStdCubature(rdiscretisation%RelementDistr(ielemGroup)%celement)
        end do

      case (CUB_GEN_AUTO_LUMPMASS)
        ! Standard cubature formula that lumps the mass matrix
        do i = 1,rcubatureInfo%ninfoBlockCount
          ielemGroup = rcubatureInfo%p_RinfoBlocks(i)%ielemGroup
          rcubatureInfo%p_RinfoBlocks(i)%ccubature = &
              spdiscr_getLumpCubature(rdiscretisation%RelementDistr(ielemGroup)%celement)
        end do

      end select

    end select

    ! Do we need a summed cubature rule?
    if (present(nlevels)) then
      do i = 1,rcubatureInfo%ninfoBlockCount
        rcubatureInfo%p_RinfoBlocks(i)%ccubature =&
            cub_getSummedCubType(rcubatureInfo%p_RinfoBlocks(i)%ccubature, nlevels)
      end do
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_releaseCubStructure (rcubatureInfo)

!<description>
  ! Cleans up an cubature information structure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_scalarCubatureInfo), intent(inout) :: rcubatureInfo
!</inputoutput>

!</subroutine>
    ! local variables
    integer :: i

    ! Loop through the element sets and release memory
    do i = 1,rcubatureInfo%ninfoBlockCount
      ! Release memory if necessary -- and if the memory belongs to us
      if ((rcubatureInfo%p_RinfoBlocks(i)%h_IelementList .ne. ST_NOHANDLE) .and.&
          rcubatureInfo%p_RinfoBlocks(i)%blocalElementList) then
        call storage_free(rcubatureInfo%p_RinfoBlocks(i)%h_IelementList)
      end if
    end do
    
    if (associated(rcubatureInfo%p_RinfoBlocks))&
        deallocate(rcubatureInfo%p_RinfoBlocks)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_getElementCubMapping (rcubatureInfo,rdiscretisation,IcubPerElement)

!<description>
  ! Determins an array which identifies for every element the corresponding
  ! cubature formula.
!</description>

!<input>
  ! Cubature information structure.
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo

  ! Associated discretisation structure.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
!</input>

!<output>
  ! Cubature identifier list. For every element i, IcubPerElement(i)
  ! identifies the cubature formula on that element.
  integer(I32), dimension(:), intent(out) :: IcubPerElement
!</output>

!</subroutine>
    ! local variables
    integer :: i,j,NEL
    integer, dimension(:), pointer :: p_IelementList
    integer(I32) :: ccubature

    ! Loop through the element sets and release memory
    do i = 1,rcubatureInfo%ninfoBlockCount

      ! For every element in the associated p_IelementList, identify
      ! its cubature formula.
      call spdiscr_getStdDiscrInfo(i,rcubatureInfo,rdiscretisation,&
          ccubature=ccubature,NEL=NEL,p_IelementList=p_IelementList)

      do j = 1,NEL
        IcubPerElement(p_IelementList(j)) = ccubature
      end do

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_getStdDiscrInfo (icubatureBlock,rcubatureInfo,rdiscretisation,&
      ielemGroup,celement,ccubature,NEL,p_IelementList,&
      rtrafoInfo,itrafoBlock,ctrafoType)

!<description>
  ! Returns typically used information for an assembly block.
!</description>

!<input>
  ! Identifier of the information block in rcubatureInfo
  ! where information should be acquired from.
  integer, intent (in) :: icubatureBlock

  ! Underlying cubature information structure.
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo

  ! OPTIONAL: Associated discretisation structure.
  type(t_spatialDiscretisation), intent(in), optional, target :: rdiscretisation

  ! OPTIONAL: Transformation information structure that defines the mapping
  ! from the reference to the real element(s). If not specified, the default
  ! transformation is chosen.
  type(t_scalarTrafoInfo), intent(in), optional :: rtrafoInfo
!</input>

!<output>
  ! Id of the underlying element group in the discretisation.
  integer, intent(out), optional :: ielemGroup

  ! Element identifier of that block
  integer(I32), intent(out), optional :: celement

  ! Cubature rule in that block
  integer(I32), intent(out), optional :: ccubature

  ! Number of elements in that block
  integer, intent(out), optional :: NEL

  ! Pointer to the list of elements of that block.
  ! =null if NEL=0.
  integer, dimension(:), pointer, optional :: p_IelementList

  ! Id of the underlying transformation block or =0
  ! if no transformation block is associated and the default transformation
  ! is to be used.
  integer, intent(out), optional :: itrafoBlock

  ! Transformation rule
  integer(I32), intent(out), optional :: ctrafoType
!</output>

!</subroutine>

    ! local variables
    integer :: ielementDistLocal,itrafoBlk
    integer(I32) :: ccub
    type(t_elementDistribution), pointer :: p_relementDistr

    ielementDistLocal = rcubatureInfo%p_RinfoBlocks(icubatureBlock)%ielemGroup
    itrafoBlk = rcubatureInfo%p_RinfoBlocks(icubatureBlock)%itrafoBlock
    ccub =  rcubatureInfo%p_RinfoBlocks(icubatureBlock)%ccubature

    ! Structure initialised?
    if (rcubatureInfo%p_RinfoBlocks(icubatureBlock)%NEL .eq. -1) then
      call output_line ("rcubatureInfo structure not initialised!", &
                        OU_CLASS_ERROR,OU_MODE_STD,"spdiscr_getStdDiscrInfo")
      call sys_halt()
    end if

    ! Return all information specified in the parameters
    if (present(ielemGroup)) then
      ielemGroup = ielementDistLocal
    end if

    if (present(rdiscretisation)) then
      p_relementDistr => rdiscretisation%RelementDistr(ielementDistLocal)

      if (present(celement)) then
        celement = p_relementDistr%celement
      end if
    else
      nullify(p_relementDistr)
      
      ! No chance
      if (present(celement)) then
        celement = EL_UNDEFINED
      end if
    end if

    if (present(ccubature)) then
      ccubature =  ccub
    end if

    if (present(NEL)) then
      NEL = rcubatureInfo%p_RinfoBlocks(icubatureBlock)%NEL
    end if

    if (present(itrafoBlock)) then
      itrafoBlock = itrafoBlk
    end if

    if (present(ctrafoType)) then
      if (present(rtrafoInfo) .and. (itrafoBlk .ne. 0)) then
        ! Transformation from the trafo structure
        ctrafoType = rtrafoInfo%p_RinfoBlocks(itrafoBlk)%ctrafoType
      else if (associated(p_relementDistr)) then
        ! Default transformation of the element
        ctrafoType = elem_igetTrafoType(p_relementDistr%celement)
      else
        ! Default transformation for that cubature formula
        ctrafoType = trafo_getDefaultTrafo(cub_igetShape(ccub))
      end if
    end if

    if (present(p_IelementList)) then
      nullify(p_IelementList)
      if (rcubatureInfo%p_RinfoBlocks(icubatureBlock)%NEL .ne. 0) then
        ! If the handle of the info block structure is not associated,
        ! take all elements of the corresponding element group.
        if (rcubatureInfo%p_RinfoBlocks(icubatureBlock)%h_IelementList .ne. ST_NOHANDLE) then
          call storage_getbase_int(&
              rcubatureInfo%p_RinfoBlocks(icubatureBlock)%h_IelementList,&
              p_IelementList)
        else if (present(rtrafoInfo) .and. (itrafoBlk .ne. 0)) then
          ! Element list from the transformation or the discretisation.
          call storage_getbase_int(&
              rtrafoInfo%p_RinfoBlocks(itrafoBlk)%h_IelementList,&
              p_IelementList)
        else
          ! Element list from the discretisation
          call storage_getbase_int(&
              rdiscretisation%RelementDistr(ielementDistLocal)%h_IelementList,&
              p_IelementList)
        end if
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_createDefTrafoStructure (rdiscretisation, rtrafoInfo)

!<description>
  ! Creates a default transformation information structure based on a discretisation.
  ! All elements in an element set are assembled with their default transformation.
!</description>

!<input>
  ! A discretisation structure, the cubature information structure should be
  ! associated to.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
!</input>

!<output>
  ! Cubature information structure to be created.
  type(t_scalarTrafoInfo), intent(out) :: rtrafoInfo
!</output>

!</subroutine>
    ! local variables
    integer :: i

    ! We take as many blocks as element sets.
    rtrafoInfo%ninfoBlockCount = rdiscretisation%inumFESpaces
    allocate(rtrafoInfo%p_RinfoBlocks(rtrafoInfo%ninfoBlockCount))

    ! Loop through the element sets and insert the default cubature formula.
    do i = 1,rtrafoInfo%ninfoBlockCount
      rtrafoInfo%p_RinfoBlocks(i)%ielemGroup = i

      ! Handle all elements in the same way...
      rtrafoInfo%p_RinfoBlocks(i)%h_IelementList = ST_NOHANDLE

      rtrafoInfo%p_RinfoBlocks(i)%NEL = rdiscretisation%RelementDistr(i)%NEL
    end do

    ! Initialise the cubature formula
    call spdiscr_defaultTrafo (rtrafoInfo,rdiscretisation)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_createDefTrafoByCubInfo (rcubatureInfo, rtrafoInfo)

!<description>
  ! Creates a default transformation information structure based on a cubature
  ! information structure. The transformation structure is structured
  ! the same way as the cubature structure and receives teh default
  ! transformation for the corresponding cells.
!</description>

!<input>
  ! Cubature information structure
  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
!</input>

!<output>
  ! Transformation structure to be created.
  type(t_scalarTrafoInfo), intent(out) :: rtrafoInfo
!</output>

!</subroutine>
    ! local variables
    integer :: i

    ! We take as many blocks as element sets.
    rtrafoInfo%ninfoBlockCount = rcubatureInfo%ninfoBlockCount
    allocate(rtrafoInfo%p_RinfoBlocks(rtrafoInfo%ninfoBlockCount))

    ! Loop through the element sets and insert the default cubature formula.
    do i = 1,rtrafoInfo%ninfoBlockCount
      rtrafoInfo%p_RinfoBlocks(i)%ielemGroup = rcubatureInfo%p_RinfoBlocks(i)%ielemGroup

      ! Handle all elements in the same way...
      rtrafoInfo%p_RinfoBlocks(i)%h_IelementList = rcubatureInfo%p_RinfoBlocks(i)%h_IelementList

      rtrafoInfo%p_RinfoBlocks(i)%NEL = rcubatureInfo%p_RinfoBlocks(i)%NEL
    end do
    
    ! Default initialisation
    call spdiscr_defaultTrafo (rtrafoInfo,rcubatureInfo=rcubatureInfo)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_defaultTrafo (rtrafoInfo,rdiscretisation,rcubatureInfo)

!<description>
  ! Initialises the transformation rule in rtrafoInfo by the default settings
  ! based on a discretisation or a cubature information structure
!</description>

!<input>
  ! OPTIONAL: A discretisation structure defining the default settings for
  ! the transformation. If this is present, it overrides the settingts
  ! of rcubatureInfo (if present as well).
  type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation

  ! OPTIONAL: Cubature information structure defining the default settings for
  ! the transformation.
  type(t_scalarCubatureInfo), intent(in), optional :: rcubatureInfo
!</input>

!<inputoutput>
  ! Transformation information structure to be modified.
  type(t_scalarTrafoInfo), intent(inout) :: rtrafoInfo
!</inputoutput>

!</subroutine>
    ! local variables
    integer :: i,ielemGroup
    integer(I32) :: ccubature

    if (present(rdiscretisation)) then

      ! Resolve using the element routines.
      do i = 1,rtrafoInfo%ninfoBlockCount
        ielemGroup = rtrafoInfo%p_RinfoBlocks(i)%ielemGroup
        rtrafoInfo%p_RinfoBlocks(i)%ctrafoType = &
            elem_igetTrafoType(rdiscretisation%RelementDistr(ielemGroup)%celement)
      end do

    else if (present(rcubatureInfo)) then

      ! Initialise the trafo formula using the default transformation
      ! valid for that cubature formula
      do i = 1,rcubatureInfo%ninfoBlockCount
        ccubature = rcubatureInfo%p_RinfoBlocks(i)%ccubature
        rtrafoInfo%p_RinfoBlocks(i)%ctrafoType = trafo_getDefaultTrafo(cub_igetShape(ccubature))
      end do

    else

      call output_line ("No default initialisation for transformation available!", &
          OU_CLASS_ERROR,OU_MODE_STD,"spdiscr_defaultTrafo")
      call sys_halt()
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_releaseTrafoStructure (rtrafoInfo)

!<description>
  ! Cleans up an transformation information structure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_scalarTrafoInfo), intent(inout) :: rtrafoInfo
!</inputoutput>

!</subroutine>
    ! local variables
    integer :: i

    ! Loop through the element sets and release memory
    do i = 1,rtrafoInfo%ninfoBlockCount
      ! Release memory if necessary -- and if the memory belongs to us
      if ((rtrafoInfo%p_RinfoBlocks(i)%h_IelementList .ne. ST_NOHANDLE) .and.&
          rtrafoInfo%p_RinfoBlocks(i)%blocalElementList) then
        call storage_free(rtrafoInfo%p_RinfoBlocks(i)%h_IelementList)
      end if
    end do

    deallocate(rtrafoInfo%p_RinfoBlocks)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_getElementTrafoMapping (rtrafoInfo,rdiscretisation,ItrafoPerElement)

!<description>
  ! Determins an array which identifies for every element the corresponding
  ! transformation formula.
!</description>

!<input>
  ! Transformation information structure.
  type(t_scalarTrafoInfo), intent(in) :: rtrafoInfo

  ! Associated discretisation structure.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
!</input>

!<output>
  ! Transformation identifier list. For every element i, ItrafoPerElement(i)
  ! identifies the transformation on that element.
  integer(I32), dimension(:), intent(out) :: ItrafoPerElement
!</output>

!</subroutine>
    ! local variables
    integer :: i,j,NEL
    integer, dimension(:), pointer :: p_IelementList
    integer(I32) :: ctrafoType

    ! Loop through the element sets and release memory
    do i = 1,rtrafoInfo%ninfoBlockCount

      ! For every element in the associated p_IelementList, identify
      ! its transformation rule.
      call spdiscr_getStdTrafoInfo(i,rtrafoInfo,rdiscretisation,&
          ctrafoType=ctrafoType,NEL=NEL,p_IelementList=p_IelementList)

      do j = 1,NEL
        ItrafoPerElement(p_IelementList(j)) = ctrafoType
      end do

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_getStdTrafoInfo (itrafoBlock,rtrafoInfo,rdiscretisation,&
      ielemGroup,celement,ctrafoType,NEL,p_IelementList)

!<description>
  ! Returns typically used information for a transformation block.
!</description>

!<input>
  ! Identifier of the information block in rtrafoInfo
  ! where information should be acquired from.
  integer, intent (in) :: itrafoBlock

  ! Underlying transformation information structure.
  type(t_scalarTrafoInfo), intent(in) :: rtrafoInfo

  ! Associated discretisation structure.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
!</input>

!<output>
  ! Id of the underlying element group in the discretisation.
  integer, intent(out), optional :: ielemGroup

  ! Element identifier of that block
  integer(I32), intent(out), optional :: celement

  ! Transformation rule in that block
  integer(I32), intent(out), optional :: ctrafoType

  ! Number of elements in that block
  integer, intent(out), optional :: NEL

  ! Pointer to the list of elements of that block.
  ! =null if NEL=0.
  integer, dimension(:), pointer, optional :: p_IelementList
!</output>

!</subroutine>

    ! local variables
    integer :: ielementDistLocal
    type(t_elementDistribution), pointer :: p_relementDistr

    ielementDistLocal = rtrafoInfo%p_RinfoBlocks(itrafoBlock)%ielemGroup
    p_relementDistr => rdiscretisation%RelementDistr(ielementDistLocal)

    ! Structure initialised?
    if (rtrafoInfo%p_RinfoBlocks(itrafoBlock)%NEL .eq. -1) then
      call output_line ("rtrafoInfo structure not initialised!", &
                        OU_CLASS_ERROR,OU_MODE_STD,"spdiscr_getStdTrafoInfo")
      call sys_halt()
    end if

    ! Return all information specified in the parameters
    if (present(ielemGroup)) then
      ielemGroup = ielementDistLocal
    end if

    if (present(celement)) then
      celement = p_relementDistr%celement
    end if

    if (present(ctrafoType)) then
      ctrafoType =  rtrafoInfo%p_RinfoBlocks(itrafoBlock)%ctrafoType
    end if

    if (present(NEL)) then
      NEL = rtrafoInfo%p_RinfoBlocks(itrafoBlock)%NEL
    end if

    if (present(p_IelementList)) then
      nullify(p_IelementList)
      if (rtrafoInfo%p_RinfoBlocks(itrafoBlock)%NEL .ne. 0) then
        ! If the handle of the info block structure is not associated,
        ! take all elements of the corresponding element group.
        if (rtrafoInfo%p_RinfoBlocks(itrafoBlock)%h_IelementList .ne. ST_NOHANDLE) then
          call storage_getbase_int(&
              rtrafoInfo%p_RinfoBlocks(itrafoBlock)%h_IelementList,&
              p_IelementList)
        else
          call storage_getbase_int(&
              rdiscretisation%RelementDistr(ielementDistLocal)%h_IelementList,&
              p_IelementList)
        end if
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_infoCubatureInfo (rcubatureInfo)

!<description>
    ! This subroutine outputs information about the cubature info structure
!</description>

!<input>
    ! cubature info structure
    type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
!</input>
!</subroutine>

    ! local variables
    integer :: iblock

    call output_lbrk()
    call output_line ('CubatureInfo:')
    call output_line ('-------------')
    call output_line ('ninfoBlockCount:  '//trim(sys_siL(rcubatureInfo%ninfoBlockCount,3)))

    if (associated(rcubatureInfo%p_RinfoBlocks)) then
      do iblock = 1, rcubatureInfo%ninfoBlockCount
        call spdiscr_infoCubatureInfoBlock(rcubatureInfo%p_RinfoBlocks(iblock))
      end do
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_infoCubatureInfoBlock (rcubatureInfoBlock)

!<description>
    ! This subroutine outputs information about the block of the
    ! cubature info structure
!</description>

!<input>
    ! block of the cubature info structure
    type(t_scalarCubatureInfoBlock), intent(in) :: rcubatureInfoBlock
!</input>
!</subroutine>

    call output_lbrk()
    call output_line ('CubatureInfoBlock:')
    call output_line ('------------------')
    call output_line ('blocalElementList: '//trim(sys_sl(rcubatureInfoBlock%blocalElementList)))
    call output_line ('ccubature:         '//trim(sys_siL(int(rcubatureInfoBlock%ccubature),15)))
    call output_line ('ielemGroup:     '//trim(sys_siL(rcubatureInfoBlock%ielemGroup,15)))
    call output_line ('itrafoBlock:       '//trim(sys_siL(rcubatureInfoBlock%itrafoBlock,15)))
    call output_line ('NEL:               '//trim(sys_siL(rcubatureInfoBlock%NEL,15)))
    call output_line ('h_IelementList:    '//trim(sys_siL(rcubatureInfoBlock%h_IelementList,15)))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_infoTrafoInfo (rtrafoInfo)

!<description>
    ! This subroutine outputs information about the trafo info structure
!</description>

!<input>
    ! trafo info structure
    type(t_scalarTrafoInfo), intent(in) :: rtrafoInfo
!</input>
!</subroutine>

    ! local variables
    integer :: iblock

    call output_lbrk()
    call output_line ('TrafoInfo:')
    call output_line ('----------')
    call output_line ('ninfoBlockCount:  '//trim(sys_siL(rtrafoInfo%ninfoBlockCount,3)))

    if (associated(rtrafoInfo%p_RinfoBlocks)) then
      do iblock = 1, rtrafoInfo%ninfoBlockCount
        call spdiscr_infoTrafoInfoBlock(rtrafoInfo%p_RinfoBlocks(iblock))
      end do
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_infoTrafoInfoBlock (rtrafoInfoBlock)

!<description>
    ! This subroutine outputs information about the block of the
    ! trafo info structure
!</description>

!<input>
    ! block of the trafo info structure
    type(t_scalarTrafoInfoBlock), intent(in) :: rtrafoInfoBlock
!</input>
!</subroutine>

    call output_lbrk()
    call output_line ('TrafoInfoBlock:')
    call output_line ('---------------')
    call output_line ('blocalElementList: '//trim(sys_sl(rtrafoInfoBlock%blocalElementList)))
    call output_line ('ctrafoType:        '//trim(sys_siL(int(rtrafoInfoBlock%ctrafoType),15)))
    call output_line ('ielemGroup:     '//trim(sys_siL(rtrafoInfoBlock%ielemGroup,15)))
    call output_line ('NEL:               '//trim(sys_siL(rtrafoInfoBlock%NEL,15)))
    call output_line ('h_IelementList:    '//trim(sys_siL(rtrafoInfoBlock%h_IelementList,15)))

  end subroutine

  ! ***************************************************************************

!<function>

  integer function spdiscr_getNelemGroups (rspatialDiscr)

!<description>
  ! Returns the number of element groups. Each element group defines
  ! a set of elements with the same element type and transformati8on.
!</description>

!<input>
  ! Discretisation structure
  type(t_spatialDiscretisation), intent(in) :: rspatialDiscr
!</input>

!<result>
  ! Number of element groups.
!</result>

!</function>

    spdiscr_getNelemGroups = rspatialDiscr%inumFESpaces

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_getElemGroupInfo (rspatialDiscr,igroup,celement,NEL,p_IelementList,ctrafoType)

!<description>
  ! Determins basic information about an element group.
!</description>

!<input>
  ! Discretisation structure
  type(t_spatialDiscretisation), intent(in) :: rspatialDiscr
  
  ! Index of the group
  integer, intent(in) :: igroup
!</input>

!<output>
  ! OPTIONAL: Element type
  integer(I32), intent(out), optional :: celement
  
  ! OPTIONAL: Number of elements in this group
  integer, intent(out), optional :: NEL
  
  ! OPTIONAL: Pointer to a list of elements beloging to this group
  integer, dimension(:), pointer, optional :: p_IelementList
  
  ! OPTIONAL: Underlying transformation
  integer(I32), intent(out), optional :: ctrafoType
!</output>

!</subroutine>

    ! Structure initialised?
    if (rspatialDiscr%inumFESpaces .le. 0) then
      call output_line ("Structure not initialised!", &
          OU_CLASS_ERROR,OU_MODE_STD,"spdiscr_getElemGroupInfo")
      call sys_halt()
    end if

    if ((igroup .le. 0) .or. (igroup .gt. rspatialDiscr%inumFESpaces)) then
      call output_line ("igroup invalid!", &
          OU_CLASS_ERROR,OU_MODE_STD,"spdiscr_getElemGroupInfo")
      call sys_halt()
    end if
    
    ! Element type in this block
    if (present(celement)) then
      celement = rspatialDiscr%RelementDistr(igroup)%celement
    end if
    
    if (present(NEL)) then
      NEL = rspatialDiscr%RelementDistr(igroup)%NEL
    end if
    
    ! List of elements
    if (present(p_IelementList)) then
      if (rspatialDiscr%RelementDistr(igroup)%h_IelementList .ne. ST_NOHANDLE) then
        call storage_getbase_int (rspatialDiscr%RelementDistr(igroup)%h_IelementList,&
            p_IelementList)
      else
        nullify(p_IelementList)
      end if
    end if
    
    if (present(ctrafoType)) then
      ctrafoType = rspatialDiscr%RelementDistr(igroup)%ctrafoType
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_getElemGroupIDs (rspatialDiscr,p_IelemGroupIDs)

!<description>
  ! Returns a pointer to a list that tells for every element the group ID,
  ! the element belongs to.
!</description>

!<input>
  ! Discretisation structure
  type(t_spatialDiscretisation), intent(in) :: rspatialDiscr
!</input>

!<output>
  ! Pointer to a list associating a group ID to each element.
  integer, dimension(:), pointer, optional :: p_IelemGroupIDs
!</output>

!</subroutine>

    ! Structure initialised?
    if (rspatialDiscr%inumFESpaces .le. 0) then
      call output_line ("Structure not initialised!", &
          OU_CLASS_ERROR,OU_MODE_STD,"spdiscr_getElemGroupIDs")
      call sys_halt()
    end if

    if (rspatialDiscr%h_IelemGroupIDs .ne. ST_NOHANDLE) then
      call storage_getbase_int(rspatialDiscr%h_IelemGroupIDs, p_IelemGroupIDs)
    else
      nullify(p_IelemGroupIDs)
    end if
    
  end subroutine
  
end module
