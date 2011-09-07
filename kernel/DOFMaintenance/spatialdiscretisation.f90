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
!#      -> Outputs information about the element distribution
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
!#      -> Checks whether two element distribution structures are compatible
!# </purpose>
!##############################################################################

module spatialdiscretisation

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
  integer(I32), parameter, public :: SPDISC_CUB_AUTOMATIC = 0

  ! Cubature formula stays unchanged.
  integer(I32), parameter, public :: SPDISC_CUB_NOCHANGE  = -1

!</constantblock>

!<constantblock description="Operator types">

  ! Mass matrix
  integer(I32), parameter, public :: SPDISC_OPTP_MASS = 0

  ! Laplace matrix
  integer(I32), parameter, public :: SPDISC_OPTP_LAPLACE = 1
  
  ! RHS
  integer(I32), parameter, public :: SPDISC_OPTP_RHS = 2
  
  ! Convection matrix
  integer(I32), parameter, public :: SPDISC_OPTP_CONVEC = 3

!</constantblock>

!</constants>

!<types>

!<typeblock>
  
  ! Element distribution structure. This structure collects for one type
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
    
    ! Cubature formula to use for the discretisation of this element pair
    ! during the evaluation of bilinear forms (matrix generation).
    ! Note: When evaluating bilinear forms, the ccubTypeBilForm
    ! constant of the test space decides about the cubature formula
    ! to be used!
    integer(I32) :: ccubTypeBilForm      = 0
    
    ! Cubature formula to use for the discretisation of this element pair
    ! during the evaluation of linear forms (RHS generation).
    integer(I32) :: ccubTypeLinForm      = 0

    ! Cubature formula to use for the evaluation of integrals over an FE
    ! function. This is used e.g. in postprocessing routines to calculate
    ! an integral to get an error to a reference solution.
    integer(I32) :: ccubTypeEval         = 0
    
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
    
    ! Handle to the element distribution identifier list.
    ! For every geometric element i, IelementDistr(i) specifies the
    ! number of the element distribution that contains that element.
    ! That way one can easily access information; e.g. retrieving the
    ! element type would be possible as follows:
    !   RelementDistr(IelementDistr(i))%itrialElement
    ! In a uniform discretisation (ccomplexity=SPDISC_UNIFORM), this
    ! handle is ST_NOHANDLE as all elements are in the
    ! element distribution 1.
    integer                          :: h_IelementDistr       = ST_NOHANDLE
    
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
    
    ! List of element distribution structures for every element type
    ! that is used in the discretisation.
    type(t_elementDistribution), dimension(:), pointer :: RelementDistr => null()
    
    ! Specifies whether the DOF-mapping is precomputed.
    logical                          :: bprecompiledDofMapping = .false.
    
    ! Number of DOF`s total. Only available if bprecompiledDofMapping=true.
    integer                          :: ndof = 0
    
    ! Specifies for every element a list of DOF`s how to map a local DOF
    ! to a global DOF. p_IelementDofIdx is a list with starting indices
    ! for every element in this list.
    ! Only available if bprecompiledDofMapping=true.
    integer :: h_IelementDofs = ST_NOHANDLE
    
    ! List of starting indices. p_IelementDofIdx(iel) apecifies the index
    ! in p_IelementDofs where the global DOF`s of element iel start.
    ! DIMENSION(nelements+1).
    ! Only available if bprecompiledDofMapping=true.
    integer :: h_IelementDofIdx = ST_NOHANDLE
  
    ! Handle to the array of coordinates of all global DOF`s. If
    ! generated, this array contains the positions of the global
    ! degrees of freedom in the order of the global DOF`s as returned
    ! by the local-to-global mapping routines
    integer :: h_DdofCoords = ST_NOHANDLE
  
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
    !   distributions, and each element distribution 
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
    integer                             :: ncomponents
    
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
    ! #combinations = #element distributions * (#element distributions+1)/2
    integer :: h_IdistPositions
    
    ! Array of dimension(2,#combinations/2).
    ! Returns for each block described by IdistPositions the number of the
    ! element distributions of the elements on each side of the edge.
    ! #combinations = #element distributions * (#element distributions+1)/2
    integer :: h_Idistributions
  end type

  public :: t_edgeBlocking

!</typeblock>

!</types>

  public :: spdiscr_initBlockDiscr
  public :: spdiscr_initDiscr_simple
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
!  ! Get from the element distribution the trial space and from that
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
      
      case DEFAULT
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
        ccubType = CUB_G3_T

      case (EL_Q0)
        ! Use Gauss 1X1
        ccubType = CUB_G1X1

      case (EL_Q1)
        ! Use trapezoidal rule
        ccubType = CUB_TRZ

      case (EL_Q1T)
        ! Use midpoint rule
        ccubType = CUB_MID
      
      case DEFAULT
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
        
      case DEFAULT
        ccubType = 0
      end select        
      
    case DEFAULT
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

      case (EL_P2)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! Use Gauss-4pt
          ccubType = CUB_VMC
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! Use Gauss-3pt
          ccubType = CUB_G3_T
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

      case (EL_Q1)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 3x3 Gauss formula
          ccubType = CUB_G3X3
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 2x2 Gauss formula
          ccubType = CUB_G2X2
        end select

      case (EL_Q2)

        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 4x4 Gauss formula
          ccubType = CUB_G4X4
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 3x3 Gauss formula
          ccubType = CUB_G3X3
        end select

      case (EL_Q3)
      
        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 4x4 Gauss formula
          ccubType = CUB_G4X4
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 3x3 Gauss formula
          ccubType = CUB_G3X3
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

      case (EL_QP1)
      
        select case (ioperation)
        case (SPDISC_OPTP_MASS)
          ! 3x3 Gauss formula
          ccubType = CUB_G3X3
        case (SPDISC_OPTP_LAPLACE,SPDISC_OPTP_RHS,SPDISC_OPTP_CONVEC)
          ! 2x2 Gauss formula
          ccubType = CUB_G2X2
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

      case default
        ccubType = 0
      end select
      
    case default
      ccubType = 0
    end select

  end function
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine spdiscr_initBlockDiscr (rblockDiscr,ncomponents,&
                                     rtriangulation, rboundary)
  
!<description>
  
  ! This routine initialises a block discretisation structure accept ncomponents
  ! solution components. Pointers to the triangulation, domain and boundary
  ! conditions are saved in the structure.
  !
  ! The routine performs only basic initialisation. The caller must
  ! separately initialise the the specific scalar discretisation structures 
  ! of each solution component (as collected in the RspatialDiscr
  ! array of the rblockDiscr structure).
  
!</description>

!<input>
  
  ! The triangulation structure underlying to the discretisation.
  type(t_triangulation), intent(in), target    :: rtriangulation
  
  ! Number of solution components maintained by the block structure
  integer, intent(in), optional                :: ncomponents
  
  ! OPTIONAL: The underlying domain.
  type(t_boundary), intent(in), target, optional :: rboundary

!</input>
  
!<output>
  
  ! The block discretisation structure to be initialised.
  type(t_blockDiscretisation), intent(out) :: rblockDiscr
  
!</output>
  
!</subroutine>

  ! Initialise the variables of the structure for the simple discretisation
  rblockDiscr%ndimension       = rtriangulation%ndim
  rblockDiscr%ccomplexity      = SPDISC_UNIFORM
  rblockDiscr%p_rtriangulation => rtriangulation
  if (present(rboundary)) then
    rblockDiscr%p_rboundary    => rboundary
  else
    nullify(rblockDiscr%p_rboundary)
  end if

  rblockDiscr%ncomponents      = ncomponents
  allocate(rblockDiscr%RspatialDiscr(ncomponents))

  ! That is it.  
  
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
    call spdiscr_initBlockDiscr (rblockDiscr, 1,&
        rspatialDiscr%p_rtriangulation, rspatialDiscr%p_rboundary)
    
    ! Put a copy of the spatial discretisation to first component
    ! of the the block discretisation. Share the data.
    call spdiscr_duplicateDiscrSc (rspatialDiscr,&
        rblockDiscr%RspatialDiscr(1), .true.)
  
  end subroutine  

  ! ***************************************************************************
  
!<subroutine>

  subroutine spdiscr_deriveBlockDiscr (rsourceDiscr, rdestDiscr, &
      ifirstBlock, ilastBlock, rtriangulation, rboundary)
  
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
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  ! Any old existing information in rdestDiscr is released if necessary.
  type(t_blockDiscretisation), intent(inout), target :: rdestDiscr
!</output>
  
!</subroutine>

    ! local variables
    integer :: ifirst, ilast, ncount, i
    
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
                                     rdestDiscr%RspatialDiscr(i), .true.)
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

  subroutine spdiscr_initDiscr_simple (rspatialDiscr,celement, ccubType,&
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
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  integer(I32), intent(in) :: ccubType
  
  ! The triangulation structure underlying to the discretisation.
  type(t_triangulation), intent(in), target :: rtriangulation
  
  ! The underlying domain.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  type(t_spatialDiscretisation), intent(inout), target :: rspatialDiscr
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  integer, dimension(:), pointer :: p_Iarray
  type(t_elementDistribution), pointer :: p_relementDistr
  integer(I32) :: ccub

  ! Automatically determine cubature formula if necessary  
  ccub = ccubType
  if (ccub .eq. SPDISC_CUB_AUTOMATIC) &
      ccub = spdiscr_getStdCubature(celement)
  
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
  rspatialDiscr%h_IelementDistr = ST_NOHANDLE
  
  ! Initialise the first element distribution
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
  
  ! Check the cubature formula against the element distribution.
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

  subroutine spdiscr_initDiscr_triquad (rspatialDiscr, ieltyptri, ieltypquad,&
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
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  integer(I32), intent(in) :: ccubTypeTri

  ! Cubature formula CUB_xxxx to use for calculating integrals on 
  ! quadrilateral elements
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  integer(I32), intent(in) :: ccubTypeQuad
  
  ! The triangulation structure underlying to the discretisation.
  type(t_triangulation), intent(in), target :: rtriangulation
  
  ! The underlying domain.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  type(t_spatialDiscretisation), intent(inout), target :: rspatialDiscr
!</output>
  
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
  if (ccubTri .eq. SPDISC_CUB_AUTOMATIC) &
      ccubTri = spdiscr_getStdCubature(ieltypTri)
  ccubQuad = ccubTypeQuad
  if (ccubQuad .eq. SPDISC_CUB_AUTOMATIC) &
      ccubQuad = spdiscr_getStdCubature(ieltypQuad)
  
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
  
  ! Allocate an array containing the element distribution for each element
  call storage_new ('spdiscr_initDiscr_triquad', 'h_ItrialElements', &
      rtriangulation%NEL, ST_INT, rspatialDiscr%h_IelementDistr,   &
      ST_NEWBLOCK_NOINIT)
  call storage_getbase_int (rspatialDiscr%h_IelementDistr,p_Iarray)

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
        ! Triangular elements are in element distribution 1
        p_Iarray(i) = 1
        
        ! This is the IelemCount(1)-th triangle
        IelemCount(1) = IelemCount(1)+1
        p_IelementCounter(i) = IelemCount(1)
      else
        ! Quad elements are in element distribution 2
        p_Iarray(i) = 2

        ! This is the IelemCount(2)-th quad
        IelemCount(2) = IelemCount(2)+1
        p_IelementCounter(i) = IelemCount(2)
      end if
    end do
  else
    ! Pure triangular mesh
    do i = 1, rtriangulation%NEL
      ! Triangular elements are in element distribution 1
      p_Iarray(i) = 1
      
      ! This is the IelemCount(1)-th triangle
      IelemCount(1) = IelemCount(1)+1
      p_IelementCounter(i) = IelemCount(1)
    end do
  end if
  
  ! Initialise the first element distribution
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
  
  ! Check the cubature formula against the element distribution.
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
  ! element distribution. j counts how many elements we found
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

  subroutine spdiscr_deriveSimpleDiscrSc (rsourceDiscr, celement, ccubType, &
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
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  ! A value SPDISC_CUB_NOCHANGE means:
  ! take the cubature formula from the source discretisation.
  integer(I32), intent(in) :: ccubType
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  ! Any old existing discretisation information in rdestDiscr
  ! is released if necessary.
  type(t_spatialDiscretisation), intent(inout), target :: rdestDiscr
!</output>
  
!</subroutine>

  ! local variables
  ! INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
  ! TYPE(t_elementDistribution), POINTER :: p_relementDistr
  integer(I32) :: ccub

  ! Automatically determine cubature formula if necessary  
  ccub = ccubType
  if (ccub .eq. SPDISC_CUB_AUTOMATIC) &
      ccub = spdiscr_getStdCubature(celement)
  
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
    ! Check the cubature formula against the element distribution.
    ! This stops the program if this is not fulfilled.
    call spdiscr_checkCubature(ccub,celement)
  end if
  
  ! Copy the source structure to the destination.
  ! This copies all handles and hence all dynamic information
  rdestDiscr = rsourceDiscr
  
  ! Allocate a new element distribution and copy content from source
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

  subroutine spdiscr_deriveDiscr_triquad (rsourceDiscr, ieltypTri, ieltypQuad,&
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
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  integer(I32), intent(in) :: ccubTypeTri

  ! Cubature formula to use for calculating integrals on quad
  ! elements in the new discretisation structure.
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  integer(I32), intent(in) :: ccubTypeQuad
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  ! Any old existing discretisation information in rdestDiscr
  ! is released if necessary.
  type(t_spatialDiscretisation), intent(inout), target :: rdestDiscr
!</output>
  
!</subroutine>

    ! local variables
    ! INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
    ! TYPE(t_elementDistribution), POINTER :: p_relementDistr
    integer(I32) :: ccubTri,ccubQuad
    integer :: idistr,nve

    ! Automatically determine cubature formula if necessary  
    ccubTri = ccubTypeTri
    if (ccubTri .eq. SPDISC_CUB_AUTOMATIC) &
        ccubTri = spdiscr_getStdCubature(ieltypTri)
    ccubQuad = ccubTypeQuad
    if (ccubQuad .eq. SPDISC_CUB_AUTOMATIC) &
        ccubQuad = spdiscr_getStdCubature(ieltypQuad)
    
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
    
    ! Check the cubature formula against the element distribution.
    ! This stops the program if this is not fulfilled.
    call spdiscr_checkCubature(ccubTri,ieltypTri)
    call spdiscr_checkCubature(ccubQuad,ieltypQuad)
    
    ! Copy the source structure to the destination.
    ! This copies all handles and hence all dynamic information
    rdestDiscr = rsourceDiscr
    
    ! Allocate a new element distribution
    allocate(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
    rdestDiscr%RelementDistr(1:rdestDiscr%inumFESpaces) = &
        rsourceDiscr%RelementDistr(1:rsourceDiscr%inumFESpaces)

    ! Loop through the element distributions...
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
    
      ! Release element distribution lists.
      if (rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
        call storage_free (rspatialDiscr%h_IelementDistr)
      end if
      
    else
      rspatialDiscr%h_IelementDistr = ST_NOHANDLE
    end if
    
    ! Loop through all element distributions
    do i = 1, rspatialDiscr%inumFESpaces
    
      p_relementDistr => rspatialDiscr%RelementDistr(i)
      
      ! If the element distribution is empty, skip it
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
    deallocate(rspatialDiscr%RelementDistr)
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
  
!<output>
  ! The new discretisation structure. Any old existing information in rdestDiscr
  ! is released if necessary.
  type(t_spatialDiscretisation), intent(inout), target :: rdestDiscr
!</output>
  
!</subroutine>

    logical :: bshr
    
    bshr = .true.
    if (present(bshare)) bshr = bshare

    ! Release old information if present
    call spdiscr_releaseDiscr(rdestDiscr)

    ! Currently, this routine supports only bshare=TRUE!
    if (bshr) then
    
      ! Copy all information
      rdestDiscr = rsourceDiscr
      
      ! Duplicate the element distribution structure
      allocate(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
      rdestDiscr%RelementDistr = rsourceDiscr%RelementDistr
    
      ! Mark the new discretisation structure as 'copy', to prevent
      ! the dynamic information to be released.
      ! The dynamic information 'belongs' to rdiscrSource and not to the
      ! newly created rdiscrDest!
      rdestDiscr%bisCopy = .true.
      
    else
      call output_line ('bshare=FALSE currently not supported!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_duplicateDiscrSc')  
      call sys_halt()
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
  
!<output>
  ! The new discretisation structure. Any old existing information in rdestDiscr
  ! is released if necessary.
  type(t_blockDiscretisation), intent(inout), target :: rdestDiscr
!</output>
  
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
     call output_line ('h_IelementDistr:        '&
         //trim(sys_siL(rspatialDiscr%h_IelementDistr,15)))
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
     ! element distribution
     type(t_elementDistribution), intent(in) :: relementDistr
!</input>
!</subroutine>

     call output_lbrk()
     call output_line ('ElementDistribution:')
     call output_line ('--------------------')
     call output_line ('ielement:        '//trim(sys_siL(int(relementDistr%celement),15)))
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

    ! Release coordinate arrays
    if (rdiscretisation%h_DdofCoords .ne. ST_NOHANDLE)&
        call storage_free (rdiscretisation%h_DdofCoords)

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
    
    ! Loop through the element distributions and calculate the maximum
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
    ! #combinations = #element distributions * (#element distributions+1)/2
    integer, dimension(:), pointer :: p_IdistPositions
    
    ! Array of dimension(2,#combinations/2).
    ! Returns for each block described by IdistPositions the number of the
    ! element distributions of the elements on each side of the edge.
    ! #combinations = #element distributions * (#element distributions+1)/2
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
    ! #combinations = #element distributions * (#element distributions+1)/2
    integer, dimension(:), intent(out) :: IdistPositions
    
    ! Array of dimension(2,#combinations/2).
    ! Returns for each block described by IdistPositions the number of the
    ! element distributions of the elements on each side of the edge.
    ! Idistributions(2,:) is =0 for boundary edges.
    ! #combinations = #element distributions * (#element distributions+1)/2
    integer, dimension(:,:), intent(out) :: Idistributions
  !</output>
  
!</subroutine>

    ! local variables
    integer :: i,j,k
    type(t_triangulation), pointer :: p_rtria
    integer, dimension(:,:), allocatable :: IsortArray
    integer, dimension(:), pointer :: p_IelementDistr
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
    call storage_getbase_int(rdiscretisation%h_IelementDistr,p_IelementDistr)
    call storage_getbase_int2d(p_rtria%h_IelementsAtEdge,p_IelementsAtEdge)
    
    ! Generate an array that contains all edge numbers and the adjacent
    ! FE space identifiers from the discretisation.
    allocate(IsortArray(3,p_rtria%NMT))
    
    ! Put the edge number to coordinate 1.
    ! Put the id if the first FE space to coordinate 2.
    ! Put the íd of the 2nd FE space to coordinate 3.
    do i=1,p_rtria%NMT
      IsortArray(1,i) = i
      j = p_IelementDistr(p_IelementsAtEdge(1,i))
      
      ! The edge may be a boundary edge.
      if (p_IelementsAtEdge(2,i) .eq. 0) then
        k = 0

        IsortArray(2,i) = j
        IsortArray(3,i) = k
      else
        k = p_IelementDistr(p_IelementsAtEdge(2,i))

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
  
!<output>
  ! The discretisation structure to be initialised.
  ! Receives the concatenated block discretisation "rsourceDiscr1 + rsourceDiscr2".
  ! Any old existing information in rdestDiscr is released if necessary.
  ! The new discretisation shares its information with rsourceDiscr1 and
  ! rsourceDiscr2.
  type(t_blockDiscretisation), intent(inout), target, optional :: rdestDiscr
!</output>
  
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
  ! - they have the same structure of element distributions
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

    ! All element distributions must be compatible.
    do i = 1, rdiscr1%inumFESpaces
      call spdiscr_isElemDistrCompatible(rdiscr1%RelementDistr(i),&
          rdiscr2%RelementDistr(i), bcompatible)
    end do

  end subroutine spdiscr_isDiscrCompatible

  ! ***************************************************************************
  
!<subroutine>

  subroutine spdiscr_isElemDistrCompatible (relemDistr1,relemDistr2,bcompatible)

!<description>
  ! Checks whether two element distributions are compatible to each other.
  ! Twoelement distributions are compatible if
  ! - they have the same type of element,
  ! - they have the same cubature formulas,
  ! - they have the same type of transformation,
  ! - they have the same number of elements
!</description>
  
!<input>
  ! The first element distribution
  type(t_elementDistribution), intent(in) :: relemDistr1
  
  ! The second spatial discretisation
  type(t_elementDistribution), intent(in) :: relemDistr2
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending
  ! on whether the element distributions are compatible or not.
  ! If not given, an error will inform the user if the two element
  ! distributions are not compatible and the program will halt.
  logical, intent(out), optional :: bcompatible
!</output>

!</subroutine>

  ! We assume that we are compatible
    if (present(bcompatible)) bcompatible = .true.
    
    ! Element distributions structures must have the same dimension,
    ! complexity, number of components and DOF-mapping must be precomputed
    ! either for both discretisations or not.
    if ((relemDistr1%celement .ne. relemDistr2%celement) .or.&
        (relemDistr1%ccubTypeBilForm .ne. relemDistr1%ccubTypeBilForm) .or.&
        (relemDistr1%ccubTypeLinForm .ne. relemDistr1%ccubTypeLinForm) .or.&
        (relemDistr1%ccubTypeEval .ne. relemDistr1%ccubTypeEval) .or.&
        (relemDistr1%ctrafoType .ne. relemDistr1%ctrafoType) .or.&
        (relemDistr1%NEL .ne. relemDistr1%NEL)) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Element distributions are not compatible!', &
            OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_isElemDistrCompatible')
        call sys_halt()
      end if
    end if
    
  end subroutine spdiscr_isElemDistrCompatible

end module
