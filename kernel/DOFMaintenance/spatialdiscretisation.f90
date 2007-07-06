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
!# 1.) spdiscr_initBlockDiscr2D
!#     -> Initialises a 2D block discretisation structure for the
!#        discretisation of multiple equations
!#
!# 2.) spdiscr_initDiscr_simple
!#     -> Initialise a scalar discretisation structure. One element type for
!#        all geometric elements, for test and trial functions.
!#
!# 3.) spdiscr_initDiscr_combined
!#     -> Initialise a scalar discretisation structure. One element type for
!#        all geometric elements in the trial space, another element type
!#        for all geometric elements in the test space
!#
!# 4.) spdiscr_deriveBlockDiscr
!#     -> Creates a block discretisation structure as a subset of
!#        another block discretisation structure
!#
!# 5.) spdiscr_deriveSimpleDiscrSc
!#     -> Based on an existing discretisation structure, derive a new
!#        discretisation structure with different trial elements
!#
!# 6.) spdiscr_releaseDiscr
!#     -> Release a scalar discretisation structure.
!#
!# 7.) spdiscr_releaseBlockDiscr
!#     -> Releases a block discretisation structure from memory
!#
!# 8.) spdiscr_checkCubature
!#     -> Checks if a cubature formula is compatible to an element
!#        distribution.
!#
!# </purpose>
!##############################################################################

MODULE spatialdiscretisation

  USE fsystem
  USE storage
  USE triangulation
  USE boundary
  USE boundarycondition
  USE transformation
  USE element
  USE cubature
  
  IMPLICIT NONE
  
!<constants>

!<constantblock description="Constants defining the complexity of the discretisation">

  ! Uniform discretisation: All elements are of the same type.
  INTEGER, PARAMETER :: SPDISC_UNIFORM   = 0
  
  ! Conformal discretisation: Elements of different FE spaces are mixed,
  ! but the DOF's 'fit together' (e.g. quads/tri elements with same DOF's,
  ! isoparametric elements on the boundary).
  INTEGER, PARAMETER :: SPDISC_CONFORMAL = 1
  
  ! Mixed discretisation: Elements of different FE spaces, arbitrary mixed.
  INTEGER, PARAMETER :: SPDISC_MIXED     = 2

!</constantblock>

!<constantblock>

  ! Maximum number of different FE spaces mixed in one discretisation.
  INTEGER, PARAMETER :: SPDISC_MAXFESPACES = 8

  ! Maximum number of solution components (i.e. scalar equation in the 
  ! PDE) that are supported simultaneously.
  INTEGER, PARAMETER :: SPDISC_MAXEQUATIONS = 16
  
!</constantblock>

!</constants>

!<types>

!<typeblock>
  
  ! Element distribution structure. This structure collects for one type
  ! of element (e.g. $Q_1$), on which geometric element primitives it is
  ! to be used. In the t_spatialDiscretisation there is a list of these
  ! element structures for each type of element. This e.g. allows to use
  ! triangular elements with $P_1$, quad elements with $Q_1$ and possible
  ! isoparametric "surface" elements to be mixed in the triangulation.
  !
  ! The structure is assigned to a triangulation by means of the 'parent'
  ! structure t_spatialDiscretisation, which contains a pointer to it.
  
  TYPE t_elementDistribution
  
    ! Element identifier for trial functions to use in this element list
    ! during the evaluation of bilinear forms (matrix generation).
    INTEGER(I32) :: itrialElement        = EL_UNDEFINED
    
    ! Element identifier for test functions to use in this element list
    ! during the evaluation of linear and bilinear forms.
    INTEGER(I32) :: itestElement         = EL_UNDEFINED
    
    ! Cubature formula to use for the discretisation of this element pair
    ! during the evaluation of bilinear forms (matrix generation).
    INTEGER :: ccubTypeBilForm      = 0
    
    ! Cubature formula to use for the discretisation of this element pair
    ! during the evaluation of linear forms (RHS generation).
    INTEGER :: ccubTypeLinForm      = 0

    ! Cubature formula to use for the evaluation of integrals over an FE
    ! function. This is used e.g. in postprocessing routines to calculate
    ! an integral to get an error to a reference solution.
    INTEGER :: ccubTypeEval         = 0
    
    ! Type of transformation to use from the reference element to
    ! the real element. One of the TRAFO_IDxxxx constants of the module 
    ! 'transformation' identifying the type of transformation.
    ! The same transformation is used for both, the trial and the test
    ! space, during the evaluation of linear as well as bilinear forms
    ! (matrix/RHS generation).
    INTEGER(I32) :: ctrafoType      = TRAFO_ID_UNKNOWN
    
    ! Number of elements in the list p_IelementList.
    ! May vary from the actual length of p_IelementList!
    INTEGER(PREC_ELEMENTIDX) :: NEL = 0
    
    ! Handle to list of element numbers that are discretised with this 
    ! combination of trial/test functions
    INTEGER :: h_IelementList       = ST_NOHANDLE

  END TYPE
  
!</typeblock>
  
!<typeblock>
  
  ! The central discretisation structure corresponding to one mesh level.
  ! Here, all information about the discretisation are collected (mesh
  ! information, trial functions, test functions,...).
  !
  ! Remark: The structure realises only the discretisation of 'a scalar
  !  equation'. For multidimensional problems, there are multiple of
  !  these structures, each one for one PDE. I this case, the structure
  !  is part of the block discretisation structure below and
  !  'hung into' each scalar matrix that discretises that equation.
  
  TYPE t_spatialDiscretisation
  
    ! Dimension of the discretisation. 0=not initialised, 
    ! 2=2D discretisation, 3=3D discretisation
    INTEGER                          :: ndimension             = 0
    
    ! Whether the discretisation structure is a copy of another discretisation
    ! structure. If set to TRUE, the structure was derived from another one
    ! and shares the same dynamic information (element lists,...).
    ! (This prevents the release-routine from releasing memory when
    ! cleaning up the structure.)
    LOGICAL                          :: bisCopy                = .FALSE.
  
    ! Pointer to the domain that is discretised
    TYPE(t_boundary), POINTER        :: p_rboundary            => NULL()
    
    ! Pointer to the underlying triangulation of the mesh (2D)
    TYPE(t_triangulation), POINTER   :: p_rtriangulation       => NULL()

    ! Complexity of the discretisation. One of the SPDISC_xxxx constants
    INTEGER                          :: ccomplexity            = SPDISC_UNIFORM
    
    ! Flag: Trial and Test functions in all element distributions
    ! are the same. If FALSE, there's at öeast one element distribution
    ! with different trial and test functions.
    LOGICAL                          :: bidenticalTrialAndTest = .TRUE.
    
    ! Handle to trial function identifier list: For every geometric element, 
    ! an element identifier about which element to use.
    ! In a uniform discretisation (ccomplexity=SPDISC_UNIFORM), this
    ! handle is ST_NOHANDLE as all elements are of the same type.
    INTEGER                          :: h_ItrialElements       = ST_NOHANDLE

    ! Handle to test function identifier list: For every geometric element, 
    ! an element identifier about which element to use.
    ! Coincides with p_DtrialElements if bidenticalTrialAndTest=true!
    ! In a uniform discretisation (ccomplexity=SPDISC_UNIFORM), this
    ! handle is ST_NOHANDLE as all elements are of the same type.
    INTEGER                          :: h_ItestElements        = ST_NOHANDLE
    
    ! Handle to an 'element counter' array. For every element of every
    ! type, there is a unique running number given to that element in the
    ! corresponding element subset.
    !
    ! Example: Think of a mixed mesh of triangles and quads, discretised
    !  with $P_1$ and $Q_1$. Then there are two disjunct element sets,
    !  one with triangles, one with quads. Every triangle gets a running
    !  number (1,2,3,...) and every quad gets a running number (1,2,3,...).
    !  These numbers are stored here, corresponding to each element.
    !
    ! Note that the handle may be =ST_NOHANDLE. Whether it's set up or not
    ! depends on the discretisation. The information is usually used in
    ! the DOFMapping-routines to compute the mapping between local and
    ! global degrees of freedom.
    INTEGER                          :: h_IelementCounter      = ST_NOHANDLE
    
    ! Number of different FE spaces mixed in this discretisation.
    ! This is the number of elements occupied in RelementDisttribution.
    ! Must be <= SPDISC_MAXFESPACES!
    INTEGER                          :: inumFESpaces           = 0
    
    ! List of element distribution structures for every element type
    ! that is used in the discretisation.
    TYPE(t_elementDistribution), DIMENSION(SPDISC_MAXFESPACES) :: RelementDistribution
    
  END TYPE
  
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
  
  TYPE t_blockDiscretisation
  
    ! Dimension of the discretisation. 0=not initialised, 
    ! 2=2D discretisation, 3=3D discretisation
    INTEGER                          :: ndimension             = 0

    ! Complexity of the discretisation. One of the SPDISC_xxxx constants.
    ! SPDISC_UNIFORM = all elements in each discretisation
    !   substructure RspatialDiscretisation(:) are the same.
    ! SPDISC_CONFORMAL = Elements of different FE spaces are mixed,
    !   but the DOF's 'fit together'. Each discretisation substructure 
    !   RspatialDiscretisation(:) has exactly the same number of element
    !   distributions, and each element distribution 
    !     RspatialDiscretisation(1)%Relementistributions(i), 
    !     RspatialDiscretisation(2)%Relementistributions(i),
    !     RspatialDiscretisation(3)%Relementistributions(i),...
    !   describe exactly the same set of elements (Same size, same type,
    !   same order in the element lists,...).
    INTEGER                          :: ccomplexity            = SPDISC_UNIFORM
  
    ! Pointer to the domain that is discretised
    TYPE(t_boundary), POINTER        :: p_rboundary            => NULL()
    
    ! Pointer to the underlying triangulation of the mesh (2D)
    TYPE(t_triangulation), POINTER   :: p_rtriangulation     => NULL()

    ! Pointer to the analytical description of the boundary conditions
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions => NULL()

    ! Number of solution components maintained by this structure.
    INTEGER                             :: ncomponents
    
    ! A list of up to ncomponents scalar spatial discretisation structures.
    ! Each structure corresponds to one solution component and defines
    ! trial-/test-functions, complexity of the discretisation etc.
    TYPE(t_spatialDiscretisation), &
      DIMENSION(SPDISC_MAXEQUATIONS)    :: RspatialDiscretisation
    
  END TYPE
  
!</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_checkCubature (ccubType,ielementType)
  
!<description>
  
  ! This routine checks if the cubature formula of type icubType can be applied
  ! to the elements of the type ielementType.
  ! If this is not possible, an error is thrown.
  
!</description>

!<input>
  ! The cubature formula to be tested
  INTEGER, INTENT(IN)                       :: ccubType
  
  ! The element type the cubature formula should be checked against
   INTEGER(I32), INTENT(IN)                       :: ielementType
!</input>
  
!</subroutine>

  INTEGER :: NVE
  LOGICAL :: bcompatible

  ! Get from the element distribution the trial space and from that
  ! the number of vertices, the element expects.
  NVE = elem_igetNVE(ielementType)
  
  bcompatible = .TRUE.
  
  ! Now we directly access the cubature constants in cubature.f90!
  ! This is the only point in the kernel where this is necessary.
  
  ! 1D?
  IF (ccubType .LE. 99) bcompatible = .FALSE.

  ! 3D?
  IF (ccubType .GE. 300) bcompatible = .FALSE.
  
  ! Quad?
  IF ((ccubType .GE. 200) .AND. (ccubType .LE. 249)) THEN
    ! Tri?
    IF (NVE .EQ. 3) bcompatible = .FALSE.
  END IF
  
  ! Tri?
  IF ((ccubType .GE. 250) .AND. (ccubType .LE. 299)) THEN
    ! Quad?
    IF (NVE .EQ. 4) bcompatible = .FALSE.
  END IF
  
  IF (.NOT. bcompatible) THEN
    PRINT *,'Element and cubature formula not compatible!'
    CALL sys_halt()
  END IF
  
  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION spdiscr_getLumpCubature (ielementType,ndim) RESULT (ccubType)
  
!<description>
  ! this routine tries to determine a cubature formula identifier according
  ! to a given element type, such that the corresponding mass matrix will
  ! get diagonal (mass lumping). If this is not possible, 0 is returned.
!</description>

!<input>
  ! An element type identifier
  INTEGER(I32), INTENT(IN)                       :: ielementType
  
  ! OPTIONAL: Dimension identifier. NDIM2D=2D, NDIM3D=3D. If not specified,
  ! 2D is assumed.
  INTEGER, INTENT(IN), OPTIONAL             :: ndim
!</input>

!<result>
  ! A cubature formula identifier that will diagonalise the mass matrix,
  ! or 0 if such an identifier is unknown / not possible.
!</result>
  
!</function>

    IF (PRESENT(ndim)) THEN
      IF (ndim .NE. NDIM2D) THEN
        PRINT *,'spdiscr_getLumpCubature: Only 2D supported.'
        CALL sys_halt()
      END IF
    END IF

    SELECT CASE (elem_getPrimaryElement(ielementType))
    CASE (EL_P0)
      ! Use Gauss 1X1
      ccubType = CUB_G1_T

    CASE (EL_P1)
      ! Use trapezoidal rule
      ccubType = CUB_TRZ_T

    CASE (EL_Q0)
      ! Use Gauss 1X1
      ccubType = CUB_G1X1

    CASE (EL_Q1)
      ! Use trapezoidal rule
      ccubType = CUB_TRZ

    CASE (EL_Q1T)
      ! Use midpoint rule
      ccubType = CUB_MID
      
    CASE DEFAULT
      ccubType = 0
    END SELECT

  END FUNCTION
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_initBlockDiscr2D (rblockDiscr,ncomponents,&
                                       rtriangulation, rboundary, rboundaryConditions)
  
!<description>
  
  ! This routine initialises a block discretisation structure accept ncomponents
  ! solution components. Pointers to the triangulation, domain and boundary
  ! conditions are saved in the structure.
  !
  ! The routine performs only basic initialisation. The caller must
  ! separately initialise the the specific scalar discretisation structures 
  ! of each solution component (as collected in the RspatialDiscretisation
  ! array of the rblockDiscr structure).
  
!</description>

!<input>
  
  ! The triangulation structure underlying to the discretisation.
  TYPE(t_triangulation), INTENT(IN), TARGET    :: rtriangulation
  
  ! The underlying domain.
  TYPE(t_boundary), INTENT(IN), TARGET         :: rboundary
  
  ! Number of solution components maintained by the block structure
  INTEGER, INTENT(IN)                          :: ncomponents
  
  ! OPTIONAL: The analytical description of the boundary conditions.
  ! Parameter can be ommitted if boundary conditions are not defined.
  TYPE(t_boundaryConditions), TARGET, OPTIONAL :: rboundaryConditions
  
!</input>
  
!<output>
  
  ! The block discretisation structure to be initialised.
  TYPE(t_blockDiscretisation), INTENT(OUT) :: rblockDiscr
  
!</output>
  
!</subroutine>

  ! Initialise the variables of the structure for the simple discretisation
  rblockDiscr%ndimension             = NDIM2D
  rblockDiscr%ccomplexity            = SPDISC_UNIFORM
  rblockDiscr%p_rtriangulation       => rtriangulation
  rblockDiscr%p_rboundary            => rboundary
  IF (PRESENT(rboundaryConditions)) THEN
    rblockDiscr%p_rboundaryConditions  => rboundaryConditions
  ELSE
    NULLIFY(rblockDiscr%p_rboundaryConditions)
  END IF
  rblockDiscr%ncomponents            = ncomponents

  ! That's it.  
  
  END SUBROUTINE  

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_deriveBlockDiscr (rsourceDiscr, rdestDiscr, &
      ifirstBlock,ilastBlock)
  
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
  TYPE(t_blockDiscretisation), INTENT(IN), TARGET :: rsourceDiscr

  ! OPTIONAL: Number of the block in rsourceDiscr that should be
  ! used as first block in rdestDiscr. Default value is =1.
  INTEGER, INTENT(IN), OPTIONAL :: ifirstBlock

  ! OPTIONAL: Number of the last block in rsourceDiscr that should be
  ! used as last block in rdestDiscr. Default value is the 
  ! number of components in rsourceDiscr.
  INTEGER, INTENT(IN), OPTIONAL :: ilastBlock
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  TYPE(t_blockDiscretisation), INTENT(INOUT), TARGET :: rdestDiscr
!</output>
  
!</subroutine>

    ! local variables
    INTEGER :: ifirst, ilast, ncount
    
    ! Check that the source discretisation structure is valid.
    IF (rsourceDiscr%ndimension .LE. 0) THEN
      PRINT *,'spdiscr_deriveBlockDiscr: Source structure invalid!'
      CALL sys_halt()
    END IF

    ! Evaluate the optional parameters
    ifirst = 1
    ilast = rsourceDiscr%ncomponents
    
    IF (PRESENT(ifirstBlock)) THEN
      ifirst = MIN(MAX(ifirst,ifirstBlock),ilast)
    END IF
    
    IF (PRESENT(ilastBlock)) THEN
      ilast = MAX(MIN(ilast,ilastBlock),ifirst)
    END IF
    
    ncount = ilast-ifirst+1
    
    ! Copy all information from the source discretisation structure

    rdestDiscr%ndimension             =  rsourceDiscr%ndimension           
    rdestDiscr%ccomplexity            =  rsourceDiscr%ccomplexity          
    rdestDiscr%p_rboundary            => rsourceDiscr%p_rboundary          
    rdestDiscr%p_rtriangulation       => rsourceDiscr%p_rtriangulation     
    rdestDiscr%p_rboundaryConditions  => rsourceDiscr%p_rboundaryConditions
    rdestDiscr%ncomponents            =  rsourceDiscr%ncomponents          

    ! Copy all substructures -- from ifirstBlock to ilastBlock
    rdestDiscr%RspatialDiscretisation(1:ncount) = &
      rsourceDiscr%RspatialDiscretisation(ifirstBlock:ilastBlock)
      
    ! Mark the scalar discretisation structures as 'being a copy of something'.
    ! This is to prevent the dynamic information to be released.
    ! The dynamic information 'belongs' to rdiscrSource and not to the
    ! newly created rdiscrDest!
    rdestDiscr%RspatialDiscretisation(1:ncount)%bisCopy = .TRUE.

    END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_releaseBlockDiscr (rblockDiscr, breleaseSubstruc)
  
!<description>
  ! This routine releases a block discretisation structure from memory.
!</description>

!<input>
  ! Release substructures.
  ! If set to TRUE, the memory of all scalar spatial discretisation structures
  !   in rblockDiscr is also released from memory.
  ! Is set to FALSE, only rblockDiscr is cleaned up, the substructures
  !   are ignored.
  LOGICAL, INTENT(IN) :: breleaseSubstruc
!</input>

!<inputoutput>
  ! The block discretisation structures to be released.
  TYPE(t_blockDiscretisation), INTENT(INOUT), TARGET :: rblockDiscr
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i

  ! Cut the connection to the other structures
  NULLIFY(rblockDiscr%p_rtriangulation)
  NULLIFY(rblockDiscr%p_rboundary)
  NULLIFY(rblockDiscr%p_rboundaryConditions)
  
  ! Release substructures?
  IF (breleaseSubstruc) THEN
    DO i=1,rblockDiscr%ncomponents
      CALL spdiscr_releaseDiscr (rblockDiscr%RspatialDiscretisation(i))
    END DO
  END IF
  rblockDiscr%ncomponents = 0

  ! Structure not initialised anymore
  rblockDiscr%ndimension = 0
  
  END SUBROUTINE  

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_initDiscr_simple (rspatialDiscr,ieltyp, ccubType,&
                                       rtriangulation, rboundary, rboundaryConditions)
  
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
  INTEGER(I32), INTENT(IN)                       :: ieltyp
  
  ! Cubature formula to use for calculating integrals
  INTEGER, INTENT(IN)                       :: ccubType
  
  ! The triangulation structure underlying to the discretisation.
  TYPE(t_triangulation), INTENT(IN), TARGET    :: rtriangulation
  
  ! The underlying domain.
  TYPE(t_boundary), INTENT(IN), TARGET         :: rboundary
  
  ! OPTIONAL: The analytical description of the boundary conditions.
  ! Parameter can be ommitted if boundary conditions are not defined.
  TYPE(t_boundaryConditions), TARGET, OPTIONAL :: rboundaryConditions
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  TYPE(t_spatialDiscretisation), INTENT(INOUT), TARGET :: rspatialDiscr
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
  TYPE(t_elementDistribution), POINTER :: p_relementDistr
  
  ! Do we have a structure?
  IF (rspatialDiscr%ndimension .NE. 0) THEN
    ! Release the old structure.
    CALL spdiscr_releaseDiscr(rspatialDiscr)
  END IF

  ! Initialise the variables of the structure for the simple discretisation
  rspatialDiscr%ndimension             = NDIM2D
  rspatialDiscr%p_rtriangulation       => rtriangulation
  rspatialDiscr%p_rboundary            => rboundary
  rspatialDiscr%ccomplexity            = SPDISC_UNIFORM
  
  ! All trial elements are ieltyp:
  
!  CALL storage_new1D ('spdiscr_initDiscr_simple', 'h_ItrialElements', &
!        rtriangulation%NEL, ST_INT, rspatialDiscr%h_ItrialElements,   &
!        ST_NEWBLOCK_NOINIT)
!  CALL storage_getbase_int (rspatialDiscr%h_ItrialElements,p_Iarray)
!  DO i=1,rtriangulation%NEL
!    p_Iarray(i) = ieltyp
!  END DO
  rspatialDiscr%h_ItrialElements = ST_NOHANDLE
  
  ! All test elements are ieltyp.
  ! Use the same handle for trial and test functions to save memory!
  rspatialDiscr%bidenticalTrialAndTest = .TRUE.
  rspatialDiscr%h_ItestElements = rspatialDiscr%h_ItrialElements  
  
  ! Initialise the first element distribution
  rspatialDiscr%inumFESpaces           = 1
  p_relementDistr => rspatialDiscr%RelementDistribution(1)
  
  ! Initialise test and trial space for that block
  p_relementDistr%itrialElement = ieltyp
  p_relementDistr%itestElement = ieltyp
  p_relementDistr%ccubTypeBilForm = ccubType
  p_relementDistr%ccubTypeLinForm = ccubType
  p_relementDistr%ccubTypeEval = ccubType
  
  ! Get the typical transformation used with the element
  p_relementDistr%ctrafoType = elem_igetTrafoType(ieltyp)
  
  ! Check the cubature formula against the element distribution.
  ! This stops the program if this is not fulfilled.
  CALL spdiscr_checkCubature(ccubType,ieltyp)

  ! Initialise an 'identity' array containing the numbers of all elements.
  ! This list defines the sequence how elements are processed, e.g. in the
  ! assembly of matrices/vectors.
  CALL storage_new1D ('spdiscr_initDiscr_simple', 'h_IelementList', &
        rtriangulation%NEL, ST_INT, p_relementDistr%h_IelementList,   &
        ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (p_relementDistr%h_IelementList,p_Iarray)
  DO i=1,rtriangulation%NEL
    p_Iarray(i) = i
  END DO
  
  ! Save the number of elements in that element list.
  p_relementDistr%NEL = rtriangulation%NEL
  
  ! This is a complete new structure, everything 'belongs' to this.
  rspatialDiscr%bisCopy = .FALSE.
  
  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_initDiscr_triquad (rspatialDiscr,&
      ieltyptri,ieltypquad,ccubTypeTri,ccubTypeQuad,&
      rtriangulation, rboundary, rboundaryConditions)
  
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
  INTEGER(I32), INTENT(IN)                       :: ieltypTri

  ! The element type identifier that is to be used for all quadrilateral elements.
  INTEGER(I32), INTENT(IN)                       :: ieltypQuad
  
  ! Cubature formula to use for calculating integrals on triangular elements
  INTEGER, INTENT(IN)                       :: ccubTypeTri

  ! Cubature formula to use for calculating integrals on quadrilateral elements
  INTEGER, INTENT(IN)                       :: ccubTypeQuad
  
  ! The triangulation structure underlying to the discretisation.
  TYPE(t_triangulation), INTENT(IN), TARGET    :: rtriangulation
  
  ! The underlying domain.
  TYPE(t_boundary), INTENT(IN), TARGET         :: rboundary
  
  ! OPTIONAL: The analytical description of the boundary conditions.
  ! Parameter can be ommitted if boundary conditions are not defined.
  TYPE(t_boundaryConditions), TARGET, OPTIONAL :: rboundaryConditions
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  TYPE(t_spatialDiscretisation), INTENT(INOUT), TARGET :: rspatialDiscr
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i,j
  INTEGER, DIMENSION(2) :: IelemCount
  INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray,p_IelementCounter
  TYPE(t_elementDistribution), POINTER :: p_relementDistrTria,p_relementDistrQuad
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
  
  ! Do we have a structure?
  IF (rspatialDiscr%ndimension .NE. 0) THEN
    ! Release the old structure.
    CALL spdiscr_releaseDiscr(rspatialDiscr)
  END IF

  ! Initialise the variables of the structure for the simple discretisation
  rspatialDiscr%ndimension             = NDIM2D
  rspatialDiscr%p_rtriangulation       => rtriangulation
  rspatialDiscr%p_rboundary            => rboundary
  rspatialDiscr%ccomplexity            = SPDISC_CONFORMAL
  
  ! Allocate an array containing the element type for each element
  CALL storage_new1D ('spdiscr_initDiscr_triquad', 'h_ItrialElements', &
        rtriangulation%NEL, ST_INT, rspatialDiscr%h_ItrialElements,   &
        ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (rspatialDiscr%h_ItrialElements,p_Iarray)

  ! Allocate an array with an element counter for every element type.
  CALL storage_new1D ('spdiscr_initDiscr_triquad', 'h_IelementCounter', &
        rtriangulation%NEL, ST_INT, rspatialDiscr%h_IelementCounter,   &
        ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (rspatialDiscr%h_IelementCounter,p_IelementCounter)
  
  ! Create both arrays simultaneously.
  CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
  
  IelemCount(:) = 0
  IF (UBOUND(p_IverticesAtElement,1) .GE. 4) THEN
    ! There are quads and probably triangles in the mesh
    DO i=1,rtriangulation%NEL
      IF (p_IverticesAtElement (4,i) .EQ. 0) THEN
        ! Triangular element
        p_Iarray(i) = ieltypTri
        
        ! This is the IelemCount(1)'th triangle
        IelemCount(1) = IelemCount(1)+1
        p_IelementCounter(i) = IelemCount(1)
      ELSE
        ! Quad element
        p_Iarray(i) = ieltypQuad

        ! This is the IelemCount(2)'th quad
        IelemCount(2) = IelemCount(2)+1
        p_IelementCounter(i) = IelemCount(2)
      END IF
    END DO
  ELSE
    ! Pure triangular mesh
    DO i=1,rtriangulation%NEL
      ! Triangular element
      p_Iarray(i) = ieltypTri
      
      ! This is the IelemCount(1)'th triangle
      IelemCount(1) = IelemCount(1)+1
      p_IelementCounter(i) = IelemCount(1)
    END DO
  END IF
  
  ! Trial and test element coincide.
  ! Use the same handle for trial and test functions to save memory!
  rspatialDiscr%bidenticalTrialAndTest = .TRUE.
  rspatialDiscr%h_ItestElements = rspatialDiscr%h_ItrialElements  
  
  ! Initialise the first element distribution
  rspatialDiscr%inumFESpaces           = 2
  p_relementDistrTria => rspatialDiscr%RelementDistribution(1)
  p_relementDistrQuad => rspatialDiscr%RelementDistribution(2)
  
  ! Initialise test and trial space for that block
  p_relementDistrTria%itrialElement = ieltypTri
  p_relementDistrTria%itestElement = ieltypTri
  p_relementDistrTria%ccubTypeBilForm = ccubTypeTri
  p_relementDistrTria%ccubTypeLinForm = ccubTypeTri
  p_relementDistrTria%ccubTypeEval = ccubTypeTri

  p_relementDistrQuad%itrialElement = ieltypQuad
  p_relementDistrQuad%itestElement = ieltypQuad
  p_relementDistrQuad%ccubTypeBilForm = ccubTypeQuad
  p_relementDistrQuad%ccubTypeLinForm = ccubTypeQuad
  p_relementDistrQuad%ccubTypeEval = ccubTypeQuad
  
  ! Get the typical transformation used with the element
  p_relementDistrTria%ctrafoType = elem_igetTrafoType(ieltypTri)
  p_relementDistrQuad%ctrafoType = elem_igetTrafoType(ieltypQuad)
  
  ! Check the cubature formula against the element distribution.
  ! This stops the program if this is not fulfilled.
  CALL spdiscr_checkCubature(ccubTypeTri,ieltypTri)
  CALL spdiscr_checkCubature(ccubTypeQuad,ieltypQuad)

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
  IF (rtriangulation%InelOfType(TRIA_NVETRI2D) .NE. 0) THEN
    
    CALL storage_new1D ('spdiscr_initDiscr_triquad', 'h_IelementList', &
          rtriangulation%InelOfType(TRIA_NVETRI2D), &
          ST_INT, p_relementDistrTria%h_IelementList,   &
          ST_NEWBLOCK_NOINIT)
          
    CALL storage_getbase_int (p_relementDistrTria%h_IelementList,p_Iarray)
    
    IF (UBOUND(p_IverticesAtElement,1) .GE. 4) THEN
      ! There are quads and probably triangles in the mesh
      DO i=1,rtriangulation%NEL
        IF (p_IverticesAtElement(4,i) .EQ. 0) THEN
          j = j+1
          p_Iarray(j) = i
        END IF
      END DO
    ELSE
      ! Pure triangular mesh
      DO i=1,rtriangulation%NEL
        j = j+1
        p_Iarray(j) = i
      END DO
    END IF
  END IF
  
  ! Collect all quads
  j = 0
  IF (rtriangulation%InelOfType(TRIA_NVEQUAD2D) .NE. 0) THEN
    
    CALL storage_new1D ('spdiscr_initDiscr_triquad', 'h_IelementList', &
          rtriangulation%InelOfType(TRIA_NVEQUAD2D), &
          ST_INT, p_relementDistrQuad%h_IelementList,   &
          ST_NEWBLOCK_NOINIT)
          
    CALL storage_getbase_int (p_relementDistrQuad%h_IelementList,p_Iarray)
    
    ! Because of the IF above, there are for sure quads in the mesh!
    DO i=1,rtriangulation%NEL
      IF (p_IverticesAtElement(4,i) .NE. 0) THEN
        j = j+1
        p_Iarray(j) = i
      END IF
    END DO
    
  END IF
  
  rspatialDiscr%bisCopy = .FALSE.
  
  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_deriveSimpleDiscrSc (rsourceDiscr, ieltyp, ccubType, &
                                          rdestDiscr)
  
!<description>
  ! This routine derives a discretisation structure from another one.
  !
  ! rsourceDiscr is a given uniform discretisation structure. The structure
  ! will be copied to rdestDiscr and changed in such a way, that the element
  ! type and cubature formula are changed according to the parameters.
  !
  ! The new discretisation will also be a uniform discretisation based
  ! on the element ieltyp. It's not a complete new structure, but a
  ! 'derived' structure, i.e. it uses the same dynamic information
  ! (element lists) as rsourceDiscr.
!</description>

!<input>
  ! A source discretisation structure that should be used as template
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rsourceDiscr

  ! The element type identifier that is to be used for all elements
  ! in the new discretisation structure
  INTEGER(I32), INTENT(IN)                       :: ieltyp
  
  ! Cubature formula to use for calculating integrals
  ! in the new discretisation structure
  INTEGER, INTENT(IN)                       :: ccubType
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  TYPE(t_spatialDiscretisation), INTENT(INOUT), TARGET :: rdestDiscr
!</output>
  
!</subroutine>

  ! local variables
  ! INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
  ! TYPE(t_elementDistribution), POINTER :: p_relementDistr
  
  ! Check that the source discretisation structure is valid.
  IF (rsourceDiscr%ndimension .LE. 0) THEN
    PRINT *,'spdiscr_deriveSimpleDiscr: Source structure invalid!'
    CALL sys_halt()
  END IF
  
  ! Check that the discretisation structure is really uniform.
  ! More complex situations are not supported by this routine.
  IF (rsourceDiscr%ccomplexity .NE. SPDISC_UNIFORM) THEN
    PRINT *,'spdiscr_deriveSimpleDiscr only supports uniform discretisations!'
    CALL sys_halt()
  END IF
  
  ! Check the cubature formula against the element distribution.
  ! This stops the program if this is not fulfilled.
  CALL spdiscr_checkCubature(ccubType,ieltyp)
  
  ! Copy the source structure to the destination.
  ! This copies all handles and hence all dynamic information
  rdestDiscr = rsourceDiscr
  
  ! Change the element type of all trial functions to ieltyp
  rdestDiscr%RelementDistribution(1)%itrialElement = ieltyp
  
  IF (rdestDiscr%bidenticalTrialAndTest) THEN
    rdestDiscr%RelementDistribution(1)%itestElement = ieltyp
  END IF
  rdestDiscr%bidenticalTrialAndTest = &
    rdestDiscr%RelementDistribution(1)%itrialElement .EQ. &
    rdestDiscr%RelementDistribution(1)%itestElement
  
  ! Init the cubature rule
  rdestDiscr%RelementDistribution(1)%ccubTypeBilForm = ccubType
  rdestDiscr%RelementDistribution(1)%ccubTypeLinForm = ccubType
  rdestDiscr%RelementDistribution(1)%ccubTypeEval = ccubType
  
  ! Get the typical transformation used with the element
  rdestDiscr%RelementDistribution(1)%ctrafoType = elem_igetTrafoType(ieltyp)
  
  ! Mark the new discretisation structure as 'copy', to prevent
  ! the dynamic information to be released.
  ! The dynamic information 'belongs' to rdiscrSource and not to the
  ! newly created rdiscrDest!
  rdestDiscr%bisCopy = .TRUE.
  
  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_initDiscr_combined (rspatialDiscr, &
                               ieltypTrial, ieltypTest, ccubType,&
                               rtriangulation, rboundary, rboundaryConditions)
  
!<description>
  ! This routine initialises a discretisation structure for a uniform
  ! discretisation with one element in the trial space and a probably
  ! different element in the test space.
  !
  ! If rspatialDiscr is NULL(), a new structure will be created. Otherwise,
  ! the existing structure is recreated/updated.
!</description>

!<input>
  ! The element type identifier that is to be used for all elements
  ! in the trial space.
  INTEGER(I32), INTENT(IN)                       :: ieltypTrial
  
  ! The element type identifier that is to be used for all elements
  ! in the test space.
  INTEGER(I32), INTENT(IN)                       :: ieltypTest

  ! Cubature formula to use for calculating integrals
  INTEGER, INTENT(IN)                       :: ccubType
  
  ! The triangulation structure underlying to the discretisation.
  TYPE(t_triangulation), INTENT(IN), TARGET    :: rtriangulation
  
  ! The underlying domain.
  TYPE(t_boundary), INTENT(IN), TARGET         :: rboundary
  
  ! OPTIONAL: The analytical description of the boundary conditions.
  ! Parameter can be ommitted if boundary conditions are not defined.
  TYPE(t_boundaryConditions), TARGET, OPTIONAL :: rboundaryConditions
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  TYPE(t_spatialDiscretisation), INTENT(INOUT), TARGET :: rspatialDiscr
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i,ctrafoTest
  INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
  TYPE(t_elementDistribution), POINTER :: p_relementDistr
  
  ! Do we have a structure?
  IF (rspatialDiscr%ndimension .NE. 0) THEN
    ! Release the old structure.
    CALL spdiscr_releaseDiscr(rspatialDiscr)
  END IF

  ! Initialise the variables of the structure for the simple discretisation
  rspatialDiscr%ndimension             = NDIM2D
  rspatialDiscr%p_rtriangulation       => rtriangulation
  rspatialDiscr%p_rboundary            => rboundary
  rspatialDiscr%ccomplexity            = SPDISC_UNIFORM
  
  rspatialDiscr%bidenticalTrialAndTest = ieltypTrial .EQ. ieltypTest
  
!  ! All trial elements are ieltypTrial:
!  CALL storage_new1D ('spdiscr_initDiscr_combined', 'h_ItrialElements', &
!        rtriangulation%NEL, ST_INT, rspatialDiscr%h_ItrialElements,   &
!        ST_NEWBLOCK_NOINIT)
!  CALL storage_getbase_int (rspatialDiscr%h_ItrialElements,p_Iarray)
!  DO i=1,rtriangulation%NEL
!    p_Iarray(i) = ieltypTrial
!  END DO
!
!  ! All test elements are ieltypTest:
!  IF (.NOT. rspatialDiscr%bidenticalTrialAndTest)
!    ! All test elements are ieltypTest.
!    CALL storage_new1D ('spdiscr_initDiscr_combined', 'h_ItestElements', &
!          rtriangulation%NEL, ST_INT, rspatialDiscr%h_ItestElements,   &
!          ST_NEWBLOCK_NOINIT)
!    CALL storage_getbase_int (rspatialDiscr%h_ItestElements,p_Iarray)
!    DO i=1,rtriangulation%NEL
!      p_Iarray(i) = ieltypTest
!    END DO
!  ELSE
!    rspatialDiscr%h_ItestElements = rspatialDiscr%h_ItrialElements
!  END IF
  rspatialDiscr%h_ItrialElements = ST_NOHANDLE
  rspatialDiscr%h_ItestElements = ST_NOHANDLE

  ! Initialise the first element distribution
  rspatialDiscr%inumFESpaces = 1
  p_relementDistr => rspatialDiscr%RelementDistribution(1)
  
  ! Initialise test and trial space for that block
  p_relementDistr%itrialElement = ieltypTrial
  p_relementDistr%itestElement = ieltypTest
  p_relementDistr%ccubTypeBilForm = ccubType
  p_relementDistr%ccubTypeLinForm = ccubType
  p_relementDistr%ccubTypeEval = ccubType
  
  ! Get the typical transformation used with the element
  p_relementDistr%ctrafoType = elem_igetTrafoType(ieltypTrial)
  ctrafoTest = elem_igetTrafoType(ieltypTrial)
  
  IF (p_relementDistr%ctrafoType .NE. ctrafoTest) THEN
    PRINT *,'spdiscr_initDiscr_combined: Elements incompatible due to different'
    PRINT *,'transformation between reference and real element!'
    CALL sys_halt()
  END IF
  
  ! Check the cubature formula against the element distribution.
  ! This stops the program if this is not fulfilled.
  CALL spdiscr_checkCubature(ccubType,ieltypTrial)
  CALL spdiscr_checkCubature(ccubType,ieltypTest)

  ! Save the number of elements in that element list.
  p_relementDistr%NEL = rtriangulation%NEL

  ! Initialise an 'identity' array containing the numbers of all elements.
  ! This list defines the sequence how elements are processed, e.g. in the
  ! assembly of matrices/vectors.
  CALL storage_new1D ('spdiscr_initDiscr_combined', 'h_IelementList', &
        rtriangulation%NEL, ST_INT, p_relementDistr%h_IelementList,   &
        ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (p_relementDistr%h_IelementList,p_Iarray)
  DO i=1,rtriangulation%NEL
    p_Iarray(i) = i
  END DO
  
  ! This is a complete new structure, everything 'belongs' to this.
  rspatialDiscr%bisCopy = .FALSE.

  END SUBROUTINE  

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_releaseDiscr (rspatialDiscr)
  
!<description>
  ! This routine releases a discretisation structure from memory.
!</description>

!<inputoutput>
  ! The discretisation structure to be released.
  TYPE(t_spatialDiscretisation), INTENT(INOUT), TARGET :: rspatialDiscr
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  TYPE(t_elementDistribution), POINTER :: p_relementDistr

  ! Cut the connection to the other structures
  NULLIFY(rspatialDiscr%p_rtriangulation)
  NULLIFY(rspatialDiscr%p_rboundary)
  
  ! Release element identifier lists.
  ! The element identifier list is never a copy of another structure!
  IF (rspatialDiscr%ccomplexity .NE. SPDISC_UNIFORM) THEN
    ! The handles may coincide, so release them only once!
    IF (rspatialDiscr%h_ItestElements .NE. rspatialDiscr%h_ItrialElements) THEN
      CALL storage_free (rspatialDiscr%h_ItestElements)
    ELSE
      rspatialDiscr%h_ItestElements = ST_NOHANDLE
    END IF
    CALL storage_free (rspatialDiscr%h_ItrialElements)
  END IF
  
  ! Loop through all element distributions
  DO i=1,rspatialDiscr%inumFESpaces
  
    p_relementDistr => rspatialDiscr%RelementDistribution(i)
    
    ! Release the element list there.
    ! Take care: If the current structure is a copy of another one, the
    ! element list 'belongs' to another structure, and so we mustn't
    ! delete it from memory!
    IF (.NOT. rspatialDiscr%bisCopy) THEN
      IF (p_relementDistr%h_IelementList .NE. ST_NOHANDLE) &
        CALL storage_free (p_relementDistr%h_IelementList)
    ELSE
      p_relementDistr%h_IelementList = ST_NOHANDLE
    END IF
    
    p_relementDistr%itrialElement = EL_UNDEFINED
    p_relementDistr%itestElement  = EL_UNDEFINED
    
  END DO
  
  IF (.NOT. rspatialDiscr%bisCopy) THEN
    IF (rspatialDiscr%h_IelementCounter .NE. ST_NOHANDLE) &
      CALL storage_free (rspatialDiscr%h_IelementCounter)
  ELSE    
    rspatialDiscr%h_IelementCounter = ST_NOHANDLE
  END IF
  
  ! No FE-spaces in here anymore...
  rspatialDiscr%inumFESpaces = 0
  
  ! Structure not initialised anymore
  rspatialDiscr%ndimension = 0
  
  END SUBROUTINE  

END MODULE
