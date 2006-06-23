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
!# 2.) spdiscr_releaseBlockDiscr
!#     -> Releases a block discretisation structure from memory
!#
!# 3.) spdiscr_initDiscr_simple
!#     -> Initialise a scalar discretisation structure. One element type for
!#        all geometric elements, for test and trial functions.
!#
!# 4.) spdiscr_initDiscr_combined
!#     -> Initialise a scalar discretisation structure. One element type for
!#        all geometric elements in the trial space, another element type
!#        for all geometric elements in the test space
!#
!# 5.) spdiscr_releaseDiscr
!#     -> Release a scalar discretisation structure.
!#
!# 6.) spdiscr_checkCubature
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
    INTEGER :: itrialElement        = EL_UNDEFINED
    
    ! Element identifier for test functions to use in this element list
    ! during the evaluation of linear and bilinear forms.
    INTEGER :: itestElement         = EL_UNDEFINED
    
    ! Cubature formula to use for the discretisation of this element pair
    ! during the evaluation of bilinear forms.
    INTEGER :: ccubType             = 0
    
    ! Cubature formula to use for the discretisation of this element pair
    ! during the evaluation of linear forms (RHS generation).
    INTEGER :: ccubTypeLin          = 0
    
    ! Type of transformation to use from the reference element to
    ! the real element. One of the TRAFO_IDxxxx constants of the module 
    ! 'transformation' identifying the type of transformation.
    ! The same transformation is used for both, the trial and the test
    ! space, during the evaluation of linear as well as bilinear forms
    ! (matrix/RHS generation).
    INTEGER :: ctrafoType           = TRAFO_IDUNKNOWN
    
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
  !  these structures, each one for one PDE.  I this case, the structure
  !  is part of the block discretisation structure below and
  !  'hung into' each scalar matrix that discretises that equation.
  
  TYPE t_spatialDiscretisation
  
    ! Dimension of the discretisation. 0=not initialised, 
    ! 2=2D discretisation, 3=3D discretisation
    INTEGER                          :: ndimension             = 0
  
    ! Pointer to the domain that is discretised
    TYPE(t_boundary), POINTER        :: p_rdomain              => NULL()
    
    ! Pointer to the underlying triangulation of the mesh (2D)
    TYPE(t_triangulation), POINTER   :: p_rtriangulation     => NULL()

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
  
    ! Pointer to the domain that is discretised
    TYPE(t_boundary), POINTER        :: p_rdomain              => NULL()
    
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
  INTEGER, INTENT(IN)                       :: ielementType
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
    STOP
  END IF
  
  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_initBlockDiscr2D (rblockDiscr,ncomponents,&
                                       rtriangulation, rdomain, rboundaryConditions)
  
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
  TYPE(t_boundary), INTENT(IN), TARGET         :: rdomain
  
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
  rblockDiscr%p_rtriangulation       => rtriangulation
  rblockDiscr%p_rdomain              => rdomain
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
  NULLIFY(rblockDiscr%p_rdomain)
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
                                       rtriangulation, rdomain, rboundaryConditions)
  
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
  INTEGER, INTENT(IN)                       :: ieltyp
  
  ! Cubature formula to use for calculating integrals
  INTEGER, INTENT(IN)                       :: ccubType
  
  ! The triangulation structure underlying to the discretisation.
  TYPE(t_triangulation), INTENT(IN), TARGET    :: rtriangulation
  
  ! The underlying domain.
  TYPE(t_boundary), INTENT(IN), TARGET         :: rdomain
  
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
  rspatialDiscr%p_rdomain              => rdomain
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
  p_relementDistr%ccubType = ccubType
  p_relementDistr%ccubTypeLin = ccubType
  
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
  
  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_initDiscr_combined (rspatialDiscr, &
                               ieltypTrial, ieltypTest, ccubType,&
                               rtriangulation, rdomain, rboundaryConditions)
  
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
  INTEGER, INTENT(IN)                       :: ieltypTrial
  
  ! The element type identifier that is to be used for all elements
  ! in the test space.
  INTEGER, INTENT(IN)                       :: ieltypTest

  ! Cubature formula to use for calculating integrals
  INTEGER, INTENT(IN)                       :: ccubType
  
  ! The triangulation structure underlying to the discretisation.
  TYPE(t_triangulation), INTENT(IN), TARGET    :: rtriangulation
  
  ! The underlying domain.
  TYPE(t_boundary), INTENT(IN), TARGET         :: rdomain
  
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
  rspatialDiscr%p_rdomain              => rdomain
  rspatialDiscr%ccomplexity            = SPDISC_UNIFORM
  
  ! All trial elements are ieltypTrial:
!  CALL storage_new1D ('spdiscr_initDiscr_combined', 'h_ItrialElements', &
!        rtriangulation%NEL, ST_INT, rspatialDiscr%h_ItrialElements,   &
!        ST_NEWBLOCK_NOINIT)
!  CALL storage_getbase_int (rspatialDiscr%h_ItrialElements,p_Iarray)
!  DO i=1,rtriangulation%NEL
!    p_Iarray(i) = ieltypTrial
!  END DO
  rspatialDiscr%h_ItrialElements = ST_NOHANDLE

  rspatialDiscr%bidenticalTrialAndTest = ieltypTrial .EQ. ieltypTest
  
  IF (rspatialDiscr%bidenticalTrialAndTest) THEN
    ! Identical trial and test space.
    ! Use the same handle for trial and test functions to save memory!
    rspatialDiscr%bidenticalTrialAndTest = .TRUE.
    rspatialDiscr%h_ItestElements = rspatialDiscr%h_ItrialElements  
  ELSE
    ! All test elements are ieltypTest.
    CALL storage_new1D ('spdiscr_initDiscr_combined', 'h_ItestElements', &
          rtriangulation%NEL, ST_INT, rspatialDiscr%h_ItestElements,   &
          ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int (rspatialDiscr%h_ItestElements,p_Iarray)
    DO i=1,rtriangulation%NEL
      p_Iarray(i) = ieltypTest
    END DO
  END IF

  ! Initialise the first element distribution
  rspatialDiscr%inumFESpaces           = 1
  p_relementDistr => rspatialDiscr%RelementDistribution(1)
  
  ! Initialise test and trial space for that block
  p_relementDistr%itrialElement = ieltypTrial
  p_relementDistr%itestElement = ieltypTest
  p_relementDistr%ccubType = ccubType
  p_relementDistr%ccubTypeLin = ccubType
  
  ! Get the typical transformation used with the element
  p_relementDistr%ctrafoType = elem_igetTrafoType(ieltypTrial)
  ctrafoTest = elem_igetTrafoType(ieltypTrial)
  
  IF (p_relementDistr%ctrafoType .NE. ctrafoTest) THEN
    PRINT *,'spdiscr_initDiscr_combined: Elements incompatible due to different'
    PRINT *,'transformation between reference and real element!'
    STOP
  END IF
  
  ! Check the cubature formula against the element distribution.
  ! This stops the program if this is not fulfilled.
  CALL spdiscr_checkCubature(ccubType,ieltypTrial)
  CALL spdiscr_checkCubature(ccubType,ieltypTest)

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
  NULLIFY(rspatialDiscr%p_rdomain)
  
  ! Release element identifier lists.
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
    
    ! Release the element list there
    CALL storage_free (p_relementDistr%h_IelementList)
    
    p_relementDistr%itrialElement = 0
    p_relementDistr%itestElement  = 0
    
  END DO
  
  ! No FE-spaces in here anymore...
  rspatialDiscr%inumFESpaces = 0
  
  ! Structure not initialised anymore
  rspatialDiscr%ndimension = 0
  
  END SUBROUTINE  

END MODULE
