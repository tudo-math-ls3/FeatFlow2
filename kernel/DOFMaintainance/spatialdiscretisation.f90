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
!# 1.) spdiscr_initDiscr_simple
!#     -> Initialise a discretisation structure. One element type for all
!#        geometric elements, for test and trial functions.
!#
!# 2.) spdiscr_releaseDiscr
!#     -> Release a discretisation structure.
!#
!# </purpose>
!##############################################################################

MODULE spatialdiscretisation

  USE fsystem
  USE storage
  USE triangulation
  USE boundary
  USE scalarbc
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
    INTEGER :: itrialElement        = 0
    
    ! Element identifier for test functions to use in this element list
    ! during the evaluation of linear and bilinear forms.
    INTEGER :: itestElement         = 0
    
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
  !  these structures, each one for one PDE. The structure is 'hung into'
  !  the scalar matrix that discretises that equation.
  
  TYPE t_spatialDiscretisation
  
    ! Dimension of the discretisation. 0=not initialised, 
    ! 2=2D discretisation, 3=3D discretisation
    INTEGER                          :: ndimension             = 0
  
    ! Pointer to the domain that is discretised
    TYPE(t_boundary), POINTER        :: p_rdomain              => NULL()
    
    ! Pointer to the underlying triangulation of the mesh (2D)
    TYPE(t_triangulation), POINTER :: p_rtriangulation     => NULL()

    ! Pointer to the analytical description of the boundary conditions
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions => NULL()
    
    ! Complexity of the discretisation. One of the SPDISC_xxxx constants
    INTEGER                          :: ccomplexity            = SPDISC_UNIFORM
    
    ! Flag: Trial and Test functions in all element distributions
    ! are the same. If FALSE, there's at öeast one element distribution
    ! with different trial and test functions.
    LOGICAL                          :: bidenticalTrialAndTest = .TRUE.
    
    ! Handle to trial function identifier list: For every geometric element, 
    ! an element identifier about which element to use.
    INTEGER                          :: h_ItrialElements       = ST_NOHANDLE

    ! Handle to test function identifier list: For every geometric element, 
    ! an element identifier about which element to use.
    ! Coincides with p_DtrialElements if bidenticalTrialAndTest=true!
    INTEGER                          :: h_ItestElements        = ST_NOHANDLE
    
    ! Number of different FE spaces mixed in this discretisation.
    ! This is the number of elements occupied in RelementDisttribution.
    ! Must be <= SPDISC_MAXFESPACES!
    INTEGER                     :: inumFESpaces           = 0
    
    ! List of element distribution structures for every element type
    ! that is used in the discretisation.
    TYPE(t_elementDistribution), DIMENSION(SPDISC_MAXFESPACES) :: RelementDistribution
    
  END TYPE
  
  !</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_releaseDiscr (p_rspatialDiscr,bkeepStructure)
  
!<description>
  ! This routine releases a discretisation structure from memory.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the structure p_rspatialDiscr is not released
  ! from memory. If set to FALSE or not existent (the usual setting), the 
  ! structure p_rspatialDiscr will also be removed from the heap after 
  ! cleaning up.
  LOGICAL, INTENT(IN), OPTIONAL :: bkeepStructure
!</input>

!<inputoutput>
  ! The discretisation structure to be released.
  TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscr
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  TYPE(t_elementDistribution), POINTER :: p_relementDistr

  IF (.NOT. ASSOCIATED(p_rspatialDiscr)) RETURN

  ! Cut the connection to the other structures
  NULLIFY(p_rspatialDiscr%p_rtriangulation)
  NULLIFY(p_rspatialDiscr%p_rdomain)
  NULLIFY(p_rspatialDiscr%p_rboundaryConditions)
  
  ! Release element identifier lists.
  ! The handles may coincide, so release them only once!
  IF (p_rspatialDiscr%h_ItestElements .NE. p_rspatialDiscr%h_ItrialElements) THEN
    CALL storage_free (p_rspatialDiscr%h_ItestElements)
  ELSE
    p_rspatialDiscr%h_ItestElements = ST_NOHANDLE
  END IF
  CALL storage_free (p_rspatialDiscr%h_ItrialElements)
  
  ! Loop through all element distributions
  DO i=1,p_rspatialDiscr%inumFESpaces
  
    p_relementDistr => p_rspatialDiscr%RelementDistribution(i)
    
    ! Release the element list there
    CALL storage_free (p_relementDistr%h_IelementList)
    
    p_relementDistr%itrialElement = 0
    p_relementDistr%itestElement  = 0
    
  END DO
  
  ! No FE-spaces in here anymore...
  p_rspatialDiscr%inumFESpaces = 0
  
  ! Structure not initialised anymore
  p_rspatialDiscr%ndimension = 0
  
  ! Deallocate the structure (if we are allowed to), finish.
  IF (.NOT. PRESENT(bkeepStructure)) THEN
    DEALLOCATE(p_rspatialDiscr)
  ELSE
    IF (.NOT. bkeepStructure) DEALLOCATE(p_rspatialDiscr)
  END IF

  END SUBROUTINE  

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_initDiscr_simple (p_rspatialDiscr,ieltyp, ccubType,&
                                       rtriangulation, rdomain, rboundaryConditions)
  
!<description>
  
  ! This routine initialises a discretisation structure for a uniform
  ! discretisation with one element for all geometric element primitives, 
  ! for trial as well as for test functions.
  !
  ! If p_rspatialDiscr is NULL(), a new structure will be created. Otherwise,
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
  TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscr
  
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
  TYPE(t_elementDistribution), POINTER :: p_relementDistr
  
  ! Do we have a structure?
  IF (.NOT. ASSOCIATED(p_rspatialDiscr)) THEN
    ALLOCATE(p_rspatialDiscr)
  ELSE
    ! Release the old structure without removing it from the heap.
    CALL spdiscr_releaseDiscr(p_rspatialDiscr,.TRUE.)
  END IF

  ! Initialise the variables of the structure for the simple discretisation
  p_rspatialDiscr%ndimension             = NDIM2D
  p_rspatialDiscr%p_rtriangulation       => rtriangulation
  p_rspatialDiscr%p_rdomain              => rdomain
  IF (PRESENT(rboundaryConditions)) THEN
    p_rspatialDiscr%p_rboundaryConditions  => rboundaryConditions
  ELSE
    NULLIFY(p_rspatialDiscr%p_rboundaryConditions)
  END IF
  p_rspatialDiscr%ccomplexity            = SPDISC_UNIFORM
  
  ! All trial elements are ieltyp:
  
  CALL storage_new1D ('spdiscr_initDiscr_simple', 'h_ItrialElements', &
        rtriangulation%NEL, ST_INT, p_rspatialDiscr%h_ItrialElements,   &
        ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (p_rspatialDiscr%h_ItrialElements,p_Iarray)
  DO i=1,rtriangulation%NEL
    p_Iarray(i) = ieltyp
  END DO
  
  ! All test elements are ieltyp.
  ! Use the same handle for trial and test functions to save memory!
  p_rspatialDiscr%bidenticalTrialAndTest = .TRUE.
  p_rspatialDiscr%h_ItestElements = p_rspatialDiscr%h_ItrialElements  
  
  ! Initialise the first element distribution
  p_rspatialDiscr%inumFESpaces           = 1
  p_relementDistr => p_rspatialDiscr%RelementDistribution(1)
  
  ! Initialise test and trial space for that block
  p_relementDistr%itrialElement = ieltyp
  p_relementDistr%itestElement = ieltyp
  p_relementDistr%ccubType = ccubType
  p_relementDistr%ccubTypeLin = ccubType
  
  ! Get the typical transformation used with the element
  p_relementDistr%ctrafoType = elem_getTrafoType(ieltyp)
  
  ! Initialise an 'identity' array containing the numbers of all elements.
  CALL storage_new1D ('spdiscr_initDiscr_simple', 'h_IelementList', &
        rtriangulation%NEL, ST_INT, p_relementDistr%h_IelementList,   &
        ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (p_relementDistr%h_IelementList,p_Iarray)
  DO i=1,rtriangulation%NEL
    p_Iarray(i) = i
  END DO
  
  END SUBROUTINE  
  
END MODULE
