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
!# </purpose>
!##############################################################################

MODULE spatialdiscretisation

  USE fsystem
  USE storage
  USE triangulation
  USE boundary
  USE transformation
  USE element
  USE cubature
  USE genoutput
  
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

!<constantblock description="Additional constants for cubature formulas">

  ! Automatically determine cubature formula for a discretisation.
  INTEGER, PARAMETER :: SPDISC_CUB_AUTOMATIC = 0

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
  
    ! Element identifier for Finite Element functions to use in this 
    ! element list.
    INTEGER(I32) :: celement        = EL_UNDEFINED
    
    ! Cubature formula to use for the discretisation of this element pair
    ! during the evaluation of bilinear forms (matrix generation).
    ! Note: When evaluating bilinear forms, the ccubTypeBilForm
    ! constant of the test space decides about the cubature formula
    ! to be used!
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
    ! combination of trial/test functions.
    ! If NEL=0, the element list is empty, i.e. h_IelementList = ST_NOHANDLE!
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
    ! 1=1D discretisation, 2=2D discretisation, 3=3D discretisation
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
    
    ! Handle to the element distribution identifier list.
    ! For every geometric element i, IelementDistr(i) specifies the
    ! number of the element distribution that contains that element.
    ! That way one can easily access information; e.g. retrieving the
    ! element type would be possible as follows:
    !   RelementDistr(IelementDistr(i))%itrialElement
    ! In a uniform discretisation (ccomplexity=SPDISC_UNIFORM), this
    ! handle is ST_NOHANDLE as all elements are in the
    ! element distribution 1.
    INTEGER                          :: h_IelementDistr       = ST_NOHANDLE
    
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
    ! This is the number of elements occupied in RelementDistibution.
    INTEGER                          :: inumFESpaces           = 0
    
    ! List of element distribution structures for every element type
    ! that is used in the discretisation.
    TYPE(t_elementDistribution), DIMENSION(:), POINTER :: RelementDistr => NULL()
    
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
    ! 1=1D discretisation, 2=2D discretisation, 3=3D discretisation
    INTEGER                          :: ndimension             = 0

    ! Complexity of the discretisation. One of the SPDISC_xxxx constants.
    ! SPDISC_UNIFORM = all elements in each discretisation
    !   substructure RspatialDiscr(:) are the same.
    ! SPDISC_CONFORMAL = Elements of different FE spaces are mixed,
    !   but the DOF's 'fit together'. Each discretisation substructure 
    !   RspatialDiscr(:) has exactly the same number of element
    !   distributions, and each element distribution 
    !     RspatialDiscr(1)%Relementistributions(i), 
    !     RspatialDiscr(2)%Relementistributions(i),
    !     RspatialDiscr(3)%Relementistributions(i),...
    !   describe exactly the same set of elements (Same size, same type,
    !   same order in the element lists,...).
    INTEGER                          :: ccomplexity            = SPDISC_UNIFORM
  
    ! Pointer to the domain that is discretised
    TYPE(t_boundary), POINTER        :: p_rboundary            => NULL()
    
    ! Pointer to the underlying triangulation of the mesh (2D)
    TYPE(t_triangulation), POINTER   :: p_rtriangulation     => NULL()

    ! Number of solution components maintained by this structure.
    INTEGER                             :: ncomponents
    
    ! A list of up to ncomponents scalar spatial discretisation structures.
    ! Each structure corresponds to one solution component and defines
    ! trial-/test-functions, complexity of the discretisation etc.
    TYPE(t_spatialDiscretisation), DIMENSION(:), POINTER :: RspatialDiscr => NULL()
    
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

  INTEGER :: NVE, idim
  LOGICAL :: bcompatible

  ! Get from the element distribution the trial space and from that
  ! the number of vertices, the element expects.
  NVE = elem_igetNVE(ielementType)
  idim = elem_igetDimension(ielementType)
  
  bcompatible = .TRUE.
  
  ! Now we directly access the cubature constants in cubature.f90!
  ! This is the only point in the kernel where this is necessary.
  
  ! 1D: Line?
  IF (ccubType .LE. 99) THEN
    IF ((NVE .NE. 2) .OR. (idim .NE. NDIM1D)) bcompatible = .FALSE.
  END IF

  ! 2D: Quad?
  IF ((ccubType .GE. 200) .AND. (ccubType .LE. 249)) THEN
    ! Tri?
    IF ((NVE .NE. 4) .OR. (idim .NE. NDIM2D)) bcompatible = .FALSE.
  END IF
  
  ! 2D: Tri?
  IF ((ccubType .GE. 250) .AND. (ccubType .LE. 299)) THEN
    ! Quad?
    IF ((NVE .NE. 3) .OR. (idim .NE. NDIM2D)) bcompatible = .FALSE.
  END IF
  
  ! 3D: Hexa?
  IF ((ccubType .GE. 300) .AND. (ccubType .LE. 349)) THEN
    IF ((NVE .NE. 8) .OR. (idim .NE. NDIM3D)) bcompatible = .FALSE.
  END IF
  
  ! 3D: Tetra?
  IF ((ccubType .GE. 350) .AND. (ccubType .LE. 499)) THEN
    IF ((NVE .NE. 4) .OR. (idim .NE. NDIM3D)) bcompatible = .FALSE.
  END IF
  
  ! Q2T with bubble does not work with G1X1, Trapezoidal rule and
  ! G2X2 -- Laplace matrices would be indefinite because of the definition
  ! of the bubble function!
  IF (ielementType .EQ. EL_Q2TB) THEN
    IF ((ccubType .GE. 201) .AND. (ccubType .LE. 204)) bcompatible = .FALSE.
  END IF

  IF (.NOT. bcompatible) THEN
    CALL output_line ('Element and cubature formula not compatible!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_checkCubature')  
    CALL sys_halt()
  END IF
  
  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION spdiscr_getLumpCubature (ielementType) RESULT (ccubType)
  
!<description>
  ! This routine tries to determine a cubature formula identifier according
  ! to a given element type, such that the corresponding mass matrix will
  ! get diagonal (mass lumping). If this is not possible, 0 is returned.
!</description>

!<input>
  ! An element type identifier
  INTEGER(I32), INTENT(IN)                       :: ielementType
!</input>

!<result>
  ! A cubature formula identifier that will diagonalise the mass matrix,
  ! or 0 if such an identifier is unknown / not possible.
!</result>
  
!</function>

    SELECT CASE (elem_igetDimension(ielementType))
    CASE (NDIM1D)
    
      SELECT CASE (elem_getPrimaryElement(ielementType))
      CASE (EL_P0_1D)
        ! Use Gauss-1
        ccubType = CUB_G1_1D

      CASE (EL_P1_1D)
        ! Use trapezoidal rule
        ccubType = CUB_TRZ_1D
      
      CASE (EL_P2_1D)
        ! Use Gauss-2
        ccubType = CUB_G2_1D

      CASE (EL_S31_1D)
        ! Use Gauss-4
        ccubType = CUB_G4_1D
      
      CASE DEFAULT
        ccubType = 0
      END SELECT

    CASE (NDIM2D)
    
      SELECT CASE (elem_getPrimaryElement(ielementType))
      CASE (EL_P0)
        ! Use Gauss 1X1
        ccubType = CUB_G1_T

      CASE (EL_P1)
        ! Use trapezoidal rule
        ccubType = CUB_TRZ_T

      CASE (EL_P1T)
        ! Use Gauss-3pt
        ccubType = CUB_G3_T

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
      
    CASE (NDIM3D)

      CALL output_line ('3D not supported.', &
        OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_getLumpCubature')  
      CALL sys_halt()
      
    CASE DEFAULT
      ccubType = 0
    END SELECT

  END FUNCTION
  
  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION spdiscr_getStdCubature (ielementType) RESULT (ccubType)
  
!<description>
  ! This routine returns a standard cubature formula for an element which
  ! can be used as default when setting up matrices/vectors.
  ! If this is not possible, 0 is returned.
!</description>

!<input>
  ! An element type identifier
  INTEGER(I32), INTENT(IN)                       :: ielementType
!</input>

!<result>
  ! A standard cubature formula for the assembly of matrices/vectors
  ! with the specified element ielementType.
!</result>
  
!</function>

    SELECT CASE (elem_igetDimension(ielementType))
    CASE (NDIM1D)
    
      SELECT CASE (elem_getPrimaryElement(ielementType))
      CASE (EL_P0_1D)
        ! 1-point Gauss
        ccubType = CUB_G1_1D

      CASE (EL_P1_1D)
        ! 2-point Gauss
        ccubType = CUB_G2_1D
      
      CASE (EL_P2_1D)
        ! 3-point Gauss
        ccubType = CUB_G3_1D

      CASE (EL_S31_1D)
        ! 4-point Gauss
        ccubType = CUB_G4_1D
      
      CASE DEFAULT
        ccubType = 0
      END SELECT

    CASE (NDIM2D)
    
      SELECT CASE (elem_getPrimaryElement(ielementType))
      CASE (EL_P0)
        ! Use Gauss 1X1
        ccubType = CUB_G1_T

      CASE (EL_P1)
        ! Use Gauss-3pt
        ccubType = CUB_G3_T

      CASE (EL_P2)
        ! Gauss-3pt
        ccubType = CUB_G3_T

      CASE (EL_P1T)
        ! Gauss-3pt
        ccubType = CUB_G3_T

      CASE (EL_Q0)
        ! 1x1 Gauss formula
        ccubType = CUB_G1X1

      CASE (EL_Q1)
        ! 2x2 Gauss formula
        ccubType = CUB_G2X2

      CASE (EL_Q2)
        ! 3x3 Gauss formula
        ccubType = CUB_G3X3

      CASE (EL_Q3)
        ! 3x3 Gauss formula
        ccubType = CUB_G3X3

      CASE (EL_Q1T)
        ! 2x2 Gauss formula
        ccubType = CUB_G2X2

      CASE (EL_QP1)
        ! 2x2 Gauss formula
        ccubType = CUB_G2X2
      
      CASE DEFAULT
        ccubType = 0
      END SELECT
      
    CASE (NDIM3D)

      CALL output_line ('3D not supported.', &
        OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_getStdCubature')  
      CALL sys_halt()
      
    CASE DEFAULT
      ccubType = 0
    END SELECT

  END FUNCTION
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_initBlockDiscr (rblockDiscr,ncomponents,&
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
  TYPE(t_triangulation), INTENT(IN), TARGET    :: rtriangulation
  
  ! Number of solution components maintained by the block structure
  INTEGER, INTENT(IN), OPTIONAL                :: ncomponents
  
  ! OPTIONAL: The underlying domain.
  TYPE(t_boundary), INTENT(IN), TARGET, OPTIONAL :: rboundary

!</input>
  
!<output>
  
  ! The block discretisation structure to be initialised.
  TYPE(t_blockDiscretisation), INTENT(OUT) :: rblockDiscr
  
!</output>
  
!</subroutine>

  ! Initialise the variables of the structure for the simple discretisation
  rblockDiscr%ndimension             = rtriangulation%ndim
  rblockDiscr%ccomplexity            = SPDISC_UNIFORM
  rblockDiscr%p_rtriangulation       => rtriangulation
  IF (PRESENT(rboundary)) THEN
    rblockDiscr%p_rboundary            => rboundary
  ELSE
    NULLIFY(rblockDiscr%p_rboundary)
  END IF

  rblockDiscr%ncomponents            = ncomponents
  ALLOCATE(rblockDiscr%RspatialDiscr(ncomponents))

  ! That's it.  
  
  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_initBlockDiscr2D (rblockDiscr,ncomponents,&
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
  TYPE(t_triangulation), INTENT(IN), TARGET      :: rtriangulation
  
  ! Number of solution components maintained by the block structure
  INTEGER, INTENT(IN)                            :: ncomponents

  ! OPTIONAL: The underlying domain.
  TYPE(t_boundary), INTENT(IN), TARGET, OPTIONAL :: rboundary
  
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
  IF (PRESENT(rboundary)) THEN
    rblockDiscr%p_rboundary          => rboundary
  ELSE
    NULLIFY(rblockDiscr%p_rboundary)
  END IF

  rblockDiscr%ncomponents            = ncomponents
  ALLOCATE(rblockDiscr%RspatialDiscr(ncomponents))

  ! That's it.  
  
  END SUBROUTINE  

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_createBlockDiscrInd (rspatialDiscr,rblockDiscr)
  
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
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rspatialDiscr
  
!</input>
  
!<output>
  
  ! The block discretisation structure to be initialised.
  TYPE(t_blockDiscretisation), INTENT(OUT) :: rblockDiscr
  
!</output>
  
!</subroutine>

    ! Initialise a new block discretisation with one component.
    CALL spdiscr_initBlockDiscr2D (rblockDiscr,1,&
        rspatialDiscr%p_rtriangulation, rspatialDiscr%p_rboundary)
    
    ! Put a copy of the spatial discretisation to first component
    ! of the the block discretisation. Share the data.
    CALL spdiscr_duplicateDiscrSc (rspatialDiscr,&
        rblockDiscr%RspatialDiscr(1),.TRUE.)
  
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
  ! Any old existing information in rdestDiscr is released if necessary.
  TYPE(t_blockDiscretisation), INTENT(INOUT), TARGET :: rdestDiscr
!</output>
  
!</subroutine>

    ! local variables
    INTEGER :: ifirst, ilast, ncount, i
    
    ! Check that the source discretisation structure is valid.
    IF (rsourceDiscr%ndimension .LE. 0) THEN
      CALL output_line ('Source structure invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveBlockDiscr')  
      CALL sys_halt()
    END IF

   ! Release old information if present
    CALL spdiscr_releaseBlockDiscr(rdestDiscr,.TRUE.)

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
    rdestDiscr%ncomponents            =  ncount       
    
    ! Copy all substructures -- from ifirstBlock to ilastBlock.
    ! Use spdiscr_duplicateDiscrSc which savely copies the scalar discretisation
    ! structures. We set bshare=.TRUE. here, so the information is shared
    ! between the source and destination structure; the dynamic information 
    ! 'belongs' to rdiscrSource and not to the newly created rdiscrDest!
    ALLOCATE(rdestDiscr%RspatialDiscr(ncount))
    DO i=1,ncount
      CALL spdiscr_duplicateDiscrSc (rsourceDiscr%RspatialDiscr(ifirst+i-1), &
          rdestDiscr%RspatialDiscr(i), .TRUE.)
    END DO
      
    END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_releaseBlockDiscr (rblockDiscr, breleaseSubstruc)
  
!<description>
  ! This routine releases a block discretisation structure from memory.
!</description>

!<input>
  ! OPTIONAL: Release substructures.
  ! If set to TRUE, the memory of all scalar spatial discretisation structures
  !   in rblockDiscr is also released from memory. This is the standard setting.
  ! Is set to FALSE, only rblockDiscr is cleaned up, the substructures
  !   are ignored.
  LOGICAL, INTENT(IN), OPTIONAL :: breleaseSubstruc
!</input>

!<inputoutput>
  ! The block discretisation structures to be released.
  TYPE(t_blockDiscretisation), INTENT(INOUT), TARGET :: rblockDiscr
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  LOGICAL :: brelsub
  
  brelsub = .TRUE.
  IF (PRESENT(breleaseSubstruc)) brelsub = breleaseSubstruc

  ! Cut the connection to the other structures
  NULLIFY(rblockDiscr%p_rtriangulation)
  NULLIFY(rblockDiscr%p_rboundary)
  
  ! Release substructures?
  IF (brelsub) THEN
    IF (ASSOCIATED(rblockDiscr%RspatialDiscr)) THEN
      DO i=1,rblockDiscr%ncomponents
        CALL spdiscr_releaseDiscr (rblockDiscr%RspatialDiscr(i))
      END DO
    END IF
  END IF
  IF (ASSOCIATED(rblockDiscr%RspatialDiscr)) &
    DEALLOCATE(rblockDiscr%RspatialDiscr)
  rblockDiscr%ncomponents = 0

  ! Structure not initialised anymore
  rblockDiscr%ndimension = 0
  
  END SUBROUTINE  

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_initDiscr_simple (rspatialDiscr,ieltyp, ccubType,&
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
  INTEGER(I32), INTENT(IN)                       :: ieltyp
  
  ! Cubature formula CUB_xxxx to use for calculating integrals.
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  INTEGER, INTENT(IN)                       :: ccubType
  
  ! The triangulation structure underlying to the discretisation.
  TYPE(t_triangulation), INTENT(IN), TARGET    :: rtriangulation
  
  ! The underlying domain.
  TYPE(t_boundary), INTENT(IN), TARGET, OPTIONAL :: rboundary
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
  INTEGER :: ccub

  ! Automatically determine cubature formula if necessary  
  ccub = ccubType
  IF (ccub .EQ. SPDISC_CUB_AUTOMATIC) &
      ccub = spdiscr_getStdCubature(ieltyp)
  
  ! Do we have a structure?
  IF (rspatialDiscr%ndimension .NE. 0) THEN
    ! Release the old structure.
    CALL spdiscr_releaseDiscr(rspatialDiscr)
  END IF

  ! Initialise the variables of the structure for the simple discretisation
  !rspatialDiscr%ndimension             = NDIM2D
  rspatialDiscr%ndimension             = rtriangulation%ndim
  rspatialDiscr%p_rtriangulation       => rtriangulation
  IF (PRESENT(rboundary)) THEN
    rspatialDiscr%p_rboundary            => rboundary
  ELSE
    NULLIFY(rspatialDiscr%p_rboundary)
  END IF
  rspatialDiscr%ccomplexity            = SPDISC_UNIFORM
  
  ! All trial elements are ieltyp:
  
!  CALL storage_new1D ('spdiscr_initDiscr_simple', 'h_ItrialElements', &
!        rtriangulation%NEL, ST_INT, rspatialDiscr%h_ItrialElements,   &
!        ST_NEWBLOCK_NOINIT)
!  CALL storage_getbase_int (rspatialDiscr%h_ItrialElements,p_Iarray)
!  DO i=1,rtriangulation%NEL
!    p_Iarray(i) = ieltyp
!  END DO
  rspatialDiscr%h_IelementDistr = ST_NOHANDLE
  
  ! Initialise the first element distribution
  rspatialDiscr%inumFESpaces           = 1
  ALLOCATE(rspatialDiscr%RelementDistr(rspatialDiscr%inumFESpaces))
  p_relementDistr => rspatialDiscr%RelementDistr(1)
  
  ! Initialise FE space for that block
  p_relementDistr%celement = ieltyp
  p_relementDistr%ccubTypeBilForm = ccub
  p_relementDistr%ccubTypeLinForm = ccub
  p_relementDistr%ccubTypeEval = ccub
  
  ! Get the typical transformation used with the element
  p_relementDistr%ctrafoType = elem_igetTrafoType(ieltyp)
  
  ! Check the cubature formula against the element distribution.
  ! This stops the program if this is not fulfilled.
  CALL spdiscr_checkCubature(ccub,ieltyp)

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
  INTEGER(I32), INTENT(IN)                       :: ieltypTri

  ! The element type identifier that is to be used for all quadrilateral elements.
  INTEGER(I32), INTENT(IN)                       :: ieltypQuad
  
  ! Cubature formula CUB_xxxx to use for calculating integrals 
  ! on triangular elements
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  INTEGER, INTENT(IN)                       :: ccubTypeTri

  ! Cubature formula CUB_xxxx to use for calculating integrals on 
  ! quadrilateral elements
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  INTEGER, INTENT(IN)                       :: ccubTypeQuad
  
  ! The triangulation structure underlying to the discretisation.
  TYPE(t_triangulation), INTENT(IN), TARGET    :: rtriangulation
  
  ! The underlying domain.
  TYPE(t_boundary), INTENT(IN), TARGET, OPTIONAL :: rboundary
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
  INTEGER :: ccubTri,ccubQuad
  
  ! Automatically determine cubature formula if necessary  
  ccubTri = ccubTypeTri
  IF (ccubTri .EQ. SPDISC_CUB_AUTOMATIC) &
      ccubTri = spdiscr_getStdCubature(ieltypTri)
  ccubQuad = ccubTypeQuad
  IF (ccubQuad .EQ. SPDISC_CUB_AUTOMATIC) &
      ccubQuad = spdiscr_getStdCubature(ieltypQuad)
  
  ! Do we have a structure?
  IF (rspatialDiscr%ndimension .NE. 0) THEN
    ! Release the old structure.
    CALL spdiscr_releaseDiscr(rspatialDiscr)
  END IF

  ! Initialise the variables of the structure for the simple discretisation
  rspatialDiscr%ndimension             = NDIM2D
  rspatialDiscr%p_rtriangulation       => rtriangulation
  IF (PRESENT(rboundary)) THEN
    rspatialDiscr%p_rboundary          => rboundary
  ELSE
    NULLIFY(rspatialDiscr%p_rboundary)
  END IF
  rspatialDiscr%ccomplexity            = SPDISC_CONFORMAL
  
  ! Allocate an array containing the element distribution for each element
  CALL storage_new1D ('spdiscr_initDiscr_triquad', 'h_ItrialElements', &
        rtriangulation%NEL, ST_INT, rspatialDiscr%h_IelementDistr,   &
        ST_NEWBLOCK_NOINIT)
  CALL storage_getbase_int (rspatialDiscr%h_IelementDistr,p_Iarray)

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
        ! Triangular elements are in element distribution 1
        p_Iarray(i) = 1
        
        ! This is the IelemCount(1)'th triangle
        IelemCount(1) = IelemCount(1)+1
        p_IelementCounter(i) = IelemCount(1)
      ELSE
        ! Quad elements are in element distribution 2
        p_Iarray(i) = 2

        ! This is the IelemCount(2)'th quad
        IelemCount(2) = IelemCount(2)+1
        p_IelementCounter(i) = IelemCount(2)
      END IF
    END DO
  ELSE
    ! Pure triangular mesh
    DO i=1,rtriangulation%NEL
      ! Triangular elements are in element distribution 1
      p_Iarray(i) = 1
      
      ! This is the IelemCount(1)'th triangle
      IelemCount(1) = IelemCount(1)+1
      p_IelementCounter(i) = IelemCount(1)
    END DO
  END IF
  
  ! Initialise the first element distribution
  rspatialDiscr%inumFESpaces           = 2
  ALLOCATE(rspatialDiscr%RelementDistr(rspatialDiscr%inumFESpaces))
  p_relementDistrTria => rspatialDiscr%RelementDistr(1)
  p_relementDistrQuad => rspatialDiscr%RelementDistr(2)
  
  ! Initialise test and trial space for that block
  p_relementDistrTria%celement = ieltypTri
  p_relementDistrTria%ccubTypeBilForm = ccubTri
  p_relementDistrTria%ccubTypeLinForm = ccubTri
  p_relementDistrTria%ccubTypeEval = ccubTri

  p_relementDistrQuad%celement = ieltypQuad
  p_relementDistrQuad%ccubTypeBilForm = ccubQuad
  p_relementDistrQuad%ccubTypeLinForm = ccubQuad
  p_relementDistrQuad%ccubTypeEval = ccubQuad
  
  ! Get the typical transformation used with the element
  p_relementDistrTria%ctrafoType = elem_igetTrafoType(ieltypTri)
  p_relementDistrQuad%ctrafoType = elem_igetTrafoType(ieltypQuad)
  
  ! Check the cubature formula against the element distribution.
  ! This stops the program if this is not fulfilled.
  CALL spdiscr_checkCubature(ccubTri,ieltypTri)
  CALL spdiscr_checkCubature(ccubQuad,ieltypQuad)

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
    
    IF (UBOUND(p_IverticesAtElement,1) .GE. TRIA_NVEQUAD2D) THEN
      ! There are quads and probably triangles in the mesh
      DO i=1,rtriangulation%NEL
        IF (p_IverticesAtElement(TRIA_NVEQUAD2D,i) .EQ. 0) THEN
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
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  INTEGER, INTENT(IN)                       :: ccubType
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  ! Any old existing discretisation information in rdestDiscr
  ! is released if necessary.
  TYPE(t_spatialDiscretisation), INTENT(INOUT), TARGET :: rdestDiscr
!</output>
  
!</subroutine>

  ! local variables
  ! INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
  ! TYPE(t_elementDistribution), POINTER :: p_relementDistr
  INTEGER :: ccub

  ! Automatically determine cubature formula if necessary  
  ccub = ccubType
  IF (ccub .EQ. SPDISC_CUB_AUTOMATIC) &
      ccub = spdiscr_getStdCubature(ieltyp)
  
  ! Check that the source discretisation structure is valid.
  IF (rsourceDiscr%ndimension .LE. 0) THEN
    CALL output_line ('Source structure invalid!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveSimpleDiscr')  
    CALL sys_halt()
  END IF
  
  ! Check that the discretisation structure is really uniform.
  ! More complex situations are not supported by this routine.
  IF (rsourceDiscr%ccomplexity .NE. SPDISC_UNIFORM) THEN
    CALL output_line ('Only uniform discretisations supported!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveSimpleDiscr')  
    CALL sys_halt()
  END IF
  
  IF (elem_igetDimension(rsourceDiscr%RelementDistr(1)%celement) .ne. &
      elem_igetDimension(ieltyp)) THEN
    CALL output_line ('Element dimension different!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_deriveSimpleDiscr')  
    CALL sys_halt()
  END IF
  
  ! Release old information if present
  CALL spdiscr_releaseDiscr(rdestDiscr)
  
  ! Check the cubature formula against the element distribution.
  ! This stops the program if this is not fulfilled.
  CALL spdiscr_checkCubature(ccub,ieltyp)
  
  ! Copy the source structure to the destination.
  ! This copies all handles and hence all dynamic information
  rdestDiscr = rsourceDiscr
  
  ! Allocate a new element distribution and copy content from source
  ALLOCATE(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
  rdestDiscr%RelementDistr(1:rdestDiscr%inumFESpaces) = &
      rsourceDiscr%RelementDistr(1:rsourceDiscr%inumFESpaces)

  ! Change the element type of all trial functions to ieltyp
  rdestDiscr%RelementDistr(1)%celement = ieltyp
  
  ! Init the cubature rule
  rdestDiscr%RelementDistr(1)%ccubTypeBilForm = ccub
  rdestDiscr%RelementDistr(1)%ccubTypeLinForm = ccub
  rdestDiscr%RelementDistr(1)%ccubTypeEval = ccub
  
  ! Get the typical transformation used with the element
  rdestDiscr%RelementDistr(1)%ctrafoType = elem_igetTrafoType(ieltyp)
  
  ! Mark the new discretisation structure as 'copy', to prevent
  ! the dynamic information to be released.
  ! The dynamic information 'belongs' to rdiscrSource and not to the
  ! newly created rdiscrDest!
  rdestDiscr%bisCopy = .TRUE.
  
  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_deriveDiscr_triquad (&
      rsourceDiscr, ieltypTri, ieltypQuad, ccubTypeTri, ccubTypeQuad,&
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

  ! The element type identifier that is to be used for all triangular
  ! elements in the new discretisation structure
  INTEGER(I32), INTENT(IN)                       :: ieltypTri

  ! The element type identifier that is to be used for all quad
  ! elements in the new discretisation structure
  INTEGER(I32), INTENT(IN)                       :: ieltypQuad
  
  ! Cubature formula to use for calculating integrals on triangular
  ! elements in the new discretisation structure.
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  INTEGER, INTENT(IN)                       :: ccubTypeTri

  ! Cubature formula to use for calculating integrals on quad
  ! elements in the new discretisation structure.
  ! Alternatively, the value SPDISC_CUB_AUTOMATIC means: 
  ! automatically determine cubature formula.
  INTEGER, INTENT(IN)                       :: ccubTypeQuad
!</input>
  
!<output>
  ! The discretisation structure to be initialised.
  ! Any old existing discretisation information in rdestDiscr
  ! is released if necessary.
  TYPE(t_spatialDiscretisation), INTENT(INOUT), TARGET :: rdestDiscr
!</output>
  
!</subroutine>

    ! local variables
    ! INTEGER(I32), DIMENSION(:), POINTER :: p_Iarray
    ! TYPE(t_elementDistribution), POINTER :: p_relementDistr
    INTEGER :: ccubTri,ccubQuad,idistr,nve

    ! Automatically determine cubature formula if necessary  
    ccubTri = ccubTypeTri
    IF (ccubTri .EQ. SPDISC_CUB_AUTOMATIC) &
        ccubTri = spdiscr_getStdCubature(ieltypTri)
    ccubQuad = ccubTypeQuad
    IF (ccubQuad .EQ. SPDISC_CUB_AUTOMATIC) &
        ccubQuad = spdiscr_getStdCubature(ieltypQuad)
    
    ! Check that the source discretisation structure is valid.
    IF (rsourceDiscr%ndimension .NE. NDIM2D) THEN
      CALL output_line ('Source structure invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdiscr_deriveSimpleDiscr_triquad')  
      CALL sys_halt()
    END IF
    
    ! Check that the discretisation structure is really uniform.
    ! More complex situations are not supported by this routine.
    IF ((rsourceDiscr%ccomplexity .NE. SPDISC_UNIFORM) .AND. &
        (rsourceDiscr%ccomplexity .NE. SPDISC_CONFORMAL)) THEN
      CALL output_line ('Only uniform discretisations supported!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdiscr_deriveSimpleDiscr_triquad')  
      CALL sys_halt()
    END IF
    
    ! Release old information if present
    CALL spdiscr_releaseDiscr(rdestDiscr)
    
    ! Check the cubature formula against the element distribution.
    ! This stops the program if this is not fulfilled.
    CALL spdiscr_checkCubature(ccubTri,ieltypTri)
    CALL spdiscr_checkCubature(ccubQuad,ieltypQuad)
    
    ! Copy the source structure to the destination.
    ! This copies all handles and hence all dynamic information
    rdestDiscr = rsourceDiscr
    
    ! Allocate a new element distribution
    ALLOCATE(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
    rdestDiscr%RelementDistr(1:rdestDiscr%inumFESpaces) = &
        rsourceDiscr%RelementDistr(1:rsourceDiscr%inumFESpaces)

    ! Loop through the element distributions...
    DO idistr = 1,rdestDiscr%inumFESpaces
    
      ! Check the element there. If it's a triangular element,
      ! change the element type to ielTypTri. If it's a quad
      ! element, change the element type to ielTypQuad.
      nve = elem_igetNVE(rsourceDiscr%RelementDistr(idistr)%celement)
      SELECT CASE (nve)
      CASE (TRIA_NVETRI2D)
        rdestDiscr%RelementDistr(idistr)%celement = ieltypTri

        ! Init the cubature rule
        rdestDiscr%RelementDistr(idistr)%ccubTypeBilForm = ccubTri
        rdestDiscr%RelementDistr(idistr)%ccubTypeLinForm = ccubTri
        rdestDiscr%RelementDistr(idistr)%ccubTypeEval = ccubTri

        ! Get the typical transformation used with the element
        rdestDiscr%RelementDistr(idistr)%ctrafoType = &
            elem_igetTrafoType(ieltypTri)
        
      CASE (TRIA_NVEQUAD2D)
        rdestDiscr%RelementDistr(idistr)%celement = ieltypQuad

        ! Init the cubature rule
        rdestDiscr%RelementDistr(idistr)%ccubTypeBilForm = ccubQuad
        rdestDiscr%RelementDistr(idistr)%ccubTypeLinForm = ccubQuad
        rdestDiscr%RelementDistr(idistr)%ccubTypeEval = ccubQuad

        ! Get the typical transformation used with the element
        rdestDiscr%RelementDistr(idistr)%ctrafoType = &
            elem_igetTrafoType(ieltypQuad)
      END SELECT
      
    END DO
  
    ! Mark the new discretisation structure as 'copy', to prevent
    ! the dynamic information to be released.
    ! The dynamic information 'belongs' to rdiscrSource and not to the
    ! newly created rdiscrDest!
    rdestDiscr%bisCopy = .TRUE.
  
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
  
  IF (.NOT. rspatialDiscr%bisCopy) THEN
  
    ! Release element distribution lists.
    IF (rspatialDiscr%ccomplexity .NE. SPDISC_UNIFORM) THEN
      CALL storage_free (rspatialDiscr%h_IelementDistr)
    END IF
    
  ELSE
    rspatialDiscr%h_IelementDistr = ST_NOHANDLE
  END IF
  
  ! Loop through all element distributions
  DO i=1,rspatialDiscr%inumFESpaces
  
    p_relementDistr => rspatialDiscr%RelementDistr(i)
    
    ! If the element distribution is empty, skip it
    IF (p_relementDistr%NEL .NE. 0) THEN
    
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
      
    END IF
    
    p_relementDistr%celement = EL_UNDEFINED
    
  END DO
  
  IF (.NOT. rspatialDiscr%bisCopy) THEN
    IF (rspatialDiscr%h_IelementCounter .NE. ST_NOHANDLE) &
      CALL storage_free (rspatialDiscr%h_IelementCounter)
  ELSE    
    rspatialDiscr%h_IelementCounter = ST_NOHANDLE
  END IF
  
  ! No FE-spaces in here anymore...
  IF (ASSOCIATED(rspatialDiscr%RelementDistr)) &
    DEALLOCATE(rspatialDiscr%RelementDistr)
  rspatialDiscr%inumFESpaces = 0
  
  ! Structure not initialised anymore
  rspatialDiscr%ndimension = 0
  
  END SUBROUTINE  

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_duplicateDiscrSc (rsourceDiscr, rdestDiscr, bshare)
  
!<description>
  ! This routine creates a copy of the discretisation structure rsourceDiscr.
  ! Depending on bshare, the destination structure rdestDiscr will either
  ! obtain a 'simple' copy (i.e. sharing all handles and all information
  ! with rsourceDiscr) or a separate copy (which costs memory for all the
  ! element information!).
!</description>

!<input>
  ! A source discretisation structure that should be used as template
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rsourceDiscr
  
  ! OPTIONAL: Whether the new discretisation structure should share its information
  ! with rsourceDiscr.
  ! =FALSE: Create a complete copy of rsourceDiscr which is independent
  !  of rsourceDiscr.
  ! =TRUE: The new discretisation will not be a complete new structure, but a
  !  'derived' structure, i.e. it uses the same dynamic information
  !  (handles and therefore element lists) as rsourceDiscr.
  ! If not specified, TRUE is assumed.
  LOGICAL, INTENT(IN), OPTIONAL :: bshare
!</input>
  
!<output>
  ! The new discretisation structure. Any old existing information in rdestDiscr
  ! is released if necessary.
  TYPE(t_spatialDiscretisation), INTENT(INOUT), TARGET :: rdestDiscr
!</output>
  
!</subroutine>

    LOGICAL :: bshr
    
    bshr = .TRUE.
    IF (PRESENT(bshare)) bshr = bshare

    ! Release old information if present
    CALL spdiscr_releaseDiscr(rdestDiscr)

    ! Currently, this routine supports only bshare=TRUE!
    IF (bshr) THEN
    
      ! Copy all information
      rdestDiscr = rsourceDiscr
      
      ! Duplicate the element distribution structure
      ALLOCATE(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
      rdestDiscr%RelementDistr = rsourceDiscr%RelementDistr
    
      ! Mark the new discretisation structure as 'copy', to prevent
      ! the dynamic information to be released.
      ! The dynamic information 'belongs' to rdiscrSource and not to the
      ! newly created rdiscrDest!
      rdestDiscr%bisCopy = .TRUE.
      
    ELSE
      CALL output_line ('bshare=FALSE currently not supported!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_duplicateDiscrSc')  
      CALL sys_halt()
    END IF
  
  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_duplicateBlockDiscr (rsourceDiscr, rdestDiscr, bshare)
  
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
  TYPE(t_blockDiscretisation), INTENT(IN) :: rsourceDiscr
  
  ! OPTIONAL: Whether the new discretisation structure should share its information
  ! with rsourceDiscr.
  ! =FALSE: Create a complete copy of rsourceDiscr which is independent
  !  of rsourceDiscr.
  ! =TRUE: The new discretisation will not be a complete new structure, but a
  !  'derived' structure, i.e. it uses the same dynamic information
  !  (handles and therefore element lists) as rsourceDiscr.
  ! If not specified, TRUE is assumed.
  LOGICAL, INTENT(IN), OPTIONAL :: bshare
!</input>
  
!<output>
  ! The new discretisation structure. Any old existing information in rdestDiscr
  ! is released if necessary.
  TYPE(t_blockDiscretisation), INTENT(INOUT), TARGET :: rdestDiscr
!</output>
  
!</subroutine>

    LOGICAL :: bshr
    INTEGER :: i
    
    bshr = .TRUE.
    IF (PRESENT(bshare)) bshr = bshare

    ! Release old information if present
    CALL spdiscr_releaseBlockDiscr(rdestDiscr,.TRUE.)

    ! Check that the source discretisation structure is valid.
    IF (rsourceDiscr%ndimension .LE. 0) THEN
      CALL output_line ('Source structure invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'spdiscr_duplicateBlockDiscr')  
      CALL sys_halt()
    END IF
    
    ! At first, derive a new block discretisation strucutr

    ! Copy all information from the source discretisation structure containing
    ! all basic information.
    CALL spdiscr_deriveBlockDiscr(rsourceDiscr,rdestDiscr)
    
    ! Concerning the substructures, at the moment we share the information.
    ! If bshare = false, we have to create copies.
    IF (.NOT. bshr) THEN
      DO i=1,rsourceDiscr%ncomponents
        CALL spdiscr_duplicateDiscrSc (rsourceDiscr%RspatialDiscr(i), &
            rdestDiscr%RspatialDiscr(i), .FALSE.)
      END DO
    END IF

  END SUBROUTINE  

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE spdiscr_infoBlockDiscr (rdiscr)
  
!<description>
    ! This subroutine outputs information about the block discretisation
!</description>

!<input>
    ! block discretisation
    TYPE(t_blockDiscretisation), INTENT(IN) :: rdiscr
!</input>
!</subroutine>

    ! local variables
    INTEGER :: icomponent

    CALL output_lbrk()
    CALL output_line ('BlockDiscretisation:')
    CALL output_line ('--------------------')
    CALL output_line ('ndimension:  '//TRIM(sys_siL(rdiscr%ndimension,1)))
    CALL output_line ('ccomplexity: '//TRIM(sys_siL(rdiscr%ccomplexity,1)))
    CALL output_line ('ncomponents: '//TRIM(sys_siL(rdiscr%ncomponents,3)))

    IF (ASSOCIATED(rdiscr%RspatialDiscr)) THEN
      DO icomponent=1,rdiscr%ncomponents
        CALL spdiscr_infoDiscr(rdiscr%RspatialDiscr(icomponent))
      END DO
    END IF

  END SUBROUTINE spdiscr_infoBlockDiscr

  ! ***************************************************************************

!<subroutine>

   SUBROUTINE spdiscr_infoDiscr (rspatialDiscr)

!<description>
     ! This subroutine outputs information about the spatial discretisation
!</description>

!<input>
     ! spatial discretisation
     TYPE(t_spatialDiscretisation), INTENT(IN) :: rspatialDiscr
!</input>
!</subroutine>

     ! local variable
     INTEGER :: inumFESpace

     CALL output_lbrk()
     CALL output_line ('SpatialDiscretisation:')
     CALL output_line ('----------------------')
     CALL output_line ('ndimension:             '&
         //TRIM(sys_siL(rspatialDiscr%ndimension,1)))
     CALL output_line ('bisCopy:                '&
         //TRIM(sys_sl(rspatialDiscr%bisCopy)))
     CALL output_line ('ccomplexity:            '&
         //TRIM(sys_siL(rspatialDiscr%ccomplexity,1)))
     CALL output_line ('inumFESpaces:           '&
         //TRIM(sys_siL(rspatialDiscr%inumFESpaces,15)))
     CALL output_line ('h_IelementDistr:        '&
         //TRIM(sys_siL(rspatialDiscr%h_IelementDistr,15)))
     CALL output_line ('h_IelementCounter:      '&
         //TRIM(sys_siL(rspatialDiscr%h_IelementCounter,15)))

     IF (ASSOCIATED(rspatialDiscr%RelementDistr)) THEN
       DO inumFESpace=1,rspatialDiscr%inumFESpaces
         CALL spdisc_infoElementDistr(rspatialDiscr%RelementDistr(inumFESpace))
       END DO
     END IF

   END SUBROUTINE spdiscr_infoDiscr

   ! ***************************************************************************

!<subroutine>
   
   SUBROUTINE spdisc_infoElementDistr (relementDistr)

!<description>
     ! This subroutine outputs information about the spatial discretisation
!</description>

!<input>
     ! element distribution
     TYPE(t_elementDistribution), INTENT(IN) :: relementDistr
!</input>
!</subroutine>

     CALL output_lbrk()
     CALL output_line ('ElementDistribution:')
     CALL output_line ('--------------------')
     CALL output_line ('ielement:        '//TRIM(sys_siL(relementDistr%celement,15)))
     CALL output_line ('ccubTypeBilForm: '//TRIM(sys_siL(relementDistr%ccubTypeBilForm,15)))
     CALL output_line ('ccubTypeLinForm: '//TRIM(sys_siL(relementDistr%ccubTypeLinForm,15)))
     CALL output_line ('ccubTypeEval:    '//TRIM(sys_siL(relementDistr%ccubTypeEval,15)))
     CALL output_line ('ctrafoType:      '//TRIM(sys_siL(relementDistr%ctrafoType,15)))
     CALL output_line ('NEL:             '//TRIM(sys_siL(relementDistr%NEL,15)))
     CALL output_line ('h_IelementList:  '//TRIM(sys_siL(relementDistr%h_IelementList,15)))

   END SUBROUTINE spdisc_infoElementDistr
END MODULE
