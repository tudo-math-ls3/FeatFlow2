!##############################################################################
!# ****************************************************************************
!# <name> groupfembase </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module contains a set of constant definitions and structures
!# which are used by the different group finite element modules.
!#
!# The following routines are available:
!#
!# 1.) gfem_initGroupFEMSet = gfem_initGFEMSetDirect /
!#                            gfem_initGFEMSetByMatrix
!#     -> Initialises a group finite element structure.
!#
!# 2.) gfem_initGroupFEMBlock
!#     -> Initialises a block of group finite element structures.
!#
!# 3.) gfem_releaseGroupFEMSet
!#     -> Releases a group finite element structure.
!#
!# 4.) gfem_releaseGroupFEMBlock
!#     -> Releases a block of group finite element structures.
!#
!# 5.) gfem_resizeGroupFEMSet = gfem_resizeGFEMSetDirect /
!#                              gfem_resizeGFEMSetByMatrix
!#     -> Resizes a group finite element structure.
!#
!# 6.) gfem_resizeGroupFEMBlock = gfem_resizeGFEMBlockDirect /
!#                                gfem_resizeGFEMBlockByMatrix
!#     -> Resizes a block of group finite element structures.
!#
!# 7.) gfem_copyGroupFEMSet
!#     -> Copies a group finite element structure to another structure
!#
!# 8.) gfem_copyGroupFEMBlock
!#     -> Copies a block of group finite element structures to another block
!#
!# 9.) gfem_duplicateGroupFEMSet
!#     -> Duplicates a group finite element structure
!#
!# 10.) gfem_duplicateGroupFEMBlock
!#      -> Duplicates a block of group finite element structures
!#
!# 11.)  gfem_initCoeffsFromMatrix
!#       -> Initialises precomputed coefficients from matrix
!#
!# 12.) gfem_isMatrixCompatible = gfem_isMatrixCompatibleSc /
!#                                gfem_isMatrixCompatibleBl
!#     -> Checks whether a matrix and a group finite element set are compatible
!#
!# 13.) gfem_isVectorCompatible = gfem_isVectorCompatibleSc /
!#                                gfem_isVectorCompatibleBl
!#     -> Checks whether a vector and a group finite element set are compatible
!#
!# 14.) gfem_getbase_IedgeListIdx
!#      -> Returns pointer to the index pointer for the edge list
!#
!# 15.) gfem_getbase_IedgeList
!#      -> Returns pointer to the edge list
!#
!# 16.) gfem_getbase_InodeListIdx
!#      -> Returns pointer to the index pointer for the node list
!#
!# 17.) gfem_getbase_InodeList
!#      -> Returns pointer to the node list
!#
!# 18.) gfem_getbase_IdiagList
!#      -> Returns pointer to the diagonal list
!#
!# 19.) gfem_getbase_DcoeffsAtNode / gfem_getbase_FcoeffsAtNode
!#      -> Returns pointer to the coefficients at nodes
!#
!# 20.) gfem_getbase_DcoeffsAtEdge / gfem_getbase_FcoeffsAtEdge
!#      -> Returns pointer to the coefficients at edges
!#
!# 21.) gfem_getbase_DcoeffsAtDiag / gfem_getbase_FcoeffsAtDiag
!#      -> Returns pointer to the coefficients at matrix diagonal
!#
!# 22.) gfem_getbase_IdofList
!#      -> Returns pointer to the list of restricted DOFs
!#
!# 23.) gfem_genNodeList
!#      -> Generates the node list for a given matrix
!#
!# 24.) gfem_genEdgeList
!#      -> Generates the edge list for a given matrix
!#
!# 25.) gfem_genDiagList
!#      -> Generates the diagonal pointer for a given matrix
!#
!# 26.) gfem_infoGroupFEMSet
!#      -> Outputs information about the group finite element set
!#
!# 27.) gfem_infoGroupFEMBlock
!#      -> Outputs information about the group finite element block
!#
!# The following auxiliary routines are available:
!#
!# 1.) gfem_copyH2D_IedgeList
!#      -> Copies the edge structure from the host memory
!#         to the device memory.
!#
!# 2.) gfem_copyD2H_IedgeList
!#      -> Copies the edge structure from the device memory
!#         to the host memory.
!#
!# 3.) gfem_copyH2D_CoeffsAtNode
!#     -> Copies the coefficients at nodes from the host memory
!#        to the device memory.
!#
!# 4.) gfem_copyD2H_CoeffsAtNode
!#     -> Copies the coefficients at nodes from the device memory
!#        to the host  memory.
!#
!# 5.) gfem_copyH2D_CoeffsAtEdge
!#     -> Copies the coefficients at edges from the host memory
!#        to the device memory.
!#
!# 6.) gfem_copyD2H_CoeffsAtEdge
!#     -> Copies the coefficients at matrix diagonals from the 
!#        device memory to the host  memory.
!#
!# 7.) gfem_copyH2D_CoeffsAtDiag
!#     -> Copies the coefficients at matrix diagonals from the
!#         host memory to the device memory.
!#
!# 8.) gfem_copyD2H_CoeffsAtDiag
!#     -> Copies the coefficients at edges from the device memory
!#        to the host  memory.
!#
!# 9.) gfem_allocCoeffs
!#     -> Initialises memory for precomputed coefficients
!#
!# </purpose>
!##############################################################################

module groupfembase

  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use storage

  implicit none

  private

  public :: t_groupFEMSet,t_groupFEMBlock
  public :: gfem_initGroupFEMSet
  public :: gfem_initGroupFEMBlock
  public :: gfem_releaseGroupFEMSet
  public :: gfem_releaseGroupFEMBlock
  public :: gfem_resizeGroupFEMSet
  public :: gfem_resizeGroupFEMBlock
  public :: gfem_copyGroupFEMSet
  public :: gfem_copyGroupFEMBlock
  public :: gfem_duplicateGroupFEMSet
  public :: gfem_duplicateGroupFEMBlock
  public :: gfem_initCoeffsFromMatrix
  public :: gfem_isMatrixCompatible
  public :: gfem_isVectorCompatible
  public :: gfem_getbase_IedgeListIdx
  public :: gfem_getbase_IedgeList
  public :: gfem_getbase_InodeListIdx
  public :: gfem_getbase_InodeList
  public :: gfem_getbase_IdiagList
  public :: gfem_getbase_DcoeffsAtNode
  public :: gfem_getbase_FcoeffsAtNode
  public :: gfem_getbase_DcoeffsAtEdge
  public :: gfem_getbase_FcoeffsAtEdge
  public :: gfem_getbase_DcoeffsAtDiag
  public :: gfem_getbase_FcoeffsAtDiag
  public :: gfem_getbase_IdofList
  public :: gfem_genNodeList
  public :: gfem_genEdgeList
  public :: gfem_genDiagList
  public :: gfem_infoGroupFEMSet
  public :: gfem_infoGroupFEMBlock

  public :: gfem_copyH2D_IedgeList
  public :: gfem_copyD2H_IedgeList
  public :: gfem_copyH2D_CoeffsAtNode
  public :: gfem_copyD2H_CoeffsAtNode
  public :: gfem_copyH2D_CoeffsAtEdge
  public :: gfem_copyD2H_CoeffsAtEdge
  public :: gfem_copyH2D_CoeffsAtDiag
  public :: gfem_copyD2H_CoeffsAtDiag
  public :: gfem_allocCoeffs

!<constants>
!<constantblock description="Global format flag for group FEM assembly">

  ! Assembly is undefined
  integer, parameter, public :: GFEM_UNDEFINED = 0

  ! Node-based assembly
  integer, parameter, public :: GFEM_NODEBASED = 1

  ! Edge-based assembly
  integer, parameter, public :: GFEM_EDGEBASED = 2
!</constantblock>


!<constantblock description="Global flag for initialisation of coefficients">

  ! Initialise coefficients for node-based assembly
  integer(I32), parameter, public :: GFEM_INITCOEFFS_NODEBASED = 2_I32**GFEM_NODEBASED

  ! Initialise coefficients for edge-based assembly
  integer(I32), parameter, public :: GFEM_INITCOEFFS_EDGEBASED = 2_I32**GFEM_EDGEBASED

  ! Initialise all coefficients
  integer(I32), parameter, public :: GFEM_INITCOEFFS_ALL       = GFEM_INITCOEFFS_NODEBASED+&
                                                                 GFEM_INITCOEFFS_EDGEBASED
!</constantblock>


!<constantblock description="Bitfield identifiers for properties of the \
!                            group finite element set">

  ! List of restricted DOFs has been computed: IdofList
  integer(I32), parameter, public :: GFEM_HAS_DOFLIST       = 2_I32**1

  ! Diagonal pointers have been computed: IdiagList
  integer(I32), parameter, public :: GFEM_HAS_DIAGLIST      = 2_I32**2

  ! Nodal structure has been computed: InodeList
  integer(I32), parameter, public :: GFEM_HAS_NODESTRUCTURE = 2_I32**3

  ! Edge-based structure has been computed: IedgeListIdx, IedgeList
  integer(I32), parameter, public :: GFEM_HAS_EDGESTRUCTURE = 2_I32**4

  ! Diagonal coefficient array has been computed: CoeffsAtDiag
  integer(I32), parameter, public :: GFEM_HAS_DIAGDATA      = 2_I32**5

  ! Nodal coefficient array has been computed: CoeffsAtNode
  integer(I32), parameter, public :: GFEM_HAS_NODEDATA      = 2_I32**6

  ! Edge-based coefficient array has been computed: CoeffsAtEdge
  integer(I32), parameter, public :: GFEM_HAS_EDGEDATA      = 2_I32**7  
!</constantblock>


!<constantblock description="Duplication flags. Specifies which information is duplicated \
!                            between group finite element structures">
  
  ! Duplicate atomic structure
  integer(I32), parameter, public :: GFEM_DUP_STRUCTURE     = 2_I32**16
  
  ! Duplicate list of restricted DOFs: IdofList
  integer(I32), parameter, public :: GFEM_DUP_DOFLIST       = GFEM_HAS_DOFLIST

  ! Duplicate diagonal pointers: IdiagList
  integer(I32), parameter, public :: GFEM_DUP_DIAGLIST      = GFEM_HAS_DIAGLIST
  
  ! Duplicate edge-nodal structure: InodeList
  integer(I32), parameter, public :: GFEM_DUP_NODESTRUCTURE = GFEM_HAS_NODESTRUCTURE

  ! Duplicate edge-based structure: IedgeListIdx, IedgeList
  integer(I32), parameter, public :: GFEM_DUP_EDGESTRUCTURE = GFEM_HAS_EDGESTRUCTURE

  ! Duplicate diagonal coefficient array: CoeffsAtDiag
  integer(I32), parameter, public :: GFEM_DUP_DIAGDATA      = GFEM_HAS_DIAGDATA

  ! Duplicate nodal coefficient array: CoeffsAtNode
  integer(I32), parameter, public :: GFEM_DUP_NODEDATA      = GFEM_HAS_NODEDATA

  ! Duplicate edge-based coefficient array: CoeffsAtEdge
  integer(I32), parameter, public :: GFEM_DUP_EDGEDATA      = GFEM_HAS_EDGEDATA

  ! Duplicate everything
  integer(I32), parameter, public :: GFEM_DUP_ALL           = GFEM_DUP_STRUCTURE+&
                                                              GFEM_DUP_DOFLIST+&
                                                              GFEM_DUP_DIAGLIST+&
                                                              GFEM_DUP_NODESTRUCTURE+&
                                                              GFEM_DUP_EDGESTRUCTURE+&
                                                              GFEM_DUP_DIAGDATA+&
                                                              GFEM_DUP_NODEDATA+&
                                                              GFEM_DUP_EDGEDATA
!</constantblock>


!<constantblock description="Duplication flags. Specifies which information is shared \
!                            between group finite element structures">
  
  ! Share list of restricted DOFs: IdofList
  integer(I32), parameter, public :: GFEM_SHARE_DOFLIST       = GFEM_HAS_DOFLIST

  ! Share diagonal pointers: IdiagList
  integer(I32), parameter, public :: GFEM_SHARE_DIAGLIST      = GFEM_HAS_DIAGLIST

  ! Share edge-nodal structure: InodeList
  integer(I32), parameter, public :: GFEM_SHARE_NODESTRUCTURE = GFEM_HAS_NODESTRUCTURE

  ! Share edge-based structure: IedgeListIdx, IedgeList
  integer(I32), parameter, public :: GFEM_SHARE_EDGESTRUCTURE = GFEM_HAS_EDGESTRUCTURE

  ! Share diagonal coefficient array: CoeffsAtDiag
  integer(I32), parameter, public :: GFEM_SHARE_DIAGDATA      = GFEM_HAS_DIAGDATA

  ! Share nodal coefficient array: CoeffsAtNode
  integer(I32), parameter, public :: GFEM_SHARE_NODEDATA      = GFEM_HAS_NODEDATA

  ! Share edge-based coefficient array: CoeffsAtEdge
  integer(I32), parameter, public :: GFEM_SHARE_EDGEDATA      = GFEM_HAS_EDGEDATA 
!</constantblock>
!</constants>


!<types>
!<typeblock>

  ! This structure holds all required information to evaluate a set of
  ! degrees of freedom by the group finite element formulation.

  type t_groupFEMSet

    ! Format Tag: Identifies the type of assembly
    integer :: cassemblyType = GFEM_UNDEFINED

    ! Format Tag: Identifies the data type
    integer :: cdataType = ST_DOUBLE

    ! Duplication Flag: This is a bitfield coming from an OR
    ! combination of different GFEM_SHARE_xxxx constants and
    ! specifies which parts of the structure are shared with
    ! another structure.
    integer(I32) :: iduplicationFlag = GFEM_UNDEFINED

    ! Specification Flag: Specifies the group finite element set. This
    ! is a bitfield coming from an OR combination of different
    ! GFEM_HAS_xxxx constants and specifies various properties of the
    ! structure.
    integer(I32) :: isetSpec = GFEM_UNDEFINED
    
    ! Number of non-zero entries of the sparsity pattern
    integer :: NA = 0

    ! Number of equations of the sparsity pattern
    integer :: NEQ = 0

    ! Number of edges of the sparsity pattern
    integer :: NEDGE = 0

    ! Number of local variables; in general scalar solution vectors of
    ! size NEQ posses NEQ entries. However, scalar vectors can be
    ! interleaved, that is, each of the NEQ entries stores NVAR local
    ! variables. In this case, NEQ remains unmodified but NVAR>1 such
    ! that the physical length of the vector is NEQ*NVAR.
    integer :: NVAR = 1

    ! Number of precomputed coefficients stored in CoeffsAtDiag.
    integer :: ncoeffsAtDiag = 0

    ! Number of precomputed coefficients stored in CoeffsAtNode.
    integer :: ncoeffsAtNode = 0

    ! Number of precomputed coefficients stored in CoeffsAtEdge.
    integer :: ncoeffsAtEdge = 0

    ! Handle to list of DOFs to which the global matrix is restricted
    integer :: h_IdofList = ST_NOHANDLE

    ! Handle to diagonal positions of the underlying matrix
    ! IdiagList(1,1:NEQ) : the node number i of the equation ieq
    ! IdiagList(2,1:NEQ) : the position of the diagoanl entry ii that
    !                      corresponds to the equation ieq
    integer :: h_IdiagList = ST_NOHANDLE

    ! Handle to index pointer for edge structure
    ! The numbers IedgeListIdx(k):IedgeListIdx(k+1)-1
    ! denote the edge numbers of the k-th group of edges.
    integer :: h_IedgeListIdx = ST_NOHANDLE

    ! Handle to edge structure
    ! IedgeList(1:2,1:NEDGE) : the two end-points i and j of the edge (ij)
    ! IedgeList(3:4,1:NEDGE) : the two matrix position ij and ji that
    !                          correspond to the edge (ij)
    ! IedgeList(5:6,1:NEDGE) : the two matrix position ii and jj that
    !                          correspond to the diagonal entries
    integer :: h_IedgeList = ST_NOHANDLE
    
    ! Handle to index pointer for node structure
    ! InodeListIdx(1:NEQ+1) : the index separator of the node list
    integer :: h_InodeListIdx = ST_NOHANDLE

    ! Handle to nodal structure
    ! InodeList(1,1:NA) : the global number of the node
    ! InodeList(2,1:NA) : the global position of the corresponding matrix entry
    integer :: h_InodeList = ST_NOHANDLE
    
    ! Handle to precomputed coefficiets at matrix diagonal
    integer :: h_CoeffsAtDiag = ST_NOHANDLE
    
    ! Handle to precomputed coefficients at edges
    integer :: h_CoeffsAtEdge = ST_NOHANDLE
    
    ! Handle to precomputed coefficients at nodes
    integer :: h_CoeffsAtNode = ST_NOHANDLE
    
    ! Internal array used for consistency checks
    ! Iconsistency(1) : NA of the matrix
    ! Iconsistency(2) : NEQ of the matrix
    ! Iconsistency(3) : NCOLS of the matrix
    integer, dimension(3) :: Iconsistency
  end type t_groupFEMSet
!</typeblock>


!<typeblock>

  ! This structure holds multiple sets of degrees of freedom
  
  type t_groupFEMBlock

    ! Number of blocks
    integer :: nblocks = 0

    ! A 1D array with sets of degrees of freedom
    type(t_groupFEMSet), dimension(:), pointer :: RgroupFEMBlock => null() 

  end type t_groupFEMBlock
!</typeblock>
!</types>


  interface gfem_initGroupFEMSet
    module procedure gfem_initGFEMSetDirect
    module procedure gfem_initGFEMSetByMatrix
  end interface gfem_initGroupFEMSet

  interface gfem_resizeGroupFEMSet
    module procedure gfem_resizeGFEMSetDirect
    module procedure gfem_resizeGFEMSetByMatrix
  end interface

  interface gfem_resizeGroupFEMBlock
    module procedure gfem_resizeGFEMBlockDirect
    module procedure gfem_resizeGFEMBlockByMatrix
  end interface

  interface gfem_isMatrixCompatible
    module procedure gfem_isMatrixCompatibleSc
    module procedure gfem_isMatrixCompatibleBl
  end interface

  interface gfem_isVectorCompatible
    module procedure gfem_isVectorCompatibleSc
    module procedure gfem_isVectorCompatibleBl
  end interface

  interface gfem_getbase_InodeList
    module procedure gfem_getbase_InodeList1D
    module procedure gfem_getbase_InodeList2D
  end interface

contains

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_initGFEMSetDirect(rgroupFEMSet, NA, NEQ, NEDGE, NVAR,&
      ncoeffsAtDiag, ncoeffsAtNode, ncoeffsAtEdge, cassembly, cdataType)

!<description>
    ! This subroutine initialises a set of degrees of freedom for
    ! using the group finite element formulation directly. By setting
    ! some parameters to zero, some internla data structures are not
    ! generated which can be postponed to other subroutine calls.
!</description>

!<input>
    ! Number of non-zero entries
    integer, intent(in) :: NA

    ! Number of equations
    integer, intent(in) :: NEQ

    ! Number of edges
    integer, intent(in) :: NEDGE
    
    ! Number of local variables
    integer, intent(in) :: NVAR

    ! Number of precomputed coefficients
    integer, intent(in) :: ncoeffsAtDiag
    integer, intent(in) :: ncoeffsAtNode
    integer, intent(in) :: ncoeffsAtEdge

    ! Type of assembly
    integer, intent(in) :: cassembly

    ! OPTIONAL: data type
    ! If not present, ST_DOUBLE will be used
    integer, intent(in), optional :: cdataType
!</input>

!<output>
    ! Group finite element set
    type(t_groupFEMSet), intent(out) :: rgroupFEMSet
!</output>
!</subroutine>

    ! Set data
    rgroupFEMSet%NA            = NA
    rgroupFEMSet%NEQ           = NEQ
    rgroupFEMSet%NVAR          = NVAR
    rgroupFEMSet%NEDGE         = NEDGE
    rgroupFEMSet%cassemblyType = cassembly

    if (present(cdataType))&
        rgroupFEMSet%cdataType = cdataType

    ! Allocate memory for precomputed coefficients; if some of the
    ! necessary data is not available, no memory will be allocated
    call gfem_allocCoeffs(rgroupFEMSet,&
        ncoeffsAtDiag, ncoeffsAtNode, ncoeffsAtEdge)

  end subroutine gfem_initGFEMSetDirect

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_initGFEMSetByMatrix(rgroupFEMSet, rmatrix, ncoeffsAtDiag,&
      ncoeffsAtNode, ncoeffsAtEdge, cassembly, cdataType, IdofList)

!<description>
    ! This subroutine initialises a set of degrees of freedom for
    ! using the group finite element formulation indirectly by
    ! deriving all information from the scalar template matrix.
    ! If IdofList is given, then its entries are used as degrees
    ! of freedom to which the group finite element set should be
    ! restricted to.
!</description>

!<input>
    ! Scalar template matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! Number of precomputed coefficients
    integer, intent(in) :: ncoeffsAtDiag
    integer, intent(in) :: ncoeffsAtNode
    integer, intent(in) :: ncoeffsAtEdge

    ! Type of assembly
    integer, intent(in) :: cassembly

    ! OPTIONAL: data type
    ! If not present, ST_DOUBLE will be used
    integer, intent(in), optional :: cdataType

    ! OPTIONAL: a list of degress of freedoms to which the
    ! group finite element set should be restricted.
    integer, dimension(:), intent(in), optional :: IdofList
!</input>

!<output>
    ! group finite element set
    type(t_groupFEMSet), intent(out) :: rgroupFEMSet
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdofList
    integer :: na,neq,ncols,nedge
    
    ! Calculate dimensions of matrix (possibly with restriction)
    call lsyssc_calcDimsFromMatrix(rmatrix, na, neq, ncols, nedge, IdofList)
    
    ! Call initialisation routine
    call gfem_initGroupFEMSet(rgroupFEMSet, na, neq, nedge, rmatrix%NVAR,&
        ncoeffsAtDiag, ncoeffsAtNode, ncoeffsAtEdge, cassembly, cdataType)

    ! Initialise array for consistency checks
    rgroupFEMSet%Iconsistency(1) = rmatrix%NA
    rgroupFEMSet%Iconsistency(2) = rmatrix%NEQ
    rgroupFEMSet%Iconsistency(3) = rmatrix%NCOLS
    
    ! Set handle to restricted DOF list
    if (present(IdofList)) then
      call storage_new('gfem_initGFEMSetByMatrix', 'IdofList',&
          size(IdofList), ST_INT, rgroupFEMSet%h_IdofList,&
          ST_NEWBLOCK_NOINIT)
      
      ! Copy list of restricted DOFs
      call gfem_getbase_IdofList(rgroupFEMSet, p_IdofList)
      call lalg_copyVector(IdofList, p_IdofList)
      
      ! Set specifier
      rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST)
    end if
    
  end subroutine gfem_initGFEMSetByMatrix

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_initGroupFEMBlock(rgroupFEMBlock, nblocks)

!<description>
    ! This subroutine initialises a block of group finite element sets.
    ! If rgroupFEMBlock has already some blocks associated, then the
    ! array is re-allocated internally and data is copied accordingly.
!</description>

!<input>
    ! Number of blocks
    integer, intent(in) :: nblocks
!</input>

!<inputoutput>
    ! Block of group finite element sets
    type(t_groupFEMBlock), intent(inout) :: rgroupFEMBlock
!</inputoutput>

    ! local variables
    type(t_groupFEMSet), dimension(:), allocatable :: p_RgroupFEMBlock
    integer :: i

    if (rgroupFEMBlock%nblocks .eq. 0) then
      rgroupFEMBlock%nblocks = nblocks
      allocate(rgroupFEMBlock%RgroupFEMBlock(nblocks))
    else
      ! Make backup of group finite element sets
      allocate(p_RgroupFEMBlock(rgroupFEMBlock%nblocks))
      do i = 1, rgroupFEMBlock%nblocks
        p_RgroupFEMBlock(i) = rgroupFEMBlock%RgroupFEMBlock(i)
      end do

      ! Re-allocate memory
      deallocate(rgroupFEMBlock%RgroupFEMBlock)
      allocate(rgroupFEMBlock%RgroupFEMBlock(nblocks))
      
      ! Copy data from backup
      do i = 1, min(nblocks,rgroupFEMBlock%nblocks)
        rgroupFEMBlock%RgroupFEMBlock(i) = p_RgroupFEMBlock(i)
      end do

      ! Deallocate backup
      deallocate(p_RgroupFEMBlock)

      ! Set new blocksize
      rgroupFEMBlock%nblocks = nblocks
    end if

  end subroutine gfem_initGroupFEMBlock

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_releaseGroupFEMSet(rgroupFEMSet)

!<description>
    ! This subroutine releases a group finite element set
!</description>

!<inputoutput>
    ! group finite element set
    type(t_groupFEMSet), intent(inout) :: rgroupFEMSet
!</inputoutput>
!</subroutine>

    ! Release list of restricted DOFs
    if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_DOFLIST)) then
      if (rgroupFEMSet%h_IdofList .ne. ST_NOHANDLE)&
          call storage_free(rgroupFEMSet%h_IdofList)
    end if
    rgroupFEMSet%h_IdofList = ST_NOHANDLE

    ! Release diagonal structure
    if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_DIAGLIST)) then
      if (rgroupFEMSet%h_IdiagList .ne. ST_NOHANDLE)&
          call storage_free(rgroupFEMSet%h_IdiagList)
    end if
    rgroupFEMSet%h_IdiagList = ST_NOHANDLE

    ! Release nodal structure
    if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_NODESTRUCTURE)) then
      if (rgroupFEMSet%h_InodeList .ne. ST_NOHANDLE)&
          call storage_free(rgroupFEMSet%h_InodeList)
      if (rgroupFEMSet%h_InodeListIdx .ne. ST_NOHANDLE)&
          call storage_free(rgroupFEMSet%h_InodeListIdx)
    end if
    rgroupFEMSet%h_InodeList = ST_NOHANDLE
    rgroupFEMSet%h_InodeListIdx = ST_NOHANDLE

    ! Release edge structure
    if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_EDGESTRUCTURE)) then
      if (rgroupFEMSet%h_IedgeList .ne. ST_NOHANDLE)&
          call storage_free(rgroupFEMSet%h_IedgeList)
      if (rgroupFEMSet%h_IedgeListIdx .ne. ST_NOHANDLE)&
          call storage_free(rgroupFEMSet%h_IedgeListIdx)
    end if
    rgroupFEMSet%h_IedgeList = ST_NOHANDLE
    rgroupFEMSet%h_IedgeListIdx = ST_NOHANDLE

    ! Release precomputed coefficients at matrix diagonals
    if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_DIAGDATA)) then
      if (rgroupFEMSet%h_CoeffsAtDiag .ne. ST_NOHANDLE)&
          call storage_free(rgroupFEMSet%h_CoeffsAtDiag)
    end if
    rgroupFEMSet%h_CoeffsAtDiag = ST_NOHANDLE

    ! Release precomputed coefficients at nodes
    if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_NODEDATA)) then
      if (rgroupFEMSet%h_CoeffsAtNode .ne. ST_NOHANDLE)&
          call storage_free(rgroupFEMSet%h_CoeffsAtNode)
    end if
    rgroupFEMSet%h_CoeffsAtNode = ST_NOHANDLE

    ! Release precomputed coefficients at edges
    if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_EDGEDATA)) then
      if (rgroupFEMSet%h_CoeffsAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rgroupFEMSet%h_CoeffsAtEdge)
    end if
    rgroupFEMSet%h_CoeffsAtEdge = ST_NOHANDLE

    ! Unset ownerships and specification
    rgroupFEMSet%iduplicationFlag = GFEM_UNDEFINED
    rgroupFEMSet%isetSpec         = GFEM_UNDEFINED

    ! Reset data
    rgroupFEMSet%cassemblyType = GFEM_UNDEFINED
    rgroupFEMSet%cdataType     = ST_DOUBLE
    rgroupFEMSet%NA    = 0
    rgroupFEMSet%NEQ   = 0
    rgroupFEMSet%NEDGE = 0
    rgroupFEMSet%NVAR  = 1
    rgroupFEMSet%ncoeffsAtDiag = 0
    rgroupFEMSet%ncoeffsAtNode = 0
    rgroupFEMSet%ncoeffsAtEdge = 0

    ! Reset array for consistency checks
    rgroupFEMSet%Iconsistency = 0

  contains

    !**************************************************************
    ! Checks if bitfield ibitfield in iflag is not set.
    
    pure function check(iflag, ibitfield)
      
      integer(I32), intent(in) :: iflag,ibitfield
      
      logical :: check
      
      check = (iand(iflag,ibitfield) .ne. ibitfield)
      
    end function check
    
  end subroutine gfem_releaseGroupFEMSet

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_releaseGroupFEMBlock(rgroupFEMBlock)

!<description>
    ! This subroutine releases a block of group finite element sets
!</description>

!<inputoutput>
    ! group finite element block
    type(t_groupFEMBlock), intent(inout) :: rgroupFEMBlock
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    if (rgroupFEMBlock%nblocks > 0) then
      do i = 1, rgroupFEMBlock%nblocks
        call gfem_releaseGroupFEMSet(rgroupFEMBlock%RgroupFEMBlock(i))
      end do
      
      deallocate(rgroupFEMBlock%RgroupFEMBlock)
    end if
    
    rgroupFEMBlock%nblocks = 0

  end subroutine gfem_releaseGroupFEMBlock

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_resizeGFEMSetDirect(rgroupFEMSet, NA, NEQ, NEDGE)

!<description>
    ! This subroutine resizes a set of degrees of freedom for
    ! using the group finite element formulation directly.
    ! Note that shared data cannot be resized by this routine!
!</description>

!<input>
    ! Number of non-zero entries
    integer, intent(in) :: NA

    ! Number of equations
    integer, intent(in) :: NEQ

    ! Number of edges
    integer, intent(in) :: NEDGE
!</input>

!<inputoutput>
    ! group finite element set
    type(t_groupFEMSet), intent(inout) :: rgroupFEMSet
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(3) :: Isize3D
    integer, dimension(2) :: Isize2D
    integer :: isize,idimension

    !---------------------------------------------------------------------------
    ! Resize non-zero entries
    !---------------------------------------------------------------------------
    if (rgroupFEMSet%NA .ne. NA) then

      ! Set new number of nonzero entries
      rgroupFEMSet%NA = NA
      
      !-----------------------------------------------------------------------
      if (rgroupFEMSet%h_InodeList .ne. ST_NOHANDLE) then
        if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_NODESTRUCTURE)) then
          call storage_getdimension(rgroupFEMSet%h_InodeList, idimension)
          if (idimension .eq. 1) then
            call storage_getsize(rgroupFEMSet%h_InodeList, isize)
            if (rgroupFEMSet%NA .ne. isize) then
              call output_line('Handle h_InodeList '//&
                  'is shared and cannot be resized!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfem_resizeGFEMSetDirect')
              call sys_halt()
            end if
          else
            call storage_getsize(rgroupFEMSet%h_InodeList, Isize2D)
            if (rgroupFEMSet%NA .ne. Isize2D(2)) then
              call output_line('Handle h_InodeList '//&
                  'is shared and cannot be resized!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfem_resizeGFEMSetDirect')
              call sys_halt()
            end if
          end if
        else
          call storage_realloc('gfem_resizeGFEMSetDirect',&
              rgroupFEMSet%NA, rgroupFEMSet%h_InodeList,&
              ST_NEWBLOCK_NOINIT, .false.)
          
          ! Reset specifier
          rgroupFEMSet%isetSpec = iand(rgroupFEMSet%isetSpec,&
                                       not(GFEM_HAS_NODESTRUCTURE))
        end if
      end if

      !-------------------------------------------------------------------------
      if (rgroupFEMSet%h_CoeffsAtNode .ne. ST_NOHANDLE) then
        if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_NODEDATA)) then
          call storage_getsize(rgroupFEMSet%h_CoeffsAtNode, Isize2D)
          if (rgroupFEMSet%NA .ne. Isize2D(2)) then
            call output_line('Handle h_CoeffsAtNode '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_resizeGFEMSetDirect')
            call sys_halt()
          end if
        else
          call storage_realloc('gfem_resizeGFEMSetDirect',&
              rgroupFEMSet%NA, rgroupFEMSet%h_CoeffsAtNode,&
              ST_NEWBLOCK_NOINIT, .false.)
          
          ! Reset specifier
          rgroupFEMSet%isetSpec = iand(rgroupFEMSet%isetSpec,&
                                       not(GFEM_HAS_NODEDATA))
        end if
      end if
      
    end if

    !---------------------------------------------------------------------------
    ! Resize nodal quantities
    !---------------------------------------------------------------------------
    if (rgroupFEMSet%NEQ .ne. NEQ) then
      
      ! Set new number of nodes
      rgroupFEMSet%NEQ = NEQ

      !-----------------------------------------------------------------------
      if (rgroupFEMSet%h_IdiagList .ne. ST_NOHANDLE) then
        if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_DIAGLIST)) then
          call storage_getsize(rgroupFEMSet%h_IdiagList, Isize2D)
          if (rgroupFEMSet%NEQ .ne. Isize2D(2)) then
            call output_line('Handle h_IdiagList '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_resizeGFEMSetDirect')
            call sys_halt()
          end if
        else
          ! Resize array
          call storage_realloc('gfem_resizeGFEMSetDirect',&
              rgroupFEMSet%NEQ, rgroupFEMSet%h_IdiagList,&
              ST_NEWBLOCK_NOINIT, .false.)
          
          ! Reset specifier
          rgroupFEMSet%isetSpec = iand(rgroupFEMSet%isetSpec,&
                                       not(GFEM_HAS_DIAGLIST))
        end if
      end if

      !-----------------------------------------------------------------------
      if (rgroupFEMSet%h_InodeListIdx .ne. ST_NOHANDLE) then
        if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_NODESTRUCTURE)) then
          call storage_getsize(rgroupFEMSet%h_InodeListIdx, isize)
          if (rgroupFEMSet%NEQ+1 .ne. isize) then
            call output_line('Handle h_InodeListIdx '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_resizeGFEMSetDirect')
            call sys_halt()
          end if
        else
          ! Resize array
          call storage_realloc('gfem_resizeGFEMSetDirect',&
              rgroupFEMSet%NEQ+1, rgroupFEMSet%h_InodeListIdx,&
              ST_NEWBLOCK_NOINIT, .false.)
          
          ! Reset specifier
          rgroupFEMSet%isetSpec = iand(rgroupFEMSet%isetSpec,&
                                       not(GFEM_HAS_NODESTRUCTURE))
        end if
      end if

      !-------------------------------------------------------------------------
      if (rgroupFEMSet%h_CoeffsAtDiag .ne. ST_NOHANDLE) then
        if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_DIAGDATA)) then
          call storage_getsize(rgroupFEMSet%h_CoeffsAtDiag, Isize2D)
          if (rgroupFEMSet%NEQ .ne. Isize2D(2)) then
            call output_line('Handle h_CoeffsAtDiag '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_resizeGFEMSetDirect')
            call sys_halt()
          end if
        else
          call storage_realloc('gfem_resizeGFEMSetDirect',&
              rgroupFEMSet%NEQ, rgroupFEMSet%h_CoeffsAtDiag,&
              ST_NEWBLOCK_NOINIT, .false.)
          
          ! Reset specifier
          rgroupFEMSet%isetSpec = iand(rgroupFEMSet%isetSpec,&
                                       not(GFEM_HAS_DIAGDATA))
        end if
      end if

    end if

    !---------------------------------------------------------------------------
    ! Resize edge-based quantities
    !---------------------------------------------------------------------------
    if (rgroupFEMSet%NEDGE .ne. NEDGE) then
      
      ! Set new number of edges
      rgroupFEMSet%NEDGE = NEDGE

      !-------------------------------------------------------------------------
      if (rgroupFEMSet%h_IedgeList .ne. ST_NOHANDLE) then
        if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_EDGESTRUCTURE)) then
          call storage_getsize(rgroupFEMSet%h_IedgeList, Isize2D)
          if (rgroupFEMSet%NEDGE .ne. Isize2D(2)) then
            call output_line('Handle h_IedgeList '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_resizeGFEMSetDirect')
            call sys_halt()
          end if
        else
          ! Resize array
          call storage_realloc('gfem_resizeGFEMSetDirect',&
              rgroupFEMSet%NEDGE, rgroupFEMSet%h_IedgeList,&
              ST_NEWBLOCK_NOINIT, .false.)
          
          ! Reset specifiert
          rgroupFEMSet%isetSpec = iand(rgroupFEMSet%isetSpec,&
                                       not(GFEM_HAS_EDGESTRUCTURE))
        end if
      end if

      !-----------------------------------------------------------------------
      if (rgroupFEMSet%h_CoeffsAtEdge .ne. ST_NOHANDLE) then
        if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_EDGEDATA)) then
          call storage_getsize(rgroupFEMSet%h_CoeffsAtEdge, Isize3D)
          if (rgroupFEMSet%NEDGE .ne. Isize3D(3)) then
            call output_line('Handle h_CoeffsAtEdge '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_resizeGFEMSetDirect')
            call sys_halt()
          end if
        else
          call storage_realloc('gfem_resizeGFEMSetDirect',&
              rgroupFEMSet%NEDGE, rgroupFEMSet%h_CoeffsAtEdge,&
              ST_NEWBLOCK_NOINIT, .false.)
          
          ! Reset specifier
          rgroupFEMSet%isetSpec = iand(rgroupFEMSet%isetSpec,&
                                       not(GFEM_HAS_EDGEDATA))
        end if
      end if

    end if

  contains

    !**************************************************************
    ! Checks if bitfield ibitfield in iflag is set.
    
    pure function check(iflag, ibitfield)

      integer(I32), intent(in) :: iflag,ibitfield
      
      logical :: check
      
      check = (iand(iflag,ibitfield) .eq. ibitfield)

    end function check

  end subroutine gfem_resizeGFEMSetDirect
  
  ! ***************************************************************************

!<subroutine>

  subroutine gfem_resizeGFEMSetByMatrix(rgroupFEMSet, rmatrix, IdofList)

!<description>
    ! This subroutine resizes a set of degrees of freedom for
    ! using the group finite element formulation indirectly by
    ! deriving all information from the scalar template matrix.
    ! If IdofList is given, then its entries are used as degrees
    ! of freedom to which the group finite elment set should be
    ! restricted to. Note that shared data cannot be resized!
    ! Moreover, specifying IdofList is only admissible if the
    ! group finite element set rgroupFEMSet does not share a list
    ! of restricted DOFs.
!</description>

!<input>
    ! Scalar template matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: a list of degress of freedoms to which the
    ! group finite element set should be restricted.
    integer, dimension(:), intent(in), optional, target :: IdofList
!</input>

!<inputoutput>
    ! group finite element set
    type(t_groupFEMSet), intent(inout) :: rgroupFEMSet
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdofList
    integer :: na,neq,ncols,nedge

    ! Do we have a list of restricted DOFs?
    if (present(IdofList)) then
      ! Check if list of restricted DOFs is shared with another set
      if (check(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_DOFLIST)) then
        call output_line('List of restricted DOFs is not owned by structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_resizeGFEMSetByMatrix')
        call sys_halt()
      else   
        ! Calculate dimensions of matrix with restriction
        call lsyssc_calcDimsFromMatrix(rmatrix, na, neq, ncols, nedge, IdofList)

        ! Set handle to restricted DOF list
        if (rgroupFEMSet%h_IdofList .eq. ST_NOHANDLE) then
          call storage_new('gfem_resizeGFEMSetByMatrix', 'IdofList',&
              size(IdofList), ST_INT, rgroupFEMSet%h_IdofList,&
              ST_NEWBLOCK_NOINIT)
        else
          call storage_realloc('gfem_resizeGFEMSetByMatrix',&
              size(IdofList), rgroupFEMSet%h_IdofList,&
              ST_NEWBLOCK_NOINIT, .false.)
        end if
        
        ! Copy list of restricted DOFs
        call gfem_getbase_IdofList(rgroupFEMSet, p_IdofList)
        call lalg_copyVector(IdofList, p_IdofList)
        
        ! Set specifier
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST)
      end if
      
    elseif (check(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST)) then
      ! Calculate dimensions of matrix with internally defined restriction
      call gfem_getbase_IdofList(rgroupFEMSet, p_IdofList)
      call lsyssc_calcDimsFromMatrix(rmatrix, na, neq, ncols, nedge, p_IdofList)
    else
      ! Calculate dimensions of matrix without restriction
      call lsyssc_calcDimsFromMatrix(rmatrix, na, neq, ncols, nedge)
    end if
    
    ! Call direct resize routine
    call gfem_resizeGroupFEMSet(rgroupFEMSet, na, neq, nedge)

    ! Initialise array for consistency checks
    rgroupFEMSet%Iconsistency(1) = rmatrix%NA
    rgroupFEMSet%Iconsistency(2) = rmatrix%NEQ
    rgroupFEMSet%Iconsistency(3) = rmatrix%NCOLS
    
  contains

    !**************************************************************
    ! Checks if bitfield ibitfield in iflag is set.
    
    pure function check(iflag, ibitfield)

      integer(I32), intent(in) :: iflag,ibitfield
      
      logical :: check
      
      check = (iand(iflag,ibitfield) .eq. ibitfield)

    end function check
    
  end subroutine gfem_resizeGFEMSetByMatrix

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_resizeGFEMBlockDirect(rgroupFEMBlock, NA, NEQ, NEDGE)

!<description>
    ! This subroutine resizes block of sets of degrees of freedom for
    ! using the group finite element formulation directly.
!</description>

!<input>
    ! Number of non-zero entries
    integer, intent(in) :: NA

    ! Number of equations
    integer, intent(in) :: NEQ

    ! Number of edges
    integer, intent(in) :: NEDGE   
!</input>

!<inputoutput>
    ! Block of group finite element sets
    type(t_groupFEMBlock), intent(inout) :: rgroupFEMBlock
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i
    
    do i = 1, rgroupFEMBlock%nblocks
      call gfem_resizeGroupFEMSet(rgroupFEMBlock%RgroupFEMBlock(i), NA, NEQ, NEDGE)
    end do
  
  end subroutine gfem_resizeGFEMBlockDirect

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_resizeGFEMBlockByMatrix(rgroupFEMBlock, rmatrix, IdofList)

!<description>
    ! This subroutine resizes a block of sets of degrees of freedom
    ! for using the group finite element formulation indirectly by
    ! deriving all information from the scalar template matrix.
    ! If IdofList is given, then its entries are used as degrees
    ! of freedom to which the group finite elment set should be
    ! restricted to. Note that IdofList is only admissible if the
    ! blocks of group finite element sets rgroupFEMBlock do not
    ! share a list of restricted DOFs.
!</description>

!<input>
    ! Scalar template matrix
    type(t_matrixScalar), intent(in) :: rmatrix
    
    ! OPTIONAL: a list of degress of freedoms to which the
    ! group finite element set should be restricted.
    integer, dimension(:), intent(in), optional :: IdofList
!</input>

!<inputoutput>
    ! Block of group finite element sets
    type(t_groupFEMBlock), intent(inout) :: rgroupFEMBlock
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i
    
    do i = 1, rgroupFEMBlock%nblocks
      call gfem_resizeGroupFEMSet(rgroupFEMBlock%RgroupFEMBlock(i),&
          rmatrix, IdofList)
    end do

  end subroutine gfem_resizeGFEMBlockByMatrix

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_copyGroupFEMSet(rgroupFEMSetSrc, rgroupFEMSetDest, idupFlag,&
      bpreserveContent)

!<description>
    ! This subroutine selectively copies data from the source
    ! structure rgroupFEMSetSrc to the destination structure
    ! rgroupFEMSetDest. If the optional flag bpreserveContent
    ! is set to TRUE, then the existing content is not removed
    ! before copying some content from the source structure.
!<description>

!<input>
    ! Source group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSetSrc

    ! Duplication flag that decides on how to set up the structure
    integer(I32), intent(in) :: idupFlag

    ! OPTIONAL: Flag to force that existing content is preserved
    logical, intent(in), optional :: bpreserveContent
!</input>

!<inputoutput>
    ! Destination group finite element set
    type(t_groupFEMSet), intent(inout) :: rgroupFEMSetDest
!</inputoutput>
!</subroutine>

    ! local variables
    logical :: bclear

    ! Clear destination structre unless preservation of existing
    ! content is explicitly enforced by the calling routine
    bclear = .true.
    if (present(bpreserveContent)) bclear=.not.bpreserveContent
    if (bclear) call gfem_releaseGroupFEMSet(rgroupFEMSetDest)

    ! Copy structural data
    if (check(idupFlag, GFEM_DUP_STRUCTURE)) then
      rgroupFEMSetDest%cassemblyType = rgroupFEMSetSrc%cassemblyType
      rgroupFEMSetDest%cdataType     = rgroupFEMSetSrc%cdataType
      rgroupFEMSetDest%NA            = rgroupFEMSetSrc%NA
      rgroupFEMSetDest%NEQ           = rgroupFEMSetSrc%NEQ
      rgroupFEMSetDest%NEDGE         = rgroupFEMSetSrc%NEDGE
      rgroupFEMSetDest%NVAR          = rgroupFEMSetSrc%NVAR
      rgroupFEMSetDest%ncoeffsAtDiag = rgroupFEMSetSrc%ncoeffsAtDiag
      rgroupFEMSetDest%ncoeffsAtNode = rgroupFEMSetSrc%ncoeffsAtNode
      rgroupFEMSetDest%ncoeffsAtEdge = rgroupFEMSetSrc%ncoeffsAtEdge
    end if
    

    ! Copy list of restricted DOFs
    if (check(idupFlag, GFEM_DUP_DOFLIST)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_DOFLIST)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_DOFLIST)) then
        call storage_free(rgroupFEMSetDest%h_IdofList)
      end if
      
      ! Copy content from source to destination structure
      if (check(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DOFLIST)) then
        call storage_copy(rgroupFEMSetSrc%h_IdofList,&
            rgroupFEMSetDest%h_IdofList)
      end if
      
      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DOFLIST))
      
      ! Set ownership
      rgroupFEMSetDest%iduplicationFlag =&
          iand(rgroupFEMSetDest%iduplicationFlag, not(GFEM_SHARE_DOFLIST))
    end if


    ! Copy diagonal pointer
    if (check(idupFlag, GFEM_DUP_DIAGLIST)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_DIAGLIST)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_DIAGLIST)) then
        call storage_free(rgroupFEMSetDest%h_IdiagList)
      end if

      ! Copy content from source to destination structure
      if (check(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DIAGLIST)) then
        call storage_copy(rgroupFEMSetSrc%h_IdiagList,&
            rgroupFEMSetDest%h_IdiagList)
      end if

      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DIAGLIST))

      ! Set ownership
      rgroupFEMSetDest%iduplicationFlag =&
          iand(rgroupFEMSetDest%iduplicationFlag, not(GFEM_SHARE_DIAGLIST))
    end if


    ! Copy node structure
    if (check(idupFlag, GFEM_DUP_NODESTRUCTURE)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_NODESTRUCTURE)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_NODESTRUCTURE)) then
        call storage_free(rgroupFEMSetDest%h_InodeListIdx)
        call storage_free(rgroupFEMSetDest%h_InodeList)
      end if

      ! Copy content from source to destination structure
      if (check(rgroupFEMSetSrc%isetSpec, GFEM_HAS_NODESTRUCTURE)) then
        call storage_copy(rgroupFEMSetSrc%h_InodeListIdx,&
            rgroupFEMSetDest%h_InodeListIdx)
        call storage_copy(rgroupFEMSetSrc%h_InodeList,&
            rgroupFEMSetDest%h_InodeList)
      end if
        
      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_NODESTRUCTURE))
      
      ! Set ownership
      rgroupFEMSetDest%iduplicationFlag =&
          iand(rgroupFEMSetDest%iduplicationFlag, not(GFEM_SHARE_NODESTRUCTURE))
    end if


    ! Copy edge structure
    if (check(idupFlag, GFEM_DUP_EDGESTRUCTURE)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_EDGESTRUCTURE)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_EDGESTRUCTURE)) then
        call storage_free(rgroupFEMSetDest%h_IedgeListIdx)
        call storage_free(rgroupFEMSetDest%h_IedgeList)
      end if

      ! Copy content from source to destination structure
      if (check(rgroupFEMSetSrc%isetSpec, GFEM_HAS_EDGESTRUCTURE)) then
        call storage_copy(rgroupFEMSetSrc%h_IedgeListIdx,&
            rgroupFEMSetDest%h_IedgeListIdx)
        call storage_copy(rgroupFEMSetSrc%h_IedgeList,&
            rgroupFEMSetDest%h_IedgeList)
      end if

      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_EDGESTRUCTURE))
      
      ! Set ownership
      rgroupFEMSetDest%iduplicationFlag =&
          iand(rgroupFEMSetDest%iduplicationFlag, not(GFEM_SHARE_EDGESTRUCTURE))
    end if


    ! Copy coefficients at matrix diagonals
    if (check(idupFlag, GFEM_DUP_DIAGDATA)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_DIAGDATA)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_DIAGDATA)) then
        call storage_free(rgroupFEMSetDest%h_CoeffsAtDiag)
      end if

      ! Copy content from source to destination structure
      if (check(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DIAGDATA)) then
        call storage_copy(rgroupFEMSetSrc%h_CoeffsAtDiag,&
            rgroupFEMSetDest%h_CoeffsAtDiag)
      end if

      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DIAGDATA))
      
      ! Set ownership
      rgroupFEMSetDest%iduplicationFlag =&
          iand(rgroupFEMSetDest%iduplicationFlag, not(GFEM_SHARE_DIAGDATA))
    end if


    ! Copy coefficients at nodes
    if (check(idupFlag, GFEM_DUP_NODEDATA)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_NODEDATA)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_NODEDATA)) then
        call storage_free(rgroupFEMSetDest%h_CoeffsAtNode)
      end if

      ! Copy content from source to destination structure
      if (check(rgroupFEMSetSrc%isetSpec, GFEM_HAS_NODEDATA)) then
        call storage_copy(rgroupFEMSetSrc%h_CoeffsAtNode,&
            rgroupFEMSetDest%h_CoeffsAtNode)
      end if

      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_NODEDATA))
      
      ! Set ownership
      rgroupFEMSetDest%iduplicationFlag =&
          iand(rgroupFEMSetDest%iduplicationFlag, not(GFEM_SHARE_NODEDATA))
    end if


    ! Copy coefficients at edges
    if (check(idupFlag, GFEM_DUP_EDGEDATA)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_EDGEDATA)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_EDGEDATA)) then
        call storage_free(rgroupFEMSetDest%h_CoeffsAtEdge)
      end if

      ! Copy content from source to destination structure
      if (check(rgroupFEMSetSrc%isetSpec, GFEM_HAS_EDGEDATA)) then
        call storage_copy(rgroupFEMSetSrc%h_CoeffsAtEdge,&
            rgroupFEMSetDest%h_CoeffsAtEdge)
      end if

      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_EDGEDATA))
      
      ! Set ownership
      rgroupFEMSetDest%iduplicationFlag =&
          iand(rgroupFEMSetDest%iduplicationFlag, not(GFEM_SHARE_EDGEDATA))
    end if
      
  contains

    !**************************************************************
    ! Checks if iflag has all bits ibitfield set.

    pure function check(iflag, ibitfield)
      
      integer(I32), intent(in) :: iflag,ibitfield
      
      logical :: check
      
      check = (iand(iflag,ibitfield) .eq. ibitfield)

    end function check

    !**************************************************************
    ! Checks if ibitfield is not set in idupFlag.

    pure function checkOwner(idupFlag, ibitfield)
      
      integer(I32), intent(in) :: idupFlag,ibitfield
      
      logical :: checkOwner
      
      checkOwner = (iand(idupFlag,ibitfield) .ne. ibitfield)

    end function checkOwner

  end subroutine gfem_copyGroupFEMSet

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_copyGroupFEMBlock(rgroupFEMBlockSrc, rgroupFEMBlockDest, idupFlag)

!<description>
    ! This subroutine selectively copies data from the source
    ! structure rgroupFEMBlockSrc to the destination structure
    ! rgroupFEMBlockDest.
!<description>

!<input>
    ! Source block of group finite element sets
    type(t_groupFEMBlock), intent(in) :: rgroupFEMBlockSrc

    ! Duplication flag that decides on how to set up the structure
    integer(I32), intent(in) :: idupFlag
!</input>

!<inputoutput>
    ! Destination block of group finite element sets
    type(t_groupFEMBlock), intent(inout) :: rgroupFEMBlockDest
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! Check if block structure needs to be (re-)allocated
    if (rgroupFEMBlockSrc%nblocks .ne. rgroupFEMBlockDest%nblocks) then
      call gfem_releaseGroupFEMBlock(rgroupFEMBlockDest)
      rgroupFEMBlockDest%nblocks = rgroupFEMBlockSrc%nblocks
      allocate(rgroupFEMBlockDest%RgroupFEMBlock(rgroupFEMBlockDest%nblocks))
    end if
    
    do i = 1, rgroupFEMBlockSrc%nblocks
      call gfem_copyGroupFEMSet(rgroupFEMBlockSrc%RgroupFEMBlock(i),&
          rgroupFEMBlockDest%RgroupFEMBlock(i), idupFlag)
    end do
    
  end subroutine gfem_copyGroupFEMBlock

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_duplicateGroupFEMSet(rgroupFEMSetSrc, rgroupFEMSetDest,&
      idupFlag, bpreserveContent)

!<description>
    ! This subroutine duplicates parts of the source structure
    ! rgroupFEMSetSrc in the destination structure rgroupFEMSetDest.
    ! Note that rgroupFEMSetScr is still the owner of the duplicated
    ! content. If the optional flag bpreserveContent is set to TRUE,
    ! then the existing content is not removed before copying some
    ! content from the source structure.
!<description>

!<input>
    ! Source group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSetSrc

    ! Duplication flag that decides on how to set up the structure
    integer(I32), intent(in) :: idupFlag

    ! OPTIONAL: Flag to force that existing content is preserved
    logical, intent(in), optional :: bpreserveContent
!</input>

!<inputoutput>
    ! Destination group finite element set
    type(t_groupFEMSet), intent(inout) :: rgroupFEMSetDest
!</inputoutput>
!</subroutine>

    ! local variables
    logical :: bclear

    ! Clear destination structre unless preservation of existing
    ! content is explicitly enforced by the calling routine
    bclear = .true.
    if (present(bpreserveContent)) bclear=.not.bpreserveContent
    if (bclear) call gfem_releaseGroupFEMSet(rgroupFEMSetDest)

    ! Copy structural data
    if (check(idupFlag, GFEM_DUP_STRUCTURE)) then
      rgroupFEMSetDest%cassemblyType = rgroupFEMSetSrc%cassemblyType
      rgroupFEMSetDest%cdataType     = rgroupFEMSetSrc%cdataType
      rgroupFEMSetDest%NA            = rgroupFEMSetSrc%NA
      rgroupFEMSetDest%NEQ           = rgroupFEMSetSrc%NEQ
      rgroupFEMSetDest%NEDGE         = rgroupFEMSetSrc%NEDGE
      rgroupFEMSetDest%NVAR          = rgroupFEMSetSrc%NVAR
      rgroupFEMSetDest%ncoeffsAtDiag = rgroupFEMSetSrc%ncoeffsAtDiag
      rgroupFEMSetDest%ncoeffsAtNode = rgroupFEMSetSrc%ncoeffsAtNode
      rgroupFEMSetDest%ncoeffsAtEdge = rgroupFEMSetSrc%ncoeffsAtEdge
    end if


    ! Duplicate list of restricted DOFs
    if (check(idupFlag, GFEM_DUP_DOFLIST)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_DOFLIST)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_DIAGLIST)) then
        call storage_free(rgroupFEMSetDest%h_IdiagList)
      end if
      
      ! Copy handle from source to destination structure
      rgroupFEMSetDest%h_IdofList = rgroupFEMSetSrc%h_IdofList

      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DOFLIST))

      ! Reset ownership
      rgroupFEMSetDest%iduplicationFlag = ior(rgroupFEMSetDest%iduplicationFlag,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DOFLIST))
    end if


    ! Duplicate diagonal pointer
    if (check(idupFlag, GFEM_DUP_DIAGLIST)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_DOFLIST)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_DOFLIST)) then
        call storage_free(rgroupFEMSetDest%h_IdofList)
      end if

      ! Copy handle from source to destination structure
      rgroupFEMSetDest%h_IdiagList = rgroupFEMSetSrc%h_IdiagList
      
      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DIAGLIST))
      
      ! Reset ownership
      rgroupFEMSetDest%iduplicationFlag = ior(rgroupFEMSetDest%iduplicationFlag,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DIAGLIST))
    end if
    

    ! Duplicate node structure
    if (check(idupFlag, GFEM_DUP_NODESTRUCTURE)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_NODESTRUCTURE)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_NODESTRUCTURE)) then
        call storage_free(rgroupFEMSetDest%h_InodeListIdx)
        call storage_free(rgroupFEMSetDest%h_InodeList)
      end if

      ! Copy handle from source to destination structure
      rgroupFEMSetDest%h_InodeListIdx = rgroupFEMSetSrc%h_InodeListIdx
      rgroupFEMSetDest%h_InodeList = rgroupFEMSetSrc%h_InodeList
      
      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_NODESTRUCTURE))
      
      ! Reset ownership
      rgroupFEMSetDest%iduplicationFlag = ior(rgroupFEMSetDest%iduplicationFlag,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_NODESTRUCTURE))
    end if


    ! Duplicate edge structure
    if (check(idupFlag, GFEM_DUP_EDGESTRUCTURE)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_EDGESTRUCTURE)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_EDGESTRUCTURE)) then
        call storage_free(rgroupFEMSetDest%h_IedgeListIdx)
        call storage_free(rgroupFEMSetDest%h_IedgeList)
      end if

      ! Copy handle from source to destination structure
      rgroupFEMSetDest%h_IedgeListIdx = rgroupFEMSetSrc%h_IedgeListIdx
      rgroupFEMSetDest%h_IedgeList = rgroupFEMSetSrc%h_IedgeList
      
      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_EDGESTRUCTURE))
      
      ! Reset ownership
      rgroupFEMSetDest%iduplicationFlag = ior(rgroupFEMSetDest%iduplicationFlag,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_EDGESTRUCTURE))
    end if


    ! Duplicate coefficients at matrix diagonals
    if (check(idupFlag, GFEM_DUP_DIAGDATA)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_DIAGDATA)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_DIAGDATA)) then
        call storage_free(rgroupFEMSetDest%h_CoeffsAtDiag)
      end if

      ! Copy handle from source to destination structure
      rgroupFEMSetDest%h_CoeffsAtDiag = rgroupFEMSetSrc%h_CoeffsAtDiag

      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DIAGDATA))
      
      ! Reset ownership
      rgroupFEMSetDest%iduplicationFlag = ior(rgroupFEMSetDest%iduplicationFlag,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_DIAGDATA))
    end if


    ! Duplicate coefficients at nodes
    if (check(idupFlag, GFEM_DUP_NODEDATA)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_NODEDATA)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_NODEDATA)) then
        call storage_free(rgroupFEMSetDest%h_CoeffsAtNode)
      end if

      ! Copy handle from source to destination structure
      rgroupFEMSetDest%h_CoeffsAtNode = rgroupFEMSetSrc%h_CoeffsAtNode

      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_NODEDATA))
      
      ! Reset ownership
      rgroupFEMSetDest%iduplicationFlag = ior(rgroupFEMSetDest%iduplicationFlag,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_NODEDATA))
    end if


    ! Duplicate coefficients at edges
    if (check(idupFlag, GFEM_DUP_EDGEDATA)) then
      ! Remove existing data owned by the destination structure
      if (checkOwner(rgroupFEMSetDest%iduplicationFlag, GFEM_SHARE_EDGEDATA)&
          .and.check(rgroupFEMSetDest%isetSpec, GFEM_HAS_EDGEDATA)) then
        call storage_free(rgroupFEMSetDest%h_CoeffsAtEdge)
      end if

      ! Copy handle from source to destination structure
      rgroupFEMSetDest%h_CoeffsAtEdge = rgroupFEMSetSrc%h_CoeffsAtEdge

      ! Adjust specifier of the destination structure
      rgroupFEMSetDest%isetSpec = ior(rgroupFEMSetDest%isetSpec,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_EDGEDATA))
      
      ! Reset ownership
      rgroupFEMSetDest%iduplicationFlag = ior(rgroupFEMSetDest%iduplicationFlag,&
          iand(rgroupFEMSetSrc%isetSpec, GFEM_HAS_EDGEDATA))
    end if

  contains

    !**************************************************************
    ! Checks if iflag has all bits ibitfield set.

    pure function check(iflag, ibitfield)
      
      integer(I32), intent(in) :: iflag,ibitfield
      
      logical :: check
      
      check = (iand(iflag,ibitfield) .eq. ibitfield)

    end function check

    !**************************************************************
    ! Checks if ibitfield is not set in idupFlag.

    pure function checkOwner(idupFlag, ibitfield)
      
      integer(I32), intent(in) :: idupFlag,ibitfield
      
      logical :: checkOwner
      
      checkOwner = (iand(idupFlag,ibitfield) .ne. ibitfield)

    end function checkOwner    

  end subroutine gfem_duplicateGroupFEMSet

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_duplicateGroupFEMBlock(rgroupFEMBlockSrc, rgroupFEMBlockDest, idupFlag)

!<description>
    ! This subroutine duplicates parts of the source structure
    ! rgroupFEMBlockSrc in the destination structure rgroupFEMBlockDest.
    ! Note that rgroupFEMBlockScr is still the owner of the duplicated content.
!<description>

!<input>
    ! Source block of group finite element sets
    type(t_groupFEMBlock), intent(in) :: rgroupFEMBlockSrc

    ! Duplication flag that decides on how to set up the structure
    integer(I32), intent(in) :: idupFlag
!</input>

!<inputoutput>
    ! Destination block of group finite element sets
    type(t_groupFEMBlock), intent(inout) :: rgroupFEMBlockDest
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! Check if block structure needs to be (re-)allocated
    if (rgroupFEMBlockSrc%nblocks .ne. rgroupFEMBlockDest%nblocks) then
      call gfem_releaseGroupFEMBlock(rgroupFEMBlockDest)
      rgroupFEMBlockDest%nblocks = rgroupFEMBlockSrc%nblocks
      allocate(rgroupFEMBlockDest%RgroupFEMBlock(rgroupFEMBlockDest%nblocks))
    end if
    
    do i = 1, rgroupFEMBlockSrc%nblocks
      call gfem_duplicateGroupFEMSet(rgroupFEMBlockSrc%RgroupFEMBlock(i),&
          rgroupFEMBlockDest%RgroupFEMBlock(i), idupFlag)
    end do
        
  end subroutine gfem_duplicateGroupFEMBlock

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_initCoeffsFromMatrix(rgroupFEMSet, Rmatrices, Iposition,&
      cinitcoeffs)

!<description>
    ! This subroutine initialises the arrays of precomputed
    ! coefficient with the data given by the matrices. Depending on
    ! the type of assembly, the coefficients are stored in
    ! CoeffsAtNode, CoefssAtDiag and/or CoeffsAtEdge.
!</description>

!<input>
    ! Array of scalar coefficient matrices
    type(t_matrixScalar), dimension(:), intent(in) :: Rmatrices
    
    ! OPTIONAL: Array of integers which indicate the positions of the
    ! given matrices. If this parameter is not given, then the
    ! matrices are stored starting at position one.
    integer, dimension(:), intent(in), optional, target :: Iposition

    ! OPTIONAL: specifier that determines which coefficient arrays
    ! sould be initialised. If not given, then the decision is based
    ! on the assembly type specified in the group finite element set
    integer, intent(in), optional :: cinitcoeffs
!</input>

!<inputoutpu>
    ! group finite element set
    type(t_groupFEMSet), intent(inout) :: rgroupFEMSet
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtDiag,p_DcoeffsAtNode
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:,:,:), pointer :: p_FcoeffsAtEdge
    real(SP), dimension(:,:), pointer :: p_FcoeffsAtDiag,p_FcoeffsAtNode
    real(SP), dimension(:), pointer :: p_Fdata
    integer, dimension(:,:), pointer :: p_IedgeList,p_InodeList,p_IdiagList
    integer, dimension(:), pointer :: p_Iposition,p_IdofList,p_InodeListIdx
    integer, dimension(2) :: Isize2D
    integer, dimension(3) :: Isize3D
    integer :: cinit,ia,iedge,ieq,ii,ij,imatrix,ipos,ji,nmatrices,nmaxpos

    ! Set number of matrices
    nmatrices = size(Rmatrices)

    ! Determine types of coefficients to be initialised
    cinit = 2_I32**rgroupFEMSet%cassemblyType
    if (present(cinitcoeffs)) cinit=cinitcoeffs
    
    ! Set pointer to matrix positions
    if (present(Iposition)) then
      p_Iposition => Iposition
    else
      allocate(p_Iposition(size(Rmatrices)))
      do imatrix = 1, nmatrices
        p_Iposition(imatrix) = imatrix
      end do
    end if
    
    ! Set maximum position
    nmaxpos = maxval(p_Iposition)

    ! Check if array Rmatrices and Iposition have the same size
    if (nmatrices .ne. size(p_Iposition)) then
      call output_line('size(Rmatrices) /= size(Iposition)!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
      call sys_halt()
    end if
    
    !---------------------------------------------------------------------------
    ! Node-based assembly
    !---------------------------------------------------------------------------

    if (iand(cinit, GFEM_INITCOEFFS_NODEBASED) .eq. GFEM_INITCOEFFS_NODEBASED) then

      ! Check if structure provides node-based structure and coefficient array
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
        call sys_halt()
      end if
      
      ! Check if coefficient array has enough positions
      call storage_getsize(rgroupFEMSet%h_CoeffsAtNode, Isize2D)
      if (Isize2D(1) .lt. nmaxpos) then
        call output_line('NMAXPOS exceeds dimension of CoeffsAtNode!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
        call sys_halt()
      end if

      ! Set pointers
      call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .ne. 0) then
        call gfem_getbase_IdofList(rgroupFEMSet, p_IdofList)
        call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
      end if
      
      ! What data type we we?
      select case(rgroupFEMSet%cdataType)
      case (ST_DOUBLE)
        call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)

        ! Associate data of each matrix separately
        do imatrix = 1, nmatrices
          
          ! Get matrix position
          ipos = p_Iposition(imatrix)
          
          ! What data type are we?
          select case(Rmatrices(imatrix)%cdataType)
          case (ST_DOUBLE)
            !-------------------------------------------------------------------
            ! Copy double-precision matrix to double-precision coefficients
            !-------------------------------------------------------------------
            
            ! Set pointer to matrix data
            call lsyssc_getbase_double(Rmatrices(imatrix), p_Ddata)
            
            if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
              ! Loop over all nodes and copy matrix entries
              do ia = 1, rgroupFEMSet%NA
                p_DcoeffsAtNode(ipos,ia) = p_Ddata(ia)
              end do
            else
              ! Loop over selected nodes and copy matrix entries
              do ia = 1, rgroupFEMSet%NA
                p_DcoeffsAtNode(ipos,ia) = p_Ddata(p_InodeList(1,ia))
              end do
            end if
            
          case (ST_SINGLE)
            !-------------------------------------------------------------------
            ! Copy single-precision matrix to double-precision coefficients
            !-------------------------------------------------------------------

            ! Set pointer to matrix data
            call lsyssc_getbase_single(Rmatrices(imatrix), p_Fdata)
            
            if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
              ! Loop over all nodes and copy matrix entries
              do ia = 1, rgroupFEMSet%NA
                p_DcoeffsAtNode(ipos,ia) = p_Fdata(ia)
              end do
            else
              ! Loop over selected nodes and copy matrix entries
              do ia = 1, rgroupFEMSet%NA
                p_DcoeffsAtNode(ipos,ia) = p_Fdata(p_InodeList(1,ia))
              end do
            end if
            
          case default
            call output_line('Unsupported data type!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
            call sys_halt()
          end select
        end do
        
        ! Set specifier for precomputed nodal coefficients
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA)
        
      case (ST_SINGLE)
        call gfem_getbase_FcoeffsAtNode(rgroupFEMSet, p_FcoeffsAtNode)

        ! Associate data of each matrix separately
        do imatrix = 1, nmatrices
          
          ! Get matrix position
          ipos = p_Iposition(imatrix)
          
          ! What data type are we?
          select case(Rmatrices(imatrix)%cdataType)
          case (ST_DOUBLE)
            !-------------------------------------------------------------------
            ! Copy double-precision matrix to single-precision coefficients
            !-------------------------------------------------------------------

            ! Set pointer to matrix data
            call lsyssc_getbase_double(Rmatrices(imatrix), p_Ddata)
            
            if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
              ! Loop over all nodes and copy matrix entries
              do ia = 1, rgroupFEMSet%NA
                p_FcoeffsAtNode(ipos,ia) = p_Ddata(ia)
              end do
            else
              ! Loop over selected nodes and copy matrix entries
              do ia = 1, rgroupFEMSet%NA
                p_FcoeffsAtNode(ipos,ia) = p_Ddata(p_InodeList(1,ia))
              end do
            end if
            
          case (ST_SINGLE)
            !-------------------------------------------------------------------
            ! Copy single-precision matrix to single-precision coefficients
            !-------------------------------------------------------------------
            
            ! Set pointer to matrix data
            call lsyssc_getbase_single(Rmatrices(imatrix), p_Fdata)
            
            if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
              ! Loop over all nodes and copy matrix entries
              do ia = 1, rgroupFEMSet%NA
                p_FcoeffsAtNode(ipos,ia) = p_Fdata(ia)
              end do
            else
              ! Loop over selected nodes and copy matrix entries
              do ia = 1, rgroupFEMSet%NA
                p_DcoeffsAtNode(ipos,ia) = p_Fdata(p_InodeList(1,ia))
              end do
            end if
            
          case default
            call output_line('Unsupported data type!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
            call sys_halt()
          end select
        end do
        
        ! Set specifier for precomputed nodal coefficients
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA)
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
        call sys_halt()
      end select
    end if
               
    
    !---------------------------------------------------------------------------
    ! Edge-based assembly
    !---------------------------------------------------------------------------

    if (iand(cinit, GFEM_INITCOEFFS_EDGEBASED) .eq. GFEM_INITCOEFFS_EDGEBASED) then
      
      ! Check if structure provides edge-based structure and coefficient array
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGESTRUCTURE) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGLIST)      .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
        call sys_halt()
      end if

      ! Check if coefficient arrays has enough positions
      call storage_getsize(rgroupFEMSet%h_CoeffsAtDiag, Isize2D)
      call storage_getsize(rgroupFEMSet%h_CoeffsAtEdge, Isize3D)
      if ((Isize2D(1) .lt. nmaxpos) .or. (Isize3D(1) .lt. nmaxpos)) then
        call output_line('NMAXPOS exceeds dimension of CoeffsAtDiag/CoefssAtEdge!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
        call sys_halt()
      end if

      ! Set pointers
      call gfem_getbase_IdiagList(rgroupFEMSet, p_IdiagList)
      call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)
     
      ! What data type we we?
      select case(rgroupFEMSet%cdataType)
      case (ST_DOUBLE)
        call gfem_getbase_DcoeffsAtDiag(rgroupFEMSet, p_DcoeffsAtDiag)
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)

        ! Associate data of each matrix separately
        do imatrix = 1, nmatrices
          
          ! Get matrix position
          ipos = p_Iposition(imatrix)
          
          ! What data type are we?
          select case(Rmatrices(imatrix)%cdataType)
          case (ST_DOUBLE)
            !-------------------------------------------------------------------
            ! Copy double-precision matrix to double-precision coefficients
            !-------------------------------------------------------------------

            ! Set pointer to matrix data
            call lsyssc_getbase_double(Rmatrices(imatrix), p_Ddata)
            
            ! Loop over all edges and copy off-diagonal entries
            do iedge = 1, rgroupFEMSet%NEDGE
              ij = p_IedgeList(3,iedge)
              ji = p_IedgeList(4,iedge)
              p_DcoeffsAtEdge(ipos,1,iedge) = p_Ddata(ij)
              p_DcoeffsAtEdge(ipos,2,iedge) = p_Ddata(ji)
            end do

            ! Loop over all equations and copy diagonal entroes
            do ieq = 1, rgroupFEMSet%NEQ
              ii = p_IdiagList(2,ieq)
              p_DcoeffsAtDiag(ipos,ieq) = p_Ddata(ii)
            end do
            
          case (ST_SINGLE)
            !-------------------------------------------------------------------
            ! Copy single-precision matrix to double-precision coefficients
            !-------------------------------------------------------------------

            ! Set pointer to matrix data
            call lsyssc_getbase_single(Rmatrices(imatrix), p_Fdata)
            
            ! Loop over all edges and copy off-diagonal entries
            do iedge = 1, rgroupFEMSet%NEDGE
              ij = p_IedgeList(3,iedge)
              ji = p_IedgeList(4,iedge)
              p_DcoeffsAtEdge(ipos,1,iedge) = p_Fdata(ij)
              p_DcoeffsAtEdge(ipos,2,iedge) = p_Fdata(ji)
            end do
            
            ! Loop over all equations and copy diagonal entroes
            do ieq = 1, rgroupFEMSet%NEQ
              ii = p_IdiagList(2,ieq)
              p_DcoeffsAtDiag(ipos,ieq) = p_Fdata(ii)
            end do

          case default
            call output_line('Unsupported data type!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
            call sys_halt()
          end select
        end do
        
        ! Set specifier for precomputed edge and node coefficients
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGDATA)
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA)
        

      case (ST_SINGLE)
        call gfem_getbase_FcoeffsAtDiag(rgroupFEMSet, p_FcoeffsAtDiag)
        call gfem_getbase_FcoeffsAtEdge(rgroupFEMSet, p_FcoeffsAtEdge)

        ! Associate data of each matrix separately
        do imatrix = 1, nmatrices
          
          ! Get matrix position
          ipos = p_Iposition(imatrix)
          
          ! What data type are we?
          select case(Rmatrices(imatrix)%cdataType)
          case (ST_DOUBLE)
            !-------------------------------------------------------------------
            ! Copy double-precision matrix to single-precision coefficients
            !-------------------------------------------------------------------

            ! Set pointer to matrix data
            call lsyssc_getbase_double(Rmatrices(imatrix), p_Ddata)
            
            ! Loop over all edges and copy off-diagonal entries
            do iedge = 1, rgroupFEMSet%NEDGE
              ij = p_IedgeList(3,iedge)
              ji = p_IedgeList(4,iedge)
              p_FcoeffsAtEdge(ipos,1,iedge) = p_Ddata(ij)
              p_FcoeffsAtEdge(ipos,2,iedge) = p_Ddata(ji)
            end do

            ! Loop over all equations and copy diagonal entroes
            do ieq = 1, rgroupFEMSet%NEQ
              ii = p_IdiagList(2,ieq)
              p_FcoeffsAtDiag(ipos,ieq) = p_Ddata(ii)
            end do
            
          case (ST_SINGLE)
            !-------------------------------------------------------------------
            ! Copy single-precision matrix to single-precision coefficients
            !-------------------------------------------------------------------

            ! Set pointer to matrix data
            call lsyssc_getbase_single(Rmatrices(imatrix), p_Fdata)
            
            ! Loop over all edges and copy off-diagonal entries
            do iedge = 1, rgroupFEMSet%NEDGE
              ij = p_IedgeList(3,iedge)
              ji = p_IedgeList(4,iedge)
              p_FcoeffsAtEdge(ipos,1,iedge) = p_Fdata(ij)
              p_FcoeffsAtEdge(ipos,2,iedge) = p_Fdata(ji)
            end do
            
            ! Loop over all equations and copy diagonal entroes
            do ieq = 1, rgroupFEMSet%NEQ
              ii = p_IdiagList(2,ieq)
              p_FcoeffsAtDiag(ipos,ieq) = p_Fdata(ii)
            end do

          case default
            call output_line('Unsupported data type!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
            call sys_halt()
          end select
        end do
        
        ! Set specifier for precomputed edge and node coefficients
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGDATA)
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA)
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_initCoeffsFromMatrix')
        call sys_halt()
      end select
    end if

  end subroutine gfem_initCoeffsFromMatrix

  !*****************************************************************************

!<subroutine>

  subroutine gfem_isMatrixCompatibleSc(rgroupFEMSet, rmatrix, bcompatible,&
      bexpensiveCheck)

!<description>
    ! This subroutine checks if a scalar matrix and a group finite
    ! element set compatible to each other, i.e. if they share the
    ! same structure, size and so on.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! OPTIONAL: If given and set to true, then potentially more
    ! expensive consistency checks are performed
    logical, intent(in), optional :: bexpensiveCheck
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE
    ! depending on whether the matrix and the group finite element set
    ! are compatible or not.  If not given, an error will inform the
    ! user if the matrix/operator are not compatible and the program
    ! will halt.
    logical, intent(out), optional :: bcompatible
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdofList
    integer :: na,ncols,nedge,neq
    logical :: bexpensive

    ! Perform expensive checks?
    bexpensive = .false.
    if (present(bexpensiveCheck)) bexpensive = bexpensiveCheck

    ! Perform expensive checks?
    if (bexpensive) then
      
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
    
        ! Matrix/group finite element set must have the same size
        if ((rgroupFEMSet%NVAR  .ne. rmatrix%NVAR) .or.&
            (rgroupFEMSet%NEDGE .ne. (rmatrix%NA-rmatrix%NEQ)/2)) then
          if (present(bcompatible)) then
            bcompatible = .false.
            return
          else
            call output_line('Matrix/group finite element set not compatible,'//&
                ' different structure!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_isMatrixCompatibleSc')
            call sys_halt()
          end if
        end if
        
      else   ! consider list of restricted DOFs
        
        call gfem_getbase_IdofList(rgroupFEMSet, p_IdofList)
        call lsyssc_calcDimsFromMatrix(rmatrix, na, neq, ncols, nedge, p_IdofList)
        
        ! Matrix/group finite element set must have the same size
        if ((rgroupFEMSet%NA    .ne. na)           .or.&
            (rgroupFEMSet%NEQ   .ne. neq)          .or.&
            (rgroupFEMSet%NVAR  .ne. rmatrix%NVAR) .or.&
            (rgroupFEMSet%NEDGE .ne. nedge)) then
          if (present(bcompatible)) then
            bcompatible = .false.
            return
          else
            call output_line('Matrix/group finite element set not compatible,'//&
                ' different structure!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_isMatrixCompatibleSc')
            call sys_halt()
          end if
        end if
        
      end if
    end if
    
    ! Perform standard checks
    if ((rgroupFEMSet%Iconsistency(1) .ne. rmatrix%NA)  .or.&
        (rgroupFEMSet%Iconsistency(2) .ne. rmatrix%NEQ) .or.&
        (rgroupFEMSet%Iconsistency(3) .ne. rmatrix%NCOLS)) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Matrix/group finite element set not compatible,'//&
            ' different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_isMatrixCompatibleSc')
        call sys_halt()
      end if
    end if
    
  end subroutine gfem_isMatrixCompatibleSc

  ! *****************************************************************************
  
!<subroutine>

  subroutine gfem_isMatrixCompatibleBl(rgroupFEMSet, rmatrix, bcompatible,&
      bexpensiveCheck)

!<description>
    ! This subroutine checks whether a block matrix and a group finite
    ! element set are compatible to each other, i.e. if they share the
    ! same structure, size and so on.
    !
    ! If the matrix has only one block, then the scalar counterpart of this
    ! subroutine is called with the corresponding scalar submatrix.
    ! Otherwise, the matrix is required to possess group structure.
!</description>

!<input>
    ! Block matrix
    type(t_matrixBlock), intent(in) :: rmatrix
    
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! OPTIONAL: If given and set to true, then potentially more
    ! expensive consistency checks are performed
    logical, intent(in), optional :: bexpensiveCheck
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE
    ! depending on whether the matrix and the group finite element set
    ! are compatible or not.  If not given, an error will inform the
    ! user if the matrix/operator are not compatible and the program
    ! will halt.
    logical, intent(out), optional :: bcompatible
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdofList
    integer :: na,ncols,nedge,neq
    logical :: bexpensive
    
    ! Check if matrix has only one block
    if ((rmatrix%nblocksPerCol .eq. 1) .and.&
        (rmatrix%nblocksPerRow .eq. 1)) then
      call gfem_isMatrixCompatible(rgroupFEMSet,&
          rmatrix%RmatrixBlock(1,1), bcompatible)
      return
    end if

    ! Check if number of columns equans number of rows
    if (rmatrix%nblocksPerCol .ne.&
        rmatrix%nblocksPerRow) then
      call output_line('Block matrix must have equal number of columns and rows!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfem_isMatrixCompatibleBl')
      call sys_halt()
    end if

    ! Check if matrix exhibits group structure
    if (rmatrix%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Block matrix must have group structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_isMatrixCompatibleBl')
        call sys_halt()
      end if
    end if

    ! Perform expensive checks?
    bexpensive = .false.
    if (present(bexpensiveCheck)) bexpensive = bexpensiveCheck
    
    ! Perform expensive checks?
    if (bexpensive) then
      
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
        
        ! Matrix/operator must have the same size
        if ((rgroupFEMSet%NVAR  .ne. rmatrix%nblocksPerCol)         .or.&
            (rgroupFEMSet%NA    .ne. rmatrix%RmatrixBlock(1,1)%NA)  .or.&
            (rgroupFEMSet%NEQ   .ne. rmatrix%RmatrixBlock(1,1)%NEQ) .or.&
            (rgroupFEMSet%NEDGE .ne. (rmatrix%RmatrixBlock(1,1)%NA-&
            rmatrix%RmatrixBlock(1,1)%NEQ)/2)) then
          if (present(bcompatible)) then
            bcompatible = .false.
            return
          else
            call output_line('Matrix/group finite element set not compatible,'//&
                ' different structure!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_isMatrixCompatibleBl')
            call sys_halt()
          end if
        end if
        
      else   ! consider list of restricted DOFs
        
        call gfem_getbase_IdofList(rgroupFEMSet, p_IdofList)
        call lsyssc_calcDimsFromMatrix(rmatrix%RmatrixBlock(1,1),&
            na, neq, ncols, nedge, p_IdofList)
        
        ! Matrix/group finite element set must have the same size
        if ((rgroupFEMSet%NVAR  .ne. rmatrix%nblocksPerCol) .or.&
            (rgroupFEMSet%NA    .ne. na)                    .or.&
            (rgroupFEMSet%NEQ   .ne. neq)                   .or.&
            (rgroupFEMSet%NEDGE .ne. nedge)) then
          if (present(bcompatible)) then
            bcompatible = .false.
            return
          else
            call output_line('Matrix/group finite element set not compatible,'//&
                ' different structure!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfem_isMatrixCompatibleBl')
            call sys_halt()
          end if
        end if

      end if
    end if

    ! Perform standard checks
    if ((rgroupFEMSet%Iconsistency(1) .ne. rmatrix%RmatrixBlock(1,1)%NA)  .or.&
        (rgroupFEMSet%Iconsistency(2) .ne. rmatrix%RmatrixBlock(1,1)%NEQ) .or.&
        (rgroupFEMSet%Iconsistency(3) .ne. rmatrix%RmatrixBlock(1,1)%NCOLS)) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Matrix/group finite element set not compatible,'//&
            ' different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_isMatrixCompatibleBl')
        call sys_halt()
      end if
    end if      

  end subroutine gfem_isMatrixCompatibleBl

  !*****************************************************************************

!<subroutine>

  subroutine gfem_isVectorCompatibleSc(rgroupFEMSet, rvector, bcompatible)

!<description>
    ! This subroutine checks if a vector and a group finite element
    ! set are compatible to each other, i.e., share the same
    ! structure, size and so on.
!</description>

!<input>
    ! Scalar vector
    type(t_vectorScalar), intent(in) :: rvector

    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE
    ! depending on whether the matrix and the group finite element set
    ! are compatible or not.  If not given, an error will inform the
    ! user if the matrix/operator are not compatible and the program
    ! will halt.
    logical, intent(out), optional :: bcompatible
!</output>
!</subroutine>

    ! Vector/group finite element set must have the same size
    if ((rgroupFEMSet%Iconsistency(2) .ne. rvector%NEQ) .or.&
        (rgroupFEMSet%NVAR            .ne. rvector%NVAR)) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Vector/group finite element set not compatible,'//&
            ' different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_isVectorCompatibleSc')
        call sys_halt()
      end if
    end if
    
  end subroutine gfem_isVectorCompatibleSc

!*****************************************************************************

!<subroutine>

  subroutine gfem_isVectorCompatibleBl(rgroupFEMSet, rvectorBlock, bcompatible)

!<description>
    ! This subroutine checks whether a block vector and a group finite
    ! element set are compatible to each other, i.e., share the same
    ! structure, size and so on.
    !
    ! If the vectors has only one block, then the scalar counterpart of
    ! this subroutine is called with the corresponding scalar subvector.
!</description>

!<input>
    ! Block vector
    type(t_vectorBlock), intent(in) :: rvectorBlock

    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and group finite element set are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    logical, intent(out), optional :: bcompatible
!</output>
!</subroutine>

    ! Check if block vectors has just one block
    if (rvectorBlock%nblocks .eq. 1) then
      call gfem_isVectorCompatible(rgroupFEMSet,&
          rvectorBlock%RvectorBlock(1), bcompatible)
      return
    end if
    
    ! Vector/operator must have the same size
    if (rgroupFEMSet%Iconsistency(2)*rgroupFEMSet%NVAR&
        .ne. rvectorBlock%NEQ) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Vector/group finite element set not compatible,'//&
            ' different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_isVectorCompatibleBl')
        call sys_halt()
      end if
    end if
    
  end subroutine gfem_isVectorCompatibleBl

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)

!<description>
    ! Returns a pointer to the index pointer for the edge structure
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the edge index structure
    ! NULL() if the structure rgroupFEMSet does not provide it.
    integer, dimension(:), pointer :: p_IedgeListIdx
!</output>
!</subroutine>

    ! Do we have an edge structure at all?
    if (rgroupFEMSet%h_IedgeListIdx .eq. ST_NOHANDLE) then
      nullify(p_IedgeListIdx)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rgroupFEMSet%h_IedgeListIdx,&
        p_IedgeListIdx)

  end subroutine gfem_getbase_IedgeListIdx
  
  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)

!<description>
    ! Returns a pointer to the edge structure
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the edge structure
    ! NULL() if the structure rgroupFEMSet does not provide it.
    integer, dimension(:,:), pointer :: p_IedgeList
!</output>
!</subroutine>

    ! Do we have an edge structure at all?
    if ((rgroupFEMSet%h_IedgeList .eq. ST_NOHANDLE) .or.&
        (rgroupFEMSet%NEDGE       .eq. 0)) then
      nullify(p_IedgeList)
      return
    end if
    
    ! Get the array
    call storage_getbase_int2D(rgroupFEMSet%h_IedgeList,&
        p_IedgeList, rgroupFEMSet%NEDGE)

  end subroutine gfem_getbase_IedgeList

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)

!<description>
    ! Returns a pointer to the index pointer for the nodal structure
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the nodal index structure
    ! NULL() if the structure rgroupFEMSet does not provide it.
    integer, dimension(:), pointer :: p_InodeListIdx
!</output>
!</subroutine>

    ! Do we have a node structure at all?
    if (rgroupFEMSet%h_InodeListIdx .eq. ST_NOHANDLE) then
      nullify(p_InodeListIdx)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rgroupFEMSet%h_InodeListIdx,&
        p_InodeListIdx, rgroupFEMSet%NEQ+1)

  end subroutine gfem_getbase_InodeListIdx

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_InodeList1D(rgroupFEMSet, p_InodeList)

!<description>
    ! Returns a pointer to the node structure
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the node structure
    ! NULL() if the structure rgroupFEMSet does not provide it.
    integer, dimension(:), pointer :: p_InodeList
!</output>
!</subroutine>

    ! Do we have a node structure at all?
    if ((rgroupFEMSet%h_InodeList .eq. ST_NOHANDLE) .or.&
        (rgroupFEMSet%NA          .eq. 0)) then
      nullify(p_InodeList)
      return
    end if
       
    ! Get the array
    call storage_getbase_int(rgroupFEMSet%h_InodeList,&
        p_InodeList, rgroupFEMSet%NA)
    
  end subroutine gfem_getbase_InodeList1D

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_InodeList2D(rgroupFEMSet, p_InodeList)

!<description>
    ! Returns a pointer to the node structure
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the node structure
    ! NULL() if the structure rgroupFEMSet does not provide it.
    integer, dimension(:,:), pointer :: p_InodeList
!</output>
!</subroutine>

    ! Do we have a node structure at all?
    if ((rgroupFEMSet%h_InodeList .eq. ST_NOHANDLE) .or.&
        (rgroupFEMSet%NA          .eq. 0)) then
      nullify(p_InodeList)
      return
    end if
    
    ! Get the array
    call storage_getbase_int2D(rgroupFEMSet%h_InodeList,&
        p_InodeList, rgroupFEMSet%NA)

  end subroutine gfem_getbase_InodeList2D

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_IdiagList(rgroupFEMSet, p_IdiagList)

!<description>
    ! Returns a pointer to the equation list
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the node structure
    ! NULL() if the structure rgroupFEMSet does not provide it.
    integer, dimension(:,:), pointer :: p_IdiagList
!</output>
!</subroutine>

    ! Do we have a node structure at all?
    if (rgroupFEMSet%h_IdiagList .eq. ST_NOHANDLE) then
      nullify(p_IdiagList)
      return
    end if
       
    ! Get the array
    call storage_getbase_int2D(rgroupFEMSet%h_IdiagList,&
        p_IdiagList, rgroupFEMSet%NEQ)
    
  end subroutine gfem_getbase_IdiagList

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)

!<description>
    ! Returns a pointer to the double-valued coefficients at nodes
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the precomputed coefficients at nodes
    ! NULL() if the structure rgroupFEMSet does not provide it.
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
!</output>
!</subroutine>

    ! Do we have coefficients at nodes at all?
    if (rgroupFEMSet%h_CoeffsAtNode .eq. ST_NOHANDLE) then
      nullify(p_DcoeffsAtNode)
      return
    end if
    
    ! Get the array
    select case(rgroupFEMSet%cassemblyType)
    case (GFEM_NODEBASED)
      if (rgroupFEMSet%NA .eq. 0) then
        nullify(p_DcoeffsAtNode)
      else
        call storage_getbase_double2D(rgroupFEMSet%h_CoeffsAtNode,&
            p_DcoeffsAtNode, rgroupFEMSet%NA)
      end if

    case (GFEM_EDGEBASED)
      if (rgroupFEMSet%NEQ .eq. 0) then
        nullify(p_DcoeffsAtNode)
      else
        call storage_getbase_double2D(rgroupFEMSet%h_CoeffsAtNode,&
            p_DcoeffsAtNode, rgroupFEMSet%NEQ)
      end if

    case default
      nullify(p_DcoeffsAtNode)
    end select

  end subroutine gfem_getbase_DcoeffsAtNode

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_FcoeffsAtNode(rgroupFEMSet, p_FcoeffsAtNode)

!<description>
    ! Returns a pointer to the single-valued coefficients at nodes
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the precomputed coefficients at nodes
    ! NULL() if the structure rgroupFEMSet does not provide it.
    real(SP), dimension(:,:), pointer :: p_FcoeffsAtNode
!</output>
!</subroutine>

    ! Do we have coefficients at nodes at all?
    if (rgroupFEMSet%h_CoeffsAtNode .eq. ST_NOHANDLE) then
      nullify(p_FcoeffsAtNode)
      return
    end if
    
    ! Get the array
    select case(rgroupFEMSet%cassemblyType)
    case (GFEM_NODEBASED)
      if (rgroupFEMSet%NA .eq. 0) then
        nullify(p_FcoeffsAtNode)
      else
        call storage_getbase_single2D(rgroupFEMSet%h_CoeffsAtNode,&
            p_FcoeffsAtNode, rgroupFEMSet%NA)
      end if

    case (GFEM_EDGEBASED)
      if (rgroupFEMSet%NEQ .eq. 0) then
        nullify(p_FcoeffsAtNode)
      else
        call storage_getbase_single2D(rgroupFEMSet%h_CoeffsAtNode,&
            p_FcoeffsAtNode, rgroupFEMSet%NEQ)
      end if

    case default
      nullify(p_FcoeffsAtNode)
    end select

  end subroutine gfem_getbase_FcoeffsAtNode

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_DcoeffsAtEdge(rgroupFEMSet,p_DcoeffsAtEdge)

!<description>
    ! Returns a pointer to the double-valued coefficients at edges
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the precomputed coefficients at edges
    ! NULL() if the structure rgroupFEMSet does not provide it.
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
!</output>
!</subroutine>

    ! Do we have edge data at all?
    if ((rgroupFEMSet%h_CoeffsAtEdge .eq. ST_NOHANDLE) .or.&
        (rgroupFEMSet%NEDGE          .eq. 0)) then
      nullify(p_DcoeffsAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_double3D(rgroupFEMSet%h_CoeffsAtEdge,&
        p_DcoeffsAtEdge, rgroupFEMSet%NEDGE)

  end subroutine gfem_getbase_DcoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_FcoeffsAtEdge(rgroupFEMSet, p_FcoeffsAtEdge)

!<description>
    ! Returns a pointer to the single-valued coefficients at edges
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the precomputed coefficients at edges
    ! NULL() if the structure rgroupFEMSet does not provide it.
    real(SP), dimension(:,:,:), pointer :: p_FcoeffsAtEdge
!</output>
!</subroutine>

    ! Do we have edge data at all?
    if ((rgroupFEMSet%h_CoeffsAtEdge .eq. ST_NOHANDLE) .or.&
        (rgroupFEMSet%NEDGE          .eq. 0)) then
      nullify(p_FcoeffsAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_single3D(rgroupFEMSet%h_CoeffsAtEdge,&
        p_FcoeffsAtEdge, rgroupFEMSet%NEDGE)

  end subroutine gfem_getbase_FcoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_DcoeffsAtDiag(rgroupFEMSet,p_DcoeffsAtDiag)

!<description>
    ! Returns a pointer to the double-valued coefficients at matrix diagonals
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the precomputed coefficients at matrix diagonales
    ! NULL() if the structure rgroupFEMSet does not provide it.
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtDiag
!</output>
!</subroutine>

    ! Do we have diagonal data at all?
    if ((rgroupFEMSet%h_CoeffsAtDiag .eq. ST_NOHANDLE) .or.&
        (rgroupFEMSet%NEQ            .eq. 0)) then
      nullify(p_DcoeffsAtDiag)
      return
    end if
    
    ! Get the array
    call storage_getbase_double2D(rgroupFEMSet%h_CoeffsAtDiag,&
        p_DcoeffsAtDiag, rgroupFEMSet%NEQ)

  end subroutine gfem_getbase_DcoeffsAtDiag

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_FcoeffsAtDiag(rgroupFEMSet, p_FcoeffsAtDiag)

!<description>
    ! Returns a pointer to the single-valued coefficients at matrix diagonals
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the precomputed coefficients at matrix diagonals
    ! NULL() if the structure rgroupFEMSet does not provide it.
    real(SP), dimension(:,:), pointer :: p_FcoeffsAtDiag
!</output>
!</subroutine>

    ! Do we have edge data at all?
    if ((rgroupFEMSet%h_CoeffsAtNode .eq. ST_NOHANDLE) .or.&
        (rgroupFEMSet%NEQ            .eq. 0)) then
      nullify(p_FcoeffsAtDiag)
      return
    end if
    
    ! Get the array
    call storage_getbase_single2D(rgroupFEMSet%h_CoeffsAtDiag,&
        p_FcoeffsAtDiag, rgroupFEMSet%NEQ)

  end subroutine gfem_getbase_FcoeffsAtDiag

  !*****************************************************************************

!<subroutine>

  subroutine gfem_getbase_IdofList(rgroupFEMSet, p_IdofList)

!<description>
    ! Returns a pointer to the list of restricted DOFs
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<output>
    ! Pointer to the list of restricted DOFs
    ! NULL() if the structure rgroupFEMSet does not provide it.
    integer, dimension(:), pointer :: p_IdofList
!</output>
!</subroutine>

    ! Do we have edge data at all?
    if (rgroupFEMSet%h_IdofList .eq. ST_NOHANDLE) then
      nullify(p_IdofList)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rgroupFEMSet%h_IdofList, p_IdofList)

  end subroutine gfem_getbase_IdofList

  !*****************************************************************************

!<subroutine>

  subroutine gfem_genNodeList(rmatrix, rgroupFEMSet)

!<description>IdofList
    ! This subroutine stores the list of degress of freedom to the
    ! group finite element set. If no list of degrees of freedom is
    ! provided, then the corresponding data structure is deallocated.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<inputoutput>
    ! Group finite element structure
    type(t_groupFEMSet), intent(inout) :: rgroupFEMSet
!</inputoutput>
!</subroutine>

    ! local variables
    logical, dimension(:), allocatable :: BisActive
    integer, dimension(:,:), pointer :: p_InodeList
    integer, dimension(:), pointer :: p_InodeListIdx,p_InodeList1D
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Kdiagonal
    integer, dimension(:), pointer :: p_IdofList
    integer, dimension(2) :: Isize2D
    integer :: idx,ieq,iidx,ij,isize,jcol
    
    ! Check if node structure is owned by the structure
    if (iand(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_NODESTRUCTURE) .eq.&
        GFEM_SHARE_NODESTRUCTURE) then
      call output_line('Node list is not owned by structure and '//&
          'therefore cannot be generated',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfem_genNodeList')
      call sys_halt()
    end if
    
    ! Check if list of DOFs has been attached
    if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then

      ! Check if matrix and group finite element set are compatible
      if ((rgroupFEMSet%NA    .ne. rmatrix%NA)  .or.&
          (rgroupFEMSet%NEQ   .ne. rmatrix%NEQ) .or.&
          (rgroupFEMSet%NEDGE .ne. (rmatrix%NA-rmatrix%NEQ)/2)) then
        call output_line('Matrix/group finite element set not compatible,'//&
            ' different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_genNodeList')
        call sys_halt()
      end if
    
      ! General node data structure
      select case(rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX1)
        
        ! Allocate new memory for node list if required
        if (rgroupFEMSet%h_InodeList .eq. ST_NOHANDLE) then
          call storage_new('gfem_genNodeList', 'InodeList',&
              rgroupFEMSet%NA, ST_INT, rgroupFEMSet%h_InodeList,&
              ST_NEWBLOCK_NOINIT)
        else
          call storage_getsize(rgroupFEMSet%h_InodeList, isize)
          if (isize .ne. rgroupFEMSet%NA) then
            call storage_realloc('gfem_genNodeList',&
                rgroupFEMSet%NA, ST_INT, rgroupFEMSet%h_InodeList,&
                ST_NEWBLOCK_NOINIT, .false.)
          end if
        end if
        
        ! Allocate memory for index pointer to node list if required
        if (rgroupFEMSet%h_InodeListIdx .eq. ST_NOHANDLE) then
          call storage_new('gfem_genNodeList', 'InodeListIdx',&
              rgroupFEMSet%NEQ+1, ST_INT, rgroupFEMSet%h_InodeListIdx,&
              ST_NEWBLOCK_NOINIT)
        else
          call storage_getsize(rgroupFEMSet%h_InodeListIdx, isize)
          if (isize .ne. rgroupFEMSet%NEQ+1) then
            call storage_realloc('gfem_genNodeList',&
                rgroupFEMSet%NEQ+1, ST_INT, rgroupFEMSet%h_InodeListIdx,&
                ST_NEWBLOCK_NOINIT)
          end if
        end if
        
        ! Set pointers
        call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
        call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList1D)

        ! Initialise counter and index pointer
        idx = 1; iidx = 1; p_InodeListIdx(1) = 1
        
        ! Loop over all equations
        do ieq = 1, rmatrix%NEQ
          
          ! Increase index counter
          iidx = iidx+1
          
          ! Loop over all matrix entries in current row
          do jcol = 1, rmatrix%NCOLS           
            
            ! Set column number and matrix position
            p_InodeList1D(idx) = jcol
            
            ! Increase counter
            idx = idx+1
          end do
          
          ! Set starting position of new node
          p_InodeListIdx(iidx) = idx
        end do

        ! Set state of structure
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec,&
                                    GFEM_HAS_NODESTRUCTURE)
        
      case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)

        ! Remove nodal structure from the group finite element set
        if (check(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE)) then
          call storage_free(rgroupFEMSet%h_InodeListIdx)
          call storage_free(rgroupFEMSet%h_InodeList)
        end if
        
        ! Set handles to matrix handles
        rgroupFEMSet%h_InodeListIdx = rmatrix%h_Kld
        rgroupFEMSet%h_InodeList    = rmatrix%h_Kcol
        
        ! Set state of structure
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec,&
                                    GFEM_HAS_NODESTRUCTURE)

        ! Reset ownership
        rgroupFEMSet%iduplicationFlag = ior(rgroupFEMSet%iduplicationFlag,&
                                            GFEM_SHARE_NODESTRUCTURE)

      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)

        ! Remove nodal structure from the group finite element set
        if (check(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE)) then
          call storage_free(rgroupFEMSet%h_InodeListIdx)
          call storage_free(rgroupFEMSet%h_InodeList)
        end if
        
        ! Set handles to matrix handles
        rgroupFEMSet%h_InodeListIdx = rmatrix%h_Kld
        rgroupFEMSet%h_InodeList    = rmatrix%h_Kcol
        
        ! Set state of structure
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec,&
                                    GFEM_HAS_NODESTRUCTURE)

        ! Reset ownership
        rgroupFEMSet%iduplicationFlag = ior(rgroupFEMSet%iduplicationFlag,&
                                            GFEM_SHARE_NODESTRUCTURE)

      end select
      
    else   ! use list of DOFs for restriction

      ! Set pointer
      call gfem_getbase_IdofList(rgroupFEMSet, p_IdofList)
      
      ! Allocate memory for node list if required
      if (rgroupFEMSet%h_InodeList .eq. ST_NOHANDLE) then
        Isize2D = (/2,rgroupFEMSet%NA/)
        call storage_new('gfem_genNodeList', 'InodeList',&
            Isize2D, ST_INT, rgroupFEMSet%h_InodeList,&
            ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(rgroupFEMSet%h_InodeList, Isize2D)
        if (Isize2D(2) .ne. rgroupFEMSet%NA) then
          call storage_free(rgroupFEMSet%h_InodeList)
          Isize2D = (/2,rgroupFEMSet%NA/)
          call storage_new('gfem_genNodeList', 'InodeList',&
              Isize2D, ST_INT, rgroupFEMSet%h_InodeList,&
              ST_NEWBLOCK_NOINIT)
        end if
      end if
      
      ! Allocate memory for index pointer to node list if required
      if (rgroupFEMSet%h_InodeListIdx .eq. ST_NOHANDLE) then
        call storage_new('gfem_genNodeList', 'InodeListIdx',&
            rgroupFEMSet%NEQ+1, ST_INT, rgroupFEMSet%h_InodeListIdx,&
            ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(rgroupFEMSet%h_InodeListIdx, isize)
        if (isize .ne. rgroupFEMSet%NEQ+1) then
          call storage_free(rgroupFEMSet%h_InodeListIdx)
          call storage_new('gfem_genNodeList', 'InodeListIdx',&
              rgroupFEMSet%NEQ+1, ST_INT, rgroupFEMSet%h_InodeListIdx,&
              ST_NEWBLOCK_NOINIT)
        end if
      end if

      ! Set pointers
      call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
      call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)

      ! General node data structure
      select case(rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX1)
        
        ! Generate set of active degrees of freedom
        allocate(BisActive(max(rmatrix%NEQ,rmatrix%NCOLS))); BisActive=.false.
        do idx = 1, size(p_IdofList)
          BisActive(p_IdofList(idx))=.true.
        end do

        ! Initialise counter and index pointer
        idx = 1; iidx = 1; p_InodeListIdx(1) = 1

        ! Loop over all equations
        do ieq = 1, rmatrix%NEQ
          
          ! Check if this row belongs to an active DOF
          if (.not.BisActive(ieq)) cycle

          ! Increase index counter
          iidx = iidx+1

          ! Loop over all matrix entries in current row
          do jcol = 1, rmatrix%NCOLS
            
            ! Check if this column belongs to an active DOF
            if (.not.BisActive(jcol)) cycle

            ! Set column number and matrix position
            p_InodeList(1,idx) = jcol
            p_InodeList(2,idx) = rmatrix%NCOLS*(ieq-1)+jcol
            
            ! Increase counter
            idx = idx+1
          end do

          ! Set starting position of new node
          p_InodeListIdx(iidx) = idx
        end do
               
        ! Deallocate temporal memory
        deallocate(BisActive)

        ! Set state of structure
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec,&
                                    GFEM_HAS_NODESTRUCTURE)

      case (LSYSSC_MATRIX7, LSYSSC_MATRIX7INTL)

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        
        ! Generate set of active degrees of freedom
        allocate(BisActive(rmatrix%NEQ)); BisActive=.false.
        do idx = 1, size(p_IdofList)
          BisActive(p_IdofList(idx))=.true.
        end do

        ! Initialise counter and index pointer
        idx = 1; iidx = 1; p_InodeListIdx(1) = 1
        
        ! Loop over all rows
        do ieq = 1, size(p_Kld)-1

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(ieq)) cycle

          ! Increase index counter
          iidx = iidx+1

          ! Loop over all matrix entries in current row
          do ij = p_Kld(ieq), p_Kld(ieq+1)-1

            ! Get column number
            jcol = p_Kcol(ij)

            ! Check if this column belongs to an active DOF
            if (.not.BisActive(jcol)) cycle

            ! Set column number and matrix position
            p_InodeList(1,idx) = jcol
            p_InodeList(2,idx) = ij
            
            ! Increase counter
            idx = idx+1
          end do

          ! Set starting position of new node
          p_InodeListIdx(iidx) = idx
        end do

        ! Deallocate temporal memory
        deallocate(BisActive)

        ! Set state of structure
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec,&
                                    GFEM_HAS_NODESTRUCTURE)

      case (LSYSSC_MATRIX9, LSYSSC_MATRIX9INTL)

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        
        ! Generate set of active degrees of freedom
        allocate(BisActive(rmatrix%NEQ)); BisActive=.false.
        do idx = 1, size(p_IdofList)
          BisActive(p_IdofList(idx))=.true.
        end do
     
        ! Initialise counter and index pointer
        idx = 1; iidx = 1; p_InodeListIdx(1) = 1
        
        ! Loop over all rows
        do ieq = 1, size(p_Kld)-1

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(ieq)) cycle

          ! Increase index counter
          iidx = iidx+1

          ! Set row number and matrix position of diagonal entry
          p_InodeList(1,idx) = p_Kcol(p_Kdiagonal(ieq))
          p_InodeList(2,idx) = p_Kdiagonal(ieq)

          ! Increase counter
          idx = idx+1
          
          ! Loop over all left off-diagonal matrix entries in current row
          do ij = p_Kld(ieq), p_Kdiagonal(ieq)-1

            ! Get column number
            jcol = p_Kcol(ij)

            ! Check if this column belongs to an active DOF
            if (.not.BisActive(jcol)) cycle

            ! Set column number and matrix position
            p_InodeList(1,idx) = jcol
            p_InodeList(2,idx) = ij

            ! Increase counter
            idx = idx+1
          end do

          ! Loop over all right off-diagonal matrix entries in current row
          do ij = p_Kdiagonal(ieq)+1, p_Kld(ieq+1)-1
            ! Get column number
            jcol = p_Kcol(ij)

            ! Check if this column belongs to an active DOF
            if (.not.BisActive(jcol)) cycle

            ! Set column number and matrix position
            p_InodeList(1,idx) = jcol
            p_InodeList(2,idx) = ij

            ! Increase counter
            idx = idx+1
          end do

          ! Set starting position of new node
          p_InodeListIdx(iidx) = idx
        end do

        ! Deallocate temporal memory
        deallocate(BisActive)

        ! Set state of structure
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec,&
                                    GFEM_HAS_NODESTRUCTURE)
   
      end select
      
    end if

  contains

    !**************************************************************
    ! Checks if idupFlag has all bits ibitfield set.
    
    pure function check(idupFlag, ibitfield)
      
      integer(I32), intent(in) :: idupFlag,ibitfield
      
      logical :: check
      
      check = (iand(idupFlag,ibitfield) .eq. ibitfield)
      
    end function check

  end subroutine gfem_genNodeList

  !*****************************************************************************

!<subroutine>

  subroutine gfem_genEdgeList(rmatrix, rgroupFEMSet)

!<description>
    ! This subroutine generates the list of edges which are
    ! characterised by their two endpoints (i,j) and the absolute
    ! position of matrix entries ij and ji.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<inputoutput>
    ! Group finite element structure
    type(t_groupFEMSet), intent(inout) :: rgroupFEMSet
!</inputoutput>
!</subroutine>

    ! local variables
    !$ integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx,p_IdofList

    ! Check if edge structure is owned by the structure
    if (iand(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_EDGESTRUCTURE) .eq.&
             GFEM_SHARE_EDGESTRUCTURE) then
      call output_line('Edge list is not owned by structure and '//&
          'therefore cannot be generated',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfem_genEdgeList')
      call sys_halt()
    end if

    ! Check if list of DOFs has been attached
    if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
      
      ! Check if matrix and group finite element set are compatible
      if ((rgroupFEMSet%NA    .ne. rmatrix%NA)  .or.&
          (rgroupFEMSet%NEQ   .ne. rmatrix%NEQ) .or.&
          (rgroupFEMSet%NEDGE .ne. (rmatrix%NA-rmatrix%NEQ)/2)) then
        call output_line('Matrix/group finite element set not compatible,'//&
            ' different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_genEdgeList')
        call sys_halt()
      end if
      
      ! Generate edge list for a matrix which is structurally symmetric,
      ! i.e. edge (i,j) exists if and only if edge (j,i) exists without
      ! storing the diagonal edges (i,i).
      call lsyssc_genEdgeList(rmatrix, rgroupFEMSet%h_IedgeList,&
          LSYSSC_EDGELIST_NODESANDPOS, .true., .true., rgroupFEMSet%NEDGE)

    else

      ! Set pointer
      call gfem_getbase_IdofList(rgroupFEMSet, p_IdofList)
      
      ! Generate edge list for a matrix which is structurally symmetric,
      ! i.e. edge (i,j) exists if and only if edge (j,i) exists without
      ! storing the diagonal edges (i,i).
      call lsyssc_genEdgeList(rmatrix, rgroupFEMSet%h_IedgeList,&
          LSYSSC_EDGELIST_NODESANDPOS, .true., .true., rgroupFEMSet%NEDGE,&
          p_IdofList)
      
    end if

    ! Allocate memory
    if (rgroupFEMSet%h_IedgeListIdx .eq. ST_NOHANDLE) then
      call storage_new('gfem_genEdgeList', 'IedgeListIdx',&
          2, ST_INT, rgroupFEMSet%h_IedgeListIdx, ST_NEWBLOCK_NOINIT)
    end if
    
    ! Set pointer to edge structure
    call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
    
    ! If no OpenMP is used, then all edges belong to the same
    ! group. Otherwise, the edges will be reordered below.
    p_IedgeListIdx    = rgroupFEMSet%NEDGE+1
    p_IedgeListIdx(1) = 1
    
    ! OpenMP-Extension: Perform edge-coloring to find groups of
    ! edges which can be processed in parallel, that is, the
    ! vertices of the edges in the group are all distinct
    !$ call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)
    !$ if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX1) then
    !$   call lsyssc_regroupEdgeList(rmatrix%NEQ, p_IedgeList,&
    !$       rgroupFEMSet%h_IedgeListIdx, 2*(rmatrix%NEQ-1))
    !$ else
    !$   call lsyssc_regroupEdgeList(rmatrix%NEQ, p_IedgeList,&
    !$       rgroupFEMSet%h_IedgeListIdx)
    !$ end if

    ! Set state of structure
    rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_EDGESTRUCTURE)
 
  end subroutine gfem_genEdgeList

  !*****************************************************************************

!<subroutine>

  subroutine gfem_genDiagList(rmatrix, rgroupFEMSet)

!<description>IdofList
    ! This subroutine stores the diagonal pointer
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<inputoutput>
    ! Group finite element structure
    type(t_groupFEMSet), intent(inout) :: rgroupFEMSet
!</inputoutput>
!</subroutine>

    ! local variables
    logical, dimension(:), allocatable :: BisActive
    integer, dimension(:,:), pointer :: p_IdiagList
    integer, dimension(:), pointer :: p_IdofList
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    integer, dimension(2) :: Isize
    integer :: idx,ieq
    
    ! Check if diagonal pointer is owned by the structure
    if (iand(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_DIAGLIST) .eq.&
        GFEM_SHARE_DIAGLIST) then
      call output_line('Diagonal pointer is not owned by structure and '//&
          'therefore cannot be generated',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfem_genDiagList')
      call sys_halt()
    end if
    
    ! Check if list of DOFs has been attached
    if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
      
      ! Check if matrix and group finite element set are compatible
      if ((rgroupFEMSet%NA    .ne. rmatrix%NA)  .or.&
          (rgroupFEMSet%NEQ   .ne. rmatrix%NEQ) .or.&
          (rgroupFEMSet%NEDGE .ne. (rmatrix%NA-rmatrix%NEQ)/2)) then
        call output_line('Matrix/group finite element set not compatible,'//&
            ' different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_genDiagList')
        call sys_halt()
      end if
      
      ! Allocate new memory for diagonal pointer if required
      if (rgroupFEMSet%h_IdiagList .eq. ST_NOHANDLE) then
        Isize = (/2,rgroupFEMSet%NEQ/)
        call storage_new('gfem_genDiagList', 'IdiagList',&
            Isize, ST_INT, rgroupFEMSet%h_IdiagList,&
            ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(rgroupFEMSet%h_IdiagList, Isize)
        if (Isize(2) .ne. rgroupFEMSet%NEQ) then
          call storage_realloc('gfem_genDiagList',&
              rgroupFEMSet%NEQ, ST_INT, rgroupFEMSet%h_IdiagList,&
              ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if

      ! Set pointer
      call gfem_getbase_IdiagList(rgroupFEMSet, p_IdiagList)
      
      ! What matrix type are we?
      select case(rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX1)
        
        ! Loop over all rows of the matrix
        p_IdiagList(1,1) = 1
        do ieq = 2, rmatrix%NEQ
          p_IdiagList(1,ieq) = ieq
          p_IdiagList(2,ieq) = p_IdiagList(2,ieq-1)+rmatrix%NCOLS
        end do
        
        ! Set specifier
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGLIST)
        
      case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        
        ! Set pointer
        call lsyssc_getbase_Kld(rmatrix, p_Kld)

        ! Loop over all rows of the matrix
        do ieq = 1, rmatrix%NEQ
          p_IdiagList(1,ieq) = ieq
          p_IdiagList(2,ieq) = p_Kld(ieq)
        end do
        
        ! Set specifier
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGLIST)
                
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
        
        ! Set pointer
        call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

        ! Loop over all rows of the matrix
        do ieq = 1, rmatrix%NEQ
          p_IdiagList(1,ieq) = ieq
          p_IdiagList(2,ieq) = p_Kdiagonal(ieq)
        end do
        
        ! Set specifier
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGLIST)
        
      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfem_initGFEMSetByMatrix')
        call sys_halt()
      end select

    else   ! use list of DOFs for restriction

      ! Allocate new memory for diagonal pointer if required
      if (rgroupFEMSet%h_IdiagList .eq. ST_NOHANDLE) then
        Isize = (/2,rgroupFEMSet%NEQ/)
        call storage_new('gfem_genDiagList', 'IdiagList',&
            Isize, ST_INT, rgroupFEMSet%h_IdiagList,&
            ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(rgroupFEMSet%h_IdiagList, Isize)
        if (Isize(2) .ne. rgroupFEMSet%NEQ) then
          call storage_realloc('gfem_genDiagList',&
              rgroupFEMSet%NEQ, ST_INT, rgroupFEMSet%h_IdiagList,&
              ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if

      ! Set pointer
      call gfem_getbase_IdiagList(rgroupFEMSet, p_IdiagList)
      call gfem_getbase_IdofList(rgroupFEMSet, p_IdofList)
      
      ! What matrix type are we?
      select case(rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX1)

        ! Generate set of active degrees of freedom
        allocate(BisActive(max(rmatrix%NEQ,rmatrix%NCOLS))); BisActive=.false.
        do idx = 1, size(p_IdofList)
          BisActive(p_IdofList(idx))=.true.
        end do

        ! Initialise index
        idx = 0

        ! Loop over all equations
        do ieq = 1, rmatrix%NEQ
          
          ! Check if this row belongs to an active DOF
          if (.not.BisActive(ieq)) cycle
          
          ! Increase index counter
          idx = idx+1
          
          ! Set position of diagonal entry
          p_IdiagList(1,idx) = ieq
          p_IdiagList(2,idx) = rmatrix%NCOLS*(ieq-1)+ieq
        end do

        ! Deallocate temporal memory
        deallocate(BisActive)

        ! Set specifier
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGLIST)

      case (LSYSSC_MATRIX7, LSYSSC_MATRIX7INTL,&
            LSYSSC_MATRIX9, LSYSSC_MATRIX9INTL)

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        if ((rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9) .or.&
            (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9INTL)) then
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        else
          call lsyssc_getbase_Kld(rmatrix, p_Kdiagonal)
        end if
        
        ! Generate set of active degrees of freedom
        allocate(BisActive(rmatrix%NEQ)); BisActive=.false.
        do idx = 1, size(p_IdofList)
          BisActive(p_IdofList(idx))=.true.
        end do

        ! Initialise index
        idx = 0
        
        ! Loop over all equations
        do ieq = 1, rmatrix%NEQ

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(ieq)) cycle

          ! Increase index counter
          idx = idx+1
          
          ! Set position of diagonal entry
          p_IdiagList(1,idx) = ieq
          p_IdiagList(2,idx) = p_Kdiagonal(ieq)
        end do
        
        ! Deallocate temporal memory
        deallocate(BisActive)
        
        ! Set specifier
        rgroupFEMSet%isetSpec = ior(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGLIST)
        
      end select
      
    end if
    
  end subroutine gfem_genDiagList

    !*****************************************************************************

!<subroutine>

  subroutine gfem_infoGroupFEMSet(rgroupFEMSet)

!<description>
    ! This subroutine prints out information about the group finite element set
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>
!</subroutine>

    call output_line('GroupFEMSet:')
    call output_line('------------')
    call output_line('cassemblyType:            '//trim(sys_siL(rgroupFEMSet%cassemblyType,15)))
    call output_line('cdataType:                '//trim(sys_siL(rgroupFEMSet%cdataType,15)))
    call output_line('iduplicationFlag:         '//trim(sys_siL(rgroupFEMSet%iduplicationFlag,15)))
    call checkAndOutput('GFEM_SHARE_DOFLIST:       ',rgroupFEMSet%iduplicationFlag,GFEM_SHARE_DOFLIST)
    call checkAndOutput('GFEM_SHARE_DIAGLIST:      ',rgroupFEMSet%iduplicationFlag,GFEM_SHARE_DIAGLIST)
    call checkAndOutput('GFEM_SHARE_NODESTRUCTURE: ',rgroupFEMSet%iduplicationFlag,GFEM_SHARE_NODESTRUCTURE)
    call checkAndOutput('GFEM_SHARE_EDGESTRUCTURE: ',rgroupFEMSet%iduplicationFlag,GFEM_SHARE_EDGESTRUCTURE)
    call checkAndOutput('GFEM_SHARE_DIAGDATA:      ',rgroupFEMSet%iduplicationFlag,GFEM_SHARE_DIAGDATA)
    call checkAndOutput('GFEM_SHARE_NODEDATA:      ',rgroupFEMSet%iduplicationFlag,GFEM_SHARE_NODEDATA)
    call checkAndOutput('GFEM_SHARE_EDGEDATA:      ',rgroupFEMSet%iduplicationFlag,GFEM_SHARE_EDGEDATA)
    call output_line('isetSpec:                 '//trim(sys_siL(rgroupFEMSet%isetSpec,15)))
    call checkAndOutput('GFEM_HAS_DOFLIST:         ',rgroupFEMSet%isetSpec,GFEM_HAS_DOFLIST)
    call checkAndOutput('GFEM_HAS_DIAGLIST:        ',rgroupFEMSet%isetSpec,GFEM_HAS_DIAGLIST)
    call checkAndOutput('GFEM_HAS_NODESTRUCTURE:   ',rgroupFEMSet%isetSpec,GFEM_HAS_NODESTRUCTURE)
    call checkAndOutput('GFEM_HAS_EDGESTRUCTURE:   ',rgroupFEMSet%isetSpec,GFEM_HAS_EDGESTRUCTURE)
    call checkAndOutput('GFEM_HAS_DIAGDATA:        ',rgroupFEMSet%isetSpec,GFEM_HAS_DIAGDATA)
    call checkAndOutput('GFEM_HAS_NODEDATA:        ',rgroupFEMSet%isetSpec,GFEM_HAS_NODEDATA)
    call checkAndOutput('GFEM_HAS_EDGEDATA:        ',rgroupFEMSet%isetSpec,GFEM_HAS_EDGEDATA)
    call output_line('NA:                       '//trim(sys_siL(rgroupFEMSet%NA,15)))
    call output_line('NEQ:                      '//trim(sys_siL(rgroupFEMSet%NEQ,15)))
    call output_line('NEDGE:                    '//trim(sys_siL(rgroupFEMSet%NEDGE,15)))
    call output_line('NVAR:                     '//trim(sys_siL(rgroupFEMSet%NVAR,15)))
    call output_line('ncoeffsAtDiag:            '//trim(sys_siL(rgroupFEMSet%ncoeffsAtDiag,15)))
    call output_line('ncoeffsAtNode:            '//trim(sys_siL(rgroupFEMSet%ncoeffsAtNode,15)))
    call output_line('ncoeffsAtEdge:            '//trim(sys_siL(rgroupFEMSet%ncoeffsAtEdge,15)))
    call checkAndOutputHandle('IdofList:                 ', rgroupFEMSet%h_IdofList)
    call checkAndOutputHandle('IdiagList:                ', rgroupFEMSet%h_IdiagList)
    call checkAndOutputHandle('IedgeListIdx:             ', rgroupFEMSet%h_IedgeListIdx)
    call checkAndOutputHandle('IedgeList:                ', rgroupFEMSet%h_IedgeList)
    call checkAndOutputHandle('InodeListIdx:             ', rgroupFEMSet%h_InodeListIdx)
    call checkAndOutputHandle('InodeList:                ', rgroupFEMSet%h_InodeList)
    call checkAndOutputHandle('CoeffsAtDiag:             ', rgroupFEMSet%h_CoeffsAtDiag)
    call checkAndOutputHandle('CoeffsAtNode:             ', rgroupFEMSet%h_CoeffsAtNode)
    call checkAndOutputHandle('CoeffsAtEdge:             ', rgroupFEMSet%h_CoeffsAtEdge)

  contains

    ! Here some working routines follow

    !***************************************************************************
    ! Check bitfield and output string

    subroutine checkAndOutput(cstring, iflag, ibitfield)

      ! input parameters
      character(len=*), intent(in) :: cstring
      integer(I32), intent(in) :: iflag, ibitfield

      if (iand(iflag, ibitfield) .eq. ibitfield) then
        call output_line (cstring//'TRUE')
      else
        call output_line (cstring//'FALSE')
      end if

    end subroutine checkAndOutput

    !***************************************************************************
    ! Check handle and output shape of array

    subroutine checkAndOutputHandle(cstring, h_handle)

      ! input parameters
      character(len=*), intent(in) :: cstring
      integer, intent(in) :: h_handle

      ! local variabels
      integer :: isize,idimension
      integer, dimension(2) :: Isize2D
      integer, dimension(3) :: Isize3D

      if (h_handle .ne. ST_NOHANDLE) then
        call storage_getdimension(h_handle, idimension)

        select case(idimension)
        case (1)
          call storage_getsize(h_handle, isize)
          call output_line (cstring//trim(sys_siL(h_handle,15))//&
              ' ('//trim(sys_siL(isize,15))//')')

        case(2)
          call storage_getsize(h_handle, Isize2D)
          call output_line (cstring//trim(sys_siL(h_handle,15))//&
              ' ('//trim(sys_siL(Isize2D(1),15))//','//trim(sys_siL(Isize2D(2),15))//')')
          
        case (3)
          call storage_getsize(h_handle, Isize3D)
          call output_line (cstring//trim(sys_siL(h_handle,15))//&
              ' ('//trim(sys_siL(Isize3D(1),15))//','//trim(sys_siL(Isize3D(2),15))//&
                    trim(sys_siL(Isize3D(3),15))//')')
        end select
      else
        call output_line (cstring//trim(sys_siL(h_handle,15)))
      end if
      
    end subroutine checkAndOutputHandle

  end subroutine gfem_infoGroupFEMSet

  !*****************************************************************************

!<subroutine>

  subroutine gfem_infoGroupFEMBlock(rgroupFEMBlock)

!<description>
    ! This subroutine prints out information about the block of group
    ! finite element sets
!</description>

!<input>
    ! Block of group finite element sets
    type(t_groupFEMBlock), intent(in) :: rgroupFEMBlock
!</input>
!</subroutine>

    ! local variables
    integer :: i

    call output_line('GroupFEMBlock:')
    call output_line('--------------')
    call output_line('nblocks: '//trim(sys_siL(rgroupFEMBlock%nblocks,15)))
    call output_lbrk()
    
    do i = 1, rgroupFEMBlock%nblocks
      call gfem_infoGroupFEMSet(rgroupFEMBlock%RgroupFEMBlock(i))
    end do
    
  end subroutine gfem_infoGroupFEMBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfem_copyH2D_IedgeList(rgroupFEMSet, btranspose)

!<description>
    ! This subroutine copies the edge structure from the host memory
    ! to the memory of the coprocessor device. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>


    if (rgroupFEMSet%h_IedgeList .ne. ST_NOHANDLE)&
        call storage_syncMemory(rgroupFEMSet%h_IedgeList,&
        ST_SYNCBLOCK_COPY_H2D, btranspose)

  end subroutine gfem_copyH2D_IedgeList

  !*****************************************************************************

!<subroutine>

  subroutine gfem_copyD2H_IedgeList(rgroupFEMSet, btranspose)

!<description>
    ! This subroutine copies the edge structure from the memory of the
    ! coprocessor device to the host memory. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>

    if (rgroupFEMSet%h_IedgeList .ne. ST_NOHANDLE)&
        call storage_syncMemory(rgroupFEMSet%h_IedgeList,&
        ST_SYNCBLOCK_COPY_D2H, btranspose)

  end subroutine gfem_copyD2H_IedgeList

  !*****************************************************************************

!<subroutine>

  subroutine gfem_copyH2D_CoeffsAtNode(rgroupFEMSet, btranspose)

!<description>
    ! This subroutine copies the coefficients at nodes from the host
    ! memory to the memory of the coprocessor device. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>


    if (rgroupFEMSet%h_CoeffsAtNode .ne. ST_NOHANDLE)&
        call storage_syncMemory(rgroupFEMSet%h_CoeffsAtNode,&
        ST_SYNCBLOCK_COPY_H2D, btranspose)

  end subroutine gfem_copyH2D_CoeffsAtNode

  !*****************************************************************************

!<subroutine>

  subroutine gfem_copyD2H_CoeffsAtNode(rgroupFEMSet, btranspose)

!<description>
    ! This subroutine copies the coefficients at nodes from the memory
    ! of the coprocessor device to the host memory. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>

    if (rgroupFEMSet%h_CoeffsAtNode .ne. ST_NOHANDLE)&
        call storage_syncMemory(rgroupFEMSet%h_CoeffsAtNode,&
        ST_SYNCBLOCK_COPY_D2H, btranspose)

  end subroutine gfem_copyD2H_CoeffsAtNode

  !*****************************************************************************

!<subroutine>

  subroutine gfem_copyH2D_CoeffsAtEdge(rgroupFEMSet, btranspose)

!<description>
    ! This subroutine copies the coefficients at edges from the host
    ! memory to the memory of the coprocessor device. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>


    if (rgroupFEMSet%h_CoeffsAtEdge .ne. ST_NOHANDLE)&
        call storage_syncMemory(rgroupFEMSet%h_CoeffsAtEdge,&
        ST_SYNCBLOCK_COPY_H2D, btranspose)
    
  end subroutine gfem_copyH2D_CoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine gfem_copyD2H_CoeffsAtEdge(rgroupFEMSet, btranspose)

!<description>
    ! This subroutine copies the coefficients at edges from the memory
    ! of the coprocessor device to the host memory. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>

    if (rgroupFEMSet%h_CoeffsAtEdge .ne. ST_NOHANDLE)&
        call storage_syncMemory(rgroupFEMSet%h_CoeffsAtEdge,&
        ST_SYNCBLOCK_COPY_D2H, btranspose)

  end subroutine gfem_copyD2H_CoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine gfem_copyH2D_CoeffsAtDiag(rgroupFEMSet, btranspose)

!<description>
    ! This subroutine copies the coefficients at matrix diagonals from
    ! the host memory to the memory of the coprocessor device. If no
    ! device is available, then an error is thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>


    if (rgroupFEMSet%h_CoeffsAtDiag .ne. ST_NOHANDLE)&
        call storage_syncMemory(rgroupFEMSet%h_CoeffsAtDiag,&
        ST_SYNCBLOCK_COPY_H2D, btranspose)
    
  end subroutine gfem_copyH2D_CoeffsAtDiag

  !*****************************************************************************

!<subroutine>

  subroutine gfem_copyD2H_CoeffsAtDiag(rgroupFEMSet, btranspose)

!<description>
    ! This subroutine copies the coefficients at matrix diagonals from
    ! the memory of the coprocessor device to the host memory. If no
    ! device is available, then an error is thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>

    if (rgroupFEMSet%h_CoeffsAtDiag .ne. ST_NOHANDLE)&
        call storage_syncMemory(rgroupFEMSet%h_CoeffsAtDiag,&
        ST_SYNCBLOCK_COPY_D2H, btranspose)

  end subroutine gfem_copyD2H_CoeffsAtDiag

  ! ***************************************************************************

!<subroutine>

  subroutine gfem_allocCoeffs(rgroupFEMSet, ncoeffsAtDiag, ncoeffsAtNode,&
      ncoeffsAtEdge, cdataType)

!<description>
    ! This subroutine initialises memory for storing precomputed
    ! coefficients at nodes and/or edges.
!</description>

!<input>
    ! Number of precomputed coefficients.
    ! If a value is negative, then the corresponding array is not modified!
    integer, intent(in) :: ncoeffsAtDiag
    integer, intent(in) :: ncoeffsAtNode
    integer, intent(in) :: ncoeffsAtEdge

    ! OPTIONAL Data type
    integer, intent(in), optional :: cdataType
!</input>

!<inputoutput>
    ! Group finite element set
    type(t_groupFEMSet), intent(inout) :: rgroupFEMSet
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize2D
    integer, dimension(3) :: Isize3D

    ! Remove auxiliary memory if not shared with others and set new dimension
    if (ncoeffsAtDiag .eq. 0) then
      if ((iand(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_DIAGDATA) .eq. 0)&
          .and.(rgroupFEMSet%h_CoeffsAtDiag .ne. ST_NOHANDLE))&
          call storage_free(rgroupFEMSet%h_CoeffsAtDiag)
      rgroupFEMSet%ncoeffsAtDiag = ncoeffsAtDiag
    end if

    if (ncoeffsAtNode .ge. 0) then
      if ((iand(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_NODEDATA) .eq. 0)&
          .and.(rgroupFEMSet%h_CoeffsAtNode .ne. ST_NOHANDLE))&
          call storage_free(rgroupFEMSet%h_CoeffsAtNode)
      rgroupFEMSet%ncoeffsAtNode = ncoeffsAtNode
    end if
      
    if (ncoeffsAtEdge .ge. 0) then
      if ((iand(rgroupFEMSet%iduplicationFlag, GFEM_SHARE_EDGEDATA) .eq. 0)&
        .and.(rgroupFEMSet%h_CoeffsAtEdge .ne. ST_NOHANDLE))&
        call storage_free(rgroupFEMSet%h_CoeffsAtEdge)
      rgroupFEMSet%ncoeffsAtEdge = ncoeffsAtEdge
    end if

    ! Set data type
    if (present(cdataType)) rgroupFEMSet%cdataType = cdataType
    
    if (ncoeffsAtDiag .gt. 0) rgroupFEMSet%ncoeffsAtDiag = ncoeffsAtDiag
    if (ncoeffsAtNode .gt. 0) rgroupFEMSet%ncoeffsAtNode = ncoeffsAtNode
    if (ncoeffsAtEdge .gt. 0) rgroupFEMSet%ncoeffsAtEdge = ncoeffsAtEdge

    ! Allocate diagonal data array
    if ((rgroupFEMSet%NEQ .gt. 0) .and. (rgroupFEMSet%ncoeffsAtDiag .gt. 0)) then
      Isize2D = (/rgroupFEMSet%ncoeffsAtDiag, rgroupFEMSet%NEQ/)
      call storage_new('gfem_allocCoeffs', 'CoeffsAtDiag',&
          Isize2D, rgroupFEMSet%cdataType, rgroupFEMSet%h_CoeffsAtDiag,&
          ST_NEWBLOCK_NOINIT)
    end if
    
    ! Allocate nodal data arary
    if ((rgroupFEMSet%NA .gt. 0) .and. (rgroupFEMSet%ncoeffsAtNode .gt. 0)) then
      Isize2D = (/rgroupFEMSet%ncoeffsAtNode, rgroupFEMSet%NA/)
      call storage_new('gfem_allocCoeffs', 'CoeffsAtNode',&
          Isize2D, rgroupFEMSet%cdataType, rgroupFEMSet%h_CoeffsAtNode,&
          ST_NEWBLOCK_NOINIT)
    end if

    ! Allocate edge-based data array
    if ((rgroupFEMSet%NEDGE .gt. 0) .and. (rgroupFEMSet%ncoeffsAtEdge .gt. 0)) then
      Isize3D = (/rgroupFEMSet%ncoeffsAtEdge, 2, rgroupFEMSet%NEDGE/)  
      call storage_new('gfem_allocCoeffs', 'CoeffsAtEdge',&
          Isize3D, rgroupFEMSet%cdataType, rgroupFEMSet%h_CoeffsAtEdge,&
          ST_NEWBLOCK_NOINIT)
    end if
    
  end subroutine gfem_allocCoeffs

end module groupfembase
