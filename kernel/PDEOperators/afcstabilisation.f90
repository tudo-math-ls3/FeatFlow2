!##############################################################################
!# ****************************************************************************
!# <name> afcstabilisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines for performing
!# discrete stabilisation by means of algebraic flux correction
!#
!# The following routines are available:
!#
!# 1.) afcstab_initFromParameterlist
!#     -> create a stabilisation structure and initialize
!#        it from the values of a given parameter list
!#
!# 2.) afcstab_releaseStabilisation
!#     -> release a stabilisation structure
!#
!# 3.) afcstab_resizeStabilisation = afcstab_resizeStabDirect /
!#                                   afcstab_resizeStabIndScalar
!#                                   afcstab_resizeStabIndBlock
!#     -> resize a stabilisation structure
!#
!# 4.) afcstab_getbase_IsupdiagEdgeIdx
!#     -> return pointer to the index pointer for the
!#        superdiagonal edge numbers
!#
!# 5.) afcstab_getbase_IverticesAtEdge
!#     -> return pointer to the vertices at edge structure
!#
!# 6.) afcstab_getbase_IsubdiagEdgeIdx
!#     -> return pointer to the index pointer for the
!#        subdiagonal edge numbers
!#
!# 7.) afcstab_getbase_IsubdiagEdge
!#     -> return pointer to the subdiagonal edge numbers
!#
!# 8.) afcstab_getbase_DcoeffsAtEdge
!#     -> return pointer to edge data
!#
!# 9.) afcstab_generateSubdiagEdges
!#      -> generate the subdiagonal edge data structure
!#
!# 10.) afcstab_generateExtSparsity
!#      -> generate the extended sparsity pattern
!#
!# </purpose>
!##############################################################################
module afcstabilisation

  use fsystem
  use genoutput
  use linearsystemscalar
  use linearsystemblock
  use paramlist
  use storage
  use triangulation

  implicit none
  
  private
  public :: t_afcstab
  public :: afcstab_initFromParameterlist
  public :: afcstab_releaseStabilisation
  public :: afcstab_resizeStabilisation
  public :: afcstab_getbase_IsupdiagEdgeIdx
  public :: afcstab_getbase_IverticesAtEdge
  public :: afcstab_getbase_DcoeffsAtEdge
  public :: afcstab_getbase_IsubdiagEdgeIdx
  public :: afcstab_getbase_IsubdiagEdge
  public :: afcstab_generateSubdiagEdges
  public :: afcstab_generateExtSparsity
  public :: afcstab_limit
 
  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>
!<constantblock description="Global format flags for AFC stabilisation">

  ! No stabilisation: use standard high-order Galerkin discretisation
  integer, parameter, public :: AFCSTAB_GALERKIN          = 0
  
  ! Stabilisation of discrete upwind type for convection operators
  integer, parameter, public :: AFCSTAB_UPWIND            = 1

  ! Stabilisation of discrete maximum principle preserving 
  ! type for anisotropic diffusion operators
  integer, parameter, public :: AFCSTAB_DMP               = 2

  ! Stabilisation of semi-implicit FEM-FCT type for convection operators
  integer, parameter, public :: AFCSTAB_FEMFCT            = 10

  ! Stabilisation of semi-explicit (classical) FEM-FCT type for convection operators
  integer, parameter, public :: AFCSTAB_FEMFCT_CLASSICAL  = 11

  ! Stabilisation of linearised FEM-FCT type for convection operators
  integer, parameter, public :: AFCSTAB_FEMFCT_LINEARIZED = 12
  
  ! Stabilisation of FEM-TVD type for convection operators
  integer, parameter, public :: AFCSTAB_FEMTVD            = 20

  ! Stabilisation of general purpose type for convection operators
  integer, parameter, public :: AFCSTAB_FEMGP             = 21
  
  ! Stabilisation of symmetric type for diffusion operators
  integer, parameter, public :: AFCSTAB_SYMMETRIC         = 30
  
!</constantblock>

!<constantblock description="Bitfield identifiers for state of stabilisation">
  
  ! Stabilisation is undefined
  integer, parameter, public :: AFCSTAB_UNDEFINED         = 2**0

  ! Stabilisation has been initialised
  integer, parameter, public :: AFCSTAB_INITIALISED       = 2**1

  ! Edge-based structure generated: KEDGE
  integer, parameter, public :: AFCSTAB_EDGESTRUCTURE     = 2**2

  ! Edge-based structure oriented: KEDGE
  integer, parameter, public :: AFCSTAB_EDGEORIENTATION   = 2**3

  ! Edge-based values computed from matrix: DEDGE
  integer, parameter, public :: AFCSTAB_EDGEVALUES        = 2**4

  ! Nodal antidiffusion: PP,PM
  integer, parameter, public :: AFCSTAB_ANTIDIFFUSION     = 2**5

  ! Nodal upper/lower bounds: QP,QM
  integer, parameter, public :: AFCSTAB_BOUNDS            = 2**6
  
  ! Nodal correction factors computed: RP,RM
  integer, parameter, public :: AFCSTAB_LIMITER           = 2**7

  ! Antidiffusive fluxes precomputed
  integer, parameter, public :: AFCSTAB_FLUXES            = 2**8
  
  ! Subdiagonal edge-based structure generated
  integer, parameter, public :: AFCSTAB_SUBDIAGONALEDGES  = 2**9
!</constantblock>

!</constants>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<types>
!<typeblock>

  ! data structure that holds all required information for stabilisation
  type t_afcstab
    
    ! Format Tag. Identifies the type of stabilisation
    integer :: ctypeAFCstabilisation                   = AFCSTAB_GALERKIN

    ! Format Tag: Specifies the stabilisation
    integer :: iSpec                                   = AFCSTAB_UNDEFINED

    ! Number of equations of the sparsity pattern
    integer(PREC_VECIDX) :: NEQ                        = 0

    ! Number of local variables; in general scalar solution vectors of
    ! size NEQ posses NEQ entries. However, scalar vectors can be interleaved,
    ! that is, each of the NEQ entries stores NVAR local variables. In this case,
    ! NEQ remains unmodified but NVAR>1 such that the physical length of the
    ! vector is NEQ*NVAR.
    integer :: NVAR                                    = 1

    ! Number of edges of the sparsity pattern
    integer(PREC_VECIDX) :: NEDGE                      = 0

    ! Maximum number of edges adjacent to one vertex. 
    ! This corresponds to the maximum number of nonzero row entries.
    integer(PREC_VECIDX) :: NNVEDGE                    = 0

    ! Handle to index pointer for superdiagonal edge numbers
    ! INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: IsuperdiagonalEdgesIdx
    ! The numbers IsuperdiagonalEdgesIdx(i):IsuperdiagonalEdgesIdx(i+1)-1
    ! denote the edge numbers of the ith vertex which are located in
    ! the upper right triangular matrix.
    integer :: h_IsuperdiagonalEdgesIdx                    = ST_NOHANDLE

    ! Handle to vertices at edge structure
    ! INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: IverticesAtEdge
    ! IverticesAtEdge(1:2,1:NEDGE) : the two end-points of the edge
    ! IverticesAtEdge(3:4,1:NEDGE) : the two matrix position that
    !                                correspond to the edge
    integer :: h_IverticesAtEdge                       = ST_NOHANDLE

    ! Handle to index pointer for subdiagonal edge numbers
    ! INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: IsubdiagonalEdgesIdx
    integer :: h_IsubdiagonalEdgesIdx                  = ST_NOHANDLE

    ! Handle to the subdiagonal edge numbers
    ! INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: IsubdiagonalEdges
    integer :: h_IsubdiagonalEdges                     = ST_NOHANDLE

    ! Handle to coefficient at edge structure
    integer :: h_DcoefficientsAtEdge                   = ST_NOHANDLE

    ! Flag whether or not the matrix is resorted.
    !  <0: Matrix is unsorted, sorting strategy is prepared in 
    !      h_IsortPermutation for a possible resorting of the entries.
    !  =0: Matrix is unsorted, no sorting strategy attached.
    !  >0: Matrix is sorted according to a sorting strategy.
    ! If <> 0, the absolute value 
    !               |isortStrategy| > 0
    ! indicates the sorting strategy to use, while the sign
    ! indicates whether the sorting strategy is active on the
    ! matrix (+) or not (-).
    ! The value is usually one of the SSTRAT_xxxx constants from
    ! the module 'sortstrategy'.
    integer :: isortStrategy                           = 0

    ! Handle to renumbering strategy for resorting the matrix.
    ! The renumbering strategy is a vector
    !   array [1..2*NEQ] of integer
    ! The first NEQ entries (1..NEQ) represent the permutation how to
    ! sort an unsorted matrix. The second NEQ entries (NEQ+1..2*NEQ)
    ! represent the inverse permutation.
    ! Looking from another viewpoint with the background of how a matrix
    ! is renumbered, one can say:
    !  p_IsortPermutation (column in sorted matrix) = column in unsorted matrix.
    !  p_IsortPermutation (NEQ+column in unsorted matrix) = column in sorted matrix.
    ! Whether or not the matrix is actually sorted depends on the
    ! flag isortStrategy!
    integer :: h_IsortPermutation                      = ST_NOHANDLE

    ! Auxiliary nodal vectors; used internally
    type(t_vectorScalar), dimension(:), pointer :: RnodalVectors      => null()

    ! Auxiliary nodal block vectors; used internally
    type(t_vectorBlock), dimension(:), pointer  :: RnodalBlockVectors => null()

    ! Auxiliary edge vectors; used internally
    type(t_vectorScalar), dimension(:), pointer :: RedgeVectors       => null()
  end type t_afcstab
!</typeblock>
!</types>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface afcstab_resizeStabilisation
    module procedure afcstab_resizeStabDirect
    module procedure afcstab_resizeStabIndScalar
    module procedure afcstab_resizeStabIndBlock
  end interface

  interface afcstab_limit
    module procedure afcstab_limit_unbounded
    module procedure afcstab_limit_bounded
  end interface

contains

  ! ***************************************************************************

!<subroutine>

  subroutine afcstab_initFromParameterlist(rparlist, ssectionName, rafcstab)

!<description>
    ! This subroutine creates a stabilisation structure and initializes
    ! its values from a given parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN)    :: rparlist

    ! Section name of the parameter list
    character(LEN=*), intent(IN)   :: ssectionName
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: istabilisation

    ! Get type of stabilisation from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        "istabilisation", istabilisation)
    
    ! Check if stabilisation should be applied
    if (istabilisation .eq. AFCSTAB_GALERKIN) then
      
      return   ! -> high-order Galerkin
      
    elseif ((istabilisation .ne. AFCSTAB_UPWIND)            .and. &
        (    istabilisation .ne. AFCSTAB_FEMFCT)            .and. &
        (    istabilisation .ne. AFCSTAB_FEMFCT_CLASSICAL)  .and. &
        (    istabilisation .ne. AFCSTAB_FEMFCT_LINEARIZED) .and. &
        (    istabilisation .ne. AFCSTAB_FEMTVD)            .and. &
        (    istabilisation .ne. AFCSTAB_FEMGP)             .and. &
        (    istabilisation .ne. AFCSTAB_DMP)               .and. &
        (    istabilisation .ne. AFCSTAB_SYMMETRIC)) then 
      
      call output_line('Invalid AFC type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_initFromParameterlist')
      call sys_halt()

    else
      ! Set type of stabilisation
      rafcstab%iSpec                 = AFCSTAB_UNDEFINED
      rafcstab%ctypeAFCstabilisation = istabilisation
      
    end if
  end subroutine afcstab_initFromParameterlist

  !*****************************************************************************

!<subroutine>
  
  subroutine afcstab_releaseStabilisation(rafcstab)

!<description>
    ! This subroutine releases a stabilisation structure
!</description>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! Free storage
    if (rafcstab%h_IsuperdiagonalEdgesIdx .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IsuperdiagonalEdgesIdx)
    if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IverticesAtEdge)
    if (rafcstab%h_IsubdiagonalEdgesIdx .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IsubdiagonalEdgesIdx)
    if (rafcstab%h_IsubdiagonalEdges .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IsubdiagonalEdges)
    if (rafcstab%h_DcoefficientsAtEdge .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_DcoefficientsAtEdge)
    
    ! Reset atomic data
    rafcstab%ctypeAFCstabilisation = AFCSTAB_GALERKIN
    rafcstab%iSpec                 = AFCSTAB_UNDEFINED
    rafcstab%NEQ                   = 0
    rafcstab%NVAR                  = 1
    rafcstab%NEDGE                 = 0
    rafcstab%NNVEDGE               = 0
    rafcstab%isortStrategy         = 0
    rafcstab%h_IsortPermutation    = ST_NOHANDLE

    ! Release auxiliary nodal vectors
    if (associated(rafcstab%RnodalVectors)) then
      do i = lbound(rafcstab%RnodalVectors,1),&
             ubound(rafcstab%RnodalVectors,1)
        call lsyssc_releaseVector(rafcstab%RnodalVectors(i))
      end do
      deallocate(rafcstab%RnodalVectors)
    end if

    ! Release auxiliary nodal block vectors
    if (associated(rafcstab%RnodalBlockVectors)) then
      do i = lbound(rafcstab%RnodalBlockVectors,1),&
             ubound(rafcstab%RnodalBlockVectors,1)
        call lsysbl_releaseVector(rafcstab%RnodalBlockVectors(i))
      end do
      deallocate(rafcstab%RnodalBlockVectors)
    end if

    ! Release auxiliary edge vectors
    if (associated(rafcstab%RedgeVectors)) then
      do i = lbound(rafcstab%RedgeVectors,1),&
             ubound(rafcstab%RedgeVectors,1)
        call lsyssc_releaseVector(rafcstab%RedgeVectors(i))
      end do
      deallocate(rafcstab%RedgeVectors)
    end if
  end subroutine afcstab_releaseStabilisation

   !*****************************************************************************

!<subroutine>

  subroutine afcstab_resizeStabDirect(rafcstab, neq, nedge)

!<description>
    ! This subroutine resizes all vectors of the stabilisation
    ! structure to the new values NEQ and NEDGE
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! number of equations
    integer(PREC_VECIDX), intent(IN) :: neq

    ! number of edges
    integer(PREC_VECIDX), intent(IN) :: nedge   
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    
    ! Resize nodal quantities
    if (rafcstab%NEQ .ne. neq) then

      ! Set new number of nodes
      rafcstab%NEQ = neq
      
      ! Resize edge index vector
      if (rafcstab%h_IsuperdiagonalEdgesIdx .ne. ST_NOHANDLE) then
        call storage_realloc('afcstab_resizeStabDirect',&
            rafcstab%NEQ+1, rafcstab%h_IsuperdiagonalEdgesIdx,&
            ST_NEWBLOCK_NOINIT, .false.)
      end if
      
      ! Resize subdiagonal edge index vector
      if (rafcstab%h_IsubdiagonalEdgesIdx .ne. ST_NOHANDLE) then
        call storage_realloc('afcstab_resizeStabDirect',&
            rafcstab%NEQ+1, rafcstab%h_IsubdiagonalEdgesIdx,&
            ST_NEWBLOCK_NOINIT, .false.)
      end if

      ! Resize auxiliary nodal vectors
      if(associated(rafcstab%RnodalVectors)) then
        do i = lbound(rafcstab%RnodalVectors,1),&
               ubound(rafcstab%RnodalVectors,1)
          call lsyssc_resizeVector(rafcstab%RnodalVectors(i),&
              rafcstab%NEQ, .false., .false.)
        end do
      end if

      ! Resize auxiliary nodal vectors
      if(associated(rafcstab%RnodalBlockVectors)) then
        do i = lbound(rafcstab%RnodalBlockVectors,1),&
               ubound(rafcstab%RnodalBlockVectors,1)
          call lsysbl_resizeVectorBlock(rafcstab%RnodalBlockVectors(i),&
              rafcstab%NEQ, .false., .false.)
        end do
      end if
    end if


    ! Resize edge quantities
    if (rafcstab%NEDGE .ne. nedge) then

      ! Set new number of edges
      rafcstab%NEDGE = nedge

      ! Resize array of edges
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE) then
        call storage_realloc('afcstab_resizeStabDirect',&
            rafcstab%NEDGE, rafcstab%h_IverticesAtEdge,&
            ST_NEWBLOCK_NOINIT, .false.)
      end if

      ! Resize array of subdiagonal edges
      if (rafcstab%h_IsubdiagonalEdges .ne. ST_NOHANDLE) then
        call storage_realloc('afcstab_resizeStabDirect',&
            rafcstab%NEDGE, rafcstab%h_IsubdiagonalEdges,&
            ST_NEWBLOCK_NOINIT, .false.)
      end if

      ! Resize array of edge data
      if (rafcstab%h_DcoefficientsAtEdge .ne. ST_NOHANDLE) then
        call storage_realloc('afcstab_resizeStabDirect',&
            rafcstab%NEDGE, rafcstab%h_DcoefficientsAtEdge,&
            ST_NEWBLOCK_NOINIT, .false.)
      end if

      ! Resize auxiliary edge vectors
      if(associated(rafcstab%RedgeVectors)) then
        do i = lbound(rafcstab%RedgeVectors,1),&
               ubound(rafcstab%RedgeVectors,1)
          call lsyssc_resizeVector(rafcstab%RedgeVectors(i),&
              rafcstab%NEDGE, .false., .false.)
        end do
      end if
    end if
  end subroutine afcstab_resizeStabDirect

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_resizeStabIndScalar(rafcstab, rmatrixTemplate)

!<description>
    ! This subroutine resizes all vectors of the stabilisation
    ! structure so that they are compatible to the template matrix.
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! template matrix
    type(t_matrixScalar), intent(IN) :: rmatrixTemplate
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: neq, nedge

    
    ! Determine number of equations and edges
    neq   = rmatrixTemplate%NEQ
    nedge = int(0.5*(rmatrixTemplate%NA-rmatrixTemplate%NEQ))

    ! Call resize routine directly
    call afcstab_resizeStabDirect(rafcstab, neq, nedge)

  end subroutine afcstab_resizeStabIndScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_resizeStabIndBlock(rafcstab, rmatrixBlockTemplate)

!<description>
    ! This subroutine resizes all vectors of the stabilisation
    ! structure so that they are compatible to the template matrix.
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! template matrix
    type(t_matrixBlock), intent(IN) :: rmatrixBlockTemplate
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(INOUT)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: neq, nedge

    
    ! Check if block matrix has only one block
    if ((rmatrixBlockTemplate%nblocksPerCol .eq. 1) .and. &
        (rmatrixBlockTemplate%nblocksPerRow .eq. 1)) then
      call afcstab_resizeStabIndScalar(rafcstab,&
          rmatrixBlockTemplate%RmatrixBlock(1,1))

      ! That's it
      return
    end if
    

    ! Determine number of equations and edges
    neq   = rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ
    nedge = int(0.5*(rmatrixBlockTemplate%RmatrixBlock(1,1)%NA-&
                     rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ))

    ! Call resize routine directly
    call afcstab_resizeStabDirect(rafcstab, neq, nedge)

  end subroutine afcstab_resizeStabIndBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)

!<description>
    ! Returns a pointer to the index pointer for vertices at edge structure
!</description>

!<input>
    ! discrete operator
    type(t_afcstab), intent(IN) :: rafcstab
!</input>

!<output>
    ! Pointer to the index pointer for superdiagonal edge numbers.
    ! NULL() if the discrete operator does not provide it.
    integer(PREC_VECIDX), dimension(:), pointer :: p_IsuperdiagonalEdgesIdx
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    if ((rafcstab%h_IsuperdiagonalEdgesIdx .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEQ                  .eq. 0)) then
      nullify(p_IsuperdiagonalEdgesIdx)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rafcstab%h_IsuperdiagonalEdgesIdx,&
        p_IsuperdiagonalEdgesIdx,rafcstab%NEQ+1)
  end subroutine afcstab_getbase_IsupdiagEdgeIdx

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

!<description>
    ! Returns a pointer to the vertices at edge structure
!</description>

!<input>
    ! discrete operator
    type(t_afcstab), intent(IN) :: rafcstab
!</input>

!<output>
    ! Pointer to the vertices at edge structure
    ! NULL() if the discrete operator does not provide it.
    integer(PREC_MATIDX), dimension(:,:), pointer :: p_IverticesAtEdge
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    if ((rafcstab%h_IverticesAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE             .eq. 0)) then
      nullify(p_IverticesAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_int2D(rafcstab%h_IverticesAtEdge,&
        p_IverticesAtEdge,rafcstab%NEDGE)
  end subroutine afcstab_getbase_IverticesAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)

!<description>
    ! Returns a pointer to the index pointer for 
    ! subdiagonal edge numbers
!</description>

!<input>
    ! discrete operator
    type(t_afcstab), intent(IN) :: rafcstab
!</input>

!<output>
    ! Pointer to the index pointer for 
    ! subdiagonal edge numbers
    ! NULL() if the discrete operator does not provide it.
    integer(PREC_MATIDX), dimension(:), pointer :: p_IsubdiagonalEdgesIdx
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    if ((rafcstab%h_IsubdiagonalEdgesIdx .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEQ                    .eq. 0)) then
      nullify(p_IsubdiagonalEdgesIdx)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rafcstab%h_IsubdiagonalEdgesIdx,&
        p_IsubdiagonalEdgesIdx,rafcstab%NEQ+1)
  end subroutine afcstab_getbase_IsubdiagEdgeIdx

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IsubdiagEdge(rafcstab,p_IsubdiagonalEdges)

!<description>
    ! Returns a pointer to the subdiagonal edge number
!</description>

!<input>
    ! discrete operator
    type(t_afcstab), intent(IN) :: rafcstab
!</input>

!<output>
    ! Pointer to the subdiagonal edge numbers
    ! NULL() if the discrete operator does not provide it.
    integer(PREC_MATIDX), dimension(:), pointer :: p_IsubdiagonalEdges
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    if ((rafcstab%h_IsubdiagonalEdges .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE               .eq. 0)) then
      nullify(p_IsubdiagonalEdges)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rafcstab%h_IsubdiagonalEdges,&
        p_IsubdiagonalEdges,rafcstab%NEDGE)
  end subroutine afcstab_getbase_IsubdiagEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_DcoeffsAtEdge(rafcstab,p_DcoefficientsAtEdge)

!<description>
    ! Returns a pointer to the double-valued edge data
!</description>

!<input>
    ! discrete operator
    type(t_afcstab), intent(IN) :: rafcstab
!</input>

!<output>
    ! Pointer to the double-valued edge data
    ! NULL() if the discrete operator does not provide it.
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    if ((rafcstab%h_DcoefficientsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE                .eq. 0)) then
      nullify(p_DcoefficientsAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_double2D(rafcstab%h_DcoefficientsAtEdge,&
        p_DcoefficientsAtEdge,rafcstab%NEDGE)
  end subroutine afcstab_getbase_DcoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_generateSubdiagEdges(rafcstab)

!<description>
    ! This subroutine generates the subdiagonal edge number
    ! structure based on a given edge data structure.
!</description>

!<inputoutput>
    ! discrete operator
    type(t_afcstab), intent(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer(PREC_MATIDX), dimension(:,:), pointer :: p_IverticesAtEdge
    integer(PREC_VECIDX), dimension(:), pointer   :: p_IsuperdiagonalEdgesIdx
    integer(PREC_MATIDX), dimension(:), pointer   :: p_IsubdiagonalEdges
    integer(PREC_VECIDX), dimension(:), pointer   :: p_IsubdiagonalEdgesIdx
    integer(PREC_MATIDX) :: iedge,nedge,istor
    integer(PREC_VECIDX) :: ieq,jeq,neq
    integer(I32)         :: isize

    ! Check if edge-based data structure is prepared
    if (iand(rafcstab%iSpec,AFCSTAB_EDGESTRUCTURE) .eq. 0) then
      call output_line('Discrete operator does not provide required &
          &edge-based data structure',OU_CLASS_ERROR,OU_MODE_STD,&
          'afcstab_generateSubdiagEdge')
      call sys_halt()
    end if

    ! store dimensions of discrete operator
    neq   = rafcstab%NEQ
    nedge = rafcstab%NEDGE

    ! Allocate memory (if required)
    if (rafcstab%h_IsubdiagonalEdgesIdx .eq. ST_NOHANDLE) then
      call storage_new('afcstab_generateSubdiagEdges','IsubdiagonalEdgesIdx',&
          neq+1,ST_INT,rafcstab%h_IsubdiagonalEdgesIdx,ST_NEWBLOCK_ZERO)
    else
      call storage_getsize(rafcstab%h_IsubdiagonalEdgesIdx,isize)
      if (isize < neq+1) then
        call storage_realloc('afcstab_generateSubdiagEdges',neq+1,&
            rafcstab%h_IsubdiagonalEdgesIdx,ST_NEWBLOCK_ZERO,.false.)
      else
        call storage_clear(rafcstab%h_IsubdiagonalEdgesIdx)
      end if
    end if

    if (rafcstab%h_IsubdiagonalEdges .eq. ST_NOHANDLE) then
      call storage_new('afcstab_generateSubdiagEdges','IsubdiagonalEdges',&
          nedge,ST_INT,rafcstab%h_IsubdiagonalEdges,ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rafcstab%h_IsubdiagonalEdges,isize)
      if (isize < nedge) then
        call storage_realloc('afcstab_generateSubdiagEdges',nedge,&
            rafcstab%h_IsubdiagonalEdges,ST_NEWBLOCK_NOINIT,.false.)
      end if
    end if
    
    ! Set pointers
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab,p_IsuperdiagonalEdgesIdx)
    call afcstab_getbase_IverticesAtEdge(rafcstab,p_IverticesAtEdge)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab,p_IsubdiagonalEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab,p_IsubdiagonalEdges)
    
    ! Count number of superdiagonal edges
    do ieq = 1, neq
      do iedge = p_IsuperdiagonalEdgesIdx(ieq),&
                 p_IsuperdiagonalEdgesIdx(ieq+1)-1

        ! Determine that end-point of the edge which is not equal to IEQ
        jeq = p_IverticesAtEdge(1,iedge) + p_IverticesAtEdge(2,iedge) - ieq + 1

        ! Increase number of edges connected to this point by one
        p_IsubdiagonalEdgesIdx(jeq) = p_IsubdiagonalEdgesIdx(jeq)+1
      end do
    end do
    
    ! Reshuffle pass 1.
    do ieq = 2, neq
      
      
      p_IsubdiagonalEdgesIdx(ieq) = p_IsubdiagonalEdgesIdx(ieq) +&
          p_IsubdiagonalEdgesIdx(ieq-1)
      
    end do
    p_IsubdiagonalEdgesIdx(neq+1) = p_IsubdiagonalEdgesIdx(neq+1) +&
        p_IsubdiagonalEdgesIdx(neq)
    
    ! Store the subdiagonal edge numbers. In addition, set the maximum
    ! number of edges adjacent to one vertex. That is, take the maximum
    ! value of the number of edges to the left of the diagonal plus 
    ! the number of edges to the right of the diagonal.

    rafcstab%NNVEDGE = 0! p_IsuperdiagonalEdgesIdx(2) - p_IsuperdiagonalEdgesIdx(1)

    do ieq = 1, neq
      
      ! For each equation, loop over the edges which are located in
      ! the upper right triangular matrix
      do iedge = p_IsuperdiagonalEdgesIdx(ieq),&
                 p_IsuperdiagonalEdgesIdx(ieq+1)-1
        
        ! Determine that end-point of the edge which is not equal to IEQ
        jeq = p_IverticesAtEdge(1,iedge) + p_IverticesAtEdge(2,iedge) - ieq
        
        ! Determine next free position
        istor = p_IsubdiagonalEdgesIdx(jeq)+1
        
        p_IsubdiagonalEdgesIdx(jeq) = istor
        p_IsubdiagonalEdges(istor)  = iedge

      end do
      
      ! Compute the maximum number of edges adjacent to one vertex
      rafcstab%NNVEDGE = max(rafcstab%NNVEDGE,&
          p_IsuperdiagonalEdgesIdx(ieq+1) - p_IsuperdiagonalEdgesIdx(ieq) +&
          p_IsubdiagonalEdgesIdx(ieq+1) - p_IsubdiagonalEdgesIdx(ieq))
    end do
    
    ! Reshuffle pass 2: Adjust the index vector which was tainted in
    ! the above loop
    do ieq = neq+1, 2, -1
      p_IsubdiagonalEdgesIdx(ieq) = p_IsubdiagonaledgesIdx(ieq-1)+1
    end do
    p_IsubdiagonalEdgesIdx(1) = 1

    ! Set specifier for extended edge structure
    rafcstab%iSpec = ior(rafcstab%iSpec,AFCSTAB_SUBDIAGONALEDGES)
  end subroutine afcstab_generateSubdiagEdges
  
  !*****************************************************************************

!<subroutine>

  subroutine afcstab_generateExtSparsity(rmatrixSrc, rmatrixExtended)

!<description>
    ! This subroutine generates the extended sparsity pattern
    ! required to assemble the Jacobian matrix from a sparse finite
    ! element matrix stored in CRS format. The basic idea is as
    ! follows. Let A ={a_ij} denote the adjacency matrix which
    ! represents the undirected connectivity graph of the finite
    ! element matrix (rmatrix). The matrix coefficients are
    ! henceforth given by
    !   a_ij=1 if there exists some edge ij, a_ij=0 otherwise
    ! Now, let Z=A^2. Then z_ij>0 if and only if there exists a path
    ! of length two connecting nodes i and j. This is due
    !   z_ij=sum_k(a_ik*a_kj)>0 <=> ex. k : a_ik=1 and a_kj=1.
!</description>

!<input>
    type(t_matrixScalar), intent(IN)    :: rmatrixSrc
!</input>

!<inputoutput>
    type(t_matrixScalar), intent(INOUT) :: rmatrixExtended
!</inputoutput>
!</subroutine>
    
    ! Clear output matrix
    if (lsyssc_hasMatrixStructure(rmatrixExtended) .or.&
        lsyssc_hasMatrixContent(rmatrixExtended)) then
      call lsyssc_releaseMatrix(rmatrixExtended)
    end if

    ! Compute Z=A*A and let the connectivity graph of Z be the
    ! extended sparsity pattern of the Jacobian matrix
    call lsyssc_multMatMat(rmatrixSrc, rmatrixSrc, rmatrixExtended,&
                           .true., .true., .false.)

  end subroutine afcstab_generateExtSparsity

  !*****************************************************************************

!<function>
  
  elemental function afcstab_limit_unbounded(p, q, default) result(r)

!<description>
    ! This function computes the ratio Q/P. If the denominator is
    ! too small, then the default value is applied.
!</description>

!<input>
    ! (de)nominator
    real(DP), intent(IN) :: p,q

    ! default value
    real(DP), intent(IN) :: default
!</input>

!<result>
    ! limited ratio
    real(DP) :: r
!</result>
!</function>

    if (p > SYS_EPSREAL) then
      r = q/p
    else
      r = default
    end if
  end function afcstab_limit_unbounded

  !*****************************************************************************

!<function>
  
  elemental function afcstab_limit_bounded(p, q, default, dbound) result(r)

!<description>
    ! This function computes the limited ratio Q/P and bounds the
    ! result by the size of dbound. If the denominator is too small
    ! then the default value is applied.
!</description>

!<input>
    ! (de)nominator
    real(DP), intent(IN) :: p,q
    
    ! default value
    real(DP), intent(IN) :: default

    ! upper bound
    real(DP), intent(IN) :: dbound
!</input>

!<result>
    ! limited ratio
    real(DP) :: r
!</result>
!</function>
    
    if (p > SYS_EPSREAL) then
      r = min(q/p, dbound)
    else
      r = default
    end if
  end function afcstab_limit_bounded
end module afcstabilisation
