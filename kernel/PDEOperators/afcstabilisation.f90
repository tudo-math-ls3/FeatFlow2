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
!#     -> Creates a stabilisation structure and initialize
!#        it from the values of a given parameter list
!#
!# 2.) afcstab_releaseStabilisation
!#     -> Releases a stabilisation structure
!#
!# 3.) afcstab_resizeStabilisation = afcstab_resizeStabDirect /
!#                                   afcstab_resizeStabIndScalar
!#                                   afcstab_resizeStabIndBlock
!#     -> Resizes a stabilisation structure
!#
!# 4.) afcstab_getbase_IverticesAtEdge
!#     -> Returns pointer to the vertices at edge structure
!#
!# 5.) afcstab_getbase_IsupdiagEdgeIdx
!#     -> Returns pointer to the index pointer for the
!#        superdiagonal edge numbers
!#
!# 6.) afcstab_getbase_IsubdiagEdgeIdx
!#     -> Returns pointer to the index pointer for the
!#        subdiagonal edge numbers
!#
!# 7.) afcstab_getbase_IsubdiagEdge
!#     -> Returns pointer to the subdiagonal edge numbers
!#
!# 8.) afcstab_getbase_DcoeffsAtEdge
!#     -> Returns pointer to edge data
!#
!# 9.) afcstab_generateVerticesAtEdge
!#     -> Generates the standard edge data structure
!#
!# 10.) afcstab_generateOffdiagEdges
!#      -> Generates the subdiagonal edge data structure
!#
!# 11.) afcstab_generateExtSparsity
!#      -> Generates the extended sparsity pattern
!#
!# The following auxiliary routines are available:
!#
!# 1.) afcstab_limit = afcstab_limit_unbounded/
!#                     afcstab_limit_bounded
!#     -> Compute the nodal correction factors, i.e., the ratio of
!#        admissible solution increments and raw antidiffusion
!# </purpose>
!##############################################################################
module afcstabilisation

  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use storage
  use triangulation

  implicit none
  
  private
  public :: t_afcstab
  public :: afcstab_initFromParameterlist
  public :: afcstab_releaseStabilisation
  public :: afcstab_resizeStabilisation
  public :: afcstab_getbase_IverticesAtEdge
  public :: afcstab_getbase_IsupdiagEdgeIdx
  public :: afcstab_getbase_IsubdiagEdgeIdx
  public :: afcstab_getbase_IsubdiagEdge
  public :: afcstab_getbase_DcoeffsAtEdge
  public :: afcstab_generateVerticesAtEdge
  public :: afcstab_generateOffdiagEdges
  public :: afcstab_generateExtSparsity
  public :: afcstab_limit
 
  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>
!<constantblock description="Global format flags for AFC stabilisation">

  ! No stabilisation: use standard high-order Galerkin discretisation
  integer, parameter, public :: AFCSTAB_GALERKIN              = 0
  
  ! Stabilisation of discrete upwind type for convection operators
  integer, parameter, public :: AFCSTAB_UPWIND                = 1

  ! Stabilisation of discrete maximum principle preserving 
  ! type for anisotropic diffusion operators
  integer, parameter, public :: AFCSTAB_DMP                   = 2

  ! Stabilisation of semi-implicit FEM-FCT type
  integer, parameter, public :: AFCSTAB_FEMFCT_IMPLICIT       = 10

  ! Stabilisation of semi-explicit (classical) FEM-FCT type
  integer, parameter, public :: AFCSTAB_FEMFCT_CLASSICAL      = 11

  ! Stabilisation of linearised FEM-FCT type
  integer, parameter, public :: AFCSTAB_FEMFCT_LINEARISED     = 12
  
  ! Stabilisation of iterative FEM-FCT type
  integer, parameter, public :: AFCSTAB_FEMFCT_ITERATIVE      = 13
  
  ! Stabilisation of characteristic FEM-FCT type
  integer, parameter, public :: AFCSTAB_FEMFCT_CHARACTERISTIC = 14

  ! Stabilisation of FEM-FCT type for mass antidiffusion
  integer, parameter, public :: AFCSTAB_FEMFCT_MASS           = 15

  ! Stabilisation of FEM-TVD type for convection operators
  integer, parameter, public :: AFCSTAB_FEMTVD                = 20

  ! Stabilisation of general purpose type for convection operators
  integer, parameter, public :: AFCSTAB_FEMGP                 = 21
  
  ! Stabilisation of symmetric type for diffusion operators
  integer, parameter, public :: AFCSTAB_SYMMETRIC             = 30
  
!</constantblock>

!<constantblock description="Bitfield identifiers for state of stabilisation">
  
  ! Stabilisation is undefined
  integer(I32), parameter, public :: AFCSTAB_UNDEFINED            = 2_I32**0

  ! Stabilisation has been initialised
  integer(I32), parameter, public :: AFCSTAB_INITIALISED          = 2_I32**1

  ! Edge-based structure generated: IverticesAtEdge
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGESTRUCTURE    = 2_I32**2

  ! Edge-based structure oriented: IverticesAtEdge
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGEORIENTATION  = 2_I32**3

  ! Edge-based values computed from matrix: DcoefficientsAtEdge
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGEVALUES       = 2_I32**4

  ! Subdiagonal edge-based structure generated
  integer(I32), parameter, public :: AFCSTAB_HAS_OFFDIAGONALEDGES = 2_I32**5
  
  ! Precomputed antidiffusive fluxes
  integer(I32), parameter, public :: AFCSTAB_HAS_ADFLUXES         = 2_I32**6
  
  ! Nodal sums of antidiffusive increments: PP, PM
  integer(I32), parameter, public :: AFCSTAB_HAS_ADINCREMENTS     = 2_I32**7

  ! Nodal upper/lower bounds: QP, QM
  integer(I32), parameter, public :: AFCSTAB_HAS_BOUNDS           = 2_I32**8
  
  ! Nodal correction factors: RP, RM
  integer(I32), parameter, public :: AFCSTAB_HAS_NODELIMITER      = 2_I32**9
  
  ! Edgewise correction factors: ALPHA
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGELIMITER      = 2_I32**10
    
!</constantblock>

!<constantblock description="Bitfield identifiers for FCT-algorithm">

  ! Initialize the edgewise correction factors by unity
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_INITALPHA    = 2_I32**0
  
  ! Compute the sums of antidiffusive increments
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_ADINCREMENTS = 2_I32**1

  ! Compute the distances to a local extremum
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_BOUNDS       = 2_I32**2

  ! Compute the nodal correction factors
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_LIMITNODAL   = 2_I32**3

  ! Compute edgewise correction factors
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_LIMITEDGE    = 2_I32**4
  
  ! Correct raw antidiffusive fluxes and apply them
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_CORRECT      = 2_I32**5

  ! Scale corrected antidiffusive fluxes by the inverse of the
  ! lumped mass matrix prior to applying it to the residual/solution
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_SCALEBYMASS  = 2_I32**6

  ! Constrain the raw antidiffusive fluxes by the size of limited
  ! limited rather than by the correction factors directly
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_CONSTRAIN    = 2_I32**7
  
  ! FEM-FCT algorithm without application of the corrected fluxes
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_PREPARE  = AFCSTAB_FCTALGO_INITALPHA +&
                                                                AFCSTAB_FCTALGO_ADINCREMENTS +&
                                                                AFCSTAB_FCTALGO_BOUNDS +&
                                                                AFCSTAB_FCTALGO_LIMITNODAL +&
                                                                AFCSTAB_FCTALGO_LIMITEDGE

  ! Standard FEM-FCT algorithm
  ! (also used in the very first iteration of all other variants)
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_STANDARD = AFCSTAB_FCTALGO_PREPARE +&
                                                                AFCSTAB_FCTALGO_CORRECT
  
!</constantblock>

!</constants>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<types>
!<typeblock>

  ! This structures holds all required information for stabilisation
  ! of algebraic flux correction type. It is applicable to scalar
  ! problems and systems of equations alike. Depending on the
  ! stabilisation strategy some internal structures will be generated.
  
  type t_afcstab
    
    ! Format Tag. Identifies the type of stabilisation
    integer :: ctypeAFCstabilisation = AFCSTAB_GALERKIN

    ! Format Tag: Specifies the stabilisation
    integer :: iSpec = AFCSTAB_UNDEFINED

    ! Number of equations of the sparsity pattern
    integer :: NEQ = 0

    ! Number of local variables; in general scalar solution vectors of
    ! size NEQ posses NEQ entries. However, scalar vectors can be
    ! interleaved, that is, each of the NEQ entries stores NVAR local
    ! variables. In this case, NEQ remains unmodified but NVAR>1 such
    ! that the physical length of the vector is NEQ*NVAR.
    integer :: NVAR = 1

    ! Number of local variables after transformation; this variable is
    ! similar to NVAR except for the fact that it corresponds to the
    ! transformed variable. Flux correction for systems of equations
    ! may require some sort of synchronisation which may be performed
    ! in terms of, e.g., primitive variables whereas the system
    ! matrices and the residual vector/right-hand side are assembled
    ! in terms of conservative variables. 
    integer :: NVARtransformed = 1

    ! Number of edges of the sparsity pattern
    integer :: NEDGE = 0

    ! Maximum number of edges adjacent to one vertex. This
    ! corresponds to the maximum number of nonzero row entries.
    integer :: NNVEDGE = 0

    ! Flag: compute auxiliary structures for flux prelimiting
    logical :: bprelimiting = .true.

    ! Handle to vertices at edge structure
    ! IverticesAtEdge(1:2,1:NEDGE) : the two end-points of the edge
    ! IverticesAtEdge(3:4,1:NEDGE) : the two matrix position that
    !                                correspond to the edge
    integer :: h_IverticesAtEdge = ST_NOHANDLE

    ! Handle to index pointer for superdiagonal edge numbers
    ! The numbers IsuperdiagEdgesIdx(i):IsuperdiagEdgesIdx(i+1)-1
    ! denote the edge numbers of the ith vertex which are located in
    ! the upper right triangular matrix.
    integer :: h_IsuperdiagEdgesIdx = ST_NOHANDLE

    ! Handle to index pointer for subdiagonal edge numbers
    integer :: h_IsubdiagEdgesIdx = ST_NOHANDLE

    ! Handle to the subdiagonal edge numbers
    integer :: h_IsubdiagEdges = ST_NOHANDLE

    ! Handle to coefficient at edge structure
    integer :: h_DcoefficientsAtEdge = ST_NOHANDLE

    ! REMARK: The following pointers are linked to some position of the
    ! auxiliary vectors (see below) during the initialisation of the
    ! stabilisation structure. You should never address any of the
    ! auxiliary vectors directly but always use one of the pointers
    ! p_rvectorXYZ below to be sure to obtain the correct data.

    ! Pointer to the vector of correction factors
    type(t_vectorScalar), pointer :: p_rvectorAlpha => null()

    ! Pointer to the vector of explicit antidiffusive fluxes
    type(t_vectorScalar), pointer :: p_rvectorFlux0 => null()

    ! Pointer to the vector of raw antidiffusive fluxes
    type(t_vectorScalar), pointer :: p_rvectorFlux => null()

    ! Pointer to the vector of prelimiting antidiffusive fluxes
    type(t_vectorScalar), pointer :: p_rvectorPrelimit => null()

    ! Pointers to the vectors of antidiffusive contributions
    type(t_vectorScalar), pointer :: p_rvectorPp => null()
    type(t_vectorScalar), pointer :: p_rvectorPm => null()

    ! Pointers to the vectors of local solution bounds
    type(t_vectorScalar), pointer :: p_rvectorQp => null()
    type(t_vectorScalar), pointer :: p_rvectorQm => null()

    ! Pointers to the vectors of nodal correction factors
    type(t_vectorScalar), pointer :: p_rvectorRp => null()
    type(t_vectorScalar), pointer :: p_rvectorRm => null()

    ! Pointer to the low-order predictor
    type(t_vectorBlock), pointer :: p_rvectorPredictor => null()


    ! Here, the auxiliary vectors follow which are allocated during the
    ! initialisation of the stabilisation structure and should neve be
    ! addressed directly by the used.

    ! Auxiliary nodal vectors; used internally
    type(t_vectorScalar), dimension(:), pointer :: RnodalVectors => null()

    ! Auxiliary nodal block vectors; used internally
    type(t_vectorBlock), dimension(:), pointer  :: RnodalBlockVectors => null()

    ! Auxiliary edge vectors; used internally
    type(t_vectorScalar), dimension(:), pointer :: RedgeVectors => null()
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
    type(t_parlist), intent(in)    :: rparlist

    ! Section name of the parameter list
    character(LEN=*), intent(in)   :: ssectionName
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: istabilisation,iprelimiting

    ! Get type of stabilisation from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        "istabilisation", istabilisation)
    
    ! Check if stabilisation should be applied
    if (istabilisation .eq. AFCSTAB_GALERKIN) then
      
      return   ! -> high-order Galerkin
      
    elseif ((istabilisation .ne. AFCSTAB_UPWIND)            .and. &
        (    istabilisation .ne. AFCSTAB_FEMFCT_CLASSICAL)  .and. &
        (    istabilisation .ne. AFCSTAB_FEMFCT_IMPLICIT)   .and. &
        (    istabilisation .ne. AFCSTAB_FEMFCT_LINEARISED) .and. &
        (    istabilisation .ne. AFCSTAB_FEMFCT_ITERATIVE)  .and. &
        (    istabilisation .ne. AFCSTAB_FEMFCT_MASS)       .and. &
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

    ! Get type of prelimiting from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        "iprelimiting", iprelimiting, 1)

    ! Set flag for prelimiting
    rafcstab%bprelimiting = (iprelimiting .ne. 0)

  end subroutine afcstab_initFromParameterlist

  !*****************************************************************************

!<subroutine>
  
  subroutine afcstab_releaseStabilisation(rafcstab)

!<description>
    ! This subroutine releases a stabilisation structure
!</description>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! Free storage
    if (rafcstab%h_IsuperdiagEdgesIdx .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IsuperdiagEdgesIdx)
    if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IverticesAtEdge)
    if (rafcstab%h_IsubdiagEdgesIdx .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IsubdiagEdgesIdx)
    if (rafcstab%h_IsubdiagEdges .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IsubdiagEdges)
    if (rafcstab%h_DcoefficientsAtEdge .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_DcoefficientsAtEdge)
    
    ! Reset atomic data
    rafcstab%ctypeAFCstabilisation = AFCSTAB_GALERKIN
    rafcstab%iSpec                 = AFCSTAB_UNDEFINED
    rafcstab%NEQ                   = 0
    rafcstab%NVAR                  = 1
    rafcstab%NVARtransformed       = 1
    rafcstab%NEDGE                 = 0
    rafcstab%NNVEDGE               = 0

    ! Nullify pointers
    nullify(rafcstab%p_rvectorAlpha)
    nullify(rafcstab%p_rvectorFlux0)
    nullify(rafcstab%p_rvectorFlux)
    nullify(rafcstab%p_rvectorPp)
    nullify(rafcstab%p_rvectorPm)
    nullify(rafcstab%p_rvectorQp)
    nullify(rafcstab%p_rvectorQm)
    nullify(rafcstab%p_rvectorRp)
    nullify(rafcstab%p_rvectorRm)
    nullify(rafcstab%p_rvectorPredictor)

    ! Release auxiliary nodal block vectors
    if (associated(rafcstab%RnodalBlockVectors)) then
      do i = lbound(rafcstab%RnodalBlockVectors,1),&
             ubound(rafcstab%RnodalBlockVectors,1)
        call lsysbl_releaseVector(rafcstab%RnodalBlockVectors(i))
      end do
      deallocate(rafcstab%RnodalBlockVectors)
    end if

    ! Release auxiliary nodal vectors
    if (associated(rafcstab%RnodalVectors)) then
      do i = lbound(rafcstab%RnodalVectors,1),&
             ubound(rafcstab%RnodalVectors,1)
        call lsyssc_releaseVector(rafcstab%RnodalVectors(i))
      end do
      deallocate(rafcstab%RnodalVectors)
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
    integer, intent(in) :: neq

    ! number of edges
    integer, intent(in) :: nedge   
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    
    ! Resize nodal quantities
    if (rafcstab%NEQ .ne. neq) then

      ! Set new number of nodes
      rafcstab%NEQ = neq
      
      ! Resize edge index vector
      if (rafcstab%h_IsuperdiagEdgesIdx .ne. ST_NOHANDLE) then
        call storage_realloc('afcstab_resizeStabDirect',&
            rafcstab%NEQ+1, rafcstab%h_IsuperdiagEdgesIdx,&
            ST_NEWBLOCK_NOINIT, .false.)
      end if
      
      ! Resize subdiagonal edge index vector
      if (rafcstab%h_IsubdiagEdgesIdx .ne. ST_NOHANDLE) then
        call storage_realloc('afcstab_resizeStabDirect',&
            rafcstab%NEQ+1, rafcstab%h_IsubdiagEdgesIdx,&
            ST_NEWBLOCK_NOINIT, .false.)
      end if

      ! Resize auxiliary nodal vectors
      if(associated(rafcstab%RnodalVectors)) then
        do i = lbound(rafcstab%RnodalVectors,1),&
               ubound(rafcstab%RnodalVectors,1)
          if (rafcstab%RnodalVectors(i)%NEQ .ne. 0)&
              call lsyssc_resizeVector(rafcstab%RnodalVectors(i),&
              rafcstab%NEQ, .false., .false.)
        end do
      end if

      ! Resize auxiliary nodal vectors
      if(associated(rafcstab%RnodalBlockVectors)) then
        do i = lbound(rafcstab%RnodalBlockVectors,1),&
               ubound(rafcstab%RnodalBlockVectors,1)
          if (rafcstab%RnodalBlockVectors(i)%NEQ .ne. 0)&
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
      if (rafcstab%h_IsubdiagEdges .ne. ST_NOHANDLE) then
        call storage_realloc('afcstab_resizeStabDirect',&
            rafcstab%NEDGE, rafcstab%h_IsubdiagEdges,&
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
          if (rafcstab%RedgeVectors(i)%NEQ .ne. 0)&
              call lsyssc_resizeVector(rafcstab%RedgeVectors(i),&
              rafcstab%NEDGE, .false., .false.)
        end do
      end if
    end if

    ! Set state of stabilisation
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      rafcstab%iSpec = AFCSTAB_UNDEFINED
    else
      rafcstab%iSpec = AFCSTAB_INITIALISED
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
    type(t_matrixScalar), intent(in) :: rmatrixTemplate
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout)   :: rafcstab
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
    type(t_matrixBlock), intent(in) :: rmatrixBlockTemplate
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: neq, nedge

    
    ! Check if block matrix has only one block
    if ((rmatrixBlockTemplate%nblocksPerCol .eq. 1) .and. &
        (rmatrixBlockTemplate%nblocksPerRow .eq. 1)) then
      call afcstab_resizeStabIndScalar(rafcstab,&
          rmatrixBlockTemplate%RmatrixBlock(1,1))

      ! That is it
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

  subroutine afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagEdgesIdx)

!<description>
    ! Returns a pointer to the index pointer for vertices at edge structure
!</description>

!<input>
    ! discrete operator
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the index pointer for superdiagonal edge numbers.
    ! NULL() if the discrete operator does not provide it.
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    if ((rafcstab%h_IsuperdiagEdgesIdx .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEQ .eq. 0)) then
      nullify(p_IsuperdiagEdgesIdx)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rafcstab%h_IsuperdiagEdgesIdx,&
        p_IsuperdiagEdgesIdx,rafcstab%NEQ+1)

  end subroutine afcstab_getbase_IsupdiagEdgeIdx

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

!<description>
    ! Returns a pointer to the vertices at edge structure
!</description>

!<input>
    ! discrete operator
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the vertices at edge structure
    ! NULL() if the discrete operator does not provide it.
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
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

  subroutine afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagEdgesIdx)

!<description>
    ! Returns a pointer to the index pointer for 
    ! subdiagonal edge numbers
!</description>

!<input>
    ! discrete operator
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the index pointer for 
    ! subdiagonal edge numbers
    ! NULL() if the discrete operator does not provide it.
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    if ((rafcstab%h_IsubdiagEdgesIdx .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEQ                    .eq. 0)) then
      nullify(p_IsubdiagEdgesIdx)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rafcstab%h_IsubdiagEdgesIdx,&
        p_IsubdiagEdgesIdx,rafcstab%NEQ+1)

  end subroutine afcstab_getbase_IsubdiagEdgeIdx

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IsubdiagEdge(rafcstab,p_IsubdiagEdges)

!<description>
    ! Returns a pointer to the subdiagonal edge number
!</description>

!<input>
    ! discrete operator
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the subdiagonal edge numbers
    ! NULL() if the discrete operator does not provide it.
    integer, dimension(:), pointer :: p_IsubdiagEdges
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    if ((rafcstab%h_IsubdiagEdges .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE               .eq. 0)) then
      nullify(p_IsubdiagEdges)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rafcstab%h_IsubdiagEdges,&
        p_IsubdiagEdges,rafcstab%NEDGE)

  end subroutine afcstab_getbase_IsubdiagEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_DcoeffsAtEdge(rafcstab,p_DcoefficientsAtEdge)

!<description>
    ! Returns a pointer to the double-valued edge data
!</description>

!<input>
    ! discrete operator
    type(t_afcstab), intent(in) :: rafcstab
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

  subroutine afcstab_generateVerticesAtEdge(rmatrixTemplate, rafcstab)

!<description>
    ! This subroutine generates the list of edges which are
    ! characterized by their two endpoints (i,j) and the absolute
    ! position of matrix entries ij and ji. 
!</description>

!<input>
    ! template matrix
    type(t_matrixScalar), intent(in) :: rmatrixTemplate
!</input>

!<inputoutput>
    ! discrete operator
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal, p_Ksep
    integer :: h_Ksep


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%iSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_generateVerticesAtEdge')
      call sys_halt()
    end if
    
    ! Set pointer to edge structure
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

    ! What kind of matrix are we?
    select case(rmatrixTemplate%cmatrixFormat)
    case(LSYSSC_MATRIX7,&
         LSYSSC_MATRIX7INTL)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrixTemplate, p_Kld)
      call lsyssc_getbase_Kcol(rmatrixTemplate, p_Kcol)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrixTemplate%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrixTemplate%NEQ+1)

      ! Generate edge structure
      call genEdgesMat7(p_Kld, p_Kcol, p_Ksep, p_IverticesAtEdge)

      ! Release diagonal separator
      call storage_free(h_Ksep)
      
      
    case(LSYSSC_MATRIX9,&
         LSYSSC_MATRIX9INTL)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrixTemplate, p_Kld)
      call lsyssc_getbase_Kcol(rmatrixTemplate, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrixTemplate, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrixTemplate%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrixTemplate%NEQ+1)

      ! Generate edge structure
      call genEdgesMat9(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_IverticesAtEdge)

      ! Release diagonal separator
      call storage_free(h_Ksep)


    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_generateVerticesAtEdge')
      call sys_halt()
    end select

    ! Set state of stabiliation
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_EDGESTRUCTURE)
         
  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Generate edge data structure for matrices in format 7

    subroutine genEdgesMat7(Kld, Kcol, Ksep, IverticesAtEdge)

      integer, dimension(:), intent(in) :: Kld, Kcol
      integer, dimension(:), intent(inout) :: Ksep
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      ! local variables
      integer :: i,j,ij,ji,iedge
      
      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, size(Kld)-1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
       
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Increatse edge counter
          iedge = iedge+1

          ! Set node numbers i and j
          IverticesAtEdge(1, iedge) = i
          IverticesAtEdge(2, iedge) = j
          
          ! Set matrix positions ij and ji
          IverticesAtEdge(3, iedge) = ij
          IverticesAtEdge(4, iedge) = ji
          
        end do
      end do

    end subroutine genEdgesMat7

    !**************************************************************
    ! Generate edge data structure for matrices in format 9

    subroutine genEdgesMat9(Kld, Kcol, Kdiagonal, Ksep, IverticesAtEdge)

      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal
      integer, dimension(:), intent(inout) :: Ksep
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      ! local variables
      integer :: i,j,ij,ji,iedge
      
      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, size(Kld)-1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
       
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Increatse edge counter
          iedge = iedge+1

          ! Set node numbers i and j
          IverticesAtEdge(1, iedge) = i
          IverticesAtEdge(2, iedge) = j
          
          ! Set matrix positions ij and ji
          IverticesAtEdge(3, iedge) = ij
          IverticesAtEdge(4, iedge) = ji
          
        end do
      end do

    end subroutine genEdgesMat9
    
  end subroutine afcstab_generateVerticesAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_generateOffdiagEdges(rafcstab)

!<description>
    ! This subroutine generates the edge data structure
    ! (superdiagonal separator and subdiagonal edges)
    ! based on a given edge data structure.
!</description>

!<inputoutput>
    ! discrete operator
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer   :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer   :: p_IsubdiagEdges
    integer, dimension(:), pointer   :: p_IsubdiagEdgesIdx
    integer :: iedge,nedge,istor
    integer :: ieq,jeq,neq
    integer :: isize

    ! Check if edge-based data structure is prepared
    if (iand(rafcstab%iSpec,AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) then
      call output_line('Discrete operator does not provide required ' // &
          'edge-based data structure',OU_CLASS_ERROR,OU_MODE_STD,&
          'afcstab_generateOffdiagEdge')
      call sys_halt()
    end if

    ! store dimensions of discrete operator
    neq   = rafcstab%NEQ
    nedge = rafcstab%NEDGE

    ! Allocate memory (if required)
    if (rafcstab%h_IsuperdiagEdgesIdx .eq. ST_NOHANDLE) then
      call storage_new('afcstab_generateOffdiagEdges','IsuperdiagEdgesIdx',&
          neq+1,ST_INT,rafcstab%h_IsuperdiagEdgesIdx,ST_NEWBLOCK_ZERO)
    else
      call storage_getsize(rafcstab%h_IsuperdiagEdgesIdx,isize)
      if (isize < neq+1) then
        call storage_realloc('afcstab_generateOffdiagEdges',neq+1,&
            rafcstab%h_IsuperdiagEdgesIdx,ST_NEWBLOCK_ZERO,.false.)
      else
        call storage_clear(rafcstab%h_IsuperdiagEdgesIdx)
      end if
    end if

    ! Allocate memory (if required)
    if (rafcstab%h_IsubdiagEdgesIdx .eq. ST_NOHANDLE) then
      call storage_new('afcstab_generateOffdiagEdges','IsubdiagEdgesIdx',&
          neq+1,ST_INT,rafcstab%h_IsubdiagEdgesIdx,ST_NEWBLOCK_ZERO)
    else
      call storage_getsize(rafcstab%h_IsubdiagEdgesIdx,isize)
      if (isize < neq+1) then
        call storage_realloc('afcstab_generateOffdiagEdges',neq+1,&
            rafcstab%h_IsubdiagEdgesIdx,ST_NEWBLOCK_ZERO,.false.)
      else
        call storage_clear(rafcstab%h_IsubdiagEdgesIdx)
      end if
    end if

    ! Allocate memory (if required)
    if (rafcstab%h_IsubdiagEdges .eq. ST_NOHANDLE) then
      call storage_new('afcstab_generateOffdiagEdges','IsubdiagEdges',&
          nedge,ST_INT,rafcstab%h_IsubdiagEdges,ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rafcstab%h_IsubdiagEdges,isize)
      if (isize < nedge) then
        call storage_realloc('afcstab_generateOffdiagEdges',nedge,&
            rafcstab%h_IsubdiagEdges,ST_NEWBLOCK_NOINIT,.false.)
      end if
    end if
    
    ! Set pointers
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab,p_IsuperdiagEdgesIdx)
    call afcstab_getbase_IverticesAtEdge(rafcstab,p_IverticesAtEdge)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab,p_IsubdiagEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab,p_IsubdiagEdges)
    
    ! Count the number of edges in each row
    do iedge = 1, nedge
      
       ! Determine the start-point of the current edge
      ieq = min(p_IverticesAtEdge(1,iedge),&
                p_IverticesAtEdge(2,iedge))

      ! Determine the end-point of the current edge
      jeq = max(p_IverticesAtEdge(1,iedge),&
                p_IverticesAtEdge(2,iedge))

      ! Increase number of edges connected to these points by one
      p_IsuperdiagEdgesIdx(ieq+1) = p_IsuperdiagEdgesIdx(ieq+1)+1
      p_IsubdiagEdgesIdx(jeq+1)   = p_IsubdiagEdgesIdx(jeq+1)+1
    end do
    
    ! Reshuffle pass 1
    rafcstab%NNVEDGE = 0
    p_IsuperdiagEdgesIdx(1) = 1
    p_IsubdiagEdgesIdx(1) = 1
    
    do ieq = 2, neq+1
      ! Compute the maximum number of edges adjacent to one vertex.
      rafcstab%NNVEDGE = max(rafcstab%NNVEDGE,&
          p_IsuperdiagEdgesIdx(ieq)+p_IsubdiagEdgesIdx(ieq))
      
      ! Compute absolute starting positions of edge numbers connected
      ! to the current node IEQ
      p_IsuperdiagEdgesIdx(ieq) =&
          p_IsuperdiagEdgesIdx(ieq) + p_IsuperdiagEdgesIdx(ieq-1)
      p_IsubdiagEdgesIdx(ieq) =&
          p_IsubdiagEdgesIdx(ieq) + p_IsubdiagEdgesIdx(ieq-1)
    end do
    
    ! Store the subdiagonal edge numbers
    do ieq = 1, neq
      ! Loop over the edges located in the upper right triangular matrix
      do iedge = p_IsuperdiagEdgesIdx(ieq),&
                 p_IsuperdiagEdgesIdx(ieq+1)-1
        
        ! Determine the end-point of the edge, i.e. 
        ! the node which is not equal to IEQ
        jeq = p_IverticesAtEdge(1,iedge) + p_IverticesAtEdge(2,iedge) - ieq
        
        ! Get and update next free position
        istor = p_IsubdiagEdgesIdx(jeq)
        p_IsubdiagEdgesIdx(jeq) = istor+1
        
        ! Store global edge number
        p_IsubdiagEdges(istor)  = iedge
      end do
    end do
    
    ! Reshuffle pass 2:
    ! Adjust the index vector which was tainted in the above loop
    do ieq = neq+1, 2, -1
      p_IsubdiagEdgesIdx(ieq) = p_IsubdiagEdgesIdx(ieq-1)
    end do
    p_IsubdiagEdgesIdx(1) = 1

    ! Set specifier for extended edge structure
    rafcstab%iSpec = ior(rafcstab%iSpec, AFCSTAB_HAS_OFFDIAGONALEDGES)

  end subroutine afcstab_generateOffdiagEdges
  
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
    ! Now, let <tex>$Z=A^2$</tex>. Then z_ij>0 if and only if there
    ! exists a path of length two connecting nodes i and j. This is due
    !   z_ij=sum_k(a_ik*a_kj)>0 <=> ex. k : a_ik=1 and a_kj=1.
!</description>

!<input>
    type(t_matrixScalar), intent(in)    :: rmatrixSrc
!</input>

!<inputoutput>
    type(t_matrixScalar), intent(inout) :: rmatrixExtended
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
    real(DP), intent(in) :: p,q

    ! default value
    real(DP), intent(in) :: default
!</input>

!<result>
    ! limited ratio
    real(DP) :: r
!</result>
!</function>

    if (abs(p) > 1e-12_DP) then
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
    real(DP), intent(in) :: p,q
    
    ! default value
    real(DP), intent(in) :: default

    ! upper bound
    real(DP), intent(in) :: dbound
!</input>

!<result>
    ! limited ratio
    real(DP) :: r
!</result>
!</function>
    
    if (abs(p) > 1e-12_DP) then
      r = min(q/p, dbound)
    else
      r = default
    end if
  end function afcstab_limit_bounded

end module afcstabilisation
