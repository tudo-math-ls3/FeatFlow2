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
!# 3.) afcstab_resizeStabilisation
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
MODULE afcstabilisation

  USE fsystem
  USE genoutput
  USE linearsystemscalar
  USE paramlist
  USE storage
  USE triangulation

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_afcstab
  PUBLIC :: afcstab_initFromParameterlist
  PUBLIC :: afcstab_releaseStabilisation
  PUBLIC :: afcstab_resizeStabilisation
  PUBLIC :: afcstab_getbase_IsupdiagEdgeIdx
  PUBLIC :: afcstab_getbase_IverticesAtEdge
  PUBLIC :: afcstab_getbase_DcoeffsAtEdge
  PUBLIC :: afcstab_getbase_IsubdiagEdgeIdx
  PUBLIC :: afcstab_getbase_IsubdiagEdge
  PUBLIC :: afcstab_generateSubdiagEdges
  PUBLIC :: afcstab_generateExtSparsity
  PUBLIC :: afcstab_limit
 
  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>
!<constantblock description="Global format flags for AFC stabilisation">

  ! No stabilisation: use standard high-order Galerkin discretisation
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_GALERKIN        = 0
  
  ! Stabilisation of discrete upwind type for convection operators
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_UPWIND          = 1

  ! Stabilisation of discrete maximum principle preserving 
  ! type for anisotropic diffusion operators
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_DMP             = 2

  ! Stabilisation of semi-implicit FEM-FCT type for convection operators
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_FEMFCT          = 10

  ! Stabilisation of semi-explicit (classical) FEM-FCT type for convection operators
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_FEMFCT_EXP      = 11

  ! Stabilisation of linearised FEM-FCT type for convection operators
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_FEMFCT_LIN      = 12
  
  ! Stabilisation of FEM-TVD type for convection operators
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_FEMTVD          = 20

  ! Stabilisation of general purpose type for convection operators
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_FEMGP           = 21
  
  ! Stabilisation of symmetric type for diffusion operators
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_SYMMETRIC       = 30
  
!</constantblock>

!<constantblock description="Global format flags for dissipation">

  ! Employ scalar dissipation (default)
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_SCALARDISSIPATION = 1

  ! Employ tensorial dissipation
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_TENSORDISSIPATION = 2

!</constantblock>

!<constantblock description="Bitfield identifiers for state of stabilisation">
  
  ! Stabilisation is undefined
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_UNDEFINED       = 2**0

  ! Stabilisation has been initialised
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_INITIALISED     = 2**1

  ! Edge-based structure generated: KEDGE
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_EDGESTRUCTURE   = 2**2

  ! Edge-based structure oriented: KEDGE
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_EDGEORIENTATION = 2**3

  ! Edge-based values computed from matrix: DEDGE
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_EDGEVALUES      = 2**4

  ! Nodal antidiffusion: PP,PM
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_ANTIDIFFUSION   = 2**5

  ! Nodal upper/lower bounds: QP,QM
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_BOUNDS          = 2**6
  
  ! Nodal correction factors computed: RP,RM
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_LIMITER         = 2**7

  ! Antidiffusive fluxes precomputed
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_FLUXES          = 2**8
  
  ! Subdiagonal edge-based structure generated
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_SUBDIAGONALEDGES= 2**9
!</constantblock>

!<constantblock description="Global type of mass matrix treatment">

  ! Adopt the lumped-mass discretisation
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_LUMPEDMASS      = 0

  ! Adopt the consistent-mass discretisation
  INTEGER, PARAMETER, PUBLIC :: AFCSTAB_CONSISTENTMASS  = 1
!</constantblock>
!</constants>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<types>
!<typeblock>

  ! data structure that holds all required information for stabilisation
  TYPE t_afcstab
    
    ! Format Tag. Identifies the type of stabilisation
    INTEGER :: ctypeAFCstabilisation                   = AFCSTAB_GALERKIN

    ! Format Tag. Specifies the type of dissipation
    INTEGER :: iDissipation                            = AFCSTAB_SCALARDISSIPATION

    ! Format Tag: Specifies the stabilisation
    INTEGER :: iSpec                                   = AFCSTAB_UNDEFINED

    ! Format Tag: Specifies whether an extended stencil should be
    ! created for the Jacobian matrix if required
    !   iextendedJacobian = 0 : use sparsity pattern of FE-matrix
    !   iextendedJacobian = 1 : extend sparsity pattern accordingly
    INTEGER :: iextendedJacobian                       = 0

    ! Format Tag: Specifies whether the consistent mass matrix should
    ! be considered in the stabilisation procedure
    !   imass = 0 : neglect consistent mass matrix
    !   imass = 1 : consider consistent mass matrix
    INTEGER :: imass                                   = 0

    ! Number of equations of the sparsity pattern
    INTEGER(PREC_VECIDX) :: NEQ                        = 0

    ! Number of local variables; in general scalar solution vectors of
    ! size NEQ posses NEQ entries. However, scalar vectors can be interleaved,
    ! that is, each of the NEQ entries stores NVAR local variables. In this case,
    ! NEQ remains unmodified but NVAR>1 such that the physical length of the
    ! vector is NEQ*NVAR.
    INTEGER :: NVAR                                    = 1

    ! Number of edges of the sparsity pattern
    INTEGER(PREC_VECIDX) :: NEDGE                      = 0

    ! Maximum number of edges adjacent to one vertex. 
    ! This corresponds to the maximum number of nonzero row entries.
    INTEGER(PREC_VECIDX) :: NNVEDGE                    = 0

    ! Handle to index pointer for superdiagonal edge numbers
    ! INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: IsuperdiagonalEdgesIdx
    ! The numbers IsuperdiagonalEdgesIdx(i):IsuperdiagonalEdgesIdx(i+1)-1
    ! denote the edge numbers of the ith vertex which are located in
    ! the upper right triangular matrix.
    INTEGER :: h_IsuperdiagonalEdgesIdx                    = ST_NOHANDLE

    ! Handle to vertices at edge structure
    ! INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: IverticesAtEdge
    ! IverticesAtEdge(1:2,1:NEDGE) : the two end-points of the edge
    ! IverticesAtEdge(3:4,1:NEDGE) : the two matrix position that
    !                                correspond to the edge
    INTEGER :: h_IverticesAtEdge                       = ST_NOHANDLE

    ! Handle to index pointer for subdiagonal edge numbers
    ! INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: IsubdiagonalEdgesIdx
    INTEGER :: h_IsubdiagonalEdgesIdx                  = ST_NOHANDLE

    ! Handle to the subdiagonal edge numbers
    ! INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: IsubdiagonalEdges
    INTEGER :: h_IsubdiagonalEdges                     = ST_NOHANDLE

    ! Handle to coefficient at edge structure
    INTEGER :: h_DcoefficientsAtEdge                   = ST_NOHANDLE

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
    INTEGER :: isortStrategy                           = 0

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
    INTEGER :: h_IsortPermutation                      = ST_NOHANDLE

    ! Auxiliary nodal vectors; used internally
    TYPE(t_vectorScalar), DIMENSION(:), POINTER :: RnodalVectors => NULL()

    ! Auxiliary edge vectors; used internally
    TYPE(t_vectorScalar), DIMENSION(:), POINTER :: RedgeVectors  => NULL()
  END TYPE t_afcstab
!</typeblock>
!</types>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  INTERFACE afcstab_limit
    MODULE PROCEDURE afcstab_limit_unbounded
    MODULE PROCEDURE afcstab_limit_bounded
  END INTERFACE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE afcstab_initFromParameterlist(rparlist, ssectionName, rafcstab)

!<description>
    ! This subroutine creates a stabilisation structure and initializes
    ! its values from a given parameter list
!</description>

!<input>
    ! parameter list
    TYPE(t_parlist), INTENT(IN)    :: rparlist

    ! Section name of the parameter list
    CHARACTER(LEN=*), INTENT(IN)   :: ssectionName
!</input>

!<inputoutput>
    ! Stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variable
    INTEGER :: istabilisation

    ! Get type of stabilisation from parameter list
    CALL parlst_getvalue_int(rparlist, ssectionName,&
        "istabilisation", istabilisation)
    
    ! Check if stabilisation should be applied
    IF (istabilisation .EQ. AFCSTAB_GALERKIN) THEN
      
      RETURN   ! -> high-order Galerkin
      
    ELSEIF ((istabilisation .NE. AFCSTAB_UPWIND)    .AND. &
        (    istabilisation .NE. AFCSTAB_FEMFCT)    .AND. &
        (    istabilisation .NE. AFCSTAB_FEMFCT_EXP).AND. &
        (    istabilisation .NE. AFCSTAB_FEMFCT_LIN).AND. &
        (    istabilisation .NE. AFCSTAB_FEMTVD)    .AND. &
        (    istabilisation .NE. AFCSTAB_FEMGP)     .AND. &
        (    istabilisation .NE. AFCSTAB_DMP)       .AND. &
        (    istabilisation .NE. AFCSTAB_SYMMETRIC)) THEN 
      
      CALL output_line('Invalid AFC type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_initFromParameterlist')
      CALL sys_halt()

    ELSE
      ! Set type of stabilisation
      rafcstab%iSpec                 = AFCSTAB_UNDEFINED
      rafcstab%ctypeAFCstabilisation = istabilisation
      
      ! Set additional parameters
      CALL parlst_getvalue_int(rparlist, ssectionName,&
          "iextendedJacobian", rafcstab%iextendedJacobian, 0)
      CALL parlst_getvalue_int(rparlist, ssectionName,&
          "idissipation", rafcstab%idissipation, AFCSTAB_SCALARDISSIPATION)
      CALL parlst_getvalue_int(rparlist, ssectionName,&
          "imass", rafcstab%imass, AFCSTAB_LUMPEDMASS)
    END IF
  END SUBROUTINE afcstab_initFromParameterlist

  !*****************************************************************************

!<subroutine>
  
  SUBROUTINE afcstab_releaseStabilisation(rafcstab)

!<description>
    ! This subroutine releases a stabilisation structure
!</description>

!<inputoutput>
    ! Stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER :: i

    ! Free storage
    IF (rafcstab%h_IsuperdiagonalEdgesIdx .NE. ST_NOHANDLE)&
        CALL storage_free(rafcstab%h_IsuperdiagonalEdgesIdx)
    IF (rafcstab%h_IverticesAtEdge .NE. ST_NOHANDLE)&
        CALL storage_free(rafcstab%h_IverticesAtEdge)
    IF (rafcstab%h_IsubdiagonalEdgesIdx .NE. ST_NOHANDLE)&
        CALL storage_free(rafcstab%h_IsubdiagonalEdgesIdx)
    IF (rafcstab%h_IsubdiagonalEdges .NE. ST_NOHANDLE)&
        CALL storage_free(rafcstab%h_IsubdiagonalEdges)
    IF (rafcstab%h_DcoefficientsAtEdge .NE. ST_NOHANDLE)&
        CALL storage_free(rafcstab%h_DcoefficientsAtEdge)
    
    ! Reset atomic data
    rafcstab%ctypeAFCstabilisation = AFCSTAB_GALERKIN
    rafcstab%iDissipation          = AFCSTAB_SCALARDISSIPATION
    rafcstab%iSpec                 = AFCSTAB_UNDEFINED
    rafcstab%imass                 = 0
    rafcstab%NEQ                   = 0
    rafcstab%NVAR                  = 1
    rafcstab%NEDGE                 = 0
    rafcstab%NNVEDGE               = 0
    rafcstab%isortStrategy         = 0
    rafcstab%h_IsortPermutation    = ST_NOHANDLE

    ! Release auxiliary nodal vectors
    IF (ASSOCIATED(rafcstab%RnodalVectors)) THEN
      DO i = LBOUND(rafcstab%RnodalVectors,1),&
             UBOUND(rafcstab%RnodalVectors,1)
        CALL lsyssc_releaseVector(rafcstab%RnodalVectors(i))
      END DO
      DEALLOCATE(rafcstab%RnodalVectors)
    END IF

    ! Release auxiliary edge vectors
    IF (ASSOCIATED(rafcstab%RedgeVectors)) THEN
      DO i = LBOUND(rafcstab%RedgeVectors,1),&
             UBOUND(rafcstab%RedgeVectors,1)
        CALL lsyssc_releaseVector(rafcstab%RedgeVectors(i))
      END DO
      DEALLOCATE(rafcstab%RedgeVectors)
    END IF
  END SUBROUTINE afcstab_releaseStabilisation

   !*****************************************************************************

!<subroutine>

  SUBROUTINE afcstab_resizeStabilisation(rafcstab, neq, nedge, nvar)

!<description>
    ! This subroutine resizes all vectors of the stabilisation structure
    ! to the new values NEQ, NEDGE and eventually NVAR.
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq

    ! number of edges
    INTEGER(PREC_VECIDX), INTENT(IN) :: nedge

    ! OPTIONAL: number of local variables
    INTEGER, INTENT(IN), OPTIONAL    :: nvar
!</input>

!<inputoutput>
    ! stabilisation structure
    TYPE(t_afcstab), INTENT(INOUT)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER :: i

    ! Check if dimension NVAR is correct
    IF (PRESENT(nvar)) THEN
      IF (rafcstab%NVAR .NE. nvar) THEN
        CALL output_line('Invalid number of variables!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabilisation')
        CALL sys_halt()
      END IF
    END IF

    ! Resize nodal quantities
    IF (rafcstab%NEQ .NE. neq) THEN

      ! Set new number of nodes
      rafcstab%NEQ = neq
      
      ! Resize edge index vector
      IF (rafcstab%h_IsuperdiagonalEdgesIdx .NE. ST_NOHANDLE) THEN
        CALL storage_realloc('afcstab_resizeStabilisation',&
            rafcstab%NEQ+1, rafcstab%h_IsuperdiagonalEdgesIdx,&
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
      
      ! Resize subdiagonal edge index vector
      IF (rafcstab%h_IsubdiagonalEdgesIdx .NE. ST_NOHANDLE) THEN
        CALL storage_realloc('afcstab_resizeStabilisation',&
            rafcstab%NEQ+1, rafcstab%h_IsubdiagonalEdgesIdx,&
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF

      ! Resize auxiliary nodal vectors
      IF(ASSOCIATED(rafcstab%RnodalVectors)) THEN
        DO i = LBOUND(rafcstab%RnodalVectors,1),&
               UBOUND(rafcstab%RnodalVectors,1)
          CALL lsyssc_resizeVector(rafcstab%RnodalVectors(i),&
              rafcstab%NEQ, .FALSE.)
        END DO
      END IF
    END IF


    ! Resize edge quantities
    IF (rafcstab%NEDGE .NE. nedge) THEN

      ! Set new number of edges
      rafcstab%NEDGE = nedge

      ! Resize array of edges
      IF (rafcstab%h_IverticesAtEdge .NE. ST_NOHANDLE) THEN
        CALL storage_realloc('afcstab_resizeStabilisation',&
            rafcstab%NEDGE, rafcstab%h_IverticesAtEdge,&
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF

      ! Resize array of subdiagonal edges
      IF (rafcstab%h_IsubdiagonalEdges .NE. ST_NOHANDLE) THEN
        CALL storage_realloc('afcstab_resizeStabilisation',&
            rafcstab%NEDGE, rafcstab%h_IsubdiagonalEdges,&
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF

      ! Resize array of edge data
      IF (rafcstab%h_DcoefficientsAtEdge .NE. ST_NOHANDLE) THEN
        CALL storage_realloc('afcstab_resizeStabilisation',&
            rafcstab%NEDGE, rafcstab%h_DcoefficientsAtEdge,&
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF

      ! Resize auxiliary edge vectors
      IF(ASSOCIATED(rafcstab%RedgeVectors)) THEN
        DO i = LBOUND(rafcstab%RedgeVectors,1),&
               UBOUND(rafcstab%RedgeVectors,1)
          CALL lsyssc_resizeVector(rafcstab%RedgeVectors(i),&
              rafcstab%NEDGE, .FALSE.)
        END DO
      END IF
    END IF
  END SUBROUTINE afcstab_resizeStabilisation

  !*****************************************************************************

!<subroutine>

  SUBROUTINE afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagonalEdgesIdx)

!<description>
    ! Returns a pointer to the index pointer for vertices at edge structure
!</description>

!<input>
    ! discrete operator
    TYPE(t_afcstab), INTENT(IN) :: rafcstab
!</input>

!<output>
    ! Pointer to the index pointer for superdiagonal edge numbers.
    ! NULL() if the discrete operator does not provide it.
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_IsuperdiagonalEdgesIdx
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    IF ((rafcstab%h_IsuperdiagonalEdgesIdx .EQ. ST_NOHANDLE) .OR.&
        (rafcstab%NEQ                  .EQ. 0)) THEN
      NULLIFY(p_IsuperdiagonalEdgesIdx)
      RETURN
    END IF
    
    ! Get the array
    CALL storage_getbase_int(rafcstab%h_IsuperdiagonalEdgesIdx,&
        p_IsuperdiagonalEdgesIdx,rafcstab%NEQ+1)
  END SUBROUTINE afcstab_getbase_IsupdiagEdgeIdx

  !*****************************************************************************

!<subroutine>

  SUBROUTINE afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

!<description>
    ! Returns a pointer to the vertices at edge structure
!</description>

!<input>
    ! discrete operator
    TYPE(t_afcstab), INTENT(IN) :: rafcstab
!</input>

!<output>
    ! Pointer to the vertices at edge structure
    ! NULL() if the discrete operator does not provide it.
    INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    IF ((rafcstab%h_IverticesAtEdge .EQ. ST_NOHANDLE) .OR.&
        (rafcstab%NEDGE             .EQ. 0)) THEN
      NULLIFY(p_IverticesAtEdge)
      RETURN
    END IF
    
    ! Get the array
    CALL storage_getbase_int2D(rafcstab%h_IverticesAtEdge,&
        p_IverticesAtEdge,rafcstab%NEDGE)
  END SUBROUTINE afcstab_getbase_IverticesAtEdge

  !*****************************************************************************

!<subroutine>

  SUBROUTINE afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagonalEdgesIdx)

!<description>
    ! Returns a pointer to the index pointer for 
    ! subdiagonal edge numbers
!</description>

!<input>
    ! discrete operator
    TYPE(t_afcstab), INTENT(IN) :: rafcstab
!</input>

!<output>
    ! Pointer to the index pointer for 
    ! subdiagonal edge numbers
    ! NULL() if the discrete operator does not provide it.
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_IsubdiagonalEdgesIdx
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    IF ((rafcstab%h_IsubdiagonalEdgesIdx .EQ. ST_NOHANDLE) .OR.&
        (rafcstab%NEQ                    .EQ. 0)) THEN
      NULLIFY(p_IsubdiagonalEdgesIdx)
      RETURN
    END IF
    
    ! Get the array
    CALL storage_getbase_int(rafcstab%h_IsubdiagonalEdgesIdx,&
        p_IsubdiagonalEdgesIdx,rafcstab%NEQ+1)
  END SUBROUTINE afcstab_getbase_IsubdiagEdgeIdx

  !*****************************************************************************

!<subroutine>

  SUBROUTINE afcstab_getbase_IsubdiagEdge(rafcstab,p_IsubdiagonalEdges)

!<description>
    ! Returns a pointer to the subdiagonal edge number
!</description>

!<input>
    ! discrete operator
    TYPE(t_afcstab), INTENT(IN) :: rafcstab
!</input>

!<output>
    ! Pointer to the subdiagonal edge numbers
    ! NULL() if the discrete operator does not provide it.
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_IsubdiagonalEdges
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    IF ((rafcstab%h_IsubdiagonalEdges .EQ. ST_NOHANDLE) .OR.&
        (rafcstab%NEDGE               .EQ. 0)) THEN
      NULLIFY(p_IsubdiagonalEdges)
      RETURN
    END IF
    
    ! Get the array
    CALL storage_getbase_int(rafcstab%h_IsubdiagonalEdges,&
        p_IsubdiagonalEdges,rafcstab%NEDGE)
  END SUBROUTINE afcstab_getbase_IsubdiagEdge

  !*****************************************************************************

!<subroutine>

  SUBROUTINE afcstab_getbase_DcoeffsAtEdge(rafcstab,p_DcoefficientsAtEdge)

!<description>
    ! Returns a pointer to the double-valued edge data
!</description>

!<input>
    ! discrete operator
    TYPE(t_afcstab), INTENT(IN) :: rafcstab
!</input>

!<output>
    ! Pointer to the double-valued edge data
    ! NULL() if the discrete operator does not provide it.
    REAL(DP), DIMENSION(:,:), POINTER :: p_DcoefficientsAtEdge
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    IF ((rafcstab%h_DcoefficientsAtEdge .EQ. ST_NOHANDLE) .OR.&
        (rafcstab%NEDGE                .EQ. 0)) THEN
      NULLIFY(p_DcoefficientsAtEdge)
      RETURN
    END IF
    
    ! Get the array
    CALL storage_getbase_double2D(rafcstab%h_DcoefficientsAtEdge,&
        p_DcoefficientsAtEdge,rafcstab%NEDGE)
  END SUBROUTINE afcstab_getbase_DcoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  SUBROUTINE afcstab_generateSubdiagEdges(rafcstab)

!<description>
    ! This subroutine generates the subdiagonal edge number
    ! structure based on a given edge data structure.
!</description>

!<inputoutput>
    ! discrete operator
    TYPE(t_afcstab), INTENT(INOUT) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_IsuperdiagonalEdgesIdx
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER   :: p_IsubdiagonalEdges
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER   :: p_IsubdiagonalEdgesIdx
    INTEGER(PREC_MATIDX) :: iedge,nedge,istor
    INTEGER(PREC_VECIDX) :: ieq,jeq,neq
    INTEGER(I32)         :: isize

    ! Check if edge-based data structure is prepared
    IF (IAND(rafcstab%iSpec,AFCSTAB_EDGESTRUCTURE) .EQ. 0) THEN
      CALL output_line('Discrete operator does not provide required &
          &edge-based data structure',OU_CLASS_ERROR,OU_MODE_STD,&
          'afcstab_generateSubdiagEdge')
      CALL sys_halt()
    END IF

    ! store dimensions of discrete operator
    neq   = rafcstab%NEQ
    nedge = rafcstab%NEDGE

    ! Allocate memory (if required)
    IF (rafcstab%h_IsubdiagonalEdgesIdx .EQ. ST_NOHANDLE) THEN
      CALL storage_new('afcstab_generateSubdiagEdges','IsubdiagonalEdgesIdx',&
          neq+1,ST_INT,rafcstab%h_IsubdiagonalEdgesIdx,ST_NEWBLOCK_ZERO)
    ELSE
      CALL storage_getsize(rafcstab%h_IsubdiagonalEdgesIdx,isize)
      IF (isize < neq+1) THEN
        CALL storage_realloc('afcstab_generateSubdiagEdges',neq+1,&
            rafcstab%h_IsubdiagonalEdgesIdx,ST_NEWBLOCK_ZERO,.FALSE.)
      ELSE
        CALL storage_clear(rafcstab%h_IsubdiagonalEdgesIdx)
      END IF
    END IF

    IF (rafcstab%h_IsubdiagonalEdges .EQ. ST_NOHANDLE) THEN
      CALL storage_new('afcstab_generateSubdiagEdges','IsubdiagonalEdges',&
          nedge,ST_INT,rafcstab%h_IsubdiagonalEdges,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(rafcstab%h_IsubdiagonalEdges,isize)
      IF (isize < nedge) THEN
        CALL storage_realloc('afcstab_generateSubdiagEdges',nedge,&
            rafcstab%h_IsubdiagonalEdges,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF
    
    ! Set pointers
    CALL afcstab_getbase_IsupdiagEdgeIdx(rafcstab,p_IsuperdiagonalEdgesIdx)
    CALL afcstab_getbase_IverticesAtEdge(rafcstab,p_IverticesAtEdge)
    CALL afcstab_getbase_IsubdiagEdgeIdx(rafcstab,p_IsubdiagonalEdgesIdx)
    CALL afcstab_getbase_IsubdiagEdge(rafcstab,p_IsubdiagonalEdges)
    
    ! Count number of superdiagonal edges
    DO ieq = 1, neq
      DO iedge = p_IsuperdiagonalEdgesIdx(ieq),&
                 p_IsuperdiagonalEdgesIdx(ieq+1)-1

        ! Determine that end-point of the edge which is not equal to IEQ
        jeq = p_IverticesAtEdge(1,iedge) + p_IverticesAtEdge(2,iedge) - ieq + 1

        ! Increase number of edges connected to this point by one
        p_IsubdiagonalEdgesIdx(jeq) = p_IsubdiagonalEdgesIdx(jeq)+1
      END DO
    END DO
    
    ! Reshuffle pass 1. In addition, set the maximum number of edges
    ! adjacent to one vertex. That is, take the maximum value of the
    ! number of edges to the left of the diagonal plus the number of
    ! edges to the right of the diagonal.
    rafcstab%NNVEDGE = p_IsuperdiagonalEdgesIdx(2) - p_IsuperdiagonalEdgesIdx(1)
    
    DO ieq = 2, neq
      rafcstab%NNVEDGE = MAX(rafcstab%NNVEDGE,&
          p_IsuperdiagonalEdgesIdx(ieq+1) - p_IsuperdiagonalEdgesIdx(ieq))
      
      p_IsubdiagonalEdgesIdx(ieq) = p_IsubdiagonalEdgesIdx(ieq) +&
          p_IsubdiagonalEdgesIdx(ieq-1)
    END DO
    p_IsubdiagonalEdgesIdx(neq+1) = p_IsubdiagonalEdgesIdx(neq+1) +&
        p_IsubdiagonalEdgesIdx(neq)
    
    ! Store the subdiagonal edge numbers
    DO ieq = 1, neq
      
      ! For each equation, loop over the edges which are located in
      ! the upper right triangular matrix
      DO iedge = p_IsuperdiagonalEdgesIdx(ieq),&
                 p_IsuperdiagonalEdgesIdx(ieq+1)-1
        
        ! Determine that end-point of the edge which is not equal to IEQ
        jeq = p_IverticesAtEdge(1,iedge) + p_IverticesAtEdge(2,iedge) - ieq
        
        ! Determine next free position
        istor = p_IsubdiagonalEdgesIdx(jeq)+1
        
        p_IsubdiagonalEdgesIdx(jeq) = istor
        p_IsubdiagonalEdges(istor)  = iedge
      END DO
    END DO
    
    ! Reshuffle pass 2: Adjust the index vector which was tainted in
    ! the above loop
    DO ieq = neq+1, 2, -1
      p_IsubdiagonalEdgesIdx(ieq) = p_IsubdiagonaledgesIdx(ieq-1)+1
    END DO
    p_IsubdiagonalEdgesIdx(1) = 1

    ! Set specifier for extended edge structure
    rafcstab%iSpec = IOR(rafcstab%iSpec,AFCSTAB_SUBDIAGONALEDGES)
  END SUBROUTINE afcstab_generateSubdiagEdges
  
  !*****************************************************************************

!<subroutine>

  SUBROUTINE afcstab_generateExtSparsity(rmatrixSrc, rmatrixExtended)

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
    TYPE(t_matrixScalar), INTENT(IN)    :: rmatrixSrc
!</input>

!<inputoutput>
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixExtended
!</inputoutput>
!</subroutine>
    
    ! Clear output matrix
    CALL lsyssc_releaseMatrix(rmatrixExtended)
    
    ! Compute Z=A*A and let the connectivity graph of Z be the
    ! extended sparsity pattern of the Jacobian matrix
    CALL lsyssc_multMatMat(rmatrixSrc, rmatrixSrc,&
        rmatrixExtended, .TRUE., .TRUE., .FALSE.)
  END SUBROUTINE afcstab_generateExtSparsity

  !*****************************************************************************

!<function>
  
  ELEMENTAL FUNCTION afcstab_limit_unbounded(p, q, default) RESULT(r)

!<description>
    ! This function computes the ratio Q/P. If the denominator is
    ! too small, then the default value is applied.
!</description>

!<input>
    ! (de)nominator
    REAL(DP), INTENT(IN) :: p,q

    ! default value
    REAL(DP), INTENT(IN) :: default
!</input>

!<result>
    ! limited ratio
    REAL(DP) :: r
!</result>
!</function>

    IF (p > SYS_EPSREAL) THEN
      r = q/p
    ELSE
      r = default
    END IF
  END FUNCTION afcstab_limit_unbounded

  !*****************************************************************************

!<function>
  
  ELEMENTAL FUNCTION afcstab_limit_bounded(p, q, default, dbound) RESULT(r)

!<description>
    ! This function computes the limited ratio Q/P and bounds the
    ! result by the size of dbound. If the denominator is too small
    ! then the default value is applied.
!</description>

!<input>
    ! (de)nominator
    REAL(DP), INTENT(IN) :: p,q
    
    ! default value
    REAL(DP), INTENT(IN) :: default

    ! upper bound
    REAL(DP), INTENT(IN) :: dbound
!</input>

!<result>
    ! limited ratio
    REAL(DP) :: r
!</result>
!</function>
    
    IF (p > SYS_EPSREAL) THEN
      r = MIN(q/p, dbound)
    ELSE
      r = default
    END IF
  END FUNCTION afcstab_limit_bounded
END MODULE afcstabilisation
