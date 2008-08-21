!##############################################################################
!# ****************************************************************************
!# <name> meshadjacency </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of algorithms which are applied onto an adjacency
!# graph which is described by two arrays. An adjacency graph is similar to the
!# sparsity pattern of a CSR matrix (a.k.a. matrix storage technique 7/9),
!# describing which nodes are adjacent to each other in respect to a specific
!# adjacency relation.
!# To make the algorithms usable for multiple purposes, no adjacency graph
!# data structure is defined in this module - the algorithms accept two
!# storage handles to the adjacency graph arrays instead.
!#
!# Remarks:
!# - Unless stated otherwise, the adjacency graph is assumed to be undirected,
!#   i.e. if a node A is adjacent to node B, then node B is also adjacent to A.
!#
!# - In contrast to the column index array of a CSR matrix, the indices of the
!#   adjacencent nodes for a node are not necessarily in ascending order.
!#   However, this can be achieved by applying the 'adj_sortAdjacencies'
!#   routine onto the adjacency graph arrays.
!#
!# The following routines can be found here:
!#
!# 1.) adj_applyPermutation
!#     -> Applies a permutation onto the nodes of an adjacency graph.
!#        This is an interface for two versions, one for undirected and one
!#        for directed graphs. The 'directed' version allows the two adjacency
!#        arrays to be permuted differently, while the 'undirected' version
!#        applies the same permutation on both adjacency arrays. 
!#
!# 2.) adj_sortAdjacencies
!#     -> Sorts the the indices of the adjacent nodes for each node to be
!#        in ascending order.
!#
!# 3.) adj_calcColouring
!#     -> Calculates a colouring for the nodes of the adjacency graph such
!#        that any two adjacent nodes have different colours.
!#
!# 4.) adj_calcCuthillMcKee
!#     -> Calculates a (reverse) Cuthill-McKee ordering for the nodes of
!#        the adjacency graph. 
!# </purpose>
!##############################################################################

MODULE adjacency

  USE fsystem
  USE storage

IMPLICIT NONE

!<constants>

!<constantblock description="flags for the Cuthill-McKee algorithm">
  
  ! If this flag is set, then the nodes on each level are sorted by their
  ! degree in ascending order.
  ! Cannot be used together with ADJ_CMK_FLAG_SORT_MAX.
  INTEGER(I32), PARAMETER :: ADJ_CMK_FLAG_SORT_MIN = 1
  
  ! If this flag is set, then the nodes on each level are sorted by their
  ! degree in descending order.
  ! Cannot be used together with ADJ_CMK_FLAG_SORT_MIN.
  INTEGER(I32), PARAMETER :: ADJ_CMK_FLAG_SORT_MAX = 2
  
  ! If this flag is set, then a root of minimum degree is chosen.
  ! Cannot be used together with ADJ_CMK_FLAG_ROOT_MAX.
  INTEGER(I32), PARAMETER :: ADJ_CMK_FLAG_ROOT_MIN = 4

  ! If this flag is set, then a root of maximum degree is chosen.
  ! Cannot be used together with ADJ_CMK_FLAG_ROOT_MIN.
  INTEGER(I32), PARAMETER :: ADJ_CMK_FLAG_ROOT_MAX = 8

  ! If this flag is set, then the Cuthill-McKee ordering is reversed.
  INTEGER(I32), PARAMETER :: ADJ_CMK_FLAG_REVERSE = 16

!</constantblock>

!</constants>

  INTERFACE adj_applyPermutation
    MODULE PROCEDURE adj_applyPermutation_undirected
    MODULE PROCEDURE adj_applyPermutation_directed
  END INTERFACE

CONTAINS

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE adj_applyPermutation_undirected(h_Iptr, h_Iidx, h_Ipermute)

!<description>
  ! This routine applies a permutation onto an (undirected) adjacency graph.
  ! Both the pointer and index arrays are permuted using the same permutation.
!</description>

!<input>
  ! A storage handle to the permutation that is to be applied onto the
  ! adjacency graph.
  INTEGER(I32), INTENT(IN) :: h_Ipermute
!</input>

!<inputoutput>
  ! A storage handle to the pointer-array of the adjacency graph.
  INTEGER(I32), INTENT(INOUT) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  INTEGER(I32), INTENT(INOUT) :: h_Iidx
!</inputoutput>

!</subroutine>

  INTEGER :: i,j,k,m, NEL, NNZE
  INTEGER, DIMENSION(:), ALLOCATABLE :: IptrOld, IidxOld,Iinv
  INTEGER, DIMENSION(:), POINTER :: p_Iptr, p_Iidx, p_Ipermute
  
    ! First of all, get the arrays from the storage
    CALL storage_getbase_int(h_Iptr, p_Iptr)
    CALL storage_getbase_int(h_Iidx, p_Iidx)
    CALL storage_getbase_int(h_Ipermute, p_Ipermute)
    
    ! Get the length of the arrays
    NEL = UBOUND(p_Iptr,1)-1
    NNZE = UBOUND(p_Iidx,1)
    
    ! Create a backup of the arrays
    ALLOCATE(IptrOld(NEL+1))
    ALLOCATE(IidxOld(NNZE))
    DO i = 1, NEL+1
      IptrOld(i) = p_Iptr(i)
    END DO
    DO i = 1, NNZE
      IidxOld(i) = p_Iidx(i)
    END DO
    
    ! Calculate the inverse permutation
    ALLOCATE(Iinv(NEL))
    DO i = 1, NEL
      Iinv(p_Ipermute(i)) = i
    END DO
    
    ! Now go through all elements
    j = 1
    DO i = 1, NEL
    
      ! Set the pointer for this element
      p_Iptr(i) = j
      
      ! Get the index of the element
      m = p_Ipermute(i)
      
      ! And copy the indices
      DO k = IptrOld(m), IptrOld(m+1)-1
        p_Iidx(j) = Iinv(IidxOld(k))
        j = j+1
      END DO
    END DO
    p_Iptr(NEL+1) = j
    
    ! Release the temporary memory
    DEALLOCATE(Iinv)
    DEALLOCATE(IidxOld)
    DEALLOCATE(IptrOld)
    
    ! That's it

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE adj_applyPermutation_directed(h_Iptr, h_Iidx, h_IpermutePtr,&
                                           h_IpermuteIdx)

!<description>
  ! This routine applies a permutation onto an (directed) adjacency graph.
  ! The pointer and index arrays may be permuted using different permutations.
!</description>

!<input>
  ! A storage handle to the permutation that is to be applied onto the pointer
  ! array of the adjacency graph. May be ST_NOHANDLE if the pointer array
  ! is not to be permuted.
  INTEGER(I32), INTENT(IN) :: h_IpermutePtr

  ! A storage handle to the permutation that is to be applied onto the index
  ! array of the adjacency graph. May be ST_NOHANDLE if the index array
  ! is not to be permuted.
  INTEGER(I32), INTENT(IN) :: h_IpermuteIdx
!</input>

!<inputoutput>
  ! A storage handle to the pointer-array of the adjacency graph.
  INTEGER(I32), INTENT(INOUT) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  INTEGER(I32), INTENT(INOUT) :: h_Iidx
!</inputoutput>

!</subroutine>

  INTEGER :: i,j,k,m, NNODE, NNZE
  INTEGER, DIMENSION(:), ALLOCATABLE :: IptrOld, IidxOld,IinvIdx
  INTEGER, DIMENSION(:), POINTER :: p_Iptr, p_Iidx, p_IpermutePtr,&
    p_IpermuteIdx
  
    ! First of all, get the arrays from the storage
    CALL storage_getbase_int(h_Iptr, p_Iptr)
    CALL storage_getbase_int(h_Iidx, p_Iidx)

    ! Get the length of the arrays
    NNODE = UBOUND(p_Iptr,1)-1
    NNZE = UBOUND(p_Iidx,1)

    ! Get the permutation arrays
    IF(h_IpermutePtr .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(h_IpermutePtr, p_IpermutePtr)
    ELSE
      NULLIFY(p_IpermutePtr)
    END IF
    IF(h_IpermuteIdx .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(h_IpermuteIdx, p_IpermuteIdx)
      
      ! Build the inverse permutation
      k = UBOUND(p_IpermuteIdx,1)
      ALLOCATE(IinvIdx(k))
      DO i = 1, k
        IinvIdx(p_IpermuteIdx(i)) = i
      END DO
    ELSE
      NULLIFY(p_IpermuteIdx)
    END IF
    
    ! What type of permutation is to be applied?
    IF(ASSOCIATED(p_IpermutePtr)) THEN
    
      ! Create a backup of the arrays
      ALLOCATE(IptrOld(NNODE+1))
      ALLOCATE(IidxOld(NNZE))
      DO i = 1, NNODE+1
        IptrOld(i) = p_Iptr(i)
      END DO
      DO i = 1, NNZE
        IidxOld(i) = p_Iidx(i)
      END DO

      IF(ASSOCIATED(p_IpermuteIdx)) THEN
        ! We permute both arrays.
        j = 1
        DO i = 1, NNODE
        
          ! Set the pointer for this element
          p_Iptr(i) = j
          
          ! Get the index of the element
          m = p_IpermutePtr(i)
          
          ! And copy the indices
          DO k = IptrOld(m), IptrOld(m+1)-1
            p_Iidx(j) = IinvIdx(IidxOld(k))
            j = j+1
          END DO
        END DO
        p_Iptr(NNODE+1) = j

      ELSE
        ! We permute the pointer array.
        j = 1
        DO i = 1, NNODE
        
          ! Set the pointer for this element
          p_Iptr(i) = j
          
          ! Get the index of the element
          m = p_IpermutePtr(i)
          
          ! And copy the indices
          DO k = IptrOld(m), IptrOld(m+1)-1
            p_Iidx(j) = IidxOld(k)
            j = j+1
          END DO
        END DO
        p_Iptr(NNODE+1) = j
      
      END IF
    
      ! Deallocate the backup
      DEALLOCATE(IidxOld)
      DEALLOCATE(IptrOld)

    ELSE IF(ASSOCIATED(p_IpermuteIdx)) THEN
      
      ! We permute the index array
      DO i = 1, NNZE
        p_Iidx(i) = IinvIdx(p_Iidx(i))
      END DO
      
    END IF
    
    ! Release the temporary memory
    IF(ALLOCATED(IinvIdx)) &
      DEALLOCATE(IinvIdx)
    
    ! That's it

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE adj_sortAdjacencies(h_Iptr, h_Iidx)

!<description>
  ! This routine sorts the index array such that the indices of the adjacent
  ! nodes are in ascending order.
!</description>

!<input>
  ! A storage handle to the pointer-array of the adjacency graph.
  INTEGER(I32), INTENT(IN) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  INTEGER(I32), INTENT(IN) :: h_Iidx
!</input>

  INTEGER, DIMENSION(:), POINTER :: p_Iptr, p_Iidx
  INTEGER :: i,j,k,m,NNODE
  
    ! Get the pointers from the storage
    CALL storage_getbase_int(h_Iptr, p_Iptr)
    CALL storage_getbase_int(h_Iidx, p_Iidx)
    
    NNODE = UBOUND(p_Iptr,1)-1
    
    ! The sorting algorithm used here is bubblesort - it might (and should)
    ! be replaced by a more efficient algorithm later...
    
    ! Go through all nodes
    DO i = 1, NNODE
      DO j = p_Iptr(i), p_Iptr(i+1)-1
        DO k = j+1, p_Iptr(i+1)-1
          IF(p_Iidx(j) .GT. p_Iidx(k)) THEN
            m = p_Iidx(j)
            p_Iidx(j) = p_Iidx(k)
            p_Iidx(k) = m
          END IF
        END DO ! k
      END DO ! j
    END DO ! i

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE adj_calcColouring(h_Iptr,h_Iidx,h_Icpt,h_Ipermute,h_InodeOrder)

!<description>
  ! This routine calculates a colouring for an adjacency graph, i.e. a mapping
  !
  !                C: {1,...,#NODES} -> {1,...,#COLORS}
  !
  ! such that the following estaminate holds:
  ! If two nodes a and b are adjacent, then C(a) != C(b).
  !
  ! Instead of returning the colour mapping C, this routine returns the handles
  ! to a colour partition table and a permutation:
  ! For a given colour i,
  !
  !       S := {p_Ipermute(p_Icpt(i)),...,p_Ipermute(p_Icpt(i+1)-1)}
  !
  ! is a set of node indices of the same colour, i.e.:
  !
  !         For all a, b in S, a != b ==> a is not adjacent to b
  !
  ! The number of colours used is implicitly given by UBOUND(p_Icpt,1)-1.
!</description>

!<input>
  ! A storage handle to the pointer-array of the adjacency graph.
  INTEGER(I32), INTENT(IN) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  INTEGER(I32), INTENT(IN) :: h_Iidx
  
  ! OPTIONAL: A storage handle to an array holding the order in which the
  ! colouring algorithm should proceed through the nodes. If not given,
  ! canonical order is used.
  INTEGER(I32), OPTIONAL, INTENT(IN) :: h_InodeOrder
!</input>

!<output>
  ! A storage handle to the permutation that has to be applied to
  ! get a coloured partitioning.
  INTEGER(I32), INTENT(OUT) :: h_Ipermute

  ! A storage handle to the colour partition table.
  ! The number of colours used is implicitly given by UBOUND(Icpt,1)-1.
  INTEGER(I32), INTENT(OUT) :: h_Icpt
!</output>

!</subroutine>

  INTEGER, DIMENSION(:), POINTER :: p_Iptr, p_Iidx, p_Icpt, p_Ipermute,&
    p_InodeOrder
  INTEGER, DIMENSION(:), ALLOCATABLE :: InodeColour, IcolourMap, InumNodes
  INTEGER :: i,j,k,inode,DEGREE,NNODE,inumColours,iminColour,iminNodes
  
    ! Get the pointers from the storage
    CALL storage_getbase_int(h_Iptr, p_Iptr)
    CALL storage_getbase_int(h_Iidx, p_Iidx)
    
    ! First of all, we need to calculate the degree of the graph, as this
    ! is the upper bound for the number of colours we will need.
    NNODE = UBOUND(p_Iptr,1)-1
    DEGREE = 1
    DO i = 1, NNODE
      DEGREE = MAX(DEGREE, p_Iptr(i+1)-p_Iptr(i))
    END DO
    
    ! Allocate the auxiliary arrays
    ALLOCATE(InodeColour(NNODE))
    ALLOCATE(IcolourMap(DEGREE))
    ALLOCATE(InumNodes(DEGREE))
    DO i = 1, DEGREE
      InumNodes(i) = 0
    END DO
    inumColours = 0
    
    ! Do we have a node order given?
    IF(PRESENT(h_InodeOrder)) THEN
    
      ! Get the node order array from the storage then
      CALL storage_getbase_int(h_InodeOrder, p_InodeOrder)
      
      ! And format the node-colour array
      DO i = 1, NNODE
        InodeColour(i) = 0
      END DO

      ! Now let's loop through the nodes
      DO i = 1, NNODE
      
        ! Get the node index
        inode = p_InodeOrder(i)
        
        ! Format the colour map
        DO j = 1, inumColours
          IcolourMap(j) = 0
        END DO
        
        ! Go through the adjacencies of this node
        DO j = p_Iptr(inode), p_Iptr(inode+1)-1
          
          ! The the node's colour
          k = InodeColour(p_Iidx(j))
          
          IF (k .GT. 0) THEN
            ! Mark this node's colour as 'used'
            IcolourMap(k) = 1
          END IF
          
        END DO ! j
        
        ! Now find a free colour with minimum nodes
        iminColour = -1
        iminNodes = NNODE + 100
        DO j = 1, inumColours
        
          ! Skip this colour if it is used
          IF(IcolourMap(j) .NE. 0) CYCLE
          
          ! Is this the colour with minimum nodes?
          IF(InumNodes(j) .LT. iminNodes) THEN
            iminNodes = InumNodes(j)
            iminColour = j
          END IF
          
        END DO ! j
        
        ! Did we find a free colour?
        IF(iminColour .GT. 0) THEN
          ! Yes, so set the node's colour to this one
          InodeColour(inode) = iminColour
          InumNodes(iminColour) = InumNodes(iminColour) + 1
        ELSE
          ! No, so add a new colour
          inumColours = inumColours+1
          InodeColour(inode) = inumColours
          InumNodes(inumColours) = 1
        END IF
      
      END DO ! i
    
    ELSE
    
      ! Now let's loop through the nodes
      DO i = 1, NNODE
        
        ! Format the colour map
        DO j = 1, inumColours
          IcolourMap(j) = 0
        END DO
        
        ! Go through the adjacencies of this node
        DO j = p_Iptr(i), p_Iptr(i+1)-1
          
          ! The the node index
          k = p_Iidx(j)
          
          IF (k .LT. i) THEN
            ! Mark this node's colour as 'used'
            IcolourMap(InodeColour(k)) = 1
          END IF
          
        END DO ! j
        
        ! Now find a free colour with minimum entries
        iminColour = -1
        iminNodes = NNODE + 100
        DO j = 1, inumColours
        
          ! Skip this colour if it is used
          IF(IcolourMap(j) .NE. 0) CYCLE
          
          ! Is this the colour with minimum nodes?
          IF(InumNodes(j) .LT. iminNodes) THEN
            iminNodes = InumNodes(j)
            iminColour = j
          END IF
          
        END DO ! j
        
        ! Did we find a free colour?
        IF(iminColour .GT. 0) THEN
          ! Yes, so set the node's colour to this one
          InodeColour(i) = iminColour
          InumNodes(iminColour) = InumNodes(iminColour) + 1
        ELSE
          ! No, so add a new colour
          inumColours = inumColours+1
          InodeColour(i) = inumColours
          InumNodes(inumColours) = 1
        END IF
      
      END DO ! i
    
    END IF
    
    ! Now IcellColour contains a valid colouring of our adjacency graph.
    ! We now need to build the colour partition table...
    CALL storage_new ('adj_calcColouring', 'p_Icpt', inumColours+1, ST_INT,&
                      h_Icpt, ST_NEWBLOCK_NOINIT)
    
    ! Get the array from the storage
    CALL storage_getbase_int(h_Icpt, p_Icpt)

    ! Go through the colours
    p_Icpt(1) = 1
    DO i = 1, inumColours
      p_Icpt(i+1) = p_Icpt(i) + InumNodes(i)
    END DO
    
    ! Now make a backup of the cpt - we will need it to build the permutation
    DO i = 1, inumColours
      IcolourMap(i) = p_Icpt(i)
    END DO
    
    ! Build the inverse permutation array - we will abuse the InodeColour array
    ! to store the inverse permutation.
    DO i = 1, NNODE
      k = InodeColour(i)
      InodeColour(i) = IcolourMap(k)
      IcolourMap(k) = IcolourMap(k)+1
    END DO

    ! Allocate an array for the permuation
    CALL storage_new ('adj_calcColouring', 'p_Ipermute', NNODE, ST_INT,&
                      h_Ipermute, ST_NEWBLOCK_NOINIT)
    
    ! Get the array from the storage
    CALL storage_getbase_int(h_Ipermute, p_Ipermute)
    
    ! Calculate the permutation
    DO i = 1, NNODE
      p_Ipermute(InodeColour(i)) = i
    END DO
    
    ! And release all temporary memory
    DEALLOCATE(InodeColour)
    DEALLOCATE(InumNodes)
    DEALLOCATE(IcolourMap)
    
    ! That's it

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE adj_calcCuthillMcKee(h_Iptr,h_Iidx,h_Ipermute,cflags)

!<description>
  ! This routine calculates a (reverse) Cuthill-McKee ordering permutation
  ! for a given adjacency graph.
!</description>

!<input>
  ! A storage handle to the pointer-array of the adjacency graph.
  INTEGER(I32), INTENT(IN) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  INTEGER(I32), INTENT(IN) :: h_Iidx
  
  ! OPTIONAL: A combination ADJ_CMK_FLAG_XXXX constants defined at the top
  ! of this module specifying the flags for the algorithm. If not given,
  ! only the ADJ_CMK_FLAG_ROOT_MIN flag is set.
  INTEGER(I32), OPTIONAL, INTENT(IN) :: cflags
!</input>

!<output>
  ! A storage handle to the (reverse) Cuthill-McKee permutation.
  INTEGER(I32), INTENT(OUT) :: h_Ipermute
!</output>

!</subroutine>

  INTEGER, DIMENSION(:), POINTER :: p_Iptr, p_Iidx, p_Ipermute
  INTEGER, DIMENSION(:), ALLOCATABLE :: Iaux
  INTEGER :: i,j,k,n,NNODE,root,lvl1,lvl2,lvl3,lvl_root,sort_choice, root_choice
  LOGICAL :: bReverse
  
    ! Reset output handle
    h_Ipermute = ST_NOHANDLE
  
    ! Set the default strategy
    sort_choice = 0
    root_choice = ADJ_CMK_FLAG_ROOT_MIN
    bReverse = .FALSE.
    IF(PRESENT(cflags)) THEN
      ! Check if the flags are valid
      IF(((IAND(cflags,ADJ_CMK_FLAG_SORT_MIN) .NE. 0) .AND. &
          (IAND(cflags,ADJ_CMK_FLAG_SORT_MAX) .NE. 0)) .OR. &
         ((IAND(cflags,ADJ_CMK_FLAG_ROOT_MIN) .NE. 0) .AND. &
          (IAND(cflags,ADJ_CMK_FLAG_ROOT_MAX) .NE. 0))) THEN
         PRINT *, 'ERROR: adj_calcCuthillMcKee'
         PRINT *, 'invalid flag combination'
         RETURN
      END IF
      sort_choice = IAND(cflags,IOR(ADJ_CMK_FLAG_SORT_MIN,ADJ_CMK_FLAG_SORT_MAX))
      root_choice = IAND(cflags,IOR(ADJ_CMK_FLAG_ROOT_MIN,ADJ_CMK_FLAG_ROOT_MAX))
      bReverse = (IAND(cflags,ADJ_CMK_FLAG_REVERSE) .NE. 0)
    END IF
    
    ! Get the pointers from the storage
    CALL storage_getbase_int(h_Iptr, p_Iptr)
    CALL storage_getbase_int(h_Iidx, p_Iidx)
    NNODE = UBOUND(p_Iptr,1)-1
    
    ! Allocate the permutation array
    CALL storage_new ('adj_calcCuthillMcKee', 'p_Ipermute', NNODE, ST_INT, &
                      h_Ipermute, ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int(h_Ipermute,p_Ipermute)
    
    ! Allocate an auxiliary array
    ALLOCATE(Iaux(NNODE))
    
    ! Calculate the degree for all nodes
    DO i = 1, NNODE
      Iaux(i) = p_Iptr(i+1) - p_Iptr(i)
    END DO

    lvl1 = 0
    lvl2 = 0
    DO WHILE(lvl2 .LT. NNODE)
    
      ! Usually, the DO-WHILE loop we're in performs just one iteration.
      ! However, if the adjacency graph is reducible, then we cannot reach all
      ! nodes by running through the adjacency levels of only one root - so we
      ! are forced to choose a new root for every disjunct subgraph...
      root = -1
      
      ! Now how do we choose our root node?
      SELECT CASE(root_choice)
      CASE (ADJ_CMK_FLAG_ROOT_MIN)
        ! Choose a free node with minimum degree as root
        n = NNODE+100
        DO i = 1, NNODE
          IF((Iaux(i) .GT. 0) .AND. (Iaux(i) .LT. n)) THEN
            root = i
            n = Iaux(i)
          END IF
        END DO

      CASE (ADJ_CMK_FLAG_ROOT_MAX)
        ! Choose a free node with maximum degree as root
        n = 0
        DO i = 1, NNODE
          IF(Iaux(i) .GT. n) THEN
            root = i
            n = Iaux(i)
          END IF
        END DO
      
      CASE DEFAULT
        ! Choose the first free node as 'root'
        DO i = 1, NNODE
          IF(Iaux(i) .GT. 0) THEN
            root = i
            EXIT
          END IF
        END DO
        
      END SELECT

      IF(root .LE. 0) THEN
      
        ! If we come out here, then something must have gone terribly wrong.
        ! We haven't processed all nodes yet, but we cannot find a new root...
        PRINT *, 'ERROR: adj_calcCuthillMcKee'
        PRINT *, 'internal error - no root found'
        CALL sys_halt()
        
      END IF
      
      ! Initialise the first level for the root node
      lvl1 = lvl1+1
      lvl2 = lvl1
      lvl3 = lvl1
      
      ! Add the root into the permutation array
      lvl_root = lvl1
      p_Ipermute(lvl1) = root
      Iaux(root) = -Iaux(root)
    
      ! Now let's go through the adjacency levels of the root, until there
      ! are no more adjacent nodes left.
      DO WHILE(lvl2 .LT. NNODE)
      
        ! Go through all nodes in the current level, and fetch the node
        ! indices for the next level
        DO i = lvl1, lvl2
        
          ! Get the index of the node
          n = p_Ipermute(i)
        
          ! Go through all nodes which are adjacent to node n
          DO j = p_Iptr(n), p_Iptr(n+1)-1
          
            ! Get the index of the adjacent node
            k = p_Iidx(j)
            
            ! Has this node already been processed?
            IF(Iaux(k) .GT. 0) THEN
            
              ! No, so let's add it to the next level
              lvl3 = lvl3+1
              p_Ipermute(lvl3) = k
              
              ! And negate it's entry in the auxiliary array to mark it
              ! as 'already processed'
              Iaux(k) = -Iaux(k)
              
            END IF
            
          END DO ! j
        
        END DO ! i
        
        ! Didn't we find any adjacent nodes anymore?
        ! Then jump out of this loop
        IF(lvl3 .LE. lvl2) EXIT

        ! Do we have to sort the permutation array for the next level?
        SELECT CASE(sort_choice)
        CASE (ADJ_CMK_FLAG_SORT_MIN)
          CALL adj_cmk_aux_sort_min(lvl2+1,lvl3,p_Ipermute,Iaux)
        CASE (ADJ_CMK_FLAG_SORT_MAX)
          CALL adj_cmk_aux_sort_max(lvl2+1,lvl3,p_Ipermute,Iaux)
        END SELECT
        
        ! Otherwise let's proceed with the next adjacency level
        lvl1 = lvl2+1
        lvl2 = lvl3
      
      END DO ! WHILE(lvl2 .LT. NNODE)
      
      ! Do we have to reverse the ordering?
      IF(bReverse) THEN
        i = lvl_root
        j = lvl2
        DO WHILE(i .LT. j)
          n = p_Ipermute(i)
          p_Ipermute(i) = p_Ipermute(j)
          p_Ipermute(j) = n
          i = i+1
          j = j-1
        END DO
      END IF
      
    END DO ! WHILE(lvl2 .LT. NNODE)
    
    ! Deallocate the auxiliary array
    DEALLOCATE(Iaux)
    
    ! That's it

  CONTAINS
  
    PURE SUBROUTINE adj_cmk_aux_sort_min(m,n,Iperm,Iaux)
    INTEGER, INTENT(IN) :: m,n
    INTEGER, DIMENSION(:), INTENT(INOUT) :: Iperm, Iaux
    
    INTEGER :: i,j,k
    
      ! Return if there's nothing to do here...
      IF(m .LE. n) RETURN
    
      ! Bubble-sort - replace with something faster later...
      
      ! Remark:
      ! Keep in mind that the Iaux array contains the negative degree of
      ! the nodes, as the nodes have already been processed! So we need
      ! to sort the nodes such that their entries in the Iaux array are
      ! in descending order instead of ascending order!
      DO i = m, n
        DO j = i+1, n
          IF(Iaux(Iperm(j)) .GT. Iaux(Iperm(i))) THEN
            k = Iperm(j)
            Iperm(j) = Iperm(i)
            Iperm(i) = k
          END IF
        END DO
      END DO
    
    END SUBROUTINE ! adj_cmk_aux_sort_min
  
    PURE SUBROUTINE adj_cmk_aux_sort_max(m,n,Iperm,Iaux)
    INTEGER, INTENT(IN) :: m,n
    INTEGER, DIMENSION(:), INTENT(INOUT) :: Iperm, Iaux
    
    INTEGER :: i,j,k
    
      ! Return if there's nothing to do here...
      IF(m .LE. n) RETURN
    
      ! Bubble-sort - replace with something faster later...
      DO i = m, n
        DO j = i+1, n
          IF(Iaux(Iperm(j)) .LT. Iaux(Iperm(i))) THEN
            k = Iperm(j)
            Iperm(j) = Iperm(i)
            Iperm(i) = k
          END IF
        END DO
      END DO
    
    END SUBROUTINE ! adj_cmk_aux_sort_max
  
  END SUBROUTINE ! adj_calcCuthillMcKee
  
END MODULE