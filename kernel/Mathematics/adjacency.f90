!##############################################################################
!# ****************************************************************************
!# <name> adjacency </name>
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

module adjacency

!$use omp_lib
  use fsystem
  use genoutput
  use storage

  implicit none
  
  private

!<constants>

!<constantblock description="flags for the Cuthill-McKee algorithm">
  
  ! If this flag is set, then the nodes on each level are sorted by their
  ! degree in ascending order.
  ! Cannot be used together with ADJ_CMK_FLAG_SORT_MAX.
  integer(I32), parameter, public :: ADJ_CMK_FLAG_SORT_MIN = 1
  
  ! If this flag is set, then the nodes on each level are sorted by their
  ! degree in descending order.
  ! Cannot be used together with ADJ_CMK_FLAG_SORT_MIN.
  integer(I32), parameter, public :: ADJ_CMK_FLAG_SORT_MAX = 2
  
  ! If this flag is set, then a root of minimum degree is chosen.
  ! Cannot be used together with ADJ_CMK_FLAG_ROOT_MAX.
  integer(I32), parameter, public :: ADJ_CMK_FLAG_ROOT_MIN = 4

  ! If this flag is set, then a root of maximum degree is chosen.
  ! Cannot be used together with ADJ_CMK_FLAG_ROOT_MIN.
  integer(I32), parameter, public :: ADJ_CMK_FLAG_ROOT_MAX = 8

  ! If this flag is set, then the Cuthill-McKee ordering is reversed.
  integer(I32), parameter, public :: ADJ_CMK_FLAG_REVERSE = 16

!</constantblock>

!</constants>

  interface adj_applyPermutation
    module procedure adj_applyPermutation_undirected
    module procedure adj_applyPermutation_directed
  end interface

  public :: adj_applyPermutation
  public :: adj_sortAdjacencies
  public :: adj_calcColouring
  public :: adj_calcCuthillMcKee

contains

  ! ***************************************************************************

!<subroutine>
  
  subroutine adj_applyPermutation_undirected(h_Iptr, h_Iidx, h_Ipermute)

!<description>
  ! This routine applies a permutation onto an (undirected) adjacency graph.
  ! Both the pointer and index arrays are permuted using the same permutation.
!</description>

!<input>
  ! A storage handle to the permutation that is to be applied onto the
  ! adjacency graph.
  integer, intent(in) :: h_Ipermute
!</input>

!<inputoutput>
  ! A storage handle to the pointer-array of the adjacency graph.
  integer, intent(in) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  integer, intent(inout) :: h_Iidx
!</inputoutput>

!</subroutine>

  integer :: i,j,k,m, NEL, NNZE
  integer, dimension(:), allocatable :: IptrOld, IidxOld,Iinv
  integer, dimension(:), pointer :: p_Iptr, p_Iidx, p_Ipermute
  
    ! First of all, get the arrays from the storage
    call storage_getbase_int(h_Iptr, p_Iptr)
    call storage_getbase_int(h_Iidx, p_Iidx)
    call storage_getbase_int(h_Ipermute, p_Ipermute)
    
    ! Get the length of the arrays
    NEL = ubound(p_Iptr,1)-1
    NNZE = ubound(p_Iidx,1)
    
    ! Create a backup of the arrays
    allocate(IptrOld(NEL+1))
    allocate(IidxOld(NNZE))
    do i = 1, NEL+1
      IptrOld(i) = p_Iptr(i)
    end do
    do i = 1, NNZE
      IidxOld(i) = p_Iidx(i)
    end do
    
    ! Calculate the inverse permutation
    allocate(Iinv(NEL))
    do i = 1, NEL
      Iinv(p_Ipermute(i)) = i
    end do
    
    ! Now go through all elements
    j = 1
    do i = 1, NEL
    
      ! Set the pointer for this element
      p_Iptr(i) = j
      
      ! Get the index of the element
      m = p_Ipermute(i)
      
      ! And copy the indices
      do k = IptrOld(m), IptrOld(m+1)-1
        p_Iidx(j) = Iinv(IidxOld(k))
        j = j+1
      end do
    end do
    p_Iptr(NEL+1) = j
    
    ! Release the temporary memory
    deallocate(Iinv)
    deallocate(IidxOld)
    deallocate(IptrOld)
    
    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine adj_applyPermutation_directed(h_Iptr, h_Iidx, h_IpermutePtr,&
                                           h_IpermuteIdx)

!<description>
  ! This routine applies a permutation onto an (directed) adjacency graph.
  ! The pointer and index arrays may be permuted using different permutations.
!</description>

!<input>
  ! A storage handle to the permutation that is to be applied onto the pointer
  ! array of the adjacency graph. May be ST_NOHANDLE if the pointer array
  ! is not to be permuted.
  integer, intent(in) :: h_IpermutePtr

  ! A storage handle to the permutation that is to be applied onto the index
  ! array of the adjacency graph. May be ST_NOHANDLE if the index array
  ! is not to be permuted.
  integer, intent(in) :: h_IpermuteIdx
!</input>

!<inputoutput>
  ! A storage handle to the pointer-array of the adjacency graph.
  integer, intent(inout) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  integer, intent(inout) :: h_Iidx
!</inputoutput>

!</subroutine>

  integer :: i,j,k,m, NNODE, NNZE
  integer, dimension(:), allocatable :: IptrOld, IidxOld,IinvIdx
  integer, dimension(:), pointer :: p_Iptr, p_Iidx, p_IpermutePtr,&
    p_IpermuteIdx
  
    ! First of all, get the arrays from the storage
    call storage_getbase_int(h_Iptr, p_Iptr)
    call storage_getbase_int(h_Iidx, p_Iidx)

    ! Get the length of the arrays
    NNODE = ubound(p_Iptr,1)-1
    NNZE = ubound(p_Iidx,1)

    ! Get the permutation arrays
    if(h_IpermutePtr .ne. ST_NOHANDLE) then
      call storage_getbase_int(h_IpermutePtr, p_IpermutePtr)
    else
      nullify(p_IpermutePtr)
    end if
    if(h_IpermuteIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(h_IpermuteIdx, p_IpermuteIdx)
      
      ! Build the inverse permutation
      k = ubound(p_IpermuteIdx,1)
      allocate(IinvIdx(k))
      do i = 1, k
        IinvIdx(p_IpermuteIdx(i)) = i
      end do
    else
      nullify(p_IpermuteIdx)
    end if
    
    ! What type of permutation is to be applied?
    if(associated(p_IpermutePtr)) then
    
      ! Create a backup of the arrays
      allocate(IptrOld(NNODE+1))
      allocate(IidxOld(NNZE))
      do i = 1, NNODE+1
        IptrOld(i) = p_Iptr(i)
      end do
      do i = 1, NNZE
        IidxOld(i) = p_Iidx(i)
      end do

      if(associated(p_IpermuteIdx)) then
        ! We permute both arrays.
        j = 1
        do i = 1, NNODE
        
          ! Set the pointer for this element
          p_Iptr(i) = j
          
          ! Get the index of the element
          m = p_IpermutePtr(i)
          
          ! And copy the indices
          do k = IptrOld(m), IptrOld(m+1)-1
            p_Iidx(j) = IinvIdx(IidxOld(k))
            j = j+1
          end do
        end do
        p_Iptr(NNODE+1) = j

      else
        ! We permute the pointer array.
        j = 1
        do i = 1, NNODE
        
          ! Set the pointer for this element
          p_Iptr(i) = j
          
          ! Get the index of the element
          m = p_IpermutePtr(i)
          
          ! And copy the indices
          do k = IptrOld(m), IptrOld(m+1)-1
            p_Iidx(j) = IidxOld(k)
            j = j+1
          end do
        end do
        p_Iptr(NNODE+1) = j
      
      end if
    
      ! Deallocate the backup
      deallocate(IidxOld)
      deallocate(IptrOld)

    else if(associated(p_IpermuteIdx)) then
      
      ! We permute the index array
      do i = 1, NNZE
        p_Iidx(i) = IinvIdx(p_Iidx(i))
      end do
      
    end if
    
    ! Release the temporary memory
    if(allocated(IinvIdx)) &
      deallocate(IinvIdx)
    
    ! That is it

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine adj_sortAdjacencies(h_Iptr, h_Iidx)

!<description>
  ! This routine sorts the index array such that the indices of the adjacent
  ! nodes are in ascending order.
!</description>

!<input>
  ! A storage handle to the pointer-array of the adjacency graph.
  integer, intent(in) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  integer, intent(in) :: h_Iidx
!</input>

  integer, dimension(:), pointer :: p_Iptr, p_Iidx
  integer :: i,j,k,m,NNODE
  
    ! Get the pointers from the storage
    call storage_getbase_int(h_Iptr, p_Iptr)
    call storage_getbase_int(h_Iidx, p_Iidx)
    
    NNODE = ubound(p_Iptr,1)-1
    
    ! The sorting algorithm used here is bubblesort - it might (and should)
    ! be replaced by a more efficient algorithm later...
    
    ! Go through all nodes
    do i = 1, NNODE
      do j = p_Iptr(i), p_Iptr(i+1)-1
        do k = j+1, p_Iptr(i+1)-1
          if(p_Iidx(j) .gt. p_Iidx(k)) then
            m = p_Iidx(j)
            p_Iidx(j) = p_Iidx(k)
            p_Iidx(k) = m
          end if
        end do ! k
      end do ! j
    end do ! i

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine adj_calcColouring(h_Iptr,h_Iidx,h_Icpt,h_Ipermute,h_InodeOrder)

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
  !         For all a, b in S, a != b  ->  a is not adjacent to b
  !
  ! The number of colours used is implicitly given by UBOUND(p_Icpt,1)-1.
!</description>

!<input>
  ! A storage handle to the pointer-array of the adjacency graph.
  integer, intent(in) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  integer, intent(in) :: h_Iidx
  
  ! OPTIONAL: A storage handle to an array holding the order in which the
  ! colouring algorithm should proceed through the nodes. If not given,
  ! canonical order is used.
  integer, optional, intent(in) :: h_InodeOrder
!</input>

!<output>
  ! A storage handle to the permutation that has to be applied to
  ! get a coloured partitioning.
  integer, intent(out) :: h_Ipermute

  ! A storage handle to the colour partition table.
  ! The number of colours used is implicitly given by UBOUND(Icpt,1)-1.
  integer, intent(out) :: h_Icpt
!</output>

!</subroutine>

  integer, dimension(:), pointer :: p_Iptr, p_Iidx, p_Icpt, p_Ipermute,&
    p_InodeOrder
  integer, dimension(:), allocatable :: InodeColour, IcolourMap, InumNodes
  integer :: i,j,k,inode,DEGREE,NNODE,inumColours,iminColour,iminNodes
  
    ! Get the pointers from the storage
    call storage_getbase_int(h_Iptr, p_Iptr)
    call storage_getbase_int(h_Iidx, p_Iidx)
    
    ! First of all, we need to calculate the degree of the graph, as this
    ! is the upper bound for the number of colours we will need.
    NNODE = ubound(p_Iptr,1)-1
    DEGREE = 1
    do i = 1, NNODE
      DEGREE = max(DEGREE, p_Iptr(i+1)-p_Iptr(i))
    end do
    
    ! Allocate the auxiliary arrays
    allocate(InodeColour(NNODE))
    allocate(IcolourMap(DEGREE))
    allocate(InumNodes(DEGREE))
    do i = 1, DEGREE
      InumNodes(i) = 0
    end do
    inumColours = 0
    
    ! Do we have a node order given?
    if(present(h_InodeOrder)) then
    
      ! Get the node order array from the storage then
      call storage_getbase_int(h_InodeOrder, p_InodeOrder)
      
      ! And format the node-colour array
      do i = 1, NNODE
        InodeColour(i) = 0
      end do

      ! Now let us loop through the nodes
      do i = 1, NNODE
      
        ! Get the node index
        inode = p_InodeOrder(i)
        
        ! Format the colour map
        do j = 1, inumColours
          IcolourMap(j) = 0
        end do
        
        ! Go through the adjacencies of this node
        do j = p_Iptr(inode), p_Iptr(inode+1)-1
          
          ! The node`s colour
          k = InodeColour(p_Iidx(j))
          
          if (k .gt. 0) then
            ! Mark this node`s colour as 'used'
            IcolourMap(k) = 1
          end if
          
        end do ! j
        
        ! Now find a free colour with minimum nodes
        iminColour = -1
        iminNodes = NNODE + 100
        do j = 1, inumColours
        
          ! Skip this colour if it is used
          if(IcolourMap(j) .ne. 0) cycle
          
          ! Is this the colour with minimum nodes?
          if(InumNodes(j) .lt. iminNodes) then
            iminNodes = InumNodes(j)
            iminColour = j
          end if
          
        end do ! j
        
        ! Did we find a free colour?
        if(iminColour .gt. 0) then
          ! Yes, so set the node`s colour to this one
          InodeColour(inode) = iminColour
          InumNodes(iminColour) = InumNodes(iminColour) + 1
        else
          ! No, so add a new colour
          inumColours = inumColours+1
          InodeColour(inode) = inumColours
          InumNodes(inumColours) = 1
        end if
      
      end do ! i
    
    else
    
      ! Now let us loop through the nodes
      do i = 1, NNODE
        
        ! Format the colour map
        do j = 1, inumColours
          IcolourMap(j) = 0
        end do
        
        ! Go through the adjacencies of this node
        do j = p_Iptr(i), p_Iptr(i+1)-1
          
          ! The the node index
          k = p_Iidx(j)
          
          if (k .lt. i) then
            ! Mark this node`s colour as 'used'
            IcolourMap(InodeColour(k)) = 1
          end if
          
        end do ! j
        
        ! Now find a free colour with minimum entries
        iminColour = -1
        iminNodes = NNODE + 100
        do j = 1, inumColours
        
          ! Skip this colour if it is used
          if(IcolourMap(j) .ne. 0) cycle
          
          ! Is this the colour with minimum nodes?
          if(InumNodes(j) .lt. iminNodes) then
            iminNodes = InumNodes(j)
            iminColour = j
          end if
          
        end do ! j
        
        ! Did we find a free colour?
        if(iminColour .gt. 0) then
          ! Yes, so set the node`s colour to this one
          InodeColour(i) = iminColour
          InumNodes(iminColour) = InumNodes(iminColour) + 1
        else
          ! No, so add a new colour
          inumColours = inumColours+1
          InodeColour(i) = inumColours
          InumNodes(inumColours) = 1
        end if
      
      end do ! i
    
    end if
    
    ! Now IcellColour contains a valid colouring of our adjacency graph.
    ! We now need to build the colour partition table...
    call storage_new ('adj_calcColouring', 'p_Icpt', inumColours+1, ST_INT,&
                      h_Icpt, ST_NEWBLOCK_NOINIT)
    
    ! Get the array from the storage
    call storage_getbase_int(h_Icpt, p_Icpt)

    ! Go through the colours
    p_Icpt(1) = 1
    do i = 1, inumColours
      p_Icpt(i+1) = p_Icpt(i) + InumNodes(i)
    end do
    
    ! Now make a backup of the cpt - we will need it to build the permutation
    do i = 1, inumColours
      IcolourMap(i) = p_Icpt(i)
    end do
    
    ! Build the inverse permutation array - we will abuse the InodeColour array
    ! to store the inverse permutation.
    do i = 1, NNODE
      k = InodeColour(i)
      InodeColour(i) = IcolourMap(k)
      IcolourMap(k) = IcolourMap(k)+1
    end do

    ! Allocate an array for the permuation
    call storage_new ('adj_calcColouring', 'p_Ipermute', NNODE, ST_INT,&
                      h_Ipermute, ST_NEWBLOCK_NOINIT)
    
    ! Get the array from the storage
    call storage_getbase_int(h_Ipermute, p_Ipermute)
    
    ! Calculate the permutation
    do i = 1, NNODE
      p_Ipermute(InodeColour(i)) = i
    end do
    
    ! And release all temporary memory
    deallocate(InodeColour)
    deallocate(InumNodes)
    deallocate(IcolourMap)
    
    ! That is it

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine adj_calcCuthillMcKee(h_Iptr,h_Iidx,h_Ipermute,cflags)

!<description>
  ! This routine calculates a (reverse) Cuthill-McKee ordering permutation
  ! for a given adjacency graph.
!</description>

!<input>
  ! A storage handle to the pointer-array of the adjacency graph.
  integer, intent(in) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  integer, intent(in) :: h_Iidx
  
  ! OPTIONAL: A combination ADJ_CMK_FLAG_XXXX constants defined at the top
  ! of this module specifying the flags for the algorithm. If not given,
  ! only the ADJ_CMK_FLAG_ROOT_MIN flag is set.
  integer(I32), optional, intent(in) :: cflags
!</input>

!<output>
  ! A storage handle to the (reverse) Cuthill-McKee permutation.
  integer, intent(out) :: h_Ipermute
!</output>

!</subroutine>

  integer, dimension(:), pointer :: p_Iptr, p_Iidx, p_Ipermute
  integer, dimension(:), allocatable :: Iaux
  integer :: i,j,k,n,NNODE,root,lvl1,lvl2,lvl3,lvl_root,sort_choice, root_choice
  logical :: bReverse
  
    ! Reset output handle
    h_Ipermute = ST_NOHANDLE
  
    ! Set the default strategy
    sort_choice = 0
    root_choice = ADJ_CMK_FLAG_ROOT_MIN
    bReverse = .false.
    if(present(cflags)) then
      ! Check if the flags are valid
      if(((iand(cflags,ADJ_CMK_FLAG_SORT_MIN) .ne. 0) .and. &
          (iand(cflags,ADJ_CMK_FLAG_SORT_MAX) .ne. 0)) .or. &
         ((iand(cflags,ADJ_CMK_FLAG_ROOT_MIN) .ne. 0) .and. &
          (iand(cflags,ADJ_CMK_FLAG_ROOT_MAX) .ne. 0))) then

        call output_line ('Invalid flag combination', &
                          OU_CLASS_ERROR,OU_MODE_STD,'adj_calcCuthillMcKee')
        call sys_halt()

      end if
      sort_choice = iand(cflags,ior(ADJ_CMK_FLAG_SORT_MIN,ADJ_CMK_FLAG_SORT_MAX))
      root_choice = iand(cflags,ior(ADJ_CMK_FLAG_ROOT_MIN,ADJ_CMK_FLAG_ROOT_MAX))
      bReverse = (iand(cflags,ADJ_CMK_FLAG_REVERSE) .ne. 0)
    end if
    
    ! Get the pointers from the storage
    call storage_getbase_int(h_Iptr, p_Iptr)
    call storage_getbase_int(h_Iidx, p_Iidx)
    NNODE = ubound(p_Iptr,1)-1
    
    ! Allocate the permutation array
    call storage_new ('adj_calcCuthillMcKee', 'p_Ipermute', NNODE, ST_INT, &
                      h_Ipermute, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(h_Ipermute,p_Ipermute)
    
    ! Allocate an auxiliary array
    allocate(Iaux(NNODE))
    
    ! Calculate the degree for all nodes
    do i = 1, NNODE
      Iaux(i) = p_Iptr(i+1) - p_Iptr(i)
    end do

    lvl1 = 0
    lvl2 = 0
    do while(lvl2 .lt. NNODE)
    
      ! Usually, the DO-WHILE loop we are in performs just one iteration.
      ! However, if the adjacency graph is reducible, then we cannot reach all
      ! nodes by running through the adjacency levels of only one root - so we
      ! are forced to choose a new root for every disjunct subgraph...
      root = -1
      
      ! Now how do we choose our root node?
      select case(root_choice)
      case (ADJ_CMK_FLAG_ROOT_MIN)
        ! Choose a free node with minimum degree as root
        n = NNODE+100
        do i = 1, NNODE
          if((Iaux(i) .gt. 0) .and. (Iaux(i) .lt. n)) then
            root = i
            n = Iaux(i)
          end if
        end do

      case (ADJ_CMK_FLAG_ROOT_MAX)
        ! Choose a free node with maximum degree as root
        n = 0
        do i = 1, NNODE
          if(Iaux(i) .gt. n) then
            root = i
            n = Iaux(i)
          end if
        end do
      
      case DEFAULT
        ! Choose the first free node as 'root'
        do i = 1, NNODE
          if(Iaux(i) .gt. 0) then
            root = i
            exit
          end if
        end do
        
      end select

      if(root .le. 0) then
      
        ! If we come out here, then something must have gone terribly wrong.
        ! We have not processed all nodes yet, but we cannot find a new root...
        call output_line ('internal error - no root found', &
                          OU_CLASS_ERROR,OU_MODE_STD,'adj_calcCuthillMcKee')

        call sys_halt()
        
      end if
      
      ! Initialise the first level for the root node
      lvl1 = lvl1+1
      lvl2 = lvl1
      lvl3 = lvl1
      
      ! Add the root into the permutation array
      lvl_root = lvl1
      p_Ipermute(lvl1) = root
      Iaux(root) = -Iaux(root)
    
      ! Now let us go through the adjacency levels of the root, until there
      ! are no more adjacent nodes left.
      do while(lvl2 .lt. NNODE)
      
        ! Go through all nodes in the current level, and fetch the node
        ! indices for the next level
        do i = lvl1, lvl2
        
          ! Get the index of the node
          n = p_Ipermute(i)
        
          ! Go through all nodes which are adjacent to node n
          do j = p_Iptr(n), p_Iptr(n+1)-1
          
            ! Get the index of the adjacent node
            k = p_Iidx(j)
            
            ! Has this node already been processed?
            if(Iaux(k) .gt. 0) then
            
              ! No, so let us add it to the next level
              lvl3 = lvl3+1
              p_Ipermute(lvl3) = k
              
              ! And negate it is entry in the auxiliary array to mark it
              ! as 'already processed'
              Iaux(k) = -Iaux(k)
              
            end if
            
          end do ! j
        
        end do ! i
        
        ! Did not we find any adjacent nodes anymore?
        ! Then jump out of this loop
        if(lvl3 .le. lvl2) exit

        ! Do we have to sort the permutation array for the next level?
        select case(sort_choice)
        case (ADJ_CMK_FLAG_SORT_MIN)
          call adj_cmk_aux_sort_min(lvl2+1,lvl3,p_Ipermute,Iaux)
        case (ADJ_CMK_FLAG_SORT_MAX)
          call adj_cmk_aux_sort_max(lvl2+1,lvl3,p_Ipermute,Iaux)
        end select
        
        ! Otherwise let us proceed with the next adjacency level
        lvl1 = lvl2+1
        lvl2 = lvl3
      
      end do ! WHILE(lvl2 .LT. NNODE)
      
      ! Do we have to reverse the ordering?
      if(bReverse) then
        i = lvl_root
        j = lvl2
        do while(i .lt. j)
          n = p_Ipermute(i)
          p_Ipermute(i) = p_Ipermute(j)
          p_Ipermute(j) = n
          i = i+1
          j = j-1
        end do
      end if
      
    end do ! WHILE(lvl2 .LT. NNODE)
    
    ! Deallocate the auxiliary array
    deallocate(Iaux)
    
    ! That is it

  contains
  
    pure subroutine adj_cmk_aux_sort_min(m,n,Iperm,Iaux)
    integer, intent(in) :: m,n
    integer, dimension(:), intent(inout) :: Iperm, Iaux
    
    integer :: i,j,k
    
      ! Return if there is nothing to do here...
      if(m .le. n) return
    
      ! Bubble-sort - replace with something faster later...
      
      ! Remark:
      ! Keep in mind that the Iaux array contains the negative degree of
      ! the nodes, as the nodes have already been processed! So we need
      ! to sort the nodes such that their entries in the Iaux array are
      ! in descending order instead of ascending order!
      do i = m, n
        do j = i+1, n
          if(Iaux(Iperm(j)) .gt. Iaux(Iperm(i))) then
            k = Iperm(j)
            Iperm(j) = Iperm(i)
            Iperm(i) = k
          end if
        end do
      end do
    
    end subroutine ! adj_cmk_aux_sort_min
  
    pure subroutine adj_cmk_aux_sort_max(m,n,Iperm,Iaux)
    integer, intent(in) :: m,n
    integer, dimension(:), intent(inout) :: Iperm, Iaux
    
    integer :: i,j,k
    
      ! Return if there is nothing to do here...
      if(m .le. n) return
    
      ! Bubble-sort - replace with something faster later...
      do i = m, n
        do j = i+1, n
          if(Iaux(Iperm(j)) .lt. Iaux(Iperm(i))) then
            k = Iperm(j)
            Iperm(j) = Iperm(i)
            Iperm(i) = k
          end if
        end do
      end do
    
    end subroutine ! adj_cmk_aux_sort_max
  
  end subroutine ! adj_calcCuthillMcKee
  
end module
