!##############################################################################
!# ****************************************************************************
!# <name> meshadjacency </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the calculation of adjacency graphs
!# for the cells of a mesh. The adjacency graph is similar to the sparsity
!# pattern of a CSR matrix (aka matrix storage technique 7/9), describing
!# which cells are adjacent to each other in respect to a specific adjacency
!# relation. The adjacency graph can be seen as a more extended version of
!# the mesh's IneighboursAtElement array, which is used e.g. for the
!# generation of a coloured partitioning of the mesh's cells.
!#
!# Currently, this module supports up to 3 adjacency types (depending on the
!# mesh's dimension):
!# 1. vertice-based element adjacency (1D/2D/3D):
!#    -> Two elements of the mesh are adjacent if they share at least one
!#       common vertice.
!#
!# 2. edge-based element adjacency (2D/3D):
!#    -> Two elements of the mesh are adjacent if they share at least one
!#       common edge.
!#
!# 3. face-based element adjacency (3D):
!#    -> Two elements of the mesh are adjacent if they share at least one
!#       common face.
!#
!# Additionally, an dimension-dependend alias called "neighbour-based element
!# adjacency" is defined, which is equivalent to "vertice-based" in 1D,
!# "edge-based" in 2D and "face-based" in 3D.
!#
!# The following routines can be found here:
!#
!# 1.) mshadj_calcElementAdjacency
!#     -> Calculates the adjacency graph for a given mesh in repsect to a
!#        specific adjacency relation.
!# </purpose>
!##############################################################################

module meshadjacency

!$use omp_lib
  use fsystem
  use genoutput
  use storage
  use basicgeometry
  use triangulation

  implicit none
  
  private
  
  public :: mshadj_calcElementAdjacency

!<constants>

!<constantblock description="mesh adjacency type identifiers">

  ! Element adjacency is calculated based on the IneighboursAtElement array
  ! from the triangulation.
  integer(I32), parameter, public :: MSHADJ_ADJ_NEIGHBOUR       = -1
  
  ! Element adjacency is calculated by vertex connectivity, i.e. two elements
  ! are adjacent to each other if they share at least one common vertice.
  integer(I32), parameter, public :: MSHADJ_ADJ_BY_VERTEX       = 1
  
  ! Element adjacency is calculated by edge connectivity, i.e. two elements
  ! are adjacent to each other if they share at least one common edge.
  ! In the case of a 1D grid, this is equal to MSHADJ_ADJ_BY_VERTEX.
  integer(I32), parameter, public :: MSHADJ_ADJ_BY_EDGE         = 2

  ! Element adjacency is calculated by face connectivity, i.e. two elements
  ! are adjacent to each other if they share at least one common face.
  ! In the case of a 1D grid, this is equal to MSHADJ_ADJ_BY_VERTEX.
  ! In the case of a 2D grid, this is equal to MSHADJ_ADJ_BY_EDGE.
  integer(I32), parameter, public :: MSHADJ_ADJ_BY_FACE         = 3

!</constantblock>

!</constants>

contains

  ! ***************************************************************************

!<subroutine>
  
  subroutine mshadj_calcElementAdjacency(rtria, h_Iptr, h_Iidx, cadjacency)

!<description>
  ! This routine calculates the element adjacency arrays for a given
  ! triangulation.
  ! The triangulation is silently assumed to be a conformal standard mesh.
!</description>

!<input>
  ! The triangulation structure for which the adjencency is to be calculated.
  type(t_triangulation), intent(in) :: rtria
  
  ! OPTIONAL: One of the MSHADJ_ADJ_XXXX identifiers which specifies how the
  ! adjacency is to be calculated. If not given, MSHADJ_ADJ_NEIGHBOUR is used.
  integer(I32), optional, intent(in) :: cadjacency
!</input>

!<output>
  ! A storage handle to the pointer-array of the adjacency graph.
  integer, intent(out) :: h_Iptr
  
  ! A storage handle to the index-array of the adjacency graph.
  integer, intent(out) :: h_Iidx
!</output>

!</subroutine>

  integer, dimension(:), pointer :: p_Iptr, p_Iidx, p_IadjIdx, p_Iadj
  integer, dimension(:,:), pointer :: p_IneighboursAtElement,&
    p_IverticesAtElement, p_IedgesAtElement, p_Iprim
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex,&
    p_IelementsAtEdgeIdx3D, p_IelementsAtEdge3D
  integer :: NEL,NNZE,NIAT,NIMT,i,j,k,iel,adj
  integer(I32) :: cadj
  logical :: bFound
  
    ! Reset the handles
    h_Iptr = ST_NOHANDLE
    h_Iidx = ST_NOHANDLE
    
    ! Get the correct adjacency type. If the given adjacency type parameter
    ! is equivalent to MSHADJ_ADJ_NEIGHBOUR, we will explicitly set it to
    ! MSHADJ_ADJ_NEIGHBOUR, as this adjacency is faster to compute than
    ! the more complex ones...
    cadj = MSHADJ_ADJ_NEIGHBOUR
    if(present(cadjacency)) then
      select case(cadjacency)
      case (MSHADJ_ADJ_NEIGHBOUR)
        cadj = MSHADJ_ADJ_NEIGHBOUR
      case (MSHADJ_ADJ_BY_VERTEX)
        if(rtria%ndim .eq. NDIM1D) then
          cadj = MSHADJ_ADJ_NEIGHBOUR
        else
          cadj = MSHADJ_ADJ_BY_VERTEX
        end if
      case (MSHADJ_ADJ_BY_EDGE)
        if((rtria%ndim .eq. NDIM1D) .or. (rtria%ndim .eq. NDIM2D)) then
          cadj = MSHADJ_ADJ_NEIGHBOUR
        else
          cadj = MSHADJ_ADJ_BY_EDGE
        end if
      case (MSHADJ_ADJ_BY_FACE)
        cadj = MSHADJ_ADJ_NEIGHBOUR
      end select
    end if
    
    ! Get all the arrays from the triangulation
    call storage_getbase_int2d(rtria%h_IneighboursAtElement,  p_IneighboursAtElement)
    call storage_getbase_int2d(rtria%h_IverticesAtElement, p_IverticesAtElement)
    if (rtria%h_IedgesAtElement .ne. ST_NOHANDLE) then
      call storage_getbase_int2d(rtria%h_IedgesAtElement, p_IedgesAtElement)
    else
      nullify(p_IedgesAtElement)
    end if
    if (rtria%h_IelementsAtVertexIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)
      call storage_getbase_int(rtria%h_IelementsAtVertex, p_IelementsAtVertex)
    else
      nullify(p_IelementsAtVertexIdx)
      nullify(p_IelementsAtVertex)
    end if
    if (rtria%h_IelementsAtEdgeIdx3D .ne. ST_NOHANDLE) then
      call storage_getbase_int(rtria%h_IelementsAtEdgeIdx3D, p_IelementsAtEdgeIdx3D)
      call storage_getbase_int(rtria%h_IelementsAtEdge3D, p_IelementsAtEdge3D)
    else
      nullify(p_IelementsAtEdgeIdx3D)
      nullify(p_IelementsAtEdge3D)
    end if

    ! Get the total number of elements from the mesh
    NEL = rtria%NEL

    ! Now comes the first interesting part - calculate the number of non-zero
    ! entries in our adjacency array...
    select case(rtria%ndim)
    case (NDIM1D)
      ! In the 1D case every element is adjacent to 3 elements - except for the
      ! two at the boundary.
      NNZE = 3*NEL - 2
    
    case (NDIM2D)
    
      ! Calculate the total number of inner edges in the mesh
      NIMT = -rtria%NMT + 3*rtria%InelOfType(TRIA_NVETRI2D)&
                        + 4*rtria%InelOfType(TRIA_NVEQUAD2D)
      
      ! What type of adjacency do we have here?
      if(cadj .eq. MSHADJ_ADJ_BY_VERTEX) then
      
        ! Element adjacency is calculated based on the vertices.
        ! This is a bit tricky, but since we assume that the mesh is conformal,
        ! we are still able to calculate the number of non-zero entries.
        NNZE = 0
        
        ! Loop through the p_IelementsAtVertexIdx array
        do j = 1, ubound(p_IelementsAtVertexIdx,1)-1
          
          ! Get the number of elements adjacent to vertice j
          k = p_IelementsAtVertexIdx(j+1) - p_IelementsAtVertexIdx(j)
          
          ! Every element adjacent to this vertice is coupled to all
          ! elements adjacent to this vertice...
          NNZE = NNZE + k*k
          
        end do
        
        ! Now every triangle in the mesh was coupled three times to itself -
        ! one time for each vertice.
        NNZE = NNZE - 2*rtria%InelOfType(TRIA_NVETRI2D)
        
        ! And every quadrilateral was coupled four times to itself.
        NNZE = NNZE - 3*rtria%InelOfType(TRIA_NVEQUAD2D)
        
        ! Now if two elements share a common edge, then we have counted the
        ! adjacency twice - one time for each vertice of the common edge.
        NNZE = NNZE - NIMT
        
      else
      
        ! Element adjacency is calculated based on the edges.
        NNZE = NEL + 2*NIMT
        
      end if
    
    case(NDIM3D)
    
      ! Calculate the total number of inner faces in the mesh
      !NIAT = -rtria%NAT + 6*rtria%InelOfType(TRIA_NVEHEXA3D)
      NIAT = -rtria%NAT + 6*NEL
      
      ! What type of adjacency do we have here?
      if(cadj .eq. MSHADJ_ADJ_BY_VERTEX) then

        call output_line ('3D vertice-based adjacency not yet implemented', &
                          OU_CLASS_ERROR,OU_MODE_STD,'mshadj_calcElementAdjacency')
        call sys_halt()
 
      else if(cadj .eq. MSHADJ_ADJ_BY_VERTEX) then

        call output_line ('3D edge-based adjacency not yet implemented', &
                          OU_CLASS_ERROR,OU_MODE_STD,'mshadj_calcElementAdjacency')
        call sys_halt()
      
      else
        
        ! Element adjacency is calculated based on the faces.
        NNZE = NEL + 2*NIAT
        
      end if
      
    end select
    
    ! Now call the storage to allocate the memory for our adjacency arrays
    call storage_new ('mshadj_calcElementAdjacency', 'p_Iptr', NEL+1, ST_INT, &
                      h_Iptr, ST_NEWBLOCK_NOINIT)
    call storage_new ('mshadj_calcElementAdjacency', 'p_Iidx', NNZE, ST_INT, &
                      h_Iidx, ST_NEWBLOCK_NOINIT)
    
    ! And get the pointers
    call storage_getbase_int(h_Iptr, p_Iptr)
    call storage_getbase_int(h_Iidx, p_Iidx)
    
    ! Do we calculate the adjacency based on the neighbours?
    if(cadj .eq. MSHADJ_ADJ_NEIGHBOUR) then
    
      ! This is the easy case - basically we only need to copy the
      ! neighbours-at-element array...
      k = 1
      do i = 1, NEL
      
        ! Set the pointer for this element
        p_Iptr(i) = k
        
        ! The element is always adjacent to itself
        p_Iidx(k) = i
        k = k+1
        
        ! Now go through the neighbourhood
        do j = 1, ubound(p_IneighboursAtElement,1)
        
          if(p_IneighboursAtElement(j,i) .gt. 0) then
            
            ! Add the neighbour to the adjacency
            p_Iidx(k) = p_IneighboursAtElement(j,i)
            k = k+1
            
          end if
          
        end do
      
      end do
      p_Iptr(NEL+1) = k
      
      ! That is it
      return
      
    end if
    
    ! Now we need to get the 3 arrays to assemble the adjacency from.
    ! These 3 arrays are:
    ! Iprim   = I****atElement
    ! IadjIdx = IelementsAt****Idx
    ! Iadj    = IelementsAt****
    ! **** must be replaced with 'Vertices' or 'Edges', depending on whether
    ! we assemble the adjacency based on vertice or edge connectivity.
    select case(cadj)
    case (MSHADJ_ADJ_BY_VERTEX)
      p_Iprim   => p_IverticesAtElement
      p_IadjIdx => p_IelementsAtVertexIdx
      p_Iadj    => p_IelementsAtVertex
    
    case (MSHADJ_ADJ_BY_EDGE)
      p_Iprim   => p_IedgesAtElement
      p_IadjIdx => p_IelementsAtEdgeIdx3D
      p_Iadj    => p_IelementsAtEdge3D
      
    end select
    
    ! Now loop through all elements in the mesh
    NNZE = 0
    p_Iptr(1) = 1
    do iel = 1, NEL
    
      ! First of all, the element is adjacent to itself
      NNZE = NNZE+1
      p_Iidx(NNZE) = iel
    
      ! Loop through all vertices/edges at that element
      do i = 1, ubound(p_Iprim,1)
      
        ! Skip this vertice/edge if it is zero (might happen for mixed
        ! triangle/quadrilateral meshes)
        if(p_Iprim(i,iel) .le. 0) cycle
        
        ! And loop through all elements adjacent to that vertice/edge
        do j = p_IadjIdx(p_Iprim(i,iel)), p_IadjIdx(p_Iprim(i,iel)+1)-1
        
          ! Get the element number
          adj = p_Iadj(j)
          
          ! We assume that we have to add the element adjacency
          bFound = .false.
          
          ! Go through the list
          do k = p_Iptr(iel), NNZE
          
            ! Is this the element we want to add the adjacency to?
            if(p_Iidx(k) .eq. adj) then
              ! Yes, so we do not need to add it anymore
              bFound = .true.
              exit
            end if
          
          end do ! k
          
          ! Do we already have this adjacency?
          if(bFound) cycle
          
          ! Okay, add the entry at the end of our list
          NNZE = NNZE+1
          p_Iidx(NNZE) = adj
          
        end do ! j
      
      end do ! i

      ! Set the offset for the next element
      p_Iptr(iel+1) = NNZE+1
    
    end do ! iel
    
    ! That is it

  end subroutine

end module
