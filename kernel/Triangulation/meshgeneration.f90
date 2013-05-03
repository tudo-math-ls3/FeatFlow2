!##############################################################################
!# ****************************************************************************
!# <name> meshgeneration </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of routines to generate standard meshes, e.g.,
!# quadratic/rectangular meshes with standard row numbering.
!#
!# The following routines can be found here:
!#
!# 1.) meshgen_rectangular2DQuadMesh
!#     -> Creates a rectantgular 2D mesh
!#
!# </purpose>
!##############################################################################

module meshgeneration

  use fsystem
  use storage
  use basicgeometry
  use triangulation

  implicit none

  private

  public :: meshgen_rectangular2DQuadMesh

contains

  ! ***************************************************************************

!<subroutine>

  subroutine meshgen_rectangular2DQuadMesh (rtriangulation,dx0,dx1,dy0,dy1,ncellsx,ncellsy)

!<description>
  ! Creates a rectangular, regular 2D mesh of the domain [x0,x1]x[y0,y1]
  ! The domain will have ncellsx cells in x-direction and ncellsy cells 
  ! in y-direction.
  !
  ! The resulting mesh will be a raw mesh and is not connected to a domain.
!</description>

!<input>
  ! Bounds of the domain
  real(DP), intent(in) :: dx0, dx1, dy0, dy1
  
  ! Number of cells in x and y
  integer, intent(in) :: ncellsx, ncellsy
!</input>

!<inputoutput>
  ! Stucture that receives the mesh.
  type(t_triangulation), intent(out) :: rtriangulation
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata2D
    integer :: ivt
    integer :: i,j
    integer, dimension(2) :: Isize
    
    ! generate the mesh topology
    call meshgen_struct2DQuadMeshTopo(rtriangulation, ncellsx, ncellsy)

    ! Initialise DvertexCoords
    Isize = (/NDIM2D,rtriangulation%NVT/)
    call storage_new ("meshgen_rectangular2DQuadMesh", "DCORVG",&
        Isize, ST_DOUBLE,rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)

    call storage_getbase_double2D(rtriangulation%h_DvertexCoords, p_Ddata2D)
    
    do i=0,ncellsy
      do j=0,ncellsx
        
        ivt = i * (ncellsx+1) + j + 1
        p_Ddata2D(1,ivt) = dx0 + real(j,DP) * (dx1-dx0) / real(ncellsx,DP)
        p_Ddata2D(2,ivt) = dy0 + real(i,DP) * (dy1-dy0) / real(ncellsy,DP)
        
      end do
    end do

    ! Create the basic boundary information
    call tria_genRawBoundary2D (rtriangulation)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine meshgen_struct2DQuadMeshTopo (rtriangulation,ncellsx,ncellsy)

!<description>
  ! Creates the topology of a structured 2D Quad mesh.
!</description>

!<input>
  ! Number of cells in x and y
  integer, intent(in) :: ncellsx, ncellsy
!</input>

!<inputoutput>
  ! Stucture that receives the mesh.
  type(t_triangulation), intent(out) :: rtriangulation
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_Idata2D
    integer, dimension(:), pointer :: p_Idata
    integer :: ivt, iel
    integer :: i,j
    integer, dimension(2) :: Isize
    
    ! Initialise basic data
    rtriangulation%ndim = NDIM2D
    rtriangulation%NEL = ncellsx * ncellsy
    rtriangulation%NVT = (ncellsx + 1) * (ncellsy + 1)
    rtriangulation%NMT = ncellsx * (ncellsy + 1)
    rtriangulation%NNVE = 4
    rtriangulation%NBCT = 1
    rtriangulation%NNEE = rtriangulation%NNVE
    
    ! Initialise KverticesAtElement
    Isize = (/rtriangulation%NNVE,rtriangulation%NEL/)
    call storage_new ("meshgen_struct2DQuadMeshTopo", "KVERT", Isize,&
        ST_INT, rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        
    call storage_getbase_int2d(rtriangulation%h_IverticesAtElement, p_Idata2D)
        
    do i=1,ncellsy
      do j=1,ncellsx
      
        iel = (i-1) * ncellsx + j
        
        ! The vertices, in counterclockwise order
        p_Idata2D(1,iel) = (i-1) * (ncellsx+1) + (j-1) + 1
        p_Idata2D(2,iel) = (i-1) * (ncellsx+1) +     j + 1
        p_Idata2D(3,iel) =     i * (ncellsx+1) +     j + 1
        p_Idata2D(4,iel) =     i * (ncellsx+1) + (j-1) + 1
      
      end do
    end do
    
    ! Number of elements per type
    rtriangulation%InelOfType(:) = 0
    rtriangulation%InelOfType(4) = rtriangulation%NEL

    ! Create InodalProperty
    call storage_new ("meshgen_struct2DQuadMeshTopo", "KNPR", &
        rtriangulation%NVT, ST_INT, &
        rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)

    call storage_getbase_int(rtriangulation%h_InodalProperty, p_Idata)
    
    ! Bottom/top row
    do i=0,ncellsx
      p_Idata(1+i) = 1
      p_Idata(1+i+ncellsy*(ncellsx+1)) = 1
    end do
    
    ! Left/right column
    do j=0,ncellsy
      p_Idata(1+j*(ncellsx+1)) = 1
      p_Idata(1+((j+1)*(ncellsx+1)-1)) = 1
    end do

  end subroutine
  
end module
