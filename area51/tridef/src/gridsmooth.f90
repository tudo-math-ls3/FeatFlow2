!##############################################################################
!# ****************************************************************************
!# <name> gridsmooth </name>
!# ****************************************************************************
!#
!# <purpose>
!# The following routines can be found here:
!#
!# 1.) gsmth_laplaceHC2d
!#     -> HC version of the laplace smoother
!# 2.) gsmth_laplace2d
!#     -> a relaxation parameter is added to soften the
!#        smoothing.
!# 3.) gsmth_umbrella2d
!#     -> implements a basic umbrella smoother
!# 4.) gsmth_laplace
!#     -> wrapper for the laplace smoother
!# 4.) gsmth_laplace
!#     -> wrapper for the laplace smoother
!# 5.) gsmth_umbrella
!#     -> wrapper for the umbrella smoother
!# </purpose>
!##############################################################################

module gridsmooth

  use triangulation
  use fsystem
  use basicgeometry
  use storage
  
  implicit none

contains

!<subroutine>
  subroutine gsmth_laplaceHC2d(rtriangulation)
  
!<description>
  !
  ! Implements a laplace smoother that uses a relaxation parameter:
  ! p_new = p + (drelax/n) * sum_{i=0}^{n-1} (q_i - p),
  ! where p: position of vertex, n: #neighbours of p, q_i: adjacent vertices
  ! of p.
  !
!</description>

!<inputoutput>
  ! triangulation.
  type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer, dimension(:,:), pointer :: p_IverticesAtEdge
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer :: iedge,ibct,isize,NMT,ivt1,ivt2,ive
  integer, dimension(:), allocatable :: Icount
  real, dimension(:,:), allocatable :: Dcoords
  integer, dimension(:), pointer   :: p_Iproperty
  real :: dx1,dx2,dy1,dy2,dn,drel1
  real :: drelax
  
  drelax = 0.5_dp
  
  call storage_getbase_int2D (rtriangulation%h_IverticesAtEdge,&
      p_IverticesAtEdge)
      
  call storage_getbase_double2D (rtriangulation%h_DvertexCoords,&
      p_DvertexCoords)
      
  call storage_getbase_int(&
      rtriangulation%h_InodalProperty,p_Iproperty)
      
  allocate(Icount(rtriangulation%NVT))
  
  allocate(Dcoords(2,rtriangulation%NVT))
  
  Icount = 0
  Dcoords = 0.0_dp
  
  ! loop over all edges
  do iedge=1,rtriangulation%NMT
    
    ! get the adjacent vertices
    ivt1 = p_IverticesAtEdge(1,iedge)
    ivt2 = p_IverticesAtEdge(2,iedge)
    
    ! get the coordinates
    dx1 = p_DvertexCoords(1,ivt1)
    dy1 = p_DvertexCoords(2,ivt1)
    
    dx2 = p_DvertexCoords(1,ivt2)
    dy2 = p_DvertexCoords(2,ivt2)
    
    if(p_Iproperty(ivt1) .eq. 0)then
      Dcoords(1,ivt1) = Dcoords(1,ivt1) + &
                        (p_DvertexCoords(1,ivt2) - p_DvertexCoords(1,ivt1))
      Dcoords(2,ivt1) = Dcoords(2,ivt1) + &
                        (p_DvertexCoords(2,ivt2) - p_DvertexCoords(2,ivt1))
      Icount(ivt1) = Icount(ivt1) + 1
    end if
    
    if(p_Iproperty(ivt2) .eq. 0)then
      Dcoords(1,ivt2) = Dcoords(1,ivt2) + &
                        (p_DvertexCoords(1,ivt1) - p_DvertexCoords(1,ivt2))
      Dcoords(2,ivt2) = Dcoords(2,ivt2) + &
                        (p_DvertexCoords(2,ivt1) - p_DvertexCoords(2,ivt2))
      Icount(ivt2) = Icount(ivt2) + 1
    end if
    
  end do ! end iedge
  
  ! loop over the vertices and correct
  do ive=1,rtriangulation%NVT
    if(p_Iproperty(ive) .ne. 0)cycle
    dn = real(Icount(ive))
    Dcoords(1,ive) = Dcoords(1,ive) * drelax/dn
    Dcoords(2,ive) = Dcoords(2,ive) * drelax/dn
  end do
  
  
  ! loop and copy
  do ive=1,rtriangulation%NVT
    if(p_Iproperty(ive) .ne. 0)cycle
    p_DvertexCoords(1,ive) = p_DvertexCoords(1,ive) + Dcoords(1,ive)
    p_DvertexCoords(2,ive) = p_DvertexCoords(2,ive) + Dcoords(2,ive)
  end do
  
  end subroutine
 
 ! ****************************************************************************************
  
!<subroutine>
  subroutine gsmth_laplace2d(rtriangulation)
  
!<description>
  !
  ! Implements a laplace smoother that uses a relaxation parameter:
  ! p_new = p + (drelax/n) * sum_{i=0}^{n-1} (q_i - p),
  ! where p: position of vertex, n: #neighbours of p, q_i: adjacent vertices
  ! of p.
  !
!</description>

!<inputoutput>
  ! triangulation.
  type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer(PREC_VERTEXIDX), dimension(:,:), pointer :: p_IverticesAtEdge
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer :: iedge,ibct,isize,NMT,ivt1,ivt2,ive
  integer, dimension(:), allocatable :: Icount
  real, dimension(:,:), allocatable :: Dcoords
  integer, dimension(:), pointer   :: p_Iproperty
  real :: dx1,dx2,dy1,dy2,dn,drel1
  real :: drelax
  
  drelax = 0.5_dp
  
  call storage_getbase_int2D (rtriangulation%h_IverticesAtEdge,&
      p_IverticesAtEdge)
      
  call storage_getbase_double2D (rtriangulation%h_DvertexCoords,&
      p_DvertexCoords)
      
  call storage_getbase_int(&
      rtriangulation%h_InodalProperty,p_Iproperty)
      
  allocate(Icount(rtriangulation%NVT))
  
  allocate(Dcoords(2,rtriangulation%NVT))
  
  Icount = 0
  Dcoords = 0.0_dp
  
  ! loop over all edges
  do iedge=1,rtriangulation%NMT
    
    ! get the adjacent vertices
    ivt1 = p_IverticesAtEdge(1,iedge)
    ivt2 = p_IverticesAtEdge(2,iedge)
    
    ! get the coordinates
    dx1 = p_DvertexCoords(1,ivt1)
    dy1 = p_DvertexCoords(2,ivt1)
    
    dx2 = p_DvertexCoords(1,ivt2)
    dy2 = p_DvertexCoords(2,ivt2)
    
    if(p_Iproperty(ivt1) .eq. 0)then
      Dcoords(1,ivt1) = Dcoords(1,ivt1) + &
                        (p_DvertexCoords(1,ivt2) - p_DvertexCoords(1,ivt1))
      Dcoords(2,ivt1) = Dcoords(2,ivt1) + &
                        (p_DvertexCoords(2,ivt2) - p_DvertexCoords(2,ivt1))
      Icount(ivt1) = Icount(ivt1) + 1
    end if
    
    if(p_Iproperty(ivt2) .eq. 0)then
      Dcoords(1,ivt2) = Dcoords(1,ivt2) + &
                        (p_DvertexCoords(1,ivt1) - p_DvertexCoords(1,ivt2))
      Dcoords(2,ivt2) = Dcoords(2,ivt2) + &
                        (p_DvertexCoords(2,ivt1) - p_DvertexCoords(2,ivt2))
      Icount(ivt2) = Icount(ivt2) + 1
    end if
    
  end do ! end iedge
  
  ! loop over the vertices and correct
  do ive=1,rtriangulation%NVT
    if(p_Iproperty(ive) .ne. 0)cycle
    dn = real(Icount(ive))
    Dcoords(1,ive) = Dcoords(1,ive) * drelax/dn
    Dcoords(2,ive) = Dcoords(2,ive) * drelax/dn
  end do
  
  
  ! loop and copy
  do ive=1,rtriangulation%NVT
    if(p_Iproperty(ive) .ne. 0)cycle
    p_DvertexCoords(1,ive) = p_DvertexCoords(1,ive) + Dcoords(1,ive)
    p_DvertexCoords(2,ive) = p_DvertexCoords(2,ive) + Dcoords(2,ive)
  end do
  
  
  end subroutine
  
 ! ****************************************************************************************
  
!<subroutine>
  subroutine gsmth_umbrella2d(rtriangulation)
  
!<description>
  !
  ! Implements a laplace smoother that uses a relaxation parameter:
  ! p_new = (1-relax) * p + (drelax/sum_{i=0}^{n-1} d_i) * sum_{i=0}^{n-1} d_i*q_i ,
  ! where p: position of vertex, n: #neighbours of p, q_i: adjacent vertices
  ! of p, d_i: length of edge(p,d_i)
  !
!</description>

!<inputoutput>
  ! triangulation.
  type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer(PREC_VERTEXIDX), dimension(:,:), pointer :: p_IverticesAtEdge
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer :: iedge,ibct,isize,NMT,ivt1,ivt2,ive
  integer, dimension(:), allocatable :: Icount
  real, dimension(:,:), allocatable :: Dcoords
  real, dimension(:), allocatable :: Dweight
  integer, dimension(:), pointer   :: p_Iproperty
  real :: dx1,dx2,dy1,dy2,dn,dw
  real :: drelax
  
  drelax = 0.5_dp
  
  call storage_getbase_int2D (rtriangulation%h_IverticesAtEdge,&
      p_IverticesAtEdge)
      
  call storage_getbase_double2D (rtriangulation%h_DvertexCoords,&
      p_DvertexCoords)
      
  call storage_getbase_int(&
      rtriangulation%h_InodalProperty,p_Iproperty)
      
  allocate(Icount(rtriangulation%NVT))
  
  allocate(Dcoords(2,rtriangulation%NVT))
  
  allocate(Dweight(rtriangulation%NVT))
  
  Icount = 0
  Dcoords = 0.0_dp
  Dweight = 0.0_dp
  
  ! loop over all edges
  do iedge=1,rtriangulation%NMT
    
    ! get the adjacent vertices
    ivt1 = p_IverticesAtEdge(1,iedge)
    ivt2 = p_IverticesAtEdge(2,iedge)
    
    ! get the coordinates
    dx1 = p_DvertexCoords(1,ivt1)
    dy1 = p_DvertexCoords(2,ivt1)
    
    dx2 = p_DvertexCoords(1,ivt2)
    dy2 = p_DvertexCoords(2,ivt2)

    dw = sqrt((p_DvertexCoords(1,ivt1)-p_DvertexCoords(1,ivt2))**2 + &
              (p_DvertexCoords(2,ivt1)-p_DvertexCoords(2,ivt2))**2)
    
    if(p_Iproperty(ivt1) .eq. 0)then
      ! add the edge weighted by its length
      Dcoords(1,ivt1) = Dcoords(1,ivt1) + dw * p_DvertexCoords(1,ivt2)
      Dcoords(2,ivt1) = Dcoords(2,ivt1) + dw * p_DvertexCoords(2,ivt2)
      Icount(ivt1) = Icount(ivt1) + 1
      Dweight(ivt1) = Dweight(ivt1) + dw
    end if
    
    if(p_Iproperty(ivt2) .eq. 0)then
      ! add the edge weighted by its length
      Dcoords(1,ivt2) = Dcoords(1,ivt2) + dw * p_DvertexCoords(1,ivt1)
      Dcoords(2,ivt2) = Dcoords(2,ivt2) + dw * p_DvertexCoords(2,ivt1)
      Icount(ivt2) = Icount(ivt2) + 1
      Dweight(ivt2) = Dweight(ivt2) + dw
    end if
    
  end do ! end iedge
  
  ! loop over the vertices and correct
  do ive=1,rtriangulation%NVT
    if(p_Iproperty(ive) .ne. 0)cycle
    dn = real(Icount(ive))
    Dcoords(1,ive) = Dcoords(1,ive) * drelax/Dweight(ive)
    Dcoords(2,ive) = Dcoords(2,ive) * drelax/Dweight(ive)
  end do
  
  dw = 1.0_dp - drelax
  ! loop and copy
  do ive=1,rtriangulation%NVT
    if(p_Iproperty(ive) .ne. 0)cycle
    p_DvertexCoords(1,ive) = dw * p_DvertexCoords(1,ive) + Dcoords(1,ive)
    p_DvertexCoords(2,ive) = dw * p_DvertexCoords(2,ive) + Dcoords(2,ive)
  end do
  
  
  end subroutine
  
 ! ****************************************************************************************
      
!<subroutine>
  subroutine gsmth_laplace(rtriangulation)
  
!<description>
  !
  ! This wrapper examines the dimension of
  ! the problem and calls the according subroutine
  !
!</description>

!<inputoutput>
  ! triangulation.
  type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>
  
  if(rtriangulation%ndim .eq. NDIM2D) then
    call gsmth_laplace2d(rtriangulation)
  else if(rtriangulation%ndim .eq. NDIM3D) then
  
  end if
  
  end subroutine
  
 ! ****************************************************************************************
  
!<subroutine>
  subroutine gsmth_umbrella(rtriangulation)
  
!<description>
  !
  ! This wrapper examines the dimension of
  ! the problem and calls the according subroutine
  !
!</description>

!<inputoutput>
  ! triangulation.
  type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>
  
  if(rtriangulation%ndim .eq. NDIM2D) then
    call gsmth_umbrella2d(rtriangulation)
  else if(rtriangulation%ndim .eq. NDIM3D) then
  
  end if
  
  end subroutine
  
end module