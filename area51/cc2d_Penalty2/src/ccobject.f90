!##############################################################################
!# ****************************************************************************
!# <name> cccallback </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module contains subroutines which deal with the penalty object/s
!# 
!# 1) cc_fracLambda
!#      -> calculates a fractional penalty parameter for elements of the 
!#         discretisation which are totally inside the penalty object or
!#         only intersects with it. For all other elements, the default 
!#         value is 0.
!# 2) cc_cutoff         
!#      -> creates a list of elements which satisfy a condition regarding
!#         the placement of DOFs with respect with penalty object
!#
!#
!#
!#
!#
!#
!#
!#
!#
!##############################################################################

module ccobject

use fsystem
use storage
use geometry
use triangulation
use ccbasic
use basicgeometry
use cubature
use derivatives
use dofmapping
use element


implicit none

contains

!**********************************************************************************************************************
!
  subroutine cc_fracLambda(ielement,p_rtriangulation,p_rgeometryObject,dlambda_frac)
!
!**********************************************************************************************************************

!<input>
  integer, intent(IN) :: ielement
  type(t_triangulation), intent(IN), pointer    :: p_rtriangulation
  type(t_geometryObject),intent(IN), pointer  :: p_rgeometryObject  
!<output>
  real(DP), intent(OUT) :: dlambda_frac
  
! local variable
  real(DP) :: ElAreea, LocAreea
  real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  real(DP), dimension(:,:), allocatable :: PolyPoints
  
  integer  :: ivert,SlopeType,iin,icount,i,j
  real(DP) :: dxmax,dymax,dxmin,dymin,dxs,dys,dxe,dye,dxi1,dxi2,dyi1,dyi2, &
              dxcenter,dycenter,dradius,ddist,dslope,da,db,dc,ddiscr,dtol, &
              dmin,dmax
  real(DP), dimension(2,5) :: p_points

!**********************************************************************************************************************
  ! Tolerance 
   dtol = 10e-7
  ! Get the triangulation array for the point coordinates
                                 
  call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,p_DvertexCoordinates)
  call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

  do i = 1,5
    j = i
    if (i .eq. 5) j = 1
    p_points(1,i) = p_DvertexCoordinates(1,p_IverticesAtElement(j,ielement))
    p_points(2,i) = p_DvertexCoordinates(2,p_IverticesAtElement(j,ielement))
  end do

! Get the parameters of the object (circle case)
  dxcenter = p_rgeometryObject%rcoord2d%dorigin(1)
  dycenter = p_rgeometryObject%rcoord2d%dorigin(2)
  dradius = p_rgeometryObject%rcircle%dradius
! Counter for polygon points
  icount = 0
! Initialise the array that saves the local polygon points
  allocate(PolyPoints(2,16))   
  PolyPoints(:,:) = 0.0    

  ! Calculate total area of the real element with a general formula
  !  ElAreea = 1/2*sum(x_i*y_j-y_i*x_j), j = i+1, i = 1:5
  ElAreea = 0.0_dp
  do i = 1,4
    ElAreea = ElAreea + p_points(1,i)*p_points(2,i+1)-p_points(2,i)*p_points(1,i+1)
  end do
  ElAreea = 0.5_dp*ElAreea

  do ivert = 1,p_rtriangulation%NNEE
    ! Pick up the values for start point and end point of an edge of the real element
    dxs = p_DvertexCoordinates(1,p_IverticesAtElement(ivert,ielement))
    dys = p_DvertexCoordinates(2,p_IverticesAtElement(ivert,ielement))
    dxe = p_DvertexCoordinates(1,p_IverticesAtElement(mod(ivert,4)+1,ielement))
    dye = p_DvertexCoordinates(2,p_IverticesAtElement(mod(ivert,4)+1,ielement))
    dmin = min(dxs,dxe) - dtol
    dmax = max(dxs,dxe) + dtol

    ! Check wheter the start point of the edge is already in the object. In affirmative case,
    ! save the point into the local polygon array and count the point. At the end we will have
    ! a counter of how many points the polygon has together with the coordinates of these points
    call geom_isInGeometry (p_rgeometryObject,(/dxs,dys/), iin)
          
    if (iin.eq.1) then
        icount=icount+1 ! add one point to the polygon
        PolyPoints(1,icount)=dxs ! x coordinate of the point
        PolyPoints(2,icount)=dys ! y coordinate of the point
    end if
              
  ! Check wheter the investigated edge is cutted by the object between start and end points
  ! We can have 3 cases: vertical edge, horizontal edge and oblic. First 2 cases are treated
  ! in here for a cartezian type mesh.
    SlopeType=0
    if (abs(dxe-dxs) .lt. dtol) then
        SlopeType=1 ! vertical edge case
    else if (abs(dye-dys) .lt. dtol) then
        SlopeType=2 ! horizontal edge case
    else
        SlopeType=3 ! oblic edge case                 
    end if
             
    select case (SlopeType)
    case (1) ! vertical edge
      ! First calculate the distance between the start point and center
      ! point of the object. if this distance is smaller then the radius
      ! then we may calculate conection points. otherwise, there is no 
      ! conection possible between the actual edge and object
      ddist=max(dxs,dxcenter)-min(dxs,dxcenter)
      if (ddist .le. dradius) then
          dyi1=dycenter-sqrt(dradius**2-(dxs-dycenter)**2)
          dyi2=dycenter+sqrt(dradius**2-(dxs-dycenter)**2)
        if ((dyi1 .lt. dmax) .and. (dyi1 .gt. dmin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxs        
             PolyPoints(2,icount)=dyi1       
        end if
        if ((dyi2 .le. dmax) .and. (dyi2 .ge. dmin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxs        
             PolyPoints(2,icount)=dyi2        
        end if
      end if
               
    case (2) ! horizontal edge
      ddist=max(dys,dycenter)-min(dys,dycenter)
      if (ddist.le.dradius) then
          dxi1=dxcenter-sqrt(dradius**2-(dys-dycenter)**2)
          dxi2=dxcenter+sqrt(dradius**2-(dys-dycenter)**2)
        if ((dxi1.le.dmax).and.(dxi1.ge.dmin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxi1        
             PolyPoints(2,icount)=dys        
        end if
        if ((dxi2.le.dmax).and.(dxi2.ge.dmin)) then
             icount=icount+1
             PolyPoints(1,icount)=dxi2        
             PolyPoints(2,icount)=dys        
        end if
      end if
               
    case (3) ! oblic edge
         dslope = (dye-dys)/(dxe-dxs)
         da = 1+dslope**2
         db = -2*dxcenter-2*dslope**2*dxs+2*dslope*(dys-dycenter)
         dc = dxcenter**2+dslope**2*dxs**2-2*dslope*dxs*(dys-dycenter)+(dys-dycenter)**2-dradius**2
         ddiscr = db**2-4*da*dc
         if (ddiscr .gt. 0.0_dp) then
           dxi1 = (-db + sqrt(ddiscr))/(2*da)
           dyi1 = dslope*(dxi1-dxs) + dys
           dxi2 = (-db - sqrt(ddiscr))/(2*da)
           dyi2 = dslope*(dxi2-dxs) + dys
           if (sqrt((dxi1-dxi2)**2+(dyi1-dyi2)**2) .ge. sqrt((dxs-dxe)**2+(dys-dye)**2)) then
             if ((dxi1 .gt. dmin) .and. (dxi1 .lt. dmax)) then
               icount=icount+1
               PolyPoints(1,icount)=dxi1        
               PolyPoints(2,icount)=dyi1
             end if
             if ((dxi2 .gt. dmin).and.(dxi2 .lt. dmax)) then
               icount=icount+1
               PolyPoints(1,icount)=dxi2        
               PolyPoints(2,icount)=dslope*(dxi2-dxs)+dys        
             end if
           end if
         end if
      
    end select
  end do    

! Calculate the local area determined by the saved polygon points 
  LocAreea = 0.0_DP
  do i = 1,icount
     LocAreea = LocAreea + PolyPoints(1,i)*PolyPoints(2,mod(i,icount)+1)-PolyPoints(2,i)*PolyPoints(1,mod(i,icount)+1)
  end do
  LocAreea = 0.5*LocAreea

  dlambda_frac = LocAreea/ElAreea
  
  deallocate (PolyPoints)              

  end subroutine

!--------------------------------------------------------------------------------------------------
!
    subroutine cc_cutoff(rgeometryObject,rmatrix,h_IelementList,listsize,nodes_in)   
!
!--------------------------------------------------------------------------------------------------
!
    type(t_geometryObject), pointer, intent(in)  ::  rgeometryObject    
    type(t_matrixscalar), intent(IN), target :: rmatrix 
    integer, intent(in), optional :: nodes_in
    integer, intent (out) :: h_IelementList
    integer, intent (out) :: listsize

 
! <local variables>
    integer :: i,j,icount,iin,iel,ivertices
    real(dp), dimension(:,:), pointer :: p_DvertexCoordinates   
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_Ielements
    real(dp), dimension(2) :: dcoords
! </local>    

    h_IelementList = ST_NOHANDLE
    
    call storage_new('cc_cutoff', 'h_IelementList', rmatrix%p_rspatialdiscrTrial%p_rtriangulation%NEL, &
                     ST_INT, h_IelementList,ST_NEWBLOCK_NOINIT)
    call storage_getbase_int (h_IelementList, p_Ielements)
        
    listsize = 0
    iel = rmatrix%p_rspatialDiscrTrial%p_rtriangulation%NEL
    do i=1,iel
      call storage_getbase_double2d (rmatrix%p_rspatialDiscrTrial%p_rtriangulation%h_DvertexCoords, &
                                     p_DvertexCoordinates)
      call storage_getbase_int2d (rmatrix%p_rspatialDiscrTrial%p_rtriangulation%h_IverticesAtElement,&
                                 p_IverticesAtElement)
    ! Check how many vertices of the element are inside the particle 
      icount = 0
      ivertices = rmatrix%p_rspatialDiscrTrial%p_rtriangulation%NNEE
      do j=1,ivertices
        dcoords(1) = p_DvertexCoordinates(1,p_IverticesAtElement(j,i))
        dcoords(2) = p_DvertexCoordinates(2,p_IverticesAtElement(j,i))                           

        call geom_isInGeometry (rgeometryObject,dcoords, iin)
        if (iin .eq. 1) then
          icount=icount+1
        end if    
      end do
      
      select case (nodes_in)
      
      case (1) ! 0 or 4 nodes inside for normal cubature formula
        if ((icount.eq.0).or.(icount .eq. 4))then
          p_Ielements(listsize+1)=i
          listsize = listsize+1
        end if
        
      case default ! 1-3 nodes inside for adaptive cubature formula
        if ((icount.gt.0).and.(icount .lt. 4))then
          p_Ielements(listsize+1)=i
          listsize = listsize+1
        end if  
        
      end select  
    end do

    if (listsize .gt. 0) then 
      call storage_realloc ("cc_cutoff", listsize, h_IelementList,ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_int(h_IelementList,p_Ielements)
    end if
  end subroutine

end module