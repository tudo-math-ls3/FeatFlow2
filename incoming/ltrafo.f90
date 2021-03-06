
!#########################################################################
!# ***********************************************************************
!# <name> transform_basicRoutines </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains basic transformation routines for lines and
!# line segments.
!#
!# </purpose>
!#########################################################################

module transform_basicRoutines

  use fsystem

  implicit none  

  contains 

!<function>  
  !REAL(DP) FUNCTION DPO2SG(dpx,dpy,Pstart,Pend)   
  real(DP) function pdist_point_StrghtLineSeg(dpx,dpy,Pstart,Pend)     
      
  !<description>
  ! Determine minimum (Euclidean) distance from a point (dpx,dpy) to 
  ! a straight line segment. 
  !
  ! Return value = minimum distance
  !</description>

  !<input>

  ! dpx - x-coodinate of point
  real(DP), intent(IN) :: dpx 

  ! dpy - y-coodinate of point
  real(DP), intent(IN) :: dpy 

  ! Pstart - starting point of the line segment.
  !          Pstart(1) = x-coordinate
  !          Pstart(2) = y-coordinate
  real(DP), dimension(:), intent(IN) :: Pstart 

  ! Pend - endpoint of the line segment.
  !        Pend(1) = x-coordinate
  !        Pend(2) = y-coordinate
  real(DP), dimension(:), intent(IN) :: Pend

  !</input>
  
!</function>

  ! local variables
  real(DP) :: dsx1,dsy1,dsx2,dsy2,dvx,dvy,dwx,dwy
  real(DP) :: dbx,dby,db
  real(DP) :: dvl,dpl

  dsx1 = Pstart(1)
  dsy1 = Pstart(2)  
  dsx2 = Pend(1)
  dsy2 = Pend(2)

  ! Vector Pstart -> Pend and its length
  dvx = dsx2-dsx1
  dvy = dsy2-dsy1
  dvl = sqrt(dvx*dvx+dvy*dvy)

  ! Vector Pstart -> (dpx,dpy) 
  dwx = dpx-dsx1
  dwy = dpy-dsy1

  ! Scalar product to calculate the length of the
  ! projected vector
  dpl = dwx*dvx+dwy*dvy
      
  ! Relative length of the projection:
  db  = dpl/(dvl*dvl)
      
  ! Relative Length <= 0 
  ! => Connection between Pstart and (dpx,dpy) is shortest distance
  !                
  ! Relative Length >= 1
  ! => Connection between Pend and (dpx,dpy) is shortest distance
  if (db.le.0D0) then
    pdist_point_StrghtLineSeg = sqrt((dpx-dsx1)**2+(dpy-dsy1)**2)
  else if (db.ge.1D0) then
    pdist_point_StrghtLineSeg = sqrt((dpx-dsx2)**2+(dpy-dsy2)**2)
  else
    ! Calculate the projection and the distance to that point.
    ! Remember to take the square root of the distances as they 
    ! are squared...
    dbx = dsx1+dvx*db
    dby = dsy1+dvy*db
    pdist_point_StrghtLineSeg = sqrt((dpx-dbx)**2+(dpy-dby)**2)
  endif

  end function pdist_point_StrghtLineSeg



!<subroutine>  
  !SUBROUTINE PPO2SG(dpx,dpy,Pstart,Pend,DPARM,DPROJ)
  subroutine pproj_point_StrghtLineSeg(dpx,dpy,Pstart,Pend,dparm,Dproj)
      
  !<description>
  ! Project a point (dpx,dpy) onto a line segment.
  !</description>

  !<input>

  ! dpx - x-coodinate of point
  real(DP), intent(IN) :: dpx 

  ! dpy - y-coodinate of point
  real(DP), intent(IN) :: dpy 

  ! Pstart - starting point of the line segment.
  !          Pstart(1) = x-coordinate
  !          Pstart(2) = y-coordinate
  real(DP), dimension(:), intent(IN) :: Pstart 

  ! Pend - endpoint of the line segment.
  !        Pend(1) = x-coordinate
  !        Pend(2) = y-coordinate
  real(DP), dimension(:), intent(IN) :: Pend

  !</input>
  
  !<output>
  
  ! Dproj - Projection of (dpx,dpy) onto the line segment.
  !         Dproj(1) = x-coordinate
  !         Dproj(2) = y-coordinate
  real(DP), dimension(:), intent(OUT) :: Dproj

  ! dparm  - Parameter value of the projected point along the line
  !          segment. Values in [0,1].
  !          Dproj = Pstart + dparm*(Pend-Pstart)
  real(DP), intent(OUT) :: dparm  

  !</output>

! </subroutine>

  ! local variables
  real(DP) :: dsx1,dsy1,dsx2,dsy2,dvx,dvy,dwx,dwy
  real(DP) :: dvl, dpl

  dsx1 = Pstart(1)
  dsy1 = Pstart(2)
  dsx2 = Pend(1)
  dsy2 = Pend(2)

  ! Vector Pstart -> Pend and its length
  dvx = dsx2-dsx1
  dvy = dsy2-dsy1
  dvl = sqrt(dvx*dvx+dvy*dvy)

  ! Vector Pstart -> (dpx,dpy) 
  dwx = dpx-dsx1
  dwy = dpy-dsy1
      
  ! Scalar product to calculate the length of the
  ! projected vector 
  dpl = dwx*dvx+dwy*dvy

  ! Relative length of our projected vector:
  dparm  = dpl/(dvl*dvl)

  ! Length <= 0 
  ! => Point is projected onto the starting point Pstart
  !                
  ! Length >= 1
  ! => Point is projected onto the endpoint Pend
  if (dparm.le.0D0) then
    dparm = 0D0
    Dproj(1)=dsx1
    Dproj(2)=dsy1
  else if (dparm.ge.1D0) then
    dparm = 1D0
    Dproj(1)=dsx2
    Dproj(2)=dsy2
  else
    ! Calculate the projection and the distance to that point.
    ! Remember to take the square root of the distances as they 
    ! are squared...
    Dproj(1) = dsx1+dvx*dparm
    Dproj(2) = dsy1+dvy*dparm
  endif

  end subroutine pproj_point_StrghtLineSeg


!<function>  
  !REAL(DP) FUNCTION DPO2SL(dpx,dpy,PPT1,PPT2)
  real(DP) function pdist_point_StrghtLine(dpx,dpy,Ppt1,Ppt2)

  !<description>
  ! Determine minimum (Euclidean) distance from a point (dpx,dpy) to 
  ! a straight line, given by two points on the line. 
  ! 
  ! In contrast to pdist_point_StrghtLineSeg the line has no starting/ending point!
  !
  ! Return value = minimum distance, >= 0.
  ! < 0 indicates an error: Ppt1=Ppt2
  !</description>

  !<input>

  ! dpx - x-coodinate of point
  real(DP), intent(IN) :: dpx 

  ! dpy - y-coodinate of point
  real(DP), intent(IN) :: dpy 

  ! Ppt1 - one point of the line segment.
  !          Ppt1(1) = x-coordinate
  !          Ppt1(2) = y-coordinate
  real(DP), dimension(:), intent(IN) :: Ppt1 

  ! Ppt2   - another point of the line segment. 
  !          Ppt2(1) = x-coordinate
  !          Ppt2(2) = y-coordinate
  real(DP), dimension(:), intent(IN) :: Ppt2

  !</input>

!</function>
       
  ! local variables
  real(DP) :: dsx1,dsy1,dsx2,dsy2,dvx,dvy,dwx,dwy
  real(DP) :: db,dbx,dby
  real(DP) :: dvl,dpl

  dsx1 = Ppt1(1)
  dsy1 = Ppt1(2)
  dsx2 = Ppt2(1)
  dsy2 = Ppt2(2)

  ! Vector Pstart -> Pend and its length DC2
  dvx = dsx2-dsx1
  dvy = dsy2-dsy1
  dvl = sqrt(dvx*dvx+dvy*dvy)
      
  ! Stop calculation if both points are the same.
  if (dvl.eq.0D0) then
    pdist_point_StrghtLine = -1
    return
  end if

  ! Vector Pstart -> (dpx,dpy) 
   dwx = dpx-dsx1
   dwy = dpy-dsy1
      
   !Calculate length of projected vector
   dpl = dwx*dvx+dwy*dvy

   ! Calculate the projection
   db  = dpl/dvl
   dbx = dsx1+dvx*db
   dby = dsy1+dvy*db
   pdist_point_StrghtLine = sqrt((dpx-dby)**2+(dpy-dby)**2)

   end function pdist_point_StrghtLine


!<subroutine>  
  !SUBROUTINE PPO2SL(dpx,dpy,Ppt1,Ppt2,dparm,Dproj)
  subroutine pproj_point_StrghtLine(dpx,dpy,Ppt1,Ppt2,dparm,Dproj)      

  !<description>
  ! In contrast to pproj_point_StrghtLineSeg the line has no starting/ending point!
  !</description>

  !<input>

  ! dpx - x-coodinate of point
  real(DP), intent(IN) :: dpx 

  ! dpy - y-coodinate of point
  real(DP), intent(IN) :: dpy 

  ! Ppt1 - one point of the line segment.
  !          Ppt1(1) = x-coordinate
  !          Ppt1(2) = y-coordinate
  real(DP), dimension(:), intent(IN) :: Ppt1 

  ! Ppt2   - another point of the line segment. 
  !          Ppt2(1) = x-coordinate
  !          Ppt2(2) = y-coordinate
  real(DP), dimension(:), intent(IN) :: Ppt2

  !</input>
  
  !<output>
  
  ! Dproj - Projection of (dpx,dpy) onto the line segment.
  !         Dproj(1) = x-coordinate
  !         Dproj(2) = y-coordinate
  real(DP), dimension(:), intent(OUT) :: Dproj

  ! dparm  - Parameter value of the projected point along the line
  !          segment. Values in [0,1].
  !          Dproj = Pstart + dparm*(Pend-Pstart)
  real(DP), intent(OUT) :: dparm  

  !</output>

! </subroutine>     

  ! local variables
  real(DP) :: dsx1,dsy1,dsx2,dsy2,dvx,dvy,dwx,dwy
  real(DP) :: dvl,dpl

  dsx1 = Ppt1(1)
  dsy1 = Ppt1(2)
  dsx2 = Ppt2(1)
  dsy2 = Ppt2(2)

  ! Vector Pstart -> Pend and its length DC2
  dvx = dsx2-dsx1
  dvy = dsy2-dsy1
  dvl = sqrt(dvx*dvx+dvy*dvy)

  ! Stop calculation if both points are the same.
  if (dvl.eq.0D0) then
    dparm = 1D99
    return
  end if

  ! Vector Pstart -> (dpx,dpy) 
  dwx = dpx-dsx1
  dwy = dpy-dsy1
      
  ! Calculate length of projected vector
  dpl = dwx*dvx+dwy*dvy

  ! Calculate the projection and the distance to that point
  dparm  = dpl/(dvl*dvl)
  Dproj(1) = dsx1+dvx*dparm
  Dproj(2) = dsy1+dvy*dparm
      
  end subroutine pproj_point_StrghtLine


!<function> 
  !INTEGER FUNCTION PTROLN (dpx,dpy,Ppt1,Ppt2)
  integer function ptest_pointRightLeftOnLine(dpx,dpy,Ppt1,Ppt2)

  !<description>
  ! Test if a point is right-side of a line
  !
  ! This routine tests, if a given point (dpx,dpy) is right-side or
  ! left-side of a line, given by two points.
  !
  ! Return value >  0, if the point is right-side of the line Ppt1->Ppt2
  !              <  0, if the point if left-side of the line
  !              =  0, if the point is on the line or if an error
  !              occurred (Ppt1=Ppt2)
  !
  ! Remark: More specifically the return value is the squared scalar
  !         product of the normal vector of the line with the vector
  !         Ppt1->point...
  !</description>

  !<input>

  ! dpx - x-coodinate of point
  real(DP), intent(IN) :: dpx 

  ! dpy - y-coodinate of point
  real(DP), intent(IN) :: dpy 

  ! Ppt1 - one point of the line segment.
  !          Ppt1(1) = x-coordinate
  !          Ppt1(2) = y-coordinate
  real(DP), dimension(:), intent(IN) :: Ppt1 

  ! Ppt2   - another point of the line segment. 
  !          Ppt2(1) = x-coordinate
  !          Ppt2(2) = y-coordinate
  real(DP), dimension(:), intent(IN) :: Ppt2

  !</input>

  ! local variables    
  real(DP) :: dnx,dny

  ! To distinguish where the point is, relative to the line, we check
  ! the vector Ppt1->point in relation to the normal vector of the line.
  !
  ! At first build the normal vector by rotating the tangential vector
  ! by -90 degrees; this points to the right side of the line.
  dnx = (Ppt2(2)-Ppt1(2))
  dny = -(Ppt2(1)-Ppt1(1))

  ! If this normal vector points into the "same" direction like the
  ! vector Ppt1->point, the point is on the right.
  ! So we build the scalar product. This is:
  !  > 0, if the point is on the right
  !  < 0, if the point is on the left
  !  = 0, if the point is on the line
  ! and thus a candidate to serve as a return value of this function!
  ptest_pointRightLeftOnLine = dnx*(dpx-Ppt1(1)) + dny*(dpy-Ppt1(2))
      
  end function ptest_pointRightLeftOnLine 
      
end module transform_basicRoutines
