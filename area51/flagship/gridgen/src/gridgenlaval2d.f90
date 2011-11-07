program gridgenlaval2d
  implicit none

  !*****************************************************************************
  ! Grid generator for SFB708 Laval nozzle based on external triangle mesher
  !
  ! Configuration
  !
  !      w1
  !    +----+
  !    |    |h1
  !    |    | w2
  !    |    +----+                                           w3
  !    |         |h2                                   +------------+
  !    |         |                                     |            |
  !    +         + - - -+ (x1,y1)                      |            |
  !    :         +  ang1                               |            |
  !    * h0       +     | r1            ---------------+            |h3
  !    : r0        +          ----------                            |
  !    +            + - +----- ang2 - - - - - - - - - -+            |
  !    |                                                            |
  !    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +
  !    |<------------------- length ------------------>|
  !
  !    h0,r0   : height of the midpoint and radius of the inlet
  !    w1,w2   : widths of the inlet chamber
  !    h1,h2   : heights of the inlet chamber
  !    length  : total length of the nozzle (without outlet chamber)
  !    (x1,y1) : position of the circle midpoints
  !    r1,ang1 : radius and angle of the circle segment
  !    ang2    : angle of the diverging chamber
  !    w3,h3   : width and height of the outlet chamber
  !
  ! The configuration is mirrored along the x-axis and starts at #.
  !
  !*****************************************************************************

  ! Definition of double precision
  integer, parameter  :: DP = selected_real_kind(15,307)

  ! Definition of mathematical constant PI
  real(DP) :: pi

  ! Predefined values
  real(DP) :: dheight0 =  7.0_DP
  real(DP) :: dheight1 =  1.0_DP
  real(DP) :: dheight2 =  1.0_DP
  real(DP) :: dheight3 = 14.0_DP
  real(DP) :: dwidth1  =  1.0_DP
  real(DP) :: dwidth2  =  1.0_DP
  real(DP) :: dwidth3  = 14.0_DP
  real(DP) :: dangle1  = 85.0_DP
  real(DP) :: dangle2  =  2.0_DP
  real(DP) :: dradius0 =  2.0_DP
  real(DP) :: dradius1 = 11.75_DP
  real(DP) :: dx1      =  0.0_DP
  real(DP) :: dy1      = 13.0_DP
  real(DP) :: dlength  = 58.85_DP

  ! Coordinates
  real(DP), dimension(2,18) :: Dpoints

  ! Definition of output filename
  character(len=80) :: cfilename = "grid"

  ! Definition of output format
  character(LEN=*), parameter :: cFormat = '(F32.15)'
 
  ! local variables
  character(len=80) :: cbuffer
  character(len=32) :: cbuffer1,cbuffer2,cbuffer3,cbuffer4,cbuffer5,cbuffer6
  integer :: i
  
  !-----------------------------------------------------------------------------
  ! Initialize mathematical constant(s)
  !-----------------------------------------------------------------------------

  pi = asin(1.0_DP)*2.0_DP

  !-----------------------------------------------------------------------------
  ! Get command line arguments
  !-----------------------------------------------------------------------------
  
  do i = 1, command_argument_count()

    call get_command_argument(i, cbuffer)
    select case(cbuffer)
    case('-h','-H','--help')
      write(*,'(A)') 'Usage: gridgenlaval2d [OPTION]'
      write(*,'(A)') 'Generate grid for laval nozzle in TRI/PRM format.'
      write(*,*)
      write(*,'(A,T30,A)') '-h,  -H,  --help','this help screen'
      write(*,'(A,T30,A)') '-a1,  --angle1' ,'angle of the circle segment'
      write(*,'(A,T30,A)') '-a2,  --angle2' ,'angle of the diverging chamber'
      write(*,'(A,T30,A)') '-h0,  --height0','height of the midpoint of the inlet'
      write(*,'(A,T30,A)') '-h1,  --height1','height of the inlet chamber'
      write(*,'(A,T30,A)') '-h2,  --height2','height of the inlet chamber'
      write(*,'(A,T30,A)') '-h3,  --height3','height of the outlet chamber'
      write(*,'(A,T30,A)') '-r0,  --radius0','radius of the inlet'
      write(*,'(A,T30,A)') '-r1,  --radius1','radius of the circle segment'
      write(*,'(A,T30,A)') '-x1'            ,'horizontal position of the circle segment'
      write(*,'(A,T30,A)') '-y1'            ,'vertical position of the circle segment'
      write(*,'(A,T30,A)') '-w1,  --width1' ,'width of the inlet chamber'
      write(*,'(A,T30,A)') '-w2,  --width2' ,'width of the inlet chamber'
      write(*,'(A,T30,A)') '-w3,  --width3' ,'width of the outlet chamber'

      write(*,'(A,T30,A)') '-f,  --filename','name of the output file'
      write(*,'(T30,A)')   'If not given, the generated grid is stored in file "grid.tri/prm".'
      write(*,*)
      write(*,'(A)') 'Report bugs to <matthias.moeller@math.tu-dortmund.de>.'
      stop

    case('-a1','--angle1')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dangle1

    case('-a2','--angle2')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dangle2
      
    case('-h0','--height0')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dheight0

    case('-h1','--height1')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dheight1

    case('-h2','--height2')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dheight2

    case('-h3','--height3')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dheight3

    case('-r0','--radius0')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dradius0

    case('-r1','--radius1')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dradius1

    case('-x1')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dx1      

    case('-y1')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dy1     

    case('-w1','--width1')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dwidth1

    case('-w2','--width2')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dwidth2

    case('-w3','--width3')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dwidth3

    case('-f','--filename')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) cfilename
    end select
  end do
  
  dangle1 = pi/180._DP*dangle1
  dangle2 = pi/180._DP*dangle2

  !-----------------------------------------------------------------------------
  ! Write statistics
  !-----------------------------------------------------------------------------

  write(*,'(A)') 'Generating grid'
  write(*,'(A)') '---------------'
  write(*,'(A,T45,A)') 'name of output file: ',trim(adjustl(cfilename))


  !-----------------------------------------------------------------------------
  ! Generate coordinates of all corners
  !-----------------------------------------------------------------------------

  Dpoints(1,1) =  dx1-dradius1-dwidth1-dwidth2
  Dpoints(2,1) = -(dy1+dheight1+dheight2)

  Dpoints(1,2) =  dx1-dradius1-dwidth2
  Dpoints(2,2) =  Dpoints(2,1)

  Dpoints(1,3) =  Dpoints(1,2)
  Dpoints(2,3) = -(dy1+dheight2)

  Dpoints(1,4) =  dx1-dradius1
  Dpoints(2,4) =  Dpoints(2,3)

  Dpoints(1,5) =  Dpoints(1,4)
  Dpoints(2,5) = -dy1

  Dpoints(1,6) =  dx1+dradius1*cos(pi-dangle1)
  Dpoints(2,6) = -dy1+dradius1*sin(pi-dangle1)

  Dpoints(1,7) =  dx1+(dlength-dradius1-dwidth1-dwidth2)
  Dpoints(2,7) =  Dpoints(2,6)-tan(dangle2)*(Dpoints(1,7)-Dpoints(1,6))

  Dpoints(1,8) =  Dpoints(1,7)
  Dpoints(2,8) = -dheight3

  Dpoints(1,9) = Dpoints(1,8)+dwidth3
  Dpoints(2,9) = Dpoints(2,8)

  Dpoints(1,10:18) =  Dpoints(1,9:1:-1)
  Dpoints(2,10:18) = -Dpoints(2,9:1:-1)

  !-----------------------------------------------------------------------------
  ! Generate PRM-file
  !-----------------------------------------------------------------------------

  open(unit=100, file=trim(adjustl(cfilename))//'.prm')

  write(100,'(A)') 'NBCT'
  write(100,'(A)') '1'
  write(100,'(A)') 'IBCT'
  write(100,'(A)') '1 '
  write(100,'(A)') 'NCOMP'
  write(cbuffer1, fmt='(I4)') 6 + merge(2,0,dheight1 > 0.0_DP)&
                                + merge(2,0,dheight2 > 0.0_DP)&
                                + merge(2,0,dheight3 > 0.0_DP)&
                                + merge(2,0,dwidth1  > 0.0_DP)&
                                + merge(2,0,dwidth2  > 0.0_DP)&
                                + merge(2,0,dwidth3  > 0.0_DP)
  write(100,'(A)') trim(adjustl(cbuffer1))
  write(100,'(A)') 'ITYP NSPLINE NPAR'
  if (dwidth1  > 0) write(100,'(A)') '1 1 2 '
  if (dheight1 > 0) write(100,'(A)') '1 1 2 '
  if (dwidth2  > 0) write(100,'(A)') '1 1 2 '
  if (dheight2 > 0) write(100,'(A)') '1 1 2 '
  write(100,'(A)') '2 1 3 '   ! lower circle segment

  write(100,'(A)') '1 1 2 '   ! lower converging chamber
  if (dwidth3  > 0) then
    write(100,'(A)') '1 1 2 '
    write(100,'(A)') '1 1 2 ' ! lower boundary outlet chamber
    write(100,'(A)') '1 1 2 ' ! right boundary outlet chamber
    write(100,'(A)') '1 1 2 ' ! upper boundary outlet chamber
    write(100,'(A)') '1 1 2 '
  else
    write(100,'(A)') '1 1 2 ' ! outlet boundary without chamber
  end if
  write(100,'(A)') '1 1 2 '   ! lower converging chamber
  write(100,'(A)') '2 1 3 '   ! lower circle segment
  if (dheight2 > 0) write(100,'(A)') '1 1 2 '
  if (dwidth2  > 0) write(100,'(A)') '1 1 2 '
  if (dheight1 > 0) write(100,'(A)') '1 1 2 '
  if (dwidth1  > 0) write(100,'(A)') '1 1 2 '
  write(100,'(A)') '1 1 2 '   ! left boundary inlet chamber

  write(100,'(A)') 'PARAMETERS'
  ! Line segment
  if (dwidth1 > 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,1)
    write(cbuffer2, cFormat) Dpoints(2,1)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,2)-Dpoints(1,1)
    write(cbuffer2, cFormat) Dpoints(2,2)-Dpoints(2,1)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dheight1 > 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,2)
    write(cbuffer2, cFormat) Dpoints(2,2)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,3)-Dpoints(1,2)
    write(cbuffer2, cFormat) Dpoints(2,3)-Dpoints(2,2)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dwidth2 > 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,3)
    write(cbuffer2, cFormat) Dpoints(2,3)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,4)-Dpoints(1,3)
    write(cbuffer2, cFormat) Dpoints(2,4)-Dpoints(2,3)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dheight2 > 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,4)
    write(cbuffer2, cFormat) Dpoints(2,4)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,5)-Dpoints(1,4)
    write(cbuffer2, cFormat) Dpoints(2,5)-Dpoints(2,4)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Circular segment
  write(cbuffer1, cFormat)  dx1
  write(cbuffer2, cFormat) -dy1
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! midpoint
  write(cbuffer1, cFormat) dradius1
  write(cbuffer2, cFormat) 0.0_DP
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! radius1
  write(cbuffer1, cFormat) pi
  write(cbuffer2, cFormat) pi-dangle1
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! angle1

  ! Line segment
  write(cbuffer1, cFormat) Dpoints(1,6)
  write(cbuffer2, cFormat) Dpoints(2,6)
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
  write(cbuffer1, cFormat) Dpoints(1,7)-Dpoints(1,6)
  write(cbuffer2, cFormat) Dpoints(2,7)-Dpoints(2,6)
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction

  ! Outlet chamber?
  if (dheight3 > 0.0_DP) then
    ! Line segments
    write(cbuffer1, cFormat) Dpoints(1,7)
    write(cbuffer2, cFormat) Dpoints(2,7)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,8)-Dpoints(1,7)
    write(cbuffer2, cFormat) Dpoints(2,8)-Dpoints(2,7)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segments
    write(cbuffer1, cFormat) Dpoints(1,8)
    write(cbuffer2, cFormat) Dpoints(2,8)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,9)-Dpoints(1,8)
    write(cbuffer2, cFormat) Dpoints(2,9)-Dpoints(2,8)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segments
    write(cbuffer1, cFormat) Dpoints(1,9)
    write(cbuffer2, cFormat) Dpoints(2,9)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,10)-Dpoints(1,9)
    write(cbuffer2, cFormat) Dpoints(2,10)-Dpoints(2,9)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segments
    write(cbuffer1, cFormat) Dpoints(1,10)
    write(cbuffer2, cFormat) Dpoints(2,10)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,11)-Dpoints(1,10)
    write(cbuffer2, cFormat) Dpoints(2,11)-Dpoints(2,10)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segments
    write(cbuffer1, cFormat) Dpoints(1,11)
    write(cbuffer2, cFormat) Dpoints(2,11)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,12)-Dpoints(1,11)
    write(cbuffer2, cFormat) Dpoints(2,12)-Dpoints(2,11)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  else
    ! Line segments
    write(cbuffer1, cFormat) Dpoints(1,7)
    write(cbuffer2, cFormat) Dpoints(2,7)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,12)-Dpoints(1,7)
    write(cbuffer2, cFormat) Dpoints(2,12)-Dpoints(2,7)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segments
  write(cbuffer1, cFormat) Dpoints(1,12)
  write(cbuffer2, cFormat) Dpoints(2,12)
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
  write(cbuffer1, cFormat) Dpoints(1,13)-Dpoints(1,12)
  write(cbuffer2, cFormat) Dpoints(2,13)-Dpoints(2,12)
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  ! Circular segment
  write(cbuffer1, cFormat) dx1
  write(cbuffer2, cFormat) dy1
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! midpoint
  write(cbuffer1, cFormat) dradius1
  write(cbuffer2, cFormat) 0.0_DP
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! radius1
  write(cbuffer1, cFormat) pi+dangle1
  write(cbuffer2, cFormat) pi
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! angle1
  ! Line segment
  if (dheight2 > 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,14)
    write(cbuffer2, cFormat) Dpoints(2,14)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,15)-Dpoints(1,14)
    write(cbuffer2, cFormat) Dpoints(2,15)-Dpoints(2,14)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dwidth2 > 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,15)
    write(cbuffer2, cFormat) Dpoints(2,15)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,16)-Dpoints(1,15)
    write(cbuffer2, cFormat) Dpoints(2,16)-Dpoints(2,15)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dheight1 > 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,16)
    write(cbuffer2, cFormat) Dpoints(2,16)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,17)-Dpoints(1,16)
    write(cbuffer2, cFormat) Dpoints(2,17)-Dpoints(2,16)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dwidth1 > 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,17)
    write(cbuffer2, cFormat) Dpoints(2,17)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,18)-Dpoints(1,17)
    write(cbuffer2, cFormat) Dpoints(2,18)-Dpoints(2,17)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  write(cbuffer1, cFormat) Dpoints(1,18)
  write(cbuffer2, cFormat) Dpoints(2,18)
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
  write(cbuffer1, cFormat) Dpoints(1,1)-Dpoints(1,18)
  write(cbuffer2, cFormat) Dpoints(2,1)-Dpoints(2,18)
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  
  close(100)


  !-----------------------------------------------------------------------------
  ! Prepare generation of TRI-file
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Generate PSLG (planar straight line graph) for the external program triangle
  !-----------------------------------------------------------------------------
  
  open(unit=100, file=trim(adjustl(cfilename))//'.poly')
  
  ! Write file header
  write(100, '(A)') "# "//trim(adjustl(cfilename))//'.poly'
  write(100, '(A)') "#"
  write(100, '(A)') "# Poly file generated by gridgenlaval2d."
  write(100, '(A)') "#"

  !-----------------------------------------------------------------------------

  ! Write number of points
  write(100, '(A)') "# Number of points, 2D, one attribute, boundary markers"
  write(cbuffer1, '(I10)') 18
  write(100, '(A)') trim(adjustl(cbuffer1))//' 2 1 1'
      
  ! Write vertex coordinates
  do i = 1, 18
    
    ! Fill buffers
    write(cbuffer1, '(I10)') i
    write(cbuffer2, cFormat) Dpoints(1,i)
    write(cbuffer3, cFormat) Dpoints(2,i)
    write(cbuffer4, cFormat) 0.0_DP
    
    ! Write line to file
    write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                      trim(adjustl(cbuffer2))//' '//&
                      trim(adjustl(cbuffer3))//' '//&
                      trim(adjustl(cbuffer4))//' 0'
  end do

  ! Write number of segments
  write(100, '(A)') "# Number of segments, boundary markers"
  write(cbuffer1, '(I10)') 18
  write(100, '(A)') trim(adjustl(cbuffer1))//' 1'

  ! Write segments
  do i = 1, 18

    ! Fill buffers
    write(cbuffer1, '(I10)') i
    write(cbuffer2, '(I10)') i
    write(cbuffer3, '(I10)') mod(i,18)+1

    write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                      trim(adjustl(cbuffer2))//' '//&
                      trim(adjustl(cbuffer3))//' 0'
  end do

  ! Write number of holes
  write(100, '(A)') "0"
  write(100, '(A)') "0"

  close(100)

  ! Generate conforming Delaunay triangulation (-D) subject to a
  ! user-defined element area (-a) without generating additional
  ! Steiner points on the boundary (-Y) based on vertex distribution
  write(cbuffer1, cFormat) 0.11
  call system('triangle -V -Y -D -e -q32.5 '//&
      ' -a'//trim(adjustl(cbuffer1))//&
      ' -p '//trim(adjustl(cfilename))//'.poly')
  

!!$
!!$    !---------------------------------------------------------------------------
!!$    ! Read data for inner triangulation and convert it into TRI format
!!$    !---------------------------------------------------------------------------
!!$    
!!$    open(unit=100, file=trim(adjustl(cfilename))//'.1.node')
!!$    read(100, fmt=*) nvt
!!$    close(100)
!!$
!!$    open(unit=100, file=trim(adjustl(cfilename))//'.1.ele')
!!$    read(100, fmt=*) nel
!!$    close(100)
!!$    
!!$    open(unit=100, file=trim(adjustl(cfilename))//'.1.edge')
!!$    read(100, fmt=*) nmt
!!$    close(100)
!!$
!!$  else
!!$   
!!$    !---------------------------------------------------------------------------
!!$    ! Determine the number of inner layers manually
!!$    !---------------------------------------------------------------------------
!!$    
!!$    ! Determine the number of regular coarsening steps to obtain the
!!$    ! number of vertices/elements/midpoints of the coarse grid.
!!$    isegment = nsegments; irefine  = 0
!!$
!!$    inner: do
!!$      if ((mod(isegment, 2) .eq. 0) .and.&
!!$          (isegment .ge. 2*nminsegments)) then
!!$        isegment = isegment/2
!!$        irefine  = irefine+1
!!$      else
!!$        exit inner
!!$      end if
!!$    end do inner
!!$    
!!$    ! Determine number of vertices/elements in the interior circle.
!!$    ! Thus, we initialize the quantities NEL, NVT and NMT by the
!!$    ! values of the initial coarse grid in perform regular subdivision.
!!$    nel = isegment
!!$    nvt = isegment+1
!!$    nmt = 2*isegment
!!$    
!!$    do i = 1, irefine
!!$      ! Each edge produces a new vertex
!!$      nvt = nvt + nmt
!!$      
!!$      ! Each edge is subdivided into two new edges and 
!!$      ! each element produces three new edges
!!$      nmt = 2*nmt + 3*nel
!!$      
!!$      ! Each element is subdivided into four new elements
!!$      nel = 4*nel
!!$    end do
!!$
!!$    ! Update number of vertices/midpoints for circle sections
!!$    if (dazimuth .lt. 2.0_DP*pi) then
!!$      nvt = nvt+2**irefine
!!$      nmt = nmt+2**(irefine+1)
!!$    end if
!!$
!!$  end if
!!$
!!$  
!!$  !-----------------------------------------------------------------------------
!!$  ! Generate TRI-file
!!$  !-----------------------------------------------------------------------------
!!$
!!$  open(unit=100, file=trim(adjustl(cfilename))//'.tri')
!!$
!!$  write(100,'(A)') 'Coarse mesh 2D'
!!$  write(100,'(A)') 'Generated by gridgencirc2d'
!!$  
!!$  write(cbuffer1, '(I10)') nel + nsegments*nlayers
!!$  write(cbuffer2, '(I10)') nvt + nsegments*nlayers +&
!!$                           merge(0, nlayers, dazimuth .ge. 2.0_DP*pi)
!!$
!!$  write(100,'(A)') trim(adjustl(cbuffer1))//' '//&
!!$                   trim(adjustl(cbuffer2))//' 0 4 1  NEL NVT NMT NVE NBCT'
!!$
!!$  !-----------------------------------------------------------------------------
!!$  write(100,'(A)') 'DCORVG'
!!$  !-----------------------------------------------------------------------------
!!$
!!$  if (bexternalDelaunay) then
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! Copy/convert coordinates of inner region from external triangulation
!!$    !---------------------------------------------------------------------------
!!$
!!$    ! Write interior vertices in the inner layer
!!$    open(unit=200, file=trim(adjustl(cfilename))//'.1.node')
!!$    read(200, fmt=*) nvt
!!$    
!!$    ! Read node file line by line and convert coordinates
!!$    do ivt = 1, nvt
!!$      read(200, fmt=*) i, x, y, bdrPar, iaux
!!$      
!!$      if (iaux .gt. 1) then
!!$        ! Write boundary parametrization
!!$        write(cbuffer1, cFormat) bdrPar
!!$        write(cbuffer2, cFormat) 0.0_DP
!!$      else
!!$        ! Write inner vertex
!!$        write(cbuffer1, cFormat) x
!!$        write(cbuffer2, cFormat) y
!!$      end if
!!$      
!!$      ! Write line to file
!!$      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
!!$    end do
!!$    
!!$    close(200)
!!$    
!!$  else
!!$    !---------------------------------------------------------------------------
!!$    ! Create coordinates of inner region manually
!!$    !---------------------------------------------------------------------------
!!$
!!$    ! Write vertex at origin
!!$    if (dazimuth .ge. 2.0_DP*pi) then
!!$      write(cbuffer1, cFormat) 0.0_DP
!!$      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer1))
!!$    else
!!$      write(cbuffer1, cFormat) 2.0_DP
!!$      write(cbuffer2, cFormat) 0.0_DP
!!$      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
!!$    end if
!!$    
!!$    ! Write vertices in the inner layer
!!$    iaux = isegment
!!$
!!$    ! Loop over all refinement steps
!!$    do j = 1, 2**irefine
!!$      
!!$      ! Compute radius of current layer
!!$      r = j*(dinnerRadius)/(2**irefine)
!!$      
!!$      ! Loop over all segments in current layer
!!$      do i = 1, iaux
!!$        
!!$        ! Compute azimuth
!!$        phi = dazimuth0+(i-1)*dazimuth/iaux
!!$        
!!$        if ((dazimuth .lt. 2.0_DP*pi) .and. (i .eq. 1)) then
!!$          ! Fill buffers
!!$          write(cbuffer1, cFormat) 2.0_DP+r/douterRadius
!!$          write(cbuffer2, cFormat) 0.0_DP
!!$        else
!!$          ! Fill buffers
!!$          write(cbuffer1, cFormat) r * cos(phi)
!!$          write(cbuffer2, cFormat) r * sin(phi)
!!$        end if
!!$        
!!$        ! Write line to file
!!$        write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
!!$      end do
!!$
!!$      ! Since the circle is not closed we have to add one extra vertex
!!$      if (dazimuth .lt. 2.0_DP*pi) then
!!$      ! Fill buffers
!!$      write(cbuffer1, cFormat) 2.0_DP-r/douterRadius
!!$        write(cbuffer2, cFormat) 0.0_DP
!!$
!!$        ! Write line to file
!!$        write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
!!$      end if
!!$      
!!$      ! Increate the number of segments
!!$      iaux = iaux + isegment
!!$    end do
!!$    
!!$  end if
!!$  
!!$  ! Compute minimal grid spacing (isotropic/anisotropic ?)
!!$  if (danisotropy .eq. 1.0_DP) then
!!$    dr0 = (douterRadius-dinnerRadius)/nlayers
!!$  else
!!$    ! Compute local anisotropy factor
!!$    danisotropy = exp(log(danisotropy)/(nlayers-1.0_DP))
!!$    dr0 = (douterRadius-dinnerRadius)*(danisotropy-1)/(danisotropy**nlayers-1)
!!$  end if
!!$  
!!$  
!!$  do i = 1, merge(nsegments, nsegments+1, dazimuth .ge. 2.0_DP*pi)
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Write interior vertices in the outer layer
!!$    !---------------------------------------------------------------------------
!!$    
!!$    do j = 1, nlayers-1
!!$      
!!$      ! Compute radius
!!$      if (danisotropy .eq. 1.0_DP) then
!!$        r = dinnerRadius + dr0*j
!!$      else
!!$        r = dinnerRadius + dr0*(danisotropy**(j)-1)/(danisotropy-1)
!!$      end if
!!$      
!!$      ! Compute azimuth
!!$      phi = dazimuth0+dazimuth*(i-1)/nsegments
!!$      
!!$      if (dazimuth .ge. 2.0_DP*pi) then
!!$        ! Fill buffers
!!$        write(cbuffer1, cFormat) r * cos(phi)
!!$        write(cbuffer2, cFormat) r * sin(phi)
!!$        
!!$        ! Write line to file
!!$        write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
!!$      else
!!$        if (i .eq. 1) then
!!$          ! Fill buffer
!!$          write(cbuffer1, cFormat) 2.0_DP+r/douterRadius
!!$          write(cbuffer2, cFormat) 0.0
!!$        elseif (i .eq. nsegments+1) then
!!$          ! Fill buffer
!!$          write(cbuffer1, cFormat) 2.0_DP-r/douterRadius
!!$          write(cbuffer2, cFormat) 0.0
!!$        else
!!$          ! Fill buffers
!!$          write(cbuffer1, cFormat) r * cos(phi)
!!$          write(cbuffer2, cFormat) r * sin(phi)
!!$        end if
!!$
!!$        ! Write line to file
!!$        write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
!!$      end if
!!$    end do
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Write parameters for boundary vertex
!!$    !---------------------------------------------------------------------------
!!$    
!!$    r = real(i-1,DP)/real(nsegments,DP)
!!$    write(cbuffer1, cFormat) r
!!$    write(cbuffer2, cFormat) 0.0
!!$    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
!!$  end do
!!$
!!$  !-----------------------------------------------------------------------------
!!$  write(100,'(A)') 'KVERT'
!!$  !-----------------------------------------------------------------------------
!!$
!!$  if (bexternalDelaunay) then
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! Copy elements of inner region from external triangulation
!!$    !---------------------------------------------------------------------------
!!$
!!$    open(unit=200, file=trim(adjustl(cfilename))//'.1.ele')
!!$    read(200, fmt=*) nel
!!$
!!$    ! Read element file line by line and convert coordinates
!!$    do i = 1, nel
!!$      read(200, fmt=*) iaux, i1, i2, i3
!!$
!!$      write(cbuffer1, '(I10)') i1
!!$      write(cbuffer2, '(I10)') i2
!!$      write(cbuffer3, '(I10)') i3
!!$      write(cbuffer4, '(I10)') 0
!!$
!!$      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
!!$          trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
!!$    end do
!!$    
!!$    close(200)
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! Process connection layer between unstructured and structured grids
!!$    !---------------------------------------------------------------------------
!!$
!!$    do i = 1, nsegments
!!$      
!!$      ! Compute base vertex number
!!$      ivt = nlayers*(i-1)
!!$      
!!$      if (dazimuth .ge. 2.0_DP*pi) then
!!$        write(cbuffer1, '(I10)') mod(i-1,         nsegments)+1
!!$        write(cbuffer2, '(I10)') mod(ivt,         nsegments*nlayers)+nvt+1
!!$        write(cbuffer3, '(I10)') mod(ivt+nlayers, nsegments*nlayers)+nvt+1
!!$        write(cbuffer4, '(I10)') mod(i,           nsegments)+1
!!$      else
!!$        write(cbuffer1, '(I10)') mod(i-1,         (nsegments+1))+1
!!$        write(cbuffer2, '(I10)') mod(ivt,         (nsegments+1)*nlayers)+nvt+1
!!$        write(cbuffer3, '(I10)') mod(ivt+nlayers, (nsegments+1)*nlayers)+nvt+1
!!$        write(cbuffer4, '(I10)') mod(i,           (nsegments+1))+1
!!$      end if
!!$      
!!$      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
!!$          trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
!!$    end do
!!$    
!!$  else
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Create triangles of inner region manually
!!$    !---------------------------------------------------------------------------
!!$    
!!$    ! Process inner most triangles
!!$    do i = 1, isegment
!!$      
!!$      if (dazimuth .ge. 2.0_DP*pi) then
!!$        write(cbuffer1, '(I10)') 1
!!$        write(cbuffer2, '(I10)') mod(i-1, isegment)+2
!!$        write(cbuffer3, '(I10)') mod(i,   isegment)+2
!!$        write(cbuffer4, '(I10)') 0
!!$      else
!!$        write(cbuffer1, '(I10)') 1
!!$        write(cbuffer2, '(I10)') mod(i-1, isegment+1)+2
!!$        write(cbuffer3, '(I10)') mod(i,   isegment+1)+2
!!$        write(cbuffer4, '(I10)') 0
!!$      end if
!!$      
!!$      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
!!$          trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
!!$    end do
!!$    
!!$    
!!$    ! Initialize counter and process triangles in interior region
!!$    iaux = 0
!!$    
!!$    ! Loop over all layers
!!$    do j = 1, 2**irefine-1
!!$      
!!$      ! Update counter
!!$      iaux = iaux + j
!!$      
!!$      ! Loop over all segments of the inner most coarse grid
!!$      do i = 1, isegment
!!$        
!!$        ! Compute number of starting vertices in current layer and the
!!$        ! next layer which is located interior to the current one
!!$        if (dazimuth .ge. 2.0_DP*pi) then
!!$          ivt = iaux*isegment+(i-1)*j+i+1
!!$          jvt = (iaux-j)*isegment+(i-1)*(j-1)+i+1
!!$        else
!!$          ivt = iaux*isegment+i*(j+1)+1 
!!$          jvt = (iaux-j)*isegment+i*j+1
!!$        end if
!!$        
!!$        ! Loop over all edges in the current layer
!!$        do k = 1, iaux*isegment+i*(j+1)+1-ivt+&
!!$                  merge(0, j, dazimuth .ge. 2.0_DP*pi)
!!$          
!!$          write(cbuffer1, '(I10)') ivt
!!$          write(cbuffer2, '(I10)') ivt+1
!!$          write(cbuffer3, '(I10)') jvt
!!$          write(cbuffer4, '(I10)') 0
!!$          
!!$          write(100,'(A)') trim(adjustl(cbuffer1))//' '//&
!!$                           trim(adjustl(cbuffer2))//' '//&
!!$                           trim(adjustl(cbuffer3))//' '//&
!!$                           trim(adjustl(cbuffer4))
!!$
!!$          ! Increase vertex number in current layer
!!$          ivt = ivt + 1
!!$          
!!$          write(cbuffer1, '(I10)') ivt
!!$          ! Check if this is the last=first vertex in the inner layer
!!$          if (jvt+1 .eq. iaux*isegment+&
!!$                         merge(2, 2+j, dazimuth .ge. 2.0_DP*pi)) then
!!$            write(cbuffer2, '(I10)') (iaux-j)*isegment+2
!!$          else
!!$            write(cbuffer2, '(I10)') jvt+1
!!$          end if
!!$          write(cbuffer3, '(I10)') jvt
!!$          write(cbuffer4, '(I10)') 0
!!$          
!!$          write(100,'(A)') trim(adjustl(cbuffer1))//' '//&
!!$                           trim(adjustl(cbuffer2))//' '//&
!!$                           trim(adjustl(cbuffer3))//' '//&
!!$                           trim(adjustl(cbuffer4))
!!$       
!!$          ! Increase vertex number in inner layer
!!$          jvt = jvt + 1
!!$          
!!$        end do
!!$        
!!$        write(cbuffer1, '(I10)') ivt
!!$        ! Check if this is the last=first vertex in the current layer
!!$        if (ivt+1 .eq. (iaux+j+1)*isegment+2) then
!!$          write(cbuffer2, '(I10)') jvt
!!$        else
!!$          write(cbuffer2, '(I10)') ivt+1
!!$        end if
!!$        ! Check if this is the last=first vertex in the inner layer
!!$        if (jvt .eq. iaux*isegment+&
!!$                     merge(2, 2+j, dazimuth .ge. 2.0_DP*pi)) then
!!$          write(cbuffer3, '(I10)') (iaux-j)*isegment+2
!!$        else
!!$          write(cbuffer3, '(I10)') jvt
!!$        end if
!!$        write(cbuffer4, '(I10)') 0
!!$        
!!$        write(100,'(A)') trim(adjustl(cbuffer1))//' '//&
!!$                         trim(adjustl(cbuffer2))//' '//&
!!$                         trim(adjustl(cbuffer3))//' '//&
!!$                         trim(adjustl(cbuffer4))
!!$      end do
!!$    end do
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Process connection layer between unstructured and structured grids
!!$    !---------------------------------------------------------------------------
!!$    do i = 1, nsegments
!!$      
!!$      ! Compute base vertex number
!!$      ivt = nlayers*(i-1)
!!$
!!$      if (dazimuth .ge. 2.0_DP*pi) then
!!$        write(cbuffer1, '(I10)') mod(i-1,         nsegments)+nvt-nsegments+1
!!$        write(cbuffer2, '(I10)') mod(ivt,         nsegments*nlayers)+nvt+1
!!$        write(cbuffer3, '(I10)') mod(ivt+nlayers, nsegments*nlayers)+nvt+1
!!$        write(cbuffer4, '(I10)') mod(i,           nsegments)+nvt-nsegments+1
!!$      else
!!$        write(cbuffer1, '(I10)') mod(i-1,         (nsegments+1))+nvt-nsegments
!!$        write(cbuffer2, '(I10)') mod(ivt,         (nsegments+1)*nlayers)+nvt+1
!!$        write(cbuffer3, '(I10)') mod(ivt+nlayers, (nsegments+1)*nlayers)+nvt+1
!!$        write(cbuffer4, '(I10)') mod(i,           (nsegments+1))+nvt-nsegments
!!$      end if
!!$      
!!$      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
!!$          trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
!!$    end do
!!$    
!!$  end if
!!$  
!!$  !-----------------------------------------------------------------------------
!!$  ! Process outer structured grid
!!$  !-----------------------------------------------------------------------------
!!$  do i = 1, nsegments
!!$    do j = 1, nlayers-1
!!$
!!$      ! Compute base vertex number
!!$      ivt = j-1 + nlayers*(i-1)
!!$      
!!$      if (dazimuth .ge. 2.0_DP*pi) then
!!$        write(cbuffer1, '(I10)') mod(ivt,           nsegments*nlayers)+nvt+1
!!$        write(cbuffer2, '(I10)') mod(ivt+1,         nsegments*nlayers)+nvt+1
!!$        write(cbuffer3, '(I10)') mod(ivt+nlayers+1, nsegments*nlayers)+nvt+1
!!$        write(cbuffer4, '(I10)') mod(ivt+nlayers,   nsegments*nlayers)+nvt+1
!!$      else
!!$        write(cbuffer1, '(I10)') mod(ivt,           (nsegments+1)*nlayers)+nvt+1
!!$        write(cbuffer2, '(I10)') mod(ivt+1,         (nsegments+1)*nlayers)+nvt+1
!!$        write(cbuffer3, '(I10)') mod(ivt+nlayers+1, (nsegments+1)*nlayers)+nvt+1
!!$        write(cbuffer4, '(I10)') mod(ivt+nlayers,   (nsegments+1)*nlayers)+nvt+1
!!$      end if
!!$      
!!$      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
!!$                       trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
!!$    end do
!!$  end do
!!$  
!!$  !-----------------------------------------------------------------------------
!!$  write(100,'(A)') 'KNPR'
!!$  !-----------------------------------------------------------------------------
!!$  
!!$  if (dazimuth .ge. 2.0_DP*pi) then
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Full circle
!!$    !---------------------------------------------------------------------------
!!$    
!!$    ! Write nodal property for interior vertices
!!$    do i = 1, nvt
!!$      write(100,'(A)') '0'
!!$    end do
!!$    
!!$    do i = 1, nsegments
!!$      ! Write nodal property for interior vertices
!!$      do j = 1, nlayers-1
!!$        write(100,'(A)') '0'
!!$      end do
!!$      ! Write nodal property for boundary vertices
!!$      write(100,'(A)') '1'
!!$    end do
!!$
!!$  else
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! Circle segment
!!$    !---------------------------------------------------------------------------
!!$
!!$    if (bexternalDelaunay) then
!!$      
!!$      ! Write nodal property for vertices in inner region
!!$      open(unit=200, file=trim(adjustl(cfilename))//'.1.node')
!!$      read(200, fmt=*) nvt
!!$      do ivt = 1, nvt
!!$        read(200, fmt=*) i, x, y, bdrPar, iaux
!!$        write(100,'(A)') merge('1', '0', iaux .gt. 1)
!!$      end do
!!$      close(200)
!!$
!!$    else ! Create nodal properties manually
!!$
!!$      ! Write nodal property for vertex at origin
!!$      write(100,'(A)') '1'
!!$
!!$      ! Write nodal properties for vertices in the inner layer
!!$      iaux = isegment
!!$
!!$      ! Loop over all refinement steps
!!$      do j = 1, 2**irefine
!!$
!!$        ! Loop over all segments in current layer
!!$        do i = 1, iaux
!!$          
!!$          if (i .eq. 1) then
!!$            write(100,'(A)') '1'
!!$          else
!!$            write(100,'(A)') '0'
!!$          end if
!!$        end do
!!$        
!!$        ! Since the circle is not closed we have to add one extra vertex
!!$        write(100,'(A)') '1'
!!$
!!$        ! Increate the number of segments
!!$        iaux = iaux + isegment
!!$      end do
!!$    
!!$    end if
!!$
!!$    ! Write nodal properties for vertices in the outer layer
!!$    do i = 1, nsegments+1
!!$      do j = 1, nlayers-1
!!$        if ((i .eq. 1) .or. (i .eq. nsegments+1)) then
!!$          ! Write nodal property for boundary vertex
!!$          write(100,'(A)') '1'
!!$        else
!!$          ! Write nodal property for interior vertex
!!$          write(100,'(A)') '0'
!!$        end if
!!$      end do
!!$      
!!$      ! Write nodal property for boundary vertex
!!$      write(100,'(A)') '1'
!!$    end do
!!$    
!!$  end if
!!$  
!!$  !-----------------------------------------------------------------------------
!!$  write(100,'(A)') 'KMM'
!!$  !-----------------------------------------------------------------------------
!!$  
!!$  if (dazimuth .ge. 2.0_DP*pi) then
!!$    write(cbuffer1, '(I10)') nlayers+nvt
!!$    write(cbuffer2, '(I10)') nlayers*nsegments+nvt
!!$  else
!!$    write(cbuffer1, '(I10)') nlayers+nvt
!!$    write(cbuffer2, '(I10)') nlayers+nvt-1
!!$  end if
!!$
!!$  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
!!$  
!!$  close(100)

end program gridgenlaval2d
