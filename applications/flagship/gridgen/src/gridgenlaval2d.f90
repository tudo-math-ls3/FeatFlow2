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
  real(DP) :: dangle1  = 85.0_DP
  real(DP) :: dangle2  =  2.0_DP
  real(DP) :: darea    =  0.5_DP
  real(DP) :: ddist    =  1.0_DP
  real(DP) :: dheight0 =  7.0_DP
  real(DP) :: dheight1 =  1.0_DP
  real(DP) :: dheight2 =  1.0_DP
  real(DP) :: dheight3 = 14.0_DP
  real(DP) :: dlength  = 58.85_DP
  real(DP) :: dradius0 =  2.0_DP
  real(DP) :: dradius1 = 11.75_DP
  real(DP) :: dwidth1  =  1.0_DP
  real(DP) :: dwidth2  =  1.0_DP
  real(DP) :: dwidth3  = 14.0_DP
  real(DP) :: dx1      =  0.0_DP
  real(DP) :: dy1      = 13.0_DP

  ! Total number of fixed points
  integer, parameter :: ncoords = 22

  ! Coordinates
  real(DP), dimension(2,ncoords) :: Dpoints

  ! Definition of output filename
  character(len=80) :: cfilename = "grid"

  ! Definition of output format
  character(LEN=*), parameter :: cFormat = '(F32.15)'
 
  ! local variables
  character(len=80) :: cbuffer
  character(len=32) :: cbuffer1,cbuffer2,cbuffer3,cbuffer4,cbuffer5,cbuffer6
  logical :: brewrite
  real(DP) :: domega,dsegmentlen,dx,dy,dpar
  integer :: i,i1,i2,i3,i4,ibdc,iel,ipoint,isubpoint,ivt
  integer :: nel,npoints,nsubsegments,nvt,npar
  
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
      write(*,'(A)') "      w1"
      write(*,'(A)') "    +----+"
      write(*,'(A)') "    |    |h1"
      write(*,'(A)') "    |    | w2"
      write(*,'(A)') "    |    +----+                                           w3"
      write(*,'(A)') "    |         |h2                                   +------------+"
      write(*,'(A)') "    |         |                                     |            |"
      write(*,'(A)') "    +         + - - -+ (x1,y1)                      |            |"
      write(*,'(A)') "    :         +  ang1                               |            |"
      write(*,'(A)') "    * h0       +     | r1            ---------------+            |h3"
      write(*,'(A)') "    : r0        +          ----------                            |"
      write(*,'(A)') "    +            + - +----- ang2 - - - - - - - - - -+            |"
      write(*,'(A)') "    |                                                            |"
      write(*,'(A)') "    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +"
      write(*,'(A)') "    |<------------------- length ------------------>|"
      write(*,*)
      write(*,'(A,T30,A)') '-h,  -H,  --help','this help screen'
      write(*,'(A,T30,A)') '-a1, --angle1' ,'angle of the circle segment'
      write(*,'(A,T30,A)') '-a2, --angle2' ,'angle of the diverging chamber'
      write(*,'(A,T30,A)') '-a,  --area' ,'average area of triangles'
      write(*,'(A,T30,A)') '-d,  --dist' ,'average distance between points along the boundary'
      write(*,'(A,T30,A)') '-f,  --filename','name of the output file'
      write(*,'(T30,A)')   'If not given, the generated grid is stored in file "grid.tri/prm".'
      write(*,'(A,T30,A)') '-h0, --height0','height of the midpoint of the inlet'
      write(*,'(A,T30,A)') '-h1, --height1','height of the inlet chamber'
      write(*,'(A,T30,A)') '-h2, --height2','height of the inlet chamber'
      write(*,'(A,T30,A)') '-h3, --height3','height of the outlet chamber'
      write(*,'(A,T30,A)') '-l,  --length','length of the nozzle (without outlet chamber)'
      write(*,'(A,T30,A)') '-r0, --radius0','radius of the inlet'
      write(*,'(A,T30,A)') '-r1, --radius1','radius of the circle segment'
      write(*,'(A,T30,A)') '-s,  --segment','number of segments of the circle segment'
      write(*,'(A,T30,A)') '-x1'            ,'horizontal position of the circle segment'
      write(*,'(A,T30,A)') '-y1'            ,'vertical position of the circle segment'
      write(*,'(A,T30,A)') '-w1, --width1' ,'width of the inlet chamber'
      write(*,'(A,T30,A)') '-w2, --width2' ,'width of the inlet chamber'
      write(*,'(A,T30,A)') '-w3, --width3' ,'width of the outlet chamber'
      write(*,*)
      write(*,'(A)') 'Report bugs to <matthias.moeller@math.tu-dortmund.de>.'
      stop

    case('-a','--area')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) darea

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

    case('-l','--length')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dlength

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
      
    case('-d','--dist')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) ddist
      
    end select
  end do
  
  dangle1 = pi/180.0_DP*dangle1
  dangle2 = pi/180.0_DP*dangle2

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

  Dpoints(1,19) = Dpoints(1,18)
  Dpoints(2,19) = min(Dpoints(2,18),&
                      0.5_DP*(Dpoints(2,18)+Dpoints(2,1))+dheight0+dradius0)

  Dpoints(1,20) = Dpoints(1,18)
  Dpoints(2,20) = min(Dpoints(2,18),&
                      0.5_DP*(Dpoints(2,18)+Dpoints(2,1))+dheight0-dradius0)

  if (dradius0 >= dheight0) then
    Dpoints(:,21) = Dpoints(:,20)
    Dpoints(:,22) = Dpoints(:,20)
  else
    Dpoints(1,21) = Dpoints(1,18)
    Dpoints(2,21) = max(Dpoints(2,1),&
                        0.5_DP*(Dpoints(2,18)+Dpoints(2,1))-dheight0+dradius0)

    Dpoints(1,22) = Dpoints(1,18)
    Dpoints(2,22) = max(Dpoints(2,1),&
                        0.5_DP*(Dpoints(2,18)+Dpoints(2,1))-dheight0-dradius0)
  end if

  !-----------------------------------------------------------------------------
  ! Generate PRM-file
  !-----------------------------------------------------------------------------

  open(unit=100, file=trim(adjustl(cfilename))//'.prm')

  write(100,'(A)') 'NBCT'
  write(100,'(A)') '1'
  write(100,'(A)') 'IBCT'
  write(100,'(A)') '1 '
  write(100,'(A)') 'NCOMP'
  write(cbuffer1, fmt='(I4)') 4 + merge(2,0,dheight1 /= 0.0_DP)&
                                + merge(2,0,dheight2 /= 0.0_DP)&
                                + merge(2,0,dwidth1  /= 0.0_DP)&
                                + merge(2,0,dwidth2  /= 0.0_DP)&
                                + merge(5,1,dwidth3  > 0.0_DP)&
                                + merge(5,merge(1,3,dradius0 >= dheight0),&
                                          dheight0 > 0.0_DP .and. dradius0 < dheight0)
  write(100,'(A)') trim(adjustl(cbuffer1))
  write(100,'(A)') 'ITYP NSPLINE NPAR'
  if (dwidth1  /= 0.0_DP) write(100,'(A)') '1 1 2 '
  if (dheight1 /= 0.0_DP) write(100,'(A)') '1 1 2 '
  if (dwidth2  /= 0.0_DP) write(100,'(A)') '1 1 2 '
  if (dheight2 /= 0.0_DP) write(100,'(A)') '1 1 2 '
  write(100,'(A)') '2 1 3 '   ! lower circle segment

  write(100,'(A)') '1 1 2 '   ! lower converging chamber
  if (dwidth3  > 0.0_DP) then
    ! Note that we do not check the special case that
    ! dwidth3 equals the height of the exist throat
    write(100,'(A)') '1 1 2 '
    write(100,'(A)') '1 1 2 ' ! lower boundary outlet chamber
    write(100,'(A)') '1 1 2 ' ! right boundary outlet chamber
    write(100,'(A)') '1 1 2 ' ! upper boundary outlet chamber
    write(100,'(A)') '1 1 2 '
  else
    write(100,'(A)') '1 1 2 ' ! outlet boundary without chamber
  end if
  write(100,'(A)') '1 1 2 '   ! upper converging chamber
  write(100,'(A)') '2 1 3 '   ! upper circle segment
  if (dheight2 /= 0.0_DP) write(100,'(A)') '1 1 2 '
  if (dwidth2  /= 0.0_DP) write(100,'(A)') '1 1 2 '
  if (dheight1 /= 0.0_DP) write(100,'(A)') '1 1 2 '
  if (dwidth1  /= 0.0_DP) write(100,'(A)') '1 1 2 '

  ! left boundary inlet chamber
  if (dheight0 > 0.0_DP .and. dradius0 < dheight0) then
    write(100,'(A)') '1 1 2 '
    write(100,'(A)') '1 1 2 '
    write(100,'(A)') '1 1 2 '
    write(100,'(A)') '1 1 2 '
    write(100,'(A)') '1 1 2 '
  else if (dradius0 >= dheight0) then
    write(100,'(A)') '1 1 2 '
  else
    write(100,'(A)') '1 1 2 '
    write(100,'(A)') '1 1 2 '
    write(100,'(A)') '1 1 2 '
  end if
  
  write(100,'(A)') 'PARAMETERS'
  ! Line segment
  if (dwidth1 /= 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,1)
    write(cbuffer2, cFormat) Dpoints(2,1)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,2)-Dpoints(1,1)
    write(cbuffer2, cFormat) Dpoints(2,2)-Dpoints(2,1)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dheight1 /= 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,2)
    write(cbuffer2, cFormat) Dpoints(2,2)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,3)-Dpoints(1,2)
    write(cbuffer2, cFormat) Dpoints(2,3)-Dpoints(2,2)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dwidth2 /= 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,3)
    write(cbuffer2, cFormat) Dpoints(2,3)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,4)-Dpoints(1,3)
    write(cbuffer2, cFormat) Dpoints(2,4)-Dpoints(2,3)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dheight2 /= 0.0_DP) then
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
  if (dwidth3 > 0.0_DP) then
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
  if (dheight2 /= 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,14)
    write(cbuffer2, cFormat) Dpoints(2,14)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,15)-Dpoints(1,14)
    write(cbuffer2, cFormat) Dpoints(2,15)-Dpoints(2,14)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dwidth2 /= 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,15)
    write(cbuffer2, cFormat) Dpoints(2,15)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,16)-Dpoints(1,15)
    write(cbuffer2, cFormat) Dpoints(2,16)-Dpoints(2,15)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dheight1 /= 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,16)
    write(cbuffer2, cFormat) Dpoints(2,16)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,17)-Dpoints(1,16)
    write(cbuffer2, cFormat) Dpoints(2,17)-Dpoints(2,16)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  ! Line segment
  if (dwidth1 /= 0.0_DP) then
    write(cbuffer1, cFormat) Dpoints(1,17)
    write(cbuffer2, cFormat) Dpoints(2,17)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,18)-Dpoints(1,17)
    write(cbuffer2, cFormat) Dpoints(2,18)-Dpoints(2,17)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  
  ! left boundary inlet chamber
  if (dheight0 > 0.0_DP .and. dradius0 < dheight0) then
    ! Line segment
    write(cbuffer1, cFormat) Dpoints(1,18)
    write(cbuffer2, cFormat) Dpoints(2,18)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,19)-Dpoints(1,18)
    write(cbuffer2, cFormat) Dpoints(2,19)-Dpoints(2,18)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segment
    write(cbuffer1, cFormat) Dpoints(1,19)
    write(cbuffer2, cFormat) Dpoints(2,19)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,20)-Dpoints(1,19)
    write(cbuffer2, cFormat) Dpoints(2,20)-Dpoints(2,19)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segment
    write(cbuffer1, cFormat) Dpoints(1,20)
    write(cbuffer2, cFormat) Dpoints(2,20)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,21)-Dpoints(1,20)
    write(cbuffer2, cFormat) Dpoints(2,21)-Dpoints(2,20)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segment
    write(cbuffer1, cFormat) Dpoints(1,21)
    write(cbuffer2, cFormat) Dpoints(2,21)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,22)-Dpoints(1,21)
    write(cbuffer2, cFormat) Dpoints(2,22)-Dpoints(2,21)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segment
    write(cbuffer1, cFormat) Dpoints(1,22)
    write(cbuffer2, cFormat) Dpoints(2,22)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,1)-Dpoints(1,22)
    write(cbuffer2, cFormat) Dpoints(2,1)-Dpoints(2,22)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  elseif (dradius0 >= dheight0) then
    ! Line segment
    write(cbuffer1, cFormat) Dpoints(1,18)
    write(cbuffer2, cFormat) Dpoints(2,18)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,1)-Dpoints(1,18)
    write(cbuffer2, cFormat) Dpoints(2,1)-Dpoints(2,18)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  else
    ! Line segment
    write(cbuffer1, cFormat) Dpoints(1,18)
    write(cbuffer2, cFormat) Dpoints(2,18)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,19)-Dpoints(1,18)
    write(cbuffer2, cFormat) Dpoints(2,19)-Dpoints(2,18)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segment
    write(cbuffer1, cFormat) Dpoints(1,19)
    write(cbuffer2, cFormat) Dpoints(2,19)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,20)-Dpoints(1,19)
    write(cbuffer2, cFormat) Dpoints(2,20)-Dpoints(2,19)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segment
    write(cbuffer1, cFormat) Dpoints(1,20)
    write(cbuffer2, cFormat) Dpoints(2,20)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) Dpoints(1,1)-Dpoints(1,20)
    write(cbuffer2, cFormat) Dpoints(2,1)-Dpoints(2,20)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
  end if
  
  close(100)

  !-----------------------------------------------------------------------------
  ! Generate PSLG (planar straight line graph) for the external program triangle
  !-----------------------------------------------------------------------------
  
  ! Initialize the number of points
  npoints = 0
  brewrite = .false.

1 open(unit=100, file=trim(adjustl(cfilename))//'.poly', form='formatted')
  
  ! Write file header
  write(100, '(A)') "# "//trim(adjustl(cfilename))//'.poly'
  write(100, '(A)') "#"
  write(100, '(A)') "# Poly file generated by gridgenlaval2d."
  write(100, '(A)') "#"

  !-----------------------------------------------------------------------------

  ! Write number of points
  write(100, '(A)') "# Number of points, 2D, one attribute, boundary markers"
  write(cbuffer1, '(I10)') npoints
  write(100, '(A)') trim(adjustl(cbuffer1))//' 2 1 1'

  ! Initialize the number of points
  npoints = 0; npar = 0

  ! Write vertex coordinates
  do ipoint = 1, ncoords

    ! Check if points coincide
    if (sqrt((Dpoints(1,ipoint)-Dpoints(1,mod(ipoint,ncoords)+1))**2 +&
             (Dpoints(2,ipoint)-Dpoints(2,mod(ipoint,ncoords)+1))**2) <= 1e-12) cycle
   
    ! Increase number of points by one
    npoints = npoints+1
 
    ! Increase number of boundary parameter
    npar = npar+1

    ! Fill buffers
    write(cbuffer1, '(I10)') npoints
    write(cbuffer2, cFormat) Dpoints(1,ipoint)
    write(cbuffer3, cFormat) Dpoints(2,ipoint)
    write(cbuffer4, cFormat) real(npar-1,DP)

    ! Write line to file
    write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                      trim(adjustl(cbuffer2))//' '//&
                      trim(adjustl(cbuffer3))//' '//&
                      trim(adjustl(cbuffer4))//' 1'

    ! Write extra points between two regular points with user-defined
    ! average distance; note that we must distinguish between points
    ! on a straight line and those on a circular segment
    select case(ipoint)
    case(1:4,6:12,14:ncoords)
      ! Compute line segment length
      dsegmentlen = sqrt((Dpoints(1,mod(ipoint,ncoords)+1)-Dpoints(1,ipoint))**2+&
                         (Dpoints(2,mod(ipoint,ncoords)+1)-Dpoints(2,ipoint))**2)

      ! Compute number of auxiliary segments
      nsubsegments = nint(dsegmentlen/ddist)

      ! Write vertex coordinates of auxiliary points
      do isubpoint = 1, nsubsegments-1
        ! Increase number of points by one
        npoints = npoints+1

        ! Compute relative position of the auxiliary point on
        ! the line segment between the points [ipoint,ipoint+1]
        domega = real(isubpoint,DP)/real(nsubsegments,DP)

        ! Fill buffer
        write(cbuffer1, '(I10)') npoints
        write(cbuffer2, cFormat) Dpoints(1,ipoint) +&
                         domega*(Dpoints(1,mod(ipoint,ncoords)+1)-Dpoints(1,ipoint))
        write(cbuffer3, cFormat) Dpoints(2,ipoint) +&
                         domega*(Dpoints(2,mod(ipoint,ncoords)+1)-Dpoints(2,ipoint))
        write(cbuffer4, cFormat) real(npar-1,DP)+domega
    
        ! Write line to file
        write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                          trim(adjustl(cbuffer2))//' '//&
                          trim(adjustl(cbuffer3))//' '//&
                          trim(adjustl(cbuffer4))//' 1'
      end do

    case(5)
      ! Compute circular segment length
      dsegmentlen = dradius1*dangle1

      ! Compute number of auxiliary segments
      nsubsegments = nint(dsegmentlen/ddist)

      ! Write vertex coordinates of auxiliary points
      do isubpoint = 1, nsubsegments-1
        ! Increase number of points by one
        npoints = npoints+1

        ! Compute relative position of the auxiliary point on
        ! the circular segment between the points ipoint and ipoint+1
        domega = real(isubpoint,DP)/real(nsubsegments,DP)
        
        ! Fill buffer
        write(cbuffer1, '(I10)')  npoints
        write(cbuffer2, cFormat)  dx1+dradius1*cos(pi-dangle1*domega)
        write(cbuffer3, cFormat) -dy1+dradius1*sin(pi-dangle1*domega)
        write(cbuffer4, cFormat)  real(npar-1,DP)+domega
    
        ! Write line to file
        write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                          trim(adjustl(cbuffer2))//' '//&
                          trim(adjustl(cbuffer3))//' '//&
                          trim(adjustl(cbuffer4))//' 1'
      end do

    case(13)
      ! Compute circular segment length
      dsegmentlen = dradius1*dangle1

      ! Compute number of auxiliary segments
      nsubsegments = nint(dsegmentlen/ddist)

      ! Write vertex coordinates of auxiliary points
      do isubpoint = nsubsegments-1, 1, -1
        ! Increase number of points by one
        npoints = npoints+1

        ! Compute relative position of the auxiliary point on
        ! the circular segment between the points ipoint and ipoint+1
        domega = real(isubpoint,DP)/real(nsubsegments,DP)
        
        ! Fill buffer
        write(cbuffer1, '(I10)')  npoints
        write(cbuffer2, cFormat)  dx1+dradius1*cos(pi-dangle1*domega)
        write(cbuffer3, cFormat)  dy1-dradius1*sin(pi-dangle1*domega)
        write(cbuffer4, cFormat)  real(npar-1,DP)+(1.0_DP-domega)
    
        ! Write line to file
        write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                          trim(adjustl(cbuffer2))//' '//&
                          trim(adjustl(cbuffer3))//' '//&
                          trim(adjustl(cbuffer4))//' 1'
      end do
    end select
  end do

  ! Write number of segments
  write(100, '(A)') "# Number of segments, boundary markers"
  write(cbuffer1, '(I10)') npoints
  write(100, '(A)') trim(adjustl(cbuffer1))//' 1'

  ! Write segments
  do ipoint = 1, npoints

    ! Fill buffers
    write(cbuffer1, '(I10)') ipoint
    write(cbuffer2, '(I10)') ipoint
    write(cbuffer3, '(I10)') mod(ipoint,npoints)+1

    write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                      trim(adjustl(cbuffer2))//' '//&
                      trim(adjustl(cbuffer3))//' 1'
  end do

  ! Write number of holes
  write(100, '(A)') "0"
  write(100, '(A)') "0"

  close(100)

  ! Now, the number of total points is known and needs to be updated
  if ((npoints /= ncoords) .and. .not.brewrite) then
    brewrite = .true.
    goto 1
  end if
  
  ! Generate conforming Delaunay triangulation (-D) subject to a
  ! user-defined element area (-a) without generating additional
  ! Steiner points on the boundary (-Y) based on vertex distribution
  write(cbuffer1, cFormat) darea
  call system('triangle -V -Y -D '//&
      ' -a'//trim(adjustl(cbuffer1))//&
      ' -p '//trim(adjustl(cfilename))//'.poly')
  
  !---------------------------------------------------------------------------
  ! Read data for inner triangulation and convert it into TRI format
  !---------------------------------------------------------------------------
  
  open(unit=100, file=trim(adjustl(cfilename))//'.1.node')
  read(100, fmt=*) nvt
  close(100)
  
  open(unit=100, file=trim(adjustl(cfilename))//'.1.ele')
  read(100, fmt=*) nel
  close(100)
  
  !-----------------------------------------------------------------------------
  ! Generate TRI-file
  !-----------------------------------------------------------------------------
  
  open(unit=100, file=trim(adjustl(cfilename))//'.tri')
  
  write(100,'(A)') 'Coarse mesh 2D'
  write(100,'(A)') 'Generated by gridgencirc2d'

  write(cbuffer1, '(I10)') nel
  write(cbuffer2, '(I10)') nvt

  write(100,'(A)') trim(adjustl(cbuffer1))//' '//&
                   trim(adjustl(cbuffer2))//' 0 3 1  NEL NVT NMT NVE NBCT'

  !-----------------------------------------------------------------------------
  write(100,'(A)') 'DCORVG'
  !-----------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Copy/convert coordinates from external triangulation
  !---------------------------------------------------------------------------

  open(unit=200, file=trim(adjustl(cfilename))//'.1.node')
  read(200, fmt=*) nvt
  
  ! Read node file line by line and convert coordinates
  do ivt = 1, nvt
    read(200, fmt=*) i, dx, dy, dpar, ibdc
    
    if (ibdc == 0) then
      ! Write inner vertex
      write(cbuffer1, cFormat) dx
      write(cbuffer2, cFormat) dy
    else
      ! Write boundary parametrization
      write(cbuffer1, cFormat) dpar
      write(cbuffer2, cFormat) 0.0_DP
    end if
    
    ! Write line to file
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
  end do

  close(200)

  !-----------------------------------------------------------------------------
  write(100,'(A)') 'KVERT'
  !-----------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Copy elements of inner region from external triangulation
  !---------------------------------------------------------------------------
  
  open(unit=200, file=trim(adjustl(cfilename))//'.1.ele')
  read(200, fmt=*) nel
  
  ! Read element file line by line and convert coordinates
  do iel = 1, nel
    read(200, fmt=*) i, i1, i2, i3
    
    write(cbuffer1, '(I10)') i1
    write(cbuffer2, '(I10)') i2
    write(cbuffer3, '(I10)') i3
    
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//&
                     trim(adjustl(cbuffer2))//' '//&
                     trim(adjustl(cbuffer3))
  end do
  
  close(200)

  !-----------------------------------------------------------------------------
  write(100,'(A)') 'KNPR'
  !-----------------------------------------------------------------------------

  ! Write nodal property for boundary vertices
  do ipoint = 1, npoints
    write(100,'(A)') '1'
  end do
  
  ! Write nodal property for boundary vertices
  do ipoint = npoints+1,nvt
    write(100,'(A)') '0'
  end do
  
  !-----------------------------------------------------------------------------
  write(100,'(A)') 'KMM'
  !-----------------------------------------------------------------------------

  write(cbuffer1, '(I10)') 1
  write(cbuffer2, '(I10)') npoints

  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))

  close(100)

end program gridgenlaval2d
