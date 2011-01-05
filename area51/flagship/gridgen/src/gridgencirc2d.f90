program gridgencirc2d
  implicit none

  ! Definition of double precision
  integer, parameter  :: DP = selected_real_kind(15,307)

  ! Definition of mathematical constant PI
  real(DP) :: pi

  ! Definition of azimuth
  real(DP) :: dazimuth

  ! Definition of the inner radius
  real(DP) :: dinnerRadius = 0.1_DP

  ! Definition of the outer radius
  real(DP) :: douterRadius = 1.0_DP

  ! Definition of the grid anisotropy
  real(DP) :: danisotropy = 1.0_DP

  ! Definition of the number of segments
  integer :: nsegments = 6

  ! Definition of the minimum number of segments in the inner circle
  integer :: nminsegments = 4

  ! Definition of the number of outer layers
  integer :: nlayers = 1

  ! Definition of output filename
  character(len=80) :: cfilename = "grid"

  ! Definition of output format
  character(LEN=*), parameter :: cFormat = '(F32.15)'
 
  ! Definition whether to use external triangulation
  logical :: bexternalDelaunay = .false.

  ! local variables
  character(len=80) :: cbuffer
  character(len=32) :: cbuffer1,cbuffer2,cbuffer3,cbuffer4,cbuffer5,cbuffer6
  real(DP) :: x,y,r,dr0,phi,darea,dlength,bdrPar
  integer :: i,j,k,ivt,jvt,iaux,nvt,nel,nmt,i1,i2,i3
  integer :: isegment,irefine,nlineSegments


  !-----------------------------------------------------------------------------
  ! Initialize mathematical constant(s)
  !-----------------------------------------------------------------------------

  pi = asin(1.0_DP)*2.0_DP; dazimuth = 2.0_DP*pi


  !-----------------------------------------------------------------------------
  ! Get command line arguments
  !-----------------------------------------------------------------------------
  
  do i = 1, command_argument_count()

    call get_command_argument(i, cbuffer)
    select case(cbuffer)
    case('-H','--help')
      write(*,'(A)') 'Usage: gridgen [OPTION]'
      write(*,'(A)') 'Generate an annular grid in TRI/PRM format.'
      write(*,*)
      write(*,'(A,T30,A)') '-A,  --azimuth','azimuth of the circular segment'
      write(*,'(T30,A)')   'If not given, a full circle is adopted.'
      write(*,'(A,T30,A)') '-a,  --anisotropy','degree of radial grid anisotropy in outer region'
      write(*,'(T30,A)')   'If not given, no anisotropy is adopted.'
      write(*,'(A,T30,A)') '-d,  --delaunay','Delaunay triangulation of the inner region'
      write(*,'(T30,A)')   'If not given, no Delaunay triangulation of the inner region is performed.'
      write(*,'(A,T30,A)') '-f,  --filename','name of the output file'
      write(*,'(T30,A)')   'If not given, the generated grid is stored in file "grid.tri/prm".'
      write(*,'(A,T30,A)') '-l,  --layers','number of layers in radial direction in outer region'
      write(*,'(T30,A)')   'If not given, a single layer is adopted.'
      write(*,'(A,T30,A)') '-ri, --innerradius','radius of the inner region'
      write(*,'(T30,A)')   'If not given, 0.1 is adopted as inner radius.'
      write(*,'(A,T30,A)') '-ro, --outerradius','radius of the outer region'
      write(*,'(T30,A)')   'If not given, 1.0 is adopted as inner radius.'
      write(*,'(A,T30,A)') '-s,  --segments','number of segments in azimuthal direction'
      write(*,'(T30,A)')   'If not given, 6 segments are adopted.'
      write(*,'(A,T30,A)') '-si, --minsegments','minimum number of segments in the inner circle'
      write(*,'(T30,A)')   'If not given, 4 segments are adopted. This parameter is not used if'
      write(*,'(T30,A)')   'Delaunay triangulation is used for the inner region.'
      write(*,*)
      write(*,'(A)') 'Report bugs to <matthias.moeller@math.tu-dortmund.de>.'
      stop

    case('-A','--azimuth')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dazimuth; dazimuth = min(2.0_DP*pi, dazimuth)

    case('-a','--anisotropy')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) danisotropy

    case('-d','--delaunay')
      bexternalDelaunay = .true.

    case('-ri','--innerradius')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dinnerRadius
      
    case('-ro','--outerradius')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) douterRadius

    case('-s','--segments')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) nsegments

    case('-si','--minsegments')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) nminsegments

    case('-l','--layers')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) nlayers

    case('-f','--filename')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) cfilename
    end select
  end do


  !-----------------------------------------------------------------------------
  ! Write statistics
  !-----------------------------------------------------------------------------

  write(*,'(A)') 'Generating annular grid'
  write(*,'(A)') '-----------------------'
  write(cbuffer, *) dazimuth
  write(*,'(A,T45,A)') 'azimuth of the circular segment:',trim(adjustl(cbuffer))
  write(cbuffer, *) dinnerRadius
  write(*,'(A,T45,A)') 'radius of inner region:',trim(adjustl(cbuffer))
  write(cbuffer, *) douterRadius
  write(*,'(A,T45,A)') 'radius of outer region:',trim(adjustl(cbuffer))
  write(cbuffer, *) nlayers
  write(*,'(A,T45,A)') 'number of layers in radial direction:',trim(adjustl(cbuffer))
  write(cbuffer, *) nsegments
  write(*,'(A,T45,A)') 'number of segments in azimuthal direction:',trim(adjustl(cbuffer))
  write(cbuffer, *) nminsegments
  write(*,'(A,T45,A)') 'minimum number of segments in inner region:',trim(adjustl(cbuffer))
  write(*,'(A,T45,A)') 'name of output file: ',trim(adjustl(cfilename))

  
  !-----------------------------------------------------------------------------
  ! Generate PRM-file
  !-----------------------------------------------------------------------------

  open(unit=100, file=trim(adjustl(cfilename))//'.prm')

  if (dazimuth .ge. 2.0_DP*pi) then
    !---------------------------------------------------------------------------
    ! Generate full circle with angle: 2*pi
    !---------------------------------------------------------------------------
    write(100,'(A)') 'NBCT'
    write(100,'(A)') '1'
    write(100,'(A)') 'IBCT'
    write(100,'(A)') '1 '
    write(100,'(A)') 'NCOMP'
    write(100,'(A)') '1'
    write(100,'(A)') 'ITYP NSPLINE NPAR'
    write(100,'(A)') '2 1 3 '
    write(100,'(A)') 'PARAMETERS'
    write(cbuffer1, cFormat) 0.0
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer1)) ! midpoint
    write(cbuffer2, cFormat) douterRadius
    write(100,'(A)') trim(adjustl(cbuffer2))//' '//trim(adjustl(cbuffer1)) ! radius
    write(cbuffer2, cFormat) 2.0_DP*pi
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! azimuth
    
  else
    !---------------------------------------------------------------------------
    ! Generate circle segment with angle: dazimuth
    !---------------------------------------------------------------------------
    write(100,'(A)') 'NBCT'
    write(100,'(A)') '1'
    write(100,'(A)') 'IBCT'
    write(100,'(A)') '1 '
    write(100,'(A)') 'NCOMP'
    write(100,'(A)') '3'
    write(100,'(A)') 'ITYP NSPLINE NPAR'
    write(100,'(A)') '2 1 3 '
    write(100,'(A)') '1 1 2 '
    write(100,'(A)') '1 1 2 '
    write(100,'(A)') 'PARAMETERS'
    ! Circular segment
    write(cbuffer1, cFormat) 0.0
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer1)) ! midpoint
    write(cbuffer2, cFormat) douterRadius
    write(100,'(A)') trim(adjustl(cbuffer2))//' '//trim(adjustl(cbuffer1)) ! radius
    write(cbuffer2, cFormat) 2.0_DP*pi
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! azimuth
    ! Line segment
    write(cbuffer1, cFormat) douterRadius*cos(dazimuth)
    write(cbuffer2, cFormat) douterRadius*sin(dazimuth)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! coordinate
    write(cbuffer1, cFormat) -douterRadius*cos(dazimuth)
    write(cbuffer2, cFormat) -douterRadius*sin(dazimuth)
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2)) ! direction
    ! Line segment
    write(cbuffer1, cFormat) 0.0
    write(cbuffer2, cFormat) douterRadius
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer1)) ! origin
    write(100,'(A)') trim(adjustl(cbuffer2))//' '//trim(adjustl(cbuffer1)) ! direction
        
  end if

  close(100)


  !-----------------------------------------------------------------------------
  ! Prepare generation of TRI-file
  !-----------------------------------------------------------------------------

  if (bexternalDelaunay) then

    !---------------------------------------------------------------------------
    ! Generate inner triangulation by the external program triangle
    !---------------------------------------------------------------------------
    
    open(unit=100, file=trim(adjustl(cfilename))//'.node')
    
    ! Write file header
    write(100, '(A)') "# "//trim(adjustl(cfilename))//'.node'
    write(100, '(A)') "#"
    write(100, '(A)') "# Node file generated by gridgencirc2d."
    write(100, '(A)') "#"

    if (dazimuth .ge. 2.0_DP*pi) then
      !-------------------------------------------------------------------------
      ! Generate full circle with angle: 2*pi
      !-------------------------------------------------------------------------

      ! Write number of points
      write(100, '(A)') "# Number of points, 2D, one attribute, boundary markers"
      write(cbuffer1, '(I10)') nsegments+1
      write(100, '(A)') trim(adjustl(cbuffer1))//' 2 1 1'
      
      ! Write vertex coordinates
      do isegment = 1, nsegments
        
        ! Compute azimuth
        phi = (isegment-1)*dazimuth/nsegments

        ! Fill buffers
        write(cbuffer1, '(I10)') isegment
        write(cbuffer2, cFormat) dinnerRadius * cos(phi)
        write(cbuffer3, cFormat) dinnerRadius * sin(phi)
        write(cbuffer4, cFormat) 0.0_DP
        
        ! Write line to file
        write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                          trim(adjustl(cbuffer2))//' '//&
                          trim(adjustl(cbuffer3))//' '//&
                          trim(adjustl(cbuffer4))//' 0'
      end do
      
      ! Write single node at the origin
      write(cbuffer1, '(I10)') nsegments+1
      write(cbuffer2, cFormat) 0.0_DP
      write(cbuffer3, cFormat) 0.0_DP
      write(cbuffer4, cFormat) 0.0_DP
      
      ! Write line to file
      write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                        trim(adjustl(cbuffer2))//' '//&
                        trim(adjustl(cbuffer3))//' '//&
                        trim(adjustl(cbuffer4))//' 0'
      
      close(100)

      ! Compute minimum edge length
      phi = dazimuth/nsegments
      darea = (dinnerradius**2) * ((cos(phi)-1)**2 + sin(phi)**2)
      
      ! Generate conforming Delaunay triangulation (-D) subject to a
      ! user-defined element area (-a) without generating additional
      ! Steiner points on the boundary (-Y) based on vertex distribution
      write(cbuffer1, cFormat) darea
      call system('triangle -V -Y -D -e -a'//&
          trim(adjustl(cbuffer1))//' '//trim(adjustl(cfilename))//'.node')

    else
      !-------------------------------------------------------------------------
      ! Generate circle segment with angle = dazimuth
      !-------------------------------------------------------------------------
      
      ! Compute number of line segments
      dlength = sqrt((dinnerRadius*cos(dazimuth/nsegments)-dinnerRadius)**2 +&
                     (dinnerRadius*sin(dazimuth/nsegments))**2)
      nlineSegments = floor(dinnerRadius/dlength)
      
      ! Write number of points
      write(100, '(A)') "# Number of points, 2D, one attribute, boundary markers"
      write(cbuffer1, '(I10)') nsegments+2*nlineSegments
      write(100, '(A)') trim(adjustl(cbuffer1))//' 2 1 1'

      ! Write vertex coordinates
      do isegment = 1, nsegments+1
        
        ! Compute azimuth
        phi = (isegment-1)*dazimuth/nsegments
        
        if (isegment .eq. 1) then
          ! Fill buffers
          write(cbuffer1, '(I10)') isegment
          write(cbuffer2, cFormat) dinnerRadius * cos(phi)
          write(cbuffer3, cFormat) dinnerRadius * sin(phi)
          write(cbuffer4, cFormat) 2.0_DP+dinnerRadius/douterRadius

          ! Write line to file
          write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                            trim(adjustl(cbuffer2))//' '//&
                            trim(adjustl(cbuffer3))//' '//&
                            trim(adjustl(cbuffer4))//' 2'

        elseif (isegment .eq. nsegments+1) then
          ! Fill buffers
          write(cbuffer1, '(I10)') isegment
          write(cbuffer2, cFormat) dinnerRadius * cos(phi)
          write(cbuffer3, cFormat) dinnerRadius * sin(phi)
          write(cbuffer4, cFormat) 2.0_DP-dinnerRadius/douterRadius

          ! Write line to file
          write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                            trim(adjustl(cbuffer2))//' '//&
                            trim(adjustl(cbuffer3))//' '//&
                            trim(adjustl(cbuffer4))//' 2'
        else
          ! Fill buffers
          write(cbuffer1, '(I10)') isegment
          write(cbuffer2, cFormat) dinnerRadius * cos(phi)
          write(cbuffer3, cFormat) dinnerRadius * sin(phi)
          write(cbuffer4, cFormat) 0.0_DP
          
          ! Write line to file
          write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                            trim(adjustl(cbuffer2))//' '//&
                            trim(adjustl(cbuffer3))//' '//&
                            trim(adjustl(cbuffer4))//' 0'
        end if
      end do
      
      ! Write line segment from last point on circle segment to origin
      x = dinnerRadius * cos(dazimuth)
      y = dinnerRadius * sin(dazimuth)

      ! Write vertex coordinates
      do isegment = 1, nlineSegments

        write(cbuffer1, '(I10)') nsegments+1+isegment
        write(cbuffer2, cFormat) (1.0_DP-isegment/real(nlineSegments,DP))*x
        write(cbuffer3, cFormat) (1.0_DP-isegment/real(nlineSegments,DP))*y
        write(cbuffer4, cFormat) (dinnerRadius*real(isegment,DP)/real(nlineSegments,DP)+&
                                  douterRadius-dinnerRadius)/douterRadius+1.0_DP

        write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                          trim(adjustl(cbuffer2))//' '//&
                          trim(adjustl(cbuffer3))//' '//&
                          trim(adjustl(cbuffer4))//' 2'
      end do

      ! Write line segment from origin to first point on circle segment
      x = dinnerRadius; y = 0.0_DP

      ! Write vertex coordinates
      do isegment = 1, nlineSegments-1

        write(cbuffer1, '(I10)') nsegments+1+nlineSegments+isegment
        write(cbuffer2, cFormat) real(isegment,DP)/real(nlineSegments,DP)*x
        write(cbuffer3, cFormat) real(isegment,DP)/real(nlineSegments,DP)*y
        write(cbuffer4, cFormat) real(isegment,DP)/real(nlineSegments,DP)*&
                                 dinnerRadius/douterRadius+2.0_DP

        write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                          trim(adjustl(cbuffer2))//' '//&
                          trim(adjustl(cbuffer3))//' '//&
                          trim(adjustl(cbuffer4))//' 2'
      end do
      
      close(100)
      
      ! Write PSLG (planar straight line graph) to file
      open(unit=100, file=trim(adjustl(cfilename))//'.poly')
      
      ! Write file header
      write(100, '(A)') "# "//trim(adjustl(cfilename))//'.poly'
      write(100, '(A)') "#"
      write(100, '(A)') "# PSLG file generated by gridgencirc2d."
      write(100, '(A)') "#"
      
      ! Write number of points
      write(100, '(A)') "# no points, 2D, no attributes, no boundary markers"
      write(100, '(A)') '0 2 0 0'
      write(100, '(A)') "# Number of segments, no boundary markers"
      write(cbuffer1, '(I10)') nsegments+1+2*nlineSegments
      write(100, '(A)') trim(adjustl(cbuffer1))
      
      ! Write segments
      do isegment = 1, nsegments+2*nlineSegments+1
        
        write(cbuffer1, '(I10)') isegment
        write(cbuffer2, '(I10)') mod(isegment,nsegments+2*nlineSegments)+1
        
        write(100, '(A)') trim(adjustl(cbuffer1))//' '//&
                          trim(adjustl(cbuffer1))//' '//&
                          trim(adjustl(cbuffer2))
      end do
      
      ! Write number of holes
      write(100, '(A)') '0'
      write(100, '(A)') '0'
      
      close(100)
      
      ! Compute minimum edge length
      phi = dazimuth/nsegments
      darea = (dinnerradius**2) * ((cos(phi)-1)**2 + sin(phi)**2)
      
      ! Generate conforming Delaunay triangulation (-D) subject to a
      ! user-defined element area (-a) without generating additional
      ! Steiner points on the boundary (-Y) based on vertex distribution
      write(cbuffer1, cFormat) darea
      call system('triangle -V -Y -D -e -a'//&
          trim(adjustl(cbuffer1))//' '//trim(adjustl(cfilename))//'.node'//&
          ' -p '//trim(adjustl(cfilename))//'.poly')
    end if

    !---------------------------------------------------------------------------
    ! Read data for inner triangulation and convert it into TRI format
    !---------------------------------------------------------------------------
    
    open(unit=100, file=trim(adjustl(cfilename))//'.1.node')
    read(100, fmt=*) nvt
    close(100)

    open(unit=100, file=trim(adjustl(cfilename))//'.1.ele')
    read(100, fmt=*) nel
    close(100)
    
    open(unit=100, file=trim(adjustl(cfilename))//'.1.edge')
    read(100, fmt=*) nmt
    close(100)

  else
   
    !---------------------------------------------------------------------------
    ! Determine the number of inner layers manually
    !---------------------------------------------------------------------------
    
    ! Determine the number of regular coarsening steps to obtain the
    ! number of vertices/elements/midpoints of the coarse grid.
    isegment = nsegments; irefine  = 0

    inner: do
      if ((mod(isegment, 2) .eq. 0) .and.&
          (isegment .ge. 2*nminsegments)) then
        isegment = isegment/2
        irefine  = irefine+1
      else
        exit inner
      end if
    end do inner
    
    ! Determine number of vertices/elements in the interior circle.
    ! Thus, we initialize the quantities NEL, NVT and NMT by the
    ! values of the initial coarse grid in perform regular subdivision.
    nel = isegment
    nvt = isegment+1
    nmt = 2*isegment
    
    do i = 1, irefine
      ! Each edge produces a new vertex
      nvt = nvt + nmt
      
      ! Each edge is subdivided into two new edges and 
      ! each element produces three new edges
      nmt = 2*nmt + 3*nel
      
      ! Each element is subdivided into four new elements
      nel = 4*nel
    end do

    ! Update number of vertices/midpoints for circle sections
    if (dazimuth .lt. 2.0_DP*pi) then
      nvt = nvt+2**irefine
      nmt = nmt+2**(irefine+1)
    end if

  end if

  
  !-----------------------------------------------------------------------------
  ! Generate TRI-file
  !-----------------------------------------------------------------------------

  open(unit=100, file=trim(adjustl(cfilename))//'.tri')

  write(100,'(A)') 'Coarse mesh 2D'
  write(100,'(A)') 'Generated by gridgencirc2d'
  
  write(cbuffer1, '(I10)') nel + nsegments*nlayers
  write(cbuffer2, '(I10)') nvt + nsegments*nlayers +&
                           merge(0, nlayers, dazimuth .ge. 2.0_DP*pi)

  write(100,'(A)') trim(adjustl(cbuffer1))//' '//&
                   trim(adjustl(cbuffer2))//' 0 4 1  NEL NVT NMT NVE NBCT'

  !-----------------------------------------------------------------------------
  write(100,'(A)') 'DCORVG'
  !-----------------------------------------------------------------------------

  if (bexternalDelaunay) then

    !---------------------------------------------------------------------------
    ! Copy/convert coordinates of inner region from external triangulation
    !---------------------------------------------------------------------------

    ! Write interior vertices in the inner layer
    open(unit=200, file=trim(adjustl(cfilename))//'.1.node')
    read(200, fmt=*) nvt
    
    ! Read node file line by line and convert coordinates
    do ivt = 1, nvt
      read(200, fmt=*) i, x, y, bdrPar, iaux
      
      if (iaux .gt. 1) then
        ! Write boundary parametrization
        write(cbuffer1, cFormat) bdrPar
        write(cbuffer2, cFormat) 0.0_DP
      else
        ! Write inner vertex
        write(cbuffer1, cFormat) x
        write(cbuffer2, cFormat) y
      end if
      
      ! Write line to file
      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
    end do
    
    close(200)
    
  else
    !---------------------------------------------------------------------------
    ! Create coordinates of inner region manually
    !---------------------------------------------------------------------------

    ! Write vertex at origin
    if (dazimuth .ge. 2.0_DP*pi) then
      write(cbuffer1, cFormat) 0.0_DP
      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer1))
    else
      write(cbuffer1, cFormat) 2.0_DP
      write(cbuffer2, cFormat) 0.0_DP
      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
    end if
    
    ! Write vertices in the inner layer
    iaux = isegment

    ! Loop over all refinement steps
    do j = 1, 2**irefine
      
      ! Compute radius of current layer
      r = j*(dinnerRadius)/(2**irefine)
      
      ! Loop over all segments in current layer
      do i = 1, iaux
        
        ! Compute azimuth
        phi = (i-1)*dazimuth/iaux
        
        if ((dazimuth .lt. 2.0_DP*pi) .and. (i .eq. 1)) then
          ! Fill buffers
          write(cbuffer1, cFormat) 2.0_DP+r/douterRadius
          write(cbuffer2, cFormat) 0.0_DP
        else
          ! Fill buffers
          write(cbuffer1, cFormat) r * cos(phi)
          write(cbuffer2, cFormat) r * sin(phi)
        end if
        
        ! Write line to file
        write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
      end do

      ! Since the circle is not closed we have to add one extra vertex
      if (dazimuth .lt. 2.0_DP*pi) then
      ! Fill buffers
      write(cbuffer1, cFormat) 2.0_DP-r/douterRadius
        write(cbuffer2, cFormat) 0.0_DP

        ! Write line to file
        write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
      end if
      
      ! Increate the number of segments
      iaux = iaux + isegment
    end do
    
  end if
  
  ! Compute minimal grid spacing (isotropic/anisotropic ?)
  if (danisotropy .eq. 1.0_DP) then
    dr0 = (douterRadius-dinnerRadius)/nlayers
  else
    ! Compute local anisotropy factor
    danisotropy = exp(log(danisotropy)/(nlayers-1.0_DP))
    dr0 = (douterRadius-dinnerRadius)*(danisotropy-1)/(danisotropy**nlayers-1)
  end if
  
  
  do i = 1, merge(nsegments, nsegments+1, dazimuth .ge. 2.0_DP*pi)
    
    !---------------------------------------------------------------------------
    ! Write interior vertices in the outer layer
    !---------------------------------------------------------------------------
    
    do j = 1, nlayers-1
      
      ! Compute radius
      if (danisotropy .eq. 1.0_DP) then
        r = dinnerRadius + dr0*j
      else
        r = dinnerRadius + dr0*(danisotropy**(j)-1)/(danisotropy-1)
      end if
      
      ! Compute azimuth
      phi = dazimuth*(i-1)/nsegments
      
      if (dazimuth .ge. 2.0_DP*pi) then
        ! Fill buffers
        write(cbuffer1, cFormat) r * cos(phi)
        write(cbuffer2, cFormat) r * sin(phi)
        
        ! Write line to file
        write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
      else
        if (i .eq. 1) then
          ! Fill buffer
          write(cbuffer1, cFormat) 2.0_DP+r/douterRadius
          write(cbuffer2, cFormat) 0.0
        elseif (i .eq. nsegments+1) then
          ! Fill buffer
          write(cbuffer1, cFormat) 2.0_DP-r/douterRadius
          write(cbuffer2, cFormat) 0.0
        else
          ! Fill buffers
          write(cbuffer1, cFormat) r * cos(phi)
          write(cbuffer2, cFormat) r * sin(phi)
        end if

        ! Write line to file
        write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
      end if
    end do
    
    !---------------------------------------------------------------------------
    ! Write parameters for boundary vertex
    !---------------------------------------------------------------------------
    
    r = real(i-1,DP)/real(nsegments,DP)
    write(cbuffer1, cFormat) r
    write(cbuffer2, cFormat) 0.0
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
  end do

  !-----------------------------------------------------------------------------
  write(100,'(A)') 'KVERT'
  !-----------------------------------------------------------------------------

  if (bexternalDelaunay) then

    !---------------------------------------------------------------------------
    ! Copy elements of inner region from external triangulation
    !---------------------------------------------------------------------------

    open(unit=200, file=trim(adjustl(cfilename))//'.1.ele')
    read(200, fmt=*) nel

    ! Read element file line by line and convert coordinates
    do i = 1, nel
      read(200, fmt=*) iaux, i1, i2, i3

      write(cbuffer1, '(I10)') i1
      write(cbuffer2, '(I10)') i2
      write(cbuffer3, '(I10)') i3
      write(cbuffer4, '(I10)') 0

      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
          trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
    end do
    
    close(200)

    !---------------------------------------------------------------------------
    ! Process connection layer between unstructured and structured grids
    !---------------------------------------------------------------------------

    do i = 1, nsegments
      
      ! Compute base vertex number
      ivt = nlayers*(i-1)
      
      if (dazimuth .ge. 2.0_DP*pi) then
        write(cbuffer1, '(I10)') mod(i-1,         nsegments)+1
        write(cbuffer2, '(I10)') mod(ivt,         nsegments*nlayers)+nvt+1
        write(cbuffer3, '(I10)') mod(ivt+nlayers, nsegments*nlayers)+nvt+1
        write(cbuffer4, '(I10)') mod(i,           nsegments)+1
      else
        write(cbuffer1, '(I10)') mod(i-1,         (nsegments+1))+1
        write(cbuffer2, '(I10)') mod(ivt,         (nsegments+1)*nlayers)+nvt+1
        write(cbuffer3, '(I10)') mod(ivt+nlayers, (nsegments+1)*nlayers)+nvt+1
        write(cbuffer4, '(I10)') mod(i,           (nsegments+1))+1
      end if
      
      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
          trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
    end do
    
  else
    
    !---------------------------------------------------------------------------
    ! Create triangles of inner region manually
    !---------------------------------------------------------------------------
    
    ! Process inner most triangles
    do i = 1, isegment
      
      if (dazimuth .ge. 2.0_DP*pi) then
        write(cbuffer1, '(I10)') 1
        write(cbuffer2, '(I10)') mod(i-1, isegment)+2
        write(cbuffer3, '(I10)') mod(i,   isegment)+2
        write(cbuffer4, '(I10)') 0
      else
        write(cbuffer1, '(I10)') 1
        write(cbuffer2, '(I10)') mod(i-1, isegment+1)+2
        write(cbuffer3, '(I10)') mod(i,   isegment+1)+2
        write(cbuffer4, '(I10)') 0
      end if
      
      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
          trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
    end do
    
    
    ! Initialize counter and process triangles in interior region
    iaux = 0
    
    ! Loop over all layers
    do j = 1, 2**irefine-1
      
      ! Update counter
      iaux = iaux + j
      
      ! Loop over all segments of the inner most coarse grid
      do i = 1, isegment
        
        ! Compute number of starting vertices in current layer and the
        ! next layer which is located interior to the current one
        if (dazimuth .ge. 2.0_DP*pi) then
          ivt = iaux*isegment+(i-1)*j+i+1
          jvt = (iaux-j)*isegment+(i-1)*(j-1)+i+1
        else
          ivt = iaux*isegment+i*(j+1)+1 
          jvt = (iaux-j)*isegment+i*j+1
        end if
        
        ! Loop over all edges in the current layer
        do k = 1, iaux*isegment+i*(j+1)+1-ivt+&
                  merge(0, j, dazimuth .ge. 2.0_DP*pi)
          
          write(cbuffer1, '(I10)') ivt
          write(cbuffer2, '(I10)') ivt+1
          write(cbuffer3, '(I10)') jvt
          write(cbuffer4, '(I10)') 0
          
          write(100,'(A)') trim(adjustl(cbuffer1))//' '//&
                           trim(adjustl(cbuffer2))//' '//&
                           trim(adjustl(cbuffer3))//' '//&
                           trim(adjustl(cbuffer4))

          ! Increase vertex number in current layer
          ivt = ivt + 1
          
          write(cbuffer1, '(I10)') ivt
          ! Check if this is the last=first vertex in the inner layer
          if (jvt+1 .eq. iaux*isegment+&
                         merge(2, 2+j, dazimuth .ge. 2.0_DP*pi)) then
            write(cbuffer2, '(I10)') (iaux-j)*isegment+2
          else
            write(cbuffer2, '(I10)') jvt+1
          end if
          write(cbuffer3, '(I10)') jvt
          write(cbuffer4, '(I10)') 0
          
          write(100,'(A)') trim(adjustl(cbuffer1))//' '//&
                           trim(adjustl(cbuffer2))//' '//&
                           trim(adjustl(cbuffer3))//' '//&
                           trim(adjustl(cbuffer4))
       
          ! Increase vertex number in inner layer
          jvt = jvt + 1
          
        end do
        
        write(cbuffer1, '(I10)') ivt
        ! Check if this is the last=first vertex in the current layer
        if (ivt+1 .eq. (iaux+j+1)*isegment+2) then
          write(cbuffer2, '(I10)') jvt
        else
          write(cbuffer2, '(I10)') ivt+1
        end if
        ! Check if this is the last=first vertex in the inner layer
        if (jvt .eq. iaux*isegment+&
                     merge(2, 2+j, dazimuth .ge. 2.0_DP*pi)) then
          write(cbuffer3, '(I10)') (iaux-j)*isegment+2
        else
          write(cbuffer3, '(I10)') jvt
        end if
        write(cbuffer4, '(I10)') 0
        
        write(100,'(A)') trim(adjustl(cbuffer1))//' '//&
                         trim(adjustl(cbuffer2))//' '//&
                         trim(adjustl(cbuffer3))//' '//&
                         trim(adjustl(cbuffer4))
      end do
    end do
    
    !---------------------------------------------------------------------------
    ! Process connection layer between unstructured and structured grids
    !---------------------------------------------------------------------------
    do i = 1, nsegments
      
      ! Compute base vertex number
      ivt = nlayers*(i-1)

      if (dazimuth .ge. 2.0_DP*pi) then
        write(cbuffer1, '(I10)') mod(i-1,         nsegments)+nvt-nsegments+1
        write(cbuffer2, '(I10)') mod(ivt,         nsegments*nlayers)+nvt+1
        write(cbuffer3, '(I10)') mod(ivt+nlayers, nsegments*nlayers)+nvt+1
        write(cbuffer4, '(I10)') mod(i,           nsegments)+nvt-nsegments+1
      else
        write(cbuffer1, '(I10)') mod(i-1,         (nsegments+1))+nvt-nsegments
        write(cbuffer2, '(I10)') mod(ivt,         (nsegments+1)*nlayers)+nvt+1
        write(cbuffer3, '(I10)') mod(ivt+nlayers, (nsegments+1)*nlayers)+nvt+1
        write(cbuffer4, '(I10)') mod(i,           (nsegments+1))+nvt-nsegments
      end if
      
      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
          trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
    end do
    
  end if
  
  !-----------------------------------------------------------------------------
  ! Process outer structured grid
  !-----------------------------------------------------------------------------
  do i = 1, nsegments
    do j = 1, nlayers-1

      ! Compute base vertex number
      ivt = j-1 + nlayers*(i-1)
      
      if (dazimuth .ge. 2.0_DP*pi) then
        write(cbuffer1, '(I10)') mod(ivt,           nsegments*nlayers)+nvt+1
        write(cbuffer2, '(I10)') mod(ivt+1,         nsegments*nlayers)+nvt+1
        write(cbuffer3, '(I10)') mod(ivt+nlayers+1, nsegments*nlayers)+nvt+1
        write(cbuffer4, '(I10)') mod(ivt+nlayers,   nsegments*nlayers)+nvt+1
      else
        write(cbuffer1, '(I10)') mod(ivt,           (nsegments+1)*nlayers)+nvt+1
        write(cbuffer2, '(I10)') mod(ivt+1,         (nsegments+1)*nlayers)+nvt+1
        write(cbuffer3, '(I10)') mod(ivt+nlayers+1, (nsegments+1)*nlayers)+nvt+1
        write(cbuffer4, '(I10)') mod(ivt+nlayers,   (nsegments+1)*nlayers)+nvt+1
      end if
      
      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
                       trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
    end do
  end do
  
  !-----------------------------------------------------------------------------
  write(100,'(A)') 'KNPR'
  !-----------------------------------------------------------------------------
  
  if (dazimuth .ge. 2.0_DP*pi) then
    
    !---------------------------------------------------------------------------
    ! Full circle
    !---------------------------------------------------------------------------
    
    ! Write nodal property for interior vertices
    do i = 1, nvt
      write(100,'(A)') '0'
    end do
    
    do i = 1, nsegments
      ! Write nodal property for interior vertices
      do j = 1, nlayers-1
        write(100,'(A)') '0'
      end do
      ! Write nodal property for boundary vertices
      write(100,'(A)') '1'
    end do

  else

    !---------------------------------------------------------------------------
    ! Circle segment
    !---------------------------------------------------------------------------

    if (bexternalDelaunay) then
      
      ! Write nodal property for vertices in inner region
      open(unit=200, file=trim(adjustl(cfilename))//'.1.node')
      read(200, fmt=*) nvt
      do ivt = 1, nvt
        read(200, fmt=*) i, x, y, bdrPar, iaux
        write(100,'(A)') merge('1', '0', iaux .gt. 1)
      end do
      close(200)

    else ! Create nodal properties manually

      ! Write nodal property for vertex at origin
      write(100,'(A)') '1'

      ! Write nodal properties for vertices in the inner layer
      iaux = isegment

      ! Loop over all refinement steps
      do j = 1, 2**irefine

        ! Loop over all segments in current layer
        do i = 1, iaux
          
          if (i .eq. 1) then
            write(100,'(A)') '1'
          else
            write(100,'(A)') '0'
          end if
        end do
        
        ! Since the circle is not closed we have to add one extra vertex
        write(100,'(A)') '1'

        ! Increate the number of segments
        iaux = iaux + isegment
      end do
    
    end if

    ! Write nodal properties for vertices in the outer layer
    do i = 1, nsegments+1
      do j = 1, nlayers-1
        if ((i .eq. 1) .or. (i .eq. nsegments+1)) then
          ! Write nodal property for boundary vertex
          write(100,'(A)') '1'
        else
          ! Write nodal property for interior vertex
          write(100,'(A)') '0'
        end if
      end do
      
      ! Write nodal property for boundary vertex
      write(100,'(A)') '1'
    end do
    
  end if
  
  !-----------------------------------------------------------------------------
  write(100,'(A)') 'KMM'
  !-----------------------------------------------------------------------------
  
  if (dazimuth .ge. 2.0_DP*pi) then
    write(cbuffer1, '(I10)') nlayers+nvt
    write(cbuffer2, '(I10)') nlayers*nsegments+nvt
  else
    write(cbuffer1, '(I10)') nlayers+nvt
    write(cbuffer2, '(I10)') nlayers+nvt-1
  end if

  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
  
  close(100)

end program gridgencirc2d
