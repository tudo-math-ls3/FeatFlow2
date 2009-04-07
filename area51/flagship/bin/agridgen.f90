program agridgen
  implicit none

  ! Definition of double precision
  integer, parameter  :: DP = selected_real_kind(15,307)

  ! Definition of mathematical constant PI
  real(DP), parameter :: pi = asin(1.0_DP)*2.0_DP

  ! Definition of the inner radius
  real(DP) :: dinnerRadius = 0.1_DP

  ! Definition of the outer radius
  real(DP) :: douterRadius = 1.0_DP

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
 
  ! local variables
  character(len=80) :: cbuffer
  character(len=32) :: cbuffer1,cbuffer2,cbuffer3,cbuffer4,cbuffer5,cbuffer6
  real(DP) :: x,y,r,phi
  integer :: i,j,k,ivt,jvt,isegment,iaux,irefine,nvt,nel,nmt

  !-----------------------------------------------------------------------------
  ! Get command line arguments
  !-----------------------------------------------------------------------------
  
  do i = 1, command_argument_count()

    call get_command_argument(i, cbuffer)
    select case(cbuffer)
    case('-H','--help')
      write(*,'(A)') 'Usage: agridgen [OPTION]'
      write(*,'(A)') 'Generate an annular grid in TRI/PRM format.'
      write(*,*)
      write(*,'(A,T30,A)') '-R0, --innerradius','radius of the inner circle'
      write(*,'(A,T30,A)') '-R1, --outerradius','radius of the outer circle'
      write(*,'(A,T30,A)') '-S,  --segments','number of segments in azimuthal direction'
      write(*,'(A,T30,A)') '-S0, --minsegments','minimum number of segments in the inner circle'
      write(*,'(A,T30,A)') '-L,  --layers','number of layers in radial direction'
      write(*,'(A,T30,A)') '-F,  --filename','name of the output file'
      write(*,*)
      write(*,'(A)') 'Report bugs to <matthias.moeller@math.tu-dortmund.de>.'
      stop

    case('-R0','--innerradius')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) dinnerRadius
      
    case('-R1','--outerradius')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) douterRadius

    case('-S','--segments')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) nsegments

    case('-S0','--minsegments')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) nminsegments

    case('-L','--layers')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) nlayers

    case('-F','--filename')
      call get_command_argument(i+1, cbuffer)
      read(cbuffer,*) cfilename
    end select
  end do

  !-----------------------------------------------------------------------------
  ! Write statistics
  !-----------------------------------------------------------------------------

  write(*,'(A)') 'Generating annular grid'
  write(*,'(A)') '-----------------------'
  write(cbuffer, *) dinnerRadius
  write(*,'(A,T45,A)') 'inner radius:',trim(adjustl(cbuffer))
  write(cbuffer, *) douterRadius
  write(*,'(A,T45,A)') 'outer radius:',trim(adjustl(cbuffer))
  write(cbuffer, *) nlayers
  write(*,'(A,T45,A)') 'number of layers in radial direction:',trim(adjustl(cbuffer))
  write(cbuffer, *) nsegments
  write(*,'(A,T45,A)') 'number of segments in angular direction:',trim(adjustl(cbuffer))
  write(cbuffer, *) nminsegments
  write(*,'(A,T45,A)') 'minimum number of segments in inner circle:',trim(adjustl(cbuffer))
    write(*,'(A,T45,A)') 'name of output file: ',trim(adjustl(cfilename))
  
  !-----------------------------------------------------------------------------
  ! Generate PRM-file
  !-----------------------------------------------------------------------------

  open(unit=100, file=trim(adjustl(cfilename))//'.prm')

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
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer1))
  write(cbuffer2, cFormat) douterRadius
  write(100,'(A)') trim(adjustl(cbuffer2))//' '//trim(adjustl(cbuffer1))
  write(cbuffer2, cFormat) 2.0_DP*pi
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
  
  close(100)
  
  !-----------------------------------------------------------------------------
  ! Determine the number of inner layers
  !-----------------------------------------------------------------------------
  
  isegment = nsegments; irefine = 0
  inner: do
    if ((mod(isegment, 2) .eq. 0) .and.&
        (isegment .ge. 2*nminsegments)) then
      isegment = isegment/2
      irefine  = irefine+1
    else
      exit inner
    end if
  end do inner
  
  ! Determine number of vertices/elements in the interior circle. To
  ! this we we initialize the quantities NEL, NVT and NMT by the
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

  !-----------------------------------------------------------------------------
  ! Generate TRI-file
  !-----------------------------------------------------------------------------
  open(unit=100, file=trim(adjustl(cfilename))//'.tri')

  write(100,'(A)') 'Coarse mesh 2D'
  write(100,'(A)') 'Generated by gridgen-otype'
  
  write(cbuffer1, '(I)') nsegments*nlayers + nel
  write(cbuffer2, '(I)') nsegments*nlayers + nvt

  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' 0 4 1  NEL NVT NMT NVE NBCT'

  !-----------------------------------------------------------------------------
  write(100,'(A)') 'DCORVG'
  !-----------------------------------------------------------------------------

  ! Write vertex at origin
  write(cbuffer1, cFormat) 0.0
  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer1))
  
  ! Write interior vertices in the inner layer
  iaux = isegment
  do j = 1, 2**irefine
    do i = 1, iaux
      
      ! Compute radius
      r = j*(dinnerRadius)/(2**irefine)

      ! Compute angle
      phi = (i-1)*2.0_DP*pi/iaux

      write(cbuffer1, '(F20.12)') r * cos(phi)
      write(cbuffer2, '(F20.12)') r * sin(phi)

      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
    end do
    
    ! Increate the number of segments
    iaux = iaux + isegment
  end do

  ! Write interior vertices in the outer layer
  do i = 1, nsegments
    do j = 1, nlayers-1

      ! Compute radius
      r = dinnerRadius + j*(douterRadius-dinnerRadius)/nlayers
      
      ! Compute angle
      phi = 2*(i-1)*pi/nsegments

      write(cbuffer1, cFormat) r * cos(phi)
      write(cbuffer2, cFormat) r * sin(phi)

      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
    end do

    ! Write parameter for boundary vertex
    r = real(i-1,DP)/real(nsegments,DP)
    write(cbuffer1, cFormat) r
    write(cbuffer2, cFormat) 0.0
    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
  end do

  !-----------------------------------------------------------------------------
  write(100,'(A)') 'KVERT'
  !-----------------------------------------------------------------------------

  ! Process inner most triangles
  do i = 1, isegment
    
    write(cbuffer1, '(I)') 1
    write(cbuffer2, '(I)') mod(i-1, isegment)+2
    write(cbuffer3, '(I)') mod(i,   isegment)+2
    write(cbuffer4, '(I)') 0

    write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
                     trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
  end do
  

  ! Initialize counter
  iaux = 0
  
  ! Loop over all layers
  do j = 1, 2**irefine-1
    
    iaux = iaux + j

    ! Loop over all segments
    do i = 1, isegment

      ! Compute number of starting vertices in current layer and the
      ! next layer which is located interior to the current one
      ivt = iaux*isegment+(i-1)*j+i+1
      jvt = (iaux-j)*isegment+(i-1)*(j-1)+i+1

      ! Loop over all edges in the current layer
      do k = 1, iaux*isegment+i*j+i+1-ivt

        write(cbuffer1, '(I)') ivt
        write(cbuffer2, '(I)') ivt+1
        write(cbuffer3, '(I)') jvt
        write(cbuffer4, '(I)') 0

        write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
                         trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
        
        ! Increase vertex number in current layer
        ivt = ivt + 1
        
        write(cbuffer1, '(I)') ivt
        ! Check if this is the last=first vertex in the inner layer
        if (jvt+1 .eq. iaux*isegment+2) then
          write(cbuffer2, '(I)') (iaux-j)*isegment+2
        else
          write(cbuffer2, '(I)') jvt+1
        end if
        write(cbuffer3, '(I)') jvt
        write(cbuffer4, '(I)') 0

        write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
                         trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
        
        ! Increase vertex number in inner layer
        jvt = jvt + 1
        
      end do

      write(cbuffer1, '(I)') ivt
      ! Check if this is the last=first vertex in the current layer
      if (ivt+1 .eq. (iaux+j+1)*isegment+2) then
        write(cbuffer2, '(I)') jvt
      else
        write(cbuffer2, '(I)') ivt+1
      end if
      ! Check if this is the last=first vertex in the inner layer
      if (jvt .eq. iaux*isegment+2) then
        write(cbuffer3, '(I)') (iaux-j)*isegment+2
      else
        write(cbuffer3, '(I)') jvt
      end if
      write(cbuffer4, '(I)') 0

      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
                       trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
    end do   
  end do

  ! Process connection layer between unstructured and structured grids
  do i = 1, nsegments

    ! Compute base vertex number
    ivt = nlayers*(i-1)

     write(cbuffer1, '(I)') mod(i-1,         nsegments)+nvt-nsegments+1
     write(cbuffer2, '(I)') mod(ivt,         nsegments*nlayers)+nvt+1
     write(cbuffer3, '(I)') mod(ivt+nlayers, nsegments*nlayers)+nvt+1
     write(cbuffer4, '(I)') mod(i,           nsegments)+nvt-nsegments+1

     write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
                      trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
  end do

  ! Process outer structured grid
  do i = 1, nsegments
    do j = 1, nlayers-1

      ! Compute base vertex number
      ivt = j-1 + nlayers*(i-1)
      
      write(cbuffer1, '(I)') mod(ivt,           nsegments*nlayers)+nvt+1
      write(cbuffer2, '(I)') mod(ivt+1,         nsegments*nlayers)+nvt+1
      write(cbuffer3, '(I)') mod(ivt+nlayers+1, nsegments*nlayers)+nvt+1
      write(cbuffer4, '(I)') mod(ivt+nlayers,   nsegments*nlayers)+nvt+1
      
      write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))//' '//&
                       trim(adjustl(cbuffer3))//' '//trim(adjustl(cbuffer4))
    end do
  end do
  
  !-----------------------------------------------------------------------------
  write(100,'(A)') 'KNPR'
  !-----------------------------------------------------------------------------
  
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
  
  !-----------------------------------------------------------------------------
  write(100,'(A)') 'KMM'
  !-----------------------------------------------------------------------------
  
  write(cbuffer1, '(I)') nlayers+nvt
  write(cbuffer2, '(I)') nlayers*nsegments+nvt

  write(100,'(A)') trim(adjustl(cbuffer1))//' '//trim(adjustl(cbuffer2))
  
  close(100)

end program agridgen
