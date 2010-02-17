!##############################################################################
!# ****************************************************************************
!# <name> euler_lagrange </name>
!# ****************************************************************************


module euler_lagrange

  use afcstabilisation
  use bilinearformevaluation
  use boundary
  use boundaryfilter
  use collection
  use derivatives
  use element
  use eulerlagrange_basic
  use eulerlagrange_callback
  use eulerlagrange_callback1d
  use eulerlagrange_callback2d
  use eulerlagrange_callback3d
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use graph
  use groupfemsystem
  use hadaptaux
  use hadaptivity
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use pprocerror
  use pprocgradients
  use pprocindicator
  use pprocsolution
  use problem
  use solveraux
  use spatialdiscretisation
  use statistics
  use stdoperators
  use storage
  use thermodynamics
  use timestep
  use timestepaux
  use ucd
  use triasearch
  use basicgeometry
  use geometryaux


   implicit none

   private
   public :: eulerlagrange_init
   public :: eulerlagrange_step
   public :: calculatebarycoords
   public :: findnewelement
   public :: wrongelement
   public :: moveparticle
   public :: checkboundary
   public :: calculatevolumepart

   type , public :: t_Particles
      !number of particles
      integer :: nPart 
      ! element
      integer(I32) :: h_element
      integer, dimension(:), pointer :: p_element
      ! position
      integer(I32) :: h_xpos, h_ypos, h_zpos
      real(DP), dimension(:), pointer :: p_xpos, p_ypos, p_zpos
      ! old position
      integer(I32) :: h_xpos_old, h_ypos_old, h_zpos_old
      real(DP), dimension(:), pointer :: p_xpos_old, p_ypos_old, p_zpos_old
      ! velocity
      integer(I32) :: h_xvelo, h_yvelo, h_zvelo
      real(DP), dimension(:), pointer :: p_xvelo, p_yvelo, p_zvelo
      ! old velocity
      integer(I32) :: h_xvelo_old, h_yvelo_old, h_zvelo_old
      real(DP), dimension(:), pointer :: p_xvelo_old, p_yvelo_old, p_zvelo_old
      ! velocity of the gas
      integer(I32) :: h_xvelo_gas, h_yvelo_gas, h_zvelo_gas
      real(DP), dimension(:), pointer :: p_xvelo_gas, p_yvelo_gas, p_zvelo_gas
      ! old velocity of the gas
      integer(I32) :: h_xvelo_gas_old, h_yvelo_gas_old, h_zvelo_gas_old
      real(DP), dimension(:), pointer :: p_xvelo_gas_old, p_yvelo_gas_old, p_zvelo_gas_old
      ! barycentric coordinates
      integer(I32) :: h_lambda1, h_lambda2, h_lambda3, h_lambda4
      real(DP), dimension(:), pointer :: p_lambda1, p_lambda2, p_lambda3, p_lambda4
      ! diameter, mass and alpha_n 
      integer(I32) :: h_diam, h_mass, h_alpha_n
      real(DP), dimension(:), pointer :: p_diam, p_mass, p_alpha_n
      ! midpoints of the element 
      integer(I32) :: h_midpoints_el
      real(DP), dimension(:,:), pointer :: p_midpoints_el
      ! volumepart of the particles 
      integer(I32) :: h_PartVol
      real(DP), dimension(:), pointer :: p_PartVol
      ! volumepart of the particles 
      integer(I32) :: h_PartVelo
      real(DP), dimension(:,:), pointer :: p_PartVelo
      ! gravity
      real(DP), dimension(2)  :: gravity
      ! viscosity of the gas
      real(DP) :: nu_g
      ! parameter for particle-wall collisions
      real(DP)  :: tang_val, norm_val
      ! variables for particle-wall collisions
      integer(I32) :: h_bdy_time, h_bdy_check
      real(DP), dimension(:), pointer :: p_bdy_time, p_bdy_check
      ! timestep for video
      integer :: iTimestep
      ! maximum value in x-direction
      real(DP):: maxvalx
    end type t_Particles


    contains

subroutine eulerlagrange_init(rparlist,p_rproblemLevel,rsolution,rtimestep,rcollection,rParticles)

    ! parameterlist
    type(t_parlist), intent(inout) :: rparlist
    
    ! collection structure
    type(t_collection), intent(inout) :: rcollection

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! primal solution vector
    type(t_vectorBlock), intent(inout), target :: rsolution

    ! particles
    type(t_Particles), intent(inout) :: rParticles

    ! local variables
    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement
 
  ! local variables
  integer :: ivt, iPart, iel

  real(DP) :: random1, random2

  ! midpoints of the elements
  integer(I32) :: h_midpoints
  real(DP), dimension(:,:), pointer :: p_midpoints_el
  integer, dimension(2) :: md_el_length

  ! startingpostions of the particles
  real(DP) :: partxmin, partxmax, partymin, partymax

  ! velocity, mass and diameter of the particles
  real(DP) :: velopartx, veloparty, particlediam, particlemass

  ! gravity
  real(DP) :: gravityx, gravityy

  ! quantity of particles
  integer :: nPart
  
  ! boundarybehaviour
  integer :: boundbehav

  ! kinematic viscosity of the gas
  real(DP) :: gas_nu
  
  ! Set pointer to triangulation
  p_rtriangulation => p_rproblemLevel%rtriangulation
  
  ! get quantity of particles
  call parlst_getvalue_int(rparlist, 'Eulerlagrange', "nPart", nPart)

  rParticles%npart = nPart
  
  ! storage_new (scall, sname, isize, ctype, ihandle, cinitNewBlock)
  call storage_new ('euler_lagrange', 'Particle:element', rParticles%npart, ST_INT, rParticles%h_element, &
                            ST_NEWBLOCK_NOINIT)

  call storage_new ('euler_lagrange', 'Particle:lambda1', rParticles%npart, ST_DOUBLE, rParticles%h_lambda1, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:lambda2', rParticles%npart, ST_DOUBLE, rParticles%h_lambda2, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:lambda3', rParticles%npart, ST_DOUBLE, rParticles%h_lambda3, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:lambda4', rParticles%npart, ST_DOUBLE, rParticles%h_lambda4, &
                            ST_NEWBLOCK_NOINIT)

  call storage_new ('euler_lagrange', 'Particle:diameter', rParticles%npart, ST_DOUBLE, rParticles%h_diam, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:mass', rParticles%npart, ST_DOUBLE, rParticles%h_mass, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:alpha_n', rParticles%npart, ST_DOUBLE, rParticles%h_alpha_n, &
                            ST_NEWBLOCK_NOINIT)

  call storage_new ('euler_lagrange', 'Particle:xpos', rParticles%npart, ST_DOUBLE, rParticles%h_xpos, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:ypos', rParticles%npart, ST_DOUBLE, rParticles%h_ypos, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:zpos', rParticles%npart, ST_DOUBLE, rParticles%h_zpos, &
                            ST_NEWBLOCK_NOINIT)

  call storage_new ('euler_lagrange', 'Particle:xpos_old', rParticles%npart, ST_DOUBLE, rParticles%h_xpos_old, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:ypos_old', rParticles%npart, ST_DOUBLE, rParticles%h_ypos_old, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:zpos_old', rParticles%npart, ST_DOUBLE, rParticles%h_zpos_old, &
                            ST_NEWBLOCK_NOINIT)

  call storage_new ('euler_lagrange', 'Particle:xvelo', rParticles%npart, ST_DOUBLE, rParticles%h_xvelo, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:yvelo', rParticles%npart, ST_DOUBLE, rParticles%h_yvelo, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:zvelo', rParticles%npart, ST_DOUBLE, rParticles%h_zvelo, &
                            ST_NEWBLOCK_NOINIT)

  call storage_new ('euler_lagrange', 'Particle:xvelo_old', rParticles%npart, ST_DOUBLE, rParticles%h_xvelo_old, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:yvelo_old', rParticles%npart, ST_DOUBLE, rParticles%h_yvelo_old, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:zvelo_old', rParticles%npart, ST_DOUBLE, rParticles%h_zvelo_old, &
                            ST_NEWBLOCK_NOINIT)

  call storage_new ('euler_lagrange', 'Particle:xvelo_gas', rParticles%npart, ST_DOUBLE, rParticles%h_xvelo_gas, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:yvelo_gas', rParticles%npart, ST_DOUBLE, rParticles%h_yvelo_gas, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:zvelo_gas', rParticles%npart, ST_DOUBLE, rParticles%h_zvelo_gas, &
                            ST_NEWBLOCK_NOINIT)

  call storage_new ('euler_lagrange', 'Particle:xvelo_gas_old', rParticles%npart, ST_DOUBLE, rParticles%h_xvelo_gas_old, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:yvelo_gas_old', rParticles%npart, ST_DOUBLE, rParticles%h_yvelo_gas_old, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:zvelo_gas_old', rParticles%npart, ST_DOUBLE, rParticles%h_zvelo_gas_old, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:bdy_time', rParticles%npart, ST_DOUBLE, rParticles%h_bdy_time, &
                            ST_NEWBLOCK_NOINIT)
  call storage_new ('euler_lagrange', 'Particle:bdy_check', rParticles%npart, ST_DOUBLE, rParticles%h_bdy_check, &
                            ST_NEWBLOCK_NOINIT)

   md_el_length(1)=2
   md_el_length(2)=p_rtriangulation%NEL

   ! midpoints of the elements
   call storage_new ('euler_lagrange', 'Elements:midpoints', md_el_length, ST_DOUBLE, rParticles%h_midpoints_el, &
                            ST_NEWBLOCK_NOINIT)
   ! volumepart of the particles
   call storage_new ('euler_lagrange', 'Elements:particlevolume', p_rtriangulation%NEL, ST_DOUBLE, rParticles%h_PartVol, &
                            ST_NEWBLOCK_NOINIT)
   ! velocity of the particles
   call storage_new ('euler_lagrange', 'Elements:particlevelocity', md_el_length, ST_DOUBLE, rParticles%h_PartVelo, &
                            ST_NEWBLOCK_NOINIT)
   
   call storage_getbase_double (rParticles%h_xpos, rParticles%p_xpos)
   call storage_getbase_double (rParticles%h_ypos, rParticles%p_ypos)
   call storage_getbase_double (rParticles%h_zpos, rParticles%p_zpos)
   call storage_getbase_double (rParticles%h_xpos_old, rParticles%p_xpos_old)
   call storage_getbase_double (rParticles%h_ypos_old, rParticles%p_ypos_old)
   call storage_getbase_double (rParticles%h_zpos_old, rParticles%p_zpos_old)
   call storage_getbase_double (rParticles%h_xvelo, rParticles%p_xvelo)
   call storage_getbase_double (rParticles%h_yvelo, rParticles%p_yvelo)
   call storage_getbase_double (rParticles%h_zvelo, rParticles%p_zvelo)
   call storage_getbase_double (rParticles%h_xvelo_old, rParticles%p_xvelo_old)
   call storage_getbase_double (rParticles%h_yvelo_old, rParticles%p_yvelo_old)
   call storage_getbase_double (rParticles%h_zvelo_old, rParticles%p_zvelo_old)
   call storage_getbase_double (rParticles%h_xvelo_gas, rParticles%p_xvelo_gas)
   call storage_getbase_double (rParticles%h_yvelo_gas, rParticles%p_yvelo_gas)
   call storage_getbase_double (rParticles%h_zvelo_gas, rParticles%p_zvelo_gas)
   call storage_getbase_double (rParticles%h_xvelo_gas_old, rParticles%p_xvelo_gas_old)
   call storage_getbase_double (rParticles%h_yvelo_gas_old, rParticles%p_yvelo_gas_old)
   call storage_getbase_double (rParticles%h_zvelo_gas_old, rParticles%p_zvelo_gas_old)
   call storage_getbase_int (rParticles%h_element, rParticles%p_element)
   call storage_getbase_double (rParticles%h_lambda1, rParticles%p_lambda1)
   call storage_getbase_double (rParticles%h_lambda2, rParticles%p_lambda2)
   call storage_getbase_double (rParticles%h_lambda3, rParticles%p_lambda3)
   call storage_getbase_double (rParticles%h_lambda4, rParticles%p_lambda4)
   call storage_getbase_double (rParticles%h_diam, rParticles%p_diam)
   call storage_getbase_double (rParticles%h_mass, rParticles%p_mass)
   call storage_getbase_double (rParticles%h_alpha_n, rParticles%p_alpha_n)
   call storage_getbase_double (rParticles%h_bdy_time, rParticles%p_bdy_time)
   call storage_getbase_double (rParticles%h_bdy_check, rParticles%p_bdy_check)

   call storage_getbase_double (rParticles%h_PartVol, rParticles%p_PartVol)
   call storage_getbase_double2D (rParticles%h_PartVelo, rParticles%p_PartVelo)
   call storage_getbase_double2D (rParticles%h_midpoints_el, rParticles%p_midpoints_el)
  
   ! Set pointer to coordinate vector
   call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
 
   ! Set pointer to vertices at element
   call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

   ! Stores the midpoint for each element
   do iel=1,p_rtriangulation%NEL

      rParticles%p_midpoints_el(1,iel)= &
                                (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
                                 p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
                                 p_DvertexCoords(1,p_IverticesAtElement(3,iel)))/3.0_dp

      rParticles%p_midpoints_el(2,iel)= &
                                (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
                                 p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
                                 p_DvertexCoords(2,p_IverticesAtElement(3,iel)))/3.0_dp

    end do

    ! get values for the startingpositions of the particles
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "xmin", partxmin)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "xmax", partxmax)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "ymin", partymin)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "ymax", partymax)
 
    ! get particlevelocity
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "velopartx", velopartx)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "veloparty", veloparty)

    ! get particle-mass and diameter
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlemass", particlemass)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "particlediam", particlediam)

    ! get values for gravity
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "gravityx", gravityx)
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "gravityy", gravityy)

    ! get value of kinematic viscosity
    call parlst_getvalue_double(rparlist, 'Eulerlagrange', "gas_nu", gas_nu)

    ! get boundarybehaviour
    call parlst_getvalue_int(rparlist, 'Eulerlagrange', "boundbehav", boundbehav)

    ! store particlesnumber, viscosity of the gas and gravity  
    rParticles%npart = nPart
    rParticles%nu_g= gas_nu
    rParticles%gravity(1)= gravityx
    rParticles%gravity(2)= gravityy
    rParticles%iTimestep= 0
    rParticles%maxvalx= maxval(p_DvertexCoords(1,:))
    rParticles%p_PartVol= 0
    rParticles%p_PartVelo= 0

    ! set boundaryconditions for the particles
    select case(boundbehav)
    case (0)
        rParticles%tang_val= 1.0d0
        rParticles%norm_val= 1.0d0
    case (1)
        rParticles%tang_val= 1.0d0
        rParticles%norm_val= 0.0d0
    case default
      call output_line('Invalid boundaryconditions!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_boundbehav')
      call sys_halt()
    end select

    ! initialize data for each particle
    do iPart=1,rParticles%npart
  
  		!Hole Zufallszahl
		call random_number(random1)
		call random_number(random2)
		
        ! set initial values for the particles
        rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
        rParticles%p_ypos(iPart)= partymin + random2*(partymax - partymin)
        rParticles%p_xpos_old(iPart)= partxmin + random1*(partxmax - partxmin)
        rParticles%p_ypos_old(iPart)= partymin + random2*(partymax - partymin)
        rParticles%p_xvelo(iPart)= velopartx
        rParticles%p_yvelo(iPart)= veloparty
        rParticles%p_xvelo_old(iPart)= velopartx
        rParticles%p_yvelo_old(iPart)= veloparty
        rParticles%p_xvelo_gas(iPart)= 0d0
        rParticles%p_yvelo_gas(iPart)= 0d0
        rParticles%p_xvelo_gas_old(iPart)= 0d0
        rParticles%p_yvelo_gas_old(iPart)= 0d0
        rParticles%p_diam(iPart)= particlediam
        rParticles%p_mass(iPart)= particlemass
        rParticles%p_alpha_n(iPart)= 0
        rParticles%p_element(iPart)= 1
        rParticles%p_bdy_time(iPart)= 0
        rParticles%p_bdy_check(iPart)= 0
        
        ! Find the start element for each particle
        call findnewelement(rparlist,p_rproblemLevel,rParticles,iPart)

        ! calculate barycentric coordinates
        call calculatebarycoords(p_rproblemLevel,rParticles,iPart)

        ! wrong element
        if ((abs(rParticles%p_lambda1(iPart))+abs(rParticles%p_lambda2(iPart))+&
                  abs(rParticles%p_lambda3(iPart))-1) .GE. 0.00001) then
            call wrongelement(rparlist,p_rproblemLevel,rParticles,iPart)
        end if

    end do


end subroutine eulerlagrange_init

subroutine eulerlagrange_step(rparlist,p_rproblemLevel,rsolution,rtimestep,rcollection,rParticles)

    ! parameterlist
    type(t_parlist), intent(inout) :: rparlist

    ! collection structure
    type(t_collection), intent(inout) :: rcollection

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! particles
    type(t_Particles), intent(inout) :: rParticles

    ! primal solution vector
    type(t_vectorBlock), intent(inout), target :: rsolution

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation
 
    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement
 
    
    ! current particlenumber
    integer :: iPart
    character(LEN=20) :: sfilename, sfilenamenew

    real(DP) :: dx,dy
    real(DP), dimension(2,4) :: DcornerCoords

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation

    ! Set pointer to coordinate vector
    call storage_getbase_double2D(&
         p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
 
    ! Set pointer to vertices at element
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)


    write(sfilename,'(i0)') rParticles%iTimestep
    sfilenamenew='particleflow'//trim(sfilename)//'.vtk'

    OPEN(20+rParticles%iTimestep,file='out/video/'//sfilenamenew)

    do iPart = 1, rParticles%nPart
      ! subroutine to find the element with the particle
      call findnewelement(rparlist,p_rproblemLevel,rParticles,iPart)
 
      dx= rParticles%p_xpos(iPart)
      dy= rParticles%p_ypos(iPart)
      
      DcornerCoords(1,1)= p_DvertexCoords(1,p_IverticesAtElement(1,rParticles%p_element(iPart)))
      DcornerCoords(1,2)= p_DvertexCoords(1,p_IverticesAtElement(2,rParticles%p_element(iPart)))
      DcornerCoords(1,3)= p_DvertexCoords(1,p_IverticesAtElement(3,rParticles%p_element(iPart)))
      DcornerCoords(2,1)= p_DvertexCoords(2,p_IverticesAtElement(1,rParticles%p_element(iPart)))
      DcornerCoords(2,2)= p_DvertexCoords(2,p_IverticesAtElement(2,rParticles%p_element(iPart)))
      DcornerCoords(2,3)= p_DvertexCoords(2,p_IverticesAtElement(3,rParticles%p_element(iPart)))
 
      ! get barycentric coordinates
      call gaux_getBarycentricCoords_tri2D (DcornerCoords,dx,dy,&
            rParticles%p_lambda1(iPart),rParticles%p_lambda2(iPart),rParticles%p_lambda3(iPart))
                   
      ! check if particle is in the wrong element
      IF ((abs(rParticles%p_lambda1(iPart))+&
           abs(rParticles%p_lambda2(iPart))+&
           abs(rParticles%p_lambda3(iPart))) .GE. 1.00001) then
        call wrongelement(rparlist,p_rproblemLevel,rParticles,iPart)
      end if

      write(20+rParticles%iTimestep,*) rParticles%p_xpos(iPart), rParticles%p_ypos(iPart)

    end do

    ! subroutine to compute the new position of the particles
    call moveparticle(rparlist,p_rproblemLevel,rsolution,rParticles)

    ! subroutine to calculate the volume part of the particles
    call calculatevolumepart(p_rproblemLevel,rParticles)

    ! subroutine to calculate the velocity of the particles 
    call calculatevelopart(p_rproblemLevel,rParticles)

    close(unit=20+rParticles%iTimestep)
    rParticles%iTimestep=rParticles%iTimestep+1

end subroutine eulerlagrange_step


!************ SUBROUTINE to calculate the barycentric coordinates ******************************************
!*

subroutine calculatebarycoords(p_rproblemLevel,rParticles,iPart)

    ! particles
    type(t_Particles), intent(inout) :: rParticles
    
    ! current number of particle
    integer, intent(inout) :: iPart

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! pointer to vertices at each element
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! pointer to the vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! current element
    integer :: currentElement

    ! coordinates of the vertices of the actual element
    real(DP), dimension(2,3) :: vert_coord

    ! determinates
    real(DP) :: det_A, det_A1, det_A2, det_A3

    ! local variables
    integer :: ivt

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
  
    ! Get vertices at element
    call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Get coordinates of the vertices
    call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! store current element
    currentElement = rParticles%p_element(iPart)

    ! store coordinates of the vertices
    do ivt=1,3

      vert_coord(1,ivt)= p_DvertexCoords(1,p_IverticesAtElement(ivt,currentElement))
      vert_coord(2,ivt)= p_DvertexCoords(2,p_IverticesAtElement(ivt,currentElement))

    end do

    ! compute determinate
    !	    1 x1 y1
    !    A=	1 x2 y2
    !		1 x3 y3
    ! detA=x2y3+x1y2+x3y1-x2y1-x3y2-x1y3
    det_A=  vert_coord(1,2) * vert_coord(2,3)+&
            vert_coord(1,1) * vert_coord(2,2)+&
            vert_coord(1,3) * vert_coord(2,1)-&
            vert_coord(1,2) * vert_coord(2,1)-&
            vert_coord(1,3) * vert_coord(2,2)-&
            vert_coord(1,1) * vert_coord(2,3)

    ! calculate barycentric coorinates (lambda1,lambda2,lambda3)
    ! lambda1
    !		1 x  y
    !   A1=	1 x2 y2
    !		1 x3 y3
    ! detA1=x2y3+xy2+x3y-x2y-x3y2-xy3
    det_A1= vert_coord(1,2)         * vert_coord(2,3)+&
            rParticles%p_xpos(iPart)* vert_coord(2,2)+&
            vert_coord(1,3)         * rParticles%p_ypos(iPart)-&
            vert_coord(1,2)         * rParticles%p_ypos(iPart)-&
            vert_coord(1,3)         * vert_coord(2,2)-&
            rParticles%p_xpos(iPart)* vert_coord(2,3)
 
    !lambda1=|det_A1/det_A|
    rParticles%p_lambda1(iPart)= abs(det_A1/det_A)

    ! lambda2
    !		1 x1 y1
    !   A2=	1 x  y
    !		1 x3 y3
    ! detA2=xy3+x1y+x3y1-xy1-x3y-x1y3
    det_A2= rParticles%p_xpos(iPart)* vert_coord(2,3)+&
            vert_coord(1,1)         * rParticles%p_ypos(iPart)+&
            vert_coord(1,3)         * vert_coord(2,1)-&
            rParticles%p_xpos(iPart)* vert_coord(2,1)-&
            vert_coord(1,3)         * rParticles%p_ypos(iPart)-&
            vert_coord(1,1)         * vert_coord(2,3)

    !lambda2=|det_A2/det_A|
    rParticles%p_lambda2(iPart)= abs(det_A2/det_A)

    ! lambda3
    !		1 x1 y1
    !   A3=	1 x2 y2
    !		1 x  y
    ! detA3=x2y+x1y2+xy1-x2y1-xy2-x1y
    det_A3= vert_coord(1,2)         * rParticles%p_ypos(iPart)+&
            vert_coord(1,1)         * vert_coord(2,2)+&
            rParticles%p_xpos(iPart)* vert_coord(2,1)-&
            vert_coord(1,2)         * vert_coord(2,1)-&
            rParticles%p_xpos(iPart)* vert_coord(2,2)-&
            vert_coord(1,1)         * rParticles%p_ypos(iPart)

    ! lambda3=|det_A3/det_A|
    rParticles%p_lambda3(iPart)= abs(det_A3/det_A)


end subroutine calculatebarycoords

    !*
    !**********************************************************************************************************************


    !******************************* SUBROUTINE to find the right element *************************************************
    !*

subroutine findnewelement(rparlist,p_rproblemLevel,rParticles,iPart)

    ! parameterlist
    type(t_parlist), intent(inout) :: rparlist
    
    ! particles
    type(t_Particles), intent(inout) :: rParticles

    ! pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation


    ! pointer to the neighbour elements adjacent to an element
    !
    ! Handle to 
    !       p_IneighboursAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! For each element, the numbers of adjacent elements
    ! in mathematically positive sense, meeting the element in an edge.
    ! p_RneighbourElement(IEL)\%Ineighbours(.) describes the elements adjacent 
    ! to IEL along the edges (p_RedgesOnElement(IEL)\%Iedges(.)-NVT).
    ! This is the old KADJ array.
    !
    ! Note:  For meshes with hanging vertices, this array is slightly
    ! modified. For 'big' elements this array contains the element
    ! numbers of the 'first' adjacent element via an edge/face.
    ! Note: To access all elements adjacent to an element via a
    ! hanging vertex, calculate the vertex number of the hanging
    ! vertex via InodalProperty and access all adjacent elements.
    ! For small 'hanging' elements, this array contains as usual
    ! the number of the 'big' adjacent element(s).
    integer, dimension(:,:), pointer :: p_IneighboursAtElement


    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

   
    ! current number of particle
    integer, intent(inout) :: iPart
  
    ! position of the element
    real(DP), dimension(2) :: particlepos

    ! element number
    integer :: iel
 
    ! particle is in element 
    logical :: binside
    
    ! search mode
    character(LEN=15) :: searchmode

    ! variables for midpoints_el
	real(DP), dimension(1:5) :: distances
	integer :: i, adj, minl, ite, mittelelement
	integer, parameter :: itemax = 10000
	real(DP) :: distToMid

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
   
    ! get vertices at element
    call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! get coordinates for elements
    call storage_getbase_double2d (&
        p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
        
    ! get neighboured elements
    call storage_getbase_int2D(&
        p_rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    
    ! set particles position
    particlepos(1)= rParticles%p_xpos(iPart)
    particlepos(2)= rParticles%p_ypos(iPart) 
   
    ! Get searchmode (brute force, raytrace, etc.)  
    call parlst_getvalue_string(rparlist, 'Eulerlagrange', "search", searchmode)

    select case(trim(searchmode))
    case('bruteforce')
        call tsrch_getElem_BruteForce(particlepos,p_DvertexCoords,p_IverticesAtElement,iel)
    case('raytrace2D')
        !call tsrch_getElem_raytrace2D(&
        !        Dpoint,rtriangulation,iel, iresult,ilastElement,ilastEdge,imaxIterations)
    case('midpoint')

	    distToMid = 10000.0d0

	    gotoNextElm: do ite = 1, itemax
	    
		    ! Calculate the distances to the midpoints
		    distances = 10000.0d0
   
		    do i = 1, 3 
		      if (p_IneighboursAtElement(i,rParticles%p_element(iPart)) > 0) then
			    adj = p_IneighboursAtElement(i,rParticles%p_element(iPart))
			    distances(i+1) = (rParticles%p_xpos(iPart)-rParticles%p_midpoints_el(1,adj))**2.0d0 +&
								 (rParticles%p_ypos(iPart)-rParticles%p_midpoints_el(2,adj))**2.0d0
			  end if
		    end do
		    
		    distances(1) = (rParticles%p_xpos(iPart)-&
		                        rParticles%p_midpoints_el(1,rParticles%p_element(iPart)))**2.0d0+&
		                   (rParticles%p_ypos(iPart)-&
		                        rParticles%p_midpoints_el(2,rParticles%p_element(iPart)))**2.0d0
		
		    ! Position with the smallest distance
		    minl = minloc(distances,1)

		    ! Check if the distance didn't change
		    if (minl .EQ. 1) exit gotoNextElm

		    ! Store element with the lowest distance
		    rParticles%p_element(iPart) = p_IneighboursAtElement(minl-1,rParticles%p_element(iPart))
		    distToMid = distances(minl)
		    
	    end do gotoNextElm
	    
    case default
        call output_line('Invalid search mode!',&
                     OU_CLASS_WARNING,OU_MODE_STD,'flagship')
        call sys_halt()
    end select

end subroutine findnewelement

!*
!**********************************************************************************************************************


!****************************************** SUBROUTINE for wrong elements *********************************************
!*

subroutine wrongelement(rparlist,p_rproblemLevel,rParticles,iPart)

    ! parameterlist
    type(t_parlist), intent(inout) :: rparlist

    ! particles
    type(t_Particles), intent(inout) :: rParticles

    ! current number of particle
    integer, intent(inout) :: iPart

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to array containing the elements adjacent to a vertex.
    !
    ! Handle to 
    !       p_IelementsAtVertex = array(1..*) of integer
    ! p_IelementsAtVertex ( p_IelementsAtVertexIdx(IVT)..p_IelementsAtVertexIdx(IVT+1)-1 )
    ! contains the number of the adjacent element in a vertex.
    ! This replaces the old KVEL array.
    !
    ! Note: For hanging vertices, this array contains only those
    ! elements which are 'corner adjacent' to a vertex (i.e. the 'smaller' elements).
    ! The 'big' elements adjacent to the edge which the hanging vertex
    ! is a midpoint of are not part of the vertex neighbourhood
    ! in this array.
    integer, dimension(:), pointer :: p_IelementsAtVertex

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Handle to 
    !       p_IelementsAtVertexIdx=array [1..NVT+1] of integer.
    ! Index array for p_IelementsAtVertex of length NVT+1 for describing the
    ! elements adjacent to a corner vertex. for vertex IVT, the array
    ! p_IelementsAtVertex contains the numbers of the elements around this
    ! vertex at indices 
    !     p_IelementsAtVertexIdx(IVT)..p_IelementsAtVertexIdx(IVT+1)-1.
    ! By subtracting
    !     p_IelementsAtVertexIdx(IVT+1)-p_IelementsAtVertexIdx(IVT)
    ! One can get the number of elements adjacent to a vertex IVT.
    integer, dimension(:), pointer ::  p_IelementsAtVertexIdx

    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

 	! local variables
 	integer :: Vert, Elm, currentelm, nVertex, Neighbour, i_NVBD, lexipunkt, i
    real(DP) :: dxi1,dxi2,dxi3
    real(DP) :: dx,dy
    real(DP), dimension(2,4) :: DcornerCoords

    logical :: binside

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
  
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
         
    call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex,&
        p_IelementsAtVertex)
        
    call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)

    call storage_getbase_double2d (&
        p_rtriangulation%h_DvertexCoords,p_DvertexCoords)


    ! store the current element
	currentelm = rParticles%p_element(iPart)
	dxi1=0d0
	dxi2=0d0
	dxi3=0d0
	dx= rParticles%p_xpos(iPart)
	dy= rParticles%p_ypos(iPart)
	
    ! Loop over the vertices of the element
	SearchVertex: do Vert = 1, 3												
		
		nVertex = p_IverticesAtElement(Vert, currentelm)
			
		! Loop over the element containing to the vertex
		SearchElement: do Elm = 1, (p_IelementsAtVertexIdx(nVertex+1)-p_IelementsAtVertexIdx(nVertex))		
											
	    	if (p_IelementsAtVertex(p_IelementsAtVertexIdx(nVertex)+Elm-1) == 0) then
				exit SearchElement
			end if

			rParticles%p_element(iPart) = p_IelementsAtVertex(p_IelementsAtVertexIdx(nVertex)+Elm-1)

            ! store coordinates of cornervertices
            DcornerCoords(1,1)= p_DvertexCoords(1,p_IverticesAtElement(1,rParticles%p_element(iPart)))
            DcornerCoords(1,2)= p_DvertexCoords(1,p_IverticesAtElement(2,rParticles%p_element(iPart)))
            DcornerCoords(1,3)= p_DvertexCoords(1,p_IverticesAtElement(3,rParticles%p_element(iPart)))
            DcornerCoords(2,1)= p_DvertexCoords(2,p_IverticesAtElement(1,rParticles%p_element(iPart)))
            DcornerCoords(2,2)= p_DvertexCoords(2,p_IverticesAtElement(2,rParticles%p_element(iPart)))
            DcornerCoords(2,3)= p_DvertexCoords(2,p_IverticesAtElement(3,rParticles%p_element(iPart)))

            ! check if the particle is in the element
            call gaux_isInElement_tri2D(dx,dy,DcornerCoords,binside)
          
            if (binside .eq. .true.) then
                exit SearchVertex
            end if 
            
		end do SearchElement

	end do SearchVertex

    ! get barycentric coordinates
    call gaux_getBarycentricCoords_tri2D (DcornerCoords,dx,dy,dxi1,dxi2,dxi3)

    if ((abs(dxi1)+abs(dxi2)+abs(dxi3)) .GE. 1.00001) then
	  rParticles%p_element(iPart) = currentelm
	  call checkboundary(rparlist,p_rproblemLevel,rParticles,iPart)
    end if


end subroutine wrongelement

!*
!**********************************************************************************************************************

!*************************************** SUBROUTINE to move the particle **********************************************
!*

subroutine moveparticle(rparlist,p_rproblemLevel,rsolutionPrimal,rParticles)

    ! parameterlist
    type(t_parlist), intent(inout) :: rparlist

    ! particles
    type(t_Particles), intent(inout) :: rParticles

     ! primal solution vector
    type(t_vectorBlock), intent(inout) :: rsolutionPrimal

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    ! pointer to elements adjacent to the boundary. 
    !
    ! Handle to 
    !       p_IelementsAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all elements on the (real) boundary
    ! in mathematically positive sense.
    ! p_IelementsAtBoundary(i) is the element adjacent to edge
    ! h_IedgesAtBoundary - therefore one element number might appear
    ! more than once in this array!
    ! The boundary elements of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KEBD array.
    integer, dimension(:), pointer :: p_IelementsAtBoundary


    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

	real(DP) :: ux1_part, uy1_part, ux2_part, uy2_part, ux3_part, uy3_part
	real(DP), dimension(3) :: rho_gas
    real(DP) :: rho_g, C_W, Re_p, Velo_rel, dt, c_pi
    
    type(t_vectorScalar) :: rvector1, rvector2, rvector3
    real(DP), dimension(:), pointer :: p_Ddata1, p_Ddata2, p_Ddata3
    integer :: iPart

    ! startingpostions of the particles
    real(DP) :: partxmin, partxmax, partymin, partymax

    ! velocity of the particles
    real(DP) :: velopartx, veloparty, random1, random2

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
 
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    rho_g= 0d0 
    C_W=0d0 
    Re_p=0d0 
    Velo_rel=0d0
    dt=0d0
    c_pi=0d0
    
    ! get values for the startingpositions of the particles
    call parlst_getvalue_double(rparlist, 'Timestepping', "dinitialStep", dt)

    ! Get data from solution
    call eulerlagrange_getVariable(rsolutionPrimal, 'velocity_x', rvector1)
    call eulerlagrange_getVariable(rsolutionPrimal, 'velocity_y', rvector2)
    call eulerlagrange_getVariable(rsolutionPrimal, 'density', rvector3)
    call lsyssc_getbase_double(rvector1, p_Ddata1)
    call lsyssc_getbase_double(rvector2, p_Ddata2)
    call lsyssc_getbase_double(rvector3, p_Ddata3)
 
    c_pi= 3.14159265358979323846264338327950288

    ! Loop over the particles
    do iPart = 1, rParticles%nPart

	! store old data
	rParticles%p_xpos_old(iPart)=	   rParticles%p_xpos(iPart)
	rParticles%p_ypos_old(iPart)=	   rParticles%p_ypos(iPart)
	rParticles%p_xvelo_old(iPart)=	   rParticles%p_xvelo(iPart)
	rParticles%p_yvelo_old(iPart)=	   rParticles%p_yvelo(iPart)
	rParticles%p_xvelo_gas_old(iPart)= rParticles%p_xvelo_gas(iPart)
	rParticles%p_yvelo_gas_old(iPart)= rParticles%p_xvelo_gas(iPart)

	! velocity and density of the gas in the first corner (in mathematically positive sense)
	ux1_part= p_Ddata1(p_IverticesAtElement(1,rParticles%p_element(iPart)))
	uy1_part= p_Ddata2(p_IverticesAtElement(1,rParticles%p_element(iPart)))
	rho_gas(1)= p_Ddata3(p_IverticesAtElement(1,rParticles%p_element(iPart)))

	! velocity and density of the gas in the second corner (in mathematically positive sense)
	ux2_part= p_Ddata1(p_IverticesAtElement(2,rParticles%p_element(iPart)))
	uy2_part= p_Ddata2(p_IverticesAtElement(2,rParticles%p_element(iPart)))
	rho_gas(2)= p_Ddata3(p_IverticesAtElement(2,rParticles%p_element(iPart)))

	! velocity and density of the gas in the third corner (in mathematically positive sense)
	ux3_part= p_Ddata1(p_IverticesAtElement(3,rParticles%p_element(iPart)))
	uy3_part= p_Ddata2(p_IverticesAtElement(3,rParticles%p_element(iPart)))
	rho_gas(3)= p_Ddata3(p_IverticesAtElement(3,rParticles%p_element(iPart)))

	! calculate velocity of the gas
	rParticles%p_xvelo_gas(iPart)= 	rParticles%p_lambda1(iPart)*ux1_part + &
									rParticles%p_lambda2(iPart)*ux2_part + &
									rParticles%p_lambda3(iPart)*ux3_part 
	rParticles%p_yvelo_gas(iPart)= 	rParticles%p_lambda1(iPart)*uy1_part + &
									rParticles%p_lambda2(iPart)*uy2_part + &
									rParticles%p_lambda3(iPart)*uy3_part

	! calculate the density of the gas in the position of the particle
	rho_g= 	rParticles%p_lambda1(iPart)*rho_gas(1) + rParticles%p_lambda2(iPart)*&
	        rho_gas(2) + rParticles%p_lambda3(iPart)*rho_gas(3) 


	! calculate particle Reynoldsnumber
	!
	! Re_p= \frac{d_p\ \left|\textbf{u}_g-\textbf{u}_p\right|}{\nu_g}
	! with \nu_g=\frac{\eta_g}{\rho_g}
	!
	Re_p= (rho_g*0.5*rParticles%p_diam(iPart)/rParticles%nu_g)*(sqrt((rParticles%p_xvelo_gas(iPart)- &
	       rParticles%p_xvelo_old(iPart))**2+&
		  (rParticles%p_yvelo_gas(iPart)-rParticles%p_yvelo_old(iPart))**2))

	! calculate the drag force coefficient
	if (Re_p<1000) then
		C_W= 24/Re_p*(1+0.15*Re_p**0.687)
	else
		C_W= 24/Re_p
	end if

	! calculate alpha_n
	rParticles%p_alpha_n(iPart)= C_W*c_pi*rho_g/8 

	! calculate the velocity
	Velo_rel= sqrt((rParticles%p_xvelo_old(iPart)-rParticles%p_xvelo_gas(iPart))**2 +&
	               (rParticles%p_yvelo_old(iPart)-rParticles%p_yvelo_gas(iPart))**2)

	! calculate new velocity of the particle
	rParticles%p_xvelo(iPart)= 	(rParticles%p_mass(iPart) * rParticles%p_xvelo_old(iPart)+&
								dt*rParticles%p_alpha_n(iPart) * Velo_rel * 0.25 * rParticles%p_diam(iPart)**2 &
								* rParticles%p_xvelo_gas(iPart)+ dt*rParticles%p_mass(iPart) * rParticles%gravity(1))/&
								(rParticles%p_mass(iPart) + dt*rParticles%p_alpha_n(iPart)*Velo_rel*&
								0.25*rParticles%p_diam(iPart)**2)
	rParticles%p_yvelo(iPart)= 	(rParticles%p_mass(iPart) * rParticles%p_yvelo_old(iPart)+&
								dt*rParticles%p_alpha_n(iPart)*Velo_rel*0.25*rParticles%p_diam(iPart)**2 &
								* rParticles%p_yvelo_gas(iPart)+ dt*rParticles%p_mass(iPart)*rParticles%gravity(2))/&
    							(rParticles%p_mass(iPart) + dt*rParticles%p_alpha_n(iPart)*Velo_rel*&
    							0.25*rParticles%p_diam(iPart)**2)

	!---------------------------------------------------------------------------------
	! Calculate the new position of the particle
    !
	! x_new= x_old + delta t * ux_old
	!---------------------------------------------------------------------------------
	
	rParticles%p_xpos(iPart) = rParticles%p_xpos_old(iPart) + dt * rParticles%p_xvelo(iPart)
	rParticles%p_ypos(iPart) = rParticles%p_ypos_old(iPart) + dt * rParticles%p_yvelo(iPart)
	
	if (rParticles%p_xpos(iPart).ge.rParticles%maxvalx) then
	
	    ! get values for the startingpositions of the particles
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "xmin", partxmin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "xmax", partxmax)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "ymin", partymin)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "ymax", partymax)
 
        ! get particlevelocity
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "velopartx", velopartx)
        call parlst_getvalue_double(rparlist, 'Eulerlagrange', "veloparty", veloparty)
	
	  	!Hole Zufallszahl
		call random_number(random1)
		call random_number(random2)
		
        ! set initial values for the particles
        rParticles%p_xpos(iPart)= partxmin + random1*(partxmax - partxmin)
        rParticles%p_ypos(iPart)= partymin + random2*(partymax - partymin)
        rParticles%p_xpos_old(iPart)= partxmin + random1*(partxmax - partxmin)
        rParticles%p_ypos_old(iPart)= partymin + random2*(partymax - partymin)
        rParticles%p_xvelo(iPart)= velopartx
        rParticles%p_yvelo(iPart)= veloparty
        rParticles%p_xvelo_old(iPart)= velopartx
        rParticles%p_yvelo_old(iPart)= veloparty
        rParticles%p_xvelo_gas(iPart)= 0d0
        rParticles%p_yvelo_gas(iPart)= 0d0
        rParticles%p_xvelo_gas_old(iPart)= 0d0
        rParticles%p_yvelo_gas_old(iPart)= 0d0
        rParticles%p_alpha_n(iPart)= 0
        rParticles%p_element(iPart)= 1
        rParticles%p_bdy_time(iPart)= 0
        rParticles%p_bdy_check(iPart)= 0
        
        ! Find the start element for each particle
        call findnewelement(rparlist,p_rproblemLevel,rParticles,iPart)

        ! calculate barycentric coordinates
        call calculatebarycoords(p_rproblemLevel,rParticles,iPart)

        ! wrong element
        if ((abs(rParticles%p_lambda1(iPart))+abs(rParticles%p_lambda2(iPart))+&
                  abs(rParticles%p_lambda3(iPart))-1) .GE. 0.00001) then
            call wrongelement(rparlist,p_rproblemLevel,rParticles,iPart)
        end if

	end if	
	
	end do
	
    ! release temporal data
    call lsyssc_releasevector(rvector1)
    call lsyssc_releasevector(rvector2)
    call lsyssc_releasevector(rvector3)

end subroutine moveparticle

!*
!**********************************************************************************************************************

!****************************************** SUBROUTINE to check the boundary ******************************************
!*

subroutine checkboundary(rparlist,p_rproblemLevel,rParticles,iPart)

    ! parameterlist
    type(t_parlist), intent(inout) :: rparlist

    ! particles
    type(t_Particles), intent(inout) :: rParticles

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel
    
    ! current number of particle
    integer, intent(inout) :: iPart

    
    ! pointer to vertices on boundary. 
    !
    ! Handle to 
    !       p_IverticesAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all vertices on the (real) boundary
    ! in mathematically positive sense.
    ! The boundary vertices of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KVBD array.
    integer, dimension(:), pointer :: p_IverticesAtBoundary

    ! pointer to edges adjacent to the boundary. 
    !
    ! Handle to 
    !       p_IedgesAtBoundary = array [1..NMBD] of integer.
    ! This array contains a list of all edges on the (real) boundary.
    ! 2D: in mathematically positive sense. 
    ! 3D: with increasing number.
    ! The boundary edges of boundary component i are saved at
    !        p_IboundaryCpEdgesIdx(i)..p_IboundaryCpEdgesIdx(i+1)-1.
    ! This is the old KMBD array.
    ! (Note: In 2D, the above index pointer coincides with
    !        p_IboundaryCpEdgesIdx(i)..p_IboundaryCpEdgesIdx(i+1)-1 ).
    integer, dimension(:), pointer :: p_IedgesAtBoundary

    ! pointer to elements adjacent to the boundary. 
    !
    ! Handle to 
    !       p_IelementsAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all elements on the (real) boundary
    ! in mathematically positive sense.
    ! p_IelementsAtBoundary(i) is the element adjacent to edge
    ! h_IedgesAtBoundary - therefore one element number might appear
    ! more than once in this array!
    ! The boundary elements of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KEBD array.
    integer, dimension(:), pointer :: p_IelementsAtBoundary

    ! pointer to the coordinates of the vertices
    !
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! local variables	
	integer :: i_NVBD, current, i, j, ij, jj
	real(DP), dimension(4) :: x,y
	real(DP) :: s
	real(DP), dimension(2) :: velo_rest
	integer, dimension(2) :: bdy_first
	integer :: change_bdy
	real(DP) :: tang_norm, proj_tang, proj_norm
	real(DP), dimension(2) :: bdy_tang, bdy_norm, bdy_move, bdy_point
	integer, dimension(2) :: bdy_Koords

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation

    call storage_getbase_int(&
         p_rtriangulation%h_IelementsAtBoundary, p_IelementsAtBoundary)
    call storage_getbase_double2D(&
         p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int(&
         p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

	x= 0
	y= 0
	j= 1
	ij= 1
	jj= 1
	tang_norm=0d0
	proj_tang=0d0
	proj_norm=0d0
	s=0d0

	!Speichern des aktuellen Elements
	current = rParticles%p_element(iPart)

	!Schleife über alle Randelemente
	Boundary_Search: do i_NVBD = 1, p_rtriangulation%NVBD

		!Wenn das aktuelle Element ein Randelement ist:
		if (current == p_IelementsAtBoundary(i_NVBD)) then

			!Bei dem Element handelt es sich um ein Randelement
			rParticles%p_bdy_check(iPart) = 1

			!Suche die beiden Eckpunkte, die auf dem Rand liegen
			e_search: do j=1,3
				x_Search: do i = 1, p_rtriangulation%NVBD
					if (p_IverticesAtBoundary(i)==p_IverticesAtElement(j,current)) then
						bdy_first(jj)= i
						jj=jj+1
						bdy_Koords(ij)=p_IverticesAtBoundary(i)
						ij= ij+1
						if (ij==3) exit e_Search
					end if
				end do x_search
			end do e_search

			!So steht in bdy_Koord(1) der erste Knotenpunkt (gegen den Uhrzeigersinn)
			if (bdy_first(1) .GE. bdy_first(2)) then
				change_bdy = bdy_Koords(1)
				bdy_Koords(1)= bdy_Koords(2)
				bdy_Koords(2)= change_bdy
			end if

			!x-Koordinate des ersten Randpunktes (gegen den Uhrzeigersinn)
			x(1)= p_DvertexCoords(1,bdy_Koords(1))
			!x-Koordinate des Strecke vom ersten zum zweiten Randpunkt
			x(2)= p_DvertexCoords(1,bdy_Koords(2)) - p_DvertexCoords(1,bdy_Koords(1)) 
			!x-Koordinate des Strecke vom ersten zum zweiten Randpunkt
			x(3)= rParticles%p_xpos_old(iPart)
			!x-Koordinate des Strecke vom ersten zum zweiten Randpunkt
			x(4)= rParticles%p_xpos(iPart) - rParticles%p_xpos_old(iPart)
			!Analog die y-Koordinaten:
			y(1)= p_DvertexCoords(2,bdy_Koords(1))
			y(2)= p_DvertexCoords(2,bdy_Koords(2)) - p_DvertexCoords(2,bdy_Koords(1))
			y(3)= rParticles%p_ypos_old(iPart)
			y(4)= rParticles%p_ypos(iPart) - rParticles%p_ypos_old(iPart)

			
			if ((x(4)*y(2)).NE.(y(4)*x(2))) then
				s=(x(4)*(y(3)-y(1))-y(4)*(x(3)-x(1)))/(x(4)*y(2)-y(4)*x(2))
				!Bestimme den Punkt, an dem das Teilchen auf den Rand trifft
				bdy_point(1)= x(1) + s*x(2)
				bdy_point(2)= y(1) + s*y(2)
				!Anteil der Bewegung am Zeitschritt
				rParticles%p_bdy_time(iPart)= (sqrt((bdy_point(1) - rParticles%p_xpos_old(iPart))**2+& 
											(bdy_point(2) - rParticles%p_ypos_old(iPart))**2))&
											/(sqrt((rParticles%p_xpos(iPart) - rParticles%p_xpos_old(iPart))**2+&
											(rParticles%p_ypos(iPart) - rParticles%p_ypos_old(iPart))**2))
			    rParticles%p_xpos(iPart) = bdy_point(1)
			    rParticles%p_ypos(iPart) = bdy_point(2)


		        !Der "Restvektor" der Geschwindigkeit, den das Teilchen nicht mehr fliegt, da es auf die Rand trifft
		        velo_rest(1)= rParticles%p_xpos(iPart) - bdy_point(1)
		        velo_rest(2)= rParticles%p_ypos(iPart) - bdy_point(2)

		        !Tangente in x-Richtung
		        bdy_tang(1)= p_DvertexCoords(1,bdy_Koords(2)) - p_DvertexCoords(1,bdy_Koords(1))
		        !Tangente in y-Richtung
		        bdy_tang(2)= p_DvertexCoords(2,bdy_Koords(2)) - p_DvertexCoords(2,bdy_Koords(1))
		        tang_norm = sqrt(bdy_tang(1)**2+bdy_tang(2)**2)
		        bdy_tang(1)= bdy_tang(1)/tang_norm
		        bdy_tang(2)= bdy_tang(2)/tang_norm
		        !Normale in x-Richtung
		        bdy_norm(1)= bdy_tang(2)
		        !Normale in y-Richtung
		        bdy_norm(2)= -bdy_tang(1)

		        !Projektion auf den Rand
		        proj_tang = velo_rest(1)*bdy_tang(1)+velo_rest(2)*bdy_tang(2)
		        proj_norm = velo_rest(1)*bdy_norm(1)+velo_rest(2)*bdy_norm(2)

			    !Berechne den Restvektor der Bewegung
			    bdy_move(1)= rParticles%tang_val*proj_tang*bdy_tang(1)+&
							 rParticles%norm_val*proj_norm*bdy_norm(1)
			    bdy_move(2)= rParticles%tang_val*proj_tang*bdy_tang(2)+&
							 rParticles%norm_val*proj_norm*bdy_norm(2)

			    !Berechne die neuen Koordinaten (ausgehend von dem Punkt, an dem das Teilchen auf den Rand trifft)
			    rParticles%p_xpos(iPart) = rParticles%p_xpos(iPart) + bdy_move(1)
			    rParticles%p_ypos(iPart) = rParticles%p_ypos(iPart) + bdy_move(2)

			    !Projektion auf den Rand
			    proj_tang = (rParticles%p_xpos_old(iPart) - rParticles%p_xpos(iPart))*bdy_tang(1)+&
						    (rParticles%p_ypos_old(iPart) - rParticles%p_ypos(iPart))*bdy_tang(2)
			    proj_norm = (rParticles%p_xpos_old(iPart) - rParticles%p_xpos(iPart))*bdy_norm(1)+&
				    		(rParticles%p_ypos_old(iPart) - rParticles%p_ypos(iPart))*bdy_norm(2)

			    rParticles%p_xvelo(iPart)= rParticles%tang_val*proj_tang*bdy_tang(1)+&
				    				       rParticles%norm_val*proj_norm*bdy_norm(1)
			    rParticles%p_yvelo(iPart)= rParticles%tang_val*proj_tang*bdy_tang(2)+&
				       			           rParticles%norm_val*proj_norm*bdy_norm(2)
			end if
			exit Boundary_Search
		end if
	end do Boundary_Search

    rParticles%p_bdy_check(iPart) = 0

    ! Find the new element for the particle
    call findnewelement(rparlist,p_rproblemLevel,rParticles,iPart)

    ! calculate barycentric coordinates
    call calculatebarycoords(p_rproblemLevel,rParticles,iPart)

    ! wrong element
    if ((abs(rParticles%p_lambda1(iPart))+abs(rParticles%p_lambda2(iPart))+&
         abs(rParticles%p_lambda3(iPart))-1) .GE. 0.00001) then
        call wrongelement(rparlist,p_rproblemLevel,rParticles,iPart)
    end if


end subroutine checkboundary

!*
!**********************************************************************************************************************

!************************************ SUBROUTINE to calculate the volumepart ******************************************
!*

subroutine calculatevolumepart(p_rproblemLevel,rParticles)

    ! particles
    type(t_Particles), intent(inout) :: rParticles

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to array of the volume for each element
    !
    ! 2D triangulation: Array with area of each element.
    ! 3D triangulation: Array with volume of each element.
    ! Handle to 
    !       p_DelementArea = array [1..NEL+1] of double.
    ! p_DelementArea [NEL+1] gives the total area/voloume of the domain.
    real(DP), dimension(:), pointer :: p_DelementVolume

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! local variables
	integer :: i, current
	real(DP) :: c_pi

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
 
    ! Get vertices at element
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Get area of each element
    call storage_getbase_double(&
         p_rtriangulation%h_DelementVolume, p_DelementVolume)

    c_pi= 3.14159265358979323846264338327950288

	do i= 1, rParticles%nPart

        ! element of the current particle
		current= rParticles%p_element(i)

        ! store the volumefraction of the particle in the gridpoints (with barycentric coordinates)
		rParticles%p_PartVol(p_IverticesAtElement(1,current))= &
		                rParticles%p_PartVol(p_IverticesAtElement(1,current)) + &
		                (rParticles%p_lambda1(i)*c_pi*0.25*rParticles%p_diam(i)**2)/&
		                p_DelementVolume(p_IverticesAtElement(1,current))
		rParticles%p_PartVol(p_IverticesAtElement(2,current))= &
		                rParticles%p_PartVol(p_IverticesAtElement(2,current)) + &
		                (rParticles%p_lambda2(i)*c_pi*0.25*rParticles%p_diam(i)**2)/&
		                p_DelementVolume(p_IverticesAtElement(2,current))
		rParticles%p_PartVol(p_IverticesAtElement(3,current))= &
		                rParticles%p_PartVol(p_IverticesAtElement(3,current)) + &
		                (rParticles%p_lambda3(i)*c_pi*0.25*rParticles%p_diam(i)**2)/&
		                p_DelementVolume(p_IverticesAtElement(3,current))

	end do


end subroutine calculatevolumepart

!*
!**********************************************************************************************************************


!********************** SUBROUTINE to calculate the velocity of the particles in the gridpoint ************************
!*

subroutine calculatevelopart(p_rproblemLevel,rParticles)

    ! particles
    type(t_Particles), intent(inout) :: rParticles

    ! Pointer to the multigrid level
    type(t_problemLevel), pointer :: p_rproblemLevel

    ! pointer to the vertices adjacent to an element
    !
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    
    ! pointer to array of the volume for each element
    !
    ! 2D triangulation: Array with area of each element.
    ! 3D triangulation: Array with volume of each element.
    ! Handle to 
    !       p_DelementArea = array [1..NEL+1] of double.
    ! p_DelementArea [NEL+1] gives the total area/voloume of the domain.
    real(DP), dimension(:), pointer :: p_DelementVolume

    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! local variables
	integer :: i, current
	real(DP) :: c_pi

    ! Set pointer to triangulation
    p_rtriangulation => p_rproblemLevel%rtriangulation
 
    ! Get vertices at element
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Get area of each element
    call storage_getbase_double(&
         p_rtriangulation%h_DelementVolume, p_DelementVolume)

	do i= 1, rParticles%nPart

        ! element of the current particle
		current= rParticles%p_element(i)

        ! store the velocity of the particle in the gridpoints (with barycentric coordinates)
		rParticles%p_PartVelo(1,p_IverticesAtElement(1,current))= &
		                rParticles%p_PartVelo(1,p_IverticesAtElement(1,current)) + &
		                (rParticles%p_xvelo(i)**2)/&
		                p_DelementVolume(p_IverticesAtElement(1,current))
		rParticles%p_PartVelo(1,p_IverticesAtElement(2,current))= &
		                rParticles%p_PartVelo(1,p_IverticesAtElement(2,current)) + &
		                (rParticles%p_xvelo(i)**2)/&
		                p_DelementVolume(p_IverticesAtElement(2,current))
		rParticles%p_PartVelo(1,p_IverticesAtElement(3,current))= &
		                rParticles%p_PartVelo(1,p_IverticesAtElement(3,current)) + &
		                (rParticles%p_xvelo(i)**2)/&
		                p_DelementVolume(p_IverticesAtElement(3,current))
		rParticles%p_PartVelo(2,p_IverticesAtElement(1,current))= &
		                rParticles%p_PartVelo(2,p_IverticesAtElement(1,current)) + &
		                (rParticles%p_yvelo(i)**2)/&
		                p_DelementVolume(p_IverticesAtElement(1,current))
		rParticles%p_PartVelo(2,p_IverticesAtElement(2,current))= &
		                rParticles%p_PartVelo(2,p_IverticesAtElement(2,current)) + &
		                (rParticles%p_yvelo(i)**2)/&
		                p_DelementVolume(p_IverticesAtElement(2,current))
		rParticles%p_PartVelo(2,p_IverticesAtElement(3,current))= &
		                rParticles%p_PartVelo(2,p_IverticesAtElement(3,current)) + &
		                (rParticles%p_yvelo(i)**2)/&
		                p_DelementVolume(p_IverticesAtElement(3,current))

	end do

end subroutine calculatevelopart

!*
!**********************************************************************************************************************


end module euler_lagrange
