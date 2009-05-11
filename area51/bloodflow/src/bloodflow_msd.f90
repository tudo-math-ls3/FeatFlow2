!##############################################################################
!# ****************************************************************************
!# <name> bloodflow_msd </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides all data structures and subroutine to solve
!# the mass-spring-damper model which is used to redistribute mesh points.
!# </purpose>
!##############################################################################

module bloodflow_msd

  use fsystem
  use storage
  use triangulation
  
  implicit none

  ! Every subroutine is declared private by default
  private

  ! Subroutine which are accessable from outside 
  ! this module are explicitly declared public
  public :: t_msd
  public :: bloodflow_createMSDSystem
  public :: bloodflow_releaseMSDSystem
  public :: bloodflow_calcForces

  !*****************************************************************************

!<types>

!<typeblock>

  ! This data structure represents a single particle

  type t_particle

    ! Mass of the particle
    real(DP) :: dmass = 1.0_DP

    ! Marker for fixed particles
    logical :: bfixed = .false.

    ! Position of the particle
    real(DP), dimension(NDIM2D) :: Dcoords = 0.0_DP

    ! Velocity of the particle
    real(DP), dimension(NDIM2D) :: Dvelocity = 0.0_DP

    ! Force acting on the particle
    real(DP), dimension(NDIM2D) :: Dforce = 0.0_DP

  end type t_particle

!</typeblock>

!<typeblock>

  ! This data structure represents a single spring

  type t_spring

    ! Elasticity constant of the spring
    real(DP) :: delasticity = 0.1_DP

    ! Damping constant of the spring
    real(DP) :: ddamping = 0.01_DP

    ! Rest length of the spring
    real(DP) :: drestLength = 0.0_DP
    
  end type t_spring

!</typeblock>

!<typeblock>

  ! This data structure represents a mass-spring-damper model

  type t_msd

    ! Number of particles
    integer :: nparticles = 0

    ! Number of springs
    integer :: nsprings = 0

    ! Data array for particles
    type(t_particle), dimension(:), pointer :: Rparticles => null()

    ! Data array for springs
    type(t_spring), dimension(:), pointer :: Rsprings => null()

  end type t_msd

!</types>

  !*****************************************************************************

contains

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_createMSDSystem(rmsd, rtriangulation)

!<description>
    
    ! This subroutine creates the mass-spring-damper system and
    ! initializes the particles and springs based on the given
    ! triangulation structure.

!</description>

!<input>

    ! Triangulation structure
    type(t_triangulation), intent(IN) :: rtriangulation

!</input>

!<output>

    ! Mass-spring-damper system
    type(t_msd), intent(OUT) :: rmsd

!</output>

!</subroutine>
    
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IneighboursAtElement, p_IverticesAtElement
    real(DP), dimension(NDIM2D) :: Dx
    integer :: i, j, k, iel, ive, ispring
    
    ! Set number of particles
    rmsd%nparticles = rtriangulation%NVT

    ! Allocate memory
    allocate(rmsd%Rparticles(rmsd%nparticles))
    
    ! Set pointer to physical vertex coordinates
    call storage_getbase_double2d(rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Loop over all particles
    do i = 1, rmsd%nparticles     

      ! Set particle position
      rmsd%Rparticles(i)%Dcoords = p_DvertexCoords(:, i)

      ! Set initial velocity to zero
      rmsd%Rparticles(i)%Dvelocity = 0.0_DP
    end do

    ! Set pointer
    call storage_getbase_int2D(rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    call storage_getbase_int2D(rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    
    ! Determine the number of edge springs
    ispring = 0   

    do iel = 1, rtriangulation%NEL
      do ive = 1, 3
        ! We only count the edge shared by elements IEL and JEL if JEL < IEL
        if (p_IneighboursAtElement(ive, iel) .gt. iel) cycle
        ispring = ispring+1
      end do
    end do
    rmsd%nsprings = ispring + 3*rtriangulation%NEL

    ! Allocate memory
    allocate(rmsd%Rsprings(rmsd%nsprings))

    ! Compute the rest length for the edge springs
    ispring = 0
    do iel = 1, rtriangulation%NEL
      do ive = 1, 3

        ! Skip edges which are addressed from the adjacent element
        if (p_IneighboursAtElement(ive, iel) .gt. iel) cycle
        ispring = ispring+1
        
        ! Get global vertex numbers
        i = p_IverticesAtElement(ive, iel)
        j = p_IverticesAtElement(mod(ive, 3), iel)

        ! Compute rest length of the spring
        Dx = p_DvertexCoords(:,i) - p_DvertexCoords(:,j)
        rmsd%Rsprings(ispring)%drestLength = sqrt(sum(Dx*Dx))
      end do
    end do

    ! Compute the rest length for the projection springs
    do iel = 1, rtriangulation%NEL
      
      ! Get global vertex numbers
      i = p_IverticesAtElement(1, iel)
      j = p_IverticesAtElement(2, iel)
      k = p_IverticesAtElement(3, iel)

      ! Compute rest length of the three projection springs
      rmsd%Rsprings(ispring+1)%drestLength = abs( (p_DvertexCoords(1,k)-p_DvertexCoords(1,j)) * &
                                                  (p_DvertexCoords(2,j)-p_DvertexCoords(2,i)) - &
                                                  (p_DvertexCoords(1,j)-p_DvertexCoords(1,i)) *&
                                                  (p_DvertexCoords(2,k)-p_DvertexCoords(2,j)) ) /&
                                            sqrt( (p_DvertexCoords(1,k)-p_DvertexCoords(1,j))**2 + &
                                                  (p_DvertexCoords(2,k)-p_DvertexCoords(2,j))**2 )

      rmsd%Rsprings(ispring+2)%drestLength = abs( (p_DvertexCoords(1,k)-p_DvertexCoords(1,i)) * &
                                                  (p_DvertexCoords(2,i)-p_DvertexCoords(2,j)) - &
                                                  (p_DvertexCoords(1,i)-p_DvertexCoords(1,j)) *&
                                                  (p_DvertexCoords(2,k)-p_DvertexCoords(2,i)) ) /&
                                            sqrt( (p_DvertexCoords(1,k)-p_DvertexCoords(1,i))**2 + &
                                                  (p_DvertexCoords(2,k)-p_DvertexCoords(2,i))**2 )

      rmsd%Rsprings(ispring+3)%drestLength = abs( (p_DvertexCoords(1,j)-p_DvertexCoords(1,i)) * &
                                                  (p_DvertexCoords(2,i)-p_DvertexCoords(2,k)) - &
                                                  (p_DvertexCoords(1,i)-p_DvertexCoords(1,k)) *&
                                                  (p_DvertexCoords(2,j)-p_DvertexCoords(2,i)) ) /&
                                            sqrt( (p_DvertexCoords(1,j)-p_DvertexCoords(1,i))**2 + &
                                                  (p_DvertexCoords(2,j)-p_DvertexCoords(2,i))**2 )
      
      ! Update spring counter
      ispring = ispring+3
      
    end do

  end subroutine bloodflow_createMSDSystem

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_releaseMSDSystem(rmsd)

!<description>

    ! This subroutine releases the mass-spring-damper system

!</description>

!<inputoutput>

    ! Mass-spring-damper system
    type(t_msd), intent(INOUT) :: rmsd

!</intputoutput>

!</subroutine>

    ! Reset number of particles
    rmsd%nparticles = 0

    ! Deallocate memory
    if (associated(rmsd%Rparticles)) then
      deallocate(rmsd%Rparticles)
      nullify(rmsd%Rparticles)
    end if

    ! Deallocate memory
    if (associated(rmsd%Rsprings)) then
      deallocate(rmsd%Rsprings)
      nullify(rmsd%Rsprings)
    end if

  end subroutine bloodflow_releaseMSDSystem

  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_calcForces(rmsd, rtriangulation)

!<description>

    ! This subroutine calculates the forces applied to each particle

!</description>

!<input>

    ! Triangulation structure
    type(t_triangulation), intent(IN) :: rtriangulation

!</input>

!<inputoutput>

    ! Mass-spring-damper system
    type(t_msd), intent(INOUT) :: rmsd

!</intputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(NDIM2D) :: Dforce,Dx
    real(DP) :: dSpringLength
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    integer :: i,j,k,iel,ive,ispring
    

    ! Set all external forces to zero
    do i = 1, rmsd%nparticles
      rmsd%Rparticles(i)%Dforce = 0.0_DP
    end do

    ! Set pointer to triangulation data
    call storage_getbase_int2d(rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2d(rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)

    
    ! Compute the accumulated forces due to edge springs
    ispring = 0
    do iel = 1, rtriangulation%NEL

      ! Loop over all edges of the element
      do ive = 1, 3

        ! Skip edges which are addressed from the adjacent element
        if (p_IneighboursAtElement(ive, iel) .gt. iel) cycle
        ispring = ispring+1
        
        ! Get global vertex numbers
        i = p_IverticesAtElement(ive, iel)
        j = p_IverticesAtElement(mod(ive,3)+1, iel)
        
        ! Compute length of the spring
        Dx = rmsd%Rparticles(i)%Dcoords - rmsd%Rparticles(j)%Dcoords
        dSpringLength = sqrt(sum(Dx*Dx))
        
        Dforce = rmsd%Rsprings(ispring)%delasticity *&
                 (dSpringLength-rmsd%Rsprings(ispring)%drestLength)
        
        Dforce = Dforce +&
                 rmsd%Rsprings(ispring)%ddamping *&
                 (rmsd%Rparticles(i)%Dvelocity - rmsd%Rparticles(j)%Dvelocity) *&
                 Dx / dSpringLength

        Dforce = -Dforce * Dx / dSpringLength

        ! Apply force to node i (if required)
        if (.not.rmsd%Rparticles(i)%bfixed)&
            rmsd%Rparticles(i)%Dforce = rmsd%Rparticles(i)%Dforce + Dforce

        ! Apply force to node j (if required)
        if (.not.rmsd%Rparticles(j)%bfixed)&
            rmsd%Rparticles(j)%Dforce = rmsd%Rparticles(j)%Dforce - Dforce
      end do
    end do

    ! Compute the accumulated forces due to projection springs
    do iel = 1, rtriangulation%NEL
      
      ! Get global vertex numbers
      i = p_IverticesAtElement(1, iel)
      j = p_IverticesAtElement(2, iel)
      k = p_IverticesAtElement(3, iel)

      ! Compute length of the spring

      ! Update spring counter
      ispring = ispring+3

    end do

  end subroutine bloodflow_calcForces
end module bloodflow_msd
