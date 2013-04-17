!##############################################################################
!# Tutorial 006k: Create a RHS vector
!##############################################################################

module tutorial006k

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use element
  use derivatives
  use spatialdiscretisation
  use linearsystemscalar

  use scalarpde
  use linearformevaluation
  
  use vectorio
  use ucd
  use collection

  implicit none
  private
  
  public :: start_tutorial006k

contains

  ! ***************************************************************************

!<subroutine>

  subroutine fcoeff_RHS (rdiscretisation, rform, &
                nelements, npointsPerElement, Dpoints, &
                IdofsTest, rdomainIntSubset, &
                Dcoefficients, rcollection)

  use basicgeometry
  use collection
  use domainintegration
  use fsystem
  use scalarpde
  use spatialdiscretisation
  use triangulation

!<description>
  ! This subroutine is called during the vector assembly. It returns
  ! the values of a RHS function f in a set of (cubature) points
  ! on a set of elements.
!</description>

!<input>
  ! Underlying discretisation structure.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! The linear form which is currently to be evaluated:
  type(t_linearForm), intent(in) :: rform

  ! Number of elements, where the coefficients must be computed.
  integer, intent(in) :: nelements

  ! Number of points per element, where the coefficients must be computed
  integer, intent(in) :: npointsPerElement

  ! Array of all points (x/y coords) on all elements where the values are needed.
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! Array with degrees of freedom of the test function.
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! A t_domainIntSubset structure specifying more detailed assembly information.
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
  ! Optional: A collection structure to provide additional information.
  type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
  ! A list of all coefficients in front of all terms in the linear form -
  ! for all given points on all given elements.
  !   DIMENSION(itermCount,npointsPerElement,nelements)
  ! with itermCount the number of terms in the linear form.
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    integer :: ielement, ipoint
    real(DP) :: dx, dy, dmult
    
    ! Get the multiplier from the collection -- will be 4.0 by default (see below).
    dmult = rcollection%DquickAccess(1)
    
    do ielement = 1,nelements
      do ipoint = 1,npointsPerElement
      
        ! x/y coordinate
        dx = Dpoints(1,ipoint,ielement)
        dy = Dpoints(2,ipoint,ielement)
        
        ! Write f(x,y)=d*16*x*(1-x)*y*(1-y) into Dcoefficients(1,:,:)
        Dcoefficients(1,ipoint,ielement) = &
            dmult * 16.0_DP * dx * (1.0_DP - dx) * dy * (1.0_DP - dy)
      
      end do
    end do
    
  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial006k

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    
    type(t_linearForm) :: rlinform
    type(t_collection) :: rcollection
    type(t_vectorScalar) :: rrhs
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006k")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 5x5-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 4, 4)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Discretise with Q1.
    ! =================================

    call spdiscr_initDiscr_simple (rdiscretisation,EL_Q1_2D,rtriangulation)

    ! =================================
    ! Create a RHS vector.
    ! =================================
    
    ! Create a vector.
    call lsyssc_createVector (rdiscretisation,rrhs)
    
    ! Prepare a linear form structure for one term in the RHS: (f, phi).
    ! The test function is just phi without derivative.
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC2D
    
    ! Build the RHS using fcoeff_RHS above. Pass a multiplier via a collection structure.
    ! rcollection%DquickAccess/IquickAccess/... can arbitrarily be used to pass values.
    rcollection%DquickAccess(1) = 4.0_DP
    
    call linf_buildVectorScalar (&
        rdiscretisation,rlinform,.true.,rrhs,fcoeff_RHS,rcollection)

    ! =================================
    ! Output to files. 
    ! =================================
    call output_line ("Writing RHS vector to text file...")

    ! Write the vector to a text file.
    call vecio_writeVectorHR (rrhs, "vector", .true., 0, &
        "post/tutorial006k_rhs.txt", "(E20.10)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
