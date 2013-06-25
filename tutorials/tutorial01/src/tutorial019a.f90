!##############################################################################
!# Tutorial 019a: Projection of an analytically given function into FEM
!##############################################################################

module tutorial019a

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use bilinearformevaluation
  
  use domainintegration
  use collection
  
  use matrixio
  use vectorio
  use ucd
  
  use analyticprojection

  implicit none
  private
  
  public :: start_tutorial019a

contains

  ! ***************************************************************************

!<subroutine>

  subroutine ffunction (cderivative, rdiscretisation, &
      nelements, npointsPerElement, Dpoints, IdofsTest, rdomainIntSubset, &
      Dvalues, rcollection)

!<description>
  ! This subroutine returns values of an analytically given function
  ! at specific points in a domain.
!</description>

!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in) :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! Number of elements, where the coefficients must be computed.
  integer, intent(in) :: nelements

  ! Number of points per element, where the coefficients must be computed
  integer, intent(in) :: npointsPerElement

  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in) :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional :: rcollection

!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: ipt, iel
    real(DP) :: dx,dy

    ! In every point on every element, return the value of our function.
    do iel = 1,nelements
      do ipt = 1,npointsPerElement
      
        ! x/y coordinates
        dx = Dpoints(1,ipt,iel)
        dy = Dpoints(2,ipt,iel)
        
        ! Value of the function. Here: u=||(x,y)||
        Dvalues(ipt,iel) = sqrt ( dx**2 + dy**2 )
      
      end do
    end do

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial019a

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    type(t_vectorScalar) :: rx
    type(t_ucdExport) :: rexport

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 019a")
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
    !
    ! Create a structure rdiscretisation
    ! which describes the discretisation.
    ! =================================

    call spdiscr_initDiscr_simple (rdiscretisation,EL_Q1_2D,rtriangulation)

    ! =================================
    ! Create a scalar vector.
    ! =================================

    ! Create a vector.
    call lsyssc_createVector (rdiscretisation,rx)
    
    ! =================================
    ! Projection of the above function into our FEM space
    ! =================================
    call anprj_discrDirect (rx, ffunction)

    ! =================================
    ! Write a VTK file with the mesh
    ! and this solution vector.
    ! =================================
    call output_line ("Writing file 'post/tutorial019a.vtk'.")

    ! Open
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "post/tutorial019a.vtk")
                       
    ! Pass the vector as solution.
    call ucd_addVectorByVertex (rexport, "x", UCD_VAR_STANDARD, rx)
          
    ! Write / close             
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the vector
    call lsyssc_releaseVector (rx)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
