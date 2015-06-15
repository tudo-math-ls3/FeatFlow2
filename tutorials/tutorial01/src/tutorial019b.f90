!##############################################################################
!# Tutorial 019b: Projection FEM spaces
!##############################################################################
!# Computes a Q0 FE function representing the cell size.
!# Projects the Q0 solution to Q1 via a projection, so one receives
!# the mean cell size in each grid point.
!##############################################################################

module tutorial019b

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput

  use boundary
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
  use spdiscprojection

  implicit none
  private

  public :: start_tutorial019b

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
    real(DP) :: dvol
    real(DP), dimension(:), pointer :: p_DcellSize
    integer, dimension(:), pointer :: p_Ielements

    ! From the triangulation, get an array with the cell size
    call storage_getbase_double (&
        rdiscretisation%p_rtriangulation%h_DelementVolume,p_DcellSize)

    ! Get the list of real element IDs which are behind our
    ! local element numbers 1..nelements
    p_Ielements => rdomainIntSubset%p_Ielements

    ! In every point on every element, return the value of our function.
    do iel = 1,nelements

      ! Volume of the element
      dvol = p_DcellSize (p_Ielements(iel))

      do ipt = 1,npointsPerElement

        ! Value of the function is the element volume
        Dvalues(ipt,iel) = dvol

      end do
    end do

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial019b

    ! Declare some variables.
    type(t_boundary) :: rboundary
    type(t_triangulation) :: rtriangulation
    character(LEN=SYS_STRLEN) :: spredir
    type(t_spatialDiscretisation) :: rdiscrQ0, rdiscrQ1
    type(t_vectorScalar) :: rxQ0, rxQ1
    type(t_ucdExport) :: rexport

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 019b")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Read the underlying domain
    ! and the mesh
    ! =================================

    if (sys_getenv_string("PREDIR",spredir)) then
      call boundary_read_prm(rboundary, trim(spredir)//"/bench1.prm")
    else
      call boundary_read_prm(rboundary, "pre/bench1.prm")
    end if
    ! The mesh must always be in "standard" format to work with it.
    ! First read, then convert to standard, based on rboundary.
    if (sys_getenv_string("PREDIR",spredir)) then
      call tria_readTriFile2D (rtriangulation, trim(spredir)//"/bench1.tri", rboundary)
    else
      call tria_readTriFile2D (rtriangulation, "pre/bench1.tri", rboundary)
    end if
    call tria_initStandardMeshFromRaw (rtriangulation, rboundary)

    ! =================================
    ! Discretise with Q0 and Q1.
    !
    ! Create a structure rdiscretisation
    ! which describes the discretisation.
    ! =================================

    call spdiscr_initDiscr_simple (rdiscrQ0,EL_Q0_2D,rtriangulation)
    call spdiscr_initDiscr_simple (rdiscrQ1,EL_Q1_2D,rtriangulation)

    ! =================================
    ! Create two vectors based on Q0 and Q1
    ! =================================

    call lsyssc_createVector (rdiscrQ0,rxQ0)
    call lsyssc_createVector (rdiscrQ1,rxQ1)

    ! =================================
    ! Projection of the cell size to Q0
    ! =================================

    call anprj_discrDirect (rxQ0, ffunction)

    ! =================================
    ! Projection Q0 to Q1
    ! =================================

    call lsyssc_clearVector (rxQ1)

    call spdp_projectSolutionScalar (rxQ0,rxQ1)

    ! =================================
    ! Write a VTK file with the mesh
    ! and this solution vector.
    ! =================================
    call output_line ("Writing file 'post/tutorial019b.vtk'.")

    ! Open
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       "post/tutorial019b.vtk")

    ! Pass the vectors as solutions.
    call ucd_addVectorByElement (rexport, "x_q0", UCD_VAR_STANDARD, rxQ0)
    call ucd_addVectorByVertex (rexport, "x_q1", UCD_VAR_STANDARD, rxQ1)

    ! Write / close
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================

    ! Release the vectors
    call lsyssc_releaseVector (rxQ1)
    call lsyssc_releaseVector (rxQ0)

    ! Release discretisation structures
    call spdiscr_releaseDiscr (rdiscrQ1)
    call spdiscr_releaseDiscr (rdiscrQ0)

    ! Release the triangulation and the boundary
    call tria_done (rtriangulation)
    call boundary_release (rboundary)

  end subroutine

end module
