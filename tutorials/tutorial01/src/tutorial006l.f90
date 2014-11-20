!##############################################################################
!# Tutorial 006l: L2 projection with simple defect correction loop.
!##############################################################################

module tutorial006l

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use linearalgebra
  use cubature
  use element
  use derivatives
  use spatialdiscretisation
  use linearsystemscalar
  use stdoperators

  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  
  use ucd
  use collection

  implicit none
  private
  
  public :: start_tutorial006l

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

  subroutine start_tutorial006l

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    type(t_ucdExport) :: rexport
    
    type(t_linearForm) :: rlinform
    type(t_collection) :: rcollection
    
    type(t_vectorScalar) :: rrhs, rdefect, rcorrection, rx
    type(t_matrixScalar) :: rmatrix
    type(t_scalarCubatureInfo), target :: rcubatureInfo_G3
    character(LEN=SYS_STRLEN) :: spostdir
    
    real(DP) :: dnorm
    integer :: ite
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006l")
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
    ! Create a mass matrix.
    ! =================================
    
    ! Create the matrix structure and content
    call bilf_createMatrixStructure (rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    call lsyssc_allocEmptyMatrix (rmatrix)

    ! Clear the matrix
    call lsyssc_clearMatrix (rmatrix)

    ! Create the matrix
    call stdop_assembleSimpleMatrix (rmatrix,DER_FUNC2D,DER_FUNC2D)

    ! =================================
    ! Create a RHS vector.
    ! =================================
    
    ! Use a 3x3-Gauss Formula for the RHS
    call spdiscr_createDefCubStructure (rdiscretisation,rcubatureInfo_G3,CUB_GEN_AUTO_G3)

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
        rlinform,.true.,rrhs,rcubatureInfo_G3,fcoeff_RHS,rcollection)
        
    ! Release the cubature formula
    call spdiscr_releaseCubStructure (rcubatureInfo_G3)

    ! =================================
    ! Create a solution, a defect vector and a correction vector.
    ! =================================
    
    call lsyssc_createVector (rdiscretisation,rx)
    call lsyssc_createVector (rdiscretisation,rdefect)
    call lsyssc_createVector (rdiscretisation,rcorrection)
    
    ! =================================
    ! Perform a defect correction loop:
    !
    !   x_n+1 = x_n + omega D^-1 (b - M x_n)
    !
    ! with D=diagonal of M.
    ! Stop if the residual drops below 1E-12.
    ! =================================

    ! Clear the solution.
    call lsyssc_clearVector (rx)

    ! Perform the iteration
    ite = 0
    do
      ite = ite + 1
      
      ! Create the defect:   d := b - M x_n
      call lsyssc_copyVector (rrhs,rdefect)
      call lsyssc_matVec (rmatrix, rx, rdefect, -1.0_DP, 1.0_DP)
      
      ! Check the norm ||d||_l2
      dnorm = lsyssc_vectorNorm (rdefect,LINALG_NORML2)
      
      call output_line ("Iteration " // trim(sys_siL(ite,10)) // &
                        ", ||RES|| = " // sys_sdEL(dnorm,4) )
      
      if (dnorm .le. 1E-12_DP) exit
      
      ! Multiply with "D^-1":     c := D^-1 (b - M x_n)
      call lsyssc_invertedDiagMatVec (rmatrix,rdefect,1.0_DP,rcorrection)
      
      ! Sum up:   x_n+1 = x_n + 0.7 c
      call lsyssc_vectorLinearComb (rcorrection, rx, 0.7_DP, 1.0_DP)
      
    end do

    ! =================================
    ! Write a VTK file with the mesh
    ! and this solution vector.
    ! =================================
    if (sys_getenv_string("POSTDIR",spostdir)) then
      call output_line ("Writing file "//trim(spostdir)//"'/tutorial006h.vtk'.")

      ! Open
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         trim(spostdir)//"/tutorial006l.vtk")
    else
      call output_line ("Writing file './post/tutorial006h.vtk'.")

      ! Open
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                         "./post/tutorial006l.vtk")
    end if
                       
    ! Pass rx as solution.
    call ucd_addVectorByVertex (rexport, "x", UCD_VAR_STANDARD, rx)
          
    ! Write / close             
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release vectors
    call lsyssc_releaseVector (rcorrection)
    call lsyssc_releaseVector (rdefect)
    call lsyssc_releaseVector (rrhs)
    call lsyssc_releaseVector (rx)
    
    ! Release matrix
    call lsyssc_releaseMatrix (rmatrix)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
