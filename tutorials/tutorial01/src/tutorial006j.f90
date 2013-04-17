!##############################################################################
!# Tutorial 006j: Create a system matrix from a mass and a Laplace matrix
!##############################################################################

module tutorial006j

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
  use bilinearformevaluation
  use stdoperators
  
  use matrixio
  use ucd

  implicit none
  private
  
  public :: start_tutorial006j

contains

  ! ***************************************************************************

  subroutine start_tutorial006j

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    
    type(t_matrixScalar) :: rmatrixTemplate
    type(t_matrixScalar) :: rmatrixMass, rmatrixLaplace
    type(t_matrixScalar) :: rmatrix
    
    logical :: bhasstruc1, bhasstruc2, bhasstruc3
    logical :: bhascont1, bhascont2, bhascont3
    logical :: bsharedstruc1, bsharedstruc2, bsharedstruc3
    logical :: bsharedcont


    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006j")
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
    ! Create a CSR-matrix for the above
    ! discretisation.
    ! =================================

    ! Create the matrix structure
    call bilf_createMatrixStructure (rdiscretisation,LSYSSC_MATRIX9,rmatrixTemplate)

    ! =================================
    ! Derive a mass and a Laplace
    ! matrix which share the same
    ! structure but have their own
    ! content.
    ! =================================

    ! ------------------------------------
    ! ----- Mass matrix ------------------

    ! Share structure, allocate memory for new content    
    call lsyssc_duplicateMatrix (rmatrixTemplate,rmatrixMass,&
        LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

    ! Clear the matrix
    call lsyssc_clearMatrix (rmatrixMass)

    ! Create the matrix
    call stdop_assembleSimpleMatrix (rmatrixMass,DER_FUNC2D,DER_FUNC2D)

    ! ------------------------------------
    ! ----- Laplace matrix ---------------
    call lsyssc_duplicateMatrix (rmatrixTemplate,rmatrixLaplace,&
        LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
        
    call lsyssc_clearMatrix (rmatrixLaplace)
    
    call stdop_assembleLaplaceMatrix (rmatrixLaplace)

    ! ------------------------------------
    ! ----- System matrix ----------------

    ! Share the structure, provide empty space for the content
    call lsyssc_duplicateMatrix (rmatrixTemplate,rmatrix,&
        LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

    ! =================================
    ! Sum up Laplace and mass matrix.
    ! =================================
    
    ! We create the system as
    !
    !    rmatrix = 1/100 Mass + Laplace

    ! Clear the destination
    call lsyssc_clearMatrix (rmatrix)
    
    ! Incorporate (1/100) Mass to the system matrix:
    !    rmatrix = 1/100 Mass + rmatrix(=0)
    call lsyssc_matrixLinearComb (rmatrixMass,rmatrix,1.0_DP/100.0_DP,1.0_DP,&
        .false.,.false.,.true.,.true.)
    
    ! Incorporate Laplace to the system matrix
    !    rmatrix = Mass + rmatrix
    call lsyssc_matrixLinearComb (rmatrixMass,rmatrix,1.0_DP,1.0_DP,&
        .false.,.false.,.true.,.true.)

    ! =================================
    ! Write the matrix to a file.
    ! =================================
    call output_line ("Writing matrix to a text file...")

    ! Write the matrix to a text file, omit nonexisting entries in the matrix.
    call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
        "post/tutorial006j_matrix.txt", "(E15.5)")

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the matrices
    call lsyssc_releaseMatrix (rmatrixMass)
    call lsyssc_releaseMatrix (rmatrixLaplace)
    
    ! Release template matrix
    call lsyssc_releaseMatrix (rmatrixTemplate)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
