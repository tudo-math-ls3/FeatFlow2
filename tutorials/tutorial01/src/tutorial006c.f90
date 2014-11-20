!##############################################################################
!# Tutorial 006c: Discretise with Q1, create a finite element matrix
!##############################################################################

module tutorial006c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use bilinearformevaluation
  
  use matrixio

  implicit none
  private
  
  public :: start_tutorial006c

contains

  ! ***************************************************************************

  subroutine start_tutorial006c

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    type(t_matrixScalar) :: rmatrix
    character(LEN=SYS_STRLEN) :: spostdir

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006c")
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
    
    ! NOTE: Based on the above discretisation, we create a matrix
    ! in CSR format (matrix format LSYSSC_MATRIX9). We only create the
    ! "structure" at first, i.e., the index arrays that describe the
    ! positions of the nonzero entries.
    ! Later on, we allocate memory for the "content", i.e., the actual
    ! entries of the matrix.

    ! Create the matrix structure
    call bilf_createMatrixStructure (rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    
    ! Allocate memory for the content.
    call lsyssc_allocEmptyMatrix (rmatrix)
    
    ! =================================
    ! Fill the matrix with 1.0.
    ! =================================
    
    call lsyssc_clearMatrix (rmatrix,1.0_DP)

    ! =================================
    ! Output of the matrix structure
    ! =================================
    call output_line ("Writing matrix to text files...")
    
    if (sys_getenv_string("POSTDIR",spostdir)) then
      ! Write the matrix to a text file, omit nonexisting entries in the matrix.
      call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
          trim(spostdir)//"/tutorial006c_matrix.txt", "(E11.2)")
      
      ! Write the matrix to a MATLAB file.
      call matio_spyMatrix(&
          trim(spostdir)//"/tutorial006c_matrix.m","matrix",rmatrix,.true.)
      
      ! Write the matrix to a MAPLE file
      call matio_writeMatrixMaple (rmatrix, "matrix", 0, &
          trim(spostdir)//"/tutorial006c_matrix.maple", "(E11.2)")
    else
      ! Write the matrix to a text file, omit nonexisting entries in the matrix.
      call matio_writeMatrixHR (rmatrix, "matrix", .true., 0, &
          "./post/tutorial006c_matrix.txt", "(E11.2)")
      
      ! Write the matrix to a MATLAB file.
      call matio_spyMatrix(&
          "./post/tutorial006c_matrix.m","matrix",rmatrix,.true.)
      
      ! Write the matrix to a MAPLE file
      call matio_writeMatrixMaple (rmatrix, "matrix", 0, &
          "./post/tutorial006c_matrix.maple", "(E11.2)")
    end if

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the matrix
    call lsyssc_releaseMatrix (rmatrix)
    
    ! Release the Q1-discretisation
    call spdiscr_releaseDiscr (rdiscretisation)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
