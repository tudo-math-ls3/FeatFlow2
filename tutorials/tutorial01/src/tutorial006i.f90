!##############################################################################
!# Tutorial 006i: Create a mass and Laplace matrix
!##############################################################################

module tutorial006i

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
  
  implicit none
  private
  
  public :: start_tutorial006i

contains

  ! ***************************************************************************

  subroutine start_tutorial006i

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialdiscretisation) :: rdiscretisation
    
    type(t_matrixScalar) :: rmatrixTemplate
    type(t_matrixScalar) :: rmatrixMass, rmatrixLaplace
    
    logical :: bhasstruc1, bhasstruc2, bhasstruc3
    logical :: bhascont1, bhascont2, bhascont3
    logical :: bsharedstruc1, bsharedstruc2, bsharedstruc3
    logical :: bsharedcont

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 006i")
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
    
    ! NOTE: "lsyssc_duplicateMatrix" allows to "derive" a matrix from another one
    ! in various ways. In the following, we use the routine to derive a mass
    ! matrix rmatrixMass from the template matrix rmatrixTemplate. The mass
    ! matrix "shares" the (CSR-) structure of the template matrix (->LSYSSC_DUP_SHARE),
    ! so no additional memory is allocated. For the "content" (the entries), 
    ! a new and empty memory block is allocated (->LSYSSC_DUP_EMPTY) which we later
    ! fill with data.

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

    ! =================================
    ! Prove that the structure is shared.
    ! =================================
    
    ! Check which matrices have data and what is shared
    bhasstruc1    = lsyssc_hasMatrixStructure(rmatrixTemplate)
    bhasstruc2    = lsyssc_hasMatrixStructure(rmatrixMass)
    bhasstruc3    = lsyssc_hasMatrixStructure(rmatrixLaplace)

    bhascont1     = lsyssc_hasMatrixContent(rmatrixTemplate)
    bhascont2     = lsyssc_hasMatrixContent(rmatrixMass)
    bhascont3     = lsyssc_hasMatrixContent(rmatrixLaplace)

    bsharedstruc1 = lsyssc_isMatrixStructureShared(rmatrixTemplate,rmatrixMass)
    bsharedstruc2 = lsyssc_isMatrixStructureShared(rmatrixTemplate,rmatrixLaplace)
    bsharedstruc3 = lsyssc_isMatrixStructureShared(rmatrixMass,rmatrixLaplace)

    bsharedcont   = lsyssc_isMatrixContentShared(rmatrixMass,rmatrixLaplace)
    
    ! Print the results
    call output_line ("rmatrixTemplate has structure    = " // &
        trim(sys_sl(bhasstruc1)))

    call output_line ("rmatrixMass has structure        = " // &
        trim(sys_sl(bhasstruc2)))

    call output_line ("rmatrixLaplace has structure     = " // &
        trim(sys_sl(bhasstruc3)))

    call output_line ("rmatrixTemplate has content      = " // &
        trim(sys_sl(bhascont1)))

    call output_line ("rmatrixMass has content          = " // &
        trim(sys_sl(bhascont2)))

    call output_line ("rmatrixLaplace has content       = " // &
        trim(sys_sl(bhascont3)))

    call output_line ("Structure shared(rmatrixTemplate,rmatrixMass)    = " // &
        trim(sys_sl(bsharedstruc1)))

    call output_line ("Structure shared(rmatrixTemplate,rmatrixLaplace) = " // &
        trim(sys_sl(bsharedstruc2)))

    call output_line ("Structure shared(rmatrixMass,rmatrixLaplace)     = " // &
        trim(sys_sl(bsharedstruc3)))

    call output_line ("Content shared(rmatrixMass,rmatrixLaplace)       = " // &
        trim(sys_sl(bsharedcont)))

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
