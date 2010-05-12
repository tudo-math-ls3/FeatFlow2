!##############################################################################
!# ****************************************************************************
!# <name> spatialoperators </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains routines for assembling and applying operators in space.
!# </purpose>
!##############################################################################

module spatialoperators

  use fsystem
  use storage
  use genoutput
  use linearalgebra
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use collection
  use vectorio
  use derivatives
  
  use physics
  
  use stdoperators
  use bilinearformevaluation
  
  implicit none

  ! Defines a couple of template matrices which are precomputed at the beginning
  ! of the program and stay unchanged as long as the associated discretisation/
  ! triangulation does not change.
  type t_matvecTemplates
  
    ! Associated discretisation structure in space
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr
  
    ! General FE template matrix.
    type(t_matrixScalar) :: rmatrixTemplateA11
  
    ! Mass matrix
    type(t_matrixScalar) :: rmatrixMassA11
    
    ! Laplace matrix
    type(t_matrixScalar) :: rmatrixLaplaceA11
  
  end type

contains

  ! ***************************************************************************

  subroutine spop_calcMatrices (rspaceDiscr, rmatvecTempl)

    ! Calculates template matrices in rmatvecTempl
    
    ! Underlying spatial discretisation structure
    type(t_blockDiscretisation), intent(in), target :: rspaceDiscr
    
    ! Template structure to be created
    type(t_matVecTemplates), intent(inout) :: rmatvecTempl
    
    ! Remember the discretisation structure
    rmatvecTempl%p_rspaceDiscr => rspaceDiscr
    
    ! Generate the underlying FE matrix structure.
    call bilf_createMatrixStructure (rspaceDiscr%RspatialDiscr(1),LSYSSC_MATRIX9,&
        rmatvecTempl%rmatrixTemplateA11)
    
    ! Calculate the Laplace matrix
    call lsyssc_duplicateMatrix(&
        rmatvecTempl%rmatrixTemplateA11,rmatvecTempl%rmatrixLaplaceA11,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    call stdop_assembleLaplaceMatrix (rmatvecTempl%rmatrixLaplaceA11,.true.,1.0_DP)    
    
    ! Calculate the mass matrix
    call lsyssc_duplicateMatrix(&
        rmatvecTempl%rmatrixTemplateA11,rmatvecTempl%rmatrixMassA11,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    call stdop_assembleSimpleMatrix (rmatvecTempl%rmatrixMassA11,DER_FUNC,DER_FUNC,1.0_DP,&
      .true.)
  
  end subroutine

  ! ***************************************************************************

  subroutine spop_releaseMatrices (rmatvecTempl)

    ! Cleans up matrices in rmatvecTempl
    
    ! Template structure to be created
    type(t_matVecTemplates), intent(inout) :: rmatvecTempl
    
    ! Release...
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixLaplaceA11)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixMassA11)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixTemplateA11)
  
  end subroutine

  ! ***************************************************************************

  subroutine spop_allocateMatrix (rmatvecTempl,rphysics,rmatrix)

    ! Allocates memory for a spatial matrix.
    
    ! Template structure 
    type(t_matVecTemplates), intent(in) :: rmatvecTempl
    
    ! Physics of the problem. Defines the matrix structure.
    type(t_physics), intent(in) :: rphysics
    
    ! Output matrix
    type(t_matrixBlock), intent(out) :: rmatrix
  
    ! Set up the corresponding block matrix.
    call lsysbl_createMatBlockByDiscr (rmatvecTempl%p_rspaceDiscr,rmatrix)

    select case (rphysics%cequation)
    case (0)
      ! Heat equation.
      !
      ! A FE-matrices at block (1..2,1..2) 
      call lsyssc_duplicateMatrix (rmatvecTempl%rmatrixTemplateA11,&
          rmatrix%RmatrixBlock(1,1), LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatvecTempl%rmatrixTemplateA11,&
          rmatrix%RmatrixBlock(2,1), LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatvecTempl%rmatrixTemplateA11,&
          rmatrix%RmatrixBlock(1,2), LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (rmatvecTempl%rmatrixTemplateA11,&
          rmatrix%RmatrixBlock(2,2), LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    end select
  
  end subroutine
  
end module
