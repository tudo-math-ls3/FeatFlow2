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
  
    ! Template for B-matrices
    type(t_matrixScalar) :: rmatrixTemplateB

    ! Template for D-matrices
    type(t_matrixScalar) :: rmatrixTemplateD

    ! Template for C-matrices
    type(t_matrixScalar) :: rmatrixTemplateC
  
    ! Mass matrix
    type(t_matrixScalar) :: rmatrixMassA11
    
    ! Laplace matrix
    type(t_matrixScalar) :: rmatrixLaplaceA11

    ! B-matrices
    type(t_matrixScalar) :: rmatrixB1
    type(t_matrixScalar) :: rmatrixB2

    ! D-matrices
    type(t_matrixScalar) :: rmatrixD1
    type(t_matrixScalar) :: rmatrixD2
  end type

contains

  ! ***************************************************************************

  subroutine spop_calcMatrices (rspaceDiscr, rphysics, rmatvecTempl)

    ! Calculates template matrices in rmatvecTempl
    
    ! Underlying spatial discretisation structure
    type(t_blockDiscretisation), intent(in), target :: rspaceDiscr

    ! Physics of the problem.
    type(t_physics), intent(in) :: rphysics
    
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
  
    select case (rphysics%cequation)
    case (0,2)
      ! Heat equation. Nothing to do.
      
    case (1)
      ! Stokes equation. Generate additional B and D matrices

      call bilf_createMatrixStructure (rspaceDiscr%RspatialDiscr(3),LSYSSC_MATRIX9,&
          rmatvecTempl%rmatrixTemplateB,rspaceDiscr%RspatialDiscr(1))
          
      call lsyssc_transposeMatrix (rmatvecTempl%rmatrixTemplateB,rmatvecTempl%rmatrixTemplateD,&
          LSYSSC_TR_STRUCTURE)

      call lsyssc_duplicateMatrix(&
          rmatvecTempl%rmatrixTemplateB,rmatvecTempl%rmatrixB1,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      call lsyssc_duplicateMatrix(&
          rmatvecTempl%rmatrixTemplateB,rmatvecTempl%rmatrixB2,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix(&
          rmatvecTempl%rmatrixTemplateD,rmatvecTempl%rmatrixD1,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      call lsyssc_duplicateMatrix(&
          rmatvecTempl%rmatrixTemplateD,rmatvecTempl%rmatrixD2,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call stdop_assembleSimpleMatrix (rmatvecTempl%rmatrixB1,DER_FUNC,DER_DERIV_X,-1.0_DP,&
        .true.)

      call stdop_assembleSimpleMatrix (rmatvecTempl%rmatrixB2,DER_FUNC,DER_DERIV_Y,-1.0_DP,&
        .true.)

      call lsyssc_transposeMatrix (rmatvecTempl%rmatrixB1,rmatvecTempl%rmatrixD1,&
          LSYSSC_TR_CONTENT)

      call lsyssc_transposeMatrix (rmatvecTempl%rmatrixB2,rmatvecTempl%rmatrixD2,&
          LSYSSC_TR_CONTENT)
      
      ! C-matrix -- diagonal in the pressure.
      call bilf_createMatrixStructure (rspaceDiscr%RspatialDiscr(3),LSYSSC_MATRIX9,&
          rmatvecTempl%rmatrixTemplateC)
          
    case default
    
      call output_line ("Equation not supported.")
      call sys_halt()

    end select
  
  end subroutine

  ! ***************************************************************************

  subroutine spop_releaseMatrices (rmatvecTempl)

    ! Cleans up matrices in rmatvecTempl
    
    ! Template structure to be created
    type(t_matVecTemplates), intent(inout) :: rmatvecTempl
    
    ! Release...
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixB1)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixB2)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixD1)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixD2)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixTemplateB)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixTemplateC)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixTemplateD)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixLaplaceA11)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixMassA11)
    call lsyssc_releaseMatrix (rmatvecTempl%rmatrixTemplateA11)
  
  end subroutine

end module
