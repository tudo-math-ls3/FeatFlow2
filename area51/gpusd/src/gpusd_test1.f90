module gpusd_test1

use fsystem
use genoutput
use storage
use boundary
use triangulation
use cubature
use element
use spatialdiscretisation
use linearsystemscalar
use stdoperators
use matrixio
use statistics

use sd_core

implicit none

contains

subroutine gpusd_1()
  
  type(t_boundary) :: rbnd
  type(t_triangulation) :: rtria
  type(t_spatialDiscretisation) :: rdiscr
  type(t_vectorScalar) :: ru1, ru2
  type(t_matrixScalar) :: rmatrix
  real(DP), dimension(:), pointer :: p_Du1, p_Du2, p_DA
  integer, dimension(:), pointer :: p_Kld, p_Kcol
  real(DP), dimension(:,:), pointer :: p_Dvtx
  integer, dimension(:,:), pointer :: p_Iverts, p_Iedges
  integer :: NLMAX, info
  real(DP) :: dnu, dalpha, dbeta, dgamma, ddelta
  type(t_timer) :: rtimer
  
    ! set up the parameters
    NLMAX = 12
    dnu = 1.0_DP
    
    dalpha = 0.0_DP     ! mass matrix
    dbeta  = 0.0_DP     ! laplace matrix
    dgamma = 1.0_DP     ! convection
    ddelta = 1.0_DP     ! stabilisation

    print *, 'FEAT2 Initialisation'
    
    ! create the triangulation
    call boundary_read_prm(rbnd, './pre/QUAD.prm')
    call tria_readTriFile2D (rtria, './pre/QUAD.tri', rbnd)
    call tria_quickRefine2LevelOrdering(NLMAX-1, rtria, rbnd)
    call tria_initStandardMeshFromRaw (rtria, rbnd)
    
    ! Fetch the necessary arrays from the triangulation
    call storage_getbase_double2D(rtria%h_DvertexCoords, p_Dvtx)
    call storage_getbase_int2D(rtria%h_IverticesAtElement, p_Iverts)
    call storage_getbase_int2D(rtria%h_IedgesAtElement, p_Iedges)
    
    ! create discretisation
    call spdiscr_initDiscr_simple (rdiscr, EL_EM30_NEW, CUB_G3_2D, rtria, rbnd)
    
    ! create two empty vectors
    call lsyssc_createVecByDiscr(rdiscr, ru1)
    call lsyssc_createVecByDiscr(rdiscr, ru2)

    ! create the matrix structure
    call bilf_createMatrixStructure (rdiscr, LSYSSC_MATRIX9, rmatrix)
    call lsyssc_allocEmptyMatrix(rmatrix, LSYSSC_SETM_UNDEFINED)
    
    ! Fetch all arrays from the system
    call lsyssc_getbase_double(ru1, p_Du1)
    call lsyssc_getbase_double(ru2, p_Du2)
    call lsyssc_getbase_double(rmatrix, p_DA)
    call lsyssc_getbase_Kld(rmatrix, p_Kld)
    call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
    
    ! Set u1 and u2 to (1,0)
    call lsyssc_clearVector(ru1, 1.0_DP)
    call lsyssc_clearVector(ru2, 0.0_DP)
    
    print *, 'FEAT2 Assembly'
    ! Assemble matrix using default SD
    
    call lsyssc_clearMatrix(rmatrix, 0.0_DP)
    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)
!    call sd_asm(p_Du1,p_Du2,rmatrix,ddelta,dnu,dalpha,dbeta,1.0_DP,dgamma,0)
!    call stdop_assembleSimpleMatrix(rmatrix, DER_FUNC, DER_FUNC,bclear=.false.)
!    call stdop_assembleLaplaceMatrix(rmatrix, .false.)
    print *, '-> done!'
    call stat_stopTimer(rtimer)
    print *, 'Time: ', stat_sgetTime_byTimer(rtimer)

    
    ! Write matrix to text file
!    call matio_writeMatrixHR (rmatrix, 'FEAT2 SD',&
!                              .true., 0, './mat_feat2.txt', '(E20.10)', 1e-5_DP)
    
    print *, 'Initialising GPU-SD...'
    ! initialise GPU-SD
    call gpusd_init(info)
    if(info .ne. 0) then
      print *, 'Failed to initialise GPU-SD!'
      stop
    end if
    print *, '-> done!'

    ! set scalar parameters
    call gpusd_setparams(dnu, dalpha, dbeta, dgamma, ddelta)
    
    ! set triangulation arrays
    call gpusd_settria(rtria%NEL, p_Dvtx, p_Iverts, p_Iedges)
    
    ! set matrix arrays
    call gpusd_setmatrix(rmatrix%NEQ, p_Kld, p_Kcol, p_DA)
    
    ! set velocity
    call gpusd_setvelo(p_Du1, p_Du2)
    
    print *, 'GPU-SD Assembly'
    ! try to assemble the matrix
    call lsyssc_clearMatrix(rmatrix, 0.0_DP)
    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)
    call gpusd_assemble(info)
    print *, '-> done! INFO: ', info
    call stat_stopTimer(rtimer)
    print *, 'Time: ', stat_sgetTime_byTimer(rtimer)
    call gpusd_statistics()

    ! release GPU-SD
    call gpusd_done()

    ! Write matrix to text file
!    call matio_writeMatrixHR (rmatrix, 'GPU-SD',&
!                              .true., 0, './mat_gpu.txt', '(E20.10)', 1e-5_DP)
    
    ! Release everything
    call lsyssc_releaseMatrix(rmatrix)
    call lsyssc_releaseVector(ru2)
    call lsyssc_releaseVector(ru1)
    call spdiscr_releaseDiscr(rdiscr)
    call tria_done(rtria)
    call boundary_release(rbnd)
  
end subroutine

end module