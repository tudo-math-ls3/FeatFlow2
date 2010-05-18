program interpdbg

  use analyticprojection
  use bilinearformevaluation
  use boundary
  use collection
  use cubature
  use derivatives
  use element
  use fsystem
  use linearsystemblock
  use linearsystemscalar
  use pprocerror
  use scalarpde
  use spatialdiscretisation
  use stdoperators
  use storage
  use triangulation
  use ucd

  use interp_callback

  implicit none
  
  ! Triangulations
  type(t_triangulation) :: rtriangulation1, rtriangulation2

  ! Boundaries
  type(t_boundary) :: rboundary1, rboundary2
  
  ! Discretisation structures
  type(t_blockDiscretisation) :: rdiscretisation1, rdiscretisation2

  ! Matrices and vectors
  type(t_matrixScalar) :: rmatrixMassLumped1, rmatrixMassLumped2
  type(t_matrixBlock) :: rmatrix1, rmatrix2
  type(t_vectorBlock), target :: rvector1, rvector2

  ! UCD output structures
  type(t_ucdExport) :: rexport1, rexport2

  ! Refinement parameters
  integer :: nref1 = 0
  integer :: nref2 = 4

  ! Summed cubature parameters
  integer :: nsumcub1 = 0
  integer :: nsumcub2 = 5

  ! Collection structure
  type(t_collection) :: rcollection

  ! local variables
  real(DP), dimension(:), pointer :: p_Ddata, p_Da
  real(DP) :: derror
  integer :: i,j


  ! Set halt behaviour
  sys_haltmode = SYS_HALT_THROWFPE

  ! Initialise Featflow and storage subsystems
  call system_init()
  call storage_init(100, 100)

  ! Read boundary parametrisation
  call boundary_read_prm(rboundary1, 'pre/triangulation1.prm')
  call boundary_read_prm(rboundary2, 'pre/triangulation2.prm')

  ! Read basic triangulations from file
  call tria_readTriFile2D (rtriangulation1, 'pre/triangulation1.tri' , rboundary1)
  call tria_readTriFile2D (rtriangulation2, 'pre/triangulation2.tri' , rboundary2)

  ! Refine triangulations
  call tria_quickRefine2LevelOrdering (nref1, rtriangulation1, rboundary1)
  call tria_quickRefine2LevelOrdering (nref2, rtriangulation2, rboundary2)

  ! Create standard meshes from raw triangulation
  call tria_initStandardMeshFromRaw (rtriangulation1, rboundary1)
  call tria_initStandardMeshFromRaw (rtriangulation2, rboundary2)

  ! Initialse discretisations
  call spdiscr_initBlockDiscr(rdiscretisation1, 1, rtriangulation1)
  call spdiscr_initBlockDiscr(rdiscretisation2, 1, rtriangulation2)

  ! Create discretisations
  call spdiscr_initDiscr_triquad(rdiscretisation1%RspatialDiscr(1),&
      EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC,&
      rtriangulation1, rboundary1)
  call spdiscr_initDiscr_triquad(rdiscretisation2%RspatialDiscr(1),&
      EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC,&
      rtriangulation2, rboundary2)
  
  ! Enforce used of summed cubature formulas
  do i = 1, rdiscretisation1%ncomponents
    do j = 1, rdiscretisation1%RspatialDiscr(i)%inumFESpaces
      rdiscretisation1%RspatialDiscr(i)%RelementDistr(j)%ccubTypeLinForm =&
          cub_getSummedCubType(rdiscretisation1%RspatialDiscr(i)&
          %RelementDistr(j)%ccubTypeLinForm, nsumcub1)
      rdiscretisation1%RspatialDiscr(i)%RelementDistr(j)%ccubTypeEval =&
          cub_getSummedCubType(rdiscretisation1%RspatialDiscr(i)&
          %RelementDistr(j)%ccubTypeEval, nsumcub1)
    end do
  end do

  do i = 1, rdiscretisation2%ncomponents
    do j = 1, rdiscretisation2%RspatialDiscr(i)%inumFESpaces
      rdiscretisation2%RspatialDiscr(i)%RelementDistr(j)%ccubTypeLinForm =&
          cub_getSummedCubType(rdiscretisation2%RspatialDiscr(i)&
          %RelementDistr(j)%ccubTypeLinForm, nsumcub2)
      rdiscretisation2%RspatialDiscr(i)%RelementDistr(j)%ccubTypeEval =&
          cub_getSummedCubType(rdiscretisation2%RspatialDiscr(i)&
          %RelementDistr(j)%ccubTypeEval, nsumcub2)
    end do
  end do
  
  ! Create matrix structures
  call lsysbl_createMatBlockByDiscr(rdiscretisation1, rmatrix1)
  call bilf_createMatrixStructure(rdiscretisation1%RspatialDiscr(1),&
      LSYSSC_MATRIX9, rmatrix1%RmatrixBlock(1,1))
  call lsysbl_createMatBlockByDiscr(rdiscretisation2, rmatrix2)
  call bilf_createMatrixStructure(rdiscretisation2%RspatialDiscr(1),&
      LSYSSC_MATRIX9, rmatrix2%RmatrixBlock(1,1))

  ! Create consistent mass matrices
  call stdop_assembleSimpleMatrix(rmatrix1%RmatrixBlock(1,1), DER_FUNC, DER_FUNC)
  call stdop_assembleSimpleMatrix(rmatrix2%RmatrixBlock(1,1), DER_FUNC, DER_FUNC)

  ! Create lumped mass matrices
  call lsyssc_duplicateMatrix(rmatrix1%RmatrixBlock(1,1), rmatrixMassLumped1,&
      LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
  call lsyssc_lumpMatrixScalar(rmatrixMassLumped1, LSYSSC_LUMP_DIAG)
  call lsyssc_duplicateMatrix(rmatrix2%RmatrixBlock(1,1), rmatrixMassLumped2,&
      LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
  call lsyssc_lumpMatrixScalar(rmatrixMassLumped2, LSYSSC_LUMP_DIAG)

  ! Initialise the solution vector1
  call lsysbl_createVectorBlock(rdiscretisation1, rvector1, .true.)
  call interp_initSolution(rvector1)
  
  ! Attach solution vector to collection structure
  call collct_init(rcollection)
  rcollection%p_rvectorQuickAccess1 => rvector1
  
  ! Project the solution onto rvector2
  call lsysbl_createVectorBlock(rdiscretisation2, rvector2, .true.)

!!$  ! Lumped L2-projection
!!$  call anprj_analytL2projectionByMass(rvector2%RvectorBlock(1), &
!!$      rmatrixMassLumped2, interp_buildVector, rcollection,&
!!$      rmatrixMassLumped = rmatrixMassLumped2)
  
!!$  ! Consistent L2-projection
!!$  call anprj_analytL2projectionByMass(rvector2%RvectorBlock(1), &
!!$      rmatrix2%RmatrixBlock(1,1), interp_buildVector, rcollection,&
!!$      rmatrixMassLumped = rmatrixMassLumped2)
 
  ! Constrained L2-projection
  call anprj_analytL2projectionConstr(rvector2%RvectorBlock(1), &
      rmatrix2%RmatrixBlock(1,1), interp_buildVector, rcollection,&
      rmatrixMassLumped = rmatrixMassLumped2)
 
  ! Output solutions
  call ucd_startGMV (rexport1, UCD_FLAG_STANDARD, rtriangulation1, 'out/out1.gmv')
  call ucd_startGMV (rexport2, UCD_FLAG_STANDARD, rtriangulation2, 'out/out2.gmv')

  call lsysbl_getbase_double (rvector1, p_Ddata)
  call ucd_addVariableVertexBased (rexport1, 'sol', UCD_VAR_STANDARD, p_Ddata)
  
  call lsysbl_getbase_double (rvector2, p_Ddata)
  call ucd_addVariableVertexBased (rexport2, 'sol', UCD_VAR_STANDARD, p_Ddata)
  
  ! Write the file to disc
  call ucd_write(rexport1)
  call ucd_write(rexport2)
  call ucd_release(rexport1)
  call ucd_release(rexport2)

  ! Compute statistical information
  call lsyssc_getbase_double(rmatrixMassLumped1, p_Da)
  call lsysbl_getbase_double(rvector1, p_Ddata)
  write (*,fmt='(A,X,G15.6)') 'Mass sol1', sum(p_Da*p_Ddata)
  write (*,fmt='(A,X,G15.6)') 'Min  sol1', minval(p_Ddata)
  write (*,fmt='(A,X,G15.6)') 'Max  sol1', maxval(p_Ddata)
  
  call pperr_scalar(rvector1%RvectorBlock(1), PPERR_L1ERROR, derror,&
      interp_refFunction, rcollection)
  write (*,fmt='(A,X,G15.6)') 'L1-error1', derror

  call pperr_scalar(rvector1%RvectorBlock(1), PPERR_L2ERROR, derror,&
      interp_refFunction, rcollection)
  write (*,fmt='(A,X,G15.6)') 'L2-error1', derror
  
  call lsyssc_getbase_double(rmatrixMassLumped2, p_Da)
  call lsysbl_getbase_double(rvector2, p_Ddata)
  write (*,fmt='(A,X,G15.6)') 'Mass sol2', sum(p_Da*p_Ddata)
  write (*,fmt='(A,X,G15.6)') 'Min  sol2', minval(p_Ddata)
  write (*,fmt='(A,X,G15.6)') 'Max  sol2', maxval(p_Ddata)

  call pperr_scalar(rvector2%RvectorBlock(1), PPERR_L1ERROR, derror,&
      interp_refFunction, rcollection)
  write (*,fmt='(A,X,G15.6)') 'L1-error2', derror

  call pperr_scalar(rvector2%RvectorBlock(1), PPERR_L2ERROR, derror,&
      interp_refFunction, rcollection)
  write (*,fmt='(A,X,G15.6)') 'L2-error2', derror
    
  ! Release all data structures
  call lsysbl_releaseVector(rvector1)
  call lsysbl_releaseVector(rvector2)
  call lsysbl_releaseMatrix(rmatrix1)
  call lsysbl_releaseMatrix(rmatrix2)
  call lsyssc_releaseMatrix(rmatrixMassLumped1)
  call lsyssc_releaseMatrix(rmatrixMassLumped2)
  
  call tria_done(rtriangulation1)
  call tria_done(rtriangulation2)

  call spdiscr_releaseBlockDiscr(rdiscretisation1)
  call spdiscr_releaseBlockDiscr(rdiscretisation2)

  call boundary_release(rboundary1)
  call boundary_release(rboundary2)

  call collct_done(rcollection)

  ! Clean up the storage management, finish
  call storage_info(.true.)
  call storage_done()

contains

  !*****************************************************************************

  subroutine interp_initSolution(rvector)

    type(t_vectorBlock), intent(inout) :: rvector
    
    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: ivt

    ! Get pointers
    p_rtriangulation => rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    call lsysbl_getbase_double(rvector, p_Ddata)
    call storage_getbase_double2d(p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
   
    do ivt = 1, size(p_DvertexCoords, 2)
      if (sqrt((p_DvertexCoords(1,ivt)-0.5_DP)**2+&
               (p_DvertexCoords(2,ivt)-0.5_DP)**2) .le. 0.3) then
        p_Ddata(ivt) = 1.0_DP
      else
        p_Ddata(ivt) = 0.01_DP
      end if
    end do
  end subroutine interp_initSolution

  
   
end program interpdbg
