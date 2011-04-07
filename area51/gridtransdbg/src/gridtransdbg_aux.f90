
module gridtransdbg_aux

  use fsystem
  use storage
  use derivatives
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation
  use stdoperators
  use multileveloperators
  use multilevelprojection
  use quicksolver

  implicit none

  private

  public :: gtaux_densify
  public :: gtaux_asmProlMat_hc
  public :: gtaux_asmProlMat_l2
  public :: gtaux_asmProlMat_sl2
  public :: gtaux_asmRestMat_hc
  public :: gtaux_asmRestMat_l2
  public :: gtaux_asmRestMat_sl2
  public :: gtaux_asmInterpMat_hc
  public :: gtaux_asmInterpMat_l2
  public :: gtaux_asmInterpMat_sl2

contains

  ! ***************************************************************************

!<subroutine>
  
  subroutine gtaux_densify(rmatrix, p_Dmatrix, btranspose)

!<description>
  ! This routine densifies a type 7/9 sparse matrix, i.e. copies the matrix
  ! entries into a dense matrix.
!</description>

!<input>
  ! The scalar matrix which is to be densified.
  type(t_matrixScalar), intent(IN) :: rmatrix
  
  ! OPTIONAL:
  logical, optional, intent(IN) :: btranspose
!</input>

!<output>
  ! The dense ouput matrix.
  real(DP), dimension(:,:), pointer :: p_Dmatrix
!</output>

!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_DA
  integer, dimension(:), pointer :: p_Kld, p_Kcol
  integer :: i,j
  logical :: btrans
  
    ! Transpose?
    if(present(btranspose)) then
      btrans = btranspose
    else
      btrans = .false.
    end if
  
    ! Get the pointers from the storage
    call storage_getbase_double(rmatrix%h_DA,p_DA)
    call storage_getbase_int(rmatrix%h_Kld, p_Kld)
    call storage_getbase_int(rmatrix%h_Kcol, p_Kcol)
    
    if(btrans) then
    
      ! Allocate matrix if necessary
      if(.not. associated(p_Dmatrix)) &
        allocate(p_Dmatrix(rmatrix%NCOLS, rmatrix%NEQ))
      
      ! Format output matrix
      p_Dmatrix = 0.0_DP
      
      ! densify
      do i = 1, rmatrix%NEQ
        do j = p_Kld(i), p_Kld(i+1)-1
          p_Dmatrix(p_Kcol(j),i) = p_DA(j)
        end do ! j
      end do ! i
    
    else
    
      ! Allocate matrix if necessary
      if(.not. associated(p_Dmatrix)) &
        allocate(p_Dmatrix(rmatrix%NEQ, rmatrix%NCOLS))
      
      ! Format output matrix
      p_Dmatrix = 0.0_DP
      
      ! densify
      do i = 1, rmatrix%NEQ
        do j = p_Kld(i), p_Kld(i+1)-1
          p_Dmatrix(i,p_Kcol(j)) = p_DA(j)
        end do ! j
      end do ! i
    
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine gtaux_asmProlMat_hc(rdiscrC, rdiscrF, p_Dmat)

!<description>
  ! Assembles a dense prolongation matrix using hard-coded operators.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
!</input>

!<output>
  ! The prolongation matrix, must be set to NULL on entry.
  real(DP), dimension(:,:), pointer :: p_Dmat
!</output>

!</subroutine>

  ! Local variables
  type(t_vectorBlock) :: rvecC, rvecF
  type(t_vectorScalar) :: rtmp
  type(t_interlevelProjectionBlock) :: rproj
  real(DP), dimension(:), pointer :: p_DvecC, p_DvecF
  integer :: i,j,m,n
  
    ! Create the projection structure
    call mlprj_initProjectionDiscr(rproj, rdiscrF)
    
    ! Create two block vectors
    call lsysbl_createVecBlockByDiscr(rdiscrC, rvecC)
    call lsysbl_createVecBlockByDiscr(rdiscrF, rvecF)
    
    ! Get the data arrays
    call storage_getbase_double(rvecC%h_Ddata, p_DvecC)
    call storage_getbase_double(rvecF%h_Ddata, p_DvecF)
    
    ! Get the dimensions
    m = rvecC%NEQ
    n = rvecF%NEQ
    
    ! Create temporary vector
    i = mlprj_getTempMemoryVec (rproj,rvecC,rvecF)
    call lsyssc_createVector (rtmp,i,.false.)

    ! Allocate matrix
    allocate(p_Dmat(n,m))

    ! Go through all coarse mesh DOFs
    do i = 1, m
    
      ! Set coarse vector
      p_DvecC = 0.0_DP
      p_DvecC(i) = 1.0_DP
      
      ! Prolongate vector
      call mlprj_performProlongation(rproj, rvecC, rvecF, rtmp)
      
      ! Store result
      do j = 1, n
        p_Dmat(j,i) = p_DvecF(j)
      end do
    end do

    ! Release temporary vector
    if(rtmp%NEQ .gt. 0) then
      call lsyssc_releaseVector(rtmp)
    end if
    
    ! Release the vectors
    call lsysbl_releaseVector(rvecF)
    call lsysbl_releaseVector(rvecC)
    
    ! Release the projection structure
    call mlprj_doneProjection(rproj)
    
    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine gtaux_asmProlMat_l2(rdiscrC, rdiscrF, p_Dmat)

!<description>
  ! Assembles a dense prolongation matrix using L2-projection.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
!</input>

!<output>
  ! The prolongation matrix, must be set to NULL on entry.
  real(DP), dimension(:,:), pointer :: p_Dmat
!</output>

!</subroutine>

  ! Local variables
  type(t_vectorScalar) :: rvecC, rvecF
  type(t_matrixScalar) :: rmat2Lvl, rmatMass
  real(DP), dimension(:), pointer :: p_DvecC, p_DvecF
  real(DP), dimension(:), allocatable :: Dwork
  integer :: i,j,m,n,cinfo,niter
  real(DP) :: dtol
  
    ! Assemble a 2-Level mass matrix
    call mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), LSYSSC_MATRIX9, rmat2Lvl)
    call mlop_build2LvlMassMatrix(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), .true., rmat2Lvl)
    
    ! Assemble a mass matrix on the fine level
    call bilf_createMatrixStructure(rdiscrF%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9, rmatMass)
    call stdop_assembleSimpleMatrix(rmatMass, DER_FUNC, DER_FUNC)
    
    ! Create two block vectors
    call lsyssc_createVecByDiscr(rdiscrC%RspatialDiscr(1), rvecC)
    call lsyssc_createVecByDiscr(rdiscrF%RspatialDiscr(1), rvecF)
    
    ! Get the data arrays
    call storage_getbase_double(rvecC%h_Ddata, p_DvecC)
    call storage_getbase_double(rvecF%h_Ddata, p_DvecF)
    
    ! Get the dimensions
    m = rvecC%NEQ
    n = rvecF%NEQ
    
    ! Allocate matrix
    allocate(p_Dmat(n,m))
    
    ! Allocate work array
    allocate(Dwork(n))

    ! Go through all coarse mesh DOFs
    do i = 1, m
    
      ! Set coarse vector
      p_DvecC = 0.0_DP
      p_DvecC(i) = 1.0_DP
      
      ! Set up fine-mesh RHS vector
      call lsyssc_scalarMatVec(rmat2Lvl, rvecC, rvecF, 1.0_DP, 0.0_DP)
      
      ! Solve using SOR
      niter = 1000
      dtol = SYS_EPSREAL_DP
      call qsol_solveSOR(rmatMass, rvecF, Dwork, cinfo, niter, dtol, 1.2_DP)
      
      ! Store result
      do j = 1, n
        p_Dmat(j,i) = p_DvecF(j)
      end do
    end do
    
    ! Release work array
    deallocate(Dwork)
    
    ! Release the vectors
    call lsyssc_releaseVector(rvecF)
    call lsyssc_releaseVector(rvecC)
    
    ! Release the matrices
    call lsyssc_releaseMatrix(rmatMass)
    call lsyssc_releaseMatrix(rmat2Lvl)
    
    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine gtaux_asmProlMat_sl2(rdiscrC, rdiscrF, p_Dmat, cavrgType)

!<description>
  ! Assembles a dense prolongation matrix using sliced L2-projection.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
  
  ! The averaging type
  integer, intent(IN) :: cavrgType
!</input>

!<output>
  ! The prolongation matrix, must be set to NULL on entry.
  real(DP), dimension(:,:), pointer :: p_Dmat
!</output>

!</subroutine>

  ! Local variables
  type(t_matrixScalar) :: rmatProj
  
    ! Assemble a 2-Level prolongation matrix
    call mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), LSYSSC_MATRIX9, rmatProj)
    call mlop_build2LvlProlMatrix(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), .true., rmatProj, cavrgType)
    
    ! Densify the sparse matrix
    call gtaux_densify(rmatProj, p_Dmat, .false.)

    ! Release the matrix
    call lsyssc_releaseMatrix(rmatProj)
    
    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine gtaux_asmRestMat_hc(rdiscrC, rdiscrF, p_Dmat)

!<description>
  ! Assembles a dense restriction matrix using hard-coded operators.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
!</input>

!<output>
  ! The restriction matrix, must be set to NULL on entry.
  real(DP), dimension(:,:), pointer :: p_Dmat
!</output>

!</subroutine>

  ! Local variables
  type(t_vectorBlock) :: rvecC, rvecF
  type(t_vectorScalar) :: rtmp
  type(t_interlevelProjectionBlock) :: rproj
  real(DP), dimension(:), pointer :: p_DvecC, p_DvecF
  integer :: i,j,m,n
  
    ! Create the projection structure
    call mlprj_initProjectionDiscr(rproj, rdiscrF)
    
    ! Create two block vectors
    call lsysbl_createVecBlockByDiscr(rdiscrC, rvecC)
    call lsysbl_createVecBlockByDiscr(rdiscrF, rvecF)
    
    ! Get the data arrays
    call storage_getbase_double(rvecC%h_Ddata, p_DvecC)
    call storage_getbase_double(rvecF%h_Ddata, p_DvecF)
    
    ! Get the dimensions
    m = rvecC%NEQ
    n = rvecF%NEQ
    
    ! Create temporary vector
    i = mlprj_getTempMemoryVec (rproj,rvecC,rvecF)
    call lsyssc_createVector (rtmp,i,.false.)

    ! Allocate matrix
    allocate(p_Dmat(m,n))

    ! Go through all fine mesh DOFs
    do i = 1, n
    
      ! Set fine vector
      p_DvecF = 0.0_DP
      p_DvecF(i) = 1.0_DP
      
      ! restrict vector
      call mlprj_performRestriction(rproj, rvecC, rvecF, rtmp)
      
      ! Store result
      do j = 1, m
        p_Dmat(j,i) = p_DvecC(j)
      end do
    end do

    ! Release temporary vector
    if(rtmp%NEQ .gt. 0) then
      call lsyssc_releaseVector(rtmp)
    end if
    
    ! Release the vectors
    call lsysbl_releaseVector(rvecF)
    call lsysbl_releaseVector(rvecC)
    
    ! Release the projection structure
    call mlprj_doneProjection(rproj)
    
    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine gtaux_asmRestMat_l2(rdiscrC, rdiscrF, p_Dmat)

!<description>
  ! Assembles a dense restriction matrix using L2-projection.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
!</input>

!<output>
  ! The restriction matrix, must be set to NULL on entry.
  real(DP), dimension(:,:), pointer :: p_Dmat
!</output>

!</subroutine>

  ! Local variables
  type(t_vectorScalar) :: rvecC, rvecF
  type(t_matrixScalar) :: rmat2Lvl, rmatMass
  real(DP), dimension(:), pointer :: p_DvecC, p_DvecF
  real(DP), dimension(:), allocatable :: Dwork
  integer :: i,j,m,n,cinfo,niter
  real(DP) :: dtol
  
    ! Assemble a 2-Level mass matrix
    call mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), LSYSSC_MATRIX9, rmat2Lvl)
    call mlop_build2LvlMassMatrix(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), .true., rmat2Lvl)
    
    ! Assemble a mass matrix on the fine level
    call bilf_createMatrixStructure(rdiscrF%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9, rmatMass)
    call stdop_assembleSimpleMatrix(rmatMass, DER_FUNC, DER_FUNC)
    
    ! Create two block vectors
    call lsyssc_createVecByDiscr(rdiscrC%RspatialDiscr(1), rvecC)
    call lsyssc_createVecByDiscr(rdiscrF%RspatialDiscr(1), rvecF)
    
    ! Get the data arrays
    call storage_getbase_double(rvecC%h_Ddata, p_DvecC)
    call storage_getbase_double(rvecF%h_Ddata, p_DvecF)
    
    ! Get the dimensions
    m = rvecC%NEQ
    n = rvecF%NEQ
    
    ! Allocate matrix
    allocate(p_Dmat(m,n))
    
    ! Allocate work array
    allocate(Dwork(n))

    ! Go through all fine mesh DOFs
    do i = 1, n
    
      ! Set fine vector
      p_DvecF = 0.0_DP
      p_DvecF(i) = 1.0_DP
      
      ! Solve using SOR
      niter = 1000
      dtol = SYS_EPSREAL_DP
      call qsol_solveSOR(rmatMass, rvecF, Dwork, cinfo, niter, dtol, 1.2_DP)

      ! ...
      call lsyssc_scalarMatVec(rmat2Lvl, rvecF, rvecC, 1.0_DP, 0.0_DP, .true.)
      
      ! Store result
      do j = 1, m
        p_Dmat(j,i) = p_DvecC(j)
      end do
    end do
    
    ! Release work array
    deallocate(Dwork)
    
    ! Release the vectors
    call lsyssc_releaseVector(rvecF)
    call lsyssc_releaseVector(rvecC)
    
    ! Release the matrices
    call lsyssc_releaseMatrix(rmatMass)
    call lsyssc_releaseMatrix(rmat2Lvl)
    
    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine gtaux_asmRestMat_sl2(rdiscrC, rdiscrF, p_Dmat, cavrgType)

!<description>
  ! Assembles a dense restriction matrix using sliced L2-projection.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
  
  ! The averaging type
  integer, intent(IN) :: cavrgType
!</input>

!<output>
  ! The restriction matrix, must be set to NULL on entry.
  real(DP), dimension(:,:), pointer :: p_Dmat
!</output>

!</subroutine>

  ! Local variables
  type(t_matrixScalar) :: rmatProj
  
    ! Assemble a 2-Level prolongation matrix
    call mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), LSYSSC_MATRIX9, rmatProj)
    call mlop_build2LvlProlMatrix(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), .true., rmatProj, cavrgType)
    
    ! Densify the transposed sparse matrix
    call gtaux_densify(rmatProj, p_Dmat, .true.)

    ! Release the matrix
    call lsyssc_releaseMatrix(rmatProj)
    
    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine gtaux_asmInterpMat_hc(rdiscrC, rdiscrF, p_Dmat)

!<description>
  ! Assembles a dense interpolation matrix using hard-coded operators.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
!</input>

!<output>
  ! The interpolation matrix, must be set to NULL on entry.
  real(DP), dimension(:,:), pointer :: p_Dmat
!</output>

!</subroutine>

  ! Local variables
  type(t_vectorBlock) :: rvecC, rvecF
  type(t_vectorScalar) :: rtmp
  type(t_interlevelProjectionBlock) :: rproj
  real(DP), dimension(:), pointer :: p_DvecC, p_DvecF
  integer :: i,j,m,n
  
    ! Create the projection structure
    call mlprj_initProjectionDiscr(rproj, rdiscrF)
    
    ! Create two block vectors
    call lsysbl_createVecBlockByDiscr(rdiscrC, rvecC)
    call lsysbl_createVecBlockByDiscr(rdiscrF, rvecF)
    
    ! Get the data arrays
    call storage_getbase_double(rvecC%h_Ddata, p_DvecC)
    call storage_getbase_double(rvecF%h_Ddata, p_DvecF)
    
    ! Get the dimensions
    m = rvecC%NEQ
    n = rvecF%NEQ
    
    ! Create temporary vector
    i = mlprj_getTempMemoryVec (rproj,rvecC,rvecF)
    call lsyssc_createVector (rtmp,i,.false.)

    ! Allocate matrix
    allocate(p_Dmat(m,n))

    ! Go through all fine mesh DOFs
    do i = 1, n
    
      ! Set fine vector
      p_DvecF = 0.0_DP
      p_DvecF(i) = 1.0_DP
      
      ! interpolate vector
      call mlprj_performInterpolation(rproj, rvecC, rvecF, rtmp)
      
      ! Store result
      do j = 1, m
        p_Dmat(j,i) = p_DvecC(j)
      end do
    end do

    ! Release temporary vector
    if(rtmp%NEQ .gt. 0) then
      call lsyssc_releaseVector(rtmp)
    end if
    
    ! Release the vectors
    call lsysbl_releaseVector(rvecF)
    call lsysbl_releaseVector(rvecC)
    
    ! Release the projection structure
    call mlprj_doneProjection(rproj)
    
    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine gtaux_asmInterpMat_l2(rdiscrC, rdiscrF, p_Dmat)

!<description>
  ! Assembles a dense interpolation matrix using L2-projection.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
!</input>

!<output>
  ! The interpolation matrix, must be set to NULL on entry.
  real(DP), dimension(:,:), pointer :: p_Dmat
!</output>

!</subroutine>

  ! Local variables
  type(t_vectorScalar) :: rvecC, rvecF
  type(t_matrixScalar) :: rmat2Lvl, rmatMass
  real(DP), dimension(:), pointer :: p_DvecC, p_DvecF
  real(DP), dimension(:), allocatable :: Dwork
  integer :: i,j,m,n,cinfo,niter
  real(DP) :: dtol
  
    ! Assemble a 2-Level mass matrix
    call mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), LSYSSC_MATRIX9, rmat2Lvl)
    call mlop_build2LvlMassMatrix(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), .true., rmat2Lvl)
    
    ! Assemble a mass matrix on the coarse level
    call bilf_createMatrixStructure(rdiscrC%RspatialDiscr(1),&
                                    LSYSSC_MATRIX9, rmatMass)
    call stdop_assembleSimpleMatrix(rmatMass, DER_FUNC, DER_FUNC)
    
    ! Create two block vectors
    call lsyssc_createVecByDiscr(rdiscrC%RspatialDiscr(1), rvecC)
    call lsyssc_createVecByDiscr(rdiscrF%RspatialDiscr(1), rvecF)
    
    ! Get the data arrays
    call storage_getbase_double(rvecC%h_Ddata, p_DvecC)
    call storage_getbase_double(rvecF%h_Ddata, p_DvecF)
    
    ! Get the dimensions
    m = rvecC%NEQ
    n = rvecF%NEQ
    
    ! Allocate matrix
    allocate(p_Dmat(m,n))
    
    ! Allocate work array
    allocate(Dwork(m))

    ! Go through all fine mesh DOFs
    do i = 1, n
    
      ! Set fine vector
      p_DvecF = 0.0_DP
      p_DvecF(i) = 1.0_DP
      
      ! Set up caorse-mesh RHS vector
      call lsyssc_scalarMatVec(rmat2Lvl, rvecF, rvecC, 1.0_DP, 0.0_DP, .true.)
      
      ! Solve using SOR
      niter = 1000
      dtol = SYS_EPSREAL_DP
      call qsol_solveSOR(rmatMass, rvecC, Dwork, cinfo, niter, dtol, 1.2_DP)
      
      ! Store result
      do j = 1, m
        p_Dmat(j,i) = p_DvecC(j)
      end do
    end do
    
    ! Release work array
    deallocate(Dwork)
    
    ! Release the vectors
    call lsyssc_releaseVector(rvecF)
    call lsyssc_releaseVector(rvecC)
    
    ! Release the matrices
    call lsyssc_releaseMatrix(rmatMass)
    call lsyssc_releaseMatrix(rmat2Lvl)
    
    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine gtaux_asmInterpMat_sl2(rdiscrC, rdiscrF, p_Dmat, cavrgType)

!<description>
  ! Assembles a dense interpolation matrix using sliced L2-projection.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
  
  ! The averaging type
  integer, intent(IN) :: cavrgType
!</input>

!<output>
  ! The interpolation matrix, must be set to NULL on entry.
  real(DP), dimension(:,:), pointer :: p_Dmat
!</output>

!</subroutine>

  ! Local variables
  type(t_matrixScalar) :: rmatProj
  
    ! Assemble a 2-Level prolongation matrix
    call mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), LSYSSC_MATRIX9, rmatProj)
    call mlop_build2LvlInterpMatrix(rdiscrC%RspatialDiscr(1),&
        rdiscrF%RspatialDiscr(1), .true., rmatProj, cavrgType)
    
    ! Densify the sparse matrix
    call gtaux_densify(rmatProj, p_Dmat, .true.)

    ! Release the matrix
    call lsyssc_releaseMatrix(rmatProj)
    
    ! That's it

  end subroutine

end module
