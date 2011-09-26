!##############################################################################
!# ****************************************************************************
!# <name> spacetimematvec </name>
!# ****************************************************************************
!#
!# <purpose>
!# Application of and with space-time matrices.
!# </purpose>
!##############################################################################

module spacetimematvec

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
  use timediscretisation
  use matrixfilters
  use matrixio
  
  use stdoperators
  use bilinearformevaluation
  
  use spatialoperators
  use spacetimevectors
  
  use physics
  use spacetimebc
  
  use fespacehierarchy
  use spacetimehierarchy
  
  use spacetimematrices
  
  implicit none

contains

  ! ***************************************************************************

  subroutine stmv_matvec (rmatrix, rx, ry, cx, cy)
  
  ! Matrix-vector multiplication with a space-time matrix:
  ! ry = cx*A*rx + cy*ry
  
  ! Underlying space-time matrix
  type(t_spaceTimeMatrix), intent(in) :: rmatrix
  
  ! Input vector x.
  type(t_spaceTimeVector), intent(in) :: rx
  
  ! Input and output vector y.
  type(t_spaceTimeVector), intent(inout) :: ry
  
  ! Weight for x.
  real(DP), intent(in) :: cx
  
  ! Weioght for y
  real(DP), intent(in) :: cy
  
    ! local variables
    integer :: istep,ispaceLevel
    type(t_vectorBlock) :: rtempVecY,rtempVecX1,rtempVecX2,rtempVecX3
    type(t_matrixBlock) :: rsubmatrix
    real(DP) :: dnormy
    real(DP), dimension(:), pointer :: p_Dx1, p_Dx2, p_Dx3, p_Dy
    type(t_feSpaceLevel), pointer :: p_rfeSpaceLevel
    
    ! Level in space of the space-time matrix
    call sth_getLevel(rmatrix%p_rspaceTimeHierarchy,rmatrix%ilevel,&
        p_rfeSpaceLevel=p_rfeSpaceLevel,ispaceLevel=ispaceLevel)
    
    ! Allocate a temp vectors and matrices
    call lsysbl_createVectorBlock(p_rfeSpaceLevel%p_rdiscretisation,rtempVecX1)
    call lsysbl_createVectorBlock(p_rfeSpaceLevel%p_rdiscretisation,rtempVecX2)
    call lsysbl_createVectorBlock(p_rfeSpaceLevel%p_rdiscretisation,rtempVecX3)
    call lsysbl_createVectorBlock(p_rfeSpaceLevel%p_rdiscretisation,rtempVecY)
    
    call lsysbl_getbase_double (rtempVecX1,p_Dx1)
    call lsysbl_getbase_double (rtempVecX2,p_Dx2)
    call lsysbl_getbase_double (rtempVecX3,p_Dx3)
    call lsysbl_getbase_double (rtempVecY,p_Dy)
    
    call stmat_allocSubmatrix (rmatrix%cmatrixType,rmatrix%p_rphysics,&
        rmatrix%p_rmatVecTempl(ispaceLevel),rsubmatrix)
    
    ! Loop over all timesteps
    do istep=1,rx%NEQtime
    
      ! Load the RHS into the temp vector. Scale it by cy; we cannot 
      ! include cy in the space-time MV below since tghe weights would
      ! be wrong for subsequent MV's otherwise.
      call sptivec_getTimestepData(ry,istep,rtempVecY)
      call lsysbl_scaleVector (rtempVecY,cy)
      
      ! There are three MV's in each timestep...
      
      if (istep .gt. 1) then
        ! Left subdiagonal
        call sptivec_getTimestepData(rx,istep-1,rtempVecX1)
        call stmat_getSubmatrix (rmatrix, ispaceLevel, istep, istep-1, rsubmatrix)
        call lsysbl_blockMatVec (rsubmatrix,rtempVecX1,rtempVecY,cx,1.0_DP)
      end if
      
      ! Diagonal
      call sptivec_getTimestepData(rx,istep,rtempVecX2)
      call stmat_getSubmatrix (rmatrix, ispaceLevel, istep, istep, rsubmatrix)
      call lsysbl_blockMatVec (rsubmatrix,rtempVecX2,rtempVecY,cx,1.0_DP)
      
!      call matio_writeBlockMatrixHR (rsubmatrix, "matrix",&
!          .true., 0, "matrix", "(E20.10)")

      if (istep .lt. rx%NEQtime) then
        ! Right subdiagonal
        call sptivec_getTimestepData(rx,istep+1,rtempVecX3)
        call stmat_getSubmatrix (rmatrix, ispaceLevel, istep, istep+1, rsubmatrix)
        call lsysbl_blockMatVec (rsubmatrix,rtempVecX3,rtempVecY,cx,1.0_DP)
      end if
      
      call sptivec_setTimestepData(ry,istep,rtempVecY)

      dnormy = lsysbl_vectorNorm(rtempVecY,LINALG_NORML2)

    end do
    
    ! Release the temp vector.
    call lsysbl_releaseMatrix (rsubmatrix)
    call lsysbl_releaseVector (rtempVecX3)
    call lsysbl_releaseVector (rtempVecX2)
    call lsysbl_releaseVector (rtempVecX1)
    call lsysbl_releaseVector (rtempVecY)
  
  end subroutine

end module
