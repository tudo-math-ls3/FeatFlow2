!##############################################################################
!# ****************************************************************************
!# <name> spacetimevanca </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

MODULE spacetimevanca

  USE fsystem
  USE genoutput
  USE externalstorage
  USE spatialdiscretisation
  USE linearsystemscalar
  USE linearsystemblock
  USE collection
  USE vectorio
  
  USE cc2dmediumm2nonlinearcoreinit
  USE spacetimediscretisation
  USE spacetimevectors
  USE spacetimelinearsystem
  
  USE matrixio
  
  IMPLICIT NONE

  ! A structure that collects the pointers to all submatrices in a matrix
  ! row of the global block matrix. Circumvents the problem
  ! that Fortran cannot set up an array of pointers...
  TYPE t_matrixRow
  
    ! Pointers for submatrices below the diagonal
    REAL(DP), DIMENSION(:), POINTER :: p_Dp11 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Dp21 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Dp12 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Dp22 => NULL()
                                              
    ! Pointers for submatrices of a diagonal bock
    REAL(DP), DIMENSION(:), POINTER :: p_Da11 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da21 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da31 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da41 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da51 => NULL()
                                              
    REAL(DP), DIMENSION(:), POINTER :: p_Da12 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da22 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da32 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da42 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da52 => NULL()
                                              
    REAL(DP), DIMENSION(:), POINTER :: p_Da13 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da23 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da33 => NULL()
                                              
    REAL(DP), DIMENSION(:), POINTER :: p_Da14 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da24 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da44 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da54 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da64 => NULL()
                                              
    REAL(DP), DIMENSION(:), POINTER :: p_Da15 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da25 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da45 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da55 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da65 => NULL()
                                              
    REAL(DP), DIMENSION(:), POINTER :: p_Da46 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da56 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Da66 => NULL()
                                              
    ! Pointers for submatrices above the diagonal
    REAL(DP), DIMENSION(:), POINTER :: p_Dd44 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Dd54 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Dd45 => NULL()
    REAL(DP), DIMENSION(:), POINTER :: p_Dd55 => NULL()

    ! Scaling factors for all the submatrices.
    ! DIMENSION(6,6,-1..1), with 6=#matrix rows/columns. The last dimension specifies the
    ! matrix column: -1=left from the diagonal, 0=diagonal, 1=right from the diagonal.
    REAL(DP), DIMENSION(6,6,-1:1) :: DscaleFactor = 0.0_DP
    
  END TYPE
  
  ! Collects the positions of the DOF's on one element in all the matrices
  TYPE t_matrixPositions
  
    ! Matrix positions in the FEM matrices.
    INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: p_IposA  => NULL()
    
    ! Matrix positions in the B-matrices.
    INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: p_IposB  => NULL()
    
    ! matrix positions in the C matrix.
    INTEGER(PREC_MATIDX), DIMENSION(:,:), POINTER :: p_IposC  => NULL()
  
  END TYPE
  
  ! A structure that collects the pointers to all subvectors 
  ! of a solution chunk. Circumvents the problem
  ! that Fortran cannot set up an array of pointers...
  TYPE t_chunkVector
    
    ! Pointer to solution and RHS vector
    REAL(DP), DIMENSION(:), POINTER :: p_Dx
    REAL(DP), DIMENSION(:), POINTER :: p_Db
    
    ! Pointer to the subvectors of the solution vector
    REAL(DP), DIMENSION(:), POINTER :: p_Dx1
    REAL(DP), DIMENSION(:), POINTER :: p_Dx2
    REAL(DP), DIMENSION(:), POINTER :: p_Dx3
    REAL(DP), DIMENSION(:), POINTER :: p_Dx4
    REAL(DP), DIMENSION(:), POINTER :: p_Dx5
    REAL(DP), DIMENSION(:), POINTER :: p_Dx6
    
  END TYPE

    
CONTAINS

  SUBROUTINE createLocalMatrixQ1TQ0 (rdiscretisation,ntimesteps,rmatrix,h_ImatrixMapping)
  
  ! Creates a local matrix containing the submatrices of ntimesteps timesteps.
    
  ! Block discretisation structure of every timestep
  TYPE(t_blockDiscretisation), INTENT(IN) :: rdiscretisation
  
  ! Number of timesteps, the local matrix should hold.
  INTEGER, INTENT(IN) :: ntimesteps

  ! The local matrix to create
  TYPE(t_matrixScalar), INTENT(OUT) :: rmatrix
  
  ! A handle to an array with DIMENSION(:,:). For every matrix block column in the
  ! local matrix, this points to the first entry of the column. This array is
  ! therefore a way to directly access the block columns in the matrix.
  INTEGER, INTENT(OUT) :: h_ImatrixMapping
  
    ! local variables
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_Iindex
    INTEGER :: ndofLocal,i,j,ielType,ndofGlobal,itimestep,ndofLocalPrimalDual
    INTEGER :: ndofVelocity,ndofPressure, iYoffset,isum1,isum2
    INTEGER, DIMENSION(2) :: Isize
    REAL(DP), DIMENSION(:), POINTER :: p_Da
    
    ! Get the number of local DOF's in every timestep. For that purpose,
    ! loop over the blocks in the discretisation structure, ask the
    ! element about its number of local DOF's and sum that stuff up.
    ndofLocal = 0
    DO i=1,rdiscretisation%ncomponents
      ielType = rdiscretisation%RspatialDiscretisation(i)%&
          RelementDistribution(1)%itrialElement
      ndofLocal = ndofLocal + elem_igetNDofLoc(ielType)
    END DO
    
    ! Get the number of local DOF's in the primal space. This coincodes with
    ! the number of local DOF's in the dual space.
    ! As they both summed up give the number of local DOF's, we simply have to
    ! divide that by 2.
    ndofLocalPrimalDual = ndofLocal / 2
    
    ! Get the number of velocity DOF's in the primal space.
    ielType = rdiscretisation%RspatialDiscretisation(1)%&
        RelementDistribution(1)%itrialElement
    ndofVelocity = elem_igetNDofLoc(ielType)

    ielType = rdiscretisation%RspatialDiscretisation(3)%&
        RelementDistribution(1)%itrialElement
    ndofPressure = elem_igetNDofLoc(ielType)
    
    ! Multiplying that with the number of timesteps (+1) gives the size of the
    ! local matrix.
    ndofGlobal = ndofLocal * (ntimesteps + 1)
    
    ! Allocate memory for the pointer array that later allows quick access
    ! to the matrix columns. Each matrix is a 6x6 block matrix.
    Isize = (/6*(ntimesteps + 1),6*3/)
    CALL storage_new2D ('createLocalMatrixQ1TQ0', 'ImatrixMapping', Isize, &
        ST_INT, h_ImatrixMapping,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int2d(h_ImatrixMapping,p_Iindex)
    
    ! A local matrix has in our case the following shape
    ! (example for two timesteps)
    !
    !   A A A A  A A A A  B  M M M M
    !   A A A A  A A A A  B  M M M M
    !   A A A A  A A A A  B  M M M M
    !   A A A A  A A A A  B  M M M M
    !
    !   A A A A  A A A A  B           M M M M
    !   A A A A  A A A A  B           M M M M
    !   A A A A  A A A A  B           M M M M
    !   A A A A  A A A A  B           M M M M
    !
    !   B B B B  B B B B  C
    !
    !   R R R R  R R R R     A A A A  A A A A  B                       M M M M  M M M M
    !   R R R R  R R R R     A A A A  A A A A  B                       M M M M  M M M M
    !   R R R R  R R R R     A A A A  A A A A  B                       M M M M  M M M M
    !   R R R R  R R R R     A A A A  A A A A  B                       M M M M  M M M M
    !                                                                 
    !   R R R R  R R R R     A A A A  A A A A  B                       M M M M  M M M M
    !   R R R R  R R R R     A A A A  A A A A  B                       M M M M  M M M M
    !   R R R R  R R R R     A A A A  A A A A  B                       M M M M  M M M M
    !   R R R R  R R R R     A A A A  A A A A  B                       M M M M  M M M M
    !                       
    !                        B B B B  B B B B  C
    !
    !   M M M M  M M M M                          A A A A  A A A A  B  M M M M   
    !   M M M M  M M M M                          A A A A  A A A A  B  M M M M
    !   M M M M  M M M M                          A A A A  A A A A  B  M M M M
    !   M M M M  M M M M                          A A A A  A A A A  B  M M M M
    !                                            
    !   M M M M  M M M M                          A A A A  A A A A  B           M M M M
    !   M M M M  M M M M                          A A A A  A A A A  B           M M M M
    !   M M M M  M M M M                          A A A A  A A A A  B           M M M M
    !   M M M M  M M M M                          A A A A  A A A A  B           M M M M
    !                                            
    !                                             B B B B  B B B B  C
    !                                            
    !                                             R R R R  R R R R     A A A A  A A A A  B 
    !                                             R R R R  R R R R     A A A A  A A A A  B 
    !                                             R R R R  R R R R     A A A A  A A A A  B 
    !                                             R R R R  R R R R     A A A A  A A A A  B 
    !                                                                 
    !                                             R R R R  R R R R     A A A A  A A A A  B 
    !                                             R R R R  R R R R     A A A A  A A A A  B 
    !                                             R R R R  R R R R     A A A A  A A A A  B 
    !                                             R R R R  R R R R     A A A A  A A A A  B 
    !                                                                 
    !                                                                  B B B B  B B B B  C
    !
    ! Allocate a ndofGlobal*ndofGlobal-matrix in structure D, set up the above
    ! structure and convert it to structure 9. In the matrix entries, we save
    ! nonzero values; the conversion routine will later remove all zeroes.
    
    CALL lsyssc_createFullMatrix (rmatrix,.TRUE.,ndofGlobal)
    
    ! Get the data array
    CALL lsyssc_getbase_double (rmatrix,p_Da)
    
    ! Loop through all diagonal blocks to fill in data.
    DO itimestep = 0,ntimesteps
      ! Primal matrix
      CALL markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal,&
          1+itimestep*ndofLocal,&
          ndofLocalPrimalDual)

      ! Dual matrix
      CALL markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual,&
          ndofLocalPrimalDual)

      ! Mass matrices above the diagonal
      CALL markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual,&
          ndofVelocity)
      CALL markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal+ndofVelocity,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual+ndofVelocity,&
          ndofVelocity)

      ! The R-block below the diagonal
      CALL markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual,&
          1+itimestep*ndofLocal,&
          2*ndofVelocity)
          
    END DO
    
    ! What is still missing are the mass matrices that couple the timesteps.
    ! Loop through all diagonal blocks to fill in data.
    DO itimestep = 0,ntimesteps-1
      ! Primal matrix
      CALL markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal+ndofLocal,&
          1+itimestep*ndofLocal,&
          2*ndofVelocity)
      !CALL markSubmatrix (p_Da,ndofGlobal,&
      !    1+itimestep*ndofLocal+ndofLocal+ndofVelocity,&
      !    1+itimestep*ndofLocal+ndofVelocity,&
      !    2*ndofVelocity)

      ! Dual matrix
      CALL markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual,&
          1+itimestep*ndofLocal+ndofLocal+ndofLocalPrimalDual,&
          2*ndofVelocity)
      !CALL markSubmatrix (p_Da,ndofGlobal,&
      !    1+itimestep*ndofLocal+ndofLocalPrimalDual+ndofVelocity,&
      !    1+itimestep*ndofLocal+ndofLocal+ndofLocalPrimalDual+ndofVelocity,&
      !    ndofVelocity)
    END DO
    
    ! Convert the matrix to format 9
    CALL lsyssc_convertMatrix (rmatrix,LSYSSC_MATRIX9)
    
    ! Set up the matrix column pointers. FOr that purpose, at first
    ! write into every block the size of the block -- if it contains data.
    DO itimestep = 0,ntimesteps
      
      iYoffset = 6*itimestep

      IF (itimestep .GT. 0) THEN
        ! Left' matrix below the diagonal
        p_Iindex(iYoffset+1,0+1) = ndofVelocity
        p_Iindex(iYoffset+1,0+2) = ndofVelocity
                                     
        p_Iindex(iYoffset+2,0+1) = ndofVelocity
        p_Iindex(iYoffset+2,0+2) = ndofVelocity
      END IF

      ! 'Middle' matrix
      p_Iindex(iYoffset+1,6+1) = ndofVelocity 
      p_Iindex(iYoffset+1,6+2) = ndofVelocity 
      p_Iindex(iYoffset+1,6+3) = ndofPressure 
      p_Iindex(iYoffset+1,6+4) = ndofVelocity 
      !p_Iindex(iYoffset+1,6+5) = ndofVelocity 
                                    
      p_Iindex(iYoffset+2,6+1) = ndofVelocity 
      p_Iindex(iYoffset+2,6+2) = ndofVelocity 
      p_Iindex(iYoffset+2,6+3) = ndofPressure 
      !p_Iindex(iYoffset+2,6+4) = ndofVelocity 
      p_Iindex(iYoffset+2,6+5) = ndofVelocity 
                                    
      p_Iindex(iYoffset+3,6+1) = ndofVelocity 
      p_Iindex(iYoffset+3,6+2) = ndofVelocity 
      p_Iindex(iYoffset+3,6+3) = ndofPressure
                                    
      p_Iindex(iYoffset+4,6+1) = ndofVelocity 
      p_Iindex(iYoffset+4,6+2) = ndofVelocity 
                                    
      p_Iindex(iYoffset+4,6+4) = ndofVelocity 
      p_Iindex(iYoffset+4,6+5) = ndofVelocity 
      p_Iindex(iYoffset+4,6+6) = ndofPressure 
                                    
      p_Iindex(iYoffset+5,6+1) = ndofVelocity 
      p_Iindex(iYoffset+5,6+2) = ndofVelocity 
                                    
      p_Iindex(iYoffset+5,6+4) = ndofVelocity 
      p_Iindex(iYoffset+5,6+5) = ndofVelocity 
      p_Iindex(iYoffset+5,6+6) = ndofPressure 
                                    
      p_Iindex(iYoffset+6,6+4) = ndofVelocity 
      p_Iindex(iYoffset+6,6+5) = ndofVelocity 
      p_Iindex(iYoffset+6,6+6) = ndofPressure 
                
      IF (itimestep .LT. ntimesteps) THEN
        ! 'Right' matrix, above the diagonal
        p_Iindex(iYoffset+4,12+4) = ndofVelocity
        p_Iindex(iYoffset+4,12+5) = ndofVelocity
                                     
        p_Iindex(iYoffset+5,12+4) = ndofVelocity
        p_Iindex(iYoffset+5,12+5) = ndofVelocity
      END IF
    END DO
    
    ! Now loop through all the indices and sum them up.
    ! This gives the matrix pointers in the local matrices.
    !
    ! First shift the indices and set the first column to 0.
    DO j=6*3,2,-1
      DO i=1,UBOUND(p_Iindex,1)
        p_Iindex(i,j) = p_Iindex(i,j-1)
      END DO
    END DO
    
    DO i=1,UBOUND(p_Iindex,1)
      p_Iindex(i,1) = 0.0_DP
    END DO
    
    ! Now sum that stuff up
    DO j=2,6*3
      DO i=1,UBOUND(p_Iindex,1)
        p_Iindex(i,j) = p_Iindex(i,j) + p_Iindex(i,j-1)
      END DO
    END DO
    
    ! And add a 1. That gives the start pointers for the blocks.
    DO j=1,6*3
      DO i=1,UBOUND(p_Iindex,1)
        p_Iindex(i,j) = 1 + p_Iindex(i,j)
      END DO
    END DO
    
  CONTAINS
  
    SUBROUTINE markSubmatrix (Da,neq,iposY,iposX,isize)
    
    ! Initialises a subblock if the matrix Da with 1.0.
    ! The block starts at position (iposY,iposX) and has size isize*isize.
    
    ! Data array
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: Da
    
    ! Size of the matrix
    INTEGER, INTENT(IN) :: neq
    
    ! Start Position where to fill the matrix with 1.
    INTEGER, INTENT(IN) :: iposY,iposX
    
    ! Size of the block that should be filled with 1.
    INTEGER,INTENT(IN) :: isize
    
      ! local variables
      INTEGER :: i,j
      
      DO j=iposX-1,iposX+isize-2
        DO i=iposY,iposY+isize-1
          Da(j*neq+i) = 1.0_DP
        END DO
      END DO
      
    END SUBROUTINE
  
  END SUBROUTINE

  ! ***************************************************************************

  SUBROUTINE calcMatrixIndices (Kcol,Kld,IdofsTrial,IdofsTest,Iidx)
  
  ! Calculates the indices inside of those entries in a matrix that 
  ! belong to the local matrix. The matrix must be given in
  ! matrix structure 9.
  
  ! The column structure of the matrix.
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol

  ! The row structure of the matrix.
  INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld
  
  ! The DOF's in the trial space that belong to the local matrix.
  ! These correspond to the columns of the local matrix.
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IdofsTrial

  ! The DOF's in the test space that belong to the local matrix.
  ! THese correspond to the rows of the local matrix.
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IdofsTest
  
  ! A 2D array containing the positions in the global matrix that
  ! belong to the local matrix. The index array is set up transposed
  ! to get quicker memory access!
  ! Iidx must be a square matrix with 
  ! DIMENSION(size(IdofsTest),size(IdofsTrial))
  INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(OUT) :: Iidx
  
    ! local variables
    INTEGER :: i,k
    INTEGER(PREC_MATIDX) :: j
    INTEGER(PREC_VECIDX) :: idoflocal,icol
    
    ! Loop over the matrix rows.
    DO i=1,SIZE(IdofsTrial)
      
      idofLocal = IdofsTrial(i)
      
      ! In the row loop over the columns to find those that belong to 
      ! local DOF's.
      DO j=Kld(idofLocal),Kld(idofLocal+1)-1
      
        ! Column number
        icol = Kcol(j)
      
        ! Find the column in the Idofs array. If we find it, transfer its
        ! position to Iidx.
        DO k=1,SIZE(IdofsTest)
        
          IF (icol .EQ. IdofsTest(k)) THEN
            ! Put the index position to the local matrix
            Iidx(k,i) = j
            EXIT
          END IF
        
        END DO
      
      END DO
    
    END DO
    
  END SUBROUTINE
  
  
  
  SUBROUTINE ccopt_precSpaceTimeVanca (rproblem,rmatrix,rx,rb,domega,nchunkSize,niterations)
  
  ! Performs niterations VANCA iterations to enhance the vector rx.
  ! nchunkSize specifies the number of timesteps that are collected
  ! to a local system.
  
  ! The main problem structure
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  
  ! The underlying space time matrix.
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN) :: rmatrix
  
  ! RHS vector
  TYPE(t_spacetimeVector), INTENT(IN) :: rb
  
  ! Damping parameter
  REAL(DP), INTENT(IN) :: domega
  
  ! Chunk size. >= 0. A chunk size >= #timesteps will apply VANCA 
  ! to all timesteps simultaneously, which is VERY memory and hard disc
  ! intensive! A chunk size of 0 will apply the standard VANCA for
  ! optimal control problems that does not combine multiple time steps.
  INTEGER, INTENT(IN) :: nchunksize
  
  ! Number of VANCA iterations to apply.
  INTEGER, INTENT(IN) :: niterations
  
  ! Solution vector to be updated
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rx
  
    ! Constants
    INTEGER, PARAMETER :: ndofLocalVelocityX = 4
    INTEGER, PARAMETER :: ndofLocalVelocityY = 4
    INTEGER, PARAMETER :: ndofLocalVelocity  = ndofLocalVelocityX + ndofLocalVelocityY
    INTEGER, PARAMETER :: ndofLocalPressure  = 1
    INTEGER, PARAMETER :: ndofLocalHalf = ndofLocalVelocity + ndofLocalPressure
    INTEGER, PARAMETER :: ndofLocal = 2*ndofLocalHalf   ! primal and dual
    
    INTEGER, PARAMETER :: iposLocalPrimalVelocityX = 1
    INTEGER, PARAMETER :: iposLocalPrimalVelocityY = &
                              iposLocalPrimalVelocityX + ndofLocalVelocityX
    INTEGER, PARAMETER :: iposLocalPrimalPressure = &
                              iposLocalPrimalVelocityY + ndofLocalVelocityY

    INTEGER, PARAMETER :: iposLocalDualVelocityX = &
                              iposLocalPrimalPressure + ndofLocalPressure
    INTEGER, PARAMETER :: iposLocalDualVelocityY = &
                              iposLocalDualVelocityX + ndofLocalVelocityX
    INTEGER, PARAMETER :: iposLocalDualPressure = &
                              iposLocalDualVelocityY + ndofLocalVelocityY

    ! local variables
    INTEGER(PREC_VECIDX) :: iposGlobalPrimalVelocityX,iposGlobalPrimalVelocityY
    INTEGER(PREC_VECIDX) :: iposGlobalPrimalPressure
    INTEGER(PREC_VECIDX) :: iposGlobalDualVelocityX,iposGlobalDualVelocityY
    INTEGER(PREC_VECIDX) :: iposGlobalDualPressure 
    INTEGER(PREC_ELEMENTIDX) :: iel,NEL
    INTEGER(PREC_VERTEXIDX) :: NVT
    INTEGER(PREC_DOFIDX) :: NEQ
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_DOFIDX), DIMENSION(ndofLocalVelocityX) :: IvelocityDofs
    INTEGER(PREC_DOFIDX), DIMENSION(ndofLocalPressure) :: IpressureDofs
    !INTEGER(PREC_MATIDX), DIMENSION(ndofLocalVelocity,ndofLocalVelocity) :: ImatPosVel
    !INTEGER(PREC_MATIDX), DIMENSION(ndofLocalPressure,ndofLocalVelocity) :: ImatPosPres
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolA,p_KcolB,p_KcolC
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldA,p_KldB,p_KldC
    REAL(DP), DIMENSION(ndofLocal) :: DxLocal, DbLocal
    
    INTEGER :: nblockSize,i,j,k,ichunk,NEQtime,ichunkpos,istep,ipos,iiteration
    INTEGER :: h_ImatrixMapping
    
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! Array of system matrices for the timesteps. 
    TYPE(t_matrixBlock), DIMENSION(:,:), ALLOCATABLE :: RsystemMatrix
    
    ! Arrays with pointers to the matrix data of the matrices in RsystemMatrix
    TYPE(t_matrixRow), DIMENSION(:), ALLOCATABLE :: RmatrixRows
    
    ! The same for the solution/rhs vector
    TYPE(t_chunkVector), DIMENSION(:), ALLOCATABLE :: Rvectors
    
    ! Global DOF's in one block of the current chunk
    INTEGER(PREC_VECIDX), DIMENSION(:), ALLOCATABLE :: IchunkDofs
    
    ! Matrix position structure.
    ! Here, the positions of the matrix elements are saved, corresponding to
    ! the current element.
    TYPE(t_matrixPositions) :: rmatPositions
    
    ! Array of vectors simultaneously to handle for all the timesteps
    TYPE(t_vectorBlock), DIMENSION(:), ALLOCATABLE :: RrhsGlobal,RsolutionGlobal
    TYPE(t_vectorBlock) :: rvectorTemp
    
    ! local matrix,vectors
    TYPE(t_matrixScalar) :: rmatrixLocal
    REAL(DP), DIMENSION(:), POINTER :: p_DaLocal
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KcolLocal
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KldLocal
    TYPE(t_vectorScalar) :: rrhsLocal
    REAL(DP), DIMENSION(:), POINTER :: p_DrhsLocal
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_ImatrixMapping
    TYPE(t_matrixBlock), DIMENSION(1) :: RmatrixLocalBlock
    TYPE(t_vectorBlock) :: rrhsBlock
    TYPE(t_linsolNode), POINTER :: p_rsolver
    INTEGER :: ierror
    
    TYPE(t_vectorBlock) :: rdbgrhs,rdbgsol
    TYPE(t_matrixBlock) :: rdbgMatrix
    TYPE(t_linsolNode), POINTER :: p_rdbgsolver
    TYPE(t_matrixBlock), DIMENSION(1) :: RdbgmatrixLocalBlock
    
    ! Fetch some information
    p_rdiscretisation => rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%rdiscretisation
    NEL = p_rdiscretisation%p_rtriangulation%NEL
    NEQ = dof_igetNDofGlobblock(p_rdiscretisation)
    NVT = p_rdiscretisation%p_rtriangulation%NVT
    CALL storage_getbase_int2d (p_rdiscretisation%p_rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    
    ! Calculate the number of block that must be handled simultaneously.
    NEQtime = rmatrix%p_rspaceTimeDiscretisation%NEQtime
    nblockSize = MIN(NEQtime,MAX(1,nchunkSize+1))

    ! Allocate memory for global matrices.
    !
    ! System matrices on the diagonal of the global matrix.
    ! RsystemMatrix(-1) describes the matrices coupling the primal velocity
    ! and roughly contains mass matrices.
    ! RsystemMatrix(0) contains the system matrices on the main diagonal
    ! of the global matrix.
    ! RsystemMatrix(1) describes the matrices coupling the dual velocity
    ! and roughly contains mass matrices.
    ALLOCATE(RsystemMatrix(-1:1,nblockSize))
    
    ! Allocate one additional element before and after the next arrays. These elements
    ! are unused but there existence simplifies the handling of the array
    ! at the beginning and the ending of a chunk.
    ALLOCATE(RrhsGlobal(0:nblockSize+1))
    ALLOCATE(RsolutionGlobal(0:nblockSize+1))
    ALLOCATE(Rvectors(0:nblockSize+1))
    
    ! RmatrixRows contains pointers to the matrix data arrays of the matrices in
    ! RsystemMatrix. That gives some speed, as we don't have to access the
    ! storage-module so often.
    ALLOCATE(RmatrixRows(nblockSize))
    
    ! IchunkDofs holds the numbers of the global DOF's in one block of the current chunk.
    ALLOCATE(IchunkDofs(1:ndofLocal))
    
    ! Allocate all the matrices/vectors we have to handle simultaneously.
    DO j=1,nblockSize
      DO i=-1,1
        CALL cc_allocSystemMatrix (rproblem,rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo,&
          RsystemMatrix(i,j))
      END DO
      
      CALL lsysbl_createVecBlockByDiscr (p_rdiscretisation,RrhsGlobal(j))
      CALL lsysbl_createVecBlockByDiscr (p_rdiscretisation,RsolutionGlobal(j))
      
      ! Get pointers to the vectors for quicker access
      CALL lsysbl_getbase_double (RrhsGlobal(j),Rvectors(j)%p_Db)
      CALL lsysbl_getbase_double (RsolutionGlobal(j),Rvectors(j)%p_Dx)

      CALL lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(1),Rvectors(j)%p_Dx1)
      CALL lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(2),Rvectors(j)%p_Dx2)
      CALL lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(3),Rvectors(j)%p_Dx3)
      CALL lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(4),Rvectors(j)%p_Dx4)
      CALL lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(5),Rvectors(j)%p_Dx5)
      CALL lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(6),Rvectors(j)%p_Dx6)
    END DO
    CALL lsysbl_createVecBlockByDiscr (p_rdiscretisation,rvectorTemp)
    
    ! Allocate memory for the matrix positions.
    ! These change with every element and allow us to quickly access for every
    ! entry in each vector the corresponding matrix entry.
    ALLOCATE(rmatPositions%p_IposA(ndofLocalVelocityX,ndofLocalVelocityY))
    ALLOCATE(rmatPositions%p_IposB(ndofLocalVelocityX,ndofLocalPressure))
    ALLOCATE(rmatPositions%p_IposC(ndofLocalPressure,ndofLocalPressure))
    
    ! Allocate memory for local matrices
    CALL createLocalMatrixQ1TQ0 (&
      p_rdiscretisation,nblockSize-1,rmatrixLocal,h_ImatrixMapping)
      
    CALL storage_getbase_int2d (h_ImatrixMapping,p_ImatrixMapping)
      
    ! Get pointers to the matrix data.
    CALL lsyssc_getbase_double (rmatrixLocal,p_DaLocal)
    CALL lsyssc_getbase_Kcol (rmatrixLocal,p_KcolLocal)
    CALL lsyssc_getbase_Kld (rmatrixLocal,p_KldLocal)
      
    ! Create appropriate RHS and solution vectors.
    CALL lsyssc_createVecIndMat (rmatrixLocal,rrhsLocal,.FALSE.)
    
    ! Get pointers to the data
    CALL lsyssc_getbase_double (rrhsLocal,p_DrhsLocal)
    
    ! Create derived 1x1 block matrices/vectors for the solver.
    CALL lsysbl_createMatFromScalar (rmatrixLocal,RmatrixLocalBlock(1))
    CALL lsysbl_createVecFromScalar (rrhsLocal,rrhsBlock)
    
    ! With the matrix, prepare an UMFPACK solver to solve our system
    CALL linsol_initUMFPACK4(p_rsolver)
    CALL linsol_setMatrices(p_rsolver,RmatrixLocalBlock)
    
    CALL lsyssc_clearMatrix (rmatrixLocal,1.0_DP)
    CALL linsol_initStructure (p_rsolver,ierror)
    IF (ierror .NE. 0) THEN
      PRINT *,'UMFPACK for local system cannot be initialised!'
      CALL sys_halt()
    END IF
    
    ! We assume all matrices to share the same structure! (I.e. all velocity
    ! and mass matrices on the one hand, all B and B^T matrices on the other hand.)
    ! Fetch the matrix structure arrays.
    CALL lsyssc_getbase_Kcol (RsystemMatrix(0,1)%RmatrixBlock(1,1),p_KcolA)
    CALL lsyssc_getbase_Kld (RsystemMatrix(0,1)%RmatrixBlock(1,1),p_KldA)
    CALL lsyssc_getbase_Kcol (RsystemMatrix(0,1)%RmatrixBlock(1,3),p_KcolB)
    CALL lsyssc_getbase_Kld (RsystemMatrix(0,1)%RmatrixBlock(1,3),p_KldB)
    CALL lsyssc_getbase_Kcol (RsystemMatrix(0,1)%RmatrixBlock(3,3),p_KcolC)
    CALL lsyssc_getbase_Kld (RsystemMatrix(0,1)%RmatrixBlock(3,3),p_KldC)

    ! Set up index positions to quickly find entries in a spatial vector.
    iposGlobalPrimalVelocityX = 1
    iposGlobalPrimalVelocityY = iposGlobalPrimalVelocityX + &
        dof_igetNDofGlob(p_rdiscretisation%RspatialDiscretisation(1))
    iposGlobalPrimalPressure = iposGlobalPrimalVelocityY + &
        dof_igetNDofGlob(p_rdiscretisation%RspatialDiscretisation(2))
    iposGlobalDualVelocityX = iposGlobalPrimalPressure + &
        dof_igetNDofGlob(p_rdiscretisation%RspatialDiscretisation(3))
    iposGlobalDualVelocityY = iposGlobalDualVelocityX + &
        dof_igetNDofGlob(p_rdiscretisation%RspatialDiscretisation(4))
    iposGlobalDualPressure = iposGlobalDualVelocityY + &
        dof_igetNDofGlob(p_rdiscretisation%RspatialDiscretisation(5))

    ! Perforn niterations iterations
    DO iiteration = 1,niterations

      ! Loop through all chunks. Every chunk consists of nchunkSize timesteps -- 
      ! rounded up; or more precisely, every chunk consists of nblocks vectors.
      DO ichunk = 1,(NEQtime+nblockSize-1)/nblockSize
      
        ! Calculate the chunk position. This is normally ichunk*nblockSize
        ! except if this doesn't fit to the end position. In the latter case,
        ! we modify the chunk position appropriately such that we have
        ! nblockSize blocks in the chunk.
        ichunkPos = (ichunk-1)*nblockSize
        IF ((ichunkPos+nblockSize) .GT. NEQtime) &
          ichunkPos = NEQtime-nblockSize
        ichunkPos = ichunkPos + 1
        
        ! For all items in the current chunk, get the global solution and 
        ! RHS vectors.
        DO ipos = 1,nblockSize
          CALL sptivec_getTimestepData (rx, ichunkPos+ipos-1, RsolutionGlobal(ipos))
          CALL sptivec_getTimestepData (rb, ichunkPos+ipos-1, RrhsGlobal(ipos))
        END DO
        
        ! Build the matrices that belong to our chunk
        CALL buildMatrices (rproblem,rmatrix,&
            ichunkPos,nblockSize,RsystemMatrix,RmatrixRows)
        
        ! Set up "b-Mx" for the submatrix 'before' the chunk and 'after' the chunk;
        ! this gives the time coupling to the previous and next chunk.
        IF (ichunkPos .GT. 1) THEN
          CALL sptivec_getTimestepData (rx, ichunkPos-1, rvectorTemp)
          !CALL matio_writeBlockMatrixHR (RsystemMatrix(-1,1), 'matleft',&
          !                          .TRUE., 0, 'matleft.txt', '(E10.2)')
          CALL lsysbl_blockMatVec (RsystemMatrix(-1,1),&
              rvectorTemp,RrhsGlobal(1),-1.0_DP,1.0_DP)
        END IF
        
        IF (ichunkPos+nblockSize-1 .LT. NEQtime) THEN
          CALL sptivec_getTimestepData (rx, ichunkPos+nblockSize-1+1, rvectorTemp)
          !CALL matio_writeBlockMatrixHR (RsystemMatrix(1,nblockSize), 'matright',&
          !                          .TRUE., 0, 'matright.txt', '(E10.2)')
          CALL lsysbl_blockMatVec (RsystemMatrix(1,nblockSize),&
              rvectorTemp,RrhsGlobal(nblockSize),-1.0_DP,1.0_DP)
        END IF
        
        
        ! DEBUG!!!
        ! globalen Defekt bilden, UMFPACK anwerfen und Lösung draufaddieren.
        ! Sollte je nach Einstellung einem globalen UMFPACK oder einem
        ! UMFPACK in jd. Zeitschritt entsprechen.
        !...
        CALL lsysbl_createEmptyMatrix (rdbgmatrix,6*nblockSize)
        DO j=1,nblockSize
          DO i=1,3
            IF (j .GT. 1) THEN
              CALL lsysbl_insertSubmatrix (&
                  RsystemMatrix(-1,j),rdbgmatrix,&
                  LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,6*(j-1)+1,6*(j-1)+1-1)
            END IF
            CALL lsysbl_insertSubmatrix (&
                RsystemMatrix(0,j),rdbgmatrix,&
                LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,6*(j-1)+1,6*(j-1)+1)
            IF (j .LT. nblockSize) THEN
              CALL lsysbl_insertSubmatrix (&
                  RsystemMatrix(1,j),rdbgmatrix,&
                  LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,6*(j-1)+1,6*(j-1)+1+1)
            END IF
          END DO
        END DO
        CALL lsysbl_updateMatStrucInfo(rdbgmatrix)
        CALL matio_writeBlockMatrixHR (rdbgmatrix, 'fullmat',&
                                       .TRUE., 0, 'matfull.txt', '(E10.2)')
                                       
        CALL linsol_initUmfpack4 (p_rdbgSolver)
        CALL lsysbl_duplicateMatrix (rdbgmatrix,RdbgmatrixLocalBlock(1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        CALL linsol_setMatrices (p_rdbgSolver,RdbgmatrixLocalBlock)
        CALL linsol_initStructure (p_rdbgSolver,ierror)
        IF (ierror .NE. 0) THEN
          PRINT *,'UMFPACK for local system cannot be initialised!'
          CALL sys_halt()
        END IF
        CALL linsol_initData (p_rdbgSolver,ierror)
        IF (ierror .NE. 0) THEN
          PRINT *,'UMFPACK for local system cannot be initialised!'
          CALL sys_halt()
        END IF
        
        CALL lsysbl_createVecBlockIndMat(rdbgmatrix,rdbgrhs,.TRUE.)
        CALL lsysbl_createVecBlockIndMat(rdbgmatrix,rdbgsol,.TRUE.)
        DO i=1,nblockSize
          DO j=1,6
            CALL lsyssc_copyVector (RrhsGlobal(i)%RvectorBlock(j),&
                rdbgrhs%RvectorBlock((i-1)*6+j))
            CALL lsyssc_copyVector (RsolutionGlobal(i)%RvectorBlock(j),&
                rdbgsol%RvectorBlock((i-1)*6+j))
          END DO
        END DO
        CALL lsysbl_blockMatVec (rdbgmatrix,rdbgsol,rdbgrhs,-1.0_DP,1.0_DP)
        
        CALL linsol_precondDefect (p_rdbgSolver,rdbgrhs)
        CALL linsol_releaseSolver(p_rdbgSolver)

        DO i=1,nblockSize
          DO j=1,6
            CALL lsyssc_vectorLinearComb (rdbgrhs%RvectorBlock((i-1)*6+j),&
                RsolutionGlobal(i)%RvectorBlock(j),0.9_DP,1.0_DP)
          END DO
        END DO
                                       
        CALL lsysbl_releaseVector (rdbgsol)
        CALL lsysbl_releaseVector (rdbgrhs)
        CALL lsysbl_releaseMatrix (rdbgmatrix)
        CALL lsysbl_releaseMatrix (RdbgmatrixLocalBlock(1))
        
        
        
        
        
!        ! Now our chunk is a global problem that looks as follows (here an example
!        ! for chunk size 3, timestep 1):
!        !
!        !   A B M  |       |                   
!        !   B C    |       |          x1       f1-Mx0
!        !   M   A B|       |                 
!        !       B C|    M  |                 
!        !   -------+-------+-------          
!        !   M      |A B M  |                 
!        !          |B C    |          x2   =   f2
!        !          |M   A B|    M            
!        !          |    B C|                 
!        !   -------+-------+-------          
!        !          |M      |A B M            
!        !          |       |B C       x3       f3-Mx4
!        !          |       |M   A B          
!        !          |       |    B C            
!        !
!        !          
!        
!        ! What we have to do now is a loop over all spatial elements.
!        ! For every spatial element, we'll simultaneously handle all
!        ! timesteps in the chunk.
!        DO iel = 1,NEL
!        
!          ! Determine the spatial DOF's on our current element that we have
!          ! to extract from the matrices.
!          ! This is more or less a replacement for the DOFMAPPING to get some
!          ! speedup, with the disadvantage that this restricts to the Q1~/Q0
!          ! discretisation!
!          !
!          ! Get the pressure-DOF -- which is the number of the element.
!          IpressureDofs(1) = iel
!          
!          ! Get the velocity-DOF's. These are the numbers of the edges.
!          ! The local DOF's on the current element are 1,2,3,4. The 
!          ! corresponding global DOF's are...
!          IvelocityDofs(1) = p_IedgesAtElement(1,iel)-NVT
!          IvelocityDofs(2) = p_IedgesAtElement(2,iel)-NVT
!          IvelocityDofs(3) = p_IedgesAtElement(3,iel)-NVT
!          IvelocityDofs(4) = p_IedgesAtElement(4,iel)-NVT
!          
!          ! Figure out the numbers of the global DOF's in the current chunk.
!          ! IchunkDofs therefore defines a mapping of the 'local' DOF's 1..18
!          ! of one block in the chunk to the 'global' DOF's in each row of the chunk.
!          IchunkDofs(1) = iposGlobalPrimalVelocityX - 1 + IvelocityDofs(1)
!          IchunkDofs(2) = iposGlobalPrimalVelocityX - 1 + IvelocityDofs(2)
!          IchunkDofs(3) = iposGlobalPrimalVelocityX - 1 + IvelocityDofs(3)
!          IchunkDofs(4) = iposGlobalPrimalVelocityX - 1 + IvelocityDofs(4)
!          
!          IchunkDofs(5) = iposGlobalPrimalVelocityY - 1 + IvelocityDofs(1)
!          IchunkDofs(6) = iposGlobalPrimalVelocityY - 1 + IvelocityDofs(2)
!          IchunkDofs(7) = iposGlobalPrimalVelocityY - 1 + IvelocityDofs(3)
!          IchunkDofs(8) = iposGlobalPrimalVelocityY - 1 + IvelocityDofs(4)
!          
!          IchunkDofs(9) = iposGlobalPrimalPressure  - 1 + IpressureDofs(1)
!          
!          IchunkDofs(10) = iposGlobalDualVelocityX - 1 + IvelocityDofs(1)
!          IchunkDofs(11) = iposGlobalDualVelocityX - 1 + IvelocityDofs(2)
!          IchunkDofs(12) = iposGlobalDualVelocityX - 1 + IvelocityDofs(3)
!          IchunkDofs(13) = iposGlobalDualVelocityX - 1 + IvelocityDofs(4)
!          
!          IchunkDofs(14) = iposGlobalDualVelocityY - 1 + IvelocityDofs(1)
!          IchunkDofs(15) = iposGlobalDualVelocityY - 1 + IvelocityDofs(2)
!          IchunkDofs(16) = iposGlobalDualVelocityY - 1 + IvelocityDofs(3)
!          IchunkDofs(17) = iposGlobalDualVelocityY - 1 + IvelocityDofs(4)
!          
!          IchunkDofs(18) = iposGlobalDualPressure  - 1 + IpressureDofs(1)
!          
!          ! Figure out the matrix positions that belong to these DOF's; in
!          ! the mass/velocity matrices as well as in the B-matrices.
!          CALL getMatrixPositions (IvelocityDofs,IpressureDofs,rmatPositions,&
!              p_KcolA,p_KldA,p_KcolB,p_KldB,p_KcolC,p_KldC)
!          !CALL calcMatrixIndices (p_KcolA,p_KldA,IvelocityDofs,IvelocityDofs,ImatPosVel)
!          !CALL calcMatrixIndices (p_KcolB,p_KldB,IvelocityDofs,IpressureDofs,ImatPosPres)
!          
!          ! Set up the local matrix. This is used for preconditioning later.
!          CALL getLocalMatrices (Rmatrixrows,IvelocityDofs,IpressureDofs,&
!              p_ImatrixMapping,rmatPositions,p_DaLocal,p_KcolLocal,p_KldLocal,&
!              p_KcolA,p_KldA,p_KcolB,p_KldB,p_KcolC,p_KldC)
!              
!          CALL matio_writeMatrixHR (rmatrixLocal, 'localmat',&
!                                    .TRUE., 0, 'lmatrix.txt', '(E10.2)')
!                                    
!          ! On the current element, treat all timesteps of the current chunk.
!          ! Each timestep corresponds to one 'row' in the local matrix.
!          DO istep = 1,nblockSize
!            ! Calculate the local defect in this block of the chunk
!            CALL getLocalDefect (p_DrhsLocal((istep-1)*SIZE(IchunkDofs)+1:(istep-1+1)*SIZE(IchunkDofs)),&
!                Rmatrixrows(istep),Rvectors(istep-1:istep+1),&
!                IvelocityDofs,IpressureDofs,IchunkDofs,&
!                rmatPositions,istep .EQ. 1,istep .EQ. nblockSize,&
!                p_KcolA,p_KldA,p_KcolB,p_KldB,p_KcolC,p_KldC)
!          END DO
!        
!          ! Initialise the solver with that matrix
!          CALL linsol_initData (p_rsolver,ierror)
!          IF (ierror .NE. 0) THEN
!            PRINT *,'UMFPACK for local system cannot be initialised!'
!            CALL sys_halt()
!          END IF
!          
!          ! Solve; more precisely, precondition our local defect
!          CALL linsol_precondDefect (p_rsolver,rrhsBlock)
!          
!          ! Release solver data
!          CALL linsol_doneData (p_rsolver)
!          
!          ! Add the local correction vector to our solution
!          DO istep = 1,nblockSize
!            CALL addLocalCorrection (p_DrhsLocal((istep-1)*SIZE(IchunkDofs)+1:(istep-1+1)*SIZE(IchunkDofs)),&
!                Rvectors(istep),IchunkDofs,domega)
!          END DO
!        
!        END DO
      
        ! Save the new solution vectors.
        DO ipos = 1,nblockSize
          CALL sptivec_setTimestepData (rx, ichunkPos+ipos-1, RsolutionGlobal(ipos))
        END DO

      END DO
      
    END DO

    ! Release the solver
    CALL linsol_releaseSolver (p_rsolver)

    ! Release memory
    DEALLOCATE(rmatPositions%p_IposA)
    DEALLOCATE(rmatPositions%p_IposB)
    DEALLOCATE(rmatPositions%p_IposC)
    
    CALL lsysbl_releaseVector (rvectorTemp)
    CALL lsysbl_releaseVector (rrhsBlock)
    CALL lsysbl_releaseMatrix (RmatrixLocalBlock(1))

    CALL lsyssc_releaseVector (rrhsLocal)

    CALL storage_free(h_ImatrixMapping)
    CALL lsyssc_releaseMatrix (rmatrixLocal)
    
    DEALLOCATE(IchunkDofs)
    
    DO j=1,nblockSize
      CALL lsysbl_releaseVector (RsolutionGlobal(j))
      CALL lsysbl_releaseVector (RrhsGlobal(j))
      DO i=-1,1
        CALL lsysbl_releaseMatrix(RsystemMatrix(i,j))
      END DO
    END DO
    DEALLOCATE(RmatrixRows)
    DEALLOCATE(Rvectors)
    DEALLOCATE(RsystemMatrix)
    DEALLOCATE(RsolutionGlobal)
    DEALLOCATE(RrhsGlobal)
    
  CONTAINS
  
    SUBROUTINE buildMatrices (rproblem,rspaceTimeMatrix,&
      ichunkPos,nblockSize,RsystemMatrix,RmatrixRows)  
    
    ! Builds the matrices in timestep ichunkPos..ichunkPos+nblockSize-1 into
    ! the RsystemMatrix and RmatrixRows arrays.
    
    ! The main problem structure
    TYPE(t_problem), INTENT(INOUT) :: rproblem
    
    ! The underlying space time matrix.
    TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN) :: rspaceTimeMatrix
    
    ! Chunk position
    INTEGER, INTENT(IN) :: ichunkPos
    
    ! Size of a chunk
    INTEGER, INTENT(IN) :: nblockSize
    
    ! Destination array for the matrices
    TYPE(t_matrixBlock), DIMENSION(-1:,:), INTENT(INOUT) :: RsystemMatrix
    
    ! Destination array for matrix pointers.
    ! Variable is declared with INTENT(OUT), so all pointers are
    ! nullified automatically.
    TYPE(t_matrixRow), DIMENSION(:), INTENT(OUT) :: RmatrixRows
    
      ! local variables
      INTEGER :: i,j,k,ichunk,NEQtime,ichunkrel
      REAL(DP) :: dtheta
      TYPE(t_ccmatrixComponents) :: rmatrixComponents
      TYPE(t_vectorBlock), DIMENSION(3) :: rtimeVector
      TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr

      ! Create temp vectors for evaluating the nonlinearity      
      CALL lsysbl_createVecBlockIndMat(RsystemMatrix(0,1),rtimeVector(1))
      CALL lsysbl_createVecBlockIndMat(RsystemMatrix(0,1),rtimeVector(2))
      CALL lsysbl_createVecBlockIndMat(RsystemMatrix(0,1),rtimeVector(3))
      
      ! Fetch some information
      NEQtime = rmatrix%p_rspaceTimeDiscretisation%NEQtime
      dtheta = rproblem%rtimedependence%dtimeStepTheta
      
      ! Basic initialisation of rmatrixComponents with the pointers to the
      ! matrices / discretisation structures on the current level.
      !
      ! The weights in the rmatrixComponents structure are later initialised
      ! according to the actual situation when the matrix is to be used.
      p_rspaceTimeDiscr => rspaceTimeMatrix%p_rspaceTimeDiscretisation
      
      rmatrixComponents%p_rdiscretisation         => &
          p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation
      rmatrixComponents%p_rmatrixStokes           => &
          p_rspaceTimeDiscr%p_rlevelInfo%rmatrixStokes          
      rmatrixComponents%p_rmatrixB1               => &
          p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB1              
      rmatrixComponents%p_rmatrixB2               => &
          p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB2              
      rmatrixComponents%p_rmatrixMass             => &
          p_rspaceTimeDiscr%p_rlevelInfo%rmatrixMass            
      rmatrixComponents%p_rmatrixIdentityPressure => &
          p_rspaceTimeDiscr%p_rlevelInfo%rmatrixIdentityPressure

      rmatrixComponents%dnu = collct_getvalue_real (rproblem%rcollection,'NU')
      rmatrixComponents%iupwind1 = collct_getvalue_int (rproblem%rcollection,'IUPWIND1')
      rmatrixComponents%dupsam1 = collct_getvalue_real (rproblem%rcollection,'UPSAM1')
      rmatrixComponents%iupwind2 = collct_getvalue_int (rproblem%rcollection,'IUPWIND2')
      rmatrixComponents%dupsam2 = collct_getvalue_real (rproblem%rcollection,'UPSAM2')

      ! Fetch evaluation vectors for the nonlinearity
      IF (ichunkPos .GT. 1) &
        CALL sptivec_getTimestepData (rspaceTimeMatrix%p_rsolution, ichunkPos-1, rtimeVector(1))
        
      CALL sptivec_getTimestepData (rspaceTimeMatrix%p_rsolution, ichunkPos, rtimeVector(2))
      
      ! When retrieving matrix pointers, don't check the scaling factor!
      ! We need the 'maximum possible' layout of the matrix and so we get
      ! all pointers we are able to get.
      
      DO ichunk = ichunkPos, ichunkPos+nblockSize-1
      
        ! Current point in time
        rproblem%rtimedependence%dtime = &
            rproblem%rtimedependence%dtimeInit + (ichunk-1)*p_rspaceTimeDiscr%rtimeDiscr%dtstep
        rproblem%rtimedependence%itimestep = ichunk-1

        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
          CALL cc_updateDiscreteBC (rproblem)
        END IF

        ! Relative position of the chunk item
        ichunkrel = ichunk-ichunkPos+1
      
        ! Get the evaluation point of the nonlinearity in the next timestep
        IF (ichunk .LT. NEQtime) THEN
          CALL sptivec_getTimestepData (rspaceTimeMatrix%p_rsolution, ichunk+1, rtimeVector(3))
        END IF
        
        ! Build 'left' matrix (if it exists):
        IF (ichunk .GT. 1) THEN
        
          ! Set up the matrix weights 
          CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
            ichunk-1,-1,rmatrixComponents)
          
          ! Set up the matrix
          CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
              RsystemMatrix(-1,ichunkrel),rmatrixComponents,rtimeVector(1),rtimeVector(2),rtimeVector(3)) 
        
          ! Include the boundary conditions into that matrix.
          RsystemMatrix(-1,ichunkrel)%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
          CALL matfil_discreteBC (RsystemMatrix(-1,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC)
          CALL matfil_discreteFBC (RsystemMatrix(-1,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC)

          ! Get the pointers to the submatrices.
          IF (lsysbl_isSubmatrixPresent (RsystemMatrix(-1,ichunkrel),1,1)) &
            CALL lsyssc_getbase_double (RsystemMatrix(-1,ichunkrel)%RmatrixBlock(1,1),&
                RmatrixRows(ichunkrel)%p_Dp11)
          
          IF (lsysbl_isSubmatrixPresent (RsystemMatrix(-1,ichunkrel),2,1)) &
            CALL lsyssc_getbase_double (RsystemMatrix(-1,ichunkrel)%RmatrixBlock(2,1),&
                RmatrixRows(ichunkrel)%p_Dp21)
          
          IF (lsysbl_isSubmatrixPresent (RsystemMatrix(-1,ichunkrel),1,2)) &
            CALL lsyssc_getbase_double (RsystemMatrix(-1,ichunkrel)%RmatrixBlock(1,2),&
                RmatrixRows(ichunkrel)%p_Dp12)
          
          IF (lsysbl_isSubmatrixPresent (RsystemMatrix(-1,ichunkrel),2,2)) &
            CALL lsyssc_getbase_double (RsystemMatrix(-1,ichunkrel)%RmatrixBlock(2,2),&
                RmatrixRows(ichunkrel)%p_Dp22)
        
        END IF
        
        ! Build the 'diagonal' matrix
        !
        ! Set up the matrix weights 
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ichunk-1,0,rmatrixComponents)
        
        ! Set up the matrix
        CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            RsystemMatrix(0,ichunkrel),rmatrixComponents,&
            rtimeVector(1),rtimeVector(2),rtimeVector(3)) 
        
        ! Include the boundary conditions into that matrix.
        CALL matfil_discreteBC (RsystemMatrix(0,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC)
        CALL matfil_discreteFBC (RsystemMatrix(0,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC)

        ! Get matrix pointers
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),1,1)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(1,1),&
              RmatrixRows(ichunkrel)%p_Da11)
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),2,1)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(2,1),&
              RmatrixRows(ichunkrel)%p_Da21)

        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),3,1)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(3,1),&
              RmatrixRows(ichunkrel)%p_Da31)

        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),4,1)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(4,1),&
              RmatrixRows(ichunkrel)%p_Da41)

        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),5,1)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(5,1),&
              RmatrixRows(ichunkrel)%p_Da51)
              
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),1,2)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(1,2),&
              RmatrixRows(ichunkrel)%p_Da12)
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),2,2)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(2,2),&
              RmatrixRows(ichunkrel)%p_Da22)

        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),3,2)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(3,2),&
              RmatrixRows(ichunkrel)%p_Da32)

        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),4,2)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(4,2),&
              RmatrixRows(ichunkrel)%p_Da42)

        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),5,2)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(5,2),&
              RmatrixRows(ichunkrel)%p_Da52)
        
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),1,3)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(1,3),&
              RmatrixRows(ichunkrel)%p_Da13)
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),2,3)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(2,3),&
              RmatrixRows(ichunkrel)%p_Da23)

        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),3,3)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(3,3),&
              RmatrixRows(ichunkrel)%p_Da33)
        
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),1,4)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(1,4),&
              RmatrixRows(ichunkrel)%p_Da14)
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),2,4)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(2,4),&
              RmatrixRows(ichunkrel)%p_Da24)
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),4,4)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(4,4),&
              RmatrixRows(ichunkrel)%p_Da44)
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),5,4)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(5,4),&
              RmatrixRows(ichunkrel)%p_Da54)

        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),6,4)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(6,4),&
              RmatrixRows(ichunkrel)%p_Da64)


        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),1,5)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(1,5),&
              RmatrixRows(ichunkrel)%p_Da15)
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),2,5)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(2,5),&
              RmatrixRows(ichunkrel)%p_Da25)
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),4,5)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(4,5),&
              RmatrixRows(ichunkrel)%p_Da45)
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),5,5)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(5,5),&
              RmatrixRows(ichunkrel)%p_Da55)

        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),6,5)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(6,5),&
              RmatrixRows(ichunkrel)%p_Da65)


        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),4,6)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(4,6),&
              RmatrixRows(ichunkrel)%p_Da46)
        
        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),5,6)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(5,6),&
              RmatrixRows(ichunkrel)%p_Da56)

        IF (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),6,6)) &
          CALL lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(6,6),&
              RmatrixRows(ichunkrel)%p_Da66)
        
        ! Build 'right' matrix (if it exists):
        IF (ichunk .LT. NEQtime) THEN
        
          ! Set up the matrix weights 
          CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
            ichunk-1,1,rmatrixComponents)
          
          ! Set up the matrix
          CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
              RsystemMatrix(1,ichunkrel),rmatrixComponents,&
              rtimeVector(1),rtimeVector(2),rtimeVector(3)) 
        
          ! Include the boundary conditions into that matrix.
          RsystemMatrix(1,ichunkrel)%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
          CALL matfil_discreteBC (RsystemMatrix(1,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC)
          CALL matfil_discreteFBC (RsystemMatrix(1,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC)

          ! Get the pointers to the submatrices
          IF (lsysbl_isSubmatrixPresent (RsystemMatrix(1,ichunkrel),4,4)) &
            CALL lsyssc_getbase_double (RsystemMatrix(1,ichunkrel)%RmatrixBlock(4,4),&
                RmatrixRows(ichunkrel)%p_Dd44)
          
          IF (lsysbl_isSubmatrixPresent (RsystemMatrix(1,ichunkrel),5,4)) &
            CALL lsyssc_getbase_double (RsystemMatrix(1,ichunkrel)%RmatrixBlock(5,4),&
                RmatrixRows(ichunkrel)%p_Dd54)
          
          IF (lsysbl_isSubmatrixPresent (RsystemMatrix(1,ichunkrel),4,5)) &
            CALL lsyssc_getbase_double (RsystemMatrix(1,ichunkrel)%RmatrixBlock(4,5),&
                RmatrixRows(ichunkrel)%p_Dd45)
          
          IF (lsysbl_isSubmatrixPresent (RsystemMatrix(1,ichunkrel),5,5)) &
            CALL lsyssc_getbase_double (RsystemMatrix(1,ichunkrel)%RmatrixBlock(5,5),&
                RmatrixRows(ichunkrel)%p_Dd55)
        END IF
        

        IF (ichunk .LT. NEQtime) THEN
          ! Cycle evaluation vectors: 1 <- 2 <- 3
          CALL lsysbl_copyVector (rtimeVector(2),rtimeVector(1))
          CALL lsysbl_copyVector (rtimeVector(3),rtimeVector(2))
        END IF
        
        ! Copy all scaling factors
        DO k=-1,1
          DO j=1,6
            DO i=1,6
              RmatrixRows(ichunkrel)%DscaleFactor(i,j,k) = &
                  RsystemMatrix(k,ichunkrel)%RmatrixBlock(i,j)%dscaleFactor
            END DO
          END DO
        END DO
      
      END DO
      
      ! Release memory
      CALL lsysbl_releaseVector (rtimeVector(3))
      CALL lsysbl_releaseVector (rtimeVector(2))
      CALL lsysbl_releaseVector (rtimeVector(1))
    
    END SUBROUTINE
    
    ! -----
    
    SUBROUTINE getMatrixPositions (IvelocityDofs,IpressureDofs,rmatrixPositions,&
        KcolA,KldA,KcolB,KldB,KcolC,KldC)
        
    ! Calculates the positions in the A, B and C matrices corresponding to the
    ! velocity and pressure DOF's.
    
    ! Array with velocity DOF's.
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IvelocityDofs
    
    ! Array with pressure DOF's.
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IpressureDofs
    
    ! A t_matrixPositions structure specifying the positions that
    ! are affected by the DOF's in the global matrices
    TYPE(t_matrixPositions), INTENT(IN) :: rmatrixPositions

    ! Structure of FEM matrices for the velocity
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA

    ! Structure of FEM matrices for the gradient
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldB
    
    ! Structure of FEM matrices for the pressure. Optional.
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: KcolC
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: KldC
    
      ! local variables
      INTEGER :: i,j
      INTEGER(PREC_MATIDX) :: k
      INTEGER(PREC_DOFIDX) :: idof
    
      ! Search for the matrix positions in the velocity and pressure matrix.
      ! IvelocityDofs specify the rows in the A- and B-matrix where we have
      ! to search.
      DO i=1,SIZE(IvelocityDofs)
        idof = IvelocityDofs(i)
        ! In the A-matrix, search for the velocity DOF's.
        DO j=1,SIZE(IvelocityDofs)
          DO k=KldA(idof),KldA(idof+1)-1
            IF (KcolA(k) .EQ. IvelocityDofs(j)) THEN
              rmatrixPositions%p_IposA(i,j) = k
              EXIT
            END IF
          END DO
        END DO
        
        ! In the B-matrices search for the pressure DOF's.
        DO j=1,SIZE(IpressureDofs)
          DO k=KldB(idof),KldB(idof+1)-1
            IF (KcolB(k) .EQ. IpressureDofs(j)) THEN
              rmatrixPositions%p_IposB(i,j) = k
              EXIT
            END IF
          END DO
        END DO
      END DO
    
      IF (PRESENT(KcolC)) THEN
        ! Search for the matrix positions pressure matrix.
        ! IpressureDofs specify the rows in the C-matrix where we have
        ! to search.
        DO i=1,SIZE(IpressureDofs)
          idof = IpressureDofs(i)
          ! In the C-matrix, search for the pressure DOF's.
          DO j=1,SIZE(IpressureDofs)
            DO k=KldC(idof),KldC(idof+1)-1
              IF (KcolC(k) .EQ. IpressureDofs(j)) THEN
                rmatrixPositions%p_IposC(i,j) = k
                EXIT
              END IF
            END DO
          END DO
        END DO
    
      END IF
    
    END SUBROUTINE
    
    ! -----
    
    SUBROUTINE getLocalMatrices (Rmatrixrows,IvelocityDofs,IpressureDofs,&
        ImatrixPositions,rmatrixPositions,Da,Kcol,Kld,KcolA,KldA,KcolB,KldB,KcolC,KldC)
    
    ! Extracts the entries of all matrices and creates a local matrix of
    ! a chunk.
    
    ! The matrices all the matrix rows in the chunk.
    TYPE(t_matrixRow), DIMENSION(:), INTENT(IN) :: Rmatrixrows
    
    ! The velocity DOF's under consideration
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IvelocityDofs
    
    ! The pressure DOF's under consideration
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IpressureDofs
    
    ! An integer array that specifies the start positions of the matrix columns
    ! in the local block matrix.
    INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: ImatrixPositions
    
    ! A t_matrixPositions structure specifying the positions that
    ! are affected by the DOF's in the global matrices
    TYPE(t_matrixPositions), INTENT(IN) :: rmatrixPositions
    
    ! Column/row structure of the mass/Laplace matrices
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA

    ! Column/row structure of the B matrices
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldB

    ! Column/row structure of the C matrices. Optional.
    ! If not present, there is no C-matrix.
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: KcolC
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: KldC
    
    ! The local matrix that is to be filled with data.
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Da
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld
    
      ! Local variables
      INTEGER :: irow,icolumn,ioffsetX,ioffsetY,irowOffsetY
    
      ! At first, initialise the output with zero.
      Da(:) = 0.0_DP
      
      ! Loop over all matrix rows in the chunk
      DO irow = 1,SIZE(Rmatrixrows)
      
        ! Now get all the submatrices and write them to the local matrix.
        ! What is a little bit tricky is the destination position where we
        ! have to write data to. Note that we know in advance where to
        ! write. The destination matrix is saved in CSR format without
        ! 0-blocks. By summing up all nonzero DOF's to a destination X-position
        ! we can directly calculate where we have to write data to.

        ! Calculate the Y-offsets
        irowOffsetY = (irow-1)*6
        ioffsetY = (irow-1)*ndofLocal
        
        ! In the first matrix row, ioffsetX is =0.
        ! In the other matrix rows, ioffsetX is set to 2*#local velocity DOF's
        ! to encounter the coupling matrix for the primal velocity.
        !
        ! ioffset is only used for the first and 2nd matrix row in the local
        ! block matrix.
        IF (irow .NE. 1) THEN
         
          ! Extract the mass matrices coupling the primal velocity.
          CALL getLocalMatrix (Rmatrixrows(irow)%p_Dp11,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(1,1,-1),&
              ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,1))
          
          CALL getLocalMatrix (Rmatrixrows(irow)%p_Dp21,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(2,1,-1),&
              ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,1))

          CALL getLocalMatrix (Rmatrixrows(irow)%p_Dp12,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(1,2,-1),&
              ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,2))
          
          CALL getLocalMatrix (Rmatrixrows(irow)%p_Dp22,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(2,2,-1),&
              ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,2))
          
        END IF
        
        ! Primal equation.
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da11,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(1,1,0),&
            ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,7))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da21,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(2,1,0),&
            ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,7))

        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da31,Da,Kld,IpressureDofs,&
            rmatrixPositions%p_IposB,.TRUE.,Rmatrixrows(irow)%DscaleFactor(3,1,0),&
            ioffsetY+iposLocalPrimalPressure,ImatrixPositions(irowOffsetY+3,7))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da41,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(4,1,0),&
            ioffsetY+iposLocalDualVelocityX,ImatrixPositions(irowOffsetY+4,7))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da51,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(5,1,0),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,7))


        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da12,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(1,2,0),&
            ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,8))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da22,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(2,2,0),&
            ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,8))

        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da32,Da,Kld,IpressureDofs,&
            rmatrixPositions%p_IposB,.TRUE.,Rmatrixrows(irow)%DscaleFactor(3,2,0),&
            ioffsetY+iposLocalPrimalPressure,ImatrixPositions(irowOffsetY+3,8))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da42,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(4,2,0),&
            ioffsetY+iposLocalDualVelocityX,ImatrixPositions(irowOffsetY+4,8))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da52,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(5,2,0),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,8))


        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da13,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposB,.FALSE.,Rmatrixrows(irow)%DscaleFactor(1,3,0),&
            ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,9))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da23,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposB,.FALSE.,Rmatrixrows(irow)%DscaleFactor(2,3,0),&
            ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,9))

        IF (PRESENT(KcolC)) &
          CALL getLocalMatrix (Rmatrixrows(irow)%p_Da33,Da,Kld,IpressureDofs,&
              rmatrixPositions%p_IposC,.FALSE.,Rmatrixrows(irow)%DscaleFactor(3,3,0),&
              ioffsetY+iposLocalPrimalPressure,ImatrixPositions(irowOffsetY+3,9))
            
            
        ! Dual equation
        !
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da14,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(1,4,0),&
            ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,10))
        
        ! Does not exist
        !CALL getLocalMatrix (Rmatrixrows(irow)%p_Da24,Da,Kld,IvelocityDofs,&
        !    rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(2,4,0),&
        !    ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,10))

        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da44,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(4,4,0),&
            ioffsetY+iposLocalDualVelocityX,ImatrixPositions(irowOffsetY+4,10))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da54,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(5,4,0),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,10))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da64,Da,Kld,IpressureDofs,&
            rmatrixPositions%p_IposB,.TRUE.,Rmatrixrows(irow)%DscaleFactor(6,4,0),&
            ioffsetY+iposLocalDualPressure,ImatrixPositions(irowOffsetY+6,10))
        

        ! Does not exist
        !CALL getLocalMatrix (Rmatrixrows(irow)%p_Da15,Da,Kld,IvelocityDofs,&
        !    rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(1,5,0),&
        !    ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,11))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da25,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(2,5,0),&
            ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,11))

        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da45,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(4,5,0),&
            ioffsetY+iposLocalDualVelocityX, ImatrixPositions(irowOffsetY+4,11))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da55,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(5,5,0),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,11))

        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da65,Da,Kld,IpressureDofs,&
            rmatrixPositions%p_IposB,.TRUE.,Rmatrixrows(irow)%DscaleFactor(6,5,0),&
            ioffsetY+iposLocalDualPressure,ImatrixPositions(irowOffsetY+6,11))


        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da46,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposB,.FALSE.,Rmatrixrows(irow)%DscaleFactor(4,6,0),&
            ioffsetY+iposLocalDualVelocityX,ImatrixPositions(irowOffsetY+4,12))
        
        CALL getLocalMatrix (Rmatrixrows(irow)%p_Da56,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposB,.FALSE.,Rmatrixrows(irow)%DscaleFactor(5,6,0),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,12))

        IF (PRESENT(KcolC)) &
          CALL getLocalMatrix (Rmatrixrows(irow)%p_Da66,Da,Kld,IpressureDofs,&
              rmatrixPositions%p_IposC,.FALSE.,Rmatrixrows(irow)%DscaleFactor(6,6,0),&
              ioffsetY+iposLocalDualPressure,ImatrixPositions(irowOffsetY+6,12))


        IF (irow .LT. SIZE(Rmatrixrows)) THEN
          ! Extract the mass matrices coupling the dual velocity.
          CALL getLocalMatrix (Rmatrixrows(irow)%p_Dd44,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(4,4,1),&
            ioffsetY+iposLocalDualVelocityX,ImatrixPositions(irowOffsetY+4,16))
          
          CALL getLocalMatrix (Rmatrixrows(irow)%p_Dd54,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(5,4,1),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,16))

          CALL getLocalMatrix (Rmatrixrows(irow)%p_Dd45,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(4,5,1),&
            ioffsetY+iposLocalDualVelocityX, ImatrixPositions(irowOffsetY+4,17))
          
          CALL getLocalMatrix (Rmatrixrows(irow)%p_Dd55,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(5,5,1),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,17))
        END IF

      END DO

    END SUBROUTINE

    ! -----

    SUBROUTINE getLocalMatrix (p_Da,Db,KldB,Idofs,&
        ImatrixPositions,btransposed,dscale,idestOffsetY,idestOffsetX)

    ! Pointer to the source matrix. May point to NULL() if there is none.
    REAL(DP), DIMENSION(:), POINTER :: p_Da

    ! Destination matrix
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Db
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldB
    
    ! Rows in the source matrix that should be extracted.
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Idofs
    
    ! Array with matrix positions. The rows in this array
    ! define for every entry in Idof the positions in the full matrix of
    ! the entries that have to be extracted and written to the local matrix.
    INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN) :: ImatrixPositions
    
    ! Scaling factor for the matrix
    REAL(DP), INTENT(IN) :: dscale
    
    ! Column offset in the destination matrix. Entries are written to
    ! position idestOffsetX..idestOffsetX+UBOUND(ImatrixPositions,2)-1 into
    ! line idestOffsetY..idestOffsetY+UBOUND(ImatrixPositions,1)-1 in
    ! the compressed destination matrix.
    ! That way, idestOffsetX/Y specifies the matrix column/row where to write
    ! the data to.
    INTEGER, INTENT(IN) :: idestOffsetX,idestOffsetY
    
    ! If set to TRUE, the routine assumes that ImatrixPosition specifies
    ! the entry positions in the transposed matrix. Used for handling B^T.
    LOGICAL, INTENT(IN) :: btransposed
    
      ! local variables
      INTEGER :: i,j
      INTEGER(PREC_MATIDX) :: k,irowPosB
      
      IF (.NOT. ASSOCIATED(p_Da)) RETURN
      
      IF (.NOT. btransposed) THEN
        ! Loop over the DOF's to extract
        DO i=1,UBOUND(ImatrixPositions,1)
          irowPosB = KldB(idestOffsetY-1+i)+idestOffsetX-1
          DO j=1,UBOUND(ImatrixPositions,2)
            
            ! Get the matrix position where to find our entry.
            k = ImatrixPositions(i,j)
            
            ! Get the matrix entry and write it to the destination matrix.
            Db(irowPosB+j-1) = dscale*p_Da(k)
          END DO
        END DO
      ELSE    
        ! Transposed B-matrix.
        ! Loop over the DOF's to extract
        DO j=1,UBOUND(ImatrixPositions,2)
          irowPosB = KldB(idestOffsetY-1+j)+idestOffsetX-1
          DO i=1,UBOUND(ImatrixPositions,1)
            
            ! Get the matrix position where to find our entry.
            k = ImatrixPositions(i,j)
            
            ! Get the matrix entry and write it to the destination matrix.
            Db(irowPosB+i-1) = dscale*p_Da(k)
          END DO
        END DO
      END IF
    
    END SUBROUTINE

    ! -----
    
    SUBROUTINE getLocalDefect (Dd,rmatrixrow,Rvectors,IvelocityDofs,IpressureDofs,&
        IchunkDofs,rmatrixPositions,bfirstTimestep,blastTimestep,&
        KcolA,KldA,KcolB,KldB,KcolC,KldC)
    
    ! Calculates the local defect in one row of a chunk.
    
    ! The matrices in the current row of the chunk.
    TYPE(t_matrixRow), INTENT(IN) :: rmatrixrow
    
    ! Pointers to the RHS and solution vectors in the current row of a chunk.
    ! the array must always have dimension 3, where Rvectors(2) identifies
    ! the RHS/solution that corresponds to the diagonal matrix.
    TYPE(t_chunkVector), DIMENSION(3), INTENT(IN) :: Rvectors
    
    ! The velocity DOF's under consideration
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IvelocityDofs
    
    ! The pressure DOF's under consideration
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IpressureDofs
    
    ! A list of all DOF's in the current chunk on the current element
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IchunkDofs
    
    ! A t_matrixPositions structure specifying the positions that
    ! are affected by the DOF's in the global matrices
    TYPE(t_matrixPositions), INTENT(IN) :: rmatrixPositions
    
    ! Flag that specifies if the current timestep is the first timestep
    ! in the chunk.
    LOGICAL, INTENT(IN) :: bfirstTimestep
    
    ! Flag that specifies if the current timestep is the last timestep
    ! in the chunk.
    LOGICAL, INTENT(IN) :: blastTimestep
    
    ! Column/row structure of the mass/Laplace matrices
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA

    ! Column/row structure of the B matrices
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolB
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldB

    ! Column/row structure of the C matrices. Optional.
    ! If not present, there is no C-matrix.
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: KcolC
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: KldC
    
    ! The local defect in one row of a chunk
    REAL(DP), DIMENSION(:), INTENT(OUT), TARGET :: Dd

      ! local variables
      INTEGER :: k
      REAL(DP), DIMENSION(:), POINTER :: p_Dd1
      REAL(DP), DIMENSION(:), POINTER :: p_Dd2
      REAL(DP), DIMENSION(:), POINTER :: p_Dd3
      REAL(DP), DIMENSION(:), POINTER :: p_Dd4
      REAL(DP), DIMENSION(:), POINTER :: p_Dd5
      REAL(DP), DIMENSION(:), POINTER :: p_Dd6

      ! Get the local RHS 
      DO k=1,SIZE(Dd)
        Dd(k) = Rvectors(2)%p_Db(IchunkDofs(k))
      END DO
      
      ! Initialise some pointers to the subvectors.
      ! Otherwise, one would be totally confused to program the following...
      p_Dd1 => Dd(iposLocalPrimalVelocityX:iposLocalPrimalVelocityY-1)
      p_Dd2 => Dd(iposLocalPrimalVelocityY:iposLocalPrimalPressure-1)
      p_Dd3 => Dd(iposLocalPrimalPressure:iposLocalDualVelocityX-1)
      p_Dd4 => Dd(iposLocalDualVelocityX:iposLocalDualVelocityY-1)
      p_Dd5 => Dd(iposLocalDualVelocityY:iposLocalDualPressure-1)
      p_Dd6 => Dd(iposLocalDualPressure:)
      
      ! That was the easy part, now the defect!
      !
      ! We do now by hand a kind of matrix-vector multiplication.
      ! Every matrix that is defined in rmatrixrow must be multiplied by
      ! the corresponding subvector and be subtracted from the RHS
      ! to get the local defect.
      !
      ! Primal solution
      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,1,0),&
          rmatrixrow%p_Da11,KcolA,KldA,p_Dd1,Rvectors(2)%p_Dx1)

      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,1,0),&
          rmatrixrow%p_Da21,KcolA,KldA,p_Dd2,Rvectors(2)%p_Dx1)
          
      CALL localMatrixVectorBT (IvelocityDofs,IpressureDofs,&
          rmatrixrow%DscaleFactor(3,1,0),rmatrixrow%p_Da31,&
          KcolB,rmatrixPositions%p_IposB,p_Dd3,Rvectors(2)%p_Dx1)          



      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,1,0),&
          rmatrixrow%p_Da41,KcolA,KldA,p_Dd4,Rvectors(2)%p_Dx1)

      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,1,0),&
          rmatrixrow%p_Da51,KcolA,KldA,p_Dd5,Rvectors(2)%p_Dx1)



      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,2,0),&
          rmatrixrow%p_Da12,KcolA,KldA,p_Dd1,Rvectors(2)%p_Dx2)

      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,2,0),&
          rmatrixrow%p_Da22,KcolA,KldA,p_Dd2,Rvectors(2)%p_Dx2)

      CALL localMatrixVectorBT (IvelocityDofs,IpressureDofs,&
          rmatrixrow%DscaleFactor(3,2,0),rmatrixrow%p_Da32,&
          KcolB,rmatrixPositions%p_IposB,p_Dd3,Rvectors(2)%p_Dx2)          


      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,2,0),&
          rmatrixrow%p_Da42,KcolA,KldA,p_Dd4,Rvectors(2)%p_Dx2)

      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,2,0),&
          rmatrixrow%p_Da52,KcolA,KldA,p_Dd5,Rvectors(2)%p_Dx2)


      CALL localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(1,3,0),&
          rmatrixrow%p_Da13,KcolB,KldB,p_Dd1,Rvectors(2)%p_Dx3)

      CALL localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(2,3,0),&
          rmatrixrow%p_Da23,KcolB,KldB,p_Dd2,Rvectors(2)%p_Dx3)

      IF (PRESENT(KcolC)) &
        CALL localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(3,3,0),&
          rmatrixrow%p_Da33,KcolC,KldC,p_Dd3,Rvectors(2)%p_Dx3)

      ! Dual solution
      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,4,0),&
          rmatrixrow%p_Da14,KcolA,KldA,p_Dd1,Rvectors(2)%p_Dx4)

      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,4,0),&
          rmatrixrow%p_Da24,KcolA,KldA,p_Dd2,Rvectors(2)%p_Dx4)



      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,4,0),&
          rmatrixrow%p_Da44,KcolA,KldA,p_Dd4,Rvectors(2)%p_Dx4)

      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,4,0),&
          rmatrixrow%p_Da54,KcolA,KldA,p_Dd5,Rvectors(2)%p_Dx4)

      CALL localMatrixVectorBT (IvelocityDofs,IpressureDofs,&
          rmatrixrow%DscaleFactor(6,4,0),&
          rmatrixrow%p_Da64,KcolB,rmatrixPositions%p_IposB,p_Dd6,Rvectors(2)%p_Dx4)          


      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,5,0),&
          rmatrixrow%p_Da15,KcolA,KldA,p_Dd1,Rvectors(2)%p_Dx5)

      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,5,0),&
          rmatrixrow%p_Da25,KcolA,KldA,p_Dd2,Rvectors(2)%p_Dx5)



      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,5,0),&
          rmatrixrow%p_Da45,KcolA,KldA,p_Dd4,Rvectors(2)%p_Dx5)

      CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,5,0),&
          rmatrixrow%p_Da55,KcolA,KldA,p_Dd5,Rvectors(2)%p_Dx5)

      CALL localMatrixVectorBT (IvelocityDofs,IpressureDofs,&
          rmatrixrow%DscaleFactor(6,5,0),&
          rmatrixrow%p_Da65,KcolB,rmatrixPositions%p_IposB,p_Dd6,Rvectors(2)%p_Dx5)          


      CALL localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(4,6,0),&
          rmatrixrow%p_Da46,KcolB,KldB,p_Dd4,Rvectors(2)%p_Dx6)

      CALL localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(5,6,0),&
          rmatrixrow%p_Da56,KcolB,KldB,p_Dd5,Rvectors(2)%p_Dx6)

      IF (PRESENT(KcolC)) &
        CALL localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(6,6,0),&
          rmatrixrow%p_Da66,KcolC,KldC,p_Dd6,Rvectors(2)%p_Dx6)
      
      ! Mass matrix coupling the primal velocity
      IF (.NOT. bfirstTimestep) THEN
        CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,1,-1),&
          rmatrixrow%p_Dp11,KcolA,KldA,p_Dd1,Rvectors(1)%p_Dx1)

        CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,1,-1),&
          rmatrixrow%p_Dp21,KcolA,KldA,p_Dd2,Rvectors(1)%p_Dx1)

        CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,2,-1),&
          rmatrixrow%p_Dp12,KcolA,KldA,p_Dd1,Rvectors(1)%p_Dx2)

        CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,2,-1),&
          rmatrixrow%p_Dp22,KcolA,KldA,p_Dd2,Rvectors(1)%p_Dx2)
      END IF

      ! Mass matrix coupling the dual velocity
      IF (.NOT. blastTimestep) THEN
        CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,4,1),&
          rmatrixrow%p_Dd44,KcolA,KldA,p_Dd4,Rvectors(3)%p_Dx4)

        CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,4,1),&
          rmatrixrow%p_Dd54,KcolA,KldA,p_Dd5,Rvectors(3)%p_Dx4)

        CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,5,1),&
          rmatrixrow%p_Dd45,KcolA,KldA,p_Dd4,Rvectors(3)%p_Dx5)

        CALL localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,5,1),&
          rmatrixrow%p_Dd55,KcolA,KldA,p_Dd5,Rvectors(3)%p_Dx5)
      END IF
      
    END SUBROUTINE


    SUBROUTINE addLocalCorrection (Dd,rvectors,IchunkDofs,domega)
    
    ! Adds a local correction vector to the global solution
    
    ! A list of all DOF's in the current chunk on the current element
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IchunkDofs
    
    ! The local defect in one row of a chunk
    REAL(DP), DIMENSION(:), INTENT(IN), TARGET :: Dd

    ! Damping parameter; standard = 1.0
    REAL(DP), INTENT(IN) :: domega
    
    ! Pointers to the RHS and solution vectors in the current row of a chunk.
    ! The solution in here is updated.
    TYPE(t_chunkVector), INTENT(INOUT) :: rvectors
    
      ! local variables
      INTEGER :: k

      ! Get the local RHS 
      DO k=1,SIZE(Dd)
        rvectors%p_Dx(IchunkDofs(k)) = rvectors%p_Dx(IchunkDofs(k)) + domega*Dd(k)
      END DO
      
    END SUBROUTINE
    
    SUBROUTINE localMatrixVector (Idofs,dscale,p_Da,Kcol,Kld,Db,Dx)
    
    ! Performs a local matrix vector multiplication
    !   b = b-Ax
    ! for the rows defined in the array Idofs.
    
    ! The rows in the matrix that correspond to the elements in Db.
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: Idofs
    
    ! Scaling factor for the matrix
    REAL(DP), INTENT(IN) :: dscale
    
    ! Matrix data, structure 9. the data array is specified as a pointer;
    ! if set to NULL(), there is no matrix, so no matrix vector multiplication
    ! will be done.
    REAL(DP), DIMENSION(:), POINTER :: p_Da
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld
    
    ! Solution vector; is multiplied by the matrix and subtracted from b.
    REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
    
    ! RHS vector. The size must be the same as Idofs.
    ! Is overwritten by the defect.
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: Db
    
      INTEGER :: i,j
      
      IF (ASSOCIATED(p_Da) .AND. (dscale .NE. 0.0_DP)) THEN
      
        ! Loop through all DOF's. These defines the lines in Db where we have
        ! to do b-Ax.
        DO i=1,SIZE(Idofs)
          ! Subtract row IDofs(i) * Dx from Db.
          DO j=Kld(Idofs(i)),Kld(Idofs(i)+1)-1
            Db(i) = Db(i) - dscale*p_Da(j)*Dx(Kcol(j))
          END DO
        END DO
        
      END IF
    
    END SUBROUTINE

    SUBROUTINE localMatrixVectorBT (IvelocityDofs,IpressureDofs,dscale,&
        p_DB,KcolB,ImatrixPos,Db,Dx)
    
    ! Performs a local matrix vector multiplication
    !   b = b - B^T x
    ! for the rows defined in the array Idofs.
    
    ! Velocity DOF's under consideration
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IvelocityDofs

    ! Pressure DOF's under consideration
    INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IpressureDofs
    
    ! Scaling factor for the matrix
    REAL(DP), INTENT(IN) :: dscale
    
    ! Matrix data, structure 9. the data array is specified as a pointer;
    ! if set to NULL(), there is no matrix, so no matrix vector multiplication
    ! will be done.
    REAL(DP), DIMENSION(:), POINTER :: p_DB
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolB
    
    ! Index array specifying the entries (rows) in the B-matrix which are affected
    ! by the DOF's in IvelocityDofs. DIMENSION(#velocity DOF's,#pressure DOF's).
    INTEGER(PREC_MATIDX), DIMENSION(:,:), INTENT(IN) :: ImatrixPos
    
    ! Solution vector; is multiplied by the matrix and subtracted from b.
    REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
    
    ! RHS vector. The size must be the same as Idofs.
    ! Is overwritten by the defect.
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: Db
    
      INTEGER :: i,j
      INTEGER(PREC_MATIDX) :: k
      
      IF (ASSOCIATED(p_DB) .AND. (dscale .NE. 0.0_DP)) THEN
      
        ! Loop through all pressure DOF's. These defines the lines in DB where we have
        ! to do b-Ax.
        DO i = 1,SIZE(IpressureDofs)
          
          ! Every pressure DOF gives us a line in the B^T matrix tat we have
          ! to process. In that line, we have to loop through all velocity DOF's.
          DO j = 1,SIZE(IvelocityDofs)
          
            ! Find the matrix position in the B matrix specifying that is affected
            ! by the current velocity DOF. This defines the entry that we have to
            ! take as the entry of B^T.
            k = ImatrixPos (j,i)
            
            ! Multiply the B^T entry with the velocity DOF and subtract.
            Db(i) = Db(i) - dscale*p_DB(k)*Dx(IvelocityDofs(j))
          
          END DO
         
        END DO
        
      END IF
      
      ! Note: the above loop looks wrong on the first look. It subtracts only
      ! those entries of B^T that are specified in the ImatrixPos array --
      ! and ImatrixPpos specifies only the entries of the local matrix!
      ! Whoever tolds one that the local matrix contains all elements in B^T?!?
      !
      ! Well the FE theory does -- at least for Q0 or PQ1! Every row in B 
      ! tells a velocity DOF how to couple to all the surrounding pressure DOF's.
      ! the pressure DOF's now exists only on one cell, so coupling only to the 
      ! DOF's on that cell and not to those on the neighbour cells!
      ! Therefore, the B matrix contains only the couplings to the velocity DOF's
      ! on the current cell, and thus the B^T matrix does as well.
      ! So the B^T does not contain any entries that couple to other DOF's than
      ! those velocity DOF's on the current cell -- and these are exactly the
      ! DOF's we have in IvelocityDofs. That implies that the local B^T matrix
      ! indeed contains all entries of a line in the B^T matrix!
    
    END SUBROUTINE

  END SUBROUTINE

END MODULE
