!##############################################################################
!# ****************************************************************************
!# <name> spacetimevanca </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module spacetimevanca

  use fsystem
  use genoutput
  use externalstorage
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use collection
  use vectorio
  
  !use spacepreconditionerinit
  !use spacetimediscretisation
  use spacetimevectors
  use spacetimelinearsystem
  
  use matrixio
  
  implicit none

  ! A structure that collects the pointers to all submatrices in a matrix
  ! row of the global block matrix. Circumvents the problem
  ! that Fortran cannot set up an array of pointers...
  type t_matrixRow
  
    ! Pointers for submatrices below the diagonal
    real(DP), dimension(:), pointer :: p_Dp11 => null()
    real(DP), dimension(:), pointer :: p_Dp21 => null()
    real(DP), dimension(:), pointer :: p_Dp12 => null()
    real(DP), dimension(:), pointer :: p_Dp22 => null()
                                              
    ! Pointers for submatrices of a diagonal bock
    real(DP), dimension(:), pointer :: p_Da11 => null()
    real(DP), dimension(:), pointer :: p_Da21 => null()
    real(DP), dimension(:), pointer :: p_Da31 => null()
    real(DP), dimension(:), pointer :: p_Da41 => null()
    real(DP), dimension(:), pointer :: p_Da51 => null()
                                              
    real(DP), dimension(:), pointer :: p_Da12 => null()
    real(DP), dimension(:), pointer :: p_Da22 => null()
    real(DP), dimension(:), pointer :: p_Da32 => null()
    real(DP), dimension(:), pointer :: p_Da42 => null()
    real(DP), dimension(:), pointer :: p_Da52 => null()
                                              
    real(DP), dimension(:), pointer :: p_Da13 => null()
    real(DP), dimension(:), pointer :: p_Da23 => null()
    real(DP), dimension(:), pointer :: p_Da33 => null()
                                              
    real(DP), dimension(:), pointer :: p_Da14 => null()
    real(DP), dimension(:), pointer :: p_Da24 => null()
    real(DP), dimension(:), pointer :: p_Da44 => null()
    real(DP), dimension(:), pointer :: p_Da54 => null()
    real(DP), dimension(:), pointer :: p_Da64 => null()
                                              
    real(DP), dimension(:), pointer :: p_Da15 => null()
    real(DP), dimension(:), pointer :: p_Da25 => null()
    real(DP), dimension(:), pointer :: p_Da45 => null()
    real(DP), dimension(:), pointer :: p_Da55 => null()
    real(DP), dimension(:), pointer :: p_Da65 => null()
                                              
    real(DP), dimension(:), pointer :: p_Da46 => null()
    real(DP), dimension(:), pointer :: p_Da56 => null()
    real(DP), dimension(:), pointer :: p_Da66 => null()
                                              
    ! Pointers for submatrices above the diagonal
    real(DP), dimension(:), pointer :: p_Dd44 => null()
    real(DP), dimension(:), pointer :: p_Dd54 => null()
    real(DP), dimension(:), pointer :: p_Dd45 => null()
    real(DP), dimension(:), pointer :: p_Dd55 => null()

    ! Scaling factors for all the submatrices.
    ! DIMENSION(6,6,-1..1), with 6=#matrix rows/columns. The last dimension specifies the
    ! matrix column: -1=left from the diagonal, 0=diagonal, 1=right from the diagonal.
    real(DP), dimension(6,6,-1:1) :: DscaleFactor = 0.0_DP
    
  end type
  
  ! Collects the positions of the DOF's on one element in all the matrices
  type t_matrixPositions
  
    ! Matrix positions in the FEM matrices.
    integer, dimension(:,:), pointer :: p_IposA  => null()
    
    ! Matrix positions in the B-matrices.
    integer, dimension(:,:), pointer :: p_IposB  => null()
    
    ! matrix positions in the C matrix.
    integer, dimension(:,:), pointer :: p_IposC  => null()
  
  end type
  
  ! A structure that collects the pointers to all subvectors
  ! of a solution chunk. Circumvents the problem
  ! that Fortran cannot set up an array of pointers...
  type t_chunkVector
    
    ! Pointer to solution and RHS vector
    real(DP), dimension(:), pointer :: p_Dx
    real(DP), dimension(:), pointer :: p_Db
    
    ! Pointer to the subvectors of the solution vector
    real(DP), dimension(:), pointer :: p_Dx1
    real(DP), dimension(:), pointer :: p_Dx2
    real(DP), dimension(:), pointer :: p_Dx3
    real(DP), dimension(:), pointer :: p_Dx4
    real(DP), dimension(:), pointer :: p_Dx5
    real(DP), dimension(:), pointer :: p_Dx6
    
  end type

    
contains

  subroutine createLocalMatrixQ1TQ0 (rdiscretisation,ntimesteps,rmatrix,h_ImatrixMapping)
  
  ! Creates a local matrix containing the submatrices of ntimesteps timesteps.
    
  ! Block discretisation structure of every timestep
  type(t_blockDiscretisation), intent(IN) :: rdiscretisation
  
  ! Number of timesteps, the local matrix should hold.
  integer, intent(IN) :: ntimesteps

  ! The local matrix to create
  type(t_matrixScalar), intent(OUT) :: rmatrix
  
  ! A handle to an array with DIMENSION(:,:). For every matrix block column in the
  ! local matrix, this points to the first entry of the column. This array is
  ! therefore a way to directly access the block columns in the matrix.
  integer, intent(OUT) :: h_ImatrixMapping
  
    ! local variables
    integer(I32), dimension(:,:), pointer :: p_Iindex
    integer :: ndofLocal,i,j,ielType,ndofGlobal,itimestep,ndofLocalPrimalDual
    integer :: ndofVelocity,ndofPressure, iYoffset,isum1,isum2
    integer, dimension(2) :: Isize
    real(DP), dimension(:), pointer :: p_Da
    
    ! Get the number of local DOF's in every timestep. For that purpose,
    ! loop over the blocks in the discretisation structure, ask the
    ! element about its number of local DOF's and sum that stuff up.
    ndofLocal = 0
    do i=1,rdiscretisation%ncomponents
      ielType = rdiscretisation%RspatialDiscr(i)%&
          RelementDistr(1)%celement
      ndofLocal = ndofLocal + elem_igetNDofLoc(ielType)
    end do
    
    ! Get the number of local DOF's in the primal space. This coincodes with
    ! the number of local DOF's in the dual space.
    ! As they both summed up give the number of local DOF's, we simply have to
    ! divide that by 2.
    ndofLocalPrimalDual = ndofLocal / 2
    
    ! Get the number of velocity DOF's in the primal space.
    ielType = rdiscretisation%RspatialDiscr(1)%&
        RelementDistr(1)%celement
    ndofVelocity = elem_igetNDofLoc(ielType)

    ielType = rdiscretisation%RspatialDiscr(3)%&
        RelementDistr(1)%celement
    ndofPressure = elem_igetNDofLoc(ielType)
    
    ! Multiplying that with the number of timesteps (+1) gives the size of the
    ! local matrix.
    ndofGlobal = ndofLocal * (ntimesteps + 1)
    
    ! Allocate memory for the pointer array that later allows quick access
    ! to the matrix columns. Each matrix is a 6x6 block matrix.
    Isize = (/6*(ntimesteps + 1),6*3/)
    call storage_new ('createLocalMatrixQ1TQ0', 'ImatrixMapping', Isize, &
        ST_INT, h_ImatrixMapping,ST_NEWBLOCK_ZERO)
    call storage_getbase_int2d(h_ImatrixMapping,p_Iindex)
    
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
    
    call lsyssc_createFullMatrix (rmatrix,.true.,ndofGlobal)
    
    ! Get the data array
    call lsyssc_getbase_double (rmatrix,p_Da)
    
    ! Loop through all diagonal blocks to fill in data.
    do itimestep = 0,ntimesteps
      ! Primal matrix
      call markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal,&
          1+itimestep*ndofLocal,&
          ndofLocalPrimalDual)

      ! Dual matrix
      call markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual,&
          ndofLocalPrimalDual)

      ! Mass matrices above the diagonal
      call markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual,&
          ndofVelocity)
      call markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal+ndofVelocity,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual+ndofVelocity,&
          ndofVelocity)

      ! The R-block below the diagonal
      call markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual,&
          1+itimestep*ndofLocal,&
          2*ndofVelocity)
          
    end do
    
    ! What is still missing are the mass matrices that couple the timesteps.
    ! Loop through all diagonal blocks to fill in data.
    do itimestep = 0,ntimesteps-1
      ! Primal matrix
      call markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal+ndofLocal,&
          1+itimestep*ndofLocal,&
          2*ndofVelocity)
      !CALL markSubmatrix (p_Da,ndofGlobal,&
      !    1+itimestep*ndofLocal+ndofLocal+ndofVelocity,&
      !    1+itimestep*ndofLocal+ndofVelocity,&
      !    2*ndofVelocity)

      ! Dual matrix
      call markSubmatrix (p_Da,ndofGlobal,&
          1+itimestep*ndofLocal+ndofLocalPrimalDual,&
          1+itimestep*ndofLocal+ndofLocal+ndofLocalPrimalDual,&
          2*ndofVelocity)
      !CALL markSubmatrix (p_Da,ndofGlobal,&
      !    1+itimestep*ndofLocal+ndofLocalPrimalDual+ndofVelocity,&
      !    1+itimestep*ndofLocal+ndofLocal+ndofLocalPrimalDual+ndofVelocity,&
      !    ndofVelocity)
    end do
    
    ! Convert the matrix to format 9
    call lsyssc_convertMatrix (rmatrix,LSYSSC_MATRIX9)
    
    ! Set up the matrix column pointers. FOr that purpose, at first
    ! write into every block the size of the block -- if it contains data.
    do itimestep = 0,ntimesteps
      
      iYoffset = 6*itimestep

      if (itimestep .gt. 0) then
        ! Left' matrix below the diagonal
        p_Iindex(iYoffset+1,0+1) = ndofVelocity
        p_Iindex(iYoffset+1,0+2) = ndofVelocity
                                     
        p_Iindex(iYoffset+2,0+1) = ndofVelocity
        p_Iindex(iYoffset+2,0+2) = ndofVelocity
      end if

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
                
      if (itimestep .lt. ntimesteps) then
        ! 'Right' matrix, above the diagonal
        p_Iindex(iYoffset+4,12+4) = ndofVelocity
        p_Iindex(iYoffset+4,12+5) = ndofVelocity
                                     
        p_Iindex(iYoffset+5,12+4) = ndofVelocity
        p_Iindex(iYoffset+5,12+5) = ndofVelocity
      end if
    end do
    
    ! Now loop through all the indices and sum them up.
    ! This gives the matrix pointers in the local matrices.
    !
    ! First shift the indices and set the first column to 0.
    do j=6*3,2,-1
      do i=1,ubound(p_Iindex,1)
        p_Iindex(i,j) = p_Iindex(i,j-1)
      end do
    end do
    
    do i=1,ubound(p_Iindex,1)
      p_Iindex(i,1) = 0.0_DP
    end do
    
    ! Now sum that stuff up
    do j=2,6*3
      do i=1,ubound(p_Iindex,1)
        p_Iindex(i,j) = p_Iindex(i,j) + p_Iindex(i,j-1)
      end do
    end do
    
    ! And add a 1. That gives the start pointers for the blocks.
    do j=1,6*3
      do i=1,ubound(p_Iindex,1)
        p_Iindex(i,j) = 1 + p_Iindex(i,j)
      end do
    end do
    
  contains
  
    subroutine markSubmatrix (Da,neq,iposY,iposX,isize)
    
    ! Initialises a subblock if the matrix Da with 1.0.
    ! The block starts at position (iposY,iposX) and has size isize*isize.
    
    ! Data array
    real(DP), dimension(:), intent(INOUT) :: Da
    
    ! Size of the matrix
    integer, intent(IN) :: neq
    
    ! Start Position where to fill the matrix with 1.
    integer, intent(IN) :: iposY,iposX
    
    ! Size of the block that should be filled with 1.
    integer,intent(IN) :: isize
    
      ! local variables
      integer :: i,j
      
      do j=iposX-1,iposX+isize-2
        do i=iposY,iposY+isize-1
          Da(j*neq+i) = 1.0_DP
        end do
      end do
      
    end subroutine
  
  end subroutine

  ! ***************************************************************************

  subroutine calcMatrixIndices (Kcol,Kld,IdofsTrial,IdofsTest,Iidx)
  
  ! Calculates the indices inside of those entries in a matrix that
  ! belong to the local matrix. The matrix must be given in
  ! matrix structure 9.
  
  ! The column structure of the matrix.
  integer, dimension(:), intent(IN) :: Kcol

  ! The row structure of the matrix.
  integer, dimension(:), intent(IN) :: Kld
  
  ! The DOF's in the trial space that belong to the local matrix.
  ! These correspond to the columns of the local matrix.
  integer, dimension(:), intent(IN) :: IdofsTrial

  ! The DOF's in the test space that belong to the local matrix.
  ! THese correspond to the rows of the local matrix.
  integer, dimension(:), intent(IN) :: IdofsTest
  
  ! A 2D array containing the positions in the global matrix that
  ! belong to the local matrix. The index array is set up transposed
  ! to get quicker memory access!
  ! Iidx must be a square matrix with
  ! DIMENSION(size(IdofsTest),size(IdofsTrial))
  integer, dimension(:,:), intent(OUT) :: Iidx
  
    ! local variables
    integer :: i,k
    integer :: j
    integer :: idoflocal,icol
    
    ! Loop over the matrix rows.
    do i=1,size(IdofsTrial)
      
      idofLocal = IdofsTrial(i)
      
      ! In the row loop over the columns to find those that belong to
      ! local DOF's.
      do j=Kld(idofLocal),Kld(idofLocal+1)-1
      
        ! Column number
        icol = Kcol(j)
      
        ! Find the column in the Idofs array. If we find it, transfer its
        ! position to Iidx.
        do k=1,size(IdofsTest)
        
          if (icol .eq. IdofsTest(k)) then
            ! Put the index position to the local matrix
            Iidx(k,i) = j
            exit
          end if
        
        end do
      
      end do
    
    end do
    
  end subroutine
  
  
  
  subroutine ccopt_precSpaceTimeVanca (rproblem,rmatrix,rx,rb,domega,nchunkSize,niterations)
  
  ! Performs niterations VANKA iterations to enhance the vector rx.
  ! nchunkSize specifies the number of timesteps that are collected
  ! to a local system.
  
  ! The main problem structure
  type(t_problem), intent(INOUT) :: rproblem
  
  ! The underlying space time matrix.
  type(t_ccoptSpaceTimeMatrix), intent(IN) :: rmatrix
  
  ! RHS vector
  type(t_spacetimeVector), intent(IN) :: rb
  
  ! Damping parameter
  real(DP), intent(IN) :: domega
  
  ! Chunk size. >= 0. A chunk size >= #timesteps will apply VANKA
  ! to all timesteps simultaneously, which is VERY memory and hard disc
  ! intensive! A chunk size of 0 will apply the standard VANKA for
  ! optimal control problems that does not combine multiple time steps.
  integer, intent(IN) :: nchunksize
  
  ! Number of VANKA iterations to apply.
  integer, intent(IN) :: niterations
  
  ! Solution vector to be updated
  type(t_spacetimeVector), intent(INOUT) :: rx
  
    ! Constants
    integer, parameter :: ndofLocalVelocityX = 4
    integer, parameter :: ndofLocalVelocityY = 4
    integer, parameter :: ndofLocalVelocity  = ndofLocalVelocityX + ndofLocalVelocityY
    integer, parameter :: ndofLocalPressure  = 1
    integer, parameter :: ndofLocalHalf = ndofLocalVelocity + ndofLocalPressure
    integer, parameter :: ndofLocal = 2*ndofLocalHalf   ! primal and dual
    
    integer, parameter :: iposLocalPrimalVelocityX = 1
    integer, parameter :: iposLocalPrimalVelocityY = &
                              iposLocalPrimalVelocityX + ndofLocalVelocityX
    integer, parameter :: iposLocalPrimalPressure = &
                              iposLocalPrimalVelocityY + ndofLocalVelocityY

    integer, parameter :: iposLocalDualVelocityX = &
                              iposLocalPrimalPressure + ndofLocalPressure
    integer, parameter :: iposLocalDualVelocityY = &
                              iposLocalDualVelocityX + ndofLocalVelocityX
    integer, parameter :: iposLocalDualPressure = &
                              iposLocalDualVelocityY + ndofLocalVelocityY

    ! local variables
    integer :: iposGlobalPrimalVelocityX,iposGlobalPrimalVelocityY
    integer :: iposGlobalPrimalPressure
    integer :: iposGlobalDualVelocityX,iposGlobalDualVelocityY
    integer :: iposGlobalDualPressure
    integer :: iel,NEL
    integer :: NVT
    integer :: NEQ
    integer(PREC_EDGEIDX), dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(ndofLocalVelocityX) :: IvelocityDofs
    integer, dimension(ndofLocalPressure) :: IpressureDofs
    !integer, DIMENSION(ndofLocalVelocity,ndofLocalVelocity) :: ImatPosVel
    !integer, DIMENSION(ndofLocalPressure,ndofLocalVelocity) :: ImatPosPres
    integer, dimension(:), pointer :: p_KcolA,p_KcolB,p_KcolC
    integer, dimension(:), pointer :: p_KldA,p_KldB,p_KldC
    real(DP), dimension(ndofLocal) :: DxLocal, DbLocal
    
    integer :: nblockSize,i,j,k,ichunk,NEQtime,ichunkpos,istep,ipos,iiteration
    integer :: h_ImatrixMapping
    
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Array of system matrices for the timesteps.
    type(t_matrixBlock), dimension(:,:), allocatable :: RsystemMatrix
    
    ! Arrays with pointers to the matrix data of the matrices in RsystemMatrix
    type(t_matrixRow), dimension(:), allocatable :: RmatrixRows
    
    ! The same for the solution/rhs vector
    type(t_chunkVector), dimension(:), allocatable :: Rvectors
    
    ! Global DOF's in one block of the current chunk
    integer, dimension(:), allocatable :: IchunkDofs
    
    ! Matrix position structure.
    ! Here, the positions of the matrix elements are saved, corresponding to
    ! the current element.
    type(t_matrixPositions) :: rmatPositions
    
    ! Array of vectors simultaneously to handle for all the timesteps
    type(t_vectorBlock), dimension(:), allocatable :: RrhsGlobal,RsolutionGlobal
    type(t_vectorBlock) :: rvectorTemp
    
    ! local matrix,vectors
    type(t_matrixScalar) :: rmatrixLocal
    real(DP), dimension(:), pointer :: p_DaLocal
    integer, dimension(:), pointer :: p_KcolLocal
    integer, dimension(:), pointer :: p_KldLocal
    type(t_vectorScalar) :: rrhsLocal
    real(DP), dimension(:), pointer :: p_DrhsLocal
    integer(I32), dimension(:,:), pointer :: p_ImatrixMapping
    type(t_matrixBlock), dimension(1) :: RmatrixLocalBlock
    type(t_vectorBlock) :: rrhsBlock
    type(t_linsolNode), pointer :: p_rsolver
    integer :: ierror
    
    type(t_vectorBlock) :: rdbgrhs,rdbgsol
    type(t_matrixBlock) :: rdbgMatrix
    type(t_linsolNode), pointer :: p_rdbgsolver
    type(t_matrixBlock), dimension(1) :: RdbgmatrixLocalBlock
    
    ! Fetch some information
    p_rdiscretisation => rmatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation
    NEL = p_rdiscretisation%p_rtriangulation%NEL
    NEQ = dof_igetNDofGlobblock(p_rdiscretisation)
    NVT = p_rdiscretisation%p_rtriangulation%NVT
    call storage_getbase_int2d (p_rdiscretisation%p_rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    
    ! Calculate the number of block that must be handled simultaneously.
    NEQtime = rmatrix%p_rspaceTimeDiscr%NEQtime
    nblockSize = min(NEQtime,max(1,nchunkSize+1))

    ! Allocate memory for global matrices.
    !
    ! System matrices on the diagonal of the global matrix.
    ! RsystemMatrix(-1) describes the matrices coupling the primal velocity
    ! and roughly contains mass matrices.
    ! RsystemMatrix(0) contains the system matrices on the main diagonal
    ! of the global matrix.
    ! RsystemMatrix(1) describes the matrices coupling the dual velocity
    ! and roughly contains mass matrices.
    allocate(RsystemMatrix(-1:1,nblockSize))
    
    ! Allocate one additional element before and after the next arrays. These elements
    ! are unused but there existence simplifies the handling of the array
    ! at the beginning and the ending of a chunk.
    allocate(RrhsGlobal(0:nblockSize+1))
    allocate(RsolutionGlobal(0:nblockSize+1))
    allocate(Rvectors(0:nblockSize+1))
    
    ! RmatrixRows contains pointers to the matrix data arrays of the matrices in
    ! RsystemMatrix. That gives some speed, as we don't have to access the
    ! storage-module so often.
    allocate(RmatrixRows(nblockSize))
    
    ! IchunkDofs holds the numbers of the global DOF's in one block of the current chunk.
    allocate(IchunkDofs(1:ndofLocal))
    
    ! Allocate all the matrices/vectors we have to handle simultaneously.
    do j=1,nblockSize
      do i=-1,1
        print *,'Not implemented!'
        stop
        !call cc_allocSystemMatrix (rproblem,rmatrix%p_rspaceTimeDiscr%p_rlevelInfo,&
        !  RsystemMatrix(i,j))
      end do
      
      call lsysbl_createVecBlockByDiscr (p_rdiscretisation,RrhsGlobal(j))
      call lsysbl_createVecBlockByDiscr (p_rdiscretisation,RsolutionGlobal(j))
      
      ! Get pointers to the vectors for quicker access
      call lsysbl_getbase_double (RrhsGlobal(j),Rvectors(j)%p_Db)
      call lsysbl_getbase_double (RsolutionGlobal(j),Rvectors(j)%p_Dx)

      call lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(1),Rvectors(j)%p_Dx1)
      call lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(2),Rvectors(j)%p_Dx2)
      call lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(3),Rvectors(j)%p_Dx3)
      call lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(4),Rvectors(j)%p_Dx4)
      call lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(5),Rvectors(j)%p_Dx5)
      call lsyssc_getbase_double (RsolutionGlobal(j)%RvectorBlock(6),Rvectors(j)%p_Dx6)
    end do
    call lsysbl_createVecBlockByDiscr (p_rdiscretisation,rvectorTemp)
    
    ! Allocate memory for the matrix positions.
    ! These change with every element and allow us to quickly access for every
    ! entry in each vector the corresponding matrix entry.
    allocate(rmatPositions%p_IposA(ndofLocalVelocityX,ndofLocalVelocityY))
    allocate(rmatPositions%p_IposB(ndofLocalVelocityX,ndofLocalPressure))
    allocate(rmatPositions%p_IposC(ndofLocalPressure,ndofLocalPressure))
    
    ! Allocate memory for local matrices
    call createLocalMatrixQ1TQ0 (&
      p_rdiscretisation,nblockSize-1,rmatrixLocal,h_ImatrixMapping)
      
    call storage_getbase_int2d (h_ImatrixMapping,p_ImatrixMapping)
      
    ! Get pointers to the matrix data.
    call lsyssc_getbase_double (rmatrixLocal,p_DaLocal)
    call lsyssc_getbase_Kcol (rmatrixLocal,p_KcolLocal)
    call lsyssc_getbase_Kld (rmatrixLocal,p_KldLocal)
      
    ! Create appropriate RHS and solution vectors.
    call lsyssc_createVecIndMat (rmatrixLocal,rrhsLocal,.false.)
    
    ! Get pointers to the data
    call lsyssc_getbase_double (rrhsLocal,p_DrhsLocal)
    
    ! Create derived 1x1 block matrices/vectors for the solver.
    call lsysbl_createMatFromScalar (rmatrixLocal,RmatrixLocalBlock(1))
    call lsysbl_createVecFromScalar (rrhsLocal,rrhsBlock)
    
    ! With the matrix, prepare an UMFPACK solver to solve our system
    call linsol_initUMFPACK4(p_rsolver)
    call linsol_setMatrices(p_rsolver,RmatrixLocalBlock)
    
    call lsyssc_clearMatrix (rmatrixLocal,1.0_DP)
    call linsol_initStructure (p_rsolver,ierror)
    if (ierror .ne. 0) then
      print *,'UMFPACK for local system cannot be initialised!'
      call sys_halt()
    end if
    
    ! We assume all matrices to share the same structure! (I.e. all velocity
    ! and mass matrices on the one hand, all B and B^T matrices on the other hand.)
    ! Fetch the matrix structure arrays.
    call lsyssc_getbase_Kcol (RsystemMatrix(0,1)%RmatrixBlock(1,1),p_KcolA)
    call lsyssc_getbase_Kld (RsystemMatrix(0,1)%RmatrixBlock(1,1),p_KldA)
    call lsyssc_getbase_Kcol (RsystemMatrix(0,1)%RmatrixBlock(1,3),p_KcolB)
    call lsyssc_getbase_Kld (RsystemMatrix(0,1)%RmatrixBlock(1,3),p_KldB)
    call lsyssc_getbase_Kcol (RsystemMatrix(0,1)%RmatrixBlock(3,3),p_KcolC)
    call lsyssc_getbase_Kld (RsystemMatrix(0,1)%RmatrixBlock(3,3),p_KldC)

    ! Set up index positions to quickly find entries in a spatial vector.
    iposGlobalPrimalVelocityX = 1
    iposGlobalPrimalVelocityY = iposGlobalPrimalVelocityX + &
        dof_igetNDofGlob(p_rdiscretisation%RspatialDiscr(1))
    iposGlobalPrimalPressure = iposGlobalPrimalVelocityY + &
        dof_igetNDofGlob(p_rdiscretisation%RspatialDiscr(2))
    iposGlobalDualVelocityX = iposGlobalPrimalPressure + &
        dof_igetNDofGlob(p_rdiscretisation%RspatialDiscr(3))
    iposGlobalDualVelocityY = iposGlobalDualVelocityX + &
        dof_igetNDofGlob(p_rdiscretisation%RspatialDiscr(4))
    iposGlobalDualPressure = iposGlobalDualVelocityY + &
        dof_igetNDofGlob(p_rdiscretisation%RspatialDiscr(5))

    ! Perforn niterations iterations
    do iiteration = 1,niterations

      ! Loop through all chunks. Every chunk consists of nchunkSize timesteps --
      ! rounded up; or more precisely, every chunk consists of nblocks vectors.
      do ichunk = 1,(NEQtime+nblockSize-1)/nblockSize
      
        ! Calculate the chunk position. This is normally ichunk*nblockSize
        ! except if this doesn't fit to the end position. In the latter case,
        ! we modify the chunk position appropriately such that we have
        ! nblockSize blocks in the chunk.
        ichunkPos = (ichunk-1)*nblockSize
        if ((ichunkPos+nblockSize) .gt. NEQtime) &
          ichunkPos = NEQtime-nblockSize
        ichunkPos = ichunkPos + 1
        
        ! For all items in the current chunk, get the global solution and
        ! RHS vectors.
        do ipos = 1,nblockSize
          call sptivec_getTimestepData (rx, ichunkPos+ipos-1, RsolutionGlobal(ipos))
          call sptivec_getTimestepData (rb, ichunkPos+ipos-1, RrhsGlobal(ipos))
        end do
        
        ! Build the matrices that belong to our chunk
        call buildMatrices (rproblem,rmatrix,&
            ichunkPos,nblockSize,RsystemMatrix,RmatrixRows)
        
        ! Set up "b-Mx" for the submatrix 'before' the chunk and 'after' the chunk;
        ! this gives the time coupling to the previous and next chunk.
        if (ichunkPos .gt. 1) then
          call sptivec_getTimestepData (rx, ichunkPos-1, rvectorTemp)
          !CALL matio_writeBlockMatrixHR (RsystemMatrix(-1,1), 'matleft',&
          !                          .TRUE., 0, 'matleft.txt', '(E10.2)')
          call lsysbl_blockMatVec (RsystemMatrix(-1,1),&
              rvectorTemp,RrhsGlobal(1),-1.0_DP,1.0_DP)
        end if
        
        if (ichunkPos+nblockSize-1 .lt. NEQtime) then
          call sptivec_getTimestepData (rx, ichunkPos+nblockSize-1+1, rvectorTemp)
          !CALL matio_writeBlockMatrixHR (RsystemMatrix(1,nblockSize), 'matright',&
          !                          .TRUE., 0, 'matright.txt', '(E10.2)')
          call lsysbl_blockMatVec (RsystemMatrix(1,nblockSize),&
              rvectorTemp,RrhsGlobal(nblockSize),-1.0_DP,1.0_DP)
        end if
        
        
        ! DEBUG!!!
        ! globalen Defekt bilden, UMFPACK anwerfen und Lösung draufaddieren.
        ! Sollte je nach Einstellung einem globalen UMFPACK oder einem
        ! UMFPACK in jd. Zeitschritt entsprechen.
        !...
        call lsysbl_createEmptyMatrix (rdbgmatrix,6*nblockSize)
        do j=1,nblockSize
          do i=1,3
            if (j .gt. 1) then
              call lsysbl_insertSubmatrix (&
                  RsystemMatrix(-1,j),rdbgmatrix,&
                  LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,6*(j-1)+1,6*(j-1)+1-1)
            end if
            call lsysbl_insertSubmatrix (&
                RsystemMatrix(0,j),rdbgmatrix,&
                LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,6*(j-1)+1,6*(j-1)+1)
            if (j .lt. nblockSize) then
              call lsysbl_insertSubmatrix (&
                  RsystemMatrix(1,j),rdbgmatrix,&
                  LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,6*(j-1)+1,6*(j-1)+1+1)
            end if
          end do
        end do
        call lsysbl_updateMatStrucInfo(rdbgmatrix)
        call matio_writeBlockMatrixHR (rdbgmatrix, 'fullmat',&
                                       .true., 0, 'matfull.txt', '(E10.2)')
                                       
        call linsol_initUmfpack4 (p_rdbgSolver)
        call lsysbl_duplicateMatrix (rdbgmatrix,RdbgmatrixLocalBlock(1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        call linsol_setMatrices (p_rdbgSolver,RdbgmatrixLocalBlock)
        call linsol_initStructure (p_rdbgSolver,ierror)
        if (ierror .ne. 0) then
          print *,'UMFPACK for local system cannot be initialised!'
          call sys_halt()
        end if
        call linsol_initData (p_rdbgSolver,ierror)
        if (ierror .ne. 0) then
          print *,'UMFPACK for local system cannot be initialised!'
          call sys_halt()
        end if
        
        call lsysbl_createVecBlockIndMat(rdbgmatrix,rdbgrhs,.true.)
        call lsysbl_createVecBlockIndMat(rdbgmatrix,rdbgsol,.true.)
        do i=1,nblockSize
          do j=1,6
            call lsyssc_copyVector (RrhsGlobal(i)%RvectorBlock(j),&
                rdbgrhs%RvectorBlock((i-1)*6+j))
            call lsyssc_copyVector (RsolutionGlobal(i)%RvectorBlock(j),&
                rdbgsol%RvectorBlock((i-1)*6+j))
          end do
        end do
        call lsysbl_blockMatVec (rdbgmatrix,rdbgsol,rdbgrhs,-1.0_DP,1.0_DP)
        
        call linsol_precondDefect (p_rdbgSolver,rdbgrhs)
        call linsol_releaseSolver(p_rdbgSolver)

        do i=1,nblockSize
          do j=1,6
            call lsyssc_vectorLinearComb (rdbgrhs%RvectorBlock((i-1)*6+j),&
                RsolutionGlobal(i)%RvectorBlock(j),0.9_DP,1.0_DP)
          end do
        end do
                                       
        call lsysbl_releaseVector (rdbgsol)
        call lsysbl_releaseVector (rdbgrhs)
        call lsysbl_releaseMatrix (rdbgmatrix)
        call lsysbl_releaseMatrix (RdbgmatrixLocalBlock(1))
        
        
        
        
        
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
        do ipos = 1,nblockSize
          call sptivec_setTimestepData (rx, ichunkPos+ipos-1, RsolutionGlobal(ipos))
        end do

      end do
      
    end do

    ! Release the solver
    call linsol_releaseSolver (p_rsolver)

    ! Release memory
    deallocate(rmatPositions%p_IposA)
    deallocate(rmatPositions%p_IposB)
    deallocate(rmatPositions%p_IposC)
    
    call lsysbl_releaseVector (rvectorTemp)
    call lsysbl_releaseVector (rrhsBlock)
    call lsysbl_releaseMatrix (RmatrixLocalBlock(1))

    call lsyssc_releaseVector (rrhsLocal)

    call storage_free(h_ImatrixMapping)
    call lsyssc_releaseMatrix (rmatrixLocal)
    
    deallocate(IchunkDofs)
    
    do j=1,nblockSize
      call lsysbl_releaseVector (RsolutionGlobal(j))
      call lsysbl_releaseVector (RrhsGlobal(j))
      do i=-1,1
        call lsysbl_releaseMatrix(RsystemMatrix(i,j))
      end do
    end do
    deallocate(RmatrixRows)
    deallocate(Rvectors)
    deallocate(RsystemMatrix)
    deallocate(RsolutionGlobal)
    deallocate(RrhsGlobal)
    
  contains
  
    subroutine buildMatrices (rproblem,rspaceTimeMatrix,&
      ichunkPos,nblockSize,RsystemMatrix,RmatrixRows)
    
    ! Builds the matrices in timestep ichunkPos..ichunkPos+nblockSize-1 into
    ! the RsystemMatrix and RmatrixRows arrays.
    
    ! The main problem structure
    type(t_problem), intent(INOUT) :: rproblem
    
    ! The underlying space time matrix.
    type(t_ccoptSpaceTimeMatrix), intent(IN) :: rspaceTimeMatrix
    
    ! Chunk position
    integer, intent(IN) :: ichunkPos
    
    ! Size of a chunk
    integer, intent(IN) :: nblockSize
    
    ! Destination array for the matrices
    type(t_matrixBlock), dimension(-1:,:), intent(INOUT) :: RsystemMatrix
    
    ! Destination array for matrix pointers.
    ! Variable is declared with INTENT(OUT), so all pointers are
    ! nullified automatically.
    type(t_matrixRow), dimension(:), intent(OUT) :: RmatrixRows
    
      ! local variables
      integer :: i,j,k,ichunk,NEQtime,ichunkrel
      real(DP) :: dtheta,dtime
      type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
      type(t_vectorBlock), dimension(3) :: rtimeVector
      type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr

      ! Create temp vectors for evaluating the nonlinearity
      call lsysbl_createVecBlockIndMat(RsystemMatrix(0,1),rtimeVector(1))
      call lsysbl_createVecBlockIndMat(RsystemMatrix(0,1),rtimeVector(2))
      call lsysbl_createVecBlockIndMat(RsystemMatrix(0,1),rtimeVector(3))
      
      ! Fetch some information
      NEQtime = rmatrix%p_rspaceTimeDiscr%NEQtime
      dtheta = rproblem%rtimedependence%dtimeStepTheta
      
      ! Basic initialisation of rnonlinearSpatialMatrix with the pointers to the
      ! matrices / discretisation structures on the current level.
      
      call smva_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
          p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
          p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
      
      ! Fetch evaluation vectors for the nonlinearity
      if (ichunkPos .gt. 1) &
        call sptivec_getTimestepData (rspaceTimeMatrix%p_rsolution, ichunkPos-1, rtimeVector(1))
        
      call sptivec_getTimestepData (rspaceTimeMatrix%p_rsolution, ichunkPos, rtimeVector(2))
      
      ! When retrieving matrix pointers, don't check the scaling factor!
      ! We need the 'maximum possible' layout of the matrix and so we get
      ! all pointers we are able to get.
      
      do ichunk = ichunkPos, ichunkPos+nblockSize-1
      
        ! Current point in time
        dtime = &
            rproblem%rtimedependence%dtimeInit + (ichunk-1)*p_rspaceTimeDiscr%rtimeDiscr%dtstep

        ! -----
        ! Discretise the boundary conditions at the new point in time --
        ! if the boundary conditions are nonconstant in time!
!        if (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .ne. 0) then
!          call cc_updateDiscreteBC (rproblem,dtime)
!        end if

        ! Relative position of the chunk item
        ichunkrel = ichunk-ichunkPos+1
      
        ! Get the evaluation point of the nonlinearity in the next timestep
        if (ichunk .lt. NEQtime) then
          call sptivec_getTimestepData (rspaceTimeMatrix%p_rsolution, ichunk+1, rtimeVector(3))
        end if
        
        ! Build 'left' matrix (if it exists):
        if (ichunk .gt. 1) then
        
          ! Set up the matrix weights
          call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
            ichunk-1,-1,rnonlinearSpatialMatrix)
          
          ! Set up the matrix
          call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
              RsystemMatrix(-1,ichunkrel),rnonlinearSpatialMatrix,&
              rtimeVector(1),rtimeVector(2),rtimeVector(3))
        
          ! Include the boundary conditions into that matrix.
          RsystemMatrix(-1,ichunkrel)%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
! not yet finished
!          call matfil_discreteBC (RsystemMatrix(-1,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC)
!          call matfil_discreteFBC (RsystemMatrix(-1,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC)

          ! Get the pointers to the submatrices.
          if (lsysbl_isSubmatrixPresent (RsystemMatrix(-1,ichunkrel),1,1)) &
            call lsyssc_getbase_double (RsystemMatrix(-1,ichunkrel)%RmatrixBlock(1,1),&
                RmatrixRows(ichunkrel)%p_Dp11)
          
          if (lsysbl_isSubmatrixPresent (RsystemMatrix(-1,ichunkrel),2,1)) &
            call lsyssc_getbase_double (RsystemMatrix(-1,ichunkrel)%RmatrixBlock(2,1),&
                RmatrixRows(ichunkrel)%p_Dp21)
          
          if (lsysbl_isSubmatrixPresent (RsystemMatrix(-1,ichunkrel),1,2)) &
            call lsyssc_getbase_double (RsystemMatrix(-1,ichunkrel)%RmatrixBlock(1,2),&
                RmatrixRows(ichunkrel)%p_Dp12)
          
          if (lsysbl_isSubmatrixPresent (RsystemMatrix(-1,ichunkrel),2,2)) &
            call lsyssc_getbase_double (RsystemMatrix(-1,ichunkrel)%RmatrixBlock(2,2),&
                RmatrixRows(ichunkrel)%p_Dp22)
        
        end if
        
        ! Build the 'diagonal' matrix
        !
        ! Set up the matrix weights
        call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ichunk-1,0,rnonlinearSpatialMatrix)
        
        ! Set up the matrix
        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            RsystemMatrix(0,ichunkrel),rnonlinearSpatialMatrix,&
            rtimeVector(1),rtimeVector(2),rtimeVector(3))
        
        ! Include the boundary conditions into that matrix.
! not yet finished
!        call matfil_discreteBC (RsystemMatrix(0,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC)
!        call matfil_discreteFBC (RsystemMatrix(0,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC)

        ! Get matrix pointers
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),1,1)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(1,1),&
              RmatrixRows(ichunkrel)%p_Da11)
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),2,1)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(2,1),&
              RmatrixRows(ichunkrel)%p_Da21)

        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),3,1)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(3,1),&
              RmatrixRows(ichunkrel)%p_Da31)

        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),4,1)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(4,1),&
              RmatrixRows(ichunkrel)%p_Da41)

        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),5,1)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(5,1),&
              RmatrixRows(ichunkrel)%p_Da51)
              
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),1,2)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(1,2),&
              RmatrixRows(ichunkrel)%p_Da12)
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),2,2)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(2,2),&
              RmatrixRows(ichunkrel)%p_Da22)

        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),3,2)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(3,2),&
              RmatrixRows(ichunkrel)%p_Da32)

        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),4,2)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(4,2),&
              RmatrixRows(ichunkrel)%p_Da42)

        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),5,2)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(5,2),&
              RmatrixRows(ichunkrel)%p_Da52)
        
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),1,3)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(1,3),&
              RmatrixRows(ichunkrel)%p_Da13)
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),2,3)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(2,3),&
              RmatrixRows(ichunkrel)%p_Da23)

        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),3,3)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(3,3),&
              RmatrixRows(ichunkrel)%p_Da33)
        
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),1,4)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(1,4),&
              RmatrixRows(ichunkrel)%p_Da14)
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),2,4)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(2,4),&
              RmatrixRows(ichunkrel)%p_Da24)
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),4,4)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(4,4),&
              RmatrixRows(ichunkrel)%p_Da44)
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),5,4)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(5,4),&
              RmatrixRows(ichunkrel)%p_Da54)

        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),6,4)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(6,4),&
              RmatrixRows(ichunkrel)%p_Da64)


        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),1,5)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(1,5),&
              RmatrixRows(ichunkrel)%p_Da15)
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),2,5)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(2,5),&
              RmatrixRows(ichunkrel)%p_Da25)
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),4,5)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(4,5),&
              RmatrixRows(ichunkrel)%p_Da45)
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),5,5)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(5,5),&
              RmatrixRows(ichunkrel)%p_Da55)

        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),6,5)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(6,5),&
              RmatrixRows(ichunkrel)%p_Da65)


        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),4,6)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(4,6),&
              RmatrixRows(ichunkrel)%p_Da46)
        
        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),5,6)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(5,6),&
              RmatrixRows(ichunkrel)%p_Da56)

        if (lsysbl_isSubmatrixPresent (RsystemMatrix(0,ichunkrel),6,6)) &
          call lsyssc_getbase_double (RsystemMatrix(0,ichunkrel)%RmatrixBlock(6,6),&
              RmatrixRows(ichunkrel)%p_Da66)
        
        ! Build 'right' matrix (if it exists):
        if (ichunk .lt. NEQtime) then
        
          ! Set up the matrix weights
          call stlin_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
            ichunk-1,1,rnonlinearSpatialMatrix)
          
          ! Set up the matrix
          call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
              RsystemMatrix(1,ichunkrel),rnonlinearSpatialMatrix,&
              rtimeVector(1),rtimeVector(2),rtimeVector(3))
        
          ! Include the boundary conditions into that matrix.
          RsystemMatrix(1,ichunkrel)%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
! not yet finished
!          call matfil_discreteBC (RsystemMatrix(1,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC)
!          call matfil_discreteFBC (RsystemMatrix(1,ichunkrel),p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC)

          ! Get the pointers to the submatrices
          if (lsysbl_isSubmatrixPresent (RsystemMatrix(1,ichunkrel),4,4)) &
            call lsyssc_getbase_double (RsystemMatrix(1,ichunkrel)%RmatrixBlock(4,4),&
                RmatrixRows(ichunkrel)%p_Dd44)
          
          if (lsysbl_isSubmatrixPresent (RsystemMatrix(1,ichunkrel),5,4)) &
            call lsyssc_getbase_double (RsystemMatrix(1,ichunkrel)%RmatrixBlock(5,4),&
                RmatrixRows(ichunkrel)%p_Dd54)
          
          if (lsysbl_isSubmatrixPresent (RsystemMatrix(1,ichunkrel),4,5)) &
            call lsyssc_getbase_double (RsystemMatrix(1,ichunkrel)%RmatrixBlock(4,5),&
                RmatrixRows(ichunkrel)%p_Dd45)
          
          if (lsysbl_isSubmatrixPresent (RsystemMatrix(1,ichunkrel),5,5)) &
            call lsyssc_getbase_double (RsystemMatrix(1,ichunkrel)%RmatrixBlock(5,5),&
                RmatrixRows(ichunkrel)%p_Dd55)
        end if
        

        if (ichunk .lt. NEQtime) then
          ! Cycle evaluation vectors: 1 <- 2 <- 3
          call lsysbl_copyVector (rtimeVector(2),rtimeVector(1))
          call lsysbl_copyVector (rtimeVector(3),rtimeVector(2))
        end if
        
        ! Copy all scaling factors
        do k=-1,1
          do j=1,6
            do i=1,6
              RmatrixRows(ichunkrel)%DscaleFactor(i,j,k) = &
                  RsystemMatrix(k,ichunkrel)%RmatrixBlock(i,j)%dscaleFactor
            end do
          end do
        end do
      
      end do
      
      ! Release memory
      call lsysbl_releaseVector (rtimeVector(3))
      call lsysbl_releaseVector (rtimeVector(2))
      call lsysbl_releaseVector (rtimeVector(1))
    
    end subroutine
    
    ! -----
    
    subroutine getMatrixPositions (IvelocityDofs,IpressureDofs,rmatrixPositions,&
        KcolA,KldA,KcolB,KldB,KcolC,KldC)
        
    ! Calculates the positions in the A, B and C matrices corresponding to the
    ! velocity and pressure DOF's.
    
    ! Array with velocity DOF's.
    integer, dimension(:), intent(IN) :: IvelocityDofs
    
    ! Array with pressure DOF's.
    integer, dimension(:), intent(IN) :: IpressureDofs
    
    ! A t_matrixPositions structure specifying the positions that
    ! are affected by the DOF's in the global matrices
    type(t_matrixPositions), intent(inout) :: rmatrixPositions

    ! Structure of FEM matrices for the velocity
    integer, dimension(:), intent(IN) :: KcolA
    integer, dimension(:), intent(IN) :: KldA

    ! Structure of FEM matrices for the gradient
    integer, dimension(:), intent(IN) :: KcolB
    integer, dimension(:), intent(IN) :: KldB
    
    ! Structure of FEM matrices for the pressure. Optional.
    integer, dimension(:), intent(IN), optional :: KcolC
    integer, dimension(:), intent(IN), optional :: KldC
    
      ! local variables
      integer :: i,j
      integer :: k
      integer :: idof
    
      ! Search for the matrix positions in the velocity and pressure matrix.
      ! IvelocityDofs specify the rows in the A- and B-matrix where we have
      ! to search.
      do i=1,size(IvelocityDofs)
        idof = IvelocityDofs(i)
        ! In the A-matrix, search for the velocity DOF's.
        do j=1,size(IvelocityDofs)
          do k=KldA(idof),KldA(idof+1)-1
            if (KcolA(k) .eq. IvelocityDofs(j)) then
              rmatrixPositions%p_IposA(i,j) = k
              exit
            end if
          end do
        end do
        
        ! In the B-matrices search for the pressure DOF's.
        do j=1,size(IpressureDofs)
          do k=KldB(idof),KldB(idof+1)-1
            if (KcolB(k) .eq. IpressureDofs(j)) then
              rmatrixPositions%p_IposB(i,j) = k
              exit
            end if
          end do
        end do
      end do
    
      if (present(KcolC)) then
        ! Search for the matrix positions pressure matrix.
        ! IpressureDofs specify the rows in the C-matrix where we have
        ! to search.
        do i=1,size(IpressureDofs)
          idof = IpressureDofs(i)
          ! In the C-matrix, search for the pressure DOF's.
          do j=1,size(IpressureDofs)
            do k=KldC(idof),KldC(idof+1)-1
              if (KcolC(k) .eq. IpressureDofs(j)) then
                rmatrixPositions%p_IposC(i,j) = k
                exit
              end if
            end do
          end do
        end do
    
      end if
    
    end subroutine
    
    ! -----
    
    subroutine getLocalMatrices (Rmatrixrows,IvelocityDofs,IpressureDofs,&
        ImatrixPositions,rmatrixPositions,Da,Kcol,Kld,KcolA,KldA,KcolB,KldB,KcolC,KldC)
    
    ! Extracts the entries of all matrices and creates a local matrix of
    ! a chunk.
    
    ! The matrices all the matrix rows in the chunk.
    type(t_matrixRow), dimension(:), intent(IN) :: Rmatrixrows
    
    ! The velocity DOF's under consideration
    integer, dimension(:), intent(IN) :: IvelocityDofs
    
    ! The pressure DOF's under consideration
    integer, dimension(:), intent(IN) :: IpressureDofs
    
    ! An integer array that specifies the start positions of the matrix columns
    ! in the local block matrix.
    integer(I32), dimension(:,:), intent(IN) :: ImatrixPositions
    
    ! A t_matrixPositions structure specifying the positions that
    ! are affected by the DOF's in the global matrices
    type(t_matrixPositions), intent(IN) :: rmatrixPositions
    
    ! Column/row structure of the mass/Laplace matrices
    integer, dimension(:), intent(IN) :: KcolA
    integer, dimension(:), intent(IN) :: KldA

    ! Column/row structure of the B matrices
    integer, dimension(:), intent(IN) :: KcolB
    integer, dimension(:), intent(IN) :: KldB

    ! Column/row structure of the C matrices. Optional.
    ! If not present, there is no C-matrix.
    integer, dimension(:), intent(IN), optional :: KcolC
    integer, dimension(:), intent(IN), optional :: KldC
    
    ! The local matrix that is to be filled with data.
    real(DP), dimension(:), intent(OUT) :: Da
    integer, dimension(:), intent(IN) :: Kcol
    integer, dimension(:), intent(IN) :: Kld
    
      ! Local variables
      integer :: irow,icolumn,ioffsetX,ioffsetY,irowOffsetY
    
      ! At first, initialise the output with zero.
      Da(:) = 0.0_DP
      
      ! Loop over all matrix rows in the chunk
      do irow = 1,size(Rmatrixrows)
      
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
        if (irow .ne. 1) then
         
          ! Extract the mass matrices coupling the primal velocity.
          call getLocalMatrix (Rmatrixrows(irow)%p_Dp11,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(1,1,-1),&
              ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,1))
          
          call getLocalMatrix (Rmatrixrows(irow)%p_Dp21,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(2,1,-1),&
              ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,1))

          call getLocalMatrix (Rmatrixrows(irow)%p_Dp12,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(1,2,-1),&
              ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,2))
          
          call getLocalMatrix (Rmatrixrows(irow)%p_Dp22,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(2,2,-1),&
              ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,2))
          
        end if
        
        ! Primal equation.
        call getLocalMatrix (Rmatrixrows(irow)%p_Da11,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(1,1,0),&
            ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,7))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da21,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(2,1,0),&
            ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,7))

        call getLocalMatrix (Rmatrixrows(irow)%p_Da31,Da,Kld,IpressureDofs,&
            rmatrixPositions%p_IposB,.true.,Rmatrixrows(irow)%DscaleFactor(3,1,0),&
            ioffsetY+iposLocalPrimalPressure,ImatrixPositions(irowOffsetY+3,7))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da41,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(4,1,0),&
            ioffsetY+iposLocalDualVelocityX,ImatrixPositions(irowOffsetY+4,7))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da51,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(5,1,0),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,7))


        call getLocalMatrix (Rmatrixrows(irow)%p_Da12,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(1,2,0),&
            ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,8))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da22,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(2,2,0),&
            ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,8))

        call getLocalMatrix (Rmatrixrows(irow)%p_Da32,Da,Kld,IpressureDofs,&
            rmatrixPositions%p_IposB,.true.,Rmatrixrows(irow)%DscaleFactor(3,2,0),&
            ioffsetY+iposLocalPrimalPressure,ImatrixPositions(irowOffsetY+3,8))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da42,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(4,2,0),&
            ioffsetY+iposLocalDualVelocityX,ImatrixPositions(irowOffsetY+4,8))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da52,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(5,2,0),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,8))


        call getLocalMatrix (Rmatrixrows(irow)%p_Da13,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposB,.false.,Rmatrixrows(irow)%DscaleFactor(1,3,0),&
            ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,9))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da23,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposB,.false.,Rmatrixrows(irow)%DscaleFactor(2,3,0),&
            ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,9))

        if (present(KcolC)) &
          call getLocalMatrix (Rmatrixrows(irow)%p_Da33,Da,Kld,IpressureDofs,&
              rmatrixPositions%p_IposC,.false.,Rmatrixrows(irow)%DscaleFactor(3,3,0),&
              ioffsetY+iposLocalPrimalPressure,ImatrixPositions(irowOffsetY+3,9))
            
            
        ! Dual equation
        !
        call getLocalMatrix (Rmatrixrows(irow)%p_Da14,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(1,4,0),&
            ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,10))
        
        ! Does not exist
        !CALL getLocalMatrix (Rmatrixrows(irow)%p_Da24,Da,Kld,IvelocityDofs,&
        !    rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(2,4,0),&
        !    ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,10))

        call getLocalMatrix (Rmatrixrows(irow)%p_Da44,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(4,4,0),&
            ioffsetY+iposLocalDualVelocityX,ImatrixPositions(irowOffsetY+4,10))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da54,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(5,4,0),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,10))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da64,Da,Kld,IpressureDofs,&
            rmatrixPositions%p_IposB,.true.,Rmatrixrows(irow)%DscaleFactor(6,4,0),&
            ioffsetY+iposLocalDualPressure,ImatrixPositions(irowOffsetY+6,10))
        

        ! Does not exist
        !CALL getLocalMatrix (Rmatrixrows(irow)%p_Da15,Da,Kld,IvelocityDofs,&
        !    rmatrixPositions%p_IposA,.FALSE.,Rmatrixrows(irow)%DscaleFactor(1,5,0),&
        !    ioffsetY+iposLocalPrimalVelocityX,ImatrixPositions(irowOffsetY+1,11))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da25,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(2,5,0),&
            ioffsetY+iposLocalPrimalVelocityY,ImatrixPositions(irowOffsetY+2,11))

        call getLocalMatrix (Rmatrixrows(irow)%p_Da45,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(4,5,0),&
            ioffsetY+iposLocalDualVelocityX, ImatrixPositions(irowOffsetY+4,11))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da55,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(5,5,0),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,11))

        call getLocalMatrix (Rmatrixrows(irow)%p_Da65,Da,Kld,IpressureDofs,&
            rmatrixPositions%p_IposB,.true.,Rmatrixrows(irow)%DscaleFactor(6,5,0),&
            ioffsetY+iposLocalDualPressure,ImatrixPositions(irowOffsetY+6,11))


        call getLocalMatrix (Rmatrixrows(irow)%p_Da46,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposB,.false.,Rmatrixrows(irow)%DscaleFactor(4,6,0),&
            ioffsetY+iposLocalDualVelocityX,ImatrixPositions(irowOffsetY+4,12))
        
        call getLocalMatrix (Rmatrixrows(irow)%p_Da56,Da,Kld,IvelocityDofs,&
            rmatrixPositions%p_IposB,.false.,Rmatrixrows(irow)%DscaleFactor(5,6,0),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,12))

        if (present(KcolC)) &
          call getLocalMatrix (Rmatrixrows(irow)%p_Da66,Da,Kld,IpressureDofs,&
              rmatrixPositions%p_IposC,.false.,Rmatrixrows(irow)%DscaleFactor(6,6,0),&
              ioffsetY+iposLocalDualPressure,ImatrixPositions(irowOffsetY+6,12))


        if (irow .lt. size(Rmatrixrows)) then
          ! Extract the mass matrices coupling the dual velocity.
          call getLocalMatrix (Rmatrixrows(irow)%p_Dd44,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(4,4,1),&
            ioffsetY+iposLocalDualVelocityX,ImatrixPositions(irowOffsetY+4,16))
          
          call getLocalMatrix (Rmatrixrows(irow)%p_Dd54,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(5,4,1),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,16))

          call getLocalMatrix (Rmatrixrows(irow)%p_Dd45,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(4,5,1),&
            ioffsetY+iposLocalDualVelocityX, ImatrixPositions(irowOffsetY+4,17))
          
          call getLocalMatrix (Rmatrixrows(irow)%p_Dd55,Da,Kld,IvelocityDofs,&
              rmatrixPositions%p_IposA,.false.,Rmatrixrows(irow)%DscaleFactor(5,5,1),&
            ioffsetY+iposLocalDualVelocityY,ImatrixPositions(irowOffsetY+5,17))
        end if

      end do

    end subroutine

    ! -----

    subroutine getLocalMatrix (p_Da,Db,KldB,Idofs,&
        ImatrixPositions,btransposed,dscale,idestOffsetY,idestOffsetX)

    ! Pointer to the source matrix. May point to NULL() if there is none.
    real(DP), dimension(:), pointer :: p_Da

    ! Destination matrix
    real(DP), dimension(:), intent(OUT) :: Db
    integer, dimension(:), intent(IN) :: KldB
    
    ! Rows in the source matrix that should be extracted.
    integer, dimension(:), intent(IN) :: Idofs
    
    ! Array with matrix positions. The rows in this array
    ! define for every entry in Idof the positions in the full matrix of
    ! the entries that have to be extracted and written to the local matrix.
    integer, dimension(:,:), intent(IN) :: ImatrixPositions
    
    ! Scaling factor for the matrix
    real(DP), intent(IN) :: dscale
    
    ! Column offset in the destination matrix. Entries are written to
    ! position idestOffsetX..idestOffsetX+UBOUND(ImatrixPositions,2)-1 into
    ! line idestOffsetY..idestOffsetY+UBOUND(ImatrixPositions,1)-1 in
    ! the compressed destination matrix.
    ! That way, idestOffsetX/Y specifies the matrix column/row where to write
    ! the data to.
    integer, intent(IN) :: idestOffsetX,idestOffsetY
    
    ! If set to TRUE, the routine assumes that ImatrixPosition specifies
    ! the entry positions in the transposed matrix. Used for handling B^T.
    logical, intent(IN) :: btransposed
    
      ! local variables
      integer :: i,j
      integer :: k,irowPosB
      
      if (.not. associated(p_Da)) return
      
      if (.not. btransposed) then
        ! Loop over the DOF's to extract
        do i=1,ubound(ImatrixPositions,1)
          irowPosB = KldB(idestOffsetY-1+i)+idestOffsetX-1
          do j=1,ubound(ImatrixPositions,2)
            
            ! Get the matrix position where to find our entry.
            k = ImatrixPositions(i,j)
            
            ! Get the matrix entry and write it to the destination matrix.
            Db(irowPosB+j-1) = dscale*p_Da(k)
          end do
        end do
      else
        ! Transposed B-matrix.
        ! Loop over the DOF's to extract
        do j=1,ubound(ImatrixPositions,2)
          irowPosB = KldB(idestOffsetY-1+j)+idestOffsetX-1
          do i=1,ubound(ImatrixPositions,1)
            
            ! Get the matrix position where to find our entry.
            k = ImatrixPositions(i,j)
            
            ! Get the matrix entry and write it to the destination matrix.
            Db(irowPosB+i-1) = dscale*p_Da(k)
          end do
        end do
      end if
    
    end subroutine

    ! -----
    
    subroutine getLocalDefect (Dd,rmatrixrow,Rvectors,IvelocityDofs,IpressureDofs,&
        IchunkDofs,rmatrixPositions,bfirstTimestep,blastTimestep,&
        KcolA,KldA,KcolB,KldB,KcolC,KldC)
    
    ! Calculates the local defect in one row of a chunk.
    
    ! The matrices in the current row of the chunk.
    type(t_matrixRow), intent(IN) :: rmatrixrow
    
    ! Pointers to the RHS and solution vectors in the current row of a chunk.
    ! the array must always have dimension 3, where Rvectors(2) identifies
    ! the RHS/solution that corresponds to the diagonal matrix.
    type(t_chunkVector), dimension(3), intent(IN) :: Rvectors
    
    ! The velocity DOF's under consideration
    integer, dimension(:), intent(IN) :: IvelocityDofs
    
    ! The pressure DOF's under consideration
    integer, dimension(:), intent(IN) :: IpressureDofs
    
    ! A list of all DOF's in the current chunk on the current element
    integer, dimension(:), intent(IN) :: IchunkDofs
    
    ! A t_matrixPositions structure specifying the positions that
    ! are affected by the DOF's in the global matrices
    type(t_matrixPositions), intent(IN) :: rmatrixPositions
    
    ! Flag that specifies if the current timestep is the first timestep
    ! in the chunk.
    logical, intent(IN) :: bfirstTimestep
    
    ! Flag that specifies if the current timestep is the last timestep
    ! in the chunk.
    logical, intent(IN) :: blastTimestep
    
    ! Column/row structure of the mass/Laplace matrices
    integer, dimension(:), intent(IN) :: KcolA
    integer, dimension(:), intent(IN) :: KldA

    ! Column/row structure of the B matrices
    integer, dimension(:), intent(IN) :: KcolB
    integer, dimension(:), intent(IN) :: KldB

    ! Column/row structure of the C matrices. Optional.
    ! If not present, there is no C-matrix.
    integer, dimension(:), intent(IN), optional :: KcolC
    integer, dimension(:), intent(IN), optional :: KldC
    
    ! The local defect in one row of a chunk
    real(DP), dimension(:), intent(OUT), target :: Dd

      ! local variables
      integer :: k
      real(DP), dimension(:), pointer :: p_Dd1
      real(DP), dimension(:), pointer :: p_Dd2
      real(DP), dimension(:), pointer :: p_Dd3
      real(DP), dimension(:), pointer :: p_Dd4
      real(DP), dimension(:), pointer :: p_Dd5
      real(DP), dimension(:), pointer :: p_Dd6

      ! Get the local RHS
      do k=1,size(Dd)
        Dd(k) = Rvectors(2)%p_Db(IchunkDofs(k))
      end do
      
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
      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,1,0),&
          rmatrixrow%p_Da11,KcolA,KldA,p_Dd1,Rvectors(2)%p_Dx1)

      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,1,0),&
          rmatrixrow%p_Da21,KcolA,KldA,p_Dd2,Rvectors(2)%p_Dx1)
          
      call localMatrixVectorBT (IvelocityDofs,IpressureDofs,&
          rmatrixrow%DscaleFactor(3,1,0),rmatrixrow%p_Da31,&
          KcolB,rmatrixPositions%p_IposB,p_Dd3,Rvectors(2)%p_Dx1)



      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,1,0),&
          rmatrixrow%p_Da41,KcolA,KldA,p_Dd4,Rvectors(2)%p_Dx1)

      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,1,0),&
          rmatrixrow%p_Da51,KcolA,KldA,p_Dd5,Rvectors(2)%p_Dx1)



      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,2,0),&
          rmatrixrow%p_Da12,KcolA,KldA,p_Dd1,Rvectors(2)%p_Dx2)

      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,2,0),&
          rmatrixrow%p_Da22,KcolA,KldA,p_Dd2,Rvectors(2)%p_Dx2)

      call localMatrixVectorBT (IvelocityDofs,IpressureDofs,&
          rmatrixrow%DscaleFactor(3,2,0),rmatrixrow%p_Da32,&
          KcolB,rmatrixPositions%p_IposB,p_Dd3,Rvectors(2)%p_Dx2)


      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,2,0),&
          rmatrixrow%p_Da42,KcolA,KldA,p_Dd4,Rvectors(2)%p_Dx2)

      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,2,0),&
          rmatrixrow%p_Da52,KcolA,KldA,p_Dd5,Rvectors(2)%p_Dx2)


      call localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(1,3,0),&
          rmatrixrow%p_Da13,KcolB,KldB,p_Dd1,Rvectors(2)%p_Dx3)

      call localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(2,3,0),&
          rmatrixrow%p_Da23,KcolB,KldB,p_Dd2,Rvectors(2)%p_Dx3)

      if (present(KcolC)) &
        call localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(3,3,0),&
          rmatrixrow%p_Da33,KcolC,KldC,p_Dd3,Rvectors(2)%p_Dx3)

      ! Dual solution
      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,4,0),&
          rmatrixrow%p_Da14,KcolA,KldA,p_Dd1,Rvectors(2)%p_Dx4)

      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,4,0),&
          rmatrixrow%p_Da24,KcolA,KldA,p_Dd2,Rvectors(2)%p_Dx4)



      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,4,0),&
          rmatrixrow%p_Da44,KcolA,KldA,p_Dd4,Rvectors(2)%p_Dx4)

      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,4,0),&
          rmatrixrow%p_Da54,KcolA,KldA,p_Dd5,Rvectors(2)%p_Dx4)

      call localMatrixVectorBT (IvelocityDofs,IpressureDofs,&
          rmatrixrow%DscaleFactor(6,4,0),&
          rmatrixrow%p_Da64,KcolB,rmatrixPositions%p_IposB,p_Dd6,Rvectors(2)%p_Dx4)


      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,5,0),&
          rmatrixrow%p_Da15,KcolA,KldA,p_Dd1,Rvectors(2)%p_Dx5)

      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,5,0),&
          rmatrixrow%p_Da25,KcolA,KldA,p_Dd2,Rvectors(2)%p_Dx5)



      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,5,0),&
          rmatrixrow%p_Da45,KcolA,KldA,p_Dd4,Rvectors(2)%p_Dx5)

      call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,5,0),&
          rmatrixrow%p_Da55,KcolA,KldA,p_Dd5,Rvectors(2)%p_Dx5)

      call localMatrixVectorBT (IvelocityDofs,IpressureDofs,&
          rmatrixrow%DscaleFactor(6,5,0),&
          rmatrixrow%p_Da65,KcolB,rmatrixPositions%p_IposB,p_Dd6,Rvectors(2)%p_Dx5)


      call localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(4,6,0),&
          rmatrixrow%p_Da46,KcolB,KldB,p_Dd4,Rvectors(2)%p_Dx6)

      call localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(5,6,0),&
          rmatrixrow%p_Da56,KcolB,KldB,p_Dd5,Rvectors(2)%p_Dx6)

      if (present(KcolC)) &
        call localMatrixVector(IpressureDofs,rmatrixrow%DscaleFactor(6,6,0),&
          rmatrixrow%p_Da66,KcolC,KldC,p_Dd6,Rvectors(2)%p_Dx6)
      
      ! Mass matrix coupling the primal velocity
      if (.not. bfirstTimestep) then
        call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,1,-1),&
          rmatrixrow%p_Dp11,KcolA,KldA,p_Dd1,Rvectors(1)%p_Dx1)

        call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,1,-1),&
          rmatrixrow%p_Dp21,KcolA,KldA,p_Dd2,Rvectors(1)%p_Dx1)

        call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(1,2,-1),&
          rmatrixrow%p_Dp12,KcolA,KldA,p_Dd1,Rvectors(1)%p_Dx2)

        call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(2,2,-1),&
          rmatrixrow%p_Dp22,KcolA,KldA,p_Dd2,Rvectors(1)%p_Dx2)
      end if

      ! Mass matrix coupling the dual velocity
      if (.not. blastTimestep) then
        call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,4,1),&
          rmatrixrow%p_Dd44,KcolA,KldA,p_Dd4,Rvectors(3)%p_Dx4)

        call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,4,1),&
          rmatrixrow%p_Dd54,KcolA,KldA,p_Dd5,Rvectors(3)%p_Dx4)

        call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(4,5,1),&
          rmatrixrow%p_Dd45,KcolA,KldA,p_Dd4,Rvectors(3)%p_Dx5)

        call localMatrixVector(IvelocityDofs,rmatrixrow%DscaleFactor(5,5,1),&
          rmatrixrow%p_Dd55,KcolA,KldA,p_Dd5,Rvectors(3)%p_Dx5)
      end if
      
    end subroutine


    subroutine addLocalCorrection (Dd,rvectors,IchunkDofs,domega)
    
    ! Adds a local correction vector to the global solution
    
    ! A list of all DOF's in the current chunk on the current element
    integer, dimension(:), intent(IN) :: IchunkDofs
    
    ! The local defect in one row of a chunk
    real(DP), dimension(:), intent(IN), target :: Dd

    ! Damping parameter; standard = 1.0
    real(DP), intent(IN) :: domega
    
    ! Pointers to the RHS and solution vectors in the current row of a chunk.
    ! The solution in here is updated.
    type(t_chunkVector), intent(INOUT) :: rvectors
    
      ! local variables
      integer :: k

      ! Get the local RHS
      do k=1,size(Dd)
        rvectors%p_Dx(IchunkDofs(k)) = rvectors%p_Dx(IchunkDofs(k)) + domega*Dd(k)
      end do
      
    end subroutine
    
    subroutine localMatrixVector (Idofs,dscale,p_Da,Kcol,Kld,Db,Dx)
    
    ! Performs a local matrix vector multiplication
    !   b = b-Ax
    ! for the rows defined in the array Idofs.
    
    ! The rows in the matrix that correspond to the elements in Db.
    integer, dimension(:), intent(IN) :: Idofs
    
    ! Scaling factor for the matrix
    real(DP), intent(IN) :: dscale
    
    ! Matrix data, structure 9. the data array is specified as a pointer;
    ! if set to NULL(), there is no matrix, so no matrix vector multiplication
    ! will be done.
    real(DP), dimension(:), pointer :: p_Da
    integer, dimension(:), intent(IN) :: Kcol
    integer, dimension(:), intent(IN) :: Kld
    
    ! Solution vector; is multiplied by the matrix and subtracted from b.
    real(DP), dimension(:), intent(IN) :: Dx
    
    ! RHS vector. The size must be the same as Idofs.
    ! Is overwritten by the defect.
    real(DP), dimension(:), intent(INOUT) :: Db
    
      integer :: i,j
      
      if (associated(p_Da) .and. (dscale .ne. 0.0_DP)) then
      
        ! Loop through all DOF's. These defines the lines in Db where we have
        ! to do b-Ax.
        do i=1,size(Idofs)
          ! Subtract row IDofs(i) * Dx from Db.
          do j=Kld(Idofs(i)),Kld(Idofs(i)+1)-1
            Db(i) = Db(i) - dscale*p_Da(j)*Dx(Kcol(j))
          end do
        end do
        
      end if
    
    end subroutine

    subroutine localMatrixVectorBT (IvelocityDofs,IpressureDofs,dscale,&
        p_DB,KcolB,ImatrixPos,Db,Dx)
    
    ! Performs a local matrix vector multiplication
    !   b = b - B^T x
    ! for the rows defined in the array Idofs.
    
    ! Velocity DOF's under consideration
    integer, dimension(:), intent(IN) :: IvelocityDofs

    ! Pressure DOF's under consideration
    integer, dimension(:), intent(IN) :: IpressureDofs
    
    ! Scaling factor for the matrix
    real(DP), intent(IN) :: dscale
    
    ! Matrix data, structure 9. the data array is specified as a pointer;
    ! if set to NULL(), there is no matrix, so no matrix vector multiplication
    ! will be done.
    real(DP), dimension(:), pointer :: p_DB
    integer, dimension(:), intent(IN) :: KcolB
    
    ! Index array specifying the entries (rows) in the B-matrix which are affected
    ! by the DOF's in IvelocityDofs. DIMENSION(#velocity DOF's,#pressure DOF's).
    integer, dimension(:,:), intent(IN) :: ImatrixPos
    
    ! Solution vector; is multiplied by the matrix and subtracted from b.
    real(DP), dimension(:), intent(IN) :: Dx
    
    ! RHS vector. The size must be the same as Idofs.
    ! Is overwritten by the defect.
    real(DP), dimension(:), intent(INOUT) :: Db
    
      integer :: i,j
      integer :: k
      
      if (associated(p_DB) .and. (dscale .ne. 0.0_DP)) then
      
        ! Loop through all pressure DOF's. These defines the lines in DB where we have
        ! to do b-Ax.
        do i = 1,size(IpressureDofs)
          
          ! Every pressure DOF gives us a line in the B^T matrix tat we have
          ! to process. In that line, we have to loop through all velocity DOF's.
          do j = 1,size(IvelocityDofs)
          
            ! Find the matrix position in the B matrix specifying that is affected
            ! by the current velocity DOF. This defines the entry that we have to
            ! take as the entry of B^T.
            k = ImatrixPos (j,i)
            
            ! Multiply the B^T entry with the velocity DOF and subtract.
            Db(i) = Db(i) - dscale*p_DB(k)*Dx(IvelocityDofs(j))
          
          end do
         
        end do
        
      end if
      
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
    
    end subroutine

  end subroutine

end module
