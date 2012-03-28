!##############################################################################
!# ****************************************************************************
!# <name> globalsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to assembla a global 1x1 system matrix from
!# a block matrix containing multiple submatrices.
!#
!# The following routines can be found here:
!#
!# 1.) glsys_assembleGlobal
!#     -> Assemble a global matrix.
!# </purpose>
!##############################################################################

module globalsystem

!$use omp_lib
  use fsystem
  use genoutput
  use storage
  use linearsystemscalar
  use linearsystemblock

  implicit none

  private

  public :: glsys_assembleGlobal

contains

  ! ***************************************************************************

!<subroutine>

  subroutine glsys_assembleGlobal (rsourceMatrix, rdestMatrix, &
                                   bstructure, bcontent, &
                                   cmatrixFormat, cdataType)

!<description>
  ! This routine assembles a 1x1 block matrix rdestMatrix from a
  ! NxM block matrix rsourceMatrix.
!</description>

!<input>
  ! The source matrix.
  type(t_matrixBlock), intent(in) :: rsourceMatrix

  ! Whether to assemble a new structure in rdestMatrix.
  ! =TRUE (Standard): Release any old data/structure from rdestMatrix
  !   and create a new one.
  ! =FALSE          : rdestMatrix is assumed to be an existing matrix
  !   of the correct dimension/size which does not have to be rebuild.
  logical, intent(in) :: bstructure

  ! Whether to assemble the matrix content.
  ! =TRUE (Standard): The content of rsourceMatrix is build in rdestMatrix.
  ! =FALSE          : No matrix content is build up in rdestMatrix.
  logical, intent(in) :: bcontent

  ! OPTIONAL: Target format of the matrix rdestMatrix. Standard is Format 9.
  integer, intent(in), optional :: cmatrixFormat

  ! OPTIONAL: Data type for the entries of rdestMatrix.
  ! Standard is double precision.
  integer, intent(in), optional :: cdataType

!</input>

!<output>
  ! The destination matrix. If this matrix contains valid handles to
  ! allocated memory blocks of the correct size for structure/entries,
  ! this data is overwritten. If not, arrays are (re-)allocated in
  ! the correct size.
  type(t_matrixBlock), intent(inout) :: rdestMatrix
!</output>

!</subroutine>

    ! local variables
    integer :: cdataTypeLocal, cmatrixFormatLocal,i,j
    logical :: balloc
    type(t_matrixBlock) :: rlocalMatrix
    type(t_matrixScalar) :: rlocalMatrixScalar
    integer, dimension(:), allocatable :: Irows
    integer, dimension(:), allocatable :: Icolumns
    integer, dimension(:), pointer :: p_Kdiagonal,p_Kld
    integer, dimension(:), pointer :: p_Kcol
    integer :: isize

    ! Initialise values for data type and matrix format.
    cdataTypeLocal = ST_DOUBLE
    cmatrixFormatLocal = LSYSSC_MATRIX9

    ! If the matrix has already data, try to fetch the data type/
    ! matrix format from the existing matrix.
    if (lsysbl_isSubmatrixPresent (rdestMatrix,1,1)) then
      if (lsyssc_hasMatrixStructure(rdestMatrix%RmatrixBlock(1,1))) &
          cmatrixFormatLocal = rdestMatrix%RmatrixBlock(1,1)%cmatrixFormat
      if (lsyssc_hasMatrixContent(rdestMatrix%RmatrixBlock(1,1))) &
          cdataTypeLocal = rdestMatrix%RmatrixBlock(1,1)%cdataType
    end if

    if (present(cdataType)) cdataTypeLocal = cdataType
    if (present(cmatrixFormat)) cmatrixFormatLocal = cmatrixFormat

    ! Up to now, we do not support everything!
    ! Cancel if the input parameters want too much of us...
    if ((rsourceMatrix%nblocksPerCol .eq. 0) .or. &
        (rsourceMatrix%nblocksPerRow .eq. 0)) return

    if (cdataTypeLocal .ne. ST_DOUBLE) then
      call output_line('Only double precision destination matrix supported!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'glsys_assembleGlobal')
      call sys_halt()
    end if

    if (cmatrixFormatLocal .ne. LSYSSC_MATRIX9) then
      call output_line('Only format 9 destination matrix supported!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'glsys_assembleGlobal')
      call sys_halt()
    end if

    do j=1,rsourceMatrix%nblocksPerRow
      do i=1,rsourceMatrix%nblocksPerCol

        if (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) then

          if (rsourceMatrix%RmatrixBlock(i,j)%cdataType .ne. ST_DOUBLE) then
            call output_line('Only double precision source matrices supported!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'glsys_assembleGlobal')
            call sys_halt()
          end if

          if (rsourceMatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX9) then
            call output_line('Only format 9 source matrices supported!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'glsys_assembleGlobal')
            call sys_halt()
          end if

        end if

      end do
    end do

  ! Now start the actual assembly. What and how to assemble depend on
  ! the boolean input variables.
  if ((.not. bstructure) .and. (.not. bcontent)) return ! nothing to do.

  ! Problem: Some of the submatrices may be virtually transposed,
  ! but the algorithms here can only handle un-transposed matrices.
  ! So at first, we create an un-transposed global matrix.

  call lsysbl_duplicateMatrix (rsourceMatrix,rlocalMatrix, &
                               LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

  do j=1,rlocalMatrix%nblocksPerRow
    do i=1,rlocalMatrix%nblocksPerCol

      if (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) then
        ! Transpose the submatrix if necessary
        if (iand(rsourceMatrix%RmatrixBlock(i,j)%imatrixSpec, &
                 LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
          ! Untranspose the source-submatrix to a local matrix
          call lsyssc_transposeMatrix (rsourceMatrix%RmatrixBlock(i,j),&
                      rlocalMatrixScalar,LSYSSC_TR_VIRTUAL)

          ! Retranspose it - not vitually, but create a real transposed copy.
          if (bcontent) then
            ! Transpose everything
            call lsyssc_releaseMatrix(rlocalMatrix%RmatrixBlock(i,j))
            call lsyssc_transposeMatrix (rlocalMatrixScalar,&
                        rlocalMatrix%RmatrixBlock(i,j),LSYSSC_TR_ALL)
          else
            ! Transpose only the structure, ignore the content
            call lsyssc_releaseMatrix(rlocalMatrix%RmatrixBlock(i,j))
            call lsyssc_transposeMatrix (rlocalMatrixScalar,&
                        rlocalMatrix%RmatrixBlock(i,j),LSYSSC_TR_STRUCTURE)
          end if

        end if
      end if
    end do
  end do

  ! Ok, rlocalMatrix is now a transposed-free source matrix.
  allocate(Irows(max(rsourceMatrix%nblocksPerCol,1)+1))
  allocate(Icolumns(max(rsourceMatrix%nblocksPerRow,1)+1))

  if (bstructure .and. bcontent) then

    ! Initialise general data of the destination matrix.
    call glmatasm_initDestination (rlocalMatrix,rdestMatrix,&
                                   cmatrixFormatLocal,cdataTypeLocal)

    ! Assemble the global matrix - structure and entries.
    !
    ! What is the destination matrix structure?
    select case (cmatrixFormatLocal)
    case (LSYSSC_MATRIX9)

      ! Get the basic row/column indices of the destination matrix.
      call glmatasm_getOffsets (rlocalMatrix,Icolumns,Irows)

      ! Allocate a KLD in the destination matrix if we do not have a previous
      ! array in the correct size.
      balloc = .true.
      if ((.not. lsyssc_isMatrixStructureShared (rdestMatrix%RmatrixBlock(1,1))) .and. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Kld .ne. ST_NOHANDLE)) then
        call storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kld,isize)
        ! Matrix is for sure not transposed!
        if (isize .eq. rlocalMatrix%NEQ+1) then
          balloc = .false.
        else
          ! Release the previous memory before allocating a new block.
          call storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Kld)
        end if
      end if

      if (balloc) call storage_new ('glsys_assembleGlobal', 'Kld',  &
                                    rlocalMatrix%NEQ+1, &
                                    ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kld,&
                                    ST_NEWBLOCK_ZERO)

      ! Set up KLD and NA of the destination matrix
      call glmatasm_KLD (rlocalMatrix,rdestMatrix%RmatrixBlock(1,1),Irows)

      ! Allocate a KCOL in the destination matrix if we do not have a previous
      ! array in the correct size.
      balloc = .true.
      if ((.not. lsyssc_isMatrixStructureShared (rdestMatrix%RmatrixBlock(1,1))) .and. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Kcol .ne. ST_NOHANDLE)) then
        call storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kcol,isize)
        ! Matrix is for sure not transposed!
        if (isize .eq. rdestMatrix%RmatrixBlock(1,1)%NA) then
          balloc = .false.
        else
          ! Release the previous memory before allocating a new block.
          call storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Kcol)
        end if
      end if
      if (balloc) call storage_new ('glsys_assembleGlobal', 'Kcol', &
                                    rdestMatrix%RmatrixBlock(1,1)%NA, &
                                    ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kcol,&
                                    ST_NEWBLOCK_NOINIT)

      ! Allocate the data array if we do not have a previous
      ! array in the correct size.
      balloc = .true.
      if ((.not. lsyssc_isMatrixContentShared (rdestMatrix%RmatrixBlock(1,1))) .and. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Da .ne. ST_NOHANDLE)) then
        call storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Da,isize)
        if (isize .eq. rdestMatrix%RmatrixBlock(1,1)%NA) then
          balloc = .false.
        else
          ! Release the previous memory before allocating a new block.
          call storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Da)
        end if
      end if
      if (balloc) call storage_new ('glsys_assembleGlobal', 'Da', &
                                    rdestMatrix%RmatrixBlock(1,1)%NA, &
                                    cdataTypeLocal, &
                                    rdestMatrix%RmatrixBlock(1,1)%h_Da,&
                                    ST_NEWBLOCK_NOINIT)

      ! Set up KCOL and the matrix entries
      call glmatasm_KcolDa99dble (rlocalMatrix,rdestMatrix%RmatrixBlock(1,1),&
                                  Icolumns, Irows)

      ! Allocate a Kdiagonal in the destination matrix if we do not have a previous
      ! array in the correct size.
      balloc = .true.
      if ((.not. lsyssc_isMatrixStructureShared (rdestMatrix%RmatrixBlock(1,1))) .and. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal .ne. ST_NOHANDLE)) then
        call storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal,isize)
        ! Matrix is for sure not transposed!
        if (isize .eq. rlocalMatrix%NEQ) then
          balloc = .false.
        else
          ! Release the previous memory before allocating a new block.
          call storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal)
        end if
      end if

      ! Allocate a Kdiagonal in the destination matrix
      if (balloc) call storage_new (&
                        'glsys_assembleGlobal', 'Kdiagonal', rlocalMatrix%NEQ, &
                        ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal,&
                        ST_NEWBLOCK_NOINIT)

      ! Rebuild Kdiagonal
      call lsyssc_getbase_Kdiagonal(rdestMatrix%RmatrixBlock(1,1),p_Kdiagonal)
      call lsyssc_getbase_Kcol(rdestMatrix%RmatrixBlock(1,1),p_Kcol)
      call lsyssc_getbase_Kld(rdestMatrix%RmatrixBlock(1,1),p_Kld)
      call lsyssc_rebuildKdiagonal (p_Kcol, p_Kld, p_Kdiagonal, rdestMatrix%NEQ)

    end select

  else if (bstructure) then

    ! Assemble only the structure of the global matrix

    ! Initialise general data of the destination matrix.
    call glmatasm_initDestination (rlocalMatrix,rdestMatrix,&
                                   cmatrixFormatLocal,cdataTypeLocal)

    ! What is the destination matrix structure?
    select case (cmatrixFormatLocal)
    case (LSYSSC_MATRIX9)

      ! Get the basic row/column indices of the destination matrix.
      call glmatasm_getOffsets (rlocalMatrix,Icolumns,Irows)

      ! Allocate a KLD in the destination matrix if we do not have a previous
      ! array in the correct size.
      balloc = .true.
      if ((.not. lsyssc_isMatrixStructureShared (rdestMatrix%RmatrixBlock(1,1))) .and. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Kld .ne. ST_NOHANDLE)) then
        call storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kld,isize)
        ! Matrix is for sure not transposed!
        if (isize .eq. rlocalMatrix%NEQ+1) then
          balloc = .false.
        else
          ! Release the previous memory before allocating a new block.
          call storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Kld)
        end if
      end if

      if (balloc) call storage_new ('glsys_assembleGlobal', 'Kld',  &
                                    rlocalMatrix%NEQ+1, &
                                    ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kld,&
                                    ST_NEWBLOCK_ZERO)

      ! Set up KLD and NA of the destination matrix
      call glmatasm_KLD (rlocalMatrix,rdestMatrix%RmatrixBlock(1,1),Irows)

      ! Allocate a KCOL in the destination matrix if we do not have a previous
      ! array in the correct size.
      balloc = .true.
      if ((.not. lsyssc_isMatrixStructureShared (rdestMatrix%RmatrixBlock(1,1))) .and. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Kcol .ne. ST_NOHANDLE)) then
        call storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kcol,isize)
        ! Matrix is for sure not transposed!
        if (isize .eq. rdestMatrix%RmatrixBlock(1,1)%NA) then
          balloc = .false.
        else
          ! Release the previous memory before allocating a new block.
          call storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Kcol)
        end if
      end if
      if (balloc) call storage_new ('glsys_assembleGlobal', 'Kcol', &
                                    rdestMatrix%RmatrixBlock(1,1)%NA, &
                                    ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kcol,&
                                    ST_NEWBLOCK_NOINIT)

      ! Set up KCOL and the matrix entries
      call glmatasm_Kcol99dble (rlocalMatrix,rdestMatrix%RmatrixBlock(1,1),&
                                Icolumns, Irows)

      ! Allocate a Kdiagonal in the destination matrix
      call storage_new ('glsys_assembleGlobal', 'Kdiagonal', rlocalMatrix%NEQ, &
                        ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal,&
                        ST_NEWBLOCK_NOINIT)

      ! Rebuild Kdiagonal
      call lsyssc_getbase_Kdiagonal(rdestMatrix%RmatrixBlock(1,1),p_Kdiagonal)
      call lsyssc_getbase_Kcol(rdestMatrix%RmatrixBlock(1,1),p_Kcol)
      call lsyssc_getbase_Kld(rdestMatrix%RmatrixBlock(1,1),p_Kld)
      call lsyssc_rebuildKdiagonal (p_Kcol, p_Kld, p_Kdiagonal, rdestMatrix%NEQ)

    end select

  else

    ! Assemble only the entries of the global matrix; the structure
    ! is already there.
    ! Initialise general data of the destination matrix.
    call glmatasm_initDestination (rlocalMatrix,rdestMatrix,&
                                   cmatrixFormatLocal,cdataTypeLocal)

    ! What is the destination matrix structure?
    select case (cmatrixFormatLocal)
    case (LSYSSC_MATRIX9)

      ! Get the basic row/column indices of the destination matrix.
      call glmatasm_getOffsets (rlocalMatrix,Icolumns,Irows)

      ! Kcol/Kld are assumed to be ok.
      !
      ! Allocate the data array if we do not have a previous
      ! array in the correct size.
      balloc = .true.
      if ((.not. lsyssc_isMatrixContentShared (rdestMatrix%RmatrixBlock(1,1))) .and. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Da .ne. ST_NOHANDLE)) then
        call storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Da,isize)
        if (isize .eq. rdestMatrix%RmatrixBlock(1,1)%NA) then
          balloc = .false.
        else
          ! Release the previous memory before allocating a new block.
          call storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Da)
        end if
      end if
      if (balloc) call storage_new ('glsys_assembleGlobal', 'Da', &
                                    rdestMatrix%RmatrixBlock(1,1)%NA, &
                                    cdataTypeLocal, &
                                    rdestMatrix%RmatrixBlock(1,1)%h_Da,&
                                    ST_NEWBLOCK_NOINIT)

      ! Set up KCOL and the matrix entries
      call glmatasm_Da99dble (rlocalMatrix,rdestMatrix%RmatrixBlock(1,1),&
                              Icolumns, Irows)

      ! Allocate a Kdiagonal in the destination matrix if necessary
      if (rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal .ne. ST_NOHANDLE) then
        call storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal,isize)
        if (isize .ne. rlocalMatrix%NEQ) then
          ! Release and reallocate
          call storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Da)
        end if
      end if

      if (rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal .eq. ST_NOHANDLE) then
        call storage_new ('glsys_assembleGlobal', 'Kdiagonal', rlocalMatrix%NEQ, &
                        ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal,&
                        ST_NEWBLOCK_NOINIT)
      end if

      ! Rebuild Kdiagonal
      call lsyssc_getbase_Kdiagonal(rdestMatrix%RmatrixBlock(1,1),p_Kdiagonal)
      call lsyssc_getbase_Kcol(rdestMatrix%RmatrixBlock(1,1),p_Kcol)
      call lsyssc_getbase_Kld(rdestMatrix%RmatrixBlock(1,1),p_Kld)
      call lsyssc_rebuildKdiagonal (p_Kcol, p_Kld, p_Kdiagonal, rdestMatrix%NEQ)

    end select

  end if

  ! Release the local matrix.
  ! Note that only those submatrices are released from the heap that we
  ! created by transposing them above!
  call lsysbl_releaseMatrix(rlocalMatrix)

  deallocate(Irows)
  deallocate(Icolumns)

  end subroutine glsys_assembleGlobal

  ! Auxiliary subroutines

  !------------------------------------------------------------------
  ! Set up general data of the destination matrix

  subroutine glmatasm_initDestination (rsourceMatrix,rdestMatrix, &
                                        cmatrixFormat,cdataType)

  ! The source block matrix
  type(t_matrixBlock), intent(in) :: rsourceMatrix

  ! The matrix to be initialised
  type(t_matrixBlock), intent(inout) :: rdestMatrix

  ! Target format of the matrix rdestMatrix. Standard is Format 9.
  integer, intent(in) :: cmatrixFormat

  ! Data type for the entries of rdestMatrix.
  ! Standard is double precision.
  integer, intent(in) :: cdataType

    if ((rdestMatrix%nblocksPerCol .lt. 1) .or. &
        (rdestMatrix%nblocksPerRow .lt. 1)) then
      ! Create a new 1x1 matrix if necessary.
      call lsysbl_releaseMatrix (rdestMatrix)
      call lsysbl_createEmptyMatrix (rdestMatrix,1)
    end if
    rdestMatrix%NEQ = rsourceMatrix%NEQ
    rdestMatrix%NCOLS = rsourceMatrix%NCOLS
    rdestMatrix%imatrixSpec = rsourceMatrix%imatrixSpec

    ! There is no appropriate discretisation or boundary condition
    ! structure for a global matrix! These only apply to
    ! scalar discretisations.
    nullify(rdestMatrix%p_rblockDiscrTest)
    nullify(rdestMatrix%p_rblockDiscrTrial)
    rdestMatrix%bidenticalTrialAndTest = .true.
    nullify(rdestMatrix%p_rdiscreteBC)
    nullify(rdestMatrix%p_rdiscreteBCfict)

    ! Initialise structural information of the submatrix
    rdestMatrix%RmatrixBlock(1,1)%NEQ = rsourceMatrix%NEQ
    rdestMatrix%RmatrixBlock(1,1)%NCOLS = rsourceMatrix%NCOLS
    rdestMatrix%RmatrixBlock(1,1)%cmatrixFormat = cmatrixFormat
    rdestMatrix%RmatrixBlock(1,1)%cdataType = cdataType
    rdestMatrix%RmatrixBlock(1,1)%imatrixSpec = 0
    rdestMatrix%RmatrixBlock(1,1)%isortStrategy = 0
    rdestMatrix%RmatrixBlock(1,1)%h_IsortPermutation = ST_NOHANDLE
    rdestMatrix%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP

  end subroutine glmatasm_initDestination

  !------------------------------------------------------------------
  ! Calculates the row-block offsets in the destination matrix
  subroutine glmatasm_getOffsets (rsourceMatrix,Icolumns,Irows)

  ! The source block matrix
  type(t_matrixBlock), intent(in), target :: rsourceMatrix

  ! Real column numbers in each block-column of the block matrix.
  ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
  integer, dimension(:), intent(out) :: Icolumns

  ! Real row numbers in each block-column of the block matrix.
  ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
  integer, dimension(:), intent(out) :: Irows

    ! local variables
    integer :: i,j

    Icolumns(:) = 0
    Irows(:) = 0

    ! Loop through all matrix blocks.
    do i=1,rsourceMatrix%nblocksPerCol

      do j=1,rsourceMatrix%nblocksPerRow
        ! When checking for the presence of the matrix, do not respect
        ! the scaling factor; we only want to get the size of the matrix
        ! columns/rows!
        if (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j,.true.)) then

          Icolumns(j+1) = rsourceMatrix%RmatrixBlock(i,j)%NCOLS
          Irows(i+1) = rsourceMatrix%RmatrixBlock(i,j)%NEQ

        end if
      end do

    end do

    ! Icolumns gives the number of real columns/rows. Add them together
    ! to get the real column/row numbers, each block column/row in the
    ! global matrix starts with.
    Icolumns(1) = 1
    Irows(1) = 1
    do i=2,rsourceMatrix%nblocksPerRow
      Icolumns(i) = Icolumns(i) + Icolumns(i-1)
    end do
    do i=2,rsourceMatrix%nblocksPerCol
      Irows(i) = Irows(i) + Irows(i-1)
    end do

  end subroutine glmatasm_getOffsets

  !------------------------------------------------------------------
  ! Calculates NA and the KLD row structure of the global matrix
  subroutine glmatasm_KLD (rsourceMatrix,rdestMatrix,Irows)

  ! The source block matrix
  type(t_matrixBlock), intent(in), target :: rsourceMatrix

  ! The scalar submatrix which KLD is to be initialised.
  ! KLD must exist and be filled with 0.
  type(t_matrixScalar), intent(inout) :: rdestMatrix

  ! Real row numbers in each block-column of the block matrix.
  ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
  integer, dimension(:), intent(in) :: Irows

    ! local variables
    integer :: i,j
    integer :: irow
    integer :: irowoffset,narow
    integer, dimension(:), pointer :: p_Kld,p_KldDest
    type(t_matrixScalar), pointer :: p_rmatrix

    ! Create a new, empty KLD.
    call lsyssc_getbase_Kld (rdestMatrix,p_KldDest)

    ! Loop through all matrix blocks.
    do i=1,rsourceMatrix%nblocksPerCol

      ! Get the starting position of this block-row in the global matrix
      irowoffset = Irows(i)-1

      narow = 0

      do j=1,rsourceMatrix%nblocksPerRow
        if (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) then

          p_rmatrix => rsourceMatrix%RmatrixBlock(i,j)
          call lsyssc_getbase_Kld (p_rmatrix,p_Kld)

          ! Loop through all lines in the matrix and add the number of
          ! entries per row to get KLD:
          do irow = 1,p_rmatrix%NEQ
            ! The first entriy in KLD is always one; add the length
            ! of the line to KLD one element shifted, so calculation
            ! ok KLD is easier later.
            p_KldDest(1+irow+irowoffset) = &
              p_KldDest(1+irow+irowoffset) + (p_Kld(irow+1)-p_Kld(irow))
          end do

          ! Calculate the number of entries in this matrix-row-block
          narow = narow + rsourceMatrix%RmatrixBlock(i,j)%NA

        end if
      end do

    end do

    ! Now we have:
    ! KLD(1) = 0,
    ! KLD(2) = Number of entries in row 1,
    ! KLD(3) = Number of entries in row 2, etc.
    ! Sum up the values to get the actual KLD.
    p_KldDest(1) = 1
    do irow = 1,rsourceMatrix%NEQ
      p_KldDest(irow+1) = p_KldDest(irow+1) + p_KldDest(irow)
    end do

    ! and so we have NA.
    rdestMatrix%NA = p_KldDest(rsourceMatrix%NEQ+1)-1

  end subroutine glmatasm_KLD

  !------------------------------------------------------------------
  ! Calculates the KCOL column structure of the global matrix
  ! The KLD/NA structure must already be present in rdestMatrix!
  !
  ! Matrix-structure-9 source -> Matrix-structure-9 destination,
  subroutine glmatasm_Kcol99dble (rsourceMatrix,rdestMatrix,&
                                  Icolumns, Irows)

  ! The source block matrix
  type(t_matrixBlock), intent(in), target :: rsourceMatrix

  ! The scalar submatrix which KLD is to be initialised
  type(t_matrixScalar), intent(inout) :: rdestMatrix

  ! Real column numbers in each block-column of the block matrix.
  ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
  integer, dimension(:), intent(in) :: Icolumns

  ! Real row numbers in each block-column of the block matrix.
  ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
  integer, dimension(:), intent(in) :: Irows

    ! local variables
    integer :: i,j,h_KldTmp
    integer :: irow,ncols,irowGlobal
    integer :: ioffsetGlobal,ioffsetLocal
    integer, dimension(:), pointer :: p_Kld,p_KldDest,p_KldTmp
    integer, dimension(:), pointer :: p_Kcol,p_KcolDest
    type(t_matrixScalar), pointer :: p_rmatrix

    ! Create a copy of KLD which we use for storing the index
    ! how many entries are already allocated in each row.
    h_KldTmp = ST_NOHANDLE
    call storage_copy (rdestMatrix%h_Kld,h_KldTmp)
    call storage_getbase_int (h_KldTmp,p_KldTmp)

    ! Get KCol,Kld
    call lsyssc_getbase_Kcol (rdestMatrix,p_KcolDest)
    call lsyssc_getbase_Kld (rdestMatrix,p_KldDest)

    ! Loop through all matrix subblocks
    do i=1,rsourceMatrix%nblocksPerCol

      ! Global row index?
      irowGlobal = Irows(i)-1

      do j=1,rsourceMatrix%nblocksPerRow

        if (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) then

          ! Get the submatrix
          p_rmatrix => rsourceMatrix%RmatrixBlock(i,j)

          ! Get the local matrix pointers / column structure
          call lsyssc_getbase_Kcol (p_rmatrix,p_Kcol)
          call lsyssc_getbase_Kld (p_rmatrix,p_Kld)

          ! Loop through all rows to append them to the current rows.
          do irow = 1,p_rmatrix%NEQ

            ! How many elements to append?
            ncols = p_Kld(irow+1)-p_Kld(irow)

            ! Position of the data?
            ioffsetGlobal = p_KldTmp(irowGlobal+irow)
            ioffsetLocal  = p_Kld(irow)

            ! Copy matrix data and the column numbers
            p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
              p_Kcol(ioffsetLocal:ioffsetLocal+ncols-1)

            ! Increase the column numbers by the global column number
            ! of that matrix-block-column
            p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
              p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) + Icolumns(j)-1

            ! Increase the counter/index position array for how
            ! many elements are added to that row.
            p_KldTmp(irowGlobal+irow) = p_KldTmp(irowGlobal+irow) + ncols

          end do ! irow

        end if ! neq != 0

      end do ! j
    end do ! i

    ! Release the temp array
    call storage_free (h_KldTmp)

  end subroutine glmatasm_Kcol99dble

  !------------------------------------------------------------------
  ! Transfers the entries of the local matrices into the
  ! global matrix.
  ! The KLD/NA structure must already be present in rdestMatrix!
  !
  ! Matrix-structure-9 source -> Matrix-structure-9 destination,
  ! double precision vection
  subroutine glmatasm_Da99dble (rsourceMatrix,rdestMatrix,&
                                Icolumns, Irows)

  ! The source block matrix
  type(t_matrixBlock), intent(in), target :: rsourceMatrix

  ! The scalar submatrix which KLD is to be initialised
  type(t_matrixScalar), intent(inout) :: rdestMatrix

  ! Real column numbers in each block-column of the block matrix.
  ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
  integer, dimension(:), intent(in) :: Icolumns

  ! Real row numbers in each block-column of the block matrix.
  ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
  integer, dimension(:), intent(in) :: Irows

    ! local variables
    integer :: i,j,h_KldTmp
    real(DP) :: dscale
    integer :: irow,ncols,irowGlobal
    integer :: ioffsetGlobal,ioffsetLocal
    integer, dimension(:), pointer :: p_Kld,p_KldDest,p_KldTmp
    real(DP), dimension(:), pointer :: p_Da,p_DaDest
    type(t_matrixScalar), pointer :: p_rmatrix

    ! Create a copy of KLD which we use for storing the index
    ! how many entries are already allocated in each row.
    h_KldTmp = ST_NOHANDLE
    call storage_copy (rdestMatrix%h_Kld,h_KldTmp)
    call storage_getbase_int (h_KldTmp,p_KldTmp)

    ! Get Kld
    call lsyssc_getbase_Kld (rdestMatrix,p_KldDest)

    ! Get destination data arrays
    call lsyssc_getbase_double (rdestMatrix,p_DaDest)

    ! Loop through all matrix subblocks
    do i=1,rsourceMatrix%nblocksPerCol

      ! Global row index?
      irowGlobal = Irows(i)-1

      do j=1,rsourceMatrix%nblocksPerRow

        if (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) then

          ! Get the submatrix
          p_rmatrix => rsourceMatrix%RmatrixBlock(i,j)

          ! Get the local matrix pointers structure
          call lsyssc_getbase_Kld (p_rmatrix,p_Kld)
          call lsyssc_getbase_double (p_rmatrix,p_Da)
          dscale = p_rmatrix%dscaleFactor

          ! Loop through all rows to append them to the current rows.
          do irow = 1,p_rmatrix%NEQ

            ! How many elements to append?
            ncols = p_Kld(irow+1)-p_Kld(irow)

            ! Position of the data?
            ioffsetGlobal = p_KldTmp(irowGlobal+irow)
            ioffsetLocal  = p_Kld(irow)

            ! Copy matrix data
            p_DaDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
              dscale * p_Da(ioffsetLocal:ioffsetLocal+ncols-1)

            ! Increase the counter/index position array for how
            ! many elements are added to that row.
            p_KldTmp(irowGlobal+irow) = p_KldTmp(irowGlobal+irow) + ncols

          end do ! irow

        end if ! neq != 0

      end do ! j
    end do ! i

    ! Release the temp array
    call storage_free (h_KldTmp)

  end subroutine glmatasm_Da99dble

  !------------------------------------------------------------------
  ! Calculates the KCOL column structure of the global matrix
  ! and transfers the entries of the local matrices into the
  ! global matrix.
  ! The KLD/NA structure must already be present in rdestMatrix!
  !
  ! Matrix-structure-9 source -> Matrix-structure-9 destination,
  ! double precision vection
  subroutine glmatasm_KcolDa99dble (rsourceMatrix,rdestMatrix,&
                                    Icolumns, Irows)

  ! The source block matrix
  type(t_matrixBlock), intent(in), target :: rsourceMatrix

  ! The scalar submatrix which KLD is to be initialised
  type(t_matrixScalar), intent(inout) :: rdestMatrix

  ! Real column numbers in each block-column of the block matrix.
  ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
  integer, dimension(:), intent(in) :: Icolumns

  ! Real row numbers in each block-column of the block matrix.
  ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
  integer, dimension(:), intent(in) :: Irows

    ! local variables
    integer :: i,j,h_KldTmp
    real(DP) :: dscale
    integer :: irow,ncols,irowGlobal
    integer :: ioffsetGlobal,ioffsetLocal
    integer, dimension(:), pointer :: p_Kld,p_KldDest,p_KldTmp
    integer, dimension(:), pointer :: p_Kcol,p_KcolDest
    real(DP), dimension(:), pointer :: p_Da,p_DaDest
    type(t_matrixScalar), pointer :: p_rmatrix

    ! Create a copy of KLD which we use for storing the index
    ! how many entries are already allocated in each row.
    h_KldTmp = ST_NOHANDLE
    call storage_copy (rdestMatrix%h_Kld,h_KldTmp)
    call storage_getbase_int (h_KldTmp,p_KldTmp)

    ! Get KCol,Kld
    call lsyssc_getbase_Kcol (rdestMatrix,p_KcolDest)
    call storage_getbase_int (rdestMatrix%h_Kld,p_KldDest)

    ! Get destination data arrays
    call lsyssc_getbase_double (rdestMatrix,p_DaDest)

    ! Loop through all matrix subblocks
    do i=1,rsourceMatrix%nblocksPerCol

      ! Global row index?
      irowGlobal = Irows(i)-1

      do j=1,rsourceMatrix%nblocksPerRow

        if (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) then

          ! Get the submatrix
          p_rmatrix => rsourceMatrix%RmatrixBlock(i,j)

          ! Get the local matrix pointers / column structure
          call lsyssc_getbase_Kcol (p_rmatrix,p_Kcol)
          call lsyssc_getbase_Kld (p_rmatrix,p_Kld)
          call lsyssc_getbase_double (p_rmatrix,p_Da)
          dscale = p_rmatrix%dscaleFactor

          ! Loop through all rows to append them to the current rows.
          do irow = 1,p_rmatrix%NEQ

            ! How many elements to append?
            ncols = p_Kld(irow+1)-p_Kld(irow)

            ! Position of the data?
            ioffsetGlobal = p_KldTmp(irowGlobal+irow)
            ioffsetLocal  = p_Kld(irow)

            ! Copy matrix data and the column numbers
            p_DaDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
              dscale * p_Da(ioffsetLocal:ioffsetLocal+ncols-1)
            p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
              p_Kcol(ioffsetLocal:ioffsetLocal+ncols-1)

            ! Increase the column numbers by the global column number
            ! of that matrix-block-column
            p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
              p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) + Icolumns(j)-1

            ! Increase the counter/index position array for how
            ! many elements are added to that row.
            p_KldTmp(irowGlobal+irow) = p_KldTmp(irowGlobal+irow) + ncols

          end do ! irow

        end if ! neq != 0

      end do ! j
    end do ! i

    ! Release the temp array
    call storage_free (h_KldTmp)

  end subroutine glmatasm_KcolDa99dble

end module globalsystem
