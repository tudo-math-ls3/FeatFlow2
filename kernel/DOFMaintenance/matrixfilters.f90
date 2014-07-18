!##############################################################################
!# ****************************************************************************
!# <name> matrixfilters </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of filters that can be applied to a scalar
!# or block matrix. They can be used e.g. to modify a matrix during a nonlinear
!# solution process to prepare it for a linear solver (like implementing
!# boundary conditions).
!#
!# Typical matrix filters are those working together with discrete
!# boundary conditions. For Dirichlet boundary conditions e.g. some lines
!# of the global matrix are usually replaced by unit vectors. The corresponding
!# matrix modification routine is a kind of "filter" for a matrix and
!# can be found in this module.
!#
!# The following routines can be found here:
!#
!#  1.) matfil_discreteBC
!#      -> Apply the "discrete boundary conditions for matrices" filter
!#         onto a given (block) matrix.
!#         This e.g. replaces lines of a matrix corresponding to Dirichlet DOF`s
!#         by unit vectors for all scalar submatrices where configured.
!#
!#  2.) matfil_discreteFBC
!#      -> Apply the "discrete fictitious boundary conditions for matrices"
!#         filter onto a given (block) matrix.
!#         This e.g. replaces lines of a matrix corresponding to Dirichlet DOF`s
!#         by unit vectors for all scalar submatrices where configured.
!#
!#  3.) matfil_discreteNLSlipBC
!#      -> Implement slip boundary conditions into a matrix.
!#         Slip boundary conditions are not implemented by matfil_discreteBC!
!#         They must be implemented separately, as they are a special type
!#         of boundary conditions.
!#
!#  4.) matfil_normaliseToL20
!#      -> Adds a row to a matrix containing "ones". This implements the
!#         explicit condition "int_{\Omega} v = 0" into a matrix.
!#
!#  5.) matfil_normaliseToL20LmassMat
!#      -> Uses a lumoed mass matrix to impose the condition
!#         "int_{\Omega} v = 0" into a matrix.
!#
!# Auxiliary routines:
!#
!#  1.) matfil_imposeDirichletBC
!#      -> Imposes Dirichlet BC`s into a scalar matrix
!#
!#  2.) matfil_imposeNLSlipBC
!#      -> Imposes nonlinear slip boundary conditions into a scalar matrix
!#
!#  3.) matfil_imposeDirichletFBC
!#      -> Imposes Difichlet BC`s for fictitious boundary components
!#         into a scalar matrix
!#
!# </purpose>
!##############################################################################

module matrixfilters

!$ use omp_lib
  use fsystem
  use storage
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use discretebc
  use discretefbc
  use dofmapping
  use matrixmodification
  use genoutput
  use sortstrategybase

  implicit none

  private

  public :: matfil_discreteBC
  public :: matfil_discreteFBC
  public :: matfil_discreteNLSlipBC
  public :: matfil_normaliseToL20
  public :: matfil_imposeDirichletBC
  public :: matfil_imposeNLSlipBC
  public :: matfil_imposeDirichletFBC

contains

! *****************************************************************************
! Scalar matrix filters
! *****************************************************************************

  ! ***************************************************************************

!<subroutine>

  subroutine matfil_imposeDirichletBC (rmatrix,boffDiag,rdbcStructure)

!<description>
  ! Implements discrete Dirichlet BC`s into a scalar matrix.
  ! This is normally done by replacing some lines of the matrix rmatrix
  ! (those belonging to Dirichlet nodes) by unit vectors or by zero-vectors
  ! (depending on whether the matrix is a "diagonal" matrix or an
  ! "off-diagonal" matrix in a larger block-system).
!</description>

!<input>

  ! The t_discreteBCDirichlet that describes the discrete Dirichlet BC`s
  type(t_discreteBCDirichlet), intent(in), target  :: rdbcStructure

  ! Off-diagonal matrix.
  ! If this is present and set to TRUE, it is assumed that the matrix is not
  ! a main, guiding system matrix, but an "off-diagonal" matrix in a
  ! system with block-matrices (e.g. a matrix at position (2,1), (3,1),...
  ! or somewhere else in a block system). This modifies the way,
  ! boundary conditions are implemented into the matrix.
  logical :: boffDiag

!</input>

!<inputoutput>

  ! The scalar matrix where the boundary conditions should be imposed.
  type(t_matrixScalar), intent(inout), target :: rmatrix

!</inputoutput>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_idx
  integer, parameter :: NBLOCKSIZE = 1000
  integer, dimension(1000) :: Idofs
  integer, dimension(:), pointer :: p_Iperm
  integer i,ilenleft

  ! If nDOF=0, there are no DOF`s the current boundary condition segment,
  ! so we do not have to do anything. Maybe the case if the user selected
  ! a boundary region that is just too small.
  if (rdbcStructure%nDOF .eq. 0) return

  ! Get pointers to the structures. For the vector, get the pointer from
  ! the storage management.

  call storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.

  if (.not.associated(p_idx)) then
    call output_line ("DBC not configured",&
        OU_CLASS_ERROR,OU_MODE_STD,"matfil_imposeDirichletBC")
    call sys_halt()
  end if

  ! Ensure that the BC and the matrix are compatible.
  if (rdbcStructure%NEQ .ne. rmatrix%NEQ) then
    call output_line("BC and vector incompatible.",&
        OU_CLASS_ERROR,OU_MODE_STD,"matfil_imposeDirichletBC")
    call sys_halt()
  end if

  ! Only handle nDOF DOF`s, not the complete array!
  ! Probably, the array is longer (e.g. has the length of the vector), but
  ! contains only some entries...
  !
  ! Is the matrix sorted?
  if ((.not. rmatrix%browsSorted) .or. (.not. associated(rmatrix%p_rsortStrategyRows))) then
    ! Use mmod_replaceLinesByUnit/mmod_replaceLinesByZero to replace the
    ! corresponding rows in the matrix by unit vectors.
    ! For more complicated FE spaces, this might have
    ! to be modified in the future...
    !
    ! We use mmod_replaceLinesByUnit for "diagonal" matrix blocks (main
    ! system matrices) and mmod_replaceLinesByZero for "off-diagonal"
    ! matrix blocks (in case of larger block systems)

    if (boffDiag) then
      call mmod_replaceLinesByZero (rmatrix,p_idx(1:rdbcStructure%nDOF))
    else
      call mmod_replaceLinesByUnit (rmatrix,p_idx(1:rdbcStructure%nDOF))
    end if
  else
    ! Ok, matrix is sorted, so we have to filter all the DOF`s through the
    ! permutation before using them for implementing boundary conditions.
    ! We do this in blocks with 1000 DOF`s each to prevent the stack
    ! from being destroyed!
    !
    ! Get the permutation from the matrix - or more precisely, the
    ! back-permutation, as we need this one for the loop below.
    call sstrat_getSortedPosInfo(rmatrix%p_rsortStrategyRows,p_Iperm)

    ! How many entries to handle in the first block?
    ilenleft = min(rdbcStructure%nDOF,NBLOCKSIZE)

    ! Loop through the DOF-blocks
    do i=0,rdbcStructure%nDOF / NBLOCKSIZE
      ! Filter the DOF`s through the permutation
      Idofs(1:ilenleft) = p_Iperm(p_idx(1+i*NBLOCKSIZE:i*NBLOCKSIZE+ilenleft))

      ! And implement the BC`s with mmod_replaceLinesByUnit/
      ! mmod_replaceLinesByZero, depending on whether the matrix is a
      ! "main" system matrix or an "off-diagonal" system matrix in a larger
      ! block system.
      if (boffDiag) then
        call mmod_replaceLinesByZero (rmatrix,Idofs(1:ilenleft))
      else
        call mmod_replaceLinesByUnit (rmatrix,Idofs(1:ilenleft))
      end if

      ! How many DOF`s are left?
      ilenleft = min(rdbcStructure%nDOF-(i+1)*NBLOCKSIZE,NBLOCKSIZE)
    end do

  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine matfil_imposeNLSlipBC (rmatrix,boffDiag,bforprec,rslipBCStructure)

!<description>
  ! Implements discrete Slip BC`s into a scalar matrix.
  ! Slip BC`s are treated like Dirichlet BC`s.
  ! This is normally done by replacing some lines of the matrix rmatrix
  ! (those belonging to Dirichlet nodes) by unit vectors or by zero-vectors
  ! (depending on whether the matrix is a "diagonal" matrix or an
  ! "off-diagonal" matrix in a larger block-system).
!</description>

!<input>
  ! Prepare matrix for preconditioning.
  ! =false: Slip-nodes are handled as Dirichlet. Rows in the
  !         matrix are replaced by unit vectors.
  ! =true : Standard. Prepare matrix for preconditioning.
  !         Only the off-diagonals of rows corresponding to
  !         the slip-nodes are cleared.
  !         The matrix therefore is prepared to act as
  !         Jacobi-preconditioner for all slip-nodes.
  logical, intent(in) :: bforprec

  ! The t_discreteBCSlip that describes the discrete Slip BC`s
  type(t_discreteBCSlip), intent(in), target  :: rslipBCStructure

  ! Off-diagonal matrix.
  ! If this is present and set to TRUE, it is assumed that the matrix is not
  ! a main, guiding system matrix, but an "off-diagonal" matrix in a
  ! system with block-matrices (e.g. a matrix at position (2,1), (3,1),...
  ! or somewhere else in a block system). This modifies the way,
  ! boundary conditions are implemented into the matrix.
  logical, intent(in) :: boffDiag
!</input>

!<inputoutput>

  ! The scalar matrix where the boundary conditions should be imposed.
  type(t_matrixScalar), intent(inout), target :: rmatrix

!</inputoutput>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_idx
  integer, parameter :: NBLOCKSIZE = 1000
  integer, dimension(1000) :: Idofs
  integer, dimension(:), pointer :: p_Iperm
  integer i,ilenleft

  ! If nDOF=0, there are no DOF`s the current boundary condition segment,
  ! so we do not have to do anything. Maybe the case if the user selected
  ! a boundary region that is just too small.
  if (rslipBCStructure%nDOF .eq. 0) return

  ! Get pointers to the structures. For the vector, get the pointer from
  ! the storage management.

  call storage_getbase_int(rslipBCStructure%h_IslipDOFs,p_idx)

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.

  if (.not.associated(p_idx)) then
    call output_line ("DBC not configured",&
        OU_CLASS_ERROR,OU_MODE_STD,"matfil_imposeNLSlipBC")
    call sys_halt()
  end if

  ! Ensure that the BC and the matrix are compatible.
  if (rslipBCStructure%NEQ .ne. rmatrix%NEQ) then
    call output_line("BC and vector incompatible.",&
        OU_CLASS_ERROR,OU_MODE_STD,"matfil_imposeNLSlipBC")
    call sys_halt()
  end if

  ! Only handle nDOF DOF`s, not the complete array!
  ! Probably, the array is longer (e.g. has the length of the vector), but
  ! contains only some entries...
  !
  ! Is the matrix sorted?
  if ((.not. rmatrix%browsSorted) .or. (.not. associated(rmatrix%p_rsortStrategyRows))) then
    ! Use mmod_clearOffdiags/mmod_replaceLinesByZero clear the
    ! offdiagonals of the system matrix.
    ! For more complicated FE spaces, this might have
    ! to be modified in the future...
    !
    ! We use mmod_clearOffdiags for "diagonal" matrix blocks (main
    ! system matrices) and mmod_replaceLinesByZero for "off-diagonal"
    ! matrix blocks (in case of larger block systems)

    if (boffDiag) then
      call mmod_replaceLinesByZero (rmatrix,p_idx(1:rslipBCStructure%nDOF))
    else
      if (bforprec) then
        call mmod_clearOffdiags (rmatrix,p_idx(1:rslipBCStructure%nDOF))
      else
        call mmod_replaceLinesByUnit (rmatrix,p_idx(1:rslipBCStructure%nDOF))
      end if
    end if
  else
    ! Ok, matrix is sorted, so we have to filter all the DOF`s through the
    ! permutation before using them for implementing boundary conditions.
    ! We do this in blocks with 1000 DOF`s each to prevent the stack
    ! from being destroyed!
    !
    ! Get the permutation from the matrix - or more precisely, the
    ! back-permutation, as we need this one for the loop below.
    call sstrat_getSortedPosInfo(rmatrix%p_rsortStrategyRows,p_Iperm)

    ! How many entries to handle in the first block?
    ilenleft = min(rslipBCStructure%nDOF,NBLOCKSIZE)

    ! Loop through the DOF-blocks
    do i=0,rslipBCStructure%nDOF / NBLOCKSIZE
      ! Filter the DOF`s through the permutation
      Idofs(1:ilenleft) = p_Iperm(p_idx(1+i*NBLOCKSIZE:i*NBLOCKSIZE+ilenleft))

      ! And implement the BC`s with mmod_clearOffdiags/
      ! mmod_replaceLinesByZero, depending on whether the matrix is a
      ! "main" system matrix or an "off-diagonal" system matrix in a larger
      ! block system.
      if (boffDiag) then
        call mmod_replaceLinesByZero (rmatrix,Idofs(1:ilenleft))
      else
        if (bforprec) then
          call mmod_clearOffdiags (rmatrix,Idofs(1:ilenleft))
        else
          call mmod_replaceLinesByUnit (rmatrix,Idofs(1:ilenleft))
        end if
      end if

      ! How many DOF`s are left?
      ilenleft = min(rslipBCStructure%nDOF-(i+1)*NBLOCKSIZE,NBLOCKSIZE)
    end do

  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine matfil_imposeDirichletFBC (rmatrix,boffDiag,rdbcStructure)

!<description>
  ! Implements discrete Dirichlet BC`s of fictitious boundary components
  ! into a scalar matrix.
  ! This is normally done by replacing some lines of the matrix rmatrix
  ! (those belonging to Dirichlet nodes) by unit vectors or by zero-vectors
  ! (depending on whether the matrix is a "diagonal" matrix or an
  ! "off-diagonal" matrix in a larger block-system).
!</description>

!<input>

  ! The t_discreteFBCDirichlet that describes the discrete Dirichlet BC`s
  type(t_discreteFBCDirichlet), intent(in), target  :: rdbcStructure

  ! Off-diagonal matrix.
  ! If this is present and set to TRUE, it is assumed that the matrix is not
  ! a main, guiding system matrix, but an "off-diagonal" matrix in a
  ! system with block-matrices (e.g. a matrix at position (2,1), (3,1),...
  ! or somewhere else in a block system). This modifies the way,
  ! boundary conditions are implemented into the matrix.
  logical :: boffDiag

!</input>

!<inputoutput>

  ! The scalar matrix where the boundary conditions should be imposed.
  type(t_matrixScalar), intent(inout), target :: rmatrix

!</inputoutput>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_idx
  integer, parameter :: NBLOCKSIZE = 1000
  integer, dimension(1000) :: Idofs
  integer, dimension(:), pointer :: p_Iperm
  integer i,ilenleft

  ! If nDOF=0, there are no DOF`s the current boundary condition segment,
  ! so we do not have to do anything. Maybe the case if the user selected
  ! a boundary region that is just too small.
  if (rdbcStructure%nDOF .eq. 0) return

  ! Get pointers to the structures. For the vector, get the pointer from
  ! the storage management.

  call storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.

  if (.not.associated(p_idx)) then
    call output_line ("DBC not configured",&
        OU_CLASS_ERROR,OU_MODE_STD,"matfil_imposeDirichletFBC")
    call sys_halt()
  end if

  ! Ensure that the BC and the matrix are compatible.
  if (rdbcStructure%NEQ .ne. rmatrix%NEQ) then
    call output_line("BC and vector incompatible.",&
        OU_CLASS_ERROR,OU_MODE_STD,"matfil_imposeDirichletFBC")
    call sys_halt()
  end if

  ! Only handle nDOF DOF`s, not the complete array!
  ! Probably, the array is longer (e.g. has the length of the vector), but
  ! contains only some entries...
  !
  ! Is the matrix sorted?
  if ((.not. rmatrix%browsSorted) .or. (.not. associated(rmatrix%p_rsortStrategyRows))) then
    ! Use mmod_replaceLinesByUnit/mmod_replaceLinesByZero to replace the
    ! corresponding rows in the matrix by unit vectors.
    ! For more complicated FE spaces, this might have
    ! to be modified in the future...
    !
    ! We use mmod_replaceLinesByUnit for "diagonal" matrix blocks (main
    ! system matrices) and mmod_replaceLinesByZero for "off-diagonal"
    ! matrix blocks (in case of larger block systems)

    if (boffDiag) then
      call mmod_replaceLinesByZero (rmatrix,p_idx(1:rdbcStructure%nDOF))
    else
      call mmod_replaceLinesByUnit (rmatrix,p_idx(1:rdbcStructure%nDOF))
    end if
  else
    ! Ok, matrix is sorted, so we have to filter all the DOF`s through the
    ! permutation before using them for implementing boundary conditions.
    ! We do this in blocks with 1000 DOF`s each to prevent the stack
    ! from being destroyed!
    !
    ! Get the permutation from the matrix - or more precisely, the
    ! back-permutation, as we need this one for the loop below.
    call sstrat_getSortedPosInfo(rmatrix%p_rsortStrategyRows,p_Iperm)

    ! How many entries to handle in the first block?
    ilenleft = min(rdbcStructure%nDOF,NBLOCKSIZE)

    ! Loop through the DOF-blocks
    do i=0,rdbcStructure%nDOF / NBLOCKSIZE
      ! Filter the DOF`s through the permutation
      Idofs(1:ilenleft) = p_Iperm(p_idx(1+i*NBLOCKSIZE:i*NBLOCKSIZE+ilenleft))

      ! And implement the BC`s with mmod_replaceLinesByUnit/
      ! mmod_replaceLinesByZero, depending on whether the matrix is a
      ! "main" system matrix or an "off-diagonal" system matrix in a larger
      ! block system.
      if (boffDiag) then
        call mmod_replaceLinesByZero (rmatrix,Idofs(1:ilenleft))
      else
        call mmod_replaceLinesByUnit (rmatrix,Idofs(1:ilenleft))
      end if
    end do

  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine matfil_normaliseToL20 (rmatrix,istartColumn,iendColumn)

!<description>
  ! Modifies a scalar matrix to add an additional equation which implements
  ! <tex>$\int_{\Omega} v = 0$</tex> for all <tex>$v$</tex>. This adds another row to the
  ! matrix which sums up all vector entries when being applied to a vector.  Note that
  ! this produces a rectangular matrix with NROW+1 rows in comparison to the original one!
  !
  ! This filter can be used e.g. to impose the condition <tex>$\int_{\Omega} v = 0$</tex>
  ! to a matrix which is passed to an external linear solver (like UMFPACK)
  ! which cannot use filtering to cope with indefiniteness!
  !
  ! WARNING: UNTESTED!!!
!</description>

!<input>
  ! OPTIONAL: Start column in the matrix.
  ! This parameter can specify where to start the vector sum.
  ! If not specified, 1 is assumed.
  integer, intent(in), optional :: istartColumn

  ! OPTIONAL: End column in the matrix.
  ! This parameter can specify where to end the vector sum.
  ! If not specified, rmatrix%NCOLS is assumed.
  integer, intent(in), optional :: iendColumn
!</input>

!<inputoutput>
  ! The matrix which is to be modified.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! Local variables
    type(t_matrixScalar) :: rmatrixTemp
    integer :: irowLen,istart,iend,i
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    integer, dimension(:), pointer :: p_Kcol

    ! Matrix must not be transposed
    if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
      call output_line (&
          "Virtually transposed matrices not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"matfil_normaliseToL20")
      call sys_halt()
    end if

    ! Only double precision supported.
    if (rmatrix%cdataType .ne. ST_DOUBLE) then
      call output_line (&
          "Only double precision matrices supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"matfil_normaliseToL20")
      call sys_halt()
    end if

    ! If structure and/or content is shared, duplicate the matrix to
    ! make the entries belong to rmatrix.
    if (lsyssc_isMatrixStructureShared(rmatrix) .and. &
        lsyssc_isMatrixContentShared(rmatrix)) then
      rmatrixTemp = rmatrix
      call lsyssc_duplicateMatrix (rmatrixTemp,rmatrix,LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
    else if (lsyssc_isMatrixStructureShared(rmatrix)) then
      rmatrixTemp = rmatrix
      call lsyssc_duplicateMatrix (rmatrixTemp,rmatrix,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
    else if (lsyssc_isMatrixContentShared(rmatrix)) then
      rmatrixTemp = rmatrix
      call lsyssc_duplicateMatrix (rmatrixTemp,rmatrix,LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
    end if

    ! Length of the row
    istart = 1
    iend = rmatrix%NCOLS

    if (present(istartColumn)) &
      istart = max(istart,min(istartColumn,rmatrix%NCOLS))

    if (present(iendColumn)) &
      iend = max(istart,min(iendColumn,rmatrix%NCOLS))

    irowLen = iend-istart+1

    ! Add another row to the matrix
    call lsyssc_resizeMatrix (rmatrix, rmatrix%NEQ+1, &
        rmatrix%NCOLS, rmatrix%NA+irowLen, .false.)

    rmatrix%NEQ = rmatrix%NEQ + 1
    rmatrix%NA = rmatrix%NA+irowLen

    ! Add a row containing "1". TOgether with a RHS "0",
    ! this implements "sum_j a_ij v_j = 0".
    call lsyssc_getbase_double (rmatrix,p_Ddata)
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX1)
      istart = istart + rmatrix%NEQ * rmatrix%NCOLS
      iend = iend + rmatrix%NEQ * rmatrix%NCOLS
      do i=istart,iend
        p_Ddata(i) = 1.0_DP
      end do

    case (LSYSSC_MATRIX7)
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)

      ! New Kld entry
      p_Kld(rmatrix%NEQ+1) = p_Kld(rmatrix%NEQ) + irowLen

      ! Insert the new row
      do i=p_Kld(rmatrix%NEQ),p_Kld(rmatrix%NEQ+1)-1
        p_Ddata(i) = 1.0_DP
        p_Kcol(i) = istart + (p_Kld(rmatrix%NEQ)-i)
      end do

    case (LSYSSC_MATRIX9)
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      call lsyssc_getbase_Kld (rmatrix,p_Kdiagonal)

      ! New Kld entry
      p_Kld(rmatrix%NEQ+1) = p_Kld(rmatrix%NEQ) + irowLen

      ! Insert the new row
      do i=p_Kld(rmatrix%NEQ),p_Kld(rmatrix%NEQ+1)-1
        p_Ddata(i) = 1.0_DP
        p_Kcol(i) = istart + (p_Kld(rmatrix%NEQ)-i)
      end do

      ! New Kdiagonal entry
      p_Kdiagonal(rmatrix%NEQ) = p_Kld(rmatrix%NEQ) + rmatrix%NCOLS-istart

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine matfil_imposeFeastMirrorBC (rmatrix,boffDiag,rfmbcStructure)

!<description>
  ! Implements discrete Feast mirror BC`s into a scalar matrix.
  ! The FEAST mirror boundary condition is basically a Neumann boundary
  ! condition which is used for domain decomposition.
  ! One assumes that there is an additional "virtual" layer of cells added to
  ! a boundary edge. This leads to a slight matrix modification for all
  ! DOF`s on that boundary edge.
  ! Example: For a 5-point stencil with <tex>$Q_1$</tex>, boundary DOF`s get matrix
  ! weights "2, 1, -1/2, -1/2" (center, left, top, bottom), while inner
  ! points get matrix weights "4, -1, -1, -1, -1" (center and all surroundings).
  ! To make boundary DOF`s behave like inner DOF`s, the entries in
  ! the matrices belonging to such an edge have to be doubled,
  ! leading to "4, -1, -1".
  ! So this filter loops through the matrix and doubles all matrix entries
  ! that belong to DOF`s on FEAST mirror boundary edges.
!</description>

!<input>

  ! The t_discreteBCfeastMirror that describes the discrete FEAST mirror BC`s
  type(t_discreteBCfeastMirror), intent(in), target  :: rfmbcStructure

  ! Off-diagonal matrix.
  ! If this is present and set to TRUE, it is assumed that the matrix is not
  ! a main, guiding system matrix, but an "off-diagonal" matrix in a
  ! system with block-matrices (e.g. a matrix at position (2,1), (3,1),...
  ! or somewhere else in a block system). This modifies the way,
  ! boundary conditions are implemented into the matrix.
  logical :: boffDiag

!</input>

!<inputoutput>

  ! The scalar matrix where the boundary conditions should be imposed.
  type(t_matrixScalar), intent(inout), target :: rmatrix

!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_ImirrorDOFs,p_ImirrorDOFsClosed
    integer :: i,j
    real(DP), dimension(:), pointer :: p_Da
    integer, dimension(:), pointer :: p_Kcol,p_Iperm,p_IpermInverse
    integer, dimension(:), pointer :: p_Kld
    integer :: ia
    integer :: idof
    real(DP) :: dmirrorWeight

    ! Offdiagonal matrices are not processed by this routine up to now.
    if (boffDiag) return

    ! Impose the DOF value directly into the vector - more precisely, into the
    ! components of the subvector that is indexed by icomponent.

    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Only matrix format 7 and 9 supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"matfil_imposeFeastMirrorBC")
      call sys_halt()
    end if

    if (rmatrix%cdataType .ne. ST_DOUBLE) then
      call output_line ("Matrix must be double precision!",&
          OU_CLASS_ERROR,OU_MODE_STD,"matfil_imposeFeastMirrorBC")
      call sys_halt()
    end if

    if (rfmbcStructure%icomponent .eq. 0) then
      call output_line ("FMBC not configured!",&
          OU_CLASS_ERROR,OU_MODE_STD,"matfil_imposeFeastMirrorBC")
      call sys_halt()
    end if
    
    if (rmatrix%bcolumnsSorted .or. rmatrix%browsSorted) then
      call output_line ("Sorted matrices not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"matfil_imposeFeastMirrorBC")
      call sys_halt()
    end if

    if (rfmbcStructure%h_ImirrorDOFs .eq. ST_NOHANDLE) then
      ! No data inside of this structure.
      ! May happen if the region is not large enough to cover at least one DOF.
      return
    end if

    ! Get the matrix data
    call lsyssc_getbase_double (rmatrix,p_Da)
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rmatrix,p_Kld)

    ! Get the weight of the entries.
    ! =2 on finest level, =1.5 on level NLMAX-1,...
    !dmirrorWeight = 1.0_DP+REAL(4**rfmbcStructure%icoarseningLevel,DP)
    dmirrorWeight = 1.0_DP+1.0_DP*real(2**rfmbcStructure%icoarseningLevel,DP)

    ! Get pointers to the list of DOF`s that belong to that region and have
    ! to be tackled.
    ! p_ImirrorDOFs is a list of all DOF`s in the region.
    ! p_ImirrorDOFsClosed is a list of all DOF`s in the closure of the region.
    ! For every DOF in the region, it is neighbours have to be found in the
    ! clusure. If that is the case, the corresponding matrix entry has to be doubled.

    call storage_getbase_int(rfmbcStructure%h_ImirrorDOFs,p_ImirrorDOFs)
    call storage_getbase_int(rfmbcStructure%h_ImirrorDOFsClosed,p_ImirrorDOFsClosed)

    ! The matrix column corresponds to the DOF. For every DOF decide on
    ! whether it is on the FEAST mirror boundary component or not.
    ! If yes, double the matrix entry.

    ! Loop through the DOF`s. Each DOF gives us a matrix row to change.
    do i=1,size(p_ImirrorDOFs)

      ! Loop through the matrix row. All DOF`s in that matrix row that
      ! belong to the closed region have to be changed.
      do ia=p_Kld(p_ImirrorDOFs(i)),p_Kld(p_ImirrorDOFs(i)+1)-1
        ! Get the DOF.
        idof = p_Kcol(ia)

        ! Search the DOF in our list. Ok, that is an n^2 algorithm.
        ! It could easily replaced by an n log n algorithm using binary
        ! search since the list of DOF`s is sorted!
        ! Probably in a later implementation...
        do j=1,size(p_ImirrorDOFsClosed)
          if (p_ImirrorDOFsClosed(j) .eq. idof) then
            p_Da(ia) = dmirrorWeight * p_Da(ia)
            exit
          end if
        end do

      end do

    end do

! DOES NOT WORK AT THE MOMENT. THE SORTING OF THE COLUMNS
! MAY BE DIFFERENT FROM THE SORTING OF THE ROWS!!!
!  else
!
!    ! Ok, matrix is sorted, so we have to filter all the DOF`s through the
!    ! permutation before using them for implementing boundary conditions.
!    !
!    ! Get the permutation/inverse permutation from the matrix to renumber the columns into
!    ! the actual DOF numbers.
!    call sstrat_getUnsortedPosInfo(rmatrix%p_rsortStrategyRows,p_IpermInverse)
!    call sstrat_getSortedPosInfo(rmatrix%p_rsortStrategyRows,p_Iperm)
!
!    ! Loop through the DOF`s. Each DOF gives us a matrix row to change.
!    do i=1,size(p_ImirrorDOFs)
!
!      ! Loop through the matrix row. All DOF`s in that matrix row that
!      ! belong to our region have to be changed.
!      do ia=p_Kld(p_Iperm(p_ImirrorDOFs(i))),&
!            p_Kld(p_Iperm(p_ImirrorDOFs(i))+1)-1
!        ! Get the DOF.
!        idof = p_IpermInverse(p_Kcol(ia))
!
!        ! Search the DOF in our list. Ok, that is an n^2 algorithm.
!        ! It could easily replaced by an n log n algorithm since the list
!        ! of DOF`s is sorted!
!        do j=1,size(p_ImirrorDOFsClosed)
!          if (p_ImirrorDOFsClosed(j) .eq. idof) then
!            p_Da(ia) = dmirrorWeight * p_Da(ia)
!            exit
!          end if
!        end do
!
!      end do
!
!    end do
!
!  end if

  end subroutine

! ***************************************************************************
! Block matrix filters
! ***************************************************************************

  ! *************************************************************************
  ! Implementation of discrete boundary conditions into block matrices
  ! *************************************************************************

!<subroutine>

  subroutine matfil_discreteBC (rmatrix,rdiscreteBC)

!<description>
  ! This routine realises the `impose discrete boundary conditions to
  ! block matrix` filter.
  ! The matrix is modified either to rdiscreteBC (if specified) or
  ! to the default boundary conditions associated to the matrix
  ! (if rdiscreteBC is not specified).
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the matrix.
  ! If not specified, the default boundary conditions associated to the
  ! matrix rmatrix are imposed to the matrix.
  type(t_discreteBC), optional, intent(in), target :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block matrix where the boundary conditions should be imposed.
  type(t_matrixBlock), intent(inout),target :: rmatrix
!</inputoutput>

!</subroutine>

    integer :: iblock,jblock,i,inumEntries !,icp
    logical :: boffdiagSubmatrix
    type(t_discreteBCEntry), dimension(:), pointer :: p_RdiscreteBCEntry

    ! Imposing boundary conditions normally changes the whole matrix!
    ! Take the given BC structure or
    ! grab the boundary condition entry list from the matrix. This
    ! is a list of all discretised boundary conditions in the system.
    if (present(rdiscreteBC)) then
      p_RdiscreteBCEntry => rdiscreteBC%p_RdiscBCList
      inumEntries = rdiscreteBC%inumEntriesUsed
    else
      if (.not. associated(rmatrix%p_RdiscreteBC)) then
        ! There are no BC`s available, so we cannot do anything!
        return
      end if
      p_RdiscreteBCEntry => rmatrix%p_RdiscreteBC%p_RdiscBCList
      inumEntries = rmatrix%p_RdiscreteBC%inumEntriesUsed
    end if

    if (.not. associated(p_RdiscreteBCEntry)) return

    ! Is the matrix a "primal" matrix or is it a submatrix of another block matrix?
    boffdiagSubmatrix = rmatrix%imatrixSpec .eq. LSYSBS_MSPEC_OFFDIAGSUBMATRIX

    ! Now loop through all entries in this list:
    !DO i=1,SIZE(p_RdiscreteBCEntry)
    do i=1, inumEntries

      ! What for BC`s do we have here?
      select case (p_RdiscreteBCEntry(i)%itype)
      case (DISCBC_TPUNDEFINED)
        ! Do-nothing

      case (DISCBC_TPDIRICHLET)
        ! Dirichlet boundary conditions.
        ! On which component are they defined? The component specifies
        ! the row of the block matrix that is to be altered.
        iblock = p_RdiscreteBCEntry(i)%rdirichletBCs%icomponent

        ! Loop through this matrix row and implement the boundary conditions
        ! into the scalar submatrices.
        ! For now, this implements unit vectors into the diagonal matrices
        ! and zero-vectors into the offdiagonal matrices.
        ! Only exception: If the matrix is a submatrix of another matrix
        ! and not on the diagonal of its parent, we must replace the rows
        ! by zero vectors!
        do jblock = 1,rmatrix%nblocksPerRow
          if (lsysbl_isSubmatrixPresent(rmatrix,iblock,jblock)) then
            call matfil_imposeDirichletBC (&
                        rmatrix%RmatrixBlock(iblock,jblock), &
                        (iblock .ne. jblock) .or. boffdiagSubmatrix,&
                        p_RdiscreteBCEntry(i)%rdirichletBCs)
          end if
        end do

      case (DISCBC_TPPRESSUREDROP)
        ! Nothing to do; pressure drop BC`s are implemented only into the RHS.

      case (DISCBC_TPSLIP)
        ! Slip boundary conditions are treated like Dirichlet for all
        ! velocity components.
        ! This is a separate filter must be called manually.
        ! Therefore, there is nothing to do here.

      case (DISCBC_TPFEASTMIRROR)
        ! FEAST mirror boundary conditions.
        ! On which component are they defined? The component specifies
        ! the row of the block matrix that is to be altered.
        iblock = p_RdiscreteBCEntry(i)%rfeastMirrorBCs%icomponent

        ! Loop through this matrix row and implement the boundary conditions
        ! into the scalar submatrices.
        ! For now, this implements unit vectors into the diagonal matrices
        ! and zero-vectors into the offdiagonal matrices.
        ! Only exception: If the matrix is a submatrix of another matrix
        ! and not on the diagonal of its parent, we must replace the rows
        ! by zero vectors!
        do jblock = 1,rmatrix%nblocksPerRow
          if (lsysbl_isSubmatrixPresent(rmatrix,iblock,jblock)) then
            call matfil_imposeFeastMirrorBC (&
                        rmatrix%RmatrixBlock(iblock,jblock), &
                        (iblock .ne. jblock) .or. boffdiagSubmatrix,&
                        p_RdiscreteBCEntry(i)%rfeastMirrorBCs)
          end if
        end do

      case default
        call output_line (&
            "Unknown boundary condition"//sys_siL(p_RdiscreteBCEntry(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,"matfil_discreteBC")
        call sys_halt()

      end select
    end do

  end subroutine

  ! *************************************************************************

!<subroutine>

  subroutine matfil_discreteNLSlipBC (rmatrix,bforprec,rdiscreteBC)

!<description>
  ! Imposes nonlinear slip boundary conditions into a given matrix.
  ! The matrix is modified either to rdiscreteBC (if specified) or
  ! to the default boundary conditions associated to the matrix
  ! (if rdiscreteBC is not specified).
  !
  ! Note that slip boundary conditions are implemented by so called
  ! "nonlinear filters" (c.f. vectorfilters.f). Similarly to the vector
  ! filters, slip boundary conditions are not automatically implemented
  ! by matfil_discreteBC. To implement them, this routine must be called
  ! separately from matfil_discreteBC.
!</description>

!<input>
  ! Prepare matrix for preconditioning.
  ! =false: Slip-nodes are handled as Dirichlet. Rows in the
  !         matrix are replaced by unit vectors.
  ! =true : Standard. Prepare matrix for preconditioning.
  !         Only the off-diagonals of rows corresponding to
  !         the slip-nodes are cleared.
  !         The matrix therefore is prepared to act as
  !         Jacobi-preconditioner for all slip-nodes.
  logical, intent(in) :: bforprec

  ! OPTIONAL: boundary conditions to impose into the matrix.
  ! If not specified, the default boundary conditions associated to the
  ! matrix rmatrix are imposed to the matrix.
  type(t_discreteBC), optional, intent(in), target :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block matrix where the boundary conditions should be imposed.
  type(t_matrixBlock), intent(inout),target :: rmatrix
!</inputoutput>

!</subroutine>

    integer :: iblock,jblock,i,icp,inumentries
    logical :: boffdiagSubmatrix
    type(t_discreteBCEntry), dimension(:), pointer :: p_RdiscreteBCEntry

    nullify(p_RdiscreteBCEntry)
    
    if (present(rdiscreteBC)) then
      p_RdiscreteBCEntry => rdiscreteBC%p_RdiscBCList
      inumEntries = rdiscreteBC%inumEntriesUsed
    else
      if (.not. associated(rmatrix%p_rdiscreteBC)) then
        ! There are no BC`s available, so we cannot do anything!
        return
      end if
      ! Grab the boundary condition entry list from the matrix. This
      ! is a list of all discretised boundary conditions in the system.
      p_RdiscreteBCEntry => rmatrix%p_rdiscreteBC%p_RdiscBCList
      inumEntries = rmatrix%p_rdiscreteBC%inumEntriesUsed
    end if

    if (.not. associated(p_RdiscreteBCEntry)) return

    ! Is the matrix a "primal" matrix or is it a submatrix of another block matrix?
    boffdiagSubmatrix = rmatrix%imatrixSpec .eq. LSYSBS_MSPEC_OFFDIAGSUBMATRIX

    ! Now loop through all entries in this list:
    !DO i=1,SIZE(p_RdiscreteBCEntry)
    do i=1, inumentries

      ! Only implement slip boundary conditions.
      if (p_RdiscreteBCEntry(i)%itype .eq. DISCBC_TPSLIP) then
        ! Slip boundary conditions are treated like Dirichlet for all
        ! velocity components.

        ! Loop through all affected components to implement the BC`s.
        do icp = 1,p_RdiscreteBCEntry(i)%rslipBCs%ncomponents
          iblock = p_RdiscreteBCEntry(i)%rslipBCs%Icomponents(icp)

          do jblock = 1,rmatrix%nblocksPerRow
            if (lsysbl_isSubmatrixPresent(rmatrix,iblock,jblock)) then
              call matfil_imposeNLSlipBC (&
                          rmatrix%RmatrixBlock(iblock,jblock), &
                          (iblock .ne. jblock) .or. boffdiagSubmatrix,bforprec,&
                          p_RdiscreteBCEntry(i)%rslipBCs)
            end if
          end do

        end do
      end if

    end do

  end subroutine

  ! ***************************************************************************
  ! Implementation of discrete fictitious boundary conditions into
  ! block matrices
  ! ***************************************************************************

!<subroutine>

  subroutine matfil_discreteFBC (rmatrix,rdiscreteFBC)

!<description>
  ! This routine realises the `impose discrete fictitious boundary
  ! conditions to block matrix` filter.
  ! The matrix is modified either to rdiscreteFBC (if specified) or
  ! to the default boundary conditions associated to the matrix
  ! (if rdiscreteFBC is not specified).
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the matrix.
  ! If not specified, the default boundary conditions associated to the
  ! matrix rmatrix are imposed to the matrix.
  type(t_discreteFBC), optional, intent(in), target :: rdiscreteFBC
!</input>

!<inputoutput>
  ! The block matrix where the boundary conditions should be imposed.
  type(t_matrixBlock), intent(inout),target :: rmatrix
!</inputoutput>

!</subroutine>

    integer :: iblock,jblock,i,j,inumentries
    logical :: boffdiagSubmatrix
    type(t_discreteFBCEntry), dimension(:), pointer :: p_RdiscreteFBCEntry

    nullify(p_RdiscreteFBCEntry)

    if (present(rdiscreteFBC)) then
      p_RdiscreteFBCEntry => rdiscreteFBC%p_RdiscFBCList
      inumEntries = rdiscreteFBC%inumEntriesUsed
    else    
      if (.not. associated(rmatrix%p_rdiscreteBCfict)) then
        ! There are no BC`s available, so we cannot do anything!
        return
      end if

      ! Grab the boundary condition entry list from the matrix. This
      ! is a list of all discretised boundary conditions in the system.
      p_RdiscreteFBCEntry => rmatrix%p_rdiscreteBCfict%p_RdiscFBCList
      inumEntries = rmatrix%p_rdiscreteBCfict%inumEntriesUsed
    end if

    if (.not. associated(p_RdiscreteFBCEntry)) return

    ! Is the matrix a "primal" matrix or is it a submatrix of another block matrix?
    boffdiagSubmatrix = rmatrix%imatrixSpec .eq. LSYSBS_MSPEC_OFFDIAGSUBMATRIX

    ! Now loop through all entries in this list:
    do i=1,inumentries

      ! What for BC`s do we have here?
      select case (p_RdiscreteFBCEntry(i)%itype)
      case (DISCFBC_TPUNDEFINED)
        ! Do-nothing

      case (DISCFBC_TPDIRICHLET)
        ! Dirichlet boundary conditions.
        !
        ! Loop through all blocks where to impose the BC`s:
        do j=1,p_RdiscreteFBCEntry(i)%rdirichletFBCs%ncomponents

          iblock = p_RdiscreteFBCEntry(i)%rdirichletFBCs%Icomponents(j)

          ! Loop through this matrix row and implement the boundary conditions
          ! into the scalar submatrices.
          ! For now, this implements unit vectors into the diagonal matrices
          ! and zero-vectors into the offdiagonal matrices.
          ! Only exception: If the matrix is a submatrix of another matrix
          ! and not on the diagonal of its parent, we must replace the rows
          ! by zero vectors!
          do jblock = 1,rmatrix%nblocksPerRow
            if (lsysbl_isSubmatrixPresent(rmatrix,iblock,jblock)) then
              call matfil_imposeDirichletFBC (&
                          rmatrix%RmatrixBlock(iblock,jblock), &
                          (iblock .ne. jblock) .or. boffdiagSubmatrix,&
                          p_RdiscreteFBCEntry(i)%rdirichletFBCs)
            end if
          end do

        end do

      case default
        call output_line (&
            "Unknown boundary condition"//sys_siL(p_RdiscreteFBCEntry(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,"matfil_discreteFBC")
        call sys_halt()

      end select
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine matfil_normaliseToL20LmassMat (rmatrix,rlmassMat,irow)

!<description>
  ! Modifies a scalar matrix to impose the equation 
  ! <tex>$\int_{\Omega} v = 0$</tex> for all <tex>$v$</tex> by using a lumped
  ! mass matrix. This replaces row irow of the matrix rmatrix by the
  ! diagonal of the lumped mass matrix, which realises a weighted sum of
  ! all degrees of freedom.
  !
  ! This filter can be used e.g. to impose the condition <tex>$\int_{\Omega} v = 0$</tex>
  ! to a matrix which is passed to an external linear solver (like UMFPACK)
  ! which cannot use filtering to cope with indefiniteness!
!</description>

!<input>
  ! A lumped mass matrix in the same FEM space as rmatrix
  type(t_matrixScalar), intent(in) :: rlMassMat

  ! Row which should be replaced
  integer, intent(in) :: irow
!</input>

!<inputoutput>
  ! The matrix which is to be modified.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Da1, p_Da2
    integer, dimension(:), pointer :: p_Kld, p_Kdiagonal
    integer :: i,j

    ! Check, if matrix is not a copy of another matrix or if resize is to be enforced
    if (iand(rmatrix%imatrixSpec, LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0 .or.&
        iand(rmatrix%imatrixSpec, LSYSSC_MSPEC_CONTENTISCOPY)   .ne. 0) then
      call output_line("A copied matrix cannot be modified!",&
          OU_CLASS_ERROR,OU_MODE_STD,"matfil_normaliseToL20LmassMat")
        call sys_halt()
    end if

    if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
      call output_line("Virtually transposed matrices not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"matfil_normaliseToL20LmassMat")
      call sys_halt()
    end if

    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9)
      ! CSR format.
      if ((rmatrix%cdataType .ne. ST_DOUBLE) .or. &
          (rlMassMat%cdataType .ne. ST_DOUBLE)) then
        call output_line("Unsupported data type",&
            OU_CLASS_ERROR,OU_MODE_STD,"matfil_normaliseToL20LmassMat")
        call sys_halt()
      end if

      ! Exoand the matrix structure if necessary.
      call mmod_expandToFullRow (rmatrix,irow)

      ! Copy the diagonal entries of the lumped mass matrix
      ! to rmatrix.
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIXD)
        ! Get the data arrays
        call lsyssc_getbase_double (rlMassMat,p_Da1)
        call lsyssc_getbase_double (rmatrix,p_Da2)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)

        ! Copy the entries
        call lalg_copyVector (p_Da1,p_Da2(p_Kld(irow):),rmatrix%NEQ)

      case (LSYSSC_MATRIX9)

        ! Get the data arrays
        call lsyssc_getbase_double (rlMassMat,p_Da1)
        call lsyssc_getbase_double (rmatrix,p_Da2)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)
        call lsyssc_getbase_Kdiagonal (rlMassMat,p_Kdiagonal)

        ! Copy the entries
        j = p_Kld(irow)-1
        do i=1,rmatrix%NEQ
          p_Da2(j+i) = p_Da1(p_Kdiagonal(i))
        end do

      end select

    case default
      call output_line("Unsupported matrix format!",&
          OU_CLASS_ERROR,OU_MODE_STD,"matfil_normaliseToL20LmassMat")
    end select

  end subroutine

end module
