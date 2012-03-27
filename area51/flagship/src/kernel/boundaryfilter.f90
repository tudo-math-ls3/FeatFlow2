!##############################################################################
!# ****************************************************************************
!# <name> boundaryfilter </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!#
!# This module provides the basic routines for imposing boundary
!# conditions in the global matrix and the corresponding vector in
!# strong sense. This is realised as boundary filters which are
!# applied to the vectors/matrices after assembly.
!#
!# An alternative strategy is to prescribe boundary conditions in a
!# weak sense by including the boundary values into the bilinear
!# and/or linear form. This must be done during the assemble so that
!# only the boundary description stored in t_boundaryCondition is used
!# from this module.
!#
!#
!# The following routines are available:
!#
!# 1.) bdrf_filterMatrix = bdrf_filterMatrixScalar /
!#                         bdrf_filterMatrixBlock
!#     -> Performs matrix filtering
!#
!# 2.) bdrf_filterVectorByValue = bdrf_filterVectorScalarValue /
!#                                bdrf_filterVectorBlockValue
!#     -> Performs vector filtering by a scalar value
!#
!# 3.) bdrf_filterVectorByVector = bdrf_filterVectorScalarVector /
!#                                 bdrf_filterVectorBlockVector
!#     -> Performs vector filtering by a vector
!#
!# 4.) bdrf_filterVectorExplicit = bdrf_filterVectorScalarExplicit /
!#                                 bdrf_filterVectorBlockExplicit
!#     -> Performs vector filtering by imposing boundary values explicitly
!#
!# 5.) bdrf_filterSolution = bdrf_filterSolutionScalar /
!#                           bdrf_filterSolutionBlock /
!#                           bdrf_filterSolutionBlockScalar
!#     -> Performs explicit filtering of the solution and the defect vector
!#
!# </purpose>
!##############################################################################

module boundaryfilter

!$use omp_lib
  use basicgeometry
  use boundary
  use boundarycondaux
  use fparser
  use fsystem
  use genoutput
  use io
  use linearsystemblock
  use linearsystemscalar
  use matrixmodification
  use storage
  use triangulation

  implicit none

  private
  public :: bdrf_filterMatrix
  public :: bdrf_filterVectorByValue
  public :: bdrf_filterVectorByVector
  public :: bdrf_filterVectorExplicit
  public :: bdrf_filterSolution

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface bdrf_filterMatrix
    module procedure bdrf_filterMatrixScalar
    module procedure bdrf_filterMatrixBlock
  end interface

  interface bdrf_filterVectorByValue
    module procedure bdrf_filterVectorScalarValue
    module procedure bdrf_filterVectorBlockValue
  end interface

  interface bdrf_filterVectorByVector
    module procedure bdrf_filterVectorScalarVector
    module procedure bdrf_filterVectorBlockVector
  end interface

  interface bdrf_filterVectorExplicit
    module procedure bdrf_filterVectorScalarExplicit
    module procedure bdrf_filterVectorBlockExplicit
  end interface

  interface bdrf_filterSolution
    module procedure bdrf_filterSolutionScalar
    module procedure bdrf_filterSolutionBlock
    module procedure bdrf_filterSolutionBlockScalar
  end interface

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

contains

  ! *****************************************************************************

!<subroutine>

  subroutine bdrf_filterMatrixBlock(rboundaryCondition, rmatrix,&
      dvalue, rtriangulation)

!<description>
    ! This subroutine modifies the matrix entries of a block matrix.
    ! All off-diagonal entries are nullified. If required by the type
    ! of boundary condition, the diagonal entries are replaced by ones
    ! or by the value DVALUE if it is present. Otherwise, the diagonal
    ! entries are kept unmodified.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! OPTIONAL: value for diagonal entry
    real(DP), intent(in), optional :: dvalue

    ! OPTIONAL: The triangulation
    type(t_triangulation), intent(in), optional :: rtriangulation
!</input>

!<inputoutput>
    type(t_matrixBlock), intent(inout) :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iblock,jblock

    ! Check if there are strong boundary conditions
    if (.not.rboundaryCondition%bStrongBdrCond) return

    ! Loop over all blocks
    do iblock = 1, rmatrix%nblocksPerCol
      do jblock = 1, rmatrix%nblocksPerRow

        if (iblock .eq. jblock) then

          ! Diagonal block: impose the prescribed value if required
          call bdrf_filterMatrixScalar(rboundaryCondition,&
              rmatrix%RmatrixBlock(iblock,jblock), dvalue, rtriangulation)

        else

          ! Off-diagonal block: impose zero values
          if (rmatrix%RmatrixBlock(iblock,jblock)%cmatrixFormat .ne.&
              LSYSSC_MATRIXUNDEFINED) then
            call bdrf_filterMatrixScalar(rboundaryCondition,&
                rmatrix%RmatrixBlock(iblock,jblock), 0.0_DP, rtriangulation)
          end if

        end if
      end do
    end do
  end subroutine bdrf_filterMatrixBlock

  ! ***************************************************************************

!<subroutine>

  subroutine bdrf_filterMatrixScalar(rboundaryCondition, rmatrix,&
      dvalue, rtriangulation)

!<description>
    ! This subroutine modifies the matrix entries of a scalar matrix.
    ! All off-diagonal entries are nullified. If required by the type
    ! of boundary condition, the diagonal entries are replaced by ones
    ! or by the value DVALUE if it is present. Otherwise, the diagonal
    ! entries are kept unmodified.
!</description>

!<input>
    ! The boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! OPTIONAL: value for diagonal entry
    real(DP), intent(in), optional :: dvalue

    ! OPTIONAL triangulation
    type(t_triangulation), intent(in), optional, target :: rtriangulation
!</input>

!<inputoutput>
    type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:), pointer :: p_DA,p_DmaxPar,p_DvertexParameterValue
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal
    integer, dimension(:), pointer :: p_IboundaryCpIdx, p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IbdrCondCpIdx, p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCompPeriodic, p_IbdrCondPeriodic
    logical, dimension(:), pointer :: p_BisSegClosed
    logical :: bisSorted

    ! Check if there are strong boundary conditions
    if (.not.rboundaryCondition%bStrongBdrCond) return

    ! Check if matrix is sorted?
    if (rmatrix%isortStrategy .gt. 0) then
      bisSorted = .true.
      call lsyssc_unsortMatrix(rmatrix, .true.)
    else
      bisSorted = .false.
    end if

    ! Get underlying triangulation structure
      if (present(rtriangulation)) then
        p_rtriangulation => rtriangulation
      else
        if (.not.associated(rmatrix%p_rspatialDiscrTrial)) then
          call output_line('No discretisation associated!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterMatrixScalar')
          call sys_halt()
        end if

        if (.not.associated(rmatrix%p_rspatialDiscrTrial%p_rtriangulation)) then
          call output_line('No triangulation associated!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterMatrixScalar')
          call sys_halt()
        end if
        p_rtriangulation => rmatrix%p_rspatialDiscrTrial%p_rtriangulation
      end if

      ! Check spatial dimensions
      if (p_rtriangulation%ndim .ne. rboundaryCondition%ndimension) then
        call output_line('Spatial dimension mismatch!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterMatrixScalar')
        call sys_halt()
      end if

    ! How many spatial dimensions do we have?
    select case (rboundaryCondition%ndimension)
    case (NDIM1D)
      ! Set pointers for triangulation
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
      call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic, p_IbdrCompPeriodic)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic, p_IbdrCondPeriodic)

      ! What kind of matrix are we?
      select case (rmatrix%cmatrixFormat)

      case (LSYSSC_MATRIXD)
        call lsyssc_getbase_double(rmatrix, p_DA)

        call filtermatrix_MatD_1D(p_IbdrCondType, p_IbdrCondCpIdx,&
            rboundaryCondition%iboundarycount, p_IverticesAtBoundary,&
            p_IboundaryCpIdx, p_DA, dvalue)

      case (LSYSSC_MATRIX7)
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        call lsyssc_getbase_double(rmatrix, p_DA)

        call filtermatrix_Mat79_1D(p_IbdrCompPeriodic,&
            p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
            rboundaryCondition%iboundarycount, p_IverticesAtBoundary,&
            p_IboundaryCpIdx, p_Kld, p_Kcol, p_Kld, p_DA, dvalue)

      case (LSYSSC_MATRIX9)
        call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        call lsyssc_getbase_double(rmatrix, p_DA)

        call filtermatrix_Mat79_1D(p_IbdrCompPeriodic,&
            p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
            rboundaryCondition%iboundarycount, p_IverticesAtBoundary,&
            p_IboundaryCpIdx, p_Kld, p_Kcol, p_Kdiagonal, p_DA, dvalue)

      case (LSYSSC_MATRIX7INTL)

        select case (rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIXD)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          call filtermatrix_Mat79IntlD_1D(p_IbdrCompPeriodic,&
              p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
              rboundaryCondition%iboundarycount, p_IverticesAtBoundary,&
              p_IboundaryCpIdx, p_Kld, p_Kcol, p_Kld, rmatrix%NVAR, p_DA, dvalue)

        case (LSYSSC_MATRIX1)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          call filtermatrix_Mat79Intl1_1D(p_IbdrCompPeriodic,&
              p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
              rboundaryCondition%iboundarycount, p_IverticesAtBoundary,&
              p_IboundaryCpIdx, p_Kld, p_Kcol, p_Kld, rmatrix%NVAR, p_DA, dvalue)

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterMatrixScalar')
          call sys_halt()
        end select

      case (LSYSSC_MATRIX9INTL)

        select case (rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIXD)
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          call filtermatrix_Mat79IntlD_1D(p_IbdrCompPeriodic,&
              p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
              rboundaryCondition%iboundarycount, p_IverticesAtBoundary,&
              p_IboundaryCpIdx, p_Kld, p_Kcol, p_Kdiagonal, rmatrix%NVAR, p_DA, dvalue)

        case (LSYSSC_MATRIX1)
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          call filtermatrix_Mat79Intl1_1D(p_IbdrCompPeriodic,&
              p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
              rboundaryCondition%iboundarycount, p_IverticesAtBoundary,&
              p_IboundaryCpIdx, p_Kld, p_Kcol, p_Kdiagonal, rmatrix%NVAR, p_DA, dvalue)

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterMatrixScalar')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterMatrixScalar')
        call sys_halt()
      end select


    case (NDIM2D)
      ! Set pointers for triangulation
      call storage_getbase_double(&
          p_rtriangulation%h_DvertexParameterValue, p_DvertexParameterValue)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_double(rboundaryCondition%h_DmaxPar, p_DmaxPar)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
      call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic, p_IbdrCompPeriodic)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic, p_IbdrCondPeriodic)
      call storage_getbase_logical(rboundaryCondition%h_BisSegClosed, p_BisSegClosed)

      ! What kind of matrix are we?
      select case (rmatrix%cmatrixFormat)

      case (LSYSSC_MATRIXD)
        call lsyssc_getbase_double(rmatrix, p_DA)

        call filtermatrix_MatD_2D(p_IbdrCondType, p_IbdrCondCpIdx,&
            p_DmaxPar, p_BisSegClosed, rboundaryCondition&
            %iboundarycount, p_IverticesAtBoundary, p_IboundaryCpIdx,&
            p_DvertexParameterValue, p_DA, dvalue)

      case (LSYSSC_MATRIX7)
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        call lsyssc_getbase_double(rmatrix, p_DA)

        call filtermatrix_Mat79_2D(p_IbdrCompPeriodic,&
            p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
            p_DmaxPar, p_BisSegClosed, rboundaryCondition%iboundarycount,&
            p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
            p_Kld, p_Kcol, p_Kld, p_DA, dvalue)

      case (LSYSSC_MATRIX9)
        call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        call lsyssc_getbase_double(rmatrix, p_DA)

        call filtermatrix_Mat79_2D(p_IbdrCompPeriodic,&
            p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
            p_DmaxPar, p_BisSegClosed, rboundaryCondition%iboundarycount,&
            p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
            p_Kld, p_Kcol, p_Kdiagonal, p_DA, dvalue)

      case (LSYSSC_MATRIX7INTL)

        select case (rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIXD)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          call filtermatrix_Mat79IntlD_2D(p_IbdrCompPeriodic,&
              p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
              p_DmaxPar, p_BisSegClosed, rboundaryCondition%iboundarycount,&
              p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
              p_Kld, p_Kcol, p_Kld, rmatrix%NVAR, p_DA, dvalue)

        case (LSYSSC_MATRIX1)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          call filtermatrix_Mat79Intl1_2D(p_IbdrCompPeriodic,&
              p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
              p_DmaxPar, p_BisSegClosed, rboundaryCondition%iboundarycount,&
              p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
              p_Kld, p_Kcol, p_Kld, rmatrix%NVAR, p_DA, dvalue)

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterMatrixScalar')
          call sys_halt()
        end select

      case (LSYSSC_MATRIX9INTL)

        select case (rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIXD)
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          call filtermatrix_Mat79IntlD_2D(p_IbdrCompPeriodic,&
              p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
              p_DmaxPar, p_BisSegClosed, rboundaryCondition%iboundarycount,&
              p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
              p_Kld, p_Kcol, p_Kdiagonal, rmatrix%NVAR, p_DA, dvalue)

        case (LSYSSC_MATRIX1)
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          call filtermatrix_Mat79Intl1_2D(p_IbdrCompPeriodic,&
              p_IbdrCondPeriodic, p_IbdrCondType, p_IbdrCondCpIdx,&
              p_DmaxPar, p_BisSegClosed, rboundaryCondition%iboundarycount,&
              p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
              p_Kld, p_Kcol, p_Kdiagonal, rmatrix%NVAR, p_DA, dvalue)

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterMatrixScalar')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterMatrixScalar')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterMatrixScalar')
      call sys_halt()
    end select

    ! Do we have to re-sort matrix?
    if (bisSorted) call lsyssc_sortMatrix(rmatrix, .true., rmatrix%isortStrategy)

  contains

    ! Here are the real filter routines

    !***************************************************************
    ! Filter matrix and replace the diagonal entry by a scalar value.
    ! Here, the matrix is given as diagonal matrix in 1D.

    subroutine filtermatrix_MatD_1D(IbdrCondType, IbdrCondCpIdx,&
        nbct, IverticesAtBoundary, IboundaryCpIdx, DA, dvalue)

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Matrix data array
      real(DP), dimension(:), intent(inout) :: DA

      ! OPTIONAL: dilter value
      real(DP), intent(in), optional :: dvalue


      ! local variables
      integer :: ivbd,ivt,ibct,isegment


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          call output_line('Unable to handle periodic boundary ' // &
              'conditions for diagonal matrix!',&
              OU_CLASS_WARNING, OU_MODE_STD, 'filtermatrix_MatD_1D')
          call sys_halt()

        case default
          if (present(dvalue)) then
            ! Replace column by [0,...,0,dvalue,0,...,0]
            ivt     = IverticesAtBoundary(ivbd)
            DA(ivt) = dvalue
          end if
            ! Otherwise, replace column by [0,...,0,diag(A),0,...,0]
            ! That is easy, since we are a diagonal matrix.

        end select
      end do
    end subroutine filtermatrix_MatD_1D


    !***************************************************************
    ! Filter matrix and replace the diagonal entry by a scalar value.
    ! Here, the matrix is given as diagonal matrix in 2D.

    subroutine filtermatrix_MatD_2D(IbdrCondType, IbdrCondCpIdx,&
        DmaxParam, BisSegClosed, nbct, IverticesAtBoundary,&
        IboundaryCpIdx, DvertexParameterValue, DA, dvalue)

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Matrix data array
      real(DP), dimension(:), intent(inout) :: DA

      ! OPTIONAL: filter value
      real(DP), intent(in), optional :: dvalue


      ! local variables
      real(DP) :: dminValue,dmaxValue,dVertexParameterPeriodic
      integer :: ivbd,ivbdFirst,ivbdLast,ivt,ibct,isegment

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParametervalue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue =&
            nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle

          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

          case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
            call output_line('Unable to handle periodic boundary' // &
                'conditions for diagonal matrix!',&
                OU_CLASS_WARNING ,OU_MODE_STD,'filtermatrix_MatD_2D')
            call sys_halt()

          case default
            if (present(dvalue)) then
              ! Replace column by [0,...,0,dvalue,0,...,0]
              ivt     = IverticesAtBoundary(ivbd)
              DA(ivt) = dvalue
            end if
              ! Otherwise, replace column by [0,...,0,diag(A),0,...,0]
              ! That is easy, since we are a diagonal matrix.

          end select
        end do
      end do
    end subroutine filtermatrix_MatD_2D


    !***************************************************************
    ! Filter matrix and replace the diagonal entry by a scalar value.
    ! Here, the matrix is given in format CSR7 or CSR9 format in 1D.

    subroutine filtermatrix_Mat79_1D(IbdrCompPeriodic,&
        IbdrCondPeriodic, IbdrCondType, IbdrCondCpIdx, nbct,&
        IverticesAtBoundary, IboundaryCpIdx, Kld, Kcol, Kdiagonal,&
        DA, dvalue)

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Matrix data array
      real(DP), dimension(:), intent(inout) :: DA

      ! OPTIONAL: filter value
      real(DP), intent(in), optional :: dvalue


      ! local variables
      integer :: ibeg,iend,ipos,jpos,idiag,jdiag
      integer :: ivbd,ivbdPeriodic,ivt,ivtPeriodic
      integer :: ibct,ibctPeriodic,isegment,iaux

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          ! Compute vertex parameter value at periodic boundary
          ibctPeriodic = IbdrCompPeriodic(isegment)
          ivbdPeriodic = IboundaryCpIdx(ibctPeriodic)

          ! Get numbers of vertices
          ivt         = IverticesAtBoundary(ivbd)
          ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

          ! Add entries from the row that corresponds to node ivt to
          ! the row that corresponds to node ivtPeriodic and set
          ! DVALUE and -DVALUE in row for node ivt (if present)
          ipos  = Kld(ivt)
          jpos  = Kld(ivtPeriodic)
          idiag = Kdiagonal(ivt)
          jdiag = Kdiagonal(ivtPeriodic)

          do while(jpos .le. Kld(ivtPeriodic+1)-1)
            if (Kcol(jpos) .eq. ivtPeriodic) then
              ! Skip the diagonal of the second row since it will be
              ! processed below
              jpos = jpos+1
            elseif (Kcol(jpos) .eq. ivt) then
              ! Add diagonal entry of the first row to the
              ! corresponding entry in the second row and set the
              ! diagonal entry of the first row to DVALUE (if present).
              Da(jpos) = Da(jpos) + Da(idiag)
              if (present(dvalue)) then
                Da(idiag) = dvalue
              else
                Da(idiag) = 1.0_DP
              end if
              jpos = jpos+1
            elseif(Kcol(ipos) .eq. ivt) then
              ! Skip the diagonal of the first row since it will be
              ! processed below
              ipos = ipos+1
            elseif (Kcol(ipos) .eq. ivtPeriodic) then
              ! Add entry of the first row to the diagonal of the
              ! second row and set the corresponding entry in the
              ! first row to -DVALUE (if present).
              Da(jdiag) = Da(jdiag) + Da(ipos)
              if (present(dvalue)) then
                Da(ipos) = -dvalue
              else
                Da(ipos) = -1.0_DP
              end if
              ipos = ipos+1
            else
              ! Add entry of the first row to the second row and
              ! nullify first one
              Da(jpos) = Da(jpos)+Da(ipos)
              Da(ipos) = 0.0_DP
              ipos     = ipos+1
              jpos     = jpos+1
            end if
          end do

        case default
          ivt  = IverticesAtBoundary(ivbd)
          ibeg = Kld(ivt)
          iend = Kld(ivt+1)-1
          ipos = Kdiagonal(ivt)

          if (present(dvalue)) then
            ! Replace column by [0,...,0,dvalue,0,...,0]
            DA(ibeg:iend) = 0.0_DP
            DA(ipos)      = dvalue
          else
            ! Replace column by [0,...,0,diag(A),0,...,0]
            forall (iaux=ibeg:iend, iaux .ne. ipos)
              DA(iaux) = 0.0_DP
            end forall
          end if

        end select
      end do
    end subroutine filtermatrix_Mat79_1D


    !***************************************************************
    ! Filter matrix and replace the diagonal entry by a scalar value.
    ! Here, the matrix is given in format CSR7 or CSR9 format in 2D.

    subroutine filtermatrix_Mat79_2D(IbdrCompPeriodic,&
        IbdrCondPeriodic, IbdrCondType, IbdrCondCpIdx, DmaxParam,&
        BisSegClosed, nbct, IverticesAtBoundary, IboundaryCpIdx,&
        DvertexParameterValue, Kld, Kcol, Kdiagonal, DA, dvalue)

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Matrix data array
      real(DP), dimension(:), intent(inout) :: DA

      ! OPTIONAL: filter value
      real(DP), intent(in), optional :: dvalue


      ! local variables
      real(DP) :: dminValue,dmaxValue,dVertexParameterPeriodic
      integer :: ibeg,iend,ipos,jpos,idiag,jdiag,iaux
      integer :: ivbd,ivbdFirst,ivbdLast,ivbdPeriodic,ivt,ivtPeriodic
      integer :: ibct,ibctPeriodic,isegment,isegmentPeriodic

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue =&
            nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle

          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))
            
          case (BDRC_PERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            if (isegmentPeriodic .eq. IbdrCondCpIdx(ibctPeriodic)) then
              dVertexParameterPeriodic = DmaxParam(isegment)-&
                  DvertexParameterValue(ivbd)
            else
              dVertexParameterPeriodic = DmaxParam(isegmentPeriodic-1)+&
                  DmaxParam(isegment)-DvertexParameterValue(ivbd)
            end if

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add entries from the row that corresponds to node ivt to
            ! the row that corresponds to node ivtPeriodic and set
            ! DVALUE and -DVALUE in row for node ivt (if present)
            ipos  = Kld(ivt)
            jpos  = Kld(ivtPeriodic)
            idiag = Kdiagonal(ivt)
            jdiag = Kdiagonal(ivtPeriodic)

            do while(jpos .le. Kld(ivtPeriodic+1)-1)
              if (Kcol(jpos) .eq. ivtPeriodic) then
                ! Skip the diagonal of the second row since it will be
                ! processed below
                jpos = jpos+1
              elseif (Kcol(jpos) .eq. ivt) then
                ! Add diagonal entry of the first row to the
                ! corresponding entry in the second row and set the
                ! diagonal entry of the first row to DVALUE (if present)
                Da(jpos) = Da(jpos) + Da(idiag)
                if (present(dvalue)) then
                  Da(idiag) = dvalue
                else
                  Da(idiag) = 1.0_DP
                end if
                jpos = jpos+1
              elseif(Kcol(ipos) .eq. ivt) then
                ! Skip the diagonal of the first row since it will be
                ! processed below
                ipos = ipos+1
              elseif (Kcol(ipos) .eq. ivtPeriodic) then
                ! Add entry of the first row to the diagonal of the
                ! second row and set the corresponding entry in the
                ! first row to -DVALUE (if present).
                Da(jdiag) = Da(jdiag) + Da(ipos)
                if (present(dvalue)) then
                  Da(ipos) = -dvalue
                else
                  Da(ipos) = -1.0_DP
                end if
                ipos = ipos+1
              else
                ! Add entry of the first row to the second row and
                ! nullify first one
                Da(jpos) = Da(jpos)+Da(ipos)
                Da(ipos) = 0.0_DP
                ipos     = ipos+1
                jpos     = jpos+1
              end if
            end do

          case (BDRC_ANTIPERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            dVertexParameterPeriodic = p_DmaxPar(isegmentPeriodic)-&
                (p_DmaxPar(isegment)-p_DvertexParameterValue(ivbd))

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add entries from the row that corresponds to node ivt to
            ! the row that corresponds to node ivtPeriodic and set
            ! DVALUE and -DVALUE in row for node ivt (if present)
            ipos  = Kld(ivt)
            jpos  = Kld(ivtPeriodic)
            idiag = Kdiagonal(ivt)
            jdiag = Kdiagonal(ivtPeriodic)

            do while(jpos .le. Kld(ivtPeriodic+1)-1)
              if (Kcol(jpos) .eq. ivtPeriodic) then
                ! Skip the diagonal of the second row since it will be
                ! processed below
                jpos = jpos+1
              elseif (Kcol(jpos) .eq. ivt) then
                ! Add diagonal entry of the first row to the
                ! corresponding entry in the second row and set the
                ! diagonal entry of the first row to DVALUE
                Da(jpos) = Da(jpos) + Da(idiag)
                if (present(dvalue)) then
                  Da(idiag) = dvalue
                else
                  Da(idiag) = 1.0_DP
                end if
                jpos = jpos+1
              elseif(Kcol(ipos) .eq. ivt) then
                ! Skip the diagonal of the first row since it will be
                ! processed below
                ipos = ipos+1
              elseif (Kcol(ipos) .eq. ivtPeriodic) then
                ! Add entry of the first row to the diagonal of the
                ! second row and set the corresponding entry in the
                ! first row to -DVALUE
                Da(jdiag) = Da(jdiag) + Da(ipos)
                if (present(dvalue)) then
                  Da(ipos) = -dvalue
                else
                  Da(ipos) = -1.0_DP
                end if
                ipos = ipos+1
              else
                ! Add entry of the first row to the second row and
                ! nullify first one
                Da(jpos) = Da(jpos)+Da(ipos)
                Da(ipos) = 0.0_DP
                ipos     = ipos+1
                jpos     = jpos+1
              end if
            end do

          case default
            ivt  = IverticesAtBoundary(ivbd)
            ibeg = Kld(ivt)
            iend = Kld(ivt+1)-1
            ipos = Kdiagonal(ivt)

            if (present(dvalue)) then
              ! Replace column by [0,...,0,dvalue,0,...,0]
              DA(ibeg:iend) = 0.0_DP
              DA(ipos)      = dvalue
            else
              ! Replace column by [0,...,0,diag(A),0,...,0]
              forall (iaux=ibeg:iend, iaux .ne. ipos)
                DA(iaux) = 0.0_DP
              end forall
            end if

          end select
        end do
      end do
    end subroutine filtermatrix_Mat79_2D


    !***************************************************************
    ! Filter matrix and replace the diagonal entry by a scalar value.
    ! Here, the matrix is given in format CSR7 or CSR9 interleave format,
    ! whereby each local matrix is a diagonal matrix in 1D.

    subroutine filtermatrix_Mat79IntlD_1D(IbdrCompPeriodic,&
        IbdrCondPeriodic, IbdrCondType, IbdrCondCpIdx, nbct,&
        IverticesAtBoundary, IboundaryCpIdx, Kld, Kcol, Kdiagonal,&
        nvar, DA, dvalue)

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Number of variables
      integer, intent(in) :: nvar

      ! Matrix data array
      real(DP), dimension(nvar,*), intent(inout) :: DA

      ! OPTIONAL: filter value
      real(DP), intent(in), optional :: dvalue


      ! local variables
      integer :: ibeg,iend,ipos,jpos,iaux,idiag,jdiag
      integer :: ivbd,ivbdFirst,ivbdLast,ivbdPeriodic,ivt,ivtPeriodic
      integer :: ibct,ibctPeriodic,isegment

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))
          
        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          ! Compute vertex parameter value at periodic boundary
          ibctPeriodic = IbdrCompPeriodic(isegment)
          ivbdPeriodic = IboundaryCpIdx(ibctPeriodic)

          ! Get numbers of vertices
          ivt         = IverticesAtBoundary(ivbd)
          ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

          ! Add entries from the row that corresponds to node ivt to
          ! the row that corresponds to node ivtPeriodic and set
          ! DVALUE and -DVALUE in row for node ivt (if present).
          ipos  = Kld(ivt)
          jpos  = Kld(ivtPeriodic)
          idiag = Kdiagonal(ivt)
          jdiag = Kdiagonal(ivtPeriodic)

          do while(jpos .le. Kld(ivtPeriodic+1)-1)
            if (Kcol(jpos) .eq. ivtPeriodic) then
              ! Skip the diagonal of the second row since it will be
              ! processed below
              jpos = jpos+1
            elseif (Kcol(jpos) .eq. ivt) then
              ! Add diagonal entry of the first row to the
              ! corresponding entry in the second row and set the
              ! diagonal entry of the first row to DVALUE
              Da(:,jpos) = Da(:,jpos) + Da(:,idiag)
              if (present(dvalue)) then
                Da(:,idiag) = dvalue
              else
                Da(:,idiag) = 1.0_DP
              end if
              jpos = jpos+1
            elseif(Kcol(ipos) .eq. ivt) then
              ! Skip the diagonal of the first row since it will be
              ! processed below
              ipos = ipos+1
            elseif (Kcol(ipos) .eq. ivtPeriodic) then
              ! Add entry of the first row to the diagonal of the
              ! second row and set the corresponding entry in the
              ! first row to -DVALUE (if present).
              Da(:,jdiag) = Da(:,jdiag) + Da(:,ipos)
              if (present(dvalue)) then
                Da(:,ipos) = -dvalue
              else
                Da(:,ipos) = -1.0_DP
              end if
              ipos = ipos+1
            else
              ! Add entry of the first row to the second row and
              ! nullify first one
              Da(:,jpos) = Da(:,jpos)+Da(:,ipos)
              Da(:,ipos) = 0.0_DP
              ipos       = ipos+1
              jpos       = jpos+1
            end if
          end do

        case default
          ivt  = IverticesAtBoundary(ivbd)
          ibeg = Kld(ivt)
          iend = Kld(ivt+1)-1
          ipos = Kdiagonal(ivt)

          if (present(dvalue)) then
            ! Replace all columns by [0,...,0,dvalue,0,...,0]
            do iaux = ibeg, iend
              DA(:,iaux) = 0.0_DP
            end do
            DA(:,ipos) = dvalue
          else
            ! Replace all columns by [0,...,0,diag(A),0,...,0]
            forall (iaux=ibeg:iend, iaux .ne. ipos)
              DA(:,iaux) = 0.0_DP
            end forall
          end if

        end select
      end do
    end subroutine filtermatrix_Mat79IntlD_1D


    !***************************************************************
    ! Filter matrix and replace the diagonal entry by a scalar value.
    ! Here, the matrix is given in format CSR7 or CSR9 interleave format,
    ! whereby each local matrix is a diagonal matrix in 2D.

    subroutine filtermatrix_Mat79IntlD_2D(IbdrCompPeriodic,&
        IbdrCondPeriodic, IbdrCondType, IbdrCondCpIdx, DmaxParam,&
        BisSegClosed, nbct, IverticesAtBoundary, IboundaryCpIdx,&
        DvertexParameterValue, Kld, Kcol, Kdiagonal, nvar, DA, dvalue)

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Number of variables
      integer, intent(in) :: nvar

      ! Matrix data array
      real(DP), dimension(nvar,*), intent(inout) :: DA

      ! OPTIONAL: filter value
      real(DP), intent(in), optional :: dvalue


      ! local variables
      real(DP) :: dminValue,dmaxValue,dVertexParameterPeriodic
      integer :: ibeg,iend,ipos,jpos,iaux,idiag,jdiag
      integer :: ivbd,ivbdFirst,ivbdLast,ivbdPeriodic,ivt,ivtPeriodic
      integer :: ibct,ibctPeriodic,isegment,isegmentPeriodic

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue =&
            nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle

          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

          case (BDRC_PERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            if (isegmentPeriodic .eq. IbdrCondCpIdx(ibctPeriodic)) then
              dVertexParameterPeriodic = DmaxParam(isegment)-&
                  DvertexParameterValue(ivbd)
            else
              dVertexParameterPeriodic = DmaxParam(isegmentPeriodic-1)+&
                  DmaxParam(isegment)-DvertexParameterValue(ivbd)
            end if

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add entries from the row that corresponds to node ivt to
            ! the row that corresponds to node ivtPeriodic and set
            ! DVALUE and -DVALUE in row for node ivt (if present).
            ipos  = Kld(ivt)
            jpos  = Kld(ivtPeriodic)
            idiag = Kdiagonal(ivt)
            jdiag = Kdiagonal(ivtPeriodic)

            do while(jpos .le. Kld(ivtPeriodic+1)-1)
              if (Kcol(jpos) .eq. ivtPeriodic) then
                ! Skip the diagonal of the second row since it will be
                ! processed below
                jpos = jpos+1
              elseif (Kcol(jpos) .eq. ivt) then
                ! Add diagonal entry of the first row to the
                ! corresponding entry in the second row and set the
                ! diagonal entry of the first row to DVALUE (if present).
                Da(:,jpos) = Da(:,jpos) + Da(:,idiag)
                if (present(dvalue)) then
                  Da(:,idiag) = dvalue
                else
                  Da(:,idiag) = 1.0_DP
                end if
                jpos = jpos+1
              elseif(Kcol(ipos) .eq. ivt) then
                ! Skip the diagonal of the first row since it will be
                ! processed below
                ipos = ipos+1
              elseif (Kcol(ipos) .eq. ivtPeriodic) then
                ! Add entry of the first row to the diagonal of the
                ! second row and set the corresponding entry in the
                ! first row to -DVALUE (if present).
                Da(:,jdiag) = Da(:,jdiag) + Da(:,ipos)
                if (present(dvalue)) then
                  Da(:,ipos) = -dvalue
                else
                  Da(:,ipos) = -1.0_DP
                end if
                ipos = ipos+1
              else
                ! Add entry of the first row to the second row and
                ! nullify first one
                Da(:,jpos) = Da(:,jpos)+Da(:,ipos)
                Da(:,ipos) = 0.0_DP
                ipos       = ipos+1
                jpos       = jpos+1
              end if
            end do

          case (BDRC_ANTIPERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            dVertexParameterPeriodic = p_DmaxPar(isegmentPeriodic)-&
                (p_DmaxPar(isegment)-p_DvertexParameterValue(ivbd))

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add entries from the row that corresponds to node ivt to
            ! the row that corresponds to node ivtPeriodic and set
            ! DVALUE and -DVALUE in row for node ivt (if present).
            ipos  = Kld(ivt)
            jpos  = Kld(ivtPeriodic)
            idiag = Kdiagonal(ivt)
            jdiag = Kdiagonal(ivtPeriodic)

            do while(jpos .le. Kld(ivtPeriodic+1)-1)
              if (Kcol(jpos) .eq. ivtPeriodic) then
                ! Skip the diagonal of the second row since it will be
                ! processed below
                jpos = jpos+1
              elseif (Kcol(jpos) .eq. ivt) then
                ! Add diagonal entry of the first row to the
                ! corresponding entry in the second row and set the
                ! diagonal entry of the first row to DVALUE (if present)
                Da(:,jpos) = Da(:,jpos) + Da(:,idiag)
                if (present(dvalue)) then
                  Da(:,idiag) = dvalue
                else
                  Da(:,idiag) = 1.0_DP
                end if
                jpos = jpos+1
              elseif(Kcol(ipos) .eq. ivt) then
                ! Skip the diagonal of the first row since it will be
                ! processed below
                ipos = ipos+1
              elseif (Kcol(ipos) .eq. ivtPeriodic) then
                ! Add entry of the first row to the diagonal of the
                ! second row and set the corresponding entry in the
                ! first row to -DVALUE (if present)
                Da(:,jdiag) = Da(:,jdiag) + Da(:,ipos)
                if (present(dvalue)) then
                  Da(:,ipos) = -dvalue
                else
                  Da(:,ipos) = -1.0_DP
                end if
                ipos = ipos+1
              else
                ! Add entry of the first row to the second row and
                ! nullify first one
                Da(:,jpos) = Da(:,jpos)+Da(:,ipos)
                Da(:,ipos) = 0.0_DP
                ipos       = ipos+1
                jpos       = jpos+1
              end if
            end do

          case default
            ivt  = IverticesAtBoundary(ivbd)
            ibeg = Kld(ivt)
            iend = Kld(ivt+1)-1
            ipos = Kdiagonal(ivt)
            
            if (present(dvalue)) then
              ! Replace all columns by [0,...,0,dvalue,0,...,0]
              do iaux = ibeg, iend
                DA(:,iaux) = 0.0_DP
              end do
              DA(:,ipos) = dvalue
            else
              ! Replace all columns by [0,...,0,diag(A),0,...,0]
              forall(iaux=ibeg:iend, iaux .ne. ipos)
                DA(:,iaux) = 0.0_DP
              end forall
            end if

          end select
        end do
      end do
    end subroutine filtermatrix_Mat79IntlD_2D


    !***************************************************************
    ! Filter matrix and replace the diagonal entry by a scalar value.
    ! Here, the matrix is given in format CSR7 or CSR9 interleave format,
    ! whereby each local matrix is a full matrix in 1D.

    subroutine filtermatrix_Mat79Intl1_1D(IbdrCompPeriodic,&
        IbdrCondPeriodic, IbdrCondType, IbdrCondCpIdx, nbct,&
        IverticesAtBoundary, IboundaryCpIdx, Kld, Kcol, Kdiagonal,&
        nvar, DA, dvalue)

       ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Number of variables
      integer, intent(in) :: nvar

      ! Matrix data array
      real(DP), dimension(nvar,nvar,*), intent(inout) :: DA

      ! OPTIONAL: filter value
      real(DP), intent(in), optional :: dvalue


      ! local variables
      integer :: ibeg,iend,ipos,jpos,iaux,idiag,jdiag
      integer :: ivbd,ivbdPeriodic,ivt,ivtPeriodic
      integer :: ibct,ibctPeriodic,ivar,jvar,isegment

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          ! Compute vertex parameter value at periodic boundary
          ibctPeriodic = IbdrCompPeriodic(isegment)
          ivbdPeriodic = IboundaryCpIdx(ibctPeriodic)

          ! Get numbers of vertices
          ivt         = IverticesAtBoundary(ivbd)
          ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

          ! Add entries from the row that corresponds to node ivt to
          ! the row that corresponds to node ivtPeriodic and set
          ! DVALUE and -DVALUE in row for node ivt (if present).
          ipos  = Kld(ivt)
          jpos  = Kld(ivtPeriodic)
          idiag = Kdiagonal(ivt)
          jdiag = Kdiagonal(ivtPeriodic)

          do while(jpos .le. Kld(ivtPeriodic+1)-1)
            if (Kcol(jpos) .eq. ivtPeriodic) then
              ! Skip the diagonal of the second row since it will be
              ! processed below
              jpos = jpos+1
            elseif (Kcol(jpos) .eq. ivt) then
              ! Add diagonal entry of the first row to the
              ! corresponding entry in the second row and set the
              ! diagonal entry of the first row to DVALUE
              Da(:,:,jpos) = Da(:,:,jpos) + Da(:,:,idiag)
              if (present(dvalue)) then
                Da(:,:,idiag) = dvalue
              else
                Da(:,:,idiag) = 1.0_DP
              end if
              jpos = jpos+1
            elseif(Kcol(ipos) .eq. ivt) then
              ! Skip the diagonal of the first row since it will be
              ! processed below
              ipos = ipos+1
            elseif (Kcol(ipos) .eq. ivtPeriodic) then
              ! Add entry of the first row to the diagonal of the
              ! second row and set the corresponding entry in the
              ! first row to -DVALUE (if present).
              Da(:,:,jdiag) = Da(:,:,jdiag) + Da(:,:,ipos)
              if (present(dvalue)) then
                Da(:,:,ipos) = -dvalue
              else
                Da(:,:,ipos) = -1.0_DP
              end if
              ipos = ipos+1
            else
              ! Add entry of the first row to the second row and
              ! nullify first one
              Da(:,:,jpos) = Da(:,:,jpos)+Da(:,:,ipos)
              Da(:,:,ipos) = 0.0_DP
              ipos         = ipos+1
              jpos         = jpos+1
            end if
          end do

        case default
          ivt  = IverticesAtBoundary(ivbd)
          ibeg = Kld(ivt)
          iend = Kld(ivt+1)-1
          ipos = Kdiagonal(ivt)

          if (present(dvalue)) then
            ! Replace all columns by [0,...,0,dvalue,0,...,0]
            do iaux = ibeg, iend
              DA(:,:,iaux) = 0.0_DP
            end do

            do ivar = 1, nvar
              DA(ivar,ivar,ipos) = dvalue
            end do
          else
            ! Replace all columns by [0,...,0,diag(A),0,...,0]
            forall (iaux=ibeg:iend, iaux .ne. ipos)
              DA(:,:,iaux) = 0.0_DP
            end forall
            
            forall (ivar=1:nvar, jvar=1:nvar, ivar .ne. jvar)
              DA(ivar,jvar,ipos) = 0.0_DP
            end forall
          end if

        end select
      end do
    end subroutine filtermatrix_Mat79Intl1_1D


    !***************************************************************
    ! Filter matrix and replace the diagonal entry by a scalar value.
    ! Here, the matrix is given in format CSR7 or CSR9 interleave format,
    ! whereby each local matrix is a full matrix in 2D.

    subroutine filtermatrix_Mat79Intl1_2D(IbdrCompPeriodic,&
        IbdrCondPeriodic, IbdrCondType, IbdrCondCpIdx, DmaxParam,&
        BisSegClosed, nbct, IverticesAtBoundary, IboundaryCpIdx,&
        DvertexParameterValue, Kld, Kcol, Kdiagonal, nvar, DA, dvalue)

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Number of variables
      integer, intent(in) :: nvar

      ! Matrix data array
      real(DP), dimension(nvar,nvar,*), intent(inout) :: DA

      ! OPTIONAL: filter value
      real(DP), intent(in), optional :: dvalue


      ! local variables
      real(DP) :: dminValue,dmaxValue,dVertexParameterPeriodic
      integer :: ibeg,iend,ipos,jpos,iaux,idiag,jdiag
      integer :: ivbd,ivbdFirst,ivbdLast,ivbdPeriodic,ivt,ivtPeriodic
      integer :: ibct,ibctPeriodic,ivar,jvar,isegment,isegmentPeriodic

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue =&
            nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle

          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))
            
          case (BDRC_PERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            if (isegmentPeriodic .eq. IbdrCondCpIdx(ibctPeriodic)) then
              dVertexParameterPeriodic = DmaxParam(isegment)-&
                  DvertexParameterValue(ivbd)
            else
              dVertexParameterPeriodic = DmaxParam(isegmentPeriodic-1)+&
                  DmaxParam(isegment)-DvertexParameterValue(ivbd)
            end if

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add entries from the row that corresponds to node ivt to
            ! the row that corresponds to node ivtPeriodic and set
            ! DVALUE and -DVALUE in row for node ivt (if present).
            ipos  = Kld(ivt)
            jpos  = Kld(ivtPeriodic)
            idiag = Kdiagonal(ivt)
            jdiag = Kdiagonal(ivtPeriodic)

            do while(jpos .le. Kld(ivtPeriodic+1)-1)
              if (Kcol(jpos) .eq. ivtPeriodic) then
                ! Skip the diagonal of the second row since it will be
                ! processed below
                jpos = jpos+1
              elseif (Kcol(jpos) .eq. ivt) then
                ! Add diagonal entry of the first row to the
                ! corresponding entry in the second row and set the
                ! diagonal entry of the first row to DVALUE (if present).
                Da(:,:,jpos) = Da(:,:,jpos) + Da(:,:,idiag)
                if (present(dvalue)) then
                  Da(:,:,idiag) = dvalue
                else
                  Da(:,:,idiag) = 1.0_DP
                end if
                jpos = jpos+1
              elseif(Kcol(ipos) .eq. ivt) then
                ! Skip the diagonal of the first row since it will be
                ! processed below
                ipos = ipos+1
              elseif (Kcol(ipos) .eq. ivtPeriodic) then
                ! Add entry of the first row to the diagonal of the
                ! second row and set the corresponding entry in the
                ! first row to -DVALUE (if present)
                Da(:,:,jdiag) = Da(:,:,jdiag) + Da(:,:,ipos)
                if (present(dvalue)) then
                  Da(:,:,ipos) = -dvalue
                else
                  Da(:,:,ipos) = -1.0_DP
                end if
                ipos = ipos+1
              else
                ! Add entry of the first row to the second row and nullify first one
                Da(:,:,jpos) = Da(:,:,jpos)+Da(:,:,ipos)
                Da(:,:,ipos) = 0.0_DP
                ipos         = ipos+1
                jpos         = jpos+1
              end if
            end do

          case (BDRC_ANTIPERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            dVertexParameterPeriodic = p_DmaxPar(isegmentPeriodic)-&
                (p_DmaxPar(isegment)-p_DvertexParameterValue(ivbd))

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add entries from the row that corresponds to node ivt to
            ! the row that corresponds to node ivtPeriodic and set
            ! DVALUE and -DVALUE in row for node ivt (if present)
            ipos  = Kld(ivt)
            jpos  = Kld(ivtPeriodic)
            idiag = Kdiagonal(ivt)
            jdiag = Kdiagonal(ivtPeriodic)

            do while(jpos .le. Kld(ivtPeriodic+1)-1)
              if (Kcol(jpos) .eq. ivtPeriodic) then
                ! Skip the diagonal of the second row since it will be
                ! processed below
                jpos = jpos+1
              elseif (Kcol(jpos) .eq. ivt) then
                ! Add diagonal entry of the first row to the
                ! corresponding entry in the second row and set the
                ! diagonal entry of the first row to DVALUE (if present)
                Da(:,:,jpos) = Da(:,:,jpos) + Da(:,:,idiag)
                if (present(dvalue)) then
                  Da(:,:,idiag) = dvalue
                else
                  Da(:,:,idiag) = 1.0_DP
                end if
                jpos = jpos+1
              elseif(Kcol(ipos) .eq. ivt) then
                ! Skip the diagonal of the first row since it will be
                ! processed below
                ipos = ipos+1
              elseif (Kcol(ipos) .eq. ivtPeriodic) then
                ! Add entry of the first row to the diagonal of the
                ! second row and set the corresponding entry in the
                ! first row to -DVALUE (if present)
                Da(:,:,jdiag) = Da(:,:,jdiag) + Da(:,:,ipos)
                if (present(dvalue)) then
                  Da(:,:,ipos) = -dvalue
                else
                  Da(:,:,ipos) = -1.0_DP
                end if
                ipos = ipos+1
              else
                ! Add entry of the first row to the second row and
                ! nullify first one
                Da(:,:,jpos) = Da(:,:,jpos)+Da(:,:,ipos)
                Da(:,:,ipos) = 0.0_DP
                ipos         = ipos+1
                jpos         = jpos+1
              end if
            end do

          case default
            ivt  = IverticesAtBoundary(ivbd)
            ibeg = Kld(ivt)
            iend = Kld(ivt+1)-1
            ipos = Kdiagonal(ivt)

            if (present(dvalue)) then
              ! Replace all columns by [0,...,0,dvalue,0,...,0]
              do iaux = ibeg, iend
                DA(:,:,iaux) = 0.0_DP
              end do
              
              do ivar = 1, nvar
                DA(ivar,ivar,ipos) = dvalue
              end do
            else
              ! Replace all columns by [0,...,0,diag(A),0,...,0]
              forall (iaux=ibeg:iend, iaux .ne. ipos)
                DA(:,:,iaux) = 0.0_DP
              end forall

              forall (ivar=1:nvar, jvar=1:nvar, ivar .ne. jvar)
                DA(ivar,jvar,ipos) = 0.0_DP
              end forall
            end if

          end select
        end do
      end do
    end subroutine filtermatrix_Mat79Intl1_2D
  end subroutine bdrf_filterMatrixScalar

  ! ***************************************************************************

!<subroutine>

  subroutine bdrf_filterVectorBlockValue(rboundaryCondition, rvector,&
      dvalue, rtriangulation)

!<description>
    ! This subroutine performs vector filtering by a constant value.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! The scalar filter value
    real(DP), intent(in) :: dvalue

    ! OPTIONAL: The triangulation
    type(t_triangulation), intent(in), optional :: rtriangulation
!</input>

!<inputoutput>
    ! The vector that should be filtered
    type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iblock

    ! Check if there are strong boundary conditions
    if (.not.rboundaryCondition%bStrongBdrCond) return

    ! Loop over all blocks
    do iblock = 1, rvector%nblocks
      call bdrf_filterVectorScalarValue(rboundaryCondition, rvector&
          %RvectorBlock(iblock), dvalue, rtriangulation)
    end do
  end subroutine bdrf_filterVectorBlockValue

  ! ***************************************************************************

!<subroutine>

  subroutine bdrf_filterVectorScalarValue(rboundaryCondition, rvector,&
      dvalue, rtriangulation)

!<description>
    ! This subroutine performs vector filtering by a constant value.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! The scalar filter value
    real(DP), intent(in) :: dvalue

    ! OPTIONAL: The triangulation
    type(t_triangulation), intent(in), optional, target :: rtriangulation
!</input>

!<inputoutput>
    ! The vector that should be filtered
    type(t_vectorScalar), intent(inout) :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:), pointer :: p_Dx
    real(DP), dimension(:), pointer :: p_DmaxPar
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer, dimension(:), pointer :: p_IbdrCondType
    logical, dimension(:), pointer :: p_BisSegClosed
    logical :: bisSorted

    ! Check if there are strong boundary conditions
    if (.not.rboundaryCondition%bStrongBdrCond) return

    ! Check if vector is sorted?
    if (rvector%isortStrategy .gt. 0) then
      bisSorted = .true.
      call lsyssc_vectorActivateSorting(rvector, .false.)
    else
      bisSorted = .false.
    end if

    ! Set pointer for vector
    call lsyssc_getbase_double (rvector, p_Dx)

    ! Get underlying triangulation structure
    if (present(rtriangulation)) then
      p_rtriangulation => rtriangulation
    else
      if (.not.associated(rvector%p_rspatialDiscr)) then
        call output_line('No discretisation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarValue')
        call sys_halt()
      end if

      if (.not.associated(rvector%p_rspatialDiscr%p_rtriangulation)) then
        call output_line('No triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarValue')
        call sys_halt()
      end if
      p_rtriangulation => rvector%p_rspatialDiscr%p_rtriangulation
    end if

    ! Check spatial dimensions
    if (p_rtriangulation%ndim .ne. rboundaryCondition%ndimension) then
      call output_line('Spatial dimension mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarValue')
      call sys_halt()
    end if

    ! How many spatial dimensions do we have?
    select case (rboundaryCondition%ndimension)
    case (NDIM1D)
      ! Set pointers for triangulation
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)

      ! Set prescribed boundary values in 1D
      call filtervector_1D(p_IbdrCondType, p_IbdrCondCpIdx,&
          rboundaryCondition%iboundarycount, p_IverticesAtBoundary,&
          p_IboundaryCpIdx, rvector%NVAR, p_Dx, dvalue)

    case (NDIM2D)
      ! Set pointers for triangulation
      call storage_getbase_double(&
          p_rtriangulation%h_DvertexParameterValue, p_DvertexParameterValue)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_double(rboundaryCondition%h_DmaxPar, p_DmaxPar)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
      call storage_getbase_logical(rboundaryCondition%h_BisSegClosed, p_BisSegClosed)

      ! Set prescribed boundary values in 2D
      call filtervector_2D(p_IbdrCondType, p_IbdrCondCpIdx,&
          p_DmaxPar, p_BisSegClosed, rboundaryCondition&
          %iboundarycount, p_IverticesAtBoundary, p_IboundaryCpIdx,&
          p_DvertexParameterValue, rvector%NVAR, p_Dx, dvalue)

    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarValue')
      call sys_halt()
    end select

    ! Do we have to re-sort vector?
    if (bisSorted) call lsyssc_vectorActivateSorting(rvector, .true.)

  contains

    ! Here are the real filter routines

    !***************************************************************
    ! Filter vector in 1D.

    subroutine filtervector_1D(IbdrCondType, IbdrCondCpIdx, nbct,&
        IverticesAtBoundary, IboundaryCpIdx, nvar, Dx, dvalue)

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Number of variables
      integer, intent(in) :: nvar

      ! Vector data array
      real(DP), dimension(nvar,*), intent(inout) :: Dx

      ! Filter value
      real(DP), intent(in) :: dvalue


      ! local variables
      integer :: ivbd,ivt,ibct,isegment


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          ! Do nothing

        case default
          ! Replace vector entry by DVALUE
          ivt = IverticesAtBoundary(ivbd)
          Dx(:,ivt) = dvalue

        end select
      end do
    end subroutine filtervector_1D


    !***************************************************************
    ! Filter vector in 2D.

    subroutine filtervector_2D(IbdrCondType, IbdrCondCpIdx,&
        DmaxParam, BisSegClosed, nbct, IverticesAtBoundary,&
        IboundaryCpIdx, DvertexParameterValue, nvar, Dx, dvalue)

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

       ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Number of variables
      integer, intent(in) :: nvar

      ! Vector data array
      real(DP), dimension(nvar,*), intent(inout) :: Dx

      ! Filter value
      real(DP), intent(in) :: dvalue


      ! local variables
      real(DP) :: dminValue,dmaxValue,dVertexParameterPeriodic
      integer :: ivbd,ivbdFirst,ivbdLast,ivt,ibct,isegment


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue =&
            nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle

          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))
            
          case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
            ! Do nothing

          case default
            ! Replace vector entry by DVALUE
            ivt = IverticesAtBoundary(ivbd)
            Dx(:,ivt) = dvalue

          end select
        end do
      end do
    end subroutine filtervector_2D
  end subroutine bdrf_filterVectorScalarValue

  ! ***************************************************************************

!<subroutine>

  subroutine bdrf_filterVectorBlockVector(rboundaryCondition,&
      rvector, rvectorFilter, rtriangulation)

!<description>
    ! This subroutine performs vector filtering by another vector.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! OPTIONAL: The triangulation
    type(t_triangulation), intent(in), optional :: rtriangulation
!</input>

!<inputoutput>
    ! The vector filter
    type(t_vectorBlock), intent(inout) :: rvectorFilter

    ! The vector that should be filtered
    type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iblock

    ! Check if there are strong boundary conditions
    if (.not.rboundaryCondition%bStrongBdrCond) return

    ! Check if both vectors have the same number of blocks
    if (rvector%nblocks .ne. rvectorFilter%nblocks) then
      call output_line('Vectors are not compatible!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorBlockVector')
      call sys_halt()
    end if

    ! Loop over all blocks
    do iblock = 1, rvector%nblocks
      call bdrf_filterVectorScalarVector(rboundaryCondition,&
          rvector%RvectorBlock(iblock),&
          rvectorFilter%RvectorBlock(iblock), rtriangulation)
    end do
  end subroutine bdrf_filterVectorBlockVector

  ! ***************************************************************************

!<subroutine>

  subroutine bdrf_filterVectorScalarVector(rboundaryCondition,&
      rvector, rvectorFilter, rtriangulation)

!<description>
    ! This subroutine performs vector filtering by a constant value.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! OPTIONAL: The triangulation
    type(t_triangulation), intent(in), optional, target :: rtriangulation
!</input>

!<inputoutput>
    ! The vector filter
    type(t_vectorScalar), intent(inout) :: rvectorFilter

    ! The vector that should be filtered
    type(t_vectorScalar), intent(inout) :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:), pointer :: p_Dx, p_Dvalue
    real(DP), dimension(:), pointer :: p_DmaxPar
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer, dimension(:), pointer :: p_IbdrCondType
    logical, dimension(:), pointer :: p_BisSegClosed
    logical :: bisSorted

    ! Check if there are strong boundary conditions
    if (.not.rboundaryCondition%bStrongBdrCond) return

    ! Check if both vectors are compatible
    call lsyssc_isVectorCompatible(rvector, rvectorFilter)

    ! Check if vector is sorted?
    if (rvector%isortStrategy .gt. 0) then
      bisSorted = .true.
      call lsyssc_vectorActivateSorting(rvector, .false.)
      call lsyssc_vectorActivateSorting(rvectorFilter, .false.)
    else
      bisSorted = .false.
    end if

    ! Set pointer for vector
    call lsyssc_getbase_double (rvector, p_Dx)
    call lsyssc_getbase_double (rvectorFilter, p_Dvalue)

    ! Get underlying triangulation structure
    if (present(rtriangulation)) then
      p_rtriangulation => rtriangulation
    else
      if (.not.associated(rvector%p_rspatialDiscr)) then
        call output_line('No discretisation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarVector')
        call sys_halt()
      end if

      if (.not.associated(rvector%p_rspatialDiscr%p_rtriangulation)) then
        call output_line('No triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarVector')
        call sys_halt()
      end if
      p_rtriangulation => rvector%p_rspatialDiscr%p_rtriangulation
    end if

    ! Check spatial dimensions
    if (p_rtriangulation%ndim .ne. rboundaryCondition%ndimension) then
      call output_line('Spatial dimension mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarVector')
      call sys_halt()
    end if

    ! How many spatial dimensions do we have?
    select case (rboundaryCondition%ndimension)
    case (NDIM1D)
      ! Set pointers for triangulation
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)

      ! Set prescribed boundary values in 1D
      call filtervector_1D(p_IbdrCondType, p_IbdrCondCpIdx,&
          rboundaryCondition%iboundarycount, p_IverticesAtBoundary,&
          p_IboundaryCpIdx, rvector%NVAR, p_Dx, p_Dvalue)

    case (NDIM2D)
      ! Set pointers for triangulation
      call storage_getbase_double(&
          p_rtriangulation%h_DvertexParameterValue, p_DvertexParameterValue)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_double(rboundaryCondition%h_DmaxPar, p_DmaxPar)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
      call storage_getbase_logical(rboundaryCondition%h_BisSegClosed, p_BisSegClosed)

      ! Set prescribed boundary values in 2D
      call filtervector_2D(p_IbdrCondType, p_IbdrCondCpIdx,&
          p_DmaxPar, p_BisSegClosed, rboundaryCondition&
          %iboundarycount, p_IverticesAtBoundary, p_IboundaryCpIdx,&
          p_DvertexParameterValue, rvector%NVAR, p_Dx, p_Dvalue)

    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarVector')
      call sys_halt()
    end select

    ! Do we have to re-sort vector?
    if (bisSorted) then
      call lsyssc_vectorActivateSorting(rvector, .true.)
      call lsyssc_vectorActivateSorting(rvectorFilter, .true.)
    end if

  contains

    ! Here are the real filter routines

    !***************************************************************
    ! Filter vector in 1D.

    subroutine filtervector_1D(IbdrCondType, IbdrCondCpIdx, nbct,&
        IverticesAtBoundary, IboundaryCpIdx, nvar, Dx, Dvalue)

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Number of variables
      integer, intent(in) :: nvar

      ! Vector data array
      real(DP), dimension(nvar,*), intent(inout) :: Dx

      ! Filter data array
      real(DP), dimension(nvar,*), intent(in) :: Dvalue

      ! local variables
      integer :: ivbd,ivt,ibct,isegment


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          ! Do nothing

        case default
          ! Replace vector entry by DVALUE
          ivt = IverticesAtBoundary(ivbd)
          Dx(:,ivt) = Dvalue(:,ivt)

        end select
      end do
    end subroutine filtervector_1D


    !***************************************************************
    ! Filter vector in 2D.

    subroutine filtervector_2D(IbdrCondType, IbdrCondCpIdx,&
        DmaxParam, BisSegClosed, nbct, IverticesAtBoundary,&
        IboundaryCpIdx, DvertexParameterValue, nvar, Dx, Dvalue)

        ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

       ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Number of variables
      integer, intent(in) :: nvar

      ! Vector data array
      real(DP), dimension(nvar,*), intent(inout) :: Dx

      ! Filter data array
      real(DP), dimension(nvar,*), intent(in) :: Dvalue


      ! local variables
      real(DP) :: dminValue,dmaxValue
      integer :: ivbd,ivbdFirst,ivbdLast,ivt,ibct,isegment


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue = nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle

          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

          case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
            ! Do nothing

          case default
            ! Replace vector entry by DVALUE
            ivt = IverticesAtBoundary(ivbd)
            Dx(:,ivt) = Dvalue(:,ivt)

          end select
        end do
      end do
    end subroutine filtervector_2D
  end subroutine bdrf_filterVectorScalarVector

  ! ***************************************************************************

!<subroutine>

  subroutine bdrf_filterVectorBlockExplicit(rboundaryCondition,&
      rvector, ttime, fcb_calcBoundaryvalues, istatus, rboundary,&
      rtriangulation)

!<description>
    ! This subroutine performs vector filtering by imposing the
    ! prescribed boundary conditions explicitly.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! The simulation time
    real(DP), intent(in) :: ttime

    ! OPTIONAL: The boundary description
    type(t_boundary), intent(in), optional, target :: rboundary

    ! OPTIONAL: The triangulation
    type(t_triangulation), intent(in), optional, target :: rtriangulation

    ! OPTIONAL: The callback function
    include 'intf_bdrcallback.inc'
    optional :: fcb_calcBoundaryvalues
!</input>

!<inputoutput>
    ! The block vector that should be filtered
    type(t_vectorBlock), intent(inout) :: rvector

    ! OPTIONAL: Status of the callback function
    integer, intent(inout), optional :: istatus
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundary), pointer :: p_rboundary
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Dx
    real(DP), dimension(:), pointer :: p_DmaxPar
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer, dimension(:), pointer :: p_IbdrCondType
    logical, dimension(:), pointer :: p_BisSegClosed
    logical, dimension(rvector%nblocks) :: BisSorted
    integer :: nvt,iblock


    ! Initialise status
    if (present(istatus)) istatus = 0

    ! Check if there are strong boundary conditions
    if (.not.rboundaryCondition%bStrongBdrCond) return

    ! Check if block vector has only one block
    if (rvector%nblocks .eq. 1) then
      call bdrf_filterVectorScalarExplicit(rboundaryCondition,&
          rvector%RvectorBlock(1), ttime, fcb_calcBoundaryvalues,&
          istatus, rboundary, rtriangulation)
      ! That is it, return
      return
    end if

    ! Check if vector and boundary are compatible
    if (rvector%nblocks .ne. rboundaryCondition%nmaxExpressions) then
      call output_line('Vector and boundary description are not compatible',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorBlockExplicit')
      call sys_halt()
    end if

    ! Get number of vertices
    nvt = int(rvector%NEQ/rvector%nblocks, I32)

    do iblock = 1, rvector%nblocks

      ! Check if vector is compatible
      if (rvector%RvectorBlock(iblock)%NEQ .ne. nvt) then
        call output_line('Subvector and boundary description are not compatible',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorBlockExplicit')
        call sys_halt()
      end if

      ! Check if vector is sorted
      if (rvector%RvectorBlock(iblock)%isortStrategy .gt. 0) then
        BisSorted(iblock) = .true.
        call lsyssc_vectorActivateSorting(rvector%RvectorBlock(iblock), .false.)
      else
        BisSorted(iblock) = .false.
      end if
    end do

    ! Set pointer for vector
    call lsysbl_getbase_double(rvector, p_Dx)

    ! Get underlying triangulation structure
    if (present(rtriangulation)) then
      p_rtriangulation => rtriangulation
    else
      if (.not.associated(rvector%RvectorBlock(1)%p_rspatialDiscr)) then
        call output_line('No discretisation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorBlockExplicit')
        call sys_halt()
      end if

      if (.not.associated(rvector%RvectorBlock(1)%p_rspatialDiscr&
          %p_rtriangulation)) then
        call output_line('No triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorBlockExplicit')
        call sys_halt()
      end if
      p_rtriangulation => rvector%RvectorBlock(1)%p_rspatialDiscr&
          %p_rtriangulation
    end if

    ! Check spatial dimensions
    if (p_rtriangulation%ndim .ne. rboundaryCondition%ndimension) then
      call output_line('Spatial dimension mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorBlockExplicit')
      call sys_halt()
    end if

    ! How many spatial dimensions do we have?
    select case (rboundaryCondition%ndimension)
    case (NDIM1D)
      ! Set pointers for triangulation
      call storage_getbase_double2d(&
          p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)

      ! Impose boundary conditions explicitly in 1D
      call filtervector_1D(rboundaryCondition%rfparser,&
          p_IbdrCondType, p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
          rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
          p_IboundaryCpIdx, p_DvertexCoords, nvt, rvector%nblocks, p_Dx,&
          fcb_calcBoundaryvalues, istatus)

    case (NDIM2D)
      ! Set pointers for triangulation
      call storage_getbase_double2d(&
          p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
      call storage_getbase_double(&
          p_rtriangulation%h_DvertexParameterValue, p_DvertexParameterValue)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_double(rboundaryCondition%h_DmaxPar, p_DmaxPar)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
      call storage_getbase_logical(rboundaryCondition%h_BisSegClosed, p_BisSegClosed)

      ! Get underlying boundary structure
      if (present(rboundary)) then
        p_rboundary => rboundary
      else
        if (.not.associated(rvector%RvectorBlock(1)%p_rspatialDiscr&
            %p_rboundary)) then
          call output_line('No boundary associated!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorBlockExplicit')
          call sys_halt()
        end if
        p_rboundary => rvector%RvectorBlock(1)%p_rspatialDiscr&
            %p_rboundary
      end if

      ! Impose boundary conditions explicitly in 2D
      call filtervector_2D(rboundaryCondition%rfparser,&
          p_IbdrCondType, p_IbdrCondCpIdx, p_DmaxPar, p_BisSegClosed,&
          rboundaryCondition%iboundarycount, rboundaryCondition%nmaxExpressions,&
          p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
          p_DvertexCoords, nvt, rvector%nblocks, p_Dx, p_rboundary,&
          fcb_calcBoundaryvalues, istatus)

    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorBlockExplicit')
      call sys_halt()
    end select

    ! Do we have to re-sort the vector?
    do iblock = 1, rvector%nblocks
      if (BisSorted(iblock))&
          call lsyssc_vectorActivateSorting(rvector%RvectorBlock(iblock), .true.)
    end do

  contains

    ! Here are the real filter routines

    !***************************************************************
    ! Filter vector in 1D.

    subroutine filtervector_1D(rfparser, IbdrCondType, IbdrCondCpIdx,&
        nbct, nexpr, IverticesAtBoundary, IboundaryCpIdx, DvertexCoords,&
        nvt, nvar, Dx, fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Total number of vertices
      integer, intent(in) :: nvt

      ! Number of variables
      integer, intent(in) :: nvar

      ! Vector data array
      real(DP), dimension(nvt,nvar), intent(inout) :: Dx

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(nvar) :: Daux,Daux0
      real(DP), dimension(NDIM1D) :: DbdrNormal
      integer :: ivbd,ivt,ibct,ivar,isegment,iexpr


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          ! Do nothing
          
        case default
          ! Get vertex number
          ivt = IverticesAtBoundary(ivbd)

          ! Initialise variable values [x,0,0,time] for parser
          DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
          DvariableValues(NDIM2D)   = 0.0_DP
          DvariableValues(NDIM3D)   = 0.0_DP
          DvariableValues(NDIM3D+1) = ttime

          ! Check for user-defined callback-function
          if (present(fcb_calcBoundaryvalues)) then
            
            ! Get desired boundary values from parser
            do iexpr = 1, nexpr
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                  DvariableValues, DbdrValues(iexpr))
            end do
            
            ! Compute the analytical normal vector
            if (mod(ibct, 2) .eq. 0) then
              DbdrNormal = 1.0_DP
            else
              DbdrNormal = -1.0_DP
            end if
            
            ! Apply callback function to determine boundary conditions
            Daux  = Dx(ivt,:)
            Daux0 = Dx(ivt,:)
            call fcb_calcBoundaryvalues(DbdrNormal, DbdrNormal,&
                DbdrValues, IbdrCondType(isegment), Daux, Daux0, istatus)
            Dx(ivt,:) = Daux

          else
            
            ! Impose prescribed boundary conditions directly
            do ivar = 1, nvar
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                  DvariableValues, Dx(ivt,ivar))
            end do
            
          end if

        end select
      end do
    end subroutine filtervector_1D


    !***************************************************************
    ! Filter vector in 2D.

    subroutine filtervector_2D(rfparser, IbdrCondType, IbdrCondCpIdx,&
        DmaxParam, BisSegClosed, nbct, nexpr, IverticesAtBoundary,&
        IboundaryCpIdx, DvertexParameterValue, DvertexCoords, nvt,&
        nvar, Dx, rboundary, fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Total number of vertices
      integer, intent(in) :: nvt

      ! Number of variables
      integer, intent(in) :: nvar

      ! Vector data array
      real(DP), dimension(nvt,nvar), intent(inout) :: Dx

      ! Boundary description
      type(t_boundary), intent(in) :: rboundary

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(nvar) :: Daux,Daux0
      real(DP), dimension(NDIM2D) :: DbdrNormal,DpointNormal
      real(DP) :: dminValue,dmaxValue,dnx,dny,dnxL,dnxR,dnyL,dnyR,dw,dwL,dwR
      integer :: ivbd,ivbdFirst,ivbdLast,ivt,ivtL,ivtR,ibct,ivar,isegment,iexpr


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue = nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle

          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

          case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
            ! Do nothing

          case default
            ! Get vertex number
            ivt = IverticesAtBoundary(ivbd)

            ! Initialise variable values [x,y,0,time] for parser
            DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
            DvariableValues(NDIM2D)   = DvertexCoords(2,ivt)
            DvariableValues(NDIM3D)   = 0.0_DP
            DvariableValues(NDIM3D+1) = ttime

            ! Check for user-defined callback-function
            if (present(fcb_calcBoundaryvalues)) then
              
              ! Get desired boundary values from parser
              do iexpr = 1, nexpr
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                    DvariableValues, DbdrValues(iexpr))
              end do
              
              ! Get vertex number of predecessor
              if (ivbd .eq. ivbdFirst) then
                ivtL = IverticesAtBoundary(ivbdLast)
              else
                ivtL = IverticesAtBoundary(ivbd-1)
              end if
              
              ! Get vertex number of sucessor
              if (ivbd .eq. ivbdLast) then
                ivtR = IverticesAtBoundary(ivbdFirst)
              else
                ivtR = IverticesAtBoundary(ivbd+1)
              end if
              
              ! Compute the analytical and approximate normal vectors
              if (DvertexParameterValue(ivbd) .eq. dminValue) then
                ! We are at the beginning of a segment that includes
                ! the starting point. Otherwise, the segment counter
                ! would not have been increased.
                
                ! Thus, compute the analytical normal vector from the right
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_RIGHT)
                
                ! Compute the approximate normal vector from the right element
                dnx = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dny = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              elseif (DvertexParameterValue(ivbd) .eq. dmaxValue) then
                ! We are at the end of a segment that includes the end
                ! point. Otherwise, the segment counter would have
                ! been increased already.
                
                ! Thus, compute the analytical normal vector from the left.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_LEFT)
                
                ! Compute the approximate normal vector from the left element
                dnx = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dny = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              else
                ! We are "inside" a segment.
                ! Thus, compute the mean value of left and right normal vector.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2))
                
                ! Compute the approximate normal vector
                dnxL = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dnyL = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                
                dnxR = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dnyR = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                
                dwL = sqrt( dnxL*dnxL + dnyL*dnyL )
                dwR = sqrt( dnxR*dnxR + dnyR*dnyR )
                
                dnx = (dnxL/dwL + dnxR/dwR)/(1/dwL + 1/dwR)
                dny = (dnyL/dwL + dnyR/dwR)/(1/dwL + 1/dwR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
              end if
              
              ! Apply callback function to determine boundary conditions
              Daux  = Dx(ivt,:)
              Daux0 = Dx(ivt,:)
              call fcb_calcBoundaryvalues(DbdrNormal, DpointNormal,&
                  DbdrValues, IbdrCondType(isegment), Daux, Daux0, istatus)
              Dx(ivt,:) = Daux
              
            else
              
              ! Impose prescribed boundary conditions directly
              do ivar = 1, nvar
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                    Dvariablevalues, Dx(ivt,ivar))
              end do

            end if
            
          end select
        end do
      end do
    end subroutine filtervector_2D

  end subroutine bdrf_filterVectorBlockExplicit

  ! ***************************************************************************

!<subroutine>

  subroutine bdrf_filterVectorScalarExplicit(rboundaryCondition,&
      rvector, ttime, fcb_calcBoundaryvalues, istatus, rboundary,&
      rtriangulation)

!<description>
    ! This subroutine performs vector filtering by imposing the
    ! prescribed boundary conditions explicitly.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! The simulation time
    real(DP), intent(in) :: ttime

    ! OPTIONAL: The boundary description
    type(t_boundary), intent(in), optional, target :: rboundary

    ! OPTIONAL: The triangulation
    type(t_triangulation), intent(in), optional, target :: rtriangulation

    ! OPTIONAL: The callback function
    include 'intf_bdrcallback.inc'
    optional :: fcb_calcBoundaryvalues
!</input>

!<inputoutput>
    ! The scalar vector that should be filtered
    type(t_vectorScalar), intent(inout) :: rvector

    ! OPTIONAL: Status of the callback function
    integer, intent(inout), optional :: istatus
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundary), pointer :: p_rboundary
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Dx
    real(DP), dimension(:), pointer :: p_DmaxPar
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer, dimension(:), pointer :: p_IbdrCondType
    logical, dimension(:), pointer :: p_BisSegClosed
    logical :: bisSorted


    ! Initialise status
    if (present(istatus)) istatus = 0

    ! Check if there are strong boundary conditions
    if (.not.rboundaryCondition%bStrongBdrCond) return

    ! Check if vector is sorted?
    if (rvector%isortStrategy .gt. 0) then
      bisSorted = .true.
      call lsyssc_vectorActivateSorting(rvector, .false.)
    else
      bisSorted = .false.
    end if

    ! Set pointer for vector
    call lsyssc_getbase_double(rvector, p_Dx)

    ! Get underlying triangulation structure
    if (present(rtriangulation)) then
      p_rtriangulation => rtriangulation
    else
      if (.not.associated(rvector%p_rspatialDiscr)) then
        call output_line('No discretisation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarExplicit')
        call sys_halt()
      end if

      if (.not.associated(rvector%p_rspatialDiscr%p_rtriangulation)) then
        call output_line('No triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarExplicit')
        call sys_halt()
      end if
      p_rtriangulation => rvector%p_rspatialDiscr%p_rtriangulation
    end if

    ! Check spatial dimensions
    if (p_rtriangulation%ndim .ne. rboundaryCondition%ndimension) then
      call output_line('Spatial dimension mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarExplicit')
      call sys_halt()
    end if

    ! How many spatial dimensions do we have?
    select case (rboundaryCondition%ndimension)
    case (NDIM1D)
      ! Set pointers for triangulation
      call storage_getbase_double2d(&
          p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)

      ! Impose boundary conditions explicitly in 1D
      call filtervector_1D(rboundaryCondition%rfparser,&
          p_IbdrCondType, p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
          rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
          p_IboundaryCpIdx, p_DvertexCoords, rvector%NVAR, p_Dx,&
          fcb_calcBoundaryvalues, istatus)

    case (NDIM2D)
      ! Set pointers for triangulation
      call storage_getbase_double2d(&
          p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
      call storage_getbase_double(&
          p_rtriangulation%h_DvertexParameterValue, p_DvertexParameterValue)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_double(rboundaryCondition%h_DmaxPar, p_DmaxPar)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
      call storage_getbase_logical(rboundaryCondition%h_BisSegClosed, p_BisSegClosed)

      ! Get underlying boundary structure
      if (present(rboundary)) then
        p_rboundary => rboundary
      else
        if (.not.associated(rvector%p_rspatialDiscr%p_rboundary)) then
          call output_line('No boundary associated!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarExplicit')
          call sys_halt()
        end if
        p_rboundary => rvector%p_rspatialDiscr%p_rboundary
      end if

      ! Impose boundary conditions explicitly in 2D
      call filtervector_2D(rboundaryCondition%rfparser,&
          p_IbdrCondType, p_IbdrCondCpIdx, p_DmaxPar, p_BisSegClosed,&
          rboundaryCondition%iboundarycount, rboundaryCondition%nmaxExpressions,&
          p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
          p_DvertexCoords, rvector%NVAR, p_Dx, p_rboundary,&
          fcb_calcBoundaryvalues, istatus)

    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterVectorScalarExplicit')
      call sys_halt()
    end select

    ! Do we have to re-sort vector
    if (bisSorted) call lsyssc_vectorActivateSorting(rvector, .true.)

  contains

    ! Here are the real filter routines

    !***************************************************************
    ! Filter vector in 1D.

    subroutine filtervector_1D(rfparser, IbdrCondType, IbdrCondCpIdx,&
        nbct, nexpr, IverticesAtBoundary, IboundaryCpIdx, DvertexCoords,&
        nvar, Dx, fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Number of variables
      integer, intent(in) :: nvar

      ! Vector data array
      real(DP), dimension(nvar,*), intent(inout) :: Dx

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(NDIM1D) :: DbdrNormal
      integer :: ivbd,ivt,ivar,ibct,isegment,iexpr


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle
        
        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          ! Do nothing

        case default
          ! Get vertex number
          ivt = IverticesAtBoundary(ivbd)

          ! Initialise variable values [x,0,0,time] for parser
          DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
          DvariableValues(NDIM2D)   = 0.0_DP
          DvariableValues(NDIM3D)   = 0.0_DP
          DvariableValues(NDIM3D+1) = ttime

          ! Check for user-defined callback-function
          if (present(fcb_calcBoundaryvalues)) then
            
            ! Get desired boundary values from parser
            do iexpr = 1, nexpr
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                  DvariableValues, DbdrValues(iexpr))
            end do
            
            ! Compute the analytical normal vector
            if (mod(ibct, 2) .eq. 0) then
              DbdrNormal = 1.0_DP
            else
              DbdrNormal = -1.0_DP
            end if
            
            ! Apply callback function to determine boundary conditions
            call fcb_calcBoundaryvalues(DbdrNormal, DbdrNormal, DbdrValues,&
                IbdrCondType(isegment), Dx(:,ivt), Dx(:,ivt), istatus)
            
          else
            
            ! Impose prescribed boundary conditions directly
            do ivar = 1, nvar
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                  DvariableValues, Dx(ivar,ivt))
            end do

          end if
          
        end select
      end do
    end subroutine filtervector_1D

    !***************************************************************
    ! Filter vector in 2D.

    subroutine filtervector_2D(rfparser, IbdrCondType, IbdrCondCpIdx,&
        DmaxParam, BisSegClosed, nbct, nexpr, IverticesAtBoundary,&
        IboundaryCpIdx, DvertexParameterValue, DvertexCoords, nvar,&
        Dx, rboundary, fcb_calcBoundaryvalues, istatus)

       ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Number of variables
      integer, intent(in) :: nvar

      ! Vector data array
      real(DP), dimension(nvar,*), intent(inout) :: Dx

      ! Boundary description
      type(t_boundary), intent(in) :: rboundary

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(NDIM2D) :: DbdrNormal,DpointNormal
      real(DP) :: dminValue,dmaxValue,dnx,dny,dnxL,dnxR,dnyL,dnyR,dw,dwL,dwR
      integer :: ivbd,ivbdFirst,ivbdLast,ivt,ivtL,ivtR,ivar,ibct,isegment,iexpr


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue = nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle

          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

          case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
            ! Do nothing

          case default
            ! Get vertex number
            ivt = IverticesAtBoundary(ivbd)

            ! Initialise variable values [x,y,0,time] for parser
            DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
            DvariableValues(NDIM2D)   = DvertexCoords(2,ivt)
            DvariableValues(NDIM3D)   = 0.0_DP
            DvariableValues(NDIM3D+1) = ttime

            ! Check for user-defined callback-function
            if (present(fcb_calcBoundaryvalues)) then

              ! Get desired boundary values from parser
              do iexpr = 1, nexpr
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                    DvariableValues, DbdrValues(iexpr))
              end do

              ! Get vertex number of predecessor
              if (ivbd .eq. ivbdFirst) then
                ivtL = IverticesAtBoundary(ivbdLast)
              else
                ivtL = IverticesAtBoundary(ivbd-1)
              end if
              
              ! Get vertex number of sucessor
              if (ivbd .eq. ivbdLast) then
                ivtR = IverticesAtBoundary(ivbdFirst)
              else
                ivtR = IverticesAtBoundary(ivbd+1)
              end if
              
              ! Compute the analytical and approximate normal vectors
              if (DvertexParameterValue(ivbd) .eq. dminValue) then
                ! We are at the beginning of a segment that includes
                ! the starting point.  Otherwise, the segment counter
                ! would not have been increased.
                
                ! Thus, compute the analytical normal vector from the right
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_RIGHT)
                
                ! Compute the approximate normal vector from the right element
                dnx = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dny = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              elseif (DvertexParameterValue(ivbd) .eq. dmaxValue) then
                ! We are at the end of a segment that includes the end
                ! point.  Otherwise, the segment counter would have
                ! been increased already.
                
                ! Thus, compute the analytical normal vector from the left.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_LEFT)
                
                ! Compute the approximate normal vector from the left element
                dnx = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dny = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              else
                ! We are "inside" a segment.
                ! Thus, compute the mean value of left and right normal vector.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2))
                
                ! Compute the approximate normal vector
                dnxL = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dnyL = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                
                dnxR = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dnyR = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                
                dwL = sqrt( dnxL*dnxL + dnyL*dnyL )
                dwR = sqrt( dnxR*dnxR + dnyR*dnyR )
                
                dnx = (dnxL/dwL + dnxR/dwR)/(1/dwL + 1/dwR)
                dny = (dnyL/dwL + dnyR/dwR)/(1/dwL + 1/dwR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
              end if
              
              ! Apply callback function to determine boundary conditions
              call fcb_calcBoundaryvalues(DbdrNormal, DpointNormal, DbdrValues,&
                  IbdrCondType(isegment), Dx(:,ivt), Dx(:,ivt), istatus)

            else

              ! Impose prescribed boundary conditions directly
              do ivar = 1, nvar
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                    DvariableValues, Dx(ivar,ivt))
              end do
              
            end if

          end select
        end do
      end do
    end subroutine filtervector_2D

  end subroutine bdrf_filterVectorScalarExplicit

  ! ***************************************************************************

!<subroutine>

  subroutine bdrf_filterSolutionBlock(rboundaryCondition, rmatrix,&
      rsolution, rdefect, rsolution0, ttime, fcb_calcBoundaryvalues,&
      istatus, rboundary, rtriangulation)

!<description>
    ! This subroutine imposes prescribed boundary conditions by
    ! filtering the solution vector in the following way:
    !
    ! 1. Compute a provisional solution by scaling the defect vector
    !    by the inverse of the diagonal matrix entries and adding
    !    the result to the solution from the last iteration:
    !
    !    $$U^*_i=U_i+A_{ii}^{-1}R_i,\qquad \mbox{whereby}\quad
    !    A_{ii}=\mbox{diag}\{a_{ii}^{kk}$$
    !
    ! 2. Modify the predicted solution values by means of
    !    an external subroutine or impose the prescribed Dirichlet
    !    boundary values
    !
    !    $$U_i^{**}=\mbox{bc}(U_i^*)$$
    !
    ! 3. Nullify all entries of the defect vector $R_i=0$
    !
    ! Then the increment $\Delta U_i=A_{ii}^{-1}R_i=0$ so that
    ! for the updated solution vector $U_i=U_i^{**}$.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! The initial solution vector from the last time step
    type(t_vectorBlock), intent(in) :: rsolution0

    ! The simulation time
    real(DP), intent(in) :: ttime

    ! OPTIONAL: The boundary description
    type(t_boundary), intent(in), optional, target :: rboundary

    ! OPTIONAL: The triangulation
    type(t_triangulation), intent(in), optional, target :: rtriangulation

    ! OPTIONAL: The callback function
    include 'intf_bdrcallback.inc'
    optional :: fcb_calcBoundaryvalues
!</input>

!<inputoutput>
    ! The matrix data structure
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! The solution vector that should be filtered
    type(t_vectorBlock), intent(inout) :: rsolution

    ! The defect vector that should be filtered
    type(t_vectorBlock), intent(inout) :: rdefect

    ! OPTIONAL: Status of the callback function
    integer, intent(inout), optional :: istatus
!</inputoutput>
!</subroutine>

    ! local type
    type t_rarray
      real(DP), dimension(:), pointer :: Da
    end type t_rarray


    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundary), pointer :: p_rboundary
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DmaxPar
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    real(DP), dimension(:), pointer :: p_Du, p_Dr, p_Du0
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    integer, dimension(:), pointer :: p_Kcol
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer, dimension(:), pointer :: p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCompPeriodic
    integer, dimension(:), pointer :: p_IbdrCondPeriodic
    logical, dimension(:), pointer :: p_BisSegClosed

    logical, dimension(rmatrix%nblocksPerCol,&
                       rmatrix%nblocksPerRow) :: BisMatrixSorted
    logical, dimension(rsolution%nblocks) :: BisSolutionSorted
    logical, dimension(rdefect%nblocks) :: BisDefectSorted

    type(t_rarray), dimension(rmatrix%nblocksPerCol,&
                              rmatrix%nblocksPerRow) :: rarray

    integer :: iblock,jblock


    ! Initialise status
    if (present(istatus)) istatus = 0

    ! Check if there are strong boundary conditions
    if (.not.rboundaryCondition%bStrongBdrCond) return

    ! Check if block matrix/vector have only one block
    if ((rmatrix%nblocksPerCol .eq. 1) .and.&
        (rmatrix%nblocksPerRow .eq. 1)) then
      if (rsolution%nblocks .eq. 1 .and.&
          rdefect%nblocks   .eq. 1) then
        call bdrf_filterSolutionScalar(rboundaryCondition,&
            rmatrix%Rmatrixblock(1,1), rsolution%RvectorBlock(1),&
            rdefect%RvectorBlock(1), rsolution0%RvectorBlock(1),&
            ttime, fcb_calcBoundaryvalues, istatus, rboundary,&
            rtriangulation)
      else
        call bdrf_filterSolutionBlockScalar(rboundaryCondition,&
            rmatrix%Rmatrixblock(1,1), rsolution, rdefect,&
            rsolution0, ttime, fcb_calcBoundaryvalues,&
            istatus, rboundary, rtriangulation)
      end if

      ! That is it, return
      return
    end if

    ! Check if matrix exhibits group structure
    if (rmatrix%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionBlock')
      call sys_halt()
    end if

    ! Check if vector and boundary are compatible
    if (rsolution%nblocks .ne. rboundaryCondition%nmaxExpressions .or.&
        rdefect%nblocks   .ne. rboundaryCondition%nmaxExpressions) then
      call output_line('Vector and boundary description are not compatible',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionBlock')
      call sys_halt()
    end if

    ! Check if matrix is sorted and prepare auxiliary data structure?
    do iblock = 1, rmatrix%nblocksPerCol
      do jblock = 1, rmatrix%nblocksPerRow

        if (rmatrix%RmatrixBlock(jblock,iblock)%isortStrategy .gt. 0) then
          BisMatrixSorted(jblock,iblock) =.true.
          call lsyssc_unsortMatrix(rmatrix%RmatrixBlock(jblock,iblock), .true.)
        else
          BisMatrixSorted(jblock,iblock) =.false.
        end if

        if (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(jblock,iblock))) then
          call lsyssc_getbase_double(rmatrix%RmatrixBlock(jblock,iblock),&
              rarray(jblock,iblock)%Da)
        else
          nullify(rarray(jblock,iblock)%Da)
        end if
      end do
    end do

    ! Check if solution vector is sorted?
    do iblock = 1, rsolution%nblocks
      if (rsolution%RvectorBlock(iblock)%isortStrategy .gt. 0) then
        BisSolutionSorted(iblock) = .true.
        call lsyssc_vectorActivateSorting(rsolution%RvectorBlock(iblock), .false.)
      else
        BisSolutionSorted(iblock) = .false.
      end if
    end do

    ! Check if defect vector is sorted?
    do iblock = 1, rdefect%nblocks
      if (rdefect%RvectorBlock(iblock)%isortStrategy .gt. 0) then
        BisDefectSorted(iblock) = .true.
        call lsyssc_vectorActivateSorting(rdefect%RvectorBlock(iblock), .false.)
      else
        BisDefectSorted(iblock) = .false.
      end if
    end do

    ! Set pointers for solution and defect
    call lsysbl_getbase_double(rsolution0, p_Du0)
    call lsysbl_getbase_double(rsolution, p_Du)
    call lsysbl_getbase_double(rdefect, p_Dr)

    ! Get underlying triangulation structure
    if (present(rtriangulation)) then
      p_rtriangulation => rtriangulation
    else
      if (.not.associated(rsolution%RvectorBlock(1)%p_rspatialDiscr)) then
        call output_line('No discretisation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionBlock')
        call sys_halt()
      end if

      if (.not.associated(rsolution%RvectorBlock(1)%p_rspatialDiscr&
          %p_rtriangulation)) then
        call output_line('No triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionBlock')
        call sys_halt()
      end if
      p_rtriangulation => rsolution%RvectorBlock(1)%p_rspatialDiscr&
          %p_rtriangulation
    end if

    ! Check spatial dimensions
    if (p_rtriangulation%ndim .ne. rboundaryCondition%ndimension) then
      call output_line('Spatial dimension mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionBlock')
      call sys_halt()
    end if

    ! How many spatial dimensions do we have?
    select case (rboundaryCondition%ndimension)
    case (NDIM1D)
      ! Set pointers for triangulation
      call storage_getbase_double2d(&
          p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
      call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic, p_IbdrCompPeriodic)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic, p_IbdrCondPeriodic)

      ! What kind of matrix are we?
      select case (rmatrix%RmatrixBlock(1,1)%cmatrixFormat)

      case (LSYSSC_MATRIX7)
        call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,1), p_Kld)
        call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1), p_Kcol)

        ! Set prescribed boundary values in 1D
        call filtersolution_Mat79_1D(rboundaryCondition%rfparser,&
            p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
            p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
            rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
            p_IboundaryCpIdx, p_DvertexCoords, p_Kld, p_Kcol, p_Kld,&
            rmatrix%RmatrixBlock(1,1)%NEQ, rmatrix%nblocksPerCol,&
            rarray, p_Du, p_Dr, p_Du0, fcb_calcBoundaryvalues, istatus)

      case (LSYSSC_MATRIX9)
        call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), p_Kdiagonal)
        call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), p_Kld)
        call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1), p_Kcol)

        ! Set prescribed boundary values in 1D
        call filtersolution_Mat79_1D(rboundaryCondition%rfparser,&
            p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
            p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
            rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
            p_IboundaryCpIdx, p_DvertexCoords, p_Kld, p_Kcol,&
            p_Kdiagonal, rmatrix%RmatrixBlock(1,1)%NEQ,&
            rmatrix%nblocksPerCol, rarray, p_Du, p_Dr, p_Du0,&
            fcb_calcBoundaryvalues, istatus)

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionBlock')
        call sys_halt()
      end select


    case (NDIM2D)
      ! Set pointers for triangulation
      call storage_getbase_double2d(&
          p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
      call storage_getbase_double(&
          p_rtriangulation%h_DvertexParameterValue, p_DvertexParameterValue)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_double(rboundaryCondition%h_DmaxPar, p_DmaxPar)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
      call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic, p_IbdrCompPeriodic)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic, p_IbdrCondPeriodic)
      call storage_getbase_logical(rboundaryCondition%h_BisSegClosed, p_BisSegClosed)

      ! Get underlying boundary structure
      if (present(rboundary)) then
        p_rboundary => rboundary
      else
        if (.not.associated(rsolution%RvectorBlock(1)%p_rspatialDiscr%p_rboundary)) then
          call output_line('No boundary associated!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionBlock')
          call sys_halt()
        end if
        p_rboundary => rsolution%RvectorBlock(1)%p_rspatialDiscr%p_rboundary
      end if

      ! What kind of matrix are we?
      select case (rmatrix%RmatrixBlock(1,1)%cmatrixFormat)

      case (LSYSSC_MATRIX7)
        call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), p_Kld)
        call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1), p_Kcol)

        ! Set prescribed boundary values in 2D
        call filtersolution_Mat79_2D(rboundaryCondition%rfparser,&
            p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
            p_IbdrCondCpIdx, p_DmaxPar, p_BisSegClosed,&
            rboundaryCondition%iboundarycount, rboundaryCondition%nmaxExpressions,&
            p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
            p_DvertexCoords, p_Kld, p_Kcol, p_Kld, rmatrix%RmatrixBlock(1,1)%NEQ,&
            rmatrix%nblocksPerCol, rarray, p_Du, p_Dr, p_Du0, p_rboundary,&
            fcb_calcBoundaryvalues, istatus)

      case (LSYSSC_MATRIX9)
        call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), p_Kdiagonal)
        call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), p_Kld)
        call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1), p_Kcol)

        ! Set prescribed boundary values in 2D
        call filtersolution_Mat79_2D(rboundaryCondition%rfparser,&
            p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
            p_IbdrCondCpIdx, p_DmaxPar, p_BisSegClosed,&
            rboundaryCondition%iboundarycount, rboundaryCondition%nmaxExpressions,&
            p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
            p_DvertexCoords, p_Kld, p_Kcol, p_Kdiagonal, rmatrix%RmatrixBlock(1,1)%NEQ,&
            rmatrix%nblocksPerCol, rarray, p_Du, p_Dr, p_Du0, p_rboundary,&
            fcb_calcBoundaryvalues, istatus)

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionBlock')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionBlock')
      call sys_halt()
    end select

    ! Do we have to re-sort matrix?
    do iblock = 1, rmatrix%nblocksPerRow
      do jblock = 1, rmatrix%nblocksPerCol
        if (BisMatrixSorted(jblock,iblock))&
            call lsyssc_sortMatrix(rmatrix%RmatrixBlock(jblock,iblock),&
            .true., rmatrix%RmatrixBlock(jblock,iblock)%isortStrategy)
      end do
    end do

    ! Do we have to re-sort solution vector?
    do iblock = 1, rsolution%nblocks
      if (BisSolutionSorted(iblock))&
          call lsyssc_vectorActivateSorting(rsolution%RvectorBlock(iblock), .true.)
    end do

    ! Do we have to re-sort defect vector
    do iblock = 1, rdefect%nblocks
      if (BisDefectSorted(iblock))&
          call lsyssc_vectorActivateSorting(rdefect%RvectorBlock(iblock), .true.)
    end do

  contains

    ! Here are the real filter routines

    !***************************************************************
    ! Filter solution and defect vector in 1D
    ! Here, the matrix is store in CSR7 and CSR9 format.

    subroutine filtersolution_Mat79_1D(rfparser, IbdrCompPeriodic,&
        IbdrCondPeriodic, IbdrCondType, IbdrCondCpIdx, nbct, nexpr,&
        IverticesAtBoundary, IboundaryCpIdx, DvertexCoords, Kld,&
        Kcol, Kdiagonal, neq, nvar, rarray, Du, Dr, Du0,&
        fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Number of equations per variable
      integer, intent(in) :: neq

      ! Number of variables
      integer, intent(in) :: nvar

      ! Auxiliary data structure representing the block matrix
      type(t_rarray), dimension(:,:), intent(in) :: rarray

      ! Solution and residual data array
      real(DP), dimension(neq,nvar), intent(inout) :: Du, Dr

      ! Initial solution data array
      real(DP), dimension(neq,nvar), intent(in) :: Du0

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus

      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(nvar) :: Daux,Daux0
      real(DP), dimension(NDIM1D) :: DbdrNormal
      integer :: ivbd,ivbdPeriodic,ivt,ivtPeriodic,ipos
      integer :: ibct,ibctPeriodic,isegment,ivar,jvar,iexpr

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          ! Compute vertex parameter value at periodic boundary
          ibctPeriodic = IbdrCompPeriodic(isegment)
          ivbdPeriodic = IboundaryCpIdx(ibctPeriodic)

          ! Get numbers of vertices
          ivt         = IverticesAtBoundary(ivbd)
          ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

          ! Add residual of the equation that corresponds to node ivt
          ! to the equation that corresponds to node ivtPeriodic and
          ! set residual of first equation equal to zero
          do ivar = 1, NVAR
            Dr(ivtPeriodic,ivar) = Dr(ivt,ivar)+Dr(ivtPeriodic,ivar)
            Dr(ivt,ivar)         = 0.0_DP
          end do
          
        case default
          ! Get vertex number
          ivt  = IverticesAtBoundary(ivbd)
          ipos = Kdiagonal(ivt)

          ! Initialise variable values [x,0,0,time] for parser
          DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
          DvariableValues(NDIM2D)   = 0.0_DP
          DvariableValues(NDIM3D)   = 0.0_DP
          DvariableValues(NDIM3D+1) = ttime

          ! Check for user-defined callback-function
          if (present(fcb_calcBoundaryvalues)) then
            
            ! Predict boundary values
            do ivar = 1, nvar
              Daux(ivar) = Du(ivt,ivar)+Dr(ivt,ivar)/rarray(ivar,ivar)%Da(ipos)
            end do
            
            ! Get values from old solution
            do ivar = 1, nvar
              Daux0(ivar) = Du0(ivt,ivar)
            end do

            ! Get desired boundary values from parser
            do iexpr = 1, nexpr
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                  DvariableValues, DbdrValues(iexpr))
            end do
            
            ! Compute the analytical normal vector
            if (mod(ibct, 2) .eq. 0) then
              DbdrNormal = 1.0_DP
            else
              DbdrNormal = -1.0_DP
            end if
            
            ! Apply callback function to determine boundary conditions
            call fcb_calcBoundaryvalues(DbdrNormal, DbdrNormal, DbdrValues,&
                IbdrCondType(isegment), Daux, Daux0, istatus)
            do ivar = 1, nvar
              Du(ivt,ivar) = Daux(ivar)
            end do
            
          else
            
            ! Impose prescribed boundary conditions directly
            do ivar = 1, nvar
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                  DvariableValues, Du(ivt,ivar))
            end do

          end if
          
          ! Nullify residual entries
          do ivar = 1, nvar
            Dr(ivt,ivar) = 0.0_DP
          end do

        end select
      end do
    end subroutine filtersolution_Mat79_1D


    !***************************************************************
    ! Filter solution and defect vector in 2D
    ! Here, the matrix is store in CSR7 and CSR9 format.

    subroutine filtersolution_Mat79_2D(rfparser, IbdrCompPeriodic,&
        IbdrCondPeriodic, IbdrCondType, IbdrCondCpIdx, DmaxParam,&
        BisSegClosed, nbct, nexpr, IverticesAtBoundary, IboundaryCpIdx,&
        DvertexParameterValue, DvertexCoords, Kld, Kcol, Kdiagonal,&
        neq, nvar, rarray, Du, Dr, Du0, rboundary,&
        fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Number of equations per variable
      integer, intent(in) :: neq

      ! Number of variables
      integer, intent(in) :: nvar

      ! Auxiliary data structure representing the block matrix
      type(t_rarray), dimension(:,:), intent(in) :: rarray

      ! Solution and residual data array
      real(DP), dimension(neq,nvar), intent(inout) :: Du, Dr

      ! Initial solution data array
      real(DP), dimension(neq,nvar), intent(in) :: Du0

      ! Boundary description
      type(t_boundary), intent(in) :: rboundary

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(nvar) :: Daux,Daux0
      real(DP), dimension(NDIM2D) :: DbdrNormal,DpointNormal
      real(DP) :: dminValue,dmaxValue,dVertexParameterPeriodic
      real(DP) :: dnx,dny,dnxL,dnxR,dnyL,dnyR,dw,dwL,dwR
      integer :: ibeg,iend,idiag,jdiag,ipos,jpos
      integer :: ivbd,ivbdFirst,ivbdLast,ivbdPeriodic,ivt,ivtPeriodic,ivtL,ivtR
      integer :: ibct,ibctPeriodic,isegment,isegmentPeriodic,ivar,jvar,iexpr

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue = nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle

          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

          case (BDRC_PERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            if (isegmentPeriodic .eq. IbdrCondCpIdx(ibctPeriodic)) then
              dVertexParameterPeriodic = DmaxParam(isegment)-&
                  DvertexParameterValue(ivbd)
            else
              dVertexParameterPeriodic = DmaxParam(isegmentPeriodic-1)+&
                  DmaxParam(isegment)-DvertexParameterValue(ivbd)
            end if

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add residual of the equation that corresponds to node
            ! ivt to the equation that corresponds to node ivtPeriodic
            ! and set residual of first equation equal to zero
            do ivar = 1, NVAR
              Dr(ivtPeriodic,ivar) = Dr(ivt,ivar)+Dr(ivtPeriodic,ivar)
              Dr(ivt,ivar)         = 0.0_DP
            end do

          case (BDRC_ANTIPERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            dVertexParameterPeriodic = p_DmaxPar(isegmentPeriodic)-&
                (p_DmaxPar(isegment)-p_DvertexParameterValue(ivbd))

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add residual of the equation that corresponds to node
            ! ivt to the equation that corresponds to node ivtPeriodic
            ! and set residual of first equation equal to zero
            do ivar = 1, NVAR
              Dr(ivtPeriodic,ivar) = Dr(ivt,ivar)+Dr(ivtPeriodic,ivar)
              Dr(ivt,ivar)         = 0.0_DP
            end do
            
          case default
            ! Get vertex number
            ivt  = IverticesAtBoundary(ivbd)
            ipos = Kdiagonal(ivt)

            ! Initialise variable values [x,y,0,time] for parser
            DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
            DvariableValues(NDIM2D)   = DvertexCoords(2,ivt)
            DvariableValues(NDIM3D)   = 0.0_DP
            DvariableValues(NDIM3D+1) = ttime

            ! Check for user-defined callback-function
            if (present(fcb_calcBoundaryvalues)) then

              ! Predict boundary values
              do ivar = 1, nvar
                Daux(ivar) = Du(ivt,ivar)+Dr(ivt,ivar)/rarray(ivar,ivar)%Da(ipos)
              end do
              
              ! Get values from old solution
              do ivar = 1, nvar
                Daux0(ivar) = Du0(ivt,ivar)
              end do

              ! Get desired boundary values from parser
              do iexpr = 1, nexpr
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                    DvariableValues, DbdrValues(iexpr))
              end do
              
              ! Get vertex number of predecessor
              if (ivbd .eq. ivbdFirst) then
                ivtL = IverticesAtBoundary(ivbdLast)
              else
                ivtL = IverticesAtBoundary(ivbd-1)
              end if
              
              ! Get vertex number of sucessor
              if (ivbd .eq. ivbdLast) then
                ivtR = IverticesAtBoundary(ivbdFirst)
              else
                ivtR = IverticesAtBoundary(ivbd+1)
              end if
              
              ! Compute the analytical and approximate normal vectors
              if (DvertexParameterValue(ivbd) .eq. dminValue) then
                ! We are at the beginning of a segment that includes
                ! the starting point.  Otherwise, the segment counter
                ! would not have been increased.
                
                ! Thus, compute the analytical normal vector from the right
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_RIGHT)
                
                ! Compute the approximate normal vector from the right element
                dnx = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dny = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              elseif (DvertexParameterValue(ivbd) .eq. dmaxValue) then
                ! We are at the end of a segment that includes the end
                ! point.  Otherwise, the segment counter would have
                ! been increased already.

                ! Thus, compute the analytical normal vector from the left.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_LEFT)
                
                ! Compute the approximate normal vector from the left element
                dnx = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dny = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              else
                ! We are "inside" a segment.  Thus, compute the mean
                ! value of left and right normal vector.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2))
                
                ! Compute the approximate normal vector
                dnxL = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dnyL = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                
                dnxR = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dnyR = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                
                dwL = sqrt( dnxL*dnxL + dnyL*dnyL )
                dwR = sqrt( dnxR*dnxR + dnyR*dnyR )
                
                dnx = (dnxL/dwL + dnxR/dwR)/(1/dwL + 1/dwR)
                dny = (dnyL/dwL + dnyR/dwR)/(1/dwL + 1/dwR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
              end if
              
              ! Apply callback function to determine boundary conditions
              call fcb_calcBoundaryvalues(DbdrNormal, DpointNormal, DbdrValues,&
                  IbdrCondType(isegment), Daux, Daux0, istatus)
              do ivar = 1, nvar
                Du(ivt,ivar) = Daux(ivar)
              end do

            else
              
              ! Impose prescribed boundary conditions directly
              do ivar = 1, nvar
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                    DvariableValues, Du(ivt,ivar))
              end do

            end if
              
            ! Nullify residual entries
            do ivar = 1, nvar
              Dr(ivt,ivar) = 0.0_DP
            end do

          end select
        end do
      end do
    end subroutine filtersolution_Mat79_2D

  end subroutine bdrf_filterSolutionBlock

  ! ***************************************************************************

!<subroutine>

  subroutine bdrf_filterSolutionBlockScalar(rboundaryCondition,&
      rmatrix, rsolution, rdefect, rsolution0, ttime,&
      fcb_calcBoundaryvalues, istatus, rboundary, rtriangulation)

!<description>
    ! This subroutine imposes prescribed boundary conditions by
    ! filtering the solution vector in the following way:
    !
    ! 1. Compute a provisional solution by scaling the defect vector
    !    by the inverse of the diagonal matrix entries and adding
    !    the result to the solution from the last iteration:
    !
    !    $$U^*_i=U_i+A_{ii}^{-1}R_i,\qquad \mbox{whereby}\quad
    !    A_{ii}=\mbox{diag}\{a_{ii}^{kk}$$
    !
    ! 2. Modify the predicted solution values by means of
    !    an external subroutine or impose the prescribed Dirichlet
    !    boundary values
    !
    !    $$U_i^{**}=\mbox{bc}(U_i^*)$$
    !
    ! 3. Nullify all entries of the defect vector $R_i=0$
    !
    ! Then the increment $\Delta U_i=A_{ii}^{-1}R_i=0$ so that
    ! for the updated solution vector $U_i=U_i^{**}$.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! The initial solution vector from the last time step
    type(t_vectorBlock), intent(in) :: rsolution0

    ! The simulation time
    real(DP), intent(in) :: ttime

    ! OPTIONAL: The boundary description
    type(t_boundary), intent(in), optional :: rboundary

    ! OPTIONAL: The triangulation
    type(t_triangulation), intent(in), optional :: rtriangulation

    ! OPTIONAL: The callback function
    include 'intf_bdrcallback.inc'
    optional :: fcb_calcBoundaryvalues
!</input>

!<inputoutput>
    ! The matrix data structure
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! The solution vector that should be filtered
    type(t_vectorBlock), intent(inout) :: rsolution

    ! The defect vector that should be filtered
    type(t_vectorBlock), intent(inout) :: rdefect

    ! OPTIONAL: Status of the callback function
    integer, intent(inout), optional :: istatus
!</inputoutput>
!</subroutine>

    ! Initialise status
    if (present(istatus)) istatus = 0

    ! Check if block vectors have only one block
    if ((rsolution%nblocks .eq. 1) .and.&
        (rdefect%nblocks   .eq. 1)) then
      call bdrf_filterSolutionScalar(rboundaryCondition, rmatrix,&
          rsolution%RvectorBlock(1), rdefect%RvectorBlock(1),&
          rsolution0%RvectorBlock(1), ttime, &
          fcb_calcBoundaryvalues, istatus, rboundary, rtriangulation)
    else
      call output_line('System matrix and vectors are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionBlockScalar')
      call sys_halt()
    end if
  end subroutine bdrf_filterSolutionBlockScalar

  ! ***************************************************************************

!<subroutine>

  subroutine bdrf_filterSolutionScalar(rboundaryCondition, rmatrix,&
      rsolution, rdefect, rsolution0, ttime, fcb_calcBoundaryvalues,&
      istatus, rboundary, rtriangulation)

!<description>
    ! This subroutine imposes prescribed boundary conditions by
    ! filtering the solution vector in the following way:
    !
    ! 1. Compute a provisional solution by scaling the defect vector
    !    by the inverse of the diagonal matrix entries and adding
    !    the result to the solution from the last iteration:
    !
    !    $$U^*_i=U_i+A_{ii}^{-1}R_i,\qquad \mbox{whereby}\quad
    !    A_{ii}=\mbox{diag}\{a_{ii}^{kk}$$
    !
    ! 2. Modify the predicted solution values by means of
    !    an external subroutine or impose the prescribed Dirichlet
    !    boundary values
    !
    !    $$U_i^{**}=\mbox{bc}(U_i^*)$$
    !
    ! 3. Nullify all entries of the defect vector $R_i=0$
    !
    ! Then the increment $\Delta U_i=A_{ii}^{-1}R_i=0$ so that
    ! for the updated solution vector $U_i=U_i^{**}$.
!</description>

!<input>
    ! The boundary conditions
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! The initial solution vector from the last time step
    type(t_vectorScalar), intent(in) :: rsolution0

    ! The simulation time
    real(DP), intent(in) :: ttime

    ! OPTIONAL: The boundary description
    type(t_boundary), intent(in), optional, target :: rboundary

    ! The triangulation
    type(t_triangulation), intent(in), optional, target :: rtriangulation

    ! OPTIONAL: The callback function
    include 'intf_bdrcallback.inc'
    optional :: fcb_calcBoundaryvalues
!</input>

!<inputoutput>
    ! The scalar matrix
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! The scalar solution vector that should be filtered
    type(t_vectorScalar), intent(inout) :: rsolution

    ! The scalar defect vector that should be filtered
    type(t_vectorScalar), intent(inout) :: rdefect

    ! OPTIONAL: Status of the callback function
    integer, intent(inout), optional :: istatus
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundary), pointer :: p_rboundary
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DmaxPar
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    real(DP), dimension(:), pointer :: p_DA, p_Du, p_Dr, p_Du0
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    integer, dimension(:), pointer :: p_Kcol
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IbdrCondCpIdx
    integer, dimension(:), pointer :: p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCompPeriodic
    integer, dimension(:), pointer :: p_IbdrCondPeriodic
    logical, dimension(:), pointer :: p_BisSegClosed

    logical :: bisMatrixSorted,bisSolutionSorted,bisDefectSorted

    ! Initialise status
    if (present(istatus)) istatus = 0

    ! Check if there are strong boundary conditions
    if (.not.rboundaryCondition%bStrongBdrCond) return

    ! Check if matrix is sorted?
    if (rmatrix%isortStrategy .gt. 0) then
      bisMatrixSorted = .true.
      call lsyssc_unsortMatrix(rmatrix, .true.)
    else
      bisMatrixSorted = .false.
    end if

    ! Check if solution vector is sorted?
    if (rsolution%isortStrategy .gt. 0) then
      bisSolutionSorted = .true.
      call lsyssc_vectorActivateSorting(rsolution,.false.)
    else
      bisSolutionSorted = .false.
    end if

    ! Check if defect vector is sorted?
    if (rdefect%isortStrategy .gt. 0) then
      bisDefectSorted = .true.
      call lsyssc_vectorActivateSorting(rdefect,.false.)
    else
      bisDefectSorted = .false.
    end if

    ! Set pointers for solution and defect
    call lsyssc_getbase_double(rsolution0, p_Du0)
    call lsyssc_getbase_double(rsolution, p_Du)
    call lsyssc_getbase_double(rdefect, p_Dr)

    ! Get underlying triangulation structure
    if (present(rtriangulation)) then
      p_rtriangulation => rtriangulation
    else
      if (.not.associated(rsolution%p_rspatialDiscr)) then
        call output_line('No discretisation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionScalar')
        call sys_halt()
      end if

      if (.not.associated(rsolution%p_rspatialDiscr%p_rtriangulation)) then
        call output_line('No triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionScalar')
        call sys_halt()
      end if
      p_rtriangulation => rsolution%p_rspatialDiscr%p_rtriangulation
    end if

    ! Check spatial dimensions
    if (p_rtriangulation%ndim .ne. rboundaryCondition%ndimension) then
      call output_line('Spatial dimension mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionScalar')
      call sys_halt()
    end if

    ! How many spatial dimensions do we have?
    select case (rboundaryCondition%ndimension)
    case (NDIM1D)
      ! Set pointers for triangulation
      call storage_getbase_double2d(&
          p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
      call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic, p_IbdrCompPeriodic)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic, p_IbdrCondPeriodic)

      ! What kind of matrix are we?
      select case (rmatrix%cmatrixFormat)

      case (LSYSSC_MATRIXD)
        call lsyssc_getbase_double(rmatrix, p_DA)

        ! Set prescribed boundary values in 1D
        call filtersolution_MatD_1D(rboundaryCondition%rfparser,&
            p_IbdrCondType, p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
            rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
            p_IboundaryCpIdx, p_DvertexCoords, 1, p_DA, p_Du, p_Dr, p_Du0,&
            fcb_calcBoundaryvalues, istatus)

      case (LSYSSC_MATRIX7)
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        call lsyssc_getbase_double(rmatrix, p_DA)

        ! Set prescribed boundary values in 1D
        call filtersolution_Mat79IntlD_1D(rboundaryCondition%rfparser,&
            p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
            p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
            rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
            p_IboundaryCpIdx, p_DvertexCoords, p_Kld, p_Kcol, p_Kld,&
            1, p_DA, p_Du, p_Dr, p_Du0, fcb_calcBoundaryvalues, istatus)

      case (LSYSSC_MATRIX9)
        call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        call lsyssc_getbase_double(rmatrix, p_DA)

        ! Set prescribed boundary values in 1D
        call filtersolution_Mat79IntlD_1D(rboundaryCondition%rfparser,&
            p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
            p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
            rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
            p_IboundaryCpIdx, p_DvertexCoords, p_Kld, p_Kcol, p_Kdiagonal,&
            1, p_DA, p_Du, p_Dr, p_Du0, fcb_calcBoundaryvalues, istatus)

      case (LSYSSC_MATRIX7INTL)
        select case (rmatrix%cinterleavematrixFormat)

        case (LSYSSC_MATRIXD)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          ! Set prescribed boundary values in 1D
          call filtersolution_Mat79IntlD_1D(rboundaryCondition%rfparser,&
              p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
              p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
              rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
              p_IboundaryCpIdx, p_DvertexCoords, p_Kld, p_Kcol, p_Kld,&
              rmatrix%NVAR, p_DA, p_Du, p_Dr, p_Du0, fcb_calcBoundaryvalues, istatus)

        case (LSYSSC_MATRIX1)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          ! Set prescribed boundary values in 1D
          call filtersolution_Mat79Intl1_1D(rboundaryCondition%rfparser,&
              p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
              p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
              rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
              p_IboundaryCpIdx, p_DvertexCoords, p_Kld, p_Kcol, p_Kld,&
              rmatrix%NVAR, p_DA, p_Du, p_Dr, p_Du0, fcb_calcBoundaryvalues, istatus)

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_systemScalar')
          call sys_halt()
        end select

      case (LSYSSC_MATRIX9INTL)
        select case (rmatrix%cinterleavematrixFormat)

        case (LSYSSC_MATRIXD)
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          ! Set prescribed boundary values in 1D
          call filtersolution_Mat79IntlD_1D(rboundaryCondition%rfparser,&
              p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
              p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
              rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
              p_IboundaryCpIdx, p_DvertexCoords, p_Kld, p_Kcol, p_Kdiagonal,&
              rmatrix%NVAR, p_DA, p_Du, p_Dr, p_Du0, fcb_calcBoundaryvalues, istatus)

        case (LSYSSC_MATRIX1)
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          ! Set prescribed boundary values in 1D
          call filtersolution_Mat79Intl1_1D(rboundaryCondition%rfparser,&
              p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
              p_IbdrCondCpIdx, rboundaryCondition%iboundarycount,&
              rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
              p_IboundaryCpIdx, p_DvertexCoords, p_Kld, p_Kcol, p_Kdiagonal,&
              rmatrix%NVAR, p_DA, p_Du, p_Dr, p_Du0, fcb_calcBoundaryvalues, istatus)

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionScalar')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionScalar')
        call sys_halt()
      end select


    case (NDIM2D)
      ! Set pointers for triangulation
      call storage_getbase_double2d(&
          p_rtriangulation%h_DvertexCoords, p_DvertexCoords)
      call storage_getbase_double(&
          p_rtriangulation%h_DvertexParameterValue, p_DvertexParameterValue)
      call storage_getbase_int(&
          p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
      call storage_getbase_int(&
          p_rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)

      ! Set pointers for boundary
      call storage_getbase_double(rboundaryCondition%h_DmaxPar, p_DmaxPar)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx, p_IbdrCondCpIdx)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
      call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic, p_IbdrCompPeriodic)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic, p_IbdrCondPeriodic)
      call storage_getbase_logical(rboundaryCondition%h_BisSegClosed, p_BisSegClosed)

      ! Get underlying boundary structure
      if (present(rboundary)) then
        p_rboundary => rboundary
      else
        if (.not.associated(rsolution%p_rspatialDiscr%p_rboundary)) then
          call output_line('No boundary associated!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionScalar')
          call sys_halt()
        end if
        p_rboundary => rsolution%p_rspatialDiscr%p_rboundary
      end if

      ! What kind of matrix are we?
      select case (rmatrix%cmatrixFormat)

      case (LSYSSC_MATRIXD)
        call lsyssc_getbase_double(rmatrix, p_DA)

        ! Set prescribed boundary values in 2D
        call filtersolution_MatD_2D(rboundaryCondition%rfparser,&
            p_IbdrCondType, p_IbdrCondCpIdx, p_DmaxPar,&
            p_BisSegClosed, rboundaryCondition%iboundarycount,&
            rboundaryCondition%nmaxExpressions, p_IverticesAtBoundary,&
            p_IboundaryCpIdx, p_DvertexParameterValue, p_DvertexCoords, 1,&
            p_DA, p_Du, p_Dr, p_Du0, rboundary, fcb_calcBoundaryvalues, istatus)

      case (LSYSSC_MATRIX7)
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        call lsyssc_getbase_double(rmatrix, p_DA)

        ! Set prescribed boundary values in 2D
        call filtersolution_Mat79IntlD_2D(rboundaryCondition%rfparser,&
            p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
            p_IbdrCondCpIdx, p_DmaxPar, p_BisSegClosed,&
            rboundaryCondition%iboundarycount, rboundaryCondition%nmaxExpressions,&
            p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
            p_DvertexCoords, p_Kld, p_Kcol, p_Kld, 1, p_DA, p_Du, p_Dr, p_Du0,&
            rboundary, fcb_calcBoundaryvalues, istatus)

      case (LSYSSC_MATRIX9)
        call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        call lsyssc_getbase_double(rmatrix, p_DA)

        ! Set prescribed boundary values in 2D
        call filtersolution_Mat79IntlD_2D(rboundaryCondition%rfparser,&
            p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
            p_IbdrCondCpIdx, p_DmaxPar, p_BisSegClosed,&
            rboundaryCondition%iboundarycount, rboundaryCondition%nmaxExpressions,&
            p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
            p_DvertexCoords, p_Kld, p_Kcol, p_Kdiagonal, 1, p_DA, p_Du, p_Dr, p_Du0,&
            rboundary, fcb_calcBoundaryvalues, istatus)

      case (LSYSSC_MATRIX7INTL)
        select case (rmatrix%cinterleavematrixFormat)

        case (LSYSSC_MATRIXD)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          ! Set prescribed boundary values in 2D
          call filtersolution_Mat79IntlD_2D(rboundaryCondition%rfparser,&
              p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
              p_IbdrCondCpIdx, p_DmaxPar, p_BisSegClosed,&
              rboundaryCondition%iboundarycount, rboundaryCondition%nmaxExpressions,&
              p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
              p_DvertexCoords, p_Kld, p_Kcol, p_Kld, rmatrix%NVAR, p_DA, p_Du,&
              p_Dr, p_Du0, p_rboundary, fcb_calcBoundaryvalues, istatus)

        case (LSYSSC_MATRIX1)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          ! Set prescribed boundary values in 2D
          call filtersolution_Mat79Intl1_2D(rboundaryCondition%rfparser,&
              p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
              p_IbdrCondCpIdx, p_DmaxPar, p_BisSegClosed,&
              rboundaryCondition%iboundarycount, rboundaryCondition%nmaxExpressions,&
              p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
              p_DvertexCoords, p_Kld, p_Kcol, p_Kld, rmatrix%NVAR, p_DA, p_Du,&
              p_Dr, p_Du0, p_rboundary, fcb_calcBoundaryvalues, istatus)

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_systemScalar')
          call sys_halt()
        end select

      case (LSYSSC_MATRIX9INTL)
        select case (rmatrix%cinterleavematrixFormat)

        case (LSYSSC_MATRIXD)
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          ! Set prescribed boundary values in 2D
          call filtersolution_Mat79IntlD_2D(rboundaryCondition%rfparser,&
              p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
              p_IbdrCondCpIdx, p_DmaxPar, p_BisSegClosed,&
              rboundaryCondition%iboundarycount, rboundaryCondition%nmaxExpressions,&
              p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
              p_DvertexCoords, p_Kld, p_Kcol, p_Kdiagonal, rmatrix%NVAR, p_DA,&
              p_Du, p_Dr, p_Du0, p_rboundary, fcb_calcBoundaryvalues, istatus)

        case (LSYSSC_MATRIX1)
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
          call lsyssc_getbase_double(rmatrix, p_DA)

          ! Set prescribed boundary values in 2D
          call filtersolution_Mat79Intl1_2D(rboundaryCondition%rfparser,&
              p_IbdrCompPeriodic, p_IbdrCondPeriodic, p_IbdrCondType,&
              p_IbdrCondCpIdx, p_DmaxPar, p_BisSegClosed,&
              rboundaryCondition%iboundarycount, rboundaryCondition%nmaxExpressions,&
              p_IverticesAtBoundary, p_IboundaryCpIdx, p_DvertexParameterValue,&
              p_DvertexCoords, p_Kld, p_Kcol, p_Kdiagonal, rmatrix%NVAR, p_DA,&
              p_Du, p_Dr, p_Du0, p_rboundary, fcb_calcBoundaryvalues, istatus)

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionScalar')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionScalar')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdrf_filterSolutionScalar')
      call sys_halt()
    end select

    ! Do we have to re-sort matrix?
    if (bisMatrixSorted) call lsyssc_sortMatrix(rmatrix,.true.,rmatrix%isortStrategy)

    ! Do we have to re-sort solution vector?
    if (bisSolutionSorted) call lsyssc_vectorActivateSorting(rsolution,.true.)

    ! Do we have to re-sort defect vector
    if (bisDefectSorted) call lsyssc_vectorActivateSorting(rdefect,.true.)

  contains

    ! Here are the real filter routines

    !***************************************************************
    ! Filter solution and defect vector in 1D.
    ! Here, the matrix is given as diagonal matrix.
    ! This subroutine can only handle Neumann and Dirichlet boundarys.

    subroutine filtersolution_MatD_1D(rfparser, IbdrCondType,&
        IbdrCondCpIdx, nbct, nexpr, IverticesAtBoundary,&
        IboundaryCpIdx, DvertexCoords, nvar, DA, Du, Dr, Du0,&
        fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Number of variables
      integer, intent(in) :: nvar

      ! Matrix data array
      real(DP), dimension(nvar,*), intent(in) :: DA

      ! Solution and residual data array
      real(DP), dimension(nvar,*), intent(inout) :: Du, Dr

      ! Initial solution data array
      real(DP), dimension(nvar,*), intent(in) :: Du0

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(NDIM1D) :: DbdrNormal
      integer :: ivbd,ivt,ibct,isegment,iexpr,ivar

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          call output_line('Unable to handle periodic boundary ' // &
              'conditions for diagonal matrix!',&
              OU_CLASS_WARNING, OU_MODE_STD, 'filtersolution_MatD_1D')
          call sys_halt()
          
        case default
          ! Get vertex number
          ivt = IverticesAtBoundary(ivbd)

          ! Initialise variable values [x,0,0,time] for parser
          DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
          DvariableValues(NDIM2D)   = 0.0_DP
          DvariableValues(NDIM3D)   = 0.0_DP
          DvariableValues(NDIM3D+1) = ttime

          ! Check for user-defined callback-function
          if (present(fcb_calcBoundaryvalues)) then
            
            ! Predict boundary values
            Du(:,ivt) = Du(:,ivt)+Dr(:,ivt)/Da(:,ivt)
            
            ! Get desired boundary values from parser
            do iexpr = 1, nexpr
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                  DvariableValues, DbdrValues(iexpr))
            end do
            
            ! Compute the analytical normal vector
            if (mod(ibct, 2) .eq. 0) then
              DbdrNormal = 1.0_DP
            else
              DbdrNormal = -1.0_DP
            end if
            
            ! Apply callback function to determine boundary conditions
            call fcb_calcBoundaryvalues(DbdrNormal, DbdrNormal, DbdrValues,&
                IbdrCondType(isegment), Du(:,ivt), Du0(:,ivt), istatus)

          else

            ! Impose prescribed boundary conditions directly
            do ivar = 1, nvar
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                  DvariableValues, Du(ivar,ivt))
            end do

          end if

          ! Nullify residual entry
          Dr(:,ivt) = 0.0_DP

        end select
      end do
    end subroutine filtersolution_MatD_1D


    !***************************************************************
    ! Filter solution and defect vector in 2D.
    ! Here, the matrix is given as diagonal matrix.
    ! This subroutine can only handle Neumann and Dirichlet boundarys.

    subroutine filtersolution_MatD_2D(rfparser, IbdrCondType,&
        IbdrCondCpIdx, DmaxParam, BisSegClosed, nbct, nexpr,&
        IverticesAtBoundary, IboundaryCpIdx, DvertexParameterValue,&
        DvertexCoords, nvar, DA, Du, Dr, Du0, rboundary,&
        fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam
      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Number of variables
      integer, intent(in) :: nvar

      ! Matrix data array
      real(DP), dimension(nvar,*), intent(in) :: DA

      ! Solution and residual data array
      real(DP), dimension(nvar,*), intent(inout) :: Du, Dr

      ! Initial solution data array
      real(DP), dimension(nvar,*), intent(in) :: Du0

      ! Boundary description
      type(t_boundary), intent(in) :: rboundary

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(NDIM2D) :: DbdrNormal,DpointNormal
      real(DP) :: dminValue,dmaxValue
      real(DP) :: dnx,dny,dnxL,dnxR,dnyL,dnyR,dw,dwL,dwR
      integer :: ivbd,ivbdFirst,ivbdLast,ivt,ibct,ivtL,ivtR,isegment,iexpr,ivar


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue = nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle
          
          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))
            
          case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
            call output_line('Unable to handle periodic boundary ' // &
                'conditions for diagonal matrix!',&
                OU_CLASS_WARNING, OU_MODE_STD, 'filtersolution_MatD_2D')
            call sys_halt()
            
          case default
            ! Get vertex number
            ivt = IverticesAtBoundary(ivbd)

            ! Initialise variable values [x,y,0,time] for parser
            DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
            DvariableValues(NDIM2D)   = DvertexCoords(2,ivt)
            DvariableValues(NDIM3D)   = 0.0_DP
            DvariableValues(NDIM3D+1) = ttime

            ! Check for user-defined callback-function
            if (present(fcb_calcBoundaryvalues)) then

              ! Predict boundary values
              Du(:,ivt) = Du(:,ivt)+Dr(:,ivt)/Da(:,ivt)
              
              ! Get desired boundary values from parser
              do iexpr = 1, nexpr
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                    DvariableValues, DbdrValues(iexpr))
              end do
              
              ! Get vertex number of predecessor
              if (ivbd .eq. ivbdFirst) then
                ivtL = IverticesAtBoundary(ivbdLast)
              else
                ivtL = IverticesAtBoundary(ivbd-1)
              end if
              
              ! Get vertex number of sucessor
              if (ivbd .eq. ivbdLast) then
                ivtR = IverticesAtBoundary(ivbdFirst)
              else
                ivtR = IverticesAtBoundary(ivbd+1)
              end if
              
              ! Compute the analytical and approximate normal vectors
              if (DvertexParameterValue(ivbd) .eq. dminValue) then
                ! We are at the beginning of a segment that includes
                ! the starting point.  Otherwise, the segment counter
                ! would not have been increased.

                ! Thus, compute the analytical normal vector from the right
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_RIGHT)
                
                ! Compute the approximate normal vector from the right element
                dnx = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dny = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              elseif (DvertexParameterValue(ivbd) .eq. dmaxValue) then
                ! We are at the end of a segment that includes the end
                ! point.  Otherwise, the segment counter would have
                ! been increased already.
                
                ! Thus, compute the analytical normal vector from the left.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_LEFT)
                
                ! Compute the approximate normal vector from the left element
                dnx = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dny = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              else
                ! We are "inside" a segment.  Thus, compute the mean
                ! value of left and right normal vector.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2))
                
                ! Compute the approximate normal vector
                dnxL = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dnyL = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                
                dnxR = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dnyR = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                
                dwL = sqrt( dnxL*dnxL + dnyL*dnyL )
                dwR = sqrt( dnxR*dnxR + dnyR*dnyR )
                
                dnx = (dnxL/dwL + dnxR/dwR)/(1/dwL + 1/dwR)
                dny = (dnyL/dwL + dnyR/dwR)/(1/dwL + 1/dwR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
              end if
              
              ! Apply callback function to determine boundary conditions
              call fcb_calcBoundaryvalues(DbdrNormal, DpointNormal, DbdrValues,&
                  IbdrCondType(isegment), Du(:,ivt), Du0(:,ivt), istatus)

            else

              ! Impose prescribed Dirichlet boundary conditions
              do ivar = 1, nvar
                call fparser_evalFunction(rfparser, nvar*(isegment-1)+ivar,&
                    DvariableValues, Du(ivar,ivt))
              end do

            end if

            ! Nullify residual entry
            Dr(:,ivt) = 0.0_DP
            
          end select
        end do
      end do
    end subroutine filtersolution_MatD_2D
    
    
    !***************************************************************
    ! Filter solution and defect vector in 1D.
    ! Here, the matrix is store in CSR7 and CSR9 interleave format,
    ! whereby each local matrix is a diagonal matrix.

    subroutine filtersolution_Mat79IntlD_1D(rfparser,&
        IbdrCompPeriodic, IbdrCondPeriodic, IbdrCondType,&
        IbdrCondCpIdx, nbct, nexpr, IverticesAtBoundary,&
        IboundaryCpIdx, DvertexCoords, Kld, Kcol, Kdiagonal, nvar,&
        DA, Du, Dr, Du0, fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Number of variables
      integer, intent(in) :: nvar

      ! Matrix data array
      real(DP), dimension(nvar,*), intent(in) :: DA

      ! Solution and residual data array
      real(DP), dimension(nvar,*), intent(inout) :: Du, Dr

      ! Initial solution data array
      real(DP), dimension(nvar,*), intent(in) :: Du0

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(NDIM1D) :: DbdrNormal
      integer :: ivbd,ivbdPeriodic,ivt,ivtPeriodic
      integer :: ipos,ibct,ibctPeriodic,ivar,isegment,iexpr

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          ! Compute vertex parameter value at periodic boundary
          ibctPeriodic = IbdrCompPeriodic(isegment)
          ivbdPeriodic = IboundaryCpIdx(ibctPeriodic)

          ! Get numbers of vertices
          ivt         = IverticesAtBoundary(ivbd)
          ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

          ! Add residual of the equation that corresponds to node ivt to the equation that
          ! corresponds to node ivtPeriodic and set residual of first equation equal to zero
          Dr(:,ivtPeriodic) = Dr(:,ivt)+Dr(:,ivtPeriodic)
          Dr(:,ivt)         = 0.0_DP
          
        case default
          ! Get vertex number
          ivt  = IverticesAtBoundary(ivbd)
          
          ! Initialise variable values [x,0,0,time] for parser
          DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
          DvariableValues(NDIM2D)   = 0.0_DP
          DvariableValues(NDIM3D)   = 0.0_DP
          DvariableValues(NDIM3D+1) = ttime

          ! Check for user-defined callback-function
          if (present(fcb_calcBoundaryvalues)) then
            
            ! Predict boundary values
            ipos = Kdiagonal(ivt)
            Du(:,ivt) = Du(:,ivt)+Dr(:,ivt)/Da(:,ipos)
            
            ! Get desired boundary values from parser
            do iexpr = 1, nexpr
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                  DvariableValues, DbdrValues(iexpr))
            end do
            
            ! Compute the analytical normal vector
            if (mod(ibct, 2) .eq. 0) then
              DbdrNormal = 1.0_DP
            else
              DbdrNormal = -1.0_DP
            end if
            
            ! Apply callback function to determine boundary conditions
            call fcb_calcBoundaryvalues(DbdrNormal, DbdrNormal, DbdrValues,&
                IbdrCondType(isegment), Du(:,ivt), Du0(:,ivt), istatus)

          else

            ! Impose prescribed boundary conditions directly
            do ivar = 1, nvar
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                  DvariableValues, Du(ivar,ivt))
            end do

          end if
          
          ! Nullify residual entry
          Dr(:,ivt) = 0.0_DP

        end select
      end do
    end subroutine filtersolution_Mat79IntlD_1D


    !***************************************************************
    ! Filter solution and defect vector in 2D.
    ! Here, the matrix is store in CSR7 and CSR9 interleave format,
    ! whereby each local matrix is a diagonal matrix.

    subroutine filtersolution_Mat79IntlD_2D(rfparser,&
        IbdrCompPeriodic, IbdrCondPeriodic, IbdrCondType,&
        IbdrCondCpIdx, DmaxParam, BisSegClosed, nbct, nexpr,&
        IverticesAtBoundary, IboundaryCpIdx, DvertexParameterValue,&
        DvertexCoords, Kld, Kcol, Kdiagonal, nvar, DA, Du, Dr, Du0,&
        rboundary, fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Number of variables
      integer, intent(in) :: nvar

      ! Matrix data array
      real(DP), dimension(nvar,*), intent(in) :: DA

      ! Solution and residual data array
      real(DP), dimension(nvar,*), intent(inout) :: Du, Dr

      ! Initial solution data array
      real(DP), dimension(nvar,*), intent(in) :: Du0

      ! Boundary description
      type(t_boundary), intent(in) :: rboundary

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(NDIM2D) :: DbdrNormal,DpointNormal
      real(DP) :: dminValue,dmaxValue,dVertexParameterPeriodic
      real(DP) :: dnx,dny,dnxL,dnxR,dnyL,dnyR,dw,dwL,dwR
      integer :: ivbd,ivbdFirst,ivbdLast,ivbdPeriodic,ivt,ivtPeriodic,ivtL,ivtR
      integer :: ipos,ibct,ibctPeriodic,ivar,isegment,isegmentPeriodic,iexpr

      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue = nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle
          
          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

          case (BDRC_PERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            if (isegmentPeriodic .eq. IbdrCondCpIdx(ibctPeriodic)) then
              dVertexParameterPeriodic = DmaxParam(isegment)-&
                  DvertexParameterValue(ivbd)
            else
              dVertexParameterPeriodic = DmaxParam(isegmentPeriodic-1)+&
                  DmaxParam(isegment)-DvertexParameterValue(ivbd)
            end if

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add residual of the equation that corresponds to node
            ! ivt to the equation that corresponds to node ivtPeriodic
            ! and set residual of first equation equal to zero
            Dr(:,ivtPeriodic) = Dr(:,ivt)+Dr(:,ivtPeriodic)
            Dr(:,ivt)         = 0.0_DP

          case (BDRC_ANTIPERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            dVertexParameterPeriodic = p_DmaxPar(isegmentPeriodic)-&
                (p_DmaxPar(isegment)-p_DvertexParameterValue(ivbd))

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add residual of the equation that corresponds to node
            ! ivt to the equation that corresponds to node ivtPeriodic
            ! and set residual of first equation equal to zero
            Dr(:,ivtPeriodic) = Dr(:,ivt)+Dr(:,ivtPeriodic)
            Dr(:,ivt)         = 0.0_DP
            
          case default
            ! Get vertex number
            ivt = IverticesAtBoundary(ivbd)

            ! Initialise variable values [x,y,0,time] for parser
            DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
            DvariableValues(NDIM2D)   = DvertexCoords(2,ivt)
            DvariableValues(NDIM3D)   = DvertexCoords(2,ivt)
            DvariableValues(NDIM3D+1) = ttime

            ! Check for user-defined callback-function
            if (present(fcb_calcBoundaryvalues)) then

              ! Predict boundary values
              ipos = Kdiagonal(ivt)
              Du(:,ivt) = Du(:,ivt)+Dr(:,ivt)/Da(:,ipos)
              
              ! Get desired boundary values from parser
              do iexpr = 1, nexpr
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                    DvariableValues, DbdrValues(iexpr))
              end do
              
              ! Get vertex number of predecessor
              if (ivbd .eq. ivbdFirst) then
                ivtL = IverticesAtBoundary(ivbdLast)
              else
                ivtL = IverticesAtBoundary(ivbd-1)
              end if
              
              ! Get vertex number of sucessor
              if (ivbd .eq. ivbdLast) then
                ivtR = IverticesAtBoundary(ivbdFirst)
              else
                ivtR = IverticesAtBoundary(ivbd+1)
              end if
              
              ! Compute the analytical and approximate normal vectors
              if (DvertexParameterValue(ivbd) .eq. dminValue) then
                ! We are at the beginning of a segment that includes
                ! the starting point.  Otherwise, the segment counter
                ! would not have been increased.

                ! Thus, compute the analytical normal vector from the right
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_RIGHT)
                
                ! Compute the approximate normal vector from the right element
                dnx = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dny = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              elseif (DvertexParameterValue(ivbd) .eq. dmaxValue) then
                ! We are at the end of a segment that includes the end
                ! point.  Otherwise, the segment counter would have
                ! been increased already.
                
                ! Thus, compute the analytical normal vector from the left.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_LEFT)
                
                ! Compute the approximate normal vector from the left element
                dnx = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dny = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              else
                ! We are "inside" a segment.  Thus, compute the mean
                ! value of left and right normal vector.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2))
                
                ! Compute the approximate normal vector
                dnxL = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dnyL = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                
                dnxR = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dnyR = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                
                dwL = sqrt( dnxL*dnxL + dnyL*dnyL )
                dwR = sqrt( dnxR*dnxR + dnyR*dnyR )
                
                dnx = (dnxL/dwL + dnxR/dwR)/(1/dwL + 1/dwR)
                dny = (dnyL/dwL + dnyR/dwR)/(1/dwL + 1/dwR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
              end if
              
              ! Apply callback function to determine boundary conditions
              call fcb_calcBoundaryvalues(DbdrNormal, DpointNormal, DbdrValues,&
                  IbdrCondType(isegment), Du(:,ivt), Du0(:,ivt), istatus)

            else

              ! Impose prescribed boundary conditions directly
              do ivar = 1, nvar
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                    DvariableValues, Du(ivar,ivt))
              end do
              
            end if
            
            ! Nullify residual entry
            Dr(:,ivt) = 0.0_DP

          end select
        end do
      end do
    end subroutine filtersolution_Mat79IntlD_2D


    !***************************************************************
    ! Filter solution and defect vector in 1D.
    ! Here, the matrix is store in CSR7 and CSR9 interleave format,
    ! whereby each local matrix is a full matrix.

    subroutine filtersolution_Mat79Intl1_1D(rfparser,&
        IbdrCompPeriodic, IbdrCondPeriodic, IbdrCondType,&
        IbdrCondCpIdx, nbct, nexpr, IverticesAtBoundary, IboundaryCpIdx,&
        DvertexCoords, Kld, Kcol, Kdiagonal, nvar, DA, Du, Dr, Du0,&
        fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Number of variables
      integer, intent(in) :: nvar

      ! Matrix data array
      real(DP), dimension(nvar,nvar,*), intent(in) :: DA

      ! Solution and residual data array
      real(DP), dimension(nvar,*), intent(inout) :: Du, Dr

      ! Initial solution data array
      real(DP), dimension(nvar,*), intent(in) :: Du0

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(NDIM1D) :: DbdrNormal,DpointNormal
      integer :: ivbd,ivbdPeriodic,ivt,ivtPeriodic
      integer :: ipos,ibct,ibctPeriodic,ivar,isegment,iexpr


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Get vertex of the boundary component
        ivbd = IboundaryCpIdx(ibct)

        ! Get first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct)

        ! Check if this segment has strong boundary conditions
        if (iand(int(p_IbdrCondType(isegment),I32),&
                 BDRC_STRONG) .ne. BDRC_STRONG) cycle

        ! What kind of boundary condition are we?
        select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))

        case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          ! Compute vertex parameter value at periodic boundary
          ibctPeriodic = IbdrCompPeriodic(isegment)
          ivbdPeriodic = IboundaryCpIdx(ibctPeriodic)

          ! Get numbers of vertices
          ivt         = IverticesAtBoundary(ivbd)
          ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

          ! Add residual of the equation that corresponds to node ivt
          ! to the equation that corresponds to node ivtPeriodic and
          ! set residual of first equation equal to zero
          Dr(:,ivtPeriodic) = Dr(:,ivt)+Dr(:,ivtPeriodic)
          Dr(:,ivt)         = 0.0_DP

        case default
          ! Get vertex number
          ivt = IverticesAtBoundary(ivbd)

          ! Initialise variable values [x,0,0,time] for parser
          DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
          DvariableValues(NDIM2D)   = 0.0_DP
          DvariableValues(NDIM3D)   = 0.0_DP
          DvariableValues(NDIM3D+1) = ttime

          ! Check for user-defined callback-function
          if (present(fcb_calcBoundaryvalues)) then

            ! Predict boundary values
            ipos = Kdiagonal(ivt)
            do ivar = 1, nvar
              Du(ivar,ivt) = Du(ivar,ivt)+Dr(ivar,ivt)/Da(ivar,ivar,ipos)
            end do
            
            ! Get desired boundary values from parser
            do iexpr = 1, nexpr
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                  DvariableValues, DbdrValues(iexpr))
            end do
            
            ! Compute the analytical normal vector
            if (mod(ibct, 2) .eq. 0) then
              DbdrNormal = 1.0_DP
            else
              DbdrNormal = -1.0_DP
            end if
            
            ! Apply callback function to determine boundary conditions
            call fcb_calcBoundaryvalues(DbdrNormal, DbdrNormal, DbdrValues,&
                IbdrCondType(isegment), Du(:,ivt), Du0(:,ivt), istatus)
            
          else

            ! Impose prescribed Dirichlet boundary conditions
            do ivar = 1, nvar
              call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                  DvariableValues, Du(ivar,ivt))
            end do

          end if

          ! Nullify residual entry
          Dr(:,ivt) = 0.0_DP

        end select
      end do
    end subroutine filtersolution_Mat79Intl1_1D
    

    !***************************************************************
    ! Filter solution and defect vector in 2D.
    ! Here, the matrix is store in CSR7 and CSR9 interleave format,
    ! whereby each local matrix is a full matrix.

    subroutine filtersolution_Mat79Intl1_2D(rfparser,&
        IbdrCompPeriodic, IbdrCondPeriodic, IbdrCondType,&
        IbdrCondCpIdx, DmaxParam, BisSegClosed, nbct, nexpr,&
        IverticesAtBoundary, IboundaryCpIdx, DvertexParameterValue,&
        DvertexCoords, Kld, Kcol, Kdiagonal, nvar, DA, Du, Dr, Du0,&
        rboundary, fcb_calcBoundaryvalues, istatus)

      ! Function parser used to evaluate Dirichlet boundary values
      type(t_fparser), intent(in) :: rfparser

      ! Array with periodic boundary components
      integer, dimension(:), intent(in) :: IbdrCompPeriodic

      ! Array with periodic boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondPeriodic

      ! Array with types of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondType

      ! Index array for type of boundary conditions
      integer, dimension(:), intent(in) :: IbdrCondCpIdx

      ! Array with maximum parameter value for each boundary component
      real(DP), dimension(:), intent(in) :: DmaxParam

      ! Array with booleans for the segment type
      logical, dimension(:), intent(in) :: BisSegClosed

      ! Number of boundary components
      integer, intent(in) :: nbct

      ! Number of maximum mathematical expressions
      integer, intent(in) :: nexpr

      ! Array with numbers of vertices at the boundary
      integer, dimension(:), intent(in) :: IverticesAtBoundary

      ! Index array for vertices at the boundary
      integer, dimension(:), intent(in) :: IboundaryCpIdx

      ! Array with parameter values for vertices at the boundary
      real(DP), dimension(:), intent(in) :: DvertexParameterValue

      ! Array with coordinates of vertices at the boundary
      real(DP), dimension(:,:), intent(in) :: DvertexCoords

      ! Array with matrix structure
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal

      ! Number of variables
      integer, intent(in) :: nvar

      ! Matrix data array
      real(DP), dimension(nvar,nvar,*), intent(in) :: DA

      ! Solution and residual data array
      real(DP), dimension(nvar,*), intent(inout) :: Du, Dr

      ! Initial solution data array
      real(DP), dimension(nvar,*), intent(in) :: Du0

      ! Boundary description
      type(t_boundary), intent(in) :: rboundary

      ! OPTIONAL: Callback function
      include 'intf_bdrcallback.inc'
      optional :: fcb_calcBoundaryvalues

      ! OPTIONAL: Status of the callback function
      integer, intent(inout), optional :: istatus


      ! local variables
      real(DP), dimension(NDIM3D+1) :: DvariableValues
      real(DP), dimension(nexpr) :: DbdrValues
      real(DP), dimension(NDIM2D) :: DbdrNormal,DpointNormal
      real(DP) :: dminValue,dmaxValue,dVertexParameterPeriodic
      real(DP) :: dnx,dny,dnxL,dnxR,dnyL,dnyR,dw,dwL,dwR
      integer :: ivbd,ivbdFirst,ivbdLast,ivbdPeriodic,ivt,ivtPeriodic,ivtL,ivtR
      integer :: ipos,ibct,ibctPeriodic,ivar,isegment,isegmentPeriodic,iexpr


      ! Loop over all boundary components
      do ibct = 1, nbct

        ! Set pointer to first/last vertex of the boundary component
        ivbdFirst = IboundaryCpIdx(ibct)
        ivbdLast  = IboundaryCpIdx(ibct+1)-1

        ! Set pointer to first region of the boundary component
        isegment  = IbdrCondCpIdx(ibct+1)-1
        dminValue = min(0._DP,&
            DmaxParam(isegment)-DvertexParameterValue(ivbdLast))
        dmaxValue = max(0._DP,&
            DvertexParameterValue(ivbdLast)-DmaxParam(isegment))

        ! Adjust endpoint parameter of segment
        if (.not.BisSegClosed(isegment)) dmaxValue = nearest(dmaxValue, -1._DP)

        ! Loop over all components of the boundary component
        do ivbd = ivbdFirst, ivbdLast

          ! Compute segment index
          do while(DvertexParameterValue(ivbd) .gt. dmaxValue)

            ! Adjust startpoint parameter of segment
            dminValue = nearest(dmaxValue, 1._DP)

            ! Increase segment index and reset it to first segment if required
            isegment = isegment+1
            if (isegment .gt. IbdrCondCpIdx(ibct+1)-1)&
                isegment = IbdrCondCpIdx(ibct)

            ! Adjust endpoint parameter of segment
            dmaxValue = DmaxParam(isegment)
            if (dmaxValue .lt. dminValue) dmaxValue =&
                ceiling(DvertexParameterValue(ivbdLast), DP)
            if (.not.BisSegClosed(isegment)) dmaxValue =&
                nearest(dmaxValue, -1._DP)
          end do

          ! Check if this segment has strong boundary conditions
          if (iand(int(p_IbdrCondType(isegment),I32),&
                   BDRC_STRONG) .ne. BDRC_STRONG) cycle

          ! What kind of boundary condition are we?
          select case(iand(int(IbdrCondType(isegment),I32), BDRC_TYPEMASK))
            
          case (BDRC_PERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            if (isegmentPeriodic .eq. IbdrCondCpIdx(ibctPeriodic)) then
              dVertexParameterPeriodic = DmaxParam(isegment)-&
                  DvertexParameterValue(ivbd)
            else
              dVertexParameterPeriodic = DmaxParam(isegmentPeriodic-1)+&
                  DmaxParam(isegment)-DvertexParameterValue(ivbd)
            end if

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add residual of the equation that corresponds to node
            ! ivt to the equation that corresponds to node ivtPeriodic
            ! and set residual of first equation equal to zero
            Dr(:,ivtPeriodic) = Dr(:,ivt)+Dr(:,ivtPeriodic)
            Dr(:,ivt)         = 0.0_DP

          case (BDRC_ANTIPERIODIC)
            ! Compute vertex parameter value at periodic boundary
            ibctPeriodic     = IbdrCompPeriodic(isegment)
            isegmentPeriodic = IbdrCondPeriodic(isegment)

            dVertexParameterPeriodic = p_DmaxPar(isegmentPeriodic)-&
                (p_DmaxPar(isegment)-p_DvertexParameterValue(ivbd))

            if (dVertexParameterPeriodic .eq.&
                DmaxParam(IbdrCondCpIdx(ibctPeriodic+1)-1))&
                dVertexParameterPeriodic = 0._DP

            ! Compute vertex number of nearest neighbor at boundary
            call bdrc_getNearestNeighbor2d(DvertexParameterValue,&
                dVertexParameterPeriodic, ivbdFirst, ivbdLast, ivbdPeriodic)

            ! Get numbers of vertices
            ivt         = IverticesAtBoundary(ivbd)
            ivtPeriodic = IverticesAtBoundary(ivbdPeriodic)

            ! Add residual of the equation that corresponds to node
            ! ivt to the equation that corresponds to node ivtPeriodic
            ! and set residual of first equation equal to zero
            Dr(:,ivtPeriodic) = Dr(:,ivt)+Dr(:,ivtPeriodic)
            Dr(:,ivt)         = 0.0_DP

          case default
            ! Get vertex number
            ivt = IverticesAtBoundary(ivbd)

            ! Initialise variable values [x,y,0,time] for parser
            DvariableValues(NDIM1D)   = DvertexCoords(1,ivt)
            DvariableValues(NDIM2D)   = DvertexCoords(2,ivt)
            DvariableValues(NDIM3D)   = 0.0_DP
            DvariableValues(NDIM3D+1) = ttime

            ! Check for user-defined callback-function
            if (present(fcb_calcBoundaryvalues)) then
              
              ! Predict boundary values
              ipos = Kdiagonal(ivt)
              do ivar = 1, nvar
                Du(ivar,ivt) = Du(ivar,ivt)+Dr(ivar,ivt)/Da(ivar,ivar,ipos)
              end do

              ! Get desired boundary values from parser
              do iexpr = 1, nexpr
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+iexpr,&
                    DvariableValues, DbdrValues(iexpr))
              end do
              
              ! Get vertex number of predecessor
              if (ivbd .eq. ivbdFirst) then
                ivtL = IverticesAtBoundary(ivbdLast)
              else
                ivtL = IverticesAtBoundary(ivbd-1)
              end if
              
              ! Get vertex number of sucessor
              if (ivbd .eq. ivbdLast) then
                ivtR = IverticesAtBoundary(ivbdFirst)
              else
                ivtR = IverticesAtBoundary(ivbd+1)
              end if

              ! Compute the analytical and approximate normal vectors
              if (DvertexParameterValue(ivbd) .eq. dminValue) then
                ! We are at the beginning of a segment that includes
                ! the starting point.  Otherwise, the segment counter
                ! would not have been increased.
                
                ! Thus, compute the analytical normal vector from the right
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_RIGHT)
                
                ! Compute the approximate normal vector from the right element
                dnx = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dny = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              elseif (DvertexParameterValue(ivbd) .eq. dmaxValue) then
                ! We are at the end of a segment that includes the end
                ! point.  Otherwise, the segment counter would have
                ! been increased already.
                
                ! Thus, compute the analytical normal vector from the left.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2), BDR_NORMAL_LEFT)
                
                ! Compute the approximate normal vector from the left element
                dnx = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dny = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
                
              else
                ! We are "inside" a segment. Thus, compute the mean
                ! value of left and right normal vector.
                call boundary_getNormalVec2D(rboundary, ibct,&
                    DvertexParameterValue(ivbd), DbdrNormal(1),&
                    DbdrNormal(2))
                
                ! Compute the approximate normal vector
                dnxL = DvertexCoords(2,ivt) -DvertexCoords(2,ivtL)
                dnyL = DvertexCoords(1,ivtL)-DvertexCoords(1,ivt)
                
                dnxR = DvertexCoords(2,ivtR)-DvertexCoords(2,ivt)
                dnyR = DvertexCoords(1,ivt) -DvertexCoords(1,ivtR)
                
                dwL = sqrt( dnxL*dnxL + dnyL*dnyL )
                dwR = sqrt( dnxR*dnxR + dnyR*dnyR )
                
                dnx = (dnxL/dwL + dnxR/dwR)/(1/dwL + 1/dwR)
                dny = (dnyL/dwL + dnyR/dwR)/(1/dwL + 1/dwR)
                dw  = sqrt( dnx*dnx + dny*dny )
                
                DpointNormal(1) = dnx/dw
                DpointNormal(2) = dny/dw
              end if
              
              ! Apply callback function to determine boundary conditions
              call fcb_calcBoundaryvalues(DbdrNormal, DpointNormal, DbdrValues,&
                  IbdrCondType(isegment), Du(:,ivt), Du0(:,ivt), istatus)
              
            else
              
              ! Impose prescribed Dirichlet boundary conditions
              do ivar = 1, nvar
                call fparser_evalFunction(rfparser, nexpr*(isegment-1)+ivar,&
                    DvariableValues, Du(ivar,ivt))
              end do
              
            end if
            
            ! Nullify residual entry
            Dr(:,ivt) = 0.0_DP
            
          end select
        end do
      end do
    end subroutine filtersolution_Mat79Intl1_2D

  end subroutine bdrf_filterSolutionScalar

end module boundaryfilter
