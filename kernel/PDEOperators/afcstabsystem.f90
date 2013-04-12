!##############################################################################
!# ****************************************************************************
!# <name> afcstabsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the algebraic
!# flux correction methodology proposed by Kuzmin, Moeller and Turek
!# in a series of publications. As a starting point for systems of
!# conservation laws, the reader is referred to the book chapter
!#
!#     D. Kuzmin and M. Moeller, Algebraic flux correction
!#     II. Compressible Euler Equations, In: D. Kuzmin et al. (eds),
!#     Flux-Corrected Transport: Principles, Algorithms, and
!#     Applications, Springer, 2005, 207-250.
!#
!# A more detailed description of the algorithms is given in the
!# comments of the subroutine implementing the corresponding
!# discretisation schemes. All methods are based on the stabilisation
!# structure t_afcstab which is defined in the underlying module
!# afcstabbase. The initialisation as a system stabilisation
!# structure is done by the routine afcsys_initStabilisation.
!#
!# The following routines are available:
!#
!# 1.) afcsys_initStabilisation = afcsys_initStabByMatrixSc /
!#                                afcsys_initStabByMatrixBl /
!#                                afcsys_initStabByGroupFEMSetSc /
!#                                afcsys_initStabByGroupFEMSetBl
!#     -> Initialises the stabilisation structure
!#
!# 2.) afcsys_initPerfConfig
!#       -> Initialises the global performance configuration
!#
!# </purpose>
!##############################################################################

module afcstabsystem

!$use omp_lib
  use afcstabbase
!  use basicgeometry
!  use collection
  use fsystem
  use genoutput
  use groupfembase
!  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use perfconfig
  use spatialdiscretisation
  use storage

  implicit none

  private

  public :: afcsys_perfconfig
  public :: afcsys_initStabilisation
  public :: afcsys_initPerfConfig

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of nodes to handle simultaneously when building matrices
#ifndef AFCSYS_NEQSIM
  integer, parameter, public :: AFCSYS_NEQSIM = 128
#endif

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of edges to handle simultaneously when building matrices
#ifndef AFCSYS_NEDGESIM
  integer, parameter, public :: AFCSYS_NEDGESIM = 64
#endif
!</constantblock>

!</constants>

  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: afcsys_perfconfig

  ! ****************************************************************************

  interface afcsys_initStabilisation
    module procedure afcsys_initStabByMatrixSc
    module procedure afcsys_initStabByMatrixBl
    module procedure afcsys_initStabByGroupFEMSetSc
    module procedure afcsys_initStabByGroupFEMSetBl
  end interface

contains

  !****************************************************************************

!<subroutine>

  subroutine afcsys_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<input>
  ! OPTIONAL: performance configuration that should be used to initialise
  ! the global performance configuration. If not present, the values of
  ! the legacy constants is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
!</subroutine>

    if (present(rperfconfig)) then
      afcsys_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(afcsys_perfconfig)
      afcsys_perfconfig%NEQSIM   = AFCSYS_NEQSIM
      afcsys_perfconfig%NEDGESIM = AFCSYS_NEDGESIM
    end if

  end subroutine afcsys_initPerfConfig

  !*****************************************************************************

!<subroutine>

  subroutine afcsys_initStabByMatrixBl(rmatrix, rafcstab,&
      rblockDiscretisation, NVARtransformed)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a system stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
    !
    ! Note that the matrix is required as block matrix. If this matrix
    ! contains only one block, then the scalar counterpart of this
    ! subroutine is called with the corresponding scalar submatrix.
!</description>

!<input>
    ! template block matrix
    type(t_matrixBlock), intent(in) :: rmatrix

    ! OPTIONAL: block discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_blockDiscretisation), intent(in), optional :: rblockDiscretisation

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVAR is taken from the template matrix
    integer, intent(in), optional :: NVARtransformed
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>


    ! Check if block matrix has only one block
    if ((rmatrix%nblocksPerCol .eq. 1) .and.&
        (rmatrix%nblocksPerRow .eq. 1)) then
      if (present(rblockDiscretisation)) then
        call afcsys_initStabByMatrixSc(rmatrix%RmatrixBlock(1,1), rafcstab,&
            rblockDiscretisation%RspatialDiscr(1), NVARtransformed)
      else
        call afcsys_initStabByMatrixSc(rmatrix%RmatrixBlock(1,1), rafcstab,&
            NVARtransformed=NVARtransformed)
      end if
      return
    end if

    ! Check that number of columns equations number of rows
    if (rmatrix%nblocksPerCol .ne.&
        rmatrix%nblocksPerRow) then
      call output_line('Block matrix must have equal number of columns and rows!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_initStabByMatrixBl')
      call sys_halt()
    end if

    ! Check if matrix exhibits group structure
    if (rmatrix%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_initStabByMatrixBl')
      call sys_halt()
    end if

    ! Set atomic data from first block
    rafcstab%NVAR    = rmatrix%nblocksPerCol
    rafcstab%NEQ     = rmatrix%RmatrixBlock(1,1)%NEQ
    rafcstab%NEDGE   = (rmatrix%RmatrixBlock(1,1)%NA-&
                        rmatrix%RmatrixBlock(1,1)%NEQ)/2
    rafcstab%NNVEDGE = 0

    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rafcstab%NVAR
    end if

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

    ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
    call afcstab_allocEdgeStructure(rafcstab,6)
    call afcstab_genEdgeList(rmatrix%RmatrixBlock(1,1), rafcstab)

    ! Allocate internal data structures and
    call afcstab_allocInternalData(rafcstab, .false.,&
        rblockDiscretisation=rblockDiscretisation)

  end subroutine afcsys_initStabByMatrixBl

  ! ****************************************************************************

!<subroutine>
  subroutine afcsys_initStabByMatrixSc(rmatrix, rafcstab,&
      rdiscretisation, NVARtransformed)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a system stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
    !
    ! Note that the matrix is required as scalar matrix. It can be
    ! stored in interleave format.
!</description>

!<input>
    ! template matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: spatial discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVAR is taken from the template matrix
    integer, intent(in), optional :: NVARtransformed
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>


    ! Set atomic data
    rafcstab%NVAR    = rmatrix%NVAR
    rafcstab%NEQ     = rmatrix%NEQ
    rafcstab%NEDGE   = (rmatrix%NA-rmatrix%NEQ)/2
    rafcstab%NNVEDGE = 0

    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rmatrix%NVAR
    end if

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

    ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
    call afcstab_allocEdgeStructure(rafcstab,6)
    call afcstab_genEdgeList(rmatrix, rafcstab)

    ! Allocate internal data structures
    call afcstab_allocInternalData(rafcstab, .false.,&
        rdiscretisation=rdiscretisation)

  end subroutine afcsys_initStabByMatrixSc

  !*****************************************************************************

!<subroutine>

  subroutine afcsys_initStabByGroupFEMSetBl(rgroupFEMSet, rafcstab,&
      rblockDiscretisation, NVARtransformed)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a system stabilisation. The group finite element set
    ! is used to determine the number of equations and the number of
    ! edges. Common data structures are shared with the group FE set.
!</description>

!<input>
    ! The group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Block discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_blockDiscretisation), intent(in) :: rblockDiscretisation

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVAR is taken from the template matrix
    integer, intent(in), optional :: NVARtransformed
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>


    ! Check if block discretisation has only one block
    if (rblockDiscretisation%ncomponents .eq. 1) then
      call afcsys_initStabByGroupFEMSetSc(rgroupFEMSet, rafcstab,&
          rblockDiscretisation%RspatialDiscr(1), NVARtransformed)
      return
    end if

    ! Set atomic data
    rafcstab%NVAR    = rgroupFEMSet%NVAR
    rafcstab%NEQ     = rgroupFEMSet%NEQ
    rafcstab%NEDGE   = rgroupFEMSet%NEDGE
    rafcstab%NNVEDGE = 0

    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rafcstab%NVAR
    end if

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

    ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
    if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
      rafcstab%h_IedgeListIdx     = rgroupFEMSet%h_IedgeListIdx
      rafcstab%h_IedgeList        = rgroupFEMSet%h_IedgeList
      rafcstab%iduplicationFlag   = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
      rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
          iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
    else
      call output_line('Group finite element set does not provide edge structure',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_initStabByGroupFEMSetBl')
      call sys_halt()
    end if

    ! Allocate internal data structures
    call afcstab_allocInternalData(rafcstab, .false.,&
        rblockDiscretisation=rblockDiscretisation)

  end subroutine afcsys_initStabByGroupFEMSetBl

!*****************************************************************************

!<subroutine>

  subroutine afcsys_initStabByGroupFEMSetSc(rgroupFEMSet, rafcstab,&
      rdiscretisation, NVARtransformed)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a system stabilisation. The group finite element set
    ! is used to determine the number of equations and the number of
    ! edges. Common data structures are shared with the group FE set.
!</description>

!<input>
    ! The group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! OPTIONAL: spatial discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVAR is taken from the template matrix
    integer, intent(in), optional :: NVARtransformed
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>


    ! Set atomic data
    rafcstab%NVAR    = rgroupFEMSet%NVAR
    rafcstab%NEQ     = rgroupFEMSet%NEQ
    rafcstab%NEDGE   = rgroupFEMSet%NEDGE
    rafcstab%NNVEDGE = 0

    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rafcstab%NVAR
    end if

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

    ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
    if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
      rafcstab%h_IedgeListIdx   = rgroupFEMSet%h_IedgeListIdx
      rafcstab%h_IedgeList      = rgroupFEMSet%h_IedgeList
      rafcstab%iduplicationFlag = ior(rafcstab%iduplicationFlag,&
                                      AFCSTAB_SHARE_EDGELIST)
      rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
          iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
    else
      call output_line('Group finite element set does not provide edge structure',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_initStabByGroupFEMSetSc')
      call sys_halt()
    end if

    ! Allocate internal data structures
    call afcstab_allocInternalData(rafcstab, .false.,&
        rdiscretisation=rdiscretisation)

  end subroutine afcsys_initStabByGroupFEMSetSc

end module afcstabsystem
