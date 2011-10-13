!##############################################################################
!# ****************************************************************************
!# <name> groupfemsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines for applying the
!# group-finite element formulation to systems of conservation laws.
!# The technique was proposed by C.A.J. Fletcher in:
!#
!#     C.A.J. Fletcher, The group finite element formulation
!#     Computer Methods in Applied Mechanics and Engineering (ISSN
!#     0045-7825), vol. 37, April 1983, p. 225-244.
!#
!# The group finite element formulation uses the same basis functions
!# for the unknown solution and the fluxes. This allows for an
!# efficient matrix assemble, whereby the constant coefficient
!# matrices can be assembled once and for all at the beginning of the
!# simulation and each time the grid is modified.
!#
!# The following routines are available:
!#
!# 1.) gfsys_buildOperator = gfsys_buildOperatorNodeScalar /
!#                           gfsys_buildOperatorNodeBlock1 /
!#                           gfsys_buildOperatorNodeBlock2 /
!#                           gfsys_buildOperatorNodeBlock3 /
!#                           gfsys_buildOperatorEdgeScalar /
!#                           gfsys_buildOperatorEdgeBlock1 /
!#                           gfsys_buildOperatorEdgeBlock2 /
!#                           gfsys_buildOperatorEdgeBlock3
!#     -> Assembles a discrete operator node-by-ndoe or
!#        edge-by-edge by the group finite element formulation
!#
!# 2.) gfsys_buildVectorNode = gfsys_buildVectorNodeScalar /
!#                             gfsys_buildVectorNodeBlock
!#     -> Assembles a discrete vector node-by-node by the
!#        group finite element formulation
!#
!# 3.) gfsys_buildVectorEdge = gfsys_buildVectorEdgeScalar /
!#                             gfsys_buildVectorEdgeBlock
!#     -> Assembles a discrete vector edge-by-edge by the
!#        group finite element formulation
!#
!# 4.) gfsys_initPerfConfig
!#      -> Initialises the global performance configuration
!#
!# </purpose>
!##############################################################################

module groupfemsystem

  use afcstabbase
  use basicgeometry
  use collection
  use fsystem
  use genoutput
  use groupfembase
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use perfconfig
  use spatialdiscretisation
  use storage
  
  implicit none

  private

  public :: gfsys_initPerfConfig
  public :: gfsys_buildOperator
  public :: gfsys_buildVectorNode
  public :: gfsys_buildVectorEdge

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of nodes to handle simultaneously when building matrices
#ifndef GFSYS_NEQSIM
  integer, parameter, public :: GFSYS_NEQSIM = 128
#endif

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of edges to handle simultaneously when building matrices
#ifndef GFSYS_NEDGESIM
  integer, parameter, public :: GFSYS_NEDGESIM = 64
#endif

  ! Number of nonzero entries to handle simultaneously when building matrices
#ifndef GFSYS_NASIM
  integer, parameter, public :: GFSYS_NASIM = 1000
#endif
!</constantblock>

!</constants>

  !************************************************************************
  
  ! global performance configuration
  type(t_perfconfig), target, save :: gfsys_perfconfig

  ! ****************************************************************************

  interface gfsys_buildOperator
    module procedure gfsys_buildOperatorNodeBlock1
    module procedure gfsys_buildOperatorNodeBlock2
    module procedure gfsys_buildOperatorNodeBlock3
    module procedure gfsys_buildOperatorNodeScalar
    module procedure gfsys_buildOperatorEdgeScalar
    module procedure gfsys_buildOperatorEdgeBlock1
    module procedure gfsys_buildOperatorEdgeBlock2
    module procedure gfsys_buildOperatorEdgeBlock3
  end interface
  
  interface gfsys_buildVectorNode
    module procedure gfsys_buildVectorNodeScalar
    module procedure gfsys_buildVectorNodeBlock
  end interface

  interface gfsys_buildVectorEdge
    module procedure gfsys_buildVectorEdgeScalar
    module procedure gfsys_buildVectorEdgeBlock
  end interface

contains

  !****************************************************************************

!<subroutine>

  subroutine gfsys_initPerfConfig(rperfconfig)

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
      gfsys_perfconfig = rperfconfig
    else
      gfsys_perfconfig%NEQSIM = GFSYS_NEQSIM
      gfsys_perfconfig%NEDGESIM = GFSYS_NEDGESIM
      gfsys_perfconfig%NASIM = GFSYS_NASIM
    end if
  
  end subroutine gfsys_initPerfConfig

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildOperatorNodeBlock1(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSys_sim, dscale, bclear, rmatrix,&
      rcollection, rafcstab, fcb_calcOperatorNodeSys, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! Note that this routine is designed for block matrices/vectors.
    ! If there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute matrix entries
    include 'intf_calcMatrixDiagSys_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorNodeSys.inc'
    optional :: fcb_calcOperatorNodeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_array), dimension(:,:), allocatable  :: rarray
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    real(DP), dimension(:), pointer :: p_Dx
    integer, dimension(:,:), pointer :: p_InodeList2D
    integer, dimension(:), pointer :: p_InodeListIdx,p_InodeList1D
    logical :: bisFullMatrix

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if user-defined assembly is provided
    if (present(fcb_calcOperatorNodeSys)) then
      
      call fcb_calcOperatorNodeSys(rgroupFEMSet, rx, rmatrix, dscale,&
          bclear, fcb_calcMatrixDiagSys_sim, rcollection, rafcstab)
      ! That`s it
      return

    elseif ((rx%nblocks            .eq. 1) .and.&
            (rmatrix%nblocksPerCol .eq. 1) .and.&
            (rmatrix%nblocksPerRow .eq. 1)) then

      ! Call scalar version of this routine
      call gfsys_buildOperatorNodeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSys_sim, dscale, bclear, rmatrix%RmatrixBlock(1,1),&
          rcollection, rafcstab, rperfconfig=rperfconfig)
      ! That`s it
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => gfsys_perfconfig
    end if

    ! Check if block matrix exhibits group structure
    if (rmatrix%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeBlock1')
      call sys_halt()
    end if

    ! Check if matrix and vector have the same data type
    if (rx%cdataType .ne. rgroupFEMSet%cdataType) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeBlock1')
      call sys_halt()
    end if

    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)
      
    case (GFEM_NODEBASED)
      !-------------------------------------------------------------------------
      ! Node-based assembly
      !-------------------------------------------------------------------------
      
      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODELIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA) .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorNodeBlock1')
        call sys_halt()
      end if

      ! Allocate temporal memory
      allocate(rarray(rx%nblocks,rx%nblocks))

      ! What data types are we?
      select case(rgroupFEMSet%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
        call gfem_getbase_array(rmatrix, rarray, bisFullMatrix)
        call lsysbl_getbase_double(rx, p_Dx)
    
        ! Check if only a subset of the matrix is required
        if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
          ! Set pointers
          call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
          call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList1D)
          
          ! Assemble operator node-by-node
          call doOperatorDble(rx%RvectorBlock(1)%NEQ, rx%nblocks,&
              p_InodeListIdx, p_InodeList1D, p_DcoeffsAtNode, p_Dx, dscale,&
              bclear, bisFullMatrix, rarray)
        else
          ! Set pointers
          call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
          call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList2D)
          
          ! Assemble selected part of the operator node-by-node
          call doOperatorDbleSel(rx%RvectorBlock(1)%NEQ, rx%nblocks,&
              p_InodeListIdx, p_InodeList2D, p_DcoeffsAtNode, p_Dx, dscale,&
              bclear, bisFullMatrix, rarray)
        end if
        
        ! Deallocate temporal memory
        deallocate(rarray)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorNodeBlock1')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorNodeBlock1')
      call sys_halt()
    end select

  contains

    !**************************************************************
    ! Assemble operator node-by-node without stabilisation

    subroutine doOperatorDble(NEQ, NVAR, InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, bisFullMatrix, rarray)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: InodeListIdx,InodeList
      integer, intent(in) :: NEQ,NVAR
      logical, intent(in) :: bisFullMatrix

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: rarray
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: Dcoefficients
      integer, dimension(:,:), pointer :: IdofsAtNode
      
      ! local variables
      integer :: idx,IAset,IAmax
      integer :: ij,j,ivar,jvar,ijpos
      
      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,IdofsAtNode,&
      !$omp         IAmax,idx,ij,j,ivar,jvar,ijpos)&
      !$omp if (size(InodeList) > p_rperfconfig%NAMIN_OMP)
      
      ! allocate temporal memory
      allocate(DdataAtNode(NVAR,p_rperfconfig%NASIM))
      if (bisFullMatrix) then
        allocate(Dcoefficients(NVAR*NVAR,1,p_rperfconfig%NASIM))
      else
        allocate(Dcoefficients(NVAR,1,p_rperfconfig%NASIM))
      end if
      allocate(IdofsAtNode(2,p_rperfconfig%NASIM))
      
      ! Loop over all nonzero matrix entries
      !$omp do schedule(static,1)
      do IAset = 1, size(InodeList), p_rperfconfig%NASIM

        ! We always handle NASIM matrix entries simultaneously.
        ! How many matrix entries have we actually here?
        ! Get the maximum position of matrix entries, such that we handle 
        ! at most NASIM matrix entries simultaneously.
        
        IAmax = min(size(InodeList), IAset-1+p_rperfconfig%NASIM)

        ! Loop through all nonzero matrix entries in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IAmax-IAset+1
          
          ! Get position of matrix entry
          ij = idx+IAset-1
          
          ! Get column number of matrix entry
          j = InodeList(ij)
          
          ! Fill auxiliary array
          IdofsAtNode(1,idx) = j
          IdofsAtNode(2,idx) = ij
          DdataAtNode(:,idx) = Dx(j,:)
        end do
        
        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IAmax-IAset+1),&
            DcoeffsAtNode(:,IAset:IAmax),&
            IdofsAtNode(:,1:IAmax-IAset+1),&
            dscale, IAmax-IAset+1,&
            Dcoefficients(:,:,1:IAmax-IAset+1), rcollection)
        
        ! Loop through all nonzero matrix entries in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then

          if (bisFullMatrix) then
            
            do idx = 1, IAmax-IAset+1
              
              ! Get position of matrix entry
              ij = idx+IAset-1
              
              ! Update the global operator
              do ivar = 1, NVAR
                do jvar = 1, NVAR
                  ijpos = NVAR*(ivar-1)+jvar
                  rarray(jvar,ivar)%p_Ddata(ij) = Dcoefficients(ijpos,1,idx)
                end do
              end do
            end do
          
          else   ! matrix is block-diagonal

            do idx = 1, IAmax-IAset+1
              
              ! Get position of matrix entry
              ij = idx+IAset-1
              
              ! Update the global operator
              do ivar = 1, NVAR
                rarray(ivar,ivar)%p_Ddata(ij) = Dcoefficients(ivar,1,idx)
              end do
            end do
            
          end if

        else   ! do not clear matrix
          
          if (bisFullMatrix) then

            do idx = 1, IAmax-IAset+1
              
              ! Get position of matrix entry
              ij = idx+IAset-1
              
              ! Update the global operator
              do ivar = 1, NVAR
                do jvar = 1, NVAR
                  ijpos = NVAR*(ivar-1)+jvar
                  rarray(jvar,ivar)%p_Ddata(ij) = rarray(jvar,ivar)%p_Ddata(ij)&
                                                + Dcoefficients(ijpos,1,idx)
                end do
              end do
            end do

          else   ! matrix is block-diagonal

            do idx = 1, IAmax-IAset+1
              
              ! Get position of matrix entry
              ij = idx+IAset-1
              
              ! Update the global operator
              do ivar = 1, NVAR
                rarray(ivar,ivar)%p_Ddata(ij) = rarray(ivar,ivar)%p_Ddata(ij)&
                                              +  Dcoefficients(ivar,1,idx)
              end do
            end do

          end if
        end if
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(IdofsAtNode)
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorDble

    !**************************************************************
    ! Assemble operator node-by-node without stabilisation

    subroutine doOperatorDbleSel(NEQ, NVAR, InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, bisFullMatrix, rarray)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: InodeList
      integer, dimension(:), intent(in) :: InodeListIdx
      integer, intent(in) :: NEQ,NVAR
      logical, intent(in) :: bisFullMatrix

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: rarray
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: Dcoefficients
      
      ! local variables
      integer :: idx,IAset,IAmax
      integer :: ij,j,ivar,jvar,ijpos

      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,IAmax,idx,ij,j,ivar,jvar,ijpos)&
      !$omp if (size(InodeList,2) > p_rperfconfig%NAMIN_OMP)
      
      ! allocate temporal memory
      allocate(DdataAtNode(NVAR,p_rperfconfig%NASIM))
      if (bisFullMatrix) then
        allocate(Dcoefficients(NVAR*NVAR,1,p_rperfconfig%NASIM))
      else
        allocate(Dcoefficients(NVAR*NVAR,1,p_rperfconfig%NASIM))
      end if
      
      ! Loop over all nonzero matrix entries
      !$omp do schedule(static,1)
      do IAset = 1, size(InodeList,2), p_rperfconfig%NASIM
        
        ! We always handle NASIM matrix entries simultaneously.
        ! How many matrix entries have we actually here?
        ! Get the maximum position of matrix entries, such that we handle 
        ! at most NASIM matrix entries simultaneously.
        
        IAmax = min(size(InodeList,2), IAset-1+p_rperfconfig%NASIM)
        
        ! Loop through all nonzero matrix entries in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IAmax-IAset+1
          
          ! Fill auxiliary array
          DdataAtNode(:,idx) = Dx(InodeList(1,idx+IAset-1),:)
        end do
        
        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IAmax-IAset+1),&
            DcoeffsAtNode(:,IAset:IAmax),&
            InodeList(:,IAset:IAmax),&
            dscale, IAmax-IAset+1,&
            Dcoefficients(:,:,1:IAmax-IAset+1), rcollection)
        
        ! Loop through all nonzero matrix entries in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then

          if (bisFullMatrix) then
            
            do idx = 1, IAmax-IAset+1
              
              ! Get position of matrix entry
              ij = InodeList(2,idx+IAset-1)
              
              ! Update the global operator
              do ivar = 1, NVAR
                do jvar = 1, NVAR
                  ijpos = NVAR*(ivar-1)+jvar
                  rarray(jvar,ivar)%p_Ddata(ij) = Dcoefficients(ijpos,1,idx)
                end do
              end do
            end do
          
          else   ! matrix is block-diagonal

            do idx = 1, IAmax-IAset+1
              
              ! Get position of matrix entry
              ij = InodeList(2,idx+IAset-1)
              
              ! Update the global operator
              do ivar = 1, NVAR
                rarray(ivar,ivar)%p_Ddata(ij) = Dcoefficients(ivar,1,idx)
              end do
            end do
            
          end if

        else   ! do not clear matrix
          
          if (bisFullMatrix) then

            do idx = 1, IAmax-IAset+1
              
              ! Get position of matrix entry
              ij = InodeList(2,idx+IAset-1)
              
              ! Update the global operator
              do ivar = 1, NVAR
                do jvar = 1, NVAR
                  ijpos = NVAR*(ivar-1)+jvar
                  rarray(jvar,ivar)%p_Ddata(ij) = rarray(jvar,ivar)%p_Ddata(ij)&
                                                + Dcoefficients(ijpos,1,idx)
                end do
              end do
            end do

          else   ! matrix is block-diagonal

            do idx = 1, IAmax-IAset+1
              
              ! Get position of matrix entry
              ij = InodeList(2,idx+IAset-1)
              
              ! Update the global operator
              do ivar = 1, NVAR
                rarray(ivar,ivar)%p_Ddata(ij) = rarray(ivar,ivar)%p_Ddata(ij)&
                                              +  Dcoefficients(ivar,1,idx)
              end do
            end do

          end if
        end if
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorDbleSel
    
  end subroutine gfsys_buildOperatorNodeBlock1

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildOperatorNodeBlock2(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSys_sim, dscale, bclear, rmatrix,&
      rcollection, rafcstab, fcb_calcOperatorNodeSys, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! Note that this routine is designed for block matrices/vectors.
    ! If there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! Vector on which the matrix entries may depend.
    type(t_vectorScalar), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute matrix entries
    include 'intf_calcMatrixDiagSys_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorNodeSys.inc'
    optional :: fcb_calcOperatorNodeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rxBlock

    ! Check if user-defined assembly is provided
    if (present(fcb_calcOperatorNodeSys)) then
      
      ! Create auxiliary 1-block vector
      call lsysbl_createVecFromScalar(rx, rxBlock)

      call fcb_calcOperatorNodeSys(rgroupFEMSet, rxBlock, rmatrix, dscale,&
          bclear, fcb_calcMatrixDiagSys_sim, rcollection, rafcstab)

      ! Release auxiliary 1-block vector
      call lsysbl_releaseVector(rxBlock)

    elseif ((rmatrix%nblocksPerCol .eq. 1) .and.&
            (rmatrix%nblocksPerRow .eq. 1)) then

      ! Call scalar version of this routine
      call gfsys_buildOperatorNodeScalar(rgroupFEMSet, rx,&
          fcb_calcMatrixDiagSys_sim, dscale, bclear,&
          rmatrix%RmatrixBlock(1,1), rcollection, rafcstab,&
          rperfconfig=rperfconfig)

    else
      call output_line('Matrix must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorNodeBlock2')
      call sys_halt()
    end if
    
  end subroutine gfsys_buildOperatorNodeBlock2

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildOperatorNodeBlock3(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSys_sim, dscale, bclear, rmatrix,&
      rcollection, rafcstab, fcb_calcOperatorNodeSys, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! Note that this routine is designed for block matrices/vectors.
    ! If there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute matrix entries
    include 'intf_calcMatrixDiagSys_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorNodeSys.inc'
    optional :: fcb_calcOperatorNodeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixBlock) :: rmatrixBlock

    ! Check if user-defined assembly is provided
    if (present(fcb_calcOperatorNodeSys)) then
      
      ! Create auxiliary 1-block matrix
      call lsysbl_createMatFromScalar(rmatrix, rmatrixBlock)

      call fcb_calcOperatorNodeSys(rgroupFEMSet, rx, rmatrixBlock, dscale,&
          bclear, fcb_calcMatrixDiagSys_sim, rcollection, rafcstab)

      ! Release auxiliary 1-block matrix
      call lsysbl_releaseMatrix(rmatrixBlock)

    elseif (rx%nblocks .eq. 1) then

      ! Call scalar version of this routine
      call gfsys_buildOperatorNodeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSys_sim, dscale, bclear, rmatrix, rcollection,&
          rafcstab, rperfconfig=rperfconfig)
      
    else
      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorNodeBlock3')
      call sys_halt()
    end if
    
  end subroutine gfsys_buildOperatorNodeBlock3

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildOperatorNodeScalar(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSys_sim, dscale, bclear, rmatrix,&
      rcollection, rafcstab, fcb_calcOperatorNodeSys, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! Vector on which the matrix entries may depend.
    type(t_vectorScalar), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute matrix entries
    include 'intf_calcMatrixDiagSys_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorNodeSys.inc'
    optional :: fcb_calcOperatorNodeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rxBlock
    type(t_matrixBlock) :: rmatrixBlock
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    real(DP), dimension(:), pointer :: p_Ddata,p_Dx
    integer, dimension(:,:), pointer :: p_InodeList2D
    integer, dimension(:), pointer :: p_InodeListIdx,p_InodeList1D

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if user-defined assembly is provided
    if (present(fcb_calcOperatorNodeSys)) then
      ! Create auxiliary 1-block vector and matrix
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createMatFromScalar(rmatrix, rmatrixBlock)
      
      ! Call user-defined assembly
      call fcb_calcOperatorNodeSys(rgroupFEMSet, rxBlock, rmatrixBlock, dscale,&
          bclear, fcb_calcMatrixDiagSys_sim, rcollection, rafcstab)
      
      ! Release auxiliary 1-block vector and matrix
      call lsysbl_releaseVector(rxBlock)
      call lsysbl_releaseMatrix(rmatrixBlock)

      ! That`s it
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => gfsys_perfconfig
    end if

    ! Check if matrix and vector have the same data type
    if ((rmatrix%cdataType .ne. rx%cdataType) .or.&
        (rmatrix%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorNodeScalar')
      call sys_halt()
    end if
    
    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_NODEBASED)
      !-------------------------------------------------------------------------
      ! Node-based assembly
      !-------------------------------------------------------------------------

      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODELIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA) .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorNodeScalar')
        call sys_halt()
      end if

      ! What data types are we?
      select case(rgroupFEMSet%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
        call lsyssc_getbase_double(rmatrix, p_Ddata)
        call lsyssc_getbase_double(rx, p_Dx)

        ! What type of matrix are we?
        select case(rmatrix%cinterleavematrixFormat)
          
        case (LSYSSC_MATRIX1)
          
          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
            ! Set pointers
            call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList1D)
            
            ! Assemble operator node-by-node
            call doOperatorDble(rx%NEQ, rmatrix%NA, rx%NVAR, rx%NVAR*rx%NVAR,&
                p_InodeListIdx, p_InodeList1D, p_DcoeffsAtNode, p_Dx, dscale,&
                bclear, p_Ddata)
          else
            ! Set pointers
            call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList2D)
            
            ! Assemble selected part of the operator node-by-node
            call doOperatorDbleSel(rx%NEQ, rmatrix%NA, rx%NVAR, rx%NVAR*rx%NVAR,&
                p_InodeListIdx, p_InodeList2D, p_DcoeffsAtNode, p_Dx, dscale,&
                bclear, p_Ddata)
          end if

        case (LSYSSC_MATRIXD)
          
          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
            ! Set pointers
            call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList1D)
            
            ! Assemble operator node-by-node
            call doOperatorDble(rx%NEQ, rmatrix%NA, rx%NVAR, rx%NVAR,&
                p_InodeListIdx, p_InodeList1D, p_DcoeffsAtNode, p_Dx, dscale,&
                bclear, p_Ddata)
          else
            ! Set pointers
            call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList2D)
            
            ! Assemble selected part of the operator node-by-node
            call doOperatorDbleSel(rx%NEQ, rmatrix%NA, rx%NVAR, rx%NVAR,&
                p_InodeListIdx, p_InodeList2D, p_DcoeffsAtNode, p_Dx, dscale,&
                bclear, p_Ddata)
          end if

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorNodeScalar')
          call sys_halt()
        end select
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorNodeScalar')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorNodeScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow

    !**************************************************************
    ! Assemble operator node-by-node without stabilisation

    subroutine doOperatorDble(NEQ, NA, NVAR, MVAR, InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: InodeListIdx,InodeList
      integer, intent(in) :: NEQ,NA,NVAR,MVAR

      ! input/output parameters
      real(DP), dimension(MVAR,NA), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: Dcoefficients
      integer, dimension(:,:), pointer :: IdofsAtNode
      
      ! local variables
      integer :: idx,IAset,IAmax
      integer :: ij,j
      
      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,IdofsAtNode,IAmax,idx,ij,j)&
      !$omp if (size(InodeList) > p_rperfconfig%NAMIN_OMP)
      
      ! allocate temporal memory
      allocate(DdataAtNode(NVAR,p_rperfconfig%NASIM))
      allocate(Dcoefficients(MVAR,1,p_rperfconfig%NASIM))
      allocate(IdofsAtNode(2,p_rperfconfig%NASIM))
      
      ! Loop over all nonzero matrix entries
      !$omp do schedule(static,1)
      do IAset = 1, size(InodeList), p_rperfconfig%NASIM

        ! We always handle NASIM matrix entries simultaneously.
        ! How many matrix entries have we actually here?
        ! Get the maximum position of matrix entries, such that we handle 
        ! at most NASIM matrix entries simultaneously.
        
        IAmax = min(size(InodeList), IAset-1+p_rperfconfig%NASIM)

        ! Loop through all nonzero matrix entries in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IAmax-IAset+1
          
          ! Get position of matrix entry
          ij = idx+IAset-1
          
          ! Get column number of matrix entry
          j = InodeList(ij)
          
          ! Fill auxiliary array
          IdofsAtNode(1,idx) = j
          IdofsAtNode(2,idx) = ij
          DdataAtNode(:,idx) = Dx(:,j)
        end do
        
        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IAmax-IAset+1),&
            DcoeffsAtNode(:,IAset:IAmax),&
            IdofsAtNode(:,1:IAmax-IAset+1),&
            dscale, IAmax-IAset+1,&
            Dcoefficients(:,:,1:IAmax-IAset+1), rcollection)
        
        ! Loop through all nonzero matrix entries in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          do idx = 1, IAmax-IAset+1
            
            ! Get position of matrix entry
            ij = idx+IAset-1
            
            ! Update the global operator
            Ddata(:,ij) = Dcoefficients(:,1,idx)
          end do
          
        else   ! do not clear matrix
          
          do idx = 1, IAmax-IAset+1
            
            ! Get position of matrix entry
            ij = idx+IAset-1
            
            ! Update the global operator
            Ddata(:,ij) = Ddata(:,ij) + Dcoefficients(:,1,idx)
          end do
        end if
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(IdofsAtNode)
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorDble

    !**************************************************************
    ! Assemble operator node-by-node without stabilisation

    subroutine doOperatorDbleSel(NEQ, NA, NVAR, MVAR, InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: InodeList
      integer, dimension(:), intent(in) :: InodeListIdx
      integer, intent(in) :: NEQ,NA,NVAR,MVAR

      ! input/output parameters
      real(DP), dimension(MVAR,NA), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: Dcoefficients
      
      ! local variables
      integer :: idx,IAset,IAmax
      integer :: ij,j

      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,IAmax,idx,ij,j)&
      !$omp if (size(InodeList,2) > p_rperfconfig%NAMIN_OMP)
      
      ! allocate temporal memory
      allocate(DdataAtNode(NVAR,p_rperfconfig%NASIM))
      allocate(Dcoefficients(MVAR,1,p_rperfconfig%NASIM))
      
      ! Loop over all nonzero matrix entries
      !$omp do schedule(static,1)
      do IAset = 1, size(InodeList,2), p_rperfconfig%NASIM
        
        ! We always handle NASIM matrix entries simultaneously.
        ! How many matrix entries have we actually here?
        ! Get the maximum position of matrix entries, such that we handle 
        ! at most NASIM matrix entries simultaneously.
        
        IAmax = min(size(InodeList,2), IAset-1+p_rperfconfig%NASIM)
        
        ! Loop through all nonzero matrix entries in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IAmax-IAset+1
          
          ! Fill auxiliary array
          DdataAtNode(:,idx) = Dx(:,InodeList(1,idx+IAset-1))
        end do
        
        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IAmax-IAset+1),&
            DcoeffsAtNode(:,IAset:IAmax),&
            InodeList(:,IAset:IAmax),&
            dscale, IAmax-IAset+1,&
            Dcoefficients(:,:,1:IAmax-IAset+1), rcollection)
        
        ! Loop through all nonzero matrix entries in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          do idx = 1, IAmax-IAset+1
            
            ! Get position of matrix entry
            ij = InodeList(2,idx+IAset-1)
            
            ! Update the global operator
            Ddata(:,ij) = Dcoefficients(:,1,idx)
          end do
          
        else   ! do not clear matrix
          
          do idx = 1, IAmax-IAset+1
            
            ! Get position of matrix entry
            ij = InodeList(2,idx+IAset-1)
            
            ! Update the global operator
            Ddata(:,ij) = Ddata(:,ij) + Dcoefficients(:,1,idx)
          end do
        end if
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorDbleSel

  end subroutine gfsys_buildOperatorNodeScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildOperatorEdgeBlock1(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      dscale, bclear, rmatrix, rcollection, rafcstab,&
      fcb_calcOperatorEdgeSys, rperfconfig)

!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only edge-by-edge assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! stabilisation of AFC-type for hyperbolic problems is performed.
    ! This technique is for instance described in the reference:
    !
    !     D. Kuzmin and M. Moeller, Algebraic flux correction II. Compressible
    !     Euler equations, In: D. Kuzmin et al. (eds), Flux-Corrected
    !     Transport: Principles,  Algorithms, and Applications,
    !     Springer, 2005, 207-250.
    !
    ! Note that this routine is designed for block matrices/vectors.
    ! If there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute local matrices
    include 'intf_calcMatrixDiagSys_sim.inc'
    include 'intf_calcMatrixSys_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorEdgeSys.inc'
    optional :: fcb_calcOperatorEdgeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_array), dimension(:,:), allocatable  :: rarray
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtDiag   
    real(DP), dimension(:), pointer :: p_Dx
    integer, dimension(:,:), pointer :: p_IdiagList,p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    logical :: bisFullMatrix

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if user-defined callback function is present
    if (present(fcb_calcOperatorEdgeSys)) then

      call fcb_calcOperatorEdgeSys(rgroupFEMSet, rx, rmatrix, dscale,&
          bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
          rcollection, rafcstab)
      ! That`s it
      return

    elseif ((rx%nblocks            .eq. 1) .and.&
            (rmatrix%nblocksPerCol .eq. 1) .and.&
            (rmatrix%nblocksPerRow .eq. 1)) then

      ! Call scalar version of this routine
      call gfsys_buildOperatorEdgeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, dscale,&
          bclear, rmatrix%RmatrixBlock(1,1), rcollection, rafcstab,&
          rperfconfig=rperfconfig)
      ! That`s it
      return
    end if
    
    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => gfsys_perfconfig
    end if

    ! Check if block matrix exhibits group structure
    if (rmatrix%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeBlock1')
      call sys_halt()
    end if

    ! Check if matrix and vector have the same data type
    if (rx%cdataType .ne. rgroupFEMSet%cdataType) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeBlock1')
      call sys_halt()
    end if

    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_EDGEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------
      
      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGLIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGDATA) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA) .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeBlock1')
        call sys_halt()
      end if
      
      ! Allocate temporal memory
      allocate(rarray(rx%nblocks,rx%nblocks))
      
      ! Set pointers
      call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
      call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)
      call gfem_getbase_IdiagList(rgroupFEMSet, p_IdiagList)
      
      ! What data types are we?
      select case(rgroupFEMSet%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtDiag(rgroupFEMSet, p_DcoeffsAtDiag)
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        call gfem_getbase_array(rmatrix, rarray, bisFullMatrix)
        call lsysbl_getbase_double(rx, p_Dx)
        
        ! Assemble matrix diagonal
        call doOperatorDiagDble(rx%RvectorBlock(1)%NEQ, rx%nblocks,&
            p_IdiagList, p_DcoeffsAtDiag, p_Dx, dscale, bclear, bisFullMatrix,&
            rarray)
        
        ! Do we have to build the stabilisation?
        if (present(rafcstab)) then
          
          ! Check if stabilisation has been prepared
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
            call output_line('Stabilisation has not been prepared!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeBlock1')
            call sys_halt()
          end if
          
          !-------------------------------------------------------------------
          ! Assemble operator with stabilisation
          !-------------------------------------------------------------------
          call doOperatorStabDble(rx%RvectorBlock(1)%NEQ, rx%nblocks,&
              p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge, p_Dx,&
              dscale, bclear, bisFullMatrix, rarray)
          
        else   ! no stabilisation structure present
          
          !-------------------------------------------------------------------
          ! Assemble operator without stabilisation
          !-------------------------------------------------------------------
          call doOperatorDble(rx%RvectorBlock(1)%NEQ, rx%nblocks,&
              p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge, p_Dx,&
              dscale, bclear, bisFullMatrix, rarray)
        end if
        
        ! Deallocate temporal memory
        deallocate(rarray)
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeBlock1')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeBlock1')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble diagonal part of the operator
    
    subroutine doOperatorDiagDble(NEQ, NVAR, IdiagList, DcoeffsAtDiag,&
        Dx, dscale, bclear, bisFullMatrix, rarray)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtDiag
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear,bisFullMatrix
      integer, dimension(:,:), intent(in) :: IdiagList
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: rarray
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: Dcoefficients

      ! local variables
      integer :: IEQmax,IEQset,i,ia,idx,ijpos,ivar,jvar

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,IEQmax,i,ia,idx,ijpos,ivar,jvar)&
      !$omp if (size(IdiagList,2) > p_rperfconfig%NEQMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtNode(NVAR,p_rperfconfig%NEQSIM))
      if (bisFullMatrix) then
        allocate(Dcoefficients(NVAR*NVAR,1,p_rperfconfig%NEQSIM))
      else
        allocate(Dcoefficients(NVAR,1,p_rperfconfig%NEQSIM))
      end if

      ! Loop over the equations
      !$omp do schedule(static,1)
      do IEQset = 1, size(IdiagList,2), p_rperfconfig%NEQSIM

        ! We always handle NEQSIM equations simultaneously.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most NEQSIM equations simultaneously.
        
        IEQmax = min(size(IdiagList,2), IEQset-1+p_rperfconfig%NEQSIM)
        
        ! Loop through all equations in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IEQmax-IEQset+1

          ! Get actual equation number
          i = IdiagList(1,idx+IEQset-1)
          
          ! Fill auxiliary array
          DdataAtNode(:,idx) = Dx(i,:)
        end do

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DcoeffsAtDiag(:,IEQset:IEQmax),&
            IdiagList(:,IEQset:IEQmax),&
            dscale, IEQmax-IEQset+1,&
            Dcoefficients(:,:,1:IEQmax-IEQset+1), rcollection)
        
        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then

          if (bisFullMatrix) then

            do idx = 1, IEQmax-IEQset+1
              
              ! Get position of diagonal entry
              ia = IdiagList(2,idx+IEQset-1)
              
              ! Update the diagonal coefficient
              do ivar = 1, NVAR
                do jvar = 1, NVAR
                  ijpos = NVAR*(ivar-1)+jvar
                  rarray(jvar,ivar)%p_Ddata(ia) = Dcoefficients(ijpos,1,idx)
                end do
              end do
            end do

          else   ! matrix is block-diagonal

            do idx = 1, IEQmax-IEQset+1
              
              ! Get position of diagonal entry
              ia = IdiagList(2,idx+IEQset-1)
              
              ! Update the diagonal coefficient
              do ivar = 1, NVAR
                rarray(ivar,ivar)%p_Ddata(ia) = Dcoefficients(ivar,1,idx)
              end do
            end do

          end if

        else   ! do not clear matrix

          if (bisFullMatrix) then

            do idx = 1, IEQmax-IEQset+1
              
              ! Get position of diagonal entry
              ia = IdiagList(2,idx+IEQset-1)
              
              ! Update the diagonal coefficient
              do ivar = 1, NVAR
                do jvar = 1, NVAR
                  ijpos = NVAR*(ivar-1)+jvar
                  rarray(jvar,ivar)%p_Ddata(ia) = rarray(jvar,ivar)%p_Ddata(ia)&
                                                + Dcoefficients(ijpos,1,idx)
                end do
              end do
            end do
            
          else   ! matrix is block-diagonal
            
            do idx = 1, IEQmax-IEQset+1
              
              ! Get position of diagonal entry
              ia = IdiagList(2,idx+IEQset-1)
              
              ! Update the diagonal coefficient
              do ivar = 1, NVAR
                rarray(ivar,ivar)%p_Ddata(ia) = rarray(ivar,ivar)%p_Ddata(ia)&
                                              + Dcoefficients(ivar,1,idx)
              end do
            end do
            
          end if
        end if
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel

    end subroutine doOperatorDiagDble

    !**************************************************************
    ! Assemble operator edge-by-edge without stabilisation
    
    subroutine doOperatorDble(NEQ, NVAR, IedgeListIdx, IedgeList,&
        DcoeffsAtEdge, Dx, dscale, bclear, bisFullMatrix, rarray)
      
      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear,bisFullMatrix
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: rarray
      
      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: Dcoefficients

      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ij,ji,ivar,jvar,ijpos
      
      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtEdge,IEDGEmax,idx,iedge,ij,ji,ivar,jvar,ijpos)&
      !$omp if(size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      if (bisFullMatrix) then
        allocate(Dcoefficients(NVAR*NVAR,2,p_rperfconfig%NEDGESIM))
      else
        allocate(Dcoefficients(NVAR,2,p_rperfconfig%NEDGESIM))
      end if

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+p_rperfconfig%NEDGESIM)
        
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dcoefficients(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global matrix
          if (bclear) then

            if (bisFullMatrix) then

              do idx = 1, IEDGEmax-IEDGEset+1
                
                ! Get actual edge number
                iedge = idx+IEDGEset-1
                
                ! Get position of off-diagonal entries
                ij = IedgeList(3,iedge)
                ji = IedgeList(4,iedge)
                
                ! Update the global operator
                do ivar = 1, NVAR
                  do jvar = 1, NVAR
                    ijpos = NVAR*(ivar-1)+jvar
                    rarray(jvar,ivar)%p_Ddata(ij) = Dcoefficients(ijpos,1,idx)
                    rarray(jvar,ivar)%p_Ddata(ji) = Dcoefficients(ijpos,2,idx)
                  end do
                end do
              end do

            else   ! matrix is block-diagonal
              
              do idx = 1, IEDGEmax-IEDGEset+1
                
                ! Get actual edge number
                iedge = idx+IEDGEset-1
                
                ! Get position of off-diagonal entries
                ij = IedgeList(3,iedge)
                ji = IedgeList(4,iedge)
                
                ! Update the global operator
                do ivar = 1, NVAR
                  rarray(ivar,ivar)%p_Ddata(ij) = Dcoefficients(ivar,1,idx)
                  rarray(ivar,ivar)%p_Ddata(ji) = Dcoefficients(ivar,2,idx)
                end do
              end do

            end if

          else   ! do not clear matrix

            if (bisFullMatrix) then

              do idx = 1, IEDGEmax-IEDGEset+1
                
                ! Get actual edge number
                iedge = idx+IEDGEset-1
                
                ! Get position of off-diagonal entries
                ij = IedgeList(3,iedge)
                ji = IedgeList(4,iedge)
                
                ! Update the global operator
                do ivar = 1, NVAR
                  do jvar = 1, NVAR
                    ijpos = NVAR*(ivar-1)+jvar
                    rarray(jvar,ivar)%p_Ddata(ij) = rarray(jvar,ivar)%p_Ddata(ij)&
                                                  + Dcoefficients(ijpos,1,idx)
                    rarray(jvar,ivar)%p_Ddata(ji) = rarray(jvar,ivar)%p_Ddata(ji)&
                                                  + Dcoefficients(ijpos,2,idx)
                  end do
                end do
              end do
              
            else   ! matrix is block-diagonal
              
              do idx = 1, IEDGEmax-IEDGEset+1
                
                ! Get actual edge number
                iedge = idx+IEDGEset-1
                
                ! Get position of off-diagonal entries
                ij = IedgeList(3,iedge)
                ji = IedgeList(4,iedge)
                
                ! Update the global operator
                do ivar = 1, NVAR
                  rarray(ivar,ivar)%p_Ddata(ij) = rarray(ivar,ivar)%p_Ddata(ij)&
                                                + Dcoefficients(ivar,1,idx)
                  rarray(ivar,ivar)%p_Ddata(ji) = rarray(ivar,ivar)%p_Ddata(ji)&
                                                + Dcoefficients(ivar,2,idx)
                end do
              end do
              
            end if
          end if
        end do
        !$omp end do
        
      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorDble

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation
    
    subroutine doOperatorStabDble(NEQ, NVAR, IedgeListIdx, IedgeList,&
        DcoeffsAtEdge, Dx, dscale, bclear, bisFullMatrix, rarray)
      
      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear,bisFullMatrix
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: rarray
      
      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: Dcoefficients

      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ii,ij,ijpos,ivar,ji,jj,jvar
      
      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtEdge,&
      !$omp         IEDGEmax,idx,iedge,ii,ij,ijpos,ivar,ji,jj,jvar)&
      !$omp if(size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      if (bisFullMatrix) then
        allocate(Dcoefficients(NVAR*NVAR,3,p_rperfconfig%NEDGESIM))
      else
        allocate(Dcoefficients(NVAR,3,p_rperfconfig%NEDGESIM))
      end if

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+p_rperfconfig%NEDGESIM)
        
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dcoefficients(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global matrix
          if (bclear) then

            if (bisFullMatrix) then

              do idx = 1, IEDGEmax-IEDGEset+1
                
                ! Get actual edge number
                iedge = idx+IEDGEset-1
                
                ! Get position of off-diagonal entries
                ij = IedgeList(3,iedge)
                ji = IedgeList(4,iedge)
                
                ! Get position of diagonal entries
                ii = IedgeList(5,iedge)
                jj = IedgeList(6,iedge)
                
                ! Update the global operator
                do ivar = 1, NVAR
                  do jvar = 1, NVAR
                    ijpos = NVAR*(ivar-1)+jvar
                    rarray(jvar,ivar)%p_Ddata(ii) = rarray(jvar,ivar)%p_Ddata(ii)&
                                                  - Dcoefficients(ijpos,1,idx)
                    rarray(jvar,ivar)%p_Ddata(jj) = rarray(jvar,ivar)%p_Ddata(jj)&
                                                  - Dcoefficients(ijpos,1,idx)
                    rarray(jvar,ivar)%p_Ddata(ij) = Dcoefficients(ijpos,2,idx)&
                                                  + Dcoefficients(ijpos,1,idx)
                    rarray(jvar,ivar)%p_Ddata(ji) = Dcoefficients(ijpos,3,idx)&
                                                  + Dcoefficients(ijpos,1,idx)
                  end do
                end do
              end do

            else   ! matrix is block-diagonal
              
              do idx = 1, IEDGEmax-IEDGEset+1
                
                ! Get actual edge number
                iedge = idx+IEDGEset-1
                
                ! Get position of off-diagonal entries
                ij = IedgeList(3,iedge)
                ji = IedgeList(4,iedge)

                ! Get position of diagonal entries
                ii = IedgeList(5,iedge)
                jj = IedgeList(6,iedge)
                
                ! Update the global operator
                do ivar = 1, NVAR
                  rarray(jvar,ivar)%p_Ddata(ii) = rarray(jvar,ivar)%p_Ddata(ii)&
                                                - Dcoefficients(ivar,1,idx)
                  rarray(jvar,ivar)%p_Ddata(jj) = rarray(jvar,ivar)%p_Ddata(jj)&
                                                - Dcoefficients(ivar,1,idx)
                  rarray(ivar,ivar)%p_Ddata(ij) = Dcoefficients(ivar,2,idx)&
                                                + Dcoefficients(ivar,1,idx)
                  rarray(ivar,ivar)%p_Ddata(ji) = Dcoefficients(ivar,3,idx)&
                                                + Dcoefficients(ivar,1,idx)
                end do
              end do

            end if

          else   ! do not clear matrix

            if (bisFullMatrix) then

              do idx = 1, IEDGEmax-IEDGEset+1
                
                ! Get actual edge number
                iedge = idx+IEDGEset-1
                
                ! Get position of off-diagonal entries
                ij = IedgeList(3,iedge)
                ji = IedgeList(4,iedge)
                
                ! Update the global operator
                do ivar = 1, NVAR
                  do jvar = 1, NVAR
                    ijpos = NVAR*(ivar-1)+jvar
                    rarray(jvar,ivar)%p_Ddata(ii) = rarray(jvar,ivar)%p_Ddata(ii)&
                                                  - Dcoefficients(ijpos,1,idx)
                    rarray(jvar,ivar)%p_Ddata(jj) = rarray(jvar,ivar)%p_Ddata(jj)&
                                                  - Dcoefficients(ijpos,1,idx)
                    rarray(jvar,ivar)%p_Ddata(ij) = rarray(jvar,ivar)%p_Ddata(ij)&
                                                  + Dcoefficients(ijpos,2,idx)&
                                                  + Dcoefficients(ijpos,1,idx)
                    rarray(jvar,ivar)%p_Ddata(ji) = rarray(jvar,ivar)%p_Ddata(ji)&
                                                  + Dcoefficients(ijpos,3,idx)&
                                                  + Dcoefficients(ijpos,1,idx)
                  end do
                end do
              end do
              
            else   ! matrix is block-diagonal
              
              do idx = 1, IEDGEmax-IEDGEset+1
                
                ! Get actual edge number
                iedge = idx+IEDGEset-1
                
                ! Get position of off-diagonal entries
                ij = IedgeList(3,iedge)
                ji = IedgeList(4,iedge)
                
                ! Update the global operator
                do ivar = 1, NVAR
                  rarray(jvar,ivar)%p_Ddata(ii) = rarray(jvar,ivar)%p_Ddata(ii)&
                                                - Dcoefficients(ivar,1,idx)
                  rarray(jvar,ivar)%p_Ddata(jj) = rarray(jvar,ivar)%p_Ddata(jj)&
                                                - Dcoefficients(ivar,1,idx)
                  rarray(ivar,ivar)%p_Ddata(ij) = rarray(ivar,ivar)%p_Ddata(ij)&
                                                + Dcoefficients(ivar,2,idx)&
                                                + Dcoefficients(ivar,1,idx)
                  rarray(ivar,ivar)%p_Ddata(ji) = rarray(ivar,ivar)%p_Ddata(ji)&
                                                + Dcoefficients(ivar,3,idx)&
                                                + Dcoefficients(ivar,1,idx)
                end do
              end do
              
            end if
          end if
        end do
        !$omp end do
        
      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorStabDble

  end subroutine gfsys_buildOperatorEdgeBlock1

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildOperatorEdgeBlock2(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      dscale, bclear, rmatrix, rcollection, rafcstab,&
      fcb_calcOperatorEdgeSys, rperfconfig)

!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only edge-by-edge assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! stabilisation of AFC-type for hyperbolic problems is performed.
    ! This technique is for instance described in the reference:
    !
    !     D. Kuzmin and M. Moeller, Algebraic flux correction II. Compressible
    !     Euler equations, In: D. Kuzmin et al. (eds), Flux-Corrected
    !     Transport: Principles,  Algorithms, and Applications,
    !     Springer, 2005, 207-250.
    !
    ! Note that this routine is designed for block matrices/vectors.
    ! If there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorScalar), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute local matrices
    include 'intf_calcMatrixDiagSys_sim.inc'
    include 'intf_calcMatrixSys_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorEdgeSys.inc'
    optional :: fcb_calcOperatorEdgeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! Local variables
    type(t_vectorBlock) :: rxBlock

    ! Check if user-defined callback function is present
    if (present(fcb_calcOperatorEdgeSys)) then

      ! Create auxiliary 1-block vector
      call lsysbl_createVecFromScalar(rx, rxBlock)

      call fcb_calcOperatorEdgeSys(rgroupFEMSet, rxBlock, rmatrix, dscale,&
          bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
          rcollection, rafcstab)

      ! Release auxiliary 1-block vector
      call lsysbl_releaseVector(rxBlock)

    elseif ((rmatrix%nblocksPerCol .eq. 1) .and.&
            (rmatrix%nblocksPerRow .eq. 1)) then

      ! Call scalar version of this routine
      call gfsys_buildOperatorEdgeScalar(rgroupFEMSet, rx,&
          fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, dscale,&
          bclear, rmatrix%RmatrixBlock(1,1), rcollection, rafcstab,&
          rperfconfig=rperfconfig)
    
    else
      call output_line('Matrix must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeBlock2')
      call sys_halt()
    end if

  end subroutine gfsys_buildOperatorEdgeBlock2

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildOperatorEdgeBlock3(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      dscale, bclear, rmatrix, rcollection, rafcstab,&
      fcb_calcOperatorEdgeSys, rperfconfig)

!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only edge-by-edge assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! stabilisation of AFC-type for hyperbolic problems is performed.
    ! This technique is for instance described in the reference:
    !
    !     D. Kuzmin and M. Moeller, Algebraic flux correction II. Compressible
    !     Euler equations, In: D. Kuzmin et al. (eds), Flux-Corrected
    !     Transport: Principles,  Algorithms, and Applications,
    !     Springer, 2005, 207-250.
    !
    ! Note that this routine is designed for block matrices/vectors.
    ! If there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute local matrices
    include 'intf_calcMatrixDiagSys_sim.inc'
    include 'intf_calcMatrixSys_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorEdgeSys.inc'
    optional :: fcb_calcOperatorEdgeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! Local variables
    type(t_matrixBlock) :: rmatrixBlock

    ! Check if user-defined callback function is present
    if (present(fcb_calcOperatorEdgeSys)) then

      ! Create auxiliary 1-block matrix
      call lsysbl_createMatFromScalar(rmatrix, rmatrixBlock)

      call fcb_calcOperatorEdgeSys(rgroupFEMSet, rx, rmatrixBlock, dscale,&
          bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
          rcollection, rafcstab)

      ! Release auxiliary 1-block matrix
      call lsysbl_releaseMatrix(rmatrixBlock)

    elseif (rx%nblocks .eq. 1) then

      ! Call scalar version of this routine
      call gfsys_buildOperatorEdgeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, dscale,&
          bclear, rmatrix, rcollection, rafcstab, rperfconfig=rperfconfig)

    else
      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeBlock3')
      call sys_halt()
    end if
    
  end subroutine gfsys_buildOperatorEdgeBlock3

  ! ****************************************************************************

!<subroutine>

  subroutine gfsys_buildOperatorEdgeScalar(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
      dscale, bclear, rmatrix, rcollection, rafcstab,&
      fcb_calcOperatorEdgeSys, rperfconfig)

!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only edge-by-edge assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! stabilisation of AFC-type for hyperbolic problems is performed.
    ! This technique is for instance described in the reference:
    !
    !     D. Kuzmin and M. Moeller, Algebraic flux correction II. Compressible
    !     Euler equations, In: D. Kuzmin et al. (eds), Flux-Corrected
    !     Transport: Principles,  Algorithms, and Applications,
    !     Springer, 2005, 207-250.
    !
    ! Note that this routine requires the scalar matrices/vectors
    ! are stored in the interleave format.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! Vector on which the matrix entries may depend.
    type(t_vectorScalar), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute local matrices
    include 'intf_calcMatrixDiagSys_sim.inc'
    include 'intf_calcMatrixSys_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorEdgeSys.inc'
    optional :: fcb_calcOperatorEdgeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rxBlock
    type(t_matrixBlock) :: rmatrixBlock
    real(DP), dimension(:), pointer :: p_Ddata,p_Dx
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtDiag   
    integer, dimension(:,:), pointer :: p_IdiagList,p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if user-defined assembly is provided
    if (present(fcb_calcOperatorEdgeSys)) then
      ! Create auxiliary 1-block vector and matrix
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createMatFromScalar(rmatrix, rmatrixBlock)
      
      ! Call user-defined assembly
      call fcb_calcOperatorEdgeSys(rgroupFEMSet, rxBlock, rmatrixBlock, dscale,&
          bclear, fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
          rcollection, rafcstab)
      
      ! Release auxiliary 1-block vector and matrix
      call lsysbl_releaseVector(rxBlock)
      call lsysbl_releaseMatrix(rmatrixBlock)

      ! That`s it
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => gfsys_perfconfig
    end if
    
    ! Check if matrix and vector have the same data type
    if ((rmatrix%cdataType .ne. rx%cdataType) .or.&
        (rmatrix%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeScalar')
      call sys_halt()
    end if

    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_EDGEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------
      
      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGLIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGDATA) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA) .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeScalar')
        call sys_halt()
      end if

      ! Set pointers
      call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
      call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)
      call gfem_getbase_IdiagList(rgroupFEMSet, p_IdiagList)

      ! What data types are we?
      select case(rgroupFEMSet%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtDiag(rgroupFEMSet, p_DcoeffsAtDiag)
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        call lsyssc_getbase_double(rmatrix, p_Ddata)
        call lsyssc_getbase_double(rx, p_Dx)
        
        ! What type of matrix are we?
        select case(rmatrix%cinterleavematrixFormat)
          
        case (LSYSSC_MATRIX1)
          ! Assemble matrix diagonal
          call doOperatorDiagDble(rx%NEQ, rmatrix%NA, rx%NVAR, rx%NVAR*rx%NVAR,&
              p_IdiagList, p_DcoeffsAtDiag, p_Dx, dscale, bclear, p_Ddata)
          
          ! Do we have to build the stabilisation?
          if (present(rafcstab)) then
            
            ! Check if stabilisation has been prepared
            if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
              call output_line('Stabilisation has not been prepared!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeScalar')
              call sys_halt()
            end if
            
            !-------------------------------------------------------------------
            ! Assemble operator with stabilisation
            !-------------------------------------------------------------------
            call doOperatorStabDble(rx%NEQ, rmatrix%NA, rx%NVAR, rx%NVAR*rx%NVAR,&
                p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge, p_Dx,&
                dscale, bclear, p_Ddata)
          
          else   ! no stabilisation structure present
            
            !-------------------------------------------------------------------
            ! Assemble operator without stabilisation
            !-------------------------------------------------------------------
            call doOperatorDble(rx%NEQ, rmatrix%NA, rx%NVAR, rx%NVAR*rx%NVAR,&
                p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge, p_Dx,&
                dscale, bclear, p_Ddata)
          end if

        case (LSYSSC_MATRIXD)
          ! Assemble matrix diagonal
          call doOperatorDiagDble(rx%NEQ, rmatrix%NA, rx%NVAR, rx%NVAR,&
              p_IdiagList, p_DcoeffsAtDiag, p_Dx, dscale, bclear, p_Ddata)
          
          ! Do we have to build the stabilisation?
          if (present(rafcstab)) then
            
            ! Check if stabilisation has been prepared
            if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
              call output_line('Stabilisation has not been prepared!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeScalar')
              call sys_halt()
            end if
            
            !-------------------------------------------------------------------
            ! Assemble operator with stabilisation
            !-------------------------------------------------------------------
            call doOperatorStabDble(rx%NEQ, rmatrix%NA, rx%NVAR, rx%NVAR,&
                p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge, p_Dx,&
                dscale, bclear, p_Ddata)
          
          else   ! no stabilisation structure present
            
            !-------------------------------------------------------------------
            ! Assemble operator without stabilisation
            !-------------------------------------------------------------------
            call doOperatorDble(rx%NEQ, rmatrix%NA, rx%NVAR, rx%NVAR,&
                p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge, p_Dx,&
                dscale, bclear, p_Ddata)
          end if

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeScalar')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeScalar')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorEdgeScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble diagonal part of the operator
    
    subroutine doOperatorDiagDble(NEQ, NA, NVAR, MVAR, IdiagList,&
        DcoeffsAtDiag, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtDiag
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IdiagList
      integer, intent(in) :: NEQ,NA,NVAR,MVAR

      ! input/output parameters
      real(DP), dimension(MVAR,NA), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: Dcoefficients

      ! local variables
      integer :: IEQmax,IEQset,i,ia,idx

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,IEQmax,i,ia,idx)&
      !$omp if (size(IdiagList,2) > p_rperfconfig%NEQMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtNode(NVAR,p_rperfconfig%NEQSIM))
      allocate(Dcoefficients(MVAR,1,p_rperfconfig%NEQSIM))

      ! Loop over the equations
      !$omp do schedule(static,1)
      do IEQset = 1, size(IdiagList,2), p_rperfconfig%NEQSIM

        ! We always handle NEQSIM equations simultaneously.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most NEQSIM equations simultaneously.
        
        IEQmax = min(size(IdiagList,2), IEQset-1+p_rperfconfig%NEQSIM)
        
        ! Loop through all equations in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IEQmax-IEQset+1

          ! Get actual equation number
          i = IdiagList(1,idx+IEQset-1)
          
          ! Fill auxiliary array
          DdataAtNode(:,idx) = Dx(:,i)
        end do

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DcoeffsAtDiag(:,IEQset:IEQmax),&
            IdiagList(:,IEQset:IEQmax),&
            dscale, IEQmax-IEQset+1,&
            Dcoefficients(:,:,1:IEQmax-IEQset+1), rcollection)

        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          
          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ia = IdiagList(2,idx+IEQset-1)
            
            ! Update the diagonal coefficient
            Ddata(:,ia) = Dcoefficients(:,1,idx)
          end do

        else   ! do not clear matrix

          do idx = 1, IEQmax-IEQset+1

            ! Get position of diagonal entry
            ia = IdiagList(2,idx+IEQset-1)
            
            ! Update the diagonal coefficient
            Ddata(:,ia) = Ddata(:,ia) + Dcoefficients(:,1,idx)
          end do

        end if
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel

    end subroutine doOperatorDiagDble

    !**************************************************************
    ! Assemble operator edge-by-edge without stabilisation
    
    subroutine doOperatorDble(NEQ, NA, NVAR, MVAR, IedgeListIdx,&
        IedgeList, DcoeffsAtEdge, Dx, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NA,NVAR,MVAR

      ! input/output parameters
      real(DP), dimension(MVAR,NA), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: Dcoefficients

      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ij,ji

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtEdge,IEDGEmax,idx,iedge,ij,ji)&
      !$omp if(size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(Dcoefficients(MVAR,2,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+p_rperfconfig%NEDGESIM)
        
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dcoefficients(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global matrix
          if (bclear) then

            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Update the global operator
              Ddata(:,ij) = Dcoefficients(:,1,idx)
              Ddata(:,ji) = Dcoefficients(:,2,idx)
            end do

          else   ! do not clear matrix

            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Update the global operator
              Ddata(:,ij) = Ddata(:,ij) + Dcoefficients(:,1,idx)
              Ddata(:,ji) = Ddata(:,ji) + Dcoefficients(:,2,idx)
            end do
            
          end if
        end do
        !$omp end do
        
      end do ! igroup
      
      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorDble

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation
    
    subroutine doOperatorStabDble(NEQ, NA, NVAR, MVAR, IedgeListIdx,&
        IedgeList, DcoeffsAtEdge, Dx, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NA,NVAR,MVAR

      ! input/output parameters
      real(DP), dimension(MVAR,NA), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: Dcoefficients

      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ij,ji,ii,jj

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtEdge,IEDGEmax,idx,iedge,ij,ji,ii,jj)&
      !$omp if(size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(Dcoefficients(MVAR,3,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+p_rperfconfig%NEDGESIM)
        
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dcoefficients(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global matrix
          if (bclear) then

            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Get position of diagonal entries
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Update the global operator
              Ddata(:,ii) = Ddata(:,ii) - Dcoefficients(:,1,idx)
              Ddata(:,jj) = Ddata(:,jj) - Dcoefficients(:,1,idx)
              Ddata(:,ij) = Dcoefficients(:,2,idx) + Dcoefficients(:,1,idx)
              Ddata(:,ji) = Dcoefficients(:,3,idx) + Dcoefficients(:,1,idx)
            end do

          else   ! do not clear matrix

            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)

              ! Get position of diagonal entries
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Update the global operator
              Ddata(:,ii) = Ddata(:,ii) - Dcoefficients(:,1,idx)
              Ddata(:,jj) = Ddata(:,jj) - Dcoefficients(:,1,idx)
              Ddata(:,ij) = Ddata(:,ij) + Dcoefficients(:,2,idx)&
                                        + Dcoefficients(:,1,idx)
              Ddata(:,ji) = Ddata(:,ji) + Dcoefficients(:,3,idx)&
                                        + Dcoefficients(:,1,idx)
            end do
            
          end if
        end do
        !$omp end do
        
      end do ! igroup
      
      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorStabDble

  end subroutine gfsys_buildOperatorEdgeScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildVectorNodeBlock(rgroupFEMSet, rx,&
      fcb_calcVectorSys_sim, dscale, bclear, rvector,&
      rcollection, rafcstab, fcb_calcVectorNodeSys, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete vector by the group
    ! finite element formulation. The vector entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling a user-defined callback function.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! This subroutine assembles the divergence vector for block vectors.
    ! If the vector contains only one block, then the scalar
    ! counterpart of this routine is called with the scalar subvector.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! Vector on which the vector entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FALSE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute vector entries
    include 'intf_calcVectorSys_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcVectorNodeSys.inc'
    optional :: fcb_calcVectorNodeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Destination vector
    type(t_VectorBlock), intent(inout) :: rvector

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    real(DP), dimension(:), pointer :: p_Dx,p_Ddata
    integer, dimension(:,:), pointer :: p_InodeList2D
    integer, dimension(:), pointer :: p_InodeListIdx,p_InodeList1D

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if user-defined callback function is present
    if (present(fcb_calcVectorNodeSys)) then

      call fcb_calcVectorNodeSys(rgroupFEMSet, rx, rvector, dscale,&
          bclear, fcb_calcVectorSys_sim, rcollection, rafcstab)
      ! That`s it
      return
      
    elseif ((rx%nblocks .eq. 1) .and. (rvector%nblocks .eq. 1)) then

      ! Call scalar version of this routine
      call gfsys_buildVectorNodeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcVectorSys_sim, dscale, bclear, rvector%RvectorBlock(1),&
          rcollection, rafcstab, rperfconfig=rperfconfig)
      ! That`s it
      return
    end if
    
    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => gfsys_perfconfig
    end if

    ! Check if vectors have the same data type double
    if ((rx%cdataType .ne. rvector%cdataType) .or.&
        (rx%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorNodeBlock')
      call sys_halt()
    end if
    
    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_NODEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------

      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODELIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA) .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorNodeBlock')
        call sys_halt()
      end if
      
      ! What data types are we?
      select case(rvector%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
        call lsysbl_getbase_double(rvector, p_Ddata)
        call lsysbl_getbase_double(rx, p_Dx)
        
        ! Check if only a subset of the vector is required
        if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
          ! Set pointers
          call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
          call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList1D)

          ! Assemble vector node-by-node
          call doVectorDble(rx%RvectorBlock(1)%NEQ, rx%nblocks,&
              p_InodeListIdx, p_InodeList1D, p_DcoeffsAtNode, p_Dx,&
              dscale, bclear, p_Ddata)
        else
          ! Set pointers
          call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
          call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList2D)

          ! Assemble selected part of the vector node-by-node
          call doVectorDbleSel(rx%RvectorBlock(1)%NEQ, rx%nblocks,&
              p_InodeListIdx, p_InodeList2D, p_DcoeffsAtNode, p_Dx,&
              dscale, bclear, p_Ddata)
        end if
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorNodeBlock')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorNodeBlock')
      call sys_halt()
    end select

  contains
    
    !**************************************************************
    ! Assemble vector node-by-node without stabilisation

    subroutine doVectorDble(NEQ, NVAR, InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: InodeListIdx,InodeList
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode,Dcoefficients
      integer, dimension(:,:), pointer  :: IdofsAtNode
      
      ! local variables
      integer :: IAmax,IApos,IAset,IEQmax,IEQset,ia,idx,ieq

      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------
      
      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,IdofsAtNode,&
      !$omp         IAmax,IApos,IAset,IEQmax,ia,idx,ieq)&
      !$omp if(size(InodeList) > p_rperfconfig%NAMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtNode(NVAR,p_rperfconfig%NEQSIM))
      allocate(Dcoefficients(NVAR,p_rperfconfig%NEQSIM))
      allocate(IdofsAtNode(2,p_rperfconfig%NEQSIM))

      ! Loop over all equations in blocks of size NEQSIM
      !$omp do schedule(static,1)
      do IEQset = 1, size(InodeListIdx)-1, p_rperfconfig%NEQSIM
        
        ! We always handle NEQSIM equations by one OpenMP thread.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most NEQSIM equations by one OpenMP thread.

        IEQmax = min(size(InodeListIdx)-1, IEQset-1+p_rperfconfig%NEQSIM)
        
        ! Since the number of nonzero entries per equation is not
        ! fixed we need to iterate over all equations in the current
        ! set of equations [IEQset:IEQmax] and process not more than
        ! NASIM nonzero entries simultaneously.
        
        ! Initialise the lower and upper bounds of nonzero entries;
        ! note that this is not the matrix position itself but the
        ! equation number to which the nonzero matrix entries belong
        IAset = IEQset
        IAmax = IEQset

        ! Also initialise the absolute position of first nonzero entry
        IApos = InodeListIdx(IEQset)
        
        ! Repeat until all equation in the current set have been processed
        do while (IAset .le. IEQmax)
          
          ! Initialise local index which will run from IAset..IAmax.
          ! Since IAset and IAmax are not fixed but they are updated
          ! step-by-step we need this complicated while-loop
          idx = 0

          ! Loop over the equations of the current group and include
          ! an equation into the current set of nonzero matrix entries
          ! if the total number of nonzero matrix entries in the set
          ! does not exceed the upper bound NASIM
          do while(IAmax .le. IEQmax)

            ! Exit if more than NASIM nonzero entries would be processed
            if (InodeListIdx(IAmax+1)-IApos .gt. p_rperfconfig%NASIM) exit
            
            ! Loop through all nonzero matrix entries in the current
            ! equation and prepare the auxiliary arrays for it
            do ia = InodeListIdx(IAmax), InodeListIdx(IAmax+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Fill auxiliary arrays
              IdofsAtNode(1,idx) = InodeList(ia)       ! absolut nodal value j
              IdofsAtNode(2,idx) = ia                  ! absolute matrix position ia
              DdataAtNode(:,idx) = Dx(InodeList(ia),:) ! solution value at node j
            end do
            
            ! Increase the upper bound for nonzero entries
            IAmax = IAmax+1
          end do
          
          ! Use callback function to compute matrix entries
          call fcb_calcVectorSys_sim(&
              DdataAtNode(:,1:idx),&
              DcoeffsAtNode(:,IApos:IApos+idx-1),&
              IdofsAtNode(:,1:idx),&
              dscale, idx, Dcoefficients(:,1:idx), rcollection)
          
          ! Initialise local index which will run from IAset..IAmax
          idx = 0

          ! Loop through all equations in the current set
          ! and scatter the entries to the global matrix
          do ieq = IAset, IAmax-1
            
            ! Loop over all contributions to this equation
            do ia = InodeListIdx(ieq), InodeListIdx(ieq+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Update the global vector
              Ddata(ieq,:) = Ddata(ieq,:) + Dcoefficients(:,idx)
            end do
          end do
          
          ! Proceed with next nonzero entries in current set
          IAset = IAmax
          IApos = InodeListIdx(IASet)
          
        end do
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(IdofsAtNode)
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel

    end subroutine doVectorDble
    
    !**************************************************************
    ! Assemble vector node-by-node without stabilisation

    subroutine doVectorDbleSel(NEQ, NVAR, InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: InodeListIdx
      integer, dimension(:,:), intent(in) :: InodeList
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode,Dcoefficients
      
      ! local variables
      integer :: IAmax,IApos,IAset,IEQmax,IEQset,i,ia,idx,ieq

      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------
      
      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,&
      !$omp         IAmax,IApos,IAset,IEQmax,i,ia,idx,ieq)&
      !$omp if(size(InodeList,2) > p_rperfconfig%NAMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtNode(NVAR,p_rperfconfig%NEQSIM))
      allocate(Dcoefficients(NVAR,p_rperfconfig%NEQSIM))

      ! Loop over all equations in blocks of size NEQSIM
      !$omp do schedule(static,1)
      do IEQset = 1, size(InodeListIdx)-1, p_rperfconfig%NEQSIM
        
        ! We always handle NEQSIM equations by one OpenMP thread.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most NEQSIM equations by one OpenMP thread.

        IEQmax = min(size(InodeListIdx)-1, IEQset-1+p_rperfconfig%NEQSIM)
        
        ! Since the number of nonzero entries per equation is not
        ! fixed we need to iterate over all equations in the current
        ! set of equations [IEQset:IEQmax] and process not more than
        ! NASIM nonzero entries simultaneously.
        
        ! Initialise the lower and upper bounds of nonzero entries;
        ! note that this is not the matrix position itself but the
        ! equation number to which the nonzero matrix entries belong
        IAset = IEQset
        IAmax = IEQset

        ! Also initialise the absolute position of first nonzero entry
        IApos = InodeListIdx(IEQset)
        
        ! Repeat until all equation in the current set have been processed
        do while (IAset .le. IEQmax)
          
          ! Initialise local index which will run from IAset..IAmax.
          ! Since IAset and IAmax are not fixed but they are updated
          ! step-by-step we need this complicated while-loop
          idx = 0

          ! Loop over the equations of the current group and include
          ! an equation into the current set of nonzero matrix entries
          ! if the total number of nonzero matrix entries in the set
          ! does not exceed the upper bound NASIM
          do while(IAmax .le. IEQmax)

            ! Exit if more than NASIM nonzero entries would be processed
            if (InodeListIdx(IAmax+1)-IApos .gt. p_rperfconfig%NASIM) exit
            
            ! Loop through all nonzero matrix entries in the current
            ! equation and prepare the auxiliary arrays for it
            do ia = InodeListIdx(IAmax), InodeListIdx(IAmax+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Fill auxiliary arrays
              DdataAtNode(:,idx) = Dx(InodeList(1,ia),:) ! solution value at node j
            end do
            
            ! Increase the upper bound for nonzero entries
            IAmax = IAmax+1
          end do
          
          ! Use callback function to compute matrix entries
          call fcb_calcVectorSys_sim(&
              DdataAtNode(:,1:idx),&
              DcoeffsAtNode(:,IApos:IApos+idx-1),&
              InodeList(:,IApos:IApos+idx-1),&
              dscale, idx, Dcoefficients(:,1:idx), rcollection)
          
          ! Initialise local index which will run from IAset..IAmax
          idx = 0

          ! Loop through all equations in the current set
          ! and scatter the entries to the global matrix
          do i = IAset, IAmax-1
            
            ! Get actual node number
            ieq = InodeList(1,InodeListIdx(i))
            
            ! Loop over all contributions to this equation
            do ia = InodeListIdx(i), InodeListIdx(i+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Update the global vector
              Ddata(ieq,:) = Ddata(ieq,:) + Dcoefficients(:,idx)
            end do
          end do
          
          ! Proceed with next nonzero entries in current set
          IAset = IAmax
          IApos = InodeListIdx(IASet)
          
        end do
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel

    end subroutine doVectorDbleSel

  end subroutine gfsys_buildVectorNodeBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildVectorNodeScalar(rgroupFEMSet, rx,&
      fcb_calcVectorSys_sim, dscale, bclear, rvector,&
      rcollection, rafcstab, fcb_calcVectorNodeSys, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete vector by the group
    ! finite element formulation. The vector entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling a user-defined callback function.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! This subroutine assembles the divergence vector for scalar vectors.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! Vector on which the vector entries may depend.
    type(t_vectorScalar), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FALSE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute vector entries
    include 'intf_calcVectorSys_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcVectorNodeSys.inc'
    optional :: fcb_calcVectorNodeSys
    
    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Destination vector
    type(t_VectorScalar), intent(inout) :: rvector

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rxBlock,rvectorBlock
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    real(DP), dimension(:), pointer :: p_Dx,p_Ddata
    integer, dimension(:,:), pointer :: p_InodeList2D
    integer, dimension(:), pointer :: p_InodeListIdx,p_InodeList1D

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if user-defined assembly is provided
    if (present(fcb_calcVectorNodeSys)) then
      ! Create auxiliary 1-block vectors
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createVecFromScalar(rvector, rvectorBlock)
      
      ! Call user-defined assembly
      call fcb_calcVectorNodeSys(rgroupFEMSet, rxBlock, rvectorBlock,&
          dscale, bclear, fcb_calcVectorSys_sim, rcollection, rafcstab)
      
      ! Release auxiliary 1-block vectors
      call lsysbl_releaseVector(rxBlock)
      call lsysbl_releaseVector(rvectorBlock)

      ! That`s it
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => gfsys_perfconfig
    end if

    ! Check if vectors have the same data type double
    if ((rx%cdataType .ne. rvector%cdataType) .or.&
        (rx%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorNodeScalar')
      call sys_halt()
    end if
    
    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_NODEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------

      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODELIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA) .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorNodeScalar')
        call sys_halt()
      end if
      
      ! What data types are we?
      select case(rvector%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
        call lsyssc_getbase_double(rvector, p_Ddata)
        call lsyssc_getbase_double(rx, p_Dx)
        
        ! Check if only a subset of the vector is required
        if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
          ! Set pointers
          call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
          call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList1D)

          ! Assemble vector node-by-node
          call doVectorDble(rx%NEQ, rx%NVAR, p_InodeListIdx, p_InodeList1D,&
              p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)
        else
          ! Set pointers
          call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx)
          call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList2D)

          ! Assemble selected part of the vector node-by-node
          call doVectorDbleSel(rx%NEQ, rx%NVAR, p_InodeListIdx, p_InodeList2D,&
              p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)
        end if
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorNodeScalar')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorNodeScalar')
      call sys_halt()
    end select

  contains
    
    !**************************************************************
    ! Assemble vector node-by-node without stabilisation

    subroutine doVectorDble(NEQ, NVAR, InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: InodeListIdx,InodeList
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode,Dcoefficients
      integer, dimension(:,:), pointer  :: IdofsAtNode
      
      ! local variables
      integer :: IAmax,IApos,IAset,IEQmax,IEQset,ia,idx,ieq

      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------
      
      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,IdofsAtNode,&
      !$omp         IAmax,IApos,IAset,IEQmax,ia,idx,ieq)&
      !$omp if(size(InodeList) > p_rperfconfig%NAMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtNode(NVAR,p_rperfconfig%NEQSIM))
      allocate(Dcoefficients(NVAR,p_rperfconfig%NEQSIM))
      allocate(IdofsAtNode(2,p_rperfconfig%NEQSIM))

      ! Loop over all equations in blocks of size NEQSIM
      !$omp do schedule(static,1)
      do IEQset = 1, size(InodeListIdx)-1, p_rperfconfig%NEQSIM
        
        ! We always handle NEQSIM equations by one OpenMP thread.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most NEQSIM equations by one OpenMP thread.

        IEQmax = min(size(InodeListIdx)-1, IEQset-1+p_rperfconfig%NEQSIM)
        
        ! Since the number of nonzero entries per equation is not
        ! fixed we need to iterate over all equations in the current
        ! set of equations [IEQset:IEQmax] and process not more than
        ! NASIM nonzero entries simultaneously.
        
        ! Initialise the lower and upper bounds of nonzero entries;
        ! note that this is not the matrix position itself but the
        ! equation number to which the nonzero matrix entries belong
        IAset = IEQset
        IAmax = IEQset

        ! Also initialise the absolute position of first nonzero entry
        IApos = InodeListIdx(IEQset)
        
        ! Repeat until all equation in the current set have been processed
        do while (IAset .le. IEQmax)
          
          ! Initialise local index which will run from IAset..IAmax.
          ! Since IAset and IAmax are not fixed but they are updated
          ! step-by-step we need this complicated while-loop
          idx = 0

          ! Loop over the equations of the current group and include
          ! an equation into the current set of nonzero matrix entries
          ! if the total number of nonzero matrix entries in the set
          ! does not exceed the upper bound NASIM
          do while(IAmax .le. IEQmax)

            ! Exit if more than NASIM nonzero entries would be processed
            if (InodeListIdx(IAmax+1)-IApos .gt. p_rperfconfig%NASIM) exit
            
            ! Loop through all nonzero matrix entries in the current
            ! equation and prepare the auxiliary arrays for it
            do ia = InodeListIdx(IAmax), InodeListIdx(IAmax+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Fill auxiliary arrays
              IdofsAtNode(1,idx) = InodeList(ia)       ! absolut nodal value j
              IdofsAtNode(2,idx) = ia                  ! absolute matrix position ia
              DdataAtNode(:,idx) = Dx(:,InodeList(ia)) ! solution value at node j
            end do
            
            ! Increase the upper bound for nonzero entries
            IAmax = IAmax+1
          end do
          
          ! Use callback function to compute matrix entries
          call fcb_calcVectorSys_sim(&
              DdataAtNode(:,1:idx),&
              DcoeffsAtNode(:,IApos:IApos+idx-1),&
              IdofsAtNode(:,1:idx),&
              dscale, idx, Dcoefficients(:,1:idx), rcollection)
          
          ! Initialise local index which will run from IAset..IAmax
          idx = 0

          ! Loop through all equations in the current set
          ! and scatter the entries to the global matrix
          do ieq = IAset, IAmax-1
            
            ! Loop over all contributions to this equation
            do ia = InodeListIdx(ieq), InodeListIdx(ieq+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Update the global vector
              Ddata(:,ieq) = Ddata(:,ieq) + Dcoefficients(:,idx)
            end do
          end do
          
          ! Proceed with next nonzero entries in current set
          IAset = IAmax
          IApos = InodeListIdx(IASet)
          
        end do
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(IdofsAtNode)
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel

    end subroutine doVectorDble
    
    !**************************************************************
    ! Assemble vector node-by-node without stabilisation

    subroutine doVectorDbleSel(NEQ, NVAR, InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: InodeListIdx
      integer, dimension(:,:), intent(in) :: InodeList
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtNode,Dcoefficients
      
      ! local variables
      integer :: IAmax,IApos,IAset,IEQmax,IEQset,i,ia,idx,ieq

      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------
      
      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,&
      !$omp         IAmax,IApos,IAset,IEQmax,i,ia,idx,ieq)&
      !$omp if(size(InodeList,2) > p_rperfconfig%NAMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtNode(NVAR,p_rperfconfig%NEQSIM))
      allocate(Dcoefficients(NVAR,p_rperfconfig%NEQSIM))

      ! Loop over all equations in blocks of size NEQSIM
      !$omp do schedule(static,1)
      do IEQset = 1, size(InodeListIdx)-1, p_rperfconfig%NEQSIM
        
        ! We always handle NEQSIM equations by one OpenMP thread.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most NEQSIM equations by one OpenMP thread.

        IEQmax = min(size(InodeListIdx)-1, IEQset-1+p_rperfconfig%NEQSIM)
        
        ! Since the number of nonzero entries per equation is not
        ! fixed we need to iterate over all equations in the current
        ! set of equations [IEQset:IEQmax] and process not more than
        ! NASIM nonzero entries simultaneously.
        
        ! Initialise the lower and upper bounds of nonzero entries;
        ! note that this is not the matrix position itself but the
        ! equation number to which the nonzero matrix entries belong
        IAset = IEQset
        IAmax = IEQset

        ! Also initialise the absolute position of first nonzero entry
        IApos = InodeListIdx(IEQset)
        
        ! Repeat until all equation in the current set have been processed
        do while (IAset .le. IEQmax)
          
          ! Initialise local index which will run from IAset..IAmax.
          ! Since IAset and IAmax are not fixed but they are updated
          ! step-by-step we need this complicated while-loop
          idx = 0

          ! Loop over the equations of the current group and include
          ! an equation into the current set of nonzero matrix entries
          ! if the total number of nonzero matrix entries in the set
          ! does not exceed the upper bound NASIM
          do while(IAmax .le. IEQmax)

            ! Exit if more than NASIM nonzero entries would be processed
            if (InodeListIdx(IAmax+1)-IApos .gt. p_rperfconfig%NASIM) exit
            
            ! Loop through all nonzero matrix entries in the current
            ! equation and prepare the auxiliary arrays for it
            do ia = InodeListIdx(IAmax), InodeListIdx(IAmax+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Fill auxiliary arrays
              DdataAtNode(:,idx) = Dx(:,InodeList(1,ia)) ! solution value at node j
            end do
            
            ! Increase the upper bound for nonzero entries
            IAmax = IAmax+1
          end do
          
          ! Use callback function to compute matrix entries
          call fcb_calcVectorSys_sim(&
              DdataAtNode(:,1:idx),&
              DcoeffsAtNode(:,IApos:IApos+idx-1),&
              InodeList(:,IApos:IApos+idx-1),&
              dscale, idx, Dcoefficients(:,1:idx), rcollection)
          
          ! Initialise local index which will run from IAset..IAmax
          idx = 0

          ! Loop through all equations in the current set
          ! and scatter the entries to the global matrix
          do i = IAset, IAmax-1
            
            ! Get actual node number
            ieq = InodeList(1,InodeListIdx(i))
            
            ! Loop over all contributions to this equation
            do ia = InodeListIdx(i), InodeListIdx(i+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Update the global vector
              Ddata(:,ieq) = Ddata(:,ieq) + Dcoefficients(:,idx)
            end do
          end do
          
          ! Proceed with next nonzero entries in current set
          IAset = IAmax
          IApos = InodeListIdx(IASet)
          
        end do
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel

    end subroutine doVectorDbleSel
    
  end subroutine gfsys_buildVectorNodeScalar

  ! ****************************************************************************

!<subroutine>

  subroutine gfsys_buildVectorEdgeBlock(rgroupFEMSet, rx,&
      fcb_calcFlux_sim, dscale, bclear, rvector,&
      rcollection, rafcstab,  fcb_calcVectorEdgeSys, rperfconfig)

!<description>
    ! This subroutine assembles a vector operator by the group
    ! finite element formulation. The vector entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling a user-defined callback function.
    !
    ! If the vector contains only one block, then the scalar
    ! counterpart of this routine is called with the scalar subvector.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the vector entries may depend.
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! Callback function to compute local fluxes
    include 'intf_calcFlux_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcVectorEdgeSys.inc'
    optional :: fcb_calcVectorEdgeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dx,p_Ddata
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if user-defined assembly is provided
    if (present(fcb_calcVectorEdgeSys)) then

      call fcb_calcVectorEdgeSys(rgroupFEMSet, rx, rvector,&
          dscale, bclear, fcb_calcFlux_sim, rcollection, rafcstab)
      ! That`s it
      return
    
    elseif ((rx%nblocks .eq. 1) .and. (rvector%nblocks .eq. 1)) then

      call gfsys_buildVectorEdgeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcFlux_sim, dscale, bclear, rvector%RvectorBlock(1),&
          rcollection, rafcstab, rperfconfig=rperfconfig)
      ! That`s it
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => gfsys_perfconfig
    end if

    ! Check if vectors have the same data type double
    if ((rx%cdataType .ne. rvector%cdataType) .or.&
        (rx%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorEdgeBlock')
      call sys_halt()
    end if

    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_EDGEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------
      
      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA) .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorEdgeBlock')
        call sys_halt()
      end if
        
      ! Set pointers
      call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
      call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)

      ! What data types are we?
      select case(rvector%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        call lsysbl_getbase_double(rx, p_Dx)
        call lsysbl_getbase_double(rvector, p_Ddata)

        ! Assemble vector
        call doVectorDble(rx%RvectorBlock(1)%NEQ, rx%nblocks, p_IedgeListIdx,&
            p_IedgeList, p_DcoeffsAtEdge, p_Dx, dscale, bclear, p_Ddata)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorEdgeBlock')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorEdgeBlock')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble vector edge-by-edge
    
    subroutine doVectorDble(NEQ, NVAR, IedgeListIdx, IedgeList,&
        DcoeffsAtEdge, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Ddata

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DfluxesAtEdge
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DfluxesAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do
          
          ! Use callback function to compute internodal fluxes
          call fcb_calcFlux_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Update the global vector
            Ddata(i,:) = Ddata(i,:)+DfluxesAtEdge(:,1,idx)
            Ddata(j,:) = Ddata(j,:)+DfluxesAtEdge(:,2,idx)
          end do
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      !$omp end parallel
      
    end subroutine doVectorDble
    
  end subroutine gfsys_buildVectorEdgeBlock

  ! ****************************************************************************

!<subroutine>

  subroutine gfsys_buildVectorEdgeScalar(rgroupFEMSet, rx,&
      fcb_calcFlux_sim, dscale, bclear, rvector,&
      rcollection, rafcstab, fcb_calcVectorEdgeSys, rperfconfig)

!<description>
    ! This subroutine assembles a vector operator by the group
    ! finite element formulation. The vector entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling a user-defined callback function.
    !
    ! Note that the vectors are required as scalar vectors which are
    ! stored in the interleave format.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the vector entries may depend.
    type(t_vectorScalar), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute local fluxes
    include 'intf_calcFlux_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcVectorEdgeSys.inc'
    optional :: fcb_calcVectorEdgeSys

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorScalar), intent(inout) :: rvector

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rxBlock,rvectorBlock
    real(DP), dimension(:), pointer :: p_Dx,p_Ddata
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if user-defined assembly is provided
    if (present(fcb_calcVectorEdgeSys)) then
      ! Create auxiliary 1-block vector
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createVecFromScalar(rvector, rvectorBlock)
      
      ! Call user-defined assembly
      call fcb_calcVectorEdgeSys(rgroupFEMSet, rxBlock, rvectorBlock,&
          dscale, bclear, fcb_calcFlux_sim, rcollection, rafcstab)
      
      ! Release auxiliary 1-block vectors
      call lsysbl_releaseVector(rxBlock)
      call lsysbl_releaseVector(rvectorBlock)

      ! That`s it
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => gfsys_perfconfig
    end if

    ! Check if vectors have the same data type double
    if ((rx%cdataType .ne. rvector%cdataType) .or.&
        (rx%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorEdgeScalar')
      call sys_halt()
    end if

    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_EDGEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------
      
      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA) .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorEdgeScalar')
        call sys_halt()
      end if
        
      ! Set pointers
      call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
      call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)

      ! What data types are we?
      select case(rvector%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        call lsyssc_getbase_double(rx, p_Dx)
        call lsyssc_getbase_double(rvector, p_Ddata)

        ! Assemble vector
        call doVectorDble(rx%NEQ, rx%NVAR, p_IedgeListIdx, p_IedgeList,&
            p_DcoeffsAtEdge, p_Dx, dscale, bclear, p_Ddata)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorEdgeScalar')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorEdgeScalar')
      call sys_halt()
    end select
    
  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble vector edge-by-edge

    subroutine doVectorDble(NEQ, NVAR, IedgeListIdx, IedgeList,&
        DcoeffsAtEdge, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Ddata

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DfluxesAtEdge
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DfluxesAtEdge,IEDGEmax,i,idx,iedge,j)&
      !$omp if(size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute internodal fluxes
          call fcb_calcFlux_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Update the global vector
            Ddata(:,i) = Ddata(:,i)+DfluxesAtEdge(:,1,idx)
            Ddata(:,j) = Ddata(:,j)+DfluxesAtEdge(:,2,idx)
          end do
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      !$omp end parallel
      
    end subroutine doVectorDble

  end subroutine gfsys_buildVectorEdgeScalar
 
end module groupfemsystem
