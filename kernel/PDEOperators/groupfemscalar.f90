!##############################################################################
!# ****************************************************************************
!# <name> groupfemscalar </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the
!# group-finite element formulation to scalar problems,
!# i.e. conservation laws. The technique was proposed by
!# C.A.J. Fletcher in:
!#
!#     C.A.J. Fletcher, The group finite element formulation
!#     Computer Methods in Applied Mechanics and Engineering (ISSN
!#     0045-7825), vol. 37, April 1983, p. 225-244.
!#
!# The group finite element formulation uses the same basis functions
!# for the unknown solution and the fluxes. This allows for an
!# efficient matrix assembly, whereby the constant coefficient
!# matrices can be assembled once and for all at the beginning of the
!# simulation and each time the grid is modified.
!#
!# The following routines are available:
!#
!# 1.) gfsc_buildOperatorNode = gfsc_buildOperatorConst /
!#                              gfsc_buildOperatorNodeScalar /
!#                              gfsc_buildOperatorNodeBlock1 /
!#                              gfsc_buildOperatorNodeBlock2 /
!#                              gfsc_buildOperatorNodeBlock3
!#     -> Assembles a discrete operator node-by-node by the
!#        group finite element formulation
!#
!# 2.) gfsc_buildOperatorEdge = gfsc_buildOperatorEdgeScalar /
!#                              gfsc_buildOperatorEdgeBlock1 /
!#                              gfsc_buildOperatorEdgeBlock2 /
!#                              gfsc_buildOperatorEdgeBlock3
!#     -> Assembles a discrete operator edge-by-edge by the
!#        group finite element formulation
!#
!# 3.) gfsc_buildVectorNode = gfsc_buildVectorNodeScalar /
!#                            gfsc_buildVectorNodeBlock
!#     -> Assembles a discrete vector node-by-node by the
!#        group finite element formulation
!#
!# 4.) gfsc_buildVectorEdge = gfsc_buildVectorEdgeScalar /
!#                            gfsc_buildVectorEdgeBlock
!#     -> Assembles a discrete vector edge-by-edge by the
!#        group finite element formulation
!#
!# 5.) gfsc_buildJacobian = gfsc_buildJacobianScalar /
!#                          gfsc_buildJacobianBlock
!#      -> Assembles the Jacobian matrix for the group finite element formulation
!#
!# 6.) gfsc_initPerfConfig
!#      -> Initialises the global performance configuration
!#
!# </purpose>
!##############################################################################

module groupfemscalar

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
  
  public :: gfsc_initPerfConfig
  public :: gfsc_buildOperatorNode
  public :: gfsc_buildOperatorEdge
  public :: gfsc_buildVectorNode
  public :: gfsc_buildVectorEdge
  public :: gfsc_buildJacobian
  
!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of nodes to handle simultaneously when building matrices
#ifndef GFSC_NEQSIM
  integer, parameter, public :: GFSC_NEQSIM = 128
#endif

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of edges to handle simultaneously when building matrices
#ifndef GFSC_NEDGESIM
  integer, parameter, public :: GFSC_NEDGESIM = 64
#endif

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of nonzero entries to handle simultaneously when building matrices
#ifndef GFSC_NASIM
  integer, parameter, public :: GFSC_NASIM = 1000
#endif
!</constantblock>

!</constants>

  !************************************************************************
  
  ! global performance configuration
  type(t_perfconfig), target, save :: gfsc_perfconfig
  
  ! ****************************************************************************

  interface gfsc_buildOperatorNode
    module procedure gfsc_buildOperatorConst
    module procedure gfsc_buildOperatorNodeScalar
    module procedure gfsc_buildOperatorNodeBlock1
    module procedure gfsc_buildOperatorNodeBlock2
    module procedure gfsc_buildOperatorNodeBlock3
  end interface

  interface gfsc_buildOperatorEdge
    module procedure gfsc_buildOperatorConst
    module procedure gfsc_buildOperatorEdgeScalar
    module procedure gfsc_buildOperatorEdgeBlock1
    module procedure gfsc_buildOperatorEdgeBlock2
    module procedure gfsc_buildOperatorEdgeBlock3
  end interface
  
  interface gfsc_buildVectorNode
    module procedure gfsc_buildVectorNodeScalar
    module procedure gfsc_buildVectorNodeBlock
  end interface

  interface gfsc_buildVectorEdge
    module procedure gfsc_buildVectorEdgeScalar
    module procedure gfsc_buildVectorEdgeBlock
  end interface

  interface gfsc_buildJacobian
    module procedure gfsc_buildJacobianScalar
    module procedure gfsc_buildJacobianBlock
  end interface

contains
 
  !****************************************************************************

!<subroutine>

  subroutine gfsc_initPerfConfig(rperfconfig)

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
      gfsc_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(gfsc_perfconfig)
      gfsc_perfconfig%NEQSIM   = GFSC_NEQSIM
      gfsc_perfconfig%NEDGESIM = GFSC_NEDGESIM
      gfsc_perfconfig%NASIM    = GFSC_NASIM
    end if
  
  end subroutine gfsc_initPerfConfig

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorConst(rgroupFEMSet, dscale, bclear, rmatrix,&
      cconstrType, rafcstab, rperfconfig)

!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries are taken from
    ! the precomputed values stored in the structure rgroupFEMSet.
    !
    ! This is the simplest of the gfsc_buildOperatorX routines. It can
    ! handle both node-by-node and edge-by-edge matrix assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! for each edge IJ the matrix entries and the artificial diffusion
    ! coefficient is stored according to the following convention:
    !
    ! For upwind-biased stabilisation
    !   Coefficients(1:3, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! For symmetric stabilisation
    !   Coefficients(1:2, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! Otherwise, the standard operator is assembled without stabilisation.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Scaling factor
    real(DP), intent(in) :: dscale
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtDiag
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge

    real(SP), dimension(:), pointer :: p_Fdata
    real(SP), dimension(:,:), pointer :: p_Fcoefficients
    real(SP), dimension(:,:), pointer :: p_FcoeffsAtDiag
    real(SP), dimension(:,:), pointer :: p_FcoeffsAtNode
    real(SP), dimension(:,:,:), pointer :: p_FcoeffsAtEdge

    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:,:), pointer :: p_InodeList
    integer, dimension(:,:), pointer :: p_IdiagList
    integer, dimension(:,:), pointer :: p_InodeListIdx2D
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer, dimension(:), pointer :: p_InodeListIdx1D
    
    integer :: ccType
    logical :: bsymm

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => gfsc_perfconfig
    end if

    ! Check if matrix and vector have the same data type
    if (rmatrix%cdataType .ne. rgroupFEMSet%cdataType) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
      call sys_halt()
    end if

    ! Set type of matrix construction method
    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType

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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
        call sys_halt()
      end if

      if (present(rafcstab)) then
        call output_line('Stabilisation is not feasible in node-by-node assembly!',&
            OU_CLASS_WARNING,OU_MODE_STD,'gfsc_buildOperatorConst')
      end if
      
      ! What data types are we?
      select case(rgroupFEMSet%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
        call lsyssc_getbase_double(rmatrix, p_Ddata)
        
        ! Check if only a subset of the matrix is required
        if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then

          select case(ccType)
          case (GFEM_MATC_CONSISTENT)
            call doOperatorNodeConsistDble(p_DcoeffsAtNode,&
                dscale, bclear, p_Ddata)

          case (GFEM_MATC_LUMPED)
            call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx1D)
            call doOperatorNodeLumpedDble(p_InodeListIdx1D,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata)
          end select

        else ! use restricted DOFs for assembly

          select case(ccType)
          case (GFEM_MATC_CONSISTENT)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOperatorNodeConsistDbleSel(p_InodeList,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata)

          case (GFEM_MATC_LUMPED)
            call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx2D)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOperatorNodeLumpedDbleSel(p_InodeListIdx2D, p_InodeList,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata)
          end select

        end if
        
      case (ST_SINGLE)
        ! Set pointers
        call gfem_getbase_FcoeffsAtNode(rgroupFEMSet, p_FcoeffsAtNode)
        call lsyssc_getbase_single(rmatrix, p_Fdata)
        
        ! Check if only a subset of the matrix is required
        if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then

          select case(ccType)
          case (GFEM_MATC_CONSISTENT)
            call doOperatorNodeConsistSngl(p_FcoeffsAtNode,&
                real(dscale,SP), bclear, p_Fdata)
            
          case (GFEM_MATC_LUMPED)
            call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx1D)
            call doOperatorNodeLumpedSngl(p_InodeListIdx1D,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata)
          end select

        else ! use restricted DOFs for assembly

          select case(ccType)
          case (GFEM_MATC_CONSISTENT)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOperatorNodeConsistSnglSel(p_InodeList,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata)

          case (GFEM_MATC_LUMPED)
            call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx2D)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOperatorNodeLumpedSnglSel(p_InodeListIdx2D, p_InodeList,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata)
          end select
            
        end if
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
        call sys_halt()
      end select

      
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
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
        
        ! Assemble matrix diagonal
        call doOperatorDiagDble(p_IdiagList, p_DcoeffsAtDiag,&
            dscale, bclear, p_Ddata)
                
        ! Do we have to build the stabilisation?
        if (present(rafcstab)) then
          
          ! Check if stabilisation has been prepared
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
            call output_line('Stabilisation has not been prepared!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
            call sys_halt()
          end if
          
          ! Symmetric artificial diffusion?
          bsymm = .not.(rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED)
          
          ! Check if coefficients should be stored in stabilisation
          if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE) then
            
            ! Check if stabilisation has the same data type
            if (rafcstab%cdataType .ne. ST_DOUBLE) then
              call output_line('Stabilisation must have double precision!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
              call sys_halt()
            end if
            
            ! Set additional pointers
            call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
            
            !-------------------------------------------------------------------
            ! Assemble operator with stabilisation and generate coefficients
            !-------------------------------------------------------------------
            call doOperatorEdgeAFCDble(p_IedgeListIdx, p_IedgeList,&
                p_DcoeffsAtEdge, dscale, bclear, bsymm, p_Ddata, p_Dcoefficients)
            
            ! Set state of stabilisation
            rafcstab%istabilisationSpec =&
                ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)
            
            ! Do we need edge orientation?
            if (rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED) then
              call afcstab_upwindOrientation(p_Dcoefficients, p_IedgeList,&
                  p_DcoeffsAtEdge, 2, 3)
              rafcstab%istabilisationSpec =&
                  ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION)
            else
              rafcstab%istabilisationSpec =&
                  iand(rafcstab%istabilisationSpec, not(AFCSTAB_HAS_EDGEORIENTATION))
            end if
            
          else
            
            !-------------------------------------------------------------------
            ! Assemble operator with stabilisation but do not generate coeffs
            !-------------------------------------------------------------------
            call doOperatorEdgeStabDble(p_IedgeListIdx, p_IedgeList,&
                p_DcoeffsAtEdge, dscale, bclear, bsymm, p_Ddata)
          end if
          
        else   ! no stabilisation structure present
          
          !---------------------------------------------------------------------
          ! Assemble operator without stabilisation
          !---------------------------------------------------------------------
          call doOperatorEdgeDble(p_IedgeListIdx, p_IedgeList,&
              p_DcoeffsAtEdge, dscale, bclear, p_Ddata)
        end if
        
      case (ST_SINGLE)
        ! Set pointers
        call gfem_getbase_FcoeffsAtDiag(rgroupFEMSet, p_FcoeffsAtDiag)
        call gfem_getbase_FcoeffsAtEdge(rgroupFEMSet, p_FcoeffsAtEdge)
        call lsyssc_getbase_single(rmatrix, p_Fdata)
        
        ! Assemble matrix diagonal
        call doOperatorDiagSngl(p_IdiagList, p_FcoeffsAtDiag,&
            real(dscale,SP), bclear, p_Fdata)
        
        ! Do we have to build the stabilisation?
        if (present(rafcstab)) then
          
          ! Check if stabilisation has been prepared
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
            call output_line('Stabilisation has not been prepared!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
            call sys_halt()
          end if
          
          ! Symmetric artificial diffusion?
          bsymm = .not.(rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED)
          
          ! Check if coefficients should be stored in stabilisation
          if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE) then
            
            ! Check if stabilisation has the same data type
            if (rafcstab%cdataType .ne. ST_SINGLE) then
              call output_line('Stabilisation must have double precision!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
              call sys_halt()
            end if
            
            ! Set additional pointers
            call afcstab_getbase_FcoeffsAtEdge(rafcstab, p_Fcoefficients)
            
            !-------------------------------------------------------------------
            ! Assemble operator with stabilisation
            !-------------------------------------------------------------------
            call doOperatorEdgeAFCSngl(p_IedgeListIdx, p_IedgeList,&
                p_FcoeffsAtEdge, real(dscale,SP), bclear, bsymm, p_Fdata, p_Fcoefficients)
            
            ! Set state of stabilisation
            rafcstab%istabilisationSpec =&
                ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)
            
            ! Do we need edge orientation?
            if (rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED) then
              call afcstab_upwindOrientation(p_Fcoefficients, p_IedgeList,&
                  p_FcoeffsAtEdge, 2, 3)
              rafcstab%istabilisationSpec =&
                  ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION)
            else
              rafcstab%istabilisationSpec =&
                  iand(rafcstab%istabilisationSpec, not(AFCSTAB_HAS_EDGEORIENTATION))
            end if
            
          else
            
            !-------------------------------------------------------------------
            ! Assemble operator with stabilisation but do not generate coeffs
            !-------------------------------------------------------------------
            call doOperatorEdgeStabSngl(p_IedgeListIdx, p_IedgeList,&
                p_FcoeffsAtEdge, real(dscale,SP), bclear, bsymm, p_Fdata)
          end if
          
        else   ! no stabilisation structure present
          
          !---------------------------------------------------------------------
          ! Assemble operator without stabilisation
          !---------------------------------------------------------------------
          call doOperatorEdgeSngl(p_IedgeListIdx, p_IedgeList,&
              p_FcoeffsAtEdge, real(dscale,SP), bclear, p_Fdata)
        end if
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow

    !**************************************************************
    ! Assemble operator node-by-node in consistent manner
    
    subroutine doOperatorNodeConsistDble(DcoeffsAtNode, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! local variables
      integer :: ia

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------
      
      if (bclear) then
        
        ! Loop over all matrix entries
        !$omp parallel do default(shared)&
        !$omp if(size(Ddata) > p_rperfconfig%NAMIN_OMP)
        do ia = 1, size(Ddata)
          
          ! Update the matrix coefficient
          Ddata(ia) = dscale*DcoeffsAtNode(1,ia)
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over all matrix entries
        !$omp parallel do default(shared)&
        !$omp if(size(Ddata) > p_rperfconfig%NAMIN_OMP)
        do ia = 1, size(Ddata)
          
          ! Update the matrix coefficient
          Ddata(ia) = Ddata(ia) + dscale*DcoeffsAtNode(1,ia)
        end do
        !$omp end parallel do
        
      end if
      
    end subroutine doOperatorNodeConsistDble

    !**************************************************************
    ! Assemble operator node-by-node in lumped manner
    
    subroutine doOperatorNodeLumpedDble(InodeListIdx,&
        DcoeffsAtNode, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      integer, dimension(:), intent(in) :: InodeListIdx
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! local variables
      real(DP) :: dtemp
      integer :: ieq,ia

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------
      
      if (bclear) then
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(dtemp,ia)&
        !$omp if(size(Ddata) > p_rperfconfig%NAMIN_OMP)
        do ieq = 1, size(InodeListIdx)-1
          
          ! Clear temporal data
          dtemp = 0.0_DP

          ! Loop over all matrix enties in current row
          do ia = InodeListIdx(ieq), InodeListIdx(ieq+1)-1
            
            ! Update the matrix coefficient
            dtemp = dtemp + DcoeffsAtNode(1,ia)
          end do
          
          ! Update the diagonal entry of the global operator
          Ddata(InodeListIdx(ieq)) = dscale*dtemp
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(dtemp,ia)&
        !$omp if(size(Ddata) > p_rperfconfig%NAMIN_OMP)
        do ieq = 1, size(InodeListIdx)-1
          
          ! Clear temporal data
          dtemp = 0.0_DP

          ! Loop over all matrix enties in current row
          do ia = InodeListIdx(ieq), InodeListIdx(ieq+1)-1
            
            ! Update the matrix coefficient
            dtemp = dtemp + DcoeffsAtNode(1,ia)
          end do
          
          ! Update the diagonal entry of the global operator
          Ddata(InodeListIdx(ieq)) = Ddata(InodeListIdx(ieq)) + dscale*dtemp
        end do
        !$omp end parallel do
        
      end if
      
    end subroutine doOperatorNodeLumpedDble

    !**************************************************************
    ! Assemble operator node-by-node in consistent manner
    
    subroutine doOperatorNodeConsistDbleSel(InodeList,&
        DcoeffsAtNode, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: InodeList

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! local variables
      integer :: ia,idx

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------

      if (bclear) then
        
        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(ia)&
        !$omp if (size(InodeList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(InodeList,2)
          
          ! Get position of matrix entry
          ia = InodeList(2,idx)
          
          ! Update the matrix coefficient
          Ddata(ia) = dscale*DcoeffsAtNode(1,idx)
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(ia)&
        !$omp if (size(InodeList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(InodeList,2)
          
          ! Get position of matrix entry
          ia = InodeList(2,idx)
          
          ! Update the matrix coefficient
          Ddata(ia) = Ddata(ia) + dscale*DcoeffsAtNode(1,idx)
        end do
        !$omp end parallel do
        
      end if

    end subroutine doOperatorNodeConsistDbleSel

    !**************************************************************
    ! Assemble operator node-by-node in lumped manner
    
    subroutine doOperatorNodeLumpedDbleSel(InodeListIdx, InodeList,&
        DcoeffsAtNode, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: InodeListIdx,InodeList

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! local variables
      real(DP) :: dtemp
      integer :: ia,idx

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------
      
      if (bclear) then
        
        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(dtemp,ia)&
        !$omp if (size(InodeList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(InodeListIdx,2)
          
          ! Clear temporal data
          dtemp = 0.0_DP

          ! Loop over all matrix enties in current row
          do ia = InodeListIdx(2,idx), InodeListIdx(2,idx+1)-1

            ! Update the matrix coefficient
            dtemp = dtemp + DcoeffsAtNode(1,idx)
          end do

          ! Update the diagonal entry of the global operator
          Ddata(InodeListIdx(2,idx)) = dscale*dtemp
        end do
        !$omp end parallel do
        
      else

        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(dtemp,ia)&
        !$omp if (size(InodeList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(InodeListIdx,2)
          
          ! Clear temporal data
          dtemp = 0.0_DP

          ! Loop over all matrix enties in current row
          do ia = InodeListIdx(2,idx), InodeListIdx(2,idx+1)-1

            ! Update the matrix coefficient
            dtemp = dtemp + DcoeffsAtNode(1,idx)
          end do

          ! Update the diagonal entry of the global operator
          Ddata(InodeListIdx(2,idx)) = Ddata(InodeListIdx(2,idx)) + dscale*dtemp
        end do
        !$omp end parallel do

      end if
      
    end subroutine doOperatorNodeLumpedDbleSel

    !**************************************************************
    ! Assemble operator node-by-node in consistent manner

    subroutine doOperatorNodeConsistSngl(FcoeffsAtNode, fscale, bclear, Fdata)

      ! input parameters
      real(SP), dimension(:,:), intent(in) :: FcoeffsAtNode
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! local variables
      integer :: ia

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------
      
      if (bclear) then
        
        ! Loop over all matrix entries
        !$omp parallel do default(shared)&
        !$omp if (size(Fdata) > p_rperfconfig%NAMIN_OMP)
        do ia = 1, size(Fdata)
          
          ! Update the matrix coefficient
          Fdata(ia) = fscale*FcoeffsAtNode(1,ia)
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over all matrix entries
        !$omp parallel do default(shared)&
        !$omp if (size(Fdata) > p_rperfconfig%NAMIN_OMP)
        do ia = 1, size(Fdata)
          
          ! Update the matrix coefficient
          Fdata(ia) = Fdata(ia) + fscale*FcoeffsAtNode(1,ia)
        end do
        !$omp end parallel do
        
      end if
      
    end subroutine doOperatorNodeConsistSngl

    !**************************************************************
    ! Assemble operator node-by-node in lumped manner
    
    subroutine doOperatorNodeLumpedSngl(InodeListIdx,&
        FcoeffsAtNode, fscale, bclear, Fdata)

      ! input parameters
      real(SP), dimension(:,:), intent(in) :: FcoeffsAtNode
      integer, dimension(:), intent(in) :: InodeListIdx
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! local variables
      real(SP) :: ftemp
      integer :: ieq,ia

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------
      
      if (bclear) then
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(ftemp,ia)&
        !$omp if(size(Fdata) > p_rperfconfig%NAMIN_OMP)
        do ieq = 1, size(InodeListIdx)-1
          
          ! Clear temporal data
          ftemp = 0.0_SP

          ! Loop over all matrix enties in current row
          do ia = InodeListIdx(ieq), InodeListIdx(ieq+1)-1
            
            ! Update the matrix coefficient
            ftemp = ftemp + FcoeffsAtNode(1,ia)
          end do
          
          ! Update the diagonal entry of the global operator
          Fdata(InodeListIdx(ieq)) = fscale*ftemp
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(ftemp,ia)&
        !$omp if(size(Fdata) > p_rperfconfig%NAMIN_OMP)
        do ieq = 1, size(InodeListIdx)-1
          
          ! Clear temporal data
          ftemp = 0.0_SP

          ! Loop over all matrix enties in current row
          do ia = InodeListIdx(ieq), InodeListIdx(ieq+1)-1
            
            ! Update the matrix coefficient
            ftemp = ftemp + FcoeffsAtNode(1,ia)
          end do
          
          ! Update the diagonal entry of the global operator
          Fdata(InodeListIdx(ieq)) = Fdata(InodeListIdx(ieq)) + fscale*ftemp
        end do
        !$omp end parallel do
        
      end if
      
    end subroutine doOperatorNodeLumpedSngl

    !**************************************************************
    ! Assemble operator node-by-node in consistent manner

    subroutine doOperatorNodeConsistSnglSel(InodeList,&
        FcoeffsAtNode, fscale, bclear, Fdata)

      ! input parameters
      real(SP), dimension(:,:), intent(in) :: FcoeffsAtNode
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: InodeList

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! local variables
      integer :: ia,idx

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------

      if (bclear) then
        
        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(ia)&
        !$omp if (size(InodeList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(InodeList,2)
          
          ! Get position of matrix entry
          ia = InodeList(2,idx)
          
          ! Update the matrix coefficient
          Fdata(ia) = fscale*FcoeffsAtNode(1,idx)
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(ia)&
        !$omp if (size(InodeList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(InodeList,2)
          
          ! Get position of matrix entry
          ia = InodeList(2,idx)
          
          ! Update the matrix coefficient
          Fdata(ia) = Fdata(ia) + fscale*FcoeffsAtNode(1,idx)
        end do
        !$omp end parallel do
        
      end if
      
    end subroutine doOperatorNodeConsistSnglSel

    !**************************************************************
    ! Assemble operator node-by-node in lumped manner
    
    subroutine doOperatorNodeLumpedSnglSel(InodeListIdx, InodeList,&
        FcoeffsAtNode, fscale, bclear, Fdata)

      ! input parameters
      real(SP), dimension(:,:), intent(in) :: FcoeffsAtNode
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: InodeListIdx,InodeList

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! local variables
      real(SP) :: ftemp
      integer :: ia,idx

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------
      
      if (bclear) then
        
        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(ftemp,ia)&
        !$omp if (size(InodeList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(InodeListIdx,2)
          
          ! Clear temporal data
          ftemp = 0.0_SP

          ! Loop over all matrix enties in current row
          do ia = InodeListIdx(2,idx), InodeListIdx(2,idx+1)-1

            ! Update the matrix coefficient
            ftemp = ftemp + FcoeffsAtNode(1,idx)
          end do

          ! Update the diagonal entry of the global operator
          Fdata(InodeListIdx(2,idx)) = fscale*ftemp
        end do
        !$omp end parallel do
        
      else

        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(ftemp,ia)&
        !$omp if (size(InodeList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(InodeListIdx,2)
          
          ! Clear temporal data
          ftemp = 0.0_SP

          ! Loop over all matrix enties in current row
          do ia = InodeListIdx(2,idx), InodeListIdx(2,idx+1)-1

            ! Update the matrix coefficient
            ftemp = ftemp + FcoeffsAtNode(1,idx)
          end do

          ! Update the diagonal entry of the global operator
          Fdata(InodeListIdx(2,idx)) = Fdata(InodeListIdx(2,idx)) + fscale*ftemp
        end do
        !$omp end parallel do

      end if
      
    end subroutine doOperatorNodeLumpedSnglSel

    !**************************************************************
    ! Assemble diagonal part of the operator

    subroutine doOperatorDiagDble(IdiagList, DcoeffsAtDiag,&
        dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtDiag
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IdiagList

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! local variables
      integer :: ia,idx

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------
      
      if (bclear) then
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(ia)&
        !$omp if (size(IdiagList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(IdiagList,2)
          
          ! Get position of diagonal entry
          ia = IdiagList(2,idx)
          
          ! Update the diagonal coefficient
          Ddata(ia) = dscale*DcoeffsAtDiag(1,idx)
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(ia)&
        !$omp if (size(IdiagList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(IdiagList,2)
          
          ! Get position of diagonal entry
          ia = IdiagList(2,idx)
          
          ! Update the diagonal coefficient
          Ddata(ia) = Ddata(ia) + dscale*DcoeffsAtDiag(1,idx)
        end do
        !$omp end parallel do
        
      end if
      
    end subroutine doOperatorDiagDble
    
    !**************************************************************
    ! Assemble diagonal part of the operator

    subroutine doOperatorDiagSngl(IdiagList, FcoeffsAtDiag,&
        fscale, bclear, Fdata)

      ! input parameters
      real(SP), dimension(:,:), intent(in) :: FcoeffsAtDiag
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IdiagList

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! local variables
      integer :: ia,idx

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(ia)&
        !$omp if (size(IdiagList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(IdiagList,2)
          
          ! Get position of diagonal entry
          ia = IdiagList(2,idx)
          
          ! Update the diagonal coefficient
          Fdata(ia) = fscale*FcoeffsAtDiag(1,idx)
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(ia)&
        !$omp if (size(IdiagList,2) > p_rperfconfig%NEQMIN_OMP)
        do idx = 1, size(IdiagList,2)
          
          ! Get position of diagonal entry
          ia = IdiagList(2,idx)
          
          ! Update the diagonal coefficient
          Fdata(ia) = Fdata(ia) + fscale*FcoeffsAtDiag(1,idx)
        end do
        !$omp end parallel do
        
      end if
      
    end subroutine doOperatorDiagSngl
        
    !**************************************************************
    ! Assemble operator edge-by-edge without stabilisation
    
    subroutine doOperatorEdgeDble(IedgeListIdx, IedgeList,&
        DcoeffsAtEdge, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata

      ! local variables
      integer :: iedge,igroup,ij,ji

      !$omp parallel default(shared) private(ij,ji)&
      !$omp if (size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)
 
      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then
        
        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronize memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
        
            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            
            ! Update the global operator
            Ddata(ij) = dscale*DcoeffsAtEdge(1,1,iedge)
            Ddata(ji) = dscale*DcoeffsAtEdge(1,2,iedge)
          end do
          !$omp end do
        end do ! igroup

      else   ! do not clear matrix

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronize memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
            
            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            
            ! Update the global operator
            Ddata(ij) = Ddata(ij) + dscale*DcoeffsAtEdge(1,1,iedge)
            Ddata(ji) = Ddata(ji) + dscale*DcoeffsAtEdge(1,2,iedge)
          end do
          !$omp end do
        end do ! igroup

      end if
      !$omp end parallel

    end subroutine doOperatorEdgeDble

    !**************************************************************
    ! Assemble operator edge-by-edge without stabilisation
    
    subroutine doOperatorEdgeSngl(IedgeListIdx, IedgeList,&
        FcoeffsAtEdge, fscale, bclear, Fdata)
      
      ! input parameters
      real(SP), dimension(:,:,:), intent(in) :: FcoeffsAtEdge
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata

      ! local variables
      integer :: iedge,igroup,ij,ji
      
      !$omp parallel default(shared) private(ij,ji)&
      !$omp if (size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then
        
        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronize memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
        
            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            
            ! Update the global operator
            Fdata(ij) = fscale*FcoeffsAtEdge(1,1,iedge)
            Fdata(ji) = fscale*FcoeffsAtEdge(1,2,iedge)
          end do
          !$omp end do
        end do ! igroup

      else   ! do not clear matrix

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronize memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
            
            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            
            ! Update the global operator
            Fdata(ij) = Fdata(ij) + fscale*FcoeffsAtEdge(1,1,iedge)
            Fdata(ji) = Fdata(ji) + fscale*FcoeffsAtEdge(1,2,iedge)
          end do
          !$omp end do
        end do ! igroup

      end if
      !$omp end parallel

    end subroutine doOperatorEdgeSngl

    !**************************************************************
    ! Assemble edge-by-edge operator with stabilisation
    
    subroutine doOperatorEdgeStabDble(IedgeListIdx, IedgeList,&
        DcoeffsAtEdge, dscale, bclear, bsymm, Ddata)
      
      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear, bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata

      ! local variables
      real(DP) :: d_ij
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij)&
      !$omp if (size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then
   
        if (bsymm) then

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Symmetric artificial diffusion coefficient
              d_ij = dscale*max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) =             dscale*DcoeffsAtEdge(1,1,iedge) + d_ij
              Ddata(ji) =             dscale*DcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Non-symmetric artificial diffusion coefficient
              d_ij = dscale*max(0.0_DP,&
                  -DcoeffsAtEdge(1,1,iedge), -DcoeffsAtEdge(1,2,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) =             dscale*DcoeffsAtEdge(1,1,iedge) + d_ij
              Ddata(ji) =             dscale*DcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup

        end if
        
      else   ! do not clear matrix
        
        if (bsymm) then

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Symmetric artificial diffusion coefficient
              d_ij = dscale*max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*DcoeffsAtEdge(1,1,iedge) + d_ij
              Ddata(ji) = Ddata(ji) + dscale*DcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup

        else
          
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Non-symmetric artificial diffusion coefficient
              d_ij = dscale*max(0.0_DP,&
                  -DcoeffsAtEdge(1,1,iedge), -DcoeffsAtEdge(1,2,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*DcoeffsAtEdge(1,1,iedge) + d_ij
              Ddata(ji) = Ddata(ji) + dscale*DcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup

        end if

      end if
      !$omp end parallel

    end subroutine doOperatorEdgeStabDble

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation
    
    subroutine doOperatorEdgeStabSngl(IedgeListIdx, IedgeList,&
        FcoeffsAtEdge, fscale, bclear, bsymm, Fdata)

      ! input parameters
      real(SP), dimension(:,:,:), intent(in) :: FcoeffsAtEdge
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear, bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata

      ! local variables
      real(SP) :: d_ij
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij)&
      !$omp if (size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then
   
        if (bsymm) then

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Symmetric artificial diffusion coefficient
              d_ij = fscale*max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) =             fscale*FcoeffsAtEdge(1,1,iedge) + d_ij
              Fdata(ji) =             fscale*FcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup

        else
          
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Non-symmetric artificial diffusion coefficient
              d_ij = fscale*max(0.0_SP,&
                  -FcoeffsAtEdge(1,1,iedge), -FcoeffsAtEdge(1,2,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) =             fscale*FcoeffsAtEdge(1,1,iedge) + d_ij
              Fdata(ji) =             fscale*FcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup
          
        end if

      else   ! do not clear matrix
        
        if (bsymm) then

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Symmetric atificial diffusion coefficient
              d_ij = fscale*max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*FcoeffsAtEdge(1,1,iedge) + d_ij
              Fdata(ji) = Fdata(ji) + fscale*FcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Non-symmetric atificial diffusion coefficient
              d_ij = fscale*max(0.0_SP,&
                  -FcoeffsAtEdge(1,1,iedge), -FcoeffsAtEdge(1,2,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*FcoeffsAtEdge(1,1,iedge) + d_ij
              Fdata(ji) = Fdata(ji) + fscale*FcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup
          
        end if

      end if
      !$omp end parallel

    end subroutine doOperatorEdgeStabSngl

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation and AFC data.
    
    subroutine doOperatorEdgeAFCDble(IedgeListIdx, IedgeList,&
        DcoeffsAtEdge, dscale, bclear, bsymm, Ddata, Dcoefficients)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear,bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! output parameters
      real(DP), dimension(:,:), intent(out) :: Dcoefficients
      
      ! local variables
      real(DP) :: d_ij,k_ij,k_ji
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij,k_ij,k_ji)&
      !$omp if (size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then

        if (bsymm) then
        
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Artificial diffusion coefficient
              d_ij = dscale*max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              k_ij = dscale*max(0.0_DP,  DcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Dcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) =             dscale*DcoeffsAtEdge(1,1,iedge) + d_ij
              Ddata(ji) =             dscale*DcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Artificial diffusion coefficient
              d_ij = dscale*max(0.0_DP,&
                  -DcoeffsAtEdge(1,1,iedge), -DcoeffsAtEdge(1,2,iedge))
              k_ij = dscale*DcoeffsAtEdge(1,1,iedge) + d_ij
              k_ji = dscale*DcoeffsAtEdge(1,2,iedge) + d_ij
              
              ! Non-symmetric AFC w/o edge orientation
              Dcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) =             k_ij
              Ddata(ji) =             k_ji
            end do
            !$omp end do
          end do ! igroup

        end if

      else   ! do not clear matrix
        
        if (bsymm) then
          
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Artificial diffusion coefficient
              d_ij = dscale*max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              k_ij = dscale*max(0.0_DP,  DcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Dcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*DcoeffsAtEdge(1,1,iedge) + d_ij
              Ddata(ji) = Ddata(ji) + dscale*DcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Artificial diffusion coefficient
              d_ij = dscale*max(0.0_DP,&
                  -DcoeffsAtEdge(1,1,iedge), -DcoeffsAtEdge(1,2,iedge))
              k_ij = dscale*DcoeffsAtEdge(1,1,iedge) + d_ij
              k_ji = dscale*DcoeffsAtEdge(1,2,iedge) + d_ij
              
              ! Non-symmetric AFC w/o edge orientation
              Dcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + k_ij
              Ddata(ji) = Ddata(ji) + k_ji
            end do
            !$omp end do
          end do ! igroup

        end if
        
      end if
      !$omp end parallel

    end subroutine doOperatorEdgeAFCDble

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation and AFC data.
    
    subroutine doOperatorEdgeAFCSngl(IedgeListIdx, IedgeList,&
        FcoeffsAtEdge, fscale, bclear, bsymm, Fdata, Fcoefficients)

      ! input parameters
      real(SP), dimension(:,:,:), intent(in) :: FcoeffsAtEdge
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear,bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! output parameters
      real(SP), dimension(:,:), intent(out) :: Fcoefficients
      
      ! local variables
      real(SP) :: d_ij,k_ij,k_ji
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij,k_ij,k_ji)&
      !$omp if (size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then

        if (bsymm) then
        
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Artificial diffusion coefficient
              d_ij = fscale*max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              k_ij = fscale*max(0.0_SP,  FcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Fcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) =             fscale*FcoeffsAtEdge(1,1,iedge) + d_ij
              Fdata(ji) =             fscale*FcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Artificial diffusion coefficient
              d_ij = fscale*max(0.0_SP,&
                  -FcoeffsAtEdge(1,1,iedge), -FcoeffsAtEdge(1,2,iedge))
              k_ij = fscale*FcoeffsAtEdge(1,1,iedge) + d_ij
              k_ji = fscale*FcoeffsAtEdge(1,2,iedge) + d_ij
              
              ! Non-symmetric AFC w/o edge orientation
              Fcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) =             k_ij
              Fdata(ji) =             k_ji
            end do
            !$omp end do
          end do ! igroup

        end if

      else   ! do not clear matrix
        
        if (bsymm) then
          
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Artificial diffusion coefficient
              d_ij = fscale*max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              k_ij = fscale*max(0.0_SP,  FcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Fcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*FcoeffsAtEdge(1,1,iedge) + d_ij
              Fdata(ji) = Fdata(ji) + fscale*FcoeffsAtEdge(1,2,iedge) + d_ij
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Artificial diffusion coefficient
              d_ij = fscale*max(0.0_SP,&
                  -FcoeffsAtEdge(1,1,iedge), -FcoeffsAtEdge(1,2,iedge))
              k_ij = fscale*FcoeffsAtEdge(1,1,iedge) + d_ij
              k_ji = fscale*FcoeffsAtEdge(1,2,iedge) + d_ij
              
              ! Non-symmetric AFC w/o edge orientation
              Fcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + k_ij
              Fdata(ji) = Fdata(ji) + k_ji
            end do
            !$omp end do
          end do ! igroup

        end if
        
      end if
      !$omp end parallel

    end subroutine doOperatorEdgeAFCSngl

  end subroutine gfsc_buildOperatorConst

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorNodeBlock1(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix, cconstrType,&
      rcollection, rafcstab, fcb_calcOperatorNodeSc, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! Note that this routine serves as a wrapper for block vectors and
    ! matrices. If there is only one block, then the corresponding
    ! scalar routine is called. Otherwise, an error is thrown.
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
    include 'intf_calcMatrixDiagSc_sim.inc'

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorNodeSc.inc'
    optional :: fcb_calcOperatorNodeSc

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
    
    ! Check if user-defined callback function is present
    if (present(fcb_calcOperatorNodeSc)) then
      
      call fcb_calcOperatorNodeSc(rgroupFEMSet, rx, rmatrix,&
          dscale, bclear, cconstrType, fcb_calcMatrixDiagSc_sim,&
          rcollection, rafcstab)
      
    elseif ((rx%nblocks            .eq. 1) .and.&
            (rmatrix%nblocksPerCol .eq. 1) .and.&
            (rmatrix%nblocksPerRow .eq. 1)) then

      ! Call scalar version of this routine
      call gfsc_buildOperatorNodeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix%RmatrixBlock(1,1),&
          cconstrType, rcollection, rafcstab, rperfconfig=rperfconfig)
      
    else
      call output_line('Matrix/Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeBlock1')
      call sys_halt()
    end if
    
  end subroutine gfsc_buildOperatorNodeBlock1

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorNodeBlock2(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, dscale, bclear, cconstrType, rmatrix,&
      rcollection, rafcstab, fcb_calcOperatorNodeSc, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! Note that this routine serves as a wrapper for block matrices. If
    ! there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, an error is thrown.
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
    include 'intf_calcMatrixDiagSc_sim.inc'

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorNodeSc.inc'
    optional :: fcb_calcOperatorNodeSc

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
    if (present(fcb_calcOperatorNodeSc)) then
      
      ! Create auxiliary 1-block vector
      call lsysbl_createVecFromScalar(rx, rxBlock)

      call fcb_calcOperatorNodeSc(rgroupFEMSet, rxBlock, rmatrix,&
          dscale, bclear, cconstrType, fcb_calcMatrixDiagSc_sim,&
          rcollection, rafcstab)
      
      ! Release auxiliary 1-block vector
      call lsysbl_releaseVector(rxBlock)

    elseif ((rmatrix%nblocksPerCol .eq. 1) .and.&
            (rmatrix%nblocksPerRow .eq. 1)) then

      ! Call scalar version of this routine
      call gfsc_buildOperatorNodeScalar(rgroupFEMSet, rx,&
          fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix%RmatrixBlock(1,1),&
          cconstrType, rcollection, rafcstab, rperfconfig=rperfconfig)
      
    else
      call output_line('Matrix must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeBlock2')
      call sys_halt()
    end if
    
  end subroutine gfsc_buildOperatorNodeBlock2

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorNodeBlock3(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix, cconstrType,&
      rcollection, rafcstab, fcb_calcOperatorNodeSc, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! Note that this routine serves as a wrapper for block vectors. If
    ! there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, an error is thrown.
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
    include 'intf_calcMatrixDiagSc_sim.inc'

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorNodeSc.inc'
    optional :: fcb_calcOperatorNodeSc

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
    if (present(fcb_calcOperatorNodeSc)) then
      
      ! Create auxiliary 1-block matrix
      call lsysbl_createMatFromScalar(rmatrix, rmatrixBlock)

      call fcb_calcOperatorNodeSc(rgroupFEMSet, rx, rmatrixBlock,&
          dscale, bclear, cconstrType, fcb_calcMatrixDiagSc_sim,&
          rcollection, rafcstab)
      
      ! Release auxiliary 1-block matrix
      call lsysbl_releaseMatrix(rmatrixBlock)

    elseif (rx%nblocks .eq. 1) then

      ! Call scalar version of this routine
      call gfsc_buildOperatorNodeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix,&
          cconstrType, rcollection, rafcstab, rperfconfig=rperfconfig)
      
    else
      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeBlock3')
      call sys_halt()
    end if
    
  end subroutine gfsc_buildOperatorNodeBlock3

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorNodeScalar(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix, cconstrType,&
      rcollection, rafcstab, fcb_calcOperatorNodeSc, rperfconfig)
    
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
    include 'intf_calcMatrixDiagSc_sim.inc'
    
    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorNodeSc.inc'
    optional :: fcb_calcOperatorNodeSc

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
    integer, dimension(:,:), pointer :: p_InodeListIdx2D,p_InodeList2D
    integer, dimension(:), pointer :: p_InodeListIdx1D,p_InodeList1D
    integer :: ccType

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    ! Check if user-defined assembly is provided
    if (present(fcb_calcOperatorNodeSc)) then
      ! Create auxiliary 1-block vector and matrix
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createMatFromScalar(rmatrix, rmatrixBlock)
      
      ! Call user-defined assembly
      call fcb_calcOperatorNodeSc(rgroupFEMSet, rxBlock, rmatrixBlock,&
          dscale, bclear, cconstrType, fcb_calcMatrixDiagSc_sim,&
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
      p_rperfconfig => gfsc_perfconfig
    end if

    ! Check if matrix and vector have the same data type
    if ((rmatrix%cdataType .ne. rx%cdataType) .or.&
        (rmatrix%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
      call sys_halt()
    end if

    ! Set type of matrix construction method
    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType

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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
        call sys_halt()
      end if

      ! What data types are we?
      select case(rgroupFEMSet%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
        call lsyssc_getbase_double(rmatrix, p_Ddata)
        call lsyssc_getbase_double(rx, p_Dx)
        
        ! Check if only a subset of the matrix is required
        if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DOFLIST) .eq. 0) then
          
          select case(ccType)
          case (GFEM_MATC_CONSISTENT)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList1D)
            call doOperatorConsistDble(p_InodeList1D,&
                p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)
            
          case (GFEM_MATC_LUMPED)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList1D)
            call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx1D)
            call doOperatorLumpedDble(p_InodeListIdx1D, p_InodeList1D,&
                p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)
          end select

        else ! use restricted DOFs for assembly

          select case(ccType)
          case (GFEM_MATC_CONSISTENT)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList2D)
            call doOperatorConsistDbleSel(p_InodeList2D,&
                p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)

          case (GFEM_MATC_LUMPED)
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList2D)
            call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx2D)
            call doOperatorLumpedDbleSel(p_InodeListIdx2D, p_InodeList2D,&
                p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)
          end select

        end if

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow

    !**************************************************************
    ! Assemble operator node-by-node in consistent manner

    subroutine doOperatorConsistDble(InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: InodeList

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:), pointer :: DdataAtNode
      real(DP), dimension(:,:), pointer :: Dcoefficients
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
      
      ! Allocate temporal memory
      allocate(DdataAtNode(p_rperfconfig%NASIM))
      allocate(Dcoefficients(1,p_rperfconfig%NASIM))
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
          DdataAtNode(idx)   = Dx(j)
        end do
        
        ! Use callback function to compute matrix entries
        call fcb_calcMatrixDiagSc_sim(&
            DdataAtNode(1:IAmax-IAset+1),&
            DcoeffsAtNode(:,IAset:IAmax),&
            IdofsAtNode(:,1:IAmax-IAset+1),&
            dscale, IAmax-IAset+1,&
            Dcoefficients(:,1:IAmax-IAset+1), rcollection)
        
        ! Loop through all nonzero matrix entries in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then

          do idx = 1, IAmax-IAset+1
            
            ! Get position of matrix entry
            ij = idx+IAset-1
            
            ! Update the global operator
            Ddata(ij) = Dcoefficients(1,idx)
          end do
          
        else   ! do not clear matrix
          
          do idx = 1, IAmax-IAset+1
            
            ! Get position of matrix entry
            ij = idx+IAset-1
            
            ! Update the global operator
            Ddata(ij) = Ddata(ij) + Dcoefficients(1,idx)
          end do
        end if
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(IdofsAtNode)
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorConsistDble

    !**************************************************************
    ! Assemble operator node-by-node in lumped manner

    subroutine doOperatorLumpedDble(InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: InodeListIdx, InodeList

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:), pointer :: DdataAtNode
      real(DP), dimension(:,:), pointer :: Dcoefficients
      integer, dimension(:,:), pointer :: IdofsAtNode
      
      ! local variables
      real(DP) :: dtemp
      integer :: IAmax,IApos,IAset,IEQmax,IEQset,ia,idx,ieq
      
      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------

      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,dtemp,IdofsAtNode,&
      !$omp         IAmax,IApos,IAset,IEQmax,ia,idx,ieq)&
      !$omp if(size(InodeList) > p_rperfconfig%NAMIN_OMP)
      
      ! Allocate temporal memory
      allocate(DdataAtNode(p_rperfconfig%NASIM))
      allocate(Dcoefficients(1,p_rperfconfig%NASIM))
      allocate(IdofsAtNode(2,p_rperfconfig%NASIM))

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
              IdofsAtNode(1,idx) = InodeList(ia)     ! absolut nodal value j
              IdofsAtNode(2,idx) = ia                ! absolute matrix position ia
              DdataAtNode(idx)   = Dx(InodeList(ia)) ! solution value at node j
            end do
            
            ! Increase the upper bound for nonzero entries
            IAmax = IAmax+1
          end do
          
          ! Use callback function to compute matrix entries
          call fcb_calcMatrixDiagSc_sim(&
              DdataAtNode(1:idx),&
              DcoeffsAtNode(:,IApos:IApos+idx-1),&
              IdofsAtNode(:,1:idx),&
              dscale, idx,&
              Dcoefficients(:,1:idx), rcollection)

          ! Initialise local index which will run from IAset..IAmax
          idx = 0

          ! Loop through all nonzero matrix entries in the current set
          ! and scatter the entries to the diagonal of the global matrix
          do ieq = IAset, IAmax-1
            
            ! Clear temporal data
            dtemp = 0.0_DP
            
            ! Loop over all contributions to this equation
            do ia = InodeListIdx(ieq), InodeListIdx(ieq+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Update temporal data
              dtemp = dtemp + Dcoefficients(1,idx)
            end do
            
            ! Update the diagonal entry of the global operator
            Ddata(InodeListIdx(ieq)) = Ddata(InodeListIdx(ieq))+dtemp
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
            
    end subroutine doOperatorLumpedDble

    !**************************************************************
    ! Assemble operator node-by-node in consistent manner

    subroutine doOperatorConsistDbleSel(InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: InodeList

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:), pointer :: DdataAtNode
      real(DP), dimension(:,:), pointer :: Dcoefficients
      
      ! local variables
      integer :: idx,IAset,IAmax
      integer :: ij,j

      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,IAmax,idx,ij,j)&
      !$omp if (size(InodeList,2) > p_rperfconfig%NAMIN_OMP)
      
      ! Allocate temporal memory
      allocate(DdataAtNode(p_rperfconfig%NASIM))
      allocate(Dcoefficients(1,p_rperfconfig%NASIM))
      
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
          DdataAtNode(idx) = Dx(InodeList(1,idx+IAset-1))
        end do
        
        ! Use callback function to compute matrix entries
        call fcb_calcMatrixDiagSc_sim(&
            DdataAtNode(1:IAmax-IAset+1),&
            DcoeffsAtNode(:,IAset:IAmax),&
            InodeList(:,IAset:IAmax),&
            dscale, IAmax-IAset+1,&
            Dcoefficients(:,1:IAmax-IAset+1), rcollection)
        
        ! Loop through all nonzero matrix entries in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then

          do idx = 1, IAmax-IAset+1
            
            ! Get position of matrix entry
            ij = InodeList(2,idx+IAset-1)
            
            ! Update the global operator
            Ddata(ij) = Dcoefficients(1,idx)
          end do
          
        else   ! do not clear matrix
          
          do idx = 1, IAmax-IAset+1
            
            ! Get position of matrix entry
            ij = InodeList(2,idx+IAset-1)

            ! Update the global operator
            Ddata(ij) = Ddata(ij) + Dcoefficients(1,idx)
          end do
        end if
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorConsistDbleSel

    !**************************************************************
    ! Assemble operator node-by-node in lumped manner

    subroutine doOperatorLumpedDbleSel(InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: InodeListIdx,InodeList

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:), pointer :: DdataAtNode
      real(DP), dimension(:,:), pointer :: Dcoefficients
      
      ! local variables
      real(DP) :: dtemp
      integer :: IAmax,IApos,IAset,IEQmax,IEQset,ia,idx,iidx

      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------

      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,dtemp,&
      !$omp         IAmax,IApos,IAset,IEQmax,ia,idx,iidx)&
      !$omp if(size(InodeList,2) > p_rperfconfig%NAMIN_OMP)
      
      ! Allocate temporal memory
      allocate(DdataAtNode(p_rperfconfig%NASIM))
      allocate(Dcoefficients(1,p_rperfconfig%NASIM))
      
      ! Loop over all equations in blocks of size NEQSIM
      !$omp do schedule(static,1)
      do IEQset = 1, size(InodeListIdx,2)-1, p_rperfconfig%NEQSIM
        
        ! We always handle NEQSIM equations by one OpenMP thread.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle
        ! at most NEQSIM equations by one OpenMP thread.

        IEQmax = min(size(InodeListIdx,2)-1, IEQset-1+p_rperfconfig%NEQSIM)
        
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
        IApos = InodeListIdx(1,IEQset)
        
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
            if (InodeListIdx(1,IAmax+1)-IApos .gt. p_rperfconfig%NASIM) exit
            
            ! Loop through all nonzero matrix entries in the current
            ! equation and prepare the auxiliary arrays for it
            do ia = InodeListIdx(1,IAmax), InodeListIdx(1,IAmax+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Fill auxiliary arrays
              DdataAtNode(idx) = Dx(InodeList(1,ia)) ! solution value at node j
            end do
            
            ! Increase the upper bound for nonzero entries
            IAmax = IAmax+1
          end do
        
          ! Use callback function to compute matrix entries
          call fcb_calcMatrixDiagSc_sim(&
              DdataAtNode(1:idx),&
              DcoeffsAtNode(:,IApos:IApos+idx-1),&
              InodeList(:,IApos:IApos+idx-1),&
              dscale, idx,&
              Dcoefficients(:,1:idx), rcollection)

          ! Initialise local index which will run from IAset..IAmax
          idx = 0
          
          ! Loop through all equations in the current set
          ! and scatter the entries to the global matrix
          do iidx = IAset, IAmax-1
                        
            ! Clear temporal date
            dtemp = 0.0_DP

            ! Loop over all contributions to this equation
            do ia = InodeListIdx(1,iidx), InodeListIdx(1,iidx+1)-1
              
              ! Update local index
              idx = idx+1

              ! Update temporal data
              dtemp = dtemp + Dcoefficients(1,idx)
            end do
            
            ! Update the diagonal entry of the global operator
            Ddata(InodeListIdx(3,iidx)) = Ddata(InodeListIdx(3,iidx))+dtemp
          end do
          
          ! Proceed with next nonzero entries in current set
          IAset = IAmax
          IApos = InodeListIdx(1,IASet)
          
        end do
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel

    end subroutine doOperatorLumpedDbleSel

  end subroutine gfsc_buildOperatorNodeScalar
  
  !*****************************************************************************

!<subroutine>


  subroutine gfsc_buildOperatorEdgeBlock1(rgroupFEMSet, rx,&
      fcb_calcMatrixSc_sim, dscale, bclear, rmatrix,&
      fcb_calcMatrixDiagSc_sim, cconstrType, rcollection,&
      rafcstab, fcb_calcOperatorEdgeSc, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only edge-by-edge assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! for each edge IJ the matrix entries and the artificial diffusion
    ! coefficient is stored according to the following convention:
    !
    ! For upwind-biased stabilisation
    !   Coefficients(1:3, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! For symmetric stabilisation
    !   Coefficients(1:2, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! Note that this routine serves as a wrapper for block vectors and
    ! matrices. If there is only one block, then the corresponding
    ! scalar routine is called. Otherwise, an error is thrown.
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
    include 'intf_calcMatrixDiagSc_sim.inc'
    include 'intf_calcMatrixSc_sim.inc'
    optional :: fcb_calcMatrixDiagSc_sim

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorEdgeSc.inc'
    optional :: fcb_calcOperatorEdgeSc

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

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

    ! Check if user-defined callback function is present
    if (present(fcb_calcOperatorEdgeSc)) then
      
      call fcb_calcOperatorEdgeSc(rgroupFEMSet, rx, rmatrix,&
          dscale, bclear, fcb_calcMatrixDiagSc_sim, fcb_calcMatrixSc_sim,&
          rcollection, rafcstab)
      
    elseif ((rx%nblocks            .eq. 1) .and.&
            (rmatrix%nblocksPerCol .eq. 1) .and.&
            (rmatrix%nblocksPerRow .eq. 1)) then

      ! Call scalar version of this routine
      call gfsc_buildOperatorEdgeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixSc_sim, dscale, bclear, rmatrix%RmatrixBlock(1,1),&
          fcb_calcMatrixDiagSc_sim, cconstrType, rcollection, rafcstab,&
          rperfconfig=rperfconfig)
      
    else
      call output_line('Matrix/Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeBlock1')
      call sys_halt()
    end if
    
  end subroutine gfsc_buildOperatorEdgeBlock1

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorEdgeBlock2(rgroupFEMSet, rx,&
      fcb_calcMatrixSc_sim, dscale, bclear, rmatrix,&
      fcb_calcMatrixDiagSc_sim, cconstrType, rcollection,&
      rafcstab, fcb_calcOperatorEdgeSc, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only edge-by-edge assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! for each edge IJ the matrix entries and the artificial diffusion
    ! coefficient is stored according to the following convention:
    !
    ! For upwind-biased stabilisation
    !   Coefficients(1:3, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! For symmetric stabilisation
    !   Coefficients(1:2, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! Note that this routine serves as a wrapper for block matrices. If
    ! there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, an error is thrown.
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
    include 'intf_calcMatrixDiagSc_sim.inc'
    include 'intf_calcMatrixSc_sim.inc'
    optional :: fcb_calcMatrixDiagSc_sim

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorEdgeSc.inc'
    optional :: fcb_calcOperatorEdgeSc

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

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
    if (present(fcb_calcOperatorEdgeSc)) then

      ! Create auxiliary 1-block vector
      call lsysbl_createVecFromScalar(rx, rxBlock)
      
      call fcb_calcOperatorEdgeSc(rgroupFEMSet, rxBlock, rmatrix,&
          dscale, bclear, fcb_calcMatrixDiagSc_sim, fcb_calcMatrixSc_sim,&
          rcollection, rafcstab)
      
      ! Release auxiliary 1-block vector
      call lsysbl_releaseVector(rxBlock)
      
    elseif ((rmatrix%nblocksPerCol .eq. 1) .and.&
            (rmatrix%nblocksPerRow .eq. 1)) then
      
      ! Call scalar version of this routine
      call gfsc_buildOperatorEdgeScalar(rgroupFEMSet, rx,&
          fcb_calcMatrixSc_sim, dscale, bclear, rmatrix%RmatrixBlock(1,1),&
          fcb_calcMatrixDiagSc_sim, cconstrType, rcollection, rafcstab,&
          rperfconfig=rperfconfig)
      
    else
      call output_line('Matrix must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeBlock2')
      call sys_halt()
    end if
    
  end subroutine gfsc_buildOperatorEdgeBlock2

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorEdgeBlock3(rgroupFEMSet, rx,&
      fcb_calcMatrixSc_sim, dscale, bclear, rmatrix,&
      fcb_calcMatrixDiagSc_sim, cconstrType, rcollection,&
      rafcstab, fcb_calcOperatorEdgeSc, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only edge-by-edge assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! for each edge IJ the matrix entries and the artificial diffusion
    ! coefficient is stored according to the following convention:
    !
    ! For upwind-biased stabilisation
    !   Coefficients(1:3, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! For symmetric stabilisation
    !   Coefficients(1:2, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! Note that this routine serves as a wrapper for block vectors. If
    ! there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, an error is thrown.
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
    include 'intf_calcMatrixDiagSc_sim.inc'
    include 'intf_calcMatrixSc_sim.inc'
    optional :: fcb_calcMatrixDiagSc_sim

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorEdgeSc.inc'
    optional :: fcb_calcOperatorEdgeSc

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

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
    if (present(fcb_calcOperatorEdgeSc)) then
      
      ! Create auxiliary 1-block matrix
      call lsysbl_createMatFromScalar(rmatrix, rmatrixBlock)
      
      call fcb_calcOperatorEdgeSc(rgroupFEMSet, rx, rmatrixBlock,&
          dscale, bclear, fcb_calcMatrixDiagSc_sim, fcb_calcMatrixSc_sim,&
          rcollection, rafcstab)
      
      ! Release auxiliary 1-block matrix
      call lsysbl_releaseMatrix(rmatrixBlock)
      
    elseif (rx%nblocks .eq. 1) then
      
      ! Call scalar version of this routine
      call gfsc_buildOperatorEdgeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixSc_sim, dscale, bclear, rmatrix,&
          fcb_calcMatrixDiagSc_sim, cconstrType, rcollection,&
          rafcstab, rperfconfig=rperfconfig)
      
    else
      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeBlock2')
      call sys_halt()
    end if
    
  end subroutine gfsc_buildOperatorEdgeBlock3
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorEdgeScalar(rgroupFEMSet, rx,&
      fcb_calcMatrixSc_sim, dscale, bclear, rmatrix,&
      fcb_calcMatrixDiagSc_sim, cconstrType, rcollection,&
      rafcstab, fcb_calcOperatorEdgeSc, rperfconfig)

!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only edge-by-edge assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! for each edge IJ the matrix entries and the artificial diffusion
    ! coefficient is stored according to the following convention:
    !
    ! For upwind-biased stabilisation
    !   Coefficients(1:3, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! For symmetric stabilisation
    !   Coefficients(1:2, IJ) = (/d_ij, k_ij, k_ji/)
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
    include 'intf_calcMatrixDiagSc_sim.inc'
    include 'intf_calcMatrixSc_sim.inc'
    optional :: fcb_calcMatrixDiagSc_sim

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcOperatorEdgeSc.inc'
    optional :: fcb_calcOperatorEdgeSc

    ! OPTIONAL: One of the GFEM_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! GFEM_MATC_CONSISTENT is used.
    integer, intent(in), optional :: cconstrType

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
    real(DP), dimension(:,:), pointer :: p_Dcoefficients
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:,:), pointer :: p_IdiagList
    integer, dimension(:), pointer :: p_IedgeListIdx

    integer :: ccType

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    ! Check if user-defined assembly is provided
    if (present(fcb_calcOperatorEdgeSc)) then
      ! Create auxiliary 1-block vector and matrix
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createMatFromScalar(rmatrix, rmatrixBlock)
      
      ! Call user-defined assembly
      call fcb_calcOperatorEdgeSc(rgroupFEMSet, rxBlock, rmatrixBlock,&
          dscale, bclear, fcb_calcMatrixDiagSc_sim, fcb_calcMatrixSc_sim,&
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
      p_rperfconfig => gfsc_perfconfig
    end if

    ! Check if matrix and vector have the same data type
    if ((rmatrix%cdataType .ne. rx%cdataType) .or.&
        (rmatrix%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
      call sys_halt()
    end if
    
    ! Set type of matrix construction method
    ccType = GFEM_MATC_CONSISTENT
    if (present(cconstrType)) ccType = cconstrType

    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_EDGEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------
      
      if ((ccType .eq. GFEM_MATC_LUMPED) .and. bclear)&
          call lsyssc_clearMatrix(rmatrix)

      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA) .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
        call sys_halt()
      end if
        
      ! Set pointers
      call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
      call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)
      
      ! What data types are we?
      select case(rgroupFEMSet%cdataType)
      case (ST_DOUBLE)
        ! Set pointers
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        call lsyssc_getbase_double(rmatrix, p_Ddata)
        call lsyssc_getbase_double(rx, p_Dx)
        
        ! Assemble matrix diagonal?
        if (present(fcb_calcMatrixDiagSc_sim)) then
          
          ! Check if group finite element set is prepared
          if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGLIST) .eq. 0) .or.&
              (iand(rgroupFEMSet%isetSpec, GFEM_HAS_DIAGDATA) .eq. 0)) then
            call output_line('Group finite element set does not provide diagonal data!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
            call sys_halt()
          end if
          
          ! Assemble diagonal entries
          call gfem_getbase_IdiagList(rgroupFEMSet, p_IdiagList)
          call gfem_getbase_DcoeffsAtDiag(rgroupFEMSet, p_DcoeffsAtDiag)
          call doOperatorDiagDble(p_IdiagList, p_DcoeffsAtDiag,&
              p_Dx, dscale, bclear, p_Ddata)
        end if
        
        ! Do we have to build the stabilisation?
        if (present(rafcstab)) then
          
          ! Check if stabilisation has been prepared
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
            call output_line('Stabilisation has not been prepared!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
            call sys_halt()
          end if
          
          ! Check if coefficients should be stored in stabilisation
          if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE) then
            
            ! Check if stabilisation has the same data type
            if (rafcstab%cdataType .ne. ST_DOUBLE) then
              call output_line('Stabilisation must have double precision!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
              call sys_halt()
            end if
            
            ! Set additional pointers
            call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
            
            !-------------------------------------------------------------------
            ! Assemble operator with stabilisation and generate coefficients
            !-------------------------------------------------------------------
            call doOperatorAFCDble(p_IedgeListIdx, p_IedgeList,&
                p_DcoeffsAtEdge, p_Dx, dscale, bclear, ccType,&
                p_Ddata, p_Dcoefficients)
            
            ! Set state of stabilisation
            rafcstab%istabilisationSpec =&
                ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)
            
            ! Do we need edge orientation?
            if (rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED) then
              call afcstab_upwindOrientation(p_Dcoefficients, p_IedgeList,&
                  p_DcoeffsAtEdge, 2, 3)
              rafcstab%istabilisationSpec =&
                  ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION)
            else
              rafcstab%istabilisationSpec =&
                  iand(rafcstab%istabilisationSpec, not(AFCSTAB_HAS_EDGEORIENTATION))
            end if
            
          else
            
            !-------------------------------------------------------------------
            ! Assemble operator with stabilisation but do not generate coeffs
            !-------------------------------------------------------------------
            call doOperatorStabDble(p_IedgeListIdx, p_IedgeList,&
                  p_DcoeffsAtEdge, p_Dx, dscale, bclear, ccType, p_Ddata)
          end if
          
        else   ! no stabilisation structure present
          
          !---------------------------------------------------------------------
          ! Assemble operator without stabilisation
          !---------------------------------------------------------------------
          call doOperatorDble(p_IedgeListIdx, p_IedgeList,&
              p_DcoeffsAtEdge, p_Dx, dscale, bclear, ccType, p_Ddata)
        end if
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
      call sys_halt()
    end select
    
  contains

    ! Here, the working routine follow
    
    !**************************************************************
    ! Assemble diagonal part of the operator
    
    subroutine doOperatorDiagDble(IdiagList, DcoeffsAtDiag, Dx,&
        dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtDiag
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IdiagList

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:), pointer :: DdataAtNode
      real(DP), dimension(:,:), pointer :: Dcoefficients

      ! local variables
      integer :: IEQmax,IEQset,i,ia,idx

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------
      
      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,IEQmax,i,idx,ia)&
      !$omp if (size(IdiagList,2) > p_rperfconfig%NEQMIN_OMP)
      
      ! Allocate temporal memory
      allocate(DdataAtNode(p_rperfconfig%NEQSIM))
      allocate(Dcoefficients(1,p_rperfconfig%NEQSIM))
      
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
          DdataAtNode(idx) = Dx(i)
        end do

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSc_sim(&
            DdataAtNode(1:IEQmax-IEQset+1),&
            DcoeffsAtDiag(:,IEQset:IEQmax),&
            IdiagList(:,IEQset:IEQmax),&
            dscale, IEQmax-IEQset+1,&
            Dcoefficients(:,1:IEQmax-IEQset+1), rcollection)
        
        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then

          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ia = IdiagList(2,idx+IEQset-1)
            
            ! Update the diagonal coefficient
            Ddata(ia) = Dcoefficients(1,idx)
          end do
          
        else   ! do not clear matrix
          
          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ia = IdiagList(2,idx+IEQset-1)
            
            ! Update the diagonal coefficient
            Ddata(ia) = Ddata(ia) + Dcoefficients(1,idx)
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
    
    subroutine doOperatorDble(IedgeListIdx, IedgeList,&
        DcoeffsAtEdge, Dx, dscale, bclear, ccType, Ddata)
      
      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: ccType

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: Dcoefficients

      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ii,jj,ij,ji
      
      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtEdge,IEDGEmax,idx,iedge,ii,jj,ij,ji)&
      !$omp if(size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(2,p_rperfconfig%NEDGESIM))
      allocate(Dcoefficients(2,p_rperfconfig%NEDGESIM))

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
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSc_sim(&
              DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dcoefficients(:,1:IEDGEmax-IEDGEset+1), rcollection)
          
          ! What type of assembly are we?
          if (ccType .eq. GFEM_MATC_LUMPED) then

            ! Loop through all edges in the current set and scatter
            ! the entries to the diagonal of the global matrix
            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)

              ! Update the global operator
              Ddata(ii) = Ddata(ii) + Dcoefficients(1,idx)
              Ddata(jj) = Ddata(jj) + Dcoefficients(2,idx)
            end do
            
          else
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
                Ddata(ij) = Dcoefficients(1,idx)
                Ddata(ji) = Dcoefficients(2,idx)
              end do
              
            else   ! do not clear matrix
              
              do idx = 1, IEDGEmax-IEDGEset+1
                
                ! Get actual edge number
                iedge = idx+IEDGEset-1
                
                ! Get position of off-diagonal entries
                ij = IedgeList(3,iedge)
                ji = IedgeList(4,iedge)
                
                ! Update the global operator
                Ddata(ij) = Ddata(ij) + Dcoefficients(1,idx)
                Ddata(ji) = Ddata(ji) + Dcoefficients(2,idx)
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
    ! Assemble operator edge-by-edge with stabilisation.
    
    subroutine doOperatorStabDble(IedgeListIdx, IedgeList,&
        DcoeffsAtEdge, Dx, dscale, bclear, ccType, Ddata)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: ccType

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: Dcoefficients
      
      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ii,ij,ji,jj

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtEdge,IEDGEmax,idx,iedge,ii,ij,ji,jj)&
      !$omp if(size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)
      
      ! Allocate temporal memory
      allocate(DdataAtEdge(2,p_rperfconfig%NEDGESIM))
      allocate(Dcoefficients(3,p_rperfconfig%NEDGESIM))
      
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
          
          ! Dcoefficients is not allocated as temporal array
          ! since it is given as global output parameter.
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSc_sim(&
              DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dcoefficients(:,1:IEDGEmax-IEDGEset+1), rcollection)
          
          ! What type of assembly are we?
          if (ccType .eq. GFEM_MATC_LUMPED) then

            ! Loop through all edges in the current set and scatter
            ! the entries to the diagonal of the global matrix
            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) + Dcoefficients(2,idx)
              Ddata(jj) = Ddata(jj) + Dcoefficients(3,idx)
            end do

          else

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
                Ddata(ii) = Ddata(ii) - Dcoefficients(1,idx)
                Ddata(jj) = Ddata(jj) - Dcoefficients(1,idx)
                Ddata(ij) =             Dcoefficients(2,idx) + Dcoefficients(1,idx)
                Ddata(ji) =             Dcoefficients(3,idx) + Dcoefficients(1,idx)
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
                Ddata(ii) = Ddata(ii) - Dcoefficients(1,idx)
                Ddata(jj) = Ddata(jj) - Dcoefficients(1,idx)
                Ddata(ij) = Ddata(ij) + Dcoefficients(2,idx) + Dcoefficients(1,idx)
                Ddata(ji) = Ddata(ji) + Dcoefficients(3,idx) + Dcoefficients(1,idx)
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

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation and AFC data.
    
    subroutine doOperatorAFCDble(IedgeListIdx, IedgeList,&
        DcoeffsAtEdge, Dx, dscale, bclear, ccType, Ddata, Dcoefficients)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: ccType

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! output parameters
      real(DP), dimension(:,:), intent(out) :: Dcoefficients
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtEdge
      
      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ii,ij,ji,jj

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,IEDGEmax,idx,iedge,ii,ij,ji,jj)&
      !$omp if(size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)
      
      ! Allocate temporal memory
      allocate(DdataAtEdge(2,p_rperfconfig%NEDGESIM))
      
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
          
          ! Dcoefficients is not allocated as temporal array
          ! since it is given as global output parameter.
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSc_sim(&
              DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dcoefficients(:,IEDGEset:IEDGEmax), rcollection)
          
          ! What type of assembly are we?
          if (ccType .eq. GFEM_MATC_LUMPED) then
            
            ! Loop through all edges in the current set and scatter
            ! the entries to the diagonal of the global matrix
            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = IedgeList(5,iedge)
              jj = IedgeList(6,iedge)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) + Dcoefficients(2,iedge)
              Ddata(jj) = Ddata(jj) + Dcoefficients(3,iedge)
            end do
            
          else

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
                
                ! Compute entries of low-order operator
                Dcoefficients(2,iedge) = Dcoefficients(2,iedge)&
                                       + Dcoefficients(1,iedge)
                Dcoefficients(3,iedge) = Dcoefficients(3,iedge)&
                                       + Dcoefficients(1,iedge)
                
                ! Update the global operator
                Ddata(ii) = Ddata(ii) - Dcoefficients(1,iedge)
                Ddata(jj) = Ddata(jj) - Dcoefficients(1,iedge)
                Ddata(ij) =             Dcoefficients(2,iedge)
                Ddata(ji) =             Dcoefficients(3,iedge)
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
                
                ! Compute entries of low-order operator
                Dcoefficients(2,iedge) = Dcoefficients(2,iedge)&
                                       + Dcoefficients(1,iedge)
                Dcoefficients(3,iedge) = Dcoefficients(3,iedge)&
                                       + Dcoefficients(1,iedge)
                
                ! Update the global operator
                Ddata(ii) = Ddata(ii) - Dcoefficients(1,iedge)
                Ddata(jj) = Ddata(jj) - Dcoefficients(1,iedge)
                Ddata(ij) = Ddata(ij) + Dcoefficients(2,iedge)
                Ddata(ji) = Ddata(ji) + Dcoefficients(3,iedge)
              end do
              
            end if
          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      !$omp end parallel
      
    end subroutine doOperatorAFCDble
        
  end subroutine gfsc_buildOperatorEdgeScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildVectorNodeBlock(rgroupFEMSet, rx,&
      fcb_calcVectorSc_sim, dscale, bclear, rvector,&
      rcollection, rafcstab, fcb_calcVectorNodeSc, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete vector by the group
    ! finite element formulation. The vector entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling a user-defined callback function.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! Note that this routine serves as a wrapper for block vectors. If
    ! there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, an error is thrown.
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
    include 'intf_calcVectorSc_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcVectorNodeSc.inc'
    optional :: fcb_calcVectorNodeSc

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

    ! Check if user-defined callback function is present
    if (present(fcb_calcVectorNodeSc)) then

      call fcb_calcVectorNodeSc(rgroupFEMSet, rx, rvector, dscale,&
          bclear, fcb_calcVectorSc_sim, rcollection, rafcstab)
      
    elseif ((rx%nblocks .eq. 1) .and. (rvector%nblocks .eq. 1)) then

      ! Call scalar version of this routine
      call gfsc_buildVectorNodeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcVectorSc_sim, dscale, bclear, rvector%RvectorBlock(1),&
          rcollection, rafcstab, rperfconfig=rperfconfig)

    else
      call output_line('Vectors must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorNodeBlock')
      call sys_halt()
    end if

  end subroutine gfsc_buildVectorNodeBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildVectorNodeScalar(rgroupFEMSet, rx,&
      fcb_calcVectorSc_sim, dscale, bclear, rvector,&
      rcollection, rafcstab, fcb_calcVectorNodeSc, rperfconfig)
    
!<description>
    ! This subroutine assembles a discrete vector by the group
    ! finite element formulation. The vector entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling a user-defined callback function.
    !
    ! This routine supports only node-by-node assembly.
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
    include 'intf_calcVectorSc_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcVectorNodeSc.inc'
    optional :: fcb_calcVectorNodeSc

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
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    real(DP), dimension(:), pointer :: p_Dx,p_Ddata
    integer, dimension(:,:), pointer :: p_InodeListIdx2D,p_InodeList2D
    integer, dimension(:), pointer :: p_InodeListIdx1D,p_InodeList1D
    
    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    ! Check if user-defined assembly is provided
    if (present(fcb_calcVectorNodeSc)) then
      ! Create auxiliary 1-block vectors
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createVecFromScalar(rvector, rvectorBlock)
      
      ! Call user-defined assembly
      call fcb_calcVectorNodeSc(rgroupFEMSet, rxBlock, rvectorBlock,&
          dscale, bclear, fcb_calcVectorSc_sim, rcollection, rafcstab)
      
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
      p_rperfconfig => gfsc_perfconfig
    end if

    ! Check if vectors have the same data type double
    if ((rx%cdataType .ne. rvector%cdataType) .or.&
        (rx%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorNodeScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorNodeScalar')
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
          call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx1D)
          call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList1D)

          ! Assemble vector node-by-node
          call doVectorDble(p_InodeListIdx1D, p_InodeList1D,&
              p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)
        else
          ! Set pointers
          call gfem_getbase_InodeListIdx(rgroupFEMSet, p_InodeListIdx2D)
          call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList2D)

          ! Assemble selected part of the vector node-by-node
          call doVectorDbleSel(p_InodeListIdx2D, p_InodeList2D,&
              p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)
        end if
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorNodeScalar')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorNodeScalar')
      call sys_halt()
    end select

  contains
    
    !**************************************************************
    ! Assemble vector node-by-node without stabilisation

    subroutine doVectorDble(InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: InodeListIdx,InodeList

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:), pointer :: DdataAtNode,Dcoefficients
      integer, dimension(:,:), pointer  :: IdofsAtNode
      
      ! local variables
      real(DP) :: dtemp
      integer :: IAmax,IApos,IAset,IEQmax,IEQset,ia,idx,ieq

      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------
      
      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,dtemp,IdofsAtNode,&
      !$omp         IAmax,IApos,IAset,IEQmax,ia,idx,ieq)&
      !$omp if(size(InodeList) > p_rperfconfig%NAMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtNode(p_rperfconfig%NEQSIM))
      allocate(Dcoefficients(p_rperfconfig%NEQSIM))
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
              IdofsAtNode(1,idx) = InodeList(ia)     ! absolut nodal value j
              IdofsAtNode(2,idx) = ia                ! absolute matrix position ia
              DdataAtNode(idx)   = Dx(InodeList(ia)) ! solution value at node j
            end do
            
            ! Increase the upper bound for nonzero entries
            IAmax = IAmax+1
          end do
          
          ! Use callback function to compute matrix entries
          call fcb_calcVectorSc_sim(&
              DdataAtNode(1:idx),&
              DcoeffsAtNode(:,IApos:IApos+idx-1),&
              IdofsAtNode(:,1:idx),&
              dscale, idx,&
              Dcoefficients(1:idx), rcollection)
          
          ! Initialise local index which will run from IAset..IAmax
          idx = 0

          ! Loop through all equations in the current set
          ! and scatter the entries to the global matrix
          do ieq = IAset, IAmax-1
            
            ! Clear temporal data
            dtemp = 0.0_DP

            ! Loop over all contributions to this equation
            do ia = InodeListIdx(ieq), InodeListIdx(ieq+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Update temporal data
              dtemp = dtemp + Dcoefficients(idx)
            end do

            ! Update the global vector
            Ddata(ieq) = Ddata(ieq)+dtemp
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

    subroutine doVectorDbleSel(InodeListIdx, InodeList,&
        DcoeffsAtNode, Dx, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: InodeListIdx,InodeList

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:), pointer :: DdataAtNode,Dcoefficients
      
      ! local variables
      real(DP) :: dtemp
      integer :: IAmax,IApos,IAset,IEQmax,IEQset,ia,idx,ieq,iidx

      !-------------------------------------------------------------------------
      ! Assemble all entries
      !-------------------------------------------------------------------------
      
      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtNode,dtemp,&
      !$omp         IAmax,IApos,IAset,IEQmax,ia,idx,ieq,iidx)&
      !$omp if(size(InodeList,2) > p_rperfconfig%NAMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtNode(p_rperfconfig%NEQSIM))
      allocate(Dcoefficients(p_rperfconfig%NEQSIM))

      ! Loop over all equations in blocks of size NEQSIM
      !$omp do schedule(static,1)
      do IEQset = 1, size(InodeListIdx,2)-1, p_rperfconfig%NEQSIM
        
        ! We always handle NEQSIM equations by one OpenMP thread.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle
        ! at most NEQSIM equations by one OpenMP thread.

        IEQmax = min(size(InodeListIdx,2)-1, IEQset-1+p_rperfconfig%NEQSIM)
        
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
        IApos = InodeListIdx(1,IEQset)
        
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
            if (InodeListIdx(1,IAmax+1)-IApos .gt. p_rperfconfig%NASIM) exit
            
            ! Loop through all nonzero matrix entries in the current
            ! equation and prepare the auxiliary arrays for it
            do ia = InodeListIdx(1,IAmax), InodeListIdx(1,IAmax+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Fill auxiliary arrays
              DdataAtNode(idx) = Dx(InodeList(1,ia)) ! solution value at node j
            end do
            
            ! Increase the upper bound for nonzero entries
            IAmax = IAmax+1
          end do
          
          ! Use callback function to compute matrix entries
          call fcb_calcVectorSc_sim(&
              DdataAtNode(1:idx),&
              DcoeffsAtNode(:,IApos:IApos+idx-1),&
              InodeList(:,IApos:IApos+idx-1),&
              dscale, idx,&
              Dcoefficients(1:idx), rcollection)
          
          ! Initialise local index which will run from IAset..IAmax
          idx = 0

          ! Loop through all equations in the current set
          ! and scatter the entries to the global matrix
          do iidx = IAset, IAmax-1
            
            ! Get actual node number
            ieq = InodeListIdx(2,iidx)
            
            ! Clear temporal date
            dtemp = 0.0_DP

            ! Loop over all contributions to this equation
            do ia = InodeListIdx(1,iidx), InodeListIdx(1,iidx+1)-1
              
              ! Update local index
              idx = idx+1
              
              ! Update temporal data
              dtemp = dtemp + Dcoefficients(idx)
            end do
            
            ! Update the global vector
            Ddata(ieq) = Ddata(ieq)+dtemp
          end do
          
          ! Proceed with next nonzero entries in current set
          IAset = IAmax
          IApos = InodeListIdx(1,IASet)
          
        end do
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(DdataAtNode)
      deallocate(Dcoefficients)
      !$omp end parallel

    end subroutine doVectorDbleSel

  end subroutine gfsc_buildVectorNodeScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildVectorEdgeBlock(rgroupFEMSet, rx,&
      fcb_calcFluxSc_sim, dscale, bclear, rvector,&
      rcollection, rafcstab, fcb_calcVectorEdgeSc, rperfconfig)
    
!<description>
    ! This subroutine assembles a vector operator by the group
    ! finite element formulation. The vector entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling a user-defined callback function.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! for each edge IJ the vector entries and the artificial diffusion
    ! coefficient is stored according to the following convention:
    !
    ! For upwind-biased stabilisation
    !   Coefficients(1:3, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! For symmetric stabilisation
    !   Coefficients(1:2, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! Note that this routine serves as a wrapper for block vectors. If
    ! there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, an error is thrown.
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

    ! Callback function to compute vector entries
    include 'intf_calcFluxSc_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcVectorEdgeSc.inc'
    optional :: fcb_calcVectorEdgeSc

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

    ! Check if user-defined callback function is present
    if (present(fcb_calcVectorEdgeSc)) then

      call fcb_calcVectorEdgeSc(rgroupFEMSet, rx, rvector, dscale,&
          bclear, fcb_calcFluxSc_sim, rcollection, rafcstab)

    elseif ((rx%nblocks .eq. 1) .and. (rvector%nblocks .eq. 1)) then

      ! Call scalar version of this routine
      call gfsc_buildVectorEdgeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcFluxSc_sim, dscale, bclear, rvector%RvectorBlock(1),&
          rcollection, rafcstab, rperfconfig=rperfconfig)

    else
      call output_line('Vectors must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorEdgeBlock')
      call sys_halt()
    end if

  end subroutine gfsc_buildVectorEdgeBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildVectorEdgeScalar(rgroupFEMSet, rx,&
      fcb_calcFluxSc_sim, dscale, bclear, rvector,&
      rcollection, rafcstab, fcb_calcVectorEdgeSc, rperfconfig)
    
!<description>
    ! This subroutine assembles a vector operator by the group
    ! finite element formulation. The vector entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling a user-defined callback function.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! for each edge IJ the vector entries and the artificial diffusion
    ! coefficient is stored according to the following convention:
    !
    ! For upwind-biased stabilisation
    !   Coefficients(1:3, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! For symmetric stabilisation
    !   Coefficients(1:2, IJ) = (/d_ij, k_ij, k_ji/)
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

    ! Callback function to compute vector entries
    include 'intf_calcFluxSc_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcVectorEdgeSc.inc'
    optional :: fcb_calcVectorEdgeSc

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
    real(DP), dimension(:,:), pointer :: p_Dcoefficients
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    ! Check if user-defined assembly is provided
    if (present(fcb_calcVectorEdgeSc)) then
      ! Create auxiliary 1-block vector
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createVecFromScalar(rvector, rvectorBlock)
      
      ! Call user-defined assembly
      call fcb_calcVectorEdgeSc(rgroupFEMSet, rxBlock, rvectorBlock,&
          dscale, bclear, fcb_calcFluxSc_sim, rcollection, rafcstab)
      
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
      p_rperfconfig => gfsc_perfconfig
    end if

    ! Check if vectors have the same data type double
    if ((rx%cdataType .ne. rvector%cdataType) .or.&
        (rx%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorEdgeScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorEdgeScalar')
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
        
        ! Do we have to build the stabilisation?
        if (present(rafcstab)) then
          
          ! Check if stabilisation has been prepared
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
            call output_line('Stabilisation has not been prepared!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorEdgeScalar')
            call sys_halt()
          end if
          
          ! Check if coefficients should be stored in stabilisation
          if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE) then
            
            ! Check if stabilisation has the same data type
            if (rafcstab%cdataType .ne. ST_DOUBLE) then
              call output_line('Stabilisation must have double precision!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorEdgeScalar')
              call sys_halt()
            end if
            
            ! Set additional pointers
            call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
            
            !-------------------------------------------------------------------
            ! Assemble vector with stabilisation and generate coefficients
            !-------------------------------------------------------------------
            call doVectorDble(p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge,&
                p_Dx, dscale, bclear, p_Ddata, p_Dcoefficients)
            
            ! Set state of stabilisation
            rafcstab%istabilisationSpec =&
                ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)
            
            ! Do we need edge orientation?
            if (rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED) then
              call afcstab_upwindOrientation(p_Dcoefficients, p_IedgeList,&
                  p_DcoeffsAtEdge, 2, 3)
              rafcstab%istabilisationSpec =&
                  ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION)
            else
              rafcstab%istabilisationSpec =&
                  iand(rafcstab%istabilisationSpec, not(AFCSTAB_HAS_EDGEORIENTATION))
            end if
            
          else
            
            !-------------------------------------------------------------------
            ! Assemble vector without stabilisation
            !-------------------------------------------------------------------
            call doVectorDble(p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge,&
                p_Dx, dscale, bclear, p_Ddata)
          end if
          
        else   ! no stabilisation structure present
          
          !---------------------------------------------------------------------
          ! Assemble vector without stabilisation
          !---------------------------------------------------------------------
          call doVectorDble(p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge,&
              p_Dx, dscale, bclear, p_Ddata)
          
        end if
        
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorEdgeScalar')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorEdgeScalar')
      call sys_halt()
    end select
      
  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble vector edge-by-edge without stabilisation

    subroutine doVectorDble(IedgeListIdx, IedgeList, DcoeffsAtEdge,&
        Dx, dscale, bclear, Ddata, Dcoefficients)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata

      ! output parameters
      real(DP), dimension(:,:), intent(out), optional :: Dcoefficients

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: DfluxesAtEdge
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      if (bclear) call lalg_clearVector(Ddata)

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DfluxesAtEdge,IEDGEmax,i,idx,iedge,j)&
      !$omp if(size(IedgeList,2) > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(2,p_rperfconfig%NEDGESIM))

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
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute internodal fluxes
          if (present(Dcoefficients)) then
            call fcb_calcFluxSc_sim(&
                DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
                DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
                IedgeList(:,IEDGEset:IEDGEmax),&
                dscale, IEDGEmax-IEDGEset+1,&
                DfluxesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
                Dcoefficients(:,IEDGEset:IEDGEmax), rcollection)
          else
            call fcb_calcFluxSc_sim(&
                DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
                DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
                IedgeList(:,IEDGEset:IEDGEmax),&
                dscale, IEDGEmax-IEDGEset+1,&
                DfluxesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
                rcollection=rcollection)
          end if
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Update the global vector
            Ddata(i) = Ddata(i)+DfluxesAtEdge(1,idx)
            Ddata(j) = Ddata(j)+DfluxesAtEdge(2,idx)
          end do
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      !$omp end parallel

    end subroutine doVectorDble

  end subroutine gfsc_buildVectorEdgeScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianBlock(RcoeffMatrices, rx, fcb_calcMatrixSc_sim,&
      hstep, dscale, bbuildStabilisation, bclear, rjacobian, rcollection)

!<description>
    ! This subroutine assembles the Jacobian matrix for the convective
    ! part of the discrete transport operator for a scalar convection
    ! equation. Note that this routine serves as a wrapper for block
    ! vectors. If there is only one block, then the corresponding
    ! scalar routine is called. Otherwise, an error is thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx
    
    ! perturbation parameter
    real(DP), intent(in) :: hstep

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(in) :: bbuildStabilisation

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute matrix entries
    include 'intf_calcMatrixSc_sim.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    if (rx%nblocks .ne. 1) then
      
      call output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianBlock')
      call sys_halt()

    else
      
      call gfsc_buildJacobianScalar(&
          RcoeffMatrices, rx%RvectorBlock(1), fcb_calcMatrixSc_sim, hstep,&
          dscale, bbuildStabilisation, bclear, rjacobian, rcollection)
      
    end if

  end subroutine gfsc_buildJacobianBlock

   !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianScalar(RcoeffMatrices, rx, fcb_calcMatrixSc_sim,&
      hstep, dscale, bbuildStabilisation, bclear, rjacobian, rcollection)

!<description>
    ! This subroutine assembles the Jacobian matrix for the convective part
    ! of the discrete transport operator for a scalar convection equation.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx
    
    ! perturbation parameter
    real(DP), intent(in) :: hstep

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(in) :: bbuildStabilisation

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute matrix entries
    include 'intf_calcMatrixSc_sim.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Jac,p_Dx
    integer :: h_Ksep,ndim
    
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)
    
    ! Set pointers
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      
    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_DcoeffZ)

    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalar')
      call sys_halt()
    end select
    
    
    ! What kind of matrix are we?
    select case(rjacobian%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      
      ! Do we have to build the upwind Jacobian?
      if (bbuildStabilisation) then
        
        select case(ndim)
        case (NDIM1D)
          call doUpwindMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doUpwindMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doUpwindMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Jac)
        end select

      else   ! bbuildStabilisation

        select case(ndim)
        case (NDIM1D)
          call doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Jac)
        end select

      end if   ! bbuildStabilisation

      ! Release diagonal separator
      call storage_free(h_Ksep)
      
      
    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      
      ! Do we have to build the upwind Jacobian?
      if (bbuildStabilisation) then
        
        select case(ndim)
        case (NDIM1D)
          call doUpwindMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doUpwindMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doUpwindMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Jac)
        end select
      
      else   ! bbuildStabilisation

        select case(ndim)
        case (NDIM1D)
          call doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Jac)
        end select

      end if   ! bbuildStabilisation

      ! Release diagonal separator
      call storage_free(h_Ksep)

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow
    
    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doGalerkinMat7_1D(Kld, Kcol, Ksep, NEQ, DcoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients k_ij and k_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij ,k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat7_1D


    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doGalerkinMat7_2D(Kld, Kcol, Ksep, NEQ, DcoeffX, DcoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients k_ij and k_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij ,k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat7_2D

    
    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doGalerkinMat7_3D(Kld, Kcol, Ksep, NEQ, DcoeffX, DcoeffY, DcoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
          C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients k_ij and k_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat7_3D

    
    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doGalerkinMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      
      
      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat9_1D


    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doGalerkinMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, DcoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat9_2D


    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doGalerkinMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, DcoeffY, DcoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
          C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat9_3D


    !**************************************************************
    ! Assemble upwind Jacobian matrix for convective
    ! operator in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_1D(Kld, Kcol, Ksep, NEQ, DcoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j
      
      
      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat7_1D
    

    !**************************************************************
    ! Assemble upwind Jacobian matrix for convective
    ! operator in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_2D(Kld, Kcol, Ksep, NEQ, DcoeffX, DcoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat7_2D


    !**************************************************************
    ! Assemble Jacobian matrix for convective
    ! operator in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_3D(Kld, Kcol, Ksep, NEQ, DcoeffX, DcoeffY, DcoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j

      
      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
          C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat7_3D
    
    
    !**************************************************************
    ! Assemble Jacobian matrix for convective
    ! operator in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j

      
      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat9_1D


    !**************************************************************
    ! Assemble Jacobian matrix for convective
    ! operator in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, DcoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat9_2D


    !**************************************************************
    ! Assemble Jacobian matrix for convective
    ! operator in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, DcoeffY, DcoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
          C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij =(a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat9_3D

  end subroutine gfsc_buildJacobianScalar

end module groupfemscalar
