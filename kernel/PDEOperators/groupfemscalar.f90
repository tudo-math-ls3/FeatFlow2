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
!# 1.) gfsc_buildOperator = gfsc_buildOperatorConst /
!#                          gfsc_buildOperatorNodeScalar /
!#                          gfsc_buildOperatorNodeBlock /
!#                          gfsc_buildOperatorEdgeScalar /
!#                          gfsc_buildOperatorEdgeBlock
!#     -> Assembles a discrete operator by the group finite element formulation
!#
!# 2.) gfsc_buildVector = gfsc_buildVectorNodeScalar /
!#                        gfsc_buildVectorNodeBlock /
!#                        gfsc_buildVectorEdgeScalar /
!#                        gfsc_buildVectorEdgeBlock
!#     -> Assembles a discrete vector by the group finite element formulation
!#
!# 3.) gfsc_buildJacobian = gfsc_buildJacobianScalar /
!#                          gfsc_buildJacobianBlock
!#      -> Assembles the Jacobian matrix for the group finite element formulation
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
  use spatialdiscretisation
  use storage

  implicit none

  private
  
  public :: gfsc_buildOperator
  public :: gfsc_buildVector
  public :: gfsc_buildJacobian
  
!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of nodes to handle simultaneously when building matrices
#ifndef GFSC_NEQSIM
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSC_NEQSIM = 128
#else
  integer, public            :: GFSC_NEQSIM = 128
#endif
#endif

  ! Number of edges to handle simultaneously when building matrices
#ifndef GFSC_NEDGESIM
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSC_NEDGESIM = 64
#else
  integer, public            :: GFSC_NEDGESIM = 64
#endif
#endif
  
  ! Minimum number of nodes for OpenMP parallelisation: If the number of
  ! nodes is below this value, then no parallelisation is performed.
#ifndef GFSC_NEQMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSC_NEQMIN_OMP = 1000
#else
  integer, public            :: GFSC_NEQMIN_OMP = 1000
#endif
#endif

  ! Minimum number of edges for OpenMP parallelisation: If the number of
  ! edges is below this value, then no parallelisation is performed.
#ifndef GFSC_NEDGEMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSC_NEDGEMIN_OMP = 1000
#else
  integer, public            :: GFSC_NEDGEMIN_OMP = 1000
#endif
#endif
!</constantblock>

!</constants>

  ! ****************************************************************************

  interface gfsc_buildOperator
    module procedure gfsc_buildOperatorConst
    module procedure gfsc_buildOperatorNodeScalar
    module procedure gfsc_buildOperatorNodeBlock
    module procedure gfsc_buildOperatorEdgeScalar
    module procedure gfsc_buildOperatorEdgeBlock
  end interface
  
  interface gfsc_buildVector
!    module procedure gfsc_buildVectorNodeScalar
!    module procedure gfsc_buildVectorNodeBlock
    module procedure gfsc_buildVectorEdgeScalar
    module procedure gfsc_buildVectorEdgeBlock
  end interface

  interface gfsc_buildJacobian
    module procedure gfsc_buildJacobianScalar
    module procedure gfsc_buildJacobianBlock
  end interface

contains
 
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorConst(rgroupFEMSet, dscale, bclear, rmatrix, rafcstab)

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
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge

    real(SP), dimension(:), pointer :: p_Fdata
    real(SP), dimension(:,:), pointer :: p_Fcoefficients
    real(SP), dimension(:,:), pointer :: p_FcoeffsAtNode
    real(SP), dimension(:,:,:), pointer :: p_FcoeffsAtEdge

    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:,:), pointer :: p_InodeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer, dimension(:), pointer :: p_Kdiagonal,p_Kld

    logical :: bsymm

    ! Check if matrix and vector have the same data type
    if (rmatrix%cdataType .ne. rgroupFEMSet%cdataType) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
      call sys_halt()
    end if

    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_NODEBASED)
      !-------------------------------------------------------------------------
      ! Node-based assembly
      !-------------------------------------------------------------------------

      ! Check if group finite element set is prepared
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA) .eq. 0) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
        call sys_halt()
      end if

      if (present(rafcstab)) then
        call output_line('Stabilisation is not feasible in node-by-node assembly!',&
            OU_CLASS_WARNING,OU_MODE_STD,'gfsc_buildOperatorConst')
      end if
      
      ! What kind of matrix are we?
      select case(rmatrix%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrix, p_Kld)

        ! What data types are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Set pointers
          call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
          call lsyssc_getbase_double(rmatrix, p_Ddata)

          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOpNodeMat79Dble(p_Kld, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOpNodeMat79Dble(p_Kld, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata, p_InodeList)
          end if

        case (ST_SINGLE)
          ! Set pointers
          call gfem_getbase_FcoeffsAtNode(rgroupFEMSet, p_FcoeffsAtNode)
          call lsyssc_getbase_single(rmatrix, p_Fdata)

          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOpNodeMat79Sngl(p_Kld, rgroupFEMSet%NEQ,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOpNodeMat79Sngl(p_Kld, rgroupFEMSet%NEQ,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata, p_InodeList)
          end if

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
        call sys_halt()
      end select

      
    case (GFEM_EDGEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------

      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGESTRUCTURE) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA)      .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA)      .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
        call sys_halt()
      end if
      
      ! Set pointers
      call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
      call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)
      
      ! What kind of matrix are we?
      select case(rmatrix%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-----------------------------------------------------------------------
        
        ! Set diagonal pointer
        if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX7) then
          call lsyssc_getbase_Kld(rmatrix, p_Kdiagonal)
        else
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        end if
        
        ! What data types are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Set pointers
          call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
          call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
          call lsyssc_getbase_double(rmatrix, p_Ddata)

          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOpDiagMat79Dble(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOpDiagMat79Dble(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata, p_InodeList)
          end if
          

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
            if (rafcstab%h_CoefficientsAtEdge .ne. ST_NOHANDLE) then
              
              ! Check if stabilisation has the same data type
              if (rafcstab%cdataType .ne. ST_DOUBLE) then
                call output_line('Stabilisation must have double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
                call sys_halt()
              end if
              
              ! Set additional pointers
              call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
              
              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation and generate coefficients
              !-----------------------------------------------------------------
              call doOpEdgeAFCMat79Dble(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, dscale, bclear, bsymm,&
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
              
              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation but do not generate coeffs
              !-----------------------------------------------------------------
              call doOpEdgeStabMat79Dble(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, dscale, bclear, bsymm, p_Ddata)
            end if

          else   ! no stabilisation structure present

            !-------------------------------------------------------------------
            ! Assemble operator without stabilisation
            !-------------------------------------------------------------------
            call doOpEdgeMat79Dble(p_IedgeListIdx, p_IedgeList,&
                rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, dscale, bclear, p_Ddata)
          end if

        case (ST_SINGLE)
          ! Set pointers
          call gfem_getbase_FcoeffsAtNode(rgroupFEMSet, p_FcoeffsAtNode)
          call gfem_getbase_FcoeffsAtEdge(rgroupFEMSet, p_FcoeffsAtEdge)
          call lsyssc_getbase_single(rmatrix, p_Fdata)
          
          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOpDiagMat79Sngl(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOpDiagMat79Sngl(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata, p_InodeList)
          end if

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
            if (rafcstab%h_CoefficientsAtEdge .ne. ST_NOHANDLE) then
              
              ! Check if stabilisation has the same data type
              if (rafcstab%cdataType .ne. ST_SINGLE) then
                call output_line('Stabilisation must have double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
                call sys_halt()
              end if
              
              ! Set additional pointers
              call afcstab_getbase_FcoeffsAtEdge(rafcstab, p_Fcoefficients)
              
              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation
              !-----------------------------------------------------------------
              call doOpEdgeAFCMat79Sngl(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_FcoeffsAtEdge, real(dscale,SP),&
                  bclear, bsymm, p_Fdata, p_Fcoefficients)
              
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
              
              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation but do not generate coeffs
              !-----------------------------------------------------------------
              call doOpEdgeStabMat79Sngl(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_FcoeffsAtEdge, real(dscale,SP),&
                  bclear, bsymm, p_Fdata)
            end if
            
          else   ! no stabilisation structure present

            !-------------------------------------------------------------------
            ! Assemble operator without stabilisation
            !-------------------------------------------------------------------
            call doOpEdgeMat79Sngl(p_IedgeListIdx, p_IedgeList, rgroupFEMSet%NEDGE,&
                p_FcoeffsAtEdge, real(dscale,SP), bclear, p_Fdata)
          end if
          
        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
          call sys_halt()
        end select
        
      case default
        call output_line('Unsupported matrix format!',&
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
    ! Assemble operator node-by-node without stabilisation
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOpNodeMat79Dble(Kld, NEQ, DcoeffsAtNode,&
        dscale, bclear, Ddata, InodeList)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! local variables
      integer :: idx,ieq,ij

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------

      if (bclear) call lalg_clearVector(Ddata)

      if (present(InodeList)) then
        
        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(ij)
        do idx = 1, size(InodeList,2)

          ! Get position of matrix entry
          ij  = InodeList(2,idx)

          ! Update the matrix coefficient
          Ddata(ij) = Ddata(ij) + dscale*DcoeffsAtNode(1,idx)
        end do
        !$omp end parallel do

      else
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(ij)
        do ieq = 1, NEQ
          do ij = Kld(ieq), Kld(ieq+1)-1

            ! Update the matrix coefficient
            Ddata(ij) = Ddata(ij) + dscale*DcoeffsAtNode(1,ij)
          end do
        end do
        !$omp end parallel do

      end if

    end subroutine doOpNodeMat79Dble

    !**************************************************************
    ! Assemble operator node-by-node without stabilisation
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOpNodeMat79Sngl(Kld, NEQ, FcoeffsAtNode,&
        fscale, bclear, Fdata, InodeList)

      ! input parameters
      real(SP), dimension(:,:), intent(in) :: FcoeffsAtNode
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! local variables
      integer :: idx,ieq,ij

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------

      if (bclear) call lalg_clearVector(Fdata)

      if (present(InodeList)) then
        
        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(ij)
        do idx = 1, size(InodeList,2)

          ! Get position of matrix entry
          ij  = InodeList(2,idx)

          ! Update the matrix coefficient
          Fdata(ij) = Fdata(ij) + fscale*FcoeffsAtNode(1,idx)
        end do
        !$omp end parallel do

      else
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(ij)
        do ieq = 1, NEQ
          do ij = Kld(ieq), Kld(ieq+1)-1

            ! Update the matrix coefficient
            Fdata(ij) = Fdata(ij) + dscale*FcoeffsAtNode(1,ij)
          end do
        end do
        !$omp end parallel do

      end if

    end subroutine doOpNodeMat79Sngl

    !**************************************************************
    ! Assemble diagonal part of the operator
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOpDiagMat79Dble(Kdiagonal, NEQ, DcoeffsAtNode,&
        dscale, bclear, Ddata, InodeList)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! local variables
      integer :: ieq,ii,idx

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then

        if (present(InodeList)) then

          ! Loop over the subset of equations
          !$omp parallel do default(shared) private(ii,ieq)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do idx = 1, NEQ

            ! Get equation number and  position of diagonal entry
            ieq = InodeList(1,idx)
            ii  = InodeList(2,idx)
            
            ! Update the diagonal coefficient
            Ddata(ii) = dscale*DcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        else

          ! Loop over all equations
          !$omp parallel do default(shared) private(ii)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do ieq = 1, NEQ
            
            ! Get position of diagonal entry
            ii = Kdiagonal(ieq)
            
            ! Update the diagonal coefficient
            Ddata(ii) = dscale*DcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do
        end if

      else   ! do not clear matrix

        if (present(InodeList)) then

          ! Loop over the subset of equations
          !$omp parallel do default(shared) private(ii,ieq)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do idx = 1, NEQ

            ! Get equation number and position of diagonal entry
            ieq = InodeList(1,idx)
            ii  = InodeList(2,idx)
            
            ! Update the diagonal coefficient
            Ddata(ii) = Ddata(ii) + dscale*DcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        else

          ! Loop over all equations
          !$omp parallel do default(shared) private(ii)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do ieq = 1, NEQ
            
            ! Get position of diagonal entry
            ii = Kdiagonal(ieq)
            
            ! Update the diagonal coefficient
            Ddata(ii) = Ddata(ii) + dscale*DcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        end if

      end if

    end subroutine doOpDiagMat79Dble

    !**************************************************************
    ! Assemble diagonal part of the operator
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOpDiagMat79Sngl(Kdiagonal, NEQ, FcoeffsAtNode,&
        fscale, bclear, Fdata, InodeList)

      ! input parameters
      real(SP), dimension(:,:), intent(in) :: FcoeffsAtNode
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! local variables
      integer :: ieq,ii,idx

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then

        if (present(InodeList)) then

          ! Loop over the subset of equations
          !$omp parallel do default(shared) private(ii,ieq)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do idx = 1, NEQ

            ! Get equation number and  position of diagonal entry
            ieq = InodeList(1,idx)
            ii  = InodeList(2,idx)
            
            ! Update the diagonal coefficient
            Fdata(ii) = fscale*FcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        else

          ! Loop over all equations
          !$omp parallel do default(shared) private(ii)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do ieq = 1, NEQ
            
            ! Get position of diagonal entry
            ii = Kdiagonal(ieq)
            
            ! Update the diagonal coefficient
            Fdata(ii) = fscale*FcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do
        end if

      else   ! do not clear matrix

        if (present(InodeList)) then

          ! Loop over the subset of equations
          !$omp parallel do default(shared) private(ii,ieq)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do idx = 1, NEQ

            ! Get equation number and position of diagonal entry
            ieq = InodeList(1,idx)
            ii  = InodeList(2,idx)
            
            ! Update the diagonal coefficient
            Fdata(ii) = Fdata(ii) + fscale*FcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        else

          ! Loop over all equations
          !$omp parallel do default(shared) private(ii)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do ieq = 1, NEQ
            
            ! Get position of diagonal entry
            ii = Kdiagonal(ieq)
            
            ! Update the diagonal coefficient
            Fdata(ii) = Fdata(ii) + fscale*FcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        end if

      end if

    end subroutine doOpDiagMat79Sngl
    
    !**************************************************************
    ! Assemble operator edge-by-edge without stabilisation
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeMat79Dble(IedgeListIdx, IedgeList, NEDGE,&
        DcoeffsAtEdge, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata

      ! local variables
      integer :: iedge,igroup,ij,ji     

      !$omp parallel default(shared) private(ij,ji)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
 
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

    end subroutine doOpEdgeMat79Dble

    !**************************************************************
    ! Assemble operator edge-by-edge without stabilisation
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeMat79Sngl(IedgeListIdx, IedgeList, NEDGE,&
        FcoeffsAtEdge, fscale, bclear, Fdata)
      
      ! input parameters
      real(SP), dimension(:,:,:), intent(in) :: FcoeffsAtEdge
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata

      ! local variables
      integer :: iedge,igroup,ij,ji
      
      !$omp parallel default(shared) private(ij,ji)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

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

    end subroutine doOpEdgeMat79Sngl

    !**************************************************************
    ! Assemble edge-by-edge operator with stabilisation
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeStabMat79Dble(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, DcoeffsAtEdge, dscale, bclear, bsymm, Ddata)
      
      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear, bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata

      ! local variables
      real(DP) :: d_ij
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Symmetric artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Non-symmetric artificial diffusion coefficient
              d_ij = max(0.0_DP,&
                  -DcoeffsAtEdge(1,1,iedge), -DcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Symmetric artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = Ddata(ji) + dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Non-symmetric artificial diffusion coefficient
              d_ij = max(0.0_DP,&
                  -DcoeffsAtEdge(1,1,iedge), -DcoeffsAtEdge(1,2,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = Ddata(ji) + dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        end if

      end if
      !$omp end parallel

    end subroutine doOpEdgeStabMat79Dble

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeStabMat79Sngl(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, FcoeffsAtEdge, fscale, bclear, bsymm, Fdata)

      ! input parameters
      real(SP), dimension(:,:,:), intent(in) :: FcoeffsAtEdge
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear, bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata

      ! local variables
      real(SP) :: d_ij
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Symmetric artificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Non-symmetric artificial diffusion coefficient
              d_ij = max(0.0_SP,&
                  -FcoeffsAtEdge(1,1,iedge), -FcoeffsAtEdge(1,2,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Symmetric atificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = Fdata(ji) + fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Non-symmetric atificial diffusion coefficient
              d_ij = max(0.0_SP,&
                  -FcoeffsAtEdge(1,1,iedge), -FcoeffsAtEdge(1,2,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = Fdata(ji) + fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup
          
        end if

      end if
      !$omp end parallel

    end subroutine doOpEdgeStabMat79Sngl

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation and AFC data.
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeAFCMat79Dble(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, DcoeffsAtEdge, dscale, bclear, bsymm, Ddata, Dcoefficients)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear,bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! output parameters
      real(DP), dimension(:,:), intent(out) :: Dcoefficients
      
      ! local variables
      real(DP) :: d_ij,k_ij,k_ji
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij,k_ij,k_ji)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_DP,  DcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Dcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_DP,  DcoeffsAtEdge(1,1,iedge))
              k_ji = max(0.0_DP,  DcoeffsAtEdge(1,2,iedge))
              
              ! Non-symmetric AFC w/o edge orientation
              Dcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_DP,  DcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Dcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = Ddata(ji) + dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_DP,  DcoeffsAtEdge(1,1,iedge))
              k_ji = max(0.0_DP,  DcoeffsAtEdge(1,2,iedge))
              
              ! Non-symmetric AFC w/o edge orientation
              Dcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = Ddata(ji) + dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        end if
        
      end if
      !$omp end parallel

    end subroutine doOpEdgeAFCMat79Dble

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation and AFC data.
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeAFCMat79Sngl(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, FcoeffsAtEdge, fscale, bclear, bsymm, Fdata, Fcoefficients)

      ! input parameters
      real(SP), dimension(:,:,:), intent(in) :: FcoeffsAtEdge
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear,bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! output parameters
      real(SP), dimension(:,:), intent(out) :: Fcoefficients
      
      ! local variables
      real(SP) :: d_ij,k_ij,k_ji
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij,k_ij,k_ji)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_SP,  FcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Fcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_SP,  FcoeffsAtEdge(1,1,iedge))
              k_ji = max(0.0_SP,  FcoeffsAtEdge(1,2,iedge))
              
              ! Non-symmetric AFC w/o edge orientation
              Fcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_SP,  FcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Fcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = Fdata(ji) + fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
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
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_SP,  FcoeffsAtEdge(1,1,iedge))
              k_ji = max(0.0_SP,  FcoeffsAtEdge(1,2,iedge))
              
              ! Non-symmetric AFC w/o edge orientation
              Fcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = Fdata(ji) + fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        end if
        
      end if
      !$omp end parallel

    end subroutine doOpEdgeAFCMat79Sngl

  end subroutine gfsc_buildOperatorConst

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorNodeBlock(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix, rcollection)
    
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
    ! is called.  Otherwise, an error is thrown.
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
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>
    
    ! Check if block vector contains exactly one block
    if (rx%nblocks .ne. 1) then

      call output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeBlock')
      call sys_halt()

    else

      call gfsc_buildOperatorNodeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix, rcollection)

    end if
    
  end subroutine gfsc_buildOperatorNodeBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorNodeScalar(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix, rcollection)
    
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
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata,p_Dx
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    
    integer, dimension(:,:), pointer :: p_InodeList
    integer, dimension(:), pointer :: p_Kld

    ! Check if matrix and vector have the same data type
    if ((rmatrix%cdataType .ne. rx%cdataType) .or.&
        (rmatrix%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
      call sys_halt()
    end if

    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_NODEBASED)
      !-------------------------------------------------------------------------
      ! Node-based assembly
      !-------------------------------------------------------------------------

      ! Check if group finite element set is prepared
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA) .eq. 0) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
        call sys_halt()
      end if

      ! What kind of matrix are we?
      select case(rmatrix%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrix, p_Kld)

        ! What data types are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Set pointers
          call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
          call lsyssc_getbase_double(rmatrix, p_Ddata)
          call lsyssc_getbase_double(rx, p_Dx)

          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOperatorMat79Dble(p_Kld, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOperatorMat79Dble(p_Kld, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata, p_InodeList)
          end if

           case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported matrix format!',&
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
    ! Assemble operator node-by-node without stabilisation
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Dble(Kld, NEQ, DcoeffsAtNode, Dx,&
        dscale, bclear, Ddata, InodeList)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
    end subroutine doOperatorMat79Dble

  end subroutine gfsc_buildOperatorNodeScalar
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorEdgeBlock(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, fcb_calcMatrixSc_sim,&
      dscale, bclear, rmatrix, rcollection, rafcstab)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
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
    ! is called.  Otherwise, an error is thrown.
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

    ! Check if block vector contains exactly one block
    if (rx%nblocks .ne. 1) then

      call output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeBlock')
      call sys_halt()

    else

      call gfsc_buildOperatorEdgeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSc_sim, fcb_calcMatrixSc_sim, dscale,&
          bclear, rmatrix, rcollection, rafcstab)

    end if

  end subroutine gfsc_buildOperatorEdgeBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorEdgeScalar(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, fcb_calcMatrixSc_sim, dscale,&
      bclear, rmatrix, rcollection, rafcstab)

!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
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
    real(DP), dimension(:), pointer :: p_Ddata,p_Dx
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    real(DP), dimension(:,:), pointer :: p_Dcoefficients
    
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:,:), pointer :: p_InodeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer, dimension(:), pointer :: p_Kdiagonal

    ! Check if matrix and vector have the same data type
    if ((rmatrix%cdataType .ne. rx%cdataType) .or.&
        (rmatrix%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
      call sys_halt()
    end if
    
    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_EDGEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------
      
      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGESTRUCTURE) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA)      .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA)      .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
        call sys_halt()
      end if
        
      ! Set pointers
      call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
      call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)
      
      ! What kind of matrix are we?
      select case(rmatrix%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-----------------------------------------------------------------------
        
        ! Set diagonal pointer
        if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX7) then
          call lsyssc_getbase_Kld(rmatrix, p_Kdiagonal)
        else
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        end if
        
        ! What data types are we?
        select case(rmatrix%cdataType)
          
        case (ST_DOUBLE)
          ! Set pointers
          call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
          call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
          call lsyssc_getbase_double(rmatrix, p_Ddata)
          call lsyssc_getbase_double(rx, p_Dx)
        
          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOpDiagMat79Dble(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOpDiagMat79Dble(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata, p_InodeList)
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
            if (rafcstab%h_CoefficientsAtEdge .ne. ST_NOHANDLE) then

              ! Check if stabilisation has the same data type
              if (rafcstab%cdataType .ne. ST_DOUBLE) then
                call output_line('Stabilisation must have double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
                call sys_halt()
              end if
                   
              ! Set additional pointers
              call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)

              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation and generate coefficients
              !-----------------------------------------------------------------
              call doOperatorAFCMat79Dble(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, p_Dx, dscale, bclear,&
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

              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation but do not generate coeffs
              !-----------------------------------------------------------------
              call doOperatorStabMat79Dble(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, p_Dx, dscale, bclear, p_Ddata)
            end if
            
          else   ! no stabilisation structure present

            !-------------------------------------------------------------------
            ! Assemble operator without stabilisation
            !-------------------------------------------------------------------
            call doOperatorMat79Dble(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, p_Dx, dscale, bclear, p_Ddata)
          end if
          
        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
          call sys_halt()
        end select
              
      case default
        call output_line('Unsupported matrix format!',&
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
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOpDiagMat79Dble(Kdiagonal, NEQ, DcoeffsAtNode, Dx,&
        dscale, bclear, Ddata, InodeList)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:), pointer :: DdataAtNode
      real(DP), dimension(:,:), pointer :: DcoefficientsAtNode
      integer, dimension(:,:), pointer  :: IdofsAtNode

      ! local variables
      integer :: idx,IEQset,IEQmax
      integer :: i,iidx,ii

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------
      
      !$omp parallel default(shared)&
      !$omp private(DcoefficientsAtNode,DdataAtNode,IEQmax,&
      !$omp         IdofsAtNode,i,idx,ii)

      ! Allocate temporal memory
      allocate(IdofsAtNode(2,GFSC_NEQSIM))
      allocate(DdataAtNode(GFSC_NEQSIM))
      allocate(DcoefficientsAtNode(1,GFSC_NEQSIM))

      ! Loop over the equations
      !$omp do schedule(static,1)
      do IEQset = 1, NEQ, GFSC_NEQSIM

        ! We always handle GFSC_NEQSIM equations simultaneously.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most GFSC_NEQSIM equations simultaneously.
        
        IEQmax = min(NEQ, IEQset-1+GFSC_NEQSIM)
        
        if (present(InodeList)) then

          ! Loop through all equations in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEQmax-IEQset+1
            
            ! Get actual index
            iidx = idx+IEQset-1
            
            ! Get equation number and position of diagonal entry
            i  = InodeList(1,iidx)
            ii = InodeList(2,iidx)
            
            ! Fill auxiliary arrays
            IdofsAtNode(1,idx) = i
            IdofsAtNode(2,idx) = ii
            DdataAtNode(idx)       = Dx(i)
          end do

        else

          ! Loop through all equations in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEQmax-IEQset+1

            ! Get actual equation number
            i = idx+IEQset-1
            
            ! Get position of diagonal entry
            ii = Kdiagonal(i)
            
            ! Fill auxiliary arrays
            IdofsAtNode(1,idx) = i
            IdofsAtNode(2,idx) = ii
            DdataAtNode(idx)       = Dx(i)
          end do

        end if

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSc_sim(&
            DdataAtNode(1:IEQmax-IEQset+1),&
            DcoeffsAtNode(:,IEQset:IEQmax),&
            IdofsAtNode(:,1:IEQmax-IEQset+1),&
            dscale, IEQmax-IEQset+1,&
            DcoefficientsAtNode(:,1:IEQmax-IEQset+1), rcollection)

        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ii = IdofsAtNode(2,idx)
            
            ! Update the diagonal coefficient
            Ddata(ii) = DcoefficientsAtNode(1,idx)
          end do

        else   ! do not clear matrix
          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ii = IdofsAtNode(2,idx)
            
            ! Update the diagonal coefficient
            Ddata(ii) = Ddata(ii) + DcoefficientsAtNode(1,idx)
          end do

        end if
      end do
      !$omp end do

      ! Deallocate temporal memory
      deallocate(IdofsAtNode)
      deallocate(DdataAtNode)
      deallocate(DcoefficientsAtNode)
      !$omp end parallel

    end subroutine doOpDiagMat79Dble

    !**************************************************************
    ! Assemble operator edge-by-edge without stabilisation
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOperatorMat79Dble(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, DcoeffsAtEdge, Dx, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: Dcoefficients

      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ij,ji
      
      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtEdge,IEDGEmax,idx,iedge,ij,ji)

      ! Allocate temporal memory
      allocate(DdataAtEdge(2,GFSC_NEDGESIM))
      allocate(Dcoefficients(2,GFSC_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSC_NEDGESIM

          ! We always handle GFSC_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSC_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+GFSC_NEDGESIM)
        
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
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorMat79Dble
    
    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation.
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOperatorStabMat79Dble(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, DcoeffsAtEdge, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

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
      !$omp private(Dcoefficients,DdataAtEdge,IEDGEmax,idx,iedge,ii,ij,ji,jj)
      
      ! Allocate temporal memory
      allocate(DdataAtEdge(2,GFSC_NEDGESIM))
      allocate(Dcoefficients(3,GFSC_NEDGESIM))
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSC_NEDGESIM

          ! We always handle GFSC_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSC_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+GFSC_NEDGESIM)
          
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
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global matrix
          if (bclear) then
            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - Dcoefficients(1,idx)
              Ddata(jj) = Ddata(jj) - Dcoefficients(1,idx)
              Ddata(ij) = Dcoefficients(2,idx) + Dcoefficients(1,idx) 
              Ddata(ji) = Dcoefficients(3,idx) + Dcoefficients(1,idx)
            end do

          else   ! do not clear matrix

            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - Dcoefficients(1,idx)
              Ddata(jj) = Ddata(jj) - Dcoefficients(1,idx)
              Ddata(ij) = Ddata(ij) + Dcoefficients(2,idx) + Dcoefficients(1,idx) 
              Ddata(ji) = Ddata(ji) + Dcoefficients(3,idx) + Dcoefficients(1,idx)
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorStabMat79Dble

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation and AFC data.
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOperatorAFCMat79Dble(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, DcoeffsAtEdge, Dx, dscale, bclear, Ddata, Dcoefficients)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

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
      !$omp private(DdataAtEdge,IEDGEmax,idx,iedge,ii,ij,ji,jj)
      
      ! Allocate temporal memory
      allocate(DdataAtEdge(2,GFSC_NEDGESIM))
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSC_NEDGESIM

          ! We always handle GFSC_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSC_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+GFSC_NEDGESIM)
          
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
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global matrix
          if (bclear) then
            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Compute entries of low-order operator
              Dcoefficients(2,iedge) = Dcoefficients(2,iedge) + Dcoefficients(1,iedge)
              Dcoefficients(3,iedge) = Dcoefficients(3,iedge) + Dcoefficients(1,iedge)

              ! Update the global operator
              Ddata(ii) = Ddata(ii) - Dcoefficients(1,iedge)
              Ddata(jj) = Ddata(jj) - Dcoefficients(1,iedge)
              Ddata(ij) = Dcoefficients(2,iedge)
              Ddata(ji) = Dcoefficients(3,iedge)
            end do

          else   ! do not clear matrix

            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Compute entries of low-order operator
              Dcoefficients(2,iedge) = Dcoefficients(2,iedge) + Dcoefficients(1,iedge)
              Dcoefficients(3,iedge) = Dcoefficients(3,iedge) + Dcoefficients(1,iedge)

              ! Update the global operator
              Ddata(ii) = Ddata(ii) - Dcoefficients(1,iedge)
              Ddata(jj) = Ddata(jj) - Dcoefficients(1,iedge)
              Ddata(ij) = Ddata(ij) + Dcoefficients(2,iedge)
              Ddata(ji) = Ddata(ji) + Dcoefficients(3,iedge)
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      !$omp end parallel
      
    end subroutine doOperatorAFCMat79Dble
        
  end subroutine gfsc_buildOperatorEdgeScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildVectorNodeBlock(rgroupFEMSet, rx,&
      fcb_calcVectorSc_sim, dscale, bclear, rvector, rcollection)
    
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
    ! is called.  Otherwise, an error is thrown.
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
!</input>

!<inputoutput>
    ! Destination vector
    type(t_VectorBlock), intent(inout) :: rvector

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! Check if both block vectors contain exactly one block
    if ((rx%nblocks .ne. 1) .or. (rvector%nblocks .ne. 1)) then

      call output_line('Vectors must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorNodeBlock')
      call sys_halt()

    else

      call gfsc_buildVectorNodeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcVectorSc_sim, dscale, bclear, rvector%RvectorBlock(1),&
          rcollection)

    end if

  end subroutine gfsc_buildVectorNodeBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildVectorNodeScalar(rgroupFEMSet, rx,&
      fcb_calcVectorSc_sim, dscale, bclear, rvector, rcollection)
    
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
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorScalar), intent(inout) :: rvector

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    print *, "Not implemented yet!"
    stop

  end subroutine gfsc_buildVectorNodeScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildVectorEdgeBlock(rgroupFEMSet, rx,&
      fcb_calcFluxSc_sim, dscale, bclear, rvector, rcollection, rafcstab)
    
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
    ! is called.  Otherwise, an error is thrown.
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

    ! Check if both block vectors contain exactly one block
    if ((rx%nblocks .ne. 1) .or. (rvector%nblocks .ne. 1)) then

      call output_line('Vectors must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVectorEdgeBlock')
      call sys_halt()

    else

      call gfsc_buildVectorEdgeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcFluxSc_sim, dscale, bclear, rvector%RvectorBlock(1),&
          rcollection, rafcstab)

    end if

  end subroutine gfsc_buildVectorEdgeBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildVectorEdgeScalar(rgroupFEMSet, rx,&
      fcb_calcFluxSc_sim, dscale, bclear, rvector, rcollection, rafcstab)
    
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
    real(DP), dimension(:), pointer :: p_Dx,p_Ddata
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    real(DP), dimension(:,:), pointer :: p_Dcoefficients
    
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

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
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGESTRUCTURE) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA)      .eq. 0)) then
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
          if (rafcstab%h_CoefficientsAtEdge .ne. ST_NOHANDLE) then
            
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
                p_Dx, dscale, p_Ddata, p_Dcoefficients)
            
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
            ! Assemble operator without stabilisation
            !-------------------------------------------------------------------
            call doVectorDble(p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge,&
                p_Dx, dscale, p_Ddata)
          end if
          
        else   ! no stabilisation structure present
          
          !---------------------------------------------------------------------
          ! Assemble operator without stabilisation
          !---------------------------------------------------------------------
          call doVectorDble(p_IedgeListIdx, p_IedgeList, p_DcoeffsAtEdge,&
              p_Dx, dscale, p_Ddata)
          
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
        Dx, dscale, Ddata, Dcoefficients)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
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

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DfluxesAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(2,GFSC_NEDGESIM))
      allocate(DfluxesAtEdge(2,GFSC_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSC_NEDGESIM

          ! We always handle GFSC_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSC_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSC_NEDGESIM)
          
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
                rcollection=rcollection)
          else
            call fcb_calcFluxSc_sim(&
                DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
                DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
                IedgeList(:,IEDGEset:IEDGEmax),&
                dscale, IEDGEmax-IEDGEset+1,&
                DfluxesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
                Dcoefficients(:,IEDGEset:IEDGEmax), rcollection)
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
    ! equation.  Note that this routine serves as a wrapper for block
    ! vectors. If there is only one block, then the corresponding
    ! scalar routine is called.  Otherwise, an error is thrown.
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

    case DEFAULT
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

    case DEFAULT
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
