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
!# 1.) gfsys_buildOperator = gfsys_buildOperatorScalar /
!#                           gfsys_buildOperatorBlock
!#     -> Assembles a discrete operator by the group finite element formulation
!#
!# 2.) gfsys_buildVector = gfsys_buildVectorScalar /
!#                         gfsys_buildVectorBlock
!#     -> Assembles a discrete vector by the group finite element formulation
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
  use spatialdiscretisation
  use storage
  
  implicit none

  private

  public :: gfsys_buildOperator
  public :: gfsys_buildVector

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of nodes to handle simultaneously when building matrices
#ifndef GFSYS_NEQSIM
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSYS_NEQSIM = 128
#else
  integer, public            :: GFSYS_NEQSIM = 128
#endif
#endif

  ! Number of edges to handle simultaneously when building matrices
#ifndef GFSYS_NEDGESIM
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSYS_NEDGESIM = 64
#else
  integer, public            :: GFSYS_NEDGESIM = 64
#endif
#endif
  
  ! Minimum number of nodes for OpenMP parallelisation: If the number of
  ! nodes is below this value, then no parallelisation is performed.
#ifndef GFSYS_NEQMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSYS_NEQMIN_OMP = 1000
#else
  integer, public            :: GFSYS_NEQMIN_OMP = 1000
#endif
#endif

  ! Minimum number of edges for OpenMP parallelisation: If the number of
  ! edges is below this value, then no parallelisation is performed.
#ifndef GFSYS_NEDGEMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSYS_NEDGEMIN_OMP = 1000
#else
  integer, public            :: GFSYS_NEDGEMIN_OMP = 1000
#endif
#endif
!</constantblock>

!</constants>

  ! ****************************************************************************

  interface gfsys_buildOperator
     module procedure gfsys_buildOperatorScalar
     module procedure gfsys_buildOperatorBlock
  end interface

  interface gfsys_buildVector
    module procedure gfsys_buildVectorScalar
    module procedure gfsys_buildVectorBlock
  end interface

contains

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildOperatorBlock(rafcstab, rx,&
      fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, dscale,&
      bclear, rdivOp, rcollection, fcb_calcDivOperator)

!<description>
    ! This subroutine assembles the discrete divergence operator which results
    ! from the group finite element formulation of the continuous problem
    !
    !   <tex> $$ \nabla\cdot{\bf F}(u) $$ </tex>
    !
    ! where ${\bf f}(U)$ is a user-defined flux function for the
    ! multi-component field $U$.
    !
    ! This routine can be used to apply the following discretisations:
    !
    ! (1) the standard Galerkin finite element method
    !     which will be referred to as high-order approximation
    !
    ! (2) discrete upwinding for hyperbolic systems which results from
    !     the conservative elimination of negative eigenvalues from
    !     the Galerkin operator. This technique is for instance
    !     described in the reference:
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
    ! The solution vector
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

    ! OPTIONAL: callback  function to overwrite the standard operation
    include 'intf_calcDivOperator.inc'
    optional :: fcb_calcDivOperator
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! The divergence operator
    type(t_matrixBlock), intent(inout) :: rdivOp

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_array), dimension(:,:), allocatable  :: rarray
    real(DP), dimension(:), pointer :: p_Dx
    real(DP), dimension(:,:), pointer :: p_DmatrixCoeffsAtNode
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx,p_Kdiagonal
    logical :: bisFullMatrix


    ! Check if block vector contains only one block, if global
    ! operator is stored in interleave format and no user-defined
    ! callback function is provided.
    if (.not.present(fcb_calcDivOperator) .and.&
        (rx%nblocks .eq. 1)               .and.&
        (rdivOp%nblocksPerCol .eq. 1) .and.&
        (rdivOp%nblocksPerRow .eq. 1)) then
      call gfsys_buildOperatorScalar(rafcstab, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
          dscale, bclear, rdivOp%RmatrixBlock(1,1), rcollection,&
          fcb_calcDivOperator)
      return
    end if

    ! Check if block matrix exhibits group structure
    if (rdivOp%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorBlock')
      call sys_halt()
    end if

    ! Check if user-defined assembly is provided
    if (present(fcb_calcDivOperator)) then
      ! Call used-defined assembly
      call fcb_calcDivOperator(rafcstab, rx, rdivOp, dscale, bclear,&
          fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, rcollection)
    else
      ! Allocate temporal memory
      allocate(rarray(rx%nblocks,rx%nblocks))
      
      ! Set pointers
      call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
      call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
      call afcstab_getbase_DmatCoeffAtNode(rafcstab, p_DmatrixCoeffsAtNode)
      call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
      call afcstab_getbase_array(rdivOp, rarray, bisFullMatrix)
      call lsysbl_getbase_double(rx, p_Dx)
      
      ! What kind of matrix are we?
      select case(rdivOp%RmatrixBlock(1,1)%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-----------------------------------------------------------------------
        
        ! Set diagonal pointer
        if (rdivOp%RmatrixBlock(1,1)%cmatrixFormat .eq. LSYSSC_MATRIX7) then
          call lsyssc_getbase_Kld(rdivOp%RmatrixBlock(1,1), p_Kdiagonal)
        else
          call lsyssc_getbase_Kdiagonal(rdivOp%RmatrixBlock(1,1), p_Kdiagonal)
        end if
        
        ! What type of matrix are we?
        if (bisFullMatrix) then
          
          call doOperatorMat79(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEQ, rx%nblocks, p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge,&
              p_Dx, dscale, bclear, rarray)
          
        else   ! bisFullMatrix == no
          
          call doOperatorMat79Diag(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEQ, rx%nblocks, p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge,&
              p_Dx, dscale, bclear, rarray)
          
        end if   ! bisFullMatrix
        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorBlock')
        call sys_halt()
      end select
      
      ! Deallocate temporal memory
      deallocate(rarray)
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble block-diagonal divergence operator K
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Diag(Kdiagonal, IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtNode, DmatrixCoeffsAtEdge, Dx, dscale, bclear, K)
      
      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: K

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtNode
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtEdge
      integer, dimension(:,:), pointer  :: IdofsAtNode
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,IEQmax,IEQset,idx,igroup
      integer :: i,iedge,ii,ij,ivar,ji,jj

      
      !$omp parallel default(shared)&
      !$omp private(DcoefficientsAtEdge,DcoefficientsAtNode,DdataAtEdge,&
      !$omp         DdataAtNode,IEDGEmax,IEQmax,IdofsAtNode,i,idx,&
      !$omp         iedge,ii,ij,ivar,ji,jj)

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(IdofsAtNode(2,GFSYS_NEQSIM))
      allocate(DdataAtNode(NVAR,GFSYS_NEQSIM))
      allocate(DcoefficientsAtNode(NVAR,1,GFSYS_NEQSIM))

      ! Loop over the equations
      !$omp do schedule(static,1)
      do IEQset = 1, NEQ, GFSYS_NEQSIM

        ! We always handle GFSYS_NEQSIM equations simultaneously.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most GFSYS_NEQSIM equations simultaneously.
        
        IEQmax = min(NEQ, IEQset-1+GFSYS_NEQSIM)
        
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
          DdataAtNode(:,idx)     = Dx(i,:)
        end do

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DmatrixCoeffsAtNode(:,IEQset:IEQmax),&
            IdofsAtNode(:,1:IEQmax-IEQset+1),&
            dscale, IEQmax-IEQset+1,&
            DcoefficientsAtNode(:,:,1:IEQmax-IEQset+1), rcollection)
        
        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          do idx = 1, IEQmax-IEQset+1

            ! Get position of diagonal entry
            ii = IdofsAtNode(2,idx)
            
            ! Update the diagonal coefficient
            do ivar = 1, NVAR
              K(ivar,ivar)%p_Ddata(ii) = DcoefficientsAtNode(ivar,1,idx)
            end do
          end do

        else   ! do not clear matrix
          do idx = 1, IEQmax-IEQset+1

            ! Get position of diagonal entry
            ii = IdofsAtNode(2,idx)
            
            ! Update the diagonal coefficient
            do ivar = 1, NVAR
              K(ivar,ivar)%p_Ddata(ii) = K(ivar,ivar)%p_Ddata(ii)+&
                  DcoefficientsAtNode(ivar,1,idx)
            end do
          end do

        end if
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(IdofsAtNode)
      deallocate(DdataAtNode)
      deallocate(DcoefficientsAtNode)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcoefficientsAtEdge(NVAR,3,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DcoefficientsAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)
          
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
              do ivar = 1, NVAR
                K(ivar,ivar)%p_Ddata(ii) = K(ivar,ivar)%p_Ddata(ii)-&
                    DcoefficientsAtEdge(ivar,1,idx)
                K(ivar,ivar)%p_Ddata(jj) = K(ivar,ivar)%p_Ddata(jj)-&
                    DcoefficientsAtEdge(ivar,1,idx)
                K(ivar,ivar)%p_Ddata(ij) = &
                    DcoefficientsAtEdge(ivar,2,idx) + DcoefficientsAtEdge(ivar,1,idx) 
                K(ivar,ivar)%p_Ddata(ji) = &
                    DcoefficientsAtEdge(ivar,3,idx) + DcoefficientsAtEdge(ivar,1,idx) 
              end do
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
              do ivar = 1, NVAR
                K(ivar,ivar)%p_Ddata(ii) = K(ivar,ivar)%p_Ddata(ii)-&
                    DcoefficientsAtEdge(ivar,1,idx)
                K(ivar,ivar)%p_Ddata(jj) = K(ivar,ivar)%p_Ddata(jj)-&
                    DcoefficientsAtEdge(ivar,1,idx)
                K(ivar,ivar)%p_Ddata(ij) = K(ivar,ivar)%p_Ddata(ij)+&
                    DcoefficientsAtEdge(ivar,2,idx) + DcoefficientsAtEdge(ivar,1,idx) 
                K(ivar,ivar)%p_Ddata(ji) = K(ivar,ivar)%p_Ddata(ji)+&
                    DcoefficientsAtEdge(ivar,3,idx) + DcoefficientsAtEdge(ivar,1,idx) 
              end do
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcoefficientsAtEdge)
      !$omp end parallel

    end subroutine doOperatorMat79Diag

    
    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79(Kdiagonal, IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtNode, DmatrixCoeffsAtEdge, Dx, dscale, bclear, K)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: K

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtNode
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtEdge
      integer, dimension(:,:), pointer  :: IdofsAtNode
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,IEQmax,IEQset,idx,igroup
      integer :: i,iedge,ii,ij,ijpos,ivar,ji,jj,jvar

      
      !$omp parallel default(shared)&
      !$omp private(DcoefficientsAtEdge,DcoefficientsAtNode,DdataAtEdge,&
      !$omp         DdataAtNode,IEDGEmax,IEQmax,IdofsAtNode,i,idx,&
      !$omp         iedge,ii,ij,ijpos,ivar,jvar,ji,jj)

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(IdofsAtNode(2,GFSYS_NEQSIM))
      allocate(DdataAtNode(NVAR,GFSYS_NEQSIM))
      allocate(DcoefficientsAtNode(NVAR*NVAR,1,GFSYS_NEQSIM))

      ! Loop over the equations
      !$omp do schedule(static,1)
      do IEQset = 1, NEQ, GFSYS_NEQSIM

        ! We always handle GFSYS_NEQSIM equations simultaneously.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most GFSYS_NEQSIM equations simultaneously.
        
        IEQmax = min(NEQ, IEQset-1+GFSYS_NEQSIM)
        
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
          DdataAtNode(:,idx)     = Dx(i,:)
        end do

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DmatrixCoeffsAtNode(:,IEQset:IEQmax),&
            IdofsAtNode(:,1:IEQmax-IEQset+1),&
            dscale, IEQmax-IEQset+1,&
            DcoefficientsAtNode(:,:,1:IEQmax-IEQset+1), rcollection)

        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ii = IdofsAtNode(2,idx)
            
            ! Update the diagonal coefficient
            do ivar = 1, NVAR
              do jvar = 1, NVAR
                ijpos = NVAR*(ivar-1)+jvar
                K(jvar,ivar)%p_Ddata(ii) = DcoefficientsAtNode(ijpos,1,idx)
              end do
            end do
          end do

        else   ! do not clear matrix
          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ii = IdofsAtNode(2,idx)
            
            ! Update the diagonal coefficient
            do ivar = 1, NVAR
              do jvar = 1, NVAR
                ijpos = NVAR*(ivar-1)+jvar
                K(jvar,ivar)%p_Ddata(ii) = K(jvar,ivar)%p_Ddata(ii)+&
                    DcoefficientsAtNode(ijpos,1,idx)
              end do
            end do
          end do

        end if
      end do
      !$omp end do

      ! Deallocate temporal memory
      deallocate(IdofsAtNode)
      deallocate(DdataAtNode)
      deallocate(DcoefficientsAtNode)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcoefficientsAtEdge(NVAR*NVAR,3,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DcoefficientsAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)
          
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
              do ivar = 1, NVAR
                do jvar = 1, NVAR
                  ijpos = NVAR*(ivar-1)+jvar
                  K(jvar,ivar)%p_Ddata(ii) = K(jvar,ivar)%p_Ddata(ii)-&
                    DcoefficientsAtEdge(ijpos,1,idx)
                  K(jvar,ivar)%p_Ddata(jj) = K(jvar,ivar)%p_Ddata(jj)-&
                      DcoefficientsAtEdge(ijpos,1,idx)
                  K(jvar,ivar)%p_Ddata(ij) = &
                      DcoefficientsAtEdge(ijpos,2,idx) + DcoefficientsAtEdge(ijpos,1,idx) 
                  K(jvar,ivar)%p_Ddata(ji) = &
                      DcoefficientsAtEdge(ijpos,3,idx) + DcoefficientsAtEdge(ijpos,1,idx) 
                end do
              end do
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
              do ivar = 1, NVAR
                do jvar = 1, NVAR
                  ijpos = NVAR*(ivar-1)+jvar
                  K(jvar,ivar)%p_Ddata(ii) = K(jvar,ivar)%p_Ddata(ii)-&
                    DcoefficientsAtEdge(ijpos,1,idx)
                  K(jvar,ivar)%p_Ddata(jj) = K(jvar,ivar)%p_Ddata(jj)-&
                      DcoefficientsAtEdge(ijpos,1,idx)
                  K(jvar,ivar)%p_Ddata(ij) = K(jvar,ivar)%p_Ddata(ij)+&
                      DcoefficientsAtEdge(ijpos,2,idx) + DcoefficientsAtEdge(ijpos,1,idx) 
                  K(jvar,ivar)%p_Ddata(ji) = K(jvar,ivar)%p_Ddata(ji)+&
                      DcoefficientsAtEdge(ijpos,3,idx) + DcoefficientsAtEdge(ijpos,1,idx) 
                end do
              end do
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcoefficientsAtEdge)
      !$omp end parallel

    end subroutine doOperatorMat79

  end subroutine gfsys_buildOperatorBlock

  ! ****************************************************************************

!<subroutine>

  subroutine gfsys_buildOperatorScalar(rafcstab, rx,&
      fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, dscale,&
      bclear, rdivOp, rcollection, fcb_calcDivOperator)

!<description>
    ! This subroutine assembles the discrete divergence operator which results
    ! from the group finite element formulation of the continuous problem
    !
    !   <tex> $$ \nabla\cdot{\bf F}(u) $$ </tex>
    !
    ! where ${\bf f}(U)$ is a user-defined flux function for the
    ! multi-component field $U$.
    !
    ! This routine can be used to apply the following discretisations:
    !
    ! (1) the standard Galerkin finite element method
    !     which will be referred to as high-order approximation
    !
    ! (2) discrete upwinding for hyperbolic systems which results from
    !     the conservative elimination of negative eigenvalues from
    !     the Galerkin operator. This technique is for instance
    !     described in the reference:
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
    ! The solution vector
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

    ! OPTIONAL: callback  function to overwrite the standard operation
    include 'intf_calcDivOperator.inc'
    optional :: fcb_calcDivOperator
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! The divergence operator
    type(t_matrixScalar), intent(inout) :: rdivOp

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rxBlock
    type(t_matrixBlock) :: rdivOpBlock
    real(DP), dimension(:), pointer :: p_DivOp,p_Dx
    real(DP), dimension(:,:), pointer :: p_DmatrixCoeffsAtNode
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_Kdiagonal,p_IedgeListIdx


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorScalar')
      call sys_halt()
    end if

    ! Check if user-defined assembly is provided
    if (present(fcb_calcDivOperator)) then
      ! Create auxiliary 1-block vector and matrix
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createMatFromScalar(rdivOp, rdivOpBlock)
      
      ! Call user-defined assembly
      call fcb_calcDivOperator(rafcstab, rxBlock, rdivOpBlock, dscale, bclear,&
          fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, rcollection)
      
      ! Release auxiliary 1-block vector and matrix
      call lsysbl_releaseVector(rxBlock)
      call lsysbl_releaseMatrix(rdivOpBlock)
    else
      ! Set pointers
      call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
      call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
      call afcstab_getbase_DmatCoeffAtNode(rafcstab, p_DmatrixCoeffsAtNode)
      call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
      call lsyssc_getbase_double(rdivOp, p_DivOp)
      call lsyssc_getbase_double(rx, p_Dx)
      
      ! What kind of matrix are we?
      select case(rdivOp%cmatrixFormat)
      case(LSYSSC_MATRIX7INTL, LSYSSC_MATRIX9INTL)
        !-------------------------------------------------------------------------
        ! Matrix format 7 and 9 interleaved
        !-------------------------------------------------------------------------
        
        ! Set diagonal pointer
        if (rdivOp%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
          call lsyssc_getbase_Kld(rdivOp, p_Kdiagonal)
        else
          call lsyssc_getbase_Kdiagonal(rdivOp, p_Kdiagonal)
        end if
        
        ! What type of matrix are we?
        select case(rdivOp%cinterleavematrixFormat)
          
        case (LSYSSC_MATRIX1)
          call doOperatorMat79(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEQ, rdivOp%NA, rx%NVAR, rx%NVAR*rx%NVAR,&
              p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge, p_Dx, dscale, bclear, p_DivOp)
          
        case (LSYSSC_MATRIXD)
          call doOperatorMat79(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEQ, rdivOp%NA, rx%NVAR, rx%NVAR,&
              p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge, p_Dx, dscale, bclear, p_DivOp)
          
        case DEFAULT
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorScalar')
          call sys_halt()
        end select
        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildOperatorScalar')
        call sys_halt()
      end select
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence operator K
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79(Kdiagonal, IedgeListIdx, IedgeList,&
        NEQ, NA, NVAR, MVAR, DmatrixCoeffsAtNode, DmatrixCoeffsAtEdge, Dx,&
        dscale, bclear, K)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal,IedgeListIdx
      integer, intent(in) :: NEQ,NA,NVAR,MVAR

      ! input/output parameters
      real(DP), dimension(MVAR,NA), intent(inout) :: K

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtNode
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtEdge
      integer, dimension(:,:), pointer  :: IdofsAtNode
      
      ! local variables
      integer :: igroup,idx,IEQset,IEQmax,IEDGEset,IEDGEmax
      integer :: i,ii,jj,ij,ji,iedge

      
      !$omp parallel default(shared)&
      !$omp private(DcoefficientsAtEdge,DcoefficientsAtNode,DdataAtEdge,&
      !$omp         DdataAtNode,IEDGEmax,IEQmax,IdofsAtNode,i,idx,&
      !$omp         iedge,ii,ij,ji,jj)

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(IdofsAtNode(2,GFSYS_NEQSIM))
      allocate(DdataAtNode(NVAR,GFSYS_NEQSIM))
      allocate(DcoefficientsAtNode(MVAR,1,GFSYS_NEQSIM))

      ! Loop over the equations
      !$omp do schedule(static,1)
      do IEQset = 1, NEQ, GFSYS_NEQSIM

        ! We always handle GFSYS_NEQSIM equations simultaneously.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most GFSYS_NEQSIM equations simultaneously.
        
        IEQmax = min(NEQ, IEQset-1+GFSYS_NEQSIM)
        
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
          DdataAtNode(:,idx)     = Dx(:,i)
        end do

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DmatrixCoeffsAtNode(:,IEQset:IEQmax),&
            IdofsAtNode(:,1:IEQmax-IEQset+1),&
            dscale, IEQmax-IEQset+1,&
            DcoefficientsAtNode(:,:,1:IEQmax-IEQset+1), rcollection)

        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          do idx = 1, IEQmax-IEQset+1

            ! Get position of diagonal entry
            ii = IdofsAtNode(2,idx)
            
            ! Update the diagonal coefficient
            K(:,ii) = DcoefficientsAtNode(:,1,idx)
          end do

        else   ! do not clear matrix
          do idx = 1, IEQmax-IEQset+1

            ! Get position of diagonal entry
            ii = IdofsAtNode(2,idx)
            
            ! Update the diagonal coefficient
            K(:,ii) = K(:,ii) + DcoefficientsAtNode(:,1,idx)
          end do

        end if
      end do
      !$omp end do

      ! Deallocate temporal memory
      deallocate(IdofsAtNode)
      deallocate(DdataAtNode)
      deallocate(DcoefficientsAtNode)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcoefficientsAtEdge(MVAR,3,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+GFSYS_NEDGESIM)

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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DcoefficientsAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)
          
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
              K(:,ii) = K(:,ii) - DcoefficientsAtEdge(:,1,idx)
              K(:,jj) = K(:,jj) - DcoefficientsAtEdge(:,1,idx)
              K(:,ij) = DcoefficientsAtEdge(:,2,idx) + DcoefficientsAtEdge(:,1,idx) 
              K(:,ji) = DcoefficientsAtEdge(:,3,idx) + DcoefficientsAtEdge(:,1,idx) 
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
              K(:,ii) = K(:,ii) - DcoefficientsAtEdge(:,1,idx)
              K(:,jj) = K(:,jj) - DcoefficientsAtEdge(:,1,idx)
              K(:,ij) = K(:,ij) + DcoefficientsAtEdge(:,2,idx) + DcoefficientsAtEdge(:,1,idx) 
              K(:,ji) = K(:,ji) + DcoefficientsAtEdge(:,3,idx) + DcoefficientsAtEdge(:,1,idx) 
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcoefficientsAtEdge)
      !$omp end parallel

    end subroutine doOperatorMat79

  end subroutine gfsys_buildOperatorScalar

  ! ****************************************************************************

!<subroutine>

  subroutine gfsys_buildVectorBlock(rafcstab, rx, fcb_calcFlux_sim,&
      dscale, bclear, ry, rcollection, fcb_calcDivVector)

!<description>
    ! This subroutine assembles the divergence vector for block vectors.
    ! If the vector contains only one block, then the scalar
    ! counterpart of this routine is called with the scalar subvector.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! callback function to compute local fluxes
    include 'intf_calcFlux_sim.inc'

    ! OPTIONAL: callback  function to overwrite the standard operation
    include 'intf_calcDivVector.inc'
    optional :: fcb_calcDivVector
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx


    ! Check if block vectors contain only one block and no
    ! user-defined callback function is provided.
    if (.not.present(fcb_calcDivVector) .and.&
        (rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1) ) then
      call gfsys_buildVectorScalar(rafcstab, rx%RvectorBlock(1),&
          fcb_calcFlux_sim, dscale, bclear, ry%RvectorBlock(1), rcollection)
      return
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorBlock')
      call sys_halt()
    end if

    ! Check if user-defined assembly is provided
    if (present(fcb_calcDivVector)) then
      ! Call used-defined assembly
      call fcb_calcDivVector(rafcstab, rx, ry, dscale, bclear,&
          fcb_calcFlux_sim, rcollection)
    else
      ! Set pointers
      call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
      call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
      call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
      call lsysbl_getbase_double(rx, p_Dx)
      call lsysbl_getbase_double(ry, p_Dy)
      
      ! Clear vector?
      if (bclear) call lsysbl_clearVector(ry)

      ! Assemble the divergence vector
      call doDivVector(p_IedgeListIdx, p_IedgeList, rafcstab%NEQ,&
          rx%nblocks, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dy)
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector
    
    subroutine doDivVector(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtEdge, Dx, dscale, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DfluxesAtEdge
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j


      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DfluxesAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
            Dy(i,:) = Dy(i,:)+DfluxesAtEdge(:,1,idx)
            Dy(j,:) = Dy(j,:)+DfluxesAtEdge(:,2,idx)
          end do
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      !$omp end parallel
      
    end subroutine doDivVector
    
  end subroutine gfsys_buildVectorBlock

  ! ****************************************************************************

!<subroutine>

  subroutine gfsys_buildVectorScalar(rafcstab, rx, fcb_calcFlux_sim,&
      dscale, bclear, ry, rcollection, fcb_calcDivVector)

!<description>
    ! This subroutine assembles the divergence vector. Note that the
    ! vectors are required as scalar vectors which are stored in the
    ! interleave format.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute local fluxes
    include 'intf_calcFlux_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcDivVector.inc'
    optional :: fcb_calcDivVector
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorScalar), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rxBlock,ryBlock
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildVectorScalar')
      call sys_halt()
    end if

    ! Check if user-defined assembly is provided
    if (present(fcb_calcDivVector)) then
      ! Create auxiliary 1-block vectors
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createVecFromScalar(ry, ryBlock)

      ! Call user-defined assembly
      call fcb_calcDivVector(rafcstab, rxBlock, ryBlock, dscale, bclear,&
          fcb_calcFlux_sim, rcollection)

      ! Release auxiliary 1-block vectors
      call lsysbl_releaseVector(rxBlock)
      call lsysbl_releaseVector(ryBlock)
    else
      ! Set pointers
      call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
      call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
      call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
      call lsyssc_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(ry, p_Dy)

      ! Clear vector?
      if (bclear) call lsyssc_clearVector(ry)
      
      ! Assemble the divergence vector
      call doDivVector(p_IedgeListIdx, p_IedgeList, rafcstab%NEQ,&
          rx%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dy)
    end if
    
  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector

    subroutine doDivVector(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtEdge, Dx, dscale, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DfluxesAtEdge
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DfluxesAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
            Dy(:,i) = Dy(:,i)+DfluxesAtEdge(:,1,idx)
            Dy(:,j) = Dy(:,j)+DfluxesAtEdge(:,2,idx)
          end do
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      !$omp end parallel
      
    end subroutine doDivVector

  end subroutine gfsys_buildVectorScalar
 
end module groupfemsystem
