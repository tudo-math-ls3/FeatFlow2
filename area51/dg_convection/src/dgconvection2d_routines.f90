module dgconvection2d_routines

    use fsystem
    use storage
    use triangulation
    use spatialdiscretisation
    use linearsystemscalar
    use linearsystemblock
    use boundary
    use bilinearformevaluation
    use genoutput
    use scalarpde
    use element
    use cubature
    use basicgeometry
    use transformation
    use dofmapping
    use elementpreprocessing
    use derivatives
    

    implicit none
    
    type t_array
		! Pointer to the double-valued matrix or vector data
    	real(DP), dimension(:), pointer :: Da
	end type t_array
	
	
	integer, parameter				:: nvar2d = 3
	

contains

        
    



!
!
!  !****************************************************************************
!
!!<subroutine>
!
!  subroutine assembleMat9DGBoundaryTerms2D (rmatrixAssembly, rmatrix,&
!      rboundaryRegion, IelementList, IelementOrientation, DedgePosition,&
!      cconstrType, fcoeff_buildMatrixScBdr2D_sim, rcollection)
!
!!<description>
!
!  ! Assembles the matrix entries for a submesh by integrating over the boundary region.
!
!!</description>
!
!!<input>
!
!  ! A boundary region where to assemble the contribution
!  type(t_boundaryRegion), intent(in) :: rboundaryRegion
!
!  ! List of elements where to assemble the bilinear form.
!  integer, dimension(:), intent(in), target :: IelementList
!
!  ! List of element orientations where to assemble the linear form.
!  integer, dimension(:), intent(in) :: IelementOrientation
!
!  ! List of start- and end-parameter values of the edges on the boundary
!  real(DP), dimension(:,:), intent(in) :: DedgePosition
!
!  ! One of the BILF_MATC_xxxx constants that allow to specify the
!  ! matrix construction method.
!  integer, intent(in) :: cconstrType
!
!  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
!  ! Must be present if the matrix has nonconstant coefficients!
!  include 'intf_coefficientMatrixScBdr2D.inc'
!  optional :: fcoeff_buildMatrixScBdr2D_sim
!
!!</input>
!
!!<inputoutput>
!
!  ! A matrix assembly structure prepared with bilf_initAssembly.
!  type(t_bilfMatrixAssembly), intent(inout), target :: rmatrixAssembly
!
!  ! A matrix where to assemble the contributions to.
!  type(t_matrixScalar), intent(inout) :: rmatrix
!
!  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
!  ! callback function for nonconstant coefficients to provide additional
!  ! information.
!  type(t_collection), intent(inout), target, optional :: rcollection
!
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables, used by all processors
!    real(DP), dimension(:), pointer :: p_DA
!    integer :: indofTest,indofTrial,ncubp
!
!    ! local data of every processor when using OpenMP
!    integer :: IELset,IELmax,ibdc,k
!    integer :: iel,icubp,ialbet,ia,ib,idofe,jdofe
!    real(DP) :: domega,daux,db,dlen
!    integer(I32) :: cevaluationTag
!    type(t_bilfMatrixAssembly), target :: rlocalMatrixAssembly
!    type(t_domainIntSubset) :: rintSubset
!    integer, dimension(:,:,:), pointer :: p_Kentry
!integer, dimension(:,:,:), pointer :: p_Kentryu1v1
!integer, dimension(:,:,:), pointer :: p_Kentryu1v2
!integer, dimension(:,:,:), pointer :: p_Kentryu2v1
!integer, dimension(:,:,:), pointer :: p_Kentryu2v2
!    real(DP), dimension(:,:,:), pointer :: p_Dentry
!real(DP), dimension(:,:,:), pointer :: p_Dentryu1v1
!real(DP), dimension(:,:,:), pointer :: p_Dentryu1v2
!real(DP), dimension(:,:,:), pointer :: p_Dentryu2v1
!real(DP), dimension(:,:,:), pointer :: p_Dentryu2v2
!    real(DP), dimension(:), pointer :: p_Domega
!    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest, p_DbasTest_1, p_DbasTest_2
!    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial, p_DbasTrial_1, p_DbasTrial_2
!    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
!    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
!    real(DP), dimension(:), pointer :: p_DcoefficientsBilf
!    integer, dimension(:,:), pointer :: p_IdofsTest
!integer, dimension(:,:), pointer :: p_IdofsTest_1, p_IdofsTest_2
!    integer, dimension(:,:), pointer :: p_IdofsTrial
!integer, dimension(:,:), pointer :: p_IdofsTrial_1, p_IdofsTrial_2
!    type(t_evalElementSet), pointer :: p_revalElementSet,p_revalElementSet_1,p_revalElementSet_2
!type(t_evalElementSet),target::r_revalElementSet_1,r_revalElementSet_2
!    integer, dimension(:,:),pointer :: p_Idescriptors
!
!    ! Arrays for cubature points 1D->2D
!    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
!    real(DP), dimension(:,:,:), allocatable :: Dxi2D_1,Dxi2D_2,DpointsRef_1,DpointsRef_2,DpointsRef,Dxi2D
!real(DP), dimension(:,:,:), allocatable, target :: Dentryu1v1, Dentryu1v2, Dentryu2v1, Dentryu2v2
!integer, dimension(:,:,:), allocatable, target :: Kentryu1v1, Kentryu1v2, Kentryu2v1, Kentryu2v2
!    real(DP), dimension(:,:), allocatable :: DpointsPar
!
!    integer(i32) :: icoordSystem
!
!INTEGER, DIMENSION(:), ALLOCATABLE :: allelements,iedgeset,ielementset_1,ielementset_2
!INTEGER, DIMENSION(:), ALLOCATABLE :: IelementOrientation_1, IelementOrientation_2
!
!integer, dimension(:,:), pointer :: p_IverticesAtEdge
!integer, dimension(:,:), pointer :: p_IelementsAtEdge, p_IedgesAtElement
!real(dp), dimension(:,:), pointer :: p_DvertexCoords
!
!real(DP), dimension(:,:,:,:), allocatable, target :: DbasTest_1, DbasTest_2
!real(DP), dimension(:,:,:,:), allocatable, target :: DbasTrial_1, DbasTrial_2
!
!real(dp) :: dxl1, dxl2, dxr1, dxr2, dyl1, dyl2, dyr1, dyr2, db_1, db_2
!
!    ! Boundary component?
!    ibdc = rboundaryRegion%iboundCompIdx
!
!    ! Get some pointers for faster access
!    call lsyssc_getbase_double (rmatrix,p_DA)
!    indofTest = rmatrixAssembly%indofTest
!    indofTrial = rmatrixAssembly%indofTrial
!    ncubp = rmatrixAssembly%ncubp
!
!    ! Copy the matrix assembly data to the local matrix assembly data,
!    ! where we can allocate memory.
!    ! For single processor machines, this is actually boring and nonsense.
!    ! But using OpenMP, here we get a local copy of the matrix
!    ! assembly structure to where we can add some local data which
!    ! is released upon return without changing the original matrix assembly
!    ! stucture or disturbing the data of the other processors.
!    rlocalMatrixAssembly = rmatrixAssembly
!    call bilf_allocAssemblyData(rlocalMatrixAssembly)
!
!    ! Get some more pointers to local data.
!    p_Kentry => rlocalMatrixAssembly%p_Kentry
!    p_Dentry => rlocalMatrixAssembly%p_Dentry
!    p_Domega => rlocalMatrixAssembly%p_Domega
!    p_DbasTest => rlocalMatrixAssembly%p_DbasTest
!    p_DbasTrial => rlocalMatrixAssembly%p_DbasTrial
!    p_Dcoefficients => rlocalMatrixAssembly%p_Dcoefficients
!    p_DcubPtsRef => rlocalMatrixAssembly%p_DcubPtsRef
!    p_Idescriptors => rlocalMatrixAssembly%rform%Idescriptors
!    p_IdofsTest => rlocalMatrixAssembly%p_IdofsTest
!    p_IdofsTrial => rlocalMatrixAssembly%p_IdofsTrial
!    p_revalElementSet => rlocalMatrixAssembly%revalElementSet
!    p_DcoefficientsBilf => rlocalMatrixAssembly%rform%Dcoefficients
!
!r_revalElementSet_1 = p_revalElementSet
!r_revalElementSet_2 = p_revalElementSet
!p_revalElementSet_1 => r_revalElementSet_1
!p_revalElementSet_2 => r_revalElementSet_2
!
!
!
!
!
!
!
!
!
!allocate(p_DbasTest_1(rmatrixAssembly%indofTest,&
!             elem_getMaxDerivative(rmatrixAssembly%celementTest),&
!             rmatrixAssembly%ncubp,rmatrixAssembly%nelementsPerBlock))
!
!allocate(p_DbasTest_2(rmatrixAssembly%indofTest,&
!             elem_getMaxDerivative(rmatrixAssembly%celementTest),&
!             rmatrixAssembly%ncubp,rmatrixAssembly%nelementsPerBlock))
!
!    ! Allocate memory for the DOF's of all the elements.
!allocate(p_IdofsTest_1(&
!        rmatrixAssembly%indofTest,rmatrixAssembly%nelementsPerBlock))
!allocate(p_IdofsTest_2(&
!        rmatrixAssembly%indofTest,rmatrixAssembly%nelementsPerBlock))
!
!
!    ! the same for the trial basis functions -- if this is not the same FE space.
!    if (rmatrixAssembly%bIdenticalTrialAndTest) then
!      p_DbasTrial_1 => p_DbasTest_1
!      p_IdofsTrial_1 => p_IdofsTest_1
!p_DbasTrial_2 => p_DbasTest_2
!      p_IdofsTrial_2 => p_IdofsTest_2
!    else
!      allocate(p_DbasTrial_1(rmatrixAssembly%indofTrial,&
!              elem_getMaxDerivative(rmatrixAssembly%celementTrial), &
!              rmatrixAssembly%ncubp,rmatrixAssembly%nelementsPerBlock))
!      allocate(p_IdofsTrial_1(&
!          rmatrixAssembly%indofTrial,rmatrixAssembly%nelementsPerBlock))
!allocate(p_DbasTrial_2(rmatrixAssembly%indofTrial,&
!              elem_getMaxDerivative(rmatrixAssembly%celementTrial), &
!              rmatrixAssembly%ncubp,rmatrixAssembly%nelementsPerBlock))
!      allocate(p_IdofsTrial_2(&
!          rmatrixAssembly%indofTrial,rmatrixAssembly%nelementsPerBlock))
!    end if
!
!    ! Allocate an array saving the local matrices for all elements
!    ! in an element set.
!    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
!    ! for this local matrix, but this would normally not fit to the cache
!    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
!    allocate(rmatrixAssembly%p_Kentry(rmatrixAssembly%indofTrial,&
!        rmatrixAssembly%indofTest,rmatrixAssembly%nelementsPerBlock))
!    allocate(rmatrixAssembly%p_Dentry(rmatrixAssembly%indofTrial,&
!        rmatrixAssembly%indofTest,rmatrixAssembly%nelementsPerBlock))
!
!
!
!
!
!
!
!
!
!
!
!
!allocate(Dentryu1v1(ubound(p_Dentry,1),ubound(p_Dentry,2),ubound(p_Dentry,3)))
!allocate(Dentryu1v2(ubound(p_Dentry,1),ubound(p_Dentry,2),ubound(p_Dentry,3)))
!allocate(Dentryu2v1(ubound(p_Dentry,1),ubound(p_Dentry,2),ubound(p_Dentry,3)))
!allocate(Dentryu2v2(ubound(p_Dentry,1),ubound(p_Dentry,2),ubound(p_Dentry,3)))
!allocate(Kentryu1v1(ubound(p_Dentry,1),ubound(p_Dentry,2),ubound(p_Dentry,3)))
!allocate(Kentryu1v2(ubound(p_Dentry,1),ubound(p_Dentry,2),ubound(p_Dentry,3)))
!allocate(Kentryu2v1(ubound(p_Dentry,1),ubound(p_Dentry,2),ubound(p_Dentry,3)))
!allocate(Kentryu2v2(ubound(p_Dentry,1),ubound(p_Dentry,2),ubound(p_Dentry,3)))
!
!p_Dentryu1v1 => Dentryu1v1
!p_Dentryu1v2 => Dentryu1v2
!p_Dentryu2v1 => Dentryu2v1
!p_Dentryu2v2 => Dentryu2v2
!p_Kentryu1v1 => Kentryu1v1
!p_Kentryu1v2 => Kentryu1v2
!p_Kentryu2v1 => Kentryu2v1
!p_Kentryu2v2 => Kentryu2v2
!
!
!CALL storage_getbase_int2D(rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
!CALL storage_getbase_double2D(rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
!CALL storage_getbase_int2D(rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_IelementsAtEdge,p_IelementsAtEdge)
!CALL storage_getbase_int2D(rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
!
!    ! Transpose the coordinate array such that we get coordinates we
!    ! can work with in the mapping between 1D and 2D.
!    do k = 1, ubound(p_DcubPtsRef,1)
!      do icubp = 1,ncubp
!        Dxi1D(icubp,k) = p_DcubPtsRef(k,icubp)
!      end do
!    end do
!
!    ! Allocate memory for the cubature points in 2D.
!    allocate(Dxi2D_1(ncubp,NDIM2D+1,rlocalMatrixAssembly%nelementsPerBlock))
!	allocate(Dxi2D_2(ncubp,NDIM2D+1,rlocalMatrixAssembly%nelementsPerBlock))
!
!    ! Allocate memory for the coordinates of the reference points
!    allocate(DpointsRef_1(NDIM2D+1,ncubp,rlocalMatrixAssembly%nelementsPerBlock))
!    allocate(DpointsRef_2(NDIM2D+1,ncubp,rlocalMatrixAssembly%nelementsPerBlock))
!
!    ! Allocate memory for the parameter values of the points on the boundary
!    allocate(DpointsPar(ncubp,rlocalMatrixAssembly%nelementsPerBlock))
!
!allocate(IelementOrientation_1(rlocalMatrixAssembly%nelementsPerBlock))
!allocate(IelementOrientation_2(rlocalMatrixAssembly%nelementsPerBlock))
!
!
!
!allocate(iedgeset(rlocalMatrixAssembly%nelementsPerBlock))
!allocate(ielementset_1(rlocalMatrixAssembly%nelementsPerBlock))
!allocate(ielementset_2(rlocalMatrixAssembly%nelementsPerBlock))
!allocate(allelements(2*rlocalMatrixAssembly%nelementsPerBlock))
!
!    ! Get the type of coordinate system
!    icoordSystem = elem_igetCoordSystem(rlocalMatrixAssembly%celementTrial)
!
!    ! Loop over the elements - blockwise.
!    !
!    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
!    ! so nelementsPerBlock local matrices are simultaneously calculated in the
!    ! inner loop(s).
!    ! The blocks have all the same size, so we can use static scheduling.
!    !
!    !%OMP do schedule(static,1)
!    do IELset = 1, size(IelementList), rlocalMatrixAssembly%nelementsPerBlock
!!******** IELset durch IEdgeset ersetzen, Ielementlist durch IedgeList
!
!      ! We always handle nelementsPerBlock elements simultaneously.
!      ! How many elements have we actually here?
!      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
!      ! elements simultaneously.
!
!      IELmax = min(size(IelementList),IELset-1+rlocalMatrixAssembly%nelementsPerBlock)
!
!	  iedgeset(:)=IelementList(IelSet:IelMax)
!	  do iel = 1,IELmax-IELset+1
!        ielementset_1(iel)=p_IelementsAtEdge(1,iedgeset(iel))
!        ielementset_2(iel)=p_IelementsAtEdge(2,iedgeset(iel))
!      end do
!
!allelements(1:iel)=ielementset_1
!allelements(iel+1:2*iel)=ielementset_2
!
!CALL bilf_findEdgeOrientation2d (ielementset_1,iedgeset,p_IedgesAtElement,IelementOrientation_1)
!CALL bilf_findEdgeOrientation2d (ielementset_2,iedgeset,p_IedgesAtElement,IelementOrientation_2)
!
!
!      ! Map the 1D cubature points to the edges in 2D.
!      do iel = 1,IELmax-IELset+1
!        call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation_1(iel), &
!            ncubp, Dxi1D, Dxi2D_1(:,:,iel))
!		call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation_2(iel), &
!            ncubp, Dxi1D, Dxi2D_2(:,:,iel))
!      end do
!
!!       ! Calculate the parameter values of the points
!!       do iel = 1,IELmax-IELset+1
!!         do icubp = 1,ncubp
!!           ! Dxi1D is in [-1,1] while the current edge has parmeter values
!!           ! [DedgePosition(1),DedgePosition(2)]. So do a linear
!!           ! transformation to transform Dxi1D into that interval, this
!!           ! gives the parameter values in length parametrisation
!!           call mprim_linearRescale(Dxi1D(icubp,1), -1.0_DP, 1.0_DP,&
!!               DedgePosition(1,IELset+iel-1), DedgePosition(2,IELset+iel-1),&
!!               DpointsPar(icubp,iel))
!!         end do
!!       end do
!
!      ! Transpose the coordinate array such that we get coordinates we
!      ! can work with.
!      do iel = 1,IELmax-IELset+1
!        do icubp = 1,ncubp
!          do k = 1,ubound(DpointsRef_1,1)
!            DpointsRef_1(k,icubp,iel) = Dxi2D_1(icubp,k,iel)
!            DpointsRef_2(k,icubp,iel) = Dxi2D_2(icubp,k,iel)
!          end do
!        end do
!      end do
!
!      ! --------------------- DOF SEARCH PHASE ------------------------
!
!      ! The outstanding feature with finite elements is: A basis
!      ! function for a DOF on one element has common support only
!      ! with the DOF's on the same element! E.g. for Q1:
!      !
!      !        #. . .#. . .#. . .#
!      !        .     .     .     .
!      !        .  *  .  *  .  *  .
!      !        #-----O-----O. . .#
!      !        |     |     |     .
!      !        |     | iel |  *  .
!      !        #-----X-----O. . .#
!      !        |     |     |     .
!      !        |     |     |  *  .
!      !        #-----#-----#. . .#
!      !
!      ! --> On element iel, the basis function at "X" only interacts
!      !     with the basis functions in "O". Elements in the
!      !     neighbourhood ("*") have no support, therefore we only have
!      !     to collect all "O" DOF's.
!      !
!      ! Calculate the global DOF's into IdofsTrial / IdofsTest.
!      !
!      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
!      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
!      call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTest, &
!          ielementset_1, p_IdofsTest_1)
!	  call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTest, &
!          ielementset_2, p_IdofsTest_2)
!
!      ! If the DOF's for the test functions are different, calculate them, too.
!      if (.not. rlocalMatrixAssembly%bIdenticalTrialAndTest) then
!        call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTrial, &
!            ielementset_1, p_IdofsTrial_1)
!        call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTrial, &
!            ielementset_2, p_IdofsTrial_2)
!      end if
!
!      ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
!
!      ! For the assembly of the global matrix, we use a "local"
!      ! approach. At first we build a "local" system matrix according
!      ! to the current element. This contains all additive
!      ! contributions of element iel, which are later added at the
!      ! right positions to the elements in the global system matrix.
!      !
!      ! We have indofTrial trial DOF's per element and
!      ! indofTest test DOF's per element. Therefore there are
!      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j)
!      ! "active" (i.e. have common support) on our current element, each
!      ! giving an additive contribution to the system matrix.
!      !
!      ! We build a quadratic indofTrial*indofTest local matrix:
!      ! Kentry(1..indofTrial,1..indofTest) receives the position
!      ! in the global system matrix, where the corresponding value
!      ! has to be added to.
!      ! (The corresponding contributions can be saved separately,
!      ! but we directly add them to the global matrix in this
!      ! approach.)
!      !
!      ! We build local matrices for all our elements
!      ! in the set simultaneously. Get the positions of the local matrices
!      ! in the global matrix.
!     call bilf_getLocalMatrixIndices (rmatrix,p_IdofsTrial,p_IdofsTest,p_Kentry,&
!          ubound(p_IdofsTrial,1),ubound(p_IdofsTest,1),IELmax-IELset+1)
!
!      call bilf_getLocalMatrixIndices (rmatrix,p_IdofsTrial_1,p_IdofsTest_1,p_Kentryu1v1,&
!          ubound(p_IdofsTrial_1,1),ubound(p_IdofsTest_1,1),IELmax-IELset+1)
!      call bilf_getLocalMatrixIndices (rmatrix,p_IdofsTrial_1,p_IdofsTest_2,p_Kentryu1v2,&
!          ubound(p_IdofsTrial_1,1),ubound(p_IdofsTest_2,1),IELmax-IELset+1)
!      call bilf_getLocalMatrixIndices (rmatrix,p_IdofsTrial_2,p_IdofsTest_1,p_Kentryu2v1,&
!          ubound(p_IdofsTrial_2,1),ubound(p_IdofsTest_1,1),IELmax-IELset+1)
!      call bilf_getLocalMatrixIndices (rmatrix,p_IdofsTrial_2,p_IdofsTest_2,p_Kentryu2v2,&
!          ubound(p_IdofsTrial_2,1),ubound(p_IdofsTest_2,1),IELmax-IELset+1)
!
!      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
!
!      ! Ok, we found the positions of the local matrix entries
!      ! that we have to change.
!      ! To calculate the matrix contributions, we have to evaluate
!      ! the elements to give us the values of the basis functions
!      ! in all the DOF's in all the elements in our set.
!
!      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
!      ! the elements later. All of them can be combined with OR, what will give
!      ! a combined evaluation tag.
!      cevaluationTag = rlocalMatrixAssembly%cevaluationTag
!
!      ! The cubature points are already initialised by 1D->2D mapping.
!      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
!
!      ! Calculate all information that is necessary to evaluate the finite element
!      ! on all cells of our subset. This includes the coordinates of the points
!      ! on the cells.
!      call elprep_prepareSetForEvaluation (p_revalElementSet,&
!          cevaluationTag, rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
!          IelementList(IELset:IELmax), rlocalMatrixAssembly%ctrafoType, &
!          DpointsRef=DpointsRef)
!call elprep_prepareSetForEvaluation (p_revalElementSet_1,&
!          cevaluationTag, rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
!          Ielementset_1(:), rlocalMatrixAssembly%ctrafoType, &
!          DpointsRef=DpointsRef_1)
!call elprep_prepareSetForEvaluation (p_revalElementSet_2,&
!          cevaluationTag, rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
!          Ielementset_2(:), rlocalMatrixAssembly%ctrafoType, &
!          DpointsRef=DpointsRef_2)
!
!!!!!!!!!!!!!!!!!!!!!!!!ielementlist + _1 _2
!!!!!!!!!!!!!!!!!!!!!!!!dDpointsref + _1 _2
!
!
!
!      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
!      if (.not. rlocalMatrixAssembly%rform%ballCoeffConstant) then
!        if (present(fcoeff_buildMatrixScBdr2d_sim)) then
!          call domint_initIntegrationByEvalSet (p_revalElementSet,rintSubset)
!          rintSubset%ielementDistribution = 0
!          rintSubset%ielementStartIdx = IELset
!          rintSubset%p_Ielements => IelementList(IELset:IELmax)
!          rintSubset%p_IdofsTrial => p_IdofsTrial
!          rintSubset%celement = rlocalMatrixAssembly%celementTrial
!          call fcoeff_buildMatrixScBdr2D_sim (rmatrix%p_rspatialDiscrTest,&
!              rmatrix%p_rspatialDiscrTrial,&
!              rlocalMatrixAssembly%rform, IELmax-IELset+1, ncubp,&
!              p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!              ibdc, DpointsPar(:,1:IELmax-IELset+1),&
!              p_IdofsTrial, p_IdofsTest, rintSubset, &
!              p_Dcoefficients(:,:,1:IELmax-IELset+1), rcollection)
!          call domint_doneIntegration (rintSubset)
!        else
!          p_Dcoefficients(:,:,1:IELmax-IELset+1) = 1.0_DP
!        end if
!      end if
!
!      ! Calculate the values of the basis functions.
!      call elem_generic_sim2 (rlocalMatrixAssembly%celementTest, &
!          p_revalElementSet, rlocalMatrixAssembly%BderTest, &
!          rlocalMatrixAssembly%p_DbasTest)
!
!call elem_generic_sim2 (rlocalMatrixAssembly%celementTest, &
!          p_revalElementSet_1, rlocalMatrixAssembly%BderTest, &
!          p_DbasTest_1)
!call elem_generic_sim2 (rlocalMatrixAssembly%celementTest, &
!          p_revalElementSet_2, rlocalMatrixAssembly%BderTest, &
!          p_DbasTest_2)
!
!      ! Omit the calculation of the trial function values if they
!      ! are identical to the test function values.
!      if (.not. rlocalMatrixAssembly%bidenticalTrialAndTest) then
!        call elem_generic_sim2 (rlocalMatrixAssembly%celementTrial, &
!            p_revalElementSet, rlocalMatrixAssembly%BderTrial, &
!            rlocalMatrixAssembly%p_DbasTrial)
!  call elem_generic_sim2 (rlocalMatrixAssembly%celementTrial, &
!            p_revalElementSet_1, rlocalMatrixAssembly%BderTrial, &
!            p_DbasTrial_1)
!  call elem_generic_sim2 (rlocalMatrixAssembly%celementTrial, &
!            p_revalElementSet_2, rlocalMatrixAssembly%BderTrial, &
!            p_DbasTrial_2)
!
!
!      end if
!
!      ! --------------------- DOF COMBINATION PHASE ------------------------
!
!      ! Values of all basis functions calculated. Now we can start
!      ! to integrate!
!
!      ! Clear the local matrices
!      p_Dentry(:,:,1:IELmax-IELset+1) = 0.0_DP
!p_Dentryu1v1(:,:,1:IELmax-IELset+1) = 0.0_DP
!p_Dentryu2v1(:,:,1:IELmax-IELset+1) = 0.0_DP
!p_Dentryu1v2(:,:,1:IELmax-IELset+1) = 0.0_DP
!p_Dentryu2v2(:,:,1:IELmax-IELset+1) = 0.0_DP
!
!      ! We have two different versions for the integration - one
!      ! with constant coefficients and one with nonconstant coefficients.
!      !
!      ! Check the bilinear form which one to use:
!
!      if (rlocalMatrixAssembly%rform%ballCoeffConstant) then
!
!        ! Constant coefficients. The coefficients are to be found in
!        ! the Dcoefficients variable of the form.
!        !
!        ! Loop over the elements in the current set.
!
!        do iel = 1,IELmax-IELset+1
!
!          ! Get the length of the edge. Let's use the parameter values
!          ! on the boundary for that purpose; this is a more general
!          ! implementation than using simple lines as it will later
!          ! support isoparametric elements.
!          !
!          ! The length of the current edge serves as a "determinant"
!          ! in the cubature, so we have to divide it by 2 as an edge on
!          ! the unit interval [-1,1] has length 2.
!          !dlen = 0.5_DP*(DedgePosition(2,IELset+iel-1)-DedgePosition(1,IELset+iel-1))
!
!
!
!dxl1=p_DvertexCoords(1,p_IverticesAtEdge(1,iel))
!dyl1=p_DvertexCoords(2,p_IverticesAtEdge(1,iel))
!dxl2=p_DvertexCoords(1,p_IverticesAtEdge(2,iel))
!dyl2=p_DvertexCoords(2,p_IverticesAtEdge(2,iel))
!dlen=0.5_DP*sqrt((dxl1-dxl2)*(dxl1-dxl2)+(dyl1-dyl2)*(dyl1-dyl2))
!
!
!
!          ! Loop over all cubature points on the current element
!          do icubp = 1, ncubp
!
!            ! Calculate the current weighting factor in the cubature formula
!            ! in that cubature point.
!
!            domega = dlen * p_Domega(icubp)
!
!            ! Loop over the additive factors in the bilinear form.
!            do ialbet = 1,rlocalMatrixAssembly%rform%itermcount
!
!              ! Get from Idescriptors the type of the derivatives for the
!              ! test and trial functions. The summand we calculate
!              ! here will be added to the matrix entry:
!              !
!              ! a_ij  =  int_... ( psi_j )_ib  *  ( phi_i )_ia
!              !
!              ! -> Ix=0: function value,
!              !      =1: first derivative, ...
!              !    as defined in the module 'derivative'.
!
!              ia = p_Idescriptors(1,ialbet)
!              ib = p_Idescriptors(2,ialbet)
!
!              ! Multiply domega with the coefficient of the form.
!              ! This gives the actual value to multiply the
!              ! function value with before summing up to the integral.
!              daux = domega * p_DcoefficientsBilf(ialbet)
!
!              ! Now loop through all possible combinations of DOF's
!              ! in the current cubature point. The outer loop
!              ! loops through the "O"'s in the above picture,
!              ! the test functions:
!
!              do idofe = 1,indofTest
!
!                ! Get the value of the (test) basis function
!                ! phi_i (our "O") in the cubature point:
!                db = p_DbasTest(idofe,ib,icubp,iel)
!db_1= p_DbasTest_1(idofe,ib,icubp,iel)
!db_2= p_DbasTest_2(idofe,ib,icubp,iel)
!
!                ! Perform an inner loop through the other DOF's
!                ! (the "X").
!
!                do jdofe = 1,indofTrial
!
!                  ! Get the value of the basis function
!                  ! psi_j (our "X") in the cubature point.
!                  ! Them multiply:
!                  !    db * dbas(..) * daux
!                  ! ~= phi_i * psi_j * coefficient * cub.weight
!                  ! Summing this up gives the integral, so the contribution
!                  ! to the global matrix.
!                  !
!                  ! Simply summing up db * dbas(..) * daux would give
!                  ! the coefficient of the local matrix. We save this
!                  ! contribution in the local matrix.
!
!                  !JCOLB = Kentry(jdofe,idofe,iel)
!                  !p_DA(JCOLB) = p_DA(JCOLB) + db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
!                  p_Dentry(jdofe,idofe,iel) = p_Dentry(jdofe,idofe,iel) + &
!                                        db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
!                  p_Dentryu1v1(jdofe,idofe,iel) = p_Dentryu1v1(jdofe,idofe,iel) + &
!                                        db_1*p_DbasTrial_1(jdofe,ia,icubp,iel)*daux
!                  p_Dentryu1v2(jdofe,idofe,iel) = p_Dentryu1v2(jdofe,idofe,iel) + &
!                                        db_2*p_DbasTrial_1(jdofe,ia,icubp,iel)*daux
!                  p_Dentryu2v1(jdofe,idofe,iel) = p_Dentryu2v1(jdofe,idofe,iel) + &
!                                        db_1*p_DbasTrial_2(jdofe,ia,icubp,iel)*daux
!                  p_Dentryu2v2(jdofe,idofe,iel) = p_Dentryu2v2(jdofe,idofe,iel) + &
!                                        db_2*p_DbasTrial_2(jdofe,ia,icubp,iel)*daux
!
!
!                end do ! jdofe
!
!              end do ! idofe
!
!            end do ! ialbet
!
!          end do ! icubp
!
!        end do ! iel
!
!      else
!
!!         ! Nonconstant coefficients. The coefficients are to be found in
!!         ! the Dcoefficients variable as computed above.
!!         !
!!         ! Loop over the elements in the current set.
!!
!!         do iel = 1,IELmax-IELset+1
!!
!!           ! Get the length of the edge. Let's use the parameter values
!!           ! on the boundary for that purpose; this is a more general
!!           ! implementation than using simple lines as it will later
!!           ! support isoparametric elements.
!!           !
!!           ! The length of the current edge serves as a "determinant"
!!           ! in the cubature, so we have to divide it by 2 as an edge on
!!           ! the unit interval [-1,1] has length 2.
!!           dlen = 0.5_DP*(DedgePosition(2,IELset+iel-1)-DedgePosition(1,IELset+iel-1))
!!
!!           ! Loop over all cubature points on the current element
!!           do icubp = 1, ncubp
!!
!!             ! calculate the current weighting factor in the cubature formula
!!             ! in that cubature point.
!!
!!             domega = dlen * p_Domega(icubp)
!!
!!             ! Loop over the additive factors in the bilinear form.
!!             do ialbet = 1,rlocalMatrixAssembly%rform%itermcount
!!
!!               ! Get from Idescriptors the type of the derivatives for the
!!               ! test and trial functions. The summand we calculate
!!               ! here will be added to the matrix entry:
!!               !
!!               ! a_ij  =  int_... ( psi_j )_ia  *  ( phi_i )_ib
!!               !
!!               ! -> Ix=0: function value,
!!               !      =1: first derivative, ...
!!               !    as defined in the module 'derivative'.
!!
!!               ia = rlocalMatrixAssembly%rform%Idescriptors(1,ialbet)
!!               ib = rlocalMatrixAssembly%rform%Idescriptors(2,ialbet)
!!
!!               ! Multiply domega with the coefficient of the form.
!!               ! This gives the actual value to multiply the
!!               ! function value with before summing up to the integral.
!!               ! Get the precalculated coefficient from the coefficient array.
!!               daux = domega * p_Dcoefficients(ialbet,icubp,iel)
!!
!!               ! Now loop through all possible combinations of DOF's
!!               ! in the current cubature point. The outer loop
!!               ! loops through the "O" in the above picture,
!!               ! the test functions:
!!
!!               do idofe = 1,indofTest
!!
!!                 ! Get the value of the (test) basis function
!!                 ! phi_i (our "O") in the cubature point:
!!                 db = p_DbasTest(idofe,ib,icubp,iel)
!!
!!                 ! Perform an inner loop through the other DOF's
!!                 ! (the "X").
!!
!!                 do jdofe = 1,indofTrial
!!
!!                   ! Get the value of the basis function
!!                   ! psi_j (our "X") in the cubature point.
!!                   ! Them multiply:
!!                   !    db * dbas(..) * daux
!!                   ! ~= phi_i * psi_j * coefficient * cub.weight
!!                   ! Summing this up gives the integral, so the contribution
!!                   ! to the global matrix.
!!                   !
!!                   ! Simply summing up db * dbas(..) * daux would give
!!                   ! the coefficient of the local matrix. We save this
!!                   ! contribution in the local matrix of element iel.
!!
!!                   !JCOLB = Kentry(jdofe,idofe,iel)
!!                   !p_DA(JCOLB) = p_DA(JCOLB) + db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
!!                   p_Dentry(jdofe,idofe,iel) = &
!!                       p_Dentry(jdofe,idofe,iel)+db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
!!
!!                 end do
!!
!!               end do ! jdofe
!!
!!             end do ! ialbet
!!
!!           end do ! icubp
!!
!!         end do ! iel
!
!      end if ! rform%ballCoeffConstant
!
!      ! Incorporate the local matrices into the global one.
!      ! Kentry gives the position of the additive contributions in Dentry.
!      !
!      ! OpenMP-Extension: This is a critical section. Only one thread is
!      ! allowed to write to the matrix, otherwise the matrix may get
!      ! messed up.
!      ! The critical section is put around both loops as indofTest/indofTrial
!      ! are usually small and quickly to handle.
!
!      if (cconstrType .eq. BILF_MATC_LUMPED) then
!
!        !%OMP CRITICAL
!        do iel = 1,IELmax-IELset+1
!
!          do idofe = 1,indofTest
!            daux = 0.0_DP
!            do jdofe = 1,indofTrial
!              daux = daux + p_Dentry(jdofe,idofe,iel)
!            end do
!            p_DA(p_Kentry(idofe,idofe,iel)) = &
!                p_DA(p_Kentry(idofe,idofe,iel)) + daux
!          end do
!
!        end do ! iel
!        !%OMP END CRITICAL
!
!      else
!
!        !%OMP CRITICAL
!        do iel = 1,IELmax-IELset+1
!
!          do idofe = 1,indofTest
!            do jdofe = 1,indofTrial
!              p_DA(p_Kentry(jdofe,idofe,iel)) = &
!                  p_DA(p_Kentry(jdofe,idofe,iel)) + p_Dentry(jdofe,idofe,iel)
!
!              p_DA(p_Kentryu1v1(jdofe,idofe,iel)) = &
!                  p_DA(p_Kentryu1v1(jdofe,idofe,iel)) + p_Dentryu1v1(jdofe,idofe,iel)
!              p_DA(p_Kentryu1v2(jdofe,idofe,iel)) = &
!                  p_DA(p_Kentryu1v2(jdofe,idofe,iel)) + p_Dentryu1v2(jdofe,idofe,iel)
!              p_DA(p_Kentryu2v1(jdofe,idofe,iel)) = &
!                  p_DA(p_Kentryu2v1(jdofe,idofe,iel)) + p_Dentryu2v1(jdofe,idofe,iel)
!              p_DA(p_Kentryu2v2(jdofe,idofe,iel)) = &
!                  p_DA(p_Kentryu2v2(jdofe,idofe,iel)) + p_Dentryu2v2(jdofe,idofe,iel)
!            end do
!          end do
!
!        end do ! iel
!        !%OMP END CRITICAL
!
!      end if
!
!    end do ! IELset
!
!    ! Release the local matrix assembly structure
!    call bilf_releaseAssemblyData(rlocalMatrixAssembly)
!
!    ! Deallocate memory
!    !deallocate(Dxi2D, DpointsRef, DpointsPar)
!deallocate(Dxi2D, DpointsRef, iedgeset, allelements,ielementset_1,ielementset_2)
!deallocate(Dentryu1v1)
!deallocate(Dentryu1v2)
!deallocate(Dentryu2v1)
!deallocate(Dentryu2v2)
!deallocate(Kentryu1v1)
!deallocate(Kentryu1v2)
!deallocate(Kentryu2v1)
!deallocate(Kentryu2v2)
!deallocate(IelementOrientation_1)
!deallocate(IelementOrientation_2)
!deallocate(p_DbasTest_1)
!deallocate(p_DbasTest_2)
!deallocate(    p_DbasTrial_1, p_DbasTrial_2)
!    deallocate(p_IdofsTest_1, p_IdofsTest_2)
!    deallocate(p_IdofsTrial_1, p_IdofsTrial_2)
!  end subroutine







  !************************************************************************
  
!<subroutine>

  subroutine generic_sim2_single (celement, revalElementSet, Bder, Dbas)

!<description>
  ! This subroutine simultaneously calculates the values of the basic
  ! functions of the finite element at multiple given points on the reference
  ! element for multiple given elements.
!</description>

!<input>
  ! t_evalElementSet-structure that contains cell-specific information and
  ! coordinates of the evaluation points. revalElementSet must be prepared
  ! for the evaluation.
  type(t_evalElementSet), intent(in) :: revalElementSet

  ! Element type identifier
  integer(I32), intent(in)  :: celement
  
  ! Derivative quantifier array. array [1..DER_MAXNDER] of boolean.
  ! If bder(DER_xxxx)=true, the corresponding derivative (identified
  ! by DER_xxxx) is computed by the element (if supported). Otherwise,
  ! the element might skip the computation of that value type, i.e.
  ! the corresponding value 'Dvalue(DER_xxxx)' is undefined.
  logical, dimension(:), intent(in) :: Bder
!</input>
  
!<output>
  ! Value/derivatives of basis functions.
  ! array [1..EL_MAXNBAS,1..DER_MAXNDER,1..npointsPerElement,nelements] of double
  ! Bder(DER_FUNC)=true  => Dbas(i,DER_FUNC,j) defines the value of the i-th
  !   basis function of the finite element in the point Dcoords(j) on the
  !   reference element,
  !   Dvalue(i,DER_DERIV_X) the value of the x-derivative of the i-th
  !   basis function,...
  ! Bder(DER_xxxx)=false => Dbas(i,DER_xxxx,.) is undefined.
  real(DP), dimension(:,:,:), intent(out) :: Dbas
!</output>

! </subroutine>

  ! local variables
  real(DP), dimension(ubound(Dbas,1)-lbound(Dbas,1)+1,ubound(Dbas,2)-lbound(Dbas,2)+1,ubound(Dbas,3)-lbound(Dbas,3)+1,1) :: Dbas_array

  call elem_generic_sim2 (celement, revalElementSet, Bder, Dbas_array)

  Dbas(:,:,:) = Dbas_array(:,:,:,1)
end subroutine









 !****************************************************************************

  !<subroutine>

  subroutine getLocalMatrixIndices_single (rmatrix,Irows,Icolumns,Kentry,&
      irowsPerElement,icolsPerElement)
  
  !<description>
  
  ! Calculates index positions of local matrices in a global matrix.
  ! For a set of elements, Icolumns and Irows define the row and column indices
  ! of local matrices which have to be accessed in a global matrix rmatrix.
  ! The routine then calculates the positions of all the corresponding matrix
  ! entries in the data array of the matrix rmatrix and saves the result
  ! to the array Kentry.
  !
  ! IMPORTANT: For performance reasons, the columns and rows in Kentry
  ! are saved *transposed* in comparison to the matrix rmatrix!
  ! That means that
  !    Kentry(j,i,:) = position of element (Irows(i,:),Icolumns(j,:))
  ! holds!
  
  !</description>
  
  !<input>
  
  ! The global matrix which has to be accessed.
  type(t_matrixScalar), intent(in) :: rmatrix
  
  ! Array identifying all rows in the global matrix which have to be
  ! accessed.
  ! DIMENSION(#rows per element, #elements).
  integer, dimension(:), intent(in) :: Irows

  ! Array identifying all columns in the global matrix which have to be
  ! accessed.
  ! DIMENSION(#columns per element, #elements).
  integer, dimension(:), intent(in) :: Icolumns
  
  ! Number of rows per element / in the local matrix
  integer, intent(in) :: irowsPerElement
  
  ! Number of columns per element / in the local matrix
  integer, intent(in) :: icolsPerElement
  
  
  !</input>
  
  !<output>
  
  ! Array receiving the positions of the local matrices in the global matrix.
  ! DIMENSION(#columns per element,#rows per element,#elements).
  ! Saved in a transposed way:
  !    Kentry(j,i,:) = position of element (Irows(i,:),Icolumns(j,:))
  integer, dimension(:,:), intent(out) :: Kentry
  
  !</output>
  
  !</subroutine>
  
  
  
  ! Local variables
  
  ! Array identifying all rows in the global matrix which have to be
  ! accessed.
  ! DIMENSION(#rows per element, #elements).
  integer, dimension(size(Irows),1) :: Irows_array

  ! Array identifying all columns in the global matrix which have to be
  ! accessed.
  ! DIMENSION(#columns per element, #elements).
  integer, dimension(size(Icolumns),1) :: Icolumns_array
  
  ! Array receiving the positions of the local matrices in the global matrix.
  ! DIMENSION(#columns per element,#rows per element,#elements).
  ! Saved in a transposed way:
  !    Kentry(j,i,:) = position of element (Irows(i,:),Icolumns(j,:))
  integer, dimension(ubound(Kentry,1)-lbound(Kentry,1)+1,ubound(Kentry,2)-lbound(Kentry,2)+1,1) :: Kentry_array
  
  
  
  Irows_array(:,1)    = Irows
  Icolumns_array(:,1) = Icolumns
  Kentry_array(:,:,1) = Kentry
  
  call bilf_getLocalMatrixIndices (rmatrix,Irows_array,Icolumns_array,Kentry_array,&
      irowsPerElement,icolsPerElement,1)

  Kentry   = Kentry_array(:,:,1)
  
  end subroutine








  !************************************************************************
  
!<subroutine>

  subroutine prepareOneForMultipleEvaluation (revalElementSet, cevaluationTag, &
      rtriangulation, Ielement, ctrafoType, DpointsRef)

!<description>
  ! This subroutine prepares a t_evalElementSet structure to be used for
  ! the evaluation of a finite element in a set of cells.
  ! Dpoints contains a list of coordinates on the reference element where to
  ! evaluate. These points are mapped onto all elements in the list
  ! IelementList from the triangulation rtriangulation.
  ! cevaluationTag specifies an 'evaluation tag' that defines which
  ! information must be prepared by this routine; that tag can be obtained
  ! by asking the finite element what it needs by calling
  ! elem_getEvaluationTag.
!</description>

!<input>
  ! Evaluation tag. This is a bitfield that specifies which information is
  ! prepared in the element set; it is a combination of EL_EVLTAG_XXXX-constants.
  !
  ! Note: If EL_EVLTAG_REFPOINTS is not specified in this tag, the coordinates on
  ! the reference element are assumed to be initialised! Dpoints is ignored
  ! in that case. This can be used to calculate this information only once
  ! in a first call while then using the same set of reference coordinates
  ! for all subsequent calls.
  integer(I32), intent(in) :: cevaluationTag

  ! Underlying triangulation of the domain
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! List of elements in the mesh where the integration is to be performed.
  integer, intent(in) :: Ielement
  
  ! Type of transformation from the reference element to the real element.
  integer(I32), intent(in) :: ctrafoType
  
  ! OPTIONAL: A set of npointsPerElement tuples (x,y) (or (x,y,z) in 3D) of the
  ! points where to evaluate. These coordinates define on the reference element
  ! the coordinates of the cubature points where to evaluate the element.
  ! DIMENSION(ndimension,npointsPerElement)
  ! This array specifies the evaluation points for exactly one element
  ! and can be used if the coordinates of the evaluation points are the
  ! same for all elements. In contrast, DpointsRef may specify different points
  ! on all elements.
  ! Ignored if EL_EVLTAG_REFPOINTS is not specified in cevaluationTag.
  ! real(DP), dimension(:,:), optional :: Dpoints
  
  ! OPTIONAL: An array with coordinates of all points for all the elements
  ! where to evaluate -- relative to the reference element.
  ! For each element, a set of points can be specified here.
  ! A pointer to this array is saved in the revalElementSet structure. The array
  ! is assumed to be maintained by the caller.
  !
  ! If not specified, the routine will automatically set up that array
  ! using Dpoints (i.e. the coordinates of Dpoints are 'distributed' to
  ! all elements described by DpointsRef).
  real(DP), dimension(:,:), target :: DpointsRef
!</input>
  
!<inputoutput>
  ! The element set that is to be initialised. If this is already initialised,
  ! previous information is overwritten.
  type(t_evalElementSet), intent(inout) :: revalElementSet
!</inputoutput>

! </subroutine>


! List of elements in the mesh where the integration is to be performed.
  integer, dimension(1) :: IelementList

  real(DP), dimension(ubound(DpointsRef,1)-lbound(DpointsRef,1)+1, &
                      ubound(DpointsRef,2)-lbound(DpointsRef,2)+1, &
                      1), target :: DpointsRef_array

  Ielementlist(1) = ielement
  DpointsRef_array (:,:,1) = DpointsRef(:,:)

  call elprep_prepareSetForEvaluation (revalElementSet,&
              cevaluationTag, rtriangulation, &
              IelementList, ctrafoType, &
              DpointsRef=DpointsRef_array)
              
  DpointsRef(:,:) = DpointsRef_array (:,:,1)

  end subroutine




  !****************************************************************************
  
!<subroutine>
  
  subroutine findEdgeOrientation2d_list (IelementList,IedgeList,p_IedgesAtElement,IelementOrientation)
  
!<description>

  !

!</description>

!<input>

  ! List of elements on which the local edgenumber has to be found
  integer, dimension(:), intent(in), target :: IelementList

  ! List of edges whose local edgenumber has to be found
  integer, dimension(:), intent(in), target :: Iedgelist

  ! Pointer to the IedgesAtElement of the underlying triangulation
  integer, dimension(:,:), intent(in), pointer :: p_IedgesAtElement
  
!</input>

!<inputoutput>
  
  ! The Output: List of element orientations
  integer, dimension(:), intent(inout) :: IelementOrientation

!</inputoutput>
  
!</subroutine>


    ! local variables
    integer :: i, iel, iedge, ilocaledge

	! First test, if all arrays have the same length or throw an error
	if ((size(IelementList).ne.size(IedgeList)).or.(size(IelementOrientation).ne.size(IedgeList))) then
	      call output_line ('Arrays are not of the same size', &
                        OU_CLASS_ERROR,OU_MODE_STD,'findEdgeOrientation2d')
      call sys_halt()
	end if

	! Loop over all elements in Ielementlist
	do i= 1, size(IelementList)

		! Get global element-number on which to find the edge
		iel = Ielementlist(i)
		
		! Get global edge-number to find on that element
		iedge = IedgeList(i)

		! Find the local edge number
		do ilocaledge = 1, ubound(p_IedgesAtElement,1)
			if (p_IedgesAtElement(ilocaledge,iel) .eq. iedge) exit
		end do ! ilocaledge

		! Set the local edge number in the output array
		Ielementorientation(i) = ilocaledge
	end do ! loop over IelementList

  end subroutine


  !****************************************************************************
  
!<subroutine>
  
  subroutine findEdgeOrientation2d (Ielement,Iedge,p_IedgesAtElement,IelementOrientation)
  
!<description>

  !

!</description>

!<input>

  ! List of elements on which the local edgenumber has to be found
  integer, intent(in) :: Ielement

  ! List of edges whose local edgenumber has to be found
  integer, intent(in):: Iedge

  ! Pointer to the IedgesAtElement of the underlying triangulation
  integer, dimension(:,:), intent(in), pointer :: p_IedgesAtElement
  
!</input>

!<inputoutput>
  
  ! The Output: List of element orientations
  integer, intent(out) :: IelementOrientation

!</inputoutput>
  
!</subroutine>


    ! local variables
    integer :: ilocaledge

		! Find the local edge number
		do ilocaledge = 1, ubound(p_IedgesAtElement,1)
			if (p_IedgesAtElement(ilocaledge,ielement) .eq. iedge) exit
		end do ! ilocaledge

		! Set the local edge number in the output array
		Ielementorientation = ilocaledge
	

  end subroutine
  
  
  
  
  
  !****************************************************************************
  
  !<subroutine>

  subroutine velocityfield (dpoint,dvelocity)
    
  use fsystem
    
  !<description>
    ! This subroutine is called during the matrix assembly of the face terms
    ! The input is the coordinates of the points dpoints
    ! The output is the velocity at this point dvelocity
  !</description>
    
  !<input>
  ! The point, where to evaluate the velocity
  real(DP), dimension(:), intent(IN) :: dpoint
  !</input>
  
  !<output>
  ! The velocity
  real(DP), dimension(size(dpoint)), intent(OUT) :: dvelocity
  !</output>
    
  !</subroutine>
  
  dvelocity(1) = -dpoint(2)
  dvelocity(2) = dpoint(1)
  
  end subroutine

  
  
  
  
  
  
  
    !****************************************************************************
  
!<subroutine>
  
  subroutine assembleFaceTerms (rmatrix, rform, ccubType, velocityfield)
  
!<description>

  ! This Routine assembles the face terms of a dg discretisation

!</description>



 
    
    
    !<input>
  ! The matrix, to which to add the entries. Has to be prepared by
  !  subroutine bilf_createMatrixStructure (rdiscretisationTrial,iformat,rmatrixScalar, &
  !                                       rdiscretisationTest,cconstrType,imemguess)
  ! with cconstrType = BILF_MATC_EDGEBASED
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! The bilinear form specifying the underlying PDE of the discretisation.
  type(t_bilinearForm), intent(in) :: rform
  
  ! The type of cubature formula used to evaluate the face terms
  integer(I32), intent(in) :: ccubType
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  include 'intf_velocityfield.inc'
  optional :: velocityfield
  
  
!</input>

!<output>
  ! A matrix assembly structure.
  
!</output>

!</subroutine>
  
    ! local variables
    integer :: ielement, NEL, ccomplexity, IelementDistr, IelementDistrN, NNVE, NNEE, iedge, iglobalEdgeNumber, ineighbour, i, j, k, i1
    integer :: ilocaledgenumber, ilocaledgenumberneighbour
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrTrial, p_rspatialDiscrTest
    integer, dimension(:), pointer :: p_IelementDistr
    integer(I32) :: celementTest, celementTrial, celementTestN, celementTrialN
    
    ! Array to tell the element which derivatives to calculate
    logical, dimension(EL_MAXNDER) :: BderTrial, BderTrialTempl
    logical, dimension(EL_MAXNDER) :: BderTest, BderTestTempl
    
    ! Allocateable arrays for the values of the basis functions -
    ! for test and trial spaces.
    real(DP), dimension(:,:,:), allocatable, target :: DbasTest, DbasTrial
    real(DP), dimension(:,:,:), allocatable, target :: DbasTestN, DbasTrialN
    real(DP), dimension(:,:,:), pointer :: p_DbasTrial, p_DbasTest
    real(DP), dimension(:,:,:), pointer :: p_DbasTrialN, p_DbasTestN
    
    ! An allocateable array accepting the DOF`s of a set of elements.
    integer, dimension(:), allocatable, target :: IdofsTest, IdofsTrial
    integer, dimension(:), allocatable, target :: IdofsTestN, IdofsTrialN
    integer, dimension(:), pointer :: p_IdofsTrial, p_IdofsTest
    integer, dimension(:), pointer :: p_IdofsTrialN, p_IdofsTestN
    
    ! Pointer to the array saving the neighbous to one element
    integer, dimension(:,:), pointer :: p_IneighboursAtElement, p_IedgesAtElement
    
    ! Pointer to the cubature weights
    real(DP), dimension(:), pointer :: p_Domega
    
    ! An array that takes coordinates of the cubature formula on the reference element
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    
    integer :: ncubp, icubp
    
    ! Whether trial and test space is identical
    logical :: bIdenticalTrialAndTest, bIdenticalTrialAndTestN
    
    ! Type of transformation
    integer(I32) :: ctrafoType, ctrafoTypeN
    
    ! Basic evaluation tag of the element spaces
    integer(I32) :: cevaluationTag, cevaluationTagN
    
    ! Arrays for cubature points 1D->2D
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    real(DP), dimension(:,:), allocatable :: Dxi2D, DpointsRef
    
    ! Type of coordinatesystem
    integer(i32) :: icoordSystem, icoordSystemN
    
    ! Local matrices
    ! E : Element
    ! N : Neighbour
    !    integer, dimension(:,:,:), pointer :: p_Kentry
    integer, dimension(:,:), pointer :: p_KentryuEvE
    integer, dimension(:,:), pointer :: p_KentryuEvN
    integer, dimension(:,:), pointer :: p_KentryuNvE
    integer, dimension(:,:), pointer :: p_KentryuNvN
    !    real(DP), dimension(:,:,:), pointer :: p_Dentry
    real(DP), dimension(:,:), pointer :: p_DentryuEvE
    real(DP), dimension(:,:), pointer :: p_DentryuEvN
    real(DP), dimension(:,:), pointer :: p_DentryuNvE
    real(DP), dimension(:,:), pointer :: p_DentryuNvN
    
    integer :: IndofsTrial, IndofsTest, IndofsTrialN, IndofsTestN
    
    ! Something for elementevaluation
    type(t_evalElementSet), target :: revalElementSet
    type(t_evalElementSet), pointer :: p_revalElementSet
    
    ! Some temp variables to calculate the length of the edge
    real(dp) :: dxl1, dxl2, dyl1, dyl2, dlen
    
    ! Pointer to the array with the vertices at edge
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    ! Pointer to the coordinates of the vertices
    real(dp), dimension(:,:), pointer :: p_DvertexCoords
    
    ! Pointer to the vertices at one element
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    
    ! The normal vector of the edge
    real(dp), dimension(ndim2d) :: normal
    
    ! Whether trial and test space is identical
    logical :: bNeighbourExists
    
    ! Scalar product of velocity and normal vector
    real(DP) :: v_n
    
    ! The velocity vector
    real(DP), dimension(ndim2d) :: velocity
    
    ! The corner coordinates of the (center) element
    real(DP), dimension(:,:), pointer :: p_DcornerCoords
    
    ! Number of vertices of the element
    integer :: NVE
    
    
    
    ! DEPRECATED: Maximum number of basic functions = maximum number of
    ! local DOF`s per element.
    ! Do not use this constant anymore - determine the number of local basis
    ! functions dynamically using the 'elem_igetNDofLoc' routine!
    ! integer, parameter, public :: EL_MAXNBAS = 27
    ! DEPRECATED: maximal size of cubature node field
    ! integer, parameter, public :: CUB_MAXCUBP = 36
    
    
            
    ! Initialisation of the element set.
    ! Well, it's an empty routine anyway...
    call elprep_init(revalElementSet)
    p_revalElementSet => revalElementSet
    
    ! Get pointers to the discretisation of test and trial functions
    p_rspatialDiscrTrial => rmatrix%p_rspatialDiscrTrial
    p_rspatialDiscrTest  => rmatrix%p_rspatialDiscrTest
    bidenticalTrialAndTest = rmatrix%bidenticalTrialAndTest
    if (.not.bidenticalTrialAndTest) then
      call output_line ('Discretisations of test and trial functions do not fit', &
                        OU_CLASS_ERROR,OU_MODE_STD,'assembleFaceTerms')
      call sys_halt()
    end if
    
    ! Number of elements and edges and vertices per element in the disrcetisation
    NEL  = p_rspatialDiscrTrial%p_rtriangulation%NEL
    NNEE = p_rspatialDiscrTrial%p_rtriangulation%NNEE
    NNVE = p_rspatialDiscrTrial%p_rtriangulation%NNVE
    
    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(ccubType)
    
    ! Allocate array for the cubature weights
    allocate(p_Domega(ncubp))
    
    ! Allocate array for the points cubature points on the reference element
    allocate(p_DcubPtsRef(ndim3d, ncubp))
    
    ! Get the cubature formula
    call cub_getCubature(ccubType, p_DcubPtsRef, p_Domega)
    
    ! Transpose the coordinate array such that we get coordinates we
    ! can work with in the mapping between 1D and 2D.
    do k = 1, ubound(p_DcubPtsRef,1)
      do icubp = 1,ncubp
        Dxi1D(icubp,k) = p_DcubPtsRef(k,icubp)
      end do
    end do
    
    ! Get complexity of the discretisation
    ccomplexity = p_rspatialDiscrTrial%ccomplexity
    if (ccomplexity .ne. SPDISC_UNIFORM) then
      ! we've got more than one type of elementdistr
      ! Get pointer to the list, which allocates
      ! the distribution number to each element
      call storage_getbase_int( &
               p_rspatialDiscrTrial%h_IelementDistr,&
               p_IelementDistr)
    end if
    
    ! Get pointer to array, which saves the neighbourelements to one element
    call storage_getbase_int2d( &
    p_rspatialDiscrTrial%p_rtriangulation%h_IneighboursAtElement, &
    p_IneighboursAtElement)
    
    ! Get pointer to array, which saves the edges of one element
    call storage_getbase_int2d( &
    p_rspatialDiscrTrial%p_rtriangulation%h_IedgesAtElement, &
    p_IedgesAtElement)
    
    ! If we wanted to walk over all edges...
    !call storage_getbase_int2d( &
    !p_rspatialDiscrTrial%p_rtriangulation%h_IelementsAtEdge, &
    !p_IelementsAtEdge)
    
    ! Get pointer to the array with the vertices at edge
    CALL storage_getbase_int2D(p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    
    ! Get pointer to the coordinates of the vertices
    CALL storage_getbase_double2D(p_rspatialDiscrTest%p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    ! Get pointer to the vertices at one element
    CALL storage_getbase_int2D(p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
        
    ! As we always only integrate over one element simultaneously,
    ! we cann allocate maximum memory
    allocate(DbasTest (EL_MAXNBAS,EL_MAXNDER,ncubp))
    allocate(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp))
    allocate(DbasTestN (EL_MAXNBAS,EL_MAXNDER,ncubp))
    allocate(DbasTrialN(EL_MAXNBAS,EL_MAXNDER,ncubp))
    ! otherwise we would just use this in the loop
    !allocate(DbasTest(indofTest,&
    !         elem_getMaxDerivative(p_relementDistrTest%celement),&
    !         ncubp))
    !allocate(DbasTrial(indofTrial,&
    !         elem_getMaxDerivative(p_relementDistrTrial%celement), &
    !         ncubp)
    p_DbasTest  => DbasTest
    p_DbasTrial => DbasTrial
    p_DbasTestN  => DbasTestN
    p_DbasTrialN => DbasTrialN
    
    ! Allocate memory for the DOF`s of the element and its neighbour
    allocate(IdofsTrial(EL_MAXNBAS),IdofsTest(EL_MAXNBAS))
    allocate(IdofsTrialN(EL_MAXNBAS),IdofsTestN(EL_MAXNBAS))
    p_IdofsTrial => IdofsTrial
    p_IdofsTest  => IdofsTest
    p_IdofsTrialN => IdofsTrialN
    p_IdofsTestN  => IdofsTestN

    ! Allocate memory for the cubature points in 2D.
    allocate(Dxi2D(ncubp,NDIM3D))

    ! Allocate memory for the coordinates of the reference points
    allocate(DpointsRef(NDIM3D,ncubp))
    
    ! Allocate memory for the local matrices
    !allocate(p_Kentry(indofTrial,indofTest))
    !allocate(p_Dentry(indofTrial,indofTest))
    allocate(p_KentryuEvE(EL_MAXNBAS,EL_MAXNBAS))
    allocate(p_DentryuEvE(EL_MAXNBAS,EL_MAXNBAS))
    allocate(p_KentryuEvN(EL_MAXNBAS,EL_MAXNBAS))
    allocate(p_DentryuEvN(EL_MAXNBAS,EL_MAXNBAS))
    allocate(p_KentryuNvE(EL_MAXNBAS,EL_MAXNBAS))
    allocate(p_DentryuNvE(EL_MAXNBAS,EL_MAXNBAS))
    allocate(p_KentryuNvN(EL_MAXNBAS,EL_MAXNBAS))
    allocate(p_DentryuNvN(EL_MAXNBAS,EL_MAXNBAS))
        
    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDERxxxx
    ! according to these.
    BderTrialTempl = .false.
    BderTestTempl = .false.
    
    ! Loop through the additive terms
    do i=1,rform%itermCount
      ! The desriptor Idescriptors gives directly the derivative
      ! which is to be computed! Build templates for BDER.
      ! We do not compute the actual BDER here, as there might be some special
      ! processing if trial/test functions are identical!
      !
      ! At first build the descriptors for the trial functions
      I1=rform%Idescriptors(1,I)
      
      if ((I1 .le.0) .or. (I1 .gt. DER_MAXNDER)) then
        call output_line ('Invalid descriptor!',&
            OU_CLASS_ERROR,OU_MODE_STD,'assembleFaceTerms')
        call sys_halt()
      endif
      
      BderTrialTempl(I1)=.true.

      ! Then those of the test functions
      I1=rform%Idescriptors(2,I)
      
      if ((I1 .le.0) .or. (I1 .gt. DER_MAXNDER)) then
        call output_line ('Invalid descriptor!',&
            OU_CLASS_ERROR,OU_MODE_STD,'assembleFaceTerms')
        call sys_halt()
      endif
      
      BderTestTempl(I1)=.true.
    end do

    if (bIdenticalTrialAndTest) then
      ! Build the actual combination of what the element should calculate.
      BderTrial(:) = BderTrialTempl(:) .or. BderTestTempl(:)
      BderTest(:) = BderTrial(:)
    else
      ! Build the actual combination of what the element should calculate.
      ! Copy BDERxxxx to BDERxxxxAct
      BderTrial(:) = BderTrialTempl(:)
      BderTest(:) = BderTestTempl(:)
    end if

    
     ! first we will loop over all the elements in the discretisation
    do ielement = 1, NEL
    
    ! Determine the type of element we are working at the moment
    ielementDistr = 1
    ! In which elementdistribution is the element?
    if (ccomplexity .ne. SPDISC_UNIFORM) then
      ! we've got only one type of elementdistr
      ielementDistr = p_IelementDistr(ielement)
    end if
    
    ! Get type of test and trial function on this element
    celementTest  = p_rspatialDiscrTest%RelementDistr(ielementDistr)%celement
    celementTrial = p_rspatialDiscrTrial%RelementDistr(ielementDistr)%celement
    
    ! Get the number of degrees of freedom of the trial and test functions
    IndofsTrial = elem_igetNDofLoc(celementTrial)
    IndofsTest  = elem_igetNDofLoc(celementTest)
    
    
    
    ! Get the number of vertices of the element, specifying the transformation
    ! form the reference to the real element.
    NVE = elem_igetNVE(celementTest)
    
    ! Allocate memory for the corner coordinates
    allocate (p_DcornerCoords(ndim2d,NVE))
    
    ! Get the corner coordinates of the element
    call trafo_getCoords (ctrafoType,p_rspatialDiscrTrial%p_rtriangulation,ielement,Dcoords)
    
    
       
      ! loop over all edges of the element
      do iedge = 1, NNEE
      
        ! Get the number of the neighbour element
        ineighbour = p_IneighboursAtElement(iedge, ielement)
        
        ! Maybe the neighbour doesn't exist (as we are at a boundary) then we don't do anything
        ! (well, we could do something - maybe later)
        bNeighbourExists = .true.
        if (ineighbour.eq.0) then
          bNeighbourExists = .false.
          cycle
        end if
        
        
        ! Find global and local edgenumber
         do i = 1, NNEE
         do j = 1, NNEE
           if ((p_IedgesAtElement(i,ielement).eq. &
                p_IedgesAtElement(j,ineighbour)).and. &
               (p_IedgesAtElement(i,ielement).ne.0)) then
                 iglobalEdgeNumber = p_IedgesAtElement(i,ielement)
                 ilocaledgenumber = i
                 ilocaledgenumberneighbour = j
               end if
         end do
         end do
         
         ! Calculate the length of the edge *0.5
         dxl1=p_DvertexCoords(1,p_IverticesAtEdge(1,iglobalEdgeNumber))
         dyl1=p_DvertexCoords(2,p_IverticesAtEdge(1,iglobalEdgeNumber))
         dxl2=p_DvertexCoords(1,p_IverticesAtEdge(2,iglobalEdgeNumber))
         dyl2=p_DvertexCoords(2,p_IverticesAtEdge(2,iglobalEdgeNumber))
         dlen=0.5_DP*sqrt((dxl1-dxl2)*(dxl1-dxl2)+(dyl1-dyl2)*(dyl1-dyl2))
         
         ! Calculate the normal vector to the element at this edge
         normal(1) = (dyl2-dyl1)/dlen*0.5_DP
         normal(2) = (dxl1-dxl2)/dlen*0.5_DP
         
         do i = 1, NNVE
           if (p_IverticesAtEdge(1,iglobalEdgeNumber).eq.p_IverticesAtElement(i,ielement)) then
           j = mod(i + 1,NNVE)
             if (p_IverticesAtEdge(2,iglobalEdgeNumber).ne.p_IverticesAtElement(j,ielement)) then
               normal(1) = -normal(1)
               normal(2) = -normal(2)
             end if
           exit
           end if
         end do
         
         ! Get from the element space the type of coordinate system
         ! that is used there:
         ctrafoType = elem_igetTrafoType(celementTest)
      
         ! Determine if trial and test space is the same.
         bIdenticalTrialAndTest = (celementTest .eq. celementTrial)

        ! Get the element evaluation tag of all FE spaces. We need it to evaluate
        ! the elements later. All of them can be combined with OR, what will give
        ! a combined evaluation tag.
        cevaluationTag = elem_getEvaluationTag(celementTest)
        cevaluationTag = ior(cevaluationTag,&
                         elem_getEvaluationTag(celementTrial))
        ! The cubature points are already initialised by 1D->2D mapping.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
         
        ! Get the type of coordinate system
        icoordSystem = elem_igetCoordSystem(celementTrial)
         
        ! Map the 1D cubature points to the edges in 2D.
        call trafo_mapCubPts1Dto2D(icoordSystem, ilocaledgenumber, &
                                   ncubp, Dxi1D, Dxi2D)
              
        ! Transpose the coordinate array such that we get coordinates we
        ! can work with.
        do icubp = 1,ncubp
          do k = 1,ubound(DpointsRef,1)
            DpointsRef(k,icubp) = Dxi2D(icubp,k)
          end do
        end do
        
        ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
        ! We call dof_locGlobMapping to calculate the
        ! global DOF`s of our element
        call dof_locGlobMapping(p_rspatialDiscrTest, &
                                Ielement, p_IdofsTest)
                                   
        p_IdofsTrial = p_IdofsTest
        ! If the DOF`s for the test functions are different, calculate them, too.
        if (.not. bIdenticalTrialAndTest) then
          call dof_locGlobMapping(p_rspatialDiscrTrial, &
                                  Ielement, p_IdofsTrial)
        end if
        
        ! We build local matrices
        ! Get the positions of the local matrices in the global matrix.
        call getLocalMatrixIndices_single (rmatrix,p_IdofsTrial,p_IdofsTest,p_KentryuEvE,&
                                           IndofsTrial,IndofsTest)
        
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call prepareOneForMultipleEvaluation (p_revalElementSet,&
            cevaluationTag, p_rspatialDiscrTest%p_rtriangulation, &
            Ielement, ctrafoType, DpointsRef)
                
        ! Calculate the values of the basis functions.
        call generic_sim2_single (celementTest, &
                                  p_revalElementSet, BderTest, &
                                  p_DbasTest)
      
        ! Omit the calculation of the trial function values if they
        ! are identical to the test function values.
        if (.not.bidenticalTrialAndTest) then
          p_DbasTrial => DbasTrial
          call generic_sim2_single (celementTrial, &
              p_revalElementSet, BderTrial, &
              p_DbasTrial)
        else
          p_DbasTrial => DbasTest
        end if
        
        ! Allocate
        
        ! Map the cubature points from the reference to the real element
        ! to evaluate the velocity field
        CALL trafo_calctrafo (ctrafoType,Dcoords,&
                              DpointRef,Djac,ddetj,DpointReal)
        
        
        
        ! Now do everything we have done for the Element for the Neighbourelement
        if (bNeighbourExists) then
        
          ! Determine the type of element we are working at the moment
          ielementDistrN = 1
          ! In which elementdistribution is the element?
          if (ccomplexity .ne. SPDISC_UNIFORM) then
            ! we've got only one type of elementdistr
            ielementDistrN = p_IelementDistr(ineighbour)
          end if
    
          ! Get type of test and trial function on this element
          celementTestN  = p_rspatialDiscrTest%RelementDistr(ielementDistrN)%celement
          celementTrialN = p_rspatialDiscrTrial%RelementDistr(ielementDistrN)%celement
    
          ! Get the number of degrees of freedom of the trial and test functions
          IndofsTrialN = elem_igetNDofLoc(celementTrialN)
          IndofsTestN  = elem_igetNDofLoc(celementTestN)
          
          ! Get from the element space the type of coordinate system
          ! that is used there:
          ctrafoTypeN = elem_igetTrafoType(celementTestN)
      
          ! Determine if trial and test space is the same.
          bIdenticalTrialAndTestN = (celementTestN .eq. celementTrialN)
          
          ! Get the element evaluation tag of all FE spaces. We need it to evaluate
          ! the elements later. All of them can be combined with OR, what will give
          ! a combined evaluation tag.
          cevaluationTagN = elem_getEvaluationTag(celementTestN)
          cevaluationTagN = ior(cevaluationTagN,&
                           elem_getEvaluationTag(celementTrialN))
          ! The cubature points are already initialised by 1D->2D mapping.
          cevaluationTagN = iand(cevaluationTagN,not(EL_EVLTAG_REFPOINTS))
     
          ! Get the type of coordinate system
          icoordSystemN = elem_igetCoordSystem(celementTrialN)
         
          ! Map the 1D cubature points to the edges in 2D.
          call trafo_mapCubPts1Dto2D(icoordSystemN, ilocaledgenumberneighbour, &
                                     ncubp, Dxi1D, Dxi2D)

          ! Transpose the coordinate array such that we get coordinates we
          ! can work with.
          do icubp = 1,ncubp
            do k = 1,ubound(DpointsRef,1)
              DpointsRef(k,icubp) = Dxi2D(icubp,k)
            end do
          end do
        
          ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
          ! We call dof_locGlobMapping to calculate the
          ! global DOF`s of our element
          call dof_locGlobMapping(p_rspatialDiscrTest, &
                                  Ineighbour, p_IdofsTestN)
                                                                     
          p_IdofsTrialN = p_IdofsTestN
          ! If the DOF`s for the test functions are different, calculate them, too.
          if (.not. bIdenticalTrialAndTestN) then
            call dof_locGlobMapping(p_rspatialDiscrTrial, &
                                    Ineighbour, p_IdofsTrialN)
          end if

          ! We build local matrices
          ! Get the positions of the local matrices in the global matrix.
          call getLocalMatrixIndices_single (rmatrix,p_IdofsTrialN,p_IdofsTestN,p_KentryuNvN,&
                                             IndofsTrialN,IndofsTestN)
          call getLocalMatrixIndices_single (rmatrix,p_IdofsTrial,p_IdofsTestN,p_KentryuEvN,&
                                             IndofsTrial,IndofsTestN)
          call getLocalMatrixIndices_single (rmatrix,p_IdofsTrialN,p_IdofsTest,p_KentryuNvE,&
                                             IndofsTrialN,IndofsTest)
                                             
        
          ! Calculate all information that is necessary to evaluate the finite element
          ! on all cells of our subset. This includes the coordinates of the points
          ! on the cells.
          call prepareOneForMultipleEvaluation (p_revalElementSet,&
              cevaluationTagN, p_rspatialDiscrTest%p_rtriangulation, &
              Ineighbour, ctrafoTypeN, DpointsRef)
        
        
          ! Calculate the values of the basis functions.
          call generic_sim2_single (celementTestN, &
                                    p_revalElementSet, BderTest, &
                                    p_DbasTestN)

          ! Omit the calculation of the trial function values if they
          ! are identical to the test function values.
          if (.not.bidenticalTrialAndTest) then
            p_DbasTrial => DbasTrial
            call generic_sim2_single (celementTrialN, &
                p_revalElementSet, BderTrial, &
                p_DbasTrialN)
          else
            p_DbasTrialN => DbasTestN
          end if
        
        end if ! Neighbour exists
        
        
        ! Loop over all cubature points on the current edge
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          domega = dlen * p_Domega(icubp)
        
          ! Calculate the velocity in that cubature point
          
        
          ! Calculate the scalar product of velocity and normal vector
          v_n = normal(1)*velocity(1) + normal(2)*velocity(2)
        
        end do ! icubp
                
        
        
        
        
        
         
      end do ! loop over all edges of the element
    
    deallocate (p_DcornerCoords)
    
    end do ! loop over all elements
    
    
    
    
    
    
    
    
    
    
    
    
    
    ! deaallocate all memory
    deallocate(DbasTest,DbasTrial)
    deallocate(DbasTestN,DbasTrialN)
    deallocate(IdofsTest,IdofsTrial)
    deallocate(IdofsTestN,IdofsTrialN)
    deallocate(p_Domega)
    deallocate(p_DcubPtsRef)
    deallocate(Dxi2D)
    deallocate(DpointsRef)
    deallocate(p_KentryuEvE)
    deallocate(p_DentryuEvE)
    deallocate(p_KentryuEvN)
    deallocate(p_DentryuEvN)
    deallocate(p_KentryuNvE)
    deallocate(p_DentryuNvE)
    deallocate(p_KentryuNvN)
    deallocate(p_DentryuNvN)
    
  end subroutine

    
end module
