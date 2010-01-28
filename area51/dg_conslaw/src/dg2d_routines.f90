module dg2d_routines

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
    use collection
    use linearformevaluation
    use domainintegration
    use linearsystemscalar
    use feevaluation
    

    implicit none
    
    type t_array
		! Pointer to the double-valued matrix or vector data
    	real(DP), dimension(:), pointer :: Da
	end type t_array
	
	
	public :: linf_dg_buildVectorScalarEdge2d

contains

        
    

!****************************************************************************

!<subroutine>

  subroutine linf_dg_buildVectorScalarEdge2d (rform, ccubType, bclear,&
                                              rvectorScalar,&
                                              rvectorScalarSol,&
                                              flux_dg_buildVectorScEdge2D_sim,&
                                              rcollection)
  
!<description>
  ! This routine assembles the entries of a vector according to the linear form
  ! \int_{elementboundary} dg_flux_function * \phi ds
  ! The flux function is given as a callback routine and is dependend on the
  ! two values of the solution on the edge, according to the left and right
  ! element
  !
  ! If bclear=TRUE, the vector is cleared before the assembly and any 
  ! sorting of the entries is switched off - the vector is set up unsorted.
  !
  ! If bclear=FALSE, the vector must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
!</description>

!<input>
  ! The linear form specifying the underlying PDE of the discretisation.
  type(t_linearForm), intent(in) :: rform

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer(I32), intent(in) :: ccubType
  
  ! Whether to clear the vector before calculating the entries.
  ! If .FALSE., the new entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! The solution vector. Used to calculate the solution on the edges.
  type(t_vectorScalar), intent(in) :: rvectorScalarSol
   
  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection
  
  ! A callback routine for the function to be discretised.
  include 'intf_flux_dg_buildVectorScEdge2D.inc'
  optional :: flux_dg_buildVectorScEdge2D_sim
!</input>

!<inputoutput>
  ! The linear form vector. Calculated entries are imposed to this vector.
  type(t_vectorScalar), intent(inout) :: rvectorScalar
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_linfVectorAssembly), dimension(2) :: rvectorAssembly
  type(t_boundary), pointer :: p_rboundary
  type(t_triangulation), pointer :: p_rtriangulation
  integer, dimension(:), pointer :: p_IedgeList
  integer, dimension(:,:), pointer :: p_IlocalEdgeNumber
  integer :: ielementDistr, NMT, NVT, iedge

  ! If the vector does not exist, stop here.
  if (rvectorScalar%h_Ddata .eq. ST_NOHANDLE) then  
    call output_line('Vector not available!',&
                     OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
  end if

  ! The vector must be unsorted, otherwise we can not set up the vector.
  if (rvectorScalar%isortStrategy .gt. 0) then
    call output_line('Vector must be unsorted!',&
                     OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalarBdr2D')
    call sys_halt()
  end if
  
  ! Clear the vector if necessary.
  if (bclear) call lsyssc_clearVector (rvectorScalar)

  ! The vector must provide a discretisation structure
  if (.not. associated(rvectorScalar%p_rspatialDiscr)) then
    call output_line('No discretisation associated!',&
        OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
    call sys_halt()
  end if

  ! The discretisation must provide a triangulation structure
  if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rtriangulation)) then
    call output_line('No triangulation associated!',&
        OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
    call sys_halt()
  end if
  
  ! The discretisation must provide a boundary structure
  if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rboundary)) then
    call output_line('No boundary associated!',&
        OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
    call sys_halt()
  end if

  ! Set pointers for quicker access
  p_rboundary => rvectorScalar%p_rspatialDiscr%p_rboundary
  p_rtriangulation => rvectorScalar%p_rspatialDiscr%p_rtriangulation
    
  ! Do we have a uniform triangulation? Would simplify a lot...
  if (rvectorScalar%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then 
    
    select case(rvectorScalar%cdataType)
      
    case(ST_DOUBLE)

      ! Get number of vertices in the triangulation
      NVT = p_rtriangulation%NVT

      ! Get number of edges belonging to elements
      NMT = p_rtriangulation%NMT
          
      ! Allocate the edgelist
      allocate(p_IedgeList(NMT))
      
      ! All edges
      forall (iedge = 1:NMT) p_IedgeList(iedge)=iedge
      
      ! Allocate space for local edge numbers
      allocate(p_IlocalEdgeNumber(2,NMT))
      
      ! Get local edge numbers
      call linf_getLocalEdgeNumbers(p_IedgeList,p_IlocalEdgeNumber,&
                                    p_rtriangulation)
                                    
     
      ! Initialise the vectorAssembly structures
      call linf_initAssembly(rvectorAssembly(1), rform,&
            rvectorScalar%p_rspatialDiscr%RelementDistr(1)%celement,&
            ccubType, LINF_NELEMSIM)
      call linf_initAssembly(rvectorAssembly(2), rform,&
            rvectorScalar%p_rspatialDiscr%RelementDistr(1)%celement,&
            ccubType, LINF_NELEMSIM)
            
      ! Assemble the data for all elements in this element distribution
      call linf_dg_assembleSubmeshVectorScalarEdge2d (rvectorAssembly,&
            rvectorScalar, rvectorScalarSol,&
            p_IedgeList(1:NMT),p_IlocalEdgeNumber&
            ,flux_dg_buildVectorScEdge2D_sim,&
            rcollection&
             )

          
      ! Release the assembly structure.
      call linf_doneAssembly(rvectorAssembly(1))
      call linf_doneAssembly(rvectorAssembly(2))

      ! Deallocate the edgelist
      deallocate(p_IedgeList)
      
      ! Deallocate space for local edge numbers
      deallocate(p_IlocalEdgeNumber)
      
    case DEFAULT
      call output_line('Single precision vectors currently not supported!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
      call sys_halt()
    end select
    
  else
    call output_line('General discretisation not implemented!',&
                     OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
    call sys_halt()
  end if
 
  end subroutine linf_dg_buildVectorScalarEdge2d















































  !****************************************************************************
  
!<subroutine>  
  
  subroutine linf_dg_assembleSubmeshVectorScalarEdge2d (rvectorAssembly,&
              rvector, rvectorSol,&
              IedgeList, IlocalEdgeNumber,&
              flux_dg_buildVectorScEdge2D_sim,&
              rcollection)

!<description>

  ! Assembles the vector entries for a submesh by integration over the given edges.

!</description>

!<input>

  ! List of edges where to assemble the linear form.
  integer, dimension(:), intent(in), target :: IedgeList
  
  ! List of local edge numbers for the two elements adjacent to the two edges.
  integer, dimension(:,:), intent(in) :: IlocalEdgeNumber
  
  ! The solution vector. Used to calculate the solution on the edges.
  type(t_vectorScalar), intent(in) :: rvectorSol

  ! A callback routine which is able to calculate the values of the
  ! function $f$ which is to be discretised.
  include 'intf_flux_dg_buildVectorScEdge2D.inc'
  optional :: flux_dg_buildVectorScEdge2D_sim 
  
!</input>

!<inputoutput>
  
  ! A vector assembly structure prepared with linf_initAssembly.
  type(t_linfVectorAssembly), dimension(2), intent(inout), target :: rvectorAssembly
  
  ! A vector where to assemble the contributions to.
  type(t_vectorScalar), intent(inout) :: rvector  
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection
!</inputoutput>
  
!</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata, p_DdataSol
    integer :: indof,ncubp
    
    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax,ibdc,k
    integer :: iel,icubp,ialbet,ia,idofe
    real(DP) :: domega,daux1,daux2,dlen
    real(DP) :: dval1, dval2
    integer(I32) :: cevaluationTag
    type(t_linfVectorAssembly), dimension(2), target :: rlocalVectorAssembly
    type(t_domainIntSubset) :: rintSubset
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    integer, dimension(:),pointer :: p_Idescriptors
    integer, dimension(:,:), pointer :: p_Idofs
    type(t_evalElementSet), pointer :: p_revalElementSet

    ! A small vector holding only the additive contributions of
    ! one element
    real(DP), dimension(2,EL_MAXNBAS) :: DlocalData
  
    ! Arrays for cubature points 1D->2D
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    real(DP), dimension(:,:,:,:), allocatable :: Dxi2D,DpointsRef
    
    ! Element list, where to assemble the form
    integer, dimension(:,:), allocatable :: IelementList
    
    integer(i32) :: icoordSystem
    
    ! Chooses the element on side ... of the edge
    integer :: iside
    
    ! Number of edges
    integer :: NMT
    
    ! Pointer to Ielementsatedge in the triangulation
    integer, dimension(:,:), pointer :: p_IelementsAtEdge
    
    ! Pointer to IverticesAtEdge in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    ! Pointer to the vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    
    ! Array for the solution values in the cubature points
    real(DP), dimension(:,:,:), allocatable :: DsolVals
    
    ! Array for the length of the edges
    real(DP), dimension(:), allocatable :: edgelength
         
    ! Array for the normal vectors
    real(DP), dimension(:,:), allocatable :: normal
    
    ! Temp variables for the coordinates of the vertices
    real(DP) :: dxl1, dxl2, dyl1, dyl2
  
  
    ! Copy the assembly data to the local assembly data,
    ! where we can allocate memory.
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the vector
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original assembly
    ! stucture or disturbing the data of the other processors.
    rlocalVectorAssembly(1) = rvectorAssembly(1)
    rlocalVectorAssembly(2) = rvectorAssembly(2)
    
    call linf_allocAssemblyData(rlocalVectorAssembly(1))
    call linf_allocAssemblyData(rlocalVectorAssembly(2))
  
    ! Get pointers to elements at edge
    call storage_getbase_int2D(&
          rvector%p_rspatialDiscr%p_rtriangulation%h_IelementsAtEdge,&
          p_IelementsAtEdge)
          
    ! Get pointers to the vertex coordinates
    call storage_getbase_double2D(&
          rvector%p_rspatialDiscr%p_rtriangulation%h_DvertexCoords,&
          p_DvertexCoords)
          
    ! Get pointers to vertices at edge
    call storage_getbase_int2D(&
          rvector%p_rspatialDiscr%p_rtriangulation%h_IverticesAtEdge,&
          p_IverticesAtEdge)          
    
    ! Get the elements adjacent to the given edges
    allocate(IelementList(3,size(IedgeList)))
    IelementList(1:2,1:size(IedgeList))=p_IelementsAtEdge(1:2,IedgeList(:))
    
    ! Copy the second component and replace 0s by 1s
    IelementList(3,size(IedgeList))=IelementList(2,size(IedgeList))
    do iel = 1,size(IedgeList)
      IelementList(3,iel)=max(IelementList(2,iel),1)
    end do
    
    ! Get some pointers for faster access
    call lsyssc_getbase_double (rvector, p_Ddata)
    indof = rvectorAssembly(1)%indof
    ncubp = rvectorAssembly(1)%ncubp
    
    ! Allocate space for the solution values in the cubature points
    allocate(DsolVals(ncubp,2,rlocalVectorAssembly(1)%nelementsPerBlock))

    ! Allocate space for normal vectors
    allocate(normal(2,min(size(IedgeList),rlocalVectorAssembly(1)%nelementsPerBlock)))
    
    ! Allocate space for edge length
    allocate(edgelength(min(size(IedgeList),rlocalVectorAssembly(1)%nelementsPerBlock)))
   
!    ! Get some more pointers to local data.
!    p_Domega => rlocalVectorAssembly%p_Domega
!    p_Dbas => rlocalVectorAssembly%p_Dbas
!    p_Dcoefficients => rlocalVectorAssembly%p_Dcoefficients
!    p_DcubPtsRef => rlocalVectorAssembly%p_DcubPtsRef
!    p_Idescriptors => rlocalVectorAssembly%rform%Idescriptors
!    p_Idofs => rlocalVectorAssembly%p_Idofs
!    p_revalElementSet => rlocalVectorAssembly%revalElementSet
    
    ! Transpose the coordinate array such that we get coordinates we
    ! can work with in the mapping between 1D and 2D.
    do k = 1, ubound(rlocalVectorAssembly(1)%p_DcubPtsRef,1)
      do icubp = 1,ncubp
        Dxi1D(icubp,k) = rlocalVectorAssembly(1)%p_DcubPtsRef(k,icubp)
      end do
    end do
 
    ! Allocate memory for the cubature points in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,2,rlocalVectorAssembly(1)%nelementsPerBlock))

    ! Allocate memory for the coordinates of the reference points
    allocate(DpointsRef(NDIM2D+1,ncubp,2,rlocalVectorAssembly(1)%nelementsPerBlock))

    ! Get the type of coordinate system
    icoordSystem = elem_igetCoordSystem(rlocalVectorAssembly(1)%celement)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !%OMP do schedule(static,1)
    do IELset = 1, size(IedgeList), rlocalVectorAssembly(1)%nelementsPerBlock

      ! We always handle nelementsPerBlock elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(IedgeList),IELset-1+rlocalVectorAssembly(1)%nelementsPerBlock)
      
      ! Map the 1D cubature points to the edges in 2D.
      do iel = 1,IELmax-IELset+1
        call trafo_mapCubPts1Dto2D(icoordSystem, IlocalEdgeNumber(1,IELset+iel-1), &
            ncubp, Dxi1D, Dxi2D(:,:,1,iel))
        call trafo_mapCubPts1Dto2D(icoordSystem, IlocalEdgeNumber(2,IELset+iel-1), &
            ncubp, Dxi1D, Dxi2D(:,:,2,iel))    
      end do
     
      ! Transpose the coordinate array such that we get coordinates we
      ! can work with.
      do iel = 1,IELmax-IELset+1
        do iside = 1,2
          do icubp = 1,ncubp
            do k = 1,ubound(DpointsRef,1)
              DpointsRef(k,icubp,iside,iel) = Dxi2D(icubp,k,iside,iel)
            end do
          end do
        end do
      end do
        
      ! --------------------- DOF SEARCH PHASE ------------------------
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF`s on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | iel |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element iel, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(rvector%p_rspatialDiscr, &
          IelementList(1,IELset:IELmax), rlocalVectorAssembly(1)%p_Idofs)
      call dof_locGlobMapping_mult(rvector%p_rspatialDiscr, &
          IelementList(3,IELset:IELmax), rlocalVectorAssembly(2)%p_Idofs)
                                   
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! To calculate the element contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = rlocalVectorAssembly(1)%cevaluationTag
      
      ! The cubature points are already initialised by 1D->2D mapping.
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

      ! Calculate all information that is necessary to evaluate the
      ! finite element on all cells of our subset. This includes the
      ! coordinates of the points on the cells.
      call elprep_prepareSetForEvaluation (&
          rlocalVectorAssembly(1)%revalElementSet,&
          cevaluationTag, rvector%p_rspatialDiscr%p_rtriangulation, &
          IelementList(1,IELset:IELmax), rlocalVectorAssembly(1)%ctrafoType, &
          DpointsRef=DpointsRef(:,:,1,:))
      call elprep_prepareSetForEvaluation (&
          rlocalVectorAssembly(2)%revalElementSet,&
          cevaluationTag, rvector%p_rspatialDiscr%p_rtriangulation, &
          IelementList(3,IELset:IELmax), rlocalVectorAssembly(2)%ctrafoType, &
          DpointsRef=DpointsRef(:,:,2,:))
          
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (rlocalVectorAssembly(1)%celement, &
          rlocalVectorAssembly(1)%revalElementSet,&
          rlocalVectorAssembly(1)%Bder, &
          rlocalVectorAssembly(1)%p_Dbas)
      call elem_generic_sim2 (rlocalVectorAssembly(2)%celement, &
          rlocalVectorAssembly(2)%revalElementSet,&
          rlocalVectorAssembly(2)%Bder, &
          rlocalVectorAssembly(2)%p_Dbas)
          
          
         
      ! ********** Get solution values in the cubature points *************
      call lsyssc_getbase_double(rvectorSol,p_DdataSol)
    
      ! Now that we have the basis functions, we want to have the function values.
      ! We get them by multiplying the FE-coefficients with the values of the
      ! basis functions and summing up.
      do iel = 1,IELmax-IELset+1      
        do icubp = 1,ncubp
          ! Calculate the value in the point
          dval1 = 0.0_DP
          dval2 = 0.0_DP
          do idofe = 1,indof
            dval1 = dval1 + &
                   p_DdataSol(rlocalVectorAssembly(1)%p_Idofs(idofe,iel)) &
                   * rlocalVectorAssembly(1)%p_Dbas(idofe,DER_FUNC,icubp,iel)
            dval2 = dval2 + &
                   p_DdataSol(rlocalVectorAssembly(2)%p_Idofs(idofe,iel)) &
                   * rlocalVectorAssembly(2)%p_Dbas(idofe,DER_FUNC,icubp,iel)
          end do
          ! Save the value in the point
          DsolVals(icubp,1,iel) = dval1
          DsolVals(icubp,2,iel) = dval2
        end do
      end do
      
     
     do iel = 1,IELmax-IELset+1
      if(IelementList(2,IELset+iel-1).eq.0) then
        DsolVals(1:ncubp,2,iel) = 1.0_DP
      end if
      
    end do
     
      
     ! --------------------- Get normal vectors ---------------------------
     
     do iel = 1,IELmax-IELset+1    
       ! Calculate the length of the edge 
       dxl1=p_DvertexCoords(1,p_IverticesAtEdge(1,IedgeList(IELset+iel-1)))
       dyl1=p_DvertexCoords(2,p_IverticesAtEdge(1,IedgeList(IELset+iel-1)))
       dxl2=p_DvertexCoords(1,p_IverticesAtEdge(2,IedgeList(IELset+iel-1)))
       dyl2=p_DvertexCoords(2,p_IverticesAtEdge(2,IedgeList(IELset+iel-1)))
       edgelength(iel)=sqrt((dxl1-dxl2)*(dxl1-dxl2)+(dyl1-dyl2)*(dyl1-dyl2))
         
       ! Calculate the normal vector to the element at this edge
       normal(1,iel) = (dyl2-dyl1)/edgelength(iel)
       normal(2,iel) = (dxl1-dxl2)/edgelength(iel)
     end do
      
      
     ! ---------------------- Get values of the flux function --------------
     
     
     !!      ! Now it is time to call our coefficient function to calculate the
!!      ! function values in the cubature points:
!!      if (present(fcoeff_buildVectorScBdr2D_sim)) then
!!        call domint_initIntegrationByEvalSet (p_revalElementSet, rintSubset)
!!        rintSubset%ielementDistribution = 0
!!        rintSubset%ielementStartIdx = IELset
!!        rintSubset%p_Ielements => IelementList(IELset:IELmax)
!!        rintSubset%p_IdofsTrial => p_Idofs
!!        rintSubset%celement = rlocalVectorAssembly%celement
!!        call fcoeff_buildVectorScBdr2D_sim (rvector%p_rspatialDiscr,&
!!            rlocalVectorAssembly%rform,  IELmax-IELset+1, ncubp,&
!!            p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!!            ibdc, DpointsPar(:,1:IELmax-IELset+1),&
!!            p_Idofs, rintSubset, &
!!            p_Dcoefficients(:,:,1:IELmax-IELset+1), rcollection)
!!        call domint_doneIntegration (rintSubset)
!!      else
!!        p_Dcoefficients(:,:,1:IELmax-IELset+1) = 1.0_DP
!!      end if
      
      
      call flux_dg_buildVectorScEdge2D_sim (&
            rlocalVectorAssembly(1)%p_Dcoefficients(1,:,1:IELmax-IELset+1),&
            DsolVals(:,:,1:IELmax-IELset+1),&
            normal(:,1:IELmax-IELset+1),&
            rcollection )
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the DOF`s and cubature points to calculate the
      ! integral:

      do iel = 1,IELmax-IELset+1
        
        ! We make a 'local' approach, i.e. we calculate the values of the
        ! integral into the vector DlocalData and add them later into
        ! the large solution vector.
        
        ! Clear the output vector.
        DlocalData(1:2,1:indof) = 0.0_DP

        ! Get the length of the edge. Let us use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit interval [-1,1] has length 2.
        dlen = 0.5_DP*edgelength(iel)

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature
          ! formula in that cubature point.

          domega = dlen * rlocalVectorAssembly(1)%p_Domega(icubp)

          
          ! Loop over the additive factors in the bilinear form.
          do ialbet = 1,rlocalVectorAssembly(1)%rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be:
            !
            ! int_...  f * ( phi_i )_IA
            !
            ! -> IA=0: function value, 
            !      =1: first derivative, 
            !      =2: 2nd derivative,...
            !    as defined in the module 'derivative'.
            
            ia = rlocalVectorAssembly(1)%rform%Idescriptors(ialbet)
            
            ! Multiply domega with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            daux1 = domega * rlocalVectorAssembly(1)%p_Dcoefficients(ialbet,icubp,iel)
            daux2 = domega * rlocalVectorAssembly(1)%p_Dcoefficients(ialbet,ncubp-icubp+1,iel) *(-1.0_dp)

            ! Now loop through all possible combinations of DOF`s
            ! in the current cubature point. 

            do idofe = 1,indof
              
              ! Get the value of the basis function 
              ! phi_o in the cubature point. 
              ! Them multiply:
              !    DBAS(..) * AUX
              ! ~= phi_i * coefficient * cub.weight
              ! Summing this up gives the integral, so the contribution
              ! to the vector. 
              !
              ! Simply summing up DBAS(..) * AUX would give
              ! the additive contribution for the vector. We save this
              ! contribution in the local array.
              
              DlocalData(1,idofe) = DlocalData(1,idofe)+&
                                    rlocalVectorAssembly(1)%p_Dbas(idofe,ia,icubp,iel)*daux1
              
              if(IelementList(2,IELset+iel-1).ne.0) then
                DlocalData(2,idofe) = DlocalData(2,idofe)+&
                                      rlocalVectorAssembly(2)%p_Dbas(idofe,ia,ncubp-icubp+1,iel)*daux2
              end if
             
            end do ! idofe
            
          end do ! ialbet

        end do ! icubp 
        
        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        do idofe = 1,indof
        
        
          p_Ddata(rlocalVectorAssembly(1)%p_Idofs(idofe,iel)) =&
                       p_Ddata(rlocalVectorAssembly(1)%p_Idofs(idofe,iel)) +&
                       DlocalData(1,idofe)
          p_Ddata(rlocalVectorAssembly(2)%p_Idofs(idofe,iel)) =&
                       p_Ddata(rlocalVectorAssembly(2)%p_Idofs(idofe,iel)) +&
                       DlocalData(2,idofe)
        end do

      end do ! iel











    end do ! IELset
    
    ! Release the local vector assembly structure
    call linf_releaseAssemblyData(rlocalVectorAssembly(1))
    call linf_releaseAssemblyData(rlocalVectorAssembly(2))

    ! Deallocate memory
    deallocate(Dxi2D,DpointsRef,IelementList,DsolVals,edgelength,normal)

  end subroutine














  !****************************************************************************
  
!<subroutine>  
  
  subroutine linf_getLocalEdgeNumbers(IedgeList,IlocalEdgeNumber,&
                                      rtriangulation)

!<description>

  ! For each global edge number in the list the local edge numbers is
  ! of the two adjacent elements are calculated.

!</description>

!<input>
 
  ! List of edges for which to find the local edge numbers
  integer, dimension(:), intent(in), pointer :: IedgeList
  
  ! The underlying triangulation
  type(t_triangulation), intent(in), pointer :: rtriangulation

  
!</input>

!<output>
  
  ! List of edges for which to find the local edge numbers
  integer, dimension(:,:), intent(out), pointer :: IlocalEdgeNumber
!</output>
  
!</subroutine>
  ! local variables
  integer :: nedge, iedge, IglobalEdgeNumber, ImaxedgePerElement, ilocal
  integer :: iel1, iel2
  integer, dimension(:,:), pointer :: p_IelementsAtEdge, p_IedgesAtElement
  
  
  ! Get pointers to the needed arrays from the triangulation
  call storage_getbase_int2D(rtriangulation%h_IelementsAtEdge,&
                             p_IelementsAtEdge)
  call storage_getbase_int2D(rtriangulation%h_IedgesAtElement,&
                             p_IedgesAtElement)
                             
  ! Get maximum number of edges per element in the triangulation
  ImaxedgePerElement = rtriangulation%NNEE
                             
  ! Get the number of edges
  nedge = size(IedgeList)
  
  ! Loop over all edges
  do iedge = 1, nedge
    
    ! Get the global edge number
    Iglobaledgenumber = IedgeList(iedge)
    
    ! Get the numbers of the elements adjacent to that edge
    iel1 = p_IelementsAtEdge(1,Iglobaledgenumber)
    iel2 = p_IelementsAtEdge(2,Iglobaledgenumber)
    
    ! Get the local edge number of the global edge in the first element adjacent to that edge
    do ilocal = 1,ImaxedgePerElement
      if (Iglobaledgenumber.eq.p_IedgesAtElement(ilocal,iel1)) exit
    end do
    IlocalEdgeNumber(1,iedge) = ilocal
    
    ! Test, if the second element exists (otherwise the edge is a boundary edge)
    ! If so, find its local edge number, too
    if (iel2.eq.0) then
      ! Second element doesn't exist
      ! Set local edge number as 0
      IlocalEdgeNumber(2,iedge) = 0
    else 
      ! Second element does exist, so do the same thing as with the first element
      do ilocal = 1,ImaxedgePerElement
        if (Iglobaledgenumber.eq.p_IedgesAtElement(ilocal,iel2)) exit
      end do
      IlocalEdgeNumber(2,iedge) = ilocal
    end if
  
  end do ! iedge



  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    !****************************************************************************
  
!<subroutine>  
  
   subroutine dg2gmv(rvector,extraNodesPerEdge)

!<description>

  ! Output a DG vector to gmv format

!</description>

!<input>
 
  ! The solution vector to output
  type(t_vectorScalar), intent(in) :: rvector
  
  ! Refinement level of the output grid (0 = No, n = n extra points on edge)
  integer, intent(in) :: extraNodesPerEdge
  
    
!</input>

!<output>
!</output>
  
!</subroutine>
  ! local variables
  
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! Space for the coordinates of the points on the reference element
  real(dp), dimension(:,:), allocatable :: drefCoords
  
  ! The corner points of the cells on the reference element
  integer, dimension(:,:), allocatable ::irefCornerNodesOfCell
  
  ! The coordinates of the real vertices (not on the reference element)
  real(dp), dimension(:,:), allocatable :: dNodeCoords
  
  ! Vertices at element of the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  
  ! Vertices at element of the triangulation
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  ! Maps element to corner vertices
  integer, dimension(:,:), allocatable :: iCornerNodesOfCell
  
  integer :: nnodesOnRef, i, j, ncellsOnRef, icell, nnodes, iel, ibaseC, ibaseN, NEL, inode, ncells, iunit
  
  real(dp) :: dx, dy
  
  integer(I32) :: ctrafoType
  
  real(DP), dimension(:,:), allocatable :: Djac
  
  real(DP), dimension(:), allocatable :: Ddetj
  
  real(DP), dimension(:), allocatable :: dnodeValues
  
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvector%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Get ointers to the data form the truangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)

  ! Get number of elements from the triangulation                               
  NEL = p_rtriangulation%NEL
  
  ! First calculate the number of nodes on the reference element
  nnodesOnRef = (2+extraNodesPerEdge)*(2+extraNodesPerEdge)

  ! Allocate space for the coordinates of the nodes on the reference element
  allocate(drefCoords(3,nnodesOnRef))

  ! Calculate the coordinates of the nodes on the reference element
  inode = 1
  do i=1,2+extraNodesPerEdge
    dy = -1.0_dp + (i-1)*2.0_dp/(1+extraNodesPerEdge)
    do j=1,2+extraNodesPerEdge
      dx = -1.0_dp + (j-1)*2.0_dp/(1+extraNodesPerEdge)
      drefCoords(1,inode) = dx
      drefCoords(2,inode) = dy
      drefCoords(3,inode) = 0.0_dp
      inode=inode+1
    end do
  end do


  ! First calculate the number of cells on the reference element
  ncellsOnRef = (1+extraNodesPerEdge)*(1+extraNodesPerEdge)

  ! Allocate space for the array taking the corner nodes of the cells
  ! on the reference element
  allocate(irefCornerNodesOfCell(4,ncellsOnRef))

  ! Calculate the array taking the corner nodes of the cells
  ! on the reference element
  icell = 1
  do i=1,1+extraNodesPerEdge
    do j=1,1+extraNodesPerEdge
      irefCornerNodesOfCell(1,icell) = (j-1)*(2+extraNodesPerEdge)+(i  )
      irefCornerNodesOfCell(2,icell) = (j-1)*(2+extraNodesPerEdge)+(i+1)
      irefCornerNodesOfCell(3,icell) = (j  )*(2+extraNodesPerEdge)+(i+1)
      irefCornerNodesOfCell(4,icell) = (j  )*(2+extraNodesPerEdge)+(i  )
      icell=icell+1
    end do
  end do


  ! Now calculate the total number of nodes to write to the gmv file
  nnodes = nnodesOnRef * NEL

  ! Allocate array for the coordinates of these nodes
  allocate(dNodeCoords(3,nnodes))
  
  ! Get type of transformation
  ctrafotype = p_rspatialDiscr%RelementDistr(1)%ctrafotype
  
  ! Allocate temp space for mapping
  allocate(Djac(4,nnodesOnRef))
  allocate(Ddetj(nnodesOnRef))

  ! Calculate the real coordinates of the nodes
  do iel = 1, NEL
  
    call trafo_calctrafo_mult (ctrafoType,nnodesOnRef,&
                             p_DvertexCoords(:,p_IverticesAtElement(:,iel)),&
                             drefCoords(1:2,:),Djac,Ddetj,&
                             dNodeCoords(1:2,(iel-1)*nnodesOnRef+1:(iel)*nnodesOnRef))

  end do
  
  ! Third coordinate is zero
  dNodeCoords(3,:) = 0.0_DP

  ! Deallocate temp space for mapping
  deallocate(Djac)
  deallocate(Ddetj)

  ! Calculate the total number of cells
  ncells = ncellsOnRef*NEL

  ! Allocate space for the array taking the corner nodes of the cells
  allocate(iCornerNodesOfCell(4,ncells))

  ! Calculate the array taking the corner nodes of the cells
  do iel=1,NEL
    ibaseC = (iel-1)*ncellsOnRef
    ibaseN = (iel-1)*nnodesOnRef
    do i=1,ncellsOnRef
      iCornerNodesOfCell(1,ibaseC+i) = irefCornerNodesOfCell(1,i)+ibaseN
      iCornerNodesOfCell(2,ibaseC+i) = irefCornerNodesOfCell(2,i)+ibaseN
      iCornerNodesOfCell(3,ibaseC+i) = irefCornerNodesOfCell(3,i)+ibaseN
      iCornerNodesOfCell(4,ibaseC+i) = irefCornerNodesOfCell(4,i)+ibaseN
    end do
  end do
  
  
  
  ! Evaluate the values of the solution vector in the points
  
  allocate(dnodeValues(nnodes))
  
  do iel = 1, NEL
    ibaseN = (iel-1)*nnodesOnRef
    call fevl_evaluate_mult1 (DER_FUNC, dnodeValues(ibaseN+1:ibaseN+nnodesOnRef),&
                              rvector, iel, drefCoords(1:2,:))
  end do
  
  
  
  
  
  
  ! ************ WRITE TO FILE PHASE *******************
  
  iunit = sys_getFreeUnit()
  
  open(iunit, file='./gmv/u2d.gmv')
  
  write(iunit,'(A)') 'gmvinput ascii'
  
  write(iunit,'(A,1X,I10)') 'nodes', nnodes
  do j=1,3
    do i=1,nnodes
      write(iunit,'(E15.8)') dNodeCoords(j,i)
    end do
  end do

  write(iunit,'(A,1X,I10)') 'cells',ncells
  do i=1,ncells
    write(iunit,'(A)')'quad 4'
    write(iunit,'(4I8)') iCornerNodesOfCell(1:4,i)
  end do

write (iunit,'(A)') 'variable'  
write (iunit,'(A,1X,I5)') trim('sol'),1
do i=1,nnodes
  write (iunit,'(E15.8)') dnodevalues(i)
end do

write (iunit,'(A)') 'endvars'


write (iunit,'(A)') 'endgmv'
close(iunit)


  deallocate(drefCoords)
  deallocate(irefCornerNodesOfCell)
  deallocate(dNodeCoords)
  deallocate(iCornerNodesOfCell)
  deallocate(dnodeValues)

end subroutine



    
end module
