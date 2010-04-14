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
    use feevaluation
    
  use stdoperators
  use genoutput
  use linearsolver
  use boundary
  use filtersupport
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use pprocerror
  use genoutput
  
  use ucd
  
    
  use dg2d_callback
    

  implicit none
    
  type t_dpPointer
		! Pointer to the double-valued matrix or vector data
    	real(DP), dimension(:), pointer :: p_Ddata
	end type 
	
	type t_additionalTriaData
	  ! Pointer to normal vectors of the triangulation (ndim,nedge)
	  real(DP), dimension(:,:), pointer :: p_Dnormals
	  ! Pointer to local edge numbers (2,nedge)
	  integer, dimension(:,:), pointer :: p_IlocalEdgeNumber
	  ! Lengths of the edges
	  real(DP), dimension(:), pointer :: p_DedgeLength
	  ! Midpoints of the elements
	  real(dp), dimension(:,:), pointer :: p_DmidPoints
	  ! Pointer to dx/dy
    real(DP), dimension(:,:), pointer:: p_Ddxdy
	end type 
	
	type t_profiler
	  integer :: ntimer, icurrenttimer
	  real(dp) :: dstarttime, dendtime, dlasttime
	  real(dp), dimension(:), allocatable :: Dtimers
	end type
	
	public :: linf_dg_buildVectorScalarEdge2d

contains

        
    

!****************************************************************************

!<subroutine>

  subroutine linf_dg_buildVectorScalarEdge2d (rform, ccubType, bclear,&
                                              rvectorScalar,&
                                              rvectorScalarSol,&
                                              raddTriaData,&
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
  
  ! Additional triangulation data
  type(t_additionalTriaData), intent(in) :: raddTriaData
   
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
  if (.not. associated(rvectorScalarSol%p_rspatialDiscr)) then
    call output_line('No discretisation associated!',&
        OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
    call sys_halt()
  end if

  ! The discretisation must provide a triangulation structure
  if (.not. associated(rvectorScalarSol%p_rspatialDiscr%p_rtriangulation)) then
    call output_line('No triangulation associated!',&
        OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
    call sys_halt()
  end if
  
  ! The discretisation must provide a boundary structure
  if (.not. associated(rvectorScalarSol%p_rspatialDiscr%p_rboundary)) then
    call output_line('No boundary associated!',&
        OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
    call sys_halt()
  end if

  ! Set pointers for quicker access
  p_rboundary => rvectorScalarSol%p_rspatialDiscr%p_rboundary
  p_rtriangulation => rvectorScalarSol%p_rspatialDiscr%p_rtriangulation
    
  ! Do we have a uniform triangulation? Would simplify a lot...
  if (rvectorScalarSol%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then 
    
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
     
      ! Initialise the vectorAssembly structures
      call linf_initAssembly(rvectorAssembly(1), rform,&
            rvectorScalarSol%p_rspatialDiscr%RelementDistr(1)%celement,&
            ccubType, LINF_NELEMSIM)
      call linf_initAssembly(rvectorAssembly(2), rform,&
            rvectorScalarSol%p_rspatialDiscr%RelementDistr(1)%celement,&
            ccubType, LINF_NELEMSIM)
            
      ! Assemble the data for all elements in this element distribution
      call linf_dg_assembleSubmeshVectorScalarEdge2d (rvectorAssembly,&
            rvectorScalar, rvectorScalarSol,&
            p_IedgeList(1:NMT),raddTriaData,&
            flux_dg_buildVectorScEdge2D_sim,&
            rcollection&
             )

          
      ! Release the assembly structure.
      call linf_doneAssembly(rvectorAssembly(1))
      call linf_doneAssembly(rvectorAssembly(2))

      ! Deallocate the edgelist
      deallocate(p_IedgeList)
            
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
              IedgeList, raddTriaData,&
              flux_dg_buildVectorScEdge2D_sim,&
              rcollection)

!<description>

  ! Assembles the vector entries for a submesh by integration over the given edges.

!</description>

!<input>

  ! List of edges where to assemble the linear form.
  integer, dimension(:), intent(in), target :: IedgeList
    
  ! The solution vector. Used to calculate the solution on the edges.
  type(t_vectorScalar), intent(in) :: rvectorSol
  
  ! Additional triangulation data
  type(t_additionalTriaData), intent(in) :: raddTriaData

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
    real(DP) :: domega1,domega2,daux1,daux2,dlen
    real(DP) :: dval1, dval2
    integer(I32) :: cevaluationTag
    type(t_linfVectorAssembly), dimension(2), target :: rlocalVectorAssembly
    type(t_domainIntSubset), dimension(2) :: rintSubset
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP), dimension(:,:,:,:), pointer :: p_Dcoefficients
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
    integer, dimension(:,:), allocatable, target :: IelementList
    
    integer(i32) :: icoordSystem
    
    ! Chooses the element on side ... of the edge
    integer :: iside
    
    ! Number of edges
    integer :: NMT
    
    ! Vertices per element
    integer :: NVE
    
    integer :: ive
    
    ! Pointer to Ielementsatedge in the triangulation
    integer, dimension(:,:), pointer :: p_IelementsAtEdge
    
    ! Pointer to IverticesAtEdge in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    ! Pointer to the vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    
    ! Pointer to IverticesAtEelement in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    
!    ! Array for the solution values in the cubature points
!    real(DP), dimension(:,:,:), allocatable :: DsolVals
         
    ! Array for the normal vectors
    real(DP), dimension(:,:), allocatable :: normal
    
    ! Temp variables for the coordinates of the vertices
    real(DP) :: dxl1, dxl2, dyl1, dyl2
    
    ! The Jacobian matrix of the mapping for each point.
    ! DIMENSION(number of entries in the matrix,npointsPerEl,nelements)
    real(DP), dimension(:,:,:), allocatable :: Djac
  
    ! Jacobian determinants of the mapping for all the points from the
    ! reference element to the real element.
    ! DIMENSION(npointsPerEl,nelements)
    real(DP), dimension(:,:), allocatable :: Ddetj
  
    ! Array receiving the coordinates of the points in DpointsRef,
    ! mapped from the reference element to the real element.
    ! If not specified, they are not computed.
    ! DIMENSION(#space dimensions,npointsPerEl,nelements)
    real(DP), dimension(:,:,:), allocatable :: DpointsReal
    
    ! The coordinates of the corner edges of the elements for the transformation
    ! DIMENSION(#space dimensions,NVE,nelements)
    real(DP), dimension(:,:,:), allocatable :: Dcoords
    
    integer(I32) :: ctrafotype
    
        ! Space for the values of the flux function
    real(DP), dimension(:,:,:), allocatable :: DfluxValues
  
  
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
    
    ! Get pointers to vertices at elements
    call storage_getbase_int2D(&
          rvector%p_rspatialDiscr%p_rtriangulation%h_IverticesAtElement,&
          p_IverticesAtElement)   
    
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
    NVE = rvectorAssembly(1)%NVE
    
!    ! Allocate space for the solution values in the cubature points
!    allocate(DsolVals(2,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock))
    
!    ! Allocate space for the jacobi matrices in the cubature points
!    allocate(Djac(4,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock))
!    
!    ! Allocate space for the determinants of the jacobi matrices in the cubature points
!    allocate(Ddetj(ncubp,rlocalVectorAssembly(1)%nelementsPerBlock))
!    
!    ! Allocate space for the integration points on the real elements
!    allocate(DpointsReal(ndim2d,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock))

!    ! Allocate space for normal vectors
!    allocate(normal(2,min(size(IedgeList),rlocalVectorAssembly(1)%nelementsPerBlock)))
!    
!    ! Allocate space for edge length
!    allocate(edgelength(min(size(IedgeList),rlocalVectorAssembly(1)%nelementsPerBlock)))
    
!    ! The coordinates of the corner edges of the elements for the transformation
!    allocate(Dcoords(ndim2d,NVE,rlocalVectorAssembly(1)%nelementsPerBlock))
      
  
  
  
  
   
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
    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock,2))
    
    ! Allocate space for the flux variables DIM(nvar,ialbet,ncubp,elementsperblock)
    allocate(DfluxValues(1,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock))

    ! Get the type of coordinate system
    icoordSystem = elem_igetCoordSystem(rlocalVectorAssembly(1)%celement)
    
    ! Get type of transformation
    ctrafotype = elem_igetTrafoType(rlocalVectorAssembly(1)%celement)

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
        call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(1,Iedgelist(IELset+iel-1)), &
            ncubp, Dxi1D, Dxi2D(:,:,1,iel))
        call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(2,Iedgelist(IELset+iel-1)), &
            ncubp, Dxi1D, Dxi2D(:,:,2,iel))
      end do
     
      ! Transpose the coordinate array such that we get coordinates we
      ! can work with.
      do iside = 1,2
        do iel = 1,IELmax-IELset+1
          do icubp = 1,ncubp
            do k = 1,ubound(DpointsRef,1)
              DpointsRef(k,icubp,iel,iside) = Dxi2D(icubp,k,iside,iel)
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
          DpointsRef=DpointsRef(:,:,:,1))
      call elprep_prepareSetForEvaluation (&
          rlocalVectorAssembly(2)%revalElementSet,&
          cevaluationTag, rvector%p_rspatialDiscr%p_rtriangulation, &
          IelementList(3,IELset:IELmax), rlocalVectorAssembly(2)%ctrafoType, &
          DpointsRef=DpointsRef(:,:,:,2))
          
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
    
!      ! Now that we have the basis functions, we want to have the function values.
!      ! We get them by multiplying the FE-coefficients with the values of the
!      ! basis functions and summing up.
!      do iel = 1,IELmax-IELset+1      
!        do icubp = 1,ncubp
!          ! Calculate the value in the point
!          dval1 = 0.0_DP
!          dval2 = 0.0_DP
!          do idofe = 1,indof
!            dval1 = dval1 + &
!                   p_DdataSol(rlocalVectorAssembly(1)%p_Idofs(idofe,iel)) &
!                   * rlocalVectorAssembly(1)%p_Dbas(idofe,DER_FUNC,icubp,iel)
!            dval2 = dval2 + &
!                   p_DdataSol(rlocalVectorAssembly(2)%p_Idofs(idofe,iel)) &
!                   * rlocalVectorAssembly(2)%p_Dbas(idofe,DER_FUNC,icubp,iel)
!          end do
!          ! Save the value in the point
!          DsolVals(1,icubp,iel) = dval1
!          DsolVals(2,icubp,iel) = dval2
!        end do
!      end do




      
!     ! Set values at boundary
!     do iel = 1,IELmax-IELset+1
!      if(IelementList(2,IELset+iel-1).eq.0) then
!        DsolVals(2,1:ncubp,iel) = 0.0_DP
!      end if
!      
!      end do
     
      
     ! ---------------------- Get values of the flux function --------------
     
     
!      ! Now it is time to call our coefficient function to calculate the
!      ! function values in the cubature points:
!      if (present(fcoeff_buildVectorScBdr2D_sim)) then
!        call domint_initIntegrationByEvalSet (p_revalElementSet, rintSubset)
!        rintSubset%ielementDistribution = 0
!        rintSubset%ielementStartIdx = IELset
!        rintSubset%p_Ielements => IelementList(IELset:IELmax)
!        rintSubset%p_IdofsTrial => p_Idofs
!        rintSubset%celement = rlocalVectorAssembly%celement
!        call fcoeff_buildVectorScBdr2D_sim (rvector%p_rspatialDiscr,&
!            rlocalVectorAssembly%rform,  IELmax-IELset+1, ncubp,&
!            p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!            ibdc, DpointsPar(:,1:IELmax-IELset+1),&
!            p_Idofs, rintSubset, &
!            p_Dcoefficients(:,:,1:IELmax-IELset+1), rcollection)
!        call domint_doneIntegration (rintSubset)
!      else
!        p_Dcoefficients(:,:,1:IELmax-IELset+1) = 1.0_DP
!      end if





!      ! Fill the corner coordinates of the elements
!      do iel = 1,IELmax-IELset+1
!        do ive = 1, NVE
!          Dcoords(1:ndim2d,ive,iel)=&
!                   p_DvertexCoords(1:ndim2d,p_IverticesAtElement(ive,IelementList(1,IELset+iel-1)))
!        end do
!      end do
!
!
!      ! The numerical flux function needs the x- and y- values
!      ! So we have to call the mapping from the reference- to the real element
!      call trafo_calctrafo_sim (ctrafoType,IELmax-IELset+1,ncubp,Dcoords,&
!                                DpointsRef(1:ndim2d,:,1:IELmax-IELset+1,1),Djac(1:4,1:ncubp,1:IELmax-IELset+1),&
!                                Ddetj(1:ncubp,1:IELmax-IELset+1),&
!                                DpointsReal(1:ndim2d,1:ncubp,1:IELmax-IELset+1))










    ! If the flux function needs other, than just the function values from the solution
    ! (for example the derivatives), we will give an evalElementSet to it
    ! This is filled here
    
    call domint_initIntegrationByEvalSet (rlocalVectorAssembly(1)%revalElementSet, rintSubset(1))
    call domint_initIntegrationByEvalSet (rlocalVectorAssembly(2)%revalElementSet, rintSubset(2))
        !rintSubset(1)%ielementDistribution = 0
        rintSubset(1)%ielementStartIdx = IELset
        rintSubset(1)%p_Ielements => IelementList(1,IELset:IELmax)
        rintSubset(1)%p_IdofsTrial => rlocalVectorAssembly(1)%p_Idofs
        rintSubset(1)%celement = rlocalVectorAssembly(1)%celement
        !rintSubset(2)%ielementDistribution = 0
        rintSubset(2)%ielementStartIdx = IELset
        rintSubset(2)%p_Ielements => IelementList(2,IELset:IELmax)
        rintSubset(2)%p_IdofsTrial => rlocalVectorAssembly(2)%p_Idofs
        rintSubset(2)%celement = rlocalVectorAssembly(2)%celement





      
      
      call flux_dg_buildVectorScEdge2D_sim (&
            DfluxValues(:,:,1:IELmax-IELset+1),&
!            DsolVals(:,:,1:IELmax-IELset+1),&
            IelementList(2,IELset:IELmax),&
            raddTriaData%p_Dnormals(:,Iedgelist(IELset:IELmax)),&
            !DpointsReal(1:ndim2d,1:ncubp,1:IELmax-IELset+1),&
            rintSubset,&
            rcollection )
            
      
      call domint_doneIntegration (rintSubset(1))
      call domint_doneIntegration (rintSubset(2))
      
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
        dlen = 0.5_DP*raddTriaData%p_Dedgelength(Iedgelist(IELset+iel-1))

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature
          ! formula in that cubature point.

          domega1 = dlen * rlocalVectorAssembly(1)%p_Domega(icubp)
          domega2 = dlen * rlocalVectorAssembly(1)%p_Domega(icubp)

          
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
            daux1 = domega1 * DfluxValues(ialbet,icubp        ,iel)
            daux2 = domega2 * DfluxValues(ialbet,ncubp-icubp+1,iel) *(-1.0_dp)
           
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
              
!              if(IelementList(2,IELset+iel-1).ne.0) then
!                DlocalData(2,idofe) = DlocalData(2,idofe)+&
!                                      rlocalVectorAssembly(2)%p_Dbas(idofe,ia,icubp,iel)*daux2
!              end if
                    
              DlocalData(2,idofe) = DlocalData(2,idofe)+&
                                    rlocalVectorAssembly(2)%p_Dbas(idofe,ia,icubp,iel)*daux2* real(min(1,IelementList(2,IELset+iel-1)))
              
             
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
    deallocate(Dxi2D,DpointsRef,IelementList)!,DsolVals,edgelength,normal,Djac,Ddetj,DpointsReal,Dcoords)
    deallocate(DfluxValues)

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
  
   subroutine dg2gmv(rvector,extraNodesPerEdge,sofile,ifilenumber)

!<description>

  ! Output a DG vector to gmv format

!</description>

!<input>
 
  ! The solution vector to output
  type(t_vectorScalar), intent(in) :: rvector
  
  ! Refinement level of the output grid (0 = No, n = n extra points on edge)
  integer, intent(in) :: extraNodesPerEdge
  
  ! Name of output file
  character (LEN=SYS_STRLEN), intent(in) :: sofile
  
  ! Filenumber
  integer, intent(in) :: ifilenumber
  
    
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
  
  character (LEN=10) :: sfilenumber
  
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvector%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Get pointers to the data form the truangulation
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
  
  write(sfilenumber,'(i0)') ifilenumber

  if (ifilenumber>-1) then
    open(iunit, file=trim(sofile) // trim(sfilenumber) // '.gmv')
  else
    open(iunit, file=trim(sofile) // '.gmv')
  end if
  
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





    !****************************************************************************
  
!<subroutine>  
  
   subroutine dg2vtk(rvector,extraNodesPerEdge,sofile,ifilenumber)

!<description>

  ! Output a DG vector to gmv format

!</description>

!<input>
 
  ! The solution vector to output
  type(t_vectorScalar), intent(in) :: rvector
  
  ! Refinement level of the output grid (0 = No, n = n extra points on edge)
  integer, intent(in) :: extraNodesPerEdge
  
  ! Name of output file
  character (LEN=SYS_STRLEN), intent(in) :: sofile
  
  ! Filenumber
  integer, intent(in) :: ifilenumber
  
    
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
  
  character (LEN=10) :: sfilenumber
  
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvector%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Get pointers to the data form the truangulation
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
  
  write(sfilenumber,'(i0)') ifilenumber

  if (ifilenumber>-1) then
    open(iunit, file=trim(sofile) // trim(sfilenumber) // '.vtk')
  else
    open(iunit, file=trim(sofile) // '.vtk')
  end if
  
  write(iunit, '(A)') "# vtk DataFile Version 2.0"
  write(iunit, '(A)') "Generated by FEATFlow 2.x"
  write(iunit, '(A)') "ASCII"

  write(iunit,'(A)') 'DATASET UNSTRUCTURED_GRID'
  write(iunit,'(A,I10,A)') "POINTS", nnodes, " double"
    do i=1,nnodes
      write(iunit,'(3E16.7)') dNodeCoords(1:3,i)
    end do
  
  write(iunit,'(A,2I10)') "CELLS",ncells,ncells*5
  do i=1,ncells
    write(iunit,'(5I8)') 4,iCornerNodesOfCell(1:4,i)-1
  end do
  
  write(iunit,'(A,I10)') "CELL_TYPES",ncells
  do i=1,ncells
    write(iunit,'(I3)') 9
  end do
  
  write (iunit,'(A,I10)') "POINT_DATA", nnodes
  write (iunit,'(A)') "SCALARS scalars double 1"
  write (iunit,'(A)') "LOOKUP_TABLE default"
  do i=1,nnodes
    write (iunit,'(ES16.8E3)') dnodevalues(i)
  end do

  close(iunit)


  deallocate(drefCoords)
  deallocate(irefCornerNodesOfCell)
  deallocate(dNodeCoords)
  deallocate(iCornerNodesOfCell)
  deallocate(dnodeValues)

end subroutine



  !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_linearLimiter (rvector)

!<description>

  ! Limits the linear part of a dg_T1 element vector.

!</description>

!<input>
!</input>

!<inputoutput>

  ! A vector to limit
  type(t_vectorScalar), intent(inout) :: rvector  
  
!</inputoutput>
  
!</subroutine>

  ! local variables, used by all processors
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: indof, NEL, iel, NVE, ivt, NVT
    
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! The coordinates of the points in which to evaluate the solution vector
  real(dp), dimension(2,17) :: Dpoints
  
  ! The list of elements, in which these points can be found
  integer, dimension(17) :: Ielements
  
  ! The values of the solution vector in the points
  real(dp), dimension(17) :: Dvalues
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  real(dp), dimension(:), allocatable :: duimax, duimin
  
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  
  real(dp) :: xc,yc
  integer :: iidx, nvert, ivert, ineighbour, ineighElm
  integer, dimension(4) :: IhomeIndex
  
  real(dp) :: dui, ddu, dalpha, dalphatemp, duc
  
  integer, dimension(3) :: IdofGlob
  
  integer :: NVBD

  ! Get pointer to the solution data
  call lsyssc_getbase_double (rvector,p_Ddata)
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvector%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Get pointer to the data of the vector
  
  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  ! Number of vertives at boundary
  NVBD = p_rtriangulation%NVBD
  
  ! Get pointers to the data form the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
                               p_IelementsAtVertex) 
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
                               p_IelementsAtVertexIdx)
  call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary ,&
                               p_IverticesAtBoundary)                               
                               
                               
   ! Set pointer to coordinate vector
   call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

   ! Set pointer to vertices at element
   call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
        
  NVT = p_rtriangulation%NVT
  
  allocate(duimax(NVT),duimin(NVT))
  
  duimax= -SYS_MAXREAL
  duimin=  SYS_MAXREAL
  
  do iel = 1, NEL
    
    ! Get global DOFs of the element
    call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    duc = p_Ddata(IdofGlob(1))
    
    ! elem_igetNVE(celement)
    NVE = 4
    
    do ivt = 1, NVE
      nvt = p_IverticesAtElement(ivt,iel)
      duimax(nvt) = max(duc,duimax(nvt))
      duimin(nvt) = min(duc,duimin(nvt))
    end do
    
  end do
  
!  do ivt = 1, NVBD
!    duimax(p_IverticesAtBoundary(ivt)) = 100000000.0_DP
!    duimin(p_IverticesAtBoundary(ivt)) = -100000000.0_DP
!  end do

                               
  do iel = 1, NEL
  
    ! Get number of corner vertices
    ! elem_igetNVE(celement)
    NVE = 4
    
    ! Get midpoint of the element
    xc = &
         (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(4,iel)))/4.0_dp

    yc = &
         (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(4,iel)))/4.0_dp
          
    ! The first point we want to evaluate the solution in, is the midpoint of the element
    Dpoints(1,1) = xc
    Dpoints(2,1) = yc
    Ielements(1) = iel
    
    ! Initialise the limiting factor
    dalpha = 1.0_dp
    
    ! Now start to set the points, where to evaluate the solution
    
    ! Loop over the vertices of the element
    do ivert = 1, NVE
		
      nvert = p_IverticesAtElement(ivert, iel)
      
      ! The second point we want to evaluate the solution in, is in the corner of the mother element
      Dpoints(1,1+ivert) = p_DvertexCoords(1,nvert)
      Dpoints(2,1+ivert) = p_DvertexCoords(2,nvert)
      Ielements(1+ivert) = iel
	  end do
	  
    ! Evaluate the solution
    call fevl_evaluate (DER_FUNC, Dvalues(1:5), rvector, Dpoints(1:2,1:5), &
                          Ielements(1:5))
      
     ! Start calculating the limiting factor
     duc = Dvalues(1)
      
    do ivert = 1, NVE  
      dui = Dvalues(1+ivert)
      ddu = dui-duc
      nvert = p_IverticesAtElement(ivert, iel)
            
      ! Find the maximum/minimum value of the solution in the centroids
      ! of all elements containing this vertex
      if (ddu > 0.0_dp) then
        dalphatemp = min(1.0_dp, (duimax(nvert)-duc)/ddu)
      elseif (ddu < 0.0_dp) then
        dalphatemp = min(1.0_dp, (duimin(nvert)-duc)/ddu)
      else ! (dui==duc)
        dalphatemp = 1.0_dp
      end if
      
      dalpha = min(dalphatemp,dalpha)
      

    end do ! ivert
    
    
    ! Now we have the limitingfactor dalpha
    ! We need to multiply it with the corresponding DOFs
    
    ! Get global DOFs of the element
    call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    ! Multiply the linear part of the solution vector with the correction factor
    p_Ddata(IdofGlob(2:3)) = p_Ddata(IdofGlob(2:3))*dalpha
    
  end do ! iel
  
  deallocate(duimax,duimin)

  
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  subroutine dg_proj2steady(rsolBlock,rtriangulation, rboundary)
  
  
  
  
  
  
  
  
  
      ! An object for saving the domain:
    type(t_boundary), intent(in) :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), intent(in) :: rtriangulation
    
    type(t_vectorBlock),intent(in), target :: rsolBlock




    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrixMC

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock) :: rsolSteadyBlock,rrhsBlock,rtempBlock
    
    type(t_vectorScalar) :: rrhs,rsolSteady,rtemp

    ! A solver node that accepts parameters for the linear solver    
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices
    
    real(DP), dimension(:), pointer :: p_Ddata
    
    type(t_collection) :: rcollection
    
    type(t_filterChain), dimension(1), target :: RfilterChain
    
    type(t_ucdExport) :: rexport
    
    integer :: ielementtype
    
    integer :: ierror
  
  
  
   ielementtype = EL_Q1
  
      ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                 rtriangulation, rboundary)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   ielementType,CUB_G5X5,rtriangulation, rboundary)
                 
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrixMC) !,&
                                     !rdiscretisation%RspatialDiscr(1),&
                                     !BILF_MATC_EDGEBASED)
    
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 2D.

    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    
    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = 1.0 
    rform%Dcoefficients(2)  = 1.0 

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical
    ! data.
    call bilf_buildMatrixScalar (rform,.true.,rmatrixMC,coeff_Laplace_2D)
    
    
    
    
    
    
    
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rrhs ,.true.,ST_DOUBLE)
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rsolSteady ,.true.,ST_DOUBLE)
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rtemp ,.true.,ST_DOUBLE)
                                 

    call lsysbl_createMatFromScalar (rmatrixMC,rmatrixBlock,rdiscretisation)
    call lsysbl_createVecFromScalar (rrhs,rrhsBlock,rdiscretisation)
    call lsysbl_createVecFromScalar (rsolSteady,rsolSteadyBlock,rdiscretisation)
    
    call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
    !call lsysbl_createVecFromScalar (rtemp,rtempBlock,rdiscretisation)
    
    
    
    
    rcollection%p_rvectorQuickAccess1 => rsolBlock
    
    
    
    
    
    ! Now set the initial conditions via L2 projection
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC2D
    call linf_buildVectorScalar2 (rlinform, .true., rrhs,&
                                  coeff_Steadyproj, rcollection)
                                  
                                  
                                  
                                  
                                  
    nullify(p_rpreconditioner)
    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 0
    
    ! Attach the system matrix to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation 
    ! routine to attach all these matrices.
    ! Remark: Do not make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
    ! This does not work on all compilers, since the compiler would have
    ! to create a temp array on the stack - which does not always work!
    Rmatrices = (/rmatrixBlock/)
    call linsol_setMatrices(p_RsolverNode,Rmatrices)
    
        ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initData (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
                                  
    call linsol_solveAdaptively (p_rsolverNode,rsolSteadyBlock,rrhsBlock,rtempBlock)
    
    
    
  
    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       './gmv/steady.gmv')
    
    call lsyssc_getbase_double (rsolSteadyBlock%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    
    
    
    
    call lsysbl_releaseVector (rrhsBlock)
    call lsysbl_releaseVector (rTempBlock)
    call lsysbl_releaseVector (rsolSteadyBlock)
    call lsyssc_releaseVector (rrhs)
    call lsyssc_releaseVector (rsolsteady)
    call lsyssc_releaseVector (rtemp)
    
    call lsysbl_releaseMatrix (rmatrixBlock)
    call lsyssc_releaseMatrix (rmatrixMC)
    
    
    call spdiscr_releaseBlockDiscr(rdiscretisation)


  end subroutine
  
  
  
  
  
  
  
  
  
  
    !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_quadraticLimiter (rvector)

!<description>

  ! Limits the linear part of a dg_T1 element vector.

!</description>

!<input>
!</input>

!<inputoutput>

  ! A vector to limit
  type(t_vectorScalar), intent(inout) :: rvector  
  
!</inputoutput>
  
!</subroutine>

  ! local variables, used by all processors
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: indof, NEL, iel, NVE, ivt, NVT
    
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! The coordinates of the points in which to evaluate the solution vector
  real(dp), dimension(2,17) :: Dpoints
  
  ! The list of elements, in which these points can be found
  integer, dimension(17) :: Ielements
  
  ! The values of the solution vector in the points
  real(dp), dimension(17) :: Dvalues
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  real(dp), dimension(:), allocatable :: duimax, duimin
  
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  
  real(dp) :: xc,yc
  integer :: iidx, nvert, ivert, ineighbour, ineighElm
  integer, dimension(4) :: IhomeIndex
  
  real(dp) :: dui, ddu, dalphatemp, duc
  
  integer, dimension(6) :: IdofGlob
  
  integer, dimension(5) :: Isep
  
  integer :: NVBD, ilim, ideriv
  
  
  real(DP), dimension(:,:), allocatable :: Dalpha

  ! Get pointer to the solution data
  call lsyssc_getbase_double (rvector,p_Ddata)
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvector%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Get pointer to the data of the vector
  
  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  ! Number of vertives at boundary
  NVBD = p_rtriangulation%NVBD
  
  ! Get pointers to the data form the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
                               p_IelementsAtVertex) 
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
                               p_IelementsAtVertexIdx)
  call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary ,&
                               p_IverticesAtBoundary)                               
                               
                               
   ! Set pointer to coordinate vector
   call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

   ! Set pointer to vertices at element
   call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
        
  NVT = p_rtriangulation%NVT
  
  allocate(Duimax(NVT),Duimin(NVT),Dalpha(3,NEL))
  
  Dalpha = 1.0_dp
  
  ! Now limit the x- and y- derivative and finally the linear part
  
  do ilim = 1,3
  
    
  
    select case (ilim)
      case (1)
        ideriv = DER_DERIV_X
      case (2)
        ideriv = DER_DERIV_Y
      case (3)
        ideriv = DER_FUNC
    end select
  
  
  duimax= -SYS_MAXREAL
  duimin=  SYS_MAXREAL
  
  do iel = 1, NEL
    
    ! Get number of corner vertices
    ! elem_igetNVE(celement)
    NVE = 4
    
    ! Get midpoint of the element
    xc = &
         (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(4,iel)))/4.0_dp

    yc = &
         (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(4,iel)))/4.0_dp
          
    ! The first point we want to evaluate the solution in, is the midpoint of the element
    Dpoints(1,1) = xc
    Dpoints(2,1) = yc
    Ielements(1) = iel
    
    ! Evaluate the solution
    call fevl_evaluate (ideriv, Dvalues(1:1), rvector, Dpoints(1:2,1:1), &
                          Ielements(1:1))
                          
    duc = Dvalues(1)
    
    if (ilim == 3) then
      call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
      duc = p_Ddata(IdofGlob(1))
    end if
    
    do ivt = 1, NVE
      nvt = p_IverticesAtElement(ivt,iel)
      duimax(nvt) = max(duc,duimax(nvt))
      duimin(nvt) = min(duc,duimin(nvt))
    end do
    
    
    
  end do
  

                               
  do iel = 1, NEL
  
    ! Get number of corner vertices
    ! elem_igetNVE(celement)
    NVE = 4
    
    ! Get midpoint of the element
    xc = &
         (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(4,iel)))/4.0_dp

    yc = &
         (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(4,iel)))/4.0_dp
          
    ! The first point we want to evaluate the solution in, is the midpoint of the element
    Dpoints(1,1) = xc
    Dpoints(2,1) = yc
    Ielements(1) = iel
   
    ! Now start to set the points, where to evaluate the solution
    
    ! Loop over the vertices of the element
    do ivert = 1, NVE
		
      nvert = p_IverticesAtElement(ivert, iel)
      
      ! The second point we want to evaluate the solution in, is in the corner of the mother element
      Dpoints(1,1+ivert) = p_DvertexCoords(1,nvert)
      Dpoints(2,1+ivert) = p_DvertexCoords(2,nvert)
      Ielements(1+ivert) = iel
	  end do
	  
    ! Evaluate the solution
    call fevl_evaluate (ideriv, Dvalues(1:5), rvector, Dpoints(1:2,1:5), &
                          Ielements(1:5))
      
     ! Start calculating the limiting factor
     duc = Dvalues(1)
     
    if (ilim == 3) then
      call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
      duc = p_Ddata(IdofGlob(1))
    end if
    
!    do ivt = 1, NVE
!      nvt = p_IverticesAtElement(ivt,iel)
!      duimax(nvt) = max(duc,duimax(nvt))
!      duimin(nvt) = min(duc,duimin(nvt))
!    end do
      
    do ivert = 1, NVE  
      dui = Dvalues(1+ivert)
      
      if (ilim == 3) then
        call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
        dui = p_Ddata(IdofGlob(1))
        
        select case (ivert)
        case(1)
        dui = dui - p_Ddata(IdofGlob(2)) - p_Ddata(IdofGlob(3))
        case(2)
        dui = dui + p_Ddata(IdofGlob(2)) - p_Ddata(IdofGlob(3))
        case(3)
        dui = dui + p_Ddata(IdofGlob(2)) + p_Ddata(IdofGlob(3))
        case(4)
        dui = dui - p_Ddata(IdofGlob(2)) + p_Ddata(IdofGlob(3))
        end select
      end if
      
      ddu = dui-duc
      nvert = p_IverticesAtElement(ivert, iel)
            
      ! Find the maximum/minimum value of the solution in the centroids
      ! of all elements containing this vertex
      if (ddu > 0.0_dp) then
        dalphatemp = min(1.0_dp, (duimax(nvert)-duc)/ddu)
      elseif (ddu < 0.0_dp) then
        dalphatemp = min(1.0_dp, (duimin(nvert)-duc)/ddu)
      else ! (dui==duc)
        dalphatemp = 1.0_dp
      end if
      
      Dalpha(ilim,iel) = min(dalphatemp,Dalpha(ilim,iel))
      

    end do ! ivert
        
  end do ! iel
  
  
  
  
  select case (ilim)
    case (2)
  
      ! Now we have the limitingfactors for the quadratic part in Dalpha(1:2,:)
      ! Get the Minimum and multiply it with the corresponding DOFs
      do iel = 1, NEL  
  
        Dalpha(1,iel) = min(Dalpha(1,iel),Dalpha(2,iel))
      
        ! Get global DOFs of the element
        call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    
        ! Multiply the linear part of the solution vector with the correction factor
        p_Ddata(IdofGlob(4:6)) = p_Ddata(IdofGlob(4:6))*Dalpha(1,iel)
    
      end do ! iel
    
    case (3)
      do iel = 1, NEL  
  
        Dalpha(3,iel) = max(Dalpha(1,iel),Dalpha(3,iel))
      
        ! Get global DOFs of the element
        call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    
        ! Multiply the linear part of the solution vector with the correction factor
        p_Ddata(IdofGlob(2:3)) = p_Ddata(IdofGlob(2:3))*Dalpha(3,iel)
    
      end do ! iel
      
    end select
  
  end do ! ilim
  
  
  
  
  
  
  deallocate(duimax,duimin,dalpha)

  
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  !****************************************************************************
  
!<subroutine>  
  
  subroutine addTriaData(rtriangulation, raddTriaData)

!<description>

  ! For each global edge number in the list the local edge numbers
  ! of the two adjacent elements are calculated.
  ! And the normal vectors of the first element adjacent to each edge.

!</description>

!<input>

  ! The underlying triangulation
  type(t_triangulation), intent(in) :: rtriangulation
  
!</input>

!<output>
  
  ! The additional triangulation data
  type(t_additionalTriaData), intent(out):: raddTriaData
!</output>
  
!</subroutine>
  ! local variables
  integer :: nedge, iedge, IglobalEdgeNumber, ImaxedgePerElement, ilocal
  integer :: iel1, iel2, iel, NEL
  integer, dimension(:,:), pointer :: p_IelementsAtEdge, p_IedgesAtElement, p_IverticesAtElement
  ! Pointer to the vertex coordinates
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  real(dp) :: dxl1, dxl2, dyl1, dyl2
  ! Pointer to IverticesAtEdge in the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtEdge
  real(dp) :: dedgeLength
  integer :: ncorners
  
  
  ! Get pointers to the needed arrays from the triangulation
  call storage_getbase_int2D(rtriangulation%h_IelementsAtEdge,&
                             p_IelementsAtEdge)
  call storage_getbase_int2D(rtriangulation%h_IedgesAtElement,&
                             p_IedgesAtElement)
  call storage_getbase_int2D(rtriangulation%h_IverticesAtElement,&
                             p_IverticesAtElement)
                             
  ! Get maximum number of edges per element in the triangulation
  ImaxedgePerElement = rtriangulation%NNEE
  
  ! Get number of elements
  NEL = rtriangulation%NEL
                             
  ! Get the number of edges
  nedge = rtriangulation%NMT
  
  ! Allocate space for local edge numbers
  allocate(raddTriaData%p_IlocalEdgeNumber(2,nedge))
  
  ! Loop over all edges
  do iedge = 1, nedge
    
    ! Get the global edge number
    Iglobaledgenumber = iedge !IedgeList(iedge)
    
    ! Get the numbers of the elements adjacent to that edge
    iel1 = p_IelementsAtEdge(1,Iglobaledgenumber)
    iel2 = p_IelementsAtEdge(2,Iglobaledgenumber)
    
    ! Get the local edge number of the global edge in the first element adjacent to that edge
    do ilocal = 1,ImaxedgePerElement
      if (Iglobaledgenumber.eq.p_IedgesAtElement(ilocal,iel1)) exit
    end do
    raddTriaData%p_IlocalEdgeNumber(1,iedge) = ilocal
    
    ! Test, if the second element exists (otherwise the edge is a boundary edge)
    ! If so, find its local edge number, too
    if (iel2.eq.0) then
      ! Second element doesn't exist
      ! Set local edge number as 0
      raddTriaData%p_IlocalEdgeNumber(2,iedge) = 0
    else 
      ! Second element does exist, so do the same thing as with the first element
      do ilocal = 1,ImaxedgePerElement
        if (Iglobaledgenumber.eq.p_IedgesAtElement(ilocal,iel2)) exit
      end do
      raddTriaData%p_IlocalEdgeNumber(2,iedge) = ilocal
    end if
  
  end do ! iedge
  
  
  
  ! *** Calculate normal vectors ***
  
  ! Get pointer to the vertex coordinates
  call storage_getbase_double2D(rtriangulation%h_DvertexCoords,&
          p_DvertexCoords)
  
  ! Get pointer to vertices at edge
    call storage_getbase_int2D(&
          rtriangulation%h_IverticesAtEdge,&
          p_IverticesAtEdge)
	  
	! Allocate space for normal vectors
  allocate(raddTriaData%p_Dnormals(2,nedge))
    
  ! Allocate space for edge length
  allocate(raddTriaData%p_DedgeLength(nedge))

  do iedge = 1, nedge
    ! Calculate the length of the edge 
    dxl1=p_DvertexCoords(1,p_IverticesAtEdge(1,iedge))
    dyl1=p_DvertexCoords(2,p_IverticesAtEdge(1,iedge))
    dxl2=p_DvertexCoords(1,p_IverticesAtEdge(2,iedge))
    dyl2=p_DvertexCoords(2,p_IverticesAtEdge(2,iedge))
    dedgelength = sqrt((dxl1-dxl2)*(dxl1-dxl2)+(dyl1-dyl2)*(dyl1-dyl2))
    raddTriaData%p_DedgeLength(iedge) = dedgeLength
         
    ! Calculate the normal vector to the element at this edge
    raddTriaData%p_Dnormals(1,iedge) = (dyl2-dyl1)/dedgeLength
    raddTriaData%p_Dnormals(2,iedge) = (dxl1-dxl2)/dedgeLength
  end do
  
  
  ! *** Calculate midpoints of the elements ***
  
  ! Allocate space for midpoints
  allocate(raddTriaData%p_DmidPoints(2,NEL))
  
  do iel = 1, NEL
    raddTriaData%p_DmidPoints(1,iel) = &
         (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(4,iel)))*0.25_dp

    raddTriaData%p_DmidPoints(2,iel) = &
         (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(4,iel)))*0.25_dp
   end do
   
   
  
  ! *** Calculate dx and dy ***
  
  ! Allocate space for dx and dy
  allocate(raddTriaData%p_Ddxdy(2,NEL))
  
  ncorners = 4
  
  do iel = 1, NEL
    raddTriaData%p_Ddxdy(1,iel) = &
          maxval(p_DvertexCoords(1,p_IverticesAtElement(1:ncorners,iel))) - &
          minval(p_DvertexCoords(1,p_IverticesAtElement(1:ncorners,iel)))
         
    raddTriaData%p_Ddxdy(2,iel) = &
          maxval(p_DvertexCoords(2,p_IverticesAtElement(1:ncorners,iel))) - &
          minval(p_DvertexCoords(2,p_IverticesAtElement(1:ncorners,iel)))
   end do
  

  end subroutine
  
  
  
  !****************************************************************************
  
!<subroutine>  
  
  subroutine releaseTriaData(raddTriaData)

!<description>

  ! Release additional triangulation data.

!</description>

!<output>
  
  ! The additional triangulation data
  type(t_additionalTriaData), intent(inout):: raddTriaData
!</output>
  
!</subroutine>

  ! Deallocate space for local edge numbers
  deallocate(raddTriaData%p_IlocalEdgeNumber)

	! Deallocate space for normal vectors
  deallocate(raddTriaData%p_Dnormals)
    
  ! Deallocate space for edge length
  deallocate(raddTriaData%p_DedgeLength)
  
  ! Deallocate the midpoints
  deallocate(raddTriaData%p_DmidPoints)
  
  ! Deallocate dx and dy
  deallocate(raddTriaData%p_Ddxdy)
  
  end subroutine
  
  
  
  






!****************************************************************************

!<subroutine>

  subroutine linf_dg_buildVectorBlockEdge2d (rform, ccubType, bclear,&
                                              rvectorBlock,&
                                              rvectorBlockSol,&
                                              raddTriaData,&
                                              flux_dg_buildVectorBlEdge2D_sim,&
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
  type(t_vectorBlock), intent(in) :: rvectorBlockSol
  
  ! Additional triangulation data
  type(t_additionalTriaData), intent(in) :: raddTriaData
   
  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection
  
  ! A callback routine for the function to be discretised.
  include 'intf_flux_dg_buildVectorBlEdge2D.inc'
  optional :: flux_dg_buildVectorBlEdge2D_sim
!</input>

!<inputoutput>
  ! The linear form vector. Calculated entries are imposed to this vector.
  type(t_vectorBlock), intent(inout) :: rvectorBlock
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_linfVectorAssembly), dimension(2) :: rvectorAssembly
  type(t_boundary), pointer :: p_rboundary
  type(t_triangulation), pointer :: p_rtriangulation
  integer, dimension(:), pointer :: p_IedgeList
  integer, dimension(:,:), pointer :: p_IlocalEdgeNumber
  integer :: ielementDistr, NMT, NVT, iedge

!  ! If the vector does not exist, stop here.
!  if (rvectorScalar%h_Ddata .eq. ST_NOHANDLE) then  
!    call output_line('Vector not available!',&
!                     OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
!  end if
!
!  ! The vector must be unsorted, otherwise we can not set up the vector.
!  if (rvectorScalar%isortStrategy .gt. 0) then
!    call output_line('Vector must be unsorted!',&
!                     OU_CLASS_ERROR,OU_MODE_STD,'linf_buildVectorScalarBdr2D')
!    call sys_halt()
!  end if
  
  ! Clear the vector if necessary.
  if (bclear) call lsysbl_clearVector (rvectorBlock)

!  ! The vector must provide a discretisation structure
!  if (.not. associated(rvectorScalar%p_rspatialDiscr)) then
!    call output_line('No discretisation associated!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
!    call sys_halt()
!  end if
!
!  ! The discretisation must provide a triangulation structure
!  if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rtriangulation)) then
!    call output_line('No triangulation associated!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
!    call sys_halt()
!  end if
!  
!  ! The discretisation must provide a boundary structure
!  if (.not. associated(rvectorScalar%p_rspatialDiscr%p_rboundary)) then
!    call output_line('No boundary associated!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'linf_dg_buildVectorScalarEdge2d')
!    call sys_halt()
!  end if

  ! Set pointers for quicker access
  p_rboundary => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr%p_rboundary
  p_rtriangulation => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
  ! Do we have a uniform triangulation? Would simplify a lot...
  if (rvectorBlock%RvectorBlock(1)%p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then 
    
    select case(rvectorBlock%cdataType)
      
    case(ST_DOUBLE)

      ! Get number of vertices in the triangulation
      NVT = p_rtriangulation%NVT

      ! Get number of edges belonging to elements
      NMT = p_rtriangulation%NMT
          
      ! Allocate the edgelist
      allocate(p_IedgeList(NMT))
      
      ! All edges
      forall (iedge = 1:NMT) p_IedgeList(iedge)=iedge
     
      ! Initialise the vectorAssembly structures
      call linf_initAssembly(rvectorAssembly(1), rform,&
            rvectorBlock%RvectorBlock(1)%p_rspatialDiscr%RelementDistr(1)%celement,&
            ccubType, LINF_NELEMSIM)
      call linf_initAssembly(rvectorAssembly(2), rform,&
            rvectorBlock%RvectorBlock(1)%p_rspatialDiscr%RelementDistr(1)%celement,&
            ccubType, LINF_NELEMSIM)
            
      ! Assemble the data for all elements in this element distribution
      call linf_dg_assembleSubmeshVectorBlockEdge2d (rvectorAssembly,&
            rvectorBlock, rvectorBlockSol,&
            p_IedgeList(1:NMT),raddTriaData,&
            flux_dg_buildVectorBlEdge2D_sim,&
            rcollection&
             )

          
      ! Release the assembly structure.
      call linf_doneAssembly(rvectorAssembly(1))
      call linf_doneAssembly(rvectorAssembly(2))

      ! Deallocate the edgelist
      deallocate(p_IedgeList)
            
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
 
  end subroutine linf_dg_buildVectorBlockEdge2d















































  !****************************************************************************
  
!<subroutine>  
  
  subroutine linf_dg_assembleSubmeshVectorBlockEdge2d (RvectorAssembly,&
              rvector, rvectorSol,&
              IedgeList, raddTriaData,&
              flux_dg_buildVectorBlEdge2D_sim,&
              rcollection)

!<description>

  ! Assembles the vector entries for a submesh by integration over the given edges.

!</description>

!<input>

  ! List of edges where to assemble the linear form.
  integer, dimension(:), intent(in), target :: IedgeList
    
  ! The solution vector. Used to calculate the solution on the edges.
  type(t_vectorBlock), intent(in) :: rvectorSol
  
  ! Additional triangulation data
  type(t_additionalTriaData), intent(in) :: raddTriaData

  ! A callback routine which is able to calculate the values of the
  ! function $f$ which is to be discretised.
  include 'intf_flux_dg_buildVectorBlEdge2D.inc'
  optional :: flux_dg_buildVectorBlEdge2D_sim 
  
!</input>

!<inputoutput>
  
  ! A vector assembly structure prepared with linf_initAssembly.
  type(t_linfVectorAssembly), dimension(2), intent(inout), target :: rvectorAssembly
  
  ! A vector where to assemble the contributions to.
  type(t_vectorBlock), intent(inout) :: rvector  
  
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
    real(DP) :: domega1,domega2,daux1,daux2,dlen
    real(DP), dimension(:,:), allocatable :: daux
    real(DP) :: dval1, dval2
    integer(I32) :: cevaluationTag
    type(t_linfVectorAssembly), dimension(2), target :: rlocalVectorAssembly
    type(t_domainIntSubset), dimension(2) :: rintSubset
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    integer, dimension(:),pointer :: p_Idescriptors
    integer, dimension(:,:), pointer :: p_Idofs
    type(t_evalElementSet), pointer :: p_revalElementSet

    ! A small vector holding only the additive contributions of
    ! one element
    real(DP), dimension(:,:,:), allocatable :: DlocalData
  
    ! Arrays for cubature points 1D->2D
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    real(DP), dimension(:,:,:,:), allocatable :: Dxi2D,DpointsRef
    
    ! Element list, where to assemble the form
    integer, dimension(:,:), allocatable, target :: IelementList
    
    integer(i32) :: icoordSystem
    
    ! Chooses the element on side ... of the edge
    integer :: iside
    
    ! Number of edges
    integer :: NMT
    
    ! Vertices per element
    integer :: NVE
    
    integer :: ive
    
    ! Pointer to Ielementsatedge in the triangulation
    integer, dimension(:,:), pointer :: p_IelementsAtEdge
    
    ! Pointer to IverticesAtEdge in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    ! Pointer to the vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    
    ! Pointer to IverticesAtEelement in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    
!    ! Array for the solution values in the cubature points
!    real(DP), dimension(:,:,:), allocatable :: DsolVals
         
    ! Array for the normal vectors
    real(DP), dimension(:,:), allocatable :: normal
    
    ! Temp variables for the coordinates of the vertices
    real(DP) :: dxl1, dxl2, dyl1, dyl2
    
    ! Space for the values of the flux function
    real(DP), dimension(:,:,:,:), allocatable :: DfluxValues
    
    ! The Jacobian matrix of the mapping for each point.
    ! DIMENSION(number of entries in the matrix,npointsPerEl,nelements)
    real(DP), dimension(:,:,:), allocatable :: Djac
  
    ! Jacobian determinants of the mapping for all the points from the
    ! reference element to the real element.
    ! DIMENSION(npointsPerEl,nelements)
    real(DP), dimension(:,:), allocatable :: Ddetj
  
    ! Array receiving the coordinates of the points in DpointsRef,
    ! mapped from the reference element to the real element.
    ! If not specified, they are not computed.
    ! DIMENSION(#space dimensions,npointsPerEl,nelements)
    real(DP), dimension(:,:,:), allocatable :: DpointsReal
    
    ! The coordinates of the corner edges of the elements for the transformation
    ! DIMENSION(#space dimensions,NVE,nelements)
    real(DP), dimension(:,:,:), allocatable :: Dcoords
    
    integer :: nvar
    
    integer(I32) :: ctrafotype
    
    ! Pointer to the data of the different blocks of the output vector
    type(t_dpPointer), dimension(:), allocatable :: p_DoutputData
    
    ! Number of elements of the triangulation
    integer :: NEL
    
    ! Number of current variable
    integer :: ivar
  
  
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
    
    ! Get number of variables of the system
    nvar = rvector%nblocks
    
    ! Get number of elements
    NEL = rvectorSol%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
  
    ! Get pointers to elements at edge
    call storage_getbase_int2D(&
          rvectorSol%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%h_IelementsAtEdge,&
          p_IelementsAtEdge)
          
    ! Get pointers to the vertex coordinates
    call storage_getbase_double2D(&
          rvectorSol%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%h_DvertexCoords,&
          p_DvertexCoords)
          
    ! Get pointers to vertices at edge
    call storage_getbase_int2D(&
          rvectorSol%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%h_IverticesAtEdge,&
          p_IverticesAtEdge)
    
    ! Get pointers to vertices at elements
    call storage_getbase_int2D(&
          rvectorSol%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%h_IverticesAtElement,&
          p_IverticesAtElement)   
    
    ! Get the elements adjacent to the given edges
    allocate(IelementList(3,size(IedgeList)))
    IelementList(1:2,1:size(IedgeList))=p_IelementsAtEdge(1:2,IedgeList(:))
    
    ! Allocate daux
    allocate(daux(nvar,2))
    
    ! Copy the second component and replace 0s by 1s
    IelementList(3,size(IedgeList))=IelementList(2,size(IedgeList))
    do iel = 1,size(IedgeList)
      IelementList(3,iel)=max(IelementList(2,iel),1)
    end do
    
    ! Get some pointers for faster access
!    call lsyssc_getbase_double (rvector, p_Ddata)
    indof = rvectorAssembly(1)%indof
    ncubp = rvectorAssembly(1)%ncubp
    NVE = rvectorAssembly(1)%NVE
    
!    ! Allocate space for the solution values in the cubature points
!    allocate(DsolVals(2,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock))
    
!    ! Allocate space for the jacobi matrices in the cubature points
!    allocate(Djac(4,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock))
!    
!    ! Allocate space for the determinants of the jacobi matrices in the cubature points
!    allocate(Ddetj(ncubp,rlocalVectorAssembly(1)%nelementsPerBlock))
!    
!    ! Allocate space for the integration points on the real elements
!    allocate(DpointsReal(ndim2d,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock))

!    ! Allocate space for normal vectors
!    allocate(normal(2,min(size(IedgeList),rlocalVectorAssembly(1)%nelementsPerBlock)))
!    
!    ! Allocate space for edge length
!    allocate(edgelength(min(size(IedgeList),rlocalVectorAssembly(1)%nelementsPerBlock)))
    
!    ! The coordinates of the corner edges of the elements for the transformation
!    allocate(Dcoords(ndim2d,NVE,rlocalVectorAssembly(1)%nelementsPerBlock))
      
  ! Allocate space for the flux variables DIM(nvar,ialbet,ncubp,elementsperblock)
  allocate(DfluxValues(nvar,1,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock))
  
  
   
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
    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock,2))
    
    ! Allocate DlocalData, which will later 
    allocate(DlocalData(nvar,2,indof))
    
    ! Allocate the space for the pointer to the Data of the different blocks of the output vector
    allocate(p_DoutputData(nvar))
    
    !do ivar = 1, nvar
    !  call lsyssc_getbase_double(rvector%RvectorBlock(ivar),p_DoutputData(ivar)%p_Ddata)
    !end do
    call lsysbl_getbase_double (rvector, p_Ddata)

    ! Get the type of coordinate system
    icoordSystem = elem_igetCoordSystem(rlocalVectorAssembly(1)%celement)
    
    ! Get type of transformation
    ctrafotype = elem_igetTrafoType(rlocalVectorAssembly(1)%celement)

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
        call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(1,Iedgelist(IELset+iel-1)), &
            ncubp, Dxi1D, Dxi2D(:,:,1,iel))
        call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(2,Iedgelist(IELset+iel-1)), &
            ncubp, Dxi1D, Dxi2D(:,:,2,iel))
      end do
     
      ! Transpose the coordinate array such that we get coordinates we
      ! can work with.
      do iside = 1,2
        do iel = 1,IELmax-IELset+1
          do icubp = 1,ncubp
            do k = 1,ubound(DpointsRef,1)
              DpointsRef(k,icubp,iel,iside) = Dxi2D(icubp,k,iside,iel)
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
      call dof_locGlobMapping_mult( rvectorSol%RvectorBlock(1)%p_rspatialDiscr, &
          IelementList(1,IELset:IELmax), rlocalVectorAssembly(1)%p_Idofs)
      call dof_locGlobMapping_mult( rvectorSol%RvectorBlock(1)%p_rspatialDiscr, &
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
          cevaluationTag,  rvectorSol%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation, &
          IelementList(1,IELset:IELmax), rlocalVectorAssembly(1)%ctrafoType, &
          DpointsRef=DpointsRef(:,:,:,1))
      call elprep_prepareSetForEvaluation (&
          rlocalVectorAssembly(2)%revalElementSet,&
          cevaluationTag,  rvectorSol%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation, &
          IelementList(3,IELset:IELmax), rlocalVectorAssembly(2)%ctrafoType, &
          DpointsRef=DpointsRef(:,:,:,2))
          
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (rlocalVectorAssembly(1)%celement, &
          rlocalVectorAssembly(1)%revalElementSet,&
          rlocalVectorAssembly(1)%Bder, &
          rlocalVectorAssembly(1)%p_Dbas)
      call elem_generic_sim2 (rlocalVectorAssembly(2)%celement, &
          rlocalVectorAssembly(2)%revalElementSet,&
          rlocalVectorAssembly(2)%Bder, &
          rlocalVectorAssembly(2)%p_Dbas)
          
          
         
!      ! ********** Get solution values in the cubature points *************
!      call lsyssc_getbase_double(rvectorSol,p_DdataSol)
!    
!      ! Now that we have the basis functions, we want to have the function values.
!      ! We get them by multiplying the FE-coefficients with the values of the
!      ! basis functions and summing up.
!      do iel = 1,IELmax-IELset+1      
!        do icubp = 1,ncubp
!          ! Calculate the value in the point
!          dval1 = 0.0_DP
!          dval2 = 0.0_DP
!          do idofe = 1,indof
!            dval1 = dval1 + &
!                   p_DdataSol(rlocalVectorAssembly(1)%p_Idofs(idofe,iel)) &
!                   * rlocalVectorAssembly(1)%p_Dbas(idofe,DER_FUNC,icubp,iel)
!            dval2 = dval2 + &
!                   p_DdataSol(rlocalVectorAssembly(2)%p_Idofs(idofe,iel)) &
!                   * rlocalVectorAssembly(2)%p_Dbas(idofe,DER_FUNC,icubp,iel)
!          end do
!          ! Save the value in the point
!          DsolVals(1,icubp,iel) = dval1
!          DsolVals(2,icubp,iel) = dval2
!        end do
!      end do




      
!     ! Set values at boundary
!     do iel = 1,IELmax-IELset+1
!      if(IelementList(2,IELset+iel-1).eq.0) then
!        DsolVals(2,1:ncubp,iel) = 0.0_DP
!      end if
!      
!      end do
     
      
     ! ---------------------- Get values of the flux function --------------
     
     
!      ! Now it is time to call our coefficient function to calculate the
!      ! function values in the cubature points:
!      if (present(fcoeff_buildVectorScBdr2D_sim)) then
!        call domint_initIntegrationByEvalSet (p_revalElementSet, rintSubset)
!        rintSubset%ielementDistribution = 0
!        rintSubset%ielementStartIdx = IELset
!        rintSubset%p_Ielements => IelementList(IELset:IELmax)
!        rintSubset%p_IdofsTrial => p_Idofs
!        rintSubset%celement = rlocalVectorAssembly%celement
!        call fcoeff_buildVectorScBdr2D_sim (rvector%p_rspatialDiscr,&
!            rlocalVectorAssembly%rform,  IELmax-IELset+1, ncubp,&
!            p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!            ibdc, DpointsPar(:,1:IELmax-IELset+1),&
!            p_Idofs, rintSubset, &
!            p_Dcoefficients(:,:,1:IELmax-IELset+1), rcollection)
!        call domint_doneIntegration (rintSubset)
!      else
!        p_Dcoefficients(:,:,1:IELmax-IELset+1) = 1.0_DP
!      end if





!      ! Fill the corner coordinates of the elements
!      do iel = 1,IELmax-IELset+1
!        do ive = 1, NVE
!          Dcoords(1:ndim2d,ive,iel)=&
!                   p_DvertexCoords(1:ndim2d,p_IverticesAtElement(ive,IelementList(1,IELset+iel-1)))
!        end do
!      end do
!
!
!      ! The numerical flux function needs the x- and y- values
!      ! So we have to call the mapping from the reference- to the real element
!      call trafo_calctrafo_sim (ctrafoType,IELmax-IELset+1,ncubp,Dcoords,&
!                                DpointsRef(1:ndim2d,:,1:IELmax-IELset+1,1),Djac(1:4,1:ncubp,1:IELmax-IELset+1),&
!                                Ddetj(1:ncubp,1:IELmax-IELset+1),&
!                                DpointsReal(1:ndim2d,1:ncubp,1:IELmax-IELset+1))










    ! If the flux function needs other, than just the function values from the solution
    ! (for example the derivatives), we will give an evalElementSet to it
    ! This is filled here
    
    call domint_initIntegrationByEvalSet (rlocalVectorAssembly(1)%revalElementSet, rintSubset(1))
    call domint_initIntegrationByEvalSet (rlocalVectorAssembly(2)%revalElementSet, rintSubset(2))
        !rintSubset(1)%ielementDistribution = 0
        rintSubset(1)%ielementStartIdx = IELset
        rintSubset(1)%p_Ielements => IelementList(1,IELset:IELmax)
        rintSubset(1)%p_IdofsTrial => rlocalVectorAssembly(1)%p_Idofs
        rintSubset(1)%celement = rlocalVectorAssembly(1)%celement
        !rintSubset(2)%ielementDistribution = 0
        rintSubset(2)%ielementStartIdx = IELset
        rintSubset(2)%p_Ielements => IelementList(2,IELset:IELmax)
        rintSubset(2)%p_IdofsTrial => rlocalVectorAssembly(2)%p_Idofs
        rintSubset(2)%celement = rlocalVectorAssembly(2)%celement





      
      
      call flux_dg_buildVectorBlEdge2D_sim (&
!            rlocalVectorAssembly(1)%p_Dcoefficients(1,:,1:IELmax-IELset+1),&
!            DsolVals(:,:,1:IELmax-IELset+1),&
            DfluxValues(:,:,:,1:IELmax-IELset+1),&
            rvectorSol,&
            IelementList(2,IELset:IELmax),&
            raddTriaData%p_Dnormals(:,Iedgelist(IELset:IELmax)),&
            !DpointsReal(1:ndim2d,1:ncubp,1:IELmax-IELset+1),&
            rintSubset,&
            rcollection )
            
      
      call domint_doneIntegration (rintSubset(1))
      call domint_doneIntegration (rintSubset(2))
      
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
        DlocalData(1:nvar,1:2,1:indof) = 0.0_DP

        ! Get the length of the edge. Let us use the parameter values
        ! on the boundary for that purpose; this is a more general
        ! implementation than using simple lines as it will later 
        ! support isoparametric elements.
        !
        ! The length of the current edge serves as a "determinant"
        ! in the cubature, so we have to divide it by 2 as an edge on 
        ! the unit interval [-1,1] has length 2.
        dlen = 0.5_DP*raddTriaData%p_Dedgelength(Iedgelist(IELset+iel-1))

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature
          ! formula in that cubature point.

          domega1 = dlen * rlocalVectorAssembly(1)%p_Domega(icubp)
          domega2 = dlen * rlocalVectorAssembly(1)%p_Domega(icubp)

          
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
            !daux1 = domega1 * rlocalVectorAssembly(1)%p_Dcoefficients(ialbet,icubp        ,iel)
            !daux2 = domega2 * rlocalVectorAssembly(1)%p_Dcoefficients(ialbet,ncubp-icubp+1,iel) *(-1.0_dp)
            Daux(:,1) = domega1 * DfluxValues(:,ialbet,icubp        ,iel)
            Daux(:,2) = domega2 * DfluxValues(:,ialbet,ncubp-icubp+1,iel) *(-1.0_dp)
            
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
              
              DlocalData(:,1,idofe) = DlocalData(:,1,idofe)+&
                                    rlocalVectorAssembly(1)%p_Dbas(idofe,ia,icubp,iel)*daux(:,1)
              
!              if(IelementList(2,IELset+iel-1).ne.0) then
!                DlocalData(2,idofe) = DlocalData(2,idofe)+&
!                                      rlocalVectorAssembly(2)%p_Dbas(idofe,ia,icubp,iel)*daux2
!              end if
                    
              DlocalData(:,2,idofe) = DlocalData(:,2,idofe)+&
                                    rlocalVectorAssembly(2)%p_Dbas(idofe,ia,icubp,iel)*daux(:,2)*&
                                    real(min(1,IelementList(2,IELset+iel-1)))
              
             
            end do ! idofe
            
          end do ! ialbet

        end do ! icubp 
        
        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        
        do ivar = 1, nvar
        
          do idofe = 1,indof
                    
            p_Ddata(rvector%RvectorBlock(ivar)%iidxFirstEntry+rlocalVectorAssembly(1)%p_Idofs(idofe,iel)-1) = &            
              p_Ddata(rvector%RvectorBlock(ivar)%iidxFirstEntry+rlocalVectorAssembly(1)%p_Idofs(idofe,iel)-1) + &
              DlocalData(ivar,1,idofe)
            
            p_Ddata(rvector%RvectorBlock(ivar)%iidxFirstEntry+rlocalVectorAssembly(2)%p_Idofs(idofe,iel)-1) = &            
              p_Ddata(rvector%RvectorBlock(ivar)%iidxFirstEntry+rlocalVectorAssembly(2)%p_Idofs(idofe,iel)-1) + &
              DlocalData(ivar,2,idofe)
        
!            p_DoutputData(ivar)%p_Ddata(rlocalVectorAssembly(1)%p_Idofs(idofe,iel)) =&
!                         p_DoutputData(ivar)%p_Ddata(rlocalVectorAssembly(1)%p_Idofs(idofe,iel)) +&
!                         DlocalData(ivar,1,idofe)
!            p_DoutputData(ivar)%p_Ddata(rlocalVectorAssembly(2)%p_Idofs(idofe,iel)) =&
!                         p_DoutputData(ivar)%p_Ddata(rlocalVectorAssembly(2)%p_Idofs(idofe,iel)) +&
!                         DlocalData(ivar,2,idofe)
          end do
          
        end do ! nvar
        
      end do ! iel

    end do ! IELset



    
    ! Release the local vector assembly structure
    call linf_releaseAssemblyData(rlocalVectorAssembly(1))
    call linf_releaseAssemblyData(rlocalVectorAssembly(2))

    ! Deallocate memory
    deallocate(Dxi2D,DpointsRef,IelementList)!,DsolVals,edgelength,normal,Djac,Ddetj,DpointsReal,Dcoords)
    deallocate(DfluxValues,daux,DlocalData,p_DoutputData)

  end subroutine
  
  
  
  
  
  
  
  
  
  
  !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_linearLimiterBlockIndicatorVar (rvectorBlock, iindicatorVar)

!<description>

  ! Limits the linear part of a dg_T1 element vector.

!</description>

!<input>
!</input>

!<inputoutput>

  ! A vector to limit
  type(t_vectorBlock), intent(inout) :: rvectorBlock
  
  ! Which component of the blockvector is the indicatorvariable
  integer , intent(in) :: iindicatorVar
  
!</inputoutput>
  
!</subroutine>

  ! local variables, used by all processors
  real(DP), dimension(:), pointer :: p_Ddata, p_DdataOut
  integer :: indof, NEL, iel, NVE, ivt, NVT
    
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! The coordinates of the points in which to evaluate the solution vector
  real(dp), dimension(2,17) :: Dpoints
  
  ! The list of elements, in which these points can be found
  integer, dimension(17) :: Ielements
  
  ! The values of the solution vector in the points
  real(dp), dimension(17) :: Dvalues
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  real(dp), dimension(:), allocatable :: duimax, duimin
  
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  
  real(dp) :: xc,yc
  integer :: iidx, nvert, ivert, ineighbour, ineighElm
  integer, dimension(4) :: IhomeIndex
  
  real(dp) :: dui, ddu, dalpha, dalphatemp, duc
  
  integer, dimension(3) :: IdofGlob
  
  integer, dimension(5) :: Isep
  
  integer :: NVBD
  
  integer :: nvar, ivar
  
  ! Pointer to the scalar vectorblock, which is the indicator variable
  ! A vector to limit
  type(t_vectorScalar), pointer :: p_rvector
  
  ! Array of pointers to the data of the blockvector to limit
  type(t_dpPointer), dimension(:), allocatable :: p_DoutputData
  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
  
  ! Allocate the space for the pointer to the Data of the different blocks of the output vector
  allocate(p_DoutputData(nvar))
    
!  do ivar = 1, nvar
!    call lsyssc_getbase_double(rvectorBlock%RvectorBlock(ivar),p_DoutputData(ivar)%p_Ddata)
!  end do
  
  call lsysbl_getbase_double (rvectorBlock, p_DdataOut)
  
  ! Point to the indicatorvariable
  p_rvector => rvectorBlock%RvectorBlock(iindicatorVar)

  ! Get pointer to the solution data
  call lsyssc_getbase_double (p_rvector,p_Ddata)
  
  ! Get pointers for quicker access
  p_rspatialDiscr => p_rvector%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Get pointer to the data of the vector
  
  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  ! Number of vertives at boundary
  NVBD = p_rtriangulation%NVBD
  
  ! Get pointers to the data form the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
                               p_IelementsAtVertex) 
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
                               p_IelementsAtVertexIdx)
  call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary ,&
                               p_IverticesAtBoundary)                               
                               
                               
   ! Set pointer to coordinate vector
   call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

   ! Set pointer to vertices at element
   call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
        
  NVT = p_rtriangulation%NVT
  
  allocate(duimax(NVT),duimin(NVT))
  
  duimax= -SYS_MAXREAL
  duimin=  SYS_MAXREAL
  
  do iel = 1, NEL
    
    ! Get global DOFs of the element
    call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    duc = p_Ddata(IdofGlob(1))
    
    ! elem_igetNVE(celement)
    NVE = 4
    
    do ivt = 1, NVE
      nvt = p_IverticesAtElement(ivt,iel)
      duimax(nvt) = max(duc,duimax(nvt))
      duimin(nvt) = min(duc,duimin(nvt))
    end do
    
    
!    NVE = 4
!    
!    do ivt = 1, NVE
!      nvt = p_IverticesAtElement(ivt,iel)
!      Dpoints(1,ivt) = p_DvertexCoords(1,nvt)
!      Dpoints(2,ivt) = p_DvertexCoords(2,nvt)
!      Ielements(ivt) = iel
!    end do
!    
!    call fevl_evaluate (DER_FUNC, Dvalues(1:4), p_rvector, Dpoints(1:2,1:4), &
!                          Ielements(1:4))
!   do ivt = 1, NVE
!      nvt = p_IverticesAtElement(ivt,iel)
!      duimax(nvt) = max(duimax(nvt),Dvalues(ivt))
!      duimin(nvt) = min(duimin(nvt),Dvalues(ivt))
!    end do
    
    
  end do
  
!  do ivt = 1, NVBD
!    duimax(p_IverticesAtBoundary(ivt)) = 100000000.0_DP
!    duimin(p_IverticesAtBoundary(ivt)) = -100000000.0_DP
!  end do

                               
  do iel = 1, NEL
  
    ! Get number of corner vertices
    ! elem_igetNVE(celement)
    NVE = 4
    
    ! Get midpoint of the element
    xc = &
         (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(4,iel)))/4.0_dp

    yc = &
         (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(4,iel)))/4.0_dp
          
    ! The first point we want to evaluate the solution in, is the midpoint of the element
    Dpoints(1,1) = xc
    Dpoints(2,1) = yc
    Ielements(1) = iel
    
    ! Initialise the limiting factor
    dalpha = 1.0_dp
    
    ! Now start to set the points, where to evaluate the solution
    
    ! Loop over the vertices of the element
    do ivert = 1, NVE
		
      nvert = p_IverticesAtElement(ivert, iel)
      
      ! The second point we want to evaluate the solution in, is in the corner of the mother element
      Dpoints(1,1+ivert) = p_DvertexCoords(1,nvert)
      Dpoints(2,1+ivert) = p_DvertexCoords(2,nvert)
      Ielements(1+ivert) = iel
	  end do
	  
    ! Evaluate the solution
    call fevl_evaluate (DER_FUNC, Dvalues(1:5), p_rvector, Dpoints(1:2,1:5), &
                          Ielements(1:5))
      
     ! Start calculating the limiting factor
     duc = Dvalues(1)
      
    do ivert = 1, NVE  
      dui = Dvalues(1+ivert)
      ddu = dui-duc
      nvert = p_IverticesAtElement(ivert, iel)
            
      ! Find the maximum/minimum value of the solution in the centroids
      ! of all elements containing this vertex
      if (ddu > 0.0_dp) then
        dalphatemp = min(1.0_dp, (duimax(nvert)-duc)/ddu)
      elseif (ddu < 0.0_dp) then
        dalphatemp = min(1.0_dp, (duimin(nvert)-duc)/ddu)
      else ! (dui==duc)
        dalphatemp = 1.0_dp
      end if
      
      dalpha = min(dalphatemp,dalpha)
      

    end do ! ivert
    
    
    ! Now we have the limitingfactor dalpha
    ! We need to multiply it with the corresponding DOFs
    
    ! Get global DOFs of the element
    call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    ! Multiply the linear part of the solution vector with the correction factor
    
    do ivar = 1, nvar
      !p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3))*dalpha
      p_DdataOut(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3)-1) = &            
        p_DdataOut(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3)-1) * &
        dalpha
    end do
    
  end do ! iel
  
  deallocate(duimax,duimin)
  deallocate(p_DoutputData)

  
  end subroutine
  
  
  
  
  !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_quadraticLimiterBlockIndicatorVar (rvectorBlock, iindicatorVar)

!<description>

  ! Limits the dg_T2 element vector.

!</description>

!<input>
!</input>

!<inputoutput>

  ! A vector to limit
  type(t_vectorBlock), intent(inout) :: rvectorBlock
  
  ! Which component of the blockvector is the indicatorvariable
  integer , intent(in) :: iindicatorVar
  
!</inputoutput>
  
!</subroutine>

  ! local variables, used by all processors
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: indof, NEL, iel, NVE, ivt, NVT
    
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! The coordinates of the points in which to evaluate the solution vector
  real(dp), dimension(2,17) :: Dpoints
  
  ! The list of elements, in which these points can be found
  integer, dimension(17) :: Ielements
  
  ! The values of the solution vector in the points
  real(dp), dimension(17) :: Dvalues
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  real(dp), dimension(:), allocatable :: duimax, duimin
  
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  
  real(dp) :: xc,yc
  integer :: iidx, nvert, ivert, ineighbour, ineighElm
  integer, dimension(4) :: IhomeIndex
  
  real(dp) :: dui, ddu, dalphatemp, duc
  
  integer, dimension(6) :: IdofGlob
  
  integer, dimension(5) :: Isep
  
  integer :: NVBD
  
  integer :: nvar, ivar
  
  integer :: ilim, ideriv
  
  real(DP), dimension(:,:), allocatable :: Dalpha
  
  ! Pointer to the scalar vectorblock, which is the indicator variable
  ! A vector to limit
  type(t_vectorScalar), pointer :: p_rvector
  
  ! Array of pointers to the data of the blockvector to limit
  type(t_dpPointer), dimension(:), allocatable :: p_DoutputData
  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
  
  ! Allocate the space for the pointer to the Data of the different blocks of the output vector
  allocate(p_DoutputData(nvar))
    
  do ivar = 1, nvar
    call lsyssc_getbase_double(rvectorBlock%RvectorBlock(ivar),p_DoutputData(ivar)%p_Ddata)
  end do
  
  ! Point to the indicatorvariable
  p_rvector => rvectorBlock%RvectorBlock(iindicatorVar)

  ! Get pointer to the solution data
  call lsyssc_getbase_double (p_rvector,p_Ddata)
  
  ! Get pointers for quicker access
  p_rspatialDiscr => p_rvector%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Get pointer to the data of the vector
  
  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  ! Number of vertives at boundary
  NVBD = p_rtriangulation%NVBD
  
  ! Get pointers to the data from the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
                               p_IelementsAtVertex) 
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
                               p_IelementsAtVertexIdx)
  call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary ,&
                               p_IverticesAtBoundary)                               
                               
                               
   ! Set pointer to coordinate vector
   call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

   ! Set pointer to vertices at element
   call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
        
  NVT = p_rtriangulation%NVT
  
  allocate(Duimax(NVT),Duimin(NVT),Dalpha(3,NEL))
  
  Dalpha = 1.0_dp
  
  ! Now limit the x- and y- derivative and finally the linear part
  
  do ilim = 1, 3
  
    
  
    select case (ilim)
      case (1)
        ideriv = DER_DERIV_X
      case (2)
        ideriv = DER_DERIV_Y
      case (3)
        ideriv = DER_FUNC
    end select
  
  
  duimax= -SYS_MAXREAL
  duimin=  SYS_MAXREAL
  
  do iel = 1, NEL
    
    ! Get number of corner vertices
    ! elem_igetNVE(celement)
    NVE = 4
    
    ! Get midpoint of the element
    xc = &
         (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(4,iel)))/4.0_dp

    yc = &
         (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(4,iel)))/4.0_dp
          
    ! The first point we want to evaluate the solution in, is the midpoint of the element
    Dpoints(1,1) = xc
    Dpoints(2,1) = yc
    Ielements(1) = iel
    
    ! Evaluate the solution
    call fevl_evaluate (ideriv, Dvalues(1:1), p_rvector, Dpoints(1:2,1:1), &
                          Ielements(1:1))
                          
    duc = Dvalues(1)
    
    do ivt = 1, NVE
      nvt = p_IverticesAtElement(ivt,iel)
      duimax(nvt) = max(duc,duimax(nvt))
      duimin(nvt) = min(duc,duimin(nvt))
    end do
    
    
    
  end do
  

                               
  do iel = 1, NEL
  
    ! Get number of corner vertices
    ! elem_igetNVE(celement)
    NVE = 4
    
    ! Get midpoint of the element
    xc = &
         (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(4,iel)))/4.0_dp

    yc = &
         (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(4,iel)))/4.0_dp
          
    ! The first point we want to evaluate the solution in, is the midpoint of the element
    Dpoints(1,1) = xc
    Dpoints(2,1) = yc
    Ielements(1) = iel
   
    ! Now start to set the points, where to evaluate the solution
    
    ! Loop over the vertices of the element
    do ivert = 1, NVE
		
      nvert = p_IverticesAtElement(ivert, iel)
      
      ! The second point we want to evaluate the solution in, is in the corner of the mother element
      Dpoints(1,1+ivert) = p_DvertexCoords(1,nvert)
      Dpoints(2,1+ivert) = p_DvertexCoords(2,nvert)
      Ielements(1+ivert) = iel
	  end do
	  
    ! Evaluate the solution
    call fevl_evaluate (ideriv, Dvalues(1:5), p_rvector, Dpoints(1:2,1:5), &
                          Ielements(1:5))
      
     ! Start calculating the limiting factor
     duc = Dvalues(1)
      
    do ivert = 1, NVE  
      dui = Dvalues(1+ivert)
      ddu = dui-duc
      nvert = p_IverticesAtElement(ivert, iel)
            
      ! Find the maximum/minimum value of the solution in the centroids
      ! of all elements containing this vertex
      if (ddu > 0.0_dp) then
        dalphatemp = min(1.0_dp, (duimax(nvert)-duc)/ddu)
      elseif (ddu < 0.0_dp) then
        dalphatemp = min(1.0_dp, (duimin(nvert)-duc)/ddu)
      else ! (dui==duc)
        dalphatemp = 1.0_dp
      end if
      
      Dalpha(ilim,iel) = min(dalphatemp,Dalpha(ilim,iel))
      

    end do ! ivert
        
  end do ! iel
  
  
  
  
  select case (ilim)
    case (2)
  
      ! Now we have the limitingfactors for the quadratic part in Dalpha(1:2,:)
      ! Get the Minimum and multiply it with the corresponding DOFs
      do iel = 1, NEL  
  
        Dalpha(1,iel) = min(Dalpha(1,iel),Dalpha(2,iel))
      
        ! Get global DOFs of the element
        call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    
        ! Multiply the quadratic part of the solution vector with the correction factor
        do ivar = 1, nvar
          p_DoutputData(ivar)%p_Ddata(IdofGlob(4:6)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(4:6))*Dalpha(1,iel)
        end do
    
      end do ! iel
    
    case (3)
      do iel = 1, NEL  
  
        Dalpha(3,iel) = max(Dalpha(1,iel),Dalpha(3,iel))
     
        ! Get global DOFs of the element
        call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    
        ! Multiply the linear part of the solution vector with the correction factor
        do ivar = 1, nvar
          p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3))*Dalpha(3,iel)
        end do
    
      end do ! iel
      
    end select
  
  end do ! ilim
  
  
  
  deallocate(duimax,duimin,dalpha)
  

  deallocate(p_DoutputData)

  
  end subroutine
  
  
  
  
  
  
  
  !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_linearLimiterBlockCharVar (rvectorBlock)

!<description>

  ! Limits the linear part of a dg_T1 element vector.

!</description>

!<input>
!</input>

!<inputoutput>

  ! A vector to limit
  type(t_vectorBlock), intent(inout) :: rvectorBlock
  
!</inputoutput>
  
!</subroutine>

  ! local variables, used by all processors
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: indof, NEL, iel, NVE, ivt, NVT
    
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! The coordinates of the points in which to evaluate the solution vector
  real(dp), dimension(2,17) :: Dpoints
  
  ! The list of elements, in which these points can be found
  integer, dimension(17) :: Ielements
  
  ! The values of the solution vector in the points
  real(dp), dimension(17) :: Dvalues
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  real(dp), dimension(:), allocatable :: duimax, duimin
  
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  
  real(dp) :: xc,yc
  integer :: iidx, nvert, ivert, ineighbour, ineighElm
  integer, dimension(4) :: IhomeIndex
  
  real(dp) :: dui, ddu, duc
  
  integer, dimension(3) :: IdofGlob
  
  integer, dimension(5) :: Isep
  
  integer :: NVBD
  
  integer :: nvar, ivar, idim
  
  real(dp), dimension(:), allocatable :: DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi
  
  real(dp), dimension(:,:), allocatable :: DLin, DtLin, Dalphaei, DL, DR, Dalpha
    
  ! Array of pointers to the data of the blockvector to limit
  type(t_dpPointer), dimension(:), allocatable :: p_DoutputData
  
  integer :: iglobVtNumber, iglobNeighNum
  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Allocate the space for the pointer to the Data of the different blocks of the output vector
  allocate(p_DoutputData(nvar))
    
!  do ivar = 1, nvar
!    call lsyssc_getbase_double(rvectorBlock%RvectorBlock(ivar),p_DoutputData(ivar)%p_Ddata)
!  end do
  call lsysbl_getbase_double (rvectorBlock, p_Ddata)
  
  
    
  
  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  ! Number of vertices per element
  !NVE = p_rtriangulation%NVE
  ! elem_igetNVE(celement)
  NVE = 4
  
  ! Number of vertives at boundary
  NVBD = p_rtriangulation%NVBD
  
  ! Number of vertices
  NVT = p_rtriangulation%NVT
  
  ! Get pointers to the data form the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
                               p_IelementsAtVertex) 
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
                               p_IelementsAtVertexIdx)
  call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary ,&
                               p_IverticesAtBoundary)                               
                              
  ! Allocate the space for solution differences, transformed solution differences,
  ! limited transformed solution differences, limited backtransformed solution differences
  ! and limiting factors
  allocate(DVec(nvar), DVei(nvar), DIi(nvar), DtIi(nvar), DtLinMax(nvar), DtLinMin(nvar), DltIi(nvar), DlIi(nvar))
  !allocate(DLin(nvar,NVE-1), DtLin(nvar,NVE-1), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
  allocate(DLin(nvar,10), DtLin(nvar,10), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
  
  do iel = 1, NEL
    
    ! Get global DOFs of the element
    call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    ! Get values in the center of the element for all variables
    do ivar = 1, nvar
      !DVec(ivar) = p_DoutputData(ivar)%p_Ddata(IdofGlob(1))
      DVec(ivar) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1)-1)
    end do
    
    ! Here we should maybe get a local value of NVE
    
    ! Initialise the correction factor
    !Dalphaei(:,:) = 1.0_dp
    Dalphaei(:,:) = 0.0_dp
    
    ! Now calculate the limiting factor for every vertex on our element
    do ivt = 1, NVE
      
      ! Get global vertex number of our local vertex
      iglobVtNumber = p_IverticesAtElement(ivt,iel)
      
      ! Get solution value in this edge in the element
      Ielements(1) = iel
      
      do ivar = 1, nvar
        call fevl_evaluate (DER_FUNC, DVei(ivar:ivar), rvectorBlock%RvectorBlock(ivar),&
                            p_DvertexCoords(1:2,iglobVtNumber:iglobVtNumber), Ielements(1:1))
      end do
      
      ! Calculate solution difference
      DIi = DVei - DVec
      !write(*,*) DIi
      
      ! Get center values of all variables in all neighbour vertices and calculate
      ! the solution differences Dlin(nvar,nneighbors)
      iidx = 0
      DLin = 0.0_dp
      do ineighbour = p_IelementsAtVertexIdx(iglobVtNumber), p_IelementsAtVertexIdx(iglobVtNumber+1)-1
        iglobNeighNum = p_IelementsAtVertex(ineighbour)
        !if (iglobNeighNum.ne.iel) then
          call dof_locGlobMapping(p_rspatialDiscr, iglobNeighNum, IdofGlob)
          iidx = iidx + 1          
          do ivar = 1, nvar
            !DLin(ivar,iidx) = p_DoutputData(ivar)%p_Ddata(IdofGlob(1)) - DVec(ivar)
            DLin(ivar,iidx) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1)-1) - DVec(ivar)
          end do

          
        !end if
      end do
          
      if (iidx>0) then
      ! Dimensional splitting
      do idim = 1, NDIM2D
        ! Now we need the trafo matrices
        DL = buildInvTrafo(DVec,idim)
        DR = buildTrafo(DVec,idim)
        
        ! Transform the solution differences
        DtIi = matmul(DL,DIi)
        do ineighbour = 1, iidx
          DtLin(:,ineighbour) = matmul(DL,Dlin(:,ineighbour))
        end do
        

        
        ! Get max and min of the transformed solution differences
        do ivar = 1, nvar
          DtLinMax(ivar) = maxval(DtLin(ivar,1:iidx))
          DtLinMin(ivar) = minval(DtLin(ivar,1:iidx))
          !DtLinMax(ivar) = max(DtLinMax(ivar),0.0_dp)
          !DtLinMin(ivar) = min(DtLinMin(ivar),0.0_dp)
        end do
       
        ! Now, as the differences are transformed, we can limit every component on its own
        do ivar = 1, nvar
          DltIi(ivar) = max(min(DtLinMax(ivar),DtIi(ivar)),DtLinMin(ivar))
        end do
        
        
        
        ! Now we can trafo back
        DlIi = matmul(DR,DltIi)
        
        ! Calculate the correction factor
        ! for this element, for this edge, for this dimension (take min of all dimensions)
        do ivar = 1, nvar
          if (abs(DIi(ivar))<SYS_EPSREAL) then
            !Dalphaei(ivar,ivt) = 1.0_dp
          else
            ! This is the one following the principles
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))
            
            ! This one is less limiting
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), min(abs(DlIi(ivar)/DIi(ivar)),1.0_dp ))
            
            ! This is the one with the max of the dimensional splitting
            !Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))
            
            ! This is the least limiting one
            Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt), min(abs(DlIi(ivar)/DIi(ivar)),1.0_dp ) )
          end if
        end do
      
      end do ! idim
      end if

      
    end do ! ivt


    
    ! Get minimum of all correction factors of all vertices on this element
    do ivar = 1, nvar
      Dalpha(ivar, iel) = minval(Dalphaei(ivar,:))
    end do
    
    
!    if (iel.eq.985) then
!      write(*,*) Dalphaei(1,:)
!      write(*,*) Dalpha(1, iel)
!    end if
    
    ! Get global DOFs of the element
    call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    ! Multiply the linear part of the solution vector with the correction factor
    do ivar = 1, nvar
      ! p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3))*Dalpha(ivar, iel)
      
      p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3)-1) = &            
        p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3)-1) * &
        Dalpha(ivar, iel)
      
    end do
    
  end do !iel

  
  deallocate(p_DoutputData)
  deallocate(DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi)
  deallocate(DLin, DtLin, Dalphaei)

  
  end subroutine
 
 
 
 
 
 
  !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_quadraticLimiterBlockIndicatorVar_2 (rvectorBlock, iindicatorVar)

!<description>

  ! Limits the dg_T2 element vector.

!</description>

!<input>
!</input>

!<inputoutput>

  ! A vector to limit
  type(t_vectorBlock), intent(inout) :: rvectorBlock
  
  ! Which component of the blockvector is the indicatorvariable
  integer , intent(in) :: iindicatorVar
  
!</inputoutput>
  
!</subroutine>

  ! local variables, used by all processors
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: indof, NEL, iel, NVE, ivt, NVT
    
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! The coordinates of the points in which to evaluate the solution vector
  real(dp), dimension(2,17) :: Dpoints
  
  ! The list of elements, in which these points can be found
  integer, dimension(17) :: Ielements
  
  ! The values of the solution vector in the points
  real(dp), dimension(17) :: Dvalues
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  real(dp), dimension(:), allocatable :: duimax, duimin
  
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  
  real(dp) :: xc,yc
  integer :: iidx, nvert, ivert, ineighbour, ineighElm
  integer, dimension(4) :: IhomeIndex
  
  real(dp) :: dui, ddu, dalphatemp, duc
  
  integer, dimension(6) :: IdofGlob
  
  integer, dimension(5) :: Isep
  
  integer :: NVBD
  
  integer :: iglobNeighNum, iglobVtNumber
  
  integer :: nvar, ivar
  
  integer :: ilim, ideriv
  
  real(DP), dimension(:,:), allocatable :: Dalpha
  
  real(dp), dimension(10) :: Dalphaei
  
  real(dp) :: ddiff, dalphae, dimax, dimin
  
  ! Pointer to the scalar vectorblock, which is the indicator variable
  ! A vector to limit
  type(t_vectorScalar), pointer :: p_rvector
  
  ! Array of pointers to the data of the blockvector to limit
  type(t_dpPointer), dimension(:), allocatable :: p_DoutputData
  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
  
  ! Allocate the space for the pointer to the Data of the different blocks of the output vector
  allocate(p_DoutputData(nvar))
    
  do ivar = 1, nvar
    call lsyssc_getbase_double(rvectorBlock%RvectorBlock(ivar),p_DoutputData(ivar)%p_Ddata)
  end do
  
  ! Point to the indicatorvariable
  p_rvector => rvectorBlock%RvectorBlock(iindicatorVar)

  ! Get pointer to the solution data
  call lsyssc_getbase_double (p_rvector,p_Ddata)
  
  ! Get pointers for quicker access
  p_rspatialDiscr => p_rvector%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Get pointer to the data of the vector
  
  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  ! Number of vertives at boundary
  NVBD = p_rtriangulation%NVBD
  
  ! Get pointers to the data from the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
                               p_IelementsAtVertex) 
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
                               p_IelementsAtVertexIdx)
  call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary ,&
                               p_IverticesAtBoundary)                               
                               
                               
   ! Set pointer to coordinate vector
   call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

   ! Set pointer to vertices at element
   call storage_getbase_int2D(&
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
        
  NVT = p_rtriangulation%NVT
  
  !allocate(Duimax(NVT),Duimin(NVT))
  allocate(Dalpha(3,NEL))
  
  Dalpha = 1.0_dp
  
  ! Now limit the x- and y- derivative and finally the linear part
  
  do ilim = 1, 3
  
    
  
    select case (ilim)
      case (1)
        ideriv = DER_DERIV_X
      case (2)
        ideriv = DER_DERIV_Y
      case (3)
        ideriv = DER_FUNC
    end select
  
  
!  duimax= -SYS_MAXREAL
!  duimin=  SYS_MAXREAL
!  
!  do iel = 1, NEL
!    
!    ! Get number of corner vertices
!    ! elem_igetNVE(celement)
!    NVE = 4
!    
!    ! Get midpoint of the element
!    xc = &
!         (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
!          p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
!          p_DvertexCoords(1,p_IverticesAtElement(3,iel))+&
!          p_DvertexCoords(1,p_IverticesAtElement(4,iel)))/4.0_dp
!
!    yc = &
!         (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
!          p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
!          p_DvertexCoords(2,p_IverticesAtElement(3,iel))+&
!          p_DvertexCoords(2,p_IverticesAtElement(4,iel)))/4.0_dp
!          
!    ! The first point we want to evaluate the solution in, is the midpoint of the element
!    Dpoints(1,1) = xc
!    Dpoints(2,1) = yc
!    Ielements(1) = iel
!    
!    ! Evaluate the solution
!    call fevl_evaluate (ideriv, Dvalues(1:1), p_rvector, Dpoints(1:2,1:1), &
!                          Ielements(1:1))
!                          
!    duc = Dvalues(1)
!    
!    do ivt = 1, NVE
!      nvt = p_IverticesAtElement(ivt,iel)
!      duimax(nvt) = max(duc,duimax(nvt))
!      duimin(nvt) = min(duc,duimin(nvt))
!    end do
!    
!    
!    
!  end do
  

                               
  do iel = 1, NEL
  
    ! Get number of corner vertices
    ! elem_igetNVE(celement)
    NVE = 4
    
    ! Get midpoint of the element
    xc = &
         (p_DvertexCoords(1,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(1,p_IverticesAtElement(4,iel)))/4.0_dp

    yc = &
         (p_DvertexCoords(2,p_IverticesAtElement(1,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(2,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(3,iel))+&
          p_DvertexCoords(2,p_IverticesAtElement(4,iel)))/4.0_dp
          
    ! The first point we want to evaluate the solution in, is the midpoint of the element
    Dpoints(1,1) = xc
    Dpoints(2,1) = yc
    Ielements(1) = iel
   
    ! Now start to set the points, where to evaluate the solution
    
    call fevl_evaluate (ideriv, Dvalues(1:1), p_rvector, Dpoints(1:2,1:1), &
                        Ielements(1:1))
    
    ! Center (derivative) value
    duc = Dvalues(1)
    
          
    ! Initialise the correction factor
    Dalphaei(:) = 1.0_dp
    
    ! Now calculate the limiting factor for every vertex on our element
    do ivt = 1, NVE
      
      ! Get global vertex number of our local vertex
      iglobVtNumber = p_IverticesAtElement(ivt,iel)
      
      ! Get solution value in this edge in the element
      Dpoints(1,2) = p_DvertexCoords(1,iglobVtNumber)
      Dpoints(2,2) = p_DvertexCoords(2,iglobVtNumber)
      Ielements(2) = iel
      
      
!      call fevl_evaluate (DER_FUNC, Dvalues(1:1), p_rvector,&
!                          p_DvertexCoords(1:2,iglobVtNumber:iglobVtNumber), Ielements(1:1))
!
!      ! Solution value in this edge                    
!      dui = Dvalues(1:1)
      
      
      ! Calculate solution difference
      !DIi = DVei - DVec
      !write(*,*) DIi
      
      ! Get center values of all variables in all neighbour vertices and calculate
      ! the solution differences Dlin(nvar,nneighbors)
      iidx = 2
      
      do ineighbour = p_IelementsAtVertexIdx(iglobVtNumber), p_IelementsAtVertexIdx(iglobVtNumber+1)-1
        iglobNeighNum = p_IelementsAtVertex(ineighbour)
        if (iglobNeighNum.ne.iel) then
          iidx = iidx + 1
          Dpoints(1,iidx) = p_DvertexCoords(1,iglobVtNumber)
          Dpoints(2,iidx) = p_DvertexCoords(2,iglobVtNumber)
          Ielements(iidx) = iglobNeighNum
        end if
      end do
      
      ! Initialise the factor for this edge
      dalphae = 1.0_dp
      
      if (iidx.ne.2) then
      
      ! Evaluate the (derivatives) of the solution
      call fevl_evaluate (ideriv, Dvalues(1:iidx), p_rvector,&
                          Dpoints(1:2,1:iidx), Ielements(1:iidx))
      
      ! Center (derivative) value
      duc = Dvalues(1)
      ! Solution value in this edge                    
      dui = Dvalues(2)
      
      ! Now calculate limiting factor
      ddiff = dui-duc
      
      dimax = maxval(Dvalues(3:iidx))
      dimin = minval(Dvalues(3:iidx))
      
      ! Loop over all elements adjacent to the edge
      do ineighbour = 3, iidx
        if (ddiff >0) then
          dalphae = min(dalphae,min(1.0_dp,(dimax-duc)/ddiff))
        elseif (ddiff <0) then
          dalphae = min(dalphae,min(1.0_dp,(dimin-duc)/ddiff))
        end if
      end do
      
      end if
      
      Dalphaei (ivt) = dalphae
        
    end do ! ivt
    
    ! Get minimum of all correction factors of all vertices on this element
    Dalpha(ilim,iel) = max(0.0_dp,minval(Dalphaei(1:ivt-1)))
    
  end do ! iel
    
  
  select case (ilim)
    case (2)
  
      ! Now we have the limitingfactors for the quadratic part in Dalpha(1:2,:)
      ! Get the Minimum and multiply it with the corresponding DOFs
      do iel = 1, NEL  
  
        Dalpha(1,iel) = min(Dalpha(1,iel),Dalpha(2,iel))
      
        ! Get global DOFs of the element
        call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    
        ! Multiply the quadratic part of the solution vector with the correction factor
        do ivar = 1, nvar
          p_DoutputData(ivar)%p_Ddata(IdofGlob(4:6)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(4:6))*Dalpha(1,iel)
        end do
    
      end do ! iel
    
    case (3)
      do iel = 1, NEL  
  
        Dalpha(3,iel) = max(Dalpha(1,iel),Dalpha(3,iel))
     
        ! Get global DOFs of the element
        call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    
        ! Multiply the linear part of the solution vector with the correction factor
        do ivar = 1, nvar
          p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3))*Dalpha(3,iel)
        end do
    
      end do ! iel
      
    end select
  
  end do ! ilim
  
  
  
  !deallocate(duimax,duimin)
  deallocate(dalpha)
  

  deallocate(p_DoutputData)

  
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_quadraticLimiterBlockCharVar (rvectorBlock, raddTriaData)

!<description>

  ! Limits the linear part of a dg_T1 element vector.

!</description>

!<input>
! The additional triangulation data
  type(t_additionalTriaData), intent(in):: raddTriaData
!</input>

!<inputoutput>

  ! A vector to limit
  type(t_vectorBlock), intent(inout) :: rvectorBlock
  
!</inputoutput>
  
!</subroutine>

  ! local variables, used by all processors
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: indof, NEL, iel, NVE, ivt, NVT
    
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! The coordinates of the points in which to evaluate the solution vector
  real(dp), dimension(2,17) :: Dpoints
  
  ! The list of elements, in which these points can be found
  integer, dimension(17) :: Ielements
  
  ! The values of the solution vector in the points
  real(dp), dimension(17) :: Dvalues
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  real(dp), dimension(:), allocatable :: duimax, duimin
  
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  
  real(dp) :: xc,yc
  integer :: iidx, nvert, ivert, ineighbour, ineighElm
  integer, dimension(4) :: IhomeIndex
  
  real(dp) :: dui, ddu, duc
  
  integer, dimension(6) :: IdofGlob
  
  integer, dimension(5) :: Isep
  
  integer :: NVBD
  
  integer :: nvar, ivar, idim
  
  real(dp), dimension(:), allocatable :: DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi
  
  real(dp), dimension(:,:), allocatable :: DLin, DtLin, Dalphaei, DL, DR, Dalpha
  
  ! Array of pointers to the data of the blockvector to limit
  type(t_dpPointer), dimension(:), allocatable :: p_DoutputData
  
  integer :: iglobVtNumber, iglobNeighNum
  
  integer :: ilim, ideriv
  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Allocate the space for the pointer to the Data of the different blocks of the output vector
  allocate(p_DoutputData(nvar))
    
  do ivar = 1, nvar
    call lsyssc_getbase_double(rvectorBlock%RvectorBlock(ivar),p_DoutputData(ivar)%p_Ddata)
  end do
    
  
  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  ! Number of vertices per element
  !NVE = p_rtriangulation%NVE
  ! elem_igetNVE(celement)
  NVE = 4
  
  ! Number of vertives at boundary
  NVBD = p_rtriangulation%NVBD
  
  ! Number of vertices
  NVT = p_rtriangulation%NVT
  
  ! Get pointers to the data form the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
                               p_IelementsAtVertex) 
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
                               p_IelementsAtVertexIdx)
  call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary ,&
                               p_IverticesAtBoundary)                               
                              
  ! Allocate the space for solution differences, transformed solution differences,
  ! limited transformed solution differences, limited backtransformed solution differences
  ! and limiting factors
  allocate(DVec(nvar), DVei(nvar), DIi(nvar), DtIi(nvar), DtLinMax(nvar), DtLinMin(nvar), DltIi(nvar), DlIi(nvar))
  !allocate(DLin(nvar,NVE-1), DtLin(nvar,NVE-1), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
  allocate(DLin(nvar,10), DtLin(nvar,10), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
  
  Dalpha = 1.0_DP
  
  do ilim = 1, 3
  
  select case (ilim)
    case (1)
      ideriv = DER_DERIV_X
    case (2)
      ideriv = DER_DERIV_Y
    case (3)
      ideriv = DER_FUNC
  end select
  
  
  do iel = 1, NEL
  
    ! Get coordinates of the center of the element
    Dpoints(1,1) = raddTriaData%p_DmidPoints(1,iel)
    Dpoints(2,1) = raddTriaData%p_DmidPoints(2,iel)
    Ielements(1) = iel
    
    ! Get values in the center of the element for all variables
    do ivar = 1, nvar
      call fevl_evaluate (ideriv, DVec(ivar:ivar), rvectorBlock%RvectorBlock(ivar),&
                          Dpoints(1:2,1:1), Ielements(1:1))
    end do
    
    ! Here we should maybe get a local value of NVE
    
    ! Initialise the correction factor
    Dalphaei(:,:) = 1.0_dp
    
    ! Now calculate the limiting factor for every vertex on our element
    do ivt = 1, NVE
      
      ! Get global vertex number of our local vertex
      iglobVtNumber = p_IverticesAtElement(ivt,iel)
      
      ! Get solution value in this edge in the element
      Ielements(1) = iel
      
      do ivar = 1, nvar
        call fevl_evaluate (DER_FUNC, DVei(ivar:ivar), rvectorBlock%RvectorBlock(ivar),&
                            p_DvertexCoords(1:2,iglobVtNumber:iglobVtNumber), Ielements(1:1))
      end do
      
      ! Calculate solution difference
      DIi = DVei - DVec
      !write(*,*) DIi
      
      ! Get center values of all variables in all neighbour vertices and calculate
      ! the solution differences Dlin(nvar,nneighbors)
      iidx = 0
      DLin = 0.0_dp
      do ineighbour = p_IelementsAtVertexIdx(iglobVtNumber), p_IelementsAtVertexIdx(iglobVtNumber+1)-1
        iglobNeighNum = p_IelementsAtVertex(ineighbour)
        if (iglobNeighNum.ne.iel) then
        
          ! Get midpoint
          Dpoints(1,1) = raddTriaData%p_DmidPoints(1,iglobNeighNum)
          Dpoints(2,1) = raddTriaData%p_DmidPoints(2,iglobNeighNum)
          Ielements(1) = iglobNeighNum
    
          iidx = iidx +1
    
          ! Get values in the center of the element for all variables
          do ivar = 1, nvar
            call fevl_evaluate (ideriv, DLin(ivar:ivar,iidx), rvectorBlock%RvectorBlock(ivar),&
                                Dpoints(1:2,1:1), Ielements(1:1))
          end do
        
          ! Calculate solution difference
          DLin(:,iidx) = DLin(:,iidx) - DVec(:)
          
          
        end if
      end do
      
      iidx=iidx+1
      
      if (iidx.ne.0) then
      ! Dimensional splitting
      do idim = 1, NDIM2D
        ! Now we need the trafo matrices
        DL = buildInvTrafo(DVec,idim)
        DR = buildTrafo(DVec,idim)
        
        ! Transform the solution differences
        DtIi = matmul(DL,DIi)
        do ineighbour = 1, iidx
          DtLin(:,ineighbour) = matmul(DL,Dlin(:,ineighbour))
        end do
        
        ! Get max and min of the transformed solution differences
        do ivar = 1, nvar
          !DtLinMax(ivar) = max(maxval(DtLin(ivar,1:iidx)),0.0_dp)
          !DtLinMin(ivar) = min(minval(DtLin(ivar,1:iidx)),0.0_dp)
          DtLinMax(ivar) = maxval(DtLin(ivar,1:iidx))
          DtLinMin(ivar) = minval(DtLin(ivar,1:iidx))
        end do
        
        ! Now, as the differences are transformed, we can limit every component on its own
        do ivar = 1, nvar
          DltIi(ivar) = max(min(DtLinMax(ivar),DtIi(ivar)),DtLinMin(ivar))
        end do
        
        
        
        ! Now we can trafo back
        DlIi = matmul(DR,DltIi)
        
        ! Calculate the correction factor
        ! for this element, for this edge, for this dimension (take min of all dimensions)
        do ivar = 1, nvar
          if (abs(DIi(ivar))<SYS_EPSREAL) then
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), 1.0_dp)
            ! That's the same as: Do nothing
          else
            Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))
          end if
        end do
      
      end do ! idim
      end if
      
    end do ! ivt
    
    
    
    select case (ilim)
    case (1)
    
      ! Get minimum of all correction factors of all vertices on this element
      do ivar = 1, nvar
        Dalpha(ivar, iel) = minval(Dalphaei(ivar,:))
      end do
      
    case (2)
      
      ! Get minimum of all correction factors of all vertices on this element
      do ivar = 1, nvar
        Dalpha(ivar, iel) = min(Dalpha(ivar, iel),minval(Dalphaei(ivar,:)))
      end do
      
    case (3)
    
      ! Get minimum of all correction factors of all vertices on this element
      do ivar = 1, nvar
        Dalpha(ivar, iel) = max(Dalpha(ivar, iel),minval(Dalphaei(ivar,:)))
      end do
      
    end select
    
  end do !iel
  
  
  ! *** Now limit the solution ***
  
  do iel = 1, NEL
  
  ! Get global DOFs of the element
  call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
  
  select case (ilim)
  case (1)
  
    ! Do nothing      
    
  case (2)
      
    ! Multiply the quadratic part of the solution vector with the correction factor
    do ivar = 1, nvar
      p_DoutputData(ivar)%p_Ddata(IdofGlob(4:6)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(4:6))*Dalpha(ivar, iel)
    end do

      
  case (3)

   ! Multiply the linear part of the solution vector with the correction factor
    do ivar = 1, nvar
      p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3))*Dalpha(ivar, iel)
    end do
    

  end select
  
    end do ! iel
  
  
  end do ! ilim

  
  deallocate(p_DoutputData)
  deallocate(DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi)
  deallocate(DLin, DtLin, Dalphaei)

  
  end subroutine
  
  
    
  
  subroutine profiler_init(rprofiler, inumtimers)
    type(t_profiler), intent(inout) :: rprofiler
    integer, intent(in) :: inumtimers
    
    allocate(rprofiler%Dtimers(inumtimers))
    rprofiler%Dtimers(:)    = 0.0_dp
    rprofiler%dendtime      = 0.0_dp
    rprofiler%dlasttime     = 0.0_dp
    rprofiler%ntimer        = inumtimers
    rprofiler%icurrenttimer = -1
    call cpu_time(rprofiler%dstarttime)
    
  end subroutine
  
  
  subroutine profiler_release(rprofiler)
    type(t_profiler), intent(inout) :: rprofiler
    
    real(dp) :: dalltime
    integer :: i
    
    call cpu_time(rprofiler%dendtime)
    
    dalltime = rprofiler%dendtime - rprofiler%dstarttime
    
    write(*,*) ''
    write(*,*) '*********************************************************************'
    write(*,*) 'Profiler statistics'
    write(*,*) 'Full time:' , dalltime
    
    do i = 1, rprofiler%ntimer
      write(*,*) i,  (rprofiler%Dtimers(i)/(dalltime+SYS_EPSREAL)*100),'%'
    end do
    write(*,*) '*********************************************************************'
    write(*,*) ''
    
    rprofiler%dstarttime    = 0.0_dp
    rprofiler%dendtime      = 0.0_dp
    rprofiler%dlasttime     = 0.0_dp
    rprofiler%ntimer        = -1
    rprofiler%icurrenttimer = -1
    deallocate(rprofiler%Dtimers)
    
  end subroutine
  
  
  subroutine profiler_measure(rprofiler,icurrtimer)
    type(t_profiler), intent(inout) :: rprofiler
    integer, intent(in) :: icurrtimer
    
    real(dp) :: dcurrtime
    integer :: i
    
    call cpu_time(dcurrtime)
    
    if (rprofiler%icurrenttimer>0) then
      rprofiler%Dtimers(rprofiler%icurrenttimer) = rprofiler%Dtimers(rprofiler%icurrenttimer) + dcurrtime-rprofiler%dlasttime
    end if
    
    rprofiler%dlasttime = dcurrtime
    
    rprofiler%icurrenttimer = icurrtimer
    
  end subroutine
  
  
  
!****************************************************************************
  
!<subroutine>  
  
  subroutine getDtByCfl (rvectorBlock, raddTriaData, dCFL, dt, dgravconst)

!<description>

  ! Gets the timestepsize dt to fit the given CFL number

!</description>

!<input>
! The additional triangulation data
  type(t_additionalTriaData), intent(in):: raddTriaData
  type(t_vectorBlock), intent(in) :: rvectorBlock
  real(dp), intent(in) :: dCFL, dgravconst
!</input>

!<output>
  real(dp), intent(out) :: dt
!</output>
  
!</subroutine>
  
  ! Local variables
  
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  integer(I32) :: celement
  
  integer :: NEL, iel
  
  integer, dimension(:), allocatable :: IelIdx
  integer, dimension(:,:), allocatable :: IdofGlob
  
  real(dp), dimension(:), pointer :: p_Ddata
  
  real(dp), dimension(3) :: DQ
  
  real(dp) :: dh, du, dv, dc, ddx, ddy

  
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation

  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  celement = p_rspatialDiscr%RelementDistr(1)%celement
  
  allocate(IelIdx(NEL))
  
  do iel = 1, NEL
    IelIdx(iel) = iel
  end do
  
  allocate(IdofGlob(elem_igetNDofLoc(celement),NEL))
  
  ! Get global DOFs of the elements  
  call dof_locGlobMapping_mult(p_rspatialDiscr, IelIdx, IdofGlob)
  
  call lsysbl_getbase_double (rvectorBlock, p_Ddata)
  
  dt = SYS_MAXREAL
  
  do iel = 1, NEL
    
    DQ = p_Ddata(rvectorBlock%RvectorBlock(:)%iidxFirstEntry+IdofGlob(1,iel)-1)
    
    dh = DQ(1)
    du = abs(DQ(2)/dh)
    dv = abs(DQ(3)/dh)
    dc = sqrt(dh*dgravconst)
    ddx = raddTriaData%p_Ddxdy(1,iel)
    ddy = raddTriaData%p_Ddxdy(2,iel)
    
    dt = min(dt,dCFL/((du+dc)/ddx+(dv+dc)/ddy))
    
  end do
  
  
  ! Deallocate memory
  deallocate(IdofGlob,IelIdx)
  
end subroutine



  
  


!****************************************************************************
  
!<subroutine>  
  
  subroutine dg_linearLimiterBlockCharVar_mixedJacobian (rvectorBlock)

!<description>

  ! Limits the linear part of a dg_T1 element vector.

!</description>

!<input>
!</input>

!<inputoutput>

  ! A vector to limit
  type(t_vectorBlock), intent(inout) :: rvectorBlock
  
!</inputoutput>
  
!</subroutine>

  ! local variables, used by all processors
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: indof, NEL, iel, NVE, ivt, NVT
    
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! The coordinates of the points in which to evaluate the solution vector
  real(dp), dimension(2,17) :: Dpoints
  
  ! The list of elements, in which these points can be found
  integer, dimension(17) :: Ielements
  
  ! The values of the solution vector in the points
  real(dp), dimension(17) :: Dvalues
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  real(dp), dimension(:), allocatable :: duimax, duimin
  
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  
  real(dp) :: xc,yc
  integer :: iidx, nvert, ivert, ineighbour, ineighElm
  integer, dimension(4) :: IhomeIndex
  
  real(dp) :: dui, ddu, duc
  
  integer, dimension(3) :: IdofGlob
  
  integer, dimension(5) :: Isep
  
  integer :: NVBD
  
  integer :: nvar, ivar, idim
  
  real(dp), dimension(:), allocatable :: DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi
  
  real(dp), dimension(:,:), allocatable :: DLin, DtLin, Dalphaei, DL, DR, Dalpha
    
  ! Array of pointers to the data of the blockvector to limit
  type(t_dpPointer), dimension(:), allocatable :: p_DoutputData
  
  integer :: iglobVtNumber, iglobNeighNum
  
  real(dp) :: da, db, dquo
  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Allocate the space for the pointer to the Data of the different blocks of the output vector
  allocate(p_DoutputData(nvar))
    
!  do ivar = 1, nvar
!    call lsyssc_getbase_double(rvectorBlock%RvectorBlock(ivar),p_DoutputData(ivar)%p_Ddata)
!  end do
  call lsysbl_getbase_double (rvectorBlock, p_Ddata)
  
  
    
  
  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  ! Number of vertices per element
  !NVE = p_rtriangulation%NVE
  ! elem_igetNVE(celement)
  NVE = 4
  
  ! Number of vertives at boundary
  NVBD = p_rtriangulation%NVBD
  
  ! Number of vertices
  NVT = p_rtriangulation%NVT
  
  ! Get pointers to the data form the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
                               p_IelementsAtVertex) 
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
                               p_IelementsAtVertexIdx)
  call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary ,&
                               p_IverticesAtBoundary)                               
                              
  ! Allocate the space for solution differences, transformed solution differences,
  ! limited transformed solution differences, limited backtransformed solution differences
  ! and limiting factors
  allocate(DVec(nvar), DVei(nvar), DIi(nvar), DtIi(nvar), DtLinMax(nvar), DtLinMin(nvar), DltIi(nvar), DlIi(nvar))
  !allocate(DLin(nvar,NVE-1), DtLin(nvar,NVE-1), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
  allocate(DLin(nvar,10), DtLin(nvar,10), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
  
  do iel = 1, NEL
    
    ! Get global DOFs of the element
    call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    ! Get values in the center of the element for all variables
    do ivar = 1, nvar
      !DVec(ivar) = p_DoutputData(ivar)%p_Ddata(IdofGlob(1))
      DVec(ivar) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1)-1)
    end do
    
    ! Here we should maybe get a local value of NVE
    
    ! Initialise the correction factor
    !Dalphaei(:,:) = 1.0_dp
    Dalphaei(:,:) = 0.0_dp
    
    ! Now calculate the limiting factor for every vertex on our element
    do ivt = 1, NVE
      
      ! Get global vertex number of our local vertex
      iglobVtNumber = p_IverticesAtElement(ivt,iel)
      
      ! Get solution value in this edge in the element
      Ielements(1) = iel
      
!      do ivar = 1, nvar
!        call fevl_evaluate (DER_FUNC, DVei(ivar:ivar), rvectorBlock%RvectorBlock(ivar),&
!                            p_DvertexCoords(1:2,iglobVtNumber:iglobVtNumber), Ielements(1:1))
!      end do
      call dg_evaluateLinearPart (rvectorBlock, iel, iglobVtNumber, DVei)
      
      ! Calculate solution difference
      DIi = DVei - DVec
      !write(*,*) DIi
      
      ! Get center values of all variables in all neighbour vertices and calculate
      ! the solution differences Dlin(nvar,nneighbors)
      iidx = 0
      DLin = 0.0_dp
      do ineighbour = p_IelementsAtVertexIdx(iglobVtNumber), p_IelementsAtVertexIdx(iglobVtNumber+1)-1
        iglobNeighNum = p_IelementsAtVertex(ineighbour)
        !if (iglobNeighNum.ne.iel) then
          call dof_locGlobMapping(p_rspatialDiscr, iglobNeighNum, IdofGlob)
          iidx = iidx + 1          
          do ivar = 1, nvar
            !DLin(ivar,iidx) = p_DoutputData(ivar)%p_Ddata(IdofGlob(1)) - DVec(ivar)
            DLin(ivar,iidx) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1)-1) - DVec(ivar)
          end do

          
        !end if
      end do
          
      if (iidx>0) then
      ! Dimensional splitting
      do idim = 3, 4
        ! Now we need the trafo matrices
        if(idim<3) then
        DL = buildInvTrafo(DVec,idim)
        DR = buildTrafo(DVec,idim)
        else if (idim==3) then
        da = DVec(2)/DVec(1)
        db = DVec(3)/DVec(1)
        dquo = da*da+db*db
        if (dquo<SYS_EPSREAL) then
          DL = buildInvTrafo(DVec,idim-2)
          DR = buildTrafo(DVec,idim-2)
        else
          da = da/dquo
          db = db/dquo
          DL = buildMixedL2(DVec,da,db)
          DR = buildMixedR2(DVec,da,db)
        end if
        
        else if (idim==4) then
        da = DVec(2)/DVec(1)
        db = DVec(3)/DVec(1)
        dquo = da*da+db*db
        if (dquo<SYS_EPSREAL) then
          DL = buildInvTrafo(DVec,idim-2)
          DR = buildTrafo(DVec,idim-2)
        else
          da = da/dquo
          db = db/dquo
          DL = buildMixedL2(DVec,-db,da)
          DR = buildMixedR2(DVec,-db,da)
        end if
        end if
        
        
        ! Transform the solution differences
        DtIi = matmul(DL,DIi)
        do ineighbour = 1, iidx
          DtLin(:,ineighbour) = matmul(DL,Dlin(:,ineighbour))
        end do
        

        
        ! Get max and min of the transformed solution differences
        do ivar = 1, nvar
          DtLinMax(ivar) = maxval(DtLin(ivar,1:iidx))
          DtLinMin(ivar) = minval(DtLin(ivar,1:iidx))
          !DtLinMax(ivar) = max(DtLinMax(ivar),0.0_dp)
          !DtLinMin(ivar) = min(DtLinMin(ivar),0.0_dp)
        end do
       
        ! Now, as the differences are transformed, we can limit every component on its own
        do ivar = 1, nvar
          DltIi(ivar) = max(min(DtLinMax(ivar),DtIi(ivar)),DtLinMin(ivar))
        end do
        
        
        
        ! Now we can trafo back
        DlIi = matmul(DR,DltIi)
        
        ! Calculate the correction factor
        ! for this element, for this edge, for this dimension (take min of all dimensions)
        do ivar = 1, nvar
          if (abs(DIi(ivar))<SYS_EPSREAL) then
            !Dalphaei(ivar,ivt) = 1.0_dp
          else
            ! This is the one following the principles
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))
            
            ! This one is less limiting
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), min(abs(DlIi(ivar)/DIi(ivar)),1.0_dp ))
            
            ! This is the one with the max of the dimensional splitting
            Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))
            
            ! This is the least limiting one
            !Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt), min(abs(DlIi(ivar)/DIi(ivar)),1.0_dp ) )
          end if
          
          ! No limiting at boundary
          if (iidx<3) Dalphaei(ivar,ivt) = 1.0_dp
          
        end do
      
      end do ! idim
      end if

      
    end do ! ivt


    
    ! Get minimum of all correction factors of all vertices on this element
    do ivar = 1, nvar
      Dalpha(ivar, iel) = minval(Dalphaei(ivar,:))
    end do
    
    
!    if (iel.eq.985) then
!      write(*,*) Dalphaei(1,:)
!      write(*,*) Dalpha(1, iel)
!    end if
    
    ! Get global DOFs of the element
    call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
    
    ! Multiply the linear part of the solution vector with the correction factor
    do ivar = 1, nvar
      ! p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3))*Dalpha(ivar, iel)
      
      p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3)-1) = &            
        p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3)-1) * &
        Dalpha(ivar, iel)
      
    end do
    
  end do !iel

  
  deallocate(p_DoutputData)
  deallocate(DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi)
  deallocate(DLin, DtLin, Dalphaei)

  
  end subroutine
  
  
  
  
!****************************************************************************
  
!<subroutine>  
  
  subroutine dg_evaluateLinearPart (rvectorBlock, iel, iglobVtNumber, Dvalues)

!<description>

  ! Evaluates the linear part of an dg_T1 or dg_T2 element in element # iel and on the vertex # iglobVtNumber.

!</description>

!<input>
type(t_vectorBlock), intent(in) :: rvectorBlock
integer, intent(in) :: iel, iglobVtNumber

!</input>

!<output>
real(dp), dimension(:) :: Dvalues
!</output>
  
!</subroutine>

  ! Local variables
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: ivt, nvar, ivar
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  integer, dimension(17) :: IdofGlob
  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Get pointer to the data of the (solution) vector
  call lsysbl_getbase_double (rvectorBlock, p_Ddata)
  
  ! Get pointer to the data form the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                             p_IverticesAtElement)

  ! Find local vertex number
  do ivt = 1, 4 ! 3 shoud be enough, if it has to be really fast
    if (p_IverticesAtElement(ivt,iel).eq.iglobVtNumber) exit
  end do
  
  ! Get global DOFs of the element
  call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
  
  select case (ivt)
    case (1)
      do ivar = 1, nvar
        Dvalues(ivar) = &
            p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1)-1) &
            - p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2)-1) &
            - p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3)-1)
      end do
      case (2)
      do ivar = 1, nvar
        Dvalues(ivar) = &
            p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1)-1) &
            + p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2)-1) &
            - p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3)-1)
      end do
      case (3)
      do ivar = 1, nvar
        Dvalues(ivar) = &
            p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1)-1) &
            + p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2)-1) &
            + p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3)-1)
      end do
      case (4) ! case default would increase speed a little ?!
      do ivar = 1, nvar
        Dvalues(ivar) = &
            p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1)-1) &
            - p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2)-1) &
            + p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3)-1)
      end do
      
  end select
  
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_quadraticLimiterBlockCharVar_mixedJacobian (rvectorBlock, raddTriaData)

!<description>

  ! Limits the linear part of a dg_T1 element vector.

!</description>

!<input>
! The additional triangulation data
  type(t_additionalTriaData), intent(in):: raddTriaData
!</input>

!<inputoutput>

  ! A vector to limit
  type(t_vectorBlock), intent(inout) :: rvectorBlock
  
!</inputoutput>
  
!</subroutine>

  ! local variables, used by all processors
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: indof, NEL, iel, NVE, ivt, NVT
    
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! The coordinates of the points in which to evaluate the solution vector
  real(dp), dimension(2,17) :: Dpoints
  
  ! The list of elements, in which these points can be found
  integer, dimension(17) :: Ielements
  
  ! The values of the solution vector in the points
  real(dp), dimension(17) :: Dvalues
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  real(dp), dimension(:), allocatable :: duimax, duimin
  
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  
  real(dp) :: xc,yc
  integer :: iidx, nvert, ivert, ineighbour, ineighElm
  integer, dimension(4) :: IhomeIndex
  
  real(dp) :: dui, ddu, duc
  
  integer, dimension(6) :: IdofGlob
  
  integer, dimension(5) :: Isep
  
  integer :: NVBD
  
  integer :: nvar, ivar, idim
  
  real(dp), dimension(:), allocatable :: DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi, DQchar
  
  real(dp), dimension(:,:), allocatable :: DLin, DtLin, Dalphaei, DL, DR, Dalpha
  
  ! Array of pointers to the data of the blockvector to limit
  type(t_dpPointer), dimension(:), allocatable :: p_DoutputData
  
  integer :: iglobVtNumber, iglobNeighNum
  
  integer :: ilim, ideriv
  
  real(dp) :: da, db, dquo
  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Allocate the space for the pointer to the Data of the different blocks of the output vector
!  allocate(p_DoutputData(nvar))
!    
!  do ivar = 1, nvar
!    call lsyssc_getbase_double(rvectorBlock%RvectorBlock(ivar),p_DoutputData(ivar)%p_Ddata)
!  end do

  ! Get pointer to the data of the (solution) vector
  call lsysbl_getbase_double (rvectorBlock, p_Ddata)
    
  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  ! Number of vertices per element
  !NVE = p_rtriangulation%NVE
  ! elem_igetNVE(celement)
  NVE = 4
  
  ! Number of vertives at boundary
  NVBD = p_rtriangulation%NVBD
  
  ! Number of vertices
  NVT = p_rtriangulation%NVT
  
  ! Get pointers to the data form the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
                               p_IelementsAtVertex) 
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
                               p_IelementsAtVertexIdx)
  call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary ,&
                               p_IverticesAtBoundary)                               
                              
  ! Allocate the space for solution differences, transformed solution differences,
  ! limited transformed solution differences, limited backtransformed solution differences
  ! and limiting factors
  allocate(DVec(nvar), DVei(nvar), DIi(nvar), DtIi(nvar), DtLinMax(nvar), DtLinMin(nvar), DltIi(nvar), DlIi(nvar),DQchar(nvar))
  !allocate(DLin(nvar,NVE-1), DtLin(nvar,NVE-1), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
  allocate(DLin(nvar,10), DtLin(nvar,10), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
  
  Dalpha = 1.0_DP
  
  do ilim = 1, 3
  
  select case (ilim)
    case (1)
      ideriv = DER_DERIV_X
    case (2)
      ideriv = DER_DERIV_Y
    case (3)
      ideriv = DER_FUNC
  end select
  
  
  do iel = 1, NEL
  
    ! Get coordinates of the center of the element
    Dpoints(1,1) = raddTriaData%p_DmidPoints(1,iel)
    Dpoints(2,1) = raddTriaData%p_DmidPoints(2,iel)
    Ielements(1) = iel
    
    ! Get values in the center of the element for all variables
    
    if (ilim<3) then
      do ivar = 1, nvar
        call fevl_evaluate (ideriv, DVec(ivar:ivar), rvectorBlock%RvectorBlock(ivar),&
                            Dpoints(1:2,1:1), Ielements(1:1))
      end do
    else
      ! Get global DOFs of the element
      call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
      
      do ivar = 1, nvar
        !DVec(ivar) = p_DoutputData(ivar)%p_Ddata(IdofGlob(1))
        DVec(ivar) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1)-1)
      end do
    end if
    
    
    ! Get global DOFs of the element
    call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
      
    do ivar = 1, nvar
      DQchar(ivar) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1)-1)
    end do
    
    ! Here we should maybe get a local value of NVE
    
    ! Initialise the correction factor
    Dalphaei(:,:) = 0.0_dp
    
    ! Now calculate the limiting factor for every vertex on our element
    do ivt = 1, NVE
      
      ! Get global vertex number of our local vertex
      iglobVtNumber = p_IverticesAtElement(ivt,iel)
      
      ! Get solution value in this edge in the element
      Ielements(1) = iel
      
      if (ilim<3) then
        do ivar = 1, nvar
          call fevl_evaluate (ideriv, DVei(ivar:ivar), rvectorBlock%RvectorBlock(ivar),&
                              p_DvertexCoords(1:2,iglobVtNumber:iglobVtNumber), Ielements(1:1))
        end do
      else
        call dg_evaluateLinearPart (rvectorBlock, iel, iglobVtNumber, DVei)
      end if
      
      ! Calculate solution difference
      DIi = DVei - DVec
      !write(*,*) DIi
      
      ! Get center values of all variables in all neighbour vertices and calculate
      ! the solution differences Dlin(nvar,nneighbors)
      iidx = 0
      DLin = 0.0_dp
      do ineighbour = p_IelementsAtVertexIdx(iglobVtNumber), p_IelementsAtVertexIdx(iglobVtNumber+1)-1
        iglobNeighNum = p_IelementsAtVertex(ineighbour)
        !if (iglobNeighNum.ne.iel) then
        
          ! Get midpoint
          Dpoints(1,1) = raddTriaData%p_DmidPoints(1,iglobNeighNum)
          Dpoints(2,1) = raddTriaData%p_DmidPoints(2,iglobNeighNum)
          Ielements(1) = iglobNeighNum
    
          iidx = iidx +1
    
          ! Get values in the center of the element for all variables
          do ivar = 1, nvar
            call fevl_evaluate (ideriv, DLin(ivar:ivar,iidx), rvectorBlock%RvectorBlock(ivar),&
                                Dpoints(1:2,1:1), Ielements(1:1))
          end do
        
          ! Calculate solution difference
          DLin(:,iidx) = DLin(:,iidx) - DVec(:)
          
        !end if
      end do
      
     
      
      if (iidx.ne.0) then
      ! Dimensional splitting
      do idim = 3, 4
        ! Now we need the trafo matrices
        if(idim<3) then
        DL = buildInvTrafo(DQchar,idim)
        DR = buildTrafo(DQchar,idim)
        else if (idim==3) then
        da = DQchar(2)/DQchar(1)
        db = DQchar(3)/DQchar(1)
        dquo = da*da+db*db
        if (dquo<SYS_EPSREAL) then
          DL = buildInvTrafo(DQchar,idim-2)
          DR = buildTrafo(DQchar,idim-2)
        else
          da = da/dquo
          db = db/dquo
          DL = buildMixedL2(DQchar,da,db)
          DR = buildMixedR2(DQchar,da,db)
        end if
        
        else if (idim==4) then
        da = DQchar(2)/DQchar(1)
        db = DQchar(3)/DQchar(1)
        dquo = da*da+db*db
        if (dquo<SYS_EPSREAL) then
          DL = buildInvTrafo(DQchar,idim-2)
          DR = buildTrafo(DQchar,idim-2)
        else
          da = da/dquo
          db = db/dquo
          DL = buildMixedL2(DQchar,-db,da)
          DR = buildMixedR2(DQchar,-db,da)
        end if
        end if
        
        ! Transform the solution differences
        DtIi = matmul(DL,DIi)
        do ineighbour = 1, iidx
          DtLin(:,ineighbour) = matmul(DL,Dlin(:,ineighbour))
        end do
        
        ! Get max and min of the transformed solution differences
        do ivar = 1, nvar
          !DtLinMax(ivar) = max(maxval(DtLin(ivar,1:iidx)),0.0_dp)
          !DtLinMin(ivar) = min(minval(DtLin(ivar,1:iidx)),0.0_dp)
          DtLinMax(ivar) = maxval(DtLin(ivar,1:iidx))
          DtLinMin(ivar) = minval(DtLin(ivar,1:iidx))
        end do
        
        ! Now, as the differences are transformed, we can limit every component on its own
        do ivar = 1, nvar
          DltIi(ivar) = max(min(DtLinMax(ivar),DtIi(ivar)),DtLinMin(ivar))
        end do
        
        
        
        ! Now we can trafo back
        DlIi = matmul(DR,DltIi)
        
        ! Calculate the correction factor
        ! for this element, for this edge, for this dimension (take min of all dimensions)
        do ivar = 1, nvar
          if (abs(DIi(ivar))<SYS_EPSREAL) then
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), 1.0_dp)
            ! That's the same as: Do nothing
          else
            ! This is the one following the principles
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))
            
            ! This one is less limiting
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), min(abs(DlIi(ivar)/DIi(ivar)),1.0_dp ))
            
            ! This is the one with the max of the dimensional splitting
            Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))
            
            ! This is the least limiting one
            !Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt), min(abs(DlIi(ivar)/DIi(ivar)),1.0_dp ) )
            
          end if
          
          ! No limiting at boundary
          if (iidx<3) then
            Dalphaei(ivar,ivt) = 1.0_dp
          end if
          
        end do
      
      end do ! idim
      end if
      
    end do ! ivt
    
    
    
    select case (ilim)
    case (1)
    
      ! Get minimum of all correction factors of all vertices on this element
      do ivar = 1, nvar
        Dalpha(ivar, iel) = minval(Dalphaei(ivar,1:NVE))
      end do
      
    case (2)
      
      ! Get minimum of all correction factors of all vertices on this element
      do ivar = 1, nvar
        Dalpha(ivar, iel) = min(Dalpha(ivar, iel),minval(Dalphaei(ivar,1:NVE)))
      end do
      
    case (3)
    
      ! Get minimum of all correction factors of all vertices on this element
      do ivar = 1, nvar
        Dalpha(ivar, iel) = max(Dalpha(ivar, iel),minval(Dalphaei(ivar,1:NVE)))
      end do
      
    end select
    
  end do !iel
  
  
  ! *** Now limit the solution ***
  do iel = 1, NEL
  
  ! Get global DOFs of the element
  call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
  
  select case (ilim)
  case (1)
  
    ! Do nothing      
    
  case (2)
      
    ! Multiply the quadratic part of the solution vector with the correction factor
    do ivar = 1, nvar
      !p_DoutputData(ivar)%p_Ddata(IdofGlob(4:6)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(4:6))*Dalpha(ivar, iel)
      p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(4:6)-1) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(4:6)-1)*Dalpha(ivar, iel)
    end do

      
  case (3)

   ! Multiply the linear part of the solution vector with the correction factor
    do ivar = 1, nvar
      !p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3))*Dalpha(ivar, iel)
      p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3)-1) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3)-1)*Dalpha(ivar, iel)
    end do
    

  end select
  
    end do ! iel
  
  
  end do ! ilim

  
  !deallocate(p_DoutputData)
  deallocate(DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi, DQchar)
  deallocate(DLin, DtLin, Dalphaei)

  
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_kuzminLimiterBlockCharVar_mixedJacobian (rvectorBlock, raddTriaData)

!<description>

  ! Limits a dg_T1 or dg_T2 element vector.

!</description>

!<input>
! The additional triangulation data
  type(t_additionalTriaData), intent(in):: raddTriaData
!</input>

!<inputoutput>

  ! A vector to limit
  type(t_vectorBlock), intent(inout) :: rvectorBlock
  
!</inputoutput>
  
!</subroutine>

  ! local variables, used by all processors
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: indof, NEL, iel, NVE, ivt, NVT
    
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! The coordinates of the points in which to evaluate the solution vector
  real(dp), dimension(2,17) :: Dpoints
  
  ! The list of elements, in which these points can be found
  integer, dimension(17) :: Ielements
  
  ! The values of the solution vector in the points
  real(dp), dimension(17) :: Dvalues
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  
  real(dp), dimension(:), allocatable :: duimax, duimin
  
  integer, dimension(:), pointer :: p_IverticesAtBoundary
  
  real(dp) :: xc,yc
  integer :: iidx, nvert, ivert, ineighbour, ineighElm
  integer, dimension(4) :: IhomeIndex
  
  real(dp) :: dui, ddu, duc
  
  integer, dimension(:,:), allocatable :: IdofGlob

  integer :: NVBD
  
  integer :: nvar, ivar, idim
  
  real(dp), dimension(:), allocatable :: DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi, DQchar
  
  real(dp), dimension(:,:), allocatable :: DLin, DtLin, Dalphaei, DL, DR, Dalpha
  
  ! Array of pointers to the data of the blockvector to limit
  type(t_dpPointer), dimension(:), allocatable :: p_DoutputData
  
  integer :: iglobVtNumber, iglobNeighNum
  
  integer :: ilim, ideriv, ilimstart
  
  integer :: nDOFloc
 
  real(dp) :: da, db, dquo
  
  integer, parameter :: nelemSim = 1000
  
  integer(I32) :: celement
  
  integer :: npoints
  
  integer, dimension(:), pointer :: p_IelIdx
  
  ! Values of the linear part of the solution, Dimension (nvar,# points per element (midpoint+corners),NEL)
  real(dp), dimension(:,:,:), allocatable :: DLinPartValues
  
  ! Dimension (nvar,# points per element (midpoint+corners),# derivatives (2, x- and y-),NEL)
  real(dp), dimension(:,:,:,:), allocatable :: DDerivativeQuadraticValues
  
  integer, parameter :: nmaxneighbours = 10
  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
  
  ! Get pointers for quicker access
  p_rspatialDiscr => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
  
  ! Allocate the space for the pointer to the Data of the different blocks of the output vector
!  allocate(p_DoutputData(nvar))
!    
!  do ivar = 1, nvar
!    call lsyssc_getbase_double(rvectorBlock%RvectorBlock(ivar),p_DoutputData(ivar)%p_Ddata)
!  end do

  ! Get pointer to the data of the (solution) vector
  call lsysbl_getbase_double (rvectorBlock, p_Ddata)
    
  ! Get number of elements
  NEL = p_rtriangulation%NEL
  
  ! What is the current element type?
  celement = p_rspatialDiscr%RelementDistr(1)%celement
  
  ! Number of vertices per element
  NVE = elem_igetNVE(celement)
  !NVE = 4
  
  ! Get number of local DOF
  nDOFloc = elem_igetNDofLoc(celement)
  
  ! Number of vertices
  NVT = p_rtriangulation%NVT
  
  ! Get pointers to the data form the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
                               p_IelementsAtVertex) 
  call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
                               p_IelementsAtVertexIdx)
  call storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary ,&
                               p_IverticesAtBoundary)                               
                              
  ! Allocate the space for solution differences, transformed solution differences,
  ! limited transformed solution differences, limited backtransformed solution differences
  ! and limiting factors
  allocate(DVec(nvar), DVei(nvar), DIi(nvar), DtIi(nvar), DtLinMax(nvar), DtLinMin(nvar), DltIi(nvar), DlIi(nvar),DQchar(nvar))
  !allocate(DLin(nvar,NVE-1), DtLin(nvar,NVE-1), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
  allocate(DLin(nvar,10), DtLin(nvar,10), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
  
!  ! Get the number of elements
!  NEL = p_rspatialDiscr%RelementDistr(1)%NEL
  
  ! Allocate the space for the global DOF
  allocate(IdofGlob(nDOFloc,NEL))
  
  ! Number of points on each element, where the solution has to be evaluated (corner vertives + midpoint)
  npoints = 1 + NVE
  
  ! Allocate space for the values of the linear part of the solution vector 
  allocate(DLinPartValues(nvar,npoints,NEL))
 
  ! Get list of (all) elements
  call storage_getbase_int(p_rspatialDiscr%RelementDistr(1)%h_IelementList, p_IelIdx)
    
  ! Get global DOFs for these (all) elements
  call dof_locGlobMapping_mult(p_rspatialDiscr, p_IelIdx, IdofGlob)
  
  ! Get the values of the linear part of the solution vector 
  call dg_evaluateLinearPart_mult(rvectorBlock, p_IelIdx, IdofGlob, DLinPartValues)
  
  ! Test, if we have quadratic elements
  if (celement == EL_DG_T2_2D) then
    ! Allocate space for the values of the derivatives of the quadratic solution vector 
    allocate(DDerivativeQuadraticValues(nvar,npoints,2,NEL))
  
    ! Get the values of the derivatives of the quadratic solution vector 
    call dg_evaluateDerivativeQuadratic_mult (rvectorBlock, p_IelIdx, DDerivativeQuadraticValues, raddTriaData)
  end if
  
  
  select case (celement)
  case (EL_DG_T2_2D) 
    ilimstart = 1
    Dalpha = 1.0_DP
  case (EL_DG_T1_2D) 
    ilimstart = 3
    Dalpha = 0.0_DP
  end select
  
    
!  ! Initialise the limiting factors
!  Dalpha = 1.0_DP
!  
!  do iel = 1, NEL
!    
!    
!    do ilim = ilimstart, 3
!      
!      ! Get solution values, that are used to evaluate the trofo matrices
!      DQchar = DLinPartValues(:,1,iel)
!      
!      do idim = 3, 4
!        ! Now we need the trafo matrices
!        if(idim<3) then
!        DL = buildInvTrafo(DQchar,idim)
!        DR = buildTrafo(DQchar,idim)
!        else if (idim==3) then
!        da = DQchar(2)/DQchar(1)
!        db = DQchar(3)/DQchar(1)
!        dquo = da*da+db*db
!        if (dquo<SYS_EPSREAL) then
!          DL = buildInvTrafo(DQchar,idim-2)
!          DR = buildTrafo(DQchar,idim-2)
!        else
!          da = da/dquo
!          db = db/dquo
!          DL = buildMixedL2(DQchar,da,db)
!          DR = buildMixedR2(DQchar,da,db)
!        end if
!        
!        else if (idim==4) then
!        da = DQchar(2)/DQchar(1)
!        db = DQchar(3)/DQchar(1)
!        dquo = da*da+db*db
!        if (dquo<SYS_EPSREAL) then
!          DL = buildInvTrafo(DQchar,idim-2)
!          DR = buildTrafo(DQchar,idim-2)
!        else
!          da = da/dquo
!          db = db/dquo
!          DL = buildMixedL2(DQchar,-db,da)
!          DR = buildMixedR2(DQchar,-db,da)
!        end if
!        end if
!        
!        
!        
!        
!      end do ! idim
!  
!    
!    
!    end do ! ilim
!  
!  
!  end do ! iel
  
  
  
  
  
  
  
  
  
  
  
!  Dalpha = 1.0_DP
  
  do ilim = ilimstart, 3
  
  do iel = 1, NEL
   
    ! Get values in the center of the element for all variables
    if (ilim<3) then
      DVec = DDerivativeQuadraticValues(:,1,ilim,iel)
    else
      DVec = DLinPartValues(:,1,iel)
    end if

    DQchar = DLinPartValues(:,1,iel)

    ! Here we should maybe get a local value of NVE

    ! Initialise the correction factor
    Dalphaei(:,:) = 0.0_dp
    
    ! Now calculate the limiting factor for every vertex on our element
    do ivt = 1, NVE
      
      ! Get global vertex number of our local vertex
      iglobVtNumber = p_IverticesAtElement(ivt,iel)
      
      if (ilim<3) then
        DVei = DDerivativeQuadraticValues(:,1+ivt,ilim,iel)
      else
        DVei = DLinPartValues(:,1+ivt,iel)
      end if
      
      ! Calculate solution difference
      DIi = DVei - DVec
      
      ! Get center values of all variables in all neighbour vertices and calculate
      ! the solution differences Dlin(nvar,nneighbors)
      iidx = 0
      DLin = 0.0_dp
      do ineighbour = p_IelementsAtVertexIdx(iglobVtNumber), p_IelementsAtVertexIdx(iglobVtNumber+1)-1
        iglobNeighNum = p_IelementsAtVertex(ineighbour)
        !if (iglobNeighNum.ne.iel) then
        
        iidx = iidx +1
    
        ! Get values in the center of the element for all variables
        

        if (ilim<3) then
          DLin(:,iidx) = DDerivativeQuadraticValues(:,1,ilim,iglobNeighNum)
        else
          DLin(:,iidx) = DLinPartValues(:,1,iglobNeighNum)
        end if
        
        ! Calculate solution difference
        DLin(:,iidx) = DLin(:,iidx) - DVec(:)
          
        !end if
      end do
      
     
      
      if (iidx.ne.0) then
      ! Dimensional splitting
      do idim = 1, 2
        ! Now we need the trafo matrices
        if(idim<3) then
        DL = buildInvTrafo(DQchar,idim)
        DR = buildTrafo(DQchar,idim)
        else if (idim==3) then
        da = DQchar(2)/DQchar(1)
        db = DQchar(3)/DQchar(1)
        dquo = da*da+db*db
        if (dquo<SYS_EPSREAL) then
          DL = buildInvTrafo(DQchar,idim-2)
          DR = buildTrafo(DQchar,idim-2)
        else
          da = da/dquo
          db = db/dquo
          DL = buildMixedL2(DQchar,da,db)
          DR = buildMixedR2(DQchar,da,db)
        end if
        
        else if (idim==4) then
        da = DQchar(2)/DQchar(1)
        db = DQchar(3)/DQchar(1)
        dquo = da*da+db*db
        if (dquo<SYS_EPSREAL) then
          DL = buildInvTrafo(DQchar,idim-2)
          DR = buildTrafo(DQchar,idim-2)
        else
          da = da/dquo
          db = db/dquo
          DL = buildMixedL2(DQchar,-db,da)
          DR = buildMixedR2(DQchar,-db,da)
        end if
        end if
        
        ! Transform the solution differences
        DtIi = matmul(DL,DIi)
        do ineighbour = 1, iidx
          DtLin(:,ineighbour) = matmul(DL,Dlin(:,ineighbour))
        end do
        
        ! Get max and min of the transformed solution differences
        do ivar = 1, nvar
          !DtLinMax(ivar) = max(maxval(DtLin(ivar,1:iidx)),0.0_dp)
          !DtLinMin(ivar) = min(minval(DtLin(ivar,1:iidx)),0.0_dp)
          DtLinMax(ivar) = maxval(DtLin(ivar,1:iidx))
          DtLinMin(ivar) = minval(DtLin(ivar,1:iidx))
        end do
        
        ! Now, as the differences are transformed, we can limit every component on its own
        do ivar = 1, nvar
          DltIi(ivar) = max(min(DtLinMax(ivar),DtIi(ivar)),DtLinMin(ivar))
        end do
        
        
        
        ! Now we can trafo back
        DlIi = matmul(DR,DltIi)
        
        ! Calculate the correction factor
        ! for this element, for this edge, for this dimension (take min of all dimensions)
        do ivar = 1, nvar
          if (abs(DIi(ivar))<SYS_EPSREAL) then
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), 1.0_dp)
            ! That's the same as: Do nothing
          else
            ! This is the one following the principles
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))
            
            ! This one is less limiting
            !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), min(abs(DlIi(ivar)/DIi(ivar)),1.0_dp ))
            
            ! This is the one with the max of the dimensional splitting
            Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))
            
            ! This is the least limiting one
            !Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt), min(abs(DlIi(ivar)/DIi(ivar)),1.0_dp ) )
            
          end if
          
          ! No limiting at boundary
          if (iidx<3) then
            Dalphaei(ivar,ivt) = 1.0_dp
          end if
          
        end do
      
      end do ! idim
      end if
      
    end do ! ivt
    
    
    
    select case (ilim)
    case (1)
    
      ! Get minimum of all correction factors of all vertices on this element
      do ivar = 1, nvar
        Dalpha(ivar, iel) = minval(Dalphaei(ivar,1:NVE))
      end do
      
    case (2)
      
      ! Get minimum of all correction factors of all vertices on this element
      do ivar = 1, nvar
        Dalpha(ivar, iel) = min(Dalpha(ivar, iel),minval(Dalphaei(ivar,1:NVE)))
        p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(4:6,iel)-1) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(4:6,iel)-1)*Dalpha(ivar, iel)
      end do
      
    case (3)
    
      ! Get minimum of all correction factors of all vertices on this element
      do ivar = 1, nvar
        Dalpha(ivar, iel) = max(Dalpha(ivar, iel),minval(Dalphaei(ivar,1:NVE)))
        p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3,iel)-1) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3,iel)-1)*Dalpha(ivar, iel)
      end do
      
    end select
    
  end do !iel
  
  
!  ! *** Now limit the solution ***
!  do iel = 1, NEL
!  
!  
!  select case (ilim)
!  case (1)
!  
!    ! Do nothing      
!    
!  case (2)
!      
!    ! Multiply the quadratic part of the solution vector with the correction factor
!    do ivar = 1, nvar
!      !p_DoutputData(ivar)%p_Ddata(IdofGlob(4:6)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(4:6))*Dalpha(ivar, iel)
!      p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(4:6,iel)-1) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(4:6,iel)-1)*Dalpha(ivar, iel)
!    end do
!
!      
!  case (3)
!
!   ! Multiply the linear part of the solution vector with the correction factor
!    do ivar = 1, nvar
!      !p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3)) = p_DoutputData(ivar)%p_Ddata(IdofGlob(2:3))*Dalpha(ivar, iel)
!      p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3,iel)-1) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3,iel)-1)*Dalpha(ivar, iel)
!    end do
!    
!
!  end select
!  
!    end do ! iel
!  
!  
  end do ! ilim

  
  !deallocate(p_DoutputData)
  deallocate(DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi, DQchar)
  deallocate(DLin, DtLin, Dalphaei)
  deallocate(IdofGlob)
  deallocate(DLinPartValues)
  if (celement == EL_DG_T2_2D) deallocate(DDerivativeQuadraticValues)
  
  end subroutine
  
  
  
  
  
  
  !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_evaluateLinearPart_mult (rvectorBlock, IelIdx, IdofGlob, Dvalues)

!<description>

  ! Evaluates the linear part of an dg_T1 or dg_T2 element in the elements listed in
  ! IelIdx in the midpoints and the corners

!</description>

!<input>
type(t_vectorBlock), intent(in) :: rvectorBlock
integer, dimension(:), intent(in) :: IelIdx
integer, dimension(:,:), intent(in) :: IdofGlob
!</input>

!<output>
real(dp), dimension(:,:,:) :: Dvalues
!</output>
  
!</subroutine>

  ! Local variables
  real(DP), dimension(:), pointer :: p_Ddata
  integer :: iel, nvar, ivar
  
  ! Get pointer to the data of the (solution) vector
  call lsysbl_getbase_double (rvectorBlock, p_Ddata)
  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
  
  do iel = 1, size(IelIdx,1)
    do ivar = 1, nvar
      ! Get value in midpoint
      Dvalues(ivar,1,iel) = &
              p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1,iel)-1)
      ! Get value in first local corner
      Dvalues(ivar,2,iel) = &
              p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1,iel)-1) &
            - p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2,iel)-1) &
            - p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3,iel)-1)
      ! Get value in second local corner
      Dvalues(ivar,3,iel) = &
              p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1,iel)-1) &
            + p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2,iel)-1) &
            - p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3,iel)-1)
      ! Get value in third local corner
      Dvalues(ivar,4,iel) = &
              p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1,iel)-1) &
            + p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2,iel)-1) &
            + p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3,iel)-1)
      ! Get value in fourth local corner
      Dvalues(ivar,5,iel) = &
              p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1,iel)-1) &
            - p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2,iel)-1) &
            + p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3,iel)-1)
    end do
  end do
  
  end subroutine
  
  
  !****************************************************************************
  
!<subroutine>  
  
  subroutine dg_evaluateDerivativeQuadratic_mult (rvectorBlock, IelIdx, Dvalues, raddTriaData)

!<description>

  ! Evaluates the x- and y-derivatives of a dg_T2 element in the elements listed in
  ! IelIdx in the midpoints and the corners

!</description>

!<input>
type(t_vectorBlock), intent(in) :: rvectorBlock
type(t_additionalTriaData), intent(in) :: raddTriaData
integer, dimension(:), intent(in) :: IelIdx
!</input>

!<output>
! Dimension (nvar,# points per elements,# derivatives, nel)
real(dp), dimension(:,:,:,:) :: Dvalues
!</output>
  
!</subroutine>

  ! Local variables
  integer :: nvar, ivar, iel, ivt, NVE, iglobVtNumber
  integer(I32) :: celement
  real(dp), dimension(:,:,:), allocatable :: DpointsRef, Dpoints
  ! The underlying triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! The underlying spatial discretisation
  type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
  
  ! Pointers to some data from the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  real(dp), dimension(:,:), pointer :: p_DvertexCoords

  
  ! Get number of variables of the system
  nvar = rvectorBlock%nblocks
    
  ! Get pointers for quicker access
  p_rspatialDiscr => rvectorBlock%RvectorBlock(1)%p_rspatialDiscr
  p_rtriangulation => p_rspatialDiscr%p_rtriangulation
    
  ! What is the current element type?
  celement = p_rspatialDiscr%RelementDistr(1)%celement
  
  ! Number of vertices per element
  NVE = elem_igetNVE(celement)
  
  ! Get pointers to the data form the triangulation
  call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
  call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                               p_DvertexCoords)

  ! Build arrays with the coordinates, where to evaluate
  allocate(Dpoints(NDIM2D, 5, size(IelIdx,1)), DpointsRef(NDIM2D, 5, size(IelIdx,1)))
  
  do iel = 1, size(IelIdx,1)
    ! Save the midpoints
    Dpoints(1,1,iel) = raddTriaData%p_DmidPoints(1,IelIdx(iel))
    Dpoints(2,1,iel) = raddTriaData%p_DmidPoints(2,IelIdx(iel))
    
    ! Now calculate the limiting factor for every vertex on our element
    do ivt = 1, NVE
      ! Get global vertex number of our local vertex
      iglobVtNumber = p_IverticesAtElement(ivt,iel)
      
      Dpoints(1:2,ivt+1,iel) = p_DvertexCoords(1:2,iglobVtNumber)
    end do ! ivt
    
    ! Set coordinates on the reference element
    DpointsRef(1,1,iel) = 0.0_dp
    DpointsRef(2,1,iel) = 0.0_dp
    
    DpointsRef(1,2,iel) = -1.0_dp
    DpointsRef(2,2,iel) = -1.0_dp
    
    DpointsRef(1,3,iel) = +1.0_dp
    DpointsRef(2,3,iel) = -1.0_dp
    
    DpointsRef(1,4,iel) = +1.0_dp
    DpointsRef(2,4,iel) = +1.0_dp
    
    DpointsRef(1,5,iel) = -1.0_dp
    DpointsRef(2,5,iel) = +1.0_dp
  end do ! iel
   
  ! Now evaluate in the points
  do ivar = 1, nvar
    call fevl_evaluate_sim1 (DER_DERIV_X, Dvalues(ivar,:,1,:), &
                             rvectorBlock%rvectorBlock(ivar), &
                             Dpoints, IelIdx, DpointsRef)
    call fevl_evaluate_sim1 (DER_DERIV_Y, Dvalues(ivar,:,2,:), &
                             rvectorBlock%rvectorBlock(ivar), &
                             Dpoints, IelIdx, DpointsRef)
  end do ! ivar
  
  deallocate(Dpoints,DpointsRef)
  
  end subroutine
  
  
  
  
  !****************************************************************************
  
!<subroutine>  
  
   subroutine saveSolutionData(rvector,sofile,ifilenumber)

!<description>

  ! Output a DG vector to gmv format

!</description>

!<input>
 
  ! The solution vector to output
  type(t_vectorScalar), intent(in) :: rvector
  
  ! Name of output file
  character (LEN=SYS_STRLEN), intent(in) :: sofile
  
  ! Filenumber
  integer, intent(in) :: ifilenumber
  
    
!</input>

!<output>
!</output>
  
!</subroutine>
  ! local variables
  
  integer :: i, iunit, ilength
  character (LEN=10) :: sfilenumber
  real(dp), dimension(:), pointer :: p_Ddata
  
  ! Get pointers to the data form the truangulation
  call lsyssc_getbase_double(rvector,p_Ddata)

  ! Get the length of the data array
  ilength = size(p_Ddata,1)  
  
  
  ! ************ WRITE TO FILE PHASE *******************
  
  iunit = sys_getFreeUnit()
  open(iunit, file=trim(sofile) // '.data')
 
  
  write(iunit,'(I10)') ilength

  do i=1, ilength
    write(iunit,'(E25.16E3)') p_Ddata(i)
  end do
  
  close(iunit)



end subroutine




!****************************************************************************
  
!<subroutine>  
  
   subroutine loadSolutionData(rvector,sofile)

!<description>

  ! Loads the DOFs of a scalar vector

!</description>

!<input>
 
  ! The solution vector to output
  type(t_vectorScalar), intent(inout) :: rvector
  
  ! Name of output file
  character (LEN=SYS_STRLEN), intent(in) :: sofile
      
!</input>

!<output>
!</output>
  
!</subroutine>
  ! local variables
  
  integer :: i, iunit, ilength
  character (LEN=10) :: sfilenumber
  real(dp), dimension(:), pointer :: p_Ddata
  
  ! Get pointers to the data form the truangulation
  call lsyssc_getbase_double(rvector,p_Ddata)

  ! Get the length of the data array
  ilength = size(p_Ddata,1)  
  
  
  ! ************ WRITE TO FILE PHASE *******************
  
  iunit = sys_getFreeUnit()
  open(iunit, file=trim(sofile) // '.data')
 
  
  read(iunit,'(I10)') ilength

  do i=1, ilength
    read(iunit,'(E25.16E3)') p_Ddata(i)
  end do
  
  close(iunit)



end subroutine
  
  
end module
