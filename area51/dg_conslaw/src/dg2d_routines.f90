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
  
  use linearsystemblock
  use linearsystemscalar


  use ucd


  use dg2d_callback


  implicit none

  type t_dpPointer
     ! Pointer to the double-valued matrix or vector data
     real(DP), dimension(:), pointer :: p_Ddata
  end type t_dpPointer

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
  end type t_additionalTriaData

  type t_profiler
     integer :: ntimer, icurrenttimer
     real(dp) :: dstarttime, dendtime, dlasttime
     real(dp), dimension(:), allocatable :: Dtimers
  end type t_profiler
  
  type t_pMultigrid
     integer :: inumDiscr
     type(t_BlockDiscretisation), pointer, dimension(:) :: p_rDiscrLev
  end type t_pMultigrid

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
                daux1 = domega1 * DfluxValues(ialbet,icubp,iel)
                daux2 = domega2 * DfluxValues(ialbet,icubp,iel) *(-1.0_dp)

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

  end subroutine linf_dg_assembleSubmeshVectorScalarEdge2d














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



  end subroutine linf_getLocalEdgeNumbers














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

  end subroutine dg2gmv





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

    real(DP), dimension(:), allocatable :: dnodeValues,dnodeValues1,dnodeValues2

    character (LEN=10) :: sfilenumber

    real(dp) :: xc,yc

    real(dp), dimension(2) :: Dvel


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

    allocate(dnodeValues(nnodes),dnodeValues1(nnodes),dnodeValues2(nnodes))

    do iel = 1, NEL
       ibaseN = (iel-1)*nnodesOnRef
       call fevl_evaluate_mult1 (DER_FUNC, dnodeValues(ibaseN+1:ibaseN+nnodesOnRef),&
            rvector, iel, drefCoords(1:2,:))


       !    
       !    
       !    call fevl_evaluate_mult1 (DER_DERIV_X, dnodeValues1(ibaseN+1:ibaseN+nnodesOnRef),&
       !                              rvector, iel, drefCoords(1:2,:))
       !    call fevl_evaluate_mult1 (DER_DERIV_Y, dnodeValues2(ibaseN+1:ibaseN+nnodesOnRef),&
       !                              rvector, iel, drefCoords(1:2,:))
       !                              
       !                              
       !                              
       !                              
       !                              
       !                              
       !                              
       !                              
       !   xc = &
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
       !     ! Steady circular convection
       !    Dvel(1)=yc
       !    Dvel(2)=1.0_DP-xc
       !    
       !    Dvel(1) = Dvel(1)/(sqrt(Dvel(1)*Dvel(1)+Dvel(2)*Dvel(2)+SYS_EPSREAL_DP))
       !    Dvel(2) = Dvel(2)/(sqrt(Dvel(1)*Dvel(1)+Dvel(2)*Dvel(2)+SYS_EPSREAL_DP))
       !                          
       !    dnodeValues(ibaseN+1:ibaseN+nnodesOnRef) = dnodeValues1(ibaseN+1:ibaseN+nnodesOnRef)*Dvel(1) + dnodeValues2(ibaseN+1:ibaseN+nnodesOnRef)*Dvel(2)
       !    


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

  end subroutine dg2vtk



  !****************************************************************************

  !<subroutine>  

  subroutine dg_linearLimiter (rvector,ralpha)

    !<description>

    ! Limits the linear part of a dg_T1 element vector.

    !</description>

    !<input>
    !</input>

    !<inputoutput>

    ! A vector to limit
    type(t_vectorScalar), intent(inout) :: rvector  

    ! The limiting factors
    type(t_vectorScalar), intent(inout) :: ralpha

    !</inputoutput>

    !</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata, p_rAlpha_Ddata
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

    real(dp) :: dui, dui2, ddu, dalpha, dalphatemp, duc

    integer, dimension(3) :: IdofGlob

    integer :: NVBD

    integer, dimension(:), pointer :: p_InodalProperty

    ! Get pointer to the solution data
    call lsyssc_getbase_double (rvector,p_Ddata)

    call lsyssc_getbase_double (ralpha,p_rAlpha_Ddata)

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
    call storage_getbase_int(p_rtriangulation%h_InodalProperty ,&
         p_InodalProperty)

    ! Set pointer to coordinate vector
    call storage_getbase_double2D(&
         p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Set pointer to vertices at element
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    NVT = p_rtriangulation%NVT

    allocate(duimax(NVT),duimin(NVT))

    duimax= -SYS_MAXREAL_DP
    duimin=  SYS_MAXREAL_DP

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

          if(p_InodalProperty(nvt)>0) then
             duimax(nvt) =  1000000.0_dp
             duimin(nvt) = -1000000.0_dp
          end if
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

             !if (abs(duimax(nvert)-duc)<abs(0.1*duc)) dalphatemp = 1.0_dp
          elseif (ddu < 0.0_dp) then
             dalphatemp = min(1.0_dp, (duimin(nvert)-duc)/ddu)

             !if (abs(duimin(nvert)-duc)<abs(0.1*duc)) dalphatemp = 1.0_dp
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

       ! Save limiting factor for output later
       p_rAlpha_Ddata(IdofGlob(1)) = dalpha
       p_rAlpha_Ddata(IdofGlob(2:3)) = 0.0_dp

    end do ! iel

    deallocate(duimax,duimin)


  end subroutine dg_linearLimiter




















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


    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
         './gmv/steady.vtk')

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


  end subroutine dg_proj2steady










  !****************************************************************************

  !<subroutine>  

  subroutine dg_quadraticLimiter (rvector,ralpha)

    !<description>

    ! Limits the linear part of a dg_T1 element vector.

    !</description>

    !<input>
    !</input>

    !<inputoutput>

    ! A vector to limit
    type(t_vectorScalar), intent(inout) :: rvector  

    ! The limiting factors
    type(t_vectorScalar), intent(inout) :: ralpha

    !</inputoutput>

    !</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata, p_rAlpha_Ddata
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

    integer :: iglobNeighNum

    integer, dimension(:), pointer :: p_InodalProperty 

    integer :: ilocalVtnumber, iglobVtNumber

    real(DP), dimension(:,:), allocatable :: Dalpha


    ! Get pointer to the solution data
    call lsyssc_getbase_double (rvector,p_Ddata)

    call lsyssc_getbase_double (ralpha,p_rAlpha_Ddata)

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
    call storage_getbase_int(p_rtriangulation%h_InodalProperty ,&
         p_InodalProperty)
    call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
         p_IelementsAtVertex) 
    call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
         p_IelementsAtVertexIdx)



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


!!! First calculate the bounds

       duimax= -SYS_MAXREAL_DP
       duimin=  SYS_MAXREAL_DP

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



!!!!!!!!!!!!!!!!!!!!!!!  NEU !!!!!!!!!!!!!!!!!!!!!!!

          if(ilim<3) then

             !    do ivert = 1, NVE
             !		
             !      nvert = p_IverticesAtElement(ivert, iel)
             !      
             !      ! The second point we want to evaluate the solution in, is in the corner of the mother element
             !      Dpoints(1,1+ivert) = p_DvertexCoords(1,nvert)
             !      Dpoints(2,1+ivert) = p_DvertexCoords(2,nvert)
             !      Ielements(1+ivert) = iel
             !	  end do
             !	  
             !    ! Evaluate the solution
             !    call fevl_evaluate (ideriv, Dvalues(1:5), rvector, Dpoints(1:2,1:5), &
             !                          Ielements(1:5))
             !    do ivt = 1, NVE
             !      nvt = p_IverticesAtElement(ivt,iel)
             !      duimax(nvt) = max(maxval(Dvalues(1:5)),duimax(nvt))
             !      duimin(nvt) = min(minval(Dvalues(1:5)),duimin(nvt))
             !    end do

          else

             !      call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
             !      duc = p_Ddata(IdofGlob(1))+abs(p_Ddata(IdofGlob(2)))+abs(p_Ddata(IdofGlob(3)))
             !      
             !    do ivt = 1, NVE
             !      nvt = p_IverticesAtElement(ivt,iel)
             !      duimax(nvt) = max(p_Ddata(IdofGlob(1))+abs(p_Ddata(IdofGlob(2)))+abs(p_Ddata(IdofGlob(3))),duimax(nvt))
             !      duimin(nvt) = min(p_Ddata(IdofGlob(1))-abs(p_Ddata(IdofGlob(2)))-abs(p_Ddata(IdofGlob(3))),duimin(nvt))
             !    end do

          end if

!!!!!!!!!!!!!!!!!!!!!!!  NEU !!!!!!!!!!!!!!!!!!!!!!!



       end do





       !  do ivt = 1, nvt
       !  iidx = 0
       !    do ineighbour = p_IelementsAtVertexIdx(ivt), p_IelementsAtVertexIdx(ivt+1)-1
       !    
       !      if(ilim<3) then
       !        iglobNeighNum = p_IelementsAtVertex(ineighbour)
       !      
       !        iidx = iidx + 1
       !      
       !        Dpoints(1,iidx) = p_DvertexCoords(1,iglobNeighNum)
       !        Dpoints(2,iidx) = p_DvertexCoords(2,iglobNeighNum)
       !        Ielements(iidx) = iglobNeighNum 
       !        
       !      else
       !        iidx = iidx + 1
       !        
       !        do ilocalVtNumber = 1, 4
       !          if (p_IverticesAtElement(ilocalVtnumber,iglobNeighNum )==ivt) exit
       !        end do
       !        
       !        call dof_locGlobMapping(p_rspatialDiscr, iglobNeighNum, IdofGlob)
       !        dui = p_Ddata(IdofGlob(1))
       !        
       !        select case (ilocalVtnumber)
       !        case(1)
       !        dui = dui - p_Ddata(IdofGlob(2)) - p_Ddata(IdofGlob(3))
       !        case(2)
       !        dui = dui + p_Ddata(IdofGlob(2)) - p_Ddata(IdofGlob(3))
       !        case(3)
       !        dui = dui + p_Ddata(IdofGlob(2)) + p_Ddata(IdofGlob(3))
       !        case(4)
       !        dui = dui - p_Ddata(IdofGlob(2)) + p_Ddata(IdofGlob(3))
       !        end select
       !        
       !        Dvalues(iidx) = dui
       !        
       !      end if
       !    end do
       !    
       !    ! Evaluate the solution
       !    if(ilim<3) call fevl_evaluate (ideriv, Dvalues(1:iidx), rvector, Dpoints(1:2,1:iidx), &
       !                                   Ielements(1:iidx))
       !      
       !      
       !    duimax(ivt) = max(maxval(Dvalues(1:iidx)),duimax(ivt))
       !    duimin(ivt) = min(minval(Dvalues(1:iidx)),duimin(ivt))
       !        
       !  
       !  end do




       ! do iel = 1, NEL
       ! 
       ! iidx = 0
       ! 
       !  do ivt = 1, 4
       !  
       !    nvt = p_IverticesAtElement(ivt,iel)
       !  
       !    do ineighbour = p_IelementsAtVertexIdx(nvt), p_IelementsAtVertexIdx(nvt+1)-1
       !    
       !    iglobNeighNum = p_IelementsAtVertex(ineighbour)
       !    iidx = iidx + 1
       !    
       !      if(ilim<3) then
       !        
       !      
       !        
       !        ! Get midpoint of the element
       !        xc = &
       !         (p_DvertexCoords(1,p_IverticesAtElement(1,iglobNeighNum))+&
       !          p_DvertexCoords(1,p_IverticesAtElement(2,iglobNeighNum))+&
       !          p_DvertexCoords(1,p_IverticesAtElement(3,iglobNeighNum))+&
       !          p_DvertexCoords(1,p_IverticesAtElement(4,iglobNeighNum)))/4.0_dp
       !
       !        yc = &
       !         (p_DvertexCoords(2,p_IverticesAtElement(1,iglobNeighNum))+&
       !          p_DvertexCoords(2,p_IverticesAtElement(2,iglobNeighNum))+&
       !          p_DvertexCoords(2,p_IverticesAtElement(3,iglobNeighNum))+&
       !          p_DvertexCoords(2,p_IverticesAtElement(4,iglobNeighNum)))/4.0_dp
       !        
       !        Dpoints(1,iidx) = xc
       !        Dpoints(2,iidx) = yc
       !        Ielements(iidx) = iglobNeighNum 
       !  
       !      else
       !      
       !        call dof_locGlobMapping(p_rspatialDiscr, iglobNeighNum, IdofGlob)
       !        Dvalues(iidx) = p_Ddata(IdofGlob(1))
       !      
       !      end if
       !    
       !    end do
       !  end do
       !    
       !  ! Evaluate the solution
       !  if(ilim<3) call fevl_evaluate (ideriv, Dvalues(1:iidx), rvector, Dpoints(1:2,1:iidx), &
       !                                 Ielements(1:iidx))
       !                                 
       !  do ivt = 1, 4
       !    iglobVtNumber = p_IverticesAtElement(ivt,iel)
       !    duimax(iglobVtNumber) = max(maxval(Dvalues(1:iidx)),duimax(iglobVtNumber))
       !    duimin(iglobVtNumber) = min(minval(Dvalues(1:iidx)),duimin(iglobVtNumber))
       !  
       !  end do
       !    
       !end do    









!!! Start limiting



       do iel = 1, NEL


          ! No limiting of elements at boundary
          if ((p_InodalProperty(p_IverticesAtElement(1, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(2, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(3, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(4, iel))>0)) cycle


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

                !        ! Extremumfix
                !        if (duimax(nvert)-duc<SYS_EPSREAL_DP) dalphatemp = 1.0_dp
                !        if ((ilim>0).and.(duimax(nvert)-duc<abs(0.0001_dp*duc))) dalphatemp = 1.0_dp

             elseif (ddu < 0.0_dp) then
                dalphatemp = min(1.0_dp, (duimin(nvert)-duc)/ddu)

                !        ! Extremumfix
                !        if (duimin(nvert)-duc>-SYS_EPSREAL_DP) dalphatemp = 1.0_dp
                !        if ((ilim>0).and.(duimin(nvert)-duc>-abs(0.0001_dp*duc))) dalphatemp = 1.0_dp

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


             ! Limiting like Kuzmin did it
             Dalpha(1,iel) = min(Dalpha(1,iel),Dalpha(2,iel))

             ! Get global DOFs of the element
             call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)


             ! Multiply the quadratic part of the solution vector with the correction factor
             p_Ddata(IdofGlob(4:6)) = p_Ddata(IdofGlob(4:6))*Dalpha(1,iel)

             p_rAlpha_Ddata(IdofGlob(1)) = Dalpha(1,iel)
             p_rAlpha_Ddata(IdofGlob(2:3)) = 0.0_dp


             !        ! Limiting using different limiting factors for the second derivatives
             !        ! Get global DOFs of the element
             !        call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
             !    
             !    
             !        ! Multiply the linear part of the solution vector with the correction factor
             !        p_Ddata(IdofGlob(4)) = p_Ddata(IdofGlob(4))*Dalpha(1,iel)
             !        p_Ddata(IdofGlob(5)) = p_Ddata(IdofGlob(5))*min(Dalpha(1,iel),Dalpha(2,iel))
             !        p_Ddata(IdofGlob(6)) = p_Ddata(IdofGlob(6))*Dalpha(2,iel)
             !        
             !        
             !        Dalpha(1,iel) = min(Dalpha(1,iel),Dalpha(2,iel))
             !        
             !        p_rAlpha_Ddata(IdofGlob(1)) = Dalpha(1,iel)
             !        p_rAlpha_Ddata(IdofGlob(2:3)) = 0.0_dp


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


  end subroutine dg_quadraticLimiter





















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


  end subroutine addTriaData



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

  end subroutine releaseTriaData










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
    integer :: ielementDistr, NMT, NVT, iedge, i

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
          call dg_initAssembly_reverseCubPoints(rvectorAssembly(2), rform,&
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
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D_1, Dxi1D_2
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

    ! Copy the second component and replace 0s by 1s
    IelementList(3,size(IedgeList))=IelementList(2,size(IedgeList))
    do iel = 1,size(IedgeList)
       IelementList(3,iel)=max(IelementList(2,iel),1)
    end do

    ! Allocate daux
    allocate(daux(nvar,2))

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

    ! Allocate space for the flux variables DIM(nvar,ialbet,ncubp,elementsPerBlock)
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
          Dxi1D_1(icubp,k) = rlocalVectorAssembly(1)%p_DcubPtsRef(k,icubp)
          Dxi1D_2(icubp,k) = rlocalVectorAssembly(2)%p_DcubPtsRef(k,icubp)
       end do
    end do

    ! Allocate memory for the cubature points in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,2,rlocalVectorAssembly(1)%nelementsPerBlock))

    ! Allocate memory for the coordinates of the reference points
    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalVectorAssembly(1)%nelementsPerBlock,2))

    ! Allocate DlocalData, which will later 
    allocate(DlocalData(nvar,2,indof))

    !    ! Allocate the space for the pointer to the Data of the different blocks of the output vector
    !    allocate(p_DoutputData(nvar))

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
               ncubp, Dxi1D_1, Dxi2D(:,:,1,iel))
          call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(2,Iedgelist(IELset+iel-1)), &
               ncubp, Dxi1D_2, Dxi2D(:,:,2,iel))
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
             domega2 = dlen * rlocalVectorAssembly(2)%p_Domega(icubp)


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
                !daux1 = domega1 * rlocalVectorAssembly(1)%p_Dcoefficients(ialbet,icubp,iel)
                !daux2 = domega2 * rlocalVectorAssembly(2)%p_Dcoefficients(ialbet,icubp,iel) *(-1.0_dp)
                Daux(:,1) = domega1 * DfluxValues(:,ialbet,icubp,iel)
                Daux(:,2) = domega2 * DfluxValues(:,ialbet,icubp,iel) *(-1.0_dp)

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
    deallocate(DfluxValues,daux,DlocalData)!,p_DoutputData)

  end subroutine linf_dg_assembleSubmeshVectorBlockEdge2d










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

    duimax= -SYS_MAXREAL_DP
    duimin=  SYS_MAXREAL_DP

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


  end subroutine dg_linearLimiterBlockIndicatorVar




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


       duimax= -SYS_MAXREAL_DP
       duimin=  SYS_MAXREAL_DP

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


  end subroutine dg_quadraticLimiterBlockIndicatorVar







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
                   if (abs(DIi(ivar))<SYS_EPSREAL_DP) then
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


  end subroutine dg_linearLimiterBlockCharVar






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


       !  duimax= -SYS_MAXREAL_DP
       !  duimin=  SYS_MAXREAL_DP
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


  end subroutine dg_quadraticLimiterBlockIndicatorVar_2

































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
                      if (abs(DIi(ivar))<SYS_EPSREAL_DP) then
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


  end subroutine dg_quadraticLimiterBlockCharVar




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

  end subroutine profiler_init


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
       write(*,*) i,  (rprofiler%Dtimers(i)/(dalltime+SYS_EPSREAL_DP)*100),'%'
    end do
    write(*,*) '*********************************************************************'
    write(*,*) ''

    rprofiler%dstarttime    = 0.0_dp
    rprofiler%dendtime      = 0.0_dp
    rprofiler%dlasttime     = 0.0_dp
    rprofiler%ntimer        = -1
    rprofiler%icurrenttimer = -1
    deallocate(rprofiler%Dtimers)

  end subroutine profiler_release


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

  end subroutine profiler_measure



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

    dt = SYS_MAXREAL_DP

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

  end subroutine getDtByCfl







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
                   if (dquo<SYS_EPSREAL_DP) then
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
                   if (dquo<SYS_EPSREAL_DP) then
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
                   if (abs(DIi(ivar))<SYS_EPSREAL_DP) then
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


  end subroutine dg_linearLimiterBlockCharVar_mixedJacobian




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

  end subroutine dg_evaluateLinearPart
















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
                      if (dquo<SYS_EPSREAL_DP) then
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
                      if (dquo<SYS_EPSREAL_DP) then
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
                      if (abs(DIi(ivar))<SYS_EPSREAL_DP) then
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


  end subroutine dg_quadraticLimiterBlockCharVar_mixedJacobian





































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

    integer, dimension(:), pointer :: p_InodalProperty 

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
    call storage_getbase_int(p_rtriangulation%h_InodalProperty ,&
         p_InodalProperty)                              

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
    !        if (dquo<SYS_EPSREAL_DP) then
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
    !        if (dquo<SYS_EPSREAL_DP) then
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

          !    ! No limiting of elements at boundary
          !    if ((p_InodalProperty(p_IverticesAtElement(1, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(2, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(3, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(4, iel))>0)) cycle

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


                   if (idim==1) then
                      DL = buildMixedL2(DQchar,1.0_dp,0.0_dp)
                      DR = buildMixedR2(DQchar,1.0_dp,0.0_dp)
                   elseif (idim==2) then
                      DL = buildMixedL2(DQchar,0.0_dp,1.0_dp)
                      DR = buildMixedR2(DQchar,0.0_dp,1.0_dp)

                   elseif (idim==3) then
                      DL = buildMixedL2(DQchar,-1.0_dp,0.0_dp)
                      DR = buildMixedR2(DQchar,-1.0_dp,0.0_dp)

                   elseif (idim==4) then
                      DL = buildMixedL2(DQchar,0.0_dp,-1.0_dp)
                      DR = buildMixedR2(DQchar,0.0_dp,-1.0_dp)

                   else if (idim==5) then
                      da = DQchar(2)/DQchar(1)
                      db = DQchar(3)/DQchar(1)
                      dquo = da*da+db*db
                      if (dquo<SYS_EPSREAL_DP) then
                         DL = buildMixedL2(DQchar,1.0_dp,0.0_dp)
                         DR = buildMixedR2(DQchar,1.0_dp,0.0_dp)
                      else
                         da = da/dquo
                         db = db/dquo
                         DL = buildMixedL2(DQchar,da,db)
                         DR = buildMixedR2(DQchar,da,db)
                      end if

                   else if (idim==6) then
                      da = DQchar(2)/DQchar(1)
                      db = DQchar(3)/DQchar(1)
                      dquo = da*da+db*db
                      if (dquo<SYS_EPSREAL_DP) then
                         DL = buildMixedL2(DQchar,0.0_dp,1.0_dp)
                         DR = buildMixedR2(DQchar,0.0_dp,1.0_dp)
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
                      if (abs(DIi(ivar))<10.0_dp*SYS_EPSREAL_DP) then
                         Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), 1.0_dp)
                      else
                         !            !This is the one following the principles (than even adjust initialisation of Dalphaei!!!)
                         !            Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))

                         ! This one is less limiting
                         !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt), min(abs(DlIi(ivar)/DIi(ivar)),1.0_dp ))

                         ! This is the one with the max of the dimensional splitting
                         Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt), max(0.0_dp, min(DlIi(ivar)/DIi(ivar),1.0_dp) ))

                         ! This is the least limiting one
                         !Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt), min(abs(DlIi(ivar)/DIi(ivar)),1.0_dp ) )

                      end if

                      !          ! No limiting at boundary
                      !          if (iidx<3) then
                      !            Dalphaei(ivar,ivt) = 1.0_dp
                      !          end if

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

  end subroutine dg_kuzminLimiterBlockCharVar_mixedJacobian






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

  end subroutine dg_evaluateLinearPart_mult


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

  end subroutine dg_evaluateDerivativeQuadratic_mult




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



  end subroutine saveSolutionData




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

    ! Get pointers to the data form the triangulation
    call lsyssc_getbase_double(rvector,p_Ddata)

    ! Get the length of the data array
    ilength = size(p_Ddata,1)

    ! write(*,*) ilength


    ! ************ WRITE TO FILE PHASE *******************

    iunit = sys_getFreeUnit()
    open(iunit, file=trim(sofile) // '.data')


    read(iunit,'(I10)') ilength

    ! write(*,*) ilength

    do i=1, ilength
       read(iunit,'(E25.16E3)') p_Ddata(i)
    end do

    close(iunit)



  end subroutine loadSolutionData



















  subroutine calc_error(rvector, derror, raddtriadata)

    use transformation

    type(t_vectorScalar), intent(in) :: rvector
    real(dp), intent(out) :: derror
    type(t_additionalTriaData), intent(in):: raddTriaData


    integer :: iel, ipoint






    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: indof, NEL, NVE, ivt, NVT

    ! The underlying triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! The underlying spatial discretisation
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr


    ! The list of elements, in which these points can be found
    integer, dimension(17) :: Ielements


    ! Pointers to some data from the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx, p_IelementsAtVertex
    real(dp), dimension(:,:), pointer :: p_DvertexCoords

    integer, dimension(:), pointer :: p_IverticesAtBoundary

    real(dp) :: xc,yc
    integer :: iidx, nvert, ivert, ineighbour, ineighElm

    real(dp) :: dui, ddu, duc

    integer, dimension(:,:), allocatable :: IdofGlob


    integer :: nvar, ivar, idim


    integer :: nDOFloc


    integer(I32) :: celement

    integer :: npoints

    integer, dimension(:), pointer :: p_IelIdx

    ! Values of the linear part of the solution, Dimension (nvar,# points per element (midpoint+corners),NEL)
    real(dp), dimension(:,:,:), allocatable :: DLinPartValues

    ! Dimension (nvar,# points per element (midpoint+corners),# derivatives (2, x- and y-),NEL)
    real(dp), dimension(:,:,:,:), allocatable :: DDerivativeQuadraticValues

    integer :: ccubtype, ctrafotype, ncubp,i ,iunit,icubp

    real(dp), dimension(:), allocatable :: Domega

    real(dp), dimension(:,:), allocatable :: Dcoords

    real(dp), dimension(:,:), pointer :: p_DcubPtsRef

    real(dp), dimension(8) :: DjacPrep

    real(dp) :: r,n,ddetj,dxreal,dyreal,drefsol

    real(dp), dimension(1)::Dvalues
    real(dp), dimension(2,1)::Dpoints

    real(dp), dimension(10001) :: Dreference


    ! Get pointers for quicker access
    p_rspatialDiscr => rvector%p_rspatialDiscr
    p_rtriangulation => p_rspatialDiscr%p_rtriangulation

    ! Allocate the space for the pointer to the Data of the different blocks of the output vector
    !  allocate(p_DoutputData(nvar))
    !    
    !  do ivar = 1, nvar
    !    call lsyssc_getbase_double(rvectorBlock%RvectorBlock(ivar),p_DoutputData(ivar)%p_Ddata)
    !  end do

    ! Get pointer to the data of the (solution) vector
    !  call lsysbl_getbase_double (rvectorBlock, p_Ddata)

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


    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(celement)

    ccubtype = CUB_G5x5

    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(ccubtype)

    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType), ncubp))
    allocate(Dcoords(trafo_igetReferenceDimension(ctrafoType),nve))

    ! Get the cubature formula
    call cub_getCubature(ccubtype, p_DcubPtsRef, Domega)

    ! Get global DOFs for these (all) elements
    call dof_locGlobMapping_mult(p_rspatialDiscr, p_IelIdx, IdofGlob)

    ! Test, if we have quadratic elements
    !if (celement == EL_DG_T2_2D) then





    ! Circular dambreak

    iunit = sys_getFreeUnit()

    open(iunit, file='h')

    do i = 1, 10000
       read(iunit,*) Dreference(i)
    end do

    close(iunit)






    derror = 0.0_dp




    do iel = 1, nel
       write(*,*)real(real(iel)/real(nel))

       call trafo_getCoords (ctrafoType,p_rtriangulation,iel,Dcoords)

       call trafo_prepJac_quad2D(Dcoords, DjacPrep)

       ddetj = raddTriaData%p_Ddxdy(1,iel)*raddTriaData%p_Ddxdy(2,iel)

       do icubp = 1, ncubp

          call trafo_calcRealCoords (DjacPrep,p_DcubPtsRef(1,icubp),p_DcubPtsRef(2,icubp),dxreal,dyreal)


          r = sqrt(dxreal*dxreal+dyreal*dyreal)
          n = 1+1000*r

          drefsol =(1.0_dp-(n-real(int(n))))* Dreference(int(n)) +(n-real(int(n)))* Dreference(int(n)+1)

          Dpoints(1,1) =dxreal
          Dpoints(2,1) =dyreal

          call fevl_evaluate (DER_FUNC, Dvalues,&
               rvector, Dpoints)




          derror = derror + ddetj*Domega(icubp)*abs(drefsol-Dvalues(1))

       end do !icubp



    end do ! iel



    deallocate(Dcoords)
    deallocate(p_DcubPtsRef)
    deallocate(Domega)





  end subroutine calc_error








  !****************************************************************************

  !<subroutine>  

  subroutine dg_realKuzmin (rvectorBlock, raddTriaData)

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

    real(dp), dimension(:), allocatable :: DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi, DQchar, DWc

    real(dp), dimension(:,:), allocatable :: DLin, DtLin, Dalphaei, DL, DR, Dalpha, DL2, DR2, DlinearGradient, Dquadraticgradient, DmAlpha

    ! Array of pointers to the data of the blockvector to limit
    type(t_dpPointer), dimension(:), allocatable :: p_DoutputData

    integer :: iglobVtNumber, iglobNeighNum

    integer :: ilim, ideriv, ilimstart

    integer :: nDOFloc

    real(dp) :: da, db, dquo

    integer, parameter :: nelemSim = 1000

    integer(I32) :: celement

    integer :: npoints, ipoint

    integer :: i,j

    integer, dimension(:), pointer :: p_IelIdx

    ! Values of the linear part of the solution, Dimension (nvar,# points per element (midpoint+corners),NEL)
    real(dp), dimension(:,:,:), allocatable :: DLinPartValues

    ! Dimension (nvar,# points per element (midpoint+corners),# derivatives (2, x- and y-),NEL)
    real(dp), dimension(:,:,:,:), allocatable :: DDerivativeQuadraticValues

    integer, parameter :: nmaxneighbours = 10

    integer, dimension(:), pointer :: p_InodalProperty 

    real(dp), dimension(3,3) :: DTemp1, DTemp2

    real(dp) :: dWstar

    real(dp) :: dl2x, dl2y

    real(dp), dimension(5) :: DQcharExt

    real(dp) :: dlf, dh




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
    call storage_getbase_int(p_rtriangulation%h_InodalProperty ,&
         p_InodalProperty)                              

    ! Allocate the space for solution differences, transformed solution differences,
    ! limited transformed solution differences, limited backtransformed solution differences
    ! and limiting factors
    allocate(DVec(nvar), DVei(nvar), DIi(nvar), DtIi(nvar), DtLinMax(nvar), DtLinMin(nvar), DltIi(nvar), DlIi(nvar),DQchar(nvar),DWc(nvar))
    !allocate(DLin(nvar,NVE-1), DtLin(nvar,NVE-1), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL))
    allocate(DLin(nvar,10), DtLin(nvar,10), Dalphaei(nvar,NVE), DL(nvar,nvar), DR(nvar,nvar), Dalpha(nvar, NEL), DL2(nvar,nvar), DR2(nvar,nvar),DlinearGradient(nvar,2), DquadraticGradient(nvar,3), DmAlpha(nvar,nvar))

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

    !  do iel = 1,nel
    !    if (p_IelIdx(iel) .ne. iel) write(*,*) 'Warning'
    !  end do

    select case (celement)
    case (EL_DG_T2_2D) 
       ilimstart = 1
       Dalpha = 1.0_DP
    case (EL_DG_T1_2D) 
       ilimstart = 3
       Dalpha = 0.0_DP
    end select




    do ilim = ilimstart, 3

       do iel = 1, NEL

          ! No limiting of elements at boundary
          !if ((p_InodalProperty(p_IverticesAtElement(1, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(2, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(3, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(4, iel))>0)) cycle

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

                iidx = iidx + 1

                ! Get values in the center of the element for all variables


                if (ilim<3) then
                   DLin(:,iidx) = DDerivativeQuadraticValues(:,1,ilim,iglobNeighNum)
                else
                   DLin(:,iidx) = DLinPartValues(:,1,iglobNeighNum)
                end if

                ! Calculate solution difference
                DLin(:,iidx) = DLin(:,iidx)

                !end if
             end do



             if (iidx.ne.0) then
                ! Dimensional splitting
                do idim = -1,-1
                   ! Now we need the trafo matrices


                   if (idim==1) then
                      DL = buildMixedL2(DQchar,1.0_dp,0.0_dp)
                      DR = buildMixedR2(DQchar,1.0_dp,0.0_dp)
                   elseif (idim==-1) then  
                      DL = 1.0_dp
                      DR = 1.0_dp
                   elseif (idim==0) then  
                      da = 1.0_dp/sqrt(2.0_dp)
                      db = 1.0_dp/sqrt(2.0_dp)
                      DL = buildMixedL2(DQchar,da,db)
                      DR = buildMixedR2(DQchar,da,db)
                   elseif (idim==2) then
                      DL = buildMixedL2(DQchar,0.0_dp,1.0_dp)
                      DR = buildMixedR2(DQchar,0.0_dp,1.0_dp)
                   elseif (idim==3) then
                      DL = buildMixedL2(DQchar,-1.0_dp,0.0_dp)
                      DR = buildMixedR2(DQchar,-1.0_dp,0.0_dp)

                   elseif (idim==4) then
                      DL = buildMixedL2(DQchar,0.0_dp,-1.0_dp)
                      DR = buildMixedR2(DQchar,0.0_dp,-1.0_dp)

                   else if (idim==5) then
                      da = DQchar(2)
                      db = DQchar(3)
                      dquo = sqrt(da*da+db*db)
                      if (dquo<SYS_EPSREAL_DP) then
                         DL = buildMixedL2(DQchar,1.0_dp,0.0_dp)
                         DR = buildMixedR2(DQchar,1.0_dp,0.0_dp)
                      else
                         da = da/dquo
                         db = db/dquo
                         DL = buildMixedL2(DQchar,da,db)
                         DR = buildMixedR2(DQchar,da,db)
                      end if

                   else if (idim==6) then
                      da = DQchar(2)
                      db = DQchar(3)
                      dquo = sqrt(da*da+db*db)
                      if (dquo<SYS_EPSREAL_DP) then
                         DL = buildMixedL2(DQchar,0.0_dp,1.0_dp)
                         DR = buildMixedR2(DQchar,0.0_dp,1.0_dp)
                      else
                         da = da/dquo
                         db = db/dquo
                         DL = buildMixedL2(DQchar,-db,da)
                         DR = buildMixedR2(DQchar,-db,da)
                      end if

                   else if (idim==7) then
                      da = DQchar(2)
                      db = DQchar(3)
                      if (abs(da)>abs(db)) then
                         DL = buildMixedL2(DQchar,1.0_dp,0.0_dp)
                         DR = buildMixedR2(DQchar,1.0_dp,0.0_dp)
                      else
                         DL = buildMixedL2(DQchar,0.0_dp,1.0_dp)
                         DR = buildMixedR2(DQchar,0.0_dp,1.0_dp)
                      end if

                   else if (idim==8) then
                      da = DQchar(2)
                      db = DQchar(3)

                      if (abs(da)>abs(db)) then
                         if(da>0.0_dp)then
                            DL = buildMixedL2(DQchar,1.0_dp,0.0_dp)
                            DR = buildMixedR2(DQchar,1.0_dp,0.0_dp)
                         else
                            DL = buildMixedL2(DQchar,-1.0_dp,0.0_dp)
                            DR = buildMixedR2(DQchar,-1.0_dp,0.0_dp)
                         end if
                      else
                         if(db>0.0_dp)then
                            DL = buildMixedL2(DQchar,0.0_dp,1.0_dp)
                            DR = buildMixedR2(DQchar,0.0_dp,1.0_dp)
                         else
                            DL = buildMixedL2(DQchar,0.0_dp,-1.0_dp)
                            DR = buildMixedR2(DQchar,0.0_dp,-1.0_dp)
                         end if
                      end if

                   else if (idim==9) then
                      da = raddTriaData%p_DmidPoints(1,iel)-0.5_dp
                      db = raddTriaData%p_DmidPoints(2,iel)-0.5_dp

                      if (abs(da)>abs(db)) then
                         if(da>0.0_dp)then
                            DL = buildMixedL2(DQchar,1.0_dp,0.0_dp)
                            DR = buildMixedR2(DQchar,1.0_dp,0.0_dp)
                         else
                            DL = buildMixedL2(DQchar,-1.0_dp,0.0_dp)
                            DR = buildMixedR2(DQchar,-1.0_dp,0.0_dp)
                         end if
                      else
                         if(db>0.0_dp)then
                            DL = buildMixedL2(DQchar,0.0_dp,1.0_dp)
                            DR = buildMixedR2(DQchar,0.0_dp,1.0_dp)
                         else
                            DL = buildMixedL2(DQchar,0.0_dp,-1.0_dp)
                            DR = buildMixedR2(DQchar,0.0_dp,-1.0_dp)
                         end if
                      end if

                   else if (idim==10) then
                      da = DQchar(2)
                      db = DQchar(3)
                      dquo = sqrt(da*da+db*db)
                      DQcharext = Euler_transformVector(DQchar)

                      if (dquo<10.0_dp*SYS_EPSREAL_DP) then
                         DL = Euler_buildMixedLcfromRoe(DQcharext,1.0_dp,0.0_dp)
                         DR = Euler_buildMixedLcfromRoe(DQcharext,1.0_dp,0.0_dp)
                      else
                         da = da/dquo
                         db = db/dquo
                         DL = Euler_buildMixedLcfromRoe(DQcharext,da,db)
                         DR = Euler_buildMixedRcfromRoe(DQcharext,da,db)
                      end if

                   end if

                   ! Transform neighbouring limits
                   DtIi = matmul(DL,DIi)
                   do ineighbour = 1, iidx

                      !          if (idim==1) then
                      !            DL2 = buildMixedL2(Dlin(:,ineighbour),1.0_dp,0.0_dp)
                      !          
                      !          elseif (idim==2) then
                      !            DL2 = buildMixedL2(Dlin(:,ineighbour),0.0_dp,1.0_dp)
                      !          
                      !          end if

                      DtLin(:,ineighbour) = matmul(DL,Dlin(:,ineighbour))
                   end do

                   !        ! or
                   !        DtLin(:,1:iidx) = matmul(DL,Dlin(:,1:iidx))


                   ! Now, as the differences are transformed, we can limit every component on its own
                   DWc=matmul(DL,DVec)

                   do ivar = 1, nvar
                      if (DtIi(ivar) > 0.0_dp) then
                         dWstar = maxval(DtLin(ivar,1:iidx))
                         !Dalphaei(ivar,ivt) = min(  Dalphaei(ivar,ivt) ,min(1.0_dp,(dwstar-DWc(ivar))/(DtIi(ivar))))
                         
                         ! This is the one I take
                         Dalphaei(ivar,ivt) = max(  Dalphaei(ivar,ivt) ,min(1.0_dp,(dwstar-DWc(ivar))/(DtIi(ivar)+SYS_EPSREAL_DP)))
                      
                         !            ! Extremumfix
                         !            if (dwstar-DWc(ivar)<10.0*SYS_EPSREAL_DP) Dalphaei(ivar,ivt) = 1.0_dp
                         !            if ((ilim==3).and.(abs(dwstar-DWc(ivar))<abs(0.001*DWc(ivar)))) Dalphaei(ivar,ivt) = 1.0_dp


                         !            ! Smoothed limiter 1
                         !            dlf = (dwstar-DWc(ivar))/(DtIi(ivar)+SYS_EPSREAL_DP)
                         !            dlf = dlf*(dlf+2.0_dp)/(dlf*(dlf+1.0_dp)+2.0_dp)
                         !            Dalphaei(ivar,ivt) = max(  Dalphaei(ivar,ivt) , dlf )

                         !            ! Smoothed limiter 2
                         !            dlf = (dwstar-DWc(ivar))/(DtIi(ivar)+SYS_EPSREAL_DP)
                         !            dlf = min(1.0_dp,dlf*dlf)
                         !            Dalphaei(ivar,ivt) = max(  Dalphaei(ivar,ivt) , dlf )


                      elseif (DtIi(ivar) < 0.0_dp) then
                         dWstar = minval(DtLin(ivar,1:iidx))
                         !Dalphaei(ivar,ivt) = min(  Dalphaei(ivar,ivt) ,min(1.0_dp,(dwstar-DWc(ivar))/(DtIi(ivar))))
                         
                         ! This is the one I take
                         Dalphaei(ivar,ivt) = max(  Dalphaei(ivar,ivt) ,min(1.0_dp,(dwstar-DWc(ivar))/(DtIi(ivar)-SYS_EPSREAL_DP)))


                         !            ! Extremeumfix
                         !            if (dwstar-DWc(ivar)<10.0*SYS_EPSREAL_DP) Dalphaei(ivar,ivt) = 1.0_dp
                         !            if ((ilim==3).and.(abs(dwstar-DWc(ivar))<abs(0.001*DWc(ivar)))) Dalphaei(ivar,ivt) = 1.0_dp

                         !            ! Smoothed limiter 1
                         !            dlf = (dwstar-DWc(ivar))/(DtIi(ivar)-SYS_EPSREAL_DP)
                         !            dlf = dlf*(dlf+2.0_dp)/(dlf*(dlf+1.0_dp)+2.0_dp)
                         !            Dalphaei(ivar,ivt) = max(  Dalphaei(ivar,ivt) , dlf )

                         !            ! Smoothed limiter 2
                         !            dlf = (dwstar-DWc(ivar))/(DtIi(ivar)-SYS_EPSREAL_DP)
                         !            dlf = min(1.0_dp,dlf*dlf)
                         !            Dalphaei(ivar,ivt) = max(  Dalphaei(ivar,ivt) , dlf )

                      else
                         !Dalphaei(ivar,ivt) = min(Dalphaei(ivar,ivt),1.0_dp)
                         
                         ! This is the one I take
                         Dalphaei(ivar,ivt) = max(Dalphaei(ivar,ivt),1.0_dp)
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


             ! Limiting all second derivatives with the same limiting factor
             DmAlpha = 0.0_dp

             ! Get minimum of all correction factors of all vertices on this element
             do ivar = 1, nvar

                Dalpha(ivar, iel) = min(Dalpha(ivar, iel),minval(Dalphaei(ivar,1:NVE)))

                DmAlpha(ivar,ivar) = Dalpha(ivar, iel) !* (1.0_dp-SYS_EPSREAL_DP)

                Dquadraticgradient(ivar,1:3) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(4:6,iel)-1)

             end do

             DquadraticGradient = matmul(matmul(matmul(DR,DmAlpha),DL),DquadraticGradient)

             do ivar = 1, nvar  
                p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(4:6,iel)-1) = Dquadraticgradient(ivar,1:3)
             end do



             !      ! Limiting the second derivatives with different limiting factors
             !      DmAlpha = 0.0_dp
             !      
             !      ! Get minimum of all correction factors of all vertices on this element
             !      do ivar = 1, nvar
             !        
             !        Dquadraticgradient(ivar,1:3) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(4:6,iel)-1)
             !
             !      end do
             !      
             !      
             !      
             !      do ivar = 1, nvar
             !      
             !        DmAlpha(ivar,ivar) = Dalpha(ivar, iel)
             !
             !      end do
             !      
             !      DquadraticGradient(1:nvar,1) = matmul(matmul(matmul(DR,DmAlpha),DL),DquadraticGradient(1:nvar,1))
             !      
             !      do ivar = 1, nvar
             !      
             !        DmAlpha(ivar,ivar) = min(Dalpha(ivar, iel),minval(Dalphaei(ivar,1:NVE)))
             !
             !      end do
             !      
             !      DquadraticGradient(1:nvar,2) = matmul(matmul(matmul(DR,DmAlpha),DL),DquadraticGradient(1:nvar,2))
             !      
             !      do ivar = 1, nvar
             !      
             !        DmAlpha(ivar,ivar) = minval(Dalphaei(ivar,1:NVE))
             !
             !      end do
             !      
             !      DquadraticGradient(1:nvar,3) = matmul(matmul(matmul(DR,DmAlpha),DL),DquadraticGradient(1:nvar,3))
             !      
             !      
             !        
             !      do ivar = 1, nvar  
             !        p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(4:6,iel)-1) = Dquadraticgradient(ivar,1:3)
             !        
             !        Dalpha(ivar, iel) = min(Dalpha(ivar, iel),minval(Dalphaei(ivar,1:NVE)))
             !      end do





          case (3)

             ! Get minimum of all correction factors of all vertices on this element
             DmAlpha = 0.0_dp
             do ivar = 1, nvar
                Dalpha(ivar, iel) = max(Dalpha(ivar, iel),minval(Dalphaei(ivar,1:NVE)))

                DmAlpha(ivar,ivar) = Dalpha(ivar, iel) !* (1.0_dp-SYS_EPSREAL_DP)

                DlinearGradient(ivar,1:2) = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3,iel)-1)

             end do

             !      do idim = 1, 1
             !        ! Now we need the trafo matrices
             !        
             !        if (idim==1) then
             !          DL = buildMixedL2(DQchar,1.0_dp,0.0_dp)
             !          DR = buildMixedR2(DQchar,1.0_dp,0.0_dp)
             !        elseif (idim==2) then
             !          DL = buildMixedL2(DQchar,0.0_dp,1.0_dp)
             !          DR = buildMixedR2(DQchar,0.0_dp,1.0_dp)
             !        end if
             !        
             !        if (idim==1) then
             !          DTemp1(1:3,1:2) = matmul(matmul(matmul(DR,DmAlpha),DL),DlinearGradient)
             !        else
             !          DTemp2(1:3,1:2) = matmul(matmul(matmul(DR,DmAlpha),DL),DlinearGradient)
             !          do i=1,3
             !          do j=1,2
             !            if (abs(DTemp2(i,j))>abs(DTemp1(i,j))) DTemp1(i,j) = DTemp2(i,j)
             !          end do
             !          end do
             !        end if
             !        
             !      end do
             !      DlinearGradient=DTemp1(1:3,1:2)

             DlinearGradient = matmul(matmul(matmul(DR,DmAlpha),DL),DlinearGradient)  

             do ivar = 1, nvar
                p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2:3,iel)-1) = DlinearGradient(ivar,1:2)
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
    deallocate(DVec, DVei, DIi, DtIi, DtLinMax, DtLinMin, DltIi, DlIi, DQchar,DWc)
    deallocate(DL,DR,DL2,DR2,DlinearGradient, Dquadraticgradient, DmAlpha)
    deallocate(DLin, DtLin, Dalphaei)
    deallocate(IdofGlob)
    deallocate(DLinPartValues)
    if (celement == EL_DG_T2_2D) deallocate(DDerivativeQuadraticValues)

  end subroutine dg_realKuzmin





























  !****************************************************************************

  !<subroutine>  

  subroutine dg_quadraticLimiter_2 (rvector,ralpha)

    !<description>

    ! Limits the linear part of a dg_T1 element vector.

    !</description>

    !<input>
    !</input>

    !<inputoutput>

    ! A vector to limit
    type(t_vectorScalar), intent(inout) :: rvector  

    ! The limiting factors
    type(t_vectorScalar), intent(inout) :: ralpha

    !</inputoutput>

    !</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata, p_rAlpha_Ddata
    integer :: indof, NEL, iel, NVE, ivt, NVT

    ! The underlying triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! The underlying spatial discretisation
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr

    ! The coordinates of the points in which to evaluate the solution vector
    real(dp), dimension(2,17) :: Dpoints

    ! The list of elements, in which these points can be found
    integer, dimension(:),allocatable :: Ielements

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

    integer, dimension(:), pointer :: p_InodalProperty 

    real(DP), dimension(:,:), allocatable :: Dalpha

    real(DP), dimension(:,:,:), allocatable :: Dallvalues, DallPoints

    real(dp), dimension(2) :: Dvel

    ! Get pointer to the solution data
    call lsyssc_getbase_double (rvector,p_Ddata)

    call lsyssc_getbase_double (ralpha,p_rAlpha_Ddata)

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
    call storage_getbase_int(p_rtriangulation%h_InodalProperty ,&
         p_InodalProperty)



    ! Set pointer to coordinate vector
    call storage_getbase_double2D(&
         p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Set pointer to vertices at element
    call storage_getbase_int2D(&
         p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    NVT = p_rtriangulation%NVT

    allocate(Duimax(NVT),Duimin(NVT),Dalpha(3,NEL))
    allocate(DallValues(5,NEL,3),DallPoints(2,5,NEL))
    allocate(Ielements(NEL))


    ! Evaluate the solution and all of its derivatives in all points

    do iel= 1,nel

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
       Dallpoints(1,1,iel) = xc
       Dallpoints(2,1,iel) = yc
       Ielements(iel) = iel

       Dallpoints(1:2,2:5,iel) = p_DvertexCoords(1:2,p_IverticesAtElement(1:4,iel))


    end do


    call fevl_evaluate_sim1 (DER_FUNC, Dallvalues(1:5,1:NEL,1), rvector, Dallpoints(1:2,1:5,1:NEL), &
         Ielements)
    call fevl_evaluate_sim1 (DER_DERIV_X, Dallvalues(1:5,1:NEL,2), rvector, Dallpoints(1:2,1:5,1:NEL), &
         Ielements)
    call fevl_evaluate_sim1 (DER_DERIV_Y, Dallvalues(1:5,1:NEL,3), rvector, Dallpoints(1:2,1:5,1:NEL), &
         Ielements)



    Dalpha = 1.0_dp

    ! Now limit the x- and y- derivative and finally the linear part

    do ilim = 2,3




       duimax= -SYS_MAXREAL_DP
       duimin=  SYS_MAXREAL_DP

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

          ! Steady circular convection
          Dvel(1)=yc
          Dvel(2)=1.0_DP-xc

          Dvel(1) = Dvel(1)/(sqrt(Dvel(1)*Dvel(1)+Dvel(2)*Dvel(2)+SYS_EPSREAL_DP))
          Dvel(2) = Dvel(2)/(sqrt(Dvel(1)*Dvel(1)+Dvel(2)*Dvel(2)+SYS_EPSREAL_DP))

          duc = Dallvalues(1,iel,2)*Dvel(1) + Dallvalues(1,iel,3)*Dvel(2)

          if (ilim == 3) then
             duc = Dallvalues(1,iel,1)
          end if

          do ivt = 1, NVE
             nvt = p_IverticesAtElement(ivt,iel)
             duimax(nvt) = max(duc,duimax(nvt))
             duimin(nvt) = min(duc,duimin(nvt))
          end do



       end do



       do iel = 1, NEL


          ! No limiting of elements at boundary
          if ((p_InodalProperty(p_IverticesAtElement(1, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(2, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(3, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(4, iel))>0)) cycle


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

          ! Steady circular convection
          Dvel(1)=yc
          Dvel(2)=1.0_DP-xc

          Dvel(1) = Dvel(1)/(sqrt(Dvel(1)*Dvel(1)+Dvel(2)*Dvel(2)+SYS_EPSREAL_DP))
          Dvel(2) = Dvel(2)/(sqrt(Dvel(1)*Dvel(1)+Dvel(2)*Dvel(2)+SYS_EPSREAL_DP))

          duc = Dallvalues(1,iel,2)*Dvel(1) + Dallvalues(1,iel,3)*Dvel(2)

          if (ilim == 3) then
             duc = Dallvalues(1,iel,1)
          end if



          do ivert = 1, NVE  

             dui = Dallvalues(ivert+1,iel,2)*Dvel(1) + Dallvalues(ivert+1,iel,3)*Dvel(2)



             if (ilim == 3) then

                dui = Dallvalues(ivert+1,iel,1)

             end if

             ddu = dui-duc
             nvert = p_IverticesAtElement(ivert, iel)

             ! Find the maximum/minimum value of the solution in the centroids
             ! of all elements containing this vertex
             if (ddu > 10.0_dp*SYS_EPSREAL_DP) then
                dalphatemp = min(1.0_dp, (duimax(nvert)-duc)/ddu)

                !        ! Extremumfix
                !        if (duimax(nvert)-duc<SYS_EPSREAL_DP) dalphatemp = 1.0_dp

             elseif (ddu < -10.0_dp*SYS_EPSREAL_DP) then
                dalphatemp = min(1.0_dp, (duimin(nvert)-duc)/ddu)

                !        ! Extremumfix
                !        if (duimin(nvert)-duc>-SYS_EPSREAL_DP) dalphatemp = 1.0_dp

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


             ! Limiting like Kuzmin did it
             Dalpha(1,iel) = Dalpha(2,iel)

             ! Get global DOFs of the element
             call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)


             ! Multiply the linear part of the solution vector with the correction factor
             p_Ddata(IdofGlob(4:6)) = p_Ddata(IdofGlob(4:6))*Dalpha(1,iel)

             p_rAlpha_Ddata(IdofGlob(1)) = Dalpha(1,iel)
             p_rAlpha_Ddata(IdofGlob(2:3)) = 0.0_dp


             !        ! Limiting using different limiting factors for the second derivatives
             !        ! Get global DOFs of the element
             !        call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
             !    
             !    
             !        ! Multiply the linear part of the solution vector with the correction factor
             !        p_Ddata(IdofGlob(4)) = p_Ddata(IdofGlob(4))*Dalpha(1,iel)
             !        p_Ddata(IdofGlob(5)) = p_Ddata(IdofGlob(5))*min(Dalpha(1,iel),Dalpha(2,iel))
             !        p_Ddata(IdofGlob(6)) = p_Ddata(IdofGlob(6))*Dalpha(2,iel)
             !        
             !        
             !        Dalpha(1,iel) = min(Dalpha(1,iel),Dalpha(2,iel))
             !        
             !        p_rAlpha_Ddata(IdofGlob(1)) = Dalpha(1,iel)
             !        p_rAlpha_Ddata(IdofGlob(2:3)) = 0.0_dp


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
    deallocate(DallValues,DallPoints)
    deallocate(Ielements)


  end subroutine dg_quadraticLimiter_2



  !****************************************************************************

  !<subroutine>  

  subroutine dg_quadraticLimiter_3 (rvector,ralpha)

    !<description>

    ! Limits the linear part of a dg_T1 element vector.

    !</description>

    !<input>
    !</input>

    !<inputoutput>

    ! A vector to limit
    type(t_vectorScalar), intent(inout) :: rvector  

    ! The limiting factors
    type(t_vectorScalar), intent(inout) :: ralpha

    !</inputoutput>

    !</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_Ddata, p_rAlpha_Ddata
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

    integer :: iglobNeighNum

    integer, dimension(:), pointer :: p_InodalProperty 

    integer :: ilocalVtnumber, iglobVtNumber

    real(DP), dimension(:,:), allocatable :: Dalpha


    ! Get pointer to the solution data
    call lsyssc_getbase_double (rvector,p_Ddata)

    call lsyssc_getbase_double (ralpha,p_rAlpha_Ddata)

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
    call storage_getbase_int(p_rtriangulation%h_InodalProperty ,&
         p_InodalProperty)
    call storage_getbase_int(p_rtriangulation%h_IelementsAtVertex ,&
         p_IelementsAtVertex) 
    call storage_getbase_int(p_rtriangulation%h_IelementsAtVertexIdx ,&
         p_IelementsAtVertexIdx)



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

    do ilim = 1,1



       select case (ilim)
       case (1)
          ideriv = DER_FUNC
       end select


!!! First calculate the bounds

       duimax= -SYS_MAXREAL_DP
       duimin=  SYS_MAXREAL_DP

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

          do ivt = 1, NVE
             nvt = p_IverticesAtElement(ivt,iel)
             duimax(nvt) = max(duc,duimax(nvt))
             duimin(nvt) = min(duc,duimin(nvt))
          end do



!!!!!!!!!!!!!!!!!!!!!!!  NEU !!!!!!!!!!!!!!!!!!!!!!!

          if(ilim<3) then

             !    do ivert = 1, NVE
             !		
             !      nvert = p_IverticesAtElement(ivert, iel)
             !      
             !      ! The second point we want to evaluate the solution in, is in the corner of the mother element
             !      Dpoints(1,1+ivert) = p_DvertexCoords(1,nvert)
             !      Dpoints(2,1+ivert) = p_DvertexCoords(2,nvert)
             !      Ielements(1+ivert) = iel
             !	  end do
             !	  
             !    ! Evaluate the solution
             !    call fevl_evaluate (ideriv, Dvalues(1:5), rvector, Dpoints(1:2,1:5), &
             !                          Ielements(1:5))
             !    do ivt = 1, NVE
             !      nvt = p_IverticesAtElement(ivt,iel)
             !      duimax(nvt) = max(maxval(Dvalues(1:5)),duimax(nvt))
             !      duimin(nvt) = min(minval(Dvalues(1:5)),duimin(nvt))
             !    end do

          else

             !      call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
             !      duc = p_Ddata(IdofGlob(1))+abs(p_Ddata(IdofGlob(2)))+abs(p_Ddata(IdofGlob(3)))
             !      
             !    do ivt = 1, NVE
             !      nvt = p_IverticesAtElement(ivt,iel)
             !      duimax(nvt) = max(p_Ddata(IdofGlob(1))+abs(p_Ddata(IdofGlob(2)))+abs(p_Ddata(IdofGlob(3))),duimax(nvt))
             !      duimin(nvt) = min(p_Ddata(IdofGlob(1))-abs(p_Ddata(IdofGlob(2)))-abs(p_Ddata(IdofGlob(3))),duimin(nvt))
             !    end do

          end if

!!!!!!!!!!!!!!!!!!!!!!!  NEU !!!!!!!!!!!!!!!!!!!!!!!



       end do





       !  do ivt = 1, nvt
       !  iidx = 0
       !    do ineighbour = p_IelementsAtVertexIdx(ivt), p_IelementsAtVertexIdx(ivt+1)-1
       !    
       !      if(ilim<3) then
       !        iglobNeighNum = p_IelementsAtVertex(ineighbour)
       !      
       !        iidx = iidx + 1
       !      
       !        Dpoints(1,iidx) = p_DvertexCoords(1,iglobNeighNum)
       !        Dpoints(2,iidx) = p_DvertexCoords(2,iglobNeighNum)
       !        Ielements(iidx) = iglobNeighNum 
       !        
       !      else
       !        iidx = iidx + 1
       !        
       !        do ilocalVtNumber = 1, 4
       !          if (p_IverticesAtElement(ilocalVtnumber,iglobNeighNum )==ivt) exit
       !        end do
       !        
       !        call dof_locGlobMapping(p_rspatialDiscr, iglobNeighNum, IdofGlob)
       !        dui = p_Ddata(IdofGlob(1))
       !        
       !        select case (ilocalVtnumber)
       !        case(1)
       !        dui = dui - p_Ddata(IdofGlob(2)) - p_Ddata(IdofGlob(3))
       !        case(2)
       !        dui = dui + p_Ddata(IdofGlob(2)) - p_Ddata(IdofGlob(3))
       !        case(3)
       !        dui = dui + p_Ddata(IdofGlob(2)) + p_Ddata(IdofGlob(3))
       !        case(4)
       !        dui = dui - p_Ddata(IdofGlob(2)) + p_Ddata(IdofGlob(3))
       !        end select
       !        
       !        Dvalues(iidx) = dui
       !        
       !      end if
       !    end do
       !    
       !    ! Evaluate the solution
       !    if(ilim<3) call fevl_evaluate (ideriv, Dvalues(1:iidx), rvector, Dpoints(1:2,1:iidx), &
       !                                   Ielements(1:iidx))
       !      
       !      
       !    duimax(ivt) = max(maxval(Dvalues(1:iidx)),duimax(ivt))
       !    duimin(ivt) = min(minval(Dvalues(1:iidx)),duimin(ivt))
       !        
       !  
       !  end do


       !
       !
       do iel = 1, NEL

          iidx = 0

          do ivt = 1, 4

             nvt = p_IverticesAtElement(ivt,iel)

             do ineighbour = p_IelementsAtVertexIdx(nvt), p_IelementsAtVertexIdx(nvt+1)-1

                iglobNeighNum = p_IelementsAtVertex(ineighbour)
                iidx = iidx + 1



                ! Get midpoint of the element
                xc = &
                     (p_DvertexCoords(1,p_IverticesAtElement(1,iglobNeighNum))+&
                     p_DvertexCoords(1,p_IverticesAtElement(2,iglobNeighNum))+&
                     p_DvertexCoords(1,p_IverticesAtElement(3,iglobNeighNum))+&
                     p_DvertexCoords(1,p_IverticesAtElement(4,iglobNeighNum)))/4.0_dp

                yc = &
                     (p_DvertexCoords(2,p_IverticesAtElement(1,iglobNeighNum))+&
                     p_DvertexCoords(2,p_IverticesAtElement(2,iglobNeighNum))+&
                     p_DvertexCoords(2,p_IverticesAtElement(3,iglobNeighNum))+&
                     p_DvertexCoords(2,p_IverticesAtElement(4,iglobNeighNum)))/4.0_dp

                Dpoints(1,iidx) = xc
                Dpoints(2,iidx) = yc
                Ielements(iidx) = iglobNeighNum 

             end do
          end do

          ! Evaluate the solution
          if(ilim<3) call fevl_evaluate (ideriv, Dvalues(1:iidx), rvector, Dpoints(1:2,1:iidx), &
               Ielements(1:iidx))

          do ivt = 1, 4
             iglobVtNumber = p_IverticesAtElement(ivt,iel)
             duimax(iglobVtNumber) = max(maxval(Dvalues(1:iidx)),duimax(iglobVtNumber))
             duimin(iglobVtNumber) = min(minval(Dvalues(1:iidx)),duimin(iglobVtNumber))

          end do

       end do









!!! Start limiting



       do iel = 1, NEL


          ! No limiting of elements at boundary
          if ((p_InodalProperty(p_IverticesAtElement(1, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(2, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(3, iel))>0).or.(p_InodalProperty(p_IverticesAtElement(4, iel))>0)) cycle


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

          do ivert = 1, NVE  
             dui = Dvalues(1+ivert)
             ddu = dui-duc
             nvert = p_IverticesAtElement(ivert, iel)

             ! Find the maximum/minimum value of the solution in the centroids
             ! of all elements containing this vertex
             if (ddu > 0.0_dp) then
                dalphatemp = min(1.0_dp, (duimax(nvert)-duc)/ddu)

                !        ! Extremumfix
                !        if (duimax(nvert)-duc<SYS_EPSREAL_DP) dalphatemp = 1.0_dp
                if ((ilim>0).and.(duimax(nvert)-duc<abs(0.0001_dp*duc))) dalphatemp = 1.0_dp

             elseif (ddu < 0.0_dp) then
                dalphatemp = min(1.0_dp, (duimin(nvert)-duc)/ddu)

                !        ! Extremumfix
                !        if (duimin(nvert)-duc>-SYS_EPSREAL_DP) dalphatemp = 1.0_dp
                if ((ilim>0).and.(duimin(nvert)-duc>-abs(0.0001_dp*duc))) dalphatemp = 1.0_dp

             else ! (dui==duc)
                dalphatemp = 1.0_dp
             end if

             Dalpha(ilim,iel) = min(dalphatemp,Dalpha(ilim,iel))


          end do ! ivert

       end do ! iel




       select case (ilim)
       case (1)

          ! Now we have the limitingfactors for the quadratic part in Dalpha(1:2,:)
          ! Get the Minimum and multiply it with the corresponding DOFs
          do iel = 1, NEL  



             ! Get global DOFs of the element
             call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)


             ! Multiply the linear part of the solution vector with the correction factor
             p_Ddata(IdofGlob(2:6)) = p_Ddata(IdofGlob(2:6))*Dalpha(1,iel)

             p_rAlpha_Ddata(IdofGlob(1)) = Dalpha(1,iel)
             p_rAlpha_Ddata(IdofGlob(2:3)) = 0.0_dp


             !        ! Limiting using different limiting factors for the second derivatives
             !        ! Get global DOFs of the element
             !        call dof_locGlobMapping(p_rspatialDiscr, iel, IdofGlob)
             !    
             !    
             !        ! Multiply the linear part of the solution vector with the correction factor
             !        p_Ddata(IdofGlob(4)) = p_Ddata(IdofGlob(4))*Dalpha(1,iel)
             !        p_Ddata(IdofGlob(5)) = p_Ddata(IdofGlob(5))*min(Dalpha(1,iel),Dalpha(2,iel))
             !        p_Ddata(IdofGlob(6)) = p_Ddata(IdofGlob(6))*Dalpha(2,iel)
             !        
             !        
             !        Dalpha(1,iel) = min(Dalpha(1,iel),Dalpha(2,iel))
             !        
             !        p_rAlpha_Ddata(IdofGlob(1)) = Dalpha(1,iel)
             !        p_rAlpha_Ddata(IdofGlob(2:3)) = 0.0_dp


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


  end subroutine dg_quadraticLimiter_3






  ! Invert a dg-Massmatrix (block diagonal with small block size)
  subroutine dg_invertMassMatrix(rmatrix)

    type(t_matrixScalar), intent(inout) :: rmatrix


!!! Local variables

    ! The underlying triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! The underlying spatial discretisation
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr

    integer :: indof, NEL, iel, i, j, igel, ig, jg, iFESpace, iloc, iout, ii

    integer, dimension(6) :: IdofGlob

    type(t_elementDistribution), pointer :: p_relemDist

    real(dp), dimension(:,:), allocatable :: DlocMat, DilocMat

    INTEGER, DIMENSION(:), POINTER :: p_KLD, p_KCOL
    REAL(DP), DIMENSION(:), POINTER :: p_DA

    real(dp) :: ddet

    integer(I32) :: celement

    integer, dimension(:), pointer :: p_IelementList



    ! Get pointers of the matrix
    call lsyssc_getbase_Kcol (rmatrix,p_KCOL)
    call lsyssc_getbase_Kld (rmatrix,p_KLD)
    call lsyssc_getbase_double (rmatrix,p_DA)

    ! Get pointers for quicker access
    p_rspatialDiscr => rmatrix%p_rspatialDiscrTest
    p_rtriangulation => p_rspatialDiscr%p_rtriangulation

    ! Get number of elements
    NEL = p_rtriangulation%NEL

    ! Loop over all element distributions
    do iFESpace = 1, p_rspatialDiscr%inumFESpaces

       ! Get the element distribution from the spatial discretisation
       p_relemDist => p_rspatialDiscr%RelementDistr(iFESpace)  

       ! Get type of element in this distribution
       celement =  p_relemDist%celement

       ! Get number of local DOF
       indof = elem_igetNDofLoc(celement)

       ! Get list of all elements in this distribution
       call storage_getbase_int(p_relemDist%h_IelementList, p_IelementList)

       ! Allocate local matrix
       allocate(DlocMat(indof,indof),DilocMat(indof,indof))

       ! Loop over all elements in this element distrubution
       do iel = 1, p_relemDist%NEL

          ! Get global element number
          igel = p_IelementList(iel)

          call dof_locGlobMapping(p_rspatialDiscr, igel, IdofGlob)

!!! Get entries from the global into the local matrix
          ! Loop over all local lines
          do i = 1, indof
             ! Get global line
             ig = IdofGlob(i)

             ! Loop over all local columns
             do j = 1, indof
                ! Get global columns
                jg = IdofGlob(j)

                do iloc = p_KLD(ig), p_KLD(ig+1)-1
                   if (jg.eq.p_KCOL(iloc)) exit
                end do

                DlocMat(i,j) = p_DA(iloc)

             end do
          end do

!!! Invert local matrix
          if (indof.eq.3) then
             ddet = DlocMat(1,1)*(DlocMat(2,2)*DlocMat(3,3)-DlocMat(3,2)*DlocMat(2,3)) - DlocMat(2,2)*(DlocMat(1,2)*DlocMat(3,3)-DlocMat(3,2)*DlocMat(1,3)) + DlocMat(3,3)*(DlocMat(1,2)*DlocMat(2,3)-DlocMat(2,2)*DlocMat(1,3))

             DilocMat(1,1) = DlocMat(2,2) * DlocMat(3,3) - DlocMat(2,3) * DlocMat(3,2)
             DilocMat(1,2) = DlocMat(1,3) * DlocMat(3,2) - DlocMat(1,2) * DlocMat(3,3)
             DilocMat(1,3) = DlocMat(1,2) * DlocMat(2,3) - DlocMat(1,3) * DlocMat(2,2)

             DilocMat(2,1) = DlocMat(2,3) * DlocMat(3,1) - DlocMat(2,1) * DlocMat(3,3)
             DilocMat(2,2) = DlocMat(1,1) * DlocMat(3,3) - DlocMat(1,3) * DlocMat(3,1)
             DilocMat(2,3) = DlocMat(1,3) * DlocMat(2,1) - DlocMat(1,1) * DlocMat(2,3)

             DilocMat(3,1) = DlocMat(2,1) * DlocMat(3,2) - DlocMat(2,2) * DlocMat(3,1)
             DilocMat(3,2) = DlocMat(1,2) * DlocMat(3,1) - DlocMat(1,1) * DlocMat(3,2)
             DilocMat(3,3) = DlocMat(1,1) * DlocMat(2,2) - DlocMat(1,2) * DlocMat(2,1)

             DilocMat = 1/ddet*DilocMat

          elseif ((indof.eq.6).or.(indof.eq.1)) then
             !ddet = DlocMat(1,1)*(DlocMat(2,2)*DlocMat(3,3)*DlocMat(4,4)*(DlocMat(5,5)*DlocMat(6,6)

             DilocMat = 0.0_dp

             do ii = 1, indof
                DilocMat(ii,ii) = 1.0_dp/DlocMat(ii,ii)
             end do
          end if

          !      ! Output local matrices
          !      write(*,*) ''
          !      write(*,*) 'MC'
          !      do iout = 1,3
          !        write(*,*) DlocMat(iout,1),DlocMat(iout,2),DlocMat(iout,3)
          !      end do
          !      write(*,*) 'iMC'
          !      do iout = 1,3
          !        write(*,*) DilocMat(iout,1),DilocMat(iout,2),DilocMat(iout,3)
          !      end do
          !      pause

          !      ! Output local matrices
          !      write(*,*) ''
          !      write(*,*) 'MC'
          !      do iout = 1,6
          !        write(*,*) DlocMat(iout,1),DlocMat(iout,2),DlocMat(iout,3),DlocMat(iout,4),DlocMat(iout,5),DlocMat(iout,6)
          !      end do
          !      write(*,*) 'iMC'
          !      do iout = 1,6
          !        write(*,*) DilocMat(iout,1),DilocMat(iout,2),DilocMat(iout,3),DilocMat(iout,4),DilocMat(iout,5),DilocMat(iout,6)
          !      end do
          !      pause


!!! Get entries from the local into the global matrix
          ! Loop over all local lines
          do i = 1, indof
             ! Get global line
             ig = IdofGlob(i)

             ! Loop over all local columns
             do j = 1, indof
                ! Get global columns
                jg = IdofGlob(j)

                do iloc = p_KLD(ig), p_KLD(ig+1)-1
                   if (jg.eq.p_KCOL(iloc)) exit
                end do

                p_DA(iloc) = DilocMat(i,j)

             end do
          end do

       end do

       ! Deallocate local matrix
       deallocate(DlocMat,DilocMat)

    end do


  end subroutine dg_invertMassMatrix











  !****************************************************************************

  !<subroutine>

  subroutine dg_initAssembly_reverseCubPoints(rvectorAssembly,rform,&
       celement,ccubType,nelementsPerBlock)

    !<description>
    ! Initialise a vector assembly structure for assembling a linear form.
    !</description>

    !<input>
    ! The bilinear form specifying the underlying PDE of the discretisation.
    type(t_linearForm), intent(in) :: rform

    ! Type of element in the test space.
    integer(I32), intent(in) :: celement

    ! Type of cubature formula to use.
    integer(I32), intent(in) :: ccubType

    ! Optional: Maximum number of elements to process simultaneously.
    ! If not specified, LINF_NELEMSIM is assumed.
    integer, intent(in), optional :: nelementsPerBlock
    !</input>

    !<output>
    ! A vector assembly structure.
    type(t_linfVectorAssembly), intent(out) :: rvectorAssembly
    !</output>

    !</subroutine>

    ! local variables
    integer :: i,i1

    real(dp) :: dtemp
    real(dp), dimension(2) :: dtemp2

    ! Initialise the structure.
    rvectorAssembly%rform = rform
    rvectorAssembly%ccubType = ccubType
    rvectorAssembly%nelementsPerBlock = LINF_NELEMSIM
    if (present(nelementsPerBlock)) &
         rvectorAssembly%nelementsPerBlock = nelementsPerBlock
    rvectorAssembly%celement = celement

    ! Get the number of local DOF`s for trial and test functions
    rvectorAssembly%indof = elem_igetNDofLoc(celement)

    ! Which derivatives of basis functions are needed?
    ! Check the descriptors of the bilinear form and set BDERxxxx
    ! according to these.
    rvectorAssembly%Bder(:) = .false.

    ! Loop through the additive terms
    do i = 1,rform%itermCount
       ! The desriptor Idescriptors gives directly the derivative
       ! which is to be computed! Build templates for BDER.
       ! We do not compute the actual BDER here, as there might be some special
       ! processing if trial/test functions are identical!
       !
       ! At first build the descriptors for the trial functions
       I1=rform%Idescriptors(i)

       if ((I1 .le.0) .or. (I1 .gt. DER_MAXNDER)) then
          call output_line ('Invalid descriptor!',&
               OU_CLASS_ERROR,OU_MODE_STD,'linf_initAssembly')
          call sys_halt()
       endif

       rvectorAssembly%Bder(I1)=.true.
    end do

    ! Get the number of vertices of the element, specifying the transformation
    ! form the reference to the real element.
    rvectorAssembly%NVE = elem_igetNVE(celement)

    ! Get from the element space the type of coordinate system
    ! that is used there:
    rvectorAssembly%ctrafoType = elem_igetTrafoType(celement)

    ! Get the number of cubature points for the cubature formula
    rvectorAssembly%ncubp = cub_igetNumPts(ccubType)

    ! Allocate two arrays for the points and the weights
    allocate(rvectorAssembly%p_Domega(rvectorAssembly%ncubp))
    allocate(rvectorAssembly%p_DcubPtsRef(&
         trafo_igetReferenceDimension(rvectorAssembly%ctrafoType),&
         rvectorAssembly%ncubp))

    ! Get the cubature formula
    call cub_getCubature(ccubType,rvectorAssembly%p_DcubPtsRef,rvectorAssembly%p_Domega)

    ! Revert the numbering of cubature points
    do i = 1, size(rvectorAssembly%p_Domega)/2
       dtemp2(:) = rvectorAssembly%p_DcubPtsRef(:,i)
       rvectorAssembly%p_DcubPtsRef(:,i) = rvectorAssembly%p_DcubPtsRef(:,size(rvectorAssembly%p_Domega)+1-i)
       rvectorAssembly%p_DcubPtsRef(:,size(rvectorAssembly%p_Domega)+1-i) = dtemp2(:)

       dtemp = rvectorAssembly%p_Domega(i)
       rvectorAssembly%p_Domega(i) = rvectorAssembly%p_Domega(size(rvectorAssembly%p_Domega)+1-i)
       rvectorAssembly%p_Domega(size(rvectorAssembly%p_Domega)+1-i) = dtemp

    end do

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    rvectorAssembly%cevaluationTag = elem_getEvaluationTag(rvectorAssembly%celement)

  end subroutine dg_initAssembly_reverseCubPoints

































  !****************************************************************************

  !<subroutine>

  subroutine bilf_dg_buildMatrixScEdge2D (rform, ccubType, bclear, rmatrix,&
       rvectorSol, raddTriaData,&
       flux_dg_buildMatrixScEdge2D_sim,&
       rcollection, cconstrType)

    !<description>
    ! This routine calculates the entries of a finite element matrix in 2D.
    ! The matrix structure must be prepared with bilf_createMatrixStructure
    ! in advance.
    ! In case the array for the matrix entries does not exist, the routine
    ! allocates memory in size of the matrix of the heap for the matrix entries.
    !
    ! For setting up the entries, the discretisation structure attached to
    ! the matrix is used (rmatrix%p_rdiscretisation). This is
    ! normally attached to the matrix by bilf_createMatrixStructure.
    !
    ! The matrix must be unsorted when this routine is called, 
    ! otherwise an error is thrown.
    !</description>

    !<input>
    ! The bilinear form specifying the underlying PDE of the discretisation.
    type(t_bilinearForm), intent(in) :: rform

    ! A line cubature formula CUB_xxxx_1D to be used for line integration.
    integer(I32), intent(in) :: ccubType

    ! The solution vector. Used to calculate the solution on the edges.
    type(t_vectorScalar), intent(in) :: rvectorSol

    ! Additional triangulation data
    type(t_additionalTriaData), intent(in) :: raddTriaData

    ! Whether to clear the matrix before calculating the entries.
    ! If .FALSE., the new matrix entries are added to the existing entries.
    logical, intent(in) :: bclear

    ! A callback routine for the flux function.
    include 'intf_flux_dg_buildMatrixScEdge2D.inc'
    optional :: flux_dg_buildMatrixScEdge2D_sim

    ! OPTIONAL: One of the BILF_MATC_xxxx constants that allow to specify
    ! the matrix construction method. If not specified,
    ! BILF_MATC_EDGEBASED is used.
    integer, intent(in), optional :: cconstrType
    !</input>

    !<inputoutput>
    ! The FE matrix. Calculated matrix entries are imposed to this matrix.
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: A collection structure. This structure is given to the
    ! callback function for nonconstant coefficients to provide additional
    ! information. 
    type(t_collection), intent(inout), target, optional :: rcollection
    !</inputoutput>

    !</subroutine>

    ! local variables
    type(t_bilfMatrixAssembly), dimension(2) :: rmatrixAssembly
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: IelementList, p_IedgeList
    integer :: ccType
    integer :: iedge, ielementDistr

    ! The matrix must be unsorted, otherwise we can not set up the matrix.
    ! Note that we cannot switch off the sorting as easy as in the case
    ! of a vector, since there is a structure behind the matrix! So the caller
    ! has to make sure, the matrix is unsorted when this routine is called.
    if (rmatrix%isortStrategy .gt. 0) then
       call output_line ('Matrix-structure must be unsorted!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ! The matrix must provide discretisation structures
    if ((.not. associated(rmatrix%p_rspatialDiscrTest)) .or. &
         (.not. associated(rmatrix%p_rspatialDiscrTrial))) then
       call output_line ('No discretisation associated!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ! The discretisation must provide a triangulation structure
    if ((.not. associated(rmatrix%p_rspatialDiscrTest%p_rtriangulation)) .or. &
         (.not. associated(rmatrix%p_rspatialDiscrTrial%p_rtriangulation))) then
       call output_line('No triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    !  ! The discretisation must provide a boundary structure
    !  if ((.not. associated(rmatrix%p_rspatialDiscrTest%p_rboundary)) .or. &
    !      (.not. associated(rmatrix%p_rspatialDiscrTrial%p_rboundary))) then
    !    call output_line('No boundary associated!',&
    !        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
    !    call sys_halt()
    !  end if

    !  ! Set pointers for quicker access
    !  p_rboundary => rmatrix%p_rspatialDiscrTest%p_rboundary
    !  if (.not.associated(p_rboundary, rmatrix%p_rspatialDiscrTrial%p_rboundary)) then
    !    call output_line('Invalid boundary associated!',&
    !        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
    !    call sys_halt()
    !  end if

    p_rtriangulation => rmatrix%p_rspatialDiscrTest%p_rtriangulation
    if (.not.associated(p_rtriangulation, rmatrix%p_rspatialDiscrTrial%p_rtriangulation)) then
       call output_line('Invalid triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ccType = BILF_MATC_EDGEBASED
    if (present(cconstrType)) ccType = cconstrType

    ! Do we have a uniform triangulation? Would simplify a lot...
    select case (rmatrix%p_rspatialDiscrTest%ccomplexity)
    case (SPDISC_UNIFORM,SPDISC_CONFORMAL) 
       ! Uniform and conformal discretisations
       select case (rmatrix%cdataType)
       case (ST_DOUBLE) 
          ! Which matrix structure do we have?
          select case (rmatrix%cmatrixFormat) 
          case (LSYSSC_MATRIX9)

             ! Probably allocate/clear the matrix
             !        if (rmatrix%h_DA .eq. ST_NOHANDLE) then
             !          call lsyssc_allocEmptyMatrix(rmatrix,LSYSSC_SETM_ZERO)
             !        else
             if (bclear) call lsyssc_clearMatrix (rmatrix)
             !        end if


             ! Allocate the edgelist
             allocate(p_IedgeList(rmatrix%p_rspatialDiscrTest%p_rtriangulation%NMT))

             ! All edges
             forall (iedge = 1:rmatrix%p_rspatialDiscrTest%p_rtriangulation%NMT) p_IedgeList(iedge)=iedge

             ! Initialise a matrix assembly structure for that element distribution
             ielementDistr = 1
             call bilf_initAssembly(rmatrixAssembly(1),rform,&
                  rmatrix%p_rspatialDiscrTest%RelementDistr(ielementDistr)%celement,&
                  rmatrix%p_rspatialDiscrTrial%RelementDistr(ielementDistr)%celement,&
                  ccubType, BILF_NELEMSIM)

             ! Do the same for the other side of the egde
             ielementDistr = 1
             call dg_bilf_initAssembly_reverseCubPoints(rmatrixAssembly(2),rform,&
                  rmatrix%p_rspatialDiscrTest%RelementDistr(ielementDistr)%celement,&
                  rmatrix%p_rspatialDiscrTrial%RelementDistr(ielementDistr)%celement,&
                  ccubType, BILF_NELEMSIM)

             ! Assemble the data for all elements in this element distribution
             call dg_bilf_assembleSubmeshMat9Bdr2D (rmatrixAssembly, rmatrix,&
                  rvectorSol, raddTriaData, p_IedgeList(1:rmatrix%p_rspatialDiscrTest%p_rtriangulation%NMT),&
                  ccType, flux_dg_buildMatrixScEdge2D_sim, rcollection)


             ! Release the assembly structure.
             call bilf_doneAssembly(rmatrixAssembly(1))
             call bilf_doneAssembly(rmatrixAssembly(2))

             ! Deallocate the edgelist
             deallocate(p_IedgeList)



          case default
             call output_line ('Not supported matrix structure!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
             call sys_halt()
          end select

       case default
          call output_line ('Single precision matrices currently not supported!', &
               OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
          call sys_halt()
       end select

    case default
       call output_line ('General discretisation not implemented!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end select

  end subroutine bilf_dg_buildMatrixScEdge2D



  !****************************************************************************

  !<subroutine>  

  subroutine dg_bilf_assembleSubmeshMat9Bdr2D (rmatrixAssembly, rmatrix,&
       rvectorSol, raddTriaData, IedgeList,&
       cconstrType, flux_dg_buildMatrixScEdge2D_sim, rcollection)

    !<description>

    ! Assembles the matrix entries for a submesh by integrating over the
    ! boundary region in 2D.

    !</description>

    !<input>

    ! List of elements where to assemble the bilinear form.
    integer, dimension(:), intent(in), target :: IedgeList

    ! One of the BILF_MATC_xxxx constants that allow to specify the
    ! matrix construction method.
    integer, intent(in) :: cconstrType

    ! The solution vector. Used to calculate the solution on the edges.
    type(t_vectorScalar), intent(in) :: rvectorSol

    ! Additional triangulation data
    type(t_additionalTriaData), intent(in) :: raddTriaData

    ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
    ! Must be present if the matrix has nonconstant coefficients!
    include 'intf_flux_dg_buildMatrixScEdge2D.inc'
    optional :: flux_dg_buildMatrixScEdge2D_sim

    !</input>

    !<inputoutput>

    ! A matrix assembly structure prepared with bilf_initAssembly.
    type(t_bilfMatrixAssembly), intent(inout), dimension(2), target :: rmatrixAssembly

    ! A matrix where to assemble the contributions to.
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: A pointer to a collection structure. This structure is given to the
    ! callback function for nonconstant coefficients to provide additional
    ! information. 
    type(t_collection), intent(inout), target, optional :: rcollection

    !</inputoutput>

    !</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_DA
    integer :: indofTest,indofTrial,ncubp

    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax,ibdc,k
    integer :: iel,icubp,ialbet,ia,ib,idofe,jdofe,nve
    real(DP) :: domega1,domega2,daux1,daux2,db1,db2,dlen
    integer(I32) :: cevaluationTag
    type(t_bilfMatrixAssembly), dimension(2), target :: rlocalMatrixAssembly
    type(t_domainIntSubset), dimension(2) :: rintSubset
    integer, dimension(:,:,:), pointer :: p_Kentryii, p_Kentryia, p_Kentryai, p_Kentryaa
    real(DP), dimension(:,:,:), pointer :: p_Dentryii, p_Dentryia, p_Dentryai, p_Dentryaa
    real(DP), dimension(:,:,:), pointer :: p_Dcoords
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_Dside
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial
    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    real(DP), dimension(:), pointer :: p_DcoefficientsBilf
    integer, dimension(:,:), pointer :: p_IdofsTest
    integer, dimension(:,:), pointer :: p_IdofsTrial
    type(t_evalElementSet), pointer :: p_revalElementSet
    integer, dimension(:,:),pointer :: p_Idescriptors
    integer, dimension(:,:), allocatable, target :: IelementList

    ! Pointer to Ielementsatedge in the triangulation
    integer, dimension(:,:), pointer :: p_IelementsAtEdge

    ! Pointer to IverticesAtEdge in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtEdge

    ! Pointer to the vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Pointer to IverticesAtEelement in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Space for the values of the flux function
    real(DP), dimension(:,:,:), allocatable :: DfluxValues

    ! Arrays for cubature points 1D->2D
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D_1, Dxi1D_2
    real(DP), dimension(:,:,:,:), allocatable :: Dxi2D,DpointsRef
    real(DP), dimension(:,:), allocatable :: DpointsPar
    real(DP), dimension(:), allocatable :: DedgeLength

    integer(i32) :: icoordSystem
    integer :: NEL
    integer :: iside
    logical :: bisLinearTrafo

    !    ! Boundary component?
    !    ibdc = rboundaryRegion%iboundCompIdx

    ! Get some pointers for faster access
    call lsyssc_getbase_double (rmatrix,p_DA)
    indofTest = rmatrixAssembly(1)%indofTest
    indofTrial = rmatrixAssembly(1)%indofTrial
    ncubp = rmatrixAssembly(1)%ncubp

    ! Open-MP-Extension: Copy the matrix assembly data to the local
    ! matrix assembly data, where we can allocate memory.
    !
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the matrix
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original matrix assembly
    ! stucture or disturbing the data of the other processors.
    !
    !$omp parallel default(shared) &
    !$omp private(DedgeLength,DpointsPar,DpointsRef,Dxi1D,Dxi2D,IELmax,bisLinearTrafo,&
    !$omp         cevaluationTag,daux,db,dlen,domega,ia,ialbet,ib,icoordSystem,icubp,&
    !$omp         idofe,iel,jdofe,k,p_DbasTest,p_DbasTrial,p_Dcoefficients,&
    !$omp         p_DcoefficientsBilf,p_Dcoords,p_DcubPtsRef,p_Dentry,p_Domega,&
    !$omp         p_Idescriptors,p_IdofsTest,p_IdofsTrial,p_Kentry,&
    !$omp         p_revalElementSet,rintSubset,rlocalMatrixAssembly)
    rlocalMatrixAssembly(1) = rmatrixAssembly(1)
    rlocalMatrixAssembly(2) = rmatrixAssembly(2)
    call bilf_allocAssemblyData(rlocalMatrixAssembly(1))
    call bilf_allocAssemblyData(rlocalMatrixAssembly(2))

    ! Allocate space for the positions of the DOFs in the matrix
    allocate(p_Kentryii(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryia(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryai(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryaa(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))

    ! Allocate space for the coefficient of the solutions DOFs on each side of the edge
    allocate(p_Dside(2,rlocalMatrixAssembly(1)%rform%itermcount,ncubp,rmatrixAssembly(1)%nelementsPerBlock))

    ! Allocate space for the entries in the local matrices
    allocate(p_Dentryii(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryia(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryai(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryaa(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))



    ! Allocate space for the flux variables DIM(nvar,ialbet,ncubp,elementsperblock)
    allocate(DfluxValues(rlocalMatrixAssembly(1)%rform%itermcount,ncubp,rlocalMatrixAssembly(1)%nelementsPerBlock))

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


    ! Get number of elements
    NEL = rmatrix%p_rspatialDiscrTest%p_rtriangulation%NEL

    ! Get pointers to elements at edge
    call storage_getbase_int2D(&
         rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_IelementsAtEdge,&
         p_IelementsAtEdge)

    ! Get pointers to the vertex coordinates
    call storage_getbase_double2D(&
         rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_DvertexCoords,&
         p_DvertexCoords)

    ! Get pointers to vertices at edge
    call storage_getbase_int2D(&
         rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtEdge,&
         p_IverticesAtEdge)

    ! Get pointers to vertices at elements
    call storage_getbase_int2D(&
         rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtElement,&
         p_IverticesAtElement)   

    ! Get the elements adjacent to the given edges
    allocate(IelementList(3,size(IedgeList)))
    IelementList(1:2,1:size(IedgeList))=p_IelementsAtEdge(1:2,IedgeList(:))

    ! Copy the second component and replace 0s by 1s
    IelementList(3,size(IedgeList))=IelementList(2,size(IedgeList))
    do iel = 1,size(IedgeList)
       IelementList(3,iel)=max(IelementList(2,iel),1)
    end do

    ! Transpose the coordinate array such that we get coordinates we
    ! can work with in the mapping between 1D and 2D.
    do k = 1, ubound(rlocalmatrixAssembly(1)%p_DcubPtsRef,1)
       do icubp = 1,ncubp
          Dxi1D_1(icubp,k) = rlocalmatrixAssembly(1)%p_DcubPtsRef(k,icubp)
          Dxi1D_2(icubp,k) = rlocalmatrixAssembly(2)%p_DcubPtsRef(k,icubp)
       end do
    end do

    ! Allocate memory for the cubature points in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,2,rlocalMatrixAssembly(1)%nelementsPerBlock))

    ! Allocate memory for the coordinates of the reference points
    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalMatrixAssembly(1)%nelementsPerBlock,2))

    !    ! Allocate memory for the parameter values of the points on the boundary
    !    allocate(DpointsPar(ncubp,rlocalMatrixAssembly%nelementsPerBlock))
    !
    !    ! Allocate memory for the length of edges on the boundary
    !    allocate(DedgeLength(rlocalMatrixAssembly%nelementsPerBlock))

    ! Get the type of coordinate system
    icoordSystem = elem_igetCoordSystem(rlocalMatrixAssembly(1)%celementTrial)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !$omp do schedule(static,1)
    do IELset = 1, size(IedgeList), rlocalMatrixAssembly(1)%nelementsPerBlock

       ! We always handle nelementsPerBlock elements simultaneously.
       ! How many elements have we actually here?
       ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
       ! elements simultaneously.

       IELmax = min(size(IedgeList),IELset-1+rlocalMatrixAssembly(1)%nelementsPerBlock)

       ! Map the 1D cubature points to the edges in 2D.
       do iel = 1,IELmax-IELset+1
          call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(1,Iedgelist(IELset+iel-1)), &
               ncubp, Dxi1D_1, Dxi2D(:,:,1,iel))
          call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(2,Iedgelist(IELset+iel-1)), &
               ncubp, Dxi1D_2, Dxi2D(:,:,2,iel))
       end do

       !      ! Calculate the parameter values of the points
       !      do iel = 1,IELmax-IELset+1
       !        do icubp = 1,ncubp
       !          ! Dxi1D is in [-1,1] while the current edge has parmeter values
       !          ! [DedgePosition(1),DedgePosition(2)]. So do a linear
       !          ! transformation to transform Dxi1D into that interval, this 
       !          ! gives the parameter values in length parametrisation
       !          call mprim_linearRescale(Dxi1D(icubp,1), -1.0_DP, 1.0_DP,&
       !              DedgePosition(1,IELset+iel-1), DedgePosition(2,IELset+iel-1),&
       !              DpointsPar(icubp,iel))
       !        end do
       !      end do

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
       ! global DOF`s of our BILF_NELEMSIM elements simultaneously.
       !      call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTest, &
       !          IelementList(IELset:IELmax), p_IdofsTest)
       call dof_locGlobMapping_mult( rmatrix%p_rspatialDiscrTest, &
            IelementList(1,IELset:IELmax), rlocalMatrixAssembly(1)%p_IdofsTest)
       call dof_locGlobMapping_mult( rmatrix%p_rspatialDiscrTest, &
            IelementList(3,IELset:IELmax), rlocalMatrixAssembly(2)%p_IdofsTest)


       !      ! If the DOF`s for the trial functions are different, calculate them, too.
       !      if (.not. rlocalMatrixAssembly%bIdenticalTrialAndTest) then
       !        call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTrial, &
       !            IelementList(IELset:IELmax), p_IdofsTrial)
       !      end if

       ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------

       ! For the assembly of the global matrix, we use a "local"
       ! approach. At first we build a "local" system matrix according
       ! to the current element. This contains all additive
       ! contributions of element iel, which are later added at the
       ! right positions to the elements in the global system matrix.
       !
       ! We have indofTrial trial DOF`s per element and
       ! indofTest test DOF`s per element. Therefore there are
       ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
       ! "active" (i.e. have common support) on our current element, each 
       ! giving an additive contribution to the system matrix.
       !
       ! We build a quadratic indofTrial*indofTest local matrix:
       ! Kentry(1..indofTrial,1..indofTest) receives the position 
       ! in the global system matrix, where the corresponding value 
       ! has to be added to.
       ! (The corresponding contributions can be saved separately, 
       ! but we directly add them to the global matrix in this 
       ! approach.)
       !
       ! We build local matrices for all our elements 
       ! in the set simultaneously. Get the positions of the local matrices
       ! in the global matrix.
       !      call bilf_getLocalMatrixIndices (rmatrix,p_IdofsTest,p_IdofsTrial,p_Kentry,&
       !          ubound(p_IdofsTest,1),ubound(p_IdofsTrial,1),IELmax-IELset+1)  
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
       !            rlocalMatrixAssembly(1)%p_IdofsTrial, p_Kentryii,&
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
       !            rlocalMatrixAssembly(2)%p_IdofsTrial, p_Kentryia,&
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
       !            rlocalMatrixAssembly(1)%p_IdofsTrial, p_Kentryai,&
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
       !            rlocalMatrixAssembly(2)%p_IdofsTrial, p_Kentryaa,&
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTrial,1), IELmax-IELset+1)



       call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
            rlocalMatrixAssembly(1)%p_IdofsTest, p_Kentryii,&
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
            rlocalMatrixAssembly(2)%p_IdofsTest, p_Kentryai,&
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
            rlocalMatrixAssembly(1)%p_IdofsTest, p_Kentryia,&
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
            rlocalMatrixAssembly(2)%p_IdofsTest, p_Kentryaa,&
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), IELmax-IELset+1)   

       ! -------------------- ELEMENT EVALUATION PHASE ----------------------

       ! Ok, we found the positions of the local matrix entries
       ! that we have to change.
       ! To calculate the matrix contributions, we have to evaluate
       ! the elements to give us the values of the basis functions
       ! in all the DOF`s in all the elements in our set.

       ! Get the element evaluation tag of all FE spaces. We need it to evaluate
       ! the elements later. All of them can be combined with OR, what will give
       ! a combined evaluation tag. 
       cevaluationTag = rlocalMatrixAssembly(1)%cevaluationTag

       ! The cubature points are already initialised by 1D->2D mapping.
       cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

       !      ! Do we have a (multi-)linear transformation?
       !      bisLinearTrafo = trafo_isLinearTrafo(rlocalMatrixAssembly%ctrafoType)
       !
       !      if (bisLinearTrafo) then
       !        ! We need the vertices of the element corners and the number
       !        ! of vertices per element to compute the length of the element
       !        ! edge at the boundary
       !        cevaluationTag = ior(cevaluationTag, EL_EVLTAG_COORDS)
       !        nve = trafo_igetNVE(rlocalMatrixAssembly%ctrafoType)
       !      end if

       ! Calculate all information that is necessary to evaluate the finite element
       ! on all cells of our subset. This includes the coordinates of the points
       ! on the cells.
       !      call elprep_prepareSetForEvaluation (p_revalElementSet,&
       !          cevaluationTag, rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
       !          IelementList(IELset:IELmax), rlocalMatrixAssembly%ctrafoType, &
       !          DpointsRef=DpointsRef)
       !      p_Dcoords => p_revalElementSet%p_Dcoords

       call elprep_prepareSetForEvaluation (&
            rlocalMatrixAssembly(1)%revalElementSet,&
            cevaluationTag,  rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
            IelementList(1,IELset:IELmax), rlocalMatrixAssembly(1)%ctrafoType, &
            DpointsRef=DpointsRef(:,:,:,1))
       call elprep_prepareSetForEvaluation (&
            rlocalMatrixAssembly(2)%revalElementSet,&
            cevaluationTag,  rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
            IelementList(3,IELset:IELmax), rlocalMatrixAssembly(2)%ctrafoType, &
            DpointsRef=DpointsRef(:,:,:,2))




       !      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
       !      if (.not. rlocalMatrixAssembly%rform%ballCoeffConstant) then
       !        if (present(fcoeff_buildMatrixScBdr2D_sim)) then
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

       ! If the flux function needs other, than just the function values from the solution
       ! (for example the derivatives), we will give an evalElementSet to it
       ! This is filled here

       call domint_initIntegrationByEvalSet (rlocalMatrixAssembly(1)%revalElementSet, rintSubset(1))
       call domint_initIntegrationByEvalSet (rlocalMatrixAssembly(2)%revalElementSet, rintSubset(2))
       !rintSubset(1)%ielementDistribution = 0
       rintSubset(1)%ielementStartIdx = IELset
       rintSubset(1)%p_Ielements => IelementList(1,IELset:IELmax)
       rintSubset(1)%p_IdofsTrial => rlocalMatrixAssembly(1)%p_IdofsTest
       rintSubset(1)%celement = rlocalMatrixAssembly(1)%celementTest
       !rintSubset(2)%ielementDistribution = 0
       rintSubset(2)%ielementStartIdx = IELset
       rintSubset(2)%p_Ielements => IelementList(2,IELset:IELmax)
       rintSubset(2)%p_IdofsTrial => rlocalMatrixAssembly(2)%p_IdofsTest
       rintSubset(2)%celement = rlocalMatrixAssembly(2)%celementTest


       call flux_dg_buildMatrixScEdge2D_sim (&
            !            rlocalVectorAssembly(1)%p_Dcoefficients(1,:,1:IELmax-IELset+1),&
            !            DsolVals(:,:,1:IELmax-IELset+1),&
       DfluxValues(:,:,1:IELmax-IELset+1),&
            rvectorSol,&
            IelementList(2,IELset:IELmax),&
            p_Dside,&
            raddTriaData%p_Dnormals(:,Iedgelist(IELset:IELmax)),&
            !DpointsReal(1:ndim2d,1:ncubp,1:IELmax-IELset+1),&
       rintSubset,&
            rcollection )


       call domint_doneIntegration (rintSubset(1))
       call domint_doneIntegration (rintSubset(2))





       ! Calculate the values of the basis functions.
       !      call elem_generic_sim2 (rlocalMatrixAssembly%celementTest, &
       !          p_revalElementSet, rlocalMatrixAssembly%BderTest, &
       !          rlocalMatrixAssembly%p_DbasTest) 
       call elem_generic_sim2 (rlocalMatrixAssembly(1)%celementTest, &
            rlocalMatrixAssembly(1)%revalElementSet,&
            rlocalMatrixAssembly(1)%BderTest, &
            rlocalMatrixAssembly(1)%p_DbasTest)
       call elem_generic_sim2 (rlocalMatrixAssembly(2)%celementTest, &
            rlocalMatrixAssembly(2)%revalElementSet,&
            rlocalMatrixAssembly(2)%BderTest, &
            rlocalMatrixAssembly(2)%p_DbasTest)

       !      ! Omit the calculation of the trial function values if they
       !      ! are identical to the test function values.
       !      if (.not. rlocalMatrixAssembly%bidenticalTrialAndTest) then
       !        call elem_generic_sim2 (rlocalMatrixAssembly%celementTrial, &
       !            p_revalElementSet, rlocalMatrixAssembly%BderTrial, &
       !            rlocalMatrixAssembly%p_DbasTrial)
       !      end if

       !      ! Calculate the length of egdes on the boundary. Depending on
       !      ! whether the transformation is (multi-)linear or not we compute
       !      ! the edge length as the distance between the two corner
       !      ! vertices of the element located on the boundary or as the real
       !      ! length of the boundary segment of the element.
       !      !
       !      ! The length of the current edge serves as a "determinant" in
       !      ! the cubature, so we have to divide it by 2 as an edge on the
       !      ! unit interval [-1,1] has length 2.
       !      if (bisLinearTrafo) then
       !        do iel = 1,IELmax-IELset+1
       !          DedgeLength(iel) = 0.5_DP*sqrt(&
       !              (p_Dcoords(1,    IelementOrientation(IELset+iel-1),iel)-&
       !               p_Dcoords(1,mod(IelementOrientation(IELset+iel-1),nve)+1,iel))**2+&
       !              (p_Dcoords(2,    IelementOrientation(IELset+iel-1),iel)-&
       !               p_Dcoords(2,mod(IelementOrientation(IELset+iel-1),nve)+1,iel))**2)
       !        end do
       !      else
       !        do iel = 1,IELmax-IELset+1
       !          DedgeLength(iel) = 0.5_DP*(DedgePosition(2,IELset+iel-1)-&
       !                                     DedgePosition(1,IELset+iel-1))
       !        end do
       !      end if

       ! --------------------- DOF COMBINATION PHASE ------------------------

       ! Values of all basis functions calculated. Now we can start 
       ! to integrate!

       ! Clear the local matrices
       p_Dentryii(:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryai(:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryia(:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryaa(:,:,1:IELmax-IELset+1) = 0.0_DP

       !      p_Dentryii = 0.0_DP
       !      p_Dentryai = 0.0_DP
       !      p_Dentryia = 0.0_DP
       !      p_Dentryaa = 0.0_DP

       ! We have two different versions for the integration - one
       ! with constant coefficients and one with nonconstant coefficients.
       !
       ! Check the bilinear form which one to use:

       !      if (rlocalMatrixAssembly%rform%ballCoeffConstant) then
       !      
       !        ! Constant coefficients. The coefficients are to be found in
       !        ! the Dcoefficients variable of the form.
       !        !
       !        ! Loop over the elements in the current set.
       !
       !        do iel = 1,IELmax-IELset+1
       !          
       !          ! Get the length of the edge.
       !          dlen = DedgeLength(iel)
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
       !              ! Now loop through all possible combinations of DOF`s
       !              ! in the current cubature point. The outer loop
       !              ! loops through the "O"`s in the above picture,
       !              ! the test functions:
       !
       !              do idofe = 1,indofTest
       !              
       !                ! Get the value of the (test) basis function 
       !                ! phi_i (our "O") in the cubature point:
       !                db = p_DbasTest(idofe,ib,icubp,iel)
       !                
       !                ! Perform an inner loop through the other DOF`s
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

       ! Nonconstant coefficients. The coefficients are to be found in
       ! the Dcoefficients variable as computed above.
       !
       ! Loop over the elements.

       do iel = 1,IELmax-IELset+1

          ! Get the length of the edge.
          !dlen = DedgeLength(iel)
          dlen = 0.5_DP*raddTriaData%p_Dedgelength(Iedgelist(IELset+iel-1))

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

             ! calculate the current weighting factor in the cubature formula
             ! in that cubature point.
             !            domega = dlen * p_Domega(icubp)
             domega1 = dlen * rlocalMatrixAssembly(1)%p_Domega(icubp)
             domega2 = dlen * rlocalMatrixAssembly(2)%p_Domega(icubp)

             ! Loop over the additive factors in the bilinear form.
             do ialbet = 1,rlocalMatrixAssembly(1)%rform%itermcount

                ! Get from Idescriptors the type of the derivatives for the 
                ! test and trial functions. The summand we calculate
                ! here will be added to the matrix entry:
                !
                ! a_ij  =  int_... ( psi_j )_ia  *  ( phi_i )_ib
                !
                ! -> Ix=0: function value, 
                !      =1: first derivative, ...
                !    as defined in the module 'derivative'.

                ia = rlocalMatrixAssembly(1)%rform%Idescriptors(1,ialbet)
                ib = rlocalMatrixAssembly(1)%rform%Idescriptors(2,ialbet)

                ! Multiply domega with the coefficient of the form.
                ! This gives the actual value to multiply the
                ! function value with before summing up to the integral.
                ! Get the precalculated coefficient from the coefficient array.
                daux1 = domega1 * DfluxValues(ialbet,icubp,iel)
                daux2 = domega2 * DfluxValues(ialbet,icubp,iel) * (-1.0_dp)

                ! Now loop through all possible combinations of DOF`s
                ! in the current cubature point. The outer loop
                ! loops through the "O" in the above picture,
                ! the test functions:

                do idofe = 1,indofTest

                   ! Get the value of the (test) basis function 
                   ! phi_i (our "O") in the cubature point:
                   db1 = rlocalMatrixAssembly(1)%p_DbasTest(idofe,ib,icubp,iel)
                   db2 = rlocalMatrixAssembly(2)%p_DbasTest(idofe,ib,icubp,iel)

                   ! Perform an inner loop through the other DOF`s
                   ! (the "X"). 

                   do jdofe = 1,indofTrial

                      ! Get the value of the basis function 
                      ! psi_j (our "X") in the cubature point. 
                      ! Them multiply:
                      !    db * dbas(..) * daux
                      ! ~= phi_i * psi_j * coefficient * cub.weight
                      ! Summing this up gives the integral, so the contribution
                      ! to the global matrix. 
                      !
                      ! Simply summing up db * dbas(..) * daux would give
                      ! the coefficient of the local matrix. We save this
                      ! contribution in the local matrix of element iel.

                      !JCOLB = Kentry(jdofe,idofe,iel)
                      !p_DA(JCOLB) = p_DA(JCOLB) + db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
                      !                  p_Dentry(jdofe,idofe,iel) = &
                      !                      p_Dentry(jdofe,idofe,iel)+db*p_DbasTrial(jdofe,ia,icubp,iel)*daux

                      !                  ! Testfunction on the 'first' (i) side
                      !                  p_Dentryii(jdofe,idofe,iel) = &
                      !                      p_Dentryii(jdofe,idofe,iel)+db1*rlocalMatrixAssembly(1)%p_DbasTrial(jdofe,ia,icubp,iel)*daux1*p_Dside(1,icubp,iel)   
                      !                  p_Dentryai(jdofe,idofe,iel) = &
                      !                      p_Dentryai(jdofe,idofe,iel)+db1*rlocalMatrixAssembly(2)%p_DbasTrial(jdofe,ia,icubp,iel)*daux1*p_Dside(2,icubp,iel)   
                      !                  
                      !                  ! Testfunction on the 'second' (a) side
                      !                  p_Dentryia(jdofe,idofe,iel) = &
                      !                      p_Dentryia(jdofe,idofe,iel)+db2*rlocalMatrixAssembly(1)%p_DbasTrial(jdofe,ia,icubp,iel)*daux2*p_Dside(1,icubp,iel)
                      !                  p_Dentryaa(jdofe,idofe,iel) = &
                      !                      p_Dentryaa(jdofe,idofe,iel)+db2*rlocalMatrixAssembly(2)%p_DbasTrial(jdofe,ia,icubp,iel)*daux2*p_Dside(2,icubp,iel)
                      !                


                      ! Testfunction on the 'first' (i) side
                      p_Dentryii(jdofe,idofe,iel) = &
                           p_Dentryii(jdofe,idofe,iel)+db1*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux1*p_Dside(1,ialbet,icubp,iel)
                      p_Dentryai(jdofe,idofe,iel) = &
                           p_Dentryai(jdofe,idofe,iel)+db1*rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)*daux1*p_Dside(2,ialbet,icubp,iel)

                      ! Testfunction on the 'second' (a) side
                      p_Dentryia(jdofe,idofe,iel) = &
                           p_Dentryia(jdofe,idofe,iel)+db2*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux2*p_Dside(1,ialbet,icubp,iel)

                      !                      if ((p_Dentryia(jdofe,idofe,iel)<-1000000000.0_dp).and.(IelementList(2,IELset+iel-1).ne.0)) then
                      !                write(*,*) 'Added', db2*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux2*p_Dside(1,iel)      
                      !                write(*,*) 'ia',ia
                      !                write(*,*) 'daux1',daux1
                      !                write(*,*) 'daux2',daux2
                      !                write(*,*) 'db1',db1
                      !                write(*,*) 'db2',db2
                      !                write(*,*) 'dside1',p_Dside(1,iel)
                      !                write(*,*) 'dside2',p_Dside(2,iel)
                      !                write(*,*) 'test1',rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)
                      !                write(*,*) 'test2',rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)
                      !                        pause
                      !                      end if

                      p_Dentryaa(jdofe,idofe,iel) = &
                           p_Dentryaa(jdofe,idofe,iel)+db2*rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)*daux2*p_Dside(2,ialbet,icubp,iel)

                      !                write(*,*) 'ia',ia
                      !                write(*,*) 'daux1',daux1
                      !                write(*,*) 'daux2',daux2
                      !                write(*,*) 'db1',db1
                      !                write(*,*) 'db2',db2
                      !                write(*,*) 'dside1',p_Dside(1,iel)
                      !                write(*,*) 'dside2',p_Dside(2,iel)
                      !                write(*,*) 'test1',rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)
                      !                write(*,*) 'test2',rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)
                      !                pause


                   end do

                end do ! jdofe

             end do ! ialbet

          end do ! icubp 

       end do ! iel

       !      end if ! rform%ballCoeffConstant

       ! Incorporate the local matrices into the global one.
       ! Kentry gives the position of the additive contributions in Dentry.
       !
       ! OpenMP-Extension: This is a critical section. Only one thread is
       ! allowed to write to the matrix, otherwise the matrix may get
       ! messed up.
       ! The critical section is put around both loops as indofTest/indofTrial
       ! are usually small and quickly to handle.

       !      if (cconstrType .eq. BILF_MATC_LUMPED) then
       !
       !        !$omp critical
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
       !        !$omp end critical
       !
       !      else

       !$omp critical
       do iel = 1,IELmax-IELset+1

          do idofe = 1,indofTest
             do jdofe = 1,indofTrial
                !              p_DA(p_Kentry(jdofe,idofe,iel)) = &
                !                  p_DA(p_Kentry(jdofe,idofe,iel)) + p_Dentry(jdofe,idofe,iel)

                p_DA(p_Kentryii(jdofe,idofe,iel)) = &
                     p_DA(p_Kentryii(jdofe,idofe,iel)) + p_Dentryii(jdofe,idofe,iel)



                if (IelementList(2,IELset+iel-1).ne.0) then

                   p_DA(p_Kentryia(jdofe,idofe,iel)) = &
                        p_DA(p_Kentryia(jdofe,idofe,iel)) + p_Dentryia(jdofe,idofe,iel)

                   p_DA(p_Kentryai(jdofe,idofe,iel)) = &
                        p_DA(p_Kentryai(jdofe,idofe,iel)) + p_Dentryai(jdofe,idofe,iel)!*real(min(1,IelementList(2,IELset+iel-1)))
                   p_DA(p_Kentryaa(jdofe,idofe,iel)) = &
                        p_DA(p_Kentryaa(jdofe,idofe,iel)) + p_Dentryaa(jdofe,idofe,iel)!*real(min(1,IelementList(2,IELset+iel-1)))
                end if

             end do
          end do

       end do ! iel
       !$omp end critical

       !      end if

    end do ! IELset
    !$omp end do

    ! Release the local matrix assembly structure
    call bilf_releaseAssemblyData(rlocalMatrixAssembly(1))
    call bilf_releaseAssemblyData(rlocalMatrixAssembly(2))

    ! Deallocate memory
    deallocate(Dxi2D, DpointsRef) !, DpointsPar, DedgeLength)
    deallocate(p_Kentryii,p_Kentryia,p_Kentryai,p_Kentryaa)
    deallocate(p_Dentryii,p_Dentryia,p_Dentryai,p_Dentryaa)
    deallocate(p_Dside)
    deallocate(IelementList)

    !$omp end parallel

  end subroutine dg_bilf_assembleSubmeshMat9Bdr2D









  !****************************************************************************

  !<subroutine>

  subroutine dg_bilf_initAssembly_reverseCubPoints(rmatrixAssembly,rform,celementTest,&
       celementTrial,ccubType,nelementsPerBlock)

    !<description>
    ! Initialise a matrix assembly structure for assembling a bilinear form.
    !</description>

    !<input>
    ! The bilinear form specifying the underlying PDE of the discretisation.
    type(t_bilinearForm), intent(in) :: rform

    ! Type of element in the test space.
    integer(I32), intent(in) :: celementTest

    ! Type of element in the trial space.
    integer(I32), intent(in) :: celementTrial

    ! Type of cubature formula to use.
    integer(I32), intent(in) :: ccubType

    ! Optional: Maximum number of elements to process simultaneously.
    ! If not specified, BILF_NELEMSIM is assumed.
    integer, intent(in), optional :: nelementsPerBlock
    !</input>

    !<output>
    ! A matrix assembly structure.
    type(t_bilfMatrixAssembly), intent(out) :: rmatrixAssembly
    !</output>

    !</subroutine>

    ! local variables
    real(dp) :: dtemp
    real(dp), dimension(2) :: dtemp2

    logical, dimension(EL_MAXNDER) :: BderTrialTempl, BderTestTempl
    integer :: i,i1

    ! Initialise the structure.
    rmatrixAssembly%rform = rform
    rmatrixAssembly%ccubType = ccubType
    rmatrixAssembly%nelementsPerBlock = BILF_NELEMSIM
    if (present(nelementsPerBlock)) &
         rmatrixAssembly%nelementsPerBlock = nelementsPerBlock
    rmatrixAssembly%celementTrial = celementTrial
    rmatrixAssembly%celementTest = celementTest

    ! Get the number of local DOF`s for trial and test functions
    rmatrixAssembly%indofTrial = elem_igetNDofLoc(celementTrial)
    rmatrixAssembly%indofTest = elem_igetNDofLoc(celementTest)

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
               OU_CLASS_ERROR,OU_MODE_STD,'bilf_initAssembly')
          call sys_halt()
       endif

       BderTrialTempl(I1)=.true.

       ! Then those of the test functions
       I1=rform%Idescriptors(2,I)

       if ((I1 .le.0) .or. (I1 .gt. DER_MAXNDER)) then
          call output_line ('Invalid descriptor!',&
               OU_CLASS_ERROR,OU_MODE_STD,'bilf_initAssembly')
          call sys_halt()
       endif

       BderTestTempl(I1)=.true.
    end do

    ! Determine if trial and test space is the same.
    rmatrixAssembly%bIdenticalTrialAndTest = (celementTest .eq. celementTrial)

    if (rmatrixAssembly%bIdenticalTrialAndTest) then
       ! Build the actual combination of what the element should calculate.
       rmatrixAssembly%BderTrial(:) = BderTrialTempl(:) .or. BderTestTempl(:)
       rmatrixAssembly%BderTest(:) = rmatrixAssembly%BderTrial(:)
    else
       ! Build the actual combination of what the element should calculate.
       ! Copy BDERxxxx to BDERxxxxAct
       rmatrixAssembly%BderTrial(:) = BderTrialTempl(:)
       rmatrixAssembly%BderTest(:) = BderTestTempl(:)
    end if

    ! Get the number of vertices of the element, specifying the transformation
    ! form the reference to the real element.
    rmatrixAssembly%NVE = elem_igetNVE(celementTest)
    if (rmatrixAssembly%NVE .ne. elem_igetNVE(celementTrial)) then
       call output_line ('Element spaces incompatible!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_initAssembly')
       call sys_halt()
    end if

    ! Get from the element space the type of coordinate system
    ! that is used there:
    rmatrixAssembly%ctrafoType = elem_igetTrafoType(celementTest)

    ! Get the number of cubature points for the cubature formula
    rmatrixAssembly%ncubp = cub_igetNumPts(ccubType)

    ! Allocate two arrays for the points and the weights
    allocate(rmatrixAssembly%p_Domega(rmatrixAssembly%ncubp))
    allocate(rmatrixAssembly%p_DcubPtsRef(&
         trafo_igetReferenceDimension(rmatrixAssembly%ctrafoType),&
         rmatrixAssembly%ncubp))

    ! Get the cubature formula
    call cub_getCubature(ccubType,rmatrixAssembly%p_DcubPtsRef,rmatrixAssembly%p_Domega)


    ! Revert the numbering of cubature points
    do i = 1, size(rmatrixAssembly%p_Domega)/2
       dtemp2(:) = rmatrixAssembly%p_DcubPtsRef(:,i)
       rmatrixAssembly%p_DcubPtsRef(:,i) = rmatrixAssembly%p_DcubPtsRef(:,size(rmatrixAssembly%p_Domega)+1-i)
       rmatrixAssembly%p_DcubPtsRef(:,size(rmatrixAssembly%p_Domega)+1-i) = dtemp2(:)

       dtemp = rmatrixAssembly%p_Domega(i)
       rmatrixAssembly%p_Domega(i) = rmatrixAssembly%p_Domega(size(rmatrixAssembly%p_Domega)+1-i)
       rmatrixAssembly%p_Domega(size(rmatrixAssembly%p_Domega)+1-i) = dtemp

    end do



    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    rmatrixAssembly%cevaluationTag = elem_getEvaluationTag(rmatrixAssembly%celementTest)
    rmatrixAssembly%cevaluationTag = ior(rmatrixAssembly%cevaluationTag,&
         elem_getEvaluationTag(rmatrixAssembly%celementTrial))

  end subroutine dg_bilf_initAssembly_reverseCubPoints










































  !****************************************************************************

  !<subroutine>

  subroutine bilf_dg_buildMatrixBlEdge2D (rform, ccubType, bclear, rmatrix,&
       rvectorSol, raddTriaData,&
       flux_dg_buildMatrixBlEdge2D_sim,&
       rcollection, cconstrType)

    !<description>
    ! This routine calculates the entries of a finite element matrix in 2D.
    ! The matrix structure must be prepared with bilf_createMatrixStructure
    ! in advance.
    ! In case the array for the matrix entries does not exist, the routine
    ! allocates memory in size of the matrix of the heap for the matrix entries.
    !
    ! For setting up the entries, the discretisation structure attached to
    ! the matrix is used (rmatrix%p_rdiscretisation). This is
    ! normally attached to the matrix by bilf_createMatrixStructure.
    !
    ! The matrix must be unsorted when this routine is called, 
    ! otherwise an error is thrown.
    !</description>

    !<input>
    ! The bilinear form specifying the underlying PDE of the discretisation.
    type(t_bilinearForm), intent(in) :: rform

    ! A line cubature formula CUB_xxxx_1D to be used for line integration.
    integer(I32), intent(in) :: ccubType

    ! The solution vector. Used to calculate the solution on the edges.
    type(t_vectorBlock), intent(in) :: rvectorSol

    ! Additional triangulation data
    type(t_additionalTriaData), intent(in) :: raddTriaData

    ! Whether to clear the matrix before calculating the entries.
    ! If .FALSE., the new matrix entries are added to the existing entries.
    logical, intent(in) :: bclear

    ! A callback routine for the flux function.
    include 'intf_flux_dg_buildMatrixBlEdge2D.inc'
    optional :: flux_dg_buildMatrixBlEdge2D_sim

    ! OPTIONAL: One of the BILF_MATC_xxxx constants that allow to specify
    ! the matrix construction method. If not specified,
    ! BILF_MATC_EDGEBASED is used.
    integer, intent(in), optional :: cconstrType
    !</input>

    !<inputoutput>
    ! The FE matrix. Calculated matrix entries are imposed to this matrix.
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: A collection structure. This structure is given to the
    ! callback function for nonconstant coefficients to provide additional
    ! information. 
    type(t_collection), intent(inout), target, optional :: rcollection
    !</inputoutput>

    !</subroutine>

    ! local variables
    type(t_bilfMatrixAssembly), dimension(2) :: rmatrixAssembly
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: IelementList, p_IedgeList
    integer :: ccType
    integer :: iedge, ielementDistr

    ! The matrix must be unsorted, otherwise we can not set up the matrix.
    ! Note that we cannot switch off the sorting as easy as in the case
    ! of a vector, since there is a structure behind the matrix! So the caller
    ! has to make sure, the matrix is unsorted when this routine is called.
    if (rmatrix%RmatrixBlock(1,1)%isortStrategy .gt. 0) then
       call output_line ('Matrix-structure must be unsorted!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ! The matrix must provide discretisation structures
    if ((.not. associated(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest)) .or. &
         (.not. associated(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial))) then
       call output_line ('No discretisation associated!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ! The discretisation must provide a triangulation structure
    if ((.not. associated(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation)) .or. &
         (.not. associated(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial%p_rtriangulation))) then
       call output_line('No triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    !  ! The discretisation must provide a boundary structure
    !  if ((.not. associated(rmatrix%p_rspatialDiscrTest%p_rboundary)) .or. &
    !      (.not. associated(rmatrix%p_rspatialDiscrTrial%p_rboundary))) then
    !    call output_line('No boundary associated!',&
    !        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
    !    call sys_halt()
    !  end if

    !  ! Set pointers for quicker access
    !  p_rboundary => rmatrix%p_rspatialDiscrTest%p_rboundary
    !  if (.not.associated(p_rboundary, rmatrix%p_rspatialDiscrTrial%p_rboundary)) then
    !    call output_line('Invalid boundary associated!',&
    !        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
    !    call sys_halt()
    !  end if

    p_rtriangulation => rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation
    if (.not.associated(p_rtriangulation, rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial%p_rtriangulation)) then
       call output_line('Invalid triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ccType = BILF_MATC_EDGEBASED
    if (present(cconstrType)) ccType = cconstrType

    ! Do we have a uniform triangulation? Would simplify a lot...
    select case (rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%ccomplexity)
    case (SPDISC_UNIFORM,SPDISC_CONFORMAL) 
       ! Uniform and conformal discretisations
       select case (rmatrix%RmatrixBlock(1,1)%cdataType)
       case (ST_DOUBLE) 
          ! Which matrix structure do we have?
          select case (rmatrix%RmatrixBlock(1,1)%cmatrixFormat) 
          case (LSYSSC_MATRIX9)

             ! Probably allocate/clear the matrix
             !        if (rmatrix%h_DA .eq. ST_NOHANDLE) then
             !          call lsyssc_allocEmptyMatrix(rmatrix,LSYSSC_SETM_ZERO)
             !        else
             if (bclear) call lsysbl_clearMatrix (rmatrix)
             !        end if


             ! Allocate the edgelist
             allocate(p_IedgeList(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%NMT))

             ! All edges
             forall (iedge = 1:rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%NMT) p_IedgeList(iedge)=iedge

             ! Initialise a matrix assembly structure for that element distribution
             ielementDistr = 1
             call bilf_initAssembly(rmatrixAssembly(1),rform,&
                  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(ielementDistr)%celement,&
                  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial%RelementDistr(ielementDistr)%celement,&
                  ccubType, BILF_NELEMSIM)

             ! Do the same for the other side of the egde
             ielementDistr = 1
             call dg_bilf_initAssembly_reverseCubPoints(rmatrixAssembly(2),rform,&
                  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(ielementDistr)%celement,&
                  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial%RelementDistr(ielementDistr)%celement,&
                  ccubType, BILF_NELEMSIM)

             ! Assemble the data for all elements in this element distribution
             call dg_bilf_assembleSubmeshMat9Bdr2D_Block (rmatrixAssembly, rmatrix,&
                  rvectorSol, raddTriaData, p_IedgeList(1:rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%NMT),&
                  ccType, flux_dg_buildMatrixBlEdge2D_sim, rcollection)


             ! Release the assembly structure.
             call bilf_doneAssembly(rmatrixAssembly(1))
             call bilf_doneAssembly(rmatrixAssembly(2))

             ! Deallocate the edgelist
             deallocate(p_IedgeList)



          case default
             call output_line ('Not supported matrix structure!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
             call sys_halt()
          end select

       case default
          call output_line ('Single precision matrices currently not supported!', &
               OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
          call sys_halt()
       end select

    case default
       call output_line ('General discretisation not implemented!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end select

  end subroutine bilf_dg_buildMatrixBlEdge2D



  !****************************************************************************

  !<subroutine>  

  subroutine dg_bilf_assembleSubmeshMat9Bdr2D_Block (rmatrixAssembly, rmatrix,&
       rvectorSol, raddTriaData, IedgeList,&
       cconstrType, flux_dg_buildMatrixBlEdge2D_sim, rcollection)

    !<description>

    ! Assembles the matrix entries for a submesh by integrating over the
    ! boundary region in 2D.

    !</description>

    !<input>

    ! List of elements where to assemble the bilinear form.
    integer, dimension(:), intent(in), target :: IedgeList

    ! One of the BILF_MATC_xxxx constants that allow to specify the
    ! matrix construction method.
    integer, intent(in) :: cconstrType

    ! The solution vector. Used to calculate the solution on the edges.
    type(t_vectorBlock), intent(in) :: rvectorSol

    ! Additional triangulation data
    type(t_additionalTriaData), intent(in) :: raddTriaData

    ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
    ! Must be present if the matrix has nonconstant coefficients!
    include 'intf_flux_dg_buildMatrixBlEdge2D.inc'
    optional :: flux_dg_buildMatrixBlEdge2D_sim

    !</input>

    !<inputoutput>

    ! A matrix assembly structure prepared with bilf_initAssembly.
    type(t_bilfMatrixAssembly), intent(inout), dimension(2), target :: rmatrixAssembly

    ! A matrix where to assemble the contributions to.
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: A pointer to a collection structure. This structure is given to the
    ! callback function for nonconstant coefficients to provide additional
    ! information. 
    type(t_collection), intent(inout), target, optional :: rcollection

    !</inputoutput>

    !</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_DA
    integer :: indofTest,indofTrial,ncubp

    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax,ibdc,k
    integer :: iel,icubp,ialbet,ia,ib,idofe,jdofe,nve
    real(DP) :: domega1,domega2,db1,db2,dlen
    real(dp), dimension(:,:), allocatable :: daux1, daux2
    integer(I32) :: cevaluationTag
    type(t_bilfMatrixAssembly), dimension(2), target :: rlocalMatrixAssembly
    type(t_domainIntSubset), dimension(2) :: rintSubset
    integer, dimension(:,:,:), pointer :: p_Kentryii, p_Kentryia, p_Kentryai, p_Kentryaa
    real(DP), dimension(:,:,:,:,:), pointer :: p_Dentryii, p_Dentryia, p_Dentryai, p_Dentryaa
    real(DP), dimension(:,:,:), pointer :: p_Dcoords
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:,:,:), pointer :: p_Dside
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial
    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    real(DP), dimension(:), pointer :: p_DcoefficientsBilf
    integer, dimension(:,:), pointer :: p_IdofsTest
    integer, dimension(:,:), pointer :: p_IdofsTrial
    type(t_evalElementSet), pointer :: p_revalElementSet
    integer, dimension(:,:),pointer :: p_Idescriptors
    integer, dimension(:,:), allocatable, target :: IelementList
    integer :: nvar

    ! Pointer to Ielementsatedge in the triangulation
    integer, dimension(:,:), pointer :: p_IelementsAtEdge

    ! Pointer to IverticesAtEdge in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtEdge

    ! Pointer to the vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Pointer to IverticesAtEelement in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Space for the values of the flux function
    real(DP), dimension(:,:,:,:,:), allocatable :: DfluxValues

    ! Arrays for cubature points 1D->2D
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D_1, Dxi1D_2
    real(DP), dimension(:,:,:,:), allocatable :: Dxi2D,DpointsRef
    real(DP), dimension(:,:), allocatable :: DpointsPar
    real(DP), dimension(:), allocatable :: DedgeLength

    integer(i32) :: icoordSystem
    integer :: NEL
    integer :: iside
    logical :: bisLinearTrafo

    integer :: iblock,jblock

    !    type t_array
    !     ! Pointer to the double-valued matrix or vector data
    !     real(DP), dimension(:), pointer :: Da
    !    end type t_array

    !    type(t_array), dimension(:,:), allocatable :: p_matrixBlockDataPointers

    !    ! Boundary component?
    !    ibdc = rboundaryRegion%iboundCompIdx

    ! Get some pointers for faster access
    ! call lsyssc_getbase_double (rmatrix,p_DA)
    indofTest = rmatrixAssembly(1)%indofTest
    indofTrial = rmatrixAssembly(1)%indofTrial
    ncubp = rmatrixAssembly(1)%ncubp
    nvar = rvectorSol%nblocks

    ! Open-MP-Extension: Copy the matrix assembly data to the local
    ! matrix assembly data, where we can allocate memory.
    !
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the matrix
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original matrix assembly
    ! stucture or disturbing the data of the other processors.
    !
    !$omp parallel default(shared) &
    !$omp private(DedgeLength,DpointsPar,DpointsRef,Dxi1D,Dxi2D,IELmax,bisLinearTrafo,&
    !$omp         cevaluationTag,daux,db,dlen,domega,ia,ialbet,ib,icoordSystem,icubp,&
    !$omp         idofe,iel,jdofe,k,p_DbasTest,p_DbasTrial,p_Dcoefficients,&
    !$omp         p_DcoefficientsBilf,p_Dcoords,p_DcubPtsRef,p_Dentry,p_Domega,&
    !$omp         p_Idescriptors,p_IdofsTest,p_IdofsTrial,p_Kentry,&
    !$omp         p_revalElementSet,rintSubset,rlocalMatrixAssembly)
    rlocalMatrixAssembly(1) = rmatrixAssembly(1)
    rlocalMatrixAssembly(2) = rmatrixAssembly(2)
    call bilf_allocAssemblyData(rlocalMatrixAssembly(1))
    call bilf_allocAssemblyData(rlocalMatrixAssembly(2))

    ! Allocate space for the positions of the DOFs in the matrix
    allocate(p_Kentryii(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryia(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryai(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryaa(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))

    ! Allocate auxiliary vectors
    allocate(daux1(nvar,nvar),daux2(nvar,nvar))

    ! Allocate space for the coefficient of the solutions DOFs on each side of the edge
    allocate(p_Dside(nvar,nvar,2,rmatrixAssembly(1)%rform%itermCount,ncubp,rmatrixAssembly(1)%nelementsPerBlock))

    ! Allocate space for the entries in the local matrices
    allocate(p_Dentryii(nvar,nvar,rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryia(nvar,nvar,rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryai(nvar,nvar,rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryaa(nvar,nvar,rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))



    ! Allocate space for the flux variables DIM(nvar,nvar,ialbet,ncubp,elementsperblock)
    allocate(DfluxValues(nvar,nvar,rmatrixAssembly(1)%rform%itermCount,ncubp,rlocalMatrixAssembly(1)%nelementsPerBlock))

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


    ! Get number of elements
    NEL = rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%NEL

    ! Get pointers to elements at edge
    call storage_getbase_int2D(&
         rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%h_IelementsAtEdge,&
         p_IelementsAtEdge)

    ! Get pointers to the vertex coordinates
    call storage_getbase_double2D(&
         rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%h_DvertexCoords,&
         p_DvertexCoords)

    ! Get pointers to vertices at edge
    call storage_getbase_int2D(&
         rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtEdge,&
         p_IverticesAtEdge)

    ! Get pointers to vertices at elements
    call storage_getbase_int2D(&
         rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtElement,&
         p_IverticesAtElement)   

    ! Get the elements adjacent to the given edges
    allocate(IelementList(3,size(IedgeList)))
    IelementList(1:2,1:size(IedgeList))=p_IelementsAtEdge(1:2,IedgeList(:))

    ! Copy the second component and replace 0s by 1s
    IelementList(3,size(IedgeList))=IelementList(2,size(IedgeList))
    do iel = 1,size(IedgeList)
       IelementList(3,iel)=max(IelementList(2,iel),1)
    end do

    ! Transpose the coordinate array such that we get coordinates we
    ! can work with in the mapping between 1D and 2D.
    do k = 1, ubound(rlocalmatrixAssembly(1)%p_DcubPtsRef,1)
       do icubp = 1,ncubp
          Dxi1D_1(icubp,k) = rlocalmatrixAssembly(1)%p_DcubPtsRef(k,icubp)
          Dxi1D_2(icubp,k) = rlocalmatrixAssembly(2)%p_DcubPtsRef(k,icubp)
       end do
    end do

    ! Allocate memory for the cubature points in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,2,rlocalMatrixAssembly(1)%nelementsPerBlock))

    ! Allocate memory for the coordinates of the reference points
    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalMatrixAssembly(1)%nelementsPerBlock,2))

    !    ! Allocate memory for the parameter values of the points on the boundary
    !    allocate(DpointsPar(ncubp,rlocalMatrixAssembly%nelementsPerBlock))
    !
    !    ! Allocate memory for the length of edges on the boundary
    !    allocate(DedgeLength(rlocalMatrixAssembly%nelementsPerBlock))

    ! Get the type of coordinate system
    icoordSystem = elem_igetCoordSystem(rlocalMatrixAssembly(1)%celementTrial)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !$omp do schedule(static,1)
    do IELset = 1, size(IedgeList), rlocalMatrixAssembly(1)%nelementsPerBlock

       ! We always handle nelementsPerBlock elements simultaneously.
       ! How many elements have we actually here?
       ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
       ! elements simultaneously.

       IELmax = min(size(IedgeList),IELset-1+rlocalMatrixAssembly(1)%nelementsPerBlock)

       ! Map the 1D cubature points to the edges in 2D.
       do iel = 1,IELmax-IELset+1
          call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(1,Iedgelist(IELset+iel-1)), &
               ncubp, Dxi1D_1, Dxi2D(:,:,1,iel))
          call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(2,Iedgelist(IELset+iel-1)), &
               ncubp, Dxi1D_2, Dxi2D(:,:,2,iel))
       end do

       !      ! Calculate the parameter values of the points
       !      do iel = 1,IELmax-IELset+1
       !        do icubp = 1,ncubp
       !          ! Dxi1D is in [-1,1] while the current edge has parmeter values
       !          ! [DedgePosition(1),DedgePosition(2)]. So do a linear
       !          ! transformation to transform Dxi1D into that interval, this 
       !          ! gives the parameter values in length parametrisation
       !          call mprim_linearRescale(Dxi1D(icubp,1), -1.0_DP, 1.0_DP,&
       !              DedgePosition(1,IELset+iel-1), DedgePosition(2,IELset+iel-1),&
       !              DpointsPar(icubp,iel))
       !        end do
       !      end do

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
       ! global DOF`s of our BILF_NELEMSIM elements simultaneously.
       !      call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTest, &
       !          IelementList(IELset:IELmax), p_IdofsTest)
       call dof_locGlobMapping_mult( rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest, &
            IelementList(1,IELset:IELmax), rlocalMatrixAssembly(1)%p_IdofsTest)
       call dof_locGlobMapping_mult( rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest, &
            IelementList(3,IELset:IELmax), rlocalMatrixAssembly(2)%p_IdofsTest)


       !      ! If the DOF`s for the trial functions are different, calculate them, too.
       !      if (.not. rlocalMatrixAssembly%bIdenticalTrialAndTest) then
       !        call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTrial, &
       !            IelementList(IELset:IELmax), p_IdofsTrial)
       !      end if

       ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------

       ! For the assembly of the global matrix, we use a "local"
       ! approach. At first we build a "local" system matrix according
       ! to the current element. This contains all additive
       ! contributions of element iel, which are later added at the
       ! right positions to the elements in the global system matrix.
       !
       ! We have indofTrial trial DOF`s per element and
       ! indofTest test DOF`s per element. Therefore there are
       ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
       ! "active" (i.e. have common support) on our current element, each 
       ! giving an additive contribution to the system matrix.
       !
       ! We build a quadratic indofTrial*indofTest local matrix:
       ! Kentry(1..indofTrial,1..indofTest) receives the position 
       ! in the global system matrix, where the corresponding value 
       ! has to be added to.
       ! (The corresponding contributions can be saved separately, 
       ! but we directly add them to the global matrix in this 
       ! approach.)
       !
       ! We build local matrices for all our elements 
       ! in the set simultaneously. Get the positions of the local matrices
       ! in the global matrix.
       !      call bilf_getLocalMatrixIndices (rmatrix,p_IdofsTest,p_IdofsTrial,p_Kentry,&
       !          ubound(p_IdofsTest,1),ubound(p_IdofsTrial,1),IELmax-IELset+1)  
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
       !            rlocalMatrixAssembly(1)%p_IdofsTrial, p_Kentryii,&
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
       !            rlocalMatrixAssembly(2)%p_IdofsTrial, p_Kentryia,&
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
       !            rlocalMatrixAssembly(1)%p_IdofsTrial, p_Kentryai,&
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
       !            rlocalMatrixAssembly(2)%p_IdofsTrial, p_Kentryaa,&
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTrial,1), IELmax-IELset+1)



       call bilf_getLocalMatrixIndices (rmatrix%RmatrixBlock(1,1),rlocalMatrixAssembly(1)%p_IdofsTest, &
            rlocalMatrixAssembly(1)%p_IdofsTest, p_Kentryii,&
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix%RmatrixBlock(1,1),rlocalMatrixAssembly(1)%p_IdofsTest, &
            rlocalMatrixAssembly(2)%p_IdofsTest, p_Kentryai,&
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix%RmatrixBlock(1,1),rlocalMatrixAssembly(2)%p_IdofsTest, &
            rlocalMatrixAssembly(1)%p_IdofsTest, p_Kentryia,&
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix%RmatrixBlock(1,1),rlocalMatrixAssembly(2)%p_IdofsTest, &
            rlocalMatrixAssembly(2)%p_IdofsTest, p_Kentryaa,&
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), IELmax-IELset+1)   

       ! -------------------- ELEMENT EVALUATION PHASE ----------------------

       ! Ok, we found the positions of the local matrix entries
       ! that we have to change.
       ! To calculate the matrix contributions, we have to evaluate
       ! the elements to give us the values of the basis functions
       ! in all the DOF`s in all the elements in our set.

       ! Get the element evaluation tag of all FE spaces. We need it to evaluate
       ! the elements later. All of them can be combined with OR, what will give
       ! a combined evaluation tag. 
       cevaluationTag = rlocalMatrixAssembly(1)%cevaluationTag

       ! The cubature points are already initialised by 1D->2D mapping.
       cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

       !      ! Do we have a (multi-)linear transformation?
       !      bisLinearTrafo = trafo_isLinearTrafo(rlocalMatrixAssembly%ctrafoType)
       !
       !      if (bisLinearTrafo) then
       !        ! We need the vertices of the element corners and the number
       !        ! of vertices per element to compute the length of the element
       !        ! edge at the boundary
       !        cevaluationTag = ior(cevaluationTag, EL_EVLTAG_COORDS)
       !        nve = trafo_igetNVE(rlocalMatrixAssembly%ctrafoType)
       !      end if

       ! Calculate all information that is necessary to evaluate the finite element
       ! on all cells of our subset. This includes the coordinates of the points
       ! on the cells.
       !      call elprep_prepareSetForEvaluation (p_revalElementSet,&
       !          cevaluationTag, rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
       !          IelementList(IELset:IELmax), rlocalMatrixAssembly%ctrafoType, &
       !          DpointsRef=DpointsRef)
       !      p_Dcoords => p_revalElementSet%p_Dcoords

       call elprep_prepareSetForEvaluation (&
            rlocalMatrixAssembly(1)%revalElementSet,&
            cevaluationTag,  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation, &
            IelementList(1,IELset:IELmax), rlocalMatrixAssembly(1)%ctrafoType, &
            DpointsRef=DpointsRef(:,:,:,1))
       call elprep_prepareSetForEvaluation (&
            rlocalMatrixAssembly(2)%revalElementSet,&
            cevaluationTag,  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation, &
            IelementList(3,IELset:IELmax), rlocalMatrixAssembly(2)%ctrafoType, &
            DpointsRef=DpointsRef(:,:,:,2))




       !      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
       !      if (.not. rlocalMatrixAssembly%rform%ballCoeffConstant) then
       !        if (present(fcoeff_buildMatrixScBdr2D_sim)) then
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

       ! If the flux function needs other, than just the function values from the solution
       ! (for example the derivatives), we will give an evalElementSet to it
       ! This is filled here

       call domint_initIntegrationByEvalSet (rlocalMatrixAssembly(1)%revalElementSet, rintSubset(1))
       call domint_initIntegrationByEvalSet (rlocalMatrixAssembly(2)%revalElementSet, rintSubset(2))
       !rintSubset(1)%ielementDistribution = 0
       rintSubset(1)%ielementStartIdx = IELset
       rintSubset(1)%p_Ielements => IelementList(1,IELset:IELmax)
       rintSubset(1)%p_IdofsTrial => rlocalMatrixAssembly(1)%p_IdofsTest
       rintSubset(1)%celement = rlocalMatrixAssembly(1)%celementTest
       !rintSubset(2)%ielementDistribution = 0
       rintSubset(2)%ielementStartIdx = IELset
       rintSubset(2)%p_Ielements => IelementList(2,IELset:IELmax)
       rintSubset(2)%p_IdofsTrial => rlocalMatrixAssembly(2)%p_IdofsTest
       rintSubset(2)%celement = rlocalMatrixAssembly(2)%celementTest


       call flux_dg_buildMatrixBlEdge2D_sim (&
            !            rlocalVectorAssembly(1)%p_Dcoefficients(1,:,1:IELmax-IELset+1),&
            !            DsolVals(:,:,1:IELmax-IELset+1),&
       DfluxValues(:,:,:,:,1:IELmax-IELset+1),&
            rvectorSol,&
            IelementList(2,IELset:IELmax),&
            p_Dside,&
            raddTriaData%p_Dnormals(:,Iedgelist(IELset:IELmax)),&
            !DpointsReal(1:ndim2d,1:ncubp,1:IELmax-IELset+1),&
       rintSubset,&
            rcollection )


       call domint_doneIntegration (rintSubset(1))
       call domint_doneIntegration (rintSubset(2))





       ! Calculate the values of the basis functions.
       !      call elem_generic_sim2 (rlocalMatrixAssembly%celementTest, &
       !          p_revalElementSet, rlocalMatrixAssembly%BderTest, &
       !          rlocalMatrixAssembly%p_DbasTest) 
       call elem_generic_sim2 (rlocalMatrixAssembly(1)%celementTest, &
            rlocalMatrixAssembly(1)%revalElementSet,&
            rlocalMatrixAssembly(1)%BderTest, &
            rlocalMatrixAssembly(1)%p_DbasTest)
       call elem_generic_sim2 (rlocalMatrixAssembly(2)%celementTest, &
            rlocalMatrixAssembly(2)%revalElementSet,&
            rlocalMatrixAssembly(2)%BderTest, &
            rlocalMatrixAssembly(2)%p_DbasTest)

       !      ! Omit the calculation of the trial function values if they
       !      ! are identical to the test function values.
       !      if (.not. rlocalMatrixAssembly%bidenticalTrialAndTest) then
       !        call elem_generic_sim2 (rlocalMatrixAssembly%celementTrial, &
       !            p_revalElementSet, rlocalMatrixAssembly%BderTrial, &
       !            rlocalMatrixAssembly%p_DbasTrial)
       !      end if

       !      ! Calculate the length of egdes on the boundary. Depending on
       !      ! whether the transformation is (multi-)linear or not we compute
       !      ! the edge length as the distance between the two corner
       !      ! vertices of the element located on the boundary or as the real
       !      ! length of the boundary segment of the element.
       !      !
       !      ! The length of the current edge serves as a "determinant" in
       !      ! the cubature, so we have to divide it by 2 as an edge on the
       !      ! unit interval [-1,1] has length 2.
       !      if (bisLinearTrafo) then
       !        do iel = 1,IELmax-IELset+1
       !          DedgeLength(iel) = 0.5_DP*sqrt(&
       !              (p_Dcoords(1,    IelementOrientation(IELset+iel-1),iel)-&
       !               p_Dcoords(1,mod(IelementOrientation(IELset+iel-1),nve)+1,iel))**2+&
       !              (p_Dcoords(2,    IelementOrientation(IELset+iel-1),iel)-&
       !               p_Dcoords(2,mod(IelementOrientation(IELset+iel-1),nve)+1,iel))**2)
       !        end do
       !      else
       !        do iel = 1,IELmax-IELset+1
       !          DedgeLength(iel) = 0.5_DP*(DedgePosition(2,IELset+iel-1)-&
       !                                     DedgePosition(1,IELset+iel-1))
       !        end do
       !      end if

       ! --------------------- DOF COMBINATION PHASE ------------------------

       ! Values of all basis functions calculated. Now we can start 
       ! to integrate!

       ! Clear the local matrices
       p_Dentryii(:,:,:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryai(:,:,:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryia(:,:,:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryaa(:,:,:,:,1:IELmax-IELset+1) = 0.0_DP

       !      p_Dentryii = 0.0_DP
       !      p_Dentryai = 0.0_DP
       !      p_Dentryia = 0.0_DP
       !      p_Dentryaa = 0.0_DP

       ! We have two different versions for the integration - one
       ! with constant coefficients and one with nonconstant coefficients.
       !
       ! Check the bilinear form which one to use:

       !      if (rlocalMatrixAssembly%rform%ballCoeffConstant) then
       !      
       !        ! Constant coefficients. The coefficients are to be found in
       !        ! the Dcoefficients variable of the form.
       !        !
       !        ! Loop over the elements in the current set.
       !
       !        do iel = 1,IELmax-IELset+1
       !          
       !          ! Get the length of the edge.
       !          dlen = DedgeLength(iel)
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
       !              ! Now loop through all possible combinations of DOF`s
       !              ! in the current cubature point. The outer loop
       !              ! loops through the "O"`s in the above picture,
       !              ! the test functions:
       !
       !              do idofe = 1,indofTest
       !              
       !                ! Get the value of the (test) basis function 
       !                ! phi_i (our "O") in the cubature point:
       !                db = p_DbasTest(idofe,ib,icubp,iel)
       !                
       !                ! Perform an inner loop through the other DOF`s
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

       ! Nonconstant coefficients. The coefficients are to be found in
       ! the Dcoefficients variable as computed above.
       !
       ! Loop over the elements.

       do iel = 1,IELmax-IELset+1

          ! Get the length of the edge.
          !dlen = DedgeLength(iel)
          dlen = 0.5_DP*raddTriaData%p_Dedgelength(Iedgelist(IELset+iel-1))

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

             ! calculate the current weighting factor in the cubature formula
             ! in that cubature point.
             !            domega = dlen * p_Domega(icubp)
             domega1 = dlen * rlocalMatrixAssembly(1)%p_Domega(icubp)
             domega2 = dlen * rlocalMatrixAssembly(2)%p_Domega(icubp)

             ! Loop over the additive factors in the bilinear form.
             do ialbet = 1,rlocalMatrixAssembly(1)%rform%itermcount

                ! Get from Idescriptors the type of the derivatives for the 
                ! test and trial functions. The summand we calculate
                ! here will be added to the matrix entry:
                !
                ! a_ij  =  int_... ( psi_j )_ia  *  ( phi_i )_ib
                !
                ! -> Ix=0: function value, 
                !      =1: first derivative, ...
                !    as defined in the module 'derivative'.

                ia = rlocalMatrixAssembly(1)%rform%Idescriptors(1,ialbet)
                ib = rlocalMatrixAssembly(1)%rform%Idescriptors(2,ialbet)

                ! Multiply domega with the coefficient of the form.
                ! This gives the actual value to multiply the
                ! function value with before summing up to the integral.
                ! Get the precalculated coefficient from the coefficient array.
                daux1 = domega1 * DfluxValues(:,:,ialbet,icubp,iel)
                daux2 = domega2 * DfluxValues(:,:,ialbet,icubp,iel) * (-1.0_dp)

                ! Now loop through all possible combinations of DOF`s
                ! in the current cubature point. The outer loop
                ! loops through the "O" in the above picture,
                ! the test functions:

                do idofe = 1,indofTest

                   ! Get the value of the (test) basis function 
                   ! phi_i (our "O") in the cubature point:
                   db1 = rlocalMatrixAssembly(1)%p_DbasTest(idofe,ib,icubp,iel)
                   db2 = rlocalMatrixAssembly(2)%p_DbasTest(idofe,ib,icubp,iel)

                   ! Perform an inner loop through the other DOF`s
                   ! (the "X"). 

                   do jdofe = 1,indofTrial

                      ! Get the value of the basis function 
                      ! psi_j (our "X") in the cubature point. 
                      ! Them multiply:
                      !    db * dbas(..) * daux
                      ! ~= phi_i * psi_j * coefficient * cub.weight
                      ! Summing this up gives the integral, so the contribution
                      ! to the global matrix. 
                      !
                      ! Simply summing up db * dbas(..) * daux would give
                      ! the coefficient of the local matrix. We save this
                      ! contribution in the local matrix of element iel.

                      !JCOLB = Kentry(jdofe,idofe,iel)
                      !p_DA(JCOLB) = p_DA(JCOLB) + db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
                      !                  p_Dentry(jdofe,idofe,iel) = &
                      !                      p_Dentry(jdofe,idofe,iel)+db*p_DbasTrial(jdofe,ia,icubp,iel)*daux

                      !                  ! Testfunction on the 'first' (i) side
                      !                  p_Dentryii(jdofe,idofe,iel) = &
                      !                      p_Dentryii(jdofe,idofe,iel)+db1*rlocalMatrixAssembly(1)%p_DbasTrial(jdofe,ia,icubp,iel)*daux1*p_Dside(1,icubp,iel)   
                      !                  p_Dentryai(jdofe,idofe,iel) = &
                      !                      p_Dentryai(jdofe,idofe,iel)+db1*rlocalMatrixAssembly(2)%p_DbasTrial(jdofe,ia,icubp,iel)*daux1*p_Dside(2,icubp,iel)   
                      !                  
                      !                  ! Testfunction on the 'second' (a) side
                      !                  p_Dentryia(jdofe,idofe,iel) = &
                      !                      p_Dentryia(jdofe,idofe,iel)+db2*rlocalMatrixAssembly(1)%p_DbasTrial(jdofe,ia,icubp,iel)*daux2*p_Dside(1,icubp,iel)
                      !                  p_Dentryaa(jdofe,idofe,iel) = &
                      !                      p_Dentryaa(jdofe,idofe,iel)+db2*rlocalMatrixAssembly(2)%p_DbasTrial(jdofe,ia,icubp,iel)*daux2*p_Dside(2,icubp,iel)
                      !                


                      ! Testfunction on the 'first' (i) side
                      do iblock = 1, nvar
                         do jblock = 1, nvar
                            p_Dentryii(iblock,jblock,jdofe,idofe,iel) = &
                                 p_Dentryii(iblock,jblock,jdofe,idofe,iel)+db1*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux1(iblock,jblock)*p_Dside(iblock,jblock,1,ialbet,icubp,iel)
                            p_Dentryai(iblock,jblock,jdofe,idofe,iel) = &
                                 p_Dentryai(iblock,jblock,jdofe,idofe,iel)+db1*rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)*daux1(iblock,jblock)*p_Dside(iblock,jblock,2,ialbet,icubp,iel)

                            ! Testfunction on the 'second' (a) side
                            p_Dentryia(iblock,jblock,jdofe,idofe,iel) = &
                                 p_Dentryia(iblock,jblock,jdofe,idofe,iel)+db2*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux2(iblock,jblock)*p_Dside(iblock,jblock,1,ialbet,icubp,iel)

                            !                      if ((p_Dentryia(jdofe,idofe,iel)<-1000000000.0_dp).and.(IelementList(2,IELset+iel-1).ne.0)) then
                            !                write(*,*) 'Added', db2*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux2*p_Dside(1,iel)      
                            !                write(*,*) 'ia',ia
                            !                write(*,*) 'daux1',daux1
                            !                write(*,*) 'daux2',daux2
                            !                write(*,*) 'db1',db1
                            !                write(*,*) 'db2',db2
                            !                write(*,*) 'dside1',p_Dside(1,iel)
                            !                write(*,*) 'dside2',p_Dside(2,iel)
                            !                write(*,*) 'test1',rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)
                            !                write(*,*) 'test2',rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)
                            !                        pause
                            !                      end if

                            p_Dentryaa(iblock,jblock,jdofe,idofe,iel) = &
                                 p_Dentryaa(iblock,jblock,jdofe,idofe,iel)+db2*rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)*daux2(iblock,jblock)*p_Dside(iblock,jblock,2,ialbet,icubp,iel)
                         end do
                      end do
                      !                write(*,*) 'ia',ia
                      !                write(*,*) 'daux1',daux1
                      !                write(*,*) 'daux2',daux2
                      !                write(*,*) 'db1',db1
                      !                write(*,*) 'db2',db2
                      !                write(*,*) 'dside1',p_Dside(1,iel)
                      !                write(*,*) 'dside2',p_Dside(2,iel)
                      !                write(*,*) 'test1',rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)
                      !                write(*,*) 'test2',rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)
                      !                pause

                   end do ! idofe

                end do ! jdofe

             end do ! ialbet

          end do ! icubp 

       end do ! iel

       !      end if ! rform%ballCoeffConstant

       ! Incorporate the local matrices into the global one.
       ! Kentry gives the position of the additive contributions in Dentry.
       !
       ! OpenMP-Extension: This is a critical section. Only one thread is
       ! allowed to write to the matrix, otherwise the matrix may get
       ! messed up.
       ! The critical section is put around both loops as indofTest/indofTrial
       ! are usually small and quickly to handle.

       !      if (cconstrType .eq. BILF_MATC_LUMPED) then
       !
       !        !$omp critical
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
       !        !$omp end critical
       !
       !      else




       !      call lsysbl_getbase_double(rmatrix,p_Da)

       !      ! Get pointers to the data entries of the block matrix
       !      do iblock = 1, nvar
       !        do jblock = 1,nvar
       !          call lsyssc_getbase_double(rmatrix%RmatrixBlock(iblock,jblock),p_matrixBlockDataPointers(iblock,jblock)%Da)
       !        end do
       !      end do

       !p_Ddata(rvector%RvectorBlock(ivar)%iidxFirstEntry+rlocalVectorAssembly(1)%p_Idofs(idofe,iel)-1) = &            
       !              p_Ddata(rvector%RvectorBlock(ivar)%iidxFirstEntry+rlocalVectorAssembly(1)%p_Idofs(idofe,iel)-1) + &
       !              DlocalData(ivar,1,idofe)






       !$omp critical
       do iblock = 1, nvar
          do jblock = 1, nvar
             call lsyssc_getbase_double(rmatrix%RmatrixBlock(iblock,jblock),p_Da)
             do iel = 1,IELmax-IELset+1

                do idofe = 1,indofTest
                   do jdofe = 1,indofTrial
                      !              p_DA(p_Kentry(jdofe,idofe,iel)) = &
                      !                  p_DA(p_Kentry(jdofe,idofe,iel)) + p_Dentry(jdofe,idofe,iel)

                      !                p_matrixBlockDataPointers(iblock,jblock)%Da(p_Kentryii(jdofe,idofe,iel)) = &
                      !                  p_matrixBlockDataPointers(iblock,jblock)%Da(p_Kentryii(jdofe,idofe,iel)) + p_Dentryii(iblock,jblock,jdofe,idofe,iel)
                      !                
                      !                if (IelementList(2,IELset+iel-1).ne.0) then
                      !                
                      !                p_matrixBlockDataPointers(iblock,jblock)%Da(p_Kentryia(jdofe,idofe,iel)) = &
                      !                  p_matrixBlockDataPointers(iblock,jblock)%Da(p_Kentryia(jdofe,idofe,iel)) + p_Dentryia(iblock,jblock,jdofe,idofe,iel)
                      !                
                      !                p_matrixBlockDataPointers(iblock,jblock)%Da(p_Kentryai(jdofe,idofe,iel)) = &
                      !                  p_matrixBlockDataPointers(iblock,jblock)%Da(p_Kentryai(jdofe,idofe,iel)) + p_Dentryai(iblock,jblock,jdofe,idofe,iel)!*real(min(1,IelementList(2,IELset+iel-1)))
                      !                p_matrixBlockDataPointers(iblock,jblock)%Da(p_Kentryaa(jdofe,idofe,iel)) = &
                      !                  p_matrixBlockDataPointers(iblock,jblock)%Da(p_Kentryaa(jdofe,idofe,iel)) + p_Dentryaa(iblock,jblock,jdofe,idofe,iel)!*real(min(1,IelementList(2,IELset+iel-1)))
                      !                end if  
                      p_Da(p_Kentryii(jdofe,idofe,iel)) = &
                           p_Da(p_Kentryii(jdofe,idofe,iel)) + p_Dentryii(iblock,jblock,jdofe,idofe,iel)

                      if (IelementList(2,IELset+iel-1).ne.0) then

                         p_Da(p_Kentryia(jdofe,idofe,iel)) = &
                              p_Da(p_Kentryia(jdofe,idofe,iel)) + p_Dentryia(iblock,jblock,jdofe,idofe,iel)

                         p_Da(p_Kentryai(jdofe,idofe,iel)) = &
                              p_Da(p_Kentryai(jdofe,idofe,iel)) + p_Dentryai(iblock,jblock,jdofe,idofe,iel)!*real(min(1,IelementList(2,IELset+iel-1)))
                         p_Da(p_Kentryaa(jdofe,idofe,iel)) = &
                              p_Da(p_Kentryaa(jdofe,idofe,iel)) + p_Dentryaa(iblock,jblock,jdofe,idofe,iel)!*real(min(1,IelementList(2,IELset+iel-1)))
                      end if

                   end do
                end do

             end do ! iel
          end do ! jblock
       end do ! iblock
       !$omp end critical

       !      end if

    end do ! IELset
    !$omp end do

    ! Release the local matrix assembly structure
    call bilf_releaseAssemblyData(rlocalMatrixAssembly(1))
    call bilf_releaseAssemblyData(rlocalMatrixAssembly(2))

    ! Deallocate memory
    deallocate(Dxi2D, DpointsRef) !, DpointsPar, DedgeLength)
    deallocate(p_Kentryii,p_Kentryia,p_Kentryai,p_Kentryaa)
    deallocate(p_Dentryii,p_Dentryia,p_Dentryai,p_Dentryaa)
    deallocate(p_Dside)
    deallocate(IelementList)
    deallocate(daux1,daux2)

    !$omp end parallel

  end subroutine dg_bilf_assembleSubmeshMat9Bdr2D_Block




















































  !****************************************************************************

  !<subroutine>

  subroutine bilf_buildMatrixBlock2 (rform, bclear, rmatrix,&
       fcoeff_buildMatrixBl_sim,rcollection,rscalarAssemblyInfo)
    use extstdassemblyinfo
    !<description>
    ! This routine calculates the entries of a finite element matrix.
    ! The matrix structure must be prepared with bilf_createMatrixStructure
    ! in advance.
    ! In case the array for the matrix entries does not exist, the routine
    ! allocates memory in size of the matrix of the heap for the matrix entries.
    !
    ! For setting up the entries, the discretisation structure attached to
    ! the matrix is used (rmatrix%p_rdiscretisation). This is
    ! normally attached to the matrix by bilf_createMatrixStructure.
    !
    ! The matrix must be unsorted when this routine is called, 
    ! otherwise an error is thrown.
    !
    ! IMPLEMENTATIONAL REMARK:
    ! This is a new implementation of the matrix assembly using element subsets.
    ! In contrast to bilf_buildMatrixScalar, this routine loops itself about
    ! the element subsets and calls bilf_initAssembly/
    ! bilf_assembleSubmeshMatrix9/bilf_doneAssembly to assemble matrix
    ! contributions of a submesh.
    ! The bilf_assembleSubmeshMatrix9 interface allows to assemble parts of a
    ! matrix based on an arbitrary element list which is not bound to an
    ! element distribution.
    !</description>

    !<input>
    ! The bilinear form specifying the underlying PDE of the discretisation.
    type(t_bilinearForm), intent(in) :: rform

    ! Whether to clear the matrix before calculating the entries.
    ! If .FALSE., the new matrix entries are added to the existing entries.
    logical, intent(in) :: bclear

    ! OPTIONAL: A collection structure. This structure is given to the
    ! callback function for nonconstant coefficients to provide additional
    ! information. 
    type(t_collection), intent(inout), target, optional :: rcollection

    ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
    ! Must be present if the matrix has nonconstant coefficients!
    include 'intf_coefficientMatrixBl.inc'
    optional :: fcoeff_buildMatrixBl_sim

    ! OPTIONAL: A scalar assembly structure that gives additional information
    ! about how to set up the matrix (e.g. cubature formula). If not specified,
    ! default settings are used.
    type(t_extScalarAssemblyInfo), intent(in), optional, target :: rscalarAssemblyInfo
    !</input>

    !<inputoutput>
    ! The FE matrix. Calculated matrix entries are imposed to this matrix.
    type(t_matrixBlock), intent(inout) :: rmatrix
    !</inputoutput>

    !</subroutine>

    ! local variables
    type(t_bilfMatrixAssembly) :: rmatrixAssembly
    integer :: ielementDistr,iinfoBlock
    integer, dimension(:), pointer :: p_IelementList
    type(t_extScalarAssemblyInfo), target :: rlocalScalarAssemblyInfo
    type(t_extScalarAssemblyInfo), pointer :: p_rscalarAssemblyInfo


    ! The matrix must be unsorted, otherwise we can not set up the matrix.
    ! Note that we cannot switch off the sorting as easy as in the case
    ! of a vector, since there is a structure behind the matrix! So the caller
    ! has to make sure, the matrix is unsorted when this routine is called.
    !  if (rmatrix%isortStrategy .gt. 0) then
    !    call output_line ('Matrix-structure must be unsorted!', &
    !        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
    !    call sys_halt()
    !  end if

    !  if ((.not. associated(rmatrix%p_rspatialDiscrTest)) .or. &
    !      (.not. associated(rmatrix%p_rspatialDiscrTrial))) then
    !    call output_line ('No discretisation associated!', &
    !        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
    !    call sys_halt()
    !  end if

    ! If we do not have it, create a scalar assembly info structure that
    ! defines how to do the assembly.
    if (.not. present(rscalarAssemblyInfo)) then
       call easminfo_createDefInfoStructure(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial,&
            rlocalScalarAssemblyInfo,0)
       p_rscalarAssemblyInfo => rlocalScalarAssemblyInfo
    else
       p_rscalarAssemblyInfo => rscalarAssemblyInfo
    end if

    ! Do we have a uniform triangulation? Would simplify a lot...
    select case (rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%ccomplexity)
    case (SPDISC_UNIFORM,SPDISC_CONFORMAL) 
       ! Uniform and conformal discretisations
       select case (rmatrix%RmatrixBlock(1,1)%cdataType)
       case (ST_DOUBLE) 
          ! Which matrix structure do we have?
          select case (rmatrix%RmatrixBlock(1,1)%cmatrixFormat) 
          case (LSYSSC_MATRIX9,LSYSSC_MATRIX9ROWC)

             !        ! Probably allocate/clear the matrix
             !        if (rmatrix%h_DA .eq. ST_NOHANDLE) then
             !          call lsyssc_allocEmptyMatrix(rmatrix,LSYSSC_SETM_ZERO)
             !        else
             if (bclear) call lsysbl_clearMatrix (rmatrix)
             !        end if

             ! Loop over the element blocks to discretise
             do iinfoBlock = 1,p_rscalarAssemblyInfo%ninfoBlockCount

                ! Get the element distribution of that block.
                ielementDistr = p_rscalarAssemblyInfo%p_RinfoBlocks(iinfoBlock)%ielementDistr

                ! Check if element distribution is empty
                if (p_rscalarAssemblyInfo%p_RinfoBlocks(iinfoBlock)%NEL .le. 0 ) cycle

                ! Get list of elements present in the element distribution.
                ! If the handle of the info block structure is not associated,
                ! take all elements of the corresponding element distribution.
                if (p_rscalarAssemblyInfo%p_RinfoBlocks(iinfoBlock)%h_IelementList .ne. ST_NOHANDLE) then
                   call storage_getbase_int(&
                        p_rscalarAssemblyInfo%p_RinfoBlocks(iinfoBlock)%h_IelementList,&
                        p_IelementList)
                else
                   call storage_getbase_int(&
                        rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial%RelementDistr(ielementDistr)%h_IelementList,&
                        p_IelementList)
                end if

                ! Initialise a matrix assembly structure for that element distribution
                call bilf_initAssembly(rmatrixAssembly,rform,&
                     rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(ielementDistr)%celement,&
                     rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial%RelementDistr(ielementDistr)%celement,&
                     p_rscalarAssemblyInfo%p_RinfoBlocks(iinfoBlock)%ccubature,&
                     min(BILF_NELEMSIM,p_rscalarAssemblyInfo%p_RinfoBlocks(iinfoBlock)%NEL))

                ! Assemble the data for all elements in this element distribution
                call bilf_assembleSubmeshMatrix9Block (rmatrixAssembly,rmatrix,&
                     p_IelementList,fcoeff_buildMatrixBl_sim,rcollection)

                ! Release the assembly structure.
                call bilf_doneAssembly(rmatrixAssembly)
             end do

          case default
             call output_line ('Not supported matrix structure!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
             call sys_halt()
          end select

       case default
          call output_line ('Single precision matrices currently not supported!', &
               OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
          call sys_halt()
       end select

    case default
       call output_line ('General discretisation not implemented!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
       call sys_halt()
    end select

    ! Release the assembly structure if necessary.
    if (.not. present(rscalarAssemblyInfo)) then
       call easminfo_releaseInfoStructure(rlocalScalarAssemblyInfo)
    end if

  end subroutine bilf_buildMatrixBlock2


  !****************************************************************************

  !<subroutine>  

  subroutine bilf_assembleSubmeshMatrix9Block (rmatrixAssembly, rmatrix, IelementList,&
       fcoeff_buildMatrixBl_sim, rcollection)

    !<description>

    ! Assembles the matrix entries for a submesh by integrating over the domain.

    !</description>

    !<input>

    ! List of elements where to assemble the bilinear form.
    integer, dimension(:), intent(in), target :: IelementList

    ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
    ! Must be present if the matrix has nonconstant coefficients!
    include 'intf_coefficientMatrixBl.inc'
    optional :: fcoeff_buildMatrixBl_sim

    !</input>

    !<inputoutput>

    ! A matrix assembly structure prepared with bilf_initAssembly.
    type(t_bilfMatrixAssembly), intent(inout), target :: rmatrixAssembly

    ! A matrix where to assemble the contributions to.
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: A pointer to a collection structure. This structure is given to the
    ! callback function for nonconstant coefficients to provide additional
    ! information. 
    type(t_collection), intent(inout), target, optional :: rcollection

    !</inputoutput>

    !</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_DA
    integer :: indofTest,indofTrial,ncubp

    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax
    integer :: iel,icubp,ialbet,ia,ib,idofe,jdofe
    real(DP) :: domega,daux,db
    integer(I32) :: cevaluationTag
    type(t_bilfMatrixAssembly), target :: rlocalMatrixAssembly
    type(t_domainIntSubset) :: rintSubset
    integer, dimension(:,:,:), pointer :: p_Kentry
    real(DP), dimension(:,:,:,:,:), pointer :: p_Dentry
    real(DP), dimension(:,:), pointer :: p_Ddetj
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial
    real(DP), dimension(:,:,:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:), pointer :: p_DcoefficientsBilf
    integer, dimension(:,:), pointer :: p_IdofsTest
    integer, dimension(:,:), pointer :: p_IdofsTrial
    type(t_evalElementSet), pointer :: p_revalElementSet
    integer, dimension(:,:),pointer :: p_Idescriptors
    integer :: nblocksi, nblocksj
    integer :: iblock, jblock

    !    type t_array
    !     ! Pointer to the double-valued matrix or vector data
    !     real(DP), dimension(:), pointer :: Da
    !    end type t_array
    !    
    !    type(t_array), dimension(:,:), allocatable:: p_matrixBlockDataPointers


    ! Get some pointers for faster access
    !    call lsyssc_getbase_double (rmatrix,p_DA)
    indofTest  = rmatrixAssembly%indofTest
    indofTrial = rmatrixAssembly%indofTrial
    ncubp      = rmatrixAssembly%ncubp

    ! Get number of blocks in the matrix and test if this equals the number of rows
    nblocksi = rmatrix%nblocksPerRow
    nblocksj = rmatrix%nblocksPerCol

    ! OpenMP-Extension: Copy the matrix assembly data to the local
    ! matrix assembly data, where we can allocate memory.
    !
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the matrix
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original matrix assembly
    ! stucture or disturbing the data of the other processors.
    !
    !$omp parallel default(shared) &
    !$omp private(IELmax,cevaluationTag,daux,db,domega,ia,ialbet,ib,icubp,&
    !$omp         idofe,iel,jdofe,p_DbasTest,p_DbasTrial,p_Dcoefficients,&
    !$omp         p_DcoefficientsBilf,p_Ddetj,p_Dentry,p_Domega,p_Idescriptors,&
    !$omp         p_IdofsTest,p_IdofsTrial,p_Kentry,p_revalElementSet,&
    !$omp         rintSubset,rlocalMatrixAssembly)
    rlocalMatrixAssembly = rmatrixAssembly
    call bilf_allocAssemblyData(rlocalMatrixAssembly)

    ! Allocate space for the local matrices, coefficients and pointers
    allocate(p_Dentry(indofTrial,indofTest,nblocksi,nblocksj,rmatrixAssembly%nelementsPerBlock))
    allocate(p_Dcoefficients(nblocksi,nblocksj,rlocalMatrixAssembly%rform%itermcount,ncubp,rmatrixAssembly%nelementsPerBlock))
    !    allocate(p_matrixBlockDataPointers(nblocksi,nblocksj))

    ! Get some more pointers to local data.
    p_Kentry            => rlocalMatrixAssembly%p_Kentry
    !    p_Dentry            => rlocalMatrixAssembly%p_Dentry
    p_Domega            => rlocalMatrixAssembly%p_Domega
    p_DbasTest          => rlocalMatrixAssembly%p_DbasTest
    p_DbasTrial         => rlocalMatrixAssembly%p_DbasTrial
    !    p_Dcoefficients     => rlocalMatrixAssembly%p_Dcoefficients
    p_Idescriptors      => rlocalMatrixAssembly%rform%Idescriptors
    p_IdofsTest         => rlocalMatrixAssembly%p_IdofsTest
    p_IdofsTrial        => rlocalMatrixAssembly%p_IdofsTrial
    p_revalElementSet   => rlocalMatrixAssembly%revalElementSet
    p_DcoefficientsBilf => rlocalMatrixAssembly%rform%Dcoefficients

    ! Loop over the elements - blockwise.
    !
    ! OpenMP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !$omp do schedule(static,1)
    do IELset = 1, size(IelementList), rmatrixAssembly%nelementsPerBlock

       ! We always handle nelementsPerBlock elements simultaneously.
       ! How many elements have we actually here?
       ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
       ! elements simultaneously.

       IELmax = min(size(IelementList),IELset-1+rmatrixAssembly%nelementsPerBlock)

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
       ! global DOF`s of our BILF_NELEMSIM elements simultaneously.
       call dof_locGlobMapping_mult(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest, &
            IelementList(IELset:IELmax), p_IdofsTest)

       ! If the DOF`s for the test functions are different, calculate them, too.
       if (.not. rlocalMatrixAssembly%bIdenticalTrialAndTest) then
          call dof_locGlobMapping_mult(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial, &
               IelementList(IELset:IELmax), p_IdofsTrial)
       end if

       ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------

       ! For the assembly of the global matrix, we use a "local"
       ! approach. At first we build a "local" system matrix according
       ! to the current element. This contains all additive
       ! contributions of element iel, which are later added at the
       ! right positions to the elements in the global system matrix.
       !
       ! We have indofTrial trial DOF`s per element and
       ! indofTest test DOF`s per element. Therefore there are
       ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
       ! "active" (i.e. have common support) on our current element, each 
       ! giving an additive contribution to the system matrix.
       !
       ! We build a quadratic indofTrial*indofTest local matrix:
       ! Kentry(1..indofTrial,1..indofTest) receives the position 
       ! in the global system matrix, where the corresponding value 
       ! has to be added to.
       ! (The corresponding contributions can be saved separately, 
       ! but we directly add them to the global matrix in this 
       ! approach.)
       !
       ! We build local matrices for all our elements 
       ! in the set simultaneously. Get the positions of the local matrices
       ! in the global matrix.
       call bilf_getLocalMatrixIndices (rmatrix%RmatrixBlock(1,1),p_IdofsTest,p_IdofsTrial,p_Kentry,&
            ubound(p_IdofsTest,1),ubound(p_IdofsTrial,1),IELmax-IELset+1)

       ! -------------------- ELEMENT EVALUATION PHASE ----------------------

       ! Ok, we found the positions of the local matrix entries
       ! that we have to change.
       ! To calculate the matrix contributions, we have to evaluate
       ! the elements to give us the values of the basis functions
       ! in all the DOF`s in all the elements in our set.

       ! Get the element evaluation tag of all FE spaces. We need it to evaluate
       ! the elements later. All of them can be combined with OR, what will give
       ! a combined evaluation tag. 
       cevaluationTag = rlocalMatrixAssembly%cevaluationTag

       ! In the first loop, calculate the coordinates on the reference element.
       ! In all later loops, use the precalculated information.
       !
       ! If the cubature points are already initialised, do not do it again.
       ! We check this by taking a look to iinitialisedElements which
       ! gives the current maximum of initialised elements.
       if (IELmax .gt. rlocalMatrixAssembly%iinitialisedElements) then

          ! (Re-)initialise!
          cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)

          ! Remember the new number of initialised elements
          rlocalMatrixAssembly%iinitialisedElements = IELmax

       else
          ! No need.
          cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
       end if

       ! Calculate all information that is necessary to evaluate the finite element
       ! on all cells of our subset. This includes the coordinates of the points
       ! on the cells.
       call elprep_prepareSetForEvaluation (p_revalElementSet,&
            cevaluationTag, rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation, &
            IelementList(IELset:IELmax), rlocalMatrixAssembly%ctrafoType, &
            rlocalMatrixAssembly%p_DcubPtsRef(:,1:ncubp))
       p_Ddetj => p_revalElementSet%p_Ddetj

       ! If the matrix has nonconstant coefficients, calculate the coefficients now.

       call domint_initIntegrationByEvalSet (p_revalElementSet,rintSubset)
       rintSubset%ielementDistribution =  1
       rintSubset%ielementStartIdx     =  IELset
       rintSubset%p_Ielements          => IelementList(IELset:IELmax)
       rintSubset%p_IdofsTrial         => p_IdofsTrial
       rintSubset%celement             =  rlocalMatrixAssembly%celementTrial
       call fcoeff_buildMatrixBl_sim (rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest,&
            rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial,&
            rlocalMatrixAssembly%rform, IELmax-IELset+1, ncubp,&
            p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
            p_IdofsTrial, p_IdofsTest, rintSubset, &
            p_Dcoefficients(:,:,:,:,1:IELmax-IELset+1), rcollection)
       call domint_doneIntegration (rintSubset)

       ! Calculate the values of the basis functions.
       call elem_generic_sim2 (rlocalMatrixAssembly%celementTest, &
            p_revalElementSet, rlocalMatrixAssembly%BderTest, &
            rlocalMatrixAssembly%p_DbasTest)

       ! Omit the calculation of the trial function values if they
       ! are identical to the test function values.
       if (.not. rlocalMatrixAssembly%bidenticalTrialAndTest) then
          call elem_generic_sim2 (rlocalMatrixAssembly%celementTrial, &
               p_revalElementSet, rlocalMatrixAssembly%BderTrial, &
               rlocalMatrixAssembly%p_DbasTrial)
       end if

       ! --------------------- DOF COMBINATION PHASE ------------------------

       ! Values of all basis functions calculated. Now we can start 
       ! to integrate!

       ! Clear the local matrices
       p_Dentry(:,:,:,:,1:IELmax-IELset+1) = 0.0_DP



       ! Nonconstant coefficients. The coefficients are to be found in
       ! the Dcoefficients variable as computed above.
       !
       ! Loop over the elements in the current set.

       do iel = 1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

             ! calculate the current weighting factor in the cubature formula
             ! in that cubature point.
             !
             ! Take the absolut value of the determinant of the mapping.
             ! In 2D, the determinant is always positive, whereas in 3D,
             ! the determinant might be negative -- that is normal!

             domega = p_Domega(icubp)*abs(p_Ddetj(icubp,iel))

             ! Loop over the additive factors in the bilinear form.
             do ialbet = 1,rlocalMatrixAssembly%rform%itermcount

                ! Get from Idescriptors the type of the derivatives for the 
                ! test and trial functions. The summand we calculate
                ! here will be added to the matrix entry:
                !
                ! a_ij  =  int_... ( psi_j )_ia  *  ( phi_i )_ib
                !
                ! -> Ix=0: function value, 
                !      =1: first derivative, ...
                !    as defined in the module 'derivative'.

                ia = rlocalMatrixAssembly%rform%Idescriptors(1,ialbet)
                ib = rlocalMatrixAssembly%rform%Idescriptors(2,ialbet)

                ! Multiply domega with the coefficient of the form.
                ! This gives the actual value to multiply the
                ! function value with before summing up to the integral.
                ! Get the precalculated coefficient from the coefficient array.
                !daux = domega * p_Dcoefficients(ialbet,icubp,iel)

                ! Now loop through all possible combinations of DOF`s
                ! in the current cubature point. The outer loop
                ! loops through the "O" in the above picture,
                ! the test functions:

                do iblock = 1, nblocksi
                   do jblock = 1, nblocksj

                      do idofe = 1,indofTest

                         ! Get the value of the (test) basis function 
                         ! phi_i (our "O") in the cubature point:
                         db = p_DbasTest(idofe,ib,icubp,iel)

                         ! Perform an inner loop through the other DOF`s
                         ! (the "X"). 

                         do jdofe = 1,indofTrial

                            ! Get the value of the basis function 
                            ! psi_j (our "X") in the cubature point. 
                            ! Them multiply:
                            !    db * dbas(..) * daux
                            ! ~= phi_i * psi_j * coefficient * cub.weight
                            ! Summing this up gives the integral, so the contribution
                            ! to the global matrix. 
                            !
                            ! Simply summing up db * dbas(..) * daux would give
                            ! the coefficient of the local matrix. We save this
                            ! contribution in the local matrix of element iel.

                            !JCOLB = Kentry(jdofe,idofe,iel)
                            !p_DA(JCOLB) = p_DA(JCOLB) + db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
                            p_Dentry(jdofe,idofe,iblock,jblock,iel) = &
                                 p_Dentry(jdofe,idofe,iblock,jblock,iel)+db*p_DbasTrial(jdofe,ia,icubp,iel)*domega*p_Dcoefficients(iblock,jblock,ialbet,icubp,iel)

                         end do

                      end do ! jdofe
                   end do
                end do

             end do ! ialbet

          end do ! icubp 

       end do ! iel

       ! Incorporate the local matrices into the global one.
       ! Kentry gives the position of the additive contributions in Dentry.
       !
       ! OpenMP-Extension: This is a critical section. Only one thread is
       ! allowed to write to the matrix, otherwise the matrix may get
       ! messed up.
       ! The critical section is put around both loops as indofTest/indofTrial
       ! are usually small and quickly to handle.
       !



       !      ! Get pointers to the data entries of the block matrix
       !      do iblock = 1, nblocksi
       !        do jblock = 1, nblocksj
       !          call lsyssc_getbase_double(rmatrix%RmatrixBlock(iblock,jblock),p_matrixBlockDataPointers(iblock,jblock)%Da)
       !        end do
       !      end do

       !$omp critical
       do iel = 1,IELmax-IELset+1
          do iblock = 1, nblocksi
             do jblock = 1, nblocksj
                call lsyssc_getbase_double(rmatrix%RmatrixBlock(iblock,jblock),p_Da)
                do idofe = 1,indofTest
                   do jdofe = 1,indofTrial
                      !            p_DA(p_Kentry(jdofe,idofe,iel)) = &
                      !                p_DA(p_Kentry(jdofe,idofe,iel)) + p_Dentry(jdofe,idofe,iel)
                      !            p_matrixBlockDataPointers(iblock,jblock)%Da(p_Kentry(jdofe,idofe,iel)) = &
                      !                  p_matrixBlockDataPointers(iblock,jblock)%Da(p_Kentry(jdofe,idofe,iel)) + p_Dentry(jdofe,idofe,iblock,jblock,iel)
                      p_Da(p_Kentry(jdofe,idofe,iel)) = &
                           p_Da(p_Kentry(jdofe,idofe,iel)) + p_Dentry(jdofe,idofe,iblock,jblock,iel)
                   end do
                end do
             end do
          end do

       end do ! iel
       !$omp end critical

    end do ! IELset
    !$omp end do

    ! Release the local matrix assembly structure
    call bilf_releaseAssemblyData(rlocalMatrixAssembly)
    !$omp end parallel

    ! Deallocate local matrices, coefficients and pointers
    deallocate(p_Dentry,p_Dcoefficients)!,p_matrixBlockDataPointers)

  end subroutine bilf_assembleSubmeshMatrix9Block







































!****************************************************************************

  !<subroutine>

  subroutine bilf_dg_buildMatrixBlEdge2D_ss (rform, ccubType, bclear, rmatrix,&
       rvectorSol, raddTriaData,&
       flux_dg_buildMatrixBlEdge2D_sim_ss,&
       rcollection, cconstrType)

    !<description>
    ! This routine calculates the entries of a finite element matrix in 2D.
    ! The matrix structure must be prepared with bilf_createMatrixStructure
    ! in advance.
    ! In case the array for the matrix entries does not exist, the routine
    ! allocates memory in size of the matrix of the heap for the matrix entries.
    !
    ! For setting up the entries, the discretisation structure attached to
    ! the matrix is used (rmatrix%p_rdiscretisation). This is
    ! normally attached to the matrix by bilf_createMatrixStructure.
    !
    ! The matrix must be unsorted when this routine is called, 
    ! otherwise an error is thrown.
    !</description>

    !<input>
    ! The bilinear form specifying the underlying PDE of the discretisation.
    type(t_bilinearForm), intent(in) :: rform

    ! A line cubature formula CUB_xxxx_1D to be used for line integration.
    integer(I32), intent(in) :: ccubType

    ! The solution vector. Used to calculate the solution on the edges.
    type(t_vectorBlock), intent(in) :: rvectorSol

    ! Additional triangulation data
    type(t_additionalTriaData), intent(in) :: raddTriaData

    ! Whether to clear the matrix before calculating the entries.
    ! If .FALSE., the new matrix entries are added to the existing entries.
    logical, intent(in) :: bclear

    ! A callback routine for the flux function.
    include 'intf_flux_dg_buildMatrixBlEdge2D_ss.inc'
    optional :: flux_dg_buildMatrixBlEdge2D_sim_ss

    ! OPTIONAL: One of the BILF_MATC_xxxx constants that allow to specify
    ! the matrix construction method. If not specified,
    ! BILF_MATC_EDGEBASED is used.
    integer, intent(in), optional :: cconstrType
    !</input>

    !<inputoutput>
    ! The FE matrix. Calculated matrix entries are imposed to this matrix.
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: A collection structure. This structure is given to the
    ! callback function for nonconstant coefficients to provide additional
    ! information. 
    type(t_collection), intent(inout), target, optional :: rcollection
    !</inputoutput>

    !</subroutine>

    ! local variables
    type(t_bilfMatrixAssembly), dimension(2) :: rmatrixAssembly
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: IelementList, p_IedgeList
    integer :: ccType
    integer :: iedge, ielementDistr

    ! The matrix must be unsorted, otherwise we can not set up the matrix.
    ! Note that we cannot switch off the sorting as easy as in the case
    ! of a vector, since there is a structure behind the matrix! So the caller
    ! has to make sure, the matrix is unsorted when this routine is called.
    if (rmatrix%RmatrixBlock(1,1)%isortStrategy .gt. 0) then
       call output_line ('Matrix-structure must be unsorted!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ! The matrix must provide discretisation structures
    if ((.not. associated(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest)) .or. &
         (.not. associated(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial))) then
       call output_line ('No discretisation associated!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ! The discretisation must provide a triangulation structure
    if ((.not. associated(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation)) .or. &
         (.not. associated(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial%p_rtriangulation))) then
       call output_line('No triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    !  ! The discretisation must provide a boundary structure
    !  if ((.not. associated(rmatrix%p_rspatialDiscrTest%p_rboundary)) .or. &
    !      (.not. associated(rmatrix%p_rspatialDiscrTrial%p_rboundary))) then
    !    call output_line('No boundary associated!',&
    !        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
    !    call sys_halt()
    !  end if

    !  ! Set pointers for quicker access
    !  p_rboundary => rmatrix%p_rspatialDiscrTest%p_rboundary
    !  if (.not.associated(p_rboundary, rmatrix%p_rspatialDiscrTrial%p_rboundary)) then
    !    call output_line('Invalid boundary associated!',&
    !        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
    !    call sys_halt()
    !  end if

    p_rtriangulation => rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation
    if (.not.associated(p_rtriangulation, rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial%p_rtriangulation)) then
       call output_line('Invalid triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ccType = BILF_MATC_EDGEBASED
    if (present(cconstrType)) ccType = cconstrType

    ! Do we have a uniform triangulation? Would simplify a lot...
    select case (rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%ccomplexity)
    case (SPDISC_UNIFORM,SPDISC_CONFORMAL) 
       ! Uniform and conformal discretisations
       select case (rmatrix%RmatrixBlock(1,1)%cdataType)
       case (ST_DOUBLE) 
          ! Which matrix structure do we have?
          select case (rmatrix%RmatrixBlock(1,1)%cmatrixFormat) 
          case (LSYSSC_MATRIX9)

             ! Probably allocate/clear the matrix
             !        if (rmatrix%h_DA .eq. ST_NOHANDLE) then
             !          call lsyssc_allocEmptyMatrix(rmatrix,LSYSSC_SETM_ZERO)
             !        else
             if (bclear) call lsysbl_clearMatrix (rmatrix)
             !        end if


             ! Allocate the edgelist
             allocate(p_IedgeList(rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%NMT))

             ! All edges
             forall (iedge = 1:rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%NMT) p_IedgeList(iedge)=iedge

             ! Initialise a matrix assembly structure for that element distribution
             ielementDistr = 1
             call bilf_initAssembly(rmatrixAssembly(1),rform,&
                  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(ielementDistr)%celement,&
                  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial%RelementDistr(ielementDistr)%celement,&
                  ccubType, BILF_NELEMSIM)

             ! Do the same for the other side of the egde
             ielementDistr = 1
             call dg_bilf_initAssembly_reverseCubPoints(rmatrixAssembly(2),rform,&
                  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(ielementDistr)%celement,&
                  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTrial%RelementDistr(ielementDistr)%celement,&
                  ccubType, BILF_NELEMSIM)

             ! Assemble the data for all elements in this element distribution
             call dg_bilf_assembleSubmeshMat9Bdr2D_Block_ss (rmatrixAssembly, rmatrix,&
                  rvectorSol, raddTriaData, p_IedgeList(1:rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%NMT),&
                  ccType, flux_dg_buildMatrixBlEdge2D_sim_ss, rcollection)


             ! Release the assembly structure.
             call bilf_doneAssembly(rmatrixAssembly(1))
             call bilf_doneAssembly(rmatrixAssembly(2))

             ! Deallocate the edgelist
             deallocate(p_IedgeList)



          case default
             call output_line ('Not supported matrix structure!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
             call sys_halt()
          end select

       case default
          call output_line ('Single precision matrices currently not supported!', &
               OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
          call sys_halt()
       end select

    case default
       call output_line ('General discretisation not implemented!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end select

  end subroutine bilf_dg_buildMatrixBlEdge2D_ss



  !****************************************************************************

  !<subroutine>  

  subroutine dg_bilf_assembleSubmeshMat9Bdr2D_Block_ss (rmatrixAssembly, rmatrix,&
       rvectorSol, raddTriaData, IedgeList,&
       cconstrType, flux_dg_buildMatrixBlEdge2D_sim_ss, rcollection)

    !<description>

    ! Assembles the matrix entries for a submesh by integrating over the
    ! boundary region in 2D.

    !</description>

    !<input>

    ! List of elements where to assemble the bilinear form.
    integer, dimension(:), intent(in), target :: IedgeList

    ! One of the BILF_MATC_xxxx constants that allow to specify the
    ! matrix construction method.
    integer, intent(in) :: cconstrType

    ! The solution vector. Used to calculate the solution on the edges.
    type(t_vectorBlock), intent(in) :: rvectorSol

    ! Additional triangulation data
    type(t_additionalTriaData), intent(in) :: raddTriaData

    ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
    ! Must be present if the matrix has nonconstant coefficients!
    include 'intf_flux_dg_buildMatrixBlEdge2D_ss.inc'
    optional :: flux_dg_buildMatrixBlEdge2D_sim_ss

    !</input>

    !<inputoutput>

    ! A matrix assembly structure prepared with bilf_initAssembly.
    type(t_bilfMatrixAssembly), intent(inout), dimension(2), target :: rmatrixAssembly

    ! A matrix where to assemble the contributions to.
    type(t_matrixBlock), intent(inout) :: rmatrix

    ! OPTIONAL: A pointer to a collection structure. This structure is given to the
    ! callback function for nonconstant coefficients to provide additional
    ! information. 
    type(t_collection), intent(inout), target, optional :: rcollection

    !</inputoutput>

    !</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_DA
    integer :: indofTest,indofTrial,ncubp

    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax,ibdc,k
    integer :: iel,icubp,ialbet,ia,ib,idofe,jdofe,nve
    real(DP) :: domega1,domega2,db1,db2,dlen
    real(dp), dimension(:,:), allocatable :: daux1, daux2
    integer(I32) :: cevaluationTag
    type(t_bilfMatrixAssembly), dimension(2), target :: rlocalMatrixAssembly
    type(t_domainIntSubset), dimension(2) :: rintSubset
    integer, dimension(:,:,:), pointer :: p_Kentryii, p_Kentryia, p_Kentryai, p_Kentryaa
    real(DP), dimension(:,:,:,:,:), pointer :: p_Dentryii, p_Dentryia, p_Dentryai, p_Dentryaa
    real(DP), dimension(:,:,:), pointer :: p_Dcoords
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:), pointer :: p_Dside
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial
    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    real(DP), dimension(:), pointer :: p_DcoefficientsBilf
    integer, dimension(:,:), pointer :: p_IdofsTest
    integer, dimension(:,:), pointer :: p_IdofsTrial
    type(t_evalElementSet), pointer :: p_revalElementSet
    integer, dimension(:,:),pointer :: p_Idescriptors
    integer, dimension(:,:), allocatable, target :: IelementList
    integer :: nvar

    ! Pointer to Ielementsatedge in the triangulation
    integer, dimension(:,:), pointer :: p_IelementsAtEdge

    ! Pointer to IverticesAtEdge in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtEdge

    ! Pointer to the vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Pointer to IverticesAtEelement in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Space for the values of the flux function
    real(DP), dimension(:,:,:,:,:), allocatable :: DfluxValues

    ! Arrays for cubature points 1D->2D
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D_1, Dxi1D_2
    real(DP), dimension(:,:,:,:), allocatable :: Dxi2D,DpointsRef
    real(DP), dimension(:,:), allocatable :: DpointsPar
    real(DP), dimension(:), allocatable :: DedgeLength

    integer(i32) :: icoordSystem
    integer :: NEL
    integer :: iside
    logical :: bisLinearTrafo

    integer :: iblock,jblock

    !    type t_array
    !     ! Pointer to the double-valued matrix or vector data
    !     real(DP), dimension(:), pointer :: Da
    !    end type t_array

    !    type(t_array), dimension(:,:), allocatable :: p_matrixBlockDataPointers

    !    ! Boundary component?
    !    ibdc = rboundaryRegion%iboundCompIdx

    ! Get some pointers for faster access
    ! call lsyssc_getbase_double (rmatrix,p_DA)
    indofTest = rmatrixAssembly(1)%indofTest
    indofTrial = rmatrixAssembly(1)%indofTrial
    ncubp = rmatrixAssembly(1)%ncubp
    nvar = rvectorSol%nblocks

    ! Open-MP-Extension: Copy the matrix assembly data to the local
    ! matrix assembly data, where we can allocate memory.
    !
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the matrix
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original matrix assembly
    ! stucture or disturbing the data of the other processors.
    !
    !$omp parallel default(shared) &
    !$omp private(DedgeLength,DpointsPar,DpointsRef,Dxi1D,Dxi2D,IELmax,bisLinearTrafo,&
    !$omp         cevaluationTag,daux,db,dlen,domega,ia,ialbet,ib,icoordSystem,icubp,&
    !$omp         idofe,iel,jdofe,k,p_DbasTest,p_DbasTrial,p_Dcoefficients,&
    !$omp         p_DcoefficientsBilf,p_Dcoords,p_DcubPtsRef,p_Dentry,p_Domega,&
    !$omp         p_Idescriptors,p_IdofsTest,p_IdofsTrial,p_Kentry,&
    !$omp         p_revalElementSet,rintSubset,rlocalMatrixAssembly)
    rlocalMatrixAssembly(1) = rmatrixAssembly(1)
    rlocalMatrixAssembly(2) = rmatrixAssembly(2)
    call bilf_allocAssemblyData(rlocalMatrixAssembly(1))
    call bilf_allocAssemblyData(rlocalMatrixAssembly(2))

    ! Allocate space for the positions of the DOFs in the matrix
    allocate(p_Kentryii(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryia(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryai(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryaa(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))

    ! Allocate auxiliary vectors
    allocate(daux1(nvar,nvar),daux2(nvar,nvar))

    ! Allocate space for the coefficient of the solutions DOFs on each side of the edge
    allocate(p_Dside(2,rmatrixAssembly(1)%rform%itermCount))

    ! Allocate space for the entries in the local matrices
    allocate(p_Dentryii(nvar,nvar,rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryia(nvar,nvar,rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryai(nvar,nvar,rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryaa(nvar,nvar,rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))



    ! Allocate space for the flux variables DIM(nvar,nvar,ialbet,ncubp,elementsperblock)
    allocate(DfluxValues(nvar,nvar,rmatrixAssembly(1)%rform%itermCount,ncubp,rlocalMatrixAssembly(1)%nelementsPerBlock))

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


    ! Get number of elements
    NEL = rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%NEL

    ! Get pointers to elements at edge
    call storage_getbase_int2D(&
         rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%h_IelementsAtEdge,&
         p_IelementsAtEdge)

    ! Get pointers to the vertex coordinates
    call storage_getbase_double2D(&
         rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%h_DvertexCoords,&
         p_DvertexCoords)

    ! Get pointers to vertices at edge
    call storage_getbase_int2D(&
         rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtEdge,&
         p_IverticesAtEdge)

    ! Get pointers to vertices at elements
    call storage_getbase_int2D(&
         rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtElement,&
         p_IverticesAtElement)   

    ! Get the elements adjacent to the given edges
    allocate(IelementList(3,size(IedgeList)))
    IelementList(1:2,1:size(IedgeList))=p_IelementsAtEdge(1:2,IedgeList(:))

    ! Copy the second component and replace 0s by 1s
    IelementList(3,size(IedgeList))=IelementList(2,size(IedgeList))
    do iel = 1,size(IedgeList)
       IelementList(3,iel)=max(IelementList(2,iel),1)
    end do

    ! Transpose the coordinate array such that we get coordinates we
    ! can work with in the mapping between 1D and 2D.
    do k = 1, ubound(rlocalmatrixAssembly(1)%p_DcubPtsRef,1)
       do icubp = 1,ncubp
          Dxi1D_1(icubp,k) = rlocalmatrixAssembly(1)%p_DcubPtsRef(k,icubp)
          Dxi1D_2(icubp,k) = rlocalmatrixAssembly(2)%p_DcubPtsRef(k,icubp)
       end do
    end do

    ! Allocate memory for the cubature points in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,2,rlocalMatrixAssembly(1)%nelementsPerBlock))

    ! Allocate memory for the coordinates of the reference points
    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalMatrixAssembly(1)%nelementsPerBlock,2))

    !    ! Allocate memory for the parameter values of the points on the boundary
    !    allocate(DpointsPar(ncubp,rlocalMatrixAssembly%nelementsPerBlock))
    !
    !    ! Allocate memory for the length of edges on the boundary
    !    allocate(DedgeLength(rlocalMatrixAssembly%nelementsPerBlock))

    ! Get the type of coordinate system
    icoordSystem = elem_igetCoordSystem(rlocalMatrixAssembly(1)%celementTrial)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !$omp do schedule(static,1)
    do IELset = 1, size(IedgeList), rlocalMatrixAssembly(1)%nelementsPerBlock

       ! We always handle nelementsPerBlock elements simultaneously.
       ! How many elements have we actually here?
       ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
       ! elements simultaneously.

       IELmax = min(size(IedgeList),IELset-1+rlocalMatrixAssembly(1)%nelementsPerBlock)

       ! Map the 1D cubature points to the edges in 2D.
       do iel = 1,IELmax-IELset+1
          call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(1,Iedgelist(IELset+iel-1)), &
               ncubp, Dxi1D_1, Dxi2D(:,:,1,iel))
          call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(2,Iedgelist(IELset+iel-1)), &
               ncubp, Dxi1D_2, Dxi2D(:,:,2,iel))
       end do

       !      ! Calculate the parameter values of the points
       !      do iel = 1,IELmax-IELset+1
       !        do icubp = 1,ncubp
       !          ! Dxi1D is in [-1,1] while the current edge has parmeter values
       !          ! [DedgePosition(1),DedgePosition(2)]. So do a linear
       !          ! transformation to transform Dxi1D into that interval, this 
       !          ! gives the parameter values in length parametrisation
       !          call mprim_linearRescale(Dxi1D(icubp,1), -1.0_DP, 1.0_DP,&
       !              DedgePosition(1,IELset+iel-1), DedgePosition(2,IELset+iel-1),&
       !              DpointsPar(icubp,iel))
       !        end do
       !      end do

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
       ! global DOF`s of our BILF_NELEMSIM elements simultaneously.
       !      call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTest, &
       !          IelementList(IELset:IELmax), p_IdofsTest)
       call dof_locGlobMapping_mult( rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest, &
            IelementList(1,IELset:IELmax), rlocalMatrixAssembly(1)%p_IdofsTest)
       call dof_locGlobMapping_mult( rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest, &
            IelementList(3,IELset:IELmax), rlocalMatrixAssembly(2)%p_IdofsTest)


       !      ! If the DOF`s for the trial functions are different, calculate them, too.
       !      if (.not. rlocalMatrixAssembly%bIdenticalTrialAndTest) then
       !        call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTrial, &
       !            IelementList(IELset:IELmax), p_IdofsTrial)
       !      end if

       ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------

       ! For the assembly of the global matrix, we use a "local"
       ! approach. At first we build a "local" system matrix according
       ! to the current element. This contains all additive
       ! contributions of element iel, which are later added at the
       ! right positions to the elements in the global system matrix.
       !
       ! We have indofTrial trial DOF`s per element and
       ! indofTest test DOF`s per element. Therefore there are
       ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
       ! "active" (i.e. have common support) on our current element, each 
       ! giving an additive contribution to the system matrix.
       !
       ! We build a quadratic indofTrial*indofTest local matrix:
       ! Kentry(1..indofTrial,1..indofTest) receives the position 
       ! in the global system matrix, where the corresponding value 
       ! has to be added to.
       ! (The corresponding contributions can be saved separately, 
       ! but we directly add them to the global matrix in this 
       ! approach.)
       !
       ! We build local matrices for all our elements 
       ! in the set simultaneously. Get the positions of the local matrices
       ! in the global matrix.
       !      call bilf_getLocalMatrixIndices (rmatrix,p_IdofsTest,p_IdofsTrial,p_Kentry,&
       !          ubound(p_IdofsTest,1),ubound(p_IdofsTrial,1),IELmax-IELset+1)  
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
       !            rlocalMatrixAssembly(1)%p_IdofsTrial, p_Kentryii,&
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
       !            rlocalMatrixAssembly(2)%p_IdofsTrial, p_Kentryia,&
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
       !            rlocalMatrixAssembly(1)%p_IdofsTrial, p_Kentryai,&
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
       !            rlocalMatrixAssembly(2)%p_IdofsTrial, p_Kentryaa,&
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTrial,1), IELmax-IELset+1)



       call bilf_getLocalMatrixIndices (rmatrix%RmatrixBlock(1,1),rlocalMatrixAssembly(1)%p_IdofsTest, &
            rlocalMatrixAssembly(1)%p_IdofsTest, p_Kentryii,&
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix%RmatrixBlock(1,1),rlocalMatrixAssembly(1)%p_IdofsTest, &
            rlocalMatrixAssembly(2)%p_IdofsTest, p_Kentryai,&
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix%RmatrixBlock(1,1),rlocalMatrixAssembly(2)%p_IdofsTest, &
            rlocalMatrixAssembly(1)%p_IdofsTest, p_Kentryia,&
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix%RmatrixBlock(1,1),rlocalMatrixAssembly(2)%p_IdofsTest, &
            rlocalMatrixAssembly(2)%p_IdofsTest, p_Kentryaa,&
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), IELmax-IELset+1)   

       ! -------------------- ELEMENT EVALUATION PHASE ----------------------

       ! Ok, we found the positions of the local matrix entries
       ! that we have to change.
       ! To calculate the matrix contributions, we have to evaluate
       ! the elements to give us the values of the basis functions
       ! in all the DOF`s in all the elements in our set.

       ! Get the element evaluation tag of all FE spaces. We need it to evaluate
       ! the elements later. All of them can be combined with OR, what will give
       ! a combined evaluation tag. 
       cevaluationTag = rlocalMatrixAssembly(1)%cevaluationTag

       ! The cubature points are already initialised by 1D->2D mapping.
       cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

       !      ! Do we have a (multi-)linear transformation?
       !      bisLinearTrafo = trafo_isLinearTrafo(rlocalMatrixAssembly%ctrafoType)
       !
       !      if (bisLinearTrafo) then
       !        ! We need the vertices of the element corners and the number
       !        ! of vertices per element to compute the length of the element
       !        ! edge at the boundary
       !        cevaluationTag = ior(cevaluationTag, EL_EVLTAG_COORDS)
       !        nve = trafo_igetNVE(rlocalMatrixAssembly%ctrafoType)
       !      end if

       ! Calculate all information that is necessary to evaluate the finite element
       ! on all cells of our subset. This includes the coordinates of the points
       ! on the cells.
       !      call elprep_prepareSetForEvaluation (p_revalElementSet,&
       !          cevaluationTag, rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
       !          IelementList(IELset:IELmax), rlocalMatrixAssembly%ctrafoType, &
       !          DpointsRef=DpointsRef)
       !      p_Dcoords => p_revalElementSet%p_Dcoords

       call elprep_prepareSetForEvaluation (&
            rlocalMatrixAssembly(1)%revalElementSet,&
            cevaluationTag,  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation, &
            IelementList(1,IELset:IELmax), rlocalMatrixAssembly(1)%ctrafoType, &
            DpointsRef=DpointsRef(:,:,:,1))
       call elprep_prepareSetForEvaluation (&
            rlocalMatrixAssembly(2)%revalElementSet,&
            cevaluationTag,  rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%p_rtriangulation, &
            IelementList(3,IELset:IELmax), rlocalMatrixAssembly(2)%ctrafoType, &
            DpointsRef=DpointsRef(:,:,:,2))




       !      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
       !      if (.not. rlocalMatrixAssembly%rform%ballCoeffConstant) then
       !        if (present(fcoeff_buildMatrixScBdr2D_sim)) then
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

       ! If the flux function needs other, than just the function values from the solution
       ! (for example the derivatives), we will give an evalElementSet to it
       ! This is filled here

       call domint_initIntegrationByEvalSet (rlocalMatrixAssembly(1)%revalElementSet, rintSubset(1))
       call domint_initIntegrationByEvalSet (rlocalMatrixAssembly(2)%revalElementSet, rintSubset(2))
       !rintSubset(1)%ielementDistribution = 0
       rintSubset(1)%ielementStartIdx = IELset
       rintSubset(1)%p_Ielements => IelementList(1,IELset:IELmax)
       rintSubset(1)%p_IdofsTrial => rlocalMatrixAssembly(1)%p_IdofsTest
       rintSubset(1)%celement = rlocalMatrixAssembly(1)%celementTest
       !rintSubset(2)%ielementDistribution = 0
       rintSubset(2)%ielementStartIdx = IELset
       rintSubset(2)%p_Ielements => IelementList(2,IELset:IELmax)
       rintSubset(2)%p_IdofsTrial => rlocalMatrixAssembly(2)%p_IdofsTest
       rintSubset(2)%celement = rlocalMatrixAssembly(2)%celementTest


       call flux_dg_buildMatrixBlEdge2D_sim_ss (&
            !            rlocalVectorAssembly(1)%p_Dcoefficients(1,:,1:IELmax-IELset+1),&
            !            DsolVals(:,:,1:IELmax-IELset+1),&
       DfluxValues(:,:,:,:,1:IELmax-IELset+1),&
            rvectorSol,&
            IelementList(2,IELset:IELmax),&
            p_Dside,&
            raddTriaData%p_Dnormals(:,Iedgelist(IELset:IELmax)),&
            !DpointsReal(1:ndim2d,1:ncubp,1:IELmax-IELset+1),&
       rintSubset,&
            rcollection )


       call domint_doneIntegration (rintSubset(1))
       call domint_doneIntegration (rintSubset(2))





       ! Calculate the values of the basis functions.
       !      call elem_generic_sim2 (rlocalMatrixAssembly%celementTest, &
       !          p_revalElementSet, rlocalMatrixAssembly%BderTest, &
       !          rlocalMatrixAssembly%p_DbasTest) 
       call elem_generic_sim2 (rlocalMatrixAssembly(1)%celementTest, &
            rlocalMatrixAssembly(1)%revalElementSet,&
            rlocalMatrixAssembly(1)%BderTest, &
            rlocalMatrixAssembly(1)%p_DbasTest)
       call elem_generic_sim2 (rlocalMatrixAssembly(2)%celementTest, &
            rlocalMatrixAssembly(2)%revalElementSet,&
            rlocalMatrixAssembly(2)%BderTest, &
            rlocalMatrixAssembly(2)%p_DbasTest)

       !      ! Omit the calculation of the trial function values if they
       !      ! are identical to the test function values.
       !      if (.not. rlocalMatrixAssembly%bidenticalTrialAndTest) then
       !        call elem_generic_sim2 (rlocalMatrixAssembly%celementTrial, &
       !            p_revalElementSet, rlocalMatrixAssembly%BderTrial, &
       !            rlocalMatrixAssembly%p_DbasTrial)
       !      end if

       !      ! Calculate the length of egdes on the boundary. Depending on
       !      ! whether the transformation is (multi-)linear or not we compute
       !      ! the edge length as the distance between the two corner
       !      ! vertices of the element located on the boundary or as the real
       !      ! length of the boundary segment of the element.
       !      !
       !      ! The length of the current edge serves as a "determinant" in
       !      ! the cubature, so we have to divide it by 2 as an edge on the
       !      ! unit interval [-1,1] has length 2.
       !      if (bisLinearTrafo) then
       !        do iel = 1,IELmax-IELset+1
       !          DedgeLength(iel) = 0.5_DP*sqrt(&
       !              (p_Dcoords(1,    IelementOrientation(IELset+iel-1),iel)-&
       !               p_Dcoords(1,mod(IelementOrientation(IELset+iel-1),nve)+1,iel))**2+&
       !              (p_Dcoords(2,    IelementOrientation(IELset+iel-1),iel)-&
       !               p_Dcoords(2,mod(IelementOrientation(IELset+iel-1),nve)+1,iel))**2)
       !        end do
       !      else
       !        do iel = 1,IELmax-IELset+1
       !          DedgeLength(iel) = 0.5_DP*(DedgePosition(2,IELset+iel-1)-&
       !                                     DedgePosition(1,IELset+iel-1))
       !        end do
       !      end if

       ! --------------------- DOF COMBINATION PHASE ------------------------

       ! Values of all basis functions calculated. Now we can start 
       ! to integrate!

       ! Clear the local matrices
       p_Dentryii(:,:,:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryai(:,:,:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryia(:,:,:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryaa(:,:,:,:,1:IELmax-IELset+1) = 0.0_DP

       !      p_Dentryii = 0.0_DP
       !      p_Dentryai = 0.0_DP
       !      p_Dentryia = 0.0_DP
       !      p_Dentryaa = 0.0_DP

       ! We have two different versions for the integration - one
       ! with constant coefficients and one with nonconstant coefficients.
       !
       ! Check the bilinear form which one to use:

       !      if (rlocalMatrixAssembly%rform%ballCoeffConstant) then
       !      
       !        ! Constant coefficients. The coefficients are to be found in
       !        ! the Dcoefficients variable of the form.
       !        !
       !        ! Loop over the elements in the current set.
       !
       !        do iel = 1,IELmax-IELset+1
       !          
       !          ! Get the length of the edge.
       !          dlen = DedgeLength(iel)
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
       !              ! Now loop through all possible combinations of DOF`s
       !              ! in the current cubature point. The outer loop
       !              ! loops through the "O"`s in the above picture,
       !              ! the test functions:
       !
       !              do idofe = 1,indofTest
       !              
       !                ! Get the value of the (test) basis function 
       !                ! phi_i (our "O") in the cubature point:
       !                db = p_DbasTest(idofe,ib,icubp,iel)
       !                
       !                ! Perform an inner loop through the other DOF`s
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

       ! Nonconstant coefficients. The coefficients are to be found in
       ! the Dcoefficients variable as computed above.
       !
       ! Loop over the elements.

       do iel = 1,IELmax-IELset+1

          ! Get the length of the edge.
          !dlen = DedgeLength(iel)
          dlen = 0.5_DP*raddTriaData%p_Dedgelength(Iedgelist(IELset+iel-1))

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

             ! calculate the current weighting factor in the cubature formula
             ! in that cubature point.
             !            domega = dlen * p_Domega(icubp)
             domega1 = dlen * rlocalMatrixAssembly(1)%p_Domega(icubp)
             domega2 = dlen * rlocalMatrixAssembly(2)%p_Domega(icubp)

             ! Loop over the additive factors in the bilinear form.
             do ialbet = 1,rlocalMatrixAssembly(1)%rform%itermcount
             
                ! Test if we only need the trial functions on one side of the edge
                if ((p_Dside(1,ialbet)==1.0_dp).and.(p_Dside(2,ialbet)==0.0_dp)) then
                ! Trialfunction only needed on first side

                    ! Get from Idescriptors the type of the derivatives for the 
                    ! test and trial functions. The summand we calculate
                    ! here will be added to the matrix entry:
                    !
                    ! a_ij  =  int_... ( psi_j )_ia  *  ( phi_i )_ib
                    !
                    ! -> Ix=0: function value, 
                    !      =1: first derivative, ...
                    !    as defined in the module 'derivative'.

                    ia = rlocalMatrixAssembly(1)%rform%Idescriptors(1,ialbet)
                    ib = rlocalMatrixAssembly(1)%rform%Idescriptors(2,ialbet)

                    ! Multiply domega with the coefficient of the form.
                    ! This gives the actual value to multiply the
                    ! function value with before summing up to the integral.
                    ! Get the precalculated coefficient from the coefficient array.
                    daux1 = domega1 * DfluxValues(:,:,ialbet,icubp,iel)
                    daux2 = domega2 * DfluxValues(:,:,ialbet,icubp,iel) * (-1.0_dp)

                    ! Now loop through all possible combinations of DOF`s
                    ! in the current cubature point. The outer loop
                    ! loops through the "O" in the above picture,
                    ! the test functions:

                    do idofe = 1,indofTest

                       ! Get the value of the (test) basis function 
                       ! phi_i (our "O") in the cubature point:
                       db1 = rlocalMatrixAssembly(1)%p_DbasTest(idofe,ib,icubp,iel)
                       db2 = rlocalMatrixAssembly(2)%p_DbasTest(idofe,ib,icubp,iel)

                       ! Perform an inner loop through the other DOF`s
                       ! (the "X"). 

                       do jdofe = 1,indofTrial

                          ! Get the value of the basis function 
                          ! psi_j (our "X") in the cubature point. 
                          ! Them multiply:
                          !    db * dbas(..) * daux
                          ! ~= phi_i * psi_j * coefficient * cub.weight
                          ! Summing this up gives the integral, so the contribution
                          ! to the global matrix. 
                          !
                          ! Simply summing up db * dbas(..) * daux would give
                          ! the coefficient of the local matrix. We save this
                          ! contribution in the local matrix of element iel.

                          
                          do iblock = 1, nvar
                             do jblock = 1, nvar
                                ! Testfunction on the 'first' (i) side
                                ! Attention: The p_DbasTest are really the trialfunctions (they are the same)
                                p_Dentryii(iblock,jblock,jdofe,idofe,iel) = &
                                     p_Dentryii(iblock,jblock,jdofe,idofe,iel)+db1*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux1(iblock,jblock)*p_Dside(1,ialbet)

                                ! Testfunction on the 'second' (a) side
                                p_Dentryia(iblock,jblock,jdofe,idofe,iel) = &
                                     p_Dentryia(iblock,jblock,jdofe,idofe,iel)+db2*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux2(iblock,jblock)*p_Dside(1,ialbet)
                             end do
                          end do
                          
                       end do ! idofe

                    end do ! jdofe                
                
                elseif ((p_Dside(1,ialbet)==0.0_dp).and.(p_Dside(2,ialbet)==1.0_dp)) then
                ! Trialfunction only needed on the second side

                    ! Get from Idescriptors the type of the derivatives for the 
                    ! test and trial functions. The summand we calculate
                    ! here will be added to the matrix entry:
                    !
                    ! a_ij  =  int_... ( psi_j )_ia  *  ( phi_i )_ib
                    !
                    ! -> Ix=0: function value, 
                    !      =1: first derivative, ...
                    !    as defined in the module 'derivative'.

                    ia = rlocalMatrixAssembly(1)%rform%Idescriptors(1,ialbet)
                    ib = rlocalMatrixAssembly(1)%rform%Idescriptors(2,ialbet)

                    ! Multiply domega with the coefficient of the form.
                    ! This gives the actual value to multiply the
                    ! function value with before summing up to the integral.
                    ! Get the precalculated coefficient from the coefficient array.
                    daux1 = domega1 * DfluxValues(:,:,ialbet,icubp,iel)
                    daux2 = domega2 * DfluxValues(:,:,ialbet,icubp,iel) * (-1.0_dp)

                    ! Now loop through all possible combinations of DOF`s
                    ! in the current cubature point. The outer loop
                    ! loops through the "O" in the above picture,
                    ! the test functions:

                    do idofe = 1,indofTest

                       ! Get the value of the (test) basis function 
                       ! phi_i (our "O") in the cubature point:
                       db1 = rlocalMatrixAssembly(1)%p_DbasTest(idofe,ib,icubp,iel)
                       db2 = rlocalMatrixAssembly(2)%p_DbasTest(idofe,ib,icubp,iel)

                       ! Perform an inner loop through the other DOF`s
                       ! (the "X"). 

                       do jdofe = 1,indofTrial

                          ! Get the value of the basis function 
                          ! psi_j (our "X") in the cubature point. 
                          ! Them multiply:
                          !    db * dbas(..) * daux
                          ! ~= phi_i * psi_j * coefficient * cub.weight
                          ! Summing this up gives the integral, so the contribution
                          ! to the global matrix. 
                          !
                          ! Simply summing up db * dbas(..) * daux would give
                          ! the coefficient of the local matrix. We save this
                          ! contribution in the local matrix of element iel.

                          
                          do iblock = 1, nvar
                             do jblock = 1, nvar
                                ! Testfunction on the 'first' (i) side
                                ! Attention: The p_DbasTest are really the trialfunctions (they are the same)
                                p_Dentryai(iblock,jblock,jdofe,idofe,iel) = &
                                     p_Dentryai(iblock,jblock,jdofe,idofe,iel)+db1*rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)*daux1(iblock,jblock)*p_Dside(2,ialbet)

                                ! Testfunction on the 'second' (a) side
                                p_Dentryaa(iblock,jblock,jdofe,idofe,iel) = &
                                     p_Dentryaa(iblock,jblock,jdofe,idofe,iel)+db2*rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)*daux2(iblock,jblock)*p_Dside(2,ialbet)
                             end do
                          end do
                          
                       end do ! idofe

                    end do ! jdofe                
                
                else ! Trialfunctions are needed on both sides
                
                    ! Get from Idescriptors the type of the derivatives for the 
                    ! test and trial functions. The summand we calculate
                    ! here will be added to the matrix entry:
                    !
                    ! a_ij  =  int_... ( psi_j )_ia  *  ( phi_i )_ib
                    !
                    ! -> Ix=0: function value, 
                    !      =1: first derivative, ...
                    !    as defined in the module 'derivative'.

                    ia = rlocalMatrixAssembly(1)%rform%Idescriptors(1,ialbet)
                    ib = rlocalMatrixAssembly(1)%rform%Idescriptors(2,ialbet)

                    ! Multiply domega with the coefficient of the form.
                    ! This gives the actual value to multiply the
                    ! function value with before summing up to the integral.
                    ! Get the precalculated coefficient from the coefficient array.
                    daux1 = domega1 * DfluxValues(:,:,ialbet,icubp,iel)
                    daux2 = domega2 * DfluxValues(:,:,ialbet,icubp,iel) * (-1.0_dp)

                    ! Now loop through all possible combinations of DOF`s
                    ! in the current cubature point. The outer loop
                    ! loops through the "O" in the above picture,
                    ! the test functions:

                    do idofe = 1,indofTest

                       ! Get the value of the (test) basis function 
                       ! phi_i (our "O") in the cubature point:
                       db1 = rlocalMatrixAssembly(1)%p_DbasTest(idofe,ib,icubp,iel)
                       db2 = rlocalMatrixAssembly(2)%p_DbasTest(idofe,ib,icubp,iel)

                       ! Perform an inner loop through the other DOF`s
                       ! (the "X"). 

                       do jdofe = 1,indofTrial

                          ! Get the value of the basis function 
                          ! psi_j (our "X") in the cubature point. 
                          ! Them multiply:
                          !    db * dbas(..) * daux
                          ! ~= phi_i * psi_j * coefficient * cub.weight
                          ! Summing this up gives the integral, so the contribution
                          ! to the global matrix. 
                          !
                          ! Simply summing up db * dbas(..) * daux would give
                          ! the coefficient of the local matrix. We save this
                          ! contribution in the local matrix of element iel.

                          
                          do iblock = 1, nvar
                             do jblock = 1, nvar
                                ! Testfunction on the 'first' (i) side
                                ! Attention: The p_DbasTest are really the trialfunctions (they are the same)
                                p_Dentryii(iblock,jblock,jdofe,idofe,iel) = &
                                     p_Dentryii(iblock,jblock,jdofe,idofe,iel)+db1*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux1(iblock,jblock)*p_Dside(1,ialbet)
                                p_Dentryai(iblock,jblock,jdofe,idofe,iel) = &
                                     p_Dentryai(iblock,jblock,jdofe,idofe,iel)+db1*rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)*daux1(iblock,jblock)*p_Dside(2,ialbet)

                                ! Testfunction on the 'second' (a) side
                                p_Dentryia(iblock,jblock,jdofe,idofe,iel) = &
                                     p_Dentryia(iblock,jblock,jdofe,idofe,iel)+db2*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux2(iblock,jblock)*p_Dside(1,ialbet)
                                p_Dentryaa(iblock,jblock,jdofe,idofe,iel) = &
                                     p_Dentryaa(iblock,jblock,jdofe,idofe,iel)+db2*rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)*daux2(iblock,jblock)*p_Dside(2,ialbet)
                             end do
                          end do
                          
                       end do ! idofe

                    end do ! jdofe
                
                end if
             
             end do ! ialbet

          end do ! icubp 

       end do ! iel

       !      end if ! rform%ballCoeffConstant

       ! Incorporate the local matrices into the global one.
       ! Kentry gives the position of the additive contributions in Dentry.
       !
       ! OpenMP-Extension: This is a critical section. Only one thread is
       ! allowed to write to the matrix, otherwise the matrix may get
       ! messed up.
       ! The critical section is put around both loops as indofTest/indofTrial
       ! are usually small and quickly to handle.

       !      if (cconstrType .eq. BILF_MATC_LUMPED) then
       !
       !        !$omp critical
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
       !        !$omp end critical
       !
       !      else




       !      call lsysbl_getbase_double(rmatrix,p_Da)

       !      ! Get pointers to the data entries of the block matrix
       !      do iblock = 1, nvar
       !        do jblock = 1,nvar
       !          call lsyssc_getbase_double(rmatrix%RmatrixBlock(iblock,jblock),p_matrixBlockDataPointers(iblock,jblock)%Da)
       !        end do
       !      end do

       !p_Ddata(rvector%RvectorBlock(ivar)%iidxFirstEntry+rlocalVectorAssembly(1)%p_Idofs(idofe,iel)-1) = &            
       !              p_Ddata(rvector%RvectorBlock(ivar)%iidxFirstEntry+rlocalVectorAssembly(1)%p_Idofs(idofe,iel)-1) + &
       !              DlocalData(ivar,1,idofe)






       !$omp critical
       do iblock = 1, nvar
          do jblock = 1, nvar
             call lsyssc_getbase_double(rmatrix%RmatrixBlock(iblock,jblock),p_Da)
             do iel = 1,IELmax-IELset+1
             
             if (IelementList(2,IELset+iel-1).ne.0) then
             ! Not at the boundary - write all DOFs
                do idofe = 1,indofTest
                   do jdofe = 1,indofTrial
                      
                      p_Da(p_Kentryii(jdofe,idofe,iel)) = &
                           p_Da(p_Kentryii(jdofe,idofe,iel)) + p_Dentryii(iblock,jblock,jdofe,idofe,iel)

                      p_Da(p_Kentryia(jdofe,idofe,iel)) = &
                           p_Da(p_Kentryia(jdofe,idofe,iel)) + p_Dentryia(iblock,jblock,jdofe,idofe,iel)

                      p_Da(p_Kentryai(jdofe,idofe,iel)) = &
                           p_Da(p_Kentryai(jdofe,idofe,iel)) + p_Dentryai(iblock,jblock,jdofe,idofe,iel)!*real(min(1,IelementList(2,IELset+iel-1)))
                      p_Da(p_Kentryaa(jdofe,idofe,iel)) = &
                           p_Da(p_Kentryaa(jdofe,idofe,iel)) + p_Dentryaa(iblock,jblock,jdofe,idofe,iel)!*real(min(1,IelementList(2,IELset+iel-1)))

                   end do
                end do
             else
             ! At the boundary - write only inner DOFs
                do idofe = 1,indofTest
                   do jdofe = 1,indofTrial
                      
                      p_Da(p_Kentryii(jdofe,idofe,iel)) = &
                           p_Da(p_Kentryii(jdofe,idofe,iel)) + p_Dentryii(iblock,jblock,jdofe,idofe,iel)

                   end do
                end do             
             end if

             end do ! iel
          end do ! jblock
       end do ! iblock
       !$omp end critical

       !      end if

    end do ! IELset
    !$omp end do

    ! Release the local matrix assembly structure
    call bilf_releaseAssemblyData(rlocalMatrixAssembly(1))
    call bilf_releaseAssemblyData(rlocalMatrixAssembly(2))

    ! Deallocate memory
    deallocate(Dxi2D, DpointsRef) !, DpointsPar, DedgeLength)
    deallocate(p_Kentryii,p_Kentryia,p_Kentryai,p_Kentryaa)
    deallocate(p_Dentryii,p_Dentryia,p_Dentryai,p_Dentryaa)
    deallocate(p_Dside)
    deallocate(IelementList)
    deallocate(daux1,daux2)

    !$omp end parallel

  end subroutine dg_bilf_assembleSubmeshMat9Bdr2D_Block_ss
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
   !****************************************************************************

  !<subroutine>  

  subroutine dg_makeSolPos (rvectorBlock)

    !<description>

    ! Checks for negative solution values and fixes the problem

    !</description>

    !<inputoutput>

    ! The solution vector. Used to calculate the solution on the edges.
    type(t_vectorBlock), intent(inout) :: rvectorBlock

    !</inputoutput>

    !</subroutine>
    
    
    
    ! Local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: iel, nvar, ivar, nlDOF
    integer, dimension(:,:), allocatable :: IdofGlob
    integer, dimension(:), allocatable :: IelList
    real(dp) :: da, db

    ! Get pointer to the data of the (solution) vector
    call lsysbl_getbase_double (rvectorBlock, p_Ddata)
    
    ! Get number of local DOFs
    nlDOF = elem_igetNDofLoc(rvectorBlock%p_rblockDiscr%RspatialDiscr(1)%RelementDistr(1)%celement)
    
    ! Allocate space for DOFs
    allocate(IdofGlob(nlDOF,rvectorBlock%p_rblockDiscr%p_rtriangulation%NEL))
    
    ! Allocate space for list of elements
    allocate(IelList(rvectorBlock%p_rblockDiscr%p_rtriangulation%NEL))
    
    forall (iel=1:rvectorBlock%p_rblockDiscr%p_rtriangulation%NEL) IelList(iel) = iel
    
    ! Get DOFs
    call dof_locGlobMapping_mult(rvectorBlock%p_rblockDiscr%RspatialDiscr(1), IelList, IdofGlob)

    ! Get number of variables of the system
    nvar = rvectorBlock%nblocks

    do ivar = 1, nvar
       if ((ivar == 2).or.(ivar==3)) cycle
       do iel = 1, rvectorBlock%p_rblockDiscr%p_rtriangulation%NEL
          ! Check value in midpoint
          if (p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1,iel)-1)<20.0_dp*SYS_EPSREAL_SP) then
             p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1,iel)-1)=20.0_dp*SYS_EPSREAL_SP
          end if
          
          
          da = (abs(p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2,iel)-1)) + &
                abs(p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3,iel)-1)))
          db = p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(1,iel)-1)
          
          if ( (da/db).ge.0.99_dp ) then
             p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2,iel)-1) = &
                p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(2,iel)-1) *db/da*0.99_dp
             p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3,iel)-1) = &
                p_Ddata(rvectorBlock%RvectorBlock(ivar)%iidxFirstEntry+IdofGlob(3,iel)-1) *db/da*0.99_dp
          end if      
          
       end do
    end do
      
    ! Deallocate space for DOFs and the elementlist
    deallocate(IdofGlob,IelList)
    
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !****************************************************************************

  !<subroutine>

  subroutine bilf_dg_buildMatrixScEdge2D_ss (rform, ccubType, bclear, rmatrix,&
       rvectorSol, raddTriaData,&
       flux_dg_buildMatrixScEdge2D_sim,&
       rcollection, cconstrType)

    !<description>
    ! This routine calculates the entries of a finite element matrix in 2D.
    ! The matrix structure must be prepared with bilf_createMatrixStructure
    ! in advance.
    ! In case the array for the matrix entries does not exist, the routine
    ! allocates memory in size of the matrix of the heap for the matrix entries.
    !
    ! For setting up the entries, the discretisation structure attached to
    ! the matrix is used (rmatrix%p_rdiscretisation). This is
    ! normally attached to the matrix by bilf_createMatrixStructure.
    !
    ! The matrix must be unsorted when this routine is called, 
    ! otherwise an error is thrown.
    !</description>

    !<input>
    ! The bilinear form specifying the underlying PDE of the discretisation.
    type(t_bilinearForm), intent(in) :: rform

    ! A line cubature formula CUB_xxxx_1D to be used for line integration.
    integer(I32), intent(in) :: ccubType

    ! The solution vector. Used to calculate the solution on the edges.
    type(t_vectorScalar), intent(in) :: rvectorSol

    ! Additional triangulation data
    type(t_additionalTriaData), intent(in) :: raddTriaData

    ! Whether to clear the matrix before calculating the entries.
    ! If .FALSE., the new matrix entries are added to the existing entries.
    logical, intent(in) :: bclear

    ! A callback routine for the flux function.
    include 'intf_flux_dg_buildMatrixScEdge2D_ss.inc'
    optional :: flux_dg_buildMatrixScEdge2D_sim

    ! OPTIONAL: One of the BILF_MATC_xxxx constants that allow to specify
    ! the matrix construction method. If not specified,
    ! BILF_MATC_EDGEBASED is used.
    integer, intent(in), optional :: cconstrType
    !</input>

    !<inputoutput>
    ! The FE matrix. Calculated matrix entries are imposed to this matrix.
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: A collection structure. This structure is given to the
    ! callback function for nonconstant coefficients to provide additional
    ! information. 
    type(t_collection), intent(inout), target, optional :: rcollection
    !</inputoutput>

    !</subroutine>

    ! local variables
    type(t_bilfMatrixAssembly), dimension(2) :: rmatrixAssembly
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: IelementList, p_IedgeList
    integer :: ccType
    integer :: iedge, ielementDistr

    ! The matrix must be unsorted, otherwise we can not set up the matrix.
    ! Note that we cannot switch off the sorting as easy as in the case
    ! of a vector, since there is a structure behind the matrix! So the caller
    ! has to make sure, the matrix is unsorted when this routine is called.
    if (rmatrix%isortStrategy .gt. 0) then
       call output_line ('Matrix-structure must be unsorted!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ! The matrix must provide discretisation structures
    if ((.not. associated(rmatrix%p_rspatialDiscrTest)) .or. &
         (.not. associated(rmatrix%p_rspatialDiscrTrial))) then
       call output_line ('No discretisation associated!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ! The discretisation must provide a triangulation structure
    if ((.not. associated(rmatrix%p_rspatialDiscrTest%p_rtriangulation)) .or. &
         (.not. associated(rmatrix%p_rspatialDiscrTrial%p_rtriangulation))) then
       call output_line('No triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    !  ! The discretisation must provide a boundary structure
    !  if ((.not. associated(rmatrix%p_rspatialDiscrTest%p_rboundary)) .or. &
    !      (.not. associated(rmatrix%p_rspatialDiscrTrial%p_rboundary))) then
    !    call output_line('No boundary associated!',&
    !        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
    !    call sys_halt()
    !  end if

    !  ! Set pointers for quicker access
    !  p_rboundary => rmatrix%p_rspatialDiscrTest%p_rboundary
    !  if (.not.associated(p_rboundary, rmatrix%p_rspatialDiscrTrial%p_rboundary)) then
    !    call output_line('Invalid boundary associated!',&
    !        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
    !    call sys_halt()
    !  end if

    p_rtriangulation => rmatrix%p_rspatialDiscrTest%p_rtriangulation
    if (.not.associated(p_rtriangulation, rmatrix%p_rspatialDiscrTrial%p_rtriangulation)) then
       call output_line('Invalid triangulation associated!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end if

    ccType = BILF_MATC_EDGEBASED
    if (present(cconstrType)) ccType = cconstrType

    ! Do we have a uniform triangulation? Would simplify a lot...
    select case (rmatrix%p_rspatialDiscrTest%ccomplexity)
    case (SPDISC_UNIFORM,SPDISC_CONFORMAL) 
       ! Uniform and conformal discretisations
       select case (rmatrix%cdataType)
       case (ST_DOUBLE) 
          ! Which matrix structure do we have?
          select case (rmatrix%cmatrixFormat) 
          case (LSYSSC_MATRIX9)

             ! Probably allocate/clear the matrix
             !        if (rmatrix%h_DA .eq. ST_NOHANDLE) then
             !          call lsyssc_allocEmptyMatrix(rmatrix,LSYSSC_SETM_ZERO)
             !        else
             if (bclear) call lsyssc_clearMatrix (rmatrix)
             !        end if


             ! Allocate the edgelist
             allocate(p_IedgeList(rmatrix%p_rspatialDiscrTest%p_rtriangulation%NMT))

             ! All edges
             forall (iedge = 1:rmatrix%p_rspatialDiscrTest%p_rtriangulation%NMT) p_IedgeList(iedge)=iedge

             ! Initialise a matrix assembly structure for that element distribution
             ielementDistr = 1
             call bilf_initAssembly(rmatrixAssembly(1),rform,&
                  rmatrix%p_rspatialDiscrTest%RelementDistr(ielementDistr)%celement,&
                  rmatrix%p_rspatialDiscrTrial%RelementDistr(ielementDistr)%celement,&
                  ccubType, BILF_NELEMSIM)

             ! Do the same for the other side of the egde
             ielementDistr = 1
             call dg_bilf_initAssembly_reverseCubPoints(rmatrixAssembly(2),rform,&
                  rmatrix%p_rspatialDiscrTest%RelementDistr(ielementDistr)%celement,&
                  rmatrix%p_rspatialDiscrTrial%RelementDistr(ielementDistr)%celement,&
                  ccubType, BILF_NELEMSIM)

             ! Assemble the data for all elements in this element distribution
             call dg_bilf_assembleSubmeshMat9Bdr2D_ss (rmatrixAssembly, rmatrix,&
                  rvectorSol, raddTriaData, p_IedgeList(1:rmatrix%p_rspatialDiscrTest%p_rtriangulation%NMT),&
                  ccType, flux_dg_buildMatrixScEdge2D_sim, rcollection)


             ! Release the assembly structure.
             call bilf_doneAssembly(rmatrixAssembly(1))
             call bilf_doneAssembly(rmatrixAssembly(2))

             ! Deallocate the edgelist
             deallocate(p_IedgeList)



          case default
             call output_line ('Not supported matrix structure!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
             call sys_halt()
          end select

       case default
          call output_line ('Single precision matrices currently not supported!', &
               OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
          call sys_halt()
       end select

    case default
       call output_line ('General discretisation not implemented!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalarBdr2D')
       call sys_halt()
    end select

  end subroutine bilf_dg_buildMatrixScEdge2D_ss



  !****************************************************************************

  !<subroutine>  

  subroutine dg_bilf_assembleSubmeshMat9Bdr2D_ss (rmatrixAssembly, rmatrix,&
       rvectorSol, raddTriaData, IedgeList,&
       cconstrType, flux_dg_buildMatrixScEdge2D_sim, rcollection)

    !<description>

    ! Assembles the matrix entries for a submesh by integrating over the
    ! boundary region in 2D.

    !</description>

    !<input>

    ! List of elements where to assemble the bilinear form.
    integer, dimension(:), intent(in), target :: IedgeList

    ! One of the BILF_MATC_xxxx constants that allow to specify the
    ! matrix construction method.
    integer, intent(in) :: cconstrType

    ! The solution vector. Used to calculate the solution on the edges.
    type(t_vectorScalar), intent(in) :: rvectorSol

    ! Additional triangulation data
    type(t_additionalTriaData), intent(in) :: raddTriaData

    ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
    ! Must be present if the matrix has nonconstant coefficients!
    include 'intf_flux_dg_buildMatrixScEdge2D_ss.inc'
    optional :: flux_dg_buildMatrixScEdge2D_sim

    !</input>

    !<inputoutput>

    ! A matrix assembly structure prepared with bilf_initAssembly.
    type(t_bilfMatrixAssembly), intent(inout), dimension(2), target :: rmatrixAssembly

    ! A matrix where to assemble the contributions to.
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: A pointer to a collection structure. This structure is given to the
    ! callback function for nonconstant coefficients to provide additional
    ! information. 
    type(t_collection), intent(inout), target, optional :: rcollection

    !</inputoutput>

    !</subroutine>

    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_DA
    integer :: indofTest,indofTrial,ncubp

    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax,ibdc,k
    integer :: iel,icubp,ialbet,ia,ib,idofe,jdofe,nve
    real(DP) :: domega1,domega2,daux1,daux2,db1,db2,dlen
    integer(I32) :: cevaluationTag
    type(t_bilfMatrixAssembly), dimension(2), target :: rlocalMatrixAssembly
    type(t_domainIntSubset), dimension(2) :: rintSubset
    integer, dimension(:,:,:), pointer :: p_Kentryii, p_Kentryia, p_Kentryai, p_Kentryaa
    real(DP), dimension(:,:,:), pointer :: p_Dentryii, p_Dentryia, p_Dentryai, p_Dentryaa
    real(DP), dimension(:,:,:), pointer :: p_Dcoords
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:), pointer :: p_Dside
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial
    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
    real(DP), dimension(:), pointer :: p_DcoefficientsBilf
    integer, dimension(:,:), pointer :: p_IdofsTest
    integer, dimension(:,:), pointer :: p_IdofsTrial
    type(t_evalElementSet), pointer :: p_revalElementSet
    integer, dimension(:,:),pointer :: p_Idescriptors
    integer, dimension(:,:), allocatable, target :: IelementList

    ! Pointer to Ielementsatedge in the triangulation
    integer, dimension(:,:), pointer :: p_IelementsAtEdge

    ! Pointer to IverticesAtEdge in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtEdge

    ! Pointer to the vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Pointer to IverticesAtEelement in the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Space for the values of the flux function
    real(DP), dimension(:,:,:), allocatable :: DfluxValues

    ! Arrays for cubature points 1D->2D
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D_1, Dxi1D_2
    real(DP), dimension(:,:,:,:), allocatable :: Dxi2D,DpointsRef
    real(DP), dimension(:,:), allocatable :: DpointsPar
    real(DP), dimension(:), allocatable :: DedgeLength

    integer(i32) :: icoordSystem
    integer :: NEL
    integer :: iside
    logical :: bisLinearTrafo

    !    ! Boundary component?
    !    ibdc = rboundaryRegion%iboundCompIdx

    ! Get some pointers for faster access
    call lsyssc_getbase_double (rmatrix,p_DA)
    indofTest = rmatrixAssembly(1)%indofTest
    indofTrial = rmatrixAssembly(1)%indofTrial
    ncubp = rmatrixAssembly(1)%ncubp

    ! Open-MP-Extension: Copy the matrix assembly data to the local
    ! matrix assembly data, where we can allocate memory.
    !
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the matrix
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original matrix assembly
    ! stucture or disturbing the data of the other processors.
    !
    !$omp parallel default(shared) &
    !$omp private(DedgeLength,DpointsPar,DpointsRef,Dxi1D,Dxi2D,IELmax,bisLinearTrafo,&
    !$omp         cevaluationTag,daux,db,dlen,domega,ia,ialbet,ib,icoordSystem,icubp,&
    !$omp         idofe,iel,jdofe,k,p_DbasTest,p_DbasTrial,p_Dcoefficients,&
    !$omp         p_DcoefficientsBilf,p_Dcoords,p_DcubPtsRef,p_Dentry,p_Domega,&
    !$omp         p_Idescriptors,p_IdofsTest,p_IdofsTrial,p_Kentry,&
    !$omp         p_revalElementSet,rintSubset,rlocalMatrixAssembly)
    rlocalMatrixAssembly(1) = rmatrixAssembly(1)
    rlocalMatrixAssembly(2) = rmatrixAssembly(2)
    call bilf_allocAssemblyData(rlocalMatrixAssembly(1))
    call bilf_allocAssemblyData(rlocalMatrixAssembly(2))

    ! Allocate space for the positions of the DOFs in the matrix
    allocate(p_Kentryii(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryia(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryai(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Kentryaa(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))

    ! Allocate space for the coefficient of the solutions DOFs on each side of the edge
    allocate(p_Dside(2,rmatrixAssembly(1)%rform%itermCount))

    ! Allocate space for the entries in the local matrices
    allocate(p_Dentryii(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryia(rmatrixAssembly(1)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryai(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(1)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))
    allocate(p_Dentryaa(rmatrixAssembly(2)%indofTrial,&
         rmatrixAssembly(2)%indofTest,rmatrixAssembly(1)%nelementsPerBlock))



    ! Allocate space for the flux variables DIM(nvar,ialbet,ncubp,elementsperblock)
    allocate(DfluxValues(rlocalMatrixAssembly(1)%rform%itermcount,ncubp,rlocalMatrixAssembly(1)%nelementsPerBlock))

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


    ! Get number of elements
    NEL = rmatrix%p_rspatialDiscrTest%p_rtriangulation%NEL

    ! Get pointers to elements at edge
    call storage_getbase_int2D(&
         rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_IelementsAtEdge,&
         p_IelementsAtEdge)

    ! Get pointers to the vertex coordinates
    call storage_getbase_double2D(&
         rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_DvertexCoords,&
         p_DvertexCoords)

    ! Get pointers to vertices at edge
    call storage_getbase_int2D(&
         rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtEdge,&
         p_IverticesAtEdge)

    ! Get pointers to vertices at elements
    call storage_getbase_int2D(&
         rmatrix%p_rspatialDiscrTest%p_rtriangulation%h_IverticesAtElement,&
         p_IverticesAtElement)   

    ! Get the elements adjacent to the given edges
    allocate(IelementList(3,size(IedgeList)))
    IelementList(1:2,1:size(IedgeList))=p_IelementsAtEdge(1:2,IedgeList(:))

    ! Copy the second component and replace 0s by 1s
    IelementList(3,size(IedgeList))=IelementList(2,size(IedgeList))
    do iel = 1,size(IedgeList)
       IelementList(3,iel)=max(IelementList(2,iel),1)
    end do

    ! Transpose the coordinate array such that we get coordinates we
    ! can work with in the mapping between 1D and 2D.
    do k = 1, ubound(rlocalmatrixAssembly(1)%p_DcubPtsRef,1)
       do icubp = 1,ncubp
          Dxi1D_1(icubp,k) = rlocalmatrixAssembly(1)%p_DcubPtsRef(k,icubp)
          Dxi1D_2(icubp,k) = rlocalmatrixAssembly(2)%p_DcubPtsRef(k,icubp)
       end do
    end do

    ! Allocate memory for the cubature points in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,2,rlocalMatrixAssembly(1)%nelementsPerBlock))

    ! Allocate memory for the coordinates of the reference points
    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalMatrixAssembly(1)%nelementsPerBlock,2))

    !    ! Allocate memory for the parameter values of the points on the boundary
    !    allocate(DpointsPar(ncubp,rlocalMatrixAssembly%nelementsPerBlock))
    !
    !    ! Allocate memory for the length of edges on the boundary
    !    allocate(DedgeLength(rlocalMatrixAssembly%nelementsPerBlock))

    ! Get the type of coordinate system
    icoordSystem = elem_igetCoordSystem(rlocalMatrixAssembly(1)%celementTrial)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !$omp do schedule(static,1)
    do IELset = 1, size(IedgeList), rlocalMatrixAssembly(1)%nelementsPerBlock

       ! We always handle nelementsPerBlock elements simultaneously.
       ! How many elements have we actually here?
       ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
       ! elements simultaneously.

       IELmax = min(size(IedgeList),IELset-1+rlocalMatrixAssembly(1)%nelementsPerBlock)

       ! Map the 1D cubature points to the edges in 2D.
       do iel = 1,IELmax-IELset+1
          call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(1,Iedgelist(IELset+iel-1)), &
               ncubp, Dxi1D_1, Dxi2D(:,:,1,iel))
          call trafo_mapCubPts1Dto2D(icoordSystem, raddTriaData%p_IlocalEdgeNumber(2,Iedgelist(IELset+iel-1)), &
               ncubp, Dxi1D_2, Dxi2D(:,:,2,iel))
       end do

       !      ! Calculate the parameter values of the points
       !      do iel = 1,IELmax-IELset+1
       !        do icubp = 1,ncubp
       !          ! Dxi1D is in [-1,1] while the current edge has parmeter values
       !          ! [DedgePosition(1),DedgePosition(2)]. So do a linear
       !          ! transformation to transform Dxi1D into that interval, this 
       !          ! gives the parameter values in length parametrisation
       !          call mprim_linearRescale(Dxi1D(icubp,1), -1.0_DP, 1.0_DP,&
       !              DedgePosition(1,IELset+iel-1), DedgePosition(2,IELset+iel-1),&
       !              DpointsPar(icubp,iel))
       !        end do
       !      end do

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
       ! global DOF`s of our BILF_NELEMSIM elements simultaneously.
       !      call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTest, &
       !          IelementList(IELset:IELmax), p_IdofsTest)
       call dof_locGlobMapping_mult( rmatrix%p_rspatialDiscrTest, &
            IelementList(1,IELset:IELmax), rlocalMatrixAssembly(1)%p_IdofsTest)
       call dof_locGlobMapping_mult( rmatrix%p_rspatialDiscrTest, &
            IelementList(3,IELset:IELmax), rlocalMatrixAssembly(2)%p_IdofsTest)


       !      ! If the DOF`s for the trial functions are different, calculate them, too.
       !      if (.not. rlocalMatrixAssembly%bIdenticalTrialAndTest) then
       !        call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTrial, &
       !            IelementList(IELset:IELmax), p_IdofsTrial)
       !      end if

       ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------

       ! For the assembly of the global matrix, we use a "local"
       ! approach. At first we build a "local" system matrix according
       ! to the current element. This contains all additive
       ! contributions of element iel, which are later added at the
       ! right positions to the elements in the global system matrix.
       !
       ! We have indofTrial trial DOF`s per element and
       ! indofTest test DOF`s per element. Therefore there are
       ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
       ! "active" (i.e. have common support) on our current element, each 
       ! giving an additive contribution to the system matrix.
       !
       ! We build a quadratic indofTrial*indofTest local matrix:
       ! Kentry(1..indofTrial,1..indofTest) receives the position 
       ! in the global system matrix, where the corresponding value 
       ! has to be added to.
       ! (The corresponding contributions can be saved separately, 
       ! but we directly add them to the global matrix in this 
       ! approach.)
       !
       ! We build local matrices for all our elements 
       ! in the set simultaneously. Get the positions of the local matrices
       ! in the global matrix.
       !      call bilf_getLocalMatrixIndices (rmatrix,p_IdofsTest,p_IdofsTrial,p_Kentry,&
       !          ubound(p_IdofsTest,1),ubound(p_IdofsTrial,1),IELmax-IELset+1)  
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
       !            rlocalMatrixAssembly(1)%p_IdofsTrial, p_Kentryii,&
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
       !            rlocalMatrixAssembly(2)%p_IdofsTrial, p_Kentryia,&
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
       !            rlocalMatrixAssembly(1)%p_IdofsTrial, p_Kentryai,&
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(1)%p_IdofsTrial,1), IELmax-IELset+1)    
       !      call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
       !            rlocalMatrixAssembly(2)%p_IdofsTrial, p_Kentryaa,&
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
       !            ubound(rlocalMatrixAssembly(2)%p_IdofsTrial,1), IELmax-IELset+1)



       call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
            rlocalMatrixAssembly(1)%p_IdofsTest, p_Kentryii,&
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(1)%p_IdofsTest, &
            rlocalMatrixAssembly(2)%p_IdofsTest, p_Kentryai,&
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
            rlocalMatrixAssembly(1)%p_IdofsTest, p_Kentryia,&
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(1)%p_IdofsTest,1), IELmax-IELset+1)    
       call bilf_getLocalMatrixIndices (rmatrix,rlocalMatrixAssembly(2)%p_IdofsTest, &
            rlocalMatrixAssembly(2)%p_IdofsTest, p_Kentryaa,&
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), &
            ubound(rlocalMatrixAssembly(2)%p_IdofsTest,1), IELmax-IELset+1)   

       ! -------------------- ELEMENT EVALUATION PHASE ----------------------

       ! Ok, we found the positions of the local matrix entries
       ! that we have to change.
       ! To calculate the matrix contributions, we have to evaluate
       ! the elements to give us the values of the basis functions
       ! in all the DOF`s in all the elements in our set.

       ! Get the element evaluation tag of all FE spaces. We need it to evaluate
       ! the elements later. All of them can be combined with OR, what will give
       ! a combined evaluation tag. 
       cevaluationTag = rlocalMatrixAssembly(1)%cevaluationTag

       ! The cubature points are already initialised by 1D->2D mapping.
       cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

       !      ! Do we have a (multi-)linear transformation?
       !      bisLinearTrafo = trafo_isLinearTrafo(rlocalMatrixAssembly%ctrafoType)
       !
       !      if (bisLinearTrafo) then
       !        ! We need the vertices of the element corners and the number
       !        ! of vertices per element to compute the length of the element
       !        ! edge at the boundary
       !        cevaluationTag = ior(cevaluationTag, EL_EVLTAG_COORDS)
       !        nve = trafo_igetNVE(rlocalMatrixAssembly%ctrafoType)
       !      end if

       ! Calculate all information that is necessary to evaluate the finite element
       ! on all cells of our subset. This includes the coordinates of the points
       ! on the cells.
       !      call elprep_prepareSetForEvaluation (p_revalElementSet,&
       !          cevaluationTag, rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
       !          IelementList(IELset:IELmax), rlocalMatrixAssembly%ctrafoType, &
       !          DpointsRef=DpointsRef)
       !      p_Dcoords => p_revalElementSet%p_Dcoords

       call elprep_prepareSetForEvaluation (&
            rlocalMatrixAssembly(1)%revalElementSet,&
            cevaluationTag,  rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
            IelementList(1,IELset:IELmax), rlocalMatrixAssembly(1)%ctrafoType, &
            DpointsRef=DpointsRef(:,:,:,1))
       call elprep_prepareSetForEvaluation (&
            rlocalMatrixAssembly(2)%revalElementSet,&
            cevaluationTag,  rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
            IelementList(3,IELset:IELmax), rlocalMatrixAssembly(2)%ctrafoType, &
            DpointsRef=DpointsRef(:,:,:,2))




       !      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
       !      if (.not. rlocalMatrixAssembly%rform%ballCoeffConstant) then
       !        if (present(fcoeff_buildMatrixScBdr2D_sim)) then
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

       ! If the flux function needs other, than just the function values from the solution
       ! (for example the derivatives), we will give an evalElementSet to it
       ! This is filled here

       call domint_initIntegrationByEvalSet (rlocalMatrixAssembly(1)%revalElementSet, rintSubset(1))
       call domint_initIntegrationByEvalSet (rlocalMatrixAssembly(2)%revalElementSet, rintSubset(2))
       !rintSubset(1)%ielementDistribution = 0
       rintSubset(1)%ielementStartIdx = IELset
       rintSubset(1)%p_Ielements => IelementList(1,IELset:IELmax)
       rintSubset(1)%p_IdofsTrial => rlocalMatrixAssembly(1)%p_IdofsTest
       rintSubset(1)%celement = rlocalMatrixAssembly(1)%celementTest
       !rintSubset(2)%ielementDistribution = 0
       rintSubset(2)%ielementStartIdx = IELset
       rintSubset(2)%p_Ielements => IelementList(2,IELset:IELmax)
       rintSubset(2)%p_IdofsTrial => rlocalMatrixAssembly(2)%p_IdofsTest
       rintSubset(2)%celement = rlocalMatrixAssembly(2)%celementTest


       call flux_dg_buildMatrixScEdge2D_sim (&
            !            rlocalVectorAssembly(1)%p_Dcoefficients(1,:,1:IELmax-IELset+1),&
            !            DsolVals(:,:,1:IELmax-IELset+1),&
       DfluxValues(:,:,1:IELmax-IELset+1),&
            rvectorSol,&
            IelementList(2,IELset:IELmax),&
            p_Dside,&
            raddTriaData%p_Dnormals(:,Iedgelist(IELset:IELmax)),&
            !DpointsReal(1:ndim2d,1:ncubp,1:IELmax-IELset+1),&
       rintSubset,&
            rcollection )


       call domint_doneIntegration (rintSubset(1))
       call domint_doneIntegration (rintSubset(2))





       ! Calculate the values of the basis functions.
       !      call elem_generic_sim2 (rlocalMatrixAssembly%celementTest, &
       !          p_revalElementSet, rlocalMatrixAssembly%BderTest, &
       !          rlocalMatrixAssembly%p_DbasTest) 
       call elem_generic_sim2 (rlocalMatrixAssembly(1)%celementTest, &
            rlocalMatrixAssembly(1)%revalElementSet,&
            rlocalMatrixAssembly(1)%BderTest, &
            rlocalMatrixAssembly(1)%p_DbasTest)
       call elem_generic_sim2 (rlocalMatrixAssembly(2)%celementTest, &
            rlocalMatrixAssembly(2)%revalElementSet,&
            rlocalMatrixAssembly(2)%BderTest, &
            rlocalMatrixAssembly(2)%p_DbasTest)

       !      ! Omit the calculation of the trial function values if they
       !      ! are identical to the test function values.
       !      if (.not. rlocalMatrixAssembly%bidenticalTrialAndTest) then
       !        call elem_generic_sim2 (rlocalMatrixAssembly%celementTrial, &
       !            p_revalElementSet, rlocalMatrixAssembly%BderTrial, &
       !            rlocalMatrixAssembly%p_DbasTrial)
       !      end if

       !      ! Calculate the length of egdes on the boundary. Depending on
       !      ! whether the transformation is (multi-)linear or not we compute
       !      ! the edge length as the distance between the two corner
       !      ! vertices of the element located on the boundary or as the real
       !      ! length of the boundary segment of the element.
       !      !
       !      ! The length of the current edge serves as a "determinant" in
       !      ! the cubature, so we have to divide it by 2 as an edge on the
       !      ! unit interval [-1,1] has length 2.
       !      if (bisLinearTrafo) then
       !        do iel = 1,IELmax-IELset+1
       !          DedgeLength(iel) = 0.5_DP*sqrt(&
       !              (p_Dcoords(1,    IelementOrientation(IELset+iel-1),iel)-&
       !               p_Dcoords(1,mod(IelementOrientation(IELset+iel-1),nve)+1,iel))**2+&
       !              (p_Dcoords(2,    IelementOrientation(IELset+iel-1),iel)-&
       !               p_Dcoords(2,mod(IelementOrientation(IELset+iel-1),nve)+1,iel))**2)
       !        end do
       !      else
       !        do iel = 1,IELmax-IELset+1
       !          DedgeLength(iel) = 0.5_DP*(DedgePosition(2,IELset+iel-1)-&
       !                                     DedgePosition(1,IELset+iel-1))
       !        end do
       !      end if

       ! --------------------- DOF COMBINATION PHASE ------------------------

       ! Values of all basis functions calculated. Now we can start 
       ! to integrate!

       ! Clear the local matrices
       p_Dentryii(:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryai(:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryia(:,:,1:IELmax-IELset+1) = 0.0_DP
       p_Dentryaa(:,:,1:IELmax-IELset+1) = 0.0_DP

       !      p_Dentryii = 0.0_DP
       !      p_Dentryai = 0.0_DP
       !      p_Dentryia = 0.0_DP
       !      p_Dentryaa = 0.0_DP

       ! We have two different versions for the integration - one
       ! with constant coefficients and one with nonconstant coefficients.
       !
       ! Check the bilinear form which one to use:

       !      if (rlocalMatrixAssembly%rform%ballCoeffConstant) then
       !      
       !        ! Constant coefficients. The coefficients are to be found in
       !        ! the Dcoefficients variable of the form.
       !        !
       !        ! Loop over the elements in the current set.
       !
       !        do iel = 1,IELmax-IELset+1
       !          
       !          ! Get the length of the edge.
       !          dlen = DedgeLength(iel)
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
       !              ! Now loop through all possible combinations of DOF`s
       !              ! in the current cubature point. The outer loop
       !              ! loops through the "O"`s in the above picture,
       !              ! the test functions:
       !
       !              do idofe = 1,indofTest
       !              
       !                ! Get the value of the (test) basis function 
       !                ! phi_i (our "O") in the cubature point:
       !                db = p_DbasTest(idofe,ib,icubp,iel)
       !                
       !                ! Perform an inner loop through the other DOF`s
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

       ! Nonconstant coefficients. The coefficients are to be found in
       ! the Dcoefficients variable as computed above.
       !
       ! Loop over the elements.

       do iel = 1,IELmax-IELset+1

          ! Get the length of the edge.
          !dlen = DedgeLength(iel)
          dlen = 0.5_DP*raddTriaData%p_Dedgelength(Iedgelist(IELset+iel-1))

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

             ! calculate the current weighting factor in the cubature formula
             ! in that cubature point.
             !            domega = dlen * p_Domega(icubp)
             domega1 = dlen * rlocalMatrixAssembly(1)%p_Domega(icubp)
             domega2 = dlen * rlocalMatrixAssembly(2)%p_Domega(icubp)

             ! Loop over the additive factors in the bilinear form.
             do ialbet = 1,rlocalMatrixAssembly(1)%rform%itermcount

                ! Get from Idescriptors the type of the derivatives for the 
                ! test and trial functions. The summand we calculate
                ! here will be added to the matrix entry:
                !
                ! a_ij  =  int_... ( psi_j )_ia  *  ( phi_i )_ib
                !
                ! -> Ix=0: function value, 
                !      =1: first derivative, ...
                !    as defined in the module 'derivative'.

                ia = rlocalMatrixAssembly(1)%rform%Idescriptors(1,ialbet)
                ib = rlocalMatrixAssembly(1)%rform%Idescriptors(2,ialbet)

                ! Multiply domega with the coefficient of the form.
                ! This gives the actual value to multiply the
                ! function value with before summing up to the integral.
                ! Get the precalculated coefficient from the coefficient array.
                daux1 = domega1 * DfluxValues(ialbet,icubp,iel)
                daux2 = domega2 * DfluxValues(ialbet,icubp,iel) * (-1.0_dp)

                ! Now loop through all possible combinations of DOF`s
                ! in the current cubature point. The outer loop
                ! loops through the "O" in the above picture,
                ! the test functions:

                do idofe = 1,indofTest

                   ! Get the value of the (test) basis function 
                   ! phi_i (our "O") in the cubature point:
                   db1 = rlocalMatrixAssembly(1)%p_DbasTest(idofe,ib,icubp,iel)
                   db2 = rlocalMatrixAssembly(2)%p_DbasTest(idofe,ib,icubp,iel)

                   ! Perform an inner loop through the other DOF`s
                   ! (the "X"). 

                   do jdofe = 1,indofTrial

                      ! Get the value of the basis function 
                      ! psi_j (our "X") in the cubature point. 
                      ! Them multiply:
                      !    db * dbas(..) * daux
                      ! ~= phi_i * psi_j * coefficient * cub.weight
                      ! Summing this up gives the integral, so the contribution
                      ! to the global matrix. 
                      !
                      ! Simply summing up db * dbas(..) * daux would give
                      ! the coefficient of the local matrix. We save this
                      ! contribution in the local matrix of element iel.

                      !JCOLB = Kentry(jdofe,idofe,iel)
                      !p_DA(JCOLB) = p_DA(JCOLB) + db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
                      !                  p_Dentry(jdofe,idofe,iel) = &
                      !                      p_Dentry(jdofe,idofe,iel)+db*p_DbasTrial(jdofe,ia,icubp,iel)*daux

                      !                  ! Testfunction on the 'first' (i) side
                      !                  p_Dentryii(jdofe,idofe,iel) = &
                      !                      p_Dentryii(jdofe,idofe,iel)+db1*rlocalMatrixAssembly(1)%p_DbasTrial(jdofe,ia,icubp,iel)*daux1*p_Dside(1,icubp,iel)   
                      !                  p_Dentryai(jdofe,idofe,iel) = &
                      !                      p_Dentryai(jdofe,idofe,iel)+db1*rlocalMatrixAssembly(2)%p_DbasTrial(jdofe,ia,icubp,iel)*daux1*p_Dside(2,icubp,iel)   
                      !                  
                      !                  ! Testfunction on the 'second' (a) side
                      !                  p_Dentryia(jdofe,idofe,iel) = &
                      !                      p_Dentryia(jdofe,idofe,iel)+db2*rlocalMatrixAssembly(1)%p_DbasTrial(jdofe,ia,icubp,iel)*daux2*p_Dside(1,icubp,iel)
                      !                  p_Dentryaa(jdofe,idofe,iel) = &
                      !                      p_Dentryaa(jdofe,idofe,iel)+db2*rlocalMatrixAssembly(2)%p_DbasTrial(jdofe,ia,icubp,iel)*daux2*p_Dside(2,icubp,iel)
                      !                


                      ! Testfunction on the 'first' (i) side
                      p_Dentryii(jdofe,idofe,iel) = &
                           p_Dentryii(jdofe,idofe,iel)+db1*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux1*p_Dside(1,ialbet)
                      p_Dentryai(jdofe,idofe,iel) = &
                           p_Dentryai(jdofe,idofe,iel)+db1*rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)*daux1*p_Dside(2,ialbet)

                      ! Testfunction on the 'second' (a) side
                      p_Dentryia(jdofe,idofe,iel) = &
                           p_Dentryia(jdofe,idofe,iel)+db2*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux2*p_Dside(1,ialbet)

                      !                      if ((p_Dentryia(jdofe,idofe,iel)<-1000000000.0_dp).and.(IelementList(2,IELset+iel-1).ne.0)) then
                      !                write(*,*) 'Added', db2*rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)*daux2*p_Dside(1,iel)      
                      !                write(*,*) 'ia',ia
                      !                write(*,*) 'daux1',daux1
                      !                write(*,*) 'daux2',daux2
                      !                write(*,*) 'db1',db1
                      !                write(*,*) 'db2',db2
                      !                write(*,*) 'dside1',p_Dside(1,iel)
                      !                write(*,*) 'dside2',p_Dside(2,iel)
                      !                write(*,*) 'test1',rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)
                      !                write(*,*) 'test2',rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)
                      !                        pause
                      !                      end if

                      p_Dentryaa(jdofe,idofe,iel) = &
                           p_Dentryaa(jdofe,idofe,iel)+db2*rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)*daux2*p_Dside(2,ialbet)

                      !                write(*,*) 'ia',ia
                      !                write(*,*) 'daux1',daux1
                      !                write(*,*) 'daux2',daux2
                      !                write(*,*) 'db1',db1
                      !                write(*,*) 'db2',db2
                      !                write(*,*) 'dside1',p_Dside(1,iel)
                      !                write(*,*) 'dside2',p_Dside(2,iel)
                      !                write(*,*) 'test1',rlocalMatrixAssembly(1)%p_DbasTest(jdofe,ia,icubp,iel)
                      !                write(*,*) 'test2',rlocalMatrixAssembly(2)%p_DbasTest(jdofe,ia,icubp,iel)
                      !                pause


                   end do

                end do ! jdofe

             end do ! ialbet

          end do ! icubp 

       end do ! iel

       !      end if ! rform%ballCoeffConstant

       ! Incorporate the local matrices into the global one.
       ! Kentry gives the position of the additive contributions in Dentry.
       !
       ! OpenMP-Extension: This is a critical section. Only one thread is
       ! allowed to write to the matrix, otherwise the matrix may get
       ! messed up.
       ! The critical section is put around both loops as indofTest/indofTrial
       ! are usually small and quickly to handle.

       !      if (cconstrType .eq. BILF_MATC_LUMPED) then
       !
       !        !$omp critical
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
       !        !$omp end critical
       !
       !      else

       !$omp critical
       do iel = 1,IELmax-IELset+1

          do idofe = 1,indofTest
             do jdofe = 1,indofTrial
                !              p_DA(p_Kentry(jdofe,idofe,iel)) = &
                !                  p_DA(p_Kentry(jdofe,idofe,iel)) + p_Dentry(jdofe,idofe,iel)

                p_DA(p_Kentryii(jdofe,idofe,iel)) = &
                     p_DA(p_Kentryii(jdofe,idofe,iel)) + p_Dentryii(jdofe,idofe,iel)



                if (IelementList(2,IELset+iel-1).ne.0) then

                   p_DA(p_Kentryia(jdofe,idofe,iel)) = &
                        p_DA(p_Kentryia(jdofe,idofe,iel)) + p_Dentryia(jdofe,idofe,iel)

                   p_DA(p_Kentryai(jdofe,idofe,iel)) = &
                        p_DA(p_Kentryai(jdofe,idofe,iel)) + p_Dentryai(jdofe,idofe,iel)!*real(min(1,IelementList(2,IELset+iel-1)))
                   p_DA(p_Kentryaa(jdofe,idofe,iel)) = &
                        p_DA(p_Kentryaa(jdofe,idofe,iel)) + p_Dentryaa(jdofe,idofe,iel)!*real(min(1,IelementList(2,IELset+iel-1)))
                end if

             end do
          end do

       end do ! iel
       !$omp end critical

       !      end if

    end do ! IELset
    !$omp end do

    ! Release the local matrix assembly structure
    call bilf_releaseAssemblyData(rlocalMatrixAssembly(1))
    call bilf_releaseAssemblyData(rlocalMatrixAssembly(2))

    ! Deallocate memory
    deallocate(Dxi2D, DpointsRef) !, DpointsPar, DedgeLength)
    deallocate(p_Kentryii,p_Kentryia,p_Kentryai,p_Kentryaa)
    deallocate(p_Dentryii,p_Dentryia,p_Dentryai,p_Dentryaa)
    deallocate(p_Dside)
    deallocate(IelementList)

    !$omp end parallel

  end subroutine dg_bilf_assembleSubmeshMat9Bdr2D_ss
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
   !****************************************************************************

  !<subroutine>  

  subroutine dg_pM_project (rvectorSource,rvectorDest)

    !<description>

    ! For p-Multigrid. Takes a dg_Tx vector and projects it to a dg_Ty vector.

    !</description>

    !<input>

    ! The lower order input vector
    type(t_vectorScalar), intent(in) :: rvectorSource

    !</input>

    !<output>

    ! The higher order output vector
    type(t_vectorScalar), intent(inout) :: rvectorDest

    !</output>

    !</subroutine>
 
    integer :: ndofloc1, ndofloc2, iel, NEL, nDOFlocSource, nDOFlocDest
    integer(I32) :: celementSource, celementDest
    real(DP), dimension(:), pointer :: p_DdataSource, p_DdataDest
    integer, dimension(6) :: IdofGlobSource, IdofGlobDest

    ! Get number of elements
    NEL = rvectorSource%p_rspatialDiscr%p_rtriangulation%NEL

    ! What is the current element type?
    celementSource = rvectorSource%p_rspatialDiscr%RelementDistr(1)%celement
    celementDest = rvectorDest%p_rspatialDiscr%RelementDistr(1)%celement

    ! Get number of local DOF
    nDOFlocSource = elem_igetNDofLoc(celementSource)
    nDOFlocDest = elem_igetNDofLoc(celementDest)
    
    ! Get Pointers to data
    call lsyssc_getbase_double(rvectorSource, p_DdataSource)
    call lsyssc_getbase_double(rvectorDest, p_DdataDest)
    
    ! Initialize destination vector DOFs
    p_DdataDest(:) = 0.0_dp
    
    ! Loop over all elements
    do iel = 1, NEL
      
      ! Get global DOFs from the elements
      call dof_locGlobMapping(rvectorSource%p_rspatialDiscr, iel, IdofGlobSource(1:nDOFlocSource))
      call dof_locGlobMapping(rvectorDest%p_rspatialDiscr, iel, IdofGlobDest(1:nDOFlocDest))
      
      ! Copy DOFs from source to dest vector
      p_DdataDest(IdofGlobDest(1:min(nDOFlocSource,nDOFlocDest))) = p_DdataSource(IdofGlobSource(1:min(nDOFlocSource,nDOFlocDest)))
      
      
    
    end do 
  
  
  
  
  end subroutine
  
  
!  !****************************************************************************
!
!!<subroutine>
!
!  subroutine dg_pperr_scalar2d_conf (cerrortype, derror, rdiscretisation,&
!                                  rvectorScalar, ffunctionReference,&
!                                  rcollection, relementError, ffunctionWeight)
!
!!<description>
!  ! This routine calculates the error of a given finite element function
!  ! in rvector to a given analytical callback function ffunctionReference.
!  ! 2D version for double-precision vectors.
!!</description>
!
!!<input> 
!  ! Type of error to compute. Bitfield. This is a combination of the
!  ! PPERR_xxxx-constants, which specifies what to compute.
!  ! Example: PPERR_L2ERROR computes the $L_2$-error.
!  integer, intent(in) :: cerrortype
!  
!  ! A discretisation structure specifying how to compute the error.
!  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
!
!  ! OPTIONAL: The FE solution vector. Represents a scalar FE function.
!  ! If omitted, the function is assumed to be constantly =0.
!  type(t_vectorScalar), intent(in), target, optional :: rvectorScalar
!
!  ! OPTIONAL: A callback function that provides the analytical reference 
!  ! function to which the error should be computed.
!  ! If not specified, the reference function is assumed to be zero!
!  include 'intf_refFunctionSc.inc'
!  optional :: ffunctionReference
!
!  ! OPTIONAL: A callback function that provides the weighting function
!  ! by which the computed error is multipled.
!  ! If not specified, the reference function is assumed to be =1!
!  optional :: ffunctionWeight
!!</input>
!
!!<inputoutput>
!  ! OPTIONAL: A collection structure to provide additional 
!  ! information for callback routines.
!  type(t_collection), intent(inout), optional :: rcollection
!
!  ! OPTIONAL: A scalar vector which holds the calculated error per element
!  type(t_vectorScalar), intent(inout), optional :: relementError
!!</inputoutput>
!
!!<output>
!  ! Array receiving the calculated error.
!  real(DP), intent(out) :: derror
!!</output>
!
!!</subroutine>
!
!    ! local variables
!    integer :: ielementDistr, ICUBP, NVE, NCOEFF
!    integer :: IEL, IELmax, IELset, IELGlobal
!    real(DP) :: OM
!    
!    ! Array to tell the element which derivatives to calculate
!    logical, dimension(EL_MAXNDER) :: Bder
!    
!    ! For every cubature point on the reference element,
!    ! the corresponding cubature weight
!    real(DP), dimension(:), allocatable :: Domega
!    
!    ! number of cubature points on the reference element
!    integer :: ncubp
!    
!    ! Number of local degees of freedom for test functions
!    integer :: indofTrial
!    
!    ! The triangulation structure - to shorten some things...
!    type(t_triangulation), pointer :: p_rtriangulation
!    
!    ! A pointer to an element-number list
!    integer, dimension(:), pointer :: p_IelementList
!    
!    ! An array receiving the coordinates of cubature points on
!    ! the reference element for all elements in a set.
!    real(DP), dimension(:,:), allocatable :: p_DcubPtsRef
!
!    ! Arrays for saving Jacobian determinants and matrices
!    real(DP), dimension(:,:), pointer :: p_Ddetj
!    
!    ! Current element distribution
!    type(t_elementDistribution), pointer :: p_relementDistribution
!    
!    ! Number of elements in the current element distribution
!    integer :: NEL, ielement1, ipoint1
!    
!    real(dp) :: dx
!
!    ! Pointer to the values of the function that are computed by the callback routine.
!    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
!    
!    ! Number of elements in a block. Normally =PPERR_NELEMSIM,
!    ! except if there are less elements in the discretisation.
!    integer :: nelementsPerBlock
!    
!    ! A t_domainIntSubset structure that is used for storing information
!    ! and passing it to callback routines.
!    type(t_domainIntSubset) :: rintSubset
!    type(t_evalElementSet) :: revalElementSet
!    
!    ! An allocateable array accepting the DOF`s of a set of elements.
!    integer, dimension(:,:), allocatable, target :: IdofsTrial
!  
!    ! Type of transformation from the reference to the real element 
!    integer(I32) :: ctrafoType
!    
!    ! Element evaluation tag; collects some information necessary for evaluating
!    ! the elements.
!    integer(I32) :: cevaluationTag
!
!    ! Pointer to the element-wise error
!    real(DP), dimension(:), pointer :: p_Derror
!
!
!    ! Which derivatives of basis functions are needed?
!    ! Check the descriptors of the bilinear form and set BDER
!    ! according to these.
!
!    Bder = .false.
!    select case (cerrortype)
!    case (PPERR_L1ERROR, PPERR_L2ERROR, PPERR_MEANERROR)
!      Bder(DER_FUNC) = .true.
!      NCOEFF = 3
!    case (PPERR_H1ERROR) 
!      Bder(DER_DERIV_X) = .true.
!      Bder(DER_DERIV_Y) = .true.
!      NCOEFF = 5
!    case default
!      NCOEFF = 2
!    end select
!        
!    ! Get a pointer to the triangulation - for easier access.
!    p_rtriangulation => rdiscretisation%p_rtriangulation
!    
!    ! For saving some memory in smaller discretisations, we calculate
!    ! the number of elements per block. For smaller triangulations,
!    ! this is NEL. If there are too many elements, it is at most
!    ! PPERR_NELEMSIM. This is only used for allocating some arrays.
!    nelementsPerBlock = min(PPERR_NELEMSIM, p_rtriangulation%NEL)
!    
!    ! Set the current error to 0 and add the error contributions of each element
!    ! to that.
!    derror = 0.0_DP
!
!    ! Set pointer to element-wise error
!    if (present(relementError)) then
!      call lsyssc_getbase_double(relementError, p_Derror)
!    end if
!
!    ! Now loop over the different element distributions (=combinations
!    ! of trial and test functions) in the discretisation.
!
!    do ielementDistr = 1, rdiscretisation%inumFESpaces
!    
!      ! Activate the current element distribution
!      p_relementDistribution => rdiscretisation%RelementDistr(ielementDistr)
!    
!      ! Cancel if this element distribution is empty.
!      if (p_relementDistribution%NEL .eq. 0) cycle
!
!      ! Get the number of local DOF`s for trial functions
!      indofTrial = elem_igetNDofLoc(p_relementDistribution%celement)
!      
!      ! Get the number of corner vertices of the element
!      NVE = elem_igetNVE(p_relementDistribution%celement)
!      
!      ! Get from the trial element space the type of coordinate system
!      ! that is used there:
!      ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
!
!      ! Get the number of cubature points for the cubature formula
!      ncubp = cub_igetNumPts(p_relementDistribution%ccubTypeEval)
!      
!      ! Allocate two arrays for the points and the weights
!      allocate(Domega(ncubp))
!      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType), ncubp))
!      
!      ! Get the cubature formula
!      call cub_getCubature(p_relementDistribution%ccubTypeEval, p_DcubPtsRef, Domega)
!      
!      ! Allocate memory for the DOF`s of all the elements.
!      allocate(IdofsTrial(indofTrial, nelementsPerBlock))
!
!      ! Allocate memory for the coefficients
!      allocate(Dcoefficients(ncubp, nelementsPerBlock, NCOEFF))
!    
!      ! Initialisation of the element set.
!      call elprep_init(revalElementSet)
!
!      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
!      ! the elements later. All of them can be combined with OR, what will give
!      ! a combined evaluation tag. 
!      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
!                      
!      if (present(ffunctionReference) .or.&
!          present(ffunctionWeight)) then
!        ! Evaluate real coordinates if necessary.
!        cevaluationTag = ior(cevaluationTag, EL_EVLTAG_REALPOINTS)
!      end if
!                      
!      ! Make sure that we have determinants.
!      cevaluationTag = ior(cevaluationTag, EL_EVLTAG_DETJ)
!
!      ! p_IelementList must point to our set of elements in the discretisation
!      ! with that combination of trial functions
!      call storage_getbase_int (p_relementDistribution%h_IelementList, &
!                                p_IelementList)
!                     
!      ! Get the number of elements there.
!      NEL = p_relementDistribution%NEL
!    
!      ! Loop over the elements - blockwise.
!      do IELset = 1, NEL, PPERR_NELEMSIM
!      
!        ! We always handle LINF_NELEMSIM elements simultaneously.
!        ! How many elements have we actually here?
!        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
!        ! elements simultaneously.
!        
!        IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
!      
!        ! Calculate the global DOF`s into IdofsTrial.
!        !
!        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
!        ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
!        call dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
!                                     IdofsTrial)
!                                     
!        ! Prepare the call to the evaluation routine of the analytic function.    
!        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
!        rintSubset%ielementDistribution = ielementDistr
!        rintSubset%ielementStartIdx = IELset
!        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
!        rintSubset%p_IdofsTrial => IdofsTrial
!        rintSubset%celement = p_relementDistribution%celement
!    
!        ! Calculate all information that is necessary to evaluate the finite element
!        ! on all cells of our subset. This includes the coordinates of the points
!        ! on the cells.
!        call elprep_prepareSetForEvaluation (revalElementSet,&
!            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
!            ctrafoType, p_DcubPtsRef(:,1:ncubp))
!        p_Ddetj => revalElementSet%p_Ddetj
!
!        ! In the next loop, we do not have to evaluate the coordinates
!        ! on the reference elements anymore.
!        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
!
!        ! At this point, we must select the correct domain integration and coefficient
!        ! calculation routine, depending which type of error we should compute!
!        
!        select case (cerrortype)
!        
!        case (PPERR_L2ERROR)
!          
!          ! L2-error uses only the values of the function.
!          
!          if (present(ffunctionReference)) then
!            ! It is time to call our coefficient function to calculate the
!            ! function values in the cubature points:  u(x,y)
!            ! The result is saved in Dcoefficients(:,:,1)
!            call ffunctionReference (DER_FUNC, rdiscretisation, &
!                        int(IELmax-IELset+1), ncubp,&
!                        revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!                        IdofsTrial, rintSubset,&
!                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP
!          end if
!          
!          if (present(rvectorScalar)) then
!            ! Calculate the values of the FE function in the
!            ! cubature points: u_h(x,y).
!            ! Save the result to Dcoefficients(:,:,2)
!            call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
!                    p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
!                    Dcoefficients(:,1:IELmax-IELset+1,2))
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,2) = 0.0_DP
!          end if
!
!          if (present(ffunctionWeight)) then
!            ! Calculate the values of the weighting function in
!            ! the cubature points: w(x).
!            ! Save the result to Dcoefficients(:,:,3)
!!            call ffunctionWeight (rdiscretisation,&
!!                int(IELmax-IELset+1), ncubp, &
!!                revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!!                IdofsTrial, rintSubset, &
!!                Dcoefficients(:,1:IELmax-IELset+1,3), rcollection)
!          else
!            
!            do ielement1 = 1,IELmax-IELset+1
!            do ipoint1 = 1,size(revalElementSet%p_DpointsReal,2)
!            dx = revalElementSet%p_DpointsReal(1,1,ielement1)
!            if ((dx <0.45_dp).and.(dx>0.25)) then
!               Dcoefficients(ipoint1,ielement1,3) = 1.0_DP
!            else
!               Dcoefficients(ipoint1,ielement1,3) = 0.0_DP
!            end if
!            end do
!            end do
!            
!          end if
!
!          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2)
!          ! and multiplication by Dcoefficients(:,:,3) yields
!          ! the error "w*[u-u_h](cubature pt.)"!
!          !        
!          ! Loop through elements in the set and for each element,
!          ! loop through the DOF`s and cubature points to calculate the
!          ! integral: int_Omega w*(u-u_h,u-u_h) dx
!          
!          if (present(relementError)) then
!
!            do IEL=1,IELmax-IELset+1
!          
!              ! Loop over all cubature points on the current element
!              do icubp = 1, ncubp
!                
!                ! calculate the current weighting factor in the cubature formula
!                ! in that cubature point.
!                !
!                ! Take the absolut value of the determinant of the mapping.
!                ! In 2D, the determinant is always positive, whereas in 3D,
!                ! the determinant might be negative -- that is normal!
!                
!                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))*Dcoefficients(icubp,IEL,3)
!                
!                ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
!
!                IELGlobal = p_IelementList(IELset+IEL-1)
!                
!                p_Derror(IELGlobal) = OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2
!
!                derror = derror + p_Derror(IELGlobal)
!                
!              end do ! ICUBP 
!              
!            end do ! IEL
!
!          else
!
!            do IEL=1,IELmax-IELset+1
!          
!              ! Loop over all cubature points on the current element
!              do icubp = 1, ncubp
!                
!                ! calculate the current weighting factor in the cubature formula
!                ! in that cubature point.
!                !
!                ! Take the absolut value of the determinant of the mapping.
!                ! In 2D, the determinant is always positive, whereas in 3D,
!                ! the determinant might be negative -- that is normal!
!                
!                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))*Dcoefficients(icubp,IEL,3)
!                
!                ! L2-error is:   int_... (u-u_h)*(u-u_h) dx
!                
!                derror = derror + &
!                         OM * (Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))**2
!                
!              end do ! ICUBP 
!              
!            end do ! IEL
!
!          end if
!
!        case (PPERR_L1ERROR)
!          
!          ! L1-error uses only the values of the function.
!          
!          if (present(ffunctionReference)) then
!            ! It is time to call our coefficient function to calculate the
!            ! function values in the cubature points:  u(x,y)
!            ! The result is saved in Dcoefficients(:,:,1)
!            call ffunctionReference (DER_FUNC, rdiscretisation, &
!                        int(IELmax-IELset+1), ncubp,&
!                        revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!                        IdofsTrial, rintSubset,&
!                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP
!          end if
!          
!          if (present(rvectorScalar)) then
!            ! Calculate the values of the FE function in the
!            ! cubature points: u_h(x,y).
!            ! Save the result to Dcoefficients(:,:,2)
!            call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
!                    p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
!                    Dcoefficients(:,1:IELmax-IELset+1,2))
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,2) = 0.0_DP
!          end if
!          
!          if (present(ffunctionWeight)) then
!            ! Calculate the values of the weighting function in
!            ! the cubature points: w(x).
!            ! Save the result to Dcoefficients(:,:,3)
!!            call ffunctionWeight (rdiscretisation,&
!!                int(IELmax-IELset+1), ncubp, &
!!                revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!!                IdofsTrial, rintSubset, &
!!                Dcoefficients(:,1:IELmax-IELset+1,3), rcollection)
!          else
!            do ielement1 = 1,IELmax-IELset+1
!            do ipoint1 = 1,size(revalElementSet%p_DpointsReal,2)
!            dx = revalElementSet%p_DpointsReal(1,1,ielement1)
!            if ((dx <0.45_dp).and.(dx>0.25)) then
!               Dcoefficients(ipoint1,ielement1,3) = 1.0_DP
!            else
!               Dcoefficients(ipoint1,ielement1,3) = 0.0_DP
!            end if
!            end do
!            end do
!          end if
!
!          ! Subtraction of Dcoefficients(:,:,1) from Dcoefficients(:,:,2)
!          ! and multiplication by Dcoefficients(:,:,3) yields
!          ! the error "w*[u-u_h](cubature pt.)"!
!          !        
!          ! Loop through elements in the set and for each element,
!          ! loop through the DOF`s and cubature points to calculate the
!          ! integral: int_Omega w*abs(u-u_h) dx
!          
!          if (present(relementError)) then
!
!            do IEL=1,IELmax-IELset+1
!              
!              ! Loop over all cubature points on the current element
!              do icubp = 1, ncubp
!                
!                ! calculate the current weighting factor in the cubature formula
!                ! in that cubature point.
!                
!                OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)*Dcoefficients(icubp,IEL,3)
!                
!                ! L1-error is:   int_... abs(u-u_h) dx
!                
!                IELGlobal = p_IelementList(IELset+IEL-1)
!
!                p_Derror(IELGlobal) = OM * abs(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))
!
!                derror = derror + p_Derror(IELGlobal)
!                
!              end do ! ICUBP 
!              
!            end do ! IEL
!
!          else
!            
!            do IEL=1,IELmax-IELset+1
!              
!              ! Loop over all cubature points on the current element
!              do icubp = 1, ncubp
!                
!                ! calculate the current weighting factor in the cubature formula
!                ! in that cubature point.
!                
!                OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)*Dcoefficients(icubp,IEL,3)
!                
!                ! L1-error is:   int_... abs(u-u_h) dx
!                
!                derror = derror + &
!                         OM * abs(Dcoefficients(icubp,IEL,2)-Dcoefficients(icubp,IEL,1))
!                
!              end do ! ICUBP 
!              
!            end do ! IEL
!
!          end if
!
!        case (PPERR_H1ERROR)
!
!          ! H1-error uses only 1st derivative of the function.
!
!          if (present(ffunctionReference)) then          
!            ! It is time to call our coefficient function to calculate the
!            ! X-derivative values in the cubature points:  u(x,y)
!            ! The result is saved in Dcoefficients(:,:,1)
!            call ffunctionReference (DER_DERIV_X, rdiscretisation, &
!                        int(IELmax-IELset+1), ncubp,&
!                        revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!                        IdofsTrial, rintSubset,&
!                        Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
!                        
!            ! Calculate the Y-derivative to Dcoefficients(:,:,2)
!            call ffunctionReference (DER_DERIV_Y,rdiscretisation, &
!                        int(IELmax-IELset+1), ncubp,&
!                        revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!                        IdofsTrial, rintSubset,&
!                        Dcoefficients(:,1:IELmax-IELset+1,2), rcollection)
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,1:2) = 0.0_DP
!          end if
!          
!          if (present(rvectorScalar)) then
!            ! Calculate the X/Y-derivative of the FE function in the
!            ! cubature points: u_h(x,y).
!            ! Save the result to Dcoefficients(:,:,3..4)
!            call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
!                    p_relementDistribution%celement, IdofsTrial, DER_DERIV_X,&
!                    Dcoefficients(:,1:IELmax-IELset+1,3))
!
!            call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
!                    p_relementDistribution%celement, IdofsTrial, DER_DERIV_Y,&
!                    Dcoefficients(:,1:IELmax-IELset+1,4))
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,3:4) = 0.0_DP
!          end if
!
!          if (present(ffunctionWeight)) then
!            ! Calculate the values of the weighting function in
!            ! the cubature points: w(x).
!            ! Save the result to Dcoefficients(:,:,5)
!            call ffunctionWeight (rdiscretisation,&
!                int(IELmax-IELset+1), ncubp, &
!                revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!                IdofsTrial, rintSubset, &
!                Dcoefficients(:,1:IELmax-IELset+1,5), rcollection)
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,5) = 1.0_DP
!          end if
!
!          ! Subtraction of Dcoefficients(:,:,1..2) from Dcoefficients(:,:,3..4) 
!          ! and multiplication by Dcoefficients(:,:,5) yields
!          ! the error "w*grad(u-u_h)(cubature pt.)"!
!          !        
!          ! Loop through elements in the set and for each element,
!          ! loop through the DOF`s and cubature points to calculate the
!          ! integral: int_Omega w*(grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
!          
!          if (present(relementError)) then
!
!            do IEL=1,IELmax-IELset+1
!              
!              ! Loop over all cubature points on the current element
!              do icubp = 1, ncubp
!                
!                ! calculate the current weighting factor in the cubature formula
!                ! in that cubature point.
!                !
!                ! Take the absolut value of the determinant of the mapping.
!                ! In 2D, the determinant is always positive, whereas in 3D,
!                ! the determinant might be negative -- that is normal!
!                
!                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))*Dcoefficients(icubp,IEL,5)
!                
!                ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
!                
!                IELGlobal = p_IelementList(IELset+IEL-1)
!
!                p_Derror(IELGlobal) = OM * ((Dcoefficients(icubp,IEL,3)-Dcoefficients(icubp,IEL,1))**2 + &
!                                            (Dcoefficients(icubp,IEL,4)-Dcoefficients(icubp,IEL,2))**2)
!
!                derror = derror + p_Derror(IELGlobal)
!
!              end do ! ICUBP 
!              
!            end do ! IEL
!
!          else
!            
!            do IEL=1,IELmax-IELset+1
!              
!              ! Loop over all cubature points on the current element
!              do icubp = 1, ncubp
!                
!                ! calculate the current weighting factor in the cubature formula
!                ! in that cubature point.
!                !
!                ! Take the absolut value of the determinant of the mapping.
!                ! In 2D, the determinant is always positive, whereas in 3D,
!                ! the determinant might be negative -- that is normal!
!                
!                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))*Dcoefficients(icubp,IEL,5)
!                
!                ! H1-error is:   int_... (grad(u)-grad(u_h),grad(u)-grad(u_h)) dx
!                
!                derror = derror + OM * &
!                         ((Dcoefficients(icubp,IEL,3)-Dcoefficients(icubp,IEL,1))**2 + &
!                          (Dcoefficients(icubp,IEL,4)-Dcoefficients(icubp,IEL,2))**2)
!
!              end do ! ICUBP 
!              
!            end do ! IEL
!
!          end if
!
!        case (PPERR_MEANERROR)
!
!          ! The integral mean value uses only the values of the function
!          
!          if (present(ffunctionReference)) then
!            ! Calculate the values of the coefficient function in the
!            ! cubature points: u(x).
!            ! Save the result to Dcoefficients(:,:,1)
!            call ffunctionReference(DER_FUNC, rdiscretisation,&
!                int(IELmax-IELset+1), ncubp, &
!                revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!                IdofsTrial, rintSubset, &
!                Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,1) = 0.0_DP
!          end if
!          
!          if (present(rvectorScalar)) then
!            ! Calculate the values of the FE function in the
!            ! cubature points: u_h(x).
!            ! Save the result to Dcoefficients(:,:,2)
!            call fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
!                p_relementDistribution%celement, IdofsTrial, DER_FUNC,&
!                Dcoefficients(:,1:IELmax-IELset+1,2))
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,2) = 0.0_DP
!          end if
!          
!          if (present(ffunctionWeight)) then
!            ! Calculate the values of the weighting function in
!            ! the cubature points: w(x).
!            ! Save the result to Dcoefficients(:,:,3)
!            call ffunctionWeight (rdiscretisation,&
!                int(IELmax-IELset+1), ncubp, &
!                revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!                IdofsTrial, rintSubset, &
!                Dcoefficients(:,1:IELmax-IELset+1,3), rcollection)
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,3) = 1.0_DP
!          end if
!
!          ! Subtraction of Dcoefficients(:,:,2) from Dcoefficients(:,:,1)
!          ! and multiplication by Dcoefficients(:,:,3) yields
!          ! the error "w*[u-u_h] (cubature pt.)"!
!          !        
!          ! Loop through elements in the set and for each element,
!          ! loop through the DOF`s and cubature points to calculate the
!          ! integral: int_Omega w*(u-u_h) dx
!
!          if (present(relementError)) then
!            
!            do IEL = 1, IELmax-IELset+1
!              
!              ! Loop over all cubature points on the current element
!              do icubp = 1, ncubp
!                
!                ! calculate the current weighting factor in the cubature formula
!                ! in that cubature point.
!                !
!                ! Take the absolut value of the determinant of the mapping.
!                ! In 2D, the determinant is always positive, whereas in 3D,
!                ! the determinant might be negative -- that is normal!
!                
!                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))*Dcoefficients(icubp,IEL,3)
!                
!                IELGlobal = p_IelementList(IELset+IEL-1)
!                
!                p_Derror(IELGlobal) = OM * ( Dcoefficients(icubp,IEL,1) - &
!                                             Dcoefficients(icubp,IEL,2) )
!
!                derror = derror + p_Derror(IELGlobal)
!                
!              end do ! ICUBP 
!              
!            end do ! IEL
!            
!          else
!            
!            do IEL = 1, IELmax-IELset+1
!              
!              ! Loop over all cubature points on the current element
!              do icubp = 1, ncubp
!                
!                ! calculate the current weighting factor in the cubature formula
!                ! in that cubature point.
!                !
!                ! Take the absolut value of the determinant of the mapping.
!                ! In 2D, the determinant is always positive, whereas in 3D,
!                ! the determinant might be negative -- that is normal!
!                
!                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))*Dcoefficients(icubp,IEL,3)
!                
!                derror = derror + &
!                         OM * ( Dcoefficients(icubp,IEL,1) - &
!                                Dcoefficients(icubp,IEL,2) )
!                
!              end do ! ICUBP 
!              
!            end do ! IEL
!
!          end if
!        
!        case default
!
!          ! This case realises
!          !
!          !  int (w u) dx
!          !  
!          ! with w being the weight and u the analytical function given by
!          ! ffunctionReference. This function must be present such that it can
!          ! be evaluated -- otherwise the user made a mistake in calling
!          ! this routine.
!
!          ! Evaluate the error with aid of the callback function
!          if (present(ffunctionReference)) then
!            ! Calculate the values of the coefficient function in the
!            ! cubature points: u(x).
!            ! Save the result to Dcoefficients(:,:,1)
!            call ffunctionReference(DER_FUNC, rdiscretisation,&
!                int(IELmax-IELset+1), ncubp, &
!                revalElementSet%p_DpointsReal,&
!                IdofsTrial, rintSubset, &
!                Dcoefficients(:,1:IELmax-IELset+1,1), rcollection)
!          else
!            call output_line('Reference function missing in user-defined error type!',&
!                OU_CLASS_ERROR,OU_MODE_STD,'pperr_scalar2d_conf')
!            call sys_halt()
!          end if
!
!          if (present(ffunctionWeight)) then
!            ! Calculate the values of the weighting function in
!            ! the cubature points: w(x).
!            ! Save the result to Dcoefficients(:,:,2)
!            call ffunctionWeight (rdiscretisation,&
!                int(IELmax-IELset+1), ncubp, &
!                revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
!                IdofsTrial, rintSubset, &
!                Dcoefficients(:,1:IELmax-IELset+1,2), rcollection)
!          else
!            Dcoefficients(:,1:IELmax-IELset+1,2) = 1.0_DP
!          end if
!
!          ! Multiplication of Dcoefficients(:,:,1) by Dcoefficients(:,:,2) yields
!          ! the error "w*f(u)(cubature pt.)"!
!          !
!          ! Loop through elements in the set and for each element,
!          ! loop through the DOF`s and cubature points to calculate the
!          ! integral: int_Omega w*f(u) dx
!
!          if (present(relementError)) then
!            
!            do IEL = 1, IELmax-IELset+1
!              
!              ! Loop over all cubature points on the current element
!              do icubp = 1, ncubp
!                
!                ! calculate the current weighting factor in the cubature formula
!                ! in that cubature point.
!                !
!                ! Take the absolut value of the determinant of the mapping.
!                ! In 2D, the determinant is always positive, whereas in 3D,
!                ! the determinant might be negative -- that is normal!
!                
!                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))*Dcoefficients(icubp,IEL,2)
!                
!                IELGlobal = p_IelementList(IELset+IEL-1)
!                
!                p_Derror(IELGlobal) = OM * Dcoefficients(icubp,IEL,1)
!
!                derror = derror + p_Derror(IELGlobal)
!                
!              end do ! ICUBP 
!              
!            end do ! IEL
!            
!          else
!            
!            do IEL = 1, IELmax-IELset+1
!              
!              ! Loop over all cubature points on the current element
!              do icubp = 1, ncubp
!                
!                ! calculate the current weighting factor in the cubature formula
!                ! in that cubature point.
!                !
!                ! Take the absolut value of the determinant of the mapping.
!                ! In 2D, the determinant is always positive, whereas in 3D,
!                ! the determinant might be negative -- that is normal!
!                
!                OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))*Dcoefficients(icubp,IEL,2)
!                
!                derror = derror + OM * Dcoefficients(icubp,IEL,1)
!                
!              end do ! ICUBP 
!              
!            end do ! IEL
!
!          end if
!
!        end select
!        
!        ! Release the temporary domain integration structure again
!        call domint_doneIntegration (rintSubset)
!    
!      end do ! IELset
!      
!      ! Release memory
!      call elprep_releaseElementSet(revalElementSet)
!
!      deallocate(p_DcubPtsRef)
!      deallocate(Dcoefficients)
!      deallocate(IdofsTrial)
!      deallocate(Domega)
!
!    end do ! ielementDistr
!
!    ! derror is ||error||^2, so take the square root at last.
!    if ((cerrortype .eq. PPERR_L2ERROR) .or.&
!        (cerrortype .eq. PPERR_H1ERROR)) then
!      derror = sqrt(derror)
!      if (present(relementError)) then
!        do IEL = 1, size(p_Derror,1)
!          p_Derror(IEL) = sqrt(p_Derror(IEL))
!        end do
!      end if
!    end if
!
!  end subroutine dg_pperr_scalar2d_conf
  
  
  


end module dg2d_routines
