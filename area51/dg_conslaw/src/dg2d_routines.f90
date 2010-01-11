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
  !                                            flux_dg_buildVectorScEdge2D_sim,&
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
   
  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  type(t_collection), intent(inout), target, optional :: rcollection
  
  ! A callback routine for the function to be discretised.
!  include 'intf_flux_dg_buildVectorScEdge2D.inc'
!  optional :: intf_flux_dg_buildVectorScEdge2D_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated entries are imposed to this vector.
  type(t_vectorScalar), intent(inout) :: rvectorScalar
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_linfVectorAssembly) :: rvectorAssembly
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
      
      ! All edges (the number of the first edge is NVT+1)
      !forall (iedge = 1:NMT) p_IedgeList(iedge)=NVT+iedge
      forall (iedge = 1:NMT) p_IedgeList(iedge)=iedge
      
      ! Allocate space for local edge numbers
      allocate(p_IlocalEdgeNumber(2,NMT))
      
      ! Get local edge numbers
      call linf_getLocalEdgeNumbers(p_IedgeList,p_IlocalEdgeNumber,&
                                    p_rtriangulation)
     
      call linf_initAssembly(rvectorAssembly, rform,&
            rvectorScalar%p_rspatialDiscr%RelementDistr(1)%celement,&
            ccubType, LINF_NELEMSIM)
       
      ! Assemble the data for all elements in this element distribution
 !     call linf_dg_assembleSubmeshVectorScalarEdge2d (rvectorAssembly, rvectorScalar,&
 !           p_IedgeList(1:NMT), flux_dg_buildVectorScEdge2D_sim, rcollection)
          
      ! Release the assembly structure.
      call linf_doneAssembly(rvectorAssembly)

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

















































!
!
!
!  !****************************************************************************
!  
!!<subroutine>  
!  
!  subroutine linf_dg_assembleSubmeshVectorScalarEdge2d (rvectorAssembly, rvector,&
!      rboundaryRegion, IelementList, IelementOrientation, DedgePosition,&
!      fcoeff_buildVectorScBdr2D_sim, rcollection)
!
!!<description>
!
!  ! Assembles the vector entries for a submesh by integration over the boundary region.
!
!!</description>
!
!!<input>
!  
!  ! A boundary region where to assemble the contribution
!  type(t_boundaryRegion), intent(in) :: rboundaryRegion
!
!  ! List of elements where to assemble the linear form.
!  integer, dimension(:), intent(in), target :: IelementList
!  
!  ! List of element orientations where to assemble the linear form.
!  integer, dimension(:), intent(in) :: IelementOrientation
!
!  ! List of start- and end-parameter values of the edges on the boundary
!  real(DP), dimension(:,:), intent(in) :: DedgePosition
!
!  ! A callback routine which is able to calculate the values of the
!  ! function $f$ which is to be discretised.
!  include 'intf_coefficientVectorScBdr2D.inc'
!  optional :: fcoeff_buildVectorScBdr2D_sim 
!  
!!</input>
!
!!<inputoutput>
!  
!  ! A vector assembly structure prepared with linf_initAssembly.
!  type(t_linfVectorAssembly), intent(inout), target :: rvectorAssembly
!  
!  ! A vector where to assemble the contributions to.
!  type(t_vectorScalar), intent(inout) :: rvector  
!  
!  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
!  ! callback function for nonconstant coefficients to provide additional
!  ! information. 
!  type(t_collection), intent(inout), target, optional :: rcollection
!!</inputoutput>
!  
!!</subroutine>
!
!    ! local variables, used by all processors
!    real(DP), dimension(:), pointer :: p_Ddata
!    integer :: indof,ncubp
!    
!    ! local data of every processor when using OpenMP
!    integer :: IELset,IELmax,ibdc,k
!    integer :: iel,icubp,ialbet,ia,idofe
!    real(DP) :: domega,daux,dlen
!    integer(I32) :: cevaluationTag
!    type(t_linfVectorAssembly), target :: rlocalVectorAssembly
!    type(t_domainIntSubset) :: rintSubset
!    real(DP), dimension(:), pointer :: p_Domega
!    real(DP), dimension(:,:,:,:), pointer :: p_Dbas
!    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
!    real(DP), dimension(:,:), pointer :: p_DcubPtsRef
!    integer, dimension(:),pointer :: p_Idescriptors
!    integer, dimension(:,:), pointer :: p_Idofs
!    type(t_evalElementSet), pointer :: p_revalElementSet
!
!    ! A small vector holding only the additive controbutions of
!    ! one element
!    real(DP), dimension(EL_MAXNBAS) :: DlocalData
!  
!    ! Arrays for cubature points 1D->2D
!    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
!    real(DP), dimension(:,:,:), allocatable :: Dxi2D,DpointsRef
!    real(DP), dimension(:,:), allocatable :: DpointsPar
!    
!    integer(i32) :: icoordSystem
!
!    ! Boundary component?
!    ibdc = rboundaryRegion%iboundCompIdx
!
!    ! Get some pointers for faster access
!    call lsyssc_getbase_double (rvector, p_Ddata)
!    indof = rvectorAssembly%indof
!    ncubp = rvectorAssembly%ncubp
!
!    ! Copy the assembly data to the local assembly data,
!    ! where we can allocate memory.
!    ! For single processor machines, this is actually boring and nonsense.
!    ! But using OpenMP, here we get a local copy of the vector
!    ! assembly structure to where we can add some local data which
!    ! is released upon return without changing the original assembly
!    ! stucture or disturbing the data of the other processors.
!    rlocalVectorAssembly = rvectorAssembly
!    call linf_allocAssemblyData(rlocalVectorAssembly)
!    
!    ! Get some more pointers to local data.
!    p_Domega => rlocalVectorAssembly%p_Domega
!    p_Dbas => rlocalVectorAssembly%p_Dbas
!    p_Dcoefficients => rlocalVectorAssembly%p_Dcoefficients
!    p_DcubPtsRef => rlocalVectorAssembly%p_DcubPtsRef
!    p_Idescriptors => rlocalVectorAssembly%rform%Idescriptors
!    p_Idofs => rlocalVectorAssembly%p_Idofs
!    p_revalElementSet => rlocalVectorAssembly%revalElementSet
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
!    allocate(Dxi2D(ncubp,NDIM2D+1,rlocalVectorAssembly%nelementsPerBlock))
!
!    ! Allocate memory for the coordinates of the reference points
!    allocate(DpointsRef(NDIM2D+1,ncubp,rlocalVectorAssembly%nelementsPerBlock))
!
!    ! Allocate memory for the parameter values of the points on the boundary
!    allocate(DpointsPar(ncubp,rlocalVectorAssembly%nelementsPerBlock))
!
!    ! Get the type of coordinate system
!    icoordSystem = elem_igetCoordSystem(rlocalVectorAssembly%celement)
!
!    ! Loop over the elements - blockwise.
!    !
!    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
!    ! so nelementsPerBlock local matrices are simultaneously calculated in the
!    ! inner loop(s).
!    ! The blocks have all the same size, so we can use static scheduling.
!    !
!    !%OMP do schedule(static,1)
!    do IELset = 1, size(IelementList), rlocalVectorAssembly%nelementsPerBlock
!    
!      ! We always handle nelementsPerBlock elements simultaneously.
!      ! How many elements have we actually here?
!      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
!      ! elements simultaneously.
!      
!      IELmax = min(size(IelementList),IELset-1+rlocalVectorAssembly%nelementsPerBlock)
!      
!      ! Map the 1D cubature points to the edges in 2D.
!      do iel = 1,IELmax-IELset+1
!        call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(IELset+iel-1), &
!            ncubp, Dxi1D, Dxi2D(:,:,iel))
!      end do
!
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
!      
!      ! Transpose the coordinate array such that we get coordinates we
!      ! can work with.
!      do iel = 1,IELmax-IELset+1
!        do icubp = 1,ncubp
!          do k = 1,ubound(DpointsRef,1)
!            DpointsRef(k,icubp,iel) = Dxi2D(icubp,k,iel)
!          end do
!        end do
!      end do
!      
!      ! --------------------- DOF SEARCH PHASE ------------------------
!    
!      ! The outstanding feature with finite elements is: A basis
!      ! function for a DOF on one element has common support only
!      ! with the DOF`s on the same element! E.g. for Q1:
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
!      !     to collect all "O" DOF`s.
!      !
!      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
!      !
!      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
!      ! global DOF`s of our LINF_NELEMSIM elements simultaneously.
!      call dof_locGlobMapping_mult(rvector%p_rspatialDiscr, &
!          IelementList(IELset:IELmax), p_Idofs)
!                                   
!      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
!      
!      ! To calculate the element contributions, we have to evaluate
!      ! the elements to give us the values of the basis functions
!      ! in all the DOF`s in all the elements in our set.
!
!      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
!      ! the elements later. All of them can be combined with OR, what will give
!      ! a combined evaluation tag. 
!      cevaluationTag = rlocalVectorAssembly%cevaluationTag
!      
!      ! The cubature points are already initialised by 1D->2D mapping.
!      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
!
!      ! Calculate all information that is necessary to evaluate the
!      ! finite element on all cells of our subset. This includes the
!      ! coordinates of the points on the cells.
!      call elprep_prepareSetForEvaluation (p_revalElementSet,&
!          cevaluationTag, rvector%p_rspatialDiscr%p_rtriangulation, &
!          IelementList(IELset:IELmax), rlocalVectorAssembly%ctrafoType, &
!          DpointsRef=DpointsRef)
!      
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
!      
!      ! Calculate the values of the basis functions.
!      call elem_generic_sim2 (rlocalVectorAssembly%celement, &
!          p_revalElementSet, rlocalVectorAssembly%Bder, &
!          rlocalVectorAssembly%p_Dbas)
!      
!      ! --------------------- DOF COMBINATION PHASE ------------------------
!      
!      ! Values of all basis functions calculated. Now we can start 
!      ! to integrate!
!      !
!      ! Loop through elements in the set and for each element,
!      ! loop through the DOF`s and cubature points to calculate the
!      ! integral:
!
!      do iel = 1,IELmax-IELset+1
!        
!        ! We make a 'local' approach, i.e. we calculate the values of the
!        ! integral into the vector DlocalData and add them later into
!        ! the large solution vector.
!        
!        ! Clear the output vector.
!        DlocalData(1:indof) = 0.0_DP
!
!        ! Get the length of the edge. Let us use the parameter values
!        ! on the boundary for that purpose; this is a more general
!        ! implementation than using simple lines as it will later 
!        ! support isoparametric elements.
!        !
!        ! The length of the current edge serves as a "determinant"
!        ! in the cubature, so we have to divide it by 2 as an edge on 
!        ! the unit interval [-1,1] has length 2.
!        dlen = 0.5_DP*(DedgePosition(2,IELset+iel-1)-DedgePosition(1,IELset+iel-1))
!
!        ! Loop over all cubature points on the current element
!        do icubp = 1, ncubp
!
!          ! Calculate the current weighting factor in the cubature
!          ! formula in that cubature point.
!
!          domega = dlen * p_Domega(icubp)
!          
!          ! Loop over the additive factors in the bilinear form.
!          do ialbet = 1,rlocalVectorAssembly%rform%itermcount
!          
!            ! Get from Idescriptors the type of the derivatives for the 
!            ! test and trial functions. The summand we calculate
!            ! here will be:
!            !
!            ! int_...  f * ( phi_i )_IA
!            !
!            ! -> IA=0: function value, 
!            !      =1: first derivative, 
!            !      =2: 2nd derivative,...
!            !    as defined in the module 'derivative'.
!            
!            ia = p_Idescriptors(ialbet)
!            
!            ! Multiply domega with the coefficient of the form.
!            ! This gives the actual value to multiply the
!            ! function value with before summing up to the integral.
!            ! Get the precalculated coefficient from the coefficient array.
!            daux = domega * p_Dcoefficients(ialbet,icubp,iel)
!          
!            ! Now loop through all possible combinations of DOF`s
!            ! in the current cubature point. 
!
!            do idofe = 1,indof
!              
!              ! Get the value of the basis function 
!              ! phi_o in the cubature point. 
!              ! Them multiply:
!              !    DBAS(..) * AUX
!              ! ~= phi_i * coefficient * cub.weight
!              ! Summing this up gives the integral, so the contribution
!              ! to the vector. 
!              !
!              ! Simply summing up DBAS(..) * AUX would give
!              ! the additive contribution for the vector. We save this
!              ! contribution in the local array.
!              
!              DlocalData(idofe) = DlocalData(idofe)+p_Dbas(idofe,ia,icubp,iel)*daux
!              
!            end do ! idofe
!            
!          end do ! ialbet
!
!        end do ! icubp 
!        
!        ! Incorporate the local vector into the global one.
!        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
!        ! the IdofsTest array.
!        do idofe = 1,indof
!          p_Ddata(p_Idofs(idofe,iel)) = p_Ddata(p_Idofs(idofe,iel)) + DlocalData(idofe)
!        end do
!
!      end do ! iel
!
!    end do ! IELset
!    
!    ! Release the local vector assembly structure
!    call linf_releaseAssemblyData(rlocalVectorAssembly)
!
!    ! Deallocate memory
!    deallocate(Dxi2D, DpointsRef, DpointsPar)
!
!  end subroutine
!













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



    
end module
