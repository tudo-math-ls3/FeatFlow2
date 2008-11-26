!#########################################################################
!# ***********************************************************************
!# <name> pprocnavierstokes </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains various routines for the postprocessing of 
!# solutions of the (Navier) Stokes system.
!#
!# The following routines can be found in this module:
!#
!# 1.) ppns2D_bdforces_uniform
!#     -> Calculate the body forces acting on a boundary component
!#        of a given parametrisation.
!#
!# 2.) ppns3D_bdforces_uniform
!#     -> Calculate the body forces acting on a boundary mesh region.
!#
!# 3.) ppns2D_streamfct_uniform
!#     -> Calculate the streamfunction of a 2D velocity field.
!# </purpose>
!#########################################################################

module pprocnavierstokes

  use fsystem
  use storage
  use boundary
  use cubature
  use triangulation
  use linearsystemscalar
  use linearsystemblock
  use spatialdiscretisation
  use derivatives
  use meshregion

  implicit none

contains

  !****************************************************************************

!<subroutine>

  subroutine ppns2D_bdforces_uniform (rvector,rregion,Dforces,ccub,df1,df2)

!<description>
  ! Calculates the drag-/lift-forces acting on a part of the real
  ! boundary for a vector rvector with the solution of the 2D
  ! (Navier-)Stokes equation. It's assumed that
  !   rvector%rvectorBlock(1) = X-velocity,
  !   rvector%rvectorBlock(2) = Y-velocity,
  !   rvector%rvectorBlock(3) = pressure.
  !
  ! rregion specifies a boundary region where to calculate
  ! the force integral. Dforces(1:NDIM2D) receives the forces in the
  ! X- and Y-direction.
  !
  ! The body forces are defined as the integrals
  !
  !    Dforces(1) = 2/df2 * int_s [df1 dut/dn n_y - p n_x] ds 
  !    Dforces(2) = 2/df2 * int_s [df1 dut/dn n_x + p n_y] ds 
  !
  ! where df1 and df2 allow to weight different integral parts
  !
  ! Usually in benchmark-like geometries there is:
  !
  !   $$ df1 = RHO*NU = (density of the fluid)*viscosity$$
  !   $$ df2 = RHO*DIST*UMEAN**2
  !          = (density of the fluid)*(length of the obstacle facing the flow)
  !           *(mean velocity of the fluid)^2 $$
  !
  ! which influences this routine to calculate the drag-/lift coefficients
  ! instead of the forces. The actual forces are calculated by
  ! setting df1=1.0 and df2=2.0, which is the standard setting when
  ! neglecting df1 and df2.
  !
  ! Double precision version for all Finite Elements without 
  ! curved boundaries.
!</description>

!<input>
  ! The FE solution vector.
  type(t_vectorBlock), intent(IN)    :: rvector
  
  ! Boundary region where to calculate the boundary forces.
  ! Can be created e.g. by boundary_createRegion.
  type(t_boundaryRegion), intent(IN) :: rregion
  
  ! 1D Cubature formula identifier to use for the line integration.
  ! One of the CUB_xxxx_1D constants in the cubature.f90.
  integer, intent(IN)                 :: ccub

  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(IN), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(IN), optional      :: df2
  
!</input>

!<output>
  ! Array receiving the forces acting on the boundary specified by rregion.
  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
  real(DP), dimension(:), intent(OUT) :: Dforces
!</output>

!</subroutine>

  ! Spatial discretisation structure of velocity and pressure
  type(t_spatialDiscretisation), pointer :: p_rdiscrU, p_rdiscrP

  ! Element type identifier for U and P
  integer(I32) :: ielemU, ielemP
  
  ! Number of local DOF's in U and P
  integer :: idoflocU, idoflocP
  
  ! Triangulation
  type (t_triangulation), pointer :: p_rtriangulation
  
  ! An accepting the DOF's of an element.
  integer, dimension(EL_MAXNBAS), target :: IdofsU, IdofsP
  
  ! Coordinates of the coordinates of an element
  real(DP), dimension(:,:), allocatable :: DCoords
  
  ! Coordinates of the cubature points on reference and real element
  real(DP), dimension(NDIM2D,CUB_MAXCUBP_1D) :: DpointsRef,DpointsReal
  
  ! Coordinate system for U and P element
  integer :: ctrafoU, ctrafoP
  
  ! U/P element parametric or nonparametric
  logical :: bnonparU,bnonparP
  
  ! Arrays for saving Jacobian determinants and matrices
  real(DP), dimension(CUB_MAXCUBP_1D) :: Ddetj
  real(DP), dimension(EL_NJACENTRIES2D,CUB_MAXCUBP_1D) :: Djac

  ! Array to tell the element which derivatives to calculate.
  logical, dimension(EL_MAXNDER) :: BderU, BderP
  
  ! Value of basis functions
  real(DP), dimension(EL_MAXNBAS,EL_MAXNDER,CUB_MAXCUBP_1D) :: DbasU, DbasP
  
  ! Pointer to vector data of solution vector
  real(DP), dimension(:), pointer :: p_DdataUX,p_DdataUY,p_DdataP

  ! Cubature point coordinates on the reference element.
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D, DXi2D

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(CUB_MAXCUBP) :: Domega
  
  ! number/index of cubature points on the reference element
  integer :: ncubp,icubp

  ! Edges, vertices and elements on the boundary
  integer, dimension(:), pointer :: p_IelementsAtBoundary
  integer, dimension(:), pointer    :: p_IedgesAtBoundary
  integer, dimension(:,:), pointer  :: p_IedgesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtEdge
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  real(DP), dimension(:), pointer                 :: p_DvertexParameterValue
  integer, dimension(:), pointer             :: p_IboundaryCpIdx
  real(DP), dimension(:,:), pointer               :: p_DvertexCoordinates

  ! other local variables
  integer :: iedgeidx,ivt1,ivt2,ilocaledge,nlocaledges,idfl,icp
  integer :: neqU,neqP
  integer :: iedge,iedgeglobal
  integer :: iel
  real(DP) :: dvtp1,dvtp2,dedgelen, dweight,dut,dpf1,dpf2
  real(DP), dimension(2) :: DintU, DintP
  real(DP) :: dpres
  real(DP), dimension(NDIM2D) :: dvt1,dvt2,dtangential,dnormal
  integer(I32), dimension(:), pointer :: p_ItwistIndex

    ! Get the vector data
    neqU = rvector%RvectorBlock(1)%NEQ
    neqP = rvector%RvectorBlock(3)%NEQ
    
    if (rvector%cdataType .ne. ST_DOUBLE) then
      print *,'ppns2D_bdforces: Unsupported vector precision.'
      call sys_halt()
    end if

    ! We support only uniform discretisation structures.
    if (.not. associated(rvector%p_rblockDiscr)) then
      print *,'ppns2D_bdforces: No discretisation structure!'
      call sys_halt()
    end if

    if (rvector%p_rblockDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
      print *,'ppns2D_bdforces_uniform: Discretisation too complex!'
      call sys_halt()
    end if
    
    ! Get pointers to the subvectors from the block vector
    call lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
    call lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
    call lsyssc_getbase_double (rvector%RvectorBlock(3),p_DdataP)
    
    if ((rvector%RvectorBlock(1)%isortStrategy > 0) .or. &
        (rvector%RvectorBlock(2)%isortStrategy > 0) .or. &
        (rvector%RvectorBlock(3)%isortStrategy > 0)) then
      print *,'ppns2D_bdforces_uniform: Resorted vectors not supported!'
      call sys_halt()
    end if
    
    ! Get pointers to the spatial discretisation structures of the
    ! veloctiy and pressure
    p_rdiscrU => rvector%RvectorBlock(1)%p_rspatialDiscr
    p_rdiscrP => rvector%RvectorBlock(3)%p_rspatialDiscr
    
    ! What is the actual element that is used for the discretisation?
    
    ielemU = p_rdiscrU%RelementDistr(1)%celement
    ielemP = p_rdiscrP%RelementDistr(1)%celement
    
    ! So far so good, we have checked that the assumptions for the integration
    ! are fulfilled. Now we can start the actual integration.
    !
    ! At first initialise a 1D cubature formula for the boundary integration.
    
    call cub_getCubPoints (ccub,ncubp,Dxi1D,Domega)
    
    ! In Dxi1D we have the 1D coordinates of the cubature points.
    ! These have to be mapped to the 2D element which is under condideration
    ! during the integration -- later.
    !
    ! Before, we need some additional information about the triangulation
    ! and about our element!
    !
    ! Number of local DOF's:
    
    idoflocU = elem_igetNDofLoc(ielemU)
    idoflocP = elem_igetNDofLoc(ielemP)
    
    ! Number of local edges
    nlocaledges = elem_igetNVE (ielemU)
    
    ! Allocate memory for the coordinates of the element
    allocate(DCoords(NDIM2D,max(elem_igetNVE(ielemU),elem_igetNVE(ielemP))))
    
    ! The triangulation - it's the same for U and P
    p_rtriangulation => p_rdiscrU%p_rtriangulation

    ! The arrays which contain the elements and edge numbers on the boundary
    ! as well as boundary index array.
    ! Fetch geometrical information
    call storage_getbase_int (p_rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement, &
        p_IverticesAtElement)
    call storage_getbase_int (p_rtriangulation%h_IelementsAtBoundary, &
        p_IelementsAtBoundary)
    call storage_getbase_int (p_rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords, &
                                   p_DvertexCoordinates)
                                   
    if (p_rtriangulation%h_DvertexParameterValue .eq. ST_NOHANDLE) then
      print *,'No boundary parameters available!'
      call sys_halt()
    end if
    
    call storage_getbase_double (p_rtriangulation%h_DvertexParameterValue, &
                                 p_DvertexParameterValue)
                                 
    ! Does the element need twist indices?
    nullify(p_ItwistIndex)
    if (p_rtriangulation%h_ItwistIndex .ne. ST_NOHANDLE) then
      call storage_getbase_int32 (p_rtriangulation%h_ItwistIndex,p_ItwistIndex)
    end if

    ! Is one of the elements nonparametric
    bnonparU = elem_isnonparametric(ielemU) 
    bnonparP = elem_isnonparametric(ielemP)

    ! Coordinate systems of U and P element
    ctrafoU = elem_igetTrafoType(ielemU)
    ctrafoP = elem_igetTrafoType(ielemP)
    
    ! Derivatives to calculate when evaluating the U and P-element, respectively.
    BderU = .false.
    BderU(DER_DERIV_X) = .true.
    BderU(DER_DERIV_Y) = .true.
    
    BderP = .false.
    BderP(DER_FUNC) = .true.
    
    ! Prepare the weighting coefficients
    dpf1 = 1.0_DP
    dpf2 = 2.0_DP
    if (present(df1)) dpf1 = df1
    if (present(df2)) dpf2 = df2

    ! We are on boundary component
    icp = rregion%iboundCompIdx
    
    ! We assemble the integral contributions separately
    DintU = 0.0_DP
    DintP = 0.0_DP
    
    ! Loop through all edges along the boundary
    do iedgeidx = p_IboundaryCpIdx(icp),p_IboundaryCpIdx(icp+1)-1
      
      ! Current element
      iel = p_IelementsAtBoundary(iedgeidx)
      
      ! Egde number
      iedgeglobal = p_IedgesAtBoundary(iedgeidx)
      iedge       = iedgeglobal
      
      ! What are the vertices adjacent to that edge?
      ivt1 = p_IverticesAtEdge(1,iedge)
      ivt2 = p_IverticesAtEdge(2,iedge)
      
      ! Parameter values of these vertices?
      ! iedgeidx is the index of the edge as well as the number of the
      ! index of the vertex preceding the edge!
      dvtp1 = p_DvertexParameterValue (iedgeidx)
      if (iedgeidx .ne. p_IboundaryCpIdx(icp+1)-1) then
        dvtp2 = p_DvertexParameterValue (iedgeidx+1)
      else
        ! Last vertex has maximum parameter value! (TMAX)
        dvtp2 = boundary_dgetMaxParVal(p_rdiscrU%p_rboundary,icp)
      end if
      
      ! Is the edge in the specified boundary region, so we are allowed 
      ! to integrate about it? We check both endpoints...
      
      if (boundary_isInRegion (rregion,icp,dvtp1) .and. &
          boundary_isInRegion (rregion,icp,dvtp2)) then
          
        ! Ok, the edge is a boundary edge in the specified region.
        !
        ! Get the coordinates of the vertices
        dvt1 = p_DvertexCoordinates(:,ivt1)
        dvt2 = p_DvertexCoordinates(:,ivt2)
        
        ! Furthermore, get the tangential and the normal vector of the edge
        dtangential = dvt2(:)-dvt1(:)
        
        dedgelen = sqrt(dtangential(1)**2 + dtangential(2)**2)
        
        dtangential = dtangential/ dedgelen
        
        dnormal(1) = -dtangential(2)
        dnormal(2) = dtangential(1)

        ! Note that this only works for a linear approximation of
        ! the boundary! In a later implementation when finite elements
        ! are allowed to have curved boundaries, this has to be
        ! changed! In that case, a normal/tangential vector has to
        ! be computed in every cubature point along the line.
        ! The length of the boundary edge can be obtained using
        ! the arc length parametrisation of the boundary...
        !
        ! Now quickly determine which edge 1..4 of the current element
        ! is the one we are treating at the moment; the edge starts with
        ! vertex IVT1!
        
        do ilocaledge = 1,nlocaledges
          if (p_IedgesAtElement(ilocaledge,iel) .eq. iedgeglobal) exit
        end do
        
        if (ilocaledge .gt. nlocaledges) then
          print *,'ppns2D_bdforces: Edge not found. KMID destroyed?'
          call sys_halt()
        end if
        
        ! The number of the edge is in ilocaledge. We have to transfer
        ! the coordinates of the cubature points from 1D to 2D depending
        ! on this edge.
        
        call trafo_mapCubPts1Dto2DRefQuad(ilocaledge, ncubp, Dxi1D, Dxi2D)

        ! Now, in Dxi2D we have all the cubature points on the 2D reference
        ! element along the edge!
        ! Map the coordinates in a proper 2D array
        DpointsRef (1:NDIM2D,1:ncubp) = transpose(Dxi2D(1:ncubp,1:NDIM2D))
                
        ! For the integration, we need the global DOF's on our element
        ! for U and P:
        
        call dof_locGlobMapping(p_rdiscrU, iel, IdofsU)
        call dof_locGlobMapping(p_rdiscrP, iel, IdofsP)
        
        ! Get the coordinates of the point on the current element
        Dcoords (1:NDIM2D,1:nlocaledges) = &
            p_DvertexCoordinates(1:NDIM2D, p_IverticesAtElement (1:nlocaledges,iel))
        
        ! Calculate the transformation for all points on the current element.
        ! If we have only parametric elements, we don't have to calculate
        ! the real coordinates of the points.
        if (bnonparU .or. bnonparP) then
          call trafo_calctrafo_mult (ctrafoU,ncubp,Dcoords,&
                                     DpointsRef,Djac,Ddetj,DpointsReal)
        else
          call trafo_calctrafo_mult (ctrafoU,ncubp,Dcoords,&
                                     DpointsRef,Djac,Ddetj)
        end if
        
        ! Evaluate the U- and P-element in all our cubature points
        if (bnonparU) then
          call elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsReal,p_ItwistIndex(iel))
        else
          call elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsRef,p_ItwistIndex(iel))
        end if

        if (bnonparP) then
          call elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsReal,p_ItwistIndex(iel))
        else
          call elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsRef,p_ItwistIndex(iel))
        end if
        
        ! Loop over the cubature points on the current element
        ! to assemble the integral
        do icubp = 1,ncubp
        
          ! Calculate the OMEGA for the integration by multiplication
          ! of the integration coefficient by the Jacobian of the
          ! mapping.
          ! The determinant of the mapping of the unit interval [-1,1]
          ! to the real line is 0.5*length of the line!
            
          dweight = Domega(icubp)*0.5_DP*dedgelen
          
          ! Loop through the DOF's on our element and calculate
          ! the tangential U as well as P.
          dut = 0.0_DP
          do idfl=1,idoflocU
            dut = dut &
            + p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)*Dtangential(1) &
                                      * Dnormal(1) &
            + p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)*Dtangential(2) &
                                      * Dnormal(1) &
            + p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)*Dtangential(1) &
                                      * Dnormal(2) &
            + p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)*Dtangential(2) &
                                      * Dnormal(2)
          end do

          dpres = 0.0_DP
          do idfl=1,idoflocP
            dpres = dpres+p_DdataP(IdofsP(idfl))*DbasP(idfl,DER_FUNC,icubp)
          end do
          
          ! Sum this up to the two integral contributions for the pressure and
          ! velocity.
          DintU(1) = DintU(1) + dweight * dut * Dnormal(2)
          DintU(2) = DintU(2) - dweight * dut * Dnormal(1)
          DintP(1) = DintP(1) - dweight * dpres * Dnormal(1)
          DintP(2) = DintP(2) - dweight * dpres * Dnormal(2)
        
        end do ! icubp
          
      end if
      
    end do ! iedgeidx
    
    ! DintU and DintP give now the contributions to the force integral:
    !
    ! DragCoeff = 2/dfp2 * (dpf1*DintU(1) + DintP(1))
    ! LiftCoeff = 2/dfp2 * (dpf1*DintU(2) + DintP(2))
    
    Dforces = 0.0_DP
    Dforces(1:NDIM2D) = 2.0_DP/dpf2 * (dpf1*DintU(:) + DintP(:))
    
    ! Deallocate memory, finish.
    deallocate(DCoords)
  
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine ppns3D_bdforces_uniform (rvector,rregion,Dforces,ccub,df1,df2)

!<description>
  ! Calculates the drag-/lift-forces acting on a part of the real
  ! boundary for a vector rvector with the solution of the 3D
  ! (Navier-)Stokes equation. It's assumed that
  !   rvector%rvectorBlock(1) = X-velocity,
  !   rvector%rvectorBlock(2) = Y-velocity,
  !   rvector%rvectorBlock(3) = Z-velocity,
  !   rvector%rvectorBlock(4) = pressure.
  !
  ! rregion specifies a boundary region where to calculate
  ! the force integral. Dforces(1:3) receives the forces in the
  ! X-, Y- and Z-direction.
  !
  ! The body forces are defined as the integrals
  !
  !    Dforces(1) = 2/df2 * int_s [df1 dut/dn n_y - p n_x] ds 
  !    Dforces(2) = 2/df2 * int_s [df1 dut/dn n_x + p n_y] ds 
  !    Dforces(3) = 2/df2 * int_s [df1 dut/dn n_z + p n_z] ds 
  !
  ! where df1 and df2 allow to weight different integral parts
  !
  ! Usually in benchmark-like geometries there is:
  !
  !   $$ df1 = RHO*NU = (density of the fluid)*viscosity$$
  !   $$ df2 = RHO*DIST*UMEAN**2
  !          = (density of the fluid)*(length of the obstacle facing the flow)
  !           *(mean velocity of the fluid)^2 $$
  !
  ! which influences this routine to calculate the drag-/lift coefficients
  ! instead of the forces. The actual forces are calculated by
  ! setting df1=1.0 and df2=2.0, which is the standard setting when
  ! neglecting df1 and df2.
  !
  ! Double precision version for all Finite Elements without 
  ! curved boundaries.
!</description>

!<input>
  ! The FE solution vector.
  type(t_vectorBlock), intent(IN)     :: rvector
  
  ! Mesh region where to calculate the boundary forces.
  type(t_meshRegion), intent(IN)      :: rregion
  
  ! 2D Cubature formula identifier to use for the quad integration.
  ! One of the CUB_xxxx constants in the cubature.f90.
  integer, intent(IN)                 :: ccub

  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(IN), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(IN), optional      :: df2
  
!</input>

!<output>
  ! Array receiving the forces acting on the boundary specified by rregion.
  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
  real(DP), dimension(:), intent(OUT) :: Dforces
!</output>

!</subroutine>

  ! Spatial discretisation structure of velocity and pressure
  type(t_spatialDiscretisation), pointer :: p_rdiscrU, p_rdiscrP

  ! Element type identifier for U and P
  integer(I32) :: ielemU, ielemP
  
  ! Number of local DOF's in U and P
  integer :: idoflocU, idoflocP
  
  ! Triangulation
  type (t_triangulation), pointer :: p_rtria
  
  ! An accepting the DOF's of an element.
  integer, dimension(EL_MAXNBAS), target :: IdofsU, IdofsP
  
  ! Coordinates of the vertices
  real(DP), dimension(NDIM3D,8) :: Dcoords
  
  ! Coordinates of the normal vectors
  real(DP), dimension(:,:), allocatable :: Dnormal

  ! Coordinates of the cubature points on reference and real element
  real(DP), dimension(NDIM3D,CUB_MAXCUBP_2D) :: DpointsRef,DpointsReal
  
  ! Coordinate system for U and P element
  integer :: ctrafoU, ctrafoP
  
  ! U/P element parametric or nonparametric
  logical :: bnonparU,bnonparP
  
  ! Arrays for saving Jacobian determinants and matrices
  real(DP), dimension(CUB_MAXCUBP_2D) :: Ddetj
  real(DP), dimension(CUB_MAXCUBP_2D) :: Ddetj_face
  real(DP), dimension(EL_NJACENTRIES3D,CUB_MAXCUBP_2D) :: Djac

  ! Array to tell the element which derivatives to calculate.
  logical, dimension(EL_MAXNDER) :: BderU, BderP
  
  ! Value of basis functions
  real(DP), dimension(EL_MAXNBAS,EL_MAXNDER,CUB_MAXCUBP_2D) :: DbasU, DbasP
  
  ! Pointer to vector data of solution vector
  real(DP), dimension(:), pointer :: p_DdataUX,p_DdataUY,p_DdataUZ,p_DdataP

  ! Cubature point coordinates on the reference element.
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi2D, Dxi3D

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(CUB_MAXCUBP) :: Domega
  
  ! number/index of cubature points on the reference element
  integer :: ncubp,icubp

  ! Edges, vertices and elements on the boundary
  integer, dimension(:), pointer     :: p_IfaceIdx
  integer, dimension(:,:), pointer   :: p_IfaceAtElem
  integer, dimension(:,:), pointer :: p_IvertAtElem
  integer, dimension(:,:), pointer :: p_IelemAtFace
  real(DP), dimension(:,:), pointer                :: p_Dvertex

  ! other local variables
  integer :: iat, iface,ilocalface,idfl,i
  integer :: neqU,neqP
  integer :: iel
  real(DP) :: dpf1,dpf2,dpres,du1,du2,du3,dweight,dv
  real(DP), dimension(3) :: DintU, DintP
  integer(I32), dimension(:), pointer :: p_ItwistIndex

    ! Get the vector data
    neqU = rvector%RvectorBlock(1)%NEQ
    neqP = rvector%RvectorBlock(4)%NEQ
    
    if (rvector%cdataType .ne. ST_DOUBLE) then
      print *,'ppns3D_bdforces: Unsupported vector precision.'
      call sys_halt()
    end if

    ! We support only uniform discretisation structures.
    if (.not. associated(rvector%p_rblockDiscr)) then
      print *,'ppns3D_bdforces: No discretisation structure!'
      call sys_halt()
    end if

    if (rvector%p_rblockDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
      print *,'ppns3D_bdforces_uniform: Discretisation too complex!'
      call sys_halt()
    end if
    
    ! Get pointers to the subvectors from the block vector
    call lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
    call lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
    call lsyssc_getbase_double (rvector%RvectorBlock(3),p_DdataUZ)
    call lsyssc_getbase_double (rvector%RvectorBlock(4),p_DdataP)
    
    if ((rvector%RvectorBlock(1)%isortStrategy > 0) .or. &
        (rvector%RvectorBlock(2)%isortStrategy > 0) .or. &
        (rvector%RvectorBlock(3)%isortStrategy > 0) .or. &
        (rvector%RvectorBlock(4)%isortStrategy > 0)) then
      print *,'ppns3D_bdforces_uniform: Resorted vectors not supported!'
      call sys_halt()
    end if
    
    ! Get pointers to the spatial discretisation structures of the
    ! veloctiy and pressure
    p_rdiscrU => rvector%RvectorBlock(1)%p_rspatialDiscr
    p_rdiscrP => rvector%RvectorBlock(4)%p_rspatialDiscr
    
    ! What is the actual element that is used for the discretisation?
    
    ielemU = p_rdiscrU%RelementDistr(1)%celement
    ielemP = p_rdiscrP%RelementDistr(1)%celement
    
    ! So far so good, we have checked that the assumptions for the integration
    ! are fulfilled. Now we can start the actual integration.
    !
    ! At first initialise a 2D cubature formula for the boundary integration.
    call cub_getCubPoints (ccub,ncubp,Dxi2D,Domega)
    
    ! In Dxi2D we have the 2D coordinates of the cubature points.
    ! These have to be mapped to the 3D element which is under condideration
    ! during the integration -- later.
    !
    ! Before, we need some additional information about the triangulation
    ! and about our element!
    !
    ! Number of local DOF's:
    idoflocU = elem_igetNDofLoc(ielemU)
    idoflocP = elem_igetNDofLoc(ielemP)
    
    ! The triangulation - it's the same for U and P
    p_rtria => p_rdiscrU%p_rtriangulation
    
    ! Get a pointer to the elements-at-face array
    call storage_getbase_int2D(p_rtria%h_IelementsAtFace, p_IelemAtFace)
       
    ! Get a pointer to the faces-at-element array
    call storage_getbase_int2D(p_rtria%h_IfacesAtElement, p_IfaceAtElem)

    ! Get a pointer to the vertices-at-element array
    call storage_getbase_int2D(p_rtria%h_IverticesAtElement, p_IvertAtElem)
    
    ! Get the vertice coordinate array
    call storage_getbase_double2D(p_rtria%h_DvertexCoords, p_Dvertex)
    
    ! And get the face index array of the mesh region
    call storage_getbase_int(rregion%h_IfaceIdx, p_IfaceIdx)
    
    ! Does the element need twist indices?
    nullify(p_ItwistIndex)
    if (p_rtria%h_ItwistIndex .ne. ST_NOHANDLE) then
      call storage_getbase_int32 (p_rtria%h_ItwistIndex,p_ItwistIndex)
    end if

    ! Is one of the elements nonparametric
    bnonparU = elem_isnonparametric(ielemU) 
    bnonparP = elem_isnonparametric(ielemP)

    ! Coordinate systems of U and P element
    ctrafoU = elem_igetTrafoType(ielemU)
    ctrafoP = elem_igetTrafoType(ielemP)
    
    ! Derivatives to calculate when evaluating the U and P-element, respectively.
    BderU = .false.
    BderU(DER_DERIV3D_X) = .true.
    BderU(DER_DERIV3D_Y) = .true.
    BderU(DER_DERIV3D_Z) = .true.
    
    BderP = .false.
    BderP(DER_FUNC3D) = .true.
    
    ! Prepare the weighting coefficients
    dpf1 = 1.0_DP
    dpf2 = 2.0_DP
    if (present(df1)) dpf1 = df1
    if (present(df2)) dpf2 = df2

    ! We assemble the integral contributions separately
    DintU = 0.0_DP
    DintP = 0.0_DP
    
    ! Calculate the normal vectors on the faces
    allocate(Dnormal(3, rregion%NAT))
    call mshreg_calcBoundaryNormals3D(rregion, Dnormal)
    
    ! Loop through all faces along the boundary
    do iat = 1, rregion%NAT
    
      ! Get the index fo the face
      iface = p_IfaceIdx(iat)
      
      ! Current element
      iel = p_IelemAtFace(1,iface)
      
      ! Get the coordinates of the corner vertices on the current element
      do i=1, 8
        Dcoords (:,i) = p_Dvertex(:, p_IvertAtElem (i,iel))
      end do

      ! Get the local index of the face on that element
      do ilocalface = 1, 6
        if(p_IfaceAtElem(ilocalface,iel) .eq. iface) exit
      end do
      
      ! We have to transfer the coordinates of the cubature points from
      ! 2D to 3D depending on this localface.
      call trafo_mapCubPts2Dto3DRefHexa(ilocalface, ncubp, Dxi2D, Dxi3D)
      
      ! And calculate the determinants for the mapping
      call ppns3D_det(Dcoords,ilocalface,Dxi2D,ncubp,Ddetj_face)

      ! Now, in Dxi3D we have all the cubature points on the 3D reference
      ! element along the face!
      ! Map the coordinates in a proper 3D array
      DpointsRef (1:NDIM3D,1:ncubp) = transpose(Dxi3D(1:ncubp,1:NDIM3D))
              
      ! For the integration, we need the global DOF's on our element
      ! for U and P:
      call dof_locGlobMapping(p_rdiscrU, iel, IdofsU)
      call dof_locGlobMapping(p_rdiscrP, iel, IdofsP)
      
      ! Calculate the transformation for all points on the current element.
      ! If we have only parametric elements, we don't have to calculate
      ! the real coordinates of the points.
      if (bnonparU .or. bnonparP) then
        call trafo_calctrafo_mult (ctrafoU,ncubp,Dcoords,&
                                   DpointsRef,Djac,Ddetj,DpointsReal)
      else
        call trafo_calctrafo_mult (ctrafoU,ncubp,Dcoords,&
                                   DpointsRef,Djac,Ddetj)
      end if
      
      if(associated(p_ItwistIndex)) then
      
        if (bnonparU) then
          call elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsReal,p_ItwistIndex(iel))
        else
          call elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsRef,p_ItwistIndex(iel))
        end if

        if (bnonparP) then
          call elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsReal,p_ItwistIndex(iel))
        else
          call elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsRef,p_ItwistIndex(iel))
        end if
        
      else
      
        if (bnonparU) then
          call elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsReal)
        else
          call elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsRef)
        end if

        if (bnonparP) then
          call elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsReal)
        else
          call elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsRef)
        end if
      
      end if
      
      ! Loop over the cubature points on the current element
      ! to assemble the integral
      do icubp = 1,ncubp
      
        ! Calculate the OMEGA for the integration by multiplication
        ! of the integration coefficient by the Jacobian of the
        ! mapping.
        dweight = Domega(icubp)*Ddetj_face(icubp)
        
        ! Loop through the DOF's on our element and calculate
        ! U as well as P.
        du1 = 0.0_DP
        du2 = 0.0_DP
        du3 = 0.0_DP
        do idfl=1,idoflocU
        
           dv = DbasU(idfl,DER_DERIV3D_X,icubp)*Dnormal(1,iat)&
              + DbasU(idfl,DER_DERIV3D_Y,icubp)*Dnormal(2,iat)&
              + DbasU(idfl,DER_DERIV3D_Z,icubp)*Dnormal(3,iat)
              
           du1 = du1 + p_DdataUX(IdofsU(idfl)) * dv
           du2 = du2 + p_DdataUY(IdofsU(idfl)) * dv
           du3 = du3 + p_DdataUZ(IdofsU(idfl)) * dv
               
        end do

        dpres = 0.0_DP
        do idfl=1,idoflocP
          dpres = dpres + p_DdataP(IdofsP(idfl))*DbasP(idfl,DER_FUNC3D,icubp)
        end do
        
        ! Sum this up to the two integral contributions for the pressure and
        ! velocity.
        DintU(1) = DintU(1) + dweight * du1
        DintU(2) = DintU(2) + dweight * du2
        DintU(3) = DintU(3) + dweight * du3
        DintP(1) = DintP(1) - dweight * dpres * Dnormal(1,iat)
        DintP(2) = DintP(2) - dweight * dpres * Dnormal(2,iat)
        DintP(3) = DintP(3) - dweight * dpres * Dnormal(3,iat)
      
      end do ! icubp
          
    end do ! iat
    
    ! DintU and DintP give now the contributions to the force integral:
    !
    ! DragCoeff = 2/dfp2 * (dpf1*DintU(1) + DintP(1))
    ! LiftCoeff = 2/dfp2 * (dpf1*DintU(2) + DintP(2))
    Dforces = 0.0_DP
    Dforces(:) = 2.0_DP/dpf2 * (dpf1*DintU(:) + DintP(:))
    
    ! Deallocate memory, finish.
    deallocate(Dnormal)
    
    ! That's it
    
    contains
    
!<subroutine>

    pure subroutine ppns3D_det(Dcoords,iface,Dcub,ncubp,Ddetj)

!<describtion>
  ! This routine calculates the "determinant" of a bilinear quadrilateral
  ! transformation from the 2D reference quadrilateral onto a hexahedron
  ! face in 3D.
!</describtion>
    
!<input>
    ! The coordinates of the eight corner vertices of the hexahedron
    real(DP), dimension(:,:), intent(IN) :: Dcoords
    
    ! The index of the face onto which the points are mapped
    integer, intent(IN) :: iface
    
    ! The 2D coordinates of the points that are to be mapped
    real(DP), dimension(:,:), intent(IN) :: Dcub
    
    ! The number of points which are to be mapped
    integer, intent(IN) :: ncubp
!</input>

!<output>
    ! The jacobian determinants of the mapping
    real(DP), dimension(:), intent(OUT) :: Ddetj
!</output>

!</subroutine>
    
    ! Local variables
    real(DP), dimension(3,4) :: Dv
    real(DP), dimension(3,3) :: Dt
    real(DP), dimension(3) :: Dx,Dy,Dn
    integer :: i,j
    
      ! Onto which face do we map the points?
      ! Note: The orientation of the vertices corresponds to the mapping
      ! routine trafo_mapCubPts2Dto3DRefHexa defined in transformation.f90.
      select case(iface)
      case (1)
        Dv(:,1) = Dcoords(:,1)
        Dv(:,2) = Dcoords(:,2)
        Dv(:,3) = Dcoords(:,3)
        Dv(:,4) = Dcoords(:,4)
      case (2)
        Dv(:,1) = Dcoords(:,1)
        Dv(:,2) = Dcoords(:,2)
        Dv(:,3) = Dcoords(:,6)
        Dv(:,4) = Dcoords(:,5)
      case (3)
        Dv(:,1) = Dcoords(:,2)
        Dv(:,2) = Dcoords(:,3)
        Dv(:,3) = Dcoords(:,7)
        Dv(:,4) = Dcoords(:,6)
      case (4)
        Dv(:,1) = Dcoords(:,3)
        Dv(:,2) = Dcoords(:,4)
        Dv(:,3) = Dcoords(:,8)
        Dv(:,4) = Dcoords(:,7)
      case (5)
        Dv(:,1) = Dcoords(:,4)
        Dv(:,2) = Dcoords(:,1)
        Dv(:,3) = Dcoords(:,5)
        Dv(:,4) = Dcoords(:,8)
      case (6)
        Dv(:,1) = Dcoords(:,5)
        Dv(:,2) = Dcoords(:,6)
        Dv(:,3) = Dcoords(:,7)
        Dv(:,4) = Dcoords(:,8)
      end select
      
      ! We have a bilinear mapping T: R^2 -> R^3, so the jacobian matrix
      ! of T is a 3x2 matrix. To get a useful replacement for the determinant
      ! we set the determinant to ||(dT/dx) X (dT/dy)||_2, where 'X' denotes
      ! the 3D cross-product.
      
      ! Calculate transformation coefficients for the jacobian matrix
      do i = 1,3
        Dt(i,1) = 0.25_DP * (-Dv(i,1) + Dv(i,2) + Dv(i,3) - Dv(i,4))
        Dt(i,2) = 0.25_DP * (-Dv(i,1) - Dv(i,2) + Dv(i,3) + Dv(i,4))
        Dt(i,3) = 0.25_DP * ( Dv(i,1) - Dv(i,2) + Dv(i,3) - Dv(i,4))
      end do
      
      ! And calculate the determinants
      do i = 1, ncubp
      
        do j = 1, 3
          ! Dx := dT / dx
          Dx(j) = Dt(j,1) + Dt(j,3)*Dcub(i,2)
          ! Dy := dT / dy
          Dy(j) = Dt(j,2) + Dt(j,3)*Dcub(i,1)
        end do
        
        ! Dn := Dx x Dy
        Dn(1) = Dx(2)*Dy(3) - Dx(3)*Dy(2)
        Dn(2) = Dx(3)*Dy(1) - Dx(1)*Dy(3)
        Dn(3) = Dx(1)*Dy(2) - Dx(2)*Dy(1)
        
        ! detj := ||Dn||_2
        Ddetj(i) = sqrt(Dn(1)**2 + Dn(2)**2 + Dn(3)**2)
        
      end do
    
    end subroutine
  
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine ppns2D_streamfct_uniform (rvector,rdestVector)

!<description>
  ! Calculates the streamfunction of a 2D velocity field rvector.
  ! It's assumed that
  !   rvector%rvectorBlock(1) = X-velocity,
  !   rvector%rvectorBlock(2) = Y-velocity
  !
  ! The vector rdestVector must be an empty vector and receives a Finite
  ! Element representation of the values of the Streamfunction. For this
  ! purpose the discretisation structure in rdestVector must specify
  ! a valid FE space.
  !
  ! Note 1: Supports only uniform discretisations, unsorted
  ! double precision vectors.
  !
  ! Note 2: Currently, this routine only supports Q1~ for the velocity
  ! subvectors in rvector. rdestVector is calculated in the Q1 space,
  ! so the discretisation structure in rdestVector must describbe a
  ! $Q_1$ discretisation!
!</description>

!<input>
  ! The FE solution vector.
  type(t_vectorBlock), intent(IN)    :: rvector
!</input>

!<inputoutput>
  ! An empty vector that is prepared with a discretisation structure to
  ! represent the streamfunction. The values in the vector are overwritten
  ! with the FE representation of the streamfunction.
  type(t_vectorScalar), intent(INOUT),target    :: rdestVector
!</inputoutput>

!</subroutine>

    ! local variables
    integer, parameter :: NVE = 4

    integer :: iel,ielaux,icurrentelement
    integer :: jve
    integer(I32) :: ieltype1,ieltype2,ieltypeDest
    integer :: haux,ive,iadj
    integer :: ilastMarked,imarkCounter,imarktmp
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to vector data of solution vector
    real(DP), dimension(:), pointer :: p_DdataUX,p_DdataUY,p_Dx
    integer, dimension(:), pointer :: p_Iind
    
    ! Stuff from the triangulation
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    if (rvector%cdataType .ne. ST_DOUBLE) then
      print *,'ppns2D_streamfct_uniform: Unsupported vector precision.'
      call sys_halt()
    end if

    ! We support only uniform discretisation structures.
    if (.not. associated(rvector%p_rblockDiscr)) then
      print *,'ppns2D_streamfct_uniform: No discretisation structure in rvector!'
      call sys_halt()
    end if

    if (rvector%p_rblockDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
      print *,'ppns2D_streamfct_uniform: Discretisation of rvector too complex!'
      call sys_halt()
    end if

    if (.not. associated(rdestVector%p_rspatialDiscr)) then
      print *,'ppns2D_streamfct_uniform: No discretisation structure in rdestVector!'
      call sys_halt()
    end if

    if (rdestVector%p_rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
      print *,'ppns2D_streamfct_uniform: Discretisation of rdestVector too complex!'
      call sys_halt()
    end if

    ieltype1 = rvector%p_rblockDiscr% &
        RspatialDiscr(1)%RelementDistr(1)%celement
    ieltype2 = rvector%p_rblockDiscr% &
        RspatialDiscr(2)%RelementDistr(1)%celement
    ieltypeDest = rdestVector%p_rspatialDiscr% &
        RelementDistr(1)%celement

    if (elem_getPrimaryElement(ieltype1) .ne. EL_Q1T) then
      print *,'ppns2D_streamfct_uniform: rvector must be discretised with Q1~!'
    end if

    if (elem_getPrimaryElement(ieltype2) .ne. EL_Q1T) then
      print *,'ppns2D_streamfct_uniform: rvector must be discretised with Q1~!'
    end if

    if (ieltypeDest .ne. EL_Q1) then
      print *,'ppns2D_streamfct_uniform: rdestVector must be discretised with Q1!'
    end if
    
    if ((rvector%RvectorBlock(1)%isortStrategy > 0) .or. &
        (rvector%RvectorBlock(2)%isortStrategy > 0) .or. &
        (rdestVector%isortStrategy > 0)) then
      print *,'ppns2D_bdforces_uniform: Resorted vectors not supported!'
      call sys_halt()
    end if

    ! Let's go. Note that we perform the same computation for all,
    ! parametric and nonparametric, point-based and integral-mean-value
    ! based Q1~ solutions. Gives a slight error, but the final Q1
    ! representation is not exact anyway!
    !
    p_rtriangulation => rdestVector%p_rspatialDiscr%p_rtriangulation
    if (.not. associated(p_rtriangulation)) then
      print *,'ppns2D_bdforces_uniform: Unknown triangulation!'
    end if
    
    ! Get pointers to the subvectors from the block vector
    call lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
    call lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
    call lsyssc_getbase_double (rdestVector,p_Dx)
    
    ! Auxiliary array
    call storage_new1D ('ppns2D_streamfct_uniform', 'aux', &
                        p_rtriangulation%NVT, ST_INT, haux,ST_NEWBLOCK_ZERO)
    call storage_getbase_int (haux,p_Iind)
    
    ! Get stuff from the triangulation
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    
    ! Clear the solution. The auxiliary array is already = 0.
    call lsyssc_clearVector (rdestVector)
    
    ! Start with element 1 and its first vertex.
    ! Set the streamfunction there to 0.0. The streamfunction
    ! in the other vertices are then relatively calculated 
    ! to this basis.
    ! p_Iind is a marker which is set to 1 for every node that 
    ! is finished.

    p_Dx(p_IverticesAtElement(1,1))   = 0.0_DP
    p_Iind(p_IverticesAtElement(1,1)) = 1

    ! Loop over the elements:

    do icurrentelement=1,p_rtriangulation%NEL

      ! We set the current element iel to icurrentelement:

      iel=icurrentelement
      
      ! On the current element, loop over the vertices of
      ! that element. Add the p_Iind-values of all vertices of the
      ! current element together. So, imarkCounter gets the number of marked
      ! vertices and ilastMarked will be the number of the "last" marked
      ! vertex of the four on the element.

      imarkCounter = 0
      do ive = 1,NVE
        jve=p_IverticesAtElement(ive,iel)
        imarkCounter=imarkCounter+p_Iind(jve)
        if (p_Iind(jve) .ge. 1) ilastMarked=ive
      end do

      ! If all four vertices are marked, there's nothing to calculate
      ! on the element. If no vertex is marked, we can't calculate
      ! anything on the current element. In both cases skip the
      ! computation and search for a better element:

      if ((imarkCounter .gt. 0) .and. (imarkCounter .lt. NVE)) then

        ! Ok, here we found an element where some of the vertices are
        ! marked and some not. Here we can calculate a part of the
        ! streamfunction.
        
        call calcSFC_Q1TQ1 (p_DvertexCoords(:,1:p_rtriangulation%NVT),&
            p_IverticesAtElement,p_IedgesAtElement,&
            iel,ilastMarked,p_Iind,p_DdataUX,p_DdataUY,p_Dx)

        ! Now on the current element iel, on all (corner) vertices the
        ! streamfunction is calculated. We go on looking to the adjacent
        ! elements of iel if there's an element where the streamfunction
        ! is not calculated in all vertices...

      end if

      ! Now we make a 'greedy' search from the current element to find
      ! as many other elements as possible where we can calculate the
      ! streamfunction.
      ! Look onto the adjacent elements of the current element if there's
      ! a suitable neighbour element where we can continue the calculation.
      
      neighbourloop: do
      
        ! Loop over the edges of the current element

        do iadj=1,NVE

          ! Get the neighbour element adjacent to the current element

          ielaux=p_IneighboursAtElement(iadj,iel)
          
          if (ielaux.ne.0) then
          
            ! Now we have the number of the neighbour element in ielaux.
            ! Loop about the vertices of the element and sum up the
            ! markers into imarkCounter.
          
            imarkCounter=0
            do ive=1,NVE
              jve=p_IverticesAtElement(ive,ielaux)
              imarkCounter=imarkCounter+p_Iind(jve)
              if (p_Iind(jve) .ge. 1) imarktmp = ive
            end do
            
            ! If there is at least one but not all markers set, the
            ! element can be used for further calculation.
            ! Switch the current element iel to that one and
            ! calculate the streamfunction here.

            if ((imarkCounter .gt. 0) .and. (imarkCounter .lt. NVE)) then

              iel = ielaux
              ilastMarked = imarktmp

              call calcSFC_Q1TQ1 (p_DvertexCoords,p_IverticesAtElement,&
                                  p_IedgesAtElement,&
                                  iel,ilastMarked,p_Iind,p_DdataUX,p_DdataUY,p_Dx)
                                  
              ! Continue the search from here
              cycle neighbourloop
            end if
            
          end if ! ielaux <> 0

        end do ! iadj

        ! We cannot continue from here. Leave the element subset we just marked/
        ! calculated and continue with the standard search to find the next element
        ! that is only partially completed.
        exit neighbourloop
        
      end do neighbourloop
    
    end do ! icurrentelement

    ! At last, normalize the streamfunction such that vertex 1
    ! has value 0.0. Remember that we assigned a value of 0.0
    ! to the first vertex of element 1, which is usually not
    ! vertex 1 of the triangulation!

    call lalg_vectorAddScalarDble (p_Dx,-p_Dx(1))
    
    ! Release temp memory, finish
    call storage_free (haux)

  contains
  
    subroutine calcSFC_Q1TQ1 (DvertexCoords,IverticesAtElement,IedgesAtElement,&
                              iel,ibaseIdx,Imarkers,Du,Dv,Dx)
    
    ! Calculates the value of the streamfunction in all vertices of
    ! element iel. Du(Dv is the X/Y-velocity in Q1~-discretisaton.
    ! Dx is the destination vector and is creatd as Q1-vector.
    
    ! Point coordinates
    real(DP), dimension(:,:), intent(IN)               :: DvertexCoords
    
    ! Vertices at the element
    integer, dimension(:,:), intent(IN) :: IverticesAtElement

    ! Edges at the element
    integer, dimension(:,:), intent(IN) :: IedgesAtElement
    
    ! Element number where to calculate the streamfunction
    integer, intent(IN)               :: iel
    
    ! Local number (1..NVE) of any of the vertices on the element iel
    ! where the streamfunction is already calculated. This will be used
    ! as 'base' to calculate the others.
    integer, intent(IN)                                :: ibaseIdx
    
    ! Marker array of length NVT. All vertices where streamfunction
    ! values are calculated are marked as 1.
    integer, intent(INOUT), dimension(:)          :: Imarkers
    
    ! X/Y-velocity.
    real(DP), dimension(:), intent(IN)                 :: Du, Dv
    
    ! Vector of size NVT; receives the values of the streamfunction
    real(DP), dimension(:), intent(INOUT)              :: Dx
  
    ! local variables
    integer, parameter :: NVE = 4
    integer :: ive,inextidx,inextvertex,imarked
    integer :: ivt,NVT
    integer :: imid
    real(DP) :: dpx1,dpx2,dpy1,dpy2,dn1,dn2
    integer :: ilastMarked
  
      ilastmarked = ibaseIdx
      NVT = ubound(DvertexCoords,2)
  
      ! Loop over the vertices on the element. We can skip the one
      ! where the SF-value is already calculated.
      do ive=1,NVE-1

        ! ibaseIdx is the index of a marked vertex. Calculate the "next"
        ! vertex following ilastMarked and its vertex number into inextvertex.

        inextidx=mod(ilastMarked,NVE)+1
        inextvertex=p_IverticesAtElement(inextidx,iel)

        ! If that vertex is not marked, the streamfunction is not
        ! calculated there. Otherwise we are just looking at two 
        ! marked neighboured vertices, so there's nothing to gain here. 

        if (Imarkers(inextvertex) .eq. 0) then
        
          ! Vertex ivt (corresponding to ilastMarked) is marked, vertex inextvertex 
          ! is not. 
          !
          !     x---x inextvertex
          !     |   |
          !     x---O ilastmarked
          !
          ! Mark vertex inextvertex to indicate that the streamfunction
          ! is not being calculated there:
        
          Imarkers(inextvertex) = 1

          ! and calculate the streamfunction in inextvertex.
          ivt  = p_IverticesAtElement(ilastMarked,iel)
          
          ! imid is now the midpoint number following ivt and thus
          ! the number of the DOF in the FE function.
          ! Calculate the normal vector of the current edge into
          ! N=(dn1,dn2) - not normalized.
          !
          !     x-------x inextvertex
          !     |       |
          !     |  imid x--> (dn1,dn2)
          !     |       |
          !     x-------O ivt

          dpx1=DvertexCoords(1,ivt)
          dpy1=DvertexCoords(2,ivt)
          dpx2=DvertexCoords(1,inextvertex)
          dpy2=DvertexCoords(2,inextvertex)
          dn1 = dpy2-dpy1
          dn2 =-dpx2+dpx1
          
          ! Get the DOF that corresponds to the current edge.
          ! The number of the edge coincides with the number of the
          ! vertex that was marked at last.
          imid = IedgesAtElement (ilastMarked,iel)
          
          ! Calculate the streamfunction in inextvertex from the value
          ! in ivt by:
          !
          ! sfc(inextvertex) = sfc(ivt) + U(imid)*N
          !
          ! which is the "amount of flow over the edge (ivt,inextvertex)".

          Dx(inextvertex)=Dx(ivt)+(Du(imid)*dn1+Dv(imid)*dn2)
        
        end if ! (Imarkers(inextvertex) == 0)

        ! Go on to the next vertex on the element to look if that one
        ! has a not marked neighbour.

        ilastMarked=inextidx
          
      end do ! ive
    
    end subroutine

  end subroutine

end module
