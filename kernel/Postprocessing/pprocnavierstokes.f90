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

MODULE pprocnavierstokes

  USE fsystem
  USE storage
  USE boundary
  USE cubature
  USE triangulation
  USE linearsystemscalar
  USE linearsystemblock
  USE spatialdiscretisation
  USE meshregion

  IMPLICIT NONE

CONTAINS

  !****************************************************************************

!<subroutine>

  SUBROUTINE ppns2D_bdforces_uniform (rvector,rregion,Dforces,ccub,df1,df2)

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
  TYPE(t_vectorBlock), INTENT(IN)    :: rvector
  
  ! Boundary region where to calculate the boundary forces.
  ! Can be created e.g. by boundary_createRegion.
  TYPE(t_boundaryRegion), INTENT(IN) :: rregion
  
  ! 1D Cubature formula identifier to use for the line integration.
  ! One of the CUB_xxxx_1D constants in the cubature.f90.
  INTEGER, INTENT(IN)                 :: ccub

  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  REAL(DP), INTENT(IN), OPTIONAL      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  REAL(DP), INTENT(IN), OPTIONAL      :: df2
  
!</input>

!<output>
  ! Array receiving the forces acting on the boundary specified by rregion.
  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dforces
!</output>

!</subroutine>

  ! Spatial discretisation structure of velocity and pressure
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrU, p_rdiscrP

  ! Element type identifier for U and P
  INTEGER(I32) :: ielemU, ielemP
  
  ! Number of local DOF's in U and P
  INTEGER :: idoflocU, idoflocP
  
  ! Triangulation
  TYPE (t_triangulation), POINTER :: p_rtriangulation
  
  ! An accepting the DOF's of an element.
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS), TARGET :: IdofsU, IdofsP
  
  ! Coordinates of the coordinates of an element
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: DCoords
  
  ! Coordinates of the cubature points on reference and real element
  REAL(DP), DIMENSION(NDIM2D,CUB_MAXCUBP_1D) :: DpointsRef,DpointsReal
  
  ! Coordinate system for U and P element
  INTEGER :: ctrafoU, ctrafoP
  
  ! U/P element parametric or nonparametric
  LOGICAL :: bnonparU,bnonparP
  
  ! Arrays for saving Jacobian determinants and matrices
  REAL(DP), DIMENSION(CUB_MAXCUBP_1D) :: Ddetj
  REAL(DP), DIMENSION(EL_NJACENTRIES2D,CUB_MAXCUBP_1D) :: Djac

  ! Array to tell the element which derivatives to calculate.
  LOGICAL, DIMENSION(EL_MAXNDER) :: BderU, BderP
  
  ! Value of basis functions
  REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,CUB_MAXCUBP_1D) :: DbasU, DbasP
  
  ! Pointer to vector data of solution vector
  REAL(DP), DIMENSION(:), POINTER :: p_DdataUX,p_DdataUY,p_DdataP

  ! Cubature point coordinates on the reference element.
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi1D, DXi2D

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! number/index of cubature points on the reference element
  INTEGER :: ncubp,icubp

  ! Edges, vertices and elements on the boundary
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtBoundary
  INTEGER(PREC_EDGEIDX), DIMENSION(:), POINTER    :: p_IedgesAtBoundary
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER  :: p_IedgesAtElement
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
  REAL(DP), DIMENSION(:), POINTER                 :: p_DvertexParameterValue
  INTEGER(I32), DIMENSION(:), POINTER             :: p_IboundaryCpIdx
  REAL(DP), DIMENSION(:,:), POINTER               :: p_DvertexCoordinates

  ! other local variables
  INTEGER :: iedgeidx,ivt1,ivt2,ilocaledge,nlocaledges,idfl,icp
  INTEGER(PREC_DOFIDX) :: neqU,neqP
  INTEGER(PREC_EDGEIDX) :: iedge,iedgeglobal
  INTEGER(PREC_ELEMENTIDX) :: iel
  REAL(DP) :: dvtp1,dvtp2,dedgelen, dweight,dut,dpf1,dpf2
  REAL(DP), DIMENSION(2) :: DintU, DintP
  REAL(DP) :: dpres
  REAL(DP), DIMENSION(NDIM2D) :: dvt1,dvt2,dtangential,dnormal
  INTEGER(I32), DIMENSION(:), POINTER :: p_ItwistIndex

    ! Get the vector data
    neqU = rvector%RvectorBlock(1)%NEQ
    neqP = rvector%RvectorBlock(3)%NEQ
    
    IF (rvector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'ppns2D_bdforces: Unsupported vector precision.'
      CALL sys_halt()
    END IF

    ! We support only uniform discretisation structures.
    IF (.NOT. ASSOCIATED(rvector%p_rblockDiscr)) THEN
      PRINT *,'ppns2D_bdforces: No discretisation structure!'
      CALL sys_halt()
    END IF

    IF (rvector%p_rblockDiscr%ccomplexity .NE. SPDISC_UNIFORM) THEN
      PRINT *,'ppns2D_bdforces_uniform: Discretisation too complex!'
      CALL sys_halt()
    END IF
    
    ! Get pointers to the subvectors from the block vector
    CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
    CALL lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
    CALL lsyssc_getbase_double (rvector%RvectorBlock(3),p_DdataP)
    
    IF ((rvector%RvectorBlock(1)%isortStrategy > 0) .OR. &
        (rvector%RvectorBlock(2)%isortStrategy > 0) .OR. &
        (rvector%RvectorBlock(3)%isortStrategy > 0)) THEN
      PRINT *,'ppns2D_bdforces_uniform: Resorted vectors not supported!'
      CALL sys_halt()
    END IF
    
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
    
    CALL cub_getCubPoints (ccub,ncubp,Dxi1D,Domega)
    
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
    ALLOCATE(DCoords(NDIM2D,MAX(elem_igetNVE(ielemU),elem_igetNVE(ielemP))))
    
    ! The triangulation - it's the same for U and P
    p_rtriangulation => p_rdiscrU%p_rtriangulation

    ! The arrays which contain the elements and edge numbers on the boundary
    ! as well as boundary index array.
    ! Fetch geometrical information
    CALL storage_getbase_int (p_rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)
    CALL storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement, &
        p_IverticesAtElement)
    CALL storage_getbase_int (p_rtriangulation%h_IelementsAtBoundary, &
        p_IelementsAtBoundary)
    CALL storage_getbase_int (p_rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    CALL storage_getbase_double2d (p_rtriangulation%h_DvertexCoords, &
                                   p_DvertexCoordinates)
                                   
    IF (p_rtriangulation%h_DvertexParameterValue .EQ. ST_NOHANDLE) THEN
      PRINT *,'No boundary parameters available!'
      CALL sys_halt()
    END IF
    
    CALL storage_getbase_double (p_rtriangulation%h_DvertexParameterValue, &
                                 p_DvertexParameterValue)
                                 
    ! Does the element need twist indices?
    NULLIFY(p_ItwistIndex)
    IF (p_rtriangulation%h_ItwistIndexEdges .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int (p_rtriangulation%h_ItwistIndexEdges,p_ItwistIndex)
    END IF

    ! Is one of the elements nonparametric
    bnonparU = elem_isnonparametric(ielemU) 
    bnonparP = elem_isnonparametric(ielemP)

    ! Coordinate systems of U and P element
    ctrafoU = elem_igetTrafoType(ielemU)
    ctrafoP = elem_igetTrafoType(ielemP)
    
    ! Derivatives to calculate when evaluating the U and P-element, respectively.
    BderU = .FALSE.
    BderU(DER_DERIV_X) = .TRUE.
    BderU(DER_DERIV_Y) = .TRUE.
    
    BderP = .FALSE.
    BderP(DER_FUNC) = .TRUE.
    
    ! Prepare the weighting coefficients
    dpf1 = 1.0_DP
    dpf2 = 2.0_DP
    IF (PRESENT(df1)) dpf1 = df1
    IF (PRESENT(df2)) dpf2 = df2

    ! We are on boundary component
    icp = rregion%iboundCompIdx
    
    ! We assemble the integral contributions separately
    DintU = 0.0_DP
    DintP = 0.0_DP
    
    ! Loop through all edges along the boundary
    DO iedgeidx = p_IboundaryCpIdx(icp),p_IboundaryCpIdx(icp+1)-1
      
      ! Current element
      iel = p_IelementsAtBoundary(iedgeidx)
      
      ! Egde number
      iedgeglobal = p_IedgesAtBoundary(iedgeidx)
      iedge       = iedgeglobal - p_rtriangulation%NVT
      
      ! What are the vertices adjacent to that edge?
      ivt1 = p_IverticesAtEdge(1,iedge)
      ivt2 = p_IverticesAtEdge(2,iedge)
      
      ! Parameter values of these vertices?
      ! iedgeidx is the index of the edge as well as the number of the
      ! index of the vertex preceding the edge!
      dvtp1 = p_DvertexParameterValue (iedgeidx)
      IF (iedgeidx .NE. p_IboundaryCpIdx(icp+1)-1) THEN
        dvtp2 = p_DvertexParameterValue (iedgeidx+1)
      ELSE
        ! Last vertex has maximum parameter value! (TMAX)
        dvtp2 = boundary_dgetMaxParVal(p_rdiscrU%p_rboundary,icp)
      END IF
      
      ! Is the edge in the specified boundary region, so we are allowed 
      ! to integrate about it? We check both endpoints...
      
      IF (boundary_isInRegion (rregion,icp,dvtp1) .AND. &
          boundary_isInRegion (rregion,icp,dvtp2)) THEN
          
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
        
        DO ilocaledge = 1,nlocaledges
          IF (p_IedgesAtElement(ilocaledge,iel) .EQ. iedgeglobal) EXIT
        END DO
        
        IF (ilocaledge .GT. nlocaledges) THEN
          PRINT *,'ppns2D_bdforces: Edge not found. KMID destroyed?'
          CALL sys_halt()
        END IF
        
        ! The number of the edge is in ilocaledge. We have to transfer
        ! the coordinates of the cubature points from 1D to 2D depending
        ! on this edge.
        
        CALL trafo_mapCubPts1Dto2DRefQuad(ilocaledge, ncubp, Dxi1D, Dxi2D)

        ! Now, in Dxi2D we have all the cubature points on the 2D reference
        ! element along the edge!
        ! Map the coordinates in a proper 2D array
        DpointsRef (1:NDIM2D,1:ncubp) = TRANSPOSE(Dxi2D(1:ncubp,1:NDIM2D))
                
        ! For the integration, we need the global DOF's on our element
        ! for U and P:
        
        CALL dof_locGlobMapping(p_rdiscrU, iel, IdofsU)
        CALL dof_locGlobMapping(p_rdiscrP, iel, IdofsP)
        
        ! Get the coordinates of the point on the current element
        Dcoords (1:NDIM2D,1:nlocaledges) = &
            p_DvertexCoordinates(1:NDIM2D, p_IverticesAtElement (1:nlocaledges,iel))
        
        ! Calculate the transformation for all points on the current element.
        ! If we have only parametric elements, we don't have to calculate
        ! the real coordinates of the points.
        IF (bnonparU .OR. bnonparP) THEN
          CALL trafo_calctrafo_mult (ctrafoU,ncubp,Dcoords,&
                                     DpointsRef,Djac,Ddetj,DpointsReal)
        ELSE
          CALL trafo_calctrafo_mult (ctrafoU,ncubp,Dcoords,&
                                     DpointsRef,Djac,Ddetj)
        END IF
        
        ! Evaluate the U- and P-element in all our cubature points
        IF (bnonparU) THEN
          CALL elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsReal,p_ItwistIndex(iel))
        ELSE
          CALL elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsRef,p_ItwistIndex(iel))
        END IF

        IF (bnonparP) THEN
          CALL elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsReal,p_ItwistIndex(iel))
        ELSE
          CALL elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsRef,p_ItwistIndex(iel))
        END IF
        
        ! Loop over the cubature points on the current element
        ! to assemble the integral
        DO icubp = 1,ncubp
        
          ! Calculate the OMEGA for the integration by multiplication
          ! of the integration coefficient by the Jacobian of the
          ! mapping.
          ! The determinant of the mapping of the unit interval [-1,1]
          ! to the real line is 0.5*length of the line!
            
          dweight = Domega(icubp)*0.5_DP*dedgelen
          
          ! Loop through the DOF's on our element and calculate
          ! the tangential U as well as P.
          dut = 0.0_DP
          DO idfl=1,idoflocU
            dut = dut &
            + p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)*Dtangential(1) &
                                      * Dnormal(1) &
            + p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)*Dtangential(2) &
                                      * Dnormal(1) &
            + p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)*Dtangential(1) &
                                      * Dnormal(2) &
            + p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)*Dtangential(2) &
                                      * Dnormal(2)
          END DO

          dpres = 0.0_DP
          DO idfl=1,idoflocP
            dpres = dpres+p_DdataP(IdofsP(idfl))*DbasP(idfl,DER_FUNC,icubp)
          END DO
          
          ! Sum this up to the two integral contributions for the pressure and
          ! velocity.
          DintU(1) = DintU(1) + dweight * dut * Dnormal(2)
          DintU(2) = DintU(2) - dweight * dut * Dnormal(1)
          DintP(1) = DintP(1) - dweight * dpres * Dnormal(1)
          DintP(2) = DintP(2) - dweight * dpres * Dnormal(2)
        
        END DO ! icubp
          
      END IF
      
    END DO ! iedgeidx
    
    ! DintU and DintP give now the contributions to the force integral:
    !
    ! DragCoeff = 2/dfp2 * (dpf1*DintU(1) + DintP(1))
    ! LiftCoeff = 2/dfp2 * (dpf1*DintU(2) + DintP(2))
    
    Dforces = 0.0_DP
    Dforces(1:NDIM2D) = 2.0_DP/dpf2 * (dpf1*DintU(:) + DintP(:))
    
    ! Deallocate memory, finish.
    DEALLOCATE(DCoords)
  
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE ppns3D_bdforces_uniform (rvector,rregion,Dforces,ccub,df1,df2)

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
  TYPE(t_vectorBlock), INTENT(IN)     :: rvector
  
  ! Mesh region where to calculate the boundary forces.
  TYPE(t_meshRegion), INTENT(IN)      :: rregion
  
  ! 2D Cubature formula identifier to use for the quad integration.
  ! One of the CUB_xxxx constants in the cubature.f90.
  INTEGER, INTENT(IN)                 :: ccub

  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  REAL(DP), INTENT(IN), OPTIONAL      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  REAL(DP), INTENT(IN), OPTIONAL      :: df2
  
!</input>

!<output>
  ! Array receiving the forces acting on the boundary specified by rregion.
  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dforces
!</output>

!</subroutine>

  ! Spatial discretisation structure of velocity and pressure
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrU, p_rdiscrP

  ! Element type identifier for U and P
  INTEGER(I32) :: ielemU, ielemP
  
  ! Number of local DOF's in U and P
  INTEGER :: idoflocU, idoflocP
  
  ! Triangulation
  TYPE (t_triangulation), POINTER :: p_rtria
  
  ! An accepting the DOF's of an element.
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS), TARGET :: IdofsU, IdofsP
  
  ! Coordinates of the vertices
  REAL(DP), DIMENSION(NDIM3D,8) :: Dcoords
  
  ! Coordinates of the normal vectors
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Dnormal

  ! Coordinates of the cubature points on reference and real element
  REAL(DP), DIMENSION(NDIM3D,CUB_MAXCUBP_2D) :: DpointsRef,DpointsReal
  
  ! Coordinate system for U and P element
  INTEGER :: ctrafoU, ctrafoP
  
  ! U/P element parametric or nonparametric
  LOGICAL :: bnonparU,bnonparP
  
  ! Arrays for saving Jacobian determinants and matrices
  REAL(DP), DIMENSION(CUB_MAXCUBP_2D) :: Ddetj
  REAL(DP), DIMENSION(CUB_MAXCUBP_2D) :: Ddetj_face
  REAL(DP), DIMENSION(EL_NJACENTRIES3D,CUB_MAXCUBP_2D) :: Djac

  ! Array to tell the element which derivatives to calculate.
  LOGICAL, DIMENSION(EL_MAXNDER) :: BderU, BderP
  
  ! Value of basis functions
  REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER,CUB_MAXCUBP_2D) :: DbasU, DbasP
  
  ! Pointer to vector data of solution vector
  REAL(DP), DIMENSION(:), POINTER :: p_DdataUX,p_DdataUY,p_DdataUZ,p_DdataP

  ! Cubature point coordinates on the reference element.
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi2D, Dxi3D

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! number/index of cubature points on the reference element
  INTEGER :: ncubp,icubp

  ! Edges, vertices and elements on the boundary
  INTEGER(PREC_FACEIDX), DIMENSION(:), POINTER     :: p_IfaceIdx
  INTEGER(PREC_FACEIDX), DIMENSION(:,:), POINTER   :: p_IfaceAtElem
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IvertAtElem
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IelemAtFace
  REAL(DP), DIMENSION(:,:), POINTER                :: p_Dvertex

  ! other local variables
  INTEGER :: iat, iface,ilocalface,idfl,i,NVT,NMT
  INTEGER(PREC_DOFIDX) :: neqU,neqP
  INTEGER(PREC_ELEMENTIDX) :: iel
  REAL(DP) :: dpf1,dpf2,dpres,du1,du2,du3,dweight,dv
  REAL(DP), DIMENSION(3) :: DintU, DintP
  INTEGER(I32), DIMENSION(:), POINTER :: p_ItwistIndex

    ! Get the vector data
    neqU = rvector%RvectorBlock(1)%NEQ
    neqP = rvector%RvectorBlock(4)%NEQ
    
    IF (rvector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'ppns3D_bdforces: Unsupported vector precision.'
      CALL sys_halt()
    END IF

    ! We support only uniform discretisation structures.
    IF (.NOT. ASSOCIATED(rvector%p_rblockDiscr)) THEN
      PRINT *,'ppns3D_bdforces: No discretisation structure!'
      CALL sys_halt()
    END IF

    IF (rvector%p_rblockDiscr%ccomplexity .NE. SPDISC_UNIFORM) THEN
      PRINT *,'ppns3D_bdforces_uniform: Discretisation too complex!'
      CALL sys_halt()
    END IF
    
    ! Get pointers to the subvectors from the block vector
    CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
    CALL lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
    CALL lsyssc_getbase_double (rvector%RvectorBlock(3),p_DdataUZ)
    CALL lsyssc_getbase_double (rvector%RvectorBlock(4),p_DdataP)
    
    IF ((rvector%RvectorBlock(1)%isortStrategy > 0) .OR. &
        (rvector%RvectorBlock(2)%isortStrategy > 0) .OR. &
        (rvector%RvectorBlock(3)%isortStrategy > 0) .OR. &
        (rvector%RvectorBlock(4)%isortStrategy > 0)) THEN
      PRINT *,'ppns3D_bdforces_uniform: Resorted vectors not supported!'
      CALL sys_halt()
    END IF
    
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
    CALL cub_getCubPoints (ccub,ncubp,Dxi2D,Domega)
    
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
    CALL storage_getbase_int2D(p_rtria%h_IelementsAtFace, p_IelemAtFace)
       
    ! Get a pointer to the faces-at-element array
    CALL storage_getbase_int2D(p_rtria%h_IfacesAtElement, p_IfaceAtElem)

    ! Get a pointer to the vertices-at-element array
    CALL storage_getbase_int2D(p_rtria%h_IverticesAtElement, p_IvertAtElem)
    
    ! Get the vertice coordinate array
    CALL storage_getbase_double2D(p_rtria%h_DvertexCoords, p_Dvertex)
    
    ! And get the face index array of the mesh region
    CALL storage_getbase_int(rregion%h_IfaceIdx, p_IfaceIdx)
    
    NVT = p_rtria%NVT
    NMT = p_rtria%NMT
                                   
    ! Does the element need twist indices?
    NULLIFY(p_ItwistIndex)
    IF (p_rtria%h_ItwistIndexEdges .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int (p_rtria%h_ItwistIndexFaces,p_ItwistIndex)
    END IF

    ! Is one of the elements nonparametric
    bnonparU = elem_isnonparametric(ielemU) 
    bnonparP = elem_isnonparametric(ielemP)

    ! Coordinate systems of U and P element
    ctrafoU = elem_igetTrafoType(ielemU)
    ctrafoP = elem_igetTrafoType(ielemP)
    
    ! Derivatives to calculate when evaluating the U and P-element, respectively.
    BderU = .FALSE.
    BderU(DER_DERIV3D_X) = .TRUE.
    BderU(DER_DERIV3D_Y) = .TRUE.
    BderU(DER_DERIV3D_Z) = .TRUE.
    
    BderP = .FALSE.
    BderP(DER_FUNC3D) = .TRUE.
    
    ! Prepare the weighting coefficients
    dpf1 = 1.0_DP
    dpf2 = 2.0_DP
    IF (PRESENT(df1)) dpf1 = df1
    IF (PRESENT(df2)) dpf2 = df2

    ! We assemble the integral contributions separately
    DintU = 0.0_DP
    DintP = 0.0_DP
    
    ! Calculate the normal vectors on the faces
    ALLOCATE(Dnormal(3, rregion%NAT))
    CALL mshreg_calcBoundaryNormals3D(rregion, Dnormal)
    
    ! Loop through all faces along the boundary
    DO iat = 1, rregion%NAT
    
      ! Get the index fo the face
      iface = p_IfaceIdx(iat)
      
      ! Current element
      iel = p_IelemAtFace(1,iface)
      
      ! Get the coordinates of the corner vertices on the current element
      DO i=1, 8
        Dcoords (:,i) = p_Dvertex(:, p_IvertAtElem (i,iel))
      END DO

      ! Get the local index of the face on that element
      DO ilocalface = 1, 6
        IF(p_IfaceAtElem(ilocalface,iel) .EQ. iface+NVT+NMT) EXIT
      END DO
      
      ! We have to transfer the coordinates of the cubature points from
      ! 2D to 3D depending on this localface.
      CALL trafo_mapCubPts2Dto3DRefHexa(ilocalface, ncubp, Dxi2D, Dxi3D)
      
      ! And calculate the determinants for the mapping
      CALL ppns3D_det(Dcoords,ilocalface,Dxi2D,ncubp,Ddetj_face)

      ! Now, in Dxi3D we have all the cubature points on the 3D reference
      ! element along the face!
      ! Map the coordinates in a proper 3D array
      DpointsRef (1:NDIM3D,1:ncubp) = TRANSPOSE(Dxi3D(1:ncubp,1:NDIM3D))
              
      ! For the integration, we need the global DOF's on our element
      ! for U and P:
      CALL dof_locGlobMapping(p_rdiscrU, iel, IdofsU)
      CALL dof_locGlobMapping(p_rdiscrP, iel, IdofsP)
      
      ! Calculate the transformation for all points on the current element.
      ! If we have only parametric elements, we don't have to calculate
      ! the real coordinates of the points.
      IF (bnonparU .OR. bnonparP) THEN
        CALL trafo_calctrafo_mult (ctrafoU,ncubp,Dcoords,&
                                   DpointsRef,Djac,Ddetj,DpointsReal)
      ELSE
        CALL trafo_calctrafo_mult (ctrafoU,ncubp,Dcoords,&
                                   DpointsRef,Djac,Ddetj)
      END IF
      
      IF(ASSOCIATED(p_ItwistIndex)) THEN
      
        IF (bnonparU) THEN
          CALL elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsReal,p_ItwistIndex(iel))
        ELSE
          CALL elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsRef,p_ItwistIndex(iel))
        END IF

        IF (bnonparP) THEN
          CALL elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsReal,p_ItwistIndex(iel))
        ELSE
          CALL elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsRef,p_ItwistIndex(iel))
        END IF
        
      ELSE
      
        IF (bnonparU) THEN
          CALL elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsReal)
        ELSE
          CALL elem_generic_mult (ielemU, Dcoords, Djac, Ddetj, &
                                  BderU, DbasU, ncubp, DpointsRef)
        END IF

        IF (bnonparP) THEN
          CALL elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsReal)
        ELSE
          CALL elem_generic_mult (ielemP, Dcoords, Djac, Ddetj, &
                                  BderP, DbasP, ncubp, DpointsRef)
        END IF
      
      END IF
      
      ! Loop over the cubature points on the current element
      ! to assemble the integral
      DO icubp = 1,ncubp
      
        ! Calculate the OMEGA for the integration by multiplication
        ! of the integration coefficient by the Jacobian of the
        ! mapping.
        dweight = Domega(icubp)*Ddetj_face(icubp)
        
        ! Loop through the DOF's on our element and calculate
        ! U as well as P.
        du1 = 0.0_DP
        du2 = 0.0_DP
        du3 = 0.0_DP
        DO idfl=1,idoflocU
        
           dv = DbasU(idfl,DER_DERIV3D_X,icubp)*Dnormal(1,iat)&
              + DbasU(idfl,DER_DERIV3D_Y,icubp)*Dnormal(2,iat)&
              + DbasU(idfl,DER_DERIV3D_Z,icubp)*Dnormal(3,iat)
              
           du1 = du1 + p_DdataUX(IdofsU(idfl)) * dv
           du2 = du2 + p_DdataUY(IdofsU(idfl)) * dv
           du3 = du3 + p_DdataUZ(IdofsU(idfl)) * dv
               
        END DO

        dpres = 0.0_DP
        DO idfl=1,idoflocP
          dpres = dpres + p_DdataP(IdofsP(idfl))*DbasP(idfl,DER_FUNC3D,icubp)
        END DO
        
        ! Sum this up to the two integral contributions for the pressure and
        ! velocity.
        DintU(1) = DintU(1) + dweight * du1
        DintU(2) = DintU(2) + dweight * du2
        DintU(3) = DintU(3) + dweight * du3
        DintP(1) = DintP(1) - dweight * dpres * Dnormal(1,iat)
        DintP(2) = DintP(2) - dweight * dpres * Dnormal(2,iat)
        DintP(3) = DintP(3) - dweight * dpres * Dnormal(3,iat)
      
      END DO ! icubp
          
    END DO ! iat
    
    ! DintU and DintP give now the contributions to the force integral:
    !
    ! DragCoeff = 2/dfp2 * (dpf1*DintU(1) + DintP(1))
    ! LiftCoeff = 2/dfp2 * (dpf1*DintU(2) + DintP(2))
    Dforces = 0.0_DP
    Dforces(:) = 2.0_DP/dpf2 * (dpf1*DintU(:) + DintP(:))
    
    ! Deallocate memory, finish.
    DEALLOCATE(Dnormal)
    
    ! That's it
    
    CONTAINS
    
!<subroutine>

    PURE SUBROUTINE ppns3D_det(Dcoords,iface,Dcub,ncubp,Ddetj)

!<describtion>
  ! This routine calculates the "determinant" of a bilinear quadrilateral
  ! transformation from the 2D reference quadrilateral onto a hexahedron
  ! face in 3D.
!</describtion>
    
!<input>
    ! The coordinates of the eight corner vertices of the hexahedron
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
    
    ! The index of the face onto which the points are mapped
    INTEGER, INTENT(IN) :: iface
    
    ! The 2D coordinates of the points that are to be mapped
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcub
    
    ! The number of points which are to be mapped
    INTEGER, INTENT(IN) :: ncubp
!</input>

!<output>
    ! The jacobian determinants of the mapping
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Ddetj
!</output>

!</subroutine>
    
    ! Local variables
    REAL(DP), DIMENSION(3,4) :: Dv
    REAL(DP), DIMENSION(3,3) :: Dt
    REAL(DP), DIMENSION(3) :: Dx,Dy,Dn
    INTEGER :: i,j
    
      ! Onto which face do we map the points?
      ! Note: The orientation of the vertices corresponds to the mapping
      ! routine trafo_mapCubPts2Dto3DRefHexa defined in transformation.f90.
      SELECT CASE(iface)
      CASE (1)
        Dv(:,1) = Dcoords(:,1)
        Dv(:,2) = Dcoords(:,2)
        Dv(:,3) = Dcoords(:,3)
        Dv(:,4) = Dcoords(:,4)
      CASE (2)
        Dv(:,1) = Dcoords(:,1)
        Dv(:,2) = Dcoords(:,2)
        Dv(:,3) = Dcoords(:,6)
        Dv(:,4) = Dcoords(:,5)
      CASE (3)
        Dv(:,1) = Dcoords(:,2)
        Dv(:,2) = Dcoords(:,3)
        Dv(:,3) = Dcoords(:,7)
        Dv(:,4) = Dcoords(:,6)
      CASE (4)
        Dv(:,1) = Dcoords(:,3)
        Dv(:,2) = Dcoords(:,4)
        Dv(:,3) = Dcoords(:,8)
        Dv(:,4) = Dcoords(:,7)
      CASE (5)
        Dv(:,1) = Dcoords(:,4)
        Dv(:,2) = Dcoords(:,1)
        Dv(:,3) = Dcoords(:,5)
        Dv(:,4) = Dcoords(:,8)
      CASE (6)
        Dv(:,1) = Dcoords(:,5)
        Dv(:,2) = Dcoords(:,6)
        Dv(:,3) = Dcoords(:,7)
        Dv(:,4) = Dcoords(:,8)
      END SELECT
      
      ! We have a bilinear mapping T: R^2 -> R^3, so the jacobian matrix
      ! of T is a 3x2 matrix. To get a useful replacement for the determinant
      ! we set the determinant to ||(dT/dx) X (dT/dy)||_2, where 'X' denotes
      ! the 3D cross-product.
      
      ! Calculate transformation coefficients for the jacobian matrix
      DO i = 1,3
        Dt(i,1) = 0.25_DP * (-Dv(i,1) + Dv(i,2) + Dv(i,3) - Dv(i,4))
        Dt(i,2) = 0.25_DP * (-Dv(i,1) - Dv(i,2) + Dv(i,3) + Dv(i,4))
        Dt(i,3) = 0.25_DP * ( Dv(i,1) - Dv(i,2) + Dv(i,3) - Dv(i,4))
      END DO
      
      ! And calculate the determinants
      DO i = 1, ncubp
      
        DO j = 1, 3
          ! Dx := dT / dx
          Dx(j) = Dt(j,1) + Dt(j,3)*Dcub(i,1)
          ! Dy := dT / dy
          Dy(j) = Dt(j,2) + Dt(j,3)*Dcub(i,2)
        END DO
        
        ! Dn := Dx x Dy
        Dn(1) = Dx(2)*Dy(3) - Dx(3)*Dy(2)
        Dn(2) = Dx(3)*Dy(1) - Dx(1)*Dy(3)
        Dn(3) = Dx(1)*Dy(2) - Dx(2)*Dy(1)
        
        ! detj := ||Dn||_2
        Ddetj(i) = SQRT(Dn(1)**2 + Dn(2)**2 + Dn(3)**2)
        
      END DO
    
    END SUBROUTINE
  
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE ppns2D_streamfct_uniform (rvector,rdestVector)

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
  TYPE(t_vectorBlock), INTENT(IN)    :: rvector
!</input>

!<inputoutput>
  ! An empty vector that is prepared with a discretisation structure to
  ! represent the streamfunction. The values in the vector are overwritten
  ! with the FE representation of the streamfunction.
  TYPE(t_vectorScalar), INTENT(INOUT),TARGET    :: rdestVector
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER, PARAMETER :: NVE = 4

    INTEGER(PREC_ELEMENTIDX) :: iel,ielaux,icurrentelement
    INTEGER(PREC_VERTEXIDX) :: jve
    INTEGER(I32) :: ieltype1,ieltype2,ieltypeDest
    INTEGER :: haux,ive,iadj
    INTEGER :: ilastMarked,imarkCounter,imarktmp
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! Pointer to vector data of solution vector
    REAL(DP), DIMENSION(:), POINTER :: p_DdataUX,p_DdataUY,p_Dx
    INTEGER(I32), DIMENSION(:), POINTER :: p_Iind
    
    ! Stuff from the triangulation
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords

    IF (rvector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'ppns2D_streamfct_uniform: Unsupported vector precision.'
      CALL sys_halt()
    END IF

    ! We support only uniform discretisation structures.
    IF (.NOT. ASSOCIATED(rvector%p_rblockDiscr)) THEN
      PRINT *,'ppns2D_streamfct_uniform: No discretisation structure in rvector!'
      CALL sys_halt()
    END IF

    IF (rvector%p_rblockDiscr%ccomplexity .NE. SPDISC_UNIFORM) THEN
      PRINT *,'ppns2D_streamfct_uniform: Discretisation of rvector too complex!'
      CALL sys_halt()
    END IF

    IF (.NOT. ASSOCIATED(rdestVector%p_rspatialDiscr)) THEN
      PRINT *,'ppns2D_streamfct_uniform: No discretisation structure in rdestVector!'
      CALL sys_halt()
    END IF

    IF (rdestVector%p_rspatialDiscr%ccomplexity .NE. SPDISC_UNIFORM) THEN
      PRINT *,'ppns2D_streamfct_uniform: Discretisation of rdestVector too complex!'
      CALL sys_halt()
    END IF

    ieltype1 = rvector%p_rblockDiscr% &
        RspatialDiscr(1)%RelementDistr(1)%celement
    ieltype2 = rvector%p_rblockDiscr% &
        RspatialDiscr(2)%RelementDistr(1)%celement
    ieltypeDest = rdestVector%p_rspatialDiscr% &
        RelementDistr(1)%celement

    IF (elem_getPrimaryElement(ieltype1) .NE. EL_Q1T) THEN
      PRINT *,'ppns2D_streamfct_uniform: rvector must be discretised with Q1~!'
    END IF

    IF (elem_getPrimaryElement(ieltype2) .NE. EL_Q1T) THEN
      PRINT *,'ppns2D_streamfct_uniform: rvector must be discretised with Q1~!'
    END IF

    IF (ieltypeDest .NE. EL_Q1) THEN
      PRINT *,'ppns2D_streamfct_uniform: rdestVector must be discretised with Q1!'
    END IF
    
    IF ((rvector%RvectorBlock(1)%isortStrategy > 0) .OR. &
        (rvector%RvectorBlock(2)%isortStrategy > 0) .OR. &
        (rdestVector%isortStrategy > 0)) THEN
      PRINT *,'ppns2D_bdforces_uniform: Resorted vectors not supported!'
      CALL sys_halt()
    END IF

    ! Let's go. Note that we perform the same computation for all,
    ! parametric and nonparametric, point-based and integral-mean-value
    ! based Q1~ solutions. Gives a slight error, but the final Q1
    ! representation is not exact anyway!
    !
    p_rtriangulation => rdestVector%p_rspatialDiscr%p_rtriangulation
    IF (.NOT. ASSOCIATED(p_rtriangulation)) THEN
      PRINT *,'ppns2D_bdforces_uniform: Unknown triangulation!'
    END IF
    
    ! Get pointers to the subvectors from the block vector
    CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
    CALL lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
    CALL lsyssc_getbase_double (rdestVector,p_Dx)
    
    ! Auxiliary array
    CALL storage_new1D ('ppns2D_streamfct_uniform', 'aux', &
                        p_rtriangulation%NVT, ST_INT, haux,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int (haux,p_Iind)
    
    ! Get stuff from the triangulation
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    CALL storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    
    ! Clear the solution. The auxiliary array is already = 0.
    CALL lsyssc_clearVector (rdestVector)
    
    ! Start with element 1 and its first vertex.
    ! Set the streamfunction there to 0.0. The streamfunction
    ! in the other vertices are then relatively calculated 
    ! to this basis.
    ! p_Iind is a marker which is set to 1 for every node that 
    ! is finished.

    p_Dx(p_IverticesAtElement(1,1))   = 0.0_DP
    p_Iind(p_IverticesAtElement(1,1)) = 1

    ! Loop over the elements:

    DO icurrentelement=1,p_rtriangulation%NEL

      ! We set the current element iel to icurrentelement:

      iel=icurrentelement
      
      ! On the current element, loop over the vertices of
      ! that element. Add the p_Iind-values of all vertices of the
      ! current element together. So, imarkCounter gets the number of marked
      ! vertices and ilastMarked will be the number of the "last" marked
      ! vertex of the four on the element.

      imarkCounter = 0
      DO ive = 1,NVE
        jve=p_IverticesAtElement(ive,iel)
        imarkCounter=imarkCounter+p_Iind(jve)
        IF (p_Iind(jve) .GE. 1) ilastMarked=ive
      END DO

      ! If all four vertices are marked, there's nothing to calculate
      ! on the element. If no vertex is marked, we can't calculate
      ! anything on the current element. In both cases skip the
      ! computation and search for a better element:

      IF ((imarkCounter .GT. 0) .AND. (imarkCounter .LT. NVE)) THEN

        ! Ok, here we found an element where some of the vertices are
        ! marked and some not. Here we can calculate a part of the
        ! streamfunction.
        
        CALL calcSFC_Q1TQ1 (p_DvertexCoords(:,1:p_rtriangulation%NVT),&
            p_IverticesAtElement,p_IedgesAtElement,&
            iel,ilastMarked,p_Iind,p_DdataUX,p_DdataUY,p_Dx)

        ! Now on the current element iel, on all (corner) vertices the
        ! streamfunction is calculated. We go on looking to the adjacent
        ! elements of iel if there's an element where the streamfunction
        ! is not calculated in all vertices...

      END IF

      ! Now we make a 'greedy' search from the current element to find
      ! as many other elements as possible where we can calculate the
      ! streamfunction.
      ! Look onto the adjacent elements of the current element if there's
      ! a suitable neighbour element where we can continue the calculation.
      
      neighbourloop: DO
      
        ! Loop over the edges of the current element

        DO iadj=1,NVE

          ! Get the neighbour element adjacent to the current element

          ielaux=p_IneighboursAtElement(iadj,iel)
          
          IF (ielaux.NE.0) THEN
          
            ! Now we have the number of the neighbour element in ielaux.
            ! Loop about the vertices of the element and sum up the
            ! markers into imarkCounter.
          
            imarkCounter=0
            DO ive=1,NVE
              jve=p_IverticesAtElement(ive,ielaux)
              imarkCounter=imarkCounter+p_Iind(jve)
              IF (p_Iind(jve) .GE. 1) imarktmp = ive
            END DO
            
            ! If there is at least one but not all markers set, the
            ! element can be used for further calculation.
            ! Switch the current element iel to that one and
            ! calculate the streamfunction here.

            IF ((imarkCounter .GT. 0) .AND. (imarkCounter .LT. NVE)) THEN

              iel = ielaux
              ilastMarked = imarktmp

              CALL calcSFC_Q1TQ1 (p_DvertexCoords,p_IverticesAtElement,&
                                  p_IedgesAtElement,&
                                  iel,ilastMarked,p_Iind,p_DdataUX,p_DdataUY,p_Dx)
                                  
              ! Continue the search from here
              CYCLE neighbourloop
            END IF
            
          END IF ! ielaux <> 0

        END DO ! iadj

        ! We cannot continue from here. Leave the element subset we just marked/
        ! calculated and continue with the standard search to find the next element
        ! that is only partially completed.
        EXIT neighbourloop
        
      END DO neighbourloop
    
    END DO ! icurrentelement

    ! At last, normalize the streamfunction such that vertex 1
    ! has value 0.0. Remember that we assigned a value of 0.0
    ! to the first vertex of element 1, which is usually not
    ! vertex 1 of the triangulation!

    CALL lalg_vectorAddScalarDble (p_Dx,-p_Dx(1))
    
    ! Release temp memory, finish
    CALL storage_free (haux)

  CONTAINS
  
    SUBROUTINE calcSFC_Q1TQ1 (DvertexCoords,IverticesAtElement,IedgesAtElement,&
                              iel,ibaseIdx,Imarkers,Du,Dv,Dx)
    
    ! Calculates the value of the streamfunction in all vertices of
    ! element iel. Du(Dv is the X/Y-velocity in Q1~-discretisaton.
    ! Dx is the destination vector and is creatd as Q1-vector.
    
    ! Point coordinates
    REAL(DP), DIMENSION(:,:), INTENT(IN)               :: DvertexCoords
    
    ! Vertices at the element
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

    ! Edges at the element
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement
    
    ! Element number where to calculate the streamfunction
    INTEGER(PREC_ELEMENTIDX), INTENT(IN)               :: iel
    
    ! Local number (1..NVE) of any of the vertices on the element iel
    ! where the streamfunction is already calculated. This will be used
    ! as 'base' to calculate the others.
    INTEGER, INTENT(IN)                                :: ibaseIdx
    
    ! Marker array of length NVT. All vertices where streamfunction
    ! values are calculated are marked as 1.
    INTEGER(I32), INTENT(INOUT), DIMENSION(:)          :: Imarkers
    
    ! X/Y-velocity.
    REAL(DP), DIMENSION(:), INTENT(IN)                 :: Du, Dv
    
    ! Vector of size NVT; receives the values of the streamfunction
    REAL(DP), DIMENSION(:), INTENT(INOUT)              :: Dx
  
    ! local variables
    INTEGER, PARAMETER :: NVE = 4
    INTEGER :: ive,inextidx,inextvertex,imarked
    INTEGER(PREC_VERTEXIDX) :: ivt,NVT
    INTEGER(PREC_EDGEIDX) :: imid
    REAL(DP) :: dpx1,dpx2,dpy1,dpy2,dn1,dn2
    INTEGER :: ilastMarked
  
      ilastmarked = ibaseIdx
      NVT = UBOUND(DvertexCoords,2)
  
      ! Loop over the vertices on the element. We can skip the one
      ! where the SF-value is already calculated.
      DO ive=1,NVE-1

        ! ibaseIdx is the index of a marked vertex. Calculate the "next"
        ! vertex following ilastMarked and its vertex number into inextvertex.

        inextidx=MOD(ilastMarked,NVE)+1
        inextvertex=p_IverticesAtElement(inextidx,iel)

        ! If that vertex is not marked, the streamfunction is not
        ! calculated there. Otherwise we are just looking at two 
        ! marked neighboured vertices, so there's nothing to gain here. 

        IF (Imarkers(inextvertex) .EQ. 0) THEN
        
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
          imid = IedgesAtElement (ilastMarked,iel)-NVT
          
          ! Calculate the streamfunction in inextvertex from the value
          ! in ivt by:
          !
          ! sfc(inextvertex) = sfc(ivt) + U(imid)*N
          !
          ! which is the "amount of flow over the edge (ivt,inextvertex)".

          Dx(inextvertex)=Dx(ivt)+(Du(imid)*dn1+Dv(imid)*dn2)
        
        END IF ! (Imarkers(inextvertex) == 0)

        ! Go on to the next vertex on the element to look if that one
        ! has a not marked neighbour.

        ilastMarked=inextidx
          
      END DO ! ive
    
    END SUBROUTINE

  END SUBROUTINE

END MODULE
