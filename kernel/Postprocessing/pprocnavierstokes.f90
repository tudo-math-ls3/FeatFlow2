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
  INTEGER :: ielemU, ielemP
  
  ! Number of local DOF's in U and P
  INTEGER :: idoflocU, idoflocP
  
  ! Triangulation
  TYPE (t_triangulation), POINTER :: p_rtriangulation
  
  ! An accepting the DOF's of an element.
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS), TARGET :: IdofsU, IdofsP
  
  ! Coordinates of the corners of an element
  REAL(DP), DIMENSION(NDIM2D,EL_MAXNVE) :: DCoords
  
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
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
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

    ! Get the vector data
    neqU = rvector%RvectorBlock(1)%NEQ
    neqP = rvector%RvectorBlock(3)%NEQ
    
    IF (rvector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'ppns2D_bdforces: Unsupported vector precision.'
      STOP
    END IF

    ! We support only uniform and conformal discretisation structures.
    IF (.NOT. ASSOCIATED(rvector%p_rblockDiscretisation)) THEN
      PRINT *,'ppns2D_bdforces: No discretisation structure!'
      STOP
    END IF

    IF (rvector%p_rblockDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) THEN
      PRINT *,'ppns2D_bdforces_uniform: Discretisation too complex!'
      STOP
    END IF
    
    ! Get pointers to the subvectors from the block vector
    CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
    CALL lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
    CALL lsyssc_getbase_double (rvector%RvectorBlock(3),p_DdataP)
    
    IF ((rvector%RvectorBlock(1)%isortStrategy > 0) .OR. &
        (rvector%RvectorBlock(2)%isortStrategy > 0) .OR. &
        (rvector%RvectorBlock(3)%isortStrategy > 0)) THEN
      PRINT *,'ppns2D_bdforces_uniform: Resorted vectors not supported!'
      STOP
    END IF
    
    ! Get pointers to the spatial discretisation structures of the
    ! veloctiy and pressure
    p_rdiscrU => rvector%RvectorBlock(1)%p_rspatialDiscretisation
    p_rdiscrP => rvector%RvectorBlock(3)%p_rspatialDiscretisation
    
    ! What is the actual element that is used for the discretisation?
    
    ielemU = p_rdiscrU%RelementDistribution(1)%itrialElement
    ielemP = p_rdiscrP%RelementDistribution(1)%itrialElement
    
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
    CALL storage_getbase_double (p_rtriangulation%h_DvertexParameterValue, &
                                 p_DvertexParameterValue)
    CALL storage_getbase_double2d (p_rtriangulation%h_DcornerCoordinates, &
                                   p_DvertexCoordinates)

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
    
    ! We assemble the integrap contributions separately
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
        dvtp2 = boundary_dgetMaxParVal(p_rdiscrU%p_rdomain,icp)
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
          STOP
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
        
        CALL dof_locGlobMapping(p_rdiscrU, iel, .FALSE., IdofsU)
        CALL dof_locGlobMapping(p_rdiscrP, iel, .FALSE., IdofsP)
        
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
            dpres = p_DdataP(IdofsP(idfl))*DbasP(idfl,DER_FUNC,icubp)
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
  
  END SUBROUTINE

END MODULE
