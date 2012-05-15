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
!#
!# 4.) ppns2D_bdforces_line
!#     -> Calculate the body forces acting on a boundary component
!#        of a given parametrisation.
!#     -> Line integral working on mixed meshes.
!#
!# 5.) ppns2D_bdforces_vol
!#     -> Calculate the body forces acting in some part of the domain
!#     -> Volume integral working on mixed meshes, based on a
!#        characteristic function of the subdomain where to integrate
!#
!# 6.) ppns2D_calcFluxThroughLine
!#     -> Calculate the flux through a line
!#
!# 7.) ppns_initPerfConfig
!#      -> Initialises the global performance configuration
!# </purpose>
!#########################################################################

module pprocnavierstokes

!$use omp_lib
  use basicgeometry
  use bcassemblybase
  use boundary
  use collection
  use cubature
  use derivatives
  use dofmapping
  use domainintegration
  use element
  use elementpreprocessing
  use feevaluation
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use meshregion
  use perfconfig
  use pprocintegrals
  use spatialdiscretisation
  use storage
  use transformation
  use triangulation

  implicit none

  private

!<constants>

!<constantblock description = "Identifies the underlying operator.">

  ! Use the simple gradient tensor formulation for computing the forces.
  ! This formulation exploits a special trick only available in 2D
  ! to compute the coefficients by
  !   Cd = 1/df1 int ( df2 * du_tx/dn - p n_x)
  !   Cl = 1/df1 int ( df2 * du_ty/dn - p n_y)
  ! In 3D, PPNAVST_GRADIENTTENSOR_SIMPLE coincides with PPNAVST_GRADIENTTENSOR.
  integer, parameter, public :: PPNAVST_GRADIENTTENSOR_SIMPLE = 0

  ! Use the gradient tensor formulation for computing the forces.
  ! This is the standard formulation
  !   C = 1/df1 int (df2 * sigma n - pI)
  ! with sigma = grad(u)) being the gradient tensor.
  integer, parameter, public :: PPNAVST_GRADIENTTENSOR = 1

  ! Use the deformation tensor formulation for computing the forces.
  ! This is the standard formulation
  !   C = 1/df1 int (df2 * sigma n - pI)
  ! with sigma = 1/2 (grad(u) + grad(u)^T) being the deformation tensor.
  integer, parameter, public :: PPNAVST_DEFORMATIONTENSOR = 2

!</constantblock>

!<constantblock description="Constants defining the blocking of the calculation.">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of elements to handle simultaneously in volume integration
#ifndef PPNS_NELEMSIM
  integer, parameter, public :: PPNS_NELEMSIM = 1000
#endif

!</constantblock>

!</constants>

  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: ppns_perfconfig

  !************************************************************************

  public :: ppns_initPerfConfig
  public :: ppns2D_bdforces_uniform
  public :: ppns3D_bdforces_uniform
  public :: ppns2D_streamfct_uniform
  public :: ppns2D_bdforces_line
  public :: ppns2D_bdforces_vol
  public :: ppns2D_calcFluxThroughLine

contains

  !****************************************************************************

!<subroutine>

  subroutine ppns_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<input>
  ! OPTIONAL: performance configuration that should be used to initialise
  ! the global performance configuration. If not present, the values of
  ! the legacy constants is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
!</subroutine>

    if (present(rperfconfig)) then
      ppns_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(ppns_perfconfig)
      ppns_perfconfig%NELEMSIM = PPNS_NELEMSIM
    end if

  end subroutine ppns_initPerfConfig

  !****************************************************************************

!<subroutine>

  subroutine ppns2D_bdforces_uniform (rvector,rregion,Dforces,ccub,df1,df2,cformulation)

!<description>
  ! Calculates the drag-/lift-forces acting on a part of the real
  ! boundary for a vector rvector with the solution of the 2D
  ! (Navier-)Stokes equation.
  !
  ! This routine is a wrapper for the following routines:
  ! 1.) ppns2D_bdforces_uni_line
  ! 2.) ppns2D_bdforces_uni_vol
  !
  ! Depending on whether ccub defined a 1D or 2D cubature rule, either the
  ! ppns2D_bdforces_uni_line (1D) routine or the ppns2D_bdforces_uni_vol (2D)
  ! routine is called.
!</description>

!<input>
  ! The FE solution vector.
  type(t_vectorBlock), intent(in)    :: rvector

  ! Boundary region where to calculate the boundary forces.
  ! Can be created e.g. by boundary_createRegion.
  type(t_boundaryRegion), intent(in) :: rregion

  ! 1D Cubature formula identifier to use for the line integration.
  ! One of the CUB_xxxx_1D constants in the cubature.f90.
  integer(I32), intent(in)            :: ccub

  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(in), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(in), optional      :: df2

  ! OPTIONAL: Type of tensor to use for the computation.
  ! May be either PPNAVST_GRADIENTTENSOR_XXXX for the gradient tensor or
  ! PPNAVST_GDEFORMATIONTENSOR for the deformation tensor formulation.
  ! If not present, PPNAVST_GRADIENTTENSOR_SIMPLE is assumed.
  integer, intent(in), optional :: cformulation

!</input>

!<output>
  ! Array receiving the forces acting on the boundary specified by rregion.
  real(DP), dimension(:), intent(out) :: Dforces
!</output>

!</subroutine>

    ! Determine the shape of the cubature rule
    select case(cub_igetShape(ccub))
    case (BGEOM_SHAPE_LINE)
      ! It is a 1D cubature formula - call line integration version.
      call ppns2D_bdforces_uni_line(rvector,rregion,Dforces,ccub,df1,df2,cformulation)

    case (BGEOM_SHAPE_TRIA, BGEOM_SHAPE_QUAD)
      ! It is a 2D cubature formula - call volume integration version.
      call ppns2D_bdforces_uni_vol(rvector,rregion,Dforces,ccub,df1,df2)

    case default
      ! Oops...
      call output_line ('Cubature formula not supported!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'ppns2D_bdforces_uniform')
      call sys_halt()
    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine ppns2D_bdforces_uni_line(rvector,rregion,Dforces,ccub,df1,df2,cformulation)

!<description>
  ! Calculates the drag-/lift-forces acting on a part of the real
  ! boundary for a vector rvector with the solution of the 2D
  ! (Navier-)Stokes equation. It is assumed that
  !   rvector%rvectorBlock(1) = X-velocity,
  !   rvector%rvectorBlock(2) = Y-velocity,
  !   rvector%rvectorBlock(3) = pressure.
  !
  ! rregion specifies a boundary region where to calculate
  ! the force integral. Dforces(1:NDIM2D) receives the forces in the
  ! X- and Y-direction.
  !
  ! The body forces are defined as the integrals
  ! <tex>
  !    $$ Dforces(1) = 2/df2 * \int_\Gamma [df1 dut/dn n_y - p n_x] ds $$
  !    $$ Dforces(2) = 2/df2 * \int_\Gamma [df1 dut/dn n_x + p n_y] ds $$
  ! </tex>
  ! where df1 and df2 allow to weight different integral parts
  !
  ! Usually in benchmark-like geometries there is:
  ! <tex>
  !   $$ df1 = RHO*NU = (\text{density of the fluid}) * \text{viscosity}$$
  !   $$ df2 = RHO*DIST*UMEAN**2
  !          = (\text{density of the fluid} *
  !             \text{length of the obstacle facing the flow})
  !           *(\text{mean velocity of the fluid})^2 $$
  ! </tex>
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
  type(t_vectorBlock), intent(in)    :: rvector

  ! Boundary region where to calculate the boundary forces.
  ! Can be created e.g. by boundary_createRegion.
  type(t_boundaryRegion), intent(in) :: rregion

  ! 1D Cubature formula identifier to use for the line integration.
  ! One of the CUB_xxxx_1D constants in the cubature.f90.
  integer(I32), intent(in)           :: ccub

  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(in), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(in), optional      :: df2

  ! OPTIONAL: Type of tensor to use for the computation.
  ! May be either PPNAVST_GRADIENTTENSOR_XXXX for the gradient tensor or
  ! PPNAVST_GDEFORMATIONTENSOR for the deformation tensor formulation.
  ! If not present, PPNAVST_GRADIENTTENSOR_SIMPLE is assumed.
  integer, intent(in), optional :: cformulation

!</input>

!<output>
  ! Array receiving the forces acting on the boundary specified by rregion.
  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
  real(DP), dimension(:), intent(out) :: Dforces
!</output>

!</subroutine>

  ! Spatial discretisation structure of velocity and pressure
  type(t_spatialDiscretisation), pointer :: p_rdiscrU, p_rdiscrP

  ! Element type identifier for U and P
  integer(I32) :: ielemU, ielemP

  ! Number of local DOF`s in U and P
  integer :: idoflocU, idoflocP

  ! Triangulation
  type (t_triangulation), pointer :: p_rtriangulation

  ! An accepting the DOF`s of an element.
  integer, dimension(EL_MAXNBAS), target :: IdofsU, IdofsP

  ! Coordinates of the coordinates of an element
  real(DP), dimension(:,:), allocatable :: DCoords

  ! Coordinates of the cubature points on reference and real element
  real(DP), dimension(NDIM2D,CUB_MAXCUBP_1D) :: DpointsRef,DpointsReal

  ! Coordinate system for U and P element
  integer(I32) :: ctrafoU, ctrafoP

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
  real(DP) :: dvtp1,dvtp2,dedgelen, dweight,dut,dpf1,dpf2,dval1,dval2
  real(DP), dimension(2) :: DintU, DintP
  real(DP) :: dpres
  real(DP), dimension(NDIM2D) :: dvt1,dvt2,dtangential,dnormal
  integer(I32), dimension(:), pointer :: p_ItwistIndex

  ! Type of formulation
  integer :: cform

    ! Get the vector data
    neqU = rvector%RvectorBlock(1)%NEQ
    neqP = rvector%RvectorBlock(3)%NEQ

    if (rvector%cdataType .ne. ST_DOUBLE) then
      print *,'ppns2D_bdforces: Unsupported vector precision.'
      call sys_halt()
    end if

    ! We support only uniform discretisation structures.
    if (.not. associated(rvector%p_rblockDiscr)) then
      print *,'ppns2D_bdforces_uni_line: No discretisation structure!'
      call sys_halt()
    end if

    if (rvector%p_rblockDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
      print *,'ppns2D_bdforces_uni_line: Discretisation too complex!'
      call sys_halt()
    end if

    ! Get pointers to the subvectors from the block vector
    call lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
    call lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
    call lsyssc_getbase_double (rvector%RvectorBlock(3),p_DdataP)

    if ((rvector%RvectorBlock(1)%isortStrategy > 0) .or. &
        (rvector%RvectorBlock(2)%isortStrategy > 0) .or. &
        (rvector%RvectorBlock(3)%isortStrategy > 0)) then
      print *,'ppns2D_bdforces_uni_line: Resorted vectors not supported!'
      call sys_halt()
    end if

    ! Get the used formulation; gradient or deformation tensor.
    cform = PPNAVST_GRADIENTTENSOR_SIMPLE
    if (present(cformulation)) cform = cformulation

    ! Get pointers to the spatial discretisation structures of the
    ! velocity and pressure
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
    ! These have to be mapped to the 2D element which is under consideration
    ! during the integration -- later.
    !
    ! Before, we need some additional information about the triangulation
    ! and about our element!
    !
    ! Number of local DOF`s:

    idoflocU = elem_igetNDofLoc(ielemU)
    idoflocP = elem_igetNDofLoc(ielemP)

    ! Number of local edges
    nlocaledges = elem_igetNVE (ielemU)

    ! Allocate memory for the coordinates of the element
    allocate(DCoords(NDIM2D,max(elem_igetNVE(ielemU),elem_igetNVE(ielemP))))

    ! The triangulation - it is the same for U and P
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
      ! to integrate along it? We check both endpoints...

      if (boundary_isInRegion (rregion,icp,dvtp1) .and. &
          boundary_isInRegion (rregion,icp,dvtp2)) then

        ! Ok, the edge is a boundary edge in the specified region.
        !
        ! Get the coordinates of the vertices
        dvt1 = p_DvertexCoordinates(:,ivt1)
        dvt2 = p_DvertexCoordinates(:,ivt2)

        ! Furthermore, get the tangential and the normal vector of the edge
        !
        ! The question of the sign of the tangential can be answered pretty easily.
        ! Consider for comparison the circle with radius r around the origin which is
        ! traversed in mathematical negative direction. Why negative? Because we are the
        ! inner circular obstacle has opposite direction... The circle can be
        ! parameterised by means of w(t) = (r cos t, -r sin t)^T, t \in [0, 2 \pi]. Then
        ! holds for the tangential:
        !    grad w(t) = (-r sin t, -r cos t)^T.
        ! Normalise it:
        !              = (-sin t, -cos t)^T.
        ! The FEAT2 convention to number the nodes of a quad is counterclockwise. From the
        ! point of view of a circle traversed clockwise the very same boundary points are
        ! traversed clockwise which is fine to determine the tangential vector. So, we
        ! don't need to reverse the points and get the formula for the tangential:
        dtangential = dvt2(:) - dvt1(:)

        dedgelen = sqrt(dtangential(1)**2 + dtangential(2)**2)

        dtangential = dtangential / dedgelen

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
          print *,'ppns2D_bdforces_uni_line: Edge not found. KMID destroyed?'
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

        ! For the integration, we need the global DOF`s on our element
        ! for U and P:

        call dof_locGlobMapping(p_rdiscrU, iel, IdofsU)
        call dof_locGlobMapping(p_rdiscrP, iel, IdofsP)

        ! Get the coordinates of the point on the current element
        Dcoords (1:NDIM2D,1:nlocaledges) = &
            p_DvertexCoordinates(1:NDIM2D, p_IverticesAtElement (1:nlocaledges,iel))

        ! Calculate the transformation for all points on the current element.
        ! If we have only parametric elements, we do not have to calculate
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

        ! Which tensor formulation do we have?
        select case (cform)
        case (PPNAVST_GRADIENTTENSOR_SIMPLE)
          ! 'Simple' gradient tensor formulation
          !
          ! Loop over the cubature points on the current element
          ! to assemble the integral
          do icubp = 1,ncubp

            ! Calculate the OMEGA for the integration by multiplication
            ! of the integration coefficient by the Jacobian of the
            ! mapping.
            ! The determinant of the mapping of the unit interval [-1,1]
            ! to the real line is 0.5*length of the line!

            dweight = Domega(icubp)*0.5_DP*dedgelen

            ! Loop through the DOF`s on our element and calculate
            ! the tangential velocity as well as pressure P.
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

        case (PPNAVST_GRADIENTTENSOR)
          ! 'Full' Gradient tensor formulation
          !
          ! We need to calculate the integral based on the following integrand:
          !   |-p         |   (du1/dx  du1/dy )   ( n_x )
          !   |       -p  | + (du2/dx  du2/dy ) * ( n_y )
          !
          ! Loop over the cubature points on the current element
          ! to assemble the integral
          do icubp = 1,ncubp

            ! Calculate the OMEGA for the integration by multiplication
            ! of the integration coefficient by the Jacobian of the
            ! mapping.
            ! The determinant of the mapping of the unit interval [-1,1]
            ! to the real line is 0.5*length of the line!

            dweight = Domega(icubp)*0.5_DP*dedgelen

            ! Loop through the DOF`s on our element and calculate
            ! the tangential U as well as P.
            dval1 = 0.0_DP
            dval2 = 0.0_DP
            do idfl=1,idoflocU
              dval1 = dval1 + (p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)) &
                              * Dnormal(1) &
                            + (p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)) &
                              * Dnormal(2)

              dval2 = dval2 + (p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)) &
                              * Dnormal(1) &
                            + (p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)) &
                              * Dnormal(2)
            end do

            dpres = 0.0_DP
            do idfl=1,idoflocP
              dpres = dpres + p_DdataP(IdofsP(idfl))*DbasP(idfl,DER_FUNC,icubp)
            end do

            ! Sum this up to the two integral contributions for the pressure and
            ! velocity.
            DintU(1) = DintU(1) + dweight * dval1
            DintU(2) = DintU(2) + dweight * dval2
            DintP(1) = DintP(1) - dweight * dpres * Dnormal(1)
            DintP(2) = DintP(2) - dweight * dpres * Dnormal(2)

          end do ! icubp

        case (PPNAVST_DEFORMATIONTENSOR)

          ! Deformation tensor formulation
          !
          ! We need to calculate the integral based on the following integrand:
          !   |-p         |   (du1/dx               1/2 (du1/dy+du2/dx) )   ( n_x )
          !   |       -p  | + (1/2 (du1/dy+du2/dx)  du2/dy              ) * ( n_y )
          !
          ! Loop over the cubature points on the current element
          ! to assemble the integral
          do icubp = 1,ncubp

            ! Calculate the OMEGA for the integration by multiplication
            ! of the integration coefficient by the Jacobian of the
            ! mapping.
            ! The determinant of the mapping of the unit interval [-1,1]
            ! to the real line is 0.5*length of the line!

            dweight = Domega(icubp)*0.5_DP*dedgelen

            ! Loop through the DOF`s on our element and calculate
            ! the tangential U as well as P.
            dval1 = 0.0_DP
            dval2 = 0.0_DP
            do idfl=1,idoflocU
              dval1 = dval1 + 2.0_DP*(p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)) &
                              * Dnormal(1) &
                            + (p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)) &
                              * Dnormal(2) &
                            + (p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)) &
                              * Dnormal(2)

              dval2 = dval2 + 2.0_DP*(p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)) &
                              * Dnormal(2) &
                            + (p_DdataUX(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_Y,icubp)) &
                              * Dnormal(1) &
                            + (p_DdataUY(IdofsU(idfl)) * DbasU(idfl,DER_DERIV_X,icubp)) &
                              * Dnormal(1)
            end do

            dpres = 0.0_DP
            do idfl=1,idoflocP
              dpres = dpres + p_DdataP(IdofsP(idfl))*DbasP(idfl,DER_FUNC,icubp)
            end do

            ! Sum this up to the two integral contributions for the pressure and
            ! velocity.
            DintU(1) = DintU(1) + dweight * dval1
            DintU(2) = DintU(2) + dweight * dval2
            DintP(1) = DintP(1) - dweight * dpres * Dnormal(1)
            DintP(2) = DintP(2) - dweight * dpres * Dnormal(2)

          end do ! icubp


        end select

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

  subroutine ppns2D_bdforces_uni_vol (rvector,rregion,Dforces,ccub,df1,df2)

!<description>
  ! Calculates the drag-/lift-forces acting on a part of the real
  ! boundary for a vector rvector with the solution of the 2D
  ! (Navier-)Stokes equation. It is assumed that
  !   rvector%rvectorBlock(1) = X-velocity,
  !   rvector%rvectorBlock(2) = Y-velocity,
  !   rvector%rvectorBlock(3) = pressure.
  !
  ! rregion specifies a boundary region where to calculate
  ! the force integral. Dforces(1:NDIM2D) receives the forces in the
  ! X- and Y-direction.
  !
  !
  ! Double precision version for all Finite Elements without
  ! curved boundaries.
!</description>

!<input>
  ! The FE solution vector.
  type(t_vectorBlock), intent(in)    :: rvector

  ! Boundary region where to calculate the boundary forces.
  ! Can be created e.g. by boundary_createRegion.
  type(t_boundaryRegion), intent(in) :: rregion

  ! 2D Cubature formula identifier to use for the volume integration.
  ! One of the CUB_xxxx constants in the cubature.f90.
  integer(I32), intent(in)            :: ccub

  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(in), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(in), optional      :: df2

!</input>

!<output>
  ! Array receiving the forces acting on the boundary specified by rregion.
  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
  real(DP), dimension(:), intent(out) :: Dforces
!</output>

!</subroutine>

  ! Spatial discretisation structure of velocity and pressure
  type(t_spatialDiscretisation), pointer :: p_rdiscrU, p_rdiscrP

  ! Element type identifier for U and P
  integer(I32) :: ielemU, ielemP

  ! Number of local DOF`s in U and P
  integer :: idoflocU, idoflocP

  ! Triangulation
  type (t_triangulation), pointer :: p_rtriangulation

  ! An accepting the DOF`s of an element.
  integer, dimension(EL_MAXNBAS), target :: IdofsU, IdofsP

  ! Coordinates of the coordinates of an element
  real(DP), dimension(:,:), allocatable :: DCoords

  ! Coordinates of the cubature points on reference and real element
  real(DP), dimension(NDIM3D,CUB_MAXCUBP) :: DpointsRef,DpointsReal

  ! Cubature point coordinates on the reference element.
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi2D

  ! Coordinate system for U and P element
  integer(I32) :: ctrafoU, ctrafoP

  ! U/P element parametric or nonparametric
  logical :: bnonparU, bnonparP

  ! Arrays for saving Jacobian determinants and matrices
  real(DP), dimension(CUB_MAXCUBP_2D) :: Ddetj
  real(DP), dimension(EL_NJACENTRIES2D,CUB_MAXCUBP_2D) :: Djac

  ! Array to tell the element which derivatives to calculate.
  logical, dimension(EL_MAXNDER) :: BderU, BderP

  ! Value of basis functions
  real(DP), dimension(EL_MAXNBAS,EL_MAXNDER,CUB_MAXCUBP_2D) :: DbasU, DbasP

  ! Value of local alpha function
  real(DP), dimension(EL_MAXNBAS) :: Dalpha

  ! Pointer to vector data of solution vector
  real(DP), dimension(:), pointer :: p_DdataUX,p_DdataUY,p_DdataP


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
  real(DP), dimension(:), pointer :: p_DvertexParameterValue
  integer, dimension(:), pointer :: p_IboundaryCpIdx
  real(DP), dimension(:,:), pointer :: p_DvertexCoordinates

  ! other local variables
  integer :: iedgeidx,ivt1,ivt2,ilocaledge,nlocaledges,idfl,icp
  integer :: neqU,neqP,i,j
  integer :: iedge,iedgeglobal
  integer :: iel
  real(DP) :: dvtp1,dvtp2, dweight,dpf1,dpf2
  real(DP) :: du1x,du1y,du1v,du2x,du2y,du2v,dpv,dax,day
  real(DP), dimension(2) :: DintU, DintP, DintV
  real(DP), dimension(NDIM2D) :: dvt1,dvt2
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
    ! velocity and pressure
    p_rdiscrU => rvector%RvectorBlock(1)%p_rspatialDiscr
    p_rdiscrP => rvector%RvectorBlock(3)%p_rspatialDiscr

    ! What is the actual element that is used for the discretisation?
    ielemU = p_rdiscrU%RelementDistr(1)%celement
    ielemP = p_rdiscrP%RelementDistr(1)%celement

    ! So far so good, we have checked that the assumptions for the integration
    ! are fulfilled. Now we can start the actual integration.
    !
    ! At first initialise a 2D cubature formula for the integration.
    call cub_getCubPoints (ccub,ncubp,Dxi2D,Domega)

    ! Map the coordinates in a proper 2D array
    do i = 1, ubound(DpointsRef,1)
      do j = 1, ubound(DpointsRef,2)
        DpointsRef(i,j) = Dxi2D(j,i)
      end do
    end do
    !DpointsRef (1:NDIM2D,1:ncubp) = transpose(Dxi2D(1:ncubp,1:NDIM2D))

    ! In Dxi1D we have the 1D coordinates of the cubature points.
    ! These have to be mapped to the 2D element which is under consideration
    ! during the integration -- later.
    !
    ! Before, we need some additional information about the triangulation
    ! and about our element!
    !
    ! Number of local DOF`s:

    idoflocU = elem_igetNDofLoc(ielemU)
    idoflocP = elem_igetNDofLoc(ielemP)

    ! Number of local edges
    nlocaledges = elem_igetNVE (ielemU)

    ! Allocate memory for the coordinates of the element
    allocate(Dcoords(NDIM2D,max(elem_igetNVE(ielemU),elem_igetNVE(ielemP))))


    ! The triangulation - it is the same for U and P
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
    BderU(DER_FUNC) = .true.
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
    DintV = 0.0_DP

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
      ! to integrate along it? We check both endpoints...

      if (boundary_isInRegion (rregion,icp,dvtp1) .and. &
          boundary_isInRegion (rregion,icp,dvtp2)) then

        ! Ok, the edge is a boundary edge in the specified region.
        !
        ! Get the coordinates of the vertices
        dvt1 = p_DvertexCoordinates(:,ivt1)
        dvt2 = p_DvertexCoordinates(:,ivt2)

        ! Now quickly determine which edge 1..4 of the current element
        ! is the one we are treating at the moment.
        do ilocaledge = 1,nlocaledges
          if (p_IedgesAtElement(ilocaledge,iel) .eq. iedgeglobal) exit
        end do

        if (ilocaledge .gt. nlocaledges) then
          print *,'ppns2D_bdforces: Edge not found. KMID destroyed?'
          call sys_halt()
        end if

!        select case(ilocaledge)
!        case(1)
!          darx =  0.0_DP
!          dary = -0.5_DP
!        case(2)
!          darx =  0.5_DP
!          dary =  0.0_DP
!        case(3)
!          darx =  0.0_DP
!          dary =  0.5_DP
!        case(4)
!          darx = -0.5_DP
!          dary =  0.0_DP
!        end select

        ! Okay, now let us build the local alpha vector
        Dalpha = 0.0_DP
        select case(elem_getPrimaryElement(ielemU))
        case(EL_P1)
          ! Set both vertices adjacent to the edge to 1
          Dalpha(ilocaledge) = 1.0_DP
          Dalpha(mod(ilocaledge,3)+1) = 1.0_DP

        case(EL_P2)
          ! Set both vertices adjacent to the edge to 1
          Dalpha(ilocaledge) = 1.0_DP
          Dalpha(mod(ilocaledge,3)+1) = 1.0_DP
          ! Set the edge value to 1
          Dalpha(3+ilocaledge) = 1.0_DP

        case(EL_Q1,EL_QPW4P1_2D)
          ! Set both vertices adjacent to the edge to 1
          Dalpha(ilocaledge) = 1.0_DP
          Dalpha(mod(ilocaledge,4)+1) = 1.0_DP

        case(EL_Q2,EL_QPW4P2_2D)
          ! Set both vertices adjacent to the edge to 1
          Dalpha(ilocaledge) = 1.0_DP
          Dalpha(mod(ilocaledge,4)+1) = 1.0_DP
          ! Set the edge value to 1
          Dalpha(4+ilocaledge) = 1.0_DP

        case(EL_P1T,EL_Q1T,EL_Q1TB,EL_Q2T,EL_Q2TB,EL_Q3T_2D)
          ! Set the first moment on the edge to 1
          Dalpha(ilocaledge) = 1.0_DP

        case default
          call output_line ('Velocity element not supported!', &
                            OU_CLASS_ERROR,OU_MODE_STD, &
                            'ppns2D_bdforces_uni_vol')
          call sys_halt()
        end select


        ! For the integration, we need the global DOF`s on our element
        ! for U and P:
        call dof_locGlobMapping(p_rdiscrU, iel, IdofsU)
        call dof_locGlobMapping(p_rdiscrP, iel, IdofsP)

        ! Get the coordinates of the point on the current element
        Dcoords (1:NDIM2D,1:nlocaledges) = &
            p_DvertexCoordinates(1:NDIM2D, p_IverticesAtElement (1:nlocaledges,iel))

        ! Calculate the transformation for all points on the current element.
        ! If we have only parametric elements, we do not have to calculate
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
          dweight = Domega(icubp)*Ddetj(icubp)

          !
          du1v = 0.0_DP
          du1x = 0.0_DP
          du1y = 0.0_DP
          du2v = 0.0_DP
          du2x = 0.0_DP
          du2y = 0.0_DP
          dax = 0.0_DP
          day = 0.0_DP
          do idfl = 1, idoflocU
            du1v = du1v + DbasU(idfl,DER_FUNC,icubp) * p_DdataUX(IdofsU(idfl))
            du1x = du1x + DbasU(idfl,DER_DERIV_X,icubp) * p_DdataUX(IdofsU(idfl))
            du1y = du1y + DbasU(idfl,DER_DERIV_Y,icubp) * p_DdataUX(IdofsU(idfl))
            du2v = du2v + DbasU(idfl,DER_FUNC,icubp) * p_DdataUY(IdofsU(idfl))
            du2x = du2x + DbasU(idfl,DER_DERIV_X,icubp) * p_DdataUY(IdofsU(idfl))
            du2y = du2y + DbasU(idfl,DER_DERIV_Y,icubp) * p_DdataUY(IdofsU(idfl))
            dax = dax + DbasU(idfl,DER_DERIV_X,icubp) * Dalpha(idfl)
            day = day + DbasU(idfl,DER_DERIV_Y,icubp) * Dalpha(idfl)
          end do

          !dav = 0.5_DP + darx*DpointsRef(1,icubp) + dary*DpointsRef(2,icubp)
          !dax =  (Djac(4,icubp)*darx - Djac(2,icubp)*dary) / Ddetj(icubp)
          !day = -(Djac(3,icubp)*darx - Djac(1,icubp)*dary) / Ddetj(icubp)

          dpv = 0.0_DP
          do idfl = 1, idoflocP
            dpv = dpv + DbasP(idfl,DER_FUNC,icubp)*p_DdataP(IdofsP(idfl))
          end do

          ! Sum this up to the two integral contributions for the pressure and
          ! velocity.
          DintU(1) = DintU(1) - dweight * (du1x * dax + du1y * day)
          DintU(2) = DintU(2) - dweight * (du2x * dax + du2y * day)
          DintP(1) = DintP(1) + dweight * dpv * dax
          DintP(2) = DintP(2) + dweight * dpv * day

          !DintV(1) = DintV(1) + dweight * dav * (du1v * du1x + du2v * du1y)
          !DintV(2) = DintV(2) + dweight * dav * (du1v * du2x + du2v * du2y)

        end do ! icubp

      end if

    end do ! iedgeidx

    ! DintU and DintP give now the contributions to the force integral:
    !
    ! DragCoeff = 2/dfp2 * (dpf1*DintU(1) + DintP(1))
    ! LiftCoeff = 2/dfp2 * (dpf1*DintU(2) + DintP(2))

    Dforces = 0.0_DP
    Dforces(1:NDIM2D) = 2.0_DP/dpf2 * (dpf1*DintU(:) + DintP(:))
    !Dforces(1:NDIM2D) = 2.0_DP/dpf2 * (dpf1*DintU(:) + DintV(:) + DintP(:))

    ! Deallocate memory, finish.
    deallocate(DCoords)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine ppns2D_bdforces_line(rvector,rregion,Dforces,ccub,df1,df2,cformulation,&
      ffunctionReference,rcollection,rperfconfig)

!<description>
  ! Calculates the drag-/lift-forces acting on a part of the real
  ! boundary for a vector rvector with the solution of the 2D
  ! (Navier-)Stokes equation. It is assumed that
  !   rvector%rvectorBlock(1) = X-velocity,
  !   rvector%rvectorBlock(2) = Y-velocity,
  !   rvector%rvectorBlock(3) = pressure.
  !
  ! rregion specifies a boundary region where to calculate
  ! the force integral. Dforces(1:NDIM2D) receives the forces in the
  ! X- and Y-direction.
  !
  ! The body forces are defined as the integrals
  ! <tex>
  !    $$ Dforces(1) = 2/df2 * \int_\Gamma [df1 dut/dn n_y - p n_x] ds $$
  !    $$ Dforces(2) = 2/df2 * \int_\Gamma [df1 dut/dn n_x + p n_y] ds $$
  ! </tex>
  ! where df1 and df2 allow to weight different integral parts
  !
  ! Usually in benchmark-like geometries there is:
  ! <tex>
  !   $$ df1 = RHO*NU = (\text{density of the fluid}) * \text{viscosity}$$
  !   $$ df2 = RHO*DIST*UMEAN**2
  !          = (\text{density of the fluid} *
  !             \text{length of the obstacle facing the flow})
  !           *(\text{mean velocity of the fluid})^2 $$
  ! </tex>
  ! which influences this routine to calculate the drag-/lift coefficients
  ! instead of the forces. The actual forces are calculated by
  ! setting df1=1.0 and df2=2.0, which is the standard setting when
  ! neglecting df1 and df2.
  !
  ! Double precision version for all Finite Elements without
  ! curved boundaries. Nonconstant df1 supported.
!</description>

!<input>
  ! The FE solution vector.
  type(t_vectorBlock), intent(in)    :: rvector

  ! Boundary region where to calculate the boundary forces.
  ! Can be created e.g. by boundary_createRegion.
  type(t_boundaryRegion), intent(in) :: rregion

  ! 1D Cubature formula identifier to use for the line integration.
  ! One of the CUB_xxxx_1D constants in the cubature.f90.
  integer(I32), intent(in)            :: ccub

  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(in), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(in), optional      :: df2

  ! OPTIONAL: Type of tensor to use for the computation.
  ! May be either PPNAVST_GRADIENTTENSOR_XXXX for the gradient tensor or
  ! PPNAVST_GDEFORMATIONTENSOR for the deformation tensor formulation.
  ! If not present, PPNAVST_GRADIENTTENSOR_SIMPLE is assumed.
  integer, intent(in), optional :: cformulation

  ! A callback routine that allows to specify a nonconstant coefficient df1.
  ! If not present, df1 is used as constant coefficient in the integral.
  ! If present, df1 is ignored and calculated by cubature point
  ! by this callback routine.
  include '../Postprocessing/intf_refFunctionSc.inc'
  optional :: ffunctionReference

  ! OPTIONAL: A collection structure that is passed to ffunctionRefSimple.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<output>
  ! Array receiving the forces acting on the boundary specified by rregion.
  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
  real(DP), dimension(:), intent(out) :: Dforces
!</output>

!</subroutine>

  ! Triangulation
  type (t_triangulation), pointer :: p_rtriangulation

  ! Elements in the boundary region
  integer :: nelementsInRegion, nelements, iel
  integer, dimension(:), allocatable, target :: IelementsInRegion, IcurrentElementSet
  integer, dimension(:), allocatable :: IelementPointer

  ! Edges in the BC-region, local index
  integer, dimension(:), allocatable, target :: IedgesLocal

  ! Element distributions and current discretisation
  integer :: ieldistr
  type(t_spatialDiscretisation), pointer :: p_rdiscr
  integer, dimension(:), pointer :: p_IelementDistr

  ! A t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  type(t_domainIntSubset) :: rintSubset

  ! Coordinate system for reference coordinate and current element type
  integer(I32) :: celement,icoordSystem,ctrafotype

  ! Cubature point coordinates on the reference element.
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D,Dxi2D

  ! Coordinates of the cubature points on reference element
  real(DP), dimension(:,:,:), allocatable :: DpointsRef

  ! An element set for evaluating FE spaces.
  type(t_evalElementSet) :: revalElementSet

  ! DOF-mapping and DOF-values on the current element set
  integer, dimension(:,:), allocatable, target :: Idofs
  integer :: ndofs,idfl,idof
  real(DP), dimension(:,:,:,:), allocatable :: Dbas
  logical, dimension(EL_MAXNDER) :: Bder

  ! variables to calculate FE values
  real(DP) :: dv1,dv2,dv3,dv4

  ! pointers to the FEM solution vector
  real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3

  ! Values in the cubature points
  real(DP), dimension(:,:,:), allocatable :: Dvalues

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(CUB_MAXCUBP) :: Domega

  ! number/index of cubature points on the reference element
  integer :: ncubp,icubp

  ! Edges, vertices and elements on the boundary
  integer, dimension(:,:), pointer  :: p_IedgesAtElement,p_IverticesAtEdge
  real(DP), dimension(:,:), pointer  :: p_DvertexCoords

  ! Variables that get the integral values of the parts
  real(DP), dimension(2) :: DintU, DintP

  ! other local variables
  integer :: iedge,ivt1,ivt2
  real(DP) :: dpf1,dpf2
  real(DP), dimension(NDIM2D) :: Dvt1,Dvt2,Dtangential,Dnormal
  real(DP) :: dedgelen,dweight,dut,dpres,dval1,dval2

  ! Type of formulation
  integer :: cform

    ! Some basic checks...

    if (rvector%cdataType .ne. ST_DOUBLE) then
      print *,'ppns2D_bdforces: Unsupported vector precision.'
      call sys_halt()
    end if

    ! We support only uniform discretisation structures.
    if (.not. associated(rvector%p_rblockDiscr)) then
      print *,'ppns2D_bdforces: No discretisation structure!'
      call sys_halt()
    end if

    if ((rvector%RvectorBlock(1)%isortStrategy > 0) .or. &
        (rvector%RvectorBlock(2)%isortStrategy > 0) .or. &
        (rvector%RvectorBlock(3)%isortStrategy > 0)) then
      print *,'ppns2D_bdforces_uniform: Resorted vectors not supported!'
      call sys_halt()
    end if

    ! Get the used formulation; gradient or deformation tensor.
    cform = PPNAVST_GRADIENTTENSOR_SIMPLE
    if (present(cformulation)) cform = cformulation

    ! Prepare the weighting coefficients
    dpf1 = 1.0_DP
    dpf2 = 2.0_DP
    if (present(df1)) dpf1 = df1
    if (present(df2)) dpf2 = df2

    ! Preparation of derivatives
    Bder(:) = .false.
    Bder(DER_FUNC) = .true.
    Bder(DER_DERIV2D_X) = .true.
    Bder(DER_DERIV2D_Y) = .true.

    ! The triangulation - it is the same for U and P
    p_rtriangulation => rvector%p_rblockDiscr%RspatialDiscr(1)%p_rtriangulation

    ! The arrays which contain the elements and edge numbers on the boundary
    ! as well as boundary index array.
    ! Fetch geometrical information
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,p_DvertexCoords)

    if (p_rtriangulation%h_DvertexParameterValue .eq. ST_NOHANDLE) then
      print *,'No boundary parameters available!'
      call sys_halt()
    end if

    ! Get the FE solution vectors
    call lsyssc_getbase_double (Rvector%RvectorBlock(1),p_Ddata1)
    call lsyssc_getbase_double (Rvector%RvectorBlock(2),p_Ddata2)
    call lsyssc_getbase_double (Rvector%RvectorBlock(3),p_Ddata3)

    ! In a very first step, figure out how many elements we have
    ! in our boundary region.
    call bcasm_getElementsInBdRegion (p_rtriangulation,rregion, nelementsInRegion)
    allocate(IelementsInRegion(nelementsInRegion))
    allocate(IedgesLocal(nelementsInRegion))
    call bcasm_getElementsInBdRegion (p_rtriangulation,rregion, nelementsInRegion, &
        IelementsInRegion,IedgeLocal=IedgesLocal)

    ! A set of elements, currently under consideration
    allocate(IcurrentElementSet(nelementsInRegion))

    ! IelementPointer points for each element in IcurrentElementSet to
    ! the original position in IelementsInRegion.
    allocate(IelementPointer(nelementsInRegion))

    ! Initialise a 1D cubature formula for the boundary integration.
    call cub_getCubPoints (ccub,ncubp,Dxi1D,Domega)

    ! Allocate memory for the values in the cubature points:
    ! Dvalues(:,:,1) = du/dx
    ! Dvalues(:,:,2) = du/dy
    ! Dvalues(:,:,3) = dv/dx
    ! Dvalues(:,:,4) = dv/dy
    ! Dvalues(:,:,5) = p
    ! Dvalues(:,:,6) = dpf1
    allocate(Dvalues(ncubp,nelementsInRegion,NDIM2D*2+1+1))

    ! We assemble the integral contributions separately
    DintU = 0.0_DP
    DintP = 0.0_DP

    ! Loop over the element distributions of the velocity vectors
    p_rdiscr => rvector%p_rblockDiscr%RspatialDiscr(1)
    do ieldistr = 1,p_rdiscr%inumFESpaces

      ! Find all elements within this element discribution.
      ! Save them in IelementsInRegion. This array will therefore receive
      ! all the elements we can handle simultaneously.
      if (p_rdiscr%h_IelementDistr .eq. ST_NOHANDLE) then
        nelements = nelementsInRegion
        call lalg_copyVectorInt(IelementsInRegion,IcurrentElementSet)
        do iel=1,nelements
          IelementPointer(iel) = iel
        end do
      else
        call storage_getbase_int(p_rdiscr%h_IelementDistr,p_IelementDistr)
        nelements = 0
        do iel = 1,nelementsInRegion
          if (p_IelementDistr(IcurrentElementSet(iel)) .eq. ieldistr) then
            nelements = nelements + 1
            IcurrentElementSet(nelements) = IcurrentElementSet(iel)
            IelementPointer(nelements) = iel
          end if
        end do
      end if

      ! What is the current element type?
      celement = p_rdiscr%RelementDistr(ieldistr)%celement

      ! Transformation
      ctrafotype = p_rdiscr%RelementDistr(ieldistr)%ctrafotype

      ! Figure out the coordinate system used in the current element set
      icoordSystem = elem_igetCoordSystem(celement)

      ! Allocate memory for the cubature points
      allocate(DpointsRef(trafo_igetReferenceDimension(ctrafoType),ncubp,nelements))

      ! Allocate memory for the points on the reference element(s)
      ! and map the 1D cubature points to the 2D reference element.
      do iel = 1,nelements

        ! Map the 1D points to 2D.
        call trafo_mapCubPts1Dto2D(icoordSystem, IedgesLocal(IelementPointer(iel)),&
            ncubp, Dxi1D, Dxi2D)
        DpointsRef (1:NDIM2D,1:ncubp,iel) = transpose(Dxi2D(1:ncubp,1:NDIM2D))

      end do

      ! Initialise an element set
      call elprep_init (revalElementSet)

      ! Allocate memory for the DOF`s on all these elements.
      ! Let us hope, there are not too many, otherwise the routine may have
      ! to be rewritten to work in chunks
      ndofs = elem_igetNDofLoc(celement)
      allocate(Idofs(ndofs,nelements))
      allocate(Dbas(ndofs,elem_getMaxDerivative(celement),ncubp,nelements))

      ! With the DOF-Mapping, get the DOF`s.
      call dof_locGlobMapping_mult(p_rdiscr, IcurrentElementSet(1:nelements), Idofs)

      ! Prepare an element set for evaluation on the current set
      call elprep_prepareSetForEvaluation (revalElementSet, &
          elem_getEvaluationTag(celement), p_rtriangulation, &
          IcurrentElementSet(1:nelements), ctrafoType, DpointsRef=DpointsRef,&
          rperfconfig=rperfconfig)

      ! In case, dpf1 is nonconstant, calculate dpf1.
      if (present(ffunctionReference)) then
        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
        !rintSubset%ielementDistribution = ieldistr
        rintSubset%ielementStartIdx = 1
        rintSubset%p_Ielements => IcurrentElementSet(1:nelements)
        rintSubset%p_IdofsTrial => Idofs
        rintSubset%celement = celement
        call ffunctionReference (DER_FUNC,p_rdiscr, &
                  nelements,ncubp,revalElementSet%p_DpointsReal, &
                  Idofs,rintSubset,Dvalues(:,:,6),rcollection)
        call domint_doneIntegration (rintSubset)
      else
        Dvalues(:,:,6) = dpf1
      end if

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (celement, revalElementSet, Bder, Dbas)

      ! Now perform a loop over the elements, cubature points and basis
      ! functions to calculate the values of the basis functions we need
      ! to Dvalues.
      do iel = 1,nelements
        do icubp = 1,ncubp
          dv1 = 0.0_DP
          dv2 = 0.0_DP
          dv3 = 0.0_DP
          dv4 = 0.0_DP
          do idfl = 1,ndofs
            idof = Idofs(idfl,iel)
            dv1 = dv1 + p_Ddata1(idof)*Dbas(idfl,DER_DERIV_X,icubp,iel)
            dv2 = dv2 + p_Ddata1(idof)*Dbas(idfl,DER_DERIV_Y,icubp,iel)
            dv3 = dv3 + p_Ddata2(idof)*Dbas(idfl,DER_DERIV_X,icubp,iel)
            dv4 = dv4 + p_Ddata2(idof)*Dbas(idfl,DER_DERIV_Y,icubp,iel)
          end do
          Dvalues(icubp,IelementPointer(iel),1) = dv1
          Dvalues(icubp,IelementPointer(iel),2) = dv2
          Dvalues(icubp,IelementPointer(iel),3) = dv3
          Dvalues(icubp,IelementPointer(iel),4) = dv4
        end do
      end do

      call elprep_releaseElementSet(revalElementSet)
      deallocate(Dbas)
      deallocate(Idofs)
      deallocate(DpointsRef)

    end do

    ! Now do the same for the pressure.

    ! Loop over the element distributions of the velocity vectors
    p_rdiscr => rvector%p_rblockDiscr%RspatialDiscr(3)
    do ieldistr = 1,p_rdiscr%inumFESpaces

      ! Find all elements within this element discribution.
      ! Save them in IelementsInRegion. This array will therefore receive
      ! all the elements we can handle simultaneously.
      if (p_rdiscr%h_IelementDistr .eq. ST_NOHANDLE) then
        nelements = nelementsInRegion
        call lalg_copyVectorInt(IelementsInRegion,IcurrentElementSet)
        do iel=1,nelements
          IelementPointer(iel) = iel
        end do
      else
        call storage_getbase_int(p_rdiscr%h_IelementDistr,p_IelementDistr)
        nelements = 0
        do iel = 1,nelementsInRegion
          if (p_IelementDistr(IcurrentElementSet(iel)) .eq. ieldistr) then
            nelements = nelements + 1
            IcurrentElementSet(nelements) = IcurrentElementSet(iel)
            IelementPointer(nelements) = iel
          end if
        end do
      end if

      ! What is the current element type?
      celement = p_rdiscr%RelementDistr(ieldistr)%celement

      ! Transformation
      ctrafotype = p_rdiscr%RelementDistr(ieldistr)%ctrafotype

      ! Figure out the coordinate system used in the current element set
      icoordSystem = elem_igetCoordSystem(celement)

      ! Allocate memory for the cubature points
      allocate(DpointsRef(trafo_igetReferenceDimension(ctrafoType),ncubp,nelements))

      ! Allocate memory for the points on the reference element(s)
      ! and map the 1D cubature points to the 2D reference element.
      do iel = 1,nelements

        ! Map the 1D points to 2D.
        call trafo_mapCubPts1Dto2D(icoordSystem, IedgesLocal(IelementPointer(iel)),&
            ncubp, Dxi1D, Dxi2D)
        DpointsRef (1:NDIM2D,1:ncubp,iel) = transpose(Dxi2D(1:ncubp,1:NDIM2D))

      end do

      ! Initialise an element set
      call elprep_init (revalElementSet)

      ! Allocate memory for the DOF`s on all these elements.
      ! Let us hope, there are not too many, otherwise the routine may have
      ! to be rewritten to work in chunks
      ndofs = elem_igetNDofLoc(celement)
      allocate(Idofs(ndofs,nelements))
      allocate(Dbas(ndofs,elem_getMaxDerivative(celement),ncubp,nelements))

      ! With the DOF-Mapping, get the DOF`s.
      call dof_locGlobMapping_mult(p_rdiscr, IcurrentElementSet(1:nelements), Idofs)

      ! Prepare an element set for evaluation on the current set
      call elprep_prepareSetForEvaluation (revalElementSet, &
          elem_getEvaluationTag(celement), p_rtriangulation, &
          IcurrentElementSet(1:nelements), ctrafoType, DpointsRef=DpointsRef,&
          rperfconfig=rperfconfig)

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (celement, revalElementSet, Bder, Dbas)

      ! Now perform a loop over the elements, cubature points and basis
      ! functions to calculate the values of the basis functions we need
      ! to Dvalues.
      do iel = 1,nelements
        do icubp = 1,ncubp
          dv1 = 0.0_DP
          do idfl = 1,ndofs
            idof = Idofs(idfl,iel)
            dv1 = dv1 + p_Ddata3(idof)*Dbas(idfl,DER_FUNC,icubp,iel)
          end do
          Dvalues(icubp,IelementPointer(iel),5) = dv1
        end do
      end do

      call elprep_releaseElementSet(revalElementSet)
      deallocate(Dbas)
      deallocate(Idofs)
      deallocate(DpointsRef)

    end do

    ! That was the most awful part, now it gets easier...

    ! Which tensor formulation do we have?
    select case (cform)
    case (PPNAVST_GRADIENTTENSOR_SIMPLE)
      ! 'Simple' gradient tensor formulation
      !
      ! Loop over the cubature points on the current element
      ! to assemble the integral
      do iel=1,nelementsInRegion

        ! Get the length of the edge on that element, the normal vector,
        ! and the tangential vector.
        iedge = p_IedgesAtElement(IedgesLocal(iel),IelementsInRegion(iel))
        ivt1 = p_IverticesAtEdge(1,iedge)
        ivt2 = p_IverticesAtEdge(2,iedge)
        Dvt1 = p_DvertexCoords(:,ivt1)
        Dvt2 = p_DvertexCoords(:,ivt2)

        Dtangential = Dvt2(:)-Dvt1(:)

        dedgelen = sqrt(dtangential(1)**2 + dtangential(2)**2)

        Dtangential = Dtangential / dedgelen

        dnormal(1) = -Dtangential(2)
        dnormal(2) = Dtangential(1)

        do icubp = 1,ncubp

          ! Get dpf1 in the cubature point. It is an additional weight.
          dpf1 = Dvalues(icubp,iel,6)

          ! Calculate the OMEGA for the integration by multiplication
          ! of the integration coefficient by the Jacobian of the
          ! mapping.
          ! The determinant of the mapping of the unit interval [-1,1]
          ! to the real line is 0.5*length of the line!

          dweight = Domega(icubp)*0.5_DP*dedgelen

          ! calculate the tangential U as well as P.
          dut = Dvalues(icubp,iel,1)*Dtangential(1) * Dnormal(1) &
              + Dvalues(icubp,iel,3)*Dtangential(2) * Dnormal(1) &
              + Dvalues(icubp,iel,2)*Dtangential(1) * Dnormal(2) &
              + Dvalues(icubp,iel,4)*Dtangential(2) * Dnormal(2)

          dpres = Dvalues(icubp,iel,5)

          ! Sum this up to the two integral contributions for the pressure and
          ! velocity.
          DintU(1) = DintU(1) + dweight * dpf1 * dut * Dnormal(2)
          DintU(2) = DintU(2) - dweight * dpf1 * dut * Dnormal(1)
          DintP(1) = DintP(1) - dweight * dpres * Dnormal(1)
          DintP(2) = DintP(2) - dweight * dpres * Dnormal(2)

        end do ! icubp

      end do ! iel

    case (PPNAVST_GRADIENTTENSOR)
      ! 'Full' Gradient tensor formulation
      !
      ! We need to calculate the integral based on the following integrand:
      !   |-p         |   (du1/dx  du1/dy )   ( n_x )
      !   |       -p  | + (du2/dx  du2/dy ) * ( n_y )
      !
      ! Loop over the cubature points on the current element
      ! to assemble the integral
      do iel=1,nelementsInRegion

        ! Get the length of the edge on that element, the normal vector,
        ! and the tangential vector.
        iedge = p_IedgesAtElement(IedgesLocal(iel),IelementsInRegion(iel))
        ivt1 = p_IverticesAtEdge(1,iedge)
        ivt2 = p_IverticesAtEdge(2,iedge)
        Dvt1 = p_DvertexCoords(:,ivt1)
        Dvt2 = p_DvertexCoords(:,ivt2)

        Dtangential = Dvt2(:)-Dvt1(:)

        dedgelen = sqrt(dtangential(1)**2 + dtangential(2)**2)

        Dtangential = Dtangential / dedgelen

        dnormal(1) = -Dtangential(2)
        dnormal(2) = Dtangential(1)

        do icubp = 1,ncubp

          ! Calculate the OMEGA for the integration by multiplication
          ! of the integration coefficient by the Jacobian of the
          ! mapping.
          ! The determinant of the mapping of the unit interval [-1,1]
          ! to the real line is 0.5*length of the line!

          dweight = Domega(icubp)*0.5_DP*dedgelen

          ! calculate the contributions in X and y direction
          dval1 =   Dvalues(icubp,iel,1) * Dnormal(1) &
                  + Dvalues(icubp,iel,2) * Dnormal(2)

          dval2 =   Dvalues(icubp,iel,3) * Dnormal(1) &
                  + Dvalues(icubp,iel,4) * Dnormal(2)

          dpres = Dvalues(icubp,iel,5)

          ! Sum this up to the two integral contributions for the pressure and
          ! velocity.
          DintU(1) = DintU(1) + dweight * dpf1 * dval1
          DintU(2) = DintU(2) + dweight * dpf1 * dval2
          DintP(1) = DintP(1) - dweight * dpres * Dnormal(1)
          DintP(2) = DintP(2) - dweight * dpres * Dnormal(2)

        end do ! icubp

      end do ! iel

    case (PPNAVST_DEFORMATIONTENSOR)

      ! Deformation tensor formulation
      !
      ! We need to calculate the integral based on the following integrand:
      !   |-p         |   (du1/dx               1/2 (du1/dy+du2/dx) )   ( n_x )
      !   |       -p  | + (1/2 (du1/dy+du2/dx)  du2/dy              ) * ( n_y )
      !
      ! Loop over the cubature points on the current element
      ! to assemble the integral
      do iel=1,nelementsInRegion

        ! Get the length of the edge on that element, the normal vector,
        ! and the tangential vector.
        iedge = p_IedgesAtElement(IedgesLocal(iel),IelementsInRegion(iel))
        ivt1 = p_IverticesAtEdge(1,iedge)
        ivt2 = p_IverticesAtEdge(2,iedge)
        Dvt1 = p_DvertexCoords(:,ivt1)
        Dvt2 = p_DvertexCoords(:,ivt2)

        Dtangential = Dvt2(:)-Dvt1(:)

        dedgelen = sqrt(dtangential(1)**2 + dtangential(2)**2)

        Dtangential = Dtangential / dedgelen

        dnormal(1) = -Dtangential(2)
        dnormal(2) = Dtangential(1)

        do icubp = 1,ncubp

          ! Calculate the OMEGA for the integration by multiplication
          ! of the integration coefficient by the Jacobian of the
          ! mapping.
          ! The determinant of the mapping of the unit interval [-1,1]
          ! to the real line is 0.5*length of the line!

          dweight = Domega(icubp)*0.5_DP*dedgelen

          ! Loop through the DOF`s on our element and calculate
          ! the tangential U as well as P.
          dval1 =   2.0_DP*Dvalues(icubp,iel,1) * Dnormal(1) &
                  + Dvalues(icubp,iel,2) * Dnormal(2) &
                  + Dvalues(icubp,iel,3) * Dnormal(2)

          dval2 =   2.0_DP*Dvalues(icubp,iel,4) * Dnormal(2) &
                  + Dvalues(icubp,iel,2) * Dnormal(1) &
                  + Dvalues(icubp,iel,3) * Dnormal(1)

          dpres = Dvalues(icubp,iel,5)

          ! Sum this up to the two integral contributions for the pressure and
          ! velocity.
          DintU(1) = DintU(1) + dweight * dpf1 * dval1
          DintU(2) = DintU(2) + dweight * dpf1 * dval2
          DintP(1) = DintP(1) - dweight * dpres * Dnormal(1)
          DintP(2) = DintP(2) - dweight * dpres * Dnormal(2)

        end do ! icubp

      end do ! iel

    end select

    ! DintU and DintP give now the contributions to the force integral:
    !
    ! DragCoeff = 2/dfp2 * (dpf1*DintU(1) + DintP(1))
    ! LiftCoeff = 2/dfp2 * (dpf1*DintU(2) + DintP(2))

    Dforces = 0.0_DP
    Dforces(1:NDIM2D) = 2.0_DP/dpf2 * (DintU(:) + DintP(:))

    ! Deallocate memory, finish.
    deallocate(Dvalues)
    deallocate(IelementsInRegion)

  end subroutine


  !****************************************************************************

!<subroutine>

  subroutine ppns3D_bdforces_uniform (rvector,rregion,Dforces,ccub,df1,df2)

!<description>
  ! Calculates the drag-/lift-forces acting on a part of the real
  ! boundary for a vector rvector with the solution of the 3D
  ! (Navier-)Stokes equation. It is assumed that
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
  ! <tex>
  !    $$ Dforces(1) = 2/df2 * \int_\Gamma [df1 dut/dn n_y - p n_x] ds $$
  !    $$ Dforces(2) = 2/df2 * \int_\Gamma [df1 dut/dn n_x + p n_y] ds $$
  !    $$ Dforces(3) = 2/df2 * \int_\Gamma [df1 dut/dn n_z + p n_z] ds $$
  ! </tex>
  ! where df1 and df2 allow to weight different integral parts
  !
  ! Usually in benchmark-like geometries there is:
  ! <tex>
  !   $$ df1 = RHO*NU = (\text{density of the fluid}) * \text{viscosity}$$
  !   $$ df2 = RHO*DIST*UMEAN**2
  !          = (\text{density of the fluid} *
  !             \text{length of the obstacle facing the flow})
  !           *(\text{mean velocity of the fluid})^2 $$
  ! </tex>
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
  type(t_vectorBlock), intent(in)     :: rvector

  ! Mesh region where to calculate the boundary forces.
  type(t_meshRegion), intent(in)      :: rregion

  ! 2D Cubature formula identifier to use for the quad integration.
  ! One of the CUB_xxxx constants in the cubature.f90.
  integer(I32), intent(in)            :: ccub

  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(in), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(in), optional      :: df2

!</input>

!<output>
  ! Array receiving the forces acting on the boundary specified by rregion.
  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
  real(DP), dimension(:), intent(out) :: Dforces
!</output>

!</subroutine>

  ! Spatial discretisation structure of velocity and pressure
  type(t_spatialDiscretisation), pointer :: p_rdiscrU, p_rdiscrP

  ! Element type identifier for U and P
  integer(I32) :: ielemU, ielemP

  ! Number of local DOF`s in U and P
  integer :: idoflocU, idoflocP

  ! Triangulation
  type (t_triangulation), pointer :: p_rtria

  ! An accepting the DOF`s of an element.
  integer, dimension(EL_MAXNBAS), target :: IdofsU, IdofsP

  ! Coordinates of the vertices
  real(DP), dimension(NDIM3D,8) :: Dcoords

  ! Coordinates of the normal vectors
  real(DP), dimension(:,:), allocatable :: Dnormal

  ! Coordinates of the cubature points on reference and real element
  real(DP), dimension(NDIM3D,CUB_MAXCUBP_2D) :: DpointsRef,DpointsReal

  ! Coordinate system for U and P element
  integer(I32) :: ctrafoU, ctrafoP

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
    ! velocity and pressure
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
    ! These have to be mapped to the 3D element which is under consideration
    ! during the integration -- later.
    !
    ! Before, we need some additional information about the triangulation
    ! and about our element!
    !
    ! Number of local DOF`s:
    idoflocU = elem_igetNDofLoc(ielemU)
    idoflocP = elem_igetNDofLoc(ielemP)

    ! The triangulation - it is the same for U and P
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

      ! For the integration, we need the global DOF`s on our element
      ! for U and P:
      call dof_locGlobMapping(p_rdiscrU, iel, IdofsU)
      call dof_locGlobMapping(p_rdiscrP, iel, IdofsP)

      ! Calculate the transformation for all points on the current element.
      ! If we have only parametric elements, we do not have to calculate
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

        ! Loop through the DOF`s on our element and calculate
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

    ! That is it

    contains

!<subroutine>

    pure subroutine ppns3D_det(Dcoords,iface,Dcub,ncubp,Ddetj)

!<description>
  ! This routine calculates the "determinant" of a bilinear quadrilateral
  ! transformation from the 2D reference quadrilateral onto a hexahedron
  ! face in 3D.
!</description>

!<input>
    ! The coordinates of the eight corner vertices of the hexahedron
    real(DP), dimension(:,:), intent(in) :: Dcoords

    ! The index of the face onto which the points are mapped
    integer, intent(in) :: iface

    ! The 2D coordinates of the points that are to be mapped
    real(DP), dimension(:,:), intent(in) :: Dcub

    ! The number of points which are to be mapped
    integer, intent(in) :: ncubp
!</input>

!<output>
    ! The jacobian determinants of the mapping
    real(DP), dimension(:), intent(out) :: Ddetj
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

  subroutine ppns2D_streamfct_uniform (rvector,rdestVector,isubvectorVel)

!<description>
  ! Calculates the streamfunction of a 2D velocity field rvector.
  ! It is assumed that
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
  ! <tex>$Q_1$</tex> discretisation!
!</description>

!<input>
  ! The FE solution vector.
  type(t_vectorBlock), intent(in)    :: rvector

  ! OPTIONAL: Defines the start index in rvector where the velocity starts.
  ! Default = 1.
  integer, intent(in), optional :: isubvectorVel
!</input>

!<inputoutput>
  ! An empty vector that is prepared with a discretisation structure to
  ! represent the streamfunction. The values in the vector are overwritten
  ! with the FE representation of the streamfunction.
  type(t_vectorScalar), intent(inout),target    :: rdestVector
!</inputoutput>

!</subroutine>

    ! local variables
    integer, parameter :: NVE = 4

    integer :: iel,ielaux,icurrentelement,ivelstart
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

    ! Let us go. Note that we perform the same computation for all,
    ! parametric and nonparametric, point-based and integral-mean-value
    ! based Q1~ solutions. Gives a slight error, but the final Q1
    ! representation is not exact anyway!
    !
    p_rtriangulation => rdestVector%p_rspatialDiscr%p_rtriangulation
    if (.not. associated(p_rtriangulation)) then
      print *,'ppns2D_bdforces_uniform: Unknown triangulation!'
    end if

    ! Get pointers to the subvectors from the block vector
    ivelstart = 1
    if (present(isubvectorVel)) ivelstart = isubvectorVel
    call lsyssc_getbase_double (rvector%RvectorBlock(ivelstart),p_DdataUX)
    call lsyssc_getbase_double (rvector%RvectorBlock(ivelstart+1),p_DdataUY)
    call lsyssc_getbase_double (rdestVector,p_Dx)

    ! Auxiliary array
    call storage_new ('ppns2D_streamfct_uniform', 'aux', &
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

      ! If all four vertices are marked, there is nothing to calculate
      ! on the element. If no vertex is marked, we can not calculate
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
        ! elements of iel if there is an element where the streamfunction
        ! is not calculated in all vertices...

      end if

      ! Now we make a 'greedy' search from the current element to find
      ! as many other elements as possible where we can calculate the
      ! streamfunction.
      ! Look onto the adjacent elements of the current element if there is
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

    ! At last, normalise the streamfunction such that vertex 1
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
    real(DP), dimension(:,:), intent(in)               :: DvertexCoords

    ! Vertices at the element
    integer, dimension(:,:), intent(in) :: IverticesAtElement

    ! Edges at the element
    integer, dimension(:,:), intent(in) :: IedgesAtElement

    ! Element number where to calculate the streamfunction
    integer, intent(in)               :: iel

    ! Local number (1..NVE) of any of the vertices on the element iel
    ! where the streamfunction is already calculated. This will be used
    ! as 'base' to calculate the others.
    integer, intent(in)                                :: ibaseIdx

    ! Marker array of length NVT. All vertices where streamfunction
    ! values are calculated are marked as 1.
    integer, intent(inout), dimension(:)          :: Imarkers

    ! X/Y-velocity.
    real(DP), dimension(:), intent(in)                 :: Du, Dv

    ! Vector of size NVT; receives the values of the streamfunction
    real(DP), dimension(:), intent(inout)              :: Dx

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
        ! marked neighboured vertices, so there is nothing to gain here.

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
          ! N=(dn1,dn2) - not normalised.
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

  !****************************************************************************

!<subroutine>

  subroutine ppns2D_bdforces_vol(rvector,rcharfct,Dforces,df1,df2,cformulation,&
      ffunctionReference,rcollection,rperfconfig)

!<description>
  ! Calculates the drag-/lift-forces of a function defined by rvector
  ! in a part of the domain identified by the characteristic function rcharfct.
  ! The characteristic function can be calculated e.g. by functions like
  ! anprj_charFctRealBdComp.
  !
  ! It is assumed that
  !   rvector%rvectorBlock(1) = X-velocity,
  !   rvector%rvectorBlock(2) = Y-velocity,
  !   rvector%rvectorBlock(3) = pressure.
  !
  ! rregion specifies a boundary region where to calculate
  ! the force integral. Dforces(1:NDIM2D) receives the forces in the
  ! X- and Y-direction.
  !
  ! The body forces are defined as the integrals
  ! <tex>
  !    $$ Dforces(1) = 2/df2 * \int_\Gamma [df1 dut/dn n_y - p n_x] ds $$
  !    $$ Dforces(2) = 2/df2 * \int_\Gamma [df1 dut/dn n_x + p n_y] ds $$
  ! </tex>
  ! where df1 and df2 allow to weight different integral parts.
  !
  ! Usually in benchmark-like geometries there is:
  ! <tex>
  !   $$ df1 = RHO*NU = (\text{density of the fluid}) * \text{viscosity}$$
  !   $$ df2 = RHO*DIST*UMEAN**2
  !          = (\text{density of the fluid} *
  !             \text{length of the obstacle facing the flow})
  !           *(\text{mean velocity of the fluid})^2 $$
  ! </tex>
  ! which influences this routine to calculate the drag-/lift coefficients
  ! instead of the forces. The actual forces are calculated by
  ! setting df1=1.0 and df2=2.0, which is the standard setting when
  ! neglecting df1 and df2.
  !
  ! Double precision version for all Finite Elements without
  ! curved boundaries. Nonconstant df1 supported.
!</description>

!<input>
  ! The FE solution vector.
  type(t_vectorBlock), intent(in)    :: rvector

  ! A characteristic function of the subdomain inside of the domain
  ! where the integraion should take place.
  ! This is a 0-1 vector, usually defined in the vertices and given
  ! as P1/Q1 vector. Is is =1 in the inner of objects e.g.
  type(t_vectorScalar), intent(in)    :: rcharfct

  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(in), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(in), optional      :: df2

  ! OPTIONAL: Type of tensor to use for the computation.
  ! May be either PPNAVST_GRADIENTTENSOR_XXXX for the gradient tensor or
  ! PPNAVST_GDEFORMATIONTENSOR for the deformation tensor formulation.
  ! If not present, PPNAVST_GRADIENTTENSOR_SIMPLE is assumed.
  integer, intent(in), optional :: cformulation

  ! A callback routine that allows to specify a nonconstant coefficient df1.
  ! If not present, df1 is used as constant coefficient in the integral.
  ! If present, df1 is ignored and calculated by cubature point
  ! by this callback routine.
  include '../Postprocessing/intf_refFunctionSc.inc'
  optional :: ffunctionReference

  ! OPTIONAL: A collection structure that is passed to ffunctionRefSimple.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
  ! Array receiving the forces acting on the boundary specified by rregion.
  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
  real(DP), dimension(:), intent(out) :: Dforces
!</output>

!</subroutine>

    ! local variables
    integer :: i,k,icurrentElementDistr, ICUBP, NVE
    integer :: IEL, IELmax, IELset
    real(DP) :: OM, DN1, DN2, dpp
    real(DP) :: ah1,ah2,du1x,du1y,du2x,du2y,dalx,daly
    real(DP), dimension(2) :: DintU, DintP

    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(CUB_MAXCUBP) :: Domega

    ! number of cubature points on the reference element
    integer :: ncubp

    ! Number of local degees of freedom for test functions
    integer :: indofTrial,indofFunc1,indofFunc2

    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation

    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList

    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj

    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistributionU
    type(t_elementDistribution), pointer :: p_relementDistributionP
    type(t_elementDistribution), pointer :: p_relementDistributionA

    ! Number of elements in the current element distribution
    integer :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients

    ! Number of elements in a block. Normally =NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock

    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_evalElementSet) :: revalElementSet
    type(t_domainIntSubset) :: rintSubset

    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsFunc1
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsFunc2


    ! Type of transformation from the reference to the real element
    integer(I32) :: ctrafoType

    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! Type of formulation
    integer :: cform

    real(dp) :: dpf2
    real(dp), dimension(:,:), pointer :: Dpf1

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => ppns_perfconfig
    end if

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rcharfct%p_rspatialDiscr%p_rtriangulation

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(p_rperfconfig%NELEMSIM,p_rtriangulation%NEL)

    ! Prepare the weighting coefficients
    dpf2 = 2.0_DP
    if (present(df2)) dpf2 = df2

    ! Get the used formulation; gradient or deformation tensor.
    cform = PPNAVST_GRADIENTTENSOR_SIMPLE
    if (present(cformulation)) cform = cformulation

    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1,rvector%p_rblockDiscr%RspatialDiscr(1)%inumFESpaces

      ! Activate the current element distribution
      p_relementDistributionU => &
      rvector%p_rblockDiscr%RspatialDiscr(1)%RelementDistr(icurrentElementDistr)

      p_relementDistributionA =>&
      rcharfct%p_rspatialDiscr%RelementDistr(icurrentElementDistr)

      p_relementDistributionP => &
      rvector%p_rblockDiscr%RspatialDiscr(3)%RelementDistr(icurrentElementDistr)

      ! Cancel if this element distribution is empty.
      if (p_relementDistributionU%NEL .eq. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistributionU%celement)
      indofFunc1 = elem_igetNDofLoc(p_relementDistributionA%celement)
      indofFunc2 = elem_igetNDofLoc(p_relementDistributionP%celement)

      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistributionU%celement)

      if (NVE .ne. elem_igetNVE(p_relementDistributionA%celement)) then
        call output_line ('Element spaces incompatible!', &
                          OU_CLASS_ERROR,OU_MODE_STD, 'ppns2D_bdforces_vol')
        call sys_halt()
      end if

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistributionU%celement)

      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      ! Now Dxi stores the point coordinates of the cubature points on the reference element
      call cub_getCubPoints(p_relementDistributionU%ccubTypeEval, ncubp, Dxi, Domega)

      ! Allocate some memory to hold the cubature points on the reference element
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

      ! Reformat the cubature points; they are in the wrong shape!
      do i=1,ncubp
        do k=1,ubound(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        end do
      end do

      ! Prepare the DF1 array
      allocate(Dpf1(ncubp,nelementsPerBlock))

      if (present(df1)) then
        Dpf1(:,:) = df1
      else
        Dpf1(:,:) = 1.0_DP
      end if

      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))
      allocate(IdofsFunc1(indofFunc1,nelementsPerBlock))
      allocate(IdofsFunc2(indofFunc2,nelementsPerBlock))

      ! Allocate memory for the coefficients
      allocate(Dcoefficients(ncubp,nelementsPerBlock,13))

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistributionU%celement)
      cevaluationTag = ior(cevaluationTag,elem_getEvaluationTag(p_relementDistributionP%celement))
      cevaluationTag = ior(cevaluationTag,elem_getEvaluationTag(p_relementDistributionA%celement))

      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REALPOINTS)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistributionU%h_IelementList, &
                                p_IelementList)

      ! Get the number of elements there.
      NEL = p_relementDistributionU%NEL

      ! We assemble the integral contributions separately
      DintU = 0.0_DP
      DintP = 0.0_DP

      ! Prepare the call to the evaluation routine of the analytic function.
      call elprep_init(revalElementSet)

      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, p_rperfconfig%NELEMSIM

        ! We always handle NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most NELEMSIM
        ! elements simultaneously.

        IELmax = min(NEL,IELset-1+p_rperfconfig%NELEMSIM)

        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our NELEMSIM elements simultaneously.

        !--------------------------------------------------------------------------------
        call dof_locGlobMapping_mult(rvector%p_rblockDiscr%RspatialDiscr(1), &
                                     p_IelementList(IELset:IELmax),IdofsTrial)
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        call dof_locGlobMapping_mult(rcharfct%p_rspatialDiscr, &
                                     p_IelementList(IELset:IELmax),IdofsFunc1)
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        call dof_locGlobMapping_mult(rvector%p_rblockDiscr%RspatialDiscr(3), &
                                     p_IelementList(IELset:IELmax),IdofsFunc2)
        !--------------------------------------------------------------------------------
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
        p_Ddetj => revalElementSet%p_Ddetj

        ! In case, dpf1 is nonconstant, calculate dpf1.
        if (present(ffunctionReference)) then
          call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
          !rintSubset%ielementDistribution = icurrentElementDistr
          rintSubset%ielementStartIdx = 1
          rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
          rintSubset%p_IdofsTrial => IdofsTrial
          rintSubset%celement = p_relementDistributionU%celement
          call ffunctionReference (DER_FUNC,rvector%p_rblockDiscr%RspatialDiscr(1), &
                    IELmax-IELset+1,ncubp,revalElementSet%p_DpointsReal, &
                    IdofsTrial,rintSubset,Dpf1(:,:),rcollection)
          call domint_doneIntegration (rintSubset)
        end if

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!

        !----------------------------------------------------------------------------
        !                         EVALUATION PHASE
        !----------------------------------------------------------------------------
        !
        ! We need to build the following system:
        !    /                                                              \
        !   | |-p(x_i)             |   |du1/dx (x_i) du1/dy (x_i) ...|       |   | -dalpha/dx (x_i) |
        !   | |       -p(x_i)      | + |...  you know this works  ...| + u^t | * | -dalpha/dy (x_i) |
        !    \                                                              /
        !
        !
        ! Get the pressure in the cubature points
        ! Save the result to Dcoefficients(:,:,1)

        ! Build the p matrix
        call fevl_evaluate_sim3 (rvector%RvectorBlock(3), revalElementSet, &
                p_relementDistributionP%celement, &
                IdofsFunc2, DER_FUNC, Dcoefficients(:,1:IELmax-IELset+1_I32,1))

        ! Build the jacobi matrix of this (u1,u2,u3)
        ! First Row -------------------------------
        ! Save the result to Dcoefficients(:,:,2:4)
        call fevl_evaluate_sim3 (rvector%RvectorBlock(1), revalElementSet, &
                p_relementDistributionU%celement, &
                IdofsTrial, DER_DERIV2D_X, Dcoefficients(:,1:IELmax-IELset+1_I32,2))

        call fevl_evaluate_sim3 (rvector%RvectorBlock(1), revalElementSet, &
                p_relementDistributionU%celement, &
                IdofsTrial, DER_DERIV2D_Y, Dcoefficients(:,1:IELmax-IELset+1_I32,3))

        ! Second Row -------------------------------
        ! Save the result to Dcoefficients(:,:,4:5)
        call fevl_evaluate_sim3 (rvector%RvectorBlock(2), revalElementSet, &
                p_relementDistributionU%celement, &
                IdofsTrial, DER_DERIV2D_X, Dcoefficients(:,1:IELmax-IELset+1_I32,4))

        call fevl_evaluate_sim3 (rvector%RvectorBlock(2), revalElementSet, &
                p_relementDistributionU%celement, &
                IdofsTrial, DER_DERIV2D_Y, Dcoefficients(:,1:IELmax-IELset+1_I32,5))

        ! Build the alpha vector
        ! Save the result to Dcoefficients(:,:,6:7)
        call fevl_evaluate_sim3 (rcharfct, revalElementSet,&
                p_relementDistributionA%celement, IdofsFunc1, DER_DERIV2D_X,&
                Dcoefficients(:,1:IELmax-IELset+1_I32,6))

        call fevl_evaluate_sim3 (rcharfct, revalElementSet,&
                p_relementDistributionA%celement, IdofsFunc1, DER_DERIV2D_Y,&
                Dcoefficients(:,1:IELmax-IELset+1_I32,7))

        select case (cform)
        case (PPNAVST_GRADIENTTENSOR_SIMPLE,PPNAVST_GRADIENTTENSOR)
          ! 'Full' Gradient tensor formulation.
          ! The 'simple' formulation is not available we do the 'full' here.
          !
          ! We need to calculate the integral based on the following integrand:
          !   |-p         |   (du1/dx  du1/dy )   ( n_x )
          !   |       -p  | + (du2/dx  du2/dy ) * ( n_y )
          !
          ! Loop over the cubature points on the current element
          ! to assemble the integral

          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (-p * I + Dj(u)) * (-grad(alpha)) dx
          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp

              OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

              ! get the pressure
              dpp  = Dcoefficients(icubp,iel,1)

              ! x and y derivative of u1
              du1x = Dcoefficients(icubp,iel,2)
              du1y = Dcoefficients(icubp,iel,3)

              ! x and y derivative of u2
              du2x = Dcoefficients(icubp,iel,4)
              du2y = Dcoefficients(icubp,iel,5)

              dalx = Dcoefficients(icubp,iel,6)
              daly = Dcoefficients(icubp,iel,7)

              ! Derivatice of the characteristic function is the (scaled) normal
              dn1  = -dalx
              dn2  = -daly

              ! gradient tensor
              ah1 = du1x*dn1 + du1y*dn2
              ah2 = du2x*dn1 + du2y*dn2

              ! Sum this up to the two integral contributions for the pressure and
              ! then velocity.
              DintU(1) = DintU(1) + om * Dpf1(icubp,iel) * ah1
              DintU(2) = DintU(2) + om * Dpf1(icubp,iel) * ah2
              DintP(1) = DintP(1) - om * dpp * dn1
              DintP(2) = DintP(2) - om * dpp * dn2

            end do ! ICUBP

          end do ! IEL

        case (PPNAVST_DEFORMATIONTENSOR)

          ! Deformation tensor formulation
          !
          ! We need to calculate the integral based on the following integrand:
          !   |-p         |   (du1/dx               1/2 (du1/dy+du2/dx) )   ( n_x )
          !   |       -p  | + (1/2 (du1/dy+du2/dx)  du2/dy              ) * ( n_y )
          !
          ! Loop over the cubature points on the current element
          ! to assemble the integral

          ! Loop through elements in the set and for each element,
          ! loop through the DOF's and cubature points to calculate the
          ! integral: int_Omega (-p * I + Dj(u)) * (-grad(alpha)) dx
          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp

              OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

              ! get the pressure
              dpp  = Dcoefficients(icubp,iel,1)

              ! x and y derivative of u1
              du1x = Dcoefficients(icubp,iel,2)
              du1y = Dcoefficients(icubp,iel,3)

              ! x and y derivative of u2
              du2x = Dcoefficients(icubp,iel,4)
              du2y = Dcoefficients(icubp,iel,5)

              dalx = Dcoefficients(icubp,iel,6)
              daly = Dcoefficients(icubp,iel,7)

              ! Derivatice of the characteristic function is  the (scaled) normal
              dn1  = -dalx
              dn2  = -daly

              ! deformation tensor
              ah1 = 2.0_dp*du1x*dn1 + (du1y+du2x)*dn2
              ah2 = (du1y+du2x)*dn1 + 2.0_dp*du2y*dn2

              ! Sum this up to the two integral contributions for the pressure and
              ! the velocity.
              DintU(1) = DintU(1) + om * Dpf1(icubp,iel) * ah1
              DintU(2) = DintU(2) + om * Dpf1(icubp,iel) * ah2
              DintP(1) = DintP(1) - om * dpp * dn1
              DintP(2) = DintP(2) - om * dpp * dn2

            end do ! ICUBP

          end do ! IEL

        end select

      end do ! IELset

      ! Release memory
      call elprep_releaseElementSet(revalElementSet)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)
      deallocate(IdofsFunc1)
      deallocate(IdofsFunc2)
      deallocate(Dpf1)

    end do ! icurrentElementDistr

    Dforces = 0.0_DP
    Dforces(1:NDIM2D) = 2.0_DP/dpf2 * (DintU(:) + DintP(:))

  end subroutine

  ! ***************************************************************************

  subroutine ffluxCallback (nelements,npointsPerElement,Ielements,&
                                  Dpoints,Dvalues,rcollection,&
                                  DpointsRef,Djac,Ddetj)

  ! Calculates the flux in a set of points on a set of elements.

  integer, intent(in)                               :: nelements
  integer, intent(in)                               :: npointsPerElement
  integer, dimension(:), intent(in)                 :: Ielements
  real(DP), dimension(:,:,:), intent(in)            :: Dpoints
  type(t_collection), intent(inout), optional       :: rcollection
  real(DP), dimension(:,:,:),intent(in), optional   :: DpointsRef
  real(DP), dimension(:,:,:),intent(in), optional   :: Djac
  real(DP), dimension(:,:),intent(in), optional     :: Ddetj

  real(DP), dimension(:,:), intent(out)             :: Dvalues

    ! local variables
    real(DP), dimension(:,:,:), allocatable :: Dvelocity
    type(t_vectorBlock), pointer :: p_rvector
    real(DP) :: dn1,dn2
    integer :: iel,ipt

    ! Calculate the velocity in all the points.
    allocate(Dvelocity(npointsPerElement,nelements,2))

    ! For that purpose, take the velocity vector from the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1

    ! Calculate the X-velocity
    call fevl_evaluate_sim (DER_FUNC, Dvelocity(:,:,1), p_rvector%RvectorBlock(1), &
        Dpoints, Ielements)

    ! Calculate the Y-velocity
    call fevl_evaluate_sim (DER_FUNC, Dvelocity(:,:,2), p_rvector%RvectorBlock(2), &
        Dpoints, Ielements)

    ! Get the normal of the line
    dn1 = rcollection%DquickAccess(1)
    dn2 = rcollection%DquickAccess(2)

    ! Calculate the flux
    do iel = 1,nelements
      do ipt = 1,npointsPerElement
        Dvalues(ipt,iel) = Dvelocity(ipt,iel,1)*dn1 + Dvelocity(ipt,iel,2)*dn2
      end do
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine ppns2D_calcFluxThroughLine (rvector,Dstart,Dend,dflux,&
      ccubature,nlevels)

!<description>
  ! Calculates the flux through a line given by two points.
!</description>

!<input>
  ! Velocity vector.
  type(t_vectorBlock), intent(in), target :: rvector

  ! Starting point of the line.
  real(dp), dimension(:), intent(in) :: Dstart

  ! Ending point of the line.
  real(dp), dimension(:), intent(in) :: Dend

  ! OPTIONAL: Basic 1D cubature rule tu use for line integration.
  integer(I32), intent(in), optional :: ccubature

  ! OPTIONAL: Refinement of the cubature rule. >= 0.
  ! The line Dstart..Dend will be divided into 2**idegree intervals and a
  ! summed cubature formula will be applied.
  integer, intent(in), optional :: nlevels

!</input>

!<output>
  ! The calculated flux.
  real(DP), intent(out) :: dflux
!</output>

!</subroutine>

    type(t_collection) :: rcollection
    real(dp) :: dnorm,dsize,dlocalh
    real(DP), dimension(NDIM2D) :: Dbox
    integer :: ireflevel,nel
    integer(I32) :: ccub

    ! Default settings
    ccub = CUB_G3_1D
    if (present(ccubature)) ccub = ccubature

    ! Initialise the collection for the calculation
    rcollection%p_rvectorQuickAccess1 => rvector

    if (.not. present(nlevels)) then
      ! Normal vector
      dnorm = sqrt((Dend(1)-Dstart(1))**2 + (Dend(2)-Dstart(2))**2)
      if (dnorm .eq. 0.0_DP) then
        ! No line, no flux...
        dflux = 0
        return
      end if
      rcollection%DquickAccess(1) = -(Dend(2)-Dstart(2)) / dnorm
      rcollection%DquickAccess(2) = (Dend(1)-Dstart(1)) / dnorm

      ! To compute the level of refinement of the cubature rule, we make a guess.
      Dbox(1:2) = rvector%p_rblockDiscr%p_rtriangulation%DboundingBoxMax(1:2) &
          - rvector%p_rblockDiscr%p_rtriangulation%DboundingBoxMin(1:2)
      dsize = min(Dbox(1),Dbox(2))

      ! In the mean, an element has this size:
      dlocalh = dsize / sqrt(real(rvector%p_rblockDiscr%p_rtriangulation%NEL,dp))

      ! Divide the length of the line by this, that's the mean number of elements
      ! on the line.
      nel = min(int(dnorm/dlocalh),1)

      ! Every refinement gives 2 elements, so the log2 of nel is the cubature
      ! refinement. We add 2 levels as 'safety' buffer.
      ireflevel = log(real(nel,dp))/log(2._DP) + 2
    else
      ! Take the given one.
      ireflevel = nlevels
    end if

    ! Go...
    call ppint_lineIntegral (dflux,rvector%p_rblockDiscr%p_rtriangulation,&
        Dstart,Dend,ccub,ireflevel,ffluxCallback,rcollection)

  end subroutine

end module
