!#########################################################################
!# ***********************************************************************
!# <name> boundaryintegral </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module allows to calculate integral values over the boundary of
!# a domain. The following routines can be found here:
!#
!# 1.) bdint_integral2D
!#     -> Calculates the boundary integral
!#            <tex> $$ \int_{\Gamma} f(x) dx $$ </tex>
!#        for an arbitrary function <tex>$ f(x) $</tex>.
!#
!# 2.) bdint_scalarBoundaryInt2D
!#     -> Calculates the boundary integral
!#            <tex> $$ \int_{\Gamma} f_h(x) dx $$ </tex>
!#        for an arbitrary scalar FE function <tex>$ f_h(x) $</tex>.
!#
!# 3.) bdint_normalDerivativeInt2D
!#     -> Calculates the integral of the normal derivative
!#            <tex> $$ \int_{\Gamma} \partial_n f_h(x) dx $$ </tex>
!#        for an arbitrary scalar FE function <tex>$ f_h(x) $</tex>.
!# </purpose>
!#########################################################################

module boundaryintegral

!$use omp_lib
  use basicgeometry
  use boundary
  use collection
  use cubature
  use derivatives
  use domainintegration
  use feevaluation
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use mprimitives
  use spatialdiscretisation
  use storage
  use transformation
  use triangulation

  implicit none

  private
  public :: bdint_integral2D
  public :: bdint_integral2D_conf
  public :: bdint_scalarBoundaryInt2D
  public :: bdint_normalDerivativeInt2D

contains

  !****************************************************************************

!<subroutine>

  subroutine bdint_integral2D (rboundary, rtriangulation, ccubType,&
      ffunction, dvalue, rboundaryRegion, rcollection)

!<description>
  ! This routine calculates a boundary integral of an arbitrary
  ! scalar function in 2D specified by the callback routine ffunction.
  !
  ! rboundaryRegion is a t_boundaryRegion object that allows to
  ! specify the boundary region where the integral should be computed.
  ! If not specified, the integral is computed over the whole boundary.
  !
  ! NOTE: The routine is SLOW but general!
!</description>

!<input>
  ! A boundary object that specifies the analytical boundary of the domain.
  type(t_boundary), intent(in) :: rboundary

  ! A triangulation object that specifies the mesh in the domain.
  type(t_triangulation), intent(in) :: rtriangulation

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer(I32), intent(in) :: ccubType

  ! A callback function that provides the analytical reference
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_functionScBdr2D.inc'

  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion

  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function to provide additional information.
  type(t_collection), intent(in), target, optional :: rcollection
!</input>

!<output>
  ! The calculated value of the integral.
  real(DP), intent(out) :: dvalue
!</output>

!</subroutine>

    ! local variables
    type(t_collection), pointer :: p_rcollection
    type(t_boundaryRegion) :: rboundaryReg
    real(DP) :: dlocalValue
    integer :: ibdc

    ! Let p_rcollection point to rcollection - or NULL if it is not
    ! given.
    if (present(rcollection)) then
      p_rcollection => rcollection
    else
      p_rcollection => null()
    end if

    if (rtriangulation%ndim .ne. NDIM2D) then
      call output_line('Only 2D discretisations allowed!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdint_integral2D')
      call sys_halt()
    end if

    ! If the boundary region is specified, call pperr_scalarBoundary2d_conf
    ! for that boundary region. Otherwise, call pperr_scalarBoundary2d_conf
    ! for all possible boundary regions and sum up the errors.
    if (present(rboundaryRegion)) then
      call bdint_integral2D_conf (rboundary, rtriangulation, ccubType,&
        ffunction, dvalue, rboundaryRegion, p_rcollection)
    else
      dvalue = 0.0_DP
      ! Create a boundary region for each boundary component and call
      ! the calculation routine for that.
      do ibdc = 1, boundary_igetNBoundComp(rboundary)
        call boundary_createRegion (rboundary, ibdc, 0, rboundaryReg)
        call bdint_integral2D_conf (rboundary, rtriangulation, ccubType,&
          ffunction, dlocalValue, rboundaryReg, p_rcollection)
        dvalue = dvalue + dlocalValue
      end do
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bdint_integral2D_conf (rboundary, rtriangulation, ccubType,&
      ffunction, dvalue, rboundaryRegion, rcollection)

!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

!<input>
  ! A boundary object that specifies the analytical boundary of the domain.
  type(t_boundary), intent(in) :: rboundary

  ! A triangulation object that specifies the mesh in the domain.
  type(t_triangulation), intent(in) :: rtriangulation

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer(I32), intent(in)                     :: ccubType

  ! A t_boundaryRegion specifying the boundary region where
  ! to calculate.
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion

  ! A callback function that provides the analytical reference
  ! function to which the error should be computed.
  ! If not specified, the reference function is assumed to be zero!
  include 'intf_functionScBdr2D.inc'

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection
!</input>

!<output>
  ! The calculated value of the integral.
  real(DP), intent(out) :: dvalue
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), allocatable :: IelementOrientation
    integer, dimension(:), allocatable :: Ielements
    real(DP), dimension(:,:), allocatable :: DedgePosition

    integer :: ibdc,ibdcoffset,iedge,ilocaledge,nve
    integer :: NEL,NELbdc,iel,i,k
    integer(I32) :: ctrafoType

    ! The triangulation structure - to shorten some things...
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:), pointer :: p_IelementsAtBoundary
    real(DP), dimension(:), pointer :: p_DedgeParameterValue
    real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Arrays for cubature points
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D
    real(DP), dimension(:,:,:), allocatable :: Dxi2D,Dpoints,DpointsRef
    real(DP), dimension(CUB_MAXCUBP) :: Domega1D
    real(DP), dimension(:,:), allocatable :: Dvalues
    real(DP), dimension(NDIM2D,TRIA_MAXNVE) :: Dcoord
    integer :: ncubp,ipoint
    integer(I32) :: icoordSystem
    real(DP) :: dlen,dpar1,dpar2

    real(DP), dimension(:,:), allocatable :: DpointPar

    ! Get some pointers and arrays for quicker access

    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_int (rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    call storage_getbase_int (rtriangulation%h_IelementsAtBoundary,&
        p_IelementsAtBoundary)
    call storage_getbase_int2d (rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,&
        p_DvertexCoordinates)
    call storage_getbase_double (rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)
    call storage_getbase_double (rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)

    ! Boundary component?
    ibdc = rboundaryRegion%iboundCompIdx

    ! Number of elements on that boundary component?
    NELbdc = p_IboundaryCpIdx(ibdc+1)-p_IboundaryCpIdx(ibdc)

    ! Position of the boundary component?
    ibdcoffset = p_IboundaryCpIdx(ibdc)

    ! In a first step, we figure out the elements on the boundary and their
    ! orientation. Allocate arrays that are large enough to hold
    ! even all elements on the boundary if necessary.
    allocate(Ielements(NELbdc), IelementOrientation(NELbdc))

    ! Allocate an array saving the start- and end-parameter values
    ! of the edges on the boundary.
    allocate(DedgePosition(2,NELbdc))

    ! Loop through the edges on the boundary component ibdc.
    ! If the edge is inside, remember the element number and figure out
    ! the orientation of the edge.
    ! NEL counts the total number of elements in the region.
    NEL = 0
    do iedge = 1,NELbdc
      if (boundary_isInRegion(rboundaryRegion,ibdc,&
          p_DedgeParameterValue(iedge))) then
        NEL = NEL + 1

        ! Element number
        Ielements(NEL) = p_IelementsAtBoundary(iedge)

        ! Element orientation; i.e. the local number of the boundary edge
        do ilocaledge = 1,ubound(p_IedgesAtElement,1)
          if (p_IedgesAtElement(ilocaledge,p_IelementsAtBoundary(iedge)) .eq. &
              p_IedgesAtBoundary(iedge)) exit
        end do
        IelementOrientation(NEL) = ilocaledge

        ! Save the start parameter value of the edge -- in length
        ! parametrisation.
        dpar1 = boundary_convertParameter(rboundary, &
            ibdc, p_DvertexParameterValue(iedge), rboundaryRegion%cparType, &
            BDR_PAR_LENGTH)

        ! Save the end parameter value. Be careful: The last edge
        ! must be treated differently!
        if (iedge .ne. NELbdc) then
          dpar2 = boundary_convertParameter(rboundary, &
              ibdc, p_DvertexParameterValue(iedge+1), rboundaryRegion%cparType, &
              BDR_PAR_LENGTH)

        else
          dpar2 = boundary_dgetMaxParVal(rboundary,ibdc,BDR_PAR_LENGTH)
        end if

        DedgePosition(1,NEL) = dpar1
        DedgePosition(2,NEL) = dpar2

      end if
    end do

    ! Get the parameter values of the 1D cubature formula
    ! as well as the number of cubature points ncubp
    call cub_getCubPoints(ccubType, ncubp, Dxi1D, Domega1D)

    ! Map the 1D cubature points to the edges in 2D.
    allocate(Dxi2D(ncubp,NDIM2D+1,NEL))

    ! The type of the coordinate system may change with every element.
    ! So we may have to switch... ItrialElements in the discretisation
    ! structure informs us about the element type.
    do iel = 1,NEL

      ! How many vertices are on the current element?
      nve = tria_getNVE(p_IverticesAtElement, iel)

      select case (nve)
      case (TRIA_NVETRI2D)
        ! Coordinates of cubature points on triangles are given in
        ! barycentric coordinates.
        icoordSystem = TRAFO_CS_BARY2DTRI
      case (TRIA_NVEQUAD2D)
        ! Coordinates of cubature points on triangles are given in
        ! standard coordinates on the reference element
        icoordSystem = TRAFO_CS_REF2DQUAD
      end select

      call trafo_mapCubPts1Dto2D(icoordSystem, IelementOrientation(iel), &
          ncubp, Dxi1D, Dxi2D(:,:,iel))
    end do

    ! Transpose the coordinate array such that we get coordinates
    ! we can work with.
    allocate(DpointsRef(NDIM2D+1,ncubp,NEL))
    do iel = 1,NEL
      do i = 1,ncubp
        do k =1,ubound(DpointsRef,1)
          DpointsRef(k,i,iel) = Dxi2D(i,k,iel)
        end do
      end do
    end do

    ! Dxi2D is not needed anymore.
    deallocate(Dxi2D)

    ! Calculate the coordinates of the points on world coordinates
    allocate(Dpoints(NDIM2D+1,ncubp,NEL))

    ! Transformation can be different for all elements
    do iel = 1,NEL

      ! Get the points forming the element
      do ipoint = 1,ubound(p_IverticesAtElement,1)
        Dcoord(1,ipoint) = &
            p_DvertexCoordinates(1,p_IverticesAtElement(ipoint,Ielements(iel)))
        Dcoord(2,ipoint) = &
            p_DvertexCoordinates(2,p_IverticesAtElement(ipoint,Ielements(iel)))
      end do

      ! Transform the cubature points
      do ipoint = 1,ncubp

        ! How many vertices are on the current element?
        nve = tria_getNVE(p_IverticesAtElement, iel)

        select case (nve)
        case (TRIA_NVETRI2D)
          ! Coordinates of cubature points on triangles are given in
          ! barycentric coordinates, 2D.
          ctrafoType = TRAFO_ID_LINSIMPLEX + TRAFO_DIM_2D
        case (TRIA_NVEQUAD2D)
          ! Coordinates of cubature points on triangles are given in
          ! standard coordinates on the reference element, 2D.
          ctrafoType = TRAFO_ID_MLINCUBE + TRAFO_DIM_2D
        end select

        call trafo_calcRealCoords (ctrafoType,Dcoord,&
            DpointsRef(:,ipoint,iel),Dpoints(:,ipoint,iel))
      end do
    end do

    ! Calculate the parameter values of the points on the boundary
    allocate (DpointPar(ncubp,NEL))
    do iel=1,NEL
      do ipoint = 1,ncubp
        ! Dxi1D is in [-1,1] while the current edge has parmeter values
        ! [DedgePosition(1),DedgePosition(2)]. So do a linear
        ! transformation to transform Dxi1D into that interval, this
        ! gives the parameter values in length parametrisation
        call mprim_linearRescale(Dxi1D(ipoint,1),-1.0_DP,1.0_DP,&
            DedgePosition(1,iel),DedgePosition(2,iel),DpointPar(ipoint,iel))
      end do
    end do

    ! So Dxi2 defines the coordinates on the reference element for all
    ! elements. Generally said, we have to evaluate the function in these
    ! points now. For that purpose, we call ffunction!

    allocate (Dvalues(ncubp,NEL))
    Dvalues = 0.0_DP

    ! Evaluate the function on the boundary
    call ffunction (DpointsRef, Dpoints, rboundaryRegion%iboundCompIdx,&
        DpointPar,Ielements(1:NEL), Dvalues(:,:),rcollection)

    ! Now, Dvalues1 contains in Dvalues1(:,:,1) the term
    ! "u(x,y)-u_h(x,y)" -- in every cubature point on every
    ! element. We are finally able to calculate the integral!
    ! That means, run over all the edges and sum up...
    ! (ok, if rvectorScalar is not specified, we have
    !  -u_h(x,y) in Dvalues1(:,:,1), but as we take the square,
    !  it does not matter if we have u_h or -u_h there!)

    dvalue = 0.0_DP
    do iel = 1,NEL

      ! Get the length of the edge. Let us use the parameter values
      ! on the boundary for that purpose; this is a more general
      ! implementation than using simple lines as it will later
      ! support isoparametric elements.
      !
      ! The length of the current edge serves as a "determinant"
      ! in the cubature, so we have to divide it by 2 as an edge on
      ! the unit interval [-1,1] has length 2.
      dlen = 0.5_DP*(DedgePosition(2,iel)-DedgePosition(1,iel))

      do ipoint = 1,ncubp
        dvalue = dvalue + dlen * Domega1D(ipoint) * (Dvalues(ipoint,iel)**2)
      end do
    end do

    deallocate(DpointPar)
    deallocate(Dvalues)

    ! Release memory

    deallocate(Dpoints)

    deallocate(DedgePosition)
    deallocate(Ielements, IelementOrientation)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bdint_scalarBoundaryInt2D (rvectorScalar, ccubType,&
      dvalue, rboundaryRegion)

!<description>
  ! This function calculates the boundary integral
  !    <tex> $$ int_{\Gamma} f_h(x) dx $$ </tex>
  ! for the scalar FE function $f_h$ specified in rvectorScalar.
  ! rboundaryRegion allows to specify a region on the boundary where
  ! the integral is computed; if not specified, the integral is
  ! computed over the whole boundary.
!</description>

!<input>
  ! A scalar FE function to compute a boundary integral of.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer(I32), intent(in) :: ccubType

  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion
!</input>

!<output>
  ! The calculated value of the integral.
  real(DP), intent(out) :: dvalue
!</output>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscr
    type(t_collection) :: rcollection

    ! Get the underlying discretisation
    p_rdiscr => rvectorScalar%p_rspatialDiscr

    if (p_rdiscr%ndimension .ne. NDIM2D) then
      call output_line('Only 2D discretisations allowed!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdint_scalarBoundaryInt2D')
      call sys_halt()
    end if

    ! Create a collection that allows us to pass the vector to
    ! the callback routine
    call collct_init(rcollection)

    ! Save the vector to the collection, so we can restore it in
    ! the callback routine
    call collct_setvalue_vecsca (rcollection, 'vector', rvectorScalar, .true.)

    call bdint_integral2D (p_rdiscr%p_rboundary, p_rdiscr%p_rtriangulation,&
      ccubType, ffunctionFEevaluation2D, dvalue, rboundaryRegion, rcollection)

    ! Release the collection
    call collct_deletevalue (rcollection,'vector')
    call collct_done(rcollection)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine ffunctionFEevaluation2D (DpointsRef, Dpoints, ibct, DpointPar, Ielements, &
        Dvalues, rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration

  !<description>
    ! This subroutine is called during the calculation of boundary
    ! integrals. It has to compute the (analytical) values of a
    ! function in a couple of points Dpoints on a couple of elements
    ! Ielements. The integral evaluation function uses these values
    ! to compute integrals.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real and reference coordinates.
    ! It has to to simultaneously compute the desired values for all these points.
  !</description>

  !<input>
    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates correspond to the reference
    ! element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: DpointsRef

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    integer, dimension(:), intent(in) :: Ielements

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
  !</input>

  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out)                      :: Dvalues
  !</output>

  !</subroutine>

    ! local variables
    type(t_vectorScalar), pointer :: p_rvector
    integer :: iel

    ! Get the vector with the FE function from the collection
    p_rvector => collct_getvalue_vecsca (rcollection, 'vector')

    ! Evaluate the FE function in the given points.
    do iel=1,size(Ielements)
      call fevl_evaluate_mult (DER_FUNC, Dvalues(:,iel), p_rvector, &
          Ielements(iel), DpointsRef(:,:,iel), Dpoints(:,:,iel))
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bdint_normalDerivativeInt2D (rvectorScalar, ccubType,&
      dvalue, rboundaryRegion)

!<description>
  ! This function calculates the boundary integral
  !    <tex> $$ int_{\Gamma} grad(f_h(x))*n dx $$ </tex>
  ! for the block FE function $f_h$ specified in rvectorScalar.
  ! rboundaryRegion allows to specify a region on the boundary where
  ! the integral is computed; if not specified, the integral is
  ! computed over the whole boundary.
!</description>

!<input>
  ! A scalar FE function to compute a boundary integral of.
  type(t_vectorScalar), intent(in) :: rvectorScalar

  ! A line cubature formula CUB_xxxx_1D to be used for line integration.
  integer(I32), intent(in) :: ccubType

  ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
  ! to calculate. If not specified, the computation is done over
  ! the whole boundary.
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion
!</input>

!<output>
  ! The calculated value of the integral.
  real(DP), intent(out) :: dvalue
!</output>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscr
    type(t_collection) :: rcollection

    ! Get the underlying discretisation
    p_rdiscr => rvectorScalar%p_rspatialDiscr

    if (p_rdiscr%ndimension .ne. NDIM2D) then
      call output_line('Only 2D discretisations allowed!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdint_normalDerivativeInt2D')
      call sys_halt()
    end if

    ! Create a collection that allows us to pass the vector to
    ! the callback routine
    call collct_init(rcollection)

    ! Save the vector to the collection, so we can restore it in
    ! the callback routine
    call collct_setvalue_vecsca (rcollection, 'vector', rvectorScalar, .true.)

    call bdint_integral2D (p_rdiscr%p_rboundary, p_rdiscr%p_rtriangulation,&
      ccubType, ffunctionNormalDeriv2D, dvalue, rboundaryRegion, rcollection)

    ! Release the collection
    call collct_deletevalue (rcollection,'vector')
    call collct_done(rcollection)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine ffunctionNormalDeriv2D (DpointsRef, Dpoints, ibct, DpointPar, Ielements, &
        Dvalues, rcollection)

    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration

  !<description>
    ! This subroutine is called during the calculation of boundary
    ! integrals. It has to compute the (analytical) values of a
    ! function in a couple of points Dpoints on a couple of elements
    ! Ielements. The integral evaluation function uses these values
    ! to compute integrals.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real and reference coordinates.
    ! It has to to simultaneously compute the desired values for all these points.
  !</description>

  !<input>
    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates correspond to the reference
    ! element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: DpointsRef

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    integer, dimension(:), intent(in) :: Ielements

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
  !</input>

  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
  !</output>

  !</subroutine>

    ! local variables
    type(t_vectorScalar), pointer :: p_rvector
    integer :: iel
    integer :: ipoint
    real(DP), dimension(:,:), allocatable :: DderivX,DderivY
    real(DP) :: dnx,dny,dminPar,dmaxPar,dt

    ! Get the vector with the FE function from the collection
    p_rvector => collct_getvalue_vecsca (rcollection, 'vector')

    ! Allocate memory for the values of the derivative
    allocate(DderivX(ubound(Dvalues,1),ubound(Dvalues,2)))
    allocate(DderivY(ubound(Dvalues,1),ubound(Dvalues,2)))

    ! Evaluate the derivative of the FE function in the given points.
    do iel=1,size(Ielements)
      call fevl_evaluate_mult (DER_DERIV_X, DderivX(:,iel), p_rvector, &
          Ielements(iel), DpointsRef(:,:,iel), Dpoints(:,:,iel))
      call fevl_evaluate_mult (DER_DERIV_Y, DderivY(:,iel), p_rvector, &
          Ielements(iel), DpointsRef(:,:,iel), Dpoints(:,:,iel))
    end do

    ! Get the minimum and maximum parameter value. The point with the minimal
    ! parameter value is the start point of the interval, the point with the
    ! maximum parameter value the endpoint.
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    do iel=1,size(Ielements)
      do ipoint = 1,ubound(Dpoints,1)
        dminPar = min(DpointPar(ipoint,iel),dminPar)
        dmaxPar = max(DpointPar(ipoint,iel),dmaxPar)
      end do
    end do

    ! Multiply the derivative with the normal in each point to get the
    ! normal derivative.
    do iel=1,size(Ielements)
      do ipoint = 1,ubound(Dpoints,1)

        dt = DpointPar(ipoint,iel)

        ! Get the normal vector in the point from the boundary.
        ! Note that the parameter value is in length parametrisation!
        ! When we are at the left or right endpoint of the interval, we
        ! calculate the normal vector based on the current edge.
        ! Without that, the behaviour of the routine may lead to some
        ! confusion if the endpoints of the interval coincide with
        ! the endpoints of a boundary edge. In such a case, the routine
        ! would normally compute the normal vector as a mean on the
        ! normal vectors of the edges adjacent to such a point!
        if (DpointPar(ipoint,iel) .eq. dminPar) then
          ! Start point
          call boundary_getNormalVec2D(p_rvector%p_rspatialDiscr%p_rboundary, ibct, &
              dt, dnx, dny,BDR_NORMAL_RIGHT,BDR_PAR_LENGTH)
        else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
          ! End point
          call boundary_getNormalVec2D(p_rvector%p_rspatialDiscr%p_rboundary, ibct, &
              dt, dnx, dny,BDR_NORMAL_LEFT,BDR_PAR_LENGTH)
        else
          ! Inner point
          call boundary_getNormalVec2D(p_rvector%p_rspatialDiscr%p_rboundary, ibct, &
              dt, dnx, dny,cparType=BDR_PAR_LENGTH)
        end if

        ! Calculate the normal derivative
        Dvalues(ipoint,iel) = DderivX(ipoint,iel)*dnx + DderivY(ipoint,iel)*dny

      end do
    end do

    ! Release memory
    deallocate(DderivX,DderivY)

  end subroutine

end module
