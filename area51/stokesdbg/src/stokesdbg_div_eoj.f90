
module stokesdbg_div_eoj

!$use omp_lib
  use basicgeometry
  use bilinearformevaluation
  use cubature
  use derivatives
  use dofmapping
  use domainintegration
  use element
  use elementpreprocessing
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use perfconfig
  use spatialdiscretisation
  use storage
  use transformation
  use triangulation

  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine stokesdbg_divEoj ( &
      rmatrixScalar,dgamma,dgammastar,deojEdgeExp,dtheta,ccubType,dnu,&
      rdiscretisation,InodeList,rperfconfig)

!<description>
  ! Unified edge oriented jump stabilisation.
  !
  ! Adds the unified edge oriented jump stabilisation to the matrix rmatrix:
  ! <tex>
  ! $$< Ju,v > = \sum_E \max(\gamma^{*} * \nu * h_E, \gamma h_E^2)
  !            \int_E [grad u] [grad v] ds $$
  ! </tex>
  ! Uniform discretisation, double precision structure-7 and 9 matrix.
  !
  ! For a rerefence about the stabilisation, see
  ! [Ouazzi, A.; Finite Element Simulation of Nonlinear Fluids, Application
  ! to Granular Material and Powder; Shaker Verlag, ISBN 3-8322-5201-0, p. 55ff]
  !
  ! WARNING: For edge oriented stabilisation, the underlying matrix rmatrix
  !   must have an extended stencil! The matrix structure must be set up with
  !   the BILF_MATC_EDGEBASED switch!!!
!</description>

!<input>
  ! Stabilisation parameter. Standard=0.01
  real(DP), intent(in) :: dgamma

  ! 2nd stabilisation parameter. Standard=dgamma=0.0
  real(DP), intent(in) :: dgammastar

  ! Exponent for edge length weight in the jump stabilisation. Standard = 2.
  ! A value of 2 corresponds to a weight h_E^2, but this can be changed here.
  real(dp), intent(in) :: deojEdgeExp

  ! Multiplication factor for the stabilisation matrix when adding
  ! it to the global matrix. Standard value = 1.0.
  real(DP), intent(in) :: dtheta

  ! 1D cubature formula to use for line integration
  ! Standard = CUB_G2_1D.
  integer(I32), intent(in) :: ccubType

  ! Viscosity parameter for the matrix if viscosity is constant.
  real(DP), intent(in) :: dnu

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), target, optional :: rdiscretisation

  ! OPTIONAL: List of edges where the operator should be computed.
  ! If not present, the operator will be computed on all edges.
  integer, dimension(:), intent(in), optional :: InodeList

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! The system matrix to be modified. Must be format 7 or 9.
  type(t_matrixScalar), dimension(:,:), intent(inout), target :: RmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: irow, jcol, idof
  integer :: IMT
  integer :: ivt1,ivt2,NVT,NMT, IMTidx
  integer :: IEL
  integer :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  real(DP) :: dedgelength,dedgeweight,dphidx,dphidy,dpsidx,dpsidy,dcoeff, dom

  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  real(DP), dimension(:), pointer :: p_Da11, p_Da12, p_Da21, p_Da22

  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsTempl
  integer, dimension(EL_MAXNBAS*2), target :: Idofs

  ! Arrays saving the local DOF numbers belonging to the global
  ! DOF numbers in Idofs
  integer, dimension(EL_MAXNBAS*2),target :: IlocalDofs

  ! Renumbering strategy for local DOF`s
  integer, dimension(EL_MAXNBAS), target :: IlocalDofRenum

  ! Number of local DOF`s on the patch
  integer :: ndof

  ! Number of local degees of freedom for trial and test functions
  integer :: indofPerElement

  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriangulation

  ! Some triangulation arrays we need frequently
  integer, dimension(:,:), pointer :: p_IneighboursAtElement
  integer, dimension(:,:), pointer :: p_IelementsAtEdge
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtEdge
  real(DP), dimension(:,:), pointer :: p_DvertexCoords

  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! Underlying discretisation structure
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation

  ! Arrays for saving Jacobian determinants and matrices
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:), allocatable :: Kentry
  real(DP), dimension(:,:,:), allocatable :: Dentry
  real(DP), dimension(2,2) :: Dval

  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoType

  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag

  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  type(t_evalElementSet) :: revalElementSet

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:,:), allocatable, target :: DcubPtsRefOnAllEdges

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:,:), pointer :: p_DcubPtsRef

  ! Cubature point weights
  real(DP), dimension(CUB_MAXCUBP) :: Domega

  ! Cubature point coordinates on the reference element
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D,Dxi2D

  ! number of cubature points on the reference element
  integer :: ncubp,icubp

  ! Derivative specifiers
  logical, dimension(EL_MAXNDER) :: Bder

  ! Whether the viscosity is constant or not
  logical :: bconstViscosity

    ! Currently we support only constant viscosity
    bconstViscosity = .true.

    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rmatrixScalar(1,1)%p_rspatialDiscrTest%p_rtriangulation
    p_rdiscretisation => rmatrixScalar(1,1)%p_rspatialDiscrTest

    ! If a discretisation structure is present, take that one.
    if (present(rdiscretisation)) &
      p_rdiscretisation => rdiscretisation

    ! Get Kvert, Kadj, Kmid, Kmel,...
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IneighboursAtElement,&
                                p_IneighboursAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IelementsAtEdge,&
                                p_IelementsAtEdge)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,&
                                p_IverticesAtEdge)
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                  p_DvertexCoords)

    NVT = p_rtriangulation%NVT
    NMT = p_rtriangulation%NMT
    if (present(InodeList)) NMT = size(InodeList)

    ! Get KLD, KCol...
    call lsyssc_getbase_Kld (rmatrixScalar(1,1),p_KLD)
    call lsyssc_getbase_Kcol (rmatrixScalar(1,1),p_Kcol)
    call lsyssc_getbase_double (rmatrixScalar(1,1),p_Da11)
    call lsyssc_getbase_double (rmatrixScalar(1,2),p_Da12)
    call lsyssc_getbase_double (rmatrixScalar(2,1),p_Da21)
    call lsyssc_getbase_double (rmatrixScalar(2,2),p_Da22)

    ! Activate the one and only element distribution
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get the number of local DOF`s for trial and test functions
    indofPerElement = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_relementDistribution%celement)

    ! Assure thath the element spaces are compatible
    if (elem_igetShape(p_relementDistribution%celement) .ne. &
        elem_igetShape(rmatrixScalar(1,1)%p_rspatialDiscrTrial%RelementDistr(1)%celement)) then
      call output_line ('Element spaces incompatible!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_uniDP')
      call sys_halt()
    end if

    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference line [-1,1]
    call cub_getCubPoints(ccubType, ncubp, Dxi1D, Domega)

    ! Allocate arrays accepting cubature point coordinates.
    ! We have ncubp cubature points on every egde.
    ! DcubPtsRef saves all possible combinations, which edge of
    ! one element might interact with one edge of one neighbour
    ! element.
    allocate(DcubPtsRefOnAllEdges(NDIM2D,ncubp,NVE))

    ! Put the 1D cubature points from to all of the edges on the
    ! reference element, so we can access them quickly.
    do iedge = 1,NVE
      call trafo_mapCubPts1Dto2DRefQuad(iedge, ncubp, Dxi1D, Dxi2D)
      do i=1,ncubp
        DcubPtsRefOnAllEdges(1,i,iedge) = Dxi2D(i,1)
        DcubPtsRefOnAllEdges(2,i,iedge) = Dxi2D(i,2)
      end do
    end do

    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

    ! Allocate some memory to hold the cubature points on the reference element
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp,2))

    ! Allocate arrays saving the local matrices for all elements
    ! in an element set. We allocate the arrays large enough...
    allocate(Kentry(indofPerElement*2*indofPerElement*2))
    allocate(Dentry(indofPerElement*2*indofPerElement*2,2,2))

    ! Allocate memory for obtaining DOF`s:
    allocate(IdofsTempl(indofPerElement,2))

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock+1))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    !
    ! We allocate space for 3 instead of 2 elements. The reason is that
    ! we later permute the values of the basis functions to get
    ! a local numbering on the patch. That is also the reason, we allocate
    ! not indofPerElement elements, but even indofPerElement,
    ! which is more than enough space to hold the values of the DOF`s of
    ! a whole element patch.

    allocate(Dbas(indofPerElement*2, &
            elem_getMaxDerivative(p_relementDistribution%celement),&
            ncubp,3))

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag.
    cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

    ! Do not calculate coordinates on the reference element -- we do this manually.
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

    ! Set up which derivatives to compute in the basis functions: X/Y-derivative
    Bder = .false.
    Bder(DER_DERIV_X) = .true.
    Bder(DER_DERIV_Y) = .true.

    ! Fill the basis function arrays with 0. Essential, as only parts
    ! of them are overwritten later.
    Dbas = 0.0_DP

    ! We loop through all edges
    do IMTidx = 1,NMT

      ! Either process all edges or the specified ones.
      if (present(InodeList)) then
        IMT = InodeList (IMTidx)
      else
        IMT = IMTidx
      end if

!      ! Check if we have 1 or 2 elements on the current edge
!      if (p_IelementsAtEdge (2,IMT) .eq. 0) then
!
!        ! This only happens if we have a boundary edge.
!        ! The boundary handling is... doing nothing! So we can skip
!        ! the handling of that edge completely!
!        cycle
!
!      end if

      ! Check how many elements we have on that edge.
      IELcount = 1
      if (p_IelementsAtEdge (2,IMT) .ne. 0) then
        IELcount = 2
      end if

      ! On an example, we now show the relationship between all the DOF`s
      ! we have to consider in the current situation. We have two elements,
      ! let us say IEL1 and IEL2, with their local and global DOF`s
      ! (example for Q1~):
      !
      !   local DOF`s on the           corresponding global
      !   element and on the patch     DOF`s
      !
      !    +----4----+----7----+       +----20---+----50---+
      !    |    4    |    4    |       |         |         |
      !    |         |         |       |         |         |
      !    1 1     3 3 1     3 6       10        60        70
      !    |         |         |       |         |         |
      !    |    2    |    2    |       |         |         |
      !    +----2----+----3----+       +----40---+----30---+
      !        IEL1      IEL2              IEL1      IEL2
      !
      !
      ! On every element, we have 4 local DOF`s (1..4). On the other hand,
      ! we have "local DOF`s on the patch" (1..7), ehich we call "patch DOF`s"
      ! from now on. To every local DOF, there belongs a global DOF
      ! (10,20,30,... or whatever), which gives the coefficient of the basis
      ! function.
      !
      ! Numbering that way, the local DOF`s of IEL1 obviously coincide
      ! with the first couple of local DOF`s of the element patch! Only
      ! the local DOF`s of element IEL2 make trouble, as they have another
      ! numbering.

      ! Get the global DOF`s of the 1 or two elements
      call dof_locGlobMapping_mult(p_rdiscretisation, &
                                  p_IelementsAtEdge (1:IELcount,IMT), &
                                  IdofsTempl)

      ! Some of the DOF`s on element 2 may coincide with DOF`s on element 1.
      ! More precisely, some must coincide! Therefore, we now have to collect the
      ! DOF`s uniquely and to figure out, which local DOF`s of element 2
      ! must renumbered to the appropriate local patch-DOF (like in the
      ! example, where local DOF 1 of element IEL2 must be renumbered
      ! to local patch DOF 3!
      !
      ! As the first couple of local DOF`s of IEL1 coincide with the local
      ! DOF`s of the patch, we can simply copy them:

      ndof = indofPerElement
      Idofs(1:ndof) = IdofsTempl(1:ndof,1)

      ! Furthermore, we read IdofsTempl and store the DOF`s there in Idofs,
      ! skipping all DOF`s we already have and setting up the renumbering strategy
      ! of local DOF`s on element IEL2 to patch DOF`s.

      if (IELcount .gt. 1) then

        ! Collect all the DOF's.
        skiploop: do IDOFE = 1,indofPerElement

          ! Do we have the DOF?
          idof = IdofsTempl(IDOFE,IELcount)
          do JDOFE=1,ndof
            if (Idofs(JDOFE) .eq. idof) then
              ! Yes, we have it.
              ! That means, the local DOF idof has to be mapped to the
              ! local patch dof...
              IlocalDofRenum (IDOFE) = JDOFE

              ! Proceed with the next one.
              cycle skiploop
            end if
          end do

          ! We do not have that DOF! Append it to Idofs.
          ndof = ndof+1
          Idofs(ndof) = idof
          IlocalDofRenum (IDOFE) = ndof

          ! Save also the number of the local DOF.
          ! Note that the global DOF`s in IdofsTempl(1..indofPerElement)
          ! belong to the local DOF`s 1..indofPerElement -- in that order!
          IlocalDofs(ndof) = IDOFE

        end do skiploop

      else

        ! Copy the DOF's.
        ndof = indofPerElement
        do IDOFE = 1,ndof

          ! We do not have that DOF! Append it to Idofs.
          Idofs(IDOFE) = IdofsTempl(IDOFE,1)
          IlocalDofRenum(IDOFE) = 0
          IlocalDofs(IDOFE) = IDOFE

        end do

      end if

      ! Now we know: Our 'local' matrix (consisting of only these DOF`s we just
      ! calculated) is a ndofsTest*ndofsTrial matrix.
      !
      ! Now extract the corresponding entries from the matrix.
      ! Kentry is an index where the entries of the local matrix Dentry
      ! (which we build up later) can be found in the global matrix.
      !
      ! The test space gives us the rows; loop through them

      do IDOFE = 0,ndof-1

        ! The 'local' DOF IDOFE corresponds in the global matrix to line...
        irow = Idofs (1+IDOFE)

        ! Loop through that line to find the columns, indexed by the local
        ! DOF`s in the trial space.
        trialspaceloop: do JDOFE = 1,ndof

          do jcol = p_Kld(irow),p_Kld(irow+1)-1

            if (p_Kcol(jcol) .eq. Idofs(JDOFE)) then

              ! Found! Put as destination pointer to Kentry
              Kentry (IDOFE*ndof+JDOFE) = jcol

              ! Clear local matrix
              Dentry (IDOFE*ndof+JDOFE,:,:) = 0.0_DP

              ! Jump out of the loop, proceed with next column
              cycle trialspaceloop

            end if

          end do

          call output_line ('Matrix invalid! Trial-DOF not found!', &
              OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_uniDP')
          call sys_halt()

        end do trialspaceloop

      end do ! JDOFE

      ! Now we can set up the local matrix in Dentry. Later, we will plug it into
      ! the global matrix using the positions in Kentry.
      !
      ! The next step is to evaluate the basis functions in the cubature
      ! points on the edge. To compute the jump, this has to be done
      ! for both elements on the edge.
      !
      ! Figure out which edge on the current element is IMT.
      ! We need this local numbering later to determine on which edge we
      ! have to place the cubature points.
      do i = 1,IELcount

        IEL = p_IelementsAtEdge(i,IMT)
        do iedge = 1,NVE
          if (p_IedgesAtElement (iedge,IEL) .eq. IMT) exit
        end do

        ! Copy the coordinates of the corresponding cubature points
        ! to DcubPtsEval. We calculated the cubature points on the
        ! reference element in advance, so we can simply take them.

        p_DcubPtsRef(:,:,i) = DcubPtsRefOnAllEdges (:,:,iedge)

      end do

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementsAtEdge (1:IELcount,IMT), &
          ctrafoType,DpointsRef=p_DcubPtsRef, rperfconfig=rperfconfig)
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

      ! Apply the permutation of the local DOF`s on the test functions
      ! on element 2. The numbers of the local DOF`s on element 1
      ! coincides with the numbers of the local DOF`s on the patch.
      ! Those on element 2, we have to renumber according to the permutation
      ! so that they are in the correct order according to the DOF`s on the patch.
      !
      ! We copy the values of the basis functions to the space in
      ! p_DcubPtsTest which is reserved for the 3rd element!
      ! That way, Dbas has:
      !
      ! local DOF on element   :       1       2       3       4
      !
      ! Global DOF on element 1:      10      40      60      20
      ! Dbas(1..4,*,*,1):    d1psi10 d1psi40 d1psi60 d1psi20
      !
      ! Global DOF on elemenr 2:      60      30      70      50
      ! Dbas(1..4,*,*,2):    d2psi60 d2psi30 d2psi70 d2psi50
      !
      ! and will be transformed to:
      !
      ! local patch-DOF:            1       2       3       4       5      6        7
      ! global DOF:                10      40      60      20      30     70       50
      ! Dbas(1..7,*,*,1): d1psi10 d1psi40 d1psi60 d1psi20     0.0    0.0      0.0
      ! Dbas(1..7,*,*,2): --------------------------unused-----------------------
      ! Dbas(1..7,*,*,3):     0.0     0.0 d2psi60     0.0 d2psi30 d2psi70 d2psi50
      !
      ! ("d1psi10" = grad(psi_10) on element 1,
      !  "psi_10" = test function on global DOF #10)
      !
      ! That space Dbas(:,:,:,3) is unused up to now, initialise by 0.0 and
      ! partially overwrite with the reordered values of the basis functions.
      Dbas (:,:,:,3) = 0.0_DP
      if (IELcount .eq. 2) then
        ! DOF values are not transferred if there is no neighbour element.
        ! In this case, the DOF values stay =0.
        Dbas (IlocalDofRenum(1:indofPerElement),:,:,3) &
          = Dbas (1:indofPerElement,:,:,2)
      end if

      ! What do we have now? When looking at the example, we have:
      !
      ! ndof = 7
      ! Idofs(1..7)             = global DOF`s in local numbering 1..7
      ! Dbas(1..7,*,1..ncubp,1) = values of basis functions on element 1
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF`s only appearing at element 2
      ! Dbas(1..7,*,1..ncubp,3) = values of basis functions on element 2
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF`s only appearing at element 1
      !
      ! Now we can start to integrate using this.

      ! Get the length of the current edge. It serves as a "determinant"
      ! in the cubature, so we have to divide it by 2 as an edge on the unit interval
      ! [-1,1] has length 2.
      ivt1 = p_IverticesAtEdge (1,IMT)
      ivt2 = p_IverticesAtEdge (2,IMT)
      dedgelength = &
        sqrt ((p_DvertexCoords(1,ivt2)-p_DvertexCoords(1,ivt1))**2 &
            + (p_DvertexCoords(2,ivt2)-p_DvertexCoords(2,ivt1))**2 )
      dedgeweight = dedgelength * 0.5_DP

      ! Compute the coefficient in front of the integral:
      ! < Ju,v > = sum_E max(gammastar*nu*h_E, gamma*h_E^2) int_E [grad u] [grad v] ds
      dcoeff = max(dgammastar * dnu * dedgelength, &
                   dgamma * dedgelength**deojEdgeExp )

      ! Now we have the values of the basis functions in all the cubature
      ! points.
      !
      ! Integrate the jump over the edges. This calculates the local matrix.
      !
      ! Loop through the test basis functions
      do IDOFE = 1,ndof

        ! Loop through the trial basis functions
        do JDOFE = 1,ndof

          Dval = 0.0_DP

          ! Loop through the cubature points to calculate the integral
          ! of the jump. Note that for computing the jump, we have to
          ! look in the inverse order to the cubature points of the neighbour
          ! element!
          ! Remember that the values of the basis functions on the first element
          ! are in p_DbasXXXX (.,.,.,1), while those of the 2nd element are
          ! in p_DbasXXXX (.,.,.,3) (rather than in p_DbasXXXX (.,.,.,2))
          ! by the above construction!
          do icubp = 1,ncubp

            ! [ grad phi ]   ( jump in the derivative of trial basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dphidx = Dbas (JDOFE,DER_DERIV_X,icubp,1) &
                  - Dbas (JDOFE,DER_DERIV_X,ncubp-icubp+1,3)

            dphidy = Dbas (JDOFE,DER_DERIV_Y,icubp,1) &
                  - Dbas (JDOFE,DER_DERIV_Y,ncubp-icubp+1,3)

            ! [ grad psi ]   ( jump in the derivative of test basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dpsidx = Dbas (IDOFE,DER_DERIV_X,icubp,1) &
                  - Dbas (IDOFE,DER_DERIV_X,ncubp-icubp+1,3)

            dpsidy = Dbas (IDOFE,DER_DERIV_Y,icubp,1) &
                  - Dbas (IDOFE,DER_DERIV_Y,ncubp-icubp+1,3)


            ! Compute int_edge ( [grad phi_i] [grad phi_j] )
            dom = Domega(icubp) * dedgeweight
            !dval = dval + Domega(icubp) * dedgeweight * &
            !              (dphidx*dpsidx + dphidy*dpsidy)
            
            Dval(1,1) = Dval(1,1) + dom * dphidx * dpsidx
            Dval(1,2) = Dval(1,2) + dom * dphidy * dpsidx
            Dval(2,1) = Dval(2,1) + dom * dphidx * dpsidy
            Dval(2,2) = Dval(2,2) + dom * dphidy * dpsidy

          end do

          ! Add the contribution to the local matrix -- weighted by the
          ! Omega from the cubature formula and the length of the edge.
          icubp = (IDOFE-1)*ndof+JDOFE
          Dentry(icubp,1,1) = Dentry(icubp,1,1) + dcoeff*dval(1,1)
          Dentry(icubp,1,2) = Dentry(icubp,1,2) + dcoeff*dval(1,2)
          Dentry(icubp,2,1) = Dentry(icubp,2,1) + dcoeff*dval(2,1)
          Dentry(icubp,2,2) = Dentry(icubp,2,2) + dcoeff*dval(2,2)

        end do

      end do

      ! Incorporate our "local" system matrix
      ! into the global matrix. The position of each entry DENTRY(X,Y)
      ! in the global matrix array A was saved in element Kentry(X,Y)
      ! before.
      ! Kentry gives the position of the additive contributions in Dentry.
      ! The entry is weighted by the current dtheta, which is usually
      ! the weighting parameter of the corresponding THETA-scheme of a
      ! nonstationary simulation. For stationary simulations, dtheta is typically
      ! 1.0 which includes the local matrix into the global one directly.)

      do IDOFE = 0,ndof-1
        do JDOFE = 1,ndof
          icubp = Kentry(IDOFE*ndof+JDOFE)
          p_Da11(icubp) = p_Da11(icubp) + dtheta * Dentry (IDOFE*ndof+JDOFE,1,1)
          p_Da12(icubp) = p_Da12(icubp) + dtheta * Dentry (IDOFE*ndof+JDOFE,1,2)
          p_Da21(icubp) = p_Da21(icubp) + dtheta * Dentry (IDOFE*ndof+JDOFE,2,1)
          p_Da22(icubp) = p_Da22(icubp) + dtheta * Dentry (IDOFE*ndof+JDOFE,2,2)

        end do
      end do

      ! Proceed with next edge

    end do ! IMTidx

    ! Clean up allocated arrays and memory.
    deallocate(DcubPtsRefOnAllEdges)

    call elprep_releaseElementSet(revalElementSet)
    deallocate(p_DcubPtsRef)

    deallocate(Dbas)

    deallocate(Kentry)
    deallocate(Dentry)

    deallocate(IdofsTempl)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine stokesdbg_flowEoj ( &
      rmatrixScalar,dgamma,ccubType,&
      rdiscretisation,InodeList,rperfconfig)

!<description>
  ! Unified edge oriented jump stabilisation.
  !
  ! Adds the unified edge oriented jump stabilisation to the matrix rmatrix:
  ! <tex>
  ! $$< Ju,v > = \sum_E \max(\gamma^{*} * \nu * h_E, \gamma h_E^2)
  !            \int_E [grad u] [grad v] ds $$
  ! </tex>
  ! Uniform discretisation, double precision structure-7 and 9 matrix.
  !
  ! For a rerefence about the stabilisation, see
  ! [Ouazzi, A.; Finite Element Simulation of Nonlinear Fluids, Application
  ! to Granular Material and Powder; Shaker Verlag, ISBN 3-8322-5201-0, p. 55ff]
  !
  ! WARNING: For edge oriented stabilisation, the underlying matrix rmatrix
  !   must have an extended stencil! The matrix structure must be set up with
  !   the BILF_MATC_EDGEBASED switch!!!
!</description>

!<input>
  ! Stabilisation parameter. Standard=0.01
  real(DP), intent(in) :: dgamma

  ! 1D cubature formula to use for line integration
  ! Standard = CUB_G2_1D.
  integer(I32), intent(in) :: ccubType

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), target, optional :: rdiscretisation

  ! OPTIONAL: List of edges where the operator should be computed.
  ! If not present, the operator will be computed on all edges.
  integer, dimension(:), intent(in), optional :: InodeList

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! The system matrix to be modified. Must be format 7 or 9.
  type(t_matrixScalar), dimension(:,:), intent(inout), target :: RmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: irow, jcol, idof
  integer :: IMT
  integer :: ivt1,ivt2,NVT,NMT, IMTidx
  integer :: IEL
  integer :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  real(DP) :: dedgelength,dedgeweight,dphidx,dphidy,dpsidx,dpsidy,dcoeff, dom, detax, detay

  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  real(DP), dimension(:), pointer :: p_Da11, p_Da12, p_Da21, p_Da22

  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsTempl
  integer, dimension(EL_MAXNBAS*2), target :: Idofs

  ! Arrays saving the local DOF numbers belonging to the global
  ! DOF numbers in Idofs
  integer, dimension(EL_MAXNBAS*2),target :: IlocalDofs

  ! Renumbering strategy for local DOF`s
  integer, dimension(EL_MAXNBAS), target :: IlocalDofRenum

  ! Number of local DOF`s on the patch
  integer :: ndof

  ! Number of local degees of freedom for trial and test functions
  integer :: indofPerElement

  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriangulation

  ! Some triangulation arrays we need frequently
  integer, dimension(:,:), pointer :: p_IneighboursAtElement
  integer, dimension(:,:), pointer :: p_IelementsAtEdge
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtEdge
  real(DP), dimension(:,:), pointer :: p_DvertexCoords

  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! Underlying discretisation structure
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation

  ! Arrays for saving Jacobian determinants and matrices
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:), allocatable :: Kentry
  real(DP), dimension(:,:,:), allocatable :: Dentry
  real(DP), dimension(2,2) :: Dval

  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoType

  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag

  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  type(t_evalElementSet) :: revalElementSet

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:,:), allocatable, target :: DcubPtsRefOnAllEdges

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:,:), pointer :: p_DcubPtsRef

  ! Cubature point weights
  real(DP), dimension(CUB_MAXCUBP) :: Domega

  ! Cubature point coordinates on the reference element
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D,Dxi2D

  ! number of cubature points on the reference element
  integer :: ncubp,icubp

  ! Derivative specifiers
  logical, dimension(EL_MAXNDER) :: Bder

  ! Whether the viscosity is constant or not
  logical :: bconstViscosity

    ! Currently we support only constant viscosity
    bconstViscosity = .true.

    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rmatrixScalar(1,1)%p_rspatialDiscrTest%p_rtriangulation
    p_rdiscretisation => rmatrixScalar(1,1)%p_rspatialDiscrTest

    ! If a discretisation structure is present, take that one.
    if (present(rdiscretisation)) &
      p_rdiscretisation => rdiscretisation

    ! Get Kvert, Kadj, Kmid, Kmel,...
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IneighboursAtElement,&
                                p_IneighboursAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IelementsAtEdge,&
                                p_IelementsAtEdge)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,&
                                p_IverticesAtEdge)
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                  p_DvertexCoords)

    NVT = p_rtriangulation%NVT
    NMT = p_rtriangulation%NMT
    if (present(InodeList)) NMT = size(InodeList)

    ! Get KLD, KCol...
    call lsyssc_getbase_Kld (rmatrixScalar(1,1),p_KLD)
    call lsyssc_getbase_Kcol (rmatrixScalar(1,1),p_Kcol)
    call lsyssc_getbase_double (rmatrixScalar(1,1),p_Da11)
    call lsyssc_getbase_double (rmatrixScalar(1,2),p_Da12)
    call lsyssc_getbase_double (rmatrixScalar(2,1),p_Da21)
    call lsyssc_getbase_double (rmatrixScalar(2,2),p_Da22)

    ! Activate the one and only element distribution
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get the number of local DOF`s for trial and test functions
    indofPerElement = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_relementDistribution%celement)

    ! Assure thath the element spaces are compatible
    if (elem_igetShape(p_relementDistribution%celement) .ne. &
        elem_igetShape(rmatrixScalar(1,1)%p_rspatialDiscrTrial%RelementDistr(1)%celement)) then
      call output_line ('Element spaces incompatible!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_uniDP')
      call sys_halt()
    end if

    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference line [-1,1]
    call cub_getCubPoints(ccubType, ncubp, Dxi1D, Domega)

    ! Allocate arrays accepting cubature point coordinates.
    ! We have ncubp cubature points on every egde.
    ! DcubPtsRef saves all possible combinations, which edge of
    ! one element might interact with one edge of one neighbour
    ! element.
    allocate(DcubPtsRefOnAllEdges(NDIM2D,ncubp,NVE))

    ! Put the 1D cubature points from to all of the edges on the
    ! reference element, so we can access them quickly.
    do iedge = 1,NVE
      call trafo_mapCubPts1Dto2DRefQuad(iedge, ncubp, Dxi1D, Dxi2D)
      do i=1,ncubp
        DcubPtsRefOnAllEdges(1,i,iedge) = Dxi2D(i,1)
        DcubPtsRefOnAllEdges(2,i,iedge) = Dxi2D(i,2)
      end do
    end do

    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

    ! Allocate some memory to hold the cubature points on the reference element
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp,2))

    ! Allocate arrays saving the local matrices for all elements
    ! in an element set. We allocate the arrays large enough...
    allocate(Kentry(indofPerElement*2*indofPerElement*2))
    allocate(Dentry(indofPerElement*2*indofPerElement*2,2,2))

    ! Allocate memory for obtaining DOF`s:
    allocate(IdofsTempl(indofPerElement,2))

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock+1))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    !
    ! We allocate space for 3 instead of 2 elements. The reason is that
    ! we later permute the values of the basis functions to get
    ! a local numbering on the patch. That is also the reason, we allocate
    ! not indofPerElement elements, but even indofPerElement,
    ! which is more than enough space to hold the values of the DOF`s of
    ! a whole element patch.

    allocate(Dbas(indofPerElement*2, &
            elem_getMaxDerivative(p_relementDistribution%celement),&
            ncubp,3))

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag.
    cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

    ! Do not calculate coordinates on the reference element -- we do this manually.
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

    ! Set up which derivatives to compute in the basis functions: X/Y-derivative
    Bder = .false.
    Bder(DER_FUNC) = .true.
    !Bder(DER_DERIV_X) = .true.
    !Bder(DER_DERIV_Y) = .true.

    ! Fill the basis function arrays with 0. Essential, as only parts
    ! of them are overwritten later.
    Dbas = 0.0_DP

    ! We loop through all edges
    do IMTidx = 1,NMT

      ! Either process all edges or the specified ones.
      if (present(InodeList)) then
        IMT = InodeList (IMTidx)
      else
        IMT = IMTidx
      end if

!      ! Check if we have 1 or 2 elements on the current edge
!      if (p_IelementsAtEdge (2,IMT) .eq. 0) then
!
!        ! This only happens if we have a boundary edge.
!        ! The boundary handling is... doing nothing! So we can skip
!        ! the handling of that edge completely!
!        cycle
!
!      end if

      detax = p_DvertexCoords(2, p_IverticesAtEdge(1,IMT)) - p_DvertexCoords(2, p_IverticesAtEdge(2,IMT))
      detay = p_DvertexCoords(1, p_IverticesAtEdge(2,IMT)) - p_DvertexCoords(1, p_IverticesAtEdge(1,IMT))
      dom = 1.0_DP / sqrt(detax**2+detay**2)
      detax = detax*dom
      detay = detay*dom

      ! Check how many elements we have on that edge.
      IELcount = 1
      if (p_IelementsAtEdge (2,IMT) .ne. 0) then
        IELcount = 2
      end if

      ! On an example, we now show the relationship between all the DOF`s
      ! we have to consider in the current situation. We have two elements,
      ! let us say IEL1 and IEL2, with their local and global DOF`s
      ! (example for Q1~):
      !
      !   local DOF`s on the           corresponding global
      !   element and on the patch     DOF`s
      !
      !    +----4----+----7----+       +----20---+----50---+
      !    |    4    |    4    |       |         |         |
      !    |         |         |       |         |         |
      !    1 1     3 3 1     3 6       10        60        70
      !    |         |         |       |         |         |
      !    |    2    |    2    |       |         |         |
      !    +----2----+----3----+       +----40---+----30---+
      !        IEL1      IEL2              IEL1      IEL2
      !
      !
      ! On every element, we have 4 local DOF`s (1..4). On the other hand,
      ! we have "local DOF`s on the patch" (1..7), ehich we call "patch DOF`s"
      ! from now on. To every local DOF, there belongs a global DOF
      ! (10,20,30,... or whatever), which gives the coefficient of the basis
      ! function.
      !
      ! Numbering that way, the local DOF`s of IEL1 obviously coincide
      ! with the first couple of local DOF`s of the element patch! Only
      ! the local DOF`s of element IEL2 make trouble, as they have another
      ! numbering.

      ! Get the global DOF`s of the 1 or two elements
      call dof_locGlobMapping_mult(p_rdiscretisation, &
                                  p_IelementsAtEdge (1:IELcount,IMT), &
                                  IdofsTempl)

      ! Some of the DOF`s on element 2 may coincide with DOF`s on element 1.
      ! More precisely, some must coincide! Therefore, we now have to collect the
      ! DOF`s uniquely and to figure out, which local DOF`s of element 2
      ! must renumbered to the appropriate local patch-DOF (like in the
      ! example, where local DOF 1 of element IEL2 must be renumbered
      ! to local patch DOF 3!
      !
      ! As the first couple of local DOF`s of IEL1 coincide with the local
      ! DOF`s of the patch, we can simply copy them:

      ndof = indofPerElement
      Idofs(1:ndof) = IdofsTempl(1:ndof,1)

      ! Furthermore, we read IdofsTempl and store the DOF`s there in Idofs,
      ! skipping all DOF`s we already have and setting up the renumbering strategy
      ! of local DOF`s on element IEL2 to patch DOF`s.

      if (IELcount .gt. 1) then

        ! Collect all the DOF's.
        skiploop: do IDOFE = 1,indofPerElement

          ! Do we have the DOF?
          idof = IdofsTempl(IDOFE,IELcount)
          do JDOFE=1,ndof
            if (Idofs(JDOFE) .eq. idof) then
              ! Yes, we have it.
              ! That means, the local DOF idof has to be mapped to the
              ! local patch dof...
              IlocalDofRenum (IDOFE) = JDOFE

              ! Proceed with the next one.
              cycle skiploop
            end if
          end do

          ! We do not have that DOF! Append it to Idofs.
          ndof = ndof+1
          Idofs(ndof) = idof
          IlocalDofRenum (IDOFE) = ndof

          ! Save also the number of the local DOF.
          ! Note that the global DOF`s in IdofsTempl(1..indofPerElement)
          ! belong to the local DOF`s 1..indofPerElement -- in that order!
          IlocalDofs(ndof) = IDOFE

        end do skiploop

      else

        ! Copy the DOF's.
        ndof = indofPerElement
        do IDOFE = 1,ndof

          ! We do not have that DOF! Append it to Idofs.
          Idofs(IDOFE) = IdofsTempl(IDOFE,1)
          IlocalDofRenum(IDOFE) = 0
          IlocalDofs(IDOFE) = IDOFE

        end do

      end if

      ! Now we know: Our 'local' matrix (consisting of only these DOF`s we just
      ! calculated) is a ndofsTest*ndofsTrial matrix.
      !
      ! Now extract the corresponding entries from the matrix.
      ! Kentry is an index where the entries of the local matrix Dentry
      ! (which we build up later) can be found in the global matrix.
      !
      ! The test space gives us the rows; loop through them

      do IDOFE = 0,ndof-1

        ! The 'local' DOF IDOFE corresponds in the global matrix to line...
        irow = Idofs (1+IDOFE)

        ! Loop through that line to find the columns, indexed by the local
        ! DOF`s in the trial space.
        trialspaceloop: do JDOFE = 1,ndof

          do jcol = p_Kld(irow),p_Kld(irow+1)-1

            if (p_Kcol(jcol) .eq. Idofs(JDOFE)) then

              ! Found! Put as destination pointer to Kentry
              Kentry (IDOFE*ndof+JDOFE) = jcol

              ! Clear local matrix
              Dentry (IDOFE*ndof+JDOFE,:,:) = 0.0_DP

              ! Jump out of the loop, proceed with next column
              cycle trialspaceloop

            end if

          end do

          call output_line ('Matrix invalid! Trial-DOF not found!', &
              OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_uniDP')
          call sys_halt()

        end do trialspaceloop

      end do ! JDOFE

      ! Now we can set up the local matrix in Dentry. Later, we will plug it into
      ! the global matrix using the positions in Kentry.
      !
      ! The next step is to evaluate the basis functions in the cubature
      ! points on the edge. To compute the jump, this has to be done
      ! for both elements on the edge.
      !
      ! Figure out which edge on the current element is IMT.
      ! We need this local numbering later to determine on which edge we
      ! have to place the cubature points.
      do i = 1,IELcount

        IEL = p_IelementsAtEdge(i,IMT)
        do iedge = 1,NVE
          if (p_IedgesAtElement (iedge,IEL) .eq. IMT) exit
        end do

        ! Copy the coordinates of the corresponding cubature points
        ! to DcubPtsEval. We calculated the cubature points on the
        ! reference element in advance, so we can simply take them.

        p_DcubPtsRef(:,:,i) = DcubPtsRefOnAllEdges (:,:,iedge)

      end do

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementsAtEdge (1:IELcount,IMT), &
          ctrafoType,DpointsRef=p_DcubPtsRef, rperfconfig=rperfconfig)
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

      ! Apply the permutation of the local DOF`s on the test functions
      ! on element 2. The numbers of the local DOF`s on element 1
      ! coincides with the numbers of the local DOF`s on the patch.
      ! Those on element 2, we have to renumber according to the permutation
      ! so that they are in the correct order according to the DOF`s on the patch.
      !
      ! We copy the values of the basis functions to the space in
      ! p_DcubPtsTest which is reserved for the 3rd element!
      ! That way, Dbas has:
      !
      ! local DOF on element   :       1       2       3       4
      !
      ! Global DOF on element 1:      10      40      60      20
      ! Dbas(1..4,*,*,1):    d1psi10 d1psi40 d1psi60 d1psi20
      !
      ! Global DOF on elemenr 2:      60      30      70      50
      ! Dbas(1..4,*,*,2):    d2psi60 d2psi30 d2psi70 d2psi50
      !
      ! and will be transformed to:
      !
      ! local patch-DOF:            1       2       3       4       5      6        7
      ! global DOF:                10      40      60      20      30     70       50
      ! Dbas(1..7,*,*,1): d1psi10 d1psi40 d1psi60 d1psi20     0.0    0.0      0.0
      ! Dbas(1..7,*,*,2): --------------------------unused-----------------------
      ! Dbas(1..7,*,*,3):     0.0     0.0 d2psi60     0.0 d2psi30 d2psi70 d2psi50
      !
      ! ("d1psi10" = grad(psi_10) on element 1,
      !  "psi_10" = test function on global DOF #10)
      !
      ! That space Dbas(:,:,:,3) is unused up to now, initialise by 0.0 and
      ! partially overwrite with the reordered values of the basis functions.
      Dbas (:,:,:,3) = 0.0_DP
      if (IELcount .eq. 2) then
        ! DOF values are not transferred if there is no neighbour element.
        ! In this case, the DOF values stay =0.
        Dbas (IlocalDofRenum(1:indofPerElement),:,:,3) &
          = Dbas (1:indofPerElement,:,:,2)
      end if

      ! What do we have now? When looking at the example, we have:
      !
      ! ndof = 7
      ! Idofs(1..7)             = global DOF`s in local numbering 1..7
      ! Dbas(1..7,*,1..ncubp,1) = values of basis functions on element 1
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF`s only appearing at element 2
      ! Dbas(1..7,*,1..ncubp,3) = values of basis functions on element 2
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF`s only appearing at element 1
      !
      ! Now we can start to integrate using this.

      ! Get the length of the current edge. It serves as a "determinant"
      ! in the cubature, so we have to divide it by 2 as an edge on the unit interval
      ! [-1,1] has length 2.
      ivt1 = p_IverticesAtEdge (1,IMT)
      ivt2 = p_IverticesAtEdge (2,IMT)
      dedgelength = &
        sqrt ((p_DvertexCoords(1,ivt2)-p_DvertexCoords(1,ivt1))**2 &
            + (p_DvertexCoords(2,ivt2)-p_DvertexCoords(2,ivt1))**2 )
      dedgeweight = dedgelength * 0.5_DP

      ! Compute the coefficient in front of the integral:
      ! < Ju,v > = sum_E max(gammastar*nu*h_E, gamma*h_E^2) int_E [grad u] [grad v] ds
      !dcoeff = max(dgammastar * dnu * dedgelength, &
      !             dgamma * dedgelength**deojEdgeExp )

      ! Now we have the values of the basis functions in all the cubature
      ! points.
      !
      ! Integrate the jump over the edges. This calculates the local matrix.
      !
      ! Loop through the test basis functions
      do IDOFE = 1,ndof

        ! Loop through the trial basis functions
        do JDOFE = 1,ndof

          Dval = 0.0_DP

          ! Loop through the cubature points to calculate the integral
          ! of the jump. Note that for computing the jump, we have to
          ! look in the inverse order to the cubature points of the neighbour
          ! element!
          ! Remember that the values of the basis functions on the first element
          ! are in p_DbasXXXX (.,.,.,1), while those of the 2nd element are
          ! in p_DbasXXXX (.,.,.,3) (rather than in p_DbasXXXX (.,.,.,2))
          ! by the above construction!
          do icubp = 1,ncubp

            ! [ grad phi ]   ( jump in the derivative of trial basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dphidx = detax*(Dbas (JDOFE,DER_FUNC,icubp,1) &
                          - Dbas (JDOFE,DER_FUNC,ncubp-icubp+1,3))

            dphidy =detay*(Dbas (JDOFE,DER_FUNC,icubp,1) &
                  - Dbas (JDOFE,DER_FUNC,ncubp-icubp+1,3))

            ! [ grad psi ]   ( jump in the derivative of test basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dpsidx = detax*(Dbas (IDOFE,DER_FUNC,icubp,1) &
                  - Dbas (IDOFE,DER_FUNC,ncubp-icubp+1,3))

            dpsidy = detay*(Dbas (IDOFE,DER_FUNC,icubp,1) &
                  - Dbas (IDOFE,DER_FUNC,ncubp-icubp+1,3))


            ! Compute int_edge ( [grad phi_i] [grad phi_j] )
            dom = Domega(icubp) * dedgeweight
            !dval = dval + Domega(icubp) * dedgeweight * &
            !              (dphidx*dpsidx + dphidy*dpsidy)
            
            Dval(1,1) = Dval(1,1) + dom * dphidx * dpsidx
            Dval(1,2) = Dval(1,2) + dom * dphidy * dpsidx
            Dval(2,1) = Dval(2,1) + dom * dphidx * dpsidy
            Dval(2,2) = Dval(2,2) + dom * dphidy * dpsidy

          end do

          ! Add the contribution to the local matrix -- weighted by the
          ! Omega from the cubature formula and the length of the edge.
          icubp = (IDOFE-1)*ndof+JDOFE
          Dentry(icubp,1,1) = Dentry(icubp,1,1) + dgamma*dval(1,1)
          Dentry(icubp,1,2) = Dentry(icubp,1,2) + dgamma*dval(1,2)
          Dentry(icubp,2,1) = Dentry(icubp,2,1) + dgamma*dval(2,1)
          Dentry(icubp,2,2) = Dentry(icubp,2,2) + dgamma*dval(2,2)

        end do

      end do

      ! Incorporate our "local" system matrix
      ! into the global matrix. The position of each entry DENTRY(X,Y)
      ! in the global matrix array A was saved in element Kentry(X,Y)
      ! before.
      ! Kentry gives the position of the additive contributions in Dentry.
      ! The entry is weighted by the current dtheta, which is usually
      ! the weighting parameter of the corresponding THETA-scheme of a
      ! nonstationary simulation. For stationary simulations, dtheta is typically
      ! 1.0 which includes the local matrix into the global one directly.)

      do IDOFE = 0,ndof-1
        do JDOFE = 1,ndof
          icubp = Kentry(IDOFE*ndof+JDOFE)
          p_Da11(icubp) = p_Da11(icubp) +  Dentry (IDOFE*ndof+JDOFE,1,1)
          p_Da12(icubp) = p_Da12(icubp) +  Dentry (IDOFE*ndof+JDOFE,1,2)
          p_Da21(icubp) = p_Da21(icubp) +  Dentry (IDOFE*ndof+JDOFE,2,1)
          p_Da22(icubp) = p_Da22(icubp) +  Dentry (IDOFE*ndof+JDOFE,2,2)

        end do
      end do

      ! Proceed with next edge

    end do ! IMTidx

    ! Clean up allocated arrays and memory.
    deallocate(DcubPtsRefOnAllEdges)

    call elprep_releaseElementSet(revalElementSet)
    deallocate(p_DcubPtsRef)

    deallocate(Dbas)

    deallocate(Kentry)
    deallocate(Dentry)

    deallocate(IdofsTempl)

  end subroutine

end module