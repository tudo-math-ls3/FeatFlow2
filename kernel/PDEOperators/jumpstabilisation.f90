!##############################################################################
!# ****************************************************************************
!# <name> jumpstabilisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines that implement the jump stabilisation.
!#
!# The following routines can be found in this module:
!#
!# 1.) jstab_calcUEOJumpStabil
!#     -> Modifies a scalar matrix to add the unified edge-oriented
!#        jump stabilisation term.
!#
!# 2.) jstab_matvecUEOJumpStabilBlk2d
!#     -> Performs a matrix vector multiplication with the jump stabilisation
!#        matrix for a 2d block vector.
!# 
!# 3.) jstab_calcReacJumpStabil
!#     -> Modifies a scalar matrix to add the reactive
!#        jump stabilisation term.
!#
!# 4.) jstab_matvecReacJumpStabilBlk2d
!#     -> Performs a matrix vector multiplication with the reactive jump
!#        stabilisation matrix for a 2d block vector.
!# 
!# Some auxiliary routines:
!#
!# 1.) jstab_ueoJumpStabil2d_m_unidble
!#     -> The actual matrix computation routine for 2D domains.
!#
!# 2.) jstab_reacJumpStabil2d_m_unidbl
!#     -> The actual matrix computation routine for 2D domains.
!#
!# 3.) jstab_ueoJumpStabil3d_m_unidble
!#     -> The actual matrix computation routine for 3D domains.
!#
!# 4.) jstab_reacJumpStabil3d_m_unidble
!#     -> The actual matrix computation routine for 3D domains.
!#
!# 5.) jstab_ueoJumpStabil1d_m_unidble
!#     -> The actual matrix computation routine for 1D domains.
!#
!# </purpose>
!##############################################################################

module jumpstabilisation

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
  
  private

  public :: jstab_calcUEOJumpStabilisation
  public :: jstab_matvecUEOJumpStabilBlk2d
  public :: jstab_calcReacJumpStabilisation
  public :: jstab_matvecReacJumpStabilBlk2d

contains

  ! ***************************************************************************

!<subroutine>

  subroutine jstab_calcUEOJumpStabilisation (&
      rmatrix,dgamma,dgammastar,deojEdgeExp,dtheta,ccubType,dnu,rdiscretisation,&
      InodeList,rperfconfig)

!<description>
  ! Edge oriented stabilisation technique. This routine incorporates the
  ! UEO stabilisation into a matrix.
!</description>

!<input>
  ! Stabilisation parameter. Standard = 0.01
  real(DP), intent(in) :: dgamma
  
  ! 2nd stabilisation parameter. Standard = 0
  real(DP), intent(in) :: dgammastar
  
  ! Exponent for edge length weight in the jump stabilisation. Standard = 2.
  ! A value of 2 corresponds to a weight h_E^2, but this can be changed here.
  real(dp), intent(in) :: deojEdgeExp

  ! Multiplication factor for the stabilisation matrix when adding
  ! it to the global matrix. Standard value = 1.0.
  real(DP), intent(in) :: dtheta
  
  ! 1D cubature formula to use for line integration.
  ! Standard = CUB_G2_1D.
  integer(I32), intent(in) :: ccubType
  
  ! Viscosity parameter for the matrix if viscosity is constant.
  real(DP), intent(in) :: dnu

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation
  
  ! OPTIONAL: List of edges/faces where the operator should be computed.
  ! If not present, the operator will be computed on all edges/faces.
  integer, dimension(:), intent(in), optional :: InodeList

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! Scalar matrix to be modified. The stabilisation term is added to rmatrix.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format 
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ('Unsupported matrix format.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_calcUEOJumpStabilisation')
      call sys_halt()
    end if

    if (rmatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM) then
      call output_line ('Unsupported discretisation.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_calcUEOJumpStabilisation')
      call sys_halt()
    end if
    
    ! Check which element we have here...
    select case(elem_igetShape(rmatrix%p_rspatialDiscrTest%RelementDistr(1)%celement))
    case (BGEOM_SHAPE_LINE)
      ! 1D line element
      call jstab_ueoJumpStabil1d_m_unidble (&
          rmatrix,dgamma,dgammastar,dtheta,dnu,rdiscretisation,rperfconfig)

    case (BGEOM_SHAPE_QUAD)
      ! 2D quadrilateral element
      call jstab_ueoJumpStabil2d_m_unidble (&
          rmatrix,dgamma,dgammastar,deojEdgeExp,dtheta,ccubType,dnu,rdiscretisation,&
          InodeList,rperfconfig)

    case (BGEOM_SHAPE_HEXA)
      ! 3D hexahedron element
      call jstab_ueoJumpStabil3d_m_unidble (&
          rmatrix,dgamma,dgammastar,dtheta,ccubType,dnu,rdiscretisation,rperfconfig)
    
    case default
      call output_line ('Unsupported element.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_calcUEOJumpStabilisation')
      call sys_halt()
      
    end select

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine jstab_ueoJumpStabil2d_m_unidble ( &
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
  type(t_matrixScalar), intent(inout), target :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: irow, jcol, idof
  integer :: IMT
  integer :: ivt1,ivt2,NVT,NMT, IMTidx
  integer :: IEL
  integer :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  real(DP) :: dval,dedgelength,dedgeweight,dphidx,dphidy,dpsidx,dpsidy,dcoeff
  
  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  real(DP), dimension(:), pointer :: p_Da
  
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
  real(DP), dimension(:), allocatable :: Dentry

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
  
  ! Whther the viscosity ís constant or not
  logical :: bconstViscosity
  
    ! Currently we support only constant viscosity
    bconstViscosity = .true.

    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rmatrixScalar%p_rspatialDiscrTest%p_rtriangulation
    p_rdiscretisation => rmatrixScalar%p_rspatialDiscrTest
    
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
    call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    call lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)
    call lsyssc_getbase_double (rmatrixScalar,p_Da)
    
    ! Activate the one and only element distribution
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get the number of local DOF`s for trial and test functions
    indofPerElement = elem_igetNDofLoc(p_relementDistribution%celement)
    
    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_relementDistribution%celement)
    
    ! Get the number of corner vertices of the element
    if (NVE .ne. elem_igetNVE(rmatrixScalar%p_rspatialDiscrTrial%&
        RelementDistr(1)%celement)) then
      call output_line ('Element spaces incompatible!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_unidble')
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
    allocate(Dentry(indofPerElement*2*indofPerElement*2))
    
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
              Dentry (IDOFE*ndof+JDOFE) = 0.0_DP
              
              ! Jump out of the loop, proceed with next column
              cycle trialspaceloop
            
            end if
            
          end do

          call output_line ('Matrix invalid! Trial-DOF not found!', &
              OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_unidble')
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
    
          dval = 0.0_DP
          
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
            dval = dval + Domega(icubp) * dedgeweight * &
                          (dphidx*dpsidx + dphidy*dpsidy)
          
          end do

          ! Add the contribution to the local matrix -- weighted by the
          ! Omega from the cubature formula and the length of the edge.
          Dentry ((IDOFE-1)*ndof+JDOFE) = &
            Dentry ((IDOFE-1)*ndof+JDOFE) + dcoeff*dval

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
        
          p_Da(Kentry(IDOFE*ndof+JDOFE)) = &
            p_Da(Kentry(IDOFE*ndof+JDOFE)) + &
            dtheta * Dentry (IDOFE*ndof+JDOFE)
        
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

    pure subroutine jstab3d_aux_calcQuadDetj(Dvtx, Dpts, Ddetj)
  
  !<description>
    ! INTERNAL AUXILIARY ROUTINE:
    ! This routine calculates the 'jacobian determinants' of a bilinear mapping
    ! from the 2D reference quadrilateral onto a face in 3D.
    ! The 'determinant' of a 3x2 jacobian matrix is defined as the euclid norm
    ! of the 3D cross-product of the jacobian matrix` columns.
  !</description>
  
  !<input>
    ! The coordinates of the face`s corner vertices.
    real(DP), dimension(3,4), intent(in) :: Dvtx
    
    ! The points for which the 'jacobian determinants' are to be calculated,
    ! given in 2D reference quadrilateral coordinates.
    real(DP), dimension(:,:), intent(in) :: Dpts
  !</input>
  
  !<output>
    ! The 'jacobian determinants' of the mapping.
    real(DP), dimension(:), intent(out) :: Ddetj
  !</output>
  
  !</subroutine>
    
    ! Some local variables
    integer :: ipt
    real(DP), dimension(3,3) :: Dtrafo
    real(DP), dimension(3,2) :: Djac
    real(DP), dimension(3) :: Dcp
    
      ! Calculate the coefficients of the bilinear mapping, but without
      ! the constant terms, as we will not need them here...
      Dtrafo(:,1) = 0.25_DP*(-Dvtx(:,1)+Dvtx(:,2)+Dvtx(:,3)-Dvtx(:,4))
      Dtrafo(:,2) = 0.25_DP*(-Dvtx(:,1)-Dvtx(:,2)+Dvtx(:,3)+Dvtx(:,4))
      Dtrafo(:,3) = 0.25_DP*( Dvtx(:,1)-Dvtx(:,2)+Dvtx(:,3)-Dvtx(:,4))
      
      ! Loop over all points
      do ipt = 1, ubound(Dpts,2)
      
        ! Calculate jacobian matrix of the bilinear mapping R^2 -> R^3
        Djac(:,1) = Dtrafo(:,1) + Dpts(2,ipt)*Dtrafo(:,3)
        Djac(:,2) = Dtrafo(:,2) + Dpts(1,ipt)*Dtrafo(:,3)
        
        ! Calculate cross-product of the columns of the jacobian matrix
        Dcp(1) = Djac(2,1)*Djac(3,2) - Djac(3,1)*Djac(2,2)
        Dcp(2) = Djac(3,1)*Djac(1,2) - Djac(1,1)*Djac(3,2)
        Dcp(3) = Djac(1,1)*Djac(2,2) - Djac(2,1)*Djac(1,2)
        
        ! Now the 'jacobian determinant' is the euclid norm of the
        ! cross-product vector we have calculated.
        Ddetj(ipt) = sqrt(Dcp(1)**2 + Dcp(2)**2 + Dcp(3)**2)
      
      end do ! ipt
    
    end subroutine jstab3d_aux_calcQuadDetj
  
  ! ***************************************************************************
    
  !<subroutine>
    
    pure subroutine jstab3d_aux_mapQuadToHexa(iat, itwist, Dpts2D, Dpts3D)
  
  !<description>
    ! INTERNAL AUXILIARY ROUTINE:
    ! Maps a set of points given on the 2D reference quadrilateral onto
    ! one of the local faces of a 3D reference hexahedron.
  !</description>
  
  !<input>
    ! The index of the local face onto which the points are to mapped.
    ! Is silently assumed to be 1 <= iat <= 6.
    integer, intent(in) :: iat
    
    ! The twist index of the hexahedron.
    integer(I32), intent(in) :: itwist
    
    ! The points which are to be mapped, given in 2D quadrilateral
    ! reference coordinates.
    real(DP), dimension(:,:), intent(in) :: Dpts2D
  !</input>
  
  !<output>
    ! The mapped points, given in 3D hexahedron reference coordinates.
    real(DP), dimension(:,:), intent(out) :: Dpts3D
  !</output>
  
  !</subroutine>
  
    real(DP), dimension(2,ubound(Dpts2D,2)) :: Dpts
    integer :: ipt
  
      ! Transform the points using the twist for this face
      select case(iand(int(ishft(itwist,-(9+3*iat))),7))
      case(0)
        Dpts(:,:) =  Dpts2D(:,:)
      case(1)
        Dpts(1,:) =  Dpts2D(2,:)
        Dpts(2,:) = -Dpts2D(1,:)
      case(2)
        Dpts(:,:) = -Dpts2D(:,:)
      case(3)
        Dpts(1,:) = -Dpts2D(2,:)
        Dpts(2,:) =  Dpts2D(1,:)
      case(4)
        Dpts(1,:) =  Dpts2D(2,:)
        Dpts(2,:) =  Dpts2D(1,:)
      case(5)
        Dpts(1,:) = -Dpts2D(1,:)
        Dpts(2,:) =  Dpts2D(2,:)
      case(6)
        Dpts(1,:) = -Dpts2D(2,:)
        Dpts(2,:) = -Dpts2D(1,:)
      case(7)
        Dpts(1,:) =  Dpts2D(1,:)
        Dpts(2,:) = -Dpts2D(2,:)
      end select
      
      ! Okay, which face to we have here?
      select case(iat)
      case(1) ! bottom face
        do ipt = 1, ubound(Dpts2D,2)
          Dpts3D(1,ipt) =  Dpts(1,ipt)
          Dpts3D(2,ipt) =  Dpts(2,ipt)
          Dpts3D(3,ipt) = -1.0_DP
        end do ! ipt
      case(2) ! front face
        do ipt = 1, ubound(Dpts2D,2)
          Dpts3D(1,ipt) =  Dpts(2,ipt)
          Dpts3D(2,ipt) = -1.0_DP
          Dpts3D(3,ipt) =  Dpts(1,ipt)
        end do ! ipt
      case(3) ! right face
        do ipt = 1, ubound(Dpts2D,2)
          Dpts3D(1,ipt) =  1.0_DP
          Dpts3D(2,ipt) = -Dpts(2,ipt)
          Dpts3D(3,ipt) = -Dpts(1,ipt)
        end do ! ipt
      case(4) ! back face
        do ipt = 1, ubound(Dpts2D,2)
          Dpts3D(1,ipt) = -Dpts(1,ipt)
          Dpts3D(2,ipt) =  1.0_DP
          Dpts3D(3,ipt) = -Dpts(2,ipt)
        end do ! ipt
      case(5) ! left face
        do ipt = 1, ubound(Dpts2D,2)
          Dpts3D(1,ipt) = -1.0_DP
          Dpts3D(2,ipt) =  Dpts(1,ipt)
          Dpts3D(3,ipt) =  Dpts(1,ipt)
        end do ! ipt
      case(6) ! top face
        do ipt = 1, ubound(Dpts2D,2)
          Dpts3D(1,ipt) = -Dpts(2,ipt)
          Dpts3D(2,ipt) = -Dpts(1,ipt)
          Dpts3D(3,ipt) =  1.0_DP
        end do ! ipt
      end select
    
    end subroutine jstab3d_aux_mapQuadToHexa

  ! ***************************************************************************
  
!<subroutine>

  subroutine jstab_ueoJumpStabil3d_m_unidble ( &
      rmatrixScalar,dgamma,dgammastar,dtheta,ccubType,dnu,rdiscretisation,rperfconfig)
      
!<description>
  ! Unified edge oriented jump stabilisation, 3D version for uniform
  ! hexahedron elements.
  !
  ! Adds the unified edge oriented jump stabilisation to the matrix rmatrix:
  ! <tex>
  ! $$< Ju,v > = \sum_E \max(\gamma^{*} * \nu * h_E, \gamma h_E^2) 
  !            \int_E [grad u] [grad v] ds$$
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
  
  ! 2nd stabilisation parameter. Standard=dgamma=0.01
  real(DP), intent(in) :: dgammastar
  
  ! Multiplication factor for the stabilisation matrix when adding
  ! it to the global matrix. Standard value = 1.0.
  real(DP), intent(in) :: dtheta
  
  ! 2D cubature formula to use for quadrilateral integration
  ! Standard = CUB_G2_2D.
  integer(I32), intent(in) :: ccubType
  
  ! Viscosity parameter for the matrix if viscosity is constant.
  real(DP), intent(in) :: dnu

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), target, optional :: rdiscretisation

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
  
!<inputoutput>
  ! The system matrix to be modified. Must be format 7 or 9.
  type(t_matrixScalar), intent(inout), target :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtria
  
  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relemDist
  
  ! Underlying discretisation structure
  type(t_spatialDiscretisation), pointer :: p_rdiscr
  
  ! Some triangulation arrays we need frequently
  integer, dimension(:,:), pointer :: p_IelementsAtFace
  integer, dimension(:,:), pointer :: p_IfacesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtFace
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer(I32), dimension(:), pointer :: p_ItwistIndex
  
  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  real(DP), dimension(:), pointer :: p_Da
  
  ! Number of DOFs - per element and per patch
  integer :: ndofs, ndofsPatch, idof,idofp,idof1,idof2
  
  ! Variables for the cubature formula
  integer :: ncubp, icubp
  real(DP), dimension(:,:), allocatable :: DcubPts2D
  real(DP), dimension(:,:,:), allocatable :: DcubPts3D
  real(DP), dimension(:), allocatable :: Domega, Dweights
  
  ! Variables for the DOF-mapping
  integer, dimension(:,:), allocatable :: Idofs
  
  ! Face corner vertices
  real(DP), dimension(3,4) :: DfaceVerts

  ! An element evaluation structure
  type(t_evalElementSet) :: reval
  
  ! Derivative specifiers
  logical, dimension(EL_MAXNDER) :: Bder
  
  ! Arrays for the evaluation of the elements
  real(DP), dimension(:,:,:,:), allocatable :: Dbas
  
  ! An array for the evaluation of the Jump
  real(DP), dimension(:,:,:), allocatable :: Djump
  
  ! An array for the local DOFs of the patch
  integer, dimension(:), allocatable :: IdofsPatch
  
  ! The local matrix
  real(DP), dimension(:,:), allocatable :: Dmatrix

  ! some other local variables
  integer(I32) :: celement, cevalTag, ctrafo
  integer :: i, j, iel1, iel2, iat, iat1, iat2, idx, irow
  integer, dimension(2) :: Iel
  real(DP) :: dquadArea, dcoeff, dvalue


    ! Get a pointer to the triangulation and discretisation.
    p_rtria => rmatrixScalar%p_rspatialDiscrTest%p_rtriangulation
    p_rdiscr => rmatrixScalar%p_rspatialDiscrTest
    
    ! If a discretisation structure is present, take that one.
    if (present(rdiscretisation)) &
      p_rdiscr => rdiscretisation

    ! Get the arrays from the triangulation
    call storage_getbase_int2d (p_rtria%h_IfacesAtElement, p_IfacesAtElement)
    call storage_getbase_int2d (p_rtria%h_IelementsAtFace, p_IelementsAtFace)
    call storage_getbase_int2d (p_rtria%h_IverticesAtFace, p_IverticesAtFace)
    call storage_getbase_double2d(p_rtria%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int32(p_rtria%h_ItwistIndex, p_ItwistIndex)
    
    ! Get KLD, KCol and Da
    call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    call lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)
    call lsyssc_getbase_double (rmatrixScalar,p_Da)
    
    ! Get the number of points for cubature formula
    ncubp = cub_igetNumPts(ccubType)
    
    ! Allocate the arrays for the 2D cubature rule
    allocate(DcubPts2D(2,ncubp))
    allocate(Domega(ncubp))
    allocate(Dweights(ncubp))
    
    ! Get the cubature formula itself
    call cub_getCubature(ccubType, DcubPts2D, Domega)

    ! Allocate an array for the 3D rules - the content is initialised later
    allocate(DcubPts3D(3,ncubp,2))

    ! Activate the one and only element distribution
    p_relemDist => p_rdiscr%RelementDistr(1)

    ! Get the number of DOFs
    celement = p_relemDist%celement
    ndofs = elem_igetNDofLoc(celement)
      
    ! Get the trafo type
    ctrafo = elem_igetTrafoType(celement)
    
    ! Get the evaluation tag of the element
    cevalTag = elem_getEvaluationTag(celement)
    
    ! Do not calculate reference coordiantes - we will do this manually.
    cevalTag = iand(cevalTag, not(EL_EVLTAG_REFPOINTS))

    ! Set up the Bder array - we will need the first derivatives
    Bder = .false.
    Bder(DER_DERIV3D_X) = .true.
    Bder(DER_DERIV3D_Y) = .true.
    Bder(DER_DERIV3D_Z) = .true.
    
    ! Allocate an array for the DOF-mapping
    allocate(Idofs(ndofs,2))
    
    ! Allocate an array for the element evaluation
    allocate(Dbas(ndofs,elem_getMaxDerivative(celement),ncubp,2))
    
    ! Allocate an array for the DOFs on the patch
    allocate(IdofsPatch(2*ndofs))
    
    ! Allocate an array for the jumps
    allocate(Djump(3,ncubp,2*ndofs))
    
    ! Allocate the local matrix
    allocate(Dmatrix(2*ndofs,2*ndofs))
    
    ! Okay, loop over all faces of the mesh
    do iat = 1, p_rtria%NAT
    
      ! Get the indices of the elements adjacent to the current face
      iel1 = p_IelementsAtFace(1,iat)
      iel2 = p_IelementsAtFace(2,iat)
      
      ! If iel2 is 0, then the face is a boundary face and we can skip it
      if(iel2 .eq. 0) cycle
      
      ! Copy (iel1,iel2) into an array
      Iel(:) = (/ iel1, iel2 /)
      
      ! Find out which of the element`s local faces corresponds
      ! to the global face we are currently processing.
      do iat1 = 1, 6
        if(p_IfacesAtElement(iat1, iel1) .eq. iat) exit
      end do
      do iat2 = 1, 6
        if(p_IfacesAtElement(iat2, iel2) .eq. iat) exit
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 1: Evaluate elements
      ! In this step, we will map the cubature points onto both hexahedra,
      ! evaluate the element on both cells and perform the DOF-mapping.
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Now let us calculate the reference coordinates for the two hexahedra
      call jstab3d_aux_mapQuadToHexa(iat1, p_ItwistIndex(iel1), DcubPts2D, &
                                     DcubPts3D(:,:,1))
      call jstab3d_aux_mapQuadToHexa(iat2, p_ItwistIndex(iel2), DcubPts2D, &
                                     DcubPts3D(:,:,2))

      ! Prepare the element evaluation structure
      call elprep_prepareSetForEvaluation (reval, cevalTag, p_rtria, Iel, &
                                           ctrafo, DpointsRef=DcubPts3D,&
                                           rperfconfig=rperfconfig)
      
      ! Evaluate the element
      call elem_generic_sim2(celement, reval, Bder, Dbas)

      ! Perform the DOF-mapping for both elements
      call dof_locGlobMapping_mult(p_rdiscr, Iel, Idofs)
    
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 2: Calculate Jumps
      ! In this step, we will calculate the gradient jumps and, at the same
      ! time, we will calculate the DOFs of the current element patch.
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Reset Djump
      Djump = 0.0_DP
      
      ! First, copy the evaluation of the first element`s gradients to the
      ! Djump array. At the same time, copy the DOF indices of the first
      ! element into IdofsPatch
      do idof = 1, ndofs
        do icubp = 1, ncubp
          Djump(1,icubp,idof) = Dbas(idof,DER_DERIV3D_X,icubp,1)
          Djump(2,icubp,idof) = Dbas(idof,DER_DERIV3D_Y,icubp,1)
          Djump(3,icubp,idof) = Dbas(idof,DER_DERIV3D_Z,icubp,1)
        end do ! icubp
        IdofsPatch(idof) = Idofs(idof,1)
      end do ! idof
      
      ! Currently, we have ndofs DOFs in our patch - these are the DOFs of
      ! the first element in the patch.
      ndofsPatch = ndofs
      
      ! Now comes the interesting part: Calculate the jumps.
      ! So loop over all local DOFs on the second element
      do idof = 1, ndofs

        ! Now let us see whether the current DOF (of the second element)
        ! is already in the patch. Please note that it is sufficient to check
        ! only the first ndofs entries in IdofsPatch, as all entries beyond
        ! ndofs belong to the second element that we are currently processing.
        do idofp = 1, ndofs
          if(IdofsPatch(idofp) .eq. Idofs(idof,2)) exit
        end do
        
        ! Update the DOF-map for the patch in the case the the current DOF was
        ! not already in the list.
        if(idofp .gt. ndofs) then
          ndofsPatch = ndofsPatch + 1
          idofp = ndofsPatch
          IdofsPatch(idofp) = Idofs(idof,2)
        end if
        
        ! Calculate the jump of the gradient in all points
        do icubp = 1, ncubp
          Djump(1,icubp,idofp) = Djump(1,icubp,idofp) &
                               - Dbas(idof,DER_DERIV3D_X,icubp,2)
          Djump(2,icubp,idofp) = Djump(2,icubp,idofp) &
                               - Dbas(idof,DER_DERIV3D_Y,icubp,2)
          Djump(3,icubp,idofp) = Djump(3,icubp,idofp) &
                               - Dbas(idof,DER_DERIV3D_Z,icubp,2)
        end do

      end do ! idof
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 3: Prepare for integration
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Now let us calculate the integration weights for the face.
      ! We will need the four corner vertices of the face for this.
      do i = 1, 4
        do j = 1, 3
          DfaceVerts(j,i) = p_DvertexCoords(j, p_IverticesAtFace(i,iat))
        end do
      end do
      
      ! Calculate the 'jacobian determinants' of the bilinear mapping
      call jstab3d_aux_calcQuadDetj(DfaceVerts, DcubPts2D, Dweights)
      
      ! And finally, multiply the 'jacobian determinants' by Domega to get
      ! the integration weights for the current face. At the same time,
      ! calculate the area of the face.
      dquadArea = 0.0_DP
      do icubp = 1, ncubp
        Dweights(icubp) = Dweights(icubp) * Domega(icubp)
        dquadArea = dquadArea + Dweights(icubp)
      end do
      
      ! Calculate the coefficient for the integral
      dcoeff = & !max(dgammastar * dnu * dquadArea, &
                !dgamma * dquadArea**2 !)
                dgamma * dquadArea

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 4: Calculate local matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Reset the local matrix for this patch
      Dmatrix = 0.0_DP
      
      ! Okay, let us loop over all DOFs in the current patch, once for the
      ! test and once for the trial space.
      do idof1 = 1, ndofsPatch
        do idof2 = 1, ndofsPatch
          
          ! Integrate the gradient jump
          dvalue = 0.0_DP
          do icubp = 1, ncubp
            dvalue = dvalue + Dweights(icubp) * ( &
                Djump(1,icubp,idof1)*Djump(1,icubp,idof2) &
              + Djump(2,icubp,idof1)*Djump(2,icubp,idof2) &
              + Djump(3,icubp,idof1)*Djump(3,icubp,idof2))
          end do ! icubp
          
          ! Store the entry in the local matrix, weighted by dcoeff
          Dmatrix(idof1,idof2) = dcoeff*dvalue

        end do ! idof2
      end do ! idof1

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 5: Update global matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Loop over all DOFs in test space
      do idof1 = 1, ndofsPatch
      
        ! Get the index of the row
        irow = IdofsPatch(idof1)
        
        ! Loop over all DOFs in trial space
        do idof2 = 1, ndofsPatch
          
          ! Try to find the entry in the global matrix
          do idx = p_Kld(irow), p_Kld(irow+1)-1
            if(p_Kcol(idx) .eq. IdofsPatch(idof2)) then
              p_Da(idx) = p_Da(idx) + dtheta*Dmatrix(idof1,idof2)
              exit
            end if 
          end do ! idx
          
        end do ! idof2
      
      end do ! idof1 
      
      ! Proceed with the next face
      
    end do ! iat
    
    ! Release the element set
    call elprep_releaseElementSet(reval)
    
    ! Deallocate everything we have allocated
    deallocate(Dmatrix)
    deallocate(Djump)
    deallocate(IdofsPatch)
    deallocate(Dbas)
    deallocate(Idofs)
    deallocate(DcubPts3D)
    deallocate(Dweights)
    deallocate(Domega)
    deallocate(DcubPts2D)
    
    ! That is it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine jstab_ueoJumpStabil1d_m_unidble ( &
      rmatrixScalar,dgamma,dgammastar,dtheta,dnu,rdiscretisation,rperfconfig)
      
!<description>
  ! Unified edge oriented jump stabilisation, 1D version for uniform
  ! hexahedron elements.
  !
  ! Adds the unified edge oriented jump stabilisation to the matrix rmatrix:
  ! <tex>
  ! $$< Ju,v > = \sum_E \max(\gamma^{*} * \nu * h_E, \gamma h_E^2) 
  !            \int_E [grad u] [grad v] ds$$
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
  
  ! 2nd stabilisation parameter. Standard=dgamma=0.01
  real(DP), intent(in) :: dgammastar
  
  ! Multiplication factor for the stabilisation matrix when adding
  ! it to the global matrix. Standard value = 1.0.
  real(DP), intent(in) :: dtheta
  
  ! Viscosity parameter for the matrix if viscosity is constant.
  real(DP), intent(in) :: dnu

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), target, optional :: rdiscretisation

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
  
!<inputoutput>
  ! The system matrix to be modified. Must be format 7 or 9.
  type(t_matrixScalar), intent(inout), target :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtria
  
  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relemDist
  
  ! Underlying discretisation structure
  type(t_spatialDiscretisation), pointer :: p_rdiscr
  
  ! Some triangulation arrays we need frequently
  integer, dimension(:), pointer :: p_IelemAtVert, p_IelemAtVertIdx
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  
  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  real(DP), dimension(:), pointer :: p_Da
  
  ! Number of DOFs - per element and per patch
  integer :: ndofs, ndofsPatch, idof,idofp,idof1,idof2
  
  ! Variables for the cubature formula
  real(DP), dimension(1,1,2) :: DcubPts1D
  
  ! Variables for the DOF-mapping
  integer, dimension(:,:), allocatable :: Idofs

  ! An element evaluation structure
  type(t_evalElementSet) :: reval
  
  ! Derivative specifiers
  logical, dimension(EL_MAXNDER) :: Bder
  
  ! Arrays for the evaluation of the elements
  real(DP), dimension(:,:,:,:), allocatable :: Dbas
  
  ! An array for the evaluation of the Jump
  real(DP), dimension(:), allocatable :: Djump
  
  ! An array for the local DOFs of the patch
  integer, dimension(:), allocatable :: IdofsPatch
  
  ! The local matrix
  real(DP), dimension(:,:), allocatable :: Dmatrix

  ! some other local variables
  integer(I32) :: celement, cevalTag, ctrafo
  integer :: iel1, iel2, ivt, ivt1, ivt2, idx, irow
  integer, dimension(2) :: Iel
  real(DP) :: dh, dcoeff


    ! Get a pointer to the triangulation and discretisation.
    p_rtria => rmatrixScalar%p_rspatialDiscrTest%p_rtriangulation
    p_rdiscr => rmatrixScalar%p_rspatialDiscrTest
    
    ! If a discretisation structure is present, take that one.
    if (present(rdiscretisation)) &
      p_rdiscr => rdiscretisation

    ! Get the arrays from the triangulation
    call storage_getbase_int2d (p_rtria%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int (p_rtria%h_IelementsAtVertex, p_IelemAtVert)
    call storage_getbase_int (p_rtria%h_IelementsAtVertexIdx, p_IelemAtVertIdx)
    call storage_getbase_double2d(p_rtria%h_DvertexCoords, p_DvertexCoords)
    
    ! Get KLD, KCol and Da
    call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    call lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)
    call lsyssc_getbase_double (rmatrixScalar,p_Da)
    
    ! Activate the one and only element distribution
    p_relemDist => p_rdiscr%RelementDistr(1)

    ! Get the number of DOFs
    celement = p_relemDist%celement
    ndofs = elem_igetNDofLoc(celement)
      
    ! Get the trafo type
    ctrafo = elem_igetTrafoType(celement)
    
    ! Get the evaluation tag of the element
    cevalTag = elem_getEvaluationTag(celement)
    
    ! Do not calculate reference coordiantes - we will do this manually.
    cevalTag = iand(cevalTag, not(EL_EVLTAG_REFPOINTS))

    ! Set up the Bder array - we will need the first derivatives
    Bder = .false.
    Bder(DER_DERIV1D_X) = .true.
    
    ! Allocate an array for the DOF-mapping
    allocate(Idofs(ndofs,2))
    
    ! Allocate an array for the element evaluation
    allocate(Dbas(ndofs,elem_getMaxDerivative(celement),2,2))
    
    ! Allocate an array for the DOFs on the patch
    allocate(IdofsPatch(2*ndofs))
    
    ! Allocate an array for the jumps
    allocate(Djump(2*ndofs))
    
    ! Allocate the local matrix
    allocate(Dmatrix(2*ndofs,2*ndofs))
    
    ! Okay, loop over all vertices of the mesh
    do ivt = 1, p_rtria%NVT
    
      ! Get the indices of the elements adjacent to the current vertice
      idx = p_IelemAtVertIdx(ivt)
      if(p_IelemAtVertIdx(ivt+1) - idx .le. 1) cycle
      iel1 = p_IelemAtVert(idx)
      iel2 = p_IelemAtVert(idx+1)
      
      ! Copy (iel1,iel2) into an array
      Iel(:) = (/ iel1, iel2 /)
      
      ! Find out which of the element`s local faces corresponds
      ! to the global face we are currently processing.
      do ivt1 = 1, 2
        if(p_IverticesAtElement(ivt1, iel1) .eq. ivt) exit
      end do
      do ivt2 = 1, 2
        if(p_IverticesAtElement(ivt2, iel2) .eq. ivt) exit
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 1: Evaluate elements
      ! Evaluate the element on both cells and perform the DOF-mapping.
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Set up the 'cubature points'
      DcubPts1D(1,1,1) = -1.0_DP + real(2*(ivt1-1),DP)
      DcubPts1D(1,1,2) = -1.0_DP + real(2*(ivt2-1),DP)

      ! Prepare the element evaluation structure
      call elprep_prepareSetForEvaluation (reval, cevalTag, p_rtria, Iel, &
                                           ctrafo, DpointsRef=DcubPts1D,&
                                           rperfconfig=rperfconfig)
      
      ! Evaluate the element
      call elem_generic_sim2(celement, reval, Bder, Dbas)

      ! Perform the DOF-mapping for both elements
      call dof_locGlobMapping_mult(p_rdiscr, Iel, Idofs)
    
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 2: Calculate Jumps
      ! In this step, we will calculate the gradient jumps and, at the same
      ! time, we will calculate the DOFs of the current element patch.
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Reset Djump
      Djump = 0.0_DP
      
      ! First, copy the evaluation of the first element`s gradients to the
      ! Djump array. At the same time, copy the DOF indices of the first
      ! element into IdofsPatch
      do idof = 1, ndofs
        Djump(idof) = Dbas(idof,DER_DERIV1D_X,1,1)
        IdofsPatch(idof) = Idofs(idof,1)
      end do ! idof
      
      ! Currently, we have ndofs DOFs in our patch - these are the DOFs of
      ! the first element in the patch.
      ndofsPatch = ndofs
      
      ! Now comes the interesting part: Calculate the jumps.
      ! So loop over all local DOFs on the second element
      do idof = 1, ndofs

        ! Now let us see whether the current DOF (of the second element)
        ! is already in the patch. Please note that it is sufficient to check
        ! only the first ndofs entries in IdofsPatch, as all entries beyond
        ! ndofs belong to the second element that we are currently processing.
        do idofp = 1, ndofs
          if(IdofsPatch(idofp) .eq. Idofs(idof,2)) exit
        end do
        
        ! Update the DOF-map for the patch in the case the the current DOF was
        ! not already in the list.
        if(idofp .gt. ndofs) then
          ndofsPatch = ndofsPatch + 1
          idofp = ndofsPatch
          IdofsPatch(idofp) = Idofs(idof,2)
        end if
        
        ! Calculate the jump of the gradient in all points
        Djump(idofp) = Djump(idofp) - Dbas(idof,DER_DERIV1D_X,1,2)

      end do ! idof
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 3: Prepare for 'integration'
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Calculate local h
      dh = 0.5_DP * abs(p_IverticesAtElement(3-ivt1,iel1) &
                       -p_IverticesAtElement(3-ivt2,iel2))
      
      ! Calculate the coefficient for the integral
      dcoeff = & !max(dgammastar * dnu * dh, &
                dgamma * dh**2

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 4: Calculate local matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Reset the local matrix for this patch
      Dmatrix = 0.0_DP
      
      ! Okay, let us loop over all DOFs in the current patch, once for the
      ! test and once for the trial space.
      do idof1 = 1, ndofsPatch
        do idof2 = 1, ndofsPatch
          
          ! Store the entry in the local matrix, weighted by dcoeff
          Dmatrix(idof1,idof2) = dcoeff*Djump(idof1)*Djump(idof2)

        end do ! idof2
      end do ! idof1

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 5: Update global matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Loop over all DOFs in test space
      do idof1 = 1, ndofsPatch
      
        ! Get the index of the row
        irow = IdofsPatch(idof1)
        
        ! Loop over all DOFs in trial space
        do idof2 = 1, ndofsPatch
          
          ! Try to find the entry in the global matrix
          do idx = p_Kld(irow), p_Kld(irow+1)-1
            if(p_Kcol(idx) .eq. IdofsPatch(idof2)) then
              p_Da(idx) = p_Da(idx) + dtheta*Dmatrix(idof1,idof2)
              exit
            end if 
          end do ! idx
          
        end do ! idof2
      
      end do ! idof1 
      
      ! Proceed with the next face
      
    end do ! iat
    
    ! Release the element set
    call elprep_releaseElementSet(reval)
    
    ! Deallocate everything we have allocated
    deallocate(Dmatrix)
    deallocate(Djump)
    deallocate(IdofsPatch)
    deallocate(Dbas)
    deallocate(Idofs)
    
    ! That is it
  
  end subroutine jstab_ueoJumpStabil1d_m_unidble
  
! ***************************************************************************

!<subroutine>

  subroutine jstab_matvecUEOJumpStabilBlk2d ( &
      dgamma,dgammastar,deojEdgeExp,ccubType,dnu,&
      rtemplateMat,rx,ry,cx,cy,rdiscretisation,&
      InodeList, rperfconfig)
!<description>
  ! Unified edge oriented jump stabilisation.
  !
  ! Performs a matrix vector multiplication
  !   y := cx * J x + cy y
  ! with J being the stabilisation matrix of the operator
  ! <tex>
  ! $$< Ju,v > = \sum_E \max(\gamma^{*} * \nu * h_E, \gamma h_E^2) 
  !            \int_E [grad u] [grad v] ds.$$
  ! </tex>
  !
  ! Uniform discretisation, double precision structure-7 and 9 matrix.
  ! The matrix rtemplateMat must contain the structure of a FE template
  ! matrix created with an extended matrix stencil: The matrix structure must 
  ! be set up with the BILF_MATC_EDGEBASED switch!!!
  !
  ! For a rerefence about the stabilisation, see
  ! [Ouazzi, A.; Finite Element Simulation of Nonlinear Fluids, Application
  ! to Granular Material and Powder; Shaker Verlag, ISBN 3-8322-5201-0, p. 55ff]
!</description>
  
!<input>

  ! Stabilisation parameter. Standard=0.01
  real(DP), intent(in) :: dgamma
  
  ! 2nd stabilisation parameter. Standard=dgamma=0.0
  real(DP), intent(in) :: dgammastar
  
  ! Exponent for edge length weight in the jump stabilisation. Standard = 2.
  ! A value of 2 corresponds to a weight h_E^2, but this can be changed here.
  real(dp), intent(in) :: deojEdgeExp

  ! 1D cubature formula to use for line integration
  ! Standard = CUB_G2_1D.
  integer(I32), intent(in) :: ccubType
  
  ! Viscosity parameter for the matrix if viscosity is constant.
  real(DP), intent(in) :: dnu

  ! Template FE matrix. Must be format 7 or 9 and created with an
  ! extended matrix stencil. This matrix is used for all dimensions.
  type(t_matrixScalar), intent(in), target :: rtemplateMat

  ! Solution vector x. Must be a 2D vector.
  type(t_vectorBlock), intent(in) :: rx

  ! Multiplication factor for x.
  real(DP), intent(in) :: cx

  ! Multiplication factor for y.
  real(DP), intent(in) :: cy

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
  ! RHS vector y where the result should be written to.
  ! Must be a 2D vector.
  type(t_vectorBlock), intent(inout) :: ry
!</inputoutput>

!</subroutine>


  ! local variables
  integer :: irow, jcol, idof
  integer :: IMT
  integer :: ivt1,ivt2,NVT,NMT,IMTidx
  integer :: IEL
  integer :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  real(DP) :: dval,dedgelength,dedgeweight,dphidx,dphidy,dpsidx,dpsidy,dcoeff
  
  ! Pointer to KLD, KCOL
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  
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
  
  ! Number of local degees of freedom for test functions
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
  real(DP), dimension(:), allocatable :: Dentry

  ! An element evaluation set for evaluating elements.
  type(t_evalElementSet) :: revalElementSet
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:,:), allocatable, target :: DcubPtsRefOnAllEdges

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:,:), allocatable :: p_DcubPtsRef

  ! Cubature point weights
  real(DP), dimension(CUB_MAXCUBP) :: Domega
  
  ! Cubature point coordinates on the reference element
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D,Dxi2D
  
  ! number of cubature points on the reference element
  integer :: ncubp,icubp
  
  ! Derivative specifiers
  logical, dimension(EL_MAXNDER) :: Bder

  ! Type of transformation from the reference to the real element 
  integer(I32) :: ctrafoType
  
  ! Data arrays for the vectors
  real(DP), dimension(:), pointer :: p_Dx1,p_Dy1
  real(DP), dimension(:), pointer :: p_Dx2,p_Dy2
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag
  
    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rtemplateMat%p_rspatialDiscrTest%p_rtriangulation
    p_rdiscretisation => rtemplateMat%p_rspatialDiscrTest
    
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
    
    ! If necessary, scale y in-advance, so we can just calculate y=cAx+y
    if (cy .eq. 0.0_DP) then
      call lsysbl_clearVector (ry)
    else if (cy .ne. 1.0_DP) then
      call lsysbl_scaleVector (ry,cy)
    end if
    
    ! Get KLD, KCol...
    call lsyssc_getbase_Kld (rtemplateMat,p_KLD)
    call lsyssc_getbase_Kcol (rtemplateMat,p_Kcol)
    
    ! Get the data vectors
    call lsyssc_getbase_double (rx%RvectorBlock(1),p_Dx1)
    call lsyssc_getbase_double (ry%RvectorBlock(1),p_Dy1)
    
    call lsyssc_getbase_double (rx%RvectorBlock(2),p_Dx2)
    call lsyssc_getbase_double (ry%RvectorBlock(2),p_Dy2)
    
    ! Activate the one and only element distribution
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get the number of local DOF`s for trial and test functions
    indofPerElement = elem_igetNDofLoc(p_relementDistribution%celement)
    
    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_relementDistribution%celement)
    
    ! Get the number of corner vertices of the element
    if (NVE .ne. elem_igetNVE(p_relementDistribution%celement)) then
      call output_line ('Element spaces incompatible!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'conv_ueoJumpStabil2d_double_uni')
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
    
    ! Get from the FE element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
    
    ! Allocate some memory to hold the cubature points on the reference element
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp,2))

    ! Allocate arrays saving the local matrices for all elements
    ! in an element set. We allocate the arrays large enough...
    allocate(Kentry(indofPerElement*2*indofPerElement*2))
    allocate(Dentry(indofPerElement*2*indofPerElement*2))
    
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
             
    ! Set up which derivatives to compute in the basis functions: X/Y-derivative
    Bder = .false.
    Bder(DER_DERIV_X) = .true.
    Bder(DER_DERIV_Y) = .true.

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

    ! Do not calculate coordinates on the reference element -- we do this manually.                    
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

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
      ! calculated) is a ndofs*ndofs matrix.
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
              Dentry (IDOFE*ndof+JDOFE) = 0.0_DP
              
              ! Jump out of the loop, proceed with next column
              cycle trialspaceloop
            
            end if
            
          end do

          call output_line ('Matrix invalid! Trial-DOF not found!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'conv_ueoJumpStabil2d_double_uni')
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
          ctrafoType, DpointsRef=p_DcubPtsRef, rperfconfig=rperfconfig)
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
    
          dval = 0.0_DP
          
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
            dval = dval + Domega(icubp) * dedgeweight * &
                          (dphidx*dpsidx + dphidy*dpsidy)
          
          end do

          ! Add the contribution to the local matrix -- weighted by the
          ! Omega from the cubature formula and the length of the edge.
          Dentry ((IDOFE-1)*ndof+JDOFE) = &
            Dentry ((IDOFE-1)*ndof+JDOFE) + dcoeff*dval

        end do
      
      end do
      
      ! Build the defect vector                     
      !     y = cAx + y
      ! This is done matrix free, only with the help of the local 
      ! matrix.                                                   
      ! In this case, D=(D1,D2) is expected to be the RHS on      
      ! entry and will be updated to be the defect vector when    
      ! this routine is left.                                     
      
      do IDOFE=0,ndof-1

        irow=Idofs(1+IDOFE)

        do JDOFE=1,ndof

          dval = cx * Dentry(IDOFE*ndof+JDOFE)         

          jcol=Idofs(JDOFE)
          p_Dy1(irow) = p_Dy1(irow) + dval*p_Dx1(jcol)
          p_Dy2(irow) = p_Dy2(irow) + dval*p_Dx2(jcol)

        end do
      end do

      ! Proceed with next edge
    
    end do ! IMT

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

  subroutine jstab_calcReacJumpStabilisation (&
      rmatrix,dgamma,dtheta,ccubType,dnu,rdiscretisation,rperfconfig)

!<description>
  ! Edge oriented stabilisation technique. This routine incorporates the
  ! reayctive jump stabilisation into a matrix.
!</description>

!<input>
  ! Stabilisation parameter. Standard=0.01
  real(DP), intent(in) :: dgamma
  
  ! Multiplication factor for the stabilisation matrix when adding
  ! it to the global matrix. Standard value = 1.0.
  real(DP), intent(in) :: dtheta
  
  ! 1D cubature formula to use for line integration.
  ! Standard = CUB_G2_1D.
  integer(I32), intent(in) :: ccubType
  
  ! Viscosity parameter for the matrix if viscosity is constant.
  real(DP), intent(in) :: dnu

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! Scalar matrix to be modified. The stabilisation term is added to rmatrix.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format 
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ('Unsupported matrix format.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_calcUEOJumpStabilisation')
      call sys_halt()
    end if

    if (rmatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM) then
      call output_line ('Unsupported discretisation.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_calcUEOJumpStabilisation')
      call sys_halt()
    end if

    ! Check which element we have here...
    select case(elem_igetShape(rmatrix%p_rspatialDiscrTest%RelementDistr(1)%celement))
    case (BGEOM_SHAPE_QUAD)
      ! 2D quadrilateral element
      call jstab_reacJumpStabil2d_m_unidbl (&
        rmatrix,dgamma,dtheta,ccubType,dnu,rdiscretisation,rperfconfig)

    case (BGEOM_SHAPE_HEXA)
      ! 3D hexahedron element
      call jstab_reacJumpStabil3d_m_unidbl (&
        rmatrix,dgamma,dtheta,ccubType,dnu,rdiscretisation,rperfconfig)
    
    case default
      call output_line ('Unsupported element.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_calcUEOJumpStabilisation')
      call sys_halt()
      
    end select

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine jstab_reacJumpStabil2d_m_unidbl ( &
      rmatrixScalar,dgamma,dtheta,ccubType,dnu,rdiscretisation,rperfconfig)
!<description>
  ! Unified edge oriented jump stabilisation.
  !
  ! Adds the reactive jump stabilisation to the matrix rmatrix:
  ! <tex>
  ! $$< Ju,v > = \sum_E \gamma \nu 1/|E| \int_E [u] [v] ds$$
  ! </tex>
  !
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

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
  
!<inputoutput>
  ! The system matrix to be modified. Must be format 7 or 9.
  type(t_matrixScalar), intent(inout), target :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: irow, jcol, idof
  integer :: IMT
  integer :: ivt1,ivt2,NVT
  integer :: IEL
  integer :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  real(DP) :: dval,dedgelength,dedgeweight,dphi,dpsi,dcoeff
  
  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  real(DP), dimension(:), pointer :: p_Da
  
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
  real(DP), dimension(:), allocatable :: Dentry

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
  
  ! Whther the viscosity ís constant or not
  logical :: bconstViscosity
  
    ! Currently we support only constant viscosity
    bconstViscosity = .true.

    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rmatrixScalar%p_rspatialDiscrTest%p_rtriangulation
    p_rdiscretisation => rmatrixScalar%p_rspatialDiscrTest
    
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
    
    ! Get KLD, KCol...
    call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    call lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)
    call lsyssc_getbase_double (rmatrixScalar,p_Da)
    
    ! Activate the one and only element distribution
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get the number of local DOF`s for trial and test functions
    indofPerElement = elem_igetNDofLoc(p_relementDistribution%celement)
    
    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_relementDistribution%celement)
    
    ! Get the number of corner vertices of the element
    if (NVE .ne. elem_igetNVE(rmatrixScalar%p_rspatialDiscrTrial%&
        RelementDistr(1)%celement)) then
      call output_line ('Element spaces incompatible!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_unidble')
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
    
    ! We have always 2 elements in each patch...
    IELcount = 2
      
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
    
    ! Allocate some memory to hold the cubature points on the reference element
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp,IELcount))

    ! Allocate arrays saving the local matrices for all elements
    ! in an element set. We allocate the arrays large enough...
    allocate(Kentry(indofPerElement*2*indofPerElement*2))
    allocate(Dentry(indofPerElement*2*indofPerElement*2))
    
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
    Bder(:) = .false.
    Bder(DER_FUNC) = .true.

      ! Fill the basis function arrays with 0. Essential, as only parts
      ! of them are overwritten later.
      Dbas = 0.0_DP
      
    ! We loop through all edges
    do IMT = 1,p_rtriangulation%NMT
    
      ! Check if we have 1 or 2 elements on the current edge
      if (p_IelementsAtEdge (2,IMT) .eq. 0) then
      
        ! This only happens if we have a boundary edge.
        ! The boundary handling is... doing nothing! So we can skip
        ! the handling of that edge completely!
        cycle
      
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
              Dentry (IDOFE*ndof+JDOFE) = 0.0_DP
              
              ! Jump out of the loop, proceed with next column
              cycle trialspaceloop
            
            end if
            
          end do

          call output_line ('Matrix invalid! Trial-DOF not found!', &
              OU_CLASS_ERROR,OU_MODE_STD,'jstab_ueoJumpStabil2d_m_unidble')
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
          ctrafoType, DpointsRef=p_DcubPtsRef, rperfconfig=rperfconfig)
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
      Dbas (IlocalDofRenum(1:indofPerElement),:,:,3) &
        = Dbas (1:indofPerElement,:,:,2)
      
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
      dedgeweight = dedgelength !* 0.5_DP
      
      ! Compute the coefficient in front of the integral:
      ! < Ju,v > = sum_E gamma*nu/|E| int_E [u] [v] ds
      dcoeff = dgamma * dnu / dedgelength
      
      ! Now we have the values of the basis functions in all the cubature 
      ! points.
      !
      ! Integrate the jump over the edges. This calculates the local matrix.
      !
      ! Loop through the test basis functions
      do IDOFE = 1,ndof
      
        ! Loop through the trial basis functions
        do JDOFE = 1,ndof
    
          dval = 0.0_DP
          
          ! Loop through the cubature points to calculate the integral
          ! of the jump. Note that for computing the jump, we have to
          ! look in the inverse order to the cubature points of the neighbour
          ! element!
          ! Remember that the values of the basis functions on the first element
          ! are in p_DbasXXXX (.,.,.,1), while those of the 2nd element are
          ! in p_DbasXXXX (.,.,.,3) (rather than in p_DbasXXXX (.,.,.,2))
          ! by the above construction!
          do icubp = 1,ncubp
          
            ! [ phi ]   ( jump in the derivative of trial basis function)
            ! = [ phi  ,  phi ]
            dphi  = Dbas (JDOFE,DER_FUNC,icubp,1) &
                  - Dbas (JDOFE,DER_FUNC,ncubp-icubp+1,3)

            ! [ grad psi ]   ( jump in the derivative of test basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dpsi  = Dbas (IDOFE,DER_FUNC,icubp,1) &
                  - Dbas (IDOFE,DER_FUNC,ncubp-icubp+1,3)
            
          
            ! Compute int_edge ( [grad phi_i] [grad phi_j] )
            dval = dval + Domega(icubp) * dedgeweight * (dphi*dpsi)
          
          end do

          ! Add the contribution to the local matrix -- weighted by the
          ! Omega from the cubature formula and the length of the edge.
          Dentry ((IDOFE-1)*ndof+JDOFE) = &
            Dentry ((IDOFE-1)*ndof+JDOFE) + dcoeff*dval

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
        
          p_Da(Kentry(IDOFE*ndof+JDOFE)) = &
            p_Da(Kentry(IDOFE*ndof+JDOFE)) + &
            dtheta * Dentry (IDOFE*ndof+JDOFE)
        
        end do
      end do
      
      ! Proceed with next edge
    
    end do ! IMT

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

  subroutine jstab_reacJumpStabil3d_m_unidbl ( &
      rmatrixScalar,dgamma,dtheta,ccubType,dnu,rdiscretisation,rperfconfig)
      
!<description>
  ! Unified edge oriented jump stabilisation, 3D version for uniform
  ! hexahedron elements.
  !
  ! Adds the reactive jump stabilisation to the matrix rmatrix:
  ! <tex>
  ! $$< Ju,v > = \sum_E \gamma \nu 1/|E| \int_E [u] [v] ds$$
  ! </tex>
  !
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
  
  ! Multiplication factor for the stabilisation matrix when adding
  ! it to the global matrix. Standard value = 1.0.
  real(DP), intent(in) :: dtheta
  
  ! 2D cubature formula to use for quadrilateral integration
  ! Standard = CUB_G2_2D.
  integer(I32), intent(in) :: ccubType
  
  ! Viscosity parameter for the matrix if viscosity is constant.
  real(DP), intent(in) :: dnu

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), target, optional :: rdiscretisation

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
  
!<inputoutput>
  ! The system matrix to be modified. Must be format 7 or 9.
  type(t_matrixScalar), intent(inout), target :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtria
  
  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relemDist
  
  ! Underlying discretisation structure
  type(t_spatialDiscretisation), pointer :: p_rdiscr
  
  ! Some triangulation arrays we need frequently
  integer, dimension(:,:), pointer :: p_IelementsAtFace
  integer, dimension(:,:), pointer :: p_IfacesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtFace
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer(I32), dimension(:), pointer :: p_ItwistIndex
  
  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  real(DP), dimension(:), pointer :: p_Da
  
  ! Number of DOFs - per element and per patch
  integer :: ndofs, ndofsPatch, idof,idofp,idof1,idof2
  
  ! Variables for the cubature formula
  integer :: ncubp, icubp
  real(DP), dimension(:,:), allocatable :: DcubPts2D
  real(DP), dimension(:,:,:), allocatable :: DcubPts3D
  real(DP), dimension(:), allocatable :: Domega, Dweights
  
  ! Variables for the DOF-mapping
  integer, dimension(:,:), allocatable :: Idofs
  
  ! Face corner vertices
  real(DP), dimension(3,4) :: DfaceVerts

  ! An element evaluation structure
  type(t_evalElementSet) :: reval
  
  ! Derivative specifiers
  logical, dimension(EL_MAXNDER) :: Bder
  
  ! Arrays for the evaluation of the elements
  real(DP), dimension(:,:,:,:), allocatable :: Dbas
  
  ! An array for the evaluation of the Jump
  real(DP), dimension(:,:), allocatable :: Djump
  
  ! An array for the local DOFs of the patch
  integer, dimension(:), allocatable :: IdofsPatch
  
  ! The local matrix
  real(DP), dimension(:,:), allocatable :: Dmatrix

  ! some other local variables
  integer(I32) :: celement, cevalTag, ctrafo
  integer :: i, j, iel1, iel2, iat, iat1, iat2, idx, irow
  integer, dimension(2) :: Iel
  real(DP) :: dquadArea, dcoeff, dvalue


    ! Get a pointer to the triangulation and discretisation.
    p_rtria => rmatrixScalar%p_rspatialDiscrTest%p_rtriangulation
    p_rdiscr => rmatrixScalar%p_rspatialDiscrTest
    
    ! If a discretisation structure is present, take that one.
    if (present(rdiscretisation)) &
      p_rdiscr => rdiscretisation

    ! Get the arrays from the triangulation
    call storage_getbase_int2d (p_rtria%h_IfacesAtElement, p_IfacesAtElement)
    call storage_getbase_int2d (p_rtria%h_IelementsAtFace, p_IelementsAtFace)
    call storage_getbase_int2d (p_rtria%h_IverticesAtFace, p_IverticesAtFace)
    call storage_getbase_double2d(p_rtria%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int32(p_rtria%h_ItwistIndex, p_ItwistIndex)
    
    ! Get KLD, KCol and Da
    call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    call lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)
    call lsyssc_getbase_double (rmatrixScalar,p_Da)
    
    ! Get the number of points for cubature formula
    ncubp = cub_igetNumPts(ccubType)
    
    ! Allocate the arrays for the 2D cubature rule
    allocate(DcubPts2D(2,ncubp))
    allocate(Domega(ncubp))
    allocate(Dweights(ncubp))
    
    ! Get the cubature formula itself
    call cub_getCubature(ccubType, DcubPts2D, Domega)

    ! Allocate an array for the 3D rules - the content is initialised later
    allocate(DcubPts3D(3,ncubp,2))

    ! Activate the one and only element distribution
    p_relemDist => p_rdiscr%RelementDistr(1)

    ! Get the number of DOFs
    celement = p_relemDist%celement
    ndofs = elem_igetNDofLoc(celement)
      
    ! Get the trafo type
    ctrafo = elem_igetTrafoType(celement)
    
    ! Get the evaluation tag of the element
    cevalTag = elem_getEvaluationTag(celement)
    
    ! Do not calculate reference coordiantes - we will do this manually.
    cevalTag = iand(cevalTag, not(EL_EVLTAG_REFPOINTS))

    ! Set up the Bder array - we will need the function values
    Bder = .false.
    Bder(DER_FUNC3D) = .true.
    
    ! Allocate an array for the DOF-mapping
    allocate(Idofs(ndofs,2))
    
    ! Allocate an array for the element evaluation
    allocate(Dbas(ndofs,elem_getMaxDerivative(celement),ncubp,2))
    
    ! Allocate an array for the DOFs on the patch
    allocate(IdofsPatch(2*ndofs))
    
    ! Allocate an array for the jumps
    allocate(Djump(ncubp,2*ndofs))
    
    ! Allocate the local matrix
    allocate(Dmatrix(2*ndofs,2*ndofs))
    
    ! Okay, loop over all faces of the mesh
    do iat = 1, p_rtria%NAT
    
      ! Get the indices of the elements adjacent to the current face
      iel1 = p_IelementsAtFace(1,iat)
      iel2 = p_IelementsAtFace(2,iat)
      
      ! If iel2 is 0, then the face is a boundary face and we can skip it
      if(iel2 .eq. 0) cycle
      
      ! Copy (iel1,iel2) into an array
      Iel(:) = (/ iel1, iel2 /)
      
      ! Find out which of the element`s local faces corresponds
      ! to the global face we are currently processing.
      do iat1 = 1, 6
        if(p_IfacesAtElement(iat1, iel1) .eq. iat) exit
      end do
      do iat2 = 1, 6
        if(p_IfacesAtElement(iat2, iel2) .eq. iat) exit
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 1: Evaluate elements
      ! In this step, we will map the cubature points onto both hexahedra,
      ! evaluate the element on both cells and perform the DOF-mapping.
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Now let us calculate the reference coordinates for the two hexahedra
      call jstab3d_aux_mapQuadToHexa(iat1, p_ItwistIndex(iel1), DcubPts2D, &
                                     DcubPts3D(:,:,1))
      call jstab3d_aux_mapQuadToHexa(iat2, p_ItwistIndex(iel2), DcubPts2D, &
                                     DcubPts3D(:,:,2))

      ! Prepare the element evaluation structure
      call elprep_prepareSetForEvaluation (reval, cevalTag, p_rtria, Iel, &
                                           ctrafo, DpointsRef=DcubPts3D,&
                                           rperfconfig=rperfconfig)
      
      ! Evaluate the element
      call elem_generic_sim2(celement, reval, Bder, Dbas)

      ! Perform the DOF-mapping for both elements
      call dof_locGlobMapping_mult(p_rdiscr, Iel, Idofs)
    
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 2: Calculate Jumps
      ! In this step, we will calculate the jumps and, at the same
      ! time, we will calculate the DOFs of the current element patch.
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Reset Djump
      Djump = 0.0_DP
      
      ! First, copy the evaluation of the first element`s values to the
      ! Djump array. At the same time, copy the DOF indices of the first
      ! element into IdofsPatch
      do idof = 1, ndofs
        do icubp = 1, ncubp
          Djump(icubp,idof) = Dbas(idof,DER_FUNC3D,icubp,1)
        end do ! icubp
        IdofsPatch(idof) = Idofs(idof,1)
      end do ! idof
      
      ! Currently, we have ndofs DOFs in our patch - these are the DOFs of
      ! the first element in the patch.
      ndofsPatch = ndofs
      
      ! Now comes the interesting part: Calculate the jumps.
      ! So loop over all local DOFs on the second element
      do idof = 1, ndofs

        ! Now let us see whether the current DOF (of the second element)
        ! is already in the patch. Please note that it is sufficient to check
        ! only the first ndofs entries in IdofsPatch, as all entries beyond
        ! ndofs belong to the second element that we are currently processing.
        do idofp = 1, ndofs
          if(IdofsPatch(idofp) .eq. Idofs(idof,2)) exit
        end do
        
        ! Update the DOF-map for the patch in the case the the current DOF was
        ! not already in the list.
        if(idofp .gt. ndofs) then
          ndofsPatch = ndofsPatch + 1
          idofp = ndofsPatch
          IdofsPatch(idofp) = Idofs(idof,2)
        end if
        
        ! Calculate the jump in all points
        do icubp = 1, ncubp
          Djump(icubp,idofp) = Djump(icubp,idofp) &
                             - Dbas(idof,DER_FUNC3D,icubp,2)
        end do

      end do ! idof
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 3: Prepare for integration
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Now let us calculate the integration weights for the face.
      ! We will need the four corner vertices of the face for this.
      do i = 1, 4
        do j = 1, 3
          DfaceVerts(j,i) = p_DvertexCoords(j, p_IverticesAtFace(i,iat))
        end do
      end do
      
      ! Calculate the 'jacobian determinants' of the bilinear mapping
      call jstab3d_aux_calcQuadDetj(DfaceVerts, DcubPts2D, Dweights)
      
      ! And finally, multiply the 'jacobian determinants' by Domega to get
      ! the integration weights for the current face. At the same time,
      ! calculate the area of the face.
      dquadArea = 0.0_DP
      do icubp = 1, ncubp
        Dweights(icubp) = Dweights(icubp) * Domega(icubp)
        dquadArea = dquadArea + Dweights(icubp)
      end do
      
      ! Calculate the coefficient for the integral
      dcoeff = dgamma * dnu / sqrt(dquadArea)

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 4: Calculate local matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Reset the local matrix for this patch
      Dmatrix = 0.0_DP
      
      ! Okay, let us loop over all DOFs in the current patch, once for the
      ! test and once for the trial space.
      do idof1 = 1, ndofsPatch
        do idof2 = 1, ndofsPatch
          
          ! Integrate the jump
          dvalue = 0.0_DP
          do icubp = 1, ncubp
            dvalue = dvalue + Dweights(icubp) *  &
                Djump(icubp,idof1)*Djump(icubp,idof2)
          end do ! icubp
          
          ! Store the entry in the local matrix, weighted by dcoeff
          Dmatrix(idof1,idof2) = dcoeff*dvalue

        end do ! idof2
      end do ! idof1

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! STEP 5: Update global matrix
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      ! Loop over all DOFs in test space
      do idof1 = 1, ndofsPatch
      
        ! Get the index of the row
        irow = IdofsPatch(idof1)
        
        ! Loop over all DOFs in trial space
        do idof2 = 1, ndofsPatch
          
          ! Try to find the entry in the global matrix
          do idx = p_Kld(irow), p_Kld(irow+1)-1
            if(p_Kcol(idx) .eq. IdofsPatch(idof2)) then
              p_Da(idx) = p_Da(idx) + dtheta*Dmatrix(idof1,idof2)
              exit
            end if 
          end do ! idx
          
        end do ! idof2
      
      end do ! idof1 
      
      ! Proceed with the next face
      
    end do ! iat
    
    ! Release the element set
    call elprep_releaseElementSet(reval)
    
    ! Deallocate everything we have allocated
    deallocate(Dmatrix)
    deallocate(Djump)
    deallocate(IdofsPatch)
    deallocate(Dbas)
    deallocate(Idofs)
    deallocate(DcubPts3D)
    deallocate(Dweights)
    deallocate(Domega)
    deallocate(DcubPts2D)
    
    ! That is it
  
  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine jstab_matvecReacJumpStabilBlk2d ( &
                  dgamma,ccubType,dnu,&
                  rtemplateMat,rx,ry,cx,cy,rdiscretisation,rperfconfig)
!<description>
  ! Unified edge oriented jump stabilisation.
  !
  ! Performs a matrix vector multiplication
  !   <tex> $$ y := cx * J x + cy y $$ </tex>
  ! with J being the stabilisation matrix of the operator
  ! <tex>
  ! $$< Ju,v > = \sum_E \gamma \nu 1/|E| \int_E [u] [v] ds.$$
  ! </tex>
  !
  ! Uniform discretisation, double precision structure-7 and 9 matrix.
  ! The matrix rtemplateMat must contain the structure of a FE template
  ! matrix created with an extended matrix stencil: The matrix structure must 
  ! be set up with the BILF_MATC_EDGEBASED switch!!!
  !
  ! For a rerefence about the stabilisation, see
  ! [Ouazzi, A.; Finite Element Simulation of Nonlinear Fluids, Application
  ! to Granular Material and Powder; Shaker Verlag, ISBN 3-8322-5201-0, p. 55ff]
!</description>
  
!<input>

  ! Stabilisation parameter. Standard=0.01
  real(DP), intent(in) :: dgamma
  
  ! 1D cubature formula to use for line integration
  ! Standard = CUB_G2_1D.
  integer(I32), intent(in) :: ccubType
  
  ! Viscosity parameter for the matrix if viscosity is constant.
  real(DP), intent(in) :: dnu

  ! Template FE matrix. Must be format 7 or 9 and created with an
  ! extended matrix stencil. This matrix is used for all dimensions.
  type(t_matrixScalar), intent(in), target :: rtemplateMat

  ! Solution vector x. Must be a 2D vector.
  type(t_vectorBlock), intent(in) :: rx

  ! Multiplication factor for x.
  real(DP), intent(in) :: cx

  ! Multiplication factor for y.
  real(DP), intent(in) :: cy

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), target, optional :: rdiscretisation

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! RHS vector y where the result should be written to.
  ! Must be a 2D vector.
  type(t_vectorBlock), intent(inout) :: ry
!</inputoutput>

!</subroutine>


  ! local variables
  integer :: irow, jcol, idof
  integer :: IMT
  integer :: ivt1,ivt2,NVT
  integer :: IEL
  integer :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  real(DP) :: dval,dedgelength,dedgeweight,dphi,dpsi,dcoeff
  
  ! Pointer to KLD, KCOL
  integer, dimension(:), pointer :: p_Kld
  integer, dimension(:), pointer :: p_Kcol
  
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
  
  ! Number of local degees of freedom for test functions
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
  real(DP), dimension(:), allocatable :: Dentry

  ! An element evaluation set for evaluating elements.
  type(t_evalElementSet) :: revalElementSet
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:,:), allocatable, target :: DcubPtsRefOnAllEdges

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:,:), allocatable :: p_DcubPtsRef

  ! Cubature point weights
  real(DP), dimension(CUB_MAXCUBP) :: Domega
  
  ! Cubature point coordinates on the reference element
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi1D,Dxi2D
  
  ! number of cubature points on the reference element
  integer :: ncubp,icubp
  
  ! Derivative specifiers
  logical, dimension(EL_MAXNDER) :: Bder

  ! Type of transformation from the reference to the real element 
  integer(I32) :: ctrafoType
  
  ! Data arrays for the vectors
  real(DP), dimension(:), pointer :: p_Dx1,p_Dy1
  real(DP), dimension(:), pointer :: p_Dx2,p_Dy2
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag
  
    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rtemplateMat%p_rspatialDiscrTest%p_rtriangulation
    p_rdiscretisation => rtemplateMat%p_rspatialDiscrTest
    
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
    
    ! If necessary, scale y in-advance, so we can just calculate y=cAx+y
    if (cy .eq. 0.0_DP) then
      call lsysbl_clearVector (ry)
    else if (cy .ne. 1.0_DP) then
      call lsysbl_scaleVector (ry,cy)
    end if
    
    ! Get KLD, KCol...
    call lsyssc_getbase_Kld (rtemplateMat,p_KLD)
    call lsyssc_getbase_Kcol (rtemplateMat,p_Kcol)
    
    ! Get the data vectors
    call lsyssc_getbase_double (rx%RvectorBlock(1),p_Dx1)
    call lsyssc_getbase_double (ry%RvectorBlock(1),p_Dy1)
    
    call lsyssc_getbase_double (rx%RvectorBlock(2),p_Dx2)
    call lsyssc_getbase_double (ry%RvectorBlock(2),p_Dy2)
    
    ! Activate the one and only element distribution
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get the number of local DOF`s for trial and test functions
    indofPerElement = elem_igetNDofLoc(p_relementDistribution%celement)
    
    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_relementDistribution%celement)
    
    ! Get the number of corner vertices of the element
    if (NVE .ne. elem_igetNVE(p_relementDistribution%celement)) then
      call output_line ('Element spaces incompatible!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'conv_ueoJumpStabil2d_double_uni')
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
    
    ! We have always 2 elements in each patch...
    IELcount = 2
      
    ! Get from the FE element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
    
    ! Allocate some memory to hold the cubature points on the reference element
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp,IELcount))

    ! Allocate arrays saving the local matrices for all elements
    ! in an element set. We allocate the arrays large enough...
    allocate(Kentry(indofPerElement*2*indofPerElement*2))
    allocate(Dentry(indofPerElement*2*indofPerElement*2))
    
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
             
    ! Set up which derivatives to compute in the basis functions: X/Y-derivative
    Bder = .false.
    Bder(DER_DERIV_X) = .true.
    Bder(DER_DERIV_Y) = .true.

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

    ! Do not calculate coordinates on the reference element -- we do this manually.                    
    cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

    ! Fill the basis function arrays with 0. Essential, as only parts
    ! of them are overwritten later.
    Dbas = 0.0_DP
    
    ! We loop through all edges
    do IMT = 1,p_rtriangulation%NMT
    
      ! Check if we have 1 or 2 elements on the current edge
      if (p_IelementsAtEdge (2,IMT) .eq. 0) then
      
        ! This only happens if we have a boundary edge.
        ! The boundary handling is... doing nothing! So we can skip
        ! the handling of that edge completely!
        cycle
      
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
        
      ! Now we know: Our 'local' matrix (consisting of only these DOF`s we just
      ! calculated) is a ndofs*ndofs matrix.
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
              Dentry (IDOFE*ndof+JDOFE) = 0.0_DP
              
              ! Jump out of the loop, proceed with next column
              cycle trialspaceloop
            
            end if
            
          end do

          call output_line ('Matrix invalid! Trial-DOF not found!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'conv_ueoJumpStabil2d_double_uni')
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
          ctrafoType, DpointsRef=p_DcubPtsRef, rperfconfig=rperfconfig)
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
      Dbas (IlocalDofRenum(1:indofPerElement),:,:,3) &
        = Dbas (1:indofPerElement,:,:,2)
      
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
      dedgeweight = dedgelength !* 0.5_DP *
      
      ! Compute the coefficient in front of the integral:
      ! < Ju,v > = sum_E max(gammastar*nu*h_E, gamma*h_E^2) int_E [grad u] [grad v] ds
      dcoeff = dgamma * dnu / dedgelength
      
      ! Now we have the values of the basis functions in all the cubature 
      ! points.
      !
      ! Integrate the jump over the edges. This calculates the local matrix.
      !
      ! Loop through the test basis functions
      do IDOFE = 1,ndof
      
        ! Loop through the trial basis functions
        do JDOFE = 1,ndof
    
          dval = 0.0_DP
          
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
            dphi  = Dbas (JDOFE,DER_FUNC,icubp,1) &
                  - Dbas (JDOFE,DER_FUNC,ncubp-icubp+1,3)

            ! [ grad psi ]   ( jump in the derivative of test basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dpsi  = Dbas (IDOFE,DER_FUNC,icubp,1) &
                  - Dbas (IDOFE,DER_FUNC,ncubp-icubp+1,3)

            ! Compute int_edge ( [grad phi_i] [grad phi_j] )
            dval = dval + Domega(icubp) * dedgeweight * (dphi*dpsi)
          
          end do

          ! Add the contribution to the local matrix -- weighted by the
          ! Omega from the cubature formula and the length of the edge.
          Dentry ((IDOFE-1)*ndof+JDOFE) = &
            Dentry ((IDOFE-1)*ndof+JDOFE) + dcoeff*dval

        end do
      
      end do
      
      ! Build the defect vector                     
      !     y = cAx + y
      ! This is done matrix free, only with the help of the local 
      ! matrix.                                                   
      ! In this case, D=(D1,D2) is expected to be the RHS on      
      ! entry and will be updated to be the defect vector when    
      ! this routine is left.                                     
      
      do IDOFE=0,ndof-1

        irow=Idofs(1+IDOFE)

        do JDOFE=1,ndof

          dval = cx * Dentry(IDOFE*ndof+JDOFE)         

          jcol=Idofs(JDOFE)
          p_Dy1(irow) = p_Dy1(irow) + dval*p_Dx1(jcol)
          p_Dy2(irow) = p_Dy2(irow) + dval*p_Dx2(jcol)

        end do
      end do

      ! Proceed with next edge
    
    end do ! IMT

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
