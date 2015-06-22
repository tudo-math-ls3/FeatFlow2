module stokesdbg_edge_div

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
use ucd
use spdiscprojection
use vectorio

implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine stokesdbg_writeEdgeDiv(rvecSol, sfilename,ccubType)
!<input>
  type(t_vectorScalar), dimension(:), intent(inout) :: rvecSol
  character(len=*), intent(in) :: sfilename
  integer(I32), intent(in) :: ccubType
!</input>

!</subroutine>

  ! local variables
  integer :: irow, jcol, idof
  integer :: IMT
  integer :: ivt1,ivt2,NVT,NMT, IMTidx
  integer :: IEL
  integer :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  real(DP) :: dedgelength,dedgeweight,dphidx,dphidy,dpsidx,dpsidy,dcoeff, dom

  ! Pointer to KLD, KCOL, DA
  real(DP), dimension(:), pointer :: p_Du1, p_Du2, p_Dvj, p_Dvm, p_Dwj, p_Dwm

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
  !integer, dimension(:), allocatable :: Kentry
  real(DP), dimension(:,:), allocatable :: Dentry
  real(DP), dimension(2) :: Dval

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

  type(t_spatialDiscretisation) :: rdiscQ1T
  type(t_vectorScalar) :: rvecJump, rvecMean
  type(t_ucdExport) :: rexport
  real(DP) :: dt

    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rvecSol(1)%p_rspatialDiscr%p_rtriangulation
    p_rdiscretisation => rvecSol(1)%p_rspatialDiscr

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

    ! allocate some vectors
    call spdiscr_initDiscr_simple(rdiscQ1T, EL_E030_2D, CUB_G1_2D, p_rtriangulation)
    call lsyssc_createVecByDiscr(rdiscQ1T, rvecJump, .true.)
    call lsyssc_createVecByDiscr(rdiscQ1T, rvecMean, .true.)
    call lsyssc_getbase_double(rvecJump, p_Dvj)
    call lsyssc_getbase_double(rvecMean, p_Dvm)
    
    ! Get KLD, KCol...
    call lsyssc_getbase_double (rvecSol(1),p_Du1)
    call lsyssc_getbase_double (rvecSol(2),p_Du2)

    ! Activate the one and only element distribution
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get the number of local DOF`s for trial and test functions
    indofPerElement = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_relementDistribution%celement)

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
        IMT = IMTidx

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
      else
        cycle
      end if


      ! Get the global DOF`s of the 1 or two elements
      call dof_locGlobMapping_mult(p_rdiscretisation, &
                                  p_IelementsAtEdge (1:IELcount,IMT), &
                                  IdofsTempl)
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
          ctrafoType,DpointsRef=p_DcubPtsRef)
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

      ! That space Dbas(:,:,:,3) is unused up to now, initialise by 0.0 and
      ! partially overwrite with the reordered values of the basis functions.
      Dbas (:,:,:,3) = 0.0_DP
      if (IELcount .eq. 2) then
        ! DOF values are not transferred if there is no neighbour element.
        ! In this case, the DOF values stay =0.
        Dbas (IlocalDofRenum(1:indofPerElement),:,:,3) &
          = Dbas (1:indofPerElement,:,:,2)
      end if

      ! Get the length of the current edge. It serves as a "determinant"
      ! in the cubature, so we have to divide it by 2 as an edge on the unit interval
      ! [-1,1] has length 2.
      ivt1 = p_IverticesAtEdge (1,IMT)
      ivt2 = p_IverticesAtEdge (2,IMT)
      dedgelength = &
        sqrt ((p_DvertexCoords(1,ivt2)-p_DvertexCoords(1,ivt1))**2 &
            + (p_DvertexCoords(2,ivt2)-p_DvertexCoords(2,ivt1))**2 )
      dedgeweight = dedgelength * 0.5_DP

      Dval = 0.0_DP
      
      ! Now we have the values of the basis functions in all the cubature
      ! points.
      !
      ! Integrate the jump over the edges. This calculates the local matrix.
      !
      ! Loop through the test basis functions
      do IDOFE = 1,ndof

          ! Loop through the cubature points to calculate the integral
          ! of the jump. Note that for computing the jump, we have to
          ! look in the inverse order to the cubature points of the neighbour
          ! element!
          ! Remember that the values of the basis functions on the first element
          ! are in p_DbasXXXX (.,.,.,1), while those of the 2nd element are
          ! in p_DbasXXXX (.,.,.,3) (rather than in p_DbasXXXX (.,.,.,2))
          ! by the above construction!
          do icubp = 1,ncubp

            ! [ grad psi ]   ( jump in the derivative of test basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dphidx = p_Du1(Idofs(IDOFE)) * Dbas (IDOFE,DER_DERIV_X,icubp,1)
            dpsidx = p_Du1(Idofs(IDOFE)) * Dbas (IDOFE,DER_DERIV_X,ncubp-icubp+1,3)

            dphidy = p_Du2(Idofs(IDOFE)) * Dbas (IDOFE,DER_DERIV_Y,icubp,1)
            dpsidy = p_Du2(Idofs(IDOFE)) * Dbas (IDOFE,DER_DERIV_Y,ncubp-icubp+1,3)

            ! Compute int_edge ( [grad phi_i] [grad phi_j] )
            dom = Domega(icubp) * dedgeweight
            !dval = dval + Domega(icubp) * dedgeweight * &
            !              (dphidx*dpsidx + dphidy*dpsidy)
            
            Dval(1) = Dval(1) + dom * ((dphidx-dpsidx)+(dphidy-dpsidy))**2
            Dval(2) = Dval(2) + dom * 0.5_DP*(dphidx+dpsidx+dphidy+dpsidy)

          end do

      end do
      
      p_Dvj(IMT) = sqrt(Dval(1))
      p_Dvm(IMT) = Dval(2)

    end do ! IMTidx

    ! Clean up allocated arrays and memory.
    deallocate(DcubPtsRefOnAllEdges)

    call elprep_releaseElementSet(revalElementSet)
    deallocate(p_DcubPtsRef)

    deallocate(Dbas)

    !deallocate(Dentry)

    deallocate(IdofsTempl)
    
    ! **********************************************************************************
    
    allocate(p_Dwj(NVT))
    allocate(p_Dwm(NVT))
    
    call spdp_projectToVertices(rvecJump, p_Dwj)
    call spdp_projectToVertices(rvecMean, p_Dwm)

    call ucd_startVTK (rexport, UCD_FLAG_STANDARD, p_rtriangulation, './ucd/' // trim(sfilename) // '.vtk')
    call ucd_addVariableVertexBased(rexport, 'div-jump', p_Dwj)
    call ucd_addVariableVertexBased(rexport, 'div-mean', p_Dwm)
    call ucd_write(rexport)
    call ucd_release(rexport)
    
    deallocate(p_Dwj)
    deallocate(p_Dwm)
    
    call vecio_writeVectorHR (rvecJump, 'div-jump', .false., 0, './ucd/' // trim(sfilename) // '_jump.txt', '(E20.12)')
    call vecio_writeVectorHR (rvecMean, 'div-mean', .false., 0, './ucd/' // trim(sfilename) // '_mean.txt', '(E20.12)')
    
    call lsyssc_releaseVector(rvecJump)
    call lsyssc_releaseVector(rvecMean)
    call spdiscr_releaseDiscr(rdiscQ1T)
    
  end subroutine

end module
