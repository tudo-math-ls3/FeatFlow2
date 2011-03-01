!##############################################################################
!# ****************************************************************************
!# <name> transport_callback1d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve scalar conservation laws in 1D.
!#
!# The following routines are available:
!#
!# 1.) transp_hadaptCallback1d
!#      -> Performs application specific tasks in the adaptation algorithm in 1D
!#
!#
!# ****************************************************************************
!#
!# The following routines for linear velocity case are available:
!#
!# 1.) transp_calcMatDiagConvP1d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 1D (primal formulation)
!#
!# 2.) transp_calcMatGalConvP1d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 1D (primal formulation)
!4
!# 3.) transp_calcMatUpwConvP1d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 1D (primal formulation)
!#
!# 4.) transp_calcMatDiagConvD1d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 1D (dual formulation)
!#
!# 5.) transp_calcMatGalConvD1d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 1D (dual formulation)
!#
!# 6.) transp_calcMatUpwConvD1d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 1D (dual formulation)
!#
!# 7.) transp_coeffVecBdrConvP1d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 1D (primal formulation)
!#
!# 8.) transp_coeffMatBdrConvP1d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 1D (primal formulation)
!#
!# 9.) transp_coeffVecBdrConvD1d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 1D (dual formulation)
!#
!# 10.) transp_coeffMatBdrConvD1d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 1D (dual formulation)
!#
!#
!# ****************************************************************************
!#
!# The following routines for Burgers` equation
!# in space-time are available:
!#
!# 1.) transp_calcMatDiagBurgP1d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for Burger`s equation in 1D (primal formulation)
!#
!# 2.) transp_calcMatGalBurgP1d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        for Burger`s equation in 1D (primal formulation)
!#
!# 3.) transp_calcMatUpwBurgP1d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for Burger`s equation in 1D (primal formulation)
!#
!# 4.) transp_coeffVecBdrBurgP1d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 2D (primal formulation)
!#
!# 5.) transp_coeffMatBdrBurgP1d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 2D (primal formulation)
!#
!#
!# ****************************************************************************
!#
!# The following routines for the Buckley-Leverett
!#  equation in space-time are available:
!#
!# 1.) transp_calcMatDiagBLevP1d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for Buckley-Leverett equation in 1D (primal formulation)
!#
!# 2.) transp_calcMatGalBLevP1d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        for Buckley-Leverett equation in 1D (primal formulation)
!#
!# 3.) transp_calcMatUpwBLevP1d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for Buckley-Leverett equation in 1D (primal formulation)
!#
!# 4.) transp_coeffVecBdrBLevP1d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 2D (primal formulation)
!#
!# 5.) transp_coeffMatBdrBLevP1d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 2D (primal formulation)
!#
!#!# </purpose>
!##############################################################################

module transport_callback1d

  use boundarycondaux
  use collection
  use derivatives
  use domainintegration
  use feevaluation
  use flagship_callback
  use fparser
  use fsystem
  use genoutput
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use scalarpde
  use spatialdiscretisation
  use storage
  use transport_basic

  implicit none

  private
  
  public :: transp_hadaptCallback1d

  public :: transp_calcMatDiagConvP1d_sim
  public :: transp_calcMatGalConvP1d_sim
  public :: transp_calcMatUpwConvP1d_sim
  public :: transp_coeffMatBdrConvP1d_sim
  public :: transp_coeffVecBdrConvP1d_sim

  public :: transp_calcMatDiagConvD1d_sim
  public :: transp_calcMatGalConvD1d_sim
  public :: transp_calcMatUpwConvD1d_sim
  public :: transp_coeffMatBdrConvD1d_sim
  public :: transp_coeffVecBdrConvD1d_sim

  public :: transp_calcMatDiagBurgP1d_sim
  public :: transp_calcMatGalBurgP1d_sim
  public :: transp_calcMatUpwBurgP1d_sim
  public :: transp_coeffVecBdrBurgP1d_sim
  public :: transp_coeffMatBdrBurgP1d_sim

  public :: transp_calcMatDiagBLevP1d_sim
  public :: transp_calcMatGalBLevP1d_sim
  public :: transp_calcMatUpwBLevP1d_sim
  public :: transp_coeffVecBdrBLevP1d_sim
  public :: transp_coeffMatBdrBLevP1d_sim

contains
  
  !*****************************************************************************

!<subroutine>

  subroutine transp_hadaptCallback1d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 1D.
!</description>

!<input>
    ! Identifier for the grid modification operation
    integer, intent(in) :: iOperation
!</input>

!<inputoutput>
    ! A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   IquickAccess(1):     NEQ or ivt
    !   IquickAccess(2:3):   ivt1,ivt2
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution


    ! What operation should be performed?
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vector from colletion and set pointer
      rsolution => rcollection%p_rvectorQuickAccess1
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_ADJUSTVERTEXDIM)
      ! Resize solution vector
      if (rsolution%NEQ .ne. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if


    case(HADAPT_OPR_INSERTVERTEXEDGE)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      p_Dsolution(rcollection%IquickAccess(1)) = &
          0.5_DP*(p_Dsolution(rcollection%IquickAccess(2))+&
                  p_Dsolution(rcollection%IquickAccess(3)))

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        p_Dsolution(rcollection%IquickAccess(1)) = &
            p_Dsolution(rcollection%IquickAccess(2))
      else
        p_Dsolution(rcollection%IquickAccess(1)) = 0.0_DP
      end if

      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback1d(iOperation, rcollection)

    end select

  end subroutine transp_hadaptCallback1d

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatDiagConvP1d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 1D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtNode
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)

    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)
#else
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = -dscale*&
          p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)
#endif
    end do
    
  end subroutine transp_calcMatDiagConvP1d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatGalConvP1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 1D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)
      ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)
#else
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)
#endif

      ! Set artificial diffusion to zero
      DcoefficientsAtEdge(1,iedge) = 0.0_DP
    end do

  end subroutine transp_calcMatGalConvP1d_sim
  
  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatUpwConvP1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 1D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its doubledata
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)
      ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)
#else
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)
#endif
      
      ! Compute artificial diffusion coefficient $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
      DcoefficientsAtEdge(1,iedge) =&
          max(-DcoefficientsAtEdge(2,iedge), 0.0_DP, -DcoefficientsAtEdge(3,iedge))
    end do

  end subroutine transp_calcMatUpwConvP1d_sim

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrConvP1d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates.  According
    ! to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! linear form the corresponding coefficients in front of the
    ! terms.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation in 1D.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for periodic boundary conditions
    !   IquickAccess(3):     number of the mirror boundary component
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution,p_rvelocity
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnx,dnv,dtime,dscale,dval
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrConvP1d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first two quick access vectors
    ! point to the solution and velocity vector (if any)
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.
      Dcoefficients = 0.0_DP

      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrConvP1d_sim')


    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form (if any).

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the Neumann values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,1).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Dcoefficients(1,ipoint,iel))

          ! Multiply by scaling coefficient
          Dcoefficients(1,ipoint,iel) = dscale * Dcoefficients(1,ipoint,iel)
        end do
      end do

      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow ({\bf v}u)\cdot{\bf n}=({\bf v}g)\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser for Dirichlet value
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)

          ! Impose Dirichlet value via penalty method
          Dcoefficients(1,ipoint,iel) = dscale * dval * BDRC_DIRICHLET_PENALTY
        end do
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form 
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.
      
      if (associated(p_rvelocity)) then

        ! Allocate temporal memory
        allocate(Daux(npointsPerElement,nelements,2))
        
        ! Evaluate the velocity field in the cubature points on the boundary
        ! and store the result in Daux(:,:,:,1)
        call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Initialize values
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime
        
        ! Evaluate the function parser for the boundary values in the
        ! cubature points on the boundary and store the result in
        ! Dcoefficients(:,:,2).
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Set values for function parser
            Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment,&
                Dvalue, Daux(ipoint,iel,2))
          end do
        end do
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
            
            ! Compute the normal velocity and impose Dirichlet boundary condition
            dnv = dnx * Daux(ipoint,iel,1)
            Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,2)
          end do
        end do
        
        ! Deallocate temporal memory
        deallocate(Daux)

      else
        
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
        
      end if
      

    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      if (associated(p_rvelocity)) then
        
        ! Allocate temporal memory
        allocate(Daux(npointsPerElement,nelements,2))
        
        ! Evaluate the velocity field in the cubature points on the boundary
        ! and store the result in Daux(:,:,:,1)
        call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Initialize values
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime
        
        ! Evaluate the function parser for the boundary values in the
        ! cubature points on the boundary and store the result in
        ! Dcoefficients(:,:,2).
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Set values for function parser
            Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment,&
                Dvalue, Daux(ipoint,iel,2))
          end do
        end do
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
            
            ! Compute the normal velocity
            dnv = dnx * Daux(ipoint,iel,1)
            
            ! Check if we are at the primal inflow boundary
            if (dnv .lt. 0.0_DP) then
              Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,2)
            else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
            end if
          end do
        end do
        
        ! Deallocate temporal memory
        deallocate(Daux)

      else
        
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
        
      end if


    case(BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Periodic/Antiperiodic boundary conditions (Flux boundary conditions):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.
      
      if (associated(p_rvelocity)) then
        
        ! Allocate temporal memory
        allocate(Daux(npointsPerElement,nelements,2))
        
        ! Evaluate the velocity field in the cubature points on the boundary
        ! and store the result in Daux(:,:,:,1)
        call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

        ! Evaluate the solution in the cubature points on the mirrored
        ! boundary and store the result in Daux(:,:,2)
        call doEvaluateAtBdr1d(DER_FUNC1D, npointsPerElement*nelements, Daux(:,:,2),&
            p_rsolution%RvectorBlock(1), rcollection%IquickAccess(3))
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
            
            ! Compute the normal velocity
            dnv = dnx * Daux(ipoint,iel,1)
            
            ! Check if we are at the primal inflow boundary
            if (dnv .lt. 0.0_DP) then
              Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,2)
            else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
            end if
          end do
        end do
        
        ! Deallocate temporal memory
        deallocate(Daux)

      else
        
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
        
      end if

      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrConvP1d_sim')
      call sys_halt()
      
    end select

  contains

    ! Here come the working routines

    !***************************************************************************
    ! Evaluate th solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    ! ***************************************************************************
    
    subroutine doEvaluateAtBdr1d(iderType, n, Dvalues, rvectorScalar, ibdc)
      
      integer, intent(in) :: iderType,ibdc,n
      type(t_vectorScalar), intent(in) :: rvectorScalar
      
      real(DP), dimension(n), intent(out) :: Dvalues

      call fevl_evaluateBdr1D(iderType, Dvalues, rvectorScalar, ibdc)
      
    end subroutine doEvaluateAtBdr1d

  end subroutine transp_coeffVecBdrConvP1d_sim

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrConvD1d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the dual problem for the
    ! convection-diffusion equation in 1D.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for periodic boundary conditions
    !   IquickAccess(3):     number of the mirror boundary component
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution,p_rvelocity
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnx,dnv,dtime,dscale,dval
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrConvD1d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first two quick access vectors
    ! point to the solution and velocity vector (if any)
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assemble of a "zero" boundary integral.
      Dcoefficients = 0.0_DP
      
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrConvD1d_sim')


    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form (if any)

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Evaluate the function parser for the Neumann values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,1).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Dcoefficients(1,ipoint,iel))
          
          ! Multiply by scaling coefficient
          Dcoefficients(1,ipoint,iel) = dscale * Dcoefficients(1,ipoint,iel)
        end do
      end do


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow ({\bf v}u)\cdot{\bf n}=({\bf v}g)\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser for Dirichlet value
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)

          ! Impose Dirichlet value via penalty method
          Dcoefficients(1,ipoint,iel) = -dscale * dval * BDRC_DIRICHLET_PENALTY
        end do
      end do
      
      
    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form 
      !
      ! $$ ({\bf v}u-d\nabla u)\cdot{\bf n}=({\bf v}g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory
        allocate(Daux(npointsPerElement, nelements, NDIM1D+1))
        
        ! Evaluate the velocity field in the cubature points on the boundary
        ! and store the result in Daux(:,:,:,1)
        call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Initialize values
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime
        
        ! Evaluate the function parser for the boundary values in the
        ! cubature points on the boundary and store the result in
        ! Dcoefficients(:,:,2).
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Set values for function parser
            Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment,&
                Dvalue, Daux(ipoint,iel,2))
          end do
        end do
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
            
            ! Compute the normal velocity and impose Dirichlet boundary condition
            dnv = dnx * Daux(ipoint,iel,1)
            Dcoefficients(1,ipoint,iel) = dscale * dnv * Daux(ipoint,iel,2)
          end do
        end do
        
        ! Deallocate temporal memory
        deallocate(Daux)
        
      else
        
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
        
      end if
        
      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ ({\bf v}u-d\nabla u)\cdot{\bf n} = ({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.
      
      if (associated(p_rvelocity)) then
        
        ! Allocate temporal memory
        allocate(Daux(npointsPerElement, nelements, NDIM1D+1))
        
        ! Evaluate the velocity field in the cubature points on the boundary
        ! and store the result in Daux(:,:,:,1)
        call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Initialize values
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime
        
        ! Evaluate the function parser for the boundary values in the
        ! cubature points on the boundary and store the result in
        ! Dcoefficients(:,:,2).
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Set values for function parser
            Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment,&
                Dvalue, Daux(ipoint,iel,2))
          end do
        end do
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
            
            ! Compute the normal velocity
            dnv = dnx * Daux(ipoint,iel,1)
            
            ! Check if we are at the dual inflow boundary
            if (dnv .gt. 0.0_DP) then
              Dcoefficients(1,ipoint,iel) = dscale * dnv * Daux(ipoint,iel,2)
            else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
            end if
          end do
        end do
        
        ! Deallocate temporal memory
        deallocate(Daux)
        
      else
        
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
        
      end if


    case(BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Periodic/Antiperiodic boundary conditions (Flux boundary conditions):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -({\bf v}u-d\nabla u)\cdot{\bf n} = -({\bf v}g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      if (associated(p_rvelocity)) then
        
        ! Allocate temporal memory
        allocate(Daux(npointsPerElement, nelements, NDIM1D+1))
        
        ! Evaluate the velocity field in the cubature points on the boundary
        ! and store the result in Daux(:,:,:,1)
        call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

        ! Evaluate the solution in the cubature points on the mirrored
        ! boundary and store the result in Daux(:,:,2)
        call doEvaluateAtBdr1d(DER_FUNC1D, npointsPerElement*nelements, Daux(:,:,2),&
            p_rsolution%RvectorBlock(1), rcollection%IquickAccess(3))
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
            
            ! Compute the normal velocity
            dnv = dnx * Daux(ipoint,iel,1)
            
            ! Check if we are at the dual inflow boundary
            if (dnv .gt. 0.0_DP) then
              Dcoefficients(1,ipoint,iel) = dscale * dnv * Daux(ipoint,iel,2)
            else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
            end if
          end do
        end do
        
        ! Deallocate temporal memory
        deallocate(Daux)

      else
        
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
        
      end if

    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrConvD1d_sim')
      call sys_halt()
      
    end select

  contains

    ! Here come the working routines

    !***************************************************************************
    ! Evaluate th solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    ! ***************************************************************************
    
    subroutine doEvaluateAtBdr1d(iderType, n, Dvalues, rvectorScalar, ibdc)
      
      integer, intent(in) :: iderType,ibdc,n
      type(t_vectorScalar), intent(in) :: rvectorScalar
      
      real(DP), dimension(n), intent(out) :: Dvalues

      call fevl_evaluateBdr1D(iderType, Dvalues, rvectorScalar, ibdc)
      
    end subroutine doEvaluateAtBdr1d
    
  end subroutine transp_coeffVecBdrConvD1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrConvP1d_sim(rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, IdofsTrial, IdofsTest, rdomainIntSubset,&
      Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation in 1D.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rsolution,p_rvelocity
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP) :: dnx,dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrConvP1d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! point to the solution and velocity vector (if any)
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral (if any)

      if (associated(p_rvelocity)) then
        
        ! Allocate temporal memory
        allocate(Daux(npointsPerElement, nelements, NDIM1D+1))
        
        ! Evaluate the velocity field in the cubature points on the boundary
        ! and store the result in Daux(:,:,:,1)
        call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
            
            ! Compute the normal velocity
            dnv = dnx * Daux(ipoint,iel,1)
            
            ! Scale normal velocity by scaling parameter
            Dcoefficients(1,ipoint,iel) = -dscale * dnv
          end do
        end do
        
        ! Free temporal memory
        deallocate(Daux)

      else

        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
        
      end if


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Impose Dirichlet boundary conditions via penalty method
          Dcoefficients(1,ipoint,iel) = -dscale * BDRC_DIRICHLET_PENALTY
        end do
      end do
      

    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form      
      Dcoefficients = 0.0_DP

      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrConvP1d_sim')

      
    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, NDIM1D+1))

      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1)
      call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
          p_rvelocity%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Multiply the velocity vector with the normal in each point
      ! to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
          
          ! Compute the normal velocity
          dnv = dnx * Daux(ipoint,iel,1)

          ! Check if we are at the primal outflow boundary
          if (dnv .gt. 0.0_DP) then
            Dcoefficients(1,ipoint,iel) = -dscale * dnv
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do

      ! Free temporal memory
      deallocate(Daux)
      
    
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrConvP1d_sim')
      call sys_halt()
      
    end select
    
  end subroutine transp_coeffMatBdrConvP1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrConvD1d_sim(rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, IdofsTrial, IdofsTest, rdomainIntSubset,&
      Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the dual problem for the
    ! convection-diffusion equation in 1D.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   rvectorQuickAccess2: velocity field
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rsolution,p_rvelocity
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP) :: dnx,dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrConvD1d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! point to the solution and velocity vector (if any)
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the boundary integral for the convective term (if any)

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory
        allocate(Daux(npointsPerElement, nelements, NDIM1D+1))
        
        ! Evaluate the velocity field in the cubature points on the boundary
        ! and store the result in Daux(:,:,:,1)
        call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Get the normal vector in the point from the boundary
            dnx = merge(1.0_DP, -1.0_DP, mod(ibct,2) .eq. 0)
            
            ! Compute the normal velocity
            dnv = dnx * Daux(ipoint,iel,1)

            ! Scale normal velocity by scaling parameter
            Dcoefficients(1,ipoint,iel) = dscale * dnv
          end do
        end do
        
        ! Free temporal memory
        deallocate(Daux)

      else

        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP
        
      end if


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Impose Dirichlet boundary conditions via penalty method
          Dcoefficients(1,ipoint,iel) = dscale * BDRC_DIRICHLET_PENALTY
        end do
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Dirichlet or Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form.
      Dcoefficients = 0.0_DP

      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrConvD1d_sim')
      

    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, NDIM1D+1))

      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1)
      call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
          p_rvelocity%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Multiply the velocity vector with the normal in each point
      ! to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)

          ! Compute the normal velocity
          dnv = dnx * Daux(ipoint,iel,1)

          ! Check if we are at the dual outflow boundary
          if (dnv .lt. 0.0_DP) then
            Dcoefficients(1,ipoint,iel) = dscale * dnv
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do

      ! Free temporal memory
      deallocate(Daux)

      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrConvD1d_sim')
      call sys_halt()

    end select

  end subroutine transp_coeffMatBdrConvD1d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatDiagConvD1d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 1D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtNode
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)

    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = -dscale*&
          p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)
#else
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          p_DvelocityX(IverticesAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)
#endif
    end do
    
  end subroutine transp_calcMatDiagConvD1d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatGalConvD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 1D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>
    
    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    
    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)
      ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)
#endif

      ! Set artificial diffusion to zero
      DcoefficientsAtEdge(1,iedge) = 0.0_DP
    end do
    
  end subroutine transp_calcMatGalConvD1d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatUpwConvD1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 1D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: velocity field
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_DvelocityX
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    
    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          p_DvelocityX(IverticesAtEdge(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)
      ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          p_DvelocityX(IverticesAtEdge(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)
#endif

      ! Compute artificial diffusion coefficient $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
      DcoefficientsAtEdge(1,iedge) = dscale*&
          max(-DcoefficientsAtEdge(2,iedge), 0.0_DP, -DcoefficientsAtEdge(3,iedge))
    end do

  end subroutine transp_calcMatUpwConvD1d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatDiagBurgP1d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for the primal Burger`s equation in 1D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtNode

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>
    
    ! local variables
    integer :: inode

    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficients  $k_{ii} = 0.5*u_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          0.5_DP*DdataAtNode(inode)*DmatrixCoeffsAtNode(1,inode)
#else
      ! Compute convective coefficients  $k_{ii} = -0.5*u_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = -dscale*&
          0.5_DP*DdataAtNode(inode)*DmatrixCoeffsAtNode(1,inode)
#endif
    end do
    
  end subroutine transp_calcMatDiagBurgP1d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalBurgP1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the primal Burger`s equation in 1D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>
    
    ! local variables
    integer :: iedge

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient  $k_{ij} = 0.5*u_j*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          0.5_DP*DdataAtEdge(2,iedge)*DmatrixCoeffsAtEdge(1,2,iedge)
          
      ! Compute convective coefficient  $k_{ji} = 0.5*u_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          0.5_DP*DdataAtEdge(1,iedge)*DmatrixCoeffsAtEdge(1,1,iedge)
#else
      ! Compute convective coefficient  $k_{ij} = -0.5*u_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          0.5_DP*DdataAtEdge(2,iedge)*DmatrixCoeffsAtEdge(1,1,iedge)
          
      ! Compute convective coefficient  $k_{ji} = -0.5*u_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          0.5_DP*DdataAtEdge(1,iedge)*DmatrixCoeffsAtEdge(1,2,iedge)
#endif
      
      ! Set artificial diffusion to zero
      DcoefficientsAtEdge(1,iedge) = 0.0_DP
    end do

  end subroutine transp_calcMatGalBurgP1d_sim
  
  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwBurgP1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the primal Burger`s equation in 1D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variables
    integer :: iedge

    do iedge  = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient  $k_{ij} = 0.5*u_j*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          0.5_DP*DdataAtEdge(2,iedge)*DmatrixCoeffsAtEdge(1,2,iedge)
          
      ! Compute convective coefficient  $k_{ji} = 0.5*u_i*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          0.5_DP*DdataAtEdge(1,iedge)*DmatrixCoeffsAtEdge(1,1,iedge)
#else
      ! Compute convective coefficient  $k_{ij} = -0.5*u_j*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          0.5_DP*DdataAtEdge(2,iedge)*DmatrixCoeffsAtEdge(1,1,iedge)

      ! Compute convective coefficient  $k_{ji} = -0.5*u_i*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          0.5_DP*DdataAtEdge(1,iedge)*DmatrixCoeffsAtEdge(1,2,iedge)
#endif
      
      ! Compute artificial diffusion coefficient 
      ! $d_{ij} = \max\{abs(v_{ij}*0.5*(c_{ij}-c_{ji})\}$,
      ! where $v_{ij} = 0.5*(f(u_j)-f(u_i))/(u_j-u_i) = 0.5*(u_i+u_j)$
      DcoefficientsAtEdge(1,iedge) =&
          abs(0.25_DP*(DdataAtEdge(1,iedge)+DdataAtEdge(2,iedge))*&
              (DmatrixCoeffsAtEdge(1,1,iedge)-DmatrixCoeffsAtEdge(1,2,iedge)))

    end do

  end subroutine transp_calcMatUpwBurgP1d_sim

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrBurgP1d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates.  According
    ! to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! linear form the corresponding coefficients in front of the
    ! terms.
    !
    ! This routine handles the primal problem for the
    ! Burgers`s equation in 1D.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for periodic boundary conditions
    !   IquickAccess(3):     number of the mirror boundary component
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnx,dnv,dtime,dscale,dval
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrBurgP1d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.
      Dcoefficients = 0.0_DP

      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrBurgP1d_sim')


    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the Neumann values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,1).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Dcoefficients(1,ipoint,iel))

          ! Multiply by scaling coefficient
          Dcoefficients(1,ipoint,iel) = dscale * Dcoefficients(1,ipoint,iel)
        end do
      end do

      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow [0.5*u]*u\cdot{\bf n}=[0.5*g]*g\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser for Dirichlet value
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)

          ! Impose Dirichlet value via penalty method
          Dcoefficients(1,ipoint,iel) = dscale * dval * BDRC_DIRICHLET_PENALTY
        end do
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form 
      !
      ! $$ -([0.5*u]*u-d\nabla u)\cdot{\bf n} = -([0.5*g]*g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.
      
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 1))

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,1).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,1))
        end do
      end do
      
      ! Multiply the velocity vector [0.5*u] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)
          
          ! Compute the normal velocity and impose Dirichlet boundary condition
          dnv = dnx * 0.5_DP*Daux(ipoint,iel,1)
          Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,1)
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux)

      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([0.5*u]*u-d\nabla u)\cdot{\bf n} = -([0.5*g]*g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 2))

      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1)
      call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,2).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,2))
        end do
      end do

      ! Multiply the velocity vector [0.5*u] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)

          ! Compute the normal velocity
          dnv = dnx * 0.5_DP*Daux(ipoint,iel,1)

          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
            Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,2)
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do
      
      ! Deallocate temporal memory
      deallocate(Daux)


    case(BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Periodic/Antiperiodic boundary conditions (Flux boundary conditions):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([0.5*u]*u-d\nabla u)\cdot{\bf n} = -([0.5*g]*g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.
      
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement,nelements,2))
      
      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1)
      call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Evaluate the solution in the cubature points on the mirrored
      ! boundary and store the result in Daux(:,:,2)
      call doEvaluateAtBdr1d(DER_FUNC1D, npointsPerElement*nelements, Daux(:,:,2),&
          p_rsolution%RvectorBlock(1), rcollection%IquickAccess(3))
      
      ! Multiply the velocity vector [0.5*u] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)
          
          ! Compute the normal velocity
          dnv = dnx * 0.5_DP*(Daux(ipoint,iel,1)+Daux(ipoint,iel,2))
          
          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
            Dcoefficients(1,ipoint,iel) = -dscale * dnx*0.5_DP*Daux(ipoint,iel,2)**2
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do
      
      ! Deallocate temporal memory
      deallocate(Daux)
      
      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrBurgP1d_sim')
      call sys_halt()
      
    end select

  contains

    ! Here come the working routines

    !***************************************************************************
    ! Evaluate th solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    ! ***************************************************************************
    
    subroutine doEvaluateAtBdr1d(iderType, n, Dvalues, rvectorScalar, ibdc)
      
      integer, intent(in) :: iderType,ibdc,n
      type(t_vectorScalar), intent(in) :: rvectorScalar
      
      real(DP), dimension(n), intent(out) :: Dvalues

      call fevl_evaluateBdr1D(iderType, Dvalues, rvectorScalar, ibdc)
      
    end subroutine doEvaluateAtBdr1d
    
  end subroutine transp_coeffVecBdrBurgP1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrBurgP1d_sim(rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, IdofsTrial, IdofsTest, rdomainIntSubset,&
      Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the primal problem for the
    ! Burger`s equation in 1D.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP) :: dnx,dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrBurgP1d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral
    
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 1))
      
      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1)
      call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Multiply the velocity vector [0.5*u] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)
          
          ! Compute the normal velocity
          dnv = dnx * 0.5_DP*Daux(ipoint,iel,1)
          
          ! Scale normal velocity by scaling parameter
          Dcoefficients(1,ipoint,iel) = -dscale * dnv
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux)


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Impose Dirichlet boundary conditions via penalty method
          Dcoefficients(1,ipoint,iel) = -dscale * BDRC_DIRICHLET_PENALTY
        end do
      end do
      

    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form      
      Dcoefficients = 0.0_DP

      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrBurgP1d_sim')

      
    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 1))

      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1)
      call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Multiply the velocity vector [0.5*u] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Get the normal vector in the point from the boundary
            dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)
          
          ! Compute the normal velocity
          dnv = dnx * 0.5_DP*Daux(ipoint,iel,1)

          ! Check if we are at the primal outflow boundary
          if (dnv .gt. 0.0_DP) then
            Dcoefficients(1,ipoint,iel) = -dscale * dnv
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do

      ! Free temporal memory
      deallocate(Daux)
      
    
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrBurgP1d_sim')
      call sys_halt()
      
    end select
    
  end subroutine transp_coeffMatBdrBurgP1d_sim
  
  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatDiagBLevP1d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IverticesAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ for the primal Buckley-Leverett equation
    ! $du/dt+df(u)/dx=0$ in 1D, whereby the flux function is given by
    ! $f(u)=u^2/(u^2+0.5*(1-u)^2)$
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtNode

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of nodes
    integer, intent(in) :: nnodes
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all nodes under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtNode
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui
    integer :: inode
    
    do inode = 1, nnodes

      ui = DdataAtNode(inode)

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient  $k_{ii} = a_i*C_{ii}$,
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtNode(1,inode) = dscale*&
          ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtNode(1,inode)
#else
      ! Compute convective coefficient  $k_{ii} = -a_i*C_{ii}$,
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtNode(1,inode) = -dscale*&
          ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtNode(1,inode)
#endif

    end do
    
  end subroutine transp_calcMatDiagBLevP1d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalBLevP1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the primal Buckley-Leverett equation
    ! $du/dt+df(u)/dx=0$ in 1D, whereby the flux function is given by
    ! $f(u)=u^2/(u^2+0.5*(1-u)^2)$
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,uj
    integer :: iedge
    
    do iedge = 1, nedges

      ui = DdataAtEdge(1,iedge); uj = DdataAtEdge(2,iedge)

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient  $k_{ij} = a_j*C_{ji}$,
      ! where $a_j=u_j/(u_^2+0.5*(1-u_j)^2)$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DmatrixCoeffsAtEdge(1,2,iedge)

      ! Compute convective coefficient  $k_{ji} = a_i*C_{ij}$,
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtEdge(1,1,iedge)
#else
      ! Compute convective coefficient  $k_{ij} = -a_j*C_{ij}$,
      ! where $a_j=u_j/(u_j^2+0.5*(1-u_j)^2)$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DmatrixCoeffsAtEdge(1,1,iedge)

      ! Compute convective coefficient  $k_{ji} = -a_i*C_{ji}$,
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtEdge(1,2,iedge)
#endif

      ! Set artificial diffusion to zero
      DcoefficientsAtEdge(1,iedge) = 0.0_DP
    end do
    
  end subroutine transp_calcMatGalBLevP1d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwBLevP1d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IverticesAtEdge, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the primal Buckley-Leverett equation
    ! $du/dt+df(u)/dx=0$ in 1D, whereby the flux function is given by
    ! $f(u)=u^2/(u^2+0.5*(1-u)^2)$. 
    !
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Number of edges
    integer, intent(in) :: nedges
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! Coefficients of the matrix for all edges under consideration
    real(DP), dimension(:,:), intent(out) :: DcoefficientsAtEdge
!</output>
!</subroutine>

    ! local variable
    real(DP) :: ui,uj,vij
    integer :: iedge
    
    do iedge = 1, nedges

      ui = DdataAtEdge(1,iedge); uj = DdataAtEdge(2,iedge)

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient  $k_{ij} = a_j*C_{ji}$,
      ! where $a_j=u_j/(u_^2+0.5*(1-u_j)^2)$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DmatrixCoeffsAtEdge(1,2,iedge)

      ! Compute convective coefficient  $k_{ji} = a_i*C_{ij}$,
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtEdge(1,1,iedge)
#else
      ! Compute convective coefficient  $k_{ij} = -a_j*C_{ij}$,
      ! where $a_j=u_j/(u_j^2+0.5*(1-u_j)^2)$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DmatrixCoeffsAtEdge(1,1,iedge)

      ! Compute convective coefficient  $k_{ji} = -a_i*C_{ji}$,
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtEdge(1,2,iedge)
#endif

      ! Calculate the characteristic speed
      if (abs(ui-uj) .gt. SYS_EPSREAL) then
        vij = (uj*uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))&
              -ui*ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui)))/(uj-ui)
      else
        vij = ui*(1.0_DP-ui)/(ui*ui-0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))**2
      end if

      ! Compute artificial diffusion coefficient 
      ! $d_{ij} = \max\{abs(0.5*v_{ij}*(c_{ij}-c_{ji})\}$      
      DcoefficientsAtEdge(1,iedge) =&
          abs(0.5_DP*vij*(DmatrixCoeffsAtEdge(1,1,iedge)-DmatrixCoeffsAtEdge(1,2,iedge)))
    end do
    
  end subroutine transp_calcMatUpwBLevP1d_sim

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrBLevP1d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates.  According
    ! to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! linear form the corresponding coefficients in front of the
    ! terms.
    !
    ! This routine handles the primal problem for the
    ! Buckley-Leverett equation in 1D.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for periodic boundary conditions
    !   IquickAccess(3):     number of the mirror boundary component
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnx,dnv,dtime,dscale,dval
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrBLevP1d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Homogeneous Neumann boundary conditions:
      !
      ! The diffusive part in the linear form vanishes since
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form.
      !
      ! Hence, this routine should not be called for homogeneous
      ! Neumann boundary conditions since it corresponds to an
      ! expensive assembly of a "zero" boundary integral.
      Dcoefficients = 0.0_DP

      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrBLevP1d_sim')


    case (BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! Inhomogeneous Neumann boundary conditions:
      !
      ! Evaluate coefficient for the diffusive part of the linear form
      !
      ! $$ d\nabla u\cdot{\bf n}=0 $$
      !
      ! The convective part is included into the bilinear form (if any).

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the Neumann values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,1).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Dcoefficients(1,ipoint,iel))

          ! Multiply by scaling coefficient
          Dcoefficients(1,ipoint,iel) = dscale * Dcoefficients(1,ipoint,iel)
        end do
      end do

      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      !
      ! Evaluate coefficient for the convective part of the linear form
      !
      ! $$ u=g \Rightarrow [a]*u\cdot{\bf n}=[a]*g\cdot{\bf n} $$
      !
      ! where $ a = u/(u^2+0.5*(1-u)^2) $ is the velocity
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser for Dirichlet value
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)

          ! Impose Dirichlet value via penalty method
          Dcoefficients(1,ipoint,iel) = dscale * dval * BDRC_DIRICHLET_PENALTY
        end do
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form 
      !
      ! $$ -([a]*u-d\nabla u)\cdot{\bf n} = -([a]*g)\cdot{\bf n} $$
      !
      ! where $ a = u/(u^2+0.5*(1-u)^2) $ is the velocity
      !
      ! and do not include any boundary integral into the bilinear form at all.
      
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 1))

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,2).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,1))
        end do
      end do
      
      ! Multiply the velocity vector [a] with the normal in each point
      ! to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)
          
          ! Compute the normal velocity and impose Dirichlet boundary condition
          dnv = dnx * Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2&
                      + 0.5_DP*(1-Daux(ipoint,iel,1))**2)
          Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,1)
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux)

      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([a]*u-d\nabla u)\cdot{\bf n} = -([a]*g)\cdot{\bf n} $$
      !
      ! where $ a = u/(u^2+0.5*(1-u)^2) $ is the velocity
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 2))

      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1)
      call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,2).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM1D) = Dpoints(1:NDIM1D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,2))
        end do
      end do

      ! Multiply the velocity vector [a] with the normal in each point
      ! to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)

          ! Compute the normal velocity
          dnv = dnx * Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2&
                      + 0.5_DP*(1-Daux(ipoint,iel,1))**2)

          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
            ! Compute the prescribed normal velocity
            dnv = dnx * Daux(ipoint,iel,2)/(Daux(ipoint,iel,2)**2&
                        + 0.5_DP*(1-Daux(ipoint,iel,2))**2)
            Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,2)
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux)


    case(BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Periodic/Antiperiodic boundary conditions (Flux boundary conditions):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([a]*u-d\nabla u)\cdot{\bf n} = -([a]*g)\cdot{\bf n} $$
      !
      ! where $ a = u/(u^2+0.5*(1-u)^2) $ is the velocity
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 2))

      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1)
      call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Evaluate the solution in the cubature points on the mirrored
      ! boundary and store the result in Daux(:,:,2)
      call doEvaluateAtBdr1d(DER_FUNC1D, npointsPerElement*nelements, Daux(:,:,2),&
          p_rsolution%RvectorBlock(1), rcollection%IquickAccess(3))
      
      ! Multiply the velocity vector [a] with the normal in each point
      ! to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)

          ! Compute the normal velocity
          dnv = dnx * Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2&
                      + 0.5_DP*(1-Daux(ipoint,iel,1))**2)

          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
            ! Compute the prescribed normal velocity
            dnv = dnx * Daux(ipoint,iel,2)/(Daux(ipoint,iel,2)**2&
                        + 0.5_DP*(1-Daux(ipoint,iel,2))**2)
            Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,2)
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux)

      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrBLevP1d_sim')
      call sys_halt()
      
    end select

  contains

    ! Here come the working routines

    !***************************************************************************
    ! Evaluate th solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    ! ***************************************************************************
    
    subroutine doEvaluateAtBdr1d(iderType, n, Dvalues, rvectorScalar, ibdc)
      
      integer, intent(in) :: iderType,ibdc,n
      type(t_vectorScalar), intent(in) :: rvectorScalar
      
      real(DP), dimension(n), intent(out) :: Dvalues

      call fevl_evaluateBdr1D(iderType, Dvalues, rvectorScalar, ibdc)
      
    end subroutine doEvaluateAtBdr1d

  end subroutine transp_coeffVecBdrBLevP1d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrBLevP1d_sim(rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, IdofsTrial, IdofsTest, rdomainIntSubset,&
      Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
    !
    ! This routine handles the primal problem for the
    ! Buckley-Leverett equation in 1D.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTrial

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in) :: rform

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTrial

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   DquickAccess(1):     simulation time
    !   DquickAccess(2):     scaling parameter
    !   IquickAccess(1):     boundary type
    !   IquickAccess(2):     segment number
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rsolution
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP) :: dnx,dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrBLevP1d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral (if any)
      
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 2))
      
      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1)
      call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Multiply the velocity vector [a] with the normal in each point
      ! to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)
          
          ! Compute the normal velocity
          dnv = dnx * Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2&
                      + 0.5_DP*(1-Daux(ipoint,iel,1))**2)
          
          ! Scale normal velocity by scaling parameter
          Dcoefficients(1,ipoint,iel) = -dscale * dnv
        end do
      end do
        
      ! Free temporal memory
      deallocate(Daux)
      

    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Impose Dirichlet boundary conditions via penalty method
          Dcoefficients(1,ipoint,iel) = -dscale * BDRC_DIRICHLET_PENALTY
        end do
      end do
      

    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form      
      Dcoefficients = 0.0_DP

      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrBLevP1d_sim')

      
    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 2))

      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1)
      call fevl_evaluate_sim(DER_FUNC1D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Multiply the velocity vector [a] with the normal in each point
      ! to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Get the normal vector in the point from the boundary
          dnx = merge(1.0, -1.0, mod(ibct,2) .eq. 0)
          
          ! Compute the normal velocity
          dnv = dnx * Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2&
                      + 0.5_DP*(1-Daux(ipoint,iel,1))**2)  

          ! Check if we are at the primal outflow boundary
          if (dnv .gt. 0.0_DP) then
            Dcoefficients(1,ipoint,iel) = -dscale * dnv
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do

      ! Free temporal memory
      deallocate(Daux)
      
    
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrBLevP1d_sim')
      call sys_halt()
      
    end select
    
  end subroutine transp_coeffMatBdrBLevP1d_sim

end module transport_callback1d
