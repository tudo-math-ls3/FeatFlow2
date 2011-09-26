!##############################################################################
!# ****************************************************************************
!# <name> transport_callback2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve scalar conservation laws in 2D.
!#
!# The following routines are available:
!#
!# 1.) transp_hadaptCallback2d
!#     -> Performs application specific tasks in the adaptation algorithm in 2D
!#
!# 2.) transp_refFuncBdrInt2d_sim
!#     -> Callback routine for the evaluation of the boundary integral
!#        of the target functional for goal-oriented error estimation
!#
!# 3.) transp_errorBdrInt2d_sim
!#     -> Callback routine for the evaluation of the boundary integral
!#        of the error in the target functional for goal-oriented
!#        error estimation
!#
!# 4.) transp_weightFuncBdrInt2d_sim
!#     -> Callback routine for the evaluation of the weights in
!#        the boundary integral of the target functional for
!#        goal-oriented error estimation
!#
!#
!# ****************************************************************************
!#
!# The following routines for linear velocity case are available:
!#
!# 1.) transp_calcMatDiagConvP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 2D (primal formulation)
!#
!# 2.) transp_calcMatGalConvP2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 2D (primal formulation)
!#
!# 3.) transp_calcMatUpwConvP2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 2D (primal formulation)
!#
!# 4.) transp_calcMatDiagConvD2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for linear convection in 2D (dual formulation)
!#
!# 5.) transp_calcMatGalConvD2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        for linear convection in 2D (dual formulation)
!#
!# 6.) transp_calcMatUpwConvD2d_sim
!#     -> Calculates the off-diagonal Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for linear convection in 2D (dual formulation)
!#
!# 7.) transp_coeffVecBdrConvP2d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 2D (primal formulation)
!#
!# 8.) transp_coeffMatBdrConvP2d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 2D (primal formulation)
!#
!# 9.) transp_coeffVecBdrConvD2d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 2D (dual formulation)
!#
!# 10.) transp_coeffMatBdrConvD2d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 2D (dual formulation)
!#
!# ****************************************************************************
!#
!# The following routines for Burgers` equation in space-time are available:
!#
!# 1.) transp_calcMatDiagSTBurgP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 2.) transp_calcMatGalSTBurgP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 3.) transp_calcMatUpwSTBurgP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 4.) transp_coeffVecBdrSTBurgP2d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 2D (primal formulation)
!#
!# 5.) transp_coeffMatBdrSTBurgP2d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 2D (primal formulation)
!#
!# ****************************************************************************
!#
!# The following routines for the Buckley-Leverett equation in
!# space-time are available:
!#
!# 1.) transp_calcMatDiagSTBLevP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for Buckley-Leverett equation in 2D (primal formulation)
!#
!# 2.) transp_calcMatGalSTBLevP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        for Buckley-Leverett equation in 2D (primal formulation)
!#
!# 3.) transp_calcMatUpwSTBLevP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for Buckley-Leverett equation in 2D (primal formulation)
!#
!# 4.) transp_coeffVecBdrSTBLevP2d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 2D (primal formulation)
!#
!# 5.) transp_coeffMatBdrSTBLevP2d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 2D (primal formulation)
!#
!# ****************************************************************************
!#
!# The following routines for the Burgers` equation in 2D are available:
!#
!# 1.) transp_calcMatDiagBurgP2d_sim
!#     -> Calculates the diagonal Galerkin transport coefficients
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 2.) transp_calcMatGalBurgP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 3.) transp_calcMatUpwBurgP2d_sim
!#     -> Calculates the Galerkin transport coefficients
!#        and applies scalar artificial diffusion (discrete upwinding)
!#        for Burger`s equation in 2D (primal formulation)
!#
!# 4.) transp_coeffVecBdrBurgP2d_sim
!#      -> Calculates the coefficients for the linear form
!#         in 2D (primal formulation)
!#
!# 5.) transp_coeffMatBdrBurgP2d_sim
!#     -> Calculates the coefficients for the bilinear form
!#        in 2D (primal formulation)
!#
!# </purpose>
!##############################################################################

module transport_callback2d

  use basicgeometry
  use boundary
  use boundarycondaux
  use collection
  use cubature
  use derivatives
  use domainintegration
  use element
  use feevaluation
  use flagship_callback
  use fparser
  use fsystem
  use genoutput
  use hadaptaux
  use linearsystemblock
  use linearsystemscalar
  use mprimitives
  use scalarpde
  use spatialdiscretisation
  use storage
  use transport_basic

  implicit none

  private

  public :: transp_hadaptCallback2d
  public :: transp_refFuncBdrInt2d_sim
  public :: transp_errorBdrInt2d_sim
  public :: transp_weightFuncBdrInt2d_sim

  public :: transp_calcMatDiagConvP2d_sim
  public :: transp_calcMatGalConvP2d_sim
  public :: transp_calcMatUpwConvP2d_sim
  public :: transp_coeffMatBdrConvP2d_sim
  public :: transp_coeffVecBdrConvP2d_sim

  public :: transp_calcMatDiagConvD2d_sim
  public :: transp_calcMatGalConvD2d_sim
  public :: transp_calcMatUpwConvD2d_sim
  public :: transp_coeffMatBdrConvD2d_sim
  public :: transp_coeffVecBdrConvD2d_sim

  public :: transp_calcMatDiagSTBurgP2d_sim
  public :: transp_calcMatGalSTBurgP2d_sim
  public :: transp_calcMatUpwSTBurgP2d_sim
  public :: transp_coeffMatBdrSTBurgP2d_sim
  public :: transp_coeffVecBdrSTBurgP2d_sim
  
  public :: transp_calcMatDiagSTBLevP2d_sim
  public :: transp_calcMatGalSTBLevP2d_sim
  public :: transp_calcMatUpwSTBLevP2d_sim
  public :: transp_coeffMatBdrSTBLevP2d_sim
  public :: transp_coeffVecBdrSTBLevP2d_sim

  public :: transp_calcMatDiagBurgP2d_sim
  public :: transp_calcMatGalBurgP2d_sim
  public :: transp_calcMatUpwBurgP2d_sim
  public :: transp_coeffMatBdrBurgP2d_sim
  public :: transp_coeffVecBdrBurgP2d_sim

contains

  !*****************************************************************************

!<subroutine>

  subroutine transp_hadaptCallback2d(iOperation, rcollection)

!<description>
    ! This callback function is used to perform postprocessing tasks
    ! such as insertion/removal of elements and or vertices in the
    ! grid adaptivity procedure in 2D.
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
    !   IquickAccess(2:5):   ivt1,...,ivt5
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer, save :: rsolution
    real(DP), dimension(:), pointer, save :: p_Dsolution


    ! What operation should be performed
    select case(iOperation)

    case(HADAPT_OPR_INITCALLBACK)
      ! Retrieve solution vector from collection
      rsolution => rcollection%p_rvectorQuickAccess1
      call lsysbl_getbase_double(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_DONECALLBACK)
      ! Nullify solution vector
      nullify(rsolution, p_Dsolution)

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


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
      p_Dsolution(rcollection%IquickAccess(1)) =&
          0.5_DP*(p_Dsolution(rcollection%IquickAccess(2))+&
                  p_Dsolution(rcollection%IquickAccess(3)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_INSERTVERTEXCENTR)
      ! Insert vertex into solution vector
      if (rsolution%NEQ .lt. rcollection%IquickAccess(1)) then
        call lsysbl_resizeVectorBlock(rsolution,&
            rcollection%IquickAccess(1), .false.)
        call lsysbl_getbase_double(rsolution, p_Dsolution)
      end if
      p_Dsolution(rcollection%IquickAccess(1)) =&
          0.25_DP*(p_Dsolution(rcollection%IquickAccess(2))+&
                   p_Dsolution(rcollection%IquickAccess(3))+&
                   p_Dsolution(rcollection%IquickAccess(4))+&
                   p_Dsolution(rcollection%IquickAccess(5)))

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case(HADAPT_OPR_REMOVEVERTEX)
      ! Remove vertex from solution
      if (rcollection%IquickAccess(2) .ne. 0) then
        p_Dsolution(rcollection%IquickAccess(1)) =&
            p_Dsolution(rcollection%IquickAccess(2))
      else
        p_Dsolution(rcollection%IquickAccess(1)) = 0.0_DP
      end if

      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)


    case default
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine transp_hadaptCallback2d

  !*****************************************************************************

!<subroutine>

  subroutine transp_refFuncBdrInt2d_sim(cderivative, rdiscretisation,&
      DpointsRef, Dpoints, ibct, DpointPar, Ielements, Dvalues,&
      rcollection)

!<description>

    ! This subroutine is called during the calculation of errors with
    ! boundary integrals. It has to compute the values of the function
    !   $$ u {\bf v}\cdot{\bf n} $$
    ! where $u$ is the exact solution and ${\bf v}$ is the exact
    ! velocity vector and ${\bf n}$ denotes the outward unit normal.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real and reference coordinates.
    ! It has to to simultaneously compute the desired values for all these points.
!</description>

!<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

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

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   DquickAccess(1): simulation time
    !   SquickAccess(1): section name in the collection
    !   SquickAccess(2): string identifying the function parser
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
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(:,:,:), pointer :: Dcoefficients
    real(DP), dimension(:,:), pointer :: Dnx,Dny
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dtime
    integer :: iel,ipoint,icompFunc,icompVelX,icompVelY

    
    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    dtime = rcollection%DquickAccess(1)

    ! This subroutine assumes that the first quick access integer
    ! value holds the number of the reference function.  Moreover,
    ! quick access interger values 3 and 4 hold the numbers of the
    ! functions to be evaluated for the x-velocity and y-velocity.
    icompFunc = rcollection%IquickAccess(1)
    icompVelX = rcollection%IquickAccess(3)
    icompVelY = rcollection%IquickAccess(4)

    ! Initialize values
    Dvalue = 0.0_DP
    Dvalue(NDIM3D+1) = dtime

    ! Allocate temporal memory
    allocate(Dcoefficients(size(Dvalues,1), size(Dvalues,2), 3))
    allocate(Dnx(size(DpointPar,1),size(DpointPar,2)))
    allocate(Dny(size(DpointPar,1),size(DpointPar,2)))

    ! Get the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
    
    ! Evaluate the reference function and the exact velocities
    do iel = 1, size(Ielements)
      do ipoint = 1, ubound(Dpoints,2)

        ! Set values for function parser
        Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)
        
        ! Evaluate function parser
        call fparser_evalFunction(p_rfparser, icompFunc,  Dvalue,&
            Dcoefficients(ipoint,iel,1))
        call fparser_evalFunction(p_rfparser, icompVelX, Dvalue,&
            Dcoefficients(ipoint,iel,2))
        call fparser_evalFunction(p_rfparser, icompVelY, Dvalue,&
            Dcoefficients(ipoint,iel,3))
        
        ! Compute the expression from the data stored in Dcoefficients
        !
        !    u * (v x n)
        !
        ! in each cubature point on each elements
        Dvalues(ipoint,iel) = Dcoefficients(ipoint,iel,1) *&
                              (Dnx(ipoint,iel) * Dcoefficients(ipoint,iel,2) +&
                               Dny(ipoint,iel) * Dcoefficients(ipoint,iel,3))
      end do
    end do

    ! Free temporal memory
    deallocate(Dcoefficients, Dnx, Dny)

  end subroutine transp_refFuncBdrInt2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_errorBdrInt2d_sim(cderivative, rdiscretisation,&
      DpointsRef, Dpoints, ibct, DpointPar, Ielements, Dvalues,&
      rcollection)

!<description>

    ! This subroutine is called during the calculation of errors with
    ! boundary integrals. It has to compute the values of the function
    !   $$ u {\bf v}\cdot{\bf n} - u_h {\bf v}_h\cdot{\bf n} $$
    ! where $u$ is the exact solution and $u_h$ is its FE approximation.
    ! Moreover, ${\bf v}$ and ${\bf v}_h$ are the exact and approximate
    ! velocity vectors and ${\bf n}$ denotes the outward unit normal.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real and reference coordinates.
    ! It has to to simultaneously compute the desired values for all these points.
!</description>

!<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

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

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   rvectorQuickAccess1: solution vector
    !   rvectorQuickAccess2: velocity field
    !   IquickAccess(1):     number of the reference function
    !   IquickAccess(3):     number of the reference x-velocity
    !   IquickAccess(4):     number of the reference y-velocity
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
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
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution, p_rvelocity
    real(DP), dimension(:,:,:), pointer :: Dcoefficients
    real(DP), dimension(:,:), pointer :: Dnx,Dny
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dtime
    integer :: iel,ipoint,icompFunc,icompVelX,icompVelY


    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access vector
    ! points to the primal solution vector and the second quick access
    ! vector points to the velocity vector
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2

    ! Evaluate the FE function in the cubature points on the boundary
    call fevl_evaluate_sim(DER_FUNC2D, Dvalues,&
        p_rsolution%RvectorBlock(1), Dpoints, Ielements, DpointsRef)

    ! Allocate temporal memory
    allocate(Dcoefficients(size(Dvalues,1), size(Dvalues,2), 5))
    allocate(Dnx(size(DpointPar,1),size(DpointPar,2)))
    allocate(Dny(size(DpointPar,1),size(DpointPar,2)))

    ! Evaluate the velocity field in the cubature points on the boundary
    ! and store the result in Dcoefficients(:,:,1:2)
    call fevl_evaluate_sim(DER_FUNC2D, Dcoefficients(:,:,1),&
        p_rvelocity%RvectorBlock(1), Dpoints, Ielements, DpointsRef)
    call fevl_evaluate_sim(DER_FUNC2D, Dcoefficients(:,:,2),&
        p_rvelocity%RvectorBlock(2), Dpoints, Ielements, DpointsRef)

    ! This subroutine assumes that the first quick access integer
    ! value holds the number of the reference function.  Moreover,
    ! quick access interger values 3 and 4 hold the numbers of the
    ! functions to be evaluated for the x-velocity and y-velocity
    icompFunc = rcollection%IquickAccess(1)
    icompVelX = rcollection%IquickAccess(3)
    icompVelY = rcollection%IquickAccess(4)

    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    dtime = rcollection%DquickAccess(1)

    ! Initialize values
    Dvalue = 0.0_DP
    Dvalue(NDIM3D+1) = dtime

    ! Get the normal vectors in the cubature points on the boundary
    call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)

    ! Evaluate the reference function and the exact velocities
    do iel = 1, size(Ielements)
      do ipoint = 1, ubound(Dpoints,2)

        ! Set values for function parser
        Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)
        
        ! Evaluate function parser
        call fparser_evalFunction(p_rfparser, icompFunc,  Dvalue,&
            Dcoefficients(ipoint,iel,3))
        call fparser_evalFunction(p_rfparser, icompVelX, Dvalue,&
            Dcoefficients(ipoint,iel,4))
        call fparser_evalFunction(p_rfparser, icompVelY, Dvalue,&
            Dcoefficients(ipoint,iel,5))
        
        ! Compute the expression from the data stored in Dcoefficients
        !
        !    u * (v x n) - u_h * (v_h x n)
        !
        ! in each cubature point on each elements
        Dvalues(ipoint,iel) = Dcoefficients(ipoint,iel,3) *&
                              (Dnx(ipoint,iel) * Dcoefficients(ipoint,iel,4) +&
                               Dny(ipoint,iel) * Dcoefficients(ipoint,iel,5))-&
                              Dvalues(ipoint,iel) *&
                              (Dnx(ipoint,iel) * Dcoefficients(ipoint,iel,1) +&
                               Dny(ipoint,iel) * Dcoefficients(ipoint,iel,2))
      end do
    end do

    ! Free temporal memory
    deallocate(Dcoefficients, Dnx, Dny)

  end subroutine transp_errorBdrInt2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_weightFuncBdrInt2d_sim(rdiscretisation, DpointsRef,&
      Dpoints, ibct, DpointPar, Ielements, Dvalues, rcollection)

!<description>
    ! This subroutine is called during the calculation of errors. It
    ! has to compute the values of a weighting function in a couple of
    ! points on a couple of elements. These values are multiplied by
    ! the calculated error.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to
    ! compute simultaneously for all these points.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

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

    ! OPTIONAL: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   IquickAccess(2): number of the weighting function
    !   SquickAccess(1): section name in the collection
    !   SquickAccess(2): string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! This array has to receive the values of the weights in all the
    ! points specified in Dpoints, or the appropriate derivative of
    ! the function, respectively, according to cderivative.
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(NDIM3D+1) :: Dvalue
    integer :: ipoint, iel, icomp


    ! Initialize values
    Dvalue = 0.0_DP

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! Moreover, this subroutine assumes that the second quick access integer
    ! value holds the number of the function to be evaluated
    icomp = rcollection%IquickAccess(2)

    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    Dvalue(NDIM3D+1) = rcollection%DquickAccess(1)

    do iel = 1, size(Ielements)
      do ipoint = 1, ubound(Dpoints,2)

        ! Set values for function parser
        Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

        ! Evaluate function parser
        call fparser_evalFunction(p_rfparser, icomp, Dvalue,&
            Dvalues(ipoint,iel))
      end do
    end do

  end subroutine transp_weightFuncBdrInt2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatDiagConvP2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IdofsAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IdofsAtNode

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
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          (p_DvelocityX(IdofsAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IdofsAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = -dscale*&
          (p_DvelocityX(IdofsAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IdofsAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagConvP2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatGalConvP2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IedgeList, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: iedge
    
    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
      DcoefficientsAtEdge(1,iedge) = dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
      DcoefficientsAtEdge(1,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
#endif
    end do

  end subroutine transp_calcMatGalConvP2d_sim
 
  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatUpwConvP2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IedgeList, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    if (dscale .gt. 0.0_DP) then

      do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
        DcoefficientsAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
        DcoefficientsAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
        DcoefficientsAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
        DcoefficientsAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient 
        !   $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
        DcoefficientsAtEdge(1,iedge) =&
            max(-DcoefficientsAtEdge(2,iedge), 0.0_DP,&
                -DcoefficientsAtEdge(3,iedge))
      end do

    else

do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = v_j*C_{ji}$
        DcoefficientsAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ij}$
        DcoefficientsAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ij}$
        DcoefficientsAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ji}$
        DcoefficientsAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient 
        !   $d_{ij} = \max\{k_{ij},0,k_{ji}\}$
        DcoefficientsAtEdge(1,iedge) =&
            max(DcoefficientsAtEdge(2,iedge), 0.0_DP,&
                DcoefficientsAtEdge(3,iedge))
      end do

    end if

  end subroutine transp_calcMatUpwConvP2d_sim
 
  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatDiagConvD2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IdofsAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IdofsAtNode

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
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: inode

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    do inode = 1, nnodes

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ii} = -v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = -dscale*&
          (p_DvelocityX(IdofsAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IdofsAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient $k_{ii} = v_i*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          (p_DvelocityX(IdofsAtNode(1,inode))*DmatrixCoeffsAtNode(1,inode)&
          +p_DvelocityY(IdofsAtNode(1,inode))*DmatrixCoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagConvD2d_sim

  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatGalConvD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IedgeList, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    do iedge = 1, nedges

#ifdef TRANSP_USE_IBP
      ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
      DcoefficientsAtEdge(1,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji] = -v_i*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
      DcoefficientsAtEdge(1,iedge) = dscale*&
          (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji] = v_i*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
#endif
    end do
    
  end subroutine transp_calcMatGalConvD2d_sim
  
  !*****************************************************************************
  
!<subroutine>

  subroutine transp_calcMatUpwConvD2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IedgeList, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the form
    ! $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

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
    real(DP), dimension(:), pointer :: p_DvelocityX,p_DvelocityY
    integer :: iedge

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity field, so lets get its double data
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(1), p_DvelocityX)
    call lsyssc_getbase_double(&
        rcollection%p_rvectorQuickAccess1%RvectorBlock(2), p_DvelocityY)

    if (dscale .gt. 0.0_DP) then

      do iedge = 1, nedges
        
#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
        DcoefficientsAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
        DcoefficientsAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
        DcoefficientsAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
        DcoefficientsAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient
        !   $d_{ij} = \max\{-k_{ij},0,-k_{ji}\}$
        DcoefficientsAtEdge(1,iedge) =&
            max(-DcoefficientsAtEdge(2,iedge), 0.0_DP,&
                -DcoefficientsAtEdge(3,iedge))
      end do

    else

      do iedge = 1, nedges
        
#ifdef TRANSP_USE_IBP
        ! Compute convective coefficient $k_{ij} = -v_j*C_{ji}$
        DcoefficientsAtEdge(2,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
        ! Compute convective coefficient $k_{ji} = -v_i*C_{ij}$
        DcoefficientsAtEdge(3,iedge) = -dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
#else
        ! Compute convective coefficient $k_{ij} = v_j*C_{ij}$
        DcoefficientsAtEdge(2,iedge) = dscale*&
            (p_DvelocityX(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(1,1,iedge)&
            +p_DvelocityY(IedgeList(2,iedge))*DmatrixCoeffsAtEdge(2,1,iedge))
        ! Compute convective coefficient $k_{ji} = v_i*C_{ji}$
        DcoefficientsAtEdge(3,iedge) = dscale*&
            (p_DvelocityX(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(1,2,iedge)&
            +p_DvelocityY(IedgeList(1,iedge))*DmatrixCoeffsAtEdge(2,2,iedge))
#endif
        
        ! Compute artificial diffusion coefficient
        !   $d_{ij} = \max\{k_{ij},0,k_{ji}\}$
        DcoefficientsAtEdge(1,iedge) =&
            max(DcoefficientsAtEdge(2,iedge), 0.0_DP,&
                DcoefficientsAtEdge(3,iedge))
      end do

    end if

  end subroutine transp_calcMatUpwConvD2d_sim

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrConvP2d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates. According
    ! to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! linear form the corresponding coefficients in front of the
    ! terms. If the code is compiled with TRANSP_USE_GFEM_AT_BOUNDARY
    ! then the boundary values are not computed directly in the
    ! cubature points. In contrast, they are computed in the degrees
    ! of freedom and their values in the cubature points it inter-
    ! polated using one-dimensional finite elements at the boundary.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation.
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

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    !   IquickAccess(3):     cubature rule
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for periodic boundary conditions
    !   DquickAccess(3):     minimim parameter of the boundary component
    !   DquickAccess(4):     maximum parameter of the boundary component
    !   DquickAccess(5):     minimim parameter of the mirror boundary component
    !   DquickAccess(6):     maximum parameter of the mirror boundary component
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
    type(t_boundaryRegion), pointer :: p_rboundaryRegionMirror
    real(DP), dimension(:), pointer :: p_Ddata, p_DvelocityX, p_DvelocityY
    real(DP), dimension(:), pointer :: Domega,DcoeffAtDOF
    real(DP), dimension(:,:), pointer :: DcubPtsRef,Dbas
    real(DP), dimension(:,:), pointer :: Dnx,Dny,DpointParMirror
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(:,:,:), pointer :: Dcoords
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dlocalData,dnv,dscale,dtime,dval
    real(DP) :: dminParam,dmaxParam,dminParamMirror,dmaxParamMirror
    integer :: ccubType,ibdrtype,icubp,iel,ipoint,isegment,ivt,npoints,nve

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrConvP2d_sim')
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

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Set pointers
    call lsysbl_getbase_double(p_rsolution, p_Ddata)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
    end if
#endif

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the cubature rule
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    ccubType = rcollection%IquickAccess(3)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Evaluate one-dimensional basis functions on the boundary edge
    if (npointsPerElement .ne. cub_igetNumPts(ccubType)) then
      call output_line('Type of cubature rule at boundary mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrConvP2d_sim')
      call sys_halt()
    else
      ! How many DOFs are located at the boundary? This should be made
      ! more flexible by checking the type of element. For the time
      ! being, only linear and bilinear finite elements are supported
      npoints = 2
      
      ! How many vertices per element do we have?
      nve = elem_igetNVE(rdomainIntSubset%celement)
      
      ! Allocate temporal memory for one-dimensional
      ! cubature along the boundary edge
      allocate(Dbas(npoints,npointsPerElement))
      allocate(Domega(npointsPerElement))
      allocate(DcubPtsRef(1,npointsPerElement))

      ! Get the coordinates of the cubature points and the
      ! corresponding weights for the given cubature rule
      call cub_getCubature(ccubType, DcubPtsRef, Domega)
      
      ! Evaluate the one-dimensional basis functions 
      ! in the cubature points on the boundary
      do icubp = 1, npointsPerElement
        Dbas(1,icubp) = 0.5_DP*(1.0_DP-DcubPtsRef(1,icubp))
        Dbas(2,icubp) = 0.5_DP*(1.0_DP+DcubPtsRef(1,icubp))
      end do

      ! Deallocate temporal memory which is no longer required
      deallocate(DcubPtsRef,Domega)
    end if    
#else
    ! Boundary values are evaluated directly at the cubature points
    npoints = npointsPerElement
#endif

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
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrConvP2d_sim')


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
      ! Dcoefficients(1,:,:).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)

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

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,1), DcoeffAtDOF(npoints))
#endif

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Loop over all elements
      do iel = 1, nelements

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,1) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,1) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
#endif

        ! Loop over all points
        do ipoint = 1, npoints

          ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,1)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif

          ! Evaluate function parser for Dirichlet value
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Impose Dirichlet value via penalty method
          DcoeffAtDOF(ipoint) = dscale * dval * BDRC_DIRICHLET_PENALTY
#else
          ! Impose Dirichlet value via penalty method
          Dcoefficients(1,ipoint,iel) = dscale * dval * BDRC_DIRICHLET_PENALTY
#endif
        end do ! ipoint
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do ! icubp
#endif
      end do ! iel

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Deallocate temporal memory
      deallocate(Dcoords, DcoeffAtDOF)
#endif


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
        allocate(Dnx(npoints,nelements), Dny(npoints,nelements),&
                 Daux(npoints,nelements,2))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Allocate temporal memory
        allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))

        ! We need the solution values and the physical coordinates of the DOFs 
        ! on the boundary, so compute them and store the result in Daux
        do iel = 1, nelements
          ! Get global DOF of first endpoints
          ipoint = rdomainIntSubset%p_IelementOrientation(iel)
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(1,iel,1) = p_DvelocityX(ivt)
          Daux(1,iel,2) = p_DvelocityY(ivt)
          
          ! Store vertex coordinate
          Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

          ! Get global DOF of second endpoints
          ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(2,iel,1) = p_DvelocityX(ivt)
          Daux(2,iel,2) = p_DvelocityY(ivt)

          ! Store vertex coordinate
          Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        end do

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
        ! Evaluate the velocity field (if any) in the cubature points on
        ! the boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Calculate the normal vectors in cubature on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
        
        ! Initialize values
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime
        
        ! Evaluate the function parser for the boundary values in the
        ! cubature points or the DOFs on the boundary and store the
        ! result in Dcoefficients(:,:,3). Multiply the velocity vector
        ! with the normal in each point to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npoints
            
            ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
            Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
            
            ! Compute the normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1) +&
                  Dny(ipoint,iel) * Daux(ipoint,iel,2)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            ! Store the Robin boundary condition at the DOF on the boundary
            DcoeffAtDOF(ipoint) = -dscale * dnv * dval
#else
            ! Store the Robin boundary condition at the cubature point on the boundary
            Dcoefficients(1,ipoint,iel) = -dscale * dnv * dval
#endif
          end do ! ipoint

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Loop over the cubature points and interpolate the Robin
          ! boundary conditions from the DOFs to the cubature points,
          ! where they are needed by the linear form assembly routine
          do icubp = 1, npointsPerElement

            dlocalData = 0.0_DP

            ! Loop over the DOFs and interpolate the Robin boundary conditions
            do ipoint = 1, npoints
              dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
            end do
            
            ! Store Robin boundary condition in the cubature points
            Dcoefficients(1,icubp,iel) = dlocalData
          end do ! icubp
#endif
        end do ! iel
        
        ! Deallocate temporal memory
        deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        deallocate(Dcoords, DcoeffAtDOF)
#endif

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
        allocate(Dnx(npoints,nelements),Dny(npoints,nelements),&
                 Daux(npoints,nelements,2))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Allocate temporal memory
        allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
        
        ! We need the solution values and the physical coordinates of the DOFs 
        ! on the boundary, so compute them and store the result in Daux
        do iel = 1, nelements
          ! Get global DOF of first endpoints
          ipoint = rdomainIntSubset%p_IelementOrientation(iel)
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(1,iel,1) = p_DvelocityX(ivt)
          Daux(1,iel,2) = p_DvelocityY(ivt)
          
          ! Store vertex coordinate
          Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

          ! Get global DOF of second endpoints
          ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(2,iel,1) = p_DvelocityX(ivt)
          Daux(2,iel,2) = p_DvelocityY(ivt)

          ! Store vertex coordinate
          Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        end do

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
        ! Evaluate the velocity field (if any) in the cubature points on
        ! the boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

        ! Get the normal vectors in the cubature points on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif     
   
        ! Initialize values
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime

        ! Evaluate the function parser for the boundary values in the
        ! cubature points or the DOFS on the boundary and store the
        ! result in Dcoefficients(:,:,3). Multiply the velocity vector
        ! with the normal in each point to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npoints
                        
            ! Compute the normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1) +&
                  Dny(ipoint,iel) * Daux(ipoint,iel,2)
            
            ! Check if we are at the primal inflow boundary
            if (dnv .lt. 0.0_DP) then
              
              ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
              Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
              
              ! Evaluate function parser
              call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = -dscale * dnv * dval
#else
              Dcoefficients(1,ipoint,iel) = -dscale * dnv * dval
#endif
            else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = 0.0_DP
#else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
            end if
          end do ! ipoint
          
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Loop over the cubature points and interpolate the Robin
          ! boundary conditions from the DOFs to the cubature points,
          ! where they are needed by the linear form assembly routine
          do icubp = 1, npointsPerElement

            dlocalData = 0.0_DP

            ! Loop over the DOFs and interpolate the Robin boundary conditions
            do ipoint = 1, npoints
              dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
            end do
            
            ! Store Robin boundary condition in the cubature points
            Dcoefficients(1,icubp,iel) = dlocalData
          end do ! icubp
#endif
        end do ! iel
        
        ! Deallocate temporal memory
        deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        deallocate(Dcoords, DcoeffAtDOF)
#endif

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
        
        ! Get mirrored boundary region from collection structure
        p_rboundaryRegionMirror => collct_getvalue_bdreg(rcollection,&
            'rboundaryRegionMirror', ssectionName=trim(rcollection%SquickAccess(1)))
        
        ! Get minimum/maximum parameter values from collection structure
        dminParam = rcollection%DquickAccess(3)
        dmaxParam = rcollection%DquickAccess(4)
        dminParamMirror = rcollection%DquickAccess(5)
        dmaxParamMirror = rcollection%DquickAccess(6)
        
        ! Allocate temporal memory
        allocate(Dnx(npoints,nelements), Dny(npoints,nelements))
        allocate(DpointParMirror(npoints,nelements), Daux(npoints,nelements,3))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Allocate temporal memory
        allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))

        ! We need the solution values and the physical coordinates of the DOFs 
        ! on the boundary, so compute them and store the result in Daux
        do iel = 1, nelements
          ! Get global DOF of first endpoints
          ipoint = rdomainIntSubset%p_IelementOrientation(iel)
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(1,iel,1) = p_DvelocityX(ivt)
          Daux(1,iel,2) = p_DvelocityY(ivt)
          
          ! Store vertex coordinate
          Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

          ! Get global DOF of second endpoints
          ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(2,iel,1) = p_DvelocityX(ivt)
          Daux(2,iel,2) = p_DvelocityY(ivt)

          ! Store vertex coordinate
          Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        end do

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else  
        ! Evaluate the velocity field (if any) in the cubature points on
        ! the boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Get the normal vectors in the cubature points on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
        
        if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Rescale parameter values DedgePosition on the boundary
          ! segment where to compute the boundary conditions into
          ! parameter values on the mirror boundary region
          call mprim_linearRescale(rdomainIntSubset%p_DedgePosition, dminParam,&
              dmaxParam, dmaxParamMirror, dminParamMirror, DpointParMirror)
#else
          ! Rescale parameter values DpointPar on the boundary segment
          ! where to compute the boundary conditions into parameter
          ! values on the mirror boundary region
          call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
              dmaxParamMirror, dminParamMirror, DpointParMirror)
#endif
        else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Rescale parameter values DedgePosition on the boundary
          ! segment where to compute the boundary conditions into
          ! parameter values on the reversely mirror boundary region
          call mprim_linearRescale(rdomainIntSubset%p_DedgePosition, dminParam,&
              dmaxParam, dmaxParamMirror, dminParamMirror, DpointParMirror)
#else
          ! Rescale parameter values DpointPar on the boundary segment
          ! where to compute the boundary conditions into parameter
          ! values on the reversely mirror boundary region
          call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
              dmaxParamMirror, dminParamMirror, DpointParMirror)
#endif
        end if
        
        ! Evaluate the solution in the cubature points on the mirrored
        ! boundary and store the result in Daux(:,:,3)
        call doEvaluateAtBdr2d(DER_FUNC, npoints*nelements,&
            Daux(:,:,3), p_rsolution%RvectorBlock(1), DpointParMirror,&
            ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)

        do iel = 1, nelements
          do ipoint = 1, npoints
                        
            ! Compute the normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1) +&
                  Dny(ipoint,iel) * Daux(ipoint,iel,2)

            ! Check if we are at the primal inflow boundary
            if (dnv .lt. 0.0_DP) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = -dscale * dnv * Daux(ipoint,iel,3)
#else
              Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,3)
#endif
            else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = 0.0_DP
#else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
            end if
          end do

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
        end do
        
        ! Deallocate temporal memory
        deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Deallocate temporal memory
      deallocate(Dcoords, DcoeffAtDOF)
#endif

      else

        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP

      end if

      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrConvP2d_sim')
      call sys_halt()
      
    end select

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Deallocate temporal memory
    deallocate(Dbas)
#endif

  contains

    ! Here come the working routines

    !***************************************************************************
    ! Evaluate th solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    ! ***************************************************************************
    
    subroutine doEvaluateAtBdr2d(iderType, n, Dvalues, rvectorScalar,&
        DpointsPar, ibdc, cparType, rboundaryRegion)
      
      integer, intent(in) :: iderType,ibdc,cparType,n
      real(DP), dimension(n), intent(in) :: DpointsPar
      type(t_vectorScalar), intent(in) :: rvectorScalar
      type(t_boundaryRegion), intent(in) :: rboundaryRegion
      
      real(DP), dimension(n), intent(out) :: Dvalues
      
      call fevl_evaluateBdr2D(iderType, Dvalues, rvectorScalar,&
          DpointsPar, ibdc, cparType, rboundaryRegion)
      
    end subroutine doEvaluateAtBdr2d
    
  end subroutine transp_coeffVecBdrConvP2d_sim

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrConvD2d_sim(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

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
    ! terms. If the code is compiled with TRANSP_USE_GFEM_AT_BOUNDARY
    ! then the boundary values are not computed directly in the
    ! cubature points. In contrast, they are computed in the degrees
    ! of freedom and their values in the cubature points it inter-
    ! polated using one-dimensional finite elements at the boundary.
    !
    ! This routine handles the dual problem for the
    ! convection-diffusion equation.
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

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    !   IquickAccess(3):     cubature rule
    !   SquickAccess(1):     section name in the collection
    !   SquickAccess(2):     string identifying the function parser
    !
    ! only for periodic boundary conditions
    !   DquickAccess(3):     minimim parameter of the boundary component
    !   DquickAccess(4):     maximum parameter of the boundary component
    !   DquickAccess(5):     minimim parameter of the mirror boundary component
    !   DquickAccess(6):     maximum parameter of the mirror boundary component
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
    type(t_vectorBlock), pointer :: p_rsolution, p_rvelocity
    type(t_boundaryRegion), pointer :: p_rboundaryRegionMirror
    real(DP), dimension(:), pointer :: p_Ddata, p_DvelocityX, p_DvelocityY
    real(DP), dimension(:), pointer :: Domega,DcoeffAtDOF
    real(DP), dimension(:,:), pointer :: DcubPtsRef,Dbas
    real(DP), dimension(:,:), pointer :: Dnx,Dny,DpointParMirror
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(:,:,:), pointer :: Dcoords
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dlocalData,dnv,dscale,dtime,dval
    real(DP) :: dminParam,dmaxParam,dminParamMirror,dmaxParamMirror
    integer :: ccubType,ibdrtype,icubp,iel,ipoint,isegment,ivt,npoints,nve
    
#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrConvD2d_sim')
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

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Set pointers
    call lsysbl_getbase_double(p_rsolution, p_Ddata)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
    end if
#endif

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the cubature rule
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    ccubType = rcollection%IquickAccess(3)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Evaluate one-dimensional basis functions on the boundary edge
    if (npointsPerElement .ne. cub_igetNumPts(ccubType)) then
      call output_line('Type of cubature rule at boundary mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrConvD2d_sim')
      call sys_halt()
    else
      ! How many DOFs are located at the boundary? This should be made
      ! more flexible by checking the type of element. For the time
      ! being, only linear and bilinear finite elements are supported
      npoints = 2
      
      ! How many vertices per element do we have?
      nve = elem_igetNVE(rdomainIntSubset%celement)
      
      ! Allocate temporal memory for one-dimensional
      ! cubature along the boundary edge
      allocate(Dbas(npoints,npointsPerElement))
      allocate(Domega(npointsPerElement))
      allocate(DcubPtsRef(1,npointsPerElement))

      ! Get the coordinates of the cubature points and the
      ! corresponding weights for the given cubature rule
      call cub_getCubature(ccubType, DcubPtsRef, Domega)
      
      ! Evaluate the one-dimensional basis functions 
      ! in the cubature points on the boundary
      do icubp = 1, npointsPerElement
        Dbas(1,icubp) = 0.5_DP*(1.0_DP-DcubPtsRef(1,icubp))
        Dbas(2,icubp) = 0.5_DP*(1.0_DP+DcubPtsRef(1,icubp))
      end do

      ! Deallocate temporal memory which is no longer required
      deallocate(DcubPtsRef,Domega)
    end if    
#else
    ! Boundary values are evaluated directly at the cubature points
    npoints = npointsPerElement
#endif

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
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrConvD2d_sim')


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
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)
          
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

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,1), DcoeffAtDOF(npoints))
#endif

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      do iel = 1, nelements

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,1) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,1) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
#endif

        do ipoint = 1, npointsPerElement

          ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,1)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif

          ! Evaluate function parser for Dirichlet value
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Impose Dirichlet value via penalty method
          DcoeffAtDOF(ipoint) = -dscale * dval * BDRC_DIRICHLET_PENALTY
#else
          ! Impose Dirichlet value via penalty method
          Dcoefficients(1,ipoint,iel) = -dscale * dval * BDRC_DIRICHLET_PENALTY
#endif
        end do

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Deallocate temporal memory
      deallocate(Dcoords, DcoeffAtDOF)
#endif

      
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
        allocate(Dnx(npoints,nelements), Dny(npoints,nelements),&
                 Daux(npoints,nelements,2))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Allocate temporal memory
        allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))

        ! We need the solution values and the physical coordinates of the DOFs 
        ! on the boundary, so compute them and store the result in Daux
        do iel = 1, nelements
          ! Get global DOF of first endpoints
          ipoint = rdomainIntSubset%p_IelementOrientation(iel)
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(1,iel,1) = p_DvelocityX(ivt)
          Daux(1,iel,2) = p_DvelocityY(ivt)
          
          ! Store vertex coordinate
          Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

          ! Get global DOF of second endpoints
          ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(2,iel,1) = p_DvelocityX(ivt)
          Daux(2,iel,2) = p_DvelocityY(ivt)

          ! Store vertex coordinate
          Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        end do

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
        ! Evaluate the velocity field (if any) in the cubature points on
        ! the boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Calculate the normal vectors in cubature on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
        
        ! Initialize values
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime
        
        ! Evaluate the function parser for the boundary values in the
        ! cubature points or the DOFs on the boundary and store the
        ! result in Dcoefficients(:,:,3). Multiply the velocity vector
        ! with the normal in each point to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npoints
            
            ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
            Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
            
            ! Compute the normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1) +&
                  Dny(ipoint,iel) * Daux(ipoint,iel,2)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            ! Store the Robin boundary condition at the DOF on the boundary
            DcoeffAtDOF(ipoint) = dscale * dnv * dval
#else
            ! Store the Robin boundary condition at the cubature point on the boundary
            Dcoefficients(1,ipoint,iel) = dscale * dnv * dval
#endif
          end do ! ipoint

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Loop over the cubature points and interpolate the Robin
          ! boundary conditions from the DOFs to the cubature points,
          ! where they are needed by the linear form assembly routine
          do icubp = 1, npointsPerElement

            dlocalData = 0.0_DP

            ! Loop over the DOFs and interpolate the Robin boundary conditions
            do ipoint = 1, npoints
              dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
            end do
            
            ! Store Robin boundary condition in the cubature points
            Dcoefficients(1,icubp,iel) = dlocalData
          end do ! icubp
#endif
        end do ! iel
        
        ! Deallocate temporal memory
        deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        deallocate(Dcoords, DcoeffAtDOF)
#endif

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
        allocate(Dnx(npoints,nelements),Dny(npoints,nelements),&
                 Daux(npoints,nelements,2))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Allocate temporal memory
        allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
        
        ! We need the solution values and the physical coordinates of the DOFs 
        ! on the boundary, so compute them and store the result in Daux
        do iel = 1, nelements
          ! Get global DOF of first endpoints
          ipoint = rdomainIntSubset%p_IelementOrientation(iel)
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(1,iel,1) = p_DvelocityX(ivt)
          Daux(1,iel,2) = p_DvelocityY(ivt)
          
          ! Store vertex coordinate
          Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

          ! Get global DOF of second endpoints
          ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(2,iel,1) = p_DvelocityX(ivt)
          Daux(2,iel,2) = p_DvelocityY(ivt)

          ! Store vertex coordinate
          Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        end do

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
        ! Evaluate the velocity field (if any) in the cubature points on
        ! the boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

        ! Get the normal vectors in the cubature points on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif     
   
        ! Initialize values
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime

        ! Evaluate the function parser for the boundary values in the
        ! cubature points or the DOFS on the boundary and store the
        ! result in Dcoefficients(:,:,3). Multiply the velocity vector
        ! with the normal in each point to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npoints
            
            ! Compute the normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1) +&
                  Dny(ipoint,iel) * Daux(ipoint,iel,2)
            
            ! Check if we are at the dual inflow boundary
            if (dnv .gt. 0.0_DP) then
              
              ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
              Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
              
              ! Evaluate function parser
              call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
              
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = dscale * dnv * dval
#else
              Dcoefficients(1,ipoint,iel) = dscale * dnv * dval
#endif
            else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = 0.0_DP
#else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
            end if
          end do ! ipoint

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Loop over the cubature points and interpolate the Robin
          ! boundary conditions from the DOFs to the cubature points,
          ! where they are needed by the linear form assembly routine
          do icubp = 1, npointsPerElement

            dlocalData = 0.0_DP

            ! Loop over the DOFs and interpolate the Robin boundary conditions
            do ipoint = 1, npoints
              dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
            end do
            
            ! Store Robin boundary condition in the cubature points
            Dcoefficients(1,icubp,iel) = dlocalData
          end do ! icubp
#endif
        end do ! iel
        
        ! Deallocate temporal memory
        deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        deallocate(Dcoords, DcoeffAtDOF)
#endif

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
        
        ! Get mirrored boundary region from collection structure
        p_rboundaryRegionMirror => collct_getvalue_bdreg(rcollection,&
            'rboundaryRegionMirror', ssectionName=trim(rcollection%SquickAccess(1)))
        
        ! Get minimum/maximum parameter values from collection structure
        dminParam = rcollection%DquickAccess(3)
        dmaxParam = rcollection%DquickAccess(4)
        dminParamMirror = rcollection%DquickAccess(5)
        dmaxParamMirror = rcollection%DquickAccess(6)
        
        ! Allocate temporal memory
        allocate(Dnx(npoints,nelements), Dny(npoints,nelements))
        allocate(DpointParMirror(npoints,nelements), Daux(npoints,nelements,3))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Allocate temporal memory
        allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))

        ! We need the solution values and the physical coordinates of the DOFs 
        ! on the boundary, so compute them and store the result in Daux
        do iel = 1, nelements
          ! Get global DOF of first endpoints
          ipoint = rdomainIntSubset%p_IelementOrientation(iel)
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(1,iel,1) = p_DvelocityX(ivt)
          Daux(1,iel,2) = p_DvelocityY(ivt)
          
          ! Store vertex coordinate
          Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

          ! Get global DOF of second endpoints
          ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(2,iel,1) = p_DvelocityX(ivt)
          Daux(2,iel,2) = p_DvelocityY(ivt)

          ! Store vertex coordinate
          Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        end do

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
        ! Evaluate the velocity field (if any) in the cubature points on
        ! the boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints, &
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Get the normal vectors in the cubature points on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
        
        if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Rescale parameter values DedgePosition on the boundary
          ! segment where to compute the boundary conditions into
          ! parameter values on the mirror boundary region
          call mprim_linearRescale(rdomainIntSubset%p_DedgePosition, dminParam,&
              dmaxParam, dmaxParamMirror, dminParamMirror, DpointParMirror)
#else
          ! Rescale parameter values DpointPar on the boundary segment
          ! where to compute the boundary conditions into parameter
          ! values on the mirror boundary region
          call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
              dmaxParamMirror, dminParamMirror, DpointParMirror)
#endif
        else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Rescale parameter values DedgePosition on the boundary
          ! segment where to compute the boundary conditions into
          ! parameter values on the reversely mirror boundary region
          call mprim_linearRescale(rdomainIntSubset%p_DedgePosition, dminParam,&
              dmaxParam, dmaxParamMirror, dminParamMirror, DpointParMirror)
#else
          ! Rescale parameter values DpointPar on the boundary segment
          ! where to compute the boundary conditions into parameter
          ! values on the reversely mirror boundary region
          call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
              dmaxParamMirror, dminParamMirror, DpointParMirror)
#endif
        end if
        
        ! Evaluate the solution in the cubature points on the mirrored
        ! boundary and store the result in Daux(:,:,3)
        call doEvaluateAtBdr2d(DER_FUNC, npointsPerElement*nelements,&
            Daux(:,:,3), p_rsolution%RvectorBlock(1), DpointParMirror,&
            ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)

        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
                        
            ! Compute the normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1) +&
                  Dny(ipoint,iel) * Daux(ipoint,iel,2)

            ! Check if we are at the dual inflow boundary
            if (dnv .gt. 0.0_DP) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = dscale * dnv * Daux(ipoint,iel,3)
#else
              Dcoefficients(1,ipoint,iel) = dscale * dnv * Daux(ipoint,iel,3)
#endif
            else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = 0.0_DP
#else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
            end if
          end do ! ipoint

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do ! icubp
#endif
      end do ! iel
      
      ! Deallocate temporal memory
      deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Deallocate temporal memory
      deallocate(Dcoords, DcoeffAtDOF)
#endif

      else

        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP

      end if


    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrConvD2d_sim')
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
    
    subroutine doEvaluateAtBdr2d(iderType, n, Dvalues, rvectorScalar,&
        DpointsPar, ibdc, cparType, rboundaryRegion)
      
      integer, intent(in) :: iderType,ibdc,cparType,n
      real(DP), dimension(n), intent(in) :: DpointsPar
      type(t_vectorScalar), intent(in) :: rvectorScalar
      type(t_boundaryRegion), intent(in) :: rboundaryRegion
      
      real(DP), dimension(n), intent(out) :: Dvalues
      
      call fevl_evaluateBdr2D(iderType, Dvalues, rvectorScalar,&
          DpointsPar, ibdc, cparType, rboundaryRegion)
      
    end subroutine doEvaluateAtBdr2d
    
  end subroutine transp_coeffVecBdrConvD2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrConvP2d_sim(rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to
    ! compute the coefficients in front of the terms of the bilinear
    ! form.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates.  According
    ! to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! bilinear form the corresponding coefficients in front of the
    ! terms. If the code is compiled with TRANSP_USE_GFEM_AT_BOUNDARY
    ! then the boundary values are not computed directly in the
    ! cubature points. In contrast, they are computed in the degrees
    ! of freedom and their values in the cubature points it inter-
    ! polated using one-dimensional finite elements at the boundary.
    !
    ! This routine handles the primal problem for the
    ! convection-diffusion equation.
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

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    !   IquickAccess(3):     cubature rule
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
    type(t_vectorBlock), pointer :: p_rsolution, p_rvelocity
    real(DP), dimension(:), pointer :: p_Ddata, p_DvelocityX, p_DvelocityY
    real(DP), dimension(:), pointer :: Domega,DcoeffAtDOF
    real(DP), dimension(:,:), pointer :: DcubPtsRef,Dbas
    real(DP), dimension(:,:), pointer :: Dnx,Dny
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(:,:,:), pointer :: Dcoords
    real(DP) :: dlocalData,dnv,dscale,dtime
    integer :: ccubType,ibdrtype,icubp,iel,ipoint,isegment,ivt,npoints,nve

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrConvP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! point to the solution and velocity vector (if any)
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Set pointers
    call lsysbl_getbase_double(p_rsolution, p_Ddata)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
    end if
#endif

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the cubature rule
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    ccubType = rcollection%IquickAccess(3)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Evaluate one-dimensional basis functions on the boundary edge
    if (npointsPerElement .ne. cub_igetNumPts(ccubType)) then
      call output_line('Type of cubature rule at boundary mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrConvP2d_sim')
      call sys_halt()
    else
      ! How many DOFs are located at the boundary? This should be made
      ! more flexible by checking the type of element. For the time
      ! being, only linear and bilinear finite elements are supported
      npoints = 2
      
      ! How many vertices per element do we have?
      nve = elem_igetNVE(rdomainIntSubset%celement)
      
      ! Allocate temporal memory for one-dimensional
      ! cubature along the boundary edge
      allocate(Dbas(npoints,npointsPerElement))
      allocate(Domega(npointsPerElement))
      allocate(DcubPtsRef(1,npointsPerElement))

      ! Get the coordinates of the cubature points and the
      ! corresponding weights for the given cubature rule
      call cub_getCubature(ccubType, DcubPtsRef, Domega)
      
      ! Evaluate the one-dimensional basis functions 
      ! in the cubature points on the boundary
      do icubp = 1, npointsPerElement
        Dbas(1,icubp) = 0.5_DP*(1.0_DP-DcubPtsRef(1,icubp))
        Dbas(2,icubp) = 0.5_DP*(1.0_DP+DcubPtsRef(1,icubp))
      end do

      ! Deallocate temporal memory which is no longer required
      deallocate(DcubPtsRef,Domega)
    end if    
#else
    ! Boundary values are evaluated directly at the cubature points
    npoints = npointsPerElement
#endif

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral (if any)

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory
        allocate(Dnx(npoints,nelements), Dny(npoints,nelements),&
                 Daux(npoints,nelements,2))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Allocate temporal memory
        allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))

        ! We need the solution values and the physical coordinates of the DOFs 
        ! on the boundary, so compute them and store the result in Daux
        do iel = 1, nelements
          ! Get global DOF of first endpoints
          ipoint = rdomainIntSubset%p_IelementOrientation(iel)
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(1,iel,1) = p_DvelocityX(ivt)
          Daux(1,iel,2) = p_DvelocityY(ivt)
          
          ! Store vertex coordinate
          Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

          ! Get global DOF of second endpoints
          ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(2,iel,1) = p_DvelocityX(ivt)
          Daux(2,iel,2) = p_DvelocityY(ivt)

          ! Store vertex coordinate
          Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        end do

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
        ! Evaluate the velocity field (if any) in the cubature points
        ! on the boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Get the normal vectors in the cubature points on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npoints
            
            ! Compute the normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1) +&
                  Dny(ipoint,iel) * Daux(ipoint,iel,2)
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY      
            ! Scale normal velocity by scaling parameter
            DcoeffAtDOF(ipoint) = -dscale * dnv
#else
            ! Scale normal velocity by scaling parameter
            Dcoefficients(1,ipoint,iel) = -dscale * dnv
#endif
          end do ! ipoint

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Loop over the cubature points and interpolate the Robin
          ! boundary conditions from the DOFs to the cubature points,
          ! where they are needed by the linear form assembly routine
          do icubp = 1, npointsPerElement
            
            dlocalData = 0.0_DP
            
            ! Loop over the DOFs and interpolate the Robin boundary conditions
            do ipoint = 1, npoints
              dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
            end do
            
            ! Store Robin boundary condition in the cubature points
            Dcoefficients(1,icubp,iel) = dlocalData
          end do ! icubp
#endif
        end do ! iel
        
        ! Free temporal memory
        deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        deallocate(Dcoords, DcoeffAtDOF)
#endif
        
      else
        
        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP

      end if


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      ! Impose penalty parameter

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
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrConvP2d_sim')

      
    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      !
      ! The convective part of the boundary integral at the outflow is
      ! likewise assembled for periodic and antiperiodic boundary conditions

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory
        allocate(Dnx(npoints,nelements), Dny(npoints,nelements),&
                 Daux(npoints,nelements,2))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Allocate temporal memory
        allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))

        ! We need the solution values and the physical coordinates of the DOFs 
        ! on the boundary, so compute them and store the result in Daux
        do iel = 1, nelements
          ! Get global DOF of first endpoints
          ipoint = rdomainIntSubset%p_IelementOrientation(iel)
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(1,iel,1) = p_DvelocityX(ivt)
          Daux(1,iel,2) = p_DvelocityY(ivt)
          
          ! Store vertex coordinate
          Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

          ! Get global DOF of second endpoints
          ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(2,iel,1) = p_DvelocityX(ivt)
          Daux(2,iel,2) = p_DvelocityY(ivt)

          ! Store vertex coordinate
          Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        end do

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
        ! Evaluate the velocity field (if any) in the cubature points
        ! on the boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

        ! Get the normal vectors in the cubature points on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npoints
            
            ! Compute the normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1) +&
                  Dny(ipoint,iel) * Daux(ipoint,iel,2)

            ! Check if we are at the primal outflow boundary
            if (dnv .gt. 0.0_DP) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = -dscale * dnv
#else
              Dcoefficients(1,ipoint,iel) = -dscale * dnv
#endif
            else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = 0.0_DP
#else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
            end if
          end do ! ipoint

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Loop over the cubature points and interpolate the Robin
          ! boundary conditions from the DOFs to the cubature points,
          ! where they are needed by the linear form assembly routine
          do icubp = 1, npointsPerElement
            
            dlocalData = 0.0_DP
            
            ! Loop over the DOFs and interpolate the Robin boundary conditions
            do ipoint = 1, npoints
              dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
            end do
            
            ! Store Robin boundary condition in the cubature points
            Dcoefficients(1,icubp,iel) = dlocalData
          end do ! icubp
#endif
        end do ! iel
        
        ! Free temporal memory
        deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        deallocate(Dcoords, DcoeffAtDOF)
#endif
        
      else

        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP

      end if
      
    
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrConvP2d_sim')
      call sys_halt()
      
    end select
    
  end subroutine transp_coeffMatBdrConvP2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrConvD2d_sim(rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the matrix assembly. It has to
    ! compute the coefficients in front of the terms of the bilinear
    ! form.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in real coordinates.  According
    ! to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the
    ! bilinear form the corresponding coefficients in front of the
    ! terms. If the code is compiled with TRANSP_USE_GFEM_AT_BOUNDARY
    ! then the boundary values are not computed directly in the
    ! cubature points. In contrast, they are computed in the degrees
    ! of freedom and their values in the cubature points it inter-
    ! polated using one-dimensional finite elements at the boundary.
    !
    ! This routine handles the dual problem for the
    ! convection-diffusion equation.
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

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    !   IquickAccess(3):     cubature rule
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
    type(t_vectorBlock), pointer :: p_rsolution, p_rvelocity
    real(DP), dimension(:), pointer :: p_Ddata, p_DvelocityX, p_DvelocityY
    real(DP), dimension(:), pointer :: Domega,DcoeffAtDOF
    real(DP), dimension(:,:), pointer :: DcubPtsRef,Dbas
    real(DP), dimension(:,:), pointer :: Dnx,Dny
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(:,:,:), pointer :: Dcoords
    real(DP) :: dlocalData,dnv,dscale,dtime
    integer :: ccubType,ibdrtype,icubp,iel,ipoint,isegment,ivt,npoints,nve

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrConvD2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first two quick access vectors
    ! point to the solution and velocity vector (if any)
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Set pointers
    call lsysbl_getbase_double(p_rsolution, p_Ddata)
    if (associated(p_rvelocity)) then
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(p_rvelocity%RvectorBlock(2), p_DvelocityY)
    end if
#endif

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the cubature rule
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    ccubType = rcollection%IquickAccess(3)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Evaluate one-dimensional basis functions on the boundary edge
    if (npointsPerElement .ne. cub_igetNumPts(ccubType)) then
      call output_line('Type of cubature rule at boundary mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrConvD2d_sim')
      call sys_halt()
    else
      ! How many DOFs are located at the boundary? This should be made
      ! more flexible by checking the type of element. For the time
      ! being, only linear and bilinear finite elements are supported
      npoints = 2
      
      ! How many vertices per element do we have?
      nve = elem_igetNVE(rdomainIntSubset%celement)
      
      ! Allocate temporal memory for one-dimensional
      ! cubature along the boundary edge
      allocate(Dbas(npoints,npointsPerElement))
      allocate(Domega(npointsPerElement))
      allocate(DcubPtsRef(1,npointsPerElement))

      ! Get the coordinates of the cubature points and the
      ! corresponding weights for the given cubature rule
      call cub_getCubature(ccubType, DcubPtsRef, Domega)
      
      ! Evaluate the one-dimensional basis functions 
      ! in the cubature points on the boundary
      do icubp = 1, npointsPerElement
        Dbas(1,icubp) = 0.5_DP*(1.0_DP-DcubPtsRef(1,icubp))
        Dbas(2,icubp) = 0.5_DP*(1.0_DP+DcubPtsRef(1,icubp))
      end do

      ! Deallocate temporal memory which is no longer required
      deallocate(DcubPtsRef,Domega)
    end if    
#else
    ! Boundary values are evaluated directly at the cubature points
    npoints = npointsPerElement
#endif

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the boundary integral for the convective term (if any)

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory
        allocate(Dnx(npoints,nelements), Dny(npoints,nelements),&
                 Daux(npoints,nelements,2))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Allocate temporal memory
        allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))

        ! We need the solution values and the physical coordinates of the DOFs 
        ! on the boundary, so compute them and store the result in Daux
        do iel = 1, nelements
          ! Get global DOF of first endpoints
          ipoint = rdomainIntSubset%p_IelementOrientation(iel)
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(1,iel,1) = p_DvelocityX(ivt)
          Daux(1,iel,2) = p_DvelocityY(ivt)
          
          ! Store vertex coordinate
          Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

          ! Get global DOF of second endpoints
          ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(2,iel,1) = p_DvelocityX(ivt)
          Daux(2,iel,2) = p_DvelocityY(ivt)

          ! Store vertex coordinate
          Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        end do

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
        ! Evaluate the velocity field (if any) in the cubature points
        ! on the boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        ! Get the normal vectors in the cubature points on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npoints
            
            ! Compute the normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1) +&
                  Dny(ipoint,iel) * Daux(ipoint,iel,2)
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY      
            ! Scale normal velocity by scaling parameter
            DcoeffAtDOF(ipoint) = dscale * dnv
#else
            ! Scale normal velocity by scaling parameter
            Dcoefficients(1,ipoint,iel) = dscale * dnv
#endif
          end do ! ipoint

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Loop over the cubature points and interpolate the Robin
          ! boundary conditions from the DOFs to the cubature points,
          ! where they are needed by the linear form assembly routine
          do icubp = 1, npointsPerElement
            
            dlocalData = 0.0_DP
            
            ! Loop over the DOFs and interpolate the Robin boundary conditions
            do ipoint = 1, npoints
              dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
            end do
            
            ! Store Robin boundary condition in the cubature points
            Dcoefficients(1,icubp,iel) = dlocalData
          end do ! icubp
#endif
        end do ! iel
        
        ! Free temporal memory
        deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        deallocate(Dcoords, DcoeffAtDOF)
#endif

      else

        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP

      end if


    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      ! Impose penalty parameter

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Impose Dirichlet boundary conditions via penalty method
          Dcoefficients(1,ipoint,iel) = dscale * BDRC_DIRICHLET_PENALTY
        end do
      end do


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      ! Do nothing since the boundary values are build into the linear form.

      Dcoefficients = 0.0_DP

      ! This routine should not be called at all for homogeneous Neumann boundary
      ! conditions since it corresponds to an expensive assembly of "zero".
      call output_line('Redundant assembly of vanishing boundary term!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrConvD2d_sim')
      

    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      !
      ! The convective part of the boundary integral at the outflow is
      ! likewise assembled for periodic and antiperiodic boundary conditions

      if (associated(p_rvelocity)) then

        ! Allocate temporal memory
        allocate(Dnx(npoints,nelements), Dny(npoints,nelements),&
                 Daux(npoints,nelements,2))

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Allocate temporal memory
        allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))

        ! We need the solution values and the physical coordinates of the DOFs 
        ! on the boundary, so compute them and store the result in Daux
        do iel = 1, nelements
          ! Get global DOF of first endpoints
          ipoint = rdomainIntSubset%p_IelementOrientation(iel)
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(1,iel,1) = p_DvelocityX(ivt)
          Daux(1,iel,2) = p_DvelocityY(ivt)
          
          ! Store vertex coordinate
          Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)

          ! Get global DOF of second endpoints
          ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
          ivt    = IdofsTest(ipoint,iel)

          ! Store velocity data
          Daux(2,iel,1) = p_DvelocityX(ivt)
          Daux(2,iel,2) = p_DvelocityY(ivt)

          ! Store vertex coordinate
          Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        end do

        ! Calculate the normal vectors in DOFs on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
        ! Evaluate the velocity field (if any) in the cubature points
        ! on the boundary and store the result in Daux(:,:,1:2)
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
            p_rvelocity%RvectorBlock(1), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
        
        call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,2),&
            p_rvelocity%RvectorBlock(2), Dpoints,&
            rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

        ! Get the normal vectors in the cubature points on the boundary
        call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
        
        ! Multiply the velocity vector with the normal in each point
        ! to get the normal velocity.
        do iel = 1, nelements
          do ipoint = 1, npoints
            
            ! Compute the normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1) +&
                  Dny(ipoint,iel) * Daux(ipoint,iel,2)

            ! Check if we are at the dual outflow boundary
            if (dnv .lt. 0.0_DP) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = dscale * dnv
#else
              Dcoefficients(1,ipoint,iel) = dscale * dnv
#endif
            else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
              DcoeffAtDOF(ipoint) = 0.0_DP
#else
              Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
            end if
          end do ! ipoint

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Loop over the cubature points and interpolate the Robin
          ! boundary conditions from the DOFs to the cubature points,
          ! where they are needed by the linear form assembly routine
          do icubp = 1, npointsPerElement
            
            dlocalData = 0.0_DP
            
            ! Loop over the DOFs and interpolate the Robin boundary conditions
            do ipoint = 1, npoints
              dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
            end do
            
            ! Store Robin boundary condition in the cubature points
            Dcoefficients(1,icubp,iel) = dlocalData
          end do ! icubp
#endif
        end do ! iel
        
        ! Free temporal memory
        deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        deallocate(Dcoords, DcoeffAtDOF)
#endif

      else

        ! Clear coefficients for zero velocity
        Dcoefficients = 0.0_DP

      end if

      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrConvD2d_sim')
      call sys_halt()

    end select

  end subroutine transp_coeffMatBdrConvD2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatDiagSTBurgP2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IdofsAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for space-time formulation of the
    ! one-dimensional primal Burgers equation $du/dt+df(u)/dx=0$,
    ! whereby the flux function is given by $f(u)=0.5*u^2$.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IdofsAtNode

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
      ! Compute convective coefficients $k_{ii} = [0.5*u_i,1]*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          (0.5*DdataAtNode(inode)*DmatrixCoeffsAtNode(1,inode)&
          +                       DmatrixCoeffsAtNode(2,inode))
#else
      ! Compute convective coefficients $k_{ii} = -[0.5*u_i,1]*C_{ii}$
      DcoefficientsAtNode(1,inode) = -dscale*&
          (0.5*DdataAtNode(inode)*DmatrixCoeffsAtNode(1,inode)&
          +                       DmatrixCoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagSTBurgP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalSTBurgP2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IedgeList, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the
    ! one-dimensional primal Burgers equation $du/dt+df(u)/dx=0$,
    ! whereby the flux function is given by $f(u)=0.5*u^2$.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

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
      ! Compute convective coefficient $k_{ij} = [0.5*u_j,1]*C_{ji}$
      DcoefficientsAtEdge(1,iedge) = dscale*&
          (0.5_DP*DdataAtEdge(2,iedge)*DmatrixCoeffsAtEdge(1,2,iedge)&
          +                            DmatrixCoeffsAtEdge(2,2,iedge))

      ! Compute convective coefficient $k_{ji} = [0.5*u_i,1]*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (0.5_DP*DdataAtEdge(1,iedge)*DmatrixCoeffsAtEdge(1,1,iedge)&
          +                            DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -[0.5*u_j,1]*C_{ij}$
      DcoefficientsAtEdge(1,iedge) = -dscale*&
          (0.5_DP*DdataAtEdge(2,iedge)*DmatrixCoeffsAtEdge(1,1,iedge)&
          +                            DmatrixCoeffsAtEdge(2,1,iedge))

      ! Compute convective coefficient $k_{ji} = -[0.5*u_i,1]*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (0.5_DP*DdataAtEdge(1,iedge)*DmatrixCoeffsAtEdge(1,2,iedge)&
          +                            DmatrixCoeffsAtEdge(2,2,iedge))
#endif
    end do

  end subroutine transp_calcMatGalSTBurgP2d_sim
  
  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwSTBurgP2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IedgeList, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the
    ! one-dimensional primal Burgers equation $du/dt+df(u)/dx=0$,
    ! whereby the flux function is given by $f(u)=0.5*u^2$.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

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
      ! Compute convective coefficient $k_{ij} = [0.5*u_j,1]*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (0.5_DP*DdataAtEdge(2,iedge)*DmatrixCoeffsAtEdge(1,2,iedge)&
          +                            DmatrixCoeffsAtEdge(2,2,iedge))

      ! Compute convective coefficient $k_{ji} = [0.5*u_i,1]*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (0.5_DP*DdataAtEdge(1,iedge)*DmatrixCoeffsAtEdge(1,1,iedge)&
          +                            DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -[0.5*u_j,1]*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (0.5_DP*DdataAtEdge(2,iedge)*DmatrixCoeffsAtEdge(1,1,iedge)&
          +                            DmatrixCoeffsAtEdge(2,1,iedge))

      ! Compute convective coefficient $k_{ji} = -[0.5*u_i,1]*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (0.5_DP*DdataAtEdge(1,iedge)*DmatrixCoeffsAtEdge(1,2,iedge)&
          +                            DmatrixCoeffsAtEdge(2,2,iedge))
#endif

      ! Compute artificial diffusion coefficient
      !   $d_{ij} = abs(v_{ij}*0.5*(cx_{ij}-cx_{ji}) + max\{cy_{ij}, cy_{ji}\}$,
      ! where 
      !   $v_{ij} = 0.5*(f(u_j)-f(u_i))/(u_j-u_i) = 0.5*(u_i+u_j)$
      DcoefficientsAtEdge(1,iedge) = dscale*&
          abs(0.25_DP*(DdataAtEdge(1,iedge)+DdataAtEdge(2,iedge))*&
              (DmatrixCoeffsAtEdge(1,1,iedge)-DmatrixCoeffsAtEdge(1,2,iedge)))+&
          max(-DmatrixCoeffsAtEdge(2,1,iedge),-DmatrixCoeffsAtEdge(2,2,iedge))
    end do

  end subroutine transp_calcMatUpwSTBurgP2d_sim

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrSTBurgP2d_sim(rdiscretisation,&
      rform, nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

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
    ! This routine handles the primal problem for the
    ! 1D Burgers equation in space and time.
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

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    !   DquickAccess(3):     minimim parameter of the boundary component
    !   DquickAccess(4):     maximum parameter of the boundary component
    !   DquickAccess(5):     minimim parameter of the mirror boundary component
    !   DquickAccess(6):     maximum parameter of the mirror boundary component
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
    real(DP), dimension(:,:), pointer :: Dnx,Dny
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnv,dtime,dscale,dval
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrSTBurgP2d_sim')
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
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrSTBurgP2d_sim')


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
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

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
      ! $$ u=g \Rightarrow [0.5*u,1]*u\cdot{\bf n}=[0.5*g,1]*g\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

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
      ! $$ -([0.5*u,1]*u-d\nabla u)\cdot{\bf n} = -([0.5*g,1]*g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 1))
      allocate(Dnx(npointsPerElement, nelements))
      allocate(Dny(npointsPerElement, nelements))

      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
      
      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,1).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,1))
          
          ! Compute the normal velocity and impose Dirichlet boundary condition
          dnv = Dnx(ipoint,iel) * 0.5*Daux(ipoint,iel,1) + Dny(ipoint,iel)
          Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,1)
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux, Dnx, Dny)


    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([0.5*u,1]*u-d\nabla u)\cdot{\bf n} = -([0.5*g,1]*g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 2))
      allocate(Dnx(npointsPerElement, nelements))
      allocate(Dny(npointsPerElement, nelements))

      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,:1)
      call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
      
      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,2). Multiply the velocity vector [0.5*u,1]
      ! with the normal in each point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,2))

          ! Compute the normal velocity
          dnv = Dnx(ipoint,iel) * 0.5*Daux(ipoint,iel,1) + Dny(ipoint,iel)

          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
            ! Compute the prescribed normal velocity
            dnv = Dnx(ipoint,iel) * 0.5*Daux(ipoint,iel,2) + Dny(ipoint,iel)
            Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,2)
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux, Dnx, Dny)


    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrSTBurgP2d_sim')
      call sys_halt()
      
    end select
    
  end subroutine transp_coeffVecBdrSTBurgP2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrSTBurgP2d_sim(&
      rdiscretisationTrial, rdiscretisationTest, rform, nelements,&
      npointsPerElement, Dpoints, ibct, DpointPar, IdofsTrial,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

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
    ! 1D Burgers equation in space and time.
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

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    real(DP), dimension(:,:), pointer :: Dnx,Dny
    real(DP) :: dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrSTBurgP2d_sim')
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
      allocate(Dnx(npointsPerElement, nelements))
      allocate(Dny(npointsPerElement, nelements))
      
      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
      
      ! Multiply the velocity vector [0.5*u,1] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Compute the normal velocity
          dnv = Dnx(ipoint,iel) * 0.5*Daux(ipoint,iel,1) + Dny(ipoint,iel)
          
          ! Scale normal velocity by scaling parameter
          Dcoefficients(1,ipoint,iel) = -dscale * dnv
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux, Dnx, Dny)

      
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
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrSTBurgP2d_sim')

      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 1))
      allocate(Dnx(npointsPerElement, nelements))
      allocate(Dny(npointsPerElement, nelements))

      ! Evaluate the solution field in the cubature points on the boundary
      ! and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)

      ! Multiply the velocity vector [0.5*u,1] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Compute the normal velocity
          dnv = Dnx(ipoint,iel) * 0.5*Daux(ipoint,iel,1) + Dny(ipoint,iel)
          
          ! Check if we are at the primal outflow boundary
          if (dnv .gt. 0.0_DP) then
            Dcoefficients(1,ipoint,iel) = -dscale * dnv
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux, Dnx, Dny)


    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrSTBurgP2d_sim')
      call sys_halt()

    end select

  end subroutine transp_coeffMatBdrSTBurgP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatDiagSTBLevP2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IdofsAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the
    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IdofsAtNode

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
      ! Compute convective coefficient $k_{ii} = [a_i,1]*C_{ii}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtNode(1,inode) = dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtNode(1,inode)&
          +                                          DmatrixCoeffsAtNode(2,inode))
#else
      ! Compute convective coefficient $k_{ii} = -[a_i,1]*C_{ii}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtNode(1,inode) = -dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtNode(1,inode)&
          +                                          DmatrixCoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagSTBLevP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalSTBLevP2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IedgeList, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the
    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

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
      ! Compute convective coefficient $k_{ij} = [a_j,1]*C_{ji}$
      ! where $a_j=u_j/(u_^2+0.5*(1-u_j)^2)$
      DcoefficientsAtEdge(1,iedge) = dscale*&
          (uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +                                          DmatrixCoeffsAtEdge(2,2,iedge))
      
      ! Compute convective coefficient $k_{ji} = [a_i,1]*C_{ij}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +                                          DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -[a_j,1]*C_{ij}$
      ! where $a_j=u_j/(u_j^2+0.5*(1-u_j)^2)$
      DcoefficientsAtEdge(1,iedge) = -dscale*&
          (uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +                                          DmatrixCoeffsAtEdge(2,1,iedge))
      
      ! Compute convective coefficient $k_{ji} = -[a_i,1]*C_{ji}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +                                          DmatrixCoeffsAtEdge(2,2,iedge))
#endif
    end do
    
  end subroutine transp_calcMatGalSTBLevP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwSTBLevP2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IedgeList, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the
    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
    !
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

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
      ! Compute convective coefficient $k_{ij} = [a_j,1]*C_{ji}$
      ! where $a_j=u_j/(u_^2+0.5*(1-u_j)^2)$
      DcoefficientsAtEdge(2,iedge) = dscale*&
          (uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +                                          DmatrixCoeffsAtEdge(2,2,iedge))
      
      ! Compute convective coefficient $k_{ji} = [a_i,1]*C_{ij}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtEdge(3,iedge) = dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +                                          DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -[a_j,1]*C_{ij}$
      ! where $a_j=u_j/(u_j^2+0.5*(1-u_j)^2)$
      DcoefficientsAtEdge(2,iedge) = -dscale*&
          (uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))*DmatrixCoeffsAtEdge(1,1,iedge)&
          +                                          DmatrixCoeffsAtEdge(2,1,iedge))
      
      ! Compute convective coefficient $k_{ji} = -[a_i,1]*C_{ji}$
      ! where $a_i=u_i/(u_i^2+0.5*(1-u_i)^2)$
      DcoefficientsAtEdge(3,iedge) = -dscale*&
          (ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))*DmatrixCoeffsAtEdge(1,2,iedge)&
          +                                          DmatrixCoeffsAtEdge(2,2,iedge))
#endif

      ! Calculate the characteristic speed
      if (abs(ui-uj) .gt. SYS_EPSREAL_DP) then
        vij = (uj*uj/(uj*uj+0.5_DP*(1.0_DP-uj)*(1.0_DP-uj))&
              -ui*ui/(ui*ui+0.5_DP*(1.0_DP-ui)*(1.0_DP-ui)))/(uj-ui)
      else
        vij = ui*(1.0_DP-ui)/(ui*ui-0.5_DP*(1.0_DP-ui)*(1.0_DP-ui))**2
      end if

      ! Compute artificial diffusion coefficient 
      ! $d_{ij} = abs(0.5*v_{ij}*(c_{ij}-c_{ji}) + max\{cy_{ij}, cy_{ji}\}$
      DcoefficientsAtEdge(1,iedge) = dscale*&
          abs(0.5_DP*vij*(DmatrixCoeffsAtEdge(1,1,iedge)-&
                          DmatrixCoeffsAtEdge(1,2,iedge)))+&
          max(-DmatrixCoeffsAtEdge(2,1,iedge),-DmatrixCoeffsAtEdge(2,2,iedge))
    end do
    
  end subroutine transp_calcMatUpwSTBLevP2d_sim

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrSTBLevP2d_sim(rdiscretisation,&
      rform, nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

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
    ! This routine handles the primal problem for the
    ! 1D Buckley-Leverett equation in space and time.
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

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    !   DquickAccess(3):     minimim parameter of the boundary component
    !   DquickAccess(4):     maximum parameter of the boundary component
    !   DquickAccess(5):     minimim parameter of the mirror boundary component
    !   DquickAccess(6):     maximum parameter of the mirror boundary component
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
    real(DP), dimension(:,:), pointer :: Dnx,Dny
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dnv,dtime,dscale,dval
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrSTBLevP2d_sim')
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
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrSTBLevP2d_sim')


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
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

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
      ! $$ u=g \Rightarrow [a,1]*u\cdot{\bf n}=[a,1]*g\cdot{\bf n} $$
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
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

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
      ! $$ -([a,1]*u-d\nabla u)\cdot{\bf n} = -([a,1]*g)\cdot{\bf n} $$
      !
      ! where $ a = u/(u^2+0.5*(1-u)^2) $ is the velocity
      !
      ! and do not include any boundary integral into the bilinear form at all.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 1))
      allocate(Dnx(npointsPerElement, nelements))
      allocate(Dny(npointsPerElement, nelements))
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,1). Multiply the velocity vector [a,1] with
      ! the normal in each point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,1))
          
          ! Compute the normal velocity and impose Dirichlet boundary condition
          dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2 +&
                0.5*(1-Daux(ipoint,iel,1))**2) + Dny(ipoint,iel)
          Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,1)
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux, Dnx, Dny)


    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([a,1]*u-d\nabla u)\cdot{\bf n} = -([a,1]*g)\cdot{\bf n} $$
      !
      ! where $ a = u/(u^2+0.5*(1-u)^2) $ is the velocity
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 2))
      allocate(Dnx(npointsPerElement, nelements))
      allocate(Dny(npointsPerElement, nelements))

      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime

      ! Evaluate the function parser for the boundary values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,2). Multiply the velocity vector [a,1] with
      ! the normal in each point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D, ipoint, iel)

          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,2))
          
          ! Compute the normal velocity
          dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2 +&
                0.5*(1-Daux(ipoint,iel,1))**2) + Dny(ipoint,iel)

          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
            ! Compute the prescribed normal velocity
            dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,2)/(Daux(ipoint,iel,2)**2 +&
                  0.5*(1-Daux(ipoint,iel,2))**2) + Dny(ipoint,iel)
            Dcoefficients(1,ipoint,iel) = -dscale * dnv * Daux(ipoint,iel,2)
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Daux, Dnx, Dny)


    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrSTBLevP2d_sim')
      call sys_halt()
      
    end select

  end subroutine transp_coeffVecBdrSTBLevP2d_sim

    !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrSTBLevP2d_sim(&
      rdiscretisationTrial, rdiscretisationTest, rform, nelements,&
      npointsPerElement, Dpoints, ibct, DpointPar, IdofsTrial,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

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
    ! 1D Buckley-Leverett equation in space and time.
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

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    real(DP), dimension(:,:), pointer :: Dnx,Dny
    real(DP) :: dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrSTBLevP2d_sim')
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
      allocate(Daux(npointsPerElement, nelements, 2))
      allocate(Dnx(npointsPerElement, nelements))
      allocate(Dny(npointsPerElement, nelements))
      
      ! Evaluate the solution in the cubature points on the boundary
      ! and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)

      ! Multiply the velocity vector [a,1] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Compute the normal velocity
          dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2 +&
                0.5*(1-Daux(ipoint,iel,1))**2) + Dny(ipoint,iel)
          
          ! Scale normal velocity by scaling parameter
          Dcoefficients(1,ipoint,iel) = -dscale * dnv
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux, Dnx, Dny)

      
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
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrSTBLevP2d_sim')

      
    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      
      ! Allocate temporal memory
      allocate(Daux(npointsPerElement, nelements, 2))
      allocate(Dnx(npointsPerElement, nelements))
      allocate(Dny(npointsPerElement, nelements))

      ! Evaluate the solution field in the cubature points on the boundary
      ! and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
      
      ! Multiply the velocity vector [a,1] with the normal in each
      ! point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Compute the normal velocity
          dnv = Dnx(ipoint,iel) * Daux(ipoint,iel,1)/(Daux(ipoint,iel,1)**2 +&
                0.5*(1-Daux(ipoint,iel,1))**2) + Dny(ipoint,iel)
          
          ! Check if we are at the primal outflow boundary
          if (dnv .gt. 0.0_DP) then
            Dcoefficients(1,ipoint,iel) = -dscale * dnv
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux, Dnx, Dny)


    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrSTBLevP2d_sim')
      call sys_halt()

    end select

  end subroutine transp_coeffMatBdrSTBLevP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatDiagBurgP2d_sim(DdataAtNode,&
      DmatrixCoeffsAtNode, IdofsAtNode, dscale, nnodes,&
      DcoefficientsAtNode, rcollection)

!<description>
    ! This subroutine computes the diagonal convective matrix
    ! coefficients $k_{ii}$ for the primal Burger`s equation in 2D.
!</description>

!<input>
    ! Nodal solution values for all nodes under consideration
    real(DP), dimension(:), intent(in) :: DdataAtNode
    
    ! Entries of the coefficient matrices for all nodes under consideration
    real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
    
    ! Numbers of vertices and matrix entries for all nodes under consideration
    integer, dimension(:,:), intent(in) :: IdofsAtNode

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
      ! Compute convective coefficients $k_{ii} = (0.5*u_i*I)*C_{ii}$
      DcoefficientsAtNode(1,inode) = dscale*&
          0.5_DP*(DdataAtNode(inode)*DmatrixCoeffsAtNode(1,inode)&
                 +DdataAtNode(inode)*DmatrixCoeffsAtNode(2,inode))
#else
      ! Compute convective coefficients $k_{ii} = -(0.5*u_i*I)*C_{ii}$
      DcoefficientsAtNode(1,inode) = -dscale*&
          0.5_DP*(DdataAtNode(inode)*DmatrixCoeffsAtNode(1,inode)&
                 +DdataAtNode(inode)*DmatrixCoeffsAtNode(2,inode))
#endif
    end do
    
  end subroutine transp_calcMatDiagBurgP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalBurgP2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IedgeList, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the primal Burger`s equation in 2D.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

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
      ! Compute convective coefficient $k_{ij} = (0.5*u_j*I)*C_{ji}$
      DcoefficientsAtEdge(1,iedge) = dscale*0.5_DP*DdataAtEdge(2,iedge)*&
          (DmatrixCoeffsAtEdge(1,2,iedge)+DmatrixCoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji} = (0.5*u_i*I)*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = dscale*0.5_DP*DdataAtEdge(1,iedge)*&
          (DmatrixCoeffsAtEdge(1,1,iedge)+DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -(0.5*u_j*I)*C_{ij}$
      DcoefficientsAtEdge(1,iedge) = -dscale*0.5_DP*DdataAtEdge(2,iedge)*&
          (DmatrixCoeffsAtEdge(1,1,iedge)+DmatrixCoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji} = -(0.5*u_i*I)*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = -dscale*0.5_DP*DdataAtEdge(1,iedge)*&
          (DmatrixCoeffsAtEdge(1,2,iedge)+DmatrixCoeffsAtEdge(2,2,iedge))
#endif
    end do

  end subroutine transp_calcMatGalBurgP2d_sim

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwBurgP2d_sim(DdataAtEdge,&
      DmatrixCoeffsAtEdge, IedgeList, dscale, nedges,&
      DcoefficientsAtEdge, rcollection)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the primal Burger`s equation in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>

!<input>
    ! Nodal solution values for all edges under consideration
    real(DP), dimension(:,:), intent(in) :: DdataAtEdge
    
    ! Entries of the coefficient matrices for all edges under consideration
    real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
    
    ! Numbers of vertices and matrix entries for all edges under consideration
    integer, dimension(:,:), intent(in) :: IedgeList

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
      ! Compute convective coefficient $k_{ij} = (0.5*u_j*I)*C_{ji}$
      DcoefficientsAtEdge(2,iedge) = dscale*0.5_DP*DdataAtEdge(2,iedge)*&
          (DmatrixCoeffsAtEdge(1,2,iedge)+DmatrixCoeffsAtEdge(2,2,iedge))
      ! Compute convective coefficient $k_{ji} = (0.5*u_i*I)*C_{ij}$
      DcoefficientsAtEdge(3,iedge) = dscale*0.5_DP*DdataAtEdge(1,iedge)*&
          (DmatrixCoeffsAtEdge(1,1,iedge)+DmatrixCoeffsAtEdge(2,1,iedge))
#else
      ! Compute convective coefficient $k_{ij} = -(0.5*u_j*I)*C_{ij}$
      DcoefficientsAtEdge(2,iedge) = -dscale*0.5_DP*DdataAtEdge(2,iedge)*&
          (DmatrixCoeffsAtEdge(1,1,iedge)+DmatrixCoeffsAtEdge(2,1,iedge))
      ! Compute convective coefficient $k_{ji} = -(0.5*u_i*I)*C_{ji}$
      DcoefficientsAtEdge(3,iedge) = -dscale*0.5_DP*DdataAtEdge(1,iedge)*&
          (DmatrixCoeffsAtEdge(1,2,iedge)+DmatrixCoeffsAtEdge(2,2,iedge))
#endif

      ! Compute artificial diffusion coefficient 
      ! $d_{ij} = abs(v_{ij}*I*0.5*(c_{ij}-c_{ji})$,
      ! where $v_{ij} = 0.5*(f(u_j)-f(u_i))/(u_j-u_i) = 0.5*(u_i+u_j)$
      DcoefficientsAtEdge(1,iedge) = dscale*&
          abs(0.25_DP*(DdataAtEdge(1,iedge)+DdataAtEdge(2,iedge))*&
          (DmatrixCoeffsAtEdge(1,1,iedge)-DmatrixCoeffsAtEdge(1,2,iedge)+&
           DmatrixCoeffsAtEdge(2,1,iedge)-DmatrixCoeffsAtEdge(2,2,iedge)))
    end do

  end subroutine transp_calcMatUpwBurgP2d_sim

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrBurgP2d_sim(rdiscretisation,&
      rform, nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

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
    ! This routine handles the constant velocities in the primal problem.
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

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    !   DquickAccess(3):     minimim parameter of the boundary component
    !   DquickAccess(4):     maximum parameter of the boundary component
    !   DquickAccess(5):     minimim parameter of the mirror boundary component
    !   DquickAccess(6):     maximum parameter of the mirror boundary component
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
    type(t_boundaryRegion), pointer :: p_rboundaryRegionMirror
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(:), pointer :: Domega,DcoeffAtDOF
    real(DP), dimension(:,:), pointer :: DcubPtsRef,Dbas
    real(DP), dimension(:,:), pointer :: Dnx,Dny,DpointParMirror
    real(DP), dimension(:,:,:), pointer :: Dcoords,Daux
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dlocalData,dnv,dscale,dtime,dval
    real(DP) :: dminParam,dmaxParam,dminParamMirror,dmaxParamMirror
    integer :: ccubType,ibdrtype,icubp,iel,ipoint,isegment,ivt,npoints,nve

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffVecBdrBurgP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access vectors
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Set pointers
    call lsysbl_getbase_double(p_rsolution, p_Ddata)
#endif

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the cubature rule
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    ccubType = rcollection%IquickAccess(3)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Evaluate one-dimensional basis functions on the boundary edge
    if (npointsPerElement .ne. cub_igetNumPts(ccubType)) then
      call output_line('Type of cubature rule at boundary mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrBurgP2d_sim')
      call sys_halt()
    else
      ! How many DOFs are located at the boundary? This should be made
      ! more flexible by checking the type of element. For the time
      ! being, only linear and bilinear finite elements are supported
      npoints = 2
      
      ! How many vertices per element do we have?
      nve = elem_igetNVE(rdomainIntSubset%celement)
      
      ! Allocate temporal memory for one-dimensional
      ! cubature along the boundary edge
      allocate(Dbas(npoints,npointsPerElement))
      allocate(Domega(npointsPerElement))
      allocate(DcubPtsRef(1,npointsPerElement))

      ! Get the coordinates of the cubature points and the
      ! corresponding weights for the given cubature rule
      call cub_getCubature(ccubType, DcubPtsRef, Domega)
      
      ! Evaluate the one-dimensional basis functions 
      ! in the cubature points on the boundary
      do icubp = 1, npointsPerElement
        Dbas(1,icubp) = 0.5_DP*(1.0_DP-DcubPtsRef(1,icubp))
        Dbas(2,icubp) = 0.5_DP*(1.0_DP+DcubPtsRef(1,icubp))
      end do

      ! Deallocate temporal memory which is no longer required
      deallocate(DcubPtsRef,Domega)
    end if    
#else
    ! Boundary values are evaluated directly at the cubature points
    npoints = npointsPerElement
#endif

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
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffVecBdrBurgP2d_sim')


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
      ! Dcoefficients(1,:,:).
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement

          ! Set values for function parser
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)

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
      ! $$ u=g \Rightarrow ([0.5*u*I]*u\cdot{\bf n}=([0.5*g*I]*g)\cdot{\bf n} $$
      !
      ! The diffusive part is included into the bilinear form.

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,1), DcoeffAtDOF(npoints))
#endif
      
      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      do iel = 1, nelements
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,1) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,1) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
#endif
        
        do ipoint = 1, npoints
          
          ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,1)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif

          ! Evaluate function parser for Dirichlet value
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Impose Dirichlet value via penalty method
          DcoeffAtDOF(ipoint) = dscale * dval * BDRC_DIRICHLET_PENALTY
#else
          ! Impose Dirichlet value via penalty method
          Dcoefficients(1,ipoint,iel) = dscale * dval * BDRC_DIRICHLET_PENALTY
#endif
        end do
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Deallocate temporal memory
      deallocate(Dcoords, DcoeffAtDOF)
#endif


    case (BDRC_ROBIN)
      !-------------------------------------------------------------------------
      ! Robin boundary conditions:
      !
      ! Evaluate coefficients for both the convective and the diffusive
      ! part of the linear form 
      !
      ! $$ -([0.5*u*I]*u-d\nabla u)\cdot{\bf n} = -([0.5*g*I]*g)\cdot{\bf n} $$
      !
      ! and do not include any boundary integral into the bilinear form at all.
      
      ! Allocate temporal memory
      allocate(Dnx(npoints,nelements), Dny(npoints,nelements))
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
      
      ! We need the physical coordinates of the DOFs on the boundary,
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
      ! Calculate the normal vectors in cubature on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
      
      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Evaluate the function parser for the boundary values in the
      ! cubature points or the DOFs on the boundary. Multiply the
      ! velocity vector with the normal in each point to get the
      ! normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, dval)
          
          ! Compute the normal velocity
          dnv = 0.5_DP * (Dnx(ipoint,iel)+Dny(ipoint,iel)) * dval

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          ! Store the Robin boundary condition at the DOF on the boundary
          DcoeffAtDOF(ipoint) = -dscale * dnv * dval
#else
          ! Store the Robin boundary condition at the cubature point on the boundary
          Dcoefficients(1,ipoint,iel) = -dscale * dnv * dval
#endif
        end do
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do
      
      ! Deallocate temporal memory
      deallocate(Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      deallocate(Dcoords, DcoeffAtDOF)
#endif


    case(BDRC_FLUX)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s prescribed at the inlet):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([0.5*u*I]*u-d\nabla u)\cdot{\bf n} = -([0.5*g*I]g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.

      ! Allocate temporal memory
      allocate(Dnx(npoints,nelements),Dny(npoints,nelements),&
               Daux(npoints,nelements,1))
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
      
      ! We need the solution values and the physical coordinates of the DOFs 
      ! on the boundary, so compute them and store the result in Daux
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(1,iel,1) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(2,iel,1) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
      ! Evaluate the solution vector in the cubature points on the
      ! boundary and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif     
      
      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Evaluate the function parser for the boundary values in the
      ! cubature points or the DOFS on the boundary. Multiply the
      ! velocity vector with the normal in each point to get the
      ! normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Set values for function parser
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
          Dvalue(1:NDIM2D) = Dcoords(1:NDIM2D,ipoint,iel)
#else
          Dvalue(1:NDIM2D) = Dpoints(1:NDIM2D,ipoint,iel)
#endif
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment, Dvalue, dval)
          
          ! Compute the normal velocity
          dnv = 0.5_DP * (Dnx(ipoint,iel) + Dny(ipoint,iel)) * Daux(ipoint,iel,1)
          
          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = -dscale * 0.5_DP*(Dnx(ipoint,iel) +&
                                                    Dny(ipoint,iel)) * dval * dval
#else
            Dcoefficients(1,ipoint,iel) = -dscale * 0.5_DP*(Dnx(ipoint,iel) +&
                                                            Dny(ipoint,iel)) * dval * dval
#endif
          else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = 0.0_DP
#else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
          end if
        end do
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do
      
      ! Deallocate temporal memory
      deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      deallocate(Dcoords, DcoeffAtDOF)
#endif


    case(BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Periodic/Antiperiodic boundary conditions (Flux boundary conditions):
      !
      ! Evaluate coefficient for both the convective and diffusive
      ! part for the linear form at the inflow boundary part.
      !
      ! $$ -([0.5*u*I]*u-d\nabla u)\cdot{\bf n} = -([0.5*g*I]*g)\cdot{\bf n} $$
      !
      ! The boundary integral at the outflow boundary is included
      ! into the bilinear form.
      
      ! Get mirrored boundary region from collection structure
      p_rboundaryRegionMirror => collct_getvalue_bdreg(rcollection,&
          'rboundaryRegionMirror', ssectionName=trim(rcollection%SquickAccess(1)))
      
      ! Get minimum/maximum parameter values from collection structure
      dminParam = rcollection%DquickAccess(3)
      dmaxParam = rcollection%DquickAccess(4)
      dminParamMirror = rcollection%DquickAccess(5)
      dmaxParamMirror = rcollection%DquickAccess(6)
      
      ! Allocate temporal memory
      allocate(Dnx(npoints,nelements), Dny(npoints,nelements))
      allocate(DpointParMirror(npoints,nelements), Daux(npoints,nelements,2))
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
      
      ! We need the solution values and the physical coordinates of the DOFs 
      ! on the boundary, so compute them and store the result in Daux
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(1,iel,1) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(2,iel,1) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else  
      ! Evaluate the solution vector in the cubature points on
      ! the boundary and store the result in Daux(:,:,1)
      call fevl_evaluate_sim(DER_FUNC2D, Daux(:,:,1),&
          p_rsolution%RvectorBlock(1), Dpoints, &
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)

      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
      
      if (iand(ibdrtype, BDRC_TYPEMASK) .eq. BDRC_PERIODIC) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Rescale parameter values DedgePosition on the boundary
        ! segment where to compute the boundary conditions into
        ! parameter values on the mirror boundary region
        call mprim_linearRescale(rdomainIntSubset%p_DedgePosition, dminParam,&
            dmaxParam, dmaxParamMirror, dminParamMirror, DpointParMirror)
#else
        ! Rescale parameter values DpointPar on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the mirror boundary region
        call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
            dmaxParamMirror, dminParamMirror, DpointParMirror)
#endif
      else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Rescale parameter values DedgePosition on the boundary
        ! segment where to compute the boundary conditions into
        ! parameter values on the reversely mirror boundary region
        call mprim_linearRescale(rdomainIntSubset%p_DedgePosition, dminParam,&
            dmaxParam, dmaxParamMirror, dminParamMirror, DpointParMirror)
#else
        ! Rescale parameter values DpointPar on the boundary segment
        ! where to compute the boundary conditions into parameter
        ! values on the reversely mirror boundary region
        call mprim_linearRescale(DpointPar, dminParam, dmaxParam,&
            dmaxParamMirror, dminParamMirror, DpointParMirror)
#endif
      end if
      
      ! Evaluate the solution in the cubature points on the mirrored
      ! boundary and store the result in Daux(:,:,2)
      call doEvaluateAtBdr2d(DER_FUNC, npoints*nelements,&
          Daux(:,:,2), p_rsolution%RvectorBlock(1), DpointParMirror,&
          ibct, BDR_PAR_LENGTH, p_rboundaryRegionMirror)

      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute the normal velocity
          dnv = 0.5_DP * (Dnx(ipoint,iel) + Dny(ipoint,iel)) * Daux(ipoint,iel,1)
          
          ! Check if we are at the primal inflow boundary
          if (dnv .lt. 0.0_DP) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = -dscale * 0.5_DP * (Dnx(ipoint,iel) +&
                                  Dny(ipoint,iel)) * Daux(ipoint,iel,2)**2
#else
            Dcoefficients(1,ipoint,iel) = -dscale * 0.5_DP * (Dnx(ipoint,iel) +&
                                          Dny(ipoint,iel)) * Daux(ipoint,iel,2)**2
#endif
          else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = 0.0_DP
#else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
          end if
        end do
        
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do
      
      ! Deallocate temporal memory
      deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Deallocate temporal memory
      deallocate(Dcoords, DcoeffAtDOF)
#endif
      

    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffVecBdrBurgP2d_sim')
      call sys_halt()
      
    end select

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Deallocate temporal memory
    deallocate(Dbas)
#endif

  contains

    ! Here come the working routines

    !***************************************************************************
    ! Evaluate th solution vector at some boundary points given in
    ! terms of their parameter values. This ugly trick is necessary
    ! since we have to pass the 2d-array Dvalues and DpointsPar to a
    ! subroutine which accepts only 1d-arrays.
    ! ***************************************************************************
    
    subroutine doEvaluateAtBdr2d(iderType, n, Dvalues, rvectorScalar,&
        DpointsPar, ibdc, cparType, rboundaryRegion)
      
      integer, intent(in) :: iderType,ibdc,cparType,n
      real(DP), dimension(n), intent(in) :: DpointsPar
      type(t_vectorScalar), intent(in) :: rvectorScalar
      type(t_boundaryRegion), intent(in) :: rboundaryRegion
      
      real(DP), dimension(n), intent(out) :: Dvalues
      
      call fevl_evaluateBdr2D(iderType, Dvalues, rvectorScalar,&
          DpointsPar, ibdc, cparType, rboundaryRegion)
      
    end subroutine doEvaluateAtBdr2d

  end subroutine transp_coeffVecBdrBurgP2d_sim

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrBurgP2d_sim( rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

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
    ! This routine handles the constant velocities in the primal problem.
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

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

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
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(:), pointer :: Domega,DcoeffAtDOF
    real(DP), dimension(:,:), pointer :: DcubPtsRef,Dbas
    real(DP), dimension(:,:), pointer :: Dnx,Dny,Daux
    real(DP), dimension(:,:,:), pointer :: Dcoords
    real(DP) :: dlocalData,dnv,dscale,dtime
    integer :: ccubType,ibdrtype,icubp,iel,ipoint,isegment,ivt,npoints,nve

    ! REMARK: This subroutine makes use of a lot of preprocessor flags
    ! to distinguish between the classical boundary integal evaluation,
    ! that is, the boundary integrals are approximated by evaluating
    ! the FE-function at the cubature points and computing the flux
    ! based on the FE-function in the cubature points.
    ! This approach is not fully compatible with the group finite
    ! element formulation which approximates the fluxes in the same
    ! way as the approxitame solution in the degrees of freedom. In
    ! this case, the fluxes are evaluated in the degrees of freedom
    ! and interpolated to the cubature points on the boundary afterwards.

#ifndef TRANSP_USE_IBP
    call output_line('Application must be compiled with flag &
        &-DTRANSP_USE_IBP if boundary conditions are imposed in weak sense',&
        OU_CLASS_ERROR, OU_MODE_STD, 'transp_coeffMatBdrBurgP2d_sim')
    call sys_halt()
#endif

    ! This subroutine assumes that the first quick access vector
    ! points to the solution vector
    p_rsolution => rcollection%p_rvectorQuickAccess1


#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Set pointers
    call lsysbl_getbase_double(p_rsolution, p_Ddata)
#endif
    
    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first three quick access integer values hold the type of
    ! boundary condition, the segment number and the cubature rule
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)
    ccubType = rcollection%IquickAccess(3)

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
    ! Evaluate one-dimensional basis functions on the boundary edge
    ! which will be used to (a) interpolate the FE-functions evaluated
    ! at the degrees of freedom into the cubature points for numerical
    ! integration and to redistribute the function values given at the
    ! cubature points to the degrees of freefom.
    if (npointsPerElement .ne. cub_igetNumPts(ccubType)) then
      call output_line('Type of cubature rule at boundary mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrBurgP2d_sim')
      call sys_halt()
    else
      ! How many DOFs are located at the boundary? This should be made
      ! more flexible by checking the type of element. For the time
      ! being, only linear and bilinear finite elements are supported
      npoints = 2
      
      ! How many vertices per element do we have?
      nve = elem_igetNVE(rdomainIntSubset%celement)
      
      ! Allocate temporal memory for one-dimensional
      ! cubature along the boundary edge
      allocate(Dbas(npoints,npointsPerElement))
      allocate(Domega(npointsPerElement))
      allocate(DcubPtsRef(1,npointsPerElement))

      ! Get the coordinates of the cubature points and the
      ! corresponding weights for the given cubature rule
      call cub_getCubature(ccubType, DcubPtsRef, Domega)
      
      ! Evaluate the one-dimensional basis functions 
      ! in the cubature points on the boundary
      do icubp = 1, npointsPerElement
        Dbas(1,icubp) = 0.5_DP*(1.0_DP-DcubPtsRef(1,icubp))
        Dbas(2,icubp) = 0.5_DP*(1.0_DP+DcubPtsRef(1,icubp))
      end do

      ! Deallocate temporal memory which is no longer required
      deallocate(DcubPtsRef,Domega)
    end if    
#else
    ! Boundary values are evaluated directly at the cubature points
    npoints = npointsPerElement
#endif

    ! What type of boundary conditions are we?
    select case(iand(ibdrtype, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN)
      !-------------------------------------------------------------------------
      ! (In-)Homogeneous Neumann boundary conditions:
      ! Assemble the convective part of the boundary integral

      ! Allocate temporal memory
      allocate(Dnx(npoints,nelements), Dny(npoints,nelements),&
               Daux(npoints,nelements))
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
      
      ! We need the solution values and the physical coordinates of the DOFs 
      ! on the boundary, so compute them and store the result in Daux
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(1,iel) = p_Ddata(ivt)
                
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(2,iel) = p_Ddata(ivt)

        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
      ! Evaluate the solution vector in the cubature points on the
      ! boundary and store the result in Daux
      call fevl_evaluate_sim(DER_FUNC2D, Daux,&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
      
      ! Multiply the velocity vector [0.5*u*I] with the normal in
      ! each point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute the normal velocity
          dnv = 0.5_DP * (Dnx(ipoint,iel) + Dny(ipoint,iel)) * Daux(ipoint,iel)
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY      
          ! Scale normal velocity by scaling parameter
          DcoeffAtDOF(ipoint) = -dscale * dnv
#else
          ! Scale normal velocity by scaling parameter
          Dcoefficients(1,ipoint,iel) = -dscale * dnv
#endif
        end do

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do
        
      ! Free temporal memory
      deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      deallocate(Dcoords, DcoeffAtDOF)
#endif
      
      
    case (BDRC_DIRICHLET)
      !-------------------------------------------------------------------------
      ! Dirichlet boundary conditions:
      ! Impose penalty parameter

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
          OU_CLASS_WARNING,OU_MODE_STD,'transp_coeffMatBdrBurgP2d_sim')


    case(BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      !-------------------------------------------------------------------------
      ! Flux boundary conditions (Robin bc`s at the outlet)
      ! Assemble the convective part of the boundary integral at the outflow
      !
      ! The convective part of the boundary integral at the outflow is
      ! likewise assembled for periodic and antiperiodic boundary conditions

      ! Allocate temporal memory
      allocate(Dnx(npoints,nelements), Dny(npoints,nelements),&
               Daux(npoints,nelements))
      
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      ! Allocate temporal memory
      allocate(Dcoords(NDIM2D,npoints,nelements), DcoeffAtDOF(npoints))
     
      ! We need the solution values and the physical coordinates of the DOFs 
      ! on the boundary, so compute them and store the result in Daux
      do iel = 1, nelements
        ! Get global DOF of first endpoints
        ipoint = rdomainIntSubset%p_IelementOrientation(iel)
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(1,iel) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,1,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
        
        ! Get global DOF of second endpoints
        ipoint = mod(rdomainIntSubset%p_IelementOrientation(iel),nve)+1
        ivt    = IdofsTest(ipoint,iel)
        
        ! Store solution data
        Daux(2,iel) = p_Ddata(ivt)
        
        ! Store vertex coordinate
        Dcoords(1:NDIM2D,2,iel) = rdomainIntSubset%p_Dcoords(1:NDIM2D,ipoint,iel)
      end do
      
      ! Calculate the normal vectors in DOFs on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dcoords, Dnx, Dny, 1)
#else
      ! Evaluate the velocity field (if any) in the cubature points on
      ! the boundary and store the result in Daux
      call fevl_evaluate_sim(DER_FUNC2D, Daux,&
          p_rsolution%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the normal vectors in the cubature points on the boundary
      call boundary_calcNormalVec2D(Dpoints, Dpoints, Dnx, Dny, 1)
#endif
      
      ! Multiply the velocity vector [0.5*u*I] with the normal in
      ! each point to get the normal velocity.
      do iel = 1, nelements
        do ipoint = 1, npoints
          
          ! Compute the normal velocity
          dnv = 0.5_DP * (Dnx(ipoint,iel) + Dny(ipoint,iel)) * Daux(ipoint,iel)

          ! Check if we are at the primal outflow boundary
          if (dnv .gt. 0.0_DP) then
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = -dscale * dnv
#else
            Dcoefficients(1,ipoint,iel) = -dscale * dnv
#endif
          else
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
            DcoeffAtDOF(ipoint) = 0.0_DP
#else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
#endif
          end if
        end do

#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
        ! Loop over the cubature points and interpolate the Robin
        ! boundary conditions from the DOFs to the cubature points,
        ! where they are needed by the linear form assembly routine
        do icubp = 1, npointsPerElement
          
          dlocalData = 0.0_DP
          
          ! Loop over the DOFs and interpolate the Robin boundary conditions
          do ipoint = 1, npoints
            dlocalData = dlocalData + Dbas(ipoint,icubp)*DcoeffAtDOF(ipoint)
          end do
          
          ! Store Robin boundary condition in the cubature points
          Dcoefficients(1,icubp,iel) = dlocalData
        end do
#endif
      end do
      
      ! Free temporal memory
      deallocate(Daux, Dnx, Dny)
#ifdef TRANSP_USE_GFEM_AT_BOUNDARY
      deallocate(Dcoords, DcoeffAtDOF)
#endif

      
    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_coeffMatBdrBurgP2d_sim')
      call sys_halt()
      
    end select

  end subroutine transp_coeffMatBdrBurgP2d_sim

end module transport_callback2d
