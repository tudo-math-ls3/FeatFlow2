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
!# 1.) transp_setVariable2d
!#     -> Sets global variables for external data, e.g., velocity fields in 2D
!#
!# 2.) transp_hadaptCallback2d
!#     -> Performs application specific tasks in the adaptation algorithm in 2D
!#
!# 3.) transp_refFuncBdrInt2d
!#     -> Callback routine for the evaluation of the boundary integral
!#        of the target functional for goal-oriented error estimation
!#
!# 4.) transp_errorBdrInt2d
!#     -> Callback routine for the evaluation of the boundary integral
!#        of the error in the target functional for goal-oriented
!#        error estimation
!#
!# 5.) transp_weightFuncBdrInt2d
!#     -> Callback routine for the evaluation of the weights in
!#        the boundary integral of the target functional for
!#        goal-oriented error estimation
!#
!#
!# ****************************************************************************
!#
!# The following routines for linear velocity case are available:
!#
!# 1.) transp_calcMatGalConvectionP2d
!#     transp_calcMatUpwConvectionP2d
!#     -> Calculates the transport coefficients for linear convection in 2D
!#
!# 2.) transp_calcMatGalConvectionD2d
!#     transp_calcMatUpwConvectionD2d
!#     -> Calculates the transport coefficients for linear convection in 2D
!#
!# 3.) transp_coeffVecBdrConvectionP2d
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 4.) transp_coeffVecBdrConvectionD2d
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 5.) transp_coeffMatBdrConvectionP2d
!#     -> Calculates the coefficients for the bilinear form in 2D
!#
!# 6.) transp_coeffMatBdrConvectionD2d
!#     -> Calculates the coefficients for the bilinear form in 2D
!#
!#
!# ****************************************************************************
!#
!# The following routines for Burgers` equation 
!# in space-time are available:
!#
!# 1.) transp_calcMatGalSTBurgersP2d
!#     transp_calcMatUpwSTBurgersP2d
!#     -> Calculates the transport coefficients for 
!#        Burgers` equation in space-time
!#
!# 2.) ...
!#
!# 3.) transp_coeffVecBdrSTBurgersP2d
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 4.) ...
!#
!# 5.) transp_coeffMatBdrSTBurgersP2d
!#     -> Calculates the coefficients for the bilinear form in 2D
!#
!# ****************************************************************************
!#
!# The following routines for the Buckley-Leverett
!#  equation in space-time are available:
!#
!# 1.) transp_calcMatGalSTBuckLevP2d
!#     transp_calcMatUpwSTBuckLevP2d
!#     -> Calculates the transport coefficients for 
!#        Buckley-Leverett equation in space-time
!#
!# 2.) ...
!#
!# 3.) transp_coeffVecBdrSTBuckLevP2d
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 4.) ...
!#
!# 5.) transp_coeffMatBdrSTBuckLevP2d
!#     -> Calculates the coefficients for the bilinear form in 2D
!#
!# ****************************************************************************
!#
!# The following routines for the Burgers` equation in 2D are available:
!#
!# 1.) transp_calcMatGalBurgersP2d
!#     transp_calcMatUpwBurgersP2d
!#     -> Calculates the transport coefficients for 
!#        Burgers` equation in 2D
!#
!# 2.) ...
!#
!# 3.) transp_coeffVecBdrBurgersP2d
!#      -> Calculates the coefficients for the linear form in 2D
!#
!# 4.) ...
!#
!# 5.) transp_coeffMatBdrBurgersP2d
!#     -> Calculates the coefficients for the bilinear form in 2D
!#
!# </purpose>
!##############################################################################

module transport_callback2d

  use collection
  use derivatives
  use fsystem
  use genoutput
  use hadaptaux
  use linearsystemscalar
  use linearsystemblock
  use storage

  use flagship_callback

  implicit none

  private
  public :: transp_calcMatGalConvectionD2d
  public :: transp_calcMatUpwConvectionD2d
  public :: transp_calcMatGalSTBuckLevP2d
  public :: transp_calcMatUpwSTBuckLevP2d
  public :: transp_calcMatGalBurgersP2d
  public :: transp_calcMatUpwBurgersP2d
  public :: transp_calcMatGalSTBurgersP2d
  public :: transp_calcMatUpwSTBurgersP2d
  public :: transp_calcMatGalConvectionP2d
  public :: transp_calcMatUpwConvectionP2d
  public :: transp_coeffMatBdrConvectionD2d
  public :: transp_coeffMatBdrSTBuckLevP2d
  public :: transp_coeffMatBdrBurgersP2d
  public :: transp_coeffMatBdrSTBurgersP2d
  public :: transp_coeffMatBdrConvectionP2d
  public :: transp_coeffVecBdrConvectionD2d
  public :: transp_coeffVecBdrSTBuckLevP2d
  public :: transp_coeffVecBdrBurgersP2d
  public :: transp_coeffVecBdrSTBurgersP2d
  public :: transp_coeffVecBdrConvectionP2d
  public :: transp_errorBdrInt2d
  public :: transp_hadaptCallback2d
  public :: transp_refFuncBdrInt2d
  public :: transp_setVariable2d
  public :: transp_weightFuncBdrInt2d

!<globals>

  !*****************************************************************
  ! Pointers to external data vectors.
  !
  ! Using global variables is not good programming style but it is the
  ! only way to allow for an efficient access to the velocity data
  ! from within the callback routines which are called repeatedly

  real(DP), dimension(:), pointer, save :: p_Dvariable1 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable2 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable3 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable4 => null()
  real(DP), dimension(:), pointer, save :: p_Dvariable5 => null()

!</globals>

contains

  !*****************************************************************************

!<subroutine>

  subroutine transp_setVariable2d(rvector, ivariable)

!<description>
    ! This subroutine sets one of the the global pointers to the given vector.
!</description>

!<input>
    ! scalar vector
    type(t_vectorScalar), intent(in) :: rvector

    ! variable number
    integer, intent(in) :: ivariable
!</input>
!</subroutine>

    select case(ivariable)
    case (1)
      call lsyssc_getbase_double(rvector, p_Dvariable1)
    case (2)
      call lsyssc_getbase_double(rvector, p_Dvariable2)
    case (3)
      call lsyssc_getbase_double(rvector, p_Dvariable3)
    case (4)
      call lsyssc_getbase_double(rvector, p_Dvariable4)
    case (5)
      call lsyssc_getbase_double(rvector, p_Dvariable5)
    case DEFAULT
      call output_line('Invalid variable number!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_setVariable2d')
      call sys_halt()
    end select
    
  end subroutine transp_setVariable2d

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
    ! Collection
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

    
    case DEFAULT
      ! Call the general callback function
      call flagship_hadaptCallback2d(iOperation, rcollection)

    end select

  end subroutine transp_hadaptCallback2d
  
  !*****************************************************************************

!<subroutine>

  subroutine transp_refFuncBdrInt2d(cderivative, rdiscretisation,&
      DpointsRef, Dpoints, ibct, DpointPar, Ielements, Dvalues,&
      rcollection)

    use basicgeometry
    use boundary
    use collection
    use domainintegration
    use feevaluation
    use fparser
    use fsystem    
    use scalarpde
    use spatialdiscretisation
    use triangulation

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
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(:,:,:), pointer :: Dcoefficients
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dt,dminPar,dmaxPar,dnx,dny,dtime
    integer :: iel,ipoint,icomp1,icomp2,icomp,ndim


    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(1)))
    
    ! This subroutine assumes that the first quick access integer
    ! value holds the number of the reference function.  Moreover,
    ! quick access interger values 3 and 4 hold the numbers of the
    ! functions to be evaluated for the x-velocity and y-velocity.
    icomp  = rcollection%IquickAccess(1)
    icomp1 = rcollection%IquickAccess(3)
    icomp2 = rcollection%IquickAccess(4)

    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    dtime = rcollection%DquickAccess(1)
    
    ! Initialize values
    Dvalue = 0.0_DP
    Dvalue(NDIM3D+1) = dtime

    ! Set number of spatial dimensions
    ndim = size(Dpoints, 1)

    ! Allocate temporal memory
    allocate(Dcoefficients(size(Dvalues,1), size(Dvalues,2), 3))

    ! Evaluate the reference function and the exact velocities in the
    ! cubature points on the boundary and store the result in
    ! Dcoefficients(:,:,1:3).
    do iel = 1, size(Ielements)
      do ipoint = 1, ubound(Dpoints,2)
        
        ! Set values for function parser
        Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

        ! Evaluate function parser
        call fparser_evalFunction(p_rfparser, icomp,  Dvalue,&
            Dcoefficients(ipoint,iel,1))
        call fparser_evalFunction(p_rfparser, icomp1, Dvalue,&
            Dcoefficients(ipoint,iel,2))
        call fparser_evalFunction(p_rfparser, icomp2, Dvalue,&
            Dcoefficients(ipoint,iel,3))
      end do
    end do

    ! Get the minimum and maximum parameter value. The point with the minimal
    ! parameter value is the start point of the interval, the point with the
    ! maximum parameter value the endpoint.
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    do iel = 1, size(Ielements)
      do ipoint = 1, ubound(Dpoints,2)
        dminPar = min(DpointPar(ipoint,iel), dminPar)
        dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
      end do
    end do

    ! Multiply the velocity vector with the normal in each point
    ! to get the normal velocity.
    do iel = 1, size(Ielements)
      do ipoint = 1, ubound(Dpoints,2)

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
          call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
              ibct, dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
        else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
          ! End point
          call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
              ibct, dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
        else
          ! Inner point
          call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
              ibct, dt, dnx, dny, cparType=BDR_PAR_LENGTH)
        end if
        
        ! Compute the expression from the data stored in Dcoefficients
        !
        !    u * (v x n)
        !
        ! in each cubature point on each elements

        Dvalues(ipoint,iel) = Dcoefficients(ipoint,iel,1) *&
                              (dnx * Dcoefficients(ipoint,iel,2) +&
                               dny * Dcoefficients(ipoint,iel,3))
      end do
    end do
    
    ! Free temporal memory
    deallocate(Dcoefficients)

  end subroutine transp_refFuncBdrInt2d

  !*****************************************************************************

!<subroutine>

  subroutine transp_errorBdrInt2d(cderivative, rdiscretisation,&
      DpointsRef, Dpoints, ibct, DpointPar, Ielements, Dvalues,&
      rcollection)

    use basicgeometry
    use boundary
    use collection
    use domainintegration
    use feevaluation
    use fparser
    use fsystem    
    use scalarpde
    use spatialdiscretisation
    use triangulation

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
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock), pointer :: p_rsolution, p_rvelocity
    real(DP), dimension(:,:,:), pointer :: Dcoefficients
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dt,dminPar,dmaxPar,dnx,dny,dtime
    integer :: iel,ipoint,icomp1,icomp2,icomp,ndim


    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access vector
    ! points to the primal solution vector and the second quick access
    ! vector points to the velocity vector
    p_rsolution => rcollection%p_rvectorQuickAccess1
    p_rvelocity => rcollection%p_rvectorQuickAccess2
    
    ! Evaluate the FE function in the cubature points on the boundary
    call fevl_evaluate_sim1(DER_FUNC, Dvalues, p_rsolution%RvectorBlock(1),&
                            Dpoints, Ielements, DpointsRef)

    ! Allocate temporal memory
    allocate(Dcoefficients(size(Dvalues,1), size(Dvalues,2), 5))
    
    ! Evaluate the velocity field in the cubature points on the boundary
    ! and store the result in Dcoefficients(:,:,1:2)
    call fevl_evaluate_sim1(DER_FUNC, Dcoefficients(:,:,1),&
        p_rvelocity%RvectorBlock(1), Dpoints, Ielements, DpointsRef)
    call fevl_evaluate_sim1(DER_FUNC, Dcoefficients(:,:,2),&
        p_rvelocity%RvectorBlock(2), Dpoints, Ielements, DpointsRef)

    ! This subroutine assumes that the first quick access integer
    ! value holds the number of the reference function.  Moreover,
    ! quick access interger values 3 and 4 hold the numbers of the
    ! functions to be evaluated for the x-velocity and y-velocity
    icomp  = rcollection%IquickAccess(1)
    icomp1 = rcollection%IquickAccess(3)
    icomp2 = rcollection%IquickAccess(4)

    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    dtime = rcollection%DquickAccess(1)
    
    ! Initialize values
    Dvalue = 0.0_DP
    Dvalue(NDIM3D+1) = dtime

    ! Set number of spatial dimensions
    ndim = size(Dpoints, 1)

    ! Evaluate the reference function and the exact velocities in the
    ! cubature points on the boundary and store the result in
    ! Dcoefficients(:,:,3:5).
    do iel = 1, size(Ielements)
      do ipoint = 1, ubound(Dpoints,2)
        
        ! Set values for function parser
        Dvalue(1:ndim) = Dpoints(:, ipoint, iel)

        ! Evaluate function parser
        call fparser_evalFunction(p_rfparser, icomp,  Dvalue,&
            Dcoefficients(ipoint,iel,3))
        call fparser_evalFunction(p_rfparser, icomp1, Dvalue,&
            Dcoefficients(ipoint,iel,4))
        call fparser_evalFunction(p_rfparser, icomp2, Dvalue,&
            Dcoefficients(ipoint,iel,5))
      end do
    end do

    ! Get the minimum and maximum parameter value. The point with the minimal
    ! parameter value is the start point of the interval, the point with the
    ! maximum parameter value the endpoint.
    dminPar = DpointPar(1,1)
    dmaxPar = DpointPar(1,1)
    do iel = 1, size(Ielements)
      do ipoint = 1, ubound(Dpoints,2)
        dminPar = min(DpointPar(ipoint,iel), dminPar)
        dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
      end do
    end do

    ! Multiply the velocity vector with the normal in each point
    ! to get the normal velocity.
    do iel = 1, size(Ielements)
      do ipoint = 1, ubound(Dpoints,2)

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
          call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
              ibct, dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
        else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
          ! End point
          call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
              ibct, dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
        else
          ! Inner point
          call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
              ibct, dt, dnx, dny, cparType=BDR_PAR_LENGTH)
        end if
        
        ! Compute the expression from the data stored in Dcoefficients
        !
        !    u * (v x n) - u_h * (v_h x n)
        !
        ! in each cubature point on each elements

        Dvalues(ipoint,iel) = Dcoefficients(ipoint,iel,3) *&
                              (dnx * Dcoefficients(ipoint,iel,4) +&
                               dny * Dcoefficients(ipoint,iel,5))-&
                              Dvalues(ipoint,iel) *&
                              (dnx * Dcoefficients(ipoint,iel,1) +&
                               dny * Dcoefficients(ipoint,iel,2))
      end do
    end do

    ! Free temporal memory
    deallocate(Dcoefficients)

  end subroutine transp_errorBdrInt2d

  !*****************************************************************************

!<subroutine>

  subroutine transp_weightFuncBdrInt2d(rdiscretisation, DpointsRef,&
      Dpoints, ibct, DpointPar, Ielements, Dvalues, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fparser
    use fsystem    
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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
    integer :: ipoint, iel, ndim, icomp


    ! Initialize values
    Dvalue = 0.0_DP
    
    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
                                       trim(rcollection%SquickAccess(1)))
   
    ! Moreover, this subroutine assumes that the second quick access integer
    ! value holds the number of the function to be evaluated
    icomp = rcollection%IquickAccess(2)
    
    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    Dvalue(NDIM3D+1) = rcollection%DquickAccess(1)
    
    ! Set number of spatial dimensions
    ndim = size(Dpoints, 1)

    do iel = 1, size(Ielements)
      do ipoint = 1, ubound(Dpoints,2)
        
        ! Set values for function parser
        Dvalue(1:ndim) = Dpoints(:, ipoint, iel)
        
        ! Evaluate function parser
        call fparser_evalFunction(p_rfparser, icomp, Dvalue,&
            Dvalues(ipoint,iel))
      end do
    end do

  end subroutine transp_weightFuncBdrInt2d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalConvectionP2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij,k_ji,d_ij
!</output>
!</subroutine>

    ! local variables
    real(DP) :: hi,hj,Ei,Ej,ui,uj,vi,vj,ci,cj

    ! Compute convective coefficients
    k_ij = -p_Dvariable1(j)*C_ij(1)-p_Dvariable2(j)*C_ij(2)
    k_ji = -p_Dvariable1(i)*C_ji(1)-p_Dvariable2(i)*C_ji(2)

    ! Set artificial diffusion to zero
    d_ij = 0.0_DP

  end subroutine transp_calcMatGalConvectionP2d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwConvectionP2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij,k_ji,d_ij
!</output>
!</subroutine>

    ! local variables
    real(DP) :: hi,hj,Ei,Ej,ui,uj,vi,vj,ci,cj

    ! Compute convective coefficients
    k_ij = -p_Dvariable1(j)*C_ij(1)-p_Dvariable2(j)*C_ij(2)
    k_ji = -p_Dvariable1(i)*C_ji(1)-p_Dvariable2(i)*C_ji(2)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)

  end subroutine transp_calcMatUpwConvectionP2d

  !*****************************************************************************

!<subroutine>

  pure subroutine transp_calcMatGalConvectionD2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D.
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij,k_ji,d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = p_Dvariable1(j)*C_ij(1)+p_Dvariable2(j)*C_ij(2)
    k_ji = p_Dvariable1(i)*C_ji(1)+p_Dvariable2(i)*C_ji(2)

    ! Set artificial diffusion to zero
    d_ij = 0.0_DP
    
  end subroutine transp_calcMatGalConvectionD2d

  !*****************************************************************************

!<subroutine>

  pure subroutine transp_calcMatUpwConvectionD2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij,k_ji,d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = p_Dvariable1(j)*C_ij(1)+p_Dvariable2(j)*C_ij(2)
    k_ji = p_Dvariable1(i)*C_ji(1)+p_Dvariable2(i)*C_ji(2)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)
    
  end subroutine transp_calcMatUpwConvectionD2d

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrConvectionP2d(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use boundary
    use boundaryfilter
    use collection
    use domainintegration
    use feevaluation
    use fparser
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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
    type(t_vectorBlock), pointer :: p_rvelocity
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dminPar,dmaxPar,dt,dnx,dny,dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint,ndim
    
    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(1)))
    
    ! This subroutine assumes that the first quick access vector
    ! points to the velocity vector
    p_rvelocity => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(ibdrtype)
      
    case(BDR_DIRICHLET_WEAK)
      
      ! Allocate temporal memory
      allocate(Daux(ubound(Dpoints,2), ubound(Dpoints,3), NDIM2D+1))
      
      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1:2)
      call fevl_evaluate_sim1(DER_FUNC2D, Daux(:,:,1),&
          p_rvelocity%RvectorBlock(1), Dpoints, rdomainIntSubset&
          %p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      call fevl_evaluate_sim1(DER_FUNC2D, Daux(:,:,2),&
          p_rvelocity%RvectorBlock(2), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Set number of spatial dimensions
      ndim = size(Dpoints, 1)

      ! Evaluate the function parser for the Dirichlet values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,3).
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          
          ! Set values for function parser
          Dvalue(1:ndim) = Dpoints(:, ipoint, iel)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,3))
        end do
      end do
      
      ! Get the minimum and maximum parameter value. The point with the minimal
      ! parameter value is the start point of the interval, the point with the
      ! maximum parameter value the endpoint.
      dminPar = DpointPar(1,1)
      dmaxPar = DpointPar(1,1)
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          dminPar = min(DpointPar(ipoint,iel), dminPar)
          dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
        end do
      end do
      
      ! Multiply the velocity vector with the normal in each point
      ! to get the normal velocity.
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          
          dt = DpointPar(ipoint,iel)
          
          ! Get the normal vector in the point from the boundary.
          ! Note that the parameter value is in length
          ! parametrisation!  When we are at the left or right
          ! endpoint of the interval, we calculate the normal vector
          ! based on the current edge.  Without that, the behaviour of
          ! the routine may lead to some confusion if the endpoints of
          ! the interval coincide with the endpoints of a boundary
          ! edge. In such a case, the routine would normally compute
          ! the normal vector as a mean on the normal vectors of the
          ! edges adjacent to such a point!
          if (DpointPar(ipoint,iel) .eq. dminPar) then
            ! Start point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
            
          else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
            ! End point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
          else
            ! Inner point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, dt, dnx, dny, cparType=BDR_PAR_LENGTH)
          end if
          
          ! Compute the normal velocity
          dnv = dnx * Daux(ipoint,iel,1) + dny * Daux(ipoint,iel,2)
          
          ! Check if we are at the primal inflow boundary
          if (dnv .lt. -SYS_EPSREAL) then
            Dcoefficients(1,ipoint,iel) = dscale * dnv * Daux(ipoint,iel,3)
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do
      
      ! Deallocate temporal memory
      deallocate(Daux)
      

    case(BDR_INHOMNEUMANN_WEAK)
      
      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Set number of spatial dimensions
      ndim = size(Dpoints, 1)
      
      ! Evaluate the function parser for the Neumann values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,1).
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          
          ! Set values for function parser
          Dvalue(1:ndim) = Dpoints(:, ipoint, iel)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Dcoefficients(1,ipoint,iel))

          ! Multiply by scaling coefficient
          Dcoefficients(1,ipoint,iel) = dscale * Dcoefficients(1,ipoint,iel)
        end do
      end do
      
    end select

  end subroutine transp_coeffVecBdrConvectionP2d

   ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrConvectionD2d(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use boundary
    use boundaryfilter
    use collection
    use domainintegration
    use feevaluation
    use fparser
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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
    ! This routine handles the constant velocities in the dual problem.
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
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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
    type(t_vectorBlock), pointer :: p_rvelocity
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dminPar,dmaxPar,dt,dnx,dny,dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint,ndim
    
    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity vector
    p_rvelocity => rcollection%p_rvectorQuickAccess1
    
    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(ibdrtype)
      
    case(BDR_DIRICHLET_WEAK)

      ! Allocate temporal memory
      allocate(Daux(ubound(Dpoints,2), ubound(Dpoints,3), NDIM2D+1))
      
      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1:2)
      call fevl_evaluate_sim1(DER_FUNC2D, Daux(:,:,1),&
          p_rvelocity%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      call fevl_evaluate_sim1(DER_FUNC2D, Daux(:,:,2),&
          p_rvelocity%RvectorBlock(2), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Set number of spatial dimensions
      ndim = size(Dpoints, 1)
      
      ! Evaluate the function parser for the Dirichlet values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,3).
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          
          ! Set values for function parser
          Dvalue(1:ndim) = Dpoints(:, ipoint, iel)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Daux(ipoint,iel,3))
        end do
      end do
      
      ! Get the minimum and maximum parameter value. The point with the minimal
      ! parameter value is the start point of the interval, the point with the
      ! maximum parameter value the endpoint.
      dminPar = DpointPar(1,1)
      dmaxPar = DpointPar(1,1)
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          dminPar = min(DpointPar(ipoint,iel), dminPar)
          dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
        end do
      end do
      
      ! Multiply the velocity vector with the normal in each point
      ! to get the normal velocity.
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          
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
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
            
          else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
            ! End point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
          else
            ! Inner point
            call boundary_getNormalVec2D(rdiscretisation%p_rboundary,&
                ibct, dt, dnx, dny, cparType=BDR_PAR_LENGTH)
          end if
          
          ! Compute the normal velocity
          dnv = dnx * Daux(ipoint,iel,1) + dny * Daux(ipoint,iel,2)
          
          ! Check if we are at the dual inflow boundary
          if (dnv .gt. SYS_EPSREAL) then
            Dcoefficients(1,ipoint,iel) = dscale * dnv * Daux(ipoint,iel,3)
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do
      
      ! Deallocate temporal memory
      deallocate(Daux)

      case(BDR_INHOMNEUMANN_WEAK)
      
      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Set number of spatial dimensions
      ndim = size(Dpoints, 1)
      
      ! Evaluate the function parser for the Neumann values in the
      ! cubature points on the boundary and store the result in
      ! Dcoefficients(:,:,1).
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          
          ! Set values for function parser
          Dvalue(1:ndim) = Dpoints(:, ipoint, iel)
          
          ! Evaluate function parser
          call fparser_evalFunction(p_rfparser, isegment,&
              Dvalue, Dcoefficients(1,ipoint,iel))

          ! Multiply by scaling coefficient
          Dcoefficients(1,ipoint,iel) = dscale * Dcoefficients(1,ipoint,iel)
        end do
      end do

    end select

  end subroutine transp_coeffVecBdrConvectionD2d  

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrConvectionP2d(rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use boundary
    use boundaryfilter
    use collection
    use domainintegration
    use feevaluation
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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
    type(t_vectorBlock), pointer :: p_rvelocity
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP) :: dminPar,dmaxPar,dt,dnx,dny,dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint,ndim

    ! The first quick access vector points to the velocity vector
    p_rvelocity => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(ibdrtype)

    case(BDR_DIRICHLET_WEAK)

      ! Allocate temporal memory
      allocate(Daux(ubound(Dpoints,2), ubound(Dpoints,3), NDIM2D+1))
      
      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1:2)
      call fevl_evaluate_sim1(DER_FUNC, Daux(:,:,1),&
          p_rvelocity%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      call fevl_evaluate_sim1(DER_FUNC, Daux(:,:,2),&
          p_rvelocity%RvectorBlock(2), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the minimum and maximum parameter value. The point with the minimal
      ! parameter value is the start point of the interval, the point with the
      ! maximum parameter value the endpoint.
      dminPar = DpointPar(1,1)
      dmaxPar = DpointPar(1,1)
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          dminPar = min(DpointPar(ipoint,iel), dminPar)
          dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
        end do
      end do

      ! Multiply the velocity vector with the normal in each point
      ! to get the normal velocity.
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          
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
            call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary,&
                ibct, dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
            
          else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
            ! End point
            call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary,&
                ibct, dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
          else
            ! Inner point
            call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary,&
                ibct, dt, dnx, dny, cparType=BDR_PAR_LENGTH)
          end if
          
          ! Compute the normal velocity
          dnv = dnx * Daux(ipoint,iel,1) + dny * Daux(ipoint,iel,2)
          
          ! Check if we are at the primal inflow boundary
          if (dnv .lt. -SYS_EPSREAL) then
            Dcoefficients(1,ipoint,iel) = dscale * dnv
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux)

    case default

      ! Do nothing

    end select
    
  end subroutine transp_coeffMatBdrConvectionP2d

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrConvectionD2d(rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use boundary
    use boundaryfilter
    use collection
    use domainintegration
    use feevaluation
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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
    ! This routine handles the constant velocities in the dual problem.
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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
    type(t_vectorBlock), pointer :: p_rvelocity
    real(DP), dimension(:,:,:), pointer :: Daux
    real(DP) :: dminPar,dmaxPar,dt,dnx,dny,dnv,dtime,dscale
    integer :: ibdrtype,isegment,iel,ipoint,ndim

    ! This subroutine assumes that the first quick access vector
    ! points to the velocity vector
    p_rvelocity => rcollection%p_rvectorQuickAccess1

    ! The first two quick access double values hold the simulation
    ! time and the scaling parameter
    dtime  = rcollection%DquickAccess(1)
    dscale = rcollection%DquickAccess(2)

    ! The first two quick access integer values hold the type of
    ! boundary condition and the segment number
    ibdrtype = rcollection%IquickAccess(1)
    isegment = rcollection%IquickAccess(2)

    ! What type of boundary conditions are we?
    select case(ibdrtype)

    case(BDR_DIRICHLET_WEAK)
      
      ! Allocate temporal memory
      allocate(Daux(ubound(Dpoints,2), ubound(Dpoints,3), NDIM2D+1))
      
      ! Evaluate the velocity field in the cubature points on the boundary
      ! and store the result in Daux(:,:,:,1:2)
      call fevl_evaluate_sim1(DER_FUNC, Daux(:,:,1),&
          p_rvelocity%RvectorBlock(1), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      call fevl_evaluate_sim1(DER_FUNC, Daux(:,:,2),&
          p_rvelocity%RvectorBlock(2), Dpoints,&
          rdomainIntSubset%p_Ielements, rdomainIntSubset%p_DcubPtsRef)
      
      ! Get the minimum and maximum parameter value. The point with the minimal
      ! parameter value is the start point of the interval, the point with the
      ! maximum parameter value the endpoint.
      dminPar = DpointPar(1,1)
      dmaxPar = DpointPar(1,1)
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          dminPar = min(DpointPar(ipoint,iel), dminPar)
          dmaxPar = max(DpointPar(ipoint,iel), dmaxPar)
        end do
      end do
      
      ! Multiply the velocity vector with the normal in each point
      ! to get the normal velocity.
      do iel = 1, size(rdomainIntSubset%p_Ielements)
        do ipoint = 1, ubound(Dpoints,2)
          
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
            call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary,&
                ibct, dt, dnx, dny, BDR_NORMAL_RIGHT, BDR_PAR_LENGTH)
            
          else if (DpointPar(ipoint,iel) .eq. dmaxPar) then
            ! End point
            call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary,&
                ibct, dt, dnx, dny, BDR_NORMAL_LEFT, BDR_PAR_LENGTH)
          else
            ! Inner point
            call boundary_getNormalVec2D(rdiscretisationTrial%p_rboundary,&
                ibct, dt, dnx, dny, cparType=BDR_PAR_LENGTH)
          end if
          
          ! Compute the normal velocity
          dnv = dnx * Daux(ipoint,iel,1) + dny * Daux(ipoint,iel,2)
          
          ! Check if we are at the dual inflow boundary
          if (dnv .gt. SYS_EPSREAL) then
            Dcoefficients(1,ipoint,iel) = dscale * dnv
          else
            Dcoefficients(1,ipoint,iel) = 0.0_DP
          end if
        end do
      end do
      
      ! Free temporal memory
      deallocate(Daux)

    case default
      
      ! Do nothing

    end select

  end subroutine transp_coeffMatBdrConvectionD2d

  !*****************************************************************************
    
!<subroutine>

  pure subroutine transp_calcMatGalSTBurgersP2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the 
    ! one-dimensional Burgers equation $du/dt+df(u)/dx=0$, whereby
    ! the flux function is given by $f(u)=0.5*u^2$.
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = -0.5_DP*(u_i+u_j)*C_ij(1)-C_ij(2)
    k_ji = -0.5_DP*(u_i+u_j)*C_ji(1)-C_ji(2)

    ! Set artificial diffusion to zero
    d_ij = 0.0_DP

  end subroutine transp_calcMatGalSTBurgersP2d

   !*****************************************************************************
    
!<subroutine>

  pure subroutine transp_calcMatUpwSTBurgersP2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the 
    ! one-dimensional Burgers equation $du/dt+df(u)/dx=0$, whereby
    ! the flux function is given by $f(u)=0.5*u^2$.
    ! Moreover, scalar artificial diffusion is applied
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = -0.5_DP*(u_i+u_j)*C_ij(1)-C_ij(2)
    k_ji = -0.5_DP*(u_i+u_j)*C_ji(1)-C_ji(2)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)

  end subroutine transp_calcMatUpwSTBurgersP2d

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrSTBurgersP2d(rdiscretisation,&
      rform, nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use boundary
    use collection
    use domainintegration
    use feevaluation
    use fparser
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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

    print *, "Weak boundary conditions are not available yet"
    stop

  end subroutine transp_coeffVecBdrSTBurgersP2d

  !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrSTBurgersP2d(&
      rdiscretisationTrial, rdiscretisationTest, rform, nelements,&
      npointsPerElement, Dpoints, ibct, DpointPar, IdofsTrial,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use boundary
    use boundaryfilter
    use collection
    use domainintegration
    use feevaluation
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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

    print *, "Weak boundary conditions are not available yet"
    stop

  end subroutine transp_coeffMatBdrSTBurgersP2d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatGalSTBuckLevP2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the 
    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
    !
    ! Here, the characteristic velocity $a(u)=f^\prime(u)$ is given
    ! by $a(u)=\frac{4u(1-u)}{(3u^2-2u+1)^2}$.
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! local variables
    real(DP) :: v_i,v_j
    
    ! Compute velocities
    v_i = 4*u_i*(1-u_i)/(3*u_i*u_i-2*u_i+1)**2
    v_j = 4*u_j*(1-u_j)/(3*u_j*u_j-2*u_j+1)**2

    ! Compute convective coefficients
    k_ij = -v_j*C_ij(1)-C_ij(2)
    k_ji = -v_i*C_ji(1)-C_ji(2)

    ! Set artificial diffusion to zero
    d_ij = 0.0_DP
        
  end subroutine transp_calcMatGalSTBuckLevP2d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine transp_calcMatUpwSTBuckLevP2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the 
    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
    !
    ! Here, the characteristic velocity $a(u)=f^\prime(u)$ is given
    ! by $a(u)=\frac{4u(1-u)}{(3u^2-2u+1)^2}$.
    ! Moreover, scalar artificial diffusion is applied.
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! local variables
    real(DP) :: v_i,v_j
    
    ! Compute velocities
    v_i = 4*u_i*(1-u_i)/(3*u_i*u_i-2*u_i+1)**2
    v_j = 4*u_j*(1-u_j)/(3*u_j*u_j-2*u_j+1)**2

    ! Compute convective coefficients
    k_ij = -v_j*C_ij(1)-C_ij(2)
    k_ji = -v_i*C_ji(1)-C_ji(2)

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)
        
  end subroutine transp_calcMatUpwSTBuckLevP2d

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrSTBuckLevP2d(rdiscretisation,&
      rform, nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use boundary
    use collection
    use domainintegration
    use feevaluation
    use fparser
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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

    print *, "Weak boundary conditions are not available yet"
    stop

  end subroutine transp_coeffVecBdrSTBuckLevP2d

    !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrSTBuckLevP2d(&
      rdiscretisationTrial, rdiscretisationTest, rform, nelements,&
      npointsPerElement, Dpoints, ibct, DpointPar, IdofsTrial,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use boundary
    use boundaryfilter
    use collection
    use domainintegration
    use feevaluation
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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

    print *, "Weak boundary conditions are not available yet"
    stop

  end subroutine transp_coeffMatBdrSTBuckLevP2d

  !*****************************************************************************
    
!<subroutine>

  pure subroutine transp_calcMatGalBurgersP2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for Burgers` equation in 2D.
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = -0.5_DP*(u_i+u_j)*(C_ij(1)+C_ij(2))
    k_ji = -0.5_DP*(u_i+u_j)*(C_ji(1)+C_ji(2))

    ! Set artificial diffusion to zero
    d_ij = 0.0_DP
    
  end subroutine transp_calcMatGalBurgersP2d

  !*****************************************************************************
    
!<subroutine>

  pure subroutine transp_calcMatUpwBurgersP2d(&
      u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji, d_ij)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for Burgers` equation in 2D.
    ! Moreover, scalar artificial diffusion is applied.
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(in) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij, C_ji

    ! nodal indices
    integer, intent(in) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(out) :: k_ij, k_ji, d_ij
!</output>
!</subroutine>

    ! Compute convective coefficients
    k_ij = -0.5_DP*(u_i+u_j)*(C_ij(1)+C_ij(2))
    k_ji = -0.5_DP*(u_i+u_j)*(C_ji(1)+C_ji(2))

    ! Compute artificial diffusion coefficient
    d_ij = max(-k_ij, 0.0_DP, -k_ji)
    
  end subroutine transp_calcMatUpwBurgersP2d

   ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVecBdrBurgersP2d(rdiscretisation,&
      rform, nelements, npointsPerElement, Dpoints, ibct, DpointPar,&
      IdofsTest, rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use boundary
    use collection
    use domainintegration
    use feevaluation
    use fparser
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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

    print *, "Weak boundary conditions are not available yet"
    stop

  end subroutine transp_coeffVecBdrBurgersP2d

    !*****************************************************************************

!<subroutine>

  subroutine transp_coeffMatBdrBurgersP2d( rdiscretisationTrial,&
      rdiscretisationTest, rform, nelements, npointsPerElement,&
      Dpoints, ibct, DpointPar, IdofsTrial, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use boundary
    use boundaryfilter
    use collection
    use domainintegration
    use feevaluation
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation
    
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

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
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

    print *, "Weak boundary conditions are not available yet"
    stop

  end subroutine transp_coeffMatBdrBurgersP2d
    
end module transport_callback2d
