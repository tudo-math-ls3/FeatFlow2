!##############################################################################
!# Tutorial 022d: Calculate an L2 and H1 error, standard integration, subdomain
!##############################################################################

module tutorial022d

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  
  use triangulation
  use meshgeneration
  
  use element
  use cubature
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation
  use derivatives
  use collection
  
  use derivatives
  use domainintegration
  use collection
  use pprocerror

  implicit none
  private
  
  public :: start_tutorial022d

contains

!****************************************************************************

!<subroutine>

  subroutine fcalc_reffunction (cderivative, rdiscretisation, &
      nelements, npointsPerElement, Dpoints, IdofsTest, rdomainIntSubset, &
      Dvalues, rcollection)

!<description>
  ! Calculates the values of the reference function in the cubature points.
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

  ! Number of elements, where the coefficients must be computed.
  integer, intent(in) :: nelements

  ! Number of points per element, where the coefficients must be computed
  integer, intent(in) :: npointsPerElement

  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: ipt, iel
    real(DP) :: dx, dy
    
    ! What do we have to calculate? Function values? Derivatives?
    select case (cderivative)
    case (DER_FUNC)
      ! Function value
      !
      ! Loop over all elements and points on the elements
      do iel = 1,nelements
        do ipt = 1,npointsPerElement
          
          ! Coordinate of the point
          dx = Dpoints(1,ipt,iel)
          dy = Dpoints(2,ipt,iel)
        
          ! Calculate the value in the point
          Dvalues(ipt,iel) = 1.0_DP/exp(dx)

        end do
      end do
    
    case (DER_DERIV2D_X)
      ! X-derivative
      !
      ! Loop over all elements and points on the elements
      do iel = 1,nelements
        do ipt = 1,npointsPerElement
          
          ! Coordinate of the point
          dx = Dpoints(1,ipt,iel)
          dy = Dpoints(2,ipt,iel)
        
          ! Calculate the value in the point
          Dvalues(ipt,iel) = -1.0_DP/exp(dx)

        end do
      end do
      
    case (DER_DERIV2D_Y)
      ! Y-derivative
      !      
      ! Loop over all elements and points on the elements
      do iel = 1,nelements
        do ipt = 1,npointsPerElement
          
          ! Coordinate of the point
          dx = Dpoints(1,ipt,iel)
          dy = Dpoints(2,ipt,iel)
        
          ! Calculate the value in the point
          Dvalues(ipt,iel) = 0.0_DP

        end do
      end do

    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fcharFunctionSubdomain (rdiscretisation, nelements,&
      npointsPerElement, Dpoints, IdofsTest, rdomainIntSubset,  Dvalues, rcollection)

!<description>
  ! Defines the characteristic function of the domain where to calculate.
  ! This is multiplied to the error/norm in the cubature point.
!</description>

!<input>
  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! Number of elements, where the coefficients must be computed.
  integer, intent(in) :: nelements

  ! Number of points per element, where the coefficients must be computed
  integer, intent(in) :: npointsPerElement

  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in) :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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
  ! This array has to receive the values of the weights in all the
  ! points specified in Dpoints, or the appropriate derivative of
  ! the function, respectively, according to cderivative.
  ! DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: iel,ipt
    real(DP) :: dx,dy

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do ipt = 1,npointsPerElement

        ! Get the coordinates of the point
        dx = Dpoints(1,ipt,iel)
        dy = Dpoints(2,ipt,iel)

        ! Characteristic function of a part of the domain.
        if (dx .le. 0.25_DP) then
          Dvalues(ipt,iel) = 1.0_DP
        else
          Dvalues(ipt,iel) = 0.0_DP
        end if
        
      end do ! icubp
    
    end do ! iel

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial022d

    ! Declare some variables.
    type(t_triangulation) :: rtriangulation
    type(t_spatialDiscretisation) :: rspatialDiscr
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_scalarCubatureInfo), target :: rcubatureInfo
    type(t_vectorBlock) :: rx

    integer :: ivt
    real(DP) :: dx, dintvalue
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 022d")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Create a brick mesh
    ! =================================

    ! The mesh must always be in "standard" format. 
    ! First create a 9x9-mesh on [0,1]x[0,1], then convert to standard.
    call meshgen_rectangular2DQuadMesh (rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 8, 8)
    call tria_initStandardMeshFromRaw (rtriangulation)

    ! =================================
    ! Discretise with Q1.
    !
    ! Create a structure rspatialDiscr
    ! which describes the discretisation.
    ! We generate a 1x1 system with Q1.
    ! =================================

    ! Create a spatial discretisation with Q1
    call spdiscr_initDiscr_simple (rspatialDiscr,EL_Q1_2D,rtriangulation)
    
    ! Create a block discretisation with 1 block Q1.
    call spdiscr_initBlockDiscr (rblockDiscr,rtriangulation)
    call spdiscr_appendBlockComponent (rblockDiscr,rspatialDiscr)
    call spdiscr_commitBlockDiscr (rblockDiscr)

    ! =================================
    ! Create a scalar vector.
    ! =================================

    ! Create a vector.
    call lsysbl_createVector (rblockDiscr,rx)
    
    ! =================================
    ! Fill the vector with data.
    ! =================================
    
    ! Get a pointer to the data.
    call lsyssc_getbase_double (rx%RvectorBlock(1),p_Ddata)
    
    ! Get a pointer to the point coordinates.
    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    ! Set the entries of the vector according to the function
    !    u(x,y) = 1/exp(x)
    do ivt=1,rx%NEQ
      dx = p_DvertexCoords(1,ivt)
      p_Ddata(ivt) = 1.0_DP / exp(dx)
    end do

    ! =================================
    ! Define cubature formula
    ! =================================

    ! Use a Gauss 3x3 formula for the discretisation.
    call spdiscr_createDefCubStructure (rspatialDiscr,rcubatureInfo,CUB_GEN_AUTO_G3)

    ! =================================
    ! Calculate the L2-error of u
    ! =================================

    ! Calculate the error.    
    call pperr_scalar (PPERR_L2ERROR, dintvalue, rx%RvectorBlock(1), &
        fcalc_reffunction,ffunctionWeight=fcharFunctionSubdomain, &
        rcubatureInfo=rcubatureInfo)
        
    call output_line ("L2-error = "//trim(sys_sdEL(dintvalue,10)))

    ! =================================
    ! Calculate the H1-(semi)error of u
    ! =================================
    
    ! Calculate the error.    
    call pperr_scalar (PPERR_H1ERROR, dintvalue, rx%RvectorBlock(1), &
        fcalc_reffunction,ffunctionWeight=fcharFunctionSubdomain, &
        rcubatureInfo=rcubatureInfo)
        
    call output_line ("H1-error = "//trim(sys_sdEL(dintvalue,10)))

    ! =================================
    ! Cleanup
    ! =================================
    
    ! Release the vector
    call lsysbl_releaseVector (rx)
    
    ! Cubature done.
    call spdiscr_releaseCubStructure (rcubatureInfo)

    ! Release the discretisation
    call spdiscr_releaseBlockDiscr (rblockDiscr)
    call spdiscr_releaseDiscr (rspatialDiscr)

    ! Release the triangulation
    call tria_done (rtriangulation)
    
  end subroutine

end module
