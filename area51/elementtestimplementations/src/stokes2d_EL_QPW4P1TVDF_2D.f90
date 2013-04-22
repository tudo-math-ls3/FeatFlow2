!##############################################################################
!# ****************************************************************************
!# <name> stokes2d_EL_QPW4P1TVDF_2D </name>
!# ****************************************************************************
!#
!# <purpose>
!# Tests the element EL_QPW4P1TVDF_2D. 
!# </purpose>
!##############################################################################

module stokes2d_EL_QPW4P1TVDF_2D

  use fsystem
  use storage
  use genoutput
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use linearalgebra
  use discretebc
  use bcassembly
  use triangulation
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use spdiscprojection
  use scalarpde
  use bilinearformevaluation
  use linearformevaluation
  use discretebc
  use filtersupport
  use coarsegridcorrection
  use linearsolver
  use ucd
  use matrixio
  use mprimitives
  use matrixmodification
  
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop
  use collection
  
  implicit none

contains
  
  !****************************************************************************

!<subroutine>

  subroutine st2d3_fcalc_Stokes(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Stokes operator.
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!<subroutine>

    real(DP) :: dnu
    type(t_collection) :: rcoll  
  
    ! Viscosity
    dnu = rcollection%DquickAccess(1)
    
    ! Laplace operator
    rcoll%DquickAccess(1) = dnu
    rcoll%IquickAccess(1) = 1
    rcoll%IquickAccess(2) = 1
    call bma_fcalc_laplace(RmatrixData,rassemblyData,rmatrixAssembly,&
        npointsPerElement,nelements,revalVectors,rcoll)

    ! Gradient/divergence operators.
    rcoll%DquickAccess(1) = 1.0_DP
    rcoll%IquickAccess(1) = 2
    rcoll%IquickAccess(2) = 1
    call bma_fcalc_gradientdiv(RmatrixData,rassemblyData,rmatrixAssembly,&
        npointsPerElement,nelements,revalVectors,rcoll)

    rcoll%DquickAccess(1) = 1.0_DP
    rcoll%IquickAccess(1) = 1
    rcoll%IquickAccess(2) = 2
    call bma_fcalc_divergence(RmatrixData,rassemblyData,rmatrixAssembly,&
        npointsPerElement,nelements,revalVectors,rcoll)

  end subroutine

  ! ***********************************************************************************************

  elemental real(DP) function funcDelta(dx, dy, dsigma, dh)
  real(DP), intent(in) :: dx, dy, dsigma, dh
  real(DP), parameter :: dr0 = 0.25_DP
  real(DP) :: dr, dn, dkappa, ddist, dw, ddelta
  
    dr = sqrt(dx*dx + dy*dy)
    dkappa = -1._DP / dr
    dn = dx / dr
    ddist = dr - dr0
    dw = ddist / dh
    if(abs(dw) < 1._DP) then
      ddelta = 35._DP/32._DP * (1._DP - 3._DP*dw**2 + 3._DP*dw**4 - dw**6) / dh
      funcDelta = ddelta*dkappa*dn*dsigma
    else
      funcDelta = 0._DP
    end if

  end function

  !****************************************************************************

!<subroutine>

  subroutine st2d3_fcalc_rhs(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Mass operator in all diagonal matrices.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: rvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!<subroutine>

    ! Local variables
    real(DP) :: dbasI1,dbasI2, dval1, dval2, dval3, dx,dy
    integer :: iel, icubp, idofe,ndimfe
    real(DP), dimension(:,:), pointer :: p_DlocalVector1,p_DlocalVector2
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest1,p_DbasTest2
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData1,p_rvectorData2
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    real(DP) :: dsigma,dh
    
    dsigma = rcollection%DquickAccess(1)
    dh = rcollection%DquickAccess(2)
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    p_rvectorData1 => RvectorData(1)
    p_rvectorData2 => RvectorData(2)
    p_DlocalVector1 => RvectorData(1)%p_Dentry
    p_DlocalVector2 => RvectorData(2)%p_Dentry
    p_DbasTest1 => RvectorData(1)%p_DbasTest
    p_DbasTest2 => RvectorData(2)%p_DbasTest
    
    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Calculate the RHS of the momentum equation
    
    ! FE space dimension
    ndimfe = RvectorData(1)%ndimfe

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Values of the velocity RHS for the X1 and X2 component
        !dval1 = funcDelta(dx, dy, dsigma, dh)
        !dval2 = funcDelta(dy, dx, dsigma, dh)
        dval1 = dsigma * 3.0_DP*dx**2
        dval2 = dsigma * 3.0_DP*dy**2
        
        ! Outer loop over the DOFs i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rvectorData1%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          ! Multiply the values of the basis functions
          ! (1st derivatives) by the cubature weight and sum up
          ! into the local vectors.
          dbasI1 = p_DbasTest1(idofe,DER_FUNC,icubp,iel)
          dbasI2 = p_DbasTest1(idofe+p_rvectorData1%ndofTest,DER_FUNC,icubp,iel)
          
          p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
              p_DcubWeight(icubp,iel) * ( dval1 * dbasI1 + dval2 * dbasI2 )
          
        end do ! jdofe
          
      end do ! icubp
    
    end do ! iel

    ! Calculate the RHS of the continuity equation
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Value of the pressure RHS
        dval3 = 0.0_DP

        ! Outer loop over the DOFs i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rvectorData2%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI1 = p_DbasTest2(idofe,DER_FUNC,icubp,iel)
          
          ! Multiply the values of the basis functions
          ! by the cubature weight and sum up into the local vectors.
          p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
              p_DcubWeight(icubp,iel) * dval3 * dbasI1
          
        end do ! jdofe

      end do ! icubp
    
    end do ! iel
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine st2d3_getBoundaryValues_2D (Icomponents,rdiscretisation,rboundaryRegion, &
      ielement,cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
  use collection
  use spatialdiscretisation
  use discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! "snapshot" of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! Component specifier.
  ! For Dirichlet boundary:
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry,
  !   2=2nd solution component, e.g. Y-velocity,...)
  integer, dimension(:), intent(in)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Boundary region that is currently being processed.
  type(t_boundaryRegion), intent(in)                          :: rboundaryRegion
  
  ! The element number on the boundary which is currently being processed
  integer, intent(in)                                         :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(in)                                         :: cinfoNeeded
  
  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC :
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDDERIV :
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN :
  !   iwhere = number of the edge where the value integral mean value
  !            should be computed
  integer, intent(in)                                          :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC :
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV :
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN :
  !   dwhere = 0 (not used)
  real(DP), intent(in)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional                 :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>
  
!</subroutine>

    real(DP) :: dx,dy
    
    call boundary_getCoords(rdiscretisation%p_rboundary, &
        rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
    
    ! Return zero Dirichlet boundary values for all situations by default.
    ! For vector-valued elements,
    ! - Dvalues(1) returns the prescribed X-velocity,
    ! - Dvalues(2) returns the prescribed Y-velocity.
    Dvalues(1) = 0.0_DP
    Dvalues(2) = 0.0_DP

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine st2d3_fcalc_L2error(dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the (squared) L2 error.
    ! The component to calculate the error for has to be provided
    ! in rcollection%IquickAccess(1). The viscosity has to be stored
    ! in rcollection%DquickAccess(1).
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dvalx,dvaly,dvalrefx,dvalrefy,dval,dvalref,dx,dy,dnu
    integer :: iel, icubp, icomp
    real(DP) :: dsigma
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
    real(DP), dimension(:,:,:), pointer :: p_DfuncVec
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Calcel if no FEM function is given.
    if (revalVectors%ncount .eq. 0) return

    dsigma = rcollection%DquickAccess(1)

    ! Skip interleaved vectors.
    if (revalVectors%p_RvectorData(1)%bisInterleaved) return

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the component from the collection
    icomp = rcollection%IquickAccess(1)

    ! Get the viscosity parameter
    dnu = rcollection%DquickAccess(1)

    dintvalue = 0.0_DP
    
    ! Which component to calculate?
    select case (icomp)
    
    ! ----------
    ! Velocity
    ! ----------
    case (1)
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_DfuncVec => revalVectors%p_RvectorData(1)%p_DdataVec(:,:,:,DER_FUNC2D)
    
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the bubble function
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          
          ! Reference functions
          dvalrefx = 0.0_DP
          dvalrefy = 0.0_DP

          ! Get the error of the FEM function to the bubble function
          dvalx = p_DfuncVec(1,icubp,iel)
          dvaly = p_DfuncVec(2,icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the (squared) L2 error:
          dintvalue = dintvalue + &
              p_DcubWeight(icubp,iel) * (dvalx - dvalrefx)**2
            
        end do ! icubp
      
      end do ! iel
      
    case (2)
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_DfuncVec => revalVectors%p_RvectorData(1)%p_DdataVec(:,:,:,DER_FUNC2D)
    
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the bubble function
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          
          ! Reference functions
          dvalrefx = 0.0_DP
          dvalrefy = 0.0_DP

          ! Get the error of the FEM function to the bubble function
          dvalx = p_DfuncVec(1,icubp,iel)
          dvaly = p_DfuncVec(2,icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the (squared) L2 error:
          dintvalue = dintvalue + &
              p_DcubWeight(icubp,iel) * (dvaly - dvalrefy)**2
            
        end do ! icubp
      
      end do ! iel

    ! ----------
    ! Pressure
    ! ----------
    case (3)
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_Dfunc => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC2D)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the bubble function
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          
          ! Reference function
          dval = (dx**3+dy**3-0.5_DP) * dsigma

          ! Get the error of the FEM function to the bubble function
          dvalref = p_Dfunc(icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the (squared) L2 error:
          dintvalue = dintvalue + &
              p_DcubWeight(icubp,iel) * (dval - dvalref)**2
            
        end do ! icubp
      
      end do ! iel

    end select
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine st2d3_fcalc_H1error(dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the (squared) H1 error.
    ! The component to calculate the error for has to be provided
    ! in rcollection%IquickAccess(1). The viscosity has to be stored
    ! in rcollection%DquickAccess(1).
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dderivX1,dderivY1,dderivX2,dderivY2,dx,dy
    integer :: iel, icubp, icomp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:,:), pointer :: p_DderivXvec,p_DderivYvec
    real(DP), dimension(:,:), pointer :: p_DderivX,p_DderivY
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Calcel if no FEM function is given.
    if (revalVectors%ncount .eq. 0) return

    ! Skip interleaved vectors.
    if (revalVectors%p_RvectorData(1)%bisInterleaved) return

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the component from the collection
    icomp = rcollection%IquickAccess(1)

    dintvalue = 0.0_DP

    ! Which component to calculate?
    select case (icomp)
    
    ! ----------
    ! X-velocity
    ! ----------
    case (1)
    
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_DderivXvec => revalVectors%p_RvectorData(1)%p_DdataVec(:,:,:,DER_DERIV2D_X)
      p_DderivYvec => revalVectors%p_RvectorData(1)%p_DdataVec(:,:,:,DER_DERIV2D_Y)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the derivatives of the bubble function in the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          
          dderivX1 = 0.0_DP
          dderivY1 = 0.0_DP

          ! Get the error of the FEM function derivatives of the bubble function
          ! in the cubature point
          dderivX2 = p_DderivXvec(1,icubp,iel)
          dderivY2 = p_DderivYvec(1,icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the (squared) H1 error:
          dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
              ( (dderivX1 - dderivX2)**2 + (dderivY1 - dderivY2)**2 )
            
        end do ! icubp
      
      end do ! iel
    
    ! ----------
    ! Y-velocity
    ! ----------
    case (2)
    
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_DderivXvec => revalVectors%p_RvectorData(1)%p_DdataVec(:,:,:,DER_DERIV2D_X)
      p_DderivYvec => revalVectors%p_RvectorData(1)%p_DdataVec(:,:,:,DER_DERIV2D_Y)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the derivatives of the bubble function in the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          
          dderivX1 = 0.0_DP
          dderivY1 = 0.0_DP

          ! Get the error of the FEM function derivatives of the bubble function
          ! in the cubature point
          dderivX2 = p_DderivXvec(2,icubp,iel)
          dderivY2 = p_DderivYvec(2,icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the (squared) H1 error:
          dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
              ( (dderivX1 - dderivX2)**2 + (dderivY1 - dderivY2)**2 )
            
        end do ! icubp
      
      end do ! iel
    
    
    ! ----------
    ! Pressure
    ! ----------
    case (3)

      ! not used here.

      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_DderivX => revalVectors%p_RvectorData(icomp)%p_Ddata(:,:,DER_DERIV2D_X)
      p_DderivY => revalVectors%p_RvectorData(icomp)%p_Ddata(:,:,DER_DERIV2D_Y)

      
    end select
            
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine test_stokes2d_EL_QPW4P1TVDF_2D
  
!<description>
  ! This is an all-in-one stokes solver for directly solving a stokes
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let us see...
    !
    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rvector,rrhs,rtempBlock,rprjVector
    type(t_matrixBlock) :: rmatrixPmass
    
    ! Cubature information structure which defines the cubature formula.
    type(t_scalarCubatureInfo) :: rcubatureInfo,rcubatureInfo2
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(10), target :: RfilterChain
    
    ! Number of filters in the filter chain
    integer :: nfilters
    
    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Collection structure for callback routines
    type(t_collection) :: rcollection
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir
    
    ! Vector evaluation structure for the calculation of errors.
    type(t_fev2Vectors) :: revalVectors

    ! Error arrays for post-processing
    real(DP), dimension(2) :: DerrorUL2, DerrorUH1
    real(DP), dimension(1) :: DerrorPL2
    real(DP) :: derror, ddiv
    real(DP), dimension(:), pointer :: p_DdataX,p_DdataP
    real(DP) :: dh, dsigma
    logical :: bpureDirichlet
    real(DP), dimension(2), parameter :: DsigmaList = (/ 1.0_DP, 1E+3_DP /)
    integer :: isigma

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    ! As we do not use a multigrid solver here, we will set the level to 5.
    do NLMAX = 5,5
    
      do isigma = 1,size(DsigmaList)
      
        ! Parameters defining the pressure function
        dh = 0.05_DP
        dsigma = DsigmaList(isigma)
        
        ! Viscosity parameter:
        dnu = 1.0_DP
        
        call output_lbrk()
        call output_line ("Level " // trim(sys_siL(NLMAX,10)) // ", sigma=" // &
            trim(sys_sdEL(dsigma,2)) )
        call output_lbrk()

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Read the domain, read the mesh, refine
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! Get the path $PREDIR from the environment, where to read .prm/.tri files
        ! from. If that does not exist, write to the directory "./pre".
        if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

        ! At first, read in the parametrisation of the boundary and save
        ! it to rboundary.
        call boundary_read_prm(rboundary, trim(spredir)//"/QUAD.prm")
            
        ! Now read in the basic triangulation.
        call tria_readTriFile2D (rtriangulation, trim(spredir)//"/QUAD.tri", rboundary)
        
        ! Refine the mesh up to the minimum level
        call tria_quickRefine2LevelOrdering(NLMAX-1,rtriangulation,rboundary)
        
        ! Create information about adjacencies and everything one needs from
        ! a triangulation. Afterwards, we have the coarse mesh.
        call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Set up a discretisation structure which tells the code which
        ! finite element to use
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! Now we can start to initialise the discretisation. At first, set up
        ! a block discretisation structure that specifies 3 blocks in the
        ! solution vector.
        call spdiscr_initBlockDiscr (rdiscretisation,2,&
                                    rtriangulation, rboundary)

        ! rdiscretisation%RspatialDiscr is a list of scalar
        ! discretisation structures for every component of the solution vector.
        ! We have a solution vector with three components:
        !  Component 1 = X-velocity
        !  Component 2 = Y-velocity
        !  Component 3 = Pressure
        ! For simplicity, we set up one discretisation structure for the
        ! velocity...
        call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1),&
                    EL_QPW4P1TVDF_2D, rtriangulation, rboundary)
                    
        ! For the pressure (3rd component), we set up a separate discretisation
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
            EL_QPW4P0_2D, rdiscretisation%RspatialDiscr(2))

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Set up an cubature info structure to tell the code which cubature
        ! formula to use
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                     
        ! Create an assembly information structure which tells the code
        ! the cubature formula to use. Standard: Gauss 3x3.
        call spdiscr_createDefCubStructure(&  
            rdiscretisation%RspatialDiscr(1),rcubatureInfo,CUB_QPW4QG5T_2D)

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Create a block matrix structure for the operator
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! Initialise the block matrix with default values based on
        ! the discretisation.
        call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrix)
        
        ! Now as the discretisation is set up, we can start to generate
        ! the structure of the system matrix which is to solve.
        ! We create that directly in the block (1,1) of the block matrix
        ! using the discretisation structure of the first block.
        !
        ! In the global system, there are two coupling matrices B1 and B2.
        ! Both have the same structure.
        !
        !    ( A         B1 )
        !    (      A    B2 )
        !    ( B1^T B2^T    )
        !
        ! Create the matrix structure of the X-velocity.
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                        LSYSSC_MATRIX9, rmatrix%RmatrixBlock(1,1))
        
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(2),&
            LSYSSC_MATRIX9, rmatrix%RmatrixBlock(1,2),rdiscretisation%RspatialDiscr(1))
        
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
            LSYSSC_MATRIX9, rmatrix%RmatrixBlock(2,1),rdiscretisation%RspatialDiscr(2))
        
        call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(2),&
            LSYSSC_MATRIX9, rmatrix%RmatrixBlock(2,2),rdiscretisation%RspatialDiscr(2))
            
        ! Now re-assign the block discretisation structure to all matrices
        call lsysbl_assignDiscrDirectMat (rmatrix,rdiscretisation)
        
        ! Allocate memory for the matrix
        call lsysbl_allocEmptyMatrix (rmatrix,LSYSSC_SETM_ZERO)

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Create a block matrix entries
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call output_line("Generating matrices...")
        
        ! Pass dnu via rcollection
        rcollection%DquickAccess(1) = dnu
        
        call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,&
              st2d3_fcalc_Stokes,rcubatureInfo=rcubatureInfo,rcollection=rcollection)
        
        call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrixPmass)
        call lsyssc_copyMatrix (rmatrix%RmatrixBlock(2,2),rmatrixPmass%RmatrixBlock(2,2))
        call bma_buildMatrix (rmatrixPmass,BMA_CALC_STANDARD,&
              bma_fcalc_massDiag,rcubatureInfo=rcubatureInfo)
        call lsyssc_lumpMatrixScalar (rmatrixPmass%RmatrixBlock(2,2),LSYSSC_LUMP_DIAG,.true.)
        call mmod_expandToFullRow (rmatrix%RmatrixBlock(2,2),1)

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Create RHS and solution vectors
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call output_line("Generating vectors...")

        ! Create a RHS and a solution vector based on the discretisation.
        ! Fill with zero.
        call lsysbl_createVectorBlock (rdiscretisation,rrhs,.true.)
        call lsysbl_createVectorBlock (rdiscretisation,rvector,.true.)

        ! Assemble the right hand side.
        rcollection%DquickAccess(1) = dsigma
        rcollection%DquickAccess(2) = dh
        call bma_buildVector (rrhs,BMA_CALC_STANDARD,&
              st2d3_fcalc_rhs,rcubatureInfo=rcubatureInfo,rcollection=rcollection)
                                    
        ! Debug
        call lsyssc_getbase_double (rrhs%RvectorBlock(1),p_DdataX)
        call lsyssc_getbase_double (rrhs%RvectorBlock(2),p_DdataP)

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Assembly of matrices/vectors finished
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call output_line("Setting up BC...")

        ! Set to TRUE for a full Dirichlet problem!!!
        bpureDirichlet = .true.

        ! For implementing boundary conditions, we use a `filter technique with
        ! discretised boundary conditions`. This means, we first have to calculate
        ! a discrete version of the analytic BC, which we can implement into the
        ! solution/RHS vectors using the corresponding filter.
        !
        ! Create a t_discreteBC structure where we store all discretised boundary
        ! conditions.
        call bcasm_initDiscreteBC(rdiscreteBC)
        
        ! We first set up the boundary conditions for the X-velocity, then those
        ! of the Y-velocity.
        !
        ! We "know" already (from the problem definition) that we have four boundary
        ! segments in the domain. Each of these, we want to use for enforcing
        ! some kind of boundary condition.
        !
        ! We ask the boundary routines to create a "boundary region" - which is
        ! simply a part of the boundary corresponding to a boundary segment.
        ! A boundary region roughly contains the type, the min/max parameter value
        ! and whether the endpoints are inside the region or not.
        call boundary_createRegion(rboundary,1,1,rboundaryRegion)
        
        ! The endpoint of this segment should also be Dirichlet. We set this by
        ! changing the region properties in rboundaryRegion.
        rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
        
        ! We use this boundary region and specify that we want to have Dirichlet
        ! boundary there. The following call does the following:
        ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
        !   We specify icomponent="1" to indicate that we set up the
        !   Dirichlet BC`s for the first (here: one and only) component in the
        !   solution vector.
        ! - Discretise the boundary condition so that the BC`s can be applied
        !   to matrices and vectors
        ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
        call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                          rboundaryRegion,rdiscreteBC,&
                                          st2d3_getBoundaryValues_2D)
                                 
        ! Edge 2 is Neumann boundary, so it is commented out.
        if (bpuredirichlet) then
          CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
          CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                            rboundaryRegion,rdiscreteBC,&
                                            st2d3_getBoundaryValues_2D)
        end if
                                 
        ! Edge 3 of boundary component 1.
        call boundary_createRegion(rboundary,1,3,rboundaryRegion)
        call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                          rboundaryRegion,rdiscreteBC,&
                                          st2d3_getBoundaryValues_2D)
        
        ! Edge 4 of boundary component 1. That is it.
        call boundary_createRegion(rboundary,1,4,rboundaryRegion)
        call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                          rboundaryRegion,rdiscreteBC,&
                                          st2d3_getBoundaryValues_2D)

        ! Edge 1 of boundary component 2. That is it.
        ! call boundary_createRegion(rboundary,2,1,rboundaryRegion)
        ! call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
        !                                    rboundaryRegion,rdiscreteBC,&
        !                                    st2d3_getBoundaryValues_2D)

        ! The pressure does not need boundary conditions.
        ! Next step is to implement boundary conditions into the RHS,
        ! solution and matrix. This is done using a vector/matrix filter
        ! for discrete boundary conditions.
        call vecfil_discreteBCrhs (rrhs,rdiscreteBC)
        call vecfil_discreteBCsol (rvector,rdiscreteBC)
        call matfil_discreteBC (rmatrix,rdiscreteBC)

        ! Special modification of the matrix/vectors in case of a pure
        ! Dirichlet problem. Implementation of the integral mean value
        ! constraint in the pressure.
        if (bpureDirichlet) then
        
          call mmod_replaceLineByLumpedMass (rmatrix%RmatrixBlock(2,2),1,rmatrixPmass%RmatrixBlock(2,2))
          call mmod_replaceLinesByZero (rmatrix%RmatrixBlock(2,1),(/1/))

          call vecfil_oneEntryZero (rrhs,2,1)
          call vecfil_oneEntryZero (rvector,2,1)

        end if

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Set up a linear solver
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! During the linear solver, the boundary conditions are also
        ! frequently imposed to the vectors. But as the linear solver
        ! does not work with the actual solution vectors but with
        ! defect vectors instead.
        ! So, set up a filter chain that filters the defect vector
        ! during the solution process to implement discrete boundary conditions.
        call filter_clearFilterChain (RfilterChain,nfilters)
        call filter_newFilterDiscBCDef (RfilterChain,nfilters,rdiscreteBC)

        if (bpureDirichlet) then
          call filter_newFilterOneEntryZero (RfilterChain,nfilters,2,1)
        end if

        ! Create a BiCGStab-solver with VANCA preconditioner.
        ! Attach the above filter chain to the solver, so that the solver
        ! automatically filters the vector during the solution process.
        nullify(p_rpreconditioner)
        call linsol_initUMFPACK4 (p_rsolverNode)

        ! Set the output level of the solver to 2 for some output
        p_rsolverNode%ioutputLevel = 2

        ! We will allow the solver to perform 200 iterations
        p_rsolverNode%nmaxIterations = 2000

        ! Attach the system matrix to the solver.
        call linsol_setMatrix(p_rsolverNode,rmatrix)
        
        ! Initialise structure/data of the solver. This allows the
        ! solver to allocate memory / perform some precalculation
        ! to the problem.
        call output_line("Symbolic factorisation...")
        call linsol_initStructure (p_rsolverNode, ierror)

        if (ierror .ne. LINSOL_ERR_NOERROR) then
          call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
          call sys_halt()
        end if
        
        call output_line("Numeric factorisation...")
        call linsol_initData (p_rsolverNode, ierror)

        if (ierror .ne. LINSOL_ERR_NOERROR) then
          call output_line("Matrix singular!",OU_CLASS_ERROR)
          call sys_halt()
        end if

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Solve the system
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! Finally solve the system. As we want to solve Ax=b with
        ! b being the real RHS and x being the real solution vector,
        ! we use linsol_solveAdaptively. If b is a defect
        ! RHS and x a defect update to be added to a solution vector,
        ! we would have to use linsol_precondDefect instead.
        call output_line("Solving...")
        call linsol_solveAdaptively (p_rsolverNode,rvector,rrhs,rtempBlock)
        
        ! Debug
        call lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataX)
        call lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataP)

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Postprocessing of the solution
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call output_line("Postprocessing...")
        
        ! That is it, rvecSol now contains our solution. We can now
        ! start the postprocessing.

        ! Project the solution to the vertices.
        nullify(p_DdataX)
        nullify(p_DdataP)
        call spdp_projectToVertices(rvector%RvectorBlock(1),p_DdataX,DER_FUNC)
        call spdp_projectToCells(rvector%RvectorBlock(2),p_DdataP,DER_FUNC)
        
        ! Get the path for writing postprocessing files from the environment variable
        ! $UCDDIR. If that does not exist, write to the directory "./gmv".
        if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

        ! Start UCD export to VTK file:
        call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                          trim(sucddir)//"/u2d_0_Q1VDF.vtk")
        
        ! Add the solution to the UCD exporter.
        ! Fir vector-valued elements, the first NVT entries contain the X-velocity, 
        ! the next NVT entries the Y-velocity.
        call ucd_addVarVertBasedVec(rexport,"velocity",&
            p_DdataX(1:rtriangulation%NVT),p_DdataX(rtriangulation%NVT+1:2*rtriangulation%NVT))
        call ucd_addVariableElementBased (rexport, "solp", UCD_VAR_STANDARD, p_DdataP)

        ! Write the file to disc, that is it.
        call ucd_write (rexport)
        call ucd_release (rexport)
        
        ! Release temporary storage.
        deallocate(p_DdataX)
        deallocate(p_DdataP)

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Projection and VTK export finished.
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Error calculation
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
        ! Store the viscosity parameter nu in the collection"s quick access array
        rcollection%DquickAccess(1) = dsigma
        rcollection%DquickAccess(2) = dh
        
        ! Set up revalVectors with the velocity/pressure vectors.
        call fev2_addVectorToEvalList(revalVectors,rvector%RvectorBlock(1),1)
        call fev2_addVectorToEvalList(revalVectors,rvector%RvectorBlock(2),0)

        ! Calculate errors of velocity and pressure against analytic solutions.
        rcollection%IquickAccess(1) = 1  ! X-component, L2-error
        call bma_buildIntegral(DerrorUL2(1),BMA_CALC_STANDARD,st2d3_fcalc_L2error,&
            rcollection=rcollection,revalVectors=revalVectors,&
            rcubatureInfo=rcubatureInfo)
            
        rcollection%IquickAccess(1) = 2  ! X-component, L2-error
        call bma_buildIntegral(DerrorUL2(2),BMA_CALC_STANDARD,st2d3_fcalc_L2error,&
            rcollection=rcollection,revalVectors=revalVectors,&
            rcubatureInfo=rcubatureInfo)

        rcollection%IquickAccess(1) = 3  ! pressure, l2-error
        call bma_buildintegral(derrorPL2(1),BMA_CALC_STANDARD,st2d3_fcalc_l2error,&
            rcollection=rcollection,revalvectors=revalvectors,&
            rcubatureinfo=rcubatureinfo)
    
        rcollection%iquickaccess(1) = 1  ! x-component, h1-error
        call bma_buildintegral(derrorUH1(1),BMA_CALC_STANDARD,st2d3_fcalc_h1error,&
            rcollection=rcollection,revalvectors=revalvectors,&
            rcubatureinfo=rcubatureinfo)
            
        rcollection%iquickaccess(1) = 2  ! y-component, h1-error
        call bma_buildintegral(derrorUH1(2),BMA_CALC_STANDARD,st2d3_fcalc_h1error,&
            rcollection=rcollection,revalvectors=revalvectors,&
            rcubatureinfo=rcubatureInfo)

        call bma_buildintegral(ddiv,BMA_CALC_STANDARD,bma_fcalc_divergenceL2norm,&
            rcollection=rcollection,revalvectors=revalvectors,&
            rcubatureinfo=rcubatureInfo)

        ! Print the errors.
        call output_lbrk()
        call output_line("|u - u_h|_L2 = " // trim(sys_sdEL(sqrt(DerrorUL2(1)), 2)) &
                                    // " " // trim(sys_sdEL(sqrt(DerrorUL2(2)), 2)) &
                                    // " => " // &
                                    trim(sys_sdEL(sqrt(DerrorUL2(1) + DerrorUL2(2)), 2)))
        call output_line("|u - u_h|_H1 = " // trim(sys_sdEL(sqrt(DerrorUH1(1)), 2)) &
                                    // " " // trim(sys_sdEL(sqrt(DerrorUH1(2)), 2)) &
                                    // " => " // &
                                    trim(sys_sdEL(sqrt(DerrorUH1(1) + DerrorUH1(2)), 2)))
        call output_line("|p - p_h|_L2 =                      " // &
            trim(sys_sdEL(sqrt(derrorPL2(1)), 2)))
            
        call output_lbrk()
        call output_line ("|err|_L2(u) / H1(u) / L2(p) / div(u) = " // &
            trim(sys_sdEL(sqrt(DerrorUL2(1) + DerrorUL2(2)), 2)) // " " // &
            trim(sys_sdEL(sqrt(DerrorUH1(1) + DerrorUH1(2)), 2)) // " " // &
            trim(sys_sdEL(sqrt(derrorPL2(1)), 2)) // " " // &
            trim(sys_sdEL(sqrt(ddiv), 2)))
        
        ! Cleanup
        call fev2_releaseVectorList(revalVectors)
        
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Clean up
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! We are finished - but not completely!
        ! Now, clean up so that all the memory is available again.
        !
        ! Release solver data and structure
        call linsol_doneData (p_rsolverNode)
        call linsol_doneStructure (p_rsolverNode)
        
        ! Release the solver node and all subnodes attached to it (if at all):
        call linsol_releaseSolver (p_rsolverNode)
        
        ! Release the block matrix/vectors
        call lsysbl_releaseVector (rtempBlock)
        call lsysbl_releaseVector (rvector)
        call lsysbl_releaseVector (rrhs)
        call lsysbl_releaseMatrix (rmatrix)
        call lsysbl_releaseMatrix (rmatrixPmass)
        
        ! Release our discrete version of the boundary conditions
        call bcasm_releaseDiscreteBC (rdiscreteBC)

        ! Release the discretisation structure and all spatial discretisation
        ! structures in it.
        call spdiscr_releaseBlockDiscr(rdiscretisation)
        
        ! Release the triangulation.
        call tria_done (rtriangulation)
        
        ! Finally release the domain, that is it.
        call boundary_release (rboundary)
      
      end do
      
    end do

  end subroutine

end module
