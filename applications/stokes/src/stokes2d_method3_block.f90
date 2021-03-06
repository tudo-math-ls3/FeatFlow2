!##############################################################################
!# ****************************************************************************
!# <name> stokes2d_method3_block </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a Stokes
!# problem on a simple domain.
!#
!# The routine uses the BiCGStab solver with a simple-VANKA preconditioner
!# for 2D saddle point problems, Jacobi-Type.
!# The discretisation uses the block assembly method to evaluate the
!# mattrix all-in-one.
!# </purpose>
!##############################################################################

module stokes2d_method3_block

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
  
  use stokes2d_callback
  
  use blockmatassemblybase
  use blockmatassembly
  use feevaluation2
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

    ! Local variables
    real(DP) :: dbasJ, dval, dbasIx, dbasIy, dbasJx, dbasJy, dnu
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA11,p_DlocalMatrixA22
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA13,p_DlocalMatrixA23
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrixA31,p_DlocalMatrixA32
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA11,p_DbasTestA11
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrialA13,p_DbasTestA13
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixDataA11,p_rmatrixDataA22,p_rmatrixDataA13
  
    ! Viscosity
    dnu = rcollection%DquickAccess(1)
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    p_rmatrixDataA11 => RmatrixData(1,1)
    p_rmatrixDataA22 => RmatrixData(2,2)
    p_rmatrixDataA13 => RmatrixData(1,3)
    p_DlocalMatrixA11 => RmatrixData(1,1)%p_Dentry
    p_DlocalMatrixA22 => RmatrixData(2,2)%p_Dentry
    p_DlocalMatrixA13 => RmatrixData(1,3)%p_Dentry
    p_DlocalMatrixA23 => RmatrixData(2,3)%p_Dentry
    p_DlocalMatrixA31 => RmatrixData(3,1)%p_Dentry
    p_DlocalMatrixA32 => RmatrixData(3,2)%p_Dentry
    p_DbasTrialA11 => RmatrixData(1,1)%p_DbasTrial
    p_DbasTestA11 => RmatrixData(1,1)%p_DbasTest
    p_DbasTrialA13 => RmatrixData(1,3)%p_DbasTrial
    p_DbasTestA13 => RmatrixData(1,3)%p_DbasTest
    
    ! Calculate Laplace in the diagonal blocks
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Outer loop over the DOFs i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rmatrixDataA11%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasIx = p_DbasTestA11(idofe,DER_DERIV2D_X,icubp,iel)
          dbasIy = p_DbasTestA11(idofe,DER_DERIV2D_Y,icubp,iel)
          
          ! Inner loop over the DOFs j=1..ndof, which corresponds to
          ! the basis function Phi_j:
          do jdofe=1,p_rmatrixDataA11%ndofTrial
            
            ! Fetch the contributions of the (trial) basis function Phi_j
            ! into dbasJ
            dbasJx = p_DbasTrialA11(jdofe,DER_DERIV2D_X,icubp,iel)
            dbasJy = p_DbasTrialA11(jdofe,DER_DERIV2D_Y,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local matrices.
            dval = dnu * p_DcubWeight(icubp,iel) * ( dbasJx*dbasIx + dbasJy*dbasIy )
            p_DlocalMatrixA11(jdofe,idofe,iel) = p_DlocalMatrixA11(jdofe,idofe,iel) + dval
            p_DlocalMatrixA22(jdofe,idofe,iel) = p_DlocalMatrixA22(jdofe,idofe,iel) + dval
                                          
          end do ! idofe
          
        end do ! jdofe

      end do ! icubp
    
    end do ! iel
    
    ! Calculate the B1/B2/B1^T/B2^T blocks
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Outer loop over the DOFs i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rmatrixDataA13%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasIx = p_DbasTestA13(idofe,DER_DERIV2D_X,icubp,iel)
          dbasIy = p_DbasTestA13(idofe,DER_DERIV2D_Y,icubp,iel)
          
          ! Inner loop over the DOFs j=1..ndof, which corresponds to
          ! the basis function Phi_j:
          do jdofe=1,p_rmatrixDataA13%ndofTrial
            
            ! Fetch the contributions of the (trial) basis function Phi_j
            ! into dbasJ
            dbasJ = p_DbasTrialA13(jdofe,DER_FUNC,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local matrices.
            
            ! B-matrices
            p_DlocalMatrixA13(jdofe,idofe,iel) = p_DlocalMatrixA13(jdofe,idofe,iel) + &
                p_DcubWeight(icubp,iel) * ( -dbasJ*dbasIx )

            p_DlocalMatrixA23(jdofe,idofe,iel) = p_DlocalMatrixA23(jdofe,idofe,iel) + &
                p_DcubWeight(icubp,iel) * ( -dbasJ*dbasIy )

            ! B^T-matrices
            p_DlocalMatrixA31(idofe,jdofe,iel) = p_DlocalMatrixA31(idofe,jdofe,iel) + &
                p_DcubWeight(icubp,iel) * ( -dbasJ*dbasIx )

            p_DlocalMatrixA32(idofe,jdofe,iel) = p_DlocalMatrixA32(idofe,jdofe,iel) + &
                p_DcubWeight(icubp,iel) * ( -dbasJ*dbasIy )
                                          
          end do ! idofe
          
        end do ! jdofe

      end do ! icubp
    
    end do ! iel

  end subroutine

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
    real(DP) :: dbasI, dval1, dval2, dval3
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector1,p_DlocalVector2,p_DlocalVector3
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest1,p_DbasTest3
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData1,p_rvectorData2,p_rvectorData3
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    p_rvectorData1 => RvectorData(1)
    p_rvectorData2 => RvectorData(2)
    p_rvectorData3 => RvectorData(3)
    p_DlocalVector1 => RvectorData(1)%p_Dentry
    p_DlocalVector2 => RvectorData(2)%p_Dentry
    p_DlocalVector3 => RvectorData(3)%p_Dentry
    p_DbasTest1 => RvectorData(1)%p_DbasTest
    p_DbasTest3 => RvectorData(3)%p_DbasTest
    
    ! Calculate the RHS of the momentum equation
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Outer loop over the DOFs i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rvectorData1%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest1(idofe,DER_FUNC,icubp,iel)
          
          ! Values of the velocity RHS for the X1 and X2 component
          dval1 = 0.0_DP
          dval2 = 0.0_DP
          
          ! Multiply the values of the basis functions
          ! (1st derivatives) by the cubature weight and sum up
          ! into the local vectors.
          p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
              p_DcubWeight(icubp,iel) * dval1 * dbasI
          p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
              p_DcubWeight(icubp,iel) * dval2 * dbasI
          
        end do ! jdofe

      end do ! icubp
    
    end do ! iel
    
    ! Calculate the RHS of the continuity equation
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Outer loop over the DOFs i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rvectorData3%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest3(idofe,DER_FUNC,icubp,iel)
          
          ! Value of the pressure RHS
          dval3 = 0.0_DP 

          ! Multiply the values of the basis functions
          ! by the cubature weight and sum up into the local vectors.
          p_DlocalVector3(idofe,iel) = p_DlocalVector3(idofe,iel) + &
              p_DcubWeight(icubp,iel) * dval3 * dbasI
          
        end do ! jdofe

      end do ! icubp
    
    end do ! iel
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine st2d3_fcalc_L2error(Dintvalue,rassemblyData,rintegralAssembly,&
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
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dval1,dval2,dx,dy,dnu
    integer :: iel, icubp, icomp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
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

    ! Get the viscosity parameter
    dnu = rcollection%DquickAccess(1)

    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => revalVectors%p_RvectorData(icomp)%p_Ddata(:,:,DER_FUNC2D)

    Dintvalue = 0.0_DP
    
    ! Which component to calculate?
    select case (icomp)
    
    ! ----------
    ! X-velocity
    ! ----------
    case (1)
    
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the bubble function
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          
          dval1 = dy * (1.0_DP-dy)

          ! Get the error of the FEM function to the bubble function
          dval2 = p_Dfunc(icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the (squared) L2 error:
          Dintvalue = Dintvalue + &
              p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
            
        end do ! icubp
      
      end do ! iel
      
    ! ----------
    ! Y-velocity
    ! ----------
    case (2)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the bubble function
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          
          dval1 = 0.0_DP

          ! Get the error of the FEM function to the bubble function
          dval2 = p_Dfunc(icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the (squared) L2 error:
          Dintvalue = Dintvalue + &
              p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
            
        end do ! icubp
      
      end do ! iel

    ! ----------
    ! Pressure
    ! ----------
    case (3)
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the bubble function
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          
          dval1 = 2.0_DP * dnu * (1.0_DP-dx)

          ! Get the error of the FEM function to the bubble function
          dval2 = p_Dfunc(icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the (squared) L2 error:
          Dintvalue = Dintvalue + &
              p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
            
        end do ! icubp
      
      end do ! iel

    end select
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine st2d3_fcalc_H1error(Dintvalue,rassemblyData,rintegralAssembly,&
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
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>    

!<subroutine>

    ! Local variables
    real(DP) :: dderivX1,dderivY1,dderivX2,dderivY2,dx,dy
    integer :: iel, icubp, icomp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
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

    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_DderivX => revalVectors%p_RvectorData(icomp)%p_Ddata(:,:,DER_DERIV2D_X)
    p_DderivY => revalVectors%p_RvectorData(icomp)%p_Ddata(:,:,DER_DERIV2D_Y)

    Dintvalue = 0.0_DP

    ! Which component to calculate?
    select case (icomp)
    
    ! ----------
    ! X-velocity
    ! ----------
    case (1)
        
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the derivatives of the bubble function in the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          
          dderivX1 = 0.0_DP
          dderivY1 = 1.0_DP-2.0_DP*dy

          ! Get the error of the FEM function derivatives of the bubble function
          ! in the cubature point
          dderivX2 = p_DderivX(icubp,iel)
          dderivY2 = p_DderivY(icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the (squared) H1 error:
          Dintvalue = Dintvalue + p_DcubWeight(icubp,iel) * &
              ( (dderivX1 - dderivX2)**2 + (dderivY1 - dderivY2)**2 )
            
        end do ! icubp
      
      end do ! iel
    
    ! ----------
    ! Y-velocity
    ! ----------
    case (2)

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
          dderivX2 = p_DderivX(icubp,iel)
          dderivY2 = p_DderivY(icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the (squared) H1 error:
          Dintvalue = Dintvalue + p_DcubWeight(icubp,iel) * &
              ( (dderivX1 - dderivX2)**2 + (dderivY1 - dderivY2)**2 )
            
        end do ! icubp
      
      end do ! iel
    
    
    ! ----------
    ! Pressure
    ! ----------
    case (3)

      ! not used here.
      
    end select
            
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine stokes2d_3_block
  
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
    type(t_blockDiscretisation) :: rdiscretisation,rprjDiscretisation
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rvector,rrhs,rtempBlock,rprjVector
    
    ! Cubature information structure which defines the cubature formula.
    type(t_scalarCubatureInfo) :: rcubatureInfo
    
    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC, rprjDiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    
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

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    ! As we do not use a multigrid solver here, we will set the level
    ! to 5 instead of 7.
    NLMAX = 5
    
    ! Viscosity parameter:
    dnu = 1.0_DP

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
    call spdiscr_initBlockDiscr (rdiscretisation,3,&
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
                EL_EM30, rtriangulation, rboundary)
                
    ! ...and copy this structure also to the discretisation structure
    ! of the 2nd component (Y-velocity). This needs no additional memory,
    ! as both structures will share the same dynamic information afterwards.
    call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
        rdiscretisation%RspatialDiscr(2))

    ! For the pressure (3rd component), we set up a separate discretisation
    ! structure, as this uses different finite elements for trial and test
    ! functions.
    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
        EL_Q0, rdiscretisation%RspatialDiscr(3))

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                 
    ! Create an assembly information structure which tells the code
    ! the cubature formula to use. Standard: Gauss 3x3.
    call spdiscr_createDefCubStructure(&  
        rdiscretisation%RspatialDiscr(1),rcubatureInfo,CUB_GEN_AUTO_G3)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a block matrix structure for the operator
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Initialise the block matrix with default values based on
    ! the discretisation.
    call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrix)
    
    ! Inform the matrix that we build a saddle-point problem.
    ! Normally, imatrixSpec has the value LSYSBS_MSPEC_GENERAL,
    ! but probably some solvers can use the special structure later.
    rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
    
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
                                     
    ! Use it for the Y-velocity
    call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
        rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)

    ! Create the structure for the B-matrices using the 3rd
    ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(3),&
                                     LSYSSC_MATRIX9, rmatrix%RmatrixBlock(1,3),&
                                     rdiscretisation%RspatialDiscr(1))

    call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,3),&
        rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
        
    ! Create the structure for the B^T matrices by transposing the structre
    ! of the B-matrices. The structure of B1^T and B2^T is again the same.
    call lsyssc_transposeMatrix (&
        rmatrix%RmatrixBlock(1,3),rmatrix%RmatrixBlock(3,1),LSYSSC_TR_STRUCTURE)

    call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(3,1),&
        rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
        
    ! Now re-assign the block discretisation structure to all matrices
    call lsysbl_assignDiscrDirectMat (rmatrix,rdiscretisation)
    
    ! Allocate memory for the matrix
    call lsysbl_allocEmptyMatrix (rmatrix,LSYSSC_SETM_ZERO)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a block matrix entries
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Pass dnu via rcollection
    rcollection%DquickAccess(1) = dnu
    
    call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,&
          st2d3_fcalc_Stokes,rcubatureInfo=rcubatureInfo,rcollection=rcollection)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create RHS and solution vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Create a RHS and a solution vector based on the discretisation.
    ! Fill with zero.
    call lsysbl_createVectorBlock (rdiscretisation,rrhs,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rvector,.true.)

    ! Assemble the right hand side.
    call bma_buildVector (rrhs,BMA_CALC_STANDARD,&
          st2d3_fcalc_rhs,rcubatureInfo=rcubatureInfo,rcollection=rcollection)
                                
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Assembly of matrices/vectors finished
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
                                       getBoundaryValues_2D)
                             
    ! Edge 2 is Neumann boundary, so it is commented out.
    ! CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues_2D)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)

    ! Now continue with defining the boundary conditions of the Y-velocity:
    !
    ! Define edge 1.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    
    ! Edge with start- and endpoint.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Edge 2 is Neumann boundary, so it is commented out.
    ! CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,2,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues_2D)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)

    ! The pressure does not need boundary conditions.

    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    call vecfil_discreteBCrhs (rrhs,rdiscreteBC)
    call vecfil_discreteBCsol (rvector,rdiscreteBC)
    call matfil_discreteBC (rmatrix,rdiscreteBC)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a linear solver
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    call filter_initFilterChain (RfilterChain,nfilters)
    call filter_newFilterDiscBCDef (RfilterChain,nfilters,rdiscreteBC)

    ! Create a BiCGStab-solver with VANCA preconditioner.
    ! Attach the above filter chain to the solver, so that the solver
    ! automatically filters the vector during the solution process.
    nullify(p_rpreconditioner)
    call linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_GENERAL)
    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)

    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! We will allow the solver to perform 200 iterations
    p_rsolverNode%nmaxIterations = 200

    ! Attach the system matrix to the solver.
    call linsol_setMatrix(p_rsolverNode,rmatrix)
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode, ierror)

    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
      call sys_halt()
    end if
    
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
    call linsol_solveAdaptively (p_rsolverNode,rvector,rrhs,rtempBlock)

    ! -------------------------------------------------------------------------
    ! Solution vector projection for projecting a solution vector of
    ! an arbitrary element to Q1 such that it can be written into a VTK file
    ! -------------------------------------------------------------------------

    ! The solution vector is probably not in the way GMV likes it!
    ! GMV for example does not understand Q1~ vectors!
    ! Therefore, we first have to convert the vector to a form that
    ! GMV understands.
    ! GMV understands only Q1 solutions! So the task is now to create
    ! a Q1 solution from p_rvector and write that out.
    !
    ! For this purpose, first create a "derived" simple discretisation
    ! structure based on Q1 by copying the main guiding block discretisation
    ! structure and modifying the discretisation structures of the
    ! two velocity subvectors:
    
    call spdiscr_duplicateBlockDiscr (rdiscretisation,rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1), &
                 EL_Q1, rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
                 EL_Q1, rprjDiscretisation%RspatialDiscr(2))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVectorBlock (rprjDiscretisation,rprjVector,.true.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q0
    ! discretisation.
    !
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rprjDiscreteBC)
    !
    ! Edge 1 of boundary component 1, X-velocity.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)

    ! Edge with start- and endpoint.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,1,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Edge 2 is Neumann boundary, so it is commented out.
    ! CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rprjDiscretisation,1,&
    !                                    rboundaryRegion,rprjDiscreteBC,&
    !                                    getBoundaryValues_2D)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,1,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,1,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)

    ! Edge 1 of boundary component 1, Y-velocity.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
  
    ! Edge with start- and endpoint.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,2,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Edge 2 is Neumann boundary, so it is commented out.
    ! CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rprjDiscretisation,2,&
    !                                    rboundaryRegion,rprjDiscreteBC,&
    !                                    getBoundaryValues_2D)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,2,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rprjDiscretisation,2,&
                                       rboundaryRegion,rprjDiscreteBC,&
                                       getBoundaryValues_2D)

    ! Send the vector to the boundary-condition implementation filter.
    ! This modifies the vector according to the discrete boundary
    ! conditions.
    call vecfil_discreteBCsol (rprjVector,rprjDiscreteBC)
    
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    ! We can now start the postprocessing.
    ! Start UCD export to GMV file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
        trim(sucddir)//"/u2d_3_block.vtk")

    ! Write velocity field
    call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X- and Y-velocity at once to the ucd module.
    call ucd_addVarVertBasedVec(rexport,"velocity",p_Ddata,p_Ddata2)

    ! If we use the GMV exporter, we might replace the line above by the
    ! following two lines:
    !CALL ucd_addVariableVertexBased (rexport,"X-vel",UCD_VAR_XVELOCITY, p_Ddata)
    !CALL ucd_addVariableVertexBased (rexport,"Y-vel",UCD_VAR_YVELOCITY, p_Ddata2)
        
    ! Write pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    call ucd_addVariableElementBased (rexport,"pressure",UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Projection and VTK export finished.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Error calculation
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Store the viscosity parameter nu in the collection"s quick access array
    rcollection%DquickAccess(1) = dnu
    
    ! Set up revalVectors with the velocity/pressure vectors.
    call fev2_addVectorToEvalList(revalVectors,rvector%RvectorBlock(1),1)
    call fev2_addVectorToEvalList(revalVectors,rvector%RvectorBlock(2),1)
    call fev2_addVectorToEvalList(revalVectors,rvector%RvectorBlock(3),0)

    ! Calculate errors of velocity and pressure against analytic solutions.
    rcollection%IquickAccess(1) = 1  ! X-component, L2-error
    call bma_buildIntegral(DerrorUL2(1),BMA_CALC_STANDARD,st2d3_fcalc_L2error,&
        rcollection=rcollection,revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo)

    rcollection%IquickAccess(1) = 2  ! Y-component, L2-error
    call bma_buildIntegral(DerrorUL2(2),BMA_CALC_STANDARD,st2d3_fcalc_L2error,&
        rcollection=rcollection,revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo)

    rcollection%IquickAccess(1) = 3  ! Pressure, L2-error
    call bma_buildIntegral(derrorPL2(1),BMA_CALC_STANDARD,st2d3_fcalc_L2error,&
        rcollection=rcollection,revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo)

    rcollection%IquickAccess(1) = 1  ! X-component, H1-error
    call bma_buildIntegral(DerrorUH1(1),BMA_CALC_STANDARD,st2d3_fcalc_H1error,&
        rcollection=rcollection,revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo)
        
    rcollection%IquickAccess(1) = 2  ! Y-component, H1-error
    call bma_buildIntegral(DerrorUH1(2),BMA_CALC_STANDARD,st2d3_fcalc_H1error,&
        rcollection=rcollection,revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo)

    ! Print the errors.
    call output_lbrk()
    call output_line("|u - u_h|_L2 = " // trim(sys_sdEL(sqrt(DerrorUL2(1)), 10)) &
                                // " " // trim(sys_sdEL(sqrt(DerrorUL2(2)), 10)))
    call output_line("|u - u_h|_H1 = " // trim(sys_sdEL(sqrt(DerrorUH1(1)), 10)) &
                                // " " // trim(sys_sdEL(sqrt(DerrorUH1(2)), 10)))
    call output_line("|p - p_h|_L2 = " // trim(sys_sdEL(sqrt(derrorPL2(1)), 10)))
    
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
    
    ! Release the filter chain
    call filter_doneFilterChain (RfilterChain,nfilters)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rprjVector)
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rrhs)
    call lsysbl_releaseMatrix (rmatrix)
    
    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rprjDiscreteBC)
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rprjDiscretisation)
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation.
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)

  end subroutine

end module
