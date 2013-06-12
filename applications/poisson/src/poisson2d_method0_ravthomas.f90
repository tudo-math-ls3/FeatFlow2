!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method0_ravthomas </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!# For the assembly, the block assembly routines are used.
!#
!# The test applies a mixed formulation of the Poisson problem using the
!# lowest order Raviart-Thomas element.
!#
!# </purpose>
!##############################################################################

module poisson2d_method0_ravthomas

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use derivatives
  use element
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use filtersupport
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use scalarpde
  use element
  use ucd
  use pprocerror
  use matrixio
  
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop
  use feevaluation2
  use collection
    
  use poisson2d_callback
  
  implicit none

! The standard mixed formulation reads
!
!    (sigma,v)     + ( grad u, v ) = 0
!    (sigma,grad q)                = ( f, q )
!
! which is calculated using RT for sigma and P1 for u
! if the following define is active:

#define INSTABLEFORMULATION

! However, the formulation is instable. One actually needs
! one degree less in the polynomial space of u.
!
! If the above define is deactivated, the following stable alternative
! formulation is calculated using P0 for the velocity:
!
!    (sigma,v)     + ( u, div v ) = 0
!    (div sigma,q)                = ( f, q )
!
! However, this will give wrong boundary values for sigma
! (as the partial integration would normally imply some boundary terms)
! and a less cute GMV/VTK picture :-)

contains

  !****************************************************************************

!<subroutine>

  subroutine fcalc_systemmat(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Laplace operator at position (x,y) in a block matrix.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! at matrix position (1,1) with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%DquickAccess(1) = multiplier in front of the matrix.
    ! rcollection%IquickAccess(1) = x-coordinate in the block matrix
    ! rcollection%IquickAccess(2) = y-coordinate in the block matrix
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

    real(DP) :: dbasIx, dbasJx, dbasIy, dbasJy, dbasIz, dbasJz, dbasI1, dbasJ1, dbasI2, dbasJ2
    real(DP) :: dval,dbasI,dbasJ
    integer :: iel, icubp, idofe, jdofe, ivar, nvar
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix,p_DlocalMatrixT
    real(DP), dimension(:,:,:,:), pointer :: p_DlocalMatrixIntl
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    integer :: ix,iy,ndimfe,idimfe
    real(DP) :: dscale

    ! Get parameters
    dscale = 1.0_DP
    iy = 1
    ix = 1

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_rmatrixData => RmatrixData(iy,ix)
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest

    ndimfe = RmatrixData(iy,ix)%ndimfeTrial

    ! A mass matrix at position (1,1)

    ! Get the matrix data      
    p_DlocalMatrix => RmatrixData(1,1)%p_Dentry

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        do idimfe = 0,ndimfe-1
        
          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Phi_i:
          do idofe=1,p_rmatrixData%ndofTest

            ! Fetch the contributions of the (test) basis functions Phi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe+idimfe*p_rmatrixData%ndofTest,DER_FUNC,icubp,iel)

            ! Inner loop over the DOF's j=1..ndof, which corresponds to
            ! the basis function Phi_j:
            do jdofe=1,p_rmatrixData%ndofTrial

              ! Fetch the contributions of the (trial) basis function Phi_j
              ! into dbasJ
              dbasJ = p_DbasTrial(jdofe+idimfe*p_rmatrixData%ndofTrial,DER_FUNC,icubp,iel)

              ! Multiply the values of the basis functions
              ! (1st derivatives) by the cubature weight and sum up
              ! into the local matrices.
              p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                  dscale * p_DcubWeight(icubp,iel) * dbasJ*dbasI

            end do ! idofe

          end do ! jdofe

        end do

      end do ! icubp

    end do ! iel

    ! A matrix (t,grad(u)) at position (1,2), its transposed to (2,1)

    ! Get the matrix data      
    ix = 2
    iy = 1
    p_rmatrixData => RmatrixData(iy,ix)
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest

    p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry
    p_DlocalMatrixT => RmatrixData(ix,iy)%p_Dentry

#ifdef INSTABLEFORMULATION
    
    ! Instable formulation with P1.

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rmatrixData%ndofTest

          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI1 = p_DbasTest(idofe+0*p_rmatrixData%ndofTest,DER_FUNC,icubp,iel)
          dbasI2 = p_DbasTest(idofe+1*p_rmatrixData%ndofTest,DER_FUNC,icubp,iel)

          ! Inner loop over the DOF's j=1..ndof, which corresponds to
          ! the basis function Phi_j:
          do jdofe=1,p_rmatrixData%ndofTrial

            ! Fetch the contributions of the (trial) basis function Phi_j
            ! into dbasJ
            dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
            dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local matrices.
            dval = (-dscale) * p_DcubWeight(icubp,iel) * ( dbasJx*dbasI1 + dbasJy*dbasI2 )
            
            p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + dval
            p_DlocalMatrixT(idofe,jdofe,iel) = p_DlocalMatrixT(idofe,jdofe,iel) + dval
                

          end do ! idofe

        end do ! jdofe
          
      end do ! icubp

    end do ! iel
    
#else

    ! Stable formulation with P0

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rmatrixData%ndofTest

          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasIx = p_DbasTest(idofe+0*p_rmatrixData%ndofTest,DER_DERIV2D_X,icubp,iel)
          dbasIy = p_DbasTest(idofe+1*p_rmatrixData%ndofTest,DER_DERIV2D_y,icubp,iel)

          ! Inner loop over the DOF's j=1..ndof, which corresponds to
          ! the basis function Phi_j:
          do jdofe=1,p_rmatrixData%ndofTrial

            ! Fetch the contributions of the (trial) basis function Phi_j
            ! into dbasJ
            dbasJ1 = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local matrices.
            dval = (-dscale) * p_DcubWeight(icubp,iel) * dbasJ1*(dbasIx+dbasIy)
            
            p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + dval
            p_DlocalMatrixT(idofe,jdofe,iel) = p_DlocalMatrixT(idofe,jdofe,iel) + dval
                

          end do ! idofe

        end do ! jdofe
          
      end do ! icubp

    end do ! iel
    
#endif

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fcalc_rhs(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=32*y*(1-y)+32*x*(1-x).
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
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
    real(DP) :: dbasI, dval, dx, dy
    integer :: icomp
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! RHS in the 2nd component.
    icomp = 2

    ! Get the data arrays of the subvector
    p_rvectorData => RvectorData(icomp)
    p_DlocalVector => RvectorData(icomp)%p_Dentry
    p_DbasTest => RvectorData(icomp)%p_DbasTest
  
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Calculate the values of the RHS using the coordinates
        ! of the cubature points.
        dval = 32.0_DP*dy*(1.0_DP-dy) + 32_DP*dx*(1.0_DP-dx)
        
        ! Outer loop over the DOF's i=1..ndof on our current element,
        ! which corresponds to the (test) basis functions Phi_i:
        do idofe=1,p_rvectorData%ndofTest
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! into dbasI
          dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
          
          ! Multiply the values of the basis functions
          ! (1st derivatives) by the cubature weight and sum up
          ! into the local vectors.
          p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
              p_DcubWeight(icubp,iel) * (-1.0_DP) * dval * dbasI
          
        end do ! jdofe

      end do ! icubp
    
    end do ! iel
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine fcalc_midpValues(Dintvalue,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the squared L2-norm of an arbitrary finite element 
    ! function: <tex> ||v||^2 </tex>.
    ! If v has multiple components, the sum of the integrals of all
    ! components is returned.
    !
    ! The FEM function(s) must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rvectorAssembly

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
    integer :: iel
    real(DP), dimension(:,:,:), pointer :: p_Dfunc
    real(DP), dimension(:), pointer :: p_DdataX, p_DdataY
    integer, dimension(:), pointer :: p_Ielements
  
    Dintvalue = 0.0_DP
  
    call storage_getbase_double (rcollection%IquickAccess(1),p_DdataX)
    call storage_getbase_double (rcollection%IquickAccess(2),p_DdataY)

    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => revalVectors%p_RvectorData(1)%p_DdataVec(:,:,:,DER_FUNC)
    
    p_Ielements => rassemblyData%p_IelementList

    ! Loop over the elements in the current set.
    do iel = 1,nelements
    
      p_DdataX(p_Ielements(iel)) = p_Dfunc(1,1,iel)
      p_DdataY(p_Ielements(iel)) = p_Dfunc(2,1,iel)
    
    end do ! iel
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_0_ravthomas
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to VTK file
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

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo,rcubatureInfo2
    
    ! A matrix, a RHS vector, a solution vector and a temporary vector. 
    ! The RHS vector accepts the RHS of the problem, the solution vector
    ! accepts the solution. All are block vectors with only one block.
    type(t_matrixBlock) :: rmatSystem
    type(t_vectorBlock) :: rvecSol,rvecRhs,rvecTmp

    ! A set of variables describing the discrete boundary conditions.
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    ! A collection for passing values
    type(t_collection) :: rcollection

    ! Vector evaluation structure for the calculation of errors.
    type(t_fev2Vectors) :: revalVectors
    
    ! Output block for UCD output to VTK file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir

    real(DP), dimension(:), pointer :: p_DdataX, p_DdataY
  
    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 5
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Read the domain, read the mesh, refine
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./pre"

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, trim(spredir)//"/TRIA.prm")
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, trim(spredir)//"/TRIA.tri", rboundary)
     
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a discretisation structure which tells the code which
    ! finite element to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,2,&
                                 rtriangulation, rboundary)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_RT0_2D,rtriangulation, rboundary)
                                   
#ifdef INSTABLEFORMULATION                                   
    ! Instable formulation with P1
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(2), &
                                   EL_P1_2D,rtriangulation, rboundary)
#else
    ! Stable formulation with P0
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(2), &
                                   EL_P0_2D,rtriangulation, rboundary)
#endif

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up an cubature info structure to tell the code which cubature
    ! formula to use
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
                 
    ! Create an assembly information structure which tells the code
    ! the cubature formula to use. Standard: Gauss 3x3.
    call spdiscr_createDefCubStructure(&  
        rdiscretisation%RspatialDiscr(1),rcubatureInfo,CUB_G3MP_T)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create a 1x1 block matrix with the operator
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! At first, create a basic 1x1 block matrix based on the discretisation.
    call lsysbl_createMatrix (rdiscretisation,rmatSystem)
    
    ! We create a scalar matrix, based on the discretisation structure.
    call bilf_createMatrixStructure (rmatSystem, 1, 1, LSYSSC_MATRIX9)
    call bilf_createMatrixStructure (rmatSystem, 2, 1, LSYSSC_MATRIX9)
    call bilf_createMatrixStructure (rmatSystem, 1, 2, LSYSSC_MATRIX9)
    call bilf_createMatrixStructure (rmatSystem, 2, 2, LSYSSC_MATRIX9)
        
    ! And now to the entries of the matrix.
    !
    ! Allocate memory for the entries, fill with zero.
    call lsysbl_allocEmptyMatrix (rmatSystem,LSYSSC_SETM_ZERO)
    
    ! Pass the constant dnu=1.0 in front of the Laplace via rcollection.
    rcollection%DquickAccess(1) = 1.0_DP
    
    ! Laplace is to be assembled at position (x,y) = (1,1) in the block matrix
    rcollection%IquickAccess(1) = 1
    rcollection%IquickAccess(2) = 1

    ! The routine "bma_fcalc_laplace" can be found in "blockmatassembly.f90" as
    ! an example for a matrix calculation routine!
    call bma_buildMatrix (rmatSystem,BMA_CALC_STANDARD,&
          fcalc_systemmat,rcubatureInfo=rcubatureInfo,rcollection=rcollection)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Create RHS and solution vectors
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
    ! Next step: Create a RHS vector, a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (rdiscretisation,rvecRhs,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rvecSol,.true.)
    call lsysbl_createVectorBlock (rdiscretisation,rvecTmp,.true.)

    ! Discretise the RHS to get a discrete version of it.
    !
    ! The routine "bma_fcalc_rhsBubble" can be found in "blockmatassembly.f90" as
    ! an example for a RHS routine!
    call bma_buildVector (rvecRhs,BMA_CALC_STANDARD,&
        fcalc_rhs,rcubatureInfo=rcubatureInfo,rcollection=rcollection)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Assembly of matrices/vectors finished
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Discretise the boundary conditions and apply them to the matrix/RHS/sol.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! Now we have the raw problem. What is missing is the definition of the boundary
    ! conditions.
    ! For implementing boundary conditions, we use a `filter technique with
    ! discretised boundary conditions`. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rdiscreteBC)
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
    
    ! We use this boundary region and specify that we want to have Dirichlet
    ! boundary there. The following call does the following:
    ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
    !   We specify icomponent="1" to indicate that we set up the
    !   Dirichlet BC`s for the first (here: one and only) component in the
    !   solution vector.
    ! - Discretise the boundary condition so that the BC`s can be applied
    !   to matrices and vectors
    ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Now to the edge 2 of boundary component 1 the domain.
    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConrealbd (rdiscretisation,2,&
                                       rboundaryregion,rdiscretebc,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)


    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    call vecfil_discreteBCrhs (rvecRhs,rdiscreteBC)
    call vecfil_discreteBCsol (rvecSol,rdiscreteBC)
    call matfil_discreteBC (rmatSystem,rdiscreteBC)

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Set up a linear solver
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initUmfpack4 (p_rsolverNode)
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Attach the system matrix to the solver.
    call linsol_setMatrix(p_RsolverNode,rmatSystem)
    
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
    call linsol_solveAdaptively (p_rsolverNode,rvecSol,rvecRhs,rvecTmp)
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Postprocessing of the solution
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
    ! That is it, rvecSol now contains our solution. We can now
    ! start the postprocessing.

    ! Small trick. Use the integral calculatiuon routine to calculate the
    ! values of the first component in the midpoints of the elements.
    ! Create some destination arrays with NEL size and save their handles
    ! to the collection. IN the callback, the arrays are available then
    ! and we can save the values of the velocity vectors into them.
    call storage_new ("Poisson", "Xvel", rtriangulation%NEL, ST_DOUBLE, &
        rcollection%IquickAccess(1), ST_NEWBLOCK_ZERO)
    call storage_new ("Poisson", "Yvel", rtriangulation%NEL, ST_DOUBLE, &
        rcollection%IquickAccess(2), ST_NEWBLOCK_ZERO)

    call fev2_addVectorToEvalList(revalVectors,rvecSol%RvectorBlock(1),0)

    ! Gauss 1-pt rule = 1 Point per element in the center.
    call spdiscr_createDefCubStructure(&  
        rdiscretisation%RspatialDiscr(1),rcubatureInfo2,CUB_G1_T)
        
    ! Calculate the midpoint values.
    call bma_buildIntegral(derror,BMA_CALC_STANDARD,fcalc_midpValues,&
        rcollection=rcollection,revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo2)
        
    call spdiscr_releaseCubStructure(rcubatureInfo2)
    call fev2_releaseVectorList(revalVectors)
    
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./gmv"

    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       trim(sucddir)//"/u2d_0_ravthomas.vtk")
    
    ! Add the solution to the UCD exporter
    call storage_getbase_double (rcollection%IquickAccess(1),p_DdataX)
    call storage_getbase_double (rcollection%IquickAccess(2),p_DdataY)
    call ucd_addVariableElementBased (rexport, "sol1", UCD_VAR_STANDARD, p_DdataX)
    call ucd_addVariableElementBased (rexport, "sol2", UCD_VAR_STANDARD, p_DdataY)

#ifdef INSTABLEFORMULATION
    ! Instable formulation with P1
    call ucd_addVectorByVertex (rexport, "sol3", UCD_VAR_STANDARD, &
        rvecSol%RvectorBlock(2))
#else
    ! Stable formulation with P0
    call ucd_addVectorByElement (rexport, "sol3", UCD_VAR_STANDARD, &
        rvecSol%RvectorBlock(2))
#endif
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Release temporary storage.
    call storage_free(rcollection%IquickAccess(1))
    call storage_free(rcollection%IquickAccess(2))

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate the error to the reference function
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Set up revalVectors with the solution vectors.
    call fev2_addVectorToEvalList(revalVectors,rvecSol%RvectorBlock(2),0)

    ! L2-error, squared
    call bma_buildIntegral(derror,BMA_CALC_STANDARD,bma_fcalc_bubbleL2error,&
        rcollection=rcollection,revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo)

    derror = sqrt(derror)

    call output_line ("L2-error: " // sys_sdEL(derror,10) )

    ! H1-error, squared
!    call bma_buildIntegral(derror,BMA_CALC_STANDARD,bma_fcalc_bubbleH1error,&
!        rcollection=rcollection,revalVectors=revalVectors,&
!        rcubatureInfo=rcubatureInfo)
!    
!    derror = sqrt(derror)
!    
!    call output_line ("H1-error: " // sys_sdEL(derror,10) )
    
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
    call lsysbl_releaseVector (rvecTmp)
    call lsysbl_releaseVector (rvecSol)
    call lsysbl_releaseVector (rvecRhs)
    call lsysbl_releaseMatrix (rmatSystem)

    ! Release the cubature info structure.
    call spdiscr_releaseCubStructure(rcubatureInfo)

    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation.
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)
    
  end subroutine

end module
