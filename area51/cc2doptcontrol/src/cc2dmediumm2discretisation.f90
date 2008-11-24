!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2discretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic spatial discretisation related routines for 
!# CC2D. Here, matrix and RHS creation routines can be found as well as
!# routines to initialise/clean up discretisation structures.
!#
!# The following routines can be found here:
!#
!# 1.) cc_initDiscretisation
!#     -> Initialise the discretisation structure inside of the problem
!#        structure using the parameters from the INI/DAT files.
!#
!# 2.) cc_allocMatVec
!#     -> Allocates memory for vectors/matrices on all levels.
!#
!# 4.) cc_generateStaticMatrices
!#     -> Assembles matrix entries of static matrices (Stokes, B)
!#        on one level
!#
!# 5.) cc_generateBasicMatrices
!#     -> Assembles the matrix entries of all static matrices on all levels.
!#
!# 6.) cc_generateBasicRHS
!#     -> Generates a general RHS vector without any boundary conditions
!#        implemented
!#
!# 7.) cc_doneMatVec
!#     -> Cleanup of matrices/vectors, release all memory
!#
!# 8.) cc_doneDiscretisation
!#     -> Cleanup of the underlying discretisation structures
!#
!# 9.) cc_initInitialSolution
!#     -> Init solution vector according to parameters in the DAT file
!#
!# 10.) cc_writeSolution
!#      -> Write solution vector as configured in the DAT file.
!#
!# </purpose>
!##############################################################################

module cc2dmediumm2discretisation

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use stdoperators
  use vectorio
  
  use collection
  use convection
    
  use cc2dmediumm2basic
  use cc2dmedium_callback
  use cc2dmediumm2nonlinearcoreinit
  
  implicit none
  
contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_get1LevelDiscretisation (rparlist,rboundary,rtriangulation,&
      nequations,rdiscretisation)
  
!<description>
  ! Sets up the discretisation structures for the (primal) velocity and pressure
  ! as specified by the parameters in the parameter list rparList.
!</description>

!<input>
  ! The parameter list with the parameters in the DAT file.
  type(t_parlist), intent(IN) :: rparlist
  
  ! Number of equations.
  ! =3 for primal or dual discretisation structure.
  ! =6 for coupled primal/dual discretisation.
  integer, intent(IN) :: nequations
  
  ! Underlying boundary.
  type(t_boundary), intent(IN), target :: rboundary
  
  ! Underlying triangulation
  type(t_triangulation), intent(IN), target :: rtriangulation
!</input>

!<inputoutput>
  ! The discretisation structure to be set up.
  type(t_blockDiscretisation), intent(INOUT) :: rdiscretisation
!</inputoutput>

!</subroutine>

    integer :: ielementType,icubA,icubB,icubF
    character(LEN=SYS_NAMELEN) :: sstr

    ! Get parameters about the discretisation from the parameter list
    call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
                              'iElementType',ielementType,3)

    call parlst_getvalue_string (rparlist,'CC-DISCRETISATION',&
                                 'scubStokes',sstr,'')
    if (sstr .eq. '') then
      call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
                                'icubStokes',icubA,CUB_G2X2)
    else
      icubA = cub_igetID(sstr)
    end if

    call parlst_getvalue_string (rparlist,'CC-DISCRETISATION',&
                                'scubB',sstr,'')
    if (sstr .eq. '') then
      call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
                                'icubB',icubB,CUB_G2X2)
    else
      icubB = cub_igetID(sstr)
    end if

    call parlst_getvalue_string (rparlist,'CC-DISCRETISATION',&
                                 'scubF',sstr,'')
    if (sstr .eq. '') then
      call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
                                'icubF',icubF,CUB_G2X2)
    else
      icubF = cub_igetID(sstr)
    end if

    call spdiscr_initBlockDiscr2D (rdiscretisation,nequations,&
                                   rtriangulation, rboundary)

    select case (ielementType)
    case (0)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_E031,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case (1)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_E030,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case (2)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM31,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case (3)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM30,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case (4)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_Q2,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_QP1, icubB, rdiscretisation%RspatialDiscr(3))
                  
    case (5)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM30_UNPIVOTED,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case (6)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM30_UNSCALED,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case DEFAULT
      print *,'Unknown discretisation: iElementType = ',ielementType
      stop
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: I,k
  
  ! Number of equations in our problem. 
  ! velocity+velocity+pressure + dual velocity+dual velocity+dual pressure = 6
  integer, parameter :: nequations = 6
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisationMass
    
    integer :: icubM
    character(LEN=SYS_NAMELEN) :: sstr
    
    ! Set up discrezisation structures on all levels:

    do i=rproblem%NLMIN,rproblem%NLMAX
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies 3 blocks in the
      ! solution vector.
      call cc_get1LevelDiscretisation (rproblem%rparamList,&
          rproblem%rboundary,p_rtriangulation,nequations,&
          rproblem%RlevelInfo(i)%rdiscretisation)

      ! Save the discretisation structure to our local LevelInfo structure
      ! for later use.
      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisation

      ! -----------------------------------------------------------------------
      ! Optimal control extension
      ! -----------------------------------------------------------------------
      ! The variables 4,5,6 are discretised like the variables 1,2,3
      ! in the discretisation structure.
      ! Variable 1 = x-velocity  -->  Variable 4 = dual x-velocity
      ! Variable 2 = y-velocity  -->  Variable 5 = dual y-velocity
      ! Variable 3 = pressure    -->  Variable 6 = dual pressure
      call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(1),&
                                     p_rdiscretisation%RspatialDiscr(4))
      call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(2),&
                                     p_rdiscretisation%RspatialDiscr(5))
      call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(3),&
                                     p_rdiscretisation%RspatialDiscr(6))
          
      ! -----------------------------------------------------------------------
      ! Separated discretisation structures for primal and dual problem
      ! -----------------------------------------------------------------------
      ! Create a separate block discretisation structure, only for the primal
      ! space. Make a copy of the one we have, i.e. reuse it.
      ! Mark the substructures as being a copy from another discretisation structure,
      ! sharing all information with part 1..3 of the global problem.
      call spdiscr_deriveBlockDiscr (rproblem%RlevelInfo(i)%rdiscretisation, &
          rproblem%RlevelInfo(i)%rdiscretisationPrimal, 1,3)
      
      ! -----------------------------------------------------------------------
      ! Mass matrices
      ! -----------------------------------------------------------------------

      ! Initialise a discretisation structure for the mass matrix.
      ! Copy the discretisation structure of the first (Stokes) block
      ! and replace the cubature-formula identifier by that which is to be
      ! used for the mass matrix.
      allocate(p_rdiscretisationMass)
      call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(1),&
          rproblem%RlevelInfo(i)%rdiscretisationMass,.true.)
      p_rdiscretisationMass => rproblem%RlevelInfo(i)%rdiscretisationMass
      
      call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                  'scubStokes',sstr,'')
      if (sstr .eq. '') then
        call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                  'icubStokes',icubM,CUB_G2X2)
      else
        icubM = cub_igetID(sstr)
      end if

      ! Initialise the cubature formula appropriately.
      do k = 1,p_rdiscretisationMass%inumFESpaces
        p_rdiscretisationMass%RelementDistr(k)%ccubTypeBilForm = icubM
      end do

    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Remove the main block discretisation structure. 
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation,.true.)
      
      ! Release the block discretisation structures of the primal and dual
      ! space.
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisationPrimal,.true.)
      
      ! -----------------------------------------------------------------------
      ! Mass matrix problem
      ! -----------------------------------------------------------------------

      ! Release the mass matrix discretisation.
      call spdiscr_releaseDiscr(rproblem%RlevelInfo(i)%rdiscretisationMass)

    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_allocMatVec (rproblem,rvector,rrhs)
  
!<description>
  ! Allocates memory for all matrices and vectors of the problem on the heap
  ! by evaluating the parameters in the problem structure.
  ! Matrices/vectors of global importance are added to the collection
  ! structure of the problem, given in rproblem. Matrix/vector entries
  ! are not calculated, they are left uninitialised.\\\\
  !
  ! The following matrices/vectors are added to the collection:
  ! Level independent data:\\
  !  'RHS'         - Right hand side vector\\
  !  'SOLUTION'    - Solution vector\\
  !
  ! On every level:\\
  !  'STOKES'    - Stokes matrix\\
  !  'SYSTEMMAT' - Global system (block) matrix\\
  !  'RTEMPVEC'  - Temporary vector
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! A vector structure for the solution vector. The structure is initialised,
  ! memory is allocated for the data entries.
  type(t_vectorBlock), intent(INOUT) :: rvector

  ! A vector structure for the RHS vector. The structure is initialised,
  ! memory is allocated for the data entries.
  type(t_vectorBlock), intent(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,cmatBuildType
  
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixScalar), pointer :: p_rmatrixStokes
    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient
    type(t_vectorBlock), pointer :: p_rtempVector

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    ! When the jump stabilisation is used, we have to create an extended
    ! matrix stencil!
    cmatBuildType = BILF_MATC_ELEMENTBASED
    
    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
                              'IUPWIND1', i)
    if (i .eq. 2) cmatBuildType = BILF_MATC_EDGEBASED

    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
                              'IUPWIND2', i)
    if (i .eq. 2) cmatBuildType = BILF_MATC_EDGEBASED
  
    ! Initialise all levels...
    do i=rproblem%NLMIN,rproblem%NLMAX

      ! -----------------------------------------------------------------------
      ! Basic (Navier-) Stokes problem
      ! -----------------------------------------------------------------------

      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisation
      
      ! The global system looks as follows:
      !
      !    ( A11       B1  M/a          )  ( u )  =  ( f1)
      !    (      A22  B2       M/a     )  ( v )     ( f2)
      !    ( B1^T B2^T                  )  ( p )     ( fp)
      !    ( -M            A44  A45  B1 )  ( l1)     ( z1)
      !    (      -M       A54  A55  B2 )  ( l2)     ( z2)
      !    (               B1^T B2^T    )  ( xi)     ( 0 )
      !
      ! with Aii = L + nonlinear Convection, a=alphaC from the given equation.
      ! We compute in advance a standard Stokes matrix L which can be added 
      ! later to the convection matrix, resulting in the nonlinear system matrix.
      !
      ! At first, we create 'template' matrices that define the structure
      ! of each of the submatrices in that global system. These matrices
      ! are later used to 'derive' the actual Laplace/Stokes matrices
      ! by sharing the structure (KCOL/KLD).
      !
      ! Get a pointer to the template FEM matrix. This is used for the
      ! Laplace/Stokes matrix and probably for the mass matrix.
      
      p_rmatrixTemplateFEM => rproblem%RlevelInfo(i)%rmatrixTemplateFEM

      ! Create the matrix structure
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixTemplateFEM,cconstrType=cmatBuildType)

      ! In the global system, there are two gradient matrices B1 and B2.
      ! Create a template matrix that defines their structure.
      p_rmatrixTemplateGradient => rproblem%RlevelInfo(i)%rmatrixTemplateGradient
      
      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(3),LSYSSC_MATRIX9,&
                p_rmatrixTemplateGradient,p_rdiscretisation%RspatialDiscr(1))
      
      ! Ok, now we use the matrices from above to create the actual submatrices
      ! that are used in the global system.
      !
      ! Get a pointer to the (scalar) Stokes matrix:
      p_rmatrixStokes => rproblem%RlevelInfo(i)%rmatrixStokes
      
      ! Connect the Stokes matrix to the template FEM matrix such that they
      ! use the same structure.
      !
      ! Don't create a content array yet, it will be created by 
      ! the assembly routines later.
      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  p_rmatrixStokes,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      
      ! Allocate memory for the entries; don't initialise the memory.
      call lsyssc_allocEmptyMatrix (p_rmatrixStokes,LSYSSC_SETM_UNDEFINED)
      
      ! In the global system, there are two coupling matrices B1 and B2.
      ! Both have the structure of the template gradient matrix.
      ! So connect the two B-matrices to the template gradient matrix
      ! such that they share the same structure.
      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
      !
      ! Don't create a content array yet, it will be created by 
      ! the assembly routines later.
      
      call lsyssc_duplicateMatrix (p_rmatrixTemplateGradient,&
                  rproblem%RlevelInfo(i)%rmatrixB1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
      call lsyssc_duplicateMatrix (p_rmatrixTemplateGradient,&
                  rproblem%RlevelInfo(i)%rmatrixB2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
      ! Allocate memory for the entries; don't initialise the memory.
      call lsyssc_allocEmptyMatrix (rproblem%RlevelInfo(i)%rmatrixB1,&
                                           LSYSSC_SETM_UNDEFINED)
      call lsyssc_allocEmptyMatrix (rproblem%RlevelInfo(i)%rmatrixB2,&
                                           LSYSSC_SETM_UNDEFINED)
                                           
      ! -----------------------------------------------------------------------
      ! Allocate memory for an identity matrix in the size of the pressure.
      ! Attach a discretisation structure that describes the pressure element
      ! as trial space.
      call lsyssc_createDiagMatrixStruc (&
          rproblem%RlevelInfo(i)%rmatrixIdentityPressure,&
          rproblem%RlevelInfo(i)%rmatrixB1%NCOLS,LSYSSC_MATRIX9)
      call lsyssc_allocEmptyMatrix (rproblem%RlevelInfo(i)%rmatrixIdentityPressure,&
          LSYSSC_SETM_UNDEFINED)
      call lsyssc_assignDiscretDirectMat (rproblem%RlevelInfo(i)%rmatrixIdentityPressure,&
          p_rdiscretisation%RspatialDiscr(3))
          
      ! -----------------------------------------------------------------------
      ! Temporary vectors
      !
      ! Now on all levels except for the maximum one, create a temporary 
      ! vector on that level, based on the matrix template.
      ! It's used for building the matrices on lower levels.
      if (i .lt. rproblem%NLMAX) then
        call lsysbl_createVecBlockByDiscr (&
            rproblem%RlevelInfo(i)%rdiscretisation,&
            rproblem%RlevelInfo(i)%rtempVector1,.true.)
        call lsysbl_createVecBlockByDiscr (&
            rproblem%RlevelInfo(i)%rdiscretisation,&
            rproblem%RlevelInfo(i)%rtempVector2,.true.)
        call lsysbl_createVecBlockByDiscr (&
            rproblem%RlevelInfo(i)%rdiscretisation,&
            rproblem%RlevelInfo(i)%rtempVector3,.true.)

        p_rtempVector => rproblem%RlevelInfo(i)%rtempVector2
            
        ! The temp vectors for the primal and dual system share their memory
        ! with that temp vector.
        call lsysbl_deriveSubvector(p_rtempVector,&
            rproblem%RlevelInfo(i)%rtempVectorPrimal,1,3,.true.)
        call lsysbl_deriveSubvector(p_rtempVector,&
            rproblem%RlevelInfo(i)%rtempVectorDual,4,6,.true.)
      end if
      
    end do
    
    ! (Only) on the finest level, we need to allocate a RHS vector
    ! and a solution vector.
    !
    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template.
    ! Initialise the vectors with 0.
    call lsysbl_createVecBlockByDiscr (&
        rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation,rrhs,.true.)
    call lsysbl_createVecBlockByDiscr (&
        rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation,rvector,.true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateStaticMatrices (rproblem,rlevelInfo)
  
!<description>
  ! Calculates entries of all static matrices (Stokes, B-matrices,...)
  ! in the specified problem structure, i.e. the entries of all matrices 
  ! that do not change during the computation or which serve as a template for
  ! generating other matrices.
  !
  ! Memory for those matrices must have been allocated before with 
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem

  ! A level-info structure. The static matrices in this structure are generated.
  type(t_problem_lvl), intent(INOUT),target :: rlevelInfo
!</inputoutput>

!</subroutine>

    ! A pointer to the Stokes and mass matrix
    type(t_matrixScalar), pointer :: p_rmatrixStokes,p_rmatrixMass

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! -----------------------------------------------------------------------
    ! Basic (Navier-) Stokes problem
    ! -----------------------------------------------------------------------
  
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%rdiscretisation
    
    ! Get a pointer to the (scalar) Stokes matrix:
    p_rmatrixStokes => rlevelInfo%rmatrixStokes
    
    call stdop_assembleLaplaceMatrix (p_rmatrixStokes,.true.,rproblem%dnu)
    
    ! In the global system, there are two coupling matrices B1 and B2.
    !
    ! Build the first pressure matrix B1.
    call stdop_assembleSimpleMatrix (rlevelInfo%rmatrixB1,&
        DER_FUNC,DER_DERIV_X,-1.0_DP)

    ! Build the second pressure matrix B2.
    call stdop_assembleSimpleMatrix (rlevelInfo%rmatrixB2,&
        DER_FUNC,DER_DERIV_Y,-1.0_DP)
                                
    ! -----------------------------------------------------------------------
    ! Mass matrices
    ! -----------------------------------------------------------------------

    p_rmatrixMass => rlevelInfo%rmatrixMass

    ! If there is an existing mass matrix, release it.
    call lsyssc_releaseMatrix (rlevelInfo%rmatrixMass)

    ! Generate mass matrix. The matrix has basically the same structure as
    ! our template FEM matrix, so we can take that.
    call lsyssc_duplicateMatrix (rlevelInfo%rmatrixTemplateFEM,&
                p_rmatrixMass,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
    ! Change the discretisation structure of the mass matrix to the
    ! correct one; at the moment it points to the discretisation structure
    ! of the Stokes matrix...
    call lsyssc_assignDiscretDirectMat (p_rmatrixMass,&
        rlevelInfo%rdiscretisationMass)

    ! Call the standard matrix setup routine to build the matrix.                    
    call stdop_assembleSimpleMatrix (p_rmatrixMass,DER_FUNC,DER_FUNC)

    ! Clean up the collection (as we are done with the assembly, that's it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

    ! -----------------------------------------------------------------------
    ! Initialise the identity matrix
    call lsyssc_initialiseIdentityMatrix (rlevelInfo%rmatrixIdentityPressure)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateBasicMatrices (rproblem)
  
!<description>
  ! Calculates the entries of all static matrices (Mass, B,...) on all levels.
  !
  ! Memory for those matrices must have been allocated before with 
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    do i=rproblem%NLMIN,rproblem%NLMAX
      call cc_generateStaticMatrices (&
          rproblem,rproblem%RlevelInfo(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateBasicRHS (rproblem,rrhs)
  
!<description>
  ! Calculates the entries of the basic right-hand-side vector on the finest
  ! level. Boundary conditions or similar things are not implemented into 
  ! the vector.
  ! Memory for the RHS vector must have been allocated in advance.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! The RHS vector which is to be filled with data.
  type(t_vectorBlock), intent(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Drhs
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get a pointer to the RHS on the finest level as well as to the
    ! block discretisation structure:
    p_rdiscretisation => rrhs%p_rblockDiscr
    
    ! The vector structure is already prepared, but the entries are missing.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the 
    ! first block.
    !
    ! We pass our collection structure as well to this routine, 
    ! so the callback routine has access to everything what is
    ! in the collection.
    !
    ! Note that the vector is unsorted after calling this routine!
    !
    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! Discretise the X-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
              rrhs%RvectorBlock(1),coeff_RHS_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
              rrhs%RvectorBlock(2),coeff_RHS_y,&
              rproblem%rcollection)
                                
    ! The third subvector must be zero initially - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(rrhs%RvectorBlock(3))
    
    ! The RHS terms for the dual equation are calculated similarly using
    ! the desired 'target' flow field.
    !
    ! Discretise the X-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
              rrhs%RvectorBlock(4),coeff_TARGET_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
              rrhs%RvectorBlock(5),coeff_TARGET_y,&
              rproblem%rcollection)
      
    ! Depending on the formulation, to get a reference dual velocity,
    ! it might be necessary to switch the sign of the target velocity field 
    ! because the RHS of the dual equation is '-z'!
    ! Remember that it this case the signs of the mass matrices that couple
    ! primal and dual velocity must be changed, too!
    
    if (rproblem%roptcontrol%ispaceTimeFormulation .eq. 0) then
      call lsyssc_scaleVector (rrhs%RvectorBlock(4),-1.0_DP)
      call lsyssc_scaleVector (rrhs%RvectorBlock(5),-1.0_DP)
    end if

    ! Dual pressure RHS is =0.
    call lsyssc_clearVector(rrhs%RvectorBlock(6))
                                
    ! Clean up the collection (as we are done with the assembly, that's it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initInitialSolution (rproblem,rvector)
  
!<description>
  ! Initialises the initial solution vector into rvector. Depending on the settings
  ! in the DAT file this is either zero or read from a file.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(IN) :: rproblem
!</input>

!<inputoutput>
  ! The solution vector to be initialised. Must be set up according to the
  ! maximum level NLMAX in rproblem!
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: istart,nblock,iblock
    type(t_vectorBlock) :: rvector1,rvector2,rtempVector
    type(t_vectorScalar) :: rvectorTemp
    character(LEN=SYS_STRLEN) :: sarray,sfile,sfileString
    integer :: ilev
    integer(PREC_VECIDX) :: NEQ
    type(t_interlevelProjectionBlock) :: rprojection 

    ! Get the parameter what to do with rvector
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'isolutionStart',istart,0)
    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'ssolutionStart',sfileString,'')

    ! Create a temp vector at level NLMAX-istart+1.
    ilev = rproblem%NLMAX-abs(istart)+1
    
    if (ilev .lt. rproblem%NLMIN) then
      print *,'Warning: Level of start vector is < NLMIN! Initialising with zero!'
      istart = 0
    end if
    
    ! Inít with zero? Read from file?
    if (istart .eq. 0) then
      ! Init with zero
      call lsysbl_clearVector (rvector)
    else
      ! Remove possible ''-characters
      read(sfileString,*) sfile

      call lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev)%rdiscretisation,rvector1,.true.)
      
      ! Read in the vector, create a new block vector.
      call vecio_readBlockVectorHR (&
          rtempVector, sarray, .true., 0, sfile, istart .gt. 0)
      
      ! Copy velocity/pressure vectors from the temp to the actual
      ! vector. If the temp vector provides dual solutions, copy
      ! them as well.
      nblock = 3
      if (rtempVector%nblocks .ge. 2*(NDIM2D+1)) then
        nblock = 6
      end if
      
      do iblock = 1,nblock
        if (rtempVector%RvectorBlock(iblock)%NEQ .ne. rvector1%RvectorBlock(iblock)%NEQ) then
          call output_line ('Start vector incompatible!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'mysubroutine')
          call sys_halt()
        end if
        
        call lsyssc_duplicateVector (rtempVector%RvectorBlock(iblock),&
            rvector1%RvectorBlock(iblock),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
      end do
      
      ! If the vector is on level < NLMAX, we have to bring it to level NLMAX
      do while (ilev .lt. rproblem%NLMAX)
        
        ! Initialise a vector for the higher level and a prolongation structure.
        call lsysbl_createVectorBlock (&
            rproblem%RlevelInfo(ilev+1)%rdiscretisation,rvector2,.false.)
        
        call mlprj_initProjectionVec (rprojection,rvector2)
        
        ! Prolongate to the next higher level.

        NEQ = mlprj_getTempMemoryVec (rprojection,rvector1,rvector2)
        if (NEQ .ne. 0) call lsyssc_createVector (rvectorTemp,NEQ,.false.)
        call mlprj_performProlongation (rprojection,rvector1, &
                                        rvector2,rvectorTemp)
        if (NEQ .ne. 0) call lsyssc_releaseVector (rvectorTemp)
        
        ! Swap rvector1 and rvector2. Release the coarse grid vector.
        call lsysbl_swapVectors (rvector1,rvector2)
        call lsysbl_releaseVector (rvector2)
        
        call mlprj_doneProjection (rprojection)
        
        ! rvector1 is now on level ilev+1
        ilev = ilev+1
        
      end do
      
      ! Copy the resulting vector rvector1 to the output.
      call lsysbl_copyVector (rvector1,rvector)
      
      ! Release the temp vector
      call lsysbl_releaseVector (rvector1)
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_writeSolution (rproblem,rvector)
  
!<description>
  ! Writes a solution vector rvector to a file as configured in the parameters
  ! in the DAT file.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(IN) :: rproblem

  ! The solution vector to be written out. Must be set up according to the
  ! maximum level NLMAX in rproblem!
  type(t_vectorBlock), intent(IN) :: rvector
!</input>

!</subroutine>

    ! local variables
    integer(I32) :: idestLevel
    type(t_vectorBlock) :: rvector1,rvector2
    type(t_vectorScalar) :: rvectorTemp
    character(LEN=SYS_STRLEN) :: sfile,sfileString
    integer :: ilev
    integer(PREC_VECIDX) :: NEQ
    type(t_interlevelProjectionBlock) :: rprojection 
    logical :: bformatted

    ! Get the parameter what to do with rvector
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iSolutionWrite',idestLevel,0)
    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'sSolutionWrite',sfileString,'')

    if (idestLevel .eq. 0) return ! nothing to do.
    
    ! Remove possible ''-characters
    read(sfileString,*) sfile
    
    bformatted = idestLevel .gt. 0
    idestLevel = rproblem%NLMAX-abs(idestLevel)+1 ! level where to write out

    if (idestLevel .lt. rproblem%NLMIN) then
      print *,'Warning: Level for solution vector is < NLMIN! &
              &Writing out at level NLMIN!'
      idestLevel = rproblem%NLMIN
    end if
    
    ! Interpolate the solution down to level istart.
    call lsysbl_copyVector (rvector,rvector1)   ! creates new rvector1!

    do ilev = rproblem%NLMAX,idestLevel+1,-1
      
      ! Initialise a vector for the lower level and a prolongation structure.
      call lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev-1)%rdiscretisation,rvector2,.false.)
      
      call mlprj_initProjectionVec (rprojection,rvector2)
      
      ! Interpolate to the next higher level.
      ! (Don't 'restrict'! Restriction would be for the dual space = RHS vectors!)

      NEQ = mlprj_getTempMemoryVec (rprojection,rvector2,rvector1)
      if (NEQ .ne. 0) call lsyssc_createVector (rvectorTemp,NEQ,.false.)
      call mlprj_performInterpolation (rprojection,rvector2,rvector1, &
                                       rvectorTemp)
      if (NEQ .ne. 0) call lsyssc_releaseVector (rvectorTemp)
      
      ! Swap rvector1 and rvector2. Release the fine grid vector.
      call lsysbl_swapVectors (rvector1,rvector2)
      call lsysbl_releaseVector (rvector2)
      
      call mlprj_doneProjection (rprojection)
      
    end do

    ! Write out the solution.
    if (bformatted) then
      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,&
                                    0, sfile, '(E22.15)')
    else
      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,0, sfile)
    end if

    ! Release temp memory.
    call lsysbl_releaseVector (rvector1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneMatVec (rproblem,rvector,rrhs)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! A vector structure for the solution vector. The structure is cleaned up,
  ! memory is released.
  type(t_vectorBlock), intent(INOUT) :: rvector

  ! A vector structure for the RHS vector. The structure is cleaned up,
  ! memory is released.
  type(t_vectorBlock), intent(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

    integer :: i

    ! Release matrices and vectors on all levels
    do i=rproblem%NLMAX,rproblem%NLMIN,-1

      ! If there is an existing mass matrix, release it.
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixMass)

      ! Release Stokes, B1 and B2 matrix
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB2)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB1)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixStokes)
      
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixIdentityPressure)
      
      ! Release the template matrices. This is the point, where the
      ! memory of the matrix structure is released.
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixTemplateGradient)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixTemplateFEM)
      
      ! Remove the temp vector that was used for interpolating the solution
      ! from higher to lower levels in the nonlinear iteration.
      if (i .lt. rproblem%NLMAX) then
        call lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVectorPrimal)
        call lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVectorDual)
        call lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector3)
        call lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector2)
        call lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector1)
      end if
      
    end do

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rrhs)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initOptControl (rproblem)
  
!<description>
  ! Initialises the optimal control problem. Reads parameters from the DAT
  ! file and evaluates them.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The roptcontrol substructure is initialised.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

    character(SYS_STRLEN) :: spar

    ! Read in the parameters
    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dalphaC',rproblem%roptcontrol%dalphaC,1.0_DP)
        
    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dgammaC',rproblem%roptcontrol%dgammaC,0.0_DP)
        
    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itypeTargetFlow',rproblem%roptcontrol%itypeTargetFlow,0)

    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ilevelTargetFlow',rproblem%roptcontrol%ilevelTargetFlow,0)

    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'smeshTargetFlow',spar,'''''')
    read(spar,*) rproblem%roptcontrol%smeshTargetFlow

    call parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'stargetFlow',spar,'')
    read(spar,*) rproblem%roptcontrol%stargetFlow

    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itargetFlowDelta',rproblem%roptcontrol%itargetFlowDelta,1)

    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itargetFlowTimesteps',rproblem%roptcontrol%itargetFlowTimesteps,-1)
    
    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ispaceTimeFormulation',rproblem%roptcontrol%ispaceTimeFormulation,0)
    
    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'iconvectionExplicit',rproblem%roptcontrol%iconvectionExplicit,0)
    
    call parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ccontrolConstraints',rproblem%roptcontrol%ccontrolConstraints,0)

    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dumin1',rproblem%roptcontrol%dumin1,-1.0E10_DP)

    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dumax1',rproblem%roptcontrol%dumax1,1.0E10_DP)

    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dumin2',rproblem%roptcontrol%dumin2,-1.0E10_DP)

    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dumax2',rproblem%roptcontrol%dumax2,1.0E10_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneOptControl (rproblem)
  
!<description>
  ! Cleans up the structure for the optimal control problem. 
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The roptcontrol substructure is initialised.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

    rproblem%roptcontrol%dalphaC = 1.0_DP
    rproblem%roptcontrol%dgammaC = 0.0_DP
    rproblem%roptcontrol%itypeTargetFlow = 0
    rproblem%roptcontrol%stargetFlow = ''

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initTargetFlow (rproblem,dstartTime,dendTime,ntimesteps)
  
!<description>
  ! Reads in the target flow of the optimal control problem.
!</description>

!<input>
  ! OPTIONAL: Start time in a nonstationary simulation.
  ! This parameter must be present in a nonstatinoary simulation. In a stationary
  ! simulation, the parameter can be skipped.
  real(DP), intent(IN), optional :: dstartTime

  ! OPTIONAL: Start time in a nonstationary simulation.
  ! This parameter must be present in a nonstatinoary simulation. In a stationary
  ! simulation, the parameter can be skipped.
  real(DP), intent(IN), optional :: dendTime

  ! OPTIONAL: Total number of timesteps in a nonstationary simulation.
  ! This parameter must be present in a nonstatinoary simulation. In a stationary
  ! simulation, the parameter can be skipped.
  integer, intent(IN), optional :: ntimesteps
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The target flow is saved to the roptcontrol substructure.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>
    
    ! local variables
    character(SYS_STRLEN) :: sarray
    integer :: iref
    integer :: nactTimesteps
    
    ! At first, what is with the triangulation and discretisation?
    ! Are we using the same mesh / discretisation as we use for the
    ! solution?
    if (rproblem%roptcontrol%smeshTargetFlow .eq. '') then
      ! Is the triangulation inside of the bounds specified by NLMIN/NLMAX
      ! or not? If yes, we can use a predefined triangulation.
      ! If not, we have to refine the mesh until we reach the destination
      ! level.
      if (rproblem%roptcontrol%ilevelTargetFlow .le. rproblem%NLMAX) then
        ! Create references to the existing triangulation
        rproblem%roptcontrol%p_rtriangulationTargetFlow => &
          rproblem%RlevelInfo(rproblem%roptcontrol%ilevelTargetFlow)%rtriangulation
          
        rproblem%roptcontrol%p_rdiscrTargetFlow => &
          rproblem%RlevelInfo(rproblem%roptcontrol%ilevelTargetFlow)%rdiscretisationPrimal
      else
        iref = rproblem%roptcontrol%ilevelTargetFlow - rproblem%NLMAX
        allocate(rproblem%roptcontrol%p_rtriangulationTargetFlow)
        call tria_duplicate (&
            rproblem%RlevelInfo(rproblem%NLMAX)%rtriangulation,&
            rproblem%roptcontrol%p_rtriangulationTargetFlow,TR_SHARE_ALL)
            
        ! Refine
        call tria_quickRefine2LevelOrdering(iref,&
            rproblem%roptcontrol%p_rtriangulationTargetFlow,rproblem%rboundary)
        call tria_initStandardMeshFromRaw (&
            rproblem%roptcontrol%p_rtriangulationTargetFlow,rproblem%rboundary)
            
        ! Create a discretisation structure corresponding to that mesh.
        allocate(rproblem%roptcontrol%p_rdiscrTargetFlow)
        call cc_get1LevelDiscretisation (rproblem%rparamList,rproblem%rboundary,&
            rproblem%roptcontrol%p_rtriangulationTargetFlow,NDIM2D+1,&
            rproblem%roptcontrol%p_rdiscrTargetFlow)
        
      end if
    else
      ! No, we use a different mesh. We must read it from discr and refine it!
      allocate(rproblem%roptcontrol%p_rtriangulationTargetFlow)
      call tria_readTriFile2D (rproblem%roptcontrol%p_rtriangulationTargetFlow, &
          rproblem%roptcontrol%smeshTargetFlow, rproblem%rboundary)

      ! Refine
      call tria_quickRefine2LevelOrdering(rproblem%roptcontrol%ilevelTargetFlow-1,&
          rproblem%roptcontrol%p_rtriangulationTargetFlow,rproblem%rboundary)
      call tria_initStandardMeshFromRaw (&
          rproblem%roptcontrol%p_rtriangulationTargetFlow,rproblem%rboundary)
      
      ! Create a discretisation structure corresponding to that mesh.
      allocate(rproblem%roptcontrol%p_rdiscrTargetFlow)
      call cc_get1LevelDiscretisation (rproblem%rparamList,rproblem%rboundary,&
          rproblem%roptcontrol%p_rtriangulationTargetFlow,NDIM2D+1,&
          rproblem%roptcontrol%p_rdiscrTargetFlow)
    end if

    ! Probably create a vector containing the target flow
    select case (rproblem%roptcontrol%itypeTargetFlow)
    case (1,3)
    
      ! Stationary target flow. Read in the vector
      call vecio_readBlockVectorHR (&
          rproblem%roptcontrol%rtargetFlow, sarray, .true., &
          0, rproblem%roptcontrol%stargetFlow, .true.)
          
      ! Connect the vector to the discretisation structure.
      call lsysbl_enforceStructureDiscr (rproblem%roptcontrol%p_rdiscrTargetFlow,&
          rproblem%roptcontrol%rtargetFlow)
          
    case (2,4)
    
      ! Nonstationary target flow given by a sequence of files.

      if (.not. present(ntimesteps)) then
        call output_line ('ntimesteps not available!')
        call sys_halt()
      end if

      ! Read the first vector. It defines the shape of all vectors in rx!
      call vecio_readBlockVectorHR (&
          rproblem%roptcontrol%rtargetFlow, sarray, .true., &
          0, trim(rproblem%roptcontrol%stargetFlow)//'.00000', .true.)

      if (dof_igetNDofGlob(&
          rproblem%roptcontrol%p_rdiscrTargetFlow%RspatialDiscr(1)) .ne. &
          rproblem%roptcontrol%rtargetFlow%RvectorBlock(1)%NEQ) then
        call output_line ('Target flow vector invalid, NEQ wrong!',&
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initTargetFlow')
        call sys_halt()
      end if

      ! Connect the vector to the discretisation structure.
      call lsysbl_enforceStructureDiscr (rproblem%roptcontrol%p_rdiscrTargetFlow,&
          rproblem%roptcontrol%rtargetFlow)
      !rproblem%roptcontrol%rtargetFlow%p_rblockDiscretisation => &
      !    rproblem%roptcontrol%p_rdiscrTargetFlow

      ! How many timesteps does the target flow have?      
      if (rproblem%roptcontrol%itargetFlowTimesteps .eq. -1) then
        ! As many files as we have timesteps
        nactTimesteps = ntimesteps
      else
        ! Exactly itargetFlowDelta timesteps, the values inbetween
        ! must be calculated by interpolation.
        nactTimesteps = rproblem%roptcontrol%itargetFlowTimesteps
      end if
      
      ! Create a time discretisation structure for the target flow and
      ! create the target flow vector. We assume the target flow to be given
      ! in the endpoints of the subintervals -- which corresponds to the
      ! implicit Euler scheme.
      call tdiscr_initTheta (dstartTime, dendTime, nacttimesteps, 1.0_DP, &
          rproblem%roptcontrol%rtargetTimeDiscr)
          
      call sptivec_initVectorDiscr (rproblem%roptcontrol%rtargetFlowNonstat,&
          rproblem%roptcontrol%rtargetTimeDiscr,&
          rproblem%roptcontrol%p_rdiscrTargetFlow)

      call sptivec_loadFromFileSequence (&
          rproblem%roptcontrol%rtargetFlowNonstat,&
          '('''//trim(rproblem%roptcontrol%stargetFlow)//'.'',I5.5)',&
          0,nacttimesteps,&
          rproblem%roptcontrol%itargetFlowDelta,.true.,.true.)

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneTargetFlow (rproblem)
  
!<description>
  ! Cleans up the structure for the optimal control problem. 
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The roptcontrol substructure is initialised.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

    ! Release memory
    select case (rproblem%roptcontrol%itypeTargetFlow) 
    case (1,3)
      call lsysbl_releaseVector (rproblem%roptcontrol%rtargetFlow)
    case (2,4)
      call lsysbl_releaseVector (rproblem%roptcontrol%rtargetFlow)
      call sptivec_releaseVector (rproblem%roptcontrol%rtargetFlowNonstat)
    end select

    if ((rproblem%roptcontrol%smeshTargetFlow .eq. '') .and. &
        (rproblem%roptcontrol%ilevelTargetFlow .le. rproblem%NLMAX)) then
      ! Create references to the existing triangulation
      nullify(rproblem%roptcontrol%p_rtriangulationTargetFlow)
      nullify(rproblem%roptcontrol%p_rdiscrTargetFlow)
    else
      ! Release the allocated triangulation/discretisation
      
      call spdiscr_releaseBlockDiscr(rproblem%roptcontrol%p_rdiscrTargetFlow,.true.)
      deallocate(rproblem%roptcontrol%p_rdiscrTargetFlow)
      
      call tria_done(rproblem%roptcontrol%p_rtriangulationTargetFlow)
      deallocate(rproblem%roptcontrol%p_rtriangulationTargetFlow)
    end if

  end subroutine

end module
