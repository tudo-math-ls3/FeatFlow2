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

MODULE cc2dmediumm2discretisation

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE coarsegridcorrection
  USE spdiscprojection
  USE nonlinearsolver
  USE paramlist
  USE stdoperators
  USE vectorio
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback
  USE cc2dmediumm2nonlinearcoreinit
  
  IMPLICIT NONE
  
CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_get1LevelDiscretisation (rparlist,rboundary,rtriangulation,&
      nequations,rdiscretisation)
  
!<description>
  ! Sets up the discretisation structures for the (primal) velocity and pressure
  ! as specified by the parameters in the parameter list rparList.
!</description>

!<input>
  ! The parameter list with the parameters in the DAT file.
  TYPE(t_parlist), INTENT(IN) :: rparlist
  
  ! Number of equations.
  ! =3 for primal or dual discretisation structure.
  ! =6 for coupled primal/dual discretisation.
  INTEGER, INTENT(IN) :: nequations
  
  ! Underlying boundary.
  TYPE(t_boundary), INTENT(IN), TARGET :: rboundary
  
  ! Underlying triangulation
  TYPE(t_triangulation), INTENT(IN), TARGET :: rtriangulation
!</input>

!<inputoutput>
  ! The discretisation structure to be set up.
  TYPE(t_blockDiscretisation), INTENT(INOUT) :: rdiscretisation
!</inputoutput>

!</subroutine>

    INTEGER :: ielementType,icubA,icubB,icubF
    CHARACTER(LEN=SYS_NAMELEN) :: sstr

    ! Get parameters about the discretisation from the parameter list
    CALL parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
                              'iElementType',ielementType,3)

    CALL parlst_getvalue_string (rparlist,'CC-DISCRETISATION',&
                                 'scubStokes',sstr,'')
    IF (sstr .EQ. '') THEN
      CALL parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
                                'icubStokes',icubA,CUB_G2X2)
    ELSE
      icubA = cub_igetID(sstr)
    END IF

    CALL parlst_getvalue_string (rparlist,'CC-DISCRETISATION',&
                                'scubB',sstr,'')
    IF (sstr .EQ. '') THEN
      CALL parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
                                'icubB',icubB,CUB_G2X2)
    ELSE
      icubB = cub_igetID(sstr)
    END IF

    CALL parlst_getvalue_string (rparlist,'CC-DISCRETISATION',&
                                 'scubF',sstr,'')
    IF (sstr .EQ. '') THEN
      CALL parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
                                'icubF',icubF,CUB_G2X2)
    ELSE
      icubF = cub_igetID(sstr)
    END IF

    CALL spdiscr_initBlockDiscr2D (rdiscretisation,nequations,&
                                   rtriangulation, rboundary)

    SELECT CASE (ielementType)
    CASE (0)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      CALL spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscretisation(1), &
                  EL_E031,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscretisation(1)% &
        RelementDistribution(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscretisation(2) = &
      !   rdiscretisation%RspatialDiscretisation(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      CALL spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscretisation(1), &
          rdiscretisation%RspatialDiscretisation(2), .TRUE.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_initDiscr_combined ( &
                  rdiscretisation%RspatialDiscretisation(3), &
                  EL_Q0,EL_E031,icubB, &
                  rtriangulation, rboundary)

    CASE (1)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      CALL spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscretisation(1), &
                  EL_E030,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscretisation(1)% &
        RelementDistribution(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscretisation(2) = &
      !   rdiscretisation%RspatialDiscretisation(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      CALL spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscretisation(1), &
          rdiscretisation%RspatialDiscretisation(2), .TRUE.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_initDiscr_combined ( &
                  rdiscretisation%RspatialDiscretisation(3), &
                  EL_Q0,EL_E030,icubB, &
                  rtriangulation, rboundary)

    CASE (2)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      CALL spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscretisation(1), &
                  EL_EM31,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscretisation(1)% &
        RelementDistribution(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscretisation(2) = &
      !   rdiscretisation%RspatialDiscretisation(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      CALL spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscretisation(1), &
          rdiscretisation%RspatialDiscretisation(2), .TRUE.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_initDiscr_combined ( &
                  rdiscretisation%RspatialDiscretisation(3), &
                  EL_Q0,EL_EM31,icubB, &
                  rtriangulation, rboundary)

    CASE (3)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      CALL spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscretisation(1), &
                  EL_EM30,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscretisation(1)% &
        RelementDistribution(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscretisation(2) = &
      !   rdiscretisation%RspatialDiscretisation(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      CALL spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscretisation(1), &
          rdiscretisation%RspatialDiscretisation(2), .TRUE.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_initDiscr_combined ( &
                  rdiscretisation%RspatialDiscretisation(3), &
                  EL_Q0,EL_EM30,icubB, &
                  rtriangulation, rboundary)

    CASE (4)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      CALL spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscretisation(1), &
                  EL_Q2,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscretisation(1)% &
        RelementDistribution(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscretisation(2) = &
      !   rdiscretisation%RspatialDiscretisation(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      CALL spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscretisation(1), &
          rdiscretisation%RspatialDiscretisation(2), .TRUE.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_initDiscr_combined ( &
                  rdiscretisation%RspatialDiscretisation(3), &
                  EL_QP1,EL_Q2,icubB, &
                  rtriangulation, rboundary)
                  
    CASE (5)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      CALL spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscretisation(1), &
                  EL_EM30_UNPIVOTED,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscretisation(1)% &
        RelementDistribution(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscretisation(2) = &
      !   rdiscretisation%RspatialDiscretisation(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      CALL spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscretisation(1), &
          rdiscretisation%RspatialDiscretisation(2), .TRUE.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_initDiscr_combined ( &
                  rdiscretisation%RspatialDiscretisation(3), &
                  EL_Q0,EL_EM30_UNPIVOTED,icubB, &
                  rtriangulation, rboundary)

    CASE (6)
      ! rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the 
      ! velocity...
      CALL spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscretisation(1), &
                  EL_EM30_UNSCALED,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscretisation(1)% &
        RelementDistribution(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscretisation(2) = &
      !   rdiscretisation%RspatialDiscretisation(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      CALL spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscretisation(1), &
          rdiscretisation%RspatialDiscretisation(2), .TRUE.)
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_initDiscr_combined ( &
                  rdiscretisation%RspatialDiscretisation(3), &
                  EL_Q0,EL_EM30_UNSCALED,icubB, &
                  rtriangulation, rboundary)

    CASE DEFAULT
      PRINT *,'Unknown discretisation: iElementType = ',ielementType
      STOP
    END SELECT
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: I,k
  
  ! Number of equations in our problem. 
  ! velocity+velocity+pressure + dual velocity+dual velocity+dual pressure = 6
  INTEGER, PARAMETER :: nequations = 6
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! An object for the block discretisation on one level
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisationMass
    
    INTEGER :: icubM
    CHARACTER(LEN=SYS_NAMELEN) :: sstr
    
    ! Set up discrezisation structures on all levels:

    DO i=rproblem%NLMIN,rproblem%NLMAX
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies 3 blocks in the
      ! solution vector.
      ALLOCATE(p_rdiscretisation)
      CALL cc_get1LevelDiscretisation (rproblem%rparamList,&
          rproblem%p_rboundary,p_rtriangulation,nequations,p_rdiscretisation)

      ! Save the discretisation structure to our local LevelInfo structure
      ! for later use.
      rproblem%RlevelInfo(i)%p_rdiscretisation => p_rdiscretisation

      ! -----------------------------------------------------------------------
      ! Optimal control extension
      ! -----------------------------------------------------------------------
      ! The variables 4,5,6 are discretised like the variables 1,2,3
      ! in the discretisation structure.
      ! Variable 1 = x-velocity  -->  Variable 4 = dual x-velocity
      ! Variable 2 = y-velocity  -->  Variable 5 = dual y-velocity
      ! Variable 3 = pressure    -->  Variable 6 = dual pressure
      CALL spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscretisation(1),&
                                     p_rdiscretisation%RspatialDiscretisation(4))
      CALL spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscretisation(2),&
                                     p_rdiscretisation%RspatialDiscretisation(5))
      CALL spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscretisation(3),&
                                     p_rdiscretisation%RspatialDiscretisation(6))
          
      ! -----------------------------------------------------------------------
      ! Separated discretisation structures for primal and dual problem
      ! -----------------------------------------------------------------------
      ! Create a separate block discretisation structure, only for the primal
      ! space. Make a copy of the one we have, i.e. reuse it.
      ! Mark the substructures as being a copy from another discretisation structure,
      ! sharing all information with part 1..3 of the global problem.
      ALLOCATE (rproblem%RlevelInfo(i)%p_rdiscretisationPrimal)
      CALL spdiscr_deriveBlockDiscr (rproblem%RlevelInfo(i)%p_rdiscretisation, &
          rproblem%RlevelInfo(i)%p_rdiscretisationPrimal, 1,3)
      
      ! -----------------------------------------------------------------------
      ! Mass matrices
      ! -----------------------------------------------------------------------

      ! Initialise a discretisation structure for the mass matrix.
      ! Copy the discretisation structure of the first (Stokes) block
      ! and replace the cubature-formula identifier by that which is to be
      ! used for the mass matrix.
      ALLOCATE(p_rdiscretisationMass)
      CALL spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscretisation(1),&
          p_rdiscretisationMass,.TRUE.)
      
      CALL parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                  'scubStokes',sstr,'')
      IF (sstr .EQ. '') THEN
        CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                  'icubStokes',icubM,CUB_G2X2)
      ELSE
        icubM = cub_igetID(sstr)
      END IF

      ! Initialise the cubature formula appropriately.
      DO k = 1,p_rdiscretisationMass%inumFESpaces
        p_rdiscretisationMass%RelementDistribution(k)%ccubTypeBilForm = icubM
      END DO

      rproblem%RlevelInfo(i)%p_rdiscretisationMass => p_rdiscretisationMass
        
    END DO
                                   
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    DO i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Remove the main block discretisation structure. 
      CALL spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%p_rdiscretisation,.TRUE.)
      
      ! Remove the discretisation from the heap.
      DEALLOCATE(rproblem%RlevelInfo(i)%p_rdiscretisation)

      ! Release the block discretisation structures of the primal and dual
      ! space.
      CALL spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%p_rdiscretisationPrimal,.TRUE.)
      DEALLOCATE(rproblem%RlevelInfo(i)%p_rdiscretisationPrimal)
      
      ! -----------------------------------------------------------------------
      ! Mass matrix problem
      ! -----------------------------------------------------------------------

      ! Release the mass matrix discretisation.
      CALL spdiscr_releaseDiscr(rproblem%RlevelInfo(i)%p_rdiscretisationMass)
      DEALLOCATE(rproblem%RlevelInfo(i)%p_rdiscretisationMass)

    END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_allocMatVec (rproblem,rvector,rrhs)
  
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
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! A vector structure for the solution vector. The structure is initialised,
  ! memory is allocated for the data entries.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector

  ! A vector structure for the RHS vector. The structure is initialised,
  ! memory is allocated for the data entries.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,cmatBuildType
  
    ! A pointer to the system matrix and the RHS/solution vectors.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes
    TYPE(t_matrixScalar), POINTER :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient
    TYPE(t_vectorBlock), POINTER :: p_rtempVector

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
  
    ! When the jump stabilisation is used, we have to create an extended
    ! matrix stencil!
    cmatBuildType = BILF_MATC_ELEMENTBASED
    
    CALL parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
                              'IUPWIND1', i)
    IF (i .EQ. 2) cmatBuildType = BILF_MATC_EDGEBASED

    CALL parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
                              'IUPWIND2', i)
    IF (i .EQ. 2) cmatBuildType = BILF_MATC_EDGEBASED
  
    ! Initialise all levels...
    DO i=rproblem%NLMIN,rproblem%NLMAX

      ! -----------------------------------------------------------------------
      ! Basic (Navier-) Stokes problem
      ! -----------------------------------------------------------------------

      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
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
      CALL bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscretisation(1),LSYSSC_MATRIX9,&
                p_rmatrixTemplateFEM,cmatBuildType)

      ! In the global system, there are two gradient matrices B1 and B2.
      ! Create a template matrix that defines their structure.
      p_rmatrixTemplateGradient => rproblem%RlevelInfo(i)%rmatrixTemplateGradient
      
      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscretisation.
      CALL bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscretisation(3),LSYSSC_MATRIX9,&
                p_rmatrixTemplateGradient)
      
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
      CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  p_rmatrixStokes,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      
      ! Allocate memory for the entries; don't initialise the memory.
      CALL lsyssc_allocEmptyMatrix (p_rmatrixStokes,LSYSSC_SETM_UNDEFINED)
      
      ! In the global system, there are two coupling matrices B1 and B2.
      ! Both have the structure of the template gradient matrix.
      ! So connect the two B-matrices to the template gradient matrix
      ! such that they share the same structure.
      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscretisation.
      !
      ! Don't create a content array yet, it will be created by 
      ! the assembly routines later.
      
      CALL lsyssc_duplicateMatrix (p_rmatrixTemplateGradient,&
                  rproblem%RlevelInfo(i)%rmatrixB1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
      CALL lsyssc_duplicateMatrix (p_rmatrixTemplateGradient,&
                  rproblem%RlevelInfo(i)%rmatrixB2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
      ! Allocate memory for the entries; don't initialise the memory.
      CALL lsyssc_allocEmptyMatrix (rproblem%RlevelInfo(i)%rmatrixB1,&
                                           LSYSSC_SETM_UNDEFINED)
      CALL lsyssc_allocEmptyMatrix (rproblem%RlevelInfo(i)%rmatrixB2,&
                                           LSYSSC_SETM_UNDEFINED)
                                           
      ! -----------------------------------------------------------------------
      ! Allocate memory for an identity matrix in the size of the pressure.
      ! Attach a discretisation structure that describes the pressure element
      ! as trial space.
      CALL lsyssc_createDiagMatrixStruc (&
          rproblem%RlevelInfo(i)%rmatrixIdentityPressure,&
          rproblem%RlevelInfo(i)%rmatrixB1%NCOLS,LSYSSC_MATRIX9)
      CALL lsyssc_allocEmptyMatrix (rproblem%RlevelInfo(i)%rmatrixIdentityPressure,&
          LSYSSC_SETM_UNDEFINED)
      rproblem%RlevelInfo(i)%rmatrixIdentityPressure%p_rspatialDiscretisation => &
        p_rdiscretisation%RspatialDiscretisation(3)

      ! -----------------------------------------------------------------------
      ! Temporary vectors
      !
      ! Now on all levels except for the maximum one, create a temporary 
      ! vector on that level, based on the matrix template.
      ! It's used for building the matrices on lower levels.
      IF (i .LT. rproblem%NLMAX) THEN
        p_rtempVector => rproblem%RlevelInfo(i)%rtempVector
        CALL lsysbl_createVecBlockByDiscr (&
            rproblem%RlevelInfo(i)%p_rdiscretisation,p_rtempVector,.TRUE.)
            
        ! The temp vectors for the primal and dual system share their memory
        ! with that temp vector.
        CALL lsysbl_deriveSubvector(p_rtempVector,&
            rproblem%RlevelInfo(i)%rtempVectorPrimal,1,3,.TRUE.)
        CALL lsysbl_deriveSubvector(p_rtempVector,&
            rproblem%RlevelInfo(i)%rtempVectorDual,4,6,.TRUE.)
      END IF
      
    END DO
    
    ! (Only) on the finest level, we need to allocate a RHS vector
    ! and a solution vector.
    !
    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template.
    ! Initialise the vectors with 0.
    CALL lsysbl_createVecBlockByDiscr (&
        rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscretisation,rrhs,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (&
        rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscretisation,rvector,.TRUE.)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_generateStaticMatrices (rproblem,rlevelInfo)
  
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
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! A level-info structure. The static matrices in this structure are generated.
  TYPE(t_problem_lvl), INTENT(INOUT),TARGET :: rlevelInfo
!</inputoutput>

!</subroutine>

    ! A pointer to the Stokes and mass matrix
    TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes,p_rmatrixMass

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    CALL cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! -----------------------------------------------------------------------
    ! Basic (Navier-) Stokes problem
    ! -----------------------------------------------------------------------
  
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%p_rdiscretisation
    
    ! Get a pointer to the (scalar) Stokes matrix:
    p_rmatrixStokes => rlevelInfo%rmatrixStokes
    
    CALL stdop_assembleLaplaceMatrix (p_rmatrixStokes,.TRUE.,rproblem%dnu)
    
    ! In the global system, there are two coupling matrices B1 and B2.
    !
    ! Build the first pressure matrix B1.
    CALL stdop_assembleSimpleMatrix (rlevelInfo%rmatrixB1,&
        DER_FUNC,DER_DERIV_X,-1.0_DP)

    ! Build the second pressure matrix B2.
    CALL stdop_assembleSimpleMatrix (rlevelInfo%rmatrixB2,&
        DER_FUNC,DER_DERIV_Y,-1.0_DP)
                                
    ! -----------------------------------------------------------------------
    ! Mass matrices
    ! -----------------------------------------------------------------------

    p_rmatrixMass => rlevelInfo%rmatrixMass

    ! If there is an existing mass matrix, release it.
    CALL lsyssc_releaseMatrix (rlevelInfo%rmatrixMass)

    ! Generate mass matrix. The matrix has basically the same structure as
    ! our template FEM matrix, so we can take that.
    CALL lsyssc_duplicateMatrix (rlevelInfo%rmatrixTemplateFEM,&
                p_rmatrixMass,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
    ! Change the discretisation structure of the mass matrix to the
    ! correct one; at the moment it points to the discretisation structure
    ! of the Stokes matrix...
    p_rmatrixMass%p_rspatialDiscretisation => &
      rlevelInfo%p_rdiscretisationMass

    ! Call the standard matrix setup routine to build the matrix.                    
    CALL stdop_assembleSimpleMatrix (p_rmatrixMass,DER_FUNC,DER_FUNC)

    ! Clean up the collection (as we are done with the assembly, that's it.
    CALL cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

    ! -----------------------------------------------------------------------
    ! Initialise the identity matrix
    CALL lsyssc_initialiseIdentityMatrix (rlevelInfo%rmatrixIdentityPressure)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_generateBasicMatrices (rproblem)
  
!<description>
  ! Calculates the entries of all static matrices (Mass, B,...) on all levels.
  !
  ! Memory for those matrices must have been allocated before with 
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i

    DO i=rproblem%NLMIN,rproblem%NLMAX
      CALL cc_generateStaticMatrices (&
          rproblem,rproblem%RlevelInfo(i))
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_generateBasicRHS (rproblem,rrhs)
  
!<description>
  ! Calculates the entries of the basic right-hand-side vector on the finest
  ! level. Boundary conditions or similar things are not implemented into 
  ! the vector.
  ! Memory for the RHS vector must have been allocated in advance.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! The RHS vector which is to be filled with data.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_linearForm) :: rlinform
    
    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Drhs
    
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rrhs,p_Drhs)

    ! Get a pointer to the RHS on the finest level as well as to the
    ! block discretisation structure:
    p_rdiscretisation => rrhs%p_rblockDiscretisation
    
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
    CALL cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! Discretise the X-velocity part:
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              rrhs%RvectorBlock(1),coeff_RHS_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(2),rlinform,.TRUE.,&
              rrhs%RvectorBlock(2),coeff_RHS_y,&
              rproblem%rcollection)
                                
    ! The third subvector must be zero initially - as it represents the RHS of
    ! the equation "div(u) = 0".
    CALL lsyssc_clearVector(rrhs%RvectorBlock(3))
    
    ! The RHS terms for the dual equation are calculated similarly using
    ! the desired 'target' flow field.
    !
    ! Discretise the X-velocity part:
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              rrhs%RvectorBlock(4),coeff_TARGET_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(2),rlinform,.TRUE.,&
              rrhs%RvectorBlock(5),coeff_TARGET_y,&
              rproblem%rcollection)
      
    ! Depending on the formulation, to get a reference dual velocity,
    ! it might be necessary to switch the sign of the target velocity field 
    ! because the RHS of the dual equation is '-z'!
    ! Remember that it this case the signs of the mass matrices that couple
    ! primal and dual velocity must be changed, too!
    
    IF (rproblem%roptcontrol%ispaceTimeFormulation .EQ. 0) THEN
      CALL lsyssc_scaleVector (rrhs%RvectorBlock(4),-1.0_DP)
      CALL lsyssc_scaleVector (rrhs%RvectorBlock(5),-1.0_DP)
    END IF

    ! Dual pressure RHS is =0.
    CALL lsyssc_clearVector(rrhs%RvectorBlock(6))
                                
    ! Clean up the collection (as we are done with the assembly, that's it.
    CALL cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initInitialSolution (rproblem,rvector)
  
!<description>
  ! Initialises the initial solution vector into rvector. Depending on the settings
  ! in the DAT file this is either zero or read from a file.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(IN) :: rproblem
!</input>

!<inputoutput>
  ! The solution vector to be initialised. Must be set up according to the
  ! maximum level NLMAX in rproblem!
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(I32) :: istart
    TYPE(t_vectorBlock) :: rvector1,rvector2
    TYPE(t_vectorScalar) :: rvectorTemp
    CHARACTER(LEN=SYS_STRLEN) :: sarray,sfile,sfileString
    INTEGER :: ilev
    INTEGER(PREC_VECIDX) :: NEQ
    TYPE(t_interlevelProjectionBlock) :: rprojection 

    ! Get the parameter what to do with rvector
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'isolutionStart',istart,0)
    CALL parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'ssolutionStart',sfileString,'')

    ! Create a temp vector at level NLMAX-istart+1.
    ilev = rproblem%NLMAX-ABS(istart)+1
    
    IF (ilev .LT. rproblem%NLMIN) THEN
      PRINT *,'Warning: Level of start vector is < NLMIN! Initialising with zero!'
      istart = 0
    END IF
    
    ! Inít with zero? Read from file?
    IF (istart .EQ. 0) THEN
      ! Init with zero
      CALL lsysbl_clearVector (rvector)
    ELSE
      ! Remove possible ''-characters
      READ(sfileString,*) sfile

      CALL lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev)%p_rdiscretisation,rvector1,.FALSE.)
      
      ! Read in the vector
      CALL vecio_readBlockVectorHR (&
          rvector1, sarray, .TRUE., 0, sfile, istart .GT. 0)
          
      ! If the vector is on level < NLMAX, we have to bring it to level NLMAX
      DO WHILE (ilev .LT. rproblem%NLMAX)
        
        ! Initialise a vector for the higher level and a prolongation structure.
        CALL lsysbl_createVectorBlock (&
            rproblem%RlevelInfo(ilev+1)%p_rdiscretisation,rvector2,.FALSE.)
        
        CALL mlprj_initProjectionVec (rprojection,rvector2)
        
        ! Prolongate to the next higher level.

        NEQ = mlprj_getTempMemoryVec (rprojection,rvector1,rvector2)
        IF (NEQ .NE. 0) CALL lsyssc_createVector (rvectorTemp,NEQ,.FALSE.)
        CALL mlprj_performProlongation (rprojection,rvector1, &
                                        rvector2,rvectorTemp)
        IF (NEQ .NE. 0) CALL lsyssc_releaseVector (rvectorTemp)
        
        ! Swap rvector1 and rvector2. Release the coarse grid vector.
        CALL lsysbl_swapVectors (rvector1,rvector2)
        CALL lsysbl_releaseVector (rvector2)
        
        CALL mlprj_doneProjection (rprojection)
        
        ! rvector1 is now on level ilev+1
        ilev = ilev+1
        
      END DO
      
      ! Copy the resulting vector rvector1 to the output.
      CALL lsysbl_copyVector (rvector1,rvector)
      
      ! Release the temp vector
      CALL lsysbl_releaseVector (rvector1)
      
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_writeSolution (rproblem,rvector)
  
!<description>
  ! Writes a solution vector rvector to a file as configured in the parameters
  ! in the DAT file.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(IN) :: rproblem

  ! The solution vector to be written out. Must be set up according to the
  ! maximum level NLMAX in rproblem!
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
!</input>

!</subroutine>

    ! local variables
    INTEGER(I32) :: idestLevel
    TYPE(t_vectorBlock) :: rvector1,rvector2
    TYPE(t_vectorScalar) :: rvectorTemp
    CHARACTER(LEN=SYS_STRLEN) :: sfile,sfileString
    INTEGER :: ilev
    INTEGER(PREC_VECIDX) :: NEQ
    TYPE(t_interlevelProjectionBlock) :: rprojection 
    LOGICAL :: bformatted

    ! Get the parameter what to do with rvector
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iSolutionWrite',idestLevel,0)
    CALL parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'sSolutionWrite',sfileString,'')

    IF (idestLevel .EQ. 0) RETURN ! nothing to do.
    
    ! Remove possible ''-characters
    READ(sfileString,*) sfile
    
    bformatted = idestLevel .GT. 0
    idestLevel = rproblem%NLMAX-ABS(idestLevel)+1 ! level where to write out

    IF (idestLevel .LT. rproblem%NLMIN) THEN
      PRINT *,'Warning: Level for solution vector is < NLMIN! &
              &Writing out at level NLMIN!'
      idestLevel = rproblem%NLMIN
    END IF
    
    ! Interpolate the solution down to level istart.
    CALL lsysbl_copyVector (rvector,rvector1)   ! creates new rvector1!

    DO ilev = rproblem%NLMAX,idestLevel+1,-1
      
      ! Initialise a vector for the lower level and a prolongation structure.
      CALL lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev-1)%p_rdiscretisation,rvector2,.FALSE.)
      
      CALL mlprj_initProjectionVec (rprojection,rvector2)
      
      ! Interpolate to the next higher level.
      ! (Don't 'restrict'! Restriction would be for the dual space = RHS vectors!)

      NEQ = mlprj_getTempMemoryVec (rprojection,rvector2,rvector1)
      IF (NEQ .NE. 0) CALL lsyssc_createVector (rvectorTemp,NEQ,.FALSE.)
      CALL mlprj_performInterpolation (rprojection,rvector2,rvector1, &
                                       rvectorTemp)
      IF (NEQ .NE. 0) CALL lsyssc_releaseVector (rvectorTemp)
      
      ! Swap rvector1 and rvector2. Release the fine grid vector.
      CALL lsysbl_swapVectors (rvector1,rvector2)
      CALL lsysbl_releaseVector (rvector2)
      
      CALL mlprj_doneProjection (rprojection)
      
    END DO

    ! Write out the solution.
    IF (bformatted) THEN
      CALL vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .TRUE.,&
                                    0, sfile, '(E22.15)')
    ELSE
      CALL vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .TRUE.,0, sfile)
    END IF

    ! Release temp memory.
    CALL lsysbl_releaseVector (rvector1)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_doneMatVec (rproblem,rvector,rrhs)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! A vector structure for the solution vector. The structure is cleaned up,
  ! memory is released.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector

  ! A vector structure for the RHS vector. The structure is cleaned up,
  ! memory is released.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

    INTEGER :: i

    ! Release matrices and vectors on all levels
    DO i=rproblem%NLMAX,rproblem%NLMIN,-1

      ! If there is an existing mass matrix, release it.
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixMass)

      ! Release Stokes, B1 and B2 matrix
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB2)
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB1)
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixStokes)
      
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixIdentityPressure)
      
      ! Release the template matrices. This is the point, where the
      ! memory of the matrix structure is released.
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixTemplateGradient)
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixTemplateFEM)
      
      ! Remove the temp vector that was used for interpolating the solution
      ! from higher to lower levels in the nonlinear iteration.
      IF (i .LT. rproblem%NLMAX) THEN
        CALL lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVectorPrimal)
        CALL lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVectorDual)
        CALL lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector)
      END IF
      
    END DO

    ! Delete solution/RHS vector
    CALL lsysbl_releaseVector (rvector)
    CALL lsysbl_releaseVector (rrhs)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initOptControl (rproblem)
  
!<description>
  ! Initialises the optimal control problem. Reads parameters from the DAT
  ! file and evaluates them.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The roptcontrol substructure is initialised.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    CHARACTER(SYS_STRLEN) :: spar

    ! Read in the parameters
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dalphaC',rproblem%roptcontrol%dalphaC,1.0_DP)
        
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
        'dgammaC',rproblem%roptcontrol%dgammaC,0.0_DP)
        
    CALL parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itypeTargetFlow',rproblem%roptcontrol%itypeTargetFlow,0)

    CALL parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ilevelTargetFlow',rproblem%roptcontrol%ilevelTargetFlow,0)

    CALL parlst_getvalue_string (rproblem%rparamList,'OPTIMALCONTROL',&
        'stargetFlow',spar,'')
    READ(spar,*) rproblem%roptcontrol%stargetFlow

    CALL parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itargetFlowDelta',rproblem%roptcontrol%itargetFlowDelta,1)

    CALL parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itargetFlowTimesteps',rproblem%roptcontrol%itargetFlowTimesteps,-1)
    
    CALL parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'ispaceTimeFormulation',rproblem%roptcontrol%ispaceTimeFormulation,0)
    
    CALL parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'iconvectionExplicit',rproblem%roptcontrol%iconvectionExplicit,0)
    
    CALL parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itypeTerminalCondition',rproblem%roptcontrol%itypeTerminalCondition,0)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_doneOptControl (rproblem)
  
!<description>
  ! Cleans up the structure for the optimal control problem. 
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The roptcontrol substructure is initialised.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    rproblem%roptcontrol%dalphaC = 1.0_DP
    rproblem%roptcontrol%dgammaC = 0.0_DP
    rproblem%roptcontrol%itypeTargetFlow = 0
    rproblem%roptcontrol%stargetFlow = ''

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initTargetFlow (rproblem,dstartTime,dendTime,ntimesteps)
  
!<description>
  ! Reads in the target flow of the optimal control problem.
!</description>

!<input>
  ! OPTIONAL: Start time in a nonstationary simulation.
  ! This parameter must be present in a nonstatinoary simulation. In a stationary
  ! simulation, the parameter can be skipped.
  REAL(DP), INTENT(IN), OPTIONAL :: dstartTime

  ! OPTIONAL: Start time in a nonstationary simulation.
  ! This parameter must be present in a nonstatinoary simulation. In a stationary
  ! simulation, the parameter can be skipped.
  REAL(DP), INTENT(IN), OPTIONAL :: dendTime

  ! OPTIONAL: Total number of timesteps in a nonstationary simulation.
  ! This parameter must be present in a nonstatinoary simulation. In a stationary
  ! simulation, the parameter can be skipped.
  INTEGER, INTENT(IN), OPTIONAL :: ntimesteps
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The target flow is saved to the roptcontrol substructure.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>
    
    ! local variables
    CHARACTER(SYS_STRLEN) :: sarray
    INTEGER :: iref
    INTEGER :: nactTimesteps
    
    ! At first, what is with the triangulation and discretisation?
    ! Is the triangulation inside of the bounds specified by NLMIN/NLMAX
    ! or not? If yes, we can use a predefined triangulation.
    ! If not, we have to refine the mesh until we reach the destination
    ! level.
    IF (rproblem%roptcontrol%ilevelTargetFlow .LE. rproblem%NLMAX) THEN
      ! Create references to the existing triangulation
      rproblem%roptcontrol%p_rtriangulation => &
        rproblem%RlevelInfo(rproblem%roptcontrol%ilevelTargetFlow)%rtriangulation
        
      rproblem%roptcontrol%p_rdiscrTargetFlow => &
        rproblem%RlevelInfo(rproblem%roptcontrol%ilevelTargetFlow)%p_rdiscretisationPrimal
    ELSE
      iref = rproblem%roptcontrol%ilevelTargetFlow - rproblem%NLMAX
      ALLOCATE(rproblem%roptcontrol%p_rtriangulation)
      CALL tria_duplicate (&
          rproblem%RlevelInfo(rproblem%NLMAX)%rtriangulation,&
          rproblem%roptcontrol%p_rtriangulation,TR_SHARE_ALL)
          
      ! Refine
      CALL tria_quickRefine2LevelOrdering(iref,&
          rproblem%roptcontrol%p_rtriangulation,rproblem%p_rboundary)
      CALL tria_initStandardMeshFromRaw (&
          rproblem%roptcontrol%p_rtriangulation,rproblem%p_rboundary)
          
      ! Create a discretisation structure corresponding to that mesh.
      ALLOCATE(rproblem%roptcontrol%p_rdiscrTargetFlow)
      CALL cc_get1LevelDiscretisation (rproblem%rparamList,rproblem%p_rboundary,&
          rproblem%roptcontrol%p_rtriangulation,NDIM2D+1,&
          rproblem%roptcontrol%p_rdiscrTargetFlow)
      
    END IF

    ! Probably create a vector containing the target flow
    SELECT CASE (rproblem%roptcontrol%itypeTargetFlow)
    CASE (1,3)
    
      ! Stationary target flow. Read in the vector
      CALL vecio_readBlockVectorHR (&
          rproblem%roptcontrol%rtargetFlow, sarray, .TRUE., &
          0, rproblem%roptcontrol%stargetFlow, .TRUE.)
          
      ! Connect the vector to the discretisation structure.
      rproblem%roptcontrol%rtargetFlow%p_rblockDiscretisation => &
          rproblem%roptcontrol%p_rdiscrTargetFlow
          
    CASE (2,4)
    
      ! Nonstationary target flow given by a sequence of files.

      IF (.NOT. PRESENT(ntimesteps)) THEN
        CALL output_line ('ntimesteps not available!')
        CALL sys_halt()
      END IF

      ! Read the first vector. It defines the shape of all vectors in rx!
      CALL vecio_readBlockVectorHR (&
          rproblem%roptcontrol%rtargetFlow, sarray, .TRUE., &
          0, TRIM(rproblem%roptcontrol%stargetFlow)//'.00000', .TRUE.)

      IF (dof_igetNDofGlob(&
          rproblem%roptcontrol%p_rdiscrTargetFlow%RspatialDiscretisation(1)) .NE. &
          rproblem%roptcontrol%rtargetFlow%RvectorBlock(1)%NEQ) THEN
        CALL output_line ('Target flow vector invalid, NEQ wrong!',&
            OU_CLASS_ERROR,OU_MODE_STD,'cc_initTargetFlow')
        CALL sys_halt()
      END IF

      ! Connect the vector to the discretisation structure.
      CALL lsysbl_enforceStructureDiscr (rproblem%roptcontrol%p_rdiscrTargetFlow,&
          rproblem%roptcontrol%rtargetFlow)
      !rproblem%roptcontrol%rtargetFlow%p_rblockDiscretisation => &
      !    rproblem%roptcontrol%p_rdiscrTargetFlow

      ! How many timesteps does the target flow have?      
      IF (rproblem%roptcontrol%itargetFlowTimesteps .EQ. -1) THEN
        ! As many files as we have timesteps
        nactTimesteps = ntimesteps
      ELSE
        ! Exactly itargetFlowDelta timesteps, the values inbetween
        ! must be calculated by interpolation.
        nactTimesteps = rproblem%roptcontrol%itargetFlowTimesteps
      END IF
      
      ! Create a time discretisation structure for the target flow and
      ! create the target flow vector. We assume the target flow to be given
      ! in the endpoints of the subintervals -- which corresponds to the
      ! implicit Euler scheme.
      CALL tdiscr_initTheta (dstartTime, dendTime, nacttimesteps, 1.0_DP, &
          rproblem%roptcontrol%rtargetTimeDiscr)
          
      CALL sptivec_initVectorDiscr (rproblem%roptcontrol%rtargetFlowNonstat,&
          rproblem%roptcontrol%rtargetTimeDiscr,&
          rproblem%roptcontrol%p_rdiscrTargetFlow)

      CALL sptivec_loadFromFileSequence (&
          rproblem%roptcontrol%rtargetFlowNonstat,&
          '('''//TRIM(rproblem%roptcontrol%stargetFlow)//'.'',I5.5)',&
          0,nacttimesteps,&
          rproblem%roptcontrol%itargetFlowDelta,.TRUE.,.TRUE.)

    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_doneTargetFlow (rproblem)
  
!<description>
  ! Cleans up the structure for the optimal control problem. 
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  ! The roptcontrol substructure is initialised.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! Release memory
    SELECT CASE (rproblem%roptcontrol%itypeTargetFlow) 
    CASE (1,3)
      CALL lsysbl_releaseVector (rproblem%roptcontrol%rtargetFlow)
    CASE (2,4)
      CALL lsysbl_releaseVector (rproblem%roptcontrol%rtargetFlow)
      CALL sptivec_releaseVector (rproblem%roptcontrol%rtargetFlowNonstat)
    END SELECT

    IF (rproblem%roptcontrol%ilevelTargetFlow .LE. rproblem%NLMAX) THEN
      ! Create references to the existing triangulation
      NULLIFY(rproblem%roptcontrol%p_rtriangulation)
      NULLIFY(rproblem%roptcontrol%p_rdiscrTargetFlow)
    ELSE
      ! Release the allocated triangulation/discretisation
      
      CALL spdiscr_releaseBlockDiscr(rproblem%roptcontrol%p_rdiscrTargetFlow,.TRUE.)
      DEALLOCATE(rproblem%roptcontrol%p_rdiscrTargetFlow)
      
      CALL tria_done(rproblem%roptcontrol%p_rtriangulation)
      DEALLOCATE(rproblem%roptcontrol%p_rtriangulation)
    END IF

  END SUBROUTINE

END MODULE
