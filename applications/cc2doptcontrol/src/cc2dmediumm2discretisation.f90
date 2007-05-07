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
!# 1.) c2d2_initDiscretisation
!#     -> Initialise the discretisation structure inside of the problem
!#        structure using the parameters from the INI/DAT files.
!#
!# 2.) c2d2_allocMatVec
!#     -> Allocates memory for vectors/matrices
!#
!# 3.) c2d2_generateStaticMatrices
!#     -> Assembles matrix entries of static matrices (Stokes, B)
!#
!# 4.) c2d2_generateStaticSystemParts
!#     -> Initialises static parts of the global system matrix with
!#        entries of template (B-) matrices
!#
!# 5.) c2d2_generateBasicRHS
!#     -> Generates a general RHS vector without any boundary conditions
!#        implemented
!#
!# 6.) c2d2_doneMatVec
!#     -> Cleanup of matrices/vectors, release all memory
!#
!# 7.) c2d2_doneDiscretisation
!#     -> Cleanup of the underlying discretisation structures
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
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback
  
  IMPLICIT NONE
  
CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initDiscretisation (rproblem)
  
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
  INTEGER :: I,k,ielementType,icubA,icubB,icubF, icubM
  CHARACTER(LEN=SYS_NAMELEN) :: sstr
  
  ! Number of equations in our problem. 
  ! velocity+velocity+pressure + dual velocity+dual velocity+dual pressure = 6
  INTEGER, PARAMETER :: nequations = 6
  
    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! An object for the block discretisation on one level
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisationMass
    
    ! Which discretisation is to use?
    ! Which cubature formula should be used?
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iElementType',ielementType,3)

    CALL parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'scubStokes',sstr,'')
    IF (sstr .EQ. '') THEN
      CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'icubStokes',icubA,CUB_G2X2)
    ELSE
      icubA = cub_igetID(sstr)
    END IF

    CALL parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                'scubB',sstr,'')
    IF (sstr .EQ. '') THEN
      CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'icubB',icubB,CUB_G2X2)
    ELSE
      icubB = cub_igetID(sstr)
    END IF

    CALL parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'scubF',sstr,'')
    IF (sstr .EQ. '') THEN
      CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'icubF',icubF,CUB_G2X2)
    ELSE
      icubF = cub_igetID(sstr)
    END IF

    ! Now set up discrezisation structures on all levels:

    DO i=rproblem%NLMIN,rproblem%NLMAX
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%p_rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies 3 blocks in the
      ! solution vector.
      ALLOCATE(p_rdiscretisation)
      CALL spdiscr_initBlockDiscr2D (p_rdiscretisation,nequations,&
                                     p_rtriangulation, p_rboundary)

      ! Save the discretisation structure to our local LevelInfo structure
      ! for later use.
      rproblem%RlevelInfo(i)%p_rdiscretisation => p_rdiscretisation

      SELECT CASE (ielementType)
      CASE (0)
        ! p_rdiscretisation%Rdiscretisations is a list of scalar 
        ! discretisation structures for every component of the solution vector.
        ! We have a solution vector with three components:
        !  Component 1 = X-velocity
        !  Component 2 = Y-velocity
        !  Component 3 = Pressure
        ! For simplicity, we set up one discretisation structure for the 
        ! velocity...
        CALL spdiscr_initDiscr_simple ( &
                    p_rdiscretisation%RspatialDiscretisation(1), &
                    EL_E031,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscretisation(1)% &
          RelementDistribution(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards,
        ! but we have to be careful when releasing the discretisation structures
        ! at the end of the program!
        p_rdiscretisation%RspatialDiscretisation(2) = &
          p_rdiscretisation%RspatialDiscretisation(1)
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscretisation(3), &
                    EL_Q0,EL_E031,icubB, &
                    p_rtriangulation, p_rboundary)

      CASE (1)
        ! p_rdiscretisation%Rdiscretisations is a list of scalar 
        ! discretisation structures for every component of the solution vector.
        ! We have a solution vector with three components:
        !  Component 1 = X-velocity
        !  Component 2 = Y-velocity
        !  Component 3 = Pressure
        ! For simplicity, we set up one discretisation structure for the 
        ! velocity...
        CALL spdiscr_initDiscr_simple ( &
                    p_rdiscretisation%RspatialDiscretisation(1), &
                    EL_E030,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscretisation(1)% &
          RelementDistribution(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards,
        ! but we have to be careful when releasing the discretisation structures
        ! at the end of the program!
        p_rdiscretisation%RspatialDiscretisation(2) = &
          p_rdiscretisation%RspatialDiscretisation(1)
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscretisation(3), &
                    EL_Q0,EL_E030,icubB, &
                    p_rtriangulation, p_rboundary)

      CASE (2)
        ! p_rdiscretisation%Rdiscretisations is a list of scalar 
        ! discretisation structures for every component of the solution vector.
        ! We have a solution vector with three components:
        !  Component 1 = X-velocity
        !  Component 2 = Y-velocity
        !  Component 3 = Pressure
        ! For simplicity, we set up one discretisation structure for the 
        ! velocity...
        CALL spdiscr_initDiscr_simple ( &
                    p_rdiscretisation%RspatialDiscretisation(1), &
                    EL_EM31,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscretisation(1)% &
          RelementDistribution(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards,
        ! but we have to be careful when releasing the discretisation structures
        ! at the end of the program!
        p_rdiscretisation%RspatialDiscretisation(2) = &
          p_rdiscretisation%RspatialDiscretisation(1)
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscretisation(3), &
                    EL_Q0,EL_EM31,icubB, &
                    p_rtriangulation, p_rboundary)

      CASE (3)
        ! p_rdiscretisation%Rdiscretisations is a list of scalar 
        ! discretisation structures for every component of the solution vector.
        ! We have a solution vector with three components:
        !  Component 1 = X-velocity
        !  Component 2 = Y-velocity
        !  Component 3 = Pressure
        ! For simplicity, we set up one discretisation structure for the 
        ! velocity...
        CALL spdiscr_initDiscr_simple ( &
                    p_rdiscretisation%RspatialDiscretisation(1), &
                    EL_EM30,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscretisation(1)% &
          RelementDistribution(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards,
        ! but we have to be careful when releasing the discretisation structures
        ! at the end of the program!
        p_rdiscretisation%RspatialDiscretisation(2) = &
          p_rdiscretisation%RspatialDiscretisation(1)
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscretisation(3), &
                    EL_Q0,EL_EM30,icubB, &
                    p_rtriangulation, p_rboundary)

      CASE (4)
        ! p_rdiscretisation%Rdiscretisations is a list of scalar 
        ! discretisation structures for every component of the solution vector.
        ! We have a solution vector with three components:
        !  Component 1 = X-velocity
        !  Component 2 = Y-velocity
        !  Component 3 = Pressure
        ! For simplicity, we set up one discretisation structure for the 
        ! velocity...
        CALL spdiscr_initDiscr_simple ( &
                    p_rdiscretisation%RspatialDiscretisation(1), &
                    EL_Q2,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscretisation(1)% &
          RelementDistribution(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards,
        ! but we have to be careful when releasing the discretisation structures
        ! at the end of the program!
        p_rdiscretisation%RspatialDiscretisation(2) = &
          p_rdiscretisation%RspatialDiscretisation(1)
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscretisation(3), &
                    EL_QP1,EL_Q2,icubB, &
                    p_rtriangulation, p_rboundary)
                    
      CASE (5)
        ! p_rdiscretisation%Rdiscretisations is a list of scalar 
        ! discretisation structures for every component of the solution vector.
        ! We have a solution vector with three components:
        !  Component 1 = X-velocity
        !  Component 2 = Y-velocity
        !  Component 3 = Pressure
        ! For simplicity, we set up one discretisation structure for the 
        ! velocity...
        CALL spdiscr_initDiscr_simple ( &
                    p_rdiscretisation%RspatialDiscretisation(1), &
                    EL_EM30_UNPIVOTED,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscretisation(1)% &
          RelementDistribution(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards,
        ! but we have to be careful when releasing the discretisation structures
        ! at the end of the program!
        p_rdiscretisation%RspatialDiscretisation(2) = &
          p_rdiscretisation%RspatialDiscretisation(1)
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscretisation(3), &
                    EL_Q0,EL_EM30_UNPIVOTED,icubB, &
                    p_rtriangulation, p_rboundary)

      CASE (6)
        ! p_rdiscretisation%Rdiscretisations is a list of scalar 
        ! discretisation structures for every component of the solution vector.
        ! We have a solution vector with three components:
        !  Component 1 = X-velocity
        !  Component 2 = Y-velocity
        !  Component 3 = Pressure
        ! For simplicity, we set up one discretisation structure for the 
        ! velocity...
        CALL spdiscr_initDiscr_simple ( &
                    p_rdiscretisation%RspatialDiscretisation(1), &
                    EL_EM30_UNSCALED,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscretisation(1)% &
          RelementDistribution(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards,
        ! but we have to be careful when releasing the discretisation structures
        ! at the end of the program!
        p_rdiscretisation%RspatialDiscretisation(2) = &
          p_rdiscretisation%RspatialDiscretisation(1)
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscretisation(3), &
                    EL_Q0,EL_EM30_UNSCALED,icubB, &
                    p_rtriangulation, p_rboundary)

      CASE DEFAULT
        PRINT *,'Unknown discretisation: iElementType = ',ielementType
        STOP
      END SELECT
      
      ! -----------------------------------------------------------------------
      ! Optimal control extension
      ! -----------------------------------------------------------------------
      ! The variables 4,5,6 are discretised like the variables 1,2,3
      ! in the discretisation structure.
      ! Variable 1 = x-velocity  -->  Variable 4 = dual x-velocity
      ! Variable 2 = y-velocity  -->  Variable 5 = dual y-velocity
      ! Variable 3 = pressure    -->  Variable 6 = dual pressure
      ! As we simply copy the structures, we have to be careful when releasing
      ! the structures at the end of the program!
      p_rdiscretisation%RspatialDiscretisation(4) = &
          p_rdiscretisation%RspatialDiscretisation(1)
      p_rdiscretisation%RspatialDiscretisation(5) = &
          p_rdiscretisation%RspatialDiscretisation(2)
      p_rdiscretisation%RspatialDiscretisation(6) = &
          p_rdiscretisation%RspatialDiscretisation(3)

      ! -----------------------------------------------------------------------
      ! Mass matrices
      ! -----------------------------------------------------------------------

      ! Initialise a discretisation structure for the mass matrix.
      ! Copy the discretisation structure of the first (Stokes) block
      ! and replace the cubature-formula identifier by that which is to be
      ! used for the mass matrix.
      ALLOCATE(p_rdiscretisationMass)
      p_rdiscretisationMass = p_rdiscretisation%RspatialDiscretisation(1)
      
      ! Mark the mass matrix discretisation structure as copy of another one -
      ! to prevent accidental deallocation of memory.
      p_rdiscretisationMass%bisCopy = .TRUE.

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

  SUBROUTINE c2d2_doneDiscretisation (rproblem)
  
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
  TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    DO i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Before we remove the block discretisation structure, remember that
      ! we copied the scalar discretisation structure for the X-velocity
      ! to the Y-velocity.
      ! To prevent errors or wrong deallocation, we manually release the
      ! spatial discretisation structures of each of the components.
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation

      ! Remove spatial discretisation structure of the velocity:
      CALL spdiscr_releaseDiscr(p_rdiscretisation%RspatialDiscretisation(1))
      
      ! Don't remove that of the Y-velocity; there is none :)
      !
      ! Remove the discretisation structure of the pressure.
      CALL spdiscr_releaseDiscr(p_rdiscretisation%RspatialDiscretisation(3))
      
      ! Don't remove the discretisation structures for the dual variables;
      ! they had been the same as for the primal variables and were
      ! released already!
      !
      ! Finally remove the block discretisation structure. Don't release
      ! the substructures again.
      CALL spdiscr_releaseBlockDiscr(p_rdiscretisation,.FALSE.)
      
      ! Remove the discretisation from the heap.
      DEALLOCATE(p_rdiscretisation)

      ! -----------------------------------------------------------------------
      ! Mass matrix problem
      ! -----------------------------------------------------------------------

      ! Release the mass matrix discretisation.
      ! Don't release the content as we created it as a copy of the Stokes
      ! discretisation structure.
      DEALLOCATE(rproblem%RlevelInfo(i)%p_rdiscretisationMass)

    END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_allocMatVec (rproblem,rvector,rrhs)
  
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
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes
    TYPE(t_matrixScalar), POINTER :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient
    TYPE(t_vectorBlock), POINTER :: p_rtempVector

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
  
    ! When the jump stabilisation is used, we have to create an extended
    ! matrix stencil!
    cmatBuildType = BILF_MATC_ELEMENTBASED
    
    CALL parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
                              'IUPWIND', i)
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
      ! Now let's come to the main system matrix, which is a block matrix.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)    
      
      ! Let's consider the global system in detail:
      !
      !    ( A11       B1  M/a          ) = ( A11  A12  A13 A14           )
      !    (      A22  B2       M/a     )   ( A21  A22  A23      A25      )
      !    ( B1^T B2^T                  )   ( A31  A32  A33               )
      !    ( -M            A44  A45  B1 )   ( A41           A44  A45  A46 )
      !    (      -M       A54  A55  B2 )   (      A52      A54  A55  A56 )
      !    (               B1^T B2^T    )   (               A64  A65  A66 )
      !

      ! The matrices A11 and A22 of the global system matrix have exactly
      ! the same structure as the original Stokes matrix from above!
      ! Initialise them with the same structure, i.e. A11, A22 and the
      ! Stokes matrix L share(!) the same structure.
      !
      ! For this purpose, use the "duplicate matrix" routine.
      ! The structure of the matrix is shared with the template FEM matrix.
      ! For the content, a new empty array is allocated which will later receive
      ! the entries.
      CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  p_rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        
      IF (.NOT. rproblem%bdecoupledXY) THEN          
        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix 
        ! A22 is identical to A11! So mirror A11 to A22 sharing the
        ! structure and the content.
        CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
                    p_rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      ELSE
        ! Otherwise, create another copy of the template matrix.
        CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    p_rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      END IF

      ! Manually change the discretisation structure of the Y-velocity 
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      p_rmatrix%RmatrixBlock(2,2)%p_rspatialDiscretisation => &
        p_rdiscretisation%RspatialDiscretisation(2)
                                  
      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create empty space for the entries. 
      ! Later, the B-matrices are copied into here and modified for boundary
      ! conditions.
      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(1,3),&
                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(2,3),&
                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      ! These matrices will not change during the whole computation,
      ! so we can put refereces to the original ones to the system matrix.
      CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(3,1),&
                                   LSYSSC_TR_VIRTUAL)

      CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(3,2),&
                                   LSYSSC_TR_VIRTUAL)
                                   
      ! -----------------------------------------------------------------------
      ! Optimal control extension
      ! -----------------------------------------------------------------------
      ! Add submatrices for the dual system to the main system matrix.
      ! The velocity submatrices have the same structure as the template
      ! FEM matrix, but they have all different content!
      
      CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  p_rmatrix%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      IF (.NOT. rproblem%bdecoupledXY) THEN          
        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix 
        ! A55 is identical to A11! So mirror A44 to A55 sharing the
        ! structure and the content.
        CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(4,4),&
                    p_rmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      ELSE
        ! Otherwise, create another copy of the template matrix.
        CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    p_rmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      END IF

      ! For a simple Stokes control equation, we don't need the matrices
      ! (4,5) and (5,4) -- they only appear if the convective operator is present.
      ! Furthermore, the presence of the convective operator forces A(5,5) to be
      ! independent of A(4,4), so overwrite the definition of A(5,5) from above.
      IF (rproblem%iequation .EQ. 0) THEN
        CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    p_rmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    p_rmatrix%RmatrixBlock(4,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                    p_rmatrix%RmatrixBlock(5,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      END IF

      ! Add the B- and B^T submatrices to the main system matrix.
      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(4,6),&
                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(5,6),&
                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
      CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(6,4),&
                                   LSYSSC_TR_VIRTUAL)

      CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(6,5),&
                                   LSYSSC_TR_VIRTUAL)
                                   
      ! Add space for the mass matrices to the system matrix.
      !
      ! They share their structure with the main mass matrix but provide
      ! empty space for the entries. (The entries miht be changed due to
      ! boundary conditions!)
      CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  p_rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  p_rmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  p_rmatrix%RmatrixBlock(4,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      CALL lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  p_rmatrix%RmatrixBlock(5,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
      ! That's it, all submatrices are basically set up.
      !
      ! Update the structural information of the block matrix, as we manually
      ! changed the submatrices:
      CALL lsysbl_updateMatStrucInfo (p_rmatrix)
      
      ! Now on all levels except for the maximum one, create a temporary 
      ! vector on that level, based on the matrix template.
      ! It's used for building the matrices on lower levels.
      IF (i .LT. rproblem%NLMAX) THEN
        p_rtempVector => rproblem%RlevelInfo(i)%rtempVector
        CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rtempVector,.FALSE.)
        
        ! Add the temp vector to the collection on level i
        ! for use in the callback routine
        CALL collct_setvalue_vec(rproblem%rcollection,PAR_TEMPVEC,p_rtempVector,&
                                .TRUE.,i)
      END IF
      

    END DO
    
    ! (Only) on the finest level, we need to have to allocate a RHS vector
    ! and a solution vector.
    !
    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template.
    ! Initialise the vectors with 0.
    CALL lsysbl_createVecBlockIndMat (p_rmatrix,rrhs, .TRUE.)
    CALL lsysbl_createVecBlockIndMat (p_rmatrix,rvector, .TRUE.)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_generateStaticMatrices (rproblem)
  
!<description>
  ! Calculates entries of all static matrices (Stokes, B-matrices,...)
  ! of the problem, i.e. of all matrices that do not change during the 
  ! computation or which serve as template for generating other matrices.
  !
  ! Memory for those matrices must have been allocated before with 
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i
  
    ! A bilinear and linear form describing the analytic problem to solve
    !TYPE(t_bilinearForm) :: rform
    
    ! A pointer to the Stokes and mass matrix
    TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes,p_rmatrixMass

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    CALL c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)

    DO i=rproblem%NLMIN,rproblem%NLMAX
    
      ! -----------------------------------------------------------------------
      ! Basic (Navier-) Stokes problem
      ! -----------------------------------------------------------------------
    
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      ! Get a pointer to the (scalar) Stokes matrix:
      p_rmatrixStokes => rproblem%RlevelInfo(i)%rmatrixStokes
      
      CALL stdop_assembleLaplaceMatrix (p_rmatrixStokes,.TRUE.,rproblem%dnu)
      
      ! In the global system, there are two coupling matrices B1 and B2.
      !
      ! Build the first pressure matrix B1.
      CALL stdop_assembleSimpleMatrix (rproblem%RlevelInfo(i)%rmatrixB1,&
          DER_FUNC,DER_DERIV_X,-1.0_DP)

      ! Build the second pressure matrix B2.
      CALL stdop_assembleSimpleMatrix (rproblem%RlevelInfo(i)%rmatrixB2,&
          DER_FUNC,DER_DERIV_Y,-1.0_DP)
                                  
      ! -----------------------------------------------------------------------
      ! Mass matrices
      ! -----------------------------------------------------------------------

      p_rmatrixMass => rproblem%RlevelInfo(i)%rmatrixMass

      ! If there is an existing mass matrix, release it.
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixMass)

      ! Generate mass matrix. The matrix has basically the same structure as
      ! our template FEM matrix, so we can take that.
      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixTemplateFEM,&
                  p_rmatrixMass,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
      ! Change the discretisation structure of the mass matrix to the
      ! correct one; at the moment it points to the discretisation structure
      ! of the Stokes matrix...
      p_rmatrixMass%p_rspatialDiscretisation => &
        rproblem%RlevelInfo(i)%p_rdiscretisationMass

      ! Call the standard matrix setup routine to build the matrix.                    
      CALL stdop_assembleSimpleMatrix (p_rmatrixMass,DER_FUNC,DER_FUNC)
                  

    END DO

    ! Clean up the collection (as we are done with the assembly, that's it.
    CALL c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_generateBasicSystemMatrix (rproblem,ilev,rmatrix)
  
!<description>
  ! Generates the basic system matrix of level ilev into the matrix rmatrix.
  ! For this purpose, rmatrix is connected to the basic underlying
  ! system matrix in the problem structure, and the entries of the
  ! static submatrices are copied to rmatrix.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! Number of the level that should be initialised.
  INTEGER, INTENT(IN) :: ilev
  
  ! OPTIONAL: A block matrix structure that receives the basic system matrix.
  ! If not present, the routine will initialise only the basic system
  ! matrix in the problem structure.
  TYPE(t_matrixBlock), INTENT(INOUT), OPTIONAL :: rmatrix
!</inputoutput>

!</subroutine>

    ! A pointer to the system matrix and the RHS/solution vectors.
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    
    REAL(DP) :: dalphaC
    
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dalphaC',dalphaC,0.1_DP)
  
    ! -----------------------------------------------------------------------
    ! Basic (Navier-) Stokes problem
    ! -----------------------------------------------------------------------

    ! Get the main system matrix, which is a block matrix.
    p_rmatrix => rproblem%RlevelInfo(ilev)%rmatrix
    
    ! Let's consider the global system in detail:
    !
    !    ( A11       B1  M/a          ) = ( A11  A12  A13 A14           )
    !    (      A22  B2       M/a     )   ( A21  A22  A23      A25      )
    !    ( B1^T B2^T                  )   ( A31  A32  A33               )
    !    ( -M            A44  A45  B1 )   ( A41           A44  A45  A46 )
    !    (      -M       A54  A55  B2 )   (      A52      A54  A55  A56 )
    !    (               B1^T B2^T    )   (               A64  A65  A66 )

    ! The purpose of this routine is to initialise the static parts of this
    ! matrix. The two A-blocks are nonlinear, ilev.e. dynamic, so we don't
    ! have to deal with them. What we do here is to initialise the
    ! B-matrices and the mass matrices of this system according to given
    ! templates!
    !
    ! The B1/B2 matrices exist up to now only in our local problem structure
    ! as scalar 'template' matrices. Put a copy of them into the block matrix.
    !
    ! Note that we share the structure of B1/B2 with those B1/B2 of the
    ! block matrix, while we create copies of the entries. The B-blocks
    ! are already prepared and memory for the entries is already allocated;
    ! so we only have to copy the entries.
    CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(ilev)%rmatrixB1, &
                                  p_rmatrix%RmatrixBlock(1,3),&
                                  LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)

    CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(ilev)%rmatrixB2, &
                                  p_rmatrix%RmatrixBlock(2,3),&
                                  LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
    
    ! Furthermore, put B1^T and B2^T to the block matrix.
    CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(ilev)%rmatrixB1, &
                                  p_rmatrix%RmatrixBlock(3,1),&
                                  LSYSSC_TR_VIRTUAL)

    CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(ilev)%rmatrixB2, &
                                  p_rmatrix%RmatrixBlock(3,2),&
                                  LSYSSC_TR_VIRTUAL)

    ! Include B1, B1^T, B2, B2^T also to the appropriate blocks for the
    ! dual system

    CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(ilev)%rmatrixB1, &
                                  p_rmatrix%RmatrixBlock(4,6),&
                                  LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)

    CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(ilev)%rmatrixB2, &
                                  p_rmatrix%RmatrixBlock(5,6),&
                                  LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
    
    ! Furthermore, put B1^T and B2^T to the block matrix.
    CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(ilev)%rmatrixB1, &
                                  p_rmatrix%RmatrixBlock(6,4),&
                                  LSYSSC_TR_VIRTUAL)

    CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(ilev)%rmatrixB2, &
                                  p_rmatrix%RmatrixBlock(6,5),&
                                  LSYSSC_TR_VIRTUAL)
                                
    ! Include the mass matrices into the global system. Share the structure,
    ! duplicate the content. It's important not to share the content as
    ! boundary implementation routines might overwrite the content!
    !
    ! The mass matrices on the upper right of the system are weighted 
    ! by 1/dalphaC.
    CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(ilev)%rmatrixMass, &
                                  p_rmatrix%RmatrixBlock(1,4),&
                                  LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
    p_rmatrix%RmatrixBlock(1,4)%dscaleFactor = 1.0_DP/dalphaC
                                  
    CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(ilev)%rmatrixMass, &
                                  p_rmatrix%RmatrixBlock(2,5),&
                                  LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
    p_rmatrix%RmatrixBlock(2,5)%dscaleFactor = 1.0_DP/dalphaC
    

    CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(ilev)%rmatrixMass, &
                                  p_rmatrix%RmatrixBlock(4,1),&
                                  LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
    p_rmatrix%RmatrixBlock(4,1)%dscaleFactor = -1.0_DP
                                  
    CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(ilev)%rmatrixMass, &
                                  p_rmatrix%RmatrixBlock(5,2),&
                                  LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
    p_rmatrix%RmatrixBlock(5,2)%dscaleFactor = -1.0_DP

    ! Make a copy of the system matrix to rmatrix --
    ! or to be more concrete: link rmatrix to the system matrix without
    ! copying and memory.
    IF (PRESENT(rmatrix)) THEN
      CALL lsysbl_duplicateMatrix (p_rmatrix,rmatrix,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_generateStaticSystemParts (rproblem)
  
!<description>
  ! Generates the static parts of the main system matrix of the coupled
  ! global system matrix. More precisely, this copies the entries of
  ! previously calculated B- and mass matrices into the global system matrix, 
  ! so that they can be modified for boundary conditions.
  !
  ! For the optimal control problem, the mass matrices will be copied into the
  ! system matrix and be scaled according to the definition of the stationary
  ! problem.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i
  
    ! For each level, call the generation routine.
    DO i=rproblem%NLMIN,rproblem%NLMAX
      CALL c2d2_generateBasicSystemMatrix (rproblem,i)
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_generateBasicRHS (rproblem,rrhs)
  
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
    CALL c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)

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
              
    ! Switch the sign of the target velocity field because the RHS of the
    ! dual equation is '-z'!
    CALL lsyssc_scaleVector (rrhs%RvectorBlock(4),-1.0_DP)
    CALL lsyssc_scaleVector (rrhs%RvectorBlock(5),-1.0_DP)

    ! Dual pressure RHS is =0.
    CALL lsyssc_clearVector(rrhs%RvectorBlock(6))
                                
    ! Clean up the collection (as we are done with the assembly, that's it.
    CALL c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_doneMatVec (rproblem,rvector,rrhs)
  
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
      ! Delete the system matrix.
      CALL lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)

      ! If there is an existing mass matrix, release it.
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixMass)

      ! Release Stokes, B1 and B2 matrix
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB2)
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB1)
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixStokes)
      
      ! Release the template matrices. This is the point, where the
      ! memory of the matrix structure is released.
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixTemplateGradient)
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixTemplateFEM)
      
      ! Remove the temp vector that was used for interpolating the solution
      ! from higher to lower levels in the nonlinear iteration.
      IF (i .LT. rproblem%NLMAX) THEN
        CALL lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector)
        CALL collct_deletevalue(rproblem%rcollection,PAR_TEMPVEC,i)
      END IF
      
    END DO

    ! Delete solution/RHS vector
    CALL lsysbl_releaseVector (rvector)
    CALL lsysbl_releaseVector (rrhs)

  END SUBROUTINE

END MODULE
