!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2discretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic spatial discretisation related routines for 
!# CC2D. Here, matrix and RHS creation routines can be found as well as
!# routines to initialise/clean up discretisation structures and routines
!# to read/write solution vectors.
!#
!# The following routines can be found here:
!#
!# 1.) c2d2_initDiscretisation
!#     -> Initialise the discretisation structure inside of the problem
!#        structure using the parameters from the INI/DAT files.
!#
!# 2.) c2d2_allocMatVec
!#     -> Allocates memory for vectors/matrices on all levels.
!#
!# 4.) c2d2_generateStaticMatrices
!#     -> Assembles matrix entries of static matrices (Stokes, B)
!#        on one level
!#
!# 5.) c2d2_generateBasicMatrices
!#     -> Assembles the matrix entries of all static matrices on all levels.
!#
!# 6.) c2d2_generateBasicRHS
!#     -> Generates a general RHS vector without any boundary conditions
!#        implemented
!#
!# 7.) c2d2_doneMatVec
!#     -> Cleanup of matrices/vectors, release all memory
!#
!# 8.) c2d2_doneDiscretisation
!#     -> Cleanup of the underlying discretisation structures
!#
!# 9.) c2d2_initInitialSolution
!#     -> Init solution vector according to parameters in the DAT file
!#
!# 10.) c2d2_writeSolution
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
  
  USE collection
  USE convection
  USE vectorio
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback
  USE cc2dmediumm2nonlinearcoreinit
  
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
  INTEGER :: I,j,k,ielementType,icubA,icubB,icubF, icubM
  CHARACTER(LEN=SYS_NAMELEN) :: sstr
  
  ! Number of equations in our problem. velocity+velocity+pressure = 3
  INTEGER, PARAMETER :: nequations = 3
  
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
                    p_rdiscretisation%RspatialDiscr(1), &
                    EL_E031,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscr(1)% &
          RelementDistr(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards.
        CALL spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscr(1),&
            p_rdiscretisation%RspatialDiscr(2))
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscr(3), &
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
                    p_rdiscretisation%RspatialDiscr(1), &
                    EL_E030,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscr(1)% &
          RelementDistr(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards.
        CALL spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscr(1),&
            p_rdiscretisation%RspatialDiscr(2))
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscr(3), &
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
                    p_rdiscretisation%RspatialDiscr(1), &
                    EL_EM31,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscr(1)% &
          RelementDistr(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards.
        CALL spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscr(1),&
            p_rdiscretisation%RspatialDiscr(2))
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscr(3), &
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
                    p_rdiscretisation%RspatialDiscr(1), &
                    EL_EM30,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscr(1)% &
          RelementDistr(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards.
        CALL spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscr(1),&
            p_rdiscretisation%RspatialDiscr(2))
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscr(3), &
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
                    p_rdiscretisation%RspatialDiscr(1), &
                    EL_Q2,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscr(1)% &
          RelementDistr(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards.
        CALL spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscr(1),&
            p_rdiscretisation%RspatialDiscr(2))
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscr(3), &
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
                    p_rdiscretisation%RspatialDiscr(1), &
                    EL_EM30_UNPIVOTED,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscr(1)% &
          RelementDistr(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards.
        CALL spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscr(1),&
            p_rdiscretisation%RspatialDiscr(2))
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscr(3), &
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
                    p_rdiscretisation%RspatialDiscr(1), &
                    EL_EM30_UNSCALED,icubA, &
                    p_rtriangulation, p_rboundary)
                    
        ! Manually set the cubature formula for the RHS as the above routine
        ! uses the same for matrix and vectors.
        p_rdiscretisation%RspatialDiscr(1)% &
          RelementDistr(1)%ccubTypeLinForm = icubF
                    
        ! ...and copy this structure also to the discretisation structure
        ! of the 2nd component (Y-velocity). This needs no additional memory, 
        ! as both structures will share the same dynamic information afterwards.
        CALL spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscr(1),&
            p_rdiscretisation%RspatialDiscr(2))
    
        ! For the pressure (3rd component), we set up a separate discretisation 
        ! structure, as this uses different finite elements for trial and test
        ! functions.
        CALL spdiscr_initDiscr_combined ( &
                    p_rdiscretisation%RspatialDiscr(3), &
                    EL_Q0,EL_EM30_UNSCALED,icubB, &
                    p_rtriangulation, p_rboundary)

      CASE DEFAULT
        PRINT *,'Unknown discretisation: iElementType = ',ielementType
        STOP
      END SELECT

      ! -----------------------------------------------------------------------
      ! Time-dependent problem
      ! -----------------------------------------------------------------------

      CALL parlst_getvalue_int_direct (rproblem%rparamList, 'TIME-DISCRETISATION', &
                                       'ITIMEDEPENDENCE', j, 0)
      IF (j .NE. 0) THEN
      
        ! Initialise a discretisation structure for the mass matrix.
        ! Copy the discretisation structure of the first (Stokes) block
        ! and replace the cubature-formula identifier by that which is to be
        ! used for the mass matrix.
        ALLOCATE(p_rdiscretisationMass)
        p_rdiscretisationMass = p_rdiscretisation%RspatialDiscr(1)
        
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
          p_rdiscretisationMass%RelementDistr(k)%ccubTypeBilForm = icubM
        END DO

        ! Should we do mass lumping?
        CALL parlst_getvalue_int_direct (rproblem%rparamList, 'CC-DISCRETISATION', &
                                        'IMASS', j, 0)
                                        
        IF (j .EQ. 0) THEN
        
          ! How to do lumping?
          CALL parlst_getvalue_int_direct (rproblem%rparamList, 'CC-DISCRETISATION', &
                                          'IMASSLUMPTYPE', j, 0)
                                          
          ! Set cubature formula for lumping. The constant from the DAT file corresponds
          ! to one of the LSYSSC_LUMP_xxxx constants for lsyssc_lumpMatrixScalar.
          ! When to do simple mass lumping, replace the cubature formula by one
          ! that is compatible with the corresponding element to generate
          ! a diagonal mass matrix.
          IF (j .EQ. LSYSSC_LUMP_STD) THEN
            
            DO k = 1,p_rdiscretisationMass%inumFESpaces
              
              j = spdiscr_getLumpCubature (&
                  p_rdiscretisationMass%RelementDistr(k)%itrialElement)
              IF (j .NE. 0) THEN
                icubM = j
              ELSE
                PRINT *,'c2d2_initDiscretisation: Unknown cubature formula for &
                        &mass lumping!'
                STOP
              END IF
              
              ! Set the cubature formula appropriately
              p_rdiscretisationMass%RelementDistr(k)%ccubTypeBilForm = icubM

            END DO
          
          END IF
        
        END IF
       
        rproblem%RlevelInfo(i)%p_rdiscretisationMass => p_rdiscretisationMass
        
      END IF
      
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
  INTEGER :: i,j

    DO i=rproblem%NLMAX,rproblem%NLMIN,-1
      
      ! Remove the block discretisation structure and all substructures.
      CALL spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%p_rdiscretisation)
      
      ! Remove the discretisation from the heap.
      DEALLOCATE(rproblem%RlevelInfo(i)%p_rdiscretisation)

      ! -----------------------------------------------------------------------
      ! Time-dependent problem
      ! -----------------------------------------------------------------------

      CALL parlst_getvalue_int_direct (rproblem%rparamList, 'TIME-DISCRETISATION', &
                                       'ITIMEDEPENDENCE', j, 0)
      IF (j .NE. 0) THEN
        ! Release the mass matrix discretisation.
        ! Don't release the content as we created it as a copy of the Stokes
        ! discretisation structure.
        DEALLOCATE(rproblem%RlevelInfo(i)%p_rdiscretisationMass)
      END IF

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
      !    ( A         B1 )
      !    (      A    B2 )
      !    ( B1^T B2^T    )
      !
      ! with A = L + nonlinear Convection. We compute in advance
      ! a standard Stokes matrix L which can be added later to the
      ! convection matrix, resulting in the nonlinear system matrix.
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
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixTemplateFEM,cmatBuildType)

      ! In the global system, there are two gradient matrices B1 and B2.
      ! Create a template matrix that defines their structure.
      p_rmatrixTemplateGradient => rproblem%RlevelInfo(i)%rmatrixTemplateGradient
      
      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
      CALL bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(3),LSYSSC_MATRIX9,&
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
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
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
      
      ! Allocate memory for that matrix with the appropriate construction routine.
      ! The modules "nonlinearcoreinit" and "nonlinearcore" are actually the
      ! only modules that 'know' the structure of the system matrix!
      CALL c2d2_allocSystemMatrix (rproblem,rproblem%RlevelInfo(i),&
          rproblem%RlevelInfo(i)%rpreallocatedSystemMatrix)

      ! -----------------------------------------------------------------------
      ! Temporary vectors
      !
      ! Now on all levels except for the maximum one, create a temporary 
      ! vector on that level, based on the matrix template.
      ! It's used for building the matrices on lower levels.
      IF (i .LT. rproblem%NLMAX) THEN
        p_rtempVector => rproblem%RlevelInfo(i)%rtempVector
        CALL lsysbl_createVecBlockIndMat (&
            rproblem%RlevelInfo(i)%rpreallocatedSystemMatrix,&
            p_rtempVector,.FALSE.)
      END IF

    END DO
    
    ! (Only) on the finest level, we need to have to allocate a RHS vector
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

  SUBROUTINE c2d2_generateStaticMatrices (rproblem,rlevelInfo)
  
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

    ! local variables
    INTEGER :: j

    ! A pointer to the Stokes and mass matrix
    TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes,p_rmatrixMass

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    CALL c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)
    
    ! -----------------------------------------------------------------------
    ! Basic (Navier-) Stokes problem
    ! -----------------------------------------------------------------------
    
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%p_rdiscretisation
    
    ! Get a pointer to the (scalar) Stokes matrix:
    p_rmatrixStokes => rlevelInfo%rmatrixStokes
    
    ! The global system looks as follows:
    !
    !    ( A         B1 )
    !    (      A    B2 )
    !    ( B1^T B2^T    )
    !
    ! with A = L + nonlinear Convection. We compute in advance
    ! a standard Stokes matrix L which can be added later to the
    ! convection matrix, resulting in the nonlinear system matrix,
    ! as well as both B-matrices.
    
!    ! For assembling of the entries, we need a bilinear form, 
!    ! which first has to be set up manually.
!    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
!    ! scalar system matrix in 2D.
!    
!    rform%itermCount = 2
!    rform%Idescriptors(1,1) = DER_DERIV_X
!    rform%Idescriptors(2,1) = DER_DERIV_X
!    rform%Idescriptors(1,2) = DER_DERIV_Y
!    rform%Idescriptors(2,2) = DER_DERIV_Y
!
!    ! In the standard case, we have constant coefficients:
!    rform%ballCoeffConstant = .TRUE.
!    rform%BconstantCoeff = .TRUE.
!    rform%Dcoefficients(1)  = rproblem%dnu
!    rform%Dcoefficients(2)  = rproblem%dnu
!
!    ! Now we can build the matrix entries.
!    ! We specify the callback function coeff_Stokes for the coefficients.
!    ! As long as we use constant coefficients, this routine is not used.
!    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
!    ! the framework will call the callback routine to get analytical data.
!    !
!    ! We pass our collection structure as well to this routine, 
!    ! so the callback routine has access to everything what is
!    ! in the collection.
!    CALL bilf_buildMatrixScalar (rform,.TRUE.,&
!                                 p_rmatrixStokes,coeff_Stokes,&
!                                 rproblem%rcollection)
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
    ! Time-dependent problem
    ! -----------------------------------------------------------------------

    CALL parlst_getvalue_int_direct (rproblem%rparamList, 'TIME-DISCRETISATION', &
                                     'ITIMEDEPENDENCE', j, 0)
    IF (j .NE. 0) THEN
    
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
      p_rmatrixMass%p_rspatialDiscretisation => rlevelInfo%p_rdiscretisationMass

      ! Call the standard matrix setup routine to build the matrix.                    
      CALL stdop_assembleSimpleMatrix (p_rmatrixMass,DER_FUNC,DER_FUNC)
                  
      ! Should we do mass lumping?
      CALL parlst_getvalue_int_direct (rproblem%rparamList, 'CC-DISCRETISATION', &
                                      'IMASS', j, 0)
                                      
      IF (j .EQ. 0) THEN
      
        ! How to do lumping?
        CALL parlst_getvalue_int_direct (rproblem%rparamList, 'CC-DISCRETISATION', &
                                        'IMASSLUMPTYPE', j, 0)
                                        
        ! Lump the mass matrix. The constant from the DAT file corresponds
        ! to one of the LSYSSC_LUMP_xxxx constants for lsyssc_lumpMatrixScalar.
        CALL lsyssc_lumpMatrixScalar (p_rmatrixMass,j)
      
      END IF
      
    END IF

    ! Clean up the collection (as we are done with the assembly, that's it.
    CALL c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_generateBasicMatrices (rproblem)
  
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
      CALL c2d2_generateStaticMatrices (&
          rproblem,rproblem%RlevelInfo(i))
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
              p_rdiscretisation%RspatialDiscr(1),rlinform,.TRUE.,&
              rrhs%RvectorBlock(1),coeff_RHS_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(2),rlinform,.TRUE.,&
              rrhs%RvectorBlock(2),coeff_RHS_y,&
              rproblem%rcollection)
                                
    ! The third subvector must be zero initially - as it represents the RHS of
    ! the equation "div(u) = 0".
    CALL lsyssc_clearVector(rrhs%RvectorBlock(3))
                                
    ! Clean up the collection (as we are done with the assembly, that's it.
    CALL c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)

  END SUBROUTINE

  ! ***************************************************************************

!!subroutine>
!
!  SUBROUTINE c2d2_initMatVec (rproblem)
!  
!!description>
!  ! Calculates entries of all static matrices (Stokes, B-matrices,...)
!  ! of the problem. 
!  ! Calculates the system matrix and RHS vector of the linear system
!  ! by discretising the problem with the default discretisation structure
!  ! in the problem structure.
!  ! Sets up a solution vector for the linear system.
!!/description>
!
!!inputoutput>
!  ! A problem structure saving problem-dependent information.
!  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!!/inputoutput>
!
!  ! local variables
!  INTEGER :: i
!  
!    ! A bilinear and linear form describing the analytic problem to solve
!    TYPE(t_bilinearForm) :: rform
!    TYPE(t_linearForm) :: rlinform
!    
!    ! A pointer to the system matrix and the RHS/solution vectors.
!    TYPE(t_matrixBlock), POINTER :: p_rmatrix
!    TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes
!    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector,p_rtempVector
!
!    ! A pointer to the discretisation structure with the data.
!    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
!  
!    DO i=rproblem%NLMIN,rproblem%NLMAX
!      ! Ask the problem structure to give us the discretisation structure
!      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
!      
!      ! The global system looks as follows:
!      !
!      !    ( A         B1 )
!      !    (      A    B2 )
!      !    ( B1^T B2^T    )
!      !
!      ! with A = L + nonlinear Convection. We compute in advance
!      ! a standard Stokes matrix L which can be added later to the
!      ! convection matrix, resulting in the nonlinear system matrix.
!      !
!      ! Get a pointer to the (scalar) Stokes matrix:
!      p_rmatrixStokes => rproblem%RlevelInfo(i)%rmatrixStokes
!      
!      ! Create the matrix structure of the Stokes matrix:
!      CALL bilf_createMatrixStructure (&
!                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
!                p_rmatrixStokes)
!      
!      ! And now to the entries of the matrix. For assembling of the entries,
!      ! we need a bilinear form, which first has to be set up manually.
!      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
!      ! scalar system matrix in 2D.
!      
!      rform%itermCount = 2
!      rform%Idescriptors(1,1) = DER_DERIV_X
!      rform%Idescriptors(2,1) = DER_DERIV_X
!      rform%Idescriptors(1,2) = DER_DERIV_Y
!      rform%Idescriptors(2,2) = DER_DERIV_Y
!
!      ! In the standard case, we have constant coefficients:
!      rform%ballCoeffConstant = .TRUE.
!      rform%BconstantCoeff = .TRUE.
!      rform%Dcoefficients(1)  = rproblem%dnu
!      rform%Dcoefficients(2)  = rproblem%dnu
!
!      ! Now we can build the matrix entries.
!      ! We specify the callback function coeff_Stokes for the coefficients.
!      ! As long as we use constant coefficients, this routine is not used.
!      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
!      ! the framework will call the callback routine to get analytical data.
!      !
!      ! We pass our collection structure as well to this routine, 
!      ! so the callback routine has access to everything what is
!      ! in the collection.
!      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
!                                   p_rmatrixStokes,coeff_Stokes,&
!                                   rproblem%rcollection)
!      
!      ! In the global system, there are two coupling matrices B1 and B2.
!      ! Both have the same structure.
!      ! Create the matrices structure of the pressure using the 3rd
!      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
!      CALL bilf_createMatrixStructure (&
!                p_rdiscretisation%RspatialDiscr(3),LSYSSC_MATRIX9,&
!                rproblem%RlevelInfo(i)%rmatrixB1)
!                
!      ! Duplicate the B1 matrix structure to the B2 matrix, so use
!      ! lsyssc_duplicateMatrix to create B2. Share the matrix 
!      ! structure between B1 and B2 (B1 is the parent and B2 the child). 
!      ! Don't create a content array yet, it will be created by 
!      ! the assembly routines later.
!      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1,&
!                  rproblem%RlevelInfo(i)%rmatrixB2,LSYSSC_DUP_COPY,LSYSSC_DUP_REMOVE)
!
!      ! Build the first pressure matrix B1.
!      ! Again first set up the bilinear form, then call the matrix assembly.
!      rform%itermCount = 1
!      rform%Idescriptors(1,1) = DER_FUNC
!      rform%Idescriptors(2,1) = DER_DERIV_X
!
!      ! In the standard case, we have constant coefficients:
!      rform%ballCoeffConstant = .TRUE.
!      rform%BconstantCoeff = .TRUE.
!      rform%Dcoefficients(1)  = -1.0_DP
!      
!      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
!                                  rproblem%RlevelInfo(i)%rmatrixB1,coeff_Pressure,&
!                                  rproblem%rcollection)
!
!      ! Build the second pressure matrix B2.
!      ! Again first set up the bilinear form, then call the matrix assembly.
!      rform%itermCount = 1
!      rform%Idescriptors(1,1) = DER_FUNC
!      rform%Idescriptors(2,1) = DER_DERIV_Y
!
!      ! In the standard case, we have constant coefficients:
!      rform%ballCoeffConstant = .TRUE.
!      rform%BconstantCoeff = .TRUE.
!      rform%Dcoefficients(1)  = -1.0_DP
!      
!      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
!                                  rproblem%RlevelInfo(i)%rmatrixB2,coeff_Pressure,&
!                                  rproblem%rcollection)
!                                  
!      ! Now let's come to the main system matrix, which is a block matrix.
!      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
!      
!      ! Initialise the block matrix with default values based on
!      ! the discretisation.
!      CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)    
!      
!      ! Inform the matrix that we build a saddle-point problem.
!      ! Normally, imatrixSpec has the value LSYSBS_MSPEC_GENERAL,
!      ! but probably some solvers can use the special structure later.
!      p_rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
!      
!      ! Let's consider the global system in detail:
!      !
!      !    ( A         B1 ) = ( A11  A12  A13 )
!      !    (      A    B2 )   ( A21  A22  A23 )
!      !    ( B1^T B2^T    )   ( A31  A32  A33 )
!      !
!      ! The matrices A11 and A22 of the global system matrix have exactly
!      ! the same structure as the original Stokes matrix from above!
!      ! Initialise them with the same structure, i.e. A11, A22 and the
!      ! Stokes matrix L share(!) the same structure.
!      !
!      ! For this purpose, use the "duplicate matric" routine.
!      ! The structure of the matrix is shared with the Stokes matrix.
!      ! For the content, a new empty array is allocated which will later receive
!      ! the entries.
!      CALL lsyssc_duplicateMatrix (p_rmatrixStokes,&
!                  p_rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!        
!      IF (.NOT. rproblem%bdecoupledXY) THEN          
!        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix 
!        ! A22 is identical to A11! So mirror A11 to A22 sharing the
!        ! structure and the content.
!        CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
!                    p_rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!        ! Save the value of bdecoupledXY to the collection.
!        CALL collct_setvalue_int(rproblem%rcollection,'DECOUPLEDXY',NO,.TRUE.)
!      ELSE
!        ! Otherwise, create another copy of the Stokes matrix.
!        CALL lsyssc_duplicateMatrix (p_rmatrixStokes,&
!                    p_rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!        ! Save the value of bdecoupledXY to the collection.
!        CALL collct_setvalue_int(rproblem%rcollection,'DECOUPLEDXY',YES,.TRUE.)
!      END IF
!
!      ! Manually change the discretisation structure of the Y-velocity 
!      ! matrix to the Y-discretisation structure.
!      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
!      ! so this is not really necessary - we do this for sure...
!      p_rmatrix%RmatrixBlock(2,2)%p_rspatialDiscretisation => &
!        p_rdiscretisation%RspatialDiscr(2)
!                                  
!      ! The B1/B2 matrices exist up to now only in our local problem structure.
!      ! Put a copy of them into the block matrix.
!      !
!      ! Note that we share the structure of B1/B2 with those B1/B2 of the
!      ! block matrix, while we create copies of the entries. The reason is
!      ! that these matrices are modified for bondary conditions later.
!      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
!                                   p_rmatrix%RmatrixBlock(1,3),&
!                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!
!      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
!                                   p_rmatrix%RmatrixBlock(2,3),&
!                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!      
!      ! Furthermore, put B1^T and B2^T to the block matrix.
!      CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
!                                   p_rmatrix%RmatrixBlock(3,1),&
!                                   LSYSSC_TR_VIRTUAL)
!
!      CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
!                                   p_rmatrix%RmatrixBlock(3,2),&
!                                   LSYSSC_TR_VIRTUAL)
!
!      ! Update the structural information of the block matrix, as we manually
!      ! changed the submatrices:
!      CALL lsysbl_updateMatStrucInfo (p_rmatrix)
!      
!      ! Now on all levels except for the maximum one, create a temporary 
!      ! vector on that level, based on the matrix template.
!      ! It's used for building the matrices on lower levels.
!      IF (i .LT. rproblem%NLMAX) THEN
!        p_rtempVector => rproblem%RlevelInfo(i)%rtempVector
!        CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rtempVector,.FALSE.)
!        
!        ! Add the temp vector to the collection on level i
!        ! for use in the callback routine
!        CALL collct_setvalue_vec(rproblem%rcollection,PAR_RTEMPVEC,p_rtempVector,&
!                                .TRUE.,i)
!      END IF
!
!    END DO
!
!    ! (Only) on the finest level, we need to calculate a RHS vector
!    ! and to allocate a solution vector.
!    
!    p_rrhs    => rproblem%rrhs   
!    p_rvector => rproblem%rvector
!
!    ! Although we could manually create the solution/RHS vector,
!    ! the easiest way to set up the vector structure is
!    ! to create it by using our matrix as template:
!    CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rrhs, .FALSE.)
!    CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rvector, .FALSE.)
!
!    ! Save the solution/RHS vector to the collection. Might be used
!    ! later (e.g. in nonlinear problems)
!    CALL collct_setvalue_vec(rproblem%rcollection,PAR_RHS,p_rrhs,.TRUE.)
!    CALL collct_setvalue_vec(rproblem%rcollection,PAR_SOLUTION,p_rvector,.TRUE.)
!    
!    ! The vector structure is ready but the entries are missing. 
!    ! So the next thing is to calculate the content of that vector.
!    !
!    ! At first set up the corresponding linear form (f,Phi_j):
!    rlinform%itermCount = 1
!    rlinform%Idescriptors(1) = DER_FUNC
!    
!    ! ... and then discretise the RHS to the first subvector of
!    ! the block vector using the discretisation structure of the 
!    ! first block.
!    !
!    ! We pass our collection structure as well to this routine, 
!    ! so the callback routine has access to everything what is
!    ! in the collection.
!    !
!    ! Note that the vector is unsorted after calling this routine!
!    CALL linf_buildVectorScalar (&
!              p_rdiscretisation%RspatialDiscr(1),rlinform,.TRUE.,&
!              p_rrhs%RvectorBlock(1),coeff_RHS_x,&
!              rproblem%rcollection)
!                                
!    ! The third subvector must be zero - as it represents the RHS of
!    ! the equation "div(u) = 0".
!    CALL lsyssc_clearVector(p_rrhs%RvectorBlock(3))
!                                
!    ! Clear the solution vector on the finest level.
!    CALL lsysbl_clearVector(rproblem%rvector)
!    
!  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initInitialSolution (rproblem,rvector)
  
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

  SUBROUTINE c2d2_writeSolution (rproblem,rvector)
  
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
      CALL lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rpreallocatedSystemMatrix)

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
      END IF
      
    END DO

    ! Delete solution/RHS vector
    CALL lsysbl_releaseVector (rvector)
    CALL lsysbl_releaseVector (rrhs)

  END SUBROUTINE

END MODULE
