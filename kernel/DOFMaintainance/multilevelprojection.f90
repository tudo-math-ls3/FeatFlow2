!##############################################################################
!# ****************************************************************************
!# <name> multilevelprojection </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains routines for prolongation, restriction and interpolation
!# of solution vectors between different levels, for scalar systems as well
!# as for block systems.
!#
!# The module contains the following routines:
!#
!# 1.) mlprj_initProjection
!#     -> Initialises a projection structure with values according to a
!#        spatial discretisation
!#
!# 2.) mlprj_doneProjection
!#     -> Cleans up a projection structure.
!#
!# 3.) mlprj_getTempMemory
!#     -> Determine the amount of temporary storage that is necessary
!#        to transport a vector from one level to another
!#
!# 4.) mlprj_performProlongation
!#     -> Interpolates a solution vector from a coarse grid to a fine grid
!#        (L2-projection in the primal space)
!#
!# 5.) mlprj_performRestriction
!#     -> Restricts a defect vector from a fine grid to a coarse grid
!#        (L2-projection in the dual space)
!#
!# 6.) mlprj_performInterpolation
!#     -> Interpolated a solution from a fine grid to a coarse grid
!#        (L2-projection in the primal space)
!# 
!# </purpose>
!##############################################################################

MODULE multilevelprojection

  USE fsystem
  USE spatialdiscretisation
  USE linearsystemscalar
  USE linearsystemblock
  USE triangulation
  USE element
  
  IMPLICIT NONE

!<types>
  
!<typeblock>
  
  ! Describes for a scalar equation which type of prolongation, restriction
  ! or interpolation to be used. The information is (of course) element
  ! dependent.
  ! One t_interlevelProjectionScalar structure is associated to one
  ! t_elementDistribution structure of a spatial discretisation.
  ! t_interlevelProjectionScalar is not part if the t_spatialDiscretisation
  ! structure, as t_interlevelProjectionScalar resides 'between' levels:
  ! If there is one t_spatialDiscretisation structure for level 3 and
  ! one for level 4, the t_interlevelProjectionScalar structure is
  ! logically arranged at level '3.5' to configure the projection
  ! between level 3 and 4. For this reason, there is one 
  ! t_interlevelProjectionScalar less than levels: As there is no level 0,
  ! the structure for level '0.5' is missing!
  
  TYPE t_interlevelProjectionScalar
  
    ! Element type that should be assumed for the prolongation. Must fit
    ! to the element type of the discretisation concerning the DOF's:
    ! If element "EM30" is used e.g., prolongation for "EM30", "EM31",
    ! "E031" and "E030" is possible, although performance will suffer
    ! in using the "wrong" prolongation. An error will be displayed if
    ! the DOF's don't fit together at all, e.g. when using $Q_1$
    ! prolongation when using $Q_2$ for the discretisation.
    !
    ! A value of EL_UNDEFINED indicates that the 'natural' prolongation
    ! (i.e. that configured by the spatial discretisation) should be used.
    INTEGER                     :: ielementTypeProlongation = EL_UNDEFINED
    
    ! Order of the prolongation to use. 
    ! -1=use default prolongation (e.g. linear for $Q_1$, quadratic for 
    !    $Q_2$,...). Some elements allow to configure the type of the 
    !    prolongation, e.g. when $Q_2$ is used, apart from 0 the following 
    !    values are allowed:
    ! 0=constant prolongation
    ! 1=linear prolongation of a once refined mesh
    ! 2=quadratic interpolation
    INTEGER                     :: iprolongationOrder = -1

    ! Element type that should be assumed for the restriction. Must fit
    ! to the element type of the discretisation concerning the DOF's:
    ! If element "EM30" is used e.g., restriction for "EM30", "EM31",
    ! "E031" and "E030" is possible, although performance will suffer
    ! in using the "wrong" restriction. An error will be displayed if
    ! the DOF's don't fit together at all, e.g. when using $Q_1$
    ! restriction when using $Q_2$ for the discretisation.
    !
    ! A value of EL_UNDEFINED indicates that the 'natural' restriction
    ! (i.e. that configured by the spatial discretisation) should be used.
    INTEGER                     :: ielementTypeRestriction = EL_UNDEFINED

    ! Order of the restriction to use. 
    ! -1=use default restriction (e.g. linear for $Q_1$, quadratic for 
    !    $Q_2$,...). Some elements allow to configure the type of the 
    !    restriction, e.g. when $Q_2$ is used, apart from 0 the following 
    !    values are allowed:
    ! 0=constant restriction
    ! 1=linear restriction of a once refined mesh
    ! 2=quadratic restriction
    INTEGER                     :: irestrictionOrder = -1

    ! Element type that should be assumed for the interpolation of a
    ! solution to a lower level. Must fit to the element type of the 
    ! discretisation concerning the DOF's.
    ! An error will be displayed if the DOF's don't fit together at 
    ! all, e.g. when using $Q_1$ interpolation when using $Q_2$ 
    ! for the discretisation.
    !
    ! A value of EL_UNDEFINED indicates that the 'natural' interpolation
    ! (i.e. that configured by the spatial discretisation) should be used.
    INTEGER                     :: ielementTypeInterpolation = EL_UNDEFINED
    
    ! Order of the interpolation to use when interpolating a solution vector
    ! to a lower level.
    ! -1=use default interpolation (e.g. linear for $Q_1$, quadratic for 
    !    $Q_2$,...). Some elements allow to configure the type of the 
    !    interpolation, e.g. when $Q_2$ is used, apart from 0 the following 
    !    values are allowed:
    ! 0=constant interpolation
    ! 1=linear interpolation of a once refined mesh
    ! 2=quadratic interpolation
    INTEGER                     :: iinterpolationOrder = -1
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! Contains a list of t_interlevelProjectionScalar structures that
  ! describes the projection of block vectors between different levels
  ! in the discretisation.
  ! The structure contains a 2-dimensional array: There is one
  ! t_interlevelProjectionScalar structure for every equation and
  ! for every element distribution in the equation.
  !
  ! The default initialisation initialises this structure with the
  ! usual values for the whole grid transfer process. If there is special
  ! need for 'reconfiguring' the grid transfer for a special element
  ! type in a special equation, the application can change the content.
  
  TYPE t_interlevelProjectionBlock
  
    ! A list of t_interlevelProjectionScalar structures for every
    ! equation and every element distribution in the discretisation.
    TYPE(t_interlevelProjectionScalar), &
      DIMENSION(SPDISC_MAXFESPACES,LSYSBL_MAXBLOCKS) :: rscalarProjection
  
  END TYPE
  
  !</typeblock>

!</types>

CONTAINS

!<subroutine>

  SUBROUTINE mlprj_initProjection (rprojection,RspatialDiscretisation)
  
!<description>
  ! This subroutine initialises an t_interlevelProjectionBlock with default
  ! values for a given PDE system. This allows the interlevel-projection
  ! routines to calculate constants that are necessary during a level
  ! change. (This is used e.g. in the general prolongation/restriction
  ! where the prolongation/restriction matrix must be calculated).
  !
  ! The calculated information is saved to rprojection and released
  ! with mlprj_doneProjection.
!</description>

!<input>
  ! An array of discretisation structures. Each discretisation structure
  ! corresponds to one scalar equation. The projection structure will be
  ! prepared to according to this discretisation list.
  TYPE(t_spatialDiscretisation), DIMENSION(:), INTENT(IN) :: RspatialDiscretisation
!</input>
  
!<output>
  ! A t_interlevelProjectionBlock structure that will be filled with data
  ! about the projection of all the equations described by RspatialDiscretisation.
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
!</output>
  
!</subroutine>

  ! local variables
  ! ...
  ! At the moment, there is nothing here to be precomputed.
  ! In a later implementation, this might change if a special prolongation/
  ! restriction needs to precalculate a prolongation/restriction matrix.
    
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_doneProjection (rprojection)
  
!<description>
  ! Cleans up a t_interlevelProjectionBlock structure. All dynamically allocated
  ! memory is released.
!</description>
  
!<inputoutput>
  ! The t_interlevelProjectionBlock structure which is to be cleaned up.
  TYPE(t_interlevelProjectionBlock), INTENT(INOUT) :: rprojection 
!</inputoutput>
  
!</subroutine>

  ! As in the current implementation there is nothing to initialise,
  ! there's also nothing to clean up!

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_getProjectionStrategy (rprojection,releDistrCoarse,releDistrFine,&
                                          ractProjection)
  
!<description>
  ! Internal subroutine. This creates a projection structure ractProjection
  ! from a template projection structure rprojection and the discretisation
  ! structures on the coarse and fine grid and checks compatibility. 
  !
  ! An error is thrown if either the DOF's on the coarse- and fine grid don't 
  ! fit together (i.e. there is no compatible projection) or if they don't fit 
  ! to the projection template (e.g. if ielementTypeProlongation=EL_Q1 when the
  ! actual discretisation structure uses itrialElement=EL_Q2).
  ! Standard values (e.g. ielementTypeProlongation=EL_UNDEFINED) are replaced
  ! by the actual values (e.g. ielementTypeProlongation=EL_Q2,...).
!</description>
  
!<input>
  ! The t_interlevelProjectionScalar structure which is used as a template
  TYPE(t_interlevelProjectionScalar), INTENT(IN) :: rprojection 
  
  ! One element distribution structure on the coarse grid
  TYPE(t_elementDistribution), INTENT(IN) :: releDistrCoarse

  ! One element distribution structure on the fine grid
  TYPE(t_elementDistribution), INTENT(IN) :: releDistrFine
!</input>

!<output>
  ! The t_interlevelProjectionScalar structure which configures the actual
  ! grid transfer between rdiscrCoarse and rdiscrFine.
  TYPE(t_interlevelProjectionScalar), INTENT(OUT) :: ractProjection 
!</output>
  
!</subroutine>

  ! Check to see if the two discretisation structures fit together
  IF (releDistrCoarse%itrialElement .NE. releDistrFine%itrialElement) THEN
    PRINT *,'Element distribution on the coarse and fine grid incompatible!'
    PRINT *,'Coarse grid: ',releDistrCoarse%itrialElement,&
            ' Fine grid: ',releDistrFine%itrialElement
    STOP
  END IF

  ! Copy the template to the actual projection structure.
  ractProjection = rprojection
  
  ! Is any element type to be replaced according to a discretisation?
  IF (ractProjection%ielementTypeProlongation .EQ. EL_UNDEFINED) THEN
    ractProjection%ielementTypeProlongation = releDistrCoarse%itrialElement
  END IF

  IF (ractProjection%ielementTypeRestriction .EQ. EL_UNDEFINED) THEN
    ractProjection%ielementTypeRestriction = releDistrCoarse%itrialElement
  END IF

  IF (ractProjection%ielementTypeInterpolation .EQ. EL_UNDEFINED) THEN
    ractProjection%ielementTypeInterpolation = releDistrCoarse%itrialElement
  END IF

  ! Check to see if the discretisation structures fit to the projection structure
  IF (ractProjection%ielementTypeProlongation .NE. releDistrCoarse%itrialElement) THEN
    ! Here, we can insert some additional code so that E030 is compatible eith EM30
    ! and so on...
    PRINT *,'Element distribution of the grid and interlevel projection incompatible!'
    PRINT *,'Grid: ',releDistrCoarse%itrialElement,&
            ' Prolongation: ',ractProjection%ielementTypeProlongation
    STOP
  END IF

  IF (ractProjection%ielementTypeRestriction .NE. releDistrCoarse%itrialElement) THEN
    ! Here, we can insert some additional code so that E030 is compatible eith EM30
    ! and so on...
    PRINT *,'Element distribution of the grid and interlevel projection incompatible!'
    PRINT *,'Grid: ',releDistrCoarse%itrialElement,&
            ' Restriction: ',ractProjection%ielementTypeRestriction
    STOP
  END IF

  IF (ractProjection%ielementTypeInterpolation .NE. releDistrCoarse%itrialElement) THEN
    ! Here, we can insert some additional code so that E030 is compatible eith EM30
    ! and so on...
    PRINT *,'Element distribution of the grid and interlevel projection incompatible!'
    PRINT *,'Grid: ',releDistrCoarse%itrialElement,&
            ' Interpolation: ',ractProjection%ielementTypeInterpolation
    STOP
  END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<function>

  INTEGER(PREC_VECIDX) FUNCTION mlprj_getTempMemoryScalar (rprojectionScalar,&
                                               rdiscrCoarse,rdiscrFine)
  
!<description>
  ! This function returns for a given projection-structure and given
  ! discretisation structures on the coarse and fine grid the amount of temporary 
  ! memory that is needed to do prolongation, restriction or interpolation.
  ! A value of 0 of course indicates that no temporary memory is needed.
!</description>
  
!<input>
  ! The t_interlevelProjectionScalar structure which configures
  ! the projection between rdiscrCoarse and rdiscrFine
  ! for each of the element distributions in rdiscrCoarse /rdiscrFine
  TYPE(t_interlevelProjectionScalar), DIMENSION(:), INTENT(IN) :: RprojectionScalar
  
  ! The element distribution structure of the equation on the coarse grid.
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscrCoarse

  ! The element distribution structure of the equation on the fine grid.
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscrFine
!</input>

!<result>
  ! Amount of temporary memory needed for prolongation/restriction/interpolation
  ! between vectors corresponding to the given combination of spatial
  ! discretisations.
  ! =0 if no memory is necessary.
!</result>
  
!</function>

  ! Currently, there is no additional memory needed.
  mlprj_getTempMemoryScalar = 0
  
  ! In case, P_0 is interpolated linearly, e.g., there might be some memory
  ! necessary!

  END FUNCTION
  
  ! ***************************************************************************
  
!<function>

  INTEGER(PREC_VECIDX) FUNCTION mlprj_getTempMemory (rprojection, &
                                                     RdiscrCoarse,RdiscrFine)
  
!<description>
  ! Returns the amount of temporary memory that is needed by the
  ! interlevel-routines to transfer vectors between two grids.
  ! RdiscrCoarse and RdiscrFine are are arrays of discretisation structures.
  ! Each discretisation structure corresponds to one scalar equation.
!</description>
  
!<input>
  ! Projection structure that configures the grid transfer for all equations
  ! and all element distributions in each equation.
  TYPE(t_interlevelProjectionBlock), INTENT(IN) :: rprojection 
  
  ! List of disretisation structures fo the equations on the Coarse grid 
  TYPE(t_spatialDiscretisation), DIMENSION(:), INTENT(IN) :: RdiscrCoarse

  ! List of disretisation structures fo the equations on the Fine grid 
  TYPE(t_spatialDiscretisation), DIMENSION(:), INTENT(IN) :: RdiscrFine
!</input>

!<result>
  ! Length of a temporary vector that is necessary for transferring vectors
  ! on the coarse grid to the fine grid and vice versa.
  ! =0 if no memory is necessary.
!</result>

  ! local variables
  INTEGER i
  INTEGER(PREC_VECIDX) :: imemmax,imemact

  IF (SIZE(RdiscrCoarse) .NE. SIZE(RdiscrFine)) THEN
    PRINT *,'mlprj_allocTempVector: Coarse and fine grid incompatible!'
    STOP
  END IF
  
  ! How much memory do we need?
  ! As we perform the grid transfer equation-by-equation and element-
  ! distribution by element-distribution, we need the maximum of all 
  ! mlprj_getTempMemoryScalar calls!
  !
  ! Loop through all blocks:
  
  imemmax = 0
  DO i=1,SIZE(RdiscrCoarse)
    ! How much memory needed for that equation?
    imemact = mlprj_getTempMemoryScalar (rprojection%rscalarProjection(:,i),&
                                         RdiscrCoarse(i), RdiscrFine(i))
            
    imemmax = MAX(imemmax,imemact)
  END DO
  
  ! Now we know how mich we need:
  mlprj_getTempMemory = imemmax
  
  END FUNCTION
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_performProlongation (rprojection,rcoarseVector, &
                                        rfineVector,rtempVector)
  
!<description>
  ! Performs a prolongation for a given block vector (i.e. a projection
  ! in the primal space where the solution lives). The vector
  ! rcoarseVector on a coarser grid is projected to the vector
  ! rfineVector on a finer grid. 
  ! rprojection configures how the grid transfer is performed.
  ! This projection structure must be build corresponding to the spatial 
  ! discretisation structures in rcoarseVector and rfineVector!
!</description>
  
!<input>
  ! The t_interlevelProjectionBlock structure that configures the grid transfer
  TYPE(t_interlevelProjectionBlock), INTENT(IN) :: rprojection 
  
  ! Coarse grid vector
  TYPE(t_vectorBlock), INTENT(IN) :: rcoarseVector
!</input>
  
!<inputoutput>
  ! Temporary vector. The vector must have the same data type as
  ! rfineVector and must be at least as long as indicated by the function
  ! mlprj_getTempMemory. If mlprj_getTempMemory was =0, the vector may
  ! be a dummy vector.
  ! The vector does not have to be connected to a discretisation structure
  ! or something similar; the content is undefined at entry and will be
  ! undefined when leaving this routine.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector
  
  ! Fine grid vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: rfineVector
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrCoarse,p_rdiscrFine
  TYPE(t_triangulation), POINTER :: p_rtriaCoarse,p_rtriaFine
  TYPE(t_interlevelProjectionScalar) :: ractProjection
  
  ! Pointers into the triangulation
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementCoarse
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementFine
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementCoarse
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementFine
  
  ! Data arrays
  REAL(DP), DIMENSION(:), POINTER :: p_DuCoarse, p_DuFine
 
    ! The vectors must be of data type DOUBLE - we don't support anything
    ! different at the moment...
    IF (rcoarseVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Coarse grid vector has unsupported data type!'
      STOP
    END IF

    IF (rfineVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Fine grid vector has unsupported data type!'
      STOP
    END IF
    
    ! Calls the correct prolongation routine for each block in the 
    ! discretisation...
    DO i=1,rcoarseVector%nblocks
    
      IF (rcoarseVector%RvectorBlock(i)%NEQ .GT. 0) THEN
        p_rdiscrCoarse => rcoarseVector%RvectorBlock(i)%p_rspatialDiscretisation
        p_rdiscrFine => rfineVector%RvectorBlock(i)%p_rspatialDiscretisation
        
        ! We need a discretisation:
        IF ((.NOT. ASSOCIATED(p_rdiscrCoarse)) .OR. &
            (.NOT. ASSOCIATED(p_rdiscrFine))) THEN
          PRINT *,'Intergrid transfer: No discretisation!'
        END IF

        ! Currently, we support only uniform triangulations.
        IF ((p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
            (p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM)) THEN
          PRINT *,'Intergrid transfer supports currently only uniform discretisations!'
          STOP
        END IF
        
        ! Get the pointers to the vectors
        CALL storage_getbase_double (rcoarseVector%RvectorBlock(1)%h_Ddata,p_DuCoarse)
        p_DuCoarse => p_DuCoarse( &
          rcoarseVector%RvectorBlock(1)%iidxFirstEntry: &
          rcoarseVector%RvectorBlock(1)%iidxFirstEntry+rcoarseVector%RvectorBlock(1)%NEQ)

        CALL storage_getbase_double (rfineVector%RvectorBlock(1)%h_Ddata,p_DuFine)
        p_DuFine => p_DuFine( &
          rfineVector%RvectorBlock(1)%iidxFirstEntry: &
          rfineVector%RvectorBlock(1)%iidxFirstEntry+rfineVector%RvectorBlock(1)%NEQ)
        
        ! Use the first projection structure as template and create
        ! the actual projection structure for our situation.
        ! Remember, in a uniform grid we only have one projection structure
        ! and one element distribution!
        CALL mlprj_getProjectionStrategy (rprojection%RscalarProjection(1,1), &
              p_rdiscrCoarse%RelementDistribution(1), &
              p_rdiscrFine%RelementDistribution(1), &
              ractProjection)
      
        ! Depending on the element type of the trial functions in the
        ! discretisation, choose the right prolongation and call it.
        p_rtriaCoarse => p_rdiscrCoarse%p_rtriangulation
        p_rtriaFine => p_rdiscrFine%p_rtriangulation
        SELECT CASE (ractProjection%ielementTypeProlongation)
        CASE (EL_P1)
          ! P1 prolongation
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
          CALL mlprj_prolUniformP1_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
               p_IneighboursAtElementCoarse,p_rtriaCoarse%NEL)
               
        CASE (EL_Q0)
          ! Q0 prolongation
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_prolUniformQ0_double (p_DuCoarse,p_DuFine, &
               p_IneighboursAtElementFine, &
               p_rtriaCoarse%NEL)
               
        CASE (EL_Q1)
          ! Q1 prolongation
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
          CALL mlprj_prolUniformQ1_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
               p_IneighboursAtElementCoarse,p_rtriaCoarse%NEL)
        CASE DEFAULT
          PRINT *,'Unsupported prolongation!'
          STOP
        END SELECT
      
      END IF

    END DO  ! i

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_performRestriction (rprojection,rcoarseVector, &
                                       rfineVector,rtempVector)
  
!<description>
  ! Performs a restriction for a given block vector (i.e. a projection
  ! in the dual space where the RHS vector lives). The vector
  ! rfineVector on a finer grid is projected to the vector
  ! rcoarseVector on a coarser grid. 
  ! rprojection configures how the grid transfer is performed.
  ! This projection structure must be build corresponding to the spatial 
  ! discretisation structures in rcoarseVector and rfineVector!
!</description>
  
!<input>
  ! The t_interlevelProjectionBlock structure that configures the grid transfer
  TYPE(t_interlevelProjectionBlock), INTENT(IN) :: rprojection 
  
  ! Fine grid vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: rfineVector
!</input>
  
!<inputoutput>
  ! Temporary vector. The vector must have the same data type as
  ! rcoarseVector and must be at least as long as indicated by the function
  ! mlprj_getTempMemory. If mlprj_getTempMemory was =0, the vector may
  ! be a dummy vector.
  ! The vector does not have to be connected to a discretisation structure
  ! or something similar; the content is undefined at entry and will be
  ! undefined when leaving this routine.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector

  ! Coarse grid vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: rcoarseVector
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrCoarse,p_rdiscrFine
  TYPE(t_interlevelProjectionScalar) :: ractProjection
  TYPE(t_triangulation), POINTER :: p_rtriaCoarse,p_rtriaFine
  
  ! Pointers into the triangulation
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementCoarse
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementFine
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementCoarse
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementFine
  
  ! Data arrays
  REAL(DP), DIMENSION(:), POINTER :: p_DuCoarse, p_DuFine
 
    ! The vectors must be of data type DOUBLE - we don't support anything
    ! different at the moment...
    IF (rcoarseVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Coarse grid vector has unsupported data type!'
      STOP
    END IF

    IF (rfineVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Fine grid vector has unsupported data type!'
      STOP
    END IF
    
    ! Calls the correct prolongation routine for each block in the 
    ! discretisation...
    DO i=1,rcoarseVector%nblocks
    
      IF (rcoarseVector%RvectorBlock(i)%NEQ .GT. 0) THEN
        p_rdiscrCoarse => rcoarseVector%RvectorBlock(i)%p_rspatialDiscretisation
        p_rdiscrFine => rfineVector%RvectorBlock(i)%p_rspatialDiscretisation
        
        ! We need a discretisation:
        IF ((.NOT. ASSOCIATED(p_rdiscrCoarse)) .OR. &
            (.NOT. ASSOCIATED(p_rdiscrFine))) THEN
          PRINT *,'Intergrid transfer: No discretisation!'
        END IF

        ! Currently, we support only uniform triangulations.
        IF ((p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
            (p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM)) THEN
          PRINT *,'Intergrid transfer supports currently only uniform discretisations!'
          STOP
        END IF
        
        ! Get the pointers to the vectors
        CALL storage_getbase_double (rcoarseVector%RvectorBlock(1)%h_Ddata,p_DuCoarse)
        p_DuCoarse => p_DuCoarse( &
          rcoarseVector%RvectorBlock(1)%iidxFirstEntry: &
          rcoarseVector%RvectorBlock(1)%iidxFirstEntry+rcoarseVector%RvectorBlock(1)%NEQ)

        CALL storage_getbase_double (rfineVector%RvectorBlock(1)%h_Ddata,p_DuFine)
        p_DuFine => p_DuFine( &
          rfineVector%RvectorBlock(1)%iidxFirstEntry: &
          rfineVector%RvectorBlock(1)%iidxFirstEntry+rfineVector%RvectorBlock(1)%NEQ)
        
        ! Use the first projection structure as template and create
        ! the actual projection structure for our situation.
        ! Remember, in a uniform grid we only have one projection structure
        ! and one element distribution!
        CALL mlprj_getProjectionStrategy (rprojection%RscalarProjection(1,1), &
              p_rdiscrCoarse%RelementDistribution(1), &
              p_rdiscrFine%RelementDistribution(1), &
              ractProjection)
      
        ! Depending on the element type of the trial functions in the
        ! discretisation, choose the right prolongation and call it.
        p_rtriaCoarse => p_rdiscrCoarse%p_rtriangulation
        p_rtriaFine => p_rdiscrFine%p_rtriangulation
        SELECT CASE (ractProjection%ielementTypeProlongation)
        CASE (EL_P1)
          ! P1 restriction
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_restUniformP1_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementFine,p_IneighboursAtElementFine, &
               p_rtriaCoarse%NEL,p_rtriaFine%NEL)
               
        CASE (EL_Q0)
          ! Q0 restriction
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_restUniformQ0_double (p_DuCoarse,p_DuFine, &
               p_IneighboursAtElementFine, &
               p_rtriaCoarse%NEL,p_rtriaFine%NEL)
               
        CASE (EL_Q1)
          ! Q1 restriction
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_restUniformQ1_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementFine,p_IneighboursAtElementFine,&
               p_rtriaCoarse%NEL,p_rtriaFine%NEL)
               
        CASE DEFAULT
          PRINT *,'Unsupported restriction!'
          STOP
        END SELECT
      
      END IF

    END DO  ! i

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_performInterpolation (rprojection,rcoarseVector, &
                                         rfineVector,rtempVector)
  
!<description>
  ! Performs an interpolation for a given block vector (i.e. a projection
  ! in the primal space where the solution vector lives). The vector
  ! rfineVector on a finer grid is projected to the vector
  ! rcoarseVector on a coarser grid. 
  ! rprojection configures how the grid transfer is performed.
  ! This projection structure must be build corresponding to the spatial 
  ! discretisation structures in rcoarseVector and rfineVector!
!</description>
  
!<input>
  
  ! The t_interlevelProjectionBlock structure that configures the grid transfer
  TYPE(t_interlevelProjectionBlock), INTENT(IN) :: rprojection 
  
  ! Fine grid vector
  TYPE(t_vectorBlock), INTENT(IN) :: rfineVector
!</input>
  
!<inputoutput>
  ! Temporary vector. The vector must have the same data type as
  ! rcoarseVector and must be at least as long as indicated by the function
  ! mlprj_getTempMemory. If mlprj_getTempMemory was =0, the vector may
  ! be a dummy vector.
  ! The vector does not have to be connected to a discretisation structure
  ! or something similar; the content is undefined at entry and will be
  ! undefined when leaving this routine.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector

  ! Coarse grid vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: rcoarseVector
!</inputoutput>
  
!</subroutine>
  
  ! local variables
  INTEGER :: i
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrCoarse,p_rdiscrFine
  TYPE(t_interlevelProjectionScalar) :: ractProjection
  TYPE(t_triangulation), POINTER :: p_rtriaCoarse,p_rtriaFine
  
  ! Pointers into the triangulation
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementCoarse
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementFine
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementCoarse
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementFine
  
  ! Data arrays
  REAL(DP), DIMENSION(:), POINTER :: p_DuCoarse, p_DuFine
 
    ! The vectors must be of data type DOUBLE - we don't support anything
    ! different at the moment...
    IF (rcoarseVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Coarse grid vector has unsupported data type!'
      STOP
    END IF

    IF (rfineVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Fine grid vector has unsupported data type!'
      STOP
    END IF
    
    ! Calls the correct prolongation routine for each block in the 
    ! discretisation...
    DO i=1,rcoarseVector%nblocks
    
      IF (rcoarseVector%RvectorBlock(i)%NEQ .GT. 0) THEN
        p_rdiscrCoarse => rcoarseVector%RvectorBlock(i)%p_rspatialDiscretisation
        p_rdiscrFine => rfineVector%RvectorBlock(i)%p_rspatialDiscretisation
        
        ! We need a discretisation:
        IF ((.NOT. ASSOCIATED(p_rdiscrCoarse)) .OR. &
            (.NOT. ASSOCIATED(p_rdiscrFine))) THEN
          PRINT *,'Intergrid transfer: No discretisation!'
        END IF

        ! Currently, we support only uniform triangulations.
        IF ((p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
            (p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM)) THEN
          PRINT *,'Intergrid transfer supports currently only uniform discretisations!'
          STOP
        END IF
        
        ! Get the pointers to the vectors
        CALL storage_getbase_double (rcoarseVector%RvectorBlock(1)%h_Ddata,p_DuCoarse)
        p_DuCoarse => p_DuCoarse( &
          rcoarseVector%RvectorBlock(1)%iidxFirstEntry: &
          rcoarseVector%RvectorBlock(1)%iidxFirstEntry+rcoarseVector%RvectorBlock(1)%NEQ)

        CALL storage_getbase_double (rfineVector%RvectorBlock(1)%h_Ddata,p_DuFine)
        p_DuFine => p_DuFine( &
          rfineVector%RvectorBlock(1)%iidxFirstEntry: &
          rfineVector%RvectorBlock(1)%iidxFirstEntry+rfineVector%RvectorBlock(1)%NEQ)
        
        ! Use the first projection structure as template and create
        ! the actual projection structure for our situation.
        ! Remember, in a uniform grid we only have one projection structure
        ! and one element distribution!
        CALL mlprj_getProjectionStrategy (rprojection%RscalarProjection(1,1), &
              p_rdiscrCoarse%RelementDistribution(1), &
              p_rdiscrFine%RelementDistribution(1), &
              ractProjection)
      
        ! Depending on the element type of the trial functions in the
        ! discretisation, choose the right prolongation and call it.
        p_rtriaCoarse => p_rdiscrCoarse%p_rtriangulation
        p_rtriaFine => p_rdiscrFine%p_rtriangulation
        SELECT CASE (ractProjection%ielementTypeProlongation)
        !CASE (EL_P1)
          ! P1 interpolation
        CASE (EL_Q0)
          ! Q0 interpolation
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_interpUniformQ0_double (p_DuCoarse,p_DuFine, &
               p_IneighboursAtElementFine, &
               p_rtriaCoarse%NEL,p_rtriaFine%NEL)
               
        !CASE (EL_Q1)
          ! Q1 interpolation
        CASE DEFAULT
          PRINT *,'Unsupported interpolation!'
          STOP
        END SELECT
      
      END IF

    END DO  ! i

  END SUBROUTINE

  ! ***************************************************************************
  ! Now the actual prolongation/restriction routines follow.
  ! For every type of element there is a set of routines that perform
  ! the actual grid transfer.
  ! We separate for
  ! - double and single precision vectors
  ! - different elements
  ! - uniform / conformal / ... discretisations
  ! ***************************************************************************

!!<subroutine>
!
!  SUBROUTINE mlprj_prolUniformQ1_double ()
!  
!!<description>
!  ! Prolongate a solution vector from a coarse grid to a fine grid.
!  ! Q1, uniform triangulation.
!!</description>
!  
!!<input>
!  
!!</input>
!  
!!<inputoutput>
!!</inputoutput>
!  
!!</subroutine>
!  
!  ! local variables
!
!  END SUBROUTINE

  ! ***************************************************************************
  ! Support for P1 element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformP1_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse,IverticesAtElementFine,&
               IneighboursAtElementCoarse,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! P1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
  
  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: Q2 = .5_DP
  INTEGER(PREC_ELEMENTIDX) :: iel
  REAL(DP) :: duh1,duh2,duh3

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    CALL lalg_vectorCopyDble (DuCoarse,DuFine(1:SIZE(DuCoarse)))

    ! Loop over the elements
    DO iel=1,NELCoarse

      duh1=DuCoarse(IverticesAtElementCoarse(1,iel))
      duh2=DuCoarse(IverticesAtElementCoarse(2,iel))
      duh3=DuCoarse(IverticesAtElementCoarse(3,iel))

      ! Now check on every of the edges, if we already computed
      ! the value in the midpoint: Compute only if the neighbour
      ! element has smaller number.
      IF (IneighboursAtElementCoarse(1,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(1,iel)) = Q2*(duh1+duh2)

      IF (IneighboursAtElementCoarse(2,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(2,iel)) = Q2*(duh2+duh3)

      IF (IneighboursAtElementCoarse(3,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(3,iel)) = Q2*(duh3+duh1)

    END DO

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformP1_double (DuCoarse,DuFine, &
               IverticesAtElementFine,IneighboursAtElementFine,&
               NELcoarse,NELfine)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! P1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
  
  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse

  ! Number of elements in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELfine
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: Q2 = .5_DP
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_POINTIDX) :: ih1
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    CALL lalg_vectorCopyDble (DuFine(1:SIZE(DuCoarse)),DuCoarse)
    
    ! Loop over the elements to collect the missing additive contributions:
    DO iel=NELcoarse+1,NELfine
      IH1 = IverticesAtElementFine(1,iel)

      ! Treat every edge only once:
      
      IF (IneighboursAtElementFine(1,iel) .LT. iel) &
        DuCoarse(IH1) = DuCoarse(IH1) + Q2*DuFine(IverticesAtElementFine(2,iel))
        
      IF (IneighboursAtElementFine(3,iel) .LT. iel) &
        DuCoarse(IH1) = DuCoarse(IH1) + Q2*DuFine(IverticesAtElementFine(3,iel))

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  ! Support for Q0 element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformQ0_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, &
               NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! Q0, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: Q2 = .5_DP
  REAL(DP), PARAMETER :: Q4 = .25_DP
  
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_POINTIDX) :: IELH1,IELH2,IELH3,IELH4
  REAL(DP) :: duh

    ! Loop over the elements
    DO iel=1,NELCoarse

      ! Get the four child elements of the coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Put the value on the coarse grid into all four child
      ! elements
      duh = DuCoarse(iel)
      DuFine(IELH1) = duh
      DuFine(IELH2) = duh
      DuFine(IELH3) = duh
      DuFine(IELH4) = duh
    END DO

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformQ0_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, &
               NELcoarse, NELfine)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! Q0, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse

  ! Number of elements in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELfine
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_ELEMENTIDX) :: IELH1,IELH2,IELH3,IELH4
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    
    ! Loop over the elements to collect the missing additive contributions:
    DO iel=NELcoarse+1,NELfine
    
      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)
      
      ! Add the values in these nodes together to get the
      ! value in the coarse grid element
      DuCoarse(iel)= DuFine(IELH1)+DuFine(IELH2)+DuFine(IELH3)+DuFine(IELH4)
      
    END DO
    
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformQ0_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse, NELfine)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! Q0, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse

  ! Number of elements in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELfine
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: Q4 = .25_DP
  
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_ELEMENTIDX) :: IELH1,IELH2,IELH3,IELH4
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    
    ! Loop over the elements to collect the missing additive contributions:
    DO iel=NELcoarse+1,NELfine
    
      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)
      
      ! From all four values, build the mean and use that as
      ! value of the coarse grid element
      DuCoarse(iel)= Q4*(DuFine(IELH1)+DuFine(IELH2)+DuFine(IELH3)+DuFine(IELH4))
      
    END DO
    
  END SUBROUTINE
  
  ! ***************************************************************************
  ! Support for Q1 element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformQ1_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse,IverticesAtElementFine,&
               IneighboursAtElementCoarse,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! Q1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
  
  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: Q2 = .5_DP
  REAL(DP), PARAMETER :: Q4 = .25_DP
  
  INTEGER(PREC_ELEMENTIDX) :: iel
  REAL(DP) :: duh1,duh2,duh3,duh4

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    CALL lalg_vectorCopyDble (DuCoarse,DuFine(1:SIZE(DuCoarse)))

    ! Loop over the elements
    DO iel=1,NELCoarse

      duh1=DuCoarse(IverticesAtElementCoarse(1,iel))
      duh2=DuCoarse(IverticesAtElementCoarse(2,iel))
      duh3=DuCoarse(IverticesAtElementCoarse(3,iel))
      duh4=DuCoarse(IverticesAtElementCoarse(4,iel))

      ! Now check on every of the edges, if we already computed
      ! the value in the midpoint: Compute only if the neighbour
      ! element has smaller number.
      IF (IneighboursAtElementCoarse(1,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(2,iel)) = Q2*(duh1+duh2)

      IF (IneighboursAtElementCoarse(2,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(2,iel)) = Q2*(duh2+duh3)

      IF (IneighboursAtElementCoarse(3,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(2,iel)) = Q2*(duh3+duh4)

      IF (IneighboursAtElementCoarse(4,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(2,iel)) = Q2*(duh4+duh1)
        
      ! Don't forget the DOF in the midpoint of the element
      DuFine(IverticesAtElementFine(3,iel)) = Q4*(duh1+duh2+duh3+duh4)

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformQ1_double (DuCoarse,DuFine, &
               IverticesAtElementFine,IneighboursAtElementFine,&
               NELcoarse,NELfine)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
  
  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse

  ! Number of elements in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELfine
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: Q2 = .5_DP
  REAL(DP), PARAMETER :: Q4 = .25_DP
  
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_POINTIDX) :: i1,i2,i3,i4
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    CALL lalg_vectorCopyDble (DuFine(1:SIZE(DuCoarse)),DuCoarse)
    
    ! Loop over the elements to collect the missing additive contributions:
    DO iel=NELcoarse+1,NELfine
      i1=IverticesAtElementFine(1,iel)
      i2=IverticesAtElementFine(2,iel)
      i3=IverticesAtElementFine(3,iel)
      i4=IverticesAtElementFine(4,iel)

      ! Additive contribution of the midpoint
      DuCoarse(i1) = DuCoarse(i1)+Q4*(DuFine(i2)+DuFine(i3)+DuFine(i4))

      ! Treat every edge only once:
      IF (IneighboursAtElementFine(1,iel) .EQ. 0) &
        DuCoarse(i1) = DuCoarse(i1)+Q4*DuFine(i2)
      IF (IneighboursAtElementFine(4,iel) .EQ. 0) &
        DuCoarse(i1) = DuCoarse(i1)+Q4*DuFine(i4)
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************
  ! Support for Q2 element
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_prolUniformQ2_double ()
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $Q_2$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  
!</input>
  
!<inputoutput>
!</inputoutput>
  
!</subroutine>
  
  ! local variables

  PRINT *,'not implemented!'
  STOP

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformQ2_double ()
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $Q_2$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  
!</input>
  
!<inputoutput>
!</inputoutput>
  
!</subroutine>
  
  ! local variables

  PRINT *,'not implemented!'
  STOP

  END SUBROUTINE
  
END MODULE
