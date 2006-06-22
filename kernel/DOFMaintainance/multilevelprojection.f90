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
!# 1.) mlprj_initProjectionDirect /
!#     mlprj_initProjectionVec /
!#     mlprj_initProjectionMat
!#     -> Initialises a projection structure with values according to a
!#        spatial discretisation.
!#        Uses a discretisation structure, a (block) vector or a (block)
!#        matrix on the fine/coarse grid as template.
!#
!# 2.) mlprj_doneProjection
!#     -> Cleans up a projection structure.
!#
!# 3.) mlprj_getTempMemoryDirect /
!#     mlprj_getTempMemoryVec /
!#     mlprj_getTempMemoryMat
!#     -> Determine the amount of temporary storage that is necessary
!#        to transport a vector from one level to another
!#        Uses a discretisation structure, a (block) vector or a (block)
!#        matrix on the fine/coarse grid as template.
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
  USE geometryaux
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
    !    prolongation, e.g. when $Q_2$ is used, apart from -1 the following 
    !    values are allowed:
    ! 0=constant prolongation
    ! 1=linear prolongation of a once refined mesh
    ! 2=quadratic interpolation
    INTEGER                     :: iprolongationOrder = -1
    
    ! Prolongation variant for nonconforming elements.
    ! Allows for the elements E030/EM30/E031/EM31 to switch to a special-type 
    ! prolongation.
    ! = 0: Use default prolongation,
    ! = 1: Use standard prolongation, equally weighted (1/2 from left, 1/2 from right),
    ! = 2: Use extended prolongation, equally weighted (1/2 from left, 1/2 from right),
    ! = 3: Use extended prolongation, weighted by element size (L2 projection),
    ! = 4: Use extended prolongation, weighted by element size of neighbour element
    ! To activate extended prolongation, set this to >= 2 after initialising the
    ! interlevel projection structure!
    INTEGER                     :: iprolEX3Yvariant = 0
    
    ! Configuration parameter for extended prolongation of E030/EM30/E031/EM31
    ! element. Only valid if iprolEX3Yvariant >= 2.
    ! Aspect-ratio indicator; controls switching to constant prolongation.
    ! <=1: switch depending on aspect ratio of current element (standard),
    !  =2: switch depending on aspect ratio of current element and
    !      neighbour element
    INTEGER                     :: iprolARIndicatorEX3Y = 1
    
    ! Configuration parameter for extended prolongation of E030/EM30/E031/EM31
    ! element. Only valid if iprolEX3Yvariant >= 2.
    ! Upper bound aspect ratio; for all elements with higher AR
    ! the prolongation is switched to constant prolongation .
    ! This is set to 20.0 by default according to the analysis in
    ! [Michael Köster, Robuste Mehrgitter-Krylowraum-Techniken für FEM-Verfahren,
    !  2004, Diploma-Theses, Chair of Mathematics, University of Dortmund]
    REAL(DP)                    :: dprolARboundEX3Y = 20.0_DP

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
    !    restriction, e.g. when $Q_2$ is used, apart from -1 the following 
    !    values are allowed:
    ! 0=constant restriction
    ! 1=linear restriction of a once refined mesh
    ! 2=quadratic restriction
    INTEGER                     :: irestrictionOrder = -1

    ! Restriction variant for nonconforming elements.
    ! Allows for the elements E030/EM30/E031/EM31 to switch to a special-type 
    ! prolongation.
    ! = 0: Use default restriction,
    ! = 1: Use standard restriction, equally weighted (1/2 from left, 1/2 from right),
    ! = 2: Use extended restriction, equally weighted (1/2 from left, 1/2 from right),
    ! = 3: Use extended restriction, weighted by element size (L2 projection),
    ! = 4: Use extended restriction, weighted by element size of neighbour element
    ! To activate extended prolongation, set this to >= 2 after initialising the
    ! interlevel projection structure!
    INTEGER                     :: irestEX3Yvariant = 0
    
    ! Configuration parameter for extended restriction of E030/EM30/E031/EM31
    ! element. Only valid if iprolEX3Yvariant >= 2.
    ! Aspect-ratio indicator; controls switching to constant prolongation.
    ! <=1: switch depending on aspect ratio of current element (standard),
    !  =2: switch depending on aspect ratio of current element and
    !      neighbour element
    INTEGER                     :: irestARIndicatorEX3Y = 1
    
    ! Configuration parameter for extended restriction of E030/EM30/E031/EM31
    ! element. Only valid if iprolEX3Yvariant >= 2.
    ! Upper bound aspect ratio; for all elements with higher AR
    ! the prolongation is switched to constant prolongation .
    ! This is set to 20.0 by default according to the analysis in
    ! [Michael Köster, Robuste Mehrgitter-Krylowraum-Techniken für FEM-Verfahren,
    !  2004, Diploma-Theses, Chair of Mathematics, University of Dortmund]
    REAL(DP)                    :: drestARboundEX3Y = 20.0_DP

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
    !    interpolation, e.g. when $Q_2$ is used, apart from -1 the following 
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

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_initProjectionDirect (rprojection,RspatialDiscretisation)
  
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

  SUBROUTINE mlprj_initProjectionVec (rprojection,rvector)
  
!<description>
  ! This subroutine initialises an t_interlevelProjectionBlock with default
  ! values for a given PDE system. This allows the interlevel-projection
  ! routines to calculate constants that are necessary during a level
  ! change. (This is used e.g. in the general prolongation/restriction
  ! where the prolongation/restriction matrix must be calculated).
  !
  ! The calculated information is saved to rprojection and released
  ! with mlprj_doneProjection.
  !
  ! The PDE system is specified by the p_rspatialDiscretisation structures
  ! saved in the given vector.
!</description>

!<input>
  ! A vector containing information about the spatial discretisation of
  ! the given PDE.
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
!</input>
  
!<output>
  ! A t_interlevelProjectionBlock structure that will be filled with data
  ! about the projection of all the equations described by RspatialDiscretisation.
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
!</output>
  
!</subroutine>

    ! local variables
    TYPE(t_spatialDiscretisation), DIMENSION(LSYSBL_MAXBLOCKS) :: Rdiscr
    INTEGER :: i

    IF (rvector%nblocks .EQ. 0) THEN
      PRINT *,'mlprj_initProjectionVec: No discretisation!'
      STOP
    END IF

    ! Set up an array of discretisation structures for all the equations
    DO i=1,rvector%nblocks
      Rdiscr(i) = &
        rvector%RvectorBlock(i)%p_rspatialDiscretisation
    END DO

    ! Call the standard initialisation routine
    CALL mlprj_initProjectionDirect (rprojection,Rdiscr)

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_initProjectionMat (rprojection,rmatrix)
  
!<description>
  ! This subroutine initialises an t_interlevelProjectionBlock with default
  ! values for a given PDE system. This allows the interlevel-projection
  ! routines to calculate constants that are necessary during a level
  ! change. (This is used e.g. in the general prolongation/restriction
  ! where the prolongation/restriction matrix must be calculated).
  !
  ! The calculated information is saved to rprojection and released
  ! with mlprj_doneProjection.
  !
  ! The PDE system is specified by the p_rspatialDiscretisation structures
  ! saved in the given matrix.
!</description>

!<input>
  ! A matrix containing information about the spatial discretisation of
  ! the given PDE.
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
!</input>
  
!<output>
  ! A t_interlevelProjectionBlock structure that will be filled with data
  ! about the projection of all the equations described by RspatialDiscretisation.
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
!</output>
  
!</subroutine>

    ! local variables
    TYPE(t_spatialDiscretisation), DIMENSION(LSYSBL_MAXBLOCKS) :: Rdiscr
    INTEGER :: i,j

    IF (rmatrix%ndiagBlocks .EQ. 0) THEN
      PRINT *,'mlprj_initProjectionMat: No discretisation!'
      STOP
    END IF

    ! Set up an array of discretisation structures for all the equations.
    ! In every 'column' of the block matrix, search for the first existing
    ! matrix and use its properties for initialisation
    DO i=1,rmatrix%ndiagBlocks
      DO j=1,rmatrix%ndiagBlocks
        IF (rmatrix%RmatrixBlock(j,i)%NEQ .NE. 0) THEN
          Rdiscr(i) = &
            rmatrix%RmatrixBlock(j,i)%p_rspatialDiscretisation
          EXIT
        END IF
      END DO
    END DO

    ! Call the standard initialisation routine
    CALL mlprj_initProjectionDirect (rprojection,Rdiscr)

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

  ! Set default prolongation/restriction for Ex3y-type elements 
  ! if not specified
  IF (ractProjection%iprolEX3Yvariant .EQ. 0) THEN
    ractProjection%iprolEX3Yvariant = 1
  END IF

  IF (ractProjection%irestEX3Yvariant .EQ. 0) THEN
    ractProjection%irestEX3Yvariant = 1
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

  INTEGER(PREC_VECIDX) FUNCTION mlprj_getTempMemoryDirect (rprojection, &
                                                     RdiscrCoarse,RdiscrFine)
  
!<description>
  ! Returns the amount of temporary memory that is needed by the
  ! interlevel-routines to transfer vectors between two grids.
  ! RdiscrCoarse and RdiscrFine are arrays of discretisation structures.
  ! Each discretisation structure corresponds to one scalar equation.
!</description>
  
!<input>
  ! Projection structure that configures the grid transfer for all equations
  ! and all element distributions in each equation.
  TYPE(t_interlevelProjectionBlock), INTENT(IN) :: rprojection 
  
  ! List of disretisation structures for the equations on the Coarse grid 
  TYPE(t_spatialDiscretisation), DIMENSION(:), INTENT(IN) :: RdiscrCoarse

  ! List of disretisation structures for the equations on the Fine grid 
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
  mlprj_getTempMemoryDirect = imemmax
  
  END FUNCTION
  
  ! ***************************************************************************
  
!<function>

  INTEGER(PREC_VECIDX) FUNCTION mlprj_getTempMemoryVec (rprojection, &
                                                        rvectorCoarse,rvectorFine)
  
!<description>
  ! Returns the amount of temporary memory that is needed by the
  ! interlevel-routines to transfer vectors between two grids.
  ! rvectorCoarse and rvectorFine define template vectors on the coarse
  ! and fine grid, respectively; the memory is calculated using
  ! the discretisation structures associated to these.
!</description>
  
!<input>
  ! Projection structure that configures the grid transfer for all equations
  ! and all element distributions in each equation.
  TYPE(t_interlevelProjectionBlock), INTENT(IN) :: rprojection 
  
  ! Coarse grid vector. Must have the discretisation structure of the
  ! coarse grid attached.
  TYPE(t_vectorBlock), INTENT(IN) :: rvectorCoarse

  ! Fine grid vector. Must have the discretisation structure of the
  ! Fine grid attached.
  TYPE(t_vectorBlock), INTENT(IN) :: rvectorFine
!</input>

!<result>
  ! Length of a temporary vector that is necessary for transferring vectors
  ! on the coarse grid to the fine grid and vice versa.
  ! =0 if no memory is necessary.
!</result>

    ! local variables
    TYPE(t_spatialDiscretisation), DIMENSION(LSYSBL_MAXBLOCKS) :: RdiscrCoarse
    TYPE(t_spatialDiscretisation), DIMENSION(LSYSBL_MAXBLOCKS) :: RdiscrFine
    INTEGER :: i

    IF ((rvectorCoarse%nblocks .EQ. 0) .OR. (rvectorFine%nblocks .EQ. 0)) THEN
      PRINT *,'mlprj_getTempMemoryVec: No discretisation!'
      STOP
    END IF

    ! Set up an array of discretisation structures for all the equations
    DO i=1,rvectorCoarse%nblocks
      RdiscrCoarse(i) = &
        rvectorCoarse%RvectorBlock(i)%p_rspatialDiscretisation
    END DO

    DO i=1,rvectorFine%nblocks
      RdiscrFine(i) = &
        rvectorFine%RvectorBlock(i)%p_rspatialDiscretisation
    END DO
      
    ! Call the standard getTempMemory routine
    mlprj_getTempMemoryVec = &
      mlprj_getTempMemoryDirect (rprojection, RdiscrCoarse,RdiscrFine)

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  INTEGER(PREC_VECIDX) FUNCTION mlprj_getTempMemoryMat (rprojection, &
                                                        rmatrixCoarse,rmatrixFine)
  
!<description>
  ! Returns the amount of temporary memory that is needed by the
  ! interlevel-routines to transfer vectors between two grids.
  ! rmatrixCoarse and rmatrixFine define template matrices on the coarse
  ! and fine grid, respectively; the memory is calculated using
  ! the discretisation structures associated to these.
!</description>
  
!<input>
  ! Projection structure that configures the grid transfer for all equations
  ! and all element distributions in each equation.
  TYPE(t_interlevelProjectionBlock), INTENT(IN) :: rprojection 
  
  ! Coarse grid matrix. Must have the discretisation structure of the
  ! coarse grid attached.
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrixCoarse

  ! Fine grid matrix. Must have the discretisation structure of the
  ! Fine grid attached.
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrixFine
!</input>

!<result>
  ! Length of a temporary vector that is necessary for transferring vectors
  ! on the coarse grid to the fine grid and vice versa.
  ! =0 if no memory is necessary.
!</result>

    ! local variables
    TYPE(t_spatialDiscretisation), DIMENSION(LSYSBL_MAXBLOCKS) :: RdiscrCoarse
    TYPE(t_spatialDiscretisation), DIMENSION(LSYSBL_MAXBLOCKS) :: RdiscrFine
    INTEGER :: i,j

    IF ((rmatrixCoarse%ndiagBlocks .EQ. 0) .OR. (rmatrixFine%ndiagBlocks .EQ. 0)) THEN
      PRINT *,'mlprj_getTempMemoryVec: No discretisation!'
      STOP
    END IF

    ! Set up an array of discretisation structures for all the equations
    DO i=1,rmatrixCoarse%ndiagBlocks
      DO j=1,rmatrixCoarse%ndiagBlocks
        IF (rmatrixCoarse%RmatrixBlock(j,i)%NEQ .NE. 0) THEN
          RdiscrCoarse(i) = &
            rmatrixCoarse%RmatrixBlock(j,i)%p_rspatialDiscretisation
          EXIT
        END IF
      END DO
    END DO

    DO i=1,rmatrixFine%ndiagBlocks
      DO j=1,rmatrixFine%ndiagBlocks
        IF (rmatrixFine%RmatrixBlock(j,i)%NEQ .NE. 0) THEN
          RdiscrFine(i) = &
            rmatrixFine%RmatrixBlock(j,i)%p_rspatialDiscretisation
        END IF
      END DO
    END DO
      
    ! Call the standard getTempMemory routine
    mlprj_getTempMemoryMat = &
      mlprj_getTempMemoryDirect (rprojection, RdiscrCoarse,RdiscrFine)

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
  TYPE(t_vectorBlock), INTENT(INOUT) :: rcoarseVector
!</input>
  
!<inputoutput>
  ! Temporary vector. The vector must have the same data type as
  ! rfineVector and must be at least as long as indicated by the function
  ! mlprj_getTempMemory. If mlprj_getTempMemory was =0, the vector may
  ! be a dummy vector.
  ! The vector does not have to be connected to a discretisation structure
  ! or something similar; the content is undefined at entry and will be
  ! undefined when leaving this routine.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rtempVector
!</inputoutput>

!<output>
  ! Fine grid vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: rfineVector
!</output>
  
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
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementCoarse
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementFine
  REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordinatesCoarse
  REAL(DP), DIMENSION(:), POINTER   :: p_DelementAreaCoarse
  
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
    
    IF (lsysbl_isVectorSorted(rfineVector) .OR. lsysbl_isVectorSorted(rcoarseVector)) THEN
      PRINT *,'Vectors must be unsorted for level change!'
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
          STOP
        END IF

        ! Currently, we support only uniform triangulations.
        IF ((p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
            (p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM)) THEN
          PRINT *,'Intergrid transfer supports currently only uniform discretisations!'
          STOP
        END IF
        
        ! Get the pointers to the vectors
        CALL lsyssc_getbase_double (rcoarseVector%RvectorBlock(i),p_DuCoarse)
        CALL lsyssc_getbase_double (rfineVector%RvectorBlock(i),p_DuFine)
        
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
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_prolUniformQ1_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
               p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,p_rtriaCoarse%NEL)
               
        CASE (EL_E030,EL_EM30)
          ! Q1~ prolongation, DOF's = integral mean value
          CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                               p_IedgesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                               p_IedgesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
          
          ! Type of prolongation? Extended or not?
          SELECT CASE (ractProjection%iprolEX3Yvariant)
          CASE (:1) ! Standard prolongation
            CALL mlprj_prolUniformEx30_double (p_DuCoarse,p_DuFine, &
                p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)
                
          CASE (2:) ! Extended prolongation; modified weights, local switch
                    ! to constant prolongation
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                                       p_IverticesAtElementCoarse)
            CALL storage_getbase_double2d(p_rtriaCoarse%h_DcornerCoordinates, &
                                          p_DcornerCoordinatesCoarse)
            CALL storage_getbase_double(p_rtriaCoarse%h_DelementArea, &
                                        p_DelementAreaCoarse)
            ! (what a nasty call...)                                       
            CALL mlprj_prolUniformEx30ext_double (p_DuCoarse,p_DuFine, &
                    p_DcornerCoordinatesCoarse,p_IverticesAtElementCoarse, &
                    p_DelementAreaCoarse,&
                    p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                    p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                    p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
                    MIN(4,ractProjection%iprolEX3Yvariant)-2, &
                    ractProjection%dprolARboundEX3Y, &
                    ractProjection%iprolARIndicatorEX3Y)
          END SELECT
               
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
  TYPE(t_vectorScalar), INTENT(INOUT) :: rtempVector

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
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementCoarse
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementFine
  REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordinatesCoarse
  REAL(DP), DIMENSION(:), POINTER   :: p_DelementAreaCoarse
  
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
    
    IF (lsysbl_isVectorSorted(rfineVector) .OR. lsysbl_isVectorSorted(rcoarseVector)) THEN
      PRINT *,'Vectors must be unsorted for level change!'
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
          STOP
        END IF

        ! Currently, we support only uniform triangulations.
        IF ((p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
            (p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM)) THEN
          PRINT *,'Intergrid transfer supports currently only uniform discretisations!'
          STOP
        END IF
        
        ! Get the pointers to the vectors
        CALL lsyssc_getbase_double (rcoarseVector%RvectorBlock(i),p_DuCoarse)
        CALL lsyssc_getbase_double (rfineVector%RvectorBlock(i),p_DuFine)
        
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
               p_rtriaCoarse%NEL)
               
        CASE (EL_Q1)
          ! Q1 restriction
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_restUniformQ1_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementFine,p_IneighboursAtElementFine,&
               p_rtriaFine%NEL)
               
        CASE (EL_E030,EL_EM30)
          ! Q1~ restriction, DOF's = integral mean values
          CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                               p_IedgesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                               p_IedgesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
                               
          ! Type of restriction? Extended or not?
          SELECT CASE (ractProjection%iprolEX3Yvariant)
          CASE (:1) ! Standard prolongation
            CALL mlprj_restUniformEx30_double (p_DuCoarse,p_DuFine, &
                p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)          
                
          CASE (2:) ! Extended prolongation; modified weights, local switch
                    ! to constant prolongation
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                                       p_IverticesAtElementCoarse)
            CALL storage_getbase_double2d(p_rtriaCoarse%h_DcornerCoordinates, &
                                          p_DcornerCoordinatesCoarse)
            CALL storage_getbase_double(p_rtriaCoarse%h_DelementArea, &
                                        p_DelementAreaCoarse)
            ! (what a nasty call...)                                       
            CALL mlprj_restUniformEx30ext_double (p_DuCoarse,p_DuFine, &
                    p_DcornerCoordinatesCoarse,p_IverticesAtElementCoarse, &
                    p_DelementAreaCoarse,&
                    p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                    p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                    p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
                    MIN(4,ractProjection%iprolEX3Yvariant)-2, &
                    ractProjection%dprolARboundEX3Y, &
                    ractProjection%iprolARIndicatorEX3Y)
          END SELECT
               
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
  TYPE(t_vectorScalar), INTENT(INOUT) :: rtempVector

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
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementCoarse
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementFine
  
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
    
    IF (lsysbl_isVectorSorted(rfineVector) .OR. lsysbl_isVectorSorted(rcoarseVector)) THEN
      PRINT *,'Vectors must be unsorted for level change!'
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
          STOP
        END IF

        ! Currently, we support only uniform triangulations.
        IF ((p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
            (p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM)) THEN
          PRINT *,'Intergrid transfer supports currently only uniform discretisations!'
          STOP
        END IF
        
        ! Get the pointers to the vectors
        CALL lsyssc_getbase_double (rcoarseVector%RvectorBlock(i),p_DuCoarse)
        CALL lsyssc_getbase_double (rfineVector%RvectorBlock(i),p_DuFine)
        
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
          ! P1 interpolation
          CALL mlprj_interpUniformP1_double (p_DuCoarse,p_DuFine,p_rtriaCoarse%NVT)
          
        CASE (EL_Q0)
          ! Q0 interpolation
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_interpUniformQ0_double (p_DuCoarse,p_DuFine, &
               p_IneighboursAtElementFine, &
               p_rtriaCoarse%NEL,p_rtriaFine%NEL)
               
        CASE (EL_Q1)
          ! Q1 interpolation
          CALL mlprj_interpUniformQ1_double (p_DuCoarse,p_DuFine, p_rtriaCoarse%NVT)
          
        CASE (EL_E030,EL_EM30,EL_E031,EL_EM31)
          ! Q1~ interpolation, DOF's = integral mean values
          ! We use the same routine also for interpolating Ex31 solutions - there's
          ! not too much difference...
          CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                               p_IedgesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                               p_IedgesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
          CALL mlprj_interpUniformEx30_double (p_DuCoarse,p_DuFine, &
               p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
               p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
               p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)          
               
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
  ! Support for P0 element
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_prolUniformP0_double ()
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_0$, uniform triangulation, double precision vector.
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

  SUBROUTINE mlprj_restUniformP0_double ()
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_0$, uniform triangulation, double precision vector.
!</description>
  
!<input>
!</input>
  
!<output>
!</output>
  
!</subroutine>
  
  PRINT *,'not implemented.'
  STOP

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformP0_double ()
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_0$, uniform triangulation, double precision vector.
!</description>
  
!<input>
!</input>
  
!<output>
!</output>
  
!</subroutine>
  
  PRINT *,'not implemented.'
  STOP
    
  END SUBROUTINE

  ! ***************************************************************************
  ! Support for P1 element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformP1_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse,IverticesAtElementFine,&
               IneighboursAtElementCoarse,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_1$, uniform triangulation, double precision vector.
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
    CALL lalg_copyVectorDble (DuCoarse,DuFine(1:SIZE(DuCoarse)))

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
  ! $P_1$, uniform triangulation, double precision vector.
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
    CALL lalg_copyVectorDble (DuFine(1:SIZE(DuCoarse)),DuCoarse)
    
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
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformP1_double (DuCoarse,DuFine,NVTcoarse)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_POINTIDX), INTENT(IN) :: NVTcoarse
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  
    ! The first coase.NVT entries of the fine grid vector define 
    ! the values on the coarse grid - because of the two-level ordering!
    CALL lalg_copyVectorDble(DUfine(1:NVTcoarse),DUcoarse(1:NVTCoarse))
    
  END SUBROUTINE

  ! ***************************************************************************
  ! Support for P2 element
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_prolUniformP2_double ()
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_2$, uniform triangulation, double precision vector.
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

  SUBROUTINE mlprj_restUniformP2_double ()
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_2$, uniform triangulation, double precision vector.
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

  SUBROUTINE mlprj_interpUniformP2_double (DuCoarse,DuFine, &
                                           NVTcoarse, NMTcoarse)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_2$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_POINTIDX), INTENT(IN) :: NVTcoarse

  ! Number of edges in the coarse grid
  INTEGER(PREC_POINTIDX), INTENT(IN) :: NMTcoarse
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  
    ! The first coase.NVT+NMT entries of the fine grid vector define 
    ! the values on the coarse grid - because of the two-level ordering!
    CALL lalg_copyVectorDble(DUfine(1:NVTcoarse+NMTcoarse),&
                             DUcoarse(1:NVTCoarse+NMTcoarse))
    
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
  ! $Q_0$, uniform triangulation, double precision vector.
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
               NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $Q_0$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
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
    DO iel=1,NELcoarse
    
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
  ! $Q_0$, uniform triangulation, double precision vector.
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
               IneighboursAtElementCoarse,IneighboursAtElementFine,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $Q_1$, uniform triangulation, double precision vector.
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
  
  INTEGER(PREC_ELEMENTIDX) :: iel,ielh1,ielh2,ielh3,ielh4
  REAL(DP) :: duh1,duh2,duh3,duh4

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    CALL lalg_copyVectorDble (DuCoarse,DuFine(1:SIZE(DuCoarse)))

    ! Loop over the elements
    DO iel=1,NELCoarse

      duh1=DuCoarse(IverticesAtElementCoarse(1,iel))
      duh2=DuCoarse(IverticesAtElementCoarse(2,iel))
      duh3=DuCoarse(IverticesAtElementCoarse(3,iel))
      duh4=DuCoarse(IverticesAtElementCoarse(4,iel))

      ielh1=iel
      ielh2=IneighboursAtElementFine(2,ielh1)
      ielh3=IneighboursAtElementFine(2,ielh2)
      ielh4=IneighboursAtElementFine(2,ielh3)

      ! Now check on every of the edges, if we already computed
      ! the value in the midpoint: Compute only if the neighbour
      ! element has smaller number.
      IF (IneighboursAtElementCoarse(1,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(2,ielh1)) = Q2*(duh1+duh2)

      IF (IneighboursAtElementCoarse(2,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(2,ielh2)) = Q2*(duh2+duh3)

      IF (IneighboursAtElementCoarse(3,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(2,ielh3)) = Q2*(duh3+duh4)

      IF (IneighboursAtElementCoarse(4,iel) .LT. iel) &
        DuFine(IverticesAtElementFine(2,ielh4)) = Q2*(duh4+duh1)
        
      ! Don't forget the DOF in the midpoint of the element
      DuFine(IverticesAtElementFine(3,iel)) = Q4*(duh1+duh2+duh3+duh4)

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformQ1_double (DuCoarse,DuFine, &
               IverticesAtElementFine,IneighboursAtElementFine,&
               NELfine)
  
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
    CALL lalg_copyVectorDble (DuFine(1:SIZE(DuCoarse)),DuCoarse)
    
    ! Loop over the elements to collect the missing additive contributions:
    DO iel=1,NELfine
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
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformQ1_double (DuCoarse,DuFine, NVTcoarse)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_POINTIDX), INTENT(IN) :: NVTcoarse
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
    ! The first coase.NVT entries of the fine grid vector define the values
    ! on the coarse grid - because of the two-level ordering!
    CALL lalg_copyVectorDble(DUfine(1:NVTcoarse),DUcoarse(1:NVTCoarse))
    
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
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformQ2_double (DuCoarse,DuFine, &
                                           NVTcoarse, NMTcoarse, NELcoarse)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $Q_2$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_POINTIDX), INTENT(IN) :: NVTcoarse

  ! Number of edges in the coarse grid
  INTEGER(PREC_POINTIDX), INTENT(IN) :: NMTcoarse

  ! Number of elements in the coarse grid
  INTEGER(PREC_POINTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  
    ! The first coase.NVT+NMT+NEL entries of the fine grid vector define 
    ! the values on the coarse grid - because of the two-level ordering!
    CALL lalg_copyVectorDble(DUfine(1:NVTcoarse+NMTcoarse+NELcoarse),&
                             DUcoarse(1:NVTCoarse+NMTcoarse+NELcoarse))
    
  END SUBROUTINE
  
  ! ***************************************************************************
  ! Support for Q1~ element, DOF's = integral mean values in the edges
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformEx30_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E030/EM30, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IedgesAtElement array on the coarse grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementCoarse
  
  ! IedgesAtElement array on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine

  ! Number of vertices in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse

  ! Number of vertices in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP) :: DUH1,DUH2,DUH3,DUH4
  INTEGER(PREC_EDGEIDX) :: IM1,IM2, IM3,IM4, IA, IB, IC
  INTEGER(PREC_ELEMENTIDX) :: iel, IELH1, IELH2, IELH3, IELH4
  
  ! Weights for the restruction; all coefficients are halfed, so dividing
  ! by 2 is not necessary in the calculation routines.
  REAL(DP), PARAMETER :: A1=0.5_DP, A2=-0.125D0, A3=0_DP, A4=0.125D0
  REAL(DP), PARAMETER :: A5=0.625_DP, A6=0.125_DP, A7=0.125_DP, A8=0.125_DP
  
    ! Clear the output vector
    CALL lalg_clearVectorDble(DuFine)
  
    ! Loop over the coarse grid elements
    DO iel=1,NELcoarse

      ! Get the DOF's of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)-NVTcoarse
      IM2 = IedgesAtElementCoarse(2,iel)-NVTcoarse
      IM3 = IedgesAtElementCoarse(3,iel)-NVTcoarse
      IM4 = IedgesAtElementCoarse(4,iel)-NVTcoarse

      ! Get the values of the corresponding DOF's
      DUH1 = DuCoarse(IM1)
      DUH2 = DuCoarse(IM2)
      DUH3 = DuCoarse(IM3)
      DUH4 = DuCoarse(IM4)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Distribute the value at the edge IM1 to the 
      ! corresponding fine inner nodes

      IF (IneighboursAtElementCoarse(1,iel).NE.0) THEN 
        ! There is a neighbour at the edge
        IA=IedgesAtElementFine(1,IELH1)-NVTfine
        IB=IedgesAtElementFine(4,IELH2)-NVTfine
        DuFine(IA)=DuFine(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
        DuFine(IB)=DuFine(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
      ELSE
        ! No neighbour; boundary element
        IA=IedgesAtElementFine(1,IELH1)-NVTfine
        IB=IedgesAtElementFine(4,IELH2)-NVTfine
        DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
        DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
      ENDIF
      IC=IedgesAtElementFine(2,IELH1)-NVTfine
      DuFine(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3

      ! Distribute the value at the edge IM2 to the 
      ! corresponding fine inner nodes

      IF (IneighboursAtElementCoarse(2,iel).NE.0) THEN 
        ! There is a neighbour at the edge
       IA=IedgesAtElementFine(1,IELH2)-NVTfine
       IB=IedgesAtElementFine(4,IELH3)-NVTfine
       DuFine(IA)=DuFine(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DuFine(IB)=DuFine(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
      ELSE
        ! No neighbour; boundary element
       IA=IedgesAtElementFine(1,IELH2)-NVTfine
       IB=IedgesAtElementFine(4,IELH3)-NVTfine
       DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
      ENDIF
      IC=IedgesAtElementFine(2,IELH2)-NVTfine
      DuFine(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4

      ! Distribute the value at the edge IM3 to the 
      ! corresponding fine inner nodes

      IF (IneighboursAtElementCoarse(3,iel).NE.0) THEN 
        ! There is a neighbour at the edge
       IA=IedgesAtElementFine(1,IELH3)-NVTfine
       IB=IedgesAtElementFine(4,IELH4)-NVTfine
       DuFine(IA)=DuFine(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DuFine(IB)=DuFine(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
      ELSE
        ! No neighbour; boundary element
       IA=IedgesAtElementFine(1,IELH3)-NVTfine
       IB=IedgesAtElementFine(4,IELH4)-NVTfine
       DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
      ENDIF
      IC=IedgesAtElementFine(2,IELH3)-NVTfine
      DuFine(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1

      ! Distribute the value at the edge IM4 to the 
      ! corresponding fine inner nodes

      IF (IneighboursAtElementCoarse(4,iel).NE.0) THEN 
        ! There is a neighbour at the edge
        IA=IedgesAtElementFine(1,IELH4)-NVTfine
        IB=IedgesAtElementFine(4,IELH1)-NVTfine
        DuFine(IA)=DuFine(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
        DuFine(IB)=DuFine(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
      ELSE
        ! No neighbour; boundary element
        IA=IedgesAtElementFine(1,IELH4)-NVTfine
        IB=IedgesAtElementFine(4,IELH1)-NVTfine
        DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
        DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
      END IF
      IC=IedgesAtElementFine(2,IELH4)-NVTfine
      DuFine(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2

    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformEx30ext_double (DuCoarse,DuFine, &
               DcornerCoordinatesCoarse,IverticesAtElementCoarse,DelementAreaCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NELcoarse, &
               iweightingType, daspectRatioBound, iarIndicator)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E030/EM30, uniform triangulation, double precision vector.
  !
  ! Extended version. Switch to constant prolongation if aspect ratio
  ! of an element is too large.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse

  ! DcornerCoordinates array on the coarse grid
  REAL(DP), DIMENSION(:,:), INTENT(IN)                :: DcornerCoordinatesCoarse

  ! IverticesAtElement array on the coarse grid
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse
  
  ! DelementArea array on the coarse grid
  REAL(DP), DIMENSION(:), INTENT(IN)                  :: DelementAreaCoarse

  ! IedgesAtElement array on the coarse grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementCoarse
  
  ! IedgesAtElement array on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse

  ! Number of vertices in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
  
  ! Type of the averaging on the element edges
  ! <=0: standard averaging of both contributions by 1/2,        
  !  =1: weighted averaging of the interpolated function values:
  !      The area of the current coarse grid element determines 
  !      the weight. (L2-projection, standard),
  !  =2: weighted averaging of the interpolated function values:
  !      The area of the neightbour element of the coarse grid 
  !      the weight. 
  INTEGER, INTENT(IN)  :: iweightingType
  
  ! Upper bound aspect ratio; for all elements with higher AR
  ! the prolongation is switched to constant prolongation 
  REAL(DP), INTENT(IN) :: daspectRatioBound
  
  ! Aspect-ratio indicator.
  ! Controls switching to constant prolongation.
  ! <=1: switch depending on aspect ratio of current element,
  !  =2: switch depending on aspect ratio of current element and
  !      neighbour element
  INTEGER, INTENT(IN)  :: iarIndicator
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP) :: DUH1,DUH2,DUH3,DUH4
  REAL(DP), DIMENSION(0:TRIA_MAXNME2D) :: daspectRatio,darea
  REAL(DP), DIMENSION(TRIA_MAXNME2D) :: dweight
  INTEGER(PREC_EDGEIDX) :: IM1,IM2, IM3,IM4, IA, IB, IC
  INTEGER(PREC_ELEMENTIDX) :: iel, IELH1, IELH2, IELH3, IELH4
  INTEGER(PREC_ELEMENTIDX), DIMENSION(0:TRIA_MAXNME2D) :: IELA
  INTEGER :: i
  INTEGER, DIMENSION(TRIA_MAXNME2D) :: idoConstant
  REAL(DP), DIMENSION(NDIM2D,TRIA_MAXNVE2D) :: dcoords
  
  ! Weights for the prolongation.
  ! PRWEIG (.,1) gives the constants for the standard prolongation,
  ! PRWEIG (.,2) gives the constants for the constant prolongation.
  REAL(DP), DIMENSION(8,2), PARAMETER :: prweight = &
      RESHAPE((/1.0_DP, -0.25_DP, 0.0_DP, 0.25_DP, &
                0.625_DP, 0.125_DP, 0.125_DP, 0.125_DP, &
                1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP/),(/8,2/))
              
    ! Clear the output vector
    CALL lalg_clearVectorDble(DuFine)
  
    ! Loop over the coarse grid elements
    DO iel=1,NELcoarse

      ! Get the numbers of the elements that are neighbours to our current coarse
      ! grid element:
      !             +--------+
      !             |        |
      !             | IELA3  |
      !             |        |
      !    +--------4--------3--------+
      !    |        |        |        |
      !    | IELA4  |  IEL   | IELA2  |
      !    |        |        |        |
      !    +--------1--------2--------+
      !             |        |
      !             | IELA1  |
      !             |        |
      !             +--------+
      IELA(0) = iel
      IELA(1) = IneighboursAtElementCoarse(1,iel)
      IELA(2) = IneighboursAtElementCoarse(2,iel)
      IELA(3) = IneighboursAtElementCoarse(3,iel)
      IELA(4) = IneighboursAtElementCoarse(4,iel)
      
      ! For these five elements, determine the aspect ratio and their area.
      !
      ! At first the element in the center, which always exists.
      !
      ! Get the aspect ratio of the current coarse grid element;
      ! if necessary, calculate the reciprocal.
      dcoords = DcornerCoordinatesCoarse(:,IverticesAtElementCoarse(:,IELA(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      IF (daspectRatio(0) .LT. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)
      
      ! and the area of that element.
      darea(0) = DelementAreaCoarse(iel)
      
      ! Then the remaining neighbours.
      DO i=1,TRIA_MAXNME2D
        IF (IELA(i) .NE. 0) THEN
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = DcornerCoordinatesCoarse(:,IverticesAtElementCoarse(:,IELA(i)))
          daspectRatio(i) = gaux_getAspectRatio_quad2D (dcoords)
          IF (daspectRatio(i) .LT. 1.0_DP) daspectRatio(i) = 1.0_DP/daspectRatio(i)
          
          ! and the area of that element.
          darea(i) = DelementAreaCoarse(IELA(i))
        ELSE
          daspectRatio(i) = 0.0_DP
          darea(i) = 0.0_DP
        END IF
      END DO
      
      ! Calculate weighting factors for the interpolation.
      ! The iweightingType parameter describes:
      ! <= 0: simple interpolation, weight both contribtutions by 1/2
      !  = 1: take the weighted mean of the interpolated function values
      !       by weighting with the area of the current coarse grid element
      ! >= 2: take the weighted mean of the interpolated function values
      !       by weighting with the area of the neighboured coarse grid element

      SELECT CASE (iweightingType)
      CASE (:0) 
        dweight = 0.5_DP
      CASE (1)
        dweight = darea(0) / (darea(0)+darea(1:TRIA_MAXNME2D))
      CASE (2:)
        dweight = darea(1:TRIA_MAXNME2D) / (darea(0)+darea(1:TRIA_MAXNME2D))
      END SELECT
      
      ! Where there is no neighbour, set the weighting factor to 1.0
!      DO i=1,TRIA_MAXNMT
!        IF (IELA(i) .EQ. 0) dweight(i) = 1.0_DP
!      END DO
      WHERE(IELA(1:TRIA_MAXNME2D) .EQ. 0) dweight(1:TRIA_MAXNME2D) = 1.0_DP

      ! Now determine on which edge to switch to constant prolongation
      ! By default, we don't use constant prolongation
      idoConstant = 1
      
      ! ... but if the caller wants us to switch in a special situation...
      IF ((iarIndicator .GE. 1) .AND. (daspectRatioBound .GE. 0.0_DP)) THEN
        
        ! ... switch to constant of our element is too large...
        IF (darea(0) .GT. daspectRatioBound) idoConstant = 2
        
        ! and if iarIndicator>2, also check the neighbour element
        IF (iarIndicator .GE. 2) THEN
!          DO i=1,TRIA_MAXNME2D
!            IF (darea(i) .GT. daspectRatioBound) idoConstant(i) = 2
!          END DO
          WHERE (darea(1:4) .GT. daspectRatioBound) idoConstant = 2
        END IF
      
      END IF
      
      ! Get the DOF's of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)-NVTcoarse
      IM2 = IedgesAtElementCoarse(2,iel)-NVTcoarse
      IM3 = IedgesAtElementCoarse(3,iel)-NVTcoarse
      IM4 = IedgesAtElementCoarse(4,iel)-NVTcoarse

      ! Get the values of the corresponding DOF's
      DUH1 = DuCoarse(IM1)
      DUH2 = DuCoarse(IM2)
      DUH3 = DuCoarse(IM3)
      DUH4 = DuCoarse(IM4)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Now let's start with the actual prolongation
      ! ---------------------------------------------
      
      ! Get the DOF's on the fine grid
      IA=IedgesAtElementFine(1,IELH1)-NVTfine
      IB=IedgesAtElementFine(4,IELH2)-NVTfine
      IC=IedgesAtElementFine(2,IELH1)-NVTfine

      ! Now we have the following situation:
      
      !   4               IM3                3
      !     ===============X================
      !     |              |               |
      !     |              |               |
      !     |    IELH4     |     IELH3     |
      !     |              |               |
      !     |                              |
      ! IM4 X----------- IEL1 -------------X IM2
      !     |                              |
      !     |              |               |
      !     |    IELH1     o IC  IELH2     |
      !     |              |               |
      !     |              |               |
      !   1 =======o=======X=======o======== 2
      !     |     IA      IM1      IB      |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |            IELA1             |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     ================================

      ! Distribute the value at the edge IM1 to the 
      ! corresponding fine inner nodes

      DuFine(IA) = DuFine(IA) &
                  + dweight(1)*(prweight(1,idoConstant(1))*DUH1 &
                              +prweight(2,idoConstant(1))*DUH2 &
                              +prweight(3,idoConstant(1))*DUH3 &
                              +prweight(4,idoConstant(1))*DUH4)
      DuFine(IB) = DuFine(IB) &
                  + dweight(1)*(prweight(1,idoConstant(1))*DUH1 &
                              +prweight(4,idoConstant(1))*DUH2 &
                              +prweight(3,idoConstant(1))*DUH3 &
                              +prweight(2,idoConstant(1))*DUH4)
      DuFine(IC) = prweight(5,idoConstant(1))*DUH1 &
                 + prweight(6,idoConstant(1))*(DUH2+DUH4) &
                 + prweight(7,idoConstant(1))*DUH3

      ! Distribute the value at the edge IM2 to the 
      ! corresponding fine inner nodes
      IA=IedgesAtElementFine(1,IELH2)-NVTfine
      IB=IedgesAtElementFine(4,IELH3)-NVTfine
      IC=IedgesAtElementFine(2,IELH2)-NVTfine

      DuFine(IA) = DuFine(IA) &
                 + dweight(2)*(prweight(1,idoConstant(2))*DUH2 &
                              +prweight(2,idoConstant(2))*DUH3 &
                              +prweight(3,idoConstant(2))*DUH4 &
                              +prweight(4,idoConstant(2))*DUH1)
      DuFine(IB) = DuFine(IB) &
                 + dweight(2)*(prweight(1,idoConstant(2))*DUH2 &
                              +prweight(4,idoConstant(2))*DUH3 &
                              +prweight(3,idoConstant(2))*DUH4 &
                              +prweight(2,idoConstant(2))*DUH1)
      DuFine(IC) = prweight(5,idoConstant(2))*DUH2 &
                 + prweight(6,idoConstant(2))*(DUH3+DUH1) &
                 + prweight(7,idoConstant(2))*DUH4

      ! Distribute the value at the edge IM3 to the 
      ! corresponding fine inner nodes
      IA=IedgesAtElementFine(1,IELH3)-NVTfine
      IB=IedgesAtElementFine(4,IELH4)-NVTfine
      IC=IedgesAtElementFine(2,IELH3)-NVTfine

      DuFine(IA) = DuFine(IA) &
                 + dweight(3)*(prweight(1,idoConstant(3))*DUH3 &
                              +prweight(2,idoConstant(3))*DUH4 &
                              +prweight(3,idoConstant(3))*DUH1 &
                              +prweight(4,idoConstant(3))*DUH2)
      DuFine(IB) = DuFine(IB) &
                 + dweight(3)*(prweight(1,idoConstant(3))*DUH3 &
                              +prweight(4,idoConstant(3))*DUH4 &
                              +prweight(3,idoConstant(3))*DUH1 &
                              +prweight(2,idoConstant(3))*DUH2)
      DuFine(IC) = prweight(5,idoConstant(3))*DUH3 &
                 + prweight(6,idoConstant(3))*(DUH4+DUH2) &
                 + prweight(7,idoConstant(3))*DUH1

      ! Distribute the value at the edge IM4 to the 
      ! corresponding fine inner nodes
      IA=IedgesAtElementFine(1,IELH4)-NVTfine
      IB=IedgesAtElementFine(4,IELH1)-NVTfine
      IC=IedgesAtElementFine(2,IELH4)-NVTfine

      DuFine(IA) = DuFine(IA) &
                 + dweight(4)*(prweight(1,idoConstant(4))*DUH4 &
                              +prweight(2,idoConstant(4))*DUH1 &
                              +prweight(3,idoConstant(4))*DUH2 &
                              +prweight(4,idoConstant(4))*DUH3)
      DuFine(IB) = DuFine(IB) &
                 + dweight(4)*(prweight(1,idoConstant(4))*DUH4 &
                              +prweight(4,idoConstant(4))*DUH1 &
                              +prweight(3,idoConstant(4))*DUH2 &
                              +prweight(2,idoConstant(4))*DUH3)
      DuFine(IC) = prweight(5,idoConstant(4))*DUH4 &
                 + prweight(6,idoConstant(4))*(DUH1+DUH3) &
                 + prweight(7,idoConstant(4))*DUH2

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformEx30_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! E030/EM30 element, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IedgesAtElement array on the coarse grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine

  ! Number of vertices in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse

  ! Number of vertices in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>

    ! local variables
    REAL(DP) :: DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8,DUH9,DUH10,DUH11,DUH12
    INTEGER(PREC_EDGEIDX) :: IM1,IM2, IM3,IM4
    INTEGER(PREC_EDGEIDX) :: I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12
    INTEGER(PREC_ELEMENTIDX) :: iel, IELH1, IELH2, IELH3, IELH4
    
    ! Weights for the restruction
    REAL(DP), PARAMETER :: A1=1.0_DP, A2=-0.125D0, A3=0.0_DP, A4=0.125D0
    REAL(DP), PARAMETER :: A5=0.625_DP, A6=0.125_DP, A7=0.125_DP, A8=0.125_DP

    ! Loop over the coarse grid elements
    DO iel=1,NELcoarse

      ! Get the DOF's of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)-NVTcoarse
      IM2 = IedgesAtElementCoarse(2,iel)-NVTcoarse
      IM3 = IedgesAtElementCoarse(3,iel)-NVTcoarse
      IM4 = IedgesAtElementCoarse(4,iel)-NVTcoarse

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Get the DOF's on the fine grid
      I1=IedgesAtElementFine(1,IELH1)-NVTfine
      I2=IedgesAtElementFine(4,IELH2)-NVTfine
      I3=IedgesAtElementFine(1,IELH2)-NVTfine
      I4=IedgesAtElementFine(4,IELH3)-NVTfine
      I5=IedgesAtElementFine(1,IELH3)-NVTfine
      I6=IedgesAtElementFine(4,IELH4)-NVTfine
      I7=IedgesAtElementFine(1,IELH4)-NVTfine
      I8=IedgesAtElementFine(4,IELH1)-NVTfine
      I9=IedgesAtElementFine(2,IELH1)-NVTfine
      I10=IedgesAtElementFine(2,IELH2)-NVTfine
      I11=IedgesAtElementFine(2,IELH3)-NVTfine
      I12=IedgesAtElementFine(2,IELH4)-NVTfine

      ! Get the values of the DOF's on the fine grid
      DUH1= DuFine(I1)
      DUH2= DuFine(I2)
      DUH3= DuFine(I3)
      DUH4= DuFine(I4)
      DUH5= DuFine(I5)
      DUH6= DuFine(I6)
      DUH7= DuFine(I7)
      DUH8= DuFine(I8)
      DUH9= DuFine(I9)
      DUH10=DuFine(I10)
      DUH11=DuFine(I11)
      DUH12=DuFine(I12)
      
      ! Now collect the information in the same way as it was
      ! 'distributed' by the prolongation routine.
      ! This realises the adjoint operator of the prolongation.
      !
      ! Calculate the value of the edge IM1

      IF (IneighboursAtElementCoarse(1,iel).NE.0) THEN
        ! inner edge
        IF (IneighboursAtElementCoarse(1,iel).GT.iel) THEN
          DuCoarse(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7) &
                       + A3*(DUH5+DUH6)+A4*(DUH3+DUH8) &
                       + A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ELSE
          DuCoarse(IM1)= DuCoarse(IM1)+A2*(DUH4+DUH7) &
                       + A3*(DUH5+DUH6)+A4*(DUH3+DUH8) &
                       + A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ENDIF
      ELSE
        ! boundary edge 
        DuCoarse(IM1)=     A1*(DUH1+DUH2)+2D0*A2*(DUH4+DUH7) &
                +2D0*A3*(DUH5+DUH6)+2D0*A4*(DUH3+DUH8) &
                +    A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
      ENDIF
 
      ! Calculate the value of the edge IM2
 
      IF (IneighboursAtElementCoarse(2,iel).NE.0) THEN 
        ! inner edge
        IF (IneighboursAtElementCoarse(2,iel).GT.iel) THEN
           DuCoarse(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1) &
                         +A3*(DUH7+DUH8)+A4*(DUH5+DUH2) &
                         +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ELSE
          DuCoarse(IM2)=DuCoarse(IM2)+A2*(DUH6+DUH1) &
                              +A3*(DUH7+DUH8)+A4*(DUH5+DUH2) &
                              +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ENDIF
      ELSE
        ! boundary edge
        DuCoarse(IM2)= A1*(DUH3+DUH4)+2D0*A2*(DUH6+DUH1) &
                +2D0*A3*(DUH7+DUH8)+2D0*A4*(DUH5+DUH2) &
                +    A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
      ENDIF
      
      ! Calculate the value of the edge IM3
      
      IF (IneighboursAtElementCoarse(3,iel).NE.0) THEN 
        ! inner edge
        IF (IneighboursAtElementCoarse(3,iel).GT.iel) THEN
           DuCoarse(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3) &
                        + A3*(DUH1+DUH2)+A4*(DUH7+DUH4) &
                        + A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ELSE
           DuCoarse(IM3)= DuCoarse(IM3)+A2*(DUH8+DUH3) &
                        + A3*(DUH1+DUH2)+A4*(DUH7+DUH4) &
                        + A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ENDIF
      ELSE
        ! boundary edge
        DuCoarse(IM3)= A1*(DUH5+DUH6)+2D0*A2*(DUH8+DUH3) &
                +2D0*A3*(DUH1+DUH2)+2D0*A4*(DUH7+DUH4) &
                +    A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
      ENDIF

      ! Calculate the value of the edge IM4
      
      IF (IneighboursAtElementCoarse(4,iel).NE.0) THEN 
        ! inner edge
        IF (IneighboursAtElementCoarse(4,iel).GT.iel) THEN
          DuCoarse(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5) &
                       + A3*(DUH3+DUH4)+A4*(DUH1+DUH6) &
                       + A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ELSE
          DuCoarse(IM4)= DuCoarse(IM4)+A2*(DUH2+DUH5) &
                       + A3*(DUH3+DUH4)+A4*(DUH1+DUH6) &
                       + A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ENDIF
      ELSE
        ! boundary edge
        DuCoarse(IM4)= A1*(DUH7+DUH8)+2D0*A2*(DUH2+DUH5) &
                +2D0*A3*(DUH3+DUH4)+2D0*A4*(DUH1+DUH6) &
                +    A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
      ENDIF

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformEx30ext_double (DuCoarse,DuFine, &
               DcornerCoordinatesCoarse,IverticesAtElementCoarse,DelementAreaCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NELcoarse, &
               iweightingType, daspectRatioBound, iarIndicator)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! E030/EM30 element, uniform triangulation, double precision vector.
  !
  ! Extended version. Switch to constant restriction if aspect ratio
  ! of an element is too large.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! DcornerCoordinates array on the coarse grid
  REAL(DP), DIMENSION(:,:), INTENT(IN)                :: DcornerCoordinatesCoarse

  ! IverticesAtElement array on the coarse grid
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse
  
  ! DelementArea array on the coarse grid
  REAL(DP), DIMENSION(:), INTENT(IN)                  :: DelementAreaCoarse

  ! IedgesAtElement array on the coarse grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine

  ! Number of vertices in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse

  ! Number of vertices in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse

  ! Type of the averaging on the element edges
  ! <=0: standard averaging of both contributions by 1/2,        
  !  =1: weighted averaging of the interpolated function values:
  !      The area of the current coarse grid element determines 
  !      the weight. (L2-projection, standard),
  !  =2: weighted averaging of the interpolated function values:
  !      The area of the neightbour element of the coarse grid 
  !      the weight. 
  INTEGER, INTENT(IN)  :: iweightingType
  
  ! Upper bound aspect ratio; for all elements with higher AR
  ! the prolongation is switched to constant prolongation 
  REAL(DP), INTENT(IN) :: daspectRatioBound
  
  ! Aspect-ratio indicator.
  ! Controls switching to constant prolongation.
  ! <=1: switch depending on aspect ratio of current element,
  !  =2: switch depending on aspect ratio of current element and
  !      neighbour element
  INTEGER, INTENT(IN)  :: iarIndicator
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>

    ! local variables
    REAL(DP) :: DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8,DUH9,DUH10,DUH11,DUH12
    INTEGER(PREC_EDGEIDX) :: IM1,IM2, IM3,IM4
    INTEGER(PREC_EDGEIDX) :: I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12
    INTEGER(PREC_ELEMENTIDX) :: iel, IELH1, IELH2, IELH3, IELH4
    INTEGER(PREC_ELEMENTIDX), DIMENSION(0:TRIA_MAXNME2D) :: IELA
    INTEGER :: i
    INTEGER, DIMENSION(TRIA_MAXNME2D) :: idoConstant
    REAL(DP), DIMENSION(0:TRIA_MAXNME2D) :: daspectRatio, darea
    REAL(DP), DIMENSION(TRIA_MAXNME2D) :: dweight
    REAL(DP), DIMENSION(NDIM2D,TRIA_MAXNVE2D) :: dcoords
    
    ! Weights for the restriction.
    ! PRWEIG (.,1) gives the constants for the standard restriction,
    ! PRWEIG (.,2) gives the constants for the constant restriction.
    REAL(DP), DIMENSION(8,2), PARAMETER :: prweight = &
        RESHAPE((/1.0_DP, -0.25_DP, 0.0_DP, 0.25_DP, &
                  0.625_DP, 0.125_DP, 0.125_DP, 0.125_DP, &
                  1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                  1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP/),(/8,2/))

    ! Clear the output vector
    CALL lalg_clearVectorDble(DuCoarse)
              
    ! Loop over the coarse grid elements
    DO iel=1,NELcoarse

      ! Get the numbers of the elements that are neighbours to our current coarse
      ! grid element:
      !             +--------+
      !             |        |
      !             | IELA3  |
      !             |        |
      !    +--------4--------3--------+
      !    |        |        |        |
      !    | IELA4  |  IEL   | IELA2  |
      !    |        |        |        |
      !    +--------1--------2--------+
      !             |        |
      !             | IELA1  |
      !             |        |
      !             +--------+
      IELA(0) = iel
      IELA(1) = IneighboursAtElementCoarse(1,iel)
      IELA(2) = IneighboursAtElementCoarse(2,iel)
      IELA(3) = IneighboursAtElementCoarse(3,iel)
      IELA(4) = IneighboursAtElementCoarse(4,iel)
      
      ! For these five elements, determine the aspect ratio and their area.
      !
      ! At first the element in the center, which always exists.
      !
      ! Get the aspect ratio of the current coarse grid element;
      ! if necessary, calculate the reciprocal.
      dcoords = DcornerCoordinatesCoarse(:,IverticesAtElementCoarse(:,IELA(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      IF (daspectRatio(0) .LT. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)
      
      ! and the area of that element.
      darea(0) = DelementAreaCoarse(iel)
      
      ! Then the remaining neighbours.
      DO i=1,TRIA_MAXNME2D
        IF (IELA(i) .NE. 0) THEN
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = DcornerCoordinatesCoarse(:,IverticesAtElementCoarse(:,IELA(i)))
          daspectRatio(i) = gaux_getAspectRatio_quad2D (dcoords)
          IF (daspectRatio(i) .LT. 1.0_DP) daspectRatio(i) = 1.0_DP/daspectRatio(i)
          
          ! and the area of that element.
          darea(i) = DelementAreaCoarse(IELA(i))
        ELSE
          daspectRatio(i) = 0.0_DP
          darea(i) = 0.0_DP
        END IF
      END DO
      
      ! Calculate weighting factors for the interpolation.
      ! The iweightingType parameter describes:
      ! <= 0: simple interpolation, weight both contribtutions by 1/2
      !  = 1: take the weighted mean of the interpolated function values
      !       by weighting with the area of the current coarse grid element
      ! >= 2: take the weighted mean of the interpolated function values
      !       by weighting with the area of the neighboured coarse grid element

      SELECT CASE (iweightingType)
      CASE (:0) 
        dweight = 0.5_DP
      CASE (1)
        dweight = darea(0) / (darea(0)+darea(1:TRIA_MAXNME2D))
      CASE (2:)
        dweight = darea(1:TRIA_MAXNME2D) / (darea(0)+darea(1:TRIA_MAXNME2D))
      END SELECT
      
      ! Where there is no neighbour, set the weighting factor to 1.0
!      DO i=1,TRIA_MAXNMT
!        IF (IELA(i) .EQ. 0) dweight(i) = 1.0_DP
!      END DO
      WHERE(IELA(1:TRIA_MAXNME2D) .EQ. 0) dweight(1:TRIA_MAXNME2D) = 1.0_DP

      ! Now determine on which edge to switch to constant prolongation
      ! By default, we don't use constant prolongation
      idoConstant = 1
      
      ! ... but if the caller wants us to switch in a special situation...
      IF ((iarIndicator .GE. 1) .AND. (daspectRatioBound .GE. 0.0_DP)) THEN
        
        ! ... switch to constant of our element is too large...
        IF (darea(0) .GT. daspectRatioBound) idoConstant = 2
        
        ! and if iarIndicator>2, also check the neighbour element
        IF (iarIndicator .GE. 2) THEN
!          DO i=1,TRIA_MAXNME2D
!            IF (darea(i) .GT. daspectRatioBound) idoConstant(i) = 2
!          END DO
          WHERE (darea(1:4) .GT. daspectRatioBound) idoConstant = 2
        END IF
      
      END IF

      ! Now let's strt with the actual restriction
      ! ------------------------------------------      

      ! Get the DOF's of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)-NVTcoarse
      IM2 = IedgesAtElementCoarse(2,iel)-NVTcoarse
      IM3 = IedgesAtElementCoarse(3,iel)-NVTcoarse
      IM4 = IedgesAtElementCoarse(4,iel)-NVTcoarse

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Get the DOF's on the fine grid
      I1=IedgesAtElementFine(1,IELH1)-NVTfine
      I2=IedgesAtElementFine(4,IELH2)-NVTfine
      I3=IedgesAtElementFine(1,IELH2)-NVTfine
      I4=IedgesAtElementFine(4,IELH3)-NVTfine
      I5=IedgesAtElementFine(1,IELH3)-NVTfine
      I6=IedgesAtElementFine(4,IELH4)-NVTfine
      I7=IedgesAtElementFine(1,IELH4)-NVTfine
      I8=IedgesAtElementFine(4,IELH1)-NVTfine
      I9=IedgesAtElementFine(2,IELH1)-NVTfine
      I10=IedgesAtElementFine(2,IELH2)-NVTfine
      I11=IedgesAtElementFine(2,IELH3)-NVTfine
      I12=IedgesAtElementFine(2,IELH4)-NVTfine

      ! Get the values of the DOF's on the fine grid
      DUH1= DuFine(I1)
      DUH2= DuFine(I2)
      DUH3= DuFine(I3)
      DUH4= DuFine(I4)
      DUH5= DuFine(I5)
      DUH6= DuFine(I6)
      DUH7= DuFine(I7)
      DUH8= DuFine(I8)
      DUH9= DuFine(I9)
      DUH10=DuFine(I10)
      DUH11=DuFine(I11)
      DUH12=DuFine(I12)
      
      ! Now we have the following situation:
      !   4       I6      IM3      I5      3
      !     =======o=======X=======o========
      !     |              |               |
      !     |              |               |
      !  I7 o    IELH4     o I11 IELH3     o I4
      !     |              |               |
      !     |                              |
      ! IM4 X------o---- IEL1 -----o-------X IM2
      !     |    I12               I10     |
      !     |              |               |
      !  I8 o    IELH1     o I9  IELH2     o I3
      !     |              |               |
      !     |              |               |
      !   1 =======o=======X=======o======== 2
      !     |     I1      IM1      I2      |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |            IELA1             |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     ================================
      
      
      ! Collect the information in the same way as it was
      ! 'distributed' by the prolongation routine.
      ! This realises the adjoint operator of the prolongation.
      !
      ! Calculate the value of the edge IM1

      DuCoarse(IM1)= DuCoarse(IM1) &
                    +dweight(1)*(prweight(1,idoConstant(1))*(DUH1+DUH2) &
                                +prweight(2,idoConstant(1))*(DUH4+DUH7) &
                                +prweight(3,idoConstant(1))*(DUH5+DUH6) &
                                +prweight(4,idoConstant(1))*(DUH3+DUH8)) &
                    +prweight(5,idoConstant(1))*DUH9 &
                    +prweight(6,idoConstant(1))*(DUH10+DUH12) &
                    +prweight(7,idoConstant(1))*DUH11
 
      ! Calculate the value of the edge IM2
      
      DuCoarse(IM2)= DuCoarse(IM2) &
                    +dweight(2)*(prweight(1,idoConstant(2))*(DUH3+DUH4) &
                                +prweight(2,idoConstant(2))*(DUH6+DUH1) &
                                +prweight(3,idoConstant(2))*(DUH7+DUH8) &
                                +prweight(4,idoConstant(2))*(DUH5+DUH2)) &
                    +prweight(5,idoConstant(2))*DUH10 &
                    +prweight(6,idoConstant(2))*(DUH11+DUH9) &
                    +prweight(7,idoConstant(2))*DUH12
      
      ! Calculate the value of the edge IM3
      
      DuCoarse(IM3)= DuCoarse(IM3) &
                    +dweight(3)*(prweight(1,idoConstant(3))*(DUH5+DUH6) &
                                +prweight(2,idoConstant(3))*(DUH8+DUH3) &
                                +prweight(3,idoConstant(3))*(DUH1+DUH2) &
                                +prweight(4,idoConstant(3))*(DUH7+DUH4)) &
                    +prweight(5,idoConstant(3))*DUH11 &
                    +prweight(6,idoConstant(3))*(DUH12+DUH10) &
                    +prweight(7,idoConstant(3))*DUH9

      ! Calculate the value of the edge IM4
      
      DuCoarse(IM4)= DuCoarse(IM4) &
                    +dweight(4)*(prweight(1,idoConstant(4))*(DUH7+DUH8) &
                                +prweight(2,idoConstant(4))*(DUH2+DUH5) &
                                +prweight(3,idoConstant(4))*(DUH3+DUH4) &
                                +prweight(4,idoConstant(4))*(DUH1+DUH6)) &
                    +prweight(5,idoConstant(4))*DUH12 &
                    +prweight(6,idoConstant(4))*(DUH9+DUH11) &
                    +prweight(7,idoConstant(4))*DUH10

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformEx30_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NELcoarse)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! Ex30, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IedgesAtElement array on the coarse grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine

  ! Number of vertices in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse

  ! Number of vertices in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>

    ! local variables
    REAL(DP) :: DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8,DUH9,DUH10,DUH11,DUH12
    INTEGER(PREC_EDGEIDX) :: IM1,IM2, IM3,IM4
    INTEGER(PREC_EDGEIDX) :: I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12
    INTEGER(PREC_ELEMENTIDX) :: iel, IELH1, IELH2, IELH3, IELH4
    
    ! Weights for the restruction
    REAL(DP), PARAMETER :: A1=0.1875_DP, A2=0.375_DP, A3=-0.0625_DP
    REAL(DP), PARAMETER :: R1=0.375_DP, R2=0.75_DP, R3=-0.125_DP

    ! Loop over the coarse grid elements
    DO iel=1,NELcoarse

      ! Get the DOF's of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)-NVTcoarse
      IM2 = IedgesAtElementCoarse(2,iel)-NVTcoarse
      IM3 = IedgesAtElementCoarse(3,iel)-NVTcoarse
      IM4 = IedgesAtElementCoarse(4,iel)-NVTcoarse

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Get the DOF's on the fine grid
      I1=IedgesAtElementFine(1,IELH1)-NVTfine
      I2=IedgesAtElementFine(4,IELH2)-NVTfine
      I3=IedgesAtElementFine(1,IELH2)-NVTfine
      I4=IedgesAtElementFine(4,IELH3)-NVTfine
      I5=IedgesAtElementFine(1,IELH3)-NVTfine
      I6=IedgesAtElementFine(4,IELH4)-NVTfine
      I7=IedgesAtElementFine(1,IELH4)-NVTfine
      I8=IedgesAtElementFine(4,IELH1)-NVTfine
      I9=IedgesAtElementFine(2,IELH1)-NVTfine
      I10=IedgesAtElementFine(2,IELH2)-NVTfine
      I11=IedgesAtElementFine(2,IELH3)-NVTfine
      I12=IedgesAtElementFine(2,IELH4)-NVTfine

      ! Get the values of the DOF's on the fine grid
      DUH1= DuFine(I1)
      DUH2= DuFine(I2)
      DUH3= DuFine(I3)
      DUH4= DuFine(I4)
      DUH5= DuFine(I5)
      DUH6= DuFine(I6)
      DUH7= DuFine(I7)
      DUH8= DuFine(I8)
      DUH9= DuFine(I9)
      DUH10=DuFine(I10)
      DUH11=DuFine(I11)
      DUH12=DuFine(I12)
      
      ! Now interpolate the fine-grid values to the coarse grid nodes.
      !
      ! Calculate the value of the edge IM1

      IF (IneighboursAtElementCoarse(1,iel).NE.0) THEN
        ! inner edge
        IF (IneighboursAtElementCoarse(1,iel).GT.iel) THEN
          DuCoarse(IM1)= A1*(DUH1+DUH2) +A2*DUH9 +A3*(DUH8+DUH3+DUH10+DUH12)
        ELSE
          DuCoarse(IM1)= DuCoarse(IM1) &
                       + A1*(DUH1+DUH2) +A2*DUH9 +A3*(DUH8+DUH3+DUH10+DUH12)
        ENDIF
      ELSE
        ! boundary edge 
        DuCoarse(IM1)= R1*(DUH1+DUH2) +R2*DUH9 +R3*(DUH8+DUH3+DUH10+DUH12)
      ENDIF
 
      ! Calculate the value of the edge IM2
 
      IF (IneighboursAtElementCoarse(2,iel).NE.0) THEN 
        ! inner edge
        IF (IneighboursAtElementCoarse(2,iel).GT.iel) THEN
           DuCoarse(IM2)= A1*(DUH3+DUH4) +A2*DUH10 +A3*(DUH2+DUH5+DUH9 +DUH11)
        ELSE
          DuCoarse(IM2)= DuCoarse(IM2) &
                       + A1*(DUH3+DUH4) +A2*DUH10 +A3*(DUH2+DUH5+DUH9 +DUH11)
        ENDIF
      ELSE
        ! boundary edge
        DuCoarse(IM2)= R1*(DUH3+DUH4) +R2*DUH10 +R3*(DUH2+DUH5+DUH9 +DUH11)
      ENDIF
      
      ! Calculate the value of the edge IM3
      
      IF (IneighboursAtElementCoarse(3,iel).NE.0) THEN 
        ! inner edge
        IF (IneighboursAtElementCoarse(3,iel).GT.iel) THEN
           DuCoarse(IM3)= A1*(DUH5+DUH6) +A2*DUH11 +A3*(DUH4+DUH7+DUH10+DUH12)
        ELSE
           DuCoarse(IM3)= DuCoarse(IM3) &
                        + A1*(DUH5+DUH6) +A2*DUH11 +A3*(DUH4+DUH7+DUH10+DUH12)
        ENDIF
      ELSE
        ! boundary edge
        DuCoarse(IM3)= R1*(DUH5+DUH6) +R2*DUH11 +R3*(DUH4+DUH7+DUH10+DUH12)
      ENDIF

      ! Calculate the value of the edge IM4
      
      IF (IneighboursAtElementCoarse(4,iel).NE.0) THEN 
        ! inner edge
        IF (IneighboursAtElementCoarse(4,iel).GT.iel) THEN
          DuCoarse(IM4)= A1*(DUH7+DUH8) +A2*DUH12 +A3*(DUH6+DUH1+DUH9 +DUH11)
        ELSE
          DuCoarse(IM4)= DuCoarse(IM4) &
                       + A1*(DUH7+DUH8) +A2*DUH12 +A3*(DUH6+DUH1+DUH9 +DUH11)
        ENDIF
      ELSE
        ! boundary edge
        DuCoarse(IM4)= R1*(DUH7+DUH8) +R2*DUH12 +R3*(DUH6+DUH1+DUH9 +DUH11)
      ENDIF

    END DO

  END SUBROUTINE

END MODULE
