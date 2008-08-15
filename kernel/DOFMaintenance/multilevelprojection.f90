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
!#     mlprj_initProjectionDiscr /
!#     mlprj_initProjectionVec /
!#     mlprj_initProjectionMat
!#     -> Initialises a projection structure with values according to a
!#        spatial discretisation.
!#        Uses a discretisation structure, a (block) vector, a (block)
!#        matrix on the fine/coarse grid or a block discretisation structure
!#        as template.
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
!# 7.) mlprj_setL2ProjMatrices
!#     -> Sets the matrices for an scalar projection structure which are
!#        needed for L2-projection.
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

!<constants>

!<constantblock description="Projection Types">
  ! Hard-coded projection operators
  INTEGER(I32), PARAMETER :: MLP_PROJ_TYPE_HARDCODED = 0
  
  ! L2-projection operators
  INTEGER(I32), PARAMETER :: MLP_PROJ_TYPE_L2 = 1

!</constantblock>

!</constants>

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
  
    ! Specifies the projection type of this projection structure.
    ! One of the MLP_PROJ_TYPE_XXXX constants defined above.
    INTEGER(I32)                :: iprojType = MLP_PROJ_TYPE_HARDCODED
  
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
    INTEGER(I32)                :: ielementTypeProlongation = EL_UNDEFINED
    
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
    ! Allows to switch to a special-type prolongation. Whether or not this
    ! has an effect depends on the discretisation.
    ! = 0: Use default prolongation.
    ! Uniform discretisation with E030/E031/EM30/EM31:
    ! = 1: Use standard prolongation, equally weighted (1/2 from left, 1/2 from right),
    ! = 2: Use extended prolongation, equally weighted (1/2 from left, 1/2 from right),
    ! = 3: Use extended prolongation, weighted by element size (L2 projection),
    ! = 4: Use extended prolongation, weighted by element size of neighbour element
    ! To activate extended prolongation, set this to >= 2 after initialising the
    ! interlevel projection structure!
    ! Uniform discretisation with Q1:
    ! = 2: FEAST mirror boundary prolongation with zero on the boundary (experimentally)
    INTEGER                     :: iprolVariant = 0
    
    ! Configuration parameter for extended prolongation of E030/EM30/E031/EM31
    ! element. Only valid if iprolVariant >= 2.
    ! Aspect-ratio indicator; controls switching to constant prolongation.
    ! <=1: switch depending on aspect ratio of current element (standard),
    !  =2: switch depending on aspect ratio of current element and
    !      neighbour element
    INTEGER                     :: iprolARIndicatorEX3Y = 1
    
    ! Configuration parameter for extended prolongation of E030/EM30/E031/EM31
    ! element. Only valid if iprolVariant >= 2.
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
    INTEGER(I32)                :: ielementTypeRestriction = EL_UNDEFINED

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
    ! Allows to switch to a special-type restriction. Whether or not this
    ! has an effect depends on the discretisation.
    ! = 0: Use default restriction.
    ! Uniform discretisation with E030/E031/EM30/EM31:
    ! = 1: Use standard restriction, equally weighted (1/2 from left, 1/2 from right),
    ! = 2: Use extended restriction, equally weighted (1/2 from left, 1/2 from right),
    ! = 3: Use extended restriction, weighted by element size (L2 projection),
    ! = 4: Use extended restriction, weighted by element size of neighbour element
    ! To activate extended prolongation, set this to >= 2 after initialising the
    ! interlevel projection structure!
    ! Uniform discretisation with Q1:
    ! = 1: Element-wise FEAST mirror boundary restriction (experimentally)
    ! = 2: FEAST mirror boundary restriction with zero on the boundary (experimentally)
    INTEGER                     :: irestVariant = 0
    
    ! Configuration parameter for extended restriction of E030/EM30/E031/EM31
    ! element. Only valid if irestVariant >= 2.
    ! Aspect-ratio indicator; controls switching to constant prolongation.
    ! <=1: switch depending on aspect ratio of current element (standard),
    !  =2: switch depending on aspect ratio of current element and
    !      neighbour element
    INTEGER                     :: irestARIndicatorEX3Y = 1
    
    ! Configuration parameter for extended restriction of E030/EM30/E031/EM31
    ! element. Only valid if irestVariant >= 2.
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
    INTEGER(I32)                :: ielementTypeInterpolation = EL_UNDEFINED
    
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
    
    ! -------------------------------------------------------------------------
    ! L2-Projection Structures
    ! -------------------------------------------------------------------------
    ! Mass matrix of the fine mesh spatial discretisation.
    TYPE(t_matrixScalar)        :: rmatrixMass
    
    ! Lumped Mass matrix of the fine mesh spatial discretisation.
    TYPE(t_matrixScalar)        :: rlumpedMass
    
    ! 2-Level-Mass matrix and its virtually transpose
    TYPE(t_matrixScalar)        :: rmatrix2LvlMass
    TYPE(t_matrixScalar)        :: rmatrix2LvlMassT
    
    ! Two temporary vectors
    TYPE(t_vectorScalar)        :: rvectorTmp
    TYPE(t_vectorScalar)        :: rvectorDef
    
    ! Number of iterations
    INTEGER                     :: imaxL2Iterations = 30
    
    ! Relative and absolute tolerance
    REAL(DP)                    :: depsRelL2 = 1E-5_DP
    REAL(DP)                    :: depsAbsL2 = 1E-10_DP
    
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
    ! DIMENSION(1..#FE-spaces, 1..#equations)
    TYPE(t_interlevelProjectionScalar), DIMENSION(:,:), POINTER :: RscalarProjection => NULL()
  
  END TYPE
  
  !</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_initProjectionDirect (rprojection,RspatialDiscr)
  
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
  TYPE(t_spatialDiscretisation), DIMENSION(:), INTENT(IN) :: RspatialDiscr
!</input>
  
!<output>
  ! A t_interlevelProjectionBlock structure that will be filled with data
  ! about the projection of all the equations described by RspatialDiscr.
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
!</output>
  
!</subroutine>

    INTEGER :: nFEspaces,nequations,i
    
    ! Get the max. number of FE-spaces and the number of equations
    nequations = SIZE(RspatialDiscr)
    nFEspaces = 0
    DO i=1,nequations
      nFEspaces = MAX(nFEspaces,RspatialDiscr(i)%inumFESpaces)
    END DO

    ! Allocate spatial discretisation structures for all the subblocks
    ALLOCATE(rprojection%RscalarProjection(nFEspaces,nequations))
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_initProjectionDiscr (rprojection,rdiscretisation)
  
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
  ! The PDE system is specified by the given block discretisation structure.
!</description>

!<input>
  ! A block discretisation structure specifying the discretisation
  TYPE(t_blockDiscretisation), INTENT(IN) :: rdiscretisation
!</input>
  
!<output>
  ! A t_interlevelProjectionBlock structure that will be filled with data
  ! about the projection of all the equations described by RspatialDiscr.
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
!</output>
  
!</subroutine>

    IF (rdiscretisation%ncomponents .EQ. 0) THEN
      PRINT *,'mlprj_initProjectionDiscr: No discretisation!'
      CALL sys_halt()
    END IF

    ! Call the standard initialisation routine
    CALL mlprj_initProjectionDirect (rprojection,&
         rdiscretisation%RspatialDiscr)

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
  ! about the projection of all the equations described by RspatialDiscr.
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
!</output>
  
!</subroutine>

    ! local variables
    TYPE(t_spatialDiscretisation), DIMENSION(rvector%nblocks) :: Rdiscr
    INTEGER :: i

    IF (rvector%nblocks .EQ. 0) THEN
      PRINT *,'mlprj_initProjectionVec: No discretisation!'
      CALL sys_halt()
    END IF

    ! Set up an array of discretisation structures for all the equations
    DO i=1,rvector%nblocks
      Rdiscr(i) = rvector%RvectorBlock(i)%p_rspatialDiscr
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
  ! about the projection of all the equations described by RspatialDiscr.
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
!</output>
  
!</subroutine>

    ! local variables;
    TYPE(t_spatialDiscretisation), DIMENSION(MAX(rmatrix%ndiagBlocks,1)) :: Rdiscr
    INTEGER :: i,j

    IF (rmatrix%ndiagBlocks .EQ. 0) THEN
      PRINT *,'mlprj_initProjectionMat: No discretisation!'
      CALL sys_halt()
    END IF

    ! Set up an array of discretisation structures for all the equations.
    ! In every 'column' of the block matrix, search for the first existing
    ! matrix and use its properties for initialisation
    DO i=1,rmatrix%ndiagBlocks
      DO j=1,rmatrix%ndiagBlocks
        IF (rmatrix%RmatrixBlock(j,i)%NEQ .NE. 0) THEN
          Rdiscr(i) = &
            rmatrix%RmatrixBlock(j,i)%p_rspatialDiscrTrial
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
  INTEGER :: i,j
  TYPE(t_interlevelProjectionScalar), POINTER :: p_rprj

    ! Release allocated memory
    IF (ASSOCIATED(rprojection%RscalarProjection)) THEN
      
      ! Go through all scalar projections
      DO i = 1, UBOUND(rprojection%RscalarProjection,1)
        DO j = 1, UBOUND(rprojection%RscalarProjection,2)
          
          p_rprj => rprojection%RscalarProjection(i,j)
        
          ! Release all matrices and vectors for L2-projection
          CALL lsyssc_releaseMatrix(p_rprj%rmatrixMass)
          CALL lsyssc_releaseMatrix(p_rprj%rlumpedMass)
          CALL lsyssc_releaseMatrix(p_rprj%rmatrix2LvlMass)
          CALL lsyssc_releaseMatrix(p_rprj%rmatrix2LvlMassT)
          IF(p_rprj%rvectorTmp%NEQ .NE. 0) THEN
            CALL lsyssc_releaseVector(p_rprj%rvectorTmp)
          END IF
          IF(p_rprj%rvectorDef%NEQ .NE. 0) THEN
            CALL lsyssc_releaseVector(p_rprj%rvectorDef)
          END IF
        
        END DO ! j
      END DO ! i
      
      ! Deallocate the scalar projection array
      DEALLOCATE(rprojection%RscalarProjection)
      
    END IF

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
  IF (releDistrCoarse%celement .NE. releDistrFine%celement) THEN
    PRINT *,'Element distribution on the coarse and fine grid incompatible!'
    PRINT *,'Coarse grid: ',releDistrCoarse%celement,&
            ' Fine grid: ',releDistrFine%celement
    CALL sys_halt()
  END IF

  ! Copy the template to the actual projection structure.
  ractProjection = rprojection
  
  ! Is any element type to be replaced according to a discretisation?
  IF (ractProjection%ielementTypeProlongation .EQ. EL_UNDEFINED) THEN
    ractProjection%ielementTypeProlongation = releDistrCoarse%celement
  END IF

  IF (ractProjection%ielementTypeRestriction .EQ. EL_UNDEFINED) THEN
    ractProjection%ielementTypeRestriction = releDistrCoarse%celement
  END IF

  IF (ractProjection%ielementTypeInterpolation .EQ. EL_UNDEFINED) THEN
    ractProjection%ielementTypeInterpolation = releDistrCoarse%celement
  END IF

  ! Check to see if the discretisation structures fit to the projection structure
  IF (ractProjection%ielementTypeProlongation .NE. releDistrCoarse%celement) THEN
    ! Here, we can insert some additional code so that E030 is compatible eith EM30
    ! and so on...
    PRINT *,'Element distribution of the grid and interlevel projection incompatible!'
    PRINT *,'Grid: ',releDistrCoarse%celement,&
            ' Prolongation: ',ractProjection%ielementTypeProlongation
    CALL sys_halt()
  END IF

  IF (ractProjection%ielementTypeRestriction .NE. releDistrCoarse%celement) THEN
    ! Here, we can insert some additional code so that E030 is compatible eith EM30
    ! and so on...
    PRINT *,'Element distribution of the grid and interlevel projection incompatible!'
    PRINT *,'Grid: ',releDistrCoarse%celement,&
            ' Restriction: ',ractProjection%ielementTypeRestriction
    CALL sys_halt()
  END IF

  IF (ractProjection%ielementTypeInterpolation .NE. releDistrCoarse%celement) THEN
    ! Here, we can insert some additional code so that E030 is compatible eith EM30
    ! and so on...
    PRINT *,'Element distribution of the grid and interlevel projection incompatible!'
    PRINT *,'Grid: ',releDistrCoarse%celement,&
            ' Interpolation: ',ractProjection%ielementTypeInterpolation
    CALL sys_halt()
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

!</function>

  ! local variables
  INTEGER :: i
  INTEGER(PREC_VECIDX) :: imemmax,imemact

  IF (SIZE(RdiscrCoarse) .NE. SIZE(RdiscrFine)) THEN
    PRINT *,'mlprj_allocTempVector: Coarse and fine grid incompatible!'
    CALL sys_halt()
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

!</function>

    ! local variables
    TYPE(t_spatialDiscretisation), DIMENSION(rvectorCoarse%nblocks) :: RdiscrCoarse
    TYPE(t_spatialDiscretisation), DIMENSION(rvectorFine%nblocks) :: RdiscrFine
    INTEGER :: i

    IF ((rvectorCoarse%nblocks .EQ. 0) .OR. (rvectorFine%nblocks .EQ. 0)) THEN
      PRINT *,'mlprj_getTempMemoryVec: No discretisation!'
      CALL sys_halt()
    END IF

    ! Set up an array of discretisation structures for all the equations
    DO i=1,rvectorCoarse%nblocks
      RdiscrCoarse(i) = &
        rvectorCoarse%RvectorBlock(i)%p_rspatialDiscr
    END DO

    DO i=1,rvectorFine%nblocks
      RdiscrFine(i) = &
        rvectorFine%RvectorBlock(i)%p_rspatialDiscr
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

!</function>

    ! local variables; 
    TYPE(t_spatialDiscretisation), DIMENSION(MAX(rmatrixCoarse%ndiagBlocks,1)) :: RdiscrCoarse
    TYPE(t_spatialDiscretisation), DIMENSION(MAX(rmatrixFine%ndiagBlocks,1)) :: RdiscrFine
    INTEGER :: i,j

    IF ((rmatrixCoarse%ndiagBlocks .EQ. 0) .OR. (rmatrixFine%ndiagBlocks .EQ. 0)) THEN
      PRINT *,'mlprj_getTempMemoryVec: No discretisation!'
      CALL sys_halt()
    END IF

    ! Set up an array of discretisation structures for all the equations
    DO i=1,rmatrixCoarse%ndiagBlocks
      DO j=1,rmatrixCoarse%ndiagBlocks
        IF (lsysbl_isSubmatrixPresent(rmatrixCoarse,j,i)) THEN
          IF (.NOT. &
              ASSOCIATED(rmatrixCoarse%RmatrixBlock(j,i)%p_rspatialDiscrTrial)) then
            PRINT *,'mlprj_getTempMemoryMat: No discretisation structure in coarse &
                  &matrix at ',i,',',j
            CALL sys_halt()
          END IF
          RdiscrCoarse(i) = &
            rmatrixCoarse%RmatrixBlock(j,i)%p_rspatialDiscrTrial
          EXIT
        END IF
      END DO
    END DO

    DO i=1,rmatrixFine%ndiagBlocks
      DO j=1,rmatrixFine%ndiagBlocks
        IF (lsysbl_isSubmatrixPresent(rmatrixFine,j,i)) THEN
          IF (.NOT. &
              ASSOCIATED(rmatrixFine%RmatrixBlock(j,i)%p_rspatialDiscrTrial)) THEN
            PRINT *,'mlprj_getTempMemoryMat: No discretisation structure in fine matrix&
                  & at ',i,',',j
            CALL sys_halt()
          END IF
          RdiscrFine(i) = &
            rmatrixFine%RmatrixBlock(j,i)%p_rspatialDiscrTrial
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
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementCoarse,&
      p_IverticesAtEdgeCoarse, p_IverticesAtElementFine, p_IverticesAtFaceCoarse
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementCoarse,&
      p_IneighboursAtElementFine
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementCoarse,&
      p_IfacesAtElementCoarse, p_IedgesAtElementFine, p_IfacesAtElementFine
  INTEGER(I32), DIMENSION(:), POINTER :: p_ItwistIndexEdgesCoarse, &
      p_ItwistIndexEdgesFine
  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoordsCoarse
  REAL(DP), DIMENSION(:), POINTER   :: p_DelementAreaCoarse
  
  ! Data arrays
  REAL(DP), DIMENSION(:), POINTER :: p_DuCoarse, p_DuFine
 
    ! The vectors must be of data type DOUBLE - we don't support anything
    ! different at the moment...
    IF (rcoarseVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Coarse grid vector has unsupported data type!'
      CALL sys_halt()
    END IF

    IF (rfineVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Fine grid vector has unsupported data type!'
      CALL sys_halt()
    END IF
    
    IF (lsysbl_isVectorSorted(rfineVector) .OR. lsysbl_isVectorSorted(rcoarseVector)) THEN
      PRINT *,'Vectors must be unsorted for level change!'
      CALL sys_halt()
    END IF

    ! Calls the correct prolongation routine for each block in the 
    ! discretisation...
    DO i=1,rcoarseVector%nblocks
    
      IF (rcoarseVector%RvectorBlock(i)%NEQ .GT. 0) THEN
      
        ! Do we use L2-Projection here?
        IF (rprojection%RscalarProjection(1,i)%iprojType .EQ. &
            MLP_PROJ_TYPE_L2) THEN
          
          ! Call scalar L2-prolongation
          CALL mlprj_prolScalarL2(rprojection%RscalarProjection(1,i), &
            rcoarseVector%RvectorBlock(i), rfineVector%RvectorBlock(i))
        
          ! Continue with next block
          CYCLE
          
        END IF
      
        p_rdiscrCoarse => rcoarseVector%RvectorBlock(i)%p_rspatialDiscr
        p_rdiscrFine => rfineVector%RvectorBlock(i)%p_rspatialDiscr
        
        ! We need a discretisation:
        IF ((.NOT. ASSOCIATED(p_rdiscrCoarse)) .OR. &
            (.NOT. ASSOCIATED(p_rdiscrFine))) THEN
          PRINT *,'Intergrid transfer: No discretisation!'
          CALL sys_halt()
        END IF

        ! Currently, we support only uniform triangulations.
        IF ((p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
            (p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM)) THEN
          PRINT *,'Intergrid transfer supports currently only uniform discretisations!'
          CALL sys_halt()
        END IF
        
        ! Get the pointers to the vectors
        CALL lsyssc_getbase_double (rcoarseVector%RvectorBlock(i),p_DuCoarse)
        CALL lsyssc_getbase_double (rfineVector%RvectorBlock(i),p_DuFine)
        
        ! Use the first projection structure as template and create
        ! the actual projection structure for our situation.
        ! Remember, in a uniform grid we only have one projection structure
        ! and one element distribution!
        CALL mlprj_getProjectionStrategy (rprojection%RscalarProjection(1,i), &
              p_rdiscrCoarse%RelementDistr(1), &
              p_rdiscrFine%RelementDistr(1), &
              ractProjection)
      
        ! Depending on the element type of the trial functions in the
        ! discretisation, choose the right prolongation and call it.
        p_rtriaCoarse => p_rdiscrCoarse%p_rtriangulation
        p_rtriaFine => p_rdiscrFine%p_rtriangulation
        SELECT CASE (elem_getPrimaryElement(ractProjection%ielementTypeProlongation))
        CASE (EL_P1_1D)
          ! P1 prolongation
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL mlprj_prolUniformP1_1D_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
               p_rtriaCoarse%NEL)

        CASE (EL_P2_1D)
          ! P2 prolongation
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL mlprj_prolUniformP2_1D_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse, p_rtriaCoarse%NVT, p_rtriaCoarse%NEL)

        CASE (EL_S31_1D)
          ! S31 prolongation
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
                               p_DvertexCoordsCoarse)
          CALL mlprj_prolUniformS31_1D_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
               p_DvertexCoordsCoarse, p_rtriaCoarse%NVT, &
               p_rtriaFine%NVT, p_rtriaCoarse%NEL)

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
               
        CASE (EL_P2)
          ! P2 prolongation
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                               p_IedgesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                               p_IedgesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_prolUniformP2_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse,&
               p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
               p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
               p_rtriaCoarse%NEL)

        CASE (EL_Q0)
          ! Q0 prolongation
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_prolUniformQ0_double (p_DuCoarse,p_DuFine, &
               p_IneighboursAtElementFine, &
               p_rtriaCoarse%NEL)
               
        CASE (EL_Q1)
          ! Q1 prolongation
          SELECT CASE (ractProjection%iprolVariant)
          CASE (:1) ! Standard prolongation
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtEdge, &
                                p_IverticesAtEdgeCoarse)
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                                p_IverticesAtElementCoarse)
            CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                                p_IverticesAtElementFine)
            CALL mlprj_prolUniformQ1_double (p_DuCoarse, p_DuFine, &
                      p_IverticesAtEdgeCoarse, p_IverticesAtElementCoarse, &
                      p_rtriaCoarse%NVT, p_rtriaCoarse%NMT, p_rtriaCoarse%NEL)
             ! 'old' implementation
!            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
!                                p_IverticesAtElementCoarse)
!            CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
!                                p_IverticesAtElementFine)
!            CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
!                                p_IneighboursAtElementCoarse)
!            CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
!                                p_IneighboursAtElementFine)
!            CALL mlprj_prolUniformQ1_double (p_DuCoarse,p_DuFine, &
!                p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
!                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,p_rtriaCoarse%NEL)

          CASE (2:) !Experimental FEAST MIRROR prolongation with zero boundary
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                                p_IverticesAtElementCoarse)
            CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                                p_IverticesAtElementFine)
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                                p_IneighboursAtElementCoarse)
            CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                                p_IneighboursAtElementFine)
            CALL mlprj_prolUnifQ1FMzero_double (p_DuCoarse,p_DuFine, &
                p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NEL,p_rtriaFine%NEL)
          END SELECT

        CASE (EL_Q2)
          ! Q2 prolongation
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                               p_IedgesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
                               
          CALL mlprj_prolUniformQ2_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementFine,p_IedgesAtElementFine,&
               p_IneighboursAtElementFine,&
               p_rtriaFine%NVT, p_rtriaFine%NMT, p_rtriaCoarse%NEL)  
               
        CASE (EL_QP1)        
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_prolUniformQP1_double (p_DuCoarse,p_DuFine, &
               p_IneighboursAtElementFine,p_rtriaCoarse%NEL,p_rtriaFine%NEL)                       
                       
        CASE (EL_Q1T)
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
          SELECT CASE (ractProjection%iprolVariant)
          CASE (:1) ! Standard prolongation
            IF (IAND(ractProjection%ielementTypeProlongation,INT(2**16,I32)) .NE. 0) THEN
              ! DOF's = integral mean values
              CALL mlprj_prolUniformEx30_double (p_DuCoarse,p_DuFine, &
                  p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                  p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                  p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)
            ELSE
              ! DOF's = edge midpoint based
              CALL mlprj_prolUniformEx31_double (p_DuCoarse,p_DuFine, &
                  p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                  p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                  p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)
            END IF
                
          CASE (2:) ! Extended prolongation; modified weights, local switch
                    ! to constant prolongation
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                                       p_IverticesAtElementCoarse)
            CALL storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
                                          p_DvertexCoordsCoarse)
            CALL storage_getbase_double(p_rtriaCoarse%h_DelementVolume, &
                                        p_DelementAreaCoarse)
            ! (what a nasty call...)                                       
            IF (IAND(ractProjection%ielementTypeProlongation,INT(2**16,I32)) .NE. 0) THEN
              ! DOF's = integral mean values
              CALL mlprj_prolUniformEx30ext_double (p_DuCoarse,p_DuFine, &
                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
                      p_DelementAreaCoarse,&
                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
                      MIN(4,ractProjection%iprolVariant)-2, &
                      ractProjection%dprolARboundEX3Y, &
                      ractProjection%iprolARIndicatorEX3Y)
            ELSE
              ! DOF's = edge midpoint based
              CALL mlprj_prolUniformEx31ext_double (p_DuCoarse,p_DuFine, &
                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
                      p_DelementAreaCoarse,&
                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
                      MIN(4,ractProjection%iprolVariant)-2, &
                      ractProjection%dprolARboundEX3Y, &
                      ractProjection%iprolARIndicatorEX3Y)
            END IF
          END SELECT
        
!        CASE (EL_E037)
!          ! Q2~ with bubble prolongation
!          CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
!                               p_IedgesAtElementFine)
!          CALL storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
!                               p_IedgesAtElementCoarse)
!          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
!                               p_IneighboursAtElementFine)
!          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
!                               p_IneighboursAtElementCoarse)
!          CALL storage_getbase_int(p_rtriaFine%h_ItwistIndexEdges, &
!                               p_ItwistIndexEdgesFine)
!          CALL storage_getbase_int(p_rtriaCoarse%h_ItwistIndexEdges, &
!                               p_ItwistIndexEdgesCoarse)
!
!          CALL mlprj_prolUniformE037_double (p_DuCoarse,p_DuFine, &
!              p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!              p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!              p_ItwistIndexEdgesCoarse,p_ItwistIndexEdgesFine,&
!              p_rtriaCoarse%NVT,p_rtriaFine%NVT,&
!              p_rtriaCoarse%NMT,p_rtriaFine%NMT,&
!              p_rtriaCoarse%NEL,p_rtriaFine%NEL)

        CASE (EL_Q0_3D)
          ! Q0 prolongation
          CALL mlprj_prolUniformQ0_3D_double (p_DuCoarse,p_DuFine, p_rtriaCoarse%NEL)

        CASE (EL_Q1_3D)
          ! Q1 prolongation
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtEdge, &
                              p_IverticesAtEdgeCoarse)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtFace, &
                              p_IverticesAtFaceCoarse)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                              p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                              p_IverticesAtElementFine)
          CALL mlprj_prolUniformQ1_3D_double (p_DuCoarse, p_DuFine, &
                    p_IverticesAtEdgeCoarse, p_IverticesAtFaceCoarse,&
                    p_IverticesAtElementCoarse, p_rtriaCoarse%NVT, &
                    p_rtriaCoarse%NMT, p_rtriaCoarse%NAT, p_rtriaCoarse%NEL)

        CASE (EL_Q1T_3D)
          ! Q1~ prolongation, DOF's = integral mean value
          CALL storage_getbase_int2d(p_rtriaFine%h_IfacesAtElement, &
                               p_IfacesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IfacesAtElement, &
                               p_IfacesAtElementCoarse)
          !CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
          !                     p_IneighboursAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
          
          ! Type of prolongation? Extended or not?
!          SELECT CASE (ractProjection%iprolVariant)
!          CASE (:1) ! Standard prolongation
!            IF (IAND(ractProjection%ielementTypeProlongation,INT(2**16,I32)) .NE. 0) THEN
!              ! DOF's = integral mean values
!              CALL mlprj_prolUniformEx30_double (p_DuCoarse,p_DuFine, &
!                  p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                  p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                  p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)
!            ELSE
              ! DOF's = face midpoint based
              CALL mlprj_prolUniformEx3x_3D_double (p_DuCoarse,p_DuFine, &
                  p_IfacesAtElementCoarse,p_IfacesAtElementFine,&
                  p_IneighboursAtElementCoarse,p_rtriaCoarse%NVT,p_rtriaFine%NVT,&
                  p_rtriaCoarse%NMT,p_rtriaFine%NMT,p_rtriaCoarse%NEL,&
                  ractProjection%ielementTypeProlongation)
!            END IF
!                
!          CASE (2:) ! Extended prolongation; modified weights, local switch
!                    ! to constant prolongation
!            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
!                                       p_IverticesAtElementCoarse)
!            CALL storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
!                                          p_DvertexCoordsCoarse)
!            CALL storage_getbase_double(p_rtriaCoarse%h_DelementVolume, &
!                                        p_DelementAreaCoarse)
!            ! (what a nasty call...)                                       
!            IF (IAND(ractProjection%ielementTypeProlongation,INT(2**16,I32)) .NE. 0) THEN
!              ! DOF's = integral mean values
!              CALL mlprj_prolUniformEx30ext_double (p_DuCoarse,p_DuFine, &
!                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
!                      p_DelementAreaCoarse,&
!                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
!                      MIN(4,ractProjection%iprolVariant)-2, &
!                      ractProjection%dprolARboundEX3Y, &
!                      ractProjection%iprolARIndicatorEX3Y)
!            ELSE
!              ! DOF's = edge midpoint based
!              CALL mlprj_prolUniformEx31ext_double (p_DuCoarse,p_DuFine, &
!                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
!                      p_DelementAreaCoarse,&
!                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
!                      MIN(4,ractProjection%iprolVariant)-2, &
!                      ractProjection%dprolARboundEX3Y, &
!                      ractProjection%iprolARIndicatorEX3Y)
!            END IF
!          END SELECT

        CASE DEFAULT
          PRINT *,'Unsupported prolongation!'
          CALL sys_halt()
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
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementCoarse,&
      p_IverticesAtElementFine,p_IverticesAtEdgeCoarse,p_IverticesAtFaceCoarse
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementCoarse,&
      p_IneighboursAtElementFine
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementCoarse,&
      p_IedgesAtElementFine, p_IfacesAtElementCoarse, p_IfacesAtElementFine
  INTEGER(I32), DIMENSION(:), POINTER :: p_ItwistIndexEdgesCoarse, &
      p_ItwistIndexEdgesFine
  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoordsCoarse
  REAL(DP), DIMENSION(:), POINTER   :: p_DelementAreaCoarse
  
  ! Data arrays
  REAL(DP), DIMENSION(:), POINTER :: p_DuCoarse, p_DuFine
 
    ! The vectors must be of data type DOUBLE - we don't support anything
    ! different at the moment...
    IF (rcoarseVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Coarse grid vector has unsupported data type!'
      CALL sys_halt()
    END IF

    IF (rfineVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Fine grid vector has unsupported data type!'
      CALL sys_halt()
    END IF
    
    IF (lsysbl_isVectorSorted(rfineVector) .OR. lsysbl_isVectorSorted(rcoarseVector)) THEN
      PRINT *,'Vectors must be unsorted for level change!'
      CALL sys_halt()
    END IF
    
    ! Calls the correct prolongation routine for each block in the 
    ! discretisation...
    DO i=1,rcoarseVector%nblocks
    
      IF (rcoarseVector%RvectorBlock(i)%NEQ .GT. 0) THEN

        ! Do we use L2-Projection here?
        IF (rprojection%RscalarProjection(1,i)%iprojType .EQ. &
            MLP_PROJ_TYPE_L2) THEN
          
          ! Call scalar L2-restriction
          CALL mlprj_restScalarL2(rprojection%RscalarProjection(1,i), &
            rcoarseVector%RvectorBlock(i), rfineVector%RvectorBlock(i))
        
          ! Continue with next block
          CYCLE
          
        END IF
      

        p_rdiscrCoarse => rcoarseVector%RvectorBlock(i)%p_rspatialDiscr
        p_rdiscrFine => rfineVector%RvectorBlock(i)%p_rspatialDiscr
        
        ! We need a discretisation:
        IF ((.NOT. ASSOCIATED(p_rdiscrCoarse)) .OR. &
            (.NOT. ASSOCIATED(p_rdiscrFine))) THEN
          PRINT *,'Intergrid transfer: No discretisation!'
          CALL sys_halt()
        END IF

        ! Currently, we support only uniform triangulations.
        IF ((p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
            (p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM)) THEN
          PRINT *,'Intergrid transfer supports currently only uniform discretisations!'
          CALL sys_halt()
        END IF
        
        ! Get the pointers to the vectors
        CALL lsyssc_getbase_double (rcoarseVector%RvectorBlock(i),p_DuCoarse)
        CALL lsyssc_getbase_double (rfineVector%RvectorBlock(i),p_DuFine)
        
        ! Use the first projection structure as template and create
        ! the actual projection structure for our situation.
        ! Remember, in a uniform grid we only have one projection structure
        ! and one element distribution!
        CALL mlprj_getProjectionStrategy (rprojection%RscalarProjection(1,i), &
              p_rdiscrCoarse%RelementDistr(1), &
              p_rdiscrFine%RelementDistr(1), &
              ractProjection)
      
        ! Depending on the element type of the trial functions in the
        ! discretisation, choose the right prolongation and call it.
        p_rtriaCoarse => p_rdiscrCoarse%p_rtriangulation
        p_rtriaFine => p_rdiscrFine%p_rtriangulation
        SELECT CASE (elem_getPrimaryElement(ractProjection%ielementTypeRestriction))
        CASE (EL_P1_1D)
          ! P1 restriction
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL mlprj_restUniformP1_1D_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse,p_IverticesAtElementFine, &
               p_rtriaCoarse%NEL)
               
        CASE (EL_P2_1D)
          ! P2 restriction
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL mlprj_restUniformP2_1D_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse,p_rtriaCoarse%NVT, p_rtriaCoarse%NEL)

        CASE (EL_S31_1D)
          ! S31 restriction
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                               p_IverticesAtElementFine)
          CALL storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
                               p_DvertexCoordsCoarse)
          CALL mlprj_restUniformS31_1D_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse,p_IverticesAtElementFine, &
               p_DvertexCoordsCoarse,p_rtriaCoarse%NVT, &
               p_rtriaFine%NVT, p_rtriaCoarse%NEL)

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

        CASE (EL_P2)
          ! P2 restriction
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                               p_IedgesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                               p_IedgesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_restUniformP2_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse,&
               p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
               p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
               p_rtriaCoarse%NEL)
          
        CASE (EL_Q0)
          ! Q0 restriction
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_restUniformQ0_double (p_DuCoarse,p_DuFine, &
               p_IneighboursAtElementFine, &
               p_rtriaCoarse%NEL)
               
        CASE (EL_Q1)
          ! Q1 restriction
          
          ! Type of restriction? 
          SELECT CASE (ractProjection%irestVariant)
          CASE (:0) ! Standard restriction
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtEdge, &
                                p_IverticesAtEdgeCoarse)
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                                p_IverticesAtElementCoarse)
            CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                                p_IverticesAtElementFine)
            CALL mlprj_restUniformQ1_double (p_DuCoarse, p_DuFine, &
                      p_IverticesAtEdgeCoarse, p_IverticesAtElementCoarse, &
                      p_rtriaCoarse%NVT, p_rtriaCoarse%NMT, p_rtriaCoarse%NEL)
             ! 'old' implementation
!            CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
!                                p_IverticesAtElementFine)
!            CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
!                                p_IneighboursAtElementFine)
!            CALL mlprj_restUniformQ1_double (p_DuCoarse,p_DuFine, &
!                p_IverticesAtElementFine,p_IneighboursAtElementFine,&
!                p_rtriaFine%NEL)
          CASE (1) ! All boundaries are FEAST mirror boundaries
            CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                                p_IverticesAtElementFine)
            CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                                p_IneighboursAtElementFine)
            CALL mlprj_restUniformQ1FM_double (p_DuCoarse,p_DuFine, &
                p_IverticesAtElementFine,p_IneighboursAtElementFine,&
                p_rtriaFine%NEL)
          CASE (2:) ! All boundaries are FEAST mirror boundaries with zero boundary
            CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                                p_IverticesAtElementFine)
            CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                                p_IneighboursAtElementFine)
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                                p_IverticesAtElementCoarse)
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                                p_IneighboursAtElementCoarse)
            CALL mlprj_restUnifQ1FMzero_double (p_DuCoarse,p_DuFine, &
                p_IverticesAtElementFine,p_IneighboursAtElementFine,&
                p_IverticesAtElementCoarse,p_IneighboursAtElementCoarse,&
                p_rtriaFine%NEL,p_rtriaCoarse%NEL)
            
          END SELECT
               
        CASE (EL_Q2)
          ! Q2 restriction
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                               p_IedgesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                               p_IedgesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
                               
          CALL mlprj_restUniformQ2_double (p_DuCoarse,p_DuFine, &
               p_IverticesAtElementCoarse, p_IedgesAtElementCoarse, &
               p_IedgesAtElementFine, p_IneighboursAtElementFine,&
               p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NMT,&
               p_rtriaFine%NMT,p_rtriaCoarse%NEL)
                       
        CASE (EL_QP1)       
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_restUniformQP1_double (p_DuCoarse,p_DuFine, &
               p_IneighboursAtElementFine, p_rtriaCoarse%NEL,p_rtriaFine%NEL)                              
                              
        CASE (EL_Q1T)
          ! Q1~ restriction
          CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                               p_IedgesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                               p_IedgesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)

          ! Type of restriction? Extended or not?
          SELECT CASE (ractProjection%irestVariant)
          CASE (:1) ! Standard prolongation
            IF (IAND(ractProjection%ielementTypeRestriction,INT(2**16,I32)) .NE. 0) THEN
              ! DOF's = integral mean values
              CALL mlprj_restUniformEx30_double (p_DuCoarse,p_DuFine, &
                  p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                  p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                  p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)          
            ELSE
              ! DOF's = edge midpoint based
              CALL mlprj_restUniformEx31_double (p_DuCoarse,p_DuFine, &
                  p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                  p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                  p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)          
            END IF
                
          CASE (2:) ! Extended prolongation; modified weights, local switch
                    ! to constant prolongation
            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                                       p_IverticesAtElementCoarse)
            CALL storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
                                          p_DvertexCoordsCoarse)
            CALL storage_getbase_double(p_rtriaCoarse%h_DelementVolume, &
                                        p_DelementAreaCoarse)
            ! (what a nasty call...)                                       
            IF (IAND(ractProjection%ielementTypeRestriction,INT(2**16,I32)) .NE. 0) THEN
              ! DOF's = integral mean values
              CALL mlprj_restUniformEx30ext_double (p_DuCoarse,p_DuFine, &
                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
                      p_DelementAreaCoarse,&
                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
                      MIN(4,ractProjection%irestVariant)-2, &
                      ractProjection%dprolARboundEX3Y, &
                      ractProjection%iprolARIndicatorEX3Y)
            ELSE
              ! DOF's = edge midpoint based
              CALL mlprj_restUniformEx31ext_double (p_DuCoarse,p_DuFine, &
                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
                      p_DelementAreaCoarse,&
                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
                      MIN(4,ractProjection%irestVariant)-2, &
                      ractProjection%dprolARboundEX3Y, &
                      ractProjection%iprolARIndicatorEX3Y)
            END IF
          END SELECT

!        CASE (EL_E037)
!          ! Q2~ with bubble prolongation
!          CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
!                               p_IedgesAtElementFine)
!          CALL storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
!                               p_IedgesAtElementCoarse)
!          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
!                               p_IneighboursAtElementFine)
!          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
!                               p_IneighboursAtElementCoarse)
!          CALL storage_getbase_int(p_rtriaFine%h_ItwistIndexEdges, &
!                               p_ItwistIndexEdgesFine)
!          CALL storage_getbase_int(p_rtriaCoarse%h_ItwistIndexEdges, &
!                               p_ItwistIndexEdgesCoarse)
!
!          CALL mlprj_restUniformE037_double (p_DuCoarse,p_DuFine, &
!              p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!              p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!              p_ItwistIndexEdgesCoarse,p_ItwistIndexEdgesFine,&
!              p_rtriaCoarse%NVT,p_rtriaFine%NVT,&
!              p_rtriaCoarse%NMT,p_rtriaFine%NMT,&
!              p_rtriaCoarse%NEL,p_rtriaFine%NEL)

        CASE (EL_Q0_3D)
          ! Q0 restriction
          CALL mlprj_restUniformQ0_3D_double (p_DuCoarse,p_DuFine, p_rtriaCoarse%NEL)

        CASE (EL_Q1_3D)
          ! Q1 prolongation
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtEdge, &
                              p_IverticesAtEdgeCoarse)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtFace, &
                              p_IverticesAtFaceCoarse)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                              p_IverticesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                              p_IverticesAtElementFine)
          CALL mlprj_restUniformQ1_3D_double (p_DuCoarse, p_DuFine, &
                    p_IverticesAtEdgeCoarse, p_IverticesAtFaceCoarse,&
                    p_IverticesAtElementCoarse, p_rtriaCoarse%NVT, &
                    p_rtriaCoarse%NMT, p_rtriaCoarse%NAT, p_rtriaCoarse%NEL)

        CASE (EL_Q1T_3D)
          ! Q1~ restriction
          CALL storage_getbase_int2d(p_rtriaFine%h_IfacesAtElement, &
                               p_IfacesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IfacesAtElement, &
                               p_IfacesAtElementCoarse)
          !CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
          !                     p_IneighboursAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)

!          ! Type of restriction? Extended or not?
!          SELECT CASE (ractProjection%irestVariant)
!          CASE (:1) ! Standard prolongation
!            IF (IAND(ractProjection%ielementTypeRestriction,INT(2**16,I32)) .NE. 0) THEN
!              ! DOF's = integral mean values
!              CALL mlprj_restUniformEx30_double (p_DuCoarse,p_DuFine, &
!                  p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                  p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                  p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)          
!            ELSE
              ! DOF's = face midpoint based
              CALL mlprj_restUniformEx3x_3D_double (p_DuCoarse,p_DuFine, &
                  p_IfacesAtElementCoarse,p_IfacesAtElementFine,&
                  p_IneighboursAtElementCoarse,p_rtriaCoarse%NVT,p_rtriaFine%NVT,&
                  p_rtriaCoarse%NMT,p_rtriaFine%NMT,p_rtriaCoarse%NEL,&
                  ractProjection%ielementTypeRestriction)          
!            END IF
!                
!          CASE (2:) ! Extended prolongation; modified weights, local switch
!                    ! to constant prolongation
!            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
!                                       p_IverticesAtElementCoarse)
!            CALL storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
!                                          p_DvertexCoordsCoarse)
!            CALL storage_getbase_double(p_rtriaCoarse%h_DelementVolume, &
!                                        p_DelementAreaCoarse)
!            ! (what a nasty call...)                                       
!            IF (IAND(ractProjection%ielementTypeRestriction,INT(2**16,I32)) .NE. 0) THEN
!              ! DOF's = integral mean values
!              CALL mlprj_restUniformEx30ext_double (p_DuCoarse,p_DuFine, &
!                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
!                      p_DelementAreaCoarse,&
!                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
!                      MIN(4,ractProjection%irestVariant)-2, &
!                      ractProjection%dprolARboundEX3Y, &
!                      ractProjection%iprolARIndicatorEX3Y)
!            ELSE
!              ! DOF's = edge midpoint based
!              CALL mlprj_restUniformEx31ext_double (p_DuCoarse,p_DuFine, &
!                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
!                      p_DelementAreaCoarse,&
!                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
!                      MIN(4,ractProjection%irestVariant)-2, &
!                      ractProjection%dprolARboundEX3Y, &
!                      ractProjection%iprolARIndicatorEX3Y)
!            END IF
!          END SELECT

        CASE DEFAULT
          PRINT *,'Unsupported restriction!'
          CALL sys_halt()
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
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementCoarse
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementFine
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementCoarse
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementFine
  INTEGER(PREC_FACEIDX), DIMENSION(:,:), POINTER :: p_IfacesAtElementCoarse
  INTEGER(PREC_FACEIDX), DIMENSION(:,:), POINTER :: p_IfacesAtElementFine
  
  ! Data arrays
  REAL(DP), DIMENSION(:), POINTER :: p_DuCoarse, p_DuFine
 
    ! The vectors must be of data type DOUBLE - we don't support anything
    ! different at the moment...
    IF (rcoarseVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Coarse grid vector has unsupported data type!'
      CALL sys_halt()
    END IF

    IF (rfineVector%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'Fine grid vector has unsupported data type!'
      CALL sys_halt()
    END IF
    
    IF (lsysbl_isVectorSorted(rfineVector) .OR. lsysbl_isVectorSorted(rcoarseVector)) THEN
      PRINT *,'Vectors must be unsorted for level change!'
      CALL sys_halt()
    END IF

    ! Calls the correct prolongation routine for each block in the 
    ! discretisation...
    DO i=1,rcoarseVector%nblocks
    
      IF (rcoarseVector%RvectorBlock(i)%NEQ .GT. 0) THEN
        p_rdiscrCoarse => rcoarseVector%RvectorBlock(i)%p_rspatialDiscr
        p_rdiscrFine => rfineVector%RvectorBlock(i)%p_rspatialDiscr
        
        ! We need a discretisation:
        IF ((.NOT. ASSOCIATED(p_rdiscrCoarse)) .OR. &
            (.NOT. ASSOCIATED(p_rdiscrFine))) THEN
          PRINT *,'Intergrid transfer: No discretisation!'
          CALL sys_halt()
        END IF

        ! Currently, we support only uniform triangulations.
        IF ((p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
            (p_rdiscrCoarse%ccomplexity .NE. SPDISC_UNIFORM)) THEN
          PRINT *,'Intergrid transfer supports currently only uniform discretisations!'
          CALL sys_halt()
        END IF
        
        ! Get the pointers to the vectors
        CALL lsyssc_getbase_double (rcoarseVector%RvectorBlock(i),p_DuCoarse)
        CALL lsyssc_getbase_double (rfineVector%RvectorBlock(i),p_DuFine)
        
        ! Use the first projection structure as template and create
        ! the actual projection structure for our situation.
        ! Remember, in a uniform grid we only have one projection structure
        ! and one element distribution!
        CALL mlprj_getProjectionStrategy (rprojection%RscalarProjection(1,i), &
              p_rdiscrCoarse%RelementDistr(1), &
              p_rdiscrFine%RelementDistr(1), &
              ractProjection)
      
        ! Depending on the element type of the trial functions in the
        ! discretisation, choose the right prolongation and call it.
        p_rtriaCoarse => p_rdiscrCoarse%p_rtriangulation
        p_rtriaFine => p_rdiscrFine%p_rtriangulation
        SELECT CASE (elem_getPrimaryElement(ractProjection%ielementTypeProlongation))
        CASE (EL_P1_1D)
          ! P1 interpolation
          CALL mlprj_interpUniformP1_1D_double (p_DuCoarse,p_DuFine,p_rtriaCoarse%NVT)

        CASE (EL_P2_1D)
          ! P2 interpolation
          CALL mlprj_interpUniformP2_1D_double (p_DuCoarse,p_DuFine,&
               p_rtriaCoarse%NVT, p_rtriaCoarse%NEL)

        CASE (EL_S31_1D)
          ! S31 interpolation
          CALL mlprj_interpUniS31_1D_double (p_DuCoarse,p_DuFine,&
               p_rtriaCoarse%NVT, p_rtriaFine%NVT)

        CASE (EL_P1)
          ! P1 interpolation
          CALL mlprj_interpUniformP1_double (p_DuCoarse,p_DuFine,p_rtriaCoarse%NVT)

        CASE (EL_P2)
          ! P2 interpolation
          CALL mlprj_interpUniformP2_double (p_DuCoarse,p_DuFine, &
                                           p_rtriaCoarse%NVT, p_rtriaCoarse%NMT)
                                           
        CASE (EL_Q0)
          ! Q0 interpolation
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_interpUniformQ0_double (p_DuCoarse,p_DuFine, &
               p_IneighboursAtElementFine, p_rtriaCoarse%NEL)
               
        CASE (EL_Q1)
          ! Q1 interpolation
          CALL mlprj_interpUniformQ1_double (p_DuCoarse,p_DuFine, p_rtriaCoarse%NVT)
          
        CASE (EL_Q2)
          ! Q2 interpolation
          CALL mlprj_interpUniformQ2_double (p_DuCoarse,p_DuFine, &
                     p_rtriaCoarse%NVT, p_rtriaCoarse%NMT, p_rtriaCoarse%NEL)          
                     
        CASE (EL_QP1)
          ! QP1 interpolation
          CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
          CALL mlprj_interpUniformQP1_double (p_DuCoarse,p_DuFine, &
                  p_IneighboursAtElementFine, p_rtriaCoarse%NEL,p_rtriaFine%NEL)          
                  
        CASE (EL_Q1T)
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

        CASE (EL_Q0_3D)
          ! Q0 interpolation
          CALL mlprj_interpUniformQ0_3D_double (p_DuCoarse,p_DuFine,p_rtriaCoarse%NEL)

        CASE (EL_Q1_3D)
          ! Q1 interpolation
          CALL mlprj_interpUniformQ1_3D_double (p_DuCoarse,p_DuFine,p_rtriaCoarse%NVT)

        CASE (EL_Q1T_3D)
          ! Q1~ interpolation, DOF's = integral mean values
          ! We use the same routine also for interpolating Ex31 solutions - there's
          ! not too much difference...
          CALL storage_getbase_int2d(p_rtriaFine%h_IfacesAtElement, &
                               p_IfacesAtElementFine)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IfacesAtElement, &
                               p_IfacesAtElementCoarse)
          CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
          CALL mlprj_interpUniformEx3x_3D_dbl (p_DuCoarse,p_DuFine, &
                  p_IfacesAtElementCoarse,p_IfacesAtElementFine,&
                  p_IneighboursAtElementCoarse,p_rtriaCoarse%NVT,p_rtriaFine%NVT,&
                  p_rtriaCoarse%NMT,p_rtriaFine%NMT,p_rtriaCoarse%NEL,&
                  ractProjection%ielementTypeInterpolation)          

        CASE DEFAULT
          PRINT *,'Unsupported interpolation!'
          CALL sys_halt()
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
  ! Support for 1D P1 element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformP1_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
                                            IvertsAtElemFine,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertsAtElemCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertsAtElemFine
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: Q2 = 0.5_DP
  INTEGER(PREC_ELEMENTIDX) :: iel
  REAL(DP) :: duh1,duh2

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    CALL lalg_copyVectorDble (DuCoarse,DuFine(1:SIZE(DuCoarse)))

    ! Loop over the elements
    DO iel=1,NELCoarse

      duh1=DuCoarse(IvertsAtElemCoarse(1,iel))
      duh2=DuCoarse(IvertsAtElemCoarse(2,iel))

      ! If IEL is the current element index of the coarse grid element,
      ! then it was refined into 2 new fine grid elements with indices
      ! IEL and NVT+IEL, where NVT is the number of coarse grid vertices.
      ! The 'new' vertice in the fine grid corresponding to the coarse grid
      ! element IEL is the second vertice of the fine grid element IEL
      ! (and the first of the fine grid element NVT+IEL).
      DuFine(IvertsAtElemFine(2,iel)) = Q2 * (duh1 + duh2)

    END DO

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformP1_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
               IvertsAtElemFine,NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine

  ! IverticesAtElement array (KVERT) on the caorse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertsAtElemCoarse
  
  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertsAtElemFine
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: Q2 = .5_DP
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_VERTEXIDX) :: ifgv, icgv1, icgv2
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    CALL lalg_copyVectorDble (DuFine(1:SIZE(DuCoarse)),DuCoarse)
    
    ! Loop over the elements to collect the missing additive contributions:
    DO iel=1,NELcoarse
    
      ! Get the 'new' fine grid vertice
      ifgv = IvertsAtElemFine(2,iel)
      
      ! Get the 'old' coarse grid vertices
      icgv1 = IvertsAtElemCoarse(1,iel)
      icgv2 = IvertsAtElemCoarse(2,iel)

      DuCoarse(icgv1) = DuCoarse(icgv1) + Q2 * DuFine(ifgv)
      DuCoarse(icgv2) = DuCoarse(icgv2) + Q2 * DuFine(ifgv)

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformP1_1D_double (DuCoarse,DuFine,NVTcoarse)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse
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
  ! Support for 1D P2 element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformP2_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
                                            NVTcoarse,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertsAtElemCoarse

  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: P_3_8 = 0.375_DP
  REAL(DP), PARAMETER :: P_1_8 = 0.125_DP
  REAL(DP), PARAMETER :: P_3_4 = 0.75_DP
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER :: NVTfine
  REAL(DP) :: dv1,dv2,del

    ! Copy all DOFs from the coarse grid into the fine grid - 
    ! the DOFs belonging to the 'new' fine grid vertices get
    ! their values from the edge midpoints in the coarse grid.
    CALL lalg_copyVectorDble (DuCoarse,DuFine(1:SIZE(DuCoarse)))
    
    ! Calculate number of vertices in fine grid
    NVTfine = NVTcoarse+NELcoarse

    ! Loop over the elements
    DO iel=1,NELCoarse
    
      ! Get the DOFs from the coarse grid vertices...
      dv1 = DuCoarse(IvertsAtElemCoarse(1,iel))
      dv2 = DuCoarse(IvertsAtElemCoarse(2,iel))
      
      ! ...and get the DOF from the coarse grid edge midpoint
      del = DuCoarse(NVTcoarse + iel)
      
      ! Perform quadratic interpolation to calculate the DOFs for the
      ! fine grid edge midpoints.
      DuFine(NVTfine + iel) = P_3_8*dv1 - P_1_8*dv2 + P_3_4*del
      DuFine(NVTfine + NELcoarse + iel) = -P_1_8*dv1 + P_3_8*dv2 + P_3_4*del

    END DO

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformP2_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
               NVTcoarse,NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine

  ! IverticesAtElement array (KVERT) on the caorse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertsAtElemCoarse
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: P_3_8 = 0.375_DP
  REAL(DP), PARAMETER :: P_1_8 = 0.125_DP
  REAL(DP), PARAMETER :: P_3_4 = 0.75_DP
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER :: NVTfine
  REAL(DP) :: dem1,dem2
  INTEGER(PREC_VERTEXIDX) :: icgv1, icgv2,icgem
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    CALL lalg_copyVectorDble (DuFine(1:SIZE(DuCoarse)),DuCoarse)
    
    ! Calculate number of vertices in fine grid
    NVTfine = NVTcoarse+NELcoarse

    ! Loop over the elements to collect the missing additive contributions:
    DO iel=1,NELcoarse
    
      ! Get the fine grid edge midpoints
      dem1 = DuFine(NVTfine+iel)
      dem2 = DuFine(NVTfine+NELcoarse+iel)

      icgv1 = IvertsAtElemCoarse(1,iel)
      icgv2 = IvertsAtElemCoarse(2,iel)
      icgem = NVTcoarse+iel
      
      ! distribute the DOFs
      DuCoarse(icgv1) = DuCoarse(icgv1) + P_3_8*dem1 - P_1_8*dem2
      DuCoarse(icgv2) = DuCoarse(icgv2) - P_1_8*dem1 + P_3_8*dem2
      DuCoarse(icgem) = DuCoarse(icgem) + P_3_4*(dem1 + dem2)

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformP2_1D_double (DuCoarse,DuFine,NVTcoarse,NELcoarse)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  INTEGER :: len
  
    len = NVTcoarse+NELcoarse
  
    ! The first coase.NVT entries of the fine grid vector define 
    ! the values on the coarse grid - because of the two-level ordering!
    CALL lalg_copyVectorDble(DuFine(1:len),DuCoarse(1:len))
    
  END SUBROUTINE
  
  ! ***************************************************************************
  ! Support for 1D S31 element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformS31_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
             IvertsAtElemFine,DvertexCoordsCoarse,NVTcoarse,NVTfine,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $S_31$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertsAtElemCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertsAtElemFine
  
  ! DvertexCoords array (DCORVG) on the coarse grid
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: DvertexCoordsCoarse
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse
  
  ! Number of vertices in the fine grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: Q12 = 0.5_DP
  REAL(DP), PARAMETER :: Q14 = 0.25_DP
  REAL(DP), PARAMETER :: Q34 = 0.75_DP
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_VERTEXIDX) :: ivt1, ivt2,ivt
  REAL(DP) :: dp1,dp2,dq1,dq2,ddet
  
    ! Copy the first NVTcoarse entries - these are the coefficients of the
    ! basis functions for the function values.
    CALL lalg_copyVectorDble(DuCoarse(1:NVTcoarse),DuFine(1:NVTcoarse))
    
    ! Copy NVTcoarse entries beginning at NVTcoarse+1 of the coarse
    ! vector into the fine vector beginning at NVTfine+1 - these are the
    ! coefficients for the function derivative values.
    CALL lalg_copyVectorDble(DuCoarse(NVTcoarse + 1 : 2*NVTcoarse), &
                  DuFine(NVTfine + 1 : NVTfine + NVTcoarse))

    ! Loop over the elements
    DO iel=1,NELCoarse

      ! Vertices of element in coarse vector
      ivt1 = IvertsAtElemCoarse(1,iel)
      ivt2 = IvertsAtElemCoarse(2,iel)
      
      ! Calculate the determinant of this line
      ddet = 0.5 * (DvertexCoordsCoarse(1,ivt2) - DvertexCoordsCoarse(1,ivt1))

      ! Function value coefficients
      dp1=DuCoarse(ivt1)
      dp2=DuCoarse(ivt2)
      ! Function derivative coefficients
      dq1=DuCoarse(ivt1 + NVTcoarse)
      dq2=DuCoarse(ivt2 + NVTcoarse)

      ! Refined Vertice in fine vector
      ivt = IvertsAtElemFine(2,iel)
      
      ! Function Value
      DuFine(ivt) = Q12 * (dp1 + dp2) + Q14 * ddet * (dq1 - dq2)
      ! Function Derivative
      DuFine(ivt+NVTfine) = Q34 * (dp2 - dp1) / ddet - Q14 * (dq1 + dq2)

    END DO


  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformS31_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
             IvertsAtElemFine,DvertexCoordsCoarse,NVTcoarse,NVTfine,NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine

  ! IverticesAtElement array (KVERT) on the caorse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertsAtElemCoarse
  
  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertsAtElemFine
  
  ! DvertexCoords array (DCORVG) on the coarse grid
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: DvertexCoordsCoarse
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse
  
  ! Number of vertices in the fine grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP), PARAMETER :: Q12 = 0.5_DP
  REAL(DP), PARAMETER :: Q14 = 0.25_DP
  REAL(DP), PARAMETER :: Q34 = 0.75_DP
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_VERTEXIDX) :: ivtp1, ivtp2, ivtq1, ivtq2,ivt
  REAL(DP) :: dpf,dqf,ddet
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.

    ! Copy the first NVTcoarse entries - these are the coefficients of the
    ! basis functions for the function values.
    CALL lalg_copyVectorDble(DuFine(1:NVTcoarse),DuCoarse(1:NVTcoarse))
    
    ! Copy NVTcoarse entries beginning at NVTfine+1 of the fine vector into
    ! the coarse vector beginning at NVTfine+1 - these are the
    ! coefficients for the function derivative values.
    CALL lalg_copyVectorDble(DuFine(NVTfine + 1 : NVTfine + NVTcoarse),&
                             DuCoarse(NVTcoarse + 1 : 2*NVTcoarse))

    ! Loop over the elements to collect the missing additive contributions:
    DO iel=1,NELcoarse
    
      ! Vertices of element in coarse vector
      ivtp1 = IvertsAtElemCoarse(1,iel)
      ivtp2 = IvertsAtElemCoarse(2,iel)
      ivtq1 = ivtp1+NVTcoarse
      ivtq2 = ivtp2+NVTcoarse

      ! Calculate the determinant of this line
      ddet = 0.5 * (DvertexCoordsCoarse(1,ivtp2) - DvertexCoordsCoarse(1,ivtp1))

      ! Refined Vertice in fine vector
      ivt = IvertsAtElemFine(2,iel)
      dpf = DuFine(ivt)
      dqf = DuFine(ivt+NVTfine)

      ! Function values
      DuCoarse(ivtp1) = DuCoarse(ivtp1) + Q12*dpf - Q34*dqf/ddet
      DuCoarse(ivtp2) = DuCoarse(ivtp2) + Q12*dpf + Q34*dqf/ddet
      
      ! Function derivatives
      DuCoarse(ivtq1) = DuCoarse(ivtq1) + Q14*(dpf*ddet - dqf)
      DuCoarse(ivtq2) = DuCoarse(ivtq2) - Q14*(dpf*ddet + dqf)

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniS31_1D_double (DuCoarse,DuFine,NVTcoarse,NVTfine)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse
  
  ! Number of vertices in the fine grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTfine
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  
    ! The first NVTcoarse entries of the fine grid vector define 
    ! the values on the coarse grid - because of the two-level ordering!
    CALL lalg_copyVectorDble(DUfine(1:NVTcoarse), DUcoarse(1:NVTCoarse))
    CALL lalg_copyVectorDble(DUfine(NVTfine + 1 : NVTfine + NVTcoarse),&
                             DUcoarse(NVTCoarse + 1 : 2*NVTCoarse))
    
  END SUBROUTINE

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
  CALL sys_halt()

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
  CALL sys_halt()

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
  CALL sys_halt()
    
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
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
  
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
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
  
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
  INTEGER(PREC_VERTEXIDX) :: ih1
  
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
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse
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

  SUBROUTINE mlprj_prolUniformP2_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_2$, uniform triangulation, double precision vector.
!</description>
  
!<input>  
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse

  ! IedgesAtElement array (KMID) on the coarse grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementCoarse
  
  ! IedgesAtElement array (KMID) on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine
  
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

  ! local variables
  INTEGER(PREC_ELEMENTIDX) :: IEL,IEL1,IEL2
  REAL(DP) :: duc1,duc2,duc3,dum1,dum2,dum3
  REAL(DP), PARAMETER :: Q8 = 0.125_DP
  REAL(DP), PARAMETER :: Q4 = 0.25_DP
  REAL(DP), PARAMETER :: Q2 = 0.5_DP

    ! First, we remember the refinement scheme to clarify the
    ! local numbering of the vertices in the coarse- and fine-
    ! grid triangles.
    ! Let a coarse grid triangle be locally numbered as:
    !
    !   2 
    !   |  \
    !   |    \
    !   | IEL  \
    !   |        \
    !   3----------1
    ! 
    ! Then the refinement process assigns the following numbers:
    !
    !   2
    !   |  \
    !   |     \
    !   |        \
    !   2*-------- 1* 
    !   | \   IEL  |   \
    !   |    \     |      \
    !   |       \  |         \
    !   3----------3*----------1
    !
    ! i.e. the element number "IEL" is put into the middle with
    ! the corner vertices numbered according to the local edge
    ! numbers of the coarse grid element.
    !
    ! To access information on the edges of the fine grid element,
    ! we have to work with adjacencies!
    !  
    ! First copy DuCoarse to DuFine. This will transfer the corner values
    ! from the coarse grid to the fine grid. More precisely, this
    ! will map:
    !   Coarse grid vertices 1..NVTC  -> Fine grid vertices 1..NVTC
    !   Coarse grid midpoints 1..NMTC -> Fine grid vertices NVTC+1..NVTC+NMTC = NVTF
    ! Afterwards, we only have to create the missing midpoint values!

    CALL lalg_copyVectorDble (DuCoarse,DuFine(1:SIZE(DuCoarse)))
    
    ! loop over the elements

    DO IEL = 1,nelCoarse
    
      ! We fetch the function values of the coarse grid element
      ! into variables following the following scheme:
      !
      !     DUC2
      !     |   \
      !     |      \
      !     |         \
      !     DUM2         DUM1
      !     |                \
      !     |                   \
      !     |                      \
      !     DUC3 --------DUM3-------  DUC1

      DUC1 = DuCoarse(IverticesAtElementCoarse(1,IEL))
      DUC2 = DuCoarse(IverticesAtElementCoarse(2,IEL))
      DUC3 = DuCoarse(IverticesAtElementCoarse(3,IEL))
      DUM1 = DuCoarse(IedgesAtElementCoarse(1,IEL))
      DUM2 = DuCoarse(IedgesAtElementCoarse(2,IEL))
      DUM3 = DuCoarse(IedgesAtElementCoarse(3,IEL))

      ! We have to calculate the function values in the new DOF's on the
      ! fine grid, which are located here:
      !
      !     DUC2
      !     |   \
      !     X     X
      !     |        \
      !     DUM2 --X-- DUM1
      !     |   \      |   \
      !     X     X    X     X
      !     |       \  |        \
      !     DUC3---X---DUM3---X---DUC1
      !
      ! On this trangle, the function is a quadratic polynomial:
      !
      !   P(X,Y) = c1 + c2*x + c3*y + c4*x^2 + c5*x*y + c6*y^2
      !
      ! Solving for the coefficients such that the polynomial takes
      ! the DUxy-values in the corners/midpoints, we obtain the 
      ! polynomial as:
      !
      !   P(X,Y) = DUC3 + 
      !             (-3*DUC3-DUC1+4*DUM3)*x + 
      !             (-3*DUC3-DUC2+4*DUM2)*y + 
      !             (4*DUC3-4*DUM3+4*DUM1-4*DUM2)*x*y + 
      !             (2*DUC3+2*DUC1-4*DUM3)*x^2 + 
      !             (2*DUC3+2*DUC2-4*DUM2)*y^2
      !
      ! This has to be evaluated in the new points, marked as "X"
      ! in the above sketch.
      !
      ! Remember, that the coarse grid element IEL is moved
      ! to the inner fine-grid element. The corners of the coars
      ! grid element are always the first vertices of the triangles
      ! on the fine grid elements.
      !
      ! First calculate the edge mitpoint values of that one!
      !
      !   |        \
      !   DUM2---X---DUM1
      !   | \   IEL  |   \
      !         X     X        
      !           \  |           
      !           --DUM3--        
      !
      ! DUF(IedgesAtElementFine(1,IEL)) = P(1/4,1/2)
      !                   = -1/8*DUC3-1/8*DUC1+1/4*DUM3+1/2*DUM2+1/2*DUM1
      ! DUF(IedgesAtElementFine(2,IEL)) = P(1/4,1/2)
      !                   = -1/8*DUC1+1/2*DUM3-1/8*DUC2+1/2*DUM2+1/4*DUM1
      ! DUF(IedgesAtElementFine(3,IEL)) = P(1/2,1/4)
      !                   = -1/8*DUC3+1/2*DUM3-1/8*DUC2+1/4*DUM2+1/2*DUM1

      DuFine(IedgesAtElementFine(1,IEL)) = -Q8*DUC3-Q8*DUC1+Q4*DUM3+Q2*DUM2+Q2*DUM1
      DuFine(IedgesAtElementFine(2,IEL)) = -Q8*DUC1+Q2*DUM3-Q8*DUC2+Q2*DUM2+Q4*DUM1
      DuFine(IedgesAtElementFine(3,IEL)) = -Q8*DUC3+Q2*DUM3-Q8*DUC2+Q4*DUM2+Q2*DUM1

      ! Now to the values on the other edges.
      ! Only calculate information on edges that are adjacent to
      ! an element with lower element number. This prevents
      ! information from being computed twice.
      !
      ! Check the first edge:

      IF (IneighboursAtElementCoarse(1,IEL).LT.IEL) THEN
      
        ! Neighbor element has a smaller number than our element.
        ! Calculate the information on the current edge.
        ! This is the edge corresponding to the edge midpoint DUM1!
        ! We have to calculate the values in the new edge midpoints.
        !
        ! DUC2
        !  |   \
        !  |     X
        !  |  IEL1  \
        ! DUM2 ----- DUM1
        !      \ IEL  |   \
        !        \    |     X
        !          \  |  IEL2  \
        !            DUM3-------DUC1
        !
        ! Use adjacencies to get the fine-grid elements IEL1 and IEL2.

        IEL1 = IneighboursAtElementFine(1,IEL)
        IEL2 = IneighboursAtElementFine(3,IEL)
        
        ! Calculate the new edge midpoints:
        !
        !  DUF(IedgesAtElementFine(3,IEL1)) = P(1/4,3/4)
        !                     = -1/8*DUC1+3/8*DUC2+3/4*DUM1
        !  DUF(IedgesAtElementFine(1,IEL2)) = P(3/4,1/4)
        !                     = 3/8*DUC1-1/8*DUC2+3/4*DUM1
        
        DuFine(IedgesAtElementFine(3,IEL1)) = -Q8*DUC1+3.0_DP*Q8*DUC2+3.0_DP*Q4*DUM1
        DuFine(IedgesAtElementFine(1,IEL2)) = 3.0_DP*Q8*DUC1-Q8*DUC2+3.0_DP*Q4*DUM1
      
      END IF
    
      ! Check the next edge:

      IF (IneighboursAtElementCoarse(2,IEL).LT.IEL) THEN
      
        ! DUC2
        !  |   \
        !  X     \
        !  |   IEL1 \
        ! DUM2 ----- DUM1
        !  |   \      |     
        !  X     \IEL |       
        !  |  IEL2 \  |          
        ! DUC3-------DUM3               
        ! 
        ! Use adjacencies to get the fine-grid elements IEL1 and IEL2.

        IEL1 = IneighboursAtElementFine(1,IEL)
        IEL2 = IneighboursAtElementFine(2,IEL)
        
        ! Calculate the new edge midpoints:
        !  
        !  DUF(IedgesAtElementFine(1,IEL1)) = P(0,3/4)
        !                     = -1/8*DUC3+3/8*DUC2+3/4*DUM2
        !  DUF(IedgesAtElementFine(3,IEL2)) = P(0,1/4)
        !                     = 3/8*DUC3-1/8*DUC2+3/4*DUM2
        
        DuFine(IedgesAtElementFine(1,IEL1)) = &
            -Q8*DUC3+3.0_DP*Q8*DUC2+3.0_DP*Q4*DUM2
        DuFine(IedgesAtElementFine(3,IEL2)) = &
            3.0_DP*Q8*DUC3-1.0_DP*Q8*DUC2+3.0_DP*Q4*DUM2
      
      END IF

      ! Check the last edge

      IF (IneighboursAtElementCoarse(3,IEL).LT.IEL) THEN
      
        ! DUM2 ----- DUM1
        !  |   \  IEL |   \
        !  | IEL2\    | IEL1\
        !  |       \  |        \
        ! DUC3---X---DUM3---X---DUC1
        !
        ! Use adjacencies to get the fine-grid elements IEL1 and IEL2.

        IEL1 = IneighboursAtElementFine(3,IEL)
        IEL2 = IneighboursAtElementFine(2,IEL)
        
        ! Calculate the new edge midpoints:
        !
        !  DUF(IedgesAtElementFine(3,IEL1)) = P(3/4,0)
        !                     = -1/8*DUC3+3/8*DUC1+3/4*DUM3
        !  DUF(IedgesAtElementFine(1,IEL2)) = P(1/4,0)
        !                     = 3/8*DUC3-1/8*DUC1+3/4*DUM3
        
        DuFine(IedgesAtElementFine(3,IEL1)) = -Q8*DUC3+3.0_DP*Q8*DUC1+3.0_DP*Q4*DUM3
        DuFine(IedgesAtElementFine(1,IEL2)) = 3.0_DP*Q8*DUC3-Q8*DUC1+3.0_DP*Q4*DUM3
      
      END IF
      
      !  That's it - next element.
      
    END DO ! IEL

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformP2_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_2$, uniform triangulation, double precision vector.
!</description>
  
!<input>  
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse

  ! IedgesAtElement array (KMID) on the coarse grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementCoarse
  
  ! IedgesAtElement array (KMID) on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine
  
  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse

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
    INTEGER(PREC_ELEMENTIDX) :: IEL,IEL2,IEL3,IEL4
    REAL(DP) :: dn1,dn2,dn3
    REAL(DP) :: duf11,duf12,duf13,duf21,duf23,duf31,duf33
    REAL(DP) :: duf41,DUF43
    REAL(DP), PARAMETER :: Q8 = 0.125_DP
    REAL(DP), PARAMETER :: Q4 = 0.25_DP
    REAL(DP), PARAMETER :: Q2 = 0.5_DP

    ! First we remember the refinement scheme to clarify the
    ! local numbering of the vertices in the coarse- and fine-
    ! grid triangles.
    ! Let a coarse grid triangle be locally numbered as:
    !
    !   2 
    !   |  \
    !   |    \
    !   | IEL  \
    !   |        \
    !   3----------1
    ! 
    ! Then the refinement process assigns the following numbers:
    !
    !   2
    !   |  \
    !   |     \
    !   |        \
    !   2*-------- 1* 
    !   | \   IEL  |   \
    !   |    \     |      \
    !   |       \  |         \
    !   3----------3*----------1
    !
    ! i.e. the element number "IEL" is put into the middle with
    ! the corner vertices numbered according to the local edge
    ! numbers of the coarse grid element.
    !
    ! To access information on the edges of the fine grid element,
    ! we have to work with adjacencies!
    !
    ! Copy the first NVTC+NMTC values from DUF to DUC. This will
    ! transfer the contribution of the values from the 
    ! fine-grid vertices that are coarse-grid vertices or midpoints
    ! as well. More precisely, this will transfer:
    !
    !   Fine grid vertices 1..NVTC -> Coarse grid vertices 1..NVTC  
    !   Fine grid vertices NVTC+1..NVTC+NMTC = NVTF -> Coarse grid midpoints 1..NMTC
    !
    ! Afterwards, we have to add only the contribution of the fine grid
    ! edge midpoints to the coarse grid vertices/midpoints.

    CALL lalg_copyVectorDble (DuFine(1:SIZE(DuCoarse)),DuCoarse)
      
    ! loop over the elements

    DO IEL = 1,nelCoarse
    
      ! The prolongation created from DUC1-3 and DUM1-3 in the 
      ! following sketch the fine grid values DUFxy:
      !
      !    DUC2
      !     |  \
      !     |    \
      !   DUF21  DUF23
      !     |        \
      !     |          \
      !    DUM2--DUF11--DUM1
      !     |  \         |  \
      !     |    \       |    \
      !   DUF33  DUF12 DUF13  DUF41
      !     |        \   |        \
      !     |          \ |          \
      !    DUC3--DUF31--DUM3--DUF43--DUC1
      !
      ! This was done by a weighted averaging. The restriction now
      ! builds the DUC1-3 and DUM1-3 values as a weighted sum
      ! of themselves (DUM1-3, DUC1-3) and the new midpoints
      ! DUFxy using the same weights.
      !
      ! We had 
      !   DUCx(fine grid) := 1*DUCx(coarse grid)
      !   DUMy(fine grid) := 1*DUCy(coarse grid)
      ! Therefore:
      !   DUCx(coarse grid) := 1*DUCx(fine grid) + ...
      !   DUCy(coarse grid) := 1*DUMy(fine grid) + ...
      !
      ! This part of the sum is already written to DUC with the
      ! above LCP1 command! now comes the rest of the sum.
      !
      ! The prolongation used the following formulas to calculate
      ! the fine grid vertices:
      !
      !  DUF11 = P(1/4,1/2)
      !        = -1/8*DUC3-1/8*DUC1+1/4*DUM3+1/2*DUM2+1/2*DUM1
      !  DUF12 = P(1/4,1/4)
      !        = -1/8*DUC1+1/2*DUM3-1/8*DUC2+1/2*DUM2+1/4*DUM1
      !  DUF13 = P(1/2,1/4)
      !        = -1/8*DUC3+1/2*DUM3-1/8*DUC2+1/4*DUM2+1/2*DUM1
      !
      !  DUF23 = P(1/4,3/4)
      !        = -1/8*DUC1+3/8*DUC2+3/4*DUM1
      !  DUF41 = P(3/4,1/4)
      !        = 3/8*DUC1-1/8*DUC2+3/4*DUM1
      !
      !  DUF21 = P(0,3/4)
      !        = -1/8*DUC3+3/8*DUC2+3/4*DUM2
      !  DUF33 = P(0,1/4)
      !        = 3/8*DUC3-1/8*DUC2+3/4*DUM2
      !
      !  DUF43 = P(3/4,0)
      !        = -1/8*DUC3+3/8*DUC1+3/4*DUM3
      !  DUF31 = P(1/4,0)
      !        = 3/8*DUC3-1/8*DUC1+3/4*DUM3
      !
      ! This is equivalent to the system
      !
      !  DUF11    [-1/8,    0, -1/8, 1/2,  1/2,  1/4]   DUC1
      !  DUF12    [-1/8, -1/8 ,   0, 1/4,  1/2,  1/2]   DUC2
      !  DUF13    [   0, -1/8, -1/8, 1/2,  1/4,  1/2]   DUC3
      !  DUF21    [   0,  3/8, -1/8,   0,  3/4,    0]   DUM1
      !  DUF23 =  [-1/8,  3/8,    0, 3/4,    0,    0] * DUM2
      !  DUF31    [-1/8,    0,  3/8,   0,    0,  3/4]   DUM3
      !  DUF33    [   0, -1/8,  3/8,   0,  3/4,    0]
      !  DUF41    [ 3/8, -1/8 ,   0, 3/4,    0,    0]
      !  DUF43    [ 3/8,    0, -1/8,   0,    0,  3/4]
      !
      ! Transposing it says what is left to add to DUC1-3/DUM1-3:
      !
      !   DUC1    [-1/8, -1/8,    0,    0, -1/8, -1/8,    0,  3/8,  3/8]   DUF11
      !   DUC2    [   0, -1/8, -1/8,  3/8,  3/8,    0, -1/8, -1/8,    0]   DUF12
      !   DUC3 += [-1/8,    0, -1/8, -1/8,    0,  3/8,  3/8,    0, -1/8] * DUF13
      !   DUM1    [ 1/2,  1/4,  1/2,    0,  3/4,    0,    0,  3/4,    0]   DUF21
      !   DUM2    [ 1/2,  1/2,  1/4, 3/4,     0,    0,  3/4,    0,    0]   DUF23
      !   DUM3    [ 1/4,  1/2,  1/2,    0,    0,  3/4,    0,    0,  3/4]   DUF31
      !                                                                    DUF33
      !                                                                    DUF41
      !                                                                    DUF43
      !
      ! Fetch the fine grid values into local variables:

      IEL2 = IneighboursAtElementFine(1,IEL)
      IEL3 = IneighboursAtElementFine(2,IEL)
      IEL4 = IneighboursAtElementFine(3,IEL)
      DUF11 = DuFine(IedgesAtElementFine(1,IEL))
      DUF12 = DuFine(IedgesAtElementFine(2,IEL))
      DUF13 = DuFine(IedgesAtElementFine(3,IEL))
      DUF21 = DuFine(IedgesAtElementFine(1,IEL2))
      DUF23 = DuFine(IedgesAtElementFine(3,IEL2))
      DUF31 = DuFine(IedgesAtElementFine(1,IEL3))
      DUF33 = DuFine(IedgesAtElementFine(3,IEL3))
      DUF41 = DuFine(IedgesAtElementFine(1,IEL4))
      DUF43 = DuFine(IedgesAtElementFine(3,IEL4))

      ! When we add the information to DUC1-3/DUM1-3 we have to take
      ! into account whether there is a neighbor element or not!
      ! If there is a neighbor, we only add half of the information
      ! of the edge with the neighbor to DUC1-3/DUM1-3. 
      ! In the loop over the elements here, we will
      ! later reach the neighbor and add another time half of the
      ! information, which that way completes that edge.
      !
      ! Set Nx=0.5 if there is a neighbor on edge x on the
      ! coarse grid or =1.0 if there is no neighbor

      dn1 = 0.5_DP
      dn2 = 0.5_DP
      dn3 = 0.5_DP
      
      IF (IneighboursAtElementCoarse(1,IEL) .EQ. 0) dn1 = 1.0_DP
      IF (IneighboursAtElementCoarse(2,IEL) .EQ. 0) dn2 = 1.0_DP
      IF (IneighboursAtElementCoarse(3,IEL) .EQ. 0) dn3 = 1.0_DP
      
      ! Now sum up the restriction.

      ! DUC1:

      DuCoarse(IverticesAtElementCoarse(1,IEL)) = DuCoarse(IverticesAtElementCoarse(1,IEL)) &
             -Q8*DUF11 -Q8*DUF12                        &
        +dn1*(-Q8*DUF23 +3.0_DP*Q8*DUF41)               &
        +dn3*(-Q8*DUF31 +3.0_DP*Q8*DUF43)               
        
      ! DUC2:

      DuCoarse(IverticesAtElementCoarse(2,IEL)) = DuCoarse(IverticesAtElementCoarse(2,IEL)) &
             -Q8*DUF12 -Q8*DUF13                        &
        +dn2*(3.0_DP*Q8*DUF21-Q8*DUF33)                 &
        +dn1*(3.0_DP*Q8*DUF23-Q8*DUF41)                    
    
      ! DUC3:

      DuCoarse(IverticesAtElementCoarse(3,IEL)) = DuCoarse(IverticesAtElementCoarse(3,IEL)) &
             -Q8*DUF11 -Q8*DUF13                        &
        +dn2*(-Q8*DUF21 +3.0_DP*Q8*DUF33)               &
        +dn3*(-Q8*DUF43 +3.0_DP*Q8*DUF31)                  
    
      ! DUM1:

      DuCoarse(IedgesAtElementCoarse(1,IEL)) = DuCoarse(IedgesAtElementCoarse(1,IEL)) &
       +     Q2*DUF11 +Q4*DUF12 +Q2*DUF13                                             &
       + dn1*(3.0_DP*Q4*DUF23 +3*Q4*DUF41)

      ! DUM2:

      DuCoarse(IedgesAtElementCoarse(2,IEL)) = DuCoarse(IedgesAtElementCoarse(2,IEL)) &
       +     Q2*DUF11 +Q2*DUF12 +Q4*DUF13                                             &
       + dn2*(3.0_DP*Q4*DUF21+3.0_DP*Q4*DUF33)

      ! DUM3:

      DuCoarse(IedgesAtElementCoarse(3,IEL)) = DuCoarse(IedgesAtElementCoarse(3,IEL)) &
       +     Q4*DUF11 +Q2*DUF12 +Q2*DUF13                                             &
       +dn3*(3.0_DP*Q4*DUF43 +3.0_DP*Q4*DUF31)

      ! That's it - next element.
      
    END DO ! IEL

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
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse

  ! Number of edges in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NMTcoarse
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
  INTEGER(PREC_VERTEXIDX) :: IELH1,IELH2,IELH3,IELH4
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
      
      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      DuCoarse(iel)= DuFine(IELH1)+DuFine(IELH2)+DuFine(IELH3)+DuFine(IELH4)
      
    END DO
    
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformQ0_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)
  
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
    DO iel=1,NELcoarse
    
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
               IverticesAtEdgeCoarse, IverticesAtElementCoarse, &
               NVTcoarse, NMTcoarse, NELcoarse)
! 'old' parameter list
!               IverticesAtElementCoarse,IverticesAtElementFine,&
!               IneighboursAtElementCoarse,IneighboursAtElementFine,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $Q_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IverticesAtEdge array on the coarse grid.
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdgeCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse

!  ! IverticesAtElement array (KVERT) on the fine grid
!  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
!  
!  ! IneighboursAtElement array on the coarse grid
!  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
!  
!  ! IneighboursAtElement array on the fine grid
!  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine

  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse

  ! Number of edges in the coarse grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTcoarse

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_EDGEIDX) :: iedge

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    CALL lalg_copyVectorDble (DuCoarse,DuFine(1:NVTcoarse))

    ! Loop over the edges
    DO iedge = 1, NMTcoarse
      ! Calculate the edge midpoint DOF
      DuFine(NVTcoarse + iedge) = 0.5_DP * (&
          DuCoarse(IverticesAtEdgeCoarse(1,iedge))+&
          DuCoarse(IverticesAtEdgeCoarse(2,iedge)))
    
    END DO

    ! Loop over the elements
    DO iel = 1, NELcoarse
      ! Calculate the quad cell midpoint DOF
      DuFine(NVTcoarse + NMTcoarse + iel) = 0.25_DP * (&
          DuCoarse(IverticesAtElementCoarse(1,iel))+&
          DuCoarse(IverticesAtElementCoarse(2,iel))+&
          DuCoarse(IverticesAtElementCoarse(3,iel))+&
          DuCoarse(IverticesAtElementCoarse(4,iel)))
    END DO

!  This is the 'old' implemention, based on Feat 1.x
!  ! local variables
!  REAL(DP), PARAMETER :: Q2 = .5_DP
!  REAL(DP), PARAMETER :: Q4 = .25_DP
!  
!  INTEGER(PREC_ELEMENTIDX) :: iel,ielh1,ielh2,ielh3,ielh4
!  REAL(DP) :: duh1,duh2,duh3,duh4
!
!    ! Copy the first NVT entries - they belong to the coarse grid vertices
!    ! that are fine grid vertices at the same time.
!    CALL lalg_copyVectorDble (DuCoarse,DuFine(1:SIZE(DuCoarse)))
!
!    ! Loop over the elements
!    DO iel=1,NELCoarse
!
!      duh1=DuCoarse(IverticesAtElementCoarse(1,iel))
!      duh2=DuCoarse(IverticesAtElementCoarse(2,iel))
!      duh3=DuCoarse(IverticesAtElementCoarse(3,iel))
!      duh4=DuCoarse(IverticesAtElementCoarse(4,iel))
!
!      ielh1=iel
!      ielh2=IneighboursAtElementFine(2,ielh1)
!      ielh3=IneighboursAtElementFine(2,ielh2)
!      ielh4=IneighboursAtElementFine(2,ielh3)
!
!      ! Now check on every of the edges, if we already computed
!      ! the value in the midpoint: Compute only if the neighbour
!      ! element has smaller number.
!      IF (IneighboursAtElementCoarse(1,iel) .LT. iel) &
!        DuFine(IverticesAtElementFine(2,ielh1)) = Q2*(duh1+duh2)
!
!      IF (IneighboursAtElementCoarse(2,iel) .LT. iel) &
!        DuFine(IverticesAtElementFine(2,ielh2)) = Q2*(duh2+duh3)
!
!      IF (IneighboursAtElementCoarse(3,iel) .LT. iel) &
!        DuFine(IverticesAtElementFine(2,ielh3)) = Q2*(duh3+duh4)
!
!      IF (IneighboursAtElementCoarse(4,iel) .LT. iel) &
!        DuFine(IverticesAtElementFine(2,ielh4)) = Q2*(duh4+duh1)
!        
!      ! Don't forget the DOF in the midpoint of the element
!      DuFine(IverticesAtElementFine(3,iel)) = Q4*(duh1+duh2+duh3+duh4)
!
!    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformQ1_double (DuCoarse,DuFine, &
             IverticesAtEdgeCoarse, IverticesAtElementCoarse, &
             NVTcoarse, NMTcoarse, NELcoarse)
!               IverticesAtElementFine,IneighboursAtElementFine,&
!               NELfine)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IverticesAtEdge array on the coarse grid.
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdgeCoarse
  
  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse

  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse

  ! Number of edges in the coarse grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTcoarse
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse

! 'old' parameters
!  ! IverticesAtElement array (KVERT) on the fine grid
!  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
!  
!  ! IneighboursAtElement array on the coarse grid
!  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
!  
!  ! Number of elements in the fine grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELfine
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_EDGEIDX) :: iedge
  INTEGER(PREC_VERTEXIDX) :: ivt
  REAL(DP) :: dx
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    CALL lalg_copyVectorDble (DuFine(1:NVTcoarse),DuCoarse)
    
    ! Loop over the edges
    DO iedge = 1, NMTcoarse
      ! get the fine grid DOF
      dx = 0.5_DP * DuFine(NVTcoarse + iedge)

      ! distribute it to the coarse grid DOFs
      ivt = IverticesAtEdgeCoarse(1,iedge)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtEdgeCoarse(2,iedge)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
    END DO
    
    ! Loop over the elements
    DO iel = 1, NELcoarse
      ! get the fine grid DOF
      dx = 0.25_DP * DuFine(NVTcoarse + NMTcoarse + iel)
      
      ! distribute it to the coarse grid DOFs
      ivt = IverticesAtElementCoarse(1, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(2, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(3, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(4, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
    END DO
    
! This is the 'old' implementation, based on Feat 1.x
!  ! local variables
!  REAL(DP), PARAMETER :: Q2 = .5_DP
!  REAL(DP), PARAMETER :: Q4 = .25_DP
!  
!  INTEGER(PREC_ELEMENTIDX) :: iel
!  INTEGER(PREC_VERTEXIDX) :: i1,i2,i3,i4
!  
!    ! The information that was 'distributed' in the prolongation has to
!    ! be 'collected'.
!    !
!    ! Copy the first NVT entries - this gives the first additive contribution.
!    CALL lalg_copyVectorDble (DuFine(1:SIZE(DuCoarse)),DuCoarse)
!    
!    ! Loop over the elements to collect the missing additive contributions:
!    DO iel=1,NELfine
!      i1=IverticesAtElementFine(1,iel)
!      i2=IverticesAtElementFine(2,iel)
!      i3=IverticesAtElementFine(3,iel)
!      i4=IverticesAtElementFine(4,iel)
!
!      ! Additive contribution of the midpoint
!      DuCoarse(i1) = DuCoarse(i1)+Q4*(DuFine(i2)+DuFine(i3)+DuFine(i4))
!
!      ! Additional contribution on the boundary:
!      IF (IneighboursAtElementFine(1,iel) .EQ. 0) &
!        DuCoarse(i1) = DuCoarse(i1)+Q4*DuFine(i2)
!      IF (IneighboursAtElementFine(4,iel) .EQ. 0) &
!        DuCoarse(i1) = DuCoarse(i1)+Q4*DuFine(i4)
!    END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUnifQ1FMzero_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse,IverticesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse,NELfine)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $Q_1$, uniform triangulation, double precision vector.
  ! Experimental FEAST MIRROR variant, zeroes the entries on the boundary.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
  
  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse

  ! Number of elements in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELfine
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
  INTEGER(PREC_POINTIDX) :: i1,i2,i3,i4
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

    ! DOF's on the boundary get value 0.0.
    DO iel=1,NELfine
      i1=IverticesAtElementFine(1,iel)
      i2=IverticesAtElementFine(2,iel)
      i3=IverticesAtElementFine(3,iel)
      i4=IverticesAtElementFine(4,iel)

      ! Additional contribution on the boundary:
      IF (IneighboursAtElementFine(1,iel) .EQ. 0) THEN
        DuFine(i1) = 0.0_DP
        DuFine(i2) = 0.0_DP
      END IF
      IF (IneighboursAtElementFine(2,iel) .EQ. 0) THEN
        DuFine(i2) = 0.0_DP
        DuFine(i3) = 0.0_DP
      END IF
      IF (IneighboursAtElementFine(3,iel) .EQ. 0) THEN
        DuFine(i3) = 0.0_DP
        DuFine(i4) = 0.0_DP
      END IF
      IF (IneighboursAtElementFine(4,iel) .EQ. 0) THEN
        DuFine(i4) = 0.0_DP
        DuFine(i1) = 0.0_DP
      END IF
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformQ1FM_double (DuCoarse,DuFine, &
               IverticesAtElementFine,IneighboursAtElementFine,&
               NELfine)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
  ! Cellwise approach for FEAST mirror boundary.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
  
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
  INTEGER(PREC_VERTEXIDX) :: i1,i2,i3,i4
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    CALL lalg_clearVectorDble (DuCoarse)
    
    ! Loop over the elements to collect the missing additive contributions:
    DO iel=1,NELfine
      i1=IverticesAtElementFine(1,iel)
      i2=IverticesAtElementFine(2,iel)
      i3=IverticesAtElementFine(3,iel)
      i4=IverticesAtElementFine(4,iel)

      ! Additive contribution of all vertices to the coarse grid vertex
      DuCoarse(i1) = DuCoarse(i1)+Q4*(DUfine(i1)+DuFine(i2)+DuFine(i3)+DuFine(i4))

      ! Additional contribution on the boundary:
      IF (IneighboursAtElementFine(1,iel) .EQ. 0) THEN
        DuCoarse(i1) = DuCoarse(i1)+2.0_DP*Q4*DuFine(i2)
        DuCoarse(i1) = DuCoarse(i1)+2.0_DP*Q4*DuFine(i1)
      END IF
      IF (IneighboursAtElementFine(4,iel) .EQ. 0) THEN
        DuCoarse(i1) = DuCoarse(i1)+2.0_DP*Q4*DuFine(i4)
        DuCoarse(i1) = DuCoarse(i1)+2.0_DP*Q4*DuFine(i1)
      END IF
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUnifQ1FMzero_double (DuCoarse,DuFine,&
               IverticesAtElementFine,IneighboursAtElementFine,&
               IverticesAtElementCoarse,IneighboursAtElementCoarse,&
               NELfine,NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
  ! Cellwise approach for FEAST mirror boundary.
  ! The DOF's on the boundary are set to 0.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
  
  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse
  
  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! Number of elements in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
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
  INTEGER(PREC_VERTEXIDX) :: i1,i2,i3,i4
  
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

      ! Additional contribution on the boundary:
      IF (IneighboursAtElementFine(1,iel) .EQ. 0) &
        DuCoarse(i1) = DuCoarse(i1)+Q4*DuFine(i2)
      IF (IneighboursAtElementFine(4,iel) .EQ. 0) &
        DuCoarse(i1) = DuCoarse(i1)+Q4*DuFine(i4)
    END DO
    
    ! DOF's on the boundary get value 0.0.
    DO iel=1,NELcoarse
      i1=IverticesAtElementCoarse(1,iel)
      i2=IverticesAtElementCoarse(2,iel)
      i3=IverticesAtElementCoarse(3,iel)
      i4=IverticesAtElementCoarse(4,iel)

      ! Additional contribution on the boundary:
      IF (IneighboursAtElementCoarse(1,iel) .EQ. 0) THEN
        DuCoarse(i1) = 0.0_DP
        DuCoarse(i2) = 0.0_DP
      END IF
      IF (IneighboursAtElementCoarse(2,iel) .EQ. 0) THEN
        DuCoarse(i2) = 0.0_DP
        DuCoarse(i3) = 0.0_DP
      END IF
      IF (IneighboursAtElementCoarse(3,iel) .EQ. 0) THEN
        DuCoarse(i3) = 0.0_DP
        DuCoarse(i4) = 0.0_DP
      END IF
      IF (IneighboursAtElementCoarse(4,iel) .EQ. 0) THEN
        DuCoarse(i4) = 0.0_DP
        DuCoarse(i1) = 0.0_DP
      END IF
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
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse
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

  SUBROUTINE mlprj_prolUniformQ2_double (DuCoarse,DuFine, &
               IverticesAtElementFine,IedgesAtElementFine,&
               IneighboursAtElementFine,&
               NVTfine, NMTfine, NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $Q_2$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IverticesAtElement array (KVERT) on the fine grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementFine

  ! IedgesAtElement array (KMID) on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine

  ! Number of vertices on the fine grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTfine

  ! Number of edges in the fine grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER :: i
  INTEGER(PREC_ELEMENTIDX), DIMENSION(4) :: IelFine

    ! Copy the first NVT+NMT+NEL entries - they belong to the coarse grid 
    ! vertices/edge midpoints/element midpoints and
    ! are fine grid vertices at the same time.
    CALL lalg_copyVectorDble (DuCoarse,DuFine(1:SIZE(DuCoarse)))

    ! Loop over the elements of the coarse grid
    DO iel=1,NELcoarse
   
      ! Obtain the numbers of the fine grid elements inside of the
      ! coarse grid element. According to the regular refinement, the element
      ! number of the coarse grid element is the element number of the
      ! first element inside the coarse grid element on the fine grid.
      IelFine(1)=iel
      IelFine(2)=IneighboursAtElementFine(2,IelFine(1))
      IelFine(3)=IneighboursAtElementFine(2,IelFine(2))
      IelFine(4)=IneighboursAtElementFine(2,IelFine(3))

      ! Loop over the fine grid elements in the coarse grid element.
      ! 'Distribute' the information from the edge midpoints and the element
      ! midpoint to the edge midpoints/element midpoints of the
      ! fine grid element.      
      DO i=1,4
        ! Distribute information on the edges of the coarse grid element
        ! to the edges of the fine grid element i inside of the coarse
        ! grid element.
        DUfine(IedgesAtElementFine(1,IelFine(i)))= &
             +(3.0/8.0)*DUfine(IverticesAtElementFine(1,IelFine(i))) &
             +(3.0/4.0)*DUfine(IverticesAtElementFine(2,IelFine(i))) &
             -(1.0/8.0)*DUfine(IverticesAtElementFine(1,IelFine(MOD(i,4)+1)))
        DUfine(IedgesAtElementFine(4,IelFine(i)))= &
             +(3.0/8.0)*DUfine(IverticesAtElementFine(1,IelFine(i))) &
             +(3.0/4.0)*DUfine(IverticesAtElementFine(4,IelFine(i))) &
             -(1.0/8.0)*DUfine(IverticesAtElementFine(1,IelFine(MOD(i+2,4)+1)))
        DUfine(IedgesAtElementFine(2,IelFine(i)))= &
             +(3.0/8.0)*DUfine(IverticesAtElementFine(2,IelFine(i))) &
             +(3.0/4.0)*DUfine(IverticesAtElementFine(3,IelFine(i))) &
             -(1.0/8.0)*DUfine(IverticesAtElementFine(4,IelFine(MOD(i+2,4)+1)))
             
        ! Distribute information of the coarse grid midpoint to
        ! the fine grid midpoint of the fine grid element i inside
        ! of the coarse grid element.
        DUfine(NVTfine+NMTfine+IelFine(i))= &
             +(9.0/64.0)*DUfine(IverticesAtElementFine(1,IelFine(i))) &
             +(18.0/64.0)*DUfine(IverticesAtElementFine(2,IelFine(i))) &
             +(36.0/64.0)*DUfine(IverticesAtElementFine(3,IelFine(i))) &
             +(18.0/64.0)*DUfine(IverticesAtElementFine(4,IelFine(i))) &
             -(3.0/64.0)*DUfine(IverticesAtElementFine(1,IelFine(MOD(i,4)+1))) &
             -(6.0/64.0)*DUfine(IverticesAtElementFine(2,IelFine(MOD(i,4)+1))) &
             -(3.0/64.0)*DUfine(IverticesAtElementFine(1,IelFine(MOD(i+2,4)+1))) &
             -(6.0/64.0)*DUfine(IverticesAtElementFine(4,IelFine(MOD(i+2,4)+1))) &
             +(1.0/64.0)*DUfine(IverticesAtElementFine(1,IelFine(MOD(i+1,4)+1)))
      ENDDO
    ENDDO

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformQ2_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse, IedgesAtElementCoarse, &
               IedgesAtElementFine, IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NMTcoarse,NMTfine,NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $Q_2$, uniform triangulation, double precision vector.
!</description>

!<input>  
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine

  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse
  
  ! IedgesAtElement array (KMID) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementCoarse
  
  ! IedgesAtElement array (KMID) on the fine grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine
  
  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse

  ! Number of vertices in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine

  ! Number of edges in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NMTcoarse

  ! Number of elements in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NMTfine

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
  INTEGER :: i
  INTEGER(PREC_ELEMENTIDX), DIMENSION(4) :: IelFine

    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    !
    ! Copy the first NVT+NMT+NEL (coarse) entries - this gives the first 
    ! additive contribution: The values in the corners/edge midpoints/
    ! element midpoints of the coarse grid stem from with the values of the
    ! corners of the fine grid.
    CALL lalg_copyVectorDble (DuFine(1:SIZE(DuCoarse)),DuCoarse)

    ! Loop over the elements to collect the missing additive contributions:
    DO iel=1,NELcoarse
    
      ! Obtain the numbers of the fine grid elements inside of the
      ! coarse grid element. According to the regular refinement, the element
      ! number of the coarse grid element is the element number of the
      ! first element inside the coarse grid element on the fine grid.
      IelFine(1)=iel
      IelFine(2)=IneighboursAtElementFine(2,IelFine(1))
      IelFine(3)=IneighboursAtElementFine(2,IelFine(2))
      IelFine(4)=IneighboursAtElementFine(2,IelFine(3))

      DO i=1,4
      
        ! Collect information from the corners of the fine grid element
        DUcoarse(IverticesAtElementCoarse(i,iel)) = &
             DUcoarse(IverticesAtElementCoarse(i,iel)) &
             +0.5*(24.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(i))) &
             +0.5*(24.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(i))) &
             -0.5*(8.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(MOD(i,4)+1))) &
             -0.5*(8.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(MOD(i+2,4)+1))) &
             +(9.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(i)) &
             -(3.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(MOD(i,4)+1)) &
             +(1.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(MOD(i+1,4)+1)) &
             -(3.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(MOD(i+2,4)+1))
             
        IF(IneighboursAtElementFine(1,IelFine(i)).EQ.0) THEN
          DUcoarse(IverticesAtElementCoarse(i,iel))= &
              DUcoarse(IverticesAtElementCoarse(i,iel)) &
              +0.5*(24.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(i)))
        END IF
        
        IF(IneighboursAtElementFine(4,IelFine(i)).EQ.0) THEN
          DUcoarse(IverticesAtElementCoarse(i,iel))= &
              DUcoarse(IverticesAtElementCoarse(i,iel)) &
              +0.5*(24.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(i)))
        END IF
        
        IF(IneighboursAtElementFine(4,IelFine(MOD(i,4)+1)).EQ.0) THEN
          DUcoarse(IverticesAtElementCoarse(i,iel))= &
              DUcoarse(IverticesAtElementCoarse(i,iel)) &
              -0.5*(8.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(MOD(i,4)+1)))
        END IF
        
        IF(IneighboursAtElementFine(1,IelFine(MOD(i+2,4)+1)).EQ.0) THEN
          DUcoarse(IverticesAtElementCoarse(i,iel))= &
              DUcoarse(IverticesAtElementCoarse(i,iel)) &
              -0.5*(8.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(MOD(i+2,4)+1)))
        END IF

        ! Collect information from the edge midpoints of the fine grid element
        DUcoarse(IEdgesAtElementCoarse(i,iel))= &
            DUcoarse(IEdgesAtElementCoarse(i,iel)) &
             +0.5*(48.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(i))) &
             +(24.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(i))) &
             +0.5*(48.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(MOD(i,4)+1))) &
             -(8.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(MOD(i+1,4)+1))) &
             +(18.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(i)) &
             +(18.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(MOD(i,4)+1)) &
             -(6.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(MOD(i+1,4)+1)) &
             -(6.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(MOD(i+2,4)+1))
             
        IF(IneighboursAtElementFine(1,IelFine(i)).EQ.0) THEN
          DUcoarse(IEdgesAtElementCoarse(i,iel))= &
              DUcoarse(IEdgesAtElementCoarse(i,iel)) &
              +0.5*(48.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(i)))
        END IF
        
        IF(IneighboursAtElementFine(4,IelFine(MOD(i,4)+1)).EQ.0) THEN
          DUcoarse(IEdgesAtElementCoarse(i,iel))= &
              DUcoarse(IEdgesAtElementCoarse(i,iel)) &
              +0.5*(48.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(MOD(i,4)+1)))
        END IF

      ENDDO
      
      ! Collect information from the midpoints of the fine grid elements
      DUcoarse(NVTcoarse+NMTCoarse+iel)= &
          DUcoarse(NVTcoarse+NMTCoarse+iel) &
           +(48.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(1))) &
           +(48.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(2))) &
           +(48.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(3))) &
           +(48.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(4))) &
           +(36.0/64.0)*DUfine(NVTfine+NMTfine+IelFine(1)) &
           +(36.0/64.0)*DUfine(NVTfine+NMTfine+IelFine(2)) &
           +(36.0/64.0)*DUfine(NVTfine+NMTfine+IelFine(3)) &
           +(36.0/64.0)*DUfine(NVTfine+NMTfine+IelFine(4))
       
    ENDDO

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
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse

  ! Number of edges in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NMTcoarse

  ! Number of elements in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NELcoarse
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
  ! Support for QP1 element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformQP1_double (DuCoarse,DuFine, &
               IneighboursAtElementFine,NELcoarse,NELfine)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! QP1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IneighboursAtElement array on the fine grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse

  ! Number of elements in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELfine
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>

  ! In the QP1 element, the DOF's in the coarse and fine grid are 
  ! organised as follows:
  !
  ! +-----------+-----------+
  ! |           |           |
  ! |           |           |
  ! |     X4->      <-X1    |
  ! |     |     ^     |     |
  ! |     v     |     v     |
  ! +--------   O->  -------+
  ! |     ^           ^     |
  ! |     |     |     |     |
  ! |     X1->  |   <-X1    |
  ! |           |           |
  ! |           |           |
  ! +-----------+-----------+
  !
  ! The function value of "O" must be transported by a linear mapping to
  ! all the X. The derivative in "O" is mapped to all "X" as it's constant
  ! in the whole element. This gives the following prolongation matrix
  ! (where O=(f,u,v), Xi = (fi,ui,vi) )
  !
  !   ( 1 -.5 -.5 )  *  ( f )  =  ( f1 )
  !   ( 1  .5 -.5 )     ( u )     ( f2 )
  !   ( 1  .5  .5 )     ( v )     ( f3 )
  !   ( 1 -.5  .5 )               ( f4 )
  !   (     1   0 )               ( u1 )
  !   (     0  -1 )               ( u2 )
  !   (    -1   0 )               ( u3 )
  !   (     1   1 )               ( u4 )
  !   (     0   1 )               ( v1 )
  !   (     1   0 )               ( v2 )
  !   (     0  -1 )               ( v3 )
  !   (    -1   0 )               ( v4 )
  !
  ! The restriction matrix is the transposed of that...
  
  ! local variables
  REAL(DP), PARAMETER :: Q2 = .5_DP
  
  INTEGER(PREC_ELEMENTIDX) :: iel,ielh1,ielh2,ielh3,ielh4
  REAL(DP) :: duh1,duh2,duh3

    ! Loop over the elements
    DO iel=1,NELCoarse

      ! Get the values of the basis functions of the coarse grid element
      duh1=DuCoarse(iel)
      duh2=DuCoarse(iel+NELcoarse)
      duh3=DuCoarse(iel+2*NELcoarse)

      ! Get fine grid element numbers
      ielh1=iel
      ielh2=IneighboursAtElementFine(2,ielh1)
      ielh3=IneighboursAtElementFine(2,ielh2)
      ielh4=IneighboursAtElementFine(2,ielh3)

      ! Apply the prolonfation matrix to the coarse grid basis functions
      ! to get the fine grid values.
      DuFine(ielh1) = duh1 - Q2*duh2 - Q2*duh3
      DuFine(ielh2) = duh1 + Q2*duh2 - Q2*duh3
      DuFine(ielh3) = duh1 + Q2*duh2 + Q2*duh3
      DuFine(ielh4) = duh1 - Q2*duh2 + Q2*duh3
      
      DuFine(NELfine+ielh1) = duh2
      DuFine(NELfine+ielh2) = -duh3
      DuFine(NELfine+ielh3) = -duh2
      DuFine(NELfine+ielh4) = duh3

      DuFine(2*NELfine+ielh1) = duh3
      DuFine(2*NELfine+ielh2) = duh2
      DuFine(2*NELfine+ielh3) = -duh3
      DuFine(2*NELfine+ielh4) = -duh2

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformQP1_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse,NELfine)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! QP1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
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
  INTEGER(PREC_ELEMENTIDX) :: ielh1,ielh2,ielh3,ielh4
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'. This means, we apply the transposed prolongation
    ! matrix to the RHS vector.
    
    ! Loop over the elements to collect the additive contributions:
    DO iel=1,NELcoarse
    
      ! Get fine grid element numbers
      ielh1=iel
      ielh2=IneighboursAtElementFine(2,ielh1)
      ielh3=IneighboursAtElementFine(2,ielh2)
      ielh4=IneighboursAtElementFine(2,ielh3)

      ! Collect the distributed values to form the coarse grid RHS:
      DuCoarse(iel) = DuFine(ielh1) + DuFine(ielh2) &
                  + DuFine(ielh3) + DuFine(ielh4)
                  
      DuCoarse(NELcoarse+iel) = &
          -Q2*DuFine(ielh1) + Q2*DuFine(ielh2)      &
            + Q2*DuFine(ielh3) - Q2*DuFine(ielh4)   &
          + DuFine(NELfine+ielh1) - DuFine(NELfine+ielh3)    &
          + DuFine(2*NELfine+ielh2) - DuFine(2*NELfine+ielh4)

      DuCoarse(2*NELcoarse+iel) = &
          -Q2*DuFine(ielh1) - Q2*DuFine(ielh2)      &
            + Q2*DuFine(ielh3) + Q2*DuFine(ielh4)   &
          - DuFine(NELfine+ielh2) + DuFine(NELfine+ielh4)    &
          + DuFine(2*NELfine+ielh1) - DuFine(2*NELfine+ielh3)
          
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformQP1_double (DuCoarse,DuFine, &
                  IneighboursAtElementFine, NELcoarse,NELfine)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! QP1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
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
  INTEGER(PREC_ELEMENTIDX) :: ielh1,ielh2,ielh3,ielh4

    ! Loop over the elements 
    DO iel=1,NELcoarse
    
      ! Get fine grid element numbers
      ielh1=iel
      ielh2=IneighboursAtElementFine(2,ielh1)
      ielh3=IneighboursAtElementFine(2,ielh2)
      ielh4=IneighboursAtElementFine(2,ielh3)
      
      ! Interpolate the solution on the fine grid to the coarse grid;
      ! take the mean of 4 values each. Take care of the orientation
      ! of the basis functions!

      DuCoarse(iel) = Q4*(DuFine(ielh1) + DuFine(ielh2) &
                         +DuFine(ielh3) + DuFine(ielh4) )
                         
      DuCoarse(NELcoarse+iel) = &
        Q4*(DuFine(ielh1+NELfine) - DuFine(ielh2+2*NELfine) &
           -DuFine(ielh3+NELfine) + DuFine(ielh4+2*NELfine) )

      DuCoarse(2*NELcoarse+iel) = &
        Q4*(DuFine(ielh1+2*NELfine) + DuFine(ielh2+NELfine) &
           -DuFine(ielh3+2*NELfine) - DuFine(ielh4+NELfine) )
          
    END DO
    
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
  
  ! Weights for the restriction; all coefficients are halfed, so dividing
  ! by 2 is not necessary in the calculation routines.
  REAL(DP), PARAMETER :: A1=0.5_DP, A2=-0.125_DP, A3=0.0_DP, A4=0.125_DP
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
               DvertexCoordsCoarse,IverticesAtElementCoarse,DelementAreaCoarse,&
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

  ! DvertexCoords array on the coarse grid
  REAL(DP), DIMENSION(:,:), INTENT(IN)                :: DvertexCoordsCoarse

  ! IverticesAtElement array on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse
  
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
      dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      IF (daspectRatio(0) .LT. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)
      
      ! and the area of that element.
      darea(0) = DelementAreaCoarse(iel)
      
      ! Then the remaining neighbours.
      DO i=1,TRIA_MAXNME2D
        IF (IELA(i) .NE. 0) THEN
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(i)))
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
        IF (daspectRatio(0) .GT. daspectRatioBound) idoConstant = 2
        
        ! and if iarIndicator>2, also check the neighbour element
        IF (iarIndicator .GE. 2) THEN
!          DO i=1,TRIA_MAXNME2D
!            IF (daspectRatio(i) .GT. daspectRatioBound) idoConstant(i) = 2
!          END DO
          WHERE (daspectRatio(1:4) .GT. daspectRatioBound) idoConstant = 2
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
    
    ! Weights for the restriction
    REAL(DP), PARAMETER :: A1=1.0_DP, A2=-0.125_DP, A3=0.0_DP, A4=0.125_DP
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
        DuCoarse(IM1)=     A1*(DUH1+DUH2)+2.0_DP*A2*(DUH4+DUH7) &
                +2.0_DP*A3*(DUH5+DUH6)+2.0_DP*A4*(DUH3+DUH8) &
                +       A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
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
        DuCoarse(IM2)= A1*(DUH3+DUH4)+2.0_DP*A2*(DUH6+DUH1) &
                +2.0_DP*A3*(DUH7+DUH8)+2.0_DP*A4*(DUH5+DUH2) &
                +       A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
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
        DuCoarse(IM3)= A1*(DUH5+DUH6)+2.0_DP*A2*(DUH8+DUH3) &
                +2.0_DP*A3*(DUH1+DUH2)+2.0_DP*A4*(DUH7+DUH4) &
                +       A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
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
        DuCoarse(IM4)= A1*(DUH7+DUH8)+2.0_DP*A2*(DUH2+DUH5) &
                +2.0_DP*A3*(DUH3+DUH4)+2.0_DP*A4*(DUH1+DUH6) &
                +       A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
      ENDIF

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformEx30ext_double (DuCoarse,DuFine, &
               DvertexCoordsCoarse,IverticesAtElementCoarse,DelementAreaCoarse,&
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
  
  ! DvertexCoords array on the coarse grid
  REAL(DP), DIMENSION(:,:), INTENT(IN)                :: DvertexCoordsCoarse

  ! IverticesAtElement array on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse
  
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
      dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      IF (daspectRatio(0) .LT. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)
      
      ! and the area of that element.
      darea(0) = DelementAreaCoarse(iel)
      
      ! Then the remaining neighbours.
      DO i=1,TRIA_MAXNME2D
        IF (IELA(i) .NE. 0) THEN
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(i)))
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
        IF (daspectRatio(0) .GT. daspectRatioBound) idoConstant = 2
        
        ! and if iarIndicator>2, also check the neighbour element
        IF (iarIndicator .GE. 2) THEN
!          DO i=1,TRIA_MAXNME2D
!            IF (daspectRatio(i) .GT. daspectRatioBound) idoConstant(i) = 2
!          END DO
          WHERE (daspectRatio(1:4) .GT. daspectRatioBound) idoConstant = 2
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
    
    ! Weights for the restriction
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

  ! ***************************************************************************
  ! Support for Q1~ element, DOF's = edge midpoints
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformEx31_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E031/EM31, uniform triangulation, double precision vector.
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
  
  ! Weights for the restriction; all coefficients are halfed, so dividing
  ! by 2 is not necessary in the calculation routines.
  REAL(DP), PARAMETER :: A1=0.46875_DP, A2=-0.09375_DP, A3=-0.03125_DP, A4=0.15625_DP
  REAL(DP), PARAMETER :: A5=0.5625_DP, A6=0.1875_DP, A7=0.0625_DP, A8=0.1875_DP
  
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

  SUBROUTINE mlprj_prolUniformEx31ext_double (DuCoarse,DuFine, &
               DvertexCoordsCoarse,IverticesAtElementCoarse,DelementAreaCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NELcoarse, &
               iweightingType, daspectRatioBound, iarIndicator)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E031/EM31, uniform triangulation, double precision vector.
  !
  ! Extended version. Switch to constant prolongation if aspect ratio
  ! of an element is too large.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse

  ! DvertexCoords array on the coarse grid
  REAL(DP), DIMENSION(:,:), INTENT(IN)                :: DvertexCoordsCoarse

  ! IverticesAtElement array on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse
  
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
      RESHAPE((/0.9375_DP, -0.1875_DP, -0.0625_DP, 0.3125_DP, &
                0.5625_DP, 0.1875_DP, 0.0625_DP, 0.1875_DP, &
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
      dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      IF (daspectRatio(0) .LT. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)
      
      ! and the area of that element.
      darea(0) = DelementAreaCoarse(iel)
      
      ! Then the remaining neighbours.
      DO i=1,TRIA_MAXNME2D
        IF (IELA(i) .NE. 0) THEN
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(i)))
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
        IF (daspectRatio(0) .GT. daspectRatioBound) idoConstant = 2
        
        ! and if iarIndicator>2, also check the neighbour element
        IF (iarIndicator .GE. 2) THEN
!          DO i=1,TRIA_MAXNME2D
!            IF (daspectRatio(i) .GT. daspectRatioBound) idoConstant(i) = 2
!          END DO
          WHERE (daspectRatio(1:4) .GT. daspectRatioBound) idoConstant = 2
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

  SUBROUTINE mlprj_restUniformEx31_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! E031/EM31 element, uniform triangulation, double precision vector.
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
    
    ! Weights for the restriction
    REAL(DP), PARAMETER :: A1=0.9375_DP, A2=-0.09375_DP, A3=-0.03125_DP, A4=0.15625_DP
    REAL(DP), PARAMETER :: A5=0.5625_DP, A6=0.1875_DP, A7=0.0625_DP, A8=0.1875_DP

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
        DuCoarse(IM1)=     A1*(DUH1+DUH2)+2.0_DP*A2*(DUH4+DUH7) &
                +2.0_DP*A3*(DUH5+DUH6)+2.0_DP*A4*(DUH3+DUH8) &
                +       A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
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
        DuCoarse(IM2)= A1*(DUH3+DUH4)+2.0_DP*A2*(DUH6+DUH1) &
                +2.0_DP*A3*(DUH7+DUH8)+2.0_DP*A4*(DUH5+DUH2) &
                +       A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
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
        DuCoarse(IM3)= A1*(DUH5+DUH6)+2.0_DP*A2*(DUH8+DUH3) &
                +2.0_DP*A3*(DUH1+DUH2)+2.0_DP*A4*(DUH7+DUH4) &
                +       A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
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
        DuCoarse(IM4)= A1*(DUH7+DUH8)+2.0_DP*A2*(DUH2+DUH5) &
                +2.0_DP*A3*(DUH3+DUH4)+2.0_DP*A4*(DUH1+DUH6) &
                +       A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
      ENDIF

    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformEx31ext_double (DuCoarse,DuFine, &
               DvertexCoordsCoarse,IverticesAtElementCoarse,DelementAreaCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NELcoarse, &
               iweightingType, daspectRatioBound, iarIndicator)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! E031/EM31 element, uniform triangulation, double precision vector.
  !
  ! Extended version. Switch to constant restriction if aspect ratio
  ! of an element is too large.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! DvertexCoords array on the coarse grid
  REAL(DP), DIMENSION(:,:), INTENT(IN)                :: DvertexCoordsCoarse

  ! IverticesAtElement array on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse
  
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
        RESHAPE((/0.9375_DP, -0.1875_DP, -0.0625_DP, 0.3125_DP, &
                  0.5625_DP, 0.1875_DP, 0.0625_DP, 0.1875_DP, &
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
      dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      IF (daspectRatio(0) .LT. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)
      
      ! and the area of that element.
      darea(0) = DelementAreaCoarse(iel)
      
      ! Then the remaining neighbours.
      DO i=1,TRIA_MAXNME2D
        IF (IELA(i) .NE. 0) THEN
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(i)))
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
        IF (daspectRatio(0) .GT. daspectRatioBound) idoConstant = 2
        
        ! and if iarIndicator>2, also check the neighbour element
        IF (iarIndicator .GE. 2) THEN
!          DO i=1,TRIA_MAXNME2D
!            IF (daspectRatio(i) .GT. daspectRatioBound) idoConstant(i) = 2
!          END DO
          WHERE (daspectRatio(1:4) .GT. daspectRatioBound) idoConstant = 2
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
  ! Support for 3D Q0 element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformQ0_3D_double (DuCoarse,DuFine,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $Q_0$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: ielf
  REAL(DP) :: duh

    ! Loop over the elements
    DO iel = 1, NELCoarse

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = NELcoarse + 7*(iel-1) + 1
      ielf(3) = ielf(2) + 1
      ielf(4) = ielf(3) + 1
      ielf(5) = ielf(4) + 1
      ielf(6) = ielf(5) + 1
      ielf(7) = ielf(6) + 1
      ielf(8) = ielf(7) + 1

      ! Put the value on the coarse grid into all four child
      ! elements
      duh = DuCoarse(iel)
      DuFine(ielf(1:8)) = duh
    END DO

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformQ0_3D_double (DuCoarse,DuFine,NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $Q_0$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
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
  INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: ielf
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    
    ! Loop over the elements to collect the missing additive contributions:
    DO iel=1,NELcoarse
    
      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = NELcoarse + 7*(iel-1) + 1
      ielf(3) = ielf(2) + 1
      ielf(4) = ielf(3) + 1
      ielf(5) = ielf(4) + 1
      ielf(6) = ielf(5) + 1
      ielf(7) = ielf(6) + 1
      ielf(8) = ielf(7) + 1
      
      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      DuCoarse(iel)= DuFine(ielf(1))+DuFine(ielf(2))+&
                     DuFine(ielf(3))+DuFine(ielf(4))+DuFine(ielf(5))+&
                     DuFine(ielf(6))+DuFine(ielf(7))+DuFine(ielf(8))
      
    END DO
    
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformQ0_3D_double (DuCoarse,DuFine,NELcoarse)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $Q_0$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
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
  INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: ielf
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    
    ! Loop over the elements to collect the missing additive contributions:
    DO iel=1,NELcoarse
    
      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = NELcoarse + 7*(iel-1) + 1
      ielf(3) = ielf(2) + 1
      ielf(4) = ielf(3) + 1
      ielf(5) = ielf(4) + 1
      ielf(6) = ielf(5) + 1
      ielf(7) = ielf(6) + 1
      ielf(8) = ielf(7) + 1
      
      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      DuCoarse(iel)= 0.125_DP * (DuFine(ielf(1))+DuFine(ielf(2))+&
                     DuFine(ielf(3))+DuFine(ielf(4))+DuFine(ielf(5))+&
                     DuFine(ielf(6))+DuFine(ielf(7))+DuFine(ielf(8)))
      
    END DO
    
  END SUBROUTINE
  
  ! ***************************************************************************
  ! Support for 3D Q1 element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformQ1_3D_double (DuCoarse, DuFine, &
             IverticesAtEdgeCoarse, IverticesAtFaceCoarse, &
             IverticesAtElementCoarse, NVTcoarse, NMTcoarse, NATcoarse,NELcoarse)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $Q_1$, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IverticesAtEdge array on the coarse grid.
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdgeCoarse
  
  ! IverticesAtFace array on the coarse grid.
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtFaceCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse

  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse

  ! Number of edges in the coarse grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTcoarse
  
  ! Number of faces in the coarse grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NATcoarse
  
  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_EDGEIDX) :: iedge, iface

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    CALL lalg_copyVectorDble (DuCoarse,DuFine(1:NVTcoarse))
    
    ! Loop over the edges
    DO iedge = 1, NMTcoarse
      ! Calculate the edge midpoint DOF
      DuFine(NVTcoarse + iedge) = 0.5_DP * (&
          DuCoarse(IverticesAtEdgeCoarse(1,iedge))+&
          DuCoarse(IverticesAtEdgeCoarse(2,iedge)))
    
    END DO
    
    ! Loop over the faces
    DO iface = 1, NATcoarse
      ! Calculate the face midpoint DOF
      DuFine(NVTCoarse + NMTCoarse + iface) = 0.25_DP * (&
          DuCoarse(IverticesAtFaceCoarse(1,iface))+&
          DuCoarse(IverticesAtFaceCoarse(2,iface))+&
          DuCoarse(IverticesAtFaceCoarse(3,iface))+&
          DuCoarse(IverticesAtFaceCoarse(4,iface)))
    END DO

    ! Loop over the elements
    DO iel = 1, NELcoarse
      ! Calculate the hexahedron cell midpoint DOF
      DuFine(NVTcoarse + NMTcoarse + NATcoarse + iel) = 0.125_DP * (&
          DuCoarse(IverticesAtElementCoarse(1,iel))+&
          DuCoarse(IverticesAtElementCoarse(2,iel))+&
          DuCoarse(IverticesAtElementCoarse(3,iel))+&
          DuCoarse(IverticesAtElementCoarse(4,iel))+&
          DuCoarse(IverticesAtElementCoarse(5,iel))+&
          DuCoarse(IverticesAtElementCoarse(6,iel))+&
          DuCoarse(IverticesAtElementCoarse(7,iel))+&
          DuCoarse(IverticesAtElementCoarse(8,iel)))
    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restUniformQ1_3D_double (DuCoarse,DuFine, &
             IverticesAtEdgeCoarse, IverticesAtFaceCoarse, &
             IverticesAtElementCoarse, NVTcoarse, NMTcoarse, NATcoarse,NELcoarse)
  
!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IverticesAtEdge array on the coarse grid.
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtEdgeCoarse
  
  ! IverticesAtFace array on the coarse grid.
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtFaceCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElementCoarse

  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse

  ! Number of edges in the coarse grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTcoarse
  
  ! Number of faces in the coarse grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NATcoarse
  
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
  INTEGER(PREC_EDGEIDX) :: iedge, iface
  INTEGER(PREC_VERTEXIDX) :: ivt
  REAL(DP) :: dx
  
    ! The information that was 'distributed' in the prolongation has to
    ! be 'collected'.
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    CALL lalg_copyVectorDble (DuFine(1:NVTcoarse),DuCoarse)
    
    ! Loop over the edges
    DO iedge = 1, NMTcoarse
      ! get the fine grid DOF
      dx = 0.5_DP * DuFine(NVTcoarse + iedge)

      ! distribute it to the coarse grid DOFs
      ivt = IverticesAtEdgeCoarse(1,iedge)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtEdgeCoarse(2,iedge)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
    END DO
    
    ! Loop over the faces
    DO iface = 1, NATcoarse
      ! get the fine grid DOF
      dx = 0.25_DP * DuFine(NVTcoarse + NMTcoarse + iface)
      
      ! distribute it to the coarse grid DOFs
      ivt = IverticesAtFaceCoarse(1,iface)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtFaceCoarse(2,iface)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtFaceCoarse(3,iface)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtFaceCoarse(4,iface)
    END DO
    
    ! Loop over the elements
    DO iel = 1, NELcoarse
      ! get the fine grid DOF
      dx = 0.125_DP * DuFine(NVTcoarse + NMTcoarse + NATcoarse + iel)
      
      ! distribute it to the coarse grid DOFs
      ivt = IverticesAtElementCoarse(1, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(2, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(3, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(4, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(5, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(6, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(7, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(8, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_interpUniformQ1_3D_double (DuCoarse,DuFine, NVTcoarse)
  
!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVTcoarse
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
  ! Support for 3D Q1~ element
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_prolUniformEx3x_3D_double (DuCoarse,DuFine, &
               IfacesAtElementCoarse,IfacesAtElementFine,&
               IneighboursAtElementCoarse,NVTcoarse,NVTfine,&
               NMTcoarse,NMTfine,NELcoarse,ielType)
  
!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E030/E031/EM30/EM31, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
  
  ! IfacesAtElement array on the coarse grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IfacesAtElementCoarse
  
  ! IfacesAtElement array on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IfacesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! IneighboursAtElement array on the fine grid
  !INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine

  ! Number of vertices in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse

  ! Number of vertices in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine

  ! Number of edges in the coarse grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTcoarse
  
  ! Number of edges in the fine grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
  
  ! Element type to use for prolongation
  INTEGER, INTENT(IN) :: ielType
!</input>
  
!<output>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP) :: dw, A1, A2, A3, A4, A5, A6, A7
  REAL(DP), DIMENSION(6) :: Dx
  INTEGER(PREC_EDGEIDX), DIMENSION(6) :: IM
  INTEGER(PREC_EDGEIDX), DIMENSION(4) :: idx
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: ielf
  
    ! Weights for the restriction; all coefficients are halfed, so dividing
    ! by 2 is not necessary in the calculation routines.
    IF (IAND(ielType, 2**16) .EQ. 0) THEN
      ! Weights for Ex31 element
      A1 =  0.458333333333333_DP   ! = 11/24
      A2 =  0.145833333333333_DP   ! =  7/48
      A3 = -0.104166666666667_DP   ! = -5/48
      A4 = -0.041666666666667_DP   ! = -1/24
      A5 =  0.458333333333333_DP   ! = 11/24
      A6 =  0.083333333333333_DP   ! =  1/12
      A7 = -0.041666666666667_DP   ! = -1/24
    ELSE
      ! Weights for Ex30 element
      A1 = 0.5_DP
      A2 = 0.125_DP
      A3 = -0.125_DP
      A4 = 0.0_DP
      A5 = 0.5_DP
      A6 = 0.0_DP
      A7 = 0.0_DP
    END IF
 
    ! Clear the output vector
    CALL lalg_clearVectorDble(DuFine)
  
    ! Loop over the coarse grid elements
    DO iel = 1, NELcoarse

      ! Get the DOF's of the coarse grid element
      IM(1) = IfacesAtElementCoarse(1,iel) - NVTcoarse - NMTcoarse
      IM(2) = IfacesAtElementCoarse(2,iel) - NVTcoarse - NMTcoarse
      IM(3) = IfacesAtElementCoarse(3,iel) - NVTcoarse - NMTcoarse
      IM(4) = IfacesAtElementCoarse(4,iel) - NVTcoarse - NMTcoarse
      IM(5) = IfacesAtElementCoarse(5,iel) - NVTcoarse - NMTcoarse
      IM(6) = IfacesAtElementCoarse(6,iel) - NVTcoarse - NMTcoarse

      ! Get the values of the corresponding DOF's
      Dx(1) = DuCoarse(IM(1))
      Dx(2) = DuCoarse(IM(2))
      Dx(3) = DuCoarse(IM(3))
      Dx(4) = DuCoarse(IM(4))
      Dx(5) = DuCoarse(IM(5))
      Dx(6) = DuCoarse(IM(6))

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = NELcoarse + 7*(iel-1) + 1
      ielf(3) = ielf(2) + 1
      ielf(4) = ielf(3) + 1
      ielf(5) = ielf(4) + 1
      ielf(6) = ielf(5) + 1
      ielf(7) = ielf(6) + 1
      ielf(8) = ielf(7) + 1

      ! First, we are going to prolongate the DOFs on the boundary of
      ! the coarse grid hexahedron.
      ! Face 1
      ! Calculate scaling factor. If this face of the hexahedron belongs
      ! to the domain's boundary (i.e. the hexahedron does not have a neighbour
      ! at this face), then the scaling factor is 2, otherwise 1.
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(1,iel) .EQ. 0) dw = 2.0_DP
      ! Calculate fine grid DOF indices for this coarse grid face
      idx(1) = IfacesAtElementFine(1, ielf(1)) - NVTfine - NMTfine
      idx(2) = IfacesAtElementFine(1, ielf(2)) - NVTfine - NMTfine
      idx(3) = IfacesAtElementFine(1, ielf(3)) - NVTfine - NMTfine
      idx(4) = IfacesAtElementFine(1, ielf(4)) - NVTfine - NMTfine
      ! distribute the DOF
      DuFine(idx(1))=DuFine(idx(1))+&
         dw*(A1*Dx(1)+A2*Dx(2)+A3*Dx(3)+A3*Dx(4)+A2*Dx(5)+A4*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
         dw*(A1*Dx(1)+A2*Dx(2)+A2*Dx(3)+A3*Dx(4)+A3*Dx(5)+A4*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
         dw*(A1*Dx(1)+A3*Dx(2)+A2*Dx(3)+A2*Dx(4)+A3*Dx(5)+A4*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
         dw*(A1*Dx(1)+A3*Dx(2)+A3*Dx(3)+A2*Dx(4)+A2*Dx(5)+A4*Dx(6))

      ! Face 2
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(2,iel) .EQ. 0) dw = 2.0_DP
      idx(1) = IfacesAtElementFine(2, ielf(1)) - NVTfine - NMTfine
      idx(2) = IfacesAtElementFine(5, ielf(2)) - NVTfine - NMTfine
      idx(3) = IfacesAtElementFine(5, ielf(6)) - NVTfine - NMTfine
      idx(4) = IfacesAtElementFine(2, ielf(5)) - NVTfine - NMTfine
      DuFine(idx(1))=DuFine(idx(1))+&
        dw*(A2*Dx(1)+A1*Dx(2)+A3*Dx(3)+A4*Dx(4)+A2*Dx(5)+A3*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
        dw*(A2*Dx(1)+A1*Dx(2)+A2*Dx(3)+A4*Dx(4)+A3*Dx(5)+A3*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
        dw*(A3*Dx(1)+A1*Dx(2)+A2*Dx(3)+A4*Dx(4)+A3*Dx(5)+A2*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
        dw*(A3*Dx(1)+A1*Dx(2)+A3*Dx(3)+A4*Dx(4)+A2*Dx(5)+A2*Dx(6))

      ! Face 3
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(3,iel) .EQ. 0) dw = 2.0_DP
      idx(1) = IfacesAtElementFine(2, ielf(2)) - NVTfine - NMTfine
      idx(2) = IfacesAtElementFine(5, ielf(3)) - NVTfine - NMTfine
      idx(3) = IfacesAtElementFine(5, ielf(7)) - NVTfine - NMTfine
      idx(4) = IfacesAtElementFine(2, ielf(6)) - NVTfine - NMTfine
      DuFine(idx(1))=DuFine(idx(1))+&
        dw*(A2*Dx(1)+A2*Dx(2)+A1*Dx(3)+A3*Dx(4)+A4*Dx(5)+A3*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
        dw*(A2*Dx(1)+A3*Dx(2)+A1*Dx(3)+A2*Dx(4)+A4*Dx(5)+A3*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
        dw*(A3*Dx(1)+A3*Dx(2)+A1*Dx(3)+A2*Dx(4)+A4*Dx(5)+A2*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
        dw*(A3*Dx(1)+A2*Dx(2)+A1*Dx(3)+A3*Dx(4)+A4*Dx(5)+A2*Dx(6))

      ! Face 4
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(4,iel) .EQ. 0) dw = 2.0_DP
      idx(1) = IfacesAtElementFine(2, ielf(3)) - NVTfine - NMTfine
      idx(2) = IfacesAtElementFine(5, ielf(4)) - NVTfine - NMTfine
      idx(3) = IfacesAtElementFine(5, ielf(8)) - NVTfine - NMTfine
      idx(4) = IfacesAtElementFine(2, ielf(7)) - NVTfine - NMTfine
      DuFine(idx(1))=DuFine(idx(1))+&
        dw*(A2*Dx(1)+A4*Dx(2)+A2*Dx(3)+A1*Dx(4)+A3*Dx(5)+A3*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
        dw*(A2*Dx(1)+A4*Dx(2)+A3*Dx(3)+A1*Dx(4)+A2*Dx(5)+A3*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
        dw*(A3*Dx(1)+A4*Dx(2)+A3*Dx(3)+A1*Dx(4)+A2*Dx(5)+A2*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
        dw*(A3*Dx(1)+A4*Dx(2)+A2*Dx(3)+A1*Dx(4)+A3*Dx(5)+A2*Dx(6))

      ! Face 5
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(5,iel) .EQ. 0) dw = 2.0_DP
      idx(1) = IfacesAtElementFine(2, ielf(4)) - NVTfine - NMTfine
      idx(2) = IfacesAtElementFine(5, ielf(1)) - NVTfine - NMTfine
      idx(3) = IfacesAtElementFine(5, ielf(5)) - NVTfine - NMTfine
      idx(4) = IfacesAtElementFine(2, ielf(8)) - NVTfine - NMTfine
      DuFine(idx(1))=DuFine(idx(1))+&
        dw*(A2*Dx(1)+A3*Dx(2)+A4*Dx(3)+A2*Dx(4)+A1*Dx(5)+A3*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
        dw*(A2*Dx(1)+A2*Dx(2)+A4*Dx(3)+A3*Dx(4)+A1*Dx(5)+A3*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
        dw*(A3*Dx(1)+A2*Dx(2)+A4*Dx(3)+A3*Dx(4)+A1*Dx(5)+A2*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
        dw*(A3*Dx(1)+A3*Dx(2)+A4*Dx(3)+A2*Dx(4)+A1*Dx(5)+A2*Dx(6))

      ! Face 6
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(6,iel) .EQ. 0) dw = 2.0_DP
      idx(1) = IfacesAtElementFine(1, ielf(5)) - NVTfine - NMTfine
      idx(2) = IfacesAtElementFine(1, ielf(6)) - NVTfine - NMTfine
      idx(3) = IfacesAtElementFine(1, ielf(7)) - NVTfine - NMTfine
      idx(4) = IfacesAtElementFine(1, ielf(8)) - NVTfine - NMTfine
      DuFine(idx(1))=DuFine(idx(1))+&
        dw*(A4*Dx(1)+A2*Dx(2)+A3*Dx(3)+A3*Dx(4)+A2*Dx(5)+A1*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
        dw*(A4*Dx(1)+A2*Dx(2)+A2*Dx(3)+A3*Dx(4)+A3*Dx(5)+A1*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
        dw*(A4*Dx(1)+A3*Dx(2)+A2*Dx(3)+A2*Dx(4)+A3*Dx(5)+A1*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
        dw*(A4*Dx(1)+A3*Dx(2)+A3*Dx(3)+A2*Dx(4)+A2*Dx(5)+A1*Dx(6))


      ! Now we need to calculate the DOFs of the fine grid which lie
      ! inside the coarse grid tetrahedron.
      idx(1) = IfacesAtElementFine(3, ielf(1)) - NVTfine - NMTfine
      idx(2) = IfacesAtElementFine(3, ielf(2)) - NVTfine - NMTfine
      idx(3) = IfacesAtElementFine(3, ielf(3)) - NVTfine - NMTfine
      idx(4) = IfacesAtElementFine(3, ielf(4)) - NVTfine - NMTfine
      DuFine(idx(1))=A5*Dx(1)+A5*Dx(2)+A6*Dx(3)+A7*Dx(4)+A6*Dx(5)+A7*Dx(6)
      DuFine(idx(2))=A5*Dx(1)+A6*Dx(2)+A5*Dx(3)+A6*Dx(4)+A7*Dx(5)+A7*Dx(6)
      DuFine(idx(3))=A5*Dx(1)+A7*Dx(2)+A6*Dx(3)+A5*Dx(4)+A6*Dx(5)+A7*Dx(6)
      DuFine(idx(4))=A5*Dx(1)+A6*Dx(2)+A7*Dx(3)+A6*Dx(4)+A5*Dx(5)+A7*Dx(6)

      idx(1) = IfacesAtElementFine(6, ielf(1)) - NVTfine - NMTfine
      idx(2) = IfacesAtElementFine(6, ielf(2)) - NVTfine - NMTfine
      idx(3) = IfacesAtElementFine(6, ielf(3)) - NVTfine - NMTfine
      idx(4) = IfacesAtElementFine(6, ielf(4)) - NVTfine - NMTfine
      DuFine(idx(1))=A6*Dx(1)+A5*Dx(2)+A7*Dx(3)+A7*Dx(4)+A5*Dx(5)+A6*Dx(6)
      DuFine(idx(2))=A6*Dx(1)+A5*Dx(2)+A5*Dx(3)+A7*Dx(4)+A7*Dx(5)+A6*Dx(6)
      DuFine(idx(3))=A6*Dx(1)+A7*Dx(2)+A5*Dx(3)+A5*Dx(4)+A7*Dx(5)+A6*Dx(6)
      DuFine(idx(4))=A6*Dx(1)+A7*Dx(2)+A7*Dx(3)+A5*Dx(4)+A5*Dx(5)+A6*Dx(6)

      idx(1) = IfacesAtElementFine(3, ielf(5)) - NVTfine - NMTfine
      idx(2) = IfacesAtElementFine(3, ielf(6)) - NVTfine - NMTfine
      idx(3) = IfacesAtElementFine(3, ielf(7)) - NVTfine - NMTfine
      idx(4) = IfacesAtElementFine(3, ielf(8)) - NVTfine - NMTfine
      DuFine(idx(1))=A7*Dx(1)+A5*Dx(2)+A6*Dx(3)+A7*Dx(4)+A6*Dx(5)+A5*Dx(6)
      DuFine(idx(2))=A7*Dx(1)+A6*Dx(2)+A5*Dx(3)+A6*Dx(4)+A7*Dx(5)+A5*Dx(6)
      DuFine(idx(3))=A7*Dx(1)+A7*Dx(2)+A6*Dx(3)+A5*Dx(4)+A6*Dx(5)+A5*Dx(6)
      DuFine(idx(4))=A7*Dx(1)+A6*Dx(2)+A7*Dx(3)+A6*Dx(4)+A5*Dx(5)+A5*Dx(6)

    END DO

  END SUBROUTINE

!<subroutine>

  SUBROUTINE mlprj_restUniformEx3x_3D_double (DuCoarse,DuFine, &
               IfacesAtElementCoarse,IfacesAtElementFine,&
               IneighboursAtElementCoarse,NVTcoarse,NVTfine,&
               NMTcoarse,NMTfine,NELcoarse,ielType)
  
!<description>
  ! Restrict a defect vector from a fine grid to a coarse grid.
  ! E031/EM31, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IfacesAtElement array on the coarse grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IfacesAtElementCoarse
  
  ! IfacesAtElement array on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IfacesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse

  ! Number of vertices in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine

  ! Number of edges in the coarse grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTcoarse
  
  ! Number of edges in the fine grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse

  ! Element type to use for prolongation
  INTEGER, INTENT(IN) :: ielType
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP) :: dw, A1, A2, A3, A4, A5, A6, A7
  REAL(DP), DIMENSION(36) :: Dx
  INTEGER(PREC_EDGEIDX), DIMENSION(6) :: IM
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: ielf
  
    ! Weights for the restriction; all coefficients are halfed, so dividing
    ! by 2 is not necessary in the calculation routines.
    IF (IAND(ielType, 2**16) .EQ. 0) THEN
      ! Weights for Ex31 element
      A1 =  0.458333333333333_DP   ! = 11/24
      A2 =  0.145833333333333_DP   ! =  7/48
      A3 = -0.104166666666667_DP   ! = -5/48
      A4 = -0.041666666666667_DP   ! = -1/24
      A5 =  0.458333333333333_DP   ! = 11/24
      A6 =  0.083333333333333_DP   ! =  1/12
      A7 = -0.041666666666667_DP   ! = -1/24
    ELSE
      ! Weights for Ex30 element
      A1 =  0.5_DP
      A2 =  0.125_DP
      A3 = -0.125_DP
      A4 =  0.0_DP
      A5 =  0.5_DP
      A6 =  0.0_DP
      A7 =  0.0_DP
    END IF
  
    ! Clear the output vector
    CALL lalg_clearVectorDble(DuCoarse)
  
    ! Loop over the coarse grid elements
    DO iel = 1, NELcoarse

      ! Get the DOF's of the coarse grid element
      IM(1) = IfacesAtElementCoarse(1,iel) - NVTcoarse - NMTcoarse
      IM(2) = IfacesAtElementCoarse(2,iel) - NVTcoarse - NMTcoarse
      IM(3) = IfacesAtElementCoarse(3,iel) - NVTcoarse - NMTcoarse
      IM(4) = IfacesAtElementCoarse(4,iel) - NVTcoarse - NMTcoarse
      IM(5) = IfacesAtElementCoarse(5,iel) - NVTcoarse - NMTcoarse
      IM(6) = IfacesAtElementCoarse(6,iel) - NVTcoarse - NMTcoarse

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = NELcoarse + 7*(iel-1) + 1
      ielf(3) = ielf(2) + 1
      ielf(4) = ielf(3) + 1
      ielf(5) = ielf(4) + 1
      ielf(6) = ielf(5) + 1
      ielf(7) = ielf(6) + 1
      ielf(8) = ielf(7) + 1

      ! Get the fine grid DOFs
      Dx( 1)=DuFine(IfacesAtElementFine(1,ielf(1))-NVTfine - NMTfine)
      Dx( 2)=DuFine(IfacesAtElementFine(1,ielf(2))-NVTfine - NMTfine)
      Dx( 3)=DuFine(IfacesAtElementFine(1,ielf(3))-NVTfine - NMTfine)
      Dx( 4)=DuFine(IfacesAtElementFine(1,ielf(4))-NVTfine - NMTfine)
      Dx( 5)=DuFine(IfacesAtElementFine(2,ielf(1))-NVTfine - NMTfine)
      Dx( 6)=DuFine(IfacesAtElementFine(5,ielf(2))-NVTfine - NMTfine)
      Dx( 7)=DuFine(IfacesAtElementFine(5,ielf(6))-NVTfine - NMTfine)
      Dx( 8)=DuFine(IfacesAtElementFine(2,ielf(5))-NVTfine - NMTfine)
      Dx( 9)=DuFine(IfacesAtElementFine(2,ielf(2))-NVTfine - NMTfine)
      Dx(10)=DuFine(IfacesAtElementFine(5,ielf(3))-NVTfine - NMTfine)
      Dx(11)=DuFine(IfacesAtElementFine(5,ielf(7))-NVTfine - NMTfine)
      Dx(12)=DuFine(IfacesAtElementFine(2,ielf(6))-NVTfine - NMTfine)
      Dx(13)=DuFine(IfacesAtElementFine(2,ielf(3))-NVTfine - NMTfine)
      Dx(14)=DuFine(IfacesAtElementFine(5,ielf(4))-NVTfine - NMTfine)
      Dx(15)=DuFine(IfacesAtElementFine(5,ielf(8))-NVTfine - NMTfine)
      Dx(16)=DuFine(IfacesAtElementFine(2,ielf(7))-NVTfine - NMTfine)
      Dx(17)=DuFine(IfacesAtElementFine(2,ielf(4))-NVTfine - NMTfine)
      Dx(18)=DuFine(IfacesAtElementFine(5,ielf(1))-NVTfine - NMTfine)
      Dx(19)=DuFine(IfacesAtElementFine(5,ielf(5))-NVTfine - NMTfine)
      Dx(20)=DuFine(IfacesAtElementFine(2,ielf(8))-NVTfine - NMTfine)
      Dx(21)=DuFine(IfacesAtElementFine(1,ielf(5))-NVTfine - NMTfine)
      Dx(22)=DuFine(IfacesAtElementFine(1,ielf(6))-NVTfine - NMTfine)
      Dx(23)=DuFine(IfacesAtElementFine(1,ielf(7))-NVTfine - NMTfine)
      Dx(24)=DuFine(IfacesAtElementFine(1,ielf(8))-NVTfine - NMTfine)
      Dx(25)=DuFine(IfacesAtElementFine(3,ielf(1))-NVTfine - NMTfine)
      Dx(26)=DuFine(IfacesAtElementFine(3,ielf(2))-NVTfine - NMTfine)
      Dx(27)=DuFine(IfacesAtElementFine(3,ielf(3))-NVTfine - NMTfine)
      Dx(28)=DuFine(IfacesAtElementFine(3,ielf(4))-NVTfine - NMTfine)
      Dx(29)=DuFine(IfacesAtElementFine(6,ielf(1))-NVTfine - NMTfine)
      Dx(30)=DuFine(IfacesAtElementFine(6,ielf(2))-NVTfine - NMTfine)
      Dx(31)=DuFine(IfacesAtElementFine(6,ielf(3))-NVTfine - NMTfine)
      Dx(32)=DuFine(IfacesAtElementFine(6,ielf(4))-NVTfine - NMTfine)
      Dx(33)=DuFine(IfacesAtElementFine(3,ielf(5))-NVTfine - NMTfine)
      Dx(34)=DuFine(IfacesAtElementFine(3,ielf(6))-NVTfine - NMTfine)
      Dx(35)=DuFine(IfacesAtElementFine(3,ielf(7))-NVTfine - NMTfine)
      Dx(36)=DuFine(IfacesAtElementFine(3,ielf(8))-NVTfine - NMTfine)
      
      ! Face 1
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(1,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(1))= DuCoarse(IM(1)) + dw * (&
        +A1*(Dx(1)+Dx(2)+Dx(3)+Dx(4))&
        +A2*(Dx(5)+Dx(6)+Dx(9)+Dx(10)+Dx(13)+Dx(14)+Dx(17)+Dx(18))&
        +A3*(Dx(7)+Dx(8)+Dx(11)+Dx(12)+Dx(15)+Dx(16)+Dx(19)+Dx(20))&
        +A4*(Dx(21)+Dx(22)+Dx(23)+Dx(24))&
        +A5*(Dx(25)+Dx(26)+Dx(27)+Dx(28))&
        +A6*(Dx(29)+Dx(30)+Dx(31)+Dx(32))&
        +A7*(Dx(33)+Dx(34)+Dx(35)+Dx(36)))    

      ! Face 2
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(2,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(2))= DuCoarse(IM(2)) + dw * (&
        +A1*(Dx(5)+Dx(6)+Dx(7)+Dx(8))&
        +A2*(Dx(1)+Dx(2)+Dx(9)+Dx(12)+Dx(21)+Dx(22)+Dx(18)+Dx(19))&
        +A3*(Dx(3)+Dx(4)+Dx(10)+Dx(11)+Dx(23)+Dx(24)+Dx(17)+Dx(20))&
        +A4*(Dx(13)+Dx(14)+Dx(15)+Dx(16))&
        +A5*(Dx(25)+Dx(29)+Dx(30)+Dx(33))&
        +A6*(Dx(26)+Dx(28)+Dx(34)+Dx(36))&
        +A7*(Dx(27)+Dx(31)+Dx(32)+Dx(35)))    

      ! Face 3
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(3,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(3))= DuCoarse(IM(3)) + dw * (&
        +A1*(Dx(9)+Dx(10)+Dx(11)+Dx(12))&
        +A2*(Dx(2)+Dx(3)+Dx(6)+Dx(7)+Dx(22)+Dx(23)+Dx(13)+Dx(16))&
        +A3*(Dx(1)+Dx(4)+Dx(5)+Dx(8)+Dx(21)+Dx(24)+Dx(14)+Dx(15))&
        +A4*(Dx(17)+Dx(18)+Dx(19)+Dx(20))&
        +A5*(Dx(26)+Dx(30)+Dx(31)+Dx(34))&
        +A6*(Dx(25)+Dx(27)+Dx(33)+Dx(35))&
        +A7*(Dx(28)+Dx(29)+Dx(32)+Dx(36)))    

      ! Face 4
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(4,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(4))= DuCoarse(IM(4)) + dw * (&
        +A1*(Dx(13)+Dx(14)+Dx(15)+Dx(16))&
        +A2*(Dx(3)+Dx(4)+Dx(10)+Dx(11)+Dx(23)+Dx(24)+Dx(17)+Dx(20))&
        +A3*(Dx(1)+Dx(2)+Dx(9)+Dx(12)+Dx(21)+Dx(22)+Dx(18)+Dx(19))&
        +A4*(Dx(5)+Dx(6)+Dx(7)+Dx(8))&
        +A5*(Dx(27)+Dx(31)+Dx(32)+Dx(35))&
        +A6*(Dx(26)+Dx(28)+Dx(34)+Dx(36))&
        +A7*(Dx(25)+Dx(29)+Dx(30)+Dx(33)))

      ! Face 5
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(5,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(5))= DuCoarse(IM(5)) + dw * (&
        +A1*(Dx(17)+Dx(18)+Dx(19)+Dx(20))&
        +A2*(Dx(1)+Dx(4)+Dx(5)+Dx(8)+Dx(21)+Dx(24)+Dx(14)+Dx(15))&
        +A3*(Dx(2)+Dx(3)+Dx(6)+Dx(7)+Dx(22)+Dx(23)+Dx(13)+Dx(16))&
        +A4*(Dx(9)+Dx(10)+Dx(11)+Dx(12))&
        +A5*(Dx(28)+Dx(29)+Dx(32)+Dx(36))&     
        +A6*(Dx(25)+Dx(27)+Dx(33)+Dx(35))&
        +A7*(Dx(26)+Dx(30)+Dx(31)+Dx(34)))

      ! Face 6
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(6,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(6))= DuCoarse(IM(6)) + dw * (&
        +A1*(Dx(21)+Dx(22)+Dx(23)+Dx(24))&
        +A2*(Dx(7)+Dx(8)+Dx(11)+Dx(12)+Dx(15)+Dx(16)+Dx(19)+Dx(20))&
        +A3*(Dx(5)+Dx(6)+Dx(9)+Dx(10)+Dx(13)+Dx(14)+Dx(17)+Dx(18))&
        +A4*(Dx(1)+Dx(2)+Dx(3)+Dx(4))&
        +A5*(Dx(33)+Dx(34)+Dx(35)+Dx(36))&
        +A6*(Dx(29)+Dx(30)+Dx(31)+Dx(32))&
        +A7*(Dx(25)+Dx(26)+Dx(27)+Dx(28)))

    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mlprj_interpUniformEx3x_3D_dbl (DuCoarse,DuFine, &
               IfacesAtElementCoarse,IfacesAtElementFine,&
               IneighboursAtElementCoarse,NVTcoarse,NVTfine,&
               NMTcoarse,NMTfine,NELcoarse,ielType)
  
!<description>
  ! Restrict a solution vector from a fine grid to a coarse grid.
  ! E031/EM31, uniform triangulation, double precision vector.
!</description>
  
!<input>
  ! Fine grid vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
  
  ! IfacesAtElement array on the coarse grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IfacesAtElementCoarse
  
  ! IfacesAtElement array on the fine grid
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IfacesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
  
  ! Number of vertices in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse

  ! Number of vertices in the fine grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine

  ! Number of edges in the coarse grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTcoarse
  
  ! Number of edges in the fine grid
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMTfine

  ! Number of elements in the coarse grid
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse

  ! Element type to use for prolongation
  INTEGER, INTENT(IN) :: ielType
!</input>
  
!<output>
  ! Coarse grid vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!</output>
  
!</subroutine>
  
  ! local variables
  REAL(DP) :: dw, A1, A2, A3
  REAL(DP), DIMENSION(36) :: Dx
  INTEGER(PREC_EDGEIDX), DIMENSION(6) :: IM
  INTEGER(PREC_ELEMENTIDX) :: iel
  INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: ielf
  
    ! Weights for the restriction; all coefficients are halfed, so dividing
    ! by 2 is not necessary in the calculation routines.
!    IF (IAND(ielType, 2**16) .EQ. 0) THEN
!      ! Weights for Ex31 element
!      A1 =  0.458333333333333_DP   ! = 11/24
!      A2 =  0.145833333333333_DP   ! =  7/48
!      A3 = -0.104166666666667_DP   ! = -5/48
!      A4 = -0.041666666666667_DP   ! = -1/24
!      A5 =  0.458333333333333_DP   ! = 11/24
!      A6 =  0.083333333333333_DP   ! =  1/12
!      A7 = -0.041666666666667_DP   ! = -1/24
!    ELSE
      ! Weights for Ex30 element
      A1 =  0.125_DP
      A2 =  0.25_DP
      A3 =  0.0833333333333333_DP
!    END IF
  
    ! Clear the output vector
    CALL lalg_clearVectorDble(DuCoarse)
  
    ! Loop over the coarse grid elements
    DO iel = 1, NELcoarse

      ! Get the DOF's of the coarse grid element
      IM(1) = IfacesAtElementCoarse(1,iel) - NVTcoarse - NMTcoarse
      IM(2) = IfacesAtElementCoarse(2,iel) - NVTcoarse - NMTcoarse
      IM(3) = IfacesAtElementCoarse(3,iel) - NVTcoarse - NMTcoarse
      IM(4) = IfacesAtElementCoarse(4,iel) - NVTcoarse - NMTcoarse
      IM(5) = IfacesAtElementCoarse(5,iel) - NVTcoarse - NMTcoarse
      IM(6) = IfacesAtElementCoarse(6,iel) - NVTcoarse - NMTcoarse

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = NELcoarse + 7*(iel-1) + 1
      ielf(3) = ielf(2) + 1
      ielf(4) = ielf(3) + 1
      ielf(5) = ielf(4) + 1
      ielf(6) = ielf(5) + 1
      ielf(7) = ielf(6) + 1
      ielf(8) = ielf(7) + 1

      ! Get the fine grid DOFs
      Dx(1)=DuFine(IfacesAtElementFine(1,ielf(1))-NVTfine - NMTfine)
      Dx(2)=DuFine(IfacesAtElementFine(2,ielf(1))-NVTfine - NMTfine)
      Dx(3)=DuFine(IfacesAtElementFine(3,ielf(1))-NVTfine - NMTfine)
      Dx(4)=DuFine(IfacesAtElementFine(4,ielf(1))-NVTfine - NMTfine)
      Dx(5)=DuFine(IfacesAtElementFine(5,ielf(1))-NVTfine - NMTfine)
      Dx(6)=DuFine(IfacesAtElementFine(6,ielf(1))-NVTfine - NMTfine)
      Dx(7)=DuFine(IfacesAtElementFine(1,ielf(2))-NVTfine - NMTfine)
      Dx(8)=DuFine(IfacesAtElementFine(2,ielf(2))-NVTfine - NMTfine)
      Dx(9)=DuFine(IfacesAtElementFine(3,ielf(2))-NVTfine - NMTfine)
      Dx(10)=DuFine(IfacesAtElementFine(5,ielf(2))-NVTfine - NMTfine)
      Dx(11)=DuFine(IfacesAtElementFine(6,ielf(2))-NVTfine - NMTfine)
      Dx(12)=DuFine(IfacesAtElementFine(1,ielf(3))-NVTfine - NMTfine)
      Dx(13)=DuFine(IfacesAtElementFine(2,ielf(3))-NVTfine - NMTfine)
      Dx(14)=DuFine(IfacesAtElementFine(3,ielf(3))-NVTfine - NMTfine)
      Dx(15)=DuFine(IfacesAtElementFine(5,ielf(3))-NVTfine - NMTfine)
      Dx(16)=DuFine(IfacesAtElementFine(6,ielf(3))-NVTfine - NMTfine)
      Dx(17)=DuFine(IfacesAtElementFine(1,ielf(4))-NVTfine - NMTfine)
      Dx(18)=DuFine(IfacesAtElementFine(2,ielf(4))-NVTfine - NMTfine)
      Dx(19)=DuFine(IfacesAtElementFine(5,ielf(4))-NVTfine - NMTfine)
      Dx(20)=DuFine(IfacesAtElementFine(6,ielf(4))-NVTfine - NMTfine)
      Dx(21)=DuFine(IfacesAtElementFine(1,ielf(5))-NVTfine - NMTfine)
      Dx(22)=DuFine(IfacesAtElementFine(2,ielf(5))-NVTfine - NMTfine)
      Dx(23)=DuFine(IfacesAtElementFine(3,ielf(5))-NVTfine - NMTfine)
      Dx(24)=DuFine(IfacesAtElementFine(4,ielf(5))-NVTfine - NMTfine)
      Dx(25)=DuFine(IfacesAtElementFine(5,ielf(5))-NVTfine - NMTfine)
      Dx(26)=DuFine(IfacesAtElementFine(1,ielf(6))-NVTfine - NMTfine)
      Dx(27)=DuFine(IfacesAtElementFine(2,ielf(6))-NVTfine - NMTfine)
      Dx(28)=DuFine(IfacesAtElementFine(3,ielf(6))-NVTfine - NMTfine)
      Dx(29)=DuFine(IfacesAtElementFine(5,ielf(6))-NVTfine - NMTfine)
      Dx(30)=DuFine(IfacesAtElementFine(1,ielf(7))-NVTfine - NMTfine)
      Dx(31)=DuFine(IfacesAtElementFine(2,ielf(7))-NVTfine - NMTfine)
      Dx(32)=DuFine(IfacesAtElementFine(3,ielf(7))-NVTfine - NMTfine)
      Dx(33)=DuFine(IfacesAtElementFine(5,ielf(7))-NVTfine - NMTfine)
      Dx(34)=DuFine(IfacesAtElementFine(1,ielf(8))-NVTfine - NMTfine)
      Dx(35)=DuFine(IfacesAtElementFine(2,ielf(8))-NVTfine - NMTfine)
      Dx(36)=DuFine(IfacesAtElementFine(5,ielf(8))-NVTfine - NMTfine)

      ! Face 1
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(1,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(1))= DuCoarse(IM(1)) + dw * (&
         A1*(Dx(1)+Dx(7)+Dx(12)+Dx(17))+&
         A2*(Dx(3)+Dx(14)+Dx(4)+Dx(9))-&
         A3*(Dx(2)+Dx(10)+Dx(8)+Dx(15)+Dx(13)+Dx(19)+&
             Dx(18)+Dx(5)+Dx(6)+Dx(11)+Dx(16)+Dx(20)))

      ! Face 2
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(2,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(2))= DuCoarse(IM(2)) + dw * (&
          A1*(Dx(2)+Dx(10)+Dx(29)+Dx(22))+&
          A2*(Dx(3)+Dx(23)+Dx(6)+Dx(11))-&
          A3*(Dx(1)+Dx(7)+Dx(8)+Dx(27)+Dx(26)+Dx(21)+&
              Dx(25)+Dx(5)+Dx(4)+Dx(9)+Dx(24)+Dx(28)))

      ! Face 3
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(3,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(3))= DuCoarse(IM(3)) + dw * (&
          A1*(Dx(8)+Dx(15)+Dx(33)+Dx(27))+&
          A2*(Dx(11)+Dx(16)+Dx(9)+Dx(28))-&
          A3*(Dx(3)+Dx(14)+Dx(32)+Dx(23)+Dx(7)+Dx(12)+&
              Dx(13)+Dx(31)+Dx(30)+Dx(26)+Dx(29)+Dx(10)))

      ! Face 4
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(4,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(4))= DuCoarse(IM(4)) + dw * (&
        A1*(Dx(13)+Dx(31)+Dx(36)+Dx(19))+&
        A2*(Dx(32)+Dx(14)+Dx(16)+Dx(20))-&
        A3*(Dx(17)+Dx(12)+Dx(15)+Dx(33)+Dx(30)+Dx(34)+&
            Dx(35)+Dx(18)+Dx(4)+Dx(9)+Dx(28)+Dx(24)))

      ! Face 5
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(5,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(5))= DuCoarse(IM(5)) + dw * (&
         A1*(Dx(5)+Dx(18)+Dx(35)+Dx(25))+&
         A2*(Dx(6)+Dx(20)+Dx(4)+Dx(24))-&
         A3*(Dx(1)+Dx(17)+Dx(19)+Dx(36)+Dx(34)+Dx(21)+&
             Dx(22)+Dx(2)+Dx(3)+Dx(14)+Dx(32)+Dx(23)))

      ! Face 6
      dw = 1.0_DP
      IF (IneighboursAtElementCoarse(6,iel) .EQ. 0) dw = 2.0_DP
      DuCoarse(IM(6))= DuCoarse(IM(6)) + dw *  (&
         A1*(Dx(26)+Dx(30)+Dx(34)+Dx(21))+&
         A2*(Dx(23)+Dx(32)+Dx(24)+Dx(28))-&
         A3*(Dx(22)+Dx(29)+Dx(27)+Dx(33)+Dx(31)+Dx(36)+&
             Dx(21)+Dx(25)+Dx(6)+Dx(11)+Dx(16)+Dx(20)))

    END DO

  END SUBROUTINE
  
!  ! ***************************************************************************
!  ! Support for Q2~ element with bubble
!  ! ***************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE mlprj_prolUniformE037_double (DuCoarse,DuFine, &
!               IedgesAtElementCoarse,IedgesAtElementFine,&
!               IneighboursAtElementCoarse,IneighboursAtElementFine,&
!               ItwistCoarse,ItwistFine,NVTcoarse,NVTfine,&
!               NMTcoarse,NMTfine,NELcoarse,NELfine)
!  
!!<description>
!  ! Prolongate a solution vector from a coarse grid to a fine grid.
!  ! E037, uniform triangulation, double precision vector.
!!</description>
!  
!!<input>
!  ! Coarse grid vector
!  REAL(DP), DIMENSION(:), INTENT(IN) :: DuCoarse
!  
!  ! IedgesAtElement array on the coarse grid
!  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementCoarse
!  
!  ! IedgesAtElement array on the fine grid
!  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine
!
!  ! IneighboursAtElement array on the coarse grid
!  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
!  
!  ! IneighboursAtElement array on the fine grid
!  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
!  
!  ! ItwistIndexEdges array on the coarse grid
!  INTEGER(I32), DIMENSION(:), INTENT(IN) :: ItwistCoarse
!
!  ! ItwistIndexEdges array on the fine grid
!  INTEGER(I32), DIMENSION(:), INTENT(IN) :: ItwistFine
!
!  ! Number of vertices in the coarse grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse
!
!  ! Number of vertices in the fine grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine
!  
!  ! Number of egdes in the coarse grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NMTcoarse
!  
!  ! Number of egdes in the fine grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NMTfine
!
!  ! Number of elements in the coarse grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!
!  ! Number of elements in the fine grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELfine
!!</input>
!  
!!<output>
!  ! Fine grid vector
!  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuFine
!!</output>
!  
!!</subroutine>
!  
!  ! local variables
!  INTEGER :: iel,i
!  INTEGER, DIMENSION(4) :: Ielf, Imc
!  
!  ! Local vectors
!  REAL(DP), DIMENSION(10) :: Dv
!  REAL(DP), DIMENSION(4:6,4) :: Dtf
!  REAL(DP), DIMENSION(4) :: Dtn
!  INTEGER, DIMENSION(8,4) :: Idf
!  
!    ! Clear the output vector
!    CALL lalg_clearVectorDble(DuFine)
!    
!    ! Loop over the coarse grid elements
!    DO iel = 1, NELcoarse
!    
!      ! Get the element numbers of the fine grid elements
!      Ielf(1) = iel
!      Ielf(2) = IneighboursAtElementFine(2,Ielf(1))
!      Ielf(3) = IneighboursAtElementFine(2,Ielf(2))
!      Ielf(4) = IneighboursAtElementFine(2,Ielf(3))
!
!      ! Get the edge indices of the coarse grid element
!      Imc(1) = IedgesAtElementCoarse(1,iel)-NVTcoarse
!      Imc(2) = IedgesAtElementCoarse(2,iel)-NVTcoarse
!      Imc(3) = IedgesAtElementCoarse(3,iel)-NVTcoarse
!      Imc(4) = IedgesAtElementCoarse(4,iel)-NVTcoarse
!      
!      Dtn = 1.0_DP
!      DO i = 1, 4
!      
!        ! Calculate the Edge-Int-Mean DOFs on the fine mesh
!        Idf(1,i) = IedgesAtElementFine(1,Ielf(i))-NVTfine
!        Idf(2,i) = IedgesAtElementFine(2,Ielf(i))-NVTfine
!        Idf(3,i) = IedgesAtElementFine(4,Ielf(i))-NVTfine
!        ! Calculate the Edge-Legendre-Int-Mean DOFs on the fine mesh
!        Idf(4,i) = Idf(1,i) + NMTfine
!        Idf(5,i) = Idf(2,i) + NMTfine
!        Idf(6,i) = Idf(3,i) + NMTfine
!        ! Calculate the Quad-Int-Mean DOF on the fine mesh
!        Idf(7,i) = 2*NMTfine + Ielf(i)
!        ! Calculate the Quad-Legendre-Int-Mean DOF on the fine mesh
!        Idf(8,i) = 2*NMTfine + NELfine + Ielf(i)
!        
!        ! Calculate the twist factors for the Edge-Legendre-Int-Mean DOFs
!        ! on the fine mesh
!        Dtf(4,i) = REAL(2*IAND(      ItwistFine(Ielf(i))   ,1)-1,DP)
!        Dtf(5,i) = REAL(2*IAND(ISHFT(ItwistFine(Ielf(i)),1),1)-1,DP)
!        Dtf(6,i) = REAL(2*IAND(ISHFT(ItwistFine(Ielf(i)),3),1)-1,DP)
!        
!        ! If there is a neighbour at the coarse mesh edge, then we
!        ! set the corresponding factor to 1/2
!        IF (IneighboursAtElementCoarse(i,iel) .NE. 0) Dtn(i) = 0.5_DP
!
!      END DO
!
!      ! Get the values of the corresponding coarse mesh DOFs
!      Dv(1:4) = DuCoarse(Imc(1:4))
!      Dv(  5) = DuCoarse(Imc(1) + NMTcoarse) &
!              * REAL(2*IAND(      ItwistCoarse(iel)   ,1)-1,DP)
!      Dv(  6) = DuCoarse(Imc(2) + NMTcoarse) &
!              * REAL(2*IAND(ISHFT(ItwistCoarse(iel),1),1)-1,DP)
!      Dv(  7) = DuCoarse(Imc(3) + NMTcoarse) &
!              * REAL(2*IAND(ISHFT(ItwistCoarse(iel),2),1)-1,DP)
!      Dv(  8) = DuCoarse(Imc(4) + NMTcoarse) &
!              * REAL(2*IAND(ISHFT(ItwistCoarse(iel),3),1)-1,DP)
!      Dv(  9) = DuCoarse(2*NMTcoarse + iel)
!      Dv( 10) = DuCoarse(2*NMTcoarse + NELcoarse + iel)
!
!      ! Dofs 1 - 12 -> Edge-Integral-Means
!      DuFine(Idf(1,1)) = DuFine(Idf(1,1)) + Dtn(1)*(Dv(1) - 1.59375_DP*Dv(5) - &
!                              0.09375_DP*(Dv(6) + Dv(7) + Dv(8)))
!      DuFine(Idf(3,2)) = DuFine(Idf(3,2)) + Dtn(1)*(Dv(1) + 1.59375_DP*Dv(5) + &
!                              0.09375_DP*(Dv(6) + Dv(7) + Dv(8)))
!      DuFine(Idf(1,2)) = DuFine(Idf(1,2)) + Dtn(2)*(Dv(2) - 1.59375_DP*Dv(6) - &
!                              0.09375_DP*(Dv(5) + Dv(7) + Dv(8)))
!      DuFine(Idf(3,3)) = DuFine(Idf(3,3)) + Dtn(2)*(Dv(2) + 1.59375_DP*Dv(6) + &
!                              0.09375_DP*(Dv(5) + Dv(7) + Dv(8)))
!      DuFine(Idf(1,3)) = DuFine(Idf(1,3)) + Dtn(3)*(Dv(3) - 1.59375_DP*Dv(7) - &
!                              0.09375_DP*(Dv(5) + Dv(6) + Dv(8)))
!      DuFine(Idf(3,4)) = DuFine(Idf(3,4)) + Dtn(3)*(Dv(3) + 1.59375_DP*Dv(7) + &
!                              0.09375_DP*(Dv(5) + Dv(6) + Dv(8)))
!      DuFine(Idf(1,4)) = DuFine(Idf(1,4)) + Dtn(4)*(Dv(4) - 1.59375_DP*Dv(8) - &
!                              0.09375_DP*(Dv(5) + Dv(6) + Dv(7)))
!      DuFine(Idf(3,1)) = DuFine(Idf(3,1)) + Dtn(4)*(Dv(4) + 1.59375_DP*Dv(8) + &
!                              0.09375_DP*(Dv(5) + Dv(6) + Dv(7)))
!      DuFine(Idf(2,1)) = DuFine(Idf(2,1)) - 0.25_DP*(Dv(2) + Dv(4)) + &
!                         1.5_DP*Dv(9) + 0.375_DP*(Dv(1) - Dv(3) + Dv(6) - Dv(8))
!      DuFine(Idf(2,2)) = DuFine(Idf(2,2)) - 0.25_DP*(Dv(1) + Dv(3)) + &
!                         1.5_DP*Dv(9) + 0.375_DP*(Dv(2) - Dv(4) - Dv(5) + Dv(7))
!      DuFine(Idf(2,3)) = DuFine(Idf(2,3)) - 0.25_DP*(Dv(2) + Dv(4)) + &
!                         1.5_DP*Dv(9) - 0.375_DP*(Dv(1) - Dv(3) + Dv(6) - Dv(8))
!      DuFine(Idf(2,4)) = DuFine(Idf(2,4)) - 0.25_DP*(Dv(1) + Dv(3)) + &
!                         1.5_DP*Dv(9) - 0.375_DP*(Dv(2) - Dv(4) - Dv(5) + Dv(7))
!
!      ! Dofs 13 - 24 -> Legendre-Int-Means
!      DuFine(Idf(4,1)) = DuFine(Idf(4,1)) + Dtn(1)*Dtf(4,1)*(&
!          0.125_DP*(Dv(1)-Dv(2)-Dv(3)-Dv(4)) + 0.25_DP*Dv(9) - 6.25_DP*Dv(10)&
!        + 0.40625_DP*Dv(5) + 0.28125_DP*Dv(6) - 0.09375_DP*Dv(7) - 0.46875_DP*Dv(8))
!      DuFine(Idf(6,2)) = DuFine(Idf(6,2)) - Dtn(1)*Dtf(6,2)*(&
!          0.125_DP*(-Dv(1)+Dv(2)+Dv(3)+Dv(4)) - 0.25_DP*Dv(9) + 6.25_DP*Dv(10)&
!        + 0.40625_DP*Dv(5) - 0.46875_DP*Dv(6) - 0.09375_DP*Dv(7) + 0.28125_DP*Dv(8))
!      DuFine(Idf(4,2)) = DuFine(Idf(4,2)) + Dtn(2)*Dtf(4,2)*(&
!          0.125_DP*(-Dv(1)+Dv(2)-Dv(3)-Dv(4)) + 0.25_DP*Dv(9) - 6.25_DP*Dv(10)&
!        - 0.46875_DP*Dv(5) + 0.40625_DP*Dv(6) + 0.28125_DP*Dv(7) - 0.09375_DP*Dv(8))
!      DuFine(Idf(6,3)) = DuFine(Idf(6,3)) - Dtn(2)*Dtf(6,3)*(&
!          0.125_DP*(Dv(1)-Dv(2)+Dv(3)+Dv(4)) - 0.25_DP*Dv(9) + 6.25_DP*Dv(10)&
!        + 0.28125_DP*Dv(5) + 0.40625_DP*Dv(6) - 0.46875_DP*Dv(7) - 0.09375_DP*Dv(8))
!      DuFine(Idf(4,3)) = DuFine(Idf(4,3)) + Dtn(3)*Dtf(4,3)*(&
!          0.125_DP*(-Dv(1)-Dv(2)+Dv(3)-Dv(4)) + 0.25_DP*Dv(9) - 6.25_DP*Dv(10)&
!        - 0.09375_DP*Dv(5) - 0.46875_DP*Dv(6) + 0.40625_DP*Dv(7) + 0.28125_DP*Dv(8))
!      DuFine(Idf(6,4)) = DuFine(Idf(6,4)) + Dtn(3)*Dtf(6,4)*(&
!          0.125_DP*(-Dv(1)-Dv(2)+Dv(3)-Dv(4)) + 0.25_DP*Dv(9) - 6.25_DP*Dv(10)&
!        + 0.09375_DP*Dv(5) - 0.28125_DP*Dv(6) - 0.40625_DP*Dv(7) + 0.46875_DP*Dv(8))
!      DuFine(Idf(4,4)) = DuFine(Idf(4,4)) + Dtn(4)*Dtf(4,4)*(&
!          0.125_DP*(-Dv(1)-Dv(2)+Dv(3)-Dv(4)) + 0.25_DP*Dv(9) - 6.25_DP*Dv(10)&
!        - 0.09375_DP*Dv(5) - 0.46875_DP*Dv(6) + 0.40625_DP*Dv(7) + 0.28125_DP*Dv(8))
!      DuFine(Idf(6,1)) = DuFine(Idf(6,1)) + Dtn(4)*Dtf(6,1)*(&
!          0.125_DP*(-Dv(1)+Dv(2)-Dv(3)-Dv(4)) + 0.25_DP*Dv(9) - 6.25_DP*Dv(10)&
!        - 0.28125_DP*Dv(5) - 0.40625_DP*Dv(6) + 0.46875_DP*Dv(7) + 0.09375_DP*Dv(8))
!      DuFine(Idf(5,1)) = DuFine(Idf(5,1)) - Dtf(5,1)*(0.25_DP*(-Dv(1) + Dv(9)) &
!                            + 0.125_DP*(-Dv(6) + Dv(8)) + 3.125_DP*Dv(10))
!      DuFine(Idf(5,2)) = DuFine(Idf(5,2)) + Dtf(5,2)*(0.25_DP*(Dv(2) - Dv(9)) &
!                            + 0.125_DP*(-Dv(5) + Dv(7)) - 3.125_DP*Dv(10))
!      DuFine(Idf(5,3)) = DuFine(Idf(5,3)) - Dtf(5,3)*(0.25_DP*(-Dv(3) + Dv(9)) &
!                            + 0.125_DP*(Dv(6) - Dv(8)) + 3.125_DP*Dv(10))
!      DuFine(Idf(5,4)) = DuFine(Idf(5,4)) + Dtf(5,4)*(0.25_DP*(-Dv(4) + Dv(9)) &
!                            + 0.125_DP*(-Dv(5) + Dv(7)) + 3.125_DP*Dv(10))
!      
!      ! Dofs 25 - 28 -> Int-Means on Quads
!      DuFine(Idf(7,1)) = DuFine(Idf(7,1)) + 0.25_DP*(Dv(1)-Dv(2)-Dv(3)+Dv(4)) &
!                + 0.1875_DP*(-Dv(5) + Dv(6) - Dv(7) + Dv(8)) + Dv(9)
!      DuFine(Idf(7,2)) = DuFine(Idf(7,2)) + 0.25_DP*(Dv(1)+Dv(2)-Dv(3)-Dv(4)) &
!                + 0.1875_DP*(Dv(5) - Dv(6) + Dv(7) - Dv(8)) + Dv(9)
!      DuFine(Idf(7,3)) = DuFine(Idf(7,3)) + 0.25_DP*(-Dv(1)+Dv(2)+Dv(3)-Dv(4)) &
!                + 0.1875_DP*(-Dv(5) + Dv(6) - Dv(7) + Dv(8)) + Dv(9)
!      DuFine(Idf(7,4)) = DuFine(Idf(7,4)) + 0.25_DP*(-Dv(1)-Dv(2)+Dv(3)+Dv(4)) &
!                + 0.1875_DP*(Dv(5) - Dv(6) + Dv(7) - Dv(8)) + Dv(9)
!
!      ! Dofs 29 - 32 -> Legendre-Int-Means on Quads
!      DuFine(Idf(8,1)) = DuFine(Idf(8,1)) + 0.0625_DP*Dv(10)
!      DuFine(Idf(8,2)) = DuFine(Idf(8,2)) + 0.0625_DP*Dv(10)
!      DuFine(Idf(8,3)) = DuFine(Idf(8,3)) + 0.0625_DP*Dv(10)
!      DuFine(Idf(8,4)) = DuFine(Idf(8,4)) + 0.0625_DP*Dv(10)
!
!      ! Go for the next element
!
!    END DO
!    
!    ! That's it
!
!  END SUBROUTINE
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE mlprj_restUniformE037_double (DuCoarse,DuFine, &
!               IedgesAtElementCoarse,IedgesAtElementFine,&
!               IneighboursAtElementCoarse,IneighboursAtElementFine,&
!               ItwistCoarse,ItwistFine,NVTcoarse,NVTfine,&
!               NMTcoarse,NMTfine,NELcoarse,NELfine)
!  
!!<description>
!  ! Restricts a defect vector from a fine grid to a coarse grid.
!  ! E037, uniform triangulation, double precision vector.
!!</description>
!  
!!<input>
!  ! Fine grid vector
!  REAL(DP), DIMENSION(:), INTENT(IN) :: DuFine
!  
!  ! IedgesAtElement array on the coarse grid
!  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementCoarse
!  
!  ! IedgesAtElement array on the fine grid
!  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElementFine
!
!  ! IneighboursAtElement array on the coarse grid
!  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementCoarse
!  
!  ! IneighboursAtElement array on the fine grid
!  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IneighboursAtElementFine
!  
!  ! ItwistIndexEdges array on the coarse grid
!  INTEGER(I32), DIMENSION(:), INTENT(IN) :: ItwistCoarse
!
!  ! ItwistIndexEdges array on the fine grid
!  INTEGER(I32), DIMENSION(:), INTENT(IN) :: ItwistFine
!
!  ! Number of vertices in the coarse grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTcoarse
!
!  ! Number of vertices in the fine grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NVTfine
!  
!  ! Number of egdes in the coarse grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NMTcoarse
!  
!  ! Number of egdes in the fine grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NMTfine
!
!  ! Number of elements in the coarse grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELcoarse
!
!  ! Number of elements in the fine grid
!  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NELfine
!!</input>
!  
!!<output>
!  ! Coarse grid vector
!  REAL(DP), DIMENSION(:), INTENT(OUT) :: DuCoarse
!!</output>
!  
!!</subroutine>
!  
!  ! local variables
!  INTEGER :: iel,i
!  INTEGER, DIMENSION(4) :: Ielf
!  
!  ! local vectors
!  REAL(DP), DIMENSION(32) :: dv
!  REAL(DP), DIMENSION(4) :: Dtn,Dtf
!  INTEGER, DIMENSION(10) :: Idf
!
!    ! Clear the output vector
!    CALL lalg_clearVectorDble(DuCoarse)
!    
!    ! Loop over the coarse grid elements
!    DO iel = 1, NELcoarse
!    
!      ! Get the element numbers of the fine grid elements
!      Ielf(1) = iel
!      Ielf(2) = IneighboursAtElementFine(2,Ielf(1))
!      Ielf(3) = IneighboursAtElementFine(2,Ielf(2))
!      Ielf(4) = IneighboursAtElementFine(2,Ielf(3))
!      
!      ! Calculate the coarse grid DOFs
!      Idf(1:4) = IedgesAtElementCoarse(1:4,iel)-NVTcoarse
!      Idf(5:8) = Idf(1:4) + NMTcoarse
!      Idf(  9) = 2*NMTcoarse + iel
!      Idf( 10) = 2*NMTcoarse + NELcoarse + iel
!      
!      ! Calculate twist index factors and neighbour scales for the
!      ! coarse grid edges
!      Dtn = 1.0_DP
!      DO i = 1, 4
!        
!        ! Twist index factor
!        Dtf(i) = REAL(2*IAND(ISHFT(ItwistCoarse(iel),i),1)-1,DP)
!        
!        ! Neighbour factor
!        IF (IneighboursAtElementCoarse(i,iel) .NE. 0) Dtn(i) = 0.5_DP
!        
!      END DO
!      
!      ! Get the values of the corresponding fine mesh DOFs
!      Dv( 1) = DuFine(IedgesAtElementFine(1,Ielf(1))-NVTfine)
!      Dv( 2) = DuFine(IedgesAtElementFine(4,Ielf(2))-NVTfine)
!      Dv( 3) = DuFine(IedgesAtElementFine(1,Ielf(2))-NVTfine)
!      Dv( 4) = DuFine(IedgesAtElementFine(4,Ielf(3))-NVTfine)
!      Dv( 5) = DuFine(IedgesAtElementFine(1,Ielf(3))-NVTfine)
!      Dv( 6) = DuFine(IedgesAtElementFine(4,Ielf(4))-NVTfine)
!      Dv( 7) = DuFine(IedgesAtElementFine(1,Ielf(4))-NVTfine)
!      Dv( 8) = DuFine(IedgesAtElementFine(4,Ielf(1))-NVTfine)
!      Dv( 9) = DuFine(IedgesAtElementFine(2,Ielf(1))-NVTfine)
!      Dv(10) = DuFine(IedgesAtElementFine(2,Ielf(2))-NVTfine)
!      Dv(11) = DuFine(IedgesAtElementFine(2,Ielf(3))-NVTfine)
!      Dv(12) = DuFine(IedgesAtElementFine(2,Ielf(4))-NVTfine)
!      Dv(13) = DuFine(IedgesAtElementFine(1,Ielf(1))-NVTfine+NMTfine)&
!             * REAL(2*IAND(      ItwistFine(Ielf(1))   ,1)-1,DP)
!      Dv(14) = DuFine(IedgesAtElementFine(4,Ielf(2))-NVTfine+NMTfine)&
!             * REAL(2*IAND(ISHFT(ItwistFine(Ielf(2)),3),1)-1,DP)
!      Dv(15) = DuFine(IedgesAtElementFine(1,Ielf(2))-NVTfine+NMTfine)&
!             * REAL(2*IAND(      ItwistFine(Ielf(2))   ,1)-1,DP)
!      Dv(16) = DuFine(IedgesAtElementFine(4,Ielf(3))-NVTfine+NMTfine)&
!             * REAL(2*IAND(ISHFT(ItwistFine(Ielf(3)),3),1)-1,DP)
!      Dv(17) = DuFine(IedgesAtElementFine(1,Ielf(3))-NVTfine+NMTfine)&
!             * REAL(2*IAND(      ItwistFine(Ielf(3))   ,1)-1,DP)
!      Dv(18) = DuFine(IedgesAtElementFine(4,Ielf(4))-NVTfine+NMTfine)&
!             * REAL(2*IAND(ISHFT(ItwistFine(Ielf(4)),3),1)-1,DP)
!      Dv(19) = DuFine(IedgesAtElementFine(1,Ielf(4))-NVTfine+NMTfine)&
!             * REAL(2*IAND(      ItwistFine(Ielf(4))   ,1)-1,DP)
!      Dv(20) = DuFine(IedgesAtElementFine(4,Ielf(1))-NVTfine+NMTfine)&
!             * REAL(2*IAND(ISHFT(ItwistFine(Ielf(1)),3),1)-1,DP)
!      Dv(21) = DuFine(IedgesAtElementFine(2,Ielf(1))-NVTfine+NMTfine)&
!             * REAL(2*IAND(ISHFT(ItwistFine(Ielf(1)),1),1)-1,DP)
!      Dv(22) = DuFine(IedgesAtElementFine(2,Ielf(1))-NVTfine+NMTfine)&
!             * REAL(2*IAND(ISHFT(ItwistFine(Ielf(1)),1),1)-1,DP)
!      Dv(23) = DuFine(IedgesAtElementFine(2,Ielf(1))-NVTfine+NMTfine)&
!             * REAL(2*IAND(ISHFT(ItwistFine(Ielf(1)),1),1)-1,DP)
!      Dv(24) = DuFine(IedgesAtElementFine(2,Ielf(1))-NVTfine+NMTfine)&
!             * REAL(2*IAND(ISHFT(ItwistFine(Ielf(1)),1),1)-1,DP)
!      Dv(25) = DuFine(2*NMTfine + Ielf(1))
!      Dv(26) = DuFine(2*NMTfine + Ielf(2))
!      Dv(27) = DuFine(2*NMTfine + Ielf(3))
!      Dv(28) = DuFine(2*NMTfine + Ielf(4))
!      Dv(29) = DuFine(2*NMTfine + NELfine + Ielf(1))
!      Dv(30) = DuFine(2*NMTfine + NELfine + Ielf(2))
!      Dv(31) = DuFine(2*NMTfine + NELfine + Ielf(3))
!      Dv(32) = DuFine(2*NMTfine + NELfine + Ielf(4))
!      
!      ! And collect the values
!      ! Dofs 1 - 4 -> Int-Means on Edges
!      DuCoarse(Idf(1)) = DuCoarse(Idf(1)) + Dtn(1)*(Dv(1)+Dv(2)) &
!        + 0.375_DP*(Dv(9) - Dv(11)) - 0.25_DP*(Dv(10) + Dv(12)) &
!        + 0.125_DP*(Dtn(1)*(Dv(13) + Dv(14)) - Dtn(2)*(Dv(15) + Dv(16)) &
!        - Dtn(3)*(Dv(17) + Dv(18)) - Dtn(4)*(Dv(19) + Dv(20)) &
!        + 0.25_DP*(-Dv(21) + Dv(25) + Dv(26) - Dv(27) - Dv(28)))
!      DuCoarse(Idf(2)) = DuCoarse(Idf(2)) + Dtn(2)*(Dv(3)+Dv(4)) &
!        + 0.375_DP*(Dv(10) - Dv(12)) - 0.25_DP*(Dv(9) + Dv(11)) &
!        + 0.125_DP*(-Dtn(1)*(Dv(13) + Dv(14)) + Dtn(2)*(Dv(15) + Dv(16)) &
!        - Dtn(3)*(Dv(17) + Dv(18)) - Dtn(4)*(Dv(19) - Dv(20)) &
!        + 0.25_DP*(Dv(22) - Dv(25) + Dv(26) + Dv(27) - Dv(28)))
!      DuCoarse(Idf(3)) = DuCoarse(Idf(3)) + Dtn(3)*(Dv(5)+Dv(6)) &
!        + 0.375_DP*(-Dv(9) - Dv(11)) - 0.25_DP*(Dv(10) + Dv(12)) &
!        + 0.125_DP*(-Dtn(1)*(Dv(13) + Dv(14)) - Dtn(2)*(Dv(15) + Dv(16)) &
!        + Dtn(3)*(Dv(17) + Dv(18)) + Dtn(4)*(Dv(19) - Dv(20)) &
!        + 0.25_DP*(Dv(23) - Dv(25) - Dv(26) + Dv(27) + Dv(28)))
!      DuCoarse(Idf(4)) = DuCoarse(Idf(4)) + Dtn(4)*(Dv(7)+Dv(8)) &
!        + 0.375_DP*(Dv(12) - Dv(10)) - 0.25_DP*(Dv(9) + Dv(11)) &
!        + 0.125_DP*(-Dtn(1)*(Dv(13) + Dv(14)) - Dtn(2)*(Dv(15) + Dv(16)) &
!        - Dtn(3)*(Dv(17) + Dv(18)) - Dtn(4)*(Dv(19) + Dv(20)) &
!        + 0.25_DP*(-Dv(24) + Dv(25) - Dv(26) - Dv(27) + Dv(28)))
!      ! Dofs 5 - 8 -> Legendre-Int-Means on Edges
!      DuCoarse(Idf(5)) = DuCoarse(Idf(5)) + Dtf(1)*(-1.59375_DP*Dtn(1)*(Dv(1) - Dv(2)) &
!        - 0.09375_DP*(Dtn(2)*(Dv(3) - Dv(4)) + Dtn(3)*(Dv(5) - Dv(6) + Dv(17) - Dv(18)) &
!        + Dtn(4)*(Dv(7) + Dv(8) + Dv(19))) - 0.375_DP*(Dv(10) - Dv(12)) &
!        + 0.40625_DP*(Dtn(1)*(Dv(13) - Dv(14))) - 0.46875_DP*Dtn(2)*Dv(15) &
!        - 0.28125_DP*(Dtn(2)*Dv(16) + Dtn(4)*Dv(20)) - 0.125_DP*(Dv(22) - Dv(24)) &
!        + 0.1875_DP*(-Dv(25) + Dv(26) - Dv(27) + Dv(28)))
!      DuCoarse(Idf(6)) = DuCoarse(Idf(6)) + Dtf(2)*(-1.59375_DP*Dtn(2)*(Dv(3) - Dv(4)) &
!        - 0.09375_DP*(Dtn(1)*(Dv(1) - Dv(2)) + Dtn(3)*(Dv(5) - Dv(6)) &
!        + Dtn(4)*(Dv(7) - Dv(8))) - 0.375_DP*(Dv(11) - Dv(9)) &
!        + 0.40625_DP*(Dtn(2)*(Dv(15) - Dv(16)) - Dtn(4)*Dv(20)) &
!        + 0.46875_DP*(Dtn(1)*Dv(14) - Dtn(3)*Dv(17) - Dtn(4)*Dv(19))&
!        - 0.28125_DP*(Dtn(1)*Dv(13) + Dtn(3)*Dv(18)) - 0.125_DP*(Dv(23) - Dv(21)) &
!        + 0.1875_DP*(+Dv(25) - Dv(26) + Dv(27) - Dv(28)))
!      DuCoarse(Idf(7)) = DuCoarse(Idf(7)) + Dtf(3)*(- 1.59375_DP*Dtn(3)*(Dv(5) - Dv(6)) &
!        - 0.09375_DP*(Dtn(1)*(Dv(1) - Dv(2) + Dv(13) - Dv(13)) + Dtn(2)*(Dv(3) - Dv(4)) &
!        + Dtn(4)*(Dv(7) - Dv(8)))- 0.375_DP*(Dv(12) - Dv(9)) &
!        + 0.40625_DP*(Dtn(3)*(Dv(17) - Dv(18)) + Dtn(4)*Dv(19)) &
!        + 0.46875_DP*(Dtn(2)*Dv(16) + Dtn(4)*Dv(20)) + 0.28125_DP*(Dtn(2)*Dv(15)) &
!        + 0.125_DP*(Dv(22) + Dv(24)) + 0.1875_DP*(-Dv(25) + Dv(26) - Dv(27) + Dv(28)))
!      DuCoarse(Idf(8)) = DuCoarse(Idf(8)) + Dtf(4)*(- 1.59375_DP*Dtn(4)*(Dv(7) - Dv(8)) &
!        - 0.09375_DP*(Dtn(1)*(Dv(1) - Dv(2)) + Dtn(2)*(Dv(3) - Dv(4) + Dv(15) + Dv(16)) &
!        + Dtn(3)*(Dv(5) - Dv(6)) - Dtn(4)*Dv(20)) - 0.375_DP*(Dv(9) - Dv(11)) &
!        - 0.46875_DP*(Dtn(1)*Dv(13) - Dtn(3)*Dv(18))&
!        - 0.28125_DP*(Dtn(1)*Dv(14) - Dtn(3)*Dv(17) - Dtn(4)*Dv(19)) &
!        - 0.125_DP*(Dv(21) - Dv(23)) + 0.1875_DP*(+Dv(25) - Dv(26) + Dv(27) - Dv(28)))
!      
!      DuCoarse(Idf(9)) = DuCoarse(Idf(9)) + 1.5_DP*(Dv(9) + Dv(10) + Dv(11) + Dv(12)) &
!        + 0.25_DP*(Dtn(1)*(Dv(13) + Dv(14)) + Dtn(2)*(Dv(15) + Dv(16)) &
!                 + Dtn(3)*(Dv(17) + Dv(18)) + Dtn(4)*(Dv(19) + Dv(20)) &
!                 - Dv(21) - Dv(22) - Dv(23) + Dv(24))
!      
!      DuCoarse(Idf(10)) = DuCoarse(Idf(10)) &
!        - 6.25_DP*(Dtn(1)*(Dv(13) + Dv(14)) + Dtn(2)*(Dv(15) + Dv(16)) &
!                 + Dtn(3)*(Dv(17) + Dv(18)) + Dtn(4)*(Dv(19) + Dv(20))) &
!        - 3.125_DP*(Dv(21) + Dv(22) + Dv(23) - Dv(24)) &
!        + 0.0625_DP*(Dv(29) + Dv(30) + Dv(31) + Dv(32))
!        
!    END DO
!    
!    ! That's it
!
!  END SUBROUTINE

  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_initL2Proj (rprojection,r2Lvlmass,rmass,rlumpedMass,&
                               rvecTemp1,rvecTemp2)

!<description>
  ! Sets the matrices and vectors which are needed for L2-projection.
!</description>
  
!<input>
  ! The 2-Level-Mass matrix
  TYPE(t_matrixScalar), INTENT(IN) :: r2LvlMass

  ! The mass matrix of the fine grid
  TYPE(t_matrixScalar), INTENT(IN) :: rmass

  ! OPTIONAL: The lumped mass matrix of the fine grid. If not given, the
  ! lumped mass matrix is created from rmassFine.
  TYPE(t_matrixScalar), OPTIONAL, INTENT(IN) :: rlumpedMass

  ! OPTIONAL: Two temporary vectors that match the structure of the fine
  ! mesh mass matrix. The vectors must not share the same data array.
  TYPE(t_vectorScalar), OPTIONAL, INTENT(IN) :: rvecTemp1
  TYPE(t_vectorScalar), OPTIONAL, INTENT(IN) :: rvecTemp2
!</input>

!<inputoutput>
  ! The scalar projection structure for which the matrices are to be set.
  TYPE(t_interlevelProjectionScalar), INTENT(INOUT) :: rprojection 
!</inputout>

!</subroutine>

     ! This is an L2-projection
     rprojection%iprojType = MLP_PROJ_TYPE_L2

     ! Create shared copies of the mass matrices
     CALL lsyssc_duplicateMatrix(r2LvlMass, rprojection%rmatrix2LvlMass, &
                                 LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
     CALL lsyssc_duplicateMatrix(rmass, rprojection%rmatrixMass, &
                                 LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
     
     ! Transpose the 2-Level-Mass matrix
     CALL lsyssc_transposeMatrix(r2LvlMass, rprojection%rmatrix2LvlMassT,&
                                 LSYSSC_TR_VIRTUAL)
     
     ! Do we have the lumped fine grid mass matrix?
     IF (PRESENT(rlumpedMass)) THEN

       ! Create a shared copy of it then.
       CALL lsyssc_duplicateMatrix(rlumpedMass, rprojection%rlumpedMass, &
                                   LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

     ELSE
     
       ! Copy mass matrix
       CALL lsyssc_duplicateMatrix(rmass, rprojection%rlumpedMass, &
                                   LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
     
       ! And lump it
       CALL lsyssc_lumpMatrixScalar(rprojection%rlumpedMass,LSYSSC_LUMP_DIAG)
     
     END IF
     
     ! Do we have the temporary vectors?
     IF (PRESENT(rvecTemp1)) THEN
       ! Create a shared copy of it
       CALL lsyssc_duplicateVector(rvecTemp1,rprojection%rvectorTmp,&
                                   LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
     ELSE
       ! Create a vector based on the matrix
       CALL lsyssc_createVecIndMat(rmass, rprojection%rvectorTmp)
     END IF
     IF (PRESENT(rvecTemp2)) THEN
       ! Create a shared copy of it
       CALL lsyssc_duplicateVector(rvecTemp2,rprojection%rvectorDef,&
                                   LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
     ELSE
       ! Create a vector based on the matrix
       CALL lsyssc_createVecIndMat(rmass, rprojection%rvectorDef)
     END IF
     
     ! That's it
     
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_prolScalarL2 (rprojection, rcoarseVector, rfineVector)
  
!<description>
  ! Performs an L2-prolongation of a solution vector on the coarse grid
  ! to a solution vector of the fine grid.
!</description>
  
!<input>
  ! The t_interlevelProjectionScalar structure that configures the grid transfer
  TYPE(t_interlevelProjectionScalar), INTENT(IN) :: rprojection 

  ! Coarse grid vector
  TYPE(t_vectorScalar), INTENT(INOUT) :: rcoarseVector
!</input>

!<output>
  ! Fine grid vector
  TYPE(t_vectorScalar), INTENT(INOUT) :: rfineVector
!</output>
  
!</subroutine>

  INTEGER :: i
  REAL(DP) :: ddefInit, ddef
  TYPE(t_vectorScalar) :: rtmp, rdef
  LOGICAL :: bcheckDef
  
    ! Get temporary vectors
    rtmp = rprojection%rvectorTmp
    rdef = rprojection%rvectorDef
  
    ! Do we check the defect?
    bcheckDef = ((rprojection%depsRelL2 .GT. 0.0_DP) .AND. &
                 (rprojection%depsAbsL2 .GT. 0.0_DP))

    ! Clear fine grid vector
    CALL lsyssc_clearVector(rfineVector)
    
    ! Multiply coarse grid vector with 2-Level-Mass
    CALL lsyssc_scalarMatVec(rprojection%rmatrix2LvlMass, rcoarseVector,&
                             rtmp, 1.0_DP, 0.0_DP)
    
    ! Calculate initial defect
    CALL lsyssc_copyVector(rtmp, rdef)
    IF (bcheckDef) THEN
      ddefInit = lsyssc_vectorNorm(rdef, LINALG_NORML2)
      IF (ddefInit .LE. rprojection%depsAbsL2) RETURN
      IF (ddefInit .LE. SYS_EPSREAL) ddefInit = 1.0_DP
    END IF
    
    ! Start the defect correction
    DO i = 1, rprojection%imaxL2Iterations
    
      ! Multiply by the inverse of the lumped mass matrix:
      ! d := M_l^-1 d
      CALL lsyssc_invertedDiagMatVec (rprojection%rlumpedMass,&
                                      rdef,1.0_DP,rdef)

      ! Add to the main vector:  x = x + omega*d
      CALL lsyssc_vectorLinearComb (rdef,rfineVector,1.0_DP,1.0_DP)
      
      ! Set up the defect: d := b-Mx
      CALL lsyssc_copyVector (rtmp,rdef)
      CALL lsyssc_scalarMatVec (rprojection%rmatrixMass, rfineVector,&
                                rdef, -1.0_DP, 1.0_DP)
      
      IF (bcheckDef) THEN
        ddef = lsyssc_vectorNorm(rdef, LINALG_NORML2)
        
        ! Are we finished?
        IF ((ddef .LE. rprojection%depsAbsL2) .AND. (ddef/ddefInit .LE.&
            rprojection%depsRelL2)) EXIT
      
      END IF

    END DO
    
    ! That's it

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlprj_restScalarL2 (rprojection,rcoarseVector,rfineVector)
  
!<description>
  ! Performs an L2-restriction of a defect vector on the fine grid
  ! to a defect vector of the coarse grid.
!</description>
  
!<input>
  ! The t_interlevelProjectionScalar structure that configures the grid transfer
  TYPE(t_interlevelProjectionScalar), INTENT(IN) :: rprojection 

  ! Fine grid vector
  TYPE(t_vectorScalar), INTENT(INOUT) :: rfineVector
!</input>

!<output>
  ! Coarse grid vector
  TYPE(t_vectorScalar), INTENT(INOUT) :: rcoarseVector
!</output>
  
!</subroutine>

  INTEGER :: i
  REAL(DP) :: ddefInit, ddef
  TYPE(t_vectorScalar) :: rtmp, rdef
  LOGICAL :: bcheckDef
  
    ! Get temporary vectors
    rtmp = rprojection%rvectorTmp
    rdef = rprojection%rvectorDef
    
    ! Do we check the defect?
    bcheckDef = ((rprojection%depsRelL2 .GT. 0.0_DP) .AND. &
                 (rprojection%depsAbsL2 .GT. 0.0_DP))
  
    ! Clear temporary vector
    CALL lsyssc_clearVector(rtmp)
        
    ! Calculate initial defect
    CALL lsyssc_copyVector(rfineVector, rdef)
    
    IF (bcheckDef) THEN
      ddefInit = lsyssc_vectorNorm(rdef, LINALG_NORML2)
      IF (ddefInit .LE. rprojection%depsAbsL2) RETURN
      IF (ddefInit .LE. SYS_EPSREAL) ddefInit = 1.0_DP
    END IF
    
    ! Start the defect correction
    DO i = 1, rprojection%imaxL2Iterations
    
      ! Multiply by the inverse of the lumped mass matrix:
      ! d := M_l^-1 d
      CALL lsyssc_invertedDiagMatVec (rprojection%rlumpedMass,&
                                      rdef,1.0_DP,rdef)

      ! Add to the main vector:  x = x + omega*d
      CALL lsyssc_vectorLinearComb (rdef,rtmp,1.0_DP,1.0_DP)
      
      ! Set up the defect: d := b-Mx
      CALL lsyssc_copyVector (rfineVector,rdef)
      CALL lsyssc_scalarMatVec (rprojection%rmatrixMass, rtmp,&
                                rdef, -1.0_DP, 1.0_DP)
      
      IF (bcheckDef) THEN
        ddef = lsyssc_vectorNorm(rdef, LINALG_NORML2)
        
        ! Are we finished?
        IF ((ddef .LE. rprojection%depsAbsL2) .AND. (ddef/ddefInit .LE.&
            rprojection%depsRelL2)) EXIT
            
      END IF

    END DO

    ! Multiply temporary vector with transposed 2-Level-Mass
    CALL lsyssc_scalarMatVec(rprojection%rmatrix2LvlMassT, rtmp,&
                             rcoarseVector, 1.0_DP, 0.0_DP)
    
    ! That's it

  END SUBROUTINE

END MODULE
