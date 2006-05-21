!##############################################################################
!# ****************************************************************************
!# <name> multilevelprojection </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains routines for prolongation, restriction and interpolation
!# of solution vectors between different levels, for scalar systems as well
!# as for block systems.
!# </purpose>
!##############################################################################

MODULE multilevelprojection

  USE fsystem
  USE spatialdiscretisation
  USE linearsystemscalar
  USE linearsystemblock
  
  IMPLICIT NONE

!<types>
  
  !<typeblock>
  
  ! Describes for a scalar equation which type of prolongation, restriction
  ! or interpolation to be used. The information is (of course) element
  ! dependent.
  
  TYPE t_interlevelProjectionScalar
  
    ! Element type that should be assumed for the prolongation. Must fit
    ! to the element type of the discretisation concerning the DOF's:
    ! If element "EM30" is used e.g., prolongation for "EM30", "EM31",
    ! "E031" and "E030" is possible, although performance will suffer
    ! in using the "wrong" prolongation. An error will be displayed if
    ! the DOF's don't fit together at all, e.g. when using $Q_1$
    ! prolongation when using $Q_2$ for the discretisation.
    INTEGER                     :: ielementTypeProlongation = 0
    
    ! Order of the prolongation to use. 
    ! 0=use default prolongation (e.g. linear for $Q_1$, quadratic for 
    ! $Q_2$,...). Some elements allow to configure the type of the 
    ! prolongation, e.g. when $Q_2$ is used, apart from 0 the following 
    ! values are allowed:
    ! 1=linear prolongation of a once refined mesh
    ! 2=usual quadratic interpolation
    INTEGER                     :: iprolongationOrder = 0

    ! Element type that should be assumed for the restriction. Must fit
    ! to the element type of the discretisation concerning the DOF's:
    ! If element "EM30" is used e.g., restriction for "EM30", "EM31",
    ! "E031" and "E030" is possible, although performance will suffer
    ! in using the "wrong" restriction. An error will be displayed if
    ! the DOF's don't fit together at all, e.g. when using $Q_1$
    ! restriction when using $Q_2$ for the discretisation.
    INTEGER                     :: ielementTypeRestriction = 0

    ! Order of the restriction to use. 
    ! 0=use default restriction (e.g. linear for $Q_1$, quadratic for 
    ! $Q_2$,...). Some elements allow to configure the type of the 
    ! restriction, e.g. when $Q_2$ is used, apart from 0 the following 
    ! values are allowed:
    ! 1=linear restriction of a once refined mesh
    ! 2=usual quadratic restriction
    INTEGER                     :: irestrictionOrder = 0

    ! Element type that should be assumed for the interpolation of a
    ! solution to a lower level. Must fit to the element type of the 
    ! discretisation concerning the DOF's.
    ! An error will be displayed if the DOF's don't fit together at 
    ! all, e.g. when using $Q_1$ interpolation when using $Q_2$ 
    ! for the discretisation.
    INTEGER                     :: ielementTypeInterpolation = 0
    
    ! Order of the interpolation to use when interpolating a solution vector
    ! to a lower level.
    ! 0=use default interpolation (e.g. linear for $Q_1$, quadratic for 
    ! $Q_2$,...). Some elements allow to configure the type of the 
    ! interpolation, e.g. when $Q_2$ is used, apart from 0 the following 
    ! values are allowed:
    ! 1=linear interpolation of a once refined mesh
    ! 2=usual quadratic interpolation
    INTEGER                     :: iinterpolationOrder = 0
    
  END TYPE
  
  !</typeblock>

  !<typeblock>
  
  ! Contains a list of t_interlevelProjectionScalar structures that
  ! describes the projection of block vectors between different levels
  ! in the discretisation
  
  TYPE t_interlevelProjectionBlock
  
    ! A list of t_interlevelProjectionScalar structures for every
    ! equation in the discretisation.
    TYPE(t_interlevelProjectionScalar), DIMENSION(:),POINTER :: p_rscalarProjection => NULL()
  
  END TYPE
  
  !</typeblock>

!</types>

CONTAINS

!<subroutine>

  SUBROUTINE mlevprj_initProjection (rprojection,DelementList)
  
!<description>
  ! This subroutine initialises an t_interlevelProjectionBlock with default
  ! values for a given PDE systen. DelementList is a list of element identifiers
  ! for all the equations in the system. 
!</description>
  
!<input>
  ! An array of element identifiers for each equation in the solution vector
  ! of the PDE. E.g. for CC2D: /elem_Q1Tnonpar,elem_Q1Tnonpar,elem_P0/
  INTEGER, DIMENSION(:), INTENT(IN) :: DelementList 
!</input>
  
!<output>
  ! A t_interlevelProjectionBlock structure that will be filled with data
  ! about the projection.
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
!</output>
  
!</subroutine>

  ! local variables
  
  INTEGER :: i
  
  ! Allocate an array that contains a t_interlevelProjectionScalar structure
  ! for every equation and initialise it
  
  ALLOCATE(rprojection%p_rscalarProjection(SIZE(DelementList)))
  
  DO i=1,SIZE(DelementList)
    rprojection%p_rscalarProjection(i)%ielementTypeProlongation  = DelementList(i)
    rprojection%p_rscalarProjection(i)%ielementTypeRestriction   = DelementList(i)
    rprojection%p_rscalarProjection(i)%ielementTypeInterpolation = DelementList(i)
  END DO
    
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlevprj_doneProjection (rprojection)
  
!<description>
  ! Cleans up a t_interlevelProjectionBlock structure. All dynamically allocated
  ! memory is released.
!</description>
  
!<inputoutput>
  ! The t_interlevelProjectionBlock structure which is to be cleaned up.
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
!</inputoutput>
  
!</subroutine>

  IF (.NOT. ASSOCIATED(rprojection%p_rscalarProjection)) THEN
    PRINT *,'mlevprj_doneProjection: Warning: trying to deallocate unallocated memory'
    RETURN
  END IF
  
  ! Deallocate allocated memory
  DEALLOCATE(rprojection%p_rscalarProjection)
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlevprj_performProlongation (rprojection,DcoarseVector, &
                     DfineVector,rdiscrCoarse,rdiscrFine)
  
!<description>
  ! Performs a prolongation for a given block vector (i.e. a projection
  ! in the primal space where the solution lives). The vector
  ! DcoarseVector on a coarser grid is projected to the vector
  ! DfineVector on a finer grid. rdiscrCoarse and rdiscrFine
  ! represent the coarse and fine grid discretisation structures 
  ! (triangulations,...). rprojection configures how the grid transfer
  ! is performed.
!</description>
  
!<input>
  ! The t_interlevelProjectionBlock structure that configures the grid transfer
  TYPE(t_interlevelProjectionBlock), INTENT(IN) :: rprojection 
  
  ! Coarse grid vector
  TYPE(t_vectorBlock), INTENT(IN) :: DcoarseVector
  
  ! Coarse grid discretisation
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscrCoarse

  ! Fine grid discretisation
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscrFine
!</input>
  
!<inputoutput>
  ! Fine grid vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: DfineVector
!</inputoutput>
  
!</subroutine>

  ! Calls the correct prolongation routines for each block in the 
  ! discretisation...

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlevprj_performRestriction (rprojection,DcoarseVector, &
                     DfineVector,rdiscrCoarse,rdiscrFine)
  
!<description>
  ! Performs a restriction for a given block vector (i.e. a projection
  ! in the dual space where the RHS vector lives). The vector
  ! DfineVector on a finer grid is projected to the vector
  ! DcoarseVector on a coarser grid. rdiscrCoarse and rdiscrFine
  ! represent the coarse and fine grid discretisation structures 
  ! (triangulations,...). rprojection configures how the grid transfer
  ! is performed.
!</description>
  
!<input>
  ! The t_interlevelProjectionBlock structure that configures the grid transfer
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
  
  ! Fine grid vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: DfineVector

  ! Fine grid discretisation
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscrFine

  ! Coarse grid discretisation
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscrCoarse
!</input>
  
!<inputoutput>
  ! Coarse grid vector
  TYPE(t_vectorBlock), INTENT(IN) :: DcoarseVector
!</inputoutput>
  
!</subroutine>

  ! Calls the correct restriction routines for each block in the 
  ! discretisation...

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mlevprj_performInterpolation (rprojection,DcoarseVector, &
                     DfineVector,rdiscrCoarse,rdiscrFine)
  
!<description>
  ! Performs an interpolation for a given block vector (i.e. a projection
  ! in the primal space where the solution vector lives). The vector
  ! DfineVector on a finer grid is projected to the vector
  ! DcoarseVector on a coarser grid. rdiscrCoarse and rdiscrFine
  ! represent the coarse and fine grid discretisation structures 
  ! (triangulations,...). rprojection configures how the grid transfer
  ! is performed.
!</description>
  
!<input>
  
  ! The t_interlevelProjectionBlock structure that configures the grid transfer
  TYPE(t_interlevelProjectionBlock), INTENT(OUT) :: rprojection 
  
  ! Fine grid vector
  TYPE(t_vectorBlock), INTENT(IN) :: DfineVector

  ! Fine grid discretisation
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscrFine

  ! Coarse grid discretisation
  TYPE(t_spatialDiscretisation), INTENT(INOUT) :: rdiscrCoarse
!</input>
  
!<inputoutput>
  ! Coarse grid vector
  TYPE(t_vectorBlock), INTENT(IN) :: DcoarseVector
!</inputoutput>
  
!</subroutine>
  
  ! Calls the correct interpolation routines for each block in the 
  ! discretisation...

  END SUBROUTINE
  
END MODULE
