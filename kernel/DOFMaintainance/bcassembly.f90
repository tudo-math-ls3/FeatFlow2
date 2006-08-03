!##############################################################################
!# ****************************************************************************
!# <name> bcassembly </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to discretise analytically given boundary
!# conditions. Analytically given boundary conditions are 'discretised', i.e.
!# a discrete version (realised by the structure t_discreteBCEntry in the case
!# of BC's on the real boundary and by the structure t_discreteFBCEntry in 
!# the case of fictitious boundary) is calculated.
!# This structure is used during the solution process to impose the boundary
!# conditions into the solution vector. Therefore, this module contains
!# the bridge how to put analytic boundary conditions into a discrete vector.
!#
!# The module works in tight relationship to the module 'vectorfilters'.
!# While bcassembly provides the functionality to *create* the structures,
!# the module 'vectorfilters' contains routines to *apply* the structure
!# to a given vector.
!# (This separation is necessary to prevent circular dependencies!)
!#
!# The following routines can be found here:
!#
!# 1.) bcasm_discretiseBC
!#     -> Discretises the analytic boundary definitions in a discretisation
!#        structure on the real bondary
!#
!# 2.) bcasm_releaseDiscreteBC
!#     -> Cleans up a structure with discrete boundary conditions and
!#        releases all memory allocated by it.
!#
!# 1.) bcasm_discretiseFBC
!#     -> Discretises the analytic boundary definitions of fictitious boundary
!#        objects in a discretisation structure
!#
!# 2.) bcasm_releaseDiscreteFBC
!#     -> Cleans up a structure with discrete boundary conditions of fictitious
!#        boundary objects and releases all memory allocated by it.
!#
!# The actual worker routines called by the above general ones are:
!#
!# 1.) bcasm_discrBCDirichlet
!#     -> Discretise Dirichlet boundary conditions on the real bondary
!#
!# 2.) bcasm_releaseDirichlet
!#     -> Release discrete Dirichlet boundary conditions  on the real bondary
!#
!# There are some routines that might be useful for finding out whether
!# edges/vertices exist in a boundary condition segment or not.
!#
!# 1.) bcasm_getVertInBCregion
!#     -> Returns indices in the boundary vertex array which vertices
!#        lie in a boundary region
!#
!# 1.) bcasm_getEdgesInBCregion
!#     -> Returns indices in the boundary vertex array which edges
!#        lie in a boundary region
!#
!# </purpose>
!##############################################################################

MODULE bcassembly

  USE fsystem
  USE boundarycondition
  USE discretebc
  USE discretefbc
  USE collection
  USE triangulation
  USE dofmapping
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Complexity of the discretised BC's">

  ! Discretise BC's for implementing them into a defect vector
  INTEGER, PARAMETER :: BCASM_DISCFORDEF = 2**0

  ! Discretise BC's for implementing them into a solution vector
  INTEGER, PARAMETER :: BCASM_DISCFORSOL = 2**1

  ! Discretise BC's for implementing them into a RHS vector
  INTEGER, PARAMETER :: BCASM_DISCFORRHS = 2**2

  ! Discretise BC's for implementing them into a matrix
  INTEGER, PARAMETER :: BCASM_DISCFORMAT = 2**3

  ! Discretise BC's for implementing them into matrix and defect vector
  INTEGER, PARAMETER :: BCASM_DISCFORDEFMAT = BCASM_DISCFORDEF + BCASM_DISCFORMAT
  
  ! Discretise BC's for implementing them into everything
  INTEGER, PARAMETER :: BCASM_DISCFORALL = BCASM_DISCFORDEF + BCASM_DISCFORSOL + &
                                           BCASM_DISCFORRHS + BCASM_DISCFORMAT

!</constantblock>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of entries to handle simultaneously (-> number of points or edges,
  ! depending on the situation what to discretise)
  INTEGER, PARAMETER :: FBCASM_MAXSIM   = 1000
  
!</constantblock>

!</constants>

CONTAINS

! *****************************************************************************
! Support for boundary conditions on real boundary
! *****************************************************************************

!<subroutine>

  SUBROUTINE bcasm_discretiseBC (rblockDiscretisation,p_rdiscreteBC,bforceRebuild, &
                                 fgetBoundaryValues,rcollection,casmComplexity)
  
!<description>
  ! This routine discretises an analytic definition of boundary conditions.
  ! The definition of the boundary conditions is taken from the discretisation
  ! structure rspatialDiscretisation. The discrete version is build up in
  ! p_rdiscreteBC. If p_rdiscreteBC is NULL(), a new structure is created,
  ! otherwise the old structure is updated (or even destroyed and recreated if
  ! necessary).
!</description>

!<input>
  
  ! The block discretisation structure of the underlying PDE. The boundary
  ! conditions inside of this structure are discretised.
  TYPE(t_blockDiscretisation), INTENT(IN) :: rblockDiscretisation
  
  ! Can be set to TRUE to force a complete rebuild of the rdiscreteBC structure.
  ! Normally, the structure is completely set up only in the first call
  ! or is there is a massive change in the boundary conditions (e.g. the
  ! number change)
  ! In later calls, the structure is recomputed only in those parts of the boundary,
  ! which have the t_bcRegion\%bisstatic flag set FALSE (the standard value) -
  ! except there , in which case the whole structure is rebuild.
  ! By setting bforceRebuild to TRUE, one can enforce a complete
  ! rebuild of the structure, independent of which regions are marked
  ! as static.
  LOGICAL                                   :: bforceRebuild
  
  ! A callback function that calculates values on the boundary.
  ! Is declared in the interface include file 'intf_bcassembly.inc'.
  INCLUDE 'intf_bcassembly.inc'
  
  ! OPTIONAL: A collection structure to inform the callback function with
  ! additional information. Can undefined if there is no
  ! information to pass.
  TYPE(t_collection), TARGET, OPTIONAL :: rcollection
  
  ! OPTIONAL: A combination of BCASM_DISCFORxxx constants that specify
  ! the complexity of the discretisation that is to perform. This allows to
  ! discretise only parts of the BC's, e.g. only setting up those
  ! information that are necessary for filtering defect vectors.
  ! If not specified, BCASM_DISCFORALL is assumed, i.e. the resulting
  ! boundary conditions can be used for everything.
  INTEGER(I32), INTENT(IN), OPTIONAL :: casmComplexity
!</input>

!<inputoutput>
  ! A discretised version of the analytic boundary conditions.
  ! This is a pointer to a t_discreteBC structures, 
  ! representing the boundary discretised in a discretisation- dependent way.
  ! If this pointer points to NULL(), a complete new structure is set up.
  TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: icurrentRegion, ccompl
  TYPE(t_bcRegion), POINTER :: p_rbcRegion
  LOGICAL :: bbuildAll
  TYPE(t_collection), POINTER :: p_rcoll
  
  ! Pointer to the boundary condition object
  TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions
  
  ! For quicker access:
  p_rboundaryConditions => rblockDiscretisation%p_rboundaryConditions
  
  ! We replace the optional parameter by NULL() if it does not exist
  IF (PRESENT(rcollection)) THEN
    p_rcoll => rcollection
  ELSE
    p_rcoll => NULL()
  END IF
  
  IF (.NOT. ASSOCIATED(p_rboundaryConditions)) THEN
    PRINT *,'Warning in bcasm_discretiseBC: Boundary conditions not associated!'
    RETURN
  END IF
  
  ! Default target for the discretisation if everything.
  IF (PRESENT(casmComplexity)) THEN
    ccompl = casmComplexity
  ELSE
    ccompl = BCASM_DISCFORALL
  END IF
  
  bbuildAll = bforceRebuild
  
  ! Is there a structure we can work with?
  IF (.NOT. ASSOCIATED(p_rdiscreteBC)) THEN
    ! Create a new structure on the heap.
    ALLOCATE(p_rdiscreteBC)
  END IF
  
  ! Now we work with the array in the structure. Does it exist?
  
  IF (ASSOCIATED(p_rdiscreteBC%p_RdiscBCList)) THEN
  
    ! Is there a massive change or do we have o rebuild everything?
    IF ((p_rboundaryConditions%iregionCount .NE. SIZE(p_rdiscreteBC%p_RdiscBCList)) &
        .OR. bforceRebuild) THEN
        
      ! Oh, we have to destroy the structure completely.
      ! Don't remove the structure itself from the heap - we want to
      ! fill it with data immediately afterwards!
      CALL bcasm_releaseDiscreteBC (p_rdiscreteBC,.TRUE.)

      ! Allocate a new structure array with iregionCount 
      ! entries for all the boundary conditions.
      ALLOCATE(p_rdiscreteBC%p_RdiscBCList(p_rboundaryConditions%iregionCount))
      
      bbuildAll = .TRUE.
    ELSE
      ! Otherwise, release only those information belonging to 
      ! non-static boundary regions
      DO icurrentRegion = 1,p_rboundaryConditions%iregionCount
        
        ! Get a pointer to it so we can deal with it more easily
        p_rbcRegion => p_rboundaryConditions%p_Rregions(icurrentRegion)
      
        IF (.NOT. p_rbcRegion%bisStatic) THEN
        
          ! Release all allocated information to this boundary region
          SELECT CASE (p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%itype)

          CASE (DISCBC_TPDIRICHLET)

            ! Discrete Dirichlet boundary conditions. Release the old structure.
            CALL bcasm_releaseDirichlet( &
                 p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%rdirichletBCs)

          CASE (DISCBC_TPPRESSUREDROP)

            ! Discrete Dirichlet boundary conditions. Release the old structure.
            CALL bcasm_releasePressureDrop( &
                 p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%rpressureDropBCs)

          END SELECT

          ! BC released, indicate this
          p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%itype = DISCBC_TPUNDEFINED
          
        END IF
      
      END DO  
    END IF
    
  ELSE
  
    ! Allocate a structure array with iregionCount entries
    ! for all the boundary conditions.
    ALLOCATE(p_rdiscreteBC%p_RdiscBCList(p_rboundaryConditions%iregionCount))
    bbuildAll = .TRUE.
  
  END IF
  
  ! Loop through the regions on the boundary
  DO icurrentRegion = 1,p_rboundaryConditions%iregionCount
    
    ! Get a pointer to it so we can deal with it more easily
    p_rbcRegion => p_rboundaryConditions%p_Rregions(icurrentRegion)
  
    ! Do we have to process this region?
    IF (bforceRebuild .OR. .NOT. p_rbcRegion%bisStatic) THEN
    
      ! Ok, let's go...
      ! What for BC do we have here?
      SELECT CASE (p_rbcRegion%ctype)
      
      CASE (BC_DIRICHLET)
        CALL bcasm_discrBCDirichlet (rblockDiscretisation, &
                   p_rbcRegion, p_rdiscreteBC%p_RdiscBCList(icurrentRegion), &
                   ccompl,fgetBoundaryValues,p_rcoll)
                   
      CASE (BC_PRESSUREDROP)

        CALL bcasm_discrBCpressureDrop (rblockDiscretisation, &
                   p_rbcRegion, p_rdiscreteBC%p_RdiscBCList(icurrentRegion), &
                   ccompl,fgetBoundaryValues,p_rcoll)
        
      END SELECT
        
    END IF
  
  END DO  

  END SUBROUTINE
      
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcasm_releaseDiscreteBC (p_rdiscreteBC,bkeepBCStructure)
  
!<description>
  ! This routine cleans up an array of t_discreteBCEntry describing discrete
  ! boundary conditions. All allocated memory is released. The structure
  ! p_rdiscreteBC is released from the heap as well.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the structure p_rdiscreteBC is not released
  ! from memory. If set to FALSE or not existent (the usual setting), the 
  ! structure p_rdiscreteBC will also be removed from the heap after 
  ! cleaning up.
  LOGICAL, INTENT(IN), OPTIONAL :: bkeepBCStructure
!</input>

!<inputoutput>
  ! A pointer to discretised boundary conditions. All memory allocated
  ! in and by this array is released from the heap. 
  TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: icurrentRegion
  
  IF (.NOT. ASSOCIATED(p_rdiscreteBC)) THEN
    PRINT *,'Warning in bcasm_releaseDiscreteBC: Nothing to cleanup!'
    RETURN
  END IF
  
  ! Destroy the content of the structure completely!
  IF (ASSOCIATED(p_rdiscreteBC%p_RdiscBCList)) THEN
  
    ! Destroy all the substructures in the array.
    DO icurrentRegion = 1,SIZE(p_rdiscreteBC%p_RdiscBCList)  
      
      ! Release all allocated information to this boundary region
      SELECT CASE (p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%itype)
      CASE (DISCBC_TPDIRICHLET)
        ! Discrete Dirichlet boundary conditions. Release the old structure.
        CALL bcasm_releaseDirichlet(&
          p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%rdirichletBCs)

      CASE (DISCBC_TPPRESSUREDROP)
        ! Discrete pressure drop boundary conditions. Release the old structure.
        CALL bcasm_releasePressureDrop(&
             p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%rpressureDropBCs)

      END SELECT

      ! BC released, indicate this
      p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%itype = DISCBC_TPUNDEFINED
      
    END DO  
    
    ! Release the array itself.
    DEALLOCATE(p_rdiscreteBC%p_RdiscBCList)
    
  END IF
  
  ! Deallocate the structure (if we are allowed to), finish.
  IF (.NOT. PRESENT(bkeepBCStructure)) THEN
    DEALLOCATE(p_rdiscreteBC)
  ELSE
    IF (.NOT. bkeepBCStructure) THEN
      DEALLOCATE(p_rdiscreteBC)
    END IF
  END IF

  END SUBROUTINE
      
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE bcasm_getVertInBCregion (rtriangulation,rboundary,rregion, &
                                      IminIndex,ImaxIndex,icount)
  
!<description>
  ! This routine receives a boundary region rregion describing a part on the
  ! boundary. According to the triangulation, this boundary region contains
  ! some vertices and edges on the boundary. The routine now figures out,
  ! which index in the IverticesAtBoundary in the triangulation structure 
  ! of the vertices that are inside of this boundary region.
  !
  ! Each boundary region is simply connected, but they may cross the
  ! maximum parameter value. Therefore, we find wither one or two 'sets'
  ! of indices in the IverticesAtBoundary array:
  ! 1.) One simly connected index set:
  !     ..... Vertex Vertex Vertex ......
  !     Here is icount=1 and IminIndex(1) / ImaxIndex(1) gives the
  !     first/last index in IverticesAtBoundary
  ! 2.) Two sets crossing the maximum parameter value:
  !     Vertex Vertex Vertex ........... Vertex Vertex Vertex
  !     Here is icount=2, IminIndex(2) / ImaxIndex(2) gives the
  !     first/last vertex of the first set (which was found at first when
  !     starting to search in the given interval) and IminIndex(1) / ImaxIndex(1)
  !     the first/last vertex if the 2nd set.
!</description>

!<input>
  ! The triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! The description of the domain boundary
  TYPE(t_boundary), INTENT(IN) :: rboundary
  
  ! A boundary region structure describing a part of a boundary
  TYPE(t_boundaryRegion), INTENT(IN) :: rregion
!</input>

!<output>
  ! The index of the first vertex in IverticesAtBoundary belonging to the
  ! bonudary region rregion. = 0, if no vertices belong to the given region.
  INTEGER, DIMENSION(2), INTENT(OUT) :: IminIndex

  ! The index of the last vertex in IverticesAtBoundary belonging to the
  ! bonudary region rregion. = -1, if no vertices belong to the given region.
  INTEGER, DIMENSION(2), INTENT(OUT) :: ImaxIndex
  
  ! Number of index sets in IverticesAtBoundary belonging
  ! to rregion. IminIndex(1)..ImaxIndex(1) is the first set,...
  ! IminIndex(icount)..ImaxIndex(icount) is the last set.
  INTEGER, INTENT(OUT) :: icount
!</output>

!</subroutine>

  ! local variables
  LOGICAL binside,bfinish
  REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
  REAL(DP), DIMENSION(2) :: dbegin,dend
  REAL(DP) :: dmaxpar
  INTEGER(I32), DIMENSION(:), POINTER :: p_IboundaryCpIdx
  INTEGER :: i, ifoundRegions
  INTEGER, DIMENSION(2) :: Iindex
  
  ! Get the parameter value array from the triangulation
  CALL storage_getbase_double (rtriangulation%h_DvertexParameterValue, &
                               p_DvertexParameterValue)
                               
  ! Get the array describing the start/end of each boundary component
  CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx, &
                            p_IboundaryCpIdx)
  
  ! Get the real min/max parameter value of the boundary region
  dmaxpar = boundary_dgetMaxParVal(rboundary, rregion%iboundCompIdx)
  dbegin(1) = rregion%dminParam
  dend(1) = rregion%dmaxParam
  IF (dend(1)-dbegin(1) .GE. dmaxpar) THEN
    ! Occupy everything
    icount = 1
    dbegin(1) = 0.0_DP
    dend(1) = dmaxpar
    dbegin(2) = 0.0_DP
    dend(2) = -1.0_DP
  ELSE
    ! Calculate the real min/max
    dbegin(1) = MOD(rregion%dminParam,dmaxpar)
    dend(1) = MOD(rregion%dmaxParam,dmaxpar)
    icount = 1
    IF (dend(1) .LT. dbegin(1)) THEN
      ! probably two sets
      dbegin(2) = dend(1)
      dend(2) = dmaxpar
      dbegin(1) = 0.0_DP
      icount = 2
    END IF
  END IF
  
  ! All regions in [dbegin(.),dend(.)] are now in asecnding order.
  ! icount indicates the maximum number of regions we expect to have
  ! vertices inside.
  ! Remark: dbegin/dend are not used anymore here - maybe in a later
  ! version if necessary for some reasons...
  !
  ! For now, use a simple linear search. 
  ! The boundary region structure tells us which boundary component
  ! to use; loop through all vertices there.
  
  ifoundRegions = 0
  IminIndex = 0
  ImaxIndex = -1
  binside = .FALSE.
  bfinish = .FALSE.
  
  ! Loop through all vertices in the current boundary component
  
  DO I = p_IboundaryCpIdx(rregion%iboundCompIdx), &
         p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
         
    ! Check if the vertex is inside the region.
    IF (boundary_isInRegion (rregion,rregion%iboundCompIdx,p_DvertexParameterValue(I))) THEN
      ! We are inside
      IF (.NOT. binside) THEN
        ! We are inside for the first time
        binside = .TRUE.
        IminIndex(ifoundRegions+1) = I
      END IF
    ELSE
      ! We are outside - for the first time?
      ! If yes, quit the loop.
      IF (binside) THEN
        binside = .FALSE.
        ImaxIndex(ifoundRegions+1) = I-1
        
        ! We completed the current region successfully
        ifoundRegions = ifoundRegions + 1
        
        ! Decrement icount. If it's still > 0, there's another region
        ! in dbegin/dend that may contain points
        icount = icount - 1
        
        IF (icount .LE. 0) THEN
          ! Finish that, we quit the search here as no more regions
          ! are expected.
          icount = ifoundRegions
          RETURN
        END IF
      END IF
    END IF
         
  END DO
  
  ! The loop is completed. Question: Are we still inside a region or not?
  ! If yes...
  
  IF (binside) THEN
    ! Save the last vertex number
    ImaxIndex(ifoundRegions+1) = p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
    ! Complete the region and finish.
    icount = ifoundRegions + 1
  ELSE 
    ! No we aren't. So we were awaiting points for another part of the region
    ! that never came! Reset the last index pair and quit.
    IminIndex(ifoundRegions+1) = 0
    ImaxIndex(ifoundRegions+1) = -1
    icount = ifoundRegions
  END IF
  
  IF (icount .EQ. 2) THEN
    ! We found the intervals in the 'wrong order'. When going through the boundary
    ! starting at a parameter value x, we would first find the 2nd segment
    ! x <= y..TMAX, then the first one 0..z. on the other hand set up IminIndex/
    ! ImaxIndex in the other order. So exchange IxxxIndex(1) with IxxxIndex(2)
    ! to get the right order again.
    Iindex(1) = IminIndex(2)
    IminIndex(2) = IminIndex(1)
    IminIndex(1) = Iindex(1)

    Iindex(2) = ImaxIndex(2)
    ImaxIndex(2) = ImaxIndex(1)
    ImaxIndex(1) = Iindex(2)
  END IF
  
  END SUBROUTINE
    
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE bcasm_getEdgesInBCregion (rtriangulation,rboundary,rregion, &
                                      IminIndex,ImaxIndex,icount)
  
!<description>
  ! This routine receives a boundary region rregion describing a part on the
  ! boundary. According to the triangulation, this boundary region contains
  ! some vertices and edges on the boundary. The routine now figures out,
  ! which  index in the IverticesAtBoundary in the triangulation structure 
  ! of the vertices that are inside of this boundary region.
  !
  ! Each boundary region is simply connected, but they may cross the
  ! maximum parameter value. Therefore, we find wither one or two 'sets'
  ! of indices in the IverticesAtBoundary array:
  ! 1.) One simly connected index set:
  !     ..... Edge Edge Edge......
  !     Here is icount=1 and IminIndex(1) / ImaxIndex(1) gives the
  !     first/last index in IverticesAtBoundary
  ! 2.) Two sets crossing the maximum parameter value:
  !     Edge Edge Edge ........... Edge Edge Edge
  !     Here is icount=2, IminIndex(2) / ImaxIndex(2) gives the
  !     first/last vertex of the first set (which was found at first when
  !     starting to search in the given interval) and IminIndex(1) / ImaxIndex(1)
  !     the first/last vertex if the 2nd set.
!</description>

!<input>
  ! The triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! The description of the domain boundary
  TYPE(t_boundary), INTENT(IN) :: rboundary
  
  ! A boundary region structure describing a part of a boundary
  TYPE(t_boundaryRegion), INTENT(IN) :: rregion
!</input>

!<output>
  ! The index of the first edge in IedgesAtBoundary belonging to the
  ! bonudary region rregion. = 0, if no edges belong to the given region.
  INTEGER, DIMENSION(2), INTENT(OUT) :: IminIndex

  ! The index of the last edges in IedgesAtBoundary belonging to the
  ! bonudary region rregion. = -1, if no edges belong to the given region.
  INTEGER, DIMENSION(2), INTENT(OUT) :: ImaxIndex
  
  ! Number of index sets in IedgesAtBoundary belonging
  ! to rregion. IminIndex(1)..ImaxIndex(1) is the first set,...
  ! IminIndex(icount)..ImaxIndex(icount) is the last set.
  INTEGER, INTENT(OUT) :: icount
!</output>

!</subroutine>

  ! local variables
  LOGICAL binside,bfinish
  REAL(DP), DIMENSION(:), POINTER :: p_DedgeParameterValue
  REAL(DP), DIMENSION(2) :: dbegin,dend
  REAL(DP) :: dmaxpar
  INTEGER(I32), DIMENSION(:), POINTER :: p_IboundaryCpIdx
  INTEGER :: i, ifoundRegions
  INTEGER, DIMENSION(2) :: Iindex
  
  ! Get the parameter value array from the triangulation
  CALL storage_getbase_double (rtriangulation%h_DedgeParameterValue, &
                               p_DedgeParameterValue)
                               
  ! Get the array describing the start/end of each boundary component
  CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx, &
                            p_IboundaryCpIdx)
  
  ! Get the real min/max parameter value of the boundary region
  dmaxpar = boundary_dgetMaxParVal(rboundary, rregion%iboundCompIdx)
  dbegin(1) = rregion%dminParam
  dend(1) = rregion%dmaxParam
  IF (dend(1)-dbegin(1) .GE. dmaxpar) THEN
    ! Occupy everything
    icount = 1
    dbegin(1) = 0.0_DP
    dend(1) = dmaxpar
    dbegin(2) = 0.0_DP
    dend(2) = -1.0_DP
  ELSE
    ! Calculate the real min/max
    dbegin(1) = MOD(rregion%dminParam,dmaxpar)
    dend(1) = MOD(rregion%dmaxParam,dmaxpar)
    icount = 1
    IF (dend(1) .LT. dbegin(1)) THEN
      ! probably two sets
      dbegin(2) = dend(1)
      dend(2) = dmaxpar
      dbegin(1) = 0.0_DP
      icount = 2
    END IF
  END IF
  
  ! All regions in [dbegin(.),dend(.)] are now in asecnding order.
  ! icount indicates the maximum number of regions we expect to have
  ! vertices inside.
  ! Remark: dbegin/dend are not used anymore here - maybe in a later
  ! version if necessary for some reasons...
  !
  ! For now, use a simple linear search. 
  ! The boundary region structure tells us which boundary component
  ! to use; loop through all vertices there.
  
  ifoundRegions = 0
  IminIndex = 0
  ImaxIndex = -1
  binside = .FALSE.
  bfinish = .FALSE.
  
  ! Loop through all vertices in the current boundary component
  
  DO I = p_IboundaryCpIdx(rregion%iboundCompIdx), &
         p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
         
    ! Check if the edge is inside the region.
    IF (boundary_isInRegion (rregion,rregion%iboundCompIdx,p_DedgeParameterValue(I))) THEN
      ! We are inside
      IF (.NOT. binside) THEN
        ! We are inside for the first time
        binside = .TRUE.
        IminIndex(ifoundRegions+1) = I
      END IF
    ELSE
      ! We are outside - for the first time?
      ! If yes, quit the loop.
      IF (binside) THEN
        binside = .FALSE.
        ImaxIndex(ifoundRegions+1) = I-1
        
        ! We completed the current region successfully
        ifoundRegions = ifoundRegions + 1
        
        ! Decrement icount. If it's still > 0, there's another region
        ! in dbegin/dend that may contain points
        icount = icount - 1
        
        IF (icount .LE. 0) THEN
          ! Finish that, we quit the search here as no more regions
          ! are expected.
          icount = ifoundRegions
          RETURN
        END IF
      END IF
    END IF
         
  END DO
  
  ! The loop is completed. Question: Are we still inside a region or not?
  ! If yes...
  
  IF (binside) THEN
    ! Save the last edge number
    ImaxIndex(ifoundRegions+1) = p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
    ! Complete the region and finish.
    icount = ifoundRegions + 1
  ELSE 
    ! No we aren't. So we were awaiting points for another part of the region
    ! that never came! Reset the last index pair and quit.
    IminIndex(ifoundRegions+1) = 0
    ImaxIndex(ifoundRegions+1) = -1
    icount = ifoundRegions
  END IF

  IF (icount .EQ. 2) THEN
    ! We found the intervals in the 'wrong order'. When going through the boundary
    ! starting at a parameter value x, we would first find the 2nd segment
    ! x <= y..TMAX, then the first one 0..z. on the other hand set up IminIndex/
    ! ImaxIndex in the other order. So exchange IxxxIndex(1) with IxxxIndex(2)
    ! to get the right order again.
    Iindex(1) = IminIndex(2)
    IminIndex(2) = IminIndex(1)
    IminIndex(1) = Iindex(1)

    Iindex(2) = ImaxIndex(2)
    ImaxIndex(2) = ImaxIndex(1)
    ImaxIndex(1) = Iindex(2)
  END IF

  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcasm_discrBCDirichlet (rblockDiscretisation, &
                                     rbcRegion, rdiscreteBC, casmComplexity, &
                                     fgetBoundaryValues,p_rcollection)
  
!<description>
  ! Creates a discrete version of Dirichlet boundary conditions.
  ! rbcRegion describes the region which is to be discretised. The discretised
  ! boundary conditions are created in rdiscreteBC, which is assumed
  ! to be undefined when entering this routine.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation. The boundary
  ! conditions inside of this structure are discretised.
  TYPE(t_blockDiscretisation), INTENT(IN), TARGET :: rblockDiscretisation

  ! The BC region which is to be discrised
  TYPE(t_bcRegion), INTENT(IN) :: rbcRegion

  ! A callback function that calculates values on the boundary.
  ! Is declared in the interface include file 'intf_bcassembly.inc'.
  INCLUDE 'intf_bcassembly.inc'
  
  ! A collection structure to inform the callback function with
  ! additional information. Can be NULL() if there is no information to pass.
  TYPE(t_collection), POINTER :: p_rcollection

  ! A combination of BCASM_DISCFORxxx constants that specify
  ! the complexity of the discretisation that is to perform. This allows to
  ! discretise only parts of the BC's, e.g. only setting up those
  ! information that are necessary for filtering defect vectors.
  ! If not specified, BCASM_DISCFORALL is assumed, i.e. the resulting
  ! boundary conditions can be used for everything.
  INTEGER(I32), INTENT(IN) :: casmComplexity
!</input>  

!<output>
  ! This structure receives the result of the discretisation of rbcRegion.
  ! When entering the routine, the content of this structure is undefined,
  ! all pointers are invalid. The routine fills everything with appropriate
  ! data.
  TYPE(t_discreteBCEntry), INTENT(OUT), TARGET :: rdiscreteBC
!</output>

!</subroutine>

  ! local variables
  INTEGER, DIMENSION(2) :: IminVertex,ImaxVertex,IminEdge,ImaxEdge,Iminidx,Imaxidx
  INTEGER :: i,ilocalEdge,ieltype,ielement,icount,icount2,ipart,j,icomponent
  INTEGER :: icountmin,icountmax
  INTEGER(I32) :: iedge,ipoint1,ipoint2,NVT
  INTEGER, DIMENSION(1) :: Icomponents
  TYPE(t_discreteBCDirichlet),POINTER         :: p_rdirichletBCs
  TYPE(t_triangulation), POINTER              :: p_rtriangulation
  TYPE(t_spatialDiscretisation), POINTER      :: p_rspatialDiscretisation
  INTEGER(I32), DIMENSION(:), POINTER         :: p_IelementsAtBoundary,p_IverticesAtBoundary
  INTEGER(I32), DIMENSION(:), POINTER         :: p_IedgesAtBoundary
  INTEGER(I32), DIMENSION(:), POINTER         :: p_ItrialElements
  INTEGER, DIMENSION(:,:), ALLOCATABLE        :: Idofs
  REAL(DP), DIMENSION(:,:), ALLOCATABLE       :: DdofValue
  REAL(DP), DIMENSION(:), POINTER             :: p_DedgeParameterValue,p_DvertexParameterValue
  INTEGER(I32), DIMENSION(:,:), POINTER       :: p_IverticesAtEdge,p_IedgesAtElement
  REAL(DP), DIMENSION(:), POINTER             :: p_IdirichletValues
  INTEGER(I32), DIMENSION(:), POINTER         :: p_IdirichletDOFs
  
  REAL(DP), DIMENSION(DER_MAXNDER)            :: Dvalues
  
  ! Which component is to be discretised?
  icomponent = rbcRegion%Iequations(1)
  p_rspatialDiscretisation => rblockDiscretisation%RspatialDiscretisation(icomponent)

  ! For easier access:
  p_rtriangulation => p_rspatialDiscretisation%p_rtriangulation
  CALL storage_getbase_int(p_rtriangulation%h_IelementsAtBoundary,p_IelementsAtBoundary)
  CALL storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
  CALL storage_getbase_int(p_rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)
  CALL storage_getbase_double(p_rtriangulation%h_DedgeParameterValue,p_DedgeParameterValue)
  CALL storage_getbase_double(p_rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
  CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
  CALL storage_getbase_int2D(p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
  NVT = p_rtriangulation%NVT

  IF (p_rspatialDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) THEN
    ! Every element can be of different type.
    CALL storage_getbase_int(p_rspatialDiscretisation%h_ItrialElements,p_ItrialElements)
  ELSE
    ! All elements are of the samne type. Get it in advance.
    ieltype = p_rspatialDiscretisation%RelementDistribution(1)%itrialElement
  END IF
  
  ! We have Dirichlet boundary conditions
  rdiscreteBC%itype = DISCBC_TPDIRICHLET
  
  ! Connect to the boundary condition structure
  rdiscreteBC%p_rboundaryConditions => rblockDiscretisation%p_rboundaryConditions
  
  ! Fill the structure for discrete Dirichlet BC's in the
  ! t_discreteBCEntry structure
  p_rdirichletBCs => rdiscreteBC%rdirichletBCs
  
  p_rdirichletBCs%icomponent = icomponent
  Icomponents(1) = icomponent
  
  ! We have to deal with all DOF's on the boundary. This is highly element
  ! dependent and therefore a little bit tricky :(
  !
  ! As we are in 2D, we can use parameter values at first to figure out,
  ! which points and which edges are on the boundary.
  ! What we have is a boundary segment. Now ask the boundary-index routine
  ! to give us the minimum and maximum index of the vertices and edges on the
  ! bondary that belong to this boundary segment.
  
  CALL bcasm_getVertInBCregion (p_rtriangulation,p_rspatialDiscretisation%p_rdomain, &
                                rbcRegion%rboundaryRegion, &
                                IminVertex,ImaxVertex,icount)
  CALL bcasm_getEdgesInBCregion (p_rtriangulation,p_rspatialDiscretisation%p_rdomain, &
                                 rbcRegion%rboundaryRegion, &
                                 IminEdge,ImaxEdge,icount2)
                                 
  ! Cancel if the set is empty!
  icountmin = MIN(icount,icount2)
  icountmax = MAX(icount,icount2)
  
  IF (icount .EQ. 0) THEN
    RETURN
  END IF
                                 
  ! The min/max values of these give the indices of the elements on the
  ! boundary that share our BC region - because it may be that an edge belongs
  ! to the BC region while the endpoints do not and vice versa, so we
  ! have to make sure to get an index which covers everything.
  
  IF (icount .EQ. icount2) THEN
    Iminidx = MIN(IminVertex,IminEdge)
    Imaxidx = MAX(ImaxVertex,ImaxEdge)
  ELSE
    IF (icount .GT. icount2) THEN
      ! More vertex sets than edge sets. Transfer the vetex set and take the
      ! min only with the first edge set
      Iminidx = IminVertex
      Imaxidx = ImaxVertex
    ELSE
      ! More edge sets than vertex sets. Transfer the vetex set and take the
      ! min only with the first vertex set
      Iminidx = IminVertex
      Imaxidx = ImaxVertex
    END IF
    ! Note: icountmin is usually =1 in this situation
    Iminidx(1:icountmin) = MIN(IminVertex(1:icountmin),IminEdge(1:icountmin))
    Imaxidx(1:icountmin) = MAX(ImaxVertex(1:icountmin),ImaxEdge(1:icountmin))
  END IF
  
  ! Reserve some memory to save temporarily all DOF's of all boundary
  ! elements.
  ! We handle all boundary elements simultaneously - let's hope that there are 
  ! never so many elements on the boundary that our memory runs out :-)
  ALLOCATE (Idofs(EL_MAXNBAS,ImaxIdx(1)-IminIdx(1)+1+ImaxIdx(2)-IminIdx(2)+1))
  
  ! If the boundary conditions are only to be discretised for modifying the
  ! defect vector, we can skip calculating the values in the boundary.
  IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
    ALLOCATE (DdofValue(EL_MAXNBAS,ImaxIdx(1)-IminIdx(1)+1+ImaxIdx(2)-IminIdx(2)+1))
  END IF
  Idofs = 0

  ! Now the elements with indices iminidx..imaxidx in the ItrialElements
  ! of the triangulation are on the boundary. Some elements may appear
  ! twice (on edges e.g.) but we don't care.
  !
  ! Ask the DOF-mapping routine to get us those DOF's belonging to elements
  ! on the boundary.
  IF (icountmax .GE. 1) THEN
    CALL dof_locGlobMapping_mult(p_rspatialDiscretisation, &
              p_IelementsAtBoundary(Iminidx(1):Imaxidx(1)), .FALSE., Idofs)
  END IF
            
  IF (icountmax .GE. 2) THEN
    CALL dof_locGlobMapping_mult(p_rspatialDiscretisation, &
              p_IelementsAtBoundary(Iminidx(2):Imaxidx(2)), .FALSE., &
              Idofs( :, Imaxidx(1)-Iminidx(1)+1+1 : ))
  END IF
                               
  ! Loop through the index sets
  DO ipart = 1,icountmax
  
    ! Now the element-dependent part. For each element type, we have to
    ! figure out which DOF's are on the boundary!
    DO i=Iminidx(ipart),Imaxidx(ipart)
    
      ! We proceed as follows: We figure out, which DOF is on the
      ! boundary. Then, we ask our computation routine to calculate
      ! the necessary value and translate them into a DOF value.
      ! All DOF values are collected later.

      ! Where are we at the boundary? Element? Edge? Adjacent vertices?  
      ielement = p_IelementsAtBoundary(I)
      iedge = p_IedgesAtBoundary(I)
      ipoint1 = p_IverticesAtEdge(1,iedge-NVT)
      ipoint2 = p_IverticesAtEdge(2,iedge-NVT)
      
      ! Maybe the points are in the wrong order...
      IF (ipoint1 .NE. p_IverticesAtBoundary(I)) THEN
        ipoint1 = p_IverticesAtEdge(2,I)
        ipoint2 = p_IverticesAtEdge(1,I)
        ! now they are not!
      END IF

      ! Which DOF's are on the current boundary segment? Find out where in 
      ! the element the edge on the boundary is oriented.
      DO ilocalEdge = 1,TRIA_MAXNME2D
        IF (p_IedgesAtElement(ilocalEdge,ielement) .EQ. iedge) EXIT
      END DO
      
      IF (ilocalEdge .GT. TRIA_MAXNME2D) THEN
        PRINT *,'Error in bcasm_discrBCDirichlet: Edge not found!'
        STOP
      END IF
      
      ! Now, ilocalEdge is the 'local' number of the edge
      ! corresponding to the 'global' edge number iedge
      !
      ! Get the element type in case we don't have a uniform triangulation.
      ! Otherwise, ieltype was set to the trial element type above.
      IF (p_rspatialDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) THEN
        ieltype = p_ItrialElements(p_IelementsAtBoundary(I))
      END IF
      
      SELECT CASE (ieltype)
      
      CASE (EL_P0,EL_Q0)
        ! Nice element, only one DOF :-)
        ! Either the edge or an adjacent vertex is on the boundary.
        
        IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
          ! Get the value at the corner point and accept it as
          ! Dirichlet value.
          CALL fgetBoundaryValues (Icomponents,p_rspatialDiscretisation,&
                                  rbcRegion,ielement, DISCBC_NEEDFUNC,&
                                  ipoint1,p_DvertexParameterValue(I), &
                                  p_rcollection, Dvalues)
                                 
          ! Dvalues(1) gives us the function value in the point. Save it
          DdofValue(1,I-Iminidx(ipart)+1) = Dvalues(1)
          
        END IF
        
        ! Set the DOF number < 0 to indicate that it is Dirichlet
        Idofs(1,I-Iminidx(ipart)+1) = -ABS(Idofs(1,I-Iminidx(ipart)+1))
        
      CASE (EL_P1,EL_Q1)

        ! Left point inside? -> Corresponding DOF must be computed
        IF ( (I .GE. IminVertex(ipart)) .AND. (I .LE. ImaxVertex(ipart)) ) THEN
        
          IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN          
            CALL fgetBoundaryValues (Icomponents,p_rspatialDiscretisation,&
                                    rbcRegion,ielement, DISCBC_NEEDFUNC,&
                                    ipoint1,p_DvertexParameterValue(I), &
                                    p_rcollection, Dvalues)
                                    
            ! Save the computed function value
            DdofValue(ilocalEdge,I-Iminidx(ipart)+1) = Dvalues(1) 
          END IF
          
          ! Set the DOF number < 0 to indicate that it is Dirichlet.
          ! ilocalEdge is the number of the local edge - and at the same
          ! time the number of the local DOF of Q1, as an edge always
          ! follows a corner vertex!
          Idofs(ilocalEdge,I-Iminidx(ipart)+1) = -ABS(Idofs(ilocalEdge,I-Iminidx(ipart)+1))
        END IF
        
        ! The right point does not have to be checked! It comes later
        ! with the next edge. The situation when an element crosses the
        ! maximum parameter value with its boundary is handled by the
        ! outer DO-LOOP:
        ! A boundary region with parameter value e.g. [3.0,TMAX]
        ! will produce two index sets: One index set for [0.0, 0.0]
        ! and one for [3.0, TMAX).
        
      CASE (EL_EM30,EL_E030)

        ! Edge inside? -> Calculate integral mean value over the edge
        IF ( (I .GE. IminEdge(ipart)) .AND. (I .LE. ImaxEdge(ipart)) ) THEN
          IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
            CALL fgetBoundaryValues (Icomponents,p_rspatialDiscretisation,&
                                    rbcRegion,ielement, DISCBC_NEEDINTMEAN,&
                                    iedge,p_DedgeParameterValue(I), &
                                    p_rcollection, Dvalues)
                                    
            ! Save the computed function value
            DdofValue(ilocalEdge,I-Iminidx(ipart)+1) = Dvalues(1) 
          END IF
          
          ! Set the DOF number < 0 to indicate that it is Dirichlet
          Idofs(ilocalEdge,I-Iminidx(ipart)+1) = -ABS(Idofs(ilocalEdge,I-Iminidx(ipart)+1))
        END IF
        
      CASE (EL_EM31,EL_E031)

        ! Edge inside? -> Calculate point value on midpoint of edge iedge
        IF ( (I .GE. IminEdge(ipart)) .AND. (I .LE. ImaxEdge(ipart)) ) THEN
          IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
            CALL fgetBoundaryValues (Icomponents,p_rspatialDiscretisation,&
                                    rbcRegion,ielement, DISCBC_NEEDFUNCMID,&
                                    iedge,p_DedgeParameterValue(I), &
                                    p_rcollection, Dvalues)
                                    
            ! Save the computed function value
            DdofValue(ilocalEdge,I-Iminidx(ipart)+1) = Dvalues(1) 
          END IF
          
          ! Set the DOF number < 0 to indicate that it is Dirichlet
          Idofs(ilocalEdge,I-Iminidx(ipart)+1) = -ABS(Idofs(ilocalEdge,I-Iminidx(ipart)+1))
        END IF
        
      CASE DEFAULT
      
        PRINT *,'bcasm_discrBCDirichlet: Unsupported element!'
        STOP
      
      END SELECT
    
    END DO
  
  END DO
  
  ! Now count how many values we actually have.
  icount = 0
  DO J=1,SIZE(Idofs,2)
    DO I=1,SIZE(Idofs,1)
      IF (Idofs(I,J) < 0) icount = icount + 1
    END DO
  END DO
  
  p_rdirichletBCs%nDOF = icount
  
  IF (icount .GT. 0) THEN
  
    ! Allocate arrays for storing these DOF's and their values - if values are
    ! computed.
    CALL storage_new('bcasm_discrBCDirichlet', 'h_IdirichletDOFs', &
                    icount, ST_INT, p_rdirichletBCs%h_IdirichletDOFs, &
                    ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int(p_rdirichletBCs%h_IdirichletDOFs,p_IdirichletDOFs)
    
    IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
      CALL storage_new('bcasm_discrBCDirichlet', 'h_DdirichletValues', & 
                      icount, ST_DOUBLE, p_rdirichletBCs%h_DdirichletValues, &
                      ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double(p_rdirichletBCs%h_DdirichletValues,p_IdirichletValues)
    END IF
    
    ! Transfer the DOF's and their values to these arrays.
    icount = 0
    DO J=1,SIZE(Idofs,2)
      DO I=1,SIZE(Idofs,1)
        IF (Idofs(I,J) < 0) THEN
          icount = icount + 1
          p_IdirichletDOFs(icount) = ABS(Idofs(I,J))
          IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
            p_IdirichletValues(icount) = DdofValue(I,J)
          END IF
        END IF
      END DO
    END DO
  
  ELSE
  
    ! Let's hope there is nothing saved here :-)
    p_rdirichletBCs%h_IdirichletDOFs = ST_NOHANDLE
    p_rdirichletBCs%h_DdirichletValues = ST_NOHANDLE
  
  END IF
  
  ! Remove temporary memory, finish.
  IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
    DEALLOCATE (DdofValue)
  END IF
  DEALLOCATE (Idofs)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcasm_releaseDirichlet (rdiscreteBCDirichlet)
  
!<description>
  ! This routine cleans up the discrete Dirichlet boundary conditions
  ! rdiscreteBCDirichlet.
!</description>

!<inputoutput>
  ! The discrete-BC structure which is to be cleaned up
  TYPE(t_discreteBCDirichlet), INTENT(INOUT) :: rdiscreteBCDirichlet
!</inputoutput>

!</subroutine>

  ! Release what is associated
  
  rdiscreteBCDirichlet%nDOF = 0
  IF (rdiscreteBCDirichlet%h_DdirichletValues .NE. ST_NOHANDLE) &
    CALL storage_free(rdiscreteBCDirichlet%h_DdirichletValues)
  IF (rdiscreteBCDirichlet%h_IdirichletDOFs .NE. ST_NOHANDLE) &
    CALL storage_free(rdiscreteBCDirichlet%h_IdirichletDOFs)

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcasm_discrBCpressureDrop (rblockDiscretisation, &
                                     rbcRegion, rdiscreteBC, casmComplexity, &
                                     fgetBoundaryValues,p_rcollection)
  
!<description>
  ! Creates a discrete version of pressure drop boundary conditions.
  ! rbcRegion describes the region which is to be discretised. The discretised
  ! boundary conditions are created in rdiscreteBC, which is assumed
  ! to be undefined when entering this routine.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation. The boundary
  ! conditions inside of this structure are discretised.
  TYPE(t_blockDiscretisation), INTENT(IN), TARGET :: rblockDiscretisation

  ! The BC region which is to be discrised
  TYPE(t_bcRegion), INTENT(IN) :: rbcRegion

  ! A callback function that calculates values on the boundary.
  ! Is declared in the interface include file 'intf_bcassembly.inc'.
  INCLUDE 'intf_bcassembly.inc'
  
  ! A collection structure to inform the callback function with
  ! additional information. Can be NULL() if there is no information to pass.
  TYPE(t_collection), POINTER :: p_rcollection

  ! A combination of BCASM_DISCFORxxx constants that specify
  ! the complexity of the discretisation that is to perform. This allows to
  ! discretise only parts of the BC's, e.g. only setting up those
  ! information that are necessary for filtering defect vectors.
  ! If not specified, BCASM_DISCFORALL is assumed, i.e. the resulting
  ! boundary conditions can be used for everything.
  INTEGER(I32), INTENT(IN) :: casmComplexity
!</input>  

!<output>
  ! This structure receives the result of the discretisation of rbcRegion.
  ! When entering the routine, the content of this structure is undefined,
  ! all pointers are invalid. The routine fills everything with appropriate
  ! data.
  TYPE(t_discreteBCEntry), INTENT(OUT), TARGET :: rdiscreteBC
!</output>

!</subroutine>

  ! local variables
  INTEGER :: i,icount,ipart,ieltype
  INTEGER, DIMENSION(2) :: IminEdge,ImaxEdge,Iminidx,Imaxidx
  REAL(DP), DIMENSION(DER_MAXNDER)            :: Dvalues
  REAL(DP),DIMENSION(NDIM2D)                  :: Dtangential,Dnormal
  INTEGER(PREC_POINTIDX)                      :: NVT,ipoint1,ipoint2
  INTEGER(PREC_ELEMENTIDX)                    :: ielement
  INTEGER(PREC_EDGEIDX)                       :: ndofs,idof,iedge
  INTEGER(I32), DIMENSION(2)                  :: ImodifierSize
  
  TYPE(t_spatialDiscretisation), POINTER      :: p_rspatialDiscretisation
  TYPE(t_triangulation), POINTER              :: p_rtriangulation
  INTEGER(I32), DIMENSION(:), POINTER         :: p_IelementsAtBoundary
  INTEGER(I32), DIMENSION(:), POINTER         :: p_IedgesAtBoundary,p_IverticesAtBoundary
  REAL(DP), DIMENSION(:), POINTER             :: p_DedgeParameterValue
  REAL(DP), DIMENSION(:,:), POINTER           :: p_DcornerCoordinates
  INTEGER(I32), DIMENSION(:,:), POINTER       :: p_IverticesAtEdge
  
  TYPE(t_discreteBCpressureDrop), POINTER     :: p_rpressureDropBCs
  INTEGER(PREC_DOFIDX), DIMENSION(:), POINTER :: p_IpressureDropDOFs
  REAL(DP), DIMENSION(:,:), POINTER           :: p_Dmodifier

  ! Pressure drop BC's only exist as modification of the RHS.
  ! If we should not compute them for the RHS, we don't have to do anything.
  IF (IAND(casmComplexity,BCASM_DISCFORRHS) .EQ. 0) RETURN

  ! Fill the structure for discrete pressure drop BC's in the
  ! t_discreteBCEntry structure
  p_rpressureDropBCs => rdiscreteBC%rpressureDropBCs
  
  ! Get the discretisation structures from one of the components of the solution
  ! vector that is to be modified.
  p_rspatialDiscretisation => &
    rblockDiscretisation%RspatialDiscretisation(rbcRegion%Iequations(1))  

  IF (p_rspatialDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) THEN
    PRINT *,'Discrete pressure drop boundary conditions currently only supported'
    PRINT *,'for uniform discretisations!'
    STOP
  END IF

  ieltype = p_rspatialDiscretisation%RelementDistribution(1)%itrialElement
  IF ((ieltype .NE. EL_E030) .AND. (ieltype .NE. EL_E031) .AND. &
      (ieltype .NE. EL_EM30) .AND. (ieltype .NE. EL_E031)) THEN
    PRINT *,'Discrete pressure drop boundary conditions currently only supported'
    PRINT *,'for Q1~ element!'
    STOP
  END IF
  
  IF (rbcRegion%nequations .NE. NDIM2D) THEN
    PRINT *,'Pressure drop boundary conditions only support 2D!'
    STOP
  END IF
  
  ! For easier access:
  p_rtriangulation => p_rspatialDiscretisation%p_rtriangulation

  ! Note: All elements are of the same type ieltyp.
  !
  ! We have pressure drop boundary conditions
  rdiscreteBC%itype = DISCBC_TPPRESSUREDROP
  
  ! Connect to the boundary condition structure
  rdiscreteBC%p_rboundaryConditions => rblockDiscretisation%p_rboundaryConditions
  
  ! Which components of the solution vector are affected by this boundary
  ! condition?
  p_rpressureDropBCs%ncomponents = rbcRegion%nequations
  p_rpressureDropBCs%Icomponents(1:NDIM2D) = rbcRegion%Iequations(1:NDIM2D)
  
  ! We have to deal with all DOF's on the boundary. This is highly element
  ! dependent and therefore a little bit tricky :(
  ! But here we restrict to Q1~ only, which makes life a little bit easier.
  !
  ! As we are in 2D, we can use parameter values at first to figure out,
  ! which edges are on the boundary.
  ! What we have is a boundary segment. Now ask the boundary-index routine
  ! to give us the minimum and maximum index of the edges on the
  ! bondary that belong to this boundary segment.
  
  CALL bcasm_getEdgesInBCregion (p_rtriangulation,p_rspatialDiscretisation%p_rdomain, &
                                 rbcRegion%rboundaryRegion, &
                                 IminEdge,ImaxEdge,icount)
                                 
  ! Cancel if the set is empty!
  IF (icount .EQ. 0) THEN
    RETURN
  END IF
                    
  ! Put IminEdge/ImaxEdge to Iminidx/ImaxIdx and continue working with these.
  ! in a later implementation, we probably have to include the indices of
  ! points on the boundary here, too, like in the Dirichlet case.
  Iminidx = IminEdge
  Imaxidx = ImaxEdge
  
  ! Total number of edges?
  ndofs = Imaxidx(1)-Iminidx(1)+1 + Imaxidx(2)-Iminidx(2)+1
  
  p_rpressureDropBCs%nDOF = ndofs

  ! Allocate memory to save the DOF's as well as all modifiers.
  CALL storage_new('bcasm_discrBCpressureDrop', 'h_IpressureDropDOFs', &
                  ndofs, ST_INT, p_rpressureDropBCs%h_IpressureDropDOFs, &
                  ST_NEWBLOCK_NOINIT)
  ImodifierSize = (/NDIM2D,ndofs/)
  CALL storage_new2D('bcasm_discrBCDirichlet', 'h_Dmodifier', & 
                    ImodifierSize, ST_DOUBLE, p_rpressureDropBCs%h_Dmodifier, &
                    ST_NEWBLOCK_NOINIT)
                    
  CALL storage_getbase_int(p_rpressureDropBCs%h_IpressureDropDOFs,p_IpressureDropDOFs)
  CALL storage_getbase_double2d(p_rpressureDropBCs%h_Dmodifier,p_Dmodifier)

  ! For easier access:
  CALL storage_getbase_double(p_rtriangulation%h_DedgeParameterValue,p_DedgeParameterValue)
  CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
  CALL storage_getbase_int(p_rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)
  CALL storage_getbase_int(p_rtriangulation%h_IelementsAtBoundary,p_IelementsAtBoundary)
  CALL storage_getbase_double2D(p_rtriangulation%h_DcornerCoordinates,p_DcornerCoordinates)
  CALL storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
  NVT = p_rtriangulation%NVT

  ! Now calculate the pressure drop integral; cf. p. 257 (235) in Turek's book:
  !
  ! The pressure drop boundary condition has to implement
  !
  !     - sum P_j  int_Sj  phi * n  ds
  !
  ! into the RHS vector. For each (velocity) DOF of the boundary,
  ! we save "P_j  int_Sj  phi_k * n  ds" as modifier for the DOF of
  ! the RHS!
  
  idof = 0
  
  ! Loop through the index sets. Normally we have only one, except when
  ! a boundary segment crosses the maximum parameter value...
  DO ipart = 1,icount
  
    ! Loop through all edges on the boundary belonging to our current 
    ! boundary segment.
    DO i=Iminidx(ipart),Imaxidx(ipart)
    
      ! Where are we at the boundary? Element? Edge? Adjacent vertices?  
      ielement = p_IelementsAtBoundary(I)
      iedge = p_IedgesAtBoundary(I)
      ipoint1 = p_IverticesAtEdge(1,iedge-NVT)
      ipoint2 = p_IverticesAtEdge(2,iedge-NVT)
      
      ! Maybe the points are in the wrong order...
      IF (ipoint1 .NE. p_IverticesAtBoundary(I)) THEN
        ipoint1 = p_IverticesAtEdge(2,I)
        ipoint2 = p_IverticesAtEdge(1,I)
        ! now they are not!
      END IF
      
      ! Get the coordinates of the endpoints to build the tangential
      ! vector of the edge:
      Dtangential(1:NDIM2D) = p_DcornerCoordinates(1:NDIM2D,ipoint2) &
                            - p_DcornerCoordinates(1:NDIM2D,ipoint1)
                            
      ! Get the inner normal vector. This compensates the '-' sign in front of
      ! the RHS in the formula on poage 269 in Turek's book where the outer
      ! normal vector is used.
      Dnormal(1) = -Dtangential(2)
      Dnormal(2) =  Dtangential(1)
      
      ! Don't scale the normal vector! The scaling factor cancels out later
      ! when calculating the integral on the edge with the midpoint rule!
      !
      ! At first, we ask the boundary-value routine to give us the
      ! weight P_j in the cubature point (for Q1~ in the edge midpoint)
      ! of the current boundary segment j we are discretising here!
      !
      ! Calculate normal stress in the current midpoint of the edge.
      ! In a later implementation, we might calculate the value in a cubature
      ! point to calculate the current integral - but here we use the midpoint
      ! rule...
      CALL fgetBoundaryValues (p_rpressureDropBCs%Icomponents,p_rspatialDiscretisation,&
                                rbcRegion,ielement, DISCBC_NEEDNORMALSTRESS,&
                                iedge,p_DedgeParameterValue(I), &
                                p_rcollection, Dvalues)
      
      ! Save the current modifier (normal stress value)*n to the structure for
      ! later implementation in to the RHS vector.
      idof = idof + 1
      
      p_IpressureDropDOFs(idof) = iedge-NVT
      p_Dmodifier(1:NDIM2D,idof) = Dvalues(1:NDIM2D)*Dnormal(1:NDIM2D)
    
    END DO ! i
  
  END DO ! ipart  
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcasm_releasePressureDrop (rdiscreteBCPD)
  
!<description>
  ! This routine cleans up the discrete Pressure Drop boundary conditions
  ! rdiscreteBCPD.
!</description>

!<inputoutput>
  ! The discrete-BC structure which is to be cleaned up
  TYPE(t_discreteBCpressureDrop), INTENT(INOUT) :: rdiscreteBCPD
!</inputoutput>

!</subroutine>

  ! Release what is associated
  
  rdiscreteBCPD%nDOF = 0
  IF (rdiscreteBCPD%h_IpressureDropDOFs .NE. ST_NOHANDLE) &
    CALL storage_free(rdiscreteBCPD%h_IpressureDropDOFs)
  IF (rdiscreteBCPD%h_Dmodifier .NE. ST_NOHANDLE) &
    CALL storage_free(rdiscreteBCPD%h_Dmodifier)

  END SUBROUTINE

! *****************************************************************************
! Support for boundary conditions on fictitious boundary components
! *****************************************************************************

!<subroutine>

  SUBROUTINE bcasm_discretiseFBC (rblockDiscretisation,p_rdiscreteFBC,bforceRebuild, &
                                  fgetBoundaryValuesFBC,rcollection,casmComplexity)
  
!<description>
  ! This routine discretises an analytic definition of boundary conditions
  ! on fictitious boundary components.
  ! The definition of the boundary conditions is taken from the discretisation
  ! structure rspatialDiscretisation. The discrete version is build up in
  ! p_rdiscreteFBC. If p_rdiscreteFBC is NULL(), a new structure is created,
  ! otherwise the old structure is updated (or even destroyed and recreated if
  ! necessary).
!</description>

!<input>
  
  ! The block discretisation structure of the underlying PDE. The boundary
  ! conditions inside of this structure are discretised.
  TYPE(t_blockDiscretisation), INTENT(IN) :: rblockDiscretisation
  
  ! Can be set to TRUE to force a complete rebuild of the rdiscreteFBC structure.
  ! Normally, the structure is completely set up only in the first call
  ! or is there is a massive change in the boundary conditions (e.g. the
  ! number change)
  ! In later calls, the structure is recomputed only in those parts of the boundary,
  ! which have the t_bcRegion\%bisstatic flag set FALSE (the standard value) -
  ! except there , in which case the whole structure is rebuild.
  ! By setting bforceRebuild to TRUE, one can enforce a complete
  ! rebuild of the structure, independent of which regions are marked
  ! as static.
  LOGICAL                                   :: bforceRebuild
  
  ! A callback function that calculates values on the boundary.
  ! Is declared in the interface include file 'intf_bcassembly.inc'.
  INCLUDE 'intf_fbcassembly.inc'
  
  ! OPTIONAL: A collection structure to inform the callback function with
  ! additional information. Can undefined if there is no
  ! information to pass.
  TYPE(t_collection), TARGET, OPTIONAL :: rcollection
  
  ! OPTIONAL: A combination of BCASM_DISCFORxxx constants that specify
  ! the complexity of the discretisation that is to perform. This allows to
  ! discretise only parts of the BC's, e.g. only setting up those
  ! information that are necessary for filtering defect vectors.
  ! If not specified, BCASM_DISCFORALL is assumed, i.e. the resulting
  ! boundary conditions can be used for everything.
  INTEGER(I32), INTENT(IN), OPTIONAL :: casmComplexity
!</input>

!<inputoutput>
  ! A discretised version of the analytic boundary conditions.
  ! This is a pointer to a t_discreteBC structures, 
  ! representing the boundary discretised in a discretisation- dependent way.
  ! If this pointer points to NULL(), a complete new structure is set up.
  TYPE(t_discreteFBC), POINTER :: p_rdiscreteFBC
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: icurrentRegion, ccompl
  TYPE(t_bcRegion), POINTER :: p_rbcRegion
  LOGICAL :: bbuildAll
  TYPE(t_collection), POINTER :: p_rcoll
  
  ! Pointer to the boundary condition object
  TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions
  
  ! For quicker access:
  p_rboundaryConditions => rblockDiscretisation%p_rboundaryConditions
  
  ! We replace the optional parameter by NULL() if it does not exist
  IF (PRESENT(rcollection)) THEN
    p_rcoll => rcollection
  ELSE
    p_rcoll => NULL()
  END IF
  
  IF (.NOT. ASSOCIATED(p_rboundaryConditions)) THEN
    PRINT *,'Warning in bcasm_discretiseBFC: Boundary conditions not associated!'
    RETURN
  END IF
  
  ! Default target for the discretisation if everything.
  IF (PRESENT(casmComplexity)) THEN
    ccompl = casmComplexity
  ELSE
    ccompl = BCASM_DISCFORALL
  END IF
  
  bbuildAll = bforceRebuild
  
  ! Is there a structure we can work with?
  IF (.NOT. ASSOCIATED(p_rdiscreteFBC)) THEN
    ! Create a new structure on the heap.
    ALLOCATE(p_rdiscreteFBC)
  END IF
  
  ! Now we work with the array in the structure. Does it exist?
  
  IF (ASSOCIATED(p_rdiscreteFBC%p_RdiscFBCList)) THEN
  
    ! Is there a massive change or do we have o rebuild everything?
    IF ((p_rboundaryConditions%iregionCountFBC .NE. &
         SIZE(p_rdiscreteFBC%p_RdiscFBCList)) &
        .OR. bforceRebuild) THEN
        
      ! Oh, we have to destroy the structure completely.
      ! Don't remove the structure itself from the heap - we want to
      ! fill it with data immediately afterwards!
      CALL bcasm_releaseDiscreteFBC (p_rdiscreteFBC,.TRUE.)

      ! Allocate a new structure array with iregionCount 
      ! entries for all the boundary conditions.
      ALLOCATE(p_rdiscreteFBC%p_RdiscFBCList(p_rboundaryConditions%iregionCountFBC))
      
      bbuildAll = .TRUE.
    ELSE
      ! Otherwise, release only those information belonging to 
      ! non-static boundary regions
      DO icurrentRegion = 1,p_rboundaryConditions%iregionCountFBC
        
        ! Get a pointer to it so we can deal with it more easily
        p_rbcRegion => p_rboundaryConditions%p_RregionsFBC(icurrentRegion)
      
        IF (.NOT. p_rbcRegion%bisStatic) THEN
        
          ! Release all allocated information to this boundary region
          SELECT CASE (p_rdiscreteFBC%p_RdiscFBCList(icurrentRegion)%itype)

          CASE (DISCBC_TPDIRICHLET)

            ! Discrete Dirichlet boundary conditions. Release the old structure.
            CALL bcasm_releaseFBCDirichlet( &
                 p_rdiscreteFBC%p_RdiscFBCList(icurrentRegion)%rdirichletFBCs)

          END SELECT

          ! BC released, indicate this
          p_rdiscreteFBC%p_RdiscFBCList(icurrentRegion)%itype = DISCBC_TPUNDEFINED
          
        END IF
      
      END DO  
    END IF
    
  ELSE
  
    ! Allocate a structure array with iregionCount entries
    ! for all the boundary conditions.
    ALLOCATE(p_rdiscreteFBC%p_RdiscFBCList(p_rboundaryConditions%iregionCountFBC))
    bbuildAll = .TRUE.
  
  END IF
  
  ! Loop through the regions on the boundary
  DO icurrentRegion = 1,p_rboundaryConditions%iregionCountFBC
    
    ! Get a pointer to it so we can deal with it more easily
    p_rbcRegion => p_rboundaryConditions%p_RregionsFBC(icurrentRegion)
  
    ! Do we have to process this region?
    IF (bforceRebuild .OR. .NOT. p_rbcRegion%bisStatic) THEN
    
      ! Ok, let's go...
      ! What for BC do we have here?
      SELECT CASE (p_rbcRegion%ctype)
      
      CASE (BC_DIRICHLET)
        CALL bcasm_discrFBCDirichlet (rblockDiscretisation, &
                  p_rbcRegion, p_rdiscreteFBC%p_RdiscFBCList(icurrentRegion), &
                  ccompl,fgetBoundaryValuesFBC,p_rcoll)
                    
      END SELECT
        
    END IF
  
  END DO  

  END SUBROUTINE
      
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcasm_releaseDiscreteFBC (p_rdiscreteFBC,bkeepBCStructure)
  
!<description>
  ! This routine cleans up an array of t_discreteFBCEntry describing discrete
  ! boundary conditions for fictitious boundary components. 
  ! All allocated memory is released. The structure
  ! p_rdiscreteFBC is released from the heap as well.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the structure p_rdiscreteFBC is not released
  ! from memory. If set to FALSE or not existent (the usual setting), the 
  ! structure p_rdiscreteFBC will also be removed from the heap after 
  ! cleaning up.
  LOGICAL, INTENT(IN), OPTIONAL :: bkeepBCStructure
!</input>

!<inputoutput>
  ! A pointer to discretised boundary conditions for fictitious boundary
  ! components. All memory allocated in and by this array is 
  ! released from the heap. 
  TYPE(t_discreteFBC), POINTER :: p_rdiscreteFBC
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: icurrentRegion
  
  IF (.NOT. ASSOCIATED(p_rdiscreteFBC)) THEN
    PRINT *,'Warning in bcasm_releaseDiscreteFBC: Nothing to cleanup!'
    RETURN
  END IF
  
  ! Destroy the content of the structure completely!
  IF (ASSOCIATED(p_rdiscreteFBC%p_RdiscFBCList)) THEN
  
    ! Destroy all the substructures in the array.
    DO icurrentRegion = 1,SIZE(p_rdiscreteFBC%p_RdiscFBCList)  
      
      ! Release all allocated information to this boundary region
      SELECT CASE (p_rdiscreteFBC%p_RdiscFBCList(icurrentRegion)%itype)
      CASE (DISCBC_TPDIRICHLET)
        ! Discrete Dirichlet boundary conditions. Release the old structure.
        CALL bcasm_releaseFBCDirichlet(&
          p_rdiscreteFBC%p_RdiscFBCList(icurrentRegion)%rdirichletFBCs)

      END SELECT

      ! BC released, indicate this
      p_rdiscreteFBC%p_RdiscFBCList(icurrentRegion)%itype = DISCBC_TPUNDEFINED
      
    END DO  
    
    ! Release the array itself.
    DEALLOCATE(p_rdiscreteFBC%p_RdiscFBCList)
    
  END IF
  
  ! Deallocate the structure (if we are allowed to), finish.
  IF (.NOT. PRESENT(bkeepBCStructure)) THEN
    DEALLOCATE(p_rdiscreteFBC)
  ELSE
    IF (.NOT. bkeepBCStructure) THEN
      DEALLOCATE(p_rdiscreteFBC)
    END IF
  END IF

  END SUBROUTINE
      
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcasm_discrFBCDirichlet (rblockDiscretisation, &
                                      rbcRegion, rdiscreteFBC, casmComplexity, &
                                      fgetBoundaryValuesFBC,p_rcollection)
  
!<description>
  ! Creates a discrete version of Dirichlet boundary conditions for a fictitious
  ! boundary component.
  ! rbcRegion describes the region which is to be discretised. The discretised
  ! boundary conditions are created in rdiscreteFBC, which is assumed
  ! to be undefined when entering this routine.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation. The boundary
  ! conditions inside of this structure are discretised.
  TYPE(t_blockDiscretisation), INTENT(IN), TARGET :: rblockDiscretisation

  ! The boundary condition region which is to be discrised. Must be configured
  ! as boundary condition region of a fictitious boundary object.
  TYPE(t_bcRegion), INTENT(IN) :: rbcRegion

  ! A callback function that calculates values in the domain.
  ! Is declared in the interface include file 'intf_fbcassembly.inc'.
  INCLUDE 'intf_fbcassembly.inc'
  
  ! A collection structure to inform the callback function with
  ! additional information. Can be NULL() if there is no information to pass.
  TYPE(t_collection), POINTER :: p_rcollection

  ! A combination of BCASM_DISCFORxxx constants that specify
  ! the complexity of the discretisation that is to perform. This allows to
  ! discretise only parts of the BC's, e.g. only setting up those
  ! information that are necessary for filtering defect vectors.
  ! If not specified, BCASM_DISCFORALL is assumed, i.e. the resulting
  ! boundary conditions can be used for everything.
  INTEGER(I32), INTENT(IN) :: casmComplexity
!</input>  

!<output>
  ! This structure receives the result of the discretisation of rbcRegion.
  ! When entering the routine, the content of this structure is undefined,
  ! all pointers are invalid. The routine fills everything with appropriate
  ! data.
  TYPE(t_discreteFBCEntry), INTENT(OUT), TARGET :: rdiscreteFBC
!</output>

!</subroutine>

    ! local variables
    TYPE(t_triangulation), POINTER              :: p_rtriangulation
    TYPE(t_spatialDiscretisation), POINTER      :: p_rspatialDiscretisation
    
    INTEGER(PREC_DOFIDX) :: nDOFs
    INTEGER :: nequations, h_Ddofs, h_Idofs, i, j, icomponent, ieltype
    INTEGER(PREC_DOFIDX), DIMENSION(2) :: IdofCount
    
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idofs
    REAL(DP), DIMENSION(:,:), POINTER   :: p_Ddofs

    TYPE(t_discreteFBCDirichlet),POINTER        :: p_rdirichletFBCs
    
    INTEGER(I32), DIMENSION(FBCASM_MAXSIM), TARGET      :: Isubset
    INTEGER, DIMENSION(FBCASM_MAXSIM), TARGET           :: Iinside
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET     :: p_Dsubset
    INTEGER(I32) :: isubsetStart, isubsetLength, icurrentDof
    
    TYPE(t_discreteFBCevaluation), DIMENSION(DISCFBC_MAXDISCBC) :: Revaluation
    
    ! Get the discretisation structure of the first component.
    ! In this first rough implementation, we only support uniform a uniform 
    ! discretisation, so check that we are able to handle the current situation.
    icomponent = rbcRegion%Iequations(1)
    p_rspatialDiscretisation => rblockDiscretisation%RspatialDiscretisation(icomponent)
    
    IF (p_rspatialDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) THEN
      PRINT *,'bcasm_discrFBCDirichlet: Can only handle uniform discretisation!'
      STOP
    END IF

    IF (p_rspatialDiscretisation%ndimension .NE. NDIM2D) THEN
      PRINT *,'bcasm_discrFBCDirichlet: Only 2D supported for now!'
      STOP
    END IF

    ! For easier access:
    p_rtriangulation => p_rspatialDiscretisation%p_rtriangulation

    ! All elements are of the samne type. Get it in advance.
    ieltype = p_rspatialDiscretisation%RelementDistribution(1)%itrialElement
    
    ! We have Dirichlet boundary conditions
    rdiscreteFBC%itype = DISCFBC_TPDIRICHLET
    
    ! Connect to the boundary condition structure
    rdiscreteFBC%p_rboundaryConditions => rblockDiscretisation%p_rboundaryConditions
    
    ! Fill the structure for discrete Dirichlet BC's in the
    ! t_discreteBCEntry structure
    p_rdirichletFBCs => rdiscreteFBC%rdirichletFBCs
    
    nequations = rbcRegion%nequations
    p_rdirichletFBCs%ncomponents = nequations
    p_rdirichletFBCs%Icomponents(1:nequations) = rbcRegion%Iequations(1:nequations)
    
    ! Allocate memory for intermediate values
    ALLOCATE(p_Dsubset(FBCASM_MAXSIM,1,nequations))
    
    ! We have to collect all DOF's and their values that belong to our current
    ! fictitious boundary object. Depending on the element type we will loop
    ! through all vertices, edges and elements in the triangulation
    ! to collect all the important values.
    !
    ! Allocate an array as large as a solution vector would be to store the DOF
    ! values. Allocate an integer array of the same size that receives a flag
    ! whether the DOF is actually used by our current FB object.
    ! More specifically, we remember the DOF value for each of the components!
    
    nDOFs = dof_igetNDofGlob(p_rspatialDiscretisation)
    
    IdofCount = (/nequations,nDOFs/)
    IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
      CALL storage_new2d ('bcasm_discrFBCDirichlet', 'DofValues', &
                        IdofCount, ST_DOUBLE, h_Ddofs, &
                        ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double2d (h_Ddofs,p_Ddofs)
    END IF
    CALL storage_new ('bcasm_discrFBCDirichlet', 'DofUsed', nDOFs, &
                      ST_INT, h_Idofs, ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int (h_Idofs,p_Idofs)
    
    ! Depending on the element type, prepare the evaluation structure and call
    ! the callback routine to calculate what we need.
    
    icurrentDof = 0
    
    ! Calculate values in the vertices for Q1,Q2,P1,P2
    IF ((ieltype .EQ. EL_P1) .OR. &
        (ieltype .EQ. EL_Q1) .OR. &
        (ieltype .EQ. EL_P2) .OR. &
        (ieltype .EQ. EL_Q2)) THEN
    
      ! Let's start to collect values. This is a rather element-dependent
      ! part. At first, loop through the vertices in case we have a
      ! P1/Q1/Q2 discretisation
      DO isubsetStart = 1,p_rtriangulation%NVT,FBCASM_MAXSIM
      
        isubsetLength = MIN(p_rtriangulation%NVT-isubsetStart+1,FBCASM_MAXSIM)
      
        ! Fill the subset with isubsetStart, isubsetStart+1,... to identify the
        ! subset we evaluate.
        CALL fillsubset (isubsetStart,isubsetLength,Isubset)
        
        ! Fill the evaluation structure with data for the callback routine
        DO i=1,nequations
          Revaluation(i)%cinfoNeeded = DISCFBC_NEEDFUNC
          Revaluation(i)%nvalues     = isubsetLength
          Revaluation(i)%p_Iwhere    => Isubset
          Revaluation(i)%p_Dvalues   => p_Dsubset(:,:,i)
          Revaluation(i)%p_Iinside   => Iinside
        END DO

        ! Clear the Iinside array
        Iinside = 0
        
        ! Call the callback routine to calculate the values.
        CALL fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
                                    rblockDiscretisation,rbcRegion, &
                                    Revaluation, p_rcollection)
                                    
        ! Transfer the DOF's that are affected

        IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
          
          DO i=1,isubsetLength
            IF (Iinside(i) .NE. 0) THEN
              icurrentDof = icurrentDof + 1
              DO j=1,nequations
                p_Ddofs(j,icurrentDof) = p_Dsubset(i,1,j)
              END DO
              p_Idofs(icurrentDof) = isubsetStart+i-1
            END IF
          END DO
          
        ELSE
        
          DO i=1,isubsetLength
            IF (Iinside(i) .NE. 0) THEN
              icurrentDof = icurrentDof + 1
              p_Idofs(icurrentDof) = isubsetStart+i-1
            END IF
          END DO
          
        END IF
      
      END DO
    
    END IF
    
    ! Calculate values in the edge midpoints / / integral mean values for Q1~
    IF ((ieltype .EQ. EL_E031) .OR. &
        (ieltype .EQ. EL_EM31) .OR. &
        (ieltype .EQ. EL_E030) .OR. &
        (ieltype .EQ. EL_EM30)) THEN
    
      ! Let's start to collect values. This is a rather element-dependent
      ! part. At first, loop through the vertices in case we have a
      ! P1/Q1/Q2 discretisation
      DO isubsetStart = 1,p_rtriangulation%NMT,FBCASM_MAXSIM
      
        isubsetLength = MIN(p_rtriangulation%NMT-isubsetStart+1,FBCASM_MAXSIM)
      
        ! Fill the subset with isubsetStart, isubsetStart+1,... to identify the
        ! subset we evaluate.
        CALL fillsubset (isubsetStart,isubsetLength,Isubset)
        
        ! Fill the evaluation structure with data for the callback routine
        DO i=1,nequations
          IF ((ieltype .EQ. EL_E031) .OR. &
              (ieltype .EQ. EL_EM31)) THEN
            Revaluation(i)%cinfoNeeded = DISCFBC_NEEDFUNCMID
          ELSE
            Revaluation(i)%cinfoNeeded = DISCFBC_NEEDINTMEAN
          END IF
          Revaluation(i)%nvalues     = isubsetLength
          Revaluation(i)%p_Iwhere    => Isubset
          Revaluation(i)%p_Dvalues   => p_Dsubset(:,:,i)
          Revaluation(i)%p_Iinside   => Iinside
        END DO
        
        ! Clear the Iinside array
        Iinside = 0
        
        ! Call the callback routine to calculate the values.
        CALL fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
                                    rblockDiscretisation,rbcRegion, &
                                    Revaluation, p_rcollection)
                                    
        ! Transfer the DOF's that are affected

        IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
          
          DO i=1,isubsetLength
            IF (Iinside(i) .NE. 0) THEN
              icurrentDof = icurrentDof + 1
              DO j=1,nequations
                p_Ddofs(j,icurrentDof) = p_Dsubset(i,1,j)
              END DO
              p_Idofs(icurrentDof) = isubsetStart+i-1
            END IF
          END DO
          
        ELSE
        
          DO i=1,isubsetLength
            IF (Iinside(i) .NE. 0) THEN
              icurrentDof = icurrentDof + 1
              p_Idofs(icurrentDof) = isubsetStart+i-1
            END IF
          END DO
          
        END IF
      
      END DO
    
    END IF
    
    ! Cancel if we did not find any DOF.
    IF (icurrentDof .GT. 0) THEN
      
      ! Reallocate to save memory. Store the final handles in the structure.
      IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
        ! In the 2D-array, the size of the 2nd dimension is changed to the
        ! number of DOF's.
        CALL storage_realloc ('bcasm_discrFBCDirichlet', icurrentDof, &
                              h_Ddofs, ST_NEWBLOCK_NOINIT)
        p_rdirichletFBCs%h_DdirichletValues = h_Ddofs
      END IF
      
      CALL storage_realloc ('bcasm_discrFBCDirichlet', icurrentDof, &
                            h_Idofs, ST_NEWBLOCK_NOINIT)
      p_rdirichletFBCs%h_IdirichletDOFs = h_Idofs
    ELSE
      ! Nothging inside; release arrays
      IF (IAND(casmComplexity,NOT(BCASM_DISCFORDEFMAT)) .NE. 0) THEN
        CALL storage_free (h_Ddofs)
      END IF
      CALL storage_free (h_Idofs)
    END IF

    p_rdirichletFBCs%nDOF = icurrentDof
    
    ! Release temporary data
    DEALLOCATE(p_Dsubset)
    
  CONTAINS
  
    PURE SUBROUTINE fillsubset (istart, ilength, Isubset)
    INTEGER(I32), INTENT(IN) :: istart, ilength
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Isubset
    INTEGER(I32) :: i
      DO i=1,ilength
        Isubset(i) = istart-1+i
      END DO
    END SUBROUTINE
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcasm_releaseFBCDirichlet (rdiscreteFBCDirichlet)
  
!<description>
  ! This routine cleans up the discrete Dirichlet conditions for fictitious 
  ! boundary objects in rdiscreteFBCDirichlet.
!</description>

!<inputoutput>
  ! The discrete-FBC structure which is to be cleaned up
  TYPE(t_discreteFBCDirichlet), INTENT(INOUT) :: rdiscreteFBCDirichlet
!</inputoutput>

!</subroutine>

  ! Release what is associated
  
  rdiscreteFBCDirichlet%nDOF = 0
  rdiscreteFBCDirichlet%Icomponents = 0
  rdiscreteFBCDirichlet%ncomponents = 0
  IF (rdiscreteFBCDirichlet%h_DdirichletValues .NE. ST_NOHANDLE) &
    CALL storage_free(rdiscreteFBCDirichlet%h_DdirichletValues)
  IF (rdiscreteFBCDirichlet%h_IdirichletDOFs .NE. ST_NOHANDLE) &
    CALL storage_free(rdiscreteFBCDirichlet%h_IdirichletDOFs)

  END SUBROUTINE
  
END MODULE
