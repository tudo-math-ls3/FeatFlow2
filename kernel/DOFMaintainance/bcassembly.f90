!##############################################################################
!# ****************************************************************************
!# <name> bcassembly </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to discretise analytically given boundary
!# conditions. Analytically given boundary conditions are 'discretised', i.e.
!# a discrete version (realised by the structure t_discreteBCEntry) is calculated.
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
!#        structure
!#
!# 2.) bcasm_releaseDiscreteBC
!#     -> Cleans up a structure with discrete boundary conditions and
!#        releases all memory allocated by it.
!#
!# </purpose>
!##############################################################################

MODULE bcassembly

  USE fsystem
  USE discretebc
  USE collection
  USE triangulation
  USE dofmapping
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcasm_discretiseBC (rspatialDiscretisation,p_rdiscreteBC,bforceRebuild, &
                                 fgetBoundaryValues,rcollection)
  
!<description>
  ! This routine discretises an analytic definition of boundary conditions.
  ! The definition of the boundary conditions is taken from the discretisation
  ! structure rspatialDiscretisation. The discrete version is build up in
  ! p_rdiscreteBC. If p_rdiscreteBC is NULL(), a new structure is created,
  ! otherwise the old structure is updated (or even destroyed and recreaded if
  ! necessary).
!</description>

!<input>
  
  ! The discretisation structure of the underlying discretisation. The boundary
  ! conditions inside of this structure are discretised.
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rspatialDiscretisation
  
  ! Can be set to TRUE to force a complete rebuild of the rdiscreteBC structure.
  ! Normally, the structure is completely set up only in the first call
  ! or is there is a massive change in the boundary conditions (e.g. the
  ! number change)
  ! In later calls, the structure is recomputed only in those parts of the boundary,
  ! which have the t_bcReigion%bisstatic flag set FALSE (the standard value) -
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
!</input>

!<inputoutput>
  ! A discretised version of the analytic boundary conditions.
  ! This is a pointer to a t_discreteBC structures, 
  ! representing the boundary discretised in a discretisation-
  ! dependent way.
  ! If this pointer points to NULL(), a complete new structure is
  ! set up.
  TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: icurrentRegion
  TYPE(t_bcRegion), POINTER :: p_rbcRegion
  LOGICAL :: bbuildAll
  TYPE(t_collection), POINTER :: p_rcoll
  TYPE(t_discreteBCDirichlet), POINTER :: p_rdiscreteBCDirichlet
  
  ! Pointer to the boundary condition object
  TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions
  
  ! For quicker access:
  p_rboundaryConditions => rspatialDiscretisation%p_rboundaryConditions
  
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
            p_rdiscreteBCDirichlet => &
              p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%rdirichletBCs
            CALL bcasm_releaseDirichlet(p_rdiscreteBCDirichlet)

            ! BC released, indicate this
            p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%itype = DISCBC_TPUNDEFINED

          END SELECT
          
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
        CALL bcasm_discrBCDirichlet (rspatialDiscretisation, &
                   p_rbcRegion, p_rdiscreteBC%p_RdiscBCList(icurrentRegion), &
                   fgetBoundaryValues,p_rcoll)
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
  TYPE(t_discreteBCDirichlet), POINTER :: p_rdiscreteBCDirichlet
  
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
        p_rdiscreteBCDirichlet => &
          p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%rdirichletBCs
        CALL bcasm_releaseDirichlet(p_rdiscreteBCDirichlet)

        ! BC released, indicate this
        p_rdiscreteBC%p_RdiscBCList(icurrentRegion)%itype = DISCBC_TPUNDEFINED
      END SELECT
      
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
  ! which  index in the IverticesAtBoundary in the triangulation structure 
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
  !     Here is icount=2, IminIndex(1) / ImaxIndex(1) gives the
  !     first/last vertex of the first set and IminIndex(2) / ImaxIndex(2)
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
  !     Here is icount=2, IminIndex(1) / ImaxIndex(1) gives the
  !     first/last edge of the first set and IminIndex(2) / ImaxIndex(2)
  !     the first/last edge if the 2nd set.
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

  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcasm_discrBCDirichlet (rspatialDiscretisation, &
                                     rbcRegion, rdiscreteBC,fgetBoundaryValues,p_rcollection)
  
!<description>
  ! Creates a discrete version of Dirichlet boundary conditions.
  ! rbcRegion describes the region which is to be discretised. The discretised
  ! boundary conditions are created in rdiscreteBC, which is assumed
  ! to be undefined when entering this routine.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation. The boundary
  ! conditions inside of this structure are discretised.
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rspatialDiscretisation

  ! The BC region which is to be discrised
  TYPE(t_bcRegion), INTENT(IN) :: rbcRegion

  ! A callback function that calculates values on the boundary.
  ! Is declared in the interface include file 'intf_bcassembly.inc'.
  INCLUDE 'intf_bcassembly.inc'
  
  ! A collection structure to inform the callback function with
  ! additional information. Can be NULL() if there is no information to pass.
  TYPE(t_collection), POINTER :: p_rcollection
!</input>  

!<output>
  ! This structure receives the result of the discretisation of rbcRegion.
  ! When entering the routine, the content of this structure is undefined,
  ! all pointers are invalid. The routine fills everything with appropriate
  ! data.
  TYPE(t_discreteBCEntry), INTENT(OUT), TARGET :: rdiscreteBC
!</output>

  ! local variables
  INTEGER, DIMENSION(2) :: IminVertex,ImaxVertex,IminEdge,ImaxEdge,Iminidx,Imaxidx
  INTEGER :: i,ilocalEdge,ieltype,ielement,icount,icount2,ipart,j
  INTEGER(I32) :: iedge,ipoint1,ipoint2,NVT
  TYPE(t_discreteBCDirichlet),POINTER         :: p_rdirichletBCs
  TYPE(t_triangulation), POINTER            :: p_rtriangulation
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
  
  ! For easier access:
  p_rtriangulation => rspatialDiscretisation%p_rtriangulation
  CALL storage_getbase_int(p_rtriangulation%h_IelementsAtBoundary,p_IelementsAtBoundary)
  CALL storage_getbase_int(p_rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
  CALL storage_getbase_int(p_rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)
  CALL storage_getbase_int(rspatialDiscretisation%h_ItrialElements,p_ItrialElements)
  CALL storage_getbase_double(p_rtriangulation%h_DedgeParameterValue,p_DedgeParameterValue)
  CALL storage_getbase_double(p_rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
  CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
  CALL storage_getbase_int2D(p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
  NVT = p_rtriangulation%NVT
  
  ! We have Dirichlet boundary conditions
  rdiscreteBC%itype = DISCBC_TPDIRICHLET
  
  ! Connect to the boundary condition structure
  rdiscreteBC%p_rboundaryConditions => rspatialDiscretisation%p_rboundaryConditions
  
  ! Fill the structure for discrete Dirichlet BC's in the
  ! t_discreteBCEntry structure
  p_rdirichletBCs => rdiscreteBC%rdirichletBCs
  
  ! We have to deal with all DOF's on the boundary. This is highly element
  ! dependent and therefore a little bit tricky :(
  !
  ! As we are in 2D, we can use parameter values at first to figure out,
  ! which points and which edges are on the boundary.
  ! What we have is a boundary segment. Now ask the boundary-index routine
  ! to give us the minimum and maximum index of the vertices and edges on the
  ! bondary that belong to this boundary segment.
  
  CALL bcasm_getVertInBCregion (p_rtriangulation,rspatialDiscretisation%p_rdomain, &
                                rbcRegion%rboundaryRegion, &
                                IminVertex,ImaxVertex,icount)
  CALL bcasm_getEdgesInBCregion (p_rtriangulation,rspatialDiscretisation%p_rdomain, &
                                 rbcRegion%rboundaryRegion, &
                                 IminEdge,ImaxEdge,icount2)
                                 
  ! Cancel if the set is empty!
  IF ((icount .EQ. 0) .AND. (icount2 .EQ. 0)) THEN
    RETURN
  END IF
                                 
  ! The min/max values of these give the indices of the elements on the
  ! boundary that share our BC region - because it may be that an edge belongs
  ! to the BC region while the endpoints do not and vice versa, so we
  ! have to make sure to get an index which covers everything.
  Iminidx = MIN(IminVertex,IminEdge)
  Imaxidx = MAX(ImaxVertex,ImaxEdge)
  
  ! Reserve some memory to save temporarily all DOF's of all boundary
  ! elements.
  ! We handle all boundary elements simultaneously - let's hope that there are 
  ! never so many elements on the boundary that our memory runs out :-)
  ALLOCATE (Idofs(EL_MAXNBAS,ImaxIdx(1)-IminIdx(1)+1+ImaxIdx(2)-IminIdx(2)+1))
  ALLOCATE (DdofValue(EL_MAXNBAS,ImaxIdx(1)-IminIdx(1)+1+ImaxIdx(2)-IminIdx(2)+1))
  Idofs = 0

  ! Now the elements with indices iminidx..imaxidx in the ItrialElements
  ! of the triangulation are on the boundary. Some elements may appear
  ! twice (on edges e.g.) but we don't care.
  !
  ! Ask the DOF-mapping routine to get us those DOF's belonging to elements
  ! on the boundary.
  IF (icount .GE. 1) THEN
    CALL dof_locGlobMapping_mult(rspatialDiscretisation, &
              p_IelementsAtBoundary(Iminidx(1):Imaxidx(1)), .FALSE., Idofs)
  END IF
            
  IF (icount .GE. 2) THEN
    CALL dof_locGlobMapping_mult(rspatialDiscretisation, &
              p_IelementsAtBoundary(Iminidx(2):Imaxidx(2)), .FALSE., &
              Idofs( :, Imaxidx(1)-Iminidx(1)+1+1 : ))
  END IF
                               
  ! Loop through the index sets
  DO ipart = 1,icount2
  
    ! Now the element-dependent part. For each element type, we have to
    ! figure out which DOF's are on the boundary!
    DO i=iminidx(ipart),imaxidx(ipart)
    
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
      
      ! Element type?
      ieltype = p_ItrialElements(p_IelementsAtBoundary(I))
      SELECT CASE (ieltype)
      
      CASE (EL_P0,EL_Q0)
        ! Nice element, only one DOF :-)
        ! Either the edge or an adjacent vertex is on the boundary.
        ! Get the value at the corner point and accept it as
        ! Dirichlet value.
        CALL fgetBoundaryValues (rspatialDiscretisation,rbcRegion,ielement, &
                                DISCBC_NEEDFUNC,ipoint1,p_DvertexParameterValue(I), &
                                p_rcollection, Dvalues)
                                 
        ! Dvalues(1) gives us the function value in the point. Save it
        DdofValue(1,I-Iminidx(ipart)+1) = Dvalues(1)
        
        ! Set the DOF number < 0 to indicate that it is Dirichlet
        Idofs(1,I-Iminidx(ipart)+1) = -ABS(Idofs(1,I-Iminidx(ipart)+1))
        
      CASE (EL_P1,EL_Q1)

        ! Left point inside? -> Corresponding DOF must be computed
        IF ( (I .GE. IminVertex(ipart)) .AND. (I .LE. ImaxVertex(ipart)) ) THEN
          CALL fgetBoundaryValues (rspatialDiscretisation,rbcRegion,ielement, &
                                  DISCBC_NEEDFUNC,ipoint1,p_DvertexParameterValue(I), &
                                  p_rcollection, Dvalues)
                                  
          ! Save the computed function value
          DdofValue(ilocalEdge,I-Iminidx(ipart)+1) = Dvalues(1) 
          
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
        
      CASE (EL_EM30,EL_EM31,EL_E030,EL_E031)

        ! Edge inside? -> Calculate integral mean value over the edge
        IF ( (I .GE. IminEdge(ipart)) .AND. (I .LE. ImaxEdge(ipart)) ) THEN
          CALL fgetBoundaryValues (rspatialDiscretisation,rbcRegion,ielement, &
                                  DISCBC_NEEDINTMEAN,ipoint1,p_DvertexParameterValue(I), &
                                  p_rcollection, Dvalues)
                                  
          ! Save the computed function value
          DdofValue(ilocalEdge,I-Iminidx(ipart)+1) = Dvalues(1) 
          
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
  
    ! Allocate arrays for storing these DOF's and their values
    CALL storage_new('bcasm_discrBCDirichlet', 'h_IdirichletDOFs', &
                    icount, ST_INT, p_rdirichletBCs%h_IdirichletDOFs, &
                    ST_NEWBLOCK_NOINIT)
    CALL storage_new('bcasm_discrBCDirichlet', 'h_IdirichletValues', & 
                    icount, ST_DOUBLE, p_rdirichletBCs%h_IdirichletValues, &
                    ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int(p_rdirichletBCs%h_IdirichletDOFs,p_IdirichletDOFs)
    CALL storage_getbase_double(p_rdirichletBCs%h_IdirichletValues,p_IdirichletValues)
    
    ! Transfer the DOF's and their values to these arrays.
    icount = 0
    DO J=1,SIZE(Idofs,2)
      DO I=1,SIZE(Idofs,1)
        IF (Idofs(I,J) < 0) THEN
          icount = icount + 1
          p_IdirichletDOFs(icount) = ABS(Idofs(I,J))
          p_IdirichletValues(icount) = DdofValue(I,J)
        END IF
      END DO
    END DO
  
  ELSE
  
    ! Let's hope there is nothing saved here :-)
    p_rdirichletBCs%h_IdirichletDOFs = ST_NOHANDLE
    p_rdirichletBCs%h_IdirichletValues = ST_NOHANDLE
  
  END IF
  
  ! Remove temporary memory, finish.
  DEALLOCATE (DdofValue)
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

  ! Release what is associated
  
  rdiscreteBCDirichlet%nDOF = 0
  IF (rdiscreteBCDirichlet%h_IdirichletValues .NE. ST_NOHANDLE) &
    CALL storage_free(rdiscreteBCDirichlet%h_IdirichletValues)
  IF (rdiscreteBCDirichlet%h_IdirichletDOFs .NE. ST_NOHANDLE) &
    CALL storage_free(rdiscreteBCDirichlet%h_IdirichletDOFs)

  END SUBROUTINE
  

END MODULE
  