!##############################################################################
!# ****************************************************************************
!# <name> triasearch </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains various search routines for triangulations,
!# E.g. to search for an element inside of a triangulation.
!#
!# The following routines can be found here:
!#
!# 1.) tsrch_getElem_BruteForce
!#     -> Search for an element containing a specific point.
!#        Brute force method.
!#
!# 2.) tsrch_getElem_raytrace2D
!#     -> Search for an element containing a specific point.
!#        Raytracing method in 2D.
!#
!# 3.) tsrch_getElem_hierarch
!#     -> Search for an element containing a specific point.
!#        Hierarchical method.
!#
!# 4.) tsrch_getNearestElem_BruteForce
!#     -> Find the element with the nearest midpoint to a specific point.
!#
!# </purpose>
!##############################################################################

MODULE triasearch

  USE triangulation
  USE geometryaux
  
  IMPLICIT NONE

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  INTERFACE tsrch_getElem_BruteForce
    MODULE PROCEDURE tsrch_getElem_BruteForce_dir
    MODULE PROCEDURE tsrch_getElem_BruteForce_ind
  END INTERFACE

  INTERFACE tsrch_getNearestElem_BruteForce
    MODULE PROCEDURE tsrch_getNearestMPElem_BF_dir
    MODULE PROCEDURE tsrch_getNearestMPElem_BF_ind
  END INTERFACE

  INTERFACE tsrch_getElem_raytrace2D
    MODULE PROCEDURE tsrch_getElem_raytrace2d_dir
    MODULE PROCEDURE tsrch_getElem_raytrace2d_ind
  END INTERFACE

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

CONTAINS

!************************************************************************

!<subroutine>

  PURE SUBROUTINE tsrch_getElem_BruteForce_dir (&
      Dpoint,DvertexCoords,IverticesAtElement,iel)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Brute force method (very slow on large triangulations!).
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dpoint
  
  ! Vertex coordinates of the triangulation
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: DvertexCoords
  
  ! Vertices Adjacent to the elements in the triangulation.
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement
!</input>
  
!<output>
  ! Number of the element containing the point Dpoint.
  ! =0 if no element contains the point.
  INTEGER(PREC_ELEMENTIDX), INTENT(OUT) :: iel
!</output>
  
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(NDIM2D,TRIA_MAXNVE) :: Dcorners
    LOGICAL :: binside
  
    SELECT CASE (UBOUND(Dpoint,1))
    CASE (NDIM2D) 
      
      ! Loop through all elements. Check if the element contains the point.
      ! If yes, quit.
      DO iel=1,UBOUND(IverticesAtElement,2)
      
        ! Triangular or quad element?
        IF (IverticesAtElement(4,iel) .NE. 0) THEN
          ! Fetch the coordinates of the element
          Dcorners(1:2,1:4) = DvertexCoords(1:2,IverticesAtElement(1:4,iel))
          
          ! Check if the point is inside. If yes, quit.
          CALL gaux_isInElement_quad2D(Dpoint(1),Dpoint(2),Dcorners,binside)
          IF (binside) RETURN
        ELSE
          ! Fetch the coordinates of the element
          Dcorners(1:2,1:3) = DvertexCoords(1:2,IverticesAtElement(1:3,iel))
          
          ! Check if the point is inside. If yes, quit.
          CALL gaux_isInElement_tri2D(Dpoint(1),Dpoint(2),Dcorners,binside)
          IF (binside) RETURN
        END IF
      
      END DO
      
    END SELECT
  
    ! No element found
    iel = 0
  
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tsrch_getElem_BruteForce_ind (Dpoint,rtriangulation,iel)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Brute force method (very slow on large triangulations!).
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dpoint
  
  ! Triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
!</input>
  
!<output>
  ! Number of the element containing the point Dpoint.
  ! =0 if no element contains the point.
  INTEGER(PREC_ELEMENTIDX), INTENT(OUT) :: iel
!</output>
  
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement

    CALL storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
    CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    
    CALL tsrch_getElem_BruteForce_dir (Dpoint,p_DvertexCoords,p_IverticesAtElement,iel)
  
  END SUBROUTINE

!************************************************************************

!<subroutine>

  PURE SUBROUTINE tsrch_getNearestMPElem_BF_dir (&
      Dpoint,DvertexCoords,IverticesAtElement,iel)
  
!<description>
  ! Searches for an element in the triangulation which element midpoint
  ! is closest to the point Dpoint.
  ! Brute force method (very slow on large triangulations!).
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dpoint
  
  ! Vertex coordinates of the triangulation
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: DvertexCoords
  
  ! Vertices Adjacent to the elements in the triangulation.
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement
!</input>
  
!<output>
  ! Number of the element which midpoint is closest to Dpoint.
  INTEGER(PREC_ELEMENTIDX), INTENT(OUT) :: iel
!</output>
  
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(UBOUND(DvertexCoords,1)) :: Dmidpoint
    INTEGER :: ive
    INTEGER(PREC_ELEMENTIDX) :: ielcur
    REAL(DP) :: ddist,dmindist
  
    SELECT CASE (UBOUND(Dpoint,1))
    CASE (NDIM2D) 
      
      dmindist = -1.0_DP
      
      ! Loop through all elements. 
      DO ielcur=1,UBOUND(IverticesAtElement,2)
      
        ! Get the element midpoint
        Dmidpoint(:) = 0.0_DP
        DO ive = 1,UBOUND(IverticesAtElement,1)
          IF (IverticesAtElement(ive,ielcur) .EQ. 0) EXIT ! Triangle in a quad mesh
          Dmidpoint(:) = Dmidpoint(:) + &
            DvertexCoords(:,IverticesAtElement(ive,ielcur))
        END DO
        Dmidpoint(:) = Dmidpoint(:) / REAL(ive-1,DP)
        
        ! Get the distance
        ddist = SQRT((Dmidpoint(1)-Dpoint(1))**2 + &
                (Dmidpoint(2)-Dpoint(2))**2)
        
        ! Chech if the distance is smaller than the current one.
        IF ((dmindist .LT. 0.0_DP) .OR. (ddist .LT. dmindist)) THEN 
          dmindist = ddist
          iel = ielcur
        END IF
        
      END DO
      
    END SELECT
  
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tsrch_getNearestMPElem_BF_ind (Dpoint,rtriangulation,iel)
  
!<description>
  ! Searches for an element in the triangulation which element midpoint
  ! is closest to the point Dpoint.
  ! Brute force method (very slow on large triangulations!).
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dpoint
  
  ! Triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
!</input>
  
!<output>
  ! Number of the element containing the point Dpoint.
  ! =0 if no element contains the point.
  INTEGER(PREC_ELEMENTIDX), INTENT(OUT) :: iel
!</output>
  
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement

    CALL storage_getbase_double2d (&
        rtriangulation%h_DvertexCoords,p_DvertexCoords)
    CALL storage_getbase_int2d (&
        rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    
    CALL tsrch_getNearestMPElem_BF_dir (&
        Dpoint,p_DvertexCoords,p_IverticesAtElement,iel)
  
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tsrch_getElem_raytrace2D_dir (&
      Dpoint,rtriangulation,iel, &
      DvertexCoords,IverticesAtElement,IedgesAtElement,IneighboursAtElement,&
      iresult,ilastElement,ilastEdge,imaxIterations)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Ray trace method for 2D meshes.
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dpoint
  
  ! Triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! OPTIONAL: Maximum number of cell layers to pass. Default value is 100.
  INTEGER, INTENT(IN), OPTIONAL :: imaxIterations

  ! Vertex coordinates of the mesh
  REAL(DP), DIMENSION(:,:) :: DvertexCoords
  
  ! Vertices adjacent to an element in the mesh
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:) :: IverticesAtElement
  
  ! Edges adjacent to an element in the mesh
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:) :: IedgesAtElement
  
  ! Neighbours around an element in the mesh
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:) :: IneighboursAtElement

!</input>

!<inputoutput>
  ! On input: Number of an element where to start the search.
  !           =0: Choose element automatically.
  ! On output: Number of the element containing the point Dpoint.
  !            =0 if the element was not found.
  INTEGER(PREC_ELEMENTIDX), INTENT(INOUT) :: iel
!</inputoutput>
  
!<output>
  ! OPTIONAL: Result of the search.
  ! =1 : The element was found successfully.
  ! =0 : The raytracing search broke down inside of the domain. 
  ! =-1: The search broke down because the domain was left.
  INTEGER, INTENT(OUT), OPTIONAL :: iresult
  
  ! OPTIONAL: Last analysed element.
  ! If iresult= 1: ilastElement = iel
  ! If iresult= 0: Number of the last analysed element before the search
  !                was stopped.
  ! If iresult=-1: Number of the element through which the
  !                domain was left. 
  INTEGER, INTENT(OUT), OPTIONAL :: ilastElement

  ! OPTIONAL: Number of the last analysed edge. Range 1..NMT.
  ! If iresult= 1: ilastEdge=0
  ! If iresult= 0: Number of the last analysed edge before the search
  !                was stopped.
  ! If iresult=-1: Number of the edge through which the domain was left. 
  INTEGER, INTENT(OUT), OPTIONAL :: ilastEdge
!</output>
  
!</subroutine>
    
    INTEGER :: imaxIter,ite
    INTEGER(PREC_ELEMENTIDX) :: ielold
    REAL(DP), DIMENSION(2,4) :: DcornerCoords
    INTEGER :: ive,nnve
    LOGICAL :: bcheck
    REAL(DP) :: dxmid,dymid
    REAL(DP) :: dx1,dy1,dx2,dy2
    INTEGER, DIMENSION(4), PARAMETER :: Inext = (/2,3,4,1/)
    
    ! Maximum number of iterations
    imaxIter = 100
    IF (PRESENT(imaxIterations)) imaxIter = imaxIterations

    ! Start element. If the start element is not known, we start at 1.
    IF (iel .LE. 0) iel = 1
  
    ielold = iel
    
    nnve = UBOUND(IverticesAtElement,1)
    
    ! If the point is in element iel, we are immediately done.
    
    DO ive = 1,nnve
      IF (IverticesAtElement(ive,iel) .EQ. 0) EXIT ! Triangle in a quad mesh
      DcornerCoords(1,ive) = DvertexCoords(1,IverticesAtElement(ive,iel))
      DcornerCoords(2,ive) = DvertexCoords(2,IverticesAtElement(ive,iel))
    END DO
    IF (ive .EQ. TRIA_NVETRI2D+1) THEN
      CALL gaux_isInElement_tri2D(Dpoint(1),Dpoint(2),DcornerCoords,bcheck)
    ELSE
      CALL gaux_isInElement_quad2D(Dpoint(1),Dpoint(2),DcornerCoords,bcheck)
    END IF
    
    IF (bcheck) THEN
      IF (PRESENT(iresult)) iresult = 1
      IF (PRESENT(ilastElement)) ilastElement = iel
      IF (PRESENT(ilastEdge)) ilastEdge = 0
      RETURN
    END IF
    
    ! We restrict our raytracing search to imaxIter neighbour cells;
    ! let's hope a point is not moved more than imaxIter elements in 
    ! one step...
    ieliteration: DO ite = 1,imaxIter
    
      ! Calculate a point inside of the element. We start the ray from there.
      ! For simplicity (assuming convexity of the element) we take half the
      ! line from vertex 1 to vertex 3.
      dxmid = 0.5_DP * (DvertexCoords(1,IverticesAtElement(1,iel)) + &
                        DvertexCoords(1,IverticesAtElement(3,iel)) )
      dymid = 0.5_DP * (DvertexCoords(2,IverticesAtElement(1,iel)) + &
                        DvertexCoords(2,IverticesAtElement(3,iel)) )
                        
      ! Check all edges to find an intersection point:
      DO ive = 1,nnve
      
        IF (IverticesAtElement(ive,iel) .EQ. 0) EXIT ! Triangle in a quad mesh
        
        ! Don't jump back to the old element
        IF (IneighboursAtElement(ive,iel) .NE. ielold) THEN
          
          ! Calculate point of intersection the edge ive.
          ! (The use of Inext saves one division by avoiding MOD)
          dx1 = DvertexCoords(1,IverticesAtElement(ive,iel))
          dy1 = DvertexCoords(2,IverticesAtElement(ive,iel))
          dx2 = DvertexCoords(1,IverticesAtElement(Inext(ive),iel))
          dy2 = DvertexCoords(2,IverticesAtElement(Inext(ive),iel))
          
          CALL gaux_isIntersection_line2D (dx1,dy1,dx2,dy2,dxmid,dymid,&
              Dpoint(1),Dpoint(2),bcheck)
          
          IF (bcheck) THEN
          
            ! Go on searching in the neighbour element
            ielold = iel
            
            ! Stop here if we leave our domain. To test that use the
            ! fact that KADJ()=0 if there is no neighbour element!
            ! The caller can use IEL to find out the information
            ! about the boundary where the domain was left...
            
            IF (IneighboursAtElement(ive,iel) .EQ. 0) THEN
            
              IF (PRESENT(iresult)) iresult = -1
              IF (PRESENT(ilastElement)) ilastElement = iel
              IF (PRESENT(ilastEdge)) THEN
                ilastEdge = IedgesAtElement(ive,iel)
              END IF
              
              iel = 0
              RETURN
            
            END IF
            
            ! Continue with the next element
            iel = IneighboursAtElement(ive,iel)
            CYCLE ieliteration
            
          END IF
          
        END IF
        
      END DO 
    
      ! No face was found through which we could left the element.
      ! So the element contains the whole line (dxmid,dymid)->(dx,dy) --
      ! and so indeed the point (dx,dy)! That's it.
      ! iel already has the element number, we only have to set the optional
      ! argiments.
      IF (PRESENT(iresult)) iresult = 0
      IF (PRESENT(ilastElement)) ilastElement = ielold
      IF (PRESENT(ilastEdge)) ilastEdge = 0
      RETURN
    
    END DO ieliteration
  
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tsrch_getElem_raytrace2D_ind (&
      Dpoint,rtriangulation,iel, iresult,ilastElement,ilastEdge,imaxIterations)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Ray trace method for 2D meshes.
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dpoint
  
  ! Triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! OPTIONAL: Maximum number of cell layers to pass. Default value is 100.
  INTEGER, INTENT(IN), OPTIONAL :: imaxIterations
!</input>

!<inputoutput>
  ! On input: Number of an element where to start the search.
  !           =0: Choose element automatically.
  ! On output: Number of the element containing the point Dpoint.
  !            =0 if the element was not found.
  INTEGER(PREC_ELEMENTIDX), INTENT(INOUT) :: iel
!</inputoutput>
  
!<output>
  ! OPTIONAL: Result of the search.
  ! =1 : The element was found successfully.
  ! =0 : The raytracing search broke down inside of the domain. 
  ! =-1: The search broke down because the domain was left.
  INTEGER, INTENT(OUT), OPTIONAL :: iresult
  
  ! OPTIONAL: Last analysed element.
  ! If iresult= 1: ilastElement = iel
  ! If iresult= 0: Number of the last analysed element before the search
  !                was stopped.
  ! If iresult=-1: Number of the element through which the
  !                domain was left. 
  INTEGER, INTENT(OUT), OPTIONAL :: ilastElement

  ! OPTIONAL: Number of the last analysed edge. Range 1..NMT.
  ! If iresult= 1: ilastEdge=0
  ! If iresult= 0: Number of the last analysed edge before the search
  !                was stopped.
  ! If iresult=-1: Number of the edge through which the domain was left. 
  INTEGER, INTENT(OUT), OPTIONAL :: ilastEdge
!</output>
  
!</subroutine>
    
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
    
    ! Get triangulation arrays
    CALL storage_getbase_double2d (rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    CALL storage_getbase_int2d (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    CALL storage_getbase_int2d (rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)

    CALL tsrch_getElem_raytrace2D_dir (&
        Dpoint,rtriangulation,iel, &
        p_DvertexCoords,p_IverticesAtElement,p_IedgesAtElement,p_IneighboursAtElement,&
        iresult,ilastElement,ilastEdge,imaxIterations)    
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tsrch_getElem_hierarch (Dpoint,Rtriangulation,iel)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Hierarchical element search for 2D meshes.
  ! The routine accepts a sequence of meshes Rtriangulation and a point
  ! Dpoint. In the finest mesh in Rtriangulation, that element is searched
  ! which contains the point Dpoint.
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dpoint
  
  ! A sequence of triangulation structures that should stem from starndard
  ! 2-level-ordered refinement. 
  TYPE(t_triangulation), DIMENSION(:), INTENT(IN) :: Rtriangulation
!</input>

!<inputoutput>
  ! Number of the element on the finest mesh containing the point Dpoint.
  ! =0 if the element was not found.
  INTEGER(PREC_ELEMENTIDX), INTENT(OUT) :: iel
!</inputoutput>
  
!</subroutine>

    INTEGER :: ilev
    
    iel = 0
    
    DO ilev = 1,SIZE(Rtriangulation)
      
      ! Try to make a raytracing search on the current level.
      SELECT CASE (Rtriangulation(1)%ndim)
      CASE (NDIM2D)
        CALL tsrch_getElem_raytrace2D (&
            Dpoint,Rtriangulation(1),iel)
      END SELECT
        
      IF (iel .LE. 0) THEN
        ! Oops, raytracing failed.
        ! Then start an exhaustive search through all elements.
        CALL tsrch_getElem_BruteForce_ind (Dpoint,Rtriangulation(ilev),iel)
        
        ! iel may be zero after that, as in "circle" domains it 
        ! may happen that the point on the higher level is in the domain 
        ! while the point on the lower level is not!
        ! Let's hope that the domain is somehow convex so that
        ! raytracing works in the next sweep...
        
      END IF

      ! Note: In the next sweep, we start raytracing at element iel which
      ! was calculated on the coarse mesh. In usual 2-level ordered meshes,
      ! the element with number iel on the fine grid is one subelement
      ! of element iel on the coarse grid, so the search should be quick!
      !
      ! If other refinement strategies are used, let's hope element iel on 
      ! the fine grid is nearby Dpoint, so raytracing is not too expensive...
    
    END DO
    
  END SUBROUTINE


END MODULE
