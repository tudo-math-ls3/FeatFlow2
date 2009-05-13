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

module triasearch

  use fsystem
  use storage
  use basicgeometry
  use triangulation
  use geometryaux
  
  implicit none
  
  private

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  interface tsrch_getElem_BruteForce
    module procedure tsrch_getElem_BruteForce_dir
    module procedure tsrch_getElem_BruteForce_ind
  end interface
  
  public :: tsrch_getElem_BruteForce

  interface tsrch_getNearestElem_BruteForce
    module procedure tsrch_getNearestMPElem_BF_dir
    module procedure tsrch_getNearestMPElem_BF_ind
  end interface
  
  public :: tsrch_getNearestElem_BruteForce

  interface tsrch_getElem_raytrace2D
    module procedure tsrch_getElem_raytrace2d_dir
    module procedure tsrch_getElem_raytrace2d_ind
  end interface

  public :: tsrch_getElem_raytrace2D

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

contains

!************************************************************************

!<subroutine>

  pure subroutine tsrch_getElem_BruteForce_dir (&
      Dpoint,DvertexCoords,IverticesAtElement,iel)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Brute force method (very slow on large triangulations!).
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  real(DP), dimension(:), intent(IN) :: Dpoint
  
  ! Vertex coordinates of the triangulation
  real(DP), dimension(:,:), intent(IN) :: DvertexCoords
  
  ! Vertices Adjacent to the elements in the triangulation.
  integer, dimension(:,:), intent(IN) :: IverticesAtElement
!</input>
  
!<output>
  ! Number of the element containing the point Dpoint.
  ! =0 if no element contains the point.
  integer, intent(OUT) :: iel
!</output>
  
!</subroutine>
    
    ! local variables
    real(DP), dimension(NDIM2D,TRIA_MAXNVE) :: Dcorners
    logical :: binside
  
    select case (ubound(Dpoint,1))
    case (NDIM2D) 
      
      ! Loop through all elements. Check if the element contains the point.
      ! If yes, quit.
      do iel=1,ubound(IverticesAtElement,2)
      
        ! Triangular or quad element?
        if (IverticesAtElement(4,iel) .ne. 0) then
          ! Fetch the coordinates of the element
          Dcorners(1:2,1:4) = DvertexCoords(1:2,IverticesAtElement(1:4,iel))
          
          ! Check if the point is inside. If yes, quit.
          call gaux_isInElement_quad2D(Dpoint(1),Dpoint(2),Dcorners,binside)
          if (binside) return
        else
          ! Fetch the coordinates of the element
          Dcorners(1:2,1:3) = DvertexCoords(1:2,IverticesAtElement(1:3,iel))
          
          ! Check if the point is inside. If yes, quit.
          call gaux_isInElement_tri2D(Dpoint(1),Dpoint(2),Dcorners,binside)
          if (binside) return
        end if
      
      end do
      
    end select
  
    ! No element found
    iel = 0
  
  end subroutine

!************************************************************************

!<subroutine>

  subroutine tsrch_getElem_BruteForce_ind (Dpoint,rtriangulation,iel)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Brute force method (very slow on large triangulations!).
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  real(DP), dimension(:), intent(IN) :: Dpoint
  
  ! Triangulation structure.
  type(t_triangulation), intent(IN) :: rtriangulation
!</input>
  
!<output>
  ! Number of the element containing the point Dpoint.
  ! =0 if no element contains the point.
  integer, intent(OUT) :: iel
!</output>
  
!</subroutine>
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    
    call tsrch_getElem_BruteForce_dir (Dpoint,p_DvertexCoords,p_IverticesAtElement,iel)
  
  end subroutine

!************************************************************************

!<subroutine>

  pure subroutine tsrch_getNearestMPElem_BF_dir (&
      Dpoint,DvertexCoords,IverticesAtElement,iel)
  
!<description>
  ! Searches for an element in the triangulation which element midpoint
  ! is closest to the point Dpoint.
  ! Brute force method (very slow on large triangulations!).
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  real(DP), dimension(:), intent(IN) :: Dpoint
  
  ! Vertex coordinates of the triangulation
  real(DP), dimension(:,:), intent(IN) :: DvertexCoords
  
  ! Vertices Adjacent to the elements in the triangulation.
  integer, dimension(:,:), intent(IN) :: IverticesAtElement
!</input>
  
!<output>
  ! Number of the element which midpoint is closest to Dpoint.
  integer, intent(OUT) :: iel
!</output>
  
!</subroutine>
    
    ! local variables
    real(DP), dimension(ubound(DvertexCoords,1)) :: Dmidpoint
    integer :: ive
    integer :: ielcur
    real(DP) :: ddist,dmindist
  
    select case (ubound(Dpoint,1))
    case (NDIM2D) 
      
      dmindist = -1.0_DP
      
      ! Loop through all elements. 
      do ielcur=1,ubound(IverticesAtElement,2)
      
        ! Get the element midpoint
        Dmidpoint(:) = 0.0_DP
        do ive = 1,ubound(IverticesAtElement,1)
          if (IverticesAtElement(ive,ielcur) .eq. 0) exit ! Triangle in a quad mesh
          Dmidpoint(:) = Dmidpoint(:) + &
            DvertexCoords(:,IverticesAtElement(ive,ielcur))
        end do
        Dmidpoint(:) = Dmidpoint(:) / real(ive-1,DP)
        
        ! Get the distance
        ddist = sqrt((Dmidpoint(1)-Dpoint(1))**2 + &
                (Dmidpoint(2)-Dpoint(2))**2)
        
        ! Chech if the distance is smaller than the current one.
        if ((dmindist .lt. 0.0_DP) .or. (ddist .lt. dmindist)) then 
          dmindist = ddist
          iel = ielcur
        end if
        
      end do
      
    end select
  
  end subroutine

!************************************************************************

!<subroutine>

  subroutine tsrch_getNearestMPElem_BF_ind (Dpoint,rtriangulation,iel)
  
!<description>
  ! Searches for an element in the triangulation which element midpoint
  ! is closest to the point Dpoint.
  ! Brute force method (very slow on large triangulations!).
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  real(DP), dimension(:), intent(IN) :: Dpoint
  
  ! Triangulation structure.
  type(t_triangulation), intent(IN) :: rtriangulation
!</input>
  
!<output>
  ! Number of the element containing the point Dpoint.
  ! =0 if no element contains the point.
  integer, intent(OUT) :: iel
!</output>
  
!</subroutine>
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    call storage_getbase_double2d (&
        rtriangulation%h_DvertexCoords,p_DvertexCoords)
    call storage_getbase_int2d (&
        rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    
    call tsrch_getNearestMPElem_BF_dir (&
        Dpoint,p_DvertexCoords,p_IverticesAtElement,iel)
  
  end subroutine

!************************************************************************

!<subroutine>

  subroutine tsrch_getElem_raytrace2D_dir (&
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
  real(DP), dimension(:), intent(IN) :: Dpoint
  
  ! Triangulation structure.
  type(t_triangulation), intent(IN) :: rtriangulation
  
  ! OPTIONAL: Maximum number of cell layers to pass. Default value is 100.
  integer, intent(IN), optional :: imaxIterations

  ! Vertex coordinates of the mesh
  real(DP), dimension(:,:) :: DvertexCoords
  
  ! Vertices adjacent to an element in the mesh
  integer, dimension(:,:) :: IverticesAtElement
  
  ! Edges adjacent to an element in the mesh
  integer, dimension(:,:) :: IedgesAtElement
  
  ! Neighbours around an element in the mesh
  integer, dimension(:,:) :: IneighboursAtElement

!</input>

!<inputoutput>
  ! On input: Number of an element where to start the search.
  !           =0: Choose element automatically.
  ! On output: Number of the element containing the point Dpoint.
  !            =0 if the element was not found.
  integer, intent(INOUT) :: iel
!</inputoutput>
  
!<output>
  ! OPTIONAL: Result of the search.
  ! =1 : The element was found successfully.
  ! =0 : The raytracing search broke down inside of the domain. 
  ! =-1: The search broke down because the domain was left.
  integer, intent(OUT), optional :: iresult
  
  ! OPTIONAL: Last analysed element.
  ! If iresult= 1: ilastElement = iel
  ! If iresult= 0: Number of the last analysed element before the search
  !                was stopped.
  ! If iresult=-1: Number of the element through which the
  !                domain was left. 
  integer, intent(OUT), optional :: ilastElement

  ! OPTIONAL: Number of the last analysed edge. Range 1..NMT.
  ! If iresult= 1: ilastEdge=0
  ! If iresult= 0: Number of the last analysed edge before the search
  !                was stopped.
  ! If iresult=-1: Number of the edge through which the domain was left. 
  integer, intent(OUT), optional :: ilastEdge
!</output>
  
!</subroutine>
    
    integer :: imaxIter,ite
    integer :: ielold
    real(DP), dimension(2,4) :: DcornerCoords
    integer :: ive,nnve
    logical :: bcheck
    real(DP) :: dxmid,dymid
    real(DP) :: dx1,dy1,dx2,dy2
    integer, dimension(4), parameter :: Inext = (/2,3,4,1/)
    
    ! Maximum number of iterations
    imaxIter = 100
    if (present(imaxIterations)) imaxIter = imaxIterations

    ! Start element. If the start element is not known, we start at 1.
    if (iel .le. 0) iel = 1
  
    ielold = iel
    
    nnve = ubound(IverticesAtElement,1)
    
    ! If the point is in element iel, we are immediately done.
    
    
    do ive = 1,nnve
      if (IverticesAtElement(ive,iel) .eq. 0) exit ! Triangle in a quad mesh
      DcornerCoords(1,ive) = DvertexCoords(1,IverticesAtElement(ive,iel))
      DcornerCoords(2,ive) = DvertexCoords(2,IverticesAtElement(ive,iel))
    end do
    
    ! We restrict our raytracing search to imaxIter neighbour cells;
    ! let's hope a point is not moved more than imaxIter elements in 
    ! one step...
    ieliteration: do ite = 1,imaxIter
    
      ! Check if the point is in the element.
      ! bcheck will be set to TRUE if that's the case.
      if (ive .eq. TRIA_NVETRI2D+1) then
        call gaux_isInElement_tri2D(Dpoint(1),Dpoint(2),DcornerCoords,bcheck)
      else
        call gaux_isInElement_quad2D(Dpoint(1),Dpoint(2),DcornerCoords,bcheck)
      end if
      
      if (bcheck) then
        ! We found the element.
        if (present(iresult)) iresult = 1
        if (present(ilastElement)) ilastElement = iel
        if (present(ilastEdge)) ilastEdge = 0
        return
      end if
      
      ! Calculate a point inside of the element. We start the ray from there.
      ! For simplicity (assuming convexity of the element) we take half the
      ! line from vertex 1 to vertex 3.
      dxmid = 0.5_DP * (DvertexCoords(1,IverticesAtElement(1,iel)) + &
                        DvertexCoords(1,IverticesAtElement(3,iel)) )
      dymid = 0.5_DP * (DvertexCoords(2,IverticesAtElement(1,iel)) + &
                        DvertexCoords(2,IverticesAtElement(3,iel)) )
                        
      ! Check all edges to find an intersection point:
      do ive = 1,nnve
      
        if (IverticesAtElement(ive,iel) .eq. 0) exit ! Triangle in a quad mesh
        
        ! Don't jump back to the old element
        if (IneighboursAtElement(ive,iel) .ne. ielold) then
          
          ! Calculate point of intersection the edge ive.
          ! (The use of Inext saves one division by avoiding MOD)
          dx1 = DvertexCoords(1,IverticesAtElement(ive,iel))
          dy1 = DvertexCoords(2,IverticesAtElement(ive,iel))
          dx2 = DvertexCoords(1,IverticesAtElement(Inext(ive),iel))
          dy2 = DvertexCoords(2,IverticesAtElement(Inext(ive),iel))
          
          call gaux_isIntersection_line2D (dx1,dy1,dx2,dy2,dxmid,dymid,&
              Dpoint(1),Dpoint(2),bcheck)
          
          if (bcheck) then
          
            ! Go on searching in the neighbour element
            ielold = iel
            
            ! Stop here if we leave our domain. To test that use the
            ! fact that KADJ()=0 if there is no neighbour element!
            ! The caller can use IEL to find out the information
            ! about the boundary where the domain was left...
            
            if (IneighboursAtElement(ive,iel) .eq. 0) then
            
              if (present(iresult)) iresult = -1
              if (present(ilastElement)) ilastElement = iel
              if (present(ilastEdge)) then
                ilastEdge = IedgesAtElement(ive,iel)
              end if
              
              iel = 0
              return
            
            end if
            
            ! Continue with the next element
            iel = IneighboursAtElement(ive,iel)
            cycle ieliteration
            
          end if
          
        end if
        
      end do 
    
      ! No face was found through which we could left the element.
      ! So the element contains the whole line (dxmid,dymid)->(dx,dy) --
      ! and so indeed the point (dx,dy)! That's it.
      ! iel already has the element number, we only have to set the optional
      ! argiments.
      if (present(iresult)) iresult = 0
      if (present(ilastElement)) ilastElement = ielold
      if (present(ilastEdge)) ilastEdge = 0
      return
    
    end do ieliteration
  
  end subroutine

!************************************************************************

!<subroutine>

  subroutine tsrch_getElem_raytrace2D_ind (&
      Dpoint,rtriangulation,iel, iresult,ilastElement,ilastEdge,imaxIterations)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Ray trace method for 2D meshes.
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  real(DP), dimension(:), intent(IN) :: Dpoint
  
  ! Triangulation structure.
  type(t_triangulation), intent(IN) :: rtriangulation
  
  ! OPTIONAL: Maximum number of cell layers to pass. Default value is 100.
  integer, intent(IN), optional :: imaxIterations
!</input>

!<inputoutput>
  ! On input: Number of an element where to start the search.
  !           =0: Choose element automatically.
  ! On output: Number of the element containing the point Dpoint.
  !            =0 if the element was not found.
  integer, intent(INOUT) :: iel
!</inputoutput>
  
!<output>
  ! OPTIONAL: Result of the search.
  ! =1 : The element was found successfully.
  ! =0 : The raytracing search broke down inside of the domain. 
  ! =-1: The search broke down because the domain was left.
  integer, intent(OUT), optional :: iresult
  
  ! OPTIONAL: Last analysed element.
  ! If iresult= 1: ilastElement = iel
  ! If iresult= 0: Number of the last analysed element before the search
  !                was stopped.
  ! If iresult=-1: Number of the element through which the
  !                domain was left. 
  integer, intent(OUT), optional :: ilastElement

  ! OPTIONAL: Number of the last analysed edge. Range 1..NMT.
  ! If iresult= 1: ilastEdge=0
  ! If iresult= 0: Number of the last analysed edge before the search
  !                was stopped.
  ! If iresult=-1: Number of the edge through which the domain was left. 
  integer, intent(OUT), optional :: ilastEdge
!</output>
  
!</subroutine>
    
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    
    ! Get triangulation arrays
    call storage_getbase_double2d (rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2d (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    call storage_getbase_int2d (rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)

    call tsrch_getElem_raytrace2D_dir (&
        Dpoint,rtriangulation,iel, &
        p_DvertexCoords,p_IverticesAtElement,p_IedgesAtElement,p_IneighboursAtElement,&
        iresult,ilastElement,ilastEdge,imaxIterations)    
    
  end subroutine

!************************************************************************

!<subroutine>

  subroutine tsrch_getElem_hierarch (Dpoint,Rtriangulation,iel)
  
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
  real(DP), dimension(:), intent(IN) :: Dpoint
  
  ! A sequence of triangulation structures that should stem from starndard
  ! 2-level-ordered refinement. 
  type(t_triangulation), dimension(:), intent(IN) :: Rtriangulation
!</input>

!<inputoutput>
  ! Number of the element on the finest mesh containing the point Dpoint.
  ! =0 if the element was not found.
  integer, intent(OUT) :: iel
!</inputoutput>
  
!</subroutine>

    integer :: ilev
    
    iel = 0
    
    do ilev = 1,size(Rtriangulation)
      
      ! Try to make a raytracing search on the current level.
      select case (Rtriangulation(1)%ndim)
      case (NDIM2D)
        call tsrch_getElem_raytrace2D (&
            Dpoint,Rtriangulation(1),iel)
      end select
        
      if (iel .le. 0) then
        ! Oops, raytracing failed.
        ! Then start an exhaustive search through all elements.
        call tsrch_getElem_BruteForce_ind (Dpoint,Rtriangulation(ilev),iel)
        
        ! iel may be zero after that, as in "circle" domains it 
        ! may happen that the point on the higher level is in the domain 
        ! while the point on the lower level is not!
        ! Let's hope that the domain is somehow convex so that
        ! raytracing works in the next sweep...
        
      end if

      ! Note: In the next sweep, we start raytracing at element iel which
      ! was calculated on the coarse mesh. In usual 2-level ordered meshes,
      ! the element with number iel on the fine grid is one subelement
      ! of element iel on the coarse grid, so the search should be quick!
      !
      ! If other refinement strategies are used, let's hope element iel on 
      ! the fine grid is nearby Dpoint, so raytracing is not too expensive...
    
    end do
    
  end subroutine


end module
