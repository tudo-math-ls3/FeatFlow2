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
!# 2.) tsrch_getElem_raytrace2D, tsrch_getElem_raytrace3D
!#     -> Search for an element containing a specific point.
!#        Raytracing method in 2D/3D.
!#
!# 3.) tsrch_getElem_hierarch
!#     -> Search for an element containing a specific point.
!#        Hierarchical method.
!#
!# 4.) tsrch_getNearestElem_BruteForce
!#     -> Find the element with the nearest midpoint to a specific point.
!#
!# 5.) tsrch_getElementsByRaytrace
!#     -> Determins all elements containing a set of points
!#        using raytracing search
!# </purpose>
!##############################################################################

module triasearch

  use fsystem
  use storage
  use basicgeometry
  use triangulation
  use geometryaux
  use genoutput
  
  ! DEBUG!!!
  ! use io
  ! use geometryoutput
  
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
  
  interface tsrch_getElem_raytrace3D
    module procedure tsrch_getElem_raytrace3d_dir
    module procedure tsrch_getElem_raytrace3d_ind
  end interface
  
  public :: tsrch_getElem_raytrace3D
  public :: tsrch_getElementsByRaytrace  

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

contains

!************************************************************************

!<subroutine>

 subroutine tsrch_getElem_BruteForce_dir (&
      Dpoint,DvertexCoords,IverticesAtElement,iel)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Brute force method (very slow on large triangulations!).
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  real(DP), dimension(:), intent(in) :: Dpoint
  
  ! Vertex coordinates of the triangulation
  real(DP), dimension(:,:), intent(in) :: DvertexCoords
  
  ! Vertices Adjacent to the elements in the triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement
!</input>
  
!<output>
  ! Number of the element containing the point Dpoint.
  ! =0 if no element contains the point.
  integer, intent(out) :: iel
!</output>
  
!</subroutine>
    
    ! local variables
    !
    ! WARNING: Use different dimensions for 2D and 3D for more performance
    ! since the brute force search methods expect arrays in a definite shape!
    real(DP), dimension(NDIM2D,TRIA_MAXNVE2D) :: Dcorners2D
    real(DP), dimension(NDIM3D,TRIA_MAXNVE3D) :: Dcorners3D
    logical :: binside
  
    select case (ubound(Dpoint,1))
    case (NDIM2D) 
      
      ! Loop through all elements. Check if the element contains the point.
      ! If yes, quit.
      do iel=1,ubound(IverticesAtElement,2)
      
      ! Check if pure triangle grid
        if (size(IverticesAtElement(:,1)) .eq. 3) then
       
            Dcorners2D(1:2,1:3) = DvertexCoords(1:2,IverticesAtElement(1:3,iel))

            ! Check if the point is inside. If yes, quit.
            call gaux_isInElement_tri2D(Dpoint(1),Dpoint(2),Dcorners2D,binside)
            if (binside) return

        else
      
        ! Triangular or quad element?
        if (IverticesAtElement(4,iel) .ne. 0) then
          ! Fetch the coordinates of the element
          Dcorners2D(1:2,1:4) = DvertexCoords(1:2,IverticesAtElement(1:4,iel))
          
          ! Check if the point is inside. If yes, quit.
          call gaux_isInElement_quad2D(Dpoint(1),Dpoint(2),Dcorners2D,binside)
          if (binside) return
        else
          ! Fetch the coordinates of the element
          Dcorners2D(1:2,1:3) = DvertexCoords(1:2,IverticesAtElement(1:3,iel))
          
          ! Check if the point is inside. If yes, quit.
          call gaux_isInElement_tri2D(Dpoint(1),Dpoint(2),Dcorners2D,binside)
          if (binside) return
        end if
      
      end if
      
      end do
    case (NDIM3D)
    
      ! Loop through all elements. Check if the element contains the point.
      ! If yes, quit.
      do iel=1,ubound(IverticesAtElement,2)
        ! Fetch the coordinates of the element
        Dcorners3D(1:3,1) = DvertexCoords(1:3,IverticesAtElement(1,iel))
        Dcorners3D(1:3,2) = DvertexCoords(1:3,IverticesAtElement(2,iel))
        Dcorners3D(1:3,3) = DvertexCoords(1:3,IverticesAtElement(3,iel))
        Dcorners3D(1:3,4) = DvertexCoords(1:3,IverticesAtElement(4,iel))
        Dcorners3D(1:3,5) = DvertexCoords(1:3,IverticesAtElement(5,iel))
        Dcorners3D(1:3,6) = DvertexCoords(1:3,IverticesAtElement(6,iel))
        Dcorners3D(1:3,7) = DvertexCoords(1:3,IverticesAtElement(7,iel))
        Dcorners3D(1:3,8) = DvertexCoords(1:3,IverticesAtElement(8,iel))

        call gaux_isInElement_hexa(Dpoint(1),Dpoint(2),Dpoint(3),Dcorners3D,binside)
        if (binside) return
      
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
  real(DP), dimension(:), intent(in) :: Dpoint
  
  ! Triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
!</input>
  
!<output>
  ! Number of the element containing the point Dpoint.
  ! =0 if no element contains the point.
  integer, intent(out) :: iel
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
  real(DP), dimension(:), intent(in) :: Dpoint
  
  ! Vertex coordinates of the triangulation
  real(DP), dimension(:,:), intent(in) :: DvertexCoords
  
  ! Vertices Adjacent to the elements in the triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement
!</input>
  
!<output>
  ! Number of the element which midpoint is closest to Dpoint.
  integer, intent(out) :: iel
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
  real(DP), dimension(:), intent(in) :: Dpoint
  
  ! Triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
!</input>
  
!<output>
  ! Number of the element containing the point Dpoint.
  ! =0 if no element contains the point.
  integer, intent(out) :: iel
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
  real(DP), dimension(:), intent(in) :: Dpoint
  
  ! Triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! OPTIONAL: Maximum number of cell layers to pass. Default value is 100.
  integer, intent(in), optional :: imaxIterations

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
  integer, intent(inout) :: iel
!</inputoutput>
  
!<output>
  ! OPTIONAL: Result of the search.
  ! =1 : The element was found successfully.
  ! =0 : The raytracing search broke down inside of the domain. 
  ! =-1: The search broke down because the domain was left.
  ! =-2: The maximum number of iterations was exceeded
  integer, intent(out), optional :: iresult
  
  ! OPTIONAL: Last analysed element.
  ! If iresult= 1: ilastElement = iel
  ! If iresult= 0: Number of the last analysed element before the search
  !                was stopped.
  ! If iresult=-1: Number of the element through which the
  !                domain was left. 
  ! If iresult=-2: Number of the last element in the raytracing search
  integer, intent(out), optional :: ilastElement

  ! OPTIONAL: Number of the last analysed edge. Range 1..NMT.
  ! If iresult= 1: ilastEdge=0
  ! If iresult= 0: Number of the last analysed edge before the search
  !                was stopped.
  ! If iresult=-1: Number of the edge through which the domain was left. 
  ! If iresult=-2: Number of the last edge 
  integer, intent(out), optional :: ilastEdge
!</output>
  
!</subroutine>
    
    integer :: imaxIter,ite
    integer :: ielold
    real(DP), dimension(2,4) :: DcornerCoords
    integer :: ive,nve,nnve
    logical :: bcheck
    real(DP) :: dxmid,dymid
    real(DP) :: dx1,dy1,dx2,dy2
    
    ! Maximum number of iterations
    imaxIter = 100
    if (present(imaxIterations)) imaxIter = imaxIterations

    ! Start element. If the start element is not known, we start at 1.
    if (iel .le. 0) iel = 1
  
    ielold = iel
    
    nnve = ubound(IverticesAtElement,1)
    
    ! DEBUG!!!
    ! call io_deleteFile('overview.dat')
    ! call geoout_writeGnuplotPoint (Dpoint,0,'overview.dat')
    
    ! We restrict our raytracing search to imaxIter neighbour cells;
    ! let us hope a point is not moved more than imaxIter elements in 
    ! one step...
    ieliteration: do ite = 1,imaxIter
    
      ! Fetch the element
      do ive = 1,nnve
        if (IverticesAtElement(ive,iel) .eq. 0) exit ! Triangle in a quad mesh
        DcornerCoords(1,ive) = DvertexCoords(1,IverticesAtElement(ive,iel))
        DcornerCoords(2,ive) = DvertexCoords(2,IverticesAtElement(ive,iel))
      end do
  
      ! Store number of vertices for current element
      nve = ive-1
    
      ! Check if the point is in the element.
      ! bcheck will be set to TRUE if that is the case.
      if (nve .eq. TRIA_NVETRI2D) then
        call gaux_isInElement_tri2D(Dpoint(1),Dpoint(2),DcornerCoords,bcheck)

        ! DEBUG!!!
        ! call geoout_writeGnuplotTria2D (DcornerCoords,0,'overview.dat')
      else
        call gaux_isInElement_quad2D(Dpoint(1),Dpoint(2),DcornerCoords,bcheck)

        ! DEBUG!!!
        ! call geoout_writeGnuplotQuad2D (DcornerCoords,0,'overview.dat')
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
        
        ! Do not jump back to the old element
        if (IneighboursAtElement(ive,iel) .ne. ielold) then
          
          ! Calculate point of intersection the edge ive.
          dx1 = DvertexCoords(1,IverticesAtElement(ive,iel))
          dy1 = DvertexCoords(2,IverticesAtElement(ive,iel))
          dx2 = DvertexCoords(1,IverticesAtElement(mod(ive, nve)+1,iel))
          dy2 = DvertexCoords(2,IverticesAtElement(mod(ive, nve)+1,iel))

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
    
    
    end do ieliteration
    
    ! No face was found through which we could left the element.
    ! So the element contains the whole line (dxmid,dymid)->(dx,dy) --
    ! and so indeed the point (dx,dy)! That is it.
    ! iel already has the element number, we only have to set the optional
    ! argiments.
    iel = 0
    if (present(iresult)) iresult = -2
    if (present(ilastElement)) ilastElement = ielold
    if (present(ilastEdge)) ilastEdge = 0
  
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
  real(DP), dimension(:), intent(in) :: Dpoint
  
  ! Triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! OPTIONAL: Maximum number of cell layers to pass. Default value is 100.
  integer, intent(in), optional :: imaxIterations
!</input>

!<inputoutput>
  ! On input: Number of an element where to start the search.
  !           =0: Choose element automatically.
  ! On output: Number of the element containing the point Dpoint.
  !            =0 if the element was not found.
  integer, intent(inout) :: iel
!</inputoutput>
  
!<output>
  ! OPTIONAL: Result of the search.
  ! =1 : The element was found successfully.
  ! =0 : The raytracing search broke down inside of the domain. 
  ! =-1: The search broke down because the domain was left.
  ! =-2: The maximum number of iterations was exceeded
  integer, intent(out), optional :: iresult
  
  ! OPTIONAL: Last analysed element.
  ! If iresult= 1: ilastElement = iel
  ! If iresult= 0: Number of the last analysed element before the search
  !                was stopped.
  ! If iresult=-1: Number of the element through which the
  !                domain was left. 
  ! If iresult=-2: Number of the last element in the raytracing search
  integer, intent(out), optional :: ilastElement

  ! OPTIONAL: Number of the last analysed edge. Range 1..NMT.
  ! If iresult= 1: ilastEdge=0
  ! If iresult= 0: Number of the last analysed edge before the search
  !                was stopped.
  ! If iresult=-1: Number of the edge through which the domain was left. 
  ! If iresult=-2: Number of the last element in the raytracing search
  integer, intent(out), optional :: ilastEdge
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
  real(DP), dimension(:), intent(in) :: Dpoint
  
  ! A sequence of triangulation structures that should stem from starndard
  ! 2-level-ordered refinement. 
  type(t_triangulation), dimension(:), intent(in) :: Rtriangulation
!</input>

!<inputoutput>
  ! Number of the element on the finest mesh containing the point Dpoint.
  ! =0 if the element was not found.
  integer, intent(out) :: iel
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
        ! Let us hope that the domain is somehow convex so that
        ! raytracing works in the next sweep...
        
      end if

      ! Note: In the next sweep, we start raytracing at element iel which
      ! was calculated on the coarse mesh. In usual 2-level ordered meshes,
      ! the element with number iel on the fine grid is one subelement
      ! of element iel on the coarse grid, so the search should be quick!
      !
      ! If other refinement strategies are used, let us hope element iel on 
      ! the fine grid is nearby Dpoint, so raytracing is not too expensive...
    
    end do
    
  end subroutine

!************************************************************************

!<subroutine>

  subroutine tsrch_getElem_raytrace3D_dir(&
      Dpoint,rtriangulation,iel, &
      DvertexCoords,IverticesAtElement,IfacesAtElement,IneighboursAtElement,&
      iresult,ilastElement,ilastFace,imaxIterations)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Ray trace method for 3D meshes. The setup is the same as in 2d
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  real(DP), dimension(:), intent(in) :: Dpoint
  
  ! Triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! OPTIONAL: Maximum number of cell layers to pass. Default value is 100.
  integer, intent(in), optional :: imaxIterations

!</input>

!<inputoutput>
  ! On input: Number of an element where to start the search.
  !           =0: Choose element automatically.
  ! On output: Number of the element containing the point Dpoint.
  !            =0 if the element was not found.
  integer, intent(inout) :: iel
!</inputoutput>
  
!<output>
  ! OPTIONAL: Result of the search.
  ! =1 : The element was found successfully.
  ! =0 : The raytracing search broke down inside of the domain. 
  ! =-1: The search broke down because the domain was left.
  ! =-2: The maximum number of iterations was exceeded
  integer, intent(out), optional :: iresult
  
  ! OPTIONAL: Last analysed element.
  ! If iresult= 1: ilastElement = iel
  ! If iresult= 0: Number of the last analysed element before the search
  !                was stopped.
  ! If iresult=-1: Number of the element through which the
  !                domain was left. 
  ! If iresult=-2: Number of the last element in the raytracing search
  integer, intent(out), optional :: ilastElement

  ! OPTIONAL: Number of the last analysed edge. Range 1..NMT.
  ! If iresult= 1: ilastFace=0
  ! If iresult= 0: Number of the last analysed face before the search
  !                was stopped.
  ! If iresult=-1: Number of the face through which the domain was left. 
  ! If iresult=-2: Number of the last face 
  integer, intent(out), optional :: ilastFace
!</output>
  
!</subroutine>
  ! Vertex coordinates of the mesh
  real(DP), dimension(:,:) :: DvertexCoords
  
  ! Vertices adjacent to an element in the mesh
  integer, dimension(:,:) :: IverticesAtElement
  
  ! faces of an element in the mesh
  integer, dimension(:,:) :: IfacesAtElement
  
  ! Neighbours around an element in the mesh
  integer, dimension(:,:) :: IneighboursAtElement
  
!  if (rtriangulation%h_IfacesAtElement .eq. ST_NOHANDLE) &
!      call tria_genFacesAtElement3D (rtriangulation)
!  
!  if (rtriangulation%h_IverticesAtFace .eq. ST_NOHANDLE) &
!      call tria_genVerticesAtFace3D (rtriangulation)          


    
    integer :: imaxIter,ite
    integer :: ielold
    real(DP), dimension(3,8) :: DcornerCoords
    integer :: ive,nnae,iae
    logical :: bcheck
    real(dp), dimension(3) :: Dmid
    real(dp), dimension(3,4) :: Dface
    integer, dimension(4,TRIA_NAEHEXA3D), parameter :: IverticesHexa =&
             reshape((/1,2,3,4, 1,5,6,2, 2,6,7,3,&
                       3,7,8,4, 1,4,8,5, 5,8,7,6/), (/4,TRIA_NAEHEXA3D/))
    
    
    ! Maximum number of iterations
    imaxIter = 100
    if (present(imaxIterations)) imaxIter = imaxIterations

    ! Start element. If the start element is not known, we start at 1.
    if (iel .le. 0) iel = 1
  
    ielold = iel
    
    nnae = ubound(IfacesAtElement,1)
    
    ! DEBUG!!!
    ! call io_deleteFile('overview.dat')
    ! call geoout_writeGnuplotPoint (Dpoint,0,'overview.dat')
    
    ! We restrict our raytracing search to imaxIter neighbour cells;
    ! let us hope a point is not moved more than imaxIter elements in 
    ! one step...
    ieliteration: do ite = 1,imaxIter
    
      ! Fetch the element
      do ive = 1,TRIA_NVEHEXA3D
        DcornerCoords(1:3,ive) = DvertexCoords(1:3,IverticesAtElement(ive,iel))
      end do
      
      ! Check if the point is in the element.
      ! bcheck will be set to TRUE if that is the case.
      call gaux_isInElement_hexa(Dpoint(1),Dpoint(2),Dpoint(3),DcornerCoords,bcheck)
      ! call geoout_writeGnuplotQuad2D (DcornerCoords,0,'overview.dat')
      if (bcheck) then
        ! We found the element.
        if (present(iresult)) iresult = 1
        if (present(ilastElement)) ilastElement = iel
        if (present(ilastFace)) ilastFace = 0
        return
      end if
      
      ! Calculate a point inside of the element. We start the ray from there.
      ! For simplicity we take the center point of the element
      Dmid(1:3) = DvertexCoords(1:3,IverticesAtElement(1,iel)) + &
                  DvertexCoords(1:3,IverticesAtElement(2,iel)) + &
                  DvertexCoords(1:3,IverticesAtElement(3,iel)) + &
                  DvertexCoords(1:3,IverticesAtElement(4,iel)) + &
                  DvertexCoords(1:3,IverticesAtElement(5,iel)) + &
                  DvertexCoords(1:3,IverticesAtElement(6,iel)) + &
                  DvertexCoords(1:3,IverticesAtElement(7,iel)) + &
                  DvertexCoords(1:3,IverticesAtElement(8,iel))
                  
      Dmid(:) = Dmid(:) * 0.125_dp                  
                        
      ! Check all faces to find an intersection point:
      do iae = 1,nnae
      
        do ive=1,TRIA_NVEQUAD2D
          Dface(1:3,ive)=DcornerCoords(1:3,IverticesHexa(ive,iae))
        end do
      
        ! Do not jump back to the old element
        if (IneighboursAtElement(iae,iel) .ne. ielold) then
          
          ! Calculate point of intersection with face
          call gaux_isIntersection_face (Dmid,&
               Dpoint,Dface, bcheck)
          
          if (bcheck) then
          
            ! Go on searching in the neighbour element
            ielold = iel
            
            ! Stop here if we leave our domain. To test that use the
            ! fact that KADJ()=0 if there is no neighbour element!
            ! The caller can use IEL to find out the information
            ! about the boundary where the domain was left...
            
            if (IneighboursAtElement(iae,iel) .eq. 0) then
              if (present(iresult)) iresult = -1
              if (present(ilastElement)) ilastElement = iel
              if (present(ilastFace)) then
                ilastFace = IfacesAtElement(iae,iel)
              end if
              iel = 0
              return
            end if ! if IneighboursAtElement(ive,iel) .eq. 0
            
            ! Continue with the next element
            iel = IneighboursAtElement(iae,iel)
            cycle ieliteration
            
          end if ! bcheck
          
        end if ! (IneighboursAtElement(ive,iel) .ne. ielold)
        
      end do ! for all faces
    
    
    end do ieliteration ! end for max search path length
    
    ! No face was found through which we could left the element.
    ! So the element contains the whole line (dxmid,dymid)->(dx,dy) --
    ! and so indeed the point (dx,dy)! That is it.
    ! iel already has the element number, we only have to set the optional
    ! argiments.
    iel = 0
    if (present(iresult)) iresult = -2
    if (present(ilastElement)) ilastElement = ielold
    if (present(ilastFace)) ilastFace = 0
  
  end subroutine
  
!************************************************************************

!<subroutine>

  subroutine tsrch_getElem_raytrace3D_ind (&
      Dpoint,rtriangulation,iel, iresult,ilastElement,ilastFace,imaxIterations)
  
!<description>
  ! Find an element in the triangulation containing the point Dpoint.
  ! Ray trace method for 3d meshes.
!</description>
  
!<input>
  ! The coordinate of a point where the routine should check which element
  ! contains the point.
  real(DP), dimension(:), intent(in) :: Dpoint
  
  ! Triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! OPTIONAL: Maximum number of cell layers to pass. Default value is 100.
  integer, intent(in), optional :: imaxIterations
!</input>

!<inputoutput>
  ! On input: Number of an element where to start the search.
  !           =0: Choose element automatically.
  ! On output: Number of the element containing the point Dpoint.
  !            =0 if the element was not found.
  integer, intent(inout) :: iel
!</inputoutput>
  
!<output>
  ! OPTIONAL: Result of the search.
  ! =1 : The element was found successfully.
  ! =0 : The raytracing search broke down inside of the domain. 
  ! =-1: The search broke down because the domain was left.
  integer, intent(out), optional :: iresult
  
  ! OPTIONAL: Last analysed element.
  ! If iresult= 1: ilastElement = iel
  ! If iresult= 0: Number of the last analysed element before the search
  !                was stopped.
  ! If iresult=-1: Number of the element through which the
  !                domain was left. 
  integer, intent(out), optional :: ilastElement

  ! OPTIONAL: Number of the last analysed Face. Range 1..NMT.
  ! If iresult= 1: ilastFace=0
  ! If iresult= 0: Number of the last analysed Face before the search
  !                was stopped.
  ! If iresult=-1: Number of the Face through which the domain was left. 
  integer, intent(out), optional :: ilastFace
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

    call tsrch_getElem_raytrace3d_dir (&
        Dpoint,rtriangulation,iel, &
        p_DvertexCoords,p_IverticesAtElement,p_IedgesAtElement,p_IneighboursAtElement,&
        iresult,ilastElement,ilastFace,imaxIterations)    
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine tsrch_getElementsByRaytrace (rtriangulation,Dpoints,Ielements,&
      ifirstUnknownPoint,ilastUnknownPoint)

!<description>
  ! Determine the elements containing the points in Dpoints using a raytracing
  ! linesearch.
!</description>

!<input>
  ! Underlying triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! List of points where to find the corresponding elements to.
  real(DP), dimension(:,:), intent(in) :: Dpoints
!</input>

!<output>  
  ! List of elements containing the points in Dpoints.
  ! If there was no element found for a point, the correpsonding element number
  ! is set to 0.
  integer, dimension(:), intent(out) :: Ielements
  
  ! Index in Dpoints to the first point where the corresponding element
  ! was not found. =0: All points found.
  ! If not all elements were found, the caller can loop from here through
  ! Ielements to get all unknown points.
  integer, intent(out) :: ifirstUnknownPoint

  ! Index in Dpoints to the last point where the corresponding element
  ! was not found. =0: All points found.
  ! If not all elements were found, the caller can loop from here through
  ! Ielements to get all unknown points.
  integer, intent(out) :: ilastUnknownPoint
!</output>

!</subroutine>

    integer :: iel,ipoint,ilastfoundel,iresult

    ! We use an accumulated raytracing search to find the points.
    ! This is very fast for lines and hopefully acceptable
    ! for arbitrary meshes.
    
    ifirstUnknownPoint = 0
    ilastUnknownPoint = 0
    
    ilastfoundel = 0
    do ipoint = 1,ubound(Dpoints,2)
    
      ! Find the element
      select case (ubound(Dpoints,1))
      case (NDIM2D)
      
        ! Try to start from the last found element.
        iel = ilastfoundel
        call tsrch_getElem_raytrace2D (Dpoints(:,ipoint),rtriangulation,iel,&
            iresult,ilastfoundel)
            
      case (NDIM3D)

        ! Try to start from the last found element.
        iel = ilastfoundel
        call tsrch_getElem_raytrace3D (Dpoints(:,ipoint),rtriangulation,iel,&
            iresult,ilastfoundel)

      case default
        call output_line ('Invalid dimension.', &
            OU_CLASS_ERROR,OU_MODE_STD,'tsrch_getElements')
        call sys_halt()
      end select
      
      if (iel .eq. 0) then
        ! Fallback to brute force search
        call tsrch_getElem_BruteForce (Dpoints(:,ipoint),rtriangulation,iel)
        
        ! Continue the next search from here.
        ! If not found, this is =0, so raytracing will start from the scratch.
        ilastfoundel = iel
      end if
      
      if (iel .eq. 0) then
        ! Oops, not found.
        if (ifirstUnknownPoint .eq. 0) then
          ifirstUnknownPoint = ipoint
        end if
        ilastUnknownPoint = ipoint
      end if
      
      ! Save the element
      Ielements(ipoint) = iel
      
    end do

  end subroutine

end module
