!##############################################################################
!# ****************************************************************************
!# <name> boundaryaux </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains auxiliary routines for handling boundary components.
!#
!# It contains the following set of routines:
!#
!# 1.) bdraux_getElementsAtBoundary
!#     -> Calculates the list of elements at the boundary
!#
!# 2.) bdraux_getElementsAtRegion
!#     -> Calculates the list of elements at a boundary region.
!#
!# 3.) bdraux_getElementsAtBdrComp
!#     -> Calculates the list of elements at a boundary component.
!#
!# 4.) bdraux_getNELAtBoundary
!#     -> Calculates the number of element at the boundary
!#
!# 5.) bdraux_getNELAtRegion
!#     -> Calculates the number of element at a boundary region.
!#
!# 6.) bdraux_getNELAtBdrComp
!#     -> Calculates the number of element at a boundary component.
!# </purpose>
!##############################################################################
module boundaryaux
  
  use basicgeometry
  use boundary
  use fsystem
  use genoutput
  use spatialdiscretisation
  use storage
  use triangulation
  
  implicit none

  private

  public :: bdraux_getElementsAtBoundary
  public :: bdraux_getElementsAtRegion
  public :: bdraux_getElementsAtBdrComp
  public :: bdraux_getNELAtBoundary
  public :: bdraux_getNELAtRegion
  public :: bdraux_getNELAtBdrComp
  
contains

  !****************************************************************************

!<subroutine>

  subroutine bdraux_getElementsAtBoundary(rdiscretisation,&
      NELbdc, IelementList, IelementOrientation, DedgePosition,&
      celement, cparType)

!<description>
    ! This subroutine calculates the list of elements which are
    ! adjacent to the boundary. If the number of elements
    ! exceeds the size of the working arrays an error is thrown.
!</description>

!<input>
    ! A discretisation structure
    type(t_spatialdiscretisation), intent(in) :: rdiscretisation
    
    ! OPTIONAL: Element type to be considered. If not present, then
    ! all elements adjacent to the boundary are inserted into the list
    integer(I32), intent(in), optional :: celement

    ! OPTIONAL: Type of parametrisation to use.
    ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
    integer, intent(in), optional :: cparType
!</input>

!<output>
    ! The number of elements at the boundary region
    integer, intent(out) :: NELbdc

    ! The list of elements adjacent to the boundary region
    integer, dimension(:), intent(out) :: IelementList

    ! The orientation of elements adjacent to the boundary region
    integer, dimension(:), intent(out) :: IelementOrientation

    ! OPTIONAL: The start- and end-parameter values of the edges on
    ! the boundary region
    real(DP), dimension(:,:), intent(out), optional :: DedgePosition
!</output>
!</subroutine>

    ! local variables
    type(t_boundary), pointer :: p_rboundary
    type(t_boundaryRegion) :: rboundaryRegion
    integer :: ibdc,NEL,cpar

    cpar = BDR_PAR_01
    if (present(cparType)) cpar = cparType

    ! The discretisation must provide a boundary structure
    if (associated(rdiscretisation%p_rboundary)) then
      p_rboundary => rdiscretisation%p_rboundary
    else
      call output_line('Discretisation does not provide boundary structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdraux_getElementsAtBoundary')
      call sys_halt()
    end if
    
    NELbdc = 0
    ! Create a boundary region for each boundary component and call
    ! the calculation routine for that.
    do ibdc = 1,boundary_igetNBoundComp(p_rboundary)
      call boundary_createRegion (p_rboundary, ibdc, 0, rboundaryRegion, cpar)

      if (present(DedgePosition)) then
        call bdraux_getElementsAtRegion(rboundaryRegion, rdiscretisation,&
            NEL, IelementList(NELbdc+1:), IelementOrientation(NELbdc+1:),&
            DedgePosition(:,NELbdc+1:), celement, cpar)
      else
        call bdraux_getElementsAtRegion(rboundaryRegion, rdiscretisation,&
            NEL, IelementList(NELbdc+1:), IelementOrientation(NELbdc+1:),&
            celement=celement, cparType=cpar)
      end if

      NELbdc = NELbdc + NEL
    end do
    
  end subroutine bdraux_getElementsAtBoundary

  !****************************************************************************

!<subroutine>

  subroutine bdraux_getElementsAtRegion(rboundaryRegion, rdiscretisation,&
      NELbdc, IelementList, IelementOrientation, DedgePosition,&
      celement, cparType)

!<description>
    ! This subroutine calculates the list of elements which are
    ! adjacent to the boundary region. If the number of elements
    ! exceeds the size of the working arrays an error is thrown.
!</description>

!<input>
    ! A t_boundaryRegion specifying the boundary region where
    ! to calculate. 
    type(t_boundaryRegion), intent(in) :: rboundaryRegion

    ! A discretisation structure
    type(t_spatialdiscretisation), intent(in) :: rdiscretisation

    ! OPTIONAL: Element type to be considered. If not present, then
    ! all elements adjacent to the boundary are inserted into the list
    integer(I32), intent(in), optional :: celement

    ! OPTIONAL: Type of parametrisation to use.
    ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
    integer, intent(in), optional :: cparType
!</input>

!<output>
    ! The number of elements at the boundary region
    integer, intent(out) :: NELbdc

    ! The list of elements adjacent to the boundary region
    integer, dimension(:), intent(out) :: IelementList

    ! The orientation of elements adjacent to the boundary region
    integer, dimension(:), intent(out) :: IelementOrientation

    ! OPTIONAL: The start- and end-parameter values of the edges on
    ! the boundary region
    real(DP), dimension(:,:), intent(out), optional :: DedgePosition
!</output>
!</subroutine>

    ! local variables
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:), pointer :: p_DedgeParameterValue
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:), pointer :: p_IelementsAtBoundary
    integer, dimension(:), pointer :: p_IelementDistr
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: ibdc,iedge,iel,ilocaledge,cpar

    cpar = BDR_PAR_01
    if (present(cparType)) cpar = cparType

    ! Get some pointers and arrays for quicker access
    p_rtriangulation => rdiscretisation%p_rtriangulation
    p_RelementDistribution => rdiscretisation%RelementDistr

    call storage_getbase_int (p_rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_int (p_rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    call storage_getbase_int (p_rtriangulation%h_IelementsAtBoundary,&
        p_IelementsAtBoundary)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
        p_DvertexCoordinates)
    call storage_getbase_double (p_rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)
    call storage_getbase_double (p_rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)
    
    ! Boundary component?
    ibdc = rboundaryRegion%iboundCompIdx

    ! Number of elements on that boundary component?
    NELbdc = p_IboundaryCpIdx(ibdc+1)-p_IboundaryCpIdx(ibdc)

    if (NELbdc .gt. ubound(IelementList,1) .or.&
        NELbdc .gt. ubound(IelementOrientation,1)) then
      call output_line('Insufficient memory',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdraux_getElementsAtRegion')
      call sys_halt()
    end if

    if (present(DedgePosition)) then
      if (NELbdc .gt. ubound(DedgePosition,2)) then
        call output_line('Insufficient memory',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdraux_getElementsAtRegion')
        call sys_halt()
      end if
    end if
    
    ! Do we have a uniform discretisation? Would simplify a lot...
    if (rdiscretisation%ccomplexity.eq. SPDISC_UNIFORM) then

      ! Check if there are elements which satisfy the type of element
      if (present(celement)) then
        if (p_RelementDistribution(1)%celement .ne. celement) then
          NELBdc = 0
          return
        end if
      end if
      
      ! Loop through the edges on the boundary component ibdc. If the
      ! edge is inside, remember the element number and figure out the
      ! orientation of the edge, whereby iel counts the total number
      ! of elements in the region. Be careful: The last edge must be
      ! treated differently!
      iel = 0
      do iedge = p_IboundaryCpIdx(ibdc),p_IboundaryCpIdx(ibdc+1)-2
        if (boundary_isInRegion(rboundaryRegion, ibdc,&
            p_DedgeParameterValue(iedge))) then
          
          iel = iel + 1
          
          ! Element number
          IelementList(iel) = p_IelementsAtBoundary(iedge)
          
          ! Element orientation; i.e. the local number of the boundary edge 
          do ilocaledge = 1,ubound(p_IedgesAtElement,1)
            if (p_IedgesAtElement(ilocaledge, p_IelementsAtBoundary(iedge)) .eq. &
                p_IedgesAtBoundary(iedge)) exit
          end do
          IelementOrientation(iel) = ilocaledge
          
          if (present(DedgePosition)) then
            ! Save the start and end parameter values of the edge
            DedgePosition(1,iel) = p_DvertexParameterValue(iedge)
            DedgePosition(2,iel) = p_DvertexParameterValue(iedge+1)
          end if

        end if
      end do

      ! Convert parameter values in length parametrisation
      if (present(DedgePosition))&
          call boundary_convertParameterList(rdiscretisation%p_rboundary,&
          ibdc, DedgePosition(:,1:iel), DedgePosition(:,1:iel),&
          rboundaryRegion%cparType, cpar)
      
      ! Handle the last edge differently
      iedge = p_IboundaryCpIdx(ibdc+1)-1
      
      if (boundary_isInRegion(rboundaryRegion, ibdc,&
          p_DedgeParameterValue(iedge))) then
          
        iel = iel + 1
          
        ! Element number
        IelementList(iel) = p_IelementsAtBoundary(iedge)

        ! Element orientation; i.e. the local number of the boundary edge 
        do ilocaledge = 1,ubound(p_IedgesAtElement,1)
          if (p_IedgesAtElement(ilocaledge, p_IelementsAtBoundary(iedge)) .eq. &
              p_IedgesAtBoundary(iedge)) exit
        end do
        IelementOrientation(iel) = ilocaledge
        
        if (present(DedgePosition)) then
          ! Save the start parameter value of the edge -- in length parametrisation.
          DedgePosition(1,iel) = boundary_convertParameter(rdiscretisation%p_rboundary,&
              ibdc, p_DvertexParameterValue(iedge), rboundaryRegion%cparType, cpar)
          
          ! Save the end parameter value of the edge -- in length parametrisation.
          DedgePosition(2,iel) = boundary_dgetMaxParVal(rdiscretisation%p_rboundary,&
              ibdc, cpar)
        end if
      end if
      
    else
      
      ! Set pointer
      call storage_getbase_int (rdiscretisation%h_IelementDistr,&
          p_IelementDistr)

      ! Loop through the edges on the boundary component ibdc. If the
      ! edge is inside, remember the element number and figure out the
      ! orientation of the edge, whereby iel counts the total number
      ! of elements in the region. Be careful: The last edge must be
      ! treated differently!
      iel = 0
      do iedge = p_IboundaryCpIdx(ibdc),p_IboundaryCpIdx(ibdc+1)-2
        if (boundary_isInRegion(rboundaryRegion, ibdc,&
            p_DedgeParameterValue(iedge))) then
          
          ! Check if we are the correct element distribution
          if (present(celement)) then
            if (celement .ne. p_RelementDistribution(&
                p_IelementDistr(p_IelementsAtBoundary(iedge)))%celement) cycle
          end if
          iel = iel + 1

          ! Element number
          IelementList(iel) = p_IelementsAtBoundary(iedge)
          
          ! Element orientation; i.e. the local number of the boundary edge 
          do ilocaledge = 1,ubound(p_IedgesAtElement,1)
            if (p_IedgesAtElement(ilocaledge, p_IelementsAtBoundary(iedge)) .eq. &
                p_IedgesAtBoundary(iedge)) exit
          end do
          IelementOrientation(iel) = ilocaledge
          
          if (present(DedgePosition)) then
            ! Save the start and end parameter values of the edge
            DedgePosition(1,iel) = p_DvertexParameterValue(iedge)
            DedgePosition(2,iel) = p_DvertexParameterValue(iedge+1)
          end if

        end if
      end do

      if (present(DedgePosition)) then
        ! Convert parameter values in length parametrisation
        call boundary_convertParameterList(rdiscretisation%p_rboundary,&
            ibdc, DedgePosition(:,1:iel), DedgePosition(:,1:iel),&
            rboundaryRegion%cparType, cpar)
      end if

      ! Handle the last edge differently
      iedge = p_IboundaryCpIdx(ibdc+1)-1
      
      if (boundary_isInRegion(rboundaryRegion, ibdc,&
          p_DedgeParameterValue(iedge))) then
        
        ! Check if we are the correct element distribution
        if (present(celement)) then
          if (celement .ne. p_RelementDistribution(&
              p_IelementDistr(p_IelementsAtBoundary(iedge)))%celement) then
            ! Adjust number of elements
            NELbdc = iel
            return
          end if
        end if
        iel = iel + 1
        
        ! Element number
        IelementList(iel) = p_IelementsAtBoundary(iedge)
        
        ! Element orientation; i.e. the local number of the boundary edge 
        do ilocaledge = 1,ubound(p_IedgesAtElement,1)
          if (p_IedgesAtElement(ilocaledge, p_IelementsAtBoundary(iedge)) .eq. &
              p_IedgesAtBoundary(iedge)) exit
        end do
        IelementOrientation(iel) = ilocaledge
        
        if (present(DedgePosition)) then
          ! Save the start parameter value of the edge -- in length parametrisation.
          DedgePosition(1,iel) = boundary_convertParameter(rdiscretisation%p_rboundary, &
              ibdc, p_DvertexParameterValue(iedge), rboundaryRegion%cparType, cpar)
          
          ! Save the end parameter value of the edge -- in length parametrisation.
          DedgePosition(2,iel) = boundary_dgetMaxParVal(rdiscretisation%p_rboundary,&
              ibdc, cpar)
        end if

      end if

    end if

    ! Adjust number of elements
    NELbdc = iel
        
  end subroutine bdraux_getElementsAtRegion

  !****************************************************************************

!<subroutine>

  subroutine bdraux_getElementsAtBdrComp(ibdc, rdiscretisation,&
      NELbdc, IelementList, IelementOrientation, DedgePosition,&
      celement, cparType)

!<description>
    ! This subroutine calculates the list of elements which are
    ! adjacent to the boundary component ibdrComp. If the number of
    ! elements exceeds the size of the working arrays an error is
    ! thrown.
!</description>

!<input>
    ! Number of the boundary component
    integer, intent(in) :: ibdc

    ! A discretisation structure
    type(t_spatialdiscretisation), intent(in) :: rdiscretisation

    ! OPTIONAL: Element type to be considered. If not present, then
    ! all elements adjacent to the boundary are inserted into the list
    integer(I32), intent(in), optional :: celement

    ! OPTIONAL: Type of parametrisation to use.
    ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
    integer, intent(in), optional :: cparType
!</input>

!<output>
    ! The number of elements at the boundary region
    integer, intent(out) :: NELbdc

    ! The list of elements adjacent to the boundary region
    integer, dimension(:), intent(out) :: IelementList

    ! The orientation of elements adjacent to the boundary region
    integer, dimension(:), intent(out) :: IelementOrientation

    ! OPTIONAL: The start- and end-parameter values of the edges on
    ! the boundary region
    real(DP), dimension(:,:), intent(out), optional :: DedgePosition
!</output>
!</subroutine>

    ! local variables
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IelementsAtBoundary
    integer, dimension(:), pointer :: p_IelementDistr
    integer, dimension(:), pointer :: p_InodalProperty
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: iel,iedge,cpar

    cpar = BDR_PAR_01
    if (present(cparType)) cpar = cparType

    select case(rdiscretisation%ndimension)
    case(NDIM1D)
      ! Get some pointers and arrays for quicker access
      p_rtriangulation => rdiscretisation%p_rtriangulation
      p_RelementDistribution => rdiscretisation%RelementDistr

      call storage_getbase_int (p_rtriangulation%h_IboundaryCpIdx,&
          p_IboundaryCpIdx)
      call storage_getbase_int (p_rtriangulation%h_IelementsAtBoundary,&
          p_IelementsAtBoundary)
      call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
          p_IverticesAtElement)
      call storage_getbase_int (p_rtriangulation%h_InodalProperty,&
        p_InodalProperty)

      ! Number of elements on that boundary component?
      NELbdc = p_IboundaryCpIdx(ibdc+1)-p_IboundaryCpIdx(ibdc)

      if (NELbdc .gt. ubound(IelementList,1) .or.&
          NELbdc .gt. ubound(IelementOrientation,1)) then
        call output_line('Insufficient memory',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdraux_getElementsAtBdrComp')
        call sys_halt()
      end if

      ! Do we have a uniform discretisation? Would simplify a lot...
      if (rdiscretisation%ccomplexity.eq. SPDISC_UNIFORM) then
        
        ! Check if there are elements which satisfy the type of element
        if (present(celement)) then
          if (p_RelementDistribution(1)%celement .ne. celement) then
            NELBdc = 0
            return
          end if
        end if

        ! Determine the element numbers and their orientation at the boundary
        iel = 0
        do iedge = p_IboundaryCpIdx(ibdc),p_IboundaryCpIdx(ibdc+1)-1
          
          iel = iel+1
          
          ! Element number
          IelementList(iel) = p_IelementsAtBoundary(iedge)
          
          ! Element orientation, i.e. the local number of the boundary vertex
          if (p_InodalProperty(&
              p_IverticesAtElement(1,IelementList(iel))).eq. ibdc) then
            IelementOrientation(iel) = 1
          elseif (p_InodalProperty(&
              p_IverticesAtElement(2,IelementList(iel))).eq. ibdc) then
            IelementOrientation(iel) = 2
          else
            call output_line('Unable to determine element orientation!',&
                OU_CLASS_ERROR,OU_MODE_STD,'bdraux_getElementsAtBdrComp')
          end if
        end do
        
      else

        ! Set pointer
        call storage_getbase_int (rdiscretisation%h_IelementDistr,&
            p_IelementDistr)
        
        ! Determine the element numbers and their orientation at the boundary
        iel = 0
        do iedge = p_IboundaryCpIdx(ibdc),p_IboundaryCpIdx(ibdc+1)-1

          ! Check if we are the correct element distribution
          if (present(celement)) then
            if (celement .ne. p_RelementDistribution(&
                p_IelementDistr(p_IelementsAtBoundary(iedge)))%celement) cycle
          end if
          iel = iel+1
          
          ! Element number
          IelementList(iel) = p_IelementsAtBoundary(iedge)
          
          ! Element orientation, i.e. the local number of the boundary vertex
          if (p_InodalProperty(&
              p_IverticesAtElement(1,IelementList(iel))).eq. ibdc) then
            IelementOrientation(iel) = 1
          elseif (p_InodalProperty(&
              p_IverticesAtElement(2,IelementList(iel))).eq. ibdc) then
            IelementOrientation(iel) = 2
          else
            call output_line('Unable to determine element orientation!',&
                OU_CLASS_ERROR,OU_MODE_STD,'bdraux_getElementsAtBdrComp')
          end if
        end do
        
      end if
      

    case(NDIM2D)
      ! The discretisation must provide a boundary structure
      if (associated(rdiscretisation%p_rboundary)) then
        p_rboundary => rdiscretisation%p_rboundary
      else
        call output_line('Discretisation does not provide boundary structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bdraux_getElementsAtBdrComp')
        call sys_halt()
      end if
      
      ! Create region for entire boundary and determine elements
      call boundary_createRegion (p_rboundary, ibdc, 0, rboundaryRegion)
      call bdraux_getElementsAtRegion(rboundaryRegion, rdiscretisation, NELbdc,&
          IelementList, IelementOrientation, DedgePosition, celement, cpar)


    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bdraux_getElementsAtBdrComp')
      call sys_halt()
    end select

  end subroutine bdraux_getElementsAtBdrComp

  !****************************************************************************

!<function>

  function bdraux_getNELAtBoundary(rboundary, rtriangulation) result (NELbdc)

!<description>
    ! This function calculates the number of elements which are
    ! adjacent to the boundary using the given triangulation.
!</description>

!<input>
    ! The boundary for which the number of elements is to be calculated
    type(t_boundary), intent(in) :: rboundary

    ! The triangulation structure
    type(t_triangulation), intent(in) :: rtriangulation
!</input>

!<result>
    ! The number of elements at the boundary region
    integer :: NELbdc
!</result>
!</function>

    ! local variables
    integer :: ibdc
    
    NELbdc = 0
    ! Create a boundary region for each boundary component and call
    ! the calculation routine for that.
    do ibdc = 1,boundary_igetNBoundComp(rboundary)
      NELbdc = NELbdc+bdraux_getNELAtBdrComp(ibdc, rtriangulation)
    end do
    
  end function bdraux_getNELAtBoundary

  !****************************************************************************

!<function>

  function bdraux_getNELAtRegion(rboundaryRegion, rtriangulation) result(NELbdc)

!<description>
    ! This function calculates an upper bound for the number of
    ! elements which are adjacent to the boundary using the given
    ! triangulation. 
!</description>

!<input>
    ! A t_boundaryRegion specifying the boundary region where
    ! to calculate. 
    type(t_boundaryRegion), intent(in) :: rboundaryRegion

    ! The triangulation structure
    type(t_triangulation), intent(in) :: rtriangulation
!</input>

!<result>
    ! The number of elements at the boundary region
    integer :: NELbdc
!</result>
!</function>

    ! local variables
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer :: ibdc

    ! Set pointer
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
    
    ! Boundary component?
    ibdc = rboundaryRegion%iboundCompIdx

    ! Number of elements on that boundary component?
    NELbdc = p_IboundaryCpIdx(ibdc+1)-p_IboundaryCpIdx(ibdc)
    
  end function bdraux_getNELAtRegion

  !****************************************************************************

!<function>

  function bdraux_getNELAtBdrComp(ibdc, rtriangulation) result(NELbdc)

!<description>
    ! This function calculates the number of elements which are
    ! adjacent to the boundary component using the given triangulation.
!</description>

!<input>
    ! Number of the boundary component where to calculate
    integer, intent(in) :: ibdc

    ! The triangulation structure
    type(t_triangulation), intent(in) :: rtriangulation
!</input>

!<result>
    ! The number of elements at the boundary region
    integer :: NELbdc
!</result>
!</function>

    ! local variables
    integer, dimension(:), pointer :: p_IboundaryCpIdx

    ! Set pointer
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
    
    ! Number of elements on that boundary component?
    NELbdc = p_IboundaryCpIdx(ibdc+1)-p_IboundaryCpIdx(ibdc)
    
  end function bdraux_getNELAtBdrComp
  
end module boundaryaux
