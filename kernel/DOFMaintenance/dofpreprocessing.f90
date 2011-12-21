!##############################################################################
!# ****************************************************************************
!# <name> dofpreprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic structures and routine for preparing a
!# subset of degrees of freedom which belong to a set of elements
!#
!# Routines in this module:
!#
!# 1.) dofprep_initDofSet
!#     -> Initialises a subset of degrees of freedoms based on a
!#        domain integration structure.
!#
!# 2.) dofprep_initDofSetAtBoundary
!#     -> Initialises a subset of degrees of freedoms based on a
!#        domain integration structure restricted to the boundary.
!#
!# 3.) dofprep_doneDofSet
!#     -> Releases a subset of degrees of freedom.
!#
!# </purpose>
!##############################################################################

module dofpreprocessing

  use element
  use fsystem
  use genoutput
  use domainintegration

  implicit none

  private

!<types>
  
!<typeblock>

  ! This structure is used to collect information about a set of
  ! degrees of freedom which belong to a subset of elements.
  type t_dofSubset

    ! Number of elements in the set
    integer :: nelements = 0
    
    ! Number of degrees of freedom per element
    integer :: ndofsPerElement = 0
    
    ! Element identifier
    integer(I32) :: celement = 0_I32

    ! Local ordering of the degrees of freedom
    !  DIMENSION(1..ndofsPerElement,1..nelements)
    integer, dimension(:,:), pointer :: p_IdofsLoc => null()

    ! Coordinates of the nodal degrees of freedom
    !  DIMENSION(1..ndim,1..ndofsPerElement,1..nelements)
    real(DP), dimension(:,:,:), pointer :: p_DdofCoords => null()
  end type t_dofSubset

  public :: t_dofSubset
  public :: dofprep_initDofSet
  public :: dofprep_initDofSetAtBoundary
  public :: dofprep_doneDofSet

!</typeblock>

!</types>

contains

  !*****************************************************************************

!<subroutine>

  subroutine dofprep_initDofSet(rdofSubset, rdomainIntSubset,&
      IelementOrientation)
    
!<description>
    ! This routine initialises a t_dofSubset structure using the
    ! information provided by the domain integration subset
    ! rdomainIntSubset.
!</description>

!<input>
    ! A domain integration subset
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: An array providing the orientation of the element.
    ! If not given, then the element orientation provided by
    ! rdomainIntSubset is adopted if available.
    integer, dimension(:), target, optional :: IelementOrientation
!</input>

!<output>
    ! A DOF subset structure, initialised according to the parameters
    type(t_dofSubset), intent(out) :: rdofSubset
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:), allocatable :: IdofsLoc
    integer, dimension(:), pointer :: p_IelementOrientation
    integer :: iel,idofe

    
    ! Get orientation of elements (if available)
    if (present(IelementOrientation)) then
      p_IelementOrientation => IelementOrientation
    elseif (associated(rdomainIntSubset%p_IelementOrientation)) then
      p_IelementOrientation => rdomainIntSubset%p_IelementOrientation
    else
      nullify(p_IelementOrientation)
    end if

    ! Initialise constants in the structure
    rdofSubset%nelements       = rdomainIntSubset%nelements
    rdofSubset%celement        = rdomainIntSubset%celement
    rdofSubset%ndofsPerElement = elem_igetNDofLoc(rdomainIntSubset%celement)

    ! Allocate temporal memory for local DOF numbers
    allocate(IdofsLoc(rdofSubset%ndofsPerElement))

    ! What type of finite element are we?
    select case(elem_getPrimaryElement(rdofSubset%celement))
    case (EL_P0,EL_Q0)
      IdofsLoc = (/1/)

    case (EL_P1,EL_P1T)
      IdofsLoc = (/1,2,3/)

    case (EL_Q1,EL_Q1T)
      IdofsLoc = (/1,2,3,4/)

    case (EL_P2)
      IdofsLoc = (/1,4,2,5,3,6/)

    case (EL_Q2)
      IdofsLoc = (/1,5,2,6,3,7,4,8,9/)

    case default
      call output_line('Unsupported element type!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'dofprep_initDofSet')
      call sys_halt()
    end select

    ! Allocate memory for the local DOF numbers per element
    allocate(rdofSubset%p_IdofsLoc(rdofSubset%ndofsPerElement,&
                                   rdofSubset%nelements))

    ! Loop over all degrees of freedoms of all elements and store
    ! the local number of the degrees of freedom per element
    do iel = 1, rdofSubset%nelements
      rdofSubset%p_IdofsLoc(:,iel) = IdofsLoc
    end do
    
    ! Deallocate temporal memory
    deallocate(IdofsLoc)

    ! Initialise coordinates of the DOFs?
    if (associated(rdomainIntSubset%p_Dcoords)) then
      
      ! Allocate memory for the coordinates of the DOFs
      if (rdofSubset%ndofsPerElement .gt. 0) then
        allocate(rdofSubset%p_DdofCoords(size(rdomainIntSubset%p_Dcoords,1),&
                 rdofSubset%ndofsPerElement,rdofSubset%nelements))
      else
        nullify(rdofSubset%p_DdofCoords)
      end if

      ! Calculate the coordinates of the DOFs.
      ! What type of finite element are we?
      select case(elem_getPrimaryElement(rdofSubset%celement))
      case (EL_P0)
        ! For P0 finite elements the degrees of freedom coincide with
        ! the centroid of the triangle, thus they need to be computed
        ! from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) =&
              (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel))/3.0_DP
        end do

      case (EL_Q0)
        ! For Q0 finite elements the degrees of freedom coincide with
        ! the center of the quadrilateral, thus they need to be computed
        ! from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) = 0.25_DP *&
              (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel))
        end do

      case (EL_P1,EL_Q1)
        ! For P1 and Q1 finite elements the degrees of freedom coincide
        ! with the vertices of the element, thus they can be copied.
        do iel = 1, rdofSubset%nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            rdofSubset%p_DdofCoords(:,idofe,iel) =&
                rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(idofe,iel),iel)
          end do
        end do

      case (EL_P1T)
        ! For P1~ finite elements the degrees of freedom coincide with the edge
        ! midpoints, thus they need to be computed from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel))
          rdofSubset%p_DdofCoords(:,2,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel))
          rdofSubset%p_DdofCoords(:,3,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel))
        end do

      case (EL_Q1T)
        ! For Q1~ finite elements the degrees of freedom coincide with the edge
        ! midpoints, thus they need to be computed from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel))
          rdofSubset%p_DdofCoords(:,2,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel))
          rdofSubset%p_DdofCoords(:,3,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel))
          rdofSubset%p_DdofCoords(:,4,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel))
        end do

      case (EL_P2)
        ! For P2 finite elements three degrees of freedom coincide
        ! with the vertices of the elements and three degrees of
        ! freedom coincide with the edge midpoints.
        do iel = 1, rdofSubset%nelements
          ! Copy the vertex coordinates
          rdofSubset%p_DdofCoords(:,1,iel) =&
              rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)
          rdofSubset%p_DdofCoords(:,3,iel) =&
              rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)
          rdofSubset%p_DdofCoords(:,5,iel) =&
              rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)
          
          ! Compute the coordinates of the edge midpoints
          rdofSubset%p_DdofCoords(:,2,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel))
          rdofSubset%p_DdofCoords(:,4,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel))
          rdofSubset%p_DdofCoords(:,6,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel))
        end do

      case (EL_Q2)
        ! For Q2 finite elements four degrees of freedom coincide with
        ! the vertices of the elements and four degrees of freedom
        ! coincide with the edge midpoints. One additional degree of
        ! freedom is located at the center of the element.
        do iel = 1, rdofSubset%nelements
          ! Copy the vertex coordinates
          rdofSubset%p_DdofCoords(:,1,iel) =&
              rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)
          rdofSubset%p_DdofCoords(:,3,iel) =&
              rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)
          rdofSubset%p_DdofCoords(:,5,iel) =&
              rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)
          rdofSubset%p_DdofCoords(:,7,iel) =&
              rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel)

          ! Compute the coordinates of the edge midpoints
          rdofSubset%p_DdofCoords(:,2,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel))
          rdofSubset%p_DdofCoords(:,4,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel))
          rdofSubset%p_DdofCoords(:,6,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel))
          rdofSubset%p_DdofCoords(:,8,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel))
        
          ! Compute the coordinates of the element center
          rdofSubset%p_DdofCoords(:,9,iel) =&
              0.25_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel))
        end do
        
      case default
        call output_line('Unsupported element type!',&
            OU_CLASS_ERROR, OU_MODE_STD, 'dofprep_initDofSet')
        call sys_halt()
      end select
      
    else
      nullify(rdofSubset%p_DdofCoords)
    end if

  end subroutine dofprep_initDofSet

  !*****************************************************************************

!<subroutine>

  subroutine dofprep_initDofSetAtBoundary(rdofSubset, rdomainIntSubset,&
      IelementOrientation)

!<description>
    ! This routine initialises a t_dofSubset structure using the
    ! information provided by the domain integration subset
    ! rdomainIntSubset at the boundary.
!</description>

!<input>
    ! A domain integration subset
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! OPTIONAL: An array providing the orientation of the element at
    ! the boundary. If not given, then rdomainIntSubset must provide
    ! this information. Otherwise, an error is thrown.
    integer, dimension(:), target, optional :: IelementOrientation
!</input>

!<output>
    ! A DOF subset structure, initialised according to the parameters
    type(t_dofSubset), intent(out) :: rdofSubset
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:,:), allocatable :: IdofsLoc
    integer, dimension(:), pointer :: p_IelementOrientation
    integer :: iel,idofe
    
    
    ! Get orientation of elements
    if (present(IelementOrientation)) then
      p_IelementOrientation => IelementOrientation
    elseif (associated(rdomainIntSubset%p_IelementOrientation)) then
      p_IelementOrientation => rdomainIntSubset%p_IelementOrientation
    else
      call output_line('Element orientation is not provided!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'dofprep_initDofSetAtBoundary')
      call sys_halt()
    end if

    ! Initialise constants in the structure
    rdofSubset%nelements = rdomainIntSubset%nelements
    rdofSubset%celement  = rdomainIntSubset%celement

    ! What type of finite element are we?
    select case(elem_getPrimaryElement(rdofSubset%celement))
    case (EL_P0,EL_Q0)
      rdofSubset%ndofsPerElement = 0

    case (EL_P1T)
      rdofSubset%ndofsPerElement = 1
      allocate(IdofsLoc(1,3))
      IdofsLoc = reshape((/1, 2, 3/),shape(IdofsLoc))

    case (EL_Q1T)
      rdofSubset%ndofsPerElement = 1
      allocate(IdofsLoc(1,4))
      IdofsLoc = reshape((/1, 2, 3, 4/),shape(IdofsLoc))
      
    case (EL_P1)
      rdofSubset%ndofsPerElement = 2
      allocate(IdofsLoc(2,3))
      IdofsLoc = reshape((/1,2, 2,3, 3,1/),shape(IdofsLoc))
      
    case (EL_Q1)
      rdofSubset%ndofsPerElement = 2
      allocate(IdofsLoc(2,4))
      IdofsLoc = reshape((/1,2, 2,3, 3,4, 4,1/),shape(IdofsLoc))
      
    case (EL_P2)
      rdofSubset%ndofsPerElement = 3
      allocate(IdofsLoc(3,3))
      IdofsLoc = reshape((/1,4,2, 2,5,3, 3,6,1/),shape(IdofsLoc))

    case (EL_Q2)
      rdofSubset%ndofsPerElement = 3
      allocate(IdofsLoc(3,4))
      IdofsLoc = reshape((/1,5,2, 2,6,3, 3,7,4, 4,8,1/),shape(IdofsLoc))

    case default
      call output_line('Unsupported element type!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'dofprep_initDofSetAtBoundary')
      call sys_halt()
    end select
    
    ! Allocate memory for the local DOF numbers per element
    allocate(rdofSubset%p_IdofsLoc(rdofSubset%ndofsPerElement,&
                                   rdofSubset%nelements))

    ! Loop over all degrees of freedoms of all elements and store
    ! the local number of the degrees of freedom per element
    do iel = 1, rdofSubset%nelements
      do idofe = 1, rdofSubset%ndofsPerElement
        rdofSubset%p_IdofsLoc(idofe,iel) =&
            IdofsLoc(idofe,p_IelementOrientation(iel))
      end do
    end do
    
    ! Deallocate temporal memory
    deallocate(IdofsLoc)
    
    
    ! Initialise coordinates of the DOFs?
    if (associated(rdomainIntSubset%p_Dcoords)) then
      
      ! Allocate memory for the coordinates of the DOFs
      if (rdofSubset%ndofsPerElement .gt. 0) then
        allocate(rdofSubset%p_DdofCoords(size(rdomainIntSubset%p_Dcoords,1),&
                 rdofSubset%ndofsPerElement,rdofSubset%nelements))
      else
        nullify(rdofSubset%p_DdofCoords)
      end if
    
      ! Calculate the coordinates of the DOFs.
      ! What type of finite element are we?
      select case(elem_getPrimaryElement(rdofSubset%celement))       
      case (EL_P1,EL_Q1)
        ! For P1 and Q1 finite elements the degrees of freedom coincide
        ! with the vertices of the element, thus they can be copied.
        do iel = 1, rdofSubset%nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            rdofSubset%p_DdofCoords(:,idofe,iel) =&
                rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(idofe,iel),iel)
          end do
        end do

      case (EL_P1T,EL_Q1T)
        ! For P1~ and Q1~ finite elements the degrees of freedom coincide with the
        ! edge midpoints, thus they need to be computed from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel))
        end do
      
      case (EL_P2,EL_Q2)
        ! For P2 and Q2 finite elements some degrees of freedom coincide
        ! with the vertices of the element which can be copied and others
        ! coincide with the edge midpoints which must be computed.
        do iel = 1, rdofSubset%nelements
          do idofe = 1, rdofSubset%ndofsPerElement, 2
            rdofSubset%p_DdofCoords(:,idofe,iel) =&
                rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(idofe,iel),iel)
          end do
          
          ! Compute the coordinate of the edge midpoint
          rdofSubset%p_DdofCoords(:,2,iel) =&
              0.5_DP * (rdofSubset%p_DdofCoords(:,1,iel)+&
                        rdofSubset%p_DdofCoords(:,3,iel))
        end do
        
      case default
        call output_line('Unsupported element type!',&
            OU_CLASS_ERROR, OU_MODE_STD, 'dofprep_initDofSetAtBoundary')
        call sys_halt()
      end select
      
    else
      nullify(rdofSubset%p_DdofCoords)
    end if
    
  end subroutine dofprep_initDofSetAtBoundary

  !*****************************************************************************

!<subroutine>

  subroutine dofprep_doneDofSet(rdofSubset)

!<description>
    ! This routine releases memory allocated in dofprep_initDofSetXXX
!</description>

!<inputoutput>
    ! A DIF subset structure
    type(t_dofSubset), intent(inout) :: rdofSubset
!</inputoutput>
!</subroutine>

    ! Clear structure
    rdofSubset%nelements       = 0
    rdofSubset%ndofsPerElement = 0

    if (associated(rdofSubset%p_IdofsLoc)) then
      deallocate(rdofSubset%p_IdofsLoc)
      nullify(rdofSubset%p_IdofsLoc)
    end if

    if (associated(rdofSubset%p_DdofCoords)) then
      deallocate(rdofSubset%p_DdofCoords)
      nullify(rdofSubset%p_DdofCoords)
    end if
    
  end subroutine dofprep_doneDofSet

end module dofpreprocessing
