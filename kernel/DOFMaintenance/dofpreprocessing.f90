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

!$use omp_lib
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

    ! For degrees of freedom located at the boundary, p_DedgePosition
    ! holds the parameter values of the degrees of freedom.
    ! Typically, this pointer is not associated since the position of
    ! degrees of freedom is not needed. At the boundary in 2D, it may
    ! be useful to know where DOFs are located at the boundary.
    real(DP), dimension(:,:), pointer :: p_DdofPosition => null()
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
    case (EL_P0_1D,EL_P0_2D,EL_Q0_2D,EL_P0_3D,EL_Q0_3D)
      IdofsLoc = (/1/)

      ! =-=-= 1D =-=-=

    case (EL_P1_1D)
      IdofsLoc = (/1,2/)

    case (EL_P2_1D)
      IdofsLoc = (/1,3,2/)

      ! =-=-= 2D =-=-=

    case (EL_P1_2D,EL_P1T_2D)
      IdofsLoc = (/1,2,3/)

    case (EL_Q1_2D,EL_Q1T_2D)
      IdofsLoc = (/1,2,3,4/)

    case (EL_P2_2D)
      IdofsLoc = (/1,4,2,5,3,6/)

    case (EL_Q2_2D)
      IdofsLoc = (/1,5,2,6,3,7,4,8,9/)

      ! =-=-= 3D =-=-=

    case (EL_P1_3D)
      IdofsLoc = (/1,2,3,4/)

    case (EL_Q1_3D)
      IdofsLoc = (/1,2,3,4,5,6,7,8/)

    case (EL_Q1T_3D)
      IdofsLoc = (/1,2,3,4,5,6/)

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
    if (associated(p_IelementOrientation)) then
      do iel = 1, rdofSubset%nelements
        do idofe = 1, rdofSubset%ndofsPerElement
          rdofSubset%p_IdofsLoc(:,iel) =&
              IdofsLoc(mod(p_IelementOrientation(iel)+idofe-2,&
                           rdofSubset%ndofsPerElement)+1)
        end do
      end do
    else
      do iel = 1, rdofSubset%nelements
        rdofSubset%p_IdofsLoc(:,iel) = IdofsLoc
      end do
    end if

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
      case (EL_P0_1D)
        ! For P0 finite elements in 1D the degrees of freedom coincide
        ! with the midpoint of the line, thus they need to be computed
        ! from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel))
        end do

      case (EL_P0_2D)
        ! For P0 finite elements in 2D the degrees of freedom coincide
        ! with the centroid of the triangle, thus they need to be
        ! computed from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) =&
              (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel))/3.0_DP
        end do

      case (EL_P0_3D)
        ! For P0 finite elements in 3D the degrees of freedom coincide
        ! with the centroid of the tetrahedra, thus they need to be
        ! computed from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) =&
              0.25_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel))
        end do

      case (EL_Q0_2D)
        ! For Q0 finite elements in 2D the degrees of freedom coincide
        ! with the center of the quadrilateral, thus they need to be
        ! computed from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) = 0.25_DP *&
              (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel))
        end do

      case (EL_Q0_3D)
        ! For Q0 finite elements in 3D the degrees of freedom coincide
        ! with the center of the hexahedra, thus they need to be
        ! computed from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) = 0.125_DP *&
              (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(5,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(6,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(7,iel),iel)+&
               rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(8,iel),iel))
        end do


      case (EL_P1_1D,EL_P1_2D,EL_P1_3D,EL_Q1_2D,EL_Q1_3D)
        ! For P1 and Q1 finite elements the degrees of freedom coincide
        ! with the vertices of the element, thus they can be copied.
        do iel = 1, rdofSubset%nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            rdofSubset%p_DdofCoords(:,idofe,iel) =&
                rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(idofe,iel),iel)
          end do
        end do

      case (EL_P1T_2D)
        ! For P1~ finite elements in 2D the degrees of freedom
        ! coincide with the edge midpoints, thus they need to be
        ! computed from the vertex coordinates
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

      case (EL_Q1T_2D)
        ! For Q1~ finite elements in 2D the degrees of freedom
        ! coincide with the edge midpoints, thus they need to be
        ! computed from the vertex coordinates
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

      case (EL_Q1T_3D)
        ! For Q1~ finite elements in 3D the degrees of freedom
        ! coincide with the face midpoints, thus they need to be
        ! computed from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) =&
              0.25_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel))
          rdofSubset%p_DdofCoords(:,2,iel) =&
              0.25_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(5,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(6,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel))
          rdofSubset%p_DdofCoords(:,3,iel) =&
              0.25_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(6,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(7,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel))
          rdofSubset%p_DdofCoords(:,4,iel) =&
              0.25_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(7,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(8,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel))
          rdofSubset%p_DdofCoords(:,5,iel) =&
              0.25_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(8,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(5,iel),iel))
          rdofSubset%p_DdofCoords(:,6,iel) =&
              0.25_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(5,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(8,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(7,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(6,iel),iel))
        end do

      case (EL_P2_2D)
        ! For P2 finite elements in 2D three degrees of freedom
        ! coincide with the vertices of the elements and three degrees
        ! of freedom coincide with the edge midpoints.
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

      case (EL_Q2_2D)
        ! For Q2 finite elements in 2D four degrees of freedom
        ! coincide with the vertices of the elements and four degrees
        ! of freedom coincide with the edge midpoints. One additional
        ! degree of freedom is located at the center of the element.
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
    case (EL_P0_1D,EL_P0_2D,EL_P0_3D,EL_Q0_2D,EL_Q0_3D,EL_Y0_3D,EL_R0_3D)
      rdofSubset%ndofsPerElement = 0

      ! =-=-= 1D =-=-=

    case (EL_P1_1D,EL_P2_1D)
      rdofSubset%ndofsPerElement = 1
      allocate(IdofsLoc(1,2))
      IdofsLoc = reshape((/1, 2/),shape(IdofsLoc))

      ! =-=-= 2D =-=-=

    case (EL_P1T_2D)
      rdofSubset%ndofsPerElement = 1
      allocate(IdofsLoc(1,3))
      IdofsLoc = reshape((/1, 2, 3/),shape(IdofsLoc))

    case (EL_Q1T_2D)
      rdofSubset%ndofsPerElement = 1
      allocate(IdofsLoc(1,4))
      IdofsLoc = reshape((/1, 2, 3, 4/),shape(IdofsLoc))

    case (EL_P1_2D)
      rdofSubset%ndofsPerElement = 2
      allocate(IdofsLoc(2,3))
      IdofsLoc = reshape((/1,2, 2,3, 3,1/),shape(IdofsLoc))

    case (EL_Q1_2D)
      rdofSubset%ndofsPerElement = 2
      allocate(IdofsLoc(2,4))
      IdofsLoc = reshape((/1,2, 2,3, 3,4, 4,1/),shape(IdofsLoc))

    case (EL_P2_2D)
      rdofSubset%ndofsPerElement = 3
      allocate(IdofsLoc(3,3))
      IdofsLoc = reshape((/1,4,2, 2,5,3, 3,6,1/),shape(IdofsLoc))

    case (EL_Q2_2D)
      rdofSubset%ndofsPerElement = 3
      allocate(IdofsLoc(3,4))
      IdofsLoc = reshape((/1,5,2, 2,6,3, 3,7,4, 4,8,1/),shape(IdofsLoc))

      ! =-=-= 3D =-=-=

    case (EL_P1_3D)
      rdofSubset%ndofsPerElement = 3
      allocate(IdofsLoc(3,4))
      IdofsLoc = reshape((/1,2,3, 1,4,2, 2,4,3, 3,4,1/),shape(IdofsLoc))

    case (EL_Q1_3D)
      rdofSubset%ndofsPerElement = 4
      allocate(IdofsLoc(4,6))
      IdofsLoc = reshape((/1,2,3,4, 1,5,6,2, 2,6,7,3,&
                           3,7,8,4, 1,4,8,5, 5,8,7,6/),shape(IdofsLoc))

    case (EL_Y1_3D)
      rdofSubset%ndofsPerElement = 4
      allocate(IdofsLoc(4,5))
      IdofsLoc = reshape((/1,2,3,4, 1,5,2,0, 2,5,3,0,&
                           3,5,4,0, 4,5,1,0/),shape(IdofsLoc))

    case (EL_R1_3D)
      rdofSubset%ndofsPerElement = 4
      allocate(IdofsLoc(4,5))
      IdofsLoc = reshape((/1,2,3,0, 1,4,5,2, 2,5,6,3,&
                           3,6,4,1, 4,6,5,0/),shape(IdofsLoc))

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
      case (EL_P1_1D,EL_P1_2D,EL_Q1_2D,EL_P1_3D,EL_Q1_3D,EL_Y1_3D,EL_R1_3D)
        ! For P1,Q1,Y1 and R1 finite elements the degrees of freedom coincide
        ! with the vertices of the element, thus they can be copied.
        do iel = 1, rdofSubset%nelements
          do idofe = 1, rdofSubset%ndofsPerElement
            rdofSubset%p_DdofCoords(:,idofe,iel) =&
                rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(idofe,iel),iel)
          end do
        end do

      case (EL_P1T_2D,EL_Q1T_2D)
        ! For P1~ and Q1~ finite elements in 2D the degrees of freedom
        ! coincide with the edge midpoints, thus they need to be
        ! computed from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) =&
              0.5_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                        rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel))
        end do

      case (EL_P2_2D,EL_Q2_2D)
        ! For P2 and Q2 finite elements in 2D some degrees of freedom
        ! coincide with the vertices of the element which can be
        ! copied and others coincide with the edge midpoints which
        ! must be computed.
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

      case (EL_Q1T_3D)
        ! For Q1~ finite elements in 3D the degrees of freedom
        ! coincide with the face centers, thus they need to be
        ! computed from the vertex coordinates
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofCoords(:,1,iel) =&
              0.25_DP * (rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(1,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(2,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(3,iel),iel)+&
                         rdomainIntSubset%p_Dcoords(:,rdofSubset%p_IdofsLoc(4,iel),iel))
        end do

      case default
        call output_line('Unsupported element type!',&
            OU_CLASS_ERROR, OU_MODE_STD, 'dofprep_initDofSetAtBoundary')
        call sys_halt()
      end select

    else
      nullify(rdofSubset%p_DdofCoords)
    end if


    ! Initialise positions of the DOFs?
    if (associated(rdomainIntSubset%p_DedgePosition)) then

      ! Allocate memory for the positions of the DOFs
      if (rdofSubset%ndofsPerElement .gt. 0) then
        allocate(rdofSubset%p_DdofPosition(rdofSubset%ndofsPerElement,&
                                           rdofSubset%nelements))
      else
        nullify(rdofSubset%p_DdofPosition)
      end if

      ! Calculate the positions of the DOFs.
      ! What type of finite element are we?
      select case(elem_getPrimaryElement(rdofSubset%celement))
      case (EL_P1_2D,EL_Q1_2D)
        ! For P1 and Q1 finite elements in 2D the degrees of freedom
        ! coincide with the vertices of the element, thus they are
        ! located at the endpoints of the reference interval [-1,1]
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofPosition(:,iel) = rdomainIntSubset%p_DedgePosition(:,iel)
        end do

      case (EL_P1T_2D,EL_Q1T_2D)
        ! For P1~ and Q1~ finite elements in 2D the degrees of freedom
        ! coincide with the edge midpoints, thus they are located at
        ! the midpoint of the reference interval [-1,1]
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofPosition(1,iel) =&
              0.5_DP * (rdomainIntSubset%p_DedgePosition(1,iel)+&
                        rdomainIntSubset%p_DedgePosition(2,iel))
        end do

      case (EL_P2_2D,EL_Q2_2D)
        ! For P2 and Q2 finite elements in 2D two degrees of freedom
        ! coincide with the vertices of the element and one degree of
        ! freedom coincides with the edge midpoints. Thus, they are
        ! located at the two endpoints and the midpoint of the
        ! reference interval [-1,1].
        do iel = 1, rdofSubset%nelements
          rdofSubset%p_DdofPosition(1,iel) = rdomainIntSubset%p_DedgePosition(1,iel)
          rdofSubset%p_DdofPosition(2,iel) =&
              0.5_DP * (rdomainIntSubset%p_DedgePosition(1,iel)+&
                        rdomainIntSubset%p_DedgePosition(2,iel))
          rdofSubset%p_DdofPosition(3,iel) = rdomainIntSubset%p_DedgePosition(2,iel)
        end do

      case default
        call output_line('Unsupported element type!',&
            OU_CLASS_ERROR, OU_MODE_STD, 'dofprep_initDofSetAtBoundary')
        call sys_halt()
      end select

    else
      nullify(rdofSubset%p_DdofPosition)
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

    if (associated(rdofSubset%p_DdofPosition)) then
      deallocate(rdofSubset%p_DdofPosition)
      nullify(rdofSubset%p_DdofPosition)
    end if

  end subroutine dofprep_doneDofSet

end module dofpreprocessing
