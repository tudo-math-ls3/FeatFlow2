module ExtFEcomparer_discretisation

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  use fparser
  use paramlist

  use triangulation
  use meshgeneration
  use basicgeometry

  use element
  use spatialdiscretisation

  use ExtFEcomparer_typedefs

  implicit none

  private

  public :: ExtFEcomparer_init_discretisation
  public :: ExtFEcomparer_doneDiscretisation

contains



subroutine ExtFEcomparer_init_discretisation(rproblem)
!<description>
    ! This routine is the one to be called from outside,
    ! it will branch it way to the right subroutine.
    ! It is going to be a bit tricky since we have to cover
    ! all cases here: Dimension and type of vector, that is
    ! if it is a solution of a system (ie a (u,v,p) - vector
    ! from cc2d) or if it is just the written out x-velocity
    ! of flagship. To make life easier, we do only 1 branching in here
    ! and call subroutines which branch again. It can be done
    ! better in terms of performance, but we want to have an
    ! easy-to-read-code
    ! In this routine we branch with respect to the dimension
    ! in the routines ExtFE_init_discretisationxD we branch
    ! w.r.t. the element-setting
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

    select case(rproblem%iDimension)
        case(ExtFE_NDIM1)
            call ExtFE_init_discretisation_1D(rproblem)
        case(ExtFE_NDIM2)
            call ExtFE_init_discretisation_2D(rproblem)
        case(ExtFE_NDIM3)
            call ExtFE_init_discretisation_3D(rproblem)
        case default
            call output_line(&
            "The dimension of you problem is not implemented, dimension = " &
                //sys_siL(rproblem%iDimension,10), &
                OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_discretisation")
            call sys_halt()
        end select

end subroutine

subroutine ExtFE_init_discretisation_1D(rproblem)
!<description>
    ! This routine is containing calls to all 1D-inits
    ! at the moment there is only support for a discretisation
    ! using one elementtype, so this is for future purpose only
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>
    select case(rproblem%elementSetting)

        case(ExtFE_OneElement)

            call ExtFE_init_discretisation_1D_OneElementType(rproblem)

        case(ExtFE_ElementList)
            call ExtFE_init_discretisation_1D_ElementList(rproblem)

        case default
            call output_line(&
            "Your element setting is not supported, setting = " &
                //sys_siL(rproblem%elementSetting,10), &
                OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_discretisation_1D")
            call sys_halt()
        end select

end subroutine


subroutine ExtFE_init_discretisation_1D_OneElementType(rproblem)
!<description>
  ! This routine initialises a discretisation structure for a vector
  ! that contains 1D-Variables, and all of them use the same element
  ! type
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

    ! local variables
    integer :: NLMAX, NVAR,i
    type (t_triangulation), pointer :: p_triangulation => NULL()
    type (t_blockDiscretisation), pointer :: p_discretisation => NULL()

    ! We need to create a block-discretisation-structure
    ! so that we can store one variable in each block
    NLMAX = rproblem%NLMAX
    NVAR = rproblem%NVAR
    p_triangulation => rproblem%rtriangulation
    p_discretisation => rproblem%rdiscretisation

    ! Init a block discretisation - 1 block per variable
    call spdiscr_initBlockDiscr(rproblem%rdiscretisation, &
                NVAR, p_triangulation)

    ! Now init a discretisation for each block
    do i=1,NVAR
        call spdiscr_initDiscr_simple(p_discretisation%RspatialDiscr(i), &
              rproblem%ielemtype,p_triangulation)
    end do



end subroutine

subroutine ExtFE_init_discretisation_1D_ElementList(rproblem)
!<description>
  ! This routine initialises a discretisation structure for a vector
  ! that contains 1D-Variables, and all of them use the same element
  ! type
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

    ! local variables
    integer :: NLMAX, NVAR,i
    type (t_triangulation), pointer :: p_triangulation => NULL()
    type (t_blockDiscretisation), pointer :: p_discretisation => NULL()

    ! We need to create a block-discretisation-structure
    ! so that we can store one variable in each block
    NLMAX = rproblem%NLMAX
    NVAR = rproblem%NVAR
    p_triangulation => rproblem%rtriangulation
    p_discretisation => rproblem%rdiscretisation

    ! Init a block discretisation - 1 block per variable
    call spdiscr_initBlockDiscr(rproblem%rdiscretisation, &
                NVAR, p_triangulation)

    ! Now init a discretisation for each block
    do i=1,NVAR
        call spdiscr_initDiscr_simple(p_discretisation%RspatialDiscr(i), &
              rproblem%iElemList(i),p_triangulation)
    end do

end subroutine

subroutine ExtFE_init_discretisation_2D(rproblem)
!<description>
    ! This routine is containing calls to all 2D-inits, that is
    ! in particular if it is a (u,v,p) solution (i.e. from cc2d)
    ! or if it is containing a "postprocessed" vector,
    ! that is i.e. an output of flagship that solves a system
    ! in the conservative variables (momentum, ...) but writes out
    ! the speed
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

    select case (rproblem%elementSetting)

        case(ExtFE_OneElement)
            call ExtFE_init_Discretisation_2D_OneELementType(rproblem)

        case(ExtFE_ElementPair)
            ! Element pair - i.e. one element type for the speed,
            ! one for the pressure
            call ExtFE_init_Discretisation_2D_elementPair(rproblem)

        case(ExtFE_ElementList)
            call ExtFE_init_Discretisation_2D_ElementList(rproblem)

        case default
            call output_line(&
            "This combination of dimension and element setting is not &
            &yet implemented, your input was " &
                //sys_siL(rproblem%elementSetting,10), &
                OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_discretisation_2D")
            call sys_halt()
    end select
end subroutine


subroutine ExtFE_init_Discretisation_2D_OneELementType(rproblem)

!<inputoutput>
    ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

    ! local variables
    integer :: NLMAX, NVAR,i
    type (t_triangulation), pointer :: p_triangulation => NULL()
    type (t_blockDiscretisation), pointer :: p_discretisation => NULL()
    type(t_boundary), pointer :: p_boundary => NULL()

    ! We need to create a block-discretisation-structure
    ! so that we can store one variable in each block
    NLMAX = rproblem%NLMAX
    NVAR = rproblem%NVAR
    p_triangulation => rproblem%rtriangulation
    p_discretisation => rproblem%rdiscretisation
    p_boundary => rproblem%rboundary

    ! Init a block discretisation - 1 block per variable
    call spdiscr_initBlockDiscr(rproblem%rdiscretisation, &
                NVAR, p_triangulation)

    ! Now init a discretisation for each block
    do i=1,NVAR
        call spdiscr_initDiscr_simple(p_discretisation%RspatialDiscr(i), &
                rproblem%ielemType,p_triangulation,p_boundary)
    end do

end subroutine


 subroutine ExtFE_init_Discretisation_2D_elementPair(rproblem)

!<description>
  ! This routine initialises a discretisation structure for a vector
  ! that contains a 2D-Speed-Pressure-Pair with different elements
  ! for speed and pressure. The order of the components is
  ! 1. X-Velocity, 2. Y-Velocity, 3. Pressure
!</description>

!<inputoutput>
  ! A problem structure saving all information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>


!</subroutine>

    ! local variables
    integer :: iElemPairID, NLMAX
    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary => NULL()

    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation => NULL()

    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => NULL()

    ! Which discretisation is to use?
    iElemPairID = rproblem%iElemPairID

    ! Now set up discrezisation structures on all levels:

      NLMAX = rproblem%NLMAX

      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%rboundary
      p_rtriangulation => rproblem%rtriangulation
      p_rdiscretisation => rproblem%rdiscretisation

      ! Initialise the block discretisation according to the element specifier.
      call ExtFEcomparer_getDiscretisation_2D_elementPair(iElemPairID, &
               p_rdiscretisation,rproblem%rtriangulation, &
               rproblem%rboundary)

  end subroutine


subroutine ExtFE_init_Discretisation_2D_ElementList(rproblem)

!<inputoutput>
    ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

    ! local variables
    integer :: NLMAX, NVAR,i
    type (t_triangulation), pointer :: p_triangulation => NULL()
    type (t_blockDiscretisation), pointer :: p_discretisation => NULL()
    type(t_boundary), pointer :: p_boundary => NULL()

    ! We need to create a block-discretisation-structure
    ! so that we can store one variable in each block
    NLMAX = rproblem%NLMAX
    NVAR = rproblem%NVAR
    p_triangulation => rproblem%rtriangulation
    p_discretisation => rproblem%rdiscretisation
    p_boundary => rproblem%rboundary

    ! Init a block discretisation - 1 block per variable
    call spdiscr_initBlockDiscr(rproblem%rdiscretisation, &
                NVAR, p_triangulation)

    ! Now init a discretisation for each block
    do i=1,NVAR
        call spdiscr_initDiscr_simple(p_discretisation%RspatialDiscr(i), &
                rproblem%iElemList(i),p_triangulation,p_boundary)
    end do

end subroutine

subroutine ExtFE_init_discretisation_3D(rproblem)
!<description>
! This is just a spaceholder for the 3D-part that is not implemented
! yet.
!</description>

!<inputoutput>
    ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

    call output_line(&
            "3D is not yet implemented", &
                OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_discretisation_3D")
            call sys_halt()

end subroutine

! ***************************************************************************


  subroutine ExtFEcomparer_getDiscretisation_2D_elementPair(ielementType, &
      rdiscretisation, rtriangulation,rboundary,rsourceDiscretisation)

!<description>
  ! Initialises a discretisation structure according to an element combination
  ! identifier.
!</description>

!<input>
  ! Element combination identifier to use. The following identifiers are supported:
  !  0 = Q1(E031) / Q1~(E031) / Q0
  !  1 = Q1(E030) / Q1~(E030) / Q0
  !  2 = Q1(EM31) / Q1~(EM31) / Q0
  !  3 = Q1(EM30) / Q1~(EM30) / Q0 = standard
  !  4 = Q2 (E013) / Q2 (E013) / QP1
  !  5 = Q1(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
  !  6 = Q1(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
  !  7 = Q1(EM30) / Q1~(EM30) / Q0 new interface implementation
  !  8 = Q1(EB30) / Q1~(EB30) / Q0
  ! 10 = Q2(EB50) / Q2~(EB50) / QP1
  ! 11 = Q2(EM50) / Q2~(EM50) / QP1
  ! 12 = Q2 (E013) / Q2 (E013) / QP1 (nonparametric)
  ! 20 = Q1        / Q1        / Q1
  ! 30 = P1(E020) / P1~(E020) / P0
  integer, intent(in) :: ielementType

  ! Boundary structure for the discretisation
  type(t_boundary), intent(in), target, optional :: rboundary

  ! Triangulation structure for the discretisation
  type(t_triangulation), intent(in), target, optional :: rtriangulation

  ! OPTIONAL: Source discretisation. If specified, rdiscretisation is created
  ! as copy of rsourceDiscretisation with changed element type.
  ! If not specified, rboundary and rtriangulation must be present to
  ! specify mesh and boundary.
  type(t_blockDiscretisation), intent(in), target, optional :: rsourceDiscretisation
!</input>

!<output>
  ! Block discretisation structure that defines the discretisation according
  ! to ielementType.
  type(t_blockDiscretisation), intent(out), target :: rdiscretisation
!</output>

!</subroutine>

  ! Number of equations in our problem. velocity+velocity+pressure = 3
  integer, parameter :: nequations = 3

  ! local variables
  integer(I32) :: ieltypeUV, ieltypeP

    ! Initialise the element type identifiers according to ielementType
    select case (ielementType)
    case (0)
      ieltypeUV = EL_E031
      ieltypeP = EL_Q0

    case (1)
      ieltypeUV = EL_E030
      ieltypeP = EL_Q0

    case (2)
      ieltypeUV = EL_EM31
      ieltypeP = EL_Q0

    case (3)
      ieltypeUV = EL_EM30
      ieltypeP = EL_Q0

    case (4)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1

    case (5)
      ieltypeUV = EL_EM30_UNPIVOTED
      ieltypeP = EL_Q0

    case (6)
      ieltypeUV = EL_EM30_UNSCALED
      ieltypeP = EL_Q0

    case (7)
      ieltypeUV = EL_EM30_NEW
      ieltypeP = EL_Q0

    case (8)
      ieltypeUV = EL_EB30
      ieltypeP = EL_Q0

    case (10)
      ieltypeUV = EL_EB50
      ieltypeP = EL_QP1

    case (11)
      ieltypeUV = EL_EM50
      ieltypeP = EL_QP1

    case (12)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1NP

    case (13)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1NPD

    case (14)
      ieltypeUV = EL_Q1TB_2D
      ieltypeP = EL_QP1

    case (15)
      ieltypeUV = EL_Q1TBNP_2D
      ieltypeP = EL_QP1NPD

    case (20)
      ieltypeUV = EL_Q1
      ieltypeP = EL_Q1

    case (30)
      ieltypeUV = EL_P1T
      ieltypeP = EL_P0

    case (50)
      ieltypeUV = EL_EN50_2D
      ieltypeP  = EL_DCQP1_2D

    case (51)
      ieltypeUV = EL_EN51_2D
      ieltypeP  = EL_DCQP2_2D

    case (52)
      ieltypeUV = EL_QPW4P2_2D
      ieltypeP  = EL_QPW4DCP1_2D

    case default
      call output_line (&
          "Unknown discretisation: iElementType = "//sys_siL(ielementType,10), &
          OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_getDiscretisation_2D_speedPressure")
      call sys_halt()
    end select

    if (present(rsourceDiscretisation)) then

      ! Copy the source discretisation structure.
      call spdiscr_deriveBlockDiscr (rsourceDiscretisation,rdiscretisation)

      ! Replace the spatial discretisation structures inside by new ones
      ! with the correct element type.
      call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(1),  &
          ieltypeUV, rdiscretisation%RspatialDiscr(1))

      call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(2),  &
          ieltypeUV, rdiscretisation%RspatialDiscr(2))

      call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(3),  &
          ieltypeP, rdiscretisation%RspatialDiscr(3))

    else
      ! Ensure that rboundary/rtriangulation is present
      if ((.not. present(rboundary)) .or. &
          (.not. present(rtriangulation))) then
        call output_line (&
            "Boundary or triangulation not specified.", &
            OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_getDiscretisation_2D")
        call sys_halt()
      end if

      ! Check the mesh to ensure that the discretisation is compatible.
      select case(elem_igetShape(ieltypeUV))
      case (BGEOM_SHAPE_TRIA)
        ! Triangular element; ensure that we don"t have quads in the mesh
        if(rtriangulation%InelOfType(TRIA_NVEQUAD2D) .ne. 0) then
          call output_line (&
              "Discretisation does not support the current mesh!", &
              OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_getDiscretisation_2D")
          call sys_halt()
        end if

      case (BGEOM_SHAPE_QUAD)
        ! Quadrilateral element; ensure that we don"t have triangles in the mesh
        if(rtriangulation%InelOfType(TRIA_NVETRI2D) .ne. 0) then
          call output_line (&
              "Discretisation does not support the current mesh!", &
              OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_getDiscretisation_2D")
          call sys_halt()
        end if

      end select

      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies 3 blocks in the
      ! solution vector.
      call spdiscr_initBlockDiscr (&
          rdiscretisation,nequations,rtriangulation,rboundary)

      ! rdiscretisation%RspatialDiscr is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple ( &
          rdiscretisation%RspatialDiscr(1), &
          ieltypeUV,rtriangulation, rboundary)

      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards.
      call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
          rdiscretisation%RspatialDiscr(2))

      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1),  &
          ieltypeP, rdiscretisation%RspatialDiscr(3))

    end if

  end subroutine

!<subroutine>

  ! ***************************************************************************


  subroutine ExtFEcomparer_doneDiscretisation (rproblem)

    !<description>
    ! Releases the discretisation from the heap.
    !</description>

    !<inputoutput>
    ! A problem structure saving problem-dependent information.
    type(t_problem), intent(inout), target :: rproblem
    !</inputoutput>

      ! Remove the block discretisation structure and all substructures.
      call spdiscr_releaseBlockDiscr(rproblem%rdiscretisation)


end subroutine



end module
