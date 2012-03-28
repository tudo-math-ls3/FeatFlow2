  ! This file describes the interface for a callback routine that is called
  ! for calculating values at the boundary for boundary conditions.

  interface

  !<subroutine>

    subroutine fgetBoundaryValues (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)

    use fsystem
    use boundary
    use collection
    use spatialdiscretisation
    use discretebc

  !<description>
    ! This subroutine is called during the discretisation of boundary
    ! conditions. It calculates a special quantity on the boundary, which is
    ! then used by the discretisation routines to generate a discrete
    ! 'snapshot' of the (actually analytic) boundary conditions.
  !</description>

  !<input>
    ! Component specifier.
    ! For Dirichlet boundary:
    !   Icomponents(1) defines the number of the solution component, the value
    !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry,
    !   2=2nd solution component, e.g. Y-velocity,...,
    !   3=3rd solution component, e.g. pressure)
    ! For pressure drop boundary / normal stress:
    !   Velocity components that are affected by the normal stress
    !   (usually "1 2" for x- and y-velocity while returned value musr specify
    !   the pressure at the boundary)
    integer, dimension(:), intent(in)                           :: Icomponents

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation

    ! Boundary region that is currently being processed.
    type(t_boundaryRegion), intent(in)                          :: rboundaryRegion

    ! The element number on the boundary which is currently being processed
    integer, intent(in)                                         :: ielement

    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    integer, intent(in)                                         :: cinfoNeeded

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC :
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDFUNCMID :
    !   iwhere = number of the edge in which midpoint the value
    !            should be computed
    ! cinfoNeeded=DISCBC_NEEDDERIV :
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN :
    !   iwhere = number of the edge where the value integral mean value
    !            should be computed
    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS :
    !   iwhere = Number of the edge where the normal stress should be computed.
    integer, intent(in) :: iwhere

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC :
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDDERIV :
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN :
    !   dwhere = 0 (not used)
    ! cinfoNeeded=DISCBC_NEEDNORMALSTRESS :
    !   dwhere = parameter value of the point on edge iwhere where the normal
    !            stress should be computed.
    real(DP), intent(in) :: dwhere
  !</input>

  !<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
  !</inputoutput>

  !<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1).
    ! If multiple values are needed, they are collected here (e.g. for
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    !
    ! The function may return SYS_INFINITY as a value. This indicates the
    ! framework to ignore the node and treat it as 'natural boundary condition'
    ! node.
    real(DP), dimension(:), intent(out) :: Dvalues
  !</output>

  !</subroutine>

      ! To get the X/Y-coordinates of the boundary point, use:
      !
      ! real(DP) :: dx,dy
      !
      ! call boundary_getCoords(rdiscretisation%p_rboundary, &
      !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    end subroutine

  end interface

