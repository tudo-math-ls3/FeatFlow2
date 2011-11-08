  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a callback routine that is called
  ! for calculating values at the boundary for boundary conditions of fictitious
  ! boundary components.
  
  interface
  
  !<subroutine>

    subroutine fgetBoundaryValuesFBC (Icomponents,rdiscretisation,&
                                      Revaluation, rcollection)
    
    use fsystem
    use collection
    use spatialdiscretisation
    use discretefbc
    
  !<description>
    ! This subroutine is called during the discretisation of boundary
    ! conditions on fictitious boundary components. It calculates a special quantity
    ! on the boundary, which is then used by the discretisation routines to
    ! generate a discrete 'snapshot' of the (actually analytic) boundary conditions.
    !
    ! The routine must calculate the values on all elements of the element
    ! list Ielements simultaneously. Iwhere is a list with vertex or edge numbers
    ! where information is to be retrieved. Dvalues is filled with function values
    ! while Binside is set to TRUE for every vertex/edge that is inside of the
    ! corresponding fictitious boundary region (identified by rbcRegion).
  !</description>
    
  !<input>
    ! Component specifier.
    ! For Dirichlet boundary:
    !   Icomponents(1..SIZE(Icomponents)) defines the number of the solution component,
    !   the value should be calculated for
    !   (e.g. 1=1st solution component, e.g. X-velocity,
    !         2=2nd solution component, e.g. Y-velocity,...,
    !         3=3rd solution component, e.g. pressure)
    !   Example: Icomponents(:) = [1,2] -> Compute velues for X- and Y-velocity
    !     (1=x, 2=y component)
    integer, dimension(:), intent(in) :: Icomponents
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_blockDiscretisation), intent(in) :: rdiscretisation
    
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection

  !</input>
  
  !<inputoutput>
    ! A t_discreteFBCevaluation structure array that defines what to evaluate,
    ! where to evaluate and which accepts the return values.
    ! This callback routine must check out the cinfoNeeded-entry in this structure
    ! to find out what to evaluate.
    ! The other entries in this structure describe where to evaluate.
    ! The result of the evaluation must be written into the p_Dvalues array entry
    ! in this structure.
    !
    ! The number of structures in this array depend on what to evaluate:
    !
    ! For Dirichlet boundary:
    !   revaluation contains as many entries as Icomponents; every entry in
    !   Icomponent corresponds to one entry in revaluation
    !   (so Icomponent(1)=1 defines to evaluate the X-velocity while the
    !    values for the X-velocity are written to revaluation(1)\%p_Dvalues;
    !    Icomponent(2)=2 defines to evaluate the Y-velocity while the values
    !    for the Y-velocity are written to revaluation(2)\%p_Dvalues, etc).
    !
    type(t_discreteFBCevaluation), dimension(:), intent(inout) :: Revaluation
  !</inputoutput>
    
  !</subroutine>
  
    end subroutine
    
  end interface
  
