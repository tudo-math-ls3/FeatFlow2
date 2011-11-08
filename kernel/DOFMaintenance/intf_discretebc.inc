  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a callback routine that is called
  ! for calculating values at the boundary for boundary conditions.
  ! In contrast to the fgetBoundaryValues interface defined in intf_bcassembly.inc,
  ! this interface is used for calculating values on a boundary which is described
  ! using a mesh region structure instead of an analytic boundary region structure.
  
  interface
  
  !<subroutine>

    subroutine fgetBoundaryValuesMR (Icomponents,rdiscretisation,rmeshRegion,&
                                     cinfoNeeded,Iwhere,Dwhere,Dcoords,Dvalues,&
                                     rcollection)
    
    use fsystem
    use collection
    use spatialdiscretisation
    use meshregion
    
  !<description>
    ! This subroutine is called during the assembly of boundary conditions which
    ! are defined on mesh regions.
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
    integer, dimension(:), intent(in) :: Icomponents

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! Mesh region that is currently being processed.
    type(t_meshRegion), intent(in) :: rmeshRegion

    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    integer, intent(in) :: cinfoNeeded
    
    ! An array holding information about what type of DOF is currently processed.
    ! The information is build up as follows:
    ! Iwhere(1) = vertice number of the DOF, if the DOF is vertice-based, otherwise 0
    ! Iwhere(2) = edge number of the DOF, if the DOF is edge-based, otherwise 0
    ! Iwhere(3) = face number of the DOF, if the DOF is face-based, otherwise 0
    ! Iwhere(4) = currently processed element number.
    ! If Iwhere(1) = Iwhere(2) = Iwhere(3) = 0, then the DOF is element based.
    integer, dimension(4), intent(in) :: Iwhere
    
    ! The coordinates of the point which is currently processed, given in
    ! reference coordinates of the currently processed cell type (edge,face,element).
    ! If the DOF is vertice-based, then Dwhere is undefined.
    ! If the DOF is edge-based or element-based in 1D, then Dwhere has dimension 1.
    ! If the DOF is face-based or element-based in 2D, then Dwhere has dimension 2.
    ! IF the DOF is element-based in 3D, then Dwhere has dimension 3.
    real(DP), dimension(:), intent(in) :: Dwhere
    
    ! The coordinates of the point for which the boundary values are to be
    ! calculated.
    real(DP), dimension(:), intent(in) :: Dcoords
  !</input>

  !<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine. May point to NULL() if not defined.
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
  
    end subroutine
    
  end interface
  
