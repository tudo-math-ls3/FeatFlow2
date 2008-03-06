  ! This file describes the interface for a callback routine that is called
  ! for calculating values at the boundary for boundary conditions.
  ! In contrast to the fgetBoundaryValues interface defined in intf_bcassembly.inc,
  ! this interface is used for calculating values on a boundary which is described
  ! using a mesh region structure instead of an analytic boundary region structure.
  
  INTERFACE
  
  !<subroutine>

    SUBROUTINE fgetBoundaryValuesMR (Icomponents,rdiscretisation,rmeshRegion,&
                                     cinfoNeeded,Dcoords,Dvalues,rcollection)
    
    USE collection
    USE spatialdiscretisation
    USE meshregion
    
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
    INTEGER, DIMENSION(:), INTENT(IN)                           :: Icomponents

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! Mesh region that is currently being processed.
    TYPE(t_meshRegion), INTENT(IN)                              :: rmeshRegion

    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    INTEGER, INTENT(IN)                                         :: cinfoNeeded
    
    ! The coordinates of the point for which the boundary values are to be
    ! calculated.
    REAL(DP), DIMENSION(:), INTENT(IN)                          :: Dcoords

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), INTENT(IN), OPTIONAL                    :: rcollection

  !</input>
  
  !<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1). 
    ! If multiple values are needed, they are collected here (e.g. for 
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    !
    ! The function may return SYS_INFINITY as a value. This indicates the
    ! framework to ignore the node and treat it as 'natural boundary condition'
    ! node.
    REAL(DP), DIMENSION(:), INTENT(OUT)                         :: Dvalues
  !</output>
    
  !</subroutine>
  
    END SUBROUTINE
    
  END INTERFACE
  