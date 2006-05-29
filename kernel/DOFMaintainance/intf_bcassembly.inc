  ! This file describes the interface for a callback routine that is called
  ! for calculating values at the boundary for boundary conditions.
  
  INTERFACE
  
  !<subroutine>

    SUBROUTINE fgetBoundaryValues (rdiscretisation,rbcRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, p_rcollection, Dvalues)
    
    USE collection
    USE spatialdiscretisation
    USE discretebc
    
  !<description>
    ! This subroutine is called during the discretisation of boundary
    ! conditions. It calculates a special quantity on the boundary, which is
    ! then used by the discretisation routines to generate a discrete
    ! 'snapshot' of the (actually analytic) boundary conditions.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! Boundary condition region that is currently being processed.
    ! (This e.g. defines the type of boundary conditions that are
    !  currently being calculated, as well as information about the current
    !  boundary segment 'where we are at the moment'.)
    TYPE(t_bcRegion), INTENT(IN)                                :: rbcRegion
    
    
    ! The element number on the boundary which is currently being processed
    INTEGER(I32), INTENT(IN)                                    :: ielement
    
    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    INTEGER, INTENT(IN)                                         :: cinfoNeeded
    
    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   iwhere = number of the edge where the value integral mean value
    !            should be computed
    INTEGER, INTENT(IN)                                         :: iwhere

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   dwhere = 0 (not used)
    REAL(DP), INTENT(IN)                                        :: dwhere
     
    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                  :: p_rcollection

  !</input>
  
  !<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1). 
    ! If multiple values are needed, they are collected here (e.g. for 
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    REAL(DP), DIMENSION(:), INTENT(OUT)                         :: Dvalues
  !</output>
    
  !</subroutine>
  
    END SUBROUTINE
    
  END INTERFACE
  