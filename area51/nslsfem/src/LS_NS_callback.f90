!##############################################################################
!# ****************************************************************************
!# <name> nslsfem_callback </name>
!# ****************************************************************************
!#  
!# <purpose>
!# This module contains callback functions for the Navier stokes problem
!# that are used during the matrix/vector assembly for specifying analytical
!# data.
!# There are some callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in 'intf_xxxx.inc' files.
!#
!# 1.) getBoundaryValues_2D
!#     -> Returns analitical values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# </purpose>
!##############################################################################

module LS_NS_callback

  use fsystem
  use storage
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use spatialdiscretisation
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  use triangulation
  use derivatives
  use linearsystemscalar
  use linearsystemblock
  use element
  use paramlist
  use collection, only: t_collection  
  
  use genoutput  
  
  implicit none

contains

! ***************************************************************************

!<subroutine>
  subroutine getBoundaryValues_2D (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  
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
  !   Icomponents(1) defines the number of the boundary component, the value
  !   should be calculated for (e.g. 1=1st solution component, e.g. X-velocitry,
  !   2=2nd solution component, e.g. Y-velocity,...)
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
  ! cinfoNeeded=DISCBC_NEEDDERIV :
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN :
  !   iwhere = number of the edge where the value integral mean value
  !            should be computed
  integer, intent(in)                                          :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC :
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV :
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN :
  !   dwhere = 0 (not used)
  real(DP), intent(in)                                        :: dwhere
    
  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional                 :: rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1).
  ! If multiple values are needed, they are collected here (e.g. for
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>
  
!</subroutine>

    ! Some variables
    real(DP) :: dx,dy, dC, cu, pi, n
    
    ! Even more
    integer :: icomponent
    real(DP) :: y
    
    call boundary_getCoords(rdiscretisation%p_rboundary, &
         rboundaryRegion%iboundCompIdx, dwhere, dx, dy)
    
    ! Get from the current component of the PDE we are discretising:
    icomponent = Icomponents(1)
    ! -> 1=X-velocity, 2=Y-velocity, 3=Pressure, 4=Vorticity
        
    ! Default zero value
    Dvalues(1) = 0.0_DP

    select case (rcollection%IquickAccess(9))
    
    case (0)
      ! Reg. Cavity Flow
      ! Now, depending on the problem, calculate the actual velocity value.
      select case (icomponent)
      case (1) ! X-velocity
        if ((dwhere .ge. 2.0_DP) .and. (dwhere .le. 3.0_DP)) then
          Dvalues(1) = 16.0_DP*dx*dx * (1.0_DP - dx)*(1.0_DP - dx)
        end if

      case (2) ! Y-velocity
        Dvalues(1) = 0.0_DP

      case (3) ! Pressure
      
      end select
      
    case (1,2)
      ! FAC, zero stress outflow
      ! FAC, Dirichlet velocity outflow
      ! http://featflow.de/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark_re20.html
       select case (rboundaryregion%iboundcompidx)
       case (1)
          select case (icomponent)
          case (1) ! x-velocity
           ! Outflow
            if ((dwhere .ge. 1.0_DP) .and. (dwhere .le. 2.0_DP)) then
              dvalues(1) = (1.2_DP*dy * (0.41_DP - dy))/(0.41_DP*0.41_DP)
            end if
           ! Inflow  
            if ((dwhere .ge. 3.0_DP) .and. (dwhere .le. 4.0_DP)) then
              dvalues(1) = (1.2_DP*dy * (0.41_DP - dy))/(0.41_DP*0.41_DP)
            end if

          case (2) ! y-velocity
            ! nothing to do here.

          case (3) ! pressure
            ! nothing to do here.
          end select
       
       case (2)
          ! nothing to do here.
       end select
    
    case (3,4)
      ! Poiseuielle Flow, zeros stress outflow
      !    OR
      ! Poiseuielle Flow, Dirichlet velocity outflow
      ! Now, depending on the problem, calculate the actual velocity value.
      select case (icomponent)
      case (1) ! X-velocity
        if ((dwhere .ge. 3.0_DP) .and. (dwhere .le. 4.0_DP)) then
          y = 4.0_DP-dwhere
          Dvalues(1) = y*(1.0_DP-y)
        end if

        if ((dwhere .ge. 1.0_DP) .and. (dwhere .le. 2.0_DP)) then
          y = 2.0_DP-dwhere
          Dvalues(1) = y*(1.0_DP-y)
        end if

      case (2) ! Y-velocity
        ! Nothing to do here.

      case (3) ! Pressure
        if ((dwhere .ge. 1.0_DP) .and. (dwhere .le. 2.0_DP)) then
          Dvalues(1) = 0.0_DP
        end if
      end select

    
    case (5)
      ! Bolton, step flow
       select case (rboundaryRegion%IBOUNDCOMPIDX)
       case (1)
          select case (icomponent)
          case (1) ! X-velocity
            if ((dwhere .ge. 5.0_DP) .and. (dwhere .le. 6.0_DP)) then
              Dvalues(1) = dy *(1.0_DP - dy)
            end if

            if ((dwhere .ge. 1.0_DP) .and. (dwhere .le. 2.0_DP)) then
              Dvalues(1) = dy *(2.0_DP - dy)/8.0_DP
            end if

          case (2) ! Y-velocity
            ! Nothing to do here.

          case (3) ! Pressure
            ! Nothing to do here.
          end select
       
       case (2)
          ! Nothing to do here.
       end select      
    
    case (6,7)
      ! Anlytic polynomial function
      cu = rcollection%DquickAccess(10)
      pi = 3.1415926535897932_DP
      dC = rcollection%DquickAccess(9)      
      select case (rboundaryRegion%IBOUNDCOMPIDX)
      case (1)
          select case (icomponent)
          case (1) ! X-velocity
!            Dvalues(1) = cu*exp(dx)*cos(cu*dy)
            Dvalues(1) = 0.0_DP
!            Dvalues(1) = cos(cu*pi*dy)
          case (2) ! Y-velocity
!            Dvalues(1) = -exp(dx)*sin(cu*dy)
            Dvalues(1) = 0.0_DP
!            Dvalues(1) = cos(cu*pi*dx)
          case (3) ! Pressure
!             Dvalues(1) = ( 1.0_DP-exp(-dC*dx) ) * sin(2*pi*dy)
             Dvalues(1) = dC*(dx**3 - dy**3)
!             Dvalues(1) = cos(dC*pi*dx)
          end select
       
      case (2)
          ! Nothing to do here.
      end select   


    case (8)
      ! Split channel
       select case (rboundaryregion%iboundcompidx)
       case (1)
          select case (icomponent)
          case (1) ! x-velocity
           ! Outflow
            if ((dwhere .ge. 0.0_DP) .and. (dwhere .le. 1.0_DP)) then
              dvalues(1) = (0.5_DP - dy)*(0.5_DP + dy)
            end if
           ! Inflow  
            if ((dwhere .ge. 5.0_DP) .and. (dwhere .le. 6.0_DP)) then
              dvalues(1) = (0.5_DP - dy)*(0.5_DP + dy)
            end if

          case (2) ! y-velocity
            ! nothing to do here.

          case (3) ! pressure
            ! nothing to do here.
          end select
       
       case (2)
          ! nothing to do here.
       end select


    case (9)
      ! Fully developed power law flow
      ! u/u_ave = 2n+1/n+1 * ( 1 - y^(n+1)/n )
      ! y = [0,1]
      ! n = r-1
      
      select case (icomponent)
      case (1) ! X-velocity
        if ((dwhere .ge. 3.0_DP) .and. (dwhere .le. 4.0_DP)) then
          n = rcollection%DquickAccess(18)-1.0_DP
          y = 4.0_DP-dwhere
          Dvalues(1) = (2.0_DP*n+1.0_DP)/(n+1.0_DP)* &
                       ( 1.0_DP - y**((n+1.0_DP)/n) )
        end if

        if (dwhere .eq. 0.0_DP) then
          n = rcollection%DquickAccess(18)-1.0_DP
          y = dwhere
          Dvalues(1) = (2.0_DP*n+1.0_DP)/(n+1.0_DP)* &
                       ( 1.0_DP - y**((n+1.0_DP)/n) )
        end if

      case (2) ! Y-velocity
        ! Nothing to do here.

      case (3) ! Pressure
        ! Nothing to do here.
      end select

    
    case default
      ! Un-known problem
      call output_line ("Unknown problem.", OU_CLASS_WARNING, OU_MODE_STD, &
                        "ls_BCs_Dirichlet_One")
      call sys_halt()
    end select

!    ! Step Flow
!     select case (rboundaryRegion%IBOUNDCOMPIDX)
!     case (1)
!        select case (icomponent)
!        case (1) ! X-velocity
!          if ((dwhere .ge. 3.0_DP) .and. (dwhere .le. 4.0_DP)) then
!            Dvalues(1) = (dy-0.5_DP) * (0.5319_DP - (dy-0.5_DP))
!          end if
!
!        case (2) ! Y-velocity
!          ! Nothing to do here.
!
!        case (3) ! Pressure
!          ! Nothing to do here.
!        end select
!     
!     case (2)
!        ! Nothing to do here.
!     end select


!    ! Proot 2006
!     select case (rboundaryregion%iboundcompidx)
!     case (1)
!        select case (icomponent)
!        case (1) ! x-velocity
!!         ! Inflow  
!!          if ((dwhere .ge. 3.0_DP) .and. (dwhere .le. 4.0_DP)) then
!!            dvalues(1) =  (dy + 0.75_DP) * (dy - 0.75_DP)
!!          end if
!          Dvalues(1) = 1.0_DP
!        case (2) ! y-velocity
!          ! nothing to do here.
!
!        case (3) ! pressure
!          ! nothing to do here.
!        end select
!     
!     case (2)
!        ! nothing to do here.
!     end select
     
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine analyt_project (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
   ! A simple circle
    Dcoefficients (1,:,:) = SQRT( (Dpoints(1,:,:) - 0.1_DP)**2 + (Dpoints(2,:,:) - 0.1_DP)**2 ) - 0.1_DP
   
!   !  The Osher-Fedkiw favorite star :D
!   Dcoefficients (1,:,:) = DSQRT(Dpoints(1,:,:)**2 + Dpoints(2,:,:)**2)- &
!   0.25_DP*( DCOS(7.0_DP*DATAN2(Dpoints(2,:,:),Dpoints(1,:,:))) + 2.0_DP)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine analyt_project_U (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
   !  Simple
   real(DP) :: dPI
   
!    Dcoefficients (1,:,:) = Dpoints(2,:,:)*(1.0_DP - Dpoints(2,:,:))
   Dcoefficients (1,:,:) = 1.0_DP !-Dpoints(2,:,:) 
   
!   dPI = 3.1415926535897932_DP
!   Dcoefficients (1,:,:) = -DSIN(dPI*Dpoints(1,:,:))*DSIN(dPI*Dpoints(1,:,:))*&
!   DSIN(2*dPI*Dpoints(2,:,:))

  end subroutine
  

  ! ***************************************************************************

!<subroutine>

  subroutine analyt_project_V (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
   !  Simple
   real(DP) :: dPI

  Dcoefficients (1,:,:) = 1.0_DP !Dpoints(1,:,:)
   
!   dPI = 3.1415926535897932_DP
!   Dcoefficients (1,:,:) = DSIN(dPI*Dpoints(2,:,:))*DSIN(dPI*Dpoints(2,:,:))*&
!   DSIN(2*dPI*Dpoints(1,:,:))

  end subroutine


  ! ***************************************************************************   


!<subroutine>

  subroutine getReferenceFunction_2D (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use scalarpde
  use domainintegration
  use collection
  
  
!<description>
  ! This subroutine is called during the calculation of errors. It has to compute
  ! the (analytical) values of a function in a couple of points on a couple
  ! of elements. These values are compared to those of a computed FE function
  ! and used to calculate an error.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)                      :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>
  
!</subroutine>
         
  select case (cderivative)
   case (DER_FUNC)
     ! Pressure
     Dvalues (:,:) = 1.0_DP*(Dpoints(1,:,:)**3 - Dpoints(2,:,:)**3 - 0.5_DP)
   case (DER_DERIV_X)

   case (DER_DERIV_Y)

   case DEFAULT
     ! Unknown. Set the result to 0.0.
     Dvalues = 0.0_DP
  end select 

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValuesFBC_2D (Icomponents,rdiscretisation,&
                                      Revaluation, rcollection)
  
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
  integer, dimension(:), intent(in)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_blockDiscretisation), intent(in)                     :: rdiscretisation
  
  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional                 :: rcollection

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

      ! local variables
      real(DP) :: dx, dy
      real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
      type(t_triangulation), pointer :: p_rtriangulation
      integer :: ipoint,idx
 
      ! Get the triangulation array for the point coordinates
      p_rtriangulation => rdiscretisation%RspatialDiscr(1)%p_rtriangulation
      call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                     p_DvertexCoordinates)

      
      ! Loop through the points where to evaluate:
      do idx = 1,Revaluation(1)%nvalues
      
        ! Get x- and y-coordinate
        dx = Revaluation(1)%p_Dwhere(1,idx)
        dy = Revaluation(1)%p_Dwhere(2,idx)
        
        ! Point inside?
        if (dx .eq. 0.0_DP .and. dy .eq. 0.0_DP) then
        
          ! Denote in the p_Iinside array that we prescribe a value here:
          Revaluation(1)%p_Iinside (idx) = 1
          
          ! We prescribe 0.0 as Dirichlet value here.
          Revaluation(1)%p_Dvalues (idx,1) = -1.0_DP*0.5_DP
        
        end if
        
      end do
    
  end subroutine
    
end module
