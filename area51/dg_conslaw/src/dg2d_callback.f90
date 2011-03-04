!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the Poisson problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in 'intf_xxxx.inc' files.
!#
!# --- 2D version ---
!#
!# 1.) coeff_Laplace_2D
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 2.) coeff_RHS_2D
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation. 2D case, Q2 bubble solution.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) coeff_RHS_Sin2D
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation. 2D case, sinus bubble solution.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 4.) getBoundaryValues_2D
!#     -> Returns analytic values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# 5.) getBoundaryValuesFBC_2D
!#     -> Returns analytic values in the inner of the domain on
!#        fictitious boundary objects
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcfassembly.inc'
!#
!# 6.) getBoundaryValuesMR_2D
!#     -> Returns discrete values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_discretebc.inc'
!#
!# 7.) getReferenceFunction_2D
!#     -> Returns the values of the analytic function and its derivatives,
!#        corresponding to coeff_RHS_2D, Q2 bubble solution.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 8.) getReferenceFunction_Sin2D
!#     -> Returns the values of the analytic function and its derivatives,
!#        corresponding to coeff_RHS_Sin2D, sinus bubble solution.
!#     -> Is only used for the postprocessing to calculate the $L_2$- and
!#        $H_1$-error of the FE function in comparison to the analytic
!#        function
!#
!# 9.) gethadaptMonitorFunction_2D
!#     -> Controls the grid adaption strategy in poisson2d_method1_hadapt.
!#
!# </purpose>
!##############################################################################

module dg2d_callback

  use fsystem
  use storage
  use genoutput
  use linearsolver
  use boundary
  use triangulation
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use derivatives
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use bcassembly
  use element
  use feevaluation
  use domainintegration
  use fparser
  use dg2d_problem
  use triasearch
  use linearsystemscalar
  
  implicit none

contains

! ***************************************************************************
  !<subroutine>

  subroutine coeff_Laplace_2D (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
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
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
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
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    Dcoefficients = 1.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_2D (rdiscretisation,rform, &
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

    !    u(x,y) = 16*x*(1-x)*y*(1-y)
    ! => f(x,y) = 32 * (y*(1-y)+x*(1-x))
    Dcoefficients (1,:,:) = 32.0_DP * &
                    ( Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:)) + &
                      Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) )

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_2D (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
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

integer :: iunit,i,j,iunit1,iunit2
real(dp),dimension(4000001) :: Dreference, Drefx
real(dp) :: r,n,r1
real(dp) :: dx,dy, dr, dt, drad, dh

integer :: ielement, ipoint



!! Euler: Isentropic vortex
!dt = 1.0_dp
!
!do i = 1, size(Dvalues,1)
!    do j = 1, size(Dvalues,2)
!       dx = Dpoints(1,i,j)
!       dy = Dpoints(2,i,j)
!       drad = sqrt((dx-dt-5.0_dp)**2.0_dp + dy*dy)
!       Dvalues(i,j) = (1.0_dp-0.4_dp/(16.0_dp*1.4_dp*SYS_pi**2.0_dp)*25.0_dp*exp(2.0_dp*(1.0_dp-drad*drad)))**(1.0_dp/0.4_dp)
!    end do
!    end do







! Euler: SODs shock tube

iunit1 = sys_getFreeUnit()
open(iunit1, file='./SODx.out')
iunit2 = sys_getFreeUnit()
open(iunit2, file='./SODrho.out')

do i = 1, 1601
  read(iunit1,*) Drefx(i)
  read(iunit2,*) Dreference(i)
end do

close(iunit1)
close(iunit2)


  select case (cderivative)
  case (DER_FUNC)
    
                    
    do ielement = 1, nelements
    do ipoint = 1, npointsPerElement
      dx = Dpoints(1,ipoint,ielement)
      do i = 1, 1601
        if (dx < Drefx(i)) then
          exit
        end if
      end do
      
      dh = (dx-Drefx(i-1))/(Drefx(i)-Drefx(i-1))
            
      Dvalues(ipoint,ielement) =(1.0_dp-dh)* Dreference(i-1) +(dh)* Dreference(i)
      
    end do
    end do
                              
  case (DER_DERIV_X)
  write(*,*) 'Error in calculating L2-error'
    
  case (DER_DERIV_Y)
  write(*,*) 'Error in calculating L2-error'
    
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select








!! sin
!
!  select case (cderivative)
!  case (DER_FUNC)
!    
!                    
!    do i = 1, size(Dvalues,1)
!    do j = 1, size(Dvalues,2)
!       Dvalues(i,j) = sin(SYS_PI*(Dpoints(1,i,j) - Dpoints(2,i,j)))
!    end do
!    end do
!                              
!  case (DER_DERIV_X)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case (DER_DERIV_Y)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case DEFAULT
!    ! Unknown. Set the result to 0.0.
!    Dvalues = 0.0_DP
!  end select


!! Source Term 1
!  select case (cderivative)
!  case (DER_FUNC)
!    
!                    
!    do i = 1, size(Dvalues,1)
!    do j = 1, size(Dvalues,2)
!       Dvalues(i,j) = 1.0_dp+Dpoints(1,i,j) + Dpoints(2,i,j)
!    end do
!    end do
!                              
!  case (DER_DERIV_X)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case (DER_DERIV_Y)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case DEFAULT
!    ! Unknown. Set the result to 0.0.
!    Dvalues = 0.0_DP
!  end select
  
  
  
  
  
  
!! Source Term 2
!  select case (cderivative)
!  case (DER_FUNC)
!    
!                    
!    do i = 1, size(Dvalues,1)
!    do j = 1, size(Dvalues,2)
!       Dvalues(i,j) = (1.0_dp+Dpoints(1,i,j)**4.0_dp+Dpoints(2,i,j)**4.0_dp)
!    end do
!    end do
!                              
!  case (DER_DERIV_X)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case (DER_DERIV_Y)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case DEFAULT
!    ! Unknown. Set the result to 0.0.
!    Dvalues = 0.0_DP
!  end select  
  
  
!  ! Quarter circular convection
!  select case (cderivative)
!  case (DER_FUNC)
!    
!                    
!    do i = 1, size(Dvalues,1)
!    do j = 1, size(Dvalues,2)
!      dx = Dpoints(1,i,j)
!      dy = Dpoints(2,i,j)
!       
!       dr = abs( sqrt((dx-1.0_dp)*(dx-1.0_dp)+(dy-0.0_dp)*(dy-0.0_dp)) -0.5_dp)
!!       dr = abs( sqrt((dx-1.0_dp-dy)*(dx-1.0_dp-dy)) -0.5_dp)
!       
!       Dvalues(i,j) = 0.3_DP*exp(-25.0_dp*dr*dr)
!          
!    end do
!    end do
!                              
!  case (DER_DERIV_X)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case (DER_DERIV_Y)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case DEFAULT
!    ! Unknown. Set the result to 0.0.
!    Dvalues = 0.0_DP
!  end select  
  
  
  
  

!    ! Compare to less refined solution
!
!      ! An object for saving the domain:
!    type(t_boundary) :: rboundary
!    
!    ! An object for saving the triangulation on the domain
!    type(t_triangulation) :: rtriangulation
!    
!    type(t_vectorScalar) :: rsolLoaded
!    
!    ! An object specifying the discretisation.
!    ! This contains also information about trial/test functions,...
!    type(t_blockDiscretisation) :: rdiscretisation2
!    
!    integer :: NLMAX
!    
!    real(DP), dimension(:), pointer :: p_Ddata
!    
!    integer :: ielementtype, iel,i1,i2,i3
!    
!    character (LEN=SYS_STRLEN) :: sofile
!    
!    integer :: ilength
!  
!  
!  
!    
!    
!    ielementtype = EL_DG_T2_2D
!    
!    NLMAX = 8
!    
!    sofile = 'l8.data'
!    
!    call boundary_read_prm(rboundary, './pre/quad.PRM')
!    
!    ! Now read in the basic triangulation.
!    call tria_readTriFile2D (rtriangulation, './pre/quad.TRI', rboundary, .true.)  
!    
!    
!      
!    
!    ! Refine it.
!    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
!    
!    ! And create information about adjacencies and everything one needs from
!    ! a triangulation.
!    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
!  
!    ! Now we can start to initialise the discretisation. At first, set up
!    ! a block discretisation structure that specifies the blocks in the
!    ! solution vector. In this simple problem, we only have one block.
!    call spdiscr_initBlockDiscr (rdiscretisation2,1,&
!                                 rtriangulation, rboundary)
!    
!    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
!    ! structures for every component of the solution vector.
!    ! Initialise the first element of the list to specify the element
!    ! and cubature rule for this solution component:
!    call spdiscr_initDiscr_simple (rdiscretisation2%RspatialDiscr(1), &
!                                   ielementType,CUB_G5X5,rtriangulation, rboundary)
!                 
!    
!    
!    call lsyssc_createVecByDiscr (rdiscretisation2%RspatialDiscr(1),rsolloaded ,.true.,ST_DOUBLE)
!    
!    
!    
!  ! Get pointers to the data form the truangulation
!  call lsyssc_getbase_double(rsolloaded,p_Ddata)
!
!  ! Get the length of the data array
!  ilength = size(p_Ddata,1)  
!  
!  
!  ! ************ WRITE TO FILE PHASE *******************
!  
!  iunit = sys_getFreeUnit()
!  open(iunit, file=trim(sofile))
! 
!  
!  read(iunit,'(I10)') ilength
!
!  do i=1, ilength
!    read(iunit,'(E25.16E3)') p_Ddata(i)
!  end do
!  
!  close(iunit)
!    
!    
!    
!    do iel = 1, nelements
!      ielement=0
!      call tsrch_getElem_raytrace2D (&
!            Dpoints(:,1,iel),rtriangulation,ielement,i1,i2,i3,100000)
!      call fevl_evaluate_mult1 (DER_FUNC, Dvalues(:,iel), rsolloaded, ielement, &
!                                Dpoints=Dpoints(:,:,iel))
!    end do
!
!            
!    
!    call lsyssc_releaseVector (rsolloaded)
!        
!    call spdiscr_releaseBlockDiscr(rdiscretisation2)
!
!    ! Release the triangulation. 
!    call tria_done (rtriangulation)
!        
!    ! Finally release the domain, that is it.
!    call boundary_release (rboundary)






















!! Circular dambreak
!
!iunit = sys_getFreeUnit()
!
!open(iunit, file='h')
!
!do i = 1, 10000
!  read(iunit,*) Dreference(i)
!end do
!
!close(iunit)
!
!
!  select case (cderivative)
!  case (DER_FUNC)
!    
!                    
!    do ielement = 1, nelements
!    do ipoint = 1, npointsPerElement
!      r = sqrt(Dpoints(1,ipoint,ielement)*Dpoints(1,ipoint,ielement)+Dpoints(2,ipoint,ielement)*Dpoints(2,ipoint,ielement))
!      n = 1+1000*r
!            
!       Dvalues(ipoint,ielement) =(1.0_dp-(n-real(int(n))))* Dreference(int(n)) +(n-real(int(n)))* Dreference(int(n)+1)
!
!      !Dvalues(i,j) =Dreference(int(n))
!      !write(*,*)Dreference(int(n))
!      !Dvalues(i,j) =1.0_dp
!      
!    end do
!    end do
!                              
!  case (DER_DERIV_X)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case (DER_DERIV_Y)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case DEFAULT
!    ! Unknown. Set the result to 0.0.
!    Dvalues = 0.0_DP
!  end select




!! Smooth solution
!
!iunit = sys_getFreeUnit()
!
!open(iunit, file='h')
!
!
!do i = 1, 400000
!  read(iunit,*) Dreference(i)
!
!end do
!
!close(iunit)
!
!  select case (cderivative)
!  case (DER_FUNC)
!    
!                   
!    do ielement = 1, nelements
!    do ipoint = 1, npointsPerElement
!      r = sqrt((Dpoints(1,ipoint,ielement)-0.5_dp)**2.0_dp+(Dpoints(2,ipoint,ielement)-0.5_dp)**2.0_dp)
!      n = 200000+1+100000*r             
!                    
!    
!      
!     
!      
!       Dvalues(ipoint,ielement) =(1.0_dp-(n-real(int(n))))* Dreference(int(n)) +(n-real(int(n)))* Dreference(int(n)+1)
!
!      !Dvalues(i,j) =Dreference(int(n))
!      !write(*,*)Dreference(int(n))
!      !Dvalues(i,j) =1.0_dp
!      
!    end do
!    end do
!                              
!  case (DER_DERIV_X)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case (DER_DERIV_Y)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case DEFAULT
!    ! Unknown. Set the result to 0.0.
!    Dvalues = 0.0_DP
!  end select



!! x^2+y^2
!
!  select case (cderivative)
!  case (DER_FUNC)
!    
!                    
!    do i = 1, size(Dvalues,1)
!    do j = 1, size(Dvalues,2)
!       Dvalues(i,j) = (Dpoints(1,i,j)-0.5_dp)**2.0_dp + (Dpoints(2,i,j)-0.5_dp)**2.0_dp
!    end do
!    end do
!                              
!  case (DER_DERIV_X)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case (DER_DERIV_Y)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case DEFAULT
!    ! Unknown. Set the result to 0.0.
!    Dvalues = 0.0_DP
!  end select
  
  
  
!  ! x^3+y^3
!
!  select case (cderivative)
!  case (DER_FUNC)
!    
!                    
!    do i = 1, size(Dvalues,1)
!    do j = 1, size(Dvalues,2)
!       Dvalues(i,j) = (Dpoints(1,i,j)-0.5_dp)**3.0_dp + (Dpoints(2,i,j)-0.5_dp)**3.0_dp
!    end do
!    end do
!                              
!  case (DER_DERIV_X)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case (DER_DERIV_Y)
!  write(*,*) 'Error in calculating L2-error'
!    
!  case DEFAULT
!    ! Unknown. Set the result to 0.0.
!    Dvalues = 0.0_DP
!  end select



!            ! Rotating hill
!            Dvalues(:,:)= 0.0_dp
!            do i = 1, npointsperelement
!            do j = 1, nelements
!            r1 = sqrt(((Dpoints(1,i,j)-0.5_dp)**2+(Dpoints(2,i,j)-0.5_dp)**2))
!            Dvalues(i,j)=0.5_dp*exp(-10.0_dp*r1)
!            end do
!            end do

!            ! Zalesak test
!            Dvalues(:,:)= 0.0_dp
!            do i = 1, npointsperelement
!            do j = 1, nelements
!            r1 = sqrt(((Dpoints(1,i,j)-0.5_dp)**2+(Dpoints(2,i,j)-0.75_dp)**2))
!            if ((r1.le.0.15_dp).and.((abs(Dpoints(1,i,j)-0.5).ge.0.025_dp).or.(Dpoints(2,i,j).ge.0.85_dp))) Dvalues(i,j)=1.0_dp
!            r1 = sqrt(((Dpoints(1,i,j)-0.5_dp)**2+(Dpoints(2,i,j)-0.25_dp)**2))
!            if (r1.le.0.15_dp) Dvalues(i,j)=1.0_dp-r1/0.15_dp
!            r1 = sqrt(((Dpoints(1,i,j)-0.25_dp)**2+(Dpoints(2,i,j)-0.5_dp)**2))
!            if (r1.le.0.15_dp) Dvalues(i,j)=0.25_dp*(1.0_dp+cos(SYS_PI*min(1.0_dp,r1/0.15_dp)))
!            end do
!            end do
            
!            ! Zalesak test steady
!            Dvalues(:,:)= 0.0_dp
!            do i = 1, npointsperelement
!            do j = 1, nelements
!            r1 = sqrt(((Dpoints(1,i,j)-0.5_dp)**2+(Dpoints(2,i,j)-0.25_dp)**2))
!            if (r1.le.0.15_dp) Dvalues(i,j)=1.0_dp-r1/0.15_dp
!            r1 = sqrt(((Dpoints(1,i,j)-0.25_dp)**2+(Dpoints(2,i,j)-0.5_dp)**2))
!            if (r1.le.0.15_dp) Dvalues(i,j)=0.25_dp*(1.0_dp+cos(SYS_PI*min(1.0_dp,r1/0.15_dp)))
!            end do
!            end do
            
!            ! Zalesak test smooth
!            Dvalues(:,:)= 0.0_dp
!            do i = 1, npointsperelement
!            do j = 1, nelements
!            r1 = sqrt(((Dpoints(1,i,j)-0.25_dp)**2+(Dpoints(2,i,j)-0.5_dp)**2))
!            if (r1.le.0.15_dp) Dvalues(i,j)=0.25_dp*(1.0_dp+cos(SYS_PI*min(1.0_dp,r1/0.15_dp)))
!            end do
!            end do
  

  end subroutine
  
  
  
  
! ***************************************************************************

!<subroutine>

  subroutine getCompareFunction_2D (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
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

integer :: iel



  select case (cderivative)
  case (DER_FUNC)
                              
      do iel = 1, size(Dvalues,2)
      
      call fevl_evaluate (DER_FUNC, Dvalues(:,iel),&
                 rcollection%p_rvectorQuickAccess1%RvectorBlock(1), Dpoints(:,:,iel))
                 
    end do
                              
  case (DER_DERIV_X)
  write(*,*) 'Error in calculating L2-error'
    
  case (DER_DERIV_Y)
  write(*,*) 'Error in calculating L2-error'
    
  case DEFAULT
    write(*,*) 'Error in calculating L2-error'
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_Sin2D (rdiscretisation,rform, &
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

    !    u(x,y) = SIN(PI * x) * SIN(PI * y)
    ! => f(x,y) = 2 * PI^2 * SIN(PI * x) * SIN(PI * y)
    Dcoefficients (1,:,:) = 2.0_DP * SYS_PI**2 &
                          * sin(SYS_PI * Dpoints(1,:,:)) &
                          * sin(SYS_PI * Dpoints(2,:,:))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getReferenceFunction_Sin2D (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
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
    ! u(x,y) = SIN(PI * x) * SIN(PI * y)
    Dvalues (:,:) = sin(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
  case (DER_DERIV_X)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_x(x,y) = PI * COS(PI * x) * SIN(PI * y)
    Dvalues (:,:) = SYS_PI * cos(SYS_PI*Dpoints(1,:,:)) * sin(SYS_PI*Dpoints(2,:,:))
  case (DER_DERIV_Y)
    !    u(x,y)   = SIN(PI * x) * SIN(PI * y)
    ! => u_y(x,y) = PI * SIN(PI * x) * COS(PI * y)
    Dvalues (:,:) = SYS_PI * sin(SYS_PI*Dpoints(1,:,:)) * cos(SYS_PI*Dpoints(2,:,:))
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  

  end subroutine

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

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    ! REAL(DP) :: dx,dy
    !
    ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
    !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP
  
  end subroutine

  ! ***************************************************************************
  ! Only for poisson2d_method1_fbc: Values in a fictitious boundary component:

!<subroutine>

  subroutine getBoundaryValuesFBC_2D (Icomponents,rdiscretisation,&
                                      Revaluation, rcollection)
  
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
      real(DP) :: ddistance, dxcenter, dycenter, dradius, dx, dy
      real(DP), dimension(:,:), pointer :: p_DvertexCoordinates
      type(t_triangulation), pointer :: p_rtriangulation
      integer :: ipoint,idx

      
      ! Just make sure we are evaluating in the corners.
      if (Revaluation(1)%cinfoNeeded .ne. DISCFBC_NEEDFUNC) then
        print *,'FBC: only corner evaluation supported at the moment!'
        stop
      end if
      
      ! Get the triangulation array for the point coordinates
      p_rtriangulation => rdiscretisation%RspatialDiscr(1)%p_rtriangulation
      call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                     p_DvertexCoordinates)

      ! Definition of the circle
      dxcenter = 0.4
      dycenter = 0.4
      dradius  = 0.25
      
      ! Loop through the points where to evaluate:
      do idx = 1,Revaluation(1)%nvalues
      
        ! Get the number of the point to process
        ipoint = Revaluation(1)%p_Iwhere(idx)
        
        ! Get x- and y-coordinate
        dx = p_DvertexCoordinates(1,ipoint)
        dy = p_DvertexCoordinates(2,ipoint)
        
        ! Get the distance to the center
        ddistance = sqrt( (dx-dxcenter)**2 + (dy-dycenter)**2 )
        
        ! Point inside?
        if (ddistance .le. dradius) then
        
          ! Denote in the p_Iinside array that we prescribe a value here:
          Revaluation(1)%p_Iinside (idx) = 1
          
          ! We prescribe 0.0 as Dirichlet value here.
          Revaluation(1)%p_Dvalues (idx,1) = 0.0_DP
        
        end if
        
      end do

      ! Definition of a 2nd circle
      dxcenter = 0.75
      dycenter = 0.75
      dradius  = 0.1
      
      ! Loop through the points where to evaluate:
      do idx = 1,Revaluation(1)%nvalues
      
        ! Get the number of the point to process
        ipoint = Revaluation(1)%p_Iwhere(idx)
        
        ! Get x- and y-coordinate
        dx = p_DvertexCoordinates(1,ipoint)
        dy = p_DvertexCoordinates(2,ipoint)
        
        ! Get the distance to the center
        ddistance = sqrt( (dx-dxcenter)**2 + (dy-dycenter)**2 )
        
        ! Point inside?
        if (ddistance .le. dradius) then
        
          ! Denote in the p_Iinside array that we prescribe a value here:
          Revaluation(1)%p_Iinside (idx) = 1
          
          ! We prescribe 0.0 as Dirichlet value here.
          Revaluation(1)%p_Dvalues (idx,1) = 0.0_DP
        
        end if
        
      end do
    
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine getBoundaryValuesMR_2D (Icomponents,rdiscretisation,rmeshRegion,&
                                      cinfoNeeded,Dcoords,Dvalues,rcollection)
  
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
  integer, dimension(:), intent(in)                           :: Icomponents

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
  
  ! Mesh region that is currently being processed.
  type(t_meshRegion), intent(in)                              :: rmeshRegion

  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  integer, intent(in)                                         :: cinfoNeeded
  
  ! The coordinates of the point for which the boundary values are to be
  ! calculated.
  real(DP), dimension(:), intent(in)                          :: Dcoords

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(inout), optional                 :: rcollection

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
  real(DP), dimension(:), intent(out)                         :: Dvalues
!</output>
  
!</subroutine>

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP

  end subroutine

  ! ***************************************************************************
  ! Only for poisson2d_method1_hadapt: Monitor function for adaptive grid refinement
  
!<subroutine>

  subroutine gethadaptMonitorFunction_2D(rtriangulation,rsolution,ieltype,&
      ierrorestimator,rindicator)
  
    use pprocgradients
    use pprocerror

!<description>
  ! This routine defines a 'monitor function' for the adaptive grid refinement
  ! with the h-adaptivity refinement strategy. rindicator is a vector with
  ! NEL entries for all the elements in the triangulation. The routine must
  ! fill each entry with a value that tells the h-adaptivity routines whether
  ! to refine that element or not.
!</descrition>
  
!<input>
    ! The triangulation structure of the underlying mesh which is to be refined
    type(t_triangulation), intent(in) :: rtriangulation

    ! The solution vector
    type(t_vectorScalar), intent(inout)  :: rsolution
    
    ! The type of element used for the FE solution
    integer(I32), intent(in) :: ieltype
    
    ! The type of error estimator
    integer, intent(in) :: ierrorestimator
!</input>
    
!</inputoutput>
    ! An indicator vector. Entry i in the vector rindicatir that tells the 
    ! mesh adaption routines whether to refine element i or to do coarsening
    ! with it. A value > 1.0 will refine element i, a value < 0.01 will result
    ! in coarsening -- as specified during the initialisation of the
    ! mesh refinement in the main program.
    type(t_vectorScalar), intent(inout) :: rindicator

!</subroutine>

    ! local variables
    type(t_vectorBlock)         :: rgradient,rgradientRef
    type(t_blockDiscretisation) :: rdiscrBlock,rdiscrBlockRef
    real(DP)                    :: dsolutionError,dgradientError,daux

    ! Initialise block discretisations
    call spdiscr_initBlockDiscr (rdiscrBlock,2,&
        rtriangulation, rsolution%p_rspatialDiscr%p_rboundary)
    call spdiscr_initBlockDiscr (rdiscrBlockRef,2,&
        rtriangulation, rsolution%p_rspatialDiscr%p_rboundary)

    ! What kind of element type is used for the FE solution
    select case(ieltype)
    case(1)
      ! Initialise spatial discretisations for gradient with P0-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E000, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E000, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with P1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E001, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E001, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(2))
      
    case(2)
      ! Initialise spatial discretisations for gradient with P1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E001, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E001, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with P2-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E002, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E002, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(2))

    case(11)
      ! Initialise spatial discretisations for gradient with Q0-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E010, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E010, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with Q1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E011, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E011, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(2))

    case(13)
      ! Initialise spatial discretisations for gradient with Q1-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E011, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E011, SPDISC_CUB_AUTOMATIC, rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with Q2-elements
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E013, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveSimpleDiscrSc (rsolution%p_rspatialDiscr,&
          EL_E013, SPDISC_CUB_AUTOMATIC, rdiscrBlockRef%RspatialDiscr(2))

    case(-1)
      ! Initialise spatial discretisations for gradient with P0/Q0-elements
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialDiscr,&
          EL_E000, EL_E010, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlock%RspatialDiscr(1))
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialDiscr,&
          EL_E000, EL_E010, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlock%RspatialDiscr(2))
      
      ! Initialise spatial discretisations for reference gradient with P1/Q1-elements
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialDiscr,&
          EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlockRef%RspatialDiscr(1))
      call spdiscr_deriveDiscr_triquad (rsolution%p_rspatialDiscr,&
          EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC, SPDISC_CUB_AUTOMATIC, &
          rdiscrBlockRef%RspatialDiscr(2))
      
    case DEFAULT
      call output_line('Unsupproted element type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'getMonitorFunction')
      call sys_halt()
    end select
    
    ! Create block vector for gradient values
    call lsysbl_createVecBlockByDiscr (rdiscrBlock,   rgradient,    .true.)
    call lsysbl_createVecBlockByDiscr (rdiscrBlockRef,rgradientRef, .true.)

    ! Recover consistent gradient
    call ppgrd_calcGradient (rsolution, rgradient)

    ! Recover smoothed gradient
    select case(ierrorestimator)
    case (1)
      call ppgrd_calcGradient (rsolution, rgradientRef)

    case (2)
      call ppgrd_calcGradSuperPatchRecov (rsolution, rgradientRef, PPGRD_NODEPATCH)

    case (3)
      call ppgrd_calcGradSuperPatchRecov (rsolution, rgradientRef, PPGRD_ELEMPATCH)

    case (4)
      call ppgrd_calcGradSuperPatchRecov (rsolution, rgradientRef, PPGRD_FACEPATCH)

    case DEFAULT
      call output_line('Invalid type of error estimator!',&
          OU_CLASS_ERROR,OU_MODE_STD,'getMonitorFunction')
      call sys_halt()
    end select

    ! Compute gradient error
    call pperr_blockErrorEstimate(rgradient,rgradientRef,PPERR_L2ERROR,&
        dgradientError,relementError=rindicator)
    print *, "!!gradient error!! = ",dgradientError

    ! Compute L2-norm of solution
    call pperr_scalar(rsolution,PPERR_L2ERROR,dsolutionError)

    ! Prepare indicator for grid refinement/coarsening
    daux=sqrt((dsolutionError**2+dgradientError**2)/real(rindicator%NEQ,DP))
    call lsyssc_scaleVector(rindicator,1._DP/daux)
    
    ! Release temporal discretisation structure
    call spdiscr_releaseBlockDiscr(rdiscrBlock)
    call spdiscr_releaseBlockDiscr(rdiscrBlockRef)
    
    ! Release vectors
    call lsysbl_releaseVector(rgradient)
    call lsysbl_releaseVector(rgradientRef)
    
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  



  subroutine flux_dg_buildVectorScEdge2D_sim (&
              Dcoefficients,&
!              DsolVals,&
              IelementList,&
              normal,&
              !DpointsReal,&
              rintSubset,&
              rcollection )
              
  use collection
  
  !<input>
!  real(DP), dimension(:,:,:), intent(inout) :: DsolVals
  real(DP), dimension(:,:), intent(in) :: normal
!  real(DP), dimension(:,:,:), intent(in) :: DpointsReal
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  type(t_collection), intent(inout), target, optional :: rcollection
  integer, dimension(:) , intent(in) :: IelementList
  !</input>
  
  !<output>
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>


  integer :: ipoint, iel
  real(dp), dimension(2) :: Dvel
  real(DP) :: dvn ,dx, dy, dr
  real(dp), dimension(:,:,:), allocatable :: Dsolutionvalues
  
  ! Function parser
  type(t_fparser) :: rfparser
    
  character(LEN=*), dimension(2), parameter ::&
       cvariables = (/ (/'x'/), (/'y'/) /)
    
!  ! Initialise function parser
!  call fparser_init()
!  call fparser_create(rfparser, 1)
!  call fparser_parseFunction(rfparser, 1, trim(adjustl(rcollection%SquickAccess(1))), cvariables)






  ! If the solution (or its derivatives) has to be evaluated
  allocate(Dsolutionvalues(2,ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
  ! Get values on the one side of the edge
  call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(1), &
                                 rIntSubset(1), DER_FUNC, Dsolutionvalues, 1)
  ! Get values on the other side of the edge                               
  call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(1), &
                                 rIntSubset(2), DER_FUNC, Dsolutionvalues, 2)
               
  ! Now we have
  ! Dsolutionvalues(1,ipoint,iel) = DsolVals(ipoint,1,iel)
  ! and
  ! Dsolutionvalues(2,ipoint,iel) = DsolVals(ipoint,2,iel)
  ! except on the elements, which are on the other side and where there is no other side (on the boundary)



  
  
  do iel = 1, ubound(Dcoefficients,3)
  
    do ipoint= 1, ubound(Dcoefficients,2)
    
      !dx = DpointsReal(1,ipoint,iel)
      !dy = DpointsReal(2,ipoint,iel)
      
      dx = rintSubset(1)%p_DcubPtsReal(1,ipoint,iel)
      dy = rintSubset(1)%p_DcubPtsReal(2,ipoint,iel)
      
!      ! Zalesak
!      Dvel(1)=0.5_DP-dy
!      Dvel(2)=dx-0.5_DP
      
!      ! Constant velocity
!      Dvel(1)=1.0_dp
!      Dvel(2)=1.0_dp
      
      ! Steady circular convection
      Dvel(1)=dy
      Dvel(2)=1.0_DP-dx
    
      dvn = Dvel(1)*normal(1,iel)+Dvel(2)*normal(2,iel)
      
      ! Set initial condition on the inflow boundary
      !call fparser_evalFunction(rfparser, 1, rintSubset(1)%p_DcubPtsReal(:,ipoint,iel), DsolVals(ubound(Dcoefficients,2)-ipoint+1,2,iel))
      
!      ! BC for sinus
!      if ((dx.le.0.0_dp)) Dsolutionvalues(2,ubound(Dcoefficients,2)-ipoint+1,iel) = sin(-SYS_Pi*dy)
!      if ((dy.le.0.0_dp)) Dsolutionvalues(2,ubound(Dcoefficients,2)-ipoint+1,iel) = sin( SYS_Pi*dx)
      
      ! BC quarter circlular convection
!      if ((dx.le.0.0_dp)) Dsolutionvalues(2,ubound(Dcoefficients,2)-ipoint+1,iel) = 0.0_dp !0.3_DP*exp(-25.0_dp*(dy-0.5_dp)*(dy-0.5_dp))
!      if ((dy.le.0.0_dp)) Dsolutionvalues(2,ubound(Dcoefficients,2)-ipoint+1,iel) = 0.3_DP*exp(-25.0_dp*(dx-0.5_dp)*(dx-0.5_dp))

      dr = abs( sqrt((dx-1.00_dp)*(dx-1.00_dp)+(dy-0.0_dp)*(dy-0.0_dp)) -0.5_dp)
      !dr = abs( sqrt((dx-1.00_dp-1.0_dp/64.0_dp)*(dx-1.00_dp-1.0_dp/64.0_dp)+(dy-0.0_dp)*(dy-0.0_dp)) -0.5_dp)
      ! if ((dy.le.0.0_dp)) 
      if (IelementList(iel)==0) Dsolutionvalues(2,ubound(Dcoefficients,2)-ipoint+1,iel) = 0.3_DP*exp(-25.0_dp*dr*dr)
      
      

!      ! BC for quad: zero
!      if ((dx.le.0.0_dp).or.(dx.ge.1.0_dp).or.(dy.le.0.0_dp).or.(dy.ge.1.0_dp)) Dsolutionvalues(2,ubound(Dcoefficients,2)-ipoint+1,iel) = 0.0_dp
      
!      if ((dx.le.1.0_dp).and.(dy.le.0.0000000000001_dp)) then
!        dr = sqrt((dx-1.0_dp)**2.0_dp+dy*dy)
!        if ((0.2_dp.le.dr).and.(dr.le.0.4_dp)) then
!          DsolVals(2,ubound(Dcoefficients,2)-ipoint+1,iel) = 1.0_dp
!        elseif ((0.5_dp.le.dr).and.(dr.le.0.8_dp)) then
!          DsolVals(2,ubound(Dcoefficients,2)-ipoint+1,iel) = 0.25_dp*(1+cos(SYS_PI*(dr-0.65_dp)/0.15_dp))
!        end if

!      ! BC for circular convection
!      if ((dx.le.1.0_dp).and.(dy.le.0.0000000000001_dp)) then
!        dr = sqrt((dx-1.0_dp)**2.0_dp+dy*dy)
!        if ((0.2_dp.le.dr).and.(dr.le.0.4_dp)) then
!          Dsolutionvalues(2,ubound(Dcoefficients,2)-ipoint+1,iel) = 1.0_dp
!        elseif ((0.5_dp.le.dr).and.(dr.le.0.8_dp)) then
!          Dsolutionvalues(2,ubound(Dcoefficients,2)-ipoint+1,iel) = 0.25_dp*(1+cos(SYS_PI*(dr-0.65_dp)/0.15_dp))
!        end if
        
!      ! BC for circular sinus convection
!      if ((dx.le.1.0_dp).and.(dy.le.0.0000000000001_dp)) then
!        dr = sqrt((dx-1.0_dp)**2.0_dp+dy*dy)
!        
!          Dsolutionvalues(2,ubound(Dcoefficients,2)-ipoint+1,iel) = sin(dr)
!        
!        end if
        

    
      ! Upwind flux
      if (dvn.ge.0) then
        Dcoefficients(1,ipoint,iel) = dvn *Dsolutionvalues(1,ipoint,iel)
      else
        Dcoefficients(1,ipoint,iel) = dvn *Dsolutionvalues(2,ubound(Dcoefficients,2)-ipoint+1,iel)
      end if
      
      ! Centered Flux
      !Dcoefficients(1,ipoint,iel) = 0.5_dp* dvn *(DsolVals(1,ipoint,iel)+DsolVals(2,ubound(Dcoefficients,2)-ipoint+1,iel))
      
    end do ! ipoint
  end do ! iel
  
!  ! Release the function parser
!  call fparser_release(rfparser)
  
  

                                 
  deallocate(Dsolutionvalues)



  end subroutine
  
  
  
  
  
  
    ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_IC (rdiscretisation,rform, &
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
  
    real(DP) :: r1, dx, dy, dr

    integer :: iel, ipoint
    
    ! Function parser
    type(t_fparser) :: rfparser
    
    character(LEN=*), dimension(2), parameter ::&
         cvariables = (/ (/'x'/), (/'y'/) /)
         
    integer :: iunit,i,j
    real(dp),dimension(400000) :: Dreference
    real(dp) :: r,n

    
    ! Initialise function parser
    call fparser_init()
    call fparser_create(rfparser, 1)
    
    call fparser_parseFunction(rfparser, 1, trim(adjustl(rcollection%SquickAccess(1))), cvariables)
    

!    iunit = sys_getFreeUnit()
!
!    open(iunit, file='h')
!
!
!    do i = 1, 400000
!      read(iunit,*) Dreference(i)
!    end do
!
!    close(iunit)

    Dcoefficients (1,:,:) = 0.0_dp

    
    select case (rcollection%IquickAccess(1))
      case (1) ! Set height variable h
    
        do iel = 1, size(Dcoefficients,3)
          do ipoint = 1, size(Dcoefficients,2)
          
          dx = Dpoints(1,ipoint,iel)
          dy = Dpoints(2,ipoint,iel)
          
!          ! If you take a look at the limited and unlimited solution on lev 3, then you could get the idea to take different bounds
!          Dcoefficients (1,ipoint,iel) = 3.0_dp
!          if ((dx .le. 0.5_dp).and.(dy .ge. 0.5_dp)) Dcoefficients (1,ipoint,iel) = -2.0_dp+4.0_dp*dx+ 5.0_dp*dy
!          if ((abs(dx-0.5)>0.25).or.( abs(dy-0.5)>0.25)) Dcoefficients (1,ipoint,iel) = 1.0_dp

          
            ! Zero
            Dcoefficients (1,ipoint,iel)=0.0_dp
            
            
            
       
!       dr = abs( sqrt((dx-1.0_dp)*(dx-1.0_dp)+(dy-0.0_dp)*(dy-0.0_dp)) -0.5_dp)
!!       dr = abs( sqrt((dx-1.0_dp-dy)*(dx-1.0_dp-dy)) -0.5_dp)
!       
!       Dcoefficients (1,ipoint,iel) = 0.3_DP*exp(-25.0_dp*dr*dr)
            
          
!            ! Steady state
!            Dcoefficients (1,ipoint,iel)=1.0_dp
            
            ! Gauss-Hügel
            !Dcoefficients (1,ipoint,iel) = 0.1_dp*&
            !         exp(-5.0_dp*(Dpoints(1,ipoint,iel)**2+Dpoints(2,ipoint,iel)**2))
            
            ! Cone 1
            !Dcoefficients(1,ipoint,iel) = max(1.0_dp-3.0_dp*sqrt((Dpoints(1,ipoint,iel)**2+Dpoints(2,ipoint,iel)**2)),0.0_dp)
            
!            ! Parabola
!            Dcoefficients(1,ipoint,iel) = (Dpoints(1,ipoint,iel)-0.5_dp)**2+(Dpoints(2,ipoint,iel)-0.5_dp)**2
            
!            ! d/dx parabola
!            Dcoefficients(1,ipoint,iel) = 2.0_dp*Dpoints(1,ipoint,iel)-1.0_dp
            
!            ! d/dy parabola
!            Dcoefficients(1,ipoint,iel) = 2.0_dp*Dpoints(2,ipoint,iel)-1.0_dp
            
            ! Cone 2
            !Dcoefficients(1,ipoint,iel) = max(1.0_dp-6.0_dp*sqrt(((Dpoints(1,ipoint,iel)-0.3_dp)**2+(Dpoints(2,ipoint,iel)-0.3_dp)**2)),0.0_dp)
            
            
!            ! Rotating hill
!            r1 = sqrt(((Dpoints(1,ipoint,iel)-0.5_dp)**2+(Dpoints(2,ipoint,iel)-0.5_dp)**2))
!            Dcoefficients(1,ipoint,iel)=0.5_dp*exp(-10.0_dp*r1)

!            ! Rotating unsteady thing
!            r1 = sqrt(((Dpoints(1,ipoint,iel)-0.5_dp)**2+(Dpoints(2,ipoint,iel)-0.5_dp)**2))
!            if ((r1.le.0.35_dp).and.(abs(Dpoints(1,ipoint,iel)-0.5).ge.0.025_dp).and.(abs(Dpoints(2,ipoint,iel)-0.5).ge.0.025_dp)) Dcoefficients(1,ipoint,iel)=1.0_dp

            
!            ! Zalesak
!            r1 = sqrt(((Dpoints(1,ipoint,iel)-0.5_dp)**2+(Dpoints(2,ipoint,iel)-0.75_dp)**2))
!            if ((r1.le.0.15_dp).and.((abs(Dpoints(1,ipoint,iel)-0.5).ge.0.025_dp).or.(Dpoints(2,ipoint,iel).ge.0.85_dp))) Dcoefficients(1,ipoint,iel)=1.0_dp
!            r1 = sqrt(((Dpoints(1,ipoint,iel)-0.5_dp)**2+(Dpoints(2,ipoint,iel)-0.25_dp)**2))
!            if (r1.le.0.15_dp) Dcoefficients(1,ipoint,iel)=1.0_dp-r1/0.15_dp
!            r1 = sqrt(((Dpoints(1,ipoint,iel)-0.25_dp)**2+(Dpoints(2,ipoint,iel)-0.5_dp)**2))
!            if (r1.le.0.15_dp) Dcoefficients(1,ipoint,iel)=0.25_dp*(1.0_dp+cos(SYS_PI*min(1.0_dp,r1/0.15_dp)))

!            ! Zalesak steady
!            r1 = sqrt(((Dpoints(1,ipoint,iel)-0.5_dp)**2+(Dpoints(2,ipoint,iel)-0.25_dp)**2))
!            if (r1.le.0.15_dp) Dcoefficients(1,ipoint,iel)=1.0_dp-r1/0.15_dp
!            r1 = sqrt(((Dpoints(1,ipoint,iel)-0.25_dp)**2+(Dpoints(2,ipoint,iel)-0.5_dp)**2))
!            if (r1.le.0.15_dp) Dcoefficients(1,ipoint,iel)=0.25_dp*(1.0_dp+cos(SYS_PI*min(1.0_dp,r1/0.15_dp)))

!            ! Zalesak smooth
!            r1 = sqrt(((Dpoints(1,ipoint,iel)-0.25_dp)**2+(Dpoints(2,ipoint,iel)-0.5_dp)**2))
!            if (r1.le.0.15_dp) Dcoefficients(1,ipoint,iel)=0.25_dp*(1.0_dp+cos(SYS_PI*min(1.0_dp,r1/0.15_dp)))
            
!             ! Circular convection
!            r1 = sqrt((Dpoints(1,ipoint,iel)-1.0_dp)**2.0_dp+Dpoints(2,ipoint,iel)*Dpoints(2,ipoint,iel))
!            if ((0.2_dp.le.r1).and.(r1.le.0.4_dp)) Dcoefficients(1,ipoint,iel) = 1.0_dp
!            if ((0.5_dp.le.r1).and.(r1.le.0.8_dp)) Dcoefficients(1,ipoint,iel) = 0.25_dp*(1.0_dp+cos(SYS_PI*(r1-0.65_dp)/0.15_dp))
            
            ! Parser from .dat-file
            !call fparser_evalFunction(rfparser, 1, rdomainIntSubset%p_DcubPtsReal(:,ipoint,iel), Dcoefficients(1,ipoint,iel))
            
            ! Water hill
            Dcoefficients (1,ipoint,iel) = 1.0_dp + 0.3_dp*&
                     exp(-40.0_dp*((Dpoints(1,ipoint,iel)-0.5_dp)**2+(Dpoints(2,ipoint,iel)-0.5_dp)**2))
                     
!            ! Source 2
!            Dcoefficients (1,ipoint,iel) = 1.0_dp + 1.00_dp*&
!                     ( (Dpoints(1,ipoint,iel))**4.0_dp+(Dpoints(2,ipoint,iel))**4.0_dp) 
                     
!            ! Water hill Poly
!            r = sqrt(((Dpoints(1,ipoint,iel)-0.5_dp)**2+(Dpoints(2,ipoint,iel)-0.5_dp)**2))
!            Dcoefficients (1,ipoint,iel) = 38.40_dp*r*r*r-14.40_dp*r*r+1.3_dp
!            if (r>0.25_dp) Dcoefficients (1,ipoint,iel) = 1.0_dp

!            ! Water hill Poly 2
!            r = sqrt(((Dpoints(1,ipoint,iel)-0.5_dp)**2+(Dpoints(2,ipoint,iel)-0.5_dp)**2))
!            Dcoefficients (1,ipoint,iel) = 1.0_dp-0.2*r*r
                     
!            ! Water hill big
!            Dcoefficients (1,ipoint,iel) = 0.5_dp + 2.0_dp*exp(-0.75_dp*((Dpoints(1,ipoint,iel))**2+(Dpoints(2,ipoint,iel))**2) )
                     
!            ! 'Analytical' solution to circular dambreak
!            r = sqrt(Dpoints(1,ipoint,iel)*Dpoints(1,ipoint,iel)+Dpoints(2,ipoint,iel)*Dpoints(2,ipoint,iel))
!            n = 1+1000*r
!            Dcoefficients (1,ipoint,iel) =(1.0_dp-(n-real(int(n))))* Dreference(int(n)) +(n-real(int(n)))* Dreference(int(n)+1)

!            ! 'Analytical' smooth solution
!            r = sqrt((Dpoints(1,ipoint,iel)-0.5_dp)**2.0_dp+(Dpoints(2,ipoint,iel)-0.5_dp)**2.0_dp)
!            n = 200000+1+100000*r
!            Dcoefficients (1,ipoint,iel) =(1.0_dp-(n-real(int(n))))* Dreference(int(n)) +(n-real(int(n)))* Dreference(int(n)+1)
              
!            ! Water hill at corner
!            Dcoefficients (1,ipoint,iel) = 1.0_dp + 0.1_dp*&
!                     exp(-40.0_dp*((Dpoints(1,ipoint,iel)-0.25_dp)**2.0_dp+(Dpoints(2,ipoint,iel)-0.25_dp)**2.0_dp))
           
!           ! x^2+y^2 
!           Dcoefficients (1,ipoint,iel) = (Dpoints(1,ipoint,iel)-0.5_dp)**2.0_dp+(Dpoints(2,ipoint,iel)-0.5_dp)**2.0_dp

!           ! x^3+y^3
!           Dcoefficients (1,ipoint,iel) = (Dpoints(1,ipoint,iel)-0.5_dp)**3.0_dp+(Dpoints(2,ipoint,iel)-0.5_dp)**3.0_dp
            
!          ! Riemann problem
!          if (Dpoints(1,ipoint,iel)<0.5_dp) then
!            Dcoefficients (1,ipoint,iel)=1.5_dp
!          else
!            Dcoefficients (1,ipoint,iel)=1.0_dp
!          end if
          

!          ! Zylinders
!          if (Dpoints(1,ipoint,iel)<15.0_dp) then
!            Dcoefficients (1,ipoint,iel)=10.0_dp
!          else
!            Dcoefficients (1,ipoint,iel)=1.0_dp
!          end if

!          ! Bart Simpson
!          if (Dpoints(1,ipoint,iel)<0.3_dp) then
!            Dcoefficients (1,ipoint,iel)=1.5_dp
!          else
!            Dcoefficients (1,ipoint,iel)=1.0_dp
!          end if

!          ! Dam Break
!          if (Dpoints(1,ipoint,iel)<100.0_dp) then
!            Dcoefficients (1,ipoint,iel)=1.5_dp
!          else
!            Dcoefficients (1,ipoint,iel)=1.0_dp
!          end if
          
!          ! Dede Obstacle
!          if (Dpoints(1,ipoint,iel)<-0.3_dp) then
!            Dcoefficients (1,ipoint,iel)=1.5_dp
!          else
!            Dcoefficients (1,ipoint,iel)=1.0_dp
!          end if

!            ! Water reservoir
!            if ( sqrt((Dpoints(1,ipoint,iel)-0.5_dp)**2+(Dpoints(2,ipoint,iel)-0.5_dp)**2)<0.25_dp) then  
!              Dcoefficients (1,ipoint,iel)=1.5_dp
!            else
!              Dcoefficients (1,ipoint,iel)=1.0_dp
!            end if

!            ! Circular Dambreak
!            if ( sqrt((Dpoints(1,ipoint,iel)-0.0_dp)**2+(Dpoints(2,ipoint,iel)-0.0_dp)**2)<2.5_dp) then  
!              Dcoefficients (1,ipoint,iel)=2.5_dp
!            else
!              Dcoefficients (1,ipoint,iel)=0.5_dp
!            end if

!            ! Test on Quad
!            if ( ( (Dpoints(1,ipoint,iel)-0.5_dp) * (Dpoints(2,ipoint,iel)-0.5_dp) ) >=0.0_dp ) then
!              Dcoefficients (1,ipoint,iel)=2.0_dp
!            else
!              Dcoefficients (1,ipoint,iel)=1.0_dp
!            end if

!            ! Circular Dambreak small
!            if ( sqrt((Dpoints(1,ipoint,iel)-0.0_dp)**2+(Dpoints(2,ipoint,iel)-0.0_dp)**2)<2.5_dp) then  
!              Dcoefficients (1,ipoint,iel)=1.5_dp
!            else
!              Dcoefficients (1,ipoint,iel)=1.0_dp
!            end if

!            ! Test for symmetry
!            if (( abs(Dpoints(1,ipoint,iel)-0.5_dp)<0.25_dp).and.(abs(Dpoints(2,ipoint,iel)-0.5_dp)<0.25_dp )) then  
!              Dcoefficients (1,ipoint,iel)=1.5_dp
!            else
!              Dcoefficients (1,ipoint,iel)=1.0_dp
!            end if
 
           
          end do
        end do
      
      case (2) ! Set x-momentum
      
        Dcoefficients (1,:,:) = 0.0_dp
      
      case (3) ! Set y-momentum
        Dcoefficients (1,:,:) = 0.0_dp
        
    end select
    
    ! Release the function parser
    call fparser_release(rfparser)
                    
    !Dcoefficients (1,:,:) = 0.0_dp

  end subroutine
  
  
  
  
  
  ! ***************************************************************************

!<subroutine>

  subroutine Euler_coeff_RHS_IC (rdiscretisation,rform, &
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
  
    real(DP) :: r1, dx, dy, dr, drho, dt, drad, du, dv
    
    real(dp), parameter :: gamma = 1.4_dp

    integer :: iel, ipoint
    
    ! Function parser
    type(t_fparser) :: rfparser
    
    character(LEN=*), dimension(2), parameter ::&
         cvariables = (/ (/'x'/), (/'y'/) /)
         
    integer :: iunit,i,j
    real(dp),dimension(400000) :: Dreference
    real(dp) :: r,n

    
    ! Initialise function parser
    call fparser_init()
    call fparser_create(rfparser, 1)
    
    call fparser_parseFunction(rfparser, 1, trim(adjustl(rcollection%SquickAccess(1))), cvariables)
    

!    iunit = sys_getFreeUnit()
!
!    open(iunit, file='h')
!
!
!    do i = 1, 400000
!      read(iunit,*) Dreference(i)
!    end do
!
!    close(iunit)

    Dcoefficients (1,:,:) = 0.0_dp

    
    select case (rcollection%IquickAccess(1))
      case (1) ! Set density rho
    
        do iel = 1, size(Dcoefficients,3)
          do ipoint = 1, size(Dcoefficients,2)

          ! Get coordinates          
          dx = Dpoints(1,ipoint,iel)
          dy = Dpoints(2,ipoint,iel)
  
!          ! Zero
!          Dcoefficients (1,ipoint,iel)=0.0_dp

          ! Parser from .dat-file
          !call fparser_evalFunction(rfparser, 1, rdomainIntSubset%p_DcubPtsReal(:,ipoint,iel), Dcoefficients(1,ipoint,iel))

!          ! Pressure hill
!          Dcoefficients (1,ipoint,iel) = 1.0_dp + 0.3_dp*&
!                   exp(-40.0_dp*((Dpoints(1,ipoint,iel)-0.5_dp)**2.0_dp+(Dpoints(2,ipoint,iel)-0.5_dp)**2.0_dp))

!          ! Circular dambreak (for Euler equations)
!          r = sqrt(2.0_dp*(dx-0.5_dp)**2.0_dp+(dy-0.5_dp)**2.0_dp)
!          if ( r<0.1_dp) then
!            Dcoefficients (1,ipoint,iel) = 1.0_dp + 0.15_dp
!          else
!            Dcoefficients (1,ipoint,iel) = 1.0_dp
!          end if
          
          ! Shock tube
          if (dx<0.5_dp) then
            Dcoefficients (1,ipoint,iel) = 1.0_dp
          else
            Dcoefficients (1,ipoint,iel) = 0.125_dp
          end if 
          
          
          ! For scalar problem
          Dcoefficients (1,ipoint,iel) = 0.5_dp


!        ! Isentropicvortex
!        dt = 0
!        drad = sqrt((dx-dt-5.0_dp)**2.0_dp + dy*dy)
!        
!        Dcoefficients (1,ipoint,iel) = (1.0_dp-0.4_dp/(16.0_dp*1.4_dp*SYS_pi**2.0_dp)*25.0_dp*exp(2.0_dp*(1.0_dp-drad*drad)))**(1.0_dp/0.4_dp)




!          ! Scalar convection
!          r = sqrt((dx-0.5_dp)**2.0_dp+(dy-0.5_dp)**2.0_dp)
!          Dcoefficients (1,ipoint,iel) = exp(-100.0_dp*r*r)

!            ! Smooth pressure hill
!            r = sqrt(2.0_dp*(dx-0.5_dp)**2.0_dp+(dy-0.5_dp)**2.0_dp)
!            Dcoefficients (1,ipoint,iel) = 1.0_dp+exp(-r)


!          Dcoefficients (1,ipoint,iel) = 1.0_dp + dx
           
          end do
        end do
      
      case (2) ! Set x-momentum
      
      Dcoefficients (1,:,:) = 0.0_dp
      
      do iel = 1, size(Dcoefficients,3)
          do ipoint = 1, size(Dcoefficients,2)

          ! Get coordinates          
          dx = Dpoints(1,ipoint,iel)
          dy = Dpoints(2,ipoint,iel)
          
          
!       ! Isentropicvortex
!        dt = 0
!        drad = sqrt((dx-dt-5.0_dp)**2.0_dp + dy*dy)
!        drho = (1.0_dp-0.4_dp/(16.0_dp*1.4_dp*SYS_pi**2.0_dp)*25.0_dp*exp(2.0_dp*(1.0_dp-drad*drad)))**(1.0_dp/0.4_dp)
!        Dcoefficients (1,ipoint,iel) = drho*(1.0_dp-5.0_dp*exp(1.0_dp-drad*drad)*(dy)/(2.0_dp*SYS_pi))

      
          end do
        end do
      case (3) ! Set y-momentum
      
        Dcoefficients (1,:,:) = 0.0_dp
        
        do iel = 1, size(Dcoefficients,3)
          do ipoint = 1, size(Dcoefficients,2)

          ! Get coordinates          
          dx = Dpoints(1,ipoint,iel)
          dy = Dpoints(2,ipoint,iel)
          
          
!       ! Isentropicvortex
!        dt = 0
!        drad = sqrt((dx-dt-5.0_dp)**2.0_dp + dy*dy)
!        drho = (1.0_dp-0.4_dp/(16.0_dp*1.4_dp*SYS_pi**2.0_dp)*25.0_dp*exp(2.0_dp*(1.0_dp-drad*drad)))**(1.0_dp/0.4_dp)
!        Dcoefficients (1,ipoint,iel) = drho*(5.0_dp*exp(1.0_dp-drad*drad)*(dx-5.0_dp)/(2.0_dp*SYS_pi))
        
      
          end do
        end do
        
      case (4) ! Set energy
      
        Dcoefficients (1,:,:) = 1.0_dp
        
        do iel = 1, size(Dcoefficients,3)
          do ipoint = 1, size(Dcoefficients,2)

            ! Get coordinates          
            dx = Dpoints(1,ipoint,iel)
            dy = Dpoints(2,ipoint,iel)
        
            ! Shock tube
            if (dx<0.5_dp) then
              Dcoefficients (1,ipoint,iel) = 2.5_dp
            else
              Dcoefficients (1,ipoint,iel) = 0.25_dp
            end if
            
            
            
!            ! Isentropic vortex
!            dt = 0.0_dp
!            drad = sqrt((dx-dt-5.0_dp)**2.0_dp + dy*dy)
!            drho = (1.0_dp-0.4_dp/(16.0_dp*1.4_dp*SYS_pi**2.0_dp)*25.0_dp*exp(2.0_dp*(1.0_dp-drad*drad)))**(1.0_dp/0.4_dp)
!            du = 1.0_dp-5.0_dp*exp(1.0_dp-drad*drad)*(dy)/(2.0_dp*SYS_pi)
!            dv = 5.0_dp*exp(1.0_dp-drad*drad)*(dx-5.0_dp)/(2.0_dp*SYS_pi)
!            Dcoefficients (1,ipoint,iel) = drho*((drho**gamma)/(0.4_dp*drho)+0.5_dp*(du*du+dv*dv))
            
!            ! Circular dambreak (for Euler equations)
!            r = sqrt(2.0_dp*(dx-0.5_dp)**2.0_dp+(dy-0.5_dp)**2.0_dp)
!            if ( r<0.1_dp) then
!              Dcoefficients (1,ipoint,iel) = 2.5_dp
!            else
!              Dcoefficients (1,ipoint,iel) = 0.25_dp
!            end if


!            ! Smooth pressure hill
!            r = sqrt(2.0_dp*(dx-0.5_dp)**2.0_dp+(dy-0.5_dp)**2.0_dp)
!            Dcoefficients (1,ipoint,iel) = 1.0_dp+exp(-r)
            
            
                        
!            Dcoefficients (1,ipoint,iel) = 1.0_dp
          
          end do
        end do
        
    end select
    
    ! Release the function parser
    call fparser_release(rfparser)
                    
    !Dcoefficients (1,:,:) = 0.0_dp

  end subroutine
  
  



! ***************************************************************************

!<subroutine>

  subroutine coeff_Steadyproj (rdiscretisation,rform, &
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


  integer :: iel, ipoint

  


    
    
    do iel = 1, size(Dcoefficients,3)
      call fevl_evaluate (DER_FUNC, Dcoefficients(1,:,iel),&
                 rcollection%p_rvectorQuickAccess1%RvectorBlock(1), Dpoints(:,:,iel))
    end do
                    
       ! Dcoefficients (1,:,:) =0.0_dp

  end subroutine
  
  
  
  
  
  
  
  
  
  
!   ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine flux_one (rdiscretisation,rform, &
!                  nelements,npointsPerElement,Dpoints, &
!                  IdofsTest,rdomainIntSubset,&
!                  Dcoefficients,rcollection)
!    
!    use basicgeometry
!    use triangulation
!    use collection
!    use scalarpde
!    use domainintegration
!    
!  !<description>
!    ! This subroutine is called during the vector assembly. It has to compute
!    ! the coefficients in front of the terms of the linear form.
!    !
!    ! The routine accepts a set of elements and a set of points on these
!    ! elements (cubature points) in real coordinates.
!    ! According to the terms in the linear form, the routine has to compute
!    ! simultaneously for all these points and all the terms in the linear form
!    ! the corresponding coefficients in front of the terms.
!  !</description>
!    
!  !<input>
!    ! The discretisation structure that defines the basic shape of the
!    ! triangulation with references to the underlying triangulation,
!    ! analytic boundary boundary description etc.
!    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
!    
!    ! The linear form which is currently to be evaluated:
!    type(t_linearForm), intent(in)                              :: rform
!    
!    ! Number of elements, where the coefficients must be computed.
!    integer, intent(in)                                         :: nelements
!    
!    ! Number of points per element, where the coefficients must be computed
!    integer, intent(in)                                         :: npointsPerElement
!    
!    ! This is an array of all points on all the elements where coefficients
!    ! are needed.
!    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!    ! DIMENSION(dimension,npointsPerElement,nelements)
!    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
!
!    ! An array accepting the DOF`s on all elements trial in the trial space.
!    ! DIMENSION(#local DOF`s in test space,nelements)
!    integer, dimension(:,:), intent(in) :: IdofsTest
!
!    ! This is a t_domainIntSubset structure specifying more detailed information
!    ! about the element set that is currently being integrated.
!    ! It is usually used in more complex situations (e.g. nonlinear matrices).
!    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset
!
!    ! Optional: A collection structure to provide additional 
!    ! information to the coefficient routine. 
!    type(t_collection), intent(inout), optional      :: rcollection
!    
!  !</input>
!  
!  !<output>
!    ! A list of all coefficients in front of all terms in the linear form -
!    ! for all given points on all given elements.
!    !   DIMENSION(itermCount,npointsPerElement,nelements)
!    ! with itermCount the number of terms in the linear form.
!    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
!  !</output>
!    
!  !</subroutine>
!
!    integer :: iel, ipoint
!    real(DP) :: dvx, dvy, dx, dy
!    
!    rdomainIntSubset%ielementDistribution = 1
!    
!    call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(1), &
!                                 rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!    
!    
!!    do iel = 1, size(Dcoefficients,3)
!!      call fevl_evaluate (DER_FUNC, Dcoefficients(1,:,iel),&
!!                 rcollection%p_rvectorQuickAccess1%RvectorBlock(1), Dpoints(:,:,iel))
!!    end do
!    
!    do iel = 1, size(Dcoefficients,3)
!      do ipoint = 1, size(Dcoefficients,2)
!      
!        dx = Dpoints(1,ipoint,iel)
!        dy = Dpoints(2,ipoint,iel)
!      
!        ! Zalesak
!        dvx=0.5_DP-dy
!        dvy=dx-0.5_DP
!        
!        ! Circular convection
!        dvx = dy
!        dvy = 1.0_dp - dx
!        
!        Dcoefficients (1,ipoint,iel) = Dcoefficients (1,ipoint,iel) * dvx
!      end do
!    end do
!    
!  end subroutine
!
!
!
!
!
!   ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine flux_two (rdiscretisation,rform, &
!                  nelements,npointsPerElement,Dpoints, &
!                  IdofsTest,rdomainIntSubset,&
!                  Dcoefficients,rcollection)
!    
!    use basicgeometry
!    use triangulation
!    use collection
!    use scalarpde
!    use domainintegration
!    
!  !<description>
!    ! This subroutine is called during the vector assembly. It has to compute
!    ! the coefficients in front of the terms of the linear form.
!    !
!    ! The routine accepts a set of elements and a set of points on these
!    ! elements (cubature points) in real coordinates.
!    ! According to the terms in the linear form, the routine has to compute
!    ! simultaneously for all these points and all the terms in the linear form
!    ! the corresponding coefficients in front of the terms.
!  !</description>
!    
!  !<input>
!    ! The discretisation structure that defines the basic shape of the
!    ! triangulation with references to the underlying triangulation,
!    ! analytic boundary boundary description etc.
!    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
!    
!    ! The linear form which is currently to be evaluated:
!    type(t_linearForm), intent(in)                              :: rform
!    
!    ! Number of elements, where the coefficients must be computed.
!    integer, intent(in)                                         :: nelements
!    
!    ! Number of points per element, where the coefficients must be computed
!    integer, intent(in)                                         :: npointsPerElement
!    
!    ! This is an array of all points on all the elements where coefficients
!    ! are needed.
!    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!    ! DIMENSION(dimension,npointsPerElement,nelements)
!    real(DP), dimension(:,:,:), intent(in)  :: Dpoints
!
!    ! An array accepting the DOF`s on all elements trial in the trial space.
!    ! DIMENSION(#local DOF`s in test space,nelements)
!    integer, dimension(:,:), intent(in) :: IdofsTest
!
!    ! This is a t_domainIntSubset structure specifying more detailed information
!    ! about the element set that is currently being integrated.
!    ! It is usually used in more complex situations (e.g. nonlinear matrices).
!    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset
!
!    ! Optional: A collection structure to provide additional 
!    ! information to the coefficient routine. 
!    type(t_collection), intent(inout), optional      :: rcollection
!    
!  !</input>
!  
!  !<output>
!    ! A list of all coefficients in front of all terms in the linear form -
!    ! for all given points on all given elements.
!    !   DIMENSION(itermCount,npointsPerElement,nelements)
!    ! with itermCount the number of terms in the linear form.
!    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
!  !</output>
!    
!  !</subroutine>
!
!    integer :: iel, ipoint
!    real(DP) :: dvx, dvy, dx, dy
!    
!    rdomainIntSubset%ielementDistribution = 1
!    
!    call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(1), &
!                                 rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
!    
!!    do iel = 1, size(Dcoefficients,3)
!!      call fevl_evaluate (DER_FUNC, Dcoefficients(1,:,iel),&
!!                 rcollection%p_rvectorQuickAccess1%RvectorBlock(1), Dpoints(:,:,iel))
!!    end do
!    
!    do iel = 1, size(Dcoefficients,3)
!      do ipoint = 1, size(Dcoefficients,2)
!      
!        dx = Dpoints(1,ipoint,iel)
!        dy = Dpoints(2,ipoint,iel)
!        
!        ! Zalesak
!        dvx=0.5_DP-dy
!        dvy=dx-0.5_DP
!        
!        ! Circular convection
!        dvx = dy
!        dvy = 1.0_dp - dx
!        
!        Dcoefficients (1,ipoint,iel) = Dcoefficients (1,ipoint,iel) * dvy
!      end do
!    end do
!    
!  end subroutine
  
  
  
  
  ! ***************************************************************************

!<subroutine>

  subroutine flux (rdiscretisation,rform, &
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
    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset

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

    integer :: iel, ipoint
    real(DP) :: dvx, dvy, dx, dy, dsol
    
    rdomainIntSubset%ielementDistribution = 1
    
    ! rform%%itermCount gives the number of additional terms, here 2
    
    ! First evaluate the solution in each point
    call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(1), &
                                 rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
                                 
    ! For the second additive term the solution doesnt change
    !Dcoefficients(2,:,:) = Dcoefficients(1,:,:)
    
    do iel = 1, size(Dcoefficients,3)
      do ipoint = 1, size(Dcoefficients,2)
      
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)
        
!        ! Zalesak
!        dvx=0.5_DP-dy
!        dvy=dx-0.5_DP

!        ! Constant velocity
!        dvx=1.0_dp
!        dvy=1.0_dp
        
        ! Circular convection
        dvx = dy
        dvy = 1.0_dp - dx

        dsol = Dcoefficients (1,ipoint,iel)
        
        Dcoefficients (1,ipoint,iel) = dsol * dvx
        !Dcoefficients (2,ipoint,iel) = Dcoefficients (2,ipoint,iel) * dvy
        Dcoefficients (2,ipoint,iel) = dsol * dvy
      end do
    end do
    
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    !<subroutine>

    subroutine flux_dg_buildVectorBlEdge2D_sim (&
!              Dcoefficients,&
!              DsolVals,&
              DfluxValues,&
              rvectorSolBlock,&
              IelementList,&
              normal,&
!              DpointsReal,&
              rintSubSet,&
              rcollection )
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
!  real(DP), dimension(:,:,:), intent(inout) :: DsolVals
  real(DP), dimension(:,:), intent(in) :: normal
!  real(DP), dimension(:,:,:), intent(in) :: DpointsReal
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  type(t_vectorBlock), intent(in) :: rvectorSolBlock
  integer, dimension(:) , intent(in) :: IelementList
  type(t_collection), intent(inout), target, optional :: rcollection
    
  !</input>
  
  !<output>
  ! DfluxValues(nvar,ialbet,ncubp,NEL)
  real(DP), dimension(:,:,:,:), intent(out) :: DfluxValues
  !</output>
    
  !</subroutine>
  
  
  real(dp), dimension(3) :: DQi, DQa, DQroe, DFi, DFa, DFlux
  real(dp), dimension(4) :: DQroec
  integer :: ivar, iel, ipoint
  ! Dsolutionvalues(2 sides, ncubp, NEL, nvar)
  real(dp), dimension(:,:,:,:), allocatable :: Dsolutionvalues
  real(dp) :: dx, dy
  real(dp) :: dh, du, dv, dnormalPart, dtangentialPart
  real(dp), dimension(3) :: DF1i, DF1a, DF2i, DF2a, DFx, DFy
  real(dp), dimension(3,3) :: DL, DR, DaLambda
  real(dp) :: dmaxEV
  
  real(dp),dimension(3,3) :: C
  
  real(dp), dimension(2) :: Dvel
  real(dp) :: dvn, drad
  
  
  
  
  
  


  ! Get solution values (or its derivatives)
  ! DfluxValues(nvar,ialbet,ncubp,NEL)
  ! Dsolutionvalues(2 sides, ncubp, NEL, nvar)
  allocate(Dsolutionvalues(2,ubound(DfluxValues,3),ubound(DfluxValues,4),ubound(DfluxValues,1)))
  
  do ivar = 1, size(DfluxValues,1)
  
    ! Get values on the one side of the edge
    call fevl_evaluate_sim4 (rvectorSolBlock%RvectorBlock(ivar), &
                             rIntSubset(1), DER_FUNC, Dsolutionvalues(:,:,:,ivar), 1)
    ! Get values on the other side of the edge                               
    call fevl_evaluate_sim4 (rvectorSolBlock%RvectorBlock(ivar), &
                             rIntSubset(2), DER_FUNC, Dsolutionvalues(:,:,:,ivar), 2)
   end do
               



  
  
  do iel = 1, ubound(DfluxValues,4)
  
    do ipoint= 1, ubound(DfluxValues,3)
      
      dx = rintSubset(1)%p_DcubPtsReal(1,ipoint,iel)
      dy = rintSubset(1)%p_DcubPtsReal(2,ipoint,iel)
      
    
      ! Set boundary conditions
      ! Test, if we are at a boundary
      if (IelementList(iel)==0) then
      
        ! Riemann BC
        !Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,:) = (/1.0_dp,0.0_dp,0.0_dp/)
        !if ((dx<0.00001).or.(dy<0.00001)) Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,:) = (/1.1_dp,0.0_dp,0.0_dp/)
        
        ! No BCs
        Dsolutionvalues(2,ipoint,iel,1) = Dsolutionvalues(1,ipoint,iel,1)
        Dsolutionvalues(2,ipoint,iel,2) = Dsolutionvalues(1,ipoint,iel,2)
        Dsolutionvalues(2,ipoint,iel,3) = Dsolutionvalues(1,ipoint,iel,3)

!        ! BC for Source term 1
!        Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,1) = (1.0_dp+dx+dy)
!        Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,2) = (1.0_dp+dx+dy)*dx
!        Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,3) = (1.0_dp+dx+dy)*dy
        
!        ! BC for Source term 2
!        Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,1) = (1.0_dp+dx**4.0_dp+dy**4.0_dp)
!        Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,2) = (1.0_dp+dx**4.0_dp+dy**4.0_dp)*dx
!        Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,3) = (1.0_dp+dx**4.0_dp+dy**4.0_dp)*dy
        
!        ! Reflecting BCs
!        ! Calculate x- and y- velocity from momentum
!        dh = Dsolutionvalues(1,ipoint,iel,1)
!        du = Dsolutionvalues(1,ipoint,iel,2)/dh
!        dv = Dsolutionvalues(1,ipoint,iel,3)/dh
!        
!        ! Calculate normal and tangential part
!        dnormalPart     = du*normal(1,iel) + dv*normal(2,iel)
!        dtangentialPart = du*normal(2,iel) - dv*normal(1,iel)
!        
!        ! Invert the normal part
!        dnormalPart     = -dnormalPart
!        dtangentialPart = +dtangentialPart
!        
!        
!!        if ((dx>0.1_DP).and.(dx<30.0_DP).and.(dy>0.1_DP).and.(dy<9.90_DP)) then
!!          write(*,*) du,dv
!!        end if
!        
!        ! Calculate new velocity
!        du = dnormalPart*normal(1,iel) + dtangentialPart*normal(2,iel)
!        dv = dnormalPart*normal(2,iel) - dtangentialPart*normal(1,iel)
!        
!!        if ((dx>0.1_DP).and.(dx<30.0_DP).and.(dy>0.1_DP).and.(dy<9.90_DP)) then
!!          write(*,*) du,dv
!!          write(*,*) ''
!!        end if
!        
!        ! Set new momentum
!        Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,1) = dh
!        Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,2) = dh * du
!        Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,3) = dh * dv
        
      end if
      
    
!      ! *** Upwind flux ***
!      
!      ! Get solution values on the in and outside
!      DQi = Dsolutionvalues(1,ipoint,iel,:)
!      DQa = Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,:)
!      
!      ! Get fluxes on the in and outside in x- and y-direction
!      DF1i = buildFlux(DQi,1)
!      DF1a = buildFlux(DQa,1)
!      DF2i = buildFlux(DQi,2)
!      DF2a = buildFlux(DQa,2)
!            
!      ! Calculate Roevalues
!      DQroe = calculateQroe(DQi,DQa)
!      
!      ! First calculate flux in x-direction
!      DFx= 0.5_dp*(DF1i+DF1a -& ! centered part
!                   normal(1,iel)*matmul(buildTrafo(DQroe,1),matmul(buildaLambda(DQroe,1),matmul(buildinvTrafo(DQroe,1),(dQa-dQi))))) ! artificial diffusion
!                   !sign(1.0_dp,normal(1,iel))*matmul(buildTrafo(DQroe,1),matmul(buildaLambda(DQroe,1), matmul(buildinvTrafo(DQroe,1),(dQa-dQi)))))! artificial diffusion
!
!                   
!      !DFx= 0.5_dp*(DF1i+DF1a - sign(1.0_dp,normal(1,iel))* maxval(abs(buildEigenvalues(DQroe,1)))*(dQa-dQi))
!      
!      ! First calculate flux in y-direction
!      DFy= 0.5_dp*(DF2i+DF2a -& ! centered part
!                   normal(2,iel)*matmul(buildTrafo(DQroe,2),matmul(buildaLambda(DQroe,2),matmul(buildinvTrafo(DQroe,2),(dQa-dQi))))) ! artificial diffusion
!                   !sign(1.0_dp,normal(2,iel))*matmul(buildTrafo(DQroe,2),matmul(buildaLambda(DQroe,2),matmul(buildinvTrafo(DQroe,2),(dQa-dQi))))) ! artificial diffusion
!      
!      !DFy= 0.5_dp*(DF2i+DF2a - sign(1.0_dp,normal(2,iel))* maxval(abs(buildEigenvalues(DQroe,2)))*(dQa-dQi))
!                   
!      ! Add the fluxes of the two dimensional directions to get Flux * normal
!      DFlux = DFx*normal(1,iel) + DFy*normal(2,iel)
!      
!      ! Save the calculated flux
!      DfluxValues(:,1,ipoint,iel) = DFlux
      
!      ! *** centered flux ***
!      ! Get solution values on the in and outside
!      DQi = Dsolutionvalues(1,ipoint,iel,:)
!      DQa = Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,:)
!      
!      ! Get fluxes on the in and outside in x- and y-direction
!      DF1i = buildFlux(DQi,1)
!      DF1a = buildFlux(DQa,1)
!      DF2i = buildFlux(DQi,2)
!      DF2a = buildFlux(DQa,2)
!            
!      DFx= 0.5_dp*(DF1i+DF1a)
!      DFy= 0.5_dp*(DF2i+DF2a)
!      
!      ! Add the fluxes of the two dimensional directions to get Flux * normal
!      DFlux = DFx*normal(1,iel) + DFy*normal(2,iel)
!      
!      ! Save the calculated flux
!      DfluxValues(:,1,ipoint,iel) = DFlux

!      if (IelementList(iel)==0) then
!      do ivar=1,3
!      DfluxValues(ivar,1,ipoint,iel) = 0.0_dp
!      end do
!      end if

!      ! *** Local Lax-Friedrichs flux ***
!      DQi = Dsolutionvalues(1,ipoint,iel,:)
!      DQa = Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,:)
!      
!      ! Get fluxes on the in and outside in x- and y-direction
!      DF1i = buildFlux(DQi,1)
!      DF1a = buildFlux(DQa,1)
!      DF2i = buildFlux(DQi,2)
!      DF2a = buildFlux(DQa,2)
!            
!      ! Calculate Roevalues
!      DQroe = calculateQroe(DQi,DQa)
!      
!      ! First calculate flux in x-direction
!      DFx= 0.5_dp*(DF1i+DF1a)
!      
!      ! First calculate flux in y-direction
!      DFy= 0.5_dp*(DF2i+DF2a)
!                   
!      ! Add the fluxes of the two dimensional directions to get Flux * normal
!      DFlux = 0.5_dp*(DFx*normal(1,iel) + DFy*normal(2,iel))
!      
!      ! Get estimate for biggest eigenvalue
!      dmaxEV = maxval(buildEigenvalues2(DQRoe,normal(1,iel),normal(2,iel)))
!      DFlux = DFlux - dmaxEV*(DQa - DQi)
!      
!      ! Save the calculated flux
!      DfluxValues(:,1,ipoint,iel) = DFlux

!      ! *** Upwind flux without dimensional splitting (Roe-Flux) ***
!      
!      ! Get solution values on the in and outside
!      DQi = Dsolutionvalues(1,ipoint,iel,:)
!      DQa = Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,:)
!      
!      ! Get fluxes on the in and outside in x- and y-direction
!      DF1i = buildFlux(DQi,1)
!      DF1a = buildFlux(DQa,1)
!      DF2i = buildFlux(DQi,2)
!      DF2a = buildFlux(DQa,2)
!            
!      ! Calculate Roevalues
!      DQroe = calculateQroe(DQi,DQa)
!      
!      ! First calculate flux in x-direction
!      DFx= 0.5_dp*(DF1i+DF1a)
!!      DFx = buildFlux(DQRoe,1)
!      
!      ! First calculate flux in y-direction
!      DFy= 0.5_dp*(DF2i+DF2a)
!!      DFy = buildFlux(DQRoe,2)
!      
!      ! Add the fluxes of the two dimensional directions to get Flux * normal
!      DFlux = DFx*normal(1,iel) + DFy*normal(2,iel)
!      
!      ! Add artificial diffusion
!      DL       = buildMixedL2       (DQRoe,normal(1,iel),normal(2,iel))
!      DaLambda = buildMixedaLambda2 (DQRoe,normal(1,iel),normal(2,iel))
!      DR       = buildMixedR2       (DQRoe,normal(1,iel),normal(2,iel))
!      
!      ! Save the calculated flux
!      DfluxValues(:,1,ipoint,iel) = DFlux + 0.5_dp*matmul(DR,matmul(DaLambda,matmul(DL,DQi - DQa)))



      ! *** Upwind flux without dimensional splitting (Roe-Flux) and new c ***
      
      ! Get solution values on the in and outside
      DQi = Dsolutionvalues(1,ipoint,iel,:)
      DQa = Dsolutionvalues(2,ipoint,iel,:)
      
      ! Get fluxes on the in and outside in x- and y-direction
      DF1i = buildFlux(DQi,1)
      DF1a = buildFlux(DQa,1)
      DF2i = buildFlux(DQi,2)
      DF2a = buildFlux(DQa,2)
            
      ! Calculate Roevalues
      DQroec = calculateQroec(DQi,DQa)
      
      ! First calculate flux in x-direction
      DFx= 0.5_dp*(DF1i+DF1a)
!      DFx = buildFlux(DQRoe,1)
      
      ! First calculate flux in y-direction
      DFy= 0.5_dp*(DF2i+DF2a)
!      DFy = buildFlux(DQRoe,2)
      
      ! Add the fluxes of the two dimensional directions to get Flux * normal
      DFlux = DFx*normal(1,iel) + DFy*normal(2,iel)
      
      ! Add artificial diffusion
      DL       = buildMixedL2c       (DQRoec,normal(1,iel),normal(2,iel))
      DaLambda = buildMixedaLambda2c (DQRoec,normal(1,iel),normal(2,iel))
      DR       = buildMixedR2c       (DQRoec,normal(1,iel),normal(2,iel))
      
      ! Save the calculated flux
      DfluxValues(:,1,ipoint,iel) = DFlux + 0.5_dp*matmul(DR,matmul(DaLambda,matmul(DL,DQi - DQa)))
      
      
      
!      ! *** Upwind flux for quarter circular convection ***
!    
!      ! BC for steady circular convection
!      drad = abs( sqrt((dx-1.0_dp)*(dx-1.0_dp)+(dy-0.0_dp)*(dy-0.0_dp)) -0.5_dp)
!      if (IelementList(iel)==0) Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,1) = 0.3_DP*exp(-25.0_dp*drad*drad)
!      
!      ! Get solution values on the in and outside
!      DQi(1) = Dsolutionvalues(1,ipoint,iel,1)
!      DQa(1) = Dsolutionvalues(2,ubound(DfluxValues,3)-ipoint+1,iel,1)
!      
!!      ! Zalesak
!!      Dvel(1)=0.5_DP-dy
!!      Dvel(2)=dx-0.5_DP
!      
!!      ! Constant velocity
!!      Dvel(1)=1.0_dp
!!      Dvel(2)=1.0_dp
!      
!      ! Steady circular convection
!      Dvel(1)=dy
!      Dvel(2)=1.0_DP-dx
!    
!      dvn = Dvel(1)*normal(1,iel)+Dvel(2)*normal(2,iel)
!      
!      ! Upwind flux
!      if (dvn.ge.0) then
!        DfluxValues(1,1,ipoint,iel) = dvn *DQi(1)
!      else
!        DfluxValues(1,1,ipoint,iel) = dvn *DQa(1)
!      end if
      
      
      
      
    end do ! ipoint
  end do ! iel
  

  
  

                                 
  deallocate(Dsolutionvalues)  
  
  
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
   ! ***************************************************************************

!<subroutine>

  subroutine flux_sys (rdiscretisation,rform, &
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
    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset

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

    integer :: iel, ipoint, ivar, nvar2d, currentvar,nterms
    real(DP) :: dvx, dvy, dx, dy, dsol
    ! (# of terms, ncubp, NEL, nvar2d)
    real(dp), dimension(:,:,:,:), allocatable :: DsolutionValues
    real(dp), dimension(3) :: dQ, DFx, DFy
    real(dp), dimension(:), pointer :: p_ddata
    
    !nvar2d = 1
    nvar2d = rcollection%p_rvectorQuickAccess1%Nblocks
    
    rdomainIntSubset%ielementDistribution = 1
    
    ! Which variable (of the system) we are at
    currentvar = rcollection%IquickAccess(1)
    
    ! rform%itermCount gives the number of additional terms, here 2
    nterms = rform%itermCount
    
    ! Allocate space for the solution values ! (# of terms, ncubp, NEL, nvar2d)
    allocate(DsolutionValues(nterms,size(Dpoints,2),size(Dpoints,3),nvar2d))
    
    ! First evaluate the solution in each point
    do ivar = 1, nvar2d
      call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(ivar), &
                                 rdomainIntSubset, DER_FUNC, DsolutionValues(:,:,:,ivar), 1)
    end do  
    
!call    lsyssc_getbase_double(rcollection%p_rvectorQuickAccess1%RvectorBlock(3),p_ddata)
!    write(*,*) p_ddata
                               
                                 

    
    do iel = 1, size(Dcoefficients,3)
      do ipoint = 1, size(Dcoefficients,2)
      
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)
        
        
        ! Shallow water        
        dQ = DsolutionValues(1,ipoint,iel,:)
        !write(*,*) dQ
        
        DFx = buildFlux(dQ,1)
        DFy = buildFlux(dQ,2)
        
        Dcoefficients (1,ipoint,iel) = DFx(currentvar)
        Dcoefficients (2,ipoint,iel) = DFy(currentvar)


!!        ! Steady circular convection
!!        ! Zalesak
!!        dvx=0.5_DP-dy
!!        dvy=dx-0.5_DP
!
!!        ! Constant velocity
!!        dvx=1.0_dp
!!        dvy=1.0_dp
!        
!        ! Circular convection
!        dvx = dy
!        dvy = 1.0_dp - dx
!
!        
!        dQ(1) = DsolutionValues(1,ipoint,iel,1)
!        
!        Dcoefficients (1,ipoint,iel) = dQ(1) * dvx
!        Dcoefficients (2,ipoint,iel) = dQ(1) * dvy
        
      end do
    end do
    
    deallocate(DsolutionValues)
    
  end subroutine
  



! ***************************************************************************

!<subroutine>

  subroutine flux_sys_block (rdiscretisation,rform, &
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
    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(nvar,itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    integer :: iel, ipoint, ivar, nvar2d, currentvar,nterms
    real(DP) :: dvx, dvy, dx, dy, dsol
    ! (# of terms, ncubp, NEL, nvar2d)
    real(dp), dimension(:,:,:,:), allocatable :: DsolutionValues
    real(dp), dimension(3) :: dQ, DFx, DFy
    real(dp), dimension(:), pointer :: p_ddata
    
    ! nvar2d = 3
    nvar2d = rcollection%p_rvectorQuickAccess1%Nblocks
    
    rdomainIntSubset%ielementDistribution = 1
    
    
    ! rform%itermCount gives the number of additional terms, here 2
    nterms = rform%itermCount
        
    ! Allocate space for the solution values ! (# of derivatives, ncubp, NEL, nvar2d)
    allocate(DsolutionValues(1,size(Dpoints,2),size(Dpoints,3),nvar2d))
    
    ! First evaluate the solution in each point
    do ivar = 1, nvar2d
      call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(ivar), &
                                 rdomainIntSubset, DER_FUNC, DsolutionValues(:,:,:,ivar), 1)
    end do  
 
 
    
    do iel = 1, size(Dpoints,3)
      do ipoint = 1, size(Dpoints,2)
      
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)

        ! Shallow Water        
        dQ = DsolutionValues(1,ipoint,iel,:)
        !write(*,*) dQ
        
        DFx = buildFlux(dQ,1)
        DFy = buildFlux(dQ,2)
        
        Dcoefficients (:,1,ipoint,iel) = DFx(:)
        Dcoefficients (:,2,ipoint,iel) = DFy(:)

       
       
!       ! Circular convection
!        dvx = dy
!        dvy = 1.0_dp - dx
!
!        
!        dQ(1) = DsolutionValues(1,ipoint,iel,1)
!        
!        Dcoefficients (1,1,ipoint,iel) = dQ(1) * dvx
!        Dcoefficients (1,2,ipoint,iel) = dQ(1) * dvy
        
      end do
    end do
    
    deallocate(DsolutionValues)
    
  end subroutine
    
    
    
    
    
! ***************************************************************************

!<subroutine>

  subroutine Callback_Sourceterm (rdiscretisation,rform, &
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
    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(nvar,itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    integer :: iel, ipoint, ivar, nvar2d, currentvar,nterms
    real(DP) :: dvx, dvy, dx, dy, dsol, dh, du, dv
    ! (# of terms, ncubp, NEL, nvar2d)
    real(dp), dimension(:,:,:,:), allocatable :: DsolutionValues
    real(dp), dimension(3) :: dQ, DFx, DFy
    real(dp), dimension(:), pointer :: p_ddata
    
    !nvar2d = 3
    nvar2d = rcollection%p_rvectorQuickAccess1%Nblocks
    
    rdomainIntSubset%ielementDistribution = 1
    
    
    ! rform%itermCount gives the number of additional terms, here 2
    nterms = rform%itermCount
        
    
    do iel = 1, size(Dpoints,3)
      do ipoint = 1, size(Dpoints,2)
      
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)
        
       
!        ! Source Term 1
!        Dcoefficients (1,1,ipoint,iel) = 2.0_dp+3.0_dp*dx+3.0_dp*dy
!        Dcoefficients (2,1,ipoint,iel) = 4.0_dp*dx*dx+4.0_dp*dx+4.0_dp*dx*dy+dy+1.0_dp
!        Dcoefficients (3,1,ipoint,iel) = 4.0_dp*dy*dy+4.0_dp*dy+4.0_dp*dx*dy+dx+1.0_dp
        
        ! Source term 2
        Dcoefficients (1,1,ipoint,iel) = 6.0_dp*(dx**4.0_dp+dy**4.0_dp)+2.0_dp
        Dcoefficients (2,1,ipoint,iel) = 5.0_dp*dy**4.0_dp*dx+dx+dx**5.0_dp+6.0_dp*dx**5.0_dp+2.0_dp*dx+2.0_dp*dy**4.0_dp*dx+4.0_dp*dx**3.0_dp+4.0_dp*dx**7.0_dp+4.0_dp*dx**3.0_dp*dy**4.0_dp
        Dcoefficients (3,1,ipoint,iel) = 5.0_dp*dx**4.0_dp*dy+dy+dy**5.0_dp+6.0_dp*dy**5.0_dp+2.0_dp*dy+2.0_dp*dx**4.0_dp*dy+4.0_dp*dy**3.0_dp+4.0_dp*dy**7.0_dp+4.0_dp*dy**3.0_dp*dx**4.0_dp
       
        
      end do
    end do
    
  end subroutine
  
  
  
  
  
  
  
  
  
  
  !<subroutine>

    subroutine Euler_flux_dg_buildVectorBlEdge2D_sim (&
!              Dcoefficients,&
!              DsolVals,&
              DfluxValues,&
              rvectorSolBlock,&
              IelementList,&
              normal,&
!              DpointsReal,&
              rintSubSet,&
              rcollection )
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
!  real(DP), dimension(:,:,:), intent(inout) :: DsolVals
  real(DP), dimension(:,:), intent(in) :: normal
!  real(DP), dimension(:,:,:), intent(in) :: DpointsReal
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  type(t_vectorBlock), intent(in) :: rvectorSolBlock
  integer, dimension(:) , intent(in) :: IelementList
  type(t_collection), intent(inout), target, optional :: rcollection
    
  !</input>
  
  !<output>
  ! DfluxValues(nvar,ialbet,ncubp,NEL)
  real(DP), dimension(:,:,:,:), intent(out) :: DfluxValues
  !</output>
    
  !</subroutine>
  
  
  real(dp), dimension(4) :: DQi, DQa, DQroe, DFi, DFa, DFlux
  real(dp), dimension(5) :: DQroec
  integer :: ivar, iedge, ipoint
  ! Dsolutionvalues(2 sides, ncubp, NEL, nvar)
  real(dp), dimension(:,:,:,:), allocatable :: Dsolutionvalues
  real(dp) :: dx, dy
  real(dp) :: dh, du, dv, dnormalPart, dtangentialPart
  real(dp), dimension(4) :: DFxi, DFxa, DFyi, DFya, DFx, DFy
  real(dp), dimension(4,4) :: DL, DR, DaLambda
  real(dp) :: dmaxEV
  
  real(dp),dimension(4,4) :: C
  
  real(dp), dimension(2) :: Dvel
  real(dp) :: dvn, drad, dvt
  real(dp) :: drho, dE, dpr, dc, dt
  real(dp) :: dW1, dW2, dW3, dW4, dW1o, dW2o, dW3o, dW4o
  real(dp) :: dlambda
  real(dp) :: gamma = 1.4_dp
  real(dp) :: El, Er, pl, pr
  
  real(dp), dimension(4) :: dtemporal
  
  

  ! Get solution values (or its derivatives)
  ! DfluxValues(nvar,ialbet,ncubp,NEL)
  ! Dsolutionvalues(2 sides, ncubp, NEL, nvar)
  allocate(Dsolutionvalues(2,ubound(DfluxValues,3),ubound(DfluxValues,4),ubound(DfluxValues,1)))
  
  do ivar = 1, size(DfluxValues,1)
  
    ! Get values on the one side of the edge
    call fevl_evaluate_sim4 (rvectorSolBlock%RvectorBlock(ivar), &
                             rIntSubset(1), DER_FUNC, Dsolutionvalues(:,:,:,ivar), 1)
    ! Get values on the other side of the edge                               
    call fevl_evaluate_sim4 (rvectorSolBlock%RvectorBlock(ivar), &
                             rIntSubset(2), DER_FUNC, Dsolutionvalues(:,:,:,ivar), 2)
   end do
               



  
  
  do iedge = 1, ubound(DfluxValues,4)
    
    do ipoint= 1, ubound(DfluxValues,3)
      
      ! Get coordinates
      dx = rintSubset(1)%p_DcubPtsReal(1,ipoint,iedge)
      dy = rintSubset(1)%p_DcubPtsReal(2,ipoint,iedge)

    
    
      ! Set boundary conditions
      ! Test, if we are at a boundary
      if (IelementList(iedge)==0) then
      
        ! No BCs
        Dsolutionvalues(2,ipoint,iedge,1) = Dsolutionvalues(1,ipoint,iedge,1)
        Dsolutionvalues(2,ipoint,iedge,2) = Dsolutionvalues(1,ipoint,iedge,2)
        Dsolutionvalues(2,ipoint,iedge,3) = Dsolutionvalues(1,ipoint,iedge,3)
        Dsolutionvalues(2,ipoint,iedge,4) = Dsolutionvalues(1,ipoint,iedge,4)


!        !!! Boundary conditions by Riemann invariants !!!
!        
!        ! Calculate Riemann invariants from outer (freestream) state
!        ! Here you can set the desired freestream state
!        drho = 1.0_dp
!        du = 0.0_dp
!        dv = 0.0_dp
!        dE = 0.0_dp
!        
!        dpr = (gamma-1)*drho*(dE-0.5_dp*( du*du + dv*dv ) )
!        dH = dE + dpr/drho
!        dc = sqrt((gamma-1)*(dH-0.5_dp*( du*du + dv*dv ) ))
!        
!        dvn =  du*normal(1,iedge) + dv*normal(2,iedge)
!        dvt = -du*normal(2,iedge) + dv*normal(1,iedge)
!        
!        dW1o = dvn - 2.0_dp*dc/(gamma-1)
!        dW2o = dpr/(drho**gamma)
!        dW3o = dvt
!        dW4o = dvn + 2.0_dp*dc/(gamma-1)        
!        
!        ! Calculate Riemann invariants from inner values
!        drho = Dsolutionvalues(1,ipoint,iedge,1)
!        du = Dsolutionvalues(1,ipoint,iedge,2)/drho
!        dv = Dsolutionvalues(1,ipoint,iedge,3)/drho
!        dE = Dsolutionvalues(1,ipoint,iedge,4)/drho
!        
!        dpr = (gamma-1)*drho*(dE-0.5_dp*( du*du + dv*dv ) )
!        dH = dE + dpr/drho
!        dc = sqrt((gamma-1)*(dH-0.5_dp*( du*du + dv*dv ) ))
!        
!        dvn =  du*normal(1,iedge) + dv*normal(2,iedge)
!        dvt = -du*normal(2,iedge) + dv*normal(1,iedge)
!        
!        dW1 = dvn - 2.0_dp*dc/(gamma-1)
!        dW2 = dpr/(drho**gamma)
!        dW3 = dvt
!        dW4 = dvn + 2.0_dp*dc/(gamma-1)
!        
!        ! Choose inner/outer state depending on the sign of the
!        ! eigenvalues
!        if ((dvn-dc)<0.0_dp) dW1 = dW1o
!        if (dvn<0.0_dp) then
!          dW2 = dW2o
!          dW3 = dW3o
!        end if
!        if ((dvn+dc)<0.0_dp) dW4 = dW4o
!        
!        ! Transform back to conservative variables and set value
!        ! in the ghost node
!        dc = 0.25_dp*(gamma-1.0_dp)*(dW4-dW1)
!        drho = (dc*dc/(gamma*dW2))**(1.0_dp/(gamma-1.0_dp))
!        dpr = dc*dc*drho/gamma
!        du = 0.5_dp*(dW1+dW4)*normal(1,iedge) - dW3*normal(2,iedge)
!        dv = 0.5_dp*(dW1+dW4)*normal(2,iedge) + dW3*normal(1,iedge)
!        
!        Dsolutionvalues(2,ipoint,iedge,1) = drho
!        Dsolutionvalues(2,ipoint,iedge,2) = drho*du
!        Dsolutionvalues(2,ipoint,iedge,3) = drho*dv
!        Dsolutionvalues(2,ipoint,iedge,4) = dpr/(gamma-1.0_dp) + drho*0.5_dp*(du*du+dv*dv)
        
        
        
        ! Reflecting BCs
        ! Calculate x- and y- velocity from momentum
        drho = Dsolutionvalues(1,ipoint,iedge,1)
        du = Dsolutionvalues(1,ipoint,iedge,2)/drho
        dv = Dsolutionvalues(1,ipoint,iedge,3)/drho
        
        ! Calculate normal and tangential part
        dvn =  du*normal(1,iedge) + dv*normal(2,iedge)
        dvt = -du*normal(2,iedge) + dv*normal(1,iedge)
        
        ! Invert the normal part
        dvn = -dvn
        dvt =  dvt
        
        ! Calculate new velocity
        du = dvn*normal(1,iedge) - dvt*normal(2,iedge)
        dv = dvn*normal(2,iedge) + dvt*normal(1,iedge)
        
        ! Set new momentum
        Dsolutionvalues(2,ipoint,iedge,1) = drho
        Dsolutionvalues(2,ipoint,iedge,2) = drho * du
        Dsolutionvalues(2,ipoint,iedge,3) = drho * dv
        Dsolutionvalues(2,ipoint,iedge,4) = Dsolutionvalues(1,ipoint,iedge,4)


!        !!! Boundary conditions by Riemann invariants for isentropic vortex !!!
!        
!        ! Calculate Riemann invariants from outer (freestream) state
!        ! Here you can set the desired freestream state
!        dt = rcollection%Dquickaccess(1)
!        drad = sqrt((dx-dt-5.0_dp)**2.0_dp + dy*dy)
!        
!        drho = (1.0_dp-0.4_dp/(16.0_dp*1.4_dp*SYS_pi**2.0_dp)*25.0_dp*exp(2.0_dp*(1.0_dp-drad*drad)))**(1.0_dp/0.4_dp)
!        du = 1.0_dp-5.0_dp*exp(1.0_dp-drad*drad)*(dy)/(2.0_dp*SYS_pi)
!        dv = 5.0_dp*exp(1.0_dp-drad*drad)*(dx-5.0_dp)/(2.0_dp*SYS_pi)
!        dE = (drho**gamma)/(0.4_dp*drho)+0.5_dp*(du*du+dv*dv)
!        
!        dpr = (gamma-1)*drho*(dE-0.5_dp*( du*du + dv*dv ) )
!        dH = dE + dpr/drho
!        dc = sqrt((gamma-1)*(dH-0.5_dp*( du*du + dv*dv ) ))
!        
!        dvn =  du*normal(1,iedge) + dv*normal(2,iedge)
!        dvt = -du*normal(2,iedge) + dv*normal(1,iedge)
!        
!        dW1o = dvn - 2.0_dp*dc/(gamma-1)
!        dW2o = dpr/(drho**gamma)
!        dW3o = dvt
!        dW4o = dvn + 2.0_dp*dc/(gamma-1)        
!        
!        ! Calculate Riemann invariants from inner values
!        drho = Dsolutionvalues(1,ipoint,iedge,1)
!        du = Dsolutionvalues(1,ipoint,iedge,2)/drho
!        dv = Dsolutionvalues(1,ipoint,iedge,3)/drho
!        dE = Dsolutionvalues(1,ipoint,iedge,4)/drho
!        
!        dpr = (gamma-1)*drho*(dE-0.5_dp*( du*du + dv*dv ) )
!        dH = dE + dpr/drho
!        dc = sqrt((gamma-1)*(dH-0.5_dp*( du*du + dv*dv ) ))
!        
!        dvn =  du*normal(1,iedge) + dv*normal(2,iedge)
!        dvt = -du*normal(2,iedge) + dv*normal(1,iedge)
!        
!        dW1 = dvn - 2.0_dp*dc/(gamma-1)
!        dW2 = dpr/(drho**gamma)
!        dW3 = dvt
!        dW4 = dvn + 2.0_dp*dc/(gamma-1)
!        
!        ! Choose inner/outer state depending on the sign of the
!        ! eigenvalues
!        if ((dvn-dc)<0.0_dp) dW1 = dW1o
!        if (dvn<0.0_dp) then
!          dW2 = dW2o
!          dW3 = dW3o
!        end if
!        if ((dvn+dc)<0.0_dp) dW4 = dW4o
!        
!        ! Transform back to conservative variables and set value
!        ! in the ghost node
!        dc = 0.25_dp*(gamma-1.0_dp)*(dW4-dW1)
!        drho = (dc*dc/(gamma*dW2))**(1.0_dp/(gamma-1.0_dp))
!        dpr = dc*dc*drho/gamma
!        du = 0.5_dp*(dW1+dW4)*normal(1,iedge) - dW3*normal(2,iedge)
!        dv = 0.5_dp*(dW1+dW4)*normal(2,iedge) + dW3*normal(1,iedge)
!        
!        Dsolutionvalues(2,ipoint,iedge,1) = drho
!        Dsolutionvalues(2,ipoint,iedge,2) = drho*du
!        Dsolutionvalues(2,ipoint,iedge,3) = drho*dv
!        Dsolutionvalues(2,ipoint,iedge,4) = dpr/(gamma-1.0_dp) + drho*0.5_dp*(du*du+dv*dv)
        
      end if
      

      ! *** Upwind flux without dimensional splitting (Roe-Flux) and new c ***
      
      ! Get solution values on the in and outside
      DQi = Dsolutionvalues(1,ipoint,iedge,:)
      DQa = Dsolutionvalues(2,ipoint,iedge,:)
      
!      ! Get fluxes on the in and outside in x- and y-direction
!      DFxi = Euler_buildFlux(DQi,1)
!      DFxa = Euler_buildFlux(DQa,1)
!      DFyi = Euler_buildFlux(DQi,2)
!      DFya = Euler_buildFlux(DQa,2)
!            
!      ! Calculate Roevalues
!      DQroec = Euler_calculateQroec(DQi,DQa)
!      
!      ! First calculate flux in x-direction
!      DFx = 0.5_dp*(DFxi+DFxa)
!      !DFx = Euler_buildFlux(DQRoe,1)
!      
!      ! First calculate flux in y-direction
!      DFy = 0.5_dp*(DFyi+DFya)
!      !DFy = Euler_buildFlux(DQRoe,2)
!      
!      ! Add the fluxes of the two dimensional directions to get Flux * normal
!      DFlux = DFx*normal(1,iedge) + DFy*normal(2,iedge)
!      
!      ! Add artificial diffusion
!      DL       = Euler_buildMixedLcfromRoe       (DQRoec,normal(1,iedge),normal(2,iedge))
!      DaLambda = Euler_buildMixedaLambdacfromRoe (DQRoec,normal(1,iedge),normal(2,iedge))
!      DR       = Euler_buildMixedRcfromRoe       (DQRoec,normal(1,iedge),normal(2,iedge))
!      
!      ! Save the calculated flux (Roe)
!      DfluxValues(:,1,ipoint,iedge) = DFlux + 0.5_dp*matmul(DR,matmul(DaLambda,matmul(DL,DQi - DQa)))
      
!      ! Save the calculated flux (Mean)
!      DfluxValues(:,1,ipoint,iedge) =  Euler_buildFlux(0.5_dp*(DQi+DQa),1)*normal(1,iedge) + Euler_buildFlux(0.5_dp*(DQi+DQa),2)*normal(2,iedge)
      
      
!      ! Lax Friedrichs #1
!      DfluxValues(:,1,ipoint,iedge) = DFlux + 0.5_dp*maxval(DaLambda)*(DQi - DQa)
      
      
!      ! Foreign Roe-Flux
!      DfluxValues(:,1,ipoint,iedge) = Roe(DQi, DQa, normal(1,iedge), normal(2,iedge))
      
!      ! Foreign Rotated_RHLL-Flux
!      DfluxValues(:,1,ipoint,iedge) = Rotated_RHLL(DQi, DQa, normal(1,iedge), normal(2,iedge))
      
            
!     ! Lax Friedrichs #2
!     dlambda = max(sqrt(DQi(2)*DQi(2)+dQi(3)*DQi(3))/DQi(1)+sqrt(gamma/DQi(1)*(gamma-1.0_dp)*DQi(1)*(DQi(4)/DQi(1)-0.5_dp*( (DQi(2)/DQi(1))**2.0_dp + (DQi(3)/DQi(1))**2.0_dp ) )) , sqrt(DQa(2)*DQa(2)+dQa(3)*DQa(3))/DQa(1)+sqrt(gamma/DQa(1)*(gamma-1.0_dp)*DQa(1)*(DQa(4)/DQa(1)-0.5_dp*( (DQa(2)/DQa(1))**2.0_dp + (DQa(3)/DQa(1))**2.0_dp ) )))
!     DfluxValues(:,1,ipoint,iedge) = DFlux + 0.5_dp*dlambda*(DQi - DQa)

!      ! Lax Friedrichs with Rusanov dissipation
!      El = dQi(4)/dQi(1)
!      Er = dQa(4)/dQa(1)
!      pl = (gamma-1.0_dp)*dQi(1)*(El-0.5_dp*( (dQi(2)/dQi(1))**2.0_dp + (dQi(3)/dQi(1))**2.0_dp ) )
!      pr = (gamma-1.0_dp)*dQa(1)*(Er-0.5_dp*( (dQa(2)/dQa(1))**2.0_dp + (dQa(3)/dQa(1))**2.0_dp ) )
!      dlambda = max( sqrt((DQi(2)/DQi(1)*normal(1,iedge))**2.0_dp + (DQi(3)/DQi(1)*normal(2,iedge))**2.0_dp)+sqrt(gamma/DQi(1)*pl) , sqrt((DQa(2)/DQa(1)*normal(1,iedge))**2.0_dp + (DQa(3)/DQa(1)*normal(2,iedge))**2.0_dp)+sqrt(gamma/DQa(1)*pr) )
!      DfluxValues(:,1,ipoint,iedge) = DFlux + 0.5_dp*dlambda*(DQi - DQa)

      ! Save the calculated flux (HLL)
      DfluxValues(:,1,ipoint,iedge) = Euler_buildFlux_HLL2D(DQi,DQa,normal(1,iedge),normal(2,iedge))
      
!      ! Save the calculated flux (HLLC)
!      DfluxValues(:,1,ipoint,iedge) = Euler_buildFlux_HLLC2D(DQi,DQa,normal(1,iedge),normal(2,iedge))

!      ! Galerkin
!      DfluxValues(:,1,ipoint,iedge) = DFlux
      
    end do ! ipoint
  end do ! iedge
  

  
  

                                 
  deallocate(Dsolutionvalues)  
  
  
  end subroutine
  
  
  
  ! ***************************************************************************

!<subroutine>

  subroutine Euler_flux_sys_block (rdiscretisation,rform, &
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
    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(nvar,itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    integer :: iel, ipoint, ivar, nvar2d, currentvar,nterms
    real(DP) :: dvx, dvy, dx, dy, dsol
    ! (# of terms, ncubp, NEL, nvar2d)
    real(dp), dimension(:,:,:,:), allocatable :: DsolutionValues
    real(dp), dimension(4) :: dQ, DFx, DFy
    real(dp), dimension(:), pointer :: p_ddata
    
    ! nvar2d = 4
    nvar2d = rcollection%p_rvectorQuickAccess1%Nblocks
    
    rdomainIntSubset%ielementDistribution = 1
    
    
    ! rform%itermCount gives the number of additional terms, here 2
    nterms = rform%itermCount
        
    ! Allocate space for the solution values ! (# of derivatives, ncubp, NEL, nvar2d)
    allocate(DsolutionValues(1,size(Dpoints,2),size(Dpoints,3),nvar2d))
    
    ! First evaluate the solution in each point
    do ivar = 1, nvar2d
      call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(ivar), &
                                 rdomainIntSubset, DER_FUNC, DsolutionValues(:,:,:,ivar), 1)
    end do  
 
    ! Loop over all elements and all points   
    do iel = 1, size(Dpoints,3)
      do ipoint = 1, size(Dpoints,2)
        
        ! Get coordinates
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)

        ! Get solution      
        dQ = DsolutionValues(1,ipoint,iel,:)
        
        ! Build fluxes in coordinate directions
        DFx = Euler_buildFlux(dQ,1)
        DFy = Euler_buildFlux(dQ,2)
        
        ! Set the coefficients
        Dcoefficients (:,1,ipoint,iel) = DFx(:)
        Dcoefficients (:,2,ipoint,iel) = DFy(:)
        
      end do
    end do
    
    ! Deallocate the solutionvalues    
    deallocate(DsolutionValues)
    
  end subroutine
  
  
  
  ! ***************************************************************************

!<subroutine>

  subroutine Euler_flux_sys (rdiscretisation,rform, &
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
    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset

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

    integer :: iel, ipoint, ivar, nvar2d, currentvar,nterms
    real(DP) :: dvx, dvy, dx, dy, dsol
    ! (# of terms, ncubp, NEL, nvar2d)
    real(dp), dimension(:,:,:,:), allocatable :: DsolutionValues
    real(dp), dimension(4) :: dQ, DFx, DFy
    real(dp), dimension(:), pointer :: p_ddata
    
    !nvar2d = 1
    nvar2d = rcollection%p_rvectorQuickAccess1%Nblocks
    
    rdomainIntSubset%ielementDistribution = 1
    
    ! Which variable (of the system) we are at
    currentvar = rcollection%IquickAccess(1)
    
    ! rform%itermCount gives the number of additional terms, here 2
    nterms = rform%itermCount
    
    ! Allocate space for the solution values ! (# of terms, ncubp, NEL, nvar2d)
    allocate(DsolutionValues(nterms,size(Dpoints,2),size(Dpoints,3),nvar2d))
    
    ! First evaluate the solution in each point
    do ivar = 1, nvar2d
      call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(ivar), &
                                 rdomainIntSubset, DER_FUNC, DsolutionValues(:,:,:,ivar), 1)
    end do  
    
!call    lsyssc_getbase_double(rcollection%p_rvectorQuickAccess1%RvectorBlock(3),p_ddata)
!    write(*,*) p_ddata
                               
                                 

    
    do iel = 1, size(Dcoefficients,3)
      do ipoint = 1, size(Dcoefficients,2)
      
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)
        
        
        ! Shallow water        
        dQ = DsolutionValues(1,ipoint,iel,:)
        !write(*,*) dQ
        
        DFx = Euler_buildFlux(dQ,1)
        DFy = Euler_buildFlux(dQ,2)
        
        Dcoefficients (1,ipoint,iel) = DFx(currentvar)
        Dcoefficients (2,ipoint,iel) = DFy(currentvar)



        
      end do
    end do
    
    deallocate(DsolutionValues)
    
  end subroutine



!<subroutine>

    subroutine Scalar_flux_dg_buildVectorBlEdge2D_sim (&
!              Dcoefficients,&
!              DsolVals,&
              DfluxValues,&
              rvectorSolBlock,&
              IelementList,&
              normal,&
!              DpointsReal,&
              rintSubSet,&
              rcollection )
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
!  real(DP), dimension(:,:,:), intent(inout) :: DsolVals
  real(DP), dimension(:,:), intent(in) :: normal
!  real(DP), dimension(:,:,:), intent(in) :: DpointsReal
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  type(t_vectorBlock), intent(in) :: rvectorSolBlock
  integer, dimension(:) , intent(in) :: IelementList
  type(t_collection), intent(inout), target, optional :: rcollection
    
  !</input>
  
  !<output>
  ! DfluxValues(nvar,ialbet,ncubp,NEL)
  real(DP), dimension(:,:,:,:), intent(out) :: DfluxValues
  !</output>
    
  !</subroutine>
  
  
  real(dp), dimension(1) :: DQi, DQa, DQroe, DFi, DFa, DFlux
  real(dp), dimension(5) :: DQroec
  integer :: ivar, iedge, ipoint
  ! Dsolutionvalues(2 sides, ncubp, NEL, nvar)
  real(dp), dimension(:,:,:,:), allocatable :: Dsolutionvalues
  real(dp) :: dx, dy
  real(dp) :: dh, du, dv, dnormalPart, dtangentialPart
  real(dp), dimension(4) :: DF1i, DF1a, DF2i, DF2a, DFx, DFy
  real(dp), dimension(4,4) :: DL, DR, DaLambda
  real(dp) :: dmaxEV
  
  real(dp),dimension(4,4) :: C
  
  real(dp), dimension(2) :: Dvel
  real(dp) :: dvn, drad, dvt
  real(dp) :: drho, dE, dpr, dc
  real(dp) :: dW1, dW2, dW3, dW4, dW1o, dW2o, dW3o, dW4o
  real(dp) :: dlambda
  real(dp) :: gamma = 1.4_dp
  real(dp) :: El, Er, pl, pr
  
  real(dp), dimension(4) :: dtemporal
  
  

  ! Get solution values (or its derivatives)
  ! DfluxValues(nvar,ialbet,ncubp,NEL)
  ! Dsolutionvalues(2 sides, ncubp, NEL, nvar)
  allocate(Dsolutionvalues(2,ubound(DfluxValues,3),ubound(DfluxValues,4),ubound(DfluxValues,1)))
  
  do ivar = 1, size(DfluxValues,1)
  
    ! Get values on the one side of the edge
    call fevl_evaluate_sim4 (rvectorSolBlock%RvectorBlock(ivar), &
                             rIntSubset(1), DER_FUNC, Dsolutionvalues(:,:,:,ivar), 1)
    ! Get values on the other side of the edge                               
    call fevl_evaluate_sim4 (rvectorSolBlock%RvectorBlock(ivar), &
                             rIntSubset(2), DER_FUNC, Dsolutionvalues(:,:,:,ivar), 2)
   end do
               



  
  
  do iedge = 1, ubound(DfluxValues,4)
  
    do ipoint= 1, ubound(DfluxValues,3)
      
      ! Get coordinates
      dx = rintSubset(1)%p_DcubPtsReal(1,ipoint,iedge)
      dy = rintSubset(1)%p_DcubPtsReal(2,ipoint,iedge)
      
    
      ! Set boundary conditions
      ! Test, if we are at a boundary
      if (IelementList(iedge)==0) then
      
        ! No BCs
        Dsolutionvalues(2,ipoint,iedge,1) = Dsolutionvalues(1,ipoint,iedge,1)
        
        ! Dirichlet
        Dsolutionvalues(2,ipoint,iedge,1) = 0.0_dp

      end if
      
    




      ! *** Upwind flux without dimensional splitting (Roe-Flux) and new c ***
      
      ! Get solution values on the in and outside
      DQi(1) = Dsolutionvalues(1,ipoint,iedge,1)
      DQa(1) = Dsolutionvalues(2,ipoint,iedge,1)

      
      DFlux(1) = (normal(1,iedge)*2.0_dp + normal(2,iedge)*1.0_dp)*0.5_dp*(DQi(1)+DQa(1))
      
      if ((normal(1,iedge)+normal(2,iedge))>0.0_dp) then
        DFlux(1) = (normal(1,iedge)*2.0_dp + normal(2,iedge)*1.0_dp)*(DQi(1))
      else
        DFlux(1) = (normal(1,iedge)*2.0_dp + normal(2,iedge)*1.0_dp)*(DQa(1))
      end if
      
      
      DFlux(1) = (normal(1,iedge)*2.0_dp + normal(2,iedge)*1.0_dp)*0.5_dp*(DQi(1)+DQa(1)) + 0.5_dp*abs(normal(1,iedge)*2.0_dp + normal(2,iedge)*1.0_dp)*(DQi(1)-DQa(1))
      
      
        
      
      ! Save the calculated flux
      DfluxValues(1,1,ipoint,iedge) = DFlux(1)
      
      
    end do ! ipoint
  end do ! iedge
                                 
  deallocate(Dsolutionvalues)  
  
  
  end subroutine






  ! ***************************************************************************

!<subroutine>

  subroutine Scalar_flux_sys_block (rdiscretisation,rform, &
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
    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(nvar,itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    integer :: iel, ipoint, ivar, nvar2d, currentvar,nterms
    real(DP) :: dvx, dvy, dx, dy, dsol
    ! (# of terms, ncubp, NEL, nvar2d)
    real(dp), dimension(:,:,:,:), allocatable :: DsolutionValues
    real(dp), dimension(1) :: dQ, DFx, DFy
    real(dp), dimension(:), pointer :: p_ddata
    
    ! nvar2d = 4
    nvar2d = rcollection%p_rvectorQuickAccess1%Nblocks
    
    rdomainIntSubset%ielementDistribution = 1
    
    
    ! rform%itermCount gives the number of additional terms, here 2
    nterms = rform%itermCount
        
    ! Allocate space for the solution values ! (# of derivatives, ncubp, NEL, nvar2d)
    allocate(DsolutionValues(1,size(Dpoints,2),size(Dpoints,3),nvar2d))
    
    ! First evaluate the solution in each point
    do ivar = 1, nvar2d
      call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(ivar), &
                                 rdomainIntSubset, DER_FUNC, DsolutionValues(:,:,:,ivar), 1)
    end do  
 
    ! Loop over all elements and all points   
    do iel = 1, size(Dpoints,3)
      do ipoint = 1, size(Dpoints,2)
        
        ! Get coordinates
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)

        ! Get solution      
        dQ(1) = DsolutionValues(1,ipoint,iel,1)
        
        
        ! Set the coefficients
        Dcoefficients (1,1,ipoint,iel) = 2.0_dp*dQ(1)
        Dcoefficients (1,2,ipoint,iel) = 1.0_dp*dQ(1)
        
      end do
    end do
    
    ! Deallocate the solutionvalues    
    deallocate(DsolutionValues)
    
  end subroutine
  
  
  
  
  
   ! ***************************************************************************

!<subroutine>

  subroutine Scalar_flux_sys_sc (rdiscretisation,rform, &
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
    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset

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

    integer :: iel, ipoint, ivar, nvar2d, currentvar,nterms
    real(DP) :: dvx, dvy, dx, dy, dsol
    ! (# of terms, ncubp, NEL, nvar2d)
    real(dp), dimension(:,:,:,:), allocatable :: DsolutionValues
    real(dp), dimension(3) :: dQ, DFx, DFy
    real(dp), dimension(:), pointer :: p_ddata
    
    !nvar2d = 1
    nvar2d = rcollection%p_rvectorQuickAccess1%Nblocks
    
    rdomainIntSubset%ielementDistribution = 1
    
    ! Which variable (of the system) we are at
    currentvar = rcollection%IquickAccess(1)
    
    ! rform%itermCount gives the number of additional terms, here 2
    nterms = rform%itermCount
    
    ! Allocate space for the solution values ! (# of terms, ncubp, NEL, nvar2d)
    allocate(DsolutionValues(nterms,size(Dpoints,2),size(Dpoints,3),nvar2d))
    
    ! First evaluate the solution in each point
    do ivar = 1, nvar2d
      call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(ivar), &
                                 rdomainIntSubset, DER_FUNC, DsolutionValues(:,:,:,ivar), 1)
    end do  
    
!call    lsyssc_getbase_double(rcollection%p_rvectorQuickAccess1%RvectorBlock(3),p_ddata)
!    write(*,*) p_ddata
                               
                                 

    
    do iel = 1, size(Dcoefficients,3)
      do ipoint = 1, size(Dcoefficients,2)
      
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)
        
        
        ! Get coordinates
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)

        ! Get solution      
        dQ(1) = DsolutionValues(1,ipoint,iel,1)
        
        
        ! Set the coefficients
        Dcoefficients (1,ipoint,iel) = 2.0_dp*dQ(1)
        Dcoefficients (2,ipoint,iel) = 1.0_dp*dQ(1)
        
        
        
      end do
    end do
    
    deallocate(DsolutionValues)
    
  end subroutine



  ! ***************************************************************************

!<subroutine>

  subroutine Euler_flux_sys_sc (rdiscretisation,rform, &
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
    type(t_domainIntSubset), intent(inout)              :: rdomainIntSubset

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

    integer :: iel, ipoint, ivar, nvar2d, currentvar,nterms
    real(DP) :: dvx, dvy, dx, dy, dsol
    ! (# of terms, ncubp, NEL, nvar2d)
    real(dp), dimension(:,:,:,:), allocatable :: DsolutionValues
    real(dp), dimension(4) :: dQ, DFx, DFy
    real(dp), dimension(:), pointer :: p_ddata
    
    !nvar2d = 1
    nvar2d = rcollection%p_rvectorQuickAccess1%Nblocks
    
    rdomainIntSubset%ielementDistribution = 1
    
    ! Which variable (of the system) we are at
    currentvar = rcollection%IquickAccess(1)
    
    ! rform%itermCount gives the number of additional terms, here 2
    nterms = rform%itermCount
    
    ! Allocate space for the solution values ! (# of terms, ncubp, NEL, nvar2d)
    allocate(DsolutionValues(nterms,size(Dpoints,2),size(Dpoints,3),nvar2d))
    
    ! First evaluate the solution in each point
    do ivar = 1, nvar2d
      call fevl_evaluate_sim4 (rcollection%p_rvectorQuickAccess1%RvectorBlock(ivar), &
                                 rdomainIntSubset, DER_FUNC, DsolutionValues(:,:,:,ivar), 1)
    end do  
    
!call    lsyssc_getbase_double(rcollection%p_rvectorQuickAccess1%RvectorBlock(3),p_ddata)
!    write(*,*) p_ddata
                               
                                 

    
    do iel = 1, size(Dcoefficients,3)
      do ipoint = 1, size(Dcoefficients,2)
      
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)
        
        
        ! Euler
        dQ = DsolutionValues(1,ipoint,iel,:)
        !write(*,*) dQ
        
        DFx = Euler_buildFlux(dQ,1)
        DFy = Euler_buildFlux(dQ,2)
        
        Dcoefficients (1,ipoint,iel) = DFx(currentvar)
        Dcoefficients (2,ipoint,iel) = DFy(currentvar)


        
      end do
    end do
    
    deallocate(DsolutionValues)
    
  end subroutine
  
  
  
  
  
  
  
  
  ! ***************************************************************************
  !<subroutine>

  subroutine coeff_implicitDGconvection (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(in)                            :: rform
    
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
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTrial
    
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
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
  integer :: iel, ipoint
  real(dp) :: dx, dy
  
    ! Loop over all elements and cubature points
    do iel = 1, size(Dpoints,3)
      do ipoint = 1, size(Dpoints,2)
      
        ! Get coordinates
        dx = Dpoints(1,ipoint,iel)
        dy = Dpoints(2,ipoint,iel)
        
        dx = 1.0_dp
        dy = 1.0_dp
        
        ! Set coefficients (the velocity vector)
        Dcoefficients(1,ipoint,iel) = -dy
        Dcoefficients(2,ipoint,iel) =  dx
        Dcoefficients(1,ipoint,iel) = dx
        Dcoefficients(2,ipoint,iel) = dy
      
      end do
    end do
    
  end subroutine
  
  
  
  
  !<subroutine>

    subroutine flux_dg_implicitConvection_sim (&
!              Dcoefficients,&
!              DsolVals,&
			  DfluxValues,&
			  rvectorSol,&
			  IelementList,&
			  Dside,&
              normal,&
!              DpointsReal,&
              rintSubSet,&
              rcollection )
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
!  real(DP), dimension(:,:,:), intent(inout) :: DsolVals
  real(DP), dimension(:,:,:), intent(out) :: DfluxValues
  real(DP), dimension(:,:), intent(in) :: normal
!  real(DP), dimension(:,:,:), intent(in) :: DpointsReal
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  type(t_collection), intent(inout), target, optional :: rcollection
  type(t_vectorScalar), intent(in) :: rvectorSol
  integer, dimension(:), intent(in) :: IelementList
  real(DP), dimension(:,:), intent(out) :: Dside
    
  !</input>
  
  !<output>
!  real(DP), dimension(:,:), intent(out) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
  integer :: iedge, ipoint
  real(dp) :: dx, dy
  
  
  do iedge = 1, size(DfluxValues,3)
    do ipoint = 1, size(DfluxValues,2)
    
      dx = rintSubset(1)%p_DcubPtsReal(1,ipoint,iedge)
      dy = rintSubset(1)%p_DcubPtsReal(2,ipoint,iedge)
      dx = 1.0_dp
      dy = 1.0_dp
      
      DfluxValues(1,ipoint,iedge) = -dy*normal(1,iedge) + dx*normal(2,iedge)
      DfluxValues(1,ipoint,iedge) = dx*normal(1,iedge) + dy*normal(2,iedge)
      
      if (DfluxValues(1,ipoint,iedge)>0.0_dp) then
        Dside(1,iedge) = 1.0_dp
        Dside(2,iedge) = 0.0_dp
      else
        Dside(1,iedge) = 0.0_dp
        Dside(2,iedge) = 1.0_dp
      end if
      
    end do
  end do
  
  
  end subroutine
  
  
  
  
  
  
  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_2D_zeros (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
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

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    ! REAL(DP) :: dx,dy
    !
    ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
    !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 0.0_DP
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine getBoundaryValues_2D_ones (Icomponents,rdiscretisation,rboundaryRegion,ielement, &
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

    ! To get the X/Y-coordinates of the boundary point, use:
    !
    ! REAL(DP) :: dx,dy
    !
    ! CALL boundary_getCoords(rdiscretisation%p_rboundary, &
    !     rboundaryRegion%iboundCompIdx, dwhere, dx, dy)

    ! Return zero Dirichlet boundary values for all situations.
    Dvalues(1) = 1.0_DP
  
  end subroutine
  
  


end module
