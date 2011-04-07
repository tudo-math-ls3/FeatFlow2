module shallowwater2d_routines

  use fsystem
  use storage
  use triangulation
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use boundary
  use shallowwater2d_callback
  use collection
  use linearformevaluation

  implicit none
  

  type t_array
     ! Pointer to the double-valued matrix or vector data
     real(DP), dimension(:), pointer :: Da
  end type t_array

integer, parameter :: scalardisstype = 1
   !integer, parameter				:: nvar2d = 3
  
  
  integer, parameter :: deltabtype = 3
  
  ! Unity matrix
  real(DP), dimension(3,3) :: Eye = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
  
  ! If the water heigth should go under this level, the node will be considered as empty
  real(dp), parameter :: clipwater = sys_epsreal


contains


  ! Limiter Function for TVD
  real(DP) function limiterfunc(z,n,limiter)
    implicit none
    real(DP), intent(IN)		:: z, n				! z = nominator, n = denominator of slope ratio
    integer, intent(IN)			:: limiter  		! choice of limiter (1 = Minmod, 4 = Superbee)
    real(DP)					:: r				! slope ratio

    limiterfunc = 0.0_DP

    if (abs(n)<1e-8) then
       limiterfunc = 0.0_DP
    else
       r=z/n
       select case(limiter)
       case(1) ! MinMod
          limiterfunc = max(0.0_DP,min(1.0_DP,r))
       case(2) ! Van Leer
          limiterfunc = (r+abs(r))/(1.0_DP+abs(r))
       case(3) ! MC
          limiterfunc = max(0.0_DP, min(2.0_DP,(1.0_DP+r)/2.0_DP,2.0_DP*r))
       case(4) ! Superbee
          limiterfunc = max(0.0_DP,min(1.0_DP,2.0_DP*r),min(2.0_DP,r))
       end select
    end if
  end function limiterfunc


  ! 2 parameter limiter function for TVD
  real(DP) function limiterfunc2(z,n,limiter)
    implicit none
    real(DP), intent(IN)		:: z, n				! z = nominator, n = denominator of slope ratio
    integer, intent(IN)			:: limiter  		! choice of limiter (1 = Minmod, 4 = Superbee)
    real(DP)					:: h1				! temporary variable

    h1 = (sign(1.0_DP,z)+sign(1.0_DP,n))/2.0d0

    select case(limiter)
    case(1) ! MinMod
       limiterfunc2 = h1*min(abs(z),abs(n))
    case(2) ! Van Leer
       if (abs(z)+abs(n)>0.001_DP) then	
          limiterfunc2 = h1*2.0_DP*abs(z*n)/(abs(n)+abs(z))
       else
          limiterfunc2 = 0.0_DP
       end if
    case(3) ! MC
       limiterfunc2=h1*min(2.0_DP*abs(z),(abs(z)+abs(n))/2.0_DP,2.0_DP*abs(n))
    case(4) ! Superbee
       limiterfunc2=h1*max(min(2.0_DP*abs(z),abs(n)),min(abs(z),2.0_DP*abs(n)))
    end select
  end function limiterfunc2



  ! This function returns the Roe mean values
  function calculateQroe(Ql, Qr) result(Qroe)

    ! The left and right Q values
    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Ql, Qr

    ! The computed Roe values
    real(DP), dimension(3)					:: Qroe

    ! temp variables
    real(DP)		:: whl, whr, denom, hl, hr

    ! Choose kind of mean value
    integer, parameter :: mvk = 1


    select case (mvk)

      case (0) ! Roe-meanvalues

        ! Set the height variables
        hl=Ql(1)
        hr=Qr(1)

        denom = sqrt(hl)+sqrt(hr)
        whl = 1.0_DP/sqrt(hl)
        whr = 1.0_DP/sqrt(hr)

        Qroe(1) = sqrt(hl*hr)
        Qroe(2) = Qroe(1)*(whl*Ql(2)+whr*Qr(2))/denom
        Qroe(3) = Qroe(1)*(whl*Ql(3)+whr*Qr(3))/denom

      case (1) ! Artihmetic mean

        Qroe = 0.5_dp*(Ql+Qr)

      end select

  end function calculateQroe



  ! This routine builds the jacobi matrix in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildJacobi(Q,d,g) result(J)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)	:: J

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: q

    ! The gravitational konstant
    real(DP), intent(IN)		        	:: g

    integer, intent(IN)                     :: d

    ! primitive variables
    real(DP)                                :: h, u, v, c

    ! Calculate primitive variables
    h=Q(1)
    if (h<clipwater) then
      h=0.0_dp
      u=0.0_dp
      v=0.0_dp
    else
      u=Q(2)/Q(1)
      v=Q(3)/Q(1)
    end if

    if (d==1) then
       ! build Jacobi matrix in x direction
       J(1,1) = 0.0_DP
       J(2,1) = -u**2.0_DP + g*h	! -u^2+gh
       J(3,1) = -u*v			! -uv
       J(1,2) = 1.0_DP
       J(2,2) = 2.0_DP*u				! 2u
       J(3,2) = v						! v
       J(1,3) = 0.0_DP
       J(2,3) = 0.0_DP
       J(3,3) = u						! u
    else
       ! build Jacobi matrix in y direction
       J(1,1) = 0.0_DP
       J(2,1) = -u*v			! -uv
       J(3,1) = -v**2.0_DP + g*h	! -v^2+gh
       J(1,2) = 0.0_DP
       J(2,2) = v						! v
       J(3,2) = 0.0_DP
       J(1,3) = 1.0_DP
       J(2,3) = u						! u
       J(3,3) = 2.0_DP*v				! 2v
    end if

  end function buildJacobi



  ! This routine builds the jacobi matrix in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildDissipation(Qi,Qj,k,d,g,cxij,cyij,cxji,cyji) result(Dij)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)    :: Dij

    ! The solution components q1 = h, q2 = uh, q3 = vh of the Roe-values
    real(DP), dimension(3), intent(IN)      :: Qi,Qj
    
    ! The kind of dissipation (scalar, tensorial, ...)
    integer, intent(IN)                     :: k

    ! The gravitational konstant
    real(DP), intent(IN)                    :: g
    
    ! The coefficients of the jacobians
    real(DP), intent(IN)                    :: cxij, cyij,cxji,cyji

    integer, intent(IN)                     :: d
    
    ! The Roe-Mean-Values
    real(DP), dimension(nvar2d)             :: Qroe, Q, Qroeij

    ! local variables
    real(DP)                                :: scalefactor1, scalefactor2, cRoe, uRoe, vRoe, lambda, scalardissipation, scalefactor
    
    real(dp), dimension(3) :: w,diff
    real(dp), dimension(3,3) :: A, Aij, Bij
    integer :: l
    real(dp) :: ci, cj, vi, vj, ui, uj
    real(dp) :: lambda1, lambda3
    real(dp) :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
    real(dp), dimension(nvar2d,nvar2d) :: JacobixRoeij, JacobiyRoeij
    real(dp)::h,u,v


    select case (k)
      case (0) ! No Dissipation

        Dij = 0

      case (1)! Here we use scalar dissipation

! !        scalefactor = sqrt(cxij**2.0_DP+cyij**2.0_DP)
! !        cRoe = sqrt(g*Q(1))    ! c_Roe=sqrt(g*h)
! !        ! cRoe = sqrt(0.5*gravconst*(Qi(1)+Qj(1)))
! !        uRoe = Q(2)/Q(1) 
! !        vRoe = Q(3)/Q(1)
! !        lambda = sqrt((cxij/scalefactor*uRoe)**2.0_DP+(cyij/scalefactor*vRoe)**2.0_DP)+cRoe
! !        scalarDissipation = scalefactor*lambda
! !
! !        Dij=scalarDissipation*Eye
!         
!         Qroe = calculateQroe(Qi, Qj)
!         cRoe = sqrt(0.5_DP*g*(Qi(1)+Qj(1)))
!         uRoe = QRoe(2)/QRoe(1) 
!         vRoe = QRoe(3)/QRoe(1)
!         
! !         if ((Qi(1)<1e-6_dp).or.(Qj(1)<1e-6_dp)) then
! !           uRoe = max(Qi(2)/Qi(1),Qj(2)/Qj(1))
! !           vRoe = max(Qi(3)/Qi(1),Qj(3)/Qj(1))
! !           cRoe = 2.0_dp*max(sqrt(g*Qi(1)),sqrt(g*Qj(1)))
! !         end if
!         
!         scalefactor = sqrt(cxij**2.0_DP+cyij**2.0_DP)
!         lambda = abs(cxij*uRoe+cyij*vRoe)/scalefactor+2.0_DP*cRoe
!         scalarDissipation = scalefactor*lambda
!         
!         scalarDissipation = abs(cxij*uRoe) + cxij*cRoe + abs(cyij*vRoe) + cyij*cRoe
!         
!         Dij=scalarDissipation*Eye
!         
!         
!         
!         
!         
!         A=cxij*buildJacobi(Qroe,1,g)+cyij*buildJacobi(Qroe,2,g)
!         v=(/1.0_dp,1.0_dp,1.0_dp/)
!         v=1.0_dp/sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))*v
!         
!         
!         findeig: do l = 1, 100
!           !w=A*v
!           w=matmul(A,v)
!           
!           
!            if (sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3)).le.1e-12_dp) exit findeig
!           
!           
!           w=1.0_dp/sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))*w
!           !diff=v-w
!           !if (sqrt(diff(1)*diff(1)+diff(2)*diff(2)+diff(3)*diff(3)).le.1e-3_dp) exit findeig
!           v=w
!         end do findeig
!         !w=A*w
!         w=matmul(A,w)
!         
!         scalarDissipation = sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
!         
! 
!         
!         Dij=scalarDissipation*Eye




        ! Test, if node i is dry and calculate the velocities
        if (Qi(1)<clipwater) then
          ci = 0.0_dp
          ui = 0.0_dp
          vi = 0.0_dp
        else
         ci = sqrt(g*Qi(1))
         ui = Qi(2)/Qi(1)
         vi = Qi(3)/Qi(1)
!          if(ui/ci>0.8_dp) write(*,*)ui/ci
!          ! Entropy Fix
!          lambda1 = ui-ci
!          if (abs(ui-ci)<0.1_dp) then
!           
!          end if
         
        end if

        ! Test, if node j is dry and calculate the velocities
        if (Qj(1)<clipwater) then
          cj = 0.0_dp
          uj = 0.0_dp
          vj = 0.0_dp
        else
          cj = sqrt(g*Qj(1))
          uj = Qj(2)/Qj(1)
          vj = Qj(3)/Qj(1)
!           if(uj/cj>0.8_dp)write(*,*)uj/cj
        end if

        scalefactor1 = sqrt(cxij**2.0_DP+cyij**2.0_DP)
        scalefactor2 = sqrt(cxji**2.0_DP+cyji**2.0_DP)
        scalarDissipation = max(abs(cxij*uj+cyij*vj)+scalefactor1*cj,&
                                abs(cxji*ui+cyji*vi)+scalefactor2*ci)
!         scalarDissipation = max(abs(cxij*uj+cyij*vj)+scalefactor1*cj,&
!                                 abs(cxji*ui+cyji*vi)+scalefactor2*ci,&
!                                 abs(cxij*ui+cyij*vi)+scalefactor1*ci,&
!                                 abs(cxji*uj+cyji*vj)+scalefactor2*cj)

    

        Dij=scalarDissipation*Eye
!         
!         
!         
!         Qroe = calculateQroe(Qi,Qj)
!         cRoe = sqrt(g/2.0_dp * ( Qi(1) + Qj(1)))
!         lambda1=(cxij/scalefactor1*Qroe(2)/Qroe(1)+cyij/scalefactor1*Qroe(3)/Qroe(1))-cRoe
!         lambda3=(cxij/scalefactor1*Qroe(2)/Qroe(1)+cyij/scalefactor1*Qroe(3)/Qroe(1))+cRoe
! !write(*,*) lambda1, lambda3
!         if(abs(lambda1)<0.1_dp) then
!         write(*,*)'*'
!          ! lambda1l=(cxij/scalefactor1*Qroe(2)/Qroe(1)+cyij/scalefactor1*Qroe(3)/Qroe(1))-cRoe
!         
!         end if
!         
        
        
        
        
        
        
        
        
        
!         
!         
!         
!               A=cxij*buildJacobi(Qroe,1,g)+cyij*buildJacobi(Qroe,2,g)
!         v=(/1.0_dp,1.0_dp,1.0_dp/)
!         v=1.0_dp/sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))*v
!         
!         
!         findeig1: do l = 1, 100
!           !w=A*v
!           w=matmul(A,v)
!           
!           
!            if (sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3)).le.1e-12_dp) exit findeig1
!           
!           
!           w=1.0_dp/sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))*w
!           !diff=v-w
!           !if (sqrt(diff(1)*diff(1)+diff(2)*diff(2)+diff(3)*diff(3)).le.1e-3_dp) exit findeig
!           v=w
!         end do findeig1
!         !w=A*w
!         w=matmul(A,w)
!         
!         scalarDissipation = sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
!         
!         A=cxji*buildJacobi(Qroe,1,g)+cyji*buildJacobi(Qroe,2,g)
!         v=(/1.0_dp,1.0_dp,1.0_dp/)
!         v=1.0_dp/sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))*v
!         
!         
!         findeig2: do l = 1, 100
!           !w=A*v
!           w=matmul(A,v)
!           
!           
!            if (sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3)).le.1e-12_dp) exit findeig2
!           
!           
!           w=1.0_dp/sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))*w
!           !diff=v-w
!           !if (sqrt(diff(1)*diff(1)+diff(2)*diff(2)+diff(3)*diff(3)).le.1e-3_dp) exit findeig
!           v=w
!         end do findeig2
!         !w=A*w
!         w=matmul(A,w)
!         
!         scalarDissipation = max(scalardissipation,sqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3)))
!         
!         
!         
!         
!         Dij=scalarDissipation*Eye
        
        
        
        
       ! This one works for non dry bed
!        Qroeij = calculateQroe(Qi,Qj)
! 
!        JcoeffxA = (cxij-cxji)/2.0_DP
!        JcoeffyA = (cyij-cyji)/2.0_DP
!        JcoeffxB = (cxij+cxji)/2.0_DP
!        JcoeffyB = (cyij+cyji)/2.0_DP
! 
!        scalarDissipation = max(scalardissipation,&
!             abs(JcoeffxA*maxval(abs(buildeigenvalues(Qroeij,1,g))))+ &
!             abs(JcoeffyA*maxval(abs(buildeigenvalues(Qroeij,2,g)))))
!        
!        scalardissipation = max(scalardissipation,SYS_EPSREAL_DP)
! 
!        Dij = scalarDissipation*Eye
       
       
       
       
       
       
       
!        ! Compute the jacobi matrices for x- and y- direction
!        JacobixRoeij = buildJacobi(Qroeij,1,g)
!        JacobiyRoeij = buildJacobi(Qroeij,2,g)
! 
!        ! Compute the coeffitients for the jacobi matrices
!        JcoeffxA = (CXij-CXji)/2.0_DP
!        JcoeffyA = (CYij-CYji)/2.0_DP
!        JcoeffxB = (CXij+CXji)/2.0_DP
!        JcoeffyB = (CYij+CYji)/2.0_DP
! 
!        ! Now we can compute Aij and Bij
!        Aij = JcoeffxA*JacobixRoeij + JcoeffyA*JacobiyRoeij
!        Bij = JcoeffxB*JacobixRoeij + JcoeffyB*JacobiyRoeij
!        
!        Dij = Dij + Bij

 Qroe = calculateQroe(Qi, Qj)
 h=Qroe(1)
 if(h<clipwater)then
  h=0.0_dp
  u=0.0_dp
  v=0.0_dp
 else
  u=Qroe(2)/Qroe(1)
  v=Qroe(3)/Qroe(1)
 end if
 cRoe = sqrt(g/2.0_dp * ( Qi(1) + Qj(1)))
 scalefactor1 = sqrt(cxij**2.0_DP+cyij**2.0_DP)
 scalefactor2 = sqrt(cxji**2.0_DP+cyji**2.0_DP)
 scalarDissipation = max(scalardissipation, &
  max( (abs(-cxij*u-cyij*v)+scalefactor1*cRoe), (abs(-cxji*u-cyji*v)+scalefactor2*cRoe) )       )
 Dij = scalarDissipation*Eye
        

      case (2)! Here we use tensorial dissipation in direction d

        Qroe = calculateQroe(Qi, Qj)
        
        Dij = matmul(buildTrafo(Qroe,d,g),&
               matmul(buildaLambda(Qroe,d,g),&
               buildinvTrafo(Qroe,d,g)))

    end select



!        ! Dissipation of Rusanov type
!        scalardissipation = &
!              max(abs(p_CXdata(ij)*Qj(2)/Qj(1)+p_Cydata(ij)*Qj(3)/Qj(1)) + sqrt(gravconst*Qj(1)),&
!                  abs(p_CXdata(ji)*Qi(2)/Qi(1)+p_Cydata(ji)*Qi(3)/Qi(1)) + sqrt(gravconst*Qi(1)) )




  end function buildDissipation



  ! This routine builds the trafo matrix Rij in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildTrafo(Q,d,g) result(Rij)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)	:: Rij

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    ! The gravitational konstant
    real(DP), intent(IN)		        	:: g

    ! the direction: d=1: x-direction, d=2: y-direction
    integer, intent(IN)                     :: d

    ! speed of gravitational waves
    real(DP)                                :: c

    ! temporary variable
    real(DP)                                :: coeff, h, u, v

    ! Calculate primitive variables
    h=Q(1)
    u=Q(2)/Q(1)
    v=Q(3)/Q(1)

    ! compute c = sqrt(g*h)
    c = sqrt(g*h)



    if (d==1) then
       ! build trafo matrix in x direction
       Rij(1,1) = 1.0_DP
       Rij(2,1) = u-c
       Rij(3,1) = v
       Rij(1,2) = 0.0_DP
       Rij(2,2) = 0.0_DP
       Rij(3,2) = 1.0_DP
       Rij(1,3) = 1.0_DP
       Rij(2,3) = u+c
       Rij(3,3) = v
    else
       ! build trafo matrix in y direction
       Rij(1,1) = 1.0_DP
       Rij(2,1) = u
       Rij(3,1) = v-c
       Rij(1,2) = 0.0_DP
       Rij(2,2) = 1.0_DP
       Rij(3,2) = 0.0_DP
       Rij(1,3) = 1.0_DP
       Rij(2,3) = u
       Rij(3,3) = v+c
    end if

  end function buildTrafo



  ! This routine builds the inv trafo matrix Rij^{-1} in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildInvTrafo(Q,d,g) result(invRij)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)	:: invRij

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    ! The gravitational konstant
    real(DP), intent(IN)		        	:: g

    ! the direction: d=1: x-direction, d=2: y-direction
    integer, intent(IN)                     :: d

    ! speed of gravitational waves
    real(DP)                                :: c

    ! temporary variable
    real(DP)                                :: coeff, h, u, v

    ! Calculate primitive variables
    h=Q(1)
    u=Q(2)/Q(1)
    v=Q(3)/Q(1)

    ! compute c = sqrt(g*h)
    c = sqrt(g*h)

    ! temporary variable
    coeff = 1.0_DP/(2.0_DP*c)


    if (d==1) then
       ! build inv trafo matrix in x direction
       invRij(1,1) = coeff*(u+c)
       invRij(2,1) = -v
       invRij(3,1) = coeff*(c-u)
       invRij(1,2) = -coeff
       invRij(2,2) = 0.0_DP
       invRij(3,2) = coeff
       invRij(1,3) = 0.0_DP
       invRij(2,3) = 1.0_DP
       invRij(3,3) = 0.0_DP
    else
       ! build inv trafo matrix in y direction
       invRij(1,1) = coeff*(v+c)
       invRij(2,1) = -u
       invRij(3,1) = coeff*(c-v)
       invRij(1,2) = 0
       invRij(2,2) = 1.0_DP
       invRij(3,2) = 0
       invRij(1,3) = -coeff
       invRij(2,3) = 0.0_DP
       invRij(3,3) = coeff
    end if

  end function buildInvTrafo



  ! This routine builds the diagonalised jacobi matrix Lambda in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildLambda(Q,d,g) result(Lambda)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)	:: Lambda

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    ! The gravitational konstant
    real(DP), intent(IN)		        	:: g

    ! the direction: d=1: x-direction, d=2: y-direction
    integer, intent(IN)                     :: d

    ! speed of gravitational waves
    real(DP)                                :: c

    ! temporary variable
    real(DP)                                :: coeff, h, u, v

    ! Calculate primitive variables
    h=Q(1)
    u=Q(2)/Q(1)
    v=Q(3)/Q(1)

    ! compute c = sqrt(g*h)
    c = sqrt(g*h)


    if (d==1) then
       ! build Lambda in x direction
       Lambda(1,1) = u-c
       Lambda(2,1) = 0.0_DP
       Lambda(3,1) = 0.0_DP
       Lambda(1,2) = 0.0_DP
       Lambda(2,2) = u
       Lambda(3,2) = 0.0_DP
       Lambda(1,3) = 0.0_DP
       Lambda(2,3) = 0.0_DP
       Lambda(3,3) = u+c
    else
       ! build Lambda in y direction
       Lambda(1,1) = v-c
       Lambda(2,1) = 0.0_DP
       Lambda(3,1) = 0.0_DP
       Lambda(1,2) = 0.0_DP
       Lambda(2,2) = v
       Lambda(3,2) = 0.0_DP
       Lambda(1,3) = 0.0_DP
       Lambda(2,3) = 0.0_DP
       Lambda(3,3) = v+c
    end if

  end function buildLambda



  ! This routine builds the absolute value of the diagonalised jacobi
  ! matrix aLambda in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildaLambda(Q,d,g) result(aLambda)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)	:: aLambda

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    ! The gravitational konstant
    real(DP), intent(IN)		        	:: g

    ! the direction: d=1: x-direction, d=2: y-direction
    integer, intent(IN)                     :: d

    ! speed of gravitational waves
    real(DP)                                :: c

    ! temporary variable
    real(DP)                                :: coeff, h, u, v

    ! Calculate primitive variables
    h=Q(1)
    u=Q(2)/Q(1)
    v=Q(3)/Q(1)

    ! compute c = sqrt(g*h)
    c = sqrt(g*h)


    if (d==1) then
       ! build aLambda in x direction
       aLambda(1,1) = abs(u-c)
       aLambda(2,1) = 0.0_DP
       aLambda(3,1) = 0.0_DP
       aLambda(1,2) = 0.0_DP
       aLambda(2,2) = abs(u)
       aLambda(3,2) = 0.0_DP
       aLambda(1,3) = 0.0_DP
       aLambda(2,3) = 0.0_DP
       aLambda(3,3) = abs(u+c)
    else
       ! build aLambda in y direction
       aLambda(1,1) = abs(v-c)
       aLambda(2,1) = 0.0_DP
       aLambda(3,1) = 0.0_DP
       aLambda(1,2) = 0.0_DP
       aLambda(2,2) = abs(v)
       aLambda(3,2) = 0.0_DP
       aLambda(1,3) = 0.0_DP
       aLambda(2,3) = 0.0_DP
       aLambda(3,3) = abs(v+c)
    end if

  end function buildaLambda






  ! This routine returns the eigenvalues of the jacobi matrix in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildEigenvalues(Q,d,g) result(Eigenvalues)

    ! The jacobi matrix in direction d
    real(DP), dimension(3)	:: Eigenvalues

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    ! The gravitational konstant
    real(DP), intent(IN)		        	:: g

    ! the direction: d=1: x-direction, d=2: y-direction
    integer, intent(IN)                     :: d

    ! speed of gravitational waves
    real(DP)                                :: c

    ! temporary variable
    real(DP)                                :: coeff, h, u, v

    ! Calculate primitive variables
    h=Q(1)
    if (h<clipwater) then
      h = 0.0_dp
      u = 0.0_dp
      v = 0.0_dp
    else
      u = Q(2)/Q(1)
      v = Q(3)/Q(1)
    end if

    ! compute c = sqrt(g*h)
    c = sqrt(g*h)

    if (d==1) then
       ! build eigenvalues in x direction
       Eigenvalues(1) = u-c
       Eigenvalues(2) = u
       Eigenvalues(3) = u+c
    else
       ! build eigenvalues in y direction
       Eigenvalues(1) = v-c
       Eigenvalues(2) = v
       Eigenvalues(3) = v+c
    end if

  end function buildEigenvalues



  ! This routine returns the eigenvalues of the jacobi matrix in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildFlux(Q,d,g) result(Flux)

    ! The flux vector in direction d at Q
    real(DP), dimension(3)	:: Flux

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    ! The gravitational konstant
    real(DP), intent(IN)		        	:: g

    ! the direction: d=1: x-direction, d=2: y-direction
    integer, intent(IN)                     :: d

    ! speed of gravitational waves
    real(DP)                                :: c

    ! temporary variable
    real(DP)                                :: coeff, h, u, v

    ! Calculate primitive variables
    !h=Q(1)
    !u=Q(2)/Q(1)
    !v=Q(3)/Q(1)

    ! compute c = sqrt(g*h)
    !c = sqrt(g*h)

    ! Test for dry bed case
     if (Q(1)<clipwater) then
        !dry bed case
        Flux=0
     else
       !wet bed case
       if (d==1) then
          ! build Flux in x direction
          Flux(1) = Q(2)
          Flux(2) = Q(2)*Q(2)/Q(1)+0.5_DP*g*Q(1)*Q(1)
          Flux(3) = Q(2)*Q(3)/Q(1)
       else
          ! build Flux in y direction
          Flux(1) = Q(3)
          Flux(2) = Q(2)*Q(3)/Q(1)
          Flux(3) = Q(3)*Q(3)/Q(1)+0.5_DP*g*Q(1)*Q(1)
       end if
     end if ! dry or wet bed
  end function buildFlux




  ! This routine walks over the solution vector and clips it if it's smaller than 1e-6
  subroutine ClipHeight (rarraySol, neq)

    ! parameter values
    type(t_array), dimension(nvar2d), intent(INOUT)    :: rarraySol
    integer, intent(in) :: neq

    ! variables
    integer                             :: i



    do i = 1, NEQ
       if(rarraySol(1)%Da(i)<clipwater) then
        rarraySol(1)%Da(i) = 0.0_dp
        rarraySol(2)%Da(i) = 0.0_dp
        rarraySol(3)%Da(i) = 0.0_dp
       end if


    end do


  end subroutine ClipHeight





  ! This routine builds the preconditioner for the shallow water system
  subroutine BuildShallowWaterPreconditioner (rmatrixBlockP, &
       rarrayP, rarraySol, p_CXdata, p_CYdata, &
       p_MLdata, p_Kdiagonal, p_kedge, &
       NEQ, nedge, theta, dt, gravconst)

    ! parameter values
    type(t_matrixBlock), intent(INOUT)              :: rmatrixBlockP
    type(t_array), dimension(nvar2d), intent(INOUT) :: rarrayP
    type(t_array), dimension(nvar2d), intent(IN)    :: rarraySol
    real(DP), dimension(:), pointer, intent(IN)     :: p_CXdata, p_CYdata, p_MLdata
    integer, dimension(:), pointer, intent(IN)      :: p_Kdiagonal
    integer, dimension(:,:), pointer                :: p_Kedge
    integer, intent(IN)                             :: NEQ,	nedge
    real(DP), intent(IN)                            :: theta, gravconst, dt

    ! variables
    integer                             :: i, j, l, ii, jj, ij, ji, iedge, ivar
    real(DP), dimension(nvar2d)         :: deltaQij, Qi, Qj, Qroeij
    real(DP), dimension(nvar2d,nvar2d)  :: JacobixRoeij, JacobiyRoeij
    real(DP)                            :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
    real(DP), dimension(nvar2d,nvar2d)  :: Kij, Kji, Dij
    real(DP)                            :: scalarDissipation, scalefactor, lambda
    real(DP)                            :: cRoe, uRoe, vRoe
    real(dp), dimension(3,3) :: Aij, Bij


    ! Assemble the preconditioner rmatrixBlockP: P = ML - theta*dt*L
    ! As we use the block jacobi method, we only need the main diagonal blocks


    ! First set all matrix blocks on the main diagonal of P equal to ML
    call lsysbl_clearMatrix(rmatrixBlockP)
    do i = 1, NEQ
       ii = p_Kdiagonal(i)
       do ivar = 1, nvar2d
          rarrayP(ivar)%Da(ii) = p_MLdata(i)
       end do
    end do

    ! Now walk over all edges ij and compute Aij, Bij and Dij
    ! Then write the entries of these Matrices to the corresponding entries
    ! of the matrices on the main diagonal blocks of P
    do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)

       ! get solution values at node i and j
       Qi = (/rarraySol(1)%Da(i),rarraySol(2)%Da(i),rarraySol(3)%Da(i)/)
       Qj = (/rarraySol(1)%Da(j),rarraySol(2)%Da(j),rarraySol(3)%Da(j)/)

       ! compute deltaQij = Qi - Qj
       deltaQij = Qi - Qj

       ! Compute the Roe-meanvalue for this edge
       Qroeij = calculateQroe(Qi, Qj)
       !Qroeij = 0.5_dp*(Qi+Qj)

       ! Compute the jacobi matrices for x- and y- direction
       JacobixRoeij = buildJacobi(Qroeij,1,gravconst)
       JacobiyRoeij = buildJacobi(Qroeij,2,gravconst)

       ! Compute the coeffitients for the jacobi matrices
!        JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
!        JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
!        JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
!        JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP

       ! Now we can compute Aij and Bij
!        Aij = JcoeffxA*JacobixRoeij + JcoeffyA*JacobiyRoeij
!        Bij = JcoeffxB*JacobixRoeij + JcoeffyB*JacobiyRoeij

       ! Compute the high order Matrizes
       Kij = p_CXdata(ij)*JacobixRoeij + p_CYdata(ij)*JacobiyRoeij
       Kji = p_CXdata(ji)*JacobixRoeij + p_CYdata(ji)*JacobiyRoeij

       ! Here we use scalar dissipation
       Dij = buildDissipation(Qi,Qj,1,0,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))

       ! Now add the entries of Aij and Dij to their corresponding entries in P
       ! P = M^L - theta*dt*L
       do l = 1, nvar2d
!           rarrayP(l)%Da(ii) = rarrayP(l)%Da(ii) - dt*theta*(+Aij(l,l)+Bij(l,l)-Dij(l,l))
!           rarrayP(l)%Da(ij) = rarrayP(l)%Da(ij) - dt*theta*(-Aij(l,l)-Bij(l,l)+Dij(l,l))
!           rarrayP(l)%Da(ji) = rarrayP(l)%Da(ji) - dt*theta*(+Aij(l,l)-Bij(l,l)+Dij(l,l))
!           rarrayP(l)%Da(jj) = rarrayP(l)%Da(jj) - dt*theta*(-Aij(l,l)+Bij(l,l)-Dij(l,l))
          
!           rarrayP(l)%Da(ii) = rarrayP(l)%Da(ii) - dt*theta*(+Aij(l,l)-Dij(l,l))
!           rarrayP(l)%Da(ij) = rarrayP(l)%Da(ij) - dt*theta*(-Aij(l,l)+Dij(l,l))
!           rarrayP(l)%Da(ji) = rarrayP(l)%Da(ji) - dt*theta*(+Aij(l,l)+Dij(l,l))
!           rarrayP(l)%Da(jj) = rarrayP(l)%Da(jj) - dt*theta*(-Aij(l,l)-Dij(l,l))

          rarrayP(l)%Da(ii) = rarrayP(l)%Da(ii) + dt*theta*(-Kij(l,l)+Dij(l,l))
          rarrayP(l)%Da(ij) = rarrayP(l)%Da(ij) + dt*theta*(+Kij(l,l)-Dij(l,l))
          rarrayP(l)%Da(ji) = rarrayP(l)%Da(ji) + dt*theta*(+Kji(l,l)-Dij(l,l))
          rarrayP(l)%Da(jj) = rarrayP(l)%Da(jj) + dt*theta*(-Kji(l,l)+Dij(l,l))
       end do

    end do

  end subroutine BuildShallowWaterPreconditioner






  ! This routine builds the RHS for the shallow water system
  subroutine BuildShallowWaterRHS (&
       rarrayRhs, rarraySol, rrhsBlock, rsolBlock, &
       rmatrixML, p_CXdata, p_CYdata, p_MLdata, &
       p_BXdata, p_BYdata, p_BSXdata, p_BSYdata, &
       h_fld1, p_fld1, p_fld2, &
       p_Kdiagonal, p_Kedge, NEQ, nedge, &
       theta, dt, gravconst, Method, limiter)

    ! parameter values
    type(t_array), dimension(nvar2d), intent(INOUT) :: rarrayRhs
    type(t_array), dimension(nvar2d), intent(IN)    :: rarraySol
    type(t_vectorBlock), intent(INOUT)              :: rrhsBlock
    type(t_vectorBlock), intent(IN)                 :: rsolBlock
    type(t_matrixScalar), intent(IN)                :: rmatrixML
    integer                                         :: h_fld1
    real(DP), dimension(:,:), pointer               :: p_fld1, p_fld2
    real(DP), dimension(:), pointer, intent(IN)     :: p_CXdata, p_CYdata, p_BXdata, p_BYdata, p_MLdata
    real(DP), dimension(:), pointer, intent(IN)     :: p_BSXdata, p_BSYdata
    integer, dimension(:), pointer, intent(IN)      :: p_Kdiagonal
    integer, dimension(:,:), pointer                :: p_Kedge
    integer, intent(IN)                             :: NEQ,	nedge
    real(DP), intent(IN)                            :: theta, gravconst, dt
    integer, intent(IN)                             :: Method, limiter


    ! variables
    integer                             :: i, j, l, ii, jj, ij, ji, d, iedge, ivar, inode
    real(DP), dimension(nvar2d)         :: deltaQij, Qi, Qj, Qroeij
    real(DP), dimension(nvar2d,nvar2d)  :: JacobixRoeij, JacobiyRoeij
    real(DP)                            :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
    real(DP), dimension(nvar2d,nvar2d)  :: Aij, Bij, Dij
    real(DP)                            :: scalarDissipation, scalefactor, lambda
    real(DP), dimension(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj, deltaBi, deltaBj
    real(DP)                            :: cRoe, uRoe, vRoe
    ! for TVD
    real(DP), dimension(nvar2d,nvar2d)  :: invRij, Rij
    real(DP), dimension(nvar2d)         :: deltaWij, eigenvalues, deltaFij
    real(DP)                            :: scalefactor1, scalefactor2, deltaak, deltaWijHat
    real(DP)                            :: deltaaplus, deltaaminus, deltabplus, deltabminus
    integer                             :: upwindnode


    ! Compute the RHS b: b=ML*Q^n
    do ivar = 1, nvar2d
       call lsyssc_scalarMatVec (rmatrixML, &
            rsolBlock%Rvectorblock(ivar), &
            rrhsBlock%Rvectorblock(ivar), &
            1._DP, 0._DP)
    end do


    ! Now walk over all edges ij and compute Aij, Bij and Dij
    ! and update the RHS b with the high/low order flux:
    ! b = b + dt*(1-theta)*K*u + dt*(1-theta)*D*u
    do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)

       ! get solution values at node i and j
       Qi = (/rarraySol(1)%Da(i),rarraySol(2)%Da(i),rarraySol(3)%Da(i)/)
       Qj = (/rarraySol(1)%Da(j),rarraySol(2)%Da(j),rarraySol(3)%Da(j)/)

       ! compute deltaQij = Qi - Qj
       deltaQij = Qi - Qj

       ! Compute the Roe-meanvalue for this edge
       Qroeij = calculateQroe(Qi, Qj)
       !Qroeij=0.5_dp*(Qi+Qj)

       ! Calculate the galerkin fluxes
       ! deltaKi = c_{ij}*(F(Q_i)-F(Q_j))
       ! deltaKj = c_{ji}*(F(Q_j)-F(Q_i))
       deltaKi = p_CXdata(ij)*(buildFlux(Qi,1,gravconst)-buildFlux(Qj,1,gravconst))+ &
                 p_CYdata(ij)*(buildFlux(Qi,2,gravconst)-buildFlux(Qj,2,gravconst))
       deltaKj = p_CXdata(ji)*(buildFlux(Qj,1,gravconst)-buildFlux(Qi,1,gravconst))+ &
                 p_CYdata(ji)*(buildFlux(Qj,2,gravconst)-buildFlux(Qi,2,gravconst))

       ! Now choose the artificial diffusion method
       select case (Method)
        case (0,3)
          ! Here we do not add artificial diffusion - so we have the high order method
          Dij = 0
        case (1,4,6)
          ! Here we use scalar dissipation
          Dij = buildDissipation(Qi,Qj,1,0,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
        case (2,5,7)
          ! Here we use tensorial dissipation
          Dij = abs(p_CXdata(ij))* buildDissipation(Qi,Qj,2,1,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))+&
                abs(p_CYdata(ij))* buildDissipation(Qi,Qj,2,2,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
       end select

!        ! deltaD
!        deltaDi = matmul(Dij,-deltaQij)
!        deltaDj = -deltaDi

!        ! Now we take care of the bottom profile
!        select case (deltabtype)
!         case(0)
!           DeltaBi = 0.0_dp
!           DeltaBj = 0.0_dp
!         case(1)
!           DeltaBi = (/ 0.0_DP, p_BXdata(ij) * Qi(1), p_BYdata(ij) * Qi(1) /)
!           DeltaBj = (/ 0.0_DP, p_BXdata(ji) * Qj(1), p_BYdata(ji) * Qj(1) /)
!         case(2)
!           DeltaBi = (/ 0.0_DP, p_BXdata(ij) * 0.5_dp*(Qi(1)+Qj(1)), p_BYdata(ij) * 0.5_dp*(Qi(1)+Qj(1)) /)
!           DeltaBj = (/ 0.0_DP, p_BXdata(ji) * 0.5_dp*(Qi(1)+Qj(1)), p_BYdata(ji) * 0.5_dp*(Qi(1)+Qj(1)) /)
!         case(3)
!           DeltaBi = (/ 0.0_DP, p_BXdata(ij) * (Qj(1)-Qi(1)), p_BYdata(ij) * (Qj(1)-Qi(1)) /)
!           DeltaBj = (/ 0.0_DP, p_BXdata(ji) * (Qi(1)-Qj(1)), p_BYdata(ji) * (Qi(1)-Qj(1)) /)
!         end select

       ! deltaD
       !deltaDi = matmul(Dij-Bij,-deltaQij)
       deltaDi = matmul(Dij,-deltaQij)
       deltaDj = -deltaDi

       ! add deltaK and deltaD to rhs
       ! rhs = rhs + (1-theta)*dt*K*u + (1-theta)*dt*D*u
       do l = 1, nvar2d
          rarrayRhs(l)%Da(i) = rarrayRhs(l)%Da(i) + (1-theta)*dt*deltaKi(l) &
               + (1-theta)*dt*deltaDi(l) 
          rarrayRhs(l)%Da(j) = rarrayRhs(l)%Da(j) + (1-theta)*dt*deltaKj(l) &
               + (1-theta)*dt*deltaDj(l)
       end do
        
        
       

!        ! add deltaK, deltaD and deltaB to rhs
!        ! rhs = rhs + (1-theta)*dt*K*u + (1-theta)*dt*D*u + (1-theta)*dt*S
!        do l = 1, nvar2d
!           rarrayRhs(l)%Da(i) = rarrayRhs(l)%Da(i) + (1-theta)*dt*deltaKi(l) &
!                + (1-theta)*dt*deltaDi(l) !+ (1-theta)*dt*deltaBi(l)
!           rarrayRhs(l)%Da(j) = rarrayRhs(l)%Da(j) + (1-theta)*dt*deltaKj(l) &
!                + (1-theta)*dt*deltaDj(l) !+ (1-theta)*dt*deltaBj(l)
!        end do

    end do


    ! Add the non-conservative part of the bottom source term
!     do i = 1, neq
!       DeltaBi = (/ 0.0_DP, p_BSXdata(i) * Qi(1), p_BSYdata(i) * Qi(1) /)
!       do l = 1, nvar2d
!          rarrayRhs(l)%Da(i) = rarrayRhs(l)%Da(i) + (1-theta)*dt*deltaBi(l)
!        end do
!     end do


    ! If we use the TVD scheme, then add the limited diffusive flux to the RHS
    if (Method == 3) then
       ! Do this for every dimension (dimensional splitting)
       do d = 1, 2

          ! first we fill the array fld2
          do iedge = 1, nedge
             i  = p_Kedge(1, iedge)
             j  = p_Kedge(2, iedge)
             ij = p_Kedge(3, iedge)
             ji = p_Kedge(4, iedge)
             ii = p_Kdiagonal(i)
             jj = p_Kdiagonal(j)

             ! get solution values at node i and j
             Qi = (/rarraySol(1)%Da(i),rarraySol(2)%Da(i),rarraySol(3)%Da(i)/)
             Qj = (/rarraySol(1)%Da(j),rarraySol(2)%Da(j),rarraySol(3)%Da(j)/)

             ! compute deltaQij = Qi - Qj
             deltaQij = Qi - Qj

             ! Compute the Roe-meanvalue for this edge
             Qroeij = calculateQroe(Qi, Qj)

             ! Rij^-1
             invRij = buildInvTrafo(Qroeij,d,gravconst)

             ! compute Rij^-1*(Qj - Qi)
             deltaWij = matmul(invRij,-deltaQij)

             ! compute eigenvalues
             eigenvalues = buildEigenvalues(Qroeij,d,gravconst)

             ! compute scalefactors (a_ij^d and bij^d)
             ! Compute the coeffitients for the jacobi matrices
             JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
             JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
             JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
             JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP

             ! calculate abs(aij^d)
             if (d==1) then
                scalefactor1 = JcoeffxA
                scalefactor2 = JcoeffxB
             else
                scalefactor1 = JcoeffyA
                scalefactor2 = JcoeffyB
             end if


             ! write entries to fld2
             do ivar = 1, nvar2d
                p_fld2(1+4*(ivar-1),iedge) = eigenvalues(ivar)               ! lambda_k
                p_fld2(2+4*(ivar-1),iedge) = deltaWij(ivar)                  ! deltaW^k
                p_fld2(3+4*(ivar-1),iedge) = -eigenvalues(ivar)*scalefactor1  ! k_ij^a
                p_fld2(4+4*(ivar-1),iedge) = -eigenvalues(ivar)*scalefactor2  ! k_ij^b
             end do
          end do


          ! next we fill the array fld1
          ! reset fld1
          call storage_clear (h_fld1)
          ! first compute the Piplus/minus, Qiplus/minus
          do iedge = 1, nedge
             i  = p_Kedge(1, iedge)
             j  = p_Kedge(2, iedge)
             ij = p_Kedge(3, iedge)
             ji = p_Kedge(4, iedge)
             ii = p_Kdiagonal(i)
             jj = p_Kdiagonal(j)


             do ivar = 1, nvar2d
                deltaaplus = max(0.0_DP, p_fld2(3+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
                deltaaminus = min(0.0_DP, p_fld2(3+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
                deltabplus = max(0.0_DP, p_fld2(4+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
                deltabminus = min(0.0_DP, p_fld2(4+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))

                ! write entries to fld1
                if (p_fld2(3+4*(ivar-1),iedge)<0.0_DP) then

                   p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + deltaaplus  ! Pi+
                   p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + deltaaminus ! Pi-

                   p_fld1(3+6*(ivar-1),j) = p_fld1(3+6*(ivar-1),j) + deltaaplus  ! Qj+
                   p_fld1(4+6*(ivar-1),j) = p_fld1(4+6*(ivar-1),j) + deltaaminus ! Qj-
                else
                   p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) + deltaaplus
                   p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) + deltaaminus

                   p_fld1(3+6*(ivar-1),i) = p_fld1(3+6*(ivar-1),i) + deltaaplus
                   p_fld1(4+6*(ivar-1),i) = p_fld1(4+6*(ivar-1),i) + deltaaminus
                end if
                if (p_fld2(4+4*(ivar-1),iedge)<0.0_DP) then

                   p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + deltabplus
                   p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + deltabminus

                   p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) - deltabplus
                   p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) - deltabminus
                else
                   p_fld1(4+6*(ivar-1),j) = p_fld1(4+6*(ivar-1),j) - deltabplus
                   p_fld1(3+6*(ivar-1),j) = p_fld1(3+6*(ivar-1),j) - deltabminus

                   p_fld1(3+6*(ivar-1),i) = p_fld1(3+6*(ivar-1),i) + deltabplus
                   p_fld1(4+6*(ivar-1),i) = p_fld1(4+6*(ivar-1),i) + deltabminus
                end if
             end do
          end do

          ! now compute the Riplus/minus
          do inode = 1, NEQ
             do ivar = 1, nvar2d
                p_fld1(5+6*(ivar-1),inode) = & ! Ri+
                     limiterfunc(p_fld1(3+6*(ivar-1),inode),p_fld1(1+6*(ivar-1),inode),limiter)
                p_fld1(6+6*(ivar-1),inode) = & ! Ri-
                     limiterfunc(p_fld1(4+6*(ivar-1),inode),p_fld1(2+6*(ivar-1),inode),limiter)
             end do
          end do

          ! next update the deltaWij = deltaWij - \hat deltaWij (apply limiting)
          do iedge = 1, nedge
             i  = p_Kedge(1, iedge)
             j  = p_Kedge(2, iedge)
             ij = p_Kedge(3, iedge)
             ji = p_Kedge(4, iedge)
             ii = p_Kdiagonal(i)
             jj = p_Kdiagonal(j)

             do ivar = 1, nvar2d
                if (p_fld2(3+4*(ivar-1),iedge)>0.0_DP) then
                   upwindnode = j
                else
                   upwindnode = i
                end if

                ! deltaak =  kij^a * deltaWij^k
                deltaak = p_fld2(3+4*(ivar-1),iedge) * p_fld2(2+4*(ivar-1),iedge)

                ! deltaWijHat = R(upwind)+/- * deltaWij^k
                if (deltaak<0.0_DP) then
                   deltaWijHat = p_fld1(6+6*(ivar-1),upwindnode)*p_fld2(2+4*(ivar-1),iedge)
                else
                   deltaWijHat = p_fld1(5+6*(ivar-1),upwindnode)*p_fld2(2+4*(ivar-1),iedge)
                end if

                ! deltaWij = deltaWij - \hat deltaWij
                p_fld2(2+4*(ivar-1),iedge) = p_fld2(2+4*(ivar-1),iedge) - deltaWijHat
             end do
          end do

          ! finally compute the fluxes Fij and update the rhs
          do iedge = 1, nedge
             i  = p_Kedge(1, iedge)
             j  = p_Kedge(2, iedge)
             ij = p_Kedge(3, iedge)
             ji = p_Kedge(4, iedge)
             ii = p_Kdiagonal(i)
             jj = p_Kdiagonal(j)

             do ivar = 1, nvar2d
                deltaWij(ivar) = p_fld2(2+4*(ivar-1),iedge)
             end do


             ! get solution values at node i and j
             Qi = (/rarraySol(1)%Da(i),rarraySol(2)%Da(i),rarraySol(3)%Da(i)/)
             Qj = (/rarraySol(1)%Da(j),rarraySol(2)%Da(j),rarraySol(3)%Da(j)/)

             ! Compute the Roe-meanvalue for this edge
             Qroeij = calculateQroe(Qi, Qj)

             ! Compute the coeffitients for the jacobi matrices
             JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
             JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP

             ! calculate abs(aij^d)
             if (d==1) then
                scalefactor = abs(JcoeffxA)
             else
                scalefactor = abs(JcoeffyA)
             end if

             ! compute the (limited) diffusive flux
             deltaFij = scalefactor*matmul(buildTrafo(Qroeij,d,gravconst), &
                  matmul(buildaLambda(Qroeij,d,gravconst),deltaWij))


             ! add the flux F to the rhs
             ! rhs(i) = rhs(i) + (1-theta)*dt*deltaF_ij
             ! rhs(j) = rhs(j) - (1-theta)*dt*deltaF_ij
             do l = 1, nvar2d
                rarrayRhs(l)%Da(i) = rarrayRhs(l)%Da(i) + (1-theta)*dt*deltaFij(l)
                rarrayRhs(l)%Da(j) = rarrayRhs(l)%Da(j) - (1-theta)*dt*deltaFij(l)
             end do

          end do

       end do
    end if
    



  end subroutine BuildShallowWaterRHS











  ! This routine builds the defect for the shallow water system
  subroutine BuildShallowWaterDefect (&
       rdefBlock, rstempBlock, rrhsBlock, rsolBlock, &
       rarrayRhs, rarraySol, rarrayRstemp, &
       p_CXdata, p_CYdata, p_MLdata, rmatrixML, &
       p_BXdata, p_BYdata,  p_BSXdata, p_BSYdata, &
       h_fld1, p_fld1, p_fld2, &
       p_Kdiagonal, p_Kedge, NEQ, nedge, &
       theta, dt, gravconst, Method, limiter)

    ! parameter values
    type(t_array), dimension(nvar2d), intent(INOUT) :: rarrayRstemp
    type(t_array), dimension(nvar2d), intent(IN)    :: rarraySol, rarrayRhs
    type(t_vectorBlock), intent(INOUT)              :: rdefBlock, rstempBlock
    type(t_vectorBlock), intent(IN)                 :: rsolBlock, rrhsBlock
    type(t_matrixScalar), intent(IN)                :: rmatrixML
    real(DP), dimension(:), pointer, intent(IN)     :: p_CXdata, p_CYdata, p_BXdata, p_BYdata, p_MLdata
    real(DP), dimension(:), pointer, intent(IN)     :: p_BSXdata, p_BSYdata
    integer                                         :: h_fld1
    real(DP), dimension(:,:), pointer               :: p_fld1, p_fld2
    integer, dimension(:), pointer, intent(IN)      :: p_Kdiagonal
    integer, dimension(:,:), pointer                :: p_Kedge
    integer, intent(IN)                             :: NEQ,	nedge
    real(DP), intent(IN)                            :: theta, gravconst, dt
    integer, intent(IN)                             :: Method, limiter

    ! variables
    integer                             :: i, j, l, ii, jj, ij, ji, d, iedge, ivar, inode
    real(DP), dimension(nvar2d)         :: deltaQij, Qi, Qj, Qroeij
    real(DP), dimension(nvar2d,nvar2d)  :: JacobixRoeij, JacobiyRoeij
    real(DP)                            :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
    real(DP), dimension(nvar2d,nvar2d)  :: Aij, Bij, Dij
    real(DP)                            :: scalarDissipation, scalefactor, lambda
    real(DP), dimension(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj, deltaBi, deltaBj
    real(DP)                            :: cRoe, uRoe, vRoe
    ! for TVD
    real(DP), dimension(nvar2d,nvar2d)  :: invRij, Rij
    real(DP), dimension(nvar2d)         :: deltaWij, eigenvalues, deltaFij
    real(DP)                            :: scalefactor1, scalefactor2, deltaak, deltaWijHat
    real(DP)                            :: deltaaplus, deltaaminus, deltabplus, deltabminus
    integer                             :: upwindnode



    ! Now we want to compute the defect vector
    ! def = rhs - (ML*Q - theta*dt*K*u - theta*dt*D*u) = rhs - rstemp

    ! calculate: rstemp = ML*Q
    do ivar = 1, nvar2d
       call lsyssc_scalarMatVec (rmatrixML, &
            rsolBlock%Rvectorblock(ivar), &
            rstempBlock%Rvectorblock(ivar), &
            1._DP, 0._DP)
    end do


    ! Now add the fluxes corresponding to the operator K and D,
    ! so we get the low order scheme:
    ! rstemp = rstemp - theta*dt*K*u -theta*dt*D*u
    ! So walk over all edges ij and compute Aij and Bij
    ! Then compute deltaKi=(Aij+Bij)*(Qi-Qj)
    ! and deltaKj=(Aij-Bij)*(Qi-Qj)
    ! deltaDi = (Dij+Bij)*(Qj-Qi)
    ! deltaDj = -deltaDi
    ! and add this to rstemp
    do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)

       ! get solution values at node i and j
       Qi = (/rarraySol(1)%Da(i),rarraySol(2)%Da(i),rarraySol(3)%Da(i)/)
       Qj = (/rarraySol(1)%Da(j),rarraySol(2)%Da(j),rarraySol(3)%Da(j)/)

       ! compute deltaQij = Qi - Qj
       deltaQij = Qi - Qj

       ! Calculate the galerkin fluxes
       ! deltaKi = c_{ij}*(F(Q_i)-F(Q_j))
       ! deltaKj = c_{ji}*(F(Q_j)-F(Q_i))
       deltaKi = p_CXdata(ij)*(buildFlux(Qi,1,gravconst)-buildFlux(Qj,1,gravconst))+ &
                 p_CYdata(ij)*(buildFlux(Qi,2,gravconst)-buildFlux(Qj,2,gravconst))
       deltaKj = p_CXdata(ji)*(buildFlux(Qj,1,gravconst)-buildFlux(Qi,1,gravconst))+ &
                 p_CYdata(ji)*(buildFlux(Qj,2,gravconst)-buildFlux(Qi,2,gravconst))

       ! Now choose the artificial diffusion method
       select case (Method)
        case (0,3)
          ! Here we do not add artificial diffusion - so we have the high order method
          Dij = 0
        case (1,4,6)
          ! Here we use scalar dissipation
          Dij = buildDissipation(Qi,Qj,1,0,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
        case (2,5,7)
          ! Here we use tensorial dissipation
          Dij = abs(p_CXdata(ij))* buildDissipation(Qi,Qj,2,1,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))+&
                abs(p_CYdata(ij))* buildDissipation(Qi,Qj,2,2,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
       end select

!        
!        ! Now we take care of the bottom profile
!        select case (deltabtype)
!         case(0)
!           DeltaBi = 0.0_dp
!           DeltaBj = 0.0_dp
!         case(1)
!           DeltaBi = (/ 0.0_DP, p_BXdata(ij) * Qi(1), p_BYdata(ij) * Qi(1) /)
!           DeltaBj = (/ 0.0_DP, p_BXdata(ji) * Qj(1), p_BYdata(ji) * Qj(1) /)
!         case(2)
!           DeltaBi = (/ 0.0_DP, p_BXdata(ij) * 0.5_dp*(Qi(1)+Qj(1)), p_BYdata(ij) * 0.5_dp*(Qi(1)+Qj(1)) /)
!           DeltaBj = (/ 0.0_DP, p_BXdata(ji) * 0.5_dp*(Qi(1)+Qj(1)), p_BYdata(ji) * 0.5_dp*(Qi(1)+Qj(1)) /)
!         case(3)
!           DeltaBi = (/ 0.0_DP, p_BXdata(ij) * (Qj(1)-Qi(1)), p_BYdata(ij) * (Qj(1)-Qi(1)) /)
!           DeltaBj = (/ 0.0_DP, p_BXdata(ji) * (Qi(1)-Qj(1)), p_BYdata(ji) * (Qi(1)-Qj(1)) /)
!         end select

       ! deltaD
       !deltaDi = matmul(Dij-Bij,-deltaQij)
       deltaDi = matmul(Dij,-deltaQij)
       deltaDj = -deltaDi

       ! add deltaK and deltaD to rstemp
       ! rstemp = rstemp - theta*dt*K*u - theta*dt*D*u
       do l = 1, nvar2d
          rarrayRstemp(l)%Da(i) = rarrayRstemp(l)%Da(i) - theta*dt*deltaKi(l) &
               - theta*dt*deltaDi(l)
          rarrayRstemp(l)%Da(j) = rarrayRstemp(l)%Da(j) - theta*dt*deltaKj(l) &
               - theta*dt*deltaDj(l)
       end do

    end do



!     ! Add the non-conservative part of the bottom source term
!     do i = 1, neq
!       DeltaBi = (/ 0.0_DP, p_BSXdata(i) * Qi(1), p_BSYdata(i) * Qi(1) /)
!       do l = 1, nvar2d
!         rarrayRstemp(l)%Da(i) = rarrayRstemp(l)%Da(i) -theta*dt*deltaBi(l)
!       end do
!    end do


    ! If we use the TVD scheme, then add the limited diffusive flux to the Rstemp (defect)
    if (Method == 3) then
       ! Do this for every dimension (dimensional splitting)
       do d = 1, 2

          ! first we fill the array fld2
          do iedge = 1, nedge
             i  = p_Kedge(1, iedge)
             j  = p_Kedge(2, iedge)
             ij = p_Kedge(3, iedge)
             ji = p_Kedge(4, iedge)
             ii = p_Kdiagonal(i)
             jj = p_Kdiagonal(j)

             ! get solution values at node i and j
             Qi = (/rarraySol(1)%Da(i),rarraySol(2)%Da(i),rarraySol(3)%Da(i)/)
             Qj = (/rarraySol(1)%Da(j),rarraySol(2)%Da(j),rarraySol(3)%Da(j)/)

             ! compute deltaQij = Qi - Qj
             deltaQij = Qi - Qj

             ! Compute the Roe-meanvalue for this edge
             Qroeij = calculateQroe(Qi, Qj)

             ! Rij^-1
             invRij = buildInvTrafo(Qroeij,d,gravconst)

             ! compute Rij^-1*(Qj - Qi)
             deltaWij = matmul(invRij,-deltaQij)

             ! compute eigenvalues
             eigenvalues = buildEigenvalues(Qroeij,d,gravconst)

             ! compute scalefactors (a_ij^d and bij^d)
             ! Compute the coeffitients for the jacobi matrices
             JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
             JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
             JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
             JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP

             ! calculate abs(aij^d)
             if (d==1) then
                scalefactor1 = JcoeffxA
                scalefactor2 = JcoeffxB
             else
                scalefactor1 = JcoeffyA
                scalefactor2 = JcoeffyB
             end if


             ! write entries to fld2
             do ivar = 1, nvar2d
                p_fld2(1+4*(ivar-1),iedge) = eigenvalues(ivar)               ! lambda_k
                p_fld2(2+4*(ivar-1),iedge) = deltaWij(ivar)                  ! deltaW^k
                p_fld2(3+4*(ivar-1),iedge) = -eigenvalues(ivar)*scalefactor1  ! k_ij^a
                p_fld2(4+4*(ivar-1),iedge) = -eigenvalues(ivar)*scalefactor2  ! k_ij^b
             end do
          end do


          ! next we fill the array fld1
          ! reset fld1
          call storage_clear (h_fld1)
          ! first compute the Piplus/minus, Qiplus/minus
          do iedge = 1, nedge
             i  = p_Kedge(1, iedge)
             j  = p_Kedge(2, iedge)
             ij = p_Kedge(3, iedge)
             ji = p_Kedge(4, iedge)
             ii = p_Kdiagonal(i)
             jj = p_Kdiagonal(j)


             do ivar = 1, nvar2d
                deltaaplus = max(0.0_DP, p_fld2(3+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
                deltaaminus = min(0.0_DP, p_fld2(3+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
                deltabplus = max(0.0_DP, p_fld2(4+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
                deltabminus = min(0.0_DP, p_fld2(4+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))

                ! write entries to fld1
                if (p_fld2(3+4*(ivar-1),iedge)<0.0_DP) then

                   p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + deltaaplus  ! Ri+
                   p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + deltaaminus ! Ri-

                   p_fld1(3+6*(ivar-1),j) = p_fld1(3+6*(ivar-1),j) + deltaaplus  ! Qj+
                   p_fld1(4+6*(ivar-1),j) = p_fld1(4+6*(ivar-1),j) + deltaaminus ! Qj-
                else
                   p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) + deltaaplus
                   p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) + deltaaminus

                   p_fld1(3+6*(ivar-1),i) = p_fld1(3+6*(ivar-1),i) + deltaaplus
                   p_fld1(4+6*(ivar-1),i) = p_fld1(4+6*(ivar-1),i) + deltaaminus
                end if
                if (p_fld2(4+4*(ivar-1),iedge)<0.0_DP) then

                   p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + deltabplus
                   p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + deltabminus

                   p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) - deltabplus
                   p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) - deltabminus
                else
                   p_fld1(4+6*(ivar-1),j) = p_fld1(4+6*(ivar-1),j) - deltabplus
                   p_fld1(3+6*(ivar-1),j) = p_fld1(3+6*(ivar-1),j) - deltabminus

                   p_fld1(3+6*(ivar-1),i) = p_fld1(3+6*(ivar-1),i) + deltabplus
                   p_fld1(4+6*(ivar-1),i) = p_fld1(4+6*(ivar-1),i) + deltabminus
                end if
             end do
          end do

          ! now compute the Riplus/minus
          do inode = 1, NEQ
             do ivar = 1, nvar2d
                p_fld1(5+6*(ivar-1),inode) = &
                     limiterfunc(p_fld1(3+6*(ivar-1),inode),p_fld1(1+6*(ivar-1),inode),limiter)
                p_fld1(6+6*(ivar-1),inode) = &
                     limiterfunc(p_fld1(4+6*(ivar-1),inode),p_fld1(2+6*(ivar-1),inode),limiter)
             end do
          end do

          ! next update the deltaWij = deltaWij - \hat deltaWij (apply limiting)
          do iedge = 1, nedge
             i  = p_Kedge(1, iedge)
             j  = p_Kedge(2, iedge)
             ij = p_Kedge(3, iedge)
             ji = p_Kedge(4, iedge)
             ii = p_Kdiagonal(i)
             jj = p_Kdiagonal(j)

             do ivar = 1, nvar2d
                if (p_fld2(3+4*(ivar-1),iedge)>0.0_DP) then
                   upwindnode = j
                else
                   upwindnode = i
                end if

                ! deltaak =  kij^a * deltaWij^k
                deltaak = p_fld2(3+4*(ivar-1),iedge) * p_fld2(2+4*(ivar-1),iedge)

                ! deltaWijHat = R(upwind)+/- * deltaWij^k
                if (deltaak<0.0_DP) then
                   deltaWijHat = p_fld1(6+6*(ivar-1),upwindnode)*p_fld2(2+4*(ivar-1),iedge)
                else
                   deltaWijHat = p_fld1(5+6*(ivar-1),upwindnode)*p_fld2(2+4*(ivar-1),iedge)
                end if

                ! deltaWij = deltaWij - \hat deltaWij
                p_fld2(2+4*(ivar-1),iedge) = p_fld2(2+4*(ivar-1),iedge) - deltaWijHat
             end do
          end do

          ! finally compute the fluxes Fij and update the rhs
          do iedge = 1, nedge
             i  = p_Kedge(1, iedge)
             j  = p_Kedge(2, iedge)
             ij = p_Kedge(3, iedge)
             ji = p_Kedge(4, iedge)
             ii = p_Kdiagonal(i)
             jj = p_Kdiagonal(j)

             do ivar = 1, nvar2d
                deltaWij(ivar) = p_fld2(2+4*(ivar-1),iedge)
             end do


             ! get solution values at node i and j
             Qi = (/rarraySol(1)%Da(i),rarraySol(2)%Da(i),rarraySol(3)%Da(i)/)
             Qj = (/rarraySol(1)%Da(j),rarraySol(2)%Da(j),rarraySol(3)%Da(j)/)

             ! Compute the Roe-meanvalue for this edge
             Qroeij = calculateQroe(Qi, Qj)

             ! Compute the coeffitients for the jacobi matrices
             JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
             JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP

             ! calculate abs(aij^d)
             if (d==1) then
                scalefactor = abs(JcoeffxA)
             else
                scalefactor = abs(JcoeffyA)
             end if

             ! compute the (limited) diffusive flux
             deltaFij = scalefactor*matmul(buildTrafo(Qroeij,d,gravconst), &
                  matmul(buildaLambda(Qroeij,d,gravconst),deltaWij))


             ! add the flux F to the Rstemp (defect)
             ! rstemp(i) = rstemp(i) - theta*dt*deltaF
             ! rstemp(j) = rstemp(j) + theta*dt*deltaF
             do l = 1, nvar2d
                rarrayRstemp(l)%Da(i) = rarrayRstemp(l)%Da(i) - theta*dt*deltaFij(l)
                rarrayRstemp(l)%Da(j) = rarrayRstemp(l)%Da(j) + theta*dt*deltaFij(l)
             end do

          end do

       end do
    end if


    ! Finally we can assemble the defect vector:
    ! def = rhs - rstemp
    call lsysbl_vectorLinearComb (rrhsBlock,rstempBlock, &
         1.0_DP,-1.0_DP,rdefBlock)

  end subroutine BuildShallowWaterDefect





  ! linearised FCT Method for Shallow Water
  ! adding limited antidiffusion to the low order predictor
  ! using characteristic variables and so dimensional splitting
  subroutine linFctShallowWaterAddLimitedAntidiffusion_characteristic(&
       rarraySol, rarraySolDot, rarrayRhs,&
       rdefBlock, rstempBlock, rsolBlock, rSolDotBlock, &
       rmatrixML, p_CXdata, p_CYdata, p_MLdata, p_MCdata, &
       h_fld1, p_fld1, p_fld2, &
       p_Kdiagonal, p_Kedge, NEQ, nedge, &
       gravconst, dt, Method, prelimiting)

    ! parameter values
    type(t_array), dimension(nvar2d), intent(INOUT) :: rarraySol, rarraySolDot
    type(t_array), dimension(nvar2d), intent(IN)    :: rarrayRhs
    type(t_vectorBlock), intent(INOUT)              :: rdefBlock, rstempBlock
    type(t_vectorBlock), intent(IN)                 :: rsolBlock
    type(t_vectorBlock), intent(INOUT)              :: rSolDotBlock
    type(t_matrixScalar), intent(IN)                :: rmatrixML
    real(DP), dimension(:), pointer, intent(IN)     :: p_CXdata, p_CYdata, p_MLdata, p_MCdata
    integer                                         :: h_fld1
    real(DP), dimension(:,:), pointer               :: p_fld1, p_fld2
    integer, dimension(:), pointer, intent(IN)      :: p_Kdiagonal
    integer, dimension(:,:), pointer                :: p_Kedge
    integer, intent(IN)                             :: NEQ,	nedge
    real(DP), intent(IN)                            :: gravconst, dt
    integer, intent(IN)                             :: Method, prelimiting


    ! variables
    integer                             :: i, j, l, ii, jj, ij, ji, d, iedge, ivar
    real(DP), dimension(nvar2d)         :: deltaQij, Qi, Qj, Qroeij
    real(DP), dimension(nvar2d,nvar2d)  :: JacobixRoeij, JacobiyRoeij
    real(DP)                            :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
    real(DP), dimension(nvar2d,nvar2d)  :: Aij, Bij, Dij
    real(DP)                            :: scalarDissipation, scalefactor, lambda
    real(DP), dimension(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj
    real(DP)                            :: cRoe, uRoe, vRoe
    real(DP), dimension(nvar2d)         :: Udoti, Udotj
    real(DP), dimension(nvar2d)         :: deltaWij, deltaGij, deltaFij
    real(DP), dimension(nvar2d,nvar2d)  :: Rij, invRij



    ! For the following use dimensional splitting
    do d = 1, 2

       ! Clear SolDot
       call lsysbl_clearVector (rSolDotBlock)

       ! First compute SolDot
       do iedge = 1, nedge
          i  = p_Kedge(1, iedge)
          j  = p_Kedge(2, iedge)
          ij = p_Kedge(3, iedge)
          ji = p_Kedge(4, iedge)
          ii = p_Kdiagonal(i)
          jj = p_Kdiagonal(j)

          ! get solution values at node i and j
          do ivar = 1, nvar2d
             Qi(ivar) = rarraySol(ivar)%Da(i)
             Qj(ivar) = rarraySol(ivar)%Da(j)
          end do

          ! compute deltaQij = Qi - Qj
          deltaQij = Qi - Qj

          ! Compute the Roe-meanvalue for this edge
          Qroeij = calculateQroe(Qi, Qj)

          if (d==1) then
             scalefactor = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
          else
             scalefactor = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
          end if

          ! Now we can compute Aij and Bij
          Aij = scalefactor*buildJacobi(Qroeij,d,gravconst)

          ! deltaK
          deltaKi = -matmul(Aij,deltaQij)
          deltaKj = deltaKi

          ! Calculate this alternatively by calculating
          ! deltaKi = c_{ij}*(F(Q_i)-F(Q_j))
          ! deltaKj = c_{ji}*(F(Q_j)-F(Q_i))
          if (d==1) then
             deltaKi = p_CXdata(ij)*buildFlux(Qi,1,gravconst)&
                  -p_CXdata(ij)*buildFlux(Qj,1,gravconst)
             deltaKj = p_CXdata(ji)*buildFlux(Qj,1,gravconst)&
                  -p_CXdata(ji)*buildFlux(Qi,1,gravconst)
          else
             deltaKi = p_CYdata(ij)*buildFlux(Qi,2,gravconst) &
                  -p_CYdata(ij)*buildFlux(Qj,2,gravconst)
             deltaKj = p_CYdata(ji)*buildFlux(Qj,2,gravconst) &
                  -p_CYdata(ji)*buildFlux(Qi,2,gravconst)
          end if


          ! compute Dij
          if (Method == 4) then
             ! Build Dij using scalar dissipation
             
             
             
       select case (scalardisstype)
       case (1)
       scalarDissipation = maxval(abs(buildEigenvalues(Qroeij,d,gravconst)))
       case (2)
       ! Dissipation of Rusanov type
       scalardissipation = &
             max(abs(p_CXdata(ij)*Qj(2)/Qj(1)) + sqrt(gravconst*Qj(1)),&
                 abs(p_Cydata(ji)*Qi(3)/Qi(1)) + sqrt(gravconst*Qi(1)) )
       end select
             
             
             
             Dij = abs(scalefactor)*scalarDissipation*Eye
          else if (Method == 5) then
             ! Build Dij using tensorial dissipation
             Dij = abs(scalefactor)*matmul(buildTrafo(Qroeij,d,gravconst),&
                  matmul(buildaLambda(Qroeij,d,gravconst),&
                  buildinvTrafo(Qroeij,d,gravconst)))
          end if

          ! deltaD
          deltaDi = matmul(Dij,-deltaQij)
          deltaDj = -deltaDi

          ! add deltaK and deltaD to SolDot
          ! SolDot = SolDot + 1/ML*(K*u + D*u)
          do l = 1, nvar2d
             rarraySolDot(l)%Da(i) = rarraySolDot(l)%Da(i) &
                  + 1.0_DP/p_MLdata(i)*&
                  (deltaKi(l) + deltaDi(l))
             rarraySolDot(l)%Da(j) = rarraySolDot(l)%Da(j) &
                  + 1.0_DP/p_MLdata(j)*&
                  (deltaKj(l) + deltaDj(l))
             ! Save Qroe
             p_fld2(1+2*(l-1),iedge) = QRoeij(l)
          end do

       end do
       ! SolDot fertig!



       ! Now calculate the limited antidiffusive flux

       call storage_clear (h_fld1)
       ! first we fill the array fld1 + fld2 (that means compute deltaWij and deltaGij and save it)
       do iedge = 1, nedge
          i  = p_Kedge(1, iedge)
          j  = p_Kedge(2, iedge)
          ij = p_Kedge(3, iedge)
          ji = p_Kedge(4, iedge)
          ii = p_Kdiagonal(i)
          jj = p_Kdiagonal(j)

          ! get solution values at node i and j
          do ivar = 1, nvar2d
             Qi(ivar) = rarraySol(ivar)%Da(i)
             Qj(ivar) = rarraySol(ivar)%Da(j)
          end do

          ! compute deltaQij = Qi - Qj
          deltaQij = Qi - Qj

          ! Compute the Roe-meanvalue for this edge
          Qroeij = calculateQroe(Qi, Qj)

          ! Rij^-1
          invRij = buildInvTrafo(Qroeij,d,gravconst)

          ! compute Rij^-1*(Qj - Qi)
          deltaWij = -matmul(invRij,deltaQij)

          if (d==1) then
             scalefactor = (p_CXdata(ij)-p_CXdata(ji))/2._DP
          else
             scalefactor = (p_CYdata(ij)-p_CYdata(ji))/2._DP
          end if

          ! compute Dij
          if (Method == 4) then
             ! Build Dij using scalar dissipation
       select case (scalardisstype)
       case (1)
       scalarDissipation = maxval(abs(buildEigenvalues(Qroeij,d,gravconst)))
       case (2)
       ! Dissipation of Rusanov type
       scalardissipation = &
             max(abs(p_CXdata(ij)*Qj(2)/Qj(1)) + sqrt(gravconst*Qj(1)),&
                 abs(p_Cydata(ji)*Qi(3)/Qi(1)) + sqrt(gravconst*Qi(1)) )
       end select
             Dij = abs(scalefactor)*scalarDissipation*Eye
          else if (Method == 5) then
             ! Build Dij using tensorial dissipation
             Dij = abs(scalefactor)*matmul(buildTrafo(Qroeij,d,gravconst),&
                  matmul(buildaLambda(Qroeij,d,gravconst),&
                  buildinvTrafo(Qroeij,d,gravconst)))
          end if

          ! get Udot at node i and j
          do ivar = 1, nvar2d
             Udoti(ivar) = rarraySolDot(ivar)%Da(i)
             Udotj(ivar) = rarraySolDot(ivar)%Da(j)
          end do

          ! compute deltaFij (73)
          deltaFij = p_MCdata(ij)*(Udoti-Udotj)+matmul(Dij,deltaQij)

          ! Now apply prelimiting
          if (prelimiting == 1) then
             if (Method==4) then
                ! For scalar dissipation apply
                ! MinMod prelimiting
                do ivar = 1, nvar2d
                   if (deltaFij(ivar)*deltaQij(ivar)<0) then
                      deltaFij(ivar)=0
                   else
                      if (abs(deltaFij(ivar))>abs(abs(scalefactor)*scalarDissipation*deltaQij(ivar))) then
                         deltaFij(ivar) = abs(scalefactor)*scalarDissipation*deltaQij(ivar)
                      end if
                   end if
                end do
             else if (Method==5) then
                ! For tensorial dissipation apply
                ! Simple prelimiting
                do ivar = 1, nvar2d
                   if (deltaFij(ivar)*deltaQij(ivar)<0) then
                      deltaFij(ivar) = 0
                   end if
                end do
             end if
          end if


          ! compute deltaGij
          deltaGij = matmul(invRij,deltaFij)

          do ivar = 1, nvar2d

             ! Save deltaGij
             p_fld2(2+2*(ivar-1),iedge) = deltaGij(ivar)

             ! Update the P/Q_i/j +/- (76+77)
             p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + max(0.0_DP,deltaGij(ivar)) ! Pi+
             p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + min(0.0_DP,deltaGij(ivar)) ! Pi-
             p_fld1(3+6*(ivar-1),i) = max(p_fld1(3+6*(ivar-1),i),deltaWij(ivar))! Qi+
             p_fld1(4+6*(ivar-1),i) = min(p_fld1(4+6*(ivar-1),i),deltaWij(ivar))! Qi-

             p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) + max(0.0_DP,-deltaGij(ivar)) ! Pj+
             p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) + min(0.0_DP,-deltaGij(ivar)) ! Pj-
             p_fld1(3+6*(ivar-1),j) = max(p_fld1(3+6*(ivar-1),j),-deltaWij(ivar))! Qj+
             p_fld1(4+6*(ivar-1),j) = min(p_fld1(4+6*(ivar-1),j),-deltaWij(ivar))! Qj-

          end do

       end do

       do i = 1, NEQ

          do ivar = 1, nvar2d
             ! Compute the R_i +/- (78)
             if (abs(p_fld1(1+6*(ivar-1),i))>1e-8) then
                p_fld1(5+6*(ivar-1),i) = min(1.0_DP,p_MLdata(i)*&
                     p_fld1(3+6*(ivar-1),i)/p_fld1(1+6*(ivar-1),i)/dt)! Ri+
             else
                p_fld1(5+6*(ivar-1),i) = 0.0_DP
             end if

             if (abs(p_fld1(2+6*(ivar-1),i))>1e-8) then
                p_fld1(6+6*(ivar-1),i) = min(1.0_DP, p_MLdata(i)*&
                     p_fld1(4+6*(ivar-1),i)/p_fld1(2+6*(ivar-1),i)/dt)! Ri-
             else
                p_fld1(6+6*(ivar-1),i) = 0.0_DP
             end if

          end do



       end do

       do iedge = 1, nedge
          i  = p_Kedge(1, iedge)
          j  = p_Kedge(2, iedge)
          ij = p_Kedge(3, iedge)
          ji = p_Kedge(4, iedge)
          ii = p_Kdiagonal(i)
          jj = p_Kdiagonal(j)

          do ivar = 1, nvar2d
             deltaGij(ivar) = p_fld2(2+2*(ivar-1),iedge)

             ! Limit the antidiffusive fluxes (79)
             if (deltaGij(ivar)>0.0_DP) then
                deltaGij(ivar) = deltaGij(ivar) * min(p_fld1(5+6*(ivar-1),i),p_fld1(6+6*(ivar-1),j))
             else
                deltaGij(ivar) = deltaGij(ivar) * min(p_fld1(5+6*(ivar-1),j),p_fld1(6+6*(ivar-1),i))
             end if

          end do

          ! get Roe-values at ij
          do ivar = 1, nvar2d
             Qroeij(ivar) = p_fld2(1+2*(ivar-1),iedge)
          end do

          ! Rij
          Rij = buildTrafo(Qroeij,d,gravconst)

          ! Finally we can transform back to deltaFij
          deltaFij = matmul(Rij,deltaGij)

          ! And save 
          do ivar = 1, nvar2d
             p_fld2(2+2*(ivar-1),iedge) = deltaFij(ivar)
          end do

       end do



       do iedge = 1, nedge
          i  = p_Kedge(1, iedge)
          j  = p_Kedge(2, iedge)
          ij = p_Kedge(3, iedge)
          ji = p_Kedge(4, iedge)
          ii = p_Kdiagonal(i)
          jj = p_Kdiagonal(j)

          do ivar = 1, nvar2d
             deltaFij(ivar) = p_fld2(2+2*(ivar-1),iedge)

             rarraySol(ivar)%Da(i) = rarraySol(ivar)%Da(i) + deltaFij(ivar)   *dt/p_MLdata(i)
             rarraySol(ivar)%Da(j) = rarraySol(ivar)%Da(j) - deltaFij(ivar)   *dt/p_MLdata(j)
          end do

       end do


    end do ! d=1, 2

  end subroutine linFctShallowWaterAddLimitedAntidiffusion_characteristic










  ! This subroutine takes care of boundary conditions
  ! for a reflacting boundary
  subroutine ImplementShallowWaterBCs (&
       rboundary, rtriangulation, &
       rarrayP, rarraySol, rarrayDef, &
       p_Kdiagonal, p_Kld, p_Kcol, &
       gravconst, boundarycorner)

    ! PARAMETER VALUES
    ! the boundary
    type(t_boundary), intent(IN) :: rboundary
    ! the triangulation
    type(t_triangulation), intent(IN) :: rtriangulation
    ! pointers to datas of preconditioner, solution and defect
    type(t_array), dimension(nvar2d), intent(INOUT) :: rarrayP, rarraySol, rarrayDef
    ! pointer to index of diagonal positions of matrices
    integer, dimension(:), pointer, intent(IN) :: p_Kdiagonal, p_Kld, p_Kcol
    real(DP), intent(IN) :: gravconst
    integer :: boundarycorner



    ! VARIABLES
    ! index and number of boundary components
    integer :: ibct, nbct
    ! Index of first and last node of this boundarycomponent
    integer :: ivbdFirst, ivbdLast
    ! Pointers to the triangultion datas
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    ! indices of boundary nodes
    integer :: ivbd
    ! the index of the boundary node to edit
    integer :: i
    ! index of diagonal element ii in matrix
    integer :: ii
    ! actual column
    integer :: icol
    ! parametervalue of boundary node to edit
    real(DP) :: parmv
    ! unit outward normal vector and tangential vector
    real(DP), dimension(2) :: normalV, tangentialV
    ! old primitive values at boundary node
    real(DP) :: h, u, v
    ! predicted primitive star-values at boundary node
    real(DP) :: hstar, ustar, vstar
    ! final primitive starstar-values at boundary node
    real(DP) :: hstarstar, ustarstar, vstarstar, hm, oldh
    ! normal and tangential projection of the predicted velocity
    real(DP) :: vnproj, vtproj
    ! for the fixpointiteration to calculate h**
    integer :: ite
    ! maximum iterations
    integer, parameter :: itemax = 50
    ! Element ij in matrix (preconditioner)
    integer :: ij
    ! explicit update for predicted solution
    real(DP), dimension(3) :: up
    ! Temp values for newton step
    real(dp) :: temp1, temp2
    ! actual variable
    integer :: ivar



    ! Initialise some values

    temp1 = sqrt(2*gravconst)
    temp2 = sqrt(gravconst/2)

    ! get number of boundary components
    nbct = rtriangulation%NBCT

    ! get pointer to IboundaryCpIdx, which saves where the boundarycomponents begin and end
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)

    ! get pointer to IverticesAtBoundary, which saves die numbers of the vertices on the boundary
    call storage_getbase_int (rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)

    ! get pointer to DvertexParameterValue, which saves the parametervalues 
    call storage_getbase_double (rtriangulation%h_DvertexParameterValue, p_DvertexParameterValue)


    ! Now start manipulation of boundary values, preconditioner and defect

    ! loop over all boundarycomponents
    do ibct = 1, nbct

       ! get Index of first and last node of this boundarycomponent
       ivbdFirst = p_IboundaryCpIdx(ibct)
       ivbdLast  = p_IboundaryCpIdx(ibct+1)-1


       do ivbd = ivbdFirst, ivbdLast
          ! the index of the boundary node to edit
          i = p_IverticesAtBoundary(ivbd)
          ! entry of diagonal element in matrix
          ii = p_Kdiagonal(i)
          ! parameter value of the boundary node to edit
          parmv = p_DvertexParameterValue(ivbd)

          ! Shall the nodes in a corner be manipulated, too?
          if (((abs(parmv-aint(PARMV))>1e-6)).or.(boundarycorner==1)) then


             ! get normal vector at boundary node
             call boundary_getNormalVec2D(rboundary,ibct,parmv,normalV(1),normalV(2))
             if (abs(parmv-aint(PARMV))<1e-6) then ! if we are at a corner
                ! get normalvector left and right from the corner
                call boundary_getNormalVec2D(rboundary,ibct,parmv+1e-6,normalV(1),normalV(2))
                call boundary_getNormalVec2D(rboundary,ibct,parmv-1e-6,tangentialv(1),tangentialv(2))
                ! add those normalvectors
                normalV = normalV + tangentialv
                ! and make it unit length
                normalv = normalv /sqrt(normalV(1)**2.0_DP + normalV(2)**2.0_DP)
             end if

             ! tangential Vector
             tangentialv(1) = -normalV(2)
             tangentialv(2) =  normalV(1)

!             ! Calculate the predicted values
!             up = 0
!             do ij = p_kld(i), p_kld(i+1)-1
!                if (ij .ne. ii) then
!                  icol = p_Kcol(ij)
!                  do ivar = 1, nvar2d
!                    up(ivar)= up(ivar) - rarrayP(ivar)%Da(ij) * rarrayDef(ivar)%Da(icol)
!                  end do
!                else
!                  do ivar = 1, nvar2d
!                    up(ivar)= up(ivar) + rarrayDef(ivar)%Da(i)
!                  end do
!                end if
!             end do
!             
!             do ivar = 1, nvar2d
!               rarraySol(ivar)%Da(i) = rarraySol(ivar)%Da(i) + up(ivar)/rarrayP(ivar)%Da(ii)
!             end do
                          

             ! update preconditioner
             do ij = p_kld(i), p_kld(i+1)-1
                if (ij .ne. ii) then
                   rarrayP(1)%Da(ij) = 0.0_dp
                   rarrayP(2)%Da(ij) = 0.0_dp
                   rarrayP(3)%Da(ij) = 0.0_dp
                end if
             end do


             ! predict values (AFC II (90))
             rarraySol(1)%Da(i) = rarraySol(1)%Da(i) + rarrayDef(1)%Da(i)/rarrayP(1)%Da(ii)
             rarraySol(2)%Da(i) = rarraySol(2)%Da(i) + rarrayDef(2)%Da(i)/rarrayP(2)%Da(ii)
             rarraySol(3)%Da(i) = rarraySol(3)%Da(i) + rarrayDef(3)%Da(i)/rarrayP(3)%Da(ii)
             

             ! get predicted primitive values
             hstar = rarraySol(1)%Da(i)
             if (hstar<1e-6) then
               hstar=0
               ustar=0
               vstar=0
             else
               ustar = rarraySol(2)%Da(i)/hstar
               vstar = rarraySol(3)%Da(i)/hstar
             end if

             ! update defect
             rarrayDef(1)%Da(i) = 0.0_dp
             rarrayDef(2)%Da(i) = 0.0_dp
             rarrayDef(3)%Da(i) = 0.0_dp


             ! solve the riemannian problem

             ! projection of velocity on outward normal vector and tangential vector
             vnproj = ustar * normalV(1) + vstar * normalV(2)
             vtproj = ustar * tangentialV(1) + vstar * tangentialV(2)

             ! Test, if dry bed is generated
             if (4.0_dp*(sqrt(hstar*gravconst))<-2.0_dp*vnproj) then
                !write (*,*) 'ERROR: Boundary Riemann Problem - Physical drybed generation! Wetbet Riemann Solver not applicable!'
                !write (*,*) hstar,vnproj
                hstarstar = 0
             else
             ! Bed should stay wet
             ! Initialisiere hstarstar
             ! entweder mit der predicted solution
             !hstarstar = hstar
             ! oder mit dem two-rarefaction ansatz (der explizit ausgerechnet werden kann)
             hstarstar = ((SQRT(gravconst*hstar)+0.5*vnproj)**2.0_dp)/gravconst

             iteh: do ite = 1, itemax
	        ! get h** as solution of the riemannian problem
	        oldh = hstarstar
	        if (hstarstar .le. hstar) then
            ! h** can be computed explicitely (2 rarefaction)
                   hstarstar = ((sqrt(gravconst*hstar)+0.5*vnproj)**2.0_dp)/gravconst
	        else
            ! 2 Shock
                   hm = hstarstar

                   ! compute h** by fixpointiteration
                   !hm = hstar + vnproj/sqrt(0.5*gravconst*(hm+hstar)/(hm*hstar))
                   !hm = hstar + vnproj/sqrt(0.5*gravconst*(hm+hstar)/(hm*hstar))
                   !hm = hstar + vnproj/sqrt(0.5*gravconst*(hm+hstar)/(hm*hstar))

                   ! or by newton method
                   hm = hm - (2*(hm-hstar)*sqrt(gravconst/2*(hm+hstar)/(hm*hstar))-2*vnproj) & ! !!!!! oder +2*vnproj? !!!!!
                        /(temp1*sqrt((hm+hstar)/(hm*hstar))+temp2*(hm-hstar)*(1/(hm*hstar)&
                          -(hm+hstar)/(hm*hm*hstar))/sqrt((hm+hstar)/(hm*hstar)))

                   hstarstar = hm

	        end if
         ! Test if the algorithm converged
	        if (abs((hstarstar - oldh)/(hstarstar+oldh))< 1e-8) then
                   exit iteh
	        end if
             end do iteh
end if !dry or wet bed

             ! Test for convergence and give out warnings
             if (ite==50) then
                write (*,*) 'ERROR! Boundary condition riemann problem did not converge'
             end if
             if (hstarstar<1e-6) then
                ! Clip the water heights if too small
                !write (*,*) 'ERROR! Boundary condition riemann problem h very small or smaller than zero'
                hstarstar=0
             end if
             if (ite>10) then
                write (*,*) 'WARNING! Boundary condition riemann problem convergence is slow'
             end if


             ! compute u** and v** as projection of u* and v* on the tangential Vector
             ustarstar = vtproj * tangentialV(1)
             vstarstar = vtproj * tangentialV(2)


             ! Finally transform back to primitive variables and update the solution
             rarraySol(1)%Da(i) = hstarstar
             rarraySol(2)%Da(i) = hstarstar*ustarstar
             rarraySol(3)%Da(i) = hstarstar*vstarstar


          end if


       end do
    end do

  end subroutine ImplementShallowWaterBCs





  ! linearised FCT Method for Shallow Water
  ! adding limited antidiffusion to the low order predictor
  ! using conservative variables variables and a syncronized limiting factor
  subroutine linFctShallowWaterAddLimitedAntidiffusion_syncronized(&
       rarraySol, rarraySolDot, rarrayRhs,&
       rdefBlock, rstempBlock, rsolBlock, rSolDotBlock, &
       rmatrixML, p_CXdata, p_CYdata, p_MLdata, p_MCdata, &
       h_fld1, p_fld1, p_fld2, &
       p_Kdiagonal, p_Kedge, NEQ, nedge, &
       gravconst, dt, Method, prelimiting, syncromethod, &
       rtriangulation)

    ! parameter values
    type(t_array), dimension(nvar2d), intent(INOUT) :: rarraySol, rarraySolDot
    type(t_array), dimension(nvar2d), intent(IN)    :: rarrayRhs
    type(t_vectorBlock), intent(INOUT)              :: rdefBlock, rstempBlock
    type(t_vectorBlock), intent(IN)                 :: rsolBlock
    type(t_vectorBlock), intent(INOUT)              :: rSolDotBlock
    type(t_matrixScalar), intent(IN)                :: rmatrixML
    real(DP), dimension(:), pointer, intent(IN)     :: p_CXdata, p_CYdata, p_MLdata, p_MCdata
    integer                                         :: h_fld1
    real(DP), dimension(:,:), pointer               :: p_fld1, p_fld2
    integer, dimension(:), pointer, intent(IN)      :: p_Kdiagonal
    integer, dimension(:,:), pointer                :: p_Kedge
    integer, intent(IN)                             :: NEQ,	nedge
    real(DP), intent(IN)                            :: gravconst, dt
    integer, intent(IN)                             :: Method, prelimiting, syncromethod
    type(t_triangulation), intent(IN)               :: rtriangulation


    ! variables
    integer                             :: i, j, l, ii, jj, ij, ji, d, iedge, ivar
    real(DP), dimension(nvar2d)         :: deltaQij, Qi, Qj, Qroeij
    real(DP), dimension(nvar2d,nvar2d)  :: JacobixRoeij, JacobiyRoeij
    real(DP)                            :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
    real(DP), dimension(nvar2d,nvar2d)  :: Aij, Bij, Dij, Tij, invTij
    real(DP)                            :: scalarDissipation, scalefactor, lambda
    real(DP), dimension(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj
    real(DP)                            :: cRoe, uRoe, vRoe
    real(DP), dimension(nvar2d)         :: Udoti, Udotj
    real(DP), dimension(nvar2d)         :: deltaWij, deltaGij, deltaFij
    real(DP), dimension(nvar2d,nvar2d)  :: Rij, invRij
    real(DP)                            :: alphaij
    integer                             :: ibct, nbct, ivbd, ivbdFirst, ivbdLast
    integer, dimension(:), pointer      :: p_IboundaryCpIdx
    integer, dimension(:), pointer      :: p_IverticesAtBoundary
    

    ! Clear SolDot
    call lsysbl_clearVector (rSolDotBlock)

    ! First compute SolDot
    do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)

       ! get solution values at node i and j
       do ivar = 1, nvar2d
          Qi(ivar) = rarraySol(ivar)%Da(i)
          Qj(ivar) = rarraySol(ivar)%Da(j)
       end do

       ! compute deltaQij = Qi - Qj
       deltaQij = Qi - Qj

       ! Compute the Roe-meanvalue for this edge
       Qroeij = calculateQroe(Qi, Qj)


       ! Now we can compute Aij and Bij
       !Aij = scalefactor*buildJacobi(Qroeij,d,gravconst)

       ! deltaK
       !deltaKi = -matmul(Aij,deltaQij)
       !deltaKj = deltaKi

       ! Calculate this alternatively by calculating
       ! deltaKi = c_{ij}*(F(Q_i)-F(Q_j))
       ! deltaKj = c_{ji}*(F(Q_j)-F(Q_i))
       deltaKi = p_CXdata(ij)*buildFlux(Qi,1,gravconst)+p_CYdata(ij)*buildFlux(Qi,2,gravconst) &
            -p_CXdata(ij)*buildFlux(Qj,1,gravconst)-p_CYdata(ij)*buildFlux(Qj,2,gravconst)
       deltaKj = p_CXdata(ji)*buildFlux(Qj,1,gravconst)+p_CYdata(ji)*buildFlux(Qj,2,gravconst) &
            -p_CXdata(ji)*buildFlux(Qi,1,gravconst)-p_CYdata(ji)*buildFlux(Qi,2,gravconst)


       ! compute Dij
       ! Now choose the artificial diffusion method
       select case (Method)
        case (6)
          ! Here we use scalar dissipation
          Dij = buildDissipation(Qi,Qj,1,0,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
        case (7)
          ! Here we use tensorial dissipation
          Dij = abs(p_CXdata(ij))* buildDissipation(Qi,Qj,2,1,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))+&
                abs(p_CYdata(ij))* buildDissipation(Qi,Qj,2,2,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
       end select

       ! deltaD
       deltaDi = matmul(Dij,-deltaQij)
       deltaDj = -deltaDi

       ! add deltaK and deltaD to SolDot
       ! SolDot = SolDot + 1/ML*(K*u + D*u)
       do l = 1, nvar2d
          rarraySolDot(l)%Da(i) = rarraySolDot(l)%Da(i) &
               + 1.0_DP/p_MLdata(i)*&
               (deltaKi(l) + deltaDi(l))
          rarraySolDot(l)%Da(j) = rarraySolDot(l)%Da(j) &
               + 1.0_DP/p_MLdata(j)*&
               (deltaKj(l) + deltaDj(l))
          ! Save Qroe
          p_fld2(1+2*(l-1),iedge) = QRoeij(l)
       end do

    end do
    ! SolDot fertig!


    ! Now calculate the limited antidiffusive flux

    call storage_clear (h_fld1)
    ! first we fill the array fld1 + fld2 (that means compute deltaWij and deltaGij and save it)
    do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)

       ! get solution values at node i and j
       do ivar = 1, nvar2d
          Qi(ivar) = rarraySol(ivar)%Da(i)
          Qj(ivar) = rarraySol(ivar)%Da(j)
       end do

       ! compute deltaQij = Qi - Qj
       deltaQij = Qi - Qj

       ! Compute the Roe-meanvalue for this edge
       Qroeij = calculateQroe(Qi, Qj)

       ! Trafomatrix Tij
       Tij = Eye

       ! compute Tij*(Qj - Qi), the transformed solution differences
       deltaWij = matmul(Tij,-deltaQij)

       ! Compute the coeffitients for the jacobi matrices
       JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
       JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
       JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
       JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP

       ! compute Dij
       ! Now choose the artificial diffusion method
       select case (Method)
        case (6)
          ! Here we use scalar dissipation
          Dij = buildDissipation(Qi,Qj,1,0,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
        case (7)
          ! Here we use tensorial dissipation
          Dij = abs(p_CXdata(ij))* buildDissipation(Qi,Qj,2,1,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))+&
                abs(p_CYdata(ij))* buildDissipation(Qi,Qj,2,2,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
       end select

       ! get Udot at node i and j
       do ivar = 1, nvar2d
          Udoti(ivar) = rarraySolDot(ivar)%Da(i)
          Udotj(ivar) = rarraySolDot(ivar)%Da(j)
       end do

       ! compute deltaFij (73)/(59) linearised fct
       !deltaFij = p_MCdata(ij)*(Udoti-Udotj)+matmul(Dij,deltaQij)
       deltaFij = matmul(Dij,deltaQij)

       ! compute deltaGij, the transformed antidiffusive fluxes
       deltaGij = matmul(Tij,deltaFij)
       
       ! Now apply prelimiting, siehe AFCII(70)
       if (prelimiting == 1) then
          do ivar = 1, nvar2d
             if (deltaGij(ivar)*deltaWij(ivar).ge.0) deltaGij(ivar)=0
          end do
       end if

       do ivar = 1, nvar2d

          ! Save deltaGij
          p_fld2(2+2*(ivar-1),iedge) = deltaGij(ivar)

          ! Update the P/Q_i/j +/- (76+77)
          p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + max(0.0_DP,deltaGij(ivar)) ! Pi+
          p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + min(0.0_DP,deltaGij(ivar)) ! Pi-
          p_fld1(3+6*(ivar-1),i) = max(p_fld1(3+6*(ivar-1),i),deltaWij(ivar))! Qi+
          p_fld1(4+6*(ivar-1),i) = min(p_fld1(4+6*(ivar-1),i),deltaWij(ivar))! Qi-

          p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) + max(0.0_DP,-deltaGij(ivar)) ! Pj+
          p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) + min(0.0_DP,-deltaGij(ivar)) ! Pj-
          p_fld1(3+6*(ivar-1),j) = max(p_fld1(3+6*(ivar-1),j),-deltaWij(ivar))! Qj+
          p_fld1(4+6*(ivar-1),j) = min(p_fld1(4+6*(ivar-1),j),-deltaWij(ivar))! Qj-

       end do

    end do

    do i = 1, NEQ

       do ivar = 1, nvar2d
          ! Compute the R_i +/- (73) AFC II
!           if (p_fld1(1+6*(ivar-1),i)>1e-3) then
!              p_fld1(5+6*(ivar-1),i) = min(1.0_DP,p_MLdata(i)*&
!                   p_fld1(3+6*(ivar-1),i)/p_fld1(1+6*(ivar-1),i)/dt)! Ri+
!           else
!              p_fld1(5+6*(ivar-1),i) = 0.0_DP
!           end if
! 
!           if (p_fld1(2+6*(ivar-1),i)<-1e-3) then
!              p_fld1(6+6*(ivar-1),i) = min(1.0_DP, p_MLdata(i)*&
!                   p_fld1(4+6*(ivar-1),i)/p_fld1(2+6*(ivar-1),i)/dt)! Ri-
!           else
!              p_fld1(6+6*(ivar-1),i) = 0.0_DP
!           end if
          
          
          
          
                  p_fld1(5+6*(ivar-1),i) = min(1.0_DP,p_MLdata(i)/dt*&
                  p_fld1(3+6*(ivar-1),i)/(p_fld1(1+6*(ivar-1),i)+SYS_EPSREAL_DP))! Ri+
                  
                  p_fld1(6+6*(ivar-1),i) = min(1.0_DP, p_MLdata(i)/dt*&
                  p_fld1(4+6*(ivar-1),i)/(p_fld1(2+6*(ivar-1),i)-SYS_EPSREAL_DP))! Ri-
          
          



          
!          ! If we are at a maximum (Qi+ = 0) then cancel negative flux (set Ri- = 0) to avoid clipping
!          if (p_fld1(3+6*(ivar-1),i)<1e-12) then
!            p_fld1(6+6*(ivar-1),i) = 0
!          end if
!          
!          ! If we are at a minimum (Qi- = 0) then cancel positive flux (set Ri+ = 0) to avoid clipping
!          if (p_fld1(4+6*(ivar-1),i)>1e-12) then
!            p_fld1(5+6*(ivar-1),i) = 0
!          end if

       end do

    end do
    
    
    
    
!!!!!! New Idea: Cancel antidiffusive fluxes at the boundary:
!!!!!! Set Ri+/- = 0 at boundary nodes
!    
!    ! get number of boundary components
!    nbct = rtriangulation%NBCT
!
!    ! get pointer to IboundaryCpIdx, which saves where the boundarycomponents begin and end
!    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
!
!    ! get pointer to IverticesAtBoundary, which saves die numbers of the vertices on the boundary
!    call storage_getbase_int (rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
!
!
!    ! loop over all boundarycomponents
!    do ibct = 1, nbct
!
!       ! get Index of first and last node of this boundarycomponent
!       ivbdFirst = p_IboundaryCpIdx(ibct)
!       ivbdLast  = p_IboundaryCpIdx(ibct+1)-1
!
!       do ivbd = ivbdFirst, ivbdLast
!          ! the index of the boundary node to edit
!          i = p_IverticesAtBoundary(ivbd)
!          p_fld1(5+6*(ivar-1),i) = 0
!          p_fld1(6+6*(ivar-1),i) = 0
!       end do
!    end do
    
    
    
    

    do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)

       do ivar = 1, nvar2d
          deltaGij(ivar) = p_fld2(2+2*(ivar-1),iedge)
       end do

       ! First calculate the limiting factor alphaij
       ! For this we have three options

       SELECT CASE (syncromethod)
       CASE (1)
          ! 1.: Take only the limiting factor of the variable h=Q(1)
          if (deltaGij(1)>0.0_DP) then
             alphaij = min(p_fld1(5,i),p_fld1(6,j))
          else
             alphaij = min(p_fld1(5,j),p_fld1(6,i))
          end if
          !Limit the antidiffusive fluxes Gij = Gij * alphaij
          deltaGij = deltaGij * alphaij
       CASE (2)
          ! 2.: Take the smallest of all the limiting factors
          alphaij = 1.0_DP
          do ivar = 1, nvar2d
             if (deltaGij(ivar)>0.0_DP) then
                alphaij = min(alphaij,p_fld1(5+6*(ivar-1),i),p_fld1(6+6*(ivar-1),j))
             else
                alphaij = min(alphaij,p_fld1(5+6*(ivar-1),j),p_fld1(6+6*(ivar-1),i))
             end if
          end do
          !Limit the antidiffusive fluxes Gij = Gij * alphaij
          deltaGij = deltaGij * alphaij
       CASE (3)
          ! Last Method: Limit every component on its own...  
          do ivar = 1, nvar2d
             if (deltaGij(ivar)>0.0_DP) then
                alphaij = min(p_fld1(5+6*(ivar-1),i),p_fld1(6+6*(ivar-1),j))
             else
                alphaij = min(p_fld1(5+6*(ivar-1),j),p_fld1(6+6*(ivar-1),i))
             end if
             deltaGij(ivar) = deltaGij(ivar) * alphaij
          end do
       END SELECT

       ! Get matrix to trafo back
       invTij = Eye

       ! Finally we can transform back to deltaFij
       deltaFij = matmul(invTij,deltaGij)

       deltaFij = deltaGij

       ! And save 
       do ivar = 1, nvar2d
          p_fld2(2+2*(ivar-1),iedge) = deltaFij(ivar)
       end do

    end do


    ! In the last step add the limited antidiffusion to the low order predictor
    do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)

       do ivar = 1, nvar2d
          deltaFij(ivar) = p_fld2(2+2*(ivar-1),iedge)

          rarraySol(ivar)%Da(i) = rarraySol(ivar)%Da(i) + deltaFij(ivar)   /p_MLdata(i)*dt
          rarraySol(ivar)%Da(j) = rarraySol(ivar)%Da(j) - deltaFij(ivar)   /p_MLdata(j)*dt
       end do

    end do




  end subroutine linFctShallowWaterAddLimitedAntidiffusion_syncronized
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    ! linearised FCT Method for Shallow Water
  ! adding limited antidiffusion to the low order predictor
  ! using conservative variables variables and a syncronized limiting factor
  subroutine new_syncronized(&
       rarraySol, rarraySolDot, rarrayRhs,&
       rdefBlock, rstempBlock, rsolBlock, rSolDotBlock, &
       rmatrixML, p_CXdata, p_CYdata, p_BXdata, p_BYdata, p_MLdata, p_MCdata, &
       h_fld1, p_fld1, p_fld2, &
       p_Kdiagonal, p_Kedge, NEQ, nedge, &
       gravconst, dt, Method, prelimiting, syncromethod, &
       rtriangulation)

    ! parameter values
    type(t_array), dimension(nvar2d), intent(INOUT) :: rarraySol, rarraySolDot
    type(t_array), dimension(nvar2d), intent(IN)    :: rarrayRhs
    type(t_vectorBlock), intent(INOUT)              :: rdefBlock, rstempBlock
    type(t_vectorBlock), intent(IN)                 :: rsolBlock
    type(t_vectorBlock), intent(INOUT)              :: rSolDotBlock
    type(t_matrixScalar), intent(IN)                :: rmatrixML
    real(DP), dimension(:), pointer, intent(IN)     :: p_CXdata, p_CYdata, p_BXdata, p_BYdata, p_MLdata, p_MCdata
    integer                                         :: h_fld1
    real(DP), dimension(:,:), pointer               :: p_fld1, p_fld2
    integer, dimension(:), pointer, intent(IN)      :: p_Kdiagonal
    integer, dimension(:,:), pointer                :: p_Kedge
    integer, intent(IN)                             :: NEQ, nedge
    real(DP), intent(IN)                            :: gravconst, dt
    integer, intent(IN)                             :: Method, prelimiting, syncromethod
    type(t_triangulation), intent(IN)               :: rtriangulation


    ! variables
    integer                             :: i, j, l, ii, jj, ij, ji, d, iedge, ivar
    real(DP), dimension(nvar2d)         :: deltaQij, Qi, Qj, Qroeij
    real(DP), dimension(nvar2d,nvar2d)  :: JacobixRoeij, JacobiyRoeij
    real(DP)                            :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
    real(DP), dimension(nvar2d,nvar2d)  :: Aij, Bij, Dij, Tij, invTij
    real(DP)                            :: scalarDissipation, scalefactor, lambda
    real(DP), dimension(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj, deltaBi, deltaBj
    real(DP)                            :: cRoe, uRoe, vRoe
    real(DP), dimension(nvar2d)         :: Udoti, Udotj
    real(DP), dimension(nvar2d)         :: deltaWij, deltaGij, deltaFij
    real(DP), dimension(nvar2d,nvar2d)  :: Rij, invRij
    real(DP)                            :: alphaij
    integer                             :: ibct, nbct, ivbd, ivbdFirst, ivbdLast
    integer, dimension(:), pointer      :: p_IboundaryCpIdx
    integer, dimension(:), pointer      :: p_IverticesAtBoundary
    real(dp), dimension(:,:), allocatable :: lPp, lQp, lRp
    real(dp), dimension(:,:), allocatable :: lPm, lQm, lRm
    real(dp), dimension(:,:), allocatable :: deltaF

    !allocate P/Q/R_i^+/-
    allocate(lPp(3,NEQ),lQp(3,NEQ),lRp(3,NEQ))
    allocate(lPm(3,NEQ),lQm(3,NEQ),lRm(3,NEQ))
    
    ! allocate array for antidiff fluxes
    allocate(deltaF(3,nedge))
    
    lPp(:,:) = 0
    lPm(:,:) = 0
    lQp(:,:) = 0
    lQm(:,:) = 0
    lRp(:,:) = 1
    lRm(:,:) = 1
    
        ! Clear SolDot
    call lsysbl_clearVector (rSolDotBlock)

    ! First compute SolDot
    do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)

       ! get solution values at node i and j
       do ivar = 1, nvar2d
          Qi(ivar) = rarraySol(ivar)%Da(i)
          Qj(ivar) = rarraySol(ivar)%Da(j)
       end do

       ! compute deltaQij = Qi - Qj
       deltaQij = Qi - Qj

       ! Compute the Roe-meanvalue for this edge
       Qroeij = calculateQroe(Qi, Qj)


       ! Calculate the galerkin fluxes
       ! deltaKi = c_{ij}*(F(Q_i)-F(Q_j))
       ! deltaKj = c_{ji}*(F(Q_j)-F(Q_i))
       deltaKi = p_CXdata(ij)*(buildFlux(Qi,1,gravconst)-buildFlux(Qj,1,gravconst))+ &
                 p_CYdata(ij)*(buildFlux(Qi,2,gravconst)-buildFlux(Qj,2,gravconst))
       deltaKj = p_CXdata(ji)*(buildFlux(Qj,1,gravconst)-buildFlux(Qi,1,gravconst))+ &
                 p_CYdata(ji)*(buildFlux(Qj,2,gravconst)-buildFlux(Qi,2,gravconst))

       ! compute Dij
        Dij = buildDissipation(Qi,Qj,1,0,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))

       ! deltaD
       deltaDi = matmul(Dij,-deltaQij)
       deltaDj = -deltaDi

!        ! Now we take care of the bottom profile
!        select case (deltabtype)
!         case(0)
!           DeltaBi = 0.0_dp
!           DeltaBj = 0.0_dp
!         case(1)
!           DeltaBi = (/ 0.0_DP, p_BXdata(ij) * Qi(1), p_BYdata(ij) * Qi(1) /)
!           DeltaBj = (/ 0.0_DP, p_BXdata(ji) * Qj(1), p_BYdata(ji) * Qj(1) /)
!         case(2)
!           DeltaBi = (/ 0.0_DP, p_BXdata(ij) * 0.5_dp*(Qi(1)+Qj(1)), p_BYdata(ij) * 0.5_dp*(Qi(1)+Qj(1)) /)
!           DeltaBj = (/ 0.0_DP, p_BXdata(ji) * 0.5_dp*(Qi(1)+Qj(1)), p_BYdata(ji) * 0.5_dp*(Qi(1)+Qj(1)) /)
!         end select

       ! add deltaK and deltaD to SolDot
       ! SolDot = SolDot + 1/ML*(K*u + D*u + Sourceterm(B))
       do l = 1, nvar2d
          rarraySolDot(l)%Da(i) = rarraySolDot(l)%Da(i) &
               + 1.0_DP/p_MLdata(i)*&
               (deltaKi(l) + deltaDi(l))! + DeltaBi(l))
          rarraySolDot(l)%Da(j) = rarraySolDot(l)%Da(j) &
               + 1.0_DP/p_MLdata(j)*&
               (deltaKj(l) + deltaDj(l))! + DeltaBj(l))
          ! Save Qroe
          p_fld2(1+2*(l-1),iedge) = QRoeij(l)
       end do

    end do
    ! SolDot fertig!
    
    
    
    do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)
       
       ! get solution values at node i and j
       Qi = (/rarraySol(1)%Da(i),rarraySol(2)%Da(i),rarraySol(3)%Da(i)/)
       Qj = (/rarraySol(1)%Da(j),rarraySol(2)%Da(j),rarraySol(3)%Da(j)/)
       
       ! Compute the Roe-meanvalue for this edge
       Qroeij = calculateQroe(Qi, Qj)
       !Qroeij = 0.5_dp*(Qi+Qj)
       
       deltaWij = Qj - Qi
       
       
       ! Build Dij using scalar dissipation
       ! Now choose the artificial diffusion method
       select case (Method)
        case (0,3)
          ! Here we do not add artificial diffusion - so we have the high order method
          Dij = 0
        case (1,4,6)
          ! Here we use scalar dissipation
          Dij = buildDissipation(Qi,Qj,1,0,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
        case (2,5,7)
          ! Here we use tensorial dissipation
          Dij = abs(p_CXdata(ij))* buildDissipation(Qi,Qj,2,1,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))+&
                abs(p_CYdata(ij))* buildDissipation(Qi,Qj,2,2,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
       end select
       
       ! get Udot at node i and j
       do ivar = 1, 3
          Udoti(ivar) = rarraySolDot(ivar)%Da(i)
          Udotj(ivar) = rarraySolDot(ivar)%Da(j)
       end do
       
       
       
!        JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
!        JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP
!        JacobixRoeij = buildJacobi(Qroeij,1,gravconst)
!        JacobiyRoeij = buildJacobi(Qroeij,2,gravconst)
!        ! Now we can compute Aij and Bij
!        Bij = JcoeffxB*JacobixRoeij + JcoeffyB*JacobiyRoeij
       
       
       ! compute antidiffusive flux
       deltaFij = p_MCdata(ij)*(Udoti-Udotj)+ matmul(Dij,(Qi-Qj))
       !deltaFij = matmul(Dij+Bij,(Qi-Qj))
!        deltaFij = matmul(Dij,(Qi-Qj))
       
       ! prelimit antidiff flux
!         do ivar = 1, 3
!            if (deltaWij(ivar)*deltaFij(ivar).ge.0) deltaFij(ivar) = 0.0_dp
!         end do
       
       ! and save
       deltaF(:,iedge) = deltaFij
       
       do ivar = 1, 3
       ! Augment the sums of positive/negative fluxes of indicator variable
       lPp(ivar,i) = lPp(ivar,i) + max(0.0_dp,deltaFij(ivar))
       lPm(ivar,i) = lPm(ivar,i) + min(0.0_dp,deltaFij(ivar))
       
       lPp(ivar,j) = lPp(ivar,j) + max(0.0_dp,-deltaFij(ivar))
       lPm(ivar,j) = lPm(ivar,j) + min(0.0_dp,-deltaFij(ivar))
       
       ! Update maximum/minimum admissible increments
       lQp(ivar,i) = max(lQp(ivar,i), deltaWij(ivar))
       lQm(ivar,i) = min(lQm(ivar,i), deltaWij(ivar))
       
       lQp(ivar,j) = max(lQp(ivar,j), -deltaWij(ivar))
       lQm(ivar,j) = min(lQm(ivar,j), -deltaWij(ivar))
       
       end do
     end do !iedge
     
    
     
     ! Compute the Ri+/-
     do i = 1, NEQ
     do ivar = 1, 3
       lRp(ivar,i) = min(1.0_DP,p_MLdata(i)/dt*lQp(ivar,i)/(lPp(ivar,i)+SYS_EPSREAL_DP))! Ri+
       lRm(ivar,i) = min(1.0_DP,p_MLdata(i)/dt*lQm(ivar,i)/(lPm(ivar,i)-SYS_EPSREAL_DP))! Ri-

     end do
     end do ! i
     

     do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)
       
       deltaFij = deltaF(:,iedge)
       
       select case (syncromethod)
       case(1)
       
       if (deltaFij(1) .ge. 0.0_dp) then
         alphaij = min(lRp(1,i),lRm(1,j))
       else
         alphaij = min(lRp(1,j),lRm(1,i))
       end if
       
       case(2)
       alphaij = 1.0_dp
       do ivar = 1, 3
       if (deltaFij(ivar) .ge. 0.0_dp) then
         alphaij = min(alphaij,lRp(ivar,i),lRm(ivar,j))
       else
         alphaij = min(alphaij,lRp(ivar,j),lRm(ivar,i))
       end if
       end do
       
       end select
       
       deltaFij = alphaij * deltaFij
       
       do ivar = 1, 3
         rarraySol(ivar)%Da(i) = rarraySol(ivar)%Da(i) + deltaFij(ivar)/p_MLdata(i)*dt
         rarraySol(ivar)%Da(j) = rarraySol(ivar)%Da(j) - deltaFij(ivar)/p_MLdata(j)*dt
       end do
    
     end do ! iedge
    
    ! deallocate P/Q/R_i^+/-
    deallocate(lPp,lQp,lRp)
    deallocate(lPm,lQm,lRm)
    
    ! deallocate array for antidiff fluxes
    deallocate(deltaF)

  end subroutine
  
  
  
  
  
  
  
  
  
  
  ! Add soure term explicitlely
  subroutine AddExplicitSourceTerm(rarraySol,dt,neq,h_DvertexCoords,gravconst)

    type(t_array), dimension(nvar2d), intent(inout)    :: rarraySol
    
    real(dp), intent(in) :: dt, gravconst
    
    integer(I32), intent(in) :: h_DvertexCoords
    
    integer, intent(in) :: neq
  
  ! Local variables
  
  ! Pointer to vertex coordinates
    real(DP), dimension(:,:), pointer   :: p_dVertexCoords
     
    real(DP) :: bx, by
    
    integer :: ieq
  
  
  
    ! get Pointer to the vertex coordinates
    call storage_getbase_double2d(h_DvertexCoords, p_dVertexCoords)
    
    bx=1
    by=0

    ! walk over all edges and add the source
    do ieq = 1, neq
    
      bx=0
      by=0
    
!       if ((p_dVertexCoords(1,ieq)>5).and.(p_dVertexCoords(1,ieq)<7)) then
!         bx=1
!       end if
    
      rarraySol(2)%Da(ieq)=rarraySol(2)%Da(ieq)-dt*gravconst*rarraySol(1)%Da(ieq)*bx
      rarraySol(3)%Da(ieq)=rarraySol(3)%Da(ieq)-dt*gravconst*rarraySol(1)%Da(ieq)*by
    end do

  end subroutine
  
  
  
  ! Add bottom profile before writing the solution
  subroutine AddBottomBeforeWrite(rarraySol,neq,h_bottom)

    type(t_array), dimension(nvar2d), intent(inout)    :: rarraySol
    
    
    integer(I32), intent(in) :: h_bottom
    
    integer, intent(in) :: neq
  
  ! Local variables
  
  ! Pointer to vertex coordinates
    real(DP), dimension(:), pointer   :: p_bottom
    
    integer :: ieq
  
  
    ! get Pointer to the bottom profile
    call storage_getbase_double(h_bottom, p_bottom)
    
    ! walk over all edges and add the source
    do ieq = 1, neq
    
      rarraySol(1)%Da(ieq)=rarraySol(1)%Da(ieq)+ p_bottom(ieq) !bottom heigth
      
    end do


  end subroutine
  
  
  ! Substract bottom profile after writing the solution
  subroutine SubstractBottomAfterWrite(rarraySol,neq,h_bottom)

    type(t_array), dimension(nvar2d), intent(inout)    :: rarraySol
    
    
    integer(I32), intent(in) :: h_bottom
    
    integer, intent(in) :: neq
  
  ! Local variables
  
  ! Pointer to vertex coordinates
    real(DP), dimension(:), pointer   :: p_bottom
    
    integer :: ieq
  
  
    ! get Pointer to the bottom profile
    call storage_getbase_double(h_bottom, p_bottom)
    
    ! walk over all edges and add the source
    do ieq = 1, neq
    
      rarraySol(1)%Da(ieq)=rarraySol(1)%Da(ieq)- p_bottom(ieq) !bottom heigth
      
    end do


  end subroutine
  
  
  
  
  
  
  
  ! This routine can be used to add the bottom term to the rhs or the defect vector
  subroutine addBottomTermToVec(rSolBlock,rVecBlock,rsourceBlock,rbottomBlock, &
                                Rcollection,coefficient,gravconst)
  
  ! Solution
  type(t_vectorBlock), intent(in), target :: rSolBlock
  
  ! Vector that is to be changed, i.e. rhs or defect
  type(t_vectorBlock), intent(inout) :: rVecBlock
  
  type(t_vectorBlock), intent(inout) :: rsourceBlock
  
  type(t_vectorBlock), intent(in), target :: rbottomBlock
  
  ! The collection to give data to the callback routine
  type(t_collection), intent(inout) :: Rcollection
  
  ! The coefficient before the source term, i.e. theta*dt
  real(dp), intent(in) :: coefficient
  
  ! The gravitation constant
  real(dp), intent(in) :: gravconst
  
  
  
  !! Local variables

  ! The linear form describing the source term
    type(t_linearForm) :: rform
  
  ! New try :-) Integration of \int \phi_i S(Q(x)) dx
  ! While S(Q(X)) = (0, -g h b_x, -g h b_y)
  ! \int \phi_i S(Q(x)) dx = \int -g (0, \phi_i_x b h, \phi_i_y b h) dx
  
  
  ! Build second component of the source term
  
  
  ! Set up the corresponding linear form
  rform%itermCount      = 1
  rform%Idescriptors(1) = DER_FUNC
  
  ! Hang into the collection the pointer to the solution vector
  Rcollection%p_rvectorQuickAccess1 => rsolBlock
  
  ! Hang into the collection the pointer to the bottom vector
  Rcollection%p_rvectorQuickAccess2 => rbottomBlock
  
  ! Derivate the bottom profile in direction ...
  Rcollection%IquickAccess(1) = DER_DERIV_X
  
  
  ! Build the discretized target functional.
  call linf_buildVectorScalar2(rform, .true., rsourceBlock&
       %RvectorBlock(2), shlw_SourceTermCB, rcollection)
  
  
  ! Build third component of the source term
  
  
! Set up the corresponding linear form
  rform%itermCount      = 1
  rform%Idescriptors(1) = DER_FUNC
  
  ! Hang into the collection the pointer to the solution vector
  Rcollection%p_rvectorQuickAccess1 => rsolBlock
  
  ! Hang into the collection the pointer to the bottom vector
  Rcollection%p_rvectorQuickAccess2 => rbottomBlock
  
    ! Derivate the bottom profile in direction ...
  Rcollection%IquickAccess(1) = DER_DERIV_Y
  
  
  ! Build the discretized target functional.
  call linf_buildVectorScalar2(rform, .true., rsourceBlock&
       %RvectorBlock(3), shlw_SourceTermCB, rcollection)



  ! Add the calculated sourceterm to rhs/defect
  call lsysbl_vectorLinearComb (rsourceBlock,rVecBlock,-gravconst*coefficient,1.0_DP)


  end subroutine
  
  
  
  
  
  subroutine lfctsync(rarraySol, rarraySolDot, rSolDotBlock, &
                      p_Kedge, p_Kdiagonal, NEQ, nedge,&
                      p_CXdata, p_CYdata, p_MLdata, p_MCdata, &
                      gravconst, dt, syncromethod)
  
  ! parameter values
  type(t_array), dimension(nvar2d), intent(INOUT) :: rarraySol, rarraySolDot
  type(t_vectorBlock), intent(INOUT)              :: rSolDotBlock
  integer, dimension(:,:), pointer                :: p_Kedge
  integer, dimension(:), pointer, intent(IN)      :: p_Kdiagonal
  integer, intent(IN)                             :: NEQ, nedge
  real(DP), dimension(:), pointer, intent(IN)     :: p_CXdata, p_CYdata, p_MLdata, p_MCdata
  real(DP), intent(IN)                            :: gravconst, dt
  integer, intent(IN)                             :: syncromethod

  ! local variables
  real(DP), dimension(:,:), allocatable           :: Pp, Pm, Qp, Qm, Rp, Rm
  real(DP), dimension(:,:), allocatable           :: Fijs
  real(DP), dimension(nvar2d)                     :: Qi, Qj, Udoti, Udotj
  integer                                         :: i, j, ij, ji, ii, jj
  integer                                         :: iedge, ivar, inode
  real(DP), dimension(nvar2d, nvar2d)             :: Dij
  real(dp), dimension(nvar2d)                     :: Fij
  real(dp), dimension(nvar2d)                     :: deltaQji, deltaQij
  real(dp)                                        :: alphaij
  real(dp), dimension(nvar2d)                     :: deltaKi, deltaKj, deltaDi, deltaDj
  
  
  ! function code
  
  ! Allocate memory for nodewise correction factors
  allocate(Pp(nvar2d,neq),Pm(nvar2d,neq),Qp(nvar2d,neq),Qm(nvar2d,neq),Rp(nvar2d,neq),Rm(nvar2d,neq))
  allocate(Fijs(nvar2d,nedge))
  
  ! Initialise the arrays
  Pp = 0.0_dp
  Pm = 0.0_dp
  Qp = 0.0_dp
  Qm = 0.0_dp
  Rp = 1.0_dp
  Rm = 1.0_dp
  
  
  
  
  
  
  
  
  
  ! Clear SolDot
    call lsysbl_clearVector (rSolDotBlock)

    ! First compute SolDot
    do iedge = 1, nedge
       i  = p_Kedge(1, iedge)
       j  = p_Kedge(2, iedge)
       ij = p_Kedge(3, iedge)
       ji = p_Kedge(4, iedge)
       ii = p_Kdiagonal(i)
       jj = p_Kdiagonal(j)

       ! get solution values at node i and j
       do ivar = 1, nvar2d
          Qi(ivar) = rarraySol(ivar)%Da(i)
          Qj(ivar) = rarraySol(ivar)%Da(j)
       end do

       ! compute deltaQij = Qi - Qj
       deltaQij = Qi - Qj

       ! Calculate the galerkin fluxes
       ! deltaKi = c_{ij}*(F(Q_i)-F(Q_j))
       ! deltaKj = c_{ji}*(F(Q_j)-F(Q_i))
       deltaKi = p_CXdata(ij)*(buildFlux(Qi,1,gravconst)-buildFlux(Qj,1,gravconst))+ &
                 p_CYdata(ij)*(buildFlux(Qi,2,gravconst)-buildFlux(Qj,2,gravconst))
       deltaKj = p_CXdata(ji)*(buildFlux(Qj,1,gravconst)-buildFlux(Qi,1,gravconst))+ &
                 p_CYdata(ji)*(buildFlux(Qj,2,gravconst)-buildFlux(Qi,2,gravconst))

       ! compute Dij
       Dij = buildDissipation(Qi,Qj,1,0,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))

       ! deltaD
       deltaDi = matmul(Dij,-deltaQij)
       deltaDj = -deltaDi


       ! add deltaK and deltaD to SolDot
       ! SolDot = SolDot + 1/ML*(K*u + D*u)
       do ivar = 1, nvar2d
          rarraySolDot(ivar)%Da(i) = rarraySolDot(ivar)%Da(i) &
               + 1.0_DP/p_MLdata(i)*&
               (deltaKi(ivar) + deltaDi(ivar))
          rarraySolDot(ivar)%Da(j) = rarraySolDot(ivar)%Da(j) &
               + 1.0_DP/p_MLdata(j)*&
               (deltaKj(ivar) + deltaDj(ivar))
       end do

    end do
    ! SolDot fertig!
  
  
  
  ! Walk over all edges
  do iedge = 1, nedge
    i  = p_Kedge(1, iedge)
    j  = p_Kedge(2, iedge)
    ij = p_Kedge(3, iedge)
    ji = p_Kedge(4, iedge)
    ii = p_Kdiagonal(i)
    jj = p_Kdiagonal(j)
    
    ! Get values of the solution in node i and j
    do ivar = 1, nvar2d
      Qi(ivar) = rarraySol(ivar)%da(i)
      Qj(ivar) = rarraySol(ivar)%da(j)
    end do
    
    ! Get the diffusion matrix
    Dij = buildDissipation(Qi,Qj,1,0,gravconst,p_CXdata(ij),p_CYdata(ij),p_CXdata(ji),p_CYdata(ji))
    
    ! get Udot at node i and j
    do ivar = 1, 3
      Udoti(ivar) = rarraySolDot(ivar)%Da(i)
      Udotj(ivar) = rarraySolDot(ivar)%Da(j)
    end do

    ! Calculate the raw antidiffusion
    Fij = matmul(Dij,Qi-Qj)
    Fij = p_MCdata(ij)*(Udoti-Udotj)+matmul(Dij,Qi-Qj)
    
    ! Save the raw antidiffusive flux
    Fijs(:,iedge) = Fij

    ! Calculate the difference of the solution values in this node
    deltaQji = Qj - Qi

    do ivar = 1, nvar2d
      ! Augment sum of positiv/negativ fluxes
      Pp(ivar,i) = Pp(ivar,i) + max( Fij(ivar),0.0_dp)
      Pp(ivar,j) = Pp(ivar,j) + max(-Fij(ivar),0.0_dp)
      Pm(ivar,i) = Pm(ivar,i) + min( Fij(ivar),0.0_dp)
      Pm(ivar,j) = Pm(ivar,j) + min(-Fij(ivar),0.0_dp)

      ! Update the maximum/minimum admissible increment
      Qp(ivar,i) = max(Qp(ivar,i), deltaQji(ivar))
      Qp(ivar,j) = max(Qp(ivar,j),-deltaQji(ivar))
      Qm(ivar,i) = min(Qm(ivar,i), deltaQji(ivar))
      Qm(ivar,j) = min(Qm(ivar,j),-deltaQji(ivar))
    end do

  end do ! iedge

  ! Loop over all nodes
  do inode = 1, NEQ
    do ivar = 1, nvar2d
      ! Calculate nodal correction factors
      if (Pp(ivar,inode)>0.0_dp) then
        Rp(ivar,inode) = min(1.0_dp, p_MLdata(inode)/dt*Qp(ivar,inode)/(Pp(ivar,inode)+sys_epsreal))
      else
        Rp(ivar,inode) = 0.0_dp
      end if
      if (Pm(ivar,inode)<0.0_dp) then
        Rm(ivar,inode) = min(1.0_dp, p_MLdata(inode)/dt*Qm(ivar,inode)/(Pm(ivar,inode)-sys_epsreal))
      else
        Rm(ivar,inode) = 0.0_dp
      end if
    end do
  end do ! inode
  
    ! Walk over all edges
  do iedge = 1, nedge
    i  = p_Kedge(1, iedge)
    j  = p_Kedge(2, iedge)
    ij = p_Kedge(3, iedge)
    ji = p_Kedge(4, iedge)
    ii = p_Kdiagonal(i)
    jj = p_Kdiagonal(j)
    
    ! Get values of the solution in node i and j
    do ivar = 1, nvar2d
      Qi(ivar) = rarraySol(ivar)%da(i)
      Qj(ivar) = rarraySol(ivar)%da(j)
    end do
    
    ! Get the raw antidiffusive flux
    Fij = Fijs(:,iedge)
    
    ! Calculate limiting factor by using water heigth as indicator variable
    if (Fij(1)>0.0_dp) then
      alphaij = min(Rp(1,i),Rm(1,j))
    elseif (Fij(1)<0.0_dp) then
      alphaij = min(Rp(1,j),Rm(1,i))
    else
      alphaij = 0.0_dp
    end if
    

    
    ! We can even use the minimal limiting factor
!     do ivar=2,nvar2d
!       if (Fij(ivar)>0.0_dp) then
!         alphaij = min(Rp(ivar,i),Rm(ivar,j),alphaij)
!       elseif (Fij(ivar)<0.0_dp) then
!         alphaij = min(Rp(ivar,j),Rm(ivar,i),alphaij)
!       else
!         alphaij = 0.0_dp
!       end if
!     end do

    ! If we are near a dry bed occurence, use the low order method
    if ((( dt/p_MLdata(i)*Fij(1))/Qi(1)<-0.1_dp).or.(Qi(1)<clipwater)) alphaij = 0.0_dp
    if (((-dt/p_MLdata(j)*Fij(1))/Qj(1)<-0.1_dp).or.(Qj(1)<clipwater)) alphaij = 0.0_dp
    
    ! Limit the antidiffusive flux
    Fij = alphaij*Fij
    
    ! Update the Solution with the limited antidiffusive flux
    do ivar = 1, nvar2d
      rarraySol(ivar)%da(i) = rarraySol(ivar)%da(i) + dt/p_MLdata(i)*Fij(ivar)
      rarraySol(ivar)%da(j) = rarraySol(ivar)%da(j) - dt/p_MLdata(j)*Fij(ivar)
    end do

  end do ! iedge
  
  
  ! Deallocate the memory used for the nodewise correction factors
  deallocate(Pp,Pm,Qp,Qm,Rp,Rm)
  deallocate(Fijs)
  
  end subroutine

end module shallowwater2d_routines
