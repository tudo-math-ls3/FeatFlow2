module dg2d_problem


  use fsystem
  use storage

  implicit none
  
  real(dp), parameter :: g = 1.0_dp

    

contains

        
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
    integer, parameter :: mvk = 0


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
        
      case (2) ! Roe-meanvalues

        ! Set the height variables
        hl=Ql(1)
        hr=Qr(1)

        denom = sqrt(hl)+sqrt(hr)
        whl = 1.0_DP/sqrt(hl)
        whr = 1.0_DP/sqrt(hr)

        Qroe(1) = 0.5_dp*(hl+hr)
        Qroe(2) = Qroe(1)*(whl*Ql(2)+whr*Qr(2))/denom
        Qroe(3) = Qroe(1)*(whl*Ql(3)+whr*Qr(3))/denom

      end select

  end function calculateQroe



  ! This routine builds the jacobi matrix DF of the flux in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildJacobi(Q,d) result(J)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)	:: J

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: q

    integer, intent(IN)                     :: d

    ! primitive variables
    real(DP)                                :: h, u, v, c

    ! Calculate primitive variables
    h=Q(1)
!    if (h<clipwater) then
!      h=0.0_dp
!      u=0.0_dp
!      v=0.0_dp
!    else
      u=Q(2)/Q(1)
      v=Q(3)/Q(1)
!    end if

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



  ! This routine builds the trafo matrix T for the jacobian DF of the flux in direction d (right eigenvalues)
  ! This means :   invTrafo * DF * Trafo    = D (diagonal matrix)
  ! or         :   Trafo    * D  * invTrafo = DF
  ! So to get the characteristic variables, the conservative variables have to be multiplied by
  ! invTrafo
  ! d=1: x-direction, d=2: y-direction
  function buildTrafo(Q,d) result(Rij)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)	:: Rij

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

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



  ! This routine builds the inverse of the trafo matrix T for the jacobian DF of the flux in direction d (left eigenvalues)
  ! This means :   invTrafo * DF * Trafo    = D (diagonal matrix)
  ! or         :   Trafo    * D  * invTrafo = DF
  ! So to get the characteristic variables, the conservative variables have to be multiplied by
  ! invTrafo
  ! d=1: x-direction, d=2: y-direction
  function buildInvTrafo(Q,d) result(invRij)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)	:: invRij

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

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
  function buildLambda(Q,d) result(Lambda)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)	:: Lambda

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

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
  function buildaLambda(Q,d) result(aLambda)

    ! The jacobi matrix in direction d
    real(DP), dimension(3,3)	:: aLambda

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

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
       aLambda(1,1) = min(0.1_dp,abs(v-c))
       aLambda(2,1) = 0.0_DP
       aLambda(3,1) = 0.0_DP
       aLambda(1,2) = 0.0_DP
       aLambda(2,2) = min(0.1_dp,abs(v))
       aLambda(3,2) = 0.0_DP
       aLambda(1,3) = 0.0_DP
       aLambda(2,3) = 0.0_DP
       aLambda(3,3) = min(0.1_dp,abs(v+c))
    end if
    
    
!    ! With entropy fix
!    if (d==1) then
!       ! build aLambda in x direction
!       aLambda(1,1) = min(0.1_dp,abs(u-c))
!       aLambda(2,1) = 0.0_DP
!       aLambda(3,1) = 0.0_DP
!       aLambda(1,2) = 0.0_DP
!       aLambda(2,2) = min(0.1_dp,abs(u))
!       aLambda(3,2) = 0.0_DP
!       aLambda(1,3) = 0.0_DP
!       aLambda(2,3) = 0.0_DP
!       aLambda(3,3) = min(0.1_dp,abs(u+c))
!    else
!       ! build aLambda in y direction
!       aLambda(1,1) = min(0.1_dp,abs(v-c))
!       aLambda(2,1) = 0.0_DP
!       aLambda(3,1) = 0.0_DP
!       aLambda(1,2) = 0.0_DP
!       aLambda(2,2) = min(0.1_dp,abs(v))
!       aLambda(3,2) = 0.0_DP
!       aLambda(1,3) = 0.0_DP
!       aLambda(2,3) = 0.0_DP
!       aLambda(3,3) = min(0.1_dp,abs(v+c))
!    end if

  end function buildaLambda






  ! This routine returns the eigenvalues of the jacobi matrix in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildEigenvalues(Q,d) result(Eigenvalues)

    ! The jacobi matrix in direction d
    real(DP), dimension(3)	:: Eigenvalues

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    ! the direction: d=1: x-direction, d=2: y-direction
    integer, intent(IN)                     :: d

    ! speed of gravitational waves
    real(DP)                                :: c

    ! temporary variable
    real(DP)                                :: coeff, h, u, v

    ! Calculate primitive variables
    h=Q(1)
!    if (h<clipwater) then
!      h = 0.0_dp
!      u = 0.0_dp
!      v = 0.0_dp
!    else
      u = Q(2)/Q(1)
      v = Q(3)/Q(1)
!    end if

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
  function buildFlux(Q,d) result(Flux)

    ! The flux vector in direction d at Q
    real(DP), dimension(3)	:: Flux

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

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
!     if (Q(1)<clipwater) then
!        !dry bed case
!        Flux=0
!     else
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
!     end if ! dry or wet bed
  end function buildFlux
  
  
  
  
  
  
  ! This routine builds the right eigenvectors for the mixed jacobian
  function buildMixedR(Q,a,b) result(R)

    ! Right eigenvectors
    real(DP), dimension(3,3)	:: R

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: c1, c2, c3
    
    ! Temp array
    real(dp), dimension(3) :: T
    
    if (abs(a)<10*SYS_EPSREAL) then
    
      R = buildTrafo(Q,2)
      
      T = R(:,2)
      R(:,2) = R(:,3)
      R(:,3) = R(:,1)
      R(:,1) = T
    
    else
    
      c1 = a*Q(2)+b*Q(3)
      c2 = Q(1)*(a*a+b*b)
      c3 = sqrt(c2*Q(1)*Q(1)*g)
    
    
      ! Build matrix of right eigenvectors
      R(1,1) = 0.0_DP
      R(2,1) = -b/a
      R(3,1) = 1.0_dp
      R(1,2) = 1.0_DP
      R(2,2) = ((c1+c3)*a+b*b*Q(2)-a*b*Q(3))/c2
      R(3,2) = ((c1+c3)*b+a*a*Q(3)-a*b*Q(2))/c2
      R(1,3) = 1.0_DP
      R(2,3) = ((c1-c3)*a+b*b*Q(2)-a*b*Q(3))/c2
      R(3,3) = ((c1-c3)*b+a*a*Q(3)-a*b*Q(2))/c2
    end if
    
  end function 
  
  
  
  ! This routine builds the diagonalmatrix of the absolut value of the eigenvalues
  function buildMixedaLambda(Q,a,b) result(aLambda)

    ! Left eigenvectors
    real(DP), dimension(3,3)	:: aLambda

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: c1, c2, c3, lambda1, lambda2, lambda3
    
    c1 = a*Q(2)+b*Q(3)
    c2 = Q(1)*(a*a+b*b)
    c3 = sqrt(c2*Q(1)*Q(1)*g)
    
    ! Calculate eigenvalues
    lambda1 = c1/Q(1)
    lambda2 = (c1+c3)/Q(1)
    lambda3 = (c1-c3)/Q(1)
    
    aLambda = 0.0_dp
    
    ! Build matrix of left eigenvectors
    aLambda(1,1) = abs(lambda1)
    aLambda(2,2) = abs(lambda2)
    aLambda(3,3) = abs(lambda3)
    
  end function 


  ! This routine builds the left eigenvectors for the mixed jacobian
  function buildMixedL(Q,a,b) result(L)

    ! Left eigenvectors
    real(DP), dimension(3,3)	:: L

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: c1, c2, c3
    
    c1 = a*Q(2)+b*Q(3)
    c2 = Q(1)*(a*a+b*b)
    c3 = sqrt(c2*Q(1)*Q(1)*g)
    
    
       ! Build matrix of left eigenvectors
       L(1,1) = -a*(a*Q(3)-b*Q(2))/c2
       L(2,1) = -0.5_dp*(c1-c3)/c3
       L(3,1) = 0.5_dp*(c1+c3)/c3
       L(1,2) = -a*b/(a*a+b*b)
       L(2,2) = 0.5*a*c2/(c3*(a*a+b*b))
       L(3,2) = -0.5*a*c2/(c3*(a*a+b*b))
       L(1,3) = a*a/(a*a+b*b)
       L(2,3) = 0.5_dp*b*c2/(c3*(a*a+b*b))
       L(3,3) = -0.5_dp*b*c2/(c3*(a*a+b*b))
    
  end function 


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ! This routine builds the right eigenvectors for the mixed jacobian
  function buildMixedR2(Q,a,b) result(R)

    ! Right eigenvectors
    real(DP), dimension(3,3)	:: R

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: c1, c2, c3, h, u, v, c
    
    ! Temp array
    real(dp), dimension(3) :: T
    
    
    h = Q(1)
    u = Q(2)/h
    v = Q(3)/h
    c = sqrt(h*g)
    
        
      ! Build matrix of right eigenvectors
      R(1,1) = 1.0_DP
      R(2,1) = u-c*a
      R(3,1) = v-c*b
      R(1,2) = 0.0_dp
      R(2,2) = b
      R(3,2) = -a
      R(1,3) = 1.0_DP
      R(2,3) = u+c*a
      R(3,3) = v+c*b

    
  end function 
  
  
  
  ! This routine builds the diagonalmatrix of the absolut value of the eigenvalues
  function buildMixedaLambda2(Q,a,b) result(aLambda)

    ! Left eigenvectors
    real(DP), dimension(3,3)	:: aLambda

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: c1, c2, c3, lambda1, lambda2, lambda3, u, v, h, c
    
    
    h = Q(1)
    u = Q(2)/h
    v = Q(3)/h
    c = sqrt(h*g)
    
    ! Calculate eigenvalues
    lambda2 = u*a+v*b
    lambda1 = lambda2-c
    lambda3 = lambda2+c
    
    aLambda = 0.0_dp
    
    ! Build matrix of left eigenvectors
    aLambda(1,1) = abs(lambda1)
    aLambda(2,2) = abs(lambda2)
    aLambda(3,3) = abs(lambda3)
    
  end function 


  ! This routine builds the left eigenvectors for the mixed jacobian
  function buildMixedL2(Q,a,b) result(L)

    ! Left eigenvectors
    real(DP), dimension(3,3)	:: L

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: c1, c2, c3, h, u, v, c
    
    h = Q(1)
    u = Q(2)/h
    v = Q(3)/h
    c = sqrt(h*g)

    c1=0.5_dp/c
    
    
       ! Build matrix of left eigenvectors
       L(1,1) = c1*(u*a+v*b)+0.5_dp
       L(2,1) = v*a-u*b
       L(3,1) = -c1*(u*a+v*b)+0.5_dp
       L(1,2) = -c1*a
       L(2,2) = b
       L(3,2) = c1*a
       L(1,3) = -c1*b
       L(2,3) = -a
       L(3,3) = c1*b
           
  end function 
  
  
  
  
  ! This routine returns the eigenvalues of the jacobi matrix in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildEigenvalues2(Q,a,b) result(Eigenvalues)

    ! The jacobi matrix in direction d
    real(DP), dimension(3)	:: Eigenvalues

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Q

    ! the direction: d=1: x-direction, d=2: y-direction
    real(dp), intent(IN)                     :: a,b

    ! speed of gravitational waves
    real(DP)                                :: c

    ! temporary variable
    real(DP)                                :: coeff, h, u, v, lambda1, lambda2, lambda3

    h = Q(1)
    u = Q(2)/h
    v = Q(3)/h
    c = sqrt(h*g)
    
    ! Calculate eigenvalues
    lambda2 = u*a+v*b
    lambda1 = lambda2-c
    lambda3 = lambda2+c
    
    
    ! build eigenvalues in y direction
    Eigenvalues(1) = lambda1
    Eigenvalues(2) = lambda2
    Eigenvalues(3) = lambda3
    

  end function buildEigenvalues2
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ! This function returns the Roe mean values
  function calculateQroec(Ql, Qr) result(Qroe)

    ! The left and right Q values
    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(3), intent(IN)		:: Ql, Qr

    ! The computed Roe values
    real(DP), dimension(4)					:: Qroe

    ! temp variables
    real(DP)		:: whl, whr, denom, hl, hr


        ! Set the height variables
        hl=Ql(1)
        hr=Qr(1)

        denom = sqrt(hl)+sqrt(hr)
        whl = 1.0_DP/sqrt(hl)
        whr = 1.0_DP/sqrt(hr)

        Qroe(1) = sqrt(hl*hr)
        Qroe(2) = Qroe(1)*(whl*Ql(2)+whr*Qr(2))/denom
        Qroe(3) = Qroe(1)*(whl*Ql(3)+whr*Qr(3))/denom
        
        Qroe(4) = sqrt(0.5_dp*g*(hl+hr))


  end function calculateQroec
  
      ! This routine builds the right eigenvectors for the mixed jacobian
  function buildMixedR2c(Q,a,b) result(R)

    ! Right eigenvectors
    real(DP), dimension(3,3)	:: R

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(4), intent(IN)		:: Q

    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: c1, c2, c3, h, u, v, c
    
    ! Temp array
    real(dp), dimension(3) :: T
    
    
    h = Q(1)
    u = Q(2)/h
    v = Q(3)/h
    c = Q(4)
    
        
      ! Build matrix of right eigenvectors
      R(1,1) = 1.0_DP
      R(2,1) = u-c*a
      R(3,1) = v-c*b
      R(1,2) = 0.0_dp
      R(2,2) = b
      R(3,2) = -a
      R(1,3) = 1.0_DP
      R(2,3) = u+c*a
      R(3,3) = v+c*b

    
  end function 
  
  
  
  ! This routine builds the diagonalmatrix of the absolut value of the eigenvalues
  function buildMixedaLambda2c(Q,a,b) result(aLambda)

    ! Left eigenvectors
    real(DP), dimension(3,3)	:: aLambda

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(4), intent(IN)		:: Q

    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: c1, c2, c3, lambda1, lambda2, lambda3, u, v, h, c
    
    
    h = Q(1)
    u = Q(2)/h
    v = Q(3)/h
    c = Q(4)
    
    ! Calculate eigenvalues
    lambda2 = u*a+v*b
    lambda1 = lambda2-c
    lambda3 = lambda2+c
    
    aLambda = 0.0_dp
    
    ! Build matrix of left eigenvectors
    aLambda(1,1) = abs(lambda1)
    aLambda(2,2) = abs(lambda2)
    aLambda(3,3) = abs(lambda3)
    
  end function 


  ! This routine builds the left eigenvectors for the mixed jacobian
  function buildMixedL2c(Q,a,b) result(L)

    ! Left eigenvectors
    real(DP), dimension(3,3)	:: L

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(4), intent(IN)		:: Q

    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: c1, c2, c3, h, u, v, c
    
    h = Q(1)
    u = Q(2)/h
    v = Q(3)/h
    c = Q(4)

    c1=0.5_dp/c
    
    
       ! Build matrix of left eigenvectors
       L(1,1) = c1*(u*a+v*b)+0.5_dp
       L(2,1) = v*a-u*b
       L(3,1) = -c1*(u*a+v*b)+0.5_dp
       L(1,2) = -c1*a
       L(2,2) = b
       L(3,2) = c1*a
       L(1,3) = -c1*b
       L(2,3) = -a
       L(3,3) = c1*b
           
  end function 
  
  
  
  
  ! This routine returns the eigenvalues of the jacobi matrix in direction d
  ! d=1: x-direction, d=2: y-direction
  function buildEigenvalues2c(Q,a,b) result(Eigenvalues)

    ! The jacobi matrix in direction d
    real(DP), dimension(3)	:: Eigenvalues

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(4), intent(IN)		:: Q

    ! the direction: d=1: x-direction, d=2: y-direction
    real(dp), intent(IN)                     :: a,b

    ! speed of gravitational waves
    real(DP)                                :: c

    ! temporary variable
    real(DP)                                :: coeff, h, u, v, lambda1, lambda2, lambda3

    h = Q(1)
    u = Q(2)/h
    v = Q(3)/h
    c = Q(4)
    
    ! Calculate eigenvalues
    lambda2 = u*a+v*b
    lambda1 = lambda2-c
    lambda3 = lambda2+c
    
    
    ! build eigenvalues in y direction
    Eigenvalues(1) = lambda1
    Eigenvalues(2) = lambda2
    Eigenvalues(3) = lambda3
    

  end function buildEigenvalues2c
    
    
    
    
    
    
    
    
    
    
end module
