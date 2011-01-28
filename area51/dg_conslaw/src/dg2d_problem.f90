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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ! This routine returns the flux vector for the 2d compressible euler
  ! equations of gas dynamics in direction d
  ! d=1: x-direction, d=2: y-direction
  function Euler_buildFlux(Q,d) result(Flux)

    ! The flux vector in direction d at Q
    real(DP), dimension(4)	:: Flux

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(4), intent(IN)		:: Q

    ! the direction: d=1: x-direction, d=2: y-direction
    integer, intent(IN)                     :: d

    ! pressure, stagnation enthalpy
    real(DP)                                :: p, H

    ! temporary variable
    real(DP)                                :: rho, u, v, E
    
    ! Constant Gamma
    real(dp) :: gamma = 1.4_dp

    ! Calculate primitive variables
    rho=Q(1)
    u=Q(2)/rho
    v=Q(3)/rho
    E=Q(4)/rho
    
    ! Compute the pressure
    p = (gamma - 1.0_dp)*rho*(E-0.5_dp*(u*u+v*v))
    
    ! Compute H, the stagnation enthalpy
    H = E + p/rho

    if (d==1) then
      ! build Flux in x direction
      Flux(1) = Q(2)
      Flux(2) = Q(2)*u+p
      Flux(3) = Q(3)*u
      Flux(4) = rho*H*u
    else
      ! build Flux in y direction
      Flux(1) = Q(3)
      Flux(2) = Q(2)*v
      Flux(3) = Q(3)*v+p
      Flux(4) = rho*H*v
    end if

  end function Euler_buildFlux
  
  
  
  ! This function returns the Roe mean values in PRIMITIVE
  ! variables + the speed of sound waves
  ! Qroe(1) = rho
  ! Qroe(2) = u
  ! Qroe(3) = v
  ! Qroe(4) = H
  ! Qroe(5) = c
  function Euler_calculateQroec(Ql, Qr) result(Qroe)

    ! The left and right solution values
    real(DP), dimension(4), intent(IN)		:: Ql, Qr

    ! The computed Roe values
    real(DP), dimension(5)					:: Qroe

    ! temp variables
    real(DP)		:: rhol, rhor, denom, Hl, Hr, El, Er, pl, pr, velnorm

    ! Gamma
    real(dp) :: gamma = 1.4_dp
    
    
    
    ! Get densities
    rhol = Ql(1)
    rhor = Qr(1)
    
    ! Calculate auxiliary variable
    denom = (sqrt(rhol)+sqrt(rhor))
    
    ! Set Roe-density
    Qroe(1) = sqrt(rhol*rhor)
    
    ! Set Roe-x-velocity
    Qroe(2) = (Ql(2)/sqrt(rhol) + Qr(2)/sqrt(rhor))/denom
    
    ! Set Roe-y-velocity
    Qroe(3) = (Ql(3)/sqrt(rhol) + Qr(3)/sqrt(rhor))/denom

    ! Get left and right energy states
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    
    ! Calculate left and right pressure
    pl = (gamma-1.0_dp)*rhol*(El-0.5_dp*( (Ql(2)/rhol)**2.0_dp + (Ql(3)/rhol)**2.0_dp ) )
    pr = (gamma-1.0_dp)*rhor*(Er-0.5_dp*( (Qr(2)/rhor)**2.0_dp + (Qr(3)/rhor)**2.0_dp ) )
    
    ! Calculate left and right stagnation enthalpy
    Hl = El + pl/rhol
    Hr = Er + pr/rhor
    
    ! Calculate Roe-stagnation enthalpy
    Qroe(4) = (sqrt(rhol)*Hl + sqrt(rhor)*Hr)/denom
    
    ! Calculate the speed of sound for the Roe-values
    Qroe(5) = sqrt( (gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ) )
    
  end function Euler_calculateQroec
    
    
    
  
  
  
  ! This routine builds the right eigenvectors for the mixed jacobian
  function Euler_buildMixedRcfromRoe(Q,a,b) result(R)

    ! Right eigenvectors
    real(DP), dimension(4,4)	:: R

    ! The solution components q1 = roe, q2 = u, q3 = v, q4 = H, q5 = c
    real(DP), dimension(5), intent(IN)		:: Q
    
    ! Components of the normal vector
    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: rho, u, v, H, c, ve, velnorm
   
    rho = Q(1)
    u = Q(2)
    v = Q(3)
    H = Q(4)
    c = Q(5)
    ve = a*u+b*v
    velnorm = 0.5_dp*sqrt(u*u+v*v)
    
        
      ! Build matrix of right eigenvectors
      R(1,1) = 1.0_dp
      R(2,1) = u-c*a
      R(3,1) = v-c*b
      R(4,1) = H-c*ve
      R(1,2) = 1.0_dp
      R(2,2) = u
      R(3,2) = v
      R(4,2) = velnorm
      R(1,3) = 1.0_dp
      R(2,3) = u+c*a
      R(3,3) = v+c*b
      R(4,3) = H+c*ve
      R(1,4) = 0.0_dp
      R(2,4) = b
      R(3,4) = -a
      R(4,4) = u*b-v*a
    
  end function Euler_buildMixedRcfromRoe
  
  
  
  ! This routine builds the diagonalmatrix of the absolut value of the eigenvalues
  function Euler_buildMixedaLambdacfromRoe(Q,a,b) result(aLambda)

    ! Left eigenvectors
    real(DP), dimension(4,4)	:: aLambda

    ! The solution components q1 = roe, q2 = u, q3 = v, q4 = H, q5 = c
    real(DP), dimension(5), intent(IN)		:: Q

    ! Components of the normal vector
    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: u, v, c, ve
    
    ! Constant Gamma
    real(dp) :: gamma = 1.4_dp

    u = Q(2)
    v = Q(3)
    c = Q(5)
    ve = a*u+b*v
    
    ! Build matrix with the absolute values of the eigenvalues
    aLambda = 0.0_dp
    aLambda(1,1) = abs(ve-c)
    aLambda(2,2) = abs(ve)
    aLambda(3,3) = abs(ve+c)
    aLambda(4,4) = abs(ve)
    
  end function 


  ! This routine builds the left eigenvectors for the mixed jacobian
  function Euler_buildMixedLcfromRoe(Q,a,b) result(L)

    ! Left eigenvectors
    real(DP), dimension(4,4)	:: L

    ! The solution components q1 = roe, q2 = u, q3 = v, q4 = H, q5 = c
    real(DP), dimension(5), intent(IN)		:: Q

    ! Components of the normal vector
    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: rho, u, v, H, c, ve, velnorm, b1, b2
    
    ! Constant Gamma
    real(dp) :: gamma = 1.4_dp
   
    rho = Q(1)
    u = Q(2)
    v = Q(3)
    H = Q(4)
    c = Q(5)
    ve = a*u+b*v
    velnorm = 0.5_dp*sqrt(u*u+v*v)
    b2 = (gamma-1.0_dp)/c/c
    b1= b2*velnorm
    
    
       ! Build matrix of left eigenvectors
       L(1,1) = 0.5_dp*(b1+ve/c)
       L(2,1) = 1.0_dp-b1
       L(3,1) = 0.5_dp*(b1-ve/c)
       L(4,1) = a*ve-b*u
       L(1,2) = 0.5_dp*(-b2*u-a/c)
       L(2,2) = b2*u
       L(3,2) = 0.5_dp*(-b2*u+a/c)
       L(4,2) = b
       L(1,3) = 0.5_dp*(-b2*v-b/c)
       L(2,3) = b2*v
       L(3,3) = 0.5_dp*(-b2*v+b/c)
       L(4,3) = -a
       L(1,4) = 0.5_dp*b2
       L(2,4) = -b2
       L(3,4) = 0.5_dp*b2
       L(4,4) = 0.0_dp
           
  end function 
  
  
  
  
  ! This routine returns the eigenvalues of the jacobi matrix
  function Euler_buildEigenvaluescfromRoe(Q,a,b) result(Eigenvalues)

    real(DP), dimension(4)	:: Eigenvalues

    ! The solution components q1 = roe, q2 = u, q3 = v, q4 = H, q5 = c
    real(DP), dimension(5), intent(IN)		:: Q

    ! Components of the normal vector
    real(dp), intent(IN)                     :: a,b
    
    ! Local variables
    real(dp) :: u, v, c, ve
    
    ! Constant Gamma
    real(dp) :: gamma = 1.4_dp

    u = Q(2)
    v = Q(3)
    c = Q(5)
    ve = a*u+b*v
    
    ! build eigenvalues
    Eigenvalues(1) = ve-c
    Eigenvalues(2) = ve
    Eigenvalues(3) = ve+c
    Eigenvalues(4) = ve

  end function Euler_buildEigenvaluescfromRoe 
  
  
  
  
  
  
  ! This function transforms a vector of conserved variables
  ! Q(1) = rho
  ! Q(2) = rho*u
  ! Q(3) = rho*v
  ! Q(4) = rho*E
  ! into a vector consisting of primitive variables + the speed of sound waves
  ! Q(1) = rho
  ! Q(2) = u
  ! Q(3) = v
  ! Q(4) = H
  ! Q(5) = c
  function Euler_transformVector(Qin) result(Qout)

    ! The left and right solution values
    real(DP), dimension(4), intent(IN)		:: Qin

    ! The computed Roe values
    real(DP), dimension(5)					:: Qout

    ! temp variables
    real(DP)		:: rho, H, E, p

    ! Gamma
    real(dp) :: gamma = 1.4_dp
    
    
    ! Get density
    rho = Qin(1)

    ! Set Roe-density
    Qout(1) = rho
    
    ! Set x-velocity
    Qout(2) = Qin(2)/rho
    
    ! Set y-velocity
    Qout(3) = Qin(3)/rho

    ! Get energy state
    E = Qin(4)/rho
    
    ! Calculate pressure
    p = (gamma-1.0_dp)*rho*(E-0.5_dp*(Qout(2)*Qout(2)+Qout(3)*Qout(3) ) )
    
    ! Calculate stagnation enthalpy
    H = E + p/rho
    
    ! Save stagnation enthalpy
    Qout(4) = H
    
    ! Calculate the speed of sound for the Roe-values
    !Qout(5) = sqrt( (gamma-1.0_dp)*(H - 0.5_dp*(Qout(2)*Qout(2)+Qout(3)*Qout(3)) ) )
    Qout(5) = sqrt( max(sys_EPSREAL*sys_EPSREAL , gamma*p/rho) )

  end function Euler_transformVector
    
    
    
    
    
  ! This routine returns the flux vector for the 2d compressible euler
  ! equations of gas dynamics in direction d
  ! d=1: x-direction, d=2: y-direction
  function Euler_buildFlux_HLL2D(Qii,Qaa,a,b) result(Flux)

    ! The flux vector in direction d at Q
    real(DP), dimension(4)	:: Flux

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(4), intent(IN)		:: Qii, Qaa

    ! the direction: d=1: x-direction, d=2: y-direction
    real(dp), intent(IN)                     :: a, b
    
    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(4)		:: Qi, Qa, FX

    ! pressure, stagnation enthalpy
    real(DP)                                :: pi, pa, Hi, Ha, Ei, Ea, ui, ua, vi, va, rhoi, rhoa, t1, t2, t3, SL, SR

    ! temporary variable
    real(DP)                                :: rho, u, v, E, ci, ca, denom, H, c
    
    ! Constant Gamma
    real(dp) :: gamma = 1.4_dp
    
    Qi=Qii
    Qa=Qaa
    
    ! Rotate inner trace
    t2 = Qi(2)
    t3 = Qi(3)
    Qi(2) =  a*t2 + b*t3
    Qi(3) = -b*t2 + a*t3
    
    ! Rotate outer trace
    t2 = Qa(2)
    t3 = Qa(3)
    Qa(2) =  a*t2 + b*t3
    Qa(3) = -b*t2 + a*t3

    ! Calculate primitive variables
    rhoi=Qi(1)
    ui=Qi(2)/rhoi
    vi=Qi(3)/rhoi
    Ei=Qi(4)/rhoi
    
    rhoa=Qa(1)
    ua=Qa(2)/rhoa
    va=Qa(3)/rhoa
    Ea=Qa(4)/rhoa

    ! Compute the pressure
    pi = (gamma - 1.0_dp)*rhoi*(Ei-0.5_dp*(ui*ui+vi*vi))
    pa = (gamma - 1.0_dp)*rhoa*(Ea-0.5_dp*(ua*ua+va*va))
    
    ! Compute H, the stagnation enthalpy
    Hi = Ei + pi/rhoi
    Ha = Ea + pa/rhoa
    
    ! Compute speed of sound
    ci = sqrt( (gamma-1.0_dp)*(Hi - 0.5_dp*(ui*ui+vi*vi) ) )
    ca = sqrt( (gamma-1.0_dp)*(Ha - 0.5_dp*(ua*ua+va*va) ) )
    
    ! Compute Roe-average variables
    
    
    
    
    ! Calculate auxiliary variable
    denom = (sqrt(rhoi)+sqrt(rhoa))
    
    ! Set Roe-density
    rho = sqrt(rhoi*rhoa)
    
    ! Set Roe-x-velocity
    u = (ui*sqrt(rhoi) + ua*sqrt(rhoa))/denom
    
    ! Set Roe-y-velocity
    v = (vi*sqrt(rhoi) + va*sqrt(rhoa))/denom
    
    ! Calculate Roe-stagnation enthalpy
    H = (sqrt(rhoi)*Hi + sqrt(rhoa)*Ha)/denom
    
    ! Calculate the speed of sound for the Roe-values
    c = sqrt( (gamma-1.0_dp)*(H - 0.5_dp*(u*u+v*v) ) )
    
    
    ! Compute estimate wave speeds
    SL = min(ui - ci,u-c)
    SR = max(ua + ca,u+c)
    
    ! Compute HLL flux
    t1 = (min(SR,0.0_dp)-min(0.0_dp,SL))/(SR-SL)
    t2 = 1.0_dp - t1
    t3 = (SR*abs(SL)-SL*abs(SR))/(2.0_dp*(SR-SL))
    
    FX = t1*Euler_buildFlux(Qa,1) + t2*Euler_buildFlux(Qi,1) - t3*(Qa-Qi)
    
    ! Rotate back
    Flux(1) = FX(1)
    Flux(2) = a*FX(2)-b*FX(3)
    Flux(3) = b*FX(2)+a*FX(3)
    Flux(4) = FX(4)
    
    

  end function Euler_buildFlux_HLL2D
    
    
    
end module
