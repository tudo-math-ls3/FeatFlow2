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

    if (abs(a)<10*SYS_EPSREAL_DP) then

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

  end function buildMixedR



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

  end function buildMixedaLambda


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

  end function buildMixedL

















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


  end function buildMixedR2



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

  end function buildMixedaLambda2


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

  end function buildMixedL2




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


  end function buildMixedR2c



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

  end function buildMixedaLambda2c


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

  end function buildMixedL2c




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

!    ! Calculate primitive variables
    rho=Q(1)
    u=Q(2)/rho
    v=Q(3)/rho
!    E=Q(4)/rho
!
!    ! Compute the pressure
!    p = (gamma - 1.0_dp)*rho*(E-0.5_dp*(u*u+v*v))
!
!    ! Compute H, the stagnation enthalpy
!    H = E + p/rho

    if (d==1) then
!       ! build Flux in x direction
!       Flux(1) = Q(2)
!       Flux(2) = Q(2)*u+p
!       Flux(3) = Q(3)*u
!       Flux(4) = Q(2)*H


       ! Stolen from Matthias :)
       Flux(1) = Q(2)
       Flux(2) = (gamma - 1.0_dp)*Q(4)-0.5_dp*(gamma-3.0_dp)*u*Q(2)-0.5_dp*(gamma-1.0_dp)*v*Q(3)
       Flux(3) = Q(3)*u
       Flux(4) = (gamma*Q(4)-0.5_dp*(gamma-1.0_dp)*(u*Q(2)+v*Q(3)))*u


    else
!       ! build Flux in y direction
!       Flux(1) = Q(3)
!       Flux(2) = Q(2)*v
!       Flux(3) = Q(3)*v+p
!       Flux(4) = Q(3)*H

       ! Stolen from Matthias :)
       Flux(1) = Q(3)
       Flux(2) = Q(3)*u
       Flux(3) = (gamma - 1.0_dp)*Q(4)-0.5_dp*(gamma-3.0_dp)*v*Q(3)-0.5_dp*(gamma-1.0_dp)*u*Q(2)
       Flux(4) = (gamma*Q(4)-0.5_dp*(gamma-1.0_dp)*(u*Q(2)+v*Q(3)))*v


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
    real(DP)		:: rhol, rhor, denom, Hl, Hr, El, Er, pl, pr, velnorm, aux, ul, ur, vl, vr

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
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )





    !    ! Stolen from Matthias :)
    !    ! Compute Roe mean values
    !    ul = Ql(2)/Ql(1)
    !    ur = Qr(2)/Qr(1)
    !    vl = Ql(3)/Ql(1)
    !    vr = Qr(3)/Qr(1)
    !    Qroe(1) = sqrt(rhol*rhor)
    !    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    !    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    !    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    !    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    !    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    !    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)






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

  end function Euler_buildMixedaLambdacfromRoe
  
  
  function Euler_buildMixedalambda (Q,a,b) result(aLambda)

    ! Absolute of max eigenvalue
    real(DP):: alambda

    ! The solution components
    real(DP), dimension(4), intent(IN)		:: Q

    ! Components of the normal vector
    real(dp), intent(IN)                     :: a,b

    ! Local variables
    real(dp) :: u, v, c, ve, rho, E

    ! Constant Gamma
    real(dp) :: gamma = 1.4_dp

    rho = Q(1)
    u = Q(2)/rho
    v = Q(3)/rho
    E = Q(4)/rho
    ve = a*u+b*v
    c = sqrt(max(1.4_dp*0.4_dp*(E-0.5_dp*(u*u+v*v)),0.0_dp))

    ! Build absolute value of largest eigenvalue
    alambda = abs(ve) + c
    !alambda = sqrt(u*u+v*v) + c
    
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
    L(4,1) = a*v-b*u
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

  end function Euler_buildMixedLcfromRoe




  ! This routine builds the Euler-Jacobian in x direction
  function Euler_buildJacobixcfromRoe(Q) result(Jx)

    ! Left eigenvectors
    real(DP), dimension(4,4)	:: Jx

    ! The solution components q1 = roe, q2 = u, q3 = v, q4 = H, q5 = c
    real(DP), dimension(5), intent(IN)		:: Q

    ! Local variables
    real(dp) :: rho, u, v, H, c, E

    ! Constant Gamma
    real(dp) :: gamma = 1.4_dp

    rho = Q(1)
    u = Q(2)
    v = Q(3)
    H = Q(4)
    c = Q(5)
    E = (H+(gamma-1.0_dp)*0.5_dp*(u*u+v*v))/gamma



    ! Build matrix
    Jx(1,1) = 0.0_dp
    Jx(2,1) = 0.5_dp*(gamma-3.0_dp)*u*u+0.5_dp*(gamma-1.0_dp)*v*v
    Jx(3,1) = -u*v
    Jx(4,1) = -gamma*u*E+(gamma-1.0_dp)*u*(u*u+v*v)
    Jx(1,2) = 1.0_dp
    Jx(2,2) = (3.0_dp-gamma)*u
    Jx(3,2) = v
    Jx(4,2) = gamma*E-0.5_dp*(gamma-1.0_dp)*(v*v+3.0_dp*u*u)
    Jx(1,3) = 0.0_dp
    Jx(2,3) = -(gamma-1.0_dp)*v
    Jx(3,3) = u
    Jx(4,3) = -(gamma-1.0_dp)*u*v
    Jx(1,4) = 0.0_dp
    Jx(2,4) = (gamma-1.0_dp)
    Jx(3,4) = 0.0_dp
    Jx(4,4) = gamma*u

  end function Euler_buildJacobixcfromRoe



  ! This routine builds the Euler-Jacobian in y direction
  function Euler_buildJacobiycfromRoe(Q) result(Jy)

    ! Left eigenvectors
    real(DP), dimension(4,4)	:: Jy

    ! The solution components q1 = roe, q2 = u, q3 = v, q4 = H, q5 = c
    real(DP), dimension(5), intent(IN)		:: Q

    ! Local variables
    real(dp) :: rho, u, v, H, c, E

    ! Constant Gamma
    real(dp) :: gamma = 1.4_dp

    rho = Q(1)
    u = Q(2)
    v = Q(3)
    H = Q(4)
    c = Q(5)
    E = (H+(gamma-1.0_dp)*0.5_dp*(u*u+v*v))/gamma


    ! Build matrix
    Jy(1,1) = 0.0_dp
    Jy(2,1) = -u*v
    Jy(3,1) = 0.5_dp*(gamma-3.0_dp)*v*v+0.5_dp*(gamma-1.0_dp)*u*u
    Jy(4,1) = -gamma*v*E+(gamma-1.0_dp)*v*(u*u+v*v)
    Jy(1,2) = 0.0_dp
    Jy(2,2) = v
    Jy(3,2) = -(gamma-1.0_dp)*u
    Jy(4,2) = -(gamma-1.0_dp)*u*v
    Jy(1,3) = 1.0_dp
    Jy(2,3) = u
    Jy(3,3) = (3.0_dp-gamma)*v
    Jy(4,3) = gamma*E-0.5_dp*(gamma-1.0_dp)*(u*u+3.0_dp*v*v)
    Jy(1,4) = 0.0_dp
    Jy(2,4) = 0.0_dp
    Jy(3,4) = (gamma-1.0_dp)
    Jy(4,4) = gamma*v

  end function Euler_buildJacobiycfromRoe



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
    Qout(5) = sqrt( max(sys_EPSREAL_DP*sys_EPSREAL_DP , gamma*p/rho) )

  end function Euler_transformVector





  ! This routine returns the HLL flux for the 2d compressible euler
  ! equations of gas dynamics
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
    ci = sqrt(max( (gamma-1.0_dp)*(Hi - 0.5_dp*(ui*ui+vi*vi) ),0.0_dp) )
    ca = sqrt(max( (gamma-1.0_dp)*(Ha - 0.5_dp*(ua*ua+va*va) ),0.0_dp) )

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
    c = sqrt(max( (gamma-1.0_dp)*(H - 0.5_dp*(u*u+v*v) ),0.0_dp) )


    ! Compute estimate wave speeds
    SL = min(ui - ci,u-c)
    SR = max(ua + ca,u+c)

    if (SL.ge.SR) then
       write(*,*) 'Warning'
    end if

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




  ! This routine returns the HLLC flux for the 2d compressible euler
  ! equations of gas dynamics
  function Euler_buildFlux_HLLC2D(Qii,Qaa,a,b) result(Flux)

    ! The flux vector in direction d at Q
    real(DP), dimension(4)	:: Flux

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(4), intent(IN)		:: Qii, Qaa

    ! the direction: d=1: x-direction, d=2: y-direction
    real(dp), intent(IN)                     :: a, b

    ! The solution components q1 = h, q2 = uh, q3 = vh
    real(DP), dimension(4)		:: Qi, Qa, FX

    ! pressure, stagnation enthalpy
    real(DP)                                :: pi, pa, Hi, Ha, Ei, Ea, ui, ua, vi, va, rhoi, rhoa, t1, t2, t3, SL, SR, SS

    ! temporary variable
    real(DP)                                :: rho, u, v, E, ci, ca, denom, H, c, ps, us, ql, qr

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
    ci = sqrt( max((gamma-1.0_dp)*(Hi - 0.5_dp*(ui*ui+vi*vi) ),0.0_dp ))
    ca = sqrt( max((gamma-1.0_dp)*(Ha - 0.5_dp*(ua*ua+va*va) ),0.0_dp ))

    ci = max(ci,1.0e-12)
    ca = max(ca,1.0e-12)



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
    c = sqrt(max( (gamma-1.0_dp)*(H - 0.5_dp*(u*u+v*v) ),0.0_dp) )




    !    !!! Type 1
    !
    !    ! Compute estimate wave speeds
    !    SL = min(ui - ci,ua-ca)
    !    SR = max(ua + ca,ui+ci)
    !    SS = 0.5_dp*(ui+ua) + (pi-pa)/(0.25_dp*(rhoi+rhoa)*(ci+ca))
    !
    !!    SS = max(SS,SL)
    !!    SS = min(SS,SR)
    !!
    !!    SS = 0.5_dp*(SL+SR)
    !
    !
    !
    !
    !
    !    if(SL>SR) then
    !    write(*,*) 'Warning**************'
    !    end if
    !    if (SL>SS) then
    !      write(*,*) 'Warning'
    !    end if
    !    if (SS>SR) then
    !      write(*,*) 'Warning'
    !    end if
    !
    !    FX = Euler_buildFlux(Qa+Qi,1)
    !
    !    if (SL.ge.0.0_dp) then
    !      FX = Euler_buildFlux(Qi,1)
    !    elseif ((SL<0.0_dp).and.(SS.ge.0.0_dp)) then
    !      FX = Euler_buildFlux(Qi,1) + SL*( rhoi*(SL-ui)/min((SL-SS),-SYS_EPSREAL_DP)*(/ 1.0_dp, SS, vi, Ei/rhoi+(SS-ui)*(SS+pi/(rhoi*(SL-ui))) /) -Qi)
    !    elseif ((SS<0.0_dp).and.(SR>0.0_dp)) then
    !      FX = Euler_buildFlux(Qa,1) + SR*( rhoa*(SR-ua)/max((SR-SS),SYS_EPSREAL_DP)*(/ 1.0_dp, SS, va, Ea/rhoa+(SS-ua)*(SS+pa/(rhoa*(SR-ua))) /) -Qa)
    !    elseif (SR.le.0.0_dp) then
    !      FX = Euler_buildFlux(Qa,1)
    !    end if






    !    !!! Type 2 !!!
    !
    !
    !    ps = 0.5_dp*(pi+pa)+0.5_dp*(ui-ua)*0.5_dp*(rhoi+rhoa)*0.5_dp*(ci+ca)
    !    us = 0.5_dp*(ui+ua)+0.5_dp*(pi-pa)/(0.5_dp*(rhoi+rhoa)*0.5_dp*(ci+ca))
    !
    !    if (ps>pi) then
    !      ql = sqrt(max(1.0_dp+(gamma+1.0_dp)/(2.0_dp*gamma)*(ps/pi-1.0_dp),0.0_dp))
    !    else
    !      ql = 1.0_dp
    !    end if
    !
    !    if (ps>pa) then
    !      qr = sqrt(max(1.0_dp+(gamma+1.0_dp)/(2.0_dp*gamma)*(ps/pa-1.0_dp),0.0_dp))
    !    else
    !      qr = 1.0_dp
    !    end if
    !
    !
    !
    !
    !    ! Compute estimate wave speeds
    !    SL = ui-ci*ql
    !    SR = ua+ca*qr
    !    SS = us
    !
    !!    SS = max(SS,SL)
    !!    SS = min(SS,SR)
    !
    !
    !
    !    if (SL.ge.0.0_dp) then
    !      FX = Euler_buildFlux(Qi,1)
    !    elseif ((SL<0.0_dp).and.(SS.ge.0.0_dp)) then
    !      FX = Euler_buildFlux(Qi,1) + SL*( rhoi*(SL-ui)/min((SL-SS),-SYS_EPSREAL_DP)*(/ 1.0_dp, SS, vi, Ei/rhoi+(SS-ui)*(SS+pi/(rhoi*(SL-ui))) /) -Qi)
    !    elseif ((SS<0.0_dp).and.(SR>0.0_dp)) then
    !      FX = Euler_buildFlux(Qa,1) + SR*( rhoa*(SR-ua)/max((SR-SS),SYS_EPSREAL_DP)*(/ 1.0_dp, SS, va, Ea/rhoa+(SS-ua)*(SS+pa/(rhoa*(SR-ua))) /) -Qa)
    !    elseif (SR.le.0.0_dp) then
    !      FX = Euler_buildFlux(Qa,1)
    !    else
    !      FX = 0.5_dp*Euler_buildFlux(0.5_dp*(Qi+Qa),1)
    !    end if
    !
    !
    !    if(SL>SR) then
    !      write(*,*) 'Warning**************'
    !    end if
    !    if (SL>SS) then
    !      write(*,*) 'Warning'
    !    end if
    !    if (SS>SR) then
    !      write(*,*) 'Warning'
    !    end if



    ! Type 3
    SL = min(ui - ci,ua-ca)
    SR = max(ua + ca,ui+ci)
    SS = ( rhoa*ua*(SR-ua)-rhoi*ui*(SL-ui)+pi-pa )/(rhoa*(SR-ua)-rhoi*(SL-ui))

    ps = rhoi*(SL-ui)*(SS-ui)+pi


    FX = 0.5_dp*(Euler_buildFlux(Qi,1)+Euler_buildFlux(Qa,1) -ui*Qi+(/ 0.0_dp,ps-pi,0.0_dp,ps*SS-pi*ui /) &
         -ua*Qa+(/ 0.0_dp,ps-pa,0.0_dp,ps*SS-pa*ua /) )






    ! Rotate back
    Flux(1) = FX(1)
    Flux(2) = a*FX(2)-b*FX(3)
    Flux(3) = b*FX(2)+a*FX(3)
    Flux(4) = FX(4)



  end function Euler_buildFlux_HLLC2D








  !*****************************************************************************
  !* -- Roe's Flux Function ---
  !*
  !* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
  !* Schemes, Journal of Computational Physics, 43, pp. 357-372.
  !*
  !* Katate Masatsuka, February 2009. http://www.cfdbooks.com
  !*****************************************************************************
  function Roe(uL, uR, nx, ny)
    real(dp) :: uL(4), uR(4) !  Input: conservative variables rho*[1, u, v, E]
    real(dp) :: nx, ny       !  Input: face normal vector, [nx, ny] (Left-to-Right)
    real(dp) :: Roe(4)       ! Output: Roe flux function (upwind)
    !Local constants
    real(dp) :: gamma                          ! Ratio of specific heat.
    real(dp) :: zero, fifth, half, one, two    ! Numbers
    !Local variables
    real(dp) :: tx, ty       ! Tangent vector (perpendicular to the face normal)
    real(dp) :: vxL, vxR, vyL, vyR             ! Velocity components.
    real(dp) :: rhoL, rhoR, pL, pR             ! Primitive variables.
    real(dp) :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
    real(dp) :: aL, aR, HL, HR                 ! Speeds of sound.
    real(dp) :: RT,rho,vx,vy,H,a,vn, vt        ! Roe-averages
    real(dp) :: drho,dvx,dvy,dvn,dvt,dpp,dV(4)  ! Wave strenghs
    real(dp) :: ws(4),dws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
    real(dp) :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
    integer :: i, j

    !Constants.
    gamma = 1.4_dp
    zero = 0.0_dp
    fifth = 0.2_dp
    half = 0.5_dp
    one = 1.0_dp
    two = 2.0_dp

    !Tangent vector (Do you like it? Actually, Roe flux can be implemented
    ! without any tangent vector. See "I do like CFD, VOL.1" for details.)
    tx = -ny
    ty = nx

    !Primitive and other variables.
    !  Left state
    rhoL = uL(1)
    vxL = uL(2)/uL(1)
    vyL = uL(3)/uL(1)
    vnL = vxL*nx+vyL*ny
    vtL = vxL*tx+vyL*ty
    pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
    aL = sqrt(gamma*pL/rhoL)
    HL = ( uL(4) + pL ) / rhoL
    !  Right state
    rhoR = uR(1)
    vxR = uR(2)/uR(1)
    vyR = uR(3)/uR(1)
    vnR = vxR*nx+vyR*ny
    vtR = vxR*tx+vyR*ty
    pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
    aR = sqrt(gamma*pR/rhoR)
    HR = ( uR(4) + pR ) / rhoR

    !First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
    rho = RT*rhoL
    vx = (vxL+RT*vxR)/(one+RT)
    vy = (vyL+RT*vyR)/(one+RT)
    H = ( HL+RT* HR)/(one+RT)
    a = sqrt( (gamma-one)*(H-half*(vx*vx+vy*vy)) )
    vn = vx*nx+vy*ny
    vt = vx*tx+vy*ty

    !Wave Strengths
    drho = rhoR - rhoL
    dpp = pR - pL
    dvn = vnR - vnL
    dvt = vtR - vtL

    dV(1) = (dpp - rho*a*dvn )/(two*a*a)
    dV(2) = rho*dvt/a
    dV(3) = drho - dpp/(a*a)
    dV(4) = (dpp + rho*a*dvn )/(two*a*a)

    !Wave Speed
    ws(1) = abs(vn-a)
    ws(2) = abs(vn)
    ws(3) = abs(vn)
    ws(4) = abs(vn+a)

    !Harten's Entropy Fix JCP(1983), 49, pp357-393:
    ! only for the nonlinear fields.
    dws(1) = fifth
    if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
    dws(4) = fifth
    if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

    !Right Eigenvectors
    Rv(1,1) = one
    Rv(2,1) = vx - a*nx
    Rv(3,1) = vy - a*ny
    Rv(4,1) =  H - vn*a

    Rv(1,2) = zero
    Rv(2,2) = a*tx
    Rv(3,2) = a*ty
    Rv(4,2) = vt*a

    Rv(1,3) = one
    Rv(2,3) = vx
    Rv(3,3) = vy
    Rv(4,3) = half*(vx*vx+vy*vy)

    Rv(1,4) = one
    Rv(2,4) = vx + a*nx
    Rv(3,4) = vy + a*ny
    Rv(4,4) =  H + vn*a

    !Dissipation Term
    diss = zero
    do i=1,4
       do j=1,4
          diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
       end do
    end do

    !Compute the flux.
    fL(1) = rhoL*vnL
    fL(2) = rhoL*vnL * vxL + pL*nx
    fL(3) = rhoL*vnL * vyL + pL*ny
    fL(4) = rhoL*vnL *  HL

    fR(1) = rhoR*vnR
    fR(2) = rhoR*vnR * vxR + pR*nx
    fR(3) = rhoR*vnR * vyR + pR*ny
    fR(4) = rhoR*vnR *  HR

    Roe = half * (fL + fR - diss)

  end function Roe
  
  
  
  function DRoe(uL, uR, nx, ny, h)
    real(dp) :: uL(4), uR(4)             !  Input: conservative variables rho*[1, u, v, E]
    real(dp) :: nx, ny                   !  Input: face normal vector, [nx, ny]
    real(dp) :: h                        !  Input: infinitesimal constant
    real(dp), dimension(4,8) :: DRoe     ! Output: Approximate differential of Roe flux function
    
    real(dp), dimension(4) :: F, U
    integer :: i
    
    F = Roe(uL, uR, nx, ny)
    
    ! First order approx
    do i = 1,4
      U = uL
      U(i) = U(i) + h
      DRoe(:,i) = (Roe(U, uR, nx, ny) - F)/h
    end do
    
    do i = 1,4
      U = uR
      U(i) = U(i) + h
      DRoe(:,i+4) = (Roe(uL, U, nx, ny) - F)/h
    end do

!    ! Second order approx
!    do i = 1,4
!      U = uL
!      U(i) = U(i) + h
!      DRoe(:,i) = Roe(U, uR, nx, ny)
!      U = uL
!      U(i) = U(i) - h
!      DRoe(:,i) = 0.5_dp*(DRoe(:,i)-Roe(U, uR, nx, ny))/h
!    end do
!    
!    do i = 1,4
!      U = uR
!      U(i) = U(i) + h
!      DRoe(:,i+4) = Roe(uL, U, nx, ny)
!      U = uR
!      U(i) = U(i) - h
!      DRoe(:,i+4) = 0.5_dp*(DRoe(:,i+4)-Roe(uL, U, nx, ny))/h
!    end do
    
  end function Droe


  !*****************************************************************************
  !* -- Rotated-Roe-HLL Flux Function ---
  !*
  !* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
  !* Resolving, Rotated-Hybrid Riemann Solvers,
  !* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
  !*
  !* The most robust Riemann solver known to the author (in terms of nonlinear
  !* instability such as carbuncle).
  !*
  !* NB: At a boundary face, need to switch to a geometric normal vector:
  !*               (nx2,ny2)=(nx, ny), (nx1,ny1)=(-ny,nx).
  !*     This is not implemented here. It requires information on whether
  !*     the geometric normal, (nx,ny), is on a boundary face or not.
  !*     It shouldn't be difficult for you to implement it.
  !*
  !* Katate Masatsuka, February 2010. http://www.cfdbooks.com
  !*****************************************************************************
  function Rotated_RHLL(uL, uR, nx, ny)
    real(dp) :: uL(4), uR(4)    !  Input: conservative variables rho*[1, u, v, E]
    real(dp) :: nx, ny          !  Input: face normal vector, [nx, ny] (Left-to-Right)
    real(dp) :: Rotated_RHLL(4) ! Output: Rotated_RHLL flux function.
    !Local constants
    real(dp) :: gamma                          ! Ratio of specific heat.
    real(dp) :: zero, fifth, half, one, two    ! Numbers
    real(dp) :: eps                            !
    !Local variables
    real(dp) :: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
    real(dp) :: tx, ty                         ! Tangent vector (taken as n1)
    real(dp) :: alpha1, alpha2                 ! Projections of the new normals
    real(dp) :: vxL, vxR, vyL, vyR             ! Velocity components.
    real(dp) :: rhoL, rhoR, pL, pR             ! Primitive variables.
    real(dp) :: vnL, vnR, vtL, vtR             ! Normal and tagent velocities
    real(dp) :: aL, aR, HL, HR                 ! Speeds of sound and total enthalpy
    real(dp) :: RT,rho,vx,vy,H,a               ! Roe-averages
    real(dp) :: vn, vt                         ! Normal and tagent velocities(Roe-average)
    real(dp) :: drho,dvx,dvy,dvn,dvt,dpp,dV(4)  ! Wave strenghs
    real(dp) :: abs_dq                         ! Magnitude of the velocity difference
    real(dp) :: abs_ws(4),ws(4),dws(4), Rv(4,4)! Wave speeds and right-eigevectors
    real(dp) :: SRp,SLm                        ! Wave speeds for the HLL part
    real(dp) :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
    real(dp) :: temp
    integer :: i, j

    !Constants.
    gamma = 1.4_dp
    zero = 0.0_dp
    fifth = 0.2_dp
    half = 0.5_dp
    one = 1.0_dp
    two = 2.0_dp
    eps = 1.0e-12_dp ! 1.0e-12 in the original paper (double precision)

    !Primitive and other variables.
    !  Left state
    rhoL = uL(1)
    vxL = uL(2)/uL(1)
    vyL = uL(3)/uL(1)
    pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
    pL = max(0.0_dp,pL)
    aL = sqrt(gamma*pL/rhoL)
    HL = ( uL(4) + pL ) / rhoL
    !  Right state
    rhoR = uR(1)
    vxR = uR(2)/uR(1)
    vyR = uR(3)/uR(1)
    pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
    pR = max(0.0_dp,pR)
    aR = sqrt(gamma*pR/rhoR)
    HR = ( uR(4) + pR ) / rhoR

    vnL = vxL*nx + vyL*ny
    vnR = vxR*nx + vyR*ny

    !Compute the flux.
    fL(1) = rhoL*vnL
    fL(2) = rhoL*vnL * vxL + pL*nx
    fL(3) = rhoL*vnL * vyL + pL*ny
    fL(4) = rhoL*vnL *  HL

    fR(1) = rhoR*vnR
    fR(2) = rhoR*vnR * vxR + pR*nx
    fR(3) = rhoR*vnR * vyR + pR*ny
    fR(4) = rhoR*vnR *  HR

    !Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
    !(NB: n1 and n2 may need to be frozen at some point during
    !     a steady calculation to fully make it converge. For time-accurate
    !     calculation, this is fine.)
    ! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    abs_dq = sqrt( (vxR-vxL)**2+(vyR-vyL)**2 )
    if ( abs_dq > eps) then
       nx1 = (vxR-vxL)/abs_dq
       ny1 = (vyR-vyL)/abs_dq
    else
       nx1 = -ny
       ny1 =  nx
    endif
    alpha1 = nx * nx1 + ny * ny1
    !   To make alpha1 always positive.
    temp = sign(one,alpha1)
    nx1 = temp * nx1
    ny1 = temp * ny1
    alpha1 = temp * alpha1

    ! Take n2 as perpendicular to n1.
    nx2 = -ny1
    ny2 =  nx1
    alpha2 = nx * nx2 + ny * ny2
    !   To make alpha2 always positive.
    temp = sign(one,alpha2)
    nx2 = temp * nx2
    ny2 = temp * ny2
    alpha2 = temp * alpha2

    !Now we are going to compute the Roe flux with n2 as the normal
    !and n1 as the tagent vector, with modified wave speeds (5.12)

    !Compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
    rho = RT*rhoL
    vx = (vxL+RT*vxR)/(one+RT)
    vy = (vyL+RT*vyR)/(one+RT)
    H = ( HL+RT* HR)/(one+RT)
    a = sqrt( max(0.0_dp,(gamma-one)*(H-half*(vx*vx+vy*vy))) )
    vn = vx*nx2+vy*ny2
    vt = vx*nx1+vy*ny1

    !Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnL = vxL*nx2 + vyL*ny2
    vnR = vxR*nx2 + vyR*ny2
    vtL = vxL*nx1 + vyL*ny1
    vtR = vxR*nx1 + vyR*ny1

    drho = rhoR - rhoL
    dpp =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL

    a = max(a,1.0e-12)

    dV(1) = (dpp - rho*a*dvn )/(two*a*a)
    dV(2) =  rho*dvt/a
    dV(3) =  drho - dpp/(a*a)
    dV(4) = (dpp + rho*a*dvn )/(two*a*a)

    !Wave Speeds for Roe flux part.
    ws(1) = vn-a
    ws(2) = vn
    ws(3) = vn
    ws(4) = vn+a
    abs_ws  = abs(ws)

    !Harten's Entropy Fix JCP(1983), 49, pp357-393:
    !only for the nonlinear fields.
    dws(1) = fifth
    if (abs_ws(1)<dws(1)) abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
    dws(4) = fifth
    if (abs_ws(4)<dws(4)) abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4))

    !HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
    SRp = max( zero, vtR + aR, vt + a)
    SLm = min( zero, vtL - aL, vt - a)

    !Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
    ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + two*alpha1*SRp*SLm )/ (SRp-SLm)

    !Right Eigenvectors: with n2 as normal and n1 as tangent.
    tx = nx1
    ty = ny1

    Rv(1,1) = one
    Rv(2,1) = vx - a*nx2
    Rv(3,1) = vy - a*ny2
    Rv(4,1) =  H - vn*a

    Rv(1,2) = zero
    Rv(2,2) = a*tx
    Rv(3,2) = a*ty
    Rv(4,2) = a*vt

    Rv(1,3) = one
    Rv(2,3) = vx
    Rv(3,3) = vy
    Rv(4,3) = half*(vx*vx+vy*vy)

    Rv(1,4) = one
    Rv(2,4) = vx + a*nx2
    Rv(3,4) = vy + a*ny2
    Rv(4,4) =  H + vn*a

    !Dissipation Term: Roe dissipation with the modified wave speeds.
    diss = zero
    do i=1,4
       do j=1,4
          diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
       end do
    end do

    !Compute the Rotated-RHLL flux.
    Rotated_RHLL = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss

  end function Rotated_RHLL
  
  
  
! This function builds the derivative of the absolute value of the largest eigenvalue

function Euler_buildDlambda(Qll, Qrr, a, b) result(Ddlambda)

    ! The left and right solution values
    real(DP), dimension(4), intent(IN)		:: Qll, Qrr
    
    ! The components of the normal vector
    real(dp), intent(in) :: a, b

    ! The computed derivatives of max abs eigenvalue
    real(DP), dimension(8)					:: Ddlambda

    ! temp variables
    real(DP)		:: rhol, rhor, denom, Hl, Hr, El, Er, pl, pr, velnorm, aux, ul, ur, vl, vr

    ! Roe-values and speed of sound waves
    real(dp), dimension(5) :: Qroe

    ! Gamma
    real(dp) :: gamma = 1.4_dp

    ! Delta
    real(dp) :: ddelta = 0.0000001_dp

    ! Lambda and temp Lambda
    real(dp) :: dlambda, dtlambda
    
    ! Temporary left and right solution values
    real(DP), dimension(4) :: Ql, Qr


    ! Copy solution values to their temporary opposites
    Ql = Qll
    Qr = Qrr

    ! Compute Roe mean values
    rhol = Ql(1)
    rhor = Qr(1)
    ul = Ql(2)/Ql(1)
    ur = Qr(2)/Qr(1)
    vl = Ql(3)/Ql(1)
    vr = Qr(3)/Qr(1)
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    Qroe(1) = sqrt(rhol*rhor)
    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )

    ! Compute the max abs eigenvalue
    dlambda = abs(a*Qroe(2)+b*Qroe(3)) + Qroe(5)





    ! Now compute the temporal lambdas and the differences


    ! Add delta
    Ql(1) = Ql(1) + ddelta

    ! Compute Roe mean values
    rhol = Ql(1)
    rhor = Qr(1)
    ul = Ql(2)/Ql(1)
    ur = Qr(2)/Qr(1)
    vl = Ql(3)/Ql(1)
    vr = Qr(3)/Qr(1)
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    Qroe(1) = sqrt(rhol*rhor)
    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )

    ! Compute the max abs eigenvalue
    dtlambda = abs(a*Qroe(2)+b*Qroe(3)) + Qroe(5)

    ! Save derivative
    Ddlambda(1) = (dtlambda-dlambda)/ddelta

    ! Substract delta
    Ql(1) = Ql(1) - ddelta





    ! Add delta
    Ql(2) = Ql(2) + ddelta

    ! Compute Roe mean values
    rhol = Ql(1)
    rhor = Qr(1)
    ul = Ql(2)/Ql(1)
    ur = Qr(2)/Qr(1)
    vl = Ql(3)/Ql(1)
    vr = Qr(3)/Qr(1)
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    Qroe(1) = sqrt(rhol*rhor)
    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )

    ! Compute the max abs eigenvalue
    dtlambda = abs(a*Qroe(2)+b*Qroe(3)) + Qroe(5)

    ! Save derivative
    Ddlambda(2) = (dtlambda-dlambda)/ddelta

    ! Substract delta
    Ql(2) = Ql(2) - ddelta




    ! Add delta
    Ql(3) = Ql(3) + ddelta

    ! Compute Roe mean values
    rhol = Ql(1)
    rhor = Qr(1)
    ul = Ql(2)/Ql(1)
    ur = Qr(2)/Qr(1)
    vl = Ql(3)/Ql(1)
    vr = Qr(3)/Qr(1)
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    Qroe(1) = sqrt(rhol*rhor)
    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )

    ! Compute the max abs eigenvalue
    dtlambda = abs(a*Qroe(2)+b*Qroe(3)) + Qroe(5)

    ! Save derivative
    Ddlambda(3) = (dtlambda-dlambda)/ddelta

    ! Substract delta
    Ql(3) = Ql(3) - ddelta




    ! Add delta
    Ql(4) = Ql(4) + ddelta

    ! Compute Roe mean values
    rhol = Ql(1)
    rhor = Qr(1)
    ul = Ql(2)/Ql(1)
    ur = Qr(2)/Qr(1)
    vl = Ql(3)/Ql(1)
    vr = Qr(3)/Qr(1)
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    Qroe(1) = sqrt(rhol*rhor)
    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )

    ! Compute the max abs eigenvalue
    dtlambda = abs(a*Qroe(2)+b*Qroe(3)) + Qroe(5)

    ! Save derivative
    Ddlambda(4) = (dtlambda-dlambda)/ddelta

    ! Substract delta
    Ql(4) = Ql(4) - ddelta



    ! Add delta
    Qr(1) = Qr(1) + ddelta

    ! Compute Roe mean values
    rhol = Ql(1)
    rhor = Qr(1)
    ul = Ql(2)/Ql(1)
    ur = Qr(2)/Qr(1)
    vl = Ql(3)/Ql(1)
    vr = Qr(3)/Qr(1)
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    Qroe(1) = sqrt(rhol*rhor)
    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )

    ! Compute the max abs eigenvalue
    dtlambda = abs(a*Qroe(2)+b*Qroe(3)) + Qroe(5)

    ! Save derivative
    Ddlambda(5) = (dtlambda-dlambda)/ddelta

    ! Substract delta
    Qr(1) = Qr(1) - ddelta







    ! Add delta
    Qr(2) = Qr(2) + ddelta

    ! Compute Roe mean values
    rhol = Ql(1)
    rhor = Qr(1)
    ul = Ql(2)/Ql(1)
    ur = Qr(2)/Qr(1)
    vl = Ql(3)/Ql(1)
    vr = Qr(3)/Qr(1)
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    Qroe(1) = sqrt(rhol*rhor)
    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )

    ! Compute the max abs eigenvalue
    dtlambda = abs(a*Qroe(2)+b*Qroe(3)) + Qroe(5)

    ! Save derivative
    Ddlambda(6) = (dtlambda-dlambda)/ddelta

    ! Substract delta
    Qr(2) = Qr(2) - ddelta





    ! Add delta
    Qr(3) = Qr(3) + ddelta

    ! Compute Roe mean values
    rhol = Ql(1)
    rhor = Qr(1)
    ul = Ql(2)/Ql(1)
    ur = Qr(2)/Qr(1)
    vl = Ql(3)/Ql(1)
    vr = Qr(3)/Qr(1)
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    Qroe(1) = sqrt(rhol*rhor)
    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )

    ! Compute the max abs eigenvalue
    dtlambda = abs(a*Qroe(2)+b*Qroe(3)) + Qroe(5)

    ! Save derivative
    Ddlambda(7) = (dtlambda-dlambda)/ddelta

    ! Substract delta
    Qr(3) = Qr(3) - ddelta





    ! Add delta
    Qr(4) = Qr(4) + ddelta

    ! Compute Roe mean values
    rhol = Ql(1)
    rhor = Qr(1)
    ul = Ql(2)/Ql(1)
    ur = Qr(2)/Qr(1)
    vl = Ql(3)/Ql(1)
    vr = Qr(3)/Qr(1)
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    Qroe(1) = sqrt(rhol*rhor)
    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )

    ! Compute the max abs eigenvalue
    dtlambda = abs(a*Qroe(2)+b*Qroe(3)) + Qroe(5)

    ! Save derivative
    Ddlambda(8) = (dtlambda-dlambda)/ddelta

    ! Substract delta
    Qr(4) = Qr(4) - ddelta

  end function 



! This function builds the absolute value of the largest eigenvalue

function Euler_buildlambda(Qll, Qrr, a, b) result(dlambda)

    ! The left and right solution values
    real(DP), dimension(4), intent(IN)		:: Qll, Qrr
    
    ! Lambda and temp Lambda
    real(dp) :: dlambda
    
    ! The components of the normal vector
    real(dp), intent(in) :: a, b

    ! temp variables
    real(DP)		:: rhol, rhor, denom, Hl, Hr, El, Er, pl, pr, velnorm, aux, ul, ur, vl, vr

    ! Roe-values and speed of sound waves
    real(dp), dimension(5) :: Qroe

    ! Gamma
    real(dp) :: gamma = 1.4_dp

    ! Temporary left and right solution values
    real(DP), dimension(4) :: Ql, Qr


    ! Copy solution values to their temporary opposites
    Ql = Qll
    Qr = Qrr

    ! Compute Roe mean values
    rhol = Ql(1)
    rhor = Qr(1)
    ul = Ql(2)/Ql(1)
    ur = Qr(2)/Qr(1)
    vl = Ql(3)/Ql(1)
    vr = Qr(3)/Qr(1)
    El = Ql(4)/rhol
    Er = Qr(4)/rhor
    Qroe(1) = sqrt(rhol*rhor)
    aux  = sqrt(max(rhol/rhor, SYS_EPSREAL_DP))
    Qroe(2) = (aux*ul+ur)/(aux+1.0_DP)
    Qroe(3) = (aux*vl+vr)/(aux+1.0_DP)
    hl   = GAMMA*El-0.5_dp*(gamma-1.0_dp)*(ul*ul+vl*vl)
    hr   = GAMMA*Er-0.5_dp*(gamma-1.0_dp)*(ur*ur+vr*vr)
    Qroe(4) = (aux*hl+hr)/(aux+1.0_DP)
    Qroe(5) = sqrt( max((gamma-1.0_dp)*(Qroe(4) - 0.5_dp*(Qroe(2)*Qroe(2)+Qroe(3)*Qroe(3)) ),0.0_dp) )

    ! Compute the max abs eigenvalue
    dlambda = abs(a*Qroe(2)+b*Qroe(3)) + Qroe(5)

  end function
  
  
  function LLF(uL, uR, nx, ny)
    real(dp) :: uL(4), uR(4) !  Input: conservative variables rho*[1, u, v, E]
    real(dp) :: nx, ny       !  Input: face normal vector, [nx, ny] (Left-to-Right)
    real(dp) :: LLF(4)       ! Output: Roe flux function (upwind)
    
    
    LLF = 0.5_dp*(nx*Euler_buildFlux(uL,1)+ny*Euler_buildFlux(uL,2) + nx*Euler_buildFlux(uR,1)+ny*Euler_buildFlux(uR,2) )&
          -0.5_dp*Euler_buildlambda(uL, uR, nx, ny)*(uR-uL)

    end function
    
  function DLLF(uL, uR, nx, ny, h)
    real(dp) :: uL(4), uR(4)             !  Input: conservative variables rho*[1, u, v, E]
    real(dp) :: nx, ny                   !  Input: face normal vector, [nx, ny]
    real(dp) :: h                        !  Input: infinitesimal constant
    real(dp), dimension(4,8) :: DLLF     ! Output: Approximate differential of Roe flux function
    
    real(dp), dimension(4) :: F, U
    integer :: i
    
    F = LLF(uL, uR, nx, ny)
    
    ! First order approx
    do i = 1,4
      U = uL
      U(i) = U(i) + h
      DLLF(:,i) = (LLF(U, uR, nx, ny) - F)/h
    end do
    
    do i = 1,4
      U = uR
      U(i) = U(i) + h
      DLLF(:,i+4) = (LLF(uL, U, nx, ny) - F)/h
    end do

!    ! Second order approx
!    do i = 1,4
!      U = uL
!      U(i) = U(i) + h
!      DLLF(:,i) = LLF(U, uR, nx, ny)
!      U = uL
!      U(i) = U(i) - h
!      DLLF(:,i) = 0.5_dp*(DLLF(:,i)-LLF(U, uR, nx, ny))/h
!    end do
!    
!    do i = 1,4
!      U = uR
!      U(i) = U(i) + h
!      DLLF(:,i+4) = LLF(uL, U, nx, ny)
!      U = uR
!      U(i) = U(i) - h
!      DLLF(:,i+4) = 0.5_dp*(DLLF(:,i+4)-LLF(uL, U, nx, ny))/h
!    end do
    
  end function DLLF
  
  
  
    function DRotated_RHLL_2nd(uL, uR, nx, ny, h)
    real(dp) :: uL(4), uR(4)             !  Input: conservative variables rho*[1, u, v, E]
    real(dp) :: nx, ny                   !  Input: face normal vector, [nx, ny]
    real(dp) :: h                        !  Input: infinitesimal constant
    real(dp), dimension(4,8) :: DRotated_RHLL_2nd     ! Output: Approximate differential of Roe flux function
    
    real(dp), dimension(4) :: F, U
    integer :: i
    
!    F = Rotated_RHLL(uL, uR, nx, ny)
!    
!    ! First order approx
!    do i = 1,4
!      U = uL
!      U(i) = U(i) + h
!      DRotated_RHLL_2nd(:,i) = (Rotated_RHLL(U, uR, nx, ny) - F)/h
!    end do
!    
!    do i = 1,4
!      U = uR
!      U(i) = U(i) + h
!      DRotated_RHLL_2nd(:,i+4) = (Rotated_RHLL(uL, U, nx, ny) - F)/h
!    end do

    ! Second order approx
    do i = 1,4
      U = uL
      U(i) = U(i) + h
      DRotated_RHLL_2nd(:,i) = Rotated_RHLL(U, uR, nx, ny)
      U = uL
      U(i) = U(i) - h
      DRotated_RHLL_2nd(:,i) = 0.5_dp*(DRotated_RHLL_2nd(:,i)-Rotated_RHLL(U, uR, nx, ny))/h
    end do
    
    do i = 1,4
      U = uR
      U(i) = U(i) + h
      DRotated_RHLL_2nd(:,i+4) = Rotated_RHLL(uL, U, nx, ny)
      U = uR
      U(i) = U(i) - h
      DRotated_RHLL_2nd(:,i+4) = 0.5_dp*(DRotated_RHLL_2nd(:,i+4)-Rotated_RHLL(uL, U, nx, ny))/h
    end do
    
  end function
  
  
  
  
  function Euler_WallFlux(Qi,nx,ny) result(Flux)

    ! The flux vector
    real(DP), dimension(4)	:: Flux

    ! The solution components q1 = rho, q2 = rho u, q3 = rho v, q4 = rho E
    real(DP), dimension(4), intent(IN)		:: Qi
    
    ! Components of the normal vector
    real(dp), intent(in) :: nx, ny

    ! pressure, stagnation enthalpy
    real(DP)                                :: p, H

    ! temporary variable
    real(DP)                                :: rho, u, v, E

    ! Constant Gamma
    real(dp) :: gamma = 1.4_dp

    ! Calculate primitive variables
    rho=Qi(1)
    u=Qi(2)/rho
    v=Qi(3)/rho
    E=Qi(4)/rho

    ! Compute the pressure
    p = (gamma - 1.0_dp)*rho*(E-0.5_dp*(u*u+v*v))

!    ! Compute H, the stagnation enthalpy
!    H = E + p/rho
    
    ! Build the flux vector
    Flux(1) = 0.0_dp
    Flux(2) = nx*p
    Flux(3) = ny*p
    Flux(4) = 0.0_dp

  end function
  
  
  
  function dEuler_WallFlux(uL, nx, ny, h) result(DFlux)
    real(dp) :: uL(4), uR(4)             !  Input: conservative variables rho*[1, u, v, E]
    real(dp) :: nx, ny                   !  Input: face normal vector, [nx, ny]
    real(dp) :: h                        !  Input: infinitesimal constant
    real(dp), dimension(4,8) :: DFlux    ! Output: Approximate differential of Roe flux function
    
    real(dp), dimension(4) :: F, U
    integer :: i
    

    ! Second order approx
    do i = 1,4
      U = uL
      U(i) = U(i) + h
      DFlux(:,i) = Euler_WallFlux(U,nx,ny)
      U = uL
      U(i) = U(i) - h
      DFlux(:,i) = 0.5_dp*(DFlux(:,i)-Euler_WallFlux(U,nx,ny))/h
    end do
    
    DFlux(:,5:8) = 0.0_dp
    
    
  end function
  
  
  
  function calculatePenalty(dpolgrad, dedgeLength) result(dpenalty)
  
  real(dp), intent(in) :: dpolgrad, dedgeLength
  
  real(dp) :: dPenalty
  
  
  dPenalty = (dpolgrad)*(dpolgrad+1.0_dp)/dedgeLength
  
  
  
  end function

end module dg2d_problem
