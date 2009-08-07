module shallowwater2d_routines

    use fsystem
    use storage
    use triangulation
    use spatialdiscretisation
    use linearsystemscalar
    use linearsystemblock
    use boundary

    implicit none
    
    type t_array
		! Pointer to the double-valued matrix or vector data
    	real(DP), dimension(:), pointer :: Da
	end type t_array
	
	
	integer, parameter				:: nvar2d = 3
	

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
	real(DP)		:: whl, whr, denom

	denom = sqrt(Ql(1))+sqrt(Qr(1))
	whl = 1.0_DP/sqrt(Ql(1))
	whr = 1.0_DP/sqrt(Qr(1))

	Qroe(1) = sqrt(Ql(1)*Qr(1))
	Qroe(2) = Qroe(1)*(whl*Ql(2)+whr*Qr(2))/denom
	Qroe(3) = Qroe(1)*(whl*Ql(3)+whr*Qr(3))/denom

	end function
	
	
    
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
	u=Q(2)/Q(1)
	v=Q(3)/Q(1)
	
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
	
	end function
    
    
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
	
	end function
    
    
    
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
	
	end function
	
	
	
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
	
	end function
    
    
    
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
	
	end function
	
	
	
	
	
	
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
	u=Q(2)/Q(1)
	v=Q(3)/Q(1)
	
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
	
	end function



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
	
    if (d==1) then
    ! build eigenvalues in x direction
	    Flux(1) = Q(2)
	    Flux(2) = Q(2)*Q(2)/Q(1)+0.5_DP*g*Q(1)*Q(1)
	    Flux(3) = Q(2)*Q(3)/Q(1)
    else
    ! build eigenvalues in y direction
	    Flux(1) = Q(3)
	    Flux(2) = Q(2)*Q(3)/Q(1)
	    Flux(3) = Q(3)*Q(3)/Q(1)+0.5_DP*g*Q(1)*Q(1)
	end if
	
	end function
	
	
	
	
	
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
	real(DP), dimension(nvar2d,nvar2d)  :: Aij, Bij, Dij, Eye
	real(DP)                            :: scalarDissipation, scalefactor, lambda
	real(DP)                            :: cRoe, uRoe, vRoe
	
	
	! Assemble the preconditioner rmatrixBlockP: P = ML - theta*dt*L
	! As we use the block jacobi method, we only need the main diagonal blocks
	
	    ! unit matrix
	    Eye = 0.0_DP
	    forall (ivar = 1: nvar2d)
            Eye(ivar,ivar) = 1.0_DP
	    end forall
	
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
			
			! Compute the jacobi matrices for x- and y- direction
			JacobixRoeij = buildJacobi(Qroeij,1,gravconst)
			JacobiyRoeij = buildJacobi(Qroeij,2,gravconst)
			
			! Compute the coeffitients for the jacobi matrices
			JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
			JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
			JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
			JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP
			
			! Now we can compute Aij and Bij
			Aij = JcoeffxA*JacobixRoeij + JcoeffyA*JacobiyRoeij
			Bij = JcoeffxB*JacobixRoeij + JcoeffyB*JacobiyRoeij
			
			! Here we use scalar dissipation
			! compute the maximum of the eigenvalues of Aij
			scalefactor = sqrt(JcoeffxA**2.0_DP+JcoeffyA**2.0_DP)
			cRoe = sqrt(gravconst*Qroeij(1))    ! c_Roe=sqrt(g*h)
			cRoe = sqrt(0.5*gravconst*(Qi(1)+Qj(1)))
			uRoe = Qroeij(2)/Qroeij(1)
			vRoe = Qroeij(3)/Qroeij(1)
			lambda = max(abs(uRoe-cRoe),abs(uRoe+cRoe),abs(vRoe-cRoe),abs(vRoe+cRoe),abs(uRoe),abs(vRoe))
			scalarDissipation = scalefactor*lambda
			
			scalarDissipation = &
			abs(JcoeffxA*maxval(abs(buildeigenvalues(Qroeij,1,gravconst))))+ &
			abs(JcoeffyA*maxval(abs(buildeigenvalues(Qroeij,2,gravconst))))
			
			! compute Dij
			Dij = scalarDissipation*Eye
			
			! Now add the entries of Aij and Dij to their corresponding entries in P
			! P = M^L - theta*dt*L
			do l = 1, nvar2d
				rarrayP(l)%Da(ii) = rarrayP(l)%Da(ii) - dt*theta*(+Aij(l,l)-Dij(l,l))
				rarrayP(l)%Da(ij) = rarrayP(l)%Da(ij) - dt*theta*(-Aij(l,l)+Dij(l,l))
				rarrayP(l)%Da(ji) = rarrayP(l)%Da(ji) - dt*theta*(+Aij(l,l)+Dij(l,l))
				rarrayP(l)%Da(jj) = rarrayP(l)%Da(jj) - dt*theta*(-Aij(l,l)-Dij(l,l))
			end do
			                        
		end do
	
	end subroutine
	
	
	
	
	
	
	! This routine builds the RHS for the shallow water system
	subroutine BuildShallowWaterRHS (&
	                rarrayRhs, rarraySol, rrhsBlock, rsolBlock, &
	                rmatrixML, p_CXdata, p_CYdata, p_MLdata, &
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
	real(DP), dimension(:), pointer, intent(IN)     :: p_CXdata, p_CYdata, p_MLdata
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
	real(DP), dimension(nvar2d,nvar2d)  :: Aij, Bij, Dij, Eye
	real(DP)                            :: scalarDissipation, scalefactor, lambda
	real(DP), dimension(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj
	real(DP)                            :: cRoe, uRoe, vRoe
	    ! for TVD
	real(DP), dimension(nvar2d,nvar2d)  :: invRij, Rij
	real(DP), dimension(nvar2d)         :: deltaWij, eigenvalues, deltaFij
	real(DP)                            :: scalefactor1, scalefactor2, deltaak, deltaWijHat
	real(DP)                            :: deltaaplus, deltaaminus, deltabplus, deltabminus
	integer                             :: upwindnode
	
	
	    ! unit matrix
	    Eye = 0.0_DP
	    forall (ivar = 1: nvar2d)
            Eye(ivar,ivar) = 1.0_DP
	    end forall
	    
	    
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
			
			! Compute the jacobi matrices for x- and y- direction
			JacobixRoeij = buildJacobi(Qroeij,1,gravconst)
			JacobiyRoeij = buildJacobi(Qroeij,2,gravconst)
			
			! Compute the coeffitients for the jacobi matrices
			JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
			JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
			JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
			JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP
			
			! Now we can compute Aij and Bij
			Aij = JcoeffxA*JacobixRoeij + JcoeffyA*JacobiyRoeij
			Bij = JcoeffxB*JacobixRoeij + JcoeffyB*JacobiyRoeij
			
			! deltaK = Aij*deltaQij
			deltaKi = matmul(Aij+Bij,deltaQij)
			deltaKj = matmul(Aij-Bij,deltaQij)
			
			! Calculate this alternatively by calculating
			! deltaKi = c_{ij}*(F(Q_i)-F(Q_j))
			! deltaKj = c_{ji}*(F(Q_j)-F(Q_i))
			deltaKi = p_CXdata(ij)*buildFlux(Qi,1,gravconst)+p_CYdata(ij)*buildFlux(Qi,2,gravconst) &
			         -p_CXdata(ij)*buildFlux(Qj,1,gravconst)-p_CYdata(ij)*buildFlux(Qj,2,gravconst)
			deltaKj = p_CXdata(ji)*buildFlux(Qj,1,gravconst)+p_CYdata(ji)*buildFlux(Qj,2,gravconst) &
		    	     -p_CXdata(ji)*buildFlux(Qi,1,gravconst)-p_CYdata(ji)*buildFlux(Qi,2,gravconst)
			
			! Now choose the artificial diffusion method
			if ((Method == 0).or.(Method == 3)) then
			    ! Here we do not add artificial diffusion - so we have the high order method
			    Dij = 0
			else if ((Method == 1).or.(Method == 4).or.(Method == 6)) then
			    ! Here we use scalar dissipation
			    ! compute the maximum of the eigenvalues of Aij
			    scalefactor = sqrt(JcoeffxA**2.0_DP+JcoeffyA**2.0_DP)
			    cRoe = sqrt(gravconst*Qroeij(1))    ! c_Roe=sqrt(g*h)
			    cRoe = sqrt(0.5*gravconst*(Qi(1)+Qj(1)))
			    uRoe = Qroeij(2)/Qroeij(1)
			    vRoe = Qroeij(3)/Qroeij(1)
			    lambda = max(abs(uRoe-cRoe),abs(uRoe+cRoe),abs(vRoe-cRoe),abs(vRoe+cRoe),abs(uRoe),abs(vRoe))
			    scalarDissipation = scalefactor*lambda
			    
			    scalarDissipation = &
			abs(JcoeffxA*maxval(abs(buildeigenvalues(Qroeij,1,gravconst))))+ &
			abs(JcoeffyA*maxval(abs(buildeigenvalues(Qroeij,2,gravconst))))
			    
			    ! compute Dij
			    Dij = scalarDissipation*Eye
			else if ((Method == 2).or.(Method == 5).or.(Method == 7)) then
			    ! Here we use tensorial dissipation
			    Dij = abs(JcoeffxA)*matmul(buildTrafo(Qroeij,1,gravconst),&
			                    matmul(buildaLambda(Qroeij,1,gravconst),&
			                           buildinvTrafo(Qroeij,1,gravconst))) +&
			          abs(JcoeffyA)*matmul(buildTrafo(Qroeij,2,gravconst),&
			                    matmul(buildaLambda(Qroeij,2,gravconst),&
			                           buildinvTrafo(Qroeij,2,gravconst)))
			end if
			
			! deltaD
			deltaDi = matmul(Dij-Bij,-deltaQij)
			deltaDj = -deltaDi
			
			! add deltaK and deltaD to rhs
			! rhs = rhs + (1-theta)*dt*K*u + (1-theta)*dt*D*u
			do l = 1, nvar2d
				rarrayRhs(l)%Da(i) = rarrayRhs(l)%Da(i) + (1-theta)*dt*deltaKi(l) &
			                        	+ (1-theta)*dt*deltaDi(l) 
				rarrayRhs(l)%Da(j) = rarrayRhs(l)%Da(j) + (1-theta)*dt*deltaKj(l) &
			                        	+ (1-theta)*dt*deltaDj(l)
			end do
			
		end do
		
		
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
			! rhs(i) = rhs(i) + (1-theta)*dt*deltaF(i)
			! rhs(j) = rhs(j) - (1-theta)*dt*deltaF(j)
			do l = 1, nvar2d
				rarrayRhs(l)%Da(i) = rarrayRhs(l)%Da(i) + (1-theta)*dt*deltaFij(l)
				rarrayRhs(l)%Da(j) = rarrayRhs(l)%Da(j) - (1-theta)*dt*deltaFij(l)
			end do
			
	    end do
        
        end do
		end if
		
	
	end subroutine
	
	
	
	
	
	
	
	
	
	
	
	! This routine builds the defect for the shallow water system
	subroutine BuildShallowWaterDefect (&
	                rdefBlock, rstempBlock, rrhsBlock, rsolBlock, &
	                rarrayRhs, rarraySol, rarrayRstemp, &
	                rmatrixML, p_CXdata, p_CYdata, p_MLdata, &
	                h_fld1, p_fld1, p_fld2, &
	                p_Kdiagonal, p_Kedge, NEQ, nedge, &
	                theta, dt, gravconst, Method, limiter)
	
	! parameter values
    type(t_array), dimension(nvar2d), intent(INOUT) :: rarrayRstemp
	type(t_array), dimension(nvar2d), intent(IN)    :: rarraySol, rarrayRhs
	type(t_vectorBlock), intent(INOUT)              :: rdefBlock, rstempBlock
	type(t_vectorBlock), intent(IN)                 :: rsolBlock, rrhsBlock
	type(t_matrixScalar), intent(IN)                :: rmatrixML
	real(DP), dimension(:), pointer, intent(IN)     :: p_CXdata, p_CYdata, p_MLdata
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
	real(DP), dimension(nvar2d,nvar2d)  :: Aij, Bij, Dij, Eye
	real(DP)                            :: scalarDissipation, scalefactor, lambda
	real(DP), dimension(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj
	real(DP)                            :: cRoe, uRoe, vRoe
	    ! for TVD
	real(DP), dimension(nvar2d,nvar2d)  :: invRij, Rij
	real(DP), dimension(nvar2d)         :: deltaWij, eigenvalues, deltaFij
	real(DP)                            :: scalefactor1, scalefactor2, deltaak, deltaWijHat
	real(DP)                            :: deltaaplus, deltaaminus, deltabplus, deltabminus
	integer                             :: upwindnode
	
	
	    ! unit matrix
	    Eye = 0.0_DP
	    forall (ivar = 1: nvar2d)
            Eye(ivar,ivar) = 1.0_DP
	    end forall
	    
	    
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
			
    			! Compute the Roe-meanvalue for this edge
	    		Qroeij = calculateQroe(Qi, Qj)
			
    			! Compute the jacobi matrices for x- and y- direction
	    		JacobixRoeij = buildJacobi(Qroeij,1,gravconst)
		    	JacobiyRoeij = buildJacobi(Qroeij,2,gravconst)
			
			    ! Compute the coeffitients for the jacobi matrices
			    JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
			    JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
			    JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
			    JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP
			
			    ! Now we can compute Aij and Bij
			    Aij = JcoeffxA*JacobixRoeij + JcoeffyA*JacobiyRoeij
			    Bij = JcoeffxB*JacobixRoeij + JcoeffyB*JacobiyRoeij
			
			    ! deltaK
			    deltaKi = matmul(Aij+Bij,deltaQij)
			    deltaKj = matmul(Aij-Bij,deltaQij)

				! Calculate this alternatively by calculating
				! deltaKi = c_{ij}*(F(Q_i)-F(Q_j))
				! deltaKj = c_{ji}*(F(Q_j)-F(Q_i))
				deltaKi = p_CXdata(ij)*buildFlux(Qi,1,gravconst)+p_CYdata(ij)*buildFlux(Qi,2,gravconst) &
				         -p_CXdata(ij)*buildFlux(Qj,1,gravconst)-p_CYdata(ij)*buildFlux(Qj,2,gravconst)
				deltaKj = p_CXdata(ji)*buildFlux(Qj,1,gravconst)+p_CYdata(ji)*buildFlux(Qj,2,gravconst) &
			    	     -p_CXdata(ji)*buildFlux(Qi,1,gravconst)-p_CYdata(ji)*buildFlux(Qi,2,gravconst)
			    
			    ! Now choose the artificial diffusion method
			    if ((Method == 0).or.(Method == 3)) then
			        ! Here we do not add artificial diffusion - so we have the high order method
			        Dij = 0
			    else if ((Method == 1).or.(Method == 4).or.(Method == 6)) then
			        ! Here we use scalar dissipation
			        ! compute the maximum of the eigenvalues of Aij
			        scalefactor = sqrt(JcoeffxA**2.0_DP+JcoeffyA**2.0_DP)
			        cRoe = sqrt(gravconst*Qroeij(1))    ! c_Roe=sqrt(g*h)
			        cRoe = sqrt(0.5*gravconst*(Qi(1)+Qj(1)))
			        uRoe = Qroeij(2)/Qroeij(1)
			        vRoe = Qroeij(3)/Qroeij(1)
			        lambda = max(abs(uRoe-cRoe),abs(uRoe+cRoe),abs(vRoe-cRoe),abs(vRoe+cRoe),abs(uRoe),abs(vRoe))
			        scalarDissipation = scalefactor*lambda
			        
			        scalarDissipation = &
			abs(JcoeffxA*maxval(abs(buildeigenvalues(Qroeij,1,gravconst))))+ &
			abs(JcoeffyA*maxval(abs(buildeigenvalues(Qroeij,2,gravconst))))
			        
			        ! compute Dij
			        Dij = scalarDissipation*Eye
			    else if ((Method == 2).or.(Method == 5).or.(Method == 7)) then
			        ! Here we use tensorial dissipation
			        Dij = abs(JcoeffxA)*matmul(buildTrafo(Qroeij,1,gravconst),&
			                        matmul(buildaLambda(Qroeij,1,gravconst),&
			                               buildinvTrafo(Qroeij,1,gravconst))) +&
			              abs(JcoeffyA)*matmul(buildTrafo(Qroeij,2,gravconst),&
			                        matmul(buildaLambda(Qroeij,2,gravconst),&
			                               buildinvTrafo(Qroeij,2,gravconst)))
			    end if
			    
! 			    ! deltaD
			    deltaDi = matmul(Dij-Bij,-deltaQij)
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
    
	end subroutine	
	
	
	
	
	
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
	real(DP), dimension(nvar2d,nvar2d)  :: Aij, Bij, Dij, Eye
	real(DP)                            :: scalarDissipation, scalefactor, lambda
	real(DP), dimension(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj
	real(DP)                            :: cRoe, uRoe, vRoe
	real(DP), dimension(nvar2d)         :: Udoti, Udotj
	real(DP), dimension(nvar2d)         :: deltaWij, deltaGij, deltaFij
	real(DP), dimension(nvar2d,nvar2d)  :: Rij, invRij
	
	
	! unit matrix
	Eye = 0.0_DP
	forall (ivar = 1: nvar2d)
        Eye(ivar,ivar) = 1.0_DP
	end forall

    
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
			    scalarDissipation = maxval(abs(buildEigenvalues(Qroeij,d,gravconst)))
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
			    scalarDissipation = maxval(abs(buildEigenvalues(Qroeij,d,gravconst)))
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
       
	end subroutine
	
	
	
	
	
	
	
	
	
	
	! This subroutine takes care of boundary conditions
	! for a reflacting boundary
	subroutine ImplementShallowWaterBCs (&
	                rboundary, rtriangulation, &
	                rarrayP, rarraySol, rarrayDef, &
	                p_Kdiagonal, p_Kld, &
	                gravconst, boundarycorner)
	                
	! PARAMETER VALUES
	! the boundary
	type(t_boundary), intent(IN) :: rboundary
	! the triangulation
	type(t_triangulation), intent(IN) :: rtriangulation
	! pointers to datas of preconditioner, solution and defect
	type(t_array), dimension(nvar2d), intent(INOUT) :: rarrayP, rarraySol, rarrayDef
	! pointer to index of diagonal positions of matrices
	integer, dimension(:), pointer, intent(IN) :: p_Kdiagonal, p_Kld
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
	
	
	! Initialise some values
	
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
	        ustar = rarraySol(2)%Da(i)/hstar
	        vstar = rarraySol(3)%Da(i)/hstar
	        
	        ! update defect
	        rarrayDef(1)%Da(i) = 0.0_dp
	        rarrayDef(2)%Da(i) = 0.0_dp
	        rarrayDef(3)%Da(i) = 0.0_dp
	        
	        
	        ! solve the riemannian problem
	        
	        ! projection of velocity on outward normal vector and tangential vector
	        vnproj = ustar * normalV(1) + vstar * normalV(2)
	        vtproj = ustar * tangentialV(1) + vstar * tangentialV(2)
	        
	        if (4.0_dp*(sqrt(hstar*gravconst))<-2.0_dp*vnproj) then
	            write (*,*) 'ERROR: Boundary Riemann Problem - Physical drybed generation! Wetbet Riemann Solver not applicable!'
	            write (*,*) hstar,vnproj
	        end if
	        
	        ! Initialisiere hstarstar
	        ! entweder mit der predicted solution
	         hstarstar = hstar
	        ! oder mit dem two-rarefaction ansatz (der explizit ausgerechnet werden kann)
	        !hstarstar = ((SQRT(gravconst*hstar)+0.5*vnproj)**2.0_dp)/gravconst
	        	        
	        iteh: do ite = 1, itemax
	        ! get h** as solution of the riemannian problem
	        oldh = hstarstar
	        if (hstarstar .le. hstar) then
	            ! h** can be computed explicitely
	            hstarstar = ((sqrt(gravconst*hstar)+0.5*vnproj)**2.0_dp)/gravconst
	        else
	            
	            hm = hstarstar
	            
	            ! compute h** by fixpointiteration
	            hm = hstar + vnproj/sqrt(0.5*gravconst*(hm+hstar)/(hm*hstar))
	            hm = hstar + vnproj/sqrt(0.5*gravconst*(hm+hstar)/(hm*hstar))
	            hm = hstar + vnproj/sqrt(0.5*gravconst*(hm+hstar)/(hm*hstar))
	            
	            ! or by newton method
	            ! hm = hm - ...
	            
	            hstarstar = hm
	        end if
	        ! Test if the algorithm converged
	        if (abs(hstarstar - oldh)< 1e-8) then
	            exit iteh
	        end if
	        end do iteh
	        
	        ! Test for convergence and give out warnings
	        if (ite==50) then
	            write (*,*) 'ERROR! Boundary condition riemann problem did not converge'
	        end if
	        if (hstarstar<0.05) then
	            write (*,*) 'ERROR! Boundary condition riemann problem h very small or smaller than zero'
	            hstarstar=hstar
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
	
    end subroutine





	! linearised FCT Method for Shallow Water
	! adding limited antidiffusion to the low order predictor
	! using conservative variables variables and a syncronized limiting factor
	subroutine linFctShallowWaterAddLimitedAntidiffusion_syncronized(&
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
	real(DP), dimension(nvar2d,nvar2d)  :: Aij, Bij, Dij, Eye, Tij, invTij
	real(DP)                            :: scalarDissipation, scalefactor, lambda
	real(DP), dimension(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj
	real(DP)                            :: cRoe, uRoe, vRoe
	real(DP), dimension(nvar2d)         :: Udoti, Udotj
	real(DP), dimension(nvar2d)         :: deltaWij, deltaGij, deltaFij
	real(DP), dimension(nvar2d,nvar2d)  :: Rij, invRij
	real(DP)                            :: alphaij
	
	
	! unit matrix
	Eye = 0.0_DP
	forall (ivar = 1: nvar2d)
        Eye(ivar,ivar) = 1.0_DP
	end forall

    

    
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

			! Compute the coeffitients for the jacobi matrices
			JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
			JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
			JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
			JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP
			    
			! compute Dij
			if (Method == 6) then
			    ! Build Dij using scalar dissipation
	    	scalarDissipation = &
			abs(JcoeffxA*maxval(abs(buildeigenvalues(Qroeij,1,gravconst))))+ &
			abs(JcoeffyA*maxval(abs(buildeigenvalues(Qroeij,2,gravconst))))
			
			! compute Dij
			Dij = scalarDissipation*Eye
			else if (Method == 7) then
			    ! Build Dij using tensorial dissipation
			    Dij = abs(JcoeffxA)*matmul(buildTrafo(Qroeij,1,gravconst),&
			                    matmul(buildaLambda(Qroeij,1,gravconst),&
			                           buildinvTrafo(Qroeij,1,gravconst))) +&
			          abs(JcoeffyA)*matmul(buildTrafo(Qroeij,2,gravconst),&
			                    matmul(buildaLambda(Qroeij,2,gravconst),&
			                           buildinvTrafo(Qroeij,2,gravconst)))
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
			
            ! Trafomatrix Tij
            !invRij = buildInvTrafo(Qroeij,d,gravconst)
			Tij = Eye
            
            ! compute Tij*(Qj - Qi), the transformed solution differences
            deltaWij = -matmul(Tij,deltaQij)
            
			! Compute the coeffitients for the jacobi matrices
			JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
			JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
			JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
			JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP
			    
			! compute Dij
			if (Method == 6) then
			    ! Build Dij using scalar dissipation
	    		scalarDissipation = &
				abs(JcoeffxA*maxval(abs(buildeigenvalues(Qroeij,1,gravconst))))+ &
				abs(JcoeffyA*maxval(abs(buildeigenvalues(Qroeij,2,gravconst))))
			
				Dij = scalarDissipation*Eye
			else if (Method == 7) then
			    ! Build Dij using tensorial dissipation
			    Dij = abs(JcoeffxA)*matmul(buildTrafo(Qroeij,1,gravconst),&
			                    matmul(buildaLambda(Qroeij,1,gravconst),&
			                           buildinvTrafo(Qroeij,1,gravconst))) +&
			          abs(JcoeffyA)*matmul(buildTrafo(Qroeij,2,gravconst),&
			                    matmul(buildaLambda(Qroeij,2,gravconst),&
			                           buildinvTrafo(Qroeij,2,gravconst)))
			end if
            
            ! get Udot at node i and j
            do ivar = 1, nvar2d
                Udoti(ivar) = rarraySolDot(ivar)%Da(i)
                Udotj(ivar) = rarraySolDot(ivar)%Da(j)
            end do
            
            ! compute deltaFij (73)/(59) linearised fct
            deltaFij = p_MCdata(ij)*(Udoti-Udotj)+matmul(Dij,deltaQij)
            
            ! compute deltaGij, the transformed antidiffusive fluxes
            deltaGij = matmul(Tij,deltaFij)

!!!!!!!!!!!!!!!!!!! Check methods
!!!!!!!!!!!!!!!!!!! < ---> >, siehe AFCII(70)
!             ! Now apply prelimiting
!             if (prelimiting == 1) then
!                 if (Method==4) then                 
!                     ! For scalar dissipation apply
!                     ! MinMod prelimiting
!                     do ivar = 1, nvar2d
!                     if (deltaGij(ivar)*deltaWij(ivar)<0) then
!                         deltaGij(ivar)=0
!                     else
!                         if (abs(deltaGij(ivar))>abs(abs(scalefactor)*scalarDissipation*deltaWij(ivar))) then
!                                 deltaGij(ivar) = abs(scalefactor)*scalarDissipation*deltaWij(ivar)
!                         end if
!                     end if
!                     end do
!                 else if (Method==5) then
!                 ! For tensorial dissipation apply
!                 ! Simple prelimiting
!                     do ivar = 1, nvar2d
!                         if (deltaGij(ivar)*deltaWij(ivar)<0) then
!                             deltaGij(ivar) = 0
!                         end if
!                     end do
!                 end if
!             end if


            
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
			end do

			! First calculate the limiting factor alphaij
			! For this we have two options
			! 1.: Take only the limiting factor of the variable h=Q(1)
			if (deltaGij(1)>0.0_DP) then
		        alphaij =  min(p_fld1(5,i),p_fld1(6,j))
		    else
		        alphaij = min(p_fld1(5,j),p_fld1(6,i))
		    end if

! 			! 2.: Take the smallest of all the limiting factors
! 			alphaij = 1
! 			do ivar = 1, nvar2d
! 		        if (deltaGij(ivar)>0.0_DP) then
! 		        	alphaij =  min(alphaij,p_fld1(5+6*(ivar-1),i),p_fld1(6+6*(ivar-1),j))
! 		    	else
! 		        	alphaij = min(alphaij,p_fld1(5+6*(ivar-1),j),p_fld1(6+6*(ivar-1),i))
! 		    	end if
! 			end do

			! Limit the antidiffusive fluxes Gij = Gij * alphaij
			deltaGij = deltaGij * alphaij
			
			! get Roe-values at ij
			!do ivar = 1, nvar2d
	    	!    Qroeij(ivar) = p_fld2(1+2*(ivar-1),iedge)
		    !end do
		    			
            ! Rij
            !Rij = buildTrafo(Qroeij,d,gravconst)
			
			! Get matrix to trafo back
			invTij = Eye
			
			! Finally we can transform back to deltaFij
			deltaFij = matmul(invTij,deltaGij)
			
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
			    
			    rarraySol(ivar)%Da(i) = rarraySol(ivar)%Da(i) + deltaFij(ivar)   *dt/p_MLdata(i)
			    rarraySol(ivar)%Da(j) = rarraySol(ivar)%Da(j) - deltaFij(ivar)   *dt/p_MLdata(j)
			end do
			
	    end do
            
    

       
	end subroutine



    
end module
