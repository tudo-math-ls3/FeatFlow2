MODULE shallowwater2d_routines

    USE fsystem
    USE linearsystemblock
    USE boundary

    IMPLICIT NONE
    
    TYPE t_array
		! Pointer to the double-valued matrix or vector data
    	REAL(DP), DIMENSION(:), POINTER :: Da
	END TYPE t_array
	
	
	INTEGER, PARAMETER				:: nvar2d = 3
	

CONTAINS

    
! Limiter Function for TVD
	REAL(DP) FUNCTION limiterfunc(z,n,limiter)
		IMPLICIT NONE
		REAL(DP), INTENT(IN)		:: z, n				! z = nominator, n = denominator of slope ratio
		INTEGER, INTENT(IN)			:: limiter  		! choice of limiter (1 = Minmod, 4 = Superbee)
		REAL(DP)					:: r				! slope ratio

			limiterfunc = 0.0_DP

		IF (ABS(n)<1e-8) THEN
			limiterfunc = 0.0_DP
		ELSE
			r=z/n
            SELECT CASE(limiter)
                CASE(1) ! MinMod
                    limiterfunc = max(0.0_DP,min(1.0_DP,r))
                CASE(2) ! Van Leer
                    limiterfunc = (r+abs(r))/(1.0_DP+abs(r))
                CASE(3) ! MC
                    limiterfunc = max(0.0_DP, min(2.0_DP,(1.0_DP+r)/2.0_DP,2.0_DP*r))
                CASE(4) ! Superbee
                    limiterfunc = max(0.0_DP,min(1.0_DP,2.0_DP*r),min(2.0_DP,r))
            END SELECT
         END IF
    END FUNCTION limiterfunc
    
    
! 2 parameter limiter function for TVD
	REAL(DP) FUNCTION limiterfunc2(z,n,limiter)
    	IMPLICIT NONE
		REAL(DP), INTENT(IN)		:: z, n				! z = nominator, n = denominator of slope ratio
		INTEGER, INTENT(IN)			:: limiter  		! choice of limiter (1 = Minmod, 4 = Superbee)
		REAL(DP)					:: h1				! temporary variable
		
		h1 = (SIGN(1.0_DP,z)+SIGN(1.0_DP,n))/2.0d0
		
		SELECT CASE(limiter)
			CASE(1) ! MinMod
					limiterfunc2 = h1*MIN(ABS(z),ABS(n))
			CASE(2) ! Van Leer
					IF (ABS(z)+ABS(n)>0.001_DP) THEN	
						limiterfunc2 = h1*2.0_DP*ABS(z*n)/(ABS(n)+ABS(z))
					ELSE
						limiterfunc2 = 0.0_DP
					END IF
			CASE(3) ! MC
					limiterfunc2=h1*MIN(2.0_DP*ABS(z),(ABS(z)+ABS(n))/2.0_DP,2.0_DP*ABS(n))
			CASE(4) ! Superbee
					limiterfunc2=h1*MAX(MIN(2.0_DP*ABS(z),ABS(n)),MIN(ABS(z),2.0_DP*ABS(n)))
		END SELECT
    END FUNCTION limiterfunc2
    
    

	! This function returns the Roe mean values
	FUNCTION calculateQroe(Ql, Qr) RESULT(Qroe)

	! The left and right Q values
	! The solution components q1 = h, q2 = uh, q3 = vh
	REAL(DP), DIMENSION(3), INTENT(IN)		:: Ql, Qr

	! The computed Roe values
	REAL(DP), DIMENSION(3)					:: Qroe
	
	! temp variables
	REAL(DP)		:: whl, whr, denom

	denom = sqrt(Ql(1))+sqrt(Qr(1))
	whl = 1.0_DP/sqrt(Ql(1))
	whr = 1.0_DP/sqrt(Qr(1))

	Qroe(1) = sqrt(Ql(1)*Qr(1))
	Qroe(2) = Qroe(1)*(whl*Ql(2)+whr*Qr(2))/denom
	Qroe(3) = Qroe(1)*(whl*Ql(3)+whr*Qr(3))/denom

	END FUNCTION
	
	
    
    ! This routine builds the jacobi matrix in direction d
    ! d=1: x-direction, d=2: y-direction
	FUNCTION buildJacobi(Q,d,g) RESULT(J)
	
	! The jacobi matrix in direction d
	REAL(DP), DIMENSION(3,3)	:: J
	
	! The solution components q1 = h, q2 = uh, q3 = vh
	REAL(DP), DIMENSION(3), INTENT(IN)		:: q

	! The gravitational konstant
	REAL(DP), INTENT(IN)		        	:: g
	
	INTEGER, INTENT(IN)                     :: d
	
	! primitive variables
	REAL(DP)                                :: h, u, v, c
	
	! Calculate primitive variables
	h=Q(1)
	u=Q(2)/Q(1)
	v=Q(3)/Q(1)
	
    IF (d==1) THEN
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
    ELSE
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
	END IF
	
	END FUNCTION
    
    
    ! This routine builds the trafo matrix Rij in direction d
    ! d=1: x-direction, d=2: y-direction
	FUNCTION buildTrafo(Q,d,g) RESULT(Rij)
	
	! The jacobi matrix in direction d
	REAL(DP), DIMENSION(3,3)	:: Rij
	
	! The solution components q1 = h, q2 = uh, q3 = vh
	REAL(DP), DIMENSION(3), INTENT(IN)		:: Q

	! The gravitational konstant
	REAL(DP), INTENT(IN)		        	:: g
	
	! the direction: d=1: x-direction, d=2: y-direction
	INTEGER, INTENT(IN)                     :: d
	
	! speed of gravitational waves
	REAL(DP)                                :: c
	
	! temporary variable
	REAL(DP)                                :: coeff, h, u, v
	
	! Calculate primitive variables
	h=Q(1)
	u=Q(2)/Q(1)
	v=Q(3)/Q(1)
	
	! compute c = sqrt(g*h)
	c = SQRT(g*h)
	
	

    IF (d==1) THEN
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
    ELSE
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
	END IF
	
	END FUNCTION
    
    
    
    ! This routine builds the inv trafo matrix Rij^{-1} in direction d
    ! d=1: x-direction, d=2: y-direction
	FUNCTION buildInvTrafo(Q,d,g) RESULT(invRij)
	
	! The jacobi matrix in direction d
	REAL(DP), DIMENSION(3,3)	:: invRij
	
	! The solution components q1 = h, q2 = uh, q3 = vh
	REAL(DP), DIMENSION(3), INTENT(IN)		:: Q

	! The gravitational konstant
	REAL(DP), INTENT(IN)		        	:: g
	
	! the direction: d=1: x-direction, d=2: y-direction
	INTEGER, INTENT(IN)                     :: d
	
	! speed of gravitational waves
	REAL(DP)                                :: c
	
	! temporary variable
	REAL(DP)                                :: coeff, h, u, v
	
	! Calculate primitive variables
	h=Q(1)
	u=Q(2)/Q(1)
	v=Q(3)/Q(1)
	
	! compute c = sqrt(g*h)
	c = SQRT(g*h)
	
	! temporary variable
	coeff = 1.0_DP/(2.0_DP*c)
	

    IF (d==1) THEN
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
    ELSE
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
	END IF
	
	END FUNCTION
	
	
	
	! This routine builds the diagonalised jacobi matrix Lambda in direction d
    ! d=1: x-direction, d=2: y-direction
	FUNCTION buildLambda(Q,d,g) RESULT(Lambda)
	
	! The jacobi matrix in direction d
	REAL(DP), DIMENSION(3,3)	:: Lambda
	
	! The solution components q1 = h, q2 = uh, q3 = vh
	REAL(DP), DIMENSION(3), INTENT(IN)		:: Q

	! The gravitational konstant
	REAL(DP), INTENT(IN)		        	:: g
	
	! the direction: d=1: x-direction, d=2: y-direction
	INTEGER, INTENT(IN)                     :: d
	
	! speed of gravitational waves
	REAL(DP)                                :: c
	
	! temporary variable
	REAL(DP)                                :: coeff, h, u, v
	
	! Calculate primitive variables
	h=Q(1)
	u=Q(2)/Q(1)
	v=Q(3)/Q(1)
	
	! compute c = sqrt(g*h)
	c = SQRT(g*h)
	

    IF (d==1) THEN
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
    ELSE
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
	END IF
	
	END FUNCTION
    
    
    
	! This routine builds the absolute value of the diagonalised jacobi
	! matrix aLambda in direction d
    ! d=1: x-direction, d=2: y-direction
	FUNCTION buildaLambda(Q,d,g) RESULT(aLambda)
	
	! The jacobi matrix in direction d
	REAL(DP), DIMENSION(3,3)	:: aLambda
	
	! The solution components q1 = h, q2 = uh, q3 = vh
	REAL(DP), DIMENSION(3), INTENT(IN)		:: Q

	! The gravitational konstant
	REAL(DP), INTENT(IN)		        	:: g
	
	! the direction: d=1: x-direction, d=2: y-direction
	INTEGER, INTENT(IN)                     :: d
	
	! speed of gravitational waves
	REAL(DP)                                :: c
	
	! temporary variable
	REAL(DP)                                :: coeff, h, u, v
	
	! Calculate primitive variables
	h=Q(1)
	u=Q(2)/Q(1)
	v=Q(3)/Q(1)
	
	! compute c = sqrt(g*h)
	c = SQRT(g*h)
	

    IF (d==1) THEN
    ! build aLambda in x direction
	    aLambda(1,1) = ABS(u-c)
	    aLambda(2,1) = 0.0_DP
	    aLambda(3,1) = 0.0_DP
	    aLambda(1,2) = 0.0_DP
	    aLambda(2,2) = ABS(u)
	    aLambda(3,2) = 0.0_DP
	    aLambda(1,3) = 0.0_DP
	    aLambda(2,3) = 0.0_DP
	    aLambda(3,3) = ABS(u+c)
    ELSE
    ! build aLambda in y direction
	    aLambda(1,1) = ABS(v-c)
	    aLambda(2,1) = 0.0_DP
	    aLambda(3,1) = 0.0_DP
	    aLambda(1,2) = 0.0_DP
	    aLambda(2,2) = ABS(v)
	    aLambda(3,2) = 0.0_DP
	    aLambda(1,3) = 0.0_DP
	    aLambda(2,3) = 0.0_DP
	    aLambda(3,3) = ABS(v+c)
	END IF
	
	END FUNCTION
	
	
	
	
	
	
    ! This routine returns the eigenvalues of the jacobi matrix in direction d
    ! d=1: x-direction, d=2: y-direction
	FUNCTION buildEigenvalues(Q,d,g) RESULT(Eigenvalues)
	
	! The jacobi matrix in direction d
	REAL(DP), DIMENSION(3)	:: Eigenvalues
	
	! The solution components q1 = h, q2 = uh, q3 = vh
	REAL(DP), DIMENSION(3), INTENT(IN)		:: Q

	! The gravitational konstant
	REAL(DP), INTENT(IN)		        	:: g
	
	! the direction: d=1: x-direction, d=2: y-direction
	INTEGER, INTENT(IN)                     :: d
	
	! speed of gravitational waves
	REAL(DP)                                :: c
	
	! temporary variable
	REAL(DP)                                :: coeff, h, u, v
	
	! Calculate primitive variables
	h=Q(1)
	u=Q(2)/Q(1)
	v=Q(3)/Q(1)
	
	! compute c = sqrt(g*h)
	c = SQRT(g*h)
	
    IF (d==1) THEN
    ! build eigenvalues in x direction
	    Eigenvalues(1) = u-c
	    Eigenvalues(2) = u
	    Eigenvalues(3) = u+c
    ELSE
    ! build eigenvalues in y direction
	    Eigenvalues(1) = v-c
	    Eigenvalues(2) = v
	    Eigenvalues(3) = v+c
	END IF
	
	END FUNCTION
	
	
	
	
	
	! This routine builds the preconditioner for the shallow water system
	SUBROUTINE BuildShallowWaterPreconditioner (rmatrixBlockP, &
	                rarrayP, rarraySol, p_CXdata, p_CYdata, &
	                p_MLdata, p_Kdiagonal, p_kedge, &
	                NEQ, nedge, theta, dt, gravconst)
	
	! parameter values
	TYPE(t_matrixBlock), INTENT(INOUT)              :: rmatrixBlockP
	TYPE(t_array), DIMENSION(nvar2d), INTENT(INOUT) :: rarrayP
	TYPE(t_array), DIMENSION(nvar2d), INTENT(IN)    :: rarraySol
	REAL(DP), DIMENSION(:), POINTER, INTENT(IN)     :: p_CXdata, p_CYdata, p_MLdata
	INTEGER, DIMENSION(:), POINTER, INTENT(IN)      :: p_Kdiagonal
	INTEGER(I32), DIMENSION(:,:), POINTER           :: p_Kedge
	INTEGER, INTENT(IN)                             :: NEQ,	nedge
	REAL(DP), INTENT(IN)                            :: theta, gravconst, dt
	
	! variables
	INTEGER                             :: i, j, l, ii, jj, ij, ji, iedge, ivar
	REAL(DP), DIMENSION(nvar2d)         :: deltaQij, Qi, Qj, Qroeij
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: JacobixRoeij, JacobiyRoeij
	REAL(DP)                            :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: Aij, Bij, Dij, Eye
	REAL(DP)                            :: scalarDissipation, scalefactor, lambda
	REAL(DP)                            :: cRoe, uRoe, vRoe
	
	
	! Assemble the preconditioner rmatrixBlockP: P = ML - theta*dt*L
	! As we use the block jacobi method, we only need the main diagonal blocks
	
	    ! unit matrix
	    Eye = 0.0_DP
	    FORALL (ivar = 1: nvar2d)
            Eye(ivar,ivar) = 1.0_DP
	    END FORALL
	
		! First set all matrix blocks on the main diagonal of P equal to ML
		CALL lsysbl_clearMatrix(rmatrixBlockP)
		DO i = 1, NEQ
			ii = p_Kdiagonal(i)
			DO ivar = 1, nvar2d
				rarrayP(ivar)%Da(ii) = p_MLdata(i)
			END DO
		END DO
		
		! Now walk over all edges ij and compute Aij, Bij and Dij
		! Then write the entries of these Matrices to the corresponding entries
		! of the matrices on the main diagonal blocks of P
		DO iedge = 1, nedge
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
			scalefactor = SQRT(JcoeffxA**2.0_DP+JcoeffyA**2.0_DP)
			cRoe = SQRT(gravconst*Qroeij(1))    ! c_Roe=sqrt(g*h)
			cRoe = SQRT(0.5*gravconst*(Qi(1)+Qj(1)))
			uRoe = Qroeij(2)/Qroeij(1)
			vRoe = Qroeij(3)/Qroeij(1)
			lambda = MAX(ABS(uRoe-cRoe),ABS(uRoe+cRoe),ABS(vRoe-cRoe),ABS(vRoe+cRoe),ABS(uRoe),ABS(vRoe))
			scalarDissipation = scalefactor*lambda
			
			scalarDissipation = &
			ABS(JcoeffxA*MAXVAL(ABS(buildeigenvalues(Qroeij,1,gravconst))))+ &
			ABS(JcoeffyA*MAXVAL(ABS(buildeigenvalues(Qroeij,2,gravconst))))
			
			! compute Dij
			Dij = scalarDissipation*Eye
			
			! Now add the entries of Aij and Dij to their corresponding entries in P
			! P = M^L - theta*dt*L
			DO l = 1, nvar2d
				rarrayP(l)%Da(ii) = rarrayP(l)%Da(ii) - dt*theta*(+Aij(l,l)-Dij(l,l))
				rarrayP(l)%Da(ij) = rarrayP(l)%Da(ij) - dt*theta*(-Aij(l,l)+Dij(l,l))
				rarrayP(l)%Da(ji) = rarrayP(l)%Da(ji) - dt*theta*(+Aij(l,l)+Dij(l,l))
				rarrayP(l)%Da(jj) = rarrayP(l)%Da(jj) - dt*theta*(-Aij(l,l)-Dij(l,l))
			END DO
			                        
		END DO
	
	END SUBROUTINE
	
	
	
	
	
	
	! This routine builds the RHS for the shallow water system
	SUBROUTINE BuildShallowWaterRHS (&
	                rarrayRhs, rarraySol, rrhsBlock, rsolBlock, &
	                rmatrixML, p_CXdata, p_CYdata, p_MLdata, &
	                h_fld1, p_fld1, p_fld2, &
	                p_Kdiagonal, p_Kedge, NEQ, nedge, &
	                theta, dt, gravconst, Method, limiter)
	
	! parameter values
    TYPE(t_array), DIMENSION(nvar2d), INTENT(INOUT) :: rarrayRhs
	TYPE(t_array), DIMENSION(nvar2d), INTENT(IN)    :: rarraySol
	TYPE(t_vectorBlock), INTENT(INOUT)              :: rrhsBlock
	TYPE(t_vectorBlock), INTENT(IN)                 :: rsolBlock
	TYPE(t_matrixScalar), INTENT(IN)                :: rmatrixML
	INTEGER                                         :: h_fld1
	REAL(DP), DIMENSION(:,:), POINTER               :: p_fld1, p_fld2
	REAL(DP), DIMENSION(:), POINTER, INTENT(IN)     :: p_CXdata, p_CYdata, p_MLdata
	INTEGER, DIMENSION(:), POINTER, INTENT(IN)      :: p_Kdiagonal
	INTEGER(I32), DIMENSION(:,:), POINTER           :: p_Kedge
	INTEGER, INTENT(IN)                             :: NEQ,	nedge
	REAL(DP), INTENT(IN)                            :: theta, gravconst, dt
	INTEGER, INTENT(IN)                             :: Method, limiter
	
	
	! variables
	INTEGER                             :: i, j, l, ii, jj, ij, ji, d, iedge, ivar, inode
	REAL(DP), DIMENSION(nvar2d)         :: deltaQij, Qi, Qj, Qroeij
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: JacobixRoeij, JacobiyRoeij
	REAL(DP)                            :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: Aij, Bij, Dij, Eye
	REAL(DP)                            :: scalarDissipation, scalefactor, lambda
	REAL(DP), DIMENSION(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj
	REAL(DP)                            :: cRoe, uRoe, vRoe
	    ! for TVD
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: invRij, Rij
	REAL(DP), DIMENSION(nvar2d)         :: deltaWij, eigenvalues, deltaFij
	REAL(DP)                            :: scalefactor1, scalefactor2, deltaak, deltaWijHat
	REAL(DP)                            :: deltaaplus, deltaaminus, deltabplus, deltabminus
	INTEGER                             :: upwindnode
	
	
	    ! unit matrix
	    Eye = 0.0_DP
	    FORALL (ivar = 1: nvar2d)
            Eye(ivar,ivar) = 1.0_DP
	    END FORALL
	    
	    
	    ! Compute the RHS b: b=ML*Q^n
        DO ivar = 1, nvar2d
            CALL lsyssc_scalarMatVec (rmatrixML, &
                                      rsolBlock%Rvectorblock(ivar), &
                                      rrhsBlock%Rvectorblock(ivar), &
                                      1._DP, 0._DP)
        END DO
	
		
		! Now walk over all edges ij and compute Aij, Bij and Dij
		! and update the RHS b with the high/low order flux:
        ! b = b + dt*(1-theta)*K*u + dt*(1-theta)*D*u
		DO iedge = 1, nedge
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
			deltaKi = MATMUL(Aij+Bij,deltaQij)
			deltaKj = MATMUL(Aij-Bij,deltaQij)
			
			! Now choose the artificial diffusion method
			IF ((Method == 0).OR.(Method == 3)) THEN
			    ! Here we do not add artificial diffusion - so we have the high order method
			    Dij = 0
			ELSE IF ((Method == 1).OR.(Method == 4)) THEN
			    ! Here we use scalar dissipation
			    ! compute the maximum of the eigenvalues of Aij
			    scalefactor = SQRT(JcoeffxA**2.0_DP+JcoeffyA**2.0_DP)
			    cRoe = SQRT(gravconst*Qroeij(1))    ! c_Roe=sqrt(g*h)
			    cRoe = SQRT(0.5*gravconst*(Qi(1)+Qj(1)))
			    uRoe = Qroeij(2)/Qroeij(1)
			    vRoe = Qroeij(3)/Qroeij(1)
			    lambda = MAX(ABS(uRoe-cRoe),ABS(uRoe+cRoe),ABS(vRoe-cRoe),ABS(vRoe+cRoe),ABS(uRoe),ABS(vRoe))
			    scalarDissipation = scalefactor*lambda
			    
			    scalarDissipation = &
			ABS(JcoeffxA*MAXVAL(ABS(buildeigenvalues(Qroeij,1,gravconst))))+ &
			ABS(JcoeffyA*MAXVAL(ABS(buildeigenvalues(Qroeij,2,gravconst))))
			    
			    ! compute Dij
			    Dij = scalarDissipation*Eye
			ELSE IF ((Method == 2).OR.(Method == 5)) THEN
			    ! Here we use tensorial dissipation
			    Dij = ABS(JcoeffxA)*MATMUL(buildTrafo(Qroeij,1,gravconst),&
			                    MATMUL(buildaLambda(Qroeij,1,gravconst),&
			                           buildinvTrafo(Qroeij,1,gravconst))) +&
			          ABS(JcoeffyA)*MATMUL(buildTrafo(Qroeij,2,gravconst),&
			                    MATMUL(buildaLambda(Qroeij,2,gravconst),&
			                           buildinvTrafo(Qroeij,2,gravconst)))
			END IF
			
			! deltaD
			deltaDi = MATMUL(Dij-Bij,-deltaQij)
			deltaDj = -deltaDi
			
			! add deltaK and deltaD to rhs
			! rhs = rhs + (1-theta)*dt*K*u + (1-theta)*dt*D*u
			DO l = 1, nvar2d
				rarrayRhs(l)%Da(i) = rarrayRhs(l)%Da(i) + (1-theta)*dt*deltaKi(l) &
			                        	+ (1-theta)*dt*deltaDi(l) 
				rarrayRhs(l)%Da(j) = rarrayRhs(l)%Da(j) + (1-theta)*dt*deltaKj(l) &
			                        	+ (1-theta)*dt*deltaDj(l)
			END DO
			
		END DO
		
		
	    ! If we use the TVD scheme, then add the limited diffusive flux to the RHS
	    IF (Method == 3) THEN
        ! Do this for every dimension (dimensional splitting)
        DO d = 1, 2
      
        ! first we fill the array fld2
        DO iedge = 1, nedge
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
            deltaWij = MATMUL(invRij,-deltaQij)
            
            ! compute eigenvalues
            eigenvalues = buildEigenvalues(Qroeij,d,gravconst)
            
            ! compute scalefactors (a_ij^d and bij^d)
            ! Compute the coeffitients for the jacobi matrices
			JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
			JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
			JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
			JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP
			
			! calculate abs(aij^d)
			IF (d==1) THEN
			    scalefactor1 = JcoeffxA
			    scalefactor2 = JcoeffxB
			ELSE
			    scalefactor1 = JcoeffyA
			    scalefactor2 = JcoeffyB
			END IF
			                    
            
            ! write entries to fld2
            DO ivar = 1, nvar2d
                p_fld2(1+4*(ivar-1),iedge) = eigenvalues(ivar)               ! lambda_k
                p_fld2(2+4*(ivar-1),iedge) = deltaWij(ivar)                  ! deltaW^k
                p_fld2(3+4*(ivar-1),iedge) = -eigenvalues(ivar)*scalefactor1  ! k_ij^a
                p_fld2(4+4*(ivar-1),iedge) = -eigenvalues(ivar)*scalefactor2  ! k_ij^b
            END DO
        END DO
        
        
        ! next we fill the array fld1
        ! reset fld1
        CALL storage_clear (h_fld1)
        ! first compute the Piplus/minus, Qiplus/minus
        DO iedge = 1, nedge
			i  = p_Kedge(1, iedge)
			j  = p_Kedge(2, iedge)
			ij = p_Kedge(3, iedge)
			ji = p_Kedge(4, iedge)
			ii = p_Kdiagonal(i)
			jj = p_Kdiagonal(j)
			
			
			DO ivar = 1, nvar2d
			    deltaaplus = MAX(0.0_DP, p_fld2(3+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
			    deltaaminus = MIN(0.0_DP, p_fld2(3+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
			    deltabplus = MAX(0.0_DP, p_fld2(4+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
			    deltabminus = MIN(0.0_DP, p_fld2(4+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
			    
			    ! write entries to fld1
			    IF (p_fld2(3+4*(ivar-1),iedge)<0.0_DP) THEN
			    
			        p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + deltaaplus  ! Pi+
                    p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + deltaaminus ! Pi-
                    
                    p_fld1(3+6*(ivar-1),j) = p_fld1(3+6*(ivar-1),j) + deltaaplus  ! Qj+
                    p_fld1(4+6*(ivar-1),j) = p_fld1(4+6*(ivar-1),j) + deltaaminus ! Qj-
			    ELSE
			        p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) + deltaaplus
                    p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) + deltaaminus
                    
                    p_fld1(3+6*(ivar-1),i) = p_fld1(3+6*(ivar-1),i) + deltaaplus
                    p_fld1(4+6*(ivar-1),i) = p_fld1(4+6*(ivar-1),i) + deltaaminus
			    END IF
			    IF (p_fld2(4+4*(ivar-1),iedge)<0.0_DP) THEN
			    
			        p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + deltabplus
                    p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + deltabminus
                    
                    p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) - deltabplus
                    p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) - deltabminus
	    		ELSE
		    	    p_fld1(4+6*(ivar-1),j) = p_fld1(4+6*(ivar-1),j) - deltabplus
                    p_fld1(3+6*(ivar-1),j) = p_fld1(3+6*(ivar-1),j) - deltabminus
                    
                    p_fld1(3+6*(ivar-1),i) = p_fld1(3+6*(ivar-1),i) + deltabplus
                    p_fld1(4+6*(ivar-1),i) = p_fld1(4+6*(ivar-1),i) + deltabminus
    			END IF
            END DO
        END DO
        
        ! now compute the Riplus/minus
        DO inode = 1, NEQ
            DO ivar = 1, nvar2d
                p_fld1(5+6*(ivar-1),inode) = & ! Ri+
                         limiterfunc(p_fld1(3+6*(ivar-1),inode),p_fld1(1+6*(ivar-1),inode),limiter)
                p_fld1(6+6*(ivar-1),inode) = & ! Ri-
                         limiterfunc(p_fld1(4+6*(ivar-1),inode),p_fld1(2+6*(ivar-1),inode),limiter)
            END DO
        END DO
        
        ! next update the deltaWij = deltaWij - \hat deltaWij (apply limiting)
        DO iedge = 1, nedge
			i  = p_Kedge(1, iedge)
			j  = p_Kedge(2, iedge)
			ij = p_Kedge(3, iedge)
			ji = p_Kedge(4, iedge)
			ii = p_Kdiagonal(i)
			jj = p_Kdiagonal(j)
			
			DO ivar = 1, nvar2d
			    IF (p_fld2(3+4*(ivar-1),iedge)>0.0_DP) THEN
			        upwindnode = j
			    ELSE
			        upwindnode = i
			    END IF
			    
			    ! deltaak =  kij^a * deltaWij^k
			    deltaak = p_fld2(3+4*(ivar-1),iedge) * p_fld2(2+4*(ivar-1),iedge)
			    
			    ! deltaWijHat = R(upwind)+/- * deltaWij^k
			    IF (deltaak<0.0_DP) THEN
			        deltaWijHat = p_fld1(6+6*(ivar-1),upwindnode)*p_fld2(2+4*(ivar-1),iedge)
			    ELSE
			        deltaWijHat = p_fld1(5+6*(ivar-1),upwindnode)*p_fld2(2+4*(ivar-1),iedge)
			    END IF
			    
			    ! deltaWij = deltaWij - \hat deltaWij
			    p_fld2(2+4*(ivar-1),iedge) = p_fld2(2+4*(ivar-1),iedge) - deltaWijHat
			END DO	
	    END DO
	    
	    ! finally compute the fluxes Fij and update the rhs
	    DO iedge = 1, nedge
			i  = p_Kedge(1, iedge)
			j  = p_Kedge(2, iedge)
			ij = p_Kedge(3, iedge)
			ji = p_Kedge(4, iedge)
			ii = p_Kdiagonal(i)
			jj = p_Kdiagonal(j)
			
			DO ivar = 1, nvar2d
			    deltaWij(ivar) = p_fld2(2+4*(ivar-1),iedge)
			END DO
			
			
			! get solution values at node i and j
	    	Qi = (/rarraySol(1)%Da(i),rarraySol(2)%Da(i),rarraySol(3)%Da(i)/)
		    Qj = (/rarraySol(1)%Da(j),rarraySol(2)%Da(j),rarraySol(3)%Da(j)/)
						
			! Compute the Roe-meanvalue for this edge
			Qroeij = calculateQroe(Qi, Qj)
		
			! Compute the coeffitients for the jacobi matrices
			JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
			JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
			
			! calculate abs(aij^d)
			IF (d==1) THEN
			    scalefactor = ABS(JcoeffxA)
			ELSE
			    scalefactor = ABS(JcoeffyA)
			END IF
			
			! compute the (limited) diffusive flux
			deltaFij = scalefactor*MATMUL(buildTrafo(Qroeij,d,gravconst), &
			         MATMUL(buildaLambda(Qroeij,d,gravconst),deltaWij))
			         
			
			! add the flux F to the rhs
			! rhs(i) = rhs(i) + (1-theta)*dt*deltaF(i)
			! rhs(j) = rhs(j) - (1-theta)*dt*deltaF(j)
			DO l = 1, nvar2d
				rarrayRhs(l)%Da(i) = rarrayRhs(l)%Da(i) + (1-theta)*dt*deltaFij(l)
				rarrayRhs(l)%Da(j) = rarrayRhs(l)%Da(j) - (1-theta)*dt*deltaFij(l)
			END DO
			
	    END DO
        
        END DO
		END IF
		
	
	END SUBROUTINE
	
	
	
	
	
	
	
	
	
	
	
	! This routine builds the defect for the shallow water system
	SUBROUTINE BuildShallowWaterDefect (&
	                rdefBlock, rstempBlock, rrhsBlock, rsolBlock, &
	                rarrayRhs, rarraySol, rarrayRstemp, &
	                rmatrixML, p_CXdata, p_CYdata, p_MLdata, &
	                h_fld1, p_fld1, p_fld2, &
	                p_Kdiagonal, p_Kedge, NEQ, nedge, &
	                theta, dt, gravconst, Method, limiter)
	
	! parameter values
    TYPE(t_array), DIMENSION(nvar2d), INTENT(INOUT) :: rarrayRstemp
	TYPE(t_array), DIMENSION(nvar2d), INTENT(IN)    :: rarraySol, rarrayRhs
	TYPE(t_vectorBlock), INTENT(INOUT)              :: rdefBlock, rstempBlock
	TYPE(t_vectorBlock), INTENT(IN)                 :: rsolBlock, rrhsBlock
	TYPE(t_matrixScalar), INTENT(IN)                :: rmatrixML
	REAL(DP), DIMENSION(:), POINTER, INTENT(IN)     :: p_CXdata, p_CYdata, p_MLdata
	INTEGER                                         :: h_fld1
	REAL(DP), DIMENSION(:,:), POINTER               :: p_fld1, p_fld2
	INTEGER, DIMENSION(:), POINTER, INTENT(IN)      :: p_Kdiagonal
	INTEGER(I32), DIMENSION(:,:), POINTER           :: p_Kedge
	INTEGER, INTENT(IN)                             :: NEQ,	nedge
	REAL(DP), INTENT(IN)                            :: theta, gravconst, dt
	INTEGER, INTENT(IN)                             :: Method, limiter
		
	! variables
	INTEGER                             :: i, j, l, ii, jj, ij, ji, d, iedge, ivar, inode
	REAL(DP), DIMENSION(nvar2d)         :: deltaQij, Qi, Qj, Qroeij
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: JacobixRoeij, JacobiyRoeij
	REAL(DP)                            :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: Aij, Bij, Dij, Eye
	REAL(DP)                            :: scalarDissipation, scalefactor, lambda
	REAL(DP), DIMENSION(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj
	REAL(DP)                            :: cRoe, uRoe, vRoe
	    ! for TVD
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: invRij, Rij
	REAL(DP), DIMENSION(nvar2d)         :: deltaWij, eigenvalues, deltaFij
	REAL(DP)                            :: scalefactor1, scalefactor2, deltaak, deltaWijHat
	REAL(DP)                            :: deltaaplus, deltaaminus, deltabplus, deltabminus
	INTEGER                             :: upwindnode
	
	
	    ! unit matrix
	    Eye = 0.0_DP
	    FORALL (ivar = 1: nvar2d)
            Eye(ivar,ivar) = 1.0_DP
	    END FORALL
	    
	    
	    ! Now we want to compute the defect vector
            ! def = rhs - (ML*Q - theta*dt*K*u - theta*dt*D*u) = rhs - rstemp
            
            ! calculate: rstemp = ML*Q
            DO ivar = 1, nvar2d
                CALL lsyssc_scalarMatVec (rmatrixML, &
                                          rsolBlock%Rvectorblock(ivar), &
                                          rstempBlock%Rvectorblock(ivar), &
                                          1._DP, 0._DP)
            END DO
            
            
            ! Now add the fluxes corresponding to the operator K and D,
            ! so we get the low order scheme:
            ! rstemp = rstemp - theta*dt*K*u -theta*dt*D*u
            ! So walk over all edges ij and compute Aij and Bij
		    ! Then compute deltaKi=(Aij+Bij)*(Qi-Qj)
		    ! and deltaKj=(Aij-Bij)*(Qi-Qj)
		    ! deltaDi = (Dij+Bij)*(Qj-Qi)
		    ! deltaDj = -deltaDi
		    ! and add this to rstemp
		    DO iedge = 1, nedge
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
			    deltaKi = MATMUL(Aij+Bij,deltaQij)
			    deltaKj = MATMUL(Aij-Bij,deltaQij)
			    
			    ! Now choose the artificial diffusion method
			    IF ((Method == 0).OR.(Method == 3)) THEN
			        ! Here we do not add artificial diffusion - so we have the high order method
			        Dij = 0
			    ELSE IF ((Method == 1).OR.(Method == 4)) THEN
			        ! Here we use scalar dissipation
			        ! compute the maximum of the eigenvalues of Aij
			        scalefactor = SQRT(JcoeffxA**2.0_DP+JcoeffyA**2.0_DP)
			        cRoe = SQRT(gravconst*Qroeij(1))    ! c_Roe=sqrt(g*h)
			        cRoe = SQRT(0.5*gravconst*(Qi(1)+Qj(1)))
			        uRoe = Qroeij(2)/Qroeij(1)
			        vRoe = Qroeij(3)/Qroeij(1)
			        lambda = MAX(ABS(uRoe-cRoe),ABS(uRoe+cRoe),ABS(vRoe-cRoe),ABS(vRoe+cRoe),ABS(uRoe),ABS(vRoe))
			        scalarDissipation = scalefactor*lambda
			        
			        scalarDissipation = &
			ABS(JcoeffxA*MAXVAL(ABS(buildeigenvalues(Qroeij,1,gravconst))))+ &
			ABS(JcoeffyA*MAXVAL(ABS(buildeigenvalues(Qroeij,2,gravconst))))
			        
			        ! compute Dij
			        Dij = scalarDissipation*Eye
			    ELSE IF ((Method == 2).OR.(Method == 5)) THEN
			        ! Here we use tensorial dissipation
			        Dij = ABS(JcoeffxA)*MATMUL(buildTrafo(Qroeij,1,gravconst),&
			                        MATMUL(buildaLambda(Qroeij,1,gravconst),&
			                               buildinvTrafo(Qroeij,1,gravconst))) +&
			              ABS(JcoeffyA)*MATMUL(buildTrafo(Qroeij,2,gravconst),&
			                        MATMUL(buildaLambda(Qroeij,2,gravconst),&
			                               buildinvTrafo(Qroeij,2,gravconst)))
			    END IF
			    
			    ! deltaD
			    deltaDi = MATMUL(Dij-Bij,-deltaQij)
			    deltaDj = -deltaDi
			    
			    ! add deltaK and deltaD to rstemp
			    ! rstemp = rstemp - theta*dt*K*u - theta*dt*D*u
				DO l = 1, nvar2d
					rarrayRstemp(l)%Da(i) = rarrayRstemp(l)%Da(i) - theta*dt*deltaKi(l) &
				 			                                	  - theta*dt*deltaDi(l)
					rarrayRstemp(l)%Da(j) = rarrayRstemp(l)%Da(j) - theta*dt*deltaKj(l) &
					 			                                  - theta*dt*deltaDj(l)
				END DO
				
			END DO
			
			
			
			
		! If we use the TVD scheme, then add the limited diffusive flux to the Rstemp (defect)
	    IF (Method == 3) THEN
        ! Do this for every dimension (dimensional splitting)
        DO d = 1, 2
      
        ! first we fill the array fld2
        DO iedge = 1, nedge
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
            deltaWij = MATMUL(invRij,-deltaQij)
            
            ! compute eigenvalues
            eigenvalues = buildEigenvalues(Qroeij,d,gravconst)
            
            ! compute scalefactors (a_ij^d and bij^d)
            ! Compute the coeffitients for the jacobi matrices
			JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
			JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
			JcoeffxB = (p_CXdata(ij)+p_CXdata(ji))/2.0_DP
			JcoeffyB = (p_CYdata(ij)+p_CYdata(ji))/2.0_DP
			
			! calculate abs(aij^d)
			IF (d==1) THEN
			    scalefactor1 = JcoeffxA
			    scalefactor2 = JcoeffxB
			ELSE
			    scalefactor1 = JcoeffyA
			    scalefactor2 = JcoeffyB
			END IF
			                    
            
            ! write entries to fld2
            DO ivar = 1, nvar2d
                p_fld2(1+4*(ivar-1),iedge) = eigenvalues(ivar)               ! lambda_k
                p_fld2(2+4*(ivar-1),iedge) = deltaWij(ivar)                  ! deltaW^k
                p_fld2(3+4*(ivar-1),iedge) = -eigenvalues(ivar)*scalefactor1  ! k_ij^a
                p_fld2(4+4*(ivar-1),iedge) = -eigenvalues(ivar)*scalefactor2  ! k_ij^b
            END DO
        END DO
        
        
        ! next we fill the array fld1
        ! reset fld1
        CALL storage_clear (h_fld1)
        ! first compute the Piplus/minus, Qiplus/minus
        DO iedge = 1, nedge
			i  = p_Kedge(1, iedge)
			j  = p_Kedge(2, iedge)
			ij = p_Kedge(3, iedge)
			ji = p_Kedge(4, iedge)
			ii = p_Kdiagonal(i)
			jj = p_Kdiagonal(j)
			
			
			DO ivar = 1, nvar2d
			    deltaaplus = MAX(0.0_DP, p_fld2(3+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
			    deltaaminus = MIN(0.0_DP, p_fld2(3+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
			    deltabplus = MAX(0.0_DP, p_fld2(4+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
			    deltabminus = MIN(0.0_DP, p_fld2(4+4*(ivar-1),iedge)*p_fld2(2+4*(ivar-1),iedge))
			    
			    ! write entries to fld1
			    IF (p_fld2(3+4*(ivar-1),iedge)<0.0_DP) THEN
			    
			        p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + deltaaplus  ! Ri+
                    p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + deltaaminus ! Ri-
                    
                    p_fld1(3+6*(ivar-1),j) = p_fld1(3+6*(ivar-1),j) + deltaaplus  ! Qj+
                    p_fld1(4+6*(ivar-1),j) = p_fld1(4+6*(ivar-1),j) + deltaaminus ! Qj-
			    ELSE
			        p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) + deltaaplus
                    p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) + deltaaminus
                    
                    p_fld1(3+6*(ivar-1),i) = p_fld1(3+6*(ivar-1),i) + deltaaplus
                    p_fld1(4+6*(ivar-1),i) = p_fld1(4+6*(ivar-1),i) + deltaaminus
			    END IF
			    IF (p_fld2(4+4*(ivar-1),iedge)<0.0_DP) THEN
			    
			        p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + deltabplus
                    p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + deltabminus
                    
                    p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) - deltabplus
                    p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) - deltabminus
	    		ELSE
		    	    p_fld1(4+6*(ivar-1),j) = p_fld1(4+6*(ivar-1),j) - deltabplus
                    p_fld1(3+6*(ivar-1),j) = p_fld1(3+6*(ivar-1),j) - deltabminus
                    
                    p_fld1(3+6*(ivar-1),i) = p_fld1(3+6*(ivar-1),i) + deltabplus
                    p_fld1(4+6*(ivar-1),i) = p_fld1(4+6*(ivar-1),i) + deltabminus
    			END IF
            END DO
        END DO
        
        ! now compute the Riplus/minus
        DO inode = 1, NEQ
            DO ivar = 1, nvar2d
                p_fld1(5+6*(ivar-1),inode) = &
                         limiterfunc(p_fld1(3+6*(ivar-1),inode),p_fld1(1+6*(ivar-1),inode),limiter)
                p_fld1(6+6*(ivar-1),inode) = &
                         limiterfunc(p_fld1(4+6*(ivar-1),inode),p_fld1(2+6*(ivar-1),inode),limiter)
            END DO
        END DO
        
        ! next update the deltaWij = deltaWij - \hat deltaWij (apply limiting)
        DO iedge = 1, nedge
			i  = p_Kedge(1, iedge)
			j  = p_Kedge(2, iedge)
			ij = p_Kedge(3, iedge)
			ji = p_Kedge(4, iedge)
			ii = p_Kdiagonal(i)
			jj = p_Kdiagonal(j)
			
			DO ivar = 1, nvar2d
			    IF (p_fld2(3+4*(ivar-1),iedge)>0.0_DP) THEN
			        upwindnode = j
			    ELSE
			        upwindnode = i
			    END IF
			    
			    ! deltaak =  kij^a * deltaWij^k
			    deltaak = p_fld2(3+4*(ivar-1),iedge) * p_fld2(2+4*(ivar-1),iedge)
			    
			    ! deltaWijHat = R(upwind)+/- * deltaWij^k
			    IF (deltaak<0.0_DP) THEN
			        deltaWijHat = p_fld1(6+6*(ivar-1),upwindnode)*p_fld2(2+4*(ivar-1),iedge)
			    ELSE
			        deltaWijHat = p_fld1(5+6*(ivar-1),upwindnode)*p_fld2(2+4*(ivar-1),iedge)
			    END IF
			    
			    ! deltaWij = deltaWij - \hat deltaWij
			    p_fld2(2+4*(ivar-1),iedge) = p_fld2(2+4*(ivar-1),iedge) - deltaWijHat
			END DO	
	    END DO
	    
	    ! finally compute the fluxes Fij and update the rhs
	    DO iedge = 1, nedge
			i  = p_Kedge(1, iedge)
			j  = p_Kedge(2, iedge)
			ij = p_Kedge(3, iedge)
			ji = p_Kedge(4, iedge)
			ii = p_Kdiagonal(i)
			jj = p_Kdiagonal(j)
			
			DO ivar = 1, nvar2d
			    deltaWij(ivar) = p_fld2(2+4*(ivar-1),iedge)
			END DO
			
			
			! get solution values at node i and j
	    	Qi = (/rarraySol(1)%Da(i),rarraySol(2)%Da(i),rarraySol(3)%Da(i)/)
		    Qj = (/rarraySol(1)%Da(j),rarraySol(2)%Da(j),rarraySol(3)%Da(j)/)
						
			! Compute the Roe-meanvalue for this edge
			Qroeij = calculateQroe(Qi, Qj)
		
			! Compute the coeffitients for the jacobi matrices
			JcoeffxA = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
			JcoeffyA = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
			
			! calculate abs(aij^d)
			IF (d==1) THEN
			    scalefactor = ABS(JcoeffxA)
			ELSE
			    scalefactor = ABS(JcoeffyA)
			END IF
			
			! compute the (limited) diffusive flux
			deltaFij = scalefactor*MATMUL(buildTrafo(Qroeij,d,gravconst), &
			         MATMUL(buildaLambda(Qroeij,d,gravconst),deltaWij))
			         
			
			! add the flux F to the Rstemp (defect)
			! rstemp(i) = rstemp(i) - theta*dt*deltaF
    		! rstemp(j) = rstemp(j) + theta*dt*deltaF
	    	DO l = 1, nvar2d
		    	rarrayRstemp(l)%Da(i) = rarrayRstemp(l)%Da(i) - theta*dt*deltaFij(l)
			   	rarrayRstemp(l)%Da(j) = rarrayRstemp(l)%Da(j) + theta*dt*deltaFij(l)
    		END DO
			
	    END DO
        
        END DO
		END IF
		
		
            ! Finally we can assemble the defect vector:
            ! def = rhs - rstemp
            CALL lsysbl_vectorLinearComb (rrhsBlock,rstempBlock, &
                                          1.0_DP,-1.0_DP,rdefBlock)
    
	END SUBROUTINE	
	
	
	
	
	
	! linearised FCT Method for Shallow Water
	! adding limited antidiffusion to the low order predictor
	SUBROUTINE FctShallowWaterAddLimitedAntidiffusion(&
	                rarraySol, rarraySolDot, rarrayRhs,&
	                rdefBlock, rstempBlock, rsolBlock, rSolDotBlock, &
	                rmatrixML, p_CXdata, p_CYdata, p_MLdata, p_MCdata, &
                    h_fld1, p_fld1, p_fld2, &
                    p_Kdiagonal, p_Kedge, NEQ, nedge, &
                    gravconst, dt, Method, prelimiting)
	
	! parameter values
    TYPE(t_array), DIMENSION(nvar2d), INTENT(INOUT) :: rarraySol, rarraySolDot
	TYPE(t_array), DIMENSION(nvar2d), INTENT(IN)    :: rarrayRhs
	TYPE(t_vectorBlock), INTENT(INOUT)              :: rdefBlock, rstempBlock
	TYPE(t_vectorBlock), INTENT(IN)                 :: rsolBlock
	TYPE(t_vectorBlock), INTENT(INOUT)              :: rSolDotBlock
	TYPE(t_matrixScalar), INTENT(IN)                :: rmatrixML
	REAL(DP), DIMENSION(:), POINTER, INTENT(IN)     :: p_CXdata, p_CYdata, p_MLdata, p_MCdata
	INTEGER                                         :: h_fld1
	REAL(DP), DIMENSION(:,:), POINTER               :: p_fld1, p_fld2
	INTEGER, DIMENSION(:), POINTER, INTENT(IN)      :: p_Kdiagonal
	INTEGER(I32), DIMENSION(:,:), POINTER           :: p_Kedge
	INTEGER, INTENT(IN)                             :: NEQ,	nedge
	REAL(DP), INTENT(IN)                            :: gravconst, dt
	INTEGER, INTENT(IN)                             :: Method, prelimiting
	
	
	! variables
	INTEGER                             :: i, j, l, ii, jj, ij, ji, d, iedge, ivar
	REAL(DP), DIMENSION(nvar2d)         :: deltaQij, Qi, Qj, Qroeij
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: JacobixRoeij, JacobiyRoeij
	REAL(DP)                            :: JcoeffxA, JcoeffyA, JcoeffxB, JcoeffyB
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: Aij, Bij, Dij, Eye
	REAL(DP)                            :: scalarDissipation, scalefactor, lambda
	REAL(DP), DIMENSION(nvar2d)         :: deltaKi, deltaKj, deltaDi, deltaDj
	REAL(DP)                            :: cRoe, uRoe, vRoe
	REAL(DP), DIMENSION(nvar2d)         :: Udoti, Udotj
	REAL(DP), DIMENSION(nvar2d)         :: deltaWij, deltaGij, deltaFij
	REAL(DP), DIMENSION(nvar2d,nvar2d)  :: Rij, invRij
	
	
	! unit matrix
	Eye = 0.0_DP
	FORALL (ivar = 1: nvar2d)
        Eye(ivar,ivar) = 1.0_DP
	END FORALL
	
    
    ! For the following use dimensional splitting
    DO d = 1, 2
    
        ! Clear SolDot
        CALL lsysbl_clearVector (rSolDotBlock)
        
        ! First compute SolDot
        DO iedge = 1, nedge
            i  = p_Kedge(1, iedge)
			j  = p_Kedge(2, iedge)
    		ij = p_Kedge(3, iedge)
	    	ji = p_Kedge(4, iedge)
		    ii = p_Kdiagonal(i)
			jj = p_Kdiagonal(j)
						
    		! get solution values at node i and j
		    DO ivar = 1, nvar2d
	    	    Qi(ivar) = rarraySol(ivar)%Da(i)
		        Qj(ivar) = rarraySol(ivar)%Da(j)
		    END DO
		    	
		    ! compute deltaQij = Qi - Qj
		    deltaQij = Qi - Qj
			
    		! Compute the Roe-meanvalue for this edge
	    	Qroeij = calculateQroe(Qi, Qj)
	    	
	    	IF (d==1) THEN
	    	    scalefactor = (p_CXdata(ij)-p_CXdata(ji))/2.0_DP
	    	ELSE
	    	    scalefactor = (p_CYdata(ij)-p_CYdata(ji))/2.0_DP
	    	END IF
			
    		! Now we can compute Aij and Bij
			Aij = scalefactor*buildJacobi(Qroeij,d,gravconst)
			
			! deltaK
			deltaKi = -MATMUL(Aij,deltaQij)
			deltaKj = deltaKi
			    
			! compute Dij
			IF (Method == 4) THEN
			    ! Build Dij using scalar dissipation
			    scalarDissipation = MAXVAL(ABS(buildEigenvalues(Qroeij,d,gravconst)))
			    Dij = abs(scalefactor)*scalarDissipation*Eye
			ELSE IF (Method == 5) THEN
			    ! Build Dij using tensorial dissipation
			    Dij = ABS(scalefactor)*MATMUL(buildTrafo(Qroeij,d,gravconst),&
			                    MATMUL(buildaLambda(Qroeij,d,gravconst),&
			                           buildinvTrafo(Qroeij,d,gravconst)))
			END IF
			    
			! deltaD
			deltaDi = MATMUL(Dij,-deltaQij)
			deltaDj = -deltaDi
			    
			! add deltaK and deltaD to SolDot
			! SolDot = SolDot + 1/ML*(K*u + D*u)
			DO l = 1, nvar2d
				rarraySolDot(l)%Da(i) = rarraySolDot(l)%Da(i) &
				                        + 1.0_DP/p_MLdata(i)*&
				                        (deltaKi(l) + deltaDi(l))
				rarraySolDot(l)%Da(j) = rarraySolDot(l)%Da(j) &
				                        + 1.0_DP/p_MLdata(j)*&
				                        (deltaKj(l) + deltaDj(l))
				! Save Qroe
                p_fld2(1+2*(l-1),iedge) = QRoeij(l)
			END DO
			
		END DO
		! SolDot fertig!
		
				
		
		! Now calculate the limited antidiffusive flux
        
        CALL storage_clear (h_fld1)
        ! first we fill the array fld1 + fld2 (that means compute deltaWij and deltaGij and save it)
        DO iedge = 1, nedge
			i  = p_Kedge(1, iedge)
			j  = p_Kedge(2, iedge)
			ij = p_Kedge(3, iedge)
			ji = p_Kedge(4, iedge)
			ii = p_Kdiagonal(i)
			jj = p_Kdiagonal(j)
			
			! get solution values at node i and j
			DO ivar = 1, nvar2d
	    	    Qi(ivar) = rarraySol(ivar)%Da(i)
		        Qj(ivar) = rarraySol(ivar)%Da(j)
		    END DO
		    
			! compute deltaQij = Qi - Qj
		    deltaQij = Qi - Qj
			
			! Compute the Roe-meanvalue for this edge
			Qroeij = calculateQroe(Qi, Qj)
			
            ! Rij^-1
            invRij = buildInvTrafo(Qroeij,d,gravconst)
            
            ! compute Rij^-1*(Qj - Qi)
            deltaWij = -MATMUL(invRij,deltaQij)
            
			IF (d==1) THEN
	    	    scalefactor = (p_CXdata(ij)-p_CXdata(ji))/2._DP
	    	ELSE
	    	    scalefactor = (p_CYdata(ij)-p_CYdata(ji))/2._DP
	    	END IF
            
            ! compute Dij
			IF (Method == 4) THEN
			    ! Build Dij using scalar dissipation
			    scalarDissipation = MAXVAL(ABS(buildEigenvalues(Qroeij,d,gravconst)))
			    Dij = abs(scalefactor)*scalarDissipation*Eye
			ELSE IF (Method == 5) THEN
			    ! Build Dij using tensorial dissipation
			    Dij = ABS(scalefactor)*MATMUL(buildTrafo(Qroeij,d,gravconst),&
			                    MATMUL(buildaLambda(Qroeij,d,gravconst),&
			                           buildinvTrafo(Qroeij,d,gravconst)))
			END IF
            
            ! get Udot at node i and j
            DO ivar = 1, nvar2d
                Udoti(ivar) = rarraySolDot(ivar)%Da(i)
                Udotj(ivar) = rarraySolDot(ivar)%Da(j)
            END DO
            
            ! compute deltaFij (73)
            deltaFij = p_MCdata(ij)*(Udoti-Udotj)+MATMUL(Dij,deltaQij)
            
            ! Now apply prelimiting
            IF (prelimiting == 1) THEN
                IF (Method==4) THEN
                    ! For scalar dissipation apply
                    ! MinMod prelimiting
                    DO ivar = 1, nvar2d
                    IF (deltaFij(ivar)*deltaQij(ivar)<0) THEN
                        deltaFij(ivar)=0
                    ELSE
                        IF (ABS(deltaFij(ivar))>ABS(ABS(scalefactor)*scalarDissipation*deltaQij(ivar))) THEN
                                deltaFij(ivar) = ABS(scalefactor)*scalarDissipation*deltaQij(ivar)
                        END IF
                    END IF
                    END DO
                ELSE IF (Method==5) THEN
                ! For tensorial dissipation apply
                ! Simple prelimiting
                    DO ivar = 1, nvar2d
                        IF (deltaFij(ivar)*deltaQij(ivar)<0) THEN
                            deltaFij(ivar) = 0
                        END IF
                    END DO
                END IF
            END IF
            
      
            ! compute deltaGij
            deltaGij = MATMUL(invRij,deltaFij)
            
            DO ivar = 1, nvar2d
                
                ! Save deltaGij
                p_fld2(2+2*(ivar-1),iedge) = deltaGij(ivar)
                               
                ! Update the P/Q_i/j +/- (76+77)
                p_fld1(1+6*(ivar-1),i) = p_fld1(1+6*(ivar-1),i) + MAX(0.0_DP,deltaGij(ivar)) ! Pi+
                p_fld1(2+6*(ivar-1),i) = p_fld1(2+6*(ivar-1),i) + MIN(0.0_DP,deltaGij(ivar)) ! Pi-
                p_fld1(3+6*(ivar-1),i) = MAX(p_fld1(3+6*(ivar-1),i),deltaWij(ivar))! Qi+
                p_fld1(4+6*(ivar-1),i) = MIN(p_fld1(4+6*(ivar-1),i),deltaWij(ivar))! Qi-
                
                p_fld1(1+6*(ivar-1),j) = p_fld1(1+6*(ivar-1),j) + MAX(0.0_DP,-deltaGij(ivar)) ! Pj+
                p_fld1(2+6*(ivar-1),j) = p_fld1(2+6*(ivar-1),j) + MIN(0.0_DP,-deltaGij(ivar)) ! Pj-
                p_fld1(3+6*(ivar-1),j) = MAX(p_fld1(3+6*(ivar-1),j),-deltaWij(ivar))! Qj+
                p_fld1(4+6*(ivar-1),j) = MIN(p_fld1(4+6*(ivar-1),j),-deltaWij(ivar))! Qj-
                
            END DO
            
        END DO
              
        DO i = 1, NEQ
        
            DO ivar = 1, nvar2d
                ! Compute the R_i +/- (78)
                IF (ABS(p_fld1(1+6*(ivar-1),i))>1e-8) THEN
                p_fld1(5+6*(ivar-1),i) = MIN(1.0_DP,p_MLdata(i)*&
                                         p_fld1(3+6*(ivar-1),i)/p_fld1(1+6*(ivar-1),i)/dt)! Ri+
                ELSE
                p_fld1(5+6*(ivar-1),i) = 0.0_DP
                END IF
                
                IF (ABS(p_fld1(2+6*(ivar-1),i))>1e-8) THEN
                p_fld1(6+6*(ivar-1),i) = MIN(1.0_DP, p_MLdata(i)*&
                                         p_fld1(4+6*(ivar-1),i)/p_fld1(2+6*(ivar-1),i)/dt)! Ri-
                ELSE
                p_fld1(6+6*(ivar-1),i) = 0.0_DP
                END IF
            
            END DO
            
            
            
        END DO
        
        DO iedge = 1, nedge
			i  = p_Kedge(1, iedge)
			j  = p_Kedge(2, iedge)
			ij = p_Kedge(3, iedge)
			ji = p_Kedge(4, iedge)
			ii = p_Kdiagonal(i)
			jj = p_Kdiagonal(j)
			
			DO ivar = 1, nvar2d
			    deltaGij(ivar) = p_fld2(2+2*(ivar-1),iedge)
			    
			    ! Limit the antidiffusive fluxes (79)
			    IF (deltaGij(ivar)>0.0_DP) THEN
			        deltaGij(ivar) = deltaGij(ivar) * MIN(p_fld1(5+6*(ivar-1),i),p_fld1(6+6*(ivar-1),j))
			    ELSE
			        deltaGij(ivar) = deltaGij(ivar) * MIN(p_fld1(5+6*(ivar-1),j),p_fld1(6+6*(ivar-1),i))
			    END IF
			    
			END DO
			
			! get Roe-values at ij
			DO ivar = 1, nvar2d
	    	    Qroeij(ivar) = p_fld2(1+2*(ivar-1),iedge)
		    END DO
		    			
            ! Rij
            Rij = buildTrafo(Qroeij,d,gravconst)
			
			! Finally we can transform back to deltaFij
			deltaFij = MATMUL(Rij,deltaGij)
			
			! And save 
			DO ivar = 1, nvar2d
			    p_fld2(2+2*(ivar-1),iedge) = deltaFij(ivar)
			END DO
			
		END DO
			
			
			
		DO iedge = 1, nedge
			i  = p_Kedge(1, iedge)
			j  = p_Kedge(2, iedge)
			ij = p_Kedge(3, iedge)
			ji = p_Kedge(4, iedge)
			ii = p_Kdiagonal(i)
			jj = p_Kdiagonal(j)
			
			DO ivar = 1, nvar2d
			    deltaFij(ivar) = p_fld2(2+2*(ivar-1),iedge)
			    
			    rarraySol(ivar)%Da(i) = rarraySol(ivar)%Da(i) + deltaFij(ivar)   *dt/p_MLdata(i)
			    rarraySol(ivar)%Da(j) = rarraySol(ivar)%Da(j) - deltaFij(ivar)   *dt/p_MLdata(j)
			END DO
			
	    END DO
            
    
    END DO ! d=1, 2
       
	END SUBROUTINE
	
	
	
	
	
	
	
	
	
	
	! This subroutine takes care of boundary conditions
	! for a reflacting boundary
	SUBROUTINE ImplementShallowWaterBCs (&
	                rboundary, rtriangulation, &
	                rarrayP, rarraySol, rarrayDef, &
	                p_Kdiagonal, p_Kld, &
	                gravconst, boundarycorner)
	                
	! PARAMETER VALUES
	! the boundary
	TYPE(t_boundary), INTENT(IN) :: rboundary
	! the triangulation
	TYPE(t_triangulation), INTENT(IN) :: rtriangulation
	! pointers to datas of preconditioner, solution and defect
	TYPE(t_array), DIMENSION(nvar2d), INTENT(INOUT) :: rarrayP, rarraySol, rarrayDef
	! pointer to index of diagonal positions of matrices
	INTEGER, DIMENSION(:), POINTER, INTENT(IN) :: p_Kdiagonal, p_Kld
	REAL(DP), INTENT(IN) :: gravconst
	INTEGER :: boundarycorner
	
	
	
	! VARIABLES
	! index and number of boundary components
	INTEGER :: ibct, nbct
	! Index of first and last node of this boundarycomponent
	INTEGER :: ivbdFirst, ivbdLast
	! Pointers to the triangultion datas
	INTEGER, DIMENSION(:), POINTER :: p_IboundaryCpIdx
	INTEGER, DIMENSION(:), POINTER :: p_IverticesAtBoundary
	REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
	! indices of boundary nodes
	INTEGER :: ivbd
	! the index of the boundary node to edit
	INTEGER :: i
	! index of diagonal element ii in matrix
	INTEGER :: ii
	! parametervalue of boundary node to edit
	REAL(DP) :: parmv
	! unit outward normal vector and tangential vector
	REAL(DP), DIMENSION(2) :: normalV, tangentialV
	! old primitive values at boundary node
	REAL(DP) :: h, u, v
	! predicted primitive star-values at boundary node
	REAL(DP) :: hstar, ustar, vstar
	! final primitive starstar-values at boundary node
	REAL(DP) :: hstarstar, ustarstar, vstarstar, hm, oldh
	! normal and tangential projection of the predicted velocity
	REAL(DP) :: vnproj, vtproj
	! for the fixpointiteration to calculate h**
	INTEGER :: ite
	! maximum iterations
	INTEGER, PARAMETER :: itemax = 50
	! Element ij in matrix (preconditioner)
	INTEGER :: ij
	
	
	! Initialise some values
	
	! get number of boundary components
	nbct = rtriangulation%NBCT
	
	! get pointer to IboundaryCpIdx, which saves where the boundarycomponents begin and end
	CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
	
	! get pointer to IverticesAtBoundary, which saves die numbers of the vertices on the boundary
	CALL storage_getbase_int (rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
	
	! get pointer to DvertexParameterValue, which saves the parametervalues 
	CALL storage_getbase_double (rtriangulation%h_DvertexParameterValue, p_DvertexParameterValue)
	
	
	! Now start manipulation of boundary values, preconditioner and defect
	
	! loop over all boundarycomponents
	DO ibct = 1, nbct
	    
	    ! get Index of first and last node of this boundarycomponent
	    ivbdFirst = p_IboundaryCpIdx(ibct)
	    ivbdLast  = p_IboundaryCpIdx(ibct+1)-1
	    
	    
	    DO ivbd = ivbdFirst, ivbdLast
	        ! the index of the boundary node to edit
	        i = p_IverticesAtBoundary(ivbd)
	        ! entry of diagonal element in matrix
	        ii = p_Kdiagonal(i)
	        ! parameter value of the boundary node to edit
	        parmv = p_DvertexParameterValue(ivbd)
	        

            ! Shall the nodes in a corner be manipulated, too?
	        IF (((ABS(parmv-AINT(PARMV))>1e-6)).OR.(boundarycorner==1)) THEN
	        
	        
	        ! get normal vector at boundary node
	        CALL boundary_getNormalVec2D(rboundary,ibct,parmv,normalV(1),normalV(2))
	        IF (ABS(parmv-AINT(PARMV))<1e-6) THEN ! if we are at a corner
	            ! get normalvector left and right from the corner
	            CALL boundary_getNormalVec2D(rboundary,ibct,parmv+1e-6,normalV(1),normalV(2))
	            CALL boundary_getNormalVec2D(rboundary,ibct,parmv-1e-6,tangentialv(1),tangentialv(2))
	            ! add those normalvectors
	            normalV = normalV + tangentialv
	            ! and make it unit length
	            normalv = normalv /SQRT(normalV(1)**2.0_DP + normalV(2)**2.0_DP)
	        END IF
	        
	        ! tangential Vector
	        tangentialv(1) = -normalV(2)
	        tangentialv(2) =  normalV(1)
	        
	        ! update preconditioner
	        DO ij = p_kld(i), p_kld(i+1)-1
	            IF (ij .NE. ii) THEN
	                rarrayP(1)%Da(ij) = 0.0_dp
	                rarrayP(2)%Da(ij) = 0.0_dp
	                rarrayP(3)%Da(ij) = 0.0_dp
	            END IF
	        END DO
	                
	        
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
	        
	        IF (4.0_dp*(SQRT(hstar*gravconst))<-2.0_dp*vnproj) THEN
	            WRITE (*,*) 'ERROR: Boundary Riemann Problem - Physical drybed generation! Wetbet Riemann Solver not applicable!'
	            WRITE (*,*) hstar,vnproj
	        END IF
	        
	        ! Initialisiere hstarstar
	        ! entweder mit der predicted solution
	         hstarstar = hstar
	        ! oder mit dem two-rarefaction ansatz (der explizit ausgerechnet werden kann)
	        !hstarstar = ((SQRT(gravconst*hstar)+0.5*vnproj)**2.0_dp)/gravconst
	        	        
	        iteh: DO ite = 1, itemax
	        ! get h** as solution of the riemannian problem
	        oldh = hstarstar
	        IF (hstarstar .LE. hstar) THEN
	            ! h** can be computed explicitely
	            hstarstar = ((SQRT(gravconst*hstar)+0.5*vnproj)**2.0_dp)/gravconst
	        ELSE
	            
	            hm = hstarstar
	            
	            ! compute h** by fixpointiteration
	            hm = hstar + vnproj/SQRT(0.5*gravconst*(hm+hstar)/(hm*hstar))
	            hm = hstar + vnproj/SQRT(0.5*gravconst*(hm+hstar)/(hm*hstar))
	            hm = hstar + vnproj/SQRT(0.5*gravconst*(hm+hstar)/(hm*hstar))
	            
	            ! or by newton method
	            ! hm = hm - ...
	            
	            hstarstar = hm
	        END IF
	        ! Test if the algorithm converged
	        IF (ABS(hstarstar - oldh)< 1e-8) THEN
	            EXIT iteh
	        END IF
	        END DO iteh
	        
	        ! Test for convergence and give out warnings
	        IF (ite==50) THEN
	            WRITE (*,*) 'ERROR! Boundary condition riemann problem did not converge'
	        END IF
	        IF (hstarstar<0.05) THEN
	            WRITE (*,*) 'ERROR! Boundary condition riemann problem h very small or smaller than zero'
	            hstarstar=hstar
	        END IF
	        IF (ite>10) THEN
	            WRITE (*,*) 'WARNING! Boundary condition riemann problem convergence is slow'
	        END IF
	        
	        	        
	        ! compute u** and v** as projection of u* and v* on the tangential Vector
	        ustarstar = vtproj * tangentialV(1)
	        vstarstar = vtproj * tangentialV(2)
	        
	        
	        ! Finally transform back to primitive variables and update the solution
	        rarraySol(1)%Da(i) = hstarstar
	        rarraySol(2)%Da(i) = hstarstar*ustarstar
	        rarraySol(3)%Da(i) = hstarstar*vstarstar
	        
        
	        END IF
	        	        
	    
        END DO	
	END DO
	
    END SUBROUTINE
    
END MODULE
