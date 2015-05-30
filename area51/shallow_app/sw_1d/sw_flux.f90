module sw_flux

        implicit none

	contains
	
	
	! Funktion F aus der Shallow Water PDE
	subroutine swF(F,q1,q2)
	
	double precision, intent(in)					:: q1, q2
	double precision, dimension(2), intent(out)		:: F
	double precision, parameter						:: g = 1d0
	
	F(1) = q2
	F(2) = q2*q2/q1+g*q1*q1/2.0d0
	
	end subroutine swF
	
	! diagonalisierte Jacobi Matrix der Funktion F aus der Shallow Water PDE
	subroutine swA(A,q1,q2)
	
	double precision, intent(in)					:: q1, q2
	double precision, dimension(2,2), intent(out)	:: A
	double precision, parameter						:: g = 1d0
	double precision								:: c
	
	c = sqrt(g*q1)
	
	A(1,1) = q2/q1-c
	A(1,2) = 0d0
	A(2,1) = 0d0
	A(2,2) = q2/q1+c
	
	
	end subroutine swA
	
	! Transformationsmatrix R
	subroutine swR(R,q1,q2)
	
	double precision, intent(in)					:: q1, q2
	double precision, dimension(2,2), intent(out)	:: R
	double precision, parameter						:: g = 1d0
	double precision								:: c
	
	c = sqrt(g*q1)
	
	R(1,1) = 1d0
	R(1,2) = 1d0
	R(2,1) = q2/q1-c
	R(2,2) = q2/q1+c
	
	end subroutine swR
	
	! Transformationsmatrix R^-1
	subroutine swRi(Ri,q1,q2)
	
	double precision, intent(in)					:: q1, q2
	double precision, dimension(2,2), intent(out)	:: Ri
	double precision, parameter						:: g = 1d0
	double precision								:: c, twoc
	
	c = sqrt(g*q1)
	twoc = 2d0*c
	
	Ri(1,1) = (q2/q1+c)/twoc
	Ri(1,2) = -1d0/twoc
	Ri(2,1) = (c-q2/q1)/twoc
	Ri(2,2) = 1d0/twoc
	
	end subroutine swRi


end module sw_flux
