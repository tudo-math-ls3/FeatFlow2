module exp_sw_schemes
! Schemata zum Lösen der Shallow Water Equations
	
	! benutzte Module hinzufügen
	use vartypes

	! Implizite Deklaration ausschalten
	implicit none
	
	contains
	
	subroutine onedswupwind(U, U0, x, dt, dx, funcF, funcA, funcR, funcRi)
	! Löst die eindim Shallow Water Equation mit Upwind Scheme
	
	!Variablendeklaration
	implicit none
	type(t_vectorBlock), intent(inout)		:: U0				! Lösungsvektoren (Achtung, man kann nur die Daten verändern, nicht nblocks)
	type(t_vectorBlock), intent(inout)		:: U				! Lösungsvektoren
	real(dp), intent(in)					:: dt, dx
	real(dp), dimension(:), intent(in)		:: x
	integer									:: i				! Laufvariable
	real(dp)								:: k1				! Hilfsvarible: dt/(2dx)
	real(dp), dimension(2)					:: deltaUrechts, deltaUlinks
	real(dp), dimension(2)					:: Frechts, Flinks
	real(dp)								:: hRoerechts, hRoelinks, uRoerechts, uRoelinks
	real(dp), dimension(2,2)				:: Adachrechts, Adachlinks
	real(dp), dimension(2,2)				:: tempR, tempA, tempRi
	real(dp), dimension(2)					:: tempU1, tempU2
	
	! Interfaces angeben
	interface
		! Funktion F aus der Shallow Water PDE
		subroutine funcF(F,q1,q2)
        use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2), intent(out)		:: F
		end subroutine funcF
		! diagonalisierte Jacobi Matrix der Funktion F aus der Shallow Water PDE
		subroutine funcA(A,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2,2), intent(out)	:: A
		end subroutine funcA
		! Transformationsmatrix R
		subroutine funcR(R,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2,2), intent(out)	:: R
		end subroutine funcR
		! Transformationsmatrix R^-1
		subroutine funcRi(Ri,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2,2), intent(out)	:: Ri
		end subroutine funcRi
	end interface
	
	! Hilfsvariable
	k1 = dt/(2.0d0*dx)
	
	! Randwerte übernehmen
	U%Dvectorblock(1)%datas(1) = U0%Dvectorblock(1)%datas(1)
	U%Dvectorblock(2)%datas(1) = U0%Dvectorblock(2)%datas(1)
	U%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq) = &
				    U0%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq)
	U%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq) = &
					U0%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq)
	
	
	! Lösungsschleife um alle inneren Werte zu berechnen
	do i=2, U0%Dvectorblock(1)%neq-1
		
		! Berechne einige Hilfsvariablen
		deltaUrechts(1) = U0%Dvectorblock(1)%datas(i+1) - U0%Dvectorblock(1)%datas(i)
		deltaUrechts(2) = U0%Dvectorblock(2)%datas(i+1) - U0%Dvectorblock(2)%datas(i)
		deltaUlinks(1)  = U0%Dvectorblock(1)%datas(i) - U0%Dvectorblock(1)%datas(i-1)
		deltaUlinks(2)  = U0%Dvectorblock(2)%datas(i) - U0%Dvectorblock(2)%datas(i-1)
		
		hRoerechts = sqrt(U0%Dvectorblock(1)%datas(i+1)*U0%Dvectorblock(1)%datas(i))
		hRoelinks  = sqrt(U0%Dvectorblock(1)%datas(i)*U0%Dvectorblock(1)%datas(i-1))
		uRoerechts = (1/sqrt(U0%Dvectorblock(1)%datas(i))*U0%Dvectorblock(2)%datas(i)+ &
					1/sqrt(U0%Dvectorblock(1)%datas(i+1))*U0%Dvectorblock(2)%datas(i+1))/ &
					(sqrt(U0%Dvectorblock(1)%datas(i))+sqrt(U0%Dvectorblock(1)%datas(i+1)))
		uRoelinks  = (1/sqrt(U0%Dvectorblock(1)%datas(i-1))*U0%Dvectorblock(2)%datas(i-1)+ &
					1/sqrt(U0%Dvectorblock(1)%datas(i))*U0%Dvectorblock(2)%datas(i))/ &
					(sqrt(U0%Dvectorblock(1)%datas(i-1))+sqrt(U0%Dvectorblock(1)%datas(i)))
		
		call funcR(tempR,hRoerechts,hRoerechts*uRoerechts)
		call funcA(tempA,hRoerechts,hRoerechts*uRoerechts)
		call funcRi(tempRi,hRoerechts,hRoerechts*uRoerechts)
		Adachrechts = matmul(tempR,matmul(abs(tempA),tempRi))
		
		call funcR(tempR,hRoelinks,hRoelinks*uRoelinks)
		call funcA(tempA,hRoelinks,hRoelinks*uRoelinks)
		call funcRi(tempRi,hRoelinks,hRoelinks*uRoelinks)
		Adachlinks = matmul(tempR,matmul(abs(tempA),tempRi))
		
		call funcF(Frechts,U0%Dvectorblock(1)%datas(i+1),U0%Dvectorblock(2)%datas(i+1))
		call funcF(Flinks,U0%Dvectorblock(1)%datas(i-1),U0%Dvectorblock(2)%datas(i-1))
		
		tempU1(1) = U0%Dvectorblock(1)%datas(i)
		tempU1(2) = U0%Dvectorblock(2)%datas(i)
		
		
		! Lösung berechnen
		tempU2 = tempU1 - k1*(Frechts-Flinks-matmul(Adachrechts,deltaUrechts)+ &
				matmul(Adachlinks,deltaUlinks))
		
		! Lösung abspeichern
		U%Dvectorblock(1)%datas(i) = tempU2(1)
		U%Dvectorblock(2)%datas(i) = tempU2(2)
		
		
	end do
	
	
	
		
	end subroutine onedswupwind

	
	
	
	
	
	
	subroutine onedswlaxwendroff(U, U0, x, dt, dx, funcF, funcA, funcR, funcRi)
	! Löst die eindim Shallow Water Equation mit Lax Wendroff Scheme
	
	!Variablendeklaration
	implicit none
	type(t_vectorBlock), intent(inout)		:: U0				! Lösungsvektoren (Achtung, man kann nur die Daten verändern, nicht nblocks)
	type(t_vectorBlock), intent(inout)		:: U				! Lösungsvektoren
	real(dp), intent(in)					:: dt, dx
	real(dp), dimension(:), intent(in)		:: x
	integer									:: i				! Laufvariable
	real(dp)								:: k1				! Hilfsvarible: dt/dx
	real(dp), dimension(2)					:: deltaUrechts, deltaUlinks
	real(dp), dimension(2)					:: Frechts, Flinks
	real(dp)								:: hRoerechts, hRoelinks, uRoerechts, uRoelinks
	real(dp), dimension(2,2)				:: Adachrechts, Adachlinks
	real(dp), dimension(2,2)				:: Rrechts, Rlinks, Rirechts, Rilinks
	real(dp), dimension(2,2)				:: Lambdarechts, Lambdalinks
	real(dp), dimension(2,2)				:: tempR, tempA, tempRi
	real(dp), dimension(2)					:: tempU1, tempU2
	
	! Interfaces angeben
	interface
		! Funktion F aus der Shallow Water PDE
		subroutine funcF(F,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2), intent(out)		:: F
		end subroutine funcF
		! diagonalisierte Jacobi Matrix der Funktion F aus der Shallow Water PDE
		subroutine funcA(A,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2,2), intent(out)	:: A
		end subroutine funcA
		! Transformationsmatrix R
		subroutine funcR(R,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2,2), intent(out)	:: R
		end subroutine funcR
		! Transformationsmatrix R^-1
		subroutine funcRi(Ri,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2,2), intent(out)	:: Ri
		end subroutine funcRi
	end interface
	
	! Hilfsvariable
	k1 = dt/(dx)
	
	! Randwerte übernehmen
	U%Dvectorblock(1)%datas(1) = U0%Dvectorblock(1)%datas(1)
	U%Dvectorblock(2)%datas(1) = U0%Dvectorblock(2)%datas(1)
	U%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq) = &
				    U0%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq)
	U%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq) = &
					U0%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq)
	
	
	! Lösungsschleife um alle inneren Werte zu berechnen
	do i=2, U0%Dvectorblock(1)%neq-1
		
		! Berechne einige Hilfsvariablen
		deltaUrechts(1) = U0%Dvectorblock(1)%datas(i+1) - U0%Dvectorblock(1)%datas(i)
		deltaUrechts(2) = U0%Dvectorblock(2)%datas(i+1) - U0%Dvectorblock(2)%datas(i)
		deltaUlinks(1)  = U0%Dvectorblock(1)%datas(i) - U0%Dvectorblock(1)%datas(i-1)
		deltaUlinks(2)  = U0%Dvectorblock(2)%datas(i) - U0%Dvectorblock(2)%datas(i-1)
		
		hRoerechts = sqrt(U0%Dvectorblock(1)%datas(i+1)*U0%Dvectorblock(1)%datas(i))
		hRoelinks  = sqrt(U0%Dvectorblock(1)%datas(i)*U0%Dvectorblock(1)%datas(i-1))
		uRoerechts = (1/sqrt(U0%Dvectorblock(1)%datas(i))*U0%Dvectorblock(2)%datas(i)+ &
					1/sqrt(U0%Dvectorblock(1)%datas(i+1))*U0%Dvectorblock(2)%datas(i+1))/ &
					(sqrt(U0%Dvectorblock(1)%datas(i))+sqrt(U0%Dvectorblock(1)%datas(i+1)))
		uRoelinks  = (1/sqrt(U0%Dvectorblock(1)%datas(i-1))*U0%Dvectorblock(2)%datas(i-1)+ &
					1/sqrt(U0%Dvectorblock(1)%datas(i))*U0%Dvectorblock(2)%datas(i))/ &
					(sqrt(U0%Dvectorblock(1)%datas(i-1))+sqrt(U0%Dvectorblock(1)%datas(i)))
		
		call funcR(Rrechts,hRoerechts,hRoerechts*uRoerechts)
		call funcA(Lambdarechts,hRoerechts,hRoerechts*uRoerechts)
		call funcRi(Rirechts,hRoerechts,hRoerechts*uRoerechts)
		Adachrechts = matmul(Rrechts,matmul(abs(Lambdarechts),Rirechts))
		
		call funcR(Rlinks,hRoelinks,hRoelinks*uRoelinks)
		call funcA(Lambdalinks,hRoelinks,hRoelinks*uRoelinks)
		call funcRi(Rilinks,hRoelinks,hRoelinks*uRoelinks)
		Adachlinks = matmul(Rlinks,matmul(abs(Lambdalinks),Rilinks))
		
		call funcF(Frechts,U0%Dvectorblock(1)%datas(i+1),U0%Dvectorblock(2)%datas(i+1))
		call funcF(Flinks,U0%Dvectorblock(1)%datas(i-1),U0%Dvectorblock(2)%datas(i-1))
		
		tempU1(1) = U0%Dvectorblock(1)%datas(i)
		tempU1(2) = U0%Dvectorblock(2)%datas(i)
		
		
		! Lösung berechnen
		tempU2 = tempU1 - k1/2.0d0*(Frechts-Flinks-matmul(Adachrechts,deltaUrechts)+ &		! Der Roe-Upwind Teil
				matmul(Adachlinks,deltaUlinks)+&
				( matmul(matmul(Rrechts,matmul((abs(Lambdarechts)-k1*matmul(Lambdarechts,Lambdarechts)),Rirechts)),deltaUrechts)& ! der zusätzliche LW-Teil
				-matmul(matmul(Rlinks,matmul((abs(Lambdalinks)-k1*matmul(Lambdalinks,Lambdalinks)),Rilinks)),deltaUlinks) ) )
		
		! Lösung abspeichern
		U%Dvectorblock(1)%datas(i) = tempU2(1)
		U%Dvectorblock(2)%datas(i) = tempU2(2)
		
		
	end do
	
	
	
		
	end subroutine onedswlaxwendroff
	
	
	
	
	
	
	
	subroutine onedswtvd(U, U0, x, dt, dx, funcF, funcA, funcR, funcRi)
	! Löst die eindim Shallow Water Equation mit TVD Scheme
	
	! Variablendeklaration
	implicit none
	type(t_vectorBlock), intent(inout)		:: U0				! Lösungsvektoren (Achtung, man kann nur die Daten verändern, nicht nblocks)
	type(t_vectorBlock), intent(inout)		:: U				! Lösungsvektoren
	real(dp), intent(in)					:: dt, dx
	real(dp), dimension(:), intent(in)		:: x
	integer									:: i,j				! Laufvariable
	real(dp)								:: k1				! Hilfsvarible: dt/dx
	real(dp), dimension(2)					:: deltaUrechts, deltaUlinks
	real(dp), dimension(2)					:: Frechts, Flinks
	real(dp)								:: hRoerechts, hRoelinks, uRoerechts, uRoelinks
	real(dp), dimension(2,2)				:: Adachrechts, Adachlinks
	real(dp), dimension(2,2)				:: Rrechts, Rlinks, Rirechts, Rilinks
	real(dp), dimension(2,2)				:: Lambdarechts, Lambdalinks
	real(dp), dimension(2,2)				:: tempRi
	real(dp)								:: hRoej, uRoej
	real(dp), dimension(2)					:: tempU1, tempU2
	real(dp), dimension(2)					:: tempdeltaW, tempdeltaU
	real(dp), dimension(2)					:: deltaWdachrechts, deltaWdachlinks
	real(dp), dimension(2)					:: deltaWrechts, deltaWlinks
		
	
	! Interfaces angeben
	interface
		! Funktion F aus der Shallow Water PDE
		subroutine funcF(F,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2), intent(out)		:: F
		end subroutine funcF
		! diagonalisierte Jacobi Matrix der Funktion F aus der Shallow Water PDE
		subroutine funcA(A,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2,2), intent(out)	:: A
		end subroutine funcA
		! Transformationsmatrix R
		subroutine funcR(R,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2,2), intent(out)	:: R
		end subroutine funcR
		! Transformationsmatrix R^-1
		subroutine funcRi(Ri,q1,q2)
		use vartypes
		real(dp), intent(in)					:: q1, q2
		real(dp), dimension(2,2), intent(out)	:: Ri
		end subroutine funcRi
	end interface
	
	! Hilfsvariable
	k1 = dt/(dx)
	
	! Randwerte übernehmen
	U%Dvectorblock(1)%datas(1) = U0%Dvectorblock(1)%datas(1)
	U%Dvectorblock(2)%datas(1) = U0%Dvectorblock(2)%datas(1)
	U%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq) = &
				    U0%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq)
	U%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq) = &
					U0%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq)
	U%Dvectorblock(1)%datas(2) = U0%Dvectorblock(1)%datas(2)
	U%Dvectorblock(2)%datas(2) = U0%Dvectorblock(2)%datas(2)
	U%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq-1) = &
				    U0%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq-1)
	U%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq-1) = &
					U0%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq-1)
	
	
	! Lösungsschleife um alle inneren Werte zu berechnen
	do i=3, U0%Dvectorblock(1)%neq-2
		
		! Berechne einige Hilfsvariablen
		deltaUrechts(1) = U0%Dvectorblock(1)%datas(i+1) - U0%Dvectorblock(1)%datas(i)
		deltaUrechts(2) = U0%Dvectorblock(2)%datas(i+1) - U0%Dvectorblock(2)%datas(i)
		deltaUlinks(1)  = U0%Dvectorblock(1)%datas(i) - U0%Dvectorblock(1)%datas(i-1)
		deltaUlinks(2)  = U0%Dvectorblock(2)%datas(i) - U0%Dvectorblock(2)%datas(i-1)
		
		hRoerechts = sqrt(U0%Dvectorblock(1)%datas(i+1)*U0%Dvectorblock(1)%datas(i))
		hRoelinks  = sqrt(U0%Dvectorblock(1)%datas(i)*U0%Dvectorblock(1)%datas(i-1))
		uRoerechts = (1/sqrt(U0%Dvectorblock(1)%datas(i))*U0%Dvectorblock(2)%datas(i)+ &
					1/sqrt(U0%Dvectorblock(1)%datas(i+1))*U0%Dvectorblock(2)%datas(i+1))/ &
					(sqrt(U0%Dvectorblock(1)%datas(i))+sqrt(U0%Dvectorblock(1)%datas(i+1)))
		uRoelinks  = (1/sqrt(U0%Dvectorblock(1)%datas(i-1))*U0%Dvectorblock(2)%datas(i-1)+ &
					1/sqrt(U0%Dvectorblock(1)%datas(i))*U0%Dvectorblock(2)%datas(i))/ &
					(sqrt(U0%Dvectorblock(1)%datas(i-1))+sqrt(U0%Dvectorblock(1)%datas(i)))
		
		call funcR(Rrechts,hRoerechts,hRoerechts*uRoerechts)
		call funcA(Lambdarechts,hRoerechts,hRoerechts*uRoerechts)
		call funcRi(Rirechts,hRoerechts,hRoerechts*uRoerechts)
		Adachrechts = matmul(Rrechts,matmul(abs(Lambdarechts),Rirechts))
		
		call funcR(Rlinks,hRoelinks,hRoelinks*uRoelinks)
		call funcA(Lambdalinks,hRoelinks,hRoelinks*uRoelinks)
		call funcRi(Rilinks,hRoelinks,hRoelinks*uRoelinks)
		Adachlinks = matmul(Rlinks,matmul(abs(Lambdalinks),Rilinks))
		
		call funcF(Frechts,U0%Dvectorblock(1)%datas(i+1),U0%Dvectorblock(2)%datas(i+1))
		call funcF(Flinks,U0%Dvectorblock(1)%datas(i-1),U0%Dvectorblock(2)%datas(i-1))
		
		! für tvd
		
		deltaWrechts = matmul(Rirechts,deltaUrechts)
		deltaWlinks  = matmul(Rilinks ,deltaUlinks )
		
		! berechne die rechten slope ratios und damit deltaWdach rechts
		j=i-nint(sign(1d0,Lambdarechts(1,1)))
		hRoej = sqrt(U0%Dvectorblock(1)%datas(j+1)*U0%Dvectorblock(1)%datas(j))
		uRoej = (1/sqrt(U0%Dvectorblock(1)%datas(j))*U0%Dvectorblock(2)%datas(j)+ &
					1/sqrt(U0%Dvectorblock(1)%datas(j+1))*U0%Dvectorblock(2)%datas(j+1))/ &
					(sqrt(U0%Dvectorblock(1)%datas(j))+sqrt(U0%Dvectorblock(1)%datas(j+1)))
		call funcRi(tempRi,hRoej,hRoej*uRoej)
		tempdeltaU(1) = U0%Dvectorblock(1)%datas(j+1) - U0%Dvectorblock(1)%datas(j)
		tempdeltaU(2) = U0%Dvectorblock(2)%datas(j+1) - U0%Dvectorblock(2)%datas(j)
		tempdeltaW = matmul(tempRi,tempdeltaU)
		!deltaWdachrechts(1)=deltaWrechts(1)*limiterfunc(tempdeltaW(1),deltaWrechts(1))
		deltaWdachrechts(1)=limiterfunc2(tempdeltaW(1),deltaWrechts(1))
		
		j=i-nint(sign(1d0,Lambdarechts(2,2)))
		hRoej = sqrt(U0%Dvectorblock(1)%datas(j+1)*U0%Dvectorblock(1)%datas(j))
		uRoej = (1/sqrt(U0%Dvectorblock(1)%datas(j))*U0%Dvectorblock(2)%datas(j)+ &
					1/sqrt(U0%Dvectorblock(1)%datas(j+1))*U0%Dvectorblock(2)%datas(j+1))/ &
					(sqrt(U0%Dvectorblock(1)%datas(j))+sqrt(U0%Dvectorblock(1)%datas(j+1)))
		call funcRi(tempRi,hRoej,hRoej*uRoej)
		tempdeltaU(1) = U0%Dvectorblock(1)%datas(j+1) - U0%Dvectorblock(1)%datas(j)
		tempdeltaU(2) = U0%Dvectorblock(2)%datas(j+1) - U0%Dvectorblock(2)%datas(j)
		tempdeltaW = matmul(tempRi,tempdeltaU)
		!deltaWdachrechts(2)=deltaWrechts(2)*limiterfunc(tempdeltaW(2),deltaWrechts(2))
		deltaWdachrechts(2)=limiterfunc2(tempdeltaW(2),deltaWrechts(2))
		
		! berechne die linken slope ratios und damit deltaWdach links
		j=i-1-nint(sign(1d0,Lambdalinks(1,1)))
		hRoej = sqrt(U0%Dvectorblock(1)%datas(j+1)*U0%Dvectorblock(1)%datas(j))
		uRoej = (1/sqrt(U0%Dvectorblock(1)%datas(j))*U0%Dvectorblock(2)%datas(j)+ &
					1/sqrt(U0%Dvectorblock(1)%datas(j+1))*U0%Dvectorblock(2)%datas(j+1))/ &
					(sqrt(U0%Dvectorblock(1)%datas(j))+sqrt(U0%Dvectorblock(1)%datas(j+1)))
		call funcRi(tempRi,hRoej,hRoej*uRoej)
		tempdeltaU(1) = U0%Dvectorblock(1)%datas(j+1) - U0%Dvectorblock(1)%datas(j)
		tempdeltaU(2) = U0%Dvectorblock(2)%datas(j+1) - U0%Dvectorblock(2)%datas(j)
		tempdeltaW = matmul(tempRi,tempdeltaU)
		!deltaWdachlinks(1)=deltaWlinks(1)*limiterfunc(tempdeltaw(1),deltaWlinks(1))
		deltaWdachlinks(1)=limiterfunc2(tempdeltaw(1),deltaWlinks(1))
		
		j=i-1-nint(sign(1d0,Lambdalinks(2,2)))
		hRoej = sqrt(U0%Dvectorblock(1)%datas(j+1)*U0%Dvectorblock(1)%datas(j))
		uRoej = (1/sqrt(U0%Dvectorblock(1)%datas(j))*U0%Dvectorblock(2)%datas(j)+ &
					1/sqrt(U0%Dvectorblock(1)%datas(j+1))*U0%Dvectorblock(2)%datas(j+1))/ &
					(sqrt(U0%Dvectorblock(1)%datas(j))+sqrt(U0%Dvectorblock(1)%datas(j+1)))
		call funcRi(tempRi,hRoej,hRoej*uRoej)
		tempdeltaU(1) = U0%Dvectorblock(1)%datas(j+1) - U0%Dvectorblock(1)%datas(j)
		tempdeltaU(2) = U0%Dvectorblock(2)%datas(j+1) - U0%Dvectorblock(2)%datas(j)
		tempdeltaW = matmul(tempRi,tempdeltaU)
		!deltaWdachlinks(2)=deltaWlinks(2)*limiterfunc(tempdeltaw(2),deltaWlinks(2))
		deltaWdachlinks(2)=limiterfunc2(tempdeltaw(2),deltaWlinks(2))
		
		
		
		tempU1(1) = U0%Dvectorblock(1)%datas(i)
		tempU1(2) = U0%Dvectorblock(2)%datas(i)
		
		
		! Lösung berechnen
		tempU2 = tempU1 - k1/2.0d0*(Frechts-Flinks-matmul(Adachrechts,deltaUrechts)+ &									! Der Roe-Upwind Teil
				matmul(Adachlinks,deltaUlinks)+&
				( matmul(matmul(Rrechts,(abs(Lambdarechts)-k1*matmul(Lambdarechts,Lambdarechts))),deltaWdachrechts)&	! der zusätzliche LW-Teil
				 -matmul(matmul(Rlinks ,(abs(Lambdalinks )-k1*matmul(Lambdalinks ,Lambdalinks ))),deltaWdachlinks) ) )	! inklusive Limiter in deltaWdach
		
		! Lösung abspeichern
		U%Dvectorblock(1)%datas(i) = tempU2(1)
		U%Dvectorblock(2)%datas(i) = tempU2(2)
		
		
	end do
	
	contains
	
	! Limiter Funktion
	real(dp) function limiterfunc(z,n)
		implicit none
		real(dp), intent(in)		:: z, n				! Zähler und Nenner des Slope Ratios
		integer, parameter			:: limiter = 4		! Wahl des Limiters (1 = Minmod, 4 = Superbee)
		real(dp)					:: r				! Slope Ratio

		if (abs(n)<1e-8) then
			limiterfunc = 1d0
		else
			r=z/n
        select case(limiter)
            case(1)
                    limiterfunc = max(0d0,min(1d0,r))
            case(2)
                    limiterfunc = (r+abs(r))/(1d0+abs(r))
            case(3)
                    limiterfunc = max(0d0, min(2d0,(1d0+r)/2d0,2d0*r))
            case(4)
                    limiterfunc = max(0d0,min(1d0,2d0*r),min(2d0,r))
           
         end select
         end if
    end function limiterfunc
	
	
    real(dp) function limiterfunc2(z,n)
    	implicit none
		real(dp), intent(in)		:: z, n				! Zähler und Nenner des Slope Ratios
		integer, parameter			:: limiter = 4		! Wahl des Limiters (1 = Minmod, 4 = Superbee)
		real(dp)					:: h1				! Hilfsvariable
		
		h1 = (sign(1.0d0,z)+sign(1.0d0,n))/2.0d0
		
		select case(limiter)
			case(1)
					limiterfunc2 = h1*min(abs(z),abs(n))									! MinMod
			case(2)
					if (abs(z)+abs(n)>0.001d0) then											! Van Leer
						limiterfunc2 = h1*2.0d0*abs(z*n)/(abs(n)+abs(z))
					else
						limiterfunc2 = 0d0
					end if
			case(3)
					limiterfunc2=h1*min(2.0d0*abs(z),(abs(z)+abs(n))/2.0d0,2.0d0*abs(n))	! MC
			case(4)
					limiterfunc2=h1*max(min(2.0d0*abs(z),abs(n)),min(abs(z),2.0d0*abs(n)))	!Superbee
		end select
    end function limiterfunc2
	
		
	end subroutine onedswtvd
	
	

	
	
end module exp_sw_schemes