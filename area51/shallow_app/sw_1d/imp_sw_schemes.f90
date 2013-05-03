module imp_sw_schemes
	
use vartypes
	
	!************** Implizite Shallow Water Schemata *************************
	
	contains
	
	subroutine imp_onedswupwind(U, U0, x, dt, dx, theta, scheme, funcF, funcA, funcR, funcRi)
	! Loest die eindim Shallow Water Equation mit implizitem Upwind Scheme
	
	!Variablendeklaration
	implicit none
	type(t_vectorBlock), intent(in)			:: U0				! Startwerte
	type(t_vectorBlock), intent(inout)		:: U				! Loesungsvektoren
	real(dp), intent(in)					:: dt, dx			! delta x und delta t
	real(dp), dimension(:), intent(in)		:: x				! Ortsgitter
	integer, intent(in)						:: scheme			! waehlt das Schema: 0: Upwind, 1: TVD
	integer									:: xdim				! Dimension des Problems
	integer									:: i				! Laufvariable
	real(dp), intent(in)					:: theta			! theta = 0:explizit, =1/2:crank nicolson, =1:implizit
	real(dp)								:: k1, k2, k3		! Hilfsvariblen: (1-theta)dt/(2dx), (theta)dt/(2dx), dt/dx
	real(dp), dimension(2)					:: deltaUrechts, deltaUlinks
	real(dp), dimension(2)					:: Frechts, Flinks
	real(dp)								:: hRoerechts, hRoelinks, uRoerechts, uRoelinks		! Die Roe-Werte
	real(dp), dimension(2,2)				:: bAdachrechts, bAdachlinks						! Betrag von Adach
	real(dp), dimension(2,2)				:: Adachrechts, Adachlinks							! Adach
	real(dp), dimension(2,2)				:: tempR, tempA, tempRi
	real(dp), dimension(2)					:: tempU1, tempB
	type(t_vectorBlock)						:: right_B			! Rechte Seite der Gleichungen
	type(t_blockMatrixTridiag)				:: koeff_A			! Koeffizientenmatrix
	real(dp), dimension(2,2)				:: C, D, E			! Hilfsmatrizen zur Berechnung von koeff_A
	real(dp), dimension(2,2)				:: Eye				! Einheitsmatrix
	integer									:: iter				! Fuer die Defekt-Korrektur
	integer, parameter						:: itemin=10,itemax=100
	real(dp), parameter						:: tol=1d-8, eps=1d-14
	type(t_vectorBlock)						:: residuum, deltaU
	type(t_vectorScalar)					:: tempvec1, tempvec2
	real(dp)								:: dUnorm, Unorm, resnorm	! Norm des Residuums + Norm von (delta)U
	real(dp), dimension(2)					:: tempdeltaW, tempdeltaU				! fuer TVD
	real(dp), dimension(2)					:: deltaWdachrechts, deltaWdachlinks
	real(dp), dimension(2)					:: deltaWrechts, deltaWlinks
	real(dp), dimension(2,2)				:: Rrechts, Rlinks, Rirechts, Rilinks
	real(dp), dimension(2,2)				:: Lambdarechts, Lambdalinks
	integer									:: j
	real(dp)								:: hRoej, uRoej
	
	
	! Interfaces angeben
	interface
		! Funktion F aus der Shallow Water PDE
		subroutine funcF(F,q1,q2)
        use vartypes
		real(dp), intent(in)							:: q1, q2
		real(dp), dimension(2), intent(out)				:: F
		end subroutine funcF
		! diagonalisierte Jacobi Matrix der Funktion F aus der Shallow Water PDE
		subroutine funcA(A,q1,q2)
		use vartypes
		real(dp), intent(in)							:: q1, q2
		real(dp), dimension(2,2), intent(out)			:: A
		end subroutine funcA
		! Transformationsmatrix R
		subroutine funcR(R,q1,q2)
		use vartypes
		real(dp), intent(in)							:: q1, q2
		real(dp), dimension(2,2), intent(out)			:: R
		end subroutine funcR
		! Transformationsmatrix R^-1
		subroutine funcRi(Ri,q1,q2)
		use vartypes
		real(dp), intent(in)							:: q1, q2
		real(dp), dimension(2,2), intent(out)			:: Ri
		end subroutine funcRi
	end interface
	
	
	! Hilfsvariablen
	k1 = (1-theta)*dt/(2.0d0*dx)									! fuer den expliziten Teil (rechte Seite)
	k2 = theta*dt/(2.0d0*dx)										! fuer den impliziten Teil (linke Seite)
	k3 = dt/dx														! fuer TVD
	Eye = reshape( (/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/) )		! Einheitsmatrix
	xdim = U0%Dvectorblock(1)%neq									! Dimension des Problems
	
	! U0 in U schreiben
	U%DvectorBlock(1)%datas = U0%DvectorBlock(1)%datas
	U%DvectorBlock(2)%datas = U0%DvectorBlock(2)%datas
	
	! Erstelle rechte Seite: 	right_B
	! und Koeffizientenmatrix:	koeff_A
	call make_vectorScalar(tempvec1,U0%Dvectorblock(1)%neq)
	call make_vectorScalar(tempvec2,U0%Dvectorblock(1)%neq)
	call make_vectorBlock(right_B,U0%nblocks,U0%Dvectorblock(1)%neq)
	call make_vectorBlock(residuum,U0%nblocks,U0%Dvectorblock(1)%neq)
	call make_vectorBlock(deltaU,U0%nblocks,U0%Dvectorblock(1)%neq)
	call make_blockMatrixTridiag(koeff_A,2*U0%nblocks,U0%Dvectorblock(1)%neq)
	
	! Berechne die rechte Seite der Gleichungen: right_B
	call create_right_side
	
	
	! testweise die Koeffizientenmatrix ausgeben
	!call create_koeff_Matrix
	!write(*,*) koeff_A%trimat(1)%upper%datas
	!write(*,*) koeff_A%trimat(1)%main%datas
	!write(*,*) koeff_A%trimat(1)%lower%datas
	!write(*,*) koeff_A%trimat(2)%upper%datas
	!write(*,*) koeff_A%trimat(2)%main%datas
	!write(*,*) koeff_A%trimat(2)%lower%datas
	
	
	
	! Jetzt loese mit Defektkorrektur
	iter=0
	defect_iter: do
		! neue Koeffizientenmatrix berechnen
		call create_koeff_Matrix
		
		! Residuum berechnen
		call create_Residuum
		
		
		! Loese das Gleichungssystem fuer deltaU
		call linSolvThomas(koeff_A%trimat(1),residuum%DvectorBlock(1),deltaU%DvectorBlock(1))
		call linSolvThomas(koeff_A%trimat(2),residuum%DvectorBlock(2),deltaU%DvectorBlock(2))
		
		! Addiere deltaU
		U%DvectorBlock(1)%datas = U%DvectorBlock(1)%datas + deltaU%DvectorBlock(1)%datas
		U%DvectorBlock(2)%datas = U%DvectorBlock(2)%datas + deltaU%DvectorBlock(2)%datas
		
		! Abbruchkriterium
		iter=iter+1
		resnorm = sum(abs(residuum%DvectorBlock(1)%datas))+ sum(abs(residuum%DvectorBlock(2)%datas))
		dUnorm  = sum(abs(deltaU%DvectorBlock(1)%datas))  + sum(abs(deltaU%DvectorBlock(2)%datas))
		Unorm   = sum(abs(U%DvectorBlock(1)%datas))       + sum(abs(U%DvectorBlock(2)%datas))
		if ((iter>itemin).and.( (resnorm<tol).or.(dUnorm/Unorm<eps) )) then
			exit defect_iter
		else
			if (iter>itemax) then
			write(*,*) 'Defekt-Korrektur: Warnung, keine Konvergenz!'
			exit defect_iter
			end if
		end if
	end do defect_iter
	
	! gebe Anzahl der Iterationen aus
	!write(*,*) iter
	
	! Gib Speicher frei
	call unmake_vectorScalar(tempvec1)
	call unmake_vectorScalar(tempvec2)
	call unmake_vectorBlock(right_B)			! rechte Seite
	call unmake_vectorBlock(residuum)
	call unmake_vectorBlock(deltaU)
	call unmake_blockMatrixTridiag(koeff_A)		! Koeffizientenmatrix
	
	
	
	contains
	
	
	
	
	subroutine create_right_side		! Berechnet die rechte Seite und speichert sie in right_B
		implicit none
		integer					:: i	! Laufvariable
		
		if (scheme==0) then		! Upwind
		
		! Eintraege fuer die Randwerte
		right_B%Dvectorblock(1)%datas(1) = U0%Dvectorblock(1)%datas(1)
		right_B%Dvectorblock(2)%datas(1) = U0%Dvectorblock(2)%datas(1)
		right_B%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq) = &
					    U0%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq)
		right_B%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq) = &
						U0%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq)
						
	
		! Die inneren Eintraege der rechten Seite berechnen
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
			bAdachrechts = matmul(tempR,matmul(abs(tempA),tempRi))
		
			call funcR(tempR,hRoelinks,hRoelinks*uRoelinks)
			call funcA(tempA,hRoelinks,hRoelinks*uRoelinks)
			call funcRi(tempRi,hRoelinks,hRoelinks*uRoelinks)
			bAdachlinks = matmul(tempR,matmul(abs(tempA),tempRi))
		
			call funcF(Frechts,U0%Dvectorblock(1)%datas(i+1),U0%Dvectorblock(2)%datas(i+1))
			call funcF(Flinks,U0%Dvectorblock(1)%datas(i-1),U0%Dvectorblock(2)%datas(i-1))
		
			tempU1(1) = U0%Dvectorblock(1)%datas(i)
			tempU1(2) = U0%Dvectorblock(2)%datas(i)
		
		
			! hier werden die Vektoren berechnet
			tempB = tempU1 - k1*(Frechts-Flinks-matmul(bAdachrechts,deltaUrechts)+ &
					matmul(bAdachlinks,deltaUlinks))
		
			! und hier abgespeichert
			right_B%Dvectorblock(1)%datas(i) = tempB(1)
			right_B%Dvectorblock(2)%datas(i) = tempB(2)
	
		end do
		! rechte Seite fertig: right_B

		
		else		! TVD
			
		! Eintraege fuer die Randwerte
		right_B%Dvectorblock(1)%datas(1) = U0%Dvectorblock(1)%datas(1)
		right_B%Dvectorblock(1)%datas(2) = U0%Dvectorblock(1)%datas(2)
		right_B%Dvectorblock(2)%datas(1) = U0%Dvectorblock(2)%datas(1)
		right_B%Dvectorblock(2)%datas(2) = U0%Dvectorblock(2)%datas(2)
		
		right_B%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq) = &
					    U0%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq)
		right_B%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq-1) = &
					    U0%Dvectorblock(1)%datas(U0%Dvectorblock(1)%neq-1)
		right_B%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq) = &
						U0%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq)
		right_B%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq-1) = &
						U0%Dvectorblock(2)%datas(U0%Dvectorblock(2)%neq-1)
	
	
		! Loesungsschleife um alle inneren Werte zu berechnen
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
			bAdachrechts = matmul(Rrechts,matmul(abs(Lambdarechts),Rirechts))
		
			call funcR(Rlinks,hRoelinks,hRoelinks*uRoelinks)
			call funcA(Lambdalinks,hRoelinks,hRoelinks*uRoelinks)
			call funcRi(Rilinks,hRoelinks,hRoelinks*uRoelinks)
			bAdachlinks = matmul(Rlinks,matmul(abs(Lambdalinks),Rilinks))
		
			call funcF(Frechts,U0%Dvectorblock(1)%datas(i+1),U0%Dvectorblock(2)%datas(i+1))
			call funcF(Flinks,U0%Dvectorblock(1)%datas(i-1),U0%Dvectorblock(2)%datas(i-1))
		
			! fuer tvd
		
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
		
		
			! rechte Seite berechnen
			tempB = tempU1 - k1*(Frechts-Flinks-matmul(bAdachrechts,deltaUrechts)+ &									! Der Roe-Upwind Teil
					matmul(bAdachlinks,deltaUlinks)+&
					!										  -k3*matmul(Lambdarechts,Lambdarechts)   wuerde noch bei expl. TVD da stehen
					( matmul(matmul(Rrechts,(abs(Lambdarechts)                                      )),deltaWdachrechts)&		! der zusaetzliche LW-Teil
					 -matmul(matmul(Rlinks ,(abs(Lambdalinks )                                      )),deltaWdachlinks) ) )	! inklusive Limiter in deltaWdach
					 !										  -k3*matmul(Lambdalinks ,Lambdalinks )   wuerde noch bei expl. TVD da stehen
		
			! rechte Seite abspeichern
			right_B%Dvectorblock(1)%datas(i) = tempB(1)
			right_B%Dvectorblock(2)%datas(i) = tempB(2)
		
	end do
		end if
	end subroutine create_right_side
	
	
	
	subroutine create_koeff_Matrix				! erstellt die aktuelle Koeffizientenmatrix und speichert sie in koeff_A
	implicit none
	integer					:: i	! Laufvariable
	
	! Zuerst schreibe die Eintraege fuer die Randwerte
	! fuer die Hauptdiagonalmatrizen
	koeff_A%trimat(1)%main%datas(1)=1.0d0
	koeff_A%trimat(2)%main%datas(1)=1.0d0
	koeff_A%trimat(1)%main%datas(xdim)=1.0d0
	koeff_A%trimat(2)%main%datas(xdim)=1.0d0
	koeff_A%trimat(1)%upper%datas(1)=0.0d0
	koeff_A%trimat(2)%upper%datas(1)=0.0d0
	koeff_A%trimat(1)%lower%datas(xdim-1)=0.0d0
	koeff_A%trimat(2)%lower%datas(xdim-1)=0.0d0
	! fuer die Nebendiagonalmatrizen
	koeff_A%trimat(3)%main%datas(1)=0.0d0
	koeff_A%trimat(4)%main%datas(1)=0.0d0
	koeff_A%trimat(3)%main%datas(xdim)=0.0d0
	koeff_A%trimat(4)%main%datas(xdim)=0.0d0
	koeff_A%trimat(3)%upper%datas(1)=0.0d0
	koeff_A%trimat(4)%upper%datas(1)=0.0d0
	koeff_A%trimat(3)%lower%datas(xdim-1)=0.0d0
	koeff_A%trimat(4)%lower%datas(xdim-1)=0.0d0
		
			
	! Nun berechne die inneren Eintraege
	do i=2, xdim-1
		
			! Berechne einige Hilfsvariablen
			hRoerechts = sqrt(U%Dvectorblock(1)%datas(i+1)*U%Dvectorblock(1)%datas(i))
			hRoelinks  = sqrt(U%Dvectorblock(1)%datas(i)*U%Dvectorblock(1)%datas(i-1))
			uRoerechts = (1/sqrt(U%Dvectorblock(1)%datas(i))*U%Dvectorblock(2)%datas(i)+ &
						1/sqrt(U%Dvectorblock(1)%datas(i+1))*U%Dvectorblock(2)%datas(i+1))/ &
						(sqrt(U%Dvectorblock(1)%datas(i))+sqrt(U%Dvectorblock(1)%datas(i+1)))
			uRoelinks  = (1/sqrt(U%Dvectorblock(1)%datas(i-1))*U%Dvectorblock(2)%datas(i-1)+ &
						1/sqrt(U%Dvectorblock(1)%datas(i))*U%Dvectorblock(2)%datas(i))/ &
						(sqrt(U%Dvectorblock(1)%datas(i-1))+sqrt(U%Dvectorblock(1)%datas(i)))
		
			call funcR(tempR,hRoerechts,hRoerechts*uRoerechts)
			call funcA(tempA,hRoerechts,hRoerechts*uRoerechts)
			call funcRi(tempRi,hRoerechts,hRoerechts*uRoerechts)
			Adachrechts  = matmul(tempR,matmul(tempA,tempRi))
			bAdachrechts = matmul(tempR,matmul(abs(tempA),tempRi))
		
			call funcR(tempR,hRoelinks,hRoelinks*uRoelinks)
			call funcA(tempA,hRoelinks,hRoelinks*uRoelinks)
			call funcRi(tempRi,hRoelinks,hRoelinks*uRoelinks)
			Adachlinks  = matmul(tempR,matmul(tempA,tempRi))
			bAdachlinks = matmul(tempR,matmul(abs(tempA),tempRi))
		
			
			! hier werden die lokalen Hilfsmatrizen C,D,E berechnet, die spaeter in koeff_A geschrieben werden
			C = k2*(-Adachlinks-bAdachlinks)
			D = Eye + k2*(-Adachrechts+Adachlinks+bAdachrechts+bAdachlinks)
			E = k2*(Adachrechts-bAdachrechts)
			
		
			! hier werden die Hilfsmatrizen in koeff_A geschrieben
			! in die Matrizen auf der Hauptdiagonalen
			koeff_A%trimat(1)%main%datas(i)    = D(1,1)			! Hauptdiagonale
			koeff_A%trimat(2)%main%datas(i)    = D(2,2)
			koeff_A%trimat(1)%upper%datas(i)   = E(1,1)			! obere Diagonale
			koeff_A%trimat(2)%upper%datas(i)   = E(2,2)
			koeff_A%trimat(1)%lower%datas(i-1) = C(1,1)			! untere Diagonale
			koeff_A%trimat(2)%lower%datas(i-1) = C(2,2)
			! in die Matrizen auf der Nebendiagonalen
			koeff_A%trimat(3)%main%datas(i)    = D(1,2)			! Hauptdiagonale
			koeff_A%trimat(4)%main%datas(i)    = D(2,1)
			koeff_A%trimat(3)%upper%datas(i)   = E(1,2)			! obere Diagonale
			koeff_A%trimat(4)%upper%datas(i)   = E(2,1)
			koeff_A%trimat(3)%lower%datas(i-1) = C(1,2)			! untere Diagonale
			koeff_A%trimat(4)%lower%datas(i-1) = C(2,1)
			! die Matrix-Bloecke sind jetzt wie folgt angeordnet:	koeff_A%trimat(1) koeff_A%trimat(3)
			!														koeff_A%trimat(4) koeff_A%trimat(2)
	
		end do
		! Koeffizientenmatrix fertig: koeff_A
	end subroutine create_koeff_Matrix
	
	
	
	! Limiter Funktion fuer TVD
	real(dp) function limiterfunc2(z,n)
    	implicit none
		real(dp), intent(in)		:: z, n				! Zaehler und Nenner des Slope Ratios
		integer, parameter			:: limiter = 3		! Wahl des Limiters (1 = Minmod, 4 = Superbee)
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
    
    
    
    ! Residuum berechnen
    subroutine create_Residuum
    implicit none
    integer					:: i
    
	if (scheme==0) then				! Upwind-Residuum
    	call MatVecMul(koeff_A%trimat(1),U%DvectorBlock(1),tempvec1)
		call MatVecMul(koeff_A%trimat(3),U%DvectorBlock(2),tempvec2)
		residuum%DvectorBlock(1)%datas = right_B%DvectorBlock(1)%datas - tempvec1%datas - tempvec2%datas
		call MatVecMul(koeff_A%trimat(2),U%DvectorBlock(2),tempvec1)
		call MatVecMul(koeff_A%trimat(4),U%DvectorBlock(1),tempvec2)
		residuum%DvectorBlock(2)%datas = right_B%DvectorBlock(2)%datas - tempvec1%datas - tempvec2%datas
		
	else							! TVD-Residuum
		! Zuerst uebernehme die Eintraege fuer die Randwerte
		residuum%DvectorBlock(1)%datas(1) = 0
		residuum%DvectorBlock(2)%datas(1) = 0
		residuum%DvectorBlock(1)%datas(2) = 0
		residuum%DvectorBlock(2)%datas(2) = 0
		residuum%DvectorBlock(1)%datas(xdim-1) = 0
		residuum%DvectorBlock(2)%datas(xdim-1) = 0
		residuum%DvectorBlock(1)%datas(xdim) = 0
		residuum%DvectorBlock(2)%datas(xdim) = 0
		
		
		! Schleife um alle inneren Werte zu berechnen
		do i=3, U%Dvectorblock(1)%neq-2
		
			! Berechne einige Hilfsvariablen
			deltaUrechts(1) = U0%Dvectorblock(1)%datas(i+1) - U%Dvectorblock(1)%datas(i)
			deltaUrechts(2) = U0%Dvectorblock(2)%datas(i+1) - U%Dvectorblock(2)%datas(i)
			deltaUlinks(1)  = U0%Dvectorblock(1)%datas(i) - U%Dvectorblock(1)%datas(i-1)
			deltaUlinks(2)  = U0%Dvectorblock(2)%datas(i) - U%Dvectorblock(2)%datas(i-1)
		
			hRoerechts = sqrt(U%Dvectorblock(1)%datas(i+1)*U%Dvectorblock(1)%datas(i))
			hRoelinks  = sqrt(U%Dvectorblock(1)%datas(i)*U%Dvectorblock(1)%datas(i-1))
			uRoerechts = (1/sqrt(U%Dvectorblock(1)%datas(i))*U%Dvectorblock(2)%datas(i)+ &
						1/sqrt(U%Dvectorblock(1)%datas(i+1))*U%Dvectorblock(2)%datas(i+1))/ &
						(sqrt(U%Dvectorblock(1)%datas(i))+sqrt(U%Dvectorblock(1)%datas(i+1)))
			uRoelinks  = (1/sqrt(U%Dvectorblock(1)%datas(i-1))*U%Dvectorblock(2)%datas(i-1)+ &
						1/sqrt(U%Dvectorblock(1)%datas(i))*U%Dvectorblock(2)%datas(i))/ &
						(sqrt(U%Dvectorblock(1)%datas(i-1))+sqrt(U%Dvectorblock(1)%datas(i)))
		
			call funcR(Rrechts,hRoerechts,hRoerechts*uRoerechts)
			call funcA(Lambdarechts,hRoerechts,hRoerechts*uRoerechts)
			call funcRi(Rirechts,hRoerechts,hRoerechts*uRoerechts)
			bAdachrechts = matmul(Rrechts,matmul(abs(Lambdarechts),Rirechts))
		
			call funcR(Rlinks,hRoelinks,hRoelinks*uRoelinks)
			call funcA(Lambdalinks,hRoelinks,hRoelinks*uRoelinks)
			call funcRi(Rilinks,hRoelinks,hRoelinks*uRoelinks)
			bAdachlinks = matmul(Rlinks,matmul(abs(Lambdalinks),Rilinks))
		
			call funcF(Frechts,U%Dvectorblock(1)%datas(i+1),U%Dvectorblock(2)%datas(i+1))
			call funcF(Flinks,U%Dvectorblock(1)%datas(i-1),U%Dvectorblock(2)%datas(i-1))
		
			! fuer tvd
		
			deltaWrechts = matmul(Rirechts,deltaUrechts)
			deltaWlinks  = matmul(Rilinks ,deltaUlinks )
		
			! berechne die rechten slope ratios und damit deltaWdach rechts
			j=i-nint(sign(1d0,Lambdarechts(1,1)))
			hRoej = sqrt(U%Dvectorblock(1)%datas(j+1)*U%Dvectorblock(1)%datas(j))
			uRoej = (1/sqrt(U%Dvectorblock(1)%datas(j))*U%Dvectorblock(2)%datas(j)+ &
						1/sqrt(U%Dvectorblock(1)%datas(j+1))*U%Dvectorblock(2)%datas(j+1))/ &
						(sqrt(U%Dvectorblock(1)%datas(j))+sqrt(U%Dvectorblock(1)%datas(j+1)))
			call funcRi(tempRi,hRoej,hRoej*uRoej)
			tempdeltaU(1) = U%Dvectorblock(1)%datas(j+1) - U%Dvectorblock(1)%datas(j)
			tempdeltaU(2) = U%Dvectorblock(2)%datas(j+1) - U%Dvectorblock(2)%datas(j)
			tempdeltaW = matmul(tempRi,tempdeltaU)
			!deltaWdachrechts(1)=deltaWrechts(1)*limiterfunc(tempdeltaW(1),deltaWrechts(1))
			deltaWdachrechts(1)=limiterfunc2(tempdeltaW(1),deltaWrechts(1))
		
			j=i-nint(sign(1d0,Lambdarechts(2,2)))
			hRoej = sqrt(U%Dvectorblock(1)%datas(j+1)*U%Dvectorblock(1)%datas(j))
			uRoej = (1/sqrt(U%Dvectorblock(1)%datas(j))*U%Dvectorblock(2)%datas(j)+ &
						1/sqrt(U%Dvectorblock(1)%datas(j+1))*U%Dvectorblock(2)%datas(j+1))/ &
						(sqrt(U%Dvectorblock(1)%datas(j))+sqrt(U%Dvectorblock(1)%datas(j+1)))
			call funcRi(tempRi,hRoej,hRoej*uRoej)
			tempdeltaU(1) = U%Dvectorblock(1)%datas(j+1) - U%Dvectorblock(1)%datas(j)
			tempdeltaU(2) = U%Dvectorblock(2)%datas(j+1) - U%Dvectorblock(2)%datas(j)
			tempdeltaW = matmul(tempRi,tempdeltaU)
			!deltaWdachrechts(2)=deltaWrechts(2)*limiterfunc(tempdeltaW(2),deltaWrechts(2))
			deltaWdachrechts(2)=limiterfunc2(tempdeltaW(2),deltaWrechts(2))
		
			! berechne die linken slope ratios und damit deltaWdach links
			j=i-1-nint(sign(1d0,Lambdalinks(1,1)))
			hRoej = sqrt(U%Dvectorblock(1)%datas(j+1)*U%Dvectorblock(1)%datas(j))
			uRoej = (1/sqrt(U%Dvectorblock(1)%datas(j))*U%Dvectorblock(2)%datas(j)+ &
						1/sqrt(U%Dvectorblock(1)%datas(j+1))*U%Dvectorblock(2)%datas(j+1))/ &
						(sqrt(U%Dvectorblock(1)%datas(j))+sqrt(U%Dvectorblock(1)%datas(j+1)))
			call funcRi(tempRi,hRoej,hRoej*uRoej)
			tempdeltaU(1) = U%Dvectorblock(1)%datas(j+1) - U%Dvectorblock(1)%datas(j)
			tempdeltaU(2) = U%Dvectorblock(2)%datas(j+1) - U%Dvectorblock(2)%datas(j)
			tempdeltaW = matmul(tempRi,tempdeltaU)
			!deltaWdachlinks(1)=deltaWlinks(1)*limiterfunc(tempdeltaw(1),deltaWlinks(1))
			deltaWdachlinks(1)=limiterfunc2(tempdeltaw(1),deltaWlinks(1))
			
			j=i-1-nint(sign(1d0,Lambdalinks(2,2)))
			hRoej = sqrt(U%Dvectorblock(1)%datas(j+1)*U%Dvectorblock(1)%datas(j))
			uRoej = (1/sqrt(U%Dvectorblock(1)%datas(j))*U%Dvectorblock(2)%datas(j)+ &
						1/sqrt(U%Dvectorblock(1)%datas(j+1))*U%Dvectorblock(2)%datas(j+1))/ &
						(sqrt(U%Dvectorblock(1)%datas(j))+sqrt(U%Dvectorblock(1)%datas(j+1)))
			call funcRi(tempRi,hRoej,hRoej*uRoej)
			tempdeltaU(1) = U%Dvectorblock(1)%datas(j+1) - U%Dvectorblock(1)%datas(j)
			tempdeltaU(2) = U%Dvectorblock(2)%datas(j+1) - U%Dvectorblock(2)%datas(j)
			tempdeltaW = matmul(tempRi,tempdeltaU)
			!deltaWdachlinks(2)=deltaWlinks(2)*limiterfunc(tempdeltaw(2),deltaWlinks(2))
			deltaWdachlinks(2)=limiterfunc2(tempdeltaw(2),deltaWlinks(2))
		
		
		
			tempU1(1) = U%Dvectorblock(1)%datas(i)
			tempU1(2) = U%Dvectorblock(2)%datas(i)
		
		
			! berechne, was von der rechten Seite abgezogen werden muss
			tempB = tempU1 + k2*(Frechts-Flinks-matmul(bAdachrechts,deltaUrechts)+ &									! Der Roe-Upwind Teil
					matmul(bAdachlinks,deltaUlinks)+&
					!										  -k3*matmul(Lambdarechts,Lambdarechts)   wuerde noch bei expl. TVD da stehen
					( matmul(matmul(Rrechts,(abs(Lambdarechts)                                      )),deltaWdachrechts)&		! der zusaetzliche LW-Teil
					 -matmul(matmul(Rlinks ,(abs(Lambdalinks )                                      )),deltaWdachlinks) ) )	! inklusive Limiter in deltaWdach
					 !										  -k3*matmul(Lambdalinks ,Lambdalinks )   wuerde noch bei expl. TVD da stehen
		
			! Residuum berechnen und abspeichern
			residuum%DvectorBlock(1)%datas(i) = right_B%DvectorBlock(1)%datas(i)-tempB(1)
			residuum%DvectorBlock(2)%datas(i) = right_B%DvectorBlock(2)%datas(i)-tempB(2)
		
	end do
	
	end if
    end subroutine create_Residuum
	
	
		
	end subroutine imp_onedswupwind
	
end module imp_sw_schemes
