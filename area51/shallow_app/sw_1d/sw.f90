! ********** HAUPTPROGRAMM **********

program shallow_water_prog

	! benutzte Module hinzufuegen
	use vartypes
	use exp_sw_schemes
	use imp_sw_schemes
	use sw_flux


	! Implizite Deklaration ausschalten
	implicit none

	! Variablendeklaration
	integer, parameter			:: nx = 10001			! Anzahl der x-Gitterpunkte
	integer						:: animate			! =1: alle Zeitschritte schreiben, =0: nur den letzten
	real(dp)					:: rightx			! rechtes x
	integer, parameter			:: numsolv = 2		! Anzahl der Loesungsvariablen
	real(dp)					:: tfinal			! Endzeit
	real(dp)					:: time				! Fuer die Zeitschleife
	type(t_vectorBlock)			:: U, U0			! Loesungsvektoren
	real(dp), dimension(nx)		:: x				! Ortsvektor
	real(dp)					:: dx, dt			! Schrittweiten
	real(dp)					:: time0, time1		! Zur Zeitmessung
	real(dp)                    :: rcoeff
	integer						:: i,j				! Laufvariablen
	type(t_vectorBlock)			:: tmpU				! Zum Tauschen von U und U0
	
	! Speicher reservieren fuer U und U0
	call make_vectorBlock(U, numsolv,nx)
	call make_vectorBlock(U0,numsolv,nx)
	
	
	! Parameter initialisieren
	rightx = 10d0						! rechtes x
	dx = 1.0d0/(nx-1)*rightx			! Ortsschrittweite
	dt = 0.0001d0							! Zeitschrittweite
	tfinal = 1d0						! Endzeit
	animate = 0							! Auch Zwischenschritte animieren
	
	
	! Ortsgitter initialisieren
	forall(i=1:nx)
		x(i)=(i-1)*dx
	end forall
	
	
	! Loesungsvektor initialisieren
	! Hoehe h initialisieren
	U0%DvectorBlock(1)%datas = 1.0d0
	where (abs(x)<2.5)
		!U0%DvectorBlock(1)%datas = 1.d0+0.4*exp(-5*(x-5.0d0)**2.0d0)
		U0%DvectorBlock(1)%datas = 2.0d0
	end where
	! Geschwindigkeit u initialisieren
	U0%DvectorBlock(2)%datas = 0
	where (x>25)
		U0%DvectorBlock(2)%datas = 0.0d0
	end where
	
	
	! Dateien oeffnen
	open(10,file='out.dat')
	open(11,file='parameter.dat')
	write(11,*) nx
	close(11)
	
	
	! Ausgabe des Startvektors
	if (animate==1) then
        do i=1,nx
           write(10,FMT='(E15.6E3,1X,E15.6E3,1X,E15.6E3)') x(i),&
                U0%Dvectorblock(1)%datas(i),&
                U0%Dvectorblock(2)%datas(i)    ! x, h, hu
        end do
	end if
    
	    
	! Zur Zeitberechnung
	call cpu_time(time0)
	
	
	! Loesungsschleife
	time = 0.0d0
 	timestepping: do
		time = time + dt
		if (time>tfinal) exit timestepping

        write(*,*) "Momentane Zeit:",time

		! Loeser aufrufen
		call onedswtvd(U, U0, x, dt, dx, swF, swA, swR, swRi)
		!call imp_onedswupwind(U, U0, x, dt, dx, 0.5d0, 1, swF, swA, swR, swRi)
		
		
		! Quellterm einfuegen
		do i=2,nx-1
		rcoeff = x(i)
		if (rcoeff < 1d-8) then
		    rcoeff = 0!-dt*1/1d-8
		else
		    rcoeff = -dt*1/rcoeff
		end if
		
            U%Dvectorblock(1)%datas(i)=U%Dvectorblock(1)%datas(i)&
               +rcoeff*U%Dvectorblock(2)%datas(i)
            U%Dvectorblock(2)%datas(i)=U%Dvectorblock(2)%datas(i)&
               +rcoeff*U%Dvectorblock(2)%datas(i)*U%Dvectorblock(2)%datas(i)&
                /U%Dvectorblock(1)%datas(i)
        end do
		
		
			
		! Ausgabe
		if (animate==1) then
                do i=1,nx
                   write(10,FMT='(E15.6E3,1X,E15.6E3,1X,E15.6E3)') x(i),&
                        U0%Dvectorblock(1)%datas(i),&
                        U0%Dvectorblock(2)%datas(i)    ! x, h, hu
                end do
        end if
		
		! U und U0 vertauschen
		tmpU = U
		U = U0
		U0 = tmpU
			
	end do timestepping
	
	
	
	! Ausgabe des Loesungsvektors
    do i=1,nx
       write(10,FMT='(E15.6E3,1X,E15.6E3,1X,E15.6E3)') x(i),&
            U0%Dvectorblock(1)%datas(i),&
            U0%Dvectorblock(2)%datas(i)    ! x, h, h*u
    end do
    
        
	! Ausgabe der benoetigten Zeit
	call cpu_time(time1)
	write(*,*) 'Benoetigte Zeit:', time1-time0,'Sekunden.'
	
	
	! Dateien schließen
	close(10)
	
	
	! Speicher von U und U0 freigeben
	call unmake_vectorBlock(U)
	call unmake_vectorBlock(U0)
	

end program shallow_water_prog
