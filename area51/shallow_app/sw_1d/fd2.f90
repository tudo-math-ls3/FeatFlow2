! ********** HILFSROUTINEN **********

module fdschemes

  implicit none

contains

  ! ***** Forward Algorithmus *****
  subroutine forward(u0,u,a,dx,dt)

    ! Eingabe:
    ! u0 = Ausgangsvektor
    ! u = Ausgabevektor
    ! a = Geschwindigkeit (>0)
    ! dx = Ortsschrittweite
    ! dt = Zeitschrittweite

    double precision, intent(in)			:: a, dx, dt
    double precision, dimension(:), intent(out)		:: u
    double precision, dimension(:), intent(in)		:: u0
    integer						:: i,j
    double precision					:: k1

    ! Hilfsvariable
    k1 = a*dt/dx

    ! Berechne den Lösungsvektor für alle außer dem linken x (parallel)
    forall(i=lbound(u0,1)+1:ubound(u0,1))
       u(i) = u0(i)-k1*(u0(i)-u0(i-1))
    end forall

    ! Berechne den Lösungswert für das linke x
    u(lbound(u0,1)) = u0(lbound(u0,1))+k1*u0(lbound(u0,1))

  end subroutine forward




  ! ***** generalisierter Upwind und Lax-Wendroff *****
  subroutine gul(u0, u, x, dt, dx, funcf, funcfstrich, alg)
	! Löst:
	! d/dt u + d/dx f(u) = 0
	! u(x,0) = u0(x)

	! Eingabeparameter:
	! u0            Startwert-Vektor (Zeilenvektor)
	! x		Ortsdiskretisierung (Zeilenvektor)
	! dt            Länge des Zeitschritts
	! dx            Schrittweite des Ortsgitters
	! funcf         Funktionhandle auf die Funktion f(u)
	! funcfstrich   Funktionhandle auf die Funktion f'(u)
	! alg           wählt den zu verwendenden Algorithmus
	!               0: generalisierter Upwind
	!               1: generalisierter Lax Wendroff

	double precision, intent(in)				:: dx, dt, alg
	double precision, dimension(:), intent(out)		:: u
	double precision, dimension(:), intent(in)		:: u0, x
	integer							:: i,j
	double precision					:: k1

	interface
	double precision function funcf(u,i,x)
		double precision, intent(in)	:: u, x
		integer, intent(in)		:: i
    	end function funcf
    	double precision function funcfstrich(u,i,x)
		double precision, intent(in)	:: u, x
		integer, intent(in)		:: i
    	end function funcfstrich
	end interface

	! Lösungsvektor: Randbedingung
	u(lbound(u0,1)) = u0(lbound(u0,1))
	u(ubound(u0,1)) = u0(ubound(u0,1))


	! Hilfsvariablen
	k1 = dt/dx
	
	! Berechne alle inneren Werte

	gulloop: do i = lbound(u0,1)+1,ubound(u0,1)-1
		u(i) = u0(i) - k1*(fl(i)-fl(i-1))
	end do gulloop

  contains

  ! Hilfsfunktionen für gul
  ! flux-function f_{i+1/2}  =f^L für alpha = 0, f^H für alpha = 1
    double precision function fl(i)
	double precision	:: k2, tempa
	integer, intent(in) 	:: i
        tempa = approxa(i);
        k2 = tempa*k1;    ! = nu
        fl = (funcf(u0(i+1),i+1,x(i+1))+funcf(u0(i),i,x(i)))/2.0d0-abs(tempa)/2.0d0*(u0(i+1)-u0(i))&
		+alg*abs(tempa)/2.0d0*(1-abs(k2))*(u0(i+1)-u0(i))
    end function fl
  ! berechnet a_{i+1/2}
    double precision function approxa(i)
	integer, intent(in)	:: i
        if (abs(u0(i+1)-u0(i))< 0.001d0) then
            approxa = funcfstrich(u0(i),i,x(i))
        else
            approxa = (funcf(u0(i+1),i+1,x(i+1))-funcf(u0(i),i,x(i)))/(u0(i+1)-u0(i))
        end if
    end function approxa

 end subroutine gul





  ! ***** generalisierter Upwind und Lax-Wendroff *****
  subroutine gtvd(u0, u, x, dt, dx, funcf, funcfstrich, limiter)
	! Löst:
	! d/dt u + d/dx f(u) = 0
	! u(x,0) = u0(x)

	! Eingabeparameter:
	! u0            Startwert-Vektor (Zeilenvektor)
	! x		Ortsdiskretisierung (Zeilenvektor)
	! dt            Länge des Zeitschritts
	! dx            Schrittweite des Ortsgitters
	! funcf         Funktionhandle auf die Funktion f(u)
	! funcfstrich   Funktionhandle auf die Funktion f'(u)
	! limiter	wählt den zu verwendenden Limiter (1=minmod, 2=mc, 3=vl, 4=superbee)

	double precision, intent(in)				:: dx, dt
	double precision, dimension(:), intent(out)		:: u
	double precision, dimension(:), intent(in)		:: u0, x
	integer, intent(in)					:: limiter
	integer							:: i,j
	double precision					:: k1

	interface
	double precision function funcf(u,i,x)
		double precision, intent(in)	:: u, x
		integer, intent(in)		:: i
    	end function funcf
    	double precision function funcfstrich(u,i,x)
		double precision, intent(in)	:: u, x
		integer, intent(in)		:: i
    	end function funcfstrich
	end interface

	! Lösungsvektor: Randbedingungen
	u(lbound(u0,1)) = u0(lbound(u0,1))
	u(lbound(u0,1)+1) = u0(lbound(u0,1)+1)
	u(ubound(u0,1)) = u0(ubound(u0,1))
	u(ubound(u0,1)-1) = u0(ubound(u0,1)-1)

	! Hilfsvariablen
	k1 = dt/dx
	
	! Berechne alle inneren Werte
	tvdloop: do i = lbound(u0,1)+2,ubound(u0,1)-2
		u(i) = u0(i) - k1*(fl(i)-fl(i-1))
	end do tvdloop

  contains

  ! Hilfsfunktionen für gtvd
  ! flux-function f_{i+1/2}  =f^L für alpha = 0, f^H für alpha = 1 (alpha = erg. des limiters)
    double precision function fl(i)
	double precision	:: k2, tempa
	integer, intent(in) 	:: i
	integer			:: j

        tempa = approxa(i)
        k2 = tempa*k1    ! = nu

	j=i-int(sign(1d0,tempa))          !???????????????????????????????????????????????????????????????
	!if (sign(tempa,1.0d0)>0d0) then
	!j=i-1
	!else
	!j=1+1
	!end if


        fl = (funcf(u0(i+1),i+1,x(i+1))+funcf(u0(i),i,x(i)))/2.0&
		-abs(tempa)/2.0*(u0(i+1)-u0(i))&
		 + limiterfunc(u0(i+1)-u0(i),u0(j+1)-u0(j))*abs(tempa)/2.0*(1-abs(k2))
    end function fl
  ! berechnet a_{i+1/2}
    double precision function approxa(i)
	integer, intent(in)	:: i
        if (abs(u0(i+1)-u0(i))< 0.001d0) then
            approxa = funcfstrich(u0(i),i,x(i))
        else
            approxa = (funcf(u0(i+1),i+1,x(i+1))-funcf(u0(i),i,x(i)))/(u0(i+1)-u0(i))
        end if
    end function approxa
  ! limiter function
    double precision function limiterfunc(a,b)
	double precision, intent(in)		:: a,b
	double precision			:: temp5,na,nb
        temp5=(sign(1.0d0,a)+sign(1.0d0,b))/2.0d0   ! ???????????????????
        na=abs(a)
        nb=abs(b)
        if (temp5.ne.0d0) then
        select case(limiter)
            case (1)
                    limiterfunc = temp5*min(na,nb)
            case (2)
                    if (na+nb < 0.00001d0) then
                        limiterfunc=1
                    else
                        limiterfunc = temp5*2.0d0*na*nb/(na+nb)
                    end if
            case (3)
                    limiterfunc = temp5*min(2.0d0*na,(na+nb)/2.0d0,2.0d0*nb)
            case (4)
                    limiterfunc = temp5*max(min(2.0d0*na,nb),min(na,2.0d0*nb))
        end select
        else
            limiterfunc=0d0
        end if
    end function limiterfunc

 end subroutine gtvd



  ! ***** Benchmarkfunctions *****
    double precision function bench1funcf(u,i,x)
	double precision, intent(in)	:: u,x
	integer, intent(in)		:: i
        bench1funcf = x*u
    end function bench1funcf
    double precision function bench1funcfstrich(u,i,x)
	double precision, intent(in)	:: u,x
	integer, intent(in)		:: i
        bench1funcfstrich = x
    end function bench1funcfstrich
    double precision function bench2funcf(u,i,x)
	double precision, intent(in)	:: u,x
	integer, intent(in)		:: i
        bench2funcf = (0.5d0-x)*u
    end function bench2funcf
    double precision function bench2funcfstrich(u,i,x)
	double precision, intent(in)	:: u,x
	integer, intent(in)		:: i
        bench2funcfstrich = 0.5d0-x
    end function bench2funcfstrich
    double precision function bench3funcf(u,i,x)
	double precision, intent(in)	:: u,x
	integer, intent(in)		:: i
        bench3funcf = u*u/2d0
    end function bench3funcf
    double precision function bench3funcfstrich(u,i,x)
	double precision, intent(in)	:: u,x
	integer, intent(in)		:: i
        bench3funcfstrich = u
    end function bench3funcfstrich


end module fdschemes


! ********** HAUPTPROGRAMM **********

program finite_differences
  use fdschemes

  implicit none

  integer, parameter				:: n = 101		! Anzahl der x-Gitterpunkte
  double precision, parameter			:: tfinal = 0.4		! Endzeit
  double precision				:: a, dx, dt, time
  double precision				:: time0, time1		! Zur Zeitmessung
  double precision, dimension(1:n)		:: u0 ,u ,x		! Startvektor, Lösungsvektor, Gittervektor
  integer					:: i			! Laufvariable

  ! Parameter initialisieren
  a = 1.0				! Welleneschwindigkeit (für die linear advection equation)
  dx = 1.0d0/(n-1)			! dx
  dt = 0.001d0				! dt

  ! Gittervektor initialisieren
  forall(i=1:n)
     x(i)=(i-1)*dx
  end forall

  ! Startwerte festlegen
  u0 = 0
  where (x>=0.1d0 .AND. x<=0.4d0)
     u0 = 1
  end where

  ! Datei öffnen
  open(10,file='ausgabe.dat')

  ! Startvektor in Datei schreiben
  write(10,*) u0

  ! Zeitschritte ausführen und Ergebnis in Datei schreiben
  ! mit Zeitmessung
  call cpu_time(time0)
  time = 0.0
  timestepping: do
     time = time + dt
     if (time>tfinal) exit timestepping
     call gtvd(u0, u, x, dt, dx, bench3funcf, bench3funcfstrich, 4)
     u0 = u
     write(10,*) u
  end do timestepping
  call cpu_time(time1)

  ! Datei schließen
  close(10)

  ! Zeitaufwand ausgeben
  print*,'Zeitaufwand: ',time1-time0,' Sekunden.'
end program finite_differences
