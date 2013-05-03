module vartypes
! Definition der Variablentypen

     implicit none

     integer, parameter :: dp = selected_real_kind(15,307)
     integer, parameter :: sp = selected_real_kind(6,37)


     ! Skalarer Vektor
     type t_vectorScalar
     integer 											:: neq			! Anzahl der Eintraege
     real(dp), dimension(:), pointer					:: datas
     end type t_vectorScalar

	! Blockvektor
     type t_vectorBlock
        integer											:: nblocks		! Anzahl der skalaren Vektoren
        type(t_vectorScalar), dimension(:), pointer		:: DvectorBlock
     end type t_vectorBlock

	! Skalare Tridiagonalmatrix
	type t_matrixTridiag
        integer											:: idimension	! idim x idim - Matrix
        type(t_vectorScalar), pointer					:: upper
        type(t_vectorScalar), pointer					:: main			! Hauptdiagonale
        type(t_vectorScalar), pointer					:: lower
     end type t_matrixTridiag

	! BlockTridiagonalmatrix
	type t_blockMatrixTridiag
        integer											:: idim			! Anzahl der Tridiag.matrizen
        type(t_matrixTridiag), dimension(:), pointer	:: trimat
     end type t_blockMatrixTridiag



contains


	! Erzeuge skalaren Vektor
	subroutine make_vectorScalar(SV,numEntries)
		type(t_vectorScalar), intent(inout)	:: SV						! der zu erstellende skalare Vector
		integer, intent(in)					:: numEntries				! Anzahl der Eintraege
    		allocate(SV%datas(numEntries))
			SV%neq = numEntries
	end subroutine make_vectorScalar


	! Gib skalaren Vektor frei
	subroutine unmake_vectorScalar(SV)
		type(t_vectorScalar), intent(inout)	:: SV						! der zu vernichtende skalare Vector
		deallocate(SV%datas)
	end subroutine unmake_vectorScalar


    ! Erzeuge Blockvektor
	subroutine make_vectorBlock(BV,numBlocks,numEntries)
		type(t_vectorBlock), intent(inout)	:: BV						! der zu erstellende Blockvector
		integer, intent(in)					:: numBlocks, numEntries	! Anzahl der Vektoren und Eintraege je Vektor
		integer								:: i						! Zaehlvariable
    	allocate(BV%DvectorBlock(numBlocks))
		BV%nblocks = numBlocks
		do i = 1, numBlocks
			call make_vectorScalar(BV%DvectorBlock(i),numEntries)
		end do
	end subroutine make_vectorBlock


	! Gib Blockvektor frei
	subroutine unmake_vectorBlock(BV)
		type(t_vectorBlock), intent(inout)	:: BV						! der zu vernichtende Blockvector
		integer								:: i						! Zaehlvariable
		do i = 1, BV%nblocks
			call unmake_vectorScalar(BV%DvectorBlock(i))
		end do
		deallocate(BV%DvectorBlock)
		BV%nblocks = 0
	end subroutine unmake_vectorBlock


	! Erzeuge Tridiagonalmatrix
	subroutine make_matrixTridiag(TM,idim)
		type(t_matrixTridiag), intent(inout)	:: TM					! die zu erstellende Tridiagonalmatrix
		integer, intent(in)						:: idim					! eine idim x idim - Matrix
		TM%idimension = idim
		allocate(TM%upper)
		allocate(TM%main)
		allocate(TM%lower)
    	allocate(TM%upper%datas(idim-1))
    	TM%upper%neq = idim-1
		allocate(TM%main%datas(idim))
		TM%main%neq = idim
		allocate(TM%lower%datas(idim-1))
		TM%lower%neq = idim-1
	end subroutine make_matrixTridiag


	! Gib Tridiagonalmatrix frei
	subroutine unmake_matrixTridiag(TM)
		type(t_matrixTridiag), intent(inout)	:: TM					! die zu vernichtende Tridiagonalmatrix
		TM%idimension = 0
    	deallocate(TM%upper%datas)
		deallocate(TM%main%datas)
		deallocate(TM%lower%datas)
		deallocate(TM%upper)
		deallocate(TM%main)
		deallocate(TM%lower)
	end subroutine unmake_matrixTridiag


	! Erzeuge Tridiagonalblockmatrix
	subroutine make_blockMatrixTridiag(bTM,numTri,idim)
		type(t_blockMatrixTridiag), intent(inout)	:: bTM					! die zu erstellende BlockTridiagonalmatrix
		integer, intent(in)							:: numTri, idim			! Anzahl der TridiagMatrizen und Dimension der Matrizen
		integer										:: i					! Zaehlvariable
		bTM%idim = numTri
		allocate(bTM%trimat(numTri))
		do i = 1, numTri
			call make_matrixTridiag(bTM%trimat(i),idim)
		end do
	end subroutine make_blockMatrixTridiag


	! Gib Tridiagonalblockmatrix frei
	subroutine unmake_blockMatrixTridiag(bTM)
		type(t_blockMatrixTridiag), intent(inout)	:: bTM					! die zu vernichtende BlockTridiagonalmatrix
		integer										:: i					! Zaehlvariable
		do i = 1, bTM%idim
			call unmake_matrixTridiag(bTM%trimat(i))
		end do
		deallocate(bTM%trimat)
	end subroutine unmake_blockMatrixTridiag


    ! Matrix-Vektor Multiplikation
    subroutine MatVecMul(matrix,vector,result)					! berechnet: result = Matrix * vector
    	implicit none
    	type(t_matrixTridiag), intent(in)			:: matrix
    	type(t_vectorScalar), intent(in)			:: vector
    	type(t_vectorScalar), intent(inout)			:: result
    	integer										:: i		! Zaehlvariable
    	integer										:: last

    	last = result%neq

    	! Teste, ob die Dimensionen stimmen
    	if ((matrix%idimension==vector%neq).and.(last==vector%neq)) then

    		! ersten Eintrag des Loesungsvektors berechnen
    		result%datas(1) = matrix%main%datas(1) * vector%datas(1) + matrix%upper%datas(1) * vector%datas(2)

    		! berechne die inneren Eintrage des Loesungsvektors
    		do i = 2, (last-1)
    			result%datas(i) = matrix%lower%datas(i-1) * vector%datas(i-1) &
    							+ matrix%main%datas(i) * vector%datas(i) &
    							+ matrix%upper%datas(i) * vector%datas(i+1)
    		end do

    		! letzten Eintrag des Loesungsvektors berechnen
    		result%datas(last) = matrix%lower%datas(last-1) * vector%datas(last-1) + matrix%main%datas(last) * vector%datas(last)

    	else ! Die Dimensionen sind falsch
    		write(*,*) 'Matrix-Vektor-Multiplikation: Die Dimensionen passen nicht!'
    		! stop
    	end if

    end subroutine MatVecMul


	! Linearer Loeser mit Thomas-Algorithmus			ACHTUNG: nur fuer diagonaldominante tridiag. Matrizen
	subroutine linSolvThomas(A,b,x)								! Loest A * x = b
		implicit none
	   	type(t_matrixTridiag), intent(in)			:: A		! Koeffizientenmatrix (tridiag.)
    	type(t_vectorScalar), intent(in)			:: b		! rechte Seite
    	type(t_vectorScalar), intent(inout)			:: x		! Loesungsvektor
    	type(t_vectorScalar)						:: c, d		! Hilfsvektoren
    	integer										:: i		! Zaehlvariable
    	integer										:: dimen	! Groesse des LGS
    	real(dp)									:: hk		! Hilfsvar

    	dimen = x%neq

    	! Teste, ob die Dimensionen stimmen
    	if ((A%idimension==b%neq).and.(dimen==b%neq)) then
    		! Speicher fuer Hilfsvariablen holen
    		call make_vectorScalar(c,dimen)
    		call make_vectorScalar(d,dimen)

    		! Hilfsvektoren berechnen
    		if (A%main%datas(1)==0) then
    			write(*,*) 'linSolvThomas: Division durch 0'
    			stop
    		end if
    		c%datas(1) = A%upper%datas(1)/A%main%datas(1)
    		d%datas(1) = b%datas(1)/A%main%datas(1)
    		do i = 2, dimen-1
    			hk = (A%main%datas(i)-(c%datas(i-1)*A%lower%datas(i-1)))
    			if (hk==0) then
    				write(*,*) 'linSolvThomas: Division durch 0'
    				stop
    			end if
    			c%datas(i) = A%upper%datas(i)/hk
    			d%datas(i) = (b%datas(i)-d%datas(i-1)*A%lower%datas(i-1))/hk
    		end do
    		hk = (A%main%datas(dimen)-(c%datas(dimen-1)*A%lower%datas(dimen-1)))
    		d%datas(dimen) = (b%datas(dimen)-d%datas(dimen-1)*A%lower%datas(dimen-1))/hk

    		! Rueckwaerts einsetzen
    		x%datas(dimen)=d%datas(dimen)
    		i = dimen-1
    		do
    			x%datas(i)=d%datas(i)-(c%datas(i)*x%datas(i+1))
    			i=i-1
    			if(i==0) exit
  			end do

    		! Speicher fuer Hilfsvariablen freigeben
    		call unmake_vectorScalar(c)
    		call unmake_vectorScalar(d)

    	else ! Die Dimensionen sind falsch
    		write(*,*) 'LinSolvThomas: Die Dimensionen passen nicht!'
    	end if

    end subroutine linSolvThomas


end module vartypes
