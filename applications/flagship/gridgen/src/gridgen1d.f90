program gridgen1d

  implicit none

  integer :: i, n=128

  open(unit=55, file='grid1d.tri')

  write(unit=55,fmt='(A)') 'Gitterbeschreibung'
  write(unit=55,fmt='(A)') 'Parametrisierung PARCX, PARCY, TMAXC'
  write(unit=55,fmt='(I3,X,I3,X,I3,X,I3,X,I3,X,A)') n,n+1,0,2,2,'NEL NVT NMT NVE NBCT'

  write(unit=55,fmt='(A)') 'DCORVG'
  write(unit=55,fmt='(F10.8)') 0.0
  write(unit=55,fmt='(F10.8)') 1.0
  do i = 1,n-1
    write(unit=55,fmt='(F10.8)') 1.0/n*i
  end do

  write(unit=55,fmt='(A)') 'KVERT'
  write(unit=55,fmt='(I3,X,I3)') 1,3
  do i = 3,n
    write(unit=55,fmt='(I3,X,I3)') i, i+1
  end do
  write(unit=55,fmt='(I3,X,I3)') n+1,2

  write(unit=55,fmt='(A)') 'KNPR'
  write(unit=55,fmt='(I1)') 1
  write(unit=55,fmt='(I1)') 2
  do i = 1,n-1
    write(unit=55,fmt='(I1)') 0
  end do

  write(unit=55,fmt='(A)') 'KMM'
  write(unit=55,fmt='(I1,X,I1,X,I1,X,I1)') 1,1,2,2

  close(unit=55)

end program gridgen1d
