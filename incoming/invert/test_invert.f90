program test_invert
  use inv
  implicit none

  integer, parameter :: N=4

  real(DP), dimension(N,N) :: A,M
  real(DP), dimension(N) :: x,y,z,f
  real(DP) :: tt0,tt1,tt2,tt3,tt4,tt5
  integer :: i,j,irun,nrun=500000


  M=0
  do i=1,N
     do j=1,N
        M(i,j)=i/j
     end do
  end do
  
  call cpu_time(tt0)
  do irun=1,nrun
     A=M; x=0; f=1
     call invert(A,f,x,N,0)
     call invert(A,f,x,N,1)
  end do
  call cpu_time(tt1)
  
  call cpu_time(tt2)
  do irun=1,nrun
     A=M; y=0; f=1
     call invert(A,f,y,N,2)
  end do
  call cpu_time(tt3)

  call cpu_time(tt4)
  do irun=1,nrun
     A=M; z=0; f=1
     call invertold(A,f,z,0)
     call invertold(A,f,z,1)
  end do
  call cpu_time(tt5)
  
  print *, sqrt(sum(abs(x-y)**2))
  print *, sqrt(sum(abs(x-z)**2))
  print *, sqrt(sum(abs(y-z)**2))

  print *, "CPU F90            ",tt1-tt0
  print *, "CPU LAPACK/explicit",tt3-tt2
  print *, "CPU F77            ",tt5-tt4

end program test_invert
