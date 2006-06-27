program test_invert
  USE inv
  implicit none

  INTEGER, PARAMETER :: N=4

  REAL(DP), DIMENSION(N,N) :: A,M
  REAL(DP), DIMENSION(N) :: x,y,z,f
  REAL(DP) :: tt0,tt1,tt2,tt3,tt4,tt5
  INTEGER :: i,j,irun,nrun=500000


  M=0
  DO i=1,N
     DO j=1,N
        M(i,j)=i/j
     END DO
  END DO
  
  CALL cpu_time(tt0)
  DO irun=1,nrun
     A=M; x=0; f=1
     CALL invert(A,f,x,N,0)
     CALL invert(A,f,x,N,1)
  END DO
  CALL cpu_time(tt1)
  
  CALL cpu_time(tt2)
  DO irun=1,nrun
     A=M; y=0; f=1
     CALL invert(A,f,y,N,2)
  END DO
  CALL cpu_time(tt3)

  CALL cpu_time(tt4)
  DO irun=1,nrun
     A=M; z=0; f=1
     CALL invertold(A,f,z,0)
     CALL invertold(A,f,z,1)
  END DO
  CALL cpu_time(tt5)
  
  PRINT *, SQRT(SUM(ABS(x-y)**2))
  PRINT *, SQRT(SUM(ABS(x-z)**2))
  PRINT *, SQRT(SUM(ABS(y-z)**2))

  PRINT *, "CPU F90            ",tt1-tt0
  PRINT *, "CPU LAPACK/explicit",tt3-tt2
  PRINT *, "CPU F77            ",tt5-tt4

end program test_invert
