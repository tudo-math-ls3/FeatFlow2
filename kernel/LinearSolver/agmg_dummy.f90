!##############################################################################
!# ****************************************************************************
!# <name> agmg_dummy </name>
!# ****************************************************************************
!#
!# <purpose>
!# This fortran file serves as dummy for the AGMG library which provies an
!# aggregation-based algebraic multigrid method by Yvan Notay.
!#
!# An free academic license of this library is available at:
!# http://homepages.ulb.ac.be/~ynotay/AGMG/
!#
!# This module provides the dummy routines  dagmg, cagmg, sagmg, zagmg which
!# serve as drivers for the algebraic multigrid solver in single, double, 
!# single complex and double complex precision.
!# </purpose>
!##############################################################################

subroutine cagmg( n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol )
  use fsystem
  use genoutput
  implicit none
  
  integer    :: n,ia(n+1),ja(*),ijob,iprint,nrest,iter
  complex (kind(0.0e0)) :: a(*),f(n),x(n)
  real (kind(0.0e0)) :: tol
  
  call output_line ("This is just a dummy for the original CAGMG routine provided by the AGMG library!", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call output_line ("Please obtain and appropriate license from http://homepages.ulb.ac.be/~ynotay/AGMG/", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call output_line ("and configure your application using option --buildlib=agmg to enable use of AGMG!", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call sys_halt()
  
end subroutine cagmg

! *****************************************************************************

subroutine dagmg( n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol )
  use fsystem
  use genoutput
  implicit none
  
  integer    :: n,ia(n+1),ja(*),ijob,iprint,nrest,iter
  real (kind(0.0d0)) :: a(*),f(n),x(n)
  real (kind(0.0d0)) :: tol
  
  call output_line ("This is just a dummy for the original DAGMG routine provided by the AGMG library!", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call output_line ("Please obtain and appropriate license from http://homepages.ulb.ac.be/~ynotay/AGMG/", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call output_line ("and configure your application using option --buildlib=agmg to enable use of AGMG!", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call sys_halt()
  
end subroutine dagmg

! *****************************************************************************

subroutine sagmg( n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol )
  use fsystem
  use genoutput
  implicit none
  
  integer    :: n,ia(n+1),ja(*),ijob,iprint,nrest,iter
  real (kind(0.0e0)) :: a(*),f(n),x(n)
  real (kind(0.0e0)) :: tol
    
  call output_line ("This is just a dummy for the original SAGMG routine provided by the AGMG library!", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call output_line ("Please obtain and appropriate license from http://homepages.ulb.ac.be/~ynotay/AGMG/", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call output_line ("and configure your application using option --buildlib=agmg to enable use of AGMG!", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call sys_halt()
  
end subroutine sagmg


! *****************************************************************************

subroutine zagmg( n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol )
  use fsystem
  use genoutput
  implicit none
  
  integer    :: n,ia(n+1),ja(*),ijob,iprint,nrest,iter
  complex (kind(0.0d0)) :: a(*),f(n),x(n)
  real (kind(0.0d0)) :: tol
    
  call output_line ("This is just a dummy for the original ZAGMG routine provided by the AGMG library!", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call output_line ("Please obtain and appropriate license from http://homepages.ulb.ac.be/~ynotay/AGMG/", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call output_line ("and configure your application using option --buildlib=agmg to enable use of AGMG!", &
      OU_CLASS_ERROR, OU_MODE_STD, "dagmg")
  call sys_halt()
  
end subroutine zagmg
