module weakDirichlet

  use fsystem
  use storage
  use genoutput
  use linearsolver
  use boundary
  use triangulation
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use derivatives
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use bcassembly
  use element
  implicit none
	
  ! mass matrices
  type(t_matrixScalar) :: rmassmatrixWeakD_u
  type(t_matrixScalar) :: rmassmatrixWeakD_c
  
  ! lambdas
  real(DP) :: lambda_c  
  real(DP) :: lambda_u  
    
end module 