#if 0
! Undefine the keyword 'mynewline' (i.e. keep it) if this file is preprocessed 
! by f90cpp and simply remove it if another preprocessor is used instead.
#endif

#ifdef USE_PREPROC_F90CPP
#undef MYNEWLINE
#else
#define MYNEWLINE
#endif


#if 0
! Compute the Euclidean norm of a vector with 8 entries
#endif

#define euclidnorm8(U)\
  (sqrt(U(1)**2+MYNEWLINE\
	U(2)**2+MYNEWLINE\
	U(3)**2+MYNEWLINE\
	U(4)**2+MYNEWLINE\
	U(5)**2+MYNEWLINE\
	U(6)**2+MYNEWLINE\
	U(7)**2+MYNEWLINE\
	U(8)**2))


#if 0
! Compute the sum of two arguments
#endif

#define combine(a,b)\
  (a+MYNEWLINE\
   b)


#if 0
! Stringify the given argument
#endif

#define stringify(arg) #arg


#if 0
! Concatenate two arguments
#endif

#define concatenate(a,b) stringify(a ## b)
