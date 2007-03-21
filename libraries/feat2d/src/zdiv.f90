! This routine returns how many integers fit into a double precision 
! variable.

  ELEMENTAL INTEGER FUNCTION ZILEND()
  
  IMPLICIT NONE
  
    DOUBLE PRECISION, PARAMETER :: D = 0
    
    ! (Ab)using KIND to determine how many integers/reals fit into a double
    ! is nonstandard, so let's hope there is no strange machine architecture
    ! using different values for KIND as the number of bytes!
    ZILEND = KIND(D)/KIND(0)
  
  END FUNCTION
  
! This routine returns how many reals fit into a double precision 
! variable.

  ELEMENTAL INTEGER FUNCTION ZVLEND()
  
  IMPLICIT NONE
  
    DOUBLE PRECISION, PARAMETER :: D = 0
    
    ! (Ab)using KIND to determine how many integers/reals fit into a double
    ! is nonstandard, so let's hope there is no strange machine architecture
    ! using different values for KIND as the number of bytes!
    ZVLEND = KIND(D)/KIND(0.0)
  
  END FUNCTION
  
