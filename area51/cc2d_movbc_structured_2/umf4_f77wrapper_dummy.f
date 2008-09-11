************************************************************************
* This file contains UMFPACK4 wrappers for legacy programs that can't
* use UMFPACK4. The routines here are simply empty dummies...
************************************************************************

      SUBROUTINE UMF4DEF ()
      END

      SUBROUTINE UMF4PCON ()
      END

      SUBROUTINE UMF4SYM ()
        WRITE (*,*) 'ERROR: UMFPACK4 not working!'
        STOP
      END

      SUBROUTINE UMF4NUM ()
        WRITE (*,*) 'ERROR: UMFPACK4 not working!'
        STOP
      END

      SUBROUTINE UMF4SOLR ()
      END

      SUBROUTINE UMF4SOL () 
      END

      SUBROUTINE UMF4SCAL () 
      END

      SUBROUTINE UMF4PINF () 
      END

      SUBROUTINE UMF4FNUM () 
      END

      SUBROUTINE UMF4FSYM () 
      END

      SUBROUTINE UMF4SNUM () 
      END

      SUBROUTINE UMF4SSYM () 
      END
      
      SUBROUTINE UMF4LNUM () 
      END

      SUBROUTINE UMF4LSYM () 
      END

