************************************************************************
* This file contains initialization routines for the parametrization.
************************************************************************

************************************************************************
* Generate/read parametrisation and triangulation
*
* This routine initializes the parametrization or reads it in from
* the file CPARM.
*
* Unit 65 is used for reading.
*
* BPAR   - Generate information about parametrisation
*          =0: don't generate parametrization, simply initialize
*              necessary variables (can be used in pre-initialization
*              process, depending on the problem)
*          =1: initialize variables, generate parametrizaion; standard
* IPAR   - Type of paramtrization.
*          =0: FEAT parametrisation (hard-coded in PARX/PARY/TMAX)
*          =1: OMEGA parametrisation from file
* CPARM  - Filename of the file containing the parametrisation
*          (if IPAR>0)
************************************************************************

      SUBROUTINE GENPAR (BPAR, IPAR, CPARM)
      
      IMPLICIT NONE

C parameters

      LOGICAL BPAR
      INTEGER IPAR
      CHARACTER*(*) CPARM
      
C externals
      
      INTEGER MFILE
      
      MFILE = 65
      
      IF (BPAR.AND.(IPAR.EQ.1)) THEN
C Read the domain parametrisation from file
        CALL RDPARM (CPARM,MFILE)
        CLOSE(MFILE)
      ENDIF

C     No further initialization of the parametrization necessary for
C     now. This routine is rather short, isn't it? But this allowes
C     a more complex implementation in the future...

      END

************************************************************************
* Release parametrization of the domain (->RDPARM)
*
* This routine releases all parametrization information from the heap.
*
* In:  -
* Out: -
************************************************************************

      SUBROUTINE DISPAR 
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cparqdata.inc'
      
      CALL ZDISP(0,LNCOMP,'NCOMP ')
      CALL ZDISP(0,LICPTR,'ICPTR ')
      CALL ZDISP(0,LITYP,'ITYP  ')
      CALL ZDISP(0,LNSPLN,'NSPLN ')
      CALL ZDISP(0,LNPAR,'NPAR  ')
      CALL ZDISP(0,LIPPTR,'IPPTR ')
      CALL ZDISP(0,LXPAR,'DXPAR ')
      CALL ZDISP(0,LYPAR,'DYPAR ')

      END
