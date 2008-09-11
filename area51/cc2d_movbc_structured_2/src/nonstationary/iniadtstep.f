************************************************************************
* This file contains routines to initialize the adaptive time stepping
* behaviour.
************************************************************************

************************************************************************
* Initialise adaptive time stepping
*
* This routine initializes the adaptive time stepping parameter block
* with standard values. If desired, COMMON block variables are
* transferred.
*
* In:
*   IC2PAR : =0: initialize IADTS/DADTS only with standard values
*            =1: initialize IADTS/DADTS with standard values and
*                transfer values of COMMON-blocks into them
*
* Out:
*   IADTS  : array [1..SZADTI] of integer
*   DADTS  : array [1..SZADTD] of double
*            Integer and double prec. parameter block that configures
*            the behaviour of the adaptive time stepping.
************************************************************************

      SUBROUTINE INIATS (IADTS,DADTS,IC2PAR)
      
      IMPLICIT NONE
      
      INCLUDE 'sadtstep.inc'
      
      INCLUDE 'ctimediscr.inc'
      
      INTEGER IC2PAR,IADTS(*)
      DOUBLE PRECISION DADTS(*)
      
C     Clear the structures

      CALL LCL3(IADTS,SZADTI)
      CALL LCL1(DADTS,SZADTD)
      
C     Only the double-prec. variables have to be initialized

      DADTS(ODTMIN ) = 0.000001
      DADTS(ODTMAX ) = 1.0
      DADTS(ODTFACT) = 9.0
      DADTS(OTIMSTP) = 0.1
      DADTS(OTIMEIN) = 0.5
      DADTS(OEPSADI) = 1.25D-1
      DADTS(OEPSADL) = 1.25D-3
      DADTS(OEPSADU) = 0.5
      DADTS(OEPSNS ) = 1D-5
      
C     if desired, transfer COMMON block variables

      IF (IC2PAR.NE.0) THEN
        IADTS(OIADTIM) = ABS(IADTIM)
        
        IF (IADTIM.GE.0) THEN
          IADTS(OIEXTIM) = 0
        ELSE
          IADTS(OIEXTIM) = 1
        END IF
        
        IADTS(OIADIN ) = IADIN
        IADTS(OIREPIT) = IREPIT
        IADTS(OIEPSAD) = IEPSAD
      
        DADTS(ODTMIN ) = DTMIN 
        DADTS(ODTMAX ) = DTMAX 
        DADTS(ODTFACT) = DTFACT
        DADTS(OTIMSTP) = TSTEP
        DADTS(OTIMEIN) = TIMEIN
        DADTS(OEPSADI) = EPSADI
        DADTS(OEPSADL) = EPSADL
        DADTS(OEPSADU) = EPSADU
        DADTS(OEPSNS ) = EPSNS 
      END IF
      
      END
      