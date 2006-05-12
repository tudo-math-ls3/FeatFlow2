************************************************************************
* This file contains initialization routines for the NSDEF2 solver.
************************************************************************

************************************************************************
* Initialize assembly structures for Navier-Stokes solver
*
* This routine initializes the assembly structures that define
* the discreziation and appearence of the matrices for the Navier
* Stokes solver. If desired, COMMON-block variables configuring
* the structure are transferred.
*
* In:
*   IC2PAR : =0: initialize IASMBL/DASMBL only with standard values
*            =1: initialize IASMBL/DASMBL with standard values and
*                transfer values of COMMON-blocks into them
*
* Out:
*   IASMBL  : array [1..SZASMI] of integer
*   DASMBL  : array [1..SZASMD] of double
*             Integer and double prec. parameter block that controls the
*             discretization.
************************************************************************

      SUBROUTINE INIASM (IASMBL,DASMBL,IC2PAR)
      
      IMPLICIT NONE
      
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'cdiscretization.inc'
      
C     parameters      
      
      INTEGER IASMBL(*),IC2PAR
      DOUBLE PRECISION DASMBL(*)
      
C     Clear the structures

      CALL LCL3(IASMBL,SZASMI)
      CALL LCL1(DASMBL,SZASMD)
      
C     Initialize standard cubature formulas

      IASMBL(OICUBM) = 4
      IASMBL(OICUBA) = 4
      IASMBL(OICUBN) = 4
      IASMBL(OICUBB) = 4
      IASMBL(OICUBF) = 4
      
C     Initialize Reynolds number

      DASMBL (ORE    ) = 20.0
      DASMBL (ONY    ) = 1.0/20.0
      
C     Initialize Upwinding

      DASMBL(OUPSAM) = 0.1D0
      
C     Initialize inflow

      IASMBL(OIINFCF) = 0
      DASMBL(ODPUMAX) = 0.3
      
C     The Theta-weight must be set to 1D0, which is standard
C     for a stationary simulation. A value of 0D0 will not work.
C     An instationary solver modifies that parameter according
C     to the time step length...

      DASMBL(OTHWEIG) = 1D0
            
C     Transfger data from COMMON blocks if desired

      IF (IC2PAR.NE.0) THEN
      
        IASMBL(OISTOK)  = ISTOK
      
        IASMBL(OIELEMT) = IELT
        IASMBL(OIAPRM ) = IAPRM 
        IASMBL(OICUBM ) = ICUBM 
        IASMBL(OICUBA ) = ICUBA 
        IASMBL(OICUBN ) = ICUBN 
        IASMBL(OICUBB ) = ICUBB 
        IASMBL(OICUBF ) = ICUBF 
        IASMBL(OIMASS ) = IMASS 
        IASMBL(OIMASSL) = IMASSL
        IASMBL(OIUPW  ) = IUPW  
        IASMBL(OIPRECA) = IPRECA
        IASMBL(OIPRECB) = IPRECB
        
        IASMBL(OIINFCF) = IINFCF
        DASMBL(ODPUMAX) = DPUMAX
        
        DASMBL(OUPSAM ) = UPSAM
        DASMBL(ODMTEP ) = DMTEP
        DASMBL(ORE    ) = RE
        DASMBL(ONY    ) = 1D0/RE
      
      END IF
      
      END
      