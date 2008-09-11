************************************************************************
* This file collects the initialization routines for NONSTEAD.
*
* The routines here build the general parameter structures for the
* call to NONST2 solver.
************************************************************************

************************************************************************
* Initialize Nonstationary solver structures
*
* This initializes the IPARAM/DPARAM structure arrays, which are
* used in the call of NONSTEAD, with standard values. If desired,
* values of the COMMON-blocks are transferred to these structure
* arrays.
*
* The routine does not allocate any memory, it just initializes the
* parameters.
*
* In:
*   IPARAM : array [1..SZNSDI] of integer
*            Integer parameter block for NSDEF
*   DPARAM : array [1..SZNSDD] of double
*            Double precision parameter block for NSDEF
*   IC2PAR : =0: initialize IPARAM/DPARAM only with standard values
*            =1: initialize IPARAM/DPARAM with standard values and
*                transfer values of COMMON-blocks into them
*
* Out:
*   IPARAM,
*   DPARAM : Initialized structures
*
* This routine does not allocate any memory and does not initialize
* the DString-handles of the basic filenames for file output. This
* has to be done by the caller if necessary!
************************************************************************

      SUBROUTINE INIISD (IPARAM,DPARAM,IC2PAR)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'stiming.inc'
      INCLUDE 'snonstead.inc'
      
      INCLUDE 'ctimediscr.inc'
      INCLUDE 'cdiscretization.inc'
      
      INTEGER IPARAM(SZISDI),IC2PAR
      DOUBLE PRECISION DPARAM(SZISDD)
      
C     Clean the arrays
      
      CALL LCL1 (DPARAM,SZISDD)
      CALL LCL3 (IPARAM,SZISDI)
      
C     Set up standard values.
C     At first call the standard initialization routine for the standard
C     solver structure.

      CALL INGSLV (IPARAM, DPARAM)

C     Some standard values are not used by NONSTEAD, set them to 0.
C     Others have different <> 0 standard values - initialize them.

      IPARAM (ONITMAX) = 20
      IPARAM (OTIM  )  = 2
      IPARAM (OIASRLN) = 0
      IPARAM (OMSGTRM) = 2
      IPARAM (OMTRMRS) = 1

      DPARAM (OEPSREL) = 0D0
      DPARAM (OEPSABS) = 0D0
      DPARAM (ODIVREL) = 0D0
      DPARAM (ODIVABS) = 0D0
      DPARAM (OVECZER) = 0D0
            
C     Initialize extended parameters which are different from 0.

      IPARAM (OSLTAG)  = 21
      
      DPARAM (OTIMEMX) = 100.0D0
      DPARAM (OSSTHETA) = 0.5D0
      
      DPARAM (ODMSMTM)  = -1D0
      
C     Should we transfer COMMON-block variables?

      IF (IC2PAR.NE.0) THEN
      
        IPARAM (ONITMAX) = NITNS
        IPARAM (OMSGTRM) = MT
        IPARAM (OIRHS)   = IRHS
        IPARAM (OIBDR)   = IBDR
        IPARAM (OIDMESH) = IDMESH
        IPARAM (OIMSCRT) = IMSCRT
        IPARAM (OIMSREP) = IMSREP
        IPARAM (OIMSSTP) = IMSSTP
        IPARAM (OFRSTP)  = IFRSTP
        
        DPARAM (OTIMEMX)  = TIMEMX
        DPARAM (OSSTHETA) = THETA
        DPARAM (ODMSMTM)  = DMSMTM
        
C       No stopping criteria in the parameter block of the solver;
C       they can be found in the parameters of the adaptive time
C       stepping.
      
      END IF
      
      END

