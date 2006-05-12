************************************************************************
* This file contains a postprocessing routine for calculating
* drag- and lift-forces. This routine stems from the default post-
* processing routine and is called in the optimization routine
* to calculate the forces without doing a complete postprocessing.
************************************************************************

************************************************************************
* Postprocessing: Calculate values of interest
*
* Calculates values of interest in the geometry. Usually this means
* calculating the forces on an arbitrary object. The method how to
* calculate can be specified by IDLMTH.
*
* In:
*   NLMIN  : minimum level 
*   NLMAX  : maximum level 
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for all levels.
*   MATDAT : array [1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialized with data
*   VECDAT : array [1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*            This structure array must specify the structure of
*            the vectors on each level. 
*
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*   IMGPAR : array [1..SZ020I+3*SZSLVI] of integer
*   DMGPAR : array [1..SZ020D+3*SZSLVD] of double
*            Integer and double parameter blocks for the multigrid 
*            sub-solver M020. 
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. 
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to fictitious boundary
*            routines. Not used in this routine.
*
*   NUVP   : Total length of solution vector
*   DUP    : array [1..NUVP] of double
*            Current solution vector
*   DRHS   : array [1..NUVP] of double
*            Right hand side for the next time step
*
*   IDLMTH : Method how to calculate forces.
*            =0: Line integration on boundary component ICMP
*            =1: Volume integration on fictitious boundary 
*                component ICMP
*            =2: Line integration on reconstructed interface
*                of fictitious boundary component ICMP
*   ICMP   : Number of the boundary component where to calculate
*            forces. If IDLMTH=1,2, ICMP=0 evaluates the sum of
*            all forces of all fictitious boundary components!
*   DPAR1,
*   DPAR2  : If ICMP is a real boundary component, DPAR1,DPAR2
*            must specify the minimum/maximum parameter calue
*            where to calculate the forces; otherwise
*            DPAR1,DPAR2 can be unspecified.
*   ICUB   : Cubature formula to use for calculating the forces.
*            Depends on the method how to calculate. If the method
*            is line integration, an identifier for a line cubature
*            formula must be given. If volume integration is used,
*            an identifier for the volume cubature formula must be
*            given.
*            Standard=2 (Trapezoidal rule) for line integration
*                    =4 (4x4 Gauss) for volume integration
*   COEF1,
*   COEF2  : Coefficients of the integral when calculating the forces
*
* Out:
*   DVALS  : array [1..*] of double precision
*            Is filled with data, depending on IDLMTH:
*            IDLMTH=0,1,2:
*              DVALS(1) = horizontal force
*              DVALS(2) = vertical force
*
*  VECDAT in that case, so the content of these vectors are treated
*  as undefined on entry and return of this routine.
*
*  The routine stops the time needed for postprocessing and adds this
*  to the timing substructure in DSTPAR!
************************************************************************
      
      SUBROUTINE DFPVIN (NLMIN,NLMAX,
     *                   TRIAS,MATDAT,VECDAT,
     *                   ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *                   NUVP,DUP,DRHS,
     *                   IDLMTH, ICMP, DPAR1,DPAR2, ICUB, COEF1,COEF2,
     *                   DVALS)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'smat2dns.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INCLUDE 'ssolvers.inc'
      INCLUDE 'stiming.inc'
      INCLUDE 'snsdef.inc'
      
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'stracking.inc'
      
C parameters

      INTEGER NUVP
      INTEGER IMGPAR(*),ISTPAR(*),IASMBL(*),IGEOM(*)
      INTEGER TRIAS(SZTRIA,NNLEV)
      INTEGER MATDAT(SZN2MI,NNLEV),VECDAT(SZN2VI,NNLEV),NLMIN,NLMAX
      DOUBLE PRECISION DMGPAR(*),DSTPAR(*),DASMBL(*)
      DOUBLE PRECISION DGEOM(*)
      
      DOUBLE PRECISION DUP(*),DRHS(*)
      
      INTEGER IDLMTH, ICMP, ICUB
      DOUBLE PRECISION COEF1,COEF2,DVALS(*),DPAR1,DPAR2
      
C     Coefficient of exact solution

      DOUBLE PRECISION UE,PE,UEX,UEY
      EXTERNAL UE,PE,UEX,UEY
      
C     finite elements

      EXTERNAL E030,E031,EM30,EM31,E010
      
C     local variables

      INTEGER NEQU,NEQP
      INTEGER KCORVG,KXNPR
      INTEGER KVERT,KMID,KVBD,KEBD,KVBDP,KBCT,KCORMG
      DOUBLE PRECISION DLEN
      
C     A backup of the solver structure for time measurements

      DOUBLE PRECISION TIMNG (SZNSDD)

C     At first stop the time and initialize some local variables

      CALL GTMAUX (TIMNG,DSTPAR,OTNLTIM-1+OTTPOST,0)

C     Fetch the vector size of velocity and pressure part on the 
C     maximum level - it's frequently used.

      NEQU = VECDAT(ONU,NLMAX)
      NEQP = VECDAT(ONP,NLMAX)

C     =================================================================
C     Calculation of the body forces.
C     Output of tracked solution values.
C
C     We calculate the forces and output tracked solution values
C     only on accepted solutions!

C     Fetch some variables from the geometry for less stuff to 
C     write :) We are working here onb the maximum level.

      KVERT  = L(TRIAS(OLVERT,NLMAX))
      KMID   = L(TRIAS(OLMID,NLMAX))
      KVBD   = L(TRIAS(OLVBD,NLMAX))
      KEBD   = L(TRIAS(OLEBD,NLMAX))
      KCORVG = L(TRIAS(OLCORVG,NLMAX))
      KCORMG = L(TRIAS(OLCORMG,NLMAX))
      KVBDP  = L(TRIAS(OLVBDP,NLMAX))
      KBCT   = L(TRIAS(OLBCT,NLMAX))
      KXNPR  = L(TRIAS(OLXNPR,NLMAX))

C     --------------------------------------------------------------      
C     At first the standard boundary forces evaluation on real 
C     boundary components.
C     Is there a cubature formula for the "real boundary evaluation"
C     given?

      IF (IDLMTH.EQ.0) THEN
      
        IF (IASMBL(OIELEMT).EQ.0) THEN
        
          CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *            TRIAS(1,NLMAX),KWORK(KVERT),
     *            KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *            DWORK(KCORVG),DWORK(KVBDP),E031,.FALSE.,
     *            ICMP,DPAR1,DPAR2,COEF1,COEF2,
     *            1,DVALS(1),DVALS(2),DLEN)   
     
        ELSE IF (IASMBL(OIELEMT).EQ.1) THEN
        
          CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *            TRIAS(1,NLMAX),KWORK(KVERT),
     *            KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *            DWORK(KCORVG),DWORK(KVBDP),E030,.FALSE.,
     *            ICMP,DPAR1,DPAR2,COEF1,COEF2,
     *            1,DVALS(1),DVALS(2),DLEN)   
     
        ELSE IF (IASMBL(OIELEMT).EQ.2) THEN
        
          CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *            TRIAS(1,NLMAX),KWORK(KVERT),
     *            KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *            DWORK(KCORVG),DWORK(KVBDP),EM31,.TRUE.,
     *            ICMP,DPAR1,DPAR2,COEF1,COEF2,
     *            1,DVALS(1),DVALS(2),DLEN)   
     
        ELSE IF (IASMBL(OIELEMT).EQ.3) THEN
        
          CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *            TRIAS(1,NLMAX),KWORK(KVERT),
     *            KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *            DWORK(KCORVG),DWORK(KVBDP),EM30,.TRUE.,
     *            ICMP,DPAR1,DPAR2,COEF1,COEF2,
     *            1,DVALS(1),DVALS(2),DLEN)   
     
        END IF 
        
      END IF ! IDLMTH=0

C     --------------------------------------------------------------      
C     Boundary forces evaluation on fictitious boundary 
C     by volume integration.

      IF (IDLMTH.GT.0) THEN
      
C       Call the evaluation routine, calculate forces on that
C       boundary component
     
        IF (IASMBL(OIELEMT).EQ.0) THEN
        
          CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *             DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *             KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *             E031,.FALSE.,DVALS(1),DVALS(2),COEF1,COEF2,
     *             ICMP,0.5D0,IGEOM,DGEOM)
        
        ELSE IF (IASMBL(OIELEMT).EQ.1) THEN
        
          CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *             DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *             KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *             E030,.FALSE.,DVALS(1),DVALS(2),COEF1,COEF2,
     *             ICMP,0.5D0,IGEOM,DGEOM)

        ELSE IF (IASMBL(OIELEMT).EQ.2) THEN
        
          CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *             DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *             KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *             EM31,.TRUE.,DVALS(1),DVALS(2),COEF1,COEF2,
     *             ICMP,0.5D0,IGEOM,DGEOM)

        ELSE IF (IASMBL(OIELEMT).EQ.3) THEN
        
          CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *             DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *             KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *             EM30,.TRUE.,DVALS(1),DVALS(2),COEF1,COEF2,
     *             ICMP,0.5D0,IGEOM,DGEOM)

        END IF 
        
      END IF ! IDLMTH=1

C     --------------------------------------------------------------      
C     Boundary forces evaluation on fictitious boundary 
C     by line integration.

      IF (IDLMTH.EQ.2) THEN
      
C       Call the evaluation routine, calculate forces on that
C       boundary component.
C       Pass the assembly data structures as INFO-block to
C       these routines. That way all called subroutines
C       can access information about the current assembly
C       status.
     
        IF (IASMBL(OIELEMT).EQ.0) THEN
        
          CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                 DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                 KWORK(KMID),TRIAS(1,NLMAX),
     *                 E031,.FALSE.,DVALS(1),DVALS(2),COEF1,
     *                 COEF2,ICMP,IGEOM,DGEOM)
        
        ELSE IF (IASMBL(OIELEMT).EQ.1) THEN
        
          CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                 DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                 KWORK(KMID),TRIAS(1,NLMAX),
     *                 E030,.FALSE.,DVALS(1),DVALS(2),COEF1,
     *                 COEF2,ICMP,IGEOM,DGEOM)

        ELSE IF (IASMBL(OIELEMT).EQ.2) THEN
        
          CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                 DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                 KWORK(KMID),TRIAS(1,NLMAX),
     *                 EM31,.TRUE.,DVALS(1),DVALS(2),COEF1,
     *                 COEF2,ICMP,IGEOM,DGEOM)

        ELSE IF (IASMBL(OIELEMT).EQ.3) THEN
        
          CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                 DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                 KWORK(KMID),TRIAS(1,NLMAX),
     *                 EM30,.TRUE.,DVALS(1),DVALS(2),COEF1,
     *                 COEF2,ICMP,IGEOM,DGEOM)

        END IF 
          
      END IF ! IDLMTH=2
      
C     Finally stop the time we need for postprocessing:

      CALL GTMAUX (TIMNG,DSTPAR,OTNLTIM-1+OTTPOST,1)
      
      END

