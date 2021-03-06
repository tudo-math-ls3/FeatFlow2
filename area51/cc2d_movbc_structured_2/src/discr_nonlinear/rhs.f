***********************************************************************
* This file encapsules the generation of the RHS vector
* for te stationary and nonstationary Navier Stokes system.
*
***********************************************************************
* Basically, the RHS is generated by calling the FDATIN routine,
* which is user provided. 
***********************************************************************

***********************************************************************
* Generate RHS vector
*
* This routine generates a RHS vector based on the information about
* discretization, geometries and triangulation. It can be used as a
* default generation routine in the nonstationary solver as well as
* for generating the RHS in the stationary solver by setting TIMENS
* to 0.
*
* In:
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation structures of the underlying mesh.
*   IPARAM : array [1..*] of integer
*   DPARAM : array [1..*] of integer
*            Integer/Double prec. parameter block of the current solver,
*            for which the RHS should be generated. This is passed to
*            the user defined callback routine FDATIN and not used here.
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization.
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*   TIMENS : Current simulational time.
*            Must be set to 0D0 for stationary simulation. 
*   VECDAT : array [1..SZN2VI] of integer 
*            TNS2DVectorParams-structure, corresponding to TRIA.
*            Defines the form of the RHS vector.
*   
* Out:
*   DRHS   : array [1..*] of doule precision
*            Generated RHS vector; usually length 2*NEQU+NEQP.
***********************************************************************
            
      SUBROUTINE GNRHSV (TRIA,
     *                   IPARAM,DPARAM,
     *                   IASMBL,DASMBL,IGEOM,DGEOM,
     *                   VECDAT,DRHS)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INCLUDE 'stria.inc'
      
C     Coefficient of stiffness matrix, right hand side, exact solution
      
      DOUBLE PRECISION RHS
      EXTERNAL RHS
      
C     definition of finite elements
      
      EXTERNAL E030,E031,EM30,EM31,E010
      
C     parameters
       
      INTEGER IPARAM(*),IASMBL(*),TRIA(*),IGEOM(*),VECDAT(SZN2VI)
      DOUBLE PRECISION DPARAM(*),DASMBL(*),DRHS(*),DGEOM(*)

C     constants for proper operation

      INTEGER NBLOCF
      PARAMETER (NBLOCF=2)

C     Names of vectors

      CHARACTER ARRDF*6
      DIMENSION ARRDF(NBLOCF)
      DATA ARRDF/'DF1   ','DF2   '/
      SAVE ARRDF

C     constants for forming the linear form

      INTEGER KFN(NBLOCF),KF(NNAB,NBLOCF)
      
      LOGICAL BSNGLF(NBLOCF)      
      DATA BSNGLF /.FALSE.,.FALSE./
      SAVE BSNGLF

      LOGICAL BCONF(NBLOCF)
      DATA BCONF /.FALSE.,.FALSE./
      SAVE BCONF
      
C     local variables

      INTEGER LF(NBLOCF),ICLRF,I,J,IELT,NEQU,NEQP
      DOUBLE PRECISION COECON(NNDER,NBLOCF)

C     Offsets if DB and DWORK/KWORK vector

      INTEGER KOFF(2),KVERT,KMID,KCORVG
      
C     Get the vector size of velocity and pressure part:

      NEQU = VECDAT(ONU)
      NEQP = VECDAT(ONP)

C     init the temp. array

      LF(1)=0
      LF(2)=0
      ICLRF=1
      
C     initialize bilinear form
      
      KFN(1)=1
      KFN(2)=1
      
      CALL LCL3(KF,NNAB*NBLOCF)
      KF(1,1)=1
      KF(1,2)=1
      
      CALL LCL1(COECON,NNDER*NBLOCF)
      
C     Initialize offsets
      
      KOFF(1) = 0
      KOFF(2) = NEQU
      KVERT   = L(TRIA(OLVERT))
      KMID    = L(TRIA(OLMID))
      KCORVG  = L(TRIA(OLCORVG))

C     generate right hand side(s)

      CALL LCL1(DRHS,2*NEQU)
      
      IELT = IASMBL(OIELEMT)
      IF (IELT.EQ.0) THEN
        CALL VB0X (TRIA,IPARAM,DPARAM,IGEOM,DGEOM,NBLOCF,DRHS,KOFF,
     *             KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),E031,.FALSE.,
     *             RHS,KFN,KF,BCONF,COECON,IASMBL(OICUBF))
      ELSE IF (IELT.EQ.1) THEN
        CALL VB0X (TRIA,IPARAM,DPARAM,IGEOM,DGEOM,NBLOCF,DRHS,KOFF,
     *             KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),E030,.FALSE.,
     *             RHS,KFN,KF,BCONF,COECON,IASMBL(OICUBF))
      ELSE IF (IELT.EQ.2) THEN
        CALL VB0X (TRIA,IPARAM,DPARAM,IGEOM,DGEOM,NBLOCF,DRHS,KOFF,
     *             KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),EM31,.TRUE.,
     *             RHS,KFN,KF,BCONF,COECON,IASMBL(OICUBF))
      ELSE IF (IELT.EQ.3) THEN
        CALL VB0X (TRIA,IPARAM,DPARAM,IGEOM,DGEOM,NBLOCF,DRHS,KOFF,
     *             KWORK(KVERT),KWORK(KMID),DWORK(KCORVG),EM30,.TRUE.,
     *             RHS,KFN,KF,BCONF,COECON,IASMBL(OICUBF))
      END IF
      
      END 
            
************************************************************************
* Right hand side yielding the exact solution UE
*
* The RHS is implemented by a callback to the user defined function
* FDATIN.
************************************************************************

      DOUBLE PRECISION FUNCTION RHS (X,Y,IA,IBLOC,BFIRST,
     *                               TRIA,IPARAM,DPARAM,IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'sassembly.inc'
      
C     parameters

      DOUBLE PRECISION X,Y
      INTEGER IA, IB, IBLOC
      LOGICAL BFIRST
      
      INTEGER IPARAM(*),TRIA(*),IGEOM(*)
      DOUBLE PRECISION DPARAM(*),DGEOM(*)

C     externals

      DOUBLE PRECISION FDATIN
      EXTERNAL FDATIN
      
C     local variables

      DOUBLE PRECISION TIMENS,RE
      
      RHS=FDATIN(5,IBLOC,X,Y,DPARAM(OTIMENS),DPARAM(ORE),
     *           0,TRIA,IPARAM,DPARAM,IGEOM,DGEOM)

      END
