************************************************************************
* Definition of the coefficients in the matrices, RHS,...
* for the grid adaption routines.
************************************************************************
      
************************************************************************
* Initialise descriptors of bilinear form for iterative
* grid adaption solver
* -> for the system matrix A
* -> for the mass matric M
* -> for the right hand side
*
* This routine initialises the provided KABA-descriptors according to
* the FEAT standard. The coefficients belonging to these descriptors
* are returned by the functions COEAGS/COEMGS.
*
* In: 
*  INSTR  - Bitfield; determines which variables should be initialised:
*           Bit 0 - initialise variables of system matrix
*           Bit 1 - initialise variables of mass matrix
*           Bit 2 - initialise variables of RHS
*           => INSTR=7 -> initialise everything
* Out:
*  BCONx  - whether the coefficients in the integral are constant
*  BSNGLx - whether to work in single precision
*  KABxN  - number of additive terms in integral
*  KABx   - descriptors of bilinear form
************************************************************************

      SUBROUTINE INDCGS (INSTR, BCONA, BCONM, BCONF,
     *                   BSNGLA, BSNGLM, BSNGLF,
     *                   KABAN, KABMN, KFN,
     *                   KABA, KABM, KF)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasicelem.inc'
      
      INTEGER INSTR
      INTEGER KABA(2,NNAB,*),KABM(2,NNAB,*),KF(NNAB,*)
      INTEGER KABAN(*),KABMN(*),KFN(*)
      LOGICAL BCONA(*),BCONM(*),BCONF(*)
      LOGICAL BSNGLA(*),BSNGLM(*),BSNGLF(*)
      
C KABAN - The bilinear form contains a maximum of NNAB=21 terms. KABAN
C describes the *real* number of additive terms in the bilinear form,
C i.e. in the array KABA. 
C
C KABMN - the same for the mass matrix
C
C KFN   - Similarly to KABAN: Number of read additive terms in the
C Integral of the right hand side.

C Coefficients are constant for the matrices, but not 
C for the right hand side.

C All structures are double precision.

C Descriptors of all matrices have to be set properly.

C System matrix:

      IF (IAND(INSTR,1).NE.0) THEN

        KABAN(1)  = 2
        
        BCONA (1) = .TRUE.

        BSNGLA (1) = .FALSE.
        
        KABA(1,1,1)=2
        KABA(2,1,1)=2
        KABA(1,2,1)=3
        KABA(2,2,1)=3
      END IF
      
C Mass matrix:      

      IF (IAND(INSTR,2).NE.0) THEN
        KABMN(1) = 1

        BCONM (1) = .TRUE.

        BSNGLM (1) = .FALSE.
      
        KABM(1,1,1)=1
        KABM(2,1,1)=1
      END IF

C Right hand side

      IF (IAND(INSTR,4).NE.0) THEN
        KFN  (1) = 1
    
        BCONF (1) = .FALSE.

        BSNGLF (1) = .FALSE.

        KF(1,1)=1
      END IF

      END 
      
C *** Coefficient function for the stiffness matrix
      
      DOUBLE PRECISION FUNCTION COEAGS(X,Y,IA,IB,IDA,BFIRST,
     *                                 IPARAM,DPARAM)
      IMPLICIT NONE
      INTEGER IA,IB,IDA,IPARAM(*)
      DOUBLE PRECISION X,Y,DPARAM(*)
      LOGICAL BFIRST
      COEAGS=0D0
      IF ((IA.EQ.2).AND.(IB.EQ.2)) COEAGS=1D0
      IF ((IA.EQ.3).AND.(IB.EQ.3)) COEAGS=1D0
      END
      
C *** Coefficient function for the mass matrix

      DOUBLE PRECISION FUNCTION COEMGS(X,Y,IA,IB,IDA,BFIRST,
     *                                 IPARAM,DPARAM)
      IMPLICIT NONE
      INTEGER IA,IB,IDA,IPARAM(*)
      DOUBLE PRECISION X,Y,DPARAM(*)
      LOGICAL BFIRST
      COEMGS=0D0
      IF ((IA.EQ.1).AND.(IB.EQ.1)) COEMGS=1D0
      END
