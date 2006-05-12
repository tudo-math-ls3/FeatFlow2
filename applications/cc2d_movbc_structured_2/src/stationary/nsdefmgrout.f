************************************************************************
* This file defines the callback-routines for the multigrid solver
* M020 when used in the stationary solver NSDEF.
*
* By convention, the IDATA parameter in the callback routines
* point to the starting address of the TNSDEFILinearSolution
* structure. This structure is prepared before the call of the solver
* to inform the callback routines about the current configuration.
* DDATA currently points to a dummy variable, as no double precision
* data is necessary for the callback routines.
*
* The following routines can be found here; Yxxx-routines are used
* as callback routines of M020:
*
*  SUBROUTINE RESTRU  - Interpolates a vector to a lower level
*
*  SUBROUTINE YAX2    - Matrix-vector-multiplication with NavSt-Mat.
*  SUBROUTINE YDBC2   - Implement Dirichlet boundary conditions
*  SUBROUTINE YEX2    - Solve with a one-level solver
*  SUBROUTINE YPROL2  - Prolongate a vector to a finer level
*  SUBROUTINE  PROLUP - Prolongate (U,V,P); used by YPROL2
*  SUBROUTINE YREST2  - Restrict RHS-vector to a coarser level
*  SUBROUTINE  RESTDP - Restrict defect (DU,DV,DP); used by YREST2
*  SUBROUTINE YSM2    - Smoothing of a solution vector
*  SUBROUTINE YSTEP2  - Step length control
*
*  SUBROUTINE YMG0C   - Preconditioning of BiCGStab with MG
************************************************************************
* Routines for prolongation/restriction of E03x/EM3x can be found in
* MGROUT30.F/MGROUT31.F !
* Routines for prolongation/restriction of E010 (const. pressure) can be
* found in MGROUT10.F !
************************************************************************

************************************************************************
* Prolongate solution vector
* Extended calling convention
*
* This routine prolongates a solution vector DUC on a coarser
* level to a solution vector DUF on a finer level.
*
* The vector must have the form DUx=(DU1,DU2,DUP) for velocity
* and pressure.
*
* In:
*   NUC   - Number of equations in each velocity vector on coarse grid
*   NUF   - Number of equations in each velocity vector on fine grid
*   DUC   - array [1..*] of double
*           Coarse grid solution vector; must have the form (DU1,DU2,DP)
*   TRIAC - array [1..SZTRIA] of integer
*           Coarse grid triangulation structure
*   TRIAF - array [1..SZTRIA] of integer
*           Fine grid triangulation structure
*   IINT  - Type of interpolation for velocity vectors
*           = 1, EM31-prolongation, standard, constant pressure
*           = 2: EM31-prolongation, standard, linear pressure
*           = 3, EM31-prolongation, extended, const pres.
*           = 4: EM31-prolongation, extended, lin. pres.
*           =-1, EM30-prolongation, standard, constant pressure
*           =-2: EM30-prolongation, standard, linear pressure
*           =-3, EM30-prolongation, extended, const pres.
*           =-4: EM30-prolongation, extended, lin. pres.
*   IAPR  - Configures extended prolongation/restriction;
*           Only if |IINT|=3,4. Bitfield.
*           Bit 0: switch to constant prol., depending on DPREP
*           Bit 2: switch prolongation depending on size of neighbour 
*                  element, too (additionally to Bit 0)
*   IAVPR - Type of averaging in extended prolongation/restriction.
*           Only if |IINT|=3,4. Bitfield.
*           0=simple averaging in both, prol. and rest.
*           Bit 0: weighted averaging by element size (L2 projection) in prolongation
*           Bit 2: weighted averaging by element size of neighbour element instead of L2-proj.
*                  in prolongation (additionally to Bit 0)
*   DPREP - Configures extended prolongation/restriction;
*           Only if |IINT|=3,4. Bitfield.
*           Defines an aspect ration when to switch to constant
*           prolongation/restriction
*                  
* Out:
*   DUF   - array [1..*] of double
*           Fine grid solution vector; in the form (DU1,DU2,DP)
************************************************************************
      
      SUBROUTINE PROLUP (NUC,NUF,DUC,DUF,TRIAC,TRIAF,IINT,
     *                   IAPR,IAVPR,DPREP)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'

C parameters

      DOUBLE PRECISION DUC(*),DUF(*)
      INTEGER TRIAC(SZTRIA),TRIAF(SZTRIA)
      INTEGER IAPR,NUC,NUF,IINT,IAVPR
      DOUBLE PRECISION DPREP
      
      INTEGER KCORV1,KVERT1,KMID1,KADJ1
      INTEGER KVERT2,KMID2,KADJ2
      INTEGER KAREA1
      
      INTEGER IPC,IPF

      INTEGER NMT1, NEL1, NVT1
      INTEGER NMT2, NEL2, NVT2
      
C local variables

      INTEGER IADPR1, IAVRT2

C     Fetch the addresses of the necessary arrays from the 
C     triangulation structure - so we don't have to write
C     so much :)

      KCORV1 = L(TRIAC(OLCORVG))
      KVERT1 = L(TRIAC(OLVERT))
      KMID1  = L(TRIAC(OLMID))
      KADJ1  = L(TRIAC(OLADJ))
      KAREA1 = L(TRIAC(OLAREA))

      KVERT2 = L(TRIAF(OLVERT))
      KMID2  = L(TRIAF(OLMID))
      KADJ2  = L(TRIAF(OLADJ))
      
      NMT1 = TRIAC(ONMT)
      NEL1 = TRIAC(ONEL)
      NVT1 = TRIAC(ONVT)

      NMT2 = TRIAF(ONMT)
      NEL2 = TRIAF(ONEL)
      NVT2 = TRIAF(ONVT)

      IPC = 1+2*NUC
      IPF = 1+2*NUF

C     At first interpolate the pressure - because we might need
C     the velocity vector space as temporary space:

      IF ((ABS(IINT).EQ.1).OR.(ABS(IINT).EQ.3)) THEN
      
C       Standard constant interpolation of the pressure:
      
        CALL MP010 (DUC(IPC),DUF(IPF),
     *              KWORK(KADJ1),KWORK(KADJ2),NEL1,NEL2)
      
      ELSE
      
C       Linear interpolation for pressure
C
C       That's a bit tricky now, because we need temporary space
C       for the linearly interpolated pressure. But as we have not yet
C       prolongated the velocity, we can use the two velocity vectors
C       for that purpose - they are later overwritten anyway.
C
C       At first use the routine C2N2DM to convert the constant
C       pressure on the coarse grid to linear pressure in the midpoints
C       on the coarse grid. The result is saved into the first velocity
C       vector on the fine grid

        CALL C2N2DM(DUC(IPC),DUF(1),KWORK(KMID1),
     *              KWORK(KADJ1),NEL1,NMT1,NVT1,0)

C       Then use the linear prolongation routine to prolongate the
C       vector linearly. The result is saved to the second velocity
C       vector:

        IF (IINT.GT.0) THEN
         CALL MP031(DUF(1),DUF(1+NUF),KWORK(KVERT1),KWORK(KVERT2),
     *              KWORK(KMID1),KWORK(KMID2),KWORK(KADJ1),KWORK(KADJ2),
     *              NVT1,NVT2,NEL1,NEL2,NMT2)
        ELSE
         CALL MP030(DUF(1),DUF(1+NUF),KWORK(KVERT1),KWORK(KVERT2),
     *              KWORK(KMID1),KWORK(KMID2),KWORK(KADJ1),KWORK(KADJ2),
     *              NVT1,NVT2,NEL1,NEL2,NMT2)
        ENDIF

C       Convert the pressure back to constant. Grab the information from
C       the second velocity vector and save the constant pressure to
C       its destination.

        CALL C2N2DM(DUF(IPF),DUF(1+NUF),KWORK(KMID2),
     *             KWORK(KADJ2),NEL2,NMT2,NVT2,1)
      
      END IF
      
C     Now we come to the velocity part. That's more or less easy now.
C     First check whether or not to perform standard prolongation:

      IF ((IINT.EQ.1).OR.(IINT.EQ.2)) THEN
        CALL MP031 (DUC,DUF,KWORK(KVERT1),KWORK(KVERT2),
     *                      KWORK(KMID1),KWORK(KMID2),
     *                      KWORK(KADJ1),KWORK(KADJ2),
     *                      NVT1,NVT2,NEL1,NEL2,NMT2)
        CALL MP031 (DUC(1+NUC),DUF(1+NUF),KWORK(KVERT1),KWORK(KVERT2),
     *                      KWORK(KMID1),KWORK(KMID2),
     *                      KWORK(KADJ1),KWORK(KADJ2),
     *                      NVT1,NVT2,NEL1,NEL2,NMT2)
     
      ELSE IF ((IINT.EQ.-1).OR.(IINT.EQ.-2)) THEN
        CALL MP030 (DUC,DUF,KWORK(KVERT1),KWORK(KVERT2),
     *                      KWORK(KMID1),KWORK(KMID2),
     *                      KWORK(KADJ1),KWORK(KADJ2),
     *                      NVT1,NVT2,NEL1,NEL2,NMT2)
        CALL MP030 (DUC(1+NUC),DUF(1+NUF),KWORK(KVERT1),KWORK(KVERT2),
     *                      KWORK(KMID1),KWORK(KMID2),
     *                      KWORK(KADJ1),KWORK(KADJ2),
     *                      NVT1,NVT2,NEL1,NEL2,NMT2)
     
      ELSE IF ((IINT.EQ.3).OR.(IINT.EQ.4)) THEN
      
C       The extended prolongation/restriction with respecting
C       the cell sizes need a little preparation...
C       Evaluate the parameters!
      
C       Evaluate bitfield for prolongation/restriction

        IADPR1=0
        IF(IAND(IAPR,1).EQ.1) IADPR1=1
        IF(IAND(IAPR,5).EQ.5) IADPR1=2
        
C       Evaluate bitfield for weighting

        IAVRT2=0
        IF(IAND(IAVPR,1).EQ.1) IAVRT2=1
        IF(IAND(IAVPR,5).EQ.5) IAVRT2=2
        
C       Prolongate both velocity vectors
        
        CALL MP031X(DUC,DUF,KWORK(KVERT1),KWORK(KVERT2),
     *              KWORK(KMID1),KWORK(KMID2),
     *              KWORK(KADJ1),KWORK(KADJ2),
     *              NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,
     *              DWORK(KCORV1),
     *              DWORK(KAREA1),IAVRT2,DPREP,IADPR1)
        CALL MP031X(DUC(1+NUC),DUF(1+NUF),KWORK(KVERT1),KWORK(KVERT2),
     *              KWORK(KMID1),KWORK(KMID2),
     *              KWORK(KADJ1),KWORK(KADJ2),
     *              NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,
     *              DWORK(KCORV1),
     *              DWORK(KAREA1),IAVRT2,DPREP,IADPR1)
     
      ELSE IF ((IINT.EQ.-3).OR.(IINT.EQ.-4)) THEN
      
C Evaluate bitfield for prolongation/restriction

        IADPR1=0
        IF(IAND(IAPR,1).EQ.1) IADPR1=1
        IF(IAND(IAPR,5).EQ.5) IADPR1=2
        
C Evaluate bitfield for weighting

        IAVRT2=0
        IF(IAND(IAVPR,1).EQ.1) IAVRT2=1
        IF(IAND(IAVPR,5).EQ.5) IAVRT2=2

C       Prolongate both velocity vectors
        
        CALL MP030X(DUC,DUF,KWORK(KVERT1),KWORK(KVERT2),
     *              KWORK(KMID1),KWORK(KMID2),
     *              KWORK(KADJ1),KWORK(KADJ2),
     *              NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,
     *              DWORK(KCORV1),
     *              DWORK(KAREA1),IAVRT2,DPREP,IADPR1)
        CALL MP030X(DUC(1+NUC),DUF(1+NUF),KWORK(KVERT1),KWORK(KVERT2),
     *              KWORK(KMID1),KWORK(KMID2),
     *              KWORK(KADJ1),KWORK(KADJ2),
     *              NVT1,NVT2,NEL1,NEL2,NMT1,NMT2,
     *              DWORK(KCORV1),
     *              DWORK(KAREA1),IAVRT2,DPREP,IADPR1)
      END IF
      
      END 

************************************************************************
* Restrict defect vector
*
* Extended calling convention
*
* This routine restricts a defect vector DDF on a finer
* level to a defect vector DDC on a coarser level.
*
* The vector must have the form DDx=(DD1,DD2,DDP) for velocity
* and pressure.
*
* In:
*   NUC   - Number of equations in each velocity vector on coarse grid
*   NUF   - Number of equations in each velocity vector on fine grid
*   DDF   - array [1..*] of double
*           Fine grid defect vector; must have the form (DD1,DD2,DDP)
*   TRIAC - array [1..SZTRIA] of integer
*           Coarse grid triangulation structure
*   TRIAF - array [1..SZTRIA] of integer
*           Fine grid triangulation structure
*   IINT  - Type of interpolation for velocity vectors
*           = 1: EM31-restriction, standard, constant pressure
*           = 2: EM31-restriction, standard, linear pressure
*           = 3: EM31-restriction, extended, const pres.
*           = 4: EM31-restriction, extended, lin. pres.
*           =-1: EM30-restriction, standard, constant pressure
*           =-2: EM30-restriction, standard, linear pressure
*           =-3: EM30-restriction, extended, const pres.
*           =-4: EM30-restriction, extended, lin. pres.
*   IAPR  - Configures extended prolongation/restriction;
*           Only if |IINT|=3,4. Bitfield.
*           Bit 1: switch to constant rest., depending on DPREP
*           Bit 3: switch restriction depending on size of neighbour 
*                  element, too (additionally to Bit 1)
*   IAVPR - Type of averaging in extended prolongation/restriction.
*           Only if |IINT|=3,4. Bitfield.
*           0=simple averaging in both, prol. and rest.
*           Bit 1: weighted averaging by element size (L2 projection) in restriction
*           Bit 3: weighted averaging by element size of neighbour element instead of L2-proj.
*                  in restriction (additionally to Bit 1)
*   DPREP - Configures extended prolongation/restriction;
*           Only if |IINT|=3,4. Bitfield.
*           Defines an aspect ration when to switch to constant
*           prolongation/restriction
*                  
* Out:
*   DDC   - array [1..*] of double
*           Coarse grid defect vector; in the form (DD1,DD2,DDP)
*
* Remark: For linear interpolation of the pressure, the defect vector
*   in the velocity components of the fine grid is overwritten with
*   auxiliary data. If the caller still needs it (is normally not the
*   case in a MG algorithm), it must save the information elsewhere!
************************************************************************

      SUBROUTINE RESTDP (NUC,NUF,DDC,DDF,TRIAC,TRIAF,IINT,
     *                   IAPR,IAVPR,DPREP)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'

C parameters

      DOUBLE PRECISION DDC(*),DDF(*)
      INTEGER TRIAC(SZTRIA),TRIAF(SZTRIA)
      INTEGER IAPR,IAVPR,NUC,NUF,IINT
      DOUBLE PRECISION DPREP
      
      INTEGER KCORV1,KVERT1,KMID1,KADJ1
      INTEGER KVERT2,KMID2,KADJ2
      INTEGER KAREA1
      
      INTEGER IPC,IPF

      INTEGER NMT1, NEL1, NVT1
      INTEGER NMT2, NEL2, NVT2
      
C local variables

      INTEGER IADPR1, IAVRT2

C     Fetch the addresses of the necessary arrays from the 
C     triangulation structure - so we don't have to write
C     so much :)

      KCORV1 = L(TRIAC(OLCORVG))
      KVERT1 = L(TRIAC(OLVERT))
      KMID1  = L(TRIAC(OLMID))
      KADJ1  = L(TRIAC(OLADJ))
      KAREA1 = L(TRIAC(OLAREA))

      KVERT2 = L(TRIAF(OLVERT))
      KMID2  = L(TRIAF(OLMID))
      KADJ2  = L(TRIAF(OLADJ))
      
      NMT1 = TRIAC(ONMT)
      NEL1 = TRIAC(ONEL)
      NVT1 = TRIAC(ONVT)

      NMT2 = TRIAF(ONMT)
      NEL2 = TRIAF(ONEL)
      NVT2 = TRIAF(ONVT)

      IPC = 1+2*NUC
      IPF = 1+2*NUF

C     At first restrict the velocity parts of the defect vectors.
C     Check whether to restrict with standard or extended restriction.

      IF ((IINT.EQ.1).OR.(IINT.EQ.2)) THEN
        CALL MR031(DDF,DDC,KWORK(KVERT2),KWORK(KVERT1),
     *                     KWORK(KMID2),KWORK(KMID1),
     *                     KWORK(KADJ2),KWORK(KADJ1),
     *                     NVT2,NVT1,NEL2,NEL1)
        CALL MR031(DDF(1+NUF),DDC(1+NUC),KWORK(KVERT2),KWORK(KVERT1),
     *                     KWORK(KMID2),KWORK(KMID1),
     *                     KWORK(KADJ2),KWORK(KADJ1),
     *                     NVT2,NVT1,NEL2,NEL1)
     
      ELSE IF ((IINT.EQ.-1).OR.(IINT.EQ.-2)) THEN
      
        CALL MR030(DDF,DDC,KWORK(KVERT2),KWORK(KVERT1),
     *                     KWORK(KMID2),KWORK(KMID1),
     *                     KWORK(KADJ2),KWORK(KADJ1),
     *                     NVT2,NVT1,NEL2,NEL1)
        CALL MR030(DDF(1+NUF),DDC(1+NUC),KWORK(KVERT2),KWORK(KVERT1),
     *                     KWORK(KMID2),KWORK(KMID1),
     *                     KWORK(KADJ2),KWORK(KADJ1),
     *                     NVT2,NVT1,NEL2,NEL1)
     
      ELSE IF ((IINT.EQ.3).OR.(IINT.EQ.4)) THEN

C       The extended restriction needs a little preparation.
      
C       Evaluate Bitfield for prol/rest

        IADPR1=0
        IF(IAND(IAPR,2).EQ.2) IADPR1=1
        IF(IAND(IAPR,10).EQ.10) IADPR1=2
        
C       Evaluate Bitfield for weighting

        IAVRT2=0
        IF(IAND(IAVPR,2).EQ.2) IAVRT2=1
        IF(IAND(IAVPR,10).EQ.10) IAVRT2=2
        
C       Call the restriction to calculate DUC from DUF
        
        CALL MR031X(DDF,DDC,KWORK(KVERT2),KWORK(KVERT1),
     *                      KWORK(KMID2),KWORK(KMID1),
     *                      KWORK(KADJ2),KWORK(KADJ1),
     *                      NVT2,NVT1,NEL2,NEL1,
     *                      NMT2,NMT1,DWORK(KCORV1),
     *                      DWORK(KAREA1),IAVRT2,DPREP,IADPR1)
        CALL MR031X(DDF(1+NUF),DDC(1+NUC),KWORK(KVERT2),KWORK(KVERT1),
     *                      KWORK(KMID2),KWORK(KMID1),
     *                      KWORK(KADJ2),KWORK(KADJ1),
     *                      NVT2,NVT1,NEL2,NEL1,
     *                      NMT2,NMT1,DWORK(KCORV1),
     *                      DWORK(KAREA1),IAVRT2,DPREP,IADPR1)
     
      ELSE IF ((IINT.EQ.-3).OR.(IINT.EQ.-4)) THEN
      
C Evaluate Bitfield for prol/rest

        IADPR1=0
        IF(IAND(IAPR,2).EQ.2) IADPR1=1
        IF(IAND(IAPR,10).EQ.10) IADPR1=2
        
C Evaluate Bitfield for weighting

        IAVRT2=0
        IF(IAND(IAVPR,2).EQ.2) IAVRT2=1
        IF(IAND(IAVPR,10).EQ.10) IAVRT2=2
        
C       Call the restriction to calculate DUC from DUF

        CALL MR030X(DDF,DDC,KWORK(KVERT2),KWORK(KVERT1),
     *                      KWORK(KMID2),KWORK(KMID1),
     *                      KWORK(KADJ2),KWORK(KADJ1),
     *                      NVT2,NVT1,NEL2,NEL1,
     *                      NMT2,NMT1,DWORK(KCORV1),
     *                      DWORK(KAREA1),IAVRT2,DPREP,IADPR1)
        CALL MR030X(DDF(1+NUF),DDC(1+NUC),KWORK(KVERT2),KWORK(KVERT1),
     *                      KWORK(KMID2),KWORK(KMID1),
     *                      KWORK(KADJ2),KWORK(KADJ1),
     *                      NVT2,NVT1,NEL2,NEL1,
     *                      NMT2,NMT1,DWORK(KCORV1),
     *                      DWORK(KAREA1),IAVRT2,DPREP,IADPR1)
     
      END IF

C     Now after the velocity part has been restricted, we come to
C     the pressure part. Check whether to restrict with constant
C     or linear restriction:

      IF ((ABS(IINT).EQ.1).OR.(ABS(IINT).EQ.3)) THEN

C       Constant restriction. That's easy, simply call the corresponding
C       routine:
      
        CALL MR010 (DDC(IPC),DDF(IPF),KWORK(KADJ1),KWORK(KADJ2),
     *              NEL1,NEL2)
        
      ELSE
      
C       Linear restriction. This is a little bit more tricky. We have
C       to convert the pressure to linear before restricting it.
C       The problem is, this needs some auxiliary memory.
C
C       As auxiliary memory, we use the "old" defect vector on the
C       fine grid. For standard MG algorithm, the old defect vector
C       is no more used - so there's no harm in overwriting it. If
C       the caller still needs it, there must be created a backup
C       for it.
C
C       At first convert the pressure to linear. Save the resulting
C       vector in the first velocity defect vector on the fine grid:

        CALL C2N2DM(DDF(IPF),DDF(1),KWORK(KMID2),KWORK(KADJ2),
     *              NEL2,NMT2,NVT2,0)

C       Now restrict the pressure. Save the result in the second
C       velocity defect vector overwriting the old defect.

        IF (IINT.GT.0) THEN
         CALL MR031(DDF(1),DDF(1+NUF),KWORK(KVERT2),KWORK(KVERT1),
     *              KWORK(KMID2),KWORK(KMID1),KWORK(KADJ2),KWORK(KADJ1),
     *              NVT2,NVT1,NEL2,NEL1)
        ELSE
         CALL MR030(DDF(1),DDF(1+NUF),KWORK(KVERT2),KWORK(KVERT1),
     *              KWORK(KMID2),KWORK(KMID1),KWORK(KADJ2),KWORK(KADJ1),
     *              NVT2,NVT1,NEL2,NEL1)
        ENDIF

C       Convert the restricted vector back to constant pressure on
C       the coarse level.

        CALL C2N2DM(DDC(IPC),DDF(1+NUF),KWORK(KMID1),KWORK(KADJ1),
     *              NEL1,NMT1,NVT1,1)
      
      END IF

      END

************************************************************************
* Restrict solution vector
*
* Extended calling convention
*
* This routine restricts a given solution (not defect!) vector from a
* fine grid to a coarse grid by linear interpolation. Only the velocity
* part is restricted.
*
* The vector must have the form DUx=(DU1,DV1) for the two velocity
* components. The vector is assumed to be assembled with E030,E031,
* EM30 or EM31, respectively.
*
* In:
*   DU2,
*   DV2   - array [1..NMT(fine)] of double precision
*           Fine grid solution vector
*
*   KVERT1,
*   KMID1,
*   KADJ1,
*   NEQ1,
*   NEL1,
*   NVT1  - coarse grid information
*
*   KVERT2,
*   KMID2,
*   KADJ2,
*   NEQ2,
*   NEL2,
*   NVT2  - fine grid information
*
* Out:
*   DU1,
*   DV1   - array [1..NMT(coarse)] of double precision
*           Coarse grid velocity solution vectors
*
* Remark: Works for all finds of boundary conditions.
************************************************************************

      SUBROUTINE RESTRU (DU1,DV1,DU2,DV2,
     *                   KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                   KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2 )

      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      
C constants
      
      DOUBLE PRECISION A1, A2, A3, R1, R2, R3

      PARAMETER (A1=0.1875D0, A2=0.375D0, A3=-0.0625D0)
      PARAMETER (R1=0.375D0, R2=0.75D0, R3=-0.125D0)

C parameters

      DOUBLE PRECISION DU1(*),DV1(*),  DU2(*),DV2(*)
      INTEGER KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*)
      INTEGER KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)

      INTEGER NEQ1, NEL1, NVT1
      INTEGER NEQ2, NEL2, NVT2

C local variables

      INTEGER IEL1, IM1, IM2, IM3, IM4, IELH1, IELH2, IELH3, IELH4
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12
      DOUBLE PRECISION DUH1, DUH2, DUH3, DUH4, DUH5, DUH6, DUH7, DUH8
      DOUBLE PRECISION DUH9, DUH10, DUH11, DUH12
      DOUBLE PRECISION DVH1, DVH2, DVH3, DVH4, DVH5, DVH6, DVH7, DVH8
      DOUBLE PRECISION DVH9, DVH10, DVH11, DVH12

C-----------------------------------------------------------------------
C
      DO IEL1=1,NEL1

        IM1=KMID1(1,IEL1)-NVT1
        IM2=KMID1(2,IEL1)-NVT1
        IM3=KMID1(3,IEL1)-NVT1
        IM4=KMID1(4,IEL1)-NVT1
C
        IELH1=IEL1
        IELH2=KADJ2(2,IELH1)
        IELH3=KADJ2(2,IELH2)
        IELH4=KADJ2(2,IELH3)
C
        I1=KMID2(1,IELH1)-NVT2
        I2=KMID2(4,IELH2)-NVT2
        I3=KMID2(1,IELH2)-NVT2
        I4=KMID2(4,IELH3)-NVT2
        I5=KMID2(1,IELH3)-NVT2
        I6=KMID2(4,IELH4)-NVT2
        I7=KMID2(1,IELH4)-NVT2
        I8=KMID2(4,IELH1)-NVT2
        I9=KMID2(2,IELH1)-NVT2
        I10=KMID2(2,IELH2)-NVT2
        I11=KMID2(2,IELH3)-NVT2
        I12=KMID2(2,IELH4)-NVT2
C
        DUH1= DU2(I1)
        DUH2= DU2(I2)
        DUH3= DU2(I3)
        DUH4= DU2(I4)
        DUH5= DU2(I5)
        DUH6= DU2(I6)
        DUH7= DU2(I7)
        DUH8= DU2(I8)
        DUH9= DU2(I9)
        DUH10=DU2(I10)
        DUH11=DU2(I11)
        DUH12=DU2(I12)
C
        DVH1= DV2(I1)
        DVH2= DV2(I2)
        DVH3= DV2(I3)
        DVH4= DV2(I4)
        DVH5= DV2(I5)
        DVH6= DV2(I6)
        DVH7= DV2(I7)
        DVH8= DV2(I8)
        DVH9= DV2(I9)
        DVH10=DV2(I10)
        DVH11=DV2(I11)
        DVH12=DV2(I12)
C
C ***   The edge IM1
C
        IF (KADJ1(1,IEL1).NE.0) THEN
C       case of an inner edge
         IF (KADJ1(1,IEL1).GT.IEL1) THEN
          DU1(IM1)=A1*(DUH1+DUH2) +A2*DUH9 +A3*(DUH8+DUH3+DUH10+DUH12)
          DV1(IM1)=A1*(DVH1+DVH2) +A2*DVH9 +A3*(DVH8+DVH3+DVH10+DVH12)
         ELSE
          DU1(IM1)=DU1(IM1) +
     *             A1*(DUH1+DUH2) +A2*DUH9 +A3*(DUH8+DUH3+DUH10+DUH12)
          DV1(IM1)=DV1(IM1) +
     *             A1*(DVH1+DVH2) +A2*DVH9 +A3*(DVH8+DVH3+DVH10+DVH12)
         ENDIF
        ELSE
C       case of a boundary edge
          DU1(IM1)=R1*(DUH1+DUH2) +R2*DUH9 +R3*(DUH8+DUH3+DUH10+DUH12)
          DV1(IM1)=R1*(DVH1+DVH2) +R2*DVH9 +R3*(DVH8+DVH3+DVH10+DVH12)
        ENDIF
C
C ***   The edge IM2
C
        IF (KADJ1(2,IEL1).NE.0) THEN
C       case of an inner edge
         IF (KADJ1(2,IEL1).GT.IEL1) THEN
          DU1(IM2)=A1*(DUH3+DUH4) +A2*DUH10 +A3*(DUH2+DUH5+DUH9 +DUH11)
          DV1(IM2)=A1*(DVH3+DVH4) +A2*DVH10 +A3*(DVH2+DVH5+DVH9 +DVH11)
         ELSE
          DU1(IM2)=DU1(IM2) +
     *             A1*(DUH3+DUH4) +A2*DUH10 +A3*(DUH2+DUH5+DUH9 +DUH11)
          DV1(IM2)=DV1(IM2) +
     *             A1*(DVH3+DVH4) +A2*DVH10 +A3*(DVH2+DVH5+DVH9 +DVH11)
         ENDIF
        ELSE
C       case of a boundary edge
          DU1(IM2)=R1*(DUH3+DUH4) +R2*DUH10 +R3*(DUH2+DUH5+DUH9 +DUH11)
          DV1(IM2)=R1*(DVH3+DVH4) +R2*DVH10 +R3*(DVH2+DVH5+DVH9 +DVH11)
        ENDIF
C
C ***   The edge IM3
C
        IF (KADJ1(3,IEL1).NE.0) THEN
C       case of an inner edge
         IF (KADJ1(3,IEL1).GT.IEL1) THEN
          DU1(IM3)=A1*(DUH5+DUH6) +A2*DUH11 +A3*(DUH4+DUH7+DUH10+DUH12)
          DV1(IM3)=A1*(DVH5+DVH6) +A2*DVH11 +A3*(DVH4+DVH7+DVH10+DVH12)
         ELSE
          DU1(IM3)=DU1(IM3) +
     *             A1*(DUH5+DUH6) +A2*DUH11 +A3*(DUH4+DUH7+DUH10+DUH12)
          DV1(IM3)=DV1(IM3) +
     *             A1*(DVH5+DVH6) +A2*DVH11 +A3*(DVH4+DVH7+DVH10+DVH12)
         ENDIF
        ELSE
C       case of a boundary edge
          DU1(IM3)=R1*(DUH5+DUH6) +R2*DUH11 +R3*(DUH4+DUH7+DUH10+DUH12)
          DV1(IM3)=R1*(DVH5+DVH6) +R2*DVH11 +R3*(DVH4+DVH7+DVH10+DVH12)
        ENDIF
C
C ***   The edge IM4
C
        IF (KADJ1(4,IEL1).NE.0) THEN 
C       case of an inner edge
         IF (KADJ1(4,IEL1).GT.IEL1) THEN
          DU1(IM4)=A1*(DUH7+DUH8) +A2*DUH12 +A3*(DUH6+DUH1+DUH9 +DUH11)
          DV1(IM4)=A1*(DVH7+DVH8) +A2*DVH12 +A3*(DVH6+DVH1+DVH9 +DVH11)
         ELSE
          DU1(IM4)=DU1(IM4) +
     *             A1*(DUH7+DUH8) +A2*DUH12 +A3*(DUH6+DUH1+DUH9 +DUH11)
          DV1(IM4)=DV1(IM4) +
     *             A1*(DVH7+DVH8) +A2*DVH12 +A3*(DVH6+DVH1+DVH9 +DVH11)
         ENDIF
        ELSE
C       case of a boundary edge
          DU1(IM4)=R1*(DUH7+DUH8) +R2*DUH12 +R3*(DUH6+DUH1+DUH9 +DUH11)
          DV1(IM4)=R1*(DVH7+DVH8) +R2*DVH12 +R3*(DVH6+DVH1+DVH9 +DVH11)
        ENDIF

      END DO

      END

************************************************************************
* Callback-routine of multigrid:
* Performs a matrix vector multiplication of the form
*                   DAX:= A1*(A*DX) + A2*DAX
* on the current level with the corresponding system matrix.
* Current system dimension: NEQ.
*
* Uses NSDEF-structures to fetch the matrix that is multiplied
* with the vector.
*
* The vectors DX,DAX have the structure  D=(D1,D2,DP).
* The Dirichlet components of the vector are forced to be =0.
*
* In: 
*   DX      - array [1..NEQ] of double
*   DAX     - array [1..NEQ] of double
*   NEQ     - Number of equations
*   A1      - double
*   A2      - double
*
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure of multigrid solver
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure of multigrid solver
*   IDATA   - array [1..*] of integer
*             TNSDEFILinearSolution structure with parameters
*             for the callback routine.
*   DDATA   - array [1..*] of double
*             Dummy variable, currently not used.
*
* Out:
*   DAX = A1*(A*DX) + A2*DAX, Dirichlet components forced to be =0
************************************************************************

      SUBROUTINE YAX2 (DX,DAX,NEQ,A1,A2,IPARAM,DPARAM,IDATA,DDATA)  
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'snsdeflinsol.inc'
      
C parameters

      INTEGER IPARAM(*),IDATA(*)
      DOUBLE PRECISION DPARAM(*),DDATA(*)
      DOUBLE PRECISION DX(*),DAX(*),A1,A2
      INTEGER NEQ

C local variables

      INTEGER ILEV,NU,KSCNPR,NBDMT,IDXTRI,IDXMT
      
C     The current level is saved in the multigrid status.
C     Look at which solver we use to decide where the MG parameter
C     block is to be found in IPARAM!
C     The routine can be called under three circumstances:
C     - Multigrid solver (SLTAG=11)
C     - BiCGStab solver 
C     - Multigrid preconditioner (SLTAG=111)

      IF ((IPARAM(OSLTAG).EQ.11).OR.(IPARAM(OSLTAG).EQ.111)) THEN
        ILEV = IPARAM(OILEV)
      ELSE
        ILEV = IPARAM(SZSLVI+ONLMAX)
      END IF
      
C     Using the level, we can access the correct matrices.
C     By convention (see nsdeflinsol.f), during the multigrid the 
C     handles of the matrices are saved in IMATS as substructure 
C     of the linear solver structure. This structure (described in
C     snsdeflinsol.inc) was prepared by our "parent" NSLINS
C     and transferred to us via IDATA.

      IDXMT = OIMATS+(ILEV-1)*SZN2MI

      NU    = IDATA(IDXMT+ONEQA-1)
      
C     Furthermore, the status part of the structure allowes us to fetch
C     necessary information about our triangulations. Our current
C     TRIA(1,ILEV) structure starts in IDATA behind index

      IDXTRI = OITRIAS-1 + (ILEV-1)*SZTRIA  
      
C     so we can access the shortcut nodal property array:
      
      KSCNPR = L(IDATA(IDXTRI+OTRIUD+1))
      NBDMT  = IDATA(IDXTRI+OTRIUD)
      
C     Then call the matrix multiplication routine to do its work.
      
      CALL MTMUL2(DX,DAX,A1,A2,IDATA(IDXMT))
     
C     Force the Dirichlet boundary components of the velocity
C     defect vector to zero:

      CALL BDRY0 (DAX(1),DAX(1+NU),KWORK(KSCNPR),NBDMT)

      END

************************************************************************
* Callback-routine of multigrid:
* Set Dirichlet boundary components of DX to 0.
*
* In:
*   DX      - array [1..NEQ] of double
*   NEQ     - Number of equations
*
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure of multigrid solver
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure of multigrid solver
*   IDATA   - array [1..*] of integer
*             TNSDEFILinearSolution structure with parameters
*             for the callback routine.
*   DDATA   - array [1..*] of double
*             Dummy variable, currently not used.
*
* Out:
*   Modified DX-vector
************************************************************************
      
      SUBROUTINE YDBC2 (DX,NEQ,IPARAM,DPARAM,IDATA,DDATA)  

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'snsdeflinsol.inc'

C parameters

      INTEGER IPARAM(*),IDATA(*)
      DOUBLE PRECISION DPARAM(*),DDATA(*)

      INTEGER NEQ
      DOUBLE PRECISION DX(NEQ)

C externals

      EXTERNAL BDRY0

C local variables

      INTEGER ILEV,NU,KSCNPR,NBDMT,IDXTRI,IDXMT

C     The current level is saved in the multigrid status:

      ILEV = IPARAM(OILEV)
      
C     Using the level, we can access the correct matrices.
C     By convention (see nsdeflinsol.f), during the multigrid the 
C     handles of the matrices are saved in IMATS as substructure 
C     of the linear solver structure. This structure (described in
C     snsdeflinsol.inc) was prepared by our "parent" NSLINS
C     and transferred to us via IDATA.

      IDXMT = OIMATS+(ILEV-1)*SZN2MI-1

C     The dimensions of the matrices are:

      NU = IDATA(IDXMT+ONEQA)

C     Furthermore, the status part of the structure allowes us to fetch
C     necessary information about our triangulations. Our current
C     TRIA(1,ILEV) structure starts in IDATA behind index

      IDXTRI = OITRIAS + (ILEV-1)*SZTRIA - 1 
      
C     so we can access the shortcut nodal property
      
      KSCNPR = L(IDATA(IDXTRI+OTRIUD+1))
      NBDMT  = IDATA(IDXTRI+OTRIUD)

C     Call the main routine that does the work:

      CALL BDRY0 (DX(1),DX(1+NU),KWORK(KSCNPR),NBDMT)

      END

************************************************************************
* Callback-routine of multigrid:
* Compute on the current level ILEV the solution of A*DX=DB
* with an iterative solver.
*
* DX,DB,DD have the structure  D=(D1,D2,DP)
*
* In:
*   NEQ     - length of arrays
*   DB      - array [1..NEQ] of double
*             RHS-vector
*   DD      - array [1..NEQ] of double
*             Auxiliary vector for the calculation
*
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure of multigrid solver
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure of multigrid solver
*   IDATA   - array [1..*] of integer
*             TNSDEFILinearSolution structure with parameters
*             for the callback routine.
*   DDATA   - array [1..*] of double
*             Dummy variable, currently not used.
*
* Out:
*   RHO     - Convergence rate
*   ITE     - Number of iterations
************************************************************************

      SUBROUTINE YEX2 (DX,DB,DD,NEQ,RHO,ITE,IPARAM,DPARAM,IDATA,DDATA)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'snsdeflinsol.inc'

C     parameters

      INTEGER NEQ,ITE
      DOUBLE PRECISION DX(NEQ),DB(NEQ),DD(NEQ),RHO
      
      INTEGER IPARAM(*),IDATA(*)
      DOUBLE PRECISION DPARAM(*),DDATA(*)

C     local variables

      INTEGER ILEV,NU,NP,IDXTRI,IDXMT
      INTEGER KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB
      INTEGER KMID,KXNPR,KAREA,NEL,NVT

      DOUBLE PRECISION DMAXU,DMAXP,RES,RESINI
      INTEGER NSL,INEUM
      DOUBLE PRECISION EPSSL,DMPSL,RESMAX,RLXSL
      INTEGER LAUX,KBCAUX,I,I2,IP

      INTEGER KU4SYM,KU4NUM,KSYS
      DOUBLE PRECISION INFOU4(90)

C     Set LAUX=0 - will only be changed if we allocate memory...
      
      LAUX = 0

C     The current level is saved in the multigrid status:

      ILEV = IPARAM(OILEV)
      
C     Using the level, we can access the correct matrices.
C     By convention (see nsdeflinsol.f), during the multigrid the 
C     handles of the matrices are saved in IMATS as substructure 
C     of the linear solver structure. This structure (described in
C     snsdeflinsol.inc) was prepared by our "parent" NSLINS
C     and transferred to us via IDATA.

      IDXMT = OIMATS+(ILEV-1)*SZN2MI-1

C     Tell me what the coarse grid solver is...

      IF (IPARAM(SZ020I+OSLTAG).EQ.13) THEN
      
C       We should solve the coarse grid problem by UMFPACK4. Our
C       parent NSLINS prepared the following information for us:
C       - User defined part of integer solver structure (1) 
C         = Handle to symbolic factorization
C       - User defined part of integer solver structure (2) 
C         = Handle to numeric factorization
C       - User defined part of double prec. solver structure (1..20)
C         = Control structure of UMFPACK4

        KU4SYM = IPARAM(SZ020I+OIUDSLA)
        KU4NUM = IPARAM(SZ020I+OIUDSLA+1)

C       Solve the system. Note that UMFPACK expects the matrix in
C       CSR format, which is transposed to our matrix format 7 --
C       so we solve the transposed system:
        
        KSYS = 1
        CALL UMF4SOL(KSYS,DX,DB,
     *               KU4NUM,DPARAM(SZ020D+ODUDSLA),INFOU4)
     
        ITE = 1
        RHO = 0D0
        RESINI = 1D0
        RES = 0D0
        
        IF (INFOU4(1).NE.0D0) THEN
          WRITE (MTERM,'(A,F7.0)') 
     *          'Warning: UMFPACK coarse grid solver. Status=',INFOU4(1)
          RHO = 1D0
          RES = 1D0
        END IF
        
      ELSE  
      
C       Use VANCA as coarse grid solver.
C
C       The dimensions of the matrices are:

        NU = IDATA(IDXMT+ONEQA)
        NP = IDATA(IDXMT+ONEQB)
        
        KA1 = L(IDATA(IDXMT+OLA1))
        KB1 = L(IDATA(IDXMT+OLB1))
        KB2 = L(IDATA(IDXMT+OLB2))
        
        KCOLA = L(IDATA(IDXMT+OLCLA1))
        KCOLB = L(IDATA(IDXMT+OLCLB1))

        KLDA = L(IDATA(IDXMT+OLLDA1))
        KLDB = L(IDATA(IDXMT+OLLDB1))
        
C       Furthermore, the status part of the structure allowes us to fetch
C       necessary information about our triangulations. Our current
C       TRIA(1,ILEV) structure starts in IDATA behind index...

        IDXTRI = OITRIAS + (ILEV-1)*SZTRIA - 1 
        
C       ...so we can access...
        
        KMID  = L(IDATA(IDXTRI+OLMID))
        KXNPR = L(IDATA(IDXTRI+OLXNPR))
        KAREA = L(IDATA(IDXTRI+OLAREA))
        NEL   = IDATA(IDXTRI+ONEL)
        NVT   = IDATA(IDXTRI+ONVT)

C       Theoretically we now have the possibility to choose a solver...
C       But we only have one: VANCA!
C       Vanca is basically a smoother - but remember that every smoother
C       is also a solver!  So we have a solver variant of VANCA at hand,
C       too, which is called VANCE2.
C
C       The second part of the multigrid structure is prepared to hold
C       all standard informations about the coarse grid solver. 
C       We use them to fetch the configuration parameters for our 
C       VANCA-solver:

        NSL    = IPARAM(SZ020I+ONITMAX)
        EPSSL  = DPARAM(SZ020D+OEPSABS)
        DMPSL  = DPARAM(SZ020D+OEPSREL)
        RESMAX = DPARAM(SZ020D+ODIVREL)
        RLXSL  = DPARAM(SZ020D+OOMEGA)

        INEUM = IDATA(OICNEUM)
        
        I2 = 1+NU
        IP = 1+2*NU

        IF ((IPARAM(SZ020I+OSLTAG).LT.50).OR.
     *      (IPARAM(SZ020I+OSLTAG).GT.53)) THEN
          WRITE (MTERM,*) 'Warning: Unknown coarse grid solver! '//
     *                    'Switching to VANCA!'
        END IF

C       For VANCA of type 52/53 we need an auxiliary vector

        IF ((IPARAM(SZ020I+OSLTAG).GE.52).AND.
     *      (IPARAM(SZ020I+OSLTAG).LE.53)) THEN
        
C         Try to use KBCAUX from the MG solver as temporary vector.
C         If it's not available, allocate a new aux. vector.

          KBCAUX = IPARAM(OKCBAUX)
          
          IF (KBCAUX.EQ.0) THEN
          
C           Don't initialise the vector with 0 - we are overwriting it
C           later anyway...
          
            CALL ZNEW(NEQ,-1,LAUX,'DAUX  ')
            KBCAUX = L(LAUX)
            
          END IF
          
C         As our "parent" in NSDEFLINSOL called the multigrid solver with
C         DD=DWORK, KBCAUX is a pointer in DWORK in any case.
        
        END IF

C       Here it goes...

        DO ITE=1,NSL
        
          IF (IPARAM(SZ020I+OSLTAG).EQ.53) THEN
        
C           Solving with the full-block VANCA preconditioner in a 
C           defect correction loop.
C
C           Calculate the residuum vector in DD. Force Dirichlet 
C           components to be =0 by using YAX2.

            CALL LCP1(DB,DD,NEQ)
            CALL YAX2(DX,DD,NEQ,-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)  
           
C           Our parent NSDEFLINSOL.F called the MG solver with DWORK 
C           being the basic address of DD - so KBCAUX is a pointer in 
C           DWORK in every case! Clear the auxiliary vector.

            CALL LCL1(DWORK(KBCAUX),NEQ)
            
C           Call the VANCA diagonal preconditioner; result is written to
C           the auxiliary vector. The vector to be preconditioned is 
C           in DD.

            CALL VANCP5(DWORK(KBCAUX),DWORK(KBCAUX+I2-1),
     *                  DWORK(KBCAUX+IP-1),DD(1),DD(I2),DD(IP),
     *                  DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                  DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *                  NU,NP,KWORK(KMID),
     *                  KWORK(KXNPR),NEL,NVT)
            
C           Incorporate the preconditioned vector into the solution
C           vector. This is the smoothing iteration
C
C             x_new = x_old  +  RLXSL * Vanca(b-Ax_old)

            CALL LLC1 (DWORK(KBCAUX),DX,NEQ,RLXSL,1D0)
            
C           Calculate the maximum norm of the residuum

            CALL LLI1(DWORK(KBCAUX),2*NU,DMAXU,I)
            CALL LLI1(DWORK(KBCAUX+IP-1),NP,DMAXP,I)            

          ELSE IF (IPARAM(SZ020I+OSLTAG).EQ.52) THEN
        
C           Solving with the diagonal VANCA preconditioner in a 
C           defect correction loop.
C
C           Calculate the residuum vector in DD. Force Dirichlet 
C           components to be =0 by using YAX2.

            CALL LCP1(DB,DD,NEQ)
            CALL YAX2(DX,DD,NEQ,-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)  
           
C           Our parent NSDEFLINSOL.F called the MG solver with DWORK 
C           being the basic address of DD - so KBCAUX is a pointer in 
C           DWORK in every case! Clear the auxiliary vector.

            CALL LCL1(DWORK(KBCAUX),NEQ)
            
C           Call the VANCA diagonal preconditioner; result is written to
C           the auxiliary vector. The vector to be preconditioned is 
C           in DD.

            CALL VANCP4(DWORK(KBCAUX),DWORK(KBCAUX+I2-1),
     *                  DWORK(KBCAUX+IP-1),DD(1),DD(I2),DD(IP),
     *                  DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                  DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *                  NU,NP,KWORK(KMID),
     *                  KWORK(KXNPR),NEL,NVT)
            
C           Incorporate the preconditioned vector into the solution
C           vector. This is the smoothing iteration
C
C             x_new = x_old  +  RLXSL * Vanca(b-Ax_old)

            CALL LLC1 (DWORK(KBCAUX),DX,NEQ,RLXSL,1D0)
            
C           Calculate the maximum norm of the residuum

            CALL LLI1(DWORK(KBCAUX),2*NU,DMAXU,I)
            CALL LLI1(DWORK(KBCAUX+IP-1),NP,DMAXP,I)            

          ELSE IF (IPARAM(SZ020I+OSLTAG).EQ.51) THEN
          
C           Solving with extended VANCA

            CALL VANCE3(DX(1),DX(I2),DX(IP),
     *                   DD(1),DD(I2),DD(IP),
     *                   DB(1),DB(I2),DB(IP),
     *                   DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                   DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *                   NU,NP,
     *                   KWORK(KMID),KWORK(KXNPR),NEL,NVT,
     *                   RLXSL,DMAXU,DMAXP)

          ELSE
          
C           Solving with simple VANCA
        
            CALL VANCE2(DX(1),DX(I2),DX(IP),
     *                   DD(1),DD(I2),DD(IP),
     *                   DB(1),DB(I2),DB(IP),
     *                   DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                   DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *                   NU,NP,
     *                   KWORK(KMID),KWORK(KXNPR),NEL,NVT,
     *                   RLXSL,DMAXU,DMAXP)
     
          END IF

C         Bring the pressure to L2_0 in case we have pure Dirichlet-
C         problem (e.g. Driven Cavity) - because then, the pressure has so 
C         solve an indefinite Neumann-problem, and by filtering it to L2_0
C         we make the problem definite again!

          IF (INEUM.EQ.0) CALL TOL20A(DX(1+2*NU),DWORK(KAREA),NP)

          RES=MAX(DMAXU,DMAXP)
          IF (ITE.EQ.1) RESINI=RES
          RHO=(RES/RESINI)**(1D0/DBLE(ITE))
     
          IF (IPARAM(OMCGTRM).GE.2) THEN
            WRITE (MTERM,'(A,I7,A,D25.16)') 
     *            'VANCA: Iteration ',ITE,',  !!RES!! = ',RES
          END IF

          IF ((RES.LT.EPSSL).AND.(RES.LT.DMPSL*RESINI)) GOTO 99999
           
          IF (RES.GT.RESMAX*RESINI) THEN
            WRITE (MTERM,*) 
     *            'Warning: coarse grid solver not convergent!'
            GOTO 99999
          END IF

        END DO
      
      END IF

99999 CONTINUE

C     Correction of ITE if DO-loop completely runs through

      IF (ITE.GT.NSL) ITE=NSL
      
C     In case we allocated memory, release it.

      IF (LAUX.NE.0) CALL ZDISP(0,LAUX,'DAUX  ')

C     Correct RHO on error - don't use RHO>1D99 since this works not
C     with NaN!
      
      IF (.NOT.(RHO.LE.1D99)) RHO=1D0
      
      IF (IPARAM(OMCGTRM).GE.2) THEN
        WRITE (MTERM,'(A)') 'Coarse grid solver statistics:'
        WRITE (MTERM,'(A)') ''
        WRITE (MTERM,'(A,I5)')     'Iterations              : ',ITE
        WRITE (MTERM,'(A,D24.12)') '!!INITIAL RES!!         : ',RESINI
        WRITE (MTERM,'(A,D24.12)') '!!RES!!                 : ',RES
        IF (DPARAM(ODEFINI).GT.DPARAM(OVECZER)) THEN     
          WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
     *          RES / RESINI
        ELSE
          WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',0D0
        END IF
        WRITE (MTERM,'(A,D24.12)') 'Rate of convergence     : ',RHO
        WRITE (MTERM,'(A)') ''
      END IF

      IF (IPARAM(OMCGTRM).EQ.1) THEN
        WRITE (MTERM,'(A,I5,A,D24.12)') 
     *        'Coarse grid slv.: Iterations/Rate of conv.: ',
     *        ITE,' /',RHO
      END IF

      END

************************************************************************
* Callback-routine of multigrid:
* Prolongation of solution vector
*
* Extended calling convention. 
*
* In:
*   DUC     - array [1..NEQ-coarse] of double
*             coarse grid vector, level ILEV-1, to be prolongated
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure for multigrid solver
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure for multigrid solver
*   IDATA   - array [1..*] of integer
*             TNSDEFILinearSolution structure with parameters
*             for the callback routine.
*   DDATA   - array [1..*] of double
*             Dummy variable, currently not used.
* Out:
*   DUF     - array [1..NEQ-fine] of double
*             fine grid vector, level ILEV
*
* Remarks:  - DUF and DUC have the structure  DU=(DU1,DU2,DP)
*  
************************************************************************

      SUBROUTINE YPROL2 (DUC,DUF,IPARAM,DPARAM,IDATA,DDATA)  

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'snsdeflinsol.inc'
      
C parameters

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
      DOUBLE PRECISION DUF(*),DUC(*)

C local variables
      
      INTEGER IDXTRI, IDXMT, ILEV
      
      INTEGER NUF, NUC
      
C     The current level is saved in the multigrid status.
C     We want to prolongate from level ILEV-1 to level ILEV.

      ILEV = IPARAM(OILEV)
      
C     Using the level, we can access the correct matrices.
C     By convention (see nsdeflinsol.f), during the multigrid the 
C     handles of the matrices are saved in IMATS as substructure 
C     of the linear solver structure. This structure (described in
C     snsdeflinsol.inc) was prepared by our "parent" NSLINS
C     and transferred to us via IDATA.

      IDXMT = OIMATS+(ILEV-1)*SZN2MI-1

C     The dimensions of the matrix on the coarser grid is at:

      NUC = IDATA(IDXMT-SZN2MI+ONEQA)
      
C     and that of the current grid can be found at

      NUF = IDATA(IDXMT+ONEQA)

C     Furthermore, the status part of the structure allowes us to fetch
C     the triangulations. TRIA(1,ILEV-1) structure starts in IDATA 
C     behind index...

      IDXTRI = OITRIAS + (ILEV-1)*SZTRIA 
      
C     and TRIA(1,ILEV-1) follows in a distance of SZTRIA.
      
C     Call PROLUP. This will prolongate our velocity- and pressure-
C     vectors. For the configuration of the prolongation, fetch the
C     parameters from the MG solver structure. The initialization
C     routine configured the user defined parameters in these
C     structures as follows:
C
C         PRRSI(1) = IINT
C         PRRSI(2) = IAPR
C         PRRSI(3) = IAVPR
C
C         PRRSD(1) = DPREP
      
      CALL PROLUP (NUC,NUF,DUC,DUF,IDATA(IDXTRI-SZTRIA),IDATA(IDXTRI),
     *             IPARAM(OPRRSI),IPARAM(OPRRSI+1),IPARAM(OPRRSI+2),
     *             DPARAM(OPRRSD))

      END

************************************************************************
* Callback-routine of multigrid:
* Restriction of defect vector
*
* Extended calling convention. 
*
* In:
*   DDF     - array [1..NEQ-coarse] of double
*             fine grid defect vector, level ILEV+1, to be restricted
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure for multigrid solver
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure for multigrid solver
*   IDATA   - array [1..*] of integer
*             TNSDEFILinearSolution structure with parameters
*             for the callback routine.
*   DDATA   - array [1..*] of double
*             Dummy variable, currently not used.
* Out:
*   DDC     - array [1..NEQ-fine] of double
*             coarse grid defect vector, level ILEV
*
* Remarks:  - DDF and DUC have the structure  DD=(DD1,DD2,DDP)
*  
************************************************************************

      SUBROUTINE YREST2 (DDF,DDC,IPARAM,DPARAM,IDATA,DDATA)  

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'snsdeflinsol.inc'
      
C parameters

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
      DOUBLE PRECISION DDF(*),DDC(*)

C local variables
      
      INTEGER IDXTRI, IDXMT, ILEV
      
      INTEGER NUF, NUC
      
C     The current level is saved in the multigrid status.
C     We want to restrict a solution on level ILEV+1 to the current
C     level ILEV.

      ILEV = IPARAM(OILEV)
      
C     Using the level, we can access the correct matrices.
C     By convention (see nsdeflinsol.f), during the multigrid the 
C     handles of the matrices are saved in IMATS as substructure 
C     of the linear solver structure. This structure (described in
C     snsdeflinsol.inc) was prepared by our "parent" NSLINS
C     and transferred to us via IDATA.
C
C     Go to the finer level.

      IDXMT = OIMATS+(ILEV-1)*SZN2MI-1

C     The dimensions of the matrix on the finer level is:

      NUF = IDATA(IDXMT+SZN2MI+ONEQA)
      
C     And the vector length on the coarser grid can be accessed via...

      NUC = IDATA(IDXMT+ONEQA)

C     Furthermore, the status part of the structure allowes us to fetch
C     the triangulations. TRIA(1,ILEV) structure starts in IDATA 
C     behind index...

      IDXTRI = OITRIAS + (ILEV-1)*SZTRIA 
      
C     and TRIA(1,ILEV+1) of the fine grid follows in a distance of SZTRIA.
      
C     Call RESTDP. This will prolongate our velocity- and pressure-
C     vectors. For the configuration of the prolongation, fetch the
C     parameters from the MG solver structure. The initialization
C     routine configured the user defined parameters in these
C     structures as follows:
C
C         PRRSI(1) = IINT
C         PRRSI(2) = IAPR
C         PRRSI(3) = IAVPR
C
C         PRRSD(1) = DPREP
      
      CALL RESTDP (NUC,NUF,DDC,DDF,IDATA(IDXTRI),IDATA(IDXTRI+SZTRIA),
     *             IPARAM(OPRRSI),IPARAM(OPRRSI+1),IPARAM(OPRRSI+2),
     *             DPARAM(OPRRSD))

      END

************************************************************************
* Callback-routine of multigrid:
* Perform NSM smoothing steps to system A*DX=DB
*
* In:
*   DX     - array [1..NEQ] of double;
*            vector to smooth; must have the form (D1,D2,DP)
*   DB     - array [1..NEQ] of double;
*            RHS of current system
*   DD     - array [1..NEQ] of double;
*            auxilary vector
*   NEQ    - length of vectors
*   NSM    - number of smoothing steps
*   IPARAM - array [1..*] of integer
*            Integer parameter structure for multigrid solver
*   DPARAM - array [1..*] of integer
*            Double precision parameter structure for multigrid solver
*   IDATA  - array [1..*] of integer
*            User defined data array - not used.
*   DDATA  - array [1..*] of double
*            User defined data array - not used.
*
* Out:
*   DX      - the smoothed vector
************************************************************************

      SUBROUTINE YSM2 (DX,DB,DD,NEQ,NSMS,IPARAM,DPARAM,IDATA,DDATA)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'snsdeflinsol.inc'

C     parameters

      INTEGER NEQ, NSMS
      DOUBLE PRECISION DX(NEQ),DB(NEQ),DD(NEQ)

      INTEGER IPARAM(*),IDATA(*)
      DOUBLE PRECISION DPARAM(*),DDATA(*)
      
C local variables
      
      INTEGER ILEV,NU,NP,IDXTRI,IDXMT,NEL,NVT
      INTEGER KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB
      INTEGER KMID,KXNPR,KAREA
      INTEGER INEUM
      INTEGER I2,IP,ITE
      INTEGER KBCAUX,LAUX

C     The current level is saved in the multigrid status:

      ILEV = IPARAM(OILEV)
      
C     Using the level, we can access the correct matrices.
C     By convention (see nsdeflinsol.f), during the multigrid the 
C     handles of the matrices are saved in IMATS as substructure 
C     of the linear solver structure. This structure (described in
C     snsdeflinsol.inc) was prepared by our "parent" NSLINS
C     and transferred to us via IDATA.

      IDXMT = OIMATS+(ILEV-1)*SZN2MI-1

C     The dimensions of the matrices are:

      NU = IDATA(IDXMT+ONEQA)
      NP = IDATA(IDXMT+ONEQB)
      
      KA1 = L(IDATA(IDXMT+OLA1))
      KB1 = L(IDATA(IDXMT+OLB1))
      KB2 = L(IDATA(IDXMT+OLB2))
      
      KCOLA = L(IDATA(IDXMT+OLCLA1))
      KCOLB = L(IDATA(IDXMT+OLCLB1))

      KLDA = L(IDATA(IDXMT+OLLDA1))
      KLDB = L(IDATA(IDXMT+OLLDB1))
      
C     Furthermore, the status part of the structure allowes us to fetch
C     necessary information about our triangulations. Our current
C     TRIA(1,ILEV) structure starts in IDATA behind index...

      IDXTRI = OITRIAS + (ILEV-1)*SZTRIA - 1 
      
C     ...so we can access...
      
      KMID  = L(IDATA(IDXTRI+OLMID))
      KXNPR = L(IDATA(IDXTRI+OLXNPR))
      KAREA = L(IDATA(IDXTRI+OLAREA))
      NEL   = IDATA(IDXTRI+ONEL)
      NVT   = IDATA(IDXTRI+ONVT)

C     From the NSDEF-data, get some additional information about 
C     the problem

      INEUM = IDATA(OICNEUM)

      I2 = 1+NU
      IP = 1+2*NU

C     Now perform a small check if we really should smooth with VANCA -
C     we don't have an alternative yet!

      IF ((IPARAM(OSMTAG).LT.50).OR.(IPARAM(OSMTAG).GT.53)) THEN
        WRITE (MTERM,*) 'Warning: Unknown smoother! '//
     *                  'Switching to VANCA!'
      END IF
      
C     For smoothers of type 52/53 we need an auxiliary vector

      LAUX = 0

      IF ((IPARAM(OSMTAG).GE.52).AND.(IPARAM(OSMTAG).LE.53)) THEN
      
C       Try to use KBCAUX from the MG solver as temporary vector.
C       If it's not available, allocate a new aux. vector.

        KBCAUX = IPARAM(OKCBAUX)
        
        IF (KBCAUX.EQ.0) THEN
        
C         Don't initialise the vector with 0 - we are overwriting it
C         later anyway...
        
          CALL ZNEW(NEQ,-1,LAUX,'DAUX  ')
          KBCAUX = L(LAUX)
          
        END IF
        
C       As our "parent" in NSDEFLINSOL called the multigrid solver with
C       DD=DWORK, KBCAUX is a pointer in DWORK in any case.
      
      END IF

C     Implement the smoothing directly.
C     Apply the VANCAS smoother NSMS times onto our solution vector -
C     that's it.
                                              
      DO ITE=1,NSMS

        IF (IPARAM(OSMTAG).EQ.53) THEN
        
C         Smoothing with the VANCA preconditioner.

C         Our parent NSDEFLINSOL.F called the MG solver with DWORK being 
C         the basic address of DD - so KBCAUX is a pointer in DWORK
C         in every case! Clear the auxiliary vector.

          CALL LCL1(DWORK(KBCAUX),NEQ)
          
C         Build the residuum vector in DD. Force Dirichlet components
C         to be =0 by using YAX2.

          CALL LCP1(DB,DD,NEQ)
          CALL YAX2 (DX,DD,NEQ,-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)  
          
C         Call the VANCA full-block preconditioner; result is written to
C         the auxiliary vector. The vector to be preconditioned is in DD.

          CALL VANCP5(DWORK(KBCAUX),DWORK(KBCAUX+I2-1),
     *                DWORK(KBCAUX+IP-1),DD(1),DD(I2),DD(IP),
     *                DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *                NU,NP,KWORK(KMID),
     *                KWORK(KXNPR),NEL,NVT)
          
C         Incorporate the preconditioned vector into the solution
C         vector. This is the smoothing iteration
C
C           x_new = x_old  +  RLXSM * Vanca(b-Ax_old)

          CALL LLC1 (DWORK(KBCAUX),DX,NEQ,DPARAM(OOMGSM),1D0)

        ELSE IF (IPARAM(OSMTAG).EQ.52) THEN
        
C         Smoothing with the VANCA preconditioner.

C         Our parent NSDEFLINSOL.F called the MG solver with DWORK being 
C         the basic address of DD - so KBCAUX is a pointer in DWORK
C         in every case! Clear the auxiliary vector.

          CALL LCL1(DWORK(KBCAUX),NEQ)
          
C         Build the residuum vector in DD. Force Dirichlet components
C         to be =0 by using YAX2.

          CALL LCP1(DB,DD,NEQ)
          CALL YAX2 (DX,DD,NEQ,-1D0,1D0,IPARAM,DPARAM,IDATA,DDATA)  
          
C         Call the VANCA diagonal preconditioner; result is written to
C         the auxiliary vector. The vector to be preconditioned is in DD.

          CALL VANCP4(DWORK(KBCAUX),DWORK(KBCAUX+I2-1),
     *                DWORK(KBCAUX+IP-1),DD(1),DD(I2),DD(IP),
     *                DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *                NU,NP,KWORK(KMID),
     *                KWORK(KXNPR),NEL,NVT)
          
C         Incorporate the preconditioned vector into the solution
C         vector. This is the smoothing iteration
C
C           x_new = x_old  +  RLXSM * Vanca(b-Ax_old)

          CALL LLC1 (DWORK(KBCAUX),DX,NEQ,DPARAM(OOMGSM),1D0)

        ELSE IF (IPARAM(OSMTAG).EQ.51) THEN

C         smoothing with extended VANCA

          CALL VANCS3(DX(1),DX(I2),DX(IP),DD(1),DD(I2),DD(IP),
     *                DB(1),DB(I2),DB(IP),
     *                DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *                NU,NP,KWORK(KMID),
     *                KWORK(KXNPR),NEL,NVT,DPARAM(OOMGSM))

        ELSE
        
C         Standard VANCA-smoothing

          CALL VANCS2(DX(1),DX(I2),DX(IP),DD(1),DD(I2),DD(IP),
     *                DB(1),DB(I2),DB(IP),
     *                DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *                NU,NP,KWORK(KMID),
     *                KWORK(KXNPR),NEL,NVT,DPARAM(OOMGSM))
     
        END IF

C       Bring the pressure to L2_0 in case we have pure Dirichlet-
C       problem (e.g. Driven Cavity) - because then, the pressure has so 
C       solve an indefinite Neumann-problem, and by filtering it to L2_0
C       we make the problem definite again!

        IF (INEUM.EQ.0) CALL TOL20A(DX(IP),DWORK(KAREA),NP)

      END DO
      
C     In case we allocated memory, release it.

      IF (LAUX.NE.0) CALL ZDISP(0,LAUX,'DAUX  ')

      END

************************************************************************
* Step length control
*
* Extended calling convention. 
*
* In:
*   DX      - array [1..NEQ] of double
*             Previous solution vector, finest grid; must have the form
*             (DU1,DU2,DP)
*   DD      - array [1..NEQ] of double
*             fine correction vector on level ILEV, finest grid
*   DB      - array [1..NEQ] of double
*             Right hand side vector, finest grid
*   NEQ     - length of fine-grid vectors
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure for multigrid solver
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure for multigrid solver
*   IDATA   - array [1..*] of integer
*             User defined data array - not used.
*   DDATA   - array [1..*] of double
*             User defined data array - not used.
* Out:
*   ALPHA   - Relaxation parameter for coarse grid correction,
*             according to some optimization criterion.
*             Will be in the range AMINMG/AMAXMG as configured
*             in the multigrid structure.
************************************************************************

      SUBROUTINE YSTEP2 (DX,DD,DB,NEQ,ALPHA,IPARAM,DPARAM,IDATA,DDATA)
      
      IMPLICIT NONE
      
C standard COMMON blocks
      
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'snsdeflinsol.inc'

C parameters

      INTEGER NEQ
      DOUBLE PRECISION DX(NEQ),DB(NEQ),DD(NEQ)
      DOUBLE PRECISION ALPHA
      INTEGER IPARAM(*),IDATA(*)
      DOUBLE PRECISION DPARAM(*),DDATA(*)

C local variables

      INTEGER ILEV, IEQ, NP, KCBAUX, LAUX
      DOUBLE PRECISION DBX, DBY
      DOUBLE PRECISION AMINMG,AMAXMG

C     Fetch some data from the multigrid structure
      
      AMINMG = DPARAM(OSTPMIN)
      AMAXMG = DPARAM(OSTPMAX)
      
C     We also need the size of the pressure vector. Determine the 
C     current level...

      ILEV = IPARAM(OILEV)

C     And get it from the IMATS substructures of the NSDEF structure,
C     as it corresponds to the number of equations in the B-matrices

      NP = IDATA( OIMATS-1 + (ILEV-1)*SZN2MI + ONEQB)

C     There's nothing to do if the step length is prescribed as fixed
      
      IF (AMINMG.EQ.AMAXMG) THEN
        ALPHA=AMINMG
        RETURN
      ENDIF
      
C     Calculation of optimal ALPHA. For the calculation we use
C     the auxiliary array KCBAUX that was provided by the multigrid
C     solver as temporary space.
C
C     Our parent NSDEFLINSOL.F called the MG solver with DWORK being 
C     the basic address of DD - so KBCAUX is a pointer in DWORK!

      KCBAUX = IPARAM(OKCBAUX)
      
C     If we are on the maximum level, we are not allowed to use that
C     temporary array - would destroy intermediate data of the multigrid
C     (see also m020.inc about KBCAUX). Only in that situation,
C     allocate an auxiliary vector and use that for computation.
C     We don't do that on lower levels, because there it's much faster
C     if we use a preallocated array we are allowed to use...

      LAUX = 0
      IF (ILEV.EQ.IPARAM(ONLMAX)) THEN
      
C       Don't initialize it - we are immediately overwriting it...
      
        CALL ZNEW (NEQ,-1,LAUX,'DAUXMG ')
        
        IF (IER.NE.0) THEN
        
C         If nothing helps, 1.0 might be a good guess...
        
          WRITE (MTERM,'(A)') 'Warning: not enough space for step '//
     *                        'length control!'
          ALPHA = 1D0
          RETURN
          
        END IF
        
        KCBAUX = L(LAUX)
        
      END IF
      
C     Calculate nominator of the fraction
      
      CALL LCP1 (DB,DWORK(KCBAUX),NEQ)
      CALL YAX2 (DX,DWORK(KCBAUX),NEQ,-1D0,1D0,
     *           IPARAM,DPARAM,IDATA,DDATA)
     
      DO IEQ=NEQ-1,NEQ-NP,-1
        DWORK(KCBAUX+IEQ) = -DWORK(KCBAUX+IEQ)
      END DO
      
      CALL LSP1(DD,DWORK(KCBAUX),NEQ,DBY)
      
C     Calculate denominator of the fraction
      
      CALL YAX2(DD,DWORK(KCBAUX),NEQ,1D0,0D0,
     *          IPARAM,DPARAM,IDATA,DDATA)
     
      DO IEQ=NEQ,NEQ-NP+1,-1
        DWORK(KCBAUX+IEQ-1)=-DWORK(KCBAUX+IEQ-1)
      END DO
      
      CALL LSP1(DD,DWORK(KCBAUX),NEQ,DBX)

C     Divide them, that's it.
C     Use ALPHA=0 in case DBX is 0; max happen if the coarse
C     grid contains only dirichlet-values on the boundary
C     and no inner vertices (like in the QUAD mesh).

      IF (DBX.NE.0D0) THEN
        ALPHA=DBY/DBX
      ELSE
        ALPHA=1D0
      END IF
      IF (ALPHA.LT.AMINMG) ALPHA=AMINMG
      IF (ALPHA.GT.AMAXMG) ALPHA=AMAXMG

C     Finally, on the highest level, release the temporary array

      IF (LAUX.NE.0) CALL ZDISP (0,LAUX,'DAUXMG ')

      END

************************************************************************
* Preconditioning with multigrid.
*
* This routine is called in an iterative solver (BiCGStab) to
* precondition a vector by using a Multigrid sweep.
*
* In:
*   DG     - array [1..NEQ] of double
*            Vector to be preconditioned
*   NEQ    - Length of the vector
*   IPARAM - array [1..SZMGRI] of integer
*            Integer parameter structure of the (BiCGStab-) solver
*   DPARAM - array [1..SZMGRI] of integer
*            Double precision parameter structure of the (BiCGStab-) 
*            solver
*   IDATA  - array [1..*] of integer
*            User defined integer array
*   DDATA  - array [1..*] of double
*            User defined double array
************************************************************************

      SUBROUTINE YMG0C (DG,NEQ,IPARAM,DPARAM,IDATA,DDATA)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'ssolvers.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'm020.inc'
      
C     parameters
      
      INTEGER NEQ, IPARAM(*), IDATA(*)
      DOUBLE PRECISION DG(NEQ),DPARAM(*),DDATA(*)
      
C     local variables

      INTEGER NLMAX,KRHS,KSOL

C     externals

      EXTERNAL YAX2,YPROL2,YREST2,YSM2,YEX2,YDBC2,YSTEP2,I000

C     Preconditioning means:   DG_new = MG^-1 DG
C                         <=>  MG DG_new = DG
C
C     so copy the vector DG to the RHS of the multigrid structure,
C     solve and return the solution as new DG!
C
C     IPARAM/DPARAM is the solver structure. The multigrid substrucure
C     is to be found in the second block, so at index position
C     SZSLVx of each of these arrays.
C
C     What's the maximum/current level?

      NLMAX = IPARAM(SZSLVI+ONLMAX)
      
C     Get the position of the RHS vector and the solution vector 
C     on the finest level from the
C     KOFFB index pointer in the MG structure.

      KRHS = IPARAM(SZSLVI+OKOFFB+NLMAX-1)+1
      KSOL = IPARAM(SZSLVI+OKOFFX+NLMAX-1)+1
      
C     Copy DG to the RHS and the start vector.

      CALL LCP1(DG,DWORK(KRHS),NEQ)
      CALL LCP1(DG,DWORK(KSOL),NEQ)
      
C     Call MG to solve the system - or at least to make a bounded
C     number of iterations. Give MG the MG structure, placed at
C     the second block in the list of solver structures in IPARAM/
C     DPARAM.

      CALL M020 (IPARAM(SZSLVI+1),DPARAM(SZSLVD+1), IDATA,DDATA, 
     *           DWORK(1),DWORK(1),DWORK(1),
     *           YAX2,YPROL2,YREST2,YSM2,YSM2,
     *           YEX2,YEX2,YDBC2,YSTEP2,I000)
    
C     Fetch back the solution to DG.

      CALL LCP1(DWORK(KSOL),DG,NEQ)    

C     Out "parent" NSLINS expects the user defined area of the
C     BiCGStab solver to be filled with the timing information
C     of the MG solver. Add the time that MG needed to there.

      DPARAM(ODUDSLA+1) = DPARAM(ODUDSLA+1) + DPARAM(SZSLVD+OTMTOT )
      DPARAM(ODUDSLA+2) = DPARAM(ODUDSLA+2) + DPARAM(SZSLVD+OTMPROL)
      DPARAM(ODUDSLA+3) = DPARAM(ODUDSLA+3) + DPARAM(SZSLVD+OTMREST)
      DPARAM(ODUDSLA+4) = DPARAM(ODUDSLA+4) + DPARAM(SZSLVD+OTMDEF )
      DPARAM(ODUDSLA+5) = DPARAM(ODUDSLA+5) + DPARAM(SZSLVD+OTMSMTH)
      DPARAM(ODUDSLA+6) = DPARAM(ODUDSLA+6) + DPARAM(SZSLVD+OTMCGSL)
      DPARAM(ODUDSLA+7) = DPARAM(ODUDSLA+7) + DPARAM(SZSLVD+OTMBC  )
      DPARAM(ODUDSLA+8) = DPARAM(ODUDSLA+8) + DPARAM(SZSLVD+OTMCGC )

C     That's it.
      
      END
      