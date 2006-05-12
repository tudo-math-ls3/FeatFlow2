************************************************************************
* This file contains (callback) routines for the main routines of
* CC2D. They are used by the main program only.
************************************************************************

************************************************************************
* Callback-routine: Proceed to next time step
*
* This routine is a replacement for the default routine for proceeding
* to the next time step in a nonstady simulation. 
* Basically, this routine directly calls the standard routine. 
* In contrast to the default routine the geometry information of the 
* fictitious boundary routines are updated to the current situation.
************************************************************************
      
      SUBROUTINE FBPNTS (NLMIN,NLMAX,
     *                   TRIAS,MATDAT,VECDAT,
     *                   IPARAM,DPARAM,ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *                   NUVP,DUP,DRHS,DAUX,
     *                   TIMEST, TIMENS, TSTEP)
      
      IMPLICIT NONE

      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'smat2dns.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INCLUDE 'ssolvers.inc'
      INCLUDE 'stiming.inc'
      INCLUDE 'snonstead.inc'
      
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'stracking.inc'
      
C parameters

      INTEGER NUVP
      INTEGER IPARAM (*),IMGPAR(*),ISTPAR(*),IASMBL(*),IGEOM(*)
      INTEGER TRIAS(SZTRIA,NNLEV)
      INTEGER MATDAT(SZN2MI,NNLEV),VECDAT(SZN2VI,NNLEV),NLMIN,NLMAX
      DOUBLE PRECISION DPARAM (*),DMGPAR(*),DSTPAR(*),DASMBL(*)
      DOUBLE PRECISION DGEOM(*)
      DOUBLE PRECISION TIMEST, TIMENS, TSTEP
      
      DOUBLE PRECISION DUP(*),DRHS(*),DAUX(*)
      
      INCLUDE 'sinigeometry.inc'
      
C     Set the time of the fictitious boudary routines to the
C     TIMENS+TSTEP. The current simulation time is TIMENS and
C     will be increased to TIMENS+TSTEP in the next time step.
C     In the very first time step, there is TSTEP=0.

      IF (MOD(IPARAM(OCRITE),2).EQ.0) THEN
        DGEOM(OGEOTIM) = TIMENS+TSTEP
      END IF

C     Call the default routine.

      CALL DFPRNT (NLMIN,NLMAX,
     *             TRIAS,MATDAT,VECDAT,
     *             IPARAM,DPARAM,ISTPAR,DSTPAR,
     *             IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *             NUVP,DUP,DRHS,DAUX,
     *             TIMEST, TIMENS, TSTEP)

      END
      
************************************************************************
* Callback-routine: grid adaption.
*
* This routine is called by the nonstationary solver to adapt
* the grid in each time step.
*
* In:
*   NLMIN  : minimum level 
*   NLMAX  : maximum level 
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for all levels of the source grids
*            that should be adapted
*   SOLTRI : array [1..SZTRIA,1..NLEV] of integer
*            If IMSCRT <> 0:
*              Triangulation structures for all levels of the grids that
*              correspond to the solution/RHS vectors in DUP/DRHS.
*              Allowes to calculate error measures.
*            If IMSCRT=0:
*              Dummy, points to TRIAS.
*   IMSCRT : Mesh adaption criterion.
*            =0: grid adaption by geometric aspects
*            =1: grid adaption using geometric aspects and H1-error;
*                Error measure can be calculated using DUP,DRHS,DAUX.
*                SOLTRI corresponds to DUP,DRHS.
*
*   IPARAM : array [1..SZISDI] of integer
*   DPARAM : array [1..SZISDI] of double
*            Integer and double prec. parameter blocks that define the
*            behaviour of the nonstationary solver. 
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
*   VECDAT : array [1..SZN2VI] of integer 
*            TNS2DVectorParams-structure, corresponding to TRIA.
*            Defines the form of the RHS vector.
*   NUVP   : Total length of solution vector
*   DUP    : array [1..NUVP] of double
*            Current solution vector
*   DRHS   : array [1..NUVP] of double
*            Right hand side for the next time step
*   DAUX   : array [1..NUVP] of double
*            Auxiliary vector.
*
*   TIMEST - Initial simulation time
*   TIMENS - Current simulation time
*   TSTEP  - Current time step size. Can vary from time-step to
*            time-step because of the Fractional-Step scheme
*            or the adaptive time stepping.
* Out:
*   TRIAS  - Adapted triangulation on all levels.
************************************************************************

       SUBROUTINE DFGGRI (NLMIN,NLMAX,TRIAS,SOLTRI,IMSCRT,
     *                    IPARAM,DPARAM,ISTPAR,DSTPAR,
     *                    IMGPAR,DMGPAR,IASMBL,DASMBL,
     *                    IGEOM,DGEOM,
     *                    VECDAT,NUVP,DUP,DRHS,DAUX,
     *                    TIMEST, TIMENS, TSTEP)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'sgridadaptmgstaticparams.inc'
      INCLUDE 'cgridadapt.inc'

      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
C     parameters
      
      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,NNLEV),SOLTRI(SZTRIA,NNLEV)
      INTEGER IPARAM(*),ISTPAR(*)
      INTEGER IMGPAR(*),IASMBL(*)
      INTEGER IGEOM(*),VECDAT(SZN2VI)
      INTEGER IMSCRT,NUVP
      
      DOUBLE PRECISION DPARAM(*),DSTPAR(*),DMGPAR(*),DASMBL(*)
      DOUBLE PRECISION DUP(*),DRHS(*),DAUX(*),DGEOM(*)
      DOUBLE PRECISION TIMEST, TIMENS, TSTEP
      
C     local variables

      INTEGER IDATA (3+SZTRIA*NNLEV),LU,LV,LP,KCORVG,KXNPR,I,NVT
      INTEGER LERR,KERR,KVOL
      INTEGER LCOUNT,LTMP,KTMP
      DOUBLE PRECISION DM
      INTEGER NEQU,NEQP      
      
      EXTERNAL UE,PE
      DOUBLE PRECISION UE,PE
      
C     Get the vector size of velocity and pressure part:

      NEQU = VECDAT(ONU)
      NEQP = VECDAT(ONP)
      
C     Adapt the grid if activated in the user

      IF (IGAMET.EQ.3) THEN
      
C       Should we use an error function as a monitor function?

        IF (IMSCRT.NE.0) THEN
        
C         When using an error distribution as monitor function,
C         we use IDATA and DDATA to pass necessary information
C         to the monitor-function callback routine FMON. IDATA has the
C         following internal structure (used in FMON):
C
C         IDATA(1) - =0: only use geometry details for monitor function
C                    =1: only use error distribution as monitor function
C                    =2: use both, geometric details and error distribution
C         IDATA(2) - Minimum level in SOLTRI
C         IDATA(3) - Maximum level in SOLTRI
C         IDATA(4..*) - This saves a copy of the mesh information in SOLTRI,
C                       which defines the underlying mesh for the error
C                       distribution
C         DDATA(1..NVT(NLMAX)) - Saves the error distribuition in all
C                                vertices of the triangulation SOLTRI.
C        
C         Switch on the handling of the error functional in the
C         monitor function by setting IDATA(1) to 2.

          IDATA(1) = IMSCRT
          
C         Incorporate the level/mesh information into IDATA

          IDATA(2) = NLMIN
          IDATA(3) = NLMAX
          CALL LCP3(SOLTRI,IDATA(4),SZTRIA*NNLEV)
          
C         Calculate the error distribution in all vertices.
C         This is done by interpolation of the solution vector
C         to Q1, followed by a H1-error calculation.
C
C         So interpolate the solution to Q1...

          LU = 0
          LV = 0
          LP = 0
          CALL XINTUV (DUP(1),DUP(1+NEQU),SOLTRI(1,NLMAX),LU,LV)
          CALL XINTPV (DUP(1+2*NEQU),TRIAS(1,NLMAX),LP)
          
C         Implement boundary conditions, since the interpolation
C         does not correctly handle boundaries:

          KCORVG = L(SOLTRI(OLCORVG,NLMAX))
          KXNPR  = L(SOLTRI(OLXNPR,NLMAX))
          CALL BDRCOR (DWORK(L(LU)),DWORK(L(LV)),
     *                 SOLTRI(1,NLMAX),DWORK(KCORVG),
     *                 KWORK(KXNPR),UE,1D0,DASMBL(ORE),
     *                 IASMBL,DASMBL,IGEOM,DGEOM)
          CALL BDRCRP (DWORK(L(LP)),
     *                 SOLTRI(1,NLMAX),DWORK(KCORVG),
     *                 KWORK(KXNPR),PE,1D0,DASMBL(ORE),
     *                 IASMBL,DASMBL,IGEOM,DGEOM)
          
C         Calculate the H1-error to LERR

          LERR = 0
          CALL ZNEW(SOLTRI(ONVT,NLMAX)*2,1,LERR,'DERR  ')

          CALL ERPQH1(DWORK(L(LU)),DWORK(L(LV)),SOLTRI(ONVT,NLMAX),
     *                SOLTRI(1,NLMAX),DWORK(L(LERR)))
          CALL ERPCH1(DWORK(L(LP)),SOLTRI(ONVT,NLMAX),
     *                SOLTRI(1,NLMAX),DWORK(L(LERR)+SOLTRI(ONVT,NLMAX)))
     
C         Release memory as far as possible

          CALL ZDISP(0,LP,'DP    ')
          CALL ZDISP(0,LV,'DV    ')
          CALL ZDISP(0,LU,'DU    ')
          
C         The error is currently much too big. We have to rescale
C         it to the range [0..1] for that thing to be a monitor
C         function!

          KERR = L(LERR)

C         At first make it smaller using the SQRT function. Find the
C         maximum - the minimum is 0.0 anyway, as DERR is computed
C         from a norm...
C         Add the error from velocity and pressure to a common value.

          NVT = SOLTRI(ONVT,NLMAX)
          DM = 0D0
          DO I=0,NVT-1
            DWORK(KERR+I) = SQRT(DWORK(KERR+I)) + 
     *                      SQRT(DWORK(KERR+I+NVT))
            DM = MAX(DM,DWORK(KERR+I))
          END DO
          
C         Then rescale to [0..1]. Subtract the value from 1.0, as
C         the DERR-value is large where the cell should get smaller - 
C         and the monitor function for the grid adaption wants
C         to have it vice versa...

          IF (DM.EQ.0) DM=1D0
          DM = 1D0/DM
          DO I=0,NVT-1
            DWORK(KERR+I) = 1D0 - DWORK(KERR+I)*DM
          END DO
          
C         Nice - but this is not the monitor function yet!
C         It's an indicator relative to the current grid!
C         The grid deformation on the other hand wants to have the
C         absolute cell size. Therefore we have to weight the DERR
C         function with the current cell size to get a new absolute
C         cell size distribution.
C
C         At first use the function MEASR2 (which can be found in the
C         grid adaption submodule) to calculate information about the
C         size of the cells.

          CALL ZNEW(SOLTRI(ONVT,NLMAX),3,LCOUNT,'KCOUNT ')
          CALL ZNEW(2*SOLTRI(ONEL,NLMAX)+SOLTRI(ONVT,NLMAX),
     *              1,LTMP,'DTMP  ')
     
C         Write the average cell size around each node to DWORK(L(LTMP)).
C         Collect information we don't use her after this array.

          CALL MEASR2(DWORK(L(SOLTRI(OLCORVG,NLMAX))),
     *                KWORK(L(SOLTRI(OLVERT ,NLMAX))),
     *                KWORK(L(SOLTRI(OLADJ  ,NLMAX))),
     *                SOLTRI(ONEL,NLMAX),SOLTRI(ONVT,NLMAX),
     *                KWORK(L(LCOUNT)),
     *                DWORK(L(LTMP)),
     *                DWORK(L(LTMP)+SOLTRI(ONVT,NLMAX)),
     *                DWORK(L(LTMP)+SOLTRI(ONVT,NLMAX)
     *                             +SOLTRI(ONEL,NLMAX)))

C          DO I=0,NVT-1
C            DWORK(KERR+I) = DBLE(AINT(DWORK(KERR+I)*100.0))/100.0
C          END DO
          
C         Multiply DERR by the average cell size around each node.
C         Note that MEASR2 calculated the normalized cell size,
C         which is restricted to [0,1]. Therefore, the multiplication
C         stays in [0,1]!

          KTMP = L(LTMP)
          DO I=0,NVT-1
            DWORK(KERR+I) = (DWORK(KERR+I)*DWORK(KTMP+I))
            DWORK(KERR+I) = (DBLE(AINT((DWORK(KERR+I)+1D-5)*1000D0))
     *                      /1000D0)
          END DO

C         Now we can use DERR as a monitor function...
C  
C         Release the arrays with the cell size - we don't need
C         them any more.

          CALL ZDISP(0,LTMP,'DTMP  ')
          CALL ZDISP(0,LCOUNT,'KCOUNT')
          
          CALL XSMMGW (TRIAS,NLMIN,NLMAX,IGEOM,DGEOM,
     *                 IDATA,DWORK(L(LERR)))
     
C         Release the rest of the memory

          CALL ZDISP(0,LERR,'DERR  ')
  
        ELSE

C         Switch off the handling of the error functional in the
C         monitor function by setting IDATA(1) to 0.
        
          IDATA(1) = 0

C         Call the grid adaption

          CALL XSMMGW (TRIAS,NLMIN,NLMAX,IGEOM,DGEOM,IDATA,0D0)
          
        END IF
      
      END IF
      
      END
      