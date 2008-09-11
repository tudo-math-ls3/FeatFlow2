************************************************************************
* This file contains the standard pre-/postprocessing routines for the
* stationary solver.
* Postprocessing routines are not called by the solver directly.
* After solving a problem, the user has to call these routines
* to perform the postprocessing.
************************************************************************

************************************************************************
* Default Material Classification for cells, for GMV output
*
* This is a callback routine for the GMV library. It can be used
* in the postprocessing to assign default material constants to
* cells.
* The handling here checks the midpoint of a cell. Depending on whether
* the midpoint is in the interior or in a fictitious boundary component,
* an appropriate material number is assigned to the cell.
*
* In:
*   TRIA   : triangulation structure
*   IENTRY : Number of the cell in the triangulation
*   IPARAM,
*   DPARAM : not used
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to fictitious boundary
*            routines. Not used in this routine.
* 
* Out:
*   MAT    : Material;
*            1 for interior cells,
*            1+#fict. bdry comp. if the cell is in a fict. bdry. comp.
************************************************************************

      SUBROUTINE DFGMMC (TRIA,IENTRY,IPARAM,DPARAM,MAT)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INTEGER TRIA(SZTRIA)
      INTEGER IENTRY,IPARAM(*),MAT
      DOUBLE PRECISION DPARAM(*)
      
C     local variables

      INTEGER KVERT,KCORVG,I,J
      DOUBLE PRECISION X,Y
      
C     externals

      INTEGER ISFBDY
      EXTERNAL ISFBDY
      
C     Fetch some variables for easier access

      KVERT  = L(TRIA(OLVERT))
      KCORVG = L(TRIA(OLCORVG))
      
C     Calculate the midpoint of the element

      X=0
      Y=0
      DO I=0,TRIA(ONVE)-1
        J = KWORK(KVERT+(IENTRY-1)*NNVE+I)-1
        X = X + DWORK(KCORVG+2*J)
        Y = Y + DWORK(KCORVG+2*J+1)
      END DO
      X=X/DBLE(TRIA(ONVE))
      Y=Y/DBLE(TRIA(ONVE))
      
C     Check what this vertex is for a vertex

      MAT = 1+ABS(ISFBDY (X,Y,0,IPARAM,DPARAM))
      
C     If we have a material, subtract the number of boundary
C     components, since fictitious boundary components
C     have numbers > NBCT!

      IF (MAT.NE.1) MAT = MAT - TRIA(ONBCT)
      
      END
      
************************************************************************
* Default Material Classification for vertices, for GMV output
*
* This is a callback routine for the GMV library. It can be used
* in the postprocessing to assign default material constants to
* cells.
* The handling here checks the coordinates of a vertex. The vertex
* can be a corner vertex, a vertex on an edge (defined by DCORMG) or
* an element midpoint.
* Depending on whether this is in the interior or in a fictitious 
* boundary component, an appropriate material number is assigned.
*
* In:
*   TRIA   : triangulation structure
*   IENTRY : Number of the vertex in the triangulation
*   IPARAM,
*   DPARAM : not used
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to fictitious boundary
*            routines. Not used in this routine.
* 
* Out:
*   MAT    : Material;
*            1 for interior vertices,
*            1+#fict. bdry comp. if the vertex is in a fict. bdry. comp.
************************************************************************

      SUBROUTINE DFGMMV (TRIA,IENTRY,IPARAM,DPARAM,MAT)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INTEGER TRIA(SZTRIA)
      INTEGER IENTRY,IPARAM(*),MAT
      DOUBLE PRECISION DPARAM(*)
      
C     local variables

      INTEGER KVERT,KCORVG
      DOUBLE PRECISION X,Y
      
C     externals

      INTEGER ISFBDY
      EXTERNAL ISFBDY
      
C     Fetch some variables for easier access

      KVERT  = L(TRIA(OLVERT))
      KCORVG = L(TRIA(OLCORVG))
      
C     Calculate the coordinates of the vertex

      CALL NDE2XY (IENTRY,TRIA,X,Y)

C     Check what this vertex is for a vertex

      MAT = 1+ABS(ISFBDY (X,Y,0,IPARAM,DPARAM))
      
C     If we have a material, subtract the number of boundary
C     components, since fictitious boundary components
C     have numbers > NBCT!

      IF (MAT.NE.1) MAT = MAT - TRIA(ONBCT)
      
      END

************************************************************************
* Default GMV output
*
* This routine is a standard routine for writing GMV files.
* It transforms a velocity/pressure solution vector from the Q1~/Q0
* space into Q1/Q0 space, writes it out and removes the temporary
* data again. It can be used outside from DPFSTA for "simply writing
* a GMV file". The routine should not be used in more advanced
* postprocessing routines, as it does not allow to return the
* interpolated solution vectors, which could be used for something
* else.
* The GMV-file must already be open and the header must be written
* to the file. When the routine returns, the file is still open and
* the footer is not written. This allowes the caller to add
* further information to the GMV file if necessary.
*
* In:
*   MGMV   - Handle of an output channel where the GMV output should
*            be written to. The channel must be open and will not be
*            closed by this routine!
*   TRIAS  - array [1..SZTRIA,1..*] of integer
*            Triangulation structures of all levels
*   NLMAX  - Maximum level in TRIAS
*   LVL    - Level, where the solution should be written out
*   NEQU   - Length of each velocity-component in DUP
*   DUP    - array [1..*] of double precision
*            Solution vector on level NLMAX.
*            Must have the form (U,V,P,...).
*   RE     - Reynolds number
*   IASMBL,
*   DASMBL - Assembly data structure.
*            Is passed to subroutines
*   IGEOM,
*   DGEOM  - Geometry data structure.
*            Is passed to subroutines.
************************************************************************
      
      SUBROUTINE AUXGMV (MGMV,TRIAS,NLMAX,LVL,NEQU,DUP,RE,
     *                   IASMBL,DASMBL,IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'
      
C     parameters

      INTEGER MGMV,TRIAS(SZTRIA,*),IGEOM(*),IASMBL(*)
      INTEGER NEQU,NEQP,LVL,NLMAX
      DOUBLE PRECISION DUP(*),DGEOM(*),DASMBL(*),RE
      
C     local variables

      INTEGER NCELLS,NVERTS,LU,LV,LP,LISO,K,I,KCORVG,KXNPR,LERR
      
      EXTERNAL DFGMMC,DFGMMV,UE
      
      K = MGMV
      I = LVL
      
      LU = 0
      LV = 0
      LP = 0
      LISO = 0
      LERR = 0
      
C     Write triangulation on level I.
C     Obtain NCELLS and NVERTS.

      CALL GMVTRI (K,TRIAS(1,I),1,NCELLS,NVERTS)
      
C     Write materials;
C     we don't give material names for now.

      CALL GMVMAT (K,TRIAS(1,I),0,NCELLS,TRIAS(ONEL,I),
     *             -16,'',DFGMMC,IGEOM,DGEOM)
      CALL GMVMAT (K,TRIAS(1,I),1,NVERTS,
     *             TRIAS(ONVT,I)+TRIAS(ONMT,I),
     *             -16,'',DFGMMV,IGEOM,DGEOM)
      
C     Calculate the velocity field, pressure and the stream 
C     function if they are not calculated already:

      IF ((LU.EQ.0).OR.(LV.EQ.0)) THEN
      
        CALL XINTUV (DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),LU,LV)
      
C       Implement boundary conditions, since the interpolation
C       does not correctly handle boundaries:

        KCORVG = L(TRIAS(OLCORVG,NLMAX))
        KXNPR  = L(TRIAS(OLXNPR,NLMAX))
        CALL BDRCOR (DWORK(L(LU)),DWORK(L(LV)),
     *               TRIAS(1,NLMAX),DWORK(KCORVG),
     *               KWORK(KXNPR),UE,1D0,RE,
     *               IASMBL,DASMBL,IGEOM,DGEOM)
      END IF

      IF (LP.EQ.0) THEN
        CALL XINTPV (DUP(1+2*NEQU),TRIAS(1,NLMAX),LP)
      END IF

      IF (LISO.EQ.0) THEN
        CALL XU2ISO (DUP(1),DUP(1+1*NEQU),TRIAS(1,NLMAX),LISO)
      END IF

C     Calculate the H1-error to LERR

      IF (LERR.EQ.0) THEN
        CALL ZNEW(TRIAS(ONVT,NLMAX),1,LERR,'DERR  ')
      END IF

      CALL ERPQH1(DWORK(L(LU)),DWORK(L(LV)),TRIAS(ONVT,NLMAX),
     *            TRIAS(1,NLMAX),DWORK(L(LERR)))
      CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *             DWORK(L(LERR)),'H1error')
      CALL ERPCH1(DWORK(L(LP)),TRIAS(ONVT,I),TRIAS(1,NLMAX),
     *            DWORK(L(LERR)))
      CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *             DWORK(L(LERR)),'H1errorP')
        
C     Write all these at the desired level to the GMV file:

      CALL GMVVEL (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *             DWORK(L(LU)),DWORK(L(LV)))
      CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *             DWORK(L(LP)),'pressure')
      CALL GMVSCA (K,TRIAS(1,I),0,NCELLS,TRIAS(ONEL,I),
     *             DUP(1+2*NEQU),'pressure')
      CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *             DWORK(L(LISO)),'streamfunction')
      
C     Write geometry:

      CALL FBDGMV (K,2,IGEOM,DGEOM)
      
C     Release memory, finish.

      IF (LERR.NE.0) CALL ZDISP(0,LERR,'DERR  ')
      IF (LISO.NE.0) CALL ZDISP(0,LISO,'DISO  ')
      IF (LU.NE.0) CALL ZDISP(0,LU,'DU    ')
      IF (LV.NE.0) CALL ZDISP(0,LV,'DV    ')
      IF (LP.NE.0) CALL ZDISP(0,LP,'DP    ')

      END
      
      
************************************************************************
* Default Pre- and postprocessing for stationary solver
*
* This routine can be called with a solution vector from the stationary
* solver to perform standard postprocessing like writing GMV-files,
* writing plain solution vectors,...
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
*   ITRACK : array {1..SZTRKI] of integer
*            TTrackingIParams structure that holds the output handles
*            for writing tracked solution values to files.
*
*   NUVP   : Total length of solution vector; usually 2*NEQU+NEQP
*   DUP    : array [1..NUVP] of double
*            Current solution vector
*   DRHS   : array [1..NUVP] of double
*            Right hand side for the next time step
*   DAUX   : array [1..NUVP] of double
*            Auxiliary vector.
*
*   IERANA : Cubature rule to use for error analysis in comparison
*            to reference solution.
*            =0: Don't use error analysis
*
*   IGMV   : Save GMV file to disc
*            =0: Don't save GMV-file to disc.
*            >0: Save GMV-file to disc on level IGMV.
*
*   IBDFBD : Cubature formula to use for calculation of body forces
*            (drag/lift) when calculating the forces by line integration
*            on the real boundary.
*            =0: Don't calculate body forces on real boundary
*            >0: Use line cubature formula IBDFBD
*
*   IBDFVI : Cubature formula to use for calculation of body forces
*            (drag/lift) when calculating the forces by volume
*            integration. This is used for the calculation of forces
*            arising on fictitious boundary components.
*            =0: Don't calculate body forces on fictitious boundary obj.
*            >0: Element cubature formula for CB2Q
*            =4: 2x2 Gauss; standard
*
*   IBDFLI : Cubature formula to use for calculation of body forces
*            (drag/lift) when calculating the forces by line
*            integration on the reconstructed interface of fictitious
*            boundary components.
*            =0: Don't calculate body forces on fictitious boundary
*                objects by line integration on the reconstructed
*                interface
*            >0: Element cubature formula for CB1
*            =2: Trapezoidal rule; standard
*
*   ISOLLV : Whether or not to write out the plain solution vector
*            to disc.
*             =0: Don't write solution vector to disc
*            <>0: Save solution vector to disc on level |ISOLLV|.
*                 >0: formatted output (standard)
*                 <0: unformatted output.
*   IFUSAV : Save vertex-based velocity solutions to disc.
*             =0: Don't write vertex-based velocity solutions to disc
*            <>0: Save intermediate velocities to disc on level |IFUSAV|.
*                 X-velocity is saved to "./film/DU"
*                 X-velocity is saved to "./film/DV"
*                 >0: formatted output
*                 <0: unformatted output
*   IFPSAV : Save vertex-based pressure solutions to disc.
*             =0: Don't write vertex-based pressure solutions
*                 to disc
*            <>0: Save pressures to disc on level IFPSAV.
*                 The pressure is saved to "./film/DP".
*                 >0: formatted output (standard)
*                 <0: unformatted output
*   IFXSAV : Save vertex-based streamline functions to disc.
*             =0: Don't write streamfunction to disc
*            <>0: Save streamfunction to disc on level IFXSAV.
*                 The pressure is saved to "./film/DISO".
*                 >0: formatted output (standard)
*                 <0: unformatted output
*   LFNAMS : array [1..6] of integer
*            Handles to DStrings that define the filenames for 
*            file output. If a handle is 0, a standard filename will be
*            chosen.
*            LFNAMS[1]  = Filename for output of solution vector
*            LFNAMS[2]  = Filename for film output, X-velocity output
*            LFNAMS[3]  = Filename for film output, Y-velocity output
*            LFNAMS[4]  = Filename for film output, pressure output
*            LFNAMS[5]  = Filename for film output, streamline 
*            LFNAMS[6]  = Filename of GMV-files

*  VECDAT in that case, so the content of these vectors are treated
*  as undefined on entry and return of this routine.
*
*  The routine stops the time needed for postprocessing and adds this
*  to the timing substructure in DSTPAR!
************************************************************************
      
      SUBROUTINE DFPSTA (NLMIN,NLMAX,
     *                   TRIAS,MATDAT,VECDAT,
     *                   ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *                   ITRACK,NUVP,DUP,DRHS,DAUX,
     *                   IGMV,IERANA,IBDFBD,IBDFVI,IBDFLI,
     *                   ISOLLV,IFUSAV,IFPSAV,IFXSAV,LFNAMS)
      
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
      INTEGER IMGPAR(*),ISTPAR(*),IASMBL(*),IGEOM(*),ITRACK(*),LFNAMS(6)
      INTEGER TRIAS(SZTRIA,NNLEV)
      INTEGER MATDAT(SZN2MI,NNLEV),VECDAT(SZN2VI,NNLEV),NLMIN,NLMAX
      DOUBLE PRECISION DMGPAR(*),DSTPAR(*),DASMBL(*)
      DOUBLE PRECISION DGEOM(*)
      
      DOUBLE PRECISION DUP(*),DRHS(*),DAUX(*)
      
      INTEGER IGMV,IERANA,IBDFBD,IBDFVI,IBDFLI,IFUSAV,IFPSAV,IFXSAV
      INTEGER ISOLLV
      
C     Coefficient of exact solution

      DOUBLE PRECISION UE,PE,UEX,UEY
      EXTERNAL UE,PE,UEX,UEY
      
C     finite elements

      EXTERNAL E030,E031,EM30,EM31,E010,E011
      
C     externals for GMV output

      EXTERNAL DFGMMC,DFGMMV
      
C     maximal number of tracked solutions
      
      INTEGER MAXTRK
      PARAMETER (MAXTRK=32)
      
C     local variables

      INTEGER I,J,K,MSHOW,IEL,NEQU,NEQP
      INTEGER LU,LV,LP,LISO,LERR
      INTEGER KCORVG,KXNPR
      INTEGER KVERT,KMID,KVBD,KEBD,KVBDP,KBCT,KCORMG,NCELLS,NVERTS
      INTEGER IINFO(MAXTRK)
      DOUBLE PRECISION DINFO1(MAXTRK),DINFO2(MAXTRK),DINFO3(MAXTRK)
      DOUBLE PRECISION D1,D2,COEF1,COEF2,DTMP1,DTMP2
      CHARACTER CNAME*(60)
      CHARACTER CFN*60
      
C     A backup of the solver structure for time measurements

      DOUBLE PRECISION TIMNG (SZNSDD)

C     Ok, this routine seems to get a little bit stretchy now, but this
C     lies in the nature of postprocessing, i.e. handling of all the 
C     cases :)

C     At first stop the time and iitialize some local variables

      CALL GTMAUX (TIMNG,DSTPAR,OTNLTIM-1+OTTPOST,0)

      MSHOW = ISTPAR(OMSGTRM)
      
C     Set LU,LV,LP,LISO to 0. These are handles to vertex-based solution
C     vectors which might be calculated during the postprocessing if
C     needed; in this case we have to release them later.

      LU = 0
      LV = 0
      LP = 0
      LISO = 0
      LERR = 0

C     Get the vector size of velocity and pressure part:

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

      IF (IBDFBD.GT.0) THEN
      
C       In how many boundary components should the body force be
C       evaluated? We store that number in J.

        CALL FPTSIN(1,0,1D0,DASMBL(ORE),
     *              IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *              J,D1,D2)
     
C       We only track a limited number of points

        J = MIN(J,MAXTRK)
     
C       Get the coefficients of the boundary integral into COEF1
C       and COEF2:

        CALL FPTSIN(5,0,1D0,DASMBL(ORE),
     *              IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *              0,COEF1,COEF2)
     
C       Get the boundary component segments where to evaluate the
C       forces. For now we restrict to a maximum of MAXTRK segments.
C       Evaluate there!

        DO I=1,J
        
C         Request data about the boundary component.
C         Its number is stored in K and the parameter values in D1/D2.

          CALL FPTSIN(2,I,1D0,DASMBL(ORE),
     *                IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                K,D1,D2)
     
          IINFO(I) = K
     
          IF (IASMBL(OIELEMT).EQ.0) THEN
          
            CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *              TRIAS(1,NLMAX),KWORK(KVERT),
     *              KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *              DWORK(KCORVG),DWORK(KVBDP),E031,.FALSE.,
     *              K,D1,D2,COEF1,COEF2,
     *              1,DINFO1(I),DINFO2(I),DINFO3(I))   
     
          ELSE IF (IASMBL(OIELEMT).EQ.1) THEN
          
            CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *              TRIAS(1,NLMAX),KWORK(KVERT),
     *              KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *              DWORK(KCORVG),DWORK(KVBDP),E030,.FALSE.,
     *              K,D1,D2,COEF1,COEF2,
     *              1,DINFO1(I),DINFO2(I),DINFO3(I))   
     
          ELSE IF (IASMBL(OIELEMT).EQ.2) THEN
          
            CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *              TRIAS(1,NLMAX),KWORK(KVERT),
     *              KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *              DWORK(KCORVG),DWORK(KVBDP),EM31,.TRUE.,
     *              K,D1,D2,COEF1,COEF2,
     *              1,DINFO1(I),DINFO2(I),DINFO3(I))   
     
          ELSE IF (IASMBL(OIELEMT).EQ.3) THEN
          
            CALL BDFORX(DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *              TRIAS(1,NLMAX),KWORK(KVERT),
     *              KWORK(KMID),KWORK(KVBD),KWORK(KEBD),KWORK(KBCT),
     *              DWORK(KCORVG),DWORK(KVBDP),EM30,.TRUE.,
     *              K,D1,D2,COEF1,COEF2,
     *              1,DINFO1(I),DINFO2(I),DINFO3(I))  
     
          END IF 
        
        END DO
        
C       Print the results to screen/file:

        CALL TV2TER (MSHOW,MFILE,J,IINFO,DINFO1,DINFO2,2,
     *    'Body forces real bd., constant pressure, bdc/horiz/vert')
      
C       Write to tracking file, if handles are given.
C       Write horizontal and vertical force.

        IF (ITRACK(OHBFB1).GT.0) THEN
          CALL AUXCFN ('BDC','',CNAME,-1)
          CALL TV2FIL (ITRACK(OHBFB1),J,IINFO,DINFO1,
     *                 .TRUE.,1D0,-1,CNAME)
        END IF

        IF (ITRACK(OHBFB2).GT.0) THEN
          CALL AUXCFN ('BDC','',CNAME,-1)
          CALL TV2FIL (ITRACK(OHBFB2),J,IINFO,DINFO2,
     *                 .TRUE.,1D0,-1,CNAME)
        END IF
      
      END IF ! IBDFBD <> 0

C     --------------------------------------------------------------      
C     Boundary forces evaluation on fictitious boundary 
C     by volume integration.
C     Is there a cubature formula for the "real boundary evaluation"
C     given?

      IF (IBDFVI.GT.0) THEN
      
C       In how many boundary components should the body force be
C       evaluated? We store that number in J.

        CALL FPTSIN(3,0,1D0,DASMBL(ORE),
     *              IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *              J,D1,D2)
     
C       We only track a limited number of points

        J = MIN(J,MAXTRK)
        
        IF (J.NE.0) THEN
        
C         Get the coefficients of the boundary integral into COEF1
C         and COEF2:

          CALL FPTSIN(5,0,1D0,DASMBL(ORE),
     *                IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                0,COEF1,COEF2)
     
C         Get the boundary component segments where to evaluate the
C         forces. For now we restrict to a maximum of MAXTRK segments.
C         Evaluate there!
C         If J is -1, we have to evaluate only once - for all fictitious
C         boundary components at the same time!

          DO I=1,MAX(1,J)
          
C           Request data about the boundary component.
C           The number of the fict. bdry component is stored to K.
C           If J=-1, we set K=0 to evaluate everywhere.

            IF (J.EQ.-1) THEN
              K = 0
            ELSE
              CALL FPTSIN(4,I,1D0,DASMBL(ORE),
     *                  IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                  K,D1,D2)
            END IF
            
C           Save the boundary component number to IINFO:              
            
            IINFO(I) = K
     
C           Call the evaluation routine, calculate forces on that
C           boundary component
     
            IF (IASMBL(OIELEMT).EQ.0) THEN
            
              CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *                 DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                 KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *                 E031,.FALSE.,DINFO1(I),DINFO2(I),COEF1,COEF2,
     *                 K,0.5D0,IGEOM,DGEOM)
            
            ELSE IF (IASMBL(OIELEMT).EQ.1) THEN
            
              CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *                 DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                 KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *                 E030,.FALSE.,DINFO1(I),DINFO2(I),COEF1,COEF2,
     *                 K,0.5D0,IGEOM,DGEOM)

            ELSE IF (IASMBL(OIELEMT).EQ.2) THEN
            
              CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *                 DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                 KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *                 EM31,.TRUE.,DINFO1(I),DINFO2(I),COEF1,COEF2,
     *                 K,0.5D0,IGEOM,DGEOM)

            ELSE IF (IASMBL(OIELEMT).EQ.3) THEN
            
              CALL BDFVLG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),NEQU,
     *                 DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                 KWORK(KMID),KWORK(KXNPR),TRIAS(1,NLMAX),
     *                 EM30,.TRUE.,DINFO1(I),DINFO2(I),COEF1,COEF2,
     *                 K,0.5D0,IGEOM,DGEOM)

            END IF 
          
          END DO
          
C         Print the results to screen/file:

          CALL TV2TER (MSHOW,MFILE,MAX(J,1),IINFO,DINFO1,DINFO2,2,
     *      'Body forces fict.bd., constant pressure, bdc/horiz/vert')
      
C         Write to tracking file, if handles are given.
C         Write horizontal and vertical force.

          IF (ITRACK(OHBFV1).GT.0) THEN
            CALL AUXCFN ('BDC','',CNAME,-1)
            CALL TV2FIL (ITRACK(OHBFV1),MAX(J,1),IINFO,DINFO1,
     *                   .TRUE.,1D0,-1,CNAME)
          END IF

          IF (ITRACK(OHBFV2).GT.0) THEN
            CALL AUXCFN ('BDC','',CNAME,-1)
            CALL TV2FIL (ITRACK(OHBFV2),MAX(J,1),IINFO,DINFO2,
     *                   .TRUE.,1D0,-1,CNAME)
          END IF

        END IF
      
      END IF ! IBDFVL <> 0
      
C     --------------------------------------------------------------      
C     Boundary forces evaluation on fictitious boundary 
C     by line integration.
C     Is there a cubature formula for the "real boundary evaluation"
C     given?

      IF (IBDFLI.GT.0) THEN
      
C       In how many boundary components should the body force be
C       evaluated? We store that number in J.

        CALL FPTSIN(3,0,1D0,DASMBL(ORE),
     *              IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *              J,D1,D2)
     
C       We only track a limited number of points

        J = MIN(J,MAXTRK)
        
        IF (J.NE.0) THEN
        
C         Get the coefficients of the boundary integral into COEF1
C         and COEF2:

          CALL FPTSIN(5,0,1D0,DASMBL(ORE),
     *                IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                0,COEF1,COEF2)
     
C         Get the boundary component segments where to evaluate the
C         forces. For now we restrict to a maximum of MAXTRK segments.
C         Evaluate there!
C         If J is -1, we have to evaluate only once - for all fictitious
C         boundary components at the same time!

          DO I=1,MAX(1,J)
          
C           Request data about the boundary component.
C           The number of the fict. bdry component is stored to K.
C           If J=-1, we set K=0 to evaluate everywhere.

            IF (J.EQ.-1) THEN
              K = 0
            ELSE
              CALL FPTSIN(4,I,1D0,DASMBL(ORE),
     *                  IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                  K,D1,D2)
            END IF
     
C           Save the boundary component number to IINFO:              
            
            IINFO(I) = K
     
C           Call the evaluation routine, calculate forces on that
C           boundary component.
C           Pass the assembly data structures as INFO-block to
C           these routines. That way all called subroutines
C           can access information about the current assembly
C           status.
     
            IF (IASMBL(OIELEMT).EQ.0) THEN
            
              CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                     DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                     KWORK(KMID),TRIAS(1,NLMAX),
     *                     E031,.FALSE.,DINFO1(I),DINFO2(I),COEF1,
     *                     COEF2,K,IGEOM,DGEOM)
            
            ELSE IF (IASMBL(OIELEMT).EQ.1) THEN
            
              CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                     DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                     KWORK(KMID),TRIAS(1,NLMAX),
     *                     E030,.FALSE.,DINFO1(I),DINFO2(I),COEF1,
     *                     COEF2,K,IGEOM,DGEOM)

            ELSE IF (IASMBL(OIELEMT).EQ.2) THEN
            
              CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                     DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                     KWORK(KMID),TRIAS(1,NLMAX),
     *                     EM31,.TRUE.,DINFO1(I),DINFO2(I),COEF1,
     *                     COEF2,K,IGEOM,DGEOM)

            ELSE IF (IASMBL(OIELEMT).EQ.3) THEN
            
              CALL BDFRIG (DUP(1),DUP(1+NEQU),DUP(1+2*NEQU),
     *                     DWORK(KCORVG),KWORK(KVERT),DWORK(KCORMG),
     *                     KWORK(KMID),TRIAS(1,NLMAX),
     *                     EM30,.TRUE.,DINFO1(I),DINFO2(I),COEF1,
     *                     COEF2,K,IGEOM,DGEOM)

            END IF 
          
          END DO
          
C         Print the results to screen/file:

          CALL TV2TER (MSHOW,MFILE,MAX(J,1),IINFO,DINFO1,DINFO2,2,
     *      'Body forces f.b.l.i., constant pressure, bdc/horiz/vert')
      
C         Write to tracking file, if handles are given.
C         Write horizontal and vertical force.

          IF (ITRACK(OHBFL1).GT.0) THEN
            CALL AUXCFN ('BDC','',CNAME,-1)
            CALL TV2FIL (ITRACK(OHBFL1),MAX(J,1),IINFO,DINFO1,
     *                   .TRUE.,1D0,-1,CNAME)
          END IF

          IF (ITRACK(OHBFL2).GT.0) THEN
            CALL AUXCFN ('BDC','',CNAME,-1)
            CALL TV2FIL (ITRACK(OHBFL2),MAX(J,1),IINFO,DINFO2,
     *                   .TRUE.,1D0,-1,CNAME)
          END IF
      
        END IF
      
      END IF ! IBDFVLI <> 0
      
C     --------------------------------------------------------------      
C     Output of tracked solution values.
C
C     At first: Is the output of tracked solution values active?
C     Check if we have a file handle where to output these values
C     to.
C
C     ----------
C     Velocities

      IF ((ITRACK(OHVEL1).GT.0).OR.(ITRACK(OHVEL2).GT.0)) THEN
      
C       Calculate the interpolated velocity, if it's not calculated

        IF ((LU.EQ.0).OR.(LV.EQ.0)) THEN
        
          CALL XINTUV (DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),LU,LV)
        
C         Implement boundary conditions, since the interpolation
C         does not correctly handle boundaries:

          KCORVG = L(TRIAS(OLCORVG,NLMAX))
          KXNPR  = L(TRIAS(OLXNPR,NLMAX))
          CALL BDRCOR (DWORK(L(LU)),DWORK(L(LV)),
     *                 TRIAS(1,NLMAX),DWORK(KCORVG),
     *                 KWORK(KXNPR),UE,1D0,DASMBL(ORE),
     *                 IASMBL,DASMBL,IGEOM,DGEOM)
        END IF

C       Ask the user defined routine in how many points we should
C       track the velocity; write the number into K.

        CALL FPTSIN(6,0,1D0,DASMBL(ORE),
     *            IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *            J,D1,D2)
     
C       We only track a limited number of points

        J = MIN(J,MAXTRK)
        
        IF (J.GT.0) THEN
        
C         Get the point numbers and the values to
C         IINFO/DINFO1/DINFO2:

          DO I=1,J
            CALL FPTSIN(7,I,1D0,DASMBL(ORE),
     *                  IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                  IINFO(I),D1,D2)
            IF (IINFO(I).NE.0) THEN
              DINFO1(I) = DWORK(L(LU)+IINFO(I)-1)
              DINFO2(I) = DWORK(L(LV)+IINFO(I)-1)
            ELSE
C             Oops, thats a point anywhere in the domain...
C             Search the element where it is and evaluate there.
              CALL PSRCH5(TRIAS,NLMIN,NLMAX,D1,D2,0,IEL)
              IF (IEL.EQ.0) THEN
                DINFO1(I) = 0.0
                DINFO2(I) = 0.0
              ELSE
                KCORVG = L(TRIAS(OLCORVG,NLMAX))
                KVERT  = L(TRIAS(OLVERT,NLMAX))
                KMID   = L(TRIAS(OLMID,NLMAX))
                CALL SCEVLQ (TRIAS(ONVT,NLMAX),DWORK(L(LU)),
     *                       E011,.FALSE.,D1,D2,IEL,
     *                       TRIAS(1,NLMAX),DWORK(KCORVG),
     *                       KWORK(KVERT),KWORK(KMID),
     *                       DINFO1(I), DTMP1, DTMP2) 
                CALL SCEVLQ (TRIAS(ONVT,NLMAX),DWORK(L(LV)),
     *                       E011,.FALSE.,D1,D2,IEL,
     *                       TRIAS(1,NLMAX),DWORK(KCORVG),
     *                       KWORK(KVERT),KWORK(KMID),
     *                       DINFO2(I), DTMP1, DTMP2) 
              END IF
            END IF
          END DO
          
C         Print the results to screen

          CALL TV2TER (MSHOW,MFILE,J,IINFO,DINFO1,DINFO2,2,
     *      'Tracked velocity point-values P(VELO)')

C         Write the values to the file:

          IF (ITRACK(OHVEL1).GT.0) THEN
            CALL AUXCFN ('PT','',CNAME,-1)
            CALL TV2FIL (ITRACK(OHVEL1),J,IINFO,DINFO1,
     *                   .TRUE.,1D0,-1,CNAME)
          END IF

          IF (ITRACK(OHVEL2).GT.0) THEN
            CALL AUXCFN ('PT','',CNAME,-1)
            CALL TV2FIL (ITRACK(OHVEL2),J,IINFO,DINFO2,
     *                   .TRUE.,1D0,-1,CNAME)
          END IF
          
        END IF ! J>0
        
      END IF ! HVEL1>0 or HVEL2 > 0

C     ----------
C     Pressure

      IF (ITRACK(OHPRES).GT.0) THEN
      
C       Calculate the interpolated pressure, if it's not calculated

        IF (LP.EQ.0) THEN
          CALL XINTPV (DUP(1+2*NEQU),TRIAS(1,NLMAX),LP)
        END IF

C       Ask the user defined routine in how many points we should
C       track the velocity; write the number into K.

        CALL FPTSIN(8,0,1D0,DASMBL(ORE),
     *            IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *            J,D1,D2)
     
C       We only track a limited number of points

        J = MIN(J,MAXTRK)
        
        IF (J.GT.0) THEN
        
C         Get the point numbers and the values to IINFO/DINFO1/DINFO2:

          DO I=1,J
            CALL FPTSIN(9,I,1D0,DASMBL(ORE),
     *                  IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                  IINFO(I),D1,D2)
            IF (IINFO(I).NE.0) THEN
              DINFO1(I) = DWORK(L(LP)+IINFO(I)-1)
            ELSE
C             Oops, thats a point anywhere in the domain...
C             Search the element where it is and evaluate there.
              CALL PSRCH5(TRIAS,NLMIN,NLMAX,D1,D2,0,IEL)
              IF (IEL.EQ.0) THEN
                DINFO1(I) = 0.0
              ELSE
                KCORVG = L(TRIAS(OLCORVG,NLMAX))
                KVERT  = L(TRIAS(OLVERT,NLMAX))
                KMID   = L(TRIAS(OLMID,NLMAX))
                CALL SCEVLQ (TRIAS(ONVT,NLMAX),DWORK(L(LP)),
     *                       E011,.FALSE.,D1,D2,IEL,
     *                       TRIAS(1,NLMAX),DWORK(KCORVG),
     *                       KWORK(KVERT),KWORK(KMID),
     *                       DINFO1(I), DTMP1, DTMP2) 
              END IF
            END IF
          END DO
          
C         Print the results to screen

          CALL TV2TER (MSHOW,MFILE,J,IINFO,DINFO1,DINFO2,1,
     *      'Tracked pressure point-values P(PRES)')

C         Write the values to the file:

          IF (ITRACK(OHPRES).GT.0) THEN
            CALL AUXCFN ('PT','',CNAME,-1)
            CALL TV2FIL (ITRACK(OHPRES),J,IINFO,DINFO1,
     *                   .TRUE.,1D0,-1,CNAME)
          END IF

        END IF ! J>0
        
      END IF ! HPRES > 0

C     --------------
C     Streamfunction

      IF (ITRACK(OHSTRM).GT.0) THEN
      
C       Calculate the interpolated streamfunction, if 
C       it's not calculated

        IF (LISO.EQ.0) THEN
          CALL XU2ISO (DUP(1),DUP(1+1*NEQU),TRIAS(1,NLMAX),LISO)
        END IF

C       Ask the user defined routine in how many points we should
C       track the velocity; write the number into K.

        CALL FPTSIN(10,0,1D0,DASMBL(ORE),
     *            IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *            J,D1,D2)
     
C       We only track a limited number of points

        J = MIN(J,MAXTRK)
        
        IF (J.GT.0) THEN
        
C         Get the point numbers and the values to IINFO/DINFO1:

          DO I=1,J
            CALL FPTSIN(11,I,1D0,DASMBL(ORE),
     *                  IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                  IINFO(I),D1,D2)
            IF (IINFO(I).NE.0) THEN
              DINFO1(I) = DWORK(L(LISO)+IINFO(I)-1)
            ELSE
C             Oops, thats a point anywhere in the domain...
C             Search the element where it is and evaluate there.
              CALL PSRCH5(TRIAS,NLMIN,NLMAX,D1,D2,0,IEL)
              IF (IEL.EQ.0) THEN
                DINFO1(I) = 0.0
              ELSE
                KCORVG = L(TRIAS(OLCORVG,NLMAX))
                KVERT  = L(TRIAS(OLVERT,NLMAX))
                KMID   = L(TRIAS(OLMID,NLMAX))
                CALL SCEVLQ (TRIAS(ONVT,NLMAX),DWORK(L(LISO)),
     *                       E011,.FALSE.,D1,D2,IEL,
     *                       TRIAS(1,NLMAX),DWORK(KCORVG),
     *                       KWORK(KVERT),KWORK(KMID),
     *                       DINFO1(I), DTMP1, DTMP2) 
              END IF
            END IF
          END DO
          
C         Print the results to screen

          CALL TV2TER (MSHOW,MFILE,J,IINFO,DINFO1,DINFO2,1,
     *      'Tracked streamfunction point-values P(FLUX)')

C         Write the values to the file:

          IF (ITRACK(OHSTRM).GT.0) THEN
            CALL AUXCFN ('PT','',CNAME,-1)
            CALL TV2FIL (ITRACK(OHSTRM),J,IINFO,DINFO1,
     *                   .TRUE.,1D0,-1,CNAME)
          END IF

        END IF ! J>0
        
      END IF ! OHSTRM > 0

C     -----------------------
C     Mean Integreal Pressure

      IF (ITRACK(OHIPRS).GT.0) THEN
      
C       In how many boundary components should the pressure integral
C       be evaluated? We store that number in J.

        CALL FPTSIN(12,0,1D0,DASMBL(ORE),
     *              IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *              J,D1,D2)
     
C       We only track a limited number of points

        J = MIN(J,MAXTRK)
        
        IF (J.GT.0) THEN
        
C         Get the boundary component segments where to evaluate the
C         forces. For now we restrict to a maximum of MAXTRK segments.
C         Evaluate there!
C         If J is -1, we have to evaluate only once - for all fictitious
C         boundary components at the same time!

          DO I=1,J
          
C           Request data about the boundary component.
C           The number of the fict. bdry component is stored to K.
C           If J=-1, we set K=0 to evaluate everywhere.

            CALL FPTSIN(13,I,1D0,DASMBL(ORE),
     *                IASMBL,DASMBL,IGEOM,DGEOM,TRIAS(1,NLMAX),
     *                K,D1,D2)
            
C           Save the boundary component number to IINFO:              
            
            IINFO(I) = K
     
C           Call the evaluation routine, calculate integral on that
C           boundary component into DINFO1. DINFO2 receives the length
C           of the boundary we integrate on.
     
            CALL BDINTX (DUP(1+2*NEQU),TRIAS(1,NLMAX),
     *               KWORK(KVERT),KWORK(KMID),KWORK(KVBD),KWORK(KEBD),
     *               KWORK(KBCT),DWORK(KCORVG), DWORK(KVBDP),
     *               E010,.FALSE.,
     *               K,D1,D2,1,DINFO1(I),DINFO2(I))
     
C           Divide the integral pressure by the length of the boundary
C           component to get the integral mean pressure

            IF (DINFO2(I).NE.0D0) THEN
              DINFO1(I) = DINFO1(I)/DINFO2(I)
            END IF
            
          END DO
          
C         Print the results to screen

          CALL TV2TER (MSHOW,MFILE,J,IINFO,DINFO1,DINFO2,1,
     *      'Tracked integral mean pressure values I(PRES)')

C         Write to tracking file, if handles are given.
C         Write horizontal and vertical force.

          IF (ITRACK(OHIPRS).GT.0) THEN
            CALL AUXCFN ('BDC','',CNAME,-1)
            CALL TV2FIL (ITRACK(OHIPRS),J,IINFO,DINFO1,
     *                   .TRUE.,1D0,-1,CNAME)
          END IF

        END IF
      
      END IF ! HPRES > 0

C     =================================================================
C     Pre-/postprocessing regarding writing solutions to disc.
C
C     We go on with handling the writing out of solutions.
C
C     -------------------------------------------------------------
C     Should we write out the plain solution vector as it is
C     on a defined level?

      IF (ISOLLV.NE.0) THEN
      
C       We want to write it to disc on level ISOLLV, i.e. at
C       level:

        I = MAX(NLMIN,MIN(NLMAX,ABS(ISOLLV)))
        
C       But before we can write, we must calculate it.
C       We can use e.g. one of the allocated memory vectors in VECDAT
C       on level I to store the vector, as these are only used when
C       solving linear equations - what we clearly not do here!
C       So we can e.g. store it to DWORK(J):

        J = L(VECDAT(OLTMP,I))
        
C       We have to copy two velocity and one pressure vector to that.
C       Because of the 2-level ordering, the first couple of entries
C       on the finest level belong to the coarser levels! So we
C       can concat together:

        K = VECDAT(ONU,I)
        CALL LCP1(DUP(1),DWORK(J),K)
        CALL LCP1(DUP(1+K),DWORK(J+NEQU),K)
        CALL LCP1(DUP(1+2*K),DWORK(J+2*NEQU),VECDAT(ONP,I))
        
C       Write it.
C       The sign of ISOLLV decides on whether the output is in
C       formatted or unformatted form.
        
        CALL STPUT (LFNAMS(1), CFN)
        CALL PPWRVC (0,VECDAT(ONEQV,I),DWORK(J),-1,
     *               ISOLLV,CFN)
     
      END IF

C     -------------------------------------------------------------
C     Should we write out interpolated velocities?

      IF (IFUSAV.NE.0) THEN
        
C       We want to write it to disc on level IFUSAV, i.e. at
C       level:

        I = MAX(NLMIN,MIN(NLMAX,ABS(IFUSAV)))

C       Calculate the interpolated velocity, if it's not calculated

        IF ((LU.EQ.0).OR.(LV.EQ.0)) THEN
        
          CALL XINTUV (DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),LU,LV)
        
C         Implement boundary conditions, since the interpolation
C         does not correctly handle boundaries:

          KCORVG = L(TRIAS(OLCORVG,NLMAX))
          KXNPR  = L(TRIAS(OLXNPR,NLMAX))
          CALL BDRCOR (DWORK(L(LU)),DWORK(L(LV)),
     *                 TRIAS(1,NLMAX),DWORK(KCORVG),
     *                 KWORK(KXNPR),UE,1D0,DASMBL(ORE),
     *                 IASMBL,DASMBL,IGEOM,DGEOM)
        END IF

C       Write out the values in the vertices for the X-velocity
C       The sign of IFUSAV decides on whether the output is in
C       formatted or unformatted form.

        CALL STPUT (LFNAMS(2), CFN)
        CALL PPWRVC (0,TRIAS(ONVT,I),DWORK(L(LU)),-1,
     *               IFUSAV,CFN)
                  
C       and the Y-velocity
        
        CALL STPUT (LFNAMS(3), CFN)
        CALL PPWRVC (0,TRIAS(ONVT,I),DWORK(L(LV)),-1,
     *               IFUSAV,CFN)

      END IF

C     -------------------------------------------------------------
C     Should we write out interpolated pressure values?

      IF (IFPSAV.NE.0) THEN
        
C       We want to write it to disc on level IFPSAV, i.e. at
C       level:

        I = MAX(NLMIN,MIN(NLMAX,ABS(IFPSAV)))

C       Calculate the interpolated pressure, if it's not calculated

        IF (LP.EQ.0) THEN
          CALL XINTPV (DUP(1+1*NEQU),TRIAS(1,NLMAX),LP)
        END IF

C       Write out the values in the vertices for the pressure.
C       Make use of the 2-level ordering to write out only level I.
C       The sign of IFPSAV decides on whether the output is in
C       formatted or unformatted form.

        CALL STPUT (LFNAMS(4), CFN)
        CALL PPWRVC (0,TRIAS(ONVT,I),DWORK(L(LP)),-1,
     *               IFPSAV,CFN)

      END IF
     
C     -------------------------------------------------------------
C     Should we write out the streamfunction?

      IF (IFXSAV.GT.0) THEN
        
C       We want to write it to disc on level IFPSAV, i.e. at
C       level:

        I = MAX(NLMIN,MIN(NLMAX,ABS(IFXSAV)))
        
C       Calculate the interpolated streamfunction, if 
C       it's not calculated

        IF (LISO.EQ.0) THEN
          CALL XU2ISO (DUP(1),DUP(1+1*NEQU),TRIAS(1,NLMAX),LISO)
        END IF

C       Write out the values in the vertices for the streamfunction.
C       Make use of the 2-level ordering to write out only level I.
C       The sign of IFXSAV decides on whether the output is in
C       formatted or unformatted form.

        CALL STPUT (LFNAMS(5), CFN)
        CALL PPWRVC (0,TRIAS(ONVT,I),DWORK(L(LP)),-1,
     *               IFXSAV,CFN)

      END IF
     
C     -------------------------------------------------------------
C     That's it, all solution vectors written out.
C
C     =================================================================
C     The writing of the solutions is completed here.
C     We go on now to check if we have to write out GMV files.
C     The way of doing this is very similar to the writing of
C     solution vectors...
C
C     At first check if we should write out GMV files at all:

      IF (IGMV.GT.0) THEN
      
C       We want to write the GMV file to disc on level IGMV, i.e. at
C       level:

        I = MAX(NLMIN,MIN(NLMAX,IGMV))

C       -------------------------------------------------------------
C       Use handle 69 for writing the GMV file
        
        K = 69
        
C       Open it

        CALL STPUT (LFNAMS(6), CFN)
        CALL AUXCFN (CFN,'',CNAME,-1)
        CALL GMVOF0 (K,-1,CNAME)
        
        IF (IER.EQ.0) THEN
        
C         Write header and triangulation on level I.
C         Obtain NCELLS and NVERTS.

          CALL GMVHEA (K)
          CALL GMVTRI (K,TRIAS(1,I),1,NCELLS,NVERTS)
        
C         Write materials;
C         we don't give material names for now.

          CALL GMVMAT (K,TRIAS(1,I),0,NCELLS,TRIAS(ONEL,I),
     *                 -16,'',DFGMMC,IGEOM,DGEOM)
          CALL GMVMAT (K,TRIAS(1,I),1,NVERTS,
     *                 TRIAS(ONVT,I)+TRIAS(ONMT,I),
     *                 -16,'',DFGMMV,IGEOM,DGEOM)
        
C         Calculate the velocity field, pressure and the stream 
C         function if they are not calculated already:

          IF ((LU.EQ.0).OR.(LV.EQ.0)) THEN
          
            CALL XINTUV (DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),LU,LV)
          
C           Implement boundary conditions, since the interpolation
C           does not correctly handle boundaries:

            KCORVG = L(TRIAS(OLCORVG,NLMAX))
            KXNPR  = L(TRIAS(OLXNPR,NLMAX))
            CALL BDRCOR (DWORK(L(LU)),DWORK(L(LV)),
     *                   TRIAS(1,NLMAX),DWORK(KCORVG),
     *                   KWORK(KXNPR),UE,1D0,DASMBL(ORE),
     *                   IASMBL,DASMBL,IGEOM,DGEOM)
          END IF

          IF (LP.EQ.0) THEN
            CALL XINTPV (DUP(1+2*NEQU),TRIAS(1,NLMAX),LP)
          END IF

          IF (LISO.EQ.0) THEN
            CALL XU2ISO (DUP(1),DUP(1+1*NEQU),TRIAS(1,NLMAX),LISO)
          END IF

C         Write all these at the desired level to the GMV file:

          CALL GMVVEL (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *                 DWORK(L(LU)),DWORK(L(LV)))
          CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *                 DWORK(L(LP)),'pressure')
          CALL GMVSCA (K,TRIAS(1,I),0,NCELLS,TRIAS(ONEL,I),
     *                 DUP(1+2*NEQU),'pressure')
          CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *                 DWORK(L(LISO)),'streamfunction')
        
C         Calculate the H1-error to LERR

          IF (LERR.EQ.0) THEN
            CALL ZNEW(TRIAS(ONVT,NLMAX),1,LERR,'DERR  ')
          END IF

          CALL ERPQH1(DWORK(L(LU)),DWORK(L(LV)),TRIAS(ONVT,NLMAX),
     *                TRIAS(1,NLMAX),DWORK(L(LERR)))
          CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *                 DWORK(L(LERR)),'H1error')
          CALL ERPCH1(DWORK(L(LP)),TRIAS(ONVT,NLMAX),TRIAS(1,NLMAX),
     *                DWORK(L(LERR)))
          CALL GMVSCA (K,TRIAS(1,I),1,NVERTS,TRIAS(ONVT,I),
     *                 DWORK(L(LERR)),'H1errorP')
        
C         Write geometry:

          CALL FBDGMV (K,2,IGEOM,DGEOM)
          
C         Write the footer, finish.

          CALL GMVFOT (K)
          CLOSE (K)
        
C         -------------------------------------------------------------
C         That's it, GMV file is written out.
      
        END IF ! IER = 0
          
      END IF ! IGMV > 0
      
C     =================================================================
C     Error analysis with analytic function.
C
C     Call the error analysis routine to print current errors
C     to the terminal. This is done only for accepted solutions.

      IF (IERANA.NE.0) THEN
        CALL ERATRM (MFILE,IERANA,VECDAT(1,NLMAX),DUP,DRHS,DAUX,
     *               NLMIN,NLMAX,TRIAS,1D0,
     *               IASMBL,DASMBL,IGEOM,DGEOM,
     *               UE,PE,UEX,UEY)
      END IF
      
C     Release any memory we might have used

      IF (LERR.NE.0) CALL ZDISP(0,LERR,'DERR  ')
      IF (LISO.NE.0) CALL ZDISP(0,LISO,'DISO  ')
      IF (LU.NE.0) CALL ZDISP(0,LU,'DU    ')
      IF (LV.NE.0) CALL ZDISP(0,LV,'DV    ')
      IF (LP.NE.0) CALL ZDISP(0,LP,'DP    ')
      
C     Finally stop the time we need for postprocessing:

      CALL GTMAUX (TIMNG,DSTPAR,OTNLTIM-1+OTTPOST,1)
      
C     Add the time needed for postprocessing to the "total time".
C     Since this routine is not called in the solver (is't no
C     callback routine of the solver), if the time for the
C     postprocessing should calculate to the local time, we
C     have to add it manually...

      DSTPAR(OTMTOT) = DSTPAR(OTMTOT) + DSTPAR(OTNLTIM-1+OTTPOST)
      
      END

************************************************************************
* Auxiliary routine: Create filename
*
* This routine performs the string manimulation:
*    DEST = NAME + AUX + string(IDX) (if IDX <> -1)
* NAME, AUX and DEST must be strings, DEST must be a string (60).
************************************************************************

      SUBROUTINE AUXCFN (NAME,EXT,DEST,IDX)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'dstrings.inc'
      
      INTEGER IDX
      CHARACTER NAME*(*),EXT*(*),DEST*(60)
      
      INTEGER LSTR
      
      LSTR = STNEWC (.TRUE.,NAME)

      CALL STCATC (LSTR,.TRUE.,EXT)
      IF (IDX.GT.-1) THEN
        CALL STCATI (LSTR,IDX,0,.FALSE.)
      END IF

      CALL STPUT (LSTR, DEST)
      CALL STDIS (LSTR)
      
      END
      
************************************************************************
* Standard time output for nonstationary solver
*
* This routine prints out the timing statistics on the nonstationary
* solver to terminal/file.
*
* In:
*   MSHOW  : Level of output
*   TMINI  : Time needed for preprocessing
*   TMTOT  : Total time of the algorithm including initialization/
*            postprocessing/everything.
*   IPARAM : array [1..SZNSDI] of integer
*   DPARAM : array [1..SZNSDD] of double
*            Integer and double prec. parameter blocks that define the
*            behaviour of the stationary solver. 
*
* Out:
*   If MSHOW>=0, the output is written to the file.
*   If MSHOW>=2, the output is written to the standard terminal
************************************************************************

      SUBROUTINE STATST (MSHOW,TMINI,TMTOT,ISTPAR,DSTPAR)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'

      INCLUDE 'ssolvers.inc'
      
      INCLUDE 'stiming.inc'
      INCLUDE 'snsdef.inc'
      
      INTEGER MSHOW,ISTPAR(*)
      DOUBLE PRECISION DSTPAR(*),TMTOT,TMINI
      
      INTEGER IOFS
      DOUBLE PRECISION TTMG,T
      CHARACTER CSTR*(255)
      
      IOFS = OTNLTIM-1
      
      IF (MSHOW.GT.0) THEN
        WRITE(CSTR,'(A)')
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A)') 'STATISTICS :'
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

        WRITE(CSTR,'(A,I12)') 'NWORK :      ',NWORK
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,I12)') 'IWORK :      ',IWORK
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,I12)') 'IWMAX :      ',IWMAX
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,1)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C       DSTPAR(TMTOT) in the solver structure shows the total
C       time of the solver + postprocessing. But as we want to
C       show the actual total time for everything, we print
C       out TMTOT from the parameter, not from DSTPAR!

        WRITE(CSTR,'(A,F20.10)') 'total time       : ',TMTOT
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
        WRITE(CSTR,'(A,F20.10)') 'init. time       : ',TMINI
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
        T = DSTPAR(IOFS+OTTPOST)+
     *      DSTPAR(IOFS+OTTADF ) + DSTPAR(IOFS+OTTUPW )+
     *      DSTPAR(IOFS+OTTBDR ) + DSTPAR(IOFS+OTTLC  )+
     *      DSTPAR(IOFS+OTTLSOL)
        
        WRITE(CSTR,'(A,F20.10)') 'appr. time       : ', T
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

        T = DSTPAR(IOFS+OTTADF )+DSTPAR(IOFS+OTTUPW )+
     *      DSTPAR(IOFS+OTTBDR )+DSTPAR(IOFS+OTTLC  )
     
        WRITE(CSTR,'(A,F20.10)') '-> lin.   time   : ', T
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
        WRITE(CSTR,'(A,F20.10)') '   -> mavec time : ',
     *                           DSTPAR(IOFS+OTTADF)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
        WRITE(CSTR,'(A,F20.10)') '   -> konv. time : ',
     *                           DSTPAR(IOFS+OTTUPW)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
        WRITE(CSTR,'(A,F20.10)') '   -> bdry  time : ',
     *                           DSTPAR(IOFS+OTTBDR)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
        WRITE(CSTR,'(A,F20.10)') '   -> LC    time : ',
     *                           DSTPAR(IOFS+OTTLC)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

        WRITE(CSTR,'(A,F20.10)') '-> l.sol. time   : ',
     *                            DSTPAR(IOFS+OTTLSOL)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

        WRITE(CSTR,'(A,F20.10)') '-> post.  time   : ', 
     *                            DSTPAR(IOFS+OTTPOST )
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

        WRITE(CSTR,'(A,I12)') '#nonlinear iterations : ', ISTPAR(OITE)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,I12)') '#linear iterations    : ', ISTPAR(ONLIN)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,1)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        
        TTMG = DSTPAR(IOFS+OTTMG)
        IF (TTMG.EQ.0D0) TTMG=1D0

        WRITE(CSTR,'(A)')
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A)') ' MULTIGRID COMPONENTS [in percent]:'
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') ' smoothing     :', 
     *                 1.D2*DSTPAR(IOFS+OTTSMTH)/TTMG
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') ' solver        :', 
     *                 1.D2*DSTPAR(IOFS+OTTCGC )/TTMG
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') ' defect calc.  :', 
     *                 1.D2*DSTPAR(IOFS+OTTDEF )/TTMG
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') ' prolongation  :', 
     *                 1.D2*DSTPAR(IOFS+OTTPROL)/TTMG
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,F20.10)') ' restriction   :', 
     *                 1.D2*DSTPAR(IOFS+OTTREST)/TTMG
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,1)
        CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      END IF        

1     FORMAT(79('-'))

      END
