************************************************************************
* (Re)initialise the Navier Stokes solver for a new run.
*
* In:
*  BPRT  - false=perform a complete reinitialisation
*          true =perform a partial reinitialisation, don't destroy
*                current solution vectors
*
* IGAMDE - perform grid adaption:
*          =0: no grid adaption
*          =1: grid adaption according to information in the DAT file;
*              should be done if the geometry has taken a major change
*          =2: "adaptive" grid adaption, only on level NLMAX-1;
*              to correct the grid if the geometry changed only slightly
*
* Generally the routine performs the following:
*  - regeneration of the mass matrices, Stokes matrices, RHS
*  - eventually grid deformation
*  - precalculate information from fictitious boundary information
*    about distances,...
************************************************************************

      SUBROUTINE ININSO (BPRT,IGAMDE)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cns.inc'
      INCLUDE 'cnsparfrac.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'ciniopt.inc'
      INCLUDE 'cinigeometry.inc'
      
      INCLUDE 'cgridadapt.inc'
      
C parameters

      LOGICAL BPRT
      INTEGER IGAMDE
      
C local variables

      INTEGER ILV
      INTEGER LF1,LF2,LF12P

C Reinit all Navier-Stokes structures

      WRITE (MTERM,9000) 
      WRITE (MTERM,'(A,F10.4)') 'Cylinder at position: DCXPOS=',DCXPOS
      WRITE (MTERM,'(A,F10.4)') '                      DCYPOS=',DCYPOS
      WRITE (MTERM,'(A,F10.4)') 'Radius of cylinder:   DCRAD =',DCRAD
      WRITE (MTERM,9000) 
      WRITE (MTERM,'(A)') 
      
      WRITE (MTERM,'(A)') 
      WRITE (MTERM,9000) 
      WRITE (MTERM,'(A,I10)') 
     *     'Number of opt. iteration:  OPTCNT = ',OPTCNT
      WRITE (MTERM,'(A,I10)') 
     *     'Type of obstacle:          ICIRTP = ',ICIRTP
      WRITE (MTERM,'(A,F10.4,A2,I5,A1)') 
     *     'Obstacle at position:      DCXPOS = ',DCXPOS,' (',IXNUM,')'
      WRITE (MTERM,'(A,F10.4,A2,I5,A1)') 
     *     '                           DCYPOS = ',DCYPOS,' (',IYNUM,')'
      WRITE (MTERM,'(A,F10.4)') 
     *     'X-Radius of obstacle:      DCRADX = ',DCRADX+DCRRLX
      WRITE (MTERM,'(A,F10.4)') 
     *     'Y-Radius of obstacle:      DCRADY = ',DCRADY+DCRRLY
      WRITE (MTERM,'(A,F10.4,A2,I5,A1)') 
     *     'Rotation of obstacle:      DCROT  = ',DCROT,' (',IRNUM,')'
      WRITE (MTERM,'9000') 
      WRITE (MTERM,'(A)') 

C Release the mass-/Stokes matrix entries,...

      IF (ISTAT.NE.0) THEN
        CALL DISMM (0)
      END IF
      
      IF (IPRECA.NE.4) THEN
        CALL DISSTM (0)
      ENDIF
      
      CALL DISMTB (0)

C Release any precalculated information from the grid, since the
C geometry might have changed for the new run:

      CALL DISMGP ()

C Perform a grid deformation/optimisation,...
C (As the information about elements meeting on each vertex is already 
C calculated in GENORS, we don't have to call XS2V again to calculate 
C these information)

      IF (IGAMDE.EQ.1) THEN
        
C Full grid adaption
        
        CALL XSMMGR (0,0,0)
        
      ELSE IF (IGAMDE.EQ.2) THEN

C "adaptive" grid adaption, adaption only on level NLMAX-1

C        CALL XSMMGR (1,-1,-1)

C ok, sometimes with the ellipse the above line makes the program crash!
C For now we switch back to complete re-adaption:

        CALL XSMMGR (0,0,0)

      END IF
      
C The geometry is fixed now.
C
C Try to precalculate information from all grids to save computational
C time later -- e.g. distance from all points to the boundary,...

C To do this, at first make a copy of the current information of the
C geometries on all levels:

      CALL GENMGT ()
        
C And then use this to precalculate what can be precalculated:

      CALL GENMGP ()

C After the grid deformation is completed, restore the KNPR-arrays
C and rebuild them. It's crucial that this is done AFTER the grid
C deformation in order to mark the correct inner nodes as fictitious
C boundary nodes!!!

      DO ILV=NLMIN,NLMAX
        CALL XLCP3(KLNPRO(ILV),KLNPR(ILV),KNVT(ILV)+KNMT(ILV))
        CALL IMPBDG (ILV)
      END DO
        
C ... and rebuild the mass/Stokes matrices, depending on our
C current grid

      IF (IPRECA.NE.4) THEN
        CALL GENSTM 
      ENDIF

      IF (ISTAT.NE.0) THEN
        CALL GENMM (IMASS.EQ.0,IMASSL.NE.0)
      END IF

      CALL GENMTB(IPRECB.GE.2)
          
C Rebuild right hand side like in INIT1 for finest level

      CALL GENRHS(NLMAX,LF1,LF2)

      LF12P = KLF12P(NLMAX)
      CALL  LCP1 (DWORK(L(LF1)), DWORK(L(LF12P)),KNEQA(NLMAX))
      CALL  LCP1 (DWORK(L(LF2)), DWORK(L(LF12P)+KNEQA(NLMAX)), 
     *                  KNEQA(NLMAX))
      CALL  ZDISP (0, LF1, 'LF1TMP' )
      CALL  ZDISP (0, LF2, 'LF2TMP' )

C Immediately implement boundary Dirichlet/Neumann conditions.

      CALL IMPSLR (NLMAX,(IBDR.EQ.1).AND.(ISTAT.EQ.0),
     *                 .TRUE.,.FALSE.)
     
C Reinitialise the current solution vectors if necessary

      IF (.NOT.BPRT) THEN
        CALL LCL1 (DWORK(L(KLUP(NLMAX))),KNUP(NLMAX))
      END IF
     
C Implement boundary conditions into solution. This is always necessary,
C as for a changed geometry we have at least different Dirichlet solution
C values. Although this would not be necessary in the first iteration
C or if the solution vectors contain precalculated solutions,
C we do this here for "safetyness" and "simplicity", since it's quickly
C be done...

      CALL IMPSLR (NLMAX,.FALSE.,.FALSE.,.TRUE.)

9000  FORMAT(79('-'))

      END
      
************************************************************************
* Navier-Stokes calculation
*
* This routine performs the following:
* - Initialisation of the Navier-Stokes-solver
* - Start of the Navier-Stokes solver
* - Postprocessing of the Navier-Stokes-solver.
*   If specified in the DAT file, GMV-files with unique name will be
*   written out. 
* - Calculation of the values of interest.
*
* In:
*  IVOI   - Values of interest to calculate:
*           1=calculate drag   -> DINFO(1)
*                       lift   -> DINFO(2)
*                       volume -> DINFO(3)
*  CNFIG - Type of reinitialisation/postprocessing:
*           Bit 0: true=reinitialise start vectors;
*                  otherwise use current solution as start vectors
*           Bit 1: true=rebuild problem specific information
*                  (matrix,rhs,...) according to current geometry
*           Bit 2: true=perform grid adaption to capture current
*                  geometry
*           Bit 3: true=perform, only partial grid adaption for
*                  fast-adaption of the grid, e.g. if the geometry
*                  has changed only slightly
*           Bit 4: Show results of postprocessing on screen/
*                  write into file (according to parameter MT in
*                  the COMMON block)
*           Bit 5: true=build unique GMV-filename for output with
*                  the help of level, IXNUM, IYNUM, IRNUM
*                  (variables in COMMON-blocks
*                  Otherwise: build filename with OPTCNT-counter
*           Bit 6: true=write out "optimisation" GMV-file additionally
*                  to standard GMV-file
*           Typical values:
*            5=complete reinitialisation with full grid adaption,
*              no screen output, no additional GMV-output
*            7=complete reinitialisation with full grid adaption,
*              no screen output, no additional GMV-output
*            21=standard computation without reinitialisation,
*               output on screen, no additional GMV-output
*            23=complete reinitialisation with full grid adaption,
*               output on screen, no additional GMV-output
*           119=complete reinitialisation with full grid adaption,
*               output on screen, unique GMV-output
*           127=standard computation without reinitialisation,
*               output on screen, unique GMV-output
*
* Out:
*  DINFO  - array [1..2] of double
*           calculated double-precision values of interest
*  IINFO  - array [1..2] of double
*           calculated integer values of interest
*  TPRE   - time needed for the preprocessing is added to this
*  TSOLV  - time needed for solving is added to this
*  TPOS   - time needed for postprocessing is added to this
************************************************************************

      SUBROUTINE NSCALC (IVOI,CNFIG,TPRE,TSOLV,TPOS,DINFO,IINFO)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cns.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      
      INCLUDE 'ciniopt.inc'
      INCLUDE 'cfiles.inc'
      INCLUDE 'dstrings.inc'
      
C parameters
      
      INTEGER IVOI,CNFIG,IINFO(*)
      DOUBLE PRECISION TPRE,TSOLV,TPOS,DINFO(*)
      
C externals: exact solution (for error analysis)
      
      DOUBLE PRECISION UE
      EXTERNAL UE

C routine for minimising the volume efficiency
      
      DOUBLE PRECISION MINVEF
      EXTERNAL MINVEF
      
C local variables

      DOUBLE PRECISION T1,T2,DDRAG,DLIFT
      INTEGER I,ITFILM,IFILEN,MSHOW
      CHARACTER*128 CFNAME

      INTEGER IEL,IOPTSS,IFPOST,IMETH
      DOUBLE PRECISION DRAG(3),LIFT(3),VOL(3)
      DOUBLE PRECISION CELLSZ,CAREA

C We only support drag/lift-calculation here for now.
C IINFO is not used for now.

      IF (IVOI.NE.1) RETURN
      
C Calculation method for drag/lift:

      IMETH=2

C Output control

      MSHOW=0
      IF (AND(CNFIG,16).NE.0) THEN
        MSHOW=MT
      END IF
      
C ---------------------------------------------------------------------
      
      IF ((IOPTCC.EQ.0).OR.(IOPTCC.EQ.1).OR.(IOPTCC.EQ.3)) THEN

C We perform drag/lift calculation with only one NS-calculation.

C Minimise the volume efficiency if necessary

        IF (IOPTCC.EQ.3) THEN
      
C Before we can do minimising the volume-efficiency, we should adapt
C our grid:

          IF (AND(CNFIG,4).NE.0) THEN
            IF (AND(CNFIG,8).NE.0) THEN
C "adaptive" grid adaption, adaption only on level NLMAX-1
              CALL XSMMGR (1,-1,-1)
            ELSE
C Full grid adaption
              CALL XSMMGR (0,0,0)
            END IF
          END IF
      
C Try to minimise the volume efficiency; minimum number >= 1.
C Modify the relaxation in 1/3 the size of a cell

          CALL CMINCZ (NLMAX,IEL,CAREA,CELLSZ)      

          DCRRLX = MINVEF (CELLSZ/100D0,CELLSZ/10D0,CELLSZ)
        END IF

C Reinitialise the solver if necessary.
C If the volume-efficiency was minimised before, don't perform grid-
C adaption, that was done previously.
        
        CALL ZTIME (T1)

        IF ((AND(CNFIG,6).NE.0).OR.(IOPTCC.EQ.3)) THEN
          I=0
          IF (IOPTCC.NE.3) THEN
            IF (AND(CNFIG,4).NE.0) I=1
            IF (AND(CNFIG,12).EQ.12) I=2
          END IF
          CALL ININSO (AND(CNFIG,1).EQ.0,I)
        END IF
        
        CALL ZTIME (T2)
        TPRE = TPRE+T2-T1

C Start the nonlinear solver to calculate the solution

        ITFILM=0
        IFILEN=0

        CALL ZTIME (T1)

        CALL NONSTL (MFILE1,MT,ITFILM,IFILEN)
        IF (IER.NE.0) RETURN

        CALL ZTIME (T2)
        TSOLV = TSOLV+T2-T1

C perform the postprocessing

        CALL ZTIME (T1)

C Build a filename for the GMV-output:

        IF (AND(CNFIG,32).NE.0) THEN
          I = STNEWC (.FALSE.,'gmv/OPT.gmv.')
          CALL STCATI (I,ILEV,4,.TRUE.)
          CALL STCATC (I,.FALSE.,'.')
          CALL STCATI (I,IXNUM,4,.TRUE.)
          CALL STCATC (I,.FALSE.,'.')
          CALL STCATI (I,IYNUM,4,.TRUE.)
          CALL STCATC (I,.FALSE.,'.')
          CALL STCATI (I,IRNUM,4,.TRUE.)
          CALL STPUT (I, CFNAME)
          CALL STDIS (I)
        ELSE
          I = STNEWC (.FALSE.,'gmv/OPT.gmv.')
          CALL STCATI (I,OPTCNT,4,.TRUE.)
          CALL STPUT (I, CFNAME)
          CALL STDIS (I)
        END IF

C Perform simple postprocessing

        IFILEN=0
        IFPOST=0
        I=0
        IF (AND(CNFIG,64).NE.0) I=1

        CALL FPOST(IFPOST,IFILEN,ITFILM,UE,MSHOW,DDRAG,DLIFT,I,CFNAME)
      
C Write the values of interest to the double array
      
        IF (IOPTCC.NE.1) THEN

C Recaulculate drag, lift and volume

          CALL DLCALC (NLMAX,IMETH,.TRUE.,DINFO(1),DINFO(2),DINFO(3))

        ELSE IF (IOPTCC.EQ.1) THEN
        
C Drag/Lift calculation with radius relaxation.
C Recalculate drag/lift by interpolation between 3 nearby radii

C Minimum cell size on fict. bdry., divided by 3

          CALL CMINCZ (NLMAX,IEL,CAREA,CELLSZ)      
          CELLSZ = CELLSZ / 3D0  
        
          DCRRLX = CELLSZ
          DCRRLY = CELLSZ
          CALL DLCALC (NLMAX,IMETH,.TRUE.,DRAG(1),LIFT(1),VOL(1))

          DCRRLX = 2*CELLSZ
          DCRRLY = 2*CELLSZ
          CALL DLCALC (NLMAX,IMETH,.TRUE.,DRAG(2),LIFT(2),VOL(2))

          DCRRLX = 0D0
          DCRRLY = 0D0
          CALL DLCALC (NLMAX,IMETH,.TRUE.,DRAG(3),LIFT(3),VOL(3))
        
C take the mean values
        
          DINFO(1) = (DRAG(1)+DRAG(2)+DRAG(3))/3D0
          DINFO(2) = (LIFT(1)+LIFT(2)+LIFT(3))/3D0
          DINFO(3) = (VOL(1)+VOL(2)+VOL(3))/3D0
          
          WRITE (MTERM,'(A,3D16.8)') 'I(FORCE) VOL CP interp: ',
     *                               DINFO(1),DINFO(2),DINFO(3)

        
        END IF

C ---------------------------------------------------------------------

      ELSE IF (IOPTCC.EQ.2) THEN

C We perform drag/lift calculation with 3 subsequent NS-calculations
C and obtain drag/lift-values by taking the mean.

C At first we adapt our grid for the rest of the computation here:

        IF (AND(CNFIG,4).NE.0) THEN
          IF (AND(CNFIG,8).NE.0) THEN
C "adaptive" grid adaption, adaption only on level NLMAX-1
            CALL XSMMGR (1,-1,-1)
          ELSE 
C Full grid adaption
            CALL XSMMGR (0,0,0)
          END IF
        END IF

C Perform 3 substeps with 3 different nearby radii without re-adapting
C the grid. Calculate drag/lift by taking the mean value.

C minimum cell size on fict. bdry., divided by 2:

        CALL CMINCZ (NLMAX,IEL,CAREA,CELLSZ)      
        CELLSZ = CELLSZ / 2D0  

        DO IOPTSS = 1,3
    
          IF (IOPTSS.EQ.1) THEN 
          
C no relaxation
          
            DCRRLX = 0D0
            DCRRLY = 0D0
    
          ELSE IF (IOPTSS.EQ.2) THEN 

C increase the radius a little bit
          
            DCRRLX = CELLSZ
            DCRRLY = CELLSZ
            
          ELSE IF (IOPTSS.EQ.3) THEN 

C decrease the radius a little bit
        
            DCRRLX = -CELLSZ
            DCRRLY = -CELLSZ
            
          END IF
        
C Reinitialise the solver if necessary, don't adapt the grid (since this
C is already done).
C Start with a 0-solution only in the first step and if the
C caller desires that. The following 2 calculations are then performed
C with the previous result as start vector.
        
          CALL ZTIME (T1)

          IF ((AND(CNFIG,6).NE.0).OR.(IOPTSS.GT.1)) THEN
            IF (IOPTSS.EQ.1) THEN
              CALL ININSO (AND(CNFIG,1).EQ.0,0)
            ELSE
              CALL ININSO (.TRUE.,0)
            END IF
          END IF
          
          CALL ZTIME (T2)
          TPRE = TPRE+T2-T1
        
C Start the solver
        
          CALL ZTIME (T1)
          CALL NONSTL (MFILE1,MSHOW,ITFILM,IFILEN)
          IF (IER.NE.0) RETURN
          CALL ZTIME (T2)
          TSOLV = TSOLV+T2-T1

C Build a filename for the GMV-output:

          IF (AND(CNFIG,32).NE.0) THEN
            I = STNEWC (.FALSE.,'gmv/OPT.gmv.')
            CALL STCATI (I,ILEV,4,.TRUE.)
            CALL STCATC (I,.FALSE.,'.')
            CALL STCATI (I,IXNUM,4,.TRUE.)
            CALL STCATC (I,.FALSE.,'.')
            CALL STCATI (I,IYNUM,4,.TRUE.)
            CALL STCATC (I,.FALSE.,'.')
            CALL STCATI (I,IRNUM,4,.TRUE.)
            CALL STCATC (I,.FALSE.,'.')
            CALL STCATI (I,IOPTSS,4,.TRUE.)
            CALL STPUT (I, CFNAME)
            CALL STDIS (I)
          ELSE
            I = STNEWC (.FALSE.,'gmv/OPT.gmv.')
            CALL STCATI (I,OPTCNT,4,.TRUE.)
            CALL STCATC (I,.FALSE.,'.')
            CALL STCATI (I,IOPTSS,4,.TRUE.)
            CALL STPUT (I, CFNAME)
            CALL STDIS (I)
          END IF

C Perform simple postprocessing

          IFILEN=0
          IFPOST=0
          ITFILM=0
          IFILEN=0
          MSHOW=0
          IF (AND(CNFIG,32).NE.0) THEN
            MSHOW=MT
          END IF
          I=0
          IF (AND(CNFIG,64).NE.0) I=1

          CALL FPOST(IFPOST,IFILEN,ITFILM,UE,MSHOW,DDRAG,DLIFT,I,
     *               CFNAME)
        
C recalculate drag/lift to additionally calculate the volume:

          CALL DLCALC (NLMAX,IMETH,.TRUE.,
     *                 DRAG(IOPTSS),LIFT(IOPTSS),VOL(IOPTSS))

          CALL ZTIME (T2)
          TPOS = TPOS+T2-T1

        END DO

C calculate the mean values and return these:

        DINFO(1) = (DRAG(1)+DRAG(2)+DRAG(3))/3D0
        DINFO(2) = (LIFT(1)+LIFT(2)+LIFT(3))/3D0
        DINFO(3) = (VOL(1)+VOL(2)+VOL(3))/3D0

        WRITE (MTERM,'(A,3D16.8)') 'I(FORCE) VOL CP interp 1: ',
     *                                    DRAG(1),LIFT(1),VOL(1)
        WRITE (MTERM,'(A,3D16.8)') 'I(FORCE) VOL CP interp 2: ',
     *                                    DRAG(2),LIFT(2),VOL(2)
        WRITE (MTERM,'(A,3D16.8)') 'I(FORCE) VOL CP interp 3: ',
     *                                    DRAG(3),LIFT(3),VOL(3)
        WRITE (MTERM,'(A,3D16.8)') 'I(FORCE) VOL CP interp. : ',
     *                                    DINFO(1),DINFO(2),DINFO(3)

      END IF

      END

***********************************************************************
* Calculate the minimum cell size
*
* This routine searches for the element with the smallest area.
* The result will be the size of the element (to be more exact:
* the minimum length of an edge of the element).
*
* In:
*  ILV   - The level
*
* Out:
*  IEL   - The number of the element
*  DAREA - The area of this element
*  CZ    - The minimum edge length
***********************************************************************

      SUBROUTINE CMINCZ (ILV,IEL,DAREA,CZ)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
C parameters
      
      INTEGER ILV, IEL
      DOUBLE PRECISION CZ, DAREA
      
C local variables

      INTEGER LAREA,I
      DOUBLE PRECISION E1,E2,E3,E4

C search for the element with the smallest area
      
      LAREA = L(KLAREA(ILV))
      
      DAREA = DWORK(LAREA)
      IEL = 1
      DO I=2,KNEL(ILV)
        IF (DWORK(LAREA+I-1).LT.DAREA) THEN
          IEL = I
          DAREA = DWORK(LAREA+I-1)
        END IF
      END DO
      
C get the length of the edges of this element

      CALL CELENS (IEL,DWORK(L(KLCVG(ILV))),KWORK(L(KLVERT(ILV))),
     *             KNEL(ILV),E1,E2,E3,E4)   
     
      CZ = MIN ( MIN (E1,E2) , MIN (E3,E4) )   
      
      END
      
***********************************************************************
* Calculate the length of the four edges belonging to
* Element IEL.
*
* In:
*  IEL    - Number of the element
*  NEL    - Total number of elements in the grid
*  DCORVG - Coordinates of the vertices
*  KVERT  - Vertices of the elements
*
* Out: 
*  E1,E2,E3,E4 - The four length's of the edges
***********************************************************************
      
      SUBROUTINE CELENS (IEL,DCORVG,KVERT,NEL,E1,E2,E3,E4)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      DOUBLE PRECISION DCORVG(2,*),E1,E2,E3,E4
      INTEGER IEL,NEL,KVERT(NNVE,NEL)
      
C parametrisation of the domain

      DOUBLE PRECISION PARX,PARY
      EXTERNAL PARX,PARY
      
C local variables

      DOUBLE PRECISION X1,X2,X3,X4,Y1,Y2,Y3,Y4
      
C get the coordinates

      X1 = DCORVG(1,KVERT(1,IEL))
      X2 = DCORVG(1,KVERT(2,IEL))
      X3 = DCORVG(1,KVERT(3,IEL))
      X4 = DCORVG(1,KVERT(4,IEL))

      Y1 = DCORVG(2,KVERT(1,IEL))
      Y2 = DCORVG(2,KVERT(2,IEL))
      Y3 = DCORVG(2,KVERT(3,IEL))
      Y4 = DCORVG(2,KVERT(4,IEL))
      
C Calculate E1,...,E4

      E1 = DSQRT ( ( X1-X2 )**2 + ( Y1-Y2 )**2 )
      E2 = DSQRT ( ( X2-X3 )**2 + ( Y2-Y3 )**2 )
      E3 = DSQRT ( ( X3-X4 )**2 + ( Y3-Y4 )**2 )
      E4 = dSQRT ( ( X4-X1 )**2 + ( Y4-Y1 )**2 )
     
      END

***********************************************************************
* Minimise volume efficiency
*
* This routine tries to modify the radius relaxation of the fictitious
* boundary obstacle until the volume efficiency index is as small
* as possible >= 1.
*
* In:
*  SSIZE  - The step size in that the radius relaxation is allowed
*           to change
*  STMI   - minimum absolute relaxation that is checked for the radius
*  STMX   - maximum absolute relaxation that for allowed in the radius
* 
* Return: The new guess for the radius relaxation.
*  The result will be in -STMX..STMX. At least the interval -STMN..STMN
*  is checked for a minimum of the volume efficiency.
***********************************************************************
      
      DOUBLE PRECISION FUNCTION MINVEF (STSZ,STMI,STMX)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'

      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'ciniopt.inc'
      INCLUDE 'cinigeometry.inc'
      
C parameters

      DOUBLE PRECISION STSZ,STMI,STMX

C externals

      DOUBLE PRECISION FBDVOL
      EXTERNAL FBDVOL

C local variables
      
      DOUBLE PRECISION DVOL, DRVOL, DVEFF1, DRAD, DTMP
      DOUBLE PRECISION DGUESS,DVEFFB,DVEFFI
      INTEGER NTESTS
      
      NTESTS = 1
      
C make a backup of the current radius relaxation

      DRAD = DCRRLX
      
C current guess; start with guess=0

      DCRRLX = 0D0
      
C guess for best relaxation

      DGUESS = 0D0
      
C calculate our initial volume; cancel if the fict. boundary
C implementation doesn't support the analytical calculation
C of the volume.

      DRVOL = FBDVOL (0)
      IF (DVOL.EQ.-1D0) GOTO 99999

C calculate the approximate volume and the volume efficiency

      CALL DLCALC(NLMAX,0,.FALSE.,DTMP,DTMP,DVOL)
      DVEFFB = DVOL/DRVOL
      DVEFF1 = DVEFFB
      DVEFFI = DVEFFB
      
C when < 1, we have to increase until we are >= 1; we do this by
C increasing the radius in hope to increase the vol.eff. as well

      IF (DVEFFB.LT.1D0) THEN

C Increase the radius and recalculate our volume efficiency,
C if we are allowed to do that.
C Cropped DO-WHILE-loop

        DCRRLX = -STMI-STSZ
        
1       IF ( (DCRRLX.LE.STMI).OR.
     *       ( (DCRRLX.LE.STMX) .AND. (DVEFF1.LT.1D0) ) ) THEN
          DCRRLX=DCRRLX+STSZ
          DRVOL = FBDVOL (0)
          CALL DLCALC(NLMAX,0,.FALSE.,DTMP,DTMP,DVOL)
          DVEFF1 = DVOL/DRVOL
          NTESTS=NTESTS+1

C If we can come closer to 1, we take this value.
C If we have been able to get the first number > 1, we take that;
C otherwise we take the smallest value < 1 because there is no other
C chance...
          IF (ABS(1D0-DVEFF1).LT.ABS(1D0-DVEFFB)) THEN
            IF ((DVEFF1.GE.1D0).OR.(DVEFFB.LT.1D0)) THEN
              DVEFFB=DVEFF1
              DGUESS=DCRRLX
            END IF
          END IF
          GOTO 1
        END IF

      ELSE

C Decrease the radius (and hopefully the volume efficiency) until we obtain
C the smallest vol.eff. >= 1

        DCRRLX = STMI+STSZ

2       IF ((DCRRLX.GE.-STMI).OR.
     *      ((DCRRLX.GE.-STMX).AND.(DVEFF1.GE.1D0)) ) THEN
          DCRRLX=DCRRLX-STSZ
          DRVOL = FBDVOL (0)
          CALL DLCALC(NLMAX,0,.FALSE.,DTMP,DTMP,DVOL)
          DVEFF1 = DVOL/DRVOL
          NTESTS=NTESTS+1

C If we can come closer to 1 (without being < 1), we take this value.

          IF ( (ABS(1D0-DVEFF1).LT.ABS(1D0-DVEFFB)) .AND.
     *         (DVEFF1.GE.1D0) ) THEN
            DVEFFB=DVEFF1
            DGUESS=DCRRLX
          END IF
          GOTO 2
        END IF

      END IF

99999 MINVEF = DGUESS
      DCRRLX = DRAD
      
      WRITE (MTERM,9000) 
      WRITE (MTERM,*) 'Init. Vol.eff.        = ',DVEFFI
      WRITE (MTERM,*) 'Relaxation            = ',DGUESS
      WRITE (MTERM,*) 'Radius                = ',DCRADX+DGUESS
      WRITE (MTERM,*) 'Vol.eff.              = ',DVEFFB
      WRITE (MTERM,*) '#Steps for minimising = ',NTESTS
      WRITE (MTERM,9000) 
      
9000  FORMAT(30('-'))

      END
      