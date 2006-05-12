************************************************************************
* (Static) Grid adaption, multpliple-level version
*
* Direct version without using old SETLEV-method for switching
* levels. Works on the parameters given by the initialization
* subroutine of the multiple-grid version of the grid adaption
* (INGAST/INGAMG/INGASM).
*
* This subroutine performs grid-adaption according to information
* in the triangulation structures. The grid adaption takes place
* on the levels configured by INGASM, so the caller has not to
* take any influence concerning levels of adaption etc.
*
* By specifying another NLMAX than the standard one, the caller
* can specify the level of grid adaption, on which the parameters
* configured by INGASM (i.e. in the standard implementation, the
* parameters of the .DAT file) act.
*
* In:
*   TRIAS  - array [1..SZTRIA,1..NLMAX] of integer
*            Triangulation structures on all levels
*   NLMIN  - Minimum level in TRIAS
*   NLMAX  - Maximum level in TRIAS
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*   IDATA  - array [1..*] of integer 
*   DDATA  - array [1..*] of double 
*            User defined integer- and double-precision parameter 
*            blocks. Passed to callback routines like monitor
*            function. Not used in this routine.
************************************************************************
      
      SUBROUTINE XSMMGW (TRIAS,NLMIN,NLMAX,IGEOM,DGEOM,IDATA,DDATA)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      
      INCLUDE 'cbasicmg.inc'

      INCLUDE 'stria.inc'
      
C     parameters

      INTEGER TRIAS(SZTRIA,NNLEV),NLMIN,NLMAX,IGEOM(*),IDATA(*)
      DOUBLE PRECISION DGEOM(*),DDATA(*)

C     Variables for the static method:
      
      INCLUDE 'sgridadaptgeneralparams.inc'
      INCLUDE 'sgridadaptsolverparams.inc'
      INCLUDE 'sgridadaptstaticparams.inc'
      INCLUDE 'sgridadaptmgstaticparams.inc'
      
      INTEGER GASPRI(SZGACI+SZGLSI+SZGASI),GASPMI(SZGAMI)
      DOUBLE PRECISION GASPRD(SZGACD+SZGLSD+SZGASD)
      
C     local variables

      INTEGER II,ILEV,ILVSLV,ICONF
      DOUBLE PRECISION TGES
      DOUBLE PRECISION TGTOT,TGLSPR,TGLS,TGODE,TGGPR,TGGSM,TGMSM
      DOUBLE PRECISION TGGAPR,TGGRC
      INTEGER NLSYS ,NLSITE,NCODE ,NMFEVL
      CHARACTER CSTR*255
      
C     Initialize the parameter arrays with the variables
C     in the COMMON blocks from the DAT file.

      CALL INGAST (GASPRI, GASPRD)
      
C     Initialize the parameter block of the multiple-grid adaption
C     method

      CALL INGAMG (GASPMI)
      
C     Clear timing information for output

      TGTOT  = 0.0
      TGLSPR = 0.0
      TGLS   = 0.0
      TGODE  = 0.0
      TGGPR  = 0.0
      TGGSM  = 0.0
      TGMSM  = 0.0
      TGGAPR = 0.0
      TGGRC  = 0.0
      
      NLSYS  = 0
      NLSITE = 0
      NCODE  = 0
      NMFEVL = 0
      
C     Loop through all adaption configurations

      DO ICONF = 1,GASPMI(OIMLGAS)

C       Switch to the configuration.
C       Obtain ILEV and ILVSLV, i.e. the information on which
C       levels to adapt and to solve the LGS that occur:

        CALL INGASM (ICONF,GASPRI,GASPMI,ILEV,ILVSLV)

C       Switch to the level where we should adapt:

        II=1
        IF (ILEV.LE.0) THEN
          ILEV = NLMAX+ILEV
        END IF
        IF (ILEV.GT.NLMAX) ILEV = NLMAX
        
        IF (ILVSLV.LE.0) THEN
          ILVSLV = NLMAX+ILVSLV
        END IF
        IF (ILVSLV.GT.NLMAX) ILVSLV = NLMAX
        IF (ILVSLV.LT.NLMIN) ILVSLV = NLMIN
      
        WRITE(CSTR,'(A)') 
     *              '=============================================='
        CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,I5)') 'Grid deformation configuration: ',ICONF
        CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,'(A,I5)') 'Adaptation on level:            ',ILEV
        CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
        WRITE(CSTR,*)
        CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)

C       No grid adaption if we drop below the minimum level.

        IF (ILEV.GE.NLMIN) THEN

C         Adapt the grid and recreate missing information on all levels

          CALL XSMMGS (TRIAS,NLMIN,NLMAX,ILEV,ILVSLV,
     *                 GASPRI,GASPRD,IGEOM,DGEOM,IDATA,DDATA)
          
          TGTOT  = TGTOT  + GASPRD(OTGTOT )
          TGLSPR = TGLSPR + GASPRD(OTGLSPR)
          TGLS   = TGLS   + GASPRD(OTGLS  )
          TGODE  = TGODE  + GASPRD(OTGODE )
          TGGPR  = TGGPR  + GASPRD(OTGGPR )
          TGGSM  = TGGSM  + GASPRD(OTGGSM )
          TGMSM  = TGMSM  + GASPRD(OTGMSM )
          TGGAPR = TGGAPR + GASPRD(OTGGAPR)
          TGGRC  = TGGRC  + GASPRD(OTGGRC )
          
          NLSYS  = NLSYS  + GASPRI(ONLSYS )
          NLSITE = NLSITE + GASPRI(ONLSITE)
          NCODE  = NCODE  + GASPRI(ONCODE )
          NMFEVL = NMFEVL + GASPRI(ONMFEVL)
          
        ELSE

          WRITE(CSTR,'(A)') 'Grid adaption skipped on this level.'
          CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
          WRITE(CSTR,'(I3,A,I3)') ILEV,' = ILEV < NLMIN = ',NLMIN
          CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)

        END IF
        
      END DO

C     Finally some output, finish.

      TGES=TGTOT/100D0
      WRITE(CSTR,'(A)') 
     *            '=============================================='
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A)') 'Grid deformation statistics.'
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,*) 
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,F16.6,F14.6)') 'Absolute / Rel. time    : ',
     *             TGTOT,100D0
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,F16.6,F14.6)') 'Grid prolongation/restr.: ',
     *             TGGPR,TGGPR/TGES
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,F16.6,F14.6)') 'Preparation of grid ad. : ',
     *             TGGAPR,TGGAPR/TGES
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,F16.6,F14.6)') 'Setting up linear sys.  : ',
     *             TGLSPR,TGLSPR/TGES
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,F16.6,F14.6)') 'Solving Linear Systems  : ',
     *             TGLS,TGLS/TGES
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,F16.6,F14.6)') 'Moving grid             : ',
     *             TGODE,TGODE/TGES
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,F16.6,F14.6)') 'Grid smoothing:         : ',
     *             TGGSM,TGGSM/TGES
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,F16.6,F14.6)') 'Monitor fct. smoothing  : ',
     *             TGMSM,TGMSM/TGES
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,F16.6,F14.6)') 'Gradient calculation    : ',
     *             TGGRC,TGGRC/TGES
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,*) 
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,I12)') '#linear systems         : ',
     *             NLSYS
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,I12)') '#iterations in lin. slv.: ',
     *             NLSITE
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,I12)') '#ODE''s solved           : ',
     *             NCODE
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,'(A,I12)') '#monitor function eval. : ',
     *             NMFEVL
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE(CSTR,*) 
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)

      WRITE (CSTR,'(A)') 'Grid deformation completed.'
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') 
     *             '=============================================='
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,*)
      CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)

      END 

