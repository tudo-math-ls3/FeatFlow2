************************************************************************
* Read parameters for postprocessing
*
* This routine reads the parameters for the standard postprocessing
* and stores them in the corresponding COMMON blockdiscretization.
*
* In:
*   MDATA  - Unit number to use for reading process
*   MSHOW  - Level of output in the initialization phase
*   CFNAME - Name of the file
************************************************************************

      SUBROUTINE RDPOST (MDATA,MSHOW,CFNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cpostproc.inc'
      
C     parameters
      
      CHARACTER CFNAME*(*)
      INTEGER MSHOW,MDATA
      
C     local variables

      INTEGER I,IFMTS
      CHARACTER CSTR*(255)

C     Open DAT file
      
      IFMTS = 1
      CALL OF0 (MDATA,CFNAME,IFMTS)
      
C     Ignore the first 6 lines

      DO I=1,6
        READ (MDATA,*)
      END DO
      
C     Check version number

      READ (MDATA,*) I
      IF (I.NE.100) THEN
        WRITE (*,'(A)') 
     *        'Error in reading time-discretization parameters'
        WRITE (*,'(A)') 'Version number incorrect'
        STOP
      END IF
      
      WRITE (CSTR,'(A)') ' Postprocessing parameters'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') '---------------------------'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Ignore separator
      
      READ (MDATA,*)
      
C     Read the names of the output files      
      
      READ (MDATA,*) CFLAUT
      WRITE (CSTR,1003) 
     *      'Filename auto-save            : CFLAUT = ', CFLAUT
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLSAV
      WRITE (CSTR,1003) 
     *      'Filename outp., step-b., sol. : CFLSAV = ', CFLSAV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLDUP
      WRITE (CSTR,1003) 
     *      'Filename outp., time-b., sol. : CFLDUP = ', CFLDUP
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLDU
      WRITE (CSTR,1003) 
     *      'Filename film, X-velocity     : CFLDU  = ', CFLDU
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLDV
      WRITE (CSTR,1003) 
     *      'Filename film, Y-velocity     : CFLDV  = ', CFLDV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLDP
      WRITE (CSTR,1003) 
     *      'Filename film, pressure       : CFLDP  = ', CFLDP
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLISO
      WRITE (CSTR,1003) 
     *      'Filename film, streamfct.     : CFLISO = ', CFLISO
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLGMV
      WRITE (CSTR,1003) 
     *      'Filename, GMV-output          : CFLGMV = ', CFLGMV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLFX1
      WRITE (CSTR,1003) 
     *      'Filename track., bd.i., X-dir.: CFLFX1 = ', CFLFX1
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLFY1
      WRITE (CSTR,1003) 
     *      'Filename track., bd.i., Y-dir.: CFLFY1 = ', CFLFY1
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLFX2
      WRITE (CSTR,1003) 
     *      'Filename track., vl.i., X-dir.: CFLFX2 = ', CFLFX2
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLFY2
      WRITE (CSTR,1003) 
     *      'Filename track., vl.i., Y-dir.: CFLFY2 = ', CFLFY2
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLFX3
      WRITE (CSTR,1003) 
     *      'Filename track., li.i., X-dir.: CFLFX3 = ', CFLFX3
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLFY3
      WRITE (CSTR,1003) 
     *      'Filename track., li.i., Y-dir.: CFLFY3 = ', CFLFY3
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLVLX
      WRITE (CSTR,1003) 
     *      'Filename track., X-vel. in pt.: CFLVLX = ', CFLVLX
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLVLY
      WRITE (CSTR,1003) 
     *      'Filename track., Y-vel. in pt.: CFLVLY = ', CFLVLY
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLPRE
      WRITE (CSTR,1003) 
     *      'Filename track., press. in pt.: CFLPRE = ', CFLPRE
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLSTR
      WRITE (CSTR,1003) 
     *      'Filename track., str.f. in pt.: CFLSTR = ', CFLSTR
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CFLIPR
      WRITE (CSTR,1003) 
     *      'Filename track., integr.press.: CFLIPR = ', CFLIPR
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Error analysis      
      
      READ (MDATA,*) IERANA
      IERANA=ABS(IERANA)
      WRITE (CSTR,1000) 
     *      'Error analysis                : IERANA = ', IERANA
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Auto-save during nonlinear iteration each IAUSAV iterations

      READ (MDATA,*) IAUSAV
      IAUSAV=ABS(IAUSAV)
      WRITE (CSTR,1000) 
     *      'MOD-factor for autosave in NLI: IAUSAV = ', IAUSAV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
            
C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Saving of intermediate solution vectors during time-dependent
C     simulation

      READ (MDATA,*) ISAV
      WRITE (CSTR,1000) 
     *      'Write nonsteady solutions:      ISAV   = ', ISAV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) INSAV
      WRITE (CSTR,1000) 
     *      'Stepsize for nonsteady savings: INSAV  = ', INSAV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) INSAVN
      INSAVN=ABS(INSAVN)
      WRITE (CSTR,1000) 
     *      'Modulo-no. for nonst. savings : INSAVN = ', INSAVN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Parameters for film output      

      READ (MDATA,*) DTFILM      
      WRITE (CSTR,1001) 
     *      'Time step for film output     : DTFILM = ', DTFILM
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IUPSAV
      WRITE (CSTR,1000) 
     *      'Level for solution film output: IUPSAV = ', IUPSAV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IFUSAV
      WRITE (CSTR,1000) 
     *      'Level for velocity film output: IFUSAV = ', IFUSAV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IFPSAV
      WRITE (CSTR,1000) 
     *      'Level for pressure film output: IFPSAV = ', IFPSAV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IFXSAV
      WRITE (CSTR,1000) 
     *      'Level for streamf. film output: IFXSAV = ', IFXSAV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IFINIT
      WRITE (CSTR,1000) 
     *      'Index of first file, film outp: IFINIT = ', IFINIT
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)

C     Parameters for GMV output      

      READ (MDATA,*) DTGMV      
      WRITE (CSTR,1001) 
     *      'Time step for GMV output      : DTGMV  = ', DTGMV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IGMV
      WRITE (CSTR,1000) 
     *      'Level for GMVoutput           : IGMV   = ', IGMV
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IGINIT
      WRITE (CSTR,1000) 
     *      'Index of first file, GMV outp.: IGINIT = ', IGINIT
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) ITRCES
      WRITE (CSTR,1000) 
     *      'Tracers per row/col/element   : ITRCES = ', ITRCES
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     That's it.
      
      CLOSE (MDATA)
      
1000  FORMAT (A,I4)
1001  FORMAT (A,E16.8)
1002  FORMAT (A,2I4)
1003  FORMAT (A,A)

9000  FORMAT(79('-'))
9001  FORMAT(60('-'))

      END
      