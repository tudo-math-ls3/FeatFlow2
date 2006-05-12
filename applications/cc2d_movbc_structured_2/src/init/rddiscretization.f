************************************************************************
* Read discretization parameters
*
* This routine reads the parameters about the discretization from
* a given file and stores them in the COMMON block for the
* discretization.
*
* In:
*   MDATA  - Unit number to use for reading process
*   MSHOW  - Level of output in the initialization phase
*   CFNAME - Name of the file
************************************************************************

      SUBROUTINE RDDISC (MDATA,MSHOW,CFNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cdiscretization.inc'
      
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
        WRITE (*,'(A)') 'Error in reading discretization parameters'
        WRITE (*,'(A)') 'Version number incorrect'
        STOP
      END IF
      
      WRITE (CSTR,'(A)') ' Discretization parameters'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') '---------------------------'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Ignore separator
      
      READ (MDATA,*)
      
C     Read information about start-input and final-output file

      READ (MDATA,*) ISTART
      WRITE (CSTR,1000) 
     *      'Read start vector             : ISTART = ', ISTART
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) CSTART
      
      WRITE (CSTR,1003) 
     *      'Filename of start vector      : CSTART = ', CSTART
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) ISOL
      
      WRITE (CSTR,1000) 
     *      'Write sol. vector             : ISOL   = ', ISOL
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) CSOL
      
      WRITE (CSTR,1003) 
     *      'Filename of sol. vector       : CSOL   = ', CSOL
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Separator
        
      READ (MDATA,*)
      
      WRITE (MTERM,9001)
      
C     Minimum/Maximum level of the computation

      READ (MDATA,*) NLMIN
      READ (MDATA,*) NLMAX
      
C     Maximum level must be greater than minimum level.
C     Minimum level may be negative; then the program must treat
C     the minimum level as (NLMAX+NLMIN)!

      NLMAX = MAX(1,ABS(NLMAX))
      IF (NLMIN.GT.NLMAX) NLMIN = NLMAX
      
      WRITE (CSTR,1000) 
     *      'Minimum mg-level              : NLMIN  = ', NLMIN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      WRITE (CSTR,1000) 
     *      'Minimum mg-level              : NLMIN  = ', NLMIN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      
      WRITE (MTERM,9001)

C     Configuration of the equation

      READ (MDATA,*) ISTOK
      IF (ISTOK.NE.1) ISTOK=0
      
      WRITE (CSTR,1000) 
     *      'Stokes calculation            : ISTOK  = ', ISTOK
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) RE
      
      RE = ABS(RE)
      
      WRITE (CSTR,1001) 
     *      'Viscosity parameter 1/NU      : RE     = ', RE
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) IRHS
      
      WRITE (CSTR,1000) 
     *      'RHS generation                : IRHS   = ',IRHS
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IBDR
      
      WRITE (CSTR,1000) 
     *      'Boundary generation           : IBDR   = ',IBDR
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Separator
        
      READ (MDATA,*)
      
      WRITE (MTERM,9001)

C     Element type
      
      READ(MDATA,*) IELT

      IF ((IELT.LT.0).OR.(IELT.GT.3)) IELT=3

      WRITE (CSTR,1000) 
     *      'Element type                  : IELT   = ',IELT
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      
      WRITE (MTERM,9001)
      
C     Parameters about the stabilization of the convective part
      
      READ (MDATA,*) IUPW

      WRITE (CSTR,1000) 
     *      'Convective part               : IUPW   = ',IUPW
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) UPSAM
      
      WRITE (CSTR,1001) 
     *      'Parameter for Samarskij-upwind: UPSAM  = ', UPSAM
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)
      
C     Parameters about the inflow
      
      READ (MDATA,*) IINFCF

      WRITE (CSTR,1000) 
     *      'Type of inflow profile        : IINFCF = ',IINFCF
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ (MDATA,*) DPUMAX
      
      WRITE (CSTR,1001) 
     *      'Inflow velocity               : DPUMAX = ', DPUMAX
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)
      
C     Discretization of mass matrix

      READ (MDATA,*) IMASS
      
      IF (IMASS.NE.1) IMASS=0
      
      WRITE( CSTR,1000) 
     *      'Mass evaluation               : IMASS  = ',IMASS
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IMASSL
      
      IF (IMASSL.NE.1) IMASSL=0
      
      WRITE (CSTR,1000) 
     *      'Lumped mass evaluation        : IMASSL = ',IMASSL
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
 
C     Discretization of Stokes/pressure matricex
      
      READ (MDATA,*) IPRECA
      
      IF ((IPRECA.LT.0).OR.(IPRECA.GT.4)) IPRECA=1
      IF ((IPRECA.EQ.4).AND.(IUPW.EQ.1))  IPRECA=1
      IF ((IPRECA.EQ.4).AND.(ISTOK.EQ.1)) IPRECA=1
      IF (IPRECA.EQ.3)                    IPRECA=2
      
      WRITE (CSTR,1000) 
     *      'Accuracy for ST               : IPRECA = ',IPRECA
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) IPRECB
      
      IF ((IPRECB.LT.0).OR.(IPRECB.GT.4)) IPRECB=2
      IF (IPRECB.EQ.1) IPRECB=0
      IF (IPRECB.EQ.2) IPRECB=3
      IF (IPRECB.EQ.4) IPRECB=3
      
      WRITE (CSTR,1000) 
     *      'Accuracy for B                : IPRECB = ',IPRECB
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
C     Cubature formulas for the matrices

      READ (MDATA,*) ICUBM
      ICUBM=ABS(ICUBM)
      WRITE (CSTR,1000) 
     *      'ICUB mass matrix              : ICUBM  = ',ICUBM
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) ICUBA
      ICUBA=ABS(ICUBA)
      WRITE (CSTR,1000) 
     *      'ICUB diff. matrix             : ICUBA  = ',ICUBA
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) ICUBN
      ICUBN=ABS(ICUBN)
      WRITE (CSTR,1000) 
     *      'ICUB conv. matrix             : ICUBN  = ',ICUBN
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) ICUBB
      ICUBB=ABS(ICUBB)
      WRITE (CSTR,1000) 
     *      'ICUB matrices B1,B2           : ICUBB  = ',ICUBB
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) ICUBF
      ICUBF=ABS(ICUBF)
      WRITE (CSTR,1000) 
     *      'ICUB right hand side          : ICUBF  = ',ICUBF
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Separator
        
      READ (MDATA,*)
      WRITE (MTERM,9001)
      
C     Adaptive matrix generation

      READ (MDATA,*) IAPRM      
      IF ((IAPRM.LT.0)) IAPRM=3
      WRITE (CSTR,1000) 
     *      'Type of adaptive matrix-gener.: IAPRM  = ', IAPRM
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      READ (MDATA,*) DMTEP      
      WRITE (CSTR,1001) 
     *      'Threshold.par. ad. mat.-gen.  : DMTEP  = ', DMTEP
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     That's it.
      
      WRITE (MTERM,*)
      
      CLOSE (MDATA)
      
1000  FORMAT (A,I4)
1001  FORMAT (A,E16.8)
1002  FORMAT (A,2I4)
1003  FORMAT (A,A)

9000  FORMAT(79('-'))
9001  FORMAT(60('-'))

      END
      