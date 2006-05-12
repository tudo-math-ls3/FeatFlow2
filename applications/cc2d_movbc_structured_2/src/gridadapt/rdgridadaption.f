***********************************************************************
* This file contains a procedure for reading the preferences
* for the grid adaption routines from a DAT file.
***********************************************************************

***********************************************************************
* Read grid adaption DAT file
*
* Reads in the file CFNAME from hard disc into optimisation
* COMMON blocks. It will open the file CFNAME, read the parameter and
* close the file.
*
* In:
*   MDATA  - IO file handle to use.
*   MSHOW  - Level of output in the initialization phase
*   CFNAME - name of the file
***********************************************************************

      SUBROUTINE RDGSMT (MDATA,MSHOW,CFNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'

      INCLUDE 'sgridadaptmgstaticparams.inc'      
      INCLUDE 'cgridadapt.inc'
      
C parameters

      CHARACTER CFNAME*(*)
      INTEGER MDATA,MSHOW
      
C local variables

      CHARACTER CTMP*(256)
      INTEGER IVERS,IVERSF,IFMTS
      
      WRITE (CTMP,*) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CTMP)
      WRITE (CTMP,'(A)') ' Grid adaption parameters'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CTMP)
      WRITE (CTMP,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CTMP)

C Open DAT file

      IFMTS = 1
      CLOSE (MDATA)
      CALL  OF0 (MDATA,CFNAME,IFMTS)

C Check version number of DAT file.
C Current program version number:

      IVERS = 112

      CALL GETINT (MDATA,MSHOW,
     *            'Version number of DAT-file (grid sm.) '//
     *            'IVRSGS = ',IVERSF)
      IF (IVERSF.NE.IVERS) THEN
        WRITE (*,'(A)') 'Fatal error! Version number of DAT-file '//
     *                  'incorrect!'
        WRITE (*,'(A,I6)') 'Version number expected: ',IVERS
        WRITE (*,'(A)') 'Program halted because of incompatibilities!'
        CLOSE (MDATA)
        CALL EXIT (255)
      END IF

C Read all parameters, take care of separation lines

      WRITE (CTMP,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CTMP)

      READ (MDATA,'(A)') CTMP
      CALL GETINT (MDATA,MSHOW,
     *              'Method of grid deformation:           '//
     *              'IGAMET = ',IGAMET)

      CALL GETINT (MDATA,MSHOW,
     *              '#grid adaption configurations:        '//
     *              'IMLGAS = ',IMLGAS)

      CALL GETINM (MDATA,MSHOW,
     *              'Levels of grid adaption:              '//
     *              'IGALVL = ',IGALVL,IMLGAS)

      CALL GETINM (MDATA,MSHOW,
     *              'Levels of LGS for grid adaption       '//
     *              'IGASLV = ',IGASLV,IMLGAS)

      CALL GETINM (MDATA,MSHOW,
     *              'Number of grid adaption steps:        '//
     *              'NGASTP = ',NGASTP,IMLGAS)

      CALL GETINM (MDATA,MSHOW,
     *              'Number of grid postcorrection steps:  '//
     *              'NGAPCS = ',NGAPCS,IMLGAS)

      CALL GETINT (MDATA,MSHOW,
     *              'Adaptive restart control:             '//
     *              'IARCTR = ',IARCTR)

      CALL GETINT (MDATA,MSHOW,
     *              'Type of monitor function:             '//
     *              'IMNREL = ',IMNREL)

      CALL GETDBL (MDATA,MSHOW,
     *              'Monitor function tolerance:           '//
     *              'DMNTOL = ',DMNTOL)

      CALL GETDBL (MDATA,MSHOW,
     *              'Mon. fct. restriction/scaling factor: '//
     *              'DMNSFC = ',DMNSFC)

      CALL GETDBL (MDATA,MSHOW,
     *              'EPS0 parameter for monitor function:  '//
     *              'EPS0GS = ',EPS0GS)

      WRITE (CTMP,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CTMP)

      READ (MDATA,'(A)') CTMP

      CALL GETINT (MDATA,MSHOW,
     *              'ODE solver:                           '//
     *              'ODESLV = ',ODESLV)

      CALL GETINT (MDATA,MSHOW,
     *              'Number of timesteps:                  '//
     *              'NITDT  = ',NITDT)

      CALL GETINT (MDATA,MSHOW,
     *              'Grid smoothing operator:              '//
     *              'ISMMTH = ',ISMMTH)

      CALL GETINM (MDATA,MSHOW,
     *              'Number of grid pre-smoothing steps:   '//
     *              'NPRGSM = ',NPRGSM,IMLGAS)

      CALL GETINM (MDATA,MSHOW,
     *              'Number of grid post-smoothing steps:  '//
     *              'NPOGSM = ',NPOGSM,IMLGAS)

      CALL GETDBL (MDATA,MSHOW,
     *              'Adaptive post-smoothing factor:       '//
     *              'DAPOSM = ',DAPOSM)

      CALL GETINM (MDATA,MSHOW,
     *              'Number of mof. fct. smoothing steps:  '//
     *              'NFSMTH = ',NFSMTH,IMLGAS)

      CALL GETDBL (MDATA,MSHOW,
     *              'Adaptive mon. fct. smoothing factor:  '//
     *              'DFAPSM = ',DFAPSM)

      CALL GETDBL (MDATA,MSHOW,
     *              'Mon. fct. smoothing reduction factor: '//
     *              'DFPSRF = ',DFPSRF)

      WRITE (CTMP,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CTMP)

      READ (MDATA,'(A)') CTMP

      CALL GETINT (MDATA,MSHOW,
     *              'Write computed grids to files:        '//
     *              'IGSOUT = ',IGSOUT)

      CALL GETINT (MDATA,MSHOW,
     *              'Colorize macros in GMV-output:        '//
     *              'ICOLMC = ',ICOLMC)

      WRITE (CTMP,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CTMP)

      READ (MDATA,'(A)') CTMP

      CALL GETINT (MDATA,MSHOW,
     *              'Cubature formula for grid sm. mat.:   '//
     *              'ICUBGS = ',ICUBGS)

      CALL GETINT (MDATA,MSHOW,
     *              'Cubature formula for grid sm. RHS :   '//
     *              'ICURGS = ',ICURGS)

      CALL GETINT (MDATA,MSHOW,
     *              'Cubature formula for normalizing  :   '//
     *              'ICUNGS = ',ICUNGS)

      CALL GETINT (MDATA,MSHOW,
     *              'Method for reconstructing gradients:  '//
     *              'IRCGRD = ',IRCGRD)
     
      WRITE (CTMP,9001) 
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CTMP)

      READ (MDATA,'(A)') CTMP
     
      CALL GETINT (MDATA,MSHOW,
     *              'Type of solver used for grid smooth.: '//
     *              'IGSSLV = ',IGSSLV)

      CALL GETDBL (MDATA,MSHOW,
     *              'Accuracy of grid smoothing solver:    '//
     *              'EPSGSS = ',EPSGSS)

      CALL GETINT (MDATA,MSHOW,
     *              'Relative/absolute error treatment:    '//
     *              'IRELGS = ',IRELGS)

      CALL GETINT (MDATA,MSHOW,
     *              'Max. number of solver iterations:     '//
     *              'NITSGS = ',NITSGS)

      CALL GETINT (MDATA,MSHOW,
     *              'Type of smoother/preconditioner:      '//
     *              'IPCSMG = ',IPCSMG)

      CALL GETDBL (MDATA,MSHOW,
     *              'Damping parameter of smoother/prec.:  '//
     *              'OMGPGS = ',OMGPGS)

      CALL GETINT (MDATA,MSHOW,
     *              'Number of pre-/postsmoothing steps:   '//
     *              'NSMSGS = ',NSMSGS)
     
      CALL GETINT (MDATA,MSHOW,
     *              'Preconditioning in MG-smoother:       '//
     *              'IPCSGS = ',IPCSGS)
     
      CALL GETINT (MDATA,MSHOW,
     *              'Type of resorting strategy:           '//
     *              'ISRTGS = ',ISRTGS)

      CALL GETINT (MDATA,MSHOW,
     *              'Type of coarse grid solver in MG:     '//
     *              'ICGSGS = ',ICGSGS)

      CALL GETINT (MDATA,MSHOW,
     *              'Level of the coarse grid slv. in MG:  '//
     *              'ICGSLV = ',ICGSLV)

      CALL GETINT (MDATA,MSHOW,
     *              'NIT of coarse grid solver in MG:      '//
     *              'NICGGS = ',NICGGS)

      CALL GETDBL (MDATA,MSHOW,
     *              'Accuracy of coarse grid solver in MG: '//
     *              'EPSCGS = ',EPSCGS)

      CALL GETINT (MDATA,MSHOW,
     *              'Relative/absolute error for C.G.S.:   '//
     *              'IRCGGS = ',IRCGGS)

      CALL GETINT (MDATA,MSHOW,
     *              'Preconditioner in MG-coarse gris s.:  '//
     *              'IPCGGS = ',IPCGGS)

      CALL GETDBL (MDATA,MSHOW,
     *              'Damping parameter in C.G.-solver:     '//
     *              'OMCGGS = ',OMCGGS)

      CLOSE (MDATA)

9000  FORMAT(79('-'))
9001  FORMAT(60('-'))

      END
      
***********************************************************************
* Transfer grid adaption parameter to data blocks of static method
*
* This routine will build, based on the parameters in the COMMON
* blocks, a TIGridAdaptionStaticParams and TDGridAdaptionStaticParams-
* parameter structure that can be used to call the static grid
* adaption method.
*
* In:
*   Variables in the grid-adaption COMMON blocks that save the
*   state of the INI-file of the grid adaption
*
* Out:
*   IPARAM - array [1..SZGACI+SZGLSI+SZGASI] of integer
*              Index [ 1 .. SZGACI ]
*                    TIGridAdaptionCommon-structure
*              Index [ SZGACI+1 .. SZGACI+SZGLSI ]
*                    TIGridAdaptionSolver-structure
*              Index [ SZGACI+SZGLSI+1 .. SZGACI+SZGLSI+SZGASI ]
*                    TIGridAdaptionStaticParams-structure
*            This is the integer-parameter block for the static
*            grid adaption routine. The input-variables of these
*            three blocks are filled with data.
*   DPARAM - array [1..*] of double precision
*              Index [ 1 .. SZGACD ]
*                    TDGridAdaptionCommon-structure
*              Index [ SZGACD+1 .. SZGACD+SZGLSD ]
*                    TDGridAdaptionSolver-structure
*              Index [ SZGACD+SZGLSD+1 .. SZGACD+SZGLSD+SZGASD ]
*                    TDGridAdaptionStaticParams-structure
*            This is the double precision-parameter block for the 
*            static grid adaption routine. The input-variables of these
*            three blocks are filled with data.
*   IPGAMG - array [1..SZGAMI] of integer
*            = TIGridAdaptionMGStaticParams-structure
*
*   These both structures contain the integer and double information
*   that define the behaviour of the static grid adaption method.
*   The structures can be modifies by the user before calling the
*   grid adaption if necessary.
***********************************************************************

      SUBROUTINE INGAST (IPARAM, DPARAM)
      
      IMPLICIT NONE
      
C     Include the three include-files that define the three structures:
      
      INCLUDE 'sgridadaptgeneralparams.inc'
      INCLUDE 'sgridadaptsolverparams.inc'
      INCLUDE 'sgridadaptstaticparams.inc'
      
      INCLUDE 'sgridadaptmgstaticparams.inc'

      INCLUDE 'cgridadapt.inc'
      
      INTEGER IPARAM (SZGACI+SZGLSI+SZGASI)
      DOUBLE PRECISION DPARAM (SZGACD+SZGLSD+SZGASD)
      
      INTEGER START
      
C     Clear the parameter blocks before doing anything else:

      CALL LCL1 (DPARAM,SZGACD+SZGLSD+SZGASD)
      CALL LCL3 (IPARAM,SZGACI+SZGLSI+SZGASI)
      
C     Transfer general integer-parameters into the first block

      IPARAM (ORESULT) = 0
      IPARAM (OIARCTR) = IARCTR
      IPARAM (OIGSOUT) = IGSOUT
      IPARAM (OICOLMC) = ICOLMC
      IPARAM (OISMMTH) = ISMMTH
      IPARAM (OIMNREL) = IMNREL
      IPARAM (OICUNGS) = ICUNGS

C     Transfer general double-parameters into the first block

      DPARAM (ODMNTOL) = DMNTOL
      DPARAM (OEPS0GS) = EPS0GS
      DPARAM (ODAPOSM) = DAPOSM
      DPARAM (ODMNSFC) = DMNSFC

C     Transfer solver integer-parameters into the second block

      START = SZGACI+1-1
      
      IPARAM (START+OIGSSLV) = IGSSLV
      IPARAM (START+OICUBGS) = ICUBGS
      IPARAM (START+OICURGS) = ICURGS
      IPARAM (START+OIRCGRD) = IRCGRD
      IPARAM (START+OIRELGS) = IRELGS
      IPARAM (START+ONITSGS) = NITSGS
      IPARAM (START+OIPCSMG) = IPCSMG
      IPARAM (START+ONSMSGS) = NSMSGS
      IPARAM (START+OIPCSGS) = IPCSGS
      IPARAM (START+OISRTGS) = ISRTGS
      IPARAM (START+OICGSGS) = ICGSGS
      IPARAM (START+OICGSLV) = ICGSLV
      IPARAM (START+ONICGGS) = NICGGS
      IPARAM (START+OIRCGGS) = IRCGGS
      IPARAM (START+OIPCGGS) = IPCGGS

C     Transfer solver double parameters into the second block
      
      START = SZGACD+1-1
      
      DPARAM (START+OEPSGSS) = EPSGSS
      DPARAM (START+OOMGPGS) = OMGPGS
      DPARAM (START+OEPSCGS) = EPSCGS
      DPARAM (START+OOMCGGS) = OMCGGS
      
C     Transfer static method parameters into the third block

      START = SZGACI+SZGLSI+1-1

      IPARAM (START+OODESLV) = ODESLV
      IPARAM (START+ONITDT)  = NITDT
      
C     Activate a standard configuration for the adaption:

      IPARAM(ONGASTP) = 3 
      IPARAM(START+ONGAPCS) = 2 
      IPARAM(ONPRGSM) = 20
      IPARAM(ONPOGSM) = 0
      IPARAM(ONFSMTH) = 0
      
      END

***********************************************************************
* Initialize parameter block of multiple-grid grid adaption
*
* This routine will build, based on the parameters in the COMMON
* blocks, a TIGridAdaptionMGStaticParams parameter structure that
* can be used to call in the multiple-grid version of the static grid
* adaption method.
*
* In:
*   Variables in the grid-adaption COMMON blocks that save the
*   state of the INI-file of the grid adaption
*
* Out:
*   IPGAMG - array [1..SZGAMI] of integer
*            = TIGridAdaptionMGStaticParams-structure
*
*   These both structures contain the integer and double information
*   that define the behaviour of the static grid adaption method.
*   The structures can be modifies by the user before calling the
*   grid adaption if necessary.
***********************************************************************

      SUBROUTINE INGAMG (IPGAMG)
      
      IMPLICIT NONE
      
C     Include the three include-files that define the three structures:
      
      INCLUDE 'sgridadaptmgstaticparams.inc'

      INCLUDE 'cgridadapt.inc'
      
      INTEGER IPGAMG (SZGAMI)
      
      INTEGER I
      
C     Clear the parameter block before doing anything else:

      CALL LCL3 (IPGAMG,SZGAMI)
      
C     Transfer the "multiple-grid" adaption parameters to the
C     parameter block

      IPGAMG(OIMLGAS) = IMLGAS
      
      DO I=1,IMLGAS
        IPGAMG(OIGALVL +I-1) = IGALVL(I)
        IPGAMG(OIGASLV +I-1) = IGASLV(I)
        IPGAMG(OMNGASTP+I-1) = NGASTP(I)
        IPGAMG(OMNGAPCS+I-1) = NGAPCS(I)
        IPGAMG(OMNPRGSM+I-1) = NPRGSM(I)
        IPGAMG(OMNPOGSM+I-1) = NPOGSM(I)
        IPGAMG(OMNFSMTH+I-1) = NFSMTH(I)
      END DO
      
      END

***********************************************************************
* Switch configuration in multiple grid adaption
*
* This routine switches the configuration of the "multiple-grid"
* adaption algorithm. Up to NGAMXL configurations are stored in the
* 4th parameter block IPGAMG. The routine changes parameters in
* the parameter blocks of the static grid adaption according
* to the INUM'th configuration in IPGAMG and returns the information
* on which level to solve the grid adaption LGS where to adapt.
*
* In:
*   INUM   - Number of the configuration to activate
*   IPARAM - array [1..SZGACI+SZGLSI+SZGASI] of integer
*              Index [ 1 .. SZGACI ]
*                    TIGridAdaptionCommon-structure
*              Index [ SZGACI+1 .. SZGACI+SZGLSI ]
*                    TIGridAdaptionSolver-structure
*              Index [ SZGACI+SZGLSI+1 .. SZGACI+SZGLSI+SZGASI ]
*                    TIGridAdaptionStaticParams-structure
*            Must have been initialized e.g. by INGAST.
*   IPGAMG - array [1..SZGAMI] of integer
*            = TIGridAdaptionMGStaticParams-structure
*            Must have been initialized e.g. with INGAMG.
*
* Out:
*   IPARAM - array [1..SZGACI+SZGLSI+SZGASI] of integer
*              Index [ 1 .. SZGACI ]
*                    TIGridAdaptionCommon-structure
*              Index [ SZGACI+1 .. SZGACI+SZGLSI ]
*                    TIGridAdaptionSolver-structure
*              Index [ SZGACI+SZGLSI+1 .. SZGACI+SZGLSI+SZGASI ]
*                    TIGridAdaptionStaticParams-structure
*              Index [ SZGACI+SZGLSI+SZGASI+1 .. SZGACI+SZGLSI+SZGASI+SZGAMI ]
*                    TIGridAdaptionMGStaticParams-structure
*             The parameters are changed to the INUM'th configuration
*             of the IPGAMG structure.
*   ILEV   - Level where to adapt
*   ILVSLV - Level where to solve the linear system
***********************************************************************

      SUBROUTINE INGASM (INUM,IPARAM,IPGAMG,ILEV,ILVSLV)
      
      IMPLICIT NONE
      
C     Include the three include-files that define the three structures:
      
      INCLUDE 'sgridadaptgeneralparams.inc'
      INCLUDE 'sgridadaptsolverparams.inc'
      INCLUDE 'sgridadaptstaticparams.inc'
      
      INCLUDE 'sgridadaptmgstaticparams.inc'

      INTEGER IPARAM (SZGACI+SZGLSI+SZGASI),ILEV,ILVSLV,INUM
      INTEGER IPGAMG (SZGAMI)
      
      INTEGER START,N
      
C     Transfer the desired configuration from the 4th block to the
C     first three (actually the first two) blocks:

      START = SZGACI+SZGLSI+1-1
      
      N=MIN(MAX(MAX(INUM,IPARAM(OIMLGAS)),1),NGAMXL)
      
      IPARAM(ONGASTP) = IPGAMG(OMNGASTP+N-1) 
      IPARAM(START+ONGAPCS) = IPGAMG(OMNGAPCS+N-1) 
      IPARAM(ONPRGSM) = IPGAMG(OMNPRGSM+N-1) 
      IPARAM(ONPOGSM) = IPGAMG(OMNPOGSM+N-1) 
      IPARAM(ONFSMTH) = IPGAMG(OMNFSMTH+N-1)
      
      ILEV   = IPGAMG(OIGALVL +N-1)
      ILVSLV = IPGAMG(OIGASLV +N-1)
      
      END
