************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT3D  (Release 1.1)               *
*                                                                      *
* Authors: J. Harig, S. Turek                                          *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XORS                                                                 *
*                                                                      *
* Purpose  Read subdivision from unit MFILE                            *
*          Data are assumed to be written by XOWS                      *
*                                                                      *
* Subroutines/functions called  OF0, XORA                              *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*15   Filename                                              *
* IFMT    I*4    0  Format free input                                  *
*                1  Formatted input                                    *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* MFILE   I*4    Unit number used                                      *
* CFILE   C*15   Filename used                                         *
* IER     I*4    Error indicator                                       *
*                -110  MFILE  exceeds maximum range 16,...,80          *
*                -110  Read error on unit MFILE                        *
*                -112  File CFILE could not be opened                  *
* The remaining output quantities are contained in                     *
* COMMON blocks /TRIAA/ and /TRIAD/                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE XORS(MFILE,CFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*15
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XORS'
      IF (ICHECK.GE.997) CALL OTRC('XORS  ','12/11/89')
      IER=0
C
C *** Open I/O-file
C
      CALL OF0(MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99999
C
      IF (IFMT.EQ.1) THEN
       READ (MFILE,'(3I8)',ERR=99997,END=99997)
     *                         LCORV0,LCORM0,LCORA0
       READ (MFILE,'(5I8)',ERR=99997,END=99997)
     *                         LVERT0,LEDGE0,LAREA0,LADJ0,LNPR0
       READ (MFILE,'(9I8)',ERR=99997,END=99997)
     *                         LVEL0,LEEL0,LAEL0,LVED0,LAED0,LVAR0,
     *                         LEAR0,LEVE0,LAVE0
       READ (MFILE,'(4I8)',ERR=99997,END=99997)
     *                         LBCT0,LVBD0,LEBD0,LABD0
       READ (MFILE,'(4I8)',ERR=99997,END=99997)
     *                         NEL,NVT,NET,NAT
       READ (MFILE,'(3I8)',ERR=99997,END=99997)
     *                         NVE,NEE,NAE
       READ (MFILE,'(5I8)',ERR=99997,END=99997)
     *                         NVEL,NEEL,NVED,NVAR,NEAR
       READ (MFILE,'(4I8)',ERR=99997,END=99997)
     *                         NBCT,NVBD,NEBD,NABD
      ELSE
       READ (MFILE,ERR=99997,END=99997) LCORV0,LCORM0,LCORA0
       READ (MFILE,ERR=99997,END=99997) LVERT0,LEDGE0,LAREA0,LADJ0,
     *                                  LNPR0
       READ (MFILE,ERR=99997,END=99997) LVEL0,LEEL0,LAEL0,LVED0,LAED0,
     *                                  LVAR0,LEAR0,LEVE0,LAVE0
       READ (MFILE,ERR=99997,END=99997) LBCT0,LVBD0,LEBD0,LABD0
       READ (MFILE,ERR=99997,END=99997) NEL,NVT,NET,NAT
       READ (MFILE,ERR=99997,END=99997) NVE,NEE,NAE
       READ (MFILE,ERR=99997,END=99997) NVEL,NEEL,NVED,NVAR,NEAR
       READ (MFILE,ERR=99997,END=99997) NBCT,NVBD,NEBD,NABD
      ENDIF
C
      IF (LCORV0.GT.0) THEN
       CALL XORA(LCORVG,'DCORVG',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LCORM0.GT.0) THEN
       CALL XORA(LCORMG,'DCOREG',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LCORA0.GT.0) THEN
       CALL XORA(LCORAG,'DCORAG',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LVERT0.GT.0) THEN
       CALL XORA(LVERT,'KVERT ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LEDGE0.GT.0) THEN
       CALL XORA(LEDGE,'KEDGE ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LAREA0.GT.0) THEN
       CALL XORA(LAREA,'KAREA ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LADJ0.GT.0) THEN
       CALL XORA(LADJ,'KADJ  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LNPR0.GT.0) THEN
       CALL XORA(LNPR,'KNPR  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LVEL0.GT.0) THEN
       CALL XORA(LVEL,'KVEL  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LEEL0.GT.0) THEN
       CALL XORA(LEEL,'KEEL  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LAEL0.GT.0) THEN
       CALL XORA(LAEL,'KAEL  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LVED0.GT.0) THEN
       CALL XORA(LVED,'KVED  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LAED0.GT.0) THEN
       CALL XORA(LAED,'KAED  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LVAR0.GT.0) THEN
       CALL XORA(LVAR,'KVAR  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LEAR0.GT.0) THEN
       CALL XORA(LEAR,'KEAR  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LEVE0.GT.0) THEN
       CALL XORA(LEVE,'KEVE  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LAVE0.GT.0) THEN
       CALL XORA(LAVE,'KAVE  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LBCT0.GT.0) THEN
       CALL XORA(LBCT,'KBCT  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LVBD0.GT.0) THEN
       CALL XORA(LVBD,'KVBD  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LEBD0.GT.0) THEN
       CALL XORA(LEBD,'KEBD  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
C
      IF (LABD0.GT.0) THEN
       CALL XORA(LABD,'KABD  ',MFILE,CFILE,IFMT)
       IF (IER.NE.0) GOTO 99997
      ENDIF
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORS  ')
C
99999 END
