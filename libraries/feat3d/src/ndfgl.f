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
* NDFGL                                                                *
*                                                                      *
* Purpose  Determination of the global degrees of freedom              *
*          on a given element                                          *
*                                                                      *
* Subroutines/functions called  NGLS                                   *
*                                                                      *
* Version from  07/22/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IEL      I*4    Number of element                                    *
* IPAR     I*4    Controls sorting of the resulting vectors            *
*                 (see below)                                          *
* IELTYP   I*4    Type of element corresponding to  Ennn               *
* KVERT    I*4    Array of vertices of the elements                    *
* KEDGE    I*4    Array of edges of the elements (if NET > 0)          *
* KAREA    I*4    Array of areas of the elements (IF NAT > 0)          *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KDFG     I*4    Global degrees of freedom on element IEL             *
*                  KDFG is sorted if  IPAR >= 0                        *
* KDFL     I*4    Local degrees od freedom corresponding to KDFG       *
*                  KDFL is determined only if  IPAR = 1                *
* IER      I*4    Error indicator                                      *
*                 -120 Wrong value of IELTYP                           *
*                                                                      *
************************************************************************
C
      SUBROUTINE NDFGL(IEL,IPAR,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=8,NNEE=12,NNAE=6,NNVEA=26)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION KDFG(*),KDFL(*),JVG(NNVEA),JVL(NNVEA)
      DIMENSION KVERT(NNVE,*),KEDGE(NNEE,*),KAREA(NNAE,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='NDFGL'
      IF (ICHECK.GE.999) CALL OTRC('NDFGL ','07/22/91')
C
      IER=0
C
C *** Constant elements with 1 dof per element
C
      IF (IELTYP.EQ.10) THEN
       KDFG(1)=IEL
       IF (IPAR.EQ.1) KDFL(1)=1
       GOTO 99999
      ENDIF
C
      DO 1 IVE=1,NVE
      JVG(IVE)=KVERT(IVE,IEL)
 1    JVL(IVE)=IVE
      NKE=NVE
      IF (NET.GT.0) THEN
       DO 2 IEE=1,NEE
       JVG(NKE+IEE)=KEDGE(IEE,IEL)+NVT
 2     JVL(NKE+IEE)=NKE+IEE
       NKE=NKE+NEE
      ENDIF
      IF (NAT.GT.0) THEN
       DO 3 IAE=1,NAE
       JVG(NKE+IAE)=KAREA(IAE,IEL)+NVT+NET
 3     JVL(NKE+IAE)=NKE+IAE
       NKE=NKE+NAE  
      ENDIF
      IF (IPAR.GE.0) CALL NGLS(JVG,JVL,NKE)
C
      GO TO (999,999,999,999,999,999,999,999,999,
     *       999,110,999,999,999,999,999,999,999,999,
     *       999,999,999,999,999,999,999,999,999,999,
     *       300,310,999,999,999,999,999,999,999,999), IELTYP
C
999   CONTINUE
      WRITE (CPARAM,'(I15)') IELTYP
      CALL WERR(-120,'NDFGL ')
      GOTO 99999
C
C *** Lagrange type elements with NKE dof
C
110   CONTINUE
      DO 111 IKE=1,NKE
111   KDFG(IKE)=JVG(IKE)
      IF (IPAR.EQ.1) THEN
       DO 112 IKE=1,NKE
112    KDFL(IKE)=JVL(IKE)
      ENDIF
      GOTO 99999
C
C *** Triquadratic and piecewise trilinear elements with NKE+1 dof
C
130   CONTINUE
      DO 131 IKE=1,NKE
131   KDFG(IKE)=JVG(IKE)
C *** Dof corresponding to the center of gravity
      KDFG(NKE+1)=NVT+NET+NAT+IEL
      IF (IPAR.EQ.1) THEN
       DO 132 IKE=1,NKE
132    KDFL(IKE)=JVL(IKE)
       KDFL(NKE+1)=NKE+1
      ENDIF
      GOTO 99999
C
C *** Nonconforming elements with NAE dof
C
300   CONTINUE
310   CONTINUE
      DO 311 IAE=1,NAE
      IF (NET.GT.0) THEN
       KDFG(IAE)=JVG(NKE-NAE+IAE)-NVT-NET
      ELSE
       KDFG(IAE)=JVG(NKE-NAE+IAE)-NVT
      ENDIF
311   CONTINUE
C
      IF (IPAR.EQ.1) THEN
       DO 312 IAE=1,NAE
       IF (NET.GT.0) THEN
        KDFL(IAE)=JVL(NKE-NAE+IAE)-NVE-NEE
       ELSE
        KDFL(IAE)=JVL(NKE-NAE+IAE)-NVE
       ENDIF
312    CONTINUE
      ENDIF
C
99999 END
C
C
C
      SUBROUTINE NGLS(KV1,KV2,IDIM)
C
C *** Piksort of the arrays KV1 and KV2
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KV1(*),KV2(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('NGLS  ','02/27/90')
C
      DO 1 J=2,IDIM
      IV1=KV1(J)
      IV2=KV2(J)
      DO 2 I=J-1,1,-1
      IF (KV1(I).LE.IV1) GOTO 3
      KV1(I+1)=KV1(I)
      KV2(I+1)=KV2(I)
2     CONTINUE
      I=0
3     KV1(I+1)=IV1
      KV2(I+1)=IV2
1     CONTINUE
C
      END
