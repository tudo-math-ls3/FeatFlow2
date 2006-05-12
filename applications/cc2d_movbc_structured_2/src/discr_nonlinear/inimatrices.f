************************************************************************
* This file contains initialization routines for the NSDEF2 solver.
************************************************************************

************************************************************************
* Initialize matrix structures for Navier-Stokes solver
*
* This routine allocates memory for necessary matrices and
* initializes the MATDAT structures on all levels. The matrix structures
* are initialized but the matrix entries are not calculated.
*
* In:
*   NLMIN   : Minimum level in TRIAS
*   NLMAX   : Maximum level in TRIAS
*   TRIAS   : array [1..SZTRIA,1..NNLEV] of integer
*             Triangulation structures on all levels.
*   IASMBL  : array [1..SZASMI] of integer
*   DASMBL  : array [1..SZASMD] of double
*             Integer and double prec. parameter block that controls the
*             discretization. This tells all assembly-routines how to 
*             set up the nonlinearity in the system matrices, which 
*             cubature formula to use, etc.
*
* Out:
*   MATDAT  : array [1..SZN2MI,1..NNLEV] of integer
*             TNS2DMatrixParams structure that defines the matrix
*             structure. MATDAT(.,NLMIN..NLMAX) is filled with data.
*
* Remark: If the B-matrices are to be build by exact quadrature directly
*   (IPRECB=0,1), their structure arrays are not allocated, as their
*   structure cannot be determined in advance.
************************************************************************

      SUBROUTINE ININSM (NLMIN,NLMAX,TRIAS,IASMBL,DASMBL,MATDAT)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'sassembly.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INCLUDE 'cbasicmg.inc'
      
C     parameters      
      
      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,NNLEV),IASMBL(*)
      INTEGER MATDAT(SZN2MI,NNLEV)
      DOUBLE PRECISION DASMBL(*)
      
C     local variables
      
      INTEGER I,NU,NP,IELT,NA,NB,ISYMM
      
C     definition of the elements

      EXTERNAL E030,E031,EM30,EM31,E010
      
C     Number of degrees of freedom

      INTEGER NDFGX
      EXTERNAL NDFGX

C     Get the element type for easier access.
C     Don't exploit symmetry when generating structures.

      IELT = IASMBL(OIELEMT)
      ISYMM = 0

C     Loop through all essential levels:      

      DO I=NLMIN,NLMAX

C       At first clear the structure for level I

        CALL LCL3(MATDAT(1,I),SZN2MI)

C       ---------------------------------------------------------------
C       At first generate the general matrix structure - which will
C       be shared between all mass-, Stokes- and system matrices.
C       The routine will also return in NA the number of entries and
C       in NU the number of unknowns in the system.

        IF (IELT.EQ.0) THEN
          CALL XAP7X(MATDAT(OLCLA1,I),MATDAT(OLLDA1,I),NA,
     *               NU,TRIAS(1,I),E030,ISYMM)
        ELSE IF (IELT.EQ.1) THEN
          CALL XAP7X(MATDAT(OLCLA1,I),MATDAT(OLLDA1,I),NA,
     *               NU,TRIAS(1,I),E031,ISYMM)
        ELSE IF (IELT.EQ.2) THEN
          CALL XAP7X(MATDAT(OLCLA1,I),MATDAT(OLLDA1,I),NA,
     *               NU,TRIAS(1,I),EM30,ISYMM)
        ELSE IF (IELT.EQ.3) THEN
          CALL XAP7X(MATDAT(OLCLA1,I),MATDAT(OLLDA1,I),NA,
     *               NU,TRIAS(1,I),EM31,ISYMM)
        END IF
        
        IF (IER.NE.0) RETURN

C       The above routine returns in NU the number of unknowns
C       in each velocity vector. Save it to the matrix structure.

        MATDAT (ONEQA,I) = NU
        
C       Also save the number of entries in the matrix to the structure

        MATDAT(ONA1,I) = NA
        
C       Don't forget to allocate memory for the matrix A1

        CALL ZNEW (NA,1,MATDAT(OLA1,I),'DA1   ')
        
C       ---------------------------------------------------------------
C       The 2nd A-matrix shares everything with the first one.
C       The matrices are the same.

        MATDAT(ONA2,I) = NA
        MATDAT(OLA2,I) = MATDAT(OLA1,I)
        MATDAT(OLCLA2,I) = MATDAT(OLCLA1,I)
        MATDAT(OLLDA2,I) = MATDAT(OLLDA1,I)
        
C       ---------------------------------------------------------------
C       Allocate memory for coupling matrices C1,C2.
C       *Not used in this version.*
        
C       ---------------------------------------------------------------
C       Allocate mass matrix.
C       Do we have a lumped mass matrix?

        MATDAT(OIMALMP,I) = IASMBL(OIMASS)

        IF (IASMBL(OIMASS).EQ.0) THEN

C         The mass matrix degenerates to a vector...

          CALL ZNEW (NU,1,MATDAT(OLM,I),'DMASS ')

C         The structure arrays and the number of entries
C         coincide with the system matrix
C         structure. This is a little bit inconsistent (as the
C         structure does not describe the structure of the entries),
C         but handling of lumped mass matrices always require
C         a special handling which only handles the diagonal
C         entries. The main program can switch to this special
C         handling by testing whether the mass matrix is lumped
C         or not.

          MATDAT(ONM,I)    = NA
          MATDAT(OLLDM,I)  = MATDAT(OLLDA1,I)
          MATDAT(OLCOLM,I) = MATDAT(OLCLA1,I)

C         The reason for this setting is that in the construction
C         phase when the mass matrix is lumped, first the full matrix
C         is build and then the diagonal entries of the lumped
C         matrix are generated. This always needs the structure of the
C         system matrix. To prevent the construction of the mass
C         matrix being dependent of a proper setting of the structure
C         of the system matrix, we let KCOL/KLD of the mass matrix
C         coincide with that of the system matrix. Then when generating
C         the entries, we can use the handles LCLM/LLDM rather that
C         LCLA1/LLDA1. Ok, they are the same because of the above
C         setting. But anyway, the construction of the entries of M
C         is decoupled from the structure of A, as theoretically we
C         can put different handles into LCLM/LLDM which do not
C         coincide with LCLA1/LLDA1!
          
        ELSE
        
C         Real mass matrix. We take the structure of A1.

          MATDAT(ONM,I)    = NA
          MATDAT(OLLDM,I)  = MATDAT(OLLDA1,I)
          MATDAT(OLCOLM,I) = MATDAT(OLCLA1,I)

C         Don't calculate the matrix. Allocate a vector in the
C         size of A1.

          CALL ZNEW(NA,1,MATDAT(OLM,I),'DMASS ')
          
        END IF

        IF (IER.NE.0) RETURN

C       ---------------------------------------------------------------
C       Allocate memory for Stokes matrix.
C       Only used if prescribed by IPRECA. Structure is the same as A1.

        IF (IASMBL(OIPRECA).NE.4) THEN
          CALL ZNEW(NA,1,MATDAT(OLST,I),'DSTOK ')
          MATDAT(ONST,I)   = NA
          MATDAT(OLLDST,I) = MATDAT(OLLDA1,I)
          MATDAT(OLCOLS,I) = MATDAT(OLCLA1,I)
        END IF
        
C       ---------------------------------------------------------------
C       Allocate memory for the B1 matrix.
C       Only generate the structure, don't calculate it.
C
C       Check the value of IPRECB. For IPRECB=0,1, we use exact
C       quadrature to evaluate the entry, which results in a
C       diagonal matrix - so we know the size of the matrix.
C       For IPRECB=3,4, the matrix structure is unknown and to
C       be allocated completely:

        IF (IASMBL(OIPRECB).GE.3) THEN

C         -------------------------------------------------------------
C         Generate the complete B-matrix structure for B1
C         This will generate also NB (number of unknnowns in B) and
C         NP (length of P-vector)

          IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
     *      CALL XAP9X(MATDAT(OLCLB1,I),MATDAT(OLLDB1,I),NB,NP,
     *                 TRIAS(1,I),E031,E010)
          IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
     *      CALL XAP9X(MATDAT(OLCLB1,I),MATDAT(OLLDB1,I),NB,NP,
     *                 TRIAS(1,I),E030,E010)

          IF (IER.NE.0) RETURN
          
C         Save NP to the matrix structure. We also want to save the
C         number of pressure equations - but this is not NP!
C         AP9 returns in NEQ the number of rows in the pressure matrix,
C         which is the number of velocity components. the number of
C         columns in the pressure matrix we can obtain with the
C         help of NDFG, asked for element 10 - which is the number
C         of elements. In a future version the pressure element
C         might be changed, so we hold it more or less general here.

          NP = NDFGX(10,TRIAS(1,I))

          MATDAT (ONEQB,I) = NP
          MATDAT (ONB1,I)  = NB
          
C         Don't forget to allocate memory for the B1

          CALL ZNEW (NB,1,MATDAT(OLB1,I),'DB1   ')
     
C         -------------------------------------------------------------
C         The structure of B2 is the same as B1:

          MATDAT (ONB2,I)   = NB
          MATDAT (OLCLB2,I) = MATDAT(OLCLB1,I)
          MATDAT (OLLDB2,I) = MATDAT(OLLDB1,I)
          
C         But the entries are different; allocate memory

          CALL ZNEW (NB,1,MATDAT(OLB2,I),'DB2   ')

        END IF

C       Otherwise the structure can't be determined in advance, so
C       it's determined in the calculation routine. So don't create
C       matrix structures in this case.
        
      END DO
      
      END
      
************************************************************************
* Release matrix structures for Navier-Stokes solver
*
* This routine releases memory that was allocated by ININSV.
*
* In:
*   NLMIN   : Minimum level in TRIAS
*   NLMAX   : Maximum level in TRIAS
*   MATDAT  : array [1..SZN2VI,1..NNLEV] of integer
*             TNS2DMatrixParams structure that defines the matrices.
*
* Out:
*  MATDAT(.,NLMIN..NLMAX) is released.
************************************************************************

      SUBROUTINE DONNSM (NLMIN,NLMAX,MATDAT)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'sassembly.inc'
      INCLUDE 'smat2dns.inc'
      
      INCLUDE 'cbasicmg.inc'
      
C     parameters      
      
      INTEGER NLMIN,NLMAX
      INTEGER MATDAT(SZN2MI,NNLEV)
      
C     local variables
      
      INTEGER I

C     Loop through all essential levels:      

      DO I=NLMAX,NLMIN,-1

C       Release B2-matrix; structure coincodes with B1

        MATDAT(OLLDB2,I) = 0
        MATDAT(OLCLB2,I) = 0
        IF (MATDAT(OLB2,I).NE.0)
     *    CALL ZDISP(0,MATDAT(OLB2,I),'LB2   ')

C       Release B1-matrix

        IF (MATDAT(OLLDB1,I).NE.0)
     *    CALL ZDISP(0,MATDAT(OLLDB1,I),'LLDB1 ')
        IF (MATDAT(OLCLB1,I).NE.0)
     *    CALL ZDISP(0,MATDAT(OLCLB1,I),'LCLB1 ')
        IF (MATDAT(OLB1,I).NE.0)
     *    CALL ZDISP(0,MATDAT(OLB1,I),'LB1   ')
        
C       Release Stokes matrix; structure coincodes with A1
        
        MATDAT(OLCOLS,I) = 0
        MATDAT(OLLDST,I) = 0
        IF (MATDAT(OLST,I).NE.0)
     *    CALL ZDISP(0,MATDAT(OLST,I),'LST   ')
        
C       Release mass matrix - if exists, structure coincides with A1

        MATDAT(OLLDM,I) = 0
        MATDAT(OLCOLM,I) = 0
        IF (MATDAT(OLM,I).NE.0)
     *    CALL ZDISP(0,MATDAT(OLM,I),'LM    ')

C       Release C2-matrix - this does not exist in this version!

        MATDAT(OLLDC2,I) = 0
        MATDAT(OLCLC2,I) = 0
        MATDAT(OLC2,I)   = 0

C       Release C1-matrix - this does not exist in this version!

        MATDAT(OLLDC1,I) = 0
        MATDAT(OLCLC1,I) = 0
        MATDAT(OLC1,I)   = 0
        
C       Release LA2-matrix - coincides with LA1

        MATDAT(OLLDA2,I) = 0
        MATDAT(OLCLA2,I) = 0
        MATDAT(OLA2,I)   = 0
        
C       Release A1-matrix

        IF (MATDAT(OLA1,I).NE.0)
     *    CALL ZDISP(0,MATDAT(OLA1,I),'LA1   ')
        IF (MATDAT(OLLDA1,I).NE.0)
     *    CALL ZDISP(0,MATDAT(OLLDA1,I),'LLDA1 ')
        IF (MATDAT(OLCLA1,I).NE.0)
     *    CALL ZDISP(0,MATDAT(OLCLA1,I),'LCLA1 ')
        
      END DO
      
      END

************************************************************************
* Generate matrices for Navier-Stokes solver
*
* This routine calculates the entries in the matrices depending on
* the input configuration. The bitfield IDMATS decides on which
* matrices are calculated. The memory for the matrices to be
* calculated must be allocated by the caller (e.g. with ININSM)
* prior to the call of this routine!
*
* In:
*   NLMIN   : Minimum level in TRIAS
*   NLMAX   : Maximum level in TRIAS
*   TRIAS   : array [1..SZTRIA,1..NNLEV] of integer
*             Triangulation structures on all levels.
*   IASMBL  : array [1..SZASMI] of integer
*   DASMBL  : array [1..SZASMD] of double
*             Integer and double prec. parameter block that controls the
*             discretization. This tells all assembly-routines how to 
*             set up the nonlinearity in the system matrices, which 
*             cubature formula to use, etc.
*   IDMATS  : Bitfield. Decides on which matrices are computed.
*             Bit0: Generate Stokes-(=Laplace) matrix
*             Bit1: Generate Mass matrix
*             Bit2: Generate B1- and B2-matrix
*             IDMATS=255 will build everything possible.
*
* Out:
*   MATDAT  : array [1..SZN2MI,1..NNLEV] of integer
*             TNS2DMatrixParams structure that defines the matrix
*             structure. The matrix entries referenced by 
*             MATDAT(.,NLMIN..NLMAX) are filled with data.
*
* Remarks:
*   The matrices are only build if their memory and structure 
*   was allocated/prepared previously (with ININSM).
*   Matrices that share their entries with other matrices
*   are only build once. Old matrices are overwritten.
*
*   If quadrature is activated for the B-matrices, the complete
*   B-matrix including its structure is rebuild.
************************************************************************

      SUBROUTINE GENNSM (NLMIN,NLMAX,TRIAS,IASMBL,DASMBL,IDMATS,MATDAT)
      
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'sassembly.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'arnstokes.inc'
      INCLUDE 'arnmassmat.inc'
      INCLUDE 'arnbmat.inc'
      
C     parameters      
      
      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,NNLEV),IASMBL(*),IDMATS
      INTEGER MATDAT(SZN2MI,NNLEV)
      DOUBLE PRECISION DASMBL(*)
      
C     local variables
      
      INTEGER I,NU,NP,NB,IELT,NM,ISYMM,ICLR,NST,ICUBA,ICUBM,LM
      INTEGER INU,KCOLM,KLDM,KM,KMX,ILD,ICUBB
      INTEGER LLB(2),LCOLB,LLDB
      INTEGER KVERT,KMID,KADJ,KCORVG
      DOUBLE PRECISION DMH
      
C     definition of the elements

      EXTERNAL E030,E031,EM30,EM31,E010

C local variables to define the bilinear form and the structure 
C of the matrices

C     configuration of the matrix blocks; only one block in
C     all matrices for current implementation

      INTEGER NBLOC
      PARAMETER (NBLOC=1)

C     Stokes matrix

      INTEGER KABSTN(NBLOC),KABST(2,NNAB,NBLOC)
      
      LOGICAL BCONST(NBLOC)      
      DATA BCONST/.TRUE./
      SAVE BCONST

C     Mass matrix

      INTEGER KABMN(NBLOC),KABM(2,NNAB,NBLOC)
      
      LOGICAL BCONM(NBLOC)
      DATA BCONM/.TRUE./
      SAVE BCONM
      
C     Pressure matrix; We generate both matrices at once

      INTEGER NBLOCB,NBLB1
      PARAMETER (NBLOCB=2,NBLB1=NBLOCB*NNLEV)

      INTEGER KLB(NBLB1)
      DATA KLB /NBLB1*0/
      SAVE KLB

      LOGICAL BCONB(NBLOCB)
      DATA BCONB /.TRUE.,.TRUE./
      SAVE BCONB
      
      INTEGER KABBN(NBLOCB),KABB(2,NNAB,NBLOCB)
      DATA KABBN/1,1/
      SAVE KABBN

C     External functions that define the bilinear forms

      EXTERNAL COEFST,COEFFM,COEFFB
      
C     Before we can do anything, initialize the descriptors of the
C     bilinear forms of the matrices:
C
C     Stokes matrix

      KABST (1,1,1)=2
      KABST (2,1,1)=2
      KABST (1,2,1)=3
      KABST (2,2,1)=3
      KABSTN(1)    =2
      
C     Real (not lumped) mass matrix (whether it's generated or not)

      KABM (1,1,1)=1
      KABM (2,1,1)=1
      KABMN(1)    =1

C     The two pressure matrices

      CALL LCL3(KABB,2*NNAB*NBLOCB)

C     ( A         B1 )
C     (      A    B2 )
C     ( B1^T B2^T 0  )

C     B1-block
      
      KABB(1,1,1)=2
      KABB(2,1,1)=1
      
C     B2 block

      KABB(1,1,2)=3
      KABB(2,1,2)=1
      
C     So far for the preparation. Now let's go...

C     Get the element type for easier access.
C     Don't exploit symmetry when generating structures.
C     Clear matrix before building it.

      IELT  = IASMBL(OIELEMT)
      ISYMM = 0
      ICLR  = 1
      ICUBA = IASMBL(OICUBA)
      ICUBM = IASMBL(OICUBM)
      ICUBB = IASMBL(OICUBB)

C     Loop through all essential levels:      

      DO I=NLMIN,NLMAX

C       Fetch some data for quicker access
      
        NST = MATDAT(ONST,I)
        NM  = MATDAT(ONM,I)
        NU  = MATDAT(ONEQA,I)
        NP  = MATDAT(ONEQB,I)

C       ---------------------------------------------------------------
C       If Bit0 is set and if we have a structure, calculate the 
C       Stokes matrix

        IF ((IAND(IDMATS,1).NE.0).AND.(MATDAT(OLCOLS,I).NE.0)) THEN
        
          IF (IELT.EQ.0) THEN
            CALL XAB7X(MATDAT(OLST,I),MATDAT(OLCOLS,I),
     *                  MATDAT(OLLDST,I),NST,NU,NBLOC,ICLR,TRIAS(1,I),
     *                  E031,.FALSE.,COEFST,BCONST,
     *                  KABST,KABSTN,ICUBA,ISYMM,CARRST(I),
     *                  IASMBL,DASMBL)
          ELSE IF (IELT.EQ.1) THEN
            CALL XAB7X(MATDAT(OLST,I),MATDAT(OLCOLS,I),
     *                  MATDAT(OLLDST,I),NST,NU,NBLOC,ICLR,TRIAS(1,I),
     *                  E030,.FALSE.,COEFST,BCONST,
     *                  KABST,KABSTN,ICUBA,ISYMM,CARRST(I),
     *                  IASMBL,DASMBL)
          ELSE IF (IELT.EQ.2) THEN
            CALL XAB7X(MATDAT(OLST,I),MATDAT(OLCOLS,I),
     *                  MATDAT(OLLDST,I),NST,NU,NBLOC,ICLR,TRIAS(1,I),
     *                  EM31,.TRUE.,COEFST,BCONST,
     *                  KABST,KABSTN,ICUBA,ISYMM,CARRST(I),
     *                  IASMBL,DASMBL)
          ELSE IF (IELT.EQ.3) THEN
            CALL XAB7X(MATDAT(OLST,I),MATDAT(OLCOLS,I),
     *                  MATDAT(OLLDST,I),NST,NU,NBLOC,ICLR,TRIAS(1,I),
     *                  EM30,.TRUE.,COEFST,BCONST,
     *                  KABST,KABSTN,ICUBA,ISYMM,CARRST(I),
     *                  IASMBL,DASMBL)
          END IF
        
        END IF ! AND(IDMATS,1)<>0
        
C       ---------------------------------------------------------------
C       If Bit1 is set and if we have a structure, calculate the 
C       Mass matrix

        IF ((IAND(IDMATS,2).NE.0).AND.(MATDAT(OLCOLM,I).NE.0)) THEN

C         Do we have a real mass matrix?

          IF (IASMBL(OIMASS).EQ.1) THEN
          
C           Directly build the mass matrix, overwriting the old
C           entries:

            IF (IELT.EQ.0) THEN
              CALL XAB7X(MATDAT(OLM,I),MATDAT(OLCOLM,I),
     *                    MATDAT(OLLDM,I),NM,NU,NBLOC,ICLR,
     *                    TRIAS(1,I),E031,.FALSE.,COEFFM,
     *                    BCONM,KABM,KABMN,ICUBM,ISYMM,CARRM(I),
     *                    IASMBL,DASMBL)
            ELSE IF (IELT.EQ.1) THEN
              CALL XAB7X(MATDAT(OLM,I),MATDAT(OLCOLM,I),
     *                    MATDAT(OLLDM,I),NM,NU,NBLOC,ICLR,
     *                    TRIAS(1,I),E030,.FALSE.,COEFFM,
     *                    BCONM,KABM,KABMN,ICUBM,ISYMM,CARRM(I),
     *                    IASMBL,DASMBL)
            ELSE IF (IELT.EQ.2) THEN
              CALL XAB7X(MATDAT(OLM,I),MATDAT(OLCOLM,I),
     *                    MATDAT(OLLDM,I),NM,NU,NBLOC,ICLR,
     *                    TRIAS(1,I),EM31,.TRUE.,COEFFM,
     *                    BCONM,KABM,KABMN,ICUBM,ISYMM,CARRM(I),
     *                    IASMBL,DASMBL)
            ELSE IF (IELT.EQ.3) THEN
              CALL XAB7X(MATDAT(OLM,I),MATDAT(OLCOLM,I),
     *                    MATDAT(OLLDM,I),NM,NU,NBLOC,ICLR,
     *                    TRIAS(1,I),EM30,.TRUE.,COEFFM,
     *                    BCONM,KABM,KABMN,ICUBM,ISYMM,CARRM(I),
     *                    IASMBL,DASMBL)
            END IF
          
          ELSE
          
C           We want to have a lumped mass matrix. We can either try to
C           calculate by hand or... calculate the real matrix and
C           extract cricial information. We choose the latter, as it's
C           easier :)
C           So calculate a (temporary) mass matrix like above:
          
            LM = 0
          
            IF (IELT.EQ.0) THEN
              CALL XAB7X(LM,MATDAT(OLCOLM,I),MATDAT(OLLDM,I),NM,NU,
     *                    NBLOC,ICLR,TRIAS(1,I),E031,.FALSE.,COEFFM,
     *                    BCONM,KABM,KABMN,ICUBM,ISYMM,CARRM(I),
     *                    IASMBL,DASMBL)
            ELSE IF (IELT.EQ.1) THEN
              CALL XAB7X(LM,MATDAT(OLCOLM,I),MATDAT(OLLDM,I),NM,NU,
     *                    NBLOC,ICLR,TRIAS(1,I),E030,.FALSE.,COEFFM,
     *                    BCONM,KABM,KABMN,ICUBM,ISYMM,CARRM(I),
     *                    IASMBL,DASMBL)
            ELSE IF (IELT.EQ.2) THEN
              CALL XAB7X(LM,MATDAT(OLCOLM,I),MATDAT(OLLDM,I),NM,NU,
     *                    NBLOC,ICLR,TRIAS(1,I),EM31,.TRUE.,COEFFM,
     *                    BCONM,KABM,KABMN,ICUBM,ISYMM,CARRM(I),
     *                    IASMBL,DASMBL)
            ELSE IF (IELT.EQ.3) THEN
              CALL XAB7X(LM,MATDAT(OLCOLM,I),MATDAT(OLLDM,I),NM,NU,
     *                    NBLOC,ICLR,TRIAS(1,I),EM30,.TRUE.,COEFFM,
     *                    BCONM,KABM,KABMN,ICUBM,ISYMM,CARRM(I),
     *                    IASMBL,DASMBL)
            END IF

C           Go through all lines of that matrix

            KCOLM = L(MATDAT(OLCOLM,I))
            KLDM  = L(MATDAT(OLLDM,I))
            KM    = L(LM)
            KMX   = L(MATDAT(OLM,I))

            DO INU=1,NU

C             If standard lumping is prescribed by the discretization
C             flag, simply take the diagonal entry. It's the only one in
C             the line if the correct cubature formula is active.
C       
C             Otherwise sum up the numbers on the line together and store
C             that value as the diagonal.

              IF (IASMBL(OIMASSL).EQ.1) THEN
                DMH=0D0
                DO ILD=KWORK(KLDM+INU-1),KWORK(KLDM+INU)-1
                  DMH=DMH+DWORK(KM+ILD-1)
                END DO
              ELSE
                ILD=KWORK(KLDM+INU-1)
                DMH=DWORK(KM+ILD-1)
              ENDIF

C             store the diagonal entry in the pre-allocated matrix

              DWORK(KMX+INU-1) = DMH
              
            END DO
            
C           That's it, we don't need the temporary mass matrix anymore

            CALL  ZDISP (0, LM,'DMASSH')
            IF (IER.NE.0) RETURN

          END IF ! (IASMBL(OIMASS).EQ.0)

        END IF ! AND(IDMATS,2)<>0

C       ---------------------------------------------------------------
C       If Bit2 is set and if we have a structure, calculate the 
C       B-matrices

        IF ((IAND(IDMATS,4).NE.0).AND.(MATDAT(OLCLB1,I).NE.0)) THEN
        
C         The pressure matrices can be evaluated in two ways:
C         Quick & Dirty quadrature rules or exact evaluation.
C         Which one should we take? Use IPRECB to decide...

          IF (IASMBL(OIPRECB).GE.3) THEN
          
C           We use exact evaluation. Call the evaluation routines
C           to calculate B1 and B2. Overwrite the preallocated arrays,
C           using the precalculated structure.

            KLB(1) = MATDAT(OLB1,I)
            KLB(2) = MATDAT(OLB2,I)
            NB  = MATDAT(ONB1,I)

C           Calculate the pressure part only with parametric elements -
C           it's accurate enough.

            IF ((IELT.EQ.0).OR.(IELT.EQ.2)) THEN
              CALL XAB09X(KLB,MATDAT(OLCLB1,I),MATDAT(OLLDB1,I),
     *                    NB,NBLOCB,ICLR,TRIAS(1,I),E031,E010,E010,
     *                    COEFFB,BCONB,KABB,KABBN,ICUBB,
     *                    CARRDB(I),IASMBL,DASMBL)
            ELSE IF ((IELT.EQ.1).OR.(IELT.EQ.3)) THEN
              CALL XAB09X(KLB,MATDAT(OLCLB1,I),MATDAT(OLLDB1,I),
     *                    NB,NBLOCB,ICLR,TRIAS(1,I),E030,E010,E010,
     *                    COEFFB,BCONB,KABB,KABBN,ICUBB,
     *                    CARRDB(I),IASMBL,DASMBL)
            END IF
          
C           handles should not have changed, but for sure...

            MATDAT(OLB1,I) = KLB(1)
            MATDAT(OLB2,I) = KLB(2)
          
          ELSE
          
C           Use "quick and (not so) dirty" exact quadrature rules for
C           constant pressure / linear velocities to calculate the 
C           B-matrices directly.
C
C           If there are any old matrices, release them!

            IF (MATDAT(OLLDB1,I).NE.0) 
     *        CALL ZDISP(0,MATDAT(OLLDB1,I),'KLDB  ')
     
            IF (MATDAT(OLCLB1,I).NE.0) 
     *        CALL ZDISP(0,MATDAT(OLCLB1,I),'KCOLB ')
     
            IF (MATDAT(OLB2,I).NE.0) 
     *        CALL ZDISP(0,MATDAT(OLB2,I),'DB2   ')
            IF (MATDAT(OLB1,I).NE.0) 
     *        CALL ZDISP(0,MATDAT(OLB1,I),'DB1   ')

C           Structure of the matrix B2 was the same as B1

            MATDAT(OLLDB2,I) = 0
            MATDAT(OLCLB2,I) = 0
     
C           Allocate vectors for the matrix structure. We don't assume
C           the structure to exist. We allocate more memory than we
C           actually need and will release unused memory later.
        
            NB = 2*NU
        
            CALL ZNEW(NB  ,1,LLB(1) ,'DB1   ')
            IF (IER.NE.0) RETURN
            CALL ZNEW(NB  ,1,LLB(2) ,'DB2   ')
            IF (IER.NE.0) RETURN
            CALL ZNEW(NB  ,3,LCOLB,'KCOLB ')
            IF (IER.NE.0) RETURN
            CALL ZNEW(NU+1,3,LLDB, 'LLDB  ')
            IF (IER.NE.0) RETURN
            
            KVERT  = L(TRIAS(OLVERT,I))
            KMID   = L(TRIAS(OLMID,I))
            KADJ   = L(TRIAS(OLADJ,I))
            KCORVG = L(TRIAS(OLCORVG,I))
            
            CALL BBUILD(KWORK(KVERT),KWORK(KMID),
     *             KWORK(KADJ),DWORK(KCORVG),
     *             DWORK(L(LLB(1))),DWORK(L(LLB(2))),
     *             KWORK(L(LCOLB)),KWORK(L(LLDB)),
     *             NB,TRIAS(ONEL,I),TRIAS(ONVT,I),TRIAS(ONMT,I))

C           This overwrites NB with the real number of entries
C           in the matrix.
C
C           Free unused memory; there are only NB elements in 
C           KLB, KLCOLB,...

            CALL ZDISP (NB,LLB(2) ,'DB2   ')
            CALL ZDISP (NB,LLB(1) ,'DB1   ')
            CALL ZDISP (NB,LCOLB, 'KCOLB ')
            IF (IER.NE.0) RETURN
            
C           Now we have everything together. Save the information
C           in the discretization structure. B1 and B2 share the
C           same structure.

            MATDAT(ONEQB,I) = TRIAS(ONMT,I)
            MATDAT(ONB1,I) = NB
            MATDAT(ONB2,I) = NB
            MATDAT(OLB1,I) = LLB(1)
            MATDAT(OLB2,I) = LLB(2)
            MATDAT(OLCLB1,I) = LCOLB
            MATDAT(OLCLB2,I) = LCOLB
            MATDAT(OLLDB1,I) = LLDB
            MATDAT(OLLDB1,I) = LLDB
            
          END IF ! (IASMBL(OIPRECB).GE.3)
        
        END IF ! AND(IDMATS,4)<>0
        
      END DO ! I
      
      END
      
***********************************************************************
* Implement matrix restriction, extended calling convention
*
* This routine capsules matrix restriction for the system matrix.
* The matrix entries belonging to cells with too high aspect ratio 
* are rebuild by restriction of the matrix of the next higher level.
*
* In:
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation structure
*   MATRX  : array [1..SZN2MI] of integer
*            TNS2DMatrixParams structure identifying the existing
*            matrices build of the triangulation TRIA
*   TRIAF  : array [1..SZTRIA] of integer
*            Triangulation structure of the finer level
*   MATRXF : array [1..SZN2MI] of integer
*            TNS2DMatrixParams structure identifying the existing
*            matrices build of the triangulation TRIAF
*   IAPRM  : configure adaptive P/R/matrix generation
*            =0: no adaptive matrix generation
*            =1: switch to constant mat. gen., depending on DMTEP
*            =2: switch mat. gen. depenting on size of neighbour 
*                element, too
*   DMTEP  : treshold parameter: switch construction of coarse grid 
*            matrices from standard finite element approach to locally
*            constant interpolation for all rows belonging to elements
*            with aspect ratio >= DMTEPS.
*            0D0=all elements
*            20D0=standard
*            -1D0=infinity=no matrix modifications
* Out:
*   The entries of the system matrix in MATRX are rebuild as necessary.
***********************************************************************

      SUBROUTINE IMPMRX (TRIA,MATRX,TRIAF,MATRXF,IAPRM,DMTEP)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      
C parameters

      INTEGER IAPRM
      DOUBLE PRECISION DMTEP
      INTEGER TRIA(SZTRIA),TRIAF(SZTRIA),MATRX(SZN2MI),MATRXF(SZN2MI)
      
C local variables

      INTEGER KA1,KA2,KCOLA1,KCOLA2,KLDA1,KLDA2,LLV2,LLV1,LLA2,LLA1
      INTEGER LLM2,LLM1,NVT2,NVT1,NEL2,NEL1,NMT2,NMT1,LLAR1,LLAR2
      INTEGER LCORVG
      
      INTEGER IADM1
      
C Call matrix restriction routine to rebuild the matrix
C entrys belonging to anisotropic cells. This procedure is applied to all
C matrices except for the finest level

C Initialize some helper variables
      KA1   =L(MATRX(OLA1))
      KCOLA1=L(MATRX(OLCLA1))
      KLDA1 =L(MATRX(OLLDA1))
      LLV1  =L(TRIA(OLVERT))
      LLA1  =L(TRIA(OLADJ))
      LLM1  =L(TRIA(OLMID))
      NVT1  =TRIA(ONVT)
      NEL1  =TRIA(ONEL)
      NMT1  =TRIA(ONMT)
      LLAR1 =L(TRIA(OLAREA))
      LCORVG=TRIA(OLCORVG)

      KA2   =L(MATRXF(OLA1))
      KCOLA2=L(MATRXF(OLCLA1))
      KLDA2 =L(MATRXF(OLLDA1))
      LLV2  =L(TRIAF(OLVERT))
      LLA2  =L(TRIAF(OLADJ))
      LLM2  =L(TRIAF(OLMID))
      NVT2  =TRIAF(ONVT)
      NEL2  =TRIAF(ONEL)
      NMT2  =TRIAF(ONMT)
      LLAR2 =L(TRIAF(OLAREA))

      IF ((IAPRM.NE.0).AND.(DMTEP.GE.0D0))
     *  CALL  MAREST(KWORK(LLV1),KWORK(LLV2),KWORK(LLM1),KWORK(LLM2),
     *               KWORK(LLA1),KWORK(LLA2),
     *               KWORK(KLDA1),KWORK(KLDA2),
     *               KWORK(KCOLA1),KWORK(KCOLA2),
     *               DWORK(KA1),DWORK(KA2),
     *               DWORK(L(LCORVG)),DWORK(LLAR1),
     *               NEL1,NEL2,NVT1,NVT2,
     *               DMTEP,IAPRM)
      
99999 END
      