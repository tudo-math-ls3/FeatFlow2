************************************************************************
* This file defines different general postprocessing routines
* for a solution vector. The routines are basic routines of various
* types. The following routines can be found here:
*
*  PPRDVC   - read a solution vector from disc
*  PPWRVC   - write a solution vector to disc
*  INTUV    - Interpolate midpoint-values to values in vertices of
*             a triangulation
*  XINTUV   - Interpolate midpoint-values to values in vertices of
*             a triangulation, allocate memory if necessary
*  BDRCOR   - Implement Dirichlet values into vertex-based
*             solution vectors
************************************************************************

************************************************************************
* Read solution from disc
*
* This reads a solution vector to disc. If INUM>=0, an additional
* file extension ".INUM" is added to the file name CNAME.
*
* In:
*   MUNIT  : Unit number that should be used for input.
*            Can be 0 to use default value.
*   CNAME  : Basic filename of the file; ignored if MUNIT<>0.
*   NEQ    : Number of equations in DUP
*   IFMT   : <=0: Unformattes input
*             >0: Formatted input
*   INUM   : Additional file extension identifier to add to the
*            filename.
*            =-1: Don't add an additional extension
*            >=0: Add ".INUM" to the filename
*
* Out:
*   DUP    : array [1..NEQ] of double
*            Solution vector
*   IERR   : =0, if reading was successfull
*            =1, if the file was not found
*            =2, if the array size does not match (not implemented yet)
*            =3, if there was a general error in reading the data
*
* Remark:
*   If MUNIT=0, Unit 69 will be opened for input and closed
*   afterwards. If MUNIT<>0, the output is written into that channel
*   without opening/closing it.
*   
************************************************************************

      SUBROUTINE PPRDVC (MUNIT,NEQ,DUP,INUM,IFMT,IERR,CNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      
      INCLUDE 'dstrings.inc'
      
      INTEGER NEQ,INUM,MUNIT,IFMT,F,IERR
      DOUBLE PRECISION DUP(NEQ)
      CHARACTER CNAME*(*)
      
      CHARACTER CFN*60,CARR*6
      INTEGER LSTR,M
      LOGICAL BEXISTS

C     Check "formatted"-tag
      
      IF (IFMT.LE.0) THEN
        F = 0
      ELSE
        F = 1
      END IF

C     If there does not exist a channel number, open a new file
C     and close it afterwards
      
      M=MUNIT
      
      IF (MUNIT.EQ.0) THEN
        LSTR = STNEWC (.TRUE.,CNAME)
        
        IF (INUM.GE.0) THEN
          CALL STCATC (LSTR,.FALSE.,'.')
          CALL STCATI (LSTR,INUM,1,.FALSE.)
        END IF
        
        CALL STPUT (LSTR, CFN)
        CALL STDIS (LSTR)
        
C       Cancel if the file does not exist

        INQUIRE (FILE=CFN,EXIST=BEXISTS)
        IF (.NOT.BEXISTS) THEN
          IER = 1
          RETURN
        END IF
        
        M = 69
        CALL OF0 (M,CFN,F)
      END IF
      
C     Read the vector
      
      CALL ORA1 (DUP,CARR,M,F)
      
      IERR = 0
      IF (IER.NE.0) IERR = 3
      
      IF (MUNIT.EQ.0) THEN
        CLOSE (M)
      END IF
      
      END
      
************************************************************************
* Write solution to disc
*
* This writes a solution vector to disc. If INUM>=0, an additional
* file extension ".INUM" is added to the file name CNAME.
*
* In:
*   MUNIT  : Unit number that should be used for output.
*            Can be 0 to use default value.
*   CNAME  : Basic filename of the file; ignored if MUNIT<>0.
*   NEQ    : Number of equations in DUP
*   DUP    : array [1..NEQ] of double
*            Solution vector
*   IFMT   : <=0: Unformattes output
*             >0: Formatted output
*   INUM   : Additional file extension identifier to add to the
*            filename.
*            =-1: Don't add an additional extension
*            >=0: Add ".INUM" to the filename
*
* Remark:
*   If MUNIT=0, Unit 69 will be opened for output and closed
*   afterwards. If MUNIT<>0, the output is written into that channel
*   without opening/closing it.
************************************************************************

      SUBROUTINE PPWRVC (MUNIT,NEQ,DUP,INUM,IFMT,CNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'dstrings.inc'
      
      INTEGER NEQ,INUM,MUNIT,IFMT,F
      DOUBLE PRECISION DUP(NEQ)
      CHARACTER CNAME*(*)
      
      CHARACTER CFN*(60)
      INTEGER LSTR,M

C     Check "formatted"-tag
      
      IF (IFMT.LE.0) THEN
        F = 0
      ELSE
        F = 1
      END IF

C     If there does not exist a channel number, open a new file
C     and close it afterwards
      
      M=MUNIT
      
      IF (MUNIT.EQ.0) THEN
        LSTR = STNEWC (.TRUE.,CNAME)
        
        IF (INUM.GE.0) THEN
          CALL STCATC (LSTR,.FALSE.,'.')
          CALL STCATI (LSTR,INUM,1,.FALSE.)
        END IF
        
        CALL STPUT (LSTR, CFN)
        CALL STDIS (LSTR)
        
        M = 69
        CALL OF0 (M,CFN,F)
      END IF
      
C     Write out that vector
      
      CALL OWA1 (DUP,'PPWRVC',NEQ,M,F)
      
      IF (MUNIT.EQ.0) THEN
        CLOSE (M)
      END IF
      
      END
      
************************************************************************
* Interpolate midpoint values to vertex values
*
* This routine interpolates a vector DUx(1..NMT) given as values in
* midpoints of edges to a vector DLx(1..NVT) of values in corners
* of elements. Bilinear interpolation is used. 
*
* After this routine is finished, the caller should correct Dirichlet
* boundary nodes by using the given analytic function to correct these
* nodes!
*
* In:
*   DU1,    
*   DU2    : array [1..NMT] of double
*            Vectors with values in the midpoints of edges
*   DAUX   : array [1..NVT] of double
*            Auxiliary array.
*   KVERT,
*   KMID,
*   NVT,
*   NEL    : Usual geometry information
*
* Out:
*   DL1,
*   DL2    : array [1..NVT] of double
*            Vector with values in the corners of elements
************************************************************************

      SUBROUTINE INTUV (DU1,DU2,DAUX,KVERT,KMID,NVT,NEL,DL1,DL2)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters      

      DOUBLE PRECISION DU1(*),DU2(*),DL1(*),DL2(*),DAUX(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER NVT,NEL
      
C local variables

      INTEGER IEL
      INTEGER IM1,IM2,IM3,IM4,IV1,IV2,IV3,IV4,IV
      DOUBLE PRECISION DUH1,DUH2,DUH3,DUH4,DVH1,DVH2,DVH3,DVH4

C     Clear the target vectors

      CALL LCL1 (DL1,NVT)
      CALL LCL1 (DL2,NVT)
      CALL LCL1 (DAUX,NVT)
      
C     In DAUX we sum up the number of adjacent edges to a vertex.
C
C     Loop over the elements:
      
      DO IEL=1,NEL
      
        IM1=KMID(1,IEL)-NVT
        IM2=KMID(2,IEL)-NVT
        IM3=KMID(3,IEL)-NVT
        IM4=KMID(4,IEL)-NVT

        IV1=KVERT(1,IEL)
        IV2=KVERT(2,IEL)
        IV3=KVERT(3,IEL)
        IV4=KVERT(4,IEL)

        DUH1=DU1(IM1)
        DUH2=DU1(IM2)
        DUH3=DU1(IM3)
        DUH4=DU1(IM4)

        DVH1=DU2(IM1)
        DVH2=DU2(IM2)
        DVH3=DU2(IM3)
        DVH4=DU2(IM4)

        DAUX(IV1)=DAUX(IV1)+1D0
        DAUX(IV2)=DAUX(IV2)+1D0
        DAUX(IV3)=DAUX(IV3)+1D0
        DAUX(IV4)=DAUX(IV4)+1D0

C       Bilinear interpolation gives what we have to add to the
C       value in a corner:

        DL1(IV1)=DL1(IV1) + 0.75D0*(DUH1+DUH4) - 0.25D0*(DUH2+DUH3)
        DL1(IV2)=DL1(IV2) + 0.75D0*(DUH2+DUH1) - 0.25D0*(DUH3+DUH4)
        DL1(IV3)=DL1(IV3) + 0.75D0*(DUH3+DUH2) - 0.25D0*(DUH4+DUH1)
        DL1(IV4)=DL1(IV4) + 0.75D0*(DUH4+DUH3) - 0.25D0*(DUH1+DUH2)

        DL2(IV1)=DL2(IV1) + 0.75D0*(DVH1+DVH4) - 0.25D0*(DVH2+DVH3)
        DL2(IV2)=DL2(IV2) + 0.75D0*(DVH2+DVH1) - 0.25D0*(DVH3+DVH4)
        DL2(IV3)=DL2(IV3) + 0.75D0*(DVH3+DVH2) - 0.25D0*(DVH4+DVH1)
        DL2(IV4)=DL2(IV4) + 0.75D0*(DVH4+DVH3) - 0.25D0*(DVH1+DVH2)

      END DO

C     Divide by the number of adjacent vertices, this results
C     in the interpolated solution.

      DO IV=1,NVT
        DL1(IV)=DL1(IV)/DAUX(IV)
        DL2(IV)=DL2(IV)/DAUX(IV)
      END DO
 
99999 END

************************************************************************
* Interpolate midpoint values to vertex values,
* allocate memory if necessary.
*
* This routine calls INTUP to interpolate midpoint-based functions
* to vertex-based functions. If necessary, memory is allocated
* before.
*
* In:
*   DU1,    
*   DU2    : array [1..NMT] of double
*            Vectors with values in the midpoints of edges
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation structure.
*   LX,
*   LY     : Handle to array [1..NVT] of double
*            Vertex-based velocity vectors.
*            = 0: allocate new vectors
*            <>0: overwrite previous vectors
*
* Out:
*   If LX=0 / LY=0:
*     LX,
*     LY     : Handle to array [1..NVT] of double
*              Vertex-based velocity vectors.
*   The vectors identified by LX/LY are filled with data.
************************************************************************

      SUBROUTINE XINTUV (DU1,DU2,TRIA,LX,LY)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      
      DOUBLE PRECISION DU1(*),DU2(*)
      INTEGER TRIA(SZTRIA),LX,LY
      
      INTEGER J
      
C     Allocate memory if necessary
      
      IF (LX.EQ.0) CALL ZNEW (TRIA(ONVT),1,LX,'DX    ')
      IF (LY.EQ.0) CALL ZNEW (TRIA(ONVT),1,LY,'DY    ')
      
      CALL ZNEW (TRIA(ONVT),1,J,'DAUX  ')
      
C     Call the calculation routine

      CALL INTUV (DU1,DU2,DWORK(L(J)),KWORK(L(TRIA(OLVERT))),
     *            KWORK(L(TRIA(OLMID))),TRIA(ONVT),TRIA(ONEL),
     *            DWORK(L(LX)),DWORK(L(LY)))
      
C     Release unused memory
      
      CALL ZDISP (0,J,'DAUX  ')
      
      END
      
************************************************************************
* Interpolate pressure values to vertex values
*
* This routine takes the mean of the pressure of elements adjacent
* to a vertex and uses that as pressure value in the vertex.
*
* In:
*   DP     : array [1..NEL] of double
*            Pressure values in all cells
*   DAUX   : array [1..NVT] of double
*            Auxiliary array
*   DAREA  : array [1..NEL] of double
*            Area of all cells
*   KVERT,
*   NVT,
*   NEL    : Usual geometry information
*
* Out:
*   DPL    : array [1..NVT] of double
*            Pressure, linearly interpolated into the vertices.
*            The pressure is calculated as weighted mean using the
*            element area.
************************************************************************

      SUBROUTINE INTPV (DP,DPL,DAUX,AREA,KVERT,NVT,NEL)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      DOUBLE PRECISION DP(*),DPL(*),DAUX(*)
      INTEGER KVERT(NNVE,*),NEL,NVT
      DOUBLE PRECISION AREA(*)

C local variables

      INTEGER IVT,IEL
      INTEGER IV1,IV2,IV3,IV4

      DOUBLE PRECISION DPIEL, DAREA

C     Clear the target/auxiliary arrays
      
      CALL LCL1 (DPL,NVT)
      CALL LCL1 (DAUX,NVT)     

C     Loop over the elements

      DO IEL=1,NEL

        DPIEL=DP(IEL)
        DAREA=AREA(IEL)

        IV1=KVERT(1,IEL)
        IV2=KVERT(2,IEL)
        IV3=KVERT(3,IEL)
        IV4=KVERT(4,IEL)

C       Calculate the weighted pressure into DPL

        DPL(IV1)=DPL(IV1)+0.25D0*DAREA*DPIEL
        DPL(IV2)=DPL(IV2)+0.25D0*DAREA*DPIEL
        DPL(IV3)=DPL(IV3)+0.25D0*DAREA*DPIEL
        DPL(IV4)=DPL(IV4)+0.25D0*DAREA*DPIEL

C       Calculate the total area of the cells arount each
C       vertex into DAUX

        DAUX(IV1)=DAUX(IV1)+0.25D0*DAREA
        DAUX(IV2)=DAUX(IV2)+0.25D0*DAREA
        DAUX(IV3)=DAUX(IV3)+0.25D0*DAREA
        DAUX(IV4)=DAUX(IV4)+0.25D0*DAREA

      END DO
      
C     Divide the calculated DPL-value by the area around each
C     cell to calculate the mean pressure there.
      
      DO IVT=1,NVT
        DPL(IVT)=DPL(IVT)/DAUX(IVT)
      END DO
      
      END
      
************************************************************************
* Interpolate pressure values to vertex values,
* Allocate memory if necessary
*
* This routine takes the mean of the pressure of elements adjacent
* to a vertex and uses that as pressure value in the vertex.
* If necessary, memory is allocated before.
*
* In:
*   DP     : array [1..NEL] of double
*            Pressure values in all cells
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation structure.
*   LP     : Handle to array [1..NVT] of double
*            Vertex-based pressure vector.
*            = 0: allocate new vectors
*            <>0: overwrite previous vectors
*
* Out:
*   If LX=0 / LY=0:
*     LP     : Handle to array [1..NVT] of double
*              Vertex-based pressure vector.
*   The vector identified by LP is filled with data.
************************************************************************

      SUBROUTINE XINTPV (DP,TRIA,LP)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'
      
C     parameters
      
      DOUBLE PRECISION DP(*)
      INTEGER TRIA(SZTRIA),LP
      
C     local variables
       
      INTEGER J      

C     Allocate memory if necessary
      
      IF (LP.EQ.0) CALL ZNEW (TRIA(ONVT),1,LP,'DP    ')
      
      CALL ZNEW (TRIA(ONVT),1,J,'DAUX  ')
      
C     Call the calculation routine

      CALL INTPV (DP,DWORK(L(LP)),DWORK(L(J)),DWORK(L(TRIA(OLAREA))),
     *            KWORK(L(TRIA(OLVERT))),TRIA(ONVT),TRIA(ONEL))
      
C     Release unused memory
      
      CALL ZDISP (0,J,'DAUX  ')
      
      END
      
************************************************************************
* Implement velocity Dirichlet values into vertices
*
* This routine calculates the Dirichlet values of all vertices
* and implements them into the vectors DL1,DL2 for the X- and
* Y-velocity. The routine handles only the corners
* of the elements. DL1,DL2 are expected to be solution vectors
* in element corners 1..NVT.
*
* In:
*   DCORVG : usual geometry information; must correspond to TRIA
*   KXNPR  : array with extended nodal property for deciding on
*            which node is Dirichlet; must correspond to TRIA
*   UE     : SUBROUTINE
*            Gives the exact solution in Dirichlet values
*   TRIA   : used triangulation
*   TIMENS,
*   RE,
*   IPARAM,
*   DPARAM : information passed to UE to evaluate the Dirichlet value
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to fictitious boundary
*            routines. Not used in this routine.
*
* Out:
*   DL1,
*   DL2    : array [1..NVT] of double
*            All vertices corresponding to Dirichlet vertices are
*            recalculated using UE and the user defined fictitious
*            boundary input routines. 
************************************************************************

      SUBROUTINE BDRCOR (DL1,DL2,TRIA,DCORVG,KXNPR,
     *                   UE,TIMENS,RE,IPARAM,DPARAM,IGEOM,DGEOM)

      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
C     parameters      

      DOUBLE PRECISION DL1(*),DL2(*),TIMENS,RE,DPARAM(*)
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      
      INTEGER TRIA(SZTRIA),IPARAM(*)
      DOUBLE PRECISION UE,FBINDT
      EXTERNAL UE,FBINDT
      
      INTEGER KXNPR(2,*)
      DOUBLE PRECISION DCORVG(2,*)

C     local variables

      INTEGER I
      DOUBLE PRECISION X,Y
      
C     Loop through all vertices

      DO I=1,TRIA(ONVT)
       
C       Is the vertex of Dirichlet-type?

        IF (IAND(KXNPR(1,I),2**12).NE.0) THEN
          
C         Get the coordinates

          X=DCORVG(1,I)
          Y=DCORVG(2,I)
          
C         Fictitious boundary vertex?

          IF (IAND(KXNPR(1,I),2**13).NE.0) THEN

            DL1(I) = FBINDT (1,X,Y,TIMENS,RE,I,TRIA,IPARAM,DPARAM,
     *                       IGEOM,DGEOM)
            DL2(I) = FBINDT (2,X,Y,TIMENS,RE,I,TRIA,IPARAM,DPARAM,
     *                       IGEOM,DGEOM)

          ELSE
          
C           No, standard vertex;   
C           Call the Dirichlet routine to calculate the value
            
            DL1(I) = UE(X,Y,1,TIMENS,RE,I,TRIA,IPARAM,DPARAM,
     *                  IGEOM,DGEOM)
            DL2(I) = UE(X,Y,2,TIMENS,RE,I,TRIA,IPARAM,DPARAM,
     *                  IGEOM,DGEOM)
     
          END IF
          
        END IF
      
      END DO ! I

      END

************************************************************************
* Implement pressure Dirichlet values into vertices
*
* This routine calculates the Dirichlet values of all vertices
* and implements them into the vector DLP pressure. 
* The routine handles only the corners
* of the elements. DLP are expected to be solution vectors
* in element corners 1..NVT.
*
* In:
*   DCORVG : usual geometry information; must correspond to TRIA
*   KXNPR  : array with extended nodal property for deciding on
*            which node is Dirichlet; must correspond to TRIA
*   PE     : SUBROUTINE
*            Gives the exact solution in Dirichlet values
*   TRIA   : used triangulation
*   TIMENS,
*   RE,
*   IPARAM,
*   DPARAM : information passed to UE to evaluate the Dirichlet value
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to fictitious boundary
*            routines. Not used in this routine.
*
* Out:
*   DLP    : array [1..NVT] of double
*            All vertices corresponding to Dirichlet vertices are
*            recalculated using PE and the user defined fictitious
*            boundary input routines. 
************************************************************************

      SUBROUTINE BDRCRP (DLP,TRIA,DCORVG,KXNPR,
     *                   PE,TIMENS,RE,IPARAM,DPARAM,IGEOM,DGEOM)

      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
C     parameters      

      DOUBLE PRECISION DLP(*),TIMENS,RE,DPARAM(*)
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      
      INTEGER TRIA(SZTRIA),IPARAM(*)
      DOUBLE PRECISION PE,FBINDT
      EXTERNAL PE,FBINDT
      
      INTEGER KXNPR(2,*)
      DOUBLE PRECISION DCORVG(2,*)

C     local variables

      INTEGER I
      DOUBLE PRECISION X,Y
      
C     Loop through all vertices

      DO I=1,TRIA(ONVT)
       
C       Is the vertex of Dirichlet-type?

        IF (IAND(KXNPR(1,I),2**12).NE.0) THEN
          
C         Get the coordinates

          X=DCORVG(1,I)
          Y=DCORVG(2,I)
          
C         Fictitious boundary vertex?

          IF (IAND(KXNPR(1,I),2**13).NE.0) THEN

            DLP(I) = FBINDT (3,X,Y,TIMENS,RE,I,TRIA,IPARAM,DPARAM,
     *                       IGEOM,DGEOM)

          ELSE
          
C           No, standard vertex;   
C           Call the Dirichlet routine to calculate the value
            
            DLP(I) = PE(X,Y,1,TIMENS,RE,I,TRIA,IPARAM,DPARAM,
     *                  IGEOM,DGEOM)
     
          END IF
          
        END IF
      
      END DO ! I

      END

************************************************************************
* Calculate streamfunction
*
* Based on a (U,V) solution vector, this routine routine calculates
* the (scalar) streamfunction in all vertices of the triangulation.
* U,V are expected as point values in the corners of the elements
* in a given triangulation.
*
* In:
*   DU1,
*   DU2    : array [1..NMT] of double
*            X- and Y-velocity field, midpoint/edge-based
*   DCORVG,
*   KVERT,
*   KMID,
*   KADJ,
*   NVT,
*   NVE,
*   NEL    : Usual geometry information; must correspond to TRIA!
*   DVIND  : array [1..NVT] of double
*            Auxiliary array.
*
* Out:
*   DX     : array [1..NVT] of double precision
*            The stream function in all (corner-) vertices.
************************************************************************

      SUBROUTINE U2ISO (DCORVG,KVERT,KMID,KADJ,NVT,NEL,NVE,
     *                  DVIND,DX,DU1,DU2)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'

C     parameters

      DOUBLE PRECISION DCORVG(2,*),DVIND(*),DX(*),DU1(*),DU2(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*)
      INTEGER NVT,NEL,NVE

C local variables

      INTEGER IVT, IEH, IEL, IH, IVE, JVE, IHV, IND, INDH, IVTH, IMID
      INTEGER IME, IELH
      DOUBLE PRECISION PX1,PY1,PX2,PY2,DN1,DN2,DXH

C     Clear the solution and auxiliary array

      CALL LCL1(DX,NVT)
      CALL LCL1(DVIND,NVT)

C     Start with element 1 and its first vertex.
C     Set the streamfunction there to 0.0. The streamfunction
C     in the other vertices are then relatively calculated 
C     to this basis.
C     DVIND is a marker which is set to 1 for every node that 
C     is finished.

      DX(KVERT(1,1))    = 0D0
      DVIND(KVERT(1,1)) = 1D0

C     Loop over the elements:

      DO IEH=1,NEL

C       We set the current element IEL to IEH:

        IEL=IEH
        IH=0

C       On the current element, loop over the vertices of
C       that element. Add the DVIND-values of all vertices of the
C       current element together. So, IH gets the number of marked
C       vertices and IND will be the number of the "last" marked
C       vertex of the four on the element.

        DO IVE = 1,NVE
          JVE=KVERT(IVE,IEL)
          IHV=INT(DVIND(JVE))
          IH=IH+IHV
          IF (IHV.GE.1) IND=IVE
        END DO

C       If all four vertices are marked, there's nothing to calculate
C       on the element. If no vertex is marked, we can't calculate
C       anything on the current element. In both cases skip the
C       computation and search for a better element:

        IF ((IH.GE.NVE).OR.(IH.EQ.0)) GOTO 20

C       Ok, here we found an element where some of the vertices are
C       marked and some not. Here we can calculate a part of the
C       streamfunction.

13      CONTINUE  
      
C       Loop over the vertices on the element:
      
        DO IVE=1,NVE-1

C         IND is the index of a marked vertex. Calculate the "next"
C         vertex following IND and its vertex number into IVH.

          INDH=MOD(IND,NVE)+1
          IVTH=KVERT(INDH,IEL)

C         If that vertex is not marked, the streamfunction is not
C         calculated there. Otherwise we are just looking at two 
C         marked neighboured vertices, so there's nothing to gain here. 

          IF (DVIND(IVTH).LT.1D0) THEN
          
C           Vertex IVT (corresponding to IND) is marked, vertex IVTH 
C           is not. 
C
C               x---x IVTH
C               |   |
C               x---O IVT
C
C           Mark vertex IVTH to indicate that the streamfunction
C           is not being calculated there:
          
            DVIND(IVTH)=1D0

C           and calculate the streamfunction in INTH.

            IVT =KVERT(IND,IEL)
            IMID=KMID (IND,IEL)-NVT
            
C           IMID is now the midpoint number following IVT and thus
C           the number of the DOF in the FE function.
C           Calculate the normal vector of the current edge into
C           N=(DN1,DN2) - not normalized.
C
C               x-------x IVTH
C               |       |
C               |  IMID x--> (DN1,DN2)
C               |       |
C               x-------O IVT

            PX1=DCORVG(1,IVT)
            PY1=DCORVG(2,IVT)
            PX2=DCORVG(1,IVTH)
            PY2=DCORVG(2,IVTH)
            DN1 = PY2-PY1
            DN2 =-PX2+PX1
            
C           Calculate the streamfunction in IVTH from the value
C           in IVT by:
C
C           sfc(IVTH) = sfc(IVT) + U(IMID)*N
C
C           which is the "amount of flow over the edge (IVT,IVTH)".

            DX(IVTH)=DX(IVT)+(DU1(IMID)*DN1+DU2(IMID)*DN2)
          
          END IF ! (DVIND(IVTH) < 1D0)

C         Go on to the next vertex on the element to look if that one
C         has a not marked neighbour.

          IND=INDH
            
        END DO ! IVE
        
C       Now on the current element IEL, on all (corner) vertices the
C       streamfunction is calculated. We go on looking to the adjacent
C       elements of IEL if there's an element where the streamfunction
C       is not calculated in all vertices...

20      CONTINUE

C       Look onto the adjacent elements of the current element if there's
C       a suitable neighbour element where we can continue the calculation.
C
C       Loop over the edges of the current element

        DO IME=1,NVE

C         Get the neighbour element adjacent to the current element

          IELH=KADJ(IME,IEL)
          
          IF (IELH.NE.0) THEN
          
C           Now we have the number of the neighbour element in IELH.
C           Loop about the vertices of the element and sum up the
C           markers into IH.
          
            IH=0
            DO IVE=1,NVE
              JVE=KVERT(IVE,IELH)
              IHV=INT(DVIND(JVE))
              IH=IH+IHV
              IF (IHV.GE.1) INDH=IVE
            END DO
            
C           If there is at least one but not all markers set, the
C           element can be used for further calculation.
C           Switch the current element IEL to that one and
C           continue the calculation there.

            IF ((IH.LT.NVE).AND.(IH.GT.0)) THEN
              IEL=IELH
              IND=INDH
              GOTO 13
            END IF
            
          END IF ! IELH <> 0

        END DO ! IME
      
      END DO ! IEH

C     At last, normalize the streamfunction such that vertex 1
C     has value 0.0. Remember that we assigned a value of 0.0
C     to the first vertex of element 1, which is usually not
C     vertex 1 of the triangulation!

      DXH=DX(1)
      DO IVT=1,NVT
        DX(IVT)=DX(IVT)-DXH
      END DO

      END

************************************************************************
* Calculate streamfunction,
* Allocate memory if necessary
*
* This routine calls U2ISO to calculate a streamfunction of a given
* solution vector. If necessary, memory is allocated before.
*
* In:
*   DU,
*   DV     : array [1..NMT] of double
*            Velocity field in the midpoints of the elements
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation structure.
*   LISO   : Handle to array [1..NVT] of double
*            Vertex-based streamfunction vector.
*            = 0: allocate new vector
*            <>0: overwrite previous vector
*
* Out:
*   If LISO=0:
*     LISO   : Handle to array [1..NVT] of double
*              Vertex-based streamfunction vector.
*   The vector identified by LISO is filled with data.
************************************************************************

      SUBROUTINE XU2ISO (DU,DV,TRIA,LISO)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'
      
      DOUBLE PRECISION DU(*),DV(*)
      INTEGER TRIA(SZTRIA),LISO
      
C     local variables
       
      INTEGER J      

C     Allocate memory if necessary
      
      IF (LISO.EQ.0) CALL ZNEW (TRIA(ONVT),1,LISO,'DISO  ')
      
      CALL ZNEW (TRIA(ONVT),1,J,'DAUX  ')
      
C     Call the calculation routine

      CALL  U2ISO (DWORK(L(TRIA(OLCORVG))),KWORK(L(TRIA(OLVERT))),
     *             KWORK(L(TRIA(OLMID))),KWORK(L(TRIA(OLADJ))),
     *             TRIA(ONVT),TRIA(ONEL),TRIA(ONVE),
     *             DWORK(L(J)),DWORK(L(LISO)),
     *             DU,DV)
      
C     Release unused memory
      
      CALL ZDISP (0,J,'DAUX  ')
      
      END
