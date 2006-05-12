************************************************************************
* This file collects basic input/output routines for GMV files.
* There are various routines to open a GMV-file, to write triangulation
* or solution data to it, to write polynoms to it and so on.
*
* How to write a GMV-file? Well, only some quick commands:
*
*   CALL GMVOF0      - open a gmv file for writing
*   CALL GMVHEA      - write the header
*   CALL GMVTRI      - write the triangulation
*
* Then (in arbitrary order and perhaps multiple times if necessary
*   CALL GMVMAT      - Write materials of cells or vertices
*   CALL GMVVEL      - Write velocity field
*   CALL GMVSCA      - Write scalar function
*   CALL GMVFLG      - Write cell/vertex flags
*   CALL GMVPOL      - Write polygon
*   CALL GMVCPL      - Write cell based polygon
*   CALL GMVTIM      - Write probe time
*
* At last:
*   CALL GMVFOT      - Write footer
*   CALL CLOSE       - Close the file
************************************************************************

************************************************************************
* Open a GMV file for writing
*
* This routine uses OF0 to open a GMV file for writing. 
* If IDX>=0, an index ".x" will be appended to the filename.
*
* In:
*   CNAME   - Name of the GMV file
*   IDX     - IDX=-2 will use the filename CFILE as it is.
*             IDX=-1 will append '.gmv' to the filename, thus naming
*             the file CFILE.gmv.
*             IDX>=0 will append '.gmv.[idx]' to the filename in
*             compliance to the auto-read feature of GMV.
*   MFILE   - The number of a file system handle that should represent
*             the GMV file.
*
* Out:
*   The GMV-file is opened as unit MFILE for writing.
*   It must be closed by the caller by using CLOSE.
************************************************************************

      SUBROUTINE GMVOF0 (MFILE,IDX,CNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'dstrings.inc'
      
      INTEGER MFILE,IDX
      CHARACTER*(*) CNAME
      
      INTEGER LSTR,IFMT
      CHARACTER FNAME*(60)
      
C     Create the filename using the string library.
C     Create the file in the GMV directory:

      LSTR = STNEWC (.TRUE.,CNAME)
      
      IF (IDX.GE.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.gmv')
      END IF
      
C     Append the suffix if there is one

      IF (IDX.GE.0) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IDX,5,.TRUE.)
      END IF
      
C     Transfer the final filename into FNAME

      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)
      
C     Open the file

      IFMT = 1
      CALL OF0 (MFILE,FNAME,IFMT)
      
      END
      
************************************************************************
* Write GMV header
*
* This writes the general header of a GMV file.
*
* In:
*   MFILE : The file handle of the GMV file.
************************************************************************
      
      SUBROUTINE GMVHEA (MFILE)
      
      IMPLICIT NONE
      
      INTEGER MFILE
      
      WRITE (MFILE,'(A)') 'gmvinput ascii'
      
      END
      
************************************************************************
* Write GMV header
*
* This writes the general footer of a GMV file.
*
* In:
*   MFILE : The file handle of the GMV file.
************************************************************************
      
      SUBROUTINE GMVFOT (MFILE)
      
      IMPLICIT NONE
      
      INTEGER MFILE
      
      WRITE (MFILE,'(A)') 'endgmv'
      
      END
      
************************************************************************
* Write GMV triangulation
*
* This writes nodes and cells of a triangulation structure to a GMV 
* file, thus establishing the basic shape of the domain.
*
* In:
*   MFILE : The file handle of the GMV file; must be open.
*   TRIA  : array [1..SZTRIA] of integer
*           Triangulation structure to be written into the GMV file.
*   IEDG  : Method of writing the edges and cells.
*           =0: Don't write information about edges
*           =1: Write points in DCORMG to GMV file as information about
*               edges; if DCORMG does not exist, construct edge
*               mitpoints to identify edges.
*           =2: Write midpoint of edges / DCORMG and midpoints of
*               elements to GMV file
* Out:
*   NCELLS: Number of cells in the GMV file.
*           Must be passed to the other output routines.
*   NVERTS: Number of vertices in the GMV file.
*           Must be passed to the other output routines.
************************************************************************
      
      SUBROUTINE GMVTRI (MFILE,TRIA,IEDG,NCELLS,NVERTS)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER MFILE,TRIA(SZTRIA)
      
      INTEGER I,J,IEDG,NVTV,NMTV,NELV,NCELLS,NVERTS
      INTEGER KVERT,KCORVG
      DOUBLE PRECISION X,Y

C     Fetch some information for quicker access:

      KVERT  = L(TRIA(OLVERT))
      KCORVG = L(TRIA(OLCORVG))

C     We start with the nodes.
C     Calculate total number of nodes and write it into the GMV file.

      NVTV = TRIA(ONVT)
      NMTV = MAX(TRIA(ONMT),TRIA(ONVEDT))
      NELV = MAX(TRIA(ONEL),TRIA(ONIELV))
      
      J = NVTV
      IF (IEDG.GE.1) J=J+NMTV
      IF (IEDG.GE.2) J=J+NELV
      
      WRITE(MFILE,*)'nodes ',J
      
C     Return the number of vertices in the GMV file to the caller

      NVERTS = J
      
C     All X-coordinates
      
      DO I=0,NVTV-1
        WRITE(MFILE,'(E15.8)') REAL(DWORK(KCORVG+2*I))
      END DO
      
C     Should we write the X-coordinates of the nodes in DCORMG
C     to the GMV file?

      IF (IEDG.GE.1) THEN
        DO I=1,NMTV
          CALL NDE2XY (I+NVTV,TRIA,X,Y)  
          WRITE(MFILE,'(E15.8)') REAL(X)
        END DO
      END IF

C     Should we write the X-coordinates of element midpoints
C     to the GMV file?

      IF (IEDG.GE.2) THEN
        DO I=1,NELV
          CALL NDE2XY (I+NVTV+NMTV,TRIA,X,Y)  
          WRITE(MFILE,'(E15.8)') REAL(X)
        END DO
      END IF

C     all Y-coordinates
      
      DO I=0,NVTV-1
        WRITE(MFILE,'(E15.8)') REAL(DWORK(KCORVG+2*I+1))
      END DO
      
C     Should we write the Y-coordinates of the nodes in DCORMG
C     to the GMV file?

      IF (IEDG.GE.1) THEN
        DO I=1,NMTV
          CALL NDE2XY (I+NVTV,TRIA,X,Y)
          WRITE(MFILE,'(E15.8)') REAL(Y)
        END DO
      END IF

C     Should we write the Y-coordinates of element midpoints
C     to the GMV file?

      IF (IEDG.GE.2) THEN
        DO I=1,NELV
          CALL NDE2XY (I+NVTV+NMTV,TRIA,X,Y)  
          WRITE(MFILE,'(E15.8)') REAL(Y)
        END DO
      END IF

C     and 0 as Z-coordinate - so we produce a 2D-mesh
      
      DO I=1,NVTV
        WRITE(MFILE,'(E15.8)') 0E0
      END DO

C     Should we write 0's for the Z-coordinates of the nodes in DCORMG
C     to the GMV file?

      IF (IEDG.GE.1) THEN
        DO I=0,NMTV-1
          WRITE(MFILE,'(E15.8)') 0E0
        END DO
      END IF

C     Should we write 0's for the Z-coordinates of element midpoints
C     to the GMV file?

      IF (IEDG.GE.2) THEN
        DO I=0,NELV-1
          WRITE(MFILE,'(E15.8)') 0E0
        END DO
      END IF

C     Write the connectivity to the mesh - i.e. the cells.

      WRITE(MFILE,*)'cells ',TRIA(ONEL)
      
      IF (TRIA(ONVE).EQ.4) THEN
        DO I=0,TRIA(ONEL)-1
          WRITE(MFILE,*)'quad 4'
          WRITE(MFILE,'(4I8)') (KWORK(KVERT+I*NNVE+J),J=0,3)
        END DO
      ELSE
        DO I=0,TRIA(ONEL)-1
          WRITE(MFILE,*)'tri 3'
          WRITE(MFILE,'(3I8)') (KWORK(KVERT+I*NNVE+J),J=0,2)
        END DO
      END IF      
      
C     That's it for the connectivity.
C      
C     Return the number of cells to the caller

      NCELLS = TRIA(ONEL)

      END
      
************************************************************************
* Default GET routine for GMV flag/material vectors, version 1
*
* For writing materials or flags, the corresponding GMV output
* routines expect a callback routine. In case all vertices/elements
* have no material, this routine can be used which simply returns
* materian number 1.
* 
* In:
*   TRIA   - triangulation structure, not used
*   IENTRY - not used
*   IPARAM - not used
*   DPARAM - not used
* 
* Out:
*   M      - IPARAM(IENTRY)
************************************************************************

      SUBROUTINE GMVDG1 (TRIA,IENTRY,IPARAM,DPARAM,M)
      
      IMPLICIT NONE
      
      INTEGER TRIA(*),IENTRY,IPARAM(*),M
      DOUBLE PRECISION DPARAM(*)
      
      M = 1
      
      END

************************************************************************
* Default GET routine for GMV flag/material vectors, version 2
*
* For writing materials or flags, the corresponding GMV output
* routines expect a callback routine. In case the result for this
* routine is already calculated as an integer vector, GMVDGT can
* be used as simple callback routine which simply returns
* the IENTRY'th element of IPARAM as material/flag.
* 
* In:
*   TRIA   - triangulation structure, not used
*   IENTRY - number of entry which value should be extracted from
*            the array IPARAM
*   IPARAM - array [1..*] of integer
*            Assumed to point to a vector which contains for each
*            entry the routine is called for the corresponding
*            property
*   DPARAM - not used
* 
* Out:
*   M      - IPARAM(IENTRY)
************************************************************************

      SUBROUTINE GMVDG2 (TRIA,IENTRY,IPARAM,DPARAM,M)
      
      IMPLICIT NONE
      
      INTEGER TRIA(*),IENTRY,IPARAM(*),M
      DOUBLE PRECISION DPARAM(*)
      
      M = IPARAM(IENTRY)
      
      END
      
************************************************************************
* Write GMV materials
*
* This routine assignes all cells or vertices in a given triangulation 
* a defined material. For that purpose, a callback function is called 
* that has to calculate for every cell a material number.
*
* If materials are to be assigned to vertices as well as cells,
* exactly the same number of materials with exactly the same material
* names in the same order have to be used!
*
* Material names assigned here are also valid for output of polygons.
*
* In:
*   MFILE : The file handle of the GMV file; must be open.
*   TRIA  : array [1..SZTRIA] of integer
*           Triangulation structure.
*   ITYP  : Type of output:
*           =0: Cell based material; GETMAT is called for every cell
*           =1: Vertex-based material; GETMAT is called for every vertex
*   NNODES: Total number of cells / vertices in the GMV-file.
*           This must be either NCELLS or NVERTS (depending on ITYP),
*           which was returned by a previous call to GMVTRI.
*   NODES : Either the number of cells or vertices, the caller wants to
*           specify the material. GETMAT is called with IENTRY=1..NODES.
*   NMAT  : Number of materials with names in CMAT; <= 1000
*           If < 0: |NMAT| = number of materials with standard names.
*                   CMAT is not respected.
*   CMAT  : If NMAT>0: array [1..NMAT] of string(32)
*           Material names for material 1..NMAT.
*           If NMAT<0: not respected.
*   GETMAT: SUBROUTINE (TRIA,IENTRY,IPARAM,DPARAM,MAT)
*           Has to return for every element or vertex IENTRY the
*           routine is called for a material number 
*           MAT(IEL/IVT)=1..NMAT.
*           Standard return value = 1.
*   IPARAM,
*   DPARAM: User defined integer/double precision parameter block,
*           passed to the callback function
*
* Remarks: 
*   If there exists an integer vector containing the material for
*   all cells/vertices, GMVDGT can be used as default callback routine
*   in conjunction with IPARAM=that integer vector.
*
*   NODES must not necessarily be the total number of cells or vertices
*   in the triangulation. In this case, some cells/vertices remain
*   unspecified.
************************************************************************
      
      SUBROUTINE GMVMAT (MFILE,TRIA,ITYP,NNODES,NODES,
     *                   NMAT,CMAT,GETMAT,IPARAM,DPARAM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER MFILE,TRIA(SZTRIA),NMAT,IPARAM(*),ITYP,NNODES,NODES
      DOUBLE PRECISION DPARAM(*)
      EXTERNAL GETMAT
      CHARACTER CMAT(NMAT)*32
      
      INTEGER I,MAT,N2
      
      INTEGER N,NVTV,NMTV,NELV
                        
C     Set N to the number of entries in the vector we are allowed
C     to write out:

      N = MIN(NODES,NNODES)
      IF (N.EQ.0) RETURN
      
C     Set N2 to the number of zeros we have to append

      N2 = MAX(0,NNODES-NODES)
                        
C     Write material header to GMV file

      WRITE (MFILE,'(A,I5,I5)') 'material ',MAX(0,MIN(1000,ABS(NMAT))),
     *                          MIN(ITYP,1)
      
C     Write the material names to the file - either the given names
C     or the standard names
      
      IF (NMAT.GT.0) THEN
        DO I=1,MAX(0,MIN(1000,NMAT))
          WRITE (MFILE,'(A32)') CMAT(I)
        END DO
      ELSE
        DO I=1,MAX(1,MIN(1000,ABS(NMAT)))
          WRITE (MFILE,'(I4,A)') I,'_Material'
        END DO
      END IF
      
C     Loop through all elements and write their material to the file

      DO I=1,N
        CALL GETMAT (TRIA,I,IPARAM,DPARAM,MAT)
        WRITE (MFILE,'(I10)') MIN(MAT,1000)
      END DO
      
C     The other nodes remain unspecified

      DO I=1,N2
        WRITE (MFILE,'(I10)') 0
      END DO

      END
      
************************************************************************
* Write 2D velocity
*
* This writes a 2D velocity field to a GMV file. The velocity can be
* either cell-based or vertex-based and is given as two vectors for
* the X- and Y-velocity.
*
* In:
*   MFILE : The file handle of the GMV file; must be open.
*   TRIA  : array [1..SZTRIA] of integer
*           Underlying triangulation.
*   ITYP  : Type of output:
*           =0: Cell based; DU / DV = array [1..NEL] of double
*           =1: Vertex based; DU / DV = array [1..NVT] of double
*   NNODES: Total number of cells / vertices in the GMV-file.
*           This must be either NCELLS or NVERTS (depending on ITYP),
*           which was returned by a previous call to GMVTRI.
*   NEQ   : Number of entries in DU and DV.
*           If NODES<NNODES, the remaining entries are assumed to
*           be 0, thus writing 0-velocity to the GMV-file for these
*           nodes.
*           Standard values: =NEL (is ITYP=0), =NVT (if ITYP=1)
*   DU    : array [1..NEQ] of double
*           X-velocity
*   DV    : array [1..NEQ] of double
*           Y-velocity
************************************************************************
      
      SUBROUTINE GMVVEL (MFILE,TRIA,ITYP,NNODES,NEQ,DU,DV)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER MFILE,TRIA(SZTRIA),ITYP,NNODES,NEQ
      DOUBLE PRECISION DU(NEQ),DV(NEQ)
      
      INTEGER I,NVTV,NMTV,NELV
      
      INTEGER N,N2
                   
C     Set N to the number of entries in the vector we are allowed
C     to write out:

      N = MIN(NEQ,NNODES)
      IF (N.EQ.0) RETURN
      
C     Set N2 to the number of zeros we have to append

      N2 = MAX(0,NNODES-NEQ)
                        
C     First the header.
C     We only support face and vertex data.

      WRITE (MFILE,'(A,I1)') 'velocity ',MIN(ITYP,1)
      
C     X-Component

      DO I=1,N
        WRITE (MFILE,'(E15.8)') REAL(DU(I))
      END DO
      DO I=1,N2
        WRITE (MFILE,'(E15.8)') 0E0
      END DO

C     Y-Component

      DO I=1,N
        WRITE (MFILE,'(E15.8)') REAL(DV(I))
      END DO
      DO I=1,N2
        WRITE (MFILE,'(E15.8)') 0E0
      END DO
      
C     Z-Component = Filled with zero

      DO I=1,N
        WRITE (MFILE,'(E15.8)') 0E0
      END DO
      DO I=1,N2
        WRITE (MFILE,'(E15.8)') 0E0
      END DO
      
      END
      
************************************************************************
* Write scalar field
*
* This writes a scalar vector to a GMV file. The field can be
* either cell-based or vertex-based.
*
* In:
*   MFILE : The file handle of the GMV file; must be open.
*   TRIA  : array [1..SZTRIA] of integer
*           Underlying triangulation.
*   ITYP  : Type of output:
*           =0: Cell based
*           =1: Vertex based
*   NNODES: Total number of cells / vertices in the GMV-file.
*           This must be either NCELLS or NVERTS (depending on ITYP),
*           which was returned by a previous call to GMVTRI.
*   NEQ   : Number of entries in DVEC.
*           If NODES<NNODES, the remaining entries are assumed to
*           be 0, thus writing 0-values to the GMV-file for these
*           nodes.
*           Standard values: =NEL (is ITYP=0), =NVT (if ITYP=1)
*   DVEC  : array [1..*] of double
*           The vector to write out
*   CNAME : Name of the function
************************************************************************
      
      SUBROUTINE GMVSCA (MFILE,TRIA,ITYP,NNODES,NEQ,DVEC,CNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER MFILE,TRIA(SZTRIA),ITYP,NNODES,NEQ
      DOUBLE PRECISION DVEC(NEQ)
      CHARACTER CNAME*(*)
      
      INTEGER I,NVTV,NMTV,NELV
      
      INTEGER N,N2
                        
C     Set N to the number of entries in the vector we are allowed
C     to write out:

      N = MIN(NEQ,NNODES)
      IF (N.EQ.0) RETURN
      
C     Set N2 to the number of zeros we have to append

      N2 = MAX(0,NNODES-NEQ)
                        
C     First the header.

      WRITE (MFILE,'(A)') 'variable'

C     We only support face and vertex data.

      WRITE (MFILE,'(A,I5)') CNAME,MIN(ITYP,1)
      
C     Output the vector

      DO I=1,N
        WRITE (MFILE,'(E15.8)') REAL(DVEC(I))
      END DO
      DO I=1,N2
        WRITE (MFILE,'(E15.8)') 0E0
      END DO
      
C     Finish that part

      WRITE (MFILE,'(A)') 'endvars'

      END
      
************************************************************************
* Write GMV cell flag
*
* Every cell and/or vertex in the GMV file might be assigned one or 
* multiple integer flags that define the cell. This routine allowes to
* write out for every cell one flag consisting of multiple flag types.
*
* The flag which is specified for every cell is named CNAME. This flag
* is an integer value 1..NTYP, where every number is assigned a name
* CNAME(ityp) (e.g. typ "true"=1 and "false"=2 corresponding to a flag
* "boolean").
* Writing multiple flags for every cell with multiple types for every
* flags is thus possible.
*
* A callback function is called for every cell that has to return 
* for every cell an integer number defining the corresponding flag.
*
* In:
*   MFILE : The file handle of the GMV file; must be open.
*   TRIA  : array [1..SZTRIA] of integer
*           Triangulation structure.
*   CFLAG : The name of the flag
*   ITYP  : Type of output:
*           =0: Cell based; GETFLG is called for every cell
*           =0: Vertex based; GETFLG is called for every vertex
*   NNODES: Total number of cells / vertices in the GMV-file.
*           This must be either NCELLS or NVERTS (depending on ITYP),
*           which was returned by a previous call to GMVTRI.
*   NODES : Either the number of cells or vertices, the caller wants to
*           specify flags for. GETFLG is called with IENTRY=1..NODES.
*   NTYP  : Number of types for this flag; <= 1000
*   CTYPS : array [1..NMAT] of string(32)
*           Names of flag-types for typ 1..NFLAG
*   GETFLG: SUBROUTINE (TRIA,IENTRY,IPARAM,DPARAM,TPE)
*           Has to return for every cell IENTRY=1..NEL or for
*           every vertex IENTRY=1..NVT) the routine is
*           called for a flag type TPE(IENTRY)=1..NTYP.
*           Standard return value = 1.
*   IPARAM,
*   DPARAM: User defined integer/double precision parameter block,
*           passed to the callback function
*
* Remark: 
*   If there exists an integer vector containing the flag types for
*   all cells/vertices, GMVDGT can be used as default callback routine
*   in conjunction with IPARAM=that integer vector.
************************************************************************
      
      SUBROUTINE GMVFLG (MFILE,TRIA,ITYP,NNODES,NODES,NTYP,CTYPS,GETFLG,
     *                   IPARAM,DPARAM,CFLAG)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER MFILE,TRIA(SZTRIA),ITYP,NTYP,NNODES,NODES,IPARAM(*)
      DOUBLE PRECISION DPARAM(*)
      CHARACTER CTYPS(NTYP)*32,CFLAG*(*)
      
      INTEGER I,FLAG,N2
      
      INTEGER N
                        
C     Set N to the number of entries in the vector we are allowed
C     to write out:

      N = MIN(NODES,NNODES)
      IF (N.EQ.0) RETURN
      
C     Set N2 to the number of zeros we have to append

      N2 = MAX(0,NNODES-NODES)
                        
C     Write material header to GMV file

      WRITE (MFILE,'(A)') 'flags'
      
      WRITE (MFILE,'(A,I5,I1)') CFLAG,MAX(0,MIN(1000,NTYP)),MIN(ITYP,1)
      
C     Write flag types to the file
      
      DO I=1,NTYP
        WRITE (MFILE,'(A32)') CTYPS(I)
      END DO
      
C     Loop through all nodes and write their flag types to the file

      DO I=1,N
        CALL GETFLG (TRIA,I,IPARAM,DPARAM,FLAG)
        WRITE (MFILE,'(I10)') FLAG
      END DO
      
C     Append zeros for the remaining nodes

      DO I=1,N2
        WRITE (MFILE,'(I10)') 0
      END DO
      
C     Footer

      WRITE (MFILE,'(A)') 'endflag'

      END
      
************************************************************************
* Default GET routine for polygon points
*
* For writing polygons to a GMV file, the output routine expects
* a callback routine that returns the coordinates if the points.
* In case the coordinates already exist as an array of points,
* this routine can be used to return the coordinates from an array.
*
* The coordinate array is to be expected in DPARAM as a list
* of (X,Y)-points. GMVDGC then returns the IENTRY'th point
* from DPARAM.
* 
* In:
*   TRIA   - triangulation structure, not used
*   IENTRY - number of entry which value should be extracted from
*            the array IPARAM
*   IPARAM - not used
*   DPARAM - array [1..2,1..*] of double precision
*            List of (X,Y)-points of the polygon.
* 
* Out:
*   (X,Y)  - ( DPARAM(1,IENTRY), DPARAM(2,IENTRY) )
************************************************************************

      SUBROUTINE GMVDGC (TRIA,IENTRY,IPARAM,DPARAM,X,Y)
      
      IMPLICIT NONE
      
      INTEGER TRIA(*),IENTRY,IPARAM
      DOUBLE PRECISION DPARAM(2,*),X,Y
      
      X = DPARAM(1,IENTRY)
      Y = DPARAM(2,IENTRY) 
      
      END
      
************************************************************************
* Write polygon
*
* This writes a closed polygon to a GMV file. A callback function is
* called for every point that has to return the X- and Y-coordinate
* of that point.
*
* In:
*   MFILE : The file handle of the GMV file; must be open.
*   TRIA  : array [1..SZTRIA] of integer
*           Underlying triangulation.
*   IMAT  : Material number.
*           Must be defined by GMVMAT.
*   NVERT : Number of vertices; <= 65535.
*   GETPOL: SUBROUTINE (TRIA,IENTRY,IPARAM,DPARAM,X,Y)
*           Has to return for every entry 1..NVERT in the polygon
*           the X- any Y-coordinate of that point.
*
* If the points exist as array [1..2,1..NVERT] of point coordinates,
* GMVDGC can be used as default callback routine to work with
* that array.
************************************************************************
      
      SUBROUTINE GMVPOL (MFILE,TRIA,IMAT,NVERT,GETPOL,IPARAM,
     *                   DPARAM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER MFILE,TRIA(SZTRIA),IMAT,NVERT
      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)
      EXTERNAL GETPOL
      
      INTEGER I
      DOUBLE PRECISION X,Y
      
      INTEGER N
                  
      N=MIN(MAX(NVERT,1),65535)
                        
C     First the header.

      WRITE (MFILE,'(A)') 'polygons'

C     WE start a new polygon

      WRITE (MFILE,'(I5,I6)') IMAT,N
      
C     Output the X-components of the points

      DO I=1,N
        CALL GETPOL (TRIA,I,IPARAM,DPARAM,X,Y)
        WRITE (MFILE,'(E15.8)') REAL(X)
      END DO
      
C     and the Y-components

      DO I=1,N
        CALL GETPOL (TRIA,I,IPARAM,DPARAM,X,Y)
        WRITE (MFILE,'(E15.8)') REAL(Y)
      END DO
      
C     We don't have Z-components

      DO I=1,N
        WRITE (MFILE,'(E15.8)') 0E0
      END DO

C     Finish that part

      WRITE (MFILE,'(A)') 'endpoly'

      END
      
************************************************************************
* Write cell based polynoms
*
* This writes a polygon to a GMV file. A callback function is
* called for cell in the triangulation and has to return data of all
* polygons hitting the cell. 
*
* In:
*   MFILE : The file handle of the GMV file; must be open.
*   TRIA  : array [1..SZTRIA] of integer
*           Underlying triangulation.
*   GETPOL: SUBROUTINE (TRIA,IEL,IPOL,IMAT,IPARAM,DPARAM,X1,Y1,X2,Y2)
*           Has to return for every element IEL in the triangulation
*           data about the polygons in that cell; see below.
*
* Remarks
*   GETPOL is called multiple times to extract the information about
*   polygons. Cell-wise polygons consist one or more lines per
*   cell and an assigned material for every line. GETPOL is now called
*   to obtain the starting- and ending point of that line as well as 
*   the material. The exact procedure there is as follows:
*   1.) GETPOL is called with IPOL=0.
*       In this case GETPOL must return in IMAT to the number of
*       polygonal lines in that element.
*   2.) GETPOL is called with IPOL=1..number of polygonal lines
*       (as calculated in 1.)). For every polygonal line,
*       GETPOL has to return the number of the material in IMAT
*       (where the materials are defined in GMVMAT) as well as
*       the starting- and ending-point of that line in (X1,Y1)-(X2,Y2).
************************************************************************
      
      SUBROUTINE GMVCPL (MFILE,TRIA,NVERT,GETPOL,IPARAM,DPARAM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER MFILE,TRIA(SZTRIA),NVERT
      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)
      EXTERNAL GETPOL
      
      INTEGER I,J
      DOUBLE PRECISION X1,Y1,X2,Y2
      
      INTEGER N,IMAT
                  
C     First the header.

      WRITE (MFILE,'(A)') 'polygons'
      
C     Loop through all cells

      DO I=1,TRIA(ONEL)
      
C       How many lines are in the cell?
      
        CALL GETPOL (TRIA,I,0,N,IPARAM,DPARAM,0D0,0D0,0D0,0D0)
        
C       Loop through all these lines

        DO J=1,N

C         Get the line

          CALL GETPOL (TRIA,I,J,IMAT,IPARAM,DPARAM,X1,Y1,X2,Y2)
          
C         And write that line out as one polygon,
C         First X- and then Y-coordinates, no Z-coordinates.

          WRITE (MFILE,'(I5,I2)') IMAT,2
          
          WRITE (MFILE,'(6E15.6)') REAL(X1),REAL(X2),REAL(Y1),REAL(Y2),
     *                             0E0,0E0

        END DO ! J
      
      END DO ! I

C     Finish that part

      WRITE (MFILE,'(A)') 'endpoly'

      END
      
************************************************************************
* Write time
*
* This writes a double prec. time identifier to the GMV file.
*
* In:
*   MFILE : The file handle of the GMV file; must be open.
*   PTIME : Time that is to write into the GMV file
************************************************************************
      
      SUBROUTINE GMVTIM (MFILE,PTIME)

      IMPLICIT NONE
      
      INTEGER MFILE
      DOUBLE PRECISION PTIME
      
      WRITE (MFILE,'(A,E15.8)') 'probtime ',REAL(PTIME)
      
      END
      
