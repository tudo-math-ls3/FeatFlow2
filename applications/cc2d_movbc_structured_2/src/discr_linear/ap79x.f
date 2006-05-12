************************************************************************
* Generate matrix 7/8 structure, wrapper routine
*
* Extended calling convention; Works for all types of Finite Elements
*
* This routine calls AP7X to generate a the structure of a matrix.
*
* In:
*   TRIA  : array [1..SZTRIA] of integer
*           Triangulation structure of the underlying mesh
*   ELE   : Finite Element callback routine.
*           Can be triangular or quadrilateral, but must fit to
*           the triangulation in TRIA!
*   ISYMM : =0: don't respect symmetry, create matrix in structure 7
*           =1: matrix is symmetric; create matrix in structure 8
* Out:
*   LCOL  : Handle to KCOL array of the matrix
*   LLD   : Handle to KLD array of the matrix
*   NA    : Number of entries in the matrix
*   NEQ   : Dimension/Number of rows in the matrix
************************************************************************

      SUBROUTINE XAP7X(LCOL,LLD,NA,NEQ,TRIA,ELE,ISYMM)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'stria.inc'
      
C     parameters

      INTEGER LCOL,LLD,NA,NEQ,ISYMM,TRIA(SZTRIA)
      
      EXTERNAL ELE
      
C     local variables
      
      INTEGER IELTYP,IWMAX0,IFREE,LCOL1,LIND
      INTEGER KVERT,KMID
      
      INTEGER NDFGX
      EXTERNAL NDFGX

      SUB='XAP7'
      IF (ICHECK.GE.997) CALL OTRC('XAP7  ','12/11/89')

C     Call the element wrapper EA00 to determine the
C     element type (switches between triangular and quadrilateral
C     element depending on NVE)

      CALL EA00(ELE,ELE,TRIA(ONVE),IELTYP)
      
C     Get the number of equations via NDFGX
      
      NEQ=NDFGX(IELTYP,TRIA)
      IF (IER.NE.0) GOTO 99999

C     Allocate KLD on DWORK ***
      CALL ZNEW(NEQ+1,-3,LLD,'KLD   ')
      IF (IER.NE.0) GOTO 99999

C     Determine free space on DWORK 

      IWMAX0=IWMAX
      CALL ZFREE(3,IFREE)
      IF (IER.NE.0) GOTO 99999

C     Allocate the whole space on the heap for the calculation
C     of KCOL - we don't know its size beforehand. AP7X needs
C     a KCOL, a temporary KCOL1 and a temporary KINDX. So divide
C     the space on the heap in three different portions.

      NA=IFREE/3-3
      CALL ZNEW(NA,-3,LCOL,'KCOL  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LCOL1,'KCOL1 ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LIND,'KINDX ')
      IF (IER.NE.0) GOTO 99999

C     Call AP7 to create KCOL

      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))

      CALL AP7X(KWORK(L(LCOL)),KWORK(L(LCOL1)),
     *         KWORK(L(LIND)),KWORK(L(LLD)),
     *         NA,NEQ,TRIA,ELE,ISYMM, KWORK(KVERT),KWORK(KMID))
      IF (IER.NE.0) GOTO 99999

C     KCOL is finished. NA is the size of the matrix.
C     KCOL1 and KINDX we don't need anymore. Release the space
C     we don't need.

      CALL ZDISP(NA,LIND,'KINDX ')
      CALL ZDISP(NA,LCOL1,'KCOL1 ')
      CALL ZDISP(NA,LCOL,'KCOL  ')

      IWMAX=MAX(IWORK,IWMAX0)

      CALL ZDISP(0,LIND,'KINDX ')
      CALL ZDISP(0,LCOL1,'KCOL1 ')
      
99999 END
      
************************************************************************
* Generate matrix 7/8 structure
*
* Extended calling convention; Works for all types of Finite Elements
*
* This routine creates the column structure of a matrix in format 7.
* It expects an array KCOL which is large enough to receive the
* full column structure. Two more arrays KCOL1 and KINDX of the same
* size are needed as auxiliary arrays.
* The routine returns the actual size of the matrix in NA.
*
* In:
*   KCOL1  : array [1..NA] of integer
*            Auxiliary array
*   KINDX  : array [1..NA] of integer
*            Auxiliary array
*   NA     : Maximum length of KCOL
*   NEQ    : Number of equations in teh matrix
*   TRIA  : array [1..SZTRIA] of integer
*           Triangulation structure of the underlying mesh
*   ELE   : Finite Element callback routine.
*           Can be triangular or quadrilateral, but must fit to
*           the triangulation in TRIA!
*   ISYMM  : =0: don't respect symmetry, create matrix in structure 7
*            =1: matrix is symmetric; create matrix in structure 8
*   KVERT,
*   KMID   : Arrays describing the underlying triangulation; must
*            be the arrays given by TRIA, but passed separately
*            for performance reasons.
*
* Out:
*   NA     : Actual size of the matrix
*   KCOL   : array [1..*] of integer
*            Column structure of the matrix; the first NA entries are
*            filled with data
*
* Out (COMMON blocks):
*   IER    : =-118: Not enough space in KCOL
************************************************************************
      
      SUBROUTINE AP7X(KCOL,KCOL1,KINDX,KLD,NA,NEQ,TRIA,ELE,ISYMM,
     *                KVERT,KMID)
      
      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'stria.inc'
      
C     parameters 

      INTEGER KCOL(*),KCOL1(*),KINDX(*),KLD(*),NA,NEQ,ISYMM
      INTEGER KVERT(*),KMID(*)
      
      EXTERNAL ELE
      
C     local variables

      INTEGER IFREE,IEQ,IEL,IDOFE,JDOFE,IELTYP,IDOFE1,IROW
      INTEGER IPOS,ICOL,JCOL,IHELP,TRIA(SZTRIA)
      LOGICAL BSYMM,BSORT
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      
      INTEGER NDFL
      EXTERNAL NDFL

      SUB='AP7'
      IF (ICHECK.GE.997) CALL OTRC('AP7   ','12/11/89')

      IER=0

C     We expect at least one entry per line (the diagonal entry),
C     so NEQ entries are always necessary.

      IF (NA.LT.NEQ) THEN
      
C       Oops, KCOL is definitely too small for us!

        WRITE (CPARAM,'(I15)') 1
        CALL WERR(-118,'AP7   ')
        GOTO 99999
        
      ENDIF
      
C     How much memory have we left?
      
      IFREE=NA-NEQ
      NA=NEQ

C     Initialise KINDX and KCOL1.
C     Initially, we have only diagonal elements in our matrix.
C
C     The basic idea behind the building of the matrix is a linked
C     list of column numbers in each row!
C     We collect all upcoming columns in the whole matrix in the
C     array KCOL1, i.e. each new entry is attached to that
C     (-> the array is resorted to KCOL at the end).
C     KINDX points for every entry in the matrix to the position
C     inside of KCOL1 of the next entry in the line.
C
C     At the beginning, we only have the diagonal entries in the
C     matrix, which form the "head" of this collection of linked
C     lists. Initialise KCOL1(IEQ)=IEQ, so we know: Line IEQ
C     of the matrix starts with column KCOL1(IEQ), which is the
C     diagonal entry IEQ itself.
C     KINDX(IEQ) is set to 0 to indicate that there is no following
C     element in each line, i.e. we only have diagonal elements.
C
C     Example: We want to add entry (1,3) to the matrix. Then
C     we enlarge the lists as follows:
C     - "2" (the column number" is added to the end of KCOL1, i.e.
C       KCOL1(NEQ+1) = 3
C     - We add a "follower" for the diagonal element by setting
C       KIND(1) = NEQ+1. Furthermore we set KIND(NEQ+1)=0
C     So we have a linked list of matrix entries for the rows:
C       KCOL1:    1   2   3   ...   NEQ     3
C       KINDX:  NEQ+1 0   0   ...    0      0
C                 |                        /:\
C                 +-------------------------|
C     i.e. row 1 can be computed as:
C       KCOL1(1)            (=1),
C       KCOL1(IND(1))       (=3),
C       KCOL1(IND(IND(1))   -> not defined, as IND(IND(1))=0, 
C                              line finished

      DO IEQ=1,NEQ
        KINDX(IEQ)=0
        KCOL1(IEQ)=IEQ
      END DO

C     Call EA00 to get the element type into IELTYP
      
      CALL EA00(ELE,ELE,TRIA(ONVE),IELTYP)

C     Get the number of local DOF's from NDFL:      

      IDFL=NDFL(IELTYP)
      IF (IER.NE.0) GOTO 99999

C     Symmetric structure-8 matrix?

      BSYMM=ISYMM.EQ.1
      
C     For the loop through the DOF's later, set IDOFE=1.
C     For symmetric matrices (structure-8), this will be
C     changed to build only the upper triangular part...
      
      IDOFE1=1
      
C     Loop over elements

      DO IEL=1,TRIA(ONEL)
      
C       The outstanding feature with finite elements is: A basis
C       function for a DOF on one element has common support only
C       with the DOF's on the same element! E.g. for Q1:
C
C              #. . .#. . .#. . .#
C              .     .     .     .
C              .  *  .  *  .  *  .
C              #-----O-----O. . .#
C              |     |     |     .
C              |     | IEL |  *  .
C              #-----X-----O. . .#
C              |     |     |     .
C              |     |     |  *  .
C              #-----#-----#. . .#
C              
C       --> On element IEL, the basis function at "X" only interacts
C           with the basis functions in "O". Elements in the 
C           neighbourhood ("*") have no support, therefore we only have
C           to collect all "O" DOF's.

C       Call NDFGLX to get the global DOF's on our current
C       element (the "X" and all "O"'s). 
C       We don't need the local DOF's, so by setting IPAR=0,
C       the call will only fill KDFG.

        CALL NDFGLX(TRIA,IEL,0,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.NE.0) GOTO 99999

C       Loop over the IDFL DOF's on our current element IEL.

        DO IDOFE=1,IDFL
        
C         The DOF IDOFE is now our "X".
C         This global DOF gives us the row we have to build.
        
          IROW=KDFG(IDOFE)

C         We now want to loop through all the other DOF's on
C         element IEL (the "O"'s). For a standard structure-7
C         matrix, this means doing an inner loop JDOFE=1,IDFL.
C         For a symmetric structure-8 matrix, we only build
C         the upper triangular part of the matrix; this is realised
C         by looping JDOFE=IDOFE,IDFL.
C         So in the symmetric case, set IDOFE1 from 1 to IDOFE
C         to realise this.
          
          IF (BSYMM) THEN
          
C           In case our current DOF is the last one, there's
C           nothing to calculate in the symmetric case (last row,
C           last column!). So skip further calculation in this case.
          
            IF (IROW.EQ.NEQ) GOTO 110
            
            IDOFE1=IDOFE
            
          END IF

C         Loop over the other DOF's on the current element (the "O"'s).
C         For standard matrices, this loops 1..IDFL, for symmetric
C         matrices, this loops IDOFE..IDFL!

          DO JDOFE=IDOFE1,IDFL

C           Don't do anything if our "O" is our "X" - this is the
C           diagonal element which is always in the matrix.

            IF (IDOFE.NE.JDOFE) THEN
            
C             Get the global DOF - our "O". This gives the column number
C             in the matrix where an entry occurs in row IROW (the line of 
C             the current global DOF "X").
            
              JCOL=KDFG(JDOFE)

C             This JCOL has to be inserted into line IROW.
C             But first check, whether the element is already in that line,
C             i.e. whether element (IROW,JCOL) is already in the matrix.
C             This may happen because of an earlier loop in a neighbour
C             element (imagine Q1, where there are two vertices on an edge
C             and the edge is shared between two elements)...
C
C             We start walking through the linked list of row IROW to 
C             look for column JCOL. IPOS is the position in the KCOL1 
C             array of the column in row IROW we want to test. 
C             KINDX(IPOS) is =0 if we reach the end of the linked list.
C
C             Start searching at the "head" of the list, which is the
C             diagonal element at KCOL1(IROW).

              IPOS=IROW
              
C             Loop through the elements in the list. 
              
              DO WHILE (KINDX(IPOS).NE.0)

C               If we find the column JCOL, cancel the element insertion

                IF (KCOL1(IPOS).EQ.JCOL) GOTO 120
                
C               Otherwise proceed to the next element in the list
                
                IPOS=KINDX(IPOS)
                
              END DO
              
C             Element JCOL was not found in the line - insert it!
C             Is there enough space?

              IF (IFREE.LE.0) THEN
                WRITE (CPARAM,'(I15)') IEL
                CALL WERR(-118,'AP7   ')
                GOTO 99999
              END IF

C             There is enough space. Increase NA, which is the actual
C             length of the matrix. 

              NA=NA+1
              IFREE=IFREE-1
              
C             Append JCOL to KCOL1
              
              KCOL1(NA)=JCOL
              
C             Let KINDX of the last element of the row point to our
C             new element. The new element is now the last in the row.
              
              KINDX(IPOS)=NA
              KINDX(NA)=0
             
            END IF ! IDOFE <> JDOFE
 
120       END DO ! JDOFE

110     END DO ! IDOFE

      END DO ! IEL

C     Ok, KCOL1 is built. Now build KCOL by collecting the
C     entries in the linear lists of each row.
C
C     Set back NA to 0 at first.

      NA=0
      
C     Loop through all of the NEQ linear lists:
      
      DO IEQ=1,NEQ
      
C       Add the diagonal entry (head of the list) to KCOL:
      
        NA=NA+1
        KCOL(NA)=IEQ

C       Set KLD appropriately:

        KLD(IEQ)=NA
        
C       We are at the head of the list, now we have to walk
C       through it to append the entries to KCOL.
        
        IPOS=IEQ
        
        DO WHILE (KINDX(IPOS).NE.0)
        
C         Get the position of the next entry in KCOL1:

          IPOS=KINDX(IPOS)
          
C         Add the column number to the row in KCOL:

          NA=NA+1
          KCOL(NA)=KCOL1(IPOS)
        
        END DO ! KINDX(IPOS) <> 0

      END DO ! IEQ
      
C     Append the final entry to KLD:
      
      KLD(NEQ+1)=NA+1

C     Sort off-diagonal entries on KCOL separately for each row.
C     This is a small bubble-sort...
C
C     Loop through all rows:

      DO IEQ=1,NEQ

C       BSORT is set to false as soon as we find a non-sorted 
C       element - this is a crippled REPEAT-UNTIL loop realised
C       with GOTO.

300     BSORT=.TRUE.

C       REPEAT
C
C         Loop through the line except for the diagonal and the
C         last element

          DO ICOL=KLD(IEQ)+1,KLD(IEQ+1)-2
          
C           If the next element is larger...
          
            IF (KCOL(ICOL).GT.KCOL(ICOL+1)) THEN
            
C             Change position of the current and next element
            
              IHELP=KCOL(ICOL)
              KCOL(ICOL)=KCOL(ICOL+1)
              KCOL(ICOL+1)=IHELP
              
C             And repeat the sorting of that line
              
              BSORT=.FALSE.
              
            ENDIF
            
          END DO ! ICOL
          
C       UNTIL sorted
          
        IF (.NOT.BSORT) GOTO 300

      END DO ! IEQ

99999 END

************************************************************************
* Generate matrix 9 structure, wrapper routine
*
* Extended calling convention; Works for all types of Finite Elements
*
* This routine calls AP9X to generate a the structure of a matrix.
*
* In:
*   TRIA  : array [1..SZTRIA] of integer
*           Triangulation structure of the underlying mesh
*   ELE1  : Finite Element callback routine of the matrix block
*           which is to generate (e.g. velocity).
*           Can be triangular or quadrilateral, but must fit to
*           the triangulation in TRIA!
*   ELE2  : Finite Element callback routine for elements that
*           are coupled to the matrix block which is to generate
*           (e.g. pressure).
*           Can be triangular or quadrilateral, but must fit to
*           the triangulation in TRIA!
*   ISYMM : =0: don't respect symmetry, create matrix in structure 7
*           =1: matrix is symmetric; create matrix in structure 8
* Out:
*   LCOL  : Handle to KCOL array of the matrix
*   LLD   : Handle to KLD array of the matrix
*   NA    : Number of entries in the matrix
*   NEQ   : Dimension/Number of rows in the matrix
************************************************************************

      SUBROUTINE XAP9X(LCOL,LLD,NA,NEQ,TRIA,ELE1,ELE2)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'stria.inc'
      
C     parameters

      INTEGER LCOL,LLD,NA,NEQ,TRIA(SZTRIA)
      
      EXTERNAL ELE1,ELE2
      
C     local variables
      
      INTEGER IWMAX0,IFREE,LCOL1,LIND,ITYPE1
      INTEGER KVERT,KMID
      
      INTEGER NDFGX
      EXTERNAL NDFGX

      SUB='XAP9'
      IF (ICHECK.GE.997) CALL OTRC('XAP9  ','12/20/89')
C     Call the element wrapper EA00 to determine the
C     element type (switches between triangular and quadrilateral
C     element depending on NVE)

      CALL EA00(ELE1,ELE1,TRIA(ONVE),ITYPE1)

C     Get the number of equations via NDFGX
      
      NEQ=NDFGX(ITYPE1,TRIA)
      IF (IER.NE.0) GOTO 99999

C     Allocate KLD on DWORK 

      CALL ZNEW(NEQ+1,-3,LLD,'KLD   ')
      IF (IER.NE.0) GOTO 99999

C     Determine free space on DWORK 

      IWMAX0=IWMAX
      CALL ZFREE(3,IFREE)
      IF (IER.NE.0) GOTO 99999

C     Allocate the whole space on the heap for the calculation
C     of KCOL - we don't know its size beforehand. AP7X needs
C     a KCOL, a temporary KCOL1 and a temporary KINDX. So divide
C     the space on the heap in three different portions.

      NA=IFREE/3-3
      CALL ZNEW(NA,-3,LCOL,'KCOL  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LCOL1,'KCOL1 ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LIND,'KINDX ')
      IF (IER.NE.0) GOTO 99999

C     Call AP7 to create KCOL

      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))

      CALL AP9X(KWORK(L(LCOL)),KWORK(L(LCOL1)),
     *         KWORK(L(LIND)),KWORK(L(LLD)),
     *         NA,NEQ,TRIA,ELE1,ELE2,KWORK(KVERT),KWORK(KMID))
      IF (IER.NE.0) GOTO 99999

C     KCOL is finished. NA is the size of the matrix.
C     KCOL1 and KINDX we don't need anymore. Release the space
C     we don't need.

      CALL ZDISP(NA,LIND,'KINDX ')
      CALL ZDISP(NA,LCOL1,'KCOL1 ')
      CALL ZDISP(NA,LCOL,'KCOL  ')

      IWMAX=MAX(IWORK,IWMAX0)

      CALL ZDISP(0,LIND,'KINDX ')
      CALL ZDISP(0,LCOL1,'KCOL1 ')

99999 END

************************************************************************
* Generate matrix 9 structure
*
* Extended calling convention; Works for all types of Finite Elements
*
* This routine creates the column structure of a matrix in format 9.
* It expects an array KCOL which is large enough to receive the
* full column structure. Two more arrays KCOL1 and KINDX of the same
* size are needed as auxiliary arrays.
* The routine returns the actual size of the matrix in NA.
*
* In:
*   KCOL1  : array [1..NA] of integer
*            Auxiliary array
*   KINDX  : array [1..NA] of integer
*            Auxiliary array
*   NA     : Maximum length of KCOL
*   NEQ    : Number of equations in teh matrix
*   TRIA  : array [1..SZTRIA] of integer
*           Triangulation structure of the underlying mesh
*   ELE1  : Finite Element callback routine of the matrix block
*           which is to generate (e.g. velocity).
*           Can be triangular or quadrilateral, but must fit to
*           the triangulation in TRIA!
*   ELE2  : Finite Element callback routine for elements that
*           are coupled to the matrix block which is to generate
*           (e.g. pressure).
*           Can be triangular or quadrilateral, but must fit to
*           the triangulation in TRIA!
*   ISYMM  : =0: don't respect symmetry, create matrix in structure 7
*            =1: matrix is symmetric; create matrix in structure 8
*   KVERT,
*   KMID   : Arrays describing the underlying triangulation; must
*            be the arrays given by TRIA, but passed separately
*            for performance reasons.
*
* Out:
*   NA     : Actual size of the matrix
*   KCOL   : array [1..*] of integer
*            Column structure of the matrix; the first NA entries are
*            filled with data
*
* Out (COMMON blocks):
*   IER    : =-118: Not enough space in KCOL
************************************************************************

      SUBROUTINE AP9X(KCOL,KCOL1,KINDX,KLD,NA,NEQ,TRIA,ELE1,ELE2,
     *                KVERT,KMID)
      
      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'stria.inc'
      
C     parameters 

      INTEGER KCOL(*),KCOL1(*),KINDX(*),KLD(*),NA,NEQ
      INTEGER KVERT(*),KMID(*),TRIA(SZTRIA)
      
      EXTERNAL ELE1,ELE2
      
C     local variables

      INTEGER IFREE,IEQ,IEL,JDOFE,IROW
      INTEGER IPOS,ICOL,JCOL,IHELP,ITYP1(3),NELE,I,JDOFP
      LOGICAL BSORT
      INTEGER KDFG1(NNBAS,3),KDFL1(NNBAS,3),IDFL1(3)
      
      INTEGER NDFL
      EXTERNAL NDFL

      SUB='AP9'
      IF (ICHECK.GE.997) CALL OTRC('AP9   ','12/20/89')

      IER=0

C     We expect at least one entry per line,
C     so NEQ entries are always necessary.

      IF (NA.LT.NEQ) THEN
      
C       Oops, KCOL is definitely too small for us!

        WRITE (CPARAM,'(I15)') 1
        CALL WERR(-118,'AP9   ')
        GOTO 99999
      END IF
      
C     How much memory have we left?
      
      IFREE=NA-NEQ
      NA=NEQ

C     Initialise KINDX and KCOL1.
C     Initially, we have no elements in our matrix, but we expect
C     at least one element per line.
C
C     The basic idea behind the building of the matrix is a linked
C     list of column numbers in each row!
C     We collect all upcoming columns in the whole matrix in the
C     array KCOL1, i.e. each new entry is attached to that
C     (-> the array is resorted to KCOL at the end).
C     KINDX points for every entry in the matrix to the position
C     inside of KCOL1 of the next entry in the line.
C
C     At the beginning, we only have the no entries in the
C     matrix. We initialise the "head" of this collection of linked
C     lists with 0 to indicate this. 
C     KINDX(IEQ) is set to 0 to indicate that there is no following
C     element in each line, i.e. we only have diagonal elements.
C     Later, we fill KCOL1(1..NEQ) with the first column number in
C     each row. When more entries appear in a row, they are appended
C     to KCOL at position NEQ+1..*. KINDX keeps track of the entries
C     in each row by storing (starting from the head in KCOL1(IEQ))
C     the positions of the corresponding next entry inside of KCOL1
C     in each row - so the column numbers in each row are to be
C     found in KCOL1 at positions IEQ,KINDX(IEQ),KINDX(KINDX(IEQ)),...
C     
C
C     Example: We want to add entry (1,3) to the matrix. Then
C     we enlarge the lists as follows:
C     - "2" (the column number" is added to the end of KCOL1, i.e.
C       KCOL1(NEQ+1) = 3
C     - We add a "follower" for the diagonal element by setting
C       KIND(1) = NEQ+1. Furthermore we set KIND(NEQ+1)=0
C     So we have a linked list of matrix entries for the rows:
C       KCOL1:    1   0   0   ...   NEQ     3
C       KINDX:  NEQ+1 0   0   ...    0      0
C                 |                        /:\
C                 +-------------------------|
C     i.e. row 1 can be computed as:
C       KCOL1(1)            (=1),
C       KCOL1(IND(1))       (=3),
C       KCOL1(IND(IND(1))   -> not defined, as IND(IND(1))=0, 
C                              line finished

      DO IEQ=1,NEQ
        KINDX(IEQ)=0
        KCOL1(IEQ)=0
      END DO

C     Call EA00 to get the element type of the two elements
C     into ITYP:
      
      CALL EA00(ELE1,ELE1,TRIA(ONVE),ITYP1(1))
      CALL EA00(ELE2,ELE2,TRIA(ONVE),ITYP1(2))
      
C     Check if it's the same element (e.g. in a Q1/Q1 discretisation):
      
      IF (ITYP1(1).EQ.ITYP1(2)) THEN
       NELE=1
      ELSE
       NELE=2
      ENDIF

C     Get the number of local DOF's from NDFL:      

      DO I=1,NELE
        IDFL1(I)=NDFL(ITYP1(I))
        IF (IER.NE.0) GOTO 99999
      END DO

C     Loop over elements

      DO IEL=1,TRIA(ONEL)

C       The outstanding feature with finite elements is: A basis
C       function for a DOF on one element has common support only
C       with the DOF's on the same element! E.g. for Q1:
C
C              #. . .#. . .#. . .#
C              .     .     .     .
C              .  *  .  *  .  *  .
C              .     .     .     .
C              #-----#-----#. . .#
C              |     | O O |     .
C              |     | IEL |  *  .
C              |     | O O |     .
C              #-----X-----#. . .#
C              |     |     |     .
C              |     |     |  *  .
C              |     |     |     .
C              #-----#-----#. . .#
C              
C       --> On element IEL, the basis function at "X" of element ELE1
C           only interacts with the basis functions in "O" of element 
C           ELE2. Elements in the neighbourhood ("*") have no support, 
C           therefore we only have to collect all "O" DOF's.
C
C       Call NDFGLX to get the global DOF's of both ELE's
C       on our current element (the "X" and all "O"'s). 
C       We don't need the local DOF's, so by setting IPAR=0,
C       the call will only fill KDFG.

        DO I=1,NELE
          CALL NDFGLX(TRIA,IEL,0,ITYP1(I),KVERT,KMID,
     *                KDFG1(1,I),KDFL1(1,I))
          IF (IER.NE.0) GOTO 99999
        END DO

C       Loop over the IDFL DOF's on our current element IEL
C       of our "main" element IEL1 (e.g. velocity DOF's):

        DO JDOFE=1,IDFL1(1)
        
C         The DOF IDOFE is now our "X".
C         This global DOF gives us the row we have to build.
        
          IROW=KDFG1(JDOFE,1)

C         We now loop through all the other DOF's on
C         element IEL, corresponding to the "coupled" element ELE2
C         (the "O"'s). 

          DO JDOFP=1,IDFL1(NELE)
          
C           The global DOF (one of the "O"'s) gives us the column
C           number we have to add to the current row IROW:
          
            JCOL=KDFG1(JDOFP,NELE)
            
C           This JCOL has to be inserted into line IROW.
C           But first check, whether the element is already in that line,
C           i.e. whether element (IROW,JCOL) is already in the matrix.
C           This may happen because of an earlier loop in a neighbour
C           element (imagine Q1, where there are two vertices on an edge
C           and the edge is shared between two elements)...
C
C           Is the list empty?
C
            IF (KCOL1(IROW).EQ.0) THEN

C             Yes, row IROW is empty at the moment. Add the column as
C             head of the list of row IROW.

              KCOL1(IROW)=JCOL

            ELSE
            
C             No, the list is not empty, we have a "head".
C
C             We start walking through the linked list of row IROW to 
C             look for column JCOL. IPOS is the position in the KCOL1 
C             array of the column in row IROW we want to test. 
C             KINDX(IPOS) is =0 if we reach the end of the linked list.
C
C             Start searching at the "head" of the list, which is the
C             diagonal element at KCOL1(IROW).

              IPOS=IROW
              
C             Loop through the elements in the list. 
              
              DO WHILE (KINDX(IPOS).NE.0)

C               If we find the column JCOL, cancel the element insertion

                IF (KCOL1(IPOS).EQ.JCOL) GOTO 120
                
C               Otherwise proceed to the next element in the list
                
                IPOS=KINDX(IPOS)
                
              END DO
              
C             Element JCOL was not found in the line - insert it!
C             Is there enough space?
              
              IF (IFREE.LE.0) THEN
                WRITE (CPARAM,'(I15)') IEL
                CALL WERR(-118,'AP9   ')
                GOTO 99999
              END IF

C             There is enough space. Increase NA, which is the actual
C             length of the matrix. 

              NA=NA+1
              IFREE=IFREE-1
              
C             Append JCOL to KCOL1
              
              KCOL1(NA)=JCOL
              
C             Let KINDX of the last element of the row point to our
C             new element. The new element is now the last in the row.
              
              KINDX(IPOS)=NA
              KINDX(NA)=0
              
            END IF ! KCOL1(IROW)=0
              
120       END DO ! JDOFE

        END DO ! IDOFE
        
      END DO ! IEL

C     Ok, KCOL1 is built. Now build KCOL by collecting the
C     entries in the linear lists of each row.
C
C     Set back NA to 0 at first.

      NA=0
      
C     Loop through all of the NEQ linear lists:
      
      DO IEQ=1,NEQ
      
C       Add the diagonal entry (head of the list) to KCOL:
      
        NA=NA+1
        KCOL(NA)=KCOL1(IEQ)

C       Set KLD appropriately:

        KLD(IEQ)=NA
        
C       We are at the head of the list, now we have to walk
C       through it to append the entries to KCOL.
        
        IPOS=IEQ
        
        DO WHILE (KINDX(IPOS).NE.0)
        
C         Get the position of the next entry in KCOL1:

          IPOS=KINDX(IPOS)
          
C         Add the column number to the row in KCOL:

          NA=NA+1
          KCOL(NA)=KCOL1(IPOS)
        
        END DO ! KINDX(IPOS) <> 0

      END DO ! IEQ
      
C     Append the final entry to KLD:
      
      KLD(NEQ+1)=NA+1

C     Sort off-diagonal entries on KCOL separately for each row.
C     This is a small bubble-sort...
C
C     Loop through all rows:

      DO IEQ=1,NEQ

C       BSORT is set to false as soon as we find a non-sorted 
C       element - this is a crippled REPEAT-UNTIL loop realised
C       with GOTO.

300     BSORT=.TRUE.

C       REPEAT
C
C         Loop through the line except for the diagonal and the
C         last element

          DO ICOL=KLD(IEQ),KLD(IEQ+1)-2
          
C           If the next element is larger...
          
            IF (KCOL(ICOL).GT.KCOL(ICOL+1)) THEN
            
C             Change position of the current and next element
            
              IHELP=KCOL(ICOL)
              KCOL(ICOL)=KCOL(ICOL+1)
              KCOL(ICOL+1)=IHELP
              
C             And repeat the sorting of that line
              
              BSORT=.FALSE.
              
            END IF
            
          END DO ! ICOL
          
C       UNTIL sorted
          
        IF (.NOT.BSORT) GOTO 300

      END DO ! IEQ

99999 END

