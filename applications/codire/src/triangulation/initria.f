***********************************************************************
* This file contains maintanance routines for the triangulation
* structures. It allowes to transfer data between the current 
* /TRIAx/ common blocks and a STRIA structure-array (see stria.inc).
* Furthermore it contains routines to read a triangulation or a
* sequence of triangulations from a file and to call refinement
* routines to refine a given coarse grid.
***********************************************************************

***********************************************************************
* Common-Block to STRIA
*
* Copies the content of the TRIAx COMMON-blocks to a TRIA 
* structure array.
*
* In:
*  The information in the COMMON blocks
* Out:
*  TRIA   - array [1..SZTRIA] of integer
*           the STRIA structure-array
***********************************************************************

      SUBROUTINE C2TRIA (TRIA)
      
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INTEGER TRIA(SZTRIA)
      
      TRIA(ONEL   ) = NEL   
      TRIA(ONVT   ) = NVT   
      TRIA(ONMT   ) = NMT   
      TRIA(ONVE   ) = NVE   
      TRIA(ONVEL  ) = NVEL  
      TRIA(ONBCT  ) = NBCT  
      TRIA(ONVBD  ) = NVBD  

      TRIA(OLCORVG) = LCORVG
      TRIA(OLCORMG) = LCORMG
      TRIA(OLVERT ) = LVERT 
      TRIA(OLMID  ) = LMID  
      TRIA(OLADJ  ) = LADJ  
      TRIA(OLVEL  ) = LVEL  
      TRIA(OLMEL  ) = LMEL  
      TRIA(OLNPR  ) = LNPR  
      TRIA(OLMM   ) = LMM   
      TRIA(OLVBD  ) = LVBD  
      TRIA(OLEBD  ) = LEBD  
      TRIA(OLBCT  ) = LBCT  
      TRIA(OLVBDP ) = LVBDP 
      TRIA(OLMBDP ) = LMBDP 
      
      END
      
***********************************************************************
* STRIA to Common-Block
*
* Copies the content of the TRIA-array to the TRIAx COMMON-blocks.
*
* In: 
*  TRIA   - array [1..SZTRIA] of integer
*           the STRIA structure-array
* Out:
*  The information in the COMMON-blocks.
***********************************************************************

      SUBROUTINE TRIA2C (TRIA)
      
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INTEGER TRIA(SZTRIA)
      
      NEL    = TRIA(ONEL   )
      NVT    = TRIA(ONVT   )
      NMT    = TRIA(ONMT   )
      NVE    = TRIA(ONVE   )
      NVEL   = TRIA(ONVEL  )
      NBCT   = TRIA(ONBCT  )
      NVBD   = TRIA(ONVBD  )
      
      LCORVG = TRIA(OLCORVG)
      LCORMG = TRIA(OLCORMG)
      LVERT  = TRIA(OLVERT )
      LMID   = TRIA(OLMID  )
      LADJ   = TRIA(OLADJ  )
      LVEL   = TRIA(OLVEL  )
      LMEL   = TRIA(OLMEL  )
      LNPR   = TRIA(OLNPR  )
      LMM    = TRIA(OLMM   )
      LVBD   = TRIA(OLVBD  )
      LEBD   = TRIA(OLEBD  )
      LBCT   = TRIA(OLBCT  )
      LVBDP  = TRIA(OLVBDP )
      LMBDP  = TRIA(OLMBDP )
      
      END
      
***********************************************************************
* Convert quad coarse grid to triangular coarse grid
*
* This routine converts a quadrilateral coarse grid to a triangular 
* coarse grid be dividing each quadrilateral element into two triangles.
*
* The routine should be called directly after XORSC to convert
* the coarse grid. It is not designed to convert a fully refined
* mesh with all information!
*
* In:
*   Data from the TRIAA/TRIAD COMMON blocks
*
* Out:
*   The TRIAA/TRIAD COMMON blocks are changed to reflect a triangular
*   coarse grid.
***********************************************************************

      SUBROUTINE Q2TRIC
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
C     local variables

      INTEGER LVERT2
      
C     Don't do anything if the mesh is already triangular

      IF (NVE.NE.4) RETURN
      
C     Allocate a new KVERT array. As every element is subdivided
C     into two, we exactly double the number of elements!
      
      CALL ZNEW(2*NEL*NNVE,-3,LVERT2,'KVERT2')
      
C     Creathe the new KVERT:

      CALL Q2TAX1 (NEL,KWORK(L(LVERT)),KWORK(L(LVERT2)))
      
C     Release the old KVERT arrays from heap

      CALL ZDISP(0,LVERT,'KVERT1')
      
C     Transfer the new handles to the triangulation.

      LVERT = LVERT2
      
C     Indicate that we now have a triangular mesh.

      NVE = 3
      NEL = 2*NEL
      
C     Release old KADJ and create a new one

      CALL ZDISP(0,LADJ,'KADJ  ')
      CALL XS2A
      
      END
      
***********************************************************************
* Generate standard triangulation
*
* This routine allowes to create a triangulation structure for a new
* triangulation by
*  - duplication of another triangulation
*  - reading in the triangulation from a file
*  - refining a given triangulation.
*
* In:
*   IMETH  - Method how to create the triangulation
*            =0: Create the new triangulation by refinement of an
*                old triangulation TRISRC. the old triangulation
*                will be lost.
*            =1: Create the triangulation by duplication+refinement
*                of an old one.
*                The triangulation structure TRISRC is used as
*                prototype for the new triangulation. The routine
*                uses the duplication flags IFLDUP as documented
*                in TRIDUB for the duplication. This allowes the new
*                triangulation to share information with the old one
*                (like the vertiex coordinates e.g.).
*                CMESH is ignored.
*            =2: Read the triangulation from a .tri coarse grid file.
*                IFLDUP and TRISRC is ignored.
*
*                If MFILE=0, the triangulation is read from the 
*                .tri coarse grid triangulation file CMESH. 
*                Unit no. 65 is used for reading and closed afterwards.
*
*                If MFILE<>0, the triangulation is read from the
*                input channel MFILE. If the input channel is closed,
*                it will be opened with filename CMESH; otherwise, 
*                CMESH is ignored. The file will not be closed after
*                reading.
*
*                The triangulation is expected to be stored in 
*                DEVISOR FEAT2D format. IFMT is ignored.
*
*            =3: Read the triangulation from a generated triangulation
*                file. IFMT decides on whether the input is in formatted
*                or unformatted form.
*                IFLDUP and TRISRC is ignored.
*
*                If MFILE=0, the triangulation is read from the 
*                triangulation file CMESH. Unit no. 65 is used for
*                reading.
*
*                If MFILE<>0, the triangulation is read from the
*                input channel MFILE. If the input channel is closed,
*                it will be opened with filename CMESH; otherwise, 
*                CMESH is ignored. The file will not be closed after
*                reading.
*
*                The triangulation is expected to be stored in 
*                DEVISOR FEAT2D format. IFMT is ignored.
*       
*   INVE   - Expected number of vertices per element.
*            Only effective if IMESH=2.
*            = 0: Don't care.
*            <>0: If the number of vertices per element in a readed 
*                 coarse grid differ from INVE (e.g. if the coarse grid
*                 is a quad mesh with NVE=4 and there is INVE=3 given),
*                 the GENGRI routine tries to convert the grid before
*                 refinement into a grid with INVE vertices per element.
*                 If this does not work, the program is haltet.
*
*   IFMT   - Formatted input if IMETH=3.
*            =0: Use unformatted input when reading the triangulation
*                from a file
*            =1: Use formatted input when reading the triangulation
*                from a file
*               
*   IREF   - Number of refinements.
*            After generating the mesh by duplication or reading from
*            a file, the mesh is IREF times refined. The resulting mesh
*            is stored in TRIA.
*
*            Note that every information that is shared between
*            two triangulations (by setting the bits in IFLDUP
*            appropriately) is changed for both grids in the refinement
*            process!
*
*   IFGCOR - Post-Correction of the fine grid nodes after refinement.
*            Only respected if IREF>0.
*            =0: no correction
*            =1: After each refinement, the element midpoints are
*                recalculated by interpolation. This helps e.g. against
*                anisotropic cells when the geometry contains
*                curved boundaries like circles. Standard.
*
*   IFLDUP - If IMETH=2, this flag is used to define which information
*            is duplicated and which is shared between TRISRC and TRIA.
*            Note that any information shared between TRIA and TRISRC
*            gets invalid in TRISRC when refining TRIA (see below)!
*
* Out:
*   TRIA   - array [1..SZTRIA] of integer
*            New triangulation, generated by duplication/reading and
*            refinement.
*
* Remarks: 
* 1.) Variables in the TRIAx COMMON blocks are undefined
*  after leaving this routine.
* 2.) The parametrization must be initialized before calling this
*  routine.
* 3.) Generated triangulation structures can be deleted with TRIDEL.
*  TRIDEL will automatically take care about if it's allowed to
*  release arrays of the structure from the heap by maintaining
*  the duplication flags. 
* 4.) A little bit of care has to be taken by the caller when creating
*  TRIA from TRISRC by duplication + refinement (IMETH=1)! As long
*  as IDPFLG=0, using duplication+refinement is ok. But as soon as
*  IDPFLG<>0, information is shared between TRIA and TRISRC upon
*  duplication. When using refinement of TRIA afterwards, the grid
*  information of arrays shared between TRIA and TRISRC gets
*  lost in TRISRC. 
*  This might even mean that handles in TRISRC get invalid! Therefore,
*  the caller should use TRICSH to correct the handles in TRISRC of
*  shared information to those in TRIA!
*  Example: DCORVG is shared between all refinement levels.
*  Then the caller can create the grid by:
*
*    CALL GENTRI (IMETH=2, filename='anything.tri')  -> read coarse gr.
*    DO I=NLMIN+1,NLMAX                              -> refine to NLMAX
*      CALL GENTRI (TRIA(I-1),TRIA(I),IMETH=1,IDPFLG=1)
*      --> generates TRIA(I), whereas destroys LCORVG on level I-1!
*    END DO
*    DO I=NLMAX-1,NLMIN,-1             -> write LCORVG of lv. n+1 into
*      CALL TRICSH (TRIA(I),TRIA(I+1))    LCORVG of level n
*    END DO
*
*  Afterwards, the information is again shared correctly between all
*  levels.
***********************************************************************
      
      SUBROUTINE GENTRI (TRIA, IMETH, INVE, IREF, IFGCOR, 
     *                   IFLDUP, TRISRC, 
     *                   MFILE, IFMT, CMESH)
      
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
C     We use "old" refinement routines from the FEAT library and
C     have to access COMMON block variables for that purpose:

      INCLUDE 'ctria.inc'
      
C     parameters      
      
      INTEGER TRIA(SZTRIA), TRISRC (SZTRIA), IMETH, IREF, IFLDUP, IFMT
      INTEGER MFILE,IFGCOR,INVE
      CHARACTER CMESH*(*)
      
C     local variables

      INTEGER MF
      INTEGER IMID,IADJ,IVEL,IDISP,IBDP
      
C     externals

      EXTERNAL S2DI0,S2DB0
      EXTERNAL PARX,PARY,TMAX

C     When calculating structures, use the following settings:
C
C     Generate subdivisions directly.
C     IMID=2 to calculate midpoints,
C     IVEL=1 to calculate the elements meeting on a vertex.
C     Important for force calculations and grid deformation, 
C     if used.

      IMID=2
      IADJ=1
      IVEL=1
      IDISP=1
      IBDP=2

C     Should we read the triangulation?

      IF (IMETH.EQ.3) THEN

C       We should read data that was written out by XOWA.
C       Is an input channel given?

        IF (MFILE.EQ.0) THEN

C         Read the triangulation from the file and create it in the
C         COMMON blocks:

          MF = 65
          CALL XORS(MFILE,CMESH,IFMT)
          CLOSE (MF)
          IF (IER.NE.0) RETURN
          
        ELSE 
        
C         Read the triangulation from the file without starting at
C         the beginning of the file. This can be used to read a mesh
C         out of a sequence of meshes stored with XOWA into one file.

          CALL XORS(MFILE,CMESH,IFMT)
          
C         Don't close the file, it comes from the caller...
        
        END IF
        
C       The triangulation is now in the COMMON block.
C       Also make a copy to the structure array:

        CALL C2TRIA (TRIA)

      ELSE IF (IMETH.EQ.2) THEN

C       We should read data in DEVISOR form.
C       Is an input channel given?

        IF (MFILE.EQ.0) THEN

C         Read the triangulation from the file and create it in the
C         COMMON blocks:

          MF = 65
          CALL  XORSC (MF,CMESH)
          CLOSE (MF)
          IF (IER.NE.0) RETURN

C         Eventually try to convert the grid 

          IF ((INVE.NE.0).AND.(NVE.NE.INVE)) THEN
            IF (INVE.EQ.3) CALL Q2TRIC
            IF (NVE.NE.INVE) THEN
              WRITE (MTERM,'(A)') 
     *          'Fatal: Cannot convert grid to meet discretisation '//
     *          'requirements!'
              STOP
            END IF
          END IF

C         create+initialize all structures without refinement

          IF (NVE.EQ.3) THEN
            CALL XSA0X(0,IMID,IADJ,IVEL,IDISP,IBDP,
     *                 S2DI0,S2DB0,PARX,PARY,TMAX)
          ELSE
            CALL XSB0X(0,IMID,IADJ,IVEL,IDISP,IBDP,
     *                 S2DI0,S2DB0,PARX,PARY,TMAX)
          END IF
          
        ELSE 
        
C         Read the triangulation from the given input channel.

          CALL  XORSC (MFILE,CMESH)

C         Don't close the file, it comes from the caller...
C
C         Eventually try to convert the grid 

          IF ((INVE.NE.0).AND.(NVE.NE.INVE)) THEN
            IF (INVE.EQ.3) CALL Q2TRIC
            IF (NVE.NE.INVE) THEN
              WRITE (MTERM,'(A)') 
     *          'Fatal: Cannot convert grid to meet discretisation '//
     *          'requirements!'
              STOP
            END IF
          END IF

C         create+initialize all structures without refinement

          IF (NVE.EQ.3) THEN
            CALL XSA0X(0,IMID,IADJ,IVEL,IDISP,IBDP,
     *                 S2DI0,S2DB0,PARX,PARY,TMAX)
          ELSE
            CALL XSB0X(0,IMID,IADJ,IVEL,IDISP,IBDP,
     *                 S2DI0,S2DB0,PARX,PARY,TMAX)
          END IF
           
        END IF

C       The triangulation is now in the COMMON block.
C       Also make a copy to the structure array:

        CALL C2TRIA (TRIA)

      ELSE IF (IMETH.EQ.1) THEN
      
C       Create the mesh by duplication of an old one.

        CALL TRIDUP (TRISRC,TRIA,0,IFLDUP)

C       Save it to the COMMON block for possible refinement

        CALL TRIA2C (TRIA)
      
      ELSE
      
C       Don't duplicate anything, take the new mesh as the old

        CALL LCP3 (TRISRC,TRIA,SZTRIA)
        
C       Save it to the COMMON block for possible refinement

        CALL TRIA2C (TRIA)
      
      END IF
      
C     Ok, we have a mesh. Should we refine it a little bit?

      IF (IREF.GT.0) THEN
      
C       For the refinement process, any extended information is
C       useless and must be recalculated.
C       Dispose extended structures without converting anything:

        CALL DISETR (TRIA,0)
      
C       Call the refinement routine; re-initialize 

        IF (TRIA(ONVE).EQ.3) THEN
          CALL XSA0X(IREF,IMID,IADJ,IVEL,IDISP,IBDP,
     *               S2DI0,S2DB0,PARX,PARY,TMAX)
        ELSE
          CALL XSB0X(IREF,IMID,IADJ,IVEL,IDISP,IBDP,
     *               S2DI0,S2DB0,PARX,PARY,TMAX)
        END IF
        IF (IER.NE.0) RETURN
        
C       Fo a post-correction of the coordinates in the finer grid.
C       This is only allowed for quads!

        IF ((IFGCOR.NE.0).AND.(TRIA(ONVE).EQ.4)) THEN
          CALL AVEMPC(DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LADJ)),
     *                NEL)
        END IF
        
C       We have a refined grid in the COMMON block.
C        
C       Get the mesh back from the COMMON block

        CALL C2TRIA (TRIA)
      
      END IF

C     If it does not exist, build an extended structure:

      IF (TRIA(OTRIFL).EQ.0) THEN
        CALL GENETR (TRIA)
      END IF
      
      END

