************************************************************************
* This file contains auxiliary routines for the post processing
* for tracking falues of functionals on a solution.
************************************************************************

************************************************************************
* Initialize file output 
*
* This routines opens output files.
* The solver will write tracked Drag/Lift/Solution/... values
* into these files.
* The file handles are saved into the parameter block. After the
* solver has been finished, the caller should call DONISF to close
* the files.
* The caller can provide an index number IIDX which is added
* to the filename of the output files.
*
* In:
*   IPARAM : array [1..SZNSDI] of integer
*            TTrackingIParams parameter block
*   IIDX   : File index number. If > -1, this number is added to the
*            filenames of the output files.
*            Standard = -1
*
* Out:
*   IPARAM : Is initialized
************************************************************************

      SUBROUTINE INITRK (IPARAM,IIDX)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'stracking.inc'
      
      INCLUDE 'dstrings.inc'
      
      INTEGER IPARAM(SZTRKI),IIDX
      
      INTEGER LSTR
      CHARACTER FNAME*60

C     Open files for...

C     ... tracking of drag/lift with boundary integration

      IPARAM(OHBFB1) = 40

      LSTR = STNEWC (.FALSE.,'points/force_bdy_x')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHBFB1),FILE=FNAME)

      IPARAM(OHBFB2) = 41

      LSTR = STNEWC (.FALSE.,'points/force_bdy_y')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHBFB2),FILE=FNAME)

C     ... tracking of drag/lift with volume integration

      IPARAM(OHBFV1) = 42

      LSTR = STNEWC (.FALSE.,'points/force_fbdy_vol_x')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHBFV1),FILE=FNAME)
      
      IPARAM(OHBFV2) = 43

      LSTR = STNEWC (.FALSE.,'points/force_fbdy_vol_y')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHBFV2),FILE=FNAME)

C     ... tracking of drag/lift with line integration

      IPARAM(OHBFL1) = 44

      LSTR = STNEWC (.FALSE.,'points/force_fbdy_line_x')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHBFL1),FILE=FNAME)

      IPARAM(OHBFL2) = 45

      LSTR = STNEWC (.FALSE.,'points/force_fbdy_line_y')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHBFL2),FILE=FNAME)

C     ... tracking of velocity values

      IPARAM(OHVEL1) = 46

      LSTR = STNEWC (.FALSE.,'points/velocity_x')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHVEL1),FILE=FNAME)

      IPARAM(OHVEL2) = 47

      LSTR = STNEWC (.FALSE.,'points/velocity_y')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHVEL2),FILE=FNAME)

C     ... tracking of pressure values

      IPARAM(OHPRES) = 48

      LSTR = STNEWC (.FALSE.,'points/pressure')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHPRES),FILE=FNAME)

C     ... tracking of streamfunction values

      IPARAM(OHSTRM) = 49

      LSTR = STNEWC (.FALSE.,'points/streamfunction')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHSTRM),FILE=FNAME)

C     ... tracking of mean integral pressure values

      IPARAM(OHIPRS) = 50

      LSTR = STNEWC (.FALSE.,'points/integral_pressure')
      IF (IIDX.GT.-1) THEN
        CALL STCATC (LSTR,.FALSE.,'.')
        CALL STCATI (LSTR,IIDX,4,.TRUE.)
      END IF
      CALL STPUT (LSTR, FNAME)
      CALL STDIS (LSTR)

      OPEN (UNIT=IPARAM(OHIPRS),FILE=FNAME)

      END
      
************************************************************************
* Close file output for Nonstationary solver 
*
* This routine closes the output files that were opened by INIISF.
*
* In:
*   IPARAM : array [1..SZNSDI] of integer
*            TTrackingIParams parameter block
*
* Out:
*   All handles in IPARAM are invalid.
************************************************************************

      SUBROUTINE DONTRK (IPARAM)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'stracking.inc'
      
      INTEGER IPARAM(SZTRKI)

      CLOSE (IPARAM(OHIPRS))
      CLOSE (IPARAM(OHSTRM)) 
      CLOSE (IPARAM(OHPRES))
      CLOSE (IPARAM(OHVEL2)) 
      CLOSE (IPARAM(OHVEL1)) 

      CLOSE (IPARAM(OHBFB1))
      CLOSE (IPARAM(OHBFB2)) 

      CLOSE (IPARAM(OHBFV1))
      CLOSE (IPARAM(OHBFV2))

      CLOSE (IPARAM(OHBFL1)) 
      CLOSE (IPARAM(OHBFL2)) 

      END
      
************************************************************************
* Write tracked values to terminal
*
* This routine prints the results saved in IINFO/DINFO1/DINFO2 to the
* terminal and into a file (depending on the output level).
* The values are printed with a headline CNAME and ordered in
* the form:
*   [IINFO(1)] [DINFO1(1)] [DINFO2(1)]
*   [IINFO(2)] [DINFO1(2)] [DINFO2(2)]
*   [IINFO(3)] [DINFO1(3)] [DINFO2(3)]
* ...
* If IDIM=1, only DINFO1 is printed. If IDIM=2, DINFO1 and DINFO2
* are printed.
*
* In:
*   MSHOW  - Output level
*   MFILE  - Handle of a log file where to write output, too
*   CNT    - Number of data entries in IINFO/DINFO
*   IINFO  - array [1..CNT] of integer
*            Identifier for the entries in DINFO
*   DINFO1 - array [1..CNT] of double
*            First data set that should be written out
*   DINFO2 - array [1..CNT] of double
*            Second data set that should be written out
*   IDIM   - >=1: print DINFO1
*            >=2: print DINFO2
*   CNAME  - Headline to print
************************************************************************

      SUBROUTINE TV2TER (MSHOW,MFILE,CNT,IINFO,DINFO1,DINFO2,IDIM,
     *                   CNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      
      INTEGER CNT,IINFO(*),IDIM,MSHOW,MFILE
      DOUBLE PRECISION DINFO1(*),DINFO2(*)
      CHARACTER CNAME*(*)
      
      INTEGER I
      
      IF (MSHOW.GE.2) THEN

C       Write the headline

        WRITE (MTERM,'(A)') CNAME
        
C       Write the data

        IF (IDIM.LE.1) THEN
          DO I=1,CNT
            WRITE (MTERM,'(I10,E20.10)') IINFO(I),DINFO1(I)
          END DO
        ELSE IF (IDIM.LE.2) THEN
          DO I=1,CNT
            WRITE (MTERM,'(I10,2E20.10)') IINFO(I),DINFO1(I),DINFO2(I)
          END DO
        END IF
      
        WRITE (MTERM,*)

      END IF

      IF (MSHOW.GE.1) THEN

C       Write the headline

        WRITE (MFILE,'(A)') CNAME
        
C       Write the data

        IF (IDIM.LE.1) THEN
          DO I=1,CNT
            WRITE (MFILE,'(I10,E20.10)') IINFO(I),DINFO1(I)
          END DO
        ELSE IF (IDIM.LE.2) THEN
          DO I=1,CNT
            WRITE (MFILE,'(I10,2E20.10)') IINFO(I),DINFO1(I),DINFO2(I)
          END DO
        END IF
      
        WRITE (MFILE,*)

      END IF

      END

************************************************************************
* Write tracked values to file
*
* This writes the first CNT values in DINFO to the file with
* the file handle H in the format
*    [TIMENS] [DINFO(1)] [DINFO(2)] [DINFO(3)] ...
* or
*    [INUM] [TIMENS] [DINFO(1)] [DINFO(2)] [DINFO(3)] ...
* depending on whether INUM=-1 or not.
*
* If BHEAD=true, a headline will be written to the file.
* This contains the type of the values in DINFO named by CNAME
* as well as the numbers in IINFO corresponding to DINFO:
*    [CNT] N[CNAME]
*    TIME/[CNAME] [IINFO(1)] [IINFO(2)] [IINFO(3)] ...
* IINFO can e.g. descrbe node numbers, while DINFO describes a function
* value in these nodes.
*
* In:
*   H      - Handle to the output file
*   CNT    - Number of data entries in IINFO/DINFO
*   IINFO  - array [1..CNT] of integer
*            Identifier for the entries in DINFO
*   DINFO  - array [1..CNT] of double
*            Data that should be written out
*   BHEAR  - whether or not to write a headline with CNAME
*   TIMENS - Time identifier for that data set
*   INUM   - Number of the data set.
*            =-1: No number
*   CNAME  - Name of the data set, e.g. "PT" or "BDC",...
************************************************************************

      SUBROUTINE TV2FIL (H,CNT,IINFO,DINFO,BHEAD,TIMENS,INUM,CNAME)
      
      IMPLICIT NONE
      
      INTEGER H,CNT,IINFO(*),INUM
      DOUBLE PRECISION DINFO(*),TIMENS
      CHARACTER CNAME*(*)
      LOGICAL BHEAD
      
      INTEGER I
      
      IF (BHEAD) THEN
        WRITE (H,'(I10,A2,A)') CNT,' N',CNAME
        WRITE (H,'(A$)') 'TIME/BDC '
        DO I=1,CNT-1
          WRITE (H,'(I10$)') IINFO(I)
        END DO
        IF (CNT.GE.1) WRITE (H,'(I10)') IINFO(CNT)
      END IF
      
      IF (INUM.NE.-1) WRITE (H,'(I10$)') INUM
      
      WRITE (H,'(E17.7$)') TIMENS
      DO I=1,CNT-1
        WRITE (H,'(E24.16$)') DINFO(I)
      END DO
      IF (CNT.GE.1) WRITE (H,'(E24.16)') DINFO(CNT)

!      CALL FLUSH(H)

      END
