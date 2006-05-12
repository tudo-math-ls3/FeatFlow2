************************************************************************
* Read start vector from disc
*
* In:
*   VECDAT : array [1..SZN2VI,1..NLMAX] of integer
*            Vector structure of all levels
*   TRIAS  : array [1..SZTRIA,1..NLMAX] of integer
*            Triangulation structures on all relevant levels
*   ISTART : Level of start vector and format.
*            |ISTART|=0: Don't read start vector
*                    =1: Read start vector on level NLMAX,
*                    =2: Read start vector on level NLMAX-1,
*                    =3,4,5,...
*            ISTART > 0: formatted
*            ISTART < 0: unformatted
*   NLMIN  : Minimum level in VECDAT
*   NLMAX  : Maximum level of computation. Start vector will be
*            prolongated to that level.
*   IINT,
*   IAPR,
*   IAVPR,
*   DPREP  : Parameters for the prolongation, in case the start vector
*            must be prolongated from a lower to a higher level
*   CSTART : Name of input file
*
* Out:
*   DUP    : array [1..NEQ] of double
*            The readed solution vector
*
* The "solution" vectors in the VECDAT structure are used for
* temporary data!
************************************************************************
      
      SUBROUTINE RDSTRT (NLMIN,NLMAX,VECDAT,TRIAS,DUP,
     *                   IINT,IAPR,IAVPR,DPREP,
     *                   ISTART,CSTART)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER NLMIN,NLMAX
      INTEGER VECDAT(SZN2VI,NNLEV),TRIAS(SZTRIA,*)
      INTEGER ISTART
      DOUBLE PRECISION DUP
      CHARACTER CSTART*(*)
      INTEGER IINT,IAPR,IAVPR
      DOUBLE PRECISION DPREP
      
      INTEGER KSTART,KSTRT1,KSTRT2,I
      
C     Don't do anything if ISTART is 0
      
      IF (ISTART.EQ.0) RETURN
      
C     ISTART=1: read data on maximum level - directly into DUP.
C     ISTART>2: read data on lower level; cancel if the calculated
C               lower level is less than NLMIN according to VECDAT
      
      IF (ISTART.EQ.1) THEN
        CALL PPRDVC (0,VECDAT(ONEQV,NLMAX),DUP,-1,ISTART,I,
     *              CSTART)
      ELSE
        IF ((NLMAX-ABS(ISTART)+1).LT.NLMIN) THEN
          WRITE (MTERM,'(A)') 'Can''t read start vector from disc!'
          WRITE (MTERM,'(A)') 'Vector on disc is on invalid level!'
          RETURN
        END IF
        
C       Read the vector
        
        KSTART = L(VECDAT(OLSOL,NLMAX-ABS(ISTART)+1))
        CALL PPRDVC (0,VECDAT(ONEQV,NLMAX-ABS(ISTART)+1),
     *               DWORK(KSTART),-1,ISTART,I,CSTART)
                    
      END IF

      IF (I.NE.0) THEN

        WRITE (MTERM,'(A)') 'Error reading start vector!'

      ELSE

C       Check if the reading routine destroyed everything...

        IF (I.GT.NLMAX) THEN
          WRITE (MTERM,'(A)') 'Start vector too long!'
          WRITE (MTERM,'(A)') 'Internal structures destroyed!'
          STOP
        END IF
        
C       Prolongate to current level, if |ISTART|<>1
        
        IF (ABS(ISTART).GT.1) THEN
          DO I=NLMAX-ABS(ISTART)+1,NLMAX-1
            KSTRT1 = L(VECDAT(OLSOL,I))
            KSTRT2 = L(VECDAT(OLSOL,I+1))
            CALL PROLUP (VECDAT(ONU,I),VECDAT(ONU,I+1),
     *                   DWORK(KSTRT1),DWORK(KSTRT2),
     *                   TRIAS(1,I),TRIAS(1,I+1),
     *                   IINT,IAPR,IAVPR,DPREP)
          END DO

C         Copy the resulting vector to the start vector

          CALL LCP1(DWORK(KSTRT2),DUP,VECDAT(ONEQV,NLMAX))
        
        END IF
        
      END IF
      

      END 
      
