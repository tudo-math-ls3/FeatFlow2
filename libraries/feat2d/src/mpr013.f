************************************************************************
* Prolongation/Restriction for the Q2 element
************************************************************************

************************************************************************
* Q2 prolongation, simple linear interpolation 
*
* Transfers a Q2-solution from a coarser level to a finer level.
* Uses only linear interpolation on a once refined grid.
*
* In:
*   DU1    : array [1..NVTC+NMTC] of double
*            Coarse grid solution vector.
*   KVERT1 : array [1..NNVE,1..NVTC] of double
*            Vertices on each element on the coarse grid
*   KVERT2 : array [1..NNVE,1..NVTF] of double
*            Vertices on each element on the fine grid
*   KMID1  : array [1..NNVE,1..NVTC] of double
*            Midpoints on each element on the coarse grid
*   KMID2  : array [1..NNVE,1..NVTF] of double
*            Midpoints on each element on the fine grid
*   KADJ1  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            coarse grid.
*   KADJ2  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            fine grid.
*   NEL1   : Number of elements on the coarse grid
*   NEL2   : Number of elements on the fine grid
*   NVT1   : Number of vertices on the coarse grid
*   NVT2   : Number of vertices on the fine grid
*   NMT1   : Number of edge midpoints on the coarse grid
*   NMT2   : Number of edge midpoints on the fine grid
*
* Out:
*   DU2    : array [1..NVTF+NMTF] of double
*            Fine grid solution vector
************************************************************************
 
      SUBROUTINE MP013L(DU1,DU2,NVT1,NMT1,NEL1,NVT2,NMT2,NEL2,
     *                  KVERT1,KMID1,KADJ1,KVERT2,KMID2,KADJ2)
     
      IMPLICIT NONE
     
      INTEGER NNVE
      PARAMETER (NNVE=4)
     
C     parameters
     
      DOUBLE PRECISION DU1(*),DU2(*)
      INTEGER NEL1,NEL2,NVT1,NVT2,NMT1,NMT2
      INTEGER KVERT1(NNVE,*),KVERT2(NNVE,*),
     *        KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),
     *        KADJ2(NNVE,*)

C     constants

      DOUBLE PRECISION Q2,Q4
      PARAMETER (Q2=.5D0,Q4=.25D0)
      
C     local variables      
      
      INTEGER IEL2(4),IEL1,I,J
      
C     Initialize the fine grid vector with 0.
      
      CALL  LCL1 (DU2,NVT2+NMT2+NEL2)

C     Loop over the elements.

      DO IEL1=1,NEL1
     
        IEL2(1)=IEL1
        IEL2(2)=KADJ2(2,IEL2(1))
        IEL2(3)=KADJ2(2,IEL2(2))
        IEL2(4)=KADJ2(2,IEL2(3))
        
        DO I=1,4
          DU2(KVERT2(1,IEL2(I)))=DU1(KVERT1(I,IEL1))
          DU2(KVERT2(2,IEL2(I)))=DU1(KMID1(I,IEL1))
          DU2(KVERT2(3,IEL2(I)))=DU1(NVT1+NMT1+IEL1)
          DU2(KVERT2(4,IEL2(I)))=DU1(KMID1(MOD(I+2,4)+1,IEL1))
        ENDDO
        
        DO I=1,4
          DO J=1,4
            DU2(KMID2(J,IEL2(I)))=
     *           +Q2*DU2(KVERT2(J,IEL2(I)))
     *           +Q2*DU2(KVERT2(MOD(J,4)+1,IEL2(I)))
          ENDDO
        ENDDO
        
        DO I=1,4
          DO J=1,4
            DU2(NVT2+NMT2+IEL2(I))=DU2(NVT2+NMT2+IEL2(I))
     *          +Q4*DU2(KVERT2(J,IEL2(I)))
          ENDDO
        ENDDO

      ENDDO
      
      END

************************************************************************
* Q2 restriction, simple linear interpolation 
*
* Transfers a Q2-solution from a finer level to a coarser level.
* Uses only linear interpolation on a once refined grid.
*
* In:
*   DU1    : array [1..NVTC+NMTC] of double
*            Fine grid solution vector.
*   KVERT1 : array [1..NNVE,1..NVTC] of double
*            Vertices on each element on the fine grid
*   KVERT2 : array [1..NNVE,1..NVTF] of double
*            Vertices on each element on the coarse grid
*   KMID1  : array [1..NNVE,1..NVTC] of double
*            Midpoints on each element on the fine grid
*   KMID2  : array [1..NNVE,1..NVTF] of double
*            Midpoints on each element on the coarse grid
*   KADJ1  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            fine grid.
*   KADJ2  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            coarse grid.
*   NEL1   : Number of elements on the fine grid
*   NEL2   : Number of elements on the coarse grid
*   NVT1   : Number of vertices on the fine grid
*   NVT2   : Number of vertices on the coarse grid
*   NMT1   : Number of edge midpoints on the fine grid
*   NMT2   : Number of edge midpoints on the coarse grid
*
* Out:
*   DU2    : array [1..NVTF+NMTF] of double
*            Coarse grid solution vector
************************************************************************
 
      SUBROUTINE MR013L(DU1,DU2,NVT1,NMT1,NEL1,NVT2,NMT2,NEL2,
     *                  KVERT1,KMID1,KADJ1,KVERT2,KMID2,KADJ2)

      IMPLICIT NONE
      
      INTEGER NNVE
      PARAMETER (NNVE=4)
     
C     parameters
     
      DOUBLE PRECISION DU1(*),DU2(*)
      INTEGER NEL1,NEL2,NVT1,NVT2,NMT1,NMT2
      INTEGER KVERT1(NNVE,*),KVERT2(NNVE,*),
     *        KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),
     *        KADJ2(NNVE,*)
     
C     constants
     
      DOUBLE PRECISION Q2,Q4
      PARAMETER (Q2=.5D0,Q4=.25D0)
      
C     local variables
      
      INTEGER IEL1(4),IEL2,I,J

C     Initialize the coarse grid vector with 0.

      CALL LCL1(DU2,NVT2+NMT2+NEL2)

      CALL LCP1(DU1,DU2,NVT2+NMT2+NEL2)

C     Loop over the elements

      DO IEL2=1,NEL2
     
        IEL1(1)=IEL2
        IEL1(2)=KADJ1(2,IEL1(1))
        IEL1(3)=KADJ1(2,IEL1(2))
        IEL1(4)=KADJ1(2,IEL1(3))

        DO I=1,4
          DU2(KVERT2(I,IEL2))=DU2(KVERT2(I,IEL2))
     *         +0.5*Q2*DU1(KMID1(1,IEL1(I)))
     *         +0.5*Q2*DU1(KMID1(4,IEL1(I)))
     *         +Q4*DU1(NVT1+NMT1+IEL1(I))
          IF(KADJ1(1,IEL1(I)).EQ.0) 
     *         DU2(KVERT2(I,IEL2))=DU2(KVERT2(I,IEL2))
     *         +0.5*Q2*DU1(KMID1(1,IEL1(I)))
          IF(KADJ1(4,IEL1(I)).EQ.0) 
     *         DU2(KVERT2(I,IEL2))=DU2(KVERT2(I,IEL2))
     *         +0.5*Q2*DU1(KMID1(4,IEL1(I)))
          DU2(KMID2(I,IEL2))=DU2(KMID2(I,IEL2))
     *         +0.5*Q2*DU1(KMID1(1,KADJ1(3,IEL1(I))))
     *         +Q2*DU1(KMID1(3,IEL1(I)))
     *         +0.5*Q2*DU1(KMID1(4,IEL1(I)))
     *         +Q4*DU1(NVT1+NMT1+IEL1(I))
     *         +Q4*DU1(NVT1+NMT1+KADJ1(3,IEL1(I)))
          IF(KADJ1(1,KADJ1(3,IEL1(I))).EQ.0) 
     *         DU2(KMID2(I,IEL2))=DU2(KMID2(I,IEL2))
     *         +0.5*Q2*DU1(KMID1(1,KADJ1(3,IEL1(I))))
          IF(KADJ1(4,IEL1(I)).EQ.0) 
     *         DU2(KMID2(I,IEL2))=DU2(KMID2(I,IEL2))
     *         +0.5*Q2*DU1(KMID1(4,IEL1(I)))
        ENDDO
        DU2(NVT2+NMT2+IEL2)=DU2(NVT2+NMT2+IEL2)
     *       +Q2*DU1(KMID1(2,IEL1(1)))
     *       +Q2*DU1(KMID1(2,IEL1(2)))
     *       +Q2*DU1(KMID1(2,IEL1(3)))
     *       +Q2*DU1(KMID1(2,IEL1(4)))
     *       +Q4*DU1(NVT1+NMT1+IEL1(1))
     *       +Q4*DU1(NVT1+NMT1+IEL1(2))
     *       +Q4*DU1(NVT1+NMT1+IEL1(3))
     *       +Q4*DU1(NVT1+NMT1+IEL1(4))
        
      ENDDO
      
      END
     
************************************************************************
* Q2 prolongation, full quadratic interpolation 
*
* Transfers a Q2-solution from a coarser level to a finer level.
*
* In:
*   DU1    : array [1..NVTC+NMTC] of double
*            Coarse grid solution vector.
*   KVERT1 : array [1..NNVE,1..NVTC] of double
*            Vertices on each element on the coarse grid
*   KVERT2 : array [1..NNVE,1..NVTF] of double
*            Vertices on each element on the fine grid
*   KMID1  : array [1..NNVE,1..NVTC] of double
*            Midpoints on each element on the coarse grid
*   KMID2  : array [1..NNVE,1..NVTF] of double
*            Midpoints on each element on the fine grid
*   KADJ1  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            coarse grid.
*   KADJ2  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            fine grid.
*   NEL1   : Number of elements on the coarse grid
*   NEL2   : Number of elements on the fine grid
*   NVT1   : Number of vertices on the coarse grid
*   NVT2   : Number of vertices on the fine grid
*   NMT1   : Number of edge midpoints on the coarse grid
*   NMT2   : Number of edge midpoints on the fine grid
*
* Out:
*   DU2    : array [1..NVTF+NMTF] of double
*            Fine grid solution vector
************************************************************************

      SUBROUTINE MP013 (DU1,DU2,NVT1,NMT1,NEL1,NVT2,NMT2,NEL2,
     *                  KVERT1,KMID1,KADJ1,KVERT2,KMID2,KADJ2)
     
      IMPLICIT NONE
     
      INTEGER NNVE
      PARAMETER (NNVE=4)
     
C     parameters
     
      DOUBLE PRECISION DU1(*),DU2(*)
      INTEGER NEL1,NEL2,NVT1,NVT2,NMT1,NMT2
      INTEGER KVERT1(NNVE,*),KVERT2(NNVE,*),
     *        KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),
     *        KADJ2(NNVE,*)

C     constants

      DOUBLE PRECISION Q2,Q4
      PARAMETER (Q2=.5D0,Q4=.25D0)
      
C     local variables      
      
      INTEGER IEL2(4),IEL1,I,J
      
C     Initialize the fine grid vector with 0.
      
      CALL LCL1(DU2,NVT2+NMT2+NEL2)

      CALL LCP1(DU1,DU2,NVT2+NMT2+NEL2)

C     Loop over the elements

      DO IEL1=1,NEL1
     
        IEL2(1)=IEL1
        IEL2(2)=KADJ2(2,IEL2(1))
        IEL2(3)=KADJ2(2,IEL2(2))
        IEL2(4)=KADJ2(2,IEL2(3))
        
        DO I=1,4
          DU2(KMID2(1,IEL2(I)))=
     *         +(3.0/8.0)*DU2(KVERT2(1,IEL2(I)))
     *         +(3.0/4.0)*DU2(KVERT2(2,IEL2(I)))
     *         -(1.0/8.0)*DU2(KVERT2(1,IEL2(MOD(I,4)+1)))
          DU2(KMID2(4,IEL2(I)))=
     *         +(3.0/8.0)*DU2(KVERT2(1,IEL2(I)))
     *         +(3.0/4.0)*DU2(KVERT2(4,IEL2(I)))
     *         -(1.0/8.0)*DU2(KVERT2(1,IEL2(MOD(I+2,4)+1)))
          DU2(KMID2(2,IEL2(I)))=
     *         +(3.0/8.0)*DU2(KVERT2(2,IEL2(I)))
     *         +(3.0/4.0)*DU2(KVERT2(3,IEL2(I)))
     *         -(1.0/8.0)*DU2(KVERT2(4,IEL2(MOD(I+2,4)+1)))
          DU2(NVT2+NMT2+IEL2(I))=
     *         +(9.0/64.0)*DU2(KVERT2(1,IEL2(I)))
     *         +(18.0/64.0)*DU2(KVERT2(2,IEL2(I)))
     *         +(36.0/64.0)*DU2(KVERT2(3,IEL2(I)))
     *         +(18.0/64.0)*DU2(KVERT2(4,IEL2(I)))
     *         -(3.0/64.0)*DU2(KVERT2(1,IEL2(MOD(I,4)+1)))
     *         -(6.0/64.0)*DU2(KVERT2(2,IEL2(MOD(I,4)+1)))
     *         -(3.0/64.0)*DU2(KVERT2(1,IEL2(MOD(I+2,4)+1)))
     *         -(6.0/64.0)*DU2(KVERT2(4,IEL2(MOD(I+2,4)+1)))
     *         +(1.0/64.0)*DU2(KVERT2(1,IEL2(MOD(I+1,4)+1)))
        ENDDO
      ENDDO
      
      END

************************************************************************
* Q2 restriction, full quadratic interpolation 
*
* Transfers a Q2-solution from a finer level to a coarser level.
*
* In:
*   DU1    : array [1..NVTC+NMTC] of double
*            Fine grid solution vector.
*   KVERT1 : array [1..NNVE,1..NVTC] of double
*            Vertices on each element on the fine grid
*   KVERT2 : array [1..NNVE,1..NVTF] of double
*            Vertices on each element on the coarse grid
*   KMID1  : array [1..NNVE,1..NVTC] of double
*            Midpoints on each element on the fine grid
*   KMID2  : array [1..NNVE,1..NVTF] of double
*            Midpoints on each element on the coarse grid
*   KADJ1  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            fine grid.
*   KADJ2  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            coarse grid.
*   NEL1   : Number of elements on the fine grid
*   NEL2   : Number of elements on the coarse grid
*   NVT1   : Number of vertices on the fine grid
*   NVT2   : Number of vertices on the coarse grid
*   NMT1   : Number of edge midpoints on the fine grid
*   NMT2   : Number of edge midpoints on the coarse grid
*
* Out:
*   DU2    : array [1..NVTF+NMTF] of double
*            Coarse grid solution vector
************************************************************************
 
      SUBROUTINE MR013 (DU1,DU2,NVT1,NMT1,NEL1,NVT2,NMT2,NEL2,
     *                  KVERT1,KMID1,KADJ1,KVERT2,KMID2,KADJ2)

      IMPLICIT NONE
      
      INTEGER NNVE
      PARAMETER (NNVE=4)
     
C     parameters
     
      DOUBLE PRECISION DU1(*),DU2(*)
      INTEGER NEL1,NEL2,NVT1,NVT2,NMT1,NMT2
      INTEGER KVERT1(NNVE,*),KVERT2(NNVE,*),
     *        KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),
     *        KADJ2(NNVE,*)
     
C     constants
     
      DOUBLE PRECISION Q2,Q4
      PARAMETER (Q2=.5D0,Q4=.25D0)
      
C     local variables
      
      INTEGER IEL1(4),IEL2,I,J

C     Initialize the coarse grid vector with 0.

      CALL LCL1(DU2,NVT2+NMT2+NEL2)
      CALL LCP1(DU1,DU2,NVT2+NMT2+NEL2)

      DO IEL2=1,NEL2
     
        IEL1(1)=IEL2
        IEL1(2)=KADJ1(2,IEL1(1))
        IEL1(3)=KADJ1(2,IEL1(2))
        IEL1(4)=KADJ1(2,IEL1(3))

        DO I=1,4
          DU2(KVERT2(I,IEL2))=DU2(KVERT2(I,IEL2))
     *         +0.5*(24.0/64.0)*DU1(KMID1(1,IEL1(I)))
     *         +0.5*(24.0/64.0)*DU1(KMID1(4,IEL1(I)))
     *         -0.5*(8.0/64.0)*DU1(KMID1(4,IEL1(MOD(I,4)+1)))
     *         -0.5*(8.0/64.0)*DU1(KMID1(1,IEL1(MOD(I+2,4)+1)))
     *         +(9.0/64.0)*DU1(NVT1+NMT1+IEL1(I))
     *         -(3.0/64.0)*DU1(NVT1+NMT1+IEL1(MOD(I,4)+1))
     *         +(1.0/64.0)*DU1(NVT1+NMT1+IEL1(MOD(I+1,4)+1))
     *         -(3.0/64.0)*DU1(NVT1+NMT1+IEL1(MOD(I+2,4)+1))
          IF(KADJ1(1,IEL1(I)).EQ.0)
     *         DU2(KVERT2(I,IEL2))=DU2(KVERT2(I,IEL2))
     *         +0.5*(24.0/64.0)*DU1(KMID1(1,IEL1(I)))
          IF(KADJ1(4,IEL1(I)).EQ.0)
     *         DU2(KVERT2(I,IEL2))=DU2(KVERT2(I,IEL2))
     *         +0.5*(24.0/64.0)*DU1(KMID1(4,IEL1(I)))
          IF(KADJ1(4,IEL1(MOD(I,4)+1)).EQ.0)
     *         DU2(KVERT2(I,IEL2))=DU2(KVERT2(I,IEL2))
     *         -0.5*(8.0/64.0)*DU1(KMID1(4,IEL1(MOD(I,4)+1)))
          IF(KADJ1(1,IEL1(MOD(I+2,4)+1)).EQ.0)
     *         DU2(KVERT2(I,IEL2))=DU2(KVERT2(I,IEL2))
     *         -0.5*(8.0/64.0)*DU1(KMID1(1,IEL1(MOD(I+2,4)+1)))


          DU2(KMID2(I,IEL2))=DU2(KMID2(I,IEL2))
     *         +0.5*(48.0/64.0)*DU1(KMID1(1,IEL1(I)))
     *         +(24.0/64.0)*DU1(KMID1(2,IEL1(I)))
     *         +0.5*(48.0/64.0)*DU1(KMID1(4,IEL1(MOD(I,4)+1)))
     *         -(8.0/64.0)*DU1(KMID1(2,IEL1(MOD(I+1,4)+1)))
     *         +(18.0/64.0)*DU1(NVT1+NMT1+IEL1(I))
     *         +(18.0/64.0)*DU1(NVT1+NMT1+IEL1(MOD(I,4)+1))
     *         -(6.0/64.0)*DU1(NVT1+NMT1+IEL1(MOD(I+1,4)+1))
     *         -(6.0/64.0)*DU1(NVT1+NMT1+IEL1(MOD(I+2,4)+1))
          IF(KADJ1(1,IEL1(I)).EQ.0)
     *         DU2(KMID2(I,IEL2))=DU2(KMID2(I,IEL2))
     *         +0.5*(48.0/64.0)*DU1(KMID1(1,IEL1(I)))
          IF(KADJ1(4,IEL1(MOD(I,4)+1)).EQ.0)
     *         DU2(KMID2(I,IEL2))=DU2(KMID2(I,IEL2))
     *         +0.5*(48.0/64.0)*DU1(KMID1(4,IEL1(MOD(I,4)+1)))
C     *         DU2(KVERT2(I,IEL2))=DU2(KVERT2(I,IEL2))
C     *         +0.5*(48.0/64.0)*DU1(KMID1(4,IEL1(MOD(I,4)+1)))

        ENDDO
        DU2(NVT2+NMT2+IEL2)=DU2(NVT2+NMT2+IEL2)
     *       +(48.0/64.0)*DU1(KMID1(2,IEL1(1)))
     *       +(48.0/64.0)*DU1(KMID1(2,IEL1(2)))
     *       +(48.0/64.0)*DU1(KMID1(2,IEL1(3)))
     *       +(48.0/64.0)*DU1(KMID1(2,IEL1(4)))
     *       +(36.0/64.0)*DU1(NVT1+NMT1+IEL1(1))
     *       +(36.0/64.0)*DU1(NVT1+NMT1+IEL1(2))
     *       +(36.0/64.0)*DU1(NVT1+NMT1+IEL1(3))
     *       +(36.0/64.0)*DU1(NVT1+NMT1+IEL1(4))
        
      ENDDO
      
      END
