************************************************************************
* Build B-matrices
*
* Build the B-matrix blocks in the system matrix by exact evaluation.
* This is possible due to the fact that the pressure basis functions
* are piecewise constant and the velocity basis functions
* piecewise linear. So there is an exact quadrature rule!
*
* Works only if Q1~ is used!
*
* In:
*   NEL,
*   NVT,
*   NMT,
*   KVERT,
*   KMID,
*   KADJ,
*   DCORVG  - usual geometry information
*   B1      - array [1..2*NMT] of double
*             Memory for first B-matrix block
*   B2      - array [1..2*NMT] of double
*             Memory for 2nd B-matrix block
*   KCOLB   - array [1..2*NMT] of double
*             Memory for column structure
*   KLDB    - array [1..NMT+1] of double
*             Memory for row structure structure
*   
* Out:
*   NB      - Number of entries in the B-matrices
*   B1      - The first NB entries are filled with data
*   B2      - The first NB entries are filled with data
*   KCOL    - The first NB entries are filled with data
*   B1      - Is filled with data
*
* After calling this routine, if the caller allocated memory for the
* B-matrices dynamically, the memory for B1,B2,KCOL can be released
* except for the actual NB entries.
************************************************************************

      SUBROUTINE BBUILD(KVERT,KMID,KADJ,DCORVG,B1,B2,KCOLB,KLDB,
     *                  NB,NEL,NVT,NMT)

      IMPLICIT NONE

C Common blocks
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cout.inc'
      
C parameters

      INTEGER KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER KCOLB(*),KLDB(*)
      
      INTEGER NB,NEL,NVT,NMT

C local variables

      INTEGER ILD, IADJ, IVEN, IMID, IVT, IVTN, IEL, IVE

      DOUBLE PRECISION B1(*),B2(*)
      DOUBLE PRECISION PX, PY, PXN, PYN, DN1, DN2

      ILD=0

C Loop through all elements

      DO IEL=1,NEL

C On each element loop through the vertices

        DO IVE=1,NNVE
      
C What's the neighbour element next to the edge
C starting with our current vertex?
      
          IADJ=KADJ(IVE,IEL)
          IF ((IADJ.LE.0).OR.(IADJ.GE.IEL)) THEN

C What's the next vertex = endpoint of our current edge?

            IVEN=IVE+1
            IF (IVEN.EQ.5) IVEN=1

C Calculate the number of our current vertex, the neighbour
C vertex and the midpoint inbetween...

            IMID=KMID (IVE ,IEL)-NVT
            IVT =KVERT(IVE ,IEL)
            IVTN=KVERT(IVEN,IEL)

C Then the coordinates of the two vertices

            PX =DCORVG(1,IVT)
            PY =DCORVG(2,IVT)
            PXN=DCORVG(1,IVTN)
            PYN=DCORVG(2,IVTN)

C and the length of the edge in X- and Y-direction, resp.

            DN1=-PYN+PY
            DN2= PXN-PX

C calculate the value in the B-matrix in the current row
C by exact quadrature 

            ILD=ILD+1
            KLDB (IMID)=ILD
            KCOLB(ILD )=IEL
            B1   (ILD )=DN1
            B2   (ILD )=DN2

            IF (IADJ.GT.0) THEN
              ILD=ILD+1
              KCOLB(ILD)= IADJ
              B1   (ILD)=-DN1
              B2   (ILD)=-DN2
            END IF

          END IF ! ((IADJ.LE.0).OR.(IADJ.GE.IEL))

        END DO ! IVE

      END DO ! IEL

      NB         =ILD
      KLDB(NMT+1)=ILD+1

      END

************************************************************************
* Matrix multiplication routine with B-matrix, matrix-free
*
* Build the matrix-vector-product with the B-matrices elementwise,
* i.e. without existing B-matrices. This realizes IPRECB=2.
*
* Can only be used with Q1~/Q0 !
*
* In:
*   NEL,NVT,NMT,
*   KVERT,
*   KMID,
*   KADJ,
*   DCORVG - usual geometry information
*   DP     - array [1..NEL] of double
*            Pressure vector
*   DF1    - array [1..NMT] of double
*            Velocity/RHS vector, 1st component
*   DF2    - array [1..NMT] of double
*            Velocity/RHS vector, 2nd component
*   A1     - multiplication factor for DF1
*   A2     - multiplication factor for DF2
*
* Out:
*   DFx = A1*Bx*DP + A2*DFx
************************************************************************
      
      SUBROUTINE BMUL1 (KVERT,KMID,KADJ,DCORVG,DP,DF1,DF2,
     *                  NEL,NVT,NMT,A1,A2)
      
      IMPLICIT NONE

C Common blocks
      
      INCLUDE 'cbasictria.inc'
      
C parameters

      INTEGER KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER NEL,NVT,NMT
      DOUBLE PRECISION DP(NEL),DF1(NMT),DF2(NMT),A1,A2

C local variables

      INTEGER IEL, IVE, IVEN, IMID, IVT, IVTN, IADJ
      DOUBLE PRECISION PX, PY, PXN, PYN, DN1, DN2

C     Loop over the elements

      DO IEL=1,NEL

C       Loop over the vertices

        DO IVE=1,NNVE
        
          IADJ=KADJ(IVE,IEL)
          
          IF ((IADJ.LE.0).OR.(IADJ.GE.IEL)) THEN

C           Build the entry of the B-matrix by exact quadrature
C           (which luckily exists for Q1~/Q0), like in BBUILD above:

            IVEN=IVE+1
            IF (IVEN.EQ.5) IVEN=1

            IMID=KMID (IVE ,IEL)-NVT
            IVT =KVERT(IVE ,IEL)
            IVTN=KVERT(IVEN,IEL)

            PX =DCORVG(1,IVT)
            PY =DCORVG(2,IVT)
            PXN=DCORVG(1,IVTN)
            PYN=DCORVG(2,IVTN)

            DN1=-PYN+PY
            DN2= PXN-PX

C           Calculate the matrix-vector product

            IF (IADJ.GT.0) THEN
              DF1(IMID)=A2*DF1(IMID)+A1*DN1*(DP(IEL)-DP(IADJ))
              DF2(IMID)=A2*DF2(IMID)+A1*DN2*(DP(IEL)-DP(IADJ))
            ELSE
              DF1(IMID)=A2*DF1(IMID)+A1*DN1*DP(IEL)
              DF2(IMID)=A2*DF2(IMID)+A1*DN2*DP(IEL)
            ENDIF
          
          END IF ! ((IADJ.LE.0).OR.(IADJ.GE.IEL))

        END DO ! IVE

      END DO ! IEL

      END

************************************************************************
* Matrix multiplication routine. Multiply with B^T, matrix-free.
*
* Build the matrix-vector-product with the B^T-matrices elementwise,
* i.e. without existing B-matrices. This realizes IPRECB=2.
*
* Can only be used with Q1~/Q0 !
*
* In:
*   NEL,NVT,NMT,
*   KVERT,
*   KMID,
*   KADJ,
*   DCORVG - usual geometry information
*   DU1    - array [1..NMT] of double
*            Velocity/RHS vector, 1st component
*   DU2    - array [1..NMT] of double
*            Velocity/RHS vector, 2nd component
*   A1     - multiplication factor for DF1
*
* Out:
*   DFP    - array [1..NEL] of double
*            DFP  =  A1*B^T*DU  =  A1 * (B1^T*DU1 + B2^T*DU2)
************************************************************************
      
      SUBROUTINE BTMUL1(KVERT,KMID,KADJ,DCORVG,DU1,DU2,DFP,
     *                  NEL,NVT,NMT,A1)
      
      IMPLICIT NONE

C Common blocks
      
      INCLUDE 'cbasictria.inc'
      
C parameters

      INTEGER KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER NEL,NVT,NMT
      DOUBLE PRECISION DU1(NMT),DU2(NMT),DFP(NEL),A1

C local variables

      INTEGER IEL, IVE, IVEN, IMID, IVT, IVTN
      DOUBLE PRECISION PX, PY, PXN, PYN, DN1, DN2

C     Loop over the elements

      DO IEL=1,NEL
      
        DFP(IEL)=0D0

C       Loop over the vertices

        DO IVE=1,NNVE
        
C         Build the matrix entry as in BBUILD:
        
          IVEN=IVE+1
          IF (IVEN.EQ.5) IVEN=1

          IMID=KMID (IVE ,IEL)-NVT
          IVT =KVERT(IVE ,IEL)
          IVTN=KVERT(IVEN,IEL)

          PX =DCORVG(1,IVT)
          PY =DCORVG(2,IVT)
          PXN=DCORVG(1,IVTN)
          PYN=DCORVG(2,IVTN)

          DN1=-PYN+PY
          DN2= PXN-PX

          DFP(IEL)=DFP(IEL)+A1*(DN1*DU1(IMID)+DN2*DU2(IMID))

        END DO ! IVE

      END DO ! IEL

      END
