************************************************************************
* Extended VANCA preconditioner, SCGS method
*
* Extended calling convention
*
* This routine performs one preconditioning step
*
*    (du,dp) = Vanca(F)
*
* with the full-block-Vanca preconditioner. (du,dp) must be =0 
* initially. VANCP4 executes a local Block-Gauss-Seidel preconditioner
* onto the vector F to calculate a preconditioned vector (du,dp).
* In contrast to VANCP4, this routine uses a stronger approach which
* involves the solution of a local 4x4 system for every DOF.
*
* In:
*   U1,
*   U2     - array [1..NU] of double precision
*   P      - array [1..NP] of double precision
*            Velocity and Pressure components of the preconditioned
*            vector; must be =0 initially!
*   F1,
*   F2     - array [1..NU] of double precision
*   FP     - array [1..NP] of double precision
*            Velocity and pressure components of the vector that should
*            be preconditioned.
*   A,
*   KCOLA,
*   KLDA   - Current system matrix for velocity components
*
*   B1,
*   B2,
*   KCOLB,
*   KLDB   - The two pressure matrices
*
*   NU     - Length of velocity vectors and dimension of A
*   NP     - Length of pressure vectors and dimension of B1, B2
*
*   KMBD,
*   KVERT,
*   KMID,
*   KXNPR,
*   NVT,
*   NEL    - Usual geometry information;
*            KXNPR(.) decides on Dirichlet/Neumann nodes.
*
*   RLXSM  - Relaxation parameter; 0 <= RLXSM < 1
*
* Out:
*   U1,U2,P - Preconditioned vector, velocity and pressure
***********************************************************************

      SUBROUTINE  VANCP5 (U1,U2,P,F1,F2,FP,
     *                    A,KCOLA,KLDA,B1,B2,KCOLB,KLDB,NU,NP,
     *                    KMID,KXNPR,NEL,NVT)

      IMPLICIT NONE
      
C common blocks      
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      
C parameters

      DOUBLE PRECISION U1(*),U2(*),P(*)
      DOUBLE PRECISION F1(*),F2(*),FP(*)
      DOUBLE PRECISION A(*),B1(*),B2(*)
      INTEGER KCOLA(*),KLDA(*),KCOLB(*),KLDB(*)
      INTEGER KMID(NNVE,*),KXNPR(2,*)
      INTEGER NU,NP
      INTEGER NEL,NVT
      DOUBLE PRECISION A1,A2
      
C     Local arrays for informations about one element
      
      DOUBLE PRECISION AA(4,4),BB1(4),BB2(4),FF1(4),FF2(4)
      INTEGER IU(4)
      LOGICAL BDBC(4)
      
C     local variables

      INTEGER IEL,ILD,IGD,IMID,IA,IB,IB1,IB2,JP1,JP2,ICOL,COL
      INTEGER IA1,IA2
      DOUBLE PRECISION FFP,PJP

C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================

C     Basic algorithm:
C
C     What are we doing here? Well, we want to perform 
C     *preconditioning*, i.e. we have to solve the problem
C
C       x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
C
C     for a "special" preconditioner C which we define in a moment.
C     This is equivalent to solving the system
C
C       C (x_new)  = x_old
C
C     C should be some approximation to A. Imagine our global system:
C
C         [ A   B ] (u) = (f)
C         [ B^t 0 ] (p)   (g)
C
C     In the Navier-Stokes equations with (u,p) being the preconditioned
C     vector, there should be g=0 - but this cannot be assumed
C     as it does not happen in general.
C     Now the algorithm for generating a new (u,p) vector from the old
C     one reads roughly as follows:
C
C     a) Restrict to a small part of the domain, in our case to one cell.
C     b) Fetch all the data (velocity, pressure) on that cell. On the
C        first cell, we have only "old" velocity entries. These values
C        are updated and then the calculation proceeds with the 2nd cell.
C
C               old                      new
C            |---X---|                |---X---|
C            |       |                |       |
C        old X   1   X old   -->  new X   1   X new
C            |       |                |       |
C            |---X---|                |---X---|
C               old                      new
C
C        From the second cell on, there might be "old" data and "new" 
C        data on that cell - the old data that has not been updated and
C        perhaps some already updated velocity data from a neighbor cell.
C        
C               new     old                   new     new
C            |---X---|---X---|             |---X---|---X---|
C            |       |       |             |       |       |
C        new X   1   X   2   X old --> new X   1   X   2   X new
C            |       |new    |             |       |newer  |
C            |---X---|---X---|             |---X---|---X---|
C               new     old                   new     new
C
C        These values are updated and then the calculation proceeds
C        with the next cell.
C        As can be seen in the above picture, the "new" node in the
C        middle is even going to be a "newer" node when handled again
C        for the 2nd cell. This is meant by "Gauss-Seldel" character:
C        Information is updated subsequently by using "old" data and
C        "new" data from a previous calculation.
C
C     So we start with a loop over all elements

      DO IEL=1,NEL
      
C       We now have the element
C                                            U3/V3
C       |---------|                       |----X----|
C       |         |                       |         |
C       |   IEL   |   with DOF's    U4/V4 X    P    X U2/V2
C       |         |                       |         |
C       |---------|                       |----X----|
C                                            U1/V1
C      
C       Fetch the pressure P on the current element into FFP
      
        FFP=FP(IEL)

C       Loop over all 4 U-nodes of that element

        DO ILD=1,4

C         Calculate the DOF that belongs to our edge ILD:

          IMID = KMID(ILD,IEL)
          IGD  = IMID-NVT
          
C         Write the number of the edge/node to IU.
          
          IU(ILD) = IGD
          
        END DO
        
C       Sort IU in ascending order so that 
C       it is sorted just like the global solution/rhs vector.
C
C       ONLY FOR DEBUGGING PURPOSES, not necessary in real
C       implementation; The vector IU defines the sorting.
C
C        DO ILD=1,4
C          DO COL=ILD+1,4
C            IF (IU(COL).LT.IU(ILD)) THEN
C              IGD = IU(ILD)
C              IU(ILD) = IU(COL)
C              IU(COL) = IGD
C            END IF
C          END DO
C        END DO
        
C       Now collect information about the nodes.
        
        DO ILD=1,4
          
          IGD = IU(ILD)
          IMID = IGD+NVT
          
C         Write into BDBC whether the edge is a Dirichlet-
C         Boundary edge:
          
          BDBC(ILD) = IAND(KXNPR(1,IMID),2**12).NE.0
          
C         Set FF1/FF2 initially to the value F1/F2 of the right hand
C         side vector that belongs to our current DOF corresponding
C         to II

          FF1(ILD) = F1(IGD)
          FF2(ILD) = F2(IGD)
          
        END DO
        
C       In a second loop build the local matrix and the RHS vectors.
C
C       Here a short description about the approach.
C       Again consider the problem:
C       
C          [ A   B ] (u) = (f)
C          [ B^t 0 ] (p)   (g)
C     
C       We assume, that all components in the vector (u,p) are
C       given - except for the four velocity unknowns and the
C       pressure unknown on the current element; these five unknowns
C       are located anywhere in the (u,p) vector. The idea is to
C       shift "everything known" to the right hand side to obtain
C       a system for only these unknowns!
C
C       Extracting all the lines of the system that correspond to
C       DOF's on our single element IEL results in a rectangular
C       system of the form
C
C          [ === A^ === B~ ] (|) = (f)   
C          [ B~^t       0  ] (u)   (g)   
C                            (|)
C                            (p)
C
C       <=>  [ MAT1 MAT2 B~ ] (u) = (fu)
C            [    ****   0  ] (p)   (fp)
C
C       with A^ being an 4 x (2*NVT) matrix for the two velocity
C       components and B being an (2*4) x 1 matrix that couples the
C       velocities to the pressure on our current element.
C       B~ is a 8 x 2 matrix: As every velocity couples with at most
C       two pressure elements on the neighbour cell, so we have
C       2 columns in the B-matrix.
C
C              IEL                              IEL
C           |--------|             |--------|--------|
C           |        |             |        |        |
C           |   P    |      or     |   Q    X   P    |
C           |        |             |        |        |
C         --|---X----|--           |--------|--------|
C
C       Now, throw all summands to the RHS vector that do
C       not correspond to the velocity/pressure DOF's on our single 
C       element IEL!
C
C       That way, A^ is reduced to a square matrix with two square
C       submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to 
C       two 4 x 1 submatrices (originally, every velocity couples with
C       two pressure elements on the neighbour cell, so we have
C       2 columns in the B-matrix).
C
C       The system roughly looks like:
C
C       [ ....        B1 ] (u1)   (f1)           (| )           (| )           (| )
C       [ .A~.        |  ] (u2)   (| )           (| )           (| )           (| )
C       [ ....        |  ] (u3)   (| )           (| )           (| )           (| )
C       [ ....        |  ] (u4)   (| )           (| )           (| )           (| )
C       [       ....  B2 ] (v1) = (f2) - sum u_i (ci) - sum v_i (di) - sum   Q (qi)
C       [       .A~.  |  ] (v2)   (| )   i>4     (| )   i>4     (| )  Q in N   (| )
C       [       ....  |  ] (v3)   (| )           (| )           (| )           (| )
C       [       ....  |  ] (v4)   (| )           (| )           (| )           (| )
C       [ B1^T  B2^T  0  ] (p )   (fp)           (0 )           (0 )           (0 )
C
C       assuming in the example that element IEL contains u1,u2,u3,u4,
C       v1,v2,v3,v4 and p as DOF's, and using the notation that ci 
C       and di are the i th column of matrix MAT1 and MAT2 above,
C       respectively. Q is the pressure on the neighbour elements in the 
C       neighbourhood N of element IEL and qi the corresponding column 
C       from the B~-matrix.
C       The above system obviously has the form
C       
C           [ S         B1 ] (ul)   (du)
C           [      S    B2 ] (vl) = (dv)
C           [ B1^t B2^t  0 ] (p )   (dp)
C
C       The vector (du,dv,dp) is called "local residuum"; we'll build it
C       in FF1/FF2. The 4x4-matrix S (which is the same for both, X- and
C       Y-velocity here) we build in AA.
C
C       Note: This only holds for inner edges! For boundary edges, we
C       have to do a slight modification.
C
C       Clear the local matrix AA:

        CALL LCL1(AA,16)

C       Loop over the four DOF's of our element.
        
        DO ILD=1,4

C         Each local DOF ILD corresponds to the global DOF IGD=IU(ILD),
C         so we have to extract the line IGD into the row ILD of
C         the matrix AA.

          IGD = IU(ILD)

C         Loop through the line of the matrix. CCNT gives the next position
C         where to write in the local matrix. This is used as an index
C         in COLIDX to find the correct position where to write in the
C         local matrix.

          IA1 = KLDA(IGD)
          IA2 = KLDA(IGD+1)-1
          DO ICOL = IA1,IA2
          
C           Does this column belong to one of our four DOF's, or
C           do we have to throw it to the right hand side?

            IA = KCOLA(ICOL)
            DO COL=1,4
            
C             Use the following line to switch back to simple VANCA:
C              IF ((IA.EQ.IU(COL)).AND.(IA.EQ.IGD)) GOTO 50
C
C             The following line implements the extended VANCA:
              
              IF (IA.EQ.IU(COL)) GOTO 50
              
            END DO
            
C           Not found - this DOF has to be thrown to the RHS.
C           Mark that by setting ICOL > 4.
C           Otherwise, COL is the horizontal index where to write
C           into the local matrix.
            
            COL = 5
            
50          CONTINUE
            
            IF (COL.LE.4) THEN
     
C             This is a matrix entry belonging to one of our
C             DOF's. So note it into the local matrix AA.

              AA(ILD,COL) = A(ICOL)
              
C             Of course, don't modify the RHS or anything else.
     
            ELSE
            
C             This matrix entry must be moved to the right hand
C             side to create the "local residuum". But only for
C             inner edges - boundary edges are handled slightly
C             different!
C
C             Do we have a boundary edge?

              IF (.NOT.BDBC(ILD)) THEN
     
C               No, inner edge.
C               Subtract the A*u part from the RHS for the local
C               residuum. What we throw here to the RHS is therefore 
C               on one hand a set of entries of the old velocities
C               as well as a set of entries which might have been 
C               calculated by a previous update of another cell in 
C               our neighbourhood.
     
                FF1(ILD) = FF1(ILD) - A(ICOL)*U1(IA)
                FF2(ILD) = FF2(ILD) - A(ICOL)*U2(IA)
              
              END IF ! not BDBC(ILD)
              
            END IF ! COL>4
                
          END DO ! ICOL

C         Ok, the velocity handling is finished. Not we come to handle
C         the pressure. 
C         Do we have a boundary edge?

          IF (.NOT.BDBC(ILD)) THEN

C           Get BB1,BB2 and modify FF1,FF2.
C
C           In every line of each B-matrix, there is either one entry
C           (in case of boundary nodes) or there are at most two entries
C           (for inner nodes, which where each velocity component couples
C           with two elements:
C
C                IEL                              IEL
C             |--------|             |--------|--------|
C             |        |             |        |        |
C             |   P1   |      or     |   P2   X   P1   |
C             |        |             |        |        |
C           --|---X----|--           |--------|--------|
C
C           Our current edge is either an inner edge with two elements
C           attached, or a boundary edge with Neumann boundary
C           conditions.
C
C           In the case that two elements meet each other, we subtract
C           the B*(pressure of the neighbour cell) from the RHS.
C
C           In the case that we have only one element, there's
C           only one entry in the row of the B-matrix. We fetch the
C           first and last entry of the line, which will coincide
C           then (IB1=IB2). As there is no neighbour cell, we must
C           not subtract anything from the RHS! Nevertheless we have to
C           note the column numbers in the B-matrices to our BB1/BB2
C           vectors.
C
C           Get the first and last entry from row I of the B-matrices:

            IB1=KLDB(IGD)
            IB2=KLDB(IGD+1)-1
            JP1=KCOLB(IB1)
            JP2=KCOLB(IB2)
            
C           Check which of these two belong to the current element.
C           The corresponding other one is to be subtracted from the RHS
C           vector.

            IF (JP1.EQ.IEL)  THEN
            
              PJP = P(JP2)
              BB1(ILD) = B1(IB1)
              BB2(ILD) = B2(IB1)
              
C             Don't subtract anything from the RHS if there's
C             no neighbour element!

              IF (IB1.NE.IB2) THEN
                FF1(ILD) = FF1(ILD)-B1(IB2)*PJP
                FF2(ILD) = FF2(ILD)-B2(IB2)*PJP
              END IF
              
            ELSE IF (JP2.EQ.IEL)  THEN
            
              PJP = P(JP1)
              BB1(ILD) = B1(IB2)
              BB2(ILD) = B2(IB2)
              
C             Don't subtract anything from the RHS if there's
C             no neighbour element!

              IF (IB1.NE.IB2) THEN
                FF1(ILD) = FF1(ILD)-B1(IB1)*PJP
                FF2(ILD) = FF2(ILD)-B2(IB1)*PJP
              END IF
              
            ELSE
              WRITE(MTERM,'(A)') 
     *          'ERROR in VANCP5: IEL entry in B not found'
              STOP
            ENDIF
          
          ELSE
          
C           Hmmm, the edge ILD is a Dirichlet boundary edge.
C           That means there is no adjacent element with a pressure 
C           that we have to subtract from the RHS vector.
C           Simply write the one and only element in the current 
C           line of the B-matrices to BB1/BB2.

            IB = KLDB(IGD)
            BB1(ILD) = B1(IB)
            BB2(ILD) = B2(IB)
          
          END IF ! not BDBC(ILD)
          
C         Go to the next DOF on our current element IEL...

        END DO ! ILD

C       Call the subroutine ELUPD3 for element update.
C
C       This routine will now "solve" the above system
C
C          [ S         B1 ] (ul)   (du)
C          [      S    B2 ] (vl) = (dv)
C          [ B1^t B2^t  0 ] (p )   (dp)
C
C       and write the new velocity/pressure values (ul,vl,p) to the
C       solution vector (U1,U2,P).
C       Note that although we called (du,dv,dp) "local defect", it's
C       not a real defect! (du,dv,dp) is the original RHS vector
C       where a lot of parts of the original system have been
C       subtracted in equivalent formulations! Therefore, the solution
C       of this local system really gives new values for velocity/
C       pressure and not "update" values like in a defect correction
C       approach!

        CALL  ELUPD3 (U1,U2,P,IEL,IU,BDBC,AA,BB1,BB2,FF1,FF2,FFP)

C       Go to the next element...

      END DO ! IEL
      
C=======================================================================
C
C     At this point we are finished calculating a the preconditioned
C     vector (du,dv,dp). The caller can use it e.g. in a defect
C     correction approach...

      END
