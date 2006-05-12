************************************************************************
* Extended VANCA smoother, SCGS method
*
* Extended calling convention
*
* This routine performs one smoothing step
*
*    x_new = x  +  RLXSM ( Vanca(x) - x )
*
* with Vanca(x) being a better solution than X to the system Ax=b
* calculated by a local Gauss-Seidel approach. In contrast to VANCS3,
* this routine uses a stronger smoothing approach which involves
* the solution of a local 4x4 system for every DOF.
*
* In:
*   U1,
*   U2     - array [1..NU] of double precision
*   P      - array [1..NP] of double precision
*            Velocity and Pressure components of the current
*            solution vector
*   U1OLD,
*   U2OLD  - array [1..NU] of double precision
*   POLD   - array [1..NP] of double precision
*            Auxiliary vectors
*   F1,
*   F2     - array [1..NU] of double precision
*   FP     - array [1..NP] of double precision
*            Right hand side of the system for velocity and pressure
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
*   U1,U2,P - Smoothed solution for velocity and pressure
************************************************************************

      SUBROUTINE  VANCS3 (U1,U2,P,U1OLD,U2OLD,POLD,F1,F2,FP,
     *                    A,KCOLA,KLDA,B1,B2,KCOLB,KLDB,NU,NP,
     *                    KMID,KXNPR,NEL,NVT,RLXSM)

      IMPLICIT NONE
      
C common blocks      
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      
C parameters

      DOUBLE PRECISION U1(*),U2(*),P(*),U1OLD(*),U2OLD(*),POLD(*)
      DOUBLE PRECISION F1(*),F2(*),FP(*)
      DOUBLE PRECISION A(*),B1(*),B2(*)
      INTEGER KCOLA(*),KLDA(*),KCOLB(*),KLDB(*)
      INTEGER KMID(NNVE,*),KXNPR(2,*)
      INTEGER NU,NP
      INTEGER NEL,NVT
      DOUBLE PRECISION RLXSM,A1,A2
      
C     Local arrays for informations about one element
      
      DOUBLE PRECISION AA(4,4),BB1(4),BB2(4),FF1(4),FF2(4)
      INTEGER IU(4)
      LOGICAL BDBC(4)
      
C     local variables

      INTEGER IEL,ILD,IGD,IMID,IA,IB,IB1,IB2,JP1,JP2,ICOL,COL
      INTEGER IA1,IA2
      DOUBLE PRECISION FFP,PJP

C     Copy UOLD:=U. We use that for updating U later:

      CALL LCP1 (U1,U1OLD,NU)
      CALL LCP1 (U2,U2OLD,NU)
      CALL LCP1 (P ,POLD ,NP)

C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================

C     Basic algorithm:
C
C     What are we doing here? Well, we want to perform *smoothing*.
C     This is an iteration to the solution vector in the form of a
C     defect correction:
C
C       x_n+1  =  x_n + Omega C^-1 (b - A x_n)
C
C     With C=A, this simplifies to the simple update formula:
C
C       x_n+1  =  x_n + Omega (x* - x_n)
C
C     where x* is the solution of the problem. In our case, we replace
C     x* by an approximated "better" solution x calculated from x_n.
C     Furthermore by replacing Omega=1-RLXSM, we obtain the smoothing 
C     formula
C
C       x_n+1  =  RLXSM x_n + (1-RLXSM) x
C
C     Now the big question: How to calculate x from x_n???
C     This is done in a "local Block-Gauss-Seidel" manner. This means,
C     we subsequently update some entries in x_n and use these in
C     later computations to update some other entries. Note that this is
C     not a Gauss-Seidel algorithm in closer sense, but more a
C     Gauss-Seidel type method: We calculate new entries from old ones
C     and from previously calculated new ones!
C
C     This is done in the following way: Imagine our global system:
C
C         [ A   B ] (u) = (f)
C         [ B^t 0 ] (p)   (g)
C
C     In the Navier-Stokes equations with (u,p) being the solution
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
     *          'ERROR in SMOOTH: IEL entry in B not found'
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
C     At this point we are finished calculating a (hopefully) better
C     solution (U1,U2,P) than (U1OLD,U2OLD,POLD).
C
C     Relaxation step.
C
C     We don't use directly (U1,U2,P) (although we could do),
C     but we use a relaxation to "blend" the old solution
C     into the new one. This allowes the user to influence the
C     amount of smoothing. RLXSM=1.0 means: Take the full newly
C     calculated vector as solution.
C
C          U := RLXSM*U + (1-RLXSM)*UOLD
C             = U + (1-RLXSM) (UOLD-U)

      IF (RLXSM.NE.1D0) THEN
        A1= 1D0-RLXSM
        A2=     RLXSM
        CALL  LLC1 (U1OLD,U1,NU,A1,A2)
        CALL  LLC1 (U2OLD,U2,NU,A1,A2)
        CALL  LLC1 (POLD ,P ,NP,A1,A2)
      ENDIF
      
C     That's it, smoothing finished.

      END

************************************************************************
* Extended VANCA solver, SCGS method
*
* Extended calling convention
*
* This routine performs one fix point iteration step
*
*    x_new = x  +  RLXSM ( Vanca(x) - x )
*
* with Vanca(x) being a better solution than X to the system Ax=b
* calculated by a local Gauss-Seidel approach.
* This is basically the same as the Vanca-smoother in VANCS3 above,
* but in difference the routine returns information about
* the maximum velocity/divergence defect to the caller.
* The routine uses a stronger solving approach which involves
* the solution of a local 4x4 system for every DOF.
*
* In:
*   U1,
*   U2     - array [1..NU] of double precision
*   P      - array [1..NP] of double precision
*            Velocity and Pressure components of the current
*            solution vector
*   U1OLD,
*   U2OLD  - array [1..NU] of double precision
*   POLD   - array [1..NP] of double precision
*            Auxiliary vectors
*   F1,
*   F2     - array [1..NU] of double precision
*   FP     - array [1..NP] of double precision
*            Right hand side of the system for velocity and pressure
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
*   KNPR,
*   NVT,
*   NEL    - Usual geometry information;
*            KNPR(.)<>0 identifies boundary edges
*
*   RLXSM  - Relaxation parameter; 0 <= RLXSM < 1
*
* Out:
*   U1,U2,P - Improved solution for velocity and pressure
*   DMAXU   - Maximum velocity defect
*   DMAXP   - Maximum defect in the divergence
************************************************************************

      SUBROUTINE  VANCE3 (U1,U2,P,U1OLD,U2OLD,POLD,F1,F2,FP,
     *                    A,KCOLA,KLDA,B1,B2,KCOLB,KLDB,NU,NP,
     *                    KMID,KXNPR,NEL,NVT,RLXPAR,DMAXU,DMAXP)

      IMPLICIT NONE
      
C common blocks      
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      
C parameters

      DOUBLE PRECISION U1(*),U2(*),P(*),U1OLD(*),U2OLD(*),POLD(*)
      DOUBLE PRECISION F1(*),F2(*),FP(*)
      DOUBLE PRECISION A(*),B1(*),B2(*)
      INTEGER KCOLA(*),KLDA(*),KCOLB(*),KLDB(*)
      INTEGER KMID(NNVE,*),KXNPR(2,*)
      INTEGER NU,NP
      INTEGER NEL,NVT
      DOUBLE PRECISION A1,A2,RLXPAR,DMAXU,DMAXP
      
C     Local arrays for informations about one element
      
      DOUBLE PRECISION AA(4,4),BB1(4),BB2(4),FF1(4),FF2(4)
      INTEGER IU(4)
      LOGICAL BDBC(4)
      
C     local variables

      INTEGER IEL,ILD,IGD,IMID,IA,IB,IB1,IB2,JP1,JP2,ICOL,COL
      INTEGER IA1,IA2
      INTEGER IB3(4),JP3(4)
      DOUBLE PRECISION FFP,PJP,VEC(8),FFPH

C     Copy UOLD:=U. We use that for updating U later:

      CALL LCP1 (U1,U1OLD,NU)
      CALL LCP1 (U2,U2OLD,NU)
      CALL LCP1 (P ,POLD ,NP)

C     DMAXU/DMAXP will be build by taking a maximum; initialize
C     with 0.

      DMAXU=0D0
      DMAXP=0D0

C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================

C     Basic algorithm:
C
C     What are we doing here? Well, we want to perform *smoothing*.
C     This is an iteration to the solution vector in the form of a
C     defect correction:
C
C       x_n+1  =  x_n + Omega C^-1 (b - A x_n)
C
C     With C=A, this simplifies to the simple update formula:
C
C       x_n+1  =  x_n + Omega (x* - x_n)
C
C     where x* is the solution of the problem. In our case, we replace
C     x* by an approximated "better" solution x calculated from x_n.
C     Furthermore by replacing Omega=1-RLXSM, we obtain the smoothing 
C     formula
C
C       x_n+1  =  RLXSM x_n + (1-RLXSM) x
C
C     Now the big question: How to calculate x from x_n???
C     This is done in a "local Block-Gauss-Seidel" manner. This means,
C     we subsequently update some entries in x_n and use these in
C     later computations to update some other entries. Note that this is
C     not a Gauss-Seidel algorithm in closer sense, but more a
C     Gauss-Seidel type method: We calculate new entries from old ones
C     and from previously calculated new ones!
C
C     This is done in the following way: Imagine our global system:
C
C         [ A   B ] (u) = (f)
C         [ B^t 0 ] (p)   (g)
C
C     In the Navier-Stokes equations with (u,p) being the solution
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
        FFPH=FFP

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
C           the B*pressure of the neighbour cell from the RHS.
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
            
              PJP      = P(JP2)
              BB1(ILD) = B1(IB1)
              BB2(ILD) = B2(IB1)
              
C             Don't subtract anything from the RHS if there's
C             no neighbour element!

              IF (IB1.NE.IB2) THEN
                FF1(ILD) = FF1(ILD)-B1(IB2)*PJP
                FF2(ILD) = FF2(ILD)-B2(IB2)*PJP
              END IF
              
              IB3(ILD) = IB1
              JP3(ILD) = JP1
              
            ELSE IF (JP2.EQ.IEL)  THEN
            
              PJP      = P(JP1)
              BB1(ILD) = B1(IB2)
              BB2(ILD) = B2(IB2)
              
C             Don't subtract anything from the RHS if there's
C             no neighbour element!

              IF (IB1.NE.IB2) THEN
                FF1(ILD) = FF1(ILD)-B1(IB1)*PJP     
                FF2(ILD) = FF2(ILD)-B2(IB1)*PJP     
              END IF
              
              IB3(ILD) = IB2
              JP3(ILD) = JP2
              
            ELSE
              WRITE(MTERM,'(A)') 
     *          'ERROR in VANCE3: IEL entry in B not found'
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
        
C       Calculate the maximum norm of the residual!
C
C       Additionally to IBx/JPx we now have:
C
C       PJP = pressure on the neighbour cell
C       IB3 = the index of the B-entry of the neighbour cell
C       JP3 = the column number of the B-entry of the neighbour cell
C
C       FF1/FF2 currently shows the value of 
C            f-(Parts of A)u-(neighbouring pressure)*B
C       By subtracting also (A*u) and the pressure contribution 
C       on the current cell, we form the residuum of the
C       velocity. The maximum of the residuum is stored to DMAXU.

        DO COL=1,4
          VEC(COL)   = FF1(COL)
          VEC(COL+4) = FF2(COL)
        END DO
        
        DO COL=1,4
        
          IF (.NOT.BDBC(COL)) THEN
        
            IGD = IU(COL)
            DO ILD=1,4
              VEC(ILD)   = VEC(ILD)   - AA(ILD,COL)*U1(IGD)
              VEC(ILD+4) = VEC(ILD+4) - AA(ILD,COL)*U2(IGD)
            END DO
            
C           Finally subtract also the (current cell pressure)*B.
C           The column number of the B-matrices corresponds to the
C           DOF in the pressure P. JP3 is the column number of
C           the entry in the B-matrix corresponding to the current
C           cell and thus the DOF of the pressure on the current cell. 

            PJP = P(JP3(COL))
            VEC(COL)   = VEC(COL)   - B1(IB3(COL))*PJP
            VEC(COL+4) = VEC(COL+4) - B2(IB3(COL))*PJP
          
          END IF
          
        END DO
        
C       VEC shows the "local residuum" for the X- and Y-velocity.
C       Build the MAX-norm of it and calculate the local residuum
C       for the cell pressure ("divergence defect") simultaneously.
        
        DO COL=1,4
        
          IGD = IU(COL)
          
          DMAXU = MAX(DMAXU,ABS(VEC(COL)))
          DMAXU = MAX(DMAXU,ABS(VEC(COL+4)))
          
C         FFPH on the other hand is the original RHS of the divergence
C         equation B^t*u=0 (was set that way above). Subtract B^t*u
C         from that to form the "divergence defect" g-B^t*u on the
C         current cell.

          FFPH = FFPH - BB1(COL)*U1(IGD)-BB2(COL)*U2(IGD)
          
        END DO ! COL

C       At last handle the "divergence defect", which is defined
C       cellwise.
C       REMARK: This is a difference to the original implementation
C       of VANCA, which does not use an ABS() here!
        
        DMAXP = MAX(DMAXP,ABS(FFPH))
        
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
C     At this point we are finished calculating a (hopefully) better
C     solution (U1,U2,P) than (U1OLD,U2OLD,POLD).
C
C     Relaxation step.
C
C     We don't use directly (U1,U2,P) (although we could do),
C     but we use a relaxation to "blend" the old solution
C     into the new one. This allowes the user to influence the
C     amount of smoothing. RLXSM=1.0 means: Take the full newly
C     calculated vector as solution.
C
C          U := RLXSM*U + (1-RLXSM)*UOLD
C             = U + (1-RLXSM) (UOLD-U)

      IF (RLXPAR.NE.1D0) THEN
        A1= 1D0-RLXPAR
        A2=     RLXPAR
        CALL  LLC1 (U1OLD,U1,NU,A1,A2)
        CALL  LLC1 (U2OLD,U2,NU,A1,A2)
        CALL  LLC1 (POLD ,P ,NP,A1,A2)
      ENDIF
      
      END

************************************************************************
* Element update for vector in Vanca smoother
*
* This routine performs a block Gauss-Seidel algorithm to update
* the solution vector (U1,U2,P) for all unknowns of element IEL.
* It expects a local 4x4 system which is solved using linear-algebra
* packages.
*
* In:
*   U1,
*   U2,
*   P      - array [1..*] of double precision
*            Velocity and pressure components of solution vector
*   IEL    - The element where the algorithm should work on
*   IU     - array [1..4] of integer
*            Numbers of the global velocity DOF's on element IEL
*   BDBC   - array [1..4] of boolean
*            =true , if the DOF is a boundary node
*            =false, otherwise
*   AA     - array [1..4,1..4] of double
*            Local 4x4 system, extracted from the system matrix 
*            corresponding to the DOF's in IU(.)
*   BB1,
*   BB2    - array [1..4] of double
*            Local pressure matrices, corresponding to the DOF's in
*            IU(.); both are 4x1 matrices!
*   FF1,    
*   FF2    - array [1..4] of double
*            Local RHS vector for velocity, corresponding to the 
*            DOF's in IU(.)
*   FFP    - Pressure RHS on element IEL
*
* Out:
*   The components of (U1,U2,P) belonging to element IEL are updated.
************************************************************************

      SUBROUTINE ELUPD3 (U1,U2,P,IEL,IU,BDBC,AA,BB1,BB2,FF1,FF2,FFP)

      IMPLICIT NONE

      INCLUDE 'cout.inc'

      DOUBLE PRECISION U1(*),U2(*),P(*)

C     parameters

      DOUBLE PRECISION AA(4),BB1(4),BB2(4),FF1(4),FF2(4),FFP
      LOGICAL BDBC(4)
      INTEGER IU(4)

C local variables

      INTEGER IPIV(4)
      INTEGER II, IEL, I
      DOUBLE PRECISION DP, PP
      DOUBLE PRECISION FV(4,4)
      DOUBLE PRECISION AI(4,4),DD1(4),DD2(4),UU(4,2)

C     This routine uses a "Block-Gauss-Seidel" strategy to update the
C     solution vector (U1,U2,UP) locally to a new iterate.
C     "Block-Gauss-Seidel" means here, that the new solution vector
C     is calculated by parts of the old one and by parts of the new
C     one, which were calculated by a previous call to this routine.
C
C     What we want to calculate are two things: 1.) a new pressure and
C     2.) a new velocity. Both can be calculated from the 
C     RHS using the Schur-Complement approach.
C
C     Assume we have a system:
C
C      [ S   B ] (u) = (f)
C      [ B^t 0 ] (p)   (g)
C
C     with >>S, B, u, p, f, g given locally on one patch of the mesh.<<
C     Here S is given as a 4x4 matrix, extracted from the global matrix.
C     We can write:
C     
C                u = S^-1 (f-Bp)
C            B^t u = g
C
C     Inserting the first equation into the second one gives:
C
C           B^t S^-1 (f-Bp) = g
C
C      <=>   -B^t S^-1 B p  =  g - B^t S^-1 f
C            ***********       **************
C               =: DP              =: FFP
C
C     Remark that DP is a 1x1-system, i.e. a scalar! Therefore
C     calculating DP^-1 to get p=DP^-1*FFP on the current cell
C     is trivial! So FFP/DP will be the pressure on the element IEL.
C
C     Calculating an update for the velocity 
C
C          u = S^-1 (f-Bp)
C
C     is then also trivial as S (and thus S^-1) is a full 4x4 matrix.

      DO II=1,4
      
C       DD1/DD2 saves the entries of the B1/B2-matrix - except for
C       if the vertex is a boundary vertex!
      
        IF (BDBC(II))  THEN
          DD1(II)=0D0
          DD2(II)=0D0
        ELSE
          DD1(II)=BB1(II)
          DD2(II)=BB2(II)
        ENDIF
        
      END DO
      
C     Factorization loop
C
C     What we at first want to calculate is p with
C
C             - B^t S^-1 B p  =  g - B^t S^-1 f
C
C     To calculate that for local B, S, f and p, consider at first
C     the dimensions in this system:
C
C     a) B is a 4x1 matrix 
C     b) S^-1 is a 4x4 matrix
C     c) B^t S^-1 B is therefore a 1x1 matrix, thus a scalar
C
C     In the factorization loop we use the representation
C
C       DP  = - (B^T S^-1 B)
C       FFP = g - (B^T S^-1 f)
C
C     Copy the 4x4-matrix AA to AI and factorize with LAPACK
C     to calculate a representation of S^-1:

      CALL LCP1(AA,AI,16)

      CALL DGETRF( 4, 4, AI, 4, IPIV, I )
      
      IF (I.NE.0) THEN
        WRITE(MTERM,*)'ERROR in ELUPD3: singular matrix'
        RETURN
      END IF
      
C     We assemble DP  = - (B^T S^-1 B)
C     and         FFP = g - B^T S^-1 f with f=(FF1,FF2).
C
C     Copy B=(B1,B2)->FV for the calculation of f_new(1/2) = S^-1 B.
C     Copy f->FV for the calculation of f_new(3/4) = S^-1 f.
C     Don't use BLAS since direct processing is faster with
C     8 doubles...

      FV(1,1) = DD1(1)
      FV(2,1) = DD1(2)
      FV(3,1) = DD1(3)
      FV(4,1) = DD1(4)

      FV(1,2) = DD2(1)
      FV(2,2) = DD2(2)
      FV(3,2) = DD2(3)
      FV(4,2) = DD2(4)
      
      FV(1,3) = FF1(1)
      FV(2,3) = FF1(2)
      FV(3,3) = FF1(3)
      FV(4,3) = FF1(4)
      
      FV(1,4) = FF2(1)
      FV(2,4) = FF2(2)
      FV(3,4) = FF2(3)
      FV(4,4) = FF2(4)

C     Solve these four systems simultaneously
      
      CALL DGETRS( 'N', 4, 4, AI, 4, IPIV, FV, 4, I )

      IF (I.NE.0) THEN
        WRITE(MTERM,*)'ERROR in ELUPD3: unable to calculate S^-1 f/B'
        RETURN
      END IF

C     Calculate DP = -B^t FV(1/2)
C     Calculate FFP = g - B^t FV(3/4).

      DP=0D0
      
      DO II = 1,4
        DP  = DP  - (BB1(II)*FV(II,1)+BB2(II)*FV(II,2))
        FFP = FFP - (BB1(II)*FV(II,3)+BB2(II)*FV(II,4))
      END DO

C     Solution "loop":

      IF (DABS(DP).LT.1.D-10)  THEN
c        WRITE(MTERM,*)'ERROR in ELUPDT: DP is nearly zero'
        RETURN
      ENDIF
      
C     At first we calculate the pressure on element IEL,
C     which is simply given by multiplying FFP with the
C     inverse "matrix" DP, i.e.:
      
      PP = FFP/DP
      
C     With the help of the pressure, calculate the velocity.
C     This can be done again by the Schur-Complement approach using
C
C           u = S^-1 (f-Bp)
C
C     locally on the current cell:
      
      DO II=1,4
        UU(II,1) = (FF1(II)-DD1(II)*PP)
        UU(II,2) = (FF2(II)-DD2(II)*PP)
      END DO
      CALL DGETRS( 'N', 4, 2, AI, 4, IPIV, UU, 4, I )

C     Update the global solution vector (U1,U2,P).
C     Overwrite the pressure on element IEL and the
C     velocity in the four DOF's on element IEL by the
C     new pressure/velocity values in these nodes:

      P(IEL)=PP
      
      DO II = 1,4

C       Don't overwrite velocities on Dirichlet-boundary vertices!

        IF (.NOT.BDBC(II)) THEN
          I=IU(II)
          U1(I)=UU(II,1)
          U2(I)=UU(II,2)
        END IF
        
      END DO

      END

