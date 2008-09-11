************************************************************************
* Classical VANCA smoother, SCGS method
*
* Extended calling convention
*
* This routine performs one smoothing step
*
*    x_new = x  +  RLXSM ( Vanca(x) - x )
*
* with Vanca(x) being a better solution than X to the system Ax=b
* calculated by a local Gauss-Seidel approach.
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

      SUBROUTINE  VANCS2 (U1,U2,P,U1OLD,U2OLD,POLD,F1,F2,FP,
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
      DOUBLE PRECISION RLXSM
      
C     Local arrays for informations about one element
      
      DOUBLE PRECISION AA(4),BB1(4),BB2(4),FF1(4),FF2(4)
      INTEGER IU(4)
      LOGICAL BDBC(4)
      
C     local variables

      INTEGER IEL,II,IMID,IA1,IA2,IA,IB,I,J,IB1,IB2,JP1,JP2
      DOUBLE PRECISION FFP,PJP,AOFF,A1,A2

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

        DO II=1,4

C         Set I to the DOF that belongs to our edge II:

          IMID = KMID(II,IEL)
          I    = IMID-NVT
          
C         Write the number of the edge/node to IU:
          
          IU(II) = I
          
C         Write into BDBC whether the edge is a Dirichlet-
C         Boundary edge:
          
          BDBC(II) = IAND(KXNPR(1,IMID),2**12).NE.0

C         Put on AA(.) the diagonal entry of matrix A

          IA1    = KLDA(I)
          AA(II) = A(IA1)

C         Set FF1/FF2 initially to the value F1/F2 of the right hand
C         side vector that belongs to our current DOF corresponding
C         to II

          FF1(II) = F1(I)
          FF2(II) = F2(I)
          
C         What do we have at this point?
C         FF1/FF2: "local" RHS vector belonging to the DOF's on the
C                  current element
C         AA     : Diagonal entries of A belonging to these DOF's
C         BDBC   : Whether the DOF's belong to Dirichlet nodes
C
C         And at the moment:
C         I      : number of current DOF on element IEL
C         II     : "local" number of DOF on element IEL, i.e.
C                  number of the edge
C
C         Now comes the crucial point with the "update": How to
C         subsequently update the vertex values, such that the whole
C         thing still converges to the solution, even if a node
C         is updated more than once? Here, we use a typical
C         matrix-decomposition approach:
C
C         Again consider the problem:
C         
C            [ A   B ] (u) = (f)
C            [ B^t 0 ] (p)   (g)
C     
C         We assume, that all components in the vector (u,p) are
C         given - except for the four velocity unknowns and the
C         pressure unknown on the current element; these five unknowns
C         are located anywhere in the (u,p) vector. The idea is to
C         shift "everything known" to the right hand side to obtain
C         a system for only these unknowns!
C
C         Extracting all the lines of the system that correspond to
C         DOF's on our single element IEL results in a rectangular
C         system of the form
C
C            [ === A^ === B~ ] (|) = (f)   
C            [ B~^t       0  ] (u)   (g)   
C                              (|)
C                              (p)
C
C         <=>  [ MAT1 MAT2 B~ ] (u) = (fu)
C              [    ****   0  ] (p)   (fp)
C
C         with A^ being an 4 x (2*NVT) matrix for the two velocity
C         components and B being an (2*4) x 1 matrix that couples the
C         velocities to the pressure on our current element.
C         B~ is a 8 x 2 matrix: As every velocity couples with at most
C         two pressure elements on the neighbour cell, so we have
C         2 columns in the B-matrix.
C
C                IEL                              IEL
C             |--------|             |--------|--------|
C             |        |             |        |        |
C             |   P    |      or     |   Q    X   P    |
C             |        |             |        |        |
C           --|---X----|--           |--------|--------|
C
C         Now, throw all summands to the RHS vector that do
C         not correspond to the velocity/pressure DOF's on our single 
C         element IEL!
C
C         That way, A^ is reduced to a square matrix with two square
C         submatrices A~ of size 4 x 4. The 8 x 2-matrix B~ reduces to 
C         two 4 x 1 submatrices (originally, every velocity couples with
C         two pressure elements on the neighbour cell, so we have
C         2 columns in the B-matrix).
C
C         The system roughly looks like:
C
C         [ ....        B1 ] (u1)   (f1)           (| )           (| )           (| )
C         [ .A~.        |  ] (u2)   (| )           (| )           (| )           (| )
C         [ ....        |  ] (u3)   (| )           (| )           (| )           (| )
C         [ ....        |  ] (u4)   (| )           (| )           (| )           (| )
C         [       ....  B2 ] (v1) = (f2) - sum u_i (ci) - sum v_i (di) - sum   Q (qi)
C         [       .A~.  |  ] (v2)   (| )   i>4     (| )   i>4     (| )  Q in N   (| )
C         [       ....  |  ] (v3)   (| )           (| )           (| )           (| )
C         [       ....  |  ] (v4)   (| )           (| )           (| )           (| )
C         [ B1^T  B2^T  0  ] (p )   (fp)           (0 )           (0 )           (0 )
C
C         assuming in the example that element IEL contains u1,u2,u3,u4,
C         v1,v2,v3,v4 and p as DOF's, and using the notation that ci 
C         and di are the i th column of matrix MAT1 and MAT2 above,
C         respectively. 
C         Q is the pressure on the neighbour elements in the 
C         neighbourhood N of element IEL and
C         qi the corresponding column from the B~-matrix.
C
C         But we do even more! Consider the 4x4-matrix A~. This can be
C         decomposed into A~ = L+S+R with L being the lower left, R 
C         being the upper right and S being the diagonal part of A~.
C         Again, throw Lu, Ru and Lv, Rv to the right hand side vector,
C         leaving alone the diagonal 4x4 matrix S on the left:
C        
C         [ *           B1 ] (u1)   (d11)
C         [  *          |  ] (u2)   (d12)
C         [   *         |  ] (u3)   (d13)
C         [    *        |  ] (u4)   (d14)       [ S         B1 ] (ul)   (du)
C         [       *     B2 ] (v1) = (d21)  <=>  [      S    B2 ] (vl) = (dv)
C         [        *    |  ] (v2)   (d22)       [ B1^t B2^t  0 ] (p )   (dp)
C         [         *   |  ] (v3)   (d23)
C         [          *  |  ] (v4)   (d24)
C         [ B1^T  B2^T  0  ] (p )   (d33)
C
C         with d1k = f1k - sum  A~(k,i) u_i.
C                          i<>k
C         d2k is calculated similar. (du, dv, dp) is called 
C         "local residuum".
C
C         Note: This only holds for inner edges!                                 
C         Do we have a boundary edge?

          IF (.NOT.BDBC(II)) THEN

C           No, we have an inner edge.
C
C           Create a "local" residuum!
C           Take the original RHS and subtract all entries of "A*u" to
C           this RHS as being given. What we throw here to the RHS is
C           therefore on one hand a set of entries of the old velocities
C           as well as a set of entries which might have been calculated
C           by a previous update of another cell in our neighbourhood.

            IA2=KLDA(I+1)-1
            DO IA=IA1+1,IA2
              J=KCOLA(IA)
              AOFF=A(IA)
              FF1(II) = FF1(II)-AOFF*U1(J)
              FF2(II) = FF2(II)-AOFF*U2(J)
            END DO

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

            IB1=KLDB(I)
            IB2=KLDB(I+1)-1
            JP1=KCOLB(IB1)
            JP2=KCOLB(IB2)
            
C           Check which of these two belong to the current element.
C           The corresponding other one is to be subtracted from the RHS
C           vector.

            IF (JP1.EQ.IEL)  THEN
            
              PJP = P(JP2)
              BB1(II) = B1(IB1)
              BB2(II) = B2(IB1)
              
C             Don't subtract anything from the RHS if there's
C             no neighbour element!

              IF (IB1.NE.IB2) THEN
                FF1(II) = FF1(II)-B1(IB2)*PJP
                FF2(II) = FF2(II)-B2(IB2)*PJP
              END IF
              
            ELSE IF (JP2.EQ.IEL)  THEN
            
              PJP = P(JP1)
              BB1(II) = B1(IB2)
              BB2(II) = B2(IB2)
              
C             Don't subtract anything from the RHS if there's
C             no neighbour element!

              IF (IB1.NE.IB2) THEN
                FF1(II) = FF1(II)-B1(IB1)*PJP
                FF2(II) = FF2(II)-B2(IB1)*PJP
              END IF
              
            ELSE
              WRITE (MTERM,*) 
     *              'ERROR in SMOOTH: IEL entry in B not found'
              STOP
            ENDIF
            
C           That's it, FF1/FF2 form the "local defect".

          ELSE

C           Hmmm, the edge II is a Dirichlet boundary edge.
C           That means there is no adjacent element with a pressure that
C           we have to subtract from the RHS vector.
C           Simply write the one and only element in the current line
C           of the B-matrices to BB1/BB2.

            IB = KLDB(I)
            BB1(II) = B1(IB)
            BB2(II) = B2(IB)
            
          END IF ! (not BDBC(II))

C         Go to the next DOF on our current element IEL...

        END DO ! II

C       Call the subroutine ELUPD2 for element update.
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

        CALL  ELUPD2 (U1,U2,P,IEL,IU,BDBC,AA,BB1,BB2,FF1,FF2,FFP)

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
* Classical VANCA solver, SCGS method
*
* Extended calling convention
*
* This routine performs one fix point iteration step
*
*    x_new = x  +  RLXSM ( Vanca(x) - x )
*
* with Vanca(x) being a better solution than X to the system Ax=b
* calculated by a local Gauss-Seidel approach.
* This is basically the same as the Vanca-smoother in VANCAS above,
* but in difference the routine returns information about
* the maximum velocity/divergence defect to the caller.
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

      SUBROUTINE  VANCE2 (U1,U2,P,U1OLD,U2OLD,POLD,F1,F2,FP,
     *                    A,KCOLA,KLDA,B1,B2,KCOLB,KLDB,NU,NP,
     *                    KMID,KXNPR,NEL,NVT,RLXPAR,DMAXU,DMAXP)

      IMPLICIT NONE
      
C common blocks      
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      
C parameters

      DOUBLE PRECISION U1(*),U2(*),P(*),U1OLD(*),U2OLD(*),POLD(*)
      DOUBLE PRECISION F1(*),F2(*),FP(*)
      DOUBLE PRECISION RLXPAR,DMAXP,FFPH
      DOUBLE PRECISION A(*),B1(*),B2(*)
      INTEGER KCOLA(*),KLDA(*),KCOLB(*),KLDB(*)
      INTEGER KMID(NNVE,*),KXNPR(2,*)
      INTEGER NU,NP,NEL,NVT
      
C     Local arrays for informations about one element
      
      DOUBLE PRECISION AA(4),BB1(4),BB2(4),FF1(4),FF2(4)
      INTEGER IU(4)
      LOGICAL BDBC(4)
      
C     local variables

      INTEGER IEL,II,IMID,IA1,IA2,IA,IB,I,J,IB1,IB2,JP,JP1,JP2
      INTEGER IB3,JP3,ICOLA1
      DOUBLE PRECISION FFP,PJP,AOFF,A1,A2,DMAXU,AH

C     DMAXU/DMAXP will be build by taking a maximum; initialize
C     with 0.

      DMAXU=0D0
      DMAXP=0D0

C     Copy U->UOLD; we need it later for the relaxation step
C     Unew = Uold + RLXSM (U-Uold).

      CALL  LCP1 (U1,U1OLD,NU)
      CALL  LCP1 (U2,U2OLD,NU)
      CALL  LCP1 (P ,POLD ,NP)
C
C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================

C     Loop over all elements to (re-)calculate the DOF's in velocity and
C     pressure on each element.

      DO IEL=1,NEL
      
C       Fetch the right hand side of the pressure in the equation
C       to FFP.
      
        FFP =FP(IEL)
        FFPH=FFP

C       Loop over all 4 U-nodes of that element

        DO II=1,4

          IMID=KMID(II,IEL)
          I=IMID-NVT
          
C         Write the global number of DOF into IU:
          
          IU(II)=I
          
C         Is this DOF of Dirichlet type?
          
          BDBC(II) = IAND(KXNPR(1,IMID),2**12).NE.0
          
C         We treat the current solution vector (u,p) as
C         given except for the five unknowns on the current element.
C         Therefore we split the system matrix for the velocities
C         A into A=L+S+U with S being diagonal.
C         In the velocity system Au=f <=> (L+S+U)u=f we shift everything
C         given to the right hand side, resulting in Su=f-(L+U)u
C         with S being the 4x4 diagonal matrix for the four velocity 
C         unknowns on the current element, with a_ii on the main
C         diagonal.
C
C         Put on AA(.) the diagonal entry of matrix A, which forms the
C         matrix S:

          IA1 = KLDA(I)
          AA(II) = A(IA1)
 
C         Form the RHS f-(L+U)u of the system.
C         Initial setting of FF1,FF2 is the right hand side f:

          FF1(II) = F1(I)
          FF2(II) = F2(I)

C         Check if the current node is a boundary node

          IF (.NOT.BDBC(II)) THEN

C           No, the edge is an inner edge or a Neumann boundary edge:
C
C                          IEL               IEL      
C             |--------|--------|         |--------|  
C             |        |        |         |        |  
C             |   Q    X   P    |         |   P1   |  
C             |        |        |         |        |  
C             |--------|--------|       --|---X----|--Neumann
C
C           Put all off-diagonal entries of A to the right hand side,
C           i.e. subtract (L+U)u (with L,U not being only the lower
C           left/upper right part of a triangular matrix, but the lower
C           left/upper right part of the rectangular matrix that arises
C           when extracting only the lines from the matrix that
C           correspond to our four DOF's!

            IA2=KLDA(I+1)-1
            DO IA=IA1+1,IA2
              J=KCOLA(IA)
              AOFF=A(IA)
              FF1(II)=FF1(II)-AOFF*U1(J)
              FF2(II)=FF2(II)-AOFF*U2(J)
            END DO

C           Get BB1,BB2 and modify FF1,FF2
C
C           The velocity part is subtracted from f, but that's not
C           enough. As the main matrix has the form
C             [ A         B1 ]
C             [      A    B2 ]
C             [ B1^t B2^t 0  ],
C           shifting everything given to the RHS vector after extracting 
C           only the crucial lines also makes it necessary to subtract
C           something regarding the pressure.
C           Our current element is connected to a neighbour element, as
C           we are looking to an inner edge. Get the pressure Q on the
C           neighbour elements and subtract Bx^t*Q from the RHS.
C
C           In case our edge is Neumann-boundary type, there's only one
C           entry in the B-Matrix, and thus the first and last entry of
C           the row in the B-matrix are the same. We fetch the
C           first and last entry of the line, which will coincide
C           then (IB1=IB2). As there is no neighbour cell, we must
C           not subtract anything from the RHS! Nevertheless we have to
C           note the column numbers in the B-matrices to our BB1/BB2
C           vectors.
C
C           Write the two 4x1-matrices of the pressure P on the current
C           element to BB1 and BB2.
C
C           IB1 = Index of 1st B-entry
C           IB2 = Index of 2nd B-entry (=IB1 if Neumann bdry. edge)
C           JP1 = Column of 1st B-entry
C           JP2 = Column of 2nd B-entry (=JP1 if Neumann bdry. edge)

            IB1=KLDB(I)
            IB2=KLDB(I+1)-1
            JP1=KCOLB(IB1)
            JP2=KCOLB(IB2)

            IF (JP1.EQ.IEL)  THEN
              PJP=P(JP2)
              BB1(II)=B1(IB1)
              BB2(II)=B2(IB1)
              
C             Don't subtract anything from the RHS if there's
C             no neighbour element!

              IF (IB1.NE.IB2) THEN
                FF1(II)=FF1(II)-B1(IB2)*PJP
                FF2(II)=FF2(II)-B2(IB2)*PJP
              END IF
              
              IB3    =IB1
              JP3    =JP1
              
            ELSE IF (JP2.EQ.IEL)  THEN
            
              PJP=P(JP1)
              BB1(II)=B1(IB2)
              BB2(II)=B2(IB2)
              
C             Don't subtract anything from the RHS if there's
C             no neighbour element!

              IF (IB1.NE.IB2) THEN
                FF1(II)=FF1(II)-B1(IB1)*PJP
                FF2(II)=FF2(II)-B2(IB1)*PJP
              END IF
              
              IB3    =IB2
              JP3    =JP2
              
            ELSE
              WRITE(MTERM,*) 
     *             'ERROR in VANCAE: IEL entry in B not found'
              STOP
            END IF

C           Additionally to IBx/JPx we now have:
C
C           PJP = pressure on the neighbour cell
C           IB3 = the index of the B-entry of the neighbour cell
C           JP3 = the column number of the B-entry of the neighbour cell
C
C           FF1/FF2 currently shows the value of 
C                f-(L+D)u-(neighbouring pressure)*B
C           By subtracting also (S*u) and the pressure contribution 
C           on the current cell, we form the residuum of the
C           velocity. The maximum of the residuum is stored to DMAXU.

            AH=A(IA1)
            ICOLA1=KCOLA(IA1)
            
C           The column number of the B-matrices corresponds to the
C           DOF in the pressure P. JP3 is the column number of
C           the entry in the B-matrix corresponding to the current
C           cell. By accessing P(JP3), fetch the pressure on the
C           current cell to PJP:

            PJP=P(JP3)
            
            DMAXU = MAX(DMAXU,ABS(FF1(II)-AH*U1(ICOLA1)-B1(IB3)*PJP))
            DMAXU = MAX(DMAXU,ABS(FF2(II)-AH*U2(ICOLA1)-B2(IB3)*PJP))

C           FFPH on the other hand is the original RHS of the divergence
C           equation B^t*u=0 (was set that way above). Subtract B^t*u
C           from that to form the "divergence defect" g-B^t*u on the
C           current cell.

            FFPH = FFPH - BB1(II)*U1(IU(II))-BB2(II)*U2(IU(II))

          ELSE

C           II is a Dirichlet boundary node. There's for sure no
C           neighbour cell attached, so don't have to subtract something
C           like a pressure part from the RHS.
C           Simply get the entries of the two 4x1 B-matrices for the
C           current element.

            IB = KLDB(I)
            JP = KCOLB(IB)
            BB1(II) = B1(IB)
            BB2(II) = B2(IB)

C           FFPH is the original RHS of the divergence equation
C           B^t*u=0 (was set that way above). Subtract B^t*u
C           from that to form the "divergence defect" g-B^t*u on the
C           current cell.

            FFPH = FFPH - BB1(II)*U1(IU(II))-BB2(II)*U2(IU(II))
          
          END IF ! (not BDBC(II))

C         Fetch the next DOF on the current element...

        END DO ! II

C       Build the maximum "divergence defect" g-B^t*u in DMAXP.
C       REMARK: This is a difference to the original implementation
C       of VANCA, which does not use an ABS() here!

        DMAXP=MAX(DMAXP,ABS(FFPH))

C       Call the subroutine ELUPD2 for element update.
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

        CALL  ELUPD2 (U1,U2,P,IEL,IU,BDBC,AA,BB1,BB2,FF1,FF2,FFP)

C       Improve the DOF's on the next element...

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
        A1 = 1D0-RLXPAR
        A2 =     RLXPAR
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
*   AA     - array [1..4] of double
*            Diagonal elements of the system matrix corresponding
*            to the DOF's in IU(.)
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

      SUBROUTINE ELUPD2 (U1,U2,P,IEL,IU,BDBC,AA,BB1,BB2,FF1,FF2,FFP)

      IMPLICIT NONE

      INCLUDE 'cout.inc'

      DOUBLE PRECISION U1(*),U2(*),P(*)

C     parameters

      DOUBLE PRECISION AA(4),BB1(4),BB2(4),FF1(4),FF2(4)
      DOUBLE PRECISION AI(4),DD1(4),DD2(4),UU1(4),UU2(4)
      LOGICAL BDBC(4)
      INTEGER IU(4)

C local variables

      INTEGER II, IEL, I
      DOUBLE PRECISION DP, FFP, PP

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
C     Here S is given as a diagonal matrix with diagonal entries of
C     the original system matrix. We can write:
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
C     is then also trivial as S (and thus S^-1) is a diagonal matrix.

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
        
        IF (DABS(AA(II)).LT.1.D-10)  THEN
C          WRITE(MTERM,*)'ERROR in ELUPD2: diagonal entry is nearly zero'
C          WRITE(MTERM,*)'II, I, BNDR(II) : ',II,IU(II),BDBC(II)
          RETURN
        ENDIF
        
C       AI(.) saves the diagonal matrix S^-1:
    
        AI(II)=1D0/AA(II)
        
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
C     b) S^-1 is a diagonal matrix, given by the 4 diagonal entries of A
C     c) B^t S^-1 B is therefore a 1x1 matrix, thus a scalar
C
C     So in the factorization loop we can calculate:
C
C       DP  = - (B^T S^-1 B)
C       FFP = g - (B^T S^-1 f)
C
C     As S and S^-1 are a diagonal matrices, we can exploit
C     B^T S^-1  =  S^-1 B^T  which saves some multiplications...

      DP=0D0
      
      DO II = 1,4
        DP  = DP  - AI(II)*(BB1(II)*DD1(II)+BB2(II)*DD2(II))
        FFP = FFP - AI(II)*(BB1(II)*FF1(II)+BB2(II)*FF2(II))
      END DO

C     Solution "loop":

      IF (DABS(DP).LT.1.D-10)  THEN
c        WRITE(MTERM,*)'ERROR in ELUPD2: DP is nearly zero'
        RETURN
      ENDIF
      
C     At first we calculate the pressure on element IEL,
C     which is simply given by multiplying FFP with the
C     inverte "matrix" DP, i.e.:
      
      PP = FFP/DP
      
C     With the help of the pressure, calculate the velocity.
C     This can be done again by the Schur-Complement approach using
C
C           u = S^-1 (f-Bp)
C
C     locally on the current cell:
      
      DO II=1,4
        UU1(II) = AI(II)*(FF1(II)-DD1(II)*PP)
        UU2(II) = AI(II)*(FF2(II)-DD2(II)*PP)
      END DO

C     Update the global solution vector (U1,U2,P).
C     Overwrite the pressure on element IEL and the
C     velocity in the four DOF's on element IEL by the
C     new pressure/velocity values in these nodes:

      P(IEL)=PP
      
      DO II = 1,4

C       Don't overwrite velocities on Dirichlet-boundary vertices!

        IF (.NOT.BDBC(II)) THEN
          I=IU(II)
          U1(I)=UU1(II)
          U2(I)=UU2(II)
        END IF
        
      END DO

      END

