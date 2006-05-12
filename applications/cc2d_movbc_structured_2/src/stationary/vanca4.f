***********************************************************************
* Classical VANCA preconditioner, SCGS method
*
* Extended calling convention
*
* This routine performs one preconditioning step
*
*    (du,dp) = Vanca(F)
*
* with the Diagonal-Vanca preconditioner. (du,dp) must be =0 initially.
* VANCP4 executes a local Block-Gauss-Seidel preconditioner onto the
* vector F to calculate a preconditioned vector (du,dp).
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

      SUBROUTINE  VANCP4 (U1,U2,P,F1,F2,FP,
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
      
C     Local arrays for informations about one element
      
      DOUBLE PRECISION AA(4),BB1(4),BB2(4),FF1(4),FF2(4)
      INTEGER IU(4)
      LOGICAL BDBC(4)
      
C     local variables

      INTEGER IEL,II,IMID,IA1,IA2,IA,IB,I,J,IB1,IB2,JP1,JP2
      DOUBLE PRECISION FFP,PJP,AOFF,A1,A2

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
     *              'ERROR in VANCP4: IEL entry in B not found'
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
C     At this point we are finished calculating a the preconditioned
C     vector (du,dv,dp). The caller can use it e.g. in a defect
C     correction approach...

      END
