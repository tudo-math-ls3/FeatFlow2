C UMFPACK demo program.
C
C Factor and solve a 5-by-5 system, Ax=b, using default parameters,
C except with complete printing of all arguments on input and output,
C where

C     [ 2+i  3+2i  0    0  0    ]      [  8+5i  ]                  [ 1 ]
C     [ 3-i  0     4-i  0  6+3i ]      [ 45+11i ]                  [ 2 ]
C A = [ 0   -1+i  -3+2i 2  0    ], b = [ -3+8i  ]. Solution is x = [ 3 ]
C     [ 0    0     1-i  0  0    ]      [  3-3i  ]                  [ 4 ]
C     [ 0    4     2    0  1+9i ]      [ 19+45i ]                  [ 5 ]
C
C Next, modify one entry (A (5,2) = 1.0D-14) and
C refactorize.  Solve Ax=b both without and with iterative refinement,
C with true solution (rounded to 4 digits past the decimal point):
C
C          [  -0.4634 + 1.1707i ]
C          [   2.8105 - 0.8330i ]
C          [   3.0000 - 0.0000i ]
C          [   3.9887 - 0.8218i ]
C          [   5.0976 - 0.8780i ]

        PROGRAM ZMAIN
        INTEGER NMAX, NEMAX, LVALUE, LINDEX
        PARAMETER (NMAX=20, NEMAX=100, LVALUE=300, LINDEX=300)
        INTEGER KEEP (20), INDEX (LINDEX), INFO (40),
     $     I, ICNTL (20), N, NE, AI (2*NEMAX)
        COMPLEX*16 
     $     B (NMAX), X (NMAX), W (4*NMAX), VALUE (LVALUE), AX (NEMAX)
        DOUBLE PRECISION
     $     CNTL (10), RINFO (20)

C Read input matrix and right-hand side.  Keep a copy of the triplet
C form in AI and AX.

        READ (5, *) N, NE
        READ (5, *) (AI (I), AI (NE+I), I = 1,NE)
        READ (5, 1) (AX (I), I = 1,NE)
1       FORMAT (2F5.1)
        READ (5, 1) (B (I), I = 1,N)
        DO 10 I = 1, NE
           INDEX (I) = AI (I)
           INDEX (NE+I) = AI (NE+I)
           VALUE (I) = AX (I)
10      CONTINUE

C Initialize controls, and change default printing control.  Note that
C this change from the default should only be used for test cases.  It
C can generate a lot of output for large matrices. 

        CALL UMZ21I (KEEP, CNTL, ICNTL)
        ICNTL (3) = 4

C Factorize A, and print the factors.  Input matrix is not preserved.

        CALL UMZ2FA (N, NE, 0, .FALSE., LVALUE, LINDEX, VALUE, INDEX,
     $               KEEP, CNTL, ICNTL, INFO, RINFO)
        IF (INFO (1) .LT. 0) STOP

C Reset default printing control (UMZ21I could be called instead)
        ICNTL (3) = 2

C Solve Ax = b and print solution.

        CALL UMZ2SO (N, 0, .FALSE., LVALUE, LINDEX, VALUE, INDEX,
     $               KEEP, B, X, W, CNTL, ICNTL, INFO, RINFO)
        WRITE (6, *) 'Solution to Ax=b:'
        WRITE (6, 2) (X (I), I = 1, N)
2       FORMAT (' real: ', E20.12, ' imag: ', E20.12)
        IF (INFO (1) .LT. 0) STOP

C Modify one entry of A, and refactorize using UMZ2RF.

        DO 20 I = 1, NE
           INDEX (I) = AI (I)
           INDEX (NE+I) = AI (NE+I)
           VALUE (I) = AX (I)
20      CONTINUE
C       A (5,2) happens to be (PAQ)_22, the second pivot entry:
        VALUE (10) = 1.0D-14

        CALL UMZ2RF (N, NE, 1, .FALSE., LVALUE, LINDEX, VALUE, INDEX,
     $               KEEP, CNTL, ICNTL, INFO, RINFO)
        IF (INFO (1) .LT. 0) STOP

C Solve Ax = b without iterative refinement, and print solution.
C This will be very inaccurate due to the tiny second pivot entry.

        CALL UMZ2SO (N, 0, .FALSE., LVALUE, LINDEX,  VALUE, INDEX,
     $               KEEP, B, X, W, CNTL, ICNTL, INFO, RINFO)
        WRITE (6, *) 'Solution to modified Ax=b, no iter. refinement:'
        WRITE (6, 2) (X (I), I = 1, N)
        IF (INFO (1) .LT. 0) STOP

C Solve Ax = b with iterative refinement, and print solution.
C This is much more accurate.

        ICNTL (8) = 10
        CALL UMZ2SO (N, 0, .FALSE., LVALUE, LINDEX,  VALUE, INDEX,
     $               KEEP, B, X, W, CNTL, ICNTL, INFO, RINFO)
        WRITE (6, *) 'Solution to modified Ax=b, with iter. refinement:'
        WRITE (6, 2) (X (I), I = 1, N)
        IF (INFO (1) .LT. 0) STOP
        STOP
        END
