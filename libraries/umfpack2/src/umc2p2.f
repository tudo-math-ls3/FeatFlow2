        SUBROUTINE UMC2P2 (WHO, ERROR, I, J, X, IO)
        INTEGER WHO, ERROR, I, J, IO
        COMPLEX
     $          X
        
C=== UMC2P2 ============================================================
C
C  Unsymmetric-pattern MultiFrontal Package (UMFPACK). Version 2.2c
C  Copyright (C) 1997, Timothy A. Davis, University of Florida, USA.
C  ALL RIGHTS RESERVED.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  July 7, 1997. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C***********************************************************************
C* NOTICE:  "The UMFPACK Package may be used SOLELY for educational,   *
C* research, and benchmarking purposes by non-profit organizations and *
C* the U.S. government.  Commercial and other organizations may make   *
C* use of UMFPACK SOLELY for benchmarking purposes only.  UMFPACK may  *
C* be modified by or on behalf of the User for such use but at no time *
C* shall UMFPACK or any such modified version of UMFPACK become the    *
C* property of the User.  UMFPACK is provided without warranty of any  *
C* kind, either expressed or implied.  Neither the Authors nor their   *
C* employers shall be liable for any direct or consequential loss or   *
C* damage whatsoever arising out of the use or misuse of UMFPACK by    *
C* the User.  UMFPACK must not be sold.  You may make copies of        *
C* UMFPACK, but this NOTICE and the Copyright notice must appear in    *
C* all copies.  Any other use of UMFPACK requires written permission.  *
C* Your use of UMFPACK is an implicit agreement to these conditions."  *
C*                                                                     *
C* The MA38 Package in Release 12 of the Harwell Subroutine Library    *
C* (HSL) has equivalent functionality (and identical calling interface)*
C* as UMFPACK (the HSL has single and double precision versions only,  *
C* however).  It is available for commercial use.   Technical reports, *
C* information on HSL, and matrices are available via the World Wide   *
C* Web at http://www.cis.rl.ac.uk/struct/ARCD/NUM.html, or by          *
C* anonymous ftp at seamus.cc.rl.ac.uk/pub.  Also contact Dr. Scott    *
C* Roberts, Harwell Subroutine Library, B 552, AEA Technology,         *
C* Harwell, Didcot, Oxon OX11 0RA, England.                            *
C* telephone (44) 1235 434988, fax (44) 1235 434136                    *
C* email Scott.Roberts@aeat.co.uk, who will provide details of price   *
C* and conditions of use.                                              *
C***********************************************************************

C=======================================================================
C  NOT USER-CALLABLE

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Print error and warning messages for for UMC2FA, UMC2RF, and UMC2SO.

C=======================================================================
C  INSTALLATION NOTE:
C=======================================================================
C
C  This routine can be deleted on installation (replaced with a dummy
C  routine that just returns without printing) in order to completely
C  disable the printing of all error and warning messages.  The error
C  and warning return flag (Info (1)) will not be affected.  To
C  completely disable all I/O, you can also replace the UMC2P1 routine
C  with a dummy subroutine.  If you make this modification, please do
C  not delete any original code - just comment it out instead.  Add a
C  comment and date to your modifications.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       who:            what user-callable routine called UMC2P2:
C                       1: UMC2FA, 2: UMC2RF, 3: UMC2SO
C       i, j, x:        the relevant offending value(s)
C       io:             I/O unit on which to print.  No printing
C                       occurs if < 0.
C       error:          the applicable error (<0) or warning (>0)
C                       Errors (<0) cause the factorization/solve to
C                       be terminated.  If an error occurs, a prior
C                       warning status is overwritten with the error
C                       status.
C
C  The following error codes are returned in Info (1) by UMC2ER.
C  These errors cause the factorization or solve to terminate:
C
C  Where**      Error   Description
C 
C  FA RF  -     -1      N < 1
C  FA RF  -     -2      NE < 1 or NE > maximum value
C  FA RF  -     -3      LINDEX too small
C  FA RF  -     -4      LVALUE too small
C  FA RF  -     -5      both LINDEX and LVALUE are too small
C   - RF  -     -6      prior pivot ordering no longer acceptable
C   - RF SO     -7      LU factors are uncomputed, or are corrupted

C  FA RF SO     -8      transposed or complex conj. transposed systems
C                       are not factorized or solved.

C
C  The following warning codes are returned in Info (1) by UMC2ER.
C  The factorization or solve was able to complete:
C
C  FA RF  -     1       invalid entries
C  FA RF  -     2       duplicate entries
C  FA RF  -     3       invalid and duplicate entries
C  FA RF  -     4       singular matrix
C  FA RF  -     5       invalid entries, singular matrix
C  FA RF  -     6       duplicate entries, singular matrix
C  FA RF  -     7       invalid and duplicate entries, singular matrix
C   -  - SO     8       iterative refinement cannot be done
C
C  The following are internal error codes (not returned in Info (1))
C  for printing specific invalid or duplicate entries.  These codes are
C  for UMC2CO, UMC2OF, and UMC2R2.  Invalid entries are ignored, and
C  duplicate entries are added together (and the factorization
C  continues).  Warning levels (1..7) will be set later by UMC2ER,
C  above.
C
C  FA RF  -     99      invalid entry, out of range 1..N
C  FA RF  -     98      duplicate entry
C   - RF  -     97      invalid entry:  within a diagonal block, but not
C                       in the pattern of the LU factors of that block.
C   - RF  -     96      invalid entry:  below the diagonal blocks.  Can
C                       only occur if the matrix has been ordered into
C                       block-upper-triangular form.
C   - RF  -     95      invalid entry:  matrix is singular.  The
C                       remaining rank 0 submatrix yet to be factorized
C                       is replaced with the identity matrix in the LU
C                       factors.  Any entry that remains is ignored.

C ** FA: UMC2FA, RF: UMC2RF, SO: UMC2SO

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C  Error or warning message printed on I/O unit

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutines:  UMC2ER, UMC2CO, UMC2R0, UMC2R2

C=======================================================================
C  EXECUTABLE STATEMENTS:
C       if (printing disabled on installation) return
C=======================================================================

        IF (IO .LT. 0) THEN 
C          printing of error / warning messages has not been requested
           RETURN
        ENDIF 

        IF (WHO .EQ. 1) THEN 

C          -------------------------------------------------------------
C          UMC2FA error messages
C          -------------------------------------------------------------

           IF (ERROR .EQ. -1) THEN 
              WRITE (IO, 1) 'UMC2FA: N less than one!'
           ELSE IF (ERROR .EQ. -2) THEN 
              WRITE (IO, 1) 'UMC2FA: NE less than one!'
           ELSE IF (ERROR .EQ. -3) THEN 
              WRITE (IO, 1)
     $        'UMC2FA: LINDEX too small!  Must be greater than ', I
           ELSE IF (ERROR .EQ. -4) THEN 
              WRITE (IO, 1)
     $        'UMC2FA: LVALUE too small!  Must be greater than ', I

C          -------------------------------------------------------------
C          UMC2FA cumulative warning messages
C          -------------------------------------------------------------

           ELSE IF (ERROR .EQ. 1) THEN 
              WRITE (IO, 1) 'UMC2FA: ', I,
     $        ' invalid entries ignored (out of range 1..N).'
           ELSE IF (ERROR .EQ. 2) THEN 
              WRITE (IO, 1) 'UMC2FA: ', I,' duplicate entries summed.'
           ELSE IF (ERROR .EQ. 4) THEN 
              WRITE (IO, 1)
     $        'UMC2FA: matrix is singular.  Only ', I, ' pivots found.'

C          -------------------------------------------------------------
C          UMC2FA non-cumulative warning messages (internal error codes)
C          -------------------------------------------------------------

           ELSE IF (ERROR .EQ. 99) THEN 
              WRITE (IO, 2)
     $        'UMC2FA: invalid entry (out of range 1..N):', I, J, X
           ELSE IF (ERROR .EQ. 98) THEN 
              WRITE (IO, 2)
     $        'UMC2FA: duplicate entry summed:', I, J, X

           ELSE IF (ERROR .EQ. -8) THEN 
              WRITE (IO, 1) 'UMC2FA: transa argument must be .false.!'

           ENDIF 

        ELSE IF (WHO .EQ. 2) THEN 

C          -------------------------------------------------------------
C          UMC2RF error messages
C          -------------------------------------------------------------

           IF (ERROR .EQ. -1) THEN 
              WRITE (IO, 1) 'UMC2RF: N less than one!'
           ELSE IF (ERROR .EQ. -2) THEN 
              IF (I .LT. 0) THEN 
                 WRITE (IO, 1) 'UMC2RF: NE less than one!'
              ELSE 
                 WRITE (IO, 1)
     $           'UMC2RF: NE too large!  Must be less than ', I
              ENDIF 
           ELSE IF (ERROR .EQ. -3) THEN 
              WRITE (IO, 1)
     $        'UMC2RF: LINDEX too small!  Must be greater than ', I
           ELSE IF (ERROR .EQ. -4) THEN 
              WRITE (IO, 1)
     $        'UMC2RF: LVALUE too small!  Must be greater than ', I
           ELSE IF (ERROR .EQ. -6) THEN 
              WRITE (IO, 1) 'UMC2RF: pivot order from UMC2FA failed!'
           ELSE IF (ERROR .EQ. -7) THEN 
              WRITE (IO, 1)
     $        'UMC2RF: LU factors uncomputed or corrupted!'

C          -------------------------------------------------------------
C          UMC2RF cumulative warning messages
C          -------------------------------------------------------------

           ELSE IF (ERROR .EQ. 1) THEN 
              IF (I .GT. 0) THEN 
                 WRITE (IO, 1) 'UMC2RF: ', I,
     $           ' invalid entries ignored (out of range 1..N).'
              ELSE 
                 WRITE (IO, 1) 'UMC2RF: ',-I,
     $           ' invalid entries ignored (not in prior pattern).'
              ENDIF 
           ELSE IF (ERROR .EQ. 2) THEN 
              WRITE (IO, 1) 'UMC2RF: ', I,' duplicate entries summed.'
           ELSE IF (ERROR .EQ. 4) THEN 
              WRITE (IO, 1) 'UMC2RF: matrix is singular.  Only ', I,
     $        ' pivots found.'

C          -------------------------------------------------------------
C          UMC2RF non-cumulative warning messages (internal error codes)
C          -------------------------------------------------------------

           ELSE IF (ERROR .EQ. 99) THEN 
              WRITE (IO, 2)
     $        'UMC2RF: invalid entry (out of range 1..N):', I, J, X
           ELSE IF (ERROR .EQ. 98) THEN 
              WRITE (IO, 2)
     $        'UMC2RF: duplicate entry summed:', I, J, X
           ELSE IF (ERROR .EQ. 97) THEN 
              WRITE (IO, 2)
     $        'UMC2RF: invalid entry (not in pattern of prior factors)',
     $        I, J, X
           ELSE IF (ERROR .EQ. 96) THEN 
              WRITE (IO, 2)
     $        'UMC2RF: invalid entry (below diagonal blocks):', I, J, X
           ELSE IF (ERROR .EQ. 95) THEN 
              WRITE (IO, 2)
     $        'UMC2RF: invalid entry (prior matrix singular):', I, J, X

           ELSE IF (ERROR .EQ. -8) THEN 
              WRITE (IO, 1) 'UMC2RF: transa argument must be .false.!'

           ENDIF 

        ELSE IF (WHO .EQ. 3) THEN 

C          -------------------------------------------------------------
C          UMC2SO error messages
C          -------------------------------------------------------------

           IF (ERROR .EQ. -7) THEN 
              WRITE (IO, 1)
     $        'UMC2SO: LU factors uncomputed or corrupted!'

C          -------------------------------------------------------------
C          UMC2SO non-cumulative warning messages
C          -------------------------------------------------------------

           ELSE IF (ERROR .EQ. 8) THEN 
              IF (I .EQ. 0) THEN 
                 WRITE (IO, 1)
     $  'UMC2SO: no iterative refinement: original matrix not preserved'
              ELSE 
                 WRITE (IO, 1)
     $  'UMC2SO: no iterative refinement: only for Ax=b or A''x=b'
              ENDIF 

           ELSE IF (ERROR .EQ. -8) THEN 
              WRITE (IO, 1) 'UMC2SO: transc argument must be .false.!'

           ENDIF 

        ENDIF 

        RETURN

C=======================================================================
C  FORMAT STATMENTS
C=======================================================================

1       FORMAT (' ', A, :, I12, :, A)
2       FORMAT (' ', A,/, '    row: ', I12, ' col: ', I12, ' ',
     $          '(real: ',E11.4,' imag: ',E11.4,')')
        END 
