!**********************************************************************
!* This file contains the public COMMON block and the definitions
!* responsible for pseudodynamic memory management.
!*
!* The parameter NNWORK is to be modified by the user. It describes
!* the total amount of memory (in doubles) that is used by the
!* main problem.
!*
!* NNARR describes the total amount of handles available to the
!* program.
!**********************************************************************

      INTEGER NNWORK, NNARR

! NNWORK = number of DOUBLE-values in the DWORK-block
!
! 64 Bit:
! --------
! 500 MB:
!      PARAMETER (NNWORK=35000000)
! 1 GB:
!      PARAMETER (NNWORK=70000000)
! 2 GB:
!      PARAMETER (NNWORK=140000000)
! 3 GB:
!      PARAMETER (NNWORK=210000000)
! 4 GB:
!      PARAMETER (NNWORK=280000000)
! fast 6 GB:
!      PARAMETER (NNWORK=400000000)

! 32 Bit:
! --------
! 12.8 MB
!      PARAMETER (NNWORK=2000000)
! 64 MB
!      PARAMETER (NNWORK=10000000)
! 128 MB
!      PARAMETER (NNWORK=20000000)
! 256 MB
      PARAMETER (NNWORK=40000000)
! 500 MB:
!      PARAMETER (NNWORK=70000000)
! 1 GB:
!      PARAMETER (NNWORK=140000000)
! 1.5 GB:
!      PARAMETER (NNWORK=210000000)
! 2 GB:
!      PARAMETER (NNWORK=280000000)
! nearly 3 GB:
!      PARAMETER (NNWORK=400000000)

! NNARR = max. number of handles that can be used:

      PARAMETER (NNARR=299)

! Now the declaration of the main global structures follow

      INTEGER NWORK, IWORK, IWMAX, L
      
      DOUBLE PRECISION, DIMENSION(NNWORK) :: DWORK
      REAL, DIMENSION(NNWORK*KIND(DWORK(1))/KIND(0.0))  :: VWORK
      INTEGER, DIMENSION(NNWORK*KIND(DWORK(1))/KIND(0)) :: KWORK
      
      COMMON  NWORK,IWORK,IWMAX,L(NNARR),DWORK

! Equivalence statement; essential for pseudodyn. memory management

      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
