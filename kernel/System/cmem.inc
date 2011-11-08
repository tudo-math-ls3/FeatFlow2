!**********************************************************************
!* This file contains the public common block and the definitions
!* responsible for pseudodynamic memory management.
!*
!* The parameter NNWORK is to be modified by the user. It describes
!* the total amount of memory (in doubles) that is used by the
!* main problem.
!*
!* NNARR describes the total amount of handles available to the
!* program.
!**********************************************************************

      integer NNWORK, NNARR

! NNWORK = number of double-values in the DWORK-block
!
! 64 Bit:
! --------
! 500 MB:
!      parameter (NNWORK=35000000)
! 1 GB:
!      parameter (NNWORK=70000000)
! 2 GB:
!      parameter (NNWORK=140000000)
! 3 GB:
!      parameter (NNWORK=210000000)
! 4 GB:
!      parameter (NNWORK=280000000)
! fast 6 GB:
!      parameter (NNWORK=400000000)

! 32 Bit:
! --------
! 12.8 MB
!      parameter (NNWORK=2000000)
! 64 MB
!      parameter (NNWORK=10000000)
! 128 MB
!      parameter (NNWORK=20000000)
! 256 MB
      parameter (NNWORK=40000000)
! 500 MB:
!      parameter (NNWORK=70000000)
! 1 GB:
!      parameter (NNWORK=140000000)
! 1.5 GB:
!      parameter (NNWORK=210000000)
! 2 GB:
!      parameter (NNWORK=280000000)
! nearly 3 GB:
!      parameter (NNWORK=400000000)

! NNARR = max. number of handles that can be used:

      parameter (NNARR=299)

! Now the declaration of the main global structures follow

      integer NWORK, IWORK, IWMAX, L

      double precision, dimension(NNWORK) :: DWORK
      real, dimension(NNWORK*kind(DWORK(1))/kind(0.0))  :: VWORK
      integer, dimension(NNWORK*kind(DWORK(1))/kind(0)) :: KWORK
      
      common  NWORK,IWORK,IWMAX,L(NNARR),DWORK

! equivalence statement; essential for pseudodyn. memory management

      equivalence (DWORK(1),VWORK(1),KWORK(1))
