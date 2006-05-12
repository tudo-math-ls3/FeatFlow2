c                                                                       

c**********************************************************************
c                   SPLIB Copyright (C) 1995                         **
c                     Indiana University                             **
c**********************************************************************

c-----------------------------------------------------------------------
c     logical unit number for printout of memory manager
c-----------------------------------------------------------------------
c                                                                       

      integer             iprint
      common / disks  /   iprint
c
c-----------------------------------------------------------------------
c arrays needed by pseudo dynamic memory allocation routine 'addtbl'
c-----------------------------------------------------------------------
c maxtbl........ maximum number of arrays that can be allocated
c                by the pseudo dynamic memory allocator
c given......... given allocated for work array.
c
      integer           given
      integer           maxtbl
      parameter       ( maxtbl=5000 )
      integer ibeg  , iend  , tblsiz, ntable, nxttbl
      character*6 table , origin
      character*7 tabtyp
      common / TABLE1 / table (maxtbl), tabtyp(maxtbl), origin(maxtbl)
      common / TABLE2 / ibeg  (maxtbl), iend  (maxtbl), tblsiz(maxtbl),
     &                  ntable        , nxttbl, given
