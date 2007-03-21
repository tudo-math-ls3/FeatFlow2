c                                                                       

c**********************************************************************
c                   SPLIB Copyright (C) 1995                         **
c                     Indiana University                             **
c**********************************************************************

c-----------------------------------------------------------------------
c     number of integer words per real variable on the target computer.
c-----------------------------------------------------------------------
c                             
      real*8 tstdbl                                          
      parameter  ( tstdbl = 0 )
      integer*4    tstint
      parameter  ( tstint = 0 )
      integer      intdbl
      parameter  ( intdbl = KIND(tstdbl)/KIND(tstint) )
c
