c                                                                       

c**********************************************************************
c                   SPLIB Copyright (C) 1995                         **
c                     Indiana University                             **
c**********************************************************************

c-----------------------------------------------------------------------
c     number of integer words per real variable on the target computer.
c-----------------------------------------------------------------------
c                             
      double precision tstdbl                                          
      parameter  ( tstdbl = 0 )
      integer      intdbl
      parameter  ( intdbl = KIND(tstdbl)/KIND(0) )
c
