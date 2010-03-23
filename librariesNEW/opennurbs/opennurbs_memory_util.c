/* $Header: /src4/opennurbs/opennurbs_memory_util.c 2     6/08/05 7:40a Dalelear $ */
/* $NoKeywords: $ */
/*
//
// Copyright (c) 1993-2001 Robert McNeel & Associates. All rights reserved.
// Rhinoceros is a registered trademark of Robert McNeel & Assoicates.
//
// THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY.
// ALL IMPLIED WARRANTIES OF FITNESS FOR ANY PARTICULAR PURPOSE AND OF
// MERCHANTABILITY ARE HEREBY DISCLAIMED.
//				
// For complete openNURBS copyright information see <http://www.opennurbs.org>.
//
////////////////////////////////////////////////////////////////
*/

#include "opennurbs_system.h"
#include "opennurbs_defines.h"
#include "opennurbs_memory.h"

// Memory utilities used by OpenNURBS.  If you are using
// custom memory managment, NOTHING in this file needs to
// be changed.


void* onmemdup( const void* src, size_t sz )
{
  void* p;
  if ( src && sz>0 ) 
  {
    p = onmalloc(sz);
    if (p)
      memcpy(p,src,sz);
  }
  else {
    p = 0;
  }
  return p;
}


char* onstrdup( const char* src )
{
  char* p;
  size_t sz;
  if ( src ) 
  {
    for ( sz=0;*src++;sz++)
      ; /* empty for body */
    sz++;
    p = (char*)onmemdup( src-sz, sz*sizeof(*src) );
  }
  else 
  {
    p = 0;
  }
  return p;
}


unsigned char* onmbsdup( const unsigned char* src )
{
  unsigned char* p;
  size_t sz; /* sz = number of bytes to dup (>=_mbclen(scr)) */
  if ( src ) 
  {
    for ( sz=0;*src++;sz++)
      ; /* empty for body */
    sz++;
    p = (unsigned char*)onmemdup( src-sz, sz*sizeof(*src) );
  }
  else 
  {
    p = 0;
  }
  return p;
}

#if defined(_WCHAR_T_DEFINED)

wchar_t* onwcsdup( const wchar_t* src )
{
  wchar_t* p;
  size_t sz;
  if ( src ) 
  {
    for ( sz=0;*src++;sz++)
      ; /* empty for body */
    sz++;
    p = (wchar_t*)onmemdup( src-sz, sz*sizeof(*src) );
  }
  else 
  {
    p = 0;
  }
  return p;
}
#endif
