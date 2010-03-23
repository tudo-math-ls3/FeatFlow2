/* $Header: /src4/opennurbs/opennurbs_memory_default.c 2     12/02/04 8:51a Dalelear $ */
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

#if !defined(OPENNURBS_MEMORY_CFILE_)
#error Do not directly compile opennurbs_memory_default.c.  See opennurbs_memory.c.
#endif

#if !defined(_GNU_SOURCE)
// The reason this is hidden from gcc is to suppress the gcc warning 
// that it is ignoring "message".  (These warnings generate support
// questions.)
#pragma message( " --- OpenNURBS using default OS memory manager" )
#endif

/*
// This file is included by opennurbs_memory.c.  Do NOT edit this file.
//
// If you want to provide custom memory management, put your 
// memory management code in opennurbs_memory_custom.c and change
// the file included in opennurbs_memory.c
//
// default memory managment provided by 
// system malloc()/calloc()/realloc()/free()/_msize()
*/

static void* onmemorypool_malloc( ON_MEMORY_POOL* pool, size_t sz )
{
  // 4 Sep 2003 - possible for sz=0 to make it to this point
  return ((sz>0)?malloc( sz ):0);
}

static void* onmemorypool_calloc( ON_MEMORY_POOL* pool, size_t num, size_t sz )
{
  return calloc(num,sz);
}

static void onmemorypool_free( void* memblock )
{
  free(memblock);
}

#if defined(_MSC_VER)
#if _MSC_VER == 1200
/*
//  (_MSC_VER is defined as 1200 for Microsoft Visual C++ 6.0)
//
//   NOTE WELL: Microsoft's VC 6.0 realloc() contains a bug that can cause
//              crashes and should be avoided.  See MSDN Knowledge Base
//              article ID Q225099 for more information.
*/
#define ON_REALLOC_BROKEN
#endif
#endif

static void* onmemorypool_realloc( ON_MEMORY_POOL* pool, void* memblock, size_t sz )
{
#if defined(ON_REALLOC_BROKEN)
  /* use malloc() and memcpy() instead of buggy realloc() */
  void* p;
  const size_t memblocksz = _msize(memblock);
  if ( sz <= memblocksz ) {
    /* shrink */
    if ( memblocksz <= 28 || 8*sz >= 7*memblocksz ) 
    {
      /* don't bother reallocating */
      p = memblock;
    }
    else {
      /* allocate smaller block */
      p = malloc(sz);
      if ( p ) 
      {
        memcpy( p, memblock, sz );
        free(memblock);
      }
    }
  }
  else if ( sz > memblocksz ) {
    /* grow */
    p = malloc(sz);
    if ( p ) {
      memcpy( p, memblock, memblocksz );
      free(memblock);
    }
  }
  return p;
#else
  return realloc( memblock, sz );
#endif
}

static size_t onmemorypool_msize( const void* p )
{
#if defined(ON_OS_WINDOWS)
  return _msize((void*)p);
#else
  // OS doesn't support _msize().
  return 0;
#endif
}

#if defined(ON_DEBUG)
  // This does nothing.  No "memory pool" is used with this
  // memory manager.
  // This is here for debugging purposes so the "default"
  // memory manager behaves like the SmartHeap version.
static ON_MEMORY_POOL glb_mainpool = {0,0,0,0};
static int glb_pool_count = 0;
#endif

static void onmemorypool_begin(void)
{
#if defined(ON_DEBUG)
  // This does nothing.  No "memory pool" is used with this
  // memory manager.
  // This is here for debugging purposes so the "default"
  // memory manager behaves like the SmartHeap version.
  // See opennurbs_memory_smartheap.c for details.
  static BOOL b = 0;
  if ( !b ) {
    b = 1;
    glb_mainpool.m_i = ++glb_pool_count;
    ON_SetMainMemoryPool(&glb_mainpool);
    ON_SetCurrentMemoryPool(&glb_mainpool);
  }
#endif

  return;
}

static void onmemorypool_end(void)
{
  return;
}

/*
static ON_MEMORY_POOL* onmemorypool_getmainpool( void )
{
  return 0;
}
*/

static ON_MEMORY_POOL* onmemorypool_createpool( void )
{
#if defined(ON_DEBUG)
  // This does nothing.  No "memory pool" is created.
  // This is here for debugging purposes so the "default"
  // memory manager behaves like the SmartHeap version.
  ON_MEMORY_POOL* pool = (ON_MEMORY_POOL*)oncalloc( 1, sizeof(*pool) );

  pool->m__s[ 0] = '0';
  pool->m__s[ 1] = 'N';
  pool->m__s[ 2] = '_';
  pool->m__s[ 3] = 'M';
  pool->m__s[ 4] = 'E';
  pool->m__s[ 5] = 'M';
  pool->m__s[ 6] = 'O';
  pool->m__s[ 7] = 'R';
  pool->m__s[ 8] = 'Y';
  pool->m__s[ 9] = '_';
  pool->m__s[10] = 'P';
  pool->m__s[11] = 'O';
  pool->m__s[12] = 'O';
  pool->m__s[13] = 'L';
  pool->m__s[14] = 0;

  pool->m_i = ++glb_pool_count;
  return pool;
#else
  return 0;
#endif
}

static void onmemorypool_destroypool( ON_MEMORY_POOL* pool )
{
#if defined(ON_DEBUG)
  // This does nothing.  No "memory pool" is destroyed.
  // This is here for debugging purposes so the "default"
  // memory manager behaves like the SmartHeap version.
  if ( pool && pool->m_i > 0 ) {
    pool->m_i = -pool->m_i;
    pool->m_p = 0;
  }
#endif
}

static void onmemorypool_compactpool( ON_MEMORY_POOL* pool )
{
  return;
}

