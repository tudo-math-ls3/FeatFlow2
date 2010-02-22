/* $Header: /src4/opennurbs/opennurbs_memory_custom.c 2     12/08/03 2:41p Dalelear $ */
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
#error Do not directly compile opennurbs_memory_custom.c.  See opennurbs_memory.c.
#endif

#pragma message( " --- OpenNURBS using custom memory manager" )

/*
// This file is included by opennurbs_memory.c.
//
// If you want to provide custom memory management, search for
// "TODO" in this file and fill in the blanks.  Then define
// ON_MEMORY_CUSTOM and recompile opennurbs_memory.c.
//
// For examples, see opennurbs_memory_default.c and
// opennurbs_memory_smartheap.c.
*/

static void* onmemorypool_malloc( ON_MEMORY_POOL* pool, size_t sz )
{
  // Input:
  //
  //   pool = value returned by onmemorypool_getcurrentpool() or 
  //          onmemorypool_getmainpool().
  //
  //   sz > 0
  //
  // Output:
  //   NULL or sz bytes of memory

  // TODO - allocate sz bytes of memory
  // return malloc( sz );
  return 0;
}

static void* onmemorypool_calloc( ON_MEMORY_POOL* pool, size_t num, size_t sz )
{
  // Input:
  //
  //   pool = value returned by onmemorypool_getcurrentpool() or 
  //          onmemorypool_getmainpool().
  //
  //   num > 0
  //   sz > 0
  //
  // Output:
  //   NULL or num*sz bytes of zeroed memory

  // TODO - allocate and zero num*sz bytes of memory
  // return calloc( num, sz );
  return 0;
}

static void onmemorypool_free( void* memblock )
{
  // Input:
  //
  //   memblock = non NULL pointer returned by onmemorypool_malloc(),
  //              onmemorypool_calloc(), or onmemorypool_realloc().

  // TODO - free memblock for reuse
  // free(memblock);
}


static void* onmemorypool_realloc( ON_MEMORY_POOL* pool, void* memblock, size_t sz )
{
  // Input:
  //
  //   pool = value returned by onmemorypool_getcurrentpool() or 
  //          onmemorypool_getmainpool().
  //
  //   memblock = non NULL pointer returned by onmemorypool_malloc(),
  //              onmemorypool_calloc(), or onmemorypool_realloc().
  //
  //   sz > 0
  //
  // Output:
  //   NULL or sz bytes of memory. The first min(sz,oldsz) bytes of
  //   returned block must have the same values as the input memblock.
  //   You may change the location of the block.

  // TODO - resize block 
  // return realloc( memblock, sz );
  return 0;
}

static size_t onmemorypool_msize( const void* memblock )
{
  // Input:
  //
  //   memblock = non NULL pointer returned by onmemorypool_malloc(),
  //              onmemorypool_calloc(), or onmemorypool_realloc().
  //
  // Output:
  //   size in bytes of the block

  // TODO - return number of usiable bytes in this block
  // return _msize((void*)p);
  return 0;
}

static void onmemorypool_begin(void)
{
  // TODO - whatever you want to do when OpenNURBS is initialized
  return;
}

static void onmemorypool_end(void)
{
  // TODO - whatever you want to do when OpenNURBS is terminated
  return;
}

static ON_MEMORY_POOL* onmemorypool_getmainpool( void )
{
  // TODO - return pointer to main pool
  return 0;
}

static ON_MEMORY_POOL* onmemorypool_createpool( void )
{
  // TODO - return pointer to a new pool
  return 0;
}

static void onmemorypool_destroypool( ON_MEMORY_POOL* pool )
{
  // TODO - destroy pool
  return;
}

static void onmemorypool_compactpool( ON_MEMORY_POOL* pool )
{
  // TODO - compact pool
  return;
}

