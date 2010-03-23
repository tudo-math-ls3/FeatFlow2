/* $Header: /src4/opennurbs/opennurbs_memory.c 4     9/22/05 2:52p Dalelear $ */
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
#include "opennurbs_error.h"

#if defined(ON_DLL_IMPORTS)
/*
// If you use OpenNURBS as a windows DLL, then define ON_DLL_IMPORTS 
// in applications that use OpenNURBS and they will get onmalloc(), etc.,
// from the DLL.
//
// If you use OpenNURBS as a static library, do not define ON_DLL_IMPORTS.
*/
#error opennurbs_memory.c must not be compiled with ON_DLL_IMPORTS defined.
#endif

// An opennurbs_memory_*.c file must provide definitions of the
// static onmemorypool_*() functions declared below.  If you want
// to provide custom memory management, then you should fill in the
// "TODO" areas in opennurbs_memory_custom.c and define ON_MEMORY_CUSTOM.

static void   onmemorypool_begin( void );
static void   onmemorypool_end( void );

static void*  onmemorypool_malloc( ON_MEMORY_POOL*, size_t );
static void*  onmemorypool_calloc( ON_MEMORY_POOL*, size_t, size_t );
static void   onmemorypool_free( void* );
static void*  onmemorypool_realloc( ON_MEMORY_POOL*, void*, size_t );
static size_t onmemorypool_msize( const void* );

//static ON_MEMORY_POOL* onmemorypool_getmainpool( void );
static ON_MEMORY_POOL* onmemorypool_createpool( void );
static void            onmemorypool_destroypool( ON_MEMORY_POOL* );
static void            onmemorypool_compactpool( ON_MEMORY_POOL* );


static int ON_memory_error(int);

/*
// Memory error handling
//
// See the declaration of ON_memory_error_register_handler() for instructions.
*/
static ON_memory_error_handler ON_memory_error_handler_func = 0;

static int ON_memory_error(int error_number)
{
  /*
  // error_number 0: malloc()/calloc()/realloc() returned 0
  //              1: invalid pointer passed into realloc()
  //              2: invalid pointer passed into free()
  //              3: msize() returned 0 or 0xFFFFFFFF
  //              4: ON_CreateMemoryPool() failed.
  //              5: ON_DestroyMemoryPool() failed.
  //              6: onmemorypool_free() failed.
  //              7: ON_MemoryCompactPool() failed.
  //              8: onmemorypool_exit_thread() failed.
  //              9: bad pool passed to onmemorypool_malloc().
  //
  // If this function returns a nonzero value, the call is repeated.
  // If this function returns zero, the call fails as gracefully as possible.
  */
  return (ON_memory_error_handler_func) ? ON_memory_error_handler_func(error_number) : 0;
}

#define OPENNURBS_MEMORY_CFILE_

/* Ucomment to use custom memory manager defined in opennurbs_memory_custom.c
#define ON_MEMORY_CUSTOM 
*/

/* Uncomment to use SmartHeap memory manager defined in opennurbs_memory_smartheap.c
// You must buy a SmartHeap license and link with the SmartHeap library.
#define ON_MEMORY_SMARTHEAP
*/

#if defined(ON_MEMORY_CUSTOM)
#include "opennurbs_memory_custom.c"
#elif defined(ON_MEMORY_SMARTHEAP)
#include "opennurbs_memory_smartheap.c"
#elif defined(ON_MEMORY_MICROSOFT_HEAP)
#include "opennurbs_memory_microsoft.c"
#else
#include "opennurbs_memory_default.c"
#endif



ON_memory_error_handler ON_memory_error_register_handler( ON_memory_error_handler f )
{
  ON_memory_error_handler oldf = ON_memory_error_handler_func;
  ON_memory_error_handler_func = f;
  return oldf;
}


#if defined(ON_DEBUG)
#define ON_COUNT_MEMORYUSE
#endif

static size_t onmemoryusecounters[4] = {0,0,0,0}; // {malloc_count,realloc_count,free_count,pool_count}

size_t onmemoryusecount(size_t* malloc_count, size_t* realloc_count, size_t* free_count, size_t* pool_count )
{
  // The counters are uses only when ON_COUNT_MEMORYUSE is defined.
  if ( malloc_count )
    *malloc_count = onmemoryusecounters[0];
  if ( realloc_count )
    *realloc_count = onmemoryusecounters[1];
  if ( free_count )
    *free_count = onmemoryusecounters[2];
  if ( pool_count )
    *pool_count = onmemoryusecounters[3];
  return (onmemoryusecounters[0]+onmemoryusecounters[1]+onmemoryusecounters[2]+onmemoryusecounters[3]);
}

/*
// Memory pool handling
*/

static ON_MEMORY_POOL* glb_pMainPool = 0;
static ON_MEMORY_POOL* glb_pCurrentPool = 0;


ON_MEMORY_POOL* ON_CreateMemoryPool( void )
{
  ON_MEMORY_POOL* pool = onmemorypool_createpool();
  if ( 0 == pool )
  {
    ON_memory_error(4);
  }
#if defined(ON_COUNT_MEMORYUSE)
  onmemoryusecounters[3]++;
#endif
  return pool;
}


void ON_DestroyMemoryPool( ON_MEMORY_POOL* pool )
{
  onmemorypool_destroypool(pool);
}


void ON_CompactMemoryPool( ON_MEMORY_POOL* pool )
{
  onmemorypool_compactpool(pool);
}


ON_MEMORY_POOL* ON_SetCurrentMemoryPool( ON_MEMORY_POOL* pool )
{
  ON_MEMORY_POOL* prevpool = glb_pCurrentPool;
  glb_pCurrentPool = pool;
  return prevpool;
}


ON_MEMORY_POOL* ON_GetCurrentMemoryPool(void)
{
  return glb_pCurrentPool;
}


ON_MEMORY_POOL* ON_SetMainMemoryPool( ON_MEMORY_POOL* pool )
{
  ON_MEMORY_POOL* prevpool = glb_pMainPool;
  glb_pMainPool = pool;
  return prevpool;
}


ON_MEMORY_POOL* ON_GetMainMemoryPool(void)
{
  return glb_pMainPool;
}


void* onmalloc_from_pool( ON_MEMORY_POOL* pool, size_t sz )
{
  void* p;
#if defined(ON_COUNT_MEMORYUSE)
  onmemoryusecounters[0]++;
#endif
  if (sz <= 0 )
  {
    // 4 Sep 2003 - using onmalloc(0) to check for cancel
    onmemorypool_malloc( pool, 0 );
    return 0;
  }
  for(;;) {
    p = onmemorypool_malloc( pool, sz );
    if (p)
      break;
    if (!ON_memory_error(0))
      break;
  }
  return p;
}


void* onmalloc( size_t sz )
{
  return onmalloc_from_pool( glb_pCurrentPool, sz );
}


void* oncalloc_from_pool( ON_MEMORY_POOL* pool, size_t num, size_t sz )
{
  void* p;
#if defined(ON_COUNT_MEMORYUSE)
  onmemoryusecounters[0]++;
#endif
  if ( num<= 0 || sz <=0 )
    return 0;
  for(;;) {
    p = onmemorypool_calloc( pool, num, sz );
    if (p)
      break;
    if (!ON_memory_error(0))
      break;
  }
  return p;
}


void* oncalloc( size_t num, size_t sz )
{
  return oncalloc_from_pool( glb_pCurrentPool, num, sz );
}


void onfree( void* memblock )
{
#if defined(ON_COUNT_MEMORYUSE)
  onmemoryusecounters[2]++;
#endif
  if ( memblock )
    onmemorypool_free( memblock );
}


void* onrealloc( void* memblock, size_t sz )
{
  return onrealloc_from_pool( glb_pCurrentPool, memblock, sz );
}


void* onrealloc_from_pool( ON_MEMORY_POOL* pool, void* memblock, size_t sz )
{
  void* p;
  if ( sz <= 0 ) {
    onfree(memblock);
    return 0;
  }
  if ( !memblock ) {
    return onmalloc_from_pool( pool, sz);
  }
#if defined(ON_COUNT_MEMORYUSE)
  onmemoryusecounters[1]++;
#endif
  for(;;) {
    p = onmemorypool_realloc( pool, memblock, sz );
    if (p)
      break;
    if (!ON_memory_error(0))
      break;
  }
  return p;
}


size_t onmsize( const void* memblock )
{
  size_t sz = 0;
  if (memblock) {
    for(;;) {
      sz = onmemorypool_msize( memblock );
      if ( sz > 0 && sz != 0xFFFFFFFF )
        break;
      if ( !ON_memory_error(3) ) {
        sz = 0;
        break;
      }
    }
  }
  return sz;
}

static int glb_OpenNURBS_Is_Running = 0;


void ON_OpenNURBS_Begin(void)
{
  if ( !glb_OpenNURBS_Is_Running ) {
    glb_OpenNURBS_Is_Running = 1;
    onmemorypool_begin();
  }
}


void ON_OpenNURBS_End(void)
{
  if (glb_OpenNURBS_Is_Running) {
    glb_OpenNURBS_Is_Running = 0;
    onmemorypool_end();
  }
}

