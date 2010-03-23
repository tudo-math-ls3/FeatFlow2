/* $Header: /src4/opennurbs/opennurbs_memory_microsoft.c 7     8/10/06 4:13p Dalelear $ */
/* $NoKeywords: $ */
/*
//
// Copyright (c) 1993-2003 Robert McNeel & Associates. All rights reserved.
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
#error Do not directly compile opennurbs_memory_microsoft.c.  See opennurbs_memory.c.
#endif

#pragma message( " --- OpenNURBS using Microsoft Heap memory manager" )

#if defined(_DEBUG)
#define DEBUG_WORKER_THREAD_MEMORY

#if defined(WATCH_ALLOCATION)

// Copied from Microsoft's VIsual Studio dbginit.h file
// and used to debug heap corruption.

// Put the range of request numbers here that are causing trouble.
#define ON___nRequest0 20436
#define ON___nRequest1 20436

#define ON___nNoMansLandSize 4

typedef struct ON___CrtMemBlockHeader
{
        struct ON___CrtMemBlockHeader * pBlockHeaderNext;
        struct ON___CrtMemBlockHeader * pBlockHeaderPrev;
        char *                      szFileName;
        int                         nLine;
#ifdef _WIN64
        /* These items are reversed on Win64 to eliminate gaps in the struct
         * and ensure that sizeof(struct)%16 == 0, so 16-byte alignment is
         * maintained in the debug heap.
         */
        int                         nBlockUse;
        size_t                      nDataSize;
#else  /* _WIN64 */
        size_t                      nDataSize;
        int                         nBlockUse;
#endif  /* _WIN64 */
        long                        lRequest;
        unsigned char               gap[ON___nNoMansLandSize];
        /* followed by:
         *  unsigned char           data[nDataSize];
         *  unsigned char           anotherGap[ON___nNoMansLandSize];
         */
};

typedef struct ON___CrtMemBlockEndGap
{
  unsigned char           anotherGap[ON___nNoMansLandSize];
};

static
ON___CheckCrtMemBlockHeader(void* ptr)
{
  int breakpoint_here = 0;
  struct ON___CrtMemBlockEndGap* pGap = 0;
  struct ON___CrtMemBlockHeader* pHead = (struct ON___CrtMemBlockHeader*)ptr;
  pHead--;

  // use these condtions to limit checking to requests you care about.
  if (     pHead->lRequest >= 20436
       &&  pHead->lRequest <= 20436
       )
  {
    char* byte_ptr = (char*)ptr;
    pGap = (struct ON___CrtMemBlockEndGap*)(byte_ptr + pHead->nDataSize);
    breakpoint_here = 1; // put a breakpoint here
  }
}

#endif


#endif

#if defined(DEBUG_WORKER_THREAD_MEMORY)
static void onmemorypool_debug_check(const ON_MEMORY_POOL* pool,int memory_error_value)
{
  if ( pool->m_i > 1 )
  {
    if ( 0 == pool->m_p || 0 != pool->m_j )
    {
      // pool->m_i>1 means pool is an active worker thread pool.  It
      // must have a pointer to a heap (m_p!=0) and m_j should be 
      // zero.  When the worker thread is running, m_h is is the
      // thread's handle.  A few allocations happen before the worker
      // thread starts because some GSLib caches need to be
      // built using memory from the worker thread pool so they
      // can be correctly cleaned up when a cancels a task.
      ON_memory_error(memory_error_value);
    }
  }
  else if ( 1 == pool->m_i )
  {
    if ( 0 != pool->m_p || 0 != pool->m_h || 0 != pool->m_j )
    {
      // 1==pool->m_i means pool is the main application global
      // pool.  
      // If cannot be deleted (m_j!=0) and should never have a thread
      // (m_h) set.
      // Since this particular memory manager uses the default MSC
      // CRT malloc/realloc/free for the main application global
      // pool, m_p should be NULL.  For SmartHeap, m_p will not be NULL.
      ON_memory_error(memory_error_value);
    }
  }
  else
  {
    // pool->m_i < 1 means pool has been deleted and shouldn't
    // be used
    ON_memory_error(memory_error_value);
  }
}
#endif


static void onmemorypool_exit_thread(ON_MEMORY_POOL* pool)
{
  /*
  // This function is called if pool is a worker thread pool
  // and pool->m_j != 0, which indicates the parent application
  // wants to exit the thread.  Some important checks are performed
  // to make sure we are exiting the correct thread, and then we
  // call ExitThread().
  */
  if ( pool->m_i > 1 && 0 != pool->m_h )
  {
    /* 
    // m_i > 1 means this is a worker pool and ExitThread has
    // has never been called.  We make sure that the
    // current thread is really the one we want to exit and
    // then we call ExitThread().
    //
    // In some situations, it is possible that the main thread is
    // active and the h == pool->m_h test will fail.  This means
    // we just need to wait a bit and try again.
    */
    unsigned long h = GetCurrentThreadId();
    if ( h == pool->m_h )
    {
      /*
      // Do NOT set pool->m_p = 0 or negate m_i here.  Otherwise
      // you will leak the memory allocated by the doomed worker 
      // thread.  After the thread exists, the main application
      // will call ON_DestroyMemoryPool(), which will call
      // onmemorypool_destroypool(), which will destroy the
      // heap, set m_p=0, and negate m_i.
      */
      pool->m_h = 0;
      ExitThread(pool->m_k);
    }
  }
  else
  {
    /*
    // If we get this far, then something is really wrong.
    */
    ON_memory_error(8);
  }
}

static void* onmemorypool_malloc( ON_MEMORY_POOL* pool, size_t sz )
{
  void* ptr = 0;
  // This function is only called from inside of opennurbs_memory.c and
  // the caller insures that 0!=pool.

#if defined(DEBUG_WORKER_THREAD_MEMORY)
  onmemorypool_debug_check(pool,9);
#endif

  if ( pool->m_j ) 
  {
    // In general, this will exit the thread and the rest
    // of this function will never get executed.  If it doesn't
    // exit the thread, then we have to do the best we can.
    // It may be that the main thread is active or we are allocating
    // memory from the main pool, and we have to wait until we are
    // back in the worker thread and using the worker pool.
    onmemorypool_exit_thread(pool);
  }


  // 4 Sep 2003 Dale Lear
  //   It is now possible for sz=0 to make it to this point when
  //   onmalloc(0) is used to make worker threads
  //   more responsive to canceling.
  if (sz>0)
  {
    if ( pool->m_i > 1 )
    {
      ON_MEMORY_POOL** pp = (ON_MEMORY_POOL**)HeapAlloc(pool->m_p,0,sz+8);
      if ( pp )
      {
        *pp = pool;
        ptr = ((char*)pp)+8;
      }
    }
    else
    {
      // Use CRT memory management
      ON_MEMORY_POOL** pp = (ON_MEMORY_POOL**)malloc(sz+8);
      if ( pp )
      {
#if defined(WATCH_ALLOCATION)
        ON___CheckCrtMemBlockHeader(pp);
#endif
        *pp = 0;
        ptr = ((char*)pp)+8;
      }
    }
  }

  return ptr;
}

static void* onmemorypool_calloc( ON_MEMORY_POOL* pool, size_t num, size_t sz )
{
  // This function is only called from inside of opennurbs_memory.c and
  // the caller insures that 0!=pool.  We use a onmemorypool_malloc()
  // and then memset because the 8 byte header we use to keep track
  // of the memory's parent heap.
  void* ptr = onmemorypool_malloc(pool,num*sz);
  if (ptr)
    memset(ptr,0,num*sz);
  return ptr;
}

static void onmemorypool_free( void* memblock )
{
  // This function is only called from inside of opennurbs_memory.c and
  // the caller insures that 0!=memblock;
  ON_MEMORY_POOL** pp = (ON_MEMORY_POOL**)(((char*)memblock)-8);
  if ( *pp )
  {
#if defined(DEBUG_WORKER_THREAD_MEMORY)
  onmemorypool_debug_check(*pp,10);
#endif

    // do not call onmemorypool_exit_thread() here
    if ( !HeapFree((*pp)->m_p, 0, pp ) )
      ON_memory_error(6);
  }
  else
  {
    // Use CRT memory management
#if defined(WATCH_ALLOCATION)
    ON___CheckCrtMemBlockHeader(pp);
#endif
    free(pp);
  }
}

static void* onmemorypool_realloc( ON_MEMORY_POOL* pool, void* memblock, size_t sz )
{
  // This function is only called from inside of opennurbs_memory.c and
  // the caller insures that 0!=pool,0!=memblock,sz>0;
  //
  // The input value of "pool" is ignored because a reallocation must come
  // from the same pool as memblock.

  ON_MEMORY_POOL** pp = (ON_MEMORY_POOL**)(((char*)memblock)-8);

#if defined(DEBUG_WORKER_THREAD_MEMORY)
  onmemorypool_debug_check(pool,11);
#endif

  if ( *pp )
  {
#if defined(DEBUG_WORKER_THREAD_MEMORY)
    // *pp can be different from pool
    onmemorypool_debug_check(*pp,12);
#endif
    // Memory belongs to a heap created by onmemorypool_createpool().
    pp = HeapReAlloc( (*pp)->m_p, 0, pp, sz+8);
  }
  else
  {
    // Use CRT memory management
#if defined(WATCH_ALLOCATION)
    ON___CheckCrtMemBlockHeader(pp);
#endif
    pp = realloc( pp, sz+8 );
  }
  return (pp) ? (((char*)pp)+8) : 0;
}

static size_t onmemorypool_msize( const void* memblock )
{
  // This function is only called from inside of opennurbs_memory.c and
  // the caller insures that 0!=memblock.
  size_t sz = 0;
  ON_MEMORY_POOL** pp = (ON_MEMORY_POOL**)(((char*)memblock)-8);

  if ( *pp )
  {
#if defined(DEBUG_WORKER_THREAD_MEMORY)
    onmemorypool_debug_check(*pp,14);
#endif
    sz = HeapSize( (*pp)->m_p, 0, pp );
  }
  else
  {
    sz = _msize(pp);
  }

  if ( 0xFFFFFFFF == sz || sz < 8 )
    sz = 0;
  else
    sz -= 8;

  return sz;
}

static ON_MEMORY_POOL glb_mainpool = {0,0,0,0};
static int glb_pool_count = 0;

static void onmemorypool_begin(void)
{
  // This function must be called before any memory allocation is attempted.
  // If you have globals or static member variables who initialization 
  // calls onmalloc()/oncalloc()/onrealloc()/onfree(), then this function
  // needs to be called before those variable are initialized.  One way
  // to accomplish this is to use OpenNURBS as a DLL and call this function
  // from DllMain();
  static BOOL b = 0;
  if ( !b ) 
  {

#if !defined(_DEBUG)
    // For debug builds, we will use the CRT malloc/realloc/free
    // so we can look for leaks in the VC OUTPUT window.
    ULONG heap_info = 2;
    b = 1;
    glb_mainpool.m_p = GetProcessHeap();
    HeapSetInformation(glb_mainpool.m_p,
                       HeapCompatibilityInformation,
                       &heap_info,
                       sizeof(heap_info)
                       );
#endif

    glb_mainpool.m_i = ++glb_pool_count;

    ON_SetMainMemoryPool(&glb_mainpool);
    ON_SetCurrentMemoryPool(&glb_mainpool);

#if defined(ON_OS_WINDOWS) && defined(_DEBUG)
    OutputDebugStringA("Initializing OpenNURBS's Microsoft Heap manager\n");
#endif

  }
  return;
}

static void onmemorypool_end(void)
{
  // Do not destroy existing pools as some memory may be 
  // deallocated after main()/WinMain() exits.
  return;
}

static ON_MEMORY_POOL* onmemorypool_createpool( void )
{
  static DWORD dwPageSize = 0;
  ON_MEMORY_POOL* pool = 0;
  DWORD flOptions = 0;
  SIZE_T dwInitialSize = 0;
  SIZE_T dwMaximumSize = 0; // heap can be grown
  HANDLE heap = 0;

  if ( 0 == dwPageSize )
  {
    SYSTEM_INFO system_info;
    memset(&system_info,0,sizeof(system_info));
    GetSystemInfo(&system_info);
    dwPageSize = system_info.dwPageSize;
  }
  dwInitialSize = 1024*dwPageSize;

  heap = HeapCreate( flOptions, dwInitialSize, dwMaximumSize );

  if ( 0 != heap )
  {
    // use a low fragmentation heap
    ULONG heap_info = 2;
    HeapSetInformation( heap,
                        HeapCompatibilityInformation,
                        &heap_info,
                        sizeof(heap_info)
                        );

    pool = (ON_MEMORY_POOL*)calloc( 1, sizeof(*pool) );

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

    pool->m_p = heap;
    pool->m_i = ++glb_pool_count;
  }

  return pool;
}

static void onmemorypool_destroypool( ON_MEMORY_POOL* pool )
{
  if ( pool && 0 != pool->m_p && pool->m_i > 1 ) 
  {
    pool->m_i = -pool->m_i;
    if (!HeapDestroy( pool->m_p ))
      ON_memory_error(5);
    pool->m_p = 0;
  }
  else
    ON_memory_error(5);
}

static void onmemorypool_compactpool( ON_MEMORY_POOL* pool )
{
  if ( 0 != pool && 0 != pool->m_p && pool->m_i > 1 )
  {
    if ( HeapLock( pool->m_p ) )
    {
      PROCESS_HEAP_ENTRY heap_entry;
      memset(&heap_entry,0,sizeof(heap_entry));
      while ( HeapWalk(pool->m_p,&heap_entry) )
      {
        if ( heap_entry.wFlags & PROCESS_HEAP_ENTRY_BUSY )
        {
          // heap has a memory that is in use - compact heap
          if( HeapUnlock( pool->m_p ) )
          {
            if (!HeapCompact(pool->m_p,0))
               ON_memory_error(7);
          }
          else
            ON_memory_error(7);
          return;
        }
      }
      if (HeapUnlock( pool->m_p ))
      {
        // no memory is in use - destroy the pool
        onmemorypool_destroypool( pool );
      }
      else
        ON_memory_error(7);
    }
  }
  else
    ON_memory_error(7);
}


