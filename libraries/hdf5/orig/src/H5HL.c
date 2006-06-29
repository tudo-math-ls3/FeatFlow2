/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdf.ncsa.uiuc.edu/HDF5/doc/Copyright.html.  If you do not have     *
 * access to either file, you may request a copy from hdfhelp@ncsa.uiuc.edu. *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*-------------------------------------------------------------------------
 *
 * Created:		H5HL.c
 *			Jul 16 1997
 *			Robb Matzke <matzke@llnl.gov>
 *
 * Purpose:		Heap functions for the local heaps used by symbol
 *			tables to store names (among other things).
 *
 * Modifications:
 *
 *	Robb Matzke, 5 Aug 1997
 *	Added calls to H5E.
 *
 *-------------------------------------------------------------------------
 */
#define H5F_PACKAGE		/* Suppress error about including H5Fpkg  */
#define H5HL_PACKAGE		/* Suppress error about including H5HLpkg */


#include "H5private.h"		/* Generic Functions			*/
#include "H5ACprivate.h"	/* Metadata cache			*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5Fpkg.h"             /* File access				*/
#include "H5FLprivate.h"	/* Free lists                           */
#include "H5HLpkg.h"		/* Local Heaps				*/
#include "H5MFprivate.h"	/* File memory management		*/

/* Private macros */
#define H5HL_FREE_NULL	1		/*end of free list on disk	*/
#define H5HL_MIN_HEAP   256             /* Minimum size to reduce heap buffer to */

/*
 * Local heap collection version.
 */
#define H5HL_VERSION	0

/* Private typedefs */

/* PRIVATE PROTOTYPES */
#ifdef NOT_YET
static void *H5HL_read(H5F_t *f, hid_t dxpl_id, haddr_t addr, size_t offset, size_t size,
			void *buf);
static herr_t H5HL_write(H5F_t *f, hid_t dxpl_id, haddr_t addr, size_t offset, size_t size,
			  const void *buf);
#endif /* NOT_YET */

static herr_t H5HL_serialize(H5F_t *f, H5HL_t *heap, uint8_t *buf);
static H5HL_free_t *H5HL_remove_free(H5HL_t *heap, H5HL_free_t *fl);
static herr_t H5HL_minimize_heap_space(H5F_t *f, hid_t dxpl_id, H5HL_t *heap);

/* Metadata cache callbacks */
static H5HL_t *H5HL_load(H5F_t *f, hid_t dxpl_id, haddr_t addr, const void *udata1,
			 void *udata2);
static herr_t H5HL_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr, H5HL_t *heap);
static herr_t H5HL_dest(H5F_t *f, H5HL_t *heap);
static herr_t H5HL_clear(H5F_t *f, H5HL_t *heap, hbool_t destroy);
static herr_t H5HL_compute_size(const H5F_t *f, const H5HL_t *heap, size_t *size_ptr);

/*
 * H5HL inherits cache-like properties from H5AC
 */
const H5AC_class_t H5AC_LHEAP[1] = {{
    H5AC_LHEAP_ID,
    (H5AC_load_func_t)H5HL_load,
    (H5AC_flush_func_t)H5HL_flush,
    (H5AC_dest_func_t)H5HL_dest,
    (H5AC_clear_func_t)H5HL_clear,
    (H5AC_size_func_t)H5HL_compute_size,
}};

/* Declare a free list to manage the H5HL_free_t struct */
H5FL_DEFINE_STATIC(H5HL_free_t);

/* Declare a free list to manage the H5HL_t struct */
H5FL_DEFINE_STATIC(H5HL_t);

/* Declare a PQ free list to manage the heap chunk information */
H5FL_BLK_DEFINE_STATIC(heap_chunk);


/*-------------------------------------------------------------------------
 * Function:	H5HL_create
 *
 * Purpose:	Creates a new heap data structure on disk and caches it
 *		in memory.  SIZE_HINT is a hint for the initial size of the
 *		data area of the heap.	If size hint is invalid then a
 *		reasonable (but probably not optimal) size will be chosen.
 *		If the heap ever has to grow, then REALLOC_HINT is the
 *		minimum amount by which the heap will grow.
 *
 * Return:	Success:	Non-negative. The file address of new heap is
 *				returned through the ADDR argument.
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *		matzke@llnl.gov
 *		Jul 16 1997
 *
 * Modifications:
 *
 *	Robb Matzke, 5 Aug 1997
 *	Takes a flag that determines the type of heap that is
 *	created.
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5HL_create(H5F_t *f, hid_t dxpl_id, size_t size_hint, haddr_t *addr_p/*out*/)
{
    H5HL_t	*heap = NULL;
    hsize_t	total_size;		/*total heap size on disk	*/
    size_t      sizeof_hdr;             /* Cache H5HL header size for file */
    herr_t	ret_value = SUCCEED;

    FUNC_ENTER_NOAPI(H5HL_create, FAIL);

    /* check arguments */
    assert(f);
    assert(addr_p);

    if (size_hint && size_hint < H5HL_SIZEOF_FREE(f))
	size_hint = H5HL_SIZEOF_FREE(f);
    size_hint = H5HL_ALIGN(size_hint);

    /* Cache this for later */
    sizeof_hdr= H5HL_SIZEOF_HDR(f);

    /* allocate file version */
    total_size = sizeof_hdr + size_hint;
    if (HADDR_UNDEF==(*addr_p=H5MF_alloc(f, H5FD_MEM_LHEAP, dxpl_id, total_size)))
	HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "unable to allocate file memory");

    /* allocate memory version */
    if (NULL==(heap = H5FL_CALLOC(H5HL_t)))
	HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");
    heap->addr = *addr_p + (hsize_t)sizeof_hdr;
    heap->disk_alloc = size_hint;
    heap->mem_alloc = size_hint;
    heap->disk_resrv = 0;
    if (NULL==(heap->chunk = H5FL_BLK_CALLOC(heap_chunk,(sizeof_hdr + size_hint))))
	HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");

    /* free list */
    if (size_hint) {
	if (NULL==(heap->freelist = H5FL_MALLOC(H5HL_free_t)))
	    HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");
	heap->freelist->offset = 0;
	heap->freelist->size = size_hint;
	heap->freelist->prev = heap->freelist->next = NULL;
    } else {
	heap->freelist = NULL;
    }

    /* add to cache */
    heap->cache_info.is_dirty = TRUE;
    if (H5AC_set(f, dxpl_id, H5AC_LHEAP, *addr_p, heap) < 0)
	HGOTO_ERROR(H5E_HEAP, H5E_CANTINIT, FAIL, "unable to cache heap");

done:
    if (ret_value<0) {
	if (H5F_addr_defined(*addr_p))
	    H5MF_xfree(f, H5FD_MEM_LHEAP, dxpl_id, *addr_p, total_size);
	if (heap) {
            if(H5HL_dest(f,heap)<0)
                HDONE_ERROR(H5E_HEAP, H5E_CANTFREE, FAIL, "unable to destroy local heap collection");
	}
    }
    FUNC_LEAVE_NOAPI(ret_value);
}


/*-------------------------------------------------------------------------
 * Function:	H5HL_load
 *
 * Purpose:	Loads a heap from disk.
 *
 * Return:	Success:	Ptr to a local heap memory data structure.
 *
 *		Failure:	NULL
 *
 * Programmer:	Robb Matzke
 *		matzke@llnl.gov
 *		Jul 17 1997
 *
 * Modifications:
 *		Robb Matzke, 1999-07-28
 *		The ADDR argument is passed by value.
 *
 *	Quincey Koziol, 2002-7-180
 *	Added dxpl parameter to allow more control over I/O from metadata
 *      cache.
 *-------------------------------------------------------------------------
 */
static H5HL_t *
H5HL_load(H5F_t *f, hid_t dxpl_id, haddr_t addr, const void UNUSED * udata1,
	  void UNUSED * udata2)
{
    uint8_t		hdr[52];
    size_t              sizeof_hdr;     /* Cache H5HL header size for file */
    const uint8_t	*p = NULL;
    H5HL_t		*heap = NULL;
    H5HL_free_t		*fl = NULL, *tail = NULL;
    size_t		free_block = H5HL_FREE_NULL;
    H5HL_t		*ret_value;

    FUNC_ENTER_NOAPI(H5HL_load, NULL);

    /* check arguments */
    assert(f);
    assert(H5F_addr_defined(addr));
    assert(!udata1);
    assert(!udata2);

    /* Cache this for later */
    sizeof_hdr= H5HL_SIZEOF_HDR(f);
    assert(sizeof_hdr <= sizeof(hdr));

    /* Get the local heap's header */
    if (H5F_block_read(f, H5FD_MEM_LHEAP, addr, sizeof_hdr, dxpl_id, hdr) < 0)
	HGOTO_ERROR(H5E_HEAP, H5E_READERROR, NULL, "unable to read heap header");
    p = hdr;

    /* Check magic number */
    if (HDmemcmp(hdr, H5HL_MAGIC, H5HL_SIZEOF_MAGIC))
	HGOTO_ERROR(H5E_HEAP, H5E_CANTLOAD, NULL, "bad heap signature");
    p += H5HL_SIZEOF_MAGIC;

    /* Version */
    if (H5HL_VERSION!=*p++)
	HGOTO_ERROR (H5E_HEAP, H5E_CANTLOAD, NULL, "wrong version number in global heap");

    /* Reserved */
    p += 3;

    /* Allocate space in memory for the heap */
    if (NULL==(heap = H5FL_CALLOC(H5HL_t)))
	HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");

    /* heap data size */
    H5F_DECODE_LENGTH(f, p, heap->disk_alloc);
    heap->mem_alloc = heap->disk_alloc;

    /* free list head */
    H5F_DECODE_LENGTH(f, p, free_block);
    if (free_block != H5HL_FREE_NULL && free_block >= heap->disk_alloc)
	HGOTO_ERROR(H5E_HEAP, H5E_CANTLOAD, NULL, "bad heap free list");

    /* data */
    H5F_addr_decode(f, &p, &(heap->addr));
    if (NULL==(heap->chunk = H5FL_BLK_CALLOC(heap_chunk,(sizeof_hdr + heap->mem_alloc))))
	HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");
    if (heap->disk_alloc &&
            H5F_block_read(f, H5FD_MEM_LHEAP, heap->addr, heap->disk_alloc, dxpl_id, heap->chunk + sizeof_hdr) < 0)
	HGOTO_ERROR(H5E_HEAP, H5E_CANTLOAD, NULL, "unable to read heap data");

    /* Build free list */
    while (H5HL_FREE_NULL != free_block) {
	if (free_block >= heap->disk_alloc)
	    HGOTO_ERROR(H5E_HEAP, H5E_CANTLOAD, NULL, "bad heap free list");
	if (NULL==(fl = H5FL_MALLOC(H5HL_free_t)))
	    HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");
	fl->offset = free_block;
	fl->prev = tail;
	fl->next = NULL;
	if (tail) tail->next = fl;
	tail = fl;
	if (!heap->freelist) heap->freelist = fl;

	p = heap->chunk + sizeof_hdr + free_block;
	H5F_DECODE_LENGTH(f, p, free_block);
	H5F_DECODE_LENGTH(f, p, fl->size);

	if (fl->offset + fl->size > heap->disk_alloc)
	    HGOTO_ERROR(H5E_HEAP, H5E_CANTLOAD, NULL, "bad heap free list");
    }

    /* Set return value */
    ret_value = heap;

done:
    if (!ret_value && heap) {
        if(H5HL_dest(f,heap)<0)
	    HDONE_ERROR(H5E_HEAP, H5E_CANTFREE, NULL, "unable to destroy local heap collection");
    }

    FUNC_LEAVE_NOAPI(ret_value);
}


/*-------------------------------------------------------------------------
 * Function:    H5HL_minimize_heap_space
 *
 * Purpose:     Go through the heap's freelist and determine if we can
 *              eliminate the free blocks at the tail of the buffer.
 *
 * Return:      Success:        SUCCEED
 *              Failure:        FAIL
 *
 * Programmer:  Bill Wendling
 *              wendling@ncsa.uiuc.edu
 *              Sept. 16, 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5HL_minimize_heap_space(H5F_t *f, hid_t dxpl_id, H5HL_t *heap)
{
    herr_t      ret_value = SUCCEED;
    size_t      sizeof_hdr;

    FUNC_ENTER_NOAPI(H5HL_minimize_heap_space, FAIL)

    /* check args */
    assert(f);
    assert(heap);

    sizeof_hdr = H5HL_SIZEOF_HDR(f);    /* cache H5HL header size for file */

    /*
     * When the heap is being flushed to disk, release the file space reserved
     * for it.
     */
    H5MF_free_reserved(f, (hsize_t)heap->disk_resrv);
    heap->disk_resrv = 0;

    /*
     * Check to see if we can reduce the size of the heap in memory by
     * eliminating free blocks at the tail of the buffer before flushing the
     * buffer out.
     */
    if (heap->freelist) {
        H5HL_free_t    *tmp_fl;
        H5HL_free_t    *last_fl = NULL;

        /* Search for a free block at the end of the buffer */
        for (tmp_fl = heap->freelist; tmp_fl; tmp_fl = tmp_fl->next)
            /* Check if the end of this free block is at the end of the buffer */
            if (tmp_fl->offset + tmp_fl->size == heap->mem_alloc) {
                last_fl = tmp_fl;
                break;
            }

        /*
         * Found free block at the end of the buffer, decide what to do
         * about it
         */
        if (last_fl) {
            size_t  new_mem_size = heap->mem_alloc; /* New size of memory buffer */

            /*
             * If the last free block's size is more than half the memory
             * buffer size (and the memory buffer is larger than the
             * minimum size), reduce or eliminate it.
             */
            if (last_fl->size >= (heap->mem_alloc / 2) && heap->mem_alloc > H5HL_MIN_HEAP) {
                /*
                 * Reduce size of buffer until it's too small or would
                 * eliminate the free block
                 */
                while (new_mem_size > H5HL_MIN_HEAP &&
                        new_mem_size >= (last_fl->offset + H5HL_SIZEOF_FREE(f)))
                    new_mem_size /= 2;

                /*
                 * Check if reducing the memory buffer size would
                 * eliminate the free list
                 */
                if (new_mem_size < (last_fl->offset + H5HL_SIZEOF_FREE(f))) {
                    /* Check if this is the only block on the free list */
                    if (last_fl->prev == NULL && last_fl->next == NULL) {
                        /* Double the new memory size */
                        new_mem_size *= 2;

                        /* Truncate the free block */
                        last_fl->size = H5HL_ALIGN(new_mem_size - last_fl->offset);
                        new_mem_size = last_fl->offset + last_fl->size;
                        assert(last_fl->size >= H5HL_SIZEOF_FREE(f));
                    } else {
                        /*
                         * Set the size of the memory buffer to the start
                         * of the free list
                         */
                        new_mem_size = last_fl->offset;

                        /* Eliminate the free block from the list */
                        last_fl = H5HL_remove_free(heap, last_fl);
                    }
                } else {
                    /* Truncate the free block */
                    last_fl->size = H5HL_ALIGN(new_mem_size - last_fl->offset);
                    new_mem_size = last_fl->offset + last_fl->size;
                    assert(last_fl->size >= H5HL_SIZEOF_FREE(f));
                    assert(last_fl->size == H5HL_ALIGN(last_fl->size));
                }

                /* Resize the memory buffer and reserved space in file */
                if (new_mem_size != heap->mem_alloc) {

                    heap->mem_alloc = new_mem_size;
                    heap->chunk = H5FL_BLK_REALLOC(heap_chunk, heap->chunk, (sizeof_hdr + new_mem_size));

                    if (!heap->chunk)
                        HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed")
                }
            }
        }
    }

    /*
     * If the heap grew larger or smaller than disk storage then move the
     * data segment of the heap to another contiguous block of disk
     * storage.
     */
    if (heap->mem_alloc != heap->disk_alloc) {
        haddr_t old_addr = heap->addr, new_addr;

        /* Release old space on disk */
        H5_CHECK_OVERFLOW(heap->disk_alloc, size_t, hsize_t);
        H5MF_xfree(f, H5FD_MEM_LHEAP, dxpl_id, old_addr, (hsize_t)heap->disk_alloc);
        H5E_clear();    /* don't really care if the free failed */

        /* Allocate new space on disk */
        H5_CHECK_OVERFLOW(heap->mem_alloc, size_t, hsize_t);

        if (HADDR_UNDEF == (new_addr = H5MF_alloc(f, H5FD_MEM_LHEAP, dxpl_id, (hsize_t)heap->mem_alloc)))
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "unable to allocate file space for heap")

        heap->addr = new_addr;

        /* Set new size of block on disk */
        heap->disk_alloc = heap->mem_alloc;
    }

done:
    FUNC_LEAVE_NOAPI(ret_value)
}


/*-------------------------------------------------------------------------
 * Function:    H5HL_serialize
 *
 * Purpose:     Serialize the heap. This function will eliminate free
 *              blocks at the tail of the heap and also split the block
 *              if it needs to be split for the file. This is so that we
 *              can serialize it correctly.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 * Programmer:  Bill Wendling
 *              wendling@ncsa.uiuc.edu
 *              Sept. 16, 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5HL_serialize(H5F_t *f, H5HL_t *heap, uint8_t *buf)
{
    H5HL_free_t	   *fl;
    uint8_t        *p;
    herr_t          ret_value = SUCCEED;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5HL_serialize)

    /* check args */
    assert(buf);
    assert(heap);

    /* serialize the header */
    p = buf;
    fl = heap->freelist;
    HDmemcpy(p, H5HL_MAGIC, H5HL_SIZEOF_MAGIC);
    p += H5HL_SIZEOF_MAGIC;
    *p++ = H5HL_VERSION;
    *p++ = 0;	/*reserved*/
    *p++ = 0;	/*reserved*/
    *p++ = 0;	/*reserved*/
    H5F_ENCODE_LENGTH(f, p, heap->mem_alloc);
    H5F_ENCODE_LENGTH(f, p, fl ? fl->offset : H5HL_FREE_NULL);
    H5F_addr_encode(f, &p, heap->addr);

    /* serialize the free list */
    for (; fl; fl = fl->next) {
        assert (fl->offset == H5HL_ALIGN (fl->offset));
        p = heap->chunk + H5HL_SIZEOF_HDR(f) + fl->offset;

        if (fl->next) {
            H5F_ENCODE_LENGTH(f, p, fl->next->offset);
        } else {
            H5F_ENCODE_LENGTH(f, p, H5HL_FREE_NULL);
        }

        H5F_ENCODE_LENGTH(f, p, fl->size);
    }

    FUNC_LEAVE_NOAPI(ret_value)
}


/*-------------------------------------------------------------------------
 * Function:	H5HL_flush
 *
 * Purpose:	Flushes a heap from memory to disk if it's dirty.  Optionally
 *		deletes the heap from memory.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		matzke@llnl.gov
 *		Jul 17 1997
 *
 * Modifications:
 *      rky, 1998-08-28
 *      Only p0 writes metadata to disk.
 *
 *      Robb Matzke, 1999-07-28
 *      The ADDR argument is passed by value.
 *
 *	Quincey Koziol, 2002-7-180
 *	Added dxpl parameter to allow more control over I/O from metadata
 *      cache.
 *
 *      Bill Wendling, 2003-09-16
 *      Separated out the bit that serializes the heap.
 *-------------------------------------------------------------------------
 */
static herr_t
H5HL_flush(H5F_t *f, hid_t dxpl_id, hbool_t destroy, haddr_t addr, H5HL_t *heap)
{
    herr_t  ret_value = SUCCEED;    /* Return value */

    FUNC_ENTER_NOAPI(H5HL_flush, FAIL);

    /* check arguments */
    assert(f);
    assert(H5F_addr_defined(addr));
    assert(heap);

    if (heap->cache_info.is_dirty) {
        haddr_t hdr_end_addr;
        size_t  sizeof_hdr = H5HL_SIZEOF_HDR(f);    /* cache H5HL header size for file */

        /* Minimize the heap space size if possible */
        if (H5HL_minimize_heap_space(f, dxpl_id, heap) < 0)
            HGOTO_ERROR(H5E_HEAP, H5E_CANTFREE, FAIL, "unable to minimize local heap space")

	/* Write the header */
        if (H5HL_serialize(f, heap, heap->chunk) < 0)
            HGOTO_ERROR(H5E_BTREE, H5E_CANTSERIALIZE, FAIL, "unable to serialize local heap")

	/* Copy buffer to disk */
	hdr_end_addr = addr + (hsize_t)sizeof_hdr;

	if (H5F_addr_eq(heap->addr, hdr_end_addr)) {
	    /* The header and data are contiguous */
	    if (H5F_block_write(f, H5FD_MEM_LHEAP, addr, (sizeof_hdr + heap->disk_alloc),
				dxpl_id, heap->chunk) < 0)
		HGOTO_ERROR(H5E_HEAP, H5E_WRITEERROR, FAIL, "unable to write heap header and data to file")
	} else {
	    if (H5F_block_write(f, H5FD_MEM_LHEAP, addr, sizeof_hdr, dxpl_id, heap->chunk) < 0)
		HGOTO_ERROR(H5E_HEAP, H5E_WRITEERROR, FAIL, "unable to write heap header to file")

	    if (H5F_block_write(f, H5FD_MEM_LHEAP, heap->addr, heap->disk_alloc,
				dxpl_id, heap->chunk + sizeof_hdr) < 0)
		HGOTO_ERROR(H5E_HEAP, H5E_WRITEERROR, FAIL, "unable to write heap data to file")
	}

	heap->cache_info.is_dirty = FALSE;
    }

    /* Should we destroy the memory version? */
    if (destroy) {
        if (H5HL_dest(f,heap) < 0)
	    HGOTO_ERROR(H5E_HEAP, H5E_CANTFREE, FAIL, "unable to destroy local heap collection")
    }

done:
    FUNC_LEAVE_NOAPI(ret_value)
}


/*-------------------------------------------------------------------------
 * Function:	H5HL_dest
 *
 * Purpose:	Destroys a heap in memory.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *		koziol@ncsa.uiuc.edu
 *		Jan 15 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5HL_dest(H5F_t UNUSED *f, H5HL_t *heap)
{
    H5HL_free_t	*fl;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5HL_dest);

    /* check arguments */
    assert(heap);

    /* Verify that node is clean */
    assert (heap->cache_info.is_dirty==FALSE);

    if(heap->chunk)
        heap->chunk = H5FL_BLK_FREE(heap_chunk,heap->chunk);
    while (heap->freelist) {
        fl = heap->freelist;
        heap->freelist = fl->next;
        H5FL_FREE(H5HL_free_t,fl);
    }
    H5FL_FREE(H5HL_t,heap);

    FUNC_LEAVE_NOAPI(SUCCEED);
} /* end H5HL_dest() */


/*-------------------------------------------------------------------------
 * Function:	H5HL_clear
 *
 * Purpose:	Mark a local heap in memory as non-dirty.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *		koziol@ncsa.uiuc.edu
 *		Mar 20 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5HL_clear(H5F_t *f, H5HL_t *heap, hbool_t destroy)
{
    herr_t ret_value = SUCCEED;

    FUNC_ENTER_NOAPI_NOINIT(H5HL_clear);

    /* check arguments */
    assert(heap);

    /* Mark heap as clean */
    heap->cache_info.is_dirty = FALSE;

    if (destroy)
        if (H5HL_dest(f, heap) < 0)
	    HGOTO_ERROR(H5E_HEAP, H5E_CANTFREE, FAIL, "unable to destroy local heap collection");

done:
    FUNC_LEAVE_NOAPI(ret_value);
} /* end H5HL_clear() */


/*-------------------------------------------------------------------------
 * Function:	H5HL_compute_size
 *
 * Purpose:	Compute the size in bytes of the specified instance of
 *              H5HL_t on disk, and return it in *len_ptr.  On failure,
 *              the value of *len_ptr is undefined.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	John Mainzer
 *		5/13/04
 *
 * Modifications:
 *-------------------------------------------------------------------------
 */
static herr_t
H5HL_compute_size(const H5F_t *f, const H5HL_t *heap, size_t *size_ptr)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5HL_compute_size);

    /* check arguments */
    HDassert(f);
    HDassert(heap);
    HDassert(size_ptr);

    *size_ptr = H5HL_SIZEOF_HDR(f) + heap->disk_alloc;

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* H5HL_compute_size() */


/*-------------------------------------------------------------------------
 * Function:	H5HL_read
 *
 * Purpose:	Reads some object (or part of an object) from the heap
 *		whose address is ADDR in file F.  OFFSET is the byte offset
 *		from the beginning of the heap at which to begin reading
 *		and SIZE is the number of bytes to read.
 *
 *		If BUF is the null pointer then a buffer is allocated by
 *		this function.
 *
 *		Attempting to read past the end of an object may cause this
 *		function to fail.
 *
 * Return:	Success:	BUF (or the allocated buffer)
 *
 *		Failure:	NULL
 *
 * Programmer:	Robb Matzke
 *		matzke@llnl.gov
 *		Jul 16 1997
 *
 * Modifications:
 *		Robb Matzke, 1999-07-28
 *		The ADDR argument is passed by value.
 *-------------------------------------------------------------------------
 */
#ifdef NOT_YET
static void *
H5HL_read(H5F_t *f, hid_t dxpl_id, haddr_t addr, size_t offset, size_t size, void *buf)
{
    H5HL_t	*heap = NULL;
    void      *ret_value;       /* Return value */

    FUNC_ENTER_NOAPI(H5HL_read, NULL);

    /* check arguments */
    assert(f);
    assert (H5F_addr_defined(addr));

    if (NULL == (heap = H5AC_protect(f, dxpl_id, H5AC_LHEAP, addr, NULL, NULL, H5AC_READ)))
	HGOTO_ERROR(H5E_HEAP, H5E_CANTLOAD, NULL, "unable to load heap");

    assert(offset < heap->mem_alloc);
    assert(offset + size <= heap->mem_alloc);

    if (!buf && NULL==(buf = H5MM_malloc(size)))
	HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");
    HDmemcpy(buf, heap->chunk + H5HL_SIZEOF_HDR(f) + offset, size);

    /* Set return value */
    ret_value=buf;

done:
    if (heap && H5AC_unprotect(f, dxpl_id, H5AC_LHEAP, addr, heap, FALSE) != SUCCEED)
        HDONE_ERROR(H5E_HEAP, H5E_PROTECT, NULL, "unable to release object header");

    FUNC_LEAVE_NOAPI(ret_value);
}
#endif /* NOT_YET */


/*-------------------------------------------------------------------------
 * Function:    H5HL_protect
 *
 * Purpose:     This function is a wrapper for the H5AC_protect call. The
 *              old H5HL_peek call (which this once was) wasn't "safe"
 *              for FPHDF5. (It'd get a read lock on an object but once
 *              it got that object, it'd release it keeping the old
 *              pointer value, which is no longer valid. This won't work
 *              since the pointer into some metdata block can become
 *              invalid.)
 *
 *              N.B.: This function is always called in conjunction with
 *              the H5HL_offset_into function. The return from that
 *              function is the proper pointer to the heap's object. This
 *              is done so that the return from this function can be sent
 *              to H5HL_unprotect.
 *
 * Return:      Success:    Ptr to the object. The pointer points to a
 *                          chunk of memory large enough to hold the
 *                          object from the specified offset (usually the
 *                          beginning of the object) to the end of the
 *                          object. Do not attempt to read past the end
 *                          of the object.
 *              Failure:    NULL
 *
 * Programmer:  Bill Wendling
 *              wendling@ncsa.uiuc.edu
 *              Sept. 17, 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
const H5HL_t *
H5HL_protect(H5F_t *f, hid_t dxpl_id, haddr_t addr)
{
    H5HL_t *ret_value = NULL;

    FUNC_ENTER_NOAPI(H5HL_protect, NULL);

    /* check arguments */
    assert(f);
    assert(H5F_addr_defined(addr));

    if (NULL == (ret_value = H5AC_protect(f, dxpl_id, H5AC_LHEAP, addr, NULL, NULL, H5AC_READ)))
        HGOTO_ERROR(H5E_HEAP, H5E_CANTLOAD, NULL, "unable to load heap");

done:
    FUNC_LEAVE_NOAPI(ret_value);
}


/*-------------------------------------------------------------------------
 * Function:    H5HL_offset_into
 *
 * Purpose:     Called directly after the call to H5HL_protect so that
 *              a pointer to the object in the heap can be got.
 *
 * Return:      Success:    Valid pointer.
 *              Failure:    NULL
 *
 * Programmer:  Bill Wendling
 *              wendling@ncsa.uiuc.edu
 *              Sept. 17, 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void *
H5HL_offset_into(H5F_t *f, const H5HL_t *heap, size_t offset)
{
    /*
     * We need to have called some other function before this to get a
     * valid heap pointer. So, this can remain "FUNC_ENTER_NOAPI_NOINIT"
     */
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5HL_offset_into)
    assert(f);
    assert(heap);
    assert(offset < heap->mem_alloc);
    FUNC_LEAVE_NOAPI(heap->chunk + H5HL_SIZEOF_HDR(f) + offset)
}


/*-------------------------------------------------------------------------
 * Function:    H5HL_unprotect
 *
 * Purpose:     Unprotect the data retrieved by the H5HL_protect call.
 *
 * Return:      Success:    SUCCEED
 *              Failure:    FAIL
 *
 * Programmer:  Bill Wendling
 *              wendling@ncsa.uiuc.edu
 *              Sept. 17, 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5HL_unprotect(H5F_t *f, hid_t dxpl_id, const H5HL_t *heap, haddr_t addr)
{
    herr_t  ret_value = SUCCEED;

    FUNC_ENTER_NOAPI(H5HL_unprotect, FAIL);

    /* check arguments */
    assert(f);
    assert(heap);
    assert(H5F_addr_defined(addr));

    if (H5AC_unprotect(f, dxpl_id, H5AC_LHEAP, addr, (void *)heap, FALSE) != SUCCEED)
        HGOTO_ERROR(H5E_HEAP, H5E_PROTECT, FAIL, "unable to release object header");

done:
    FUNC_LEAVE_NOAPI(ret_value);
}


/*-------------------------------------------------------------------------
 * Function:	H5HL_remove_free
 *
 * Purpose:	Removes free list element FL from the specified heap and
 *		frees it.
 *
 * Return:	NULL
 *
 * Programmer:	Robb Matzke
 *		matzke@llnl.gov
 *		Jul 17 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static H5HL_free_t *
H5HL_remove_free(H5HL_t *heap, H5HL_free_t *fl)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5HL_remove_free);

    if (fl->prev) fl->prev->next = fl->next;
    if (fl->next) fl->next->prev = fl->prev;

    if (!fl->prev) heap->freelist = fl->next;

    FUNC_LEAVE_NOAPI(H5FL_FREE(H5HL_free_t,fl));
}


/*-------------------------------------------------------------------------
 * Function:	H5HL_insert
 *
 * Purpose:	Inserts a new item into the heap.
 *
 * Return:	Success:	Offset of new item within heap.
 *
 *		Failure:	(size_t)(-1)
 *
 * Programmer:	Robb Matzke
 *		matzke@llnl.gov
 *		Jul 17 1997
 *
 * Modifications:
 *		Robb Matzke, 1999-07-28
 *		The ADDR argument is passed by value.
 *-------------------------------------------------------------------------
 */
size_t
H5HL_insert(H5F_t *f, hid_t dxpl_id, haddr_t addr, size_t buf_size, const void *buf)
{
    H5HL_t	*heap = NULL;
    H5HL_free_t	*fl = NULL, *max_fl = NULL;
    size_t	offset = 0;
    size_t	need_size, old_size, need_more;
    hbool_t	found;
    size_t	disk_resrv;	/* Amount of additional space to reserve in file */
    size_t      sizeof_hdr;     /* Cache H5HL header size for file */
    size_t	ret_value;      /* Return value */

    FUNC_ENTER_NOAPI(H5HL_insert, (size_t)(-1));

    /* check arguments */
    assert(f);
    assert(H5F_addr_defined(addr));
    assert(buf_size > 0);
    assert(buf);

    if (0==(f->intent & H5F_ACC_RDWR))
	HGOTO_ERROR (H5E_HEAP, H5E_WRITEERROR, (size_t)(-1), "no write intent on file");

    if (NULL == (heap = H5AC_protect(f, dxpl_id, H5AC_LHEAP, addr, NULL, NULL, H5AC_WRITE)))
	HGOTO_ERROR(H5E_HEAP, H5E_PROTECT, (size_t)(-1), "unable to load heap");

    heap->cache_info.is_dirty=TRUE;

    /* Cache this for later */
    sizeof_hdr= H5HL_SIZEOF_HDR(f);

    /*
     * In order to keep the free list descriptors aligned on word boundaries,
     * whatever that might mean, we round the size up to the next multiple of
     * a word.
     */
    need_size = H5HL_ALIGN(buf_size);

    /*
     * Look for a free slot large enough for this object and which would
     * leave zero or at least H5G_SIZEOF_FREE bytes left over.
     */
    for (fl=heap->freelist, found=FALSE; fl; fl=fl->next) {
	if (fl->size > need_size &&
	    fl->size - need_size >= H5HL_SIZEOF_FREE(f)) {
	    /* a bigger free block was found */
	    offset = fl->offset;
	    fl->offset += need_size;
	    fl->size -= need_size;
	    assert (fl->offset==H5HL_ALIGN (fl->offset));
	    assert (fl->size==H5HL_ALIGN (fl->size));
	    found = TRUE;
	    break;
	} else if (fl->size == need_size) {
	    /* free block of exact size found */
	    offset = fl->offset;
	    fl = H5HL_remove_free(heap, fl);
	    found = TRUE;
	    break;
	} else if (!max_fl || max_fl->offset < fl->offset) {
	    /* use worst fit */
	    max_fl = fl;
	}
    }

    /*
     * If no free chunk was large enough, then allocate more space and
     * add it to the free list.	 If the heap ends with a free chunk, we
     * can extend that free chunk.  Otherwise we'll have to make another
     * free chunk.  If the heap must expand, we double its size.
     */
    if (found==FALSE) {
	need_more = MAX3(need_size, heap->mem_alloc, H5HL_SIZEOF_FREE(f));

	/* Reserve space in the file to hold the increased heap size
	 */
	if( heap->disk_resrv == heap->mem_alloc)
		disk_resrv = need_more;
	else
		disk_resrv = heap->mem_alloc + need_more - heap->disk_resrv;

	if( H5MF_reserve(f, (hsize_t)disk_resrv) < 0 )
		HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, (size_t)(-1), "unable to reserve space in file");

	/* Update heap's record of how much space it has reserved */
	heap->disk_resrv += disk_resrv;

	if (max_fl && max_fl->offset + max_fl->size == heap->mem_alloc) {
	    /*
	     * Increase the size of the maximum free block.
	     */
	    offset = max_fl->offset;
	    max_fl->offset += need_size;
	    max_fl->size += need_more - need_size;
	    assert (max_fl->offset==H5HL_ALIGN (max_fl->offset));
	    assert (max_fl->size==H5HL_ALIGN (max_fl->size));

	    if (max_fl->size < H5HL_SIZEOF_FREE(f)) {
#ifdef H5HL_DEBUG
		if (H5DEBUG(HL) && max_fl->size) {
		    fprintf(H5DEBUG(HL), "H5HL: lost %lu bytes at line %d\n",
			    (unsigned long)(max_fl->size), __LINE__);
		}
#endif
		max_fl = H5HL_remove_free(heap, max_fl);
	    }
	} else {
	    /*
	     * Create a new free list element large enough that we can
	     * take some space out of it right away.
	     */
	    offset = heap->mem_alloc;
	    if (need_more - need_size >= H5HL_SIZEOF_FREE(f)) {
		if (NULL==(fl = H5FL_MALLOC(H5HL_free_t)))
		    HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, (size_t)(-1), "memory allocation failed");
		fl->offset = heap->mem_alloc + need_size;
		fl->size = need_more - need_size;
		assert (fl->offset==H5HL_ALIGN (fl->offset));
		assert (fl->size==H5HL_ALIGN (fl->size));
		fl->prev = NULL;
		fl->next = heap->freelist;
		if (heap->freelist) heap->freelist->prev = fl;
		heap->freelist = fl;
#ifdef H5HL_DEBUG
	    } else if (H5DEBUG(HL) && need_more > need_size) {
		fprintf(H5DEBUG(HL),
			"H5HL_insert: lost %lu bytes at line %d\n",
			(unsigned long)(need_more - need_size), __LINE__);
#endif
	    }
	}

#ifdef H5HL_DEBUG
	if (H5DEBUG(HL)) {
	    fprintf(H5DEBUG(HL),
		    "H5HL: resize mem buf from %lu to %lu bytes\n",
		    (unsigned long)(heap->mem_alloc),
		    (unsigned long)(heap->mem_alloc + need_more));
	}
#endif
	old_size = heap->mem_alloc;
	heap->mem_alloc += need_more;
	heap->chunk = H5FL_BLK_REALLOC(heap_chunk,heap->chunk,
				   (sizeof_hdr + heap->mem_alloc));
	if (NULL==heap->chunk)
	    HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, (size_t)(-1), "memory allocation failed");

	/* clear new section so junk doesn't appear in the file */
	HDmemset(heap->chunk + sizeof_hdr + old_size, 0, need_more);
    }
    /*
     * Copy the data into the heap
     */
    HDmemcpy(heap->chunk + sizeof_hdr + offset, buf, buf_size);

    /* Set return value */
    ret_value=offset;

done:
    if (heap && H5AC_unprotect(f, dxpl_id, H5AC_LHEAP, addr, heap, FALSE) != SUCCEED)
        HDONE_ERROR(H5E_HEAP, H5E_PROTECT, (size_t)(-1), "unable to release object header");

    FUNC_LEAVE_NOAPI(ret_value);
}

#ifdef NOT_YET

/*-------------------------------------------------------------------------
 * Function:	H5HL_write
 *
 * Purpose:	Writes (overwrites) the object (or part of object) stored
 *		in BUF to the heap at file address ADDR in file F.  The
 *		writing begins at byte offset OFFSET from the beginning of
 *		the heap and continues for SIZE bytes.
 *
 *		Do not partially write an object to create it;	the first
 *		write for an object must be for the entire object.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		matzke@llnl.gov
 *		Jul 16 1997
 *
 * Modifications:
 *		Robb Matzke, 1999-07-28
 *		The ADDR argument is passed by value.
 *-------------------------------------------------------------------------
 */
static herr_t
H5HL_write(H5F_t *f, hid_t dxpl_id, haddr_t addr, size_t offset, size_t size, const void *buf)
{
    H5HL_t *heap = NULL;
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5HL_write, FAIL);

    /* check arguments */
    assert(f);
    assert(H5F_addr_defined(addr));
    assert(buf);
    assert (offset==H5HL_ALIGN (offset));

    if (0==(f->intent & H5F_ACC_RDWR))
	HGOTO_ERROR (H5E_HEAP, H5E_WRITEERROR, FAIL, "no write intent on file");

    if (NULL == (heap = H5AC_protect(f, dxpl_id, H5AC_LHEAP, addr, NULL, NULL, H5AC_WRITE)))
	HGOTO_ERROR(H5E_HEAP, H5E_PROTECT, FAIL, "unable to load heap");

    assert(offset < heap->mem_alloc);
    assert(offset + size <= heap->mem_alloc);

    heap->cache_info.is_dirty=TRUE;
    HDmemcpy(heap->chunk + H5HL_SIZEOF_HDR(f) + offset, buf, size);

done:
    if (heap && H5AC_unprotect(f, dxpl_id, H5AC_LHEAP, addr, heap, FALSE) != SUCCEED && ret_value != FAIL)
        HDONE_ERROR(H5E_HEAP, H5E_PROTECT, FAIL, "unable to release object header");

    FUNC_LEAVE_NOAPI(ret_value);
}
#endif /* NOT_YET */


/*-------------------------------------------------------------------------
 * Function:	H5HL_remove
 *
 * Purpose:	Removes an object or part of an object from the heap at
 *		address ADDR of file F.	 The object (or part) to remove
 *		begins at byte OFFSET from the beginning of the heap and
 *		continues for SIZE bytes.
 *
 *		Once part of an object is removed, one must not attempt
 *		to access that part.  Removing the beginning of an object
 *		results in the object OFFSET increasing by the amount
 *		truncated.  Removing the end of an object results in
 *		object truncation.  Removing the middle of an object results
 *		in two separate objects, one at the original offset and
 *		one at the first offset past the removed portion.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		matzke@llnl.gov
 *		Jul 16 1997
 *
 * Modifications:
 *		Robb Matzke, 1999-07-28
 *		The ADDR argument is passed by value.
 *-------------------------------------------------------------------------
 */
herr_t
H5HL_remove(H5F_t *f, hid_t dxpl_id, haddr_t addr, size_t offset, size_t size)
{
    H5HL_t		*heap = NULL;
    H5HL_free_t		*fl = NULL, *fl2 = NULL;
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5HL_remove, FAIL);

    /* check arguments */
    assert(f);
    assert(H5F_addr_defined(addr));
    assert(size > 0);
    assert (offset==H5HL_ALIGN (offset));

    if (0==(f->intent & H5F_ACC_RDWR))
	HGOTO_ERROR (H5E_HEAP, H5E_WRITEERROR, FAIL, "no write intent on file");

    size = H5HL_ALIGN (size);

    if (NULL == (heap = H5AC_protect(f, dxpl_id, H5AC_LHEAP, addr, NULL, NULL, H5AC_WRITE)))
	HGOTO_ERROR(H5E_HEAP, H5E_PROTECT, FAIL, "unable to load heap");

    assert(offset < heap->mem_alloc);
    assert(offset + size <= heap->mem_alloc);

    fl = heap->freelist;
    heap->cache_info.is_dirty=TRUE;

    /*
     * Check if this chunk can be prepended or appended to an already
     * free chunk.  It might also fall between two chunks in such a way
     * that all three chunks can be combined into one.
     */
    while (fl) {
	if (offset + size == fl->offset) {
	    fl->offset = offset;
	    fl->size += size;
	    assert (fl->offset==H5HL_ALIGN (fl->offset));
	    assert (fl->size==H5HL_ALIGN (fl->size));
	    fl2 = fl->next;
	    while (fl2) {
		if (fl2->offset + fl2->size == fl->offset) {
		    fl->offset = fl2->offset;
		    fl->size += fl2->size;
		    assert (fl->offset==H5HL_ALIGN (fl->offset));
		    assert (fl->size==H5HL_ALIGN (fl->size));
		    fl2 = H5HL_remove_free(heap, fl2);
		    HGOTO_DONE(SUCCEED);
		}
		fl2 = fl2->next;
	    }
	    HGOTO_DONE(SUCCEED);

	} else if (fl->offset + fl->size == offset) {
	    fl->size += size;
	    fl2 = fl->next;
	    assert (fl->size==H5HL_ALIGN (fl->size));
	    while (fl2) {
		if (fl->offset + fl->size == fl2->offset) {
		    fl->size += fl2->size;
		    assert (fl->size==H5HL_ALIGN (fl->size));
		    fl2 = H5HL_remove_free(heap, fl2);
		    HGOTO_DONE(SUCCEED);
		}
		fl2 = fl2->next;
	    }
	    HGOTO_DONE(SUCCEED);
	}
	fl = fl->next;
    }

    /*
     * The amount which is being removed must be large enough to
     * hold the free list data.	 If not, the freed chunk is forever
     * lost.
     */
    if (size < H5HL_SIZEOF_FREE(f)) {
#ifdef H5HL_DEBUG
	if (H5DEBUG(HL)) {
	    fprintf(H5DEBUG(HL), "H5HL: lost %lu bytes\n",
		    (unsigned long) size);
	}
#endif
	HGOTO_DONE(SUCCEED);
    }
    /*
     * Add an entry to the free list.
     */
    if (NULL==(fl = H5FL_MALLOC(H5HL_free_t)))
	HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");
    fl->offset = offset;
    fl->size = size;
    assert (fl->offset==H5HL_ALIGN (fl->offset));
    assert (fl->size==H5HL_ALIGN (fl->size));
    fl->prev = NULL;
    fl->next = heap->freelist;
    if (heap->freelist)
        heap->freelist->prev = fl;
    heap->freelist = fl;

done:
    if (heap && H5AC_unprotect(f, dxpl_id, H5AC_LHEAP, addr, heap, FALSE) != SUCCEED)
        HDONE_ERROR(H5E_HEAP, H5E_PROTECT, FAIL, "unable to release object header");

    FUNC_LEAVE_NOAPI(ret_value);
}


/*-------------------------------------------------------------------------
 * Function:	H5HL_delete
 *
 * Purpose:	Deletes a local heap from disk, freeing disk space used.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *		koziol@ncsa.uiuc.edu
 *		Mar 22 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5HL_delete(H5F_t *f, hid_t dxpl_id, haddr_t addr)
{
    H5HL_t	*heap = NULL;
    size_t      sizeof_hdr;     /* Cache H5HL header size for file */
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5HL_delete, FAIL);

    /* check arguments */
    assert(f);
    assert(H5F_addr_defined(addr));

    /* Check for write access */
    if (0==(f->intent & H5F_ACC_RDWR))
	HGOTO_ERROR (H5E_HEAP, H5E_WRITEERROR, FAIL, "no write intent on file");

    /* Cache this for later */
    sizeof_hdr= H5HL_SIZEOF_HDR(f);

    /* Get heap pointer */
    if (NULL == (heap = H5AC_protect(f, dxpl_id, H5AC_LHEAP, addr, NULL, NULL, H5AC_WRITE)))
	HGOTO_ERROR(H5E_HEAP, H5E_CANTLOAD, FAIL, "unable to load heap");

    /* Check if the heap is contiguous on disk */
    assert(!H5F_addr_overflow(addr,sizeof_hdr));
    if(H5F_addr_eq(heap->addr,addr+sizeof_hdr)) {
        /* Free the contiguous local heap in one call */
        H5_CHECK_OVERFLOW(sizeof_hdr+heap->disk_alloc,size_t,hsize_t);
        if (H5MF_xfree(f, H5FD_MEM_LHEAP, dxpl_id, addr, (hsize_t)(sizeof_hdr+heap->disk_alloc))<0)
            HGOTO_ERROR(H5E_HEAP, H5E_CANTFREE, FAIL, "unable to free contiguous local heap");
    } /* end if */
    else {
        /* Free the local heap's header */
        H5_CHECK_OVERFLOW(sizeof_hdr,size_t,hsize_t);
        if (H5MF_xfree(f, H5FD_MEM_LHEAP, dxpl_id, addr, (hsize_t)sizeof_hdr)<0)
            HGOTO_ERROR(H5E_HEAP, H5E_CANTFREE, FAIL, "unable to free local heap header");

        /* Free the local heap's data */
        H5_CHECK_OVERFLOW(heap->disk_alloc,size_t,hsize_t);
        if (H5MF_xfree(f, H5FD_MEM_LHEAP, dxpl_id, heap->addr, (hsize_t)heap->disk_alloc)<0)
            HGOTO_ERROR(H5E_HEAP, H5E_CANTFREE, FAIL, "unable to free local heap data");
    } /* end else */

    /* Release the local heap metadata from the cache */
    if (H5AC_unprotect(f, dxpl_id, H5AC_LHEAP, addr, heap, TRUE)<0) {
        heap = NULL;
        HGOTO_ERROR(H5E_HEAP, H5E_PROTECT, FAIL, "unable to release local heap");
    }
    heap = NULL;

done:
    if (heap && H5AC_unprotect(f, dxpl_id, H5AC_LHEAP, addr, heap, FALSE)<0)
	HDONE_ERROR(H5E_HEAP, H5E_PROTECT, FAIL, "unable to release local heap");

    FUNC_LEAVE_NOAPI(ret_value);
} /* end H5HL_delete() */

