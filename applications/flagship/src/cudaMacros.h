#pragma once

// Minimum/maximum alignments
#define MAXALIGN(a) MIN(a,16)
#define MINALIGN(a) MAX(a,4)
#define MAX(a,b)    ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a) < (b)) ? (a) : (b) )

// Define short-hand IDX-macros for diagonal list
#ifdef DIAGLIST_DEVICE
#if (DIAGLIST_DEVICE == AOS)
#define IDX2_DIAGLIST IDX2
#else
#define IDX2_DIAGLIST IDX2T
#endif
#endif

// Define short-hand IDX-macros for edge list
#ifdef EDGELIST_DEVICE
#if (EDGELIST_DEVICE == AOS)
#define IDX2_EDGELIST IDX2
#else
#define IDX2_EDGELIST(a,i1,i2,n1,n2) a[((i1-1)/2) * 2*n2 + 2*(i2-1) + ((i1&1) == 0 ? 1 : 0)]
#endif
#endif

// Define short-hand IDX-macros for coefficients at diagonal
#ifdef COEFFSATDIAG_DEVICE
#if (COEFFSATDIAG_DEVICE == AOS)
#define IDX2_COEFFSATDIAG IDX2
#else
#define IDX2_COEFFSATDIAG IDX2T
#endif
#endif

// Define short-hand IDX-macros for coefficients at edges
#ifdef COEFFSATEDGE_DEVICE
#if (COEFFSATEDGE_DEVICE == AOS)
#define IDX3_COEFFSATEDGE IDX3
#else
#define IDX3_COEFFSATEDGE IDX3T
#endif
#endif

inline
static void prepare_baseline(const cudaDeviceProp *devProp,
			     int nitems,
			     int *nitems_per_thread,
			     int threads_per_cta,
			     int *blocks,
			     int *threads,
			     int *nitems_processed) {
  
  if (threads_per_cta > 0) {
    // Adopt number of threads per kernel
    *threads = threads_per_cta;
    
    // Compute number of items per kernel
    *nitems_per_thread = MAX((*nitems_per_thread),
			     nitems/(threads_per_cta*devProp->maxGridSize[0])+1);
    
    // Compute number of items per kernel
    *blocks =(nitems/(*nitems_per_thread)+(*threads)-1)/(*threads);
    
    // Compute total number of items processed
    *nitems_processed = (*nitems_per_thread)*threads_per_cta*(*blocks);
  }
  else {
    *threads = 0;
    *blocks  = 0;
    *nitems_processed = 0;
  }
};

inline
static void prepare_cudaDMA(const cudaDeviceProp *devProp,
			    int nitems,
			    int *nitems_per_thread,
			    int compute_threads_per_cta,
			    int dma_threads_per_ld,
			    int dma_lds,
			    int *blocks,
			    int *threads,
			    int *nitems_processed) {
  
  if (compute_threads_per_cta > 0) {
    // Compute total number of threads per kernel
    *threads = compute_threads_per_cta + dma_lds*dma_threads_per_ld;
    
    // Compute number of items per kernel
    *nitems_per_thread = MAX((*nitems_per_thread),
			     nitems/(compute_threads_per_cta*devProp->maxGridSize[0])+1);
    
    // Compute total number of threads
    *blocks = (nitems/(*nitems_per_thread))/compute_threads_per_cta;
    
    // Compute total number of items processed
    *nitems_processed = *(nitems_per_thread)*compute_threads_per_cta*(*blocks);
  }
  else {
    *threads = 0;
    *blocks  = 0;
    *nitems_processed = 0;
  }
};

