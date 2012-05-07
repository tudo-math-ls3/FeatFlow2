#pragma once

// Minimum/maximum alignments
#define MAXALIGN(a) MIN(a,16)
#define MINALIGN(a) MAX(a,4)
#define MAX(a,b)    ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a) < (b)) ? (a) : (b) )

// Define short-hand IDX-macros for diagonal list
#ifdef DIAGLIST_DEVICE
#if (DIAGLIST_DEVICE == AOS)
#define IDX1_DIAGLIST IDX1
#define IDX2_DIAGLIST IDX2
#define IDX3_DIAGLIST IDX3
#else
#define IDX1_DIAGLIST IDX1T
#define IDX2_DIAGLIST IDX2T
#define IDX3_DIAGLIST IDX3T
#endif
#endif

// Define short-hand IDX-macros for edge list
#ifdef EDGELIST_DEVICE
#if (EDGELIST_DEVICE == AOS)
#define IDX1_EDGELIST IDX1
#define IDX2_EDGELIST IDX2
#define IDX3_EDGELIST IDX3
#else
#define IDX1_EDGELIST IDX1T
#define IDX2_EDGELIST IDX2T
#define IDX3_EDGELIST IDX3T
#endif
#endif

// Define short-hand IDX-macros for coefficients at diagonal
#ifdef COEFFSATDIAG_DEVICE
#if (COEFFSATDIAG_DEVICE == AOS)
#define IDX1_COEFFSATDIAG IDX1
#define IDX2_COEFFSATDIAG IDX2
#define IDX3_COEFFSATDIAG IDX3
#else
#define IDX1_COEFFSATDIAG IDX1T
#define IDX2_COEFFSATDIAG IDX2T
#define IDX3_COEFFSATDIAG IDX3T
#endif
#endif

// Define short-hand IDX-macros for coefficients at edges
#ifdef COEFFSATEDGE_DEVICE
#if (COEFFSATEDGE_DEVICE == AOS)
#define IDX1_COEFFSATEDGE IDX1
#define IDX2_COEFFSATEDGE IDX2
#define IDX3_COEFFSATEDGE IDX3
#else
#define IDX1_COEFFSATEDGE IDX1T
#define IDX2_COEFFSATEDGE IDX2T
#define IDX3_COEFFSATEDGE IDX3T
#endif
#endif

inline
static void prepare_baseline(int nitems,
			     int nitems_per_thread,
			     int threads_per_cta,
			     int *blocks,
			     int *threads,
			     int *nitems_processed) {
  
  if (threads_per_cta > 0) {
    // Adopt number of threads per kernel
    *threads = threads_per_cta;
    
    // Compute number of items per kernel
    const int n_per_thread = MAX(nitems_per_thread,
				 nitems/(threads_per_cta*(65535))+1);
    
    // Compute number of items per kernel
    *blocks =(nitems/nitems_per_thread+(*threads)-1)/(*threads);
    
    // Compute total number of items processed
    *nitems_processed = n_per_thread*threads_per_cta*(*blocks);
  }
  else {
    *threads = 0;
    *blocks  = 0;
    *nitems_processed = 0;
  }
};

inline
static void prepare_cudaDMA(int nitems,
			    int nitems_per_thread,
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
    int n_per_thread = MAX(nitems_per_thread,
			   nitems/(compute_threads_per_cta*(65535))+1);
    
    // Compute total number of threads
    *blocks = (nitems/n_per_thread)/compute_threads_per_cta;
    
    // Compute total number of items processed
    *nitems_processed = n_per_thread*compute_threads_per_cta*(*blocks);
  }
  else {
    *threads = 0;
    *blocks  = 0;
    *nitems_processed = 0;
  }
};

