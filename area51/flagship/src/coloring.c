#include <stdio.h>
#include <stdlib.h>

/*******************************************************************************
 * The macro FNAME converts a C function name to the form that is
 * required when it is called by a FORTRAN program. On many UNIX
 * systems, subprograms called from FORTRAN have an implied underscore
 * character at the ends of their names. This macro takes care of this
 * operating system quirk.
 *******************************************************************************/
#ifdef VMS
#define FNAME(name)	name
#else
#ifdef __APPLE__
#define FNAME(name)	name
#else
#ifdef __STDC__
#define FNAME(name)	name##_
#else
#define FNAME(name)	name/**/_
#endif
#endif
#endif

typedef struct {
  int color;
  int node[2];
  int data[4];
} Edge;

/*
 * Wrapper
 */
int edgeColoring (int neq, int nedge, Edge *edgelist)
{
  return 1;
}

/*
 * Wrapper routine that can be called from Fortran
 *
 * neq  :        number of equations (=vertices)
 * nedge:        number of edges
 * ncolor:       maximum number of colors on input,
 *               actual number of colors on output
 * IedgeList:    pointer to the edge List
 * IedgeListIdx: index pointer separating groups of edges with the same color
 * 
 * Example:
 *   IedgeListIdx = [1,5,10,18];
 *   IedgeList    = [(1,2),(3,4),(5,7),(8,3),...]
 *
 * Edges  1..4  have color C1,
 * edges  5..9  have color C2,
 * edges 10..17 have color C3,...
 */
void FNAME(regroupedgelist)(int *neq, int *nedge, int *ncolor,
			    int **IedgeList, int **IedgeListIdx)
{
  // Allocate list of edges in C-format
  Edge *edgelist = (Edge*) malloc(*nedge*sizeof(Edge));
  
  int iedge,icolor;
  int *d_IedgeList = (int*) IedgeList;
  int *d_IedgeListIdx = (int*) IedgeListIdx;
  
  // Fill list of edges
  for (iedge=0; iedge<(*nedge); ++iedge) {
    // Empty color
    edgelist[iedge].color = 0;
    
    // 0-based vertex numbers
    edgelist[iedge].node[0] = d_IedgeList[6*iedge]-1;
    edgelist[iedge].node[1] = d_IedgeList[6*iedge+1]-1;

    // auxiliary data (do not care about)
    edgelist[iedge].data[0] = d_IedgeList[6*iedge+2];
    edgelist[iedge].data[1] = d_IedgeList[6*iedge+3];
    edgelist[iedge].data[2] = d_IedgeList[6*iedge+4];
    edgelist[iedge].data[3] = d_IedgeList[6*iedge+5];
  }
  
  // Apply edge coloring algorithm
  int ncolors = edgeColoring(*neq, *nedge, edgelist);
  if (ncolors > ncolor) {
    printf("Error: number of colors exceds maximum number of colors!\n");
    return;
  }
  
  // Clear index array
  for (icolor=0; icolor<(*ncolor); ++icolor)
    d_IedgeListIdx[icolor] = 0;

  // Loop over all color groups and over all edges
  int icount=0;
  for (icolor=0; icolor<ncolors; icolor++) {
    d_IedgeListIdx[icolor] = icount+1; // we are 1-based in Fortran

    for (iedge=0; iedge<(*nedge); ++iedge) {
      if (edgelist[iedge].color == icolor) {
	
	// 1-based vertex numbers
	d_IedgeList[6*icount  ] = edgelist[iedge].node[0]+1;
	d_IedgeList[6*icount+1] = edgelist[iedge].node[1]+1;

	// auxiliary data
	d_IedgeList[6*icount+2] = edgelist[iedge].data[0];
	d_IedgeList[6*icount+3] = edgelist[iedge].data[1];
	d_IedgeList[6*icount+4] = edgelist[iedge].data[2];
	d_IedgeList[6*icount+5] = edgelist[iedge].data[3];

	icount++;
      }
    }
  }
  
  // Fill unused groups
  for (icolor=ncolors; icolor<(*ncolor); ++icolor)
    d_IedgeListIdx[icolor] = icount+1;
  
  // Deallocate list of edges
  free(edgelist);
}
