#include <stdio.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif

void mmaux_dumpmesh(
  char *filename, /* filename (not null terminated!) */
  int *num_verts, /* number of vertices in mesh */
  int *num_edges, /* number of edges in mesh */
  int *num_elems, /* number of elements (quads) in mesh */
  double *coords, /* vertex coordinate array; size: 2*(*num_verts) */
  int *vert_edge, /* vertex-at-edge index array; size: 2*(*num_edges) */
  int *vert_elem, /* vertex-at-element index array; size: 4*(*num_elems) */
  int *edge_elem, /* edge-at-element index array; size: 4*(*num_elems) */
  int num_chars)  /* length of 'filename' parameter */
{
  /* binary mesh magic: "BMSH" in little endian */
  static const int magic_mesh = 0x48534D42;
  char name[512];
  FILE* file;
  int i,n,*tmp;

  /* copy filename and append null terminator */
  for(i = 0; i < num_chars; ++i)
    name[i] = filename[i];
  name[num_chars] = 0;

  /* allocate temporary buffer */
  n = (2*(*num_elems) > *num_edges) ? 4 * (*num_elems) : 2 * (*num_edges);
  tmp = (int*)malloc(4*n);

  /* try to open the file */
  if(file = fopen(name, "wb"))
  {
    /* write magic */
    fwrite(&magic_mesh, 4, 1, file);

    /* write counts */
    fwrite(num_verts, 4, 1, file);
    fwrite(num_edges, 4, 1, file);
    fwrite(num_elems, 4, 1, file);

    /* write vertex coords */
    fwrite(coords, 16, *num_verts, file);

    /* write vertex-at-edge indices */
    for(i = 0; i < 2*(*num_edges); ++i)
      tmp[i] = vert_edge[i] - 1;
    fwrite(tmp, 8, *num_edges, file);

    /* write vertex-at-elem indices */
    for(i = 0; i < 4*(*num_elems); ++i)
      tmp[i] = vert_elem[i] - 1;
    fwrite(tmp, 16, *num_elems, file);

    /* write edge-at-elem indices */
    for(i = 0; i < 4*(*num_elems); ++i)
      tmp[i] = edge_elem[i] - 1;
    fwrite(tmp, 16, *num_elems, file);

    /* close file */
    fclose(file);
  }

  /* free temporary buffer */
  free(tmp);
}

void mmaux_dumpmatrix(
  char *filename,  /* filename (not null terminated!) */
  int *neq,        /* number of rows/columns in the matrix */
  int *nnze,       /* number of non-zero entries in the matrix */
  float *asm_time, /* matrix assembly time in seconds */
  int *row_ptr,    /* CSR row-pointer array; size: *neq + 1 */
  int *col_idx,    /* CSR column-index array; size: *nnze */
  double *data,    /* CSR data array; size: *nnze */
  int num_chars)   /* length of 'filename' parameter */
{
  /* binary matrix magic: "BMTX" in little endian */
  static const int magic_matrix = 0x58544D42;
  char name[512];
  FILE* file;
  int i,n,*tmp;

  /* copy filename and append null terminator */
  for(i = 0; i < num_chars; ++i)
    name[i] = filename[i];
  name[num_chars] = 0;

  /* allocate temporary buffer */
  n = (*neq >= *nnze) ? (*neq) + 1 : (*nnze);
  tmp = (int*)malloc(4*n);

  /* try to open the file */
  if(file = fopen(name, "wb"))
  {
    /* write magic */
    fwrite(&magic_matrix, 4, 1, file);

    /* write counts */
    fwrite(neq, 4, 1, file);
    fwrite(nnze, 4, 1, file);

    /* write assembly time */
    fwrite(asm_time, 4, 1, file);

    /* write row pointer */
    for(i = 0; i <= *neq; ++i)
      tmp[i] = row_ptr[i] - 1;
    fwrite(tmp, 4, (*neq) + 1, file);

    /* write column indices */
    for(i = 0; i < *nnze; ++i)
      tmp[i] = col_idx[i] - 1;
    fwrite(tmp, 4, *nnze, file);

    /* write matrix data */
    fwrite(data, 8, *nnze, file);

    /* close file */
    fclose(file);
  }

  /* free temporary buffer */
  free(tmp);
}

void mmaux_dumpcolor(
  char *filename,  /* filename (not null terminated!) */
  int *ncol,       /* number of colors */
  int *nadj,       /* number of adjacencies in the graph */
  int *row_ptr,    /* row-pointer array; size: *ncol + 1 */
  int *col_idx,    /* column-index array; size: *nadj */
  int num_chars)   /* length of 'filename' parameter */
{
  /* binary matrix magic: "BCPT" in little endian */
  static const int magic_matrix = 0x54504342;
  char name[512];
  FILE* file;
  int i,n,*tmp;

  /* copy filename and append null terminator */
  for(i = 0; i < num_chars; ++i)
    name[i] = filename[i];
  name[num_chars] = 0;

  /* allocate temporary buffer */
  n = (*ncol >= *nadj) ? (*ncol) + 1 : (*nadj);
  tmp = (int*)malloc(4*n);

  /* try to open the file */
  if(file = fopen(name, "wb"))
  {
    /* write magic */
    fwrite(&magic_matrix, 4, 1, file);

    /* write counts */
    fwrite(ncol, 4, 1, file);
    fwrite(nadj, 4, 1, file);

    /* write row pointer */
    for(i = 0; i <= *ncol; ++i)
      tmp[i] = row_ptr[i] - 1;
    fwrite(tmp, 4, (*ncol) + 1, file);

    /* write column indices */
    for(i = 0; i < *nadj; ++i)
      tmp[i] = col_idx[i] - 1;
    fwrite(tmp, 4, *nadj, file);

    /* close file */
    fclose(file);
  }

  /* free temporary buffer */
  free(tmp);
}

/* Fortran wrappers with trailing underscores */
void mmaux_dumpmesh_(
  char *filename,
  int *num_verts,
  int *num_edges,
  int *num_elems,
  double *coords,
  int *vert_edge,
  int *vert_elem,
  int *edge_elem,
  int num_chars)
{
  mmaux_dumpmesh(
    filename,
    num_verts,
    num_edges,
    num_elems,
    coords,
    vert_edge,
    vert_elem,
    edge_elem,
    num_chars);
}

void mmaux_dumpmatrix_(
  char *filename,
  int *neq,
  int *nnze,
  float *asm_time,
  int *row_ptr,
  int *col_idx,
  double *data,
  int num_chars)
{
  mmaux_dumpmatrix(
    filename,
    neq,
    nnze,
    asm_time,
    row_ptr,
    col_idx,
    data,
    num_chars);
}

void mmaux_dumpcolor_(
  char *filename,
  int *ncol,
  int *nadj,
  int *row_ptr,
  int *col_idx,
  int num_chars)
{
  mmaux_dumpcolor(
    filename,
    ncol,
    nadj,
    row_ptr,
    col_idx,
    num_chars);
}

