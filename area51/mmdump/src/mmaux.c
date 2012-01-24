#include <stdio.h>
#include <malloc.h>

void mmaux_dumpmesh(
  char *filename,
  int *num_verts,
  int *num_edges,
  int *num_quads,
  double *coords,
  int *vert_idx,
  int *edge_idx,
  int num_chars)
{
  char name[512];
  FILE* file;
  int i,n,*tmp;

  /* copy filename and append null terminator */
  for(i = 0; i < num_chars; ++i)
    name[i] = filename[i];
  name[num_chars] = 0;

  /* allocate temporary buffer */
  n = 4 * (*num_quads);
  tmp = (int*)malloc(4*n);

  /* try to open the file */
  if(file = fopen(name, "wb"))
  {
    /* write magic */
    i = 0x4853454D; /* = "MESH" (little endian) */
    fwrite(&i, 4, 1, file);

    /* write counts */
    fwrite(num_verts, 4, 1, file);
    fwrite(num_edges, 4, 1, file);
    fwrite(num_quads, 4, 1, file);

    /* write vertex coords */
    fwrite(coords, 16, *num_verts, file);

    /* write vertex indices */
    for(i = 0; i < n; ++i)
      tmp[i] = vert_idx[i] - 1;
    fwrite(tmp, 16, *num_quads, file);

    /* write edge indices */
    for(i = 0; i < n; ++i)
      tmp[i] = edge_idx[i] - 1;
    fwrite(tmp, 16, *num_quads, file);

    /* close file */
    fclose(file);
  }

  /* free temporary buffer */
  free(tmp);
}

void mmaux_dumpmatrix(
  char *filename,
  int *neq,
  int *nnze,
  float *asm_time,
  int *row_ptr,
  int *col_idx,
  double *data,
  int num_chars)
{
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
    i = 0x5852544D; /* = "MTRX" (little endian) */
    fwrite(&i, 4, 1, file);

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

/* Fortran wrappers (trailing underscore) */
void mmaux_dumpmesh_(
  char *filename, int *num_verts, int *num_edges, int *num_quads,
  double *coords, int *vert_idx, int *edge_idx, int num_chars)
{
  mmaux_dumpmesh(filename,num_verts,num_edges,num_quads,coords,vert_idx,
    edge_idx,num_chars);
}

void mmaux_dumpmatrix_(
  char *filename, int *neq, int *nnze, float *asm_time, int *row_ptr,
  int *col_idx, double *data, int num_chars)
{
  mmaux_dumpmatrix(filename,neq,nnze,asm_time,row_ptr,col_idx,data,num_chars);
}
