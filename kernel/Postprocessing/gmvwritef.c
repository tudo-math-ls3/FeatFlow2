/*
  Copyright, 1991, The Regents of the University of
  California.  This software was produced under a U. S.
  Government contract (W-7405-ENG-36) by the Los Alamos
  National Laboratory, which is operated by the
  University of California for the U. S. Department of
  Energy.  The U. S. Government is licensed to use,
  reproduce, and distribute this software.  Permission
  is granted to the public to copy and use this software
  without charge, provided that this Notice and any statement
  of authorship are reproduced on all copies.  Neither the
  Government nor the University makes any warranty, express
  or implied, or assumes any liability or responsibility for
  the use of this software.


  Software Author: Kevin L. Bolling

  Numerous changes made by Jeff Hinrichs
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "gmvwrite.c"

char charbuffer[90];

void *fortchartoc(char* filnam);
void *fortchartocfield(char* filnam,int length);

/* ---------------------------------------------------------------- */

void fgmvwrite_openfile(char* filnam)
{
  gmvwrite_openfile(fortchartoc(filnam));
}

void FGMVWRITE_OPENFILE(char* filnam)
{
  gmvwrite_openfile(fortchartoc(filnam));
}

fgmvwrite_openfile_(char* filnam)
{
  gmvwrite_openfile(fortchartoc(filnam));
}

fgmvwrite_openfile__(char* filnam)
{
  gmvwrite_openfile(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_openfile_ascii(char* filnam)
{
  gmvwrite_openfile_ascii(fortchartoc(filnam));
}

void FGMVWRITE_OPENFILE_ASCII(char* filnam)
{
  gmvwrite_openfile_ascii(fortchartoc(filnam));
}

void fgmvwrite_openfile_ascii_(char* filnam)
{
  gmvwrite_openfile_ascii(fortchartoc(filnam));
}

void fgmvwrite_openfile_ascii__(char* filnam)
{
  gmvwrite_openfile_ascii(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_openfile_ir(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_ir(fortchartoc(filnam), *isize, *rsize);
}

void FGMVWRITE_OPENFILE_ir(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_ir(fortchartoc(filnam), *isize, *rsize);
}

void fgmvwrite_openfile_ir_(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_ir(fortchartoc(filnam), *isize, *rsize);
}

void fgmvwrite_openfile_ir__(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_ir(fortchartoc(filnam), *isize, *rsize);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_openfile_cxir(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_cxir(fortchartoc(filnam), *isize, *rsize);
}

void FGMVWRITE_OPENFILE_cxir(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_cxir(fortchartoc(filnam), *isize, *rsize);
}

void fgmvwrite_openfile_cxir_(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_cxir(fortchartoc(filnam), *isize, *rsize);
}

void fgmvwrite_openfile_cxir__(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_cxir(fortchartoc(filnam), *isize, *rsize);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_openfile_ir_ascii(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_ir_ascii(fortchartoc(filnam), *isize, *rsize);
}

void FGMVWRITE_OPENFILE_IR_ASCII(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_ir_ascii(fortchartoc(filnam), *isize, *rsize);
}

void fgmvwrite_openfile_ir_ascii_(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_ir_ascii(fortchartoc(filnam), *isize, *rsize);
}

void fgmvwrite_openfile_ir_ascii__(char *filnam, int *isize, int *rsize)
{
  gmvwrite_openfile_ir_ascii(fortchartoc(filnam), *isize, *rsize);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_closefile()
{
  gmvwrite_closefile();
}

void FGMVWRITE_CLOSEFILE()
{
  gmvwrite_closefile();
}

void fgmvwrite_closefile_()
{
  gmvwrite_closefile();
}

void fgmvwrite_closefile__()
{
  gmvwrite_closefile();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_nodes_fromfile(char filnam[], long nndes)
{
  gmvwrite_nodes_fromfile(fortchartoc(filnam), nndes);
}

void FGMVWRITE_NODES_FROMFILE(char filnam[], long nndes)
{
  gmvwrite_nodes_fromfile(fortchartoc(filnam), nndes);
}

void fgmvwrite_nodes_fromfile_(char filnam[], long nndes)
{
  gmvwrite_nodes_fromfile(fortchartoc(filnam), nndes);
}

void fgmvwrite_nodes_fromfile__(char filnam[], long nndes)
{
  gmvwrite_nodes_fromfile(fortchartoc(filnam), nndes);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_node_data(char *nndes, char *x, char *y, char *z)
{
  gmvwrite_node_data(nndes, x, y, z);
}

void FGMVWRITE_NODE_DATA(char *nndes, char *x, char *y, char *z)
{
  gmvwrite_node_data(nndes, x, y, z);
}

void fgmvwrite_node_data_(char *nndes, char *x, char *y, char *z)
{
  gmvwrite_node_data(nndes, x, y, z);
}

void fgmvwrite_node_data__(char *nndes, char *x, char *y, char *z)
{
  gmvwrite_node_data(nndes, x, y, z);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_node_data_struct(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_node_data_struct(nxv, nyv, nzv, x, y, z);
}

void FGMVWRITE_NODE_DATA_STRUCT(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_node_data_struct(nxv, nyv, nzv, x, y, z);
}

void fgmvwrite_node_data_struct_(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_node_data_struct(nxv, nyv, nzv, x, y, z);
}

void fgmvwrite_node_data_struct__(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_node_data_struct(nxv, nyv, nzv, x, y, z);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_node_data_lstruct(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_node_data_lstruct(nxv, nyv, nzv, x, y, z);
}  

void FGMVWRITE_NODE_DATA_LSTRUCT(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_node_data_lstruct(nxv, nyv, nzv, x, y, z);
}  

void fgmvwrite_node_data_lstruct_(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_node_data_lstruct(nxv, nyv, nzv, x, y, z);
}   

void fgmvwrite_node_data_lstruct__(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_node_data_lstruct(nxv, nyv, nzv, x, y, z);
} 


/* ---------------------------------------------------------------- */

void fgmvwrite_node_data_amr(int *nxc, int *nyc, int *nzc, char *x0, char *y0,
			     char *z0, char *dx, char *dy, char *dz)
{
  gmvwrite_node_data_amr(*nxc, *nyc, *nzc, x0, y0, z0, dx, dy, dz);
}

void FGMVWRITE_NODE_DATA_AMR(int *nxc, int *nyc, int *nzc, char *x0, char *y0,
			     char *z0, char *dx, char *dy, char *dz)
{
  gmvwrite_node_data_amr(*nxc, *nyc, *nzc, x0, y0, z0, dx, dy, dz);
}

void fgmvwrite_node_data_amr_(int *nxc, int *nyc, int *nzc, char *x0, char *y0,
			     char *z0, char *dx, char *dy, char *dz)
{
  gmvwrite_node_data_amr(*nxc, *nyc, *nzc, x0, y0, z0, dx, dy, dz);
}

void fgmvwrite_node_data_amr__(int *nxc, int *nyc, int *nzc, char *x0, char *y0,
			     char *z0, char *dx, char *dy, char *dz)
{
  gmvwrite_node_data_amr(*nxc, *nyc, *nzc, x0, y0, z0, dx, dy, dz);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_nodev_fromfile(char filnam[], long nndes)
{
  gmvwrite_nodev_fromfile(fortchartoc(filnam), nndes);
}

void FGMVWRITE_NODEv_FROMFILE(char filnam[], long nndes)
{
  gmvwrite_nodev_fromfile(fortchartoc(filnam), nndes);
}

void fgmvwrite_nodev_fromfile_(char filnam[], long nndes)
{
  gmvwrite_nodev_fromfile(fortchartoc(filnam), nndes);
}

void fgmvwrite_nodev_fromfile__(char filnam[], long nndes)
{
  gmvwrite_nodev_fromfile(fortchartoc(filnam), nndes);
}

/* ---------------------------------------------------------------- */

void fgmvwrite_nodev_data(char *nndes, char *x, char *y, char *z)
{
  gmvwrite_nodev_data(nndes, x, y, z);
}

void FGMVWRITE_NODEV_DATA(char *nndes, char *x, char *y, char *z)
{
  gmvwrite_nodev_data(nndes, x, y, z);
}

void fgmvwrite_nodev_data_(char *nndes, char *x, char *y, char *z)
{
  gmvwrite_nodev_data(nndes, x, y, z);
}

void fgmvwrite_nodev_data__(char *nndes, char *x, char *y, char *z)
{
  gmvwrite_nodev_data(nndes, x, y, z);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_nodev_data_lstruct(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_nodev_data_lstruct(nxv, nyv, nzv, x, y, z);
}  

void FGMVWRITE_NODEV_DATA_LSTRUCT(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_nodev_data_lstruct(nxv, nyv, nzv, x, y, z);
}  

void fgmvwrite_nodev_data_lstruct_(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_nodev_data_lstruct(nxv, nyv, nzv, x, y, z);
}   

void fgmvwrite_nodev_data_lstruct__(char *nxv, char *nyv, char *nzv, char *x, char *y, char *z)
{
  gmvwrite_nodev_data_lstruct(nxv, nyv, nzv, x, y, z);
}  


/* ---------------------------------------------------------------- */

void fgmvwrite_cells_amr(char *numcells, char *numtop, char *daughters)
{
  gmvwrite_cells_amr(numcells, numtop, daughters);
}

void  FGMVWRITE_CELLS_AMR(char *numcells, char *numtop, char *daughters)
{
  gmvwrite_cells_amr(numcells, numtop, daughters);
}

void fgmvwrite_cells_amr_(char *numcells, char *numtop, char *daughters)
{
  gmvwrite_cells_amr(numcells, numtop, daughters);
}

void fgmvwrite_cells_amr__(char *numcells, char *numtop, char *daughters)
{
  gmvwrite_cells_amr(numcells, numtop, daughters);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_cell_header(char *ncells)
{
  gmvwrite_cell_header(ncells);
}

void FGMVWRITE_CELL_HEADER(char *ncells)
{
  gmvwrite_cell_header(ncells);
}

void fgmvwrite_cell_header_(char *ncells)
{
  gmvwrite_cell_header(ncells);
}

void fgmvwrite_cell_header__(char *ncells)
{
  gmvwrite_cell_header(ncells);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_cells_fromfile(char filnam[], long ncells)
{
  gmvwrite_cells_fromfile(fortchartoc(filnam),ncells);
}

void FGMVWRITE_CELLS_FROMFILE(char filnam[], long ncells)
{
  gmvwrite_cells_fromfile(fortchartoc(filnam),ncells);
}

void fgmvwrite_cells_fromfile_(char filnam[], long ncells)
{
  gmvwrite_cells_fromfile(fortchartoc(filnam),ncells);
}

void fgmvwrite_cells_fromfile__(char filnam[], long ncells)
{
  gmvwrite_cells_fromfile(fortchartoc(filnam),ncells);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_cell_type(char cell_type[], int *nverts, char *nodes)
{
  gmvwrite_cell_type(fortchartocfield(cell_type,8), *nverts, nodes);
}

void FGMVWRITE_CELL_TYPE(char cell_type[], int *nverts, char *nodes)
{
  gmvwrite_cell_type(fortchartocfield(cell_type,8), *nverts, nodes);
}

void fgmvwrite_cell_type_(char cell_type[], int *nverts, char *nodes)
{
  gmvwrite_cell_type(fortchartocfield(cell_type,8), *nverts, nodes);
}

void fgmvwrite_cell_type__(char cell_type[], int *nverts, char *nodes)
{
  gmvwrite_cell_type(fortchartocfield(cell_type,8), *nverts, nodes);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_general_cell_type(char cell_type[], int nverts[], int *nfaces, char *nodeids)
{
  gmvwrite_general_cell_type(fortchartocfield(cell_type,8), nverts, *nfaces, nodeids);
} 

void FGMVWRITE_GENERAL_CELL_TYPE(char cell_type[], int nverts[], int *nfaces, char *nodeids)
{
  gmvwrite_general_cell_type(fortchartocfield(cell_type,8), nverts, *nfaces, nodeids);
} 

void fgmvwrite_general_cell_type_(char cell_type[], int nverts[], int *nfaces, char *nodeids)
{
  gmvwrite_general_cell_type(fortchartocfield(cell_type,8), nverts, *nfaces, nodeids);
} 

void fgmvwrite_general_cell_type__(char cell_type[], int nverts[], int *nfaces, char *nodeids)
{
  gmvwrite_general_cell_type(fortchartocfield(cell_type,8), nverts, *nfaces, nodeids);
} 


/* ---------------------------------------------------------------- */

void fgmvwrite_faces_fromfile(char filnam[], long nfaces)
{
  gmvwrite_faces_fromfile(fortchartoc(filnam), nfaces);
}

void FGMVWRITE_FACES_FROMFILE(char filnam[], long nfaces)
{
  gmvwrite_faces_fromfile(fortchartoc(filnam), nfaces);
}

void fgmvwrite_faces_fromfile_(char filnam[], long nfaces)
{
  gmvwrite_faces_fromfile(fortchartoc(filnam), nfaces);
}

void fgmvwrite_faces_fromfile__(char filnam[], long nfaces)
{
  gmvwrite_faces_fromfile(fortchartoc(filnam), nfaces);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_face_header(char *nfaces, char *ncells)
{
  gmvwrite_face_header(nfaces, ncells);
}

void FGMVWRITE_FACE_HEADER(char *nfaces, char *ncells)
{
  gmvwrite_face_header(nfaces, ncells);
}

void fgmvwrite_face_header_(char *nfaces, char *ncells)
{
  gmvwrite_face_header(nfaces, ncells);
}

void fgmvwrite_face_header__(char *nfaces, char *ncells)
{
  gmvwrite_face_header(nfaces, ncells);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_face_data(int *nverts, char *nodeids, char *cellid1, char *cellid2)
{
  gmvwrite_face_data(*nverts, nodeids, cellid1, cellid2);
}  

void FGMVWRITE_FACE_DATA(int *nverts, char *nodeids, char *cellid1, char *cellid2)
{
  gmvwrite_face_data(*nverts, nodeids, cellid1, cellid2);
}  

void fgmvwrite_face_data_(int *nverts, char *nodeids, char *cellid1, char *cellid2)
{
  gmvwrite_face_data(*nverts, nodeids, cellid1, cellid2);
}     

void fgmvwrite_face_data__(int *nverts, char *nodeids, char *cellid1, char *cellid2)
{
  gmvwrite_face_data(*nverts, nodeids, cellid1, cellid2);
}   


/* ---------------------------------------------------------------- */

void fgmvwrite_vfaces_fromfile(char filnam[], long nfaces)
{
  gmvwrite_vfaces_fromfile(fortchartoc(filnam), nfaces);
}

void FGMVWRITE_VFACES_FROMFILE(char filnam[], long nfaces)
{
  gmvwrite_vfaces_fromfile(fortchartoc(filnam), nfaces);
}

void fgmvwrite_vfaces_fromfile_(char filnam[], long nfaces)
{
  gmvwrite_vfaces_fromfile(fortchartoc(filnam), nfaces);
}

void fgmvwrite_vfaces_fromfile__(char filnam[], long nfaces)
{
  gmvwrite_vfaces_fromfile(fortchartoc(filnam), nfaces);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_vface_header(char *nfaces)
{
  gmvwrite_vface_header(nfaces);
}

void FGMVWRITE_VFACE_HEADER(char *nfaces)
{
  gmvwrite_vface_header(nfaces);
}

void fgmvwrite_vface_header_(char *nfaces)
{
  gmvwrite_vface_header(nfaces);
}

void fgmvwrite_vface_header__(char *nfaces)
{
  gmvwrite_vface_header(nfaces);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_vface_data(int *nverts, int *facepe, char *oppface, 
			  int *oppfacepe, char *cellid, char *nodeids)
{
  gmvwrite_vface_data(*nverts, *facepe, oppface, *oppfacepe, cellid, nodeids);
}

void FGMVWRITE_VFACE_DATA(int *nverts, int *facepe, char *oppface, 
			  int *oppfacepe, char *cellid, char *nodeids)
{
  gmvwrite_vface_data(*nverts, *facepe, oppface, *oppfacepe, cellid, nodeids);
}

void fgmvwrite_vface_data_(int *nverts, int *facepe, char *oppface, 
			   int *oppfacepe, char *cellid, char *nodeids)
{
  gmvwrite_vface_data(*nverts, *facepe, oppface, *oppfacepe, cellid, nodeids);
}

void fgmvwrite_vface_data__(int *nverts, int *facepe, char *oppface, 
			    int *oppfacepe, char *cellid, char *nodeids)
{
  gmvwrite_vface_data(*nverts, *facepe, oppface, *oppfacepe, cellid, nodeids);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_material_fromfile(char filnam[])
{
  gmvwrite_material_fromfile(fortchartoc(filnam));
}

void FGMVWRITE_MATERIAL_FROMFILE(char filnam[])
{
  gmvwrite_material_fromfile(fortchartoc(filnam));
}

void fgmvwrite_material_fromfile_(char filnam[])
{
  gmvwrite_material_fromfile(fortchartoc(filnam));
}

void fgmvwrite_material_fromfile__(char filnam[])
{
  gmvwrite_material_fromfile(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_material_header(int *nmats, int *data_type)
{
  gmvwrite_material_header(*nmats, *data_type);
}

void FGMVWRITE_MATERIAL_HEADER(int *nmats, int *data_type)
{
  gmvwrite_material_header(*nmats, *data_type);
}

void fgmvwrite_material_header_(int *nmats, int *data_type)
{
  gmvwrite_material_header(*nmats, *data_type);
}

void fgmvwrite_material_header__(int *nmats, int *data_type)
{
  gmvwrite_material_header(*nmats, *data_type);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_material_name(char matname[])
{
  char name[10];
  /* teminate matname with null character */
  strncpy(name, matname,8);
  *(name + 8) = '\0';
  gmvwrite_material_name(name);
}

void FGMVWRITE_MATERIAL_NAME(char matname[])
{
  char name[10];
  /* teminate matname with null character */
  strncpy(name, matname,8);
  *(name + 8) = '\0';
  gmvwrite_material_name(name);
}

void fgmvwrite_material_name_(char matname[])
{
  char name[10];
  /* teminate matname with null character */
  strncpy(name, matname,8);
  *(name + 8) = '\0';
  gmvwrite_material_name(name);
}

void fgmvwrite_material_name__(char matname[])
{
  char name[10];
  /* teminate matname with null character */
  strncpy(name, matname,8);
  *(name + 8) = '\0';
  gmvwrite_material_name(name);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_material_ids(int matids[], int *data_type)
{
  gmvwrite_material_ids(matids, *data_type);
}

void FGMVWRITE_MATERIAL_IDS(int matids[], int *data_type)
{
  gmvwrite_material_ids(matids, *data_type);
}

void fgmvwrite_material_ids_(int matids[], int *data_type)
{
  gmvwrite_material_ids(matids, *data_type);
}

void fgmvwrite_material_ids__(int matids[], int *data_type)
{
  gmvwrite_material_ids(matids, *data_type);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_velocity_data(int *data_type, char *u, char *v, char *w)
{
  gmvwrite_velocity_data(*data_type, u, v, w);
}

void FGMVWRITE_VELOCITY_DATA(int *data_type, char *u, char *v, char *w)
{
  gmvwrite_velocity_data(*data_type, u, v, w);
}

void fgmvwrite_velocity_data_(int *data_type, char *u, char *v, char *w)
{
  gmvwrite_velocity_data(*data_type, u, v, w);
}

void fgmvwrite_velocity_data__(int *data_type, char *u, char *v, char *w)
{
  gmvwrite_velocity_data(*data_type, u, v, w);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_variable_header()
{
  gmvwrite_variable_header();
}

void FGMVWRITE_VARIABLE_HEADER()
{
  gmvwrite_variable_header();
}

void fgmvwrite_variable_header_()
{
  gmvwrite_variable_header();
}

void fgmvwrite_variable_header__()
{
  gmvwrite_variable_header();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_variable_name_data(int *data_type, char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_variable_name_data(*data_type, name, vids);
}

void FGMVWRITE_VARIABLE_NAME_DATA(int *data_type, char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_variable_name_data(*data_type, name, vids);
}

void fgmvwrite_variable_name_data_(int *data_type, char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_variable_name_data(*data_type, name, vids);
}

void fgmvwrite_variable_name_data__(int *data_type, char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_variable_name_data(*data_type, name, vids);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_variable_endvars()
{
  gmvwrite_variable_endvars();
}

void FGMVWRITE_VARIABLE_ENDVARS()
{
  gmvwrite_variable_endvars();
}

void fgmvwrite_variable_endvars_()
{
  gmvwrite_variable_endvars();
}

void fgmvwrite_variable_endvars__()
{
  gmvwrite_variable_endvars();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_flag_fromfile(char filnam[])
{
  gmvwrite_flag_fromfile(fortchartoc(filnam));
}

void FGMVWRITE_FLAG_FROMFILE(char filnam[])
{
  gmvwrite_flag_fromfile(fortchartoc(filnam));
}

void fgmvwrite_flag_fromfile_(char filnam[])
{
  gmvwrite_flag_fromfile(fortchartoc(filnam));
}

void fgmvwrite_flag_fromfile__(char filnam[])
{
  gmvwrite_flag_fromfile(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_flag_header()
{
  gmvwrite_flag_header();
}

void FGMVWRITE_FLAG_HEADER()
{
  gmvwrite_flag_header();
}

void fgmvwrite_flag_header_()
{
  gmvwrite_flag_header();
}

void fgmvwrite_flag_header__()
{
  gmvwrite_flag_header();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_flag_name(char flagname[], int *numtypes, int *data_type)
{
  char name[10];
  strncpy(name, flagname, 8);
  /*  terminate flagname with null character  */
  *(name+8) = '\0';
  gmvwrite_flag_name(name, *numtypes, *data_type);
}

void FGMVWRITE_FLAG_NAME(char flagname[], int *numtypes, int *data_type)
{
  char name[10];
  strncpy(name, flagname, 8);
  /*  terminate flagname with null character  */
  *(name+8) = '\0';
  gmvwrite_flag_name(name, *numtypes, *data_type);
}

void fgmvwrite_flag_name_(char flagname[], int *numtypes, int *data_type)
{
  char name[10];
  strncpy(name, flagname, 8);
  /*  terminate flagname with null character  */
  *(name+8) = '\0';
  gmvwrite_flag_name(name, *numtypes, *data_type);
}

void fgmvwrite_flag_name__(char flagname[], int *numtypes, int *data_type)
{
  char name[10];
  strncpy(name, flagname, 8);
  /*  terminate flagname with null character  */
  *(name+8) = '\0';
  gmvwrite_flag_name(name, *numtypes, *data_type);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_flag_subname(char subname[])
{
  char type[10];
  strncpy(type, subname, 8);
  /*  terminate subname with null character  */
  *(type+8) = '\0';
  gmvwrite_flag_subname(type);
}

void FGMVWRITE_FLAG_SUBNAME(char subname[])
{
  char type[10];
  strncpy(type, subname, 8);
  /*  terminate subname with null character  */
  *(type+8) = '\0';
  gmvwrite_flag_subname(type);
}

void fgmvwrite_flag_subname_(char subname[])
{
  char type[10];
  strncpy(type, subname, 8);
  /*  terminate subname with null character  */
  *(type+8) = '\0';
  gmvwrite_flag_subname(type);
}

void fgmvwrite_flag_subname__(char subname[])
{
  char type[10];
  strncpy(type, subname, 8);
  /*  terminate subname with null character  */
  *(type+8) = '\0';
  gmvwrite_flag_subname(type);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_flag_data(int *data_type,int flagdata[])
{
  gmvwrite_flag_data(*data_type, flagdata);
}

void FGMVWRITE_FLAG_DATA(int *data_type,int flagdata[])
{
  gmvwrite_flag_data(*data_type, flagdata);
}

void fgmvwrite_flag_data_(int *data_type, int flagdata[])
{
  gmvwrite_flag_data(*data_type, flagdata);
}

void fgmvwrite_flag_data__(int *data_type, int flagdata[])
{
  gmvwrite_flag_data(*data_type, flagdata);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_flag_endflag()
{
  gmvwrite_flag_endflag();
}

void FGMVWRITE_FLAG_ENDFLAG()
{
  gmvwrite_flag_endflag();
}

void fgmvwrite_flag_endflag_()
{
  gmvwrite_flag_endflag();
}

void fgmvwrite_flag_endflag__()
{
  gmvwrite_flag_endflag();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_polygons_fromfile(char filnam[])
{
  gmvwrite_polygons_fromfile(fortchartoc(filnam));
}

void FGMVWRITE_POLYGONS_FROMFILE(char filnam[])
{
  gmvwrite_polygons_fromfile(fortchartoc(filnam));
}

void fgmvwrite_polygons_fromfile_(char filnam[])
{
  gmvwrite_polygons_fromfile(fortchartoc(filnam));
}

void fgmvwrite_polygons_fromfile__(char filnam[])
{
  gmvwrite_polygons_fromfile(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_polygons_header()
{
  gmvwrite_polygons_header();
}

void FGMVWRITE_POLYGONS_HEADER()
{
  gmvwrite_polygons_header();
}

void fgmvwrite_polygons_header_()
{
  gmvwrite_polygons_header();
}

void fgmvwrite_polygons_header__()
{
  gmvwrite_polygons_header();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_polygons_data(int *nverts, int *matnum, char *x, char *y, char *z)
{  
  gmvwrite_polygons_data(*nverts, *matnum, x, y, z);
}

void FGMVWRITE_POLYGONS_DATA(int *nverts, int *matnum, char *x, char *y, char *z)
{  
  gmvwrite_polygons_data(*nverts, *matnum, x, y, z);
}

void fgmvwrite_polygons_data_(int *nverts, int *matnum, char *x, char *y, char *z)
{  
  gmvwrite_polygons_data(*nverts, *matnum, x, y, z);
}

void fgmvwrite_polygons_data__(int *nverts, int *matnum, char *x, char *y, char *z)
{  
  gmvwrite_polygons_data(*nverts, *matnum, x, y, z);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_polygons_endpoly()
{
  gmvwrite_polygons_endpoly();
}

void FGMVWRITE_POLYGONS_ENDPOLY()
{
  gmvwrite_polygons_endpoly();
}

void fgmvwrite_polygons_endpoly_()
{
  gmvwrite_polygons_endpoly();
}

void fgmvwrite_polygons_endpoly__()
{
  gmvwrite_polygons_endpoly();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_tracers_header(int *ntracers, char *u, char *y, char *z)
{
  gmvwrite_tracers_header(*ntracers, u, y, z);
}

void FGMVWRITE_TRACERS_HEADER(int *ntracers, char *u, char *y, char *z)
{
  gmvwrite_tracers_header(*ntracers, u, y, z);
}

void fgmvwrite_tracers_header_(int *ntracers, char *u, char *y, char *z)
{
  gmvwrite_tracers_header(*ntracers, u, y, z);
}

void fgmvwrite_tracers_header__(int *ntracers, char *u, char *y, char *z)
{
  gmvwrite_tracers_header(*ntracers, u, y, z);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_tracers_name_data(int *ntracers, char tracername[], char *data)
{
  char name[10];
  strncpy(name, tracername,8);
  /*  terminate tracername with null character  */
  *(name+8) = '\0';
  gmvwrite_tracers_name_data(*ntracers, name, data);
}

void FGMVWRITE_TRACERS_NAME_DATA(int *ntracers, char tracername[], char *data)
{
  char name[10];
  strncpy(name, tracername,8);
  /*  terminate tracername with null character  */
  *(name+8) = '\0';
  gmvwrite_tracers_name_data(*ntracers, name, data);
}

void fgmvwrite_tracers_name_data_(int *ntracers, char tracername[], char *data)
{
  char name[10];
  strncpy(name, tracername,8);
  /*  terminate tracername with null character  */
  *(name+8) = '\0';
  gmvwrite_tracers_name_data(*ntracers, name, data);
}

void fgmvwrite_tracers_name_data__(int *ntracers, char tracername[], char *data)
{
  char name[10];
  strncpy(name, tracername,8);
  /*  terminate tracername with null character  */
  *(name+8) = '\0';
  gmvwrite_tracers_name_data(*ntracers, name, data);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_tracers_endtrace()
{
  gmvwrite_tracers_endtrace();
}

void FGMVWRITE_TRACERS_ENDTRACE()
{
  gmvwrite_tracers_endtrace();
}

void fgmvwrite_tracers_endtrace_()
{
  gmvwrite_tracers_endtrace();
}

void fgmvwrite_tracers_endtrace__()
{
  gmvwrite_tracers_endtrace();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_probtime(double *ptime)
{
  gmvwrite_probtime(*ptime);
}

void FGMVWRITE_PROBTIME(double *ptime)
{
  gmvwrite_probtime(*ptime);
}

void fgmvwrite_probtime_(double *ptime)
{
  gmvwrite_probtime(*ptime);
}

void fgmvwrite_probtime__(double *ptime)
{
  gmvwrite_probtime(*ptime);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_cycleno(int *number)
{
  gmvwrite_cycleno(*number);
}

void FGMVWRITE_CYCLENO(int *number)
{
  gmvwrite_cycleno(*number);
}

void fgmvwrite_cycleno_(int *number)
{
  gmvwrite_cycleno(*number);
}

void fgmvwrite_cycleno__(int *number)
{
  gmvwrite_cycleno(*number);
}


/* ---------------------------------------------------------------- */

void *fortchartoc(char* filnam){

  int i;


  /*  teminate fname with null character */

  i = 0;
  while ( *(filnam+i) != ' ')
    {
      charbuffer[i] = *(filnam+i);
      i++;
    }
  charbuffer[i] = '\0';

  return charbuffer;

}

/* ---------------------------------------------------------------- */

void *fortchartocfield(char* filnam,int length){

  int i;


  /*  teminate fname with null character */

  for(i = 0; (i < length) && (*(filnam+i) != ' '); i++)
    charbuffer[i] = *(filnam+i);

  for(; i < length; i++)
    charbuffer[i] = ' ';

  charbuffer[++i] = '\0';

  return charbuffer;

}

/* ---------------------------------------------------------------- */

void fgmvwrite_nodeids(char *nodeids)
{
  gmvwrite_nodeids(nodeids);
}  

void FGMVWRITE_NODEIDS(char *nodeids)
{
  gmvwrite_nodeids(nodeids);
}  

void fgmvwrite_nodeids_(char *nodeids)
{
  gmvwrite_nodeids(nodeids);
}     

void fgmvwrite_nodeids__(char *nodeids)
{
  gmvwrite_nodeids(nodeids);
}  

/* ---------------------------------------------------------------- */

void fgmvwrite_nodeids_fromfile(char filnam[])
{
  gmvwrite_nodeids_fromfile(fortchartoc(filnam));
}

void FGMVWRITE_NODEIDS_FROMFILE(char filnam[])
{
  gmvwrite_nodeids_fromfile(fortchartoc(filnam));
}

void fgmvwrite_nodeids_fromfile_(char filnam[])
{
  gmvwrite_nodeids_fromfile(fortchartoc(filnam));
}

void fgmvwrite_nodeids_fromfile__(char filnam[])
{
  gmvwrite_nodeids_fromfile(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_cellids(char *cellids)
{
  gmvwrite_cellids(cellids);
}   

void FGMVWRITE_CELLIDS(char *cellids)
{
  gmvwrite_cellids(cellids);
}   

void fgmvwrite_cellids_(char *cellids)
{
  gmvwrite_cellids(cellids);
}  

void fgmvwrite_cellids__(char *cellids)
{
  gmvwrite_cellids(cellids);
}  


/* ---------------------------------------------------------------- */

void fgmvwrite_cellids_fromfile(char filnam[])
{
  gmvwrite_cellids_fromfile(fortchartoc(filnam));
}

void FGMVWRITE_CELLIDS_FROMFILE(char filnam[])
{
  gmvwrite_cellids_fromfile(fortchartoc(filnam));
}

void fgmvwrite_cellids_fromfile_(char filnam[])
{
  gmvwrite_cellids_fromfile(fortchartoc(filnam));
}

void fgmvwrite_cellids_fromfile__(char filnam[])
{
  gmvwrite_cellids_fromfile(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surface_header(char *nsurf)
{
  gmvwrite_surface_header(nsurf);
}

void FGMVWRITE_SURFACE_HEADER(char *nsurf)
{
  gmvwrite_surface_header(nsurf);
}

void fgmvwrite_surface_header_(char *nsurf)
{
  gmvwrite_surface_header(nsurf);
}

void fgmvwrite_surface_header__(char *nsurf)
{
  gmvwrite_surface_header(nsurf);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surface_data(int *nverts, char *nodeids)
{
  gmvwrite_surface_data(*nverts, nodeids);
} 

void FGMVWRITE_SURFACE_DATA(int *nverts, char *nodeids)
{
  gmvwrite_surface_data(*nverts, nodeids);
} 

void fgmvwrite_surface_data_(int *nverts, char *nodeids)
{
  gmvwrite_surface_data(*nverts, nodeids);
}    

void fgmvwrite_surface_data__(int *nverts, char *nodeids)
{
  gmvwrite_surface_data(*nverts, nodeids);
}   


/* ---------------------------------------------------------------- */

void fgmvwrite_surface_fromfile(char filnam[], long nsrf)
{
  gmvwrite_surface_fromfile(fortchartoc(filnam), nsrf);
}

void FGMVWRITE_surface_FROMFILE(char filnam[], long nsrf)
{
  gmvwrite_surface_fromfile(fortchartoc(filnam), nsrf);
}

void fgmvwrite_surface_fromfile_(char filnam[], long nsrf)
{
  gmvwrite_surface_fromfile(fortchartoc(filnam), nsrf);
}

void fgmvwrite_surface_fromfile__(char filnam[], long nsrf)
{
  gmvwrite_surface_fromfile(fortchartoc(filnam), nsrf);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surfmats(int matids[])
{
  gmvwrite_surfmats(matids);
}

void FGMVWRITE_SURFMATS(int matids[])
{
  gmvwrite_surfmats(matids);
}

void fgmvwrite_surfmats_(int matids[])
{
  gmvwrite_surfmats(matids);
}

void fgmvwrite_surfmats__(int matids[])
{
  gmvwrite_surfmats(matids);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surfvel(char *u, char *v, char *w)
{
  gmvwrite_surfvel(u, v, w);
}

void FGMVWRITE_SURFVEL(char *u, char *v, char *w)
{
  gmvwrite_surfvel(u, v, w);
}

void fgmvwrite_surfvel_(char *u, char *v, char *w)
{
  gmvwrite_surfvel(u, v, w);
}

void fgmvwrite_surfvel__(char *u, char *v, char *w)
{
  gmvwrite_surfvel(u, v, w);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surfvars_header()
{
  gmvwrite_surfvars_header();
}

void FGMVWRITE_SURFVARS_HEADER()
{
  gmvwrite_surfvars_header();
}

void fgmvwrite_surfvars_header_()
{
  gmvwrite_surfvars_header();
}

void fgmvwrite_surfvars_header__()
{
  gmvwrite_surfvars_header();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surfvars_name_data(char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_surfvars_name_data(name, vids);
}

void FGMVWRITE_SURFVARS_NAME_DATA(char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_surfvars_name_data(name, vids);
}

void fgmvwrite_surfvars_name_data_(char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_surfvars_name_data(name, vids);
}

void fgmvwrite_surfvars_name_data__(char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_surfvars_name_data(name, vids);
}

/* ---------------------------------------------------------------- */

void fgmvwrite_surfvars_endsvar()
{
  gmvwrite_surfvars_endsvar();
}

void FGMVWRITE_SURFVARS_ENDSVAR()
{
  gmvwrite_surfvars_endsvar();
}

void fgmvwrite_surfvars_endsvar_()
{
  gmvwrite_surfvars_endsvar();
}

void fgmvwrite_surfvars_endsvar__()
{
  gmvwrite_surfvars_endsvar();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surfflag_header()
{
  gmvwrite_surfflag_header();
}

void FGMVWRITE_SURFFLAG_HEADER()
{
  gmvwrite_surfflag_header();
}

void fgmvwrite_surfflag_header_()
{
  gmvwrite_surfflag_header();
}

void fgmvwrite_surfflag_header__()
{
  gmvwrite_surfflag_header();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surfflag_name(char flagname[], int *numtypes)
{
  char name[10];
  strncpy(name, flagname, 8);
  /*  terminate flagname with null character  */
  *(name+8) = '\0';
  gmvwrite_surfflag_name(name, *numtypes);
}

void FGMVWRITE_SURFFLAG_NAME(char flagname[], int *numtypes)
{
  char name[10];
  strncpy(name, flagname, 8);
  /*  terminate flagname with null character  */
  *(name+8) = '\0';
  gmvwrite_surfflag_name(name, *numtypes);
}

void fgmvwrite_surfflag_name_(char flagname[], int *numtypes)
{
  char name[10];
  strncpy(name, flagname, 8);
  /*  terminate flagname with null character  */
  *(name+8) = '\0';
  gmvwrite_surfflag_name(name, *numtypes);
}

void fgmvwrite_surfflag_name__(char flagname[], int *numtypes)
{
  char name[10];
  strncpy(name, flagname, 8);
  /*  terminate flagname with null character  */
  *(name+8) = '\0';
  gmvwrite_surfflag_name(name, *numtypes);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surfflag_subname(char subname[])
{
  char type[10];
  strncpy(type, subname, 8);
  /*  terminate subname with null character  */
  *(type+8) = '\0';
  gmvwrite_surfflag_subname(type);
}

void FGMVWRITE_SURFFLAG_SUBNAME(char subname[])
{
  char type[10];
  strncpy(type, subname, 8);
  /*  terminate subname with null character  */
  *(type+8) = '\0';
  gmvwrite_surfflag_subname(type);
}

void fgmvwrite_surfflag_subname_(char subname[])
{
  char type[10];
  strncpy(type, subname, 8);
  /*  terminate subname with null character  */
  *(type+8) = '\0';
  gmvwrite_surfflag_subname(type);
}

void fgmvwrite_surfflag_subname__(char subname[])
{
  char type[10];
  strncpy(type, subname, 8);
  /*  terminate subname with null character  */
  *(type+8) = '\0';
  gmvwrite_surfflag_subname(type);
}

/* ---------------------------------------------------------------- */

void fgmvwrite_surfflag_data(int flagdata[])
{
  gmvwrite_surfflag_data(flagdata);
}

void FGMVWRITE_SURFFLAG_DATA(int flagdata[])
{
  gmvwrite_surfflag_data(flagdata);
}

void fgmvwrite_surfflag_data_(int flagdata[])
{
  gmvwrite_surfflag_data(flagdata);
}

void fgmvwrite_surfflag_data__(int flagdata[])
{
  gmvwrite_surfflag_data(flagdata);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surfflag_endflag()
{
  gmvwrite_surfflag_endflag();
}

void FGMVWRITE_SURFFLAG_ENDFLAG()
{
  gmvwrite_surfflag_endflag();
}

void fgmvwrite_surfflag_endflag_()
{
  gmvwrite_surfflag_endflag();
}

void fgmvwrite_surfflag_endflag__()
{
  gmvwrite_surfflag_endflag();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_units_fromfile(char filnam[])
{
  gmvwrite_units_fromfile(fortchartoc(filnam));
}

void FGMVWRITE_UNITS_FROMFILE(char filnam[])
{
  gmvwrite_units_fromfile(fortchartoc(filnam));
}

void fgmvwrite_units_fromfile_(char filnam[])
{
  gmvwrite_units_fromfile(fortchartoc(filnam));
}

void fgmvwrite_units_fromfile__(char filnam[])
{
  gmvwrite_units_fromfile(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_units_header()
{
  gmvwrite_units_header();
}

void FGMVWRITE_UNITS_HEADER()
{
  gmvwrite_units_header();
}

void fgmvwrite_units_header_()
{
  gmvwrite_units_header();
}
 void fgmvwrite_units_header__()
{
  gmvwrite_units_header();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_units_typehdr(int *data_type, int numtypes)
{
  gmvwrite_units_typehdr(*data_type, numtypes);
}     

void FGMVWRITE_UNITS_TYPEHDR(int *data_type, int numtypes)
{
  gmvwrite_units_typehdr(*data_type, numtypes);
}     

void fgmvwrite_units_typehdr_(int *data_type, int numtypes)
{
  gmvwrite_units_typehdr(*data_type, numtypes);
}     

void fgmvwrite_units_typehdr__(int *data_type, int numtypes)
{
  gmvwrite_units_typehdr(*data_type, numtypes);
}     


/* ---------------------------------------------------------------- */

void fgmvwrite_units_name(char fldname[], char unitname[])
{
  char type1[10], type2[20];
  strncpy(type1, fldname, 8);
  strncpy(type2, unitname, 16);
  /*  terminate with null character  */
  *(type1+8) = '\0';
  *(type2+16) = '\0';
  gmvwrite_units_name(type1,type2);
}

void FGMVWRITE_UNITS_NAME(char fldname[], char unitname[])
{
  char type1[10], type2[20];
  strncpy(type1, fldname, 8);
  strncpy(type2, unitname, 16);
  /*  terminate with null character  */
  *(type1+8) = '\0';
  *(type2+16) = '\0';
  gmvwrite_units_name(type1,type2);
}

void fgmvwrite_units_name_(char fldname[], char unitname[])
{
  char type1[10], type2[20];
  strncpy(type1, fldname, 8);
  strncpy(type2, unitname, 16);
  /*  terminate with null character  */
  *(type1+8) = '\0';
  *(type2+16) = '\0';
  gmvwrite_units_name(type1,type2);
}

void fgmvwrite_units_name__(char fldname[], char unitname[])
{
  char type1[10], type2[20];
  strncpy(type1, fldname, 8);
  strncpy(type2, unitname, 16);
  /*  terminate with null character  */
  *(type1+8) = '\0';
  *(type2+16) = '\0';
  gmvwrite_units_name(type1,type2);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_units_endunit()
{
  gmvwrite_units_endunit();
}

void FGMVWRITE_UNITS_ENDUNIT()
{
  gmvwrite_units_endunit();
}

void fgmvwrite_units_endunit_()
{
  gmvwrite_units_endunit();
}

void fgmvwrite_units_endunit__()
{
  gmvwrite_units_endunit();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_vinfo_header()
{
  gmvwrite_vinfo_header();
}

void FGMVWRITE_VINFO_HEADER()
{
  gmvwrite_vinfo_header();
}

void fgmvwrite_vinfo_header_()
{
  gmvwrite_vinfo_header();
}

void fgmvwrite_vinfo_header__()
{
  gmvwrite_vinfo_header();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_vinfo_name_data(int nelem, int nlines, char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_vinfo_name_data(nelem, nlines, name, vids);
}

void FGMVWRITE_VINFO_NAMEDATA(int nelem, int nlines, char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_vinfo_name_data(nelem, nlines, name, vids);
}

void fgmvwrite_vinfo_name_data_(int nelem, int nlines, char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_vinfo_name_data(nelem, nlines, name, vids);
}

void fgmvwrite_vinfo_name_data__(int nelem, int nlines, char varname[], char *vids)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_vinfo_name_data(nelem, nlines, name, vids);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_vinfo_endvinfo()
{
  gmvwrite_vinfo_endvinfo();
}

void FGMVWRITE_VINFO_ENDVINFO()
{
  gmvwrite_vinfo_endvinfo();
}

void fgmvwrite_vinfo_endvinfo_()
{
  gmvwrite_vinfo_endvinfo();
}

void fgmvwrite_vinfo_endvinfo__()
{
  gmvwrite_vinfo_endvinfo();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_traceids(int ntracers, char *traceids)
{
  gmvwrite_traceids(ntracers, traceids);
}     

void FGMVWRITE_TRACEIDS(int ntracers, char *traceids)
{
  gmvwrite_traceids(ntracers, traceids);
}     

void fgmvwrite_traceids_(int ntracers, char *traceids)
{
  gmvwrite_traceids(ntracers, traceids);
}     

void fgmvwrite_traceids__(int ntracers, char *traceids)
{
  gmvwrite_traceids(ntracers, traceids);
}     


/* ---------------------------------------------------------------- */

void fgmvwrite_traceids_fromfile(char filnam[])
{
  gmvwrite_traceids_fromfile(fortchartoc(filnam));
}

void FGMVWRITE_TRACEIDS_FROMFILE(char filnam[])
{
  gmvwrite_traceids_fromfile(fortchartoc(filnam));
}

void fgmvwrite_traceids_fromfile_(char filnam[])
{
  gmvwrite_traceids_fromfile(fortchartoc(filnam));
}

void fgmvwrite_traceids_fromfile__(char filnam[])
{
  gmvwrite_traceids_fromfile(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_faceids(char *faceids)
{
  gmvwrite_faceids(faceids);
}     

void FGMVWRITE_FACEIDS(char *faceids)
{
  gmvwrite_faceids(faceids);
}     

void fgmvwrite_faceids_(char *faceids)
{
  gmvwrite_faceids(faceids);
}     

void fgmvwrite_faceids__(char *faceids)
{
  gmvwrite_faceids(faceids);
}     


/* ---------------------------------------------------------------- */

void fgmvwrite_faceids_fromfile(char filnam[])
{
  gmvwrite_faceids_fromfile(fortchartoc(filnam));
}

void FGMVWRITE_FACEIDS_FROMFILE(char filnam[])
{
  gmvwrite_faceids_fromfile(fortchartoc(filnam));
}

void fgmvwrite_faceids_fromfile_(char filnam[])
{
  gmvwrite_faceids_fromfile(fortchartoc(filnam));
}

void fgmvwrite_faceids_fromfile__(char filnam[])
{
  gmvwrite_faceids_fromfile(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_group_fromfile(char filnam[])
{
  gmvwrite_group_fromfile(fortchartoc(filnam));
}

void FGMVWRITE_GROUP_FROMFILE(char filnam[])
{
  gmvwrite_group_fromfile(fortchartoc(filnam));
}

void fgmvwrite_group_fromfile_(char filnam[])
{
  gmvwrite_group_fromfile(fortchartoc(filnam));
}

void fgmvwrite_group_fromfile__(char filnam[])
{
  gmvwrite_group_fromfile(fortchartoc(filnam));
}


/* ---------------------------------------------------------------- */

void fgmvwrite_group_header()
{
  gmvwrite_group_header();
}

void FGMVWRITE_GROUP_HEADER()
{
  gmvwrite_group_header();
}

void fgmvwrite_group_header_()
{
  gmvwrite_group_header();
}

void fgmvwrite_group_header__()
{
  gmvwrite_group_header();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_group_data(char groupname[], int data_type, int numgrp, 
			  char *group_data)
{
  char name[10];
  /*  terminate groupname with null character  */
  strncpy(name, groupname, 8);
  *(name+8) = '\0';
  gmvwrite_group_data(name, data_type, numgrp, group_data);
}

void FGMVWRITE_GROUP_DATA(char groupname[], int data_type, int numgrp, 
		       char *group_data)
{
  char name[10];
  /*  terminate groupname with null character  */
  strncpy(name, groupname, 8);
  *(name+8) = '\0';
  gmvwrite_group_data(name, data_type, numgrp, group_data);
}

void fgmvwrite_group_data_(char groupname[], int data_type, int numgrp, 
			  char *group_data)
{
  char name[10];
  /*  terminate groupname with null character  */
  strncpy(name, groupname, 8);
  *(name+8) = '\0';
  gmvwrite_group_data(name, data_type, numgrp, group_data);
}

void fgmvwrite_group_data__(char groupname[], int data_type, int numgrp, 
			  char *group_data)
{
  char name[10];
  /*  terminate groupname with null character  */
  strncpy(name, groupname, 8);
  *(name+8) = '\0';
  gmvwrite_group_data(name, data_type, numgrp, group_data);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_group_endgroup()
{
  gmvwrite_group_endgroup();
}

void FGMVWRITE_GROUP_ENDGROUP()
{
  gmvwrite_group_endgroup();
}

void fgmvwrite_group_endgroup_()
{
  gmvwrite_group_endgroup();
}

void fgmvwrite_group_endgroup__()
{
  gmvwrite_group_endgroup();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_surfids(char *surfids)
{
  gmvwrite_surfids(surfids);
} 

void FGMVWRITE_SURFIDS(char *surfids)
{
  gmvwrite_surfids(surfids);
} 

void fgmvwrite_surfids_(char *surfids)
{
  gmvwrite_surfids(surfids);
}     

void fgmvwrite_surfids__(char *surfids)
{
  gmvwrite_surfids(surfids);
} 


/* ---------------------------------------------------------------- */

void fgmvwrite_codename(char codename[])
{
  char name[10];
  /*  terminate codename with null character  */
  strncpy(name, codename, 8);
  *(name+8) = '\0';
  gmvwrite_codename(name);
}

void FGMVWRITE_CODENAME(char codename[])
{
  char name[10];
  /*  terminate codename with null character  */
  strncpy(name, codename, 8);
  *(name+8) = '\0';
  gmvwrite_codename(name);
}

void fgmvwrite_codename_(char codename[])
{
  char name[10];
  /*  terminate codename with null character  */
  strncpy(name, codename, 8);
  *(name+8) = '\0';
  gmvwrite_codename(name);
}

void fgmvwrite_codename__(char codename[])
{
  char name[10];
  /*  terminate codename with null character  */
  strncpy(name, codename, 8);
  *(name+8) = '\0';
  gmvwrite_codename(name);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_codever(char codever[])
{
  char name[10];
  /*  terminate codever with null character  */
  strncpy(name, codever, 8);
  *(name+8) = '\0';
  gmvwrite_codever(name);
}

void FGMVWRITE_CODEVER(char codever[])
{
  char name[10];
  /*  terminate codever with null character  */
  strncpy(name, codever, 8);
  *(name+8) = '\0';
  gmvwrite_codever(name);
}

void fgmvwrite_codever_(char codever[])
{
  char name[10];
  /*  terminate codever with null character  */
  strncpy(name, codever, 8);
  *(name+8) = '\0';
  gmvwrite_codever(name);
}

void fgmvwrite_codever__(char codever[])
{
  char name[10];
  /*  terminate codever with null character  */
  strncpy(name, codever, 8);
  *(name+8) = '\0';
  gmvwrite_codever(name);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_simdate(char simdate[])
{
  char name[10];
  /*  terminate simdate with null character  */
  strncpy(name, simdate, 8);
  *(name+8) = '\0';
  gmvwrite_simdate(name);
}

void FGMVWRITE_SIMDATE(char simdate[])
{
  char name[10];
  /*  terminate simdate with null character  */
  strncpy(name, simdate, 8);
  *(name+8) = '\0';
  gmvwrite_simdate(name);
}

void fgmvwrite_simdate_(char simdate[])
{
  char name[10];
  /*  terminate simdate with null character  */
  strncpy(name, simdate, 8);
  *(name+8) = '\0';
  gmvwrite_simdate(name);
}

void fgmvwrite_simdate__(char simdate[])
{
  char name[10];
  /*  terminate simdate with null character  */
  strncpy(name, simdate, 8);
  *(name+8) = '\0';
  gmvwrite_simdate(name);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_subvars_header()
{
  gmvwrite_subvars_header();
}

void FGMVWRITE_SUBVARS_HEADER()
{
  gmvwrite_subvars_header();
}

void fgmvwrite_subvars_header_()
{
  gmvwrite_subvars_header();
}

void fgmvwrite_subvars_header__()
{
  gmvwrite_subvars_header();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_subvars_name_data(int data_type, int numelem, char varname[], 
				 void *vids, void *vdata)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_subvars_name_data(data_type, numelem, name, vids, vdata);
}

void FGMVWRITE_SUBVARS_NAME_DATA(int data_type, int numelem, char varname[], 
				 void *vids, void *vdata)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_subvars_name_data(data_type, numelem, name, vids, vdata);
}

void fgmvwrite_subvars_name_data_(int data_type, int numelem, char varname[], 
				 void *vids, void *vdata)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_subvars_name_data(data_type, numelem, name, vids, vdata);
}

void fgmvwrite_subvars_name_data__(int data_type, int numelem, char varname[], 
				 void *vids, void *vdata)
{
  char name[10];
  /*  terminate varname with null character  */
  strncpy(name, varname, 8);
  *(name+8) = '\0';
  gmvwrite_subvars_name_data(data_type, numelem, name, vids, vdata);
}


/* ---------------------------------------------------------------- */

void fgmvwrite_subvars_endsubv()
{
  gmvwrite_subvars_endsubv();
}

void FGMVWRITE_SUBVARS_ENDSUBV()
{
  gmvwrite_subvars_endsubv();
}

void fgmvwrite_subvars_endsubv_()
{
  gmvwrite_subvars_endsubv();
}

void fgmvwrite_subvars_endsubv__()
{
  gmvwrite_subvars_endsubv();
}


/* ---------------------------------------------------------------- */

void fgmvwrite_ghosts(int data_type, int numghst, char *ghost_data)
{
  gmvwrite_ghosts(data_type, numghst, ghost_data);
}

void FGMVWRITE_GHOSTS(int data_type, int numghst, char *ghost_data)
{
  gmvwrite_ghosts(data_type, numghst, ghost_data);
}

void fgmvwrite_ghosts_(int data_type, int numghst, char *ghost_data)
{
  gmvwrite_ghosts(data_type, numghst, ghost_data);
}

void fgmvwrite_ghosts__(int data_type, int numghst, char *ghost_data)
{
  gmvwrite_ghosts(data_type, numghst, ghost_data);
}


/* ---------------------------------------------------------------- */
