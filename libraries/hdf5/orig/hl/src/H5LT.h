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

#ifndef _H5LT_H
#define _H5LT_H

#include <hdf5.h>

#define TESTING(WHAT)	{printf("%-70s", "Testing " WHAT); fflush(stdout);}
#define PASSED()	{puts(" PASSED");fflush(stdout);}
#define H5_FAILED()	{puts("*FAILED*");fflush(stdout);}
#define SKIPPED()	{puts(" -SKIP-");fflush(stdout);}
#define EXAMPLE(WHAT)	{printf("%-70s", "Example " WHAT); fflush(stdout);}


#ifdef __cplusplus
extern "C" {
#endif


/*-------------------------------------------------------------------------
 *
 * Make dataset functions
 *
 *-------------------------------------------------------------------------
 */

herr_t H5LTmake_dataset( hid_t loc_id,
                         const char *dset_name,
                         int rank,
                         const hsize_t *dims,
                         hid_t type_id,
                         const void *buffer );

herr_t H5LTmake_dataset_char( hid_t loc_id,
                              const char *dset_name,
                              int rank,
                              const hsize_t *dims,
                              const char *buffer );

herr_t H5LTmake_dataset_short( hid_t loc_id,
                               const char *dset_name,
                               int rank,
                               const hsize_t *dims,
                               const short *buffer );

herr_t H5LTmake_dataset_int( hid_t loc_id,
                             const char *dset_name,
                             int rank,
                             const hsize_t *dims,
                             const int *buffer );

herr_t H5LTmake_dataset_long( hid_t loc_id,
                              const char *dset_name,
                              int rank,
                              const hsize_t *dims,
                              const long *buffer );

herr_t H5LTmake_dataset_float( hid_t loc_id,
                               const char *dset_name,
                               int rank,
                               const hsize_t *dims,
                               const float *buffer );

herr_t H5LTmake_dataset_double( hid_t loc_id,
                                const char *dset_name,
                                int rank,
                                const hsize_t *dims,
                                const double *buffer );

herr_t H5LTmake_dataset_string(hid_t loc_id,
                               const char *dset_name,
                               const char *buf );


/*-------------------------------------------------------------------------
 *
 * Read dataset functions
 *
 *-------------------------------------------------------------------------
 */

herr_t H5LTread_dataset( hid_t loc_id,
                         const char *dset_name,
                         hid_t type_id,
                         void *buffer );

herr_t H5LTread_dataset_char( hid_t loc_id,
                              const char *dset_name,
                              char *buffer );

herr_t H5LTread_dataset_short( hid_t loc_id,
                               const char *dset_name,
                               short *buffer );

herr_t H5LTread_dataset_int( hid_t loc_id,
                             const char *dset_name,
                             int *buffer );

herr_t H5LTread_dataset_long( hid_t loc_id,
                              const char *dset_name,
                              long *buffer );

herr_t H5LTread_dataset_float( hid_t loc_id,
                               const char *dset_name,
                               float *buffer );

herr_t H5LTread_dataset_double( hid_t loc_id,
                                const char *dset_name,
                                double *buffer );

herr_t H5LTread_dataset_string( hid_t loc_id,
                                const char *dset_name,
                                char *buf );

/*-------------------------------------------------------------------------
 *
 * Query dataset functions
 *
 *-------------------------------------------------------------------------
 */


herr_t H5LTget_dataset_ndims( hid_t loc_id,
                             const char *dset_name,
                             int *rank );

herr_t H5LTget_dataset_info( hid_t loc_id,
                             const char *dset_name,
                             hsize_t *dims,
                             H5T_class_t *type_class,
                             size_t *type_size );

herr_t H5LTfind_dataset( hid_t loc_id, const char *name );



/*-------------------------------------------------------------------------
 *
 * Set attribute functions
 *
 *-------------------------------------------------------------------------
 */


herr_t H5LTset_attribute_string( hid_t loc_id,
                                 const char *obj_name,
                                 const char *attr_name,
                                 const char *attr_data );

herr_t H5LTset_attribute_char( hid_t loc_id,
                               const char *obj_name,
                               const char *attr_name,
                               const char *buffer,
                               size_t size );

herr_t H5LTset_attribute_uchar( hid_t loc_id,
                               const char *obj_name,
                               const char *attr_name,
                               const unsigned char *buffer,
                               size_t size );

herr_t H5LTset_attribute_short( hid_t loc_id,
                              const char *obj_name,
                              const char *attr_name,
                              const short *buffer,
                              size_t size );

herr_t H5LTset_attribute_ushort( hid_t loc_id,
                              const char *obj_name,
                              const char *attr_name,
                              const unsigned short *buffer,
                              size_t size );

herr_t H5LTset_attribute_int( hid_t loc_id,
                              const char *obj_name,
                              const char *attr_name,
                              const int *buffer,
                              size_t size );

herr_t H5LTset_attribute_uint( hid_t loc_id,
                              const char *obj_name,
                              const char *attr_name,
                              const unsigned int *buffer,
                              size_t size );

herr_t H5LTset_attribute_long( hid_t loc_id,
                               const char *obj_name,
                               const char *attr_name,
                               const long *buffer,
                               size_t size );

herr_t H5LTset_attribute_ulong( hid_t loc_id,
                               const char *obj_name,
                               const char *attr_name,
                               const unsigned long *buffer,
                               size_t size );

herr_t H5LTset_attribute_float( hid_t loc_id,
                                const char *obj_name,
                                const char *attr_name,
                                const float *buffer,
                                size_t size );

herr_t H5LTset_attribute_double( hid_t loc_id,
                                 const char *obj_name,
                                 const char *attr_name,
                                 const double *buffer,
                                 size_t size );

/*-------------------------------------------------------------------------
 *
 * Get attribute functions
 *
 *-------------------------------------------------------------------------
 */

herr_t H5LTget_attribute( hid_t loc_id,
                          const char *obj_name,
                          const char *attr_name,
                          hid_t mem_type_id,
                          void *data );

herr_t H5LTget_attribute_string( hid_t loc_id,
                                 const char *obj_name,
                                 const char *attr_name,
                                 char *data );

herr_t H5LTget_attribute_char( hid_t loc_id,
                               const char *obj_name,
                               const char *attr_name,
                               char *data );

herr_t H5LTget_attribute_uchar( hid_t loc_id,
                               const char *obj_name,
                               const char *attr_name,
                               unsigned char *data );

herr_t H5LTget_attribute_short( hid_t loc_id,
                                const char *obj_name,
                                const char *attr_name,
                                short *data );

herr_t H5LTget_attribute_ushort( hid_t loc_id,
                                const char *obj_name,
                                const char *attr_name,
                                unsigned short *data );

herr_t H5LTget_attribute_int( hid_t loc_id,
                              const char *obj_name,
                              const char *attr_name,
                              int *data );

herr_t H5LTget_attribute_uint( hid_t loc_id,
                              const char *obj_name,
                              const char *attr_name,
                              unsigned int *data );

herr_t H5LTget_attribute_long( hid_t loc_id,
                               const char *obj_name,
                               const char *attr_name,
                               long *data );

herr_t H5LTget_attribute_ulong( hid_t loc_id,
                               const char *obj_name,
                               const char *attr_name,
                               unsigned long *data );

herr_t H5LTget_attribute_float( hid_t loc_id,
                                const char *obj_name,
                                const char *attr_name,
                                float *data );

herr_t H5LTget_attribute_double( hid_t loc_id,
                                 const char *obj_name,
                                 const char *attr_name,
                                 double *data );


/*-------------------------------------------------------------------------
 *
 * Query attribute functions
 *
 *-------------------------------------------------------------------------
 */


herr_t H5LTget_attribute_ndims( hid_t loc_id,
                                const char *obj_name,
                                const char *attr_name,
                                int *rank );

herr_t H5LTget_attribute_info( hid_t loc_id,
                               const char *obj_name,
                               const char *attr_name,
                               hsize_t *dims,
                               H5T_class_t *type_class,
                               size_t *type_size );





/*-------------------------------------------------------------------------
 *
 * General functions
 *
 *-------------------------------------------------------------------------
 */


hid_t H5LTcreate_compound_type( hsize_t nfields, size_t size, const char *field_names[],
                                const size_t *field_offset, const hid_t *field_types );


herr_t H5LTrepack( hsize_t nfields,
                   hsize_t nrecords,
                   size_t src_size,
                   const size_t *src_offset,
                   const size_t *src_sizes,
                   size_t dst_size,
                   const size_t *dst_offset,
                   const size_t *dst_sizes,
                   unsigned char *src_buf,
                   unsigned char *dst_buf );



/*-------------------------------------------------------------------------
 *
 * Private functions
 *
 *-------------------------------------------------------------------------
 */


herr_t H5LT_get_attribute_mem( hid_t obj_id,
                           const char *attr_name,
                           hid_t mem_type_id,
                           void *data );

herr_t H5LT_get_attribute_disk( hid_t obj_id,
                           const char *attr_name,
                           void *data );

herr_t H5LT_find_attribute( hid_t loc_id, const char *name );


herr_t H5LT_set_attribute_numerical( hid_t loc_id,
                                     const char *obj_name,
                                     const char *attr_name,
                                     size_t size,
                                     hid_t type_id,
                                     const void *data );




#ifdef __cplusplus
}
#endif

#endif
