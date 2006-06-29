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

#include "hdf5.h"
#include "h5test.h"
#include "h5repack.h"


/*-------------------------------------------------------------------------
 * Function: aux_find_obj
 *
 * Purpose: find the object name NAME (got from the traverse list)
 *  in the repack options list
 *
 *-------------------------------------------------------------------------
 */
static
int aux_find_obj(const char* name,          /* object name from traverse list */
                 pack_opt_t *options,       /* repack options */
																	pack_info_t *obj /*OUT*/)  /* info about object to filter */
{
 char *pdest;
 int  result;
	int  i;

 for ( i=0; i<options->op_tbl->nelems; i++)
 {
  if (strcmp(options->op_tbl->objs[i].path,name)==0)
		{
			*obj =  options->op_tbl->objs[i];
			return i;
		}

		pdest  = strstr(name,options->op_tbl->objs[i].path);
		result = (int)(pdest - name);

		/* found at position 1, meaning without '/' */
		if( pdest != NULL && result==1 )
		{
			*obj =  options->op_tbl->objs[i];
			return i;
		}
 }/*i*/

 return -1;
}


/*-------------------------------------------------------------------------
 * Function: aux_assign_obj
 *
 * Purpose: find the object name NAME (got from the traverse list)
 *  in the repack options list; assign the filter information OBJ
 *
 * Return: 0 not found, 1 found
 *
 *-------------------------------------------------------------------------
 */
static
int aux_assign_obj(const char* name,            /* object name from traverse list */
                   pack_opt_t *options,         /* repack options */
                   pack_info_t *obj /*OUT*/)    /* info about object to filter */
{

 int  idx, i;
	pack_info_t tmp;

	init_packobject(&tmp);

	idx = aux_find_obj(name,options,&tmp);

	/* name was on input */
	if (idx>=0)
	{

		/* applying to all objects */
  if (options->all_layout)
		{
   /* assign the global layout info to the OBJ info */
   tmp.layout=options->layout_g;
			switch (options->layout_g)
			{
			case H5D_CHUNKED:
				tmp.chunk.rank=options->chunk_g.rank;
				for ( i=0; i<tmp.chunk.rank; i++)
					tmp.chunk.chunk_lengths[i]=options->chunk_g.chunk_lengths[i];
				break;
			default:
				break;
			}/*switch*/
		}
		else
		{
			tmp.layout = options->op_tbl->objs[idx].layout;
			switch (tmp.layout)
			{
			case H5D_CHUNKED:
				tmp.chunk.rank = options->op_tbl->objs[idx].chunk.rank;
				for ( i=0; i<tmp.chunk.rank; i++)
					tmp.chunk.chunk_lengths[i]=options->op_tbl->objs[idx].chunk.chunk_lengths[i];
				break;
			default:
				break;
			}/*switch*/

		}

		/* applying to all objects */
		if (options->all_filter)
		{
			/* assign the global filter */
			tmp.nfilters=1;
			tmp.filter[0]=options->filter_g;
		} /* if all */
		else
		{
   tmp.nfilters=options->op_tbl->objs[idx].nfilters;
			for ( i=0; i<tmp.nfilters; i++)
			{
    tmp.filter[i] = options->op_tbl->objs[idx].filter[i];
			}
		}


	} /* if idx */


 /* no input name */

	else
	{

		if (options->all_filter)
		{
   /* assign the global filter */
			tmp.nfilters=1;
			tmp.filter[0]=options->filter_g;
		}


		if (options->all_layout)
		{
			/* assign the global layout info to the OBJ info */
			tmp.layout=options->layout_g;
			switch (options->layout_g)
			{
			case H5D_CHUNKED:
				tmp.chunk.rank=options->chunk_g.rank;
				for ( i=0; i<tmp.chunk.rank; i++)
					tmp.chunk.chunk_lengths[i]=options->chunk_g.chunk_lengths[i];
				break;
			default:
				break;
			}/*switch*/

		}


	}

 *obj = tmp;
 return 1;

}



/*-------------------------------------------------------------------------
 * Function: apply_filters
 *
 * Purpose: apply the filters in the object to the property list;
 *  do extra checking in the case of SZIP; delete all filters in the case
 *  of H5Z_FILTER_NONE present in the PACK_INFO_T filter array
 *
 * Return: 0 success, -1 an error occured
 *
 * Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
 *
 * Date: December 19, 2003
 *
 *-------------------------------------------------------------------------
 */
int apply_filters(const char* name,    /* object name from traverse list */
                  int rank,            /* rank of dataset */
                  hsize_t *dims,       /* dimensions of dataset */
                  hid_t dcpl_id,       /* dataset creation property list */
                  hid_t type_id,       /* dataset datatype */
                  pack_opt_t *options) /* repack options */
{
	int          nfilters;       /* number of filters in DCPL */
 unsigned     aggression;     /* the deflate level */
 hsize_t      nelmts;         /* number of elements in dataset */
 size_t       size;           /* size of datatype in bytes */
 hsize_t      chsize[64];     /* chunk size in elements */
 H5D_layout_t layout;
 int          i;
	pack_info_t  obj;

 if (rank==0)
  goto out;

		/*-------------------------------------------------------------------------
		* initialize the assigment object
		*-------------------------------------------------------------------------
 */
	init_packobject(&obj);


	/*-------------------------------------------------------------------------
 * find options
 *-------------------------------------------------------------------------
 */
 if (aux_assign_obj(name,options,&obj)==0)
  return 0;


 /* check for datasets too small */
	if ((size=H5Tget_size(type_id))==0)
		return 0;
	nelmts=1;
	for (i=0; i<rank; i++)
		nelmts*=dims[i];
	if (nelmts*size < options->threshold )
	{
		if (options->verbose)
			printf("Warning: Filter not applied to <%s>. Dataset smaller than <%d> bytes\n",
			name,(int)options->threshold);
		return 0;
	}

 /* get information about input filters */
 if ((nfilters = H5Pget_nfilters(dcpl_id))<0)
  return -1;
		/*-------------------------------------------------------------------------
		* check if we have filters in the pipeline
		* we want to replace them with the input filters
		* only remove if we are inserting new ones
		*-------------------------------------------------------------------------
 */
 if (nfilters && obj.nfilters ) {
  if (H5Premove_filter(dcpl_id,H5Z_FILTER_ALL)<0)
   return -1;
 }


 /*-------------------------------------------------------------------------
		* check if there is an existent chunk
  * read it only if there is not a requested layout
		*-------------------------------------------------------------------------
  */
 if (obj.layout == -1 )
 {
  if ((layout = H5Pget_layout(dcpl_id))<0)
   return -1;

  if (layout==H5D_CHUNKED)
  {
   if ((rank = H5Pget_chunk(dcpl_id,NELMTS(chsize),chsize/*out*/))<0)
    return -1;
   obj.layout=H5D_CHUNKED;
   obj.chunk.rank=rank;
   for ( i=0; i<rank; i++)
    obj.chunk.chunk_lengths[i] = chsize[i];
  }
 }


	/*-------------------------------------------------------------------------
 * the type of filter and additional parameter
 * type can be one of the filters
 * H5Z_FILTER_NONE       0,  uncompress if compressed
 * H5Z_FILTER_DEFLATE	   1 , deflation like gzip
 * H5Z_FILTER_SHUFFLE    2 , shuffle the data
 * H5Z_FILTER_FLETCHER32 3 , fletcher32 checksum of EDC
 * H5Z_FILTER_SZIP       4 , szip compression
 *-------------------------------------------------------------------------
 */

	if (obj.nfilters)
	{

	/*-------------------------------------------------------------------------
 * filters require CHUNK layout; if we do not have one define a default
 *-------------------------------------------------------------------------
		*/
		if (obj.layout==-1)
		{
			obj.chunk.rank=rank;
			for (i=0; i<rank; i++)
				obj.chunk.chunk_lengths[i] = dims[i];
		}

		for ( i=0; i<obj.nfilters; i++)
		{
			switch (obj.filter[i].filtn)
			{
			default:
				break;

				/*-------------------------------------------------------------------------
				* H5Z_FILTER_DEFLATE	   1 , deflation like gzip
				*-------------------------------------------------------------------------
				*/
			case H5Z_FILTER_DEFLATE:
				aggression=obj.filter[i].cd_values[0];
				/* set up for deflated data */
				if(H5Pset_chunk(dcpl_id, obj.chunk.rank, obj.chunk.chunk_lengths)<0)
					return -1;
				if(H5Pset_deflate(dcpl_id,aggression)<0)
					return -1;
				break;

				/*-------------------------------------------------------------------------
				* H5Z_FILTER_SZIP       4 , szip compression
				*-------------------------------------------------------------------------
				*/
			case H5Z_FILTER_SZIP:
				{
					unsigned  options_mask;
					unsigned  pixels_per_block;

					pixels_per_block=obj.filter[i].cd_values[0];
					if (obj.filter[i].szip_coding==0)
						options_mask=H5_SZIP_NN_OPTION_MASK;
					else
						options_mask=H5_SZIP_EC_OPTION_MASK;

					/* set up for szip data */
					if(H5Pset_chunk(dcpl_id,obj.chunk.rank,obj.chunk.chunk_lengths)<0)
						return -1;
					if (H5Pset_szip(dcpl_id,options_mask,pixels_per_block)<0)
						return -1;

				}
				break;

				/*-------------------------------------------------------------------------
				* H5Z_FILTER_SHUFFLE    2 , shuffle the data
				*-------------------------------------------------------------------------
				*/
			case H5Z_FILTER_SHUFFLE:
				if(H5Pset_chunk(dcpl_id, obj.chunk.rank, obj.chunk.chunk_lengths)<0)
					return -1;
				if (H5Pset_shuffle(dcpl_id)<0)
					return -1;
				break;

				/*-------------------------------------------------------------------------
				* H5Z_FILTER_FLETCHER32 3 , fletcher32 checksum of EDC
				*-------------------------------------------------------------------------
				*/
			case H5Z_FILTER_FLETCHER32:
				if(H5Pset_chunk(dcpl_id, obj.chunk.rank, obj.chunk.chunk_lengths)<0)
					return -1;
				if (H5Pset_fletcher32(dcpl_id)<0)
					return -1;
				break;
			} /* switch */
		}/*i*/

	}
	/*obj.nfilters*/

	/*-------------------------------------------------------------------------
 * layout
 *-------------------------------------------------------------------------
 */

	if (obj.layout>=0)
	{
		/* a layout was defined */
		if (H5Pset_layout(dcpl_id, obj.layout)<0)
			return -1;

		if (H5D_CHUNKED==obj.layout) { /* set up chunk */
			if(H5Pset_chunk(dcpl_id, obj.chunk.rank, obj.chunk.chunk_lengths)<0)
				return -1;
		}
		else if (H5D_COMPACT==obj.layout) {
			if (H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_EARLY)<0)
				return -1;
		}

	}

 return 0;

out:
 if (options->verbose)
  printf("Warning: Filter could not be applied to <%s>\n",name);
 return 0;
}



/*-------------------------------------------------------------------------
 * Function: print_filters
 *
 * Purpose: print the filters in DCPL
 *
 * Return: 0, ok, -1 no
 *
 *-------------------------------------------------------------------------
 */

int print_filters(hid_t dcpl_id)
{
 int          nfilters;       /* number of filters */
 unsigned     filt_flags;     /* filter flags */
 H5Z_filter_t filtn;          /* filter identification number */
 unsigned     cd_values[20];  /* filter client data values */
 size_t       cd_nelmts;      /* filter client number of values */
 size_t       cd_num;         /* filter client data counter */
 char         f_name[256];    /* filter name */
 char         s[64];          /* temporary string buffer */
 int          i;

 /* get information about filters */
 if ((nfilters = H5Pget_nfilters(dcpl_id))<0)
  return -1;

 for (i=0; i<nfilters; i++)
 {
  cd_nelmts = NELMTS(cd_values);
  filtn = H5Pget_filter(dcpl_id,
   (unsigned)i,
   &filt_flags,
   &cd_nelmts,
   cd_values,
   sizeof(f_name),
   f_name);

  f_name[sizeof(f_name)-1] = '\0';
  sprintf(s, "Filter-%d:", i);
  printf("    %-10s %s-%u %s {", s,
   f_name[0]?f_name:"method",
   (unsigned)filtn,
   filt_flags & H5Z_FLAG_OPTIONAL?"OPT":"");
  for (cd_num=0; cd_num<cd_nelmts; cd_num++) {
   printf("%s%u", cd_num?", ":"", cd_values[cd_num]);
  }
  printf("}\n");
 }

 return 0;
}



