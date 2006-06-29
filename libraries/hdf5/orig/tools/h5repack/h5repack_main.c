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

#include <stdlib.h>
#include "h5repack.h"

static void usage(void);

/*-------------------------------------------------------------------------
 * Function: main
 *
 * Purpose: h5repack main program
 *
 * Return: 1, error, 0, no error
 *
 * Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
 *
 * Date: May 9, 2003
 *
 * Comments:
 *
 * Modifications: July 2004
 *  Introduced the extra EC or NN option for SZIP
 *
 *-------------------------------------------------------------------------
 */


int main(int argc, char **argv)
{
 char          *infile  = NULL;
 char          *outfile = NULL;
 pack_opt_t    options;            /*the global options */
 int           i, ret;

 /* initialize options  */
 h5repack_init (&options,0);

 for ( i = 1; i < argc; i++)
 {
  if (strcmp(argv[i], "-h") == 0) {
    usage();
    exit(0);
  }
  if (strcmp(argv[i], "-i") == 0) {
   infile = argv[++i];
  }
  else if (strcmp(argv[i], "-o") == 0) {
   outfile = argv[++i];
  }
  else if (strcmp(argv[i], "-v") == 0) {
   options.verbose = 1;
  }
  else if (strcmp(argv[i], "-f") == 0) {

   /* add the -f filter option */
   if (h5repack_addfilter(argv[i+1],&options)<0)
    exit(1);

   /* jump to next */
   ++i;
  }
  else if (strcmp(argv[i], "-l") == 0) {

   /* parse the -l layout option */
   if (h5repack_addlayout(argv[i+1],&options)<0)
    exit(1);

   /* jump to next */
   ++i;
  }

  else if (strcmp(argv[i], "-m") == 0) {
   options.threshold = parse_number(argv[i+1]);
   if ((int)options.threshold==-1) {
    printf("Error: Invalid treshold size <%s>\n",argv[i+1]);
    exit(1);
   }
   ++i;
  }

  else if (strcmp(argv[i], "-e") == 0) {
   read_info(argv[++i],&options);
  }

  else if (argv[i][0] == '-') {
   usage();
   exit(1);
  }
 }

 if (infile == NULL || outfile == NULL)
 {
  usage();
  exit(1);
 }

 /* pack it */
 ret=h5repack(infile,outfile,&options);

 /* free tables */
 h5repack_end(&options);

 if (ret==-1)
  return 1;
 else
  return 0;
}

/*-------------------------------------------------------------------------
 * Function: usage
 *
 * Purpose: print usage
 *
 * Return: void
 *
 *-------------------------------------------------------------------------
 */

static
void usage(void)
{
 printf("h5repack -i input -o output [-h] [-v] [-f 'filter'] [-l 'layout'][-m number][-e file] \n");
 printf("\n");
 printf("-i input          Input HDF5 File\n");
 printf("-o output         Output HDF5 File\n");
 printf("[-h]              Print usage message\n");
 printf("[-v]              Verbose mode. Print output (list of objects, filters, warnings)\n");
 printf("[-f 'filter']     Filter type: 'filter' is a string with the format\n");
 printf("     <list of objects> : <name of filter> = <filter parameters>\n");
 printf("     <list of objects> is a comma separated list of object names\n");
 printf("       meaning apply compression only to those objects.\n");
 printf("       if no object names are specified, the filter is applied to all objects\n");
 printf("     <name of filter> can be:\n");
 printf("       GZIP, to apply the HDF5 GZIP filter (GZIP compression)\n");
 printf("       SZIP, to apply the HDF5 SZIP filter (SZIP compression)\n");
 printf("       SHUF, to apply the HDF5 shuffle filter\n");
 printf("       FLET, to apply the HDF5 checksum filter\n");
 printf("       NONE, to remove the filter\n");
 printf("     <filter parameters> is optional compression info\n");
 printf("       SHUF (no parameter)\n");
 printf("       FLET (no parameter)\n");
 printf("       GZIP=<deflation level> from 1-9\n");
 printf("       SZIP=<pixels per block,coding> (pixels per block is a even number in 2-32 and coding method is 'EC' or 'NN')\n");
 printf("[-l 'layout']     Layout type. 'layout' is a string with the format\n");
 printf("     <list of objects> : <layout type>\n");
 printf("     <list of objects> is a comma separated list of object names,\n");
 printf("       meaning that layout information is supplied for those objects.\n");
 printf("       if no object names are specified, the layout is applied to all objects\n");
 printf("     <layout type> can be:\n");
 printf("       CHUNK, to apply chunking layout\n");
 printf("       COMPA, to apply compact layout\n");
 printf("       CONTI, to apply continuous layout\n");
 printf("     <layout parameters> is present for the chunk case only\n");
 printf("       it is the chunk size of each dimension:\n");
 printf("       <dim_1 x dim_2 x ... dim_n>\n");
 printf("\n");
 printf("-e file           File with the -f and -l options (only filter and layout flags)\n");
 printf("-m number         Do not apply the filter to objects which size in bytes is smaller than number.\n");
 printf("                  If no size is specified a minimum of 1024 bytes is assumed.\n");
 printf("\n");
 printf("Examples of use:\n");
 printf("\n");
 printf("1) h5repack -i file1 -o file2 -f GZIP=1 -v\n");
 printf("\n");
 printf("   Applies GZIP compression to all objects in file1 and saves the output in file2\n");
 printf("\n");
 printf("2) h5repack -i file1 -o file2 -f dset1:SZIP=8,NN -v\n");
 printf("\n");
 printf("   Applies SZIP compression only to object 'dset1'\n");
 printf("\n");
 printf("3) h5repack -i file1 -o file2 -l dset1,dset2:CHUNK=20x10 -v\n");
 printf("\n");
 printf("   Applies chunked layout to objects 'dset1' and 'dset2'\n");
 printf("\n");
}


