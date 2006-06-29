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

/*
 * Programmer:  Robb Matzke <matzke@llnl.gov>
 *              Thursday, July 23, 1998
 *
 * Purpose:	A library for displaying the values of a dataset in a human
 *		readable format.
 */

#include <stdio.h>
#include <stdlib.h>

#include "h5tools.h"
#include "h5tools_ref.h"
#include "h5tools_str.h"
#include "h5tools_utils.h"
#include "H5private.h"

/*
 * The output functions need a temporary buffer to hold a piece of the
 * dataset while it's being printed. This constant sets the limit on the
 * size of that temporary buffer in bytes. For efficiency's sake, choose the
 * largest value suitable for your machine (for testing use a small value).
 */
#if 1
#define H5TOOLS_BUFSIZE         (1024 * 1024)
#else
#define H5TOOLS_BUFSIZE         (1024)
#endif  /* 1 */

#define ALIGN(A,Z)		((((A) + (Z) - 1) / (Z)) * (Z))

/* global variables */
int         indent;
int         compound_data;
FILE       *rawdatastream;	/* should initialize to stdout but gcc moans about it */

/* module-scoped variables */
static int  h5tools_init_g;     /* if h5tools lib has been initialized */
#ifdef H5_HAVE_PARALLEL
static int  h5tools_mpi_init_g; /* if MPI_Init() has been called */
#endif /* H5_HAVE_PARALLEL */

/* Names of VFDs */
static const char *drivernames[]={
    "sec2",
    "family",
    "split",
    "multi",
#ifdef H5_HAVE_STREAM
    "stream",
#endif	/* H5_HAVE_STREAM */
#ifdef H5_HAVE_PARALLEL
    "mpio",
    "mpiposix"
#endif /* H5_HAVE_PARALLEL */
};

/* This enum should match the entries in the above drivers_list since they
 * are indexes into the drivers_list array. */
enum {
    SEC2_IDX = 0
   ,FAMILY_IDX
   ,SPLIT_IDX
   ,MULTI_IDX
#ifdef H5_HAVE_STREAM
   ,STREAM_IDX
#endif	/* H5_HAVE_STREAM */
#ifdef H5_HAVE_PARALLEL
   ,MPIO_IDX
   ,MPIPOSIX_IDX
#endif /* H5_HAVE_PARALLEL */
} driver_idx;
#define NUM_DRIVERS     (sizeof(drivernames) / sizeof(drivernames[0]))

/*-------------------------------------------------------------------------
 * Audience:    Public
 * Chapter:     H5Tools Library
 * Purpose:     Initialize the H5 Tools library
 * Description:
 *      This should be called before any other h5tools function is called.
 *      Effect of any h5tools function called before this has been called is
 *      undetermined.
 * Return:
 *      None
 * Programmer:
 *      Albert Cheng, 2000-10-31
 * Modifications:
 *-------------------------------------------------------------------------
 */
void
h5tools_init(void)
{
    if (!h5tools_init_g) {
	if (!rawdatastream)
	    rawdatastream = stdout;

	h5tools_init_g++;
    }
}

/*-------------------------------------------------------------------------
 * Audience:    Public
 * Chapter:     H5Tools Library
 * Purpose:	Close the H5 Tools library
 * Description:
 *      Close or release resources such as files opened by the library. This
 *      should be called after all other h5tools functions have been called.
 *      Effect of any h5tools function called after this has been called is
 *      undetermined.
 * Return:
 *      None
 * Programmer:
 *      Albert Cheng, 2000-10-31
 * Modifications:
 *-------------------------------------------------------------------------
 */
void
h5tools_close(void)
{
    if (h5tools_init_g) {
	if (rawdatastream && rawdatastream != stdout) {
	    if (fclose(rawdatastream))
		perror("closing rawdatastream");
	    else
		rawdatastream = NULL;
	}

        /* Clean up the reference path table, if it's been used */
        term_ref_path_table();

        /* Shut down the library */
        H5close();

#ifdef H5_HAVE_PARALLEL
        /* Check if we need to shut down MPI */
        if(h5tools_mpi_init_g) {
            MPI_Finalize();
            h5tools_mpi_init_g=0;
        } /* end if */
#endif

	h5tools_init_g = 0;
    }
}

/*-------------------------------------------------------------------------
 * Audience:    Private
 * Chapter:     H5Tools Library
 * Purpose:	Get a FAPL for a driver
 * Description:
 *      Get a FAPL for a given VFL driver name.
 * Return:
 *      None
 * Programmer:
 *      Quincey Koziol, 2004-02-04
 * Modifications:
 *-------------------------------------------------------------------------
 */
static hid_t
h5tools_get_fapl(const char *driver, unsigned *drivernum,
#ifdef H5_HAVE_PARALLEL
int argc, const char *argv[]
#else /* H5_HAVE_PARALLEL */
int UNUSED argc, const char UNUSED *argv[]
#endif /* H5_HAVE_PARALLEL */
)
{
    hid_t               fapl = H5P_DEFAULT;

    /* Determine which driver the user wants to open the file with. Try
     * that driver. If it can't open it, then fail. */
    if (!strcmp(driver, drivernames[SEC2_IDX])) {
        if(drivernum)
            *drivernum = SEC2_IDX;
    } else if (!strcmp(driver, drivernames[FAMILY_IDX])) {
        /* FAMILY Driver */
        if((fapl = H5Pcreate(H5P_FILE_ACCESS))>=0) {
            /* Set member size to be 0 to indicate the current first member size
             * is the member size.
             */
            H5Pset_fapl_family(fapl, (hsize_t)0, H5P_DEFAULT);

            if(drivernum)
                *drivernum = FAMILY_IDX;
        } /* end if */
    } else if (!strcmp(driver, drivernames[SPLIT_IDX])) {
        /* SPLIT Driver */
        if((fapl = H5Pcreate(H5P_FILE_ACCESS))>=0) {
            H5Pset_fapl_split(fapl, "-m.h5", H5P_DEFAULT, "-r.h5", H5P_DEFAULT);

            if(drivernum)
                *drivernum = SPLIT_IDX;
        } /* end if */
    } else if (!strcmp(driver, drivernames[MULTI_IDX])) {
        /* MULTI Driver */
        if((fapl = H5Pcreate(H5P_FILE_ACCESS))>=0) {
            H5Pset_fapl_multi(fapl, NULL, NULL, NULL, NULL, TRUE);

            if(drivernum)
                *drivernum = MULTI_IDX;
        } /* end if */
#ifdef H5_HAVE_STREAM
    } else if (!strcmp(driver, drivernames[STREAM_IDX])) {
        /* STREAM Driver */
        if((fapl = H5Pcreate(H5P_FILE_ACCESS))>=0) {
            H5Pset_fapl_stream(fapl, NULL);

            if(drivernum)
                *drivernum = STREAM_IDX;
        } /* end if */
#endif	/* H5_HAVE_STREAM */
#ifdef H5_HAVE_PARALLEL
    } else if (!strcmp(driver, drivernames[MPIO_IDX])) {
        /* MPI-I/O Driver */
        if((fapl = H5Pcreate(H5P_FILE_ACCESS))>=0) {
            H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);

            /* Initialize the MPI library, if it wasn't already */
            if(!h5tools_mpi_init_g) {
                MPI_Init(&argc, (char ***)&argv);

                h5tools_mpi_init_g=1;
            } /* end if */

            if(drivernum)
                *drivernum = MPIO_IDX;
        } /* end if */
    } else if (!strcmp(driver, drivernames[MPIPOSIX_IDX])) {
        /* MPI-I/O Driver */
        if((fapl = H5Pcreate(H5P_FILE_ACCESS))>=0) {
            H5Pset_fapl_mpiposix(fapl, MPI_COMM_WORLD, TRUE);

            /* Initialize the MPI library, if it wasn't already */
            if(!h5tools_mpi_init_g) {
                MPI_Init(&argc, (char ***)&argv);

                h5tools_mpi_init_g=1;
            } /* end if */

            if(drivernum)
                *drivernum = MPIPOSIX_IDX;
        } /* end if */
#endif /* H5_HAVE_PARALLEL */
    } else {
        fapl=(-1);
    }

    return(fapl);
}

/*-------------------------------------------------------------------------
 * Audience:    Public
 * Chapter:     H5Tools Library
 * Purpose:     Open a file with various VFL drivers.
 * Description:
 *      Loop through the various types of VFL drivers trying to open FNAME.
 *      If the HDF5 library is version 1.2 or less, then we have only the SEC2
 *      driver to try out. If the HDF5 library is greater than version 1.2,
 *      then we have the FAMILY, SPLIT, and MULTI drivers to play with (and
 *      the STREAM driver if H5_HAVE_STREAM is defined, that is).
 *
 *      If DRIVER is non-NULL, then it will try to open the file with that
 *      driver first. We assume that the user knows what they are doing so, if
 *      we fail, then we won't try other file drivers.
 * Return:
 *      On success, returns a file id for the opened file. If DRIVERNAME is
 *      non-null then the first DRIVERNAME_SIZE-1 characters of the driver
 *      name are copied into the DRIVERNAME array and null terminated.
 *
 *      Otherwise, the function returns FAIL. If DRIVERNAME is non-null then
 *      the first byte is set to the null terminator.
 * Programmer:
 *      Lost in the mists of time.
 * Modifications:
 *      Robb Matzke, 2000-06-23
 *      We only have to initialize driver[] on the first call, thereby
 *      preventing memory leaks from repeated calls to H5Pcreate().
 *
 *      Robb Matzke, 2000-06-23
 *      Added DRIVERNAME_SIZE arg to prevent overflows when writing to
 *      DRIVERNAME.
 *
 *      Robb Matzke, 2000-06-23
 *      Added test to prevent coredump when the file could not be opened by
 *      any driver.
 *
 *      Robb Matzke, 2000-06-23
 *      Changed name from H5ToolsFopen() so it jives better with the names we
 *      already have at the top of this source file.
 *
 *      Thomas Radke, 2000-09-12
 *      Added Stream VFD to the driver[] array.
 *
 *      Bill Wendling, 2001-01-10
 *      Changed macro behavior so that if we have a version other than 1.2.x
 *      (i.e., > 1.2), then we do the drivers check.
 *
 *      Bill Wendling, 2001-07-30
 *      Added DRIVER parameter so that the user can specify "try this driver"
 *      instead of the default behaviour. If it fails to open the file with
 *      that driver, this will fail completely (i.e., we won't try the other
 *      drivers). We're assuming the user knows what they're doing. How UNIX
 *      of us.
 *-------------------------------------------------------------------------
 */
hid_t
h5tools_fopen(const char *fname, const char *driver, char *drivername,
              size_t drivername_size, int argc, const char *argv[])
{
    unsigned    drivernum;
    hid_t       fid = FAIL;
    hid_t       fapl = H5P_DEFAULT;

    if (driver && *driver) {
        /* Get the correct FAPL for the given driver */
        if((fapl=h5tools_get_fapl(driver,&drivernum,argc,argv))<0)
            goto done;

        H5E_BEGIN_TRY {
            fid = H5Fopen(fname, H5F_ACC_RDONLY, fapl);
        } H5E_END_TRY;

        if (fid == FAIL)
            goto done;

    } else {
        /* Try to open the file using each of the drivers */
        for (drivernum = 0; drivernum < NUM_DRIVERS; drivernum++) {
            /* Get the correct FAPL for the given driver */
            if((fapl=h5tools_get_fapl(drivernames[drivernum],NULL,argc,argv))<0)
                goto done;

            H5E_BEGIN_TRY {
                fid = H5Fopen(fname, H5F_ACC_RDONLY, fapl);
            } H5E_END_TRY;

            if (fid != FAIL)
                break;
            else {
                /* Close the FAPL */
                H5Pclose(fapl);
                fapl=H5P_DEFAULT;
            } /* end else */
        }
    }

    /* Save the driver name */
    if (drivername && drivername_size) {
        if (fid != FAIL) {
            strncpy(drivername, drivernames[drivernum], drivername_size);
            drivername[drivername_size - 1] = '\0';
        } else {
            /*no file opened*/
            drivername[0] = '\0';
        }
    }

done:
    if(fapl!=H5P_DEFAULT)
        H5Pclose(fapl);
    return fid;
}

/*-------------------------------------------------------------------------
 * Audience:    Public
 * Chapter:     H5Tools Library
 * Purpose:     Count the number of columns in a string.
 * Description:
 *      Count the number of columns in a string. This is the number of
 *      characters in the string not counting line-control characters.
 * Return:
 *      On success, returns the width of the string. Otherwise this function
 *      returns 0.
 * Programmer:
 *       Robb Matzke, Tuesday, April 27, 1999
 * Modifications:
 *-------------------------------------------------------------------------
 */
static size_t
h5tools_ncols(const char *s)
{
    register size_t i;

    for (i = 0; *s; s++)
        if (*s >= ' ')
            i++;

    return i;
}

/*-------------------------------------------------------------------------
 * Audience:    Public
 * Chapter:     H5Tools Library
 * Purpose:     Emit a simple prefix to STREAM.
 * Description:
 *      If /ctx->need_prefix/ is set then terminate the current line (if
 *      applicable), calculate the prefix string, and display it at the start
 *      of a line.
 * Return:
 *      None
 * Programmer:
 *      Robb Matzke, Monday, April 26, 1999
 * Modifications:
 *      Robb Matzke, 1999-09-29
 *	If a new prefix is printed then the current element number is set back
 *	to zero.
 *      pvn, 2004-07-08
 * Added support for printing array indices:
 *  the indentation is printed before the prefix (printed one indentation
 *  level before)
 *-------------------------------------------------------------------------
 */
static void
h5tools_simple_prefix(FILE *stream, const h5dump_t *info,
                      h5tools_context_t *ctx, hsize_t elmtno, int secnum)
{
    h5tools_str_t prefix;
    h5tools_str_t str; /*temporary for indentation */
    size_t templength = 0;
    int i, indentlevel = 0;

    if (!ctx->need_prefix)
	return;

    memset(&prefix, 0, sizeof(h5tools_str_t));
    memset(&str, 0, sizeof(h5tools_str_t));

    /* Terminate previous line, if any */
    if (ctx->cur_column) {
	    fputs(OPT(info->line_suf, ""), stream);
     putc('\n', stream);
	    fputs(OPT(info->line_sep, ""), stream);
    }

    /* Calculate new prefix */
    h5tools_str_prefix(&prefix, info, elmtno, ctx->ndims, ctx->p_min_idx,
                       ctx->p_max_idx, ctx);

    /* Write new prefix to output */
    if (ctx->indent_level >= 0) {
        indentlevel = ctx->indent_level;
    } else {
        /*
         * This is because sometimes we don't print out all the header
         * info for the data (like the tattr-2.ddl example). If that happens
         * the ctx->indent_level is negative so we need to skip the above and
         * just print out the default indent levels.
         */
	indentlevel = ctx->default_indent_level;
    }

    /* when printing array indices, print the indentation before the prefix
       the prefix is printed one indentation level before */
    if (info->pindex) {
     for (i = 0; i < indentlevel-1; i++){
      fputs(h5tools_str_fmt(&str, 0, info->line_indent), stream);
     }
    }

    if (elmtno == 0 && secnum == 0 && info->line_1st)
        fputs(h5tools_str_fmt(&prefix, 0, info->line_1st), stream);
    else if (secnum && info->line_cont)
        fputs(h5tools_str_fmt(&prefix, 0, info->line_cont), stream);
    else
        fputs(h5tools_str_fmt(&prefix, 0, info->line_pre), stream);

    templength = h5tools_str_len(&prefix);

    for (i = 0; i < indentlevel; i++){
        /*we already made the indent for the array indices case */
     if (!info->pindex) {
        fputs(h5tools_str_fmt(&prefix, 0, info->line_indent), stream);
        templength += h5tools_str_len(&prefix);
     }
     else {
      /*we cannot count the prefix for the array indices case */
      templength += h5tools_str_len(&str);
     }
    }

    ctx->cur_column = ctx->prev_prefix_len = templength;
    ctx->cur_elmt = 0;
    ctx->need_prefix = 0;

    /* Free string */
    h5tools_str_close(&prefix);
    h5tools_str_close(&str);
}

/*-------------------------------------------------------------------------
 * Audience:    Public
 * Chapter:     H5Tools Library
 * Purpose:     Prints NELMTS data elements to output STREAM.
 * Description:
 *      Prints some (NELMTS) data elements to output STREAM. The elements are
 *      stored in _MEM as type TYPE and are printed according to the format
 *      described in INFO. The CTX struct contains context information shared
 *      between calls to this function. The FLAGS is a bit field that
 *      indicates whether the data supplied in this call falls at the
 *      beginning or end of the total data to be printed (START_OF_DATA and
 *      END_OF_DATA).
 * Return:
 *      None
 * Programmer:
 *      Robb Matzke, Monday, April 26, 1999
 * Modifications:
 * 	Robb Matzke, 1999-06-04
 *	The `container' argument is the optional dataset for reference types.
 *
 * 	Robb Matzke, 1999-09-29
 *	Understands the `per_line' property which indicates that every Nth
 *	element should begin a new line.
 *
 *      Robb Matzke, LLNL, 2003-06-05
 *      Do not dereference the memory for a variable-length string here.
 *      Deref in h5tools_str_sprint() instead so recursive types are
 *      handled correctly.
	*
	*      Pedro Vicente Nunes, THG, 2005-10-19
 *        pass to the prefix in h5tools_simple_prefix the total position
	*								instead of the current stripmine position i; this is necessary 
	*								to print the array indices 
	*        new field sm_pos in h5tools_context_t, the current stripmine element position 
 *-------------------------------------------------------------------------
 */
void
h5tools_dump_simple_data(FILE *stream, const h5dump_t *info, hid_t container,
                         h5tools_context_t *ctx/*in,out*/, unsigned flags,
                         hsize_t nelmts, hid_t type, void *_mem)
{
    unsigned char	*mem = (unsigned char*)_mem;
    hsize_t		i;		/*element counter		*/
    char		*s, *section;	/*a section of output		*/
    int			secnum;		/*section sequence number	*/
    size_t		size;		/*size of each datum		*/
    size_t		ncols = 80;	/*available output width	*/
    h5tools_str_t	buffer;		/*string into which to render	*/
    int			multiline;	/*datum was multiline		*/
				hsize_t		curr_pos;		/* total data element position 		*/
    int                 elmt_counter = 0;/*counts the # elements printed.
                                          *I (ptl?) needed something that
                                          *isn't going to get reset when a new
                                          *line is formed. I'm going to use
                                          *this var to count elements and
                                          *break after we see a number equal
                                          *to the ctx->size_last_dim.   */

    /* Setup */
    memset(&buffer, 0, sizeof(h5tools_str_t));
    size = H5Tget_size(type);

    if (info->line_ncols > 0)
	ncols = info->line_ncols;

    h5tools_simple_prefix(stream, info, ctx, (hsize_t)0, 0);

    for (i = 0; i < nelmts; i++, ctx->cur_elmt++, elmt_counter++) {
        /* Render the element */
        h5tools_str_reset(&buffer);
        h5tools_str_sprint(&buffer, info, container, type, mem + i * size, ctx);

        if (i + 1 < nelmts || (flags & END_OF_DATA) == 0)
            h5tools_str_append(&buffer, "%s", OPT(info->elmt_suf1, ","));

        s = h5tools_str_fmt(&buffer, 0, "%s");

        /*
         * If the element would split on multiple lines if printed at our
         * current location...
         */
        if (info->line_multi_new == 1 &&
                (ctx->cur_column + h5tools_ncols(s) +
                 strlen(OPT(info->elmt_suf2, " ")) +
                 strlen(OPT(info->line_suf, ""))) > ncols) {
            if (ctx->prev_multiline) {
                /*
                 * ... and the previous element also occupied more than one
                 * line, then start this element at the beginning of a line.
                 */
                ctx->need_prefix = TRUE;
            } else if ((ctx->prev_prefix_len + h5tools_ncols(s) +
                    strlen(OPT(info->elmt_suf2, " ")) +
                    strlen(OPT(info->line_suf, ""))) <= ncols) {
                /*
                 * ...but *could* fit on one line otherwise, then we
                 * should end the current line and start this element on its
                 * own line.
                 */
                ctx->need_prefix = TRUE;
            }
        }

        /*
         * We need to break after each row of a dimension---> we should
         * break at the end of the each last dimension well that is the
         * way the dumper did it before
         */
        if (info->arr_linebreak && ctx->cur_elmt) {
            if (ctx->size_last_dim && (ctx->cur_elmt % ctx->size_last_dim) == 0)
                ctx->need_prefix = TRUE;

            if ((hsize_t)elmt_counter == ctx->size_last_dim) {
                ctx->need_prefix = TRUE;
                elmt_counter = 0;
            }
        }

        /*
         * If the previous element occupied multiple lines and this element
         * is too long to fit on a line then start this element at the
         * beginning of the line.
         */
        if (info->line_multi_new == 1 && ctx->prev_multiline &&
                (ctx->cur_column + h5tools_ncols(s) +
                 strlen(OPT(info->elmt_suf2, " ")) +
                 strlen(OPT(info->line_suf, ""))) > ncols)
            ctx->need_prefix = TRUE;

        /*
         * If too many elements have already been printed then we need to
         * start a new line.
         */
        if (info->line_per_line > 0 && ctx->cur_elmt >= info->line_per_line)
            ctx->need_prefix = TRUE;

        /*
         * Each OPTIONAL_LINE_BREAK embedded in the rendered string can cause
         * the data to split across multiple lines.  We display the sections
         * one-at a time.
         */
        for (secnum = 0, multiline = 0;
                 (section = strtok(secnum ? NULL : s, OPTIONAL_LINE_BREAK));
                 secnum++) {
            /*
             * If the current section plus possible suffix and end-of-line
             * information would cause the output to wrap then we need to
             * start a new line.
             */

            /*
             * Added the info->skip_first because the dumper does not want
             * this check to happen for the first line
             */
            if ((!info->skip_first || i) &&
                    (ctx->cur_column + strlen(section) +
                     strlen(OPT(info->elmt_suf2, " ")) +
                     strlen(OPT(info->line_suf, ""))) > ncols)
                ctx->need_prefix = 1;

            /*
             * Print the prefix or separate the beginning of this element
             * from the previous element.
             */
            if (ctx->need_prefix) {
                if (secnum)
                    multiline++;

																  /* pass to the prefix the total position instead of the current
                   stripmine position i; this is necessary to print the array
                   indices */
                curr_pos = ctx->sm_pos + i;

                h5tools_simple_prefix(stream, info, ctx, curr_pos, secnum);
            } else if ((i || ctx->continuation) && secnum == 0) {
                fputs(OPT(info->elmt_suf2, " "), stream);
                ctx->cur_column += strlen(OPT(info->elmt_suf2, " "));
            }

            /* Print the section */
            fputs(section, stream);
            ctx->cur_column += strlen(section);
        }

        ctx->prev_multiline = multiline;
    }

    h5tools_str_close(&buffer);
}

/*-------------------------------------------------------------------------
 * Audience:    Public
 * Chapter:     H5Tools Library
 * Purpose:     Dump out a subset of a dataset.
 * Description:
 *      Select a hyperslab from the dataset DSET using the parameters
 *      specified in SSET. Dump this out to STREAM.
 * Return:
 *      On success, return SUCCEED. Otherwise, the function returns FAIL.
 * Programmer:
 *      Bill Wendling, Wednesday, 07. March 2001
 * Modifications:
 *-------------------------------------------------------------------------
 */
static herr_t
h5tools_dump_simple_subset(FILE *stream, const h5dump_t *info, hid_t dset,
                           hid_t p_type, struct subset_t *sset,
                           int indentlevel)
{
    herr_t              ret;                    /*the value to return   */
    hid_t		f_space;		/*file data space	*/
    int			i;                      /*counters		*/
    hsize_t		zero = 0;               /*vector of zeros	*/
    unsigned int	flags;			/*buffer extent flags	*/
    hsize_t		total_size[H5S_MAX_RANK];/*total size of dataset*/

    /* Print info */
    h5tools_context_t	ctx;			/*print context		*/
    size_t		p_type_nbytes;		/*size of memory type	*/

    /* Stripmine info */
    hsize_t		sm_size[H5S_MAX_RANK];	/*stripmine size	*/
    hsize_t		sm_nbytes;		/*bytes per stripmine	*/
    hsize_t		sm_nelmts;		/*elements per stripmine*/
    unsigned char      *sm_buf = NULL;		/*buffer for raw data	*/
    hid_t		sm_space;		/*stripmine data space	*/

    hsize_t             count;

    ret = FAIL;     /* be pessimistic */
    f_space = H5Dget_space(dset);

    if (f_space == FAIL)
        goto done;

    /*
     * check that everything looks okay. the dimensionality must not be too
     * great and the dimensionality of the items selected for printing must
     * match the dimensionality of the dataset.
     */
    memset(&ctx, 0, sizeof(ctx));
    ctx.indent_level = indentlevel;
    ctx.need_prefix = 1;
    ctx.ndims = H5Sget_simple_extent_ndims(f_space);

    if ((size_t)ctx.ndims > NELMTS(sm_size))
        goto done_close;

    /* assume entire data space to be printed */
    if (ctx.ndims > 0)
        for (i = 0; i < ctx.ndims; i++)
            ctx.p_min_idx[i] = 0;

    H5Sget_simple_extent_dims(f_space, total_size, NULL);
	assert(total_size[ctx.ndims - 1]==(hsize_t)((int)(total_size[ctx.ndims - 1])));
    ctx.size_last_dim = (int)(total_size[ctx.ndims - 1]);

    count = sset->count[ctx.ndims - 1];
    sset->count[ctx.ndims - 1] = 1;

    if(ctx.ndims>0)
     init_acc_pos(&ctx,total_size);

    for (; count > 0; sset->start[ctx.ndims - 1] += sset->stride[ctx.ndims - 1],
                      count--) {
        /* calculate the potential number of elements we're going to print */
        H5Sselect_hyperslab(f_space, H5S_SELECT_SET,
            sset->start,
            sset->stride,
            sset->count,
            sset->block);
        sm_nelmts = H5Sget_select_npoints(f_space);

        /*
         * start (0, 0)
         * block (2, 2)
         * stride (15, 5)
         * count (4, 3)
         *
         * make:
         *
         * for up to "count" times.
         *
         * start (0, += stride[last_dim])
         * block (2, 2)
         * stride (15, 5)
         * count (4, 1)
         */

        if (sm_nelmts == 0) {
            /* nothing to print */
            ret = SUCCEED;
            goto done_close;
        }

        /*
         * determine the strip mine size and allocate a buffer. the strip mine is
         * a hyperslab whose size is manageable.
         */
        sm_nbytes = p_type_nbytes = H5Tget_size(p_type);

        if (ctx.ndims > 0)
            for (i = ctx.ndims; i > 0; --i) {
                sm_size[i - 1] = MIN(total_size[i - 1], H5TOOLS_BUFSIZE / sm_nbytes);
                sm_nbytes *= sm_size[i - 1];
                assert(sm_nbytes > 0);
            }

        assert(sm_nbytes == (hsize_t)((size_t)sm_nbytes)); /*check for overflow*/
        sm_buf = malloc((size_t)sm_nelmts * p_type_nbytes);
        sm_space = H5Screate_simple(1, &sm_nelmts, NULL);

        H5Sselect_hyperslab(sm_space, H5S_SELECT_SET, &zero, NULL, &sm_nelmts, NULL);

        /* Read the data */
        if (H5Dread(dset, p_type, sm_space, f_space, H5P_DEFAULT, sm_buf) < 0) {
            H5Sclose(f_space);
            H5Sclose(sm_space);
            free(sm_buf);
            return FAIL;
        }

        /* Print the data */
        flags = START_OF_DATA;

        if (count == 1)
            flags |= END_OF_DATA;

        for (i = 0; i < ctx.ndims; i++) {
            ctx.p_max_idx[i] = ctx.p_min_idx[i] + MIN(total_size[i], sm_size[i]);
        }

        h5tools_dump_simple_data(stream, info, dset, &ctx, flags, sm_nelmts,
                                 p_type, sm_buf);
        free(sm_buf);
        ctx.continuation++;
    }

    /* Terminate the output */
    if (ctx.cur_column) {
        fputs(OPT(info->line_suf, ""), stream);
        putc('\n', stream);
        fputs(OPT(info->line_sep, ""), stream);
    }

    ret = SUCCEED;

done_close:
    H5Sclose(f_space);
done:
    return ret;
}

/*-------------------------------------------------------------------------
 * Audience:    Public
 * Chapter:     H5Tools Library
 * Purpose:	Print some values from a dataset with a simple data space.
 * Description:
 *      This is a special case of h5tools_dump_dset(). This function only
 *      intended for dumping datasets -- it does strip mining and some other
 *      things which are unnecessary for smaller objects such as attributes
 *      (to print small objects like attributes simply read the attribute and
 *      call h5tools_dump_simple_mem()).
 * Return:
 *      On success, the function returns SUCCEED. Otherwise, the function
 *      returns FAIL.
 * Programmer:
 *      Robb Matzke, Thursday, July 23, 1998
 * Modifications:
 *-------------------------------------------------------------------------
 */
static int
h5tools_dump_simple_dset(FILE *stream, const h5dump_t *info, hid_t dset,
                         hid_t p_type, int indentlevel)
{
    hid_t		f_space;		/*file data space	*/
    hsize_t		elmtno;			/*counter		*/
    int			i;			/*counter		*/
    int			carry;			/*counter carry value	*/
    hsize_t		zero[8];		/*vector of zeros	*/
    unsigned int	flags;			/*buffer extent flags	*/
    hsize_t		total_size[H5S_MAX_RANK];/*total size of dataset*/

    /* Print info */
    h5tools_context_t	ctx;			/*print context		*/
    size_t		p_type_nbytes;		/*size of memory type	*/
    hsize_t		p_nelmts;		/*total selected elmts	*/

    /* Stripmine info */
    hsize_t		sm_size[H5S_MAX_RANK];	/*stripmine size	*/
    hsize_t		sm_nbytes;		/*bytes per stripmine	*/
    hsize_t		sm_nelmts;		/*elements per stripmine*/
    unsigned char      *sm_buf = NULL;		/*buffer for raw data	*/
    hid_t		sm_space;		/*stripmine data space	*/

    /* Hyperslab info */
    hsize_t		hs_offset[H5S_MAX_RANK];/*starting offset	*/
    hsize_t		hs_size[H5S_MAX_RANK];	/*size this pass	*/
    hsize_t		hs_nelmts;		/*elements in request	*/

    /* VL data special information */
    unsigned int        vl_data = 0;            /*contains VL datatypes */

    f_space = H5Dget_space(dset);

    if (f_space == FAIL)
        return FAIL;

    /*
     * Check that everything looks okay. The dimensionality must not be too
     * great and the dimensionality of the items selected for printing must
     * match the dimensionality of the dataset.
     */
    memset(&ctx, 0, sizeof(ctx));
    ctx.indent_level = indentlevel;
    ctx.need_prefix = 1;
    ctx.ndims = H5Sget_simple_extent_ndims(f_space);

    if ((size_t)ctx.ndims > NELMTS(sm_size)) {
        H5Sclose(f_space);
        return FAIL;
    }

    /* Assume entire data space to be printed */
    if (ctx.ndims > 0)
        for (i = 0; i < ctx.ndims; i++)
            ctx.p_min_idx[i] = 0;

    H5Sget_simple_extent_dims(f_space, total_size, NULL);

    /* calculate the number of elements we're going to print */
    p_nelmts = 1;

    if (ctx.ndims > 0) {
        for (i = 0; i < ctx.ndims; i++)
            p_nelmts *= total_size[i];
        ctx.size_last_dim = (total_size[ctx.ndims - 1]);
    } /* end if */
    else
        ctx.size_last_dim = 0;

    if (p_nelmts == 0) {
        /* nothing to print */
        H5Sclose(f_space);
        return SUCCEED;
    }

    /* Check if we have VL data in the dataset's datatype */
    if (H5Tdetect_class(p_type, H5T_VLEN) == TRUE)
        vl_data = TRUE;

    /*
     * Determine the strip mine size and allocate a buffer. The strip mine is
     * a hyperslab whose size is manageable.
     */
    sm_nbytes = p_type_nbytes = H5Tget_size(p_type);

    if (ctx.ndims > 0) {
        for (i = ctx.ndims; i > 0; --i) {
            sm_size[i - 1] = MIN(total_size[i - 1], H5TOOLS_BUFSIZE / sm_nbytes);
            sm_nbytes *= sm_size[i - 1];
            assert(sm_nbytes > 0);
        }
    }

    assert(sm_nbytes == (hsize_t)((size_t)sm_nbytes)); /*check for overflow*/
    sm_buf = malloc((size_t)sm_nbytes);

    sm_nelmts = sm_nbytes / p_type_nbytes;
    sm_space = H5Screate_simple(1, &sm_nelmts, NULL);

    if(ctx.ndims>0)
     init_acc_pos(&ctx,total_size);

    /* The stripmine loop */
    memset(hs_offset, 0, sizeof hs_offset);
    memset(zero, 0, sizeof zero);

    for (elmtno = 0; elmtno < p_nelmts; elmtno += hs_nelmts) {
        /* Calculate the hyperslab size */
        if (ctx.ndims > 0) {
            for (i = 0, hs_nelmts = 1; i < ctx.ndims; i++) {
                hs_size[i] = MIN(total_size[i] - hs_offset[i], sm_size[i]);
                ctx.p_max_idx[i] = ctx.p_min_idx[i] + hs_size[i];
                hs_nelmts *= hs_size[i];
            }

            H5Sselect_hyperslab(f_space, H5S_SELECT_SET, hs_offset, NULL,
                                hs_size, NULL);
            H5Sselect_hyperslab(sm_space, H5S_SELECT_SET, zero, NULL,
                                &hs_nelmts, NULL);
        } else {
            H5Sselect_all(f_space);
            H5Sselect_all(sm_space);
            hs_nelmts = 1;
        }

        /* Read the data */
        if (H5Dread(dset, p_type, sm_space, f_space, H5P_DEFAULT, sm_buf) < 0) {
            H5Sclose(f_space);
            H5Sclose(sm_space);
            free(sm_buf);
            return FAIL;
        }

        /* Print the data */
        flags = (elmtno == 0) ? START_OF_DATA : 0;
        flags |= ((elmtno + hs_nelmts) >= p_nelmts) ? END_OF_DATA : 0;

         /* initialize the current stripmine position; this is necessary to print the array
           indices */
        ctx.sm_pos = elmtno;

        h5tools_dump_simple_data(stream, info, dset, &ctx, flags, hs_nelmts,
                             p_type, sm_buf);

        /* Reclaim any VL memory, if necessary */
        if(vl_data)
            H5Dvlen_reclaim(p_type, sm_space, H5P_DEFAULT, sm_buf);

        /* Calculate the next hyperslab offset */
        for (i = ctx.ndims, carry = 1; i > 0 && carry; --i) {
            ctx.p_min_idx[i - 1] = ctx.p_max_idx[i - 1];
            hs_offset[i - 1] += hs_size[i - 1];

            if (hs_offset[i - 1] == total_size[i - 1])
                hs_offset[i - 1] = 0;
            else
                carry = 0;
        }

        ctx.continuation++;
    }

    /* Terminate the output */
    if (ctx.cur_column) {
        fputs(OPT(info->line_suf, ""), stream);
        putc('\n', stream);
        fputs(OPT(info->line_sep, ""), stream);
    }

    H5Sclose(sm_space);
    H5Sclose(f_space);

    free(sm_buf);

    return SUCCEED;
}

/*-------------------------------------------------------------------------
 * Function:	h5tools_dump_simple_mem
 *
 * Purpose:	Print some values from memory with a simple data space.
 *		This is a special case of h5tools_dump_mem().
 *
 * Return:	Success:    SUCCEED
 *
 *		Failure:    FAIL
 *
 * Programmer:	Robb Matzke
 *              Thursday, July 23, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static int
h5tools_dump_simple_mem(FILE *stream, const h5dump_t *info, hid_t obj_id,
                        hid_t type, hid_t space, void *mem, int indentlevel)
{
    int			i;			/*counters		*/
    hsize_t		nelmts;			/*total selected elmts	*/
    h5tools_context_t	ctx;			/*printing context	*/

    /*
     * Check that everything looks okay.  The dimensionality must not be too
     * great and the dimensionality of the items selected for printing must
     * match the dimensionality of the dataset.
     */
    memset(&ctx, 0, sizeof(ctx));
    ctx.ndims = H5Sget_simple_extent_ndims(space);

    if ((size_t)ctx.ndims > NELMTS(ctx.p_min_idx))
        return FAIL;

    ctx.indent_level = indentlevel;
    ctx.need_prefix = 1;

    /* Assume entire data space to be printed */
    for (i = 0; i < ctx.ndims; i++)
        ctx.p_min_idx[i] = 0;

    H5Sget_simple_extent_dims(space, ctx.p_max_idx, NULL);

    for (i = 0, nelmts = 1; ctx.ndims != 0 && i < ctx.ndims; i++)
        nelmts *= ctx.p_max_idx[i] - ctx.p_min_idx[i];

    if (nelmts == 0)
        return SUCCEED; /*nothing to print*/
    if(ctx.ndims>0) {
        assert(ctx.p_max_idx[ctx.ndims - 1]==(hsize_t)((int)ctx.p_max_idx[ctx.ndims - 1]));
        ctx.size_last_dim = (int)(ctx.p_max_idx[ctx.ndims - 1]);
    } /* end if */
    else
        ctx.size_last_dim = 0;

    if(ctx.ndims>0)
     init_acc_pos(&ctx,ctx.p_max_idx);

    /* Print it */
    h5tools_dump_simple_data(stream, info, obj_id, &ctx,
                             START_OF_DATA | END_OF_DATA, nelmts, type, mem);

    /* Terminate the output */
    if (ctx.cur_column) {
        fputs(OPT(info->line_suf, ""), stream);
        putc('\n', stream);
        fputs(OPT(info->line_sep, ""), stream);
    }

    return SUCCEED;
}

/*-------------------------------------------------------------------------
 * Function:	h5tools_dump_dset
 *
 * Purpose:	Print some values from a dataset DSET to the file STREAM
 *		after converting all types to P_TYPE (which should be a
 *		native type).  If P_TYPE is a negative value then it will be
 *		computed from the dataset type using only native types.
 *
 * Note:	This function is intended only for datasets since it does
 *		some things like strip mining which are unnecessary for
 *		smaller objects such as attributes. The easiest way to print
 *		small objects is to read the object into memory and call
 *		h5tools_dump_mem().
 *
 * Return:	Success:    SUCCEED
 *
 *		Failure:    FAIL
 *
 * Programmer:	Robb Matzke
 *              Thursday, July 23, 1998
 *
 * Modifications:
 * 		Robb Matzke, 1999-06-07
 *		If info->raw is set then the memory datatype will be the same
 *		as the file datatype.
 *
 *		Bill Wendling, 2001-02-27
 *		Renamed to ``h5tools_dump_dset'' and added the subsetting
 *		parameter.
 *
 *-------------------------------------------------------------------------
 */
int
h5tools_dump_dset(FILE *stream, const h5dump_t *info, hid_t dset, hid_t _p_type,
                  struct subset_t *sset, int indentlevel)
{
    hid_t     f_space;
    hid_t     p_type = _p_type;
    hid_t     f_type;
    H5S_class_t space_type;
    int       status = FAIL;
    h5dump_t  info_dflt;

    /* Use default values */
    if (!stream)
        stream = stdout;

    if (!info) {
        memset(&info_dflt, 0, sizeof info_dflt);
        info = &info_dflt;
    }

    if (p_type < 0) {
        f_type = H5Dget_type(dset);

        if (info->raw)
            p_type = H5Tcopy(f_type);
        else
            p_type = h5tools_get_native_type(f_type);

        H5Tclose(f_type);

        if (p_type < 0)
            goto done;
    }

    /* Check the data space */
    f_space = H5Dget_space(dset);

    space_type = H5Sget_simple_extent_type(f_space);

    /* Print the data */
    if (space_type == H5S_SIMPLE || space_type == H5S_SCALAR) {
        if (!sset)
            status = h5tools_dump_simple_dset(rawdatastream, info, dset, p_type,
                                              indentlevel);
        else
            status = h5tools_dump_simple_subset(rawdatastream, info, dset, p_type,
                                                sset, indentlevel);
    } else /* space is H5S_NULL */
        status = SUCCEED;

    /* Close the dataspace */
    H5Sclose(f_space);

done:
    if (p_type != _p_type)
        H5Tclose(p_type);

    return status;
}

/*-------------------------------------------------------------------------
 * Function:	h5tools_dump_mem
 *
 * Purpose:	Displays the data contained in MEM. MEM must have the
 *		specified data TYPE and SPACE.  Currently only simple data
 *		spaces are allowed and only the `all' selection.
 *
 * Return:	Success:    SUCCEED
 *
 *		Failure:    FAIL
 *
 * Programmer:	Robb Matzke
 *              Wednesday, January 20, 1999
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
int
h5tools_dump_mem(FILE *stream, const h5dump_t *info, hid_t obj_id, hid_t type,
                 hid_t space, void *mem, int indentlevel)
{
    h5dump_t    info_dflt;

    /* Use default values */
    if (!stream)
	stream = stdout;

    if (!info) {
	memset(&info_dflt, 0, sizeof(info_dflt));
	info = &info_dflt;
    }

    /* Check the data space */
    if (H5Sis_simple(space) <= 0)
	return -1;

    return h5tools_dump_simple_mem(stream, info, obj_id, type, space, mem,
                                   indentlevel);
}

/*-------------------------------------------------------------------------
 * Function:    init_acc_pos
 *
 * Purpose:     initialize accumulator and matrix position
 *
 * Return:      void
 *
 * Programmer:  pvn
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void init_acc_pos(h5tools_context_t	*ctx, hsize_t *dims)
{
 int i;

 assert(ctx->ndims);

 ctx->acc[ctx->ndims-1]=1;
 for(i=(ctx->ndims-2); i>=0; i--)
 {
  ctx->acc[i]=ctx->acc[i+1] * dims[i+1];
 }
 for ( i = 0; i < ctx->ndims; i++)
  ctx->pos[i]=0;
}

