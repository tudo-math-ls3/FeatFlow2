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
 * MPIO independent overlapping writes.
 *
 * First n-1 processes open 1 file.
 * Each of the n-1 process writes chunks of data to the file in round-robin
 * fashion, in a interleaved but not overlapped fashion.  Using increasing
 * chunk sizes for the benefits of testing different write sizes and also
 * reducing the numbers of writes.
 *
 * Last process (n-1) just waits.
 * First n-1 processes finish writing and cloose the file.
 * Last process opens the same file and verifies the data.
 */

#include "testphdf5.h"

/* FILENAME and filenames must have the same number of names */
const char *FILENAME[2]={
	    "MPItest",
	    NULL};
char	filenames[2][200];
int	nerrors = 0;
hid_t	fapl;				/* file access property list */

/* protocols */
static int errors_sum(int nerrs);

#define MPIO_TEST_WRITE_SIZE 1024*1024     /* 1 MB */

static int
test_mpio_overlap_writes(char *filename)
{
    int mpi_size, mpi_rank;
    MPI_Comm comm;
    MPI_Info info = MPI_INFO_NULL;
    int color, mrc;
    MPI_File	fh;
    int i;
    int vrfyerrs, nerrs;
    unsigned char  buf[4093];		/* use some prime number for size */
    int bufsize = sizeof(buf);
    MPI_Offset  stride;
    MPI_Offset  mpi_off;
    MPI_Status  mpi_stat;


    if (VERBOSE_MED)
	printf("MPIO independent overlapping writes test on file %s\n",
	    filename);

    nerrs = 0;
    /* set up MPI parameters */
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

    /* Need at least 2 processes */
    if (mpi_size < 2) {
	if (MAINPROCESS)
	    printf("Need at least 2 processes to run MPIO test.\n");
	    printf(" -SKIP- \n");
	return 0;
    }

    /* splits processes 0 to n-2 into one comm. and the last one into another */
    color = ((mpi_rank < (mpi_size - 1)) ? 0 : 1);
    mrc = MPI_Comm_split (MPI_COMM_WORLD, color, mpi_rank, &comm);
    VRFY((mrc==MPI_SUCCESS), "Comm_split succeeded");

    if (color==0){
	/* First n-1 processes (color==0) open a file and write it */
	mrc = MPI_File_open(comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR,
		info, &fh);
	VRFY((mrc==MPI_SUCCESS), "");

	stride = 1;
	mpi_off = mpi_rank*stride;
	while (mpi_off < MPIO_TEST_WRITE_SIZE){
	    /* make sure the write does not exceed the TEST_WRITE_SIZE */
	    if (mpi_off+stride > MPIO_TEST_WRITE_SIZE)
		stride = MPIO_TEST_WRITE_SIZE - mpi_off;

	    /* set data to some trivial pattern for easy verification */
	    for (i=0; i<stride; i++)
		buf[i] = (unsigned char)(mpi_off+i);
	    mrc = MPI_File_write_at(fh, mpi_off, buf, (int)stride, MPI_BYTE,
		    &mpi_stat);
	    VRFY((mrc==MPI_SUCCESS), "");

	    /* move the offset pointer to last byte written by all processes */
	    mpi_off += (mpi_size - 1 - mpi_rank) * stride;

	    /* Increase chunk size without exceeding buffer size. */
	    /* Then move the starting offset for next write. */
	    stride *= 2;
	    if (stride > bufsize)
		stride = bufsize;
	    mpi_off += mpi_rank*stride;
	}

	/* close file and free the communicator */
	mrc = MPI_File_close(&fh);
	VRFY((mrc==MPI_SUCCESS), "MPI_FILE_CLOSE");
	mrc = MPI_Comm_free(&comm);
	VRFY((mrc==MPI_SUCCESS), "MPI_Comm_free");

	/* sync with the other waiting processes */
	mrc = MPI_Barrier(MPI_COMM_WORLD);
	VRFY((mrc==MPI_SUCCESS), "Sync after writes");
    }else{
	/* last process waits till writes are done,
	 * then opens file to verify data.
	 */
	mrc = MPI_Barrier(MPI_COMM_WORLD);
	VRFY((mrc==MPI_SUCCESS), "Sync after writes");

	mrc = MPI_File_open(comm, filename, MPI_MODE_RDONLY,
		info, &fh);
	VRFY((mrc==MPI_SUCCESS), "");

	stride = bufsize;
	for (mpi_off=0; mpi_off < MPIO_TEST_WRITE_SIZE; mpi_off += bufsize){
	    /* make sure it does not read beyond end of data */
	    if (mpi_off+stride > MPIO_TEST_WRITE_SIZE)
		stride = MPIO_TEST_WRITE_SIZE - mpi_off;
	    mrc = MPI_File_read_at(fh, mpi_off, buf, (int)stride, MPI_BYTE,
		    &mpi_stat);
	    VRFY((mrc==MPI_SUCCESS), "");
	    vrfyerrs=0;
	    for (i=0; i<stride; i++){
		unsigned char expected;
		expected = (unsigned char)(mpi_off+i);
		if ((expected != buf[i]) &&
		    (vrfyerrs++ < MAX_ERR_REPORT || VERBOSE_MED)) {
			printf("proc %d: found data error at [%ld], expect %u, got %u\n",
			    mpi_rank, (long)(mpi_off+i), expected, buf[i]);
		}
	    }
	    if (vrfyerrs > MAX_ERR_REPORT && !VERBOSE_MED)
		printf("proc %d: [more errors ...]\n", mpi_rank);

	    nerrs += vrfyerrs;
	}

	/* close file and free the communicator */
	mrc = MPI_File_close(&fh);
	VRFY((mrc==MPI_SUCCESS), "MPI_FILE_CLOSE");
	mrc = MPI_Comm_free(&comm);
	VRFY((mrc==MPI_SUCCESS), "MPI_Comm_free");
    }

    /*
     * one more sync to ensure all processes have done reading
     * before ending this test.
     */
    mrc = MPI_Barrier(MPI_COMM_WORLD);
    VRFY((mrc==MPI_SUCCESS), "Sync before leaving test");
    return (nerrs);
}


#define MB      1048576        /* 1024*1024 == 2**20 */
#define GB      1073741824     /* 1024**3 == 2**30 */
#define TWO_GB_LESS1    2147483647     /* 2**31 - 1 */
#define FOUR_GB_LESS1   4294967295L     /* 2**32 - 1 */
/*
 * Verify that MPI_Offset exceeding 2**31 can be computed correctly.
 * Print any failure as information only, not as an error so that this
 * won't abort the remaining test or other separated tests.
 *
 * Test if MPIO can write file from under 2GB to over 2GB and then
 * from under 4GB to over 4GB.
 * Each process writes 1MB in round robin fashion.
 * Then reads the file back in by reverse order, that is process 0
 * reads the data of process n-1 and vice versa.
 */
static int
test_mpio_gb_file(char *filename)
{
    int mpi_size, mpi_rank;
    MPI_Info info = MPI_INFO_NULL;
    int mrc;
    MPI_File	fh;
    int i, j, n;
    int vrfyerrs;
    int writerrs;		/* write errors */
    int nerrs;
    int ntimes;			/* how many times */
    char  *buf = NULL;
    char  expected;
    MPI_Offset  mpi_off;
    MPI_Offset  mpi_off_old;
    MPI_Status  mpi_stat;
    int is_signed, sizeof_mpi_offset;

    nerrs = 0;
    /* set up MPI parameters */
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

    if (VERBOSE_MED)
        printf("MPI_Offset range test\n");

    /* figure out the signness and sizeof MPI_Offset */
    mpi_off = 0;
    is_signed = ((MPI_Offset)(mpi_off - 1)) < 0;
    sizeof_mpi_offset = (int)(sizeof(MPI_Offset));

    /*
     * Verify the sizeof MPI_Offset and correctness of handling multiple GB
     * sizes.
     */
    if (MAINPROCESS){			/* only process 0 needs to check it*/
	printf("MPI_Offset is %s %d bytes integeral type\n",
	    is_signed ? "signed" : "unsigned", (int)sizeof(MPI_Offset));
	if (sizeof_mpi_offset <= 4 && is_signed){
	    printf("Skipped 2GB range test "
		    "because MPI_Offset cannot support it\n");
	}else {
	    /* verify correctness of assigning 2GB sizes */
	    mpi_off = 2 * 1024 * (MPI_Offset)MB;
	    INFO((mpi_off>0), "2GB OFFSET assignment no overflow");
	    INFO((mpi_off-1)==TWO_GB_LESS1, "2GB OFFSET assignment succeed");

	    /* verify correctness of increasing from below 2 GB to above 2GB */
	    mpi_off = TWO_GB_LESS1;
	    for (i=0; i < 3; i++){
		mpi_off_old = mpi_off;
		mpi_off = mpi_off + 1;
		/* no overflow */
		INFO((mpi_off>0), "2GB OFFSET increment no overflow");
		/* correct inc. */
		INFO((mpi_off-1)==mpi_off_old, "2GB OFFSET increment succeed");
	    }
	}

	if (sizeof_mpi_offset <= 4){
	    printf("Skipped 4GB range test "
		    "because MPI_Offset cannot support it\n");
	}else {
	    /* verify correctness of assigning 4GB sizes */
	    mpi_off = 4 * 1024 * (MPI_Offset)MB;
	    INFO((mpi_off>0), "4GB OFFSET assignment no overflow");
	    INFO((mpi_off-1)==FOUR_GB_LESS1, "4GB OFFSET assignment succeed");

	    /* verify correctness of increasing from below 4 GB to above 4 GB */
	    mpi_off = FOUR_GB_LESS1;
	    for (i=0; i < 3; i++){
		mpi_off_old = mpi_off;
		mpi_off = mpi_off + 1;
		/* no overflow */
		INFO((mpi_off>0), "4GB OFFSET increment no overflow");
		/* correct inc. */
		INFO((mpi_off-1)==mpi_off_old, "4GB OFFSET increment succeed");
	    }
	}
    }

    /*
     * Verify if we can write to a file of multiple GB sizes.
     */
    if (VERBOSE_MED)
	printf("MPIO GB file test %s\n", filename);

    if (sizeof_mpi_offset <= 4){
	printf("Skipped GB file range test "
		"because MPI_Offset cannot support it\n");
    }else{
	buf = malloc(MB);
	VRFY((buf!=NULL), "malloc succeed");

	/* open a new file. Remove it first in case it exists. */
	if (MAINPROCESS)
	    remove(filename);
	MPI_Barrier(MPI_COMM_WORLD);	/* prevent racing condition */

	mrc = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR,
		    info, &fh);
	VRFY((mrc==MPI_SUCCESS), "MPI_FILE_OPEN");

	printf("MPIO GB file write test %s\n", filename);

	/* instead of writing every bytes of the file, we will just write
	 * some data around the 2 and 4 GB boundaries.  That should cover
	 * potential integer overflow and filesystem size limits.
	 */
	writerrs = 0;
	for (n=2; n <= 4; n+=2){
	    ntimes = GB/MB*n/mpi_size + 1;
	    for (i=ntimes-2; i <= ntimes; i++){
		mpi_off = (i*mpi_size + mpi_rank)*(MPI_Offset)MB;
		if (VERBOSE_MED)
		    HDfprintf(stdout,"proc %d: write to mpi_off=%016llx, %lld\n",
			mpi_rank, mpi_off, mpi_off);
		/* set data to some trivial pattern for easy verification */
		for (j=0; j<MB; j++)
		    *(buf+j) = i*mpi_size + mpi_rank;
		if (VERBOSE_MED)
		    HDfprintf(stdout,"proc %d: writing %d bytes at offset %lld\n",
			mpi_rank, MB, mpi_off);
		mrc = MPI_File_write_at(fh, mpi_off, buf, MB, MPI_BYTE, &mpi_stat);
		INFO((mrc==MPI_SUCCESS), "GB size file write");
		if (mrc!=MPI_SUCCESS)
		    writerrs++;
	    }
	}

	/* close file and free the communicator */
	mrc = MPI_File_close(&fh);
	VRFY((mrc==MPI_SUCCESS), "MPI_FILE_CLOSE");

	mrc = MPI_Barrier(MPI_COMM_WORLD);
	VRFY((mrc==MPI_SUCCESS), "Sync after writes");

	/*
	 * Verify if we can read the multiple GB file just created.
	 */
	/* open it again to verify the data written */
	/* but only if there was no write errors */
	printf("MPIO GB file read test %s\n", filename);
	if (errors_sum(writerrs)>0){
	    printf("proc %d: Skip read test due to previous write errors\n",
		mpi_rank);
	    goto finish;
	}
	mrc = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, info, &fh);
	VRFY((mrc==MPI_SUCCESS), "");

	/* Only read back parts of the file that have been written. */
	for (n=2; n <= 4; n+=2){
	    ntimes = GB/MB*n/mpi_size + 1;
	    for (i=ntimes-2; i <= ntimes; i++){
		mpi_off = (i*mpi_size + (mpi_size - mpi_rank - 1))*(MPI_Offset)MB;
		if (VERBOSE_MED)
		    HDfprintf(stdout,"proc %d: read from mpi_off=%016llx, %lld\n",
			mpi_rank, mpi_off, mpi_off);
		mrc = MPI_File_read_at(fh, mpi_off, buf, MB, MPI_BYTE, &mpi_stat);
		INFO((mrc==MPI_SUCCESS), "GB size file read");
		expected = i*mpi_size + (mpi_size - mpi_rank - 1);
		vrfyerrs=0;
		for (j=0; j<MB; j++){
		    if ((*(buf+j) != expected) &&
			(vrfyerrs++ < MAX_ERR_REPORT || VERBOSE_MED)){
			    printf("proc %d: found data error at [%ld+%d], expect %d, got %d\n",
				mpi_rank, (long)mpi_off, j, expected, *(buf+j));
		    }
		}
		if (vrfyerrs > MAX_ERR_REPORT && !VERBOSE_MED)
		    printf("proc %d: [more errors ...]\n", mpi_rank);

		nerrs += vrfyerrs;
	    }
	}

	/* close file and free the communicator */
	mrc = MPI_File_close(&fh);
	VRFY((mrc==MPI_SUCCESS), "MPI_FILE_CLOSE");

	/*
	 * one more sync to ensure all processes have done reading
	 * before ending this test.
	 */
	mrc = MPI_Barrier(MPI_COMM_WORLD);
	VRFY((mrc==MPI_SUCCESS), "Sync before leaving test");
    }

finish:
    if (buf)
	HDfree(buf);
    return (nerrs);
}


/*
 * MPI-IO Test: One writes, Many reads.
 * Verify if only one process writes some data and then all other
 * processes can read them back correctly. This tests if the
 * underlaying parallel I/O and file system supports parallel I/O
 * correctly.
 *
 * Algorithm: Only one process (e.g., process 0) writes some data.
 * Then all processes, including the writing process, read the data
 * back and verify them against the original values.
 */

/*
 * Default filename can be specified via first program argument.
 * Each process writes something, then reads all data back.
 */

#define DIMSIZE	32		/* Dimension size. */
#define PRINTID printf("Proc %d: ", mpi_rank)
#define USENONE 0
#define USEATOM 1		/* request atomic I/O */
#define USEFSYNC 2		/* request file_sync */


static int
test_mpio_1wMr(char *filename, int special_request)
{
    char hostname[128];
    int  mpi_size, mpi_rank;
    MPI_File fh;
    char mpi_err_str[MPI_MAX_ERROR_STRING];
    int  mpi_err_strlen;
    int  mpi_err;
    unsigned char writedata[DIMSIZE], readdata[DIMSIZE];
    unsigned char expect_val;
    int  i, irank;
    int  nerrs = 0;		/* number of errors */
    int  atomicity;
    MPI_Offset  mpi_off;
    MPI_Status  mpi_stat;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (MAINPROCESS && VERBOSE_MED){
        printf("Testing one process writes, all processes read.\n");
	printf("Using %d processes accessing file %s\n", mpi_size, filename);
        printf("    (Filename can be specified via program argument)\n");
    }

    /* show the hostname so that we can tell where the processes are running */
    if (VERBOSE_DEF){
	if (gethostname(hostname, 128) < 0){
	    PRINTID;
	    printf("gethostname failed\n");
	    return 1;
	}
	PRINTID;
	printf("hostname=%s\n", hostname);
    }

    /* Delete any old file in order to start anew. */
    /* Must delete because MPI_File_open does not have a Truncate mode. */
    /* Don't care if it has error. */
    MPI_File_delete(filename, MPI_INFO_NULL);

    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, filename,
	    MPI_MODE_RDWR | MPI_MODE_CREATE ,
	    MPI_INFO_NULL, &fh))
	    != MPI_SUCCESS){
	MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	PRINTID;
	printf("MPI_File_open failed (%s)\n", mpi_err_str);
	return 1;
    }

if (special_request & USEATOM){
    /* ==================================================
     * Set atomcity to true (1).  A POSIX compliant filesystem
     * should not need this.
     * ==================================================*/
    if ((mpi_err = MPI_File_get_atomicity(fh, &atomicity)) != MPI_SUCCESS){
	MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	PRINTID;
	printf("MPI_File_get_atomicity failed (%s)\n", mpi_err_str);
    }
    if (VERBOSE_HI)
	printf("Initial atomicity = %d\n", atomicity);
    if ((mpi_err = MPI_File_set_atomicity(fh, 1)) != MPI_SUCCESS){
	MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	PRINTID;
	printf("MPI_File_set_atomicity failed (%s)\n", mpi_err_str);
    }
    if ((mpi_err = MPI_File_get_atomicity(fh, &atomicity)) != MPI_SUCCESS){
	MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	PRINTID;
	printf("MPI_File_get_atomicity failed (%s)\n", mpi_err_str);
    }
    if (VERBOSE_HI)
	printf("After set_atomicity atomicity = %d\n", atomicity);
}

    /* This barrier is not necessary but do it anyway. */
    MPI_Barrier(MPI_COMM_WORLD);
    if (VERBOSE_HI){
	PRINTID;
	printf("between MPI_Barrier and MPI_File_write_at\n");
    }

    /* ==================================================
     * Each process calculates what to write but
     * only process irank(0) writes.
     * ==================================================*/
    irank=0;
    for (i=0; i < DIMSIZE; i++)
	writedata[i] = irank*DIMSIZE + i;
    mpi_off = irank*DIMSIZE;

    /* Only one process writes */
    if (mpi_rank==irank){
	if (VERBOSE_HI){
	    PRINTID; printf("wrote %d bytes at %ld\n", DIMSIZE, (long)mpi_off);
	}
	if ((mpi_err = MPI_File_write_at(fh, mpi_off, writedata, DIMSIZE,
			MPI_BYTE, &mpi_stat))
		!= MPI_SUCCESS){
	    MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	    PRINTID;
	    printf("MPI_File_write_at offset(%ld), bytes (%d), failed (%s)\n",
		    (long) mpi_off, DIMSIZE, mpi_err_str);
	    return 1;
	};
    };

    /* Bcast the return code and */
    /* make sure all writing are done before reading. */
    MPI_Bcast(&mpi_err, 1, MPI_INT, irank, MPI_COMM_WORLD);
    if (VERBOSE_HI){
	PRINTID;
	printf("MPI_Bcast: mpi_err = %d\n", mpi_err);
    }

if (special_request & USEFSYNC){
    /* ==================================================
     * Do a file sync.  A POSIX compliant filesystem
     * should not need this.
     * ==================================================*/
    if (VERBOSE_HI)
	printf("Apply MPI_File_sync\n");
    /* call file_sync to force the write out */
    if ((mpi_err = MPI_File_sync(fh)) != MPI_SUCCESS){
	MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	PRINTID;
	printf("MPI_File_sync failed (%s)\n", mpi_err_str);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /* call file_sync to force the write out */
    if ((mpi_err = MPI_File_sync(fh)) != MPI_SUCCESS){
	MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	PRINTID;
	printf("MPI_File_sync failed (%s)\n", mpi_err_str);
    }
}

    /* This barrier is not necessary because the Bcase or File_sync above */
    /* should take care of it.  Do it anyway. */
    MPI_Barrier(MPI_COMM_WORLD);
    if (VERBOSE_HI){
	PRINTID;
	printf("after MPI_Barrier\n");
    }

    /* ==================================================
     * Each process reads what process 0 wrote and verify.
     * ==================================================*/
    irank=0;
    mpi_off = irank*DIMSIZE;
    if ((mpi_err = MPI_File_read_at(fh, mpi_off, readdata, DIMSIZE, MPI_BYTE,
	    &mpi_stat))
	    != MPI_SUCCESS){
	MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	PRINTID;
	printf("MPI_File_read_at offset(%ld), bytes (%d), failed (%s)\n",
		(long) mpi_off, DIMSIZE, mpi_err_str);
	return 1;
    };
    for (i=0; i < DIMSIZE; i++){
	expect_val = irank*DIMSIZE + i;
	if (readdata[i] != expect_val){
	    PRINTID;
	    printf("read data[%d:%d] got %02x, expect %02x\n", irank, i,
		    readdata[i], expect_val);
	    nerrs++;
	}
    }

    MPI_File_close(&fh);

    if (VERBOSE_HI){
	PRINTID;
	printf("%d data errors detected\n", nerrs);
    }

    mpi_err = MPI_Barrier(MPI_COMM_WORLD);
    return nerrs;
}


/*
 * parse the command line options
 */
static int
parse_options(int argc, char **argv)
{
    while (--argc){
	if (**(++argv) != '-'){
	    break;
	}else{
	    switch(*(*argv+1)){
		case 'v':   if (*((*argv+1)+1))
				ParseTestVerbosity((*argv+1)+1);
			    else
				SetTestVerbosity(VERBO_MED);
			    break;
		case 'f':   if (--argc < 1) {
				nerrors++;
				return(1);
			    }
			    if (**(++argv) == '-') {
				nerrors++;
				return(1);
			    }
			    paraprefix = *argv;
			    break;
		case 'h':   /* print help message--return with nerrors set */
			    return(1);
		default:    nerrors++;
			    return(1);
	    }
	}
    } /*while*/

    /* compose the test filenames */
    {
	int i, n;
	hid_t plist;

	plist = H5Pcreate (H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
	n = sizeof(FILENAME)/sizeof(FILENAME[0]) - 1;	/* exclude the NULL */

	for (i=0; i < n; i++)
	    if (h5_fixname(FILENAME[i],plist,filenames[i],sizeof(filenames[i]))
		== NULL){
		printf("h5_fixname failed\n");
		nerrors++;
		return(1);
	    }
	H5Pclose(plist);
	if (VERBOSE_MED){
	    printf("Test filenames are:\n");
	    for (i=0; i < n; i++)
		printf("    %s\n", filenames[i]);
	}
    }

    return(0);
}


/*
 * Show command usage
 */
static void
usage(void)
{
    printf("Usage: t_mpi [-v<verbosity>] [-f <prefix>]\n");
    printf("\t-v<verbosity>\tset verbose level (0-9,l,m,h)\n");
    printf("\t-f <prefix>\tfilename prefix\n");
    printf("\n");
}

/*
 * return the sum of all errors.
 */
static int
errors_sum(int nerrs)
{
    int temp;
    MPI_Allreduce(&nerrs, &temp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return(temp);
}


int
main(int argc, char **argv)
{
    int mpi_size, mpi_rank;				/* mpi variables */
    int ret_code;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    H5open();
    if (parse_options(argc, argv) != 0){
	if (MAINPROCESS)
	    usage();
	goto finish;
    }

    if (MAINPROCESS){
	printf("===================================\n");
	printf("MPI functionality tests\n");
	printf("===================================\n");
    }

    if (VERBOSE_MED)
	h5_show_hostname();

    fapl = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);

    MPI_BANNER("MPIO 1 write Many read test...");
    ret_code = test_mpio_1wMr(filenames[0], USENONE);
    ret_code = errors_sum(ret_code);
    if (mpi_rank==0 && ret_code > 0){
	printf("***FAILED with %d total errors\n", ret_code);
	nerrors += ret_code;
    }

    /* test atomicity and file sync in high verbose mode only         */
    /* since they often hang when broken and PHDF5 does not use them. */
    if (VERBOSE_HI){
	MPI_BANNER("MPIO 1 write Many read test with atomicity...");
	ret_code = test_mpio_1wMr(filenames[0], USEATOM);
	ret_code = errors_sum(ret_code);
	if (mpi_rank==0 && ret_code > 0){
	    printf("***FAILED with %d total errors\n", ret_code);
	    nerrors += ret_code;
	}

	MPI_BANNER("MPIO 1 write Many read test with file sync...");
	ret_code = test_mpio_1wMr(filenames[0], USEFSYNC);
	ret_code = errors_sum(ret_code);
	if (mpi_rank==0 && ret_code > 0){
	    printf("***FAILED with %d total errors\n", ret_code);
	    nerrors += ret_code;
	}
    }

    MPI_BANNER("MPIO File size range test...");
    ret_code = test_mpio_gb_file(filenames[0]);
    ret_code = errors_sum(ret_code);
    if (mpi_rank==0 && ret_code > 0){
	printf("***FAILED with %d total errors\n", ret_code);
	nerrors += ret_code;
    }

    MPI_BANNER("MPIO independent overlapping writes...");
    ret_code = test_mpio_overlap_writes(filenames[0]);
    ret_code = errors_sum(ret_code);
    if (mpi_rank==0 && ret_code > 0){
	printf("***FAILED with %d total errors\n", ret_code);
	nerrors += ret_code;
    }

finish:
    /* make sure all processes are finished before final report, cleanup
     * and exit.
     */
    MPI_Barrier(MPI_COMM_WORLD);
    if (MAINPROCESS){		/* only process 0 reports */
	printf("===================================\n");
	if (nerrors){
	    printf("***MPI tests detected %d errors***\n", nerrors);
	}
	else{
	    printf("MPI tests finished with no errors\n");
	}
	printf("===================================\n");
    }

    h5_cleanup(FILENAME, fapl);
    H5close();

    /* MPI_Finalize must be called AFTER H5close which may use MPI calls */
    MPI_Finalize();

    /* cannot just return (nerrors) because exit code is limited to 1byte */
    return(nerrors!=0);
}

