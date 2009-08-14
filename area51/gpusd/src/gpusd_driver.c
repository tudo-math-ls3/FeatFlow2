#include "gpusd_driver.h"
#include <math.h>
#include <cuda.h>
//#include <stdio.h>

#if defined(WIN32) || defined(WIN64)
// we'll need this for the high-resolution counter
#include <windows.h>
#define TIMESTAMP LARGE_INTEGER
#define GET_TIMESTAMP(x) QueryPerformanceCounter(&x)
#endif

#define BATCH_SIZE		512

#define _min(i,j) ((i) < (j) ? (i) : (j))

typedef float real;

/**
 * \brief
 */
typedef struct _t_gpusd
{
	// scalar parameters
	double		dnu, dalpha, dbeta, dgamma, ddelta;

	// triangulation arrays
	int			nel;
	int			*Iverts, *Iedges;
	double		*Dvtx;

	// matrix arrays
	int			neq;
	int			*Irowptr, *Icolidx;
	double		*Ddata;

	// velocity field arrays
	double		*Du1, *Du2;

	// time for memcopy: host -> device
	long long	lltimeMemCpyHtoD_CPU;
	long long	lltimeMemCpyHtoD_GPU;
	long long	lltimeMemCpyDtoH_CPU;
	long long	lltimeMemCpyDtoH_GPU;

	// time for kernel
	long long	lltimeKernel_CPU;
	long long	lltimeKernel_GPU;

	// time for encoding/decoding
	long long	lltimeEncode_CPU;
	long long	lltimeDecode_CPU;

	// total assembly time
	long long	lltimeTotal_CPU;
	long long	lltimeTotal_GPU;

	// CUDA device
	CUdevice	cuDevice;

	// CUDA context
	CUcontext	cuContext;

	// CUDA module
	CUmodule	cuModule;

	// CUDA function handle for SD kernel
	CUfunction	cuSDkernel;

	// CUDA device pointer for "g_upsam"
	CUdeviceptr cuGupsam;

	// CUDA device pointer for "g_nu"
	CUdeviceptr cuGnu;

	// CUDA device pointer for "g_maxu"
	CUdeviceptr cuGmaxu;

	// CUDA event for timing
	// 0..1: Memcopy H->D
	// 1..2: kernel launch
	// 2..3: Memcopy D->H
	CUevent		cuEvent0, cuEvent1, cuEvent2, cuEvent3;

	// CPU time stamps
	TIMESTAMP	stamp0, stamp1, stamp2, stamp3;
	TIMESTAMP	stampBegin, stampEnd;

} t_gpusd;

static t_gpusd g_data;

typedef struct _t_buffers
{
	// input/output sizes
	int				nel, batchsize;
	int				neltodo_in,neltodo_out;
	int				neldone_in,neldone_out;

	// index of front buffer
	int				ifront;

	// size of buffers in bytes
	int				bufsize;

	// host memory buffers (double-buffered)
	real			*io[2];

	// device memory buffer
	CUdeviceptr		dev_io;
} t_buffers;

static t_buffers g_buf;

CUresult rtn;

void memformat(void *buf, size_t len)
{
	char *p = (char*)buf;
	size_t i;

	for(i = 0; i < len; i++)
		p[i] = 0;
};


void gpusd_init(int* info)
{
	int deviceCount = 0;

	if(info == NULL)
		return;

	// format the data structure
	memformat(&g_data, sizeof(t_gpusd));

	// try to initialize CUDA
	if(cuInit(0) != CUDA_SUCCESS)
	{
		*info = GPUSD_INFO_CUDA_FAIL;
		return;
	}

	// let's see whether we have a CUDA device
	cuDeviceGetCount(&deviceCount);
	if(deviceCount == 0)
	{
		*info = GPUSD_INFO_NO_DEVICE;
		return;
	}

	// get a handle for device 0
	rtn = cuDeviceGet(&g_data.cuDevice, 0);

	// get the maximum number of threads in x-direction
	//cuDeviceGetAttribute(&g_data.imaxthreads, CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X,
	//	g_data.cuDevice);

	// create a context
	rtn = cuCtxCreate(&g_data.cuContext, CU_CTX_SCHED_YIELD, g_data.cuDevice);
	//rtn = cuCtxCreate(&g_data.cuContext, CU_CTX_BLOCKING_SYNC, g_data.cuDevice);

	// create the module
	//rtn = cuModuleLoad(&g_data.cuModule, "./ptx/gpusd.ptx");
	rtn = cuModuleLoad(&g_data.cuModule, "./ptx/gpusd_sp.ptx");

	// get a handle to the kernel
	rtn = cuModuleGetFunction(&g_data.cuSDkernel, g_data.cuModule, "asm_kernel");

	// That' it!
	*info = GPUSD_INFO_SUCCESS;
}

void gpusd_done()
{
	// unload the module
	cuModuleUnload(g_data.cuModule);
	// destroy the context
	cuCtxDestroy(g_data.cuContext);
}

void gpusd_setparams(double *dnu, double *dalpha, double *dbeta,
					 double *dgamma, double *ddelta)
{
	g_data.dnu    = *dnu;
	g_data.dalpha = *dalpha;
	g_data.dbeta  = *dbeta;
	g_data.dgamma = *dgamma;
	g_data.ddelta = *ddelta;
}

void gpusd_settria(int *nel, double *Dvtx, int *Iverts, int *Iedges)
{
	g_data.nel = *nel;
	g_data.Dvtx = Dvtx;
	g_data.Iverts = Iverts;
	g_data.Iedges = Iedges;
}

void gpusd_setmatrix(int *neq, int *Irowptr, int *Icolidx, double *Ddata)
{
	g_data.neq = *neq;
	g_data.Irowptr = Irowptr;
	g_data.Icolidx = Icolidx;
	g_data.Ddata = Ddata;
}

void gpusd_setvelo(double *Du1, double *Du2)
{
	g_data.Du1 = Du1;
	g_data.Du2 = Du2;
}

long long calcCPUTime(TIMESTAMP tBegin, TIMESTAMP tEnd)
{
#if defined(WIN32) || defined(WIN64)
	LARGE_INTEGER freq;
	long long t;
	
	// get the timer frequency
	if(QueryPerformanceFrequency(&freq) == FALSE)
		return 0L;

	// calculate time difference in microseconds
	t = tEnd.QuadPart - tBegin.QuadPart;
	return (1000000L*t)/freq.QuadPart;
#else
	return 0L;
#endif
};

double calc_maxu()
{
	double dt,du;
	int i;

	du = 0.0;
	for(i = 0; i < g_data.neq; i++)
	{
		dt = g_data.Du1[i]*g_data.Du1[i] + g_data.Du2[i]*g_data.Du2[i];
		if(dt > du)
			du = dt;
	}

	return sqrt(dt);
};

void set_constants(double dmaxu)
{
	real ft;
	unsigned int bytes;
	CUdeviceptr devptr;

	ft = (real)g_data.dnu;
	rtn = cuModuleGetGlobal(&devptr, &bytes, g_data.cuModule, "g_nu");
	rtn = cuMemcpyHtoD(devptr, &ft, bytes);

	ft = (real)dmaxu;
	rtn = cuModuleGetGlobal(&devptr, &bytes, g_data.cuModule, "g_maxu");
	rtn = cuMemcpyHtoD(devptr, &ft, bytes);

	ft = (real)g_data.dalpha;
	rtn = cuModuleGetGlobal(&devptr, &bytes, g_data.cuModule, "g_alpha");
	rtn = cuMemcpyHtoD(devptr, &ft, bytes);

	ft = (real)g_data.dbeta;
	rtn = cuModuleGetGlobal(&devptr, &bytes, g_data.cuModule, "g_beta");
	rtn = cuMemcpyHtoD(devptr, &ft, bytes);

	ft = (real)g_data.dgamma;
	rtn = cuModuleGetGlobal(&devptr, &bytes, g_data.cuModule, "g_gamma");
	rtn = cuMemcpyHtoD(devptr, &ft, bytes);

	ft = (real)g_data.ddelta;
	rtn = cuModuleGetGlobal(&devptr, &bytes, g_data.cuModule, "g_delta");
	rtn = cuMemcpyHtoD(devptr, &ft, bytes);
};

void buf_init()
{
	memformat(&g_buf, sizeof(t_buffers));

	// set up the batch size
	g_buf.nel = g_data.nel;
	//g_buf.batchsize = __min(g_data.imaxthreads, g_buf.nel);
	//g_buf.batchsize = __min(g_buf.batchsize, BATCH_SIZE);
	g_buf.batchsize = __min(g_buf.nel, BATCH_SIZE);

	// calculate buffer size
	g_buf.bufsize = sizeof(real) * 16 * BATCH_SIZE;

	// try to allocate some host memory
	rtn = cuMemAllocHost(&g_buf.io[0], g_buf.bufsize);
	rtn = cuMemAllocHost(&g_buf.io[1], g_buf.bufsize);

	// also, try to allocate some memory on the device
	rtn = cuMemAlloc(&g_buf.dev_io, g_buf.bufsize);
};

void buf_done()
{
	// release device memory
	if(g_buf.dev_io != 0)
		cuMemFree(g_buf.dev_io);

	// release host memory
	if(g_buf.io[1] != NULL)
		cuMemFreeHost(g_buf.io[1]);
	if(g_buf.io[0] != NULL)
		cuMemFreeHost(g_buf.io[0]);
};

void encode_input()
{
	int i,j,k,l;
	real* buf;

	// create CPU timestamp
	GET_TIMESTAMP(g_data.stamp0);

	// calculate the number of elements to process
	g_buf.neltodo_in = __min(g_buf.batchsize, g_buf.nel - g_buf.neldone_in);

	// fetch the front buffer
	buf = g_buf.io[g_buf.ifront];

	// set up input buffers
	for(i = 0; i < g_buf.neltodo_in; i++)
	{
		for(j = 0; j < 4; j++)
		{
			/*
			l = 16*i+2*j;
			// vertex coordinates
			k = 2*(g_data.Iverts[4*(g_buf.neldone_in+i)+j]-1);
			buf[l+0] = (real)g_data.Dvtx[k+0];
			buf[l+1] = (real)g_data.Dvtx[k+1];
			// dof values
			k = g_data.Iedges[4*(g_buf.neldone_in+i)+j]-1;
			buf[l+8] = (real)g_data.Du1[k];
			buf[l+9] = (real)g_data.Du2[k];
			//*/
			//*
			l = 256*(i>>4) + (i&0xF) + 32*j;
			// vertex coordinates
			k = 2*(g_data.Iverts[4*(g_buf.neldone_in+i)+j]-1);
			buf[l+  0] = (real)g_data.Dvtx[k+0];
			buf[l+ 16] = (real)g_data.Dvtx[k+1];
			// dof values
			k = g_data.Iedges[4*(g_buf.neldone_in+i)+j]-1;
			buf[l+128] = (real)g_data.Du1[k];
			buf[l+144] = (real)g_data.Du2[k];
			//*/
		}
	}

	// increase counter
	g_buf.neldone_in += g_buf.neltodo_in;

	// remember how much elements we need to decode next time
	g_buf.neltodo_out = g_buf.neltodo_in;

	// create CPU timespamp
	GET_TIMESTAMP(g_data.stamp1);
	g_data.lltimeEncode_CPU += calcCPUTime(g_data.stamp0, g_data.stamp1);
};

void decode_output()
{
	int i,j,k,l,*d;
	real *buf,*A;

	// fetch the front buffer
	buf = g_buf.io[g_buf.ifront];

	// create CPU timestamp
	GET_TIMESTAMP(g_data.stamp0);

	// loop over all elements
	for(k = 0; k < g_buf.neltodo_out; k++)
	{
		// get a pointer to the local matrix
		//A = &buf[16*k];
		
		// get a pointer to the local dofs
		d = &g_data.Iedges[4*(g_buf.neldone_out+k)];

		// loop over the test functions
		for(i = 0; i < 4; i++)
		{
			A = &buf[256*(k>>4) + (k&0xF) + 64*i];

			// loop over the trial functions
			for(j = 0; j < 4; j++)
			{
				// loop over the global matrix row
				for(l = g_data.Irowptr[d[i]-1]; l < g_data.Irowptr[d[i]]; l++)
				{
					// is this our trial function?
					if(g_data.Icolidx[l-1] == d[j])
					{
						// write to output
						//g_data.Ddata[l-1] += (double)A[4*i+j];
						g_data.Ddata[l-1] += (double)A[j*16];
						break;
					}
				}
			}
		}
	}
	// create CPU timespamp
	GET_TIMESTAMP(g_data.stamp1);
	g_data.lltimeDecode_CPU += calcCPUTime(g_data.stamp0, g_data.stamp1);

	// update counter
	g_buf.neldone_out += g_buf.neltodo_out;
};

void wait_for_sync()
{
	float ft;

	rtn = cuCtxSynchronize();

	cuEventElapsedTime(&ft, g_data.cuEvent0, g_data.cuEvent1);
	g_data.lltimeMemCpyHtoD_GPU += (long long)(ft*1000.0f);
	cuEventElapsedTime(&ft, g_data.cuEvent1, g_data.cuEvent2);
	g_data.lltimeKernel_GPU += (long long)(ft*1000.0f);
	cuEventElapsedTime(&ft, g_data.cuEvent2, g_data.cuEvent3);
	g_data.lltimeMemCpyDtoH_GPU += (long long)(ft*1000.0f);
	cuEventElapsedTime(&ft, g_data.cuEvent0, g_data.cuEvent3);
	g_data.lltimeTotal_GPU += (long long)(ft*1000.0f);
};

void launch_kernel()
{
	int iback;
	int nthreads, nblocks;

	// switch front and back buffers
	iback = g_buf.ifront;
	g_buf.ifront = iback ^ 0x1;

	// set up the parameters
	rtn = cuParamSeti(g_data.cuSDkernel, 0, (unsigned int)g_buf.dev_io);
	rtn = cuParamSetSize(g_data.cuSDkernel, 32);

	// calculate sizes
	/*if(g_buf.neltodo_in <= 512)
	{
		nthreads = g_buf.neltodo_in;
		nblocks = 1;
	}
	else
	{
		nthreads = 512;
		nblocks = (g_buf.neltodo_in / 512) + 
			((g_buf.neltodo_in % 512) > 0 ? 1 : 0);
	}*/

	// set block shape
	rtn = cuFuncSetBlockShape(g_data.cuSDkernel, g_buf.neltodo_in, 1, 1);

	GET_TIMESTAMP(g_data.stamp0);

	// perform asynchronous memcopy: H -> D
	rtn = cuEventRecord(g_data.cuEvent0, 0);
	rtn = cuMemcpyHtoDAsync(g_buf.dev_io, g_buf.io[iback], g_buf.bufsize,0);
	GET_TIMESTAMP(g_data.stamp1);

	// launch
	rtn = cuEventRecord(g_data.cuEvent1, 0);
	rtn = cuLaunchGridAsync(g_data.cuSDkernel, 1, 1, 0);
	rtn = cuEventRecord(g_data.cuEvent2, 0);
	GET_TIMESTAMP(g_data.stamp2);

	// perform asynchronous memcopy: D -> H
	rtn = cuMemcpyDtoHAsync(g_buf.io[iback], g_buf.dev_io, g_buf.bufsize,0);
	rtn = cuEventRecord(g_data.cuEvent3, 0);
	GET_TIMESTAMP(g_data.stamp3);

	g_data.lltimeMemCpyHtoD_CPU += calcCPUTime(g_data.stamp0, g_data.stamp1);
	g_data.lltimeKernel_CPU += calcCPUTime(g_data.stamp1, g_data.stamp2);
	g_data.lltimeMemCpyDtoH_CPU += calcCPUTime(g_data.stamp2, g_data.stamp3);
};

void gpusd_assemble(int *info)
{
	float ft;
	if(info == NULL)
		return;

	// clear all timing results
	g_data.lltimeTotal_CPU = g_data.lltimeTotal_GPU = 0L;
	g_data.lltimeEncode_CPU = g_data.lltimeDecode_CPU = 0L;
	g_data.lltimeKernel_CPU = g_data.lltimeKernel_GPU= 0L;
	g_data.lltimeMemCpyHtoD_CPU = g_data.lltimeMemCpyHtoD_GPU = 0L;
	g_data.lltimeMemCpyDtoH_CPU = g_data.lltimeMemCpyDtoH_GPU = 0L;

	// first of all, let's see whether we have all necessary data
	*info = GPUSD_INFO_NOT_READY;
	if(g_data.dnu == 0.0)
		return;
	if((g_data.nel <= 0) || (g_data.Dvtx == NULL) || (g_data.Iverts == NULL) 
		|| (g_data.Iedges == NULL))
		return;
	if((g_data.neq <= 0) || (g_data.Irowptr == NULL) || (g_data.Icolidx == NULL)
		|| (g_data.Ddata == NULL))
		return;
	if((g_data.Du1 == NULL) || (g_data.Du2 == NULL))
		return;
	//*
	// create CUDA events for timing
	cuEventCreate(&g_data.cuEvent0, CU_EVENT_DEFAULT);
	cuEventCreate(&g_data.cuEvent1, CU_EVENT_DEFAULT);
	cuEventCreate(&g_data.cuEvent2, CU_EVENT_DEFAULT);
	cuEventCreate(&g_data.cuEvent3, CU_EVENT_DEFAULT);

	GET_TIMESTAMP(g_data.stampBegin);

	// Okay, if theta is zero, then we don't have to do anything
	*info = GPUSD_INFO_SUCCESS;

	// set the constants
	set_constants(calc_maxu());

	// set up the buffers
	buf_init();

	// load the initial batch onto the GPU
	encode_input();

	// launch the kernel
	launch_kernel();

	// continue until all elements have been processed
	while(g_buf.neldone_in < g_buf.nel)
	{
		// encode the next input buffer
		encode_input();

		// sync the buffers (read and write)
		wait_for_sync();

		// launch the kernel
		launch_kernel();

		// decode output buffer
		decode_output();
	}

	// sync final (read only)
	wait_for_sync();
	// switch buffers manually
	g_buf.ifront = g_buf.ifront ^ 0x1;

	// decode final output buffer
	decode_output();

	// release the buffer
	buf_done();

	cuEventDestroy(g_data.cuEvent3);
	cuEventDestroy(g_data.cuEvent2);
	cuEventDestroy(g_data.cuEvent1);
	cuEventDestroy(g_data.cuEvent0);

	// create final time stamp
	GET_TIMESTAMP(g_data.stampEnd);

	// calculate total time
	g_data.lltimeTotal_CPU = calcCPUTime(g_data.stampBegin, g_data.stampEnd);

	// That's it
	*info = GPUSD_INFO_SUCCESS;
}

void gpusd_statistics()
{
	//*
	double dGFLOPs,dGBs_HtoD_CPU,dGBs_HtoD_GPU,dGBs_DtoH_CPU,dGBs_DtoH_GPU;

	// calculate kernel GFLOP/s
	// GFLOPs = 2260*NEL / (1000 * timeKernel)
	dGFLOPs = (2.26*(double)g_data.nel) / (double)g_data.lltimeKernel_GPU;
	//dGBs_HtoD_CPU = (64.0*(double)g_data.nel) / (double)g_data.lltimeMemCpyHtoD_CPU;
	//dGBs_HtoD_CPU = (64.0*(double)g_data.nel) / (double)g_data.lltimeMemCpyHtoD_GPU;
	//dGBs_DtoH_CPU = (64.0*(double)g_data.nel) / (double)g_data.lltimeMemCpyDtoH_CPU;
	//dGBs_DtoH_CPU = (64.0*(double)g_data.nel) / (double)g_data.lltimeMemCpyDtoH_GPU;

	//      123456789-123456789-123456789-123456789-123456789-123456789
	printf("                                CPU Time             GPU Time\n");
	printf("MemCopy (H->D)....: %20I64i %20I64i\n",
		g_data.lltimeMemCpyHtoD_CPU, g_data.lltimeMemCpyHtoD_GPU);
	printf("MemCopy (D->H)....: %20I64i %20I64i\n",
		g_data.lltimeMemCpyDtoH_CPU, g_data.lltimeMemCpyDtoH_GPU);
	printf("Encode Input......: %20I64i                    -\n",
		g_data.lltimeEncode_CPU);
	printf("Decode Output.....: %20I64i                    -\n",
		g_data.lltimeDecode_CPU);
	printf("Kernel............: %20I64i %20I64i\n",
		g_data.lltimeKernel_CPU, g_data.lltimeKernel_GPU);
	printf("Total Time........: %20I64i %20I64i\n", 
		g_data.lltimeTotal_CPU, g_data.lltimeTotal_GPU);
	printf("GFLOP/s...........:                    - %20.4f\n",
		dGFLOPs);
	//printf("GB/s (H->D).......: %20.4f %20.4f\n",
	//	dGBs_HtoD_CPU, dGBs_HtoD_GPU);
	//*/
}