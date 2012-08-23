#pragma once

#if defined(__cplusplus)
extern "C" {
#endif

// Error codes
#define GPUSD_INFO_SUCCESS		0	// everything's fine
#define GPUSD_INFO_NOT_READY	1	// necessary data missing
#define GPUSD_INFO_CUDA_FAIL	2	// failed to init CUDA
#define GPUSD_INFO_NO_DEVICE	3	// no CUDA device found
#define GPUSD_INFO_MEM_ALLOC	4	// memory allocation failure

/**
 * \brief Initializes the GPU-SD driver.
 *
 * \param[out] info
 * One of the GPUSD_INFO_* error codes
 */
void gpusd_init(int* info);

/**
 * \brief Releases the GPU-SD driver.
 */
void gpusd_done();

/**
 * \brief Sets the scalar parameters.
 *
 * \param[in] dnu
 * Viscosity parameter nu = 1/RE.
 *
 * \param[in] dalpha
 * Scaling factor for the mass matrix,
 *
 * \param[in] dbeta
 * Scaling factor for the laplace matrix, pre-multiplied by nu.
 *
 * \param[in] dgamma
 * Scaling factor for the convection.
 *
 * \param[in] ddelta
 * Scaling factor for the stabilisation.
 */
void gpusd_setparams(double *dnu, double *dalpha, double *dbeta,
					 double *dgamma, double *ddelta);

/**
 * \brief Sets the arrays from the triangulation.
 *
 * \param[in] nel
 * Specifies the number of elements in the triangulation.
 *
 * \param[in] Dvtx
 * The DvertexCoords array from the FEAT2 triangulation.\n
 * Is silently assumed to be have the dimensions Dvtx[nel][2].
 *
 * \param[in] Iverts
 * The IverticesAtElement array from the FEAT2 triangulation.
 * Is silently assumed to be have the dimensions Iverts[nel][4].
 *
 * \param[in] Iedges
 * The IedgesAtElement array from the FEAT2 triangulation.
 * Is silently assumed to be have the dimensions Iedges[nel][4].
 */
void gpusd_settria(int *nel, double *Dvtx, int *Iverts, int *Iedges);

/**
 * \brief Sets the arrays from the matrix.
 *
 * \param[in] neq
 * The number of equations (=rows) of the matrix.
 *
 * \param[in] rowptr
 * The row-pointer (Kld) array of the matrix.
 *
 * \param[in] colidx
 * The column-index (Kcol) array of the matrix.
 *
 * \param[in] data
 * The data (DA) array of the matrix.
 */
void gpusd_setmatrix(int *neq, int *Irowptr, int *Icolidx, double *Ddata);

/**
 * \brief Sets the velocity field coefficient vectors.
 *
 * \param[in] Du1, Du2
 * The velocity vector field coefficient vectors.
 */
void gpusd_setvelo(double *Du1, double *Du2);

/**
 *
 */
void gpusd_assemble(int *info);

/**
 * \brief Prints the timing results of the last assembly onto the standard
 * output stream
 */
void gpusd_statistics();

#if defined(__cplusplus)
}
#endif
