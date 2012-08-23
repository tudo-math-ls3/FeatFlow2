
// This is just for the MSVC compiler, so that it stops complaining
#if defined(_SHUT_UP)
	#define __inline__
	#define __device__
	#define __global__
	#define __constant__
	struct {int x,y,z;} threadIdx;
	#define unroll
	#define sqrtf(x) x
	#define sqrt(x) x
#endif

// single or double precision?
#if defined(__DOUBLE_PRECISION__)
	// define the real type
	typedef double		real;

	// some constants
	#define ZERO		0.0
	#define QUARTER		0.25
	#define HALF		0.5
	#define ONE			1.0
	#define TWO			2.0
	#define EPS			1e-8f

	// coordinates and weights for 3-point Gauss rule
	#define PT0		   -0.7745966692414834
	#define PT1			0.0
	#define PT2			0.7745966692414834
	#define OM0			0.5555555555555556
	#define OM1			0.8888888888888888
	#define OM2			0.5555555555555556

	// coordinates of 2-point Gauss rule
	#define PS0			0.7886751345948129
	#define PS1			0.2113248654051871

	// square root
	#define SQRT(x)		sqrt(x)
#else
	// define the real type
	typedef float		real;

	// some constants
	#define ZERO		0.0f
	#define QUARTER		0.25f
	#define HALF		0.5f
	#define ONE			1.0f
	#define TWO			2.0f
	#define EPS			1e-4f

	// coordinates and weights for 3-point Gauss rule
	#define PT0		   -0.7745967f
	#define PT1			0.0f
	#define PT2			0.7745967f
	#define OM0			0.5555556f
	#define OM1			0.8888888f
	#define OM2			0.5555556f

	// coordinates of 2-point Gauss rule
	#define PS0			0.7886751f
	#define PS1			0.2113249f

	// square root
	#define SQRT(x)		sqrtf(x)
#endif // !defined(__DOUBLE_PRECISION__)

// constant: maximal velocity norm
__constant__ real g_maxu = ONE;

// constant: viscosity parameter nu = 1/RE
__constant__ real g_nu = ONE;

// constant: mass matrix scaling factor
__constant__ real g_alpha = ZERO;

// constant: laplace matrix scaling factor (pre-multiplied by nu)
__constant__ real g_beta = ZERO;

// constant: convection scaling factor
__constant__ real g_gamma = ONE;

// constant: stabilisation scaling factor
__constant__ real g_delta = ONE;


/**
 * \brief Calculates transformation coefficients and integration weights.
 *
 * \param[in] v
 * The coordinates of the corner vertices.
 *
 * \param[out] r
 * Recieves the coordinates of the quadrilateral midpoint.
 *
 * \param[out] s
 * Recieves the 2x2 inverse linear transformation matrix.
 *
 * \param[out] t
 * Revieces the product of \p s and the non-linear coefficients of the
 * bilinear transformation.
 *
 * \param[out] omega
 * Recieves the integration weights in the cubature points.
 */
__inline__ __device__ void calcTrafo(real *v, real *r, real *s, real *t, real *omega)
{
	real dets;

	// assemble linear transformation, with s_00 and s_11 exchanged
	s[0] = QUARTER * (-v[0] + v[2] + v[4] - v[6]);
	s[1] = QUARTER * (-v[0] - v[2] + v[4] + v[6]);
	s[2] = QUARTER * (-v[1] + v[3] + v[5] - v[7]);
	s[3] = QUARTER * (-v[1] - v[3] + v[5] + v[7]);

	// assemble non-linear coefficients of the bilinear trafo
	r[0] = QUARTER * (v[0] - v[2] + v[4] - v[6]);
	r[1] = QUARTER * (v[1] - v[3] + v[5] - v[7]);

	// our bilinear transformation (modulo translation) is now given as:
	// ( x ) |--> ( s_11  s_01 ) * ( x ) + ( r_0 ) * ( x*y )
	// ( y )      ( s_10  s_00 )   ( y )   ( r_1 )

	// calculate the integration weights
	#define DETJ(x,y)	((s[0]+r[0]*y)*(s[3]+r[1]*x)\
						-(s[1]+r[0]*x)*(s[2]+r[1]*y))
	omega[0] = DETJ(PT0,PT0)*OM0*OM0;
	omega[1] = DETJ(PT1,PT0)*OM1*OM0;
	omega[2] = DETJ(PT2,PT0)*OM2*OM0;
	omega[3] = DETJ(PT0,PT1)*OM0*OM1;
	omega[4] = DETJ(PT1,PT1)*OM1*OM1;
	omega[5] = DETJ(PT2,PT1)*OM2*OM1;
	omega[6] = DETJ(PT0,PT2)*OM0*OM2;
	omega[7] = DETJ(PT1,PT2)*OM1*OM2;
	omega[8] = DETJ(PT2,PT2)*OM2*OM2;
	#undef DETJ

	// calculate determinant of S
	dets = (s[0]*s[3] - s[1]*s[2]);

	// calculate the inverse of S
	// Note: Now you know why we stored S with s_00 and s_11 exchanged - assuming
	//       that you know the direct inversion formula for 2x2 matrices ^_^
	s[0] *=  dets;
	s[1] *= -dets;
	s[2] *= -dets;
	s[3] *=  dets;

	// calculate t := S^-1 * r
	t[0] = s[0]*r[0] + s[1]*r[1];
	t[1] = s[2]*r[0] + s[3]*r[1];

	// calculate quadrilateral midpoint
	r[0] = QUARTER * (v[0] + v[2] + v[4] + v[6]);
	r[1] = QUARTER * (v[1] + v[3] + v[5] + v[7]);
}

/**
 * \brief Calculates the stabilisation parameter delta.
 *
 * \param[in] v
 * The coordinates of the corner vertices.
 *
 * \param[in] d
 * The values of the DOFs.
 *
 * \param[out] delta
 * The stabilisation parameter delta.
 */
__inline__ __device__ void calcDelta(real *v, real *d, real *delta)
{
	real u[2], unorm, reloc, localH;

	// calculate mean velocity vector
	u[0] = d[0] + d[2] + d[4] + d[6];
	u[1] = d[1] + d[3] + d[5] + d[7];

	// calculate norm of mean velocity vector
	unorm = SQRT(u[0]*u[0] + u[1]*u[1]);

	// calculate delta
	if(unorm > EPS)
	{
		// calculate local H
		// todo: replace calculation of local H by "ray" version later
		localH = QUARTER * SQRT((v[4]-v[0])*(v[7]-v[3])
							   -(v[6]-v[2])*(v[5]-v[1]));

		// Samarskji style upwind
		//reloc = (unorm*localH) / g_nu;
		//*delta = (g_delta*g_maxu*localH*TWO*reloc)/(ONE+reloc);
		*delta = ZERO;
	}
	else
		*delta = ZERO;
}

/**
 * \brief Inverts a 4x4 matrix.
 *
 * \param[in] A
 * The 4x4 matrix that is to be inverted.
 *
 * \param[out] B
 * Recieves the inverse of \p A.
 */
__inline__ __device__ void invMatrix(real *A, real *B)
{
	real W[7];
	W[0]=A[8]*A[13]-A[9]*A[12];W[1]=A[8]*A[14]-A[10]*A[12];
	W[2]=A[8]*A[15]-A[11]*A[12];W[3]=A[9]*A[14]-A[10]*A[13];
	W[4]=A[9]*A[15]-A[11]*A[13];W[5]=A[10]*A[15]-A[11]*A[14];
	B[0] = A[5]*W[5]-A[6]*W[4]+A[7]*W[3];
	B[4] =-A[4]*W[5]+A[6]*W[2]-A[7]*W[1];
	B[8] = A[4]*W[4]-A[5]*W[2]+A[7]*W[0];
	B[12]=-A[4]*W[3]+A[5]*W[1]-A[6]*W[0];
	W[6]=ONE/(A[0]*B[0]+A[1]*B[4]+A[2]*B[8]+A[3]*B[12]);
	B[0]*=W[6];B[4]*=W[6];B[8]*=W[6];B[12]*=W[6];
	B[1] =W[6]*(-A[1]*W[5]+A[2]*W[4]-A[3]*W[3]);
	B[5] =W[6]*( A[0]*W[5]-A[2]*W[2]+A[3]*W[1]);
	B[9] =W[6]*(-A[0]*W[4]+A[1]*W[2]-A[3]*W[0]);
	B[13]=W[6]*( A[0]*W[3]-A[1]*W[1]+A[2]*W[0]);
	W[0]=A[0]*A[5]-A[1]*A[4];W[1]=A[0]*A[6]-A[2]*A[4];
	W[2]=A[0]*A[7]-A[3]*A[4];W[3]=A[1]*A[6]-A[2]*A[5];
	W[4]=A[1]*A[7]-A[3]*A[5];W[5]=A[2]*A[7]-A[3]*A[6];
	B[2] =W[6]*( A[13]*W[5]-A[14]*W[4]+A[15]*W[3]);
	B[6] =W[6]*(-A[12]*W[5]+A[14]*W[2]-A[15]*W[1]);
	B[10]=W[6]*( A[12]*W[4]-A[13]*W[2]+A[15]*W[0]);
	B[14]=W[6]*(-A[12]*W[3]+A[13]*W[1]-A[14]*W[0]);
	B[3] =W[6]*(-A[9]*W[5]+A[10]*W[4]-A[11]*W[3]);
	B[7] =W[6]*( A[8]*W[5]-A[10]*W[2]+A[11]*W[1]);
	B[11]=W[6]*(-A[8]*W[4]+A[9]*W[2]-A[11]*W[0]);
	B[15]=W[6]*( A[8]*W[3]-A[9]*W[1]+A[10]*W[0]);
}

/**
 * \brief Assembles coefficent matrix row
 *
 * \param[in] r
 * The coordinates of the quadrilateral midpoint.
 *
 * \param[in] s
 * The 2x2 matrix of the inverse affine transformation.
 *
 * \param[in] v1
 * The coordinates of the start-point of the edge.
 *
 * \param[in] v2
 * The coordinates of the end-point of the edge.
 *
 * \param[out] A
 * The 1x4 row of the coefficient matrix that is to be assembled.
 */
__inline__ __device__ void aux_acmr(real *r, real *s, real *v1, real *v2, real *A)
{
	real p[2],x[2],y[2],omega;

	// calculate integration weight = 0.5 * edge_len
	p[0] = v2[0]-v1[0];
	p[1] = v2[1]-v1[1];
	omega = HALF * SQRT(p[0]*p[0] + p[1]*p[1]);

	// map point 1 onto quad
	p[0] = v1[0]*PS0 + v2[0]*PS1 - r[0];
	p[1] = v1[1]*PS0 + v2[1]*PS1 - r[1];
	x[0] = s[0]*p[0] + s[1]*p[1];
	y[0] = s[2]*p[0] + s[3]*p[1];

	// map point 2 onto quad
	p[0] = v1[0]*PS1 + v2[0]*PS0 - r[0];
	p[1] = v1[1]*PS1 + v2[1]*PS0 - r[1];
	x[1] = s[0]*p[0] + s[1]*p[1];
	y[1] = s[2]*p[0] + s[3]*p[1];

	// evaluate basis functions
	A[0] = omega*TWO;
	A[1] = omega*(x[0]+x[1]);
	A[2] = omega*(y[0]+y[1]);
	A[3] = omega*(x[0]*x[0] + x[1]*x[1] - y[0]*y[0] - y[1]*y[1]);
}

/**
 * \brief Evaluates the element and updates the output matrix.
 *
 * \param[in] s, t
 * The transformation coefficients, as returned by calcTrafo().
 *
 * \param[in] d
 * The dof values on the current element.
 *
 * \param[in] C
 * The inverted coefficient matrix of the element.
 *
 * \param[in] x,y
 * The coordinates of the cubature point where to evaluate.
 *
 * \param[in] omega
 * The integration weight for the cubature point.
 *
 * \param[in] delta
 * The stabilisation parameter.
 *
 * \param[in,out] A
 * The output matrix that is to be updated.
 */
__inline__ __device__ void aux_eval(real *s, real *t, real *d, real *C,
	real x, real y, real omega, real delta, real *A)
{
	int i,j;
	real rx,ry,u1,u2,z;
	real f[4];				// function values in point
	real dx[4], dy[4];		// derivates in point

	// todo: explain this...
	rx = x * (ONE + y*t[0]);
	ry = y * (ONE + x*t[1]);

	// evaluate basis functions
	#pragma unroll
	for(i = 0; i < 4; i++)
	{
		// evaluate function value
		f[i] = C[i] + rx*(C[4+i] + rx*C[12+i])
					+ ry*(C[8+i] - ry*C[12+i]);

		// evaluate 'reference' derivatives
		rx = C[4+i] + TWO*C[12+i]*rx;
		ry = C[8+i] - TWO*C[12+i]*ry;

		// evaluate 'real' derivatives
		dx[i] = s[0]*rx + s[1]*ry;
		dy[i] = s[2]*rx + s[3]*ry;
	}

	// calculate velocity field in current point
	u1 = u2 = ZERO;
	#pragma unroll
	for(i = 0; i < 4; i++)
	{
		u1 += f[i]*d[2*i+0];
		u2 += f[i]*d[2*i+1];
	}

	// loop over all test functions
	#pragma unroll
	for(i = 0; i < 4; i++)
	{
		// pre-calculate test function factor
		z = omega*(f[i] + delta*(u1*dx[i] + u2*dy[i]));

		// loop over all trial functions
		#pragma unroll
		for(j = 0; j < 4; j++)
			A[4*i+j] += (u1*dx[j] + u2*dy[j])*z;
	}
}

/**
 * \brief Streamline-Diffusion assembly kernel
 */
//__global__ void sd_kernel(real *vtx, real *dof, real *mat)
__global__ void sd_kernel(real *io)
{
	int i,j;
	real v[8];				// local vertex coordinates
	real d[8];				// local dof values
	real r[2],s[4],t[2];	// transformation coefficients
	real omega[9];			// integration weights
	real delta;				// stabilisation parameter delta
	real A[16];				// 4x4 temporary / output matrix
	real C[16];				// 4x4 coefficient matrix

	// fetch the thread index, multiplied by 16
	j = 16*threadIdx.x;

	// Step 1: fetch the vertex coordinates and dof values for this thread
	#pragma unroll
	for(i = 0; i < 8; i++)
		v[i] = io[j++];
	#pragma unroll
	for(i = 0; i < 8; i++)
		d[i] = io[j++];

	// Step 2: calculate transformation coefficients and integration weights
	calcTrafo(v, r, s, t, omega);

	// Step 3: calculate stabilisation parameter delta
	calcDelta(v, d, &delta);

	// Step 4: assemble rows of coefficient matrix
	aux_acmr(r, s, &v[0*2], &v[1*2], &A[0*4]);
	aux_acmr(r, s, &v[1*2], &v[2*2], &A[1*4]);
	aux_acmr(r, s, &v[2*2], &v[3*2], &A[2*4]);
	aux_acmr(r, s, &v[3*2], &v[0*2], &A[3*4]);

	// Step 5: invert coefficient matrix
	invMatrix(A,C);

	// Step 6: format output matrix
	#pragma unroll
	for(i = 0; i < 16; i++)
		A[i] = ZERO;

	// Step 7: assemble matrix
	aux_eval(s, t, d, C, PT0, PT0, OM0*OM0, delta, A);
	aux_eval(s, t, d, C, PT1, PT0, OM1*OM0, delta, A);
	aux_eval(s, t, d, C, PT2, PT0, OM2*OM0, delta, A);
	aux_eval(s, t, d, C, PT0, PT1, OM0*OM1, delta, A);
	aux_eval(s, t, d, C, PT1, PT1, OM1*OM1, delta, A);
	aux_eval(s, t, d, C, PT2, PT1, OM2*OM1, delta, A);
	aux_eval(s, t, d, C, PT0, PT2, OM0*OM2, delta, A);
	aux_eval(s, t, d, C, PT1, PT2, OM1*OM2, delta, A);
	aux_eval(s, t, d, C, PT2, PT2, OM2*OM2, delta, A);

	// Step 8: write matrix into output buffer
	j = 16*threadIdx.x;
	#pragma unroll
	for(i = 0; i < 16; i++)
		io[j++] = A[i];

	// That's it!
}
