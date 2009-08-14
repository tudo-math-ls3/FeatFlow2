.version 1.4
#if defined(__DOUBLE_PRECISION__)
.target sm_13
#else
.target sm_10
#endif //  defined(__DOUBLE_PRECISION__)

/* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 * =-> PTX Streamline-Diffusion Kernel for the 2D Rannacher-Turek Element  <-=
 * -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 * =-> Copyright (C) 2009 by Peter Zajac                                   <-=
 * -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 * 
 * The kernel in this file assembles the following differential operator:
 *
 *              alpha * M + beta * L + gamma * N + delta * S
 *
 * where:
 * M is the mass matrix:              M := (          phi ,      psi  )
 * L is the stiffness matrix:         L := (     grad(phi), grad(psi) )
 * N is the gradient tensor matrix:   N := ( u * grad(phi),      psi  )
 * S is the streamline diffusion stabilisation matrix:
 *          S := sum_T delta_T * ( u * grad(phi), u * grad(psi) )_T
 *
 *
 * This file needs to be preprocessed using a C preprocessor before it can
 * be loaded into a CUDA module.
 */

#if defined(__DOUBLE_PRECISION__)
	#define _FP				f64
	#define DATASIZE		8					// 8 bytes per float
	#define THREADSIZE		2048
	#define BLOCKSIZE		65536

	#define _ZERO_			0d0000000000000000	// = 0.0
	#define _ONE_			0d3ff0000000000000	// = 1.0
	#define _TWO_			0d4000000000000000	// = 2.0
	#define _HALF_			0d3fe0000000000000	// = 0.5
	#define _QUARTER_		0d3fd0000000000000	// = 0.25
	#define _EPS_			0d3e45798ee2308c3a	// = 1e-8

	#define _PT0_			0dbfe8c97ef43f7248	// = -sqrt(2/3)
	#define _PT1_			0d0000000000000000	// = 0
	#define _PT2_			0d3fe8c97ef43f7248	// =  sqrt(2/3)
	#define _OM0_			0d3fd3c0ca4587e6b7	// = 25/81
	#define _OM1_			0d3fdf9add3c0ca458	// = 40/81
	#define _OM2_			0d3fe948b0fcd6e9e0	// = 64/81

	#define _PS0_			0d3fe93cd3a2c8198e	// = (1 + sqrt(1/3))/2
	#define _PS1_			0d3fcb0cb174df99c7	// = (1 - sqrt(1/3))/2


	#define FADD			add.rn.f64
	#define FSUB			sub.rn.f64
	#define FNEG			neg.f64
	#define FMUL			mul.rn.f64
	#define FMAD			mad.rn.f64
	#define FABS			abs.f64
	#define FRCP			rcp.rn.f64
	#define FDIV			div.rn.f64
	#define FSQRT			sqrt.rn.f64
	#define FMOV			mov.f64
	#define FLDL			ld.local.f64
	#define FLDG			ld.global.f64
	#define FLDC			ld.const.f64
	#define FSTL			st.local.f64
	#define FSTG			st.global.f64
	#define FSETP(op)		setp.op.f64

#else
	#define _FP				f32
	#define DATASIZE		4					// 4 bytes per float
	#define THREADSIZE		1024
	#define BLOCKSIZE		32768

	#define _ZERO_			0f00000000			// = 0.0
	#define _ONE_			0f3f800000			// = 1.0
	#define _TWO_			0f40000000			// = 2.0
	#define _HALF_			0f3f000000			// = 0.5
	#define _QUARTER_		0f3e800000			// = 0.25
	#define _EPS_			0f3a83126f			// = 1e-3

	#define _PT0_			0fbf464bf8			// = -sqrt(2/3)
	#define _PT1_			0f00000000			// = 0
	#define _PT2_			0f3f464bf8			// =  sqrt(2/3)
	#define _OM0_			0f3e9e0652			// = 25/81
	#define _OM1_			0f3efcd6ea			// = 40/81
	#define _OM2_			0f3f4a4588			// = 64/81

	#define _PS0_			0f3f49e69d			// = (1 + sqrt(1/3))/2
	#define _PS1_			0f3e58658e			// = (1 - sqrt(1/3))/2

	#define FADD			add.f32
	#define FSUB			sub.f32
	#define FNEG			neg.f32
	#define FMUL			mul.f32
	#define FMAD			mad.f32
	#define FABS			abs.f32
	#define FRCP			rcp.approx.f32
	#define FDIV			div.approx.f32
	#define FSQRT			sqrt.approx.f32
	#define FMOV			mov.f32
	#define FLDL			ld.local.f32
	#define FLDG			ld.global.f32
	#define FLDC			ld.const.f32
	#define FSTL			st.local.f32
	#define FSTG			st.global.f32
	#define FSETP(op)		setp.op.f32
#endif // defined(__DOUBLE_PRECISION__)

// weighting factor for mass matrix
.const ._FP g_alpha = _ZERO_;

// weighting factor for laplace matrix, pre-multiplied by nu
.const ._FP g_beta  = _ZERO_;

// weighting factor for convection
.const ._FP g_gamma = _ONE_;

// weighting factor for stabilisation
.const ._FP g_delta = _ZERO_;

// maximum norm of velocity field coefficient vector
.const ._FP g_maxu  = _ZERO_;

// viscosity parameter nu; must be > 0 if delta != 0 !!!
.const ._FP g_nu    = _ONE_;

// assembly kernel entrypoint
.entry asm_kernel(.param .u32 pbuffer)
{
	.reg .pred p;
	.reg .u32 %ri<2>;		// 32 bit integer registers
	.reg ._FP %rf<12>;		// floating point registers
	.local ._FP %lf<64>;	// floating point local memory

	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	// Step 1: initialisation mumbo-jumbo
	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	// calculate address of this thread's block in the io buffer
	ld.param.u32		%ri0, [pbuffer];
	//cvt.u32.u16			%ri1, %ctaid.x;
	//mad.lo.u32			%ri0, %ri1, BLOCKSIZE, %ri0;
	cvt.u32.u16			%ri1, %tid.x;
	shr.u32				%ri1, %ri1, 4;
	mad.lo.u32			%ri0, %ri1, THREADSIZE, %ri0;
	cvt.u32.u16			%ri1, %tid.x;
	and.b32				%ri1, %ri1, 0x0000000f;
	mad.lo.u32			%ri0, %ri1, DATASIZE, %ri0;

#if defined(__DOUBLE_PRECISION__)
	// define vertex position addresses
	#define VX0			[%ri0+   0]
	#define VY0			[%ri0+ 128]
	#define VX1			[%ri0+ 256]
	#define VY1			[%ri0+ 384]
	#define VX2			[%ri0+ 512]
	#define VY2			[%ri0+ 640]
	#define VX3			[%ri0+ 768]
	#define VY3			[%ri0+ 896]

	// define DOF value addresses
	#define D10			[%ri0+1024]
	#define D20			[%ri0+1152]
	#define D11			[%ri0+1280]
	#define D21			[%ri0+1408]
	#define D12			[%ri0+1536]
	#define D22			[%ri0+1664]
	#define D13			[%ri0+1792]
	#define D23			[%ri0+1920]
#else
	// define vertex position addresses
	#define VX0			[%ri0+  0]
	#define VY0			[%ri0+ 64] 
	#define VX1			[%ri0+128]
	#define VY1			[%ri0+192]
	#define VX2			[%ri0+256]
	#define VY2			[%ri0+320]
	#define VX3			[%ri0+384]
	#define VY3			[%ri0+448]

	// define DOF value addresses
	#define D10			[%ri0+512]
	#define D20			[%ri0+576]
	#define D11			[%ri0+640]
	#define D21			[%ri0+704]
	#define D12			[%ri0+768]
	#define D22			[%ri0+832]
	#define D13			[%ri0+896]
	#define D23			[%ri0+960]
#endif // defined(__DOUBLE_PRECISION__)

	// define output matrix addresses (overlap with input addresses)
	#define B00			VX0
	#define B01			VY0
	#define B02			VX1
	#define B03			VY1
	#define B10			VX2
	#define B11			VY2
	#define B12			VX3
	#define B13			VY3
	#define B20			D10
	#define B21			D20
	#define B22			D11
	#define B23			D21
	#define B30			D12
	#define B31			D22
	#define B32			D13
	#define B33			D23

	// define transformation coefficient variables
	#define S00			[%lf0]
	#define S01			[%lf1]
	#define S10			[%lf2]
	#define S11			[%lf3]
	#define T0			[%lf4]
	#define T1			[%lf5]
	#define R0			[%lf6]
	#define R1			[%lf7]

	// define integration weight variables
	#define OM00		[%lf8]
	#define OM01		[%lf9]
	#define OM02		[%lf10]
	#define OM10		[%lf11]
	#define OM11		[%lf12]
	#define OM12		[%lf13]
	#define OM20		[%lf14]
	#define OM21		[%lf15]
	#define OM22		[%lf16]

	// define temporary matrix
	#define A00			[%lf17]
	#define A01			[%lf18]
	#define A02			[%lf19]
	#define A03			[%lf20]
	#define A10			[%lf21]
	#define A11			[%lf22]
	#define A12			[%lf23]
	#define A13			[%lf24]
	#define A20			[%lf25]
	#define A21			[%lf26]
	#define A22			[%lf27]
	#define A23			[%lf28]
	#define A30			[%lf29]
	#define A31			[%lf30]
	#define A32			[%lf31]
	#define A33			[%lf32]

	// define coefficient matrix
	#define C00			[%lf43]
	#define C01			[%lf44]
	#define C02			[%lf45]
	#define C03			[%lf46]
	#define C10			[%lf47]
	#define C11			[%lf48]
	#define C12			[%lf49]
	#define C13			[%lf50]
	#define C20			[%lf51]
	#define C21			[%lf52]
	#define C22			[%lf53]
	#define C23			[%lf54]
	#define C30			[%lf55]
	#define C31			[%lf56]
	#define C32			[%lf57]
	#define C33			[%lf58]

	// local stabilisation parameter variable
	#define LOCDELTA	[%lf59]

	// define dof value variable (overlaps with A)
	#define DV10		A00
	#define DV11		A01
	#define DV12		A02
	#define DV13		A03
	#define DV20		A10
	#define DV21		A11
	#define DV22		A12
	#define DV23		A13

	// define derivative variables (overlaps with A)
	#define DX0			A20
	#define DX1			A21
	#define DX2			A22
	#define DX3			A23
	#define DY0			A30
	#define DY1			A31
	#define DY2			A32
	#define DY3			A33

	// define function value variables
	#define	FV0			[%lf60]
	#define	FV1			[%lf61]
	#define	FV2			[%lf62]
	#define	FV3			[%lf63]


	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	// Step 2: calculate transformation coefficients
	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	// define temporary register names
	#define Z0			%rf0
	#define Z1			%rf1
	#define Z2			%rf2
	#define Z3			%rf3
	#define NL0			%rf4
	#define NL1			%rf5
	#define ZS00		%rf6
	#define ZS01		%rf7
	#define ZS10		%rf8
	#define ZS11		%rf9
	#define ZR0			%rf10
	#define ZR1			%rf11

	// load VX* into Z*
	FLDG		Z0, VX0;
	FLDG		Z1, VX1;
	FLDG		Z2, VX2;
	FLDG		Z3, VX3;
	// calculate S00
	FNEG		ZS00, Z0;
	FADD		ZS00, ZS00, Z1;
	FADD		ZS00, ZS00, Z2;
	FSUB		ZS00, ZS00, Z3;
	FMUL		ZS00, ZS00, _QUARTER_;
	// calculate S01
	FNEG		ZS01, Z0;
	FSUB		ZS01, ZS01, Z1;
	FADD		ZS01, ZS01, Z2;
	FADD		ZS01, ZS01, Z3;
	FMUL		ZS01, ZS01, _QUARTER_;
	// calculate quadrilateral midpoint
	FMOV		ZR0, Z0;
	FADD		ZR0, ZR0, Z1;
	FADD		ZR0, ZR0, Z2;
	FADD		ZR0, ZR0, Z3;
	FMUL		ZR0, ZR0, _QUARTER_;
	FSTL		R0, ZR0;
	// calculate NL0
	FMOV		NL0, Z0;
	FSUB		NL0, NL0, Z1;
	FADD		NL0, NL0, Z2;
	FSUB		NL0, NL0, Z3;
	FMUL		NL0, NL0, _QUARTER_;
	// load VY* into Z*
	FLDG		Z0, VY0;
	FLDG		Z1, VY1;
	FLDG		Z2, VY2;
	FLDG		Z3, VY3;
	// calculate S10
	FNEG		ZS10, Z0;
	FADD		ZS10, ZS10, Z1;
	FADD		ZS10, ZS10, Z2;
	FSUB		ZS10, ZS10, Z3;
	FMUL		ZS10, ZS10, _QUARTER_;
	// calculate S11
	FNEG		ZS11, Z0;
	FSUB		ZS11, ZS11, Z1;
	FADD		ZS11, ZS11, Z2;
	FADD		ZS11, ZS11, Z3;
	FMUL		ZS11, ZS11, _QUARTER_;
	// calculate quadrilateral midpoint
	FMOV		ZR1, Z0;
	FADD		ZR1, ZR1, Z1;
	FADD		ZR1, ZR1, Z2;
	FADD		ZR1, ZR1, Z3;
	FMUL		ZR1, ZR1, _QUARTER_;
	FSTL		R1, ZR1;
	// calculate NL1
	FMOV		NL1, Z0;
	FSUB		NL1, NL1, Z1;
	FADD		NL1, NL1, Z2;
	FSUB		NL1, NL1, Z3;
	FMUL		NL1, NL1, _QUARTER_;

	// calculate integration weights
	#define __DETJ(z,x,y,o)	\
		FMAD	Z0, NL0, y, ZS00;\
		FMAD	Z1, NL1, x, ZS11;\
		FMAD	Z2, NL0, x, ZS01;\
		FMAD	Z3, NL1, y, ZS10;\
		FMUL	Z0, Z0, Z1;\
		FMUL	Z2, Z2, Z3;\
		FSUB	Z0, Z0, Z2;\
		FABS	Z0, Z0;\
		FMUL	Z0, Z0, o;\
		FSTL	z, Z0;

	__DETJ(OM00,_PT0_,_PT0_,_OM0_)
	__DETJ(OM10,_PT1_,_PT0_,_OM1_)
	__DETJ(OM20,_PT2_,_PT0_,_OM0_)
	__DETJ(OM01,_PT0_,_PT1_,_OM1_)
	__DETJ(OM11,_PT1_,_PT1_,_OM2_)
	__DETJ(OM21,_PT2_,_PT1_,_OM1_)
	__DETJ(OM02,_PT0_,_PT2_,_OM0_)
	__DETJ(OM12,_PT1_,_PT2_,_OM1_)
	__DETJ(OM22,_PT2_,_PT2_,_OM0_)

	#undef __DETJ
	
	// calculate inverse determinant of S
	FMUL		Z0, ZS00, ZS11;
	FMUL		Z1, ZS01, ZS10;
	FSUB		Z0, Z0, Z1;
	FRCP		Z0, Z0;
	FNEG		Z1, Z0;
	
	// invert S
	FMUL		ZS00, ZS00, Z0;
	FMUL		ZS01, ZS01, Z1;
	FMUL		ZS10, ZS10, Z1;
	FMUL		ZS11, ZS11, Z0;

	// calculate t := S * NL
	FMUL		Z0, ZS00, NL0;
	FMAD		Z0, ZS01, NL1, Z0;
	FMUL		Z1, ZS10, NL0;
	FMAD		Z1, ZS11, NL1, Z1;
	
	// store S
	FSTL		S00, ZS00;
	FSTL		S01, ZS01;
	FSTL		S10, ZS10;
	FSTL		S11, ZS11;

	// store t
	FSTL		T0, Z0;	
	FSTL		T1, Z1;

	// undefine temporary register names
	#undef Z0
	#undef Z1
	#undef Z2
	#undef Z3
	#undef NL0
	#undef NL1
	#undef ZS00
	#undef ZS01
	#undef ZS10
	#undef ZS11
	#undef ZT0
	#undef ZT1

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// write r,s,t and return
	FLDL		%rf0, R0;
	FSTG		B00, %rf0;
	FLDL		%rf0, R1;
	FSTG		B01, %rf0;
	FLDL		%rf0, S00;
	FSTG		B02, %rf0;
	FLDL		%rf0, S01;
	FSTG		B03, %rf0;
	FLDL		%rf0, S10;
	FSTG		B10, %rf0;
	FLDL		%rf0, S11;
	FSTG		B11, %rf0;
	FLDL		%rf0, T0;
	FSTG		B12, %rf0;
	FLDL		%rf0, T1;
	FSTG		B13, %rf0;
	ret.uni;
	//*/

	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// write omega and return
	FLDL		%rf0, OM00;
	FSTG		B00, %rf0;
	FLDL		%rf0, OM01;
	FSTG		B01, %rf0;
	FLDL		%rf0, OM02;
	FSTG		B02, %rf0;
	FLDL		%rf0, OM10;
	FSTG		B10, %rf0;
	FLDL		%rf0, OM11;
	FSTG		B11, %rf0;
	FLDL		%rf0, OM12;
	FSTG		B12, %rf0;
	FLDL		%rf0, OM20;
	FSTG		B20, %rf0;
	FLDL		%rf0, OM21;
	FSTG		B21, %rf0;
	FLDL		%rf0, OM22;
	FSTG		B22, %rf0;
	ret.uni;
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	// Step 3: assemble rows of coefficient matrix
	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	// define temporary register names
	#define X0			%rf0
	#define X1			%rf1
	#define Y0			%rf2
	#define Y1			%rf3
	#define P0			%rf4
	#define P1			%rf5
	#define C0			%rf6
	#define C1			%rf7
	#define ZS			%rf8

	#define __ACMR(i,j) \
		FLDG	C0, VX##i;\
		FLDG	C1, VX##j;\
		FLDL	P0, R0;\
		FNEG	P0, P0;\
		FMOV	P1, P0;\
		FMAD	P0, _PS0_, C0, P0;\
		FMAD	P0, _PS1_, C1, P0;\
		FMAD	P1, _PS1_, C0, P1;\
		FMAD	P1, _PS0_, C1, P1;\
		FLDL	ZS, S00;\
		FMUL	X0, ZS, P0;\
		FMUL	X1, ZS, P1;\
		FLDL	ZS, S10;\
		FMUL	Y0, ZS, P0;\
		FMUL	Y1, ZS, P1;\
		FLDG	C0, VY##i;\
		FLDG	C1, VY##j;\
		FLDL	P0, R1;\
		FNEG	P0, P0;\
		FMOV	P1, P0;\
		FMAD	P0, _PS0_, C0, P0;\
		FMAD	P0, _PS1_, C1, P0;\
		FMAD	P1, _PS1_, C0, P1;\
		FMAD	P1, _PS0_, C1, P1;\
		FLDL	ZS, S01;\
		FMAD	X0, ZS, P0, X0;\
		FMAD	X1, ZS, P1, X1;\
		FLDL	ZS, S11;\
		FMAD	Y0, ZS, P0, Y0;\
		FMAD	Y1, ZS, P1, Y1;\
		FMOV	ZS, _ONE_;\
		FSTL	A##i##0, ZS;\
		FADD	C0, X0, X1;\
		FMUL	C0, _HALF_, C0;\
		FSTL	A##i##1, C0;\
		FADD	C0, Y0, Y1;\
		FMUL	C0, _HALF_, C0;\
		FSTL	A##i##2, C0;\
		FMUL	C0, Y0, Y0;\
		FMAD	C0, Y1, Y1, C0;\
		FNEG	C0, C0;\
		FMAD	C0, X0, X0, C0;\
		FMAD	C0, X1, X1, C0;\
		FMUL	C0, _HALF_, C0;\
		FSTL	A##i##3, C0;

	__ACMR(0,1)
	__ACMR(1,2)
	__ACMR(2,3)
	__ACMR(3,0)

	#undef __ACMR

	// undefine temporary register names
	#undef X0
	#undef X1
	#undef Y0
	#undef Y1
	#undef P0
	#undef P1
	#undef C0
	#undef C1
	#undef ZS

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// transpose A
	FLDL		%rf0, A01;
	FLDL		%rf1, A10;
	FSTL		A01, %rf1;
	FSTL		A10, %rf0;
	FLDL		%rf0, A02;
	FLDL		%rf1, A20;
	FSTL		A02, %rf1;
	FSTL		A20, %rf0;
	FLDL		%rf0, A03;
	FLDL		%rf1, A30;
	FSTL		A03, %rf1;
	FSTL		A30, %rf0;
	FLDL		%rf0, A12;
	FLDL		%rf1, A21;
	FSTL		A12, %rf1;
	FSTL		A21, %rf0;
	FLDL		%rf0, A13;
	FLDL		%rf1, A31;
	FSTL		A13, %rf1;
	FSTL		A31, %rf0;
	FLDL		%rf0, A23;
	FLDL		%rf1, A32;
	FSTL		A23, %rf1;
	FSTL		A32, %rf0;
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// write out A and return
	FLDL		%rf0, A00;
	FSTG		B00, %rf0;
	FLDL		%rf0, A01;
	FSTG		B01, %rf0;
	FLDL		%rf0, A02;
	FSTG		B02, %rf0;
	FLDL		%rf0, A03;
	FSTG		B03, %rf0;
	FLDL		%rf0, A10;
	FSTG		B10, %rf0;
	FLDL		%rf0, A11;
	FSTG		B11, %rf0;
	FLDL		%rf0, A12;
	FSTG		B12, %rf0;
	FLDL		%rf0, A13;
	FSTG		B13, %rf0;
	FLDL		%rf0, A20;
	FSTG		B20, %rf0;
	FLDL		%rf0, A21;
	FSTG		B21, %rf0;
	FLDL		%rf0, A22;
	FSTG		B22, %rf0;
	FLDL		%rf0, A23;
	FSTG		B23, %rf0;
	FLDL		%rf0, A30;
	FSTG		B30, %rf0;
	FLDL		%rf0, A31;
	FSTG		B31, %rf0;
	FLDL		%rf0, A32;
	FSTG		B32, %rf0;
	FLDL		%rf0, A33;
	FSTG		B33, %rf0;
	ret.uni;
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	// Step 4: Invert element coefficient matrix
	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	#define W0			%rf0
	#define W1			%rf1
	#define W2			%rf2
	#define W3			%rf3
	#define W4			%rf4
	#define W5			%rf5
	#define WD			%rf6
	#define Z0			%rf7
	#define Z1			%rf8
	#define Z2			%rf9
	#define Z3			%rf10
	#define ZA			%rf11

	#define __CALC_WORK(i,j) \
		FLDL	Z0, A##i##0;\
		FLDL	Z1, A##i##1;\
		FLDL	Z2, A##i##2;\
		FLDL	Z3, A##i##3;\
		FLDL	ZA, A##j##0;\
		FMUL	W0, ZA, Z1;\
		FNEG	W0, W0;\
		FMUL	W1, ZA, Z2;\
		FNEG	W1, W1;\
		FMUL	W2, ZA, Z3;\
		FNEG	W2, W2;\
		FLDL	ZA, A##j##1;\
		FMUL	W3, ZA, Z2;\
		FNEG	W3, W3;\
		FMUL	W4, ZA, Z3;\
		FNEG	W4, W4;\
		FMAD	W0, ZA, Z0, W0;\
		FLDL	ZA, A##j##2;\
		FMUL	W5, ZA, Z3;\
		FNEG	W5, W5;\
		FMAD	W1, ZA, Z0, W1;\
		FMAD	W3, ZA, Z1, W3;\
		FLDL	ZA, A##j##3;\
		FMAD	W2, ZA, Z0, W2;\
		FMAD	W4, ZA, Z1, W4;\
		FMAD	W5, ZA, Z2, W5;

	#define __ELIMINATE_ROW(j) \
		FLDL	ZA, A##j##0;\
		FMUL	Z3, ZA, W3;\
		FNEG	Z3, Z3;\
		FMUL	Z2, ZA, W4;\
		FNEG	Z2, Z2;\
		FMUL	Z1, ZA, W5;\
		FNEG	Z1, Z1;\
		FLDL	ZA, A##j##1;\
		FMUL	Z0, ZA, W5;\
		FNEG	Z0, Z0;\
		FMAD	Z3, ZA, W1, Z3;\
		FNEG	Z3, Z3;\
		FMAD	Z2, ZA, W2, Z2;\
		FNEG	Z2, Z2;\
		FLDL	ZA, A##j##2;\
		FMAD	Z1, ZA, W2, Z1;\
		FNEG	Z1, Z1;\
		FMAD	Z0, ZA, W4, Z0;\
		FNEG	Z0, Z0;\
		FMAD	Z3, ZA, W0, Z3;\
		FLDL	ZA, A##j##3;\
		FMAD	Z2, ZA, W0, Z2;\
		FMAD	Z1, ZA, W1, Z1;\
		FMAD	Z0, ZA, W3, Z0;

	#define __STORE_COL_02(i) \
		FMUL	Z0, Z0, WD;\
		FMUL	Z1, Z1, WD;\
		FMUL	Z2, Z2, WD;\
		FMUL	Z3, Z3, WD;\
		FNEG	Z1, Z1;\
		FNEG	Z3, Z3;\
		FSTL	C0##i, Z0;\
		FSTL	C1##i, Z1;\
		FSTL	C2##i, Z2;\
		FSTL	C3##i, Z3;

	#define __STORE_COL_13(i) \
		FMUL	Z0, Z0, WD;\
		FMUL	Z1, Z1, WD;\
		FMUL	Z2, Z2, WD;\
		FMUL	Z3, Z3, WD;\
		FNEG	Z0, Z0;\
		FNEG	Z2, Z2;\
		FSTL	C0##i, Z0;\
		FSTL	C1##i, Z1;\
		FSTL	C2##i, Z2;\
		FSTL	C3##i, Z3;


	__CALC_WORK(2,3)
	__ELIMINATE_ROW(1)	// => column 0

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// write out W and return
	FSTG		B00, W0;
	FSTG		B01, W1;
	FSTG		B02, W2;
	FSTG		B03, W3;
	FSTG		B10, W4;
	FSTG		B11, W5;
	ret.uni;
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	// calculate determinant of A
	FLDL		ZA, A00;
	FMUL		WD, ZA, Z0;
	FNEG		WD, WD;
	FLDL		ZA, A01;
	FMAD		WD, ZA, Z1, WD;
	FNEG		WD, WD;
	FLDL		ZA, A02;
	FMAD		WD, ZA, Z2, WD;
	FNEG		WD, WD;
	FLDL		ZA, A03;
	FMAD		WD, ZA, Z3, WD;
	FNEG		WD, WD;
	FRCP		WD, WD;

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// write out determinant
	FSTG		B00, WD;
	ret.uni;
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	__STORE_COL_02(0)
	__ELIMINATE_ROW(0)	// => column 1
	__STORE_COL_13(1)
	__CALC_WORK(0,1)
	__ELIMINATE_ROW(3)	// => column 2
	__STORE_COL_02(2)
	__ELIMINATE_ROW(2)	// => column 3
	__STORE_COL_13(3)

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// write out W and return
	FSTG		B00, W0;
	FSTG		B01, W1;
	FSTG		B02, W2;
	FSTG		B03, W3;
	FSTG		B10, W4;
	FSTG		B11, W5;
	ret.uni;
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	// That's it - matrix inverted (hopefully) !
	#undef __CALC_WORK
	#undef __CALC_COL
	#undef __STORE_COL_02
	#undef __STORE_COL_13

	#undef W0
	#undef W1
	#undef W2
	#undef W3
	#undef W4
	#undef W5
	#undef WD
	#undef Z0
	#undef Z1
	#undef Z2
	#undef Z3
	#undef ZA

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// write out C and return
	FLDL		%rf0, C00;
	FSTG		B00, %rf0;
	FLDL		%rf0, C01;
	FSTG		B01, %rf0;
	FLDL		%rf0, C02;
	FSTG		B02, %rf0;
	FLDL		%rf0, C03;
	FSTG		B03, %rf0;
	FLDL		%rf0, C10;
	FSTG		B10, %rf0;
	FLDL		%rf0, C11;
	FSTG		B11, %rf0;
	FLDL		%rf0, C12;
	FSTG		B12, %rf0;
	FLDL		%rf0, C13;
	FSTG		B13, %rf0;
	FLDL		%rf0, C20;
	FSTG		B20, %rf0;
	FLDL		%rf0, C21;
	FSTG		B21, %rf0;
	FLDL		%rf0, C22;
	FSTG		B22, %rf0;
	FLDL		%rf0, C23;
	FSTG		B23, %rf0;
	FLDL		%rf0, C30;
	FSTG		B30, %rf0;
	FLDL		%rf0, C31;
	FSTG		B31, %rf0;
	FLDL		%rf0, C32;
	FSTG		B32, %rf0;
	FLDL		%rf0, C33;
	FSTG		B33, %rf0;
	ret.uni;
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	// Step 5: Calculate local delta
	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	#define U1			%rf0
	#define U2			%rf1
	#define LOCALH		%rf2
	#define UNORM		%rf3
	#define Z0			%rf4
	#define Z1			%rf5
	#define X0			%rf6
	#define Y0			%rf7
	#define X1			%rf8
	#define Y1			%rf9
	#define RELOC		%rf10

	// We can skip the calculation of the local delta if one of the following
	// conditions is true:
	// 1.) g_maxu == 0
	// 2.) g_delta == 0

	FLDC		Z0, [g_maxu];
	FLDC		Z1, [g_delta];
	FMUL		Z1, Z0, Z1;
	FSETP(eq)	p, Z1, _ZERO_;
	FMOV		Z0, _ZERO_;		// set local delta to 0
@p	bra.uni		L1;

	// calculate mean velocity on this element
	FLDG		U1, D10;
	FLDG		Z0, D11;
	FADD		U1, U1, Z0;
	FLDG		Z0, D12;
	FADD		U1, U1, Z0;
	FLDG		Z0, D13;
	FADD		U1, U1, Z0;
	FLDG		U2, D20;
	FLDG		Z0, D21;
	FADD		U2, U2, Z0;
	FLDG		Z0, D22;
	FADD		U2, U2, Z0;
	FLDG		Z0, D23;
	FADD		U2, U2, Z0;
	FMUL		Z0, U1, U1;
	FMAD		Z0, U2, U2, Z0;
	FSQRT		UNORM, Z0;
	FMUL		UNORM, _QUARTER_, UNORM;

	// check UNORM against EPS now

	FMOV		Z0, _ZERO_;		// set local delta to 0
	FSETP(lt)	p, UNORM, _EPS_;
@p	bra			L1;

	// calculate local H
	FLDG		X0, VX0;
	FLDG		X1, VX2;
	FLDG		Y0, VY1;
	FLDG		Y1, VY3;
	FSUB		X1, X1, X0;
	FSUB		Y1, Y1, Y0;
	FMUL		Z0, X1, Y1;
	FLDG		X0, VX1;
	FLDG		X1, VX3;
	FLDG		Y0, VY0;
	FLDG		Y1, VY2;
	FSUB		X1, X1, X0;
	FSUB		Y1, Y1, Y0;
	FMUL		Z1, X1, Y1;
	FSUB		Z0, Z0, Z1;
	FMUL		Z0, Z0, _HALF_;
	FSQRT		LOCALH, Z0;

	// calculate reloc
	FMUL		RELOC, UNORM, LOCALH;
	FLDC		Z0, [g_nu];
	FDIV		RELOC, RELOC, Z0;

	// calculate local delta
	FLDC		Z0, [g_delta];
	FLDC		Z1, [g_maxu];
	FMUL		Z0, Z0, Z1;
	FMUL		Z0, Z0, LOCALH;
	FMUL		Z0, Z0, _TWO_;
	FMUL		Z0, Z0, RELOC;
	FADD		RELOC, RELOC, _ONE_;
	FDIV		Z0, Z0, RELOC;

L1:
	// store delta
	FSTL		LOCDELTA, Z0;

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// write out local H and delta
	FSTG		B00, Z0;
	FSTG		B01, LOCALH;
	FSTG		B02, UNORM;
	ret.uni;
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	#undef U1
	#undef U2
	#undef LOCALH
	#undef UNORM
	#undef Z0
	#undef Z1
	#undef X0
	#undef Y0
	#undef X1
	#undef Y1
	#undef RELOC

	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	// Step 6: assemble output matrix
	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

	// copy dof values from global to local memory
	FLDG		%rf0, D10;
	FSTL		DV10, %rf0;
	FLDG		%rf0, D11;
	FSTL		DV11, %rf0;
	FLDG		%rf0, D12;
	FSTL		DV12, %rf0;
	FLDG		%rf0, D13;
	FSTL		DV13, %rf0;
	FLDG		%rf0, D20;
	FSTL		DV20, %rf0;
	FLDG		%rf0, D21;
	FSTL		DV21, %rf0;
	FLDG		%rf0, D22;
	FSTL		DV22, %rf0;
	FLDG		%rf0, D23;
	FSTL		DV23, %rf0;

	// format the output matrix
	FMOV		%rf0, _ZERO_;
	FSTG		B00, %rf0;
	FSTG		B01, %rf0;
	FSTG		B02, %rf0;
	FSTG		B03, %rf0;
	FSTG		B10, %rf0;
	FSTG		B11, %rf0;
	FSTG		B12, %rf0;
	FSTG		B13, %rf0;
	FSTG		B20, %rf0;
	FSTG		B21, %rf0;
	FSTG		B22, %rf0;
	FSTG		B23, %rf0;
	FSTG		B30, %rf0;
	FSTG		B31, %rf0;
	FSTG		B32, %rf0;
	FSTG		B33, %rf0;

	#define U1			%rf0
	#define U2			%rf1
	#define X			%rf2
	#define Y			%rf3
	#define ZS00		%rf4
	#define ZS01		%rf5
	#define ZS10		%rf6
	#define ZS11		%rf7
	#define Z0			%rf8
	#define Z1			%rf9
	#define Z2			%rf10
	#define Z3			%rf11
	#define OMEGA		X
	#define	RES			Y
	#define ALPHA		ZS00
	#define BETA		ZS01
	#define GAMMA		ZS10
	#define DELTA		ZS11

	// load S into registers
	#define __LOAD_LIN_TRAFO \
		FLDL	ZS00, S00;\
		FLDL	ZS01, S01;\
		FLDL	ZS10, S10;\
		FLDL	ZS11, S11;

	#define __BILIN_DISTO(x,y) \
		FLDL	X, T0;\
		FMAD	X, X, y, _ONE_;\
		FMUL	X, X, x;\
		FLDL	Y, T1;\
		FMAD	Y, Y, x, _ONE_;\
		FMUL	Y, Y, y;

	#define __EVAL_BASIS(i) \
		FLDL	Z0, C0##i;\
		FLDL	Z1, C1##i;\
		FLDL	Z2, C2##i;\
		FLDL	Z3, C3##i;\
		FMAD	U1, X, Z3, Z1;\
		FMUL	U2, Y, Z3;\
		FSUB	U2, Z2, U2;\
		FMAD	Z0, X, U1, Z0;\
		FMAD	Z0, Y, U2, Z0;\
		FSTL	FV##i, Z0;\
		FMUL	U1, X, Z3;\
		FMAD	U1, U1, _TWO_, Z1;\
		FMUL	U2, Y, Z3;\
		FNEG	U2, U2;\
		FMAD	U2, U2, _TWO_, Z2;\
		FMUL	Z2, ZS00, U1;\
		FMAD	Z2, ZS10, U2, Z2;\
		FMUL	Z3, ZS01, U1;\
		FMAD	Z3, ZS11, U2, Z3;\
		FSTL	DX##i, Z2;\
		FSTL	DY##i, Z3;

	#define __CALC_VELO \
		FLDL	Z0, FV0;\
		FLDL	Z2, DV10;\
		FLDL	Z3, DV20;\
		FMUL	U1, Z0, Z2;\
		FMUL	U2, Z0, Z3;\
		FLDL	Z0, FV1;\
		FLDL	Z2, DV11;\
		FLDL	Z3, DV21;\
		FMAD	U1, Z0, Z2, U1;\
		FMAD	U2, Z0, Z3, U2;\
		FLDL	Z0, FV2;\
		FLDL	Z2, DV12;\
		FLDL	Z3, DV22;\
		FMAD	U1, Z0, Z2, U1;\
		FMAD	U2, Z0, Z3, U2;\
		FLDL	Z0, FV3;\
		FLDL	Z2, DV13;\
		FLDL	Z3, DV23;\
		FMAD	U1, Z0, Z2, U1;\
		FMAD	U2, Z0, Z3, U2;

	#define __LOAD_SCALARS(o) \
		FLDC	ALPHA, [g_alpha];\
		FLDC	BETA, [g_beta];\
		FLDC	GAMMA, [g_gamma];\
		FLDL	DELTA, LOCDELTA;\
		FLDL	OMEGA, o;

	#define __ASM_TRIAL(i,j) \
		FMOV	RES, _ZERO_;\
		FLDL	Z0, DX##j;\
		FLDL	Z1, DY##j;\
		FMUL	Z2, Z0, U1;\
		FMAD	Z2, Z1, U2, Z2;\
		FMAD	RES, Z2, Z3, RES;\
		FLDL	Z2, DX##i;\
		FMUL	Z0, Z0, Z2;\
		FLDL	Z2, DY##i;\
		FMAD	Z0, Z1, Z2, Z0;\
		FMAD	RES, Z0, BETA, RES;\
		FLDL	Z0, FV##i;\
		FLDL	Z1, FV##j;\
		FMUL	Z0, Z0, Z1;\
		FMAD	RES, Z0, ALPHA, RES;\
		FLDG	Z0, B##i##j;\
		FMAD	RES, RES, OMEGA, Z0;\
		FSTG	B##i##j, RES;

	#define __ASM_TEST(i) \
		FLDL	Z0, DX##i;\
		FLDL	Z1, DY##i;\
		FMUL	Z3, U1, Z0;\
		FMAD	Z3, U2, Z1, Z3;\
		FMUL	Z3, Z3, DELTA;\
		FLDL	Z0, FV##i;\
		FMAD	Z3, Z0, GAMMA, Z3;\
		__ASM_TRIAL(i,0)\
		__ASM_TRIAL(i,1)\
		__ASM_TRIAL(i,2)\
		__ASM_TRIAL(i,3)

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// set up a mass matrix
	#undef __ASM_TRIAL
	#define __ASM_TRIAL(i,j) \
		FLDL	Z1, FV##j;\
		FLDG	RES, B##i##j;\
		FMAD	RES, Z1, Z0, RES;\
		FSTG	B##i##j, RES;
	
	#undef __ASM_TEST
	#define __ASM_TEST(i) \
		FLDL	Z0, FV##i;\
		FMUL	Z0, Z0, OMEGA;\
		__ASM_TRIAL(i,0)\
		__ASM_TRIAL(i,1)\
		__ASM_TRIAL(i,2)\
		__ASM_TRIAL(i,3)
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// set up a laplace matrix
	#undef __ASM_TRIAL
	#define __ASM_TRIAL(i,j) \
		FLDL	Z1, DX##j;\
		FLDL	Z3, DY##j;\
		FLDG	RES, B##i##j;\
		FMAD	RES, Z1, Z0, RES;\
		FMAD	RES, Z3, Z2, RES;\
		FSTG	B##i##j, RES;
	
	#undef __ASM_TEST
	#define __ASM_TEST(i) \
		FLDL	Z0, DX##i;\
		FLDL	Z2, DY##i;\
		FMUL	Z0, Z0, OMEGA;\
		FMUL	Z2, Z2, OMEGA;\
		__ASM_TRIAL(i,0)\
		__ASM_TRIAL(i,1)\
		__ASM_TRIAL(i,2)\
		__ASM_TRIAL(i,3)
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	//>>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>> DEBUG >>>
	/*
	// evaluate & write function values and exit
	__LOAD_LIN_TRAFO
	__BILIN_DISTO(_ZERO_,_ZERO_)
	__EVAL_BASIS(0)
	FSTG		B00, Z0;
	FSTG		B10, Z2;
	FSTG		B20, Z3;
	__EVAL_BASIS(1)
	FSTG		B01, Z0;
	FSTG		B11, Z2;
	FSTG		B21, Z3;
	__EVAL_BASIS(2)
	FSTG		B02, Z0;
	FSTG		B12, Z2;
	FSTG		B22, Z3;
	__EVAL_BASIS(3)
	FSTG		B03, Z0;
	FSTG		B13, Z2;
	FSTG		B23, Z3;
	ret.uni;
	//*/
	//<<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<< DEBUG <<<

	#define __ASM_ENTRIES(x,y,o) \
		__LOAD_LIN_TRAFO\
		__BILIN_DISTO(x,y)\
		__EVAL_BASIS(0)\
		__EVAL_BASIS(1)\
		__EVAL_BASIS(2)\
		__EVAL_BASIS(3)\
		__CALC_VELO\
		__LOAD_SCALARS(o)\
		__ASM_TEST(0)\
		__ASM_TEST(1)\
		__ASM_TEST(2)\
		__ASM_TEST(3)
	
	// okay, let's process all cubature points
	__ASM_ENTRIES(_PT0_, _PT0_, OM00)
	__ASM_ENTRIES(_PT0_, _PT1_, OM01)
	__ASM_ENTRIES(_PT0_, _PT2_, OM02)
	__ASM_ENTRIES(_PT1_, _PT0_, OM10)
	__ASM_ENTRIES(_PT1_, _PT1_, OM11)
	__ASM_ENTRIES(_PT1_, _PT2_, OM12)
	__ASM_ENTRIES(_PT2_, _PT0_, OM20)
	__ASM_ENTRIES(_PT2_, _PT1_, OM21)
	__ASM_ENTRIES(_PT2_, _PT2_, OM22)

} // .entry asm_kernel()
