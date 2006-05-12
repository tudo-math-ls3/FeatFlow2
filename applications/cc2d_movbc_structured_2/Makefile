#!/usr/bin/env gmake

FEATFLOW=../..

include $(FEATFLOW)/Globals.mk

INCDIR= -I./src/geometry      -I./src/nonstationary    -I./src/postprocessing  \
	-I./src/triangulation -I./src/discr_linear     -I./src/gridadapt  \
	-I./src/stationary    -I./src/core             -I./src/init \
	-I./src/optimization  -I./src/solvers          -I./src/discr_nonlinear  \
	-I./src/misc	      -I./src/parametrization  -I./include

# User-specified files

SRC1= fictbdry.f  fictbdrypc.f  indat2d.f  parq2d.f

# Sorce files of core module

SRC2= optdlcalc.f  zvalue1.f core_util.f  rdstrt.f
	#cc2d_structured.f  

# Source files for low level discretisation

SRC3= arcalc.f  bmul.f  bndry2.f  em30.f  em31.f   xmrout.f  \
	ab79x.f ap79x.f mgrout10.f  mgrout30.f  mgrout31.f mattransp.f \
	ndfglx.f interpol.f va0x.f vb0x.f
     
# Files for high-level discretisation

SRC4= coeff.f      iniassembly.f  inivectors.f  rhs.f \
	gupwdx.f   inimatrices.f  matrest.f     supwdgx.f   \
	glmatasm.f matmul.f	  normvectors.f
     
# Source files for geometry support

SRC5= compgeometries.f  linedistance.f  pccompgeometries.f  qtrafo.f \
	ttrafo.f        geometries.f    ltrafo.f      pcgeometries.f \
	rdgeometry.f
     
# Source files for grid adaption module

SRC6= bdproj.f	  gsoutput.f       prgridadaptionstatic.f  \
	bndprj.f      rcgrid.f           xmon.f    \
	domeasures.f  mon.f      rdgridadaption.f   xsmmgr2.f \
	gaassembly.f  neumannfilter.f  umbrella2.f			     \
	gradient2.f   normmon.f  umf4sort.f        xsmvgr2.f     
     
# Source files for initialization process
     
SRC7= rdlinsol_cc2d.f rdparamtriang.f    \
	init2.f       rdnonlinsol_cc2d.f  rdpostprocessing.f \
	rddiscretization.f  rdoutput.f    rdtimediscr.f
	
# Miscellaneous source files

SRC8= datutil.f     dstrings2.f  dstringsaux.f  invert.f \
	condoutput.f  dstrings1.f  dstrings3.f    generalutil.f  \
	heapsort.f    random.f
	
# Nonstationary solver

SRC9= dfktinst2.f   ininonst2.f    mgstp2.f	 stdnexttimestep.f   \
	tstepping.f iniadtstep.f   ini_nonstationary.f  nonstead2.f  \
	fpost_nonstationary.f      stopcrit.f
	
# Stationary solver

SRC10= dfkt2.f	ini_stationary.f  nsdeflinsol.f  optcnl2.f  \
	ininsdef2.f nsdef2.f	  nsdefmgrout.f  fpost_stationary.f  \
	vanca2.f    vanca3.f      vanca4.f       vanca5.f 
	 
# Optimization modules

SRC11= optalgs.f  rdopt.f

# Parametrization

SRC12= inipar.f  parq2domega.f  parq2dwrap.f  rdparm.f  recbdry.f

# Postprocessing routines

SRC13= bdintegral.f  draglift3.f  error2.f  \
	tracking.f draglift2.f  \
	post2.f    aposteriorierror.f owm.f
     
# Extended solvers
     
SRC14= gsolvers.f  m020.f  xii41x.f  yj7.f
     
# Extended triangulation support

SRC15= srchquad.f      triutil2d.f \
	gmv2doutput.f  initria.f   scalareval.f  structria.f \
	quad2tri.f     tracers.f umf4_f77wrapper_dummy.f
	
SRCF90=

include $(FEATFLOW)/kernel/kernel.inc
	
SRCF90:=$(SRCF90) testcode.f90

# put all files together

SRC=$(SRC1) $(SRC2) $(SRC3) $(SRC4) $(SRC5) $(SRC6) $(SRC7) $(SRC8) \
	$(SRC9) $(SRC10) $(SRC11) $(SRC12) $(SRC13) $(SRC14) $(SRC15) 

# path for the make where to look for which files

vpath %.f src/geometry	  src/nonstationary src/postprocessing  \
	src/triangulation src/discr_linear  src/gridadapt  \
	src/stationary    src/core          src/init	 \
	src/optimization  src/solvers	    src/discr_nonlinear \
	src/misc src/parametrization \
	kernel/System kernel/BasicGeometry \
	kernel/LinearAlgebra kernel/ElementCubature \
	kernel/Triangulation kernel/ContinuousFormulation \
	kernel/Boundary kernel/DOFMaintainance \
	kernel/LinearSystem kernel/LinearSolver \
	kernel/ProblemSupport kernel

vpath %.f90 \
	kernel/System kernel/BasicGeometry \
	kernel/LinearAlgebra kernel/ElementCubature \
	kernel/Triangulation kernel/ContinuousFormulation \
	kernel/Boundary kernel/DOFMaintainance \
	kernel/LinearSystem kernel/LinearSolver \
	kernel/ProblemSupport kernel

vpath %.for src/geometry	  src/nonstationary src/postprocessing  \
	src/triangulation src/discr_linear  src/gridadapt  \
	src/stationary    src/core          src/init	 \
	src/optimization  src/solvers	    src/discr_nonlinear \
	src/misc src/parametrization \
	kernel/System kernel/BasicGeometry \
	kernel/LinearAlgebra kernel/ElementCubature \
	kernel/Triangulation kernel/ContinuousFormulation \
	kernel/Boundary kernel/DOFMaintainance \
	kernel/LinearSystem kernel/LinearSolver \
	kernel/ProblemSupport kernel

vpath %.inc src/geometry	  src/nonstationary src/postprocessing  \
	src/triangulation src/discr_linear  src/gridadapt  \
	src/stationary    src/core          src/init	 \
	src/optimization  src/solvers	    src/discr_nonlinear \
	src/misc src/parametrization \
	include \
	kernel/System kernel/BasicGeometry \
	kernel/LinearAlgebra kernel/ElementCubature \
	kernel/Triangulation kernel/ContinuousFormulation \
	kernel/Boundary kernel/DOFMaintainance \
	kernel/LinearSystem kernel/LinearSolver \
	kernel/ProblemSupport kernel

EXEC=cc2d_movbc-$(ID)

FEATLIB= feat2d sysutils #umfpack4 amd

include $(FEATFLOW)/Rules_apps_f90.mk
-include Deps.mk

test: clean_exec run

