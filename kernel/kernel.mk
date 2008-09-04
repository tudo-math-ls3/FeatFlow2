# This makefile is included by the application makefile. It adds
# the kernel source files to the SRCF90 variable. The application
# has to add missing .F90 files to it and pass it to the
# compiler.

KERNEL=$(FEATFLOW)/kernel
INCOMING=$(FEATFLOW)/incoming

INCDIR:= $(INCDIR) -I$(KERNEL)/System -I$(KERNEL)/Triangulation \
	-I$(KERNEL)/DOFMaintenance -I$(KERNEL)/Adaptivity

KERNELSRC:=fsystem.f90 basicgeometry.f90 geometryaux.f90 afcutil.f90 \
	genoutput.f90 linearalgebra.f90 storage.f90 derivatives.f90 \
	externalstorage.f90 cubature.f90 transformation.f90 triangulation.f90\
	scalarpde.f90 error.f90 element.f90 io.f90 boundary.f90 \
	fictitiousboundary.f90 boundarycondition.f90\
	paramlist.f90 spatialdiscretisation.f90 dofmapping.f90 \
	discretebc.f90 discretefbc.f90 linearsystemscalar.f90 \
	linearsystemblock.f90 globalsystem.f90 \
	matrixmodification.f90 matrixfilters.f90 sortstrategy.f90 \
	spdiscprojection.f90 multilevelprojection.f90 \
	domainintegration.f90 coarsegridcorrection.f90 feevaluation.f90\
	vectorfilters.f90 filtersupport.f90 vanka.f90 linearsolver.f90 \
	linearsolverautoinitialise.f90 collection.f90 \
	bilinearformevaluation.f90 linearformevaluation.f90 \
	bcassembly.f90 matrixio.f90 vectorio.f90 convection.f90 \
	nonlinearsolver.f90 matrixrestriction.f90 fparser.f90 \
	stack.f90 pprocnavierstokes.f90 ucd.f90 signal.f90 signal_ccode.c \
	timestepping.f90 pprocerror.f90 trilinearformevaluation.f90 \
	stdoperators.f90 sort.f90 quadtree.f90 l2projection.f90 \
	pprocgradients.f90 geometry.f90 mprimitives.f90 \
	hadaptivity.f90 list.f90 binarytree.f90 arraylist.f90 \
	graph.f90 octree.f90 triasearch.f90 statistics.f90 \
	meshmodification.f90 mapleio.f90 jumpstabilisation.f90 \
	boundaryintegral.f90 afcstabilisation.f90 groupfemscalar.f90 \
	hadaptaux.f90 hadaptaux2d.f90 hadaptaux3d.f90 groupfemsystem.f90 \
	meshregion.f90 pprocsolution.f90 elementpreprocessing.f90 iluk.f90 \
	gmvwritef.c multileveloperators.f90 adjacency.f90 meshadjacency.f90 \
	vanka_aux.f90

# Include HDF5 subsystem if required
ifeq ($(HDF5),YES)
KERNELSRC:=$(KERNELSRC) h5lite.f90 linearsystemh5io.f90
endif

# path for the make where to look for which files

vpath %.f $(KERNEL)/System $(KERNEL)/BasicGeometry \
	$(KERNEL)/LinearAlgebra $(KERNEL)/ElementCubature \
	$(KERNEL)/Triangulation $(KERNEL)/ContinuousFormulation \
	$(KERNEL)/Boundary $(KERNEL)/DOFMaintenance \
	$(KERNEL)/LinearSystem $(KERNEL)/LinearSolver \
	$(KERNEL)/ProblemSupport $(KERNEL)/Postprocessing \
	$(KERNEL)/NonlinearSolver $(KERNEL)/Projection \
	$(KERNEL)/Mathematics $(KERNEL)/PDEOperators \
	$(KERNEL)/TimeDependence $(KERNEL)/Adaptivity \
	$(KERNEL)/DataStructures $(KERNEL)/Preprocessing \
	$(KERNEL) $(INCOMING)

vpath %.f90 $(KERNEL)/System $(KERNEL)/BasicGeometry \
	$(KERNEL)/LinearAlgebra $(KERNEL)/ElementCubature \
	$(KERNEL)/Triangulation $(KERNEL)/ContinuousFormulation \
	$(KERNEL)/Boundary $(KERNEL)/DOFMaintenance \
	$(KERNEL)/LinearSystem $(KERNEL)/LinearSolver \
	$(KERNEL)/ProblemSupport $(KERNEL)/Postprocessing \
	$(KERNEL)/NonlinearSolver $(KERNEL)/Projection \
	$(KERNEL)/Mathematics $(KERNEL)/PDEOperators \
	$(KERNEL)/TimeDependence $(KERNEL)/Adaptivity \
	$(KERNEL)/DataStructures $(KERNEL)/Preprocessing \
	$(KERNEL) $(INCOMING)

vpath %.c $(KERNEL)/System $(KERNEL)/BasicGeometry \
	$(KERNEL)/LinearAlgebra $(KERNEL)/ElementCubature \
	$(KERNEL)/Triangulation $(KERNEL)/ContinuousFormulation \
	$(KERNEL)/Boundary $(KERNEL)/DOFMaintenance \
	$(KERNEL)/LinearSystem $(KERNEL)/LinearSolver \
	$(KERNEL)/ProblemSupport $(KERNEL)/Postprocessing \
	$(KERNEL)/NonlinearSolver $(KERNEL)/Projection \
	$(KERNEL)/Mathematics $(KERNEL)/PDEOperators \
	$(KERNEL)/TimeDependence $(KERNEL)/Adaptivity \
	$(KERNEL)/DataStructures $(KERNEL)/Preprocessing \
	$(KERNEL) $(INCOMING)

vpath %.inc $(KERNEL)/System $(KERNEL)/BasicGeometry \
	$(KERNEL)/LinearAlgebra $(KERNEL)/ElementCubature \
	$(KERNEL)/Triangulation $(KERNEL)/ContinuousFormulation \
	$(KERNEL)/Boundary $(KERNEL)/DOFMaintenance \
	$(KERNEL)/LinearSystem $(KERNEL)/LinearSolver \
	$(KERNEL)/ProblemSupport $(KERNEL)/Postprocessing \
	$(KERNEL)/NonlinearSolver $(KERNEL)/Projection \
	$(KERNEL)/Mathematics $(KERNEL)/PDEOperators \
	$(KERNEL)/TimeDependence $(KERNEL)/Adaptivity \
	$(KERNEL)/DataStructures $(KERNEL)/Preprocessing \
	$(KERNEL) $(INCOMING)

