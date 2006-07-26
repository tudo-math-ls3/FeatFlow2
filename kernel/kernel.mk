# This makefile is included by the application makefile. It adds
# the kernel source files to the SRCF90 variable. The application
# has to add missing .F90 files to it and pass it to the
# compiler.

KERNEL=$(FEATFLOW)/kernel
INCOMING=$(FEATFLOW)/incoming

KERNELSRC:=fsystem.f90 basicgeometry.f90 geometryaux.f90 afcutil.f90  \
	mprimitives.f90\
	genoutput.f90 linearalgebra.f90 storage.f90 derivatives.f90 \
	cubature.f90 transformation.f90 triangulation.f90 scalarpde.f90\
	error.f90 element.f90 io.f90 boundary.f90 \
	fictitiousboundary.f90 boundarycondition.f90\
	paramlist.f90 spatialdiscretisation.f90 dofmapping.f90 \
	discretebc.f90 discretefbc.f90 linearsystemscalar.f90 \
	linearsystemblock.f90 globalsystem.f90 \
	matrixmodification.f90 matrixfilters.f90 sortstrategy.f90 \
	spdiscprojection.f90 multilevelprojection.f90 \
	domainintegration.f90 coarsegridcorrection.f90 feevaluation.f90\
	vectorfilters.f90 filtersupport.f90 vanca.f90 linearsolver.f90 \
	linearsolverautoinitialise.f90 collection.f90 \
	bilinearformevaluation.f90 linearformevaluation.f90 \
	bcassembly.f90 matrixio.f90 vectorio.f90 convection.f90 \
	nonlinearsolver.f90 matrixrestriction.f90

# path for the make where to look for which files

vpath %.f $(KERNEL)/System $(KERNEL)/BasicGeometry \
	$(KERNEL)/LinearAlgebra $(KERNEL)/ElementCubature \
	$(KERNEL)/Triangulation $(KERNEL)/ContinuousFormulation \
	$(KERNEL)/Boundary $(KERNEL)/DOFMaintainance \
	$(KERNEL)/LinearSystem $(KERNEL)/LinearSolver \
	$(KERNEL)/ProblemSupport $(KERNEL)/Postprocessing \
	$(KERNEL)/NonlinearSolver $(KERNEL)/Projection\
	$(KERNEL)/Mathmatics $(KERNEL)/PDEOperators\
	$(KERNEL) $(INCOMING)

vpath %.f90 $(KERNEL)/System $(KERNEL)/BasicGeometry \
	$(KERNEL)/LinearAlgebra $(KERNEL)/ElementCubature \
	$(KERNEL)/Triangulation $(KERNEL)/ContinuousFormulation \
	$(KERNEL)/Boundary $(KERNEL)/DOFMaintainance \
	$(KERNEL)/LinearSystem $(KERNEL)/LinearSolver \
	$(KERNEL)/ProblemSupport $(KERNEL)/Postprocessing \
	$(KERNEL)/NonlinearSolver $(KERNEL)/Projection\
	$(KERNEL)/Mathmatics $(KERNEL)/PDEOperators\
	$(KERNEL) $(INCOMING)

vpath %.inc $(KERNEL)/System $(KERNEL)/BasicGeometry \
	$(KERNEL)/LinearAlgebra $(KERNEL)/ElementCubature \
	$(KERNEL)/Triangulation $(KERNEL)/ContinuousFormulation \
	$(KERNEL)/Boundary $(KERNEL)/DOFMaintainance \
	$(KERNEL)/LinearSystem $(KERNEL)/LinearSolver \
	$(KERNEL)/ProblemSupport $(KERNEL)/Postprocessing \
	$(KERNEL)/NonlinearSolver $(KERNEL)/Projection\
	$(KERNEL)/Mathmatics $(KERNEL)/PDEOperators\
	$(KERNEL) $(INCOMING)

