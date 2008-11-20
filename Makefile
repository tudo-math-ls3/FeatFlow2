#!/usr/bin/env gmake

FEATFLOW=.
include $(FEATFLOW)/Globals.mk

MAKEFLAGS += --no-print-directory
###MAKEFLAGS += -s #for silent build

intro:
	@echo "________________________________________________________________________"
	@echo "                               FEATFLOW $(FFVER)                        "
	@echo "________________________________________________________________________"
	@echo 
	@echo "To check wether your system is recognized by the installation, type:"
	@echo 
	@echo " make id           (see Globals.mk in the case of unrecognized system)"
	@echo 
	@echo "To install Featflow (=compile everything), type:"
	@echo 
	@echo " make build        (may require 10-60 mins, 50MB disk storage per ID)"
	@echo 
	@echo "To run a short test calculation for every application, type:"
	@echo 
	@echo " make test         (may require 1-10 mins, 128MB memory)"
	@echo 
	@echo "To run the Featflow benchmark calculations, type:"
	@echo 
	@echo " make bench        (may require 20-120 mins, 256MB memory)"
	@echo 
	@echo "Other options are:"
	@echo 
	@echo " make sysinfo   -  print detailed information about current system"
	@echo 
	@echo " make all       -  compile and test everything, run benchmark"
	@echo " make libs      -  compile the libraries only"
	@echo " make apps      -  compile the applications (+libraries)"
	@echo 
	@echo " make test      -  test the applications (makes libs apps)"
	@echo " make bench     -  compile and run the benchmark test"
	@echo "                   (compiles the needed libraries and apps)"
	@echo 
	@echo " make help      -  for additional help and further options"
	@echo 
	@echo " make clean_all -  clean the installation for current architecture"
	@echo " make purge_all -  clean the installation for all architectures"
	@echo

finish:
	@echo
	@echo "________________________________________________________________________"
	@echo "         Installation finished"
	@echo "________________________________________________________________________"
	@echo
	@echo "To run a short test calculation for every application, type:"
	@echo 
	@echo " make test         (may require 1-10 mins, 128MB memory)"
	@echo 
	@echo "To run the Featflow benchmark calculations, type:"
	@echo 
	@echo " make bench        (may require 20-120 mins, 256MB memory)"
	@echo 

	@echo "To start an example application:"
	@echo " => Switch to the application directory (applications/cc2d, pp2d,...)"
	@echo "    and start the executable (cc2d, pp2d,...)"
	@echo ""
	@echo "To adapt the application to your own needs:"
	@echo " => Duplicate the directory of the application (cc2d, pp2d,...)"
	@echo "    inside of the \"applications\"-directory"
	@echo "    (or anywhere outside of the Featflow directory - in this case"
	@echo "     modify the Makefile such that the variable FEATFLOW points to:"
	@echo "     FEATFLOW=" `pwd` ")"
	@echo " => Modify the files to your own needs (indatXd.f, parqXd.f, etc.)"
	@echo " => Type \"make\" in that directory to compile the application"
	@echo " => Modify the input data file (cc2d.dat, pp2d.dat,...) inside of"
	@echo "    the #data subdirectory to your own needs"
	@echo " => Start the executable"
	@echo

.NOTPARALLEL: greet id intro finish test bench
build: greet id libs apps finish
all: greet intro id libs apps test bench finish

greet: 
	@echo "Featflow $(FFVER) starting installation ..."
	@date

#FEATLIB+= $(if $(BLASLIB),,lapack blas)
apps: libs
	@date
	@$(foreach i, $(APPS), (cd applications/$i && $(MAKE) all ); )
	@date

libs: 
	@date
	@mkdir -p $(FEATFLOW)/object/libraries
	@$(foreach i, $(BUILDLIB), (cd libraries/$i && $(MAKE) all ); )
	@date

docs: html ps pdf

html:
	@date
	@$(foreach i, $(APPS), (cd applications/$i && $(MAKE) html ); )
	@date

ps:
	@date
	@$(foreach i, $(APPS), (cd applications/$i && $(MAKE) ps ); )
	@date

pdf:
	@date
	@$(foreach i, $(APPS), (cd applications/$i && $(MAKE) pdf ); )
	@date

bench: 
	@echo "Running Featflow benchmarks..."
	@date
	@(cd applications/benchmark && $(MAKE) bench)
	@date

report: 
	@(if [ ! -f "$(BENCHLOG)" ] ; then echo 'Benchmark protocol file' $(BENCHLOG) ; echo 'for this ID:' $(ID) 'and host:' $(HOST); echo 'is not available!' ; echo 'Please run the benchmark computation (make bench) first.'; exit 1 ; fi)
	@echo "mail featflow@featflow.de <$(BENCHLOG)"

clean: clean_libs clean_apps
clean_apps:
	@$(foreach i, $(APPS), (echo "Cleaning up "$(i) ; cd applications/$i && $(MAKE) clean ); )
clean_libs: 
	@$(foreach i, $(BUILDLIB), (echo "Cleaning up "$(i) ; cd libraries/$i && $(MAKE) clean ); )

purge: purge_libs purge_apps
purge_apps:
	@$(foreach i, $(APPS), (echo "Cleaning up "$(i) ; cd applications/$i && $(MAKE) purge ); )
purge_libs:
	@$(foreach i, $(LIBS), (echo "Cleaning up "$(i) ; cd libraries/$i && $(MAKE) purge ); )

purge_all:
	@$(foreach i, $(APPS), (echo "Cleaning up "$(i) ; cd applications/$i && $(MAKE) purge_all ); )
	@$(foreach i, $(LIBS), (echo "Cleaning up "$(i) ; cd libraries/$i && $(MAKE) purge_all ); )
	-rm -f -r object/libraries/lib-*

clean_all:
	@$(foreach i, $(APPS), (echo "Cleaning up "$(i) ; cd applications/$i && $(MAKE) clean_all ); )
	@$(foreach i, $(LIBS), (echo "Cleaning up "$(i) ; cd libraries/$i && $(MAKE) clean_all ); )

debug: debug_libs debug_apps
debug_apps:
	@$(foreach i, $(APPS), (echo "Cleaning up "$(i) ; cd applications/$i && $(MAKE) debug ); )
debug_libs:
	@$(foreach i, $(BUILDLIB), (echo "Cleaning up "$(i) ; cd libraries/$i && $(MAKE) debug ); )

loc: all
	(bin/loc)

$(LIBS):
	@(cd libraries/$@ && $(MAKE) all)

$(APPS):
	@(cd applications/$@ && $(MAKE) all)

test: libs
	@date
	@$(foreach i, $(APPS), (cd applications/$i && $(MAKE) test ); )
	@date

# "make dist" creates a TARball distribution file of the FEATFLOW
# package! The version number defined in Globals.mk is used to
# name the file properly.

dist: fixperms purge_all
	@bin/clear_tildefiles.sh
	(cd .. && tar --exclude .svn --exclude .rem --exclude area51 --exclude docs --exclude documentation --exclude incoming --exclude 'libgoto*' --exclude cc2doptcontrol --exclude Deps.mk -cvzf Featflow2_$(FFVER).tar.gz Featflow2)
fixperms:
	@(cd ./bin ; chmod ug+x guess_id f77mkdep.sh clear_tildefiles.sh info_cpu info_f77)

DONEIDS=$(shell (cd object/libraries 2> /dev/null && (ls | grep lib)))
DONELIBS= \
  $(if $(DONEIDS), \
    $(foreach i, $(DONEIDS), \
      echo ' ' $(i:lib-%=%): \
      $(patsubst lib%.a,%,$(shell (cd object/libraries/$i; ls -C lib*.a)));),\
    echo 'none')

id: .id

status: .id
	@echo 'Featflow version:' $(FFVER)
	@echo 
	@echo 'Featflow modules:' $(APPS)
	@echo 
	@echo 'Available libraries:' $(LIBS)
	@echo 
	@echo 'Libraries compiled:'
	@$(DONELIBS)
	@echo 

sysinfo:
	@echo 
	@echo "  Featflow system information page"
	@echo "------------------------------------"
	@(uname -a)
	@($(MAKE) id)
	-@($(FEATFLOW)/bin/info_cpu $(ID))
	-@($(FEATFLOW)/bin/info_f77 $(ID) $(FC))

help: .help
	@echo 
	@echo 'Available targets for the top level Makefile:'
	@echo ' all           - compile all libraries and application modules'
	@echo ' build         - the same as "make all", but designed for first installation'
	@echo ' apps          - compile all application modules'
	@echo ' libs          - compile all libraries'
	@echo ' debug_libs    - compile all libraries with debug infos'
	@echo ' debug_apps    - compile all application modules with debug infos'
	@echo 
	@echo ' test          - runs an application tests for all applications'
	@echo ' bench         - runs a featflow benchmark'
	@echo 
	@echo ' status        - print id and list already compiled libs'
	@echo 
	@echo ' clean         - remove all for current ID not needed for run (object files)'
	@echo ' clean_libs    - remove all object files of libraries'
	@echo ' clean_apps    - remove all object files of applications'
	@echo ' clean_all     - remove all object files of libs and apps'
	@echo 
	@echo ' purge         - remove all that can be removed for current ID'
	@echo ' purge_libs    - remove all libraries + object files for current ID'
	@echo ' purge_apps    - remove all applications + object files for current ID'
	@echo 
	@echo ' purge_all     - remove all that can be removed for all IDs'
	@echo 

.PHONY: docs