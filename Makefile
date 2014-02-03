###################################################################
# Author:  Christian Jungreuthmayer                               #
# Date:    2. August 2002                                         #
# Licence: GPL                                                    #
###################################################################

###################################################################
# Attention: this Makefile is for gnumake                        #
# If you get an error message using it with 'make', then try to   #
# invoke 'gmake'                                                  #
###################################################################

###################################################################
# compiler flags                                                  #
###################################################################

ALTHOST =
FCOMP       = gfortran   ## ifort
#FFLAGS      = $(NAGFMOD) 
F_LNK_FLAGS = $(NAGFLIB) -lpthread

####################################################################
# site-independent files                                           #
####################################################################
FILES = simla3_0.o simla_subroutines.o qedroutines.o  

NAGFILES = simla2_3.o simla_subroutines.o qedroutines.o besselNAG.o egen_sub.o d2pdXdt.o  #invariants.o chigamma.o 



###################################################################
# FORTRAN compilation rule                                        #
###################################################################
.SUFFIXES: .f .o

# how to make a .o-file from a .f-file:
.f.o:
	$(FCOMP) $(FFLAGS) $(CORD) -c $*.f

.SUFFIXES: .f90 .o

# how to make a .o-file from a .f-file:
.f90.o:
	$(FCOMP) $(FFLAGS) $(CORD) -c $*.f90

####################################################################
# link chigamma
####################################################################
simla: include.i $(FILES) 
	$(FCOMP)  -o simla3 $(FILES) $(F_LNK_FLAGS)  $(LIBRARIES)

simlaNAG: include.i $(NAGFILES) 
	$(FCOMP)  -o simla2 $(NAGFILES) $(F_LNK_FLAGS)  $(LIBRARIES)




####################################################################

analyse: 
	$(FCOMP) $(FFLAGS_A) -c $*.f90
	include.i $(FILES)
	$(FCOMP)  -o simla2 $(F_LNK_FLAGS_A) $(FILES) $(LIBRARIES_A)
####################################################################
# print info
####################################################################
info:	
	@echo "SimlaQED routines"
	@echo "######################################"
	@echo "# DG Green     			    #"
	@echo "# Queen's University Belfast         #"
	@echo "# Copyright - August 2012            #"
	@echo "######################################"

####################################################################
# print info about set variable
####################################################################
variables_info:
	@echo "HOSTNAME    =" $(INT_HOST)
	@echo "FCOMP       =" $(FCOMP)
	@echo "FFLAGS      =" $(FFLAGS)
	@echo "F_LNK_FLAGS =" $(F_LNK_FLAGS)
	@echo "CORD        =" $(CORD)
	@echo "LIBRARIES   =" $(LIBRARIES)

####################################################################
# make all: chigamma
####################################################################
all: simla, photongen, testbessel, photonbin, gammatable

####################################################################
# recompile files when include files changed                       # 
# (using dummy "include.i"):                                       #
####################################################################

include.i: simla3_0.f90 simla_subroutines.f90 qedroutines.f90 
	touch include.i 
	rm -f simla3_0.o simla_subroutines.o qedroutines.o egen_sub.o d2pdXdt.o bessel.o
include.N: simla2_3.f90 simla_subroutines.f90 qedroutines.f90 besselNAG.f90 egen_sub.f90 d2pdXdt.f90
	touch include.N 
	rm -f simla2_3.o simla_subroutines.o qedroutines.o egen_sub.o d2pdXdt.o besselNAG.o
include.g: gammatable.f90 d2pdXdt.f90 besselNAG.f90 
	touch include.g 
	rm -f gammatable.o d2pdXdt.o besselNAG.o
include.t: photongen.f90 d2pdXdt.f90 bessel.f90 
	touch include.t 
	rm -f photongen.o d2pdXdt.o bessel.o
include.u: photonbin.f90 
	touch include.u 
	rm -f photonbin.o
include.b: testbessel.f90 besselNAG.f90 d2pdXdt.f90
	touch include.b 
	rm -f testbessel.o besselNAG.o d2pdXdt.o


####################################################################
# remove all object files and executable of simla                  +
####################################################################
clean:
	@rm -f simla $(FILES)
	@rm -f photongen $(TOOLSFILES)
	@rm -f gammatable $(TABLEFILES)
	@rm -f testbessel $(BTOOLSFILES)
	@rm -f photonbin $(BINFILES)
	@rm -f *.mod

vclean:
	@rm -f simla $(FILES)
	@rm -f photongen $(TOOLSFILES)
	@rm -f gammatable $(TABLEFILES)
	@rm -f testbessel $(BTOOLSFILES)
	@rm -f photonbin $(BINFILES)
	@rm -f *.mod
	@rm -f *.dat
####################################################################
# clean include.i                                                  #
####################################################################
clean_include:
	@rm -f include.i
	@rm -f include.t
	@rm -f include.b
	@rm -f include.u
	@rm -f include.N
	@rm -f include.g

####################################################################
# clean all
####################################################################
clean_all: clean clean_include




