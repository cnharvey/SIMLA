###################################################################
# SIMLA Makefile                               			  		  #                                         	  #
###################################################################

###################################################################
# compiler flags                                                  #
###################################################################

LIBS = -lgsl -lgslcblas #GSL libraries (C)
CWRAP = synchrotrongateway.o spectrum_subroutine.o # Object files for c-wrapper
CWRAPEMPTY = spectrum_empty.o # Empth spectrum subroutine (used to avoid 
# dependence on external libraries)
CCOMP = gcc #Erik C-compiler
ALTHOST =
FCOMP = gfortran   ## ifort


FILES = simla.o simla_subroutines.o simla_qedroutines.o  


###################################################################
# FORTRAN compilation rule                                        #
###################################################################

.SUFFIXES: .f90 .o

# how to make a .o-file from a .f-file:
.f90.o:
	$(FCOMP)  -c $*.f90  

####################################################################
# link simla
####################################################################
simla: include.i $(FILES) $(CWRAPEMPTY) 
	$(FCOMP)  -o simla $(FILES) $(CWRAPEMPTY)


####################################################################
# link simla with synchrotron wrapper # Erik
####################################################################
synchrotrongateway.o: synchrotrongateway.c
	$(CCOMP) -c synchrotrongateway.c

simla_synchrotron: include.i $(FILES) $(CWRAP)
	$(FCOMP)  -o simla $(FILES) $(CWRAP) $(LIBS) 



####################################################################
# recompile files when include files changed                       # 
# (using dummy "include.i"):                                       #
####################################################################

include.i: simla.f90 simla_subroutines.f90 simla_qedroutines.f90 
	touch include.i 
	rm -f simla.o simla_subroutines.o simla_qedroutines.o 


####################################################################
# clean up commands                  							   #
####################################################################
clean:
	@rm -f simla $(FILES)
	@rm -f *.mod
	@rm -f spectrum_empty.o spectrum_subroutine.o synchrotrongateway.o

vclean:
	@rm -f simla $(FILES)
	@rm -f *.mod
	@rm -f *.dat
	@rm -f include.i





