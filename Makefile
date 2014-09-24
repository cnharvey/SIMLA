###################################################################
# SIMLA Makefile                               			  		  #                                         	  #
###################################################################

###################################################################
# compiler flags                                                  #
###################################################################

ALTHOST =
FCOMP       = gfortran   ## ifort


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
simla: include.i $(FILES) 
	$(FCOMP)  -o simla $(FILES) 


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

vclean:
	@rm -f simla $(FILES)
	@rm -f *.mod
	@rm -f *.dat
	@rm -f include.i





