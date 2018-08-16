FC             = ifort
FFLAGS         =  -g -check all -check bounds -check uninit -ftrapuv -debug all -traceback
LDFLAGS        = -g -check all -check bounds -check uninit -ftrapuv -debug all -traceback
PROG_RM        = rm -f
PROGRAM        = DECIPHeR_v1.exe

.SUFFIXES: 
.SUFFIXES: .f90 .o

$(PROGRAM): .depend

all: $(PROGRAM)

%.o: %.f90
	-$(FC) -c $(FFLAGS) \
	$<

%.o: ../%.f90
	-$(FC) -c $(FFLAGS) \
	$<

include .depend

clean: /;
	$(PROG_RM) $(PROGRAM) *.o *.mod

clean_obj: /;
	$(PROG_RM) *.o *.mod


# Create a dependency list using makedepf90.  All files  
#  # that needs to be compiled to build the program,  
# i.e all source files except include files, should  
#  # be given on the command line to makedepf90.   
# 
#  # The argument to the '-o' option will be the name of the 
# resulting program when running 'make', in this case  
#  # 'foobar' 
depend .depend: 
	makedepf90 -b ./ -o $(PROGRAM) ../*.f90 > .depend

# dependencies file is recreated if any .f90 file changes
.depend: ../*.f90

# DO NOT DELETE
