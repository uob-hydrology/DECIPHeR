# ============================================================================
# Name        : Makefile
# Author      : Toby Dunne
# Version     :
# Copyright   : 
# Description : Make all the dta utilities
#             : Supports debug and release builds
#             : Release builds have .e extension, debug builds have .d extension
# ============================================================================
FC             = gfortran
FFLAGS         = -O0 -std=f2003 -Wall -Wextra -fbounds-check -g -fbacktrace -mcmodel=large
FFLAGS_RELEASE = -std=f2003 -Wall -mcmodel=large
LDFLAGS        = 
PROG_RM        = rm -f

RIVER_THRESH_OBJS = river_thresholds.o dta_utility.o
RIVER_FIND_HW_OBJS = river_find_hw.o dta_rivers.o dta_queue.o dta_stack.o dta_utility.o
RIVER_RUN_OBJS = river_run.o dta_rivers.o dta_queue.o dta_stack.o dta_utility.o
CATCH_FIND_OBJS = catch_find.o dta_utility.o dta_catch_cut.o dta_queue.o dta_stack.o dta_riv_tree_node.o
CATCH_CUT_OBJS = catch_cut.o dta_utility.o dta_catch_cut.o dta_queue.o dta_stack.o dta_riv_tree_node.o
CATCH_MASK_OBJS = catch_mask.o dta_catch_mask.o dta_utility.o dta_riv_tree_node.o dta_queue.o dta_stack.o
ROUTE_TREE_OBJS = route_tree.o dta_route_tree.o dta_rivers.o dta_queue.o dta_stack.o dta_utility.o dta_riv_tree_node.o dta_catch_cut.o
ATB_OBJS        = atb_wfp.o dta_utility.o dta_atb_wfp.o dta_stack.o dta_qsort.o
ROUTE_RIVER_FILE_OBJS = route_river_file.o dta_routing_file.o dta_rivers.o dta_queue.o dta_stack.o dta_utility.o dta_riv_tree_node.o
ROUTE_PROC_STD_OBJS = route_run_standalone.o dta_route_processing.o dta_utility.o dta_riv_tree_node.o
MASK_CHECK_OBJS = mask_check.o dta_utility.o
PREPROCESS_OBJS = preprocess.o dta_utility.o dta_preprocess.o
CALC_HRUS_OBJS = calculate_hrus.o dta_utility.o dta_calculate_hrus.o dta_qsort.o dta_preprocess.o

.SUFFIXES: 
.SUFFIXES: .f90 .o

$(shell mkdir -p debug release)

.PHONY: all clean debug release

RELEASE_TARGETS = river_thresholds.e river_find_hw.e river_run.e \
                  catch_cut.e catch_find.e catch_mask.e \
                  mask_check.e \
			      route_tree.e \
			      atb_wfp.e \
			      route_river_file.e \
			      route_run_standalone.e \
                  preprocess.e \
                  calculate_hrus.e
			      
DEBUG_TARGETS = $(RELEASE_TARGETS:%.e=%.d)

all: debug release

debug: $(DEBUG_TARGETS)
	echo $(DEBUG_TARGETS)
	
release: $(RELEASE_TARGETS)
	echo $(RELEASE_TARGETS)

debug/%.o: %.f90
	-$(FC) -c $(FFLAGS) -Jdebug \
	$< -o $(<:%.f90=debug/%.o)
	
release/%.o: %.f90
	-$(FC) -c $(FFLAGS_RELEASE) -Jrelease \
	$< -o $(<:%.f90=release/%.o)

clean: /;
	-$(PROG_RM) *.d *.e *.mod *.o
	-$(PROG_RM) .depend
	-$(PROG_RM) debug/*
	-$(PROG_RM) release/*

clean_obj: /;
	-$(PROG_RM) *.mod *.o
	-$(PROG_RM) .depend
	-$(PROG_RM) debug/*
	-$(PROG_RM) release/*


# Create a dependency list using makedepf90.  All files  
#  # that needs to be compiled to build the program,  
# i.d all source files except include files, should  
#  # be given on the command line to makedepf90.   
# 
#  # The argument to the '-o' option will be the name of the 
# resulting program when running 'make', in this case  
#  # 'foobar'
depend: .depend
 
.depend: 
	makedepf90 -b debug *.f90 > .depend
	makedepf90 -b release *.f90 >> .depend

# dependencies file is recreated if any .f90 file changes
.depend: *.f90

##   DEBUG
preprocess.d: $(PREPROCESS_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(PREPROCESS_OBJS:%.o=debug/%.o) $(LIBS)

calculate_hrus.d: $(CALC_HRUS_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(CALC_HRUS_OBJS:%.o=debug/%.o) $(LIBS)


route_run_standalone.d : $(ROUTE_PROC_STD_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(ROUTE_PROC_STD_OBJS:%.o=debug/%.o) $(LIBS)
	
route_river_file.d: $(ROUTE_RIVER_FILE_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(ROUTE_RIVER_FILE_OBJS:%.o=debug/%.o) $(LIBS)
	
route_tree.d: $(ROUTE_TREE_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(ROUTE_TREE_OBJS:%.o=debug/%.o) $(LIBS)

river_thresholds.d: $(RIVER_THRESH_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(RIVER_THRESH_OBJS:%.o=debug/%.o) $(LIBS)

river_find_hw.d: $(RIVER_FIND_HW_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(RIVER_FIND_HW_OBJS:%.o=debug/%.o) $(LIBS)

river_run.d: $(RIVER_RUN_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(RIVER_RUN_OBJS:%.o=debug/%.o) $(LIBS)

atb_wfp.d: $(ATB_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(ATB_OBJS:%.o=debug/%.o) $(LIBS)
	
catch_find.d: $(CATCH_FIND_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(CATCH_FIND_OBJS:%.o=debug/%.o) $(LIBS)	
	
catch_cut.d: $(CATCH_CUT_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(CATCH_CUT_OBJS:%.o=debug/%.o) $(LIBS)
	
catch_mask.d: $(CATCH_MASK_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(CATCH_MASK_OBJS:%.o=debug/%.o) $(LIBS)

mask_check.d: $(MASK_CHECK_OBJS:%.o=debug/%.o)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(MASK_CHECK_OBJS:%.o=debug/%.o) $(LIBS)

##   RELEASE
preprocess.e: $(PREPROCESS_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(PREPROCESS_OBJS:%.o=release/%.o) $(LIBS)

calculate_hrus.e: $(CALC_HRUS_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(CALC_HRUS_OBJS:%.o=release/%.o) $(LIBS)

route_run_standalone.e: $(ROUTE_PROC_STD_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(ROUTE_PROC_STD_OBJS:%.o=release/%.o) $(LIBS)

route_river_file.e: $(ROUTE_RIVER_FILE_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(ROUTE_RIVER_FILE_OBJS:%.o=release/%.o) $(LIBS)
	
route_tree.e: $(ROUTE_TREE_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(ROUTE_TREE_OBJS:%.o=release/%.o) $(LIBS)

river_thresholds.e: $(RIVER_THRESH_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(RIVER_THRESH_OBJS:%.o=release/%.o) $(LIBS)

river_find_hw.e: $(RIVER_FIND_HW_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(RIVER_FIND_HW_OBJS:%.o=release/%.o) $(LIBS)

river_run.e: $(RIVER_RUN_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(RIVER_RUN_OBJS:%.o=release/%.o) $(LIBS)

atb_wfp.e: $(ATB_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(ATB_OBJS:%.o=release/%.o) $(LIBS)

catch_find.e: $(CATCH_FIND_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(CATCH_FIND_OBJS:%.o=release/%.o) $(LIBS)

catch_cut.e: $(CATCH_CUT_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(CATCH_CUT_OBJS:%.o=release/%.o) $(LIBS)

catch_mask.e: $(CATCH_MASK_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(CATCH_MASK_OBJS:%.o=release/%.o) $(LIBS)	

mask_check.e: $(MASK_CHECK_OBJS:%.o=release/%.o)
	$(FC) -o $@ $(FFLAGS_RELEASE) $(LDFLAGS) $(MASK_CHECK_OBJS:%.o=release/%.o) $(LIBS)

include .depend
