#!/bin/bash
## This script shows how to run the DTA from scratch.

# Before running these commands:
#  1. Make sure you have all the necessary input files to run the DTA 
#  2. Copy this script and the source code onto your workspace 
#  3. Change the MAIN_DIR to point to your working folder 
#  4. Change the CODE_DIR to point to your source code folder
#  5. Change the names of the DEM, gauge list and river file names

echo ' '
echo '#=============================================================================='
echo '## STEP 0: DEFINE PATHWAYS TO ALL INPUTS'
echo '#=============================================================================='
echo ' '
echo 'Make sure you have all neccessary files and have changed path directories in the bash script'
read -p 'Press enter if you are happy to continue'
echo ' '

# -----------------------------------------------------------------------------
# Define paths to directories and filenames
# -----------------------------------------------------------------------------

# Define path to MAIN_DIR
# This is your working folder which must contain all input files. 
# THIS MUST BE CHANGED.
MAIN_DIR=/data/DYMOND/DTA

# Full directory to folder containing source code
# THIS MUST BE CHANGED.
CODE_DIR=/data/DYMOND/DTA/SOURCE_CODE

#necessary files:
#THESE MUST BE CHANGED TO REFER TO YOUR DEM AND GAUGE LIST
#Currently looking for a file called DEM.asc and gauge_list.txt
DEM=DEM
GAUGES=gauge_list.txt

#Optional files:
#THIS MUST BE CHANGED TO REFER TO YOUR RIVER NETWORK
#Currently looking for a file called RIVER.asc 
RIV=RIVER

# -----------------------------------------------------------------------------
# Change directory to main directory - this is where input and output files are
# -----------------------------------------------------------------------------

# All lines of code should be run from your main directory. 
cd ${MAIN_DIR}
echo MAIN_DIR = ${MAIN_DIR} 
echo CODE_DIR = ${CODE_DIR}
echo DEM = ${MAIN_DIR}/${DEM}
echo GAUGES = ${MAIN_DIR}/${GAUGES}
read -p 'Check above file locations are correct, and press Enter to continue'

echo ' '
echo '#=============================================================================='
echo 'STEP 1: CREATING THE RIVER NETWORK'
echo '#=============================================================================='

# -----------------------------------------------------------------------------
# 2 OPTIONS - IF STATEMENT ASKS USER IF THEY WANT TO USE AN EXTERNAL RIVER NETWORK
# -----------------------------------------------------------------------------
echo ''
echo 'Option 1: a river map is created using the DEM - no external map required'
echo 'Option 2: External river network is used'
read -p 'Would you like to use Option 1 or Option 2? [1/2]' option_river

if [ "$option_river" -eq "1" ] 
then

	echo ''
	echo '# ---------------------------------------------------------------------'
	echo '# OPTION 1: river network will be created using DEM data'
	echo '# ----------------------------------------------------------------------'

	echo '${CODE_DIR}/atb_wfp.e -dem ${DEM}.asc'
	#read -p 'Press enter to run the above command and calculate slope/area and topographic index'
	${CODE_DIR}/atb_wfp.e -dem ${DEM}.asc

	echo ' '
	echo ${CODE_DIR}/river_thresholds.e -dem ${DEM}.asc -atb ${DEM}_atb.asc -area ${DEM}_area.asc
	#read -p 'Press enter to run above command and get river thresholds'
	${CODE_DIR}/river_thresholds.e -dem ${DEM}.asc -atb ${DEM}_atb.asc -area ${DEM}_area.asc

	RIV=river_thresholds

else
	echo ''
	echo '# ---------------------------------------------------------------------'
	echo '# OPTION 2: Using an external river network'
	echo '# ----------------------------------------------------------------------'
	echo river network = ${MAIN_DIR}/${GAUGES}
fi

echo ' '
echo '# -----------------------------------------------------------------------------'
echo 'river_find_hw.e is used to find the headwaters from the river network'
echo '# -----------------------------------------------------------------------------'

echo ${CODE_DIR}/river_find_hw.e -dem ${DEM}.asc -river ${RIV}.asc
${CODE_DIR}/river_find_hw.e -dem ${DEM}.asc -river ${RIV}.asc

echo ' '
echo '# -----------------------------------------------------------------------------'
echo '# river_run.e creates river layers for Dynamic TOPMODEL'
echo '# ----------------------------------------------------------------------------- '

echo ${CODE_DIR}/river_run.e -dem ${DEM}.asc -headwater ${RIV}_HW_500m_100m.txt
${CODE_DIR}/river_run.e -dem ${DEM}.asc -headwater ${RIV}_HW_500m_100m.txt

echo ' '
echo '#=============================================================================='
echo 'STEP 2: CATCHMENT IDENTIFICATION AND MASKS'
echo '#=============================================================================='

echo ' '
echo '# -----------------------------------------------------------------------------'
echo '# catch_find.e: Generates catchment masks from XY gauge list'
echo '# -----------------------------------------------------------------------------'

echo ${CODE_DIR}/catch_find.e -dem ${DEM}.asc -river ${DEM}_riv.asc -stations ${GAUGES}
${CODE_DIR}/catch_find.e -dem ${DEM}.asc -river ${DEM}_riv.asc -stations ${GAUGES}

echo ' '
echo '# -----------------------------------------------------------------------------'
echo '# route_tree.e links all gauge locations on the river network'
echo '# -----------------------------------------------------------------------------'

echo ${CODE_DIR}/route_tree.e -dem ${DEM}.asc -river ${DEM}_riv.asc -points ${DEM}_station_match.txt
${CODE_DIR}/route_tree.e -dem ${DEM}.asc -river ${DEM}_riv.asc -points ${DEM}_station_match.txt

echo ' '
echo '# -----------------------------------------------------------------------------'
echo '# catch_cut.e and catch_mask.e create individual and combined catchment masks'
echo '# -----------------------------------------------------------------------------'

echo '# builds all the catchment masks that are used in the river tree'
mkdir masks
echo mkdir masks
echo ${CODE_DIR}/catch_cut.e -dem ${DEM}.asc -points ${DEM}_flow_point.txt -out masks
${CODE_DIR}/catch_cut.e -dem ${DEM}.asc -points ${DEM}_flow_point.txt -out masks

echo ' '
echo '#combine catchment masks based on river network tree order from flow_conn.txt'
echo ${CODE_DIR}/catch_mask.e -base ${DEM}.asc -tree ${DEM}_flow_conn.txt -mask_dirs masks
${CODE_DIR}/catch_mask.e -base ${DEM}.asc -tree ${DEM}_flow_conn.txt -mask_dirs masks 

echo ' '
echo ${CODE_DIR}/mask_check.e -mask ${DEM}_mask.asc -riv_id ${DEM}_riv_id.asc
${CODE_DIR}/mask_check.e -mask ${DEM}_mask.asc -riv_id ${DEM}_riv_id.asc

echo ' '
echo '#=============================================================================='
echo 'STEP 3: TOPOGRAPHIC INDEX'
echo '#=============================================================================='

echo ' ' 
echo '# -----------------------------------------------------------------------------'
echo '# atb_wfp.e is used to calculate slope, accumulated area and topographic index'
echo '# -----------------------------------------------------------------------------'

echo ${CODE_DIR}/atb_wfp.e -dem ${DEM}.asc -river ${DEM}_riv_id_check.asc -mask ${DEM}_mask.asc
${CODE_DIR}/atb_wfp.e -dem ${DEM}.asc -river ${DEM}_riv_id_check.asc -mask ${DEM}_mask.asc

echo ' '
echo '#=============================================================================='
echo 'STEP 4: RIVER ROUTING'
echo '#=============================================================================='

echo ' '
echo '# -----------------------------------------------------------------------------'
echo '# route_river_file.e outputs river cell information for the routing algorithm'
echo '# -----------------------------------------------------------------------------'

echo ${CODE_DIR}/route_river_file.e -dem ${DEM}.asc -river ${DEM}_riv_id_check.asc
${CODE_DIR}/route_river_file.e -dem ${DEM}.asc -river ${DEM}_riv_id_check.asc



