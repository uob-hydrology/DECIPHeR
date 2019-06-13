## This is a follow-on script from Run_DTA_Example_Script_1.sh 
#  before following this script:
#	1. Run Run_DTA_Example_Script_1.sh 
#	2. Create an output folder in your home directory
#	3. Ensure you have a HRU classifiers file and a list of gauges you wish to process
#	4. Modify the neccessary lines in this script (paths and input variables for step 2). 

echo ''
echo '#=============================================================================='
echo '# STEP 0: DEFINE PATHWAYS AND KEY VARIABLES'
echo '#=============================================================================='
echo ' '

# -----------------------------------------------------------------------------
# Define paths to directories and filenames
# -----------------------------------------------------------------------------

# Full directory to folder containing source code
# THIS MUST BE CHANGED.
CODE_DIR=/data/DYMOND/DTA/SOURCE_CODE

# Your working folder for the first part of the DTA analysis and the root file name for the DEM
# THIS MUST BE CHANGED
ROOT_FN=/data/DYMOND/DTA/DEM

#The output folder you want to put the results in 
# THIS MUST BE CHANGED
OUTPUT_DIR=/data/DYMOND/DTA/HRU

#Text file listing NRFA Gauge IDs that you want to model
GAUGELIST=${OUTPUT_DIR}/gauge_list.txt

#HRU Class file listing which classifiers you want to use to split up the catchment
HRU_CLASS_FILE=${OUTPUT_DIR}/hru_class.dat

#Catchment mask file for this gauge ID (full filepath)
# THIS MUST BE CHANGED
CATCHMASK_DIR=/data/DYMOND/DTA/masks/

# print these out to screen to prompt user to check them. 
echo CODE_DIR = ${CODE_DIR}
echo ROOT_FN = ${ROOT_FN} 
echo OUTPUT_DIR = ${OUTPUT_DIR}
echo GAUGELIST = ${GAUGELIST}
echo HRU_CLASS_FILE = ${HRU_CLASS_FILE}
echo CATCHMASK_DIR = ${CATCHMASK_DIR}

read -p 'Check above file locations are correct, and press Enter to continue'

echo ''
echo '#=============================================================================='
echo '# STEP 1: PREPROCESSING CODE'
echo '#=============================================================================='
echo ''
echo 'preprocess.e crops relevant data for specified catchment(s) only'
echo ''

echo ${CODE_DIR}/./preprocess.e -gaugelist ${GAUGELIST} -hru_class_file ${HRU_CLASS_FILE} -root_fn ${ROOT_FN} -output_folder ${OUTPUT_DIR} -catchmask_folder ${CATCHMASK_DIR}
${CODE_DIR}/./preprocess.e -gaugelist ${GAUGELIST} -hru_class_file ${HRU_CLASS_FILE} -root_fn ${ROOT_FN} -output_folder ${OUTPUT_DIR} -catchmask_folder ${CATCHMASK_DIR}


echo ''
echo '#=============================================================================='
echo '## STEP 2: CALCULATE HYDROLOGICAL RESPONSE UNITS'
echo '#=============================================================================='
echo ''

# -----------------------------------------------------------------------------
# calculate HRUs
# -----------------------------------------------------------------------------

echo ${CODE_DIR}/./calculate_hrus.e -gaugelist ${GAUGELIST} -hru_class_file ${HRU_CLASS_FILE} -output_folder ${OUTPUT_DIR}/
${CODE_DIR}/./calculate_hrus.e -gaugelist ${GAUGELIST} -hru_class_file ${HRU_CLASS_FILE} -output_folder ${OUTPUT_DIR}/
