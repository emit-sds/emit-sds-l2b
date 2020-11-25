 
TETRA_OUT_DIR=${1}
AGGREGATED_OUT_FILE=${2}
REFL_FILE=${3}
REFL_UNCERT_FILE=${4}

echo 'starting'
date

TETRA_OUT_DIR_ABS=`readlink -f ${TETRA_OUT_DIR}`
TETRA_SETUP=/shared/tetracorder5.26a.cmds/cmd-setup-tetrun

REFL_ABS_FILE=`readlink -f ${REFL_FILE}`

echo 'tetracorder'
date
#AVCL
#${TETRA_SETUP} ${TETRA_OUT_DIR} aviris_2018 cube ${REFL_FILE} -20 80 C .5 1.5 bar
#ANG
${TETRA_SETUP} ${TETRA_OUT_DIR} avirisng_2014b cube ${REFL_FILE} 1 -T -20 80 C -P .5 1.5 bar

current_dir=${PWD}
cd ${TETRA_OUT_DIR}
time ./cmd.runtet cube ${REFL_ABS_FILE} >& cmd.runtet.out
cd ${current_dir}
 

echo 'aggregation'
date
python aggregator.py ${TETRA_OUT_DIR_ABS}  ${AGGREGATED_OUT_FILE} -calculate_uncertainty 1 -reflectance_file ${REFL_FILE} -reflectance_uncertainty_file ${REFL_UNCERT_FILE}


