 
TETRA_OUT_DIR=${1}
AGGREGATED_OUT_FILE=${2}
REFL_FILE=${3}
REFL_UNCERT_FILE=${4}
EMIT_UTILS_PATH=${5}
EMIT_SDS_L2B_PATH=${6}
TETRACORDER_CMD_BASE=${7}

echo 'starting'
date

TETRA_OUT_DIR_ABS=`readlink -f ${TETRA_OUT_DIR}`
TETRA_SETUP=${TETRACORDER_CMD_BASE}/cmd-setup-tetrun

REFL_ABS_FILE=`readlink -f ${REFL_FILE}`
OUTPUT_ABS_DIR=`readlink -f ${TETRA_OUT_DIR}`
AGGREGATED_OUT_FILE_ABS=`readlink -f ${AGGREGATED_OUT_FILE}`
local_refl=`basename ${REFL_FILE}`
local_output=`basename ${TETRA_OUT_DIR}`
current_dir=${PWD}

cd /tmp

cp ${REFL_ABS_FILE} .
cp ${REFL_ABS_FILE}.hdr .

echo 'coppied'
date
${TETRA_SETUP} /tmp/${local_output} aviris_2018 cube /tmp/${local_refl} 1 -T -20 80 C -P .5 1.5 bar
cd ${local_output}
time ./cmd.runtet cube /tmp/${local_refl}  >& cmd.runtet.out
echo 'tetracorder finished'
date

cp /tmp/${local_output} ${OUTPUT_ABS_DIR} -r

rm -rf /tmp/${local_output}
rm -rf /tmp/${local_refl}

cd ${current_dir}

echo 'cleanup finished'
date

export PYTHONPATH=$PYTHONPATH:${EMIT_SDS_L2B_PATH}
export PYTHONPATH=$PYTHONPATH:${EMIT_UTILS_PATH}

cd ${EMIT_SDS_L2B_PATH}
python aggregator.py ${TETRA_OUT_DIR_ABS}  ${AGGREGATED_OUT_FILE_ABS} --calculate_uncertainty 0 --reflectance_file ${REFL_FILE} --reflectance_uncertainty_file ${REFL_UNCERT_FILE}

cd ${current_dir}
echo 'aggregation finished'
date



