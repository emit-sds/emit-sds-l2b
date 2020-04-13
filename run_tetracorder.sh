 
REFL_FILE=${1}
SCALED_REFL=${2}
SCALED_REFL_TETRA_REF=${3}
TETRA_OUT_DIR=${4}
AGGREGATED_OUT_FILE=${5}

TETRA_OUT_DIR_ABS=`readlink -f ${TETRA_OUT_DIR}`

TETRA_EXPERT_SYSTEM=/t1/tetracorder.cmds/tetracorder5.2a.cmds/cmd.lib.setup.t5.2a2
TETRA_ENVI_SL=/beegfs/scratch/brodrick/emit/emit-sds-l2b/convert_tetracorder/Spectral-Library-Reader-master/output_spectral_library

#/t1 - Tetracorder directories and commands
export GDAL_CACHEMAX=1
python scale_and_cast_refl.py ${REFL_FILE} ${SCALED_REFL}
${TETRA_EXPERT_SYSTEM} ${TETRA_OUT_DIR} aviris_2013 cube ${SCALED_REFL_TETRA_REF} -20 80 C .5 1.5 bar
echo 20 > ${TETRA_OUT_DIR}/TETNCPU.txt; 

cd ${TETRA_OUT_DIR}
time ./cmd.runtet cube >& cmd.runtet.out 
 

python /beegfs/scratch/brodrick/emit/emit-sds-l2b/convert_tetracorder/aggregator.py ${TETRA_EXPERT_SYSTEM} ${TETRA_ENVI_SL} ${TETRA_OUT_DIR_ABS} /beegfs/scratch/brodrick/emit/emit-sds-l2b/convert_tetracorder/mineral_fractions_s06an14a.csv ${AGGREGATED_OUT_FILE}
#TETRA_EXPERT_SYSTEM, TETRA_LIBRARY_FILE, TETRA_OUTPUT_DIR, MINERAL_FRACTIONS_FILE, OUTPUT
