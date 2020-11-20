 
REFL_FILE=${1}
SCALED_REFL=${2}
SCALED_REFL_TETRA_REF=${3}
TETRA_OUT_DIR=${4}
AGGREGATED_OUT_FILE=${5}

echo 'starting'
date

TETRA_OUT_DIR_ABS=`readlink -f ${TETRA_OUT_DIR}`

TETRA_EXE=/t1/tetracorder.cmds/tetracorder5.2a.cmds/cmd-setup-tetrun
TETRA_EXPERT_SYSTEM=cmd.lib.setup.t5.2a2

#AVCL
#TETRA_ENVI_SL=/beegfs/scratch/brodrick/emit/emit-sds-l2b/convert_tetracorder/Spectral-Library-Reader-master/output_spectral_library_avcl
#ANG
TETRA_ENVI_SL=/beegfs/scratch/brodrick/emit/emit-sds-l2b/convert_tetracorder/Spectral-Library-Reader-master/output_spectral_library

#/t1 - Tetracorder directories and commands
echo 'scaling'
date
export GDAL_CACHEMAX=1
python scale_and_cast_refl.py ${REFL_FILE} ${SCALED_REFL}

echo 'tetracorder'
date
#AVCL
#${TETRA_EXE} ${TETRA_OUT_DIR} aviris_2013 cube ${SCALED_REFL_TETRA_REF} -20 80 C .5 1.5 bar
#ANG
${TETRA_EXE} ${TETRA_OUT_DIR} avirisng_2014b cube ${SCALED_REFL_TETRA_REF} -20 80 C .5 1.5 bar
echo 20 > ${TETRA_OUT_DIR}/TETNCPU.txt; 

current_dir=${PWD}
cd ${TETRA_OUT_DIR}
time ./cmd.runtet cube >& cmd.runtet.out 
cd ${current_dir}
 

echo 'aggregation'
date
python /beegfs/scratch/brodrick/emit/emit-sds-l2b/convert_tetracorder/aggregator.py ${TETRA_EXPERT_SYSTEM} ${TETRA_ENVI_SL} ${TETRA_OUT_DIR_ABS} /beegfs/scratch/brodrick/emit/emit-sds-l2b/convert_tetracorder/mineral_fractions_s06an14a.csv ${AGGREGATED_OUT_FILE} -calculate_uncertainty 1 -reflectance_file ${REFL_FILE} -reflectance_uncertainty_file ${6}
#TETRA_EXPERT_SYSTEM, TETRA_LIBRARY_FILE, TETRA_OUTPUT_DIR, MINERAL_FRACTIONS_FILE, OUTPUT


# Spectral library build command
#./specpr2envi /sl1/usgs/library06.conv/s06av13a output_spectral_library_avcl
