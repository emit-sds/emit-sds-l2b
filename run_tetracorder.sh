 
REFL_FILE=${1}
SCALED_REFL=${2}
TETRA_OUT_DIR=${3}

#/t1 - Tetracorder directories and commands

python scale_and_cast_refl.py ${REFL_FILE} ${SCALED_REFL}
/t1/tetracorder.cmds/tetracorder5.2a.cmds/cmd-setup-tetrun ${TETRA_OUT_DIR} avirisng_2014b cube ${SCALED_REFL} -20 80 C .5 1.5 bar
 
