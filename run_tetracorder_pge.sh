#!/bin/bash

SCRATCH_DIR=$1
RFL_PATH=$2
LOCAL_OUTPUT=`basename ${SCRATCH_DIR}`
TMP_OUTPUT_PATH=/tmp/${LOCAL_OUTPUT}
RFL_NAME=`basename ${RFL_PATH}`
TMP_RFL_PATH=/tmp/${RFL_NAME}
TETRA_SETUP=/shared/tetracorder5.26a.cmds/cmd-setup-tetrun

# Create rfl symlinks
ln -s ${RFL_PATH} ${TMP_RFL_PATH}
ln -s ${RFL_PATH}.hdr ${TMP_RFL_PATH}.hdr

date
${TETRA_SETUP} ${TMP_OUTPUT_PATH} aviris_2018 cube ${TMP_RFL_PATH} 1 -T -20 80 C -P .5 1.5 bar
cd ${TMP_OUTPUT_PATH}
time ./cmd.runtet cube ${TMP_RFL_PATH}  >& cmd.runtet.out
echo 'tetracorder finished'
date

# Trailing slash is needed here to copy and rename the directory
cp -r ${TMP_OUTPUT_PATH}/ ${SCRATCH_DIR}/tetra_out

rm -rf ${TMP_OUTPUT_PATH}
rm -f ${TMP_RFL_PATH}
rm -f ${TMP_RFL_PATH}.hdr
