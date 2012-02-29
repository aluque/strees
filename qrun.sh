#! /bin/bash
# Script to start a simulation in the qsub queue
# Alejandro Luque - 2012

PYTHONVER=2.7
PYTHON_EXEC=python${PYTHONVER}
UNAME_FULL=`uname -a`
DATE=`date`
export LD_LIBRARY_PATH=${HOME}/lib

echo "@ [${DATE}]"
echo "@ ${UNAME_FULL}"

if [ -z $FMM_PATH ]; then
    FMM_PATH=~/fmm
fi

MAIN=${FMM_PATH}/grow_tree.py

${PYTHON_EXEC} ${MAIN} ${FMM_INPUT_FILE}

