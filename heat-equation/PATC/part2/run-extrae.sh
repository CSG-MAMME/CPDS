#!/bin/sh

export OMP_NUM_THREADS=8
export EXTRAE_CONFIG_FILE=extrae.xml
#export EXTRAE_HOME=/apps/CEPBATOOLS/extrae/2.3/bullmpi/64
source ${EXTRAE_HOME}/etc/extrae.sh

${EXTRAE_HOME}/bin/extrae -v ./heatomp test.dat
