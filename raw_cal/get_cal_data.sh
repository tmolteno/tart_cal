#!/bin/sh
DIR=cal_data
mkdir -p ${DIR}
API=https://tart.elec.ac.nz/signal
tart_set_mode --raw --api ${API} --pw ${PW}
sleep 60
tart_download_data --raw --n 1
FNAME=$(date --utc +"data_%Y-%m-%d_%H_%M")
echo $FNAME
tart_set_mode --vis --api ${API} --pw ${PW}
tart_calibration_data --api ${API} --n 1 --file ${DIR}/${FNAME}.json
mv ${FNAME}*.hdf ${DIR}/${FNAME}.hdf
