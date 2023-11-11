#!/bin/bash
# while [ 0 ]; do ./test.sh ; done
#
TART_API="https://tart.elec.ac.nz/signal"
TART_NCAL=3
TART_CAL_INT=13
DATA_DIR=./work
OUTPUT_DIR=./out 
TART_CAL_ELEVATION=50
TART_CAL_ITERATIONS=500
TART_CAL_POINTING=0
TART_CAL_POINTING_RANGE=4
TART_LOGIN_PW=sharkbait

rm -rf ${OUTPUT_DIR}
mkdir -p ${DATA_DIR}
rm -rf ${DATA_DIR}/*
python3 raw_cal/get_cal_data.py --api ${TART_API} --pw ${TART_LOGIN_PW} --n ${TART_NCAL} --interval ${TART_CAL_INT} --dir ${DATA_DIR}

# Perform optimization
# python3 raw_cal/pos_from_gps.py --api ${TART_API}  --elevation ${TART_CAL_ELEVATION}  --data ${DATA_DIR}  --dir ${OUTPUT_DIR}
python3 raw_cal/tart_cal.py --api ${TART_API}  \
                            --iterations ${TART_CAL_ITERATIONS} \
                            --pointing-range ${TART_CAL_POINTING_RANGE} \
                            --pointing ${TART_CAL_POINTING} \
                            --elevation ${TART_CAL_ELEVATION}  \
                            --data ${DATA_DIR}  \
                            --phases --dir ${OUTPUT_DIR}
# 
CALIB_OUTPUT=${OUTPUT_DIR}/BH_opt_json.json
echo "Calibration output is in ${CALIB_OUTPUT}"

tart_upload_gains --api ${TART_API} --gains ${CALIB_OUTPUT} --pw ${TART_LOGIN_PW}
tart_upload_antenna_positions --api ${TART_API} --file ${CALIB_OUTPUT} --pw ${TART_LOGIN_PW}
