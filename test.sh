#!/bin/bash
# while [ 0 ]; do ./test.sh ; done
#

echo "Calibrating ${TARGET} pw=${TART_LOGIN_PW}"



TART_API="https://tart.elec.ac.nz/${TARGET}"

: "${TART_CAL_ITERATIONS:=500}"
: "${TART_CAL_ELEVATION:=45}"
: "${TART_NCAL:=3}"
: "${TART_CAL_INT:=30}"
: "${TART_CAL_POINTING:=0}"
: "${TART_CAL_POINTING_RANGE:=3}"
: "${TART_CAL_ARGS:=""}"

DATA_DIR=./work_${TARGET}
OUTPUT_DIR=./out_${TARGET} 

# 

if [ 0 ]
then
    rm -rf ${OUTPUT_DIR}
    mkdir -p ${DATA_DIR}
    rm -rf ${DATA_DIR}/*
    python3 raw_cal/get_cal_data.py --api ${TART_API} --pw ${TART_LOGIN_PW} --n ${TART_NCAL} --interval ${TART_CAL_INT} --dir ${DATA_DIR}
fi

# Perform optimization
# python3 raw_cal/pos_from_gps.py --api ${TART_API}  --elevation ${TART_CAL_ELEVATION}  --data ${DATA_DIR}  --dir ${OUTPUT_DIR}
python3 raw_cal/tart_cal.py --api ${TART_API}  \
                            --iterations ${TART_CAL_ITERATIONS} \
                            --pointing-range ${TART_CAL_POINTING_RANGE} \
                            --pointing ${TART_CAL_POINTING} \
                            --elevation ${TART_CAL_ELEVATION}  \
                            --data ${DATA_DIR}  \
                            --get-gains --phases \
                            --dir ${OUTPUT_DIR}
# 
CALIB_OUTPUT=${OUTPUT_DIR}/BH_opt_json.json
echo "Calibration output is in ${CALIB_OUTPUT}"

tart_upload_gains --api ${TART_API} --gains ${CALIB_OUTPUT} --pw ${TART_LOGIN_PW}
echo "Gains Uploaded"
tart_upload_antenna_positions --api ${TART_API} --file ${CALIB_OUTPUT} --pw ${TART_LOGIN_PW}
echo "Positions Uploaded"
