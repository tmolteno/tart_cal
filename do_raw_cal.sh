#!/bin/bash
# while [ 0 ]; do ./test.sh ; done
#
: "${TART_CAL_UPLOAD:=0}"
: "${CAL_CODE:='.'}"

if [ ${TART_CAL_UPLOAD} != 1 ]; then
    echo "########################################################"
    echo "##                                                    ##"
    echo "## WARNING! Calibration results will NOT BE UPLOADED  ##"
    echo "## set TART_CAL_UPLOAD=1 to upload                    ##"
    echo "##                                                    ##"
    echo "########################################################"
fi


TART_API="https://api.elec.ac.nz/tart/${TARGET}"

: "${TART_CAL_ITERATIONS:=500}"
: "${TART_CAL_ELEVATION:=45}"
: "${TART_NCAL:=3}"
: "${TART_CAL_INT:=30}"
: "${TART_CAL_POINTING:=0}"
: "${TART_CAL_POINTING_RANGE:=3}"
: "${TART_CAL_ARGS:=""}"
: "${TART_GET_DATA:=1}"
: "${TART_CAL_WORK_DIR:="."}"


DATA_DIR=${TART_CAL_WORK_DIR}/work_${TARGET}
OUTPUT_DIR=${TART_CAL_WORK_DIR}/out_${TARGET}

echo "##"
echo "##"
echo "## Calibrating ${TARGET} pw=${TART_LOGIN_PW}"
echo "## Pointing:   ${TART_CAL_POINTING} +/ ${TART_CAL_POINTING_RANGE}"
echo "## Iterations: ${TART_CAL_ITERATIONS}"
echo "## Elevation:  ${TART_CAL_ELEVATION}"
echo "## Args:       ${TART_CAL_ARGS}"
echo "## get_data:   ${TART_GET_DATA}"
echo "## Work Dir:   ${TART_CAL_WORK_DIR}"
echo "## Data Dir:   ${DATA_DIR}"
echo "## Output Dir: ${OUTPUT_DIR}"
echo "##"
echo "##"



if [ ${TART_GET_DATA} == 1 ]; then
    echo "## Downloading Data from ${TART_API}"
    echo "to ${DATA_DIR} output ${OUTPUT_DIR}:"

    rm -rf ${OUTPUT_DIR}
    mkdir -p ${OUTPUT_DIR}
    mkdir -p ${DATA_DIR}
    rm -rf ${DATA_DIR}/*
    get_cal_data --api ${TART_API} --pw ${TART_LOGIN_PW} --n ${TART_NCAL} --interval ${TART_CAL_INT} --dir ${DATA_DIR}
fi


# Perform optimization
# python3 raw_cal/pos_from_gps.py --api ${TART_API}  --elevation ${TART_CAL_ELEVATION}  --data ${DATA_DIR}  --dir ${OUTPUT_DIR}
raw_calibrate --api ${TART_API}  \
            --iterations ${TART_CAL_ITERATIONS} \
            --pointing-range ${TART_CAL_POINTING_RANGE} \
            --pointing ${TART_CAL_POINTING} \
            --elevation ${TART_CAL_ELEVATION}  \
            --data ${DATA_DIR}  \
            ${TART_CAL_ARGS} \
            --dir ${OUTPUT_DIR}
# 
CALIB_OUTPUT=${OUTPUT_DIR}/BH_opt_json.json
echo "##"
echo "## Calibration output is in ${CALIB_OUTPUT}"

if [ ${TART_CAL_UPLOAD} == 1 ]; then
    tart_upload_gains --api ${TART_API} --gains ${CALIB_OUTPUT} --pw ${TART_LOGIN_PW}
    echo "## Gains Uploaded"
    tart_upload_antenna_positions --api ${TART_API} --file ${CALIB_OUTPUT} --pw ${TART_LOGIN_PW}
    echo "## Positions Uploaded"
fi
