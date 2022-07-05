#!/bin/sh
#
# Continuous Calibration Script. This should be run at regular intervals.
#
# Author Tim Molteno. tim@elec.ac.nz
#
# Reference: Molteno et al. Continuous Calibration of the TART using GPS satellites. EnzCon 2017.
#
WORKING_DIR=/work
TART_CAL_METHOD=BH # DE is really quick but does not produce good results.
DATESTR=`date "+%Y_%m_%d_%H_%M_%S"`
DIR=${WORKING_DIR}/calibration_${DATESTR}
mkdir -p ${DIR}

TART_API=https://tart.elec.ac.nz/${TARGET}
DATA_DIR=${DIR}/cal_data_${TARGET}
OUTPUT_DIR=${DIR}/cal_out_${TARGET}
CALIB_INPUT=${DIR}/calib_data.json
CALIB_OUTPUT=${DIR}/${TART_CAL_METHOD}_opt_json.json

POINTING=0

echo "Working directory: ${DIR}"
echo "Data directory: ${DATA_DIR}"

# Get calibration data
python3 /app/get_cal_data.py --api ${TART_API} --pw ${TART_LOGIN_PW} --n ${TART_NCAL} --interval ${TART_CAL_INT} --dir ${DATA_DIR}

# /usr/local/bin/tart_calibration_data --api ${TART_API} --n ${TART_NCAL} --i ${TART_CAL_INT} --file ${CALIB_INPUT}

# Perform optimization
python3 /app/tart_cal.py --api ${TART_API} --use-phases --elevation 40 --cold-start --pointing ${POINTING} --data ${DATA_DIR} --iterations 500 --dir ${OUTPUT_DIR}
# /usr/local/bin/tart_calibrate --api ${TART_API} --file ${CALIB_INPUT} --get-gains --method ${TART_CAL_METHOD} --dir ${DIR}

# Log outputs
CAL_OUTPUT_FILE=${WORKING_DIR}/cal_${DATESTR}.json
mv ${CALIB_OUTPUT} ${CAL_OUTPUT_FILE}
echo "Calibration output is in ${CAL_OUTPUT_FILE}"
/usr/local/bin/tart_upload_gains --api ${TART_API} --gains ${CAL_OUTPUT_FILE} --pw ${TART_LOGIN_PW}

echo "Uploading new antenna positions"
/usr/local/bin/tart_upload_antenna_positions --api ${TART_API} --file ${CAL_OUTPUT_FILE} --pw ${TART_LOGIN_PW}

# Clean up
rm ${DIR}/*
rmdir ${DIR}
