#!/bin/bash
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
CALIB_OUTPUT=${OUTPUT_DIR}/${TART_CAL_METHOD}_opt_json.json

# Set up default calibration parameters. These can be overridden by setting the env variables externally.
: "${TART_CAL_ITERATIONS:=700}"
: "${TART_CAL_ELEVATION:=45}"
: "${TART_CAL_INT:=30}"
: "${TART_CAL_POINTING:=0}"
: "${TART_CAL_ARGS:=""}"

echo "Working directory: ${DIR}"
echo "Data directory: ${DATA_DIR}"

# Get calibration data
python3 /app/get_cal_data.py --api ${TART_API} --pw ${TART_LOGIN_PW} --n ${TART_NCAL} --interval ${TART_CAL_INT} --dir ${DATA_DIR}

# Perform optimization
python3 /app/tart_cal.py --api ${TART_API} ${TART_CAL_ARGS} --elevation ${TART_CAL_ELEVATION} --iterations ${TART_CAL_ITERATIONS} --pointing ${TART_CAL_POINTING} --data ${DATA_DIR}  --dir ${OUTPUT_DIR}

# Log outputs
CAL_OUTPUT_FILE=${WORKING_DIR}/cal_${DATESTR}.json
mv ${CALIB_OUTPUT} ${CAL_OUTPUT_FILE}
echo "Calibration output is in ${CAL_OUTPUT_FILE}"


if [ ${TART_UPLOAD} == 1 ]; then
    echo "Uploading gains and antenna positions"

    tart_upload_gains --api ${TART_API} --gains ${CAL_OUTPUT_FILE} --pw ${TART_LOGIN_PW}
    tart_upload_antenna_positions --api ${TART_API} --file ${CAL_OUTPUT_FILE} --pw ${TART_LOGIN_PW}
fi
# Clean up
rm -r ${DIR}/*
rmdir ${DIR}
