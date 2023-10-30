
TART_API="https://tart.elec.ac.nz/rhodes"
TART_NCAL=3
TART_CAL_INT=15
DATA_DIR=./work
OUTPUT_DIR=./out 
TART_CAL_ELEVATION=55
TART_CAL_ITERATIONS=800
TART_CAL_POINTING=0

rm -rf ${OUTPUT_DIR}
mkdir -p ${DATA_DIR}
rm -rf ${DATA_DIR}/*
python3 raw_cal/get_cal_data.py --api ${TART_API} --pw sharkbait --n ${TART_NCAL} --interval ${TART_CAL_INT} --dir ${DATA_DIR}

# Perform optimization
# python3 raw_cal/pos_from_gps.py --api ${TART_API}  --elevation ${TART_CAL_ELEVATION}  --data ${DATA_DIR}  --dir ${OUTPUT_DIR}
python3 raw_cal/tart_cal.py --api ${TART_API} --phases --iterations ${TART_CAL_ITERATIONS} --pointing ${TART_CAL_POINTING} --elevation ${TART_CAL_ELEVATION}  --data ${DATA_DIR}  --dir ${OUTPUT_DIR}
