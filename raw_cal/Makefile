# Full raw data calibration using GPS correlations to flag antennas.
#
# Example usage
#
# make data TARGET=rhodes N=3 PW=foobar
# make cal TARGET=rhodes POINTING=11

TARGET=signal
PW=<password>
N=1
INTERVAL=10
POINTING=0

API=https://tart.elec.ac.nz/${TARGET}
DATA_DIR=cal_data_${TARGET}

plot:
	python3 plot_antenna_array.py --data ${DATA_DIR}

cal:
	rm -f cal_out_${TARGET}/*
	python3 tart_cal.py --api ${API} --use-phase --cold-start --pointing ${POINTING} --data ${DATA_DIR} --iterations 500 --dir cal_out_${TARGET}

	
image:
	python3 image_from_raw.py --file cal_raw_data.hdf


data:
	rm -f ${DATA_DIR}/*
	python3 get_cal_data.py --api ${API} --pw ${PW} --n ${N} --interval ${INTERVAL} --dir ${DATA_DIR}
