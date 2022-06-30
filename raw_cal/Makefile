TARGET=signal
API=https://tart.elec.ac.nz/${TARGET}
PW=<password>
N=1
INTERVAL=10
DATA_DIR=cal_data_${TARGET}
cal:
	python3 tart_cal.py --api ${API} --data ${DATA_DIR} --iterations 500 --ignore 14 --dir cal_out_${TARGET}

	
image:
	python3 image_from_raw.py --file cal_raw_data.hdf


data:
	rm -f ${DATA_DIR}/*
	python3 get_cal_data.py --api ${API} --pw ${PW} --n ${N} --interval ${INTERVAL} --dir ${DATA_DIR}
