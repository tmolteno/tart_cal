TARGET=signal
API=https://tart.elec.ac.nz/${TARGET}
PW=<password>

DIR=cal_data_${TARGET}
cal:
	python3 tart_cal.py --api ${API} --data ${DIR}
	
	
image:
	python3 image_from_raw.py --file cal_raw_data.hdf


data:
	rm -f ${DIR}/*
	python3 get_cal_data.py --api ${API} --pw ${PW} --n 1 --interval 5 --dir ${DIR}
