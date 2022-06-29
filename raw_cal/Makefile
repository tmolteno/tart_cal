TARGET=signal
API=https://tart.elec.ac.nz/${TARGET}
PW=<password>
N=1
DIR=cal_data_${TARGET}
cal:
	python3 tart_cal.py --api ${API} --data ${DIR} --iterations 500 --ignore 14 --dir cal_out_${TARGET}

	
image:
	python3 image_from_raw.py --file cal_raw_data.hdf


data:
	rm -f ${DIR}/*
	python3 get_cal_data.py --api ${API} --pw ${PW} --n ${N} --interval 5 --dir ${DIR}
