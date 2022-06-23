API=https://tart.elec.ac.nz/signal
PW=<password>

cal:
	python3 tart_cal.py --raw cal_raw_data.hdf --file cal_data.json --elevation 20
	
	
image:
	python3 image_from_raw.py --file cal_raw_data.hdf

data:
	rm -f *.hdf
	tart_set_mode --raw --api ${API} --pw ${PW}
	sleep 60
	tart_download_data --raw --n 1
	tart_set_mode --vis --api ${API} --pw ${PW}
	tart_calibration_data --api ${API} --n 1 --file cal_data.json
	mv *.hdf cal_raw_data.hdf

