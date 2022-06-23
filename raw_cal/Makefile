API=https://tart.elec.ac.nz/signal
PW=<password>

cal:
	python3 tart_cal.py --data cal_data
	
	
image:
	python3 image_from_raw.py --file cal_raw_data.hdf


data:
	PW=${PW} sh get_cal_data.sh
