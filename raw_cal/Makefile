API=https://tart.elec.ac.nz/signal
PW=<password>

cal:
	python3 tart_cal.py --raw cal_raw_data.hdf --file cal_data.json --elevation 20
	
	
image:
	python3 image_from_raw.py --file cal_raw_data.hdf


data:
	PW=${PW} sh get_cal_data.sh
