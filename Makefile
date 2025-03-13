
test:
	 TART_API=https://api.elec.ac.nz/tart/tart-kenya \
		TART_NCAL=2 TART_CAL_INT=7 \
		TART_CAL_ITER=500 \
		TART_LOGIN_PW=sharkbait \
		./tart_calibrate.sh
