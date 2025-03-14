
test:
	 TART_API=https://api.elec.ac.nz/tart/za-rhodes \
		TART_NCAL=3 TART_CAL_INT=11 \
		TART_CAL_ITER=700 \
		TART_LOGIN_PW=sharkbait \
		./tart_calibrate.sh
