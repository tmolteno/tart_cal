TART=bw-biust

test:
	 TART_API=https://api.elec.ac.nz/tart/${TART} \
		TART_NCAL=3 TART_CAL_INT=11 \
		TART_CAL_ITER=700 \
		TART_LOGIN_PW=sharkbait \
		./tart_calibrate.sh

gps:
	TARGET=za-rhodes TART_LOGIN_PW=sharkbait ./gps_pos.sh

docker:
	docker build -t tart_cal --progress=plain .

testd:
	docker run -it tart_cal TART_API=https://api.elec.ac.nz/tart/${TART} \
		TART_NCAL=3 TART_CAL_INT=11 \
		TART_CAL_ITER=700 \
		TART_LOGIN_PW=sharkbait \
		./raw_calibrate.sh
