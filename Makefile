TART=nz-elec

test:
	 TART_API=https://api.elec.ac.nz/tart/${TART} \
		TART_NCAL=1 TART_CAL_INT=0 \
		TART_CAL_ITER=700 \
		TART_LOGIN_PW=sharkbait \
		./tart_calibrate.sh

gps:
	TARGET=za-rhodes TART_LOGIN_PW=sharkbait ./gps_pos.sh

docker:
	docker build -t tart_cal --progress=plain .

testd:
	mkdir -p work
	docker run --mount type=bind,source=./work,target=/work \
		-e TARGET=${TART} \
		-e TART_NCAL="1" \
		-e TART_CAL_INT="20" \
		-e TART_CAL_ITERATIONS=100 \
		-e TART_CAL_ELEVATION=20 \
		-e TART_CAL_POINTING=0 \
		-e TART_CAL_ARGS="--phases --corr-only" \
		-e TART_LOGIN_PW=sharkbait \
		-e TART_CAL_WORK_DIR="./work" \
		-e TART_UPLOAD="0" \
		 -i -t tart_cal /raw_calibrate.sh
