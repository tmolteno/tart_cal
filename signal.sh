#!/bin/bash
while [ 0 ]
do
        TARGET=signal \
        TART_CAL_ARGS="--get-gains --phases" \
        TART_CAL_INT=13 \
        TART_NCAL=3 \
        TART_CAL_ITERATIONS=1000 \
        TART_CAL_ELEVATION=45 \
        TART_CAL_POINTING_RANGE=3 \
        TART_CAL_UPLOAD=1 \
        TART_CAL_POINTING=0 \
        TART_LOGIN_PW=sharkbait ./do_raw_cal.sh
done
