#!/bin/bash
while [ 0 ]
do
        TARGET=signal \
        TART_CAL_ARGS="--get-gains --corr-only" \
        TART_GET_DATA=1 \
        TART_NCAL=2 \
        TART_CAL_INT=13 \
        TART_CAL_ITERATIONS=800 \
        TART_CAL_ELEVATION=45 \
        TART_CAL_POINTING=0 \
        TART_CAL_POINTING_RANGE=3 \
        TART_CAL_UPLOAD=1 \
        TART_LOGIN_PW=sharkbait ./do_raw_cal.sh
done
