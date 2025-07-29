#!/bin/bash
while [ 0 ]
do
        TARGET=mu-udm \
        TART_CAL_ARGS="--get-gains --corr-only" \
        TART_GET_DATA=1 \
        TART_NCAL=3 \
        TART_CAL_INT=10 \
        TART_CAL_ITERATIONS=200 \
        TART_CAL_ELEVATION=45 \
        TART_CAL_POINTING=0 \
        TART_CAL_POINTING_RANGE=1 \
        TART_CAL_UPLOAD=1 \
        TART_LOGIN_PW=sharkbait ./do_raw_cal.sh
done
