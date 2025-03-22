#!/bin/bash
while [ 0 ]
do
        TARGET=bw-biust \
        TART_CAL_ARGS="--get-gains --corr-only" \
        TART_GET_DATA=1 \
        TART_NCAL=3 \
        TART_CAL_INT=22 \
        TART_CAL_ITERATIONS=1000 \
        TART_CAL_ELEVATION=40 \
        TART_CAL_POINTING=0 \
        TART_CAL_POINTING_RANGE=3 \
        TART_CAL_UPLOAD=1 \
        TART_LOGIN_PW=sharkbait ./do_raw_cal.sh
done
