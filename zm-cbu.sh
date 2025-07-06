#!/bin/bash
while [ 0 ]
do
        TARGET=zm-cbu \
        TART_CAL_ARGS="--get-gains --corr-only" \
        TART_GET_DATA=1 \
        TART_NCAL=2 \
        TART_CAL_INT=11 \
        TART_CAL_ITERATIONS=1000 \
        TART_CAL_ELEVATION=20 \
        TART_CAL_POINTING=0 \
        TART_CAL_POINTING_RANGE=1.5 \
        TART_CAL_UPLOAD=1 \
        TART_LOGIN_PW=sharkbait ./do_raw_cal.sh
done
