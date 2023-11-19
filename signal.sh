#!/bin/bash
while [ 0 ]
do
TARGET=signal \
        TART_CAL_ARGS="--get-gains" \
        TART_CAL_INT=13 \
        TART_NCAL=3 \
        TART_CAL_ITERATIONS=1000 \
        TART_CAL_ELEVATION=50 \
        TART_CAL_POINTING_RANGE=10 \
        TART_CAL_POINTING=0 \
        TART_LOGIN_PW=sharkbait ./test.sh
done
