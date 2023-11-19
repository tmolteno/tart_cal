#!/bin/bash
while [ 0 ]
do
    TARGET=rhodes \
        TART_CAL_ARGS="--get-gains" \
        TART_CAL_INT=13 \
        TART_NCAL=2 \
        TART_CAL_ITERATIONS=1000 \
        TART_CAL_ELEVATION=50 \
        TART_CAL_POINTING_RANGE=5 \
        TART_LOGIN_PW=sharkbait ./test.sh
done
