#!/bin/bash
while [ 0 ]
do
        TARGET=tart-kenya \
        TART_CAL_ARGS="" \
        TART_GET_DATA=1 \
        TART_NCAL=2 \
        TART_CAL_INT=7 \
        TART_CAL_ITERATIONS=800 \
        TART_CAL_ELEVATION=45 \
        TART_CAL_POINTING=0 \
        TART_CAL_POINTING_RANGE=2.5 \
        TART_CAL_UPLOAD=1 \
        TART_LOGIN_PW=nairobi ./do_raw_cal.sh
done
