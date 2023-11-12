#!/bin/bash
while [ 0 ]
do
    TARGET=signal \
        TART_CAL_ARGS="" \
        TART_CAL_INT=5 \
        TART_NCAL=3 \
        TART_CAL_ELEVATION=55 \
        TART_LOGIN_PW=sharkbait ./test.sh
done
