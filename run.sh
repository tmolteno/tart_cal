#!/bin/sh
# Run the docker image for the telescope calibration cron job (run every two hours)
# Bind Mount the directory ~/calibration_results as persistent storage for the results of calibration
#
# Author. Tim Molteno tim@elec.ac.nz. (c) 2018.
# TODO. Pass the tart login password as the first command line parameter. e.g.
#    sh run.sh <passwd>
#
# Make sure the output directory exists.
mkdir -p ~/calibration_results
# Now run the container

# For production use:
# TART_NCAL = 6 # number of measurements for calibration
# TART_CAL_INT = 30 # minutes between measurements

docker run --rm \
    -e TART_LOGIN_PW=$2 \
    -e TART_API=https://tart.elec.ac.nz/$1/ \
    -e TART_NCAL=2 \
    -e TART_CAL_INT=10 \
    -v ~/calibration_results:/work \
    --name=cal ghcr.io/tmolteno/tart_cal /tart_calibrate.sh
