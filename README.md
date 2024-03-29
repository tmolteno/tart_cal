## TART Telescope Calibration Software

Author: Tim Molteno

This docker file executes a process that uses known satellite positions to calibrate the telescope. The method of calibration is described in [1].

This instance requires a lot of processing power (approx 1 hour of CPU time on a single core high-end i7) and calibrates the TART telescope using a catalog of known L1-band objects (see object_position_server). These known objects are used as 'guide stars' to work out the gains and phases for each of the telescope antennas.

The process spends the first hour capturing 3 sets of observations at 25 minute intervals. After this, it uses these observations and the positions of the known radio sources at the time of observation. The entire process will take approximately two hours to complete (with only the second hour being very CPU intensive).

## Raw Cal : New Calibration using GPS-based gains

There are two scripts tart_calibrate.sh, and raw_calibrate.sh. The latter is the new one which is now recommended.

This container is automatically built and placed in the github container registry [https://ghcr.io/tmolteno/tart_cal]

    docker run --rm \
        -e TART_LOGIN_PW=replaceme \
        -e TART_API=https://tart.elec.ac.nz/signal/ \
        -e TARGET=signal \
        -e TART_NCAL=2 \
        -e TART_CAL_INT=20 \
        -v ~/calibration_results:/work \
        --name=cal ghcr.io/tmolteno/tart_cal /raw_calibrate.sh


#### Configuration

Configure environment variables in docker-compose.yml

#### Build and run
```bash
docker-compose build && docker-compose up -d
```

#### Follow progress
```bash
docker-compose logs --follow
```

####
```bash
docker-compose ps  # status is calibration still running?
docker-compose stop  # stop calibration run?
```



## Docker Usage

The easiest way to build this is to use docker. The continer is in the github container registry, so the following command will work (where "passwd" is the password to your telescope web api server):

    sh run.sh rhodes passwd

This will calibrate the telescope at the API endpoint 'rhodes' with password 'passwd'

### Stoping or Killing the instance

    docker ps -a
    docker stop cal
    docker rm cal


### Running this regularly

Add this command as a cron job. Modify the TART_LOGIN_PW and TART_API to refer to the URL of your TART telescope. Recommended interval is every two hours.

    docker run --rm \
        -e TART_LOGIN_PW=<your password here> \
        -e TART_API=https://tart.elec.ac.nz/signal/ \
        -e TART_NCAL=3 \
        -e TART_CAL_INT=30 \
        -v ~/calibration_results:/work \
        --name=cal ghcr.io/tmolteno/tart_cal /tart_calibrate.sh
    
### Debugging

This creates a running instance called 'cal'. You can check the logs using 

    docker attach cal

To exit type Ctrl-p Ctrl-q


To run a bash shell on the container. The logs are contained in the file /var/log/cron.log

    docker exec -it cal bash
    
Then you can manually run a calibration using

    sh /tart_calibrate.sh

The calibration normally takes 50 minutes to download data, and approximately 1.5 hours to run the optimization.

[1] Molteno et al. "Continuous Calibration of TART using GPS satellites". ENZCon2017.
