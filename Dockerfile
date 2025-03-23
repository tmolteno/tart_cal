# Create a docker file for the calibration process.. 
# NOTE. This script executes the calibration process only once.
#
#  Author. Tim Molteno. tim@elec.ac.nz (c) 2018-2022.
#
FROM debian:bookworm
LABEL Maintainer="Tim Molteno tim@elec.ac.nz"
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y &&  apt-get install -y \
    python3-venv

RUN apt-get clean -y
RUN rm -rf /var/lib/apt/lists/*

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv --system-site-packages $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install tart python packages
RUN pip3 install --no-cache-dir --no-compile poetry

WORKDIR /app
COPY raw_cal raw_cal
COPY README.md .
COPY pyproject.toml .
COPY poetry.lock .
RUN ls -rl
RUN poetry install
# setup working directory
WORKDIR /work

# Add the calibrate script.
ADD tart_calibrate.sh /tart_calibrate.sh
ADD raw_calibrate.sh /raw_calibrate.sh

# CMD sh /tart_calibrate.sh
