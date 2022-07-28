#!/bin/bash
make data TARGET=rhodes PW=sharkbait LABEL=_A$1
python3 check_sv_strength.py --data cal_data_rhodes_A$1
