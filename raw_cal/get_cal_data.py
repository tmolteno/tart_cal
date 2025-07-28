#!/usr/bin/env python
#
# Get Calibration Data from the TART telescope
#
# Copyright (c) Tim Molteno 2017-2022.

import argparse
import json
import logging
import os
import time
import urllib.parse

import numpy as np

from tart.imaging import visibility
from tart.operation import settings
from tart.util import utc

from tart_tools import api_handler, api_imaging
from tart_tools.archive_handler import handle_archive_request

logger = logging.getLogger()


def set_vis_mode(api):
    try:
        logger.info("Setting vis mode")
        api.post_with_token("mode/vis")
    except Exception as e:
        logger.exception(e)
        logger.error("Error in setting vis_mode")


def get_catalogue_json(api, config, ts):
    cat_url = api.catalog_url(lon=config.get_lon(),
                              lat=config.get_lat(),
                              datestr=utc.to_string(ts))
    logger.info(f"Getting catalog from {cat_url}")
    return api.get_url(cat_url)


def load_data(api, config):
    logger.info(f"Loading new data from {api.root}")
    set_vis_mode(api)
    time.sleep(3)  # Wait for vis mode...
    vis_json = api.get("imaging/vis")

    logger.info(f"Vis Json timestamp {vis_json['timestamp']}")
    ts = api_imaging.vis_json_timestamp(vis_json)
    logger.info(f"Timestamp {ts}")
    logger.info(f"utcnow = {utc.now()}")
    logger.info(f"URL ts {utc.to_string(ts)}")

    src_json = get_catalogue_json(api, config, ts)
    logger.info("Loading Complete")
    return vis_json, src_json


def get_raw_data(api, config, vis_json, directory='.'):
    try:
        ts = api_imaging.vis_json_timestamp(vis_json)

        logger.info(f"Get raw data to match {ts}")
        logger.info("Setting acquire raw to 2**16 to match")
        resp = api.put("acquire/raw/num_samples_exp/16")

        best_dt = 999999.0
        while best_dt > 36:
            logger.info("Setting raw mode")
            resp = api.post_with_token("mode/raw")
            interval = 2
            logger.info(f"Sleeping {interval} seconds to wait for raw data...")
            time.sleep(interval)

            set_vis_mode(api)
            time.sleep(3)
            resp_raw = api.get("raw/data")

            # Find the entry closest to the timestamp
            best_dt = 999999
            entry = None
            for e in resp_raw:
                ts_string = e['timestamp']
                ts_string = ts_string.replace("+00Z", " GMT")  # Remove old +00Z format
                e_ts = utc.from_string(ts_string)

                dt = np.abs((e_ts - ts).total_seconds())
                logger.info(f"   raw obs.ts = {e_ts} dt={dt}")
                if dt < best_dt:
                    best_dt = dt
                    entry = e
                    logger.info(f"   best ts = {e_ts} dt={dt}")
        # entry = resp_raw[0]
        data_url = urllib.parse.urljoin(f"{api.root}/", entry["filename"])

        file_name = data_url.split("/")[-1]

        file_path = os.path.join(directory, file_name)

        if os.path.isfile(file_path):
            if api_handler.sha256_checksum(file_path) == entry["checksum"]:
                logger.info(f"Skipping {file_path}")
            else:
                logger.info(f"Corrupted File: {file_path}")
                os.remove(file_path)
                api_handler.download_file(data_url, entry["checksum"], file_path)
        else:
            api_handler.download_file(data_url, entry["checksum"], file_path)
    except Exception as e:
        logger.exception(e)
    finally:
        set_vis_mode(api)


def main():
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)

    PARSER = argparse.ArgumentParser(
        description="Generate Calibration Data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    PARSER.add_argument(
        "--target",
        required=True,
        help="Telescope TART name.",
    )
    PARSER.add_argument(
        "--file",
        required=False,
        default="calibration_data.json",
        help="Calibration Output JSON file.",
    )
    PARSER.add_argument(
        "--dir",
        required=False,
        default=".",
        help="Output file directory.",
    )
    PARSER.add_argument("--pw", default="password",
                        type=str, help="API password")
    PARSER.add_argument(
        "--n", type=int, default=5,
        help="Number of samples in the JSON data file."
    )
    PARSER.add_argument(
        "--interval",
        type=int,
        default=5,
        help="Interval (in minutes) between samples in the JSON data file."
    )
    PARSER.add_argument(
        "--catalog",
        required=False,
        default="https://tart.elec.ac.nz/catalog",
        help="Catalog API URL."
    )
    PARSER.add_argument('--raw', action='store_true', help="Download raw data")
    PARSER.add_argument('--archive', action='store_true', help="Get data from the archive")

    ARGS = PARSER.parse_args()

    os.makedirs(ARGS.dir, exist_ok=True)

    api_url = f"https://api.elec.ac.nz/tart/{ARGS.target}"
    api = api_handler.APIhandler(api_url)

    ret = {}
    if ARGS.archive:
        duration = ARGS.n * ARGS.interval
        handle_archive_request(target=ARGS.target,
                               num_observations=ARGS.n,
                               output_dir=ARGS.dir,
                               start_str=f"-{duration}",
                               duration_str=f"{duration}")

        info = api.get("info")
        ant_pos = api.get("imaging/antenna_positions")
        config = settings.from_api_json(info["info"], ant_pos)
        gains_json = api.get("calibration/gain")
        ret = {"info": info, "ant_pos": ant_pos, "gains": gains_json}

        data = []
        for i in range(ARGS.n):
            fname = f"obs_{i:05d}.hdf"
            fpath = os.path.join(ARGS.dir, fname)
            # Load some vis from the HDF file
            vis = visibility.from_hdf5(fpath)
            # ret["vis_list"] = vis_list
            # ret["config_json"] = config_json
            # ret["config"] = config
            # ret["ant_pos"] = ant_pos
            # ret["timestamps"] = timestamps
            # ret["baselines"] = hdf_baselines
            # ret["gain"] = h5f["gains"][:]
            # ret["phase_offset"] = h5f["phases"][:]
            if i == 0:
                ant_pos = vis['ant_pos']
                print(f"ant_pos = {ant_pos.shape}")
                # config = vis['config']
                # gains_json = {'gain': vis['gain'].tolist(),
                #               'phase_offset': vis['phase_offset'].tolist()}
                # ret["ant_pos"] = ant_pos.tolist(),
                # ret["gains"] = gains_json

            ts = vis['timestamps'][0]
            vis_json = vis["vis_list"][0].to_json()
            src_json = get_catalogue_json(api, config, ts)
            data.append([vis_json, src_json])
        ret["data"] = data
    else:

        info = api.get("info")
        ant_pos = api.get("imaging/antenna_positions")
        config = settings.from_api_json(info["info"], ant_pos)
        gains_json = api.get("calibration/gain")

        ret = {"info": info, "ant_pos": ant_pos, "gains": gains_json}

        data = []
        for i in range(ARGS.n):
            api = api_handler.AuthorizedAPIhandler(api_url, ARGS.pw)

            vis_json, src_json = load_data(api, config)  # Get Calibration Data
            if ARGS.raw:
                get_raw_data(api, config, vis_json, ARGS.dir)
            data.append([vis_json, src_json])
            if i != ARGS.n - 1:
                logger.info(f"Sleeping {ARGS.interval} minutes")
                time.sleep(ARGS.interval * 60)

        ret["data"] = data

    with open(f"{ARGS.dir}/{ARGS.file}", "w") as fp:
        json.dump(ret, fp, indent=4, separators=(",", ": "))
