
#!/usr/bin/env python
#
# Get Calibration Data from the TART telescope
#
# Copyright (c) Tim Molteno 2017-2022.

import argparse
import time
import json
import logging
import os

import urllib.parse

from tart.operation import settings
from tart.imaging import visibility
from tart.imaging import calibration

from tart_tools import api_handler
from tart_tools import api_imaging

logger = logging.getLogger()

def split_param(x):
    rot_degrees = x[0]
    re = np.concatenate(([1], x[1:24]))
    im = np.concatenate(([0], x[24:47]))
    gains = np.sqrt(re * re + im * im)
    phase_offsets = np.arctan2(im, re)
    return rot_degrees, gains, phase_offsets


def load_data(api, config):
    logger.info(f"Loading new data from {api.root}")
    logger.info("Setting vis mode")
    resp = api.post_payload_with_token("mode/vis", {})
    vis_json = api.get("imaging/vis")
    ts = api_imaging.vis_json_timestamp(vis_json)
    logger.info(f"Getting catalog from {api.catalog_url(config, datestr=ts.isoformat())}")
    src_json = api.get_url(api.catalog_url(config, datestr=ts.isoformat()))
    logger.info("Loading Complete")
    return vis_json, src_json



def get_raw_data(api, config):
    global ARGS
    try:
        logger.info("Setting raw mode")
        resp = api.post_payload_with_token("mode/raw", {})
        logger.info("Sleeping 60 seconds to wait for raw data...")
        time.sleep(60)
        resp_raw = api.get("raw/data")
        entry = resp_raw[0]
        data_url = urllib.parse.urljoin(f"{api.root}/", entry["filename"])

        file_name = data_url.split("/")[-1]

        file_path = os.path.join(ARGS.dir, file_name)

        if os.path.isfile(file_path):
            if sha256_checksum(file_path) == entry["checksum"]:
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
        logger.info("Setting vis mode")
        resp = api.post_payload_with_token("mode/vis", {})

if __name__ == "__main__":

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
        "--api",
        required=False,
        default="https://tart.elec.ac.nz/signal",
        help="Telescope API server URL.",
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
    PARSER.add_argument("--pw", default="password", type=str, help="API password")
    PARSER.add_argument(
        "--n", type=int, default=5, help="Number of samples in the JSON data file."
    )
    PARSER.add_argument(
        "--interval",
        type=int,
        default=5,
        help="Interval (in minutes) between samples in the JSON data file.",
    )
    PARSER.add_argument(
        "--catalog",
        required=False,
        default="https://tart.elec.ac.nz/catalog",
        help="Catalog API URL.",
    )

    ARGS = PARSER.parse_args()

    os.makedirs(ARGS.dir, exist_ok=True)

    api = api_handler.AuthorizedAPIhandler(ARGS.api, ARGS.pw)

    info = api.get("info")
    ant_pos = api.get("imaging/antenna_positions")
    config = settings.from_api_json(info["info"], ant_pos)
    gains_json = api.get("calibration/gain")

    ret = {"info": info, "ant_pos": ant_pos, "gains": gains_json}
    data = []
    for i in range(ARGS.n):
        vis_json, src_json = load_data(api, config)  # Get Calibration Data
        get_raw_data(api, config)
        data.append([vis_json, src_json])
        if i != ARGS.n - 1:
            logger.info("Sleeping {} minutes".format(ARGS.interval))
            time.sleep(ARGS.interval * 60)
    ret["data"] = data
    
    
    with open(f"{ARGS.dir}/{ARGS.file}", "w") as fp:
        json.dump(ret, fp, indent=4, separators=(",", ": "))

