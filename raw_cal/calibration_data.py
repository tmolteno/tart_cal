import logging
import re

import numpy as np

from tart_tools import api_imaging
from tart.imaging import elaz

from tart.operation import observation
from tart.imaging import correlator

logger = logging.getLogger()


def load_data_from_json(vis_json, src_json,
                        config, gains, phases, flag_list, el_threshold):
    logger.info("Loading data from JSON")
    logger.info(f"    vis_json = {vis_json['timestamp']}")
    cv, ts = api_imaging.vis_calibrated(vis_json, config,
                                        gains, phases, flag_list)
    logger.info(f"    ts = {ts}")
    _src_list = elaz.from_json(src_json, el_threshold)
    return cv, ts, _src_list


def get_prn_from_name(name):
    x = re.search(r"\(([A-Za-z0-9 /]*)\)", name)
    prn_txt = x.group(1)

    # Now (PRN 123)
    x = re.search(r"[A-Za-z0-9 /]*PRN ([0-9]+)", prn_txt)
    if x is not None:
        prn_num_txt = x.group(1)
        return int(prn_num_txt), 'GPS'

    # Now (PRN E123)
    x = re.search(r"[A-Za-z0-9 /]*PRN (E[0-9]+)", prn_txt)
    if x is not None:
        prn_num_txt = x.group(1)
        print(f"    prn_num_txt={prn_num_txt}")
        return int(prn_num_txt), 'GALILEO'

    print(f"  SV: {name}")
    return None, None


def load_cal_files(raw_files,
                   calib_info,
                   config,
                   elevation_threshold):

    # Load the raw observations
    raw_obs = []
    for raw_file in raw_files:
        raw_obs.append(observation.Observation_Load(raw_file))

    n = config.get_num_antenna()

    gains = np.ones(n)
    phase_offsets = np.zeros(n)
    flag_list = []

    measurements = []
    for d in calib_info["data"]:

        vis_json, src_json = d
        cv, ts, src_list = load_data_from_json(
            vis_json,
            src_json,
            config,
            gains,
            phase_offsets,
            flag_list,
            el_threshold=elevation_threshold
        )

        prn_list = []
        for sv in src_json:
            try:
                prn, gnss_type = get_prn_from_name(sv['name'])
                if gnss_type == 'GPS':
                    prn_list.append((prn, sv))
            except:
                print(prn)

        # Find best raw observation, with the closest timestamp
        obs = None
        best_dt = 9e99
        print(f"Vis timestamp {cv.get_timestamp()}")
        for o in raw_obs:
            dt = np.abs((o.timestamp - cv.get_timestamp()).total_seconds())
            print(f"   raw obs.ts = {o.timestamp} dt={dt}")
            if dt < best_dt:
                best_dt = dt
                obs = o

        if (best_dt > 72):
            raise RuntimeError(f"Broken timestamps dt={best_dt} obs={obs.timestamp} vc={cv.get_timestamp()}")

        corr = correlator.Correlator()
        vis = corr.correlate(obs)
        print(f"Timestamp: {vis.timestamp}")
        print(f"Config: {vis.config.Dict}")

        measurements.append([cv, ts, src_list, prn_list, obs])

    return measurements


def find_good_satellites(full_acquisition_data):

    NANT = 24
    best_acq = np.zeros(NANT)
    sv_noise = np.zeros(NANT) + 0

    best_score = -999

    good_satellites = []
    sv_noise = 0
    n = 0
    for i, acquisition_data in enumerate(full_acquisition_data['satellites']):
        n_obs = 0
        good_satellites.append([])
        for d in acquisition_data:
            acq = acquisition_data[d]
            ph = np.array(acq['phases'])
            st = np.array(acq['strengths'])

            sv = acq['sv']

            if sv['el'] < 15:
                sv_noise += st

            # mean_str = np.median(st)
            good_ants = np.where(st > 7)
            std_ph = np.std(ph[good_ants])

            if (np.mean(st) > 7) and (std_ph < 0.01):
                good_satellites[i].append(sv)

                good = np.where(st > 7.0)

                best_acq[good] += st[good]
                n_obs += 1
            print(f"    Source: {int(d):02d}, stability: {std_ph:06.5f}, {np.mean(st):05.2f} {acq['sv']}")

        if n_obs == 0:
            raise RuntimeError(f"No satellites visible in obs[{i}]")

        print("Good Satellites obs[{i}]")
        for s in good_satellites:
            print(f"    {s}")

    n += n_obs

    best_acq = best_acq / n
    sv_noise = sv_noise / n
    s_n = best_acq / sv_noise
    print(f"Best acw {best_acq}")
    print(f"Noise level {s_n}")
    return good_satellites, best_acq


def load_raw_data(raw_files):
    raw_obs = []
    for raw_file in raw_files:
        raw_obs.append(observation.Observation_Load(raw_file))

    for d in calib_info["data"]:
        vis_json, src_json = d
        cv, ts, src_list = load_data_from_json(
            vis_json,
            src_json,
            config,
            gains,
            phase_offsets,
            flag_list,
            el_threshold=ARGS.elevation,
        )
        print(f"ts = {ts}")
        prn_list = []
        for sv in src_json:
            prn = sv['name'].split('PRN ')
            if len(prn) < 2:
                continue

            prn = prn[1].split(')')[0]

            try:
                prn_list.append((int(prn), sv))
            except:
                logger.info(f"Ignoring: {prn}")

        # Find best raw observation, with the closest timestamp
        obs = None
        best_dt = 9e99
        logger.info(f"Vis timestamp {cv.get_timestamp()}")
        for o in raw_obs:
            dt = np.abs((o.timestamp - cv.get_timestamp()).total_seconds())
            logger.info(f"   raw obs.ts = {o.timestamp} dt={dt}")
            if dt < best_dt:
                best_dt = dt
                obs = o

        if (best_dt > 100):
            raise RuntimeError(f"Broken timestamps dt={best_dt} obs={obs.timestamp} vc={cv.get_timestamp()}")

        corr = correlator.Correlator()
        vis = corr.correlate(obs)
        logger.info(f"Timestamp: {vis.timestamp}")
        logger.info(f"Config: {vis.config.Dict}")

        measurements.append([cv, ts, src_list, prn_list, obs])
        inv_masks.append(None)
        mask_sums.append(None)

        if len(measurements) >= ARGS.num_meas:
            break
