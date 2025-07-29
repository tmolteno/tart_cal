import logging
import re

import numpy as np
from tart.imaging import correlator, elaz
from tart.operation import observation
from tart_tools import api_imaging

logger = logging.getLogger()


def check_source(src, el_deg, jy_limit):
    print(f"check_source({src}, {el_deg})")
    if 'LUCH' in src['name']:
        return False
    if 'EGNOS' in src['name']:
        return False
    if 'EUTELSAT' in src['name']:
        return False
    if 'GAGAN' in src['name']:
        return False
    if 'IGSO' in src['name']:
        return False
    if 'BEIDOU-2 G' in src['name']:
        return False
    if 'BEIDOU-3 G' in src['name']:
        return False
    if src['el'] < el_deg:
        return False
    if src["jy"] < jy_limit:
        return False
    return True


def src_list_from_catalog_json(source_json, el_limit_deg=0.0, jy_limit=1e5):
    src_list = []
    for src in source_json:
        try:
            if check_source(src, el_deg=el_limit_deg, jy_limit=jy_limit):
                src_list.append(elaz.ElAz(src["el"], src["az"]))
        except Exception as e:
            print(f"ERROR in catalog src={src}")
            print(f"{e}")

    return src_list


def load_data_from_json(vis_json, src_json,
                        config, gains, phases, flag_list, el_threshold):
    logger.info("Loading data from JSON")
    logger.info(f"    vis_json = {vis_json['timestamp']}")
    cv, ts = api_imaging.vis_calibrated(vis_json, config,
                                        gains, phases, flag_list)
    logger.info(f"    ts = {ts}")
    _src_list = src_list_from_catalog_json(src_json, el_threshold)

    return cv, ts, _src_list


def get_prn_from_name(name):
    x = re.search(r"\(([A-Za-z0-9 /]*)\)", name)
    if x is None:
        print(f"PARSE ERROR: get_prn_from_name({name})")
        return None, None

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

    # Load the raw observations\

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
        if len(raw_files) > 0:
            raw_obs = [observation.Observation_Load(f) for f in raw_files]
            best_dt = 9e99
            logger.info(f"Vis timestamp {cv.get_timestamp()}")
            for o in raw_obs:
                dt = np.abs((o.timestamp - cv.get_timestamp()).total_seconds())
                logger.info(f"   raw obs.ts = {o.timestamp} dt={dt}")
                if dt < best_dt:
                    best_dt = dt
                    obs = o

            if (best_dt > 72):
                raise RuntimeError(f"Broken timestamps dt={best_dt} obs={obs.timestamp} vc={cv.get_timestamp()}")

            corr = correlator.Correlator()
            vis = corr.correlate(obs)
            logger.info(f"Timestamp: {vis.timestamp}")
            logger.info(f"Config: {vis.config.Dict}")

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
            logger.info(f"    Source: {int(d):02d}, stability: {std_ph:06.5f}, {np.mean(st):05.2f} {acq['sv']}")

        if n_obs == 0:
            raise RuntimeError(f"No satellites visible in obs[{i}]")

        logger.info("Good Satellites obs[{i}]")
        for s in good_satellites:
            logger.info(f"    {s}")

    n += n_obs

    best_acq = best_acq / n
    sv_noise = sv_noise / n
    s_n = best_acq / sv_noise
    logger.info(f"Best acw {best_acq}")
    logger.info(f"Noise level {s_n}")
    return good_satellites, best_acq


