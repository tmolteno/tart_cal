import logging
import numpy as np

from tart_tools import api_imaging
from tart.imaging import elaz

from tart.operation import observation
from tart.imaging import correlator

logger = logging.getLogger()


def load_data_from_json(
    vis_json, src_json, config, gains, phases, flag_list, el_threshold
):
    logger.info("Loading data from JSON")
    logger.info(f"    vis_json = {vis_json['timestamp']}")
    cv, ts = api_imaging.vis_calibrated(vis_json, config,
                                        gains, phases, flag_list)
    logger.info(f"    ts = {ts}")
    _src_list = elaz.from_json(src_json, el_threshold)
    return cv, ts, _src_list


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
            prn = sv['name'].split('PRN ')
            if len(prn) < 2:
                continue

            prn = prn[1].split(')')[0]

            try:
                prn_list.append((int(prn), sv))
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
