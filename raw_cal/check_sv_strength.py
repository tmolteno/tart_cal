import argparse
import json
import glob

import numpy as np

from tart.operation import settings
from tart.operation import observation
from tart.imaging import correlator

from tart_cal import load_data_from_json
from acquisition import acquire

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calibrate the tart telescope from downloaded data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--data",
        required=False,
        default="cal_data",
        help="Calibration Input Data Directory.",
    )

    ARGS = parser.parse_args()
    data_dir = ARGS.data
    json_files = [f for f in glob.glob(f"{data_dir}/cal*.json")]
    raw_files = [f for f in glob.glob(f"{data_dir}/*.hdf")]

    print(json_files)
    with open(json_files[0], "r") as json_file:
        calib_info = json.load(json_file)

    info = calib_info["info"]
    ant_pos = calib_info["ant_pos"]
    config = settings.from_api_json(info["info"], ant_pos)
    gains_json = calib_info["gains"]
    gains = np.asarray(gains_json["gain"])
    phase_offsets = np.asarray(gains_json["phase_offset"])
    flag_list = []
    # Now deal with the measurements.

    masks = []
    mask_sums = []
    inv_masks = []
    measurements = []
    for d, raw_file in zip(calib_info["data"], raw_files):

        vis_json, src_json = d
        cv, ts, src_list = load_data_from_json(
            vis_json,
            src_json,
            config,
            gains,
            phase_offsets,
            flag_list,
            el_threshold=0,
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
                pass
                #print(prn)

        # Load the data here from the raw file
        obs = observation.Observation_Load(raw_file)

        print(obs.timestamp , cv.get_timestamp())
        corr = correlator.Correlator()
        vis = corr.correlate(obs)
        print(f"Timestamp: {vis.timestamp}")
        print(f"Config: {vis.config.Dict}")

        measurements.append([cv, ts, src_list, prn_list, obs])
        masks.append(None)
        inv_masks.append(None)
        mask_sums.append(None)



    try:
        with open(f"{data_dir}/gps_acquisition.json", "r") as fp:
            full_acquisition_data = json.load(fp)
    except:
        full_acquisition_data = []
        for m in measurements:
            cv, ts, src_list, prn_list, obs = m

            acquisition_data = {}
            num_antenna = obs.config.get_num_antenna()
            sampling_freq = obs.get_sampling_rate()
            for svinfo in prn_list:
                prn_i, sv = svinfo
                if (prn_i <= 32):
                    acquisition_data[f"{prn_i}"] = {}
                    print(f"acquiring {svinfo}")
                    acquisition_data[f"{prn_i}"]['PRN'] = prn_i

                    strengths = []
                    phases = []
                    freqs = []
                    for i in range(num_antenna):
                        ant_i = obs.get_antenna(i)
                        mean_i = np.mean(ant_i)

                        raw_data = ant_i - mean_i

                        num_samples_per_ms = sampling_freq // 1000
                        num_samples = int(2*num_samples_per_ms)
                        [prn, strength, phase, freq] = acquire(raw_data[0:num_samples],
                                sampling_freq=sampling_freq,
                                center_freq=4.092e6, searchBand=6000, PRN=prn_i, debug=False)

                        strengths.append(strength)
                        phases.append(phase)
                        freqs.append(freq)

                    acquisition_data[f"{prn_i}"]['strengths'] = strengths
                    acquisition_data[f"{prn_i}"]['phases'] = phases
                    acquisition_data[f"{prn_i}"]['freqs'] = freqs
                    acquisition_data[f"{prn_i}"]['sv'] = sv


                    print(acquisition_data[f"{prn_i}"])
            full_acquisition_data.append(acquisition_data)

        with open(f"{data_dir}/gps_acquisition.json", "w") as fp:
            json.dump(full_acquisition_data, fp, indent=4, separators=(",", ": "))


    print("Finding visible satellites")

    best_acq = np.zeros(24)
    n = 0
    best_score = -999
    for acquisition_data in full_acquisition_data:
        for d in acquisition_data:
            acq = acquisition_data[d]
            ph = np.array(acq['phases'])
            st = np.array(acq['strengths'])

            [print(p) for p in ph]
            mean_str = np.median(st)

            if mean_str > 7.0:
                best_acq += st
                n = n + 1

            print(f"Source: {int(d):02d}, stability: {np.std(ph):06.5f}, {np.mean(st):05.2f} {acq['sv']}")

    if n == 0:
        raise RuntimeError("No satellites visible")

    best_acq = best_acq / n
    print(f"Calculating best SV signal strengths per antenna:")
    for i in range(24):
        print(f"   {i:2d} {best_acq[i]:4.2f}")

