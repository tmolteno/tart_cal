
import argparse
import glob
import json

import numpy as np
from tart.operation import settings
from tqdm import tqdm

from .acquisition import acquire_full
from .calibration_data import load_cal_files


def acquire_all_gnss(measurements):
    full_acquisition_data = {}
    full_acquisition_data['antennas'] = None
    full_acquisition_data['satellites'] = []

    for m in measurements:
        cv, ts, src_list, prn_list, obs = m

        acquisition_data = {}
        if full_acquisition_data['antennas'] is None:
            full_acquisition_data['antennas'] = \
                [i for i in range(cv.get_config().get_num_antenna())]

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
                for i in tqdm(range(num_antenna)):
                    ant_i = obs.get_antenna(i)
                    mean_i = np.mean(ant_i)

                    raw_data = ant_i - mean_i

                    num_samples_per_ms = sampling_freq // 1000
                    num_samples = int(2*num_samples_per_ms)
                    [prn, strength, phase, freq] = \
                        acquire_full(raw_data[0:num_samples],
                                     sampling_freq=sampling_freq,
                                     center_freq=4.092e6, searchBand=4000,
                                     PRN=prn_i, debug=False)

                    strengths.append(float(np.round(strength, 3)))
                    phases.append(float(phase))
                    freqs.append(float(np.round(freq, 1)))

                acquisition_data[f"{prn_i}"]['strengths'] = strengths
                acquisition_data[f"{prn_i}"]['phases'] = phases
                acquisition_data[f"{prn_i}"]['freqs'] = freqs
                acquisition_data[f"{prn_i}"]['sv'] = sv

                print(acquisition_data[f"{prn_i}"])
        full_acquisition_data['satellites'].append(acquisition_data)
    return full_acquisition_data


def get_gnss_data(measurements, fname):
    try:
        with open(fname) as fp:
            full_acquisition_data = json.load(fp)
    except:
        full_acquisition_data = acquire_all_gnss(measurements)
        with open(fname, "w") as fp:
            json.dump(full_acquisition_data, fp,
                      indent=4, separators=(",", ": "))
    return full_acquisition_data


def find_vis_satellites(full_acquisition_data, elevation):

    # Use the standard deviation of the phases to determine whether the SV is visible.
    print("Finding visible satellites")
    print(full_acquisition_data)
    antennas = full_acquisition_data['antennas']
    nant = len(antennas)
    best_acq = np.zeros(nant)  # Keep track of number of sv's visible
    n = 0
    best_score = -999

    equations = {'antennas': antennas}
    equations['sv'] = []

    for acquisition_data in full_acquisition_data['satellites']:
        for d in acquisition_data:
            acq = acquisition_data[d]
            ph = np.array(acq['phases'])
            st = np.array(acq['strengths'])

            sv = acq['sv']

            mean_str = np.median(st)
            med_ph = np.median(ph)
            good = np.where(st > 7.5)

            best_acq[good] += 1

            eqn = {}
            eqn['sv'] = sv
            eqn['cond'] = []
            for i in np.array(antennas)[good]:
                eqn['cond'].append({'codephase': ph[i],  'antenna': int(i), 'correlation': st[i]})

            #print(f"Source: {int(d):02d}, stability: {np.std(ph[good]):06.5f}, {np.mean(st):05.2f} {acq['sv']}")
            equations['sv'].append(eqn)

    if np.sum(best_acq) == 0:
        raise RuntimeError("No satellites visible")

    equations['best'] = best_acq.tolist()

    return equations


def main():
    global ARGS

    parser = argparse.ArgumentParser(
        description="Locate positions of the tart telescope from SV data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--api",
        required=False,
        default="https://tart.elec.ac.nz/signal",
        help="Telescope API server URL.",
    )
    parser.add_argument(
        "--data",
        required=False,
        default="cal_data",
        help="Calibration Input Data Directory.",
    )

    parser.add_argument("--dir", required=False, default=".", help="Output directory.")

    parser.add_argument(
        "--iterations",
        required=False,
        type=int,
        default=300,
        help="Number of iterations for basinhopping",
    )

    parser.add_argument(
        "--num-meas",
        required=False,
        type=int,
        default=999,
        help="Number of measurements to use",
    )

    parser.add_argument(
        "--elevation", type=float, default=30.0, help="Elevation threshold for sources")

    ARGS = parser.parse_args()

    # Load calibration data from the data directory

    data_dir = ARGS.data
    json_files = [f for f in glob.glob(f"{data_dir}/cal*.json")]
    raw_files = [f for f in glob.glob(f"{data_dir}/*.hdf")]

    print(json_files)
    with open(json_files[0]) as json_file:
        calib_info = json.load(json_file)

    info = calib_info["info"]
    ant_pos = calib_info["ant_pos"]
    config = settings.from_api_json(info["info"], ant_pos)

    measurements = load_cal_files(raw_files, calib_info,
                                  config, elevation_threshold=ARGS.elevation)

    fname = f"{data_dir}/gps_acquisition.json"

    full_acquisition_data = get_gnss_data(measurements, fname)

    eqn = find_vis_satellites(full_acquisition_data, ARGS.elevation)

    # Now solve for the relative positions
    print(json.dumps(eqn, indent=4))

    # Find an antenna that is present in all data.

    for b, a in zip(eqn['best'], eqn['antennas'], strict=False):
        if b < 3:
            print(a, b)

    zero_antenna = np.argmax(eqn['best'])
    print(f"Zero Antenna: {zero_antenna}")

    # Build matrix
    for sveqn in eqn['sv']:
        sv = sveqn['sv']
        print(sv)
        el = np.radians(sv['el'])
        az = np.radians(sv['az'])
        zero_phase = None
        cond = sveqn['cond']
        for c in cond:
            if (c['antenna'] == zero_antenna):
                zero_phase = c['codephase']
                zero_str = c['correlation']
                print(f"     {c}")

        if (zero_phase is None):
            break
        for c in cond:
            if (c['antenna'] != zero_antenna):
                dph = c['codephase'] - zero_phase
                weight = c['correlation']*zero_str
                phase_to_meters = 300000  # distance light travels in one ms
                dx = dph*np.cos(el)*np.sin(az)*phase_to_meters
                dy = dph*np.cos(el)*np.cos(az)*phase_to_meters
                print(f"     d[{zero_antenna},{c['antenna']}] {dph}, {dx},{dy} {weight}")


if __name__ == "__main__":
    main()
