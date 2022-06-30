#!/usr/bin/env python
"""
    Calibrate the Telescope from the RESTful API

    Copyright (c) Tim Molteno 2017-2022.
    
    This tool uses  high-dimensional optimisation to calculate the gains and phases of the 24 antennas
    of the telescope.
"""
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt

import argparse
import numpy as np
import time
import os

from copy import deepcopy

import itertools

from tart.operation import settings
from tart.operation import observation

from tart.imaging import visibility
from tart.imaging import calibration
from tart.imaging import synthesis
from tart.imaging import elaz
from tart.imaging import correlator

from tart.util import constants
from tart.util.angle import from_rad

from tart_tools import api_imaging
from tart_tools import api_handler

triplets = None
ij_index = None
jk_index = None
ik_index = None

from acquisition import acquire

REIM = False
def split_param(x):
    rot_rad = x[0]
    if REIM:
        re = np.concatenate(([1], x[1:24]))
        im = np.concatenate(([0], x[24:47]))
        gains = np.sqrt(re * re + im * im)
        phase_offsets = np.arctan2(im, re)
    else:
        gains = np.concatenate(([1], x[1:24]))
        phase_offsets = np.concatenate(([0], x[24:47]))

    return rot_rad, gains, phase_offsets


def join_param(rot_rad, gains, phase_offsets):
    ret = np.zeros(47)
    ret[0] = rot_rad
    if REIM:
        z = gains[1:24] * np.exp(phase_offsets[1:24] * 1j)
        ret[1:24] = z.real
        ret[24:47] = z.imag
    else:
        ret[1:24] = gains[1:24]
        ret[24:47] = phase_offsets[1:24]
    return ret


def param_to_json(x):
    rot_rad, gains, phase_offsets = split_param(x)
    ret = {
        "gain": np.round(gains, 4).tolist(),
        "rot_degrees": np.degrees(rot_rad),
        "phase_offset": np.round(phase_offsets, 4).tolist(),
    }
    return ret


def output_param(x, fp=None):
    ret = param_to_json(x)
    if fp is None:
        print(json.dumps(ret, indent=4, separators=(",", ": ")))
    else:
        json.dump(ret, fp, indent=4, separators=(",", ": "))


def calc_score_aux(opt_parameters, measurements, window_deg, original_positions):
    global triplets, ij_index, jk_index, ik_index, masks
    rot_rad, gains, phase_offsets = split_param(opt_parameters)

    ret_zone = 0.0
    ret_std = 0.0

    ant_idxs = np.arange(24)

    for i, m in enumerate(measurements):
        cv, ts, src_list, prn_list, obs = m

        cv.set_phase_offset(ant_idxs, phase_offsets)
        cv.set_gain(ant_idxs, gains)
        api_imaging.rotate_vis(np.degrees(rot_rad), cv, original_positions)

        n_bin = 2 ** 7
        cal_ift, cal_extent, n_fft, bin_width = api_imaging.image_from_calibrated_vis(
            cv, nw=n_bin / 4, num_bin=n_bin
        )

        abs_ift = np.abs(cal_ift)
        ift_std = np.std(abs_ift)
        ift_scaled = abs_ift / ift_std

        ret_std += -np.sqrt(ift_scaled.max())  # Peak signal to noise.

        if masks[i] is None:
            print("Creating mask")
            mask = np.zeros_like(ift_scaled)

            for s in src_list:
                x0,y0 = s.get_px(n_fft)
                d = 2*s.deg_to_pix(n_fft, window_deg)
                for y in range(mask.shape[0]):
                    for x in range(mask.shape[1]):
                        r2 = (y - y0)**2 + (x - x0)**2
                        p = np.exp(-r2/d)
                        mask[y, x] += p
            masks[i] = mask / np.sum(mask)

        mask = masks[i]
        zone_score = -np.sum(mask * ift_scaled)

        ret_zone += 4 * zone_score

    if N_IT % 100 == 0:
        print(f"S/N {ret_std:04.2f}, ZONE: {ret_zone:04.2f}")

    return (
        (ret_std + ret_zone) / len(measurements),
        ift_scaled,
        src_list,
        n_fft,
        bin_width,
        mask,
    )



def load_data_from_json(
    vis_json, src_json, config, gains, phases, flag_list, el_threshold
):

    cv, ts = api_imaging.vis_calibrated(vis_json, config, gains, phases, flag_list)
    src_list = elaz.from_json(src_json, el_threshold)
    return cv, ts, src_list



def calc_score(
    opt_parameters,
    config,
    measurements,
    window_deg,
    original_positions,
    update=False,
    show=False,
):
    global N_IT, method, output_directory, f_vs_iteration

    ret, ift_scaled, src_list, n_fft, bin_width, mask = calc_score_aux(
        opt_parameters, measurements, window_deg, original_positions
    )

    if N_IT % 100 == 0:
        print(f"Iteration {N_IT}, score={ret:04.2f}")
        f_vs_iteration.append(ret)

    if N_IT % 1000 == 0:

        #ift_sel = np.zeros_like(ift_scaled)

        #for s in src_list:
            #x_min, x_max, y_min, y_max, area = s.get_px_window(n_fft, window_deg)
            #ift_sel[y_min:y_max, x_min:x_max] = ift_scaled[y_min:y_max, x_min:x_max]

        ift_sel = ift_scaled*mask
        x_list, y_list = elaz.get_source_coordinates(src_list)

        plt.figure()
        plt.imshow(
            ift_sel,
            extent=[-1, 1, -1, 1],
            vmin=0,
        )  # vmax=8
        plt.colorbar()
        plt.xlim(1, -1)
        plt.ylim(-1, 1)
        plt.scatter(x_list, y_list, c="red", s=5)
        plt.xlabel("East-West")
        plt.ylabel("North-South")
        plt.tight_layout()
        plt.savefig("{}/opt_slice_{:05d}.png".format(output_directory, N_IT))
        plt.close()

        plt.figure()
        plt.imshow(
            ift_scaled,
            extent=[-1, 1, -1, 1],
            vmin=0,
        )  # vmax=8
        plt.colorbar()
        plt.xlim(1, -1)
        plt.ylim(-1, 1)
        plt.title(ret)
        plt.scatter(x_list, y_list, c="red", s=5)
        plt.xlabel("East-West")
        plt.ylabel("North-South")
        plt.tight_layout()
        plt.savefig(
            "{}/{}_{:5.3f}_opt_full_{:05d}.png".format(
                output_directory, method, ret, N_IT
            )
        )
        plt.close()
    N_IT += 1
    return ret

from scipy import optimize
import json


class MyTakeStep(object):
    def __init__(self, stepsize, pointing_rad):
        self.stepsize = stepsize
        self.pointing_rad = pointing_rad
        
    def __call__(self, x):
        s = self.stepsize
        
        pnt = self.pointing_rad*s
        x[0] += np.random.uniform(-pnt, pnt)
        x[1:] += np.random.uniform(-s, s, x[1:].shape)
        return x


def bh_callback(x, f, accepted):
    global output_directory, bh_basin_progress, N_IT
    print("BH f={} accepted {}".format(f, int(accepted)))
    output_param(x)
    if accepted:
        bh_basin_progress.append([N_IT, f])
        with open("{}/bh_basin_progress.json".format(output_directory), "w") as fp:
            json.dump(bh_basin_progress, fp, indent=4, separators=(",", ": "))

        with open(
            "{}/BH_basin_{:5.3f}_{}.json".format(output_directory, f, N_IT), "w"
        ) as fp:
            output_param(x, fp)


def de_callback(xk, convergence):
    print("DE at {} conv={}".format(xk, convergence))
    output_param(xk)

import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calibrate the tart telescope from API.",
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
        help="Calibration Output Data Directory.",
    )
    
    parser.add_argument("--show", action="store_true", help="show instead of save.")
    parser.add_argument(
        "--cold-start", action="store_true", help="Start from zero knowledge of gains."
    )
    parser.add_argument(
        "--get-gains",
        action="store_true",
        help="Start from current knowledge of gains.",
    )
    parser.add_argument("--dir", required=False, default=".", help="Output directory.")
    parser.add_argument(
        "--method",
        required=False,
        default="BH",
        help="Optimization Method [NM, LB, DE, BH]",
    )
    parser.add_argument(
        "--iterations",
        required=False,
        type=int,
        default=300,
        help="Number of iterations for basinhopping",
    )
    parser.add_argument(
        "--elevation", type=float, default=30.0, help="Elevation threshold for sources")
    
    parser.add_argument(
        "--pointing", type=float, default=0.0, help="Initial estimate of pointing offset (degrees)")

    parser.add_argument(
        '--ignore', nargs='+', type=int, default=[], help="Specify the list of antennas to zero out.")


    ARGS = parser.parse_args()

    # Load calibration data from the data directory
    
    data_dir = ARGS.data
    json_files = [f for f in glob.glob(f"{data_dir}/cal*.json")]
    raw_files = [f for f in glob.glob(f"{data_dir}/*.hdf")]

    print(json_files)
    with open(json_files[0], "r") as json_file:
        calib_info = json.load(json_file)
        
    info = calib_info["info"]
    ant_pos = calib_info["ant_pos"]
    config = settings.from_api_json(info["info"], ant_pos)
    
    flag_list = ARGS.ignore  # [4, 5, 14, 22] # optionally remove some unused antennas
    print(f"Flagging {flag_list}")

    method = ARGS.method
    output_directory = ARGS.dir
    os.makedirs(output_directory, exist_ok=True)
    f_vs_iteration = []

    original_positions = deepcopy(config.get_antenna_positions())

    gains_json = calib_info["gains"]
    print(gains_json["gain"])
    gains = np.asarray(gains_json["gain"])
    phase_offsets = np.asarray(gains_json["phase_offset"])
    print(gains)
    print(phase_offsets)

    if ARGS.cold_start and ARGS.get_gains:
        raise Exception("ERROR: Cannot Have both cold_start and get-gains specified")

    if ARGS.cold_start:
        gains = np.ones(len(gains_json["gain"]))
        phase_offsets = np.zeros(len(gains_json["phase_offset"]))

    config = settings.from_api_json(info["info"], ant_pos)

    init_parameters = join_param(0.0, gains, phase_offsets)
    output_param(init_parameters)

    masks = []
    measurements = []
    for d, raw_file in zip(calib_info["data"], raw_files):
            
        print(d)
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

        # Load the data here from the raw file
        obs = observation.Observation_Load(raw_file)
        corr = correlator.Correlator()
        vis = corr.correlate(obs)
        print(f"Timestamp: {vis.timestamp}")
        print(f"Config: {vis.config.Dict}")

        measurements.append([cv, ts, src_list, prn_list, obs])
        masks.append(None)

    # Acquisition to get expected list of SV's

    
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

    # Use the standard deviation of the phases to determine whether the SV is visible.
    print("Finding visible satellites")
    
    best_acq = -1
    best_score = -999
    for acquisition_data in full_acquisition_data:
        print(acquisition_data.keys())
        for d in acquisition_data:
            acq = acquisition_data[d]
            ph = np.array(acq['phases'])
            st = np.array(acq['strengths'])

            mean_str = np.median(st)
            score = np.max(mean_str - st)

            print(f"score={score}")
            if score > best_score:
                print(f"Best {d} score={score}")
                best_score = score
                best_acq = st / 6.7

            print(f"Source: {int(d):02d}, stability: {np.std(ph):06.5f}, {np.mean(st):05.2f} {acq['sv']}")

    # Now remove satellites from the catalog that we can't see.
    # https://github.com/JasonNg91/GNSS-SDR-Python/tree/master/gnsstools
        
    N_IT = 0
    window_deg = 4.0

    s = calc_score(
        init_parameters,
        config,
        measurements,
        window_deg,
        original_positions,
        update=False,
        show=False,
    )

    f = lambda param: calc_score(
        param,
        config,
        measurements,
        window_deg,
        original_positions,
        update=False,
        show=False,
    )

    pointing_error = np.radians(3)
    pointing_center = np.radians(ARGS.pointing)

    print(f"Calculating which antennas to ignore {best_acq}")
    test_gains = best_acq / best_acq[0]
    print(f"Test gains: {test_gains}")

    bounds = [0] * 47
    bounds[0] = (pointing_center-pointing_error, pointing_center + pointing_error)  # Bounds for the rotation parameter
    if REIM:
        for i in range(1,47):
            bounds[i] = (-2, 2) # Bounds for all other parameters (real and imaginary components)
    else:
        for i in range(1,24):
            tg = test_gains[i]
            bounds[i] = (max(0,tg - 0.1), tg + 0.1) # Bounds for all other parameters (real and imaginary components)
        for i in range(24,47):
            bounds[i] = (-np.pi, np.pi) # Bounds for all other parameters (real and imaginary components)




    zero_list = ARGS.ignore
    if zero_list is not None:
        print(f"Ignoring antennas {zero_list}")
        for i,a in enumerate(best_acq):
            if a < 0.5:
                print(a,i)
                bounds[i] = (0, 0.0001)

    print(f"Bounds {bounds}")
    np.random.seed(555)  # Seeded to allow replication.

    if method == "NM":
        ret = optimize.minimize(f, init_parameters, method="Nelder-Mead", tol=1e-5)
    if method == "LB":
        ret = optimize.minimize(f, init_parameters, method="L-BFGS-B", bounds=bounds)
    if method == "DE":
        ret = optimize.differential_evolution(f, bounds, disp=True)
    if method == "BH":
        bh_basin_progress = [[0, s]]
        minimizer_kwargs = {
            "method": "L-BFGS-B",
            "jac": False,
            "bounds": bounds,
            "tol": 1e-5,
            "options": {"maxcor": 48},
        }
        ret = optimize.basinhopping(
            f,
            init_parameters,
            niter=ARGS.iterations,
            T=0.5,
            take_step=MyTakeStep(1.0, pointing_error),
            disp=True,
            minimizer_kwargs=minimizer_kwargs,
            callback=bh_callback,
        )
        with open("{}/bh_basin_progress.json".format(output_directory), "w") as fp:
            json.dump(bh_basin_progress, fp, indent=4, separators=(",", ": "))

    rot_rad = ret.x[0]
    output_json = param_to_json(ret.x)
    output_json["message"] = ret.message
    output_json["optimum"] = ret.fun
    output_json["iterations"] = ret.nit

    new_positions = settings.rotate_location(
        np.degrees(rot_rad), np.array(original_positions).T
    )
    pos_list = (np.array(new_positions).T).tolist()
    output_json["antenna_positions"] = pos_list

    with open("{}/{}_opt_json.json".format(output_directory, method), "w") as fp:
        json.dump(output_json, fp, indent=4, separators=(",", ": "))

    f_history_json = {}
    f_history_json["start"] = s
    f_history_json["history"] = f_vs_iteration

    with open("{}/{}_history.json".format(output_directory, method), "w") as fp:
        json.dump(f_history_json, fp, indent=4, separators=(",", ": "))
