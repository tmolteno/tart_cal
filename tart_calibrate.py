#!/usr/bin/env python
"""
    Calibrate the Telescope from the RESTful API

    Copyright (c) Tim Molteno 2017-2025.

    This tool uses  high-dimensional optimisation to calculate the gains and phases of the 24 antennas
    of the telescope.
"""
# import matplotlib
# matplotlib.use("agg")
import argparse
import json
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from tart.imaging import elaz
from tart.imaging import imaging
from tart.operation import settings
from tart_tools import api_handler, api_imaging

triplets = None
ij_index = None
jk_index = None
ik_index = None


def split_param(x):
    rot_degrees = x[0]
    re = np.concatenate(([1], x[1:24]))
    im = np.concatenate(([0], x[24:47]))
    gains = np.sqrt(re * re + im * im)
    phase_offsets = np.arctan2(im, re)
    return rot_degrees, gains, phase_offsets


def join_param(rot_degrees, gains, phase_offsets):
    ret = np.zeros(47)
    ret[0] = rot_degrees
    z = gains[1:24] * np.exp(phase_offsets[1:24] * 1j)
    ret[1:24] = z.real
    ret[24:47] = z.imag
    return ret


def param_to_json(x):
    rot_degrees, gains, phase_offsets = split_param(x)
    ret = {
        "gain": np.round(gains, 4).tolist(),
        "rot_degrees": float(rot_degrees),
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
    rot_degrees, gains, phase_offsets = split_param(opt_parameters)

    ret_zone = 0.0
    ret_std = 0.0

    ant_idxs = np.arange(24)

    for i, m in enumerate(measurements):
        cv, ts, src_list = m

        cv.set_phase_offset(ant_idxs, phase_offsets)
        cv.set_gain(ant_idxs, gains)
        api_imaging.rotate_vis(rot_degrees, cv, original_positions)

        n_bin = 2 ** 7
        cal_ift, cal_extent, n_fft, bin_width = \
            api_imaging.image_from_calibrated_vis(
                cv, nw=n_bin / 4, num_bin=n_bin)

        abs_ift = np.abs(cal_ift)
        ift_std = np.std(abs_ift)
        ift_scaled = abs_ift / ift_std

        ret_std += -np.sqrt(ift_scaled.max())  # Peak signal to noise.

        if masks[i] is None:
            print("Creating mask")
            mask = np.zeros_like(ift_scaled)

            for s in src_list:
                x0, y0 = s.get_px(n_fft)
                d = 2*imaging.deg_to_pix(n_fft, window_deg)
                for y in range(mask.shape[0]):
                    for x in range(mask.shape[1]):
                        r2 = (y - y0)**2 + (x - x0)**2
                        p = np.exp(-r2/d)
                        mask[y, x] += p
            masks[i] = mask / np.sum(mask)

        mask = masks[i]
        zone_score = -np.sum(mask * ift_scaled)
        #zones = []
        #for s in src_list:
            ## Get a pixel window around each source...
            #x_min, x_max, y_min, y_max, area = s.get_px_window(n_fft, window_deg)
            ## Integrate the brightness in a window around that source
            #s_px = ift_scaled[y_min:y_max, x_min:x_max]
            #zones.append(np.sum(s_px) / area)
        #zones = np.array(zones)
        #zone_score = -np.mean(zones)

        ret_zone += 4 * zone_score

    ret_score = float((ret_std + ret_zone) / len(measurements))
    # if N_IT % 100 == 0:
    #     print(f"S/N {ret_std:04.2f}, ZONE: {ret_zone:04.2f}", end='\r')

    return (
        ret_score,
        ift_scaled,
        src_list,
        n_fft,
        bin_width,
        mask,
    )


def load_data_from_json(vis_json, src_json,
                        config, gains, phases, flag_list,
                        el_threshold):
    cv, ts = api_imaging.vis_calibrated(vis_json, config,
                                        gains, phases, flag_list)
    src_list = elaz.from_json(src_json, el_threshold)
    return cv, ts, src_list


def get_window(x_i, y_i, d):
    x_min = int(np.floor(x_i - d))
    x_max = int(np.ceil(x_i + d))
    y_min = int(np.floor(y_i - d))
    y_max = int(np.ceil(y_i + d))
    return x_min, x_max, y_min, y_max


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
        # print(f"Iteration {N_IT}, score={ret:04.2f}")
        f_vs_iteration.append(ret)

    if N_IT % 1000 == 0:
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
        plt.savefig(f"{output_directory}/opt_slice_{N_IT:05d}.png")
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
            f"{output_directory}/{method}_{ret:5.3f}_opt_full_{N_IT:05d}.png"
        )
        plt.close()
    N_IT += 1
    return ret


class MyTakeStep:
    def __init__(self, stepsize=0.5):
        self.stepsize = stepsize

    def __call__(self, x):
        s = self.stepsize
        x[0] += np.random.uniform(-2.0 * s, 2.0 * s)
        x[1:] += np.random.uniform(-s, s, x[1:].shape)
        return x


def bh_callback(x, f, accepted):
    global output_directory, bh_basin_progress, N_IT
    print(f"BH f={f} accepted {int(accepted)}")
    if accepted:
        output_param(x)
        bh_basin_progress.append([N_IT, float(f)])
        fname = f"{output_directory}/bh_basin_progress.json"
        with open(fname, "w") as fp:
            json.dump(bh_basin_progress, fp, indent=4, separators=(",", ": "))

        fname = f"{output_directory}/BH_basin_{float(f):5.3f}_{N_IT}.json"
        with open(fname, "w") as fp:
            output_param(x, fp)


def de_callback(xk, convergence):
    print(f"DE at {xk} conv={convergence}")
    output_param(xk)


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
        "--file",
        required=False,
        default="calibration_data.json",
        help="Calibration Output JSON file.",
    )

    #parser.add_argument(
        #"--influx",
        #required=False,
        #default=None,
        #help="Use an influxdb  (https://tartinflux.max.ac.nz/) for cal data.",
    #)

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
        default="LB",
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
        "--elevation", type=float, default=40.0, help="Elevation threshold for sources]"
    )

    parser.add_argument(
        '--ignore', nargs='+', type=int, help="Specify the list of antennas to zero out.")


    ARGS = parser.parse_args()

    # Load calibration data
    with open(ARGS.file) as json_file:
        calib_info = json.load(json_file)

    info = calib_info["info"]
    ant_pos = calib_info["ant_pos"]
    config = settings.from_api_json(info["info"], ant_pos)

    flag_list = []  # [4, 5, 14, 22] # optionally remove some unused antennas

    method = ARGS.method
    output_directory = ARGS.dir
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
    if ARGS.get_gains:
        api = api_handler.APIhandler(ARGS.api)
        gains_json = api.get("calibration/gain")
        gains = np.asarray(gains_json["gain"])
        phase_offsets = np.asarray(gains_json["phase_offset"])

    config = settings.from_api_json(info["info"], ant_pos)

    init_parameters = join_param(0.0, gains, phase_offsets)
    output_param(init_parameters)

    masks = []
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
            el_threshold=ARGS.elevation,
        )
        measurements.append([cv, ts, src_list])
        masks.append(None)

    N_IT = 0

    # Rayleigh criterion 1.2 lambda / d_max
    blmax = 0
    ant_pos_arr = np.array(ant_pos)
    n_ant = ant_pos_arr.shape[0]
    for i in range(n_ant):
        for j in range(n_ant):
            dr = ant_pos_arr[i] - ant_pos_arr[j]
            r = np.sqrt(dr[0]*dr[0] + dr[1]*dr[1])
            if r > blmax:
                blmax = r

    window_deg = np.degrees(1.2*0.2/blmax)
    print(f"window_deg: {window_deg}")
    # window_deg = 5.0

    s = calc_score(
        init_parameters,
        config,
        measurements,
        window_deg,
        original_positions,
        update=False,
        show=False,
    )
    print(f"Score from initial parameters = {s}")

    f = lambda param: calc_score(
        param,
        config,
        measurements,
        window_deg,
        original_positions,
        update=False,
        show=False,
    )

    bounds = [(-3, 3)]  # Bounds for the rotation parameter
    for i in range(46):
        bounds.append((-1.2, 1.2))
        # Bounds for all other parameters (real and imaginary components)

    zero_list = ARGS.ignore
    if zero_list is not None:
        print(f"Ignoring antennas {zero_list}")
        for i in range(24):
            if i in zero_list:
                bounds[i] = (-0.01, 0.01)
                bounds[i+23] = (-0.01, 0.01)

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
            T=0.4,
            stepsize=0.5,
            disp=True,
            minimizer_kwargs=minimizer_kwargs,
            callback=bh_callback,
        )
        with open(f"{output_directory}/bh_basin_progress.json", "w") as fp:
            json.dump(bh_basin_progress, fp, indent=4, separators=(",", ": "))

    rot_degrees = ret.x[0]
    output_json = param_to_json(ret.x)
    output_json["message"] = ret.message
    output_json["optimum"] = ret.fun
    output_json["iterations"] = ret.nit

    new_positions = settings.rotate_location(
        rot_degrees, np.array(original_positions).T
    )
    print(new_positions)
    pos_list = (np.array(new_positions).T).tolist()
    print(pos_list)
    output_json["antenna_positions"] = pos_list

    with open(f"{output_directory}/{method}_opt_json.json", "w") as fp:
        json.dump(output_json, fp, indent=4, separators=(",", ": "))

    f_history_json = {}
    f_history_json["start"] = s
    f_history_json["history"] = f_vs_iteration

    with open(f"{output_directory}/{method}_history.json", "w") as fp:
        json.dump(f_history_json, fp, indent=4, separators=(",", ": "))
