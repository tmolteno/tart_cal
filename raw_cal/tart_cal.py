#!/usr/bin/env python
"""
    Calibrate the Telescope from the RESTful API

    Copyright (c) Tim Molteno 2017-2025.

    This tool uses  high-dimensional optimisation to calculate the gains and
    phases of the 24 antennas
    of the telescope.
"""
import argparse
import glob
import json
import logging
import os
from copy import deepcopy

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from tart.imaging import elaz, imaging
from tart.operation import settings

from .mask_image import add_source as mask_add_source

from .calibration_data import find_good_satellites, load_cal_files
from .pos_from_gps import get_gnss_data

matplotlib.use("agg")
logger = logging.getLogger()

triplets = None
ij_index = None
jk_index = None
ik_index = None


NANT = 24


class Param:
    def __init__(self, nant, pointing_center, pointing_error):
        self.nant = nant
        self.gains = None
        self.phase_offsets = None
        self.rot_rad = None
        self.pointing_error = pointing_error
        self.pointing_center = pointing_center

    def take_step(self, stepsize):
        pass

    def from_vector(self, x):
        raise RuntimeError("Implement in subclass")

    def to_json(self):
        ret = {
            "gain": np.round(self.gains, 4).tolist(),
            "rot_degrees": float(np.degrees(self.rot_rad)),
            "phase_offset": np.round(self.phase_offsets, 4).tolist(),
        }
        return ret

    def output(self, fp=None):
        ret = self.to_json()
        if fp is None:
            print(json.dumps(ret, indent=4, separators=(",", ": ")))
        else:
            json.dump(ret, fp, indent=4, separators=(",", ": "))

    def pointing_bounds(self):
        return (self.pointing_center-self.pointing_error, self.pointing_center + self.pointing_error)


class ParamReIm(Param):

    def __init__(self, nant, pointing_center, pointing_error):
        super().__init__(nant, pointing_center, pointing_error)
        self.nend = int(2*self.nant-1)
        self.free_antennas = slice(1, self.nant)
        self.re_indices = slice(1, self.nant)
        self.im_indices = slice(self.nant, self.nend)

    def from_vector(self, x):
        self.rot_rad = x[0]
        re = np.concatenate(([1], x[self.re_indices]))
        im = np.concatenate(([0], x[self.im_indices]))
        self.gains = np.sqrt(re * re + im * im)
        self.phase_offsets = np.arctan2(im, re)

    def to_vector(self):
        ret = np.zeros(self.nend)
        ret[0] = self.rot_rad
        z = self.gains[self.free_antennas] * \
            np.exp(self.phase_offsets[self.free_antennas] * 1j)
        ret[self.re_indices] = z.real
        ret[self.im_indices] = z.imag
        return ret

    def bounds(self, test_gains):

        bounds = np.empty(self.nend, dtype=(float,2))

        self.blim = np.zeros(self.nant-1)
        self.blim[:] = test_gains[1:]
        print(f"Blim: {self.blim}")
        bounds[0] = self.pointing_bounds()
        bounds[self.re_indices] = (-2,2)
        bounds[self.im_indices] = (-2,2)

        for i in range(1, self.nant):
            tg = test_gains[i]
            if tg < 0.01:
                blim = 0
            else:
                blim = 2  # tg * 1.2

            bounds[i] = (-blim, blim)
            bounds[i + self.nant-1] = (-blim, blim)

        self.bounds = bounds
        return bounds

    def take_step(self, x, stepsize):
        pnt = self.pointing_error*stepsize

        bounded_step = stepsize*np.ones(self.nant-1)*self.blim

        rot_step = np.random.uniform(-pnt, pnt)
        re_step = np.random.uniform(-bounded_step, bounded_step)   # Re
        im_step = np.random.uniform(-bounded_step, bounded_step)   # Im

        ret = x + np.concatenate(([rot_step], re_step, im_step))

        return ret


class ParamPhase(Param):

    def __init__(self, nant, pointing_center, pointing_error, gains):
        super().__init__(nant, pointing_center, pointing_error)
        self.gains = gains
        self.nend = nant
        self.free_antennas = slice(1, self.nant)
        self.phase_indices = slice(1, self.nant)
        self.phase_offsets = np.zeros_like(self.gains)

    def take_step(self, x, stepsize):
        ret = np.zeros_like(x)

        pnt = self.pointing_error*stepsize

        phase_step = stepsize * 2*np.pi

        rot_step = np.random.normal(0, pnt)
        phase_steps = np.random.normal(0, phase_step, self.nant-1)

        new_rot = x[0] + rot_step
        new_phase = x[self.phase_indices] + phase_steps
        new_phase[new_phase < -np.pi] += 2*np.pi
        new_phase[new_phase > np.pi] -= 2*np.pi
        ret = np.concatenate(([new_rot], new_phase))
        return ret

    def from_vector(self, x):
        self.rot_rad = x[0]
        self.phase_offsets[self.phase_indices] = x[self.phase_indices]

    def to_vector(self):
        ret = np.zeros(self.nend)
        ret[0] = self.rot_rad
        ret[self.phase_indices] = self.phase_offsets[self.free_antennas]
        return ret

    def bounds(self, test_gains):
        bounds = [0] * self.nend
        bounds[0] = self.pointing_bounds()
        for i in range(1, self.nant):
            tg = test_gains[i]
            if tg < 0.01:
                bounds[i] = (0, 0)
            else:
                bounds[i] = (-np.pi*2, np.pi*2)  # Bounds for phases
                # bounds[i] = (-np.inf, np.inf)  # Bounds for phases

        return bounds


class ParamGainPhase(Param):

    def __init__(self, nant, pointing_center, pointing_error, gains):
        super().__init__(nant, pointing_center, pointing_error)
        self.gains = gains
        self.nend = int(2*self.nant-1)
        self.free_antennas = slice(1, self.nant)
        self.gain_indices = slice(1, self.nant)
        self.phase_indices = slice(self.nant, self.nend)
        self.phase_offsets = np.zeros_like(self.gains)

    def take_step(self, x, stepsize):
        ret = np.zeros_like(x)

        pnt = self.pointing_error*stepsize

        phase_step = stepsize * np.pi
        gain_step = stepsize * 0.1

        rot_step    = np.random.normal(0, pnt)
        phase_steps = np.random.normal(0, phase_step, self.nant-1)
        gain_steps  = np.random.normal(0, gain_step, self.nant-1)

        new_rot = x[0] + rot_step

        new_gains = x[self.gain_indices] + gain_steps
        new_phase = x[self.phase_indices] + phase_steps
        new_phase[new_phase < -np.pi] += 2*np.pi
        new_phase[new_phase > np.pi] -= 2*np.pi
        ret = np.concatenate(([new_rot], new_gains, new_phase))
        return ret

    def from_vector(self, x):
        self.rot_rad = x[0]
        self.phase_offsets[self.free_antennas] = x[self.phase_indices]
        self.gains[self.free_antennas] = x[self.gain_indices]

    def to_vector(self):
        ret = np.zeros(self.nend)
        ret[0] = self.rot_rad
        ret[self.phase_indices] = self.phase_offsets[self.free_antennas]
        ret[self.gain_indices] = self.gains[self.free_antennas]
        return ret

    def bounds(self, test_gains):
        bounds = [0] * self.nend
        bounds[0] = self.pointing_bounds()
        for i in range(1, self.nant):
            tg = test_gains[i]
            if tg < 0.01:
                bounds[i] = (0, 1e-3)  # Gain bounds
                bounds[i + self.nant - 1] = (0, 1e-3)  # Phase bounds
            else:
                bounds_err = 0.3
                bounds[i] = (test_gains[i] * (1-bounds_err), test_gains[i] * 1+(bounds_err)) # Bounds for gains
                bounds[i + self.nant - 1] = (-np.inf, np.inf) # Bounds for phases

        return bounds


def calc_score_aux(opt_parameters, measurements, window_radius_deg, original_positions):
    global triplets, ij_index, jk_index, ik_index, mask_list, myParam
    global mask_sums, ret_std, ret_zone, full_sky_masks
    global ift_scaled_list

    myParam.from_vector(opt_parameters)
    rot_rad = myParam.rot_rad
    gains = myParam.gains
    phase_offsets = myParam.phase_offsets

    ret_zone = 0.0
    ret_std = 0.0

    ant_idxs = np.arange(NANT)
    ift_scaled_list = []

    for i, m in enumerate(measurements):
        cv, ts, src_list, prn_list, obs = m

        cv.set_phase_offset(ant_idxs, phase_offsets)
        cv.set_gain(ant_idxs, gains)
        imaging.rotate_vis(np.degrees(rot_rad), cv, original_positions)

        n_bin = 2 ** 7
        cal_ift, cal_extent, n_fft, bin_width = \
            imaging.image_from_calibrated_vis(cv, nw=n_bin / 4,
                                              num_bin=n_bin)

        abs_ift = np.abs(cal_ift)
        ift_std = np.median(abs_ift)
        ift_scaled = abs_ift / ift_std
        ift_scaled_list.append(ift_scaled)

        if full_sky_masks[i] is None:
            logger.info(f"Creating full sky mask {ift_scaled.shape}")
            full_sky_mask = np.zeros_like(ift_scaled)
            src = elaz.ElAz(90, 0)
            mask_add_source(full_sky_mask, src, radius_deg=90)
            full_sky_masks[i] = full_sky_mask

        if mask_list[i] is None:
            logger.info(f"Creating mask {i}")
            mask = np.zeros_like(ift_scaled)

            for s in good_source_lists[i]:
                mask_add_source(mask, s, radius_deg=window_radius_deg)

            # Zone outside of mask
            print(f"Max mask {i} {np.max(mask)}, {np.min(mask)}")
            negative_mask = (-mask + 1)
            negative_mask[negative_mask < 0] = 0

            mask_list[i] = mask
            inv_masks[i] = negative_mask

            # mask[mask>0.4] = 1
            mask_list[i] = mask
            mask_sums[i] = np.sum(mask) + 0.01

            plt.figure()
            plt.imshow(
                mask,
                extent=[-1, 1, -1, 1],
                vmin=0,
            )  # vmax=8
            plt.colorbar()
            plt.title(f"Mask {i} {ts}")
            plt.xlim(1, -1)
            plt.ylim(-1, 1)
            plt.xlabel("East-West")
            plt.ylabel("North-South")
            plt.tight_layout()
            plt.savefig(f"{output_directory}/mask_{i}.png")
            plt.close()

        masked_sky = (ift_scaled*full_sky_masks[i])
        max_sky = np.sqrt(masked_sky.max())
        sn_score = max_sky / np.std(masked_sky)  # Peak signal to noise
        ret_std += -sn_score

        #
        # Clip the image at 0.5 maximum, and then calculated a score
        # showing how much clipped power there is in the known SV regions
        #
        # This is an attempt to avoid bright unknown sources from skewing
        # the phases towards it.
        masked_img = mask_list[i]*ift_scaled
        in_zone = np.sum((masked_img)) / mask_sums[i]
        # outmask_img = inv_masks[i]*ift_scaled
        # out_zone = np.median(outmask_img)

        zone_score = (in_zone)**3

        ret_zone += -zone_score

    ret_std = ret_std / len(measurements)
    ret_zone = ret_zone / len(measurements)

    # if N_IT % 100 == 0:
    #     print(f"S/N {ret_std:04.2f}, ZONE: {ret_zone:04.2f}, in: {in_zone:04.2f}", end='\r')

    return (
        ret_zone, ret_std,
        ift_scaled_list,
        n_fft,
        bin_width
    )


def calc_score(
    opt_parameters,
    config,
    measurements,
    window_radius_deg,
    original_positions,
    update=False,
    show=False,
):
    global N_IT, method, output_directory, f_vs_iteration, ret_zone, ret_std

    ret_zone, ret_std, ift_scaled_list, n_fft, bin_width = calc_score_aux(
        opt_parameters, measurements, window_radius_deg, original_positions
    )
    ret = float(ret_zone + ret_std)
    if N_IT % 1000 == 0:
        # print(f"Iteration {N_IT}, score={ret:04.2f}")
        f_vs_iteration.append(ret)

    N_IT += 1
    return ret


class MyTakeStep:
    def __init__(self, stepsize, pointing_rad):
        self.stepsize = stepsize
        self.pointing_rad = pointing_rad

    def __call__(self, x):
        global myParam
        return myParam.take_step(x, self.stepsize)


def bh_callback(x, f, accepted):
    global output_directory, bh_basin_progress, N_IT, ift_scaled_list, mask_list
    global ret_zone, ret_std, in_zone, myParam
    global good_source_lists, method
    # print(f"BH f={f:5.3f} accepted: {accepted}")
    if accepted:
        print(f"   S/N {ret_std:04.2f}, ZONE: {ret_zone:04.2f}")
        myParam.from_vector(x)
        myParam.output()
        bh_basin_progress.append([N_IT, float(f)])
        with open(f"{output_directory}/bh_basin_progress.json", "w") as fp:
            json.dump(bh_basin_progress, fp, indent=4, separators=(",", ": "))

        fname = f"BH_basin_{f:5.3f}_{N_IT}.json"
        with open(os.path.join(output_directory, fname), "w") as fp:
            myParam.output(fp)

        for i, mask in enumerate(mask_list):
            ift_sel = ift_scaled_list[i] * mask
            x_list, y_list = elaz.get_source_coordinates(good_source_lists[i])

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
            fname = f"mask_{-f:5.3f}_{N_IT:05d}_{i}.png"
            plt.savefig(os.path.join(output_directory, fname))
            plt.close()

            plt.figure()
            plt.imshow(
                ift_scaled_list[i]*full_sky_masks[i],
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
            fname = f"full_sky_mask_{-f:5.3f}_{N_IT:05d}_{i}.png"
            plt.savefig(os.path.join(output_directory, fname))
            plt.close()

            plt.figure()
            plt.imshow(
                ift_scaled_list[i],
                extent=[-1, 1, -1, 1],
                vmin=0,
            )  # vmax=8
            plt.colorbar()
            plt.xlim(1, -1)
            plt.ylim(-1, 1)
            plt.title(f"f={f:4.2f} i={i} r={np.degrees(myParam.rot_rad):4.2f}")
            plt.scatter(x_list, y_list, c="red", s=5)
            plt.xlabel("East-West")
            plt.ylabel("North-South")
            plt.tight_layout()
            fname = f"{method}_{f:5.3f}_accepted_{N_IT:05d}_{i}.png"
            plt.savefig(os.path.join(output_directory, fname))
            plt.close()


def de_callback(xk, convergence):
    print(f"DE at {xk} conv={convergence}")
    output_param(xk)


myParam = None
mask_list = None
full_sky_masks = None
mask_sums = []
inv_masks = []
measurements = []
good_source_lists = []
output_directory = None
N_IT = None
f_vs_iteration = None
bh_basin_progress = None
ift_scaled_list = None
method = None
ret_zone = None
ret_std = None
in_zone = None


def main():
    global myParam, output_directory, N_IT, f_vs_iteration
    global bh_basin_progress, method
    global full_sky_masks, mask_list, mask_sums, inv_masks

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

    parser = argparse.ArgumentParser(
        description="Calibrate the tart telescope from downloaded data.",
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

    parser.add_argument("--show", action="store_true", help="show instead of save.")
    parser.add_argument(
        "--cold-start", action="store_true", help="Start from zero knowledge of gains."
    )
    parser.add_argument(
        "--get-gains",
        action="store_true",
        help="Start from current knowledge of gains.",
    )
    parser.add_argument(
        "--phases",
        action="store_true",
        help="Use fixed gains from GPS and modify phases.",
    )
    parser.add_argument(
        "--gains-phases",
        action="store_true",
        help="Use nearly fixed gains from GPS and modify gains and phases.",
    )

    parser.add_argument(
        "--corr-only",
        action="store_true",
        help="Use only satellites that we can correlate against for calibration.",
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
        "--num-meas",
        required=False,
        type=int,
        default=999,
        help="Number of measurements to use",
    )

    parser.add_argument(
        "--elevation", type=float, default=30.0, help="Elevation threshold for sources")

    parser.add_argument(
        "--pointing", type=float, default=0.0, help="Initial estimate of pointing offset (degrees)")

    parser.add_argument(
        "--pointing-range", type=float, default=3.0, help="Pointing search range (degrees)")

    parser.add_argument(
        '--ignore', nargs='+', type=int, default=[], help="Specify the list of antennas to zero out.")

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

    flag_list = ARGS.ignore
    logger.info(f"Flagging {flag_list}")

    method = ARGS.method
    output_directory = ARGS.dir
    os.makedirs(output_directory, exist_ok=True)
    f_vs_iteration = []

    original_positions = deepcopy(config.get_antenna_positions())

    gains_json = calib_info["gains"]

    g_json = gains_json["gain"]
    logger.info(f"Gains JSON: {g_json}")
    gains = np.asarray(g_json)

    phase_offsets = np.asarray(gains_json["phase_offset"])
    logger.info(f"Gains {gains}")
    logger.info(f"Phases {phase_offsets}")

    if ARGS.cold_start and ARGS.get_gains:
        raise Exception("ERROR: Cannot Have both cold-start and get-gains specified")

    config = settings.from_api_json(info["info"], ant_pos)

    pointing_error = np.radians(ARGS.pointing_range)
    pointing_center = np.radians(ARGS.pointing)

    # Now deal with the measurements.
    measurements = load_cal_files(raw_files, calib_info,
                                  config, elevation_threshold=ARGS.elevation)
    inv_masks = [None] * len(measurements)
    mask_sums = [None] * len(measurements)
    mask_list = [None] * len(measurements)
    full_sky_masks = [None] * len(measurements)

    # Acquisition to get expected list of SV's
    fname = f"{data_dir}/gps_acquisition.json"
    full_acquisition_data = get_gnss_data(measurements, fname)

    # Use the standard deviation of the phases to determine whether the SV is visible.
    print("Finding visible satellites")
    good_satellites, best_acq = find_good_satellites(full_acquisition_data)

    # Find the best SV and record its index.
    best_sv = np.argmax(best_acq)

    test_gains = best_acq / best_acq[best_sv]
    print(f"test_gains = {test_gains}")
    test_gains[test_gains < 0.1] = 0.1

    test_gains = 1.0 / (test_gains)

    test_gains[test_gains > 3] = 3

    # These factors would make all SV appear equal brightness.
    # test_gains = np.ones(NANT)
    print(f"Estimated gains: {test_gains}")

    # Now remove satellites from the catalog that we can't see.
    # https://github.com/JasonNg91/GNSS-SDR-Python/tree/master/gnsstools

    for i, m in enumerate(measurements):
        good_src_list = []
        cv, ts, src_list, prn_list, obs = m
        if ARGS.corr_only:
            print(f"Source List[{i}]:")
            for s in src_list:
                good = False
                for g in good_satellites[i]:
                    el_match = (np.abs(np.degrees(s.el_r) - g['el']) < 0.1)
                    az_match = (np.abs(np.degrees(s.az_r) - g['az']) < 0.1)
                    if el_match and az_match:
                        good = True
                        good_src_list.append(s)
                        break

                print(f"    {s} {good}")
        else:
            # Remove satellites below elevation
            for s in src_list:
                if np.degrees(s.el_r) > ARGS.elevation:
                    good_src_list.append(s)

        good_source_lists.append(good_src_list)
        print(f"Source List [{i}]")
        for s in good_src_list:
            print(f"    {s}")

    if ARGS.cold_start:
        test_gains = np.ones(len(gains_json["gain"]))
        phase_offsets = np.zeros(len(gains_json["phase_offset"]))

    if ARGS.phases:
        print("Using PHASES")
        myParam = ParamPhase(NANT, pointing_center, pointing_error, test_gains)
        myParam.rot_rad = pointing_center
        myParam.phase_offsets = phase_offsets
        bh_stepsize = 1
        bh_T = 1.0

    elif ARGS.gains_phases:
        print("Using GAINS_PHASES")
        myParam = ParamGainPhase(NANT, pointing_center,
                                 pointing_error, test_gains)
        myParam.rot_rad = pointing_center
        myParam.gains = test_gains
        myParam.phase_offsets = phase_offsets
        bh_stepsize = 0.2
        bh_T = 0.4

    else:
        print("Using REIM")
        myParam = ParamReIm(NANT, pointing_center, pointing_error)
        myParam.rot_rad = pointing_center
        myParam.gains = gains
        myParam.phase_offsets = phase_offsets
        bh_stepsize = 1
        bh_T = 0.5

    bounds = myParam.bounds(test_gains)
    print(f"Bounds {bounds}")

    myParam.output()

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

    window_diameter_deg = np.degrees(1.2*0.2/blmax)
    window_radius_deg = window_diameter_deg / 2.0
    print(f"window_radius_deg: {window_radius_deg}")

    s = calc_score(
        myParam.to_vector(),
        config,
        measurements,
        window_radius_deg,
        original_positions,
        update=False,
        show=False,
    )
    s = float(s)

    print(f"Score from initial parameters = {s}")
    # bh_T = np.abs(s/40)
    # print(f"Basinhopping T = {bh_T}")

    f = lambda param: calc_score(
        param,
        config,
        measurements,
        window_radius_deg,
        original_positions,
        update=False,
        show=False,
    )

    init_parameters = myParam.to_vector()

    zero_list = ARGS.ignore

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
            "tol": 1e-3,
            "options": {"maxcor": 48},
        }
        mytakestep = MyTakeStep(bh_stepsize, pointing_error)

        ret = optimize.basinhopping(
            f,
            init_parameters,
            niter=ARGS.iterations,
            T=bh_T,
            take_step=mytakestep,
            disp=True,
            minimizer_kwargs=minimizer_kwargs,
            callback=bh_callback,
        )
        with open(f"{output_directory}/bh_basin_progress.json", "w") as fp:
            json.dump(bh_basin_progress, fp, indent=4, separators=(",", ": "))

    print(f"Success = {ret}")
    print(f"Final Score = {ret.fun}")

    myParam.from_vector(ret.x)
    rot_rad = myParam.rot_rad
    output_json = myParam.to_json()
    output_json["message"] = ret.message
    output_json["optimum"] = float(ret.fun)
    output_json["iterations"] = ret.nit

    new_positions = settings.rotate_location(
        np.degrees(rot_rad), np.array(original_positions).T
    )
    pos_list = (np.array(new_positions).T).tolist()
    output_json["antenna_positions"] = pos_list
    print(output_json)
    json_fname = f"{output_directory}/{method}_opt_json.json"
    with open(json_fname, "w") as fp:
        json.dump(output_json, fp, indent=4, separators=(",", ": "))
    print(f"Optimal solution: rot_degrees={output_json['rot_degrees']}")
    print(f"Optimal solution written to file {json_fname}")

    f_history_json = {}
    f_history_json["start"] = s
    f_history_json["history"] = f_vs_iteration

    hist_fname = f"{output_directory}/{method}_history.json"
    with open(hist_fname, "w") as fp:
        json.dump(f_history_json, fp, indent=4, separators=(",", ": "))
    print(f"Basin history written to file {hist_fname}")

    print(f"tart_upload_gains --api {ARGS.api} --gains {json_fname} --pw xxxx")
