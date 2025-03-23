#!/usr/bin/env python
"""
    Calibrate the Telescope from the RESTful API

    Copyright (c) Tim Molteno 2017-2025.

    This tool uses  high-dimensional optimisation to calculate the gains and
    phases of the 24 antennas
    of the telescope.
"""
import matplotlib

import matplotlib.pyplot as plt

import argparse
import numpy as np
import os
import logging
from scipy import optimize
import json

from copy import deepcopy

import acquisition

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


matplotlib.use("agg")
logger = logging.getLogger()

triplets = None
ij_index = None
jk_index = None
ik_index = None


ift_scaled = None


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
        self.free_antennas = slice(1,self.nant)
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
        z = self.gains[self.free_antennas] * np.exp(self.phase_offsets[self.free_antennas] * 1j)
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
        # self.phase_indices=slice(1, self.nant)
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


def calc_score_aux(opt_parameters, measurements, window_deg, original_positions):
    global triplets, ij_index, jk_index, ik_index, masks, ift_scaled, myParam, mask_sums, full_sky_mask

    myParam.from_vector(opt_parameters)
    rot_rad = myParam.rot_rad
    gains = myParam.gains
    phase_offsets = myParam.phase_offsets

    ret_zone = 0.0
    ret_std = 0.0

    ant_idxs = np.arange(NANT)

    for i, m in enumerate(measurements):
        cv, ts, src_list, prn_list, obs = m

        cv.set_phase_offset(ant_idxs, phase_offsets)
        cv.set_gain(ant_idxs, gains)
        api_imaging.rotate_vis(np.degrees(rot_rad), cv, original_positions)

        n_bin = 2 ** 8
        cal_ift, cal_extent, n_fft, bin_width = api_imaging.image_from_calibrated_vis(
            cv, nw=n_bin / 4, num_bin=n_bin
        )

        abs_ift = np.abs(cal_ift)
        ift_std = np.median(abs_ift)
        ift_scaled = abs_ift / ift_std

        if full_sky_mask is None:
            x0, y0 = n_fft // 2, n_fft // 2
            d = (n_fft // 3)**2
            full_sky_mask = np.zeros_like(ift_scaled)
            for y in range(full_sky_mask.shape[0]):
                for x in range(full_sky_mask.shape[1]):
                    r2 = (y - y0)**2 + (x - x0)**2
                    p = np.exp(-(r2/d))
                    if p > 0.2:
                        p = 1
                    full_sky_mask[y, x] = p

        if masks[i] is None:
            print("Creating mask")
            mask = np.zeros_like(ift_scaled)

            for s in src_list:
                x0, y0 = s.get_px(n_fft)
                d = s.deg_to_pix(n_fft, window_deg)
                for y in range(mask.shape[0]):
                    for x in range(mask.shape[1]):
                        r2 = (y - y0)**2 + (x - x0)**2
                        p = np.exp(-(r2/d))
                        mask[y, x] += p

            # Zone outside of mask
            print(f"Max mask {np.max(mask)}, {np.min(mask)}")
            negative_mask = (-mask + 1)
            negative_mask[negative_mask < 0] = 0

            inv_masks[i] = negative_mask

            # mask[mask>0.4] = 1
            masks[i] = mask
            mask_sums[i] = np.sum(mask)

            plt.figure()
            plt.imshow(
                mask,
                extent=[-1, 1, -1, 1],
                vmin=0,
            )  # vmax=8
            plt.colorbar()
            plt.title(f"Mask {ts}")
            plt.xlim(1, -1)
            plt.ylim(-1, 1)
            plt.xlabel("East-West")
            plt.ylabel("North-South")
            plt.tight_layout()
            plt.savefig(f"{output_directory}/mask_{i}.png")
            plt.close()

        sn_score = (ift_scaled*full_sky_mask).max()  # Peak signal to noise.
        ret_std += -np.sqrt(sn_score)

        mask = masks[i]

        #
        # Clip the image at 0.5 maximum, and then calculated a score
        # showing how much clipped power there is in the known SV regions
        #
        # This is an attempt to avoid bright unknown sources from skewing
        # the phases towards it.
        masked_img = masks[i]*np.clip(ift_scaled, a_min=0, a_max=sn_score/2)
        in_zone = np.sum(np.sqrt(masked_img)) / mask_sums[i]
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
        ift_scaled,
        src_list,
        n_fft,
        bin_width,
        mask,
    )


def load_data_from_json(
    vis_json, src_json, config, gains, phases, flag_list, el_threshold
):
    logger.info("Loading data from JSON")
    logger.info(f"    vis_json = {vis_json['timestamp']}")
    cv, ts = api_imaging.vis_calibrated(vis_json, config,
                                        gains, phases, flag_list)
    logger.info(f"    ts = {ts}")
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
    global N_IT, method, output_directory, f_vs_iteration, ret_zone, ret_std

    ret_zone, ret_std, ift_scaled, src_list, n_fft, bin_width, mask = calc_score_aux(
        opt_parameters, measurements, window_deg, original_positions
    )
    ret = float(ret_zone + ret_std)
    if N_IT % 1000 == 0:
        # print(f"Iteration {N_IT}, score={ret:04.2f}")
        f_vs_iteration.append(ret)

    N_IT += 1
    return ret


class MyTakeStep(object):
    def __init__(self, stepsize, pointing_rad):
        self.stepsize = stepsize
        self.pointing_rad = pointing_rad

    def __call__(self, x):
        return myParam.take_step(x, self.stepsize)


def bh_callback(x, f, accepted):
    global output_directory, bh_basin_progress, N_IT, ift_scaled, masks, method
    global full_sky_mask, ret_zone, ret_std, in_zone
    # print(f"BH f={f:5.3f} accepted: {accepted}")
    if accepted:
        print(f"   S/N {ret_std:04.2f}, ZONE: {ret_zone:04.2f}")
        myParam.from_vector(x)
        myParam.output()
        bh_basin_progress.append([N_IT, float(f)])
        with open(f"{output_directory}/bh_basin_progress.json", "w") as fp:
            json.dump(bh_basin_progress, fp, indent=4, separators=(",", ": "))

        with open(f"{output_directory}/BH_basin_{f:5.3f}_{N_IT}.json", "w") as fp:
            myParam.output(fp)

        mask = masks[-1]
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
        plt.savefig(f"{output_directory}/mask_{f:5.3f}_{N_IT:05d}.png")
        plt.close()

        plt.figure()
        plt.imshow(
            ift_scaled*full_sky_mask,
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
        plt.savefig(f"{output_directory}/full_sky_mask_{f:5.3f}_{N_IT:05d}.png")
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
        plt.title(f"f={f:4.2f} r={np.degrees(myParam.rot_rad):4.2f}")
        plt.scatter(x_list, y_list, c="red", s=5)
        plt.xlabel("East-West")
        plt.ylabel("North-South")
        plt.tight_layout()
        plt.savefig(f"{output_directory}/{method}_{f:5.3f}_accepted_{N_IT:05d}.png")
        plt.close()


def de_callback(xk, convergence):
    print("DE at {} conv={}".format(xk, convergence))
    output_param(xk)


import glob

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
    with open(json_files[0], "r") as json_file:
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
        raise Exception("ERROR: Cannot Have both cold_start and get-gains specified")

    if ARGS.cold_start:
        gains = np.ones(len(gains_json["gain"]))
        phase_offsets = np.zeros(len(gains_json["phase_offset"]))

    config = settings.from_api_json(info["info"], ant_pos)

    pointing_error = np.radians(ARGS.pointing_range)
    pointing_center = np.radians(ARGS.pointing)

    # Now deal with the measurements.

    # Load the raw observations
    raw_obs = []
    for raw_file in raw_files:
        raw_obs.append(observation.Observation_Load(raw_file))

    masks = []
    full_sky_mask = None
    mask_sums = []
    inv_masks = []
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
        masks.append(None)
        inv_masks.append(None)
        mask_sums.append(None)

        if len(measurements) >= ARGS.num_meas:
            break

    # Acquisition to get expected list of SV's

    try:
        with open(f"{data_dir}/gps_acquisition.json", "r") as fp:
            full_acquisition_data = json.load(fp)
            calculate = False
    except:
        calculate = True

    if calculate:
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
                        print(f"    Antenna {i}")
                        ant_i = obs.get_antenna(i)
                        mean_i = np.mean(ant_i)

                        raw_data = ant_i - mean_i

                        num_samples_per_ms = sampling_freq // 1000
                        num_samples = int(2*num_samples_per_ms)

                        [prn, strength, phase, freq] = acquisition.acquire_full(
                            raw_data,
                            sampling_freq=sampling_freq,
                            center_freq=4.092e6, searchBand=4000, PRN=prn_i,
                            debug=True)

                        strengths.append(float(np.round(strength, 3)))
                        phases.append(float(np.round(phase, 3)))
                        freqs.append(float(np.round(freq, 1)))

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

    best_acq = np.zeros(NANT)
    sv_noise = np.zeros(NANT) + 0
    n = 0
    best_score = -999

    good_satellites = []

    for acquisition_data in full_acquisition_data:
        for d in acquisition_data:
            acq = acquisition_data[d]
            ph = np.array(acq['phases'])
            st = np.array(acq['strengths'])

            sv = acq['sv']

            if sv['el'] < 15:
                sv_noise += st

            mean_str = np.median(st)
            if (mean_str > 7.0):
                good_satellites.append(sv)

                good = np.where(st > 7.0)

                best_acq[good] += st[good]
                n = n + 1
            print(f"    Source: {int(d):02d}, stability: {np.std(ph):06.5f}, {np.mean(st):05.2f} {acq['sv']}")

    if n == 0:
        raise RuntimeError("No satellites visible")

    print("Good Satellites")
    for s in good_satellites:
        print(f"    {s}")

    best_acq = best_acq / n
    sv_noise = sv_noise / n
    s_n = best_acq / sv_noise
    print(f"Best acw {best_acq}")
    print(f"Noise level {s_n}")
    test_gains = best_acq / best_acq[0]
    test_gains = 1.0 / (test_gains)
    # These factors would make all SV appear equal brightness.
    # test_gains = np.ones(NANT)
    test_gains[best_acq < 0.1] = 0
    print(f"Estimated gains: {test_gains}")

    # Now remove satellites from the catalog that we can't see.
    # https://github.com/JasonNg91/GNSS-SDR-Python/tree/master/gnsstools

    good_src_list = []
    if ARGS.corr_only:
        for m in measurements:
            cv, ts, src_list, prn_list, obs = m
            print("Source List:")
            for s in src_list:
                good = False
                for g in good_satellites:
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

    src_list = good_src_list
    print("Source List")
    for s in src_list:
        print(f"    {s}")

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

    window_deg = np.degrees(1.2*0.2/blmax)
    print(f"window_deg: {window_deg}")

    s = calc_score(
        myParam.to_vector(),
        config,
        measurements,
        window_deg,
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
        window_deg,
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
        with open("{}/bh_basin_progress.json".format(output_directory), "w") as fp:
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
