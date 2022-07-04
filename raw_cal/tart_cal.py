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

ift_scaled = None
#REIM = True
NANT=24


class Param:
    def __init__(self, nant, pointing_error):
        self.nant = nant
        self.gains = None
        self.phase_offsets = None
        self.rot_rad = None
        self.pointing_error = pointing_error

    def take_step(self, stepsize):
        pass

    def from_vector(self, x):
        raise RuntimeError("Implement in subclass")

    def to_json(self):
        ret = {
            "gain": np.round(self.gains, 4).tolist(),
            "rot_degrees": np.degrees(self.rot_rad),
            "phase_offset": np.round(self.phase_offsets, 4).tolist(),
        }
        return ret


    def output(self, fp=None):
        ret = self.to_json()
        if fp is None:
            print(json.dumps(ret, indent=4, separators=(",", ": ")))
        else:
            json.dump(ret, fp, indent=4, separators=(",", ": "))

    def pointing_bounds(self, pointing_center):
        return (pointing_center-self.pointing_error, pointing_center + self.pointing_error)


class ParamReIm(Param):

    def __init__(self, nant, pointing_error):
        super().__init__(nant, pointing_error)
        self.nend=int(2*self.nant-1)
        self.free_antennas=slice(1,self.nant)
        self.re_indices=slice(1, self.nant)
        self.im_indices=slice(self.nant, self.nend)

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

    def bounds(self, pointing_center, test_gains):
        bounds = np.empty(self.nend, dtype=(float,2))
        bounds[0] = self.pointing_bounds(pointing_center)
        bounds[self.re_indices] = (-2,2)
        bounds[self.im_indices] = (-2,2)

        for i in range(1,self.nant):
            tg = test_gains[i]
            if tg < 0.01:
                bounds[i] = (0, 1e-3)
                bounds[i + self.nant] = (0, 1e-3)

        return bounds

    def take_step(self, x, stepsize):
        pnt = self.pointing_error*stepsize

        rot_step = np.random.uniform(-pnt, pnt)
        re_step = np.random.uniform(-stepsize, stepsize, self.nant-1)   # Re
        im_step = np.random.uniform(-stepsize, stepsize, self.nant-1)   # Im

        ret = x + np.concatenate(([rot_step], re_step, im_step))

        return ret

class ParamPhase(Param):

    def __init__(self, nant, pointing_error, gains):
        super().__init__(nant, pointing_error)
        self.gains = gains
        self.nend = nant
        self.free_antennas=slice(1, self.nant)
        self.phase_indices=slice(1, self.nant)
        self.phase_offsets = np.zeros_like(self.gains)

    def take_step(self, x, stepsize):
        ret = np.zeros_like(x)

        pnt = self.pointing_error*stepsize

        phase_step = stepsize * np.pi

        rot_step = np.random.uniform(-pnt, pnt)
        phase_steps  = np.random.uniform(-phase_step, phase_step, self.nant-1)

        ret = x + np.concatenate(([rot_step], phase_steps))
        return ret

    def from_vector(self, x):
        self.rot_rad = x[0]
        self.phase_offsets[self.phase_indices] = x[self.phase_indices]

    def to_vector(self):
        ret = np.zeros(self.nend)
        ret[0] = self.rot_rad
        ret[self.phase_indices] = self.phase_offsets[self.free_antennas]
        return ret

    def bounds(self, pointing_error, test_gains):
        bounds = [0] * self.nend
        bounds[0] = self.pointing_bounds(pointing_center)
        for i in range(1,self.nant):
            tg = test_gains[i]
            if tg < 0.01:
                bounds[i] = (0, 1e-3)
            else:
                bounds[i] = (-np.pi*2, np.pi*2) # Bounds for phases

        return bounds

#def split_param(x):
    #global GAIN_INDICES, PHASE_INDICES
    #rot_rad = x[0]
    #if REIM:
        #re = np.concatenate(([1], x[GAIN_INDICES]))
        #im = np.concatenate(([0], x[PHASE_INDICES]))
        #gains = np.sqrt(re * re + im * im)
        #phase_offsets = np.arctan2(im, re)
    #else:
        #gains = np.concatenate(([1], x[GAIN_INDICES]))
        #phase_offsets = np.concatenate(([0], x[PHASE_INDICES]))

    #return rot_rad, gains, phase_offsets



#def join_param(rot_rad, gains, phase_offsets):
    #ret = np.zeros(NEND)
    #ret[0] = rot_rad
    #if REIM:
        #z = gains[FREE_ANTENNAS] * np.exp(phase_offsets[FREE_ANTENNAS] * 1j)
        #ret[GAIN_INDICES] = z.real
        #ret[PHASE_INDICES] = z.imag
    #else:
        #ret[GAIN_INDICES] = gains[FREE_ANTENNAS]
        #ret[PHASE_INDICES] = phase_offsets[FREE_ANTENNAS]
    #return ret

##def setup_indices():
    ##global REIM, NEND, FREE_ANTENNASGAIN_INDICES, PHASE_INDICES, NANT

    ##if REIM:
        ##NEND=int(2*NANT-1)
        ##FREE_ANTENNAS=slice(1,NANT)
        ##GAIN_INDICES=slice(1,NANT)
        ##PHASE_INDICES=slice(NANT, NEND)
    ##else:
        ##NEND=NANT
        ##FREE_ANTENNAS=slice(1,NANT)
        ##GAIN_INDICES=None
        ##PHASE_INDICES=slice(1, NEND)


#def param_to_json(x):
    #rot_rad, gains, phase_offsets = split_param(x)
    #ret = {
        #"gain": np.round(gains, 4).tolist(),
        #"rot_degrees": np.degrees(rot_rad),
        #"phase_offset": np.round(phase_offsets, 4).tolist(),
    #}
    #return ret


#def output_param(x, fp=None):
    #ret = param_to_json(x)
    #if fp is None:
        #print(json.dumps(ret, indent=4, separators=(",", ": ")))
    #else:
        #json.dump(ret, fp, indent=4, separators=(",", ": "))


def calc_score_aux(opt_parameters, measurements, window_deg, original_positions):
    global triplets, ij_index, jk_index, ik_index, masks, ift_scaled, myParam

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

            # Zone outside of mask
            print(f"Max mask {np.max(mask)}, {np.min(mask)}")
            negative_mask = (-mask + 1)
            negative_mask[negative_mask < 0] = 0

            inv_masks[i] = negative_mask
            masks[i] = mask

        mask = masks[i]

        masked_img = masks[i]*ift_scaled
        #outmask_img = inv_masks[i]*ift_scaled

        in_zone = np.sum(np.sqrt(masked_img))  / np.sum(masks[i])
        out_zone = 0 # np.std(outmask_img)

        zone_score = (in_zone)**2
        ret_zone += -zone_score

    ret_std = ret_std / len(measurements)
    ret_zone = ret_zone / len(measurements)

    if N_IT % 100 == 0:
        print(f"S/N {ret_std:04.2f}, ZONE: {ret_zone:04.2f}, in: {in_zone:04.2f} out: {out_zone:04.2f}", end='\r')

    return (
        (ret_zone + ret_std),
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

    if N_IT % 1000 == 0:
        #print(f"Iteration {N_IT}, score={ret:04.2f}")
        f_vs_iteration.append(ret)

    N_IT += 1
    return ret

from scipy import optimize
import json


class MyTakeStep(object):
    def __init__(self, stepsize, pointing_rad):
        self.stepsize = stepsize
        self.pointing_rad = pointing_rad

    def __call__(self, x):
        return myParam.take_step(x, self.stepsize)


def bh_callback(x, f, accepted):
    global output_directory, bh_basin_progress, N_IT, ift_scaled, masks, method
    print(f"BH f={f} accepted {accepted}")
    myParam.from_vector(x)
    myParam.output()
    if accepted:
        bh_basin_progress.append([N_IT, f])
        with open(f"{output_directory}/bh_basin_progress.json", "w") as fp:
            json.dump(bh_basin_progress, fp, indent=4, separators=(",", ": "))

        with open(f"{output_directory}/BH_basin_{f:5.3f}_{N_IT}.json", "w") as fp:
            myParam.output(fp)

        mask = masks[0]
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
        "--use-phases",
        action="store_true",
        help="Use Real/Imaginary components rather than gains and phases.",
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
        "--max-delay", type=float, default=1, help="Maximum delay in wavelengths")

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

    pointing_error = np.radians(3)
    pointing_center = np.radians(ARGS.pointing)

    if ARGS.use_phases:
        myParam = ParamPhase(NANT, pointing_error, gains)
        myParam.rot_rad = 0
    else:
        myParam = ParamReIm(NANT, pointing_error)
        myParam.rot_rad = 0
        myParam.gains = gains
        myParam.phase_offsets = phase_offsets

    #init_parameters = join_param(0.0, gains, phase_offsets)
    #output_param(init_parameters)
    myParam.output()

    masks = []
    inv_masks = []
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
        inv_masks.append(None)

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
    
    best_acq = np.zeros(NANT)
    n = 0
    best_score = -999
    for acquisition_data in full_acquisition_data:
        print(acquisition_data.keys())
        for d in acquisition_data:
            acq = acquisition_data[d]
            ph = np.array(acq['phases'])
            st = np.array(acq['strengths'])

            mean_str = np.median(st)
            
            if mean_str > 7.0:
                best_acq += st
                n = n + 1

            print(f"Source: {int(d):02d}, stability: {np.std(ph):06.5f}, {np.mean(st):05.2f} {acq['sv']}")

    if n == 0:
        raise RuntimeError("No satellites visible")
    
    best_acq = best_acq / n
    
    # Now remove satellites from the catalog that we can't see.
    # https://github.com/JasonNg91/GNSS-SDR-Python/tree/master/gnsstools
        
    N_IT = 0
    window_deg = 8.0

    s = calc_score(
        myParam.to_vector(),
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


    print(f"Calculating which antennas to ignore {best_acq}")
    test_gains = best_acq / best_acq[0]
    print(f"Estimated gains: {test_gains}")

    bounds = myParam.bounds(pointing_center, test_gains)

    init_parameters = myParam.to_vector()

    zero_list = ARGS.ignore

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

    myParam.from_vector(ret.x)
    rot_rad = myParam.rot_rad
    output_json = myParam.to_json()
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

    print(f"Optimal solution: {output_json}")
    f_history_json = {}
    f_history_json["start"] = s
    f_history_json["history"] = f_vs_iteration

    with open("{}/{}_history.json".format(output_directory, method), "w") as fp:
        json.dump(f_history_json, fp, indent=4, separators=(",", ": "))
