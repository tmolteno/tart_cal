import argparse
import json
import glob

import numpy as np
import matplotlib.pyplot as plt

from tart.operation import settings


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Visualize calibration data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--data",
        required=False,
        default="cal_data",
        help="Calibration Input Data Directory.",
    )

    parser.add_argument(
        "--rot",
        required=False,
        type=float,
        default=30,
        help="Rotated 30 degrees.",
    )

    ARGS = parser.parse_args()

    # Load calibration data from the data directory

    data_dir = ARGS.data
    json_files = [f for f in glob.glob(f"{data_dir}/cal*.json")]

    with open(json_files[0], "r") as json_file:
        calib_info = json.load(json_file)

    info = calib_info["info"]
    ant_pos = calib_info["ant_pos"]

    ant_pos = np.array(ant_pos)
    x = ant_pos.T[0]
    y = ant_pos.T[1]

    fig, ax = plt.subplots(figsize=(5,5))

    ax.plot(x,y,'.', label='original positions')
    
    rot_deg = ARGS.rot
    
    print(ant_pos)
    new_positions = settings.rotate_location(
        rot_deg, np.array(ant_pos).T
    )

    ax.grid(True)

    new_pos = np.array(new_positions).T
    print(new_pos.T)
    
    x = new_pos.T[0]
    y = new_pos.T[1]
    ax.plot(x,y,'o', label=f"Rotated {rot_deg} deg")
    plt.legend()
    plt.show()
