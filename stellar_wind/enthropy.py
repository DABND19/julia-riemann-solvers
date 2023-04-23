import json
import math
import sys

import jmespath
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


def interpolate(x: np.ndarray, y: np.ndarray, knots_count: int) -> np.ndarray:
    x_ = np.linspace(x.min(), x.max(), knots_count)
    return np.array([
        interp1d(x, row, assume_sorted=True)(x_)
        for row in y
    ])


def main():
    _, input_file, output_file = sys.argv
    with open(input_file) as f:
        records = list(map(json.loads, f))

    x = np.array(jmespath.search("[0].x", records))
    t = np.array(jmespath.search("[*].time", records))
    density = np.array(jmespath.search("[*].solution[*].density", records))
    pressure = np.array(jmespath.search("[*].solution[*].pressure", records))
    gamma = 5 / 3
    enthropy = np.log(pressure / np.power(density, gamma))

    i_min, i_max = 10700, 11750
    fig, ax = plt.subplots(1, 1)
    color_map = ax.imshow(
        interpolate(x[i_min:i_max], enthropy[::-1, i_min:i_max], 1000),
        # norm=LogNorm(
        #     vmin=density[::-1, i_min:i_max].min(), 
        #     vmax=density[::-1, i_min:i_max].max()
        # ),
        interpolation='bilinear',
        extent=[x[i_min], x[i_max], t[0], t[-1]],
        aspect="auto",
    )
    plt.colorbar(color_map, ax=ax)

    fig.savefig(output_file)


if __name__ == "__main__":
    main()
