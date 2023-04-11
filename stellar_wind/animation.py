import json
import math
import sys

import jmespath
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def main():
    _, input_file, output_file = sys.argv
    with open(input_file) as f:
        records = list(map(json.loads, f))

    x = np.array(jmespath.search("[0].x", records))
    t = np.array(jmespath.search("[*].time", records))
    density = np.array(jmespath.search("[*].solution[*].density", records))

    fig, ax = plt.subplots()

    i_min, i_max = 5600, 11000
    [line] = ax.plot(x[i_min:i_max], density[0][i_min:i_max])
    # ax.set_xlim(*get_xlimit(x))
    ax.set_ylabel("Density")
    ax.set_xlabel("Radius")
    ax.set_ylim(
        ymin=density[:, i_min:i_max].min(),
        ymax=density[:, i_min:i_max].max(),
    )
    ax.grid(True)
    # ax.set_ylim(
    #     ymin=math.pow(10, math.log10(density[:, i_min:i_max].min()) - 1),
    #     ymax=math.pow(10, math.log10(density[:, i_min:i_max].max()) + 1),
    # )
    # ax.set_yscale("log")

    def update(frame: int):
        line.set_ydata(density[frame][i_min:i_max])
        line.set_label(f"t={t[frame]:.2f}")
        ax.legend(loc="lower right")
        return (line,)

    anim = animation.FuncAnimation(
        fig, update, interval=50, frames=range(len(records))
    )
    anim.save(output_file)


if __name__ == "__main__":
    main()
