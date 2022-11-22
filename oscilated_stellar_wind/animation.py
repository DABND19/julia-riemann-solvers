from itertools import islice
import json
import math

import jmespath
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


def get_ylimit(records: np.matrix, *, gap: float = 0.2) -> tuple[float, float]:
    y_min = records.min()
    y_max = records.max()
    interval_len = y_max - y_min
    return y_min, y_max


def get_xlimit(values: np.ndarray) -> tuple[float, float]:
    x_min = values.min()
    x_max = values.max()
    return x_min, x_max


if __name__ == "__main__":
    with open("data/test_3.jsonl") as f:
        records = list(map(json.loads, islice(f, 0, 5000, 5)))
        x = np.array(jmespath.search("[0].x", records)[7500:12001])
        density = np.array(jmespath.search("[*].solution[*].density", records))

        fig, ax = plt.subplots()

        ax.set_title("Oscilations of density when A=33%.")
        [line] = ax.plot(x, density[0][7500:12001])
        # ax.set_xlim(*get_xlimit(x))
        ax.set_ylabel("Density")
        ax.set_xlabel("Radius")
        ax.set_ylim(
            ymin=density.min(),
            ymax=density.max(),
        )
        # ax.set_ylim(
        #     ymin=math.pow(10, math.log10(density.min()) - 1),
        #     ymax=math.pow(10, math.log10(density.max()) + 1),
        # )
        # ax.set_yscale("log")

        def update(frame: int):
            line.set_ydata(density[frame][7500:12001])
            return (line,)

        anim = animation.FuncAnimation(
            fig, update, interval=50, frames=range(len(records))
        )
        anim.save("oscilated_33.gif")
