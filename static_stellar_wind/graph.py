import json
import sys

import jmespath
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    _, *data_paths, picture_path = sys.argv

    fig, ((ax_1, ax_2), (ax_3, ax_4)) = plt.subplots(2, 2, figsize=(10, 8), dpi=200)
    for data_path in data_paths:
        with open(data_path) as f:
            data = json.load(f)

        x = np.array(data["x"])
        pressure = np.array(jmespath.search("solution[*].pressure", data))
        density = np.array(jmespath.search("solution[*].density", data))
        velocity = np.array(jmespath.search("solution[*].velocity", data))
        mach_number = np.array(jmespath.search("solution[*].mach_number", data))

        ax_1.plot(x, pressure, "--")
        ax_1.set_yscale('log')
        ax_1.set_title("Pressure")
        ax_1.set_ylabel("$p$")
        ax_1.set_xlabel("$r$")
        ax_1.grid(True)

        ax_2.plot(x, density, "--")
        ax_2.set_yscale('log')
        ax_2.set_title("Density")
        ax_2.set_ylabel("$\\rho$")
        ax_2.set_xlabel("$r$")
        ax_2.grid(True)

        ax_3.plot(x, velocity, "--")
        ax_3.set_yscale('log')
        ax_3.set_title("Velocity")
        ax_3.set_ylabel("$u$")
        ax_3.set_xlabel("$r$")
        ax_3.grid(True)

        ax_4.plot(x, mach_number, "--")
        ax_4.set_yscale('log')
        ax_4.set_title("Mach number")
        ax_4.set_ylabel("$M$")
        ax_4.set_xlabel("$r$")
        ax_4.grid(True)

    fig.tight_layout(pad=1.3)
    fig.savefig(picture_path)
