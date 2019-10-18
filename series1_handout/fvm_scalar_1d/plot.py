import json
import glob

import numpy as np
import matplotlib.pyplot as plt


def load_json(filename):
    with open(filename, "r") as f:
        return json.load(f)

def load_snapshot(filename):
    j = load_json(filename)
    return j["time"], np.array(j["cell_centers"]), np.array(j["data"])

def plot_snapshot(t, x, u):
    plt.plot(x, u)
    plt.title(f"t = {t:0.3f}")
    plt.show()

if __name__ == "__main__":

    files = sorted(glob.glob("testing*.json"))
    for fn in files:
        t, x, u = load_snapshot(fn)
        plot_snapshot(t, x, u)
