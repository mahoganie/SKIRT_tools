from functools import lru_cache
from importlib.resources import as_file, files
import math
import random

import numpy as np
import pandas as pd

PACKAGE_DIR = files(__package__)


@lru_cache(maxsize=1)
def _load_predefined_bands():
    band_file = PACKAGE_DIR / "PredefinedBandWavelengths.txt"
    with as_file(band_file) as path:
        return pd.read_csv(
            path,
            sep=r"\s+",
            engine="python",
            comment="#",
            names=["band", "wavelength_micron"],
        )

def get_predefined_band_name(wavelength):
    df = _load_predefined_bands()
    idx = np.argmin(np.abs(df["wavelength_micron"].to_numpy() - wavelength))
    return str(df.loc[idx, "band"])


def calculate_abs_mag(flux_Jy, distance_mpc=50):
    distance_pc = distance_mpc * 1e6  # distance in pc
    zero_point = 3631  # in Jy, zero-point for AB magnitudes
    AB_mag = -2.5 * np.log10(flux_Jy / zero_point)
    abs_mag = AB_mag - (5 * np.log10(distance_pc)) + 5
    return abs_mag


def random_viewing_angles():
    # Azimuthal angle: uniform in [0, 360)
    azimuthal = random.uniform(0, 360)
    # Inclination: sample cos(theta) uniformly
    u = random.uniform(-1, 1)
    inclination = math.degrees(math.acos(u))
    # Roll: uniform in [-360, 360]
    roll = random.uniform(-360, 360)
    return inclination, azimuthal, roll
    


# From mockobservation_tools (Courtney Klein)
# def calculate_M_obs(g_r, Lsun, pixelscale=60/400):
#     # L = 10**(0.4*(4.65-r_band_mag))
#     M_obs = obs_mass2(Lsun, g_r)
#     return M_obs