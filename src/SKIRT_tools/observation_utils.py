import math
import random
import pandas as pd
import numpy as np
from pathlib import Path

SKIRT_tools_dir = Path(__file__).resolve().parent

def get_predefined_band_name(wavelength):
    df = pd.read_csv(
        str(SKIRT_tools_dir / "PredefinedBandWavelengths.txt"),
        sep='\s+',engine='python',comment='#',names=['band','wavelength_micron'])
    idx = np.argmin(np.abs(df['wavelength_micron'].values - wavelength))
    return df.loc[idx,'band']

def calculate_abs_mag(flux_Jy, distance_mpc=50):
    distance_pc = distance_mpc * 1e6 # distance in pc
    zero_point = 3631 # in Jy, zero-point for AB magnitudes
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
