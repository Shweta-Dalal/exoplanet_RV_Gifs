# exoplanet_plottingRV.py
# --------------------------------------
# Compute and animate exoplanet orbits & radial velocity
# Requires: numpy, matplotlib, mplcyberpunk, imageio, requests, pandas
# --------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import mplcyberpunk
import imageio
import requests
import pandas as pd
import warnings
from io import StringIO

# --------------------------------------
# 1. Keplerian Radial Velocity Function
# --------------------------------------
def radvelocity(t, P, K, e, omega, T0):
    """
    Compute radial velocity (m/s) for a planet at times t.

    Parameters:
        t (array): Time array [same units as P and T0]
        P (float): Orbital period (days)
        K (float): Velocity semi-amplitude
        e (float): Eccentricity [0,1)
        omega (float): Argument of periastron (radians)
        T0 (float): Time of periastron passage

    Returns:
        vr (ndarray): Radial velocity time series
    """
    t = np.asarray(t, dtype=float)
    # Validate inputs
    if P <= 0:
        warnings.warn("Period P must be > 0; defaulting to P=1 day.")
        P = 1.0
    if K < 0:
        warnings.warn("Semi-amplitude K should be non-negative; taking absolute.")
        K = abs(K)
    if not (0 <= e < 1):
        warnings.warn("Eccentricity e must be in [0,1); setting e=0.")
        e = 0.0

    # Mean anomaly
    M = 2*np.pi*((t - T0) % P)/P
    # Solve Kepler's equation iteratively for eccentric anomaly E
    E = M.copy()
    for _ in range(50):
        delta = (M + e*np.sin(E) - E) / (1 - e*np.cos(E))
        E += delta
        if np.all(np.abs(delta) < 1e-8):
            break
    # True anomaly
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    # Radial velocity
    return K * (np.cos(nu + omega) + e*np.cos(omega))

# ------------------------------
# 2. Semi-major Axis Function
# ------------------------------
def semi_major_axis(P, M_star):
    """
    Compute semi-major axis in AU for circular orbit.
    P: orbital period (days)
    M_star: stellar mass (Msun)
    """
    if P <= 0:
        raise ValueError("Period P must be positive.")
    G = 6.67408e-11  # gravitational constant
    M = M_star * 1.989e30  # convert Msun to kg
    P_sec = P * 24*3600
    a_m = (G * M * P_sec**2 / (4*np.pi**2))**(1/3)
    return a_m / 1.496e11  # convert meters to AU

# ------------------------------
# 3. Ellipse Coordinates
# ------------------------------
def ellipse(a, e, n_points=500):
    """
    Return (x,y) for an ellipse of semi-major axis a and eccentricity e.
    """
    if a <= 0:
        raise ValueError("Semi-major axis a must be positive.")
    if not (0 <= e < 1):
        warnings.warn("Eccentricity e must be in [0,1); setting e=0.")
        e = 0.0
    theta = np.linspace(0, 2*np.pi, n_points)
    r = a*(1-e**2)/(1 + e*np.cos(theta))
    return r*np.cos(theta), r*np.sin(theta)

# --------------------------------------
# 4. Fetch Exoplanet Archive Data
# --------------------------------------

def fetch_planet_data():
    """
    Download CSV from NASA Exoplanet Archive (default_flag=1).
    Returns a pandas DataFrame.
    """
    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
    query = "SELECT * FROM ps WHERE default_flag=1"
    resp = requests.get(url, params={"query": query, "format": "csv"})
    resp.raise_for_status()
    # Avoid DtypeWarning by disabling low_memory
    return pd.read_csv(StringIO(resp.text), low_memory=False)

# --------------------------------------
# 5. Star Lookup Utilities
# --------------------------------------
def find_star_records(name, df):
    """
    Find rows matching `name`, ignoring spaces & case.
    """
    norm_input = name.replace(' ', '').upper()
    # normalize hostnames in DataFrame
    host_norm = df['hostname'].str.replace(' ', '').str.upper()
    mask = host_norm == norm_input
    if not mask.any():
        raise LookupError(f"Star '{name}' not found. Try including spaces or correct casing.")
    return df[mask]

# --------------------------------------
# 6. Plot & Animate
# --------------------------------------
frames = []
def plot_orbits_and_rv(time, rv_list, a_list, e_list, colors, star_name):
    """
    Build frames for orbit animation & radial velocity curve.
    Appends each frame to global `frames` list.
    """
    total_rv = np.sum(rv_list, axis=0)
    for i, t0 in enumerate(time):
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,4),
                                       gridspec_kw={'width_ratios':[1,3]})
        fig.patch.set_facecolor('black')
        ax1.axis('off'); ax2.axis('off')

        # Plot star
        ax1.scatter(0,0, s=150, c=colors[0], marker='*')
        ax1.set_title(star_name, color=colors[-1], pad=20)

        # Orbits & planet positions
        for P, a, e, c, rv in zip(periods, a_list, e_list, colors[1:], rv_list):
            x_e, y_e = ellipse(a, e)
            ax1.plot(x_e, y_e, c=c, lw=0.7)
            # current anomaly
            theta = 2*np.pi*(t0 % P)/P
            r = a*(1-e**2)/(1+e*np.cos(theta))
            ax1.scatter(r*np.cos(theta), r*np.sin(theta), c=colors[-1])
        ax1.set_aspect('equal')

        # Radial velocity curve
        ax2.plot(time, total_rv, lw=0.8)
        ax2.scatter(t0, total_rv[i], c=colors[-1])
        mplcyberpunk.make_lines_glow(ax2)
        ax2.set_xlabel('Time [days]')
        ax2.set_ylabel('RV [m/s]')

        # capture frame
        fig.canvas.draw()
        img = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        img = img.reshape(fig.canvas.get_width_height()[::-1]+(3,))
        frames.append(img)
        plt.close(fig)

# --------------------------------------
# Usage Example & Main
# --------------------------------------
if __name__ == '__main__':
    # Load data
    planets_df = fetch_planet_data()

    # Prompt user
    star_name = input("Enter star name (e.g. HD189733 or HD 189733): ")
    try:
        star_df = find_star_records(star_name, planets_df)
    except LookupError as e:
        print(e)
        exit(1)

    # Extract parameters
    Teff = star_df['st_teff'].iloc[0]
    M_star = star_df['st_mass'].iloc[0]
    periods = star_df['pl_orbper'].values
    eccs    = star_df['pl_orbeccen'].fillna(0).values
    Ks      = star_df['pl_rvamp'].fillna(0).values
    omegas  = np.deg2rad(star_df['pl_orblper'].fillna(0).values)
    T0s     = star_df['pl_tranmid'].fillna(0).values

    # Time grid
    time = np.linspace(0, periods.max(), 300)
    # Compute semi-major axes
    a_list = [semi_major_axis(P, M_star) for P in periods]
    # Compute RV for each planet
    rv_list = [radvelocity(time, P, K, e, w, T0)
               for P,K,e,w,T0 in zip(periods, Ks, eccs, omegas, T0s)]

    # Plot & animate
    plot_orbits_and_rv(time, rv_list, a_list, eccs,
                        colors=['#F5D300','#08F7FE','#FE53BB','#00ff41'],
                        star_name=star_df['hostname'].iloc[0])

    # Save GIF
    gif_file = 'orbit_rv.gif'
    imageio.mimsave(gif_file, frames, fps=10)
    print(f"Animation saved to {gif_file}")
