# Exoplanet Radial Velocity GIFs

A Python script to fetch exoplanet data from NASA’s Exoplanet Archive, model Keplerian orbits and radial velocity curves, and produce animated GIFs visualizing orbital motion alongside the star’s wobbles.
---

## Usage

Run the script from the command line: python exoplanet_analysis.py


You will be prompted to enter a star name (e.g. `HD189733` or `HD 189733`). The script then:

1. Downloads the exoplanet table.
2. Finds the matching star and its planets.
3. Computes orbital parameters and radial velocities.
4. Outputs `orbit_rv.gif` in the working directory.

---

## Examples

The `examples/` folder contains a sample animation: `orbit_rv.gif`.

