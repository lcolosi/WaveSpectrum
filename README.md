# Source Code for: 

Luke Colosi, Nicholas Pizzo, Laurent Grare, and Luc Lenain. Observations of Surface Gravity Wave Spectra from Moving Platforms. Journal of Atmospheric and Oceanic Technology, in prep. 

# Abstract 

Surface waves play an important role in the ocean-atmosphere coupled climate system by mediating exchanges of momentum, heat, and gas between the atmosphere and the ocean. Pseudo Lagrangian autonomous platforms (e.g. Boeing Liquid Robotics Wave Gliders) have been used to investigate the underlying physical dynamics involved in these processes to better parameterize the air-sea exchange occurring at the scale of the surface waves. This requires accurate measurements of directional surface waves down to short scales (O(1)m), as these shorter waves support most of the stress between the atmosphere and ocean.  A challenge to overcome for pseudo Lagrangian autonomous vehicles is that the platform velocity causes the observed frequency of the waves to be doppler shifted.  This leads to a modulation of the wave spectrum, in particular at high frequencies, that depends on the platform’s speed, the wave frequency, and the relative angle between the direction of wave propagation and the platform heading. In this work, we propose a novel approach to correcting this effect, building upon the work of Longuet-Higgins (1986) and more recently by Collins et al. (2017). The method is validated using a unique dataset collected from a fleet of Wave Gliders off the coast of Southern California in September 2019 operating on a tight square (500 m edge length) track over a 3-day deployment. This technique can be used to estimate wave spectra derived from other slow-moving surface vehicles such as Saildrones that use platform motion to characterize the surface wave field. 

# Plain Language Summary

# Authors 
* [Luke Colosi](https://lcolosi.github.io/)<<lcolosi@ucsd.edu>>
* [Nicholas Pizzo](https://scripps.ucsd.edu/profiles/npizzo) <<npizzo@ucsd.edu>>
* [Laurent Grare](https://airsea.ucsd.edu/people/) <<lgrare@ucsd.edu>>
* [Luc Lenain](https://scripps.ucsd.edu/profiles/llenain) <<lgrare@ucsd.edu>>

# Data

## Experiments

Data from two experiments are investigated. Data were collected off the coast of Del Mar, California (see Figure 1) during a 3-day deployment from September 9th to 11th, 2020 (referred hereon as DELMAR2020). Two Wave Gliders moved along two tight square tracks: large square with 1000 m edge length and small square with 500 m edge length. Data were collected during the submesoscale ocean dynamics experiment pilot program off the coast of San Francisco from October 29th to November 4th, 2021 (referred hereon as SMODE2021). A wave glider deployed by the Woods Hole Oceanographic team moved along in rectangular trajectories. Processed [data](insert collection url) from DELMAR2020 and SMODE2021 may be found in UCSD library digital collections.  

## Instrumentation

Instrumented Boeing Liquid Robotics SV3 wave gliders are used to collect ocean and atmospheric observations at the air-sea interface (Grare et al. 2021). Atmospheric measurements include wind speed and direction from the Gill 3D sonic Anemometer. Ocean wave measurements are derived from the motion of the platform (Lenain and Melville 2014, Thomson et al. 2018, Grare et al, 2021) recorded by Novatel's geo-positioning system (GPS) and inertial motion unit (IMU) for DELMAR2020 and VectorNav for SMODE2021.

# Funding
This work was supported by the NASA  (award [put ID# here]) and the Office of Naval Research (award [put ID# here]).

# How to use this repository

All figures in Colosi et al. (2022) can be reproduced using the MatLab scripts from this repository and processed [data](insert collection url) published to the UCSD library digital collections. To do so, follow these steps:

1. Make a local copy of this repository by either cloning or downloading it.

2. Download the processed [data](insert collection url), untar the files, and move all directories to `data` in the project root. After doing so, your directory tree should look like this:

```
WaveSpectrum/
├── data
│   ├── DelMAR2020
│   ├── SMODE2021
│   └── BATHY
├── figs
├── src
└── tools
```

3. Make sure that you have downloaded the R2022a version of MATLAB. Other versions may work, but small difference in keyword argument and outputs of functions may arise.   

4. If you follow the steps above you should be able to reproduce all figures, by running `figXX.m` from the `src` directory without having to adjust any paths.

# How to cite this code

If you wish to use the code from this repository, you may cite it as: 

[put citation from zenodo]. 
