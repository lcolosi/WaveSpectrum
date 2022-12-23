# Source Code for: 

Luke Colosi, Nicholas Pizzo, Laurent Grare, and Luc Lenain. Observations of Surface Gravity Wave Spectra from Moving Platforms. Journal of Atmospheric and Oceanic Technology, in prep. 

# Abstract 

Surface waves play an important role in the ocean-atmosphere coupled climate system by mediating exchanges of momentum, heat, and gas between the atmosphere and the ocean. Pseudo Lagrangian autonomous platforms (e.g., Boeing Liquid Robotics Wave Gliders) have been used to investigate the underlying physical dynamics involved in these processes to parameterize better the air-sea exchange occurring at the scale of the surface waves. This requires accurate measurements of directional surface waves down to short scales (O(1) meter), as these shorter waves support most of the stress between the atmosphere and ocean.  A challenge to overcome for pseudo Lagrangian autonomous vehicles is that the platform velocity causes the observed frequency of the waves to be Doppler shifted.  This leads to a modulation of the wave spectrum, particularly at high frequencies, that depends on the platform’s speed, the wave frequency, and the relative angle between the direction of wave propagation and the platform heading. In this work, we propose a novel approach correcting this effect, building upon the work of Longuet-Higgins (1986) and more recently by Collins et al. (2017) and Amador et al. (2022). The method is validated using a unique dataset collected from a fleet of Wave Gliders off the coast of Southern California in September 2019 operating on a tight square (500 m edge length) track over a 3-day deployment. This technique can be used to estimate wave spectra derived from other slow-moving surface vehicles such as Saildrones that use platform motion to characterize the surface wave field.

# Plain Language Summary

Waves at the ocean surface play an essential role in the Earth's climate system by influencing how the ocean and atmosphere interact. These interactions occur between winds, waves, and currents as they exchange properties such as energy, heat, and gases. Uncrewed vehicles that glide across the ocean surface have been used to investigate ocean-atmosphere interactions at the space and time scales of surface waves with the goal of better understanding ocean physics. This information can be applied in climate models to improve forecasting capabilities. Achieving this requires accurate measurements of short wavelength waves because these waves are crucial for transferring energy from winds to waves. A challenge for moving vehicles is that the platform's speed and direction relative to the incoming waves influence the wave period observed aboard the vehicle. This causes the energy associated with each distinct wave period present in the wave field to be shifted to a different period. Here, we introduce a new approach to correct this effect so that wave properties that are dependent on wave period are not affected by the platform's motion. This novel approach may be applied to any uncrewed surface vehicle. The method is tested using a unique data set collected from a fleet of Wave Gliders, an uncrewed surface vehicle.

# Authors 
* [Luke Colosi](https://lcolosi.github.io/)<<lcolosi@ucsd.edu>>
* [Nicholas Pizzo](https://scripps.ucsd.edu/profiles/npizzo) <<npizzo@ucsd.edu>>
* [Laurent Grare](https://airsea.ucsd.edu/people/) <<lgrare@ucsd.edu>>
* [Luc Lenain](https://scripps.ucsd.edu/profiles/llenain) <<llenain@ucsd.edu>>

# Data

## Experiments

In this study, observations collected from two experiments are investigated. The first was conducted as part of a project funded through the ONR Task Force Ocean (TFO) initiative. Observations from two Wave Gliders, Planck and Stokes, were collected off the coast of Del Mar, California (see Figure 1) during a 3-day deployment from September 9th to 11th, 2020 (referred hereon as DELMAR2020). The two Wave Gliders moved along two tight square tracks: a large square with a 1000 meter edge length and a small square with a 500 meter edge length. Here, we focus on data from the small square where Doppler effects on the frequency spectrum are more apparent for preliminary analysis. 

The second was conducted as part of the NASA Submesocale Ocean Dynamics Experiment (S-MODE) program, a project that aims to determine whether submesoscale ocean dynamics make important contributions to the vertical exchange of climate and biological variables in the upper ocean using a combination of aircraft-based remote sensing, research vessel, autonomous oceanographic platforms measurements along with numerical modeling (Farrar et al. 2020). The pilot experiment considered here was conducted off the coast of San Francisco, California (see Figure 2) from October 29th to November 4th, 2021 (referred hereon as SMODE2021). One wave glider, WHOI43, was deployed by the Upper Ocean Physics Laboratory at Woods Hole Oceanographic Institution led by Tom Farrar. WHOI43 drove in a combination of long transects and in rectangular trajectories.  

## Instrumentation

Instrumented Boeing Liquid Robotics SV3 wave gliders are used to collect ocean and atmospheric observations at the air-sea interface (Grare et al. 2021). Atmospheric measurements include wind speed and direction from the Vaisala WXT (model 530) mounted $\sim 1$ meter above the ocean surface.The wind speed 10 meters above the ocean surface is computed assuming a logarithmic wind speed profile and a neutrally stable atmosphere (Charnock 1955). Ocean wave measurements are derived from the motion of the platform (Lenain and Melville 2014, Thomson et al. 2018, Grare et al, 2021) which is collected by Novatel's Synchronous Position, Altitude, and Navigational (SPAN) technology providing a combined solution from a Novatel OEM7720 global positioning system (GPS) receiver and a Epson EG320N inertial motion unit (IMU) for DELMAR2020 and a VectorNav VN-300 GPS/IMU for SMODE2021.

# Funding
This work was supported by the NASA  (award [put ID# here]) and the Office of Naval Research (award [put ID# here]).

# How to use this repository

All figures in Colosi et al. (2023) can be reproduced using the MatLab scripts from this repository and processed [data](insert collection url) published to the UCSD library digital collections. To do so, follow these steps:

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

**Note on MatLab software version:** Code was developed using MatLab R2022a. Version released before and after this may run into some errors. 

# How to cite this code

If you wish to use the code from this repository, you may cite it as: 

[put citation from zenodo]. 
