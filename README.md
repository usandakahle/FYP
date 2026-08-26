# FYP
Final year project development files for **SHARC- Regime aware Directional waves measurement pipelin**. This repo contains the MATLAB simulation pipeline used to develop and validate a directional wave-measurement capable SHARC BUOY, building on from Noyce's thesis on the HSARC buoy hardware and firmware.

# First Iteration is on MATLAB
## Modeling the Open Ocean Wave Spectrum
Using parameters taken from ERA5 reanalysis for the period **11 July 2022- 31 July 2022** (the SCALE Winter Cruise 2022 deployment window), an initial synthetic sea-surface spectrum was developed to match a target theoretical spectrum. The developed **eta(t)** is modelled using a random-phase/Rayleigh-amplitude representation of a JONSWAP spectrum parameterized by ERA5's significant wave height and peak wave period.

Validation is performed using **Welch's windowing Method** to confirm that the generated time series recovered power spectral density matches the target spectrum **E(f)**.

## SHARC BUOY measurement
From **Noyce's thesis** on the SHARC buoy, parameters relating to the onboard sensors (ICM20649 IMU) and the existing DSP/firmware pipeline were extracted, so tha a MATLAB script could simulate what the SHARC buoy would actually measure and report when presented with the synthetic spectrum from **ocean_wave_spectrum.m**, including sensor noise, decimation, filtering and the buoy's onboard spectral estimation method (Bartlett-averaged PSD).

This was then extended to simulate a **directional** measuremnt capability (heave-pitch-roll method), using the buoy's existing gyroscope in combitaion with the accelerometer, since the buoy has no continuous horizontal position source suitable for direction recovery.

##Repository structure
**ocean_wave_spectrum.m** %Synthetic target speactrum & surface generation
**sharc_meas.m** % SHARC sensing + DSP pipeline simulation
**directional_validation.m** %Checks on the directional generation
**Era5.nc** % ERA5 wvae variables


