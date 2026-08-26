# FYP
Final year project development files

# First Iteration is on MATLAB
## Modeling the Open Ocean Wave Spectrum
Using parameters taken from the ERA5 for the period 11 July 2022 to 31 July 2022, an initial spectrum was developed to match the target spectrum. The developed spectrum eta, was modelled using the oceanic waves description in Holthuijsen 2007. 

Validation was performed using Welch's windowing method to validate that the model matched the target. 

## SHARC BUOY measurement
From Noyce's thesis on the SHARC buoy, parameters relating to the sensors  (IMU) and what they read were extracted such that a matlab script could be written to simulate what the SHARC buoy measures when presented by the spectrum in the ocean_wave_spectrum.m file. 
