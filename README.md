# eddy-characteristics

Code and figures used in my 2020 CICOES internship with Dr. Mariona Claret (CICOES) and Dr. Pascale Lelong (NWRA) on the differences in dynamical features and interactions between near-inertial waves, cyclonic and anticyclonic eddies in the northwestern Mediterranean.

**EddyTracks** holds functions used in data processing and helper functions for statistics

**eddy_map** animates eddy centers and shapes at a depth of 500m over a map of the Gulf of Lion contoured by Rossby number for the month of July

**eddy_stats** creates histograms for eddy lifetime (days), propagation (km), Rossby number, and max radius (km)

**spectra_area** creates rotary spectra of the open ocean in the Gulf of Lion at the surface, 500m, and 1500m

**spectra_eddies** creates rotary spectra following eddy centers in the Gulf of Lion at the surface, 500m, and 1500m

## Dependencies
* numpy
* matplotlib
* h5py
* scipy
