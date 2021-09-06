# MBRWSpectralDimensionPublic
As part of my dissertation, I am investigating the usefulness of spectral dimension as a measure of network complexity.
We estimate spectral dimension by fitting a power law to the mean number of nodes encountered in a segment of the walk (segment mass) as a function of the length of the segment.
The spectral dimension is 2 times the exponent of the power law.
In addition to conventional spectral dimension, our lab has developed a spectral dimension based on a memory-biased random walk (MBRW) as opposed to a simple random walk and a measure of multi-spectrality analogous to existing measures of multi-fractality.

This repository contains code from the project that I am making public with the consent of my faculty advisor, Dr. Uri Hershberg.
If you use this code, the MBRW spectral dimension, or multi-spectrality defined here, please cite
Adam Craig, Mesut Yücel, Lev Muchnik, Uri Hershberg (2021)
Comparing network complexity using random walk and memory biased random walk (MBRW) spectral dimension and multi-spectrality
presented at the 2021 NHLBI Systems Biology Symposium virtual-only conference

To run the demo, install the g++ compiler and MATLAB 2018a or later, download the Boost C++ library, place the Boost files in a folder in the parent folder of this folder, and run the MATLAB script mbrw_spectral_dimension_demo.m.
html/mbrw_spectral_dimension_demo.pdf shows the expected output of the script.

You can also compile and run the C++ executable that performs MBRW and saves log generalized means of walk segment masses.
Follow the format used in the commands passed to system() in the MATLAB script.