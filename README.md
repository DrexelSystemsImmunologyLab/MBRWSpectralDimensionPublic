# MBRWSpectralDimensionPublic
This repository is for the publicly released versions of source code used in our research into memory-biased random walk (MBRW) generalized spectral dimension as a tool for studying the community structure and heterogeneity of biological networks.
Our definition of spectral dimension is similar to the formulation in
Rammal, R., & Toulouse, G. (1983). Random walks on fractal structures and percolation clusters. Journal de Physique Lettres, 44(1), 13-22.
We estimate spectral dimension by fitting a power law to the mean number of nodes encountered in a segment of the walk (segment mass) as a function of the length of the segment.
The spectral dimension is 2 times the exponent of the power law.
However, instead of using a simple random walk, we have developed a version based on MBRW, which increases its sensitivity to meso-scale community structure and robustness against small-scale randomness.
Furthermore, we replace the expected value (simple arithmetic mean) of the segment mass for a given segment length with the generalized mean (also known as the power mean) to create a generalized spectral dimension.
We can then tune the power of the generalized mean to adjust the measure to give more weight to the most or least densely clustered parts of the network.
We can also use the range of the generalized spectral dimensions over a given domain of powers as a measure of multi-spectrality analogous to existing measures of multi-fractality.

This repository contains code I have written for the project that I am making public with the consent of my faculty advisor, Dr. Uri Hershberg.
If you use this code, please cite
Craig, A. G., Yücel, M., Muchnik, L., & Hershberg, U. (2022). Impact of Finite Size Effect on Applicability of Generalized Fractal and Spectral Dimensions to Biological Networks. Available at SSRN 4097638.

To run the demo, install the g++ compiler and MATLAB 2018a or later, download the Boost C++ library, place the Boost files in a folder in the parent folder of this folder, and run the MATLAB script mbrw_spectral_dimension_demo.m.
html/mbrw_spectral_dimension_demo.pdf shows the expected output of the script.

You can also compile and run the C++ command line utility that performs MBRW and saves log-base-2 generalized means of walk segment masses independently of the demo.
When running the utility, follow the format used in the commands passed to system() in the MATLAB script.

I have also included the source code for the utilities used to generate rescaled and fully or partially randomized networks as used in the paper, though these are not currently part of the demo.

All C++ code here uses the Boost C++ library but does not require full installation of Boost, just the presence of the library files.
To compile, use the following: g++ -std=c++0x -I [path to boost folder] [source file name minus extension].cpp -o [source file name minus extension] -O2

Files:
extract_full_gunsalus_et_al_2005.m MATLAB function that reads in the data for the example network in its original format from the supplementary materials in Gunsalus et al. 2005 and extracts a list of the edges that represent direct protein-protein interactions
kcore.m MATLAB function that takes the k-core of a Graph object
LICENSE: standard MIT License file
make_mbrw_able.m MATLAB function that extracts the largest connected component of the 2-core of a network (MBRW can only traverse networks that are connected and have a minimum degree of at least 2.)
make_similar_scalable.cpp: C++ source file for a command line utility that creates a network with mean degree and community structure similar to those of the input network but with a node count close to an input target value
mbrw_and_save_log_multimeans.cpp: C++ source file for a command line utility that performs an MBRW on an input network and outputs base-2 logs of the generalized mean segment masses for the specified power-of-2 segment lengths and the powers -10, -9, ..., 9 , 10 of the generalized mean plus the -inf-mean (minimum value) and +inf-mean (maximum value).
mbrw_and_save_segment_mass_log_multimeans_2.cpp: similar to mbrw_and_save_log_multimeans.cpp but takes in a list of powers of mean to use in place of the -10, -9, ..., 9, 10 (I have included both versions because performance is slightly better for the version where the powers are hard-coded.)
mbrw_spectral_dimension_demo.m: MATLAB script demonstrating the process of downloading the network data, extracting the desired subgraph, compiling the MBRW C++ source code, running the C++ command line utility on the network, loading the segment mass output file, using the segment masses to compute the generalized MBRW spectral dimension, and computing the multi-spectrality as the range of the spectral dimensions. 
orders.txt: a text file that mbrw_and_save_segment_mass_log_multimeans_2.cpp can take as input to specify the powers to use for the generalized means, should have one real-valued decimal number per line (Numbers with large absolute values may lead to Inf or NaN output values due to overflow or underflow errors.)
orders_pm20.txt: a second example of a list of powers to use as input to mbrw_and_save_segment_mass_log_multimeans_2.cpp
randomize_to_eq.cpp: C++ source file for a command line utility that fully or partially randomizes a network (See the paper for details of the randomization types.)
read_graph.m: MATLAB function that reads data from an edge list text file and a node name list text file and returns a Graph object
README.md: the file you are now reading
rmfileextension.m: a MATLAB function that removes the last '.' and any other characters following it from a string
spectral_dimensions_v2.m: a function that takes in a matrix of log-base-2 generalized mean segment masses (assumes that the i-th row contains masses corresponding to walk segment length 2^(i-1) and that each column corresponds to a different power of mean, though it does not need to know the power) and returns a row vector of generalized spectral dimensions, can use any of several methods to select the range of segment lengths to use

-Adam Craig, 2022-08-23