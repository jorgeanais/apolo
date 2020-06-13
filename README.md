Apolo
=====

> Status: In development 

This code was made as part of my master thesis in Astronomy. The idea is
to develop an algorithm to find stellar-clusters, focusing in the region of
the far end of the Milky Way.

Our approach consist in apply unsupervised machine learning techniques in a large
data-set of near-infrared photometry from VVV project (Saito et al. 2012). 
This data is specially suitable for this kind
of search due to range of wavelengths are less affected by interstellar extinction.
The PSF photometry catalogs came from the work by Alonso-García et al. 2018. Unfortunately,
this data is not publicly available at the present time. 

We complement this photometric 
catalog with proper-motions from VIRAC (Smith et. al 2018) and 
we also decontaminate the original sample by removing objects that are closer than expected
for our region of interest, using astrometric-data from Gaia DR2.

The original idea for this code is not just that it serves as a tool for my work, but also
to have a record of the things that I have been trying, so it is not intended for public use.
 

How it works
------------

The software is divided in the following parts:

 1. `catalog_proc`: this module performs the pre-processing of the catalogs, 
 including matching, cleaning and feature extraction.
 
 2. `clustering`: this module has the tools to perform the clustering, setup data, and do some benchmarks.

Each module contains scripts with step-by-step description of how I obtain the results.
