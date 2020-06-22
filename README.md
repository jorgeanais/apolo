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
 

Module Structure
----------------

The module is divided in the following parts:

 1. `catalog_proc`: this module performs the pre-processing of the catalogs,  including matching,
  cleaning and feature extraction.
 
 2. `clustering`: this module has the tools to perform the clustering and setup data.
 
 3. `data`: this module include data relevant to tiles in the region of interest and stellar cluster to do some test.
 
 4. `test_tools`: useful functions used to perform tests


Folder Structure
----------------

Catalogs are created/saved in a folder structure for their posterior usage. The folder structure is
defined in `data/dirconfig.py` file. 

```
DATA/
├── cross
│   ├── x_vvv-2mass
│   ├── x_vvv-2mass-combis
│   ├── x_vvv-2mass-combis-gaia
│   ├── x_vvv-2mass-combis-gaia_cont
│   ├── x_vvv-combis
│   ├── x_vvv-combis-gaia
│   ├── x_vvv-combis-gaia_cont
│   ├── x_vvv-gaia
│   └── x_vvv-gaia_cont
├── proc
│   ├── p_2mass
│   ├── p_combis
│   ├── p_gaia
│   └── p_vvv
├── raw
│   ├── 2mass
│   ├── combis
│   ├── gaia
│   └── vvv
└── test
    └── known_clusters

```

Installation
------------

1. To install Apolo, move to the folder where the setup.py is and then just execute

```
pip install -e .
```

Usage
-----

In the `scripts` folder, there are some codes that demonstrate how use apolo:

1. `script_proc.py`: It pre-process and match raw catalogs in order to be used by the clustering algorithm.
2. `basic_clustering.py`: A simple example of clustering. All the parameters are set by hand.
3. `test_known_clusters`: This script automatically perform our current methodology to a set of clusters in parallel.



