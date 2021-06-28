Apolo
=====

> Status: In development 

This code is being developed as part of my master thesis in Astronomy. The idea is
to develop an algorithm to find stellar-clusters, focusing in the region of
the far end of the Milky Way.

Our approach consist in apply unsupervised machine learning techniques in a large
data-set of near-infrared photometry from VVV project (Saito et al. 2012). 
This data is specially suitable for this kind
of search due to range of wavelengths are less affected by interstellar extinction.
The PSF photometry catalogs came from the work by Alonso-García et al. 2018.

We complement this photometric 
catalog with proper-motions from VIRAC (Smith et. al 2018) and 
we also decontaminate the original sample by removing objects that are closer than expected
for our region of interest, using astrometric-data from Gaia DR2.

The original idea for this code is not just that it serves as a tool for my work, but also
to have a record of the things that I have been trying, so it is not intended for public use.
 

Module Structure
----------------

The module is divided into the following parts:

 1. `catalog_proc`: this module contains tools to perform the pre-processing stage. The objective 
 of this stage is to get, transform and combine raw catalogs into a set of new catalogs in a common
 format that can be used in the clustering stage. 
  
 2. `clustering`: this module has the tools to perform clustering itself.
 
 3. `data`: this module include data relevant to tiles in the region of interest and stellar clusters.
 
 4. `test_tools`: functions used to test different clustering approaches, before launch the search stage.


Folder Structure
----------------

Catalogs are created/saved in a folder structure for their posterior usage. The folder structure is
defined in `apolo/data/dirconfig.py` file. 

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
    └── tiling

```

Installation
------------

1. First, download the code
 
2. Then, edit file `apolo/data/dirconfig.py` in line 8 to indicate where your base directory is.
```
...
# Base path definition
base_data_path = '/path/to/DATA'
...
```
2. To install Apolo, move to the folder where the setup.py is and then just execute

```
pip install -e .
```


Usage
-----

In the `scripts` folder, there are some codes that demonstrate how use apolo:

1. `script_proc.py`: It pre-process and match raw catalogs in order to be used by the clustering algorithm.
2. `basic_clustering.py`: A simple example of clustering. All the parameters are set by hand.
3. `test_known_clusters.pu`: This script automatically perform our current methodology to a set of clusters in parallel.
4. `join_tiles.py`: It spatially merge all the catalogs in one (in order to avoid spatial cuts in our region of interest).
5. `tiling.py`: It spatially separate the region of interest in tiles using a density-aware approach (to make it computationally approachable).
6. `perform_clustering_tiles.py`: It perform the search in the region of interest in a parallel fashion.
7. `candidates_selection.py`: Extract some stats from the resulting gruops and select candidates based on filtering some features.


Results
-------

We were able to characterize for first time, using near-infrared photometry and proper motions,
the massive stellar population of three clusters currently reported at the far end of the Galactic
bar: [MCM2005b]81, VVV CL074, and VVV CL086. The results from this study suggest that cluster VVV CL086
is closer than previously reported and therefore not in the region of the bar's far end. From the search
using Apolo, we identify three new young massive cluster candidates and recover one known cluster.
The new star cluster candidates are _Apolo 2_ at Galactic coordinates l=336.546°, b=-0.018°;
_Apolo 3_ at l=336.425°, b=0.187°; and _Apolo 4_ at l=336.580°, b=0.021°, with near-infrared
color J - K_s ~ 3.5 mag and mean proper motions, <μ<sub>α*</sub>>= -4.5±0.8 mas/yr and <μ<sub>ẟ</sub>> = -5.2±0.7 mas/yr,
consistent with the expected proper motions for objects at the far end of the bar according to the Besançon
model of stellar population synthesis in our Galaxy, and with the proper motions of the known cluster
at this place.

Download my thesis [here!](https://drive.google.com/file/d/1OcrNnFa7S-o3oVHwSqg-SB1IkH94lpzt/view)

![image](refimage.png)
Cluster candidate Apolo 2.

