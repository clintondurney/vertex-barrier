# vertex-barrier

## About
This repository contains code that was used in a 2D Vertex Model simulation of tissue mechanics. This work was included in the manuscript titled, "Zasp52 strengthens whole embryo tissue integrity through supracellular actomyosin networks", Ashour et. al. To be submitted in 2022. 

The following code has been tested on macOS Monterey 12.5 and Ubunutu 20.04.  Earlier versions of each OS should work, but have not been tested.  Windows users should be able to get the code working by mirroring the steps below for Windows.

## Author:
* [Clinton H. Durney](https://clintondurney.github.io/) (clinton.durney@jic.ac.uk)

## Dependencies 
For running the numerical simulations (no plotting), the .yml should install everything needed (see below).

For plotting snapshots of the tissue configuration, [POV-Ray](http://www.povray.org/) is needed. 

  > Ubuntu: https://www.povray.org/download/linux.php

  > macOS: https://wiki.povray.org/content/HowTo:Install_POV

For creating a seamless video from the snapshots, ffmpeg is needed. 

  > Ubuntu/macOS: https://ffmpeg.org/download.html

## To use:
__Installation__
Use [Conda](https://docs.conda.io/en/latest/) to create the environment from the .yml file:
```
conda env create -f vertex.yaml python=3.7
```

__Running__
1. Activate the conda environment
```
conda activate vertex
```
2. Run main.py 
```
python main.py
```

This will populate the directory with .pickle files that contain the attributes of the network state and .npy files that save the nodes of the cells circumference oriented counter clockwise. These files can be analyzed using a Jupyter Notebook or visualized using the repository referenced above.

3. (Optional) Create pictures
```
python plot.py
```

4. (Optional) Create video
```
python make-vid.py
```

## Citation:
If you find the code useful, please consider citing either or both of the following (to be updated upon publication of this work):
1. [Durney et. al. (2018)](https://www.sciencedirect.com/science/article/pii/S0006349518311615)
2. [Durney et. al. (2021)](https://iopscience.iop.org/article/10.1088/1478-3975/abfa69/meta)








