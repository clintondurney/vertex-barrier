# vertex-barrier

## About
This repository contains code that was used in a 2D Vertex Model simulation of tissue mechanics. This work was included in the manuscript titled, "Zasp52 strengthens whole embryo tissue integrity through supracellular actomyosin networks", Ashour et. al. To be submitted in 2022. 

The following code has been tested on macOS Monterey 12.5 and Ubunutu 20.04.  Earlier versions of each OS should work, but have not been tested.  Windows should be able to get the code working with a little effort. 

## Dependencies 
Python 3.7


## Author:
* [Clinton H. Durney](https://clintondurney.github.io/) (clinton.durney@jic.ac.uk)

## To use:
__Installation__
Use [Conda](https://docs.conda.io/en/latest/) to create the environment from the .yaml file:
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

## Citation:
If you find the code useful, please consider citing either or both of the following:
1. [Durney et. al. (2018)](https://www.sciencedirect.com/science/article/pii/S0006349518311615)
2. [Durney et. al. (2021)](https://iopscience.iop.org/article/10.1088/1478-3975/abfa69/meta)








