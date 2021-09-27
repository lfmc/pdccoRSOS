# pdccoRSOS

This is repository contains C++ source code and associated SLURM bash scripts used in the paper "Parallel generation of extensive vascular networks with application to an archetypal human kidney model".
It also contains python scripts used for the analysis and images of this same work.
For the actual input data and code output please check the Dryad repository.

## Section 3a

This section analyses how the number of segments in the baseline tree affects the final network.

## Section 3b

This section analyses if the subdomain's shape functional affects the final network.

## Section 3c

This section analyses how the algorithm copes when the domain is supplied by multiple inlets.

## Section 3d

This section analyses the scalability of the PDCCO approach.

## Section 3e

This section illustrates the application of the PDCCO to the vascularization of a prototypical human kidney.

# Needed software

## C++11

## CMake 

Used for building **VTK** and **VItA** libraries. 
Available at https://cmake.org/download/.

## VTK 8.1.2

[The Visualization Toolkit v8.1.2](https://gitlab.kitware.com/vtk/vtk/-/tree/v8.1.2) is used by the **VItA** library.

## VItA 

[The Virtual ITerative Angiogenesis](https://github.com/lfmc/VItA/tree/Cury2021_et_al_PDCCO) implements the PDCCO and [DCCO](https://doi.org/10.1038/s41598-021-85434-9) methods. For the main version please visit https://github.com/GonzaloMaso/VItA.

## Plotting and related libraries
[Python 3.8.5](python.org/downloads/).
[Pandas 1.1.3](https://pypi.org/project/pandas/1.1.3/#files).
[Numpy 1.19.2](https://numpy.org/install/).
[Matplotlib 3.3.2](https://matplotlib.org/stable/users/installing.html).
