# PINT
Comprehensive Algorithm for Tracking Convective Features in Meteorological Data

PINT is an R implementation of the efficient and adaptive algorithm for tracking convection in radar, satellite, and simulated data. It uses the Fourier phase shift method to compute the first guess of the motion and the Hungarian method for optimized matching. The algorithm also checks merging and splitting of the objects based on their size and overlap with the other objects. The algorithm works with bi-level images. Therefore,  the tracks are estimated independent of the physical properties of the objects, such as brightness temperature or reflectivity. TINT is not Titan is the python versin based on the PINT.
https://github.com/openradar/TINT


The package is made from the R script located at https://github.com/RBhupi/Darwin-Rscripts/blob/master/echo_tracking.R

