# tracR
Adoptive Tracking Algorithm for Convection Features in Meteorological Data (in R)

tracR is an implementation of an efficient and adaptive algorithm for tracking convection in radar,
satellite, and simulated data. It uses the Fourier phase shift method to compute the first guess of the
motion and the Hungarian method for optimized matching. The algorithm checks for merging and
splitting of the objects based on their size and overlapping with the other objects. The algorithm
works with bi-level images; hence the tracks are estimated independent of the objects’ physical prop-
erties, such as brightness temperature or reflectivity. "TINT is not Titan" is a python version located
at https://github.com/openradar/TINT. Both tracR and TINT packages are an independent adap-
tion of the prototype (echo_tracking.R) located at https://github.com/RBhupi/Darwin-Rscripts/

Reference for citing the tracking algorithm.
Raut, B. A., Jackson R., Picel M., Collis S., Bergemann M. & Jakob C.
(2021). A tracking algorithm for convection in simulated and remote
sensing data. *J. Appl. Meteorol. Clim.*, 60(4), 513-526. https://doi.org/10.1175/JAMC-D-20-0119.1
