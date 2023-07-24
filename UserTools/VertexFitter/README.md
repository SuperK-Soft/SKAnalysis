# Vertex fitter

Single-event fitter replicates lf_allfit with the option to use direct bonsai functions or via the bonsaifit function in bscalls.cc. 

Can choose hit data from tqreal branch, or SK common blocks so there is some additional flexibility. Sample config file in configfiles/Fitter. 

Also calculates the rest of the lowe variables but still using the Fortran functions - note that some of these do not save correctly to the LOWE.linfo array.

