Fix cernlib 64-bit address space errors.
The fix requires: lib/liblowfit_sk4_stripped.so to be linked into 'main', together with          
a "suitable" ordering of arguments to the corresponding g++ line for 'main'.
(i.e. it is not sufficient to add `-llowfit_sk4_stripped` to the end of the line, 
although the specifics of which order is required has not been fully explored).
This library provides the function `do_lowfit`. It is not necessary to invoke this function,    
but something about its presence when linked early on appears to alter the way in which 
cernlib itself is linked in as the error does not manifest. Since the function is unused 
the combination of CXXFLAGS `-O3 -Wl,--as-needed` will result in it being stripped out,
and a resurfacing of the problem. The result is that 'main' will compile, but calls to
fortran routines produce an errors of the form:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOCB/LOCF: address 0x7fb881fd3200 exceeds the 32 bit address space
or is not in the data segments
This may result in program crash or incorrect results
Therefore we will stop here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
and the applcation immediately terminates.
