# Configure files

***********************
#Description
**********************

This is a simple demonstration toolchain to show usage of the TreeReader tool when operating on SKROOT files, in order to call lf_allfit.
The toolchain consists of three tools - the LoadFileList tool handles generic specification of a set of input files. This can include a single file, a text file containing a list of files, or a glob or regex pattern. The TreeReader tool provides an interface to the ROOT file, handling entry loading and the population of fortran common blocks. Finally the lf_allfit_new tool performs lf_allfit.

************************
#Usage
************************

In LoadFileListConfig specify the input files. In the TreeReaderConfig give the name of the tree to read, the name to give to the TreeReader, and the various other options for handling the file. Many options have suitable defaults, but it is worth familiarising yourself with the full set of controls particularly for SKROOT files. In lf_allfit_newConfig specify the name of the TreeReader, which is used to retreive the logic unit number associated with the file, needed for the various skroot function calls within this Tool.
