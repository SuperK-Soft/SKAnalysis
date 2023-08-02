# Configure files

***********************
#Description
**********************

This is a simple demonstration toolchain to show usage of the TreeReader tool.
The toolchain consists of two tools - the TreeReader tool that performs interfacing with the ROOT file, and a TreeReaderDemo tool that shows a downstream tool retrieving objects from the TreeReader.

************************
#Usage
************************

In TreeReaderConfig specify the input path, the name of the tree to read, and the name to give to the TreeReader.
The input path may either be a file name, or a TChain glob pattern. The TreeReader name is used by downstream
tools to retrieve the reader from the DataModel.
Additional optional configuration variables include:
* A list of branches to read. If given, all other branches will be de-activated to improve performance, otherwise all branches will be left active.
* A starting TTree entry number. If not given, start from the first entry.
* The maximum number of entries to read. If given the toolchain will stop after reading at most that many entries, or upon reaching the end of the TTree/TChain, whichever is first. If not given (or negative), all entries in the TTree or TChain will be processed.
