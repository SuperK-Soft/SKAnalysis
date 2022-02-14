# PurewaterSpallAbundanceCuts
This tool applies a series of selection cuts on muon-lowe pair events to produce a sample of spallation events, which will be analysed by subsequent tools to produce distributions of spallation observables, and ultimately to measure the rate of production of different spallation isotopes. This is working toward a reproduction of the results from the 2015 paper on measuring cosmic spallation isotope yields.

## Data
The tool accepts its data from an MTreeReader accessing ROOT files containing muon-lowe event pairs, produced during the 2020 SRN analysis chain. The output is stored as a ROOT file describing the cuts applied in this tool and which TTree entry numbers (and branch indices, where branches store arrays) of those events that passed each cut.
Downstream tools may also access this information during the toolchain by accessing the 'SpallAbundanceSelection' MTreeSelection in the datamodel CStore.

## Configuration
Configuration variables include various cut criteria applied during selection.
```
# Tool use
# --------
verbosity             # how verbose to be
treeReaderName        # name of the input TreeReader
outputFile            # name of output ROOT file

# selection variables
# -------------------
max_closest_muon_dt   # reject events where the closest muon (preceding or following) is within this time [s]
max_closest_lowe_dx   # reject events where the next closest lowe event within 60s within this distance [cm]
ntag_FOM_threshold    # minimum neutron tagging FOM to identify an aftertrigger as a neutron capture
run_min               # minimum run number
run_max               # maximum run number
```
