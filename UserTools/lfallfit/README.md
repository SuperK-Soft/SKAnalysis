# lfallfit

This is a Tool wrapper for the various different lfallfit functions.
It chooses an appropriate version based on the SK geometry (phase), and whether the input is MC or data.
For SK-IV, it currently uses `lfallfit_sk4_final_qe43_mc_` for MC and `lfallfit_sk4_final_qe43_` for data.
These are currently hard-coded, so comment out the one you want (or add a config variable to control it).



## Data
Hits from the `skt_`, `skq_` and `skchnl_` common block, which should be populated for any kind of SK file.
Also relies on the water transparency, which should be automatically updated by the TreeReader. For MC, you should provide a a suitable value for the TreeReader config file variable `referenceWaterRun`.

## Configuration
* `verbosity`:
	how verbose to be
* `readerName`:
	name of a TreeReader Tool which is reading the input events to reconstruct.
* `hitLimitForClusfit` (aka NHITCUT):
	Maximum number of hits for clusfit to be attempted. Above this, lfallfit skips the clusfit pre-fit step and goes straight to bonsai.
* `StepsToPerform` (aka flag_skip in lfallfit):
	0 (default): do everything (vertex, direction, distance to wall, energy, clik pattern likelihood, etc.
	1: just recalculate energy and related variables
	2: recalculate energy, clik, msg and related variables

