The RunEnergyCut Tool places a run-wise cut on the bsenergy.
Cuts should be specified in the config file in the following way:

	cut runmin runmax energymin energymax

* The first word should be 'cut' and identifies a cut line.
* runmin and runmax specify the range of run numbers over which the cut is active.
* energymin and energymax specify the range of energies to ACCEPT.

i.e. reject any event for which the following is true:

	(runmin <= run <= runmax) && NOT (energymin < bsenergy < energymax)

Use a value of -1 to indicate a value of 'any'. For example, to specify a cut which is active for any runs up to and including run 68670, that would reject any events with bsenergy less than 10MeV (but with no upper limit), use a line of:

	cut -1 68670 10 -1

* Any line in the config file matching this format will be added to the set of cuts applied.
* Any number of cuts can be applied.
* Any other configuration variables (.e.g verbosity) can also be included and will be used as normal.
