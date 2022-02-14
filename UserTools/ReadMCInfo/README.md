This Tool appears to have no dependancies on upstream tools.
It populates the eventPrimaries, eventSecondaries, and eventTrueCaptures members of the DataModel,
by reading either the MCVERTEX, MCVECT and SCNDPRT and SCNDPRTVC data banks for zbs files,
or the SECONDARY branch's atmpd-style arrays for SKROOT files - Note these arrays are only populated
when generating files with SKG4, and explicitly given the XXX flag. (skdetsim does not produce this branch).
The WriteOutput Tool will then put these into the 'mc' TTree in the output file, although the WriteOutput tool also writes out other output, so has other Tool dependencies.

T0 is obtained by: `trginfo_(&T0)` and used for all primary particle times, and added to all secondary particle times.

primaries are retrieved by calling `skgetv_();` ($SKOFL_ROOT/src/sklib/skgetv.F).
This works for both zbs and ROOT, provided 'SK_FILE_FORMAT' is set first (0=zbs, 1=root).
For zbs:
	populates the `skvect_` common block with data from MCVERTEX and MCVECT banks (skdetsim primaries).
	n.b. it only returns the last vertex in MCVERTEX in skvect_.pos, no others.
	number of initial particles is skvect_.nvect (capped at 50).
	particle codes are returned in skvect_.ip[0..NVECT-1], initial positions in skvect_.pin[0..NVECT-1][0..2],
	and initial momentum in skvect_.pabs[0..NVECT-1].
	
	populates the `skvect_add_` common block with additional data from the MCPARMCONV bank
	(e.g. dark rate, random seeds, trigger config etc. - i.e. DARKDS, TCNVSK, QCNVSK, DTHRSK etc.)
For root:
	uses fortran interface skroot_get_mc(). Only last primary vertex and primary particles populated.

secondaries are retreieved by either:
for root:
	copying data from the SECONDARY branch
	in particular: iprtscnd (pdg), tscnd (time), vtxscnd (vertex), pscnd (momentum), lmecscnd (creation process), iprntprt (parent pdg) branches are used to construct Particle class instances
N.B. only some secondaries are copied: neutrons and *some hydrogen capture products*:
		-> deuterons, gammas, and electrons over ckv threshold "from an interaction other than multiple scattering" (lmecscnd!=2)
	Note also: particles where the secondary creation vertex was outside the ID (radius > RINTK, height > ZPINTK) or where it was within a PMT (determined via fortran function inpmt_()) are skipped.
	electrons and deuterons are not associated to their respective capture, only gammas are - by their creation time being within 1e-7 ns  of an existing capture....!
for zbs:
	apflscndprt_(); // note the following are not all the same. Not sure which is being used!
	$ATMPD_ROOT/src/analysis/ndecay/analyses/p2muk0/official_ntuple_pdkfit/apflscndprt.F
	$ATMPD_ROOT/src/analysis/neutron/ntag/apflscndprt.F
	$ATMPD_ROOT/src/analysis/neutron/ntag_gd/apflscndprt.F
	$ATMPD_ROOT/src/analysis/official_ntuple/apflscndprt.F
	$ATMPD_ROOT/src/recon/fitqun/apflscndprt.F
	
	based on usage this populates secndprt_ with the same members are detailed below (and possibly more...)

Note: the 'variables' tree 'true_neutron_count' branch is set as the number of neutron captures, as identified by the secondaries scan. A neutron capture is identified by finding a secondary product created by process lmec=18, at a unique time.
