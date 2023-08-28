# SkipTriggers
A simple Tool to reject a set of triggers. Specify a set of allowed triggers with:

	allowedTriggers 1,2,3

to only allow events where trigger bits 1, 2 or 3 are set and reject any events where none of these bits are set. Alternatively specify a set of rejected triggers with:

	skippedTriggers 11,12,13

to skip events where trigger bits 11, 12 or 13 are set. Trigger bits are retrieved from skhead_.idtgsk.

Triggers may be specified by their bit number, or by their name:
-	LE=0
-	HE=1
-	SLE=2
-	OD_or_Fission=3
-	Periodic=4
-	AFT_or_Cal=5
-	Veto_Start=6
-	Veto_Stop=7
-	unknown_8=8
-	unknown_9=9
-	unknown_10=10
-	Random_Wide=11
-	ID_Laser=12
-	LED=13
-	Ni=14
-	OD_Laser=15
-	LE_hitsum=16
-	HE_hitsum=17
-	SLE_hitsum=18
-	OD_hitsum=19
-	unknown_20=20
-	unknown_21=21
-	SN_Burst=22
-	mue_Decay=23
-	LINAC=24
-	LINAC_RF=25
-	unknown_26=26
-	Periodic_simple=27
-	SHE=28
-	AFT=29
-	Pedestal=30
-	T2K=31

Note these options may also be given directly to the TreeReader Tool. Using this Tool permits more flexibility in the placement of such selections. For example, a simple ToolChain may contain:

	TreeReader -> PlotQismsk -> SkipTriggers -> PlotQismsk

The hypothetical PlotQismsk Tool makes a plot of qismsk.
The SkipTriggers rejects all events where the SHE trigger bit is not set.
This toolchain would then allow one to compare the distribution of qismsk for all events vs those for which the SHE trigger bit is set.

