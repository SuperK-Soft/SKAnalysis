# CombinedFitter

Need to divide up the code and have some classes within the fitter to:

1. Get the MC info where applicable
2. Get the TQ info for the SHE
3. Add the neutron-capture events (find all triggers after SHE and add them, but have them related to the SHE) (this should be a separate bit of code. Can I adapt what the sk2p2mev code does?
4. Do the single-event fit
5. Do the combined fit for each possible combination of SHE+AFT
6. Select the SHE+AFT combo with the best likelihood or otherwise?
