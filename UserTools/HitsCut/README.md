# HitsCut

## Initialise
Gets the only config variable, verbosity.

## Execute
Pulls the total hits for the event from the `skq` common block. Note that the `skq` common block needs to be populated for this tool to work. Currently it does not search to find a populate `nhits` variable.

If the number of hits in the inner detector for the event is larger than a config variable `hitLimit` then the "Skip" variable is set in the datamodel so that the rest of the toolchain is skipped for this event.

The config variable `hitLimit` is set to 999 by default.

## Finalise
Nothing is done here. 
