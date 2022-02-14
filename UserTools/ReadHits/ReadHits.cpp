#include <skheadC.h>
#include <skparmC.h>
#include <sktqC.h>

#include "ReadHits.h"

bool ReadHits::Initialise(std::string configfile, DataModel &data)
{
	if(configfile!="")  m_variables.Initialise(configfile);
	m_data= &data;
	m_data->tool_configs[name] = &m_variables;
    return true;
}

bool ReadHits::Execute()
{
    m_data->eventVariables.Set("nevsk", skhead_.nevsk);
    
    prevEvTrigType = currentEvTrigType;
    SetTriggerType();
    
    PMTHitCluster* eventHits = &(m_data->eventPMTHits);
    Log("Clear event hits");
    eventHits->Clear();
    
    // Check if sktqz is filled
    if (!sktqz_.nqiskz) {
        Log("ReadHits: The SK common block sktqz is empty, skipping event", pWARNING, m_verbose);
        //m_data->vars.Set("Skip",true);
        return true;
    }
    
    // Copy hits from sktqz_ common block to eventHits in the DataModel
    // --------------------------------------------------------------------
    Log("Recording prompt event hits");
    for (int iHit = 0; iHit < sktqz_.nqiskz; iHit++) {
        
        PMTHit hit{ /*T*/ sktqz_.tiskz[iHit],
                    /*Q*/ sktqz_.qiskz[iHit],
                    /*I*/ sktqz_.icabiz[iHit]
                  };
        
        // Use in-gate hits only
        if (sktqz_.ihtiflz[iHit] & (1<<1)) {
            //hit.Dump();
            eventHits->Append(hit);
        }
    }
    
    // if this is an SHE with AFT and we loaded both together,
    // also append the AFT hits
    if(currentEvTrigType == tSHE && m_data->HasAFT()) {
        // load the common  blocks with the AFT hits
        m_data->LoadAFT();
        
        // it seems what we do here is try to find the last SHE hit also in the set of AFT hits?
        // That hit is then used as a common point of reference to align the time axis of the events.
        PMTHit lastHit(0, 0, 0);
        lastHit = eventHits->GetLastHit();   // get the last SHE hit
        float tOffset = 0.;
        int iHit;  // carry this over, so we only start copying from the next hit after the matching one
        bool coincidenceFound=false;
        for (iHit = 0; iHit < sktqz_.nqiskz; iHit++) {
            // match the last hit by same cable number and charge
            if (sktqz_.qiskz[iHit] == lastHit.q() && sktqz_.icabiz[iHit] == lastHit.i()) {
                tOffset = lastHit.t() - sktqz_.tiskz[iHit];
                coincidenceFound = true;
                Log(Form("coinciding hit: t: %f ns q: %f i: %d", sktqz_.tiskz[iHit], sktqz_.qiskz[iHit], sktqz_.icabiz[iHit]));
                Log(Form("Coincidence found: t = %f ns, (offset: %f ns)", lastHit.t(), tOffset));
                break;
            }
        }
        if(!coincidenceFound){
            Log("Error! ReadHits could not find a coincident hit between SHE and following AFT events! "
                "Skipping the addition of AFT hits!",0,1);
            prevEvTrigType = currentEvTrigType;  // make sure this gets set before we bail
            return false;
        } else {
            //Log("Before appending hits");
            //eventHits->DumpAllElements();
            
            Log("Appending AFT hits");
            for (iHit = iHit; iHit < sktqz_.nqiskz; iHit++) {
                PMTHit hit{ /*T*/ sktqz_.tiskz[iHit] + tOffset,
                            /*Q*/ sktqz_.qiskz[iHit],
                            /*I*/ sktqz_.icabiz[iHit]
                          };
                // Use in-gate hits only
                if (sktqz_.ihtiflz[iHit] & (1<<1)) {
                    //hit.Dump();
                    eventHits->Append(hit);
                }
            }
        }
        
        eventHits->Sort();
        
        //Log("After appending hits");
        //eventHits->DumpAllElements();
        
    }  // else not SHE with AFT, no need to add AFT hits.
    
    std::cout << "\n";
    Log("Finalised hits");
    
    return true;
}

bool ReadHits::Finalise()
{
    return true;
}

void ReadHits::SetTriggerType()
{
    // SHE
    if (skhead_.idtgsk & 1<<28)
        currentEvTrigType = tSHE;
    // AFT
    else if (skhead_.idtgsk & 1<<29)
        currentEvTrigType = tAFT;
    // Else
    else
        currentEvTrigType = tELSE;
    
    m_data->eventVariables.Set("trigger_type", (int)(currentEvTrigType));
    Log(Form("Previous trigger type: %d Current trigger type: %d", prevEvTrigType, currentEvTrigType));
}
