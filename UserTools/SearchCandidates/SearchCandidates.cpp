#include "SearchCandidates.h"
#include "Candidate.h"


bool SearchCandidates::Initialise(std::string configfile, DataModel &data)
{
	if(configfile!="")  m_variables.Initialise(configfile);
	m_data= &data;
	m_data->tool_configs[name] = &m_variables;
	m_variables.Get("m_verbose",m_verbose);
	
    // read ntag options
    if (!m_variables.Get("T0TH", T0TH)) T0TH = 2000;
    if (!m_variables.Get("T0MX", T0MX)) T0MX = 535000;
    if (!m_variables.Get("TWIDTH", TWIDTH)) TWIDTH = 14;
    if (!m_variables.Get("TMINPEAKSEP", TMINPEAKSEP)) TMINPEAKSEP = 14;
    if (!m_variables.Get("NHITSTH", NHITSTH)) NHITSTH = 7;
    if (!m_variables.Get("NHITSMX", NHITSMX)) NHITSMX = 70;
    if (!m_variables.Get("N200TH", N200TH)) N200TH = 0;
    if (!m_variables.Get("N200MX", N200MX)) N200MX = 200;

    return true;
}

bool SearchCandidates::Execute()
{
    m_data->eventCandidates.Clear();
    
    // eventPMTHits and prompt vertex must be filled
    if (m_data->eventPMTHits.IsEmpty()
        || !(m_data->eventVariables.Get("prompt_vertex", promptVertex))) {
        Log("EventPMTHits is empty, or prompt vertex is not set.", pWARNING, m_verbose);
        //m_data->vars.Set("Skip",true);
        return true;
    }
    
    //Log("Dumping hits");
    //m_data->eventPMTHits.DumpAllElements();

    int   iHitPrevious    = 0;
    int   NHitsNew        = 0;
    int   NHitsPrevious   = 0;
    int   N200Previous    = 0;
    float t0Previous      = -1e6;

    PMTHitCluster* eventHits = &(m_data->eventPMTHits);
    unsigned long nEventHits = eventHits->GetSize();

    // Loop over the saved TQ hit array from current event
    for (int iHit = 0; iHit < nEventHits; iHit++) {

        PMTHitCluster hitsInTWIDTH = eventHits->Slice(iHit, TWIDTH);

        // If (ToF-subtracted) hit comes earlier than T0TH or later than T0MX, skip:
        float firstHitTime = hitsInTWIDTH[0].t();
        if (firstHitTime < T0TH || firstHitTime > T0MX) continue;

        // Calculate NHitsNew:
        // number of hits within TWIDTH (ns) from the i-th hit
        int NHits_iHit = hitsInTWIDTH.GetSize();

        // Pass only if NHITSTH <= NHits_iHit <= NHITSMX:
        if (NHits_iHit < NHITSTH || NHits_iHit > NHITSMX) continue;

        // We've found a new peak.
        NHitsNew = NHits_iHit;
        float t0New = firstHitTime;

        // Calculate N200
        PMTHitCluster hitsIn200ns = eventHits->Slice(iHit, TWIDTH/2.-100, TWIDTH/2.+100);
        int N200New = hitsIn200ns.GetSize();

        // If peak t0 diff = t0New - t0Previous > TMINPEAKSEP, save the previous peak.
        // Also check if N200Previous is below N200 cut and if t0Previous is over t0 threshold
        if (t0New - t0Previous > TMINPEAKSEP) {
            if (N200Previous < N200MX && t0Previous > T0TH) {
                Candidate candidate(iHitPrevious);
                m_data->eventCandidates.Append(candidate);
            }
            // Reset NHitsPrevious,
            // if peaks are separated enough
            NHitsPrevious = 0;
        }

        // If NHits is not greater than previous, skip
        if ( NHitsNew <= NHitsPrevious ) continue;

        iHitPrevious  = iHit;
        t0Previous    = t0New;
        NHitsPrevious = NHitsNew;
        N200Previous  = N200New;
    }
    // Save the last peak
    if (NHitsPrevious >= NHITSTH) {
        Candidate candidate(iHitPrevious);
        m_data->eventCandidates.Append(candidate);
    }

    return true;
}

bool SearchCandidates::Finalise()
{
    return true;
}
