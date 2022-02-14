#include <skparmC.h>
#include <sktqC.h>

#include "SubtractToF.h"
#include "TVector3.h"

bool SubtractToF::Initialise(std::string configfile, DataModel &data)
{
	if(configfile!="")  m_variables.Initialise(configfile);
	m_data= &data;
	m_data->tool_configs[name] = &m_variables;
	
    return true;
}

bool SubtractToF::Execute()
{
    TVector3 promptVertex;
    m_data->eventVariables.Get("prompt_vertex", promptVertex);

    // ToF is automatically subtracted by setting vertex
    m_data->eventPMTHits.SetVertex(promptVertex);
    m_data->eventPMTHits.Sort();

    return true;
}

bool SubtractToF::Finalise()
{
    return true;
}
