#include "SpallCandidates.h"
#include "Constants.h"

#include "fortran_routines.h"
#include "softtrg_tbl.h"
#include "SK_helper_functions.h"

#include <bitset>

SpallCandidates::SpallCandidates():Tool(){}

// Selects LE events within +-60 s of a HE event. Creates a spallation selection by minusing the pre sample from
// the post sample.

bool SpallCandidates::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("treeReaderName",treeReaderName);
	
	// if getting data from TTree, check the TreeReader
	if(m_data->Trees.count(treeReaderName)==0){
		Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
		return false;
	} else {
		myTreeReader = m_data->Trees.at(treeReaderName);
	}
	
	return true;
}


bool SpallCandidates::Execute(){
	
	/*
		TODO List
		The current reduction places a +-30 s window around the candidates. Marcus mentioned a +-60 s window instead.
		Searches through all the event windows for muons
		Fit muons using lfmufit_sk4_hz. This is something that could be ported into tool.
		Calculated dE/dx using various methods
		Replaces those with dE/dx results from the fits. I don't know why it calculates these in the first place,
		maybe for the fitters? Seems to be an issue in the fitters that the ionisation energy is not being calculated
		correctly.
	*/
	
	// get HEADER branch
	myTreeReader->Get("HEADER", myHeader);
	triggerID = myHeader->idtgsk;
	
	int idetector [32], ithr [32], it0_offset [32],ipret0 [32],ipostt0 [32];
	
	softtrg_get_cond_(idetector,ithr,it0_offset,ipret0,ipostt0);
	
	
	for(int i = 0; i < 32; i++){
		if(i != 1 && i != 3){
			ithr[i] = 100000;
			it0_offset[i]=0;
			ipret0[i]=0;
			ipostt0[i]=0;
		}
	}
	softtrg_set_cond_(idetector,ithr,it0_offset,ipret0,ipostt0);
	
	// Need to call softtrg_inittrgtbl_ to populate the swtrgtbl_ common block. This requires int* arguments.
	int normalrunnum = myHeader->nrunsk;
	int* runnum = &normalrunnum;
	int tempo = 1280;
	int* max_qb = &tempo;
	int one = 1;
	int* onepointer = &one;
	int zero = 0;
	int* zeropointer = &zero; 
	
	int nhwtrg = softtrg_inittrgtbl_(runnum, zeropointer, onepointer, max_qb);
	
	int ntrg = softtrg_inittrgtbl(myHeader->nrunsk, 0, 1, 1280);
	
	int NMuons = 0;
	std::vector<float> untaggedMuonTime;
	
	//check for muons (tagged and untagged)
	for(int i = 0; i < ntrg; i++){
		if(swtrgtbl_.swtrgtype[i] == 1){
			for(int j = 0; j < ntrg; j++){
				if(swtrgtbl_.swtrgtype[j] == 3){
					if(abs(swtrgtbl_.swtrgt0ctr[j] - swtrgtbl_.swtrgt0ctr[i]) < 100){
						NMuons++;
						untaggedMuonTime.push_back(swtrgtbl_.swtrgt0ctr[i]);
						m_data->vars.Set("newMuon", true);
					}
				}
			}
		}
	}
	
	// if HE (1) and OD (3) trigger bits are already set then the event is an identified muon.
	if(triggerID.test(1) && triggerID.test(3)){
		m_data->vars.Set("newMuon", true);
		if(NMuons >= 2 && untaggedMuonTime[NMuons-1] > 400){
			std::cout << "Warning: DiMuon event found!" << std::endl;
		}
	}else if(NMuons >= 1){
		std::cout << "Untagged muon found." << std::endl;
	}
	
	return true;
	
}


bool SpallCandidates::Finalise(){
	
	return true;
}
