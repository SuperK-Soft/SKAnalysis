#ifndef ReconstructMatchedMuons_H
#define ReconstructMatchedMuons_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "MTreeReader.h"
#include "ParticleCand.h"

class ReconstructMatchedMuons: public Tool {
	
	
	public:
	
	ReconstructMatchedMuons(); ///< Simple constructor
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resorces. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Executre function used to perform Tool perpose. 
	bool Finalise(); ///< Finalise funciton used to clean up resorces.
	
	bool ReconstructNextMuon();
	bool WriteEventsOut(std::vector<ParticleCand>& eventsToWrite, int outLUN, EventType eventType);
	bool AddAftHits(const rawtqinfo_common& rawtqinfo_aft);
	
	private:
	
	MTreeReader* rfmReader = nullptr;
	std::string rfmReaderName;
	int rfmReaderLUN=-1;        // where to put output results (skroot_set_mu)
	int relicWriterLUN=-1;
	int muWriterLUN=-1;
	bool noBFF=false;  // veto fallback to BFF
	
	std::vector<int> MatchedEvNums;
	std::vector<int> MatchedEntryNums;
	std::vector<float> MatchedTimeDiff;   // [ns]
	std::vector<float> MatchedParticleE;  // [MeV]
	std::vector<bool> MatchedHasAFTs;
	uint64_t HwClockTicks;
	int NumRollovers;
	
	int relics_to_write=0;
	int muons_to_write=0;
	int relics_written=0;
	int muons_written=0;
	int muons_written_wmuboysplit=0;
	
	int current_badch_masking;
	std::vector<skroot_mu_common> reco_muons;
	
};


#endif
