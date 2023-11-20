#ifndef RelicMuonPlots_H
#define RelicMuonPlots_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "HistogramBuilder.h"
#include "MTreeReader.h"

/**
* \class RelicMuonPlots
*
* Plot distributions for muon-relic pairs
*
* $Author: M.O'Flaherty $
* $Date: 2023/27/08 $
* $Contact: moflaher@km.icrr.u-tokyo.ac.jp
*/

class RelicMuonPlots: public Tool {
	
	public:
	
	RelicMuonPlots();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	bool MakeHists(int step);
	bool MakePairVariables();
	bool GetRelicEvt();
	bool GetMuonEvt(int entrynum);
	double CalculateTrackLen(float* muon_entrypoint, float* muon_direction, double* exitpt=nullptr);
	bool GetPeTable();
	
	// for histogramming
	HistogramBuilder hb;
	HistogramBuilder relic_hb;
	HistogramBuilder mu_hb;
	std::string outputFile="relicmuplots.root";
	std::string relicFile="relicplots.root";
	std::string muFile="muonplots.root";
	std::map<int,int> muon_plotted;
	std::map<int,int> relic_nevsks; // FIXME this shouldn't be needed!
	std::map<std::pair<int,int>, int> muon_pair_plotted;
	
	// table that converts coulombs to photoelectrons
	// normally in $SKOFL_ROOT/const/lowe/petable.dat
	std::string pe_table_filename = "petable.dat";
	// rows read from this file held in memory
	std::vector<std::array<float,3>> pe_to_coulombs;
	std::vector<int> petable_startrun, petable_endrun;
	
	// for reading Trees
	int relicentrynum=0;
	std::string muReaderName="";
	std::string relicReaderName="";
	MTreeReader* relicReader=nullptr;
	MTreeReader* muReader=nullptr;
	
	// for getting branches
	// relic:
	Header* relicHeader=nullptr;
	LoweInfo* relicLowe=nullptr;
	TQReal* relicTQReal=nullptr;
	TQReal* relicTQAReal=nullptr;
	int64_t relicClockTicks=0;
	int relicRollovers=0;
	std::vector<int>* relicMatchedEntryNums=nullptr;
	std::vector<float>* relicTimeDiffs=nullptr;
	
	// muon:
	Header* muHeader=nullptr;
	MuInfo* muMu=nullptr;
	int64_t muClockTicks=0;
	int muRollovers=0;
	TQReal* muTQReal=nullptr;
	TQReal* muTQAReal=nullptr;
	
	// for tracking rates:
	int64_t lastmuticks=0;
	int lastmu_nevsk=0;
	int64_t lastrelicticks=0;
	
	// for making spallation observables
	float* muon_dedx;
	float* muon_entrypoint;
	float* muon_direction;
	double muon_tracklen;
	
	int pe_table_index = 0;
	
	int mu_before_relic;
	float dt;
	float dll;
	float dlt;
	
};


#endif
