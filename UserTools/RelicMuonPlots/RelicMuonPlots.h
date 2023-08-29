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
	
	HistogramBuilder hb;
	std::string outputFile="relicmuplots.root";
	
	std::string muReaderName="";
	std::string relicReaderName="";
	MTreeReader* relicReader=nullptr;
	MTreeReader muReader;
	Header* relicHeader=nullptr;
	Header* muHeader=nullptr;
	MuInfo* muMuInfo=nullptr;
	LoweInfo* relicLowe=nullptr;
	std::vector<int>* muEvNums=nullptr;
	float dt;
	float dll;
	float dlt;
	
};


#endif
