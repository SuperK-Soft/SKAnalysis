#include "RelicMuonMatching.h"
#include "fortran_routines.h"
#include "Constants.h"
#include "skheadC.h"
#include "ParticleCand.h"
#include <inttypes.h>
#include <algorithm>

RelicMuonMatching::RelicMuonMatching():Tool(){}


bool RelicMuonMatching::Initialise(std::string configfile, DataModel &data){
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;
	
	m_variables.Get("treeReaderName", treeReaderName);
	
	if(m_data->Trees.count(treeReaderName)==0){
	Log("Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,verbosity);
	return false;
	} else {
		myTreeReader = m_data->Trees.at(treeReaderName);
	}
	
	lun = 20;
	
	TreeManager* mgr = skroot_get_mgr(&lun);
	WriteTree = mgr->GetOTree();
	
	MatchedEvNumsBranch = WriteTree->Branch("MatchedEvNums", &MatchedEvNums);
	MatchedTimeDiffBranch = WriteTree->Branch("MatchedTimeDiff", &MatchedTimeDiff);
	PIDBranch = WriteTree->Branch("MuonTag", &PID);
	
	return true;
}


bool RelicMuonMatching::Execute(){
	
	myTreeReader->Get("HEADER", myHeader);
	myTreeReader->Get("LOWE", myLowe);
	
	m_data->vars.Get("newMuon", muonFlag);
	m_data->vars.Get("newRelic", relicFlag);
	
	if(myHeader->nsubsk != currentSubRun){
		currentSubRun = myHeader->nsubsk;
		std::cout << "Subrun number:            " << myHeader->nsubsk << std::endl;
	}
	
	/*
	Event is a muon candidate! Now need to:
	1)  Check if there have not been any relic candidates yet (i.e. the first event is a muon candidate) then 
	    save the muon candidate's event # and time to the deque muonCands
	2)  If there are already relic candidates written to the relicCands vector:
	2a) Check the time difference between the muon and the oldest relic candidate - if the time difference
	    is larger than 30s and if the relic already has muons tagged to it then set the relic to be written
	    out as we have found all of the muons within +-30s of it.then set the relic, and its set of muons as
	    they have all been tagged to a relic candidate, to be written out.
	
	ParticleCand is a struct, defined in ParticleCand.h, with structure:
		struct ParticleCand {
			int EventNumber;
			int EventTime;
			std::vector<int> matchedParticleEvNum;
		};
	*/
	
	if(muonFlag){
		RelicMuonMatch("newMuon", m_data->muonCandDeque, m_data->relicCandDeque);
		m_data->vars.Set("newMuon", false);
	}else if(relicFlag){
		RelicMuonMatch("newRelic", m_data->relicCandDeque, m_data->muonCandDeque);
		m_data->vars.Set("newRelic", false);
	}else{
		Log("Event is not a relic candidate or a muon candidate! Something has gone wrong further up in the chain"
		, v_error, verbosity);
		return false;
	}
	
	if(m_data->writeOutRelics.size()){
		WriteRelicInfo();
	}
	if(m_data->muonsToRec.size()){
		WriteMuonInfo();
	}
	
	if(muonsToRemove.size() > 0){
		RemoveFromDeque(muonsToRemove, m_data->muonCandDeque);
	}
	if(relicsToRemove.size() > 0){
		RemoveFromDeque(relicsToRemove, m_data->relicCandDeque);
	}
	
	return true;
}


bool RelicMuonMatching::Finalise(){
	
	return true;
}

bool RelicMuonMatching::RemoveFromDeque(std::vector<int>& particlesToRemove, std::deque<ParticleCand>& particleDeque){
	for(int j = particlesToRemove.size() - 1; j >= 0; j--){
		for(int i = particleDeque.size() - 1; i >= 0; i--){
			if(particleDeque[i].EntryNumber == particlesToRemove[j]){
				particleDeque.erase(particleDeque.begin() + i);
			}
		}
	}
	particlesToRemove.clear();
	return true;
}

bool RelicMuonMatching::RelicMuonMatch(std::string particleType, std::deque<ParticleCand>& currentDeque, std::deque<ParticleCand>& targetDeque){
	float currentTime = myHeader->counter_32 * 32768./1.92;
	ParticleCand currentParticle;
	currentParticle.EventNumber = myHeader->nevsk;
	currentParticle.EventTime = currentTime;
	currentParticle.EntryNumber = myTreeReader->GetEntryNumber();
	currentParticle.LowECommon = skroot_lowe_;
	
	for(int i = 0; i < targetDeque.size(); i++){
		ParticleCand targetCand = targetDeque[i];
		//find time to each target candidate (if the event is a muon cand then the target is a relic cand) oldest 
		//to newest.
		//If the time difference between the two events is less than 60 seconds then "match" the particles.
		timeDiff = (currentTime - targetCand.EventTime);
		if(abs(timeDiff) < 60. * pow(10,9)){
			//add the event # of the current event to the target particle's "matched particle" list and add the
			//if time difference is less than +-30s then the events need investigating to see if they are correlated
			//event # of the target particle to the current particle's "matched particle" list
			currentParticle.matchedParticleEvNum.push_back(targetCand.EventNumber);
			currentParticle.matchedParticleTimeDiff.push_back(timeDiff * -1.);
			currentParticle.matchedParticleBSEnergy.push_back(targetCand.LowECommon.bsenergy);
			
			targetCand.matchedParticleEvNum.push_back(currentParticle.EventNumber);
			targetCand.matchedParticleTimeDiff.push_back(timeDiff);
			targetCand.matchedParticleBSEnergy.push_back(currentParticle.LowECommon.bsenergy);
			
			//If a LowE event matches to a muon that has previously not had any muons match to it then mark
			//that muon for reconstruction
		}else if(timeDiff > 60 * pow(10,9)){
			//if timeDiff is larger than 60 seconds
			if(particleType == "newMuon"){
				int evPID = 10;
				targetCand.PID = evPID;
				relicsToRemove.push_back(targetCand.EntryNumber);
				//if the current event comes more than 0 seconds after any of the other events then no other
				//events can be tagged to it. So check that the older event has been tagged to other events and
				//set it to be written out.
				if(targetCand.matchedParticleEvNum.size()){
					m_data->writeOutRelics.push_back(targetCand);
				}
			}
			if(particleType == "newRelic"){
				int evPID = 1;
				targetCand.PID = evPID;
				muonsToRemove.push_back(targetCand.EntryNumber);
				if(targetCand.matchedParticleEvNum.size()){
					m_data->muonsToRec.push_back(targetCand);
				}
			}
		}
		targetDeque[i] = targetCand;
	}
	
	currentDeque.push_back(currentParticle);
	
	//There are ~2.5 cosmic ray muons interating in SK per second. So if the muon deque becomes larger than 140
	//then is likely that we will get muons that are over 60 seconds old and have not been matched to any
	//particles. If this is the case then the muon does not need to be kept in the deque and can be removed.
	if(particleType == "newMuon" && currentDeque.size() >= 140 && ! targetDeque.size()){
		for(int i = 0; i < currentDeque.size() - 1; i++){
			timeDiff = (currentTime - currentDeque[i].EventTime);
			if(timeDiff > 60.e+9 && ! currentDeque[i].matchedParticleEvNum.size()){
				muonsToRemove.push_back(currentDeque[i].EntryNumber);
			}else{
				break;
			}
		}
	}
	if(! targetDeque.size()){
		return true;
	}
	
	return true;
}


bool RelicMuonMatching::WriteRelicInfo(){
	
	int originalEntry = myTreeReader->GetEntryNumber();
	
	for(int writeEvent = 0; writeEvent < m_data->writeOutRelics.size(); writeEvent++){
		int currentEntry = myTreeReader->GetEntryNumber();
		if(m_data->writeOutRelics[writeEvent].EntryNumber != currentEntry){
			// don't use this anymore. m_data->GetTreeEntry(treeReaderName, entryNum) has now been added and 
			// can be used instead
			//myTreeReader->GetEntry(m_data->writeOutRelics[writeEvent]);
			m_data->getTreeEntry(treeReaderName, m_data->writeOutRelics[writeEvent].EntryNumber);
		}
		int io = 1;
		
		for(int i; i < branchestoSkip.size(); i++){
			skroot_zero_branch_(&lun, &io, branchestoSkip[i].c_str(), branchestoSkip[i].size());
		}
		
		skroot_lowe_ = m_data->writeOutRelics[writeEvent].LowECommon;
		
		
		skroot_set_lowe_(&lun,
						 skroot_lowe_.bsvertex,
						 skroot_lowe_.bsresult,
						 skroot_lowe_.bsdir,
						 skroot_lowe_.bsgood,
						 &skroot_lowe_.bsdirks,
						 skroot_lowe_.bseffhit,
						 &skroot_lowe_.bsenergy,
						 &skroot_lowe_.bsn50,
						 &skroot_lowe_.bscossun,
						 skroot_lowe_.clvertex,
						 skroot_lowe_.clresult,
						 skroot_lowe_.cldir,
						 &skroot_lowe_.clgoodness,
						 &skroot_lowe_.cldirks,
						 skroot_lowe_.cleffhit,
						 &skroot_lowe_.clenergy,
						 &skroot_lowe_.cln50,
						 &skroot_lowe_.clcossun,
						 &skroot_lowe_.latmnum,
						 &skroot_lowe_.latmh,
						 &skroot_lowe_.lmx24,
						 &skroot_lowe_.ltimediff,
						 &skroot_lowe_.lnsratio,
						 skroot_lowe_.lsdir,
						 &skroot_lowe_.spaevnum,
						 &skroot_lowe_.spaloglike,
						 &skroot_lowe_.sparesq,
						 &skroot_lowe_.spadt,
						 &skroot_lowe_.spadll,
						 &skroot_lowe_.spadlt,
						 &skroot_lowe_.spamuyn,
						 &skroot_lowe_.spamugdn,
						 skroot_lowe_.posmc,
						 skroot_lowe_.dirmc,
						 skroot_lowe_.pabsmc,
						 skroot_lowe_.energymc,
						 &skroot_lowe_.darkmc,
						 &skroot_lowe_.islekeep,
						 &skroot_lowe_.bspatlik,
						 &skroot_lowe_.clpatlik,
						 &skroot_lowe_.lwatert,
						 &skroot_lowe_.lninfo,
						 skroot_lowe_.linfo);
						 
		//delete_outside_hits_();
		skroot_set_tree_(&lun);
		WriteInfo(m_data->writeOutRelics[writeEvent]);
		skroot_fill_tree_(&lun);
	}
	
	
	if(originalEntry != myTreeReader->GetEntryNumber()){
		m_data->getTreeEntry(treeReaderName, originalEntry);
	}
	m_data->writeOutRelics.clear();
	
return true;
}

bool RelicMuonMatching::WriteMuonInfo(){
	std::vector<ParticleCand> muonsToRec = m_data->muonsToRec;
	
	if(! muonsToRec.size()) return true;
	
	myTreeReader->Get("HEADER", myHeader);
	
	int currentEntry = myTreeReader->GetEntryNumber();
	
	for(int i = 0; i < muonsToRec.size(); i++){
		if(muonsToRec[i].EntryNumber != currentEntry){
			m_data->getTreeEntry(treeReaderName, muonsToRec[i].EntryNumber);
		}
		
		
		int muyn_org, muynf;
		
		float mbentry [4];
		float mm_entry [36];
		int n_left;
		
		int currentEventNum = myHeader->nevsk;
		
		//Muon reconstruction developed by Tomoeda and Yamaguchi
		mfmuselect_(&skroot_mu_.muentpoint, &skroot_mu_.mudir, &skroot_mu_.mugoodness, &muyn_org);
		
		//muyn == 1 - good fit
		//muyn == 0 - bad fit
		
		if(muyn_org > 0){
			skroot_mu_.muyn = 1;
		}else if(muyn_org < 0){
			skroot_mu_.muyn = 0;
		}else{
			Log("Muyn_org returning as == 0. Not supported yet", v_error, verbosity);
			return false;
		}
		
		//Apply fast fit if mfmuselect has returned a bad fit
		if(skroot_mu_.muyn == 0){
			mffastfast_(&skroot_mu_.muentpoint, &skroot_mu_.mudir, &muynf);
			skroot_mu_.mufast_flag = 1;
		}else{
			skroot_mu_.mufast_flag = 0;
		}
		
		skroot_mu_.muyn = muyn_org;
		if(skroot_mu_.muyn == 0){
			skroot_mu_.muyn = muynf;
		}
		
		int muboy_class, muboy_numtracks, muboy_numleft;
		float muboy_entry [4], muboy_dir [3], muboy_otherentry [36];
		float muboy_tracklength, muboy_goodness;
		
		
		//Apply muboy
		muboy_zbs_(&currentEventNum,
			  &skroot_mu_.muboy_status,
			  &mbentry,
			  &skroot_mu_.muboy_dir,
			  &skroot_mu_.muboy_length,
			  &skroot_mu_.muboy_goodness,
			  &skroot_mu_.muboy_ntrack,
			  &mm_entry,
			  &n_left);
			  
/*		muboy_zbs_(&skhead_.nevsk, &muboy_class, &muboy_entry,
					&muboy_dir, &muboy_tracklength, 
					&muboy_goodness, &muboy_numtracks,
					&muboy_otherentry, &muboy_numleft);
		
		std::cout << "muboy: " << muboy_class << " , " << muboy_dir[0] << " " << muboy_dir[1] << " " << muboy_dir[2] << " , " << muboy_entry[0] << " " << muboy_entry[1] << " " << muboy_entry[2] << " , " << muboy_tracklength << std::endl;*/
		
		for(int track = 0; track < skroot_mu_.muboy_ntrack; track++){
			if(track == 0){
				skroot_mu_.muboy_entpos[0][track] = 3;
				skroot_mu_.muboy_entpos[1][track] = mbentry[1];
				skroot_mu_.muboy_entpos[2][track] = mbentry[2];
				skroot_mu_.muboy_entpos[3][track] = mbentry[3];
			}else{
				skroot_mu_.muboy_entpos[0][track] = mm_entry[4*track - 8];
				skroot_mu_.muboy_entpos[1][track] = mm_entry[4*track - 7];
				skroot_mu_.muboy_entpos[2][track] = mm_entry[4*track - 6];
				skroot_mu_.muboy_entpos[3][track] = mm_entry[4*track - 5];
			}
		}
		
		float muentry[4];
		for(int j = 0; j < 4; j++){
			muentry[j] = skroot_mu_.muboy_entpos[j][0];
		}
		
		
		float watert;
		
		if(myHeader->nrunsk != lastRun){
			int days_to_run_start = skday_data_.relapse[skhead_.nrunsk];
			lfwater_(&days_to_run_start, &watert);
			lastRun = myHeader->nrunsk;
		}
		
		Makededx(muentry,
				  skroot_mu_.muboy_dir,
				  skchnl_.ihcab,
				  skq_.qisk,
				  skt_.tisk,
				  geopmt_.xyzpm,
				  skq_.nqisk,
				  skroot_mu_.muboy_dedx);
		
		makededx_intg_(&muentry,
						&skroot_mu_.muboy_dir,
						&skroot_mu_.muboy_length,
						&skchnl_.ihcab,
						&skq_.qisk,
						&skt_.tisk,
						&geopmt_.xyzpm,
						&sktqz_.nqiskz,
						&skhead_.nrunsk,
						&skroot_mu_.muboy_dedx,
						&sktqz_.ihtiflz,
						&skhead_.nevsk);
		
		if(skroot_mu_.muboy_status == 1 && skroot_mu_.muboy_goodness < 0.4 && *std::max_element(muonsToRec[i].matchedParticleBSEnergy.begin(), muonsToRec[i].matchedParticleBSEnergy.end()) > 12.){
			std::cout << "Starting muon BFF" << std::endl;
			newmufit_(&bffpos, &hpos, &bffgood);
			modd = sqrt( pow((hpos[0]-bffpos[0]),2) + pow((hpos[1]-bffpos[1]),2) + pow((hpos[2]-bffpos[2]),2) );
			for(int j = 0; j < 3; j++){
				skroot_mu_.mubff_entpos[j] = bffpos[j];
				skroot_mu_.mubff_dir[j] = (hpos[j] - bffpos[j])/modd;
			}
			skroot_mu_.mubff_goodness = bffgood;
			std::cout << "Finished muon BFF" << std::endl;
			
			if(skroot_mu_.mubff_goodness > 0.3){
				for(int j = 0; j < 3; j++){
					muentry[j] = skroot_mu_.mubff_entpos[j];
				}
			}
			Makededx(muentry,
					skroot_mu_.muboy_dir,
					skchnl_.ihcab,
					skq_.qisk,
					skt_.tisk,
					geopmt_.xyzpm,
					skq_.nqisk,
					skroot_mu_.muboy_dedx);
					
			makededx_intg_(&muentry,
							&skroot_mu_.muboy_dir,
							&skroot_mu_.muboy_length,
							&skchnl_.ihcab,
							&skq_.qisk,
							&skt_.tisk,
							&geopmt_.xyzpm,
							&sktqz_.nqiskz,
							&skhead_.nrunsk,
							&skroot_mu_.muboy_dedx,
							&sktqz_.ihtiflz,
							&skhead_.nevsk);
		}
		
		
		skroot_set_mu_(&lun, 
			skroot_mu_.muentpoint, 
			skroot_mu_.mudir, 
			&skroot_mu_.mutimediff, 
			&skroot_mu_.mugoodness,
			&skroot_mu_.muqismsk, 
			&skroot_mu_.muyn, 
			&skroot_mu_.mufast_flag, 
			&skroot_mu_.muboy_status,
			&skroot_mu_.muboy_ntrack, 
			skroot_mu_.muboy_entpos, 
			skroot_mu_.muboy_dir,
			&skroot_mu_.muboy_goodness, 
			&skroot_mu_.muboy_length, 
			skroot_mu_.muboy_dedx,
			skroot_mu_.mubff_entpos,
			skroot_mu_.mubff_dir, 
			&skroot_mu_.mubff_goodness, 
			&skroot_mu_.muninfo, 
			skroot_mu_.muinfo);
			//delete_outside_hits_();
			skroot_set_tree_(&lun);
			WriteInfo(muonsToRec[i]);
			skroot_fill_tree_(&lun);
	}
	

	
	//return the previous entry so that no issues are caused with other tools
	if(currentEntry != myTreeReader->GetEntryNumber()){
		m_data->getTreeEntry(treeReaderName, currentEntry);
	}
	
	m_data->muonsToRec.clear();
return true;
}

bool RelicMuonMatching::WriteInfo(ParticleCand Event){
	MatchedEvNums.clear();
	PID = 0;
	MatchedTimeDiff.clear();
	
	MatchedEvNums = Event.matchedParticleEvNum;
	
	//MatchedEvNumsBranch->Fill();
	
	MatchedTimeDiff = Event.matchedParticleTimeDiff;
	//MatchedTimeDiffBranch->Fill();
	
	PID = Event.PID;
	
	//PIDBranch->Fill();
return true;
}

bool RelicMuonMatching::Makededx(float (&muentry)[4], float (&mdir)[3], int (&ihcab)[11146], float (&qisk)[11146], float (&tisk)[11146], float (&xyzpm)[11146][3], int& nqisk, float (&dedx)[200]){
	float length;
	int mujresult;
	float kpmt[4], kentry[4], kdir[3], cosang[2];
	double kdist[4];
	int ti, stop_bin, index;
	float sum, save, temp, percent, zerobins;
	int lbins, i, j;
	
	for(int i = 0; i <200; i++){
		dedx[i] = 0;
	}
	
	for(int i = 0; i < 4; i++){
		kentry[i] = muentry[i];
	}
	
	for(int i = 0; i < 3; i++){
		kdir[i] = mdir[i];
	}
	
	for(int i = 0; i < nqisk; i++){
		ti = ihcab[i];
		
		if(qisk[i] > 10){
			kdist[0] = 0;
			kdist[1] = 0;
			
			kpmt[0] = xyzpm[0][ti];
			kpmt[1] = xyzpm[1][ti];
			kpmt[2] = xyzpm[2][ti];
			kpmt[3] = tisk[ti];
			
			mujresult = mujecttime(kentry, kdir, kpmt, kdist, cosang);
			for(int j = 0; j < mujresult; j++){
				stop_bin = (kdist[j]/50.)+1.;
				if((stop_bin > 0) && (stop_bin < 199)){
					dedx[stop_bin] = dedx[stop_bin] + qisk[ti];
				}
			}
		}
	}
	save = 0;
	
	for(int i = 0; i < 190; i++){
		sum = 0;
		for(int j = 0; j < 9; j++){
			sum = sum + dedx[i+j-1];
		}
		if(sum > save){
			index = i;
			save = sum;
		}
	}
	lbins = (length/50.) + 1;
	j = lbins;
	zerobins = 0;
	for(int i = 0; i < j; i++){
		if(dedx[i] == 0){
			zerobins++;
		}
	}
	
	percent = zerobins/lbins;
	temp = 50. * (index + 4);
	
return true;
}


int RelicMuonMatching::mujecttime(float (&v)[4], float (&d)[3], float (&p)[4], double (&dist)[4], float (&cosang)[2]){
	float dx=p[0]-v[0],dy=p[1]-v[1],dz=p[2]-v[2];
	float sprod=dx*d[0]+dy*d[1]+dz*d[2];
	float dr2=dx*dx+dy*dy+dz*dz,dt=p[3]-v[3];
	
	float CVAC = 29.97926;
	float CMED = 21.58333;
	float INDEX = 1.33;
	float CVAC2 = CVAC*CVAC;
	float CMED2 = CMED*CMED;

	dist[0]=CVAC*(sprod*CVAC-dt*CMED2)/(CVAC2-CMED2);
	
	double rad = CVAC2*(CMED2*pow(dt, 2)-dr2)/(CVAC2-CMED2) + pow(dist[0], 2);
	
	if ((dt<0) || (rad<0)){
		dist[0]=dist[1]=-1e10;
		cosang[0]=cosang[1]=-2;
		return(0);
	}
	
	if (rad==0){
		dist[1]=dist[0];
		cosang[0]=cosang[1]=(sprod-dist[0])/sqrt(dr2-2*dist[0]*sprod+dist[0]*dist[0]);
		return(1);
	}
	
	rad=pow(rad, 0.5);
	
	dist[1]=dist[0]+rad;
	dist[0]-=rad;
	
	cosang[0]=(sprod-dist[0])/sqrt(dr2-2*dist[0]*sprod+dist[0]*dist[0]);
	cosang[1]=(sprod-dist[1])/sqrt(dr2-2*dist[1]*sprod+dist[1]*dist[1]);
	dt*=CVAC;
	if (dist[0]>dt) return(0);
	if (dist[1]>dt) return(1);
	return(2);
}

/*int RelicMuonMatching::mujectangle_(float *v,float *d,float *p,float *dist,float *time)
	{
	float dx=p[0]-v[0],dy=p[1]-v[1],dz=p[2]-v[2];
	float sprod=dx*d[0]+dy*d[1]+dz*d[2];
	float dr2=dx*dx+dy*dy+dz*dz,dt=p[3]-v[3],rad=(dr2-sprod*sprod)/(INDEX*INDEX-1);

	if (rad<0)
	{
		*dist=-1e10;
		*time=v[3]-1e10;
		return(0);
	}
	if (rad>0) rad=sqrt(rad);
	*dist=sprod-rad;
	*time=v[3]+(sprod-rad)/CVAC+rad*INDEX/CMED;
return(1);  
}*/







// DON'T LOOK AT THE CODE BELOW HERE SHHHHHHHHHHHHHHH

/*float RelicMuonMatching::rollOver(unsigned long long int currentTime, unsigned long long int targetTime){
	unsigned long long int bitOne = 1;
	unsigned long long int tDiff;
	tDiff = currentTime - targetTime;
	tDiff = (bitOne << 47) + tDiff;
	
}*/

/*bool RelicMuonMatching::AddParticletoDeque(std::deque<ParticleCand>& addToThisDeque){
	unsigned long long int newTime = bitshiftTime(skheadqb_.it0xsk, skheadqb_.nevhwsk);
	ParticleCand newParticle;
	newParticle.EventNumber = myHeader->nevsk;
	newParticle.EventTime = newTime;
	newParticle.EntryNumber = myTreeReader->GetEntryNumber();
	if(particleType == "LOWE"){
		newParticle.ReconEnergy = skroot_lowe_.bsenergy;
	} else{
		newParticle.ReconEnergy = 0.0;
	}
	addToThisDeque.push_back(newParticle);
	
	return true;
	}*/
	
/*unsigned long long int RelicMuonMatching::bitshiftTime(unsigned long long int t0Time, unsigned long long int hardwareTime){
	
	unsigned long long int shiftedt0Time, shiftedhardwareTime, oneint;
	
	// equivalent to 00000000000000001111111111111111 in binary
	oneint = 65535;
	
	shiftedt0Time = t0Time >> 16;
	shiftedt0Time = shiftedt0Time << 16;
	
	shiftedt0Time = shiftedt0Time | (t0Time & oneint);
	
	shiftedhardwareTime = hardwareTime >> 17;
	shiftedhardwareTime = shiftedhardwareTime << 32;
	
	shiftedt0Time = shiftedt0Time + shiftedhardwareTime;
	
	t0Time = shiftedt0Time;
	
	
	return t0Time;
}*/