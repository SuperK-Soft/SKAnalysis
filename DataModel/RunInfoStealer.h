class RunInfoStealer : public TreeManager {
	public:
	RunInfoStealer(){};
	
	bool StealRunInfo(std::string file){
		// get reference tree
		TFile fref(file.c_str(),"READ");
		TTree* tref=(TTree*)fref.Get("data");
		
		// copy over user info.
		// technically this isn't needed for this hack any more.
		for(int i=0; i<tref->GetUserInfo()->GetEntries(); ++i){
			theTree->GetUserInfo()->Add(tref->GetUserInfo()->At(i)->Clone());
		}
		
		// if we need to call skruninf_ now, it will be too late, as this just queries
		// the RINFO member of the TreeManager, which is set on Initialisation.
		// So what we need to do now is set the RunINFO member of the TreeManager
		// get the RunInfo from reference file
		TList* rare = (TList*)tref->GetUserInfo()->FindObject("Rare");
		RunInfo* ri_input = (RunInfo*)rare->At(0);
		
		// copy to current tree
		RINFO = (RunInfo*)ri_input->Clone();
		return (RINFO!=nullptr);
		
	}
};

/* usage demo:
	std::string readerName;
	m_variables.Get("reader",readerName);
	int LUN = m_data->GetLUN(readerName);
	TreeManager* mgr = skroot_get_mgr(&LUN);
	RunInfoStealer* h=(RunInfoStealer*)mgr;
	h->StealRunInfo("/disk2/disk02/data7/sk6/run/0866/086614/rfm_run086614.000001.root");
	
	// this will now get the run info from the reference file above
	runinfsk_();
*/
