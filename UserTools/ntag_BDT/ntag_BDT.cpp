#if defined PYTHON
#include "ntag_BDT.h"

ntag_BDT::ntag_BDT():Tool(){
	// get the name of the tool from its class name
	toolName=type_name<decltype(this)>(); toolName.pop_back();
	Py_Initialize(); // XXX THIS MUST BE CALLED IN THE CONSTRUCTOR TO USE PYBIND XXX
}

bool ntag_BDT::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="") m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	m_verbose=1;
	
	m_variables.Get("verbosity",m_verbose);
	
	// get the reader for inputs
	std::string treeReaderName;
	get_ok = m_variables.Get("treeReaderName",treeReaderName);
	if(m_data->Trees.count(treeReaderName)==0){
		Log(toolName+": Failed to find TreeReader "+treeReaderName+" in DataModel!",v_error,m_verbose);
		return false;
	}
	myTreeReader = m_data->Trees.at(treeReaderName);
	
	// intermittently write to disk every N events, so we don't lose everything in case of a crash
	m_variables.Get("writeFrequency",WRITE_FREQUENCY);
	
	// BDT variables
	NLOWINDEX = 5;
	m_variables.Get("NLOWINDEX",NLOWINDEX);
	N10TH = 6;
	m_variables.Get("n10_threshold",N10TH);
	
	// BDT model
	std::string BDT_model= "051_10M.joblib";
	m_variables.Get("BDT_model",BDT_model);
	
	// Import Python modules
	// when do we need to call this?
//	py::scoped_interpreter guard{};  // XXX
	py::object numpy = py::module::import("numpy");
	py::object joblib = py::module::import("joblib");
	py::object jobload = joblib.attr("load");
	
	// load pre-trained BDT
	py::object bdt5 = jobload(BDT_model);
	predict_proba5 = bdt5.attr("predict_proba");
	
	// make output file
	// FIXME move output to datamodel
	std::string outfilepath = "ntag_BDT_out.root";
	m_variables.Get("outfile",outfilepath);
	outfile = new TFile(outfilepath.c_str(), "RECREATE");
	
	// make output tree with same branch structure as input tree
	treeout = myTreeReader->GetTree()->CloneTree(0);
	
	// check if we have the 'type' or 'smearedvertex' branches (old relic analysis MC)
	got_type = myTreeReader->Get("type",type);
	got_smeared_vtx = myTreeReader->Get("smearedvertex", smearedvertex);
	
	// make output branch arrays
	neutron5 = new float[MAX_EVENTS];
	nlow = new int[MAX_EVENTS];
	treeout->Branch("neutron5",  neutron5, "neutron5[np]/F");
	treeout->Branch("nlow",      nlow,     "nlow[np]/I");
	
	// So that downstream Tools can access our results without needing to run a second
	// ToolChain that reads the results off disk, make a TreeReader that they can use to
	// retrieve results as they're generated, equivalently to with an upstream TreeReader Tool
	outTreeReader.Load(treeout);
	m_data->RegisterReader("ntag_BDT_OutTree", &outTreeReader);
	
	return true;
}


bool ntag_BDT::Execute(){
	
	Log(toolName+": Getting branch values",v_debug,m_verbose);
	get_ok = GetBranchValues();
	if(not get_ok){
		Log(toolName+": Error getting branch values!",v_error,m_verbose);
		//return false;
	}
	Log(toolName+": "+toString(np)+" candidates this entry",v_debug,m_verbose);
	
	// unlikely to have >500 neutron candidates in an event,
	// but still better not to segfault if we can avoid it
	if(np > MAX_EVENTS){
		Log(toolName+": expanding output arrays",v_debug,m_verbose);
		delete[] neutron5;
		delete[] nlow;
		MAX_EVENTS = np;
		neutron5 = new float[MAX_EVENTS];
		nlow = new int[MAX_EVENTS];
		treeout->SetBranchAddress("neutron5",  neutron5);
		treeout->SetBranchAddress("nlow",      nlow);
		outTreeReader.UpdateBranchPointer("neutron5");
		outTreeReader.UpdateBranchPointer("nlow");
	}
	
	// Compute position in detector for Nlow calculation (neutron tagging)
	//double rsqred = vx[0]*vx[0] + vy[0]*vy[0];
	//rsqred = rsqred/10000.; // cm^2 to m^2
	//double z = vz[0]/100.;
	////double z = fabs(vz[0])/100.;
	//int id = GetNlowIndex (rsqred, z, NLOWINDEX);
	int id = 0;
	Log("WARNING!!!!!! This is always using Nlow1!!!!",v_debug,m_verbose);
	
	// vector of indices passing preselection
	std::vector<int> passing_indices;
	
	// a vector to store a flattened array of all data; i.e.
	// {{candidate1},{candidate2}} where each {candidate} is an array of BDT input variable values
	std::vector<double> neutronvars;
	
	// number of variables in each candidate array
	int NTAG_VARS;
	
	// loop over neutron capture candidates
	for(int j=0; j<np; j++){
		Log(toolName+": Checking candidate "+toString(j)+", N10="+toString(n10[j]),v_debug,m_verbose);
		
		// initialize output metric for this candidate
		neutron5[j] = -10;
		
		// apply pre-selection cut
		if ( n10[j] < N10TH ) continue;
		
		// it passes preselection
		passing_indices.push_back(j);
		
		// choose which 'Nlow' branch to propagate to the output file
		nlow[j] = Nlow[id].at(j);
		
		// append the set of input variables for this neutron capture candidate
		// to the flattened array to be passed to the BDT
		neutronvars.push_back(n10[j]);        // N10
		// geometrical variables
		neutronvars.push_back(theta[j]);      // θmean                  ?
		neutronvars.push_back(dthetarms[j]);  // θrms?
		neutronvars.push_back(phi[j]);        // ϕrms
		neutronvars.push_back(nlow[j]);       // Nlow
		neutronvars.push_back(nc[j]);         // Ncluster               ?
		neutronvars.push_back(nnlowtheta[j]); // Nlowθ
		neutronvars.push_back(nnback[j]);     // Nback
		// PMT noise variables
		neutronvars.push_back(n300[j]);       // N300
		neutronvars.push_back(nnhighq[j]);    // NhighQ
		neutronvars.push_back(dqmean[j]);     // Qmean
		neutronvars.push_back(dqrms[j]);      // Qrms
		neutronvars.push_back(trmsold[j]);    // Trms
		neutronvars.push_back(mintrms3[j]);   // minTrms(3)
		neutronvars.push_back(mintrms6[j]);   // minTrms(6)
		// Neut-Fit variables
		neutronvars.push_back(fwall[j]);      // NFwall                 ?
		neutronvars.push_back(n10d[j]);       // δN10                   ?
		neutronvars.push_back(trmsdiff[j]);   // δTrms                  ?
		// Bonsai vertex variables
		neutronvars.push_back(bswall[j]);     // BSwall
		neutronvars.push_back(bse[j]);        // BSenergy
		// fit agreement variables
		neutronvars.push_back(bfdist[j]);     // BFdist
		neutronvars.push_back(fpdist[j]);     // FPdist
		
		// if this is the first passing candidate, note the number of variables
		if(j==0) NTAG_VARS = neutronvars.size();
		
	}
	Log(toolName+": "+toString(passing_indices.size())+" candidates passed preselection",v_debug,m_verbose);
	
	// only need to do prediction if any candidates passed preselection
	if(passing_indices.size()){
		
		// when do we need to call this?
//		py::scoped_interpreter guard{};  // XXX ???
		
		// Make a Numpy array with variable info for current event
		// the constructor constructs a 2D Numpy array from a pointer and the dimensions
		const int ncount = neutronvars.size()/NTAG_VARS;
		Log(toolName+": Building pyarray from candidate data",v_debug,m_verbose);
		if(m_verbose>v_debug){
			std::cout<<"ncount = "<<ncount<<", NTAG_VARS="<<NTAG_VARS
				     <<", neutronvars.data()="<<neutronvars.data()
				     <<", neutronvars.size()="<<neutronvars.size()<<std::endl;
		}
		auto arr = py::array_t<double>{{ncount,NTAG_VARS}, neutronvars.data()};
		
		// Make predictions
		Log(toolName+": Calling predict on candidates",v_debug,m_verbose);
		py::array_t<double, py::array::c_style | py::array::forcecast> probas5 = predict_proba5(arr);
		Log(toolName+": Prediction done",v_debug,m_verbose);
		auto probas_u5 = probas5.unchecked<2>();
		
		// map the output metric values back onto the array of all neutron candidates
		// i.e. if neutrons 1,3,5 passed preselection, BDT metrics 0,1,2 should map
		// to output branch array indexes 1,3,5
		for(int k = 0; k < passing_indices.size(); k++){
			neutron5[passing_indices[k]] = probas_u5(k,1);
		}
	}
	
	// fill output tree with BDT metrics
	Log(toolName+": Filling output branches",v_debug,m_verbose);
	treeout->Fill();
	
	unsigned long num_entries = treeout->GetEntries();
	if(num_entries%WRITE_FREQUENCY){
		// write out intermittently for safety
		Log(toolName+": Updating output TTree",v_message,m_verbose);
		outfile->Write("*",TObject::kOverwrite);
	}
	
	return true;
}


bool ntag_BDT::Finalise(){
	
	if(outfile){
		outfile->Write("*",TObject::kOverwrite);
		outfile->Close();
		delete outfile;
		outfile = nullptr;
	}
	
	// inform the TreeReader not to close the file when the destructor is called.
	outTreeReader.SetClosed();
	
	return true;
}

bool ntag_BDT::GetBranchValues(){
	
	get_ok  = (myTreeReader->Get( "HEADER",          HEADER        ));
	get_ok &= (myTreeReader->Get( "LOWE",            LOWE          ));
	get_ok &= (myTreeReader->Get( "np",              np            ));
	get_ok &= (myTreeReader->Get( "nhits",           nnhits        ));
	get_ok &= (myTreeReader->Get( "N10",             n10           ));
	get_ok &= (myTreeReader->Get( "N200M",           N200M         ));
	get_ok &= (myTreeReader->Get( "T200M",           T200M         ));
	get_ok &= (myTreeReader->Get( "Nc",              nc            ));
	get_ok &= (myTreeReader->Get( "Nback",           nnback        ));
	get_ok &= (myTreeReader->Get( "N300",            n300          ));
	get_ok &= (myTreeReader->Get( "NhighQ",          nnhighq       ));
	get_ok &= (myTreeReader->Get( "NLowtheta",       nnlowtheta    ));
	get_ok &= (myTreeReader->Get( "trms",            trmsold       ));
	get_ok &= (myTreeReader->Get( "phirms",          phi           ));
	get_ok &= (myTreeReader->Get( "thetam",          theta         ));
	get_ok &= (myTreeReader->Get( "thetarms",        dthetarms     ));
	get_ok &= (myTreeReader->Get( "Qrms",            dqrms         ));
	get_ok &= (myTreeReader->Get( "Qmean",           dqmean        ));
	get_ok &= (myTreeReader->Get( "trmsdiff",        trmsdiff      ));
	get_ok &= (myTreeReader->Get( "mintrms_6",       mintrms6      ));
	get_ok &= (myTreeReader->Get( "mintrms_3",       mintrms3      ));
	get_ok &= (myTreeReader->Get( "bwall",           bswall        ));
	get_ok &= (myTreeReader->Get( "bse",             bse           ));
	get_ok &= (myTreeReader->Get( "fpdist",          fpdist        ));
	get_ok &= (myTreeReader->Get( "bpdist",          bfdist        ));
	get_ok &= (myTreeReader->Get( "fwall",           fwall         ));
	get_ok &= (myTreeReader->Get( "N10d",            n10d          ));
	get_ok &= (myTreeReader->Get( "dt",              dt            ));
	//get_ok &= (myTreeReader->Get("neutron5",       neutron5      ));
	//get_ok &= (myTreeReader->Get("pvx",            vx            ));
	//get_ok &= (myTreeReader->Get("pvy",            vy            ));
	//get_ok &= (myTreeReader->Get("pvz",            vz            ));
	
	// extra optional branches, exist only for MC w/ old relic analysis
	if(got_smeared_vtx) get_ok &= (myTreeReader->Get( "smearedvertex", smearedvertex ));
	if(got_type) get_ok &= (myTreeReader->Get( "type",type ));
	
	// we have a variable number of 'Nlow' branches.
	// the current apply_ntag.C code defined the number of such branches
	// as a preprocessor macro, but let's make it a bit more flexible.
	// just scan from Nlow1, Nlow2... until we don't find the branch.
	int i=0;
	TTree* t = myTreeReader->GetTree();
	Nlow.clear();
	Log(toolName+": Scanning for Nlow branches",v_debug+1,m_verbose);
	while(true){
		++i;
		std::string nextbranchname = std::string("Nlow")+std::to_string(i);
		// check if branch exists
		if(t->FindBranch(nextbranchname.c_str())!=nullptr){
			// branch exists, add it to the array of pointers
			basic_array<int*> nextnlowearr;
			bool add_ok = (myTreeReader->Get(nextbranchname, nextnlowearr));
			if(not add_ok) break;
			Nlow.push_back(nextnlowearr);
		} else {
			// end of NLow branches
			break;
		}
	}
	// the current code also only ever uses Nlow1....
	get_ok &= (Nlow.size());
	
	return get_ok;
}

// Get suitable threshold for Nlow as a function of the position in the detector
Int_t ntag_BDT::GetNlowIndex(Float_t rsqred, Float_t z, const Int_t init){
	Int_t id = init;
	if ( rsqred >  80 ) id --;
	if ( rsqred > 130 ) id --;
	if ( rsqred > 160 ) id --;
	if ( rsqred > 190 ) id --;
	if ( rsqred > 210 ) id --;
	
	if ( z > 12 ) id --;
	if ( z > 14 ) id --;
	if ( z > 15 ) id --;
	
	if ( id < 0 ) id = 0;
	
	return id;
}

#endif
