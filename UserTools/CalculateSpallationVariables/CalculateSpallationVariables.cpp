#include "CalculateSpallationVariables.h"

#include "MTreeReader.h"

#include "TH1D.h"
#include "geotnkC.h"  // for SK tank geometric constants

CalculateSpallationVariables::CalculateSpallationVariables():Tool(){}


bool CalculateSpallationVariables::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
  
  GetReaders();

  m_variables.Get("run_type", run_type_str);
  if (run_type_str != "cut" && run_type_str != "calculate"){
    throw std::runtime_error("CalculateSpallationVariables::Initialise - no valid run type (calculate / cut) specified in the config file!");
  }
  return true;
}


bool CalculateSpallationVariables::Execute(){

  GetMuonBranchValues();

  // do I need to return for badly reconstructed muons? yeah, probably - only when we use `all'
  if (MU_ptr->muboy_status == 0){
    return true;
  }
  
  // "pre" means dt > 0, "post" is opposite - remember these are the times from the muon to the relic, so it's the opposite to the white paper.
  
  std::cout << "CalculateSpallationVariables: looping through " << MatchedTimeDiff_ptr->size() << " matched relics" << std::endl;
  for (int relic_idx = 0; relic_idx < MatchedTimeDiff_ptr->size(); ++relic_idx){
    
    /* dt - time difference between muon and relic candidate */
    /* mu - relic */
    float dt = MatchedTimeDiff_ptr->at(relic_idx)/pow(10,9);
    dt > 0 ? pre_dt_hist.Fill(abs(dt)) : post_dt_hist.Fill(abs(dt));
    
    relic_tree_ptr->GetEntry(MatchedOutEntryNums_ptr->at(relic_idx));
    GetRelicBranchValues();

    if (LOWE_ptr->bsenergy > 1000){
      //bad reconstruction, skipping
      continue;
    }

    /* dlt - transverse distance between muon and relic candidate */
    float dlt = 0, appr = 0;

    float* muon_entrypoint = nullptr;
    float* muon_direction = nullptr;
    
    bool did_bff = MU_ptr->muinfo[6];
    double muon_tracklen = 0;
    
    if (did_bff){
      basic_array<float> bff_entrypoint(MU_ptr->mubff_entpos);
      basic_array<float> bff_dir = MU_ptr->mubff_dir;

      muon_entrypoint = const_cast<float*>(bff_entrypoint.data());
      muon_direction = const_cast<float*>(bff_dir.data());

      muon_tracklen = CalculateTrackLen(muon_entrypoint, muon_direction);
    } else {
      basic_array<float[10][4]> muboy_entrypoint(MU_ptr->muboy_entpos);
      basic_array<float> muboy_dir = MU_ptr->muboy_dir;

      int muboy_idx = MU_ptr->muinfo[7];
      muon_entrypoint = const_cast<float*>(muboy_entrypoint[muboy_idx].data());
      muon_direction = const_cast<float*>(muboy_dir.data());

      muon_tracklen = muboy_idx == 0 ? MU_ptr->muboy_length : CalculateTrackLen(muon_entrypoint, muon_direction);
      
    }
    
    float* relic_pos = const_cast<float*>(LOWE_ptr->bsvertex);
    
    getdl_(muon_direction,
	   &relic_pos[0],
	   &relic_pos[1],
	   &relic_pos[2],
	   muon_entrypoint,
	   &dlt,
	   &appr);

    if (dlt == 0.0){
      std::cout << "CalculateSpallationVariables::Execute - dlt = 0, dumping args of getdl_" << std::endl;
      std::cout<<"calling getdl_ with:\n"
	       <<"\trelic po: ("<<relic_pos[0]<<", "<<relic_pos[1]<<", "<<relic_pos[2]<<")\n"
	       <<"\tmuon entry point: ("<<muon_entrypoint[0]<<", "<<muon_entrypoint[1]<<", "<<muon_entrypoint[2]<<")\n"
	       <<"\tmuon entry dir: ("<<muon_direction[0]<<", "<<muon_direction[1]<<", "<<muon_direction[2]<<")"<<std::endl;
      throw std::runtime_error("CalculateSpallationVariables:: bad dlt");
    }
    
    dt > 0 ? pre_dlt_hist.Fill(dlt) : post_dlt_hist.Fill(dlt);
    
    /* dll - longitudinal distance between the muon and relic candidate*/
    basic_array<float> scott_dedx(MU_ptr->muboy_dedx);
    float* muon_dedx = const_cast<float*>(scott_dedx.data());
    
    double max_edep = 0;
    int max_edep_bin=0;
    for(int i=0;i<111;i++){
      double e_dep_in_window = 0 ;
      for(int j=0;j<9;j++){
	e_dep_in_window = e_dep_in_window + muon_dedx[i+j];
      }
      if(e_dep_in_window > max_edep){
	max_edep_bin = i+4;
	max_edep = e_dep_in_window;
      }
    }
    double max_energy_dep_pos = 50.*max_edep_bin;
    float dll = max_energy_dep_pos - appr;

    if (dll == 0.0){
      std::cout << "CalculateSpallationVariables::Execute - dll = 0, dumping args of getdl_" << std::endl;
      std::cout<<"calling getdl_ with:\n"
	       <<"\trelic po: ("<<relic_pos[0]<<", "<<relic_pos[1]<<", "<<relic_pos[2]<<")\n"
	       <<"\tmuon entry point: ("<<muon_entrypoint[0]<<", "<<muon_entrypoint[1]<<", "<<muon_entrypoint[2]<<")\n"
	       <<"\tmuon entry dir: ("<<muon_direction[0]<<", "<<muon_direction[1]<<", "<<muon_direction[2]<<")"<<")\n"
	       <<"\tmax_energy_dp_pos: "<< max_energy_dep_pos<< ", appr: "<<appr<<std::endl;
      throw std::runtime_error("CalculateSpallationVariables:: bad dll");
    }


    dt > 0 ? pre_dll_hist.Fill(dll) : post_dll_hist.Fill(dll);

    /* muqismsk - max charge deposited in the detector by the muon*/
    float muqismsk = (MU_ptr->muqismsk);
    dt > 0 ? pre_muqismsk_hist.Fill(muqismsk) : post_muqismsk_hist.Fill(muqismsk);
    
    /* 
       resQ - residual charge deposited by the muon compared to the value expected from the min ionization. 
       reQ = muqismsk - q_MI * L, where q_MI is the number of photoelectrons per cm expected from the min ionization
       and L is the track length. L = sum_{i}(Li) for multiple tracks. Does that not make sense to you? yeah well get in line
    */

    double pe_per_coulomb = 30;     // this might come back to bite me but break in the sun till the sun breaks down old boy
    const float pe_per_cm = 26.78;
    double pe_from_muon = MU_ptr->muqismsk * (pe_per_cm / pe_per_coulomb);  // pe*cm^-1 / pe*C^-1 = C/cm
    double pe_from_MIP = muon_tracklen * pe_per_cm;
    float resQ = pe_from_muon - pe_from_MIP;
    dt > 0 ? pre_resQ_hist.Fill(resQ) : post_resQ_hist.Fill(resQ);

    /* lastly, we need the type of muon event: misfit=0, single_through=1, single_stopping=2, multi=3,4, corner=5 */
    int muon_type = MU_ptr->muboy_status;

    /* oh we also need to the bsenergy: */
    float bse = LOWE_ptr->bsenergy;
    
    PairingInfo p = {dt, dlt, dll, muqismsk, resQ, muon_type, bse};
    std::string p_str = GetPairingString(p);
    pairings[p_str].push_back(p);
  }
  
  return true;
}


bool CalculateSpallationVariables::Finalise(){
  
  std::string outputfile_str = "";
  m_variables.Get("outputfile_str", outputfile_str);
  if (outputfile_str.empty()){throw std::runtime_error("CalculateSpallationVariables::Finalise(): no output file specified!");}

  TFile output_file = TFile(outputfile_str.c_str(), "RECREATE");
  if (output_file.IsZombie()){throw std::runtime_error("CalculateSpallationVariables::Finalise(): couldn't open output file!");}
  output_file.cd();

  pre_dt_hist.Write();
  post_dt_hist.Write();
  pre_dlt_hist.Write();
  post_dlt_hist.Write();
  pre_dll_hist.Write();
  post_dll_hist.Write();
  pre_muqismsk_hist.Write();
  post_muqismsk_hist.Write();
  pre_resQ_hist.Write();
  post_resQ_hist.Write();

  std::cout << "number of hists " << pairings.size() << std::endl;
  
  for (const auto& [name, v_pairing] : pairings){
    std::cout << "name: " << name << std::endl;
    CreateLikelihood(name, v_pairing);
    
  }
  
  output_file.Write();
  
  return true;
}

void CalculateSpallationVariables::GetReaders(){
  std::string relic_reader_name = "";
  if (!m_variables.Get("relic_reader_name", relic_reader_name) || relic_reader_name.empty()){
    throw std::runtime_error("CalculateSpallationVariables::GetReader - no relic_reader_name specified!");
  }
  if (m_data->Trees.count(relic_reader_name) == 0){
    throw std::runtime_error("CalculateSpallationVariables::GetReader - relic reader not found!");
  }
  relic_tree_ptr = m_data->Trees.at(relic_reader_name);

  std::string muon_reader_name = "";
  if (!m_variables.Get("muon_reader_name", muon_reader_name) || muon_reader_name.empty()){
    throw std::runtime_error("CalculateSpallationVariables::GetReader - no muon_reader_name specified!");
  }
  if (m_data->Trees.count(muon_reader_name) == 0){
    throw std::runtime_error("CalculateSpallationVariables::GetReader - muon reader not found!");
  }
  muon_tree_ptr = m_data->Trees.at(muon_reader_name);

}

void CalculateSpallationVariables::GetMuonBranchValues(){
  bool ok = muon_tree_ptr->Get("MatchedTimeDiff", MatchedTimeDiff_ptr);
  if (!ok){throw std::runtime_error("CalculateSpallationVariables::GetMuonBranchValues: Couldn't retrieve MatchedTimeDiff!");}
  ok = muon_tree_ptr->Get("MatchedOutEntryNums", MatchedOutEntryNums_ptr);
  if (!ok){throw std::runtime_error("CalculateSpallationVariables::GetMuonBranchValues: Couldn't retrieve MatchedOutEntryNums!");}
  ok = muon_tree_ptr->Get("MU", MU_ptr);
  if (!ok){throw std::runtime_error("CalculateSpallationVariables::GetMuonBranchValues: Couldn't retrieve MU!");}
  return;
}

void CalculateSpallationVariables::GetRelicBranchValues(){
  LOWE_ptr = nullptr;
  bool ok = relic_tree_ptr->Get("LOWE", LOWE_ptr);
  return;
}

std::string CalculateSpallationVariables::GetPairingString(PairingInfo p) const {
  std::string hist_str = "";

  return "all";
  
  const std::vector<std::string> type_strs = {"misfit", "singlethru", "singlestop", "multi", "multi", "corner"};
  hist_str+=type_strs.at(p.muon_type);

  //  return hist_str; //debug for now
  
  hist_str+="_";
  float bse = p.bse;
  if (bse > 8 && bse < 10){hist_str+="bse8-10";}
  else if (bse > 10 && bse < 12){hist_str+="bse10-12";}
  else if (bse > 12 && bse < 14){hist_str+="bse12-14";}
  else if (bse > 14 && bse < 16){hist_str+="bse14-16";}
  else if (bse > 16 && bse < 18){hist_str+="bse16-18";}
  else if (bse > 18 && bse < 20){hist_str+="bse18-20";}
  else if (bse > 20 && bse < 24){hist_str+="bse20-24";}

  hist_str+="_";
  const float dt = abs(p.dt);
  if (dt > 0 && dt <= 0.05){hist_str+="dt-short";
  } else if (dt > 0.05 && dt <= 0.5){hist_str+="dt-med";
  } else if (dt > 0.5/* && dt < 30*/){hist_str+="dt-long";}

  if (p.muon_type == 0){return hist_str;}

  const float dlt = p.dlt;
  if (dlt > 0 && dlt <= 300){hist_str+="_dlt-short";}
  else if (dlt > 300 && dlt <=  1000){hist_str+="_dlt-med";}
  else if (dlt > 1000/* && dlt < 1000*/){hist_str+="_dlt-long";}
  return hist_str;
}

void CalculateSpallationVariables::CreateLikelihood(const std::string& name, const std::vector<PairingInfo>& vp) const {
  TH1D pre_dt("pre_dt", "pre_dt", nbins, 0, 60), post_dt("post_dt", "post_dt", nbins, 0, 60),
    pre_dlt("pre_dlt", "pre_dlt", nbins, 0, 5000), post_dlt("post_dlt", "post_dlt", nbins, 0, 5000),
    pre_dll("pre_dll", "pre_dll", nbins, -5000, 5000), post_dll("post_dll", "post_dll", nbins, -5000, 5000),
    pre_muqismsk("pre_muqismsk", "pre_muqismsk", nbins, 0, 250000), post_muqismsk("post_muqismsk", "post_muqismsk", nbins, 0, 250000),
    pre_resQ("pre_resQ", "pre_resQ", nbins, -100000, 1000000), post_resQ("post_resQ", "post_resQ", nbins, -100000, 1000000);
  for (const auto& p : vp){
    if (p.dt > 0){
      pre_dt.Fill(abs(p.dt));
      pre_dlt.Fill(p.dlt);
      pre_dll.Fill(p.dll);
      pre_muqismsk.Fill(p.muqismsk);
      pre_resQ.Fill(p.resQ);
    } else {
      post_dt.Fill(abs(p.dt));
      post_dlt.Fill(p.dlt);
      post_dll.Fill(p.dll);
      post_muqismsk.Fill(p.muqismsk);
      post_resQ.Fill(p.resQ);
    }
  }

  TH1D dt_spall(("dt_spall_"+name).c_str(), "dt_spall;dt", nbins, 0, 60);
  TH1D dt_rand(("dt_rand_"+name).c_str(), "dt_rand;dt", nbins, 0, 60);
  for (int i = 0; i <= nbins; ++i){
    dt_spall.SetBinContent(i,  pre_dt.GetBinContent(i) - post_dt.GetBinContent(i));
    dt_rand.SetBinContent(i, post_dt.GetBinContent(i));
  }

  dt_spall.Scale(1/dt_spall.Integral());
  dt_rand.Scale(1/dt_rand.Integral());

  dt_spall.Write();
  dt_rand.Write();

  TH1D dlt_spall(("dlt_spall_"+name).c_str(), "dlt_spall;dlt", nbins, 0, 5000);
  TH1D dlt_rand(("dlt_rand_"+name).c_str(), "dlt_rand;dlt", nbins, 0, 5000);
  for (int i = 0; i <= nbins; ++i){
    dlt_spall.SetBinContent(i, pre_dlt.GetBinContent(i) - post_dlt.GetBinContent(i));
    dlt_rand.SetBinContent(i, post_dlt.GetBinContent(i));
  }

  dlt_spall.Scale(1/dlt_spall.Integral());
  dlt_rand.Scale(1/dlt_rand.Integral());

  dlt_spall.Write();
  dlt_rand.Write();

  TH1D dll_spall(("dll_spall_"+name).c_str(), "dll_spall;dll", nbins, -5000, 5000);
  TH1D dll_rand(("dll_rand_"+name).c_str(), "dll_rand;dll", nbins, -5000, 5000);
  for (int i = 0; i <= nbins; ++i){
    dll_spall.SetBinContent(i, pre_dll.GetBinContent(i) - post_dll.GetBinContent(i));
    dll_rand.SetBinContent(i, post_dll.GetBinContent(i));
  }

  dll_spall.Scale(1/dll_spall.Integral());
  dll_rand.Scale(1/dll_rand.Integral());

  dll_spall.Write();
  dll_rand.Write();

  TH1D muqismsk_spall(("muqismsk_spall_"+name).c_str(), "muqismsk_spall;muqismsk", nbins, 0, 250000);
  TH1D muqismsk_rand(("muqismsk_rand_"+name).c_str(), "muqismsk_rand;muqismsk", nbins, 0, 250000);
  for (int i = 0; i <= nbins; ++i){
    muqismsk_spall.SetBinContent(i, pre_muqismsk.GetBinContent(i) - post_muqismsk.GetBinContent(i));
    muqismsk_rand.SetBinContent(i, post_muqismsk.GetBinContent(i));
  }

  muqismsk_spall.Scale(1/muqismsk_spall.Integral());
  muqismsk_rand.Scale(1/muqismsk_rand.Integral());

  muqismsk_spall.Write();
  muqismsk_rand.Write();

  TH1D resQ_spall(("resQ_spall_"+name).c_str(), "resQ_spall;resQ", nbins, -100000, 100000 );
  TH1D resQ_rand(("resQ_rand_"+name).c_str(), "resQ_rand;resQ", nbins, -100000, 100000);
  for (int i = 0; i <= nbins; ++i){
    resQ_spall.SetBinContent(i, pre_resQ.GetBinContent(i) - post_resQ.GetBinContent(i));
    resQ_rand.SetBinContent(i, post_resQ.GetBinContent(i));
  }

  resQ_spall.Scale(1/resQ_spall.Integral());
  resQ_rand.Scale(1/resQ_rand.Integral());

  resQ_spall.Write();
  resQ_rand.Write();

  TH1D likepre(("likepre_"+name).c_str(), "likepre;L_{spall}", nbins, -30, 30);
  TH1D likepost(("likepost_"+name).c_str(), "likepost;L_{spall}", nbins, -30, 30);

  for (const auto& p : vp){

    const int dt_bin = nbins * ((abs(p.dt) / (60)));
    const int dlt_bin = nbins * (p.dlt / 5000);
    const int dll_bin = nbins * ((p.dll + 5000)/(10000));
    const int muqismsk_bin = nbins * (p.muqismsk / 250000);
    const int resQ_bin = nbins * ((p.resQ + 100000)/(200000));

    //std::cout << dt_bin << dlt_bin << dll_bin << muqismsk_bin << resQ_bin << std::endl;
    
    double likelihood = std::log10(
				   dt_spall.GetBinContent(dt_bin) / dt_rand.GetBinContent(dt_bin) *
				   dlt_spall.GetBinContent(dlt_bin) / dlt_rand.GetBinContent(dlt_bin) *
				   dll_spall.GetBinContent(dll_bin) / dll_rand.GetBinContent(dll_bin) *
				   muqismsk_spall.GetBinContent(muqismsk_bin) / muqismsk_rand.GetBinContent(muqismsk_bin) *
				   resQ_spall.GetBinContent(resQ_bin) / resQ_rand.GetBinContent(resQ_bin));

    p.dt > 0 ? likepre.Fill(likelihood) : likepost.Fill(likelihood);
    
  }

  likepre.Write();
  likepost.Write();
  
}

double CalculateSpallationVariables::CalculateTrackLen(float* muon_entrypoint, float* muon_direction, double* exitpt){

  
  /*
    If you're reading this Marcus, Hello - you might be thinking that Matthew has copied this function from the RelicMuonPlots verbatim and you would indeed be correct. 
    You see, when I originally conceived of this tool I thought it wouldn't overlap quite so heavily with what you had done already - unforunately I was wrong. But nice thing for you, you won't need to check this when I push it because, hey ho, you wrote it didn't you. Feel free to bring up this sloppiness of mine in the form of a slack message but know that the only response you're going to get is: "that's terribly interesting, would you like a brew". 
  */
  
  // for reference HITKTK is the water volume height and DITKTK is its diameter,
  // these are #defined constants in geotnkC.h
	
  // sanity check muon is inward going or it's not going to go through the tank
  if( ((-muon_entrypoint[0]*muon_direction[0] + -muon_entrypoint[1]*muon_direction[1])<0) &&
      (muon_entrypoint[1]>=DITKTK/2.) ){
    // not radially inwards
    Log(m_unique_name+": Muon trajectory is not into tank!",v_error,m_verbose);
    return 0;
  }
  if( (muon_entrypoint[2] >= ( HITKTK/2.) && muon_direction[2]>0) ||
      (muon_entrypoint[2] <= (-HITKTK/2.) && muon_direction[2]<0) ){
    Log(m_unique_name+": Muon track points out of endcaps!",v_error,m_verbose);
    // pointing out of barrel
    return 0;
  }
	
  // calculate track length under assumption of a through-going muon,
  // based on its entry point and direction
  // first check for muons directed in the x-y plane
  if(std::abs(muon_direction[2]) > 0.1){
    double dist_to_endcap = (HITKTK/2.) - std::abs(muon_entrypoint[2]);
    if(std::signbit(muon_entrypoint[2]) != std::signbit(muon_direction[2])){
      dist_to_endcap = (HITKTK - dist_to_endcap);
    }
    // start with the simple case; project to the plane of the appropriate endcap
    double dxdz = muon_direction[0]/std::abs(muon_direction[2]);
    double dydz = muon_direction[1]/std::abs(muon_direction[2]);
    double proj_x = muon_entrypoint[0] + dxdz*dist_to_endcap;
    double proj_y = muon_entrypoint[1] + dydz*dist_to_endcap;
    if(((proj_x*proj_x)+(proj_y*proj_y)) < std::pow(DITKTK/2.,2.)){
      // if the projected point is within the tank radius, this is where it exits the tank
      double xtravel = proj_x-muon_entrypoint[0];
      double ytravel = proj_y-muon_entrypoint[1];
      double tracklen = std::sqrt(std::pow(xtravel,2)+
				  std::pow(ytravel,2.)+
				  std::pow(dist_to_endcap,2.));
			
      if(exitpt){
	exitpt[0] = proj_x;
	exitpt[1] = proj_y;
	exitpt[2] = muon_entrypoint[2] + dist_to_endcap*(muon_direction[2]>0 ? 1 : -1);
      }
      return tracklen;
    }
  }
  // if the muon is in the x-y plane, or the projected point is outside the tank,
  // then the track left through the barrel
	
  // get the angle of the track in the x-y plane
  double trackangle = std::atan2(muon_direction[1],muon_direction[0]);
	
  // get angle the track makes from entry point to the centre of the tank
  double entrypointpolarangle = std::atan2(-muon_entrypoint[1],-muon_entrypoint[0]);
	
  double chordlen;
  if(std::sqrt(std::pow(muon_entrypoint[0],2.)+std::pow(muon_entrypoint[1],2.))==DITKTK/2.){
    // with these we can calculate the angle subtended by the chord the muon track makes with:
    double angle_subtended = M_PI - 2.*(trackangle - entrypointpolarangle);
    // and from this we can work out the length of the chord:
    chordlen = std::abs(DITKTK * std::sin(angle_subtended/2.));
  } else {
    // ah, but if the muon entered through an endcap, the chord is truncated
    double a = std::sqrt(std::pow(muon_entrypoint[0],2.)+std::pow(muon_entrypoint[1],2.));
    // a is the distance from entry point to centre of endcap
    double C = entrypointpolarangle - trackangle;
    // C is the angle from the trajectory to tank centre
    // from a/SinA = c/SinC; -> sinA = a/c*SinC
    double A = std::asin(a/(DITKTK/2.) * std::sin(C));
    // A is the opening angle of the chord
    chordlen = a*std::cos(C) + (DITKTK/2.)*std::cos(A);
  }
	
  // amount of z travel is then chord length times z gradient
  double dzdr = muon_direction[2]/std::sqrt(std::pow(muon_direction[0],2.)+
					    std::pow(muon_direction[1],2.));
  double zdist = dzdr * chordlen;
  // sum is track length
  double tracklen = std::sqrt(std::pow(zdist,2.)+std::pow(chordlen,2.));
	
  if(exitpt){
    exitpt[0] = muon_entrypoint[0] + (chordlen * std::cos(trackangle));
    exitpt[1] = muon_entrypoint[1] + (chordlen * std::sin(trackangle));
    exitpt[2] = muon_entrypoint[2] + zdist;
  }
	
  return tracklen;
	
}
