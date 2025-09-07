#include "muechk.h"

#include "MTreeReader.h"

#include "fortran_routines.h"
extern "C" {
	void ritofcut_(float*, float*, float*);
}

#include "TH1D.h"

muechk::muechk():Tool(){}

bool muechk::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  m_log= m_data->Log;

  if(!m_variables.Get("verbose",m_verbose)) m_verbose=1;

  std::string reader_name = "";
  m_variables.Get("reader_name", reader_name);
  if (m_data->Trees.count(reader_name) != 1){
    throw std::runtime_error("muechk::Initialise: couldn't get treereader");
  }
  tree_reader_ptr = m_data->Trees.at(reader_name);
  lun = m_data->GetLUN(reader_name);

  nmue_plot = TH1D("nmue_plot", "nmue_plot", 20, 0, 0);
  nmue_times = TH1D("nmue_times", "nmue_times", 100, 0, 0);
  
  return true;
}

bool muechk::Execute(){
  
  /*
  // if we're not calling apfit, hopefully this hack isn't needed.
  // the TreeReader should set run-wise water transparency already anyway...  
  if (skwaterlen_.skwaterlen == 0.0){
    //throw std::runtime_error("skwaterlen_.skwaterlen has value 0! apfit will exit with \"ripecorr.F : APABSPT iregal 0.00000000\" ");
    
    //HACK FOR THE TIME BEING:
    skwaterlen_.skwaterlen = 10000.;
  }
  */
  
  /*
    apfit: apply reconstruction algorithms to one event:
    look in /usr/local/sklib_gcc8/atmpd-trunk/src/recon/aplib/apfit.F for description of input arguments 
    - 0 refers to "apply all reconstruction algorithms for fully contained events (FC)" 

    apfit calls apreset and muechk

    apreset: reset common blocks, calls:
    1) apclrall to clear apmring.h, apfitinf.h, appatsp.h, apmue,h
    2) appclrsep to clear apringsp.h
    3) skread(-1) to clear sktq.h <-- MAKE SURE YOU DON'T NEED THIS

    muechk: counts decay-e event , fills APMUE common
    1st arg: vertex of parent event, can be bsvertex
    2nd arg: verbosity, silent if == 1
    
    yes this tool should be called "apfit" in hindsight, but in hindsight lots of things seem simple don't they... [insert joke about doing a PhD here]
    edit: nope! apfit never worked so we DO now just call muechk directly. All according to plan.

  */

  int dummy_silent = 1;
  //apreset_();
  // apclrsep_();
  //apfit_(&dummy_silent); // << errors with COSANG NaN... 
  
  // ok instead the following process extracted from superscan
  apclrall_();
  int apversion=4;
  aprstbnk_(&apversion);
  
  /* number of sub-events (approximately decay electrons) */
  // (so it says, but always seems to be 0...)
  //kzmnum_(&nsubevts,&isubevt);
  //printf("nsubevts: %d, isubevt: %d\n",nsubevts,isubevt);
  
  /*
  float twin = (skq_.qismsk*0.65 + 1400.)/6.4/40. + 5.;
  if (twin < 30.) twin = 30.;
  float t0;
  float the_vertex[4];
  ritofcut_(the_vertex, &twin, &t0);
  the_vertex[3] = t0;
  bool read_rawdata=false;
  if (read_rawdata) {
  	int minus_lun = -std::abs(lun);
  	int ierr=0;
  	skcrawread_(&minus_lun, &ierr); // reset tube data
  } else {
  	int minus_one=-1;
  	skread_(&minus_one); // reset tube data
  }
  
  // primary vertex will be atmpd fit
  the_vertex[0] = apcommul_.appos[0];
  the_vertex[1] = apcommul_.appos[1];
  the_vertex[2] = apcommul_.appos[2];
  basic_array<float> apv(the_vertex);
  std::cout << "muechk::Execute: apfit vertex: (" << apv.at(0)<<","<<apv.at(1)<<","<<apv.at(2)<<")"<<std::endl;
  
  // these seem to invariably return t0=0 and appos (0,0,0)...
  */
  
  /*
  // even though we're using muechk, the point of this check is to exclude very low energy muons
  // that look like relic/IBD candidates - so they're not going to have been flagged as muons,
  // and won't be reconstructed with muboy.
  MuInfo* mu_ptr = nullptr;
  bool ok = tree_reader_ptr->Get("MU",mu_ptr);
  if(mu_ptr->muboy_ntrack==0){
  	std::cout<<"muechk::Execute: no muboy muons, skipping muechk"<<std::endl;
  	return true;
  }
  muboy_class muboy_status(muboy_class{mu_ptr->muboy_status});
  if(muboy_status!=muboy_class::single_stopping){ // at least exclude 0 (misfit) and 5 (corner clip)
  	std::cout<<"muechk::Execute: muboy_status="<<muboy_status<<" not single stopping, skipping muechk"<<std::endl;
  	return true;
  }
  if(mu_ptr->muboy_goodness < 0.4){
  	std::cout<<"muechk::Execute: muboy_goodness="<<mu_ptr->muboy_goodness<<" poor fit, skipping muechk"<<std::endl;
  	return true;
  }
  basic_array<float> muon_entrypoint(mu_ptr->muboy_entrypoint); // just take first muon?
  */
  
  LoweInfo* lowe_ptr = nullptr;
  bool ok = tree_reader_ptr->Get("LOWE", lowe_ptr);
  if (!ok || lowe_ptr == nullptr){
      throw std::runtime_error("couldn't get lowe branch");
  }
  float bse = lowe_ptr->bsenergy;
  basic_array<float> bsv(lowe_ptr->bsvertex);

  // or if running lfallfit upstream on an RFM file, get from common blocks
  //float bse = skroot_lowe_.bsenergy;
  //basic_array<float> bsv(skroot_lowe_.bsvertex);
  
  Log(m_unique_name+"::Execute: bsenergy: "+std::to_string(bse),v_debug,m_verbose);
  Log(m_unique_name+"::Execute: bsvertex: ("+std::to_string(bsv.at(0))+","+std::to_string(bsv.at(1))
      +","+std::to_string(bsv.at(2))+")"+", time: "+std::to_string(bsv.at(3)),v_debug,m_verbose);
  if(bse==0 || bse>9000){
  	std::cout<<"no or bad lowe reconstruction, skipping muechk"<<std::endl;
  	return true;
  } else {
  	std::cout<<"bsenergy "<<bse<<" is not 0"<<std::endl;
  }
  
  /* think this might have just been part of superscan, not reconstruction.
  // timingGateStart is a GUI control, not set from apfit as far as i can tell.
  if (skheadg_.sk_geometry >= 4) {
  	//set_timing_gate_nsec_(&timingGateStart);
  }
  */
  
  //muechk_(the_vertex, &dummy_silent);
  //muechk_(mu_ptr->muboy_entrypoint, &dummy_silent);
  muechk_(lowe_ptr->bsvertex, &dummy_silent);
  
  // only for zbs?
  //int iatmpd = -1; 
  //apgetbnk_( iatmpd ); 

  int nmue = apmue_.apnmue;
  nmue_plot.Fill(nmue);
  m_data->CStore.Set("nmue", nmue);
  std::vector<double> mue_times(apmue_.apmuetime,apmue_.apmuetime+nmue);
  m_data->CStore.Set("mue_times", mue_times);
  
  Log(m_unique_name+"::Execute: nmue = "+std::to_string(nmue),v_debug,m_verbose);
  for (int i = 0; i < nmue; ++i){
    Log(m_unique_name+"::Execute: found a decay electron with time: "+std::to_string(apmue_.apmuetime[i]),v_debug,m_verbose);
    nmue_times.Fill(apmue_.apmuetime[i]);
  }
  
  
  /*
  // from relic_sk4_ana/data_reduc/third/src/leaf.cc
  int nring = apcommul_.apnring;  
  cout << "nring: " << nring << endl;
  if (nring < 2) thirdred->ring_angle = 0;
  else{
      float ringdir[2][3];
      for (int i = 0; i < 2; i++){
            double judge = (double)appatsp_.approb[ i ][ 1 ] 
                         - (double)appatsp_.approb[ i ][ 2 ]; 
            if ( judge <= 0 ) {  // e-like
              ringdir[ i ][ 0 ] = (double)appatsp_.apmsdir[ i ][ 0 ][ 1 ]; 
              ringdir[ i ][ 1 ] = (double)appatsp_.apmsdir[ i ][ 1 ][ 1 ]; 
              ringdir[ i ][ 2 ] = (double)appatsp_.apmsdir[ i ][ 2 ][ 1 ]; 
            } 
            else {  // mu-like
              ringdir[ i ][ 0 ] = (double)appatsp_.apmsdir[ i ][ 0 ][ 2 ]; 
              ringdir[ i ][ 1 ] = (double)appatsp_.apmsdir[ i ][ 1 ][ 2 ]; 
              ringdir[ i ][ 2 ] = (double)appatsp_.apmsdir[ i ][ 2 ][ 2 ]; 
            }  
      }
      float ringnorm1 = 0;
      float ringnorm2 = 0;
      float ringprod = 0;
      for(int k = 0; k < 3; k++){
          ringnorm1 += pow(ringdir[0][k], 2);
          ringnorm2 += pow(ringdir[1][k], 2);
          ringprod += ringdir[0][k] * ringdir[1][k];
      }
      ringprod /= sqrt(ringnorm1 * ringnorm2);
      thirdred->ring_angle = acos(ringprod) * 180/3.1415926;
      cout << "nring " << nring << " with angle " << thirdred->ring_angle << endl;
  }
  */
  
  
  return true;
}

bool muechk::Finalise(){

  nmue_plot.SaveAs("nmue_plot.root");
  nmue_times.SaveAs("mue_times.root");
  
  return true;
}

extern "C" void lfhits_(int* lfflag, float* sg, float* bg, float* tspread,
		       float* signif, int* lfhit, float* tlf, float* qlf,
		       float* xyzlf, float* ttrg, int* n200, int* nbgbef,
		       int* nbgaft, int* nnba);


extern "C" void lfhits(int* lfflag, float* sg, float* bg, float* tspread,
	    float* signif, int* lfhit, float* tlf, float* qlf,
	    float* xyzlf, float* ttrg, int* n200, int* nbgbef,
	    int* nbgaft){

  std::cout << "calling decoy lfhits" << std::endl;
  
  int my_nnba = 0;
  lfhits_(lfflag, sg, bg, tspread,
	 signif, lfhit, tlf, qlf,
	 xyzlf, ttrg, n200, nbgbef,
	 nbgaft, &my_nnba);
  return;
  
}
