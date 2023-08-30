#include "ReconstructMatchedMuons.h"
#include "fortran_routines.h"
#include "Constants.h"
#include "Algorithms.h"

#include <algorithm>

ReconstructMatchedMuons::ReconstructMatchedMuons():Tool(){}


bool ReconstructMatchedMuons::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	m_variables.Get("noBFF", noBFF);
	std::string treeReaderName;
	m_variables.Get("treeReaderName",treeReaderName);
	lun = m_data->GetLUN(treeReaderName);
	if(lun<0){
		Log(m_unique_name+" Error! Could not find TreeReader '"+treeReaderName+"' in DataModel"
		    ,v_error,m_verbose);
		return false;
	}
	
	return true;
}


bool ReconstructMatchedMuons::Execute(){
	
	// for relic spallation checks we only bother with BFF for the muon
	// if the corresponding relic bsenergy is > 12... not sure the logic of that,
	// but we may wish to do a similar thing.
	bool tryBFF=false;
	m_data->vars.Get("tryBFF",tryBFF);
	
	// =================== //
	// Muon Reconstruction //
	// =================== //
	
	// store charge ranges before fix_maxqisk
	skroot_mu_.muinfo[0] = skq_.qismsk;
	skroot_mu_.muinfo[2] = skq_.qimxsk;
	
	// supposedly this undoes an upstream charge saturation correction
	// which was required for SKI-III but is no longer applicable for SKIV+
	fix_maxqisk_();
	
	// save updated charge metrics
	skroot_mu_.muqismsk = skq_.qismsk;
	skroot_mu_.muinfo[3] = skq_.qimxsk;  // should we not update the value in skq_.?
	if(skroot_mu_.muninfo < 4) skroot_mu_.muninfo = 4;
	
	int muyn_org, muynf;
	
	//Muon reconstruction developed by Tomoeda and Yamaguchi
	mfmuselect_(&skroot_mu_.muentpoint, &skroot_mu_.mudir, &skroot_mu_.mugoodness, &muyn_org);
	
	//muyn == 1 - good fit
	//muyn == 0 - bad fit
	
	if(muyn_org > 0){
		skroot_mu_.muyn = 1;
	} else if(muyn_org < 0) {
		skroot_mu_.muyn = 0;
	} else {
		Log("Muyn_org returning as == 0. Not supported yet", v_error, m_verbose);
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
	
	//Apply muboy
	int n_left;
	float  muentry[4], muboy_otherentry[36];
	// $ATMPD_ROOT/src/recon/fit/muboy.F
	muboy_zbs_(&skhead_.nevsk,
	           &skroot_mu_.muboy_status,    // stopping, through-going, corner-clipper, etc. 0=fit failed.
	           &muentry,                     // [0-2]: pos of PMT closest to entry point?, [3]: entry time
	           &skroot_mu_.muboy_dir,       // primary direction at entry point, unit normalised
	           &skroot_mu_.muboy_length,    // track length [cm]
	           &skroot_mu_.muboy_goodness,  // 0-1, higher is better*
	           &skroot_mu_.muboy_ntrack,    // num tracks ("can be 1 if multiple muons" 🤦)
	           &muboy_otherentry,            // additional entry points for tracks 2-5
	           &n_left);                    // num hit PMTs left after cluster finding
	
	// get muon track entry position(s)
	for(int track = 0; track < skroot_mu_.muboy_ntrack; track++){
		if(track == 0){
			skroot_mu_.muboy_entpos[0][track] = muentry[0];
			skroot_mu_.muboy_entpos[1][track] = muentry[1];
			skroot_mu_.muboy_entpos[2][track] = muentry[2];
			skroot_mu_.muboy_entpos[3][track] = muentry[3];
		}else{
			skroot_mu_.muboy_entpos[0][track] = muboy_otherentry[4*track - 7];
			skroot_mu_.muboy_entpos[1][track] = muboy_otherentry[4*track - 6];
			skroot_mu_.muboy_entpos[2][track] = muboy_otherentry[4*track - 5];
			skroot_mu_.muboy_entpos[3][track] = muboy_otherentry[4*track - 4];
		}
	}
	
	if(m_verbose > v_debug+2){
		std::cout << "muboy result:\n"
		          <<"\tgoodness: "<< skroot_mu_.muboy_goodness << "\n"
		          <<"\tntracks: " << skroot_mu_.muboy_ntrack   << "\n"
		          << "\tclass: "  << skroot_mu_.muboy_status   << "\n"
		          << "\tdir: ("   << skroot_mu_.muboy_dir[0]   << ", "
		                          << skroot_mu_.muboy_dir[1]   << ", "
		                          << skroot_mu_.muboy_dir[2]   << ")\n"
		          << "\tlength: " << skroot_mu_.muboy_length   << "\n"
		          << "\tfirst track entry point: ("
		                          << skroot_mu_.muboy_entpos[0] << ", "
		                          << skroot_mu_.muboy_entpos[1] << ","
		                          << skroot_mu_.muboy_entpos[2] << ")\n"
		          << std::endl;
	}
	
	// makededx needs primary entry position and direction.
	// copy muboy values to a set of variables which may overridden by BFF
	// (we already have muentry for position)
	float mudir[3];
	for(int i=0; i<3; ++i){
		mudir[i] = skroot_mu_.muboy_dir[i];
	}
	
	// according to Scott's lowe school slides:
	// for single through-going and stopping muons goodness of 0.4-0.6 is a good fit.
	// for large showers goodness of >0.15 are sometimes ok.
	// goodness <0.1 are always bad.
	
	// if muboy failed, try BFF... but as this can take up to 30 minutes per muon,
	// we have two flags: one Tool config which vetos the use of BFF on anything,
	// and one in `m_data->vars` which may depend on the details of the event.
	if(!noBFF && tryBFF && skroot_mu_.muboy_status == 1 && skroot_mu_.muboy_goodness < 0.4){
		
		Log(m_unique_name+": Muboy failed, trying BFF",v_error,m_verbose);
		
		float bffpos[3];
		float hpos[3];
		newmufit_(&bffpos, &hpos, &skroot_mu_.mubff_goodness);
		
		Log(m_unique_name+": Finished BFF",v_error,m_verbose);
		
		// copy out result
		float modd = sqrt( pow((hpos[0]-bffpos[0]),2) + pow((hpos[1]-bffpos[1]),2) + pow((hpos[2]-bffpos[2]),2) );
		for(int j = 0; j < 3; j++){
			skroot_mu_.mubff_entpos[j] = bffpos[j];
			skroot_mu_.mubff_dir[j] = (hpos[j] - bffpos[j])/modd;
		}
		
		// if BFF succeeded, update the primary muon entry position and direction for makededx
		if(skroot_mu_.mubff_goodness > 0.3){
			for(int j = 0; j < 3; j++){
				muentry[j] = skroot_mu_.mubff_entpos[j];
				mudir[j] = skroot_mu_.mubff_dir[j];
			}
		}
	}
	
	// calculate rate of energy loss along track
	/*
	// lomu_gd.F, which makes lomu_gd files
	makededx_(&muentry,
	          &mudir,
	          &skchnl_.ihcab,
	          &skq_.qisk,
	          &skt_.tisk,
	          &geopmt_.xyzpm,
	          &skq_.nqisk,
	          &skroot_mu_.muboy_dedx);
	*/
	
	// apparently this dedx is somehow 'not good'...?
	// the previous SRN analyses recalculate it with two other methods: kirk's dedx and scott's dedx.
	// FIXME understand what the issue is here and if, or why, we need to call all of them...
	
	// Kirk's method.
	// the relic_sk4_ana repo slightly modified kirk's dedx, accepting water transparency
	// and applying it in an expl function. In kirk's original version the water transparency
	// is looked up internally with lfwatert, but is never used....
	// not sure which to use tbh.
	/*
	float watert;
	get_ok = m_data->vars.Get("watert",watert);
	if(!get_ok){
		// the TreeReader should maintain a suitable value in m_data->vars...
		// is the file MC, but no reference run was provided?
		Log(m_unique_name+" Error! Kirk's `makededx` requires water transparency, "
		    "but none found in m_data->vars!",v_error,m_verbose);
		return false;
	}
	*/
	
	// recalculate dedx
	float (*muinfo)[200] = (float(*)[200])(skroot_mu_.muinfo+10);
	makededx_(&muentry,
	          &skroot_mu_.muboy_dir,
	          &skchnl_.ihcab,
	          &skq_.qisk,
	          &skt_.tisk,
	          &geopmt_.xyzpm,
	          &skq_.nqisk,
	          &skhead_.nrunsk,
//	          &watert,
	          muinfo);
	// here we put kirk's dedx array (200 elements) in muinfo[10]:muinfo[210]
	// previous relic sk4 code saved it in muinfo, but that overwrites:
	// * muinfo[0]: qismsk before fix_maxqisk
	// * mufino[1]: subtrigger number
	// * muinfo[2]: original qimxsk before fix_maxqisk
	// * mufino[3]: qimxsk after fix_maxqisk (redundant with tq_.qimxsk?)
	// * muinfo[4]: parent muon event number (parent muon event for a muon event?)
	// * muinfo[5]: subtrigger number in AFT (how does this differ from rmuinfo[1]?)
	// i can't find anything that seems to use the remaining elements, so just start from 10.
	
	// Scott's method.
	// seems like we're just overwriting the original results with this.
	// if we're just going to overwrite the result, why even call the original makededx?
	// still, this is what lomufit_gd does for the "official" lomugd files,
	// so presumably there's some reason. Maybe the original value is a prerequisite?
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
	
	// pass reconstructed muon info to output Tree branch variables
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
	
	// ======================= //
	// End Muon Reconstruction //
	// ======================= //
	
	
	
	return true;
}

bool ReconstructMatchedMuons::Finalise(){
	
	return true;
}

