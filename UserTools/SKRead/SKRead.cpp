#include "SKRead.h"

#include <TString.h>

#include <skheadC.h>
#include <fortran_interface.h>

#include "SKLibs.h"
#include "PathGetter.h"

bool SKRead::Initialise(std::string configfile, DataModel &data)
{
	if(configfile!="")  m_variables.Initialise(configfile);
	m_data= &data;
	m_data->tool_configs[name] = &m_variables;
	m_variables.Get("m_verbose",m_verbose);

    // Read options from config file
    TString inFilePath, skOptions;
    int skBadChOption, skGeometry, refRunNo;
    m_variables.Get("input_file_path", inFilePath);
    m_variables.Get("sk_options", skOptions);
    m_variables.Get("sk_bad_channel_option", skBadChOption);
    m_variables.Get("sk_geometry", skGeometry);
    m_variables.Get("reference_run", refRunNo);

    if (!fileExists(inFilePath.Data())) {
        Log("Input file does not exist. Aborting NTag!", pERROR,m_verbose);
        return false;
    }

    inputIsSKROOT = inFilePath.EndsWith(".root");
    readStatus = 0;

    //////////////////////////////////////////
    // Set read options before opening file //
    //////////////////////////////////////////

    skoptn_(const_cast<char*>(skOptions.Data()), skOptions.Length());
    m_data->GeoSet(skGeometry);
    
    // Custom bad-channel masking (M. Harada)
    // 25: mask bad channel
    // 26: read bad channel from input file
    if ( !skOptions.Contains("25") && !skOptions.Contains("26") ) {
        Log(Form("Masking bad channels with reference run %d", refRunNo));
        int refSubRunNo = 0;
        int outputErrorStatus = 0;
        skbadch_(&refRunNo, &refSubRunNo, &outputErrorStatus);
        skbadopt_(&skBadChOption);
    }

    ///////////////
    // Open file //
    ///////////////

    int lun = 10;
    // SKROOT
    if (inputIsSKROOT) {
        skroot_open_read_(&lun);
        skroot_set_input_file_(&lun, inFilePath.Data(), inFilePath.Length());
        skroot_init_(&lun);
    }
    // ZEBRA
    else {
        m_data->KZInit(); // Initialise ZEBRA

        // Set rflist and open file
        int fileIndex = 1; // 1st file in rflist
        int openError;

        set_rflist_(&lun, inFilePath.Data(),
                    "LOCAL", "", "RED", "", "", "recl=5670 status=old", "", "",
                    inFilePath.Length(), 5, 0, 3, 0, 0, 20, 0, 0);
        int ihndl=0;
        skopenf_(&lun, &fileIndex, "Z", &openError, &ihndl);

        if (openError) {
            Log("Error occurred while opening the input ZEBRA file: " + inFilePath, pERROR,m_verbose);
            return false;
        }
    }
    
    // the ExtractFeatures Tool evaluates inputIsMC in Initialise,
    // but the value was previously not set until the first Execute loop!
    // scan through the file until we have a successful read so we may
    // do the check.
    readStatus = readError;
    while(readStatus!=readOK && readStatus!=readEOF){ // keep reading until good entry or EOF
        readStatus = skread_(&lun);
        if(readStatus==readOK){
            if (skhead_.nrunsk == 999999) {
                m_data->vars.Set("inputIsMC",true);
            }
            else {
                m_data->vars.Set("inputIsMC",true);
            }
            break;
        } else if(readStatus==readEOF){
            Log("Error! SKRead hit EOF before finding any good entries!",pERROR,m_verbose);
            return false;
        }
    }
    
    return true;
}

bool SKRead::Execute()
{
    int lun = 10;
    
    // read loop until we either find a valid event or encounter the end of the toolchain.
    while(true){
        
        // skip first call to skread_ since it was already done in Initialise()
        if(exeCounter>0){ readStatus = skread_(&lun); ++exeCounter; }
        else { readStatus = readOK; ++exeCounter; }

        switch (readStatus) {
            case readOK: {
                Log(Form("Read event #%d successfully.", exeCounter));
                // pointless if checks...?
                // SHE
                if (skhead_.idtgsk & 1<<28);
                // AFT
                else if (skhead_.idtgsk & 1<<29);
                // Else
                else;
                break;
            }
            case readError: {
                Log("Read-error occured!", pWARNING,m_verbose);
                // assume this is a recoverable error and continue to next entry
                //m_data->vars.Set("Skip",true);
                break;
            }
            case readEOF: {
                Log("Reached end-of-file.");
                // skip rest of toolchain as we do not have valid data loaded
                m_data->vars.Set("Skip",true);
                // terminate toolchain
                m_data->vars.Set("StopLoop",true);
                break;
            }
        }
        
    }
    // we may choose to return false if there was a readError
    return true;
}

bool SKRead::Finalise()
{
    int lun = 10;

    if (inputIsSKROOT) {
        skroot_close_(&lun);
        skroot_end_();
    }
    else {
        skclosef_(&lun);
    }

    return true;
}
