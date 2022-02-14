// Show a demo event, using only in-gate ID hits
TChain* ts = new TChain("data");
ts->Add("*_simple.root");
ts->SetMarkerStyle(20);
ts->SetMarkerSize(0.8);
gStyle->SetPalette(1,0);
ts->Draw("hit_z:hit_x:hit_y:hit_q","Entry$==4 && hit_loc<3 && (hit_ingate&1)==1");
