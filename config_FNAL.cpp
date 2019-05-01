#include <iostream>
#include <fstream>
#include <array> 
#include "config.h"

#include "TChain.h"
#include "TString.h"
#include "TLatex.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"


namespace po = boost::program_options;
using std::vector; 

//Simple container for related histograms, etc 
struct analysis_o {

  TLatex *latex;
  TH1F *h_amp;
  TH1F *h_amp_cut;
  TH1F *h_time; 
  TProfile2D* p2_amp_vs_XY;
  TH1F *h_deltat;
  TH2F *h2_deltat_vs_amp;
  TProfile *p_deltat_vs_amp; 
  TF1 *fitAmpCorr; 
};



int main(int argc, char **argv) {


    po::options_description odesc ("Allowed Options");
    
    odesc.add_options()
      ("help", "help message")
      ("config", po::value<std::string>()->default_value("simple.ini"), "config file to parse");

    po::variables_map vmap;


    
    try { 
      po::store(po::parse_command_line(argc, argv, odesc), vmap);
      po::notify(vmap);
    }
    catch(po::error& e) {
      std::cout << e.what() << std::endl; 
    return 0;
    }
    
    if (vmap.count("help")) {
      std::cout << odesc << std::endl;
      return 0; 
    }


    std::string cfgname = vmap["config"].as<std::string>();

    std::cout << "Reading configuration values from:" << cfgname << std::endl;

    config cfg(cfgname);

    if (!(cfg.success)) {
      std::cout << "Failed to read config file..." << std::endl;
      exit(-1); 
    }


    std::cout << "Successfully read configuration file" << std::endl; 

    TLegend* leg;
    TH1F* hPad;
    //I have some code for reading from the directory and parsing out the root files
    //I'll switch to that in the next version or two, for now I'm using
    //Andrea's method 
    TChain *tree = new TChain("pulse", "pulse");
  
  
    //Make sure we end our path in a /
    if (cfg.path.back() != '/') 
      cfg.path += '/';
      
    for (int i = cfg.start; i < cfg.stop; i++) {
      TString base; 
      base.Form((cfg.path+"RawDataSaver0CMSVMETiming_Run%d_0_Raw.root").c_str(), i);
      std::cout << "Adding " << base << " to tree." << std::endl; 
      tree->Add(base); 
    }

    int nEntries = tree->GetEntries();
    std::cout << ">>> Events read: " << nEntries << std::endl;

    tree->SetBranchStatus("*", 0);
    //I'm not quite ready to switch over to a config file driven branch setup
    //it won't take much more work, but I haven't finished it yet
    //For now I'm using the manual setup

    tree->SetBranchAddress("amp", &cfg.treeData.amp);
    tree->SetBranchAddress("LP2_50", &cfg.treeData.time);
    tree->SetBranchAddress("gaus_mean", &cfg.treeData.gaus_mean);
    tree->SetBranchAddress("x_dut", &cfg.treeData.x_dut);
    tree->SetBranchAddress("y_dut", &cfg.treeData.y_dut);


    // Mostly modified form of Andrea's code follows from this point on 


    //Build all our histograms, TLatexes, etc
    selections &s = cfg.cuts;

    TH2F* h_beamXY = new TH2F("h_beamXY","h_beamXY",100,s.minX,s.maxX,100,s.minY,s.maxY);
    TH1F* h_beamX  = new TH1F("h_beamX", "h_beamX", 100,s.minX,s.maxX);
    TH1F* h_beamY  = new TH1F("h_beamY", "h_beamY", 100,s.minY,s.maxY);

    std::map<std::string, analysis_o> rootstuff;
    
    for (auto &chan : cfg.channels) {
      auto &ch = chan.second; // Get the channelSettings object
      analysis_o ana;
      ana.latex = new TLatex(0.16,0.96,Form("%s",ch.name.c_str()));
      ana.h_amp = new TH1F(Form("h_amp_%s",ch.name.c_str()), "", s.nAmpBins, s.ampMin, s.ampMax);
      
      ana.h_amp_cut = new TH1F(Form("h_amp_cut_%s", ch.name.c_str()), "", s.nAmpBins, s.ampMin, s.ampMax);
      ana.p2_amp_vs_XY = new TProfile2D(Form("p2_amp_vs_XY_%s", ch.name.c_str()), "", 100, s.minX, s.maxX, 100, s.minY, s.maxY);
      ana.h_time = new TH1F(Form("h_time_%s", ch.name.c_str()), "", s.nTimeBins, s.minTime, s.maxTime);

      //Delta T stuff -- unsure if the ampMin/Max/nBins have the same values as their
      //corresponding channels. If so, should replace these with the channel value 
      ana.h_deltat = new TH1F(Form("h_deltat_%s",ch.name.c_str()),"",s.nDeltatBins,
			  s.minDeltat, s.maxDeltat); 
      ana.h2_deltat_vs_amp = new TH2F(Form("h2_deltat_vs_amp_%s", ch.name.c_str()), "",
				  s.nAmpBins, s.ampMin, s.ampMax, s.nDeltatBins,
				  s.minDeltat, s.maxDeltat);

      ana.p_deltat_vs_amp = new TProfile(Form("p_deltat_vs_amp_%s", ch.name.c_str()), "",
					 s.nAmpBins, s.ampMin, s.ampMax); 


      ana.fitAmpCorr = new TF1(Form("fitAmpCorr_%s", ch.name.c_str()), "pol4",0.,1000.); 
      
      rootstuff[ch.name] = ana; 
      
    }

    //Event Loop
    int selectedEntries = 0;
    for (int entry = 0; entry < nEntries; entry++) {
      std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;


      tree->GetEntry(entry);
      auto &td = cfg.treeData; 
      float myX = td.x_dut[cfg.channels["photek"].ampid];
      float myY = td.y_dut[cfg.channels["photek"].ampid];

      if( myX == -999 || myY == -999 ) continue;

      h_beamX  -> Fill( myX );
      h_beamY  -> Fill( myY );
      h_beamXY -> Fill( myX,myY );


      for (auto &chan : cfg.channels)  {
	auto &ch = chan.second;

	rootstuff[ch.name].h_amp->Fill(td.amp[ch.ampid]);
	rootstuff[ch.name].p2_amp_vs_XY->Fill(myX, myY, td.amp[ch.ampid]);
	// cut on BS
	if( myX == -999 || myY == -999 ) continue;
	if( fabs(myX-s.centerX) > s.half_beam_spot_x ||
	    fabs(myY-s.centerY) > s.half_beam_spot_y ) continue;

	if (td.amp[ch.ampid] > ch.ampmin &&
	    td.amp[ch.ampid] < ch.ampmax )
	  {
	    rootstuff[ch.name].h_amp_cut->Fill(td.amp[ch.ampid]);
	    rootstuff[ch.name].h_time->Fill(td.time[ch.timeid]);
	  }

	++selectedEntries; 

      }

    }
    std::cout << "\n>>> 1st loop: selected entries " << selectedEntries << std::endl;

    //Draw Beam Histograms

    TCanvas* c_beamXY = new TCanvas("c_beamXY","c_beamXY",500,500);
    c_beamXY->cd();
    h_beamXY->SetStats(0);
    h_beamXY->SetTitle(";X [mm]; Y[mm];entries");
    h_beamXY->Draw("COLZ");


    TLine * x_BS_min = new TLine(s.centerX-s.half_beam_spot_x, s.minY,
				 s.centerX-s.half_beam_spot_x, s.maxY);

    TLine * x_BS_max = new TLine(s.centerX+s.half_beam_spot_x, s.minY,
				 s.centerX+s.half_beam_spot_x, s.maxY);

    TLine * y_BS_min = new TLine(s.minX, s.centerY-s.half_beam_spot_y,
				 s.maxX, s.centerY-s.half_beam_spot_y);

    TLine * y_BS_max = new TLine(s.minX, s.centerY+s.half_beam_spot_y,
				 s.maxX, s.centerY+s.half_beam_spot_y);


    x_BS_min->SetLineColor(kRed);
    x_BS_min->SetLineWidth(2);
    x_BS_max->SetLineColor(kRed);
    x_BS_max->SetLineWidth(2);
    y_BS_min->SetLineColor(kRed);
    y_BS_min->SetLineWidth(2);
    y_BS_max->SetLineColor(kRed);
    y_BS_max->SetLineWidth(2);

    c_beamXY -> Print(Form("c_beamXY_%s.png",cfg.label.c_str()));
    int nch = cfg.channels.size(); 
    TCanvas* c_amp_vs_XY = new TCanvas("c_amp_vs_XY","c_amp_vs_XY",1000,500*((nch-1))/2);
    c_amp_vs_XY -> Divide(2,(nch-1)/2);

    //Draw for each non-photek channel
    //Note, left / right may get swapped, will fix in future
    int canvas = 1; 
    for (auto &chan : cfg.channels) {
      auto &ch = chan.second;
      if (ch.name == "photek")
	continue;

      c_amp_vs_XY->cd(canvas);
      canvas++; 
      auto *hist = rootstuff[ch.name].p2_amp_vs_XY; 
      hist->Draw("COLZ"); 
      hist->SetStats(0);
      hist->SetTitle(";X [mm];Y [mm];amplitude [mV]");
      gPad->SetLogz();

      x_BS_min->Draw("same");
      x_BS_max->Draw("same");
      y_BS_min->Draw("same");
      y_BS_max->Draw("same");
      rootstuff[ch.name].latex->Draw("same");
    }
    c_amp_vs_XY -> Print(Form("c_amp_vs_XY_%s.png",cfg.label.c_str()));

    //Fitting and Drawing MIP Peak
    TCanvas* c_amp = new TCanvas("c_amp","c_amp",1000,500*(nch+1)/2);    
    c_amp -> Divide(2,(nch+1)/2);
    canvas = 0; 
    for(auto &chan : cfg.channels) {
      auto &ch = chan.second;
      c_amp->cd(canvas+1);
      gPad->SetLogy();
      canvas++;
      
      auto *h_amp = rootstuff[ch.name].h_amp; 
      auto *h_amp_cut = rootstuff[ch.name].h_amp_cut; 
      h_amp-> SetStats(0);
      h_amp-> SetLineColor(kBlack);
      h_amp-> SetLineWidth(2);
      h_amp-> Draw();
      h_amp-> SetTitle(";max. amp [mV];entries");

      h_amp_cut-> SetFillColor(kOrange-9);
      h_amp_cut-> SetLineColor(kBlack);
      h_amp_cut-> SetLineWidth(2);
      h_amp_cut-> Draw("same");
      if (ch.name != "photek") {
	TF1* fitMipPeak = new TF1 ("fitMipPeak","landau",ch.ampmin,ch.ampmax);
	fitMipPeak -> SetParameter(1,h_amp_cut->GetBinCenter(h_amp_cut->GetMaximum()));
	fitMipPeak->SetNpx(10000);
	fitMipPeak->SetLineColor(kRed);
	fitMipPeak->SetLineWidth(2);
	h_amp_cut->Fit(fitMipPeak,"SQR");
	ch.mip_peak = fitMipPeak->GetParameter(1);
	std::cout << "peak[" << ch.name << "] = " << ch.mip_peak << " mV" << std::endl;

	TLine* lowcut = new TLine(std::max(s.MIP_low_cut*ch.mip_peak,ch.ampmin),0.,
				  std::max(s.MIP_low_cut*ch.mip_peak,ch.ampmin),
				  h_amp->GetMaximum());
	
	TLine* higcut = new TLine(std::min(s.MIP_high_cut*ch.mip_peak,ch.ampmax),0.,
				  std::min(s.MIP_high_cut*ch.mip_peak,ch.ampmax),
				  h_amp->GetMaximum()); 
	lowcut->Draw("same");
	higcut->Draw("same");
      }
      rootstuff[ch.name].latex->Draw("same");
    }

    c_amp->Print(Form("c_amp_%s.png",cfg.label.c_str()));


    //Fitting and Drawing the Time Peak

    TF1* fitTimePeak = new TF1("fitTimePeak","gaus",s.low_time_cut,s.high_time_cut);
    TCanvas* c_time = new TCanvas("c_time","c_time",1000,500*(nch+1)/2);
    c_time -> Divide(2,(nch+1)/2);


    canvas = 0; 
    for (auto &chan : cfg.channels)  {
      auto &ch = chan.second;
      c_time -> cd(canvas+1);
      canvas++; 
      gPad -> SetLogy();
      auto *h_time = rootstuff[ch.name].h_time;
      h_time->SetStats(0);
      h_time->SetFillColor(kOrange-9);
      h_time->SetLineColor(kBlack);
      h_time->Draw();
      h_time->SetTitle(";time [ns];entries");
    
      fitTimePeak->SetParameters(1.,20.,5.);
      fitTimePeak->SetNpx(10000);
      fitTimePeak->SetLineColor(kRed);
      fitTimePeak->SetLineWidth(2);
      h_time->Fit(fitTimePeak,"QRS");
      ch.time_peak = fitTimePeak->GetParameter(1);
      ch.time_sigma = fitTimePeak->GetParameter(2);

      std::cout << "time_peak[" << ch.name << "] = " << ch.time_peak << " :: sigma_time_peak = " << ch.time_sigma << std::endl;            

      TLine* lowcut = new TLine(std::max(ch.time_peak-ch.time_sigma*s.sigma_time_cut,
					 s.low_time_cut),0,
				std::max(ch.time_peak-ch.time_sigma*s.sigma_time_cut,
					 s.low_time_cut),h_time->GetMaximum());
      
      TLine* higcut = new TLine(std::min(ch.time_peak+ch.time_sigma*s.sigma_time_cut,
					 s.high_time_cut),0,
				std::min(ch.time_peak+ch.time_sigma*s.sigma_time_cut,
					 s.high_time_cut),h_time->GetMaximum());
    
      lowcut->Draw("same");
      higcut->Draw("same");
      
      rootstuff[ch.name].latex->Draw("same");

    }
    c_time -> Print(Form("c_time_%s.png",cfg.label.c_str()));

    //-----------------------------------------------------------------------------
    //                   (2) second loop events -  to calculate amp-walk correction
    //-----------------------------------------------------------------------------

    selectedEntries = 0;
    for (int entry = 0; entry < nEntries; ++entry) {
      std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      tree->GetEntry(entry);

      auto &td = cfg.treeData;
      auto &phtk = cfg.channels["photek"];
      
      float myX = td.x_dut[phtk.ampid]; 
      float myY = td.y_dut[phtk.ampid]; 
      // cut on BS
      if( myX == -999 || myY == -999 ) continue;
      if( fabs(myX-s.centerX) > s.half_beam_spot_x ||
	  fabs(myY-s.centerY) > s.half_beam_spot_y ) continue;
    
      float time_ref = td.gaus_mean[phtk.timeid]; 
      
      //cut on MCP amp

      if ((td.amp[phtk.ampid] < phtk.ampmin) ||
	  (td.amp[phtk.ampid] > phtk.ampmax))
	continue;

      //cut on MCP time
      if ( time_ref < std::max(phtk.time_peak-phtk.time_sigma*s.sigma_time_cut,
			       s.low_time_cut) ||
	   time_ref > std::min(phtk.time_peak+phtk.time_sigma*s.sigma_time_cut,
			       s.high_time_cut))
	continue;

      for (auto &chan : cfg.channels) {
	auto &ch = chan.second;

	if ((td.amp[ch.ampid] > std::max(ch.mip_peak*s.MIP_low_cut, ch.ampmin)) &&
	    (td.amp[ch.ampid] < std::min(ch.mip_peak*s.MIP_high_cut, ch.ampmax)) &&
	    (td.time[ch.timeid] > std::max(ch.time_peak-ch.time_sigma*s.sigma_time_cut,
					   s.low_time_cut)) &&
	    (td.time[ch.timeid] < std::min(ch.time_peak+ch.time_sigma*s.sigma_time_cut,
					   s.high_time_cut))) { 
	
	  rootstuff[ch.name].h_deltat->Fill(td.time[ch.timeid]-time_ref);
	  rootstuff[ch.name].h2_deltat_vs_amp->Fill(td.amp[ch.ampid],
						    td.time[ch.timeid]-time_ref);
	  rootstuff[ch.name].p_deltat_vs_amp->Fill(td.amp[ch.ampid],
						 td.time[ch.timeid]-time_ref);
	}
	
      }

      ++selectedEntries;
      
    }

    std::cout << "\n>>> 2nd loop: selected entries " << selectedEntries << std::endl;

    //Fitting and drawing time walk corrections 


    TCanvas* c_time_vs_amp = new TCanvas("c_time_vs_amp","c_time_vs_amp",1000.,500.*((nch-1))/2);
    c_time_vs_amp -> Divide(2,(nch-1)/2);

    canvas = 1; 
    for (auto &chan : cfg.channels) {
      auto &ch = chan.second;
      if (ch.name == "photek")
	continue;

      c_time_vs_amp->cd(canvas);
      canvas++; 
      auto *h2_deltat_vs_amp = rootstuff[ch.name].h2_deltat_vs_amp;
      auto *p_deltat_vs_amp = rootstuff[ch.name].p_deltat_vs_amp; 
      auto *h_deltat = rootstuff[ch.name].h_deltat;
      auto *fitAmpCorr = rootstuff[ch.name].fitAmpCorr;
      
      h2_deltat_vs_amp->SetStats(0);
      h2_deltat_vs_amp->GetYaxis()->SetRangeUser(h_deltat->GetMean()-5.*h_deltat->GetRMS(),h_deltat->GetMean()+5.*h_deltat->GetRMS());
      h2_deltat_vs_amp->SetTitle(";max. amplitude [mV];#Deltat [ns]");
      h2_deltat_vs_amp->Draw("COLZ");
      p_deltat_vs_amp->SetMarkerStyle(20);
      p_deltat_vs_amp->SetMarkerSize(0.7);
      p_deltat_vs_amp->SetMarkerColor(kMagenta);
      p_deltat_vs_amp->Draw("same");
      p_deltat_vs_amp->Fit(fitAmpCorr, "QNRS+");

      fitAmpCorr->SetLineColor(kMagenta);
      fitAmpCorr->Draw("same");

      rootstuff[ch.name].latex->Draw("same");

    }

    c_time_vs_amp -> Print(Form("c_time_vs_amp_%s.png",cfg.label.c_str()));
    


    //--------------------------------------------------------------------------
    //                    (3) third loop events --> to apply amp-walk correction
    //--------------------------------------------------------------------------

      // define histograms
    TH1F* h_deltat_avg = new TH1F("h_deltat_avg","",6000, s.minDeltat, s.maxDeltat);
    TH1F* h_deltat_avg_ampCorr = new TH1F("h_deltat_avg_ampCorr","",6000, s.minDeltat, s.maxDeltat);
  
    // TH1F* h_deltat_avg_ampCorr_comb = new TH1F("h_deltat_avg_ampCorr_comb","",6000, minDeltat, maxDeltat);
  
    TH1F* h_deltat_left  = new TH1F("h_deltat_left","",6000, s.minDeltat, s.maxDeltat);
    TH1F* h_deltat_right = new TH1F("h_deltat_right","",6000, s.minDeltat, s.maxDeltat);
    TH1F* h_deltat_diff  = new TH1F("h_deltat_diff","",6000, s.minDeltat, s.maxDeltat);
  
    TH1F* h_deltat_left_ampCorr  = new TH1F("h_deltat_left_ampCorr","",6000, s.minDeltat, s.maxDeltat);
    TH1F* h_deltat_right_ampCorr = new TH1F("h_deltat_right_ampCorr","",6000, s.minDeltat, s.maxDeltat);
    TH1F* h_deltat_diff_ampCorr  = new TH1F("h_deltat_diff_ampCorr","",6000, s.minDeltat, s.maxDeltat);
  
    TProfile* p_left_vs_X  = new TProfile("p_left_vs_X","",400,s.minX,s.maxX);
    TProfile* p_right_vs_X = new TProfile("p_right_vs_X","",400,s.minX,s.maxX);  
    TProfile* p_avg_vs_X   = new TProfile("p_avg_vs_X","",400,s.minX,s.maxX);
    TProfile* p_diff_vs_X  = new TProfile("p_diff_vs_X","",400,s.minX,s.maxX);
  
    TProfile* p_left_vs_diff  = new TProfile("p_left_vs_diff","",6000,-10,10);
    TProfile* p_right_vs_diff = new TProfile("p_right_vs_diff","",6000,-10,10);  
    TProfile* p_avg_vs_diff   = new TProfile("p_avg_vs_diff","",6000,-10,10);

    //May want to wrap these in an object sometime as well
    
    vector<TH1F *> h_deltat_avg_ampCorr_posCutX(s.n_pos_cuts_x);
    vector<TH1F *> h_deltat_left_ampCorr_posCutX(s.n_pos_cuts_x);
    vector<TH1F *> h_deltat_right_ampCorr_posCutX(s.n_pos_cuts_x);
    vector<TH1F *> h_deltat_avg_ampCorr_posCutY(s.n_pos_cuts_y);
    vector<TH1F *> h_deltat_left_ampCorr_posCutY(s.n_pos_cuts_y);
    vector<TH1F *> h_deltat_right_ampCorr_posCutY(s.n_pos_cuts_y);
    vector< vector<TH1F *> > h_deltat_avg_ampCorr_posCutXY(s.n_pos_cuts_x); 
    
    for (int i = 0; i < s.n_pos_cuts_x; i++) {
      h_deltat_avg_ampCorr_posCutX[i]   = new TH1F(Form("h_deltat_ampCorr_posCut_%d", i), "",3000,s.minDeltat, s.maxDeltat);
      h_deltat_left_ampCorr_posCutX[i]  = new TH1F(Form("h_deltat_left_ampCorr_posCut_%d",i), "",3000,s.minDeltat, s.maxDeltat);
      h_deltat_right_ampCorr_posCutX[i]  = new TH1F(Form("h_deltat_right_ampCorr_posCut_%d",i), "",3000,s.minDeltat, s.maxDeltat);

    }


    for(int i = 0; i < s.n_pos_cuts_y; ++i)
      {
	h_deltat_avg_ampCorr_posCutY[i]   = new TH1F(Form("h_deltat_ampCorr_posCutY_%d", i), "",3000,s.minDeltat, s.maxDeltat);
	h_deltat_left_ampCorr_posCutY[i]  = new TH1F(Form("h_deltat_left_ampCorr_posCutY_%d",i), "",3000,s.minDeltat, s.maxDeltat);
	h_deltat_right_ampCorr_posCutY[i]  = new TH1F(Form("h_deltat_right_ampCorr_posCutY_%d",i), "",3000,s.minDeltat, s.maxDeltat);
      }

    for (int x = 0; x < s.n_pos_cuts_x; x++) {
      h_deltat_avg_ampCorr_posCutXY[x] =  std::vector< TH1F *>(s.n_pos_cuts_y); 
      for (int y = 0; y < s.n_pos_cuts_y; y++) {
	h_deltat_avg_ampCorr_posCutXY[x][y]  = new TH1F(Form("h_deltat_ampCorr_posCutXY_%d_%d", x,y), "",3000,s.minDeltat, s.maxDeltat);

      }
    }

    std::vector<TH1F*> h_deltat_avg_ampCorr_BSCut(s.beam_spot_bins);
    std::vector<TH1F*> h_deltat_left_ampCorr_BSCut(s.beam_spot_bins);
    std::vector<TH1F*> h_deltat_right_ampCorr_BSCut(s.beam_spot_bins);

    //BSCuts should get moved to a config driven form
    //something like a list of numbers
    std::vector<float> BScut(s.beam_spot_bins);
    BScut[0] = 10; // in mm around the center
    BScut[1] = 7;
    BScut[2] = 5;
    BScut[3] = 3;
    BScut[4] = 2;
    BScut[5] = 1;
  

  for (int iBSCut = 0; iBSCut < s.beam_spot_bins; iBSCut++)
  {
    h_deltat_avg_ampCorr_BSCut[iBSCut]  = new TH1F(Form("h_deltat_avg_ampCorr_BSCut_%d", iBSCut), "",s.nDeltatBins, s.minDeltat, s.maxDeltat);
    h_deltat_left_ampCorr_BSCut[iBSCut]  = new TH1F(Form("h_deltat_left_ampCorr_BSCut_%d", iBSCut), "",s.nDeltatBins, s.minDeltat, s.maxDeltat);
    h_deltat_right_ampCorr_BSCut[iBSCut]  = new TH1F(Form("h_deltat_right_ampCorr_BSCut_%d", iBSCut), "",s.nDeltatBins, s.minDeltat, s.maxDeltat);
  }



  //read over events
  selectedEntries = 0;

  //Grr...not a fan of this approach, will have to think about how I want to re-org
  //this to make it a bit less brittle
  

  channelSettings *l;
  channelSettings *r;
  
  for (auto &chan : cfg.channels) {
    if (chan.second.name.find("left") != std::string::npos)
      l = &chan.second;
    else if (chan.second.name.find("right") != std::string::npos)
      r = &chan.second;
  }
  channelSettings &left = *l; 
  channelSettings &right = *r; 

  float lowerPosCutX = s.centerX-s.half_beam_spot_x;
  float upperPosCutX = s.centerX+s.half_beam_spot_x;  
  float lowerPosCutY = s.centerY-s.half_beam_spot_y;
  float upperPosCutY = s.centerY+s.half_beam_spot_y;

  
  
  for (int entry = 0; entry < nEntries; ++entry)
    {

      if( entry%100 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;

      auto &td = cfg.treeData;
      auto &phtk = cfg.channels["photek"];
      
      tree->GetEntry(entry);
      
      float myX = td.x_dut[phtk.ampid];
      float myY = td.y_dut[phtk.ampid];
  
      if( myX == -999 || myY == -999 ) continue;
      if( fabs(myX-s.centerX) > s.half_beam_spot_x ||
	  fabs(myY-s.centerY) > s.half_beam_spot_y ) continue;

      float time_ref = td.gaus_mean[phtk.timeid];

      //cut on MCP amp
      if( td.amp[phtk.ampid] < phtk.ampmin ||
	  td.amp[phtk.ampid] > phtk.ampmax )
	continue;
      //cut on MCP time
      if ( time_ref < std::max(phtk.time_peak-phtk.time_sigma*s.sigma_time_cut,s.low_time_cut) ||
	   time_ref > std::min(phtk.time_peak+phtk.time_sigma*s.sigma_time_cut,s.high_time_cut) )
	continue;

      if (td.amp[left.ampid] < std::max(left.mip_peak*s.MIP_low_cut,
					left.ampmin)) continue;
      if (td.amp[left.ampid] > std::min(left.mip_peak*s.MIP_high_cut,
					left.ampmax)) continue;

      if (td.amp[right.ampid] < std::max(right.mip_peak*s.MIP_low_cut,
					 right.ampmin)) continue;
      if (td.amp[right.ampid] > std::min(right.mip_peak*s.MIP_high_cut,
					 right.ampmax)) continue;

      if (td.time[left.timeid] < std::max(left.time_peak-
					  left.time_sigma*s.sigma_time_cut,
					  s.low_time_cut)) continue;

      if (td.time[left.timeid] > std::min(left.time_peak+
					  left.time_sigma*s.sigma_time_cut,
					  s.high_time_cut)) continue;

      if (td.time[right.timeid] < std::max(right.time_peak-
					   right.time_sigma*s.sigma_time_cut,
					   s.low_time_cut)) continue;

      if (td.time[right.timeid] > std::min(right.time_peak+
					   right.time_sigma*s.sigma_time_cut,
					   s.high_time_cut)) continue;
    
      float amp1 = td.amp[left.ampid];
      float amp2 = td.amp[right.ampid];
      float time1 = td.time[left.timeid]; 
      float time2 = td.time[right.timeid];
      auto &fitAmpCorr = rootstuff[left.name].fitAmpCorr;
      auto &h_amp_cut = rootstuff[left.name].h_amp_cut; 
      float time1_ampCorr = time1 - fitAmpCorr->Eval(amp1) + fitAmpCorr->Eval(h_amp_cut->GetMean());
      fitAmpCorr = rootstuff[right.name].fitAmpCorr;
      h_amp_cut = rootstuff[right.name].h_amp_cut; 
      float time2_ampCorr = time2 - fitAmpCorr->Eval(amp2) + fitAmpCorr->Eval(h_amp_cut->GetMean());
    
      float deltat_avg = 0.5*(time1+time2) - time_ref;
      float deltat_avg_ampCorr = 0.5*(time1_ampCorr+time2_ampCorr) - time_ref;
    
      h_deltat_left  -> Fill( time1 - time_ref );
      h_deltat_right -> Fill( time2 - time_ref );
      h_deltat_diff  -> Fill( time2 - time1 );
      h_deltat_avg   -> Fill( deltat_avg );
    
      h_deltat_left_ampCorr  -> Fill( time1_ampCorr - time_ref );
      h_deltat_right_ampCorr -> Fill( time2_ampCorr - time_ref );      
      h_deltat_diff_ampCorr  -> Fill( time2_ampCorr - time1_ampCorr );    
      h_deltat_avg_ampCorr   -> Fill( deltat_avg_ampCorr );

      // filling plots vs position

      p_left_vs_X  -> Fill( myX,time1_ampCorr-time_ref );
      p_right_vs_X -> Fill( myX,time2_ampCorr-time_ref );
      p_avg_vs_X   -> Fill( myX,0.5*(time1_ampCorr+time2_ampCorr)-time_ref );
      p_diff_vs_X  -> Fill( myX,time2_ampCorr-time1_ampCorr );
    
    
      // filling plots vd t_diff
      p_left_vs_diff  -> Fill( time2_ampCorr-time1_ampCorr,time1_ampCorr-time_ref);
      p_right_vs_diff -> Fill( time2_ampCorr-time1_ampCorr,time2_ampCorr-time_ref);
      p_avg_vs_diff   -> Fill( time2_ampCorr-time1_ampCorr,0.5*(time1_ampCorr+time2_ampCorr)-time_ref);





      // X dependency  
      for(int iPosCut = 0; iPosCut < s.n_pos_cuts_x; ++iPosCut)
	{
	  float step = (upperPosCutX-lowerPosCutX) / s.n_pos_cuts_x;
	  if( myX > lowerPosCutX + iPosCut*step &&
	      myX < lowerPosCutX + (iPosCut+1)*step )
	    {
	      h_deltat_avg_ampCorr_posCutX[iPosCut]  -> Fill( deltat_avg_ampCorr );
	      h_deltat_left_ampCorr_posCutX[iPosCut] -> Fill( time1_ampCorr - time_ref );
	      h_deltat_right_ampCorr_posCutX[iPosCut] -> Fill( time2_ampCorr - time_ref );
	    }
	}
    
      // Y dependency
      for(int iPosCut = 0; iPosCut < s.n_pos_cuts_y; ++iPosCut)
	{
	  float step = (upperPosCutY-lowerPosCutY) / s.n_pos_cuts_y;
	  if( myY >lowerPosCutY + iPosCut*step &&
	      myY < lowerPosCutY+(iPosCut+1)*step )
	    {
	      h_deltat_avg_ampCorr_posCutY[iPosCut]  -> Fill( deltat_avg_ampCorr );
	      h_deltat_left_ampCorr_posCutY[iPosCut] -> Fill( time1_ampCorr - time_ref );
	      h_deltat_right_ampCorr_posCutY[iPosCut] -> Fill( time2_ampCorr - time_ref );
	    }
	}
 
    
      for (int iPosCutX = 0; iPosCutX < s.n_pos_cuts_x; ++iPosCutX)
	{
	  float stepX = (upperPosCutX-lowerPosCutX) / s.n_pos_cuts_x;
      
	  for(int iPosCutY = 0; iPosCutY < s.n_pos_cuts_y; ++iPosCutY)
	    {
	      float stepY = (upperPosCutY-lowerPosCutY) / s.n_pos_cuts_y;
        
	      if( myX > lowerPosCutX + iPosCutX*stepX &&
		  myX < lowerPosCutX + (iPosCutX+1)*stepX &&
		  myY > lowerPosCutY + iPosCutY*stepY &&
		  myY < lowerPosCutY + (iPosCutY+1)*stepY )
		{
		  h_deltat_avg_ampCorr_posCutXY[iPosCutX][iPosCutY]  -> Fill( deltat_avg_ampCorr );
		}
	    }
	}

    
      for(int iBSCut = 0; iBSCut< s.beam_spot_bins; ++iBSCut)
	{
	  if ( fabs(myX-s.centerX) < BScut[iBSCut] )
	    {
	      h_deltat_avg_ampCorr_BSCut[iBSCut] -> Fill( 0.5*(time1_ampCorr+time2_ampCorr) - time_ref );
	      h_deltat_left_ampCorr_BSCut[iBSCut] -> Fill( time1_ampCorr - time_ref );
	      h_deltat_right_ampCorr_BSCut[iBSCut] -> Fill( time2_ampCorr - time_ref );
	    }
	}
    
      ++selectedEntries; 
    }
  
  std::cout << "\n>>> 3rd loop: selected entries " << selectedEntries << std::endl;  

  // Plots, plots, plots, and....moar plots!

  TCanvas* c_time_vs_X = new TCanvas("c_time_vs_X","c_time_vs_X",1000,500);
  c_time_vs_X -> Divide(2,1);
  
  c_time_vs_X -> cd(1);
  gPad->SetGridy();
  
  TF1* fitFunc_corrX = new TF1("fitFunc_corrX","pol2",s.centerX-2.*s.half_beam_spot_x,s.centerX+2.*s.half_beam_spot_x);
  p_avg_vs_X -> Fit(fitFunc_corrX,"QNR");

  p_avg_vs_X->GetXaxis()->SetRangeUser(s.centerX-2.*s.half_beam_spot_x,s.centerX+2.*s.half_beam_spot_x);
  p_avg_vs_X->GetYaxis()->SetRangeUser(h_deltat_avg->GetMean()-15.*h_deltat_avg->GetRMS(),
                                       h_deltat_avg->GetMean()+15.*h_deltat_avg->GetRMS());
  p_avg_vs_X->Draw();
  p_avg_vs_X->SetStats(0);
  p_avg_vs_X ->SetTitle(";X [mm];#Deltat [ns]");
  fitFunc_corrX -> Draw("same");
  p_avg_vs_X-> SetMarkerStyle(21);
  p_left_vs_X-> SetLineColor(kRed+1);
  p_left_vs_X-> SetMarkerColor(kRed+1);
  p_left_vs_X-> SetMarkerStyle(20);
  p_left_vs_X->Draw("same");
  p_right_vs_X-> SetLineColor(kBlue+1);
  p_right_vs_X-> SetMarkerColor(kBlue+1);
  p_right_vs_X-> SetMarkerStyle(20);
  p_right_vs_X->Draw("same");
  p_diff_vs_X-> SetLineColor(kYellow+2);
  p_diff_vs_X-> SetMarkerColor(kYellow+2);
  p_diff_vs_X-> SetMarkerStyle(22);
  p_diff_vs_X-> SetLineStyle(7);
  p_diff_vs_X->Draw("same");
  
  leg = new TLegend(0.71,0.73,0.86,0.93,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);  
  leg->AddEntry(p_left_vs_X,  "t_{left} - t_{MCP}", "lpe");     
  leg->AddEntry(p_right_vs_X, "t_{right} - t_{MCP}", "lpe");     
  leg->AddEntry(p_avg_vs_X,   "t_{avg} - t_{MCP}", "lpe");     
  leg->AddEntry(p_diff_vs_X,  "t_{left} - t_{right}", "lpe");     
  leg->Draw("same");
  
  c_time_vs_X->cd(2);
  gPad->SetGridy();
  
  TF1* fitFunc_corrDiff = new TF1("fitFunc_corrDiff","pol2",h_deltat_diff->GetMean()-5.*h_deltat_diff->GetRMS(),h_deltat_diff->GetMean()+5.*h_deltat_diff->GetRMS());
  p_avg_vs_diff -> Fit(fitFunc_corrDiff,"QNR");
  
  p_avg_vs_diff->GetXaxis()->SetRangeUser(h_deltat_diff->GetMean()-5.*h_deltat_diff->GetRMS(),
                                          h_deltat_diff->GetMean()+5.*h_deltat_diff->GetRMS());
  p_avg_vs_diff->GetYaxis()->SetRangeUser(h_deltat_avg->GetMean()-15.*h_deltat_avg->GetRMS(),
                                          h_deltat_avg->GetMean()+15.*h_deltat_avg->GetRMS());
  p_avg_vs_diff->Draw();
  p_avg_vs_diff->SetStats(0);
  p_avg_vs_diff -> SetTitle(";t_{left} - t_{right} [ns];#Deltat [ns]");
  fitFunc_corrDiff -> Draw("same");
  p_avg_vs_diff-> SetMarkerStyle(20);
  p_left_vs_diff-> SetLineColor(kRed+1);
  p_left_vs_diff-> SetMarkerColor(kRed+1);
  p_left_vs_diff-> SetMarkerStyle(20);
  p_left_vs_diff->Draw("same");
  p_right_vs_diff-> SetLineColor(kBlue+1);
  p_right_vs_diff-> SetMarkerColor(kBlue+1);
  p_right_vs_diff-> SetMarkerStyle(20);
  p_right_vs_diff->Draw("same");
  
  leg = new TLegend(0.71,0.78,0.86,0.93,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->AddEntry(p_left_vs_diff,  "t_{left} - t_{MCP}", "lpe");     
  leg->AddEntry(p_right_vs_diff, "t_{right} - t_{MCP}", "lpe");     
  leg->AddEntry(p_avg_vs_diff,   "t_{avg} - t_{MCP}", "lpe");     
  leg->Draw("same");
  
  c_time_vs_X -> Print(Form("c_time_vs_X_%s.png",cfg.label.c_str()));
  
  
  // ---------------------------- compare left, right, sum time resolution ----------------------------
  TCanvas* c_timeRes_comp = new TCanvas("c_timeRes_comp","c_timeRes_comp",1000,500);
  c_timeRes_comp->Divide(2,1);
  
  c_timeRes_comp->cd(1);
  h_deltat_avg -> SetStats(0);
  h_deltat_avg -> SetTitle(";#Deltat (no amp-walk corr.) [ns];entries");
  h_deltat_avg -> SetLineColor(kBlack);
  h_deltat_avg -> SetLineWidth(2);
  h_deltat_avg -> GetXaxis() -> SetRangeUser(h_deltat_avg->GetMean()-15.*h_deltat_avg->GetRMS(),
                                             h_deltat_avg->GetMean()+15.*h_deltat_avg->GetRMS());
  h_deltat_avg -> Draw();
  h_deltat_left -> Draw("same");
  h_deltat_left -> SetLineColor(kRed+1);
  h_deltat_left -> SetLineWidth(2);
  h_deltat_right -> Draw("same");
  h_deltat_right -> SetLineColor(kBlue+1);
  h_deltat_right -> SetLineWidth(2);
  
  TF1* fitdeltat_left = new TF1("fitdeltat_left", "gaus", h_deltat_left->GetMean()-h_deltat_left->GetRMS()*2, h_deltat_left->GetMean()+h_deltat_left->GetRMS()*2);
  fitdeltat_left->SetLineColor(kRed+1);
  h_deltat_left->Fit(fitdeltat_left, "QR");
  TF1* fitdeltat_right = new TF1("fitdeltat_right", "gaus", h_deltat_right->GetMean()-h_deltat_right->GetRMS()*2, h_deltat_right->GetMean()+h_deltat_right->GetRMS()*2);
  fitdeltat_right->SetLineColor(kBlue+1);
  h_deltat_right->Fit(fitdeltat_right, "QR");
  TF1* fitdeltat_avg = new TF1("fitdeltat_avg", "gaus", h_deltat_avg->GetMean()-h_deltat_avg->GetRMS()*2, h_deltat_avg->GetMean()+h_deltat_avg->GetRMS()*2);
  fitdeltat_avg->SetLineColor(kBlack);
  h_deltat_avg->Fit(fitdeltat_avg, "QR");

  float sigmaLeft  = sqrt( pow(fitdeltat_left->GetParameter(2),2)  - pow(s.mcp_sigma,2) );    
  float sigmaRight = sqrt( pow(fitdeltat_right->GetParameter(2),2) - pow(s.mcp_sigma,2) );    
  float sigmaAvg   = sqrt( pow(fitdeltat_avg->GetParameter(2),2)   - pow(s.mcp_sigma,2) );    

  leg = new TLegend(0.65,0.78,0.80,0.93,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1000);
  leg->AddEntry(h_deltat_left,  Form("#sigma^{left}_{t} = %.1f ps", sigmaLeft*1000), "l");
  leg->AddEntry(h_deltat_right, Form("#sigma^{right}_{t} = %.1f ps", sigmaRight*1000), "l");
  leg->AddEntry(h_deltat_avg,   Form("#sigma^{avg}_{t} = %.1f ps", sigmaAvg*1000), "l");
  leg->Draw("same");
  
  c_timeRes_comp->cd(2);
  h_deltat_avg_ampCorr -> SetStats(0);
  h_deltat_avg_ampCorr -> SetTitle(";#Deltat (amp-walk corr.) [ns];entries");
  h_deltat_avg_ampCorr -> SetLineColor(kBlack);
  h_deltat_avg_ampCorr -> SetLineWidth(2);
  h_deltat_avg_ampCorr -> GetXaxis() -> SetRangeUser(h_deltat_avg->GetMean()-15.*h_deltat_avg->GetRMS(),
                                                     h_deltat_avg->GetMean()+15.*h_deltat_avg->GetRMS());
  h_deltat_avg_ampCorr -> Draw();
  h_deltat_left_ampCorr -> SetLineColor(kRed+1);
  h_deltat_left_ampCorr -> SetLineWidth(2);
  h_deltat_left_ampCorr -> Draw("same");
  h_deltat_right_ampCorr -> SetLineColor(kBlue+1);
  h_deltat_right_ampCorr -> SetLineWidth(2);
  h_deltat_right_ampCorr -> Draw("same");
  
  TF1* fitdeltat_left_ampCorr = new TF1("fitdeltat_left_ampCorr", "gaus", h_deltat_left_ampCorr->GetMean()-h_deltat_left_ampCorr->GetRMS()*2, h_deltat_left_ampCorr->GetMean()+h_deltat_left_ampCorr->GetRMS()*2);
  fitdeltat_left_ampCorr->SetLineColor(kRed+1);
  h_deltat_left_ampCorr->Fit(fitdeltat_left_ampCorr, "QR");
  TF1* fitdeltat_right_ampCorr = new TF1("fitdeltat_right_ampCorr", "gaus", h_deltat_right_ampCorr->GetMean()-h_deltat_right_ampCorr->GetRMS()*2, h_deltat_right_ampCorr->GetMean()+h_deltat_right_ampCorr->GetRMS()*2);
  fitdeltat_right_ampCorr->SetLineColor(kBlue+1);
  h_deltat_right_ampCorr->Fit(fitdeltat_right_ampCorr, "QR");
  TF1* fitdeltat_avg_ampCorr = new TF1("fitdeltat_avg_ampCorr", "gaus", h_deltat_avg_ampCorr->GetMean()-h_deltat_avg_ampCorr->GetRMS()*2, h_deltat_avg_ampCorr->GetMean()+h_deltat_avg_ampCorr->GetRMS()*2);
  fitdeltat_avg_ampCorr->SetLineColor(kBlack);
  h_deltat_avg_ampCorr->Fit(fitdeltat_avg_ampCorr, "QR");
  
  float sigmaLeftCorr  = sqrt(pow(fitdeltat_left_ampCorr->GetParameter(2),2)  - pow(s.mcp_sigma,2) );    
  float sigmaRightCorr = sqrt(pow(fitdeltat_right_ampCorr->GetParameter(2),2) - pow(s.mcp_sigma,2) );    
  float sigmaAvgCorr   = sqrt(pow(fitdeltat_avg_ampCorr->GetParameter(2),2)   - pow(s.mcp_sigma,2) );
  
  
  // ---------------------------- BS cut plots ----------------------------
  TF1* fitdeltat_ampCorr_BSCut_L = new TF1("fitdeltat_ampCorr_BSCut_L", "gaus", h_deltat_left_ampCorr->GetMean() -h_deltat_left_ampCorr->GetRMS()*2,  h_deltat_left_ampCorr->GetMean() +h_deltat_left_ampCorr->GetRMS()*2);
  TF1* fitdeltat_ampCorr_BSCut_R = new TF1("fitdeltat_ampCorr_BSCut_R", "gaus", h_deltat_right_ampCorr->GetMean()-h_deltat_right_ampCorr->GetRMS()*2, h_deltat_right_ampCorr->GetMean()+h_deltat_right_ampCorr->GetRMS()*2);
  TF1* fitdeltat_ampCorr_BSCut   = new TF1("fitdeltat_ampCorr_BSCut",   "gaus", h_deltat_avg_ampCorr->GetMean()  -h_deltat_avg_ampCorr->GetRMS()*2,   h_deltat_avg_ampCorr->GetMean()  +h_deltat_avg_ampCorr->GetRMS()*2);
  
  TGraphErrors * gdeltat_vs_BS_L = new TGraphErrors ();
  TGraphErrors * gdeltat_vs_BS_R = new TGraphErrors ();
  TGraphErrors * gdeltat_vs_BS   = new TGraphErrors ();
  
  int myPoint = 0;
  for (int iBSCut = 0; iBSCut< s.beam_spot_bins; iBSCut++)
  {
    //left
    h_deltat_left_ampCorr_BSCut[iBSCut]->Fit(fitdeltat_ampCorr_BSCut_L, "QNR");
    float tempSigma_L =  sqrt(pow(fitdeltat_ampCorr_BSCut_L->GetParameter(2),2) - pow(s.mcp_sigma, 2));  
    float sigmaErr_L = fitdeltat_ampCorr_BSCut_L->GetParError(2);
    
    //right
    h_deltat_right_ampCorr_BSCut[iBSCut]->Fit(fitdeltat_ampCorr_BSCut_R, "QNR");      
    float tempSigma_R =  sqrt(pow(fitdeltat_ampCorr_BSCut_R->GetParameter(2),2) - pow(s.mcp_sigma, 2));  
    float sigmaErr_R = fitdeltat_ampCorr_BSCut_R->GetParError(2);
    
    //avg
    h_deltat_avg_ampCorr_BSCut[iBSCut]->Fit(fitdeltat_ampCorr_BSCut, "QNR");
    float tempSigma =  sqrt(pow(fitdeltat_ampCorr_BSCut->GetParameter(2),2) - pow(s.mcp_sigma, 2));  
    float sigmaErr = fitdeltat_ampCorr_BSCut->GetParError(2);
    
    if (tempSigma>0 && tempSigma<0.5 &&  h_deltat_avg_ampCorr_BSCut[iBSCut]->GetEntries()>30) 
    {
      gdeltat_vs_BS_L->SetPoint(myPoint, BScut[iBSCut]*2, tempSigma_L);            
      gdeltat_vs_BS_L->SetPointError(myPoint,0, sigmaErr_L);    
      
      gdeltat_vs_BS_R->SetPoint(myPoint, BScut[iBSCut]*2, tempSigma_R);            
      gdeltat_vs_BS_R->SetPointError(myPoint,0, sigmaErr_R);    
      
      gdeltat_vs_BS->SetPoint(myPoint, BScut[iBSCut]*2, tempSigma);            
      gdeltat_vs_BS->SetPointError(myPoint,0, sigmaErr);    
      
      ++myPoint;
    }
  }
  
  TCanvas* c_timeRes_vs_BS = new TCanvas("c_timeRes_vs_BS","c_timeRes_vs_BS",500,500);
  c_timeRes_vs_BS->cd();
  gPad->SetGridy();
  gPad->SetLogx();
  hPad = (TH1F*)( gPad->DrawFrame(BScut[s.beam_spot_bins-1],0.02,3.*BScut[0],0.1) );
  hPad -> SetTitle(";beam spot width [mm];#sigma_{t} [ns]");
  hPad -> Draw();
  gdeltat_vs_BS->Draw("PLE,same");
  gdeltat_vs_BS->SetMarkerStyle(20);
  gdeltat_vs_BS_L->SetLineColor(kRed+1);
  gdeltat_vs_BS_L->SetMarkerColor(kRed+1);
  gdeltat_vs_BS_L->SetMarkerStyle(21);
  gdeltat_vs_BS_L->Draw("same LPE");
  gdeltat_vs_BS_R->SetLineColor(kBlue+1);
  gdeltat_vs_BS_R->SetMarkerColor(kBlue+1);
  gdeltat_vs_BS_R->SetMarkerStyle(21);
  gdeltat_vs_BS_R->Draw("same LPE");
  
  leg = new TLegend(0.71,0.78,0.86,0.93,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);  
  leg->AddEntry(gdeltat_vs_BS_L, Form("left only"), "lpe");     
  leg->AddEntry(gdeltat_vs_BS_R, Form("right only"), "lpe");     
  leg->AddEntry(gdeltat_vs_BS,   Form("avgrage"), "lpe");     
  leg->Draw("same");
  
  c_timeRes_vs_BS -> Print(Form("c_timeRes_vs_BS_%s.png",cfg.label.c_str()));
  
  // ---------------------------- plots in position bins ----------------------------
  TGraphErrors* g_timeRes_left_vs_X = new TGraphErrors ();
  TGraphErrors* g_timeRes_right_vs_X = new TGraphErrors ();
  TGraphErrors* g_timeRes_avg_vs_X = new TGraphErrors ();
  
  TGraphErrors* g_timeRes_left_vs_Y = new TGraphErrors ();
  TGraphErrors* g_timeRes_right_vs_Y = new TGraphErrors ();
  TGraphErrors* g_timeRes_avg_vs_Y   = new TGraphErrors ();
  
  TH2F* h2_timeRes_avg_vs_XY = new TH2F("h2_timeRes_avg_vs_XY","",s.n_pos_cuts_x,lowerPosCutX,upperPosCutX,s.n_pos_cuts_y,lowerPosCutY,upperPosCutY);
  
  // X position
  // avgrage
  myPoint = 0;
  for(int iPosCut = 0; iPosCut < s.n_pos_cuts_x; ++iPosCut)
  {
    float step = ( upperPosCutX - lowerPosCutX ) / s.n_pos_cuts_x;

    TF1* fitdeltat_ampCorr_posCut = new TF1("fitdeltat_avg_ampCorr_posCut","gaus",
                                            h_deltat_avg_ampCorr_posCutX[iPosCut]->GetMean()-h_deltat_avg_ampCorr_posCutX[iPosCut]->GetRMS()*2.,
                                            h_deltat_avg_ampCorr_posCutX[iPosCut]->GetMean()+h_deltat_avg_ampCorr_posCutX[iPosCut]->GetRMS()*2.);
    h_deltat_avg_ampCorr_posCutX[iPosCut]->Fit(fitdeltat_ampCorr_posCut, "QNR");
    float selPos = lowerPosCutX+(step*iPosCut);
    float tempSigma =  sqrt( pow(fitdeltat_ampCorr_posCut->GetParameter(2),2) - pow(s.mcp_sigma,2) );  
    float sigmaErr = fitdeltat_ampCorr_posCut->GetParError(2);      
    if( tempSigma > 0. && tempSigma < 0.5 && h_deltat_avg_ampCorr_posCutX[iPosCut]->GetEntries() > 20 ) 
    {
      g_timeRes_avg_vs_X->SetPoint(myPoint,selPos,tempSigma);            
      g_timeRes_avg_vs_X->SetPointError(myPoint,step/2/sqrt(12),sigmaErr);    
      ++myPoint;
    }

    delete fitdeltat_ampCorr_posCut;
  }

  // left only
  myPoint = 0;
  for(int iPosCut = 0; iPosCut < s.n_pos_cuts_x; ++iPosCut)
    {
      float step = ( upperPosCutX - lowerPosCutX ) / s.n_pos_cuts_x;

      TF1* fitdeltat_ampCorr_posCut = new TF1("fitdeltat_left_ampCorr_posCut","gaus",
                                              h_deltat_left_ampCorr->GetMean()-h_deltat_left_ampCorr->GetRMS()*2.,
                                              h_deltat_left_ampCorr->GetMean()+h_deltat_left_ampCorr->GetRMS()*2.);
      h_deltat_left_ampCorr_posCutX[iPosCut]->Fit(fitdeltat_ampCorr_posCut, "QNR");
      float selPos = lowerPosCutX+(step*iPosCut);
      float tempSigma =  sqrt(pow(fitdeltat_ampCorr_posCut->GetParameter(2),2) - pow(s.mcp_sigma, 2));
      float sigmaErr = fitdeltat_ampCorr_posCut->GetParError(2);
      if( tempSigma > 0. && tempSigma < 0.5 && h_deltat_left_ampCorr_posCutX[iPosCut]->GetEntries() > 20 )
      {
        g_timeRes_left_vs_X->SetPoint(myPoint, selPos, tempSigma);
        g_timeRes_left_vs_X->SetPointError(myPoint,step/2, sigmaErr);
        ++myPoint;
      }
      delete fitdeltat_ampCorr_posCut;
    }

  // right only
  myPoint = 0;
  for(int iPosCut = 0; iPosCut < s.n_pos_cuts_x; ++iPosCut)
  {
    float step = ( upperPosCutX - lowerPosCutX ) / s.n_pos_cuts_x;

    TF1* fitdeltat_ampCorr_posCut = new TF1("fitdeltat_right_ampCorr_posCut","gaus",
                                            h_deltat_right_ampCorr->GetMean()-h_deltat_right_ampCorr->GetRMS()*2.,
                                            h_deltat_right_ampCorr->GetMean()+h_deltat_right_ampCorr->GetRMS()*2.);
    h_deltat_right_ampCorr_posCutX[iPosCut]->Fit(fitdeltat_ampCorr_posCut, "QNR");
    float selPos = lowerPosCutX+(step*iPosCut);
    float tempSigma =  sqrt(pow(fitdeltat_ampCorr_posCut->GetParameter(2),2) - pow(s.mcp_sigma, 2));
    float sigmaErr = fitdeltat_ampCorr_posCut->GetParError(2);
    if (tempSigma > 0. && tempSigma < 0.5 &&  h_deltat_right_ampCorr_posCutX[iPosCut]->GetEntries() > 20 )
    {
      g_timeRes_right_vs_X->SetPoint(myPoint,selPos,tempSigma);
      g_timeRes_right_vs_X->SetPointError(myPoint,step/2,sigmaErr);
      ++myPoint;
    }
    
    delete fitdeltat_ampCorr_posCut;
  }
  
  TCanvas* c_timeRes_vs_X_Y = new TCanvas("c_timeRes_vs_X_Y","c_timeRes_vs_X_Y",1000,500);
  c_timeRes_vs_X_Y -> Divide(2,1);
  c_timeRes_vs_X_Y->cd(1);
  gPad->SetGridy();
  g_timeRes_avg_vs_X->GetYaxis()->SetRangeUser(0, 0.2);
  g_timeRes_avg_vs_X->SetTitle(";X [mm];#sigma_{t} [ns]");
  g_timeRes_avg_vs_X->SetMarkerStyle(20);
  g_timeRes_avg_vs_X->Draw("ALPE");
  g_timeRes_left_vs_X->SetLineColor(kRed+1);
  g_timeRes_left_vs_X->SetMarkerColor(kRed+1);
  g_timeRes_right_vs_X->SetLineColor(kBlue+1);
  g_timeRes_right_vs_X->SetMarkerColor(kBlue+1);
  g_timeRes_left_vs_X->Draw("same LPE");
  g_timeRes_right_vs_X->Draw("same LPE");
  
  leg = new TLegend(0.70,0.78,0.85,0.93,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1000);
  leg->AddEntry(g_timeRes_left_vs_X,  Form("#sigma^{left}_{t}"), "el");
  leg->AddEntry(g_timeRes_right_vs_X, Form("#sigma^{right}_{t}"),"el");
  leg->AddEntry(g_timeRes_avg_vs_X,   Form("#sigma^{avg}_{t}"),  "pl");
  leg->Draw("same");
  
  // Y position
  // avgrage
  myPoint = 0;
  for(int iPosCut = 0; iPosCut < s.n_pos_cuts_y; ++iPosCut)
  {
    float step = ( upperPosCutY - lowerPosCutY ) / s.n_pos_cuts_y;

    TF1* fitdeltat_ampCorr_posCut = new TF1("fitdeltat_avg_ampCorr_posCut","gaus",
                                            h_deltat_avg_ampCorr_posCutY[iPosCut]->GetMean()-h_deltat_avg_ampCorr_posCutY[iPosCut]->GetRMS()*2.,
                                            h_deltat_avg_ampCorr_posCutY[iPosCut]->GetMean()+h_deltat_avg_ampCorr_posCutY[iPosCut]->GetRMS()*2.);
    h_deltat_avg_ampCorr_posCutY[iPosCut]->Fit(fitdeltat_ampCorr_posCut, "QNR");
    float selPos = lowerPosCutY+(step*iPosCut);
    float tempSigma =  sqrt( pow(fitdeltat_ampCorr_posCut->GetParameter(2),2) - pow(s.mcp_sigma,2) );  
    float sigmaErr = fitdeltat_ampCorr_posCut->GetParError(2);      
    if( tempSigma > 0. && tempSigma < 0.5 && h_deltat_avg_ampCorr_posCutY[iPosCut]->GetEntries() > 20 ) 
    {
      g_timeRes_avg_vs_Y->SetPoint(myPoint,selPos,tempSigma);            
      g_timeRes_avg_vs_Y->SetPointError(myPoint,step/2/sqrt(12),sigmaErr);    
      ++myPoint;
    }

    delete fitdeltat_ampCorr_posCut;
  }

  // left only
  myPoint = 0;
  for(int iPosCut = 0; iPosCut < s.n_pos_cuts_y; ++iPosCut)
    {
      float step = ( upperPosCutY - lowerPosCutY ) / s.n_pos_cuts_y;

      TF1* fitdeltat_ampCorr_posCut = new TF1("fitdeltat_left_ampCorr_posCut","gaus",
                                              h_deltat_left_ampCorr->GetMean()-h_deltat_left_ampCorr->GetRMS()*2.,
                                              h_deltat_left_ampCorr->GetMean()+h_deltat_left_ampCorr->GetRMS()*2.);
      h_deltat_left_ampCorr_posCutY[iPosCut]->Fit(fitdeltat_ampCorr_posCut, "QNR");
      float selPos = lowerPosCutY+(step*iPosCut);
      float tempSigma =  sqrt(pow(fitdeltat_ampCorr_posCut->GetParameter(2),2) - pow(s.mcp_sigma, 2));
      float sigmaErr = fitdeltat_ampCorr_posCut->GetParError(2);
      if( tempSigma > 0. && tempSigma < 0.5 && h_deltat_left_ampCorr_posCutY[iPosCut]->GetEntries() > 20 )
      {
        g_timeRes_left_vs_Y->SetPoint(myPoint, selPos, tempSigma);
        g_timeRes_left_vs_Y->SetPointError(myPoint,step/2, sigmaErr);
        ++myPoint;
      }
      delete fitdeltat_ampCorr_posCut;
    }

  // right only
  myPoint = 0;
  for(int iPosCut = 0; iPosCut < s.n_pos_cuts_y; ++iPosCut)
  {
    float step = ( upperPosCutY - lowerPosCutY ) / s.n_pos_cuts_y;

    TF1* fitdeltat_ampCorr_posCut = new TF1("fitdeltat_right_ampCorr_posCut","gaus",
                                            h_deltat_right_ampCorr->GetMean()-h_deltat_right_ampCorr->GetRMS()*2.,
                                            h_deltat_right_ampCorr->GetMean()+h_deltat_right_ampCorr->GetRMS()*2.);
    h_deltat_right_ampCorr_posCutY[iPosCut]->Fit(fitdeltat_ampCorr_posCut, "QNR");
    float selPos = lowerPosCutY+(step*iPosCut);
    float tempSigma =  sqrt(pow(fitdeltat_ampCorr_posCut->GetParameter(2),2) - pow(s.mcp_sigma, 2));
    float sigmaErr = fitdeltat_ampCorr_posCut->GetParError(2);
    if (tempSigma > 0. && tempSigma < 0.5 &&  h_deltat_right_ampCorr_posCutY[iPosCut]->GetEntries() > 20 )
    {
      g_timeRes_right_vs_Y->SetPoint(myPoint,selPos,tempSigma);
      g_timeRes_right_vs_Y->SetPointError(myPoint,step/2,sigmaErr);
      ++myPoint;
    }
    
    delete fitdeltat_ampCorr_posCut;
  }
  
  c_timeRes_vs_X_Y->cd(2);
  gPad->SetGridy();
  g_timeRes_avg_vs_Y->GetYaxis()->SetRangeUser(0, 0.2);
  g_timeRes_avg_vs_Y->SetTitle(";Y [mm];#sigma_{t} [ns]");
  g_timeRes_avg_vs_Y->SetMarkerStyle(20);
  g_timeRes_avg_vs_Y->Draw("ALPE");
  g_timeRes_left_vs_Y->SetLineColor(kRed+1);
  g_timeRes_left_vs_Y->SetMarkerColor(kRed+1);
  g_timeRes_right_vs_Y->SetLineColor(kBlue+1);
  g_timeRes_right_vs_Y->SetMarkerColor(kBlue+1);
  g_timeRes_left_vs_Y->Draw("same LPE");
  g_timeRes_right_vs_Y->Draw("same LPE");
  
  leg = new TLegend(0.70,0.78,0.85,0.93,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1000);
  leg->AddEntry(g_timeRes_left_vs_Y,  Form("#sigma^{left}_{t}"), "el");
  leg->AddEntry(g_timeRes_right_vs_Y, Form("#sigma^{right}_{t}"),"el");
  leg->AddEntry(g_timeRes_avg_vs_Y,   Form("#sigma^{avg}_{t}"),  "pl");
  leg->Draw("same");
  
  c_timeRes_vs_X_Y -> Print(Form("c_timeRes_vs_X_Y_%s.png",cfg.label.c_str()));
  
  
  // XY position
  //right only
  for(int iPosCutX = 0; iPosCutX < s.n_pos_cuts_x; ++iPosCutX)
  {
    float stepX = (upperPosCutX - lowerPosCutX) / s.n_pos_cuts_x;
    
    for(int iPosCutY = 0; iPosCutY < s.n_pos_cuts_y; ++iPosCutY)
    {
      float stepY = (upperPosCutY - lowerPosCutY) / s.n_pos_cuts_y;
      
      if( h_deltat_avg_ampCorr_posCutXY[iPosCutX][iPosCutY]->GetEntries() < 20 ) continue;

      TF1* fitdeltat_ampCorr_posCut = new TF1("fitdeltat_right_ampCorr_posCut","gaus",
                                              h_deltat_avg_ampCorr_posCutXY[iPosCutX][iPosCutY]->GetMean()-h_deltat_avg_ampCorr_posCutXY[iPosCutX][iPosCutY]->GetRMS()*2.,
                                              h_deltat_avg_ampCorr_posCutXY[iPosCutX][iPosCutY]->GetMean()+h_deltat_avg_ampCorr_posCutXY[iPosCutX][iPosCutY]->GetRMS()*2.);
      h_deltat_avg_ampCorr_posCutXY[iPosCutX][iPosCutY]->Fit("fitdeltat_right_ampCorr_posCut","QNR");
      float selPosX = lowerPosCutX+(stepX*iPosCutX);
      float selPosY = lowerPosCutY+(stepY*iPosCutY);
      float tempSigma = sqrt( pow(fitdeltat_ampCorr_posCut->GetParameter(2),2) - pow(s.mcp_sigma,2) );
      // float sigmaErr = fitdeltat_ampCorr_posCut->GetParError(2);
      
      if (tempSigma > 0. && tempSigma < 0.5 && h_deltat_avg_ampCorr_posCutXY[iPosCutX][iPosCutY]->GetEntries() > 20 )
      {
        h2_timeRes_avg_vs_XY -> Fill(selPosX,selPosY,tempSigma);
      }

      delete fitdeltat_ampCorr_posCut;
    }
  }
  
  TCanvas* c_timeRes_vs_XY = new TCanvas ("c_timeRes_vs_XY","c_timeRes_vs_XY",500,500);
  gStyle -> SetPaintTextFormat(".3f");
  h2_timeRes_avg_vs_XY->SetStats(0);
  h2_timeRes_avg_vs_XY->Draw("COLZ,text");
  h2_timeRes_avg_vs_XY->GetZaxis()->SetRangeUser(0, 0.08);
  h2_timeRes_avg_vs_XY->SetTitle(";x [mm];y [mm];#sigma_{t} [ns]");
  c_timeRes_vs_XY -> Print(Form("c_timeRes_vs_XY_%s.png",cfg.label.c_str()));
  

  //----------------------------------------------------------------------
  //                    (4) fourth loop events --> to apply pos correction
  //----------------------------------------------------------------------
  // define histograms
  TH1F* h_deltat_wei_ampCorr = new TH1F("h_deltat_wei_ampCorr","",6000, s.minDeltat, s.maxDeltat);

  TH1F* h_deltat_avg_ampCorr_diffCorr = new TH1F("h_deltat_avg_ampCorr_diffCorr","",6000, s.minDeltat, s.maxDeltat);
  TH1F* h_deltat_avg_ampCorr_posCorr = new TH1F("h_deltat_avg_ampCorr_posCorr","",6000, s.minDeltat, s.maxDeltat);

  selectedEntries = 0;
  for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%100 == 0 ) std::cout << ">>> 4th loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;

      auto &td = cfg.treeData;
      auto &phtk = cfg.channels["photek"];

      tree->GetEntry(entry);
      float myX = td.x_dut[phtk.ampid];
      float myY = td.y_dut[phtk.ampid];
  
      if( myX == -999 || myY == -999 ) continue;
      if( fabs(myX-s.centerX) > s.half_beam_spot_x ||
	  fabs(myY-s.centerY) > s.half_beam_spot_y ) continue;

      float time_ref = td.gaus_mean[phtk.timeid];

      //cut on MCP amp
      if( td.amp[phtk.ampid] < phtk.ampmin ||
	  td.amp[phtk.ampid] > phtk.ampmax )
	continue;
      //cut on MCP time
      if ( time_ref < std::max(phtk.time_peak-phtk.time_sigma*s.sigma_time_cut,s.low_time_cut) ||
	   time_ref > std::min(phtk.time_peak+phtk.time_sigma*s.sigma_time_cut,s.high_time_cut) )
	continue;


      //selected bar 
      if (td.amp[left.ampid] < std::max(left.mip_peak*s.MIP_low_cut,
					left.ampmin)) continue;
      if (td.amp[left.ampid] > std::min(left.mip_peak*s.MIP_high_cut,
					left.ampmax)) continue;

      if (td.amp[right.ampid] < std::max(right.mip_peak*s.MIP_low_cut,
					 right.ampmin)) continue;
      if (td.amp[right.ampid] > std::min(right.mip_peak*s.MIP_high_cut,
					 right.ampmax)) continue;

      if (td.time[left.timeid] < std::max(left.time_peak-
					  left.time_sigma*s.sigma_time_cut,
					  s.low_time_cut)) continue;

      if (td.time[left.timeid] > std::min(left.time_peak+
					  left.time_sigma*s.sigma_time_cut,
					  s.high_time_cut)) continue;

      if (td.time[right.timeid] < std::max(right.time_peak-
					   right.time_sigma*s.sigma_time_cut,
					   s.low_time_cut)) continue;

      if (td.time[right.timeid] > std::min(right.time_peak+
					   right.time_sigma*s.sigma_time_cut,
					   s.high_time_cut)) continue;
    
      float amp1 = td.amp[left.ampid];
      float amp2 = td.amp[right.ampid];
      float time1 = td.time[left.timeid]; 
      float time2 = td.time[right.timeid];


      auto &fitAmpCorr = rootstuff[left.name].fitAmpCorr;
      auto &h_amp_cut = rootstuff[left.name].h_amp_cut; 
      float time1_ampCorr = time1 - fitAmpCorr->Eval(amp1) + fitAmpCorr->Eval(h_amp_cut->GetMean());
      fitAmpCorr = rootstuff[right.name].fitAmpCorr;
      h_amp_cut = rootstuff[right.name].h_amp_cut; 
      float time2_ampCorr = time2 - fitAmpCorr->Eval(amp2) + fitAmpCorr->Eval(h_amp_cut->GetMean());

      float deltat_avg_ampCorr = 0.5*(time1_ampCorr+time2_ampCorr) - time_ref;
      float deltat_wei_ampCorr = ( time1_ampCorr/pow(sigmaLeft,2) + time2_ampCorr/pow(sigmaRight,2) ) / ( 1/pow(sigmaLeft,2) + 1/pow(sigmaRight,2) ) - time_ref;
      float posCorr = -1.*fitFunc_corrX->Eval(myX) + fitFunc_corrX->Eval(s.centerX);
      float diffCorr = -1.*fitFunc_corrDiff->Eval(time2_ampCorr-time1_ampCorr) + fitFunc_corrDiff->Eval(h_deltat_diff->GetMean());
      
      
      h_deltat_wei_ampCorr -> Fill( deltat_wei_ampCorr );
      h_deltat_avg_ampCorr_posCorr -> Fill( deltat_avg_ampCorr + posCorr );
      h_deltat_avg_ampCorr_diffCorr-> Fill( deltat_avg_ampCorr + diffCorr );
    
      ++selectedEntries;
    }

  std::cout << "\n>>> 4th loop: selected entries " << selectedEntries << std::endl;

  c_timeRes_comp->cd(2);
  h_deltat_wei_ampCorr -> SetLineColor(kOrange+1);
  h_deltat_wei_ampCorr -> SetLineWidth(2);
  h_deltat_wei_ampCorr -> Draw("same");
  h_deltat_avg_ampCorr_posCorr -> SetLineColor(kViolet+1);
  h_deltat_avg_ampCorr_posCorr -> SetLineWidth(2);
  h_deltat_avg_ampCorr_posCorr -> Draw("same");  
  
  TF1* fitdeltat_wei_ampCorr = new TF1("fitdeltat_wei_ampCorr", "gaus", h_deltat_wei_ampCorr->GetMean()-h_deltat_wei_ampCorr->GetRMS()*2, h_deltat_wei_ampCorr->GetMean()+h_deltat_wei_ampCorr->GetRMS()*2);
  fitdeltat_wei_ampCorr->SetLineColor(kOrange+1);
  h_deltat_wei_ampCorr->Fit(fitdeltat_wei_ampCorr, "QR");
  
  TF1* fitdeltat_avg_ampCorr_posCorr = new TF1("fitdeltat_avg_ampCorr_posCorr", "gaus", h_deltat_avg_ampCorr_posCorr->GetMean()-h_deltat_avg_ampCorr_posCorr->GetRMS()*2, h_deltat_avg_ampCorr_posCorr->GetMean()+h_deltat_avg_ampCorr_posCorr->GetRMS()*2);
  fitdeltat_avg_ampCorr_posCorr->SetLineColor(kViolet+1);
  h_deltat_avg_ampCorr_posCorr->Fit(fitdeltat_avg_ampCorr_posCorr, "QR");
  
  float sigmaWeiCorr    = sqrt(pow(fitdeltat_wei_ampCorr->GetParameter(2),2)         - pow(s.mcp_sigma,2) );
  float sigmaAvgCorrPos = sqrt(pow(fitdeltat_avg_ampCorr_posCorr->GetParameter(2),2) - pow(s.mcp_sigma,2) );
  
  leg = new TLegend(0.65,0.68,0.80,0.93,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);  
  leg->AddEntry(h_deltat_left_ampCorr,       Form("#sigma^{left}_{t} = %.1f ps", sigmaLeftCorr*1000), "l");     
  leg->AddEntry(h_deltat_right_ampCorr,      Form("#sigma^{right}_{t} = %.1f ps", sigmaRightCorr*1000), "l");
  leg->AddEntry(h_deltat_avg_ampCorr,        Form("#sigma^{avg}_{t} = %.1f ps", sigmaAvgCorr*1000), "l");
  leg->AddEntry(h_deltat_wei_ampCorr,        Form("#sigma^{wei}_{t} = %.1f ps", sigmaWeiCorr*1000), "l");
  leg->AddEntry(h_deltat_avg_ampCorr_posCorr,Form("#sigma^{avg+pos}_{t} = %.1f ps", sigmaAvgCorrPos*1000), "l");
  leg->Draw("same");
  
  c_timeRes_comp -> Print(Form("c_timeRes_comp_%s.png",cfg.label.c_str()));
  
  return 0;

}
