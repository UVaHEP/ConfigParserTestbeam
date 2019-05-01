//Config File processor using boost

#include <iostream>
#include <fstream>
#include <algorithm>
#include "config.h"

namespace po = boost::program_options;

using std::string; 
using boost::algorithm::split; 
  

void channelSettings::print() {
  std::cout << this->name << std::endl
	    << this->ampid << std::endl
	    << this->timeid << std::endl
	    << this->ampmin << std::endl
	    << this->ampmax << std::endl;
}



void selections::print() {

  std::cout << "MIP_low_cut:" << this->MIP_low_cut << std::endl
	    << "MIP_high_cut:" << this->MIP_high_cut << std::endl
	    << "low_time_cut:" << this->low_time_cut << std::endl
	    << "high_time_cut:" << this->high_time_cut << std::endl
	    << "sigma_time_cut:" << this->sigma_time_cut << std::endl
	    << "minX:" << this->minX << std::endl
	    << "maxX:" << this->maxX << std::endl
	    << "centerX:" << this->centerX << std::endl
	    << "half_beam_spot_x:" << this->half_beam_spot_x << std::endl
	    << "n_pos_cuts_x:" << this->n_pos_cuts_x << std::endl
	    << "minY:" << this->minY << std::endl
	    << "maxY:" << this->maxY << std::endl
	    << "centerY:" << this->centerY << std::endl
	    << "half_beam_spot_y:" << this->half_beam_spot_y << std::endl
	    << "n_pos_cuts_y:" << this->n_pos_cuts_y << std::endl
	    << "beam_spot_bins:" << this->beam_spot_bins << std::endl
	    << "mcp_sigma:" << this->mcp_sigma << std::endl
	    << "nAmpBins:" << this->nAmpBins << std::endl
	    << "ampMin:" << this->ampMin << std::endl
	    << "ampMax:" << this->ampMax << std::endl
	    << "nTimeBins:" << this->nTimeBins << std::endl
	    << "minTime:" << this->minTime << std::endl
	    << "maxTime:" << this->maxTime << std::endl
	    << "nDeltatBins:" << this->nDeltatBins << std::endl
	    << "minDeltat:" << this->minDeltat << std::endl
	    << "maxDeltat:" << this->maxDeltat << std::endl; 
}

config::config(string filename) {

  success = false;
  
  po::options_description cdesc ("Config File Options");
  cdesc.add_options()
    ("general.branches", po::value<string>(), "branches list")
    ("general.startrun", po::value<int>(), "run to begin at")
    ("general.stoprun", po::value<int>(), "run to stop at")
    ("general.label", po::value<string>(), "label for job")
    ("general.path", po::value<string>(), "ROOT path to get files from")
    ("channel.names", po::value<string>(), "channels")
    ("channel.ampids", po::value<string>(), "ampids")
    ("channel.timeids", po::value<string>(), "timeids")
    ("channel.ampmin", po::value<string>(), "ampmin")
    ("channel.ampmax", po::value<string>(), "ampmax")
    ("processing.algorithm", po::value<string>(), "processing algo to use")
    ("selections.MIP_low_cut", po::value<float>(), "low amp cut in fraction of MIP Peak")
    ("selections.MIP_high_cut", po::value<float>(), "high amp cut in fraction of MIP peak")
    ("selections.low_time_cut", po::value<float>(), "low time cut in ns")
    ("selections.high_time_cut", po::value<float>(), "high time cut in ns")
    ("selections.sigma_time_cut", po::value<float>(), "n of sigma on time cut")
    ("selections.minX", po::value<float>(), "range of X in mm")
    ("selections.maxX", po::value<float>(), "range of X in mm")
    ("selections.centerX", po::value<float>(), "hodoscope X coordinate of crystal center in mm")
    ("selections.half_beam_spot_x", po::value<float>(), "half-size of beam spot selection around the center (mm)")
    ("selections.n_pos_cuts_x", po::value<int>(), "number of bins along X for binned time resolution")
    ("selections.minY", po::value<float>(), "range of Y in mm")
    ("selections.maxY", po::value<float>(), "range of Y in mm")
    ("selections.centerY", po::value<float>(), "hodoscope Y coordinate of crystal center in mm")
    ("selections.half_beam_spot_y", po::value<float>(), "half-size of beam spot selection around the center in mm")
    ("selections.n_pos_cuts_y", po::value<int>(), "number of bins along Y for binned time resolution")
    ("selections.beam_spot_bins", po::value<int>(), "number of beam spot bins")
    ("selections.mcp_sigma", po::value<float>(), "MCP time resolution")
    ("selections.nAmpBins", po::value<int>(), "")
    ("selections.ampMin", po::value<float>(), "")
    ("selections.ampMax", po::value<float>(), "")
    ("selections.nTimeBins", po::value<int>(), "")
    ("selections.minTime", po::value<float>(), "")
    ("selections.maxTime", po::value<float>(), "")
    ("selections.nDeltatBins", po::value<int>(), "")
    ("selections.minDeltat", po::value<float>(), "")
    ("selections.maxDeltat", po::value<float>(), "")
    
    ;
    
  po::variables_map cmap; //Config file map 

  
  try {

    std::fstream fin(filename);
    po::store(po::parse_config_file(fin, cdesc, true), cmap);
    po::notify(cmap);

    
    brlst = cmap["general.branches"].as<string>();
    path =  cmap["general.path"].as<string>();
    label = cmap["general.label"].as<string>();
    start = cmap["general.startrun"].as<int>();
    stop = cmap["general.stoprun"].as<int>();
    string chlst = cmap["channel.names"].as<string>();
    string ampids = cmap["channel.ampids"].as<string>();
    string timeids = cmap["channel.timeids"].as<string>();
    string ampmins = cmap["channel.ampmin"].as<string>();
    string ampmaxes = cmap["channel.ampmax"].as<string>();
    

    std::vector<float> amplst;
    std::vector<float> amaxlst;
    std::vector<float> aminlst;
    std::vector<float> timelst;

    std::vector<string> names;
    std::vector<string> lst;
    split(names, chlst, boost::is_any_of(";"));
    split(lst, ampids, boost::is_any_of(",")); 
    std::for_each(lst.begin(), lst.end(), [&amplst](string s) { amplst.push_back(atof(s.c_str())); });
    lst.clear();

    split(lst, ampmaxes, boost::is_any_of(",")); 
    std::for_each(lst.begin(), lst.end(), [&amaxlst](string s) { amaxlst.push_back(atof(s.c_str())); });
    lst.clear();

    split(lst, ampmins, boost::is_any_of(",")); 
    std::for_each(lst.begin(), lst.end(), [&aminlst](string s) { aminlst.push_back(atof(s.c_str())); });
    lst.clear();

    split(lst, timeids, boost::is_any_of(",")); 
    std::for_each(lst.begin(), lst.end(), [&timelst](string s) { timelst.push_back(atof(s.c_str())); });
    lst.clear();


    int count = names.size();
    

    if ((count != amaxlst.size()) ||
	(count != aminlst.size()) ||
	(count != timelst.size()) ||
	(count != amplst.size())) {
	std::cout << "Error! In our config file the length for names, ampids, etc, differs" << std::endl;
	std::cout << "Sizes" << std::endl << "--------" << std::endl
		  << "Amp IDs:" << "\t" << amplst.size() << std::endl
		  << "Amp Max:" <<  "\t\t" <<  amaxlst.size() << std::endl
		  << "Amp Min:" << "\t" << aminlst.size() << std::endl
		  << "Time Ids:" <<"\t" << timelst.size() << std::endl
		  << "Names: " << "\t\t" << names.size() << std::endl << std::endl;
	std::vector<unsigned long> szs = {names.size(), amaxlst.size(), aminlst.size(), timelst.size()}; 
	count = *std::min_element(szs.begin(), szs.end()); 
      }


    std::vector<channelSettings> chSettings; 
    for (int i = 0; i < count; i++) {
      channelSettings ch;
      ch.name = names[i];
      ch.ampid = amplst[i];
      ch.timeid = timelst[i];
      ch.ampmin = aminlst[i];
      ch.ampmax = amaxlst[i]; 
      channels[ch.name] = ch; 
    }

	
    for (auto &ch : chSettings) {
      ch.print();
    }

    // Reading and Processing selection criteria will get  moved to its own function

    //Note, I generated most of the following, didn't feel like typing it out
    //I'll try to write a script in python to autogenerate it from the ini file
    //'Cause ain't nobody got time to write all this out! 
    
    cuts.MIP_low_cut = cmap["selections.MIP_low_cut"].as<float>();
    cuts.MIP_high_cut = cmap["selections.MIP_high_cut"].as<float>();
    cuts.low_time_cut = cmap["selections.low_time_cut"].as<float>();
    cuts.high_time_cut = cmap["selections.high_time_cut"].as<float>();
    cuts.sigma_time_cut = cmap["selections.sigma_time_cut"].as<float>();
    cuts.minX = cmap["selections.minX"].as<float>();
    cuts.maxX = cmap["selections.maxX"].as<float>();
    cuts.centerX = cmap["selections.centerX"].as<float>();
    cuts.half_beam_spot_x = cmap["selections.half_beam_spot_x"].as<float>();
    cuts.n_pos_cuts_x = cmap["selections.n_pos_cuts_x"].as<int>();
    cuts.minY = cmap["selections.minY"].as<float>();
    cuts.maxY = cmap["selections.maxY"].as<float>();
    cuts.centerY = cmap["selections.centerY"].as<float>();
    cuts.half_beam_spot_y = cmap["selections.half_beam_spot_y"].as<float>();
    cuts.n_pos_cuts_y = cmap["selections.n_pos_cuts_y"].as<int>();
    cuts.beam_spot_bins = cmap["selections.beam_spot_bins"].as<int>();
    cuts.mcp_sigma = cmap["selections.mcp_sigma"].as<float>();
    cuts.nAmpBins = cmap["selections.nAmpBins"].as<int>();
    cuts.ampMin = cmap["selections.ampMin"].as<float>();
    cuts.ampMax = cmap["selections.ampMax"].as<float>();
    cuts.nTimeBins = cmap["selections.nTimeBins"].as<int>();
    cuts.minTime = cmap["selections.minTime"].as<float>();
    cuts.maxTime = cmap["selections.maxTime"].as<float>();
    cuts.nDeltatBins = cmap["selections.nDeltatBins"].as<int>();
    cuts.minDeltat = cmap["selections.minDeltat"].as<float>();
    cuts.maxDeltat = cmap["selections.maxDeltat"].as<float>();


  }
  catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    return; 
  }
  
  success = true; 
}


