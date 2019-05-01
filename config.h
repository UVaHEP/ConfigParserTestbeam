#ifndef _CONFIG_H
#define _CONFIG_H


#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <map>
#include <string> 

struct selections {


  void print();
  
  float MIP_low_cut;
  float MIP_high_cut;
  float low_time_cut;
  float high_time_cut;
  float sigma_time_cut;
  float minX;
  float maxX;
  float centerX;
  float half_beam_spot_x;
  int n_pos_cuts_x;

  float minY;
  float maxY;
  float centerY;
  float half_beam_spot_y;

  int n_pos_cuts_y;
  float beam_spot_bins;
  float mcp_sigma;

  int nAmpBins;
  float ampMin;
  float ampMax;

  int nTimeBins;
  float minTime;
  float maxTime;

  int nDeltatBins;
  float minDeltat;
  float maxDeltat;

  
  std::vector<float> BScut;


};


struct channelSettings {

  void print(); 
  
  std::string name;
  unsigned int ampid;
  unsigned int timeid;
  float ampmin;
  float ampmax;
  float mip_peak;
  float time_peak;
  float time_sigma; 
}; 



//Will probably want to move this out and into a seperate object
//As it doesn't fit the name "config" 
struct treeBranches {

  float amp[36];
  float time[36];
  float gaus_mean[36];
  float x_dut[4];
  float y_dut[4];

}; 

typedef std::map< std::string, channelSettings > channelMap; 

struct config {

 public:
  
  config(std::string filename);

  std::string path;
  std::string label; 
  std::string brlst;
  int start;
  int stop; 
  channelMap channels; 
  treeBranches treeData;
  selections cuts; 
  bool success; 

}; 



#endif 
