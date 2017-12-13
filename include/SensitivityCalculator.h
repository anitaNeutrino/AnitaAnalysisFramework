#ifndef SENSITIVITY_CALCULATOR_HH
#define SENSITIVITY_CALCULATOR_HH

#ifdef USE_ROOSTATS
#include "RooWorkspace.h" 
#endif 
 

class TGraphErrors; 

/** Sensitivity calculator ...
 *
 * Poisson signal + gaussian efficiency + separate poisson backgrounds using RooStats 
 *
 * This can easily be generalized to different efficiency / background models pretty easily (although, eventually you'll just be implementing an interface to RooStats) . 
 *
 * This can either use the ``large-sample'' profile likelihood method as described in: 
 * https://arxiv.org/abs/1505.07027 or a profile-likelihood-weighted FC-like method (which requires a large amount of MC) 
 * 
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 */


class SensitivityCalculator 
{


  public: 

    SensitivityCalculator(double CL = 0.9,bool FC = true) ;
      

    void setEfficiency(double mean, double sigma); 
    void setThermalBackground(int nobserved, double bg_to_sig_ratio); 
    void setAnthroBackground(int nobserved, double bg_to_sig_ratio); 
    void getLimit(int nobs, double * lower = 0, double * upper = 0); 
    TGraphErrors * confidenceBands(int start = 0, int stop = 20); 
  private: 

    double cl;
    bool use_fc;

#ifdef USE_ROOSTATS
    RooWorkspace w; 
#endif

};



#endif
