#ifndef ANITA_ANALYSIS_TAPER_FILTER_H
#define ANITA_ANALYSIS_TAPER_FILTER_H

#include "FilterOperation.h" 
#include "TString.h" 

/* These are filters that taper (window?) the waveform. */

class GaussianTaper : public UniformFilterOperation
{
  public: 

    const char * tag () const { return "GaussianTaper"; } 

    const char * description() const { return descStr.Data(); }

    /** Filter everything outside the pass band. Numbers are given in GHz. */ 
    GaussianTaper(double distance_ns = 5, double nsigma = 2) 
    : mean(distance_ns), sigma(distance_ns/nsigma)
    {
        descStr = TString::Format("GaussianTaper(%g,%g)", distance_ns,nsigma); 
    }

    virtual void processOne(AnalysisWaveform *) ;

  private: 
    TString descStr; 
    double mean; 
    double sigma;



}; 


#endif

