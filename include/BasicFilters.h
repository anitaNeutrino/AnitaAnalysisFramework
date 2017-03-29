#ifndef ANITA_ANALYSIS_BASIC_FILTERS
#define ANITA_ANALYSIS_BASIC_FILTERS

#include "FilterOperation.h"
#include "TString.h" 

/** A set of basic filter operations that serve as examples and also should be quite useful */ 


/** SimplePassBandFilter just cuts everything outside the pass band in
 * fourier space like a brick wall, with all the ensuing acausal goodness.
 *
 */ 


namespace FFTtools
{
  class DigitalFilter; 
}

class SimplePassBandFilter : public UniformFilterOperation 
{ 
  public: 
    const char * tag() const { return "SimplePassBandFilter"; } 

    const char * description() const { return descStr.Data(); }

    /** Filter everything outside the pass band. Numbers are given in GHz. */ 
    SimplePassBandFilter(double low = 0.2, double high = 1.2)  
      : low(low), high(high)
    {
        descStr = TString::Format("SimplePassbandFilter(%g,%g)", low,high); 
    }

    virtual void processOne(AnalysisWaveform *) ;

  private: 
    TString descStr; 
    double low; 
    double high;

}; 

/** Brick wall notch filter.  */ 
class SimpleNotchFilter : public UniformFilterOperation
{
  public: 

    /** Ghz*/ 
    SimpleNotchFilter(double minfreq, double maxfreq) 
      : min(minfreq), max(maxfreq) 
    {
      desc.Form("SimpleNotchFilter(%g,%g)",min,max); 
    }

    const char * tag() const { return "SimpleNotchFilter"; } 
    const char * description() const { return desc.Data(); } 
    virtual void processOne(AnalysisWaveform *); 
  private: 
    TString desc;  
    double min,max; 
}; 


/** ALFA filter */ 


class HybridFilter : public FilterOperation
{

  public:

    const char * tag () const { return "HybridFilter"; } 
    const char * description () const { return "Hybrid Filter"; } 
    virtual void process(FilteredAnitaEvent * event); 
};

class SumDifferenceFilter : public FilterOperation
{

  public:

    const char * tag () const { return "SumDifferenceFilter"; } 
    const char * description () const { return "SumDifference Filter"; } 
    virtual void process(FilteredAnitaEvent * event); 
};


class DigitalFilterOperation : public UniformFilterOperation 
{
  public: 
    DigitalFilterOperation(const FFTtools::DigitalFilter *digi) : digi(digi) {;}
    const char * tag () const { return "DigitalFilter"; } 
    const char * description () const { return "DigitalFilter"; } 
    virtual void processOne(AnalysisWaveform* wf); 

  private: 
    const FFTtools::DigitalFilter * digi; 


}; 

class ALFAFilter : public FilterOperation
{
  public: 
    ALFAFilter (double cutoff = 0.55); /** The cutoff is scaled by 1.3, so if you oversampled, this won't be right anymore */ 
    virtual ~ALFAFilter(); 

    virtual void process(FilteredAnitaEvent * event) ; 

    const char * tag() const { return "ALFAFilter"; } 
    const char * description() const { return descStr.Data(); } 
  private:
    DigitalFilterOperation *pb; 
    FFTtools::DigitalFilter * filt;
    double power_before; 
    double power_after;  
    TString descStr;

}; 


#endif 
