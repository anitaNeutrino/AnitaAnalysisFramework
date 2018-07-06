#ifndef ANITA_ANALYSIS_DIODE_FILTER_H
#define ANITA_ANALYSIS_DIODE_FILTER_H

#include "FilterOperation.h"
#include "TString.h" 
#include "AnalysisWaveform.h" 

/* Implements the time domain diode model */ 


class DiodeFilter : public UniformFilterOperation
{

  public: 

    DiodeFilter(); 

    const char * tag () const { return "DiodeFilter"; } 
    const char * description () const { return "DiodeFilter"; } 

    virtual void processOne(AnalysisWaveform *, const RawAnitaHeader *, int, int) ;

  private: 
    AnalysisWaveform response; 





}; 



#endif


