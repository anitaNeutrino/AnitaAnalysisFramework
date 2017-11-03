#ifndef _POLARIMETRY_HH
#define _POLARIMETRY_HH

/* some polarimetry utilities */ 

#include "TMultiGraph.h" 

class AnalysisWaveform; 

namespace polarimetry 
{

  class StokesAnalysis 
  {

    public:

      StokesAnalysis(const AnalysisWaveform * H, const AnalysisWaveform *V); 
      ~StokesAnalysis() { ; } 

      /** This computes the windowed averages over the window around Imax where I/Imax >= minIfrac */ 
      int computeWindowedAverage(double minIfrac, double * I = 0, double * Q = 0, double * U = 0 ,double * V = 0) const; 

      TMultiGraph & instantaneousGraph() { return instantaneous; } 
      TMultiGraph & cumulativeGraph() { return cumulative; } 

    private: 
      TGraph * dI, *dQ, *dU, *dV; 
      TGraph * cI, *cQ, *cU, *cV; 
      TMultiGraph instantaneous; 
      TMultiGraph cumulative; 


  }; 



}

#endif
