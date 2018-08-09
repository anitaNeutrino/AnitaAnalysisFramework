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
      StokesAnalysis(const StokesAnalysis & other); 
      ~StokesAnalysis() { ; } 

      /** This computes the windowed averages over the window around Imax where I/Imax >= minIfrac 
       *  It returns the number of points used. 
       * */ 
      int computeWindowedAverage(double minIfrac, double * I = 0, double * Q = 0, double * U = 0 ,double * V = 0) const; 


      TGraph  & instI() { return *dI; } 
      TGraph  & instQ() { return *dQ; } 
      TGraph  & instU() { return *dU; } 
      TGraph  & instV() { return *dV; } 

      TGraph  & cumuI() { return *cI; } 
      TGraph  & cumuQ() { return *cQ; } 
      TGraph  & cumuU() { return *cU; } 
      TGraph  & cumuV() { return *cV; } 

      double getAvgI() const { return avgI; } 
      double getAvgQ() const { return avgQ; } 
      double getAvgU() const { return avgU; } 
      double getAvgV() const { return avgV; } 

      void getAvgs(double * I, double * Q, double * U, double * V) const
      {
        *I = avgI; *Q=avgQ; *U=avgU; *V=avgV;
      }

      TMultiGraph & instGraphs() { return instantaneous; } 
      TMultiGraph & cumuGraphs() { return cumulative; } 

    private: 
      TGraph * dI, *dQ, *dU, *dV; 
      TGraph * cI, *cQ, *cU, *cV; 
      double avgI, avgQ, avgU, avgV; 
      TMultiGraph instantaneous; 
      TMultiGraph cumulative; 


  }; 



}

#endif
