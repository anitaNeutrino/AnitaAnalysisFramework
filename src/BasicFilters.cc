#include "BasicFilters.h"
#include "AnalysisWaveform.h" 
#include "FFTWComplex.h" 

void SimplePassBandFilter::processOne(AnalysisWaveform* g) 
{
  int nfreq = g->Nfreq(); 
  double df = g->deltaF(); 

//  printf("SimplePassBandFilter::processOne!\n"); 
  FFTWComplex * vals = g->updateFreq(); 
  for (int i = 0; i < nfreq; i++) 
  {
    if (i * df < low ||  i * df > high) 
    {
      vals[i].re = 0; 
      vals[i].im = 0; 
    }
  }
}


void SimpleNotchFilter::processOne(AnalysisWaveform* g) 
{
//  printf("SimpleNotchFilter::processOne!\n"); 
  int nfreq = g->Nfreq(); 
  double df = g->deltaF(); 

  FFTWComplex * vals = g->updateFreq(); 
  for (int i = 0; i < nfreq; i++) 
  {
    if (i * df >= min &&  i * df <= max) 
    {
      vals[i].re = 0; 
      vals[i].im = 0; 
    }
  }
}



