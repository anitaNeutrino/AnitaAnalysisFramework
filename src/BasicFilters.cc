#include "BasicFilters.h"
#include "AnalysisWaveform.h" 
#include "DigitalFilter.h" 
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


void HybridFilter::process(FilteredAnitaEvent * event) 
{
#ifdef USE_OMP
#pragma omp  parallel for 
#endif
  for (int i = 0; i < NUM_SEAVEYS; i++) 
  {
    AnalysisWaveform::basisChange(getWf(event,i, AnitaPol::kHorizontal), getWf(event,i,AnitaPol::kVertical)); 
  }


}

void SumDifferenceFilter::process(FilteredAnitaEvent * event) 
{
#ifdef USE_OMP
#pragma omp  parallel for 
#endif
  for (int i = 0; i < NUM_SEAVEYS; i++) 
  {
    AnalysisWaveform::sumDifference(getWf(event,i, AnitaPol::kHorizontal), getWf(event,i,AnitaPol::kVertical)); 
  }


}


void DigitalFilterOperation::processOne(AnalysisWaveform * g) 
{
  digi->filterGraph(g->updateEven()); 
}


ALFAFilter::ALFAFilter(double cutoff)
{
  // cutoff is a fraction of the Nyquist frequency.
  // (that's where the factor of 1.3 here comes from)  
  // So this won't work properly for interpolated or padded waveforms
  filt = new FFTtools::ButterworthFilter(FFTtools::LOWPASS, 4, cutoff/1.3); 
  pb = new DigitalFilterOperation(filt);
  descStr = TString::Format("ALFA filter - Butterworth low pass at %3.0lf MHz for 5TH and 13TH (assumes f_{Nyquist} = 1300 MHz)", 1e3*cutoff);
}

void ALFAFilter::process(FilteredAnitaEvent *event)
{
   if (AnitaVersion::get() != 3 ) return; 

   pb->processOne( getWf(event, 4, AnitaPol::kHorizontal) ); 
   //cross talk is strong in this one 
   pb->processOne( getWf(event, 12, AnitaPol::kHorizontal) ); 
}


ALFAFilter::~ALFAFilter()
{
  delete pb; 
  delete filt;
}



