#ifndef IMPULSIVITY_MEASURE_H
#define IMPULSIVITY_MEASURE_H


/** 

  In this file are a number of attempts to define the impulsivity of a signal 
  Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 

*/ 

class TH2; 
class AnalysisWaveform; 
class TGraph; 

namespace impulsivity
{

  /** Computes what I call the Impulsivity Measure of a signal
   *
   *  This is based on something I call a distance CDF, which starts at the
   *  peak of the signal and integrates normalized power outward (on both sides). Points
   *  outside the waveform are considered to be zero. This makes something sort of akin to an ROC curve for impulsivity.  
   *
   *  If all the power is contained within a single sample, the distance CDF 
   *  will just be flat with a value of 1. 
   *  
   *  If all the power is evenly distributed, the distance CDF will be linear with a y-mean of 0.5. 
   *
   *  Therefore, the impulsivity measure is defined as twice the y-mean of the distance CDF minus one. 
   *
   *  
   *  @param wf The waveform to consider
   *  @param distance_cdf If you want to look at the distance cdf, pass a pointer to a TGraph here. It will be resized appropriately. 
   *  @param pt The point to integrate from, otherwise the max abs(value) will be used. 
   *  @param use_envelope True to use the Hilbert Envelope of the waveform to integrate (recommended). Otherwise uses evenly-sampled versior. 
   *  @returns The impulsivity measure 
   */ 
  double impulsivityMeasure(const AnalysisWaveform *wf, TGraph * distance_cdf =
      0, int pt = -1, bool use_envelope = true); 
  
  

  /** The envelopogram computes the RMS envelope of the waveform with different averaging widths 
   *
   *  @param wf the input waveform
   *  @param out use this TH2 instead of allocating a new one 
   *  @param min_width the minimum averaging width to consider
   *  @param max_width the maximum averaging width to considert
   *  @param step The step size of the averaging width
   *  @param use_hilbert Use the hilbert envelope rather than the evenly sampled waveform. 
   *  @return the envelopogram
   *
   */
  TH2 * envelopogram(const AnalysisWaveform *wf , TH2 * out = 0,  int min_width = 1, int max_width = 15, 
                     int step= 1, bool use_hilbert_envelope = true); 

}




#endif
