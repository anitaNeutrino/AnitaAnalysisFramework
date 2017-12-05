#ifndef BANDWIDTH_MEASURE_H
#define BANDWIDTH_MEASURE_H

#include <vector>


/** 

  In this file are a few tries at defining a measure of the bandwidth of a signal
	Andrew Ludwig (  abl@uchicago.edu )

*/ 

class AnalysisWaveform; 
class TGraph;
class TGraphAligned; 

namespace bandwidth 
{
 /** Computes the bandwidth measure.  
  *  Was originally based on the impulsivity measure Cosmin wrote, and is calculated in nearly the same way.
  *  Turns out to be essentially the same thing as a Gini coefficient/index in economics.
  *  Takes the normalized power spectrum of a signal and sorts it from highest to lowest power, then makes a CDF of that.
  *  Returns abs(cdf_mean - 1) * 2
  *  This is just so all power in a single frequency gives a bandwidth measure of 0
  *  A flat power spectrum gives a bandwidth measure of 1
  */
  double bandwidthMeasure(const AnalysisWaveform * wf, int timeCheck = 0, TGraph* gTest = 0);
  /** This takes the normalized power spectrum from the input wave form and compares it to the normalized power spectrum of the average impulse response.
   *  Lower number (minimum zero) means more similar to the impulse response.
   *  Also runs from 0->1 with higher being more similar to an impulse response (more broadband)
   */
  double differenceFromImpulse(const AnalysisWaveform * wf, int timeCheck = 0, TGraph* gTest = 0);
  /** This takes the normalized power spectrum from the input wave form and compares it to the normalized power spectrum of the average impulse response.
   * Returns the maximum difference from an impulse response.
   * The idea is to use this to find either strong cw contamination that isnt filtered out or payload blasts
   */
  double maxDifferenceFromImpulse(const AnalysisWaveform * wf, int timeCheck = 0, TGraph* gTest = 0);
  /** Returns the Hoover index of a signal.
   *  Also from economics, is basically a k-s test vs a perfectly distributed frequency spectrum.
   *  Between 0 and 1 with 0 less broadband, 1 more so.
   */
  double hooverIndex(const AnalysisWaveform * wf, int timeCheck = 0);
  /** Returns the Theil index of a signal.
   *  Another economics metric, originally devised as a measure of how ordered a distribution is.
   *  Between 0 and 1 with 0 less broadband, 1 more so.
   */
  double theilIndex(const AnalysisWaveform * wf, int timeCheck = 0);

  /** This stuff is all used by the various metrics, not actually useful for anything else */
  void checkNotches(int timeCheck, double& notch0, double& notch1, double& notch2);
  double fillPowers(const TGraphAligned* powd, std::vector<double> &powers, double notch0, double notch1, double notch2);
  TGraph* loadImpulsePower(int timeCheck);
  void normalizePower(TGraph* g);
  TGraph* downsampleImpulse(TGraph* imp, const TGraphAligned* examp);

}




#endif
