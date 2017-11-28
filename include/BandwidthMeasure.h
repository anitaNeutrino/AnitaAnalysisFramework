#ifndef BANDWIDTH_MEASURE_H
#define BANDWIDTH_MEASURE_H

#include <vector>


/** 

  In this file is an attempt to define a measure of the bandwidth of a signal
	Andrew Ludwig (  abl@uchicago.edu )

*/ 

class AnalysisWaveform; 
class TGraphAligned; 

namespace bandwidth 
{

  void checkNotches(int timeCheck, double& notch0, double& notch1, double& notch2);
  double fillPowers(const TGraphAligned* powd, std::vector<double> &powers, double notch0, double notch1, double notch2);

  double bandwidthMeasure(const AnalysisWaveform * wf, int timeCheck = 0);
  double giniIndex(const AnalysisWaveform * wf, int timeCheck = 0);
  double hooverIndex(const AnalysisWaveform * wf, int timeCheck = 0);
  double theilIndex(const AnalysisWaveform * wf, int timeCheck = 0);
  double alternateBandwidthMeasure(const AnalysisWaveform * wf, int timeCheck = 0, double powerThreshold = .5);

}




#endif
