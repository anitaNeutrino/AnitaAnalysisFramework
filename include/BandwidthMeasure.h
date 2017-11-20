#ifndef BANDWIDTH_MEASURE_H
#define BANDWIDTH_MEASURE_H


/** 

  In this file is an attempt to define a measure of the bandwidth of a signal
	Andrew Ludwig (  abl@uchicago.edu )

*/ 

class AnalysisWaveform; 
class TGraph; 

namespace bandwidth 
{

  double bandwidthMeasure(const AnalysisWaveform * wf, TGraph* testGraph = 0);
  double alternateBandwidthMeasure(const AnalysisWaveform * wf, double powerThreshold = .5, TGraph* testGraph = 0);

}




#endif
