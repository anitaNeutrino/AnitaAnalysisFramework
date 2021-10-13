#ifndef ANITA_ANALYSIS_INTERPOLATED_NOTCH_H
#define ANITA_ANALYSIS_INTERPOLATED_NOTCH_H
?
#include "FilterOperation.h"
#include "TString.h"
#include "AnalysisWaveform.h"
?
class interpolatedNotchFilter : public UniformFilterOperation
{
?
 public:
?
  interpolatedNotchFilter();
?
  const char * tag () const { return "interpolatedNotch"; }
  const char * description () const { return "interpolatedNotchFilter"; }
?
  virtual void processOne(AnalysisWaveform *, const RawAnitaHeader *,int unused , int unsused2) ;
  virtual void process(FilteredAnitaEvent* ev);
  virtual void cutParams(bool notchLower,Double_t lowBinCenter,Double_t lowBinSigma,bool notchHigher,
			 Double_t highBinCenter,Double_t highBinSigma);
  TGraph * interpolatedFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq);
  
 private:
  Double_t lowBin,highBin,lowSigma,highSigma;
  bool notchLower,notchHigher;
?
  void init();
};
?
#endif 
