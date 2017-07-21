

#ifndef PRETTYANALYSISWAVEFORM_H
#define PRETTYANALYSISWAVEFORM_H
#include "UsefulAnitaEvent.h"
#include "CorrelationSummary.h"
#include "CorrelationSummaryAnita4.h"
#include "AnalysisWaveform.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <string>

#include <iostream>
#include <fstream>


//ANITA Includes
#include "CalibratedAnitaEvent.h"
#include "FilteredAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "FFTtools.h"
#include "FFTWComplex.h"


//ROOT Includes
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TVirtualFFT.h"
#include "TMinuit.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"





class PrettyAnalysisWaveform: public FilteredAnitaEvent
{


 public:

   PrettyAnalysisWaveform(UsefulAnitaEvent* uae, FilterStrategy* strategy, Adu5Pat* pat, RawAnitaHeader* header);


  //Putative Analysis methods
  //! Calls getMaxAntennaCorrelation
  /*!
    \param pol Which polarisation to use?
    \param peakPtr An optional pointer to a double in which to store the peak/rms value.
    \return The antenna number of the antenna with the largest signal.
  */
  int getMaxAntenna(AnitaPol::AnitaPol_t pol, Double_t *peakPtr=0);
  //! Select the antenna with the maximum voltage squared.
  /*!
    \param pol Which polarisation to use?
    \param peakPtr An optional pointer to a double in which to store the peak V^2 value.
    \return The antenna number of the antenna with the largest V^2.
  */
  int getMaxAntennaVSquared(AnitaPol::AnitaPol_t pol, Double_t *peakPtr=0);
  //! Select the upper antenna with the maximum correlation (defined as peak/rms of the correlation) with it's pair in the lower ring.
  /*!
    \param pol Which polarisation to use?
    \param peakPtr An optional pointer to a double in which to store the peak/rms correlation value.
    \return The antenna number of the antenna with the largest peak/rms correlation value.
  */  
  int getMaxAntennaCorrelation(AnitaPol::AnitaPol_t pol, Double_t *peakPtr=0);
  //! Select the upper antenna with the maximum correlation (defined as peak/rms of the correlation) with it's pair in the lower ring.
  /*!
    \param pol Which polarisation to use?
    \param peakPtr An optional pointer to a double in which to store the peak/rms correlation value.
    \return The antenna number of the antenna with the largest peak/rms correlation value.
  */  
  int getMaxAntennaCorrelationRollingAvg(AnitaPol::AnitaPol_t pol, Double_t *peakPtr=0);

  //! Generates a CorrelationSummaryAnita3 object for a set of 15 antennas.
  /*!
    \param centreAnt The number of one of the antennas in the centre of the set of 10.
    \param pol Which polarisation to use?
    \param deltaT An optional value to use if interpolation is required. This value is taking as being the desired sampling period of the interpolated waveforms.
    \return A pointer to the CorrelationSummaryAnita3 object that is created.
  */  
  CorrelationSummaryAnita4 *getCorrelationSummaryAnita4(Int_t centreAnt,AnitaPol::AnitaPol_t pol,Int_t deltaT, Int_t eventNumber);


  AnalysisWaveform *getCorrelation(int chanIndex1, int chanIndex2); ///< Wrapper around FFTtools::getCorrelationGraph
  AnalysisWaveform *getCorrelationInterpolated(int chanIndex1, int chanIndex2, AnitaPol::AnitaPol_t pol, Int_t npadfreq=8 ); ///< Wrapper around FFTtools::getInterpolatedCorrelationGraph
  
  void fillSixAntArrays(int ant, int topAnts[3], int bottomAnts[3]); ///< Utility to get neighbouring antenna numbers
  void fillNextFourAntArrays(int ant, int nextFourAnts[4]);///< Utility to get next to neighbouring antenna numbers
  void fillNadirArrays(int ant, int nadirAnts[9]);

  void fillNineAntArrays(int ant, int nineAnts[9]); ///< Utility to get neighbouring antenna numbers ( Top 0-2, Middle 3-5, Bottom 6-8)
  void fillNextSixAntArrays(int ant, int nextFourAnts[4]);///< Utility to get next to neighbouring antenna numbers
	double getHighestSnr(int centreAntenna);


#ifndef PLEASE_IGNORE_ME_DOXYGEN
  ClassDef(PrettyAnalysisWaveform,1); ///< ROOT's magic macro
#endif

  void setPassBandFilterFlag( int flag) { fPassBandFilter=flag;}
  void setNotchFilterFlag( int numNotches) { fNotchFilter=numNotches;}
  void setPassBandLimits(Double_t low, Double_t high)
     { fLowPassEdge=low; fHighPassEdge=high;}
  void setNotchBandLimits(Int_t notchNum, Double_t low, Double_t high)
     { fLowNotchEdge[notchNum]=low; fHighNotchEdge[notchNum]=high;}


 private:

  Double_t fDeltaT; ///< The interpolated sampling rate.
  Double_t fWaveOffset; ///< The difference in T0 of two channels.
  Int_t fPassBandFilter; ///< Whether or not to pass band filter the interpolated waves;
  Int_t fNotchFilter; ///< Whether or not to notch filter, and how many notches
  Double_t fLowPassEdge; ///< The lower edge of the pass band
  Double_t fHighPassEdge; ///< The higher edge of the pass band
  Double_t fLowNotchEdge[10]; ///< The lower edge of the notch band
  Double_t fHighNotchEdge[10]; ///< The higher edge of the notch band
  

};


#endif //PRETTYANALYSISWAVEFORM_H
