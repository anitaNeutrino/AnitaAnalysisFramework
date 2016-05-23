#ifndef _ANITA_EVENT_SUMMARY_H_
#define _ANITA_EVENT_SUMMARY_H_

#include "TObject.h" 
#include "Adu5Pat.h"
#include "UsefulAdu5Pat.h" 
#include "RawAnitaHeader.h"
#include <iostream>

/** Common analysis output format */ 


class AnitaEventSummary 
{
public: 

  /** The maximum number of hypotheses storable per polarization */ 
  static const Int_t maxDirectionsPerPol = 5; 

  /** A pointing hypothesis stores the results of interferometric pointing */ 
  class PointingHypothesis 
  {
  public: 
    Double_t phi;  /// peak phi, degrees
    Double_t theta; /// peak theta, degrees
    Double_t value; /// peak value
    Double_t snr; /// snr of peak
    Double_t hwAngle; /// angle with respect to triggering phi sector

    Double_t latitude;/// on continent, or -9999 if doesn't intersect
    Double_t longitude;/// on continent, or -9999 if doesn't intersect
    Double_t altitude;/// on continent, or -9999 if doesn't intersect
    
    Double_t sigma_theta;  ///error on theta
    Double_t sigma_phi;  /// error on phi
    Double_t rho;  ///correlation coefficient between theta and phi

    Bool_t triggered; /// was this in a triggered phi sector? 
    Bool_t masked; /// was this in a masked phi sector? 

    virtual ~PointingHypothesis() {;}
    ClassDef(PointingHypothesis,5); 
  }; 

  /** Stores information about a waveform (coherent or deconvolve) */ 
  class WaveformInfo
  {

  public: 
    Double_t snr; ///Signal to Noise of waveform 
    Double_t peakHilbert; /// peak of hilbert envelope
    Double_t peakVal;  /// peak value
    Double_t xPolPeakVal;  // Peak of xpol trace
    Double_t xPolPeakHilbert;  // Peak of xpol hilbert Envelope
    Double_t I,Q,U,V;  // Stokes Parameters
    Double_t bandwidth;  /// bandwidth of power spectrum 
    Int_t numAntennasInCoherent; // number of antennas used to make this 
    virtual ~WaveformInfo() {; } 

    ClassDef(WaveformInfo, 2); 
  }; 


  /** Stores various event flags */
  class EventFlags
  {
  public: 
    /** Is this event from a cal pulser? */ 
    enum CalPulser 
      {
        NONE, 
        WAIS, 
        LDB, 
        SIPLE,
        TD
      }; 

    Int_t isGood;
    Int_t isRF;
    Int_t isAdu5Trigger;
    Int_t isG12Trigger;
    Int_t isSoftwareTrigger;
    Int_t isMinBiasTrigger;
    Int_t isPayloadBlast;
    Int_t nadirFlag;
    Int_t strongCWFlag;
    Int_t isHPolTrigger;
    Int_t isVPolTrigger;

    CalPulser pulser;
    Bool_t isVarner;
    Bool_t isVarner2;
    virtual ~EventFlags() {;}

    ClassDef(EventFlags,3); 
  };

  /** A Source Hypothesis tells us about different potential sources of signals (e.g. calibration pulser) */ 
  class SourceHypothesis
  {
    public:
      Double_t theta;
      Double_t phi;
      Double_t distance;
      
      virtual ~SourceHypothesis(){;}

      ClassDef(SourceHypothesis,1);
  };

 

 
  Int_t run;
  UInt_t eventNumber;
  
  
  Int_t nPeaks[AnitaPol::kNotAPol]; ///Number of peaks actually found; this might be less than maxDirectionsPerPol 

  PointingHypothesis peak[AnitaPol::kNotAPol][maxDirectionsPerPol]; 
  WaveformInfo coherent[AnitaPol::kNotAPol][maxDirectionsPerPol]; 
  WaveformInfo deconvolved[AnitaPol::kNotAPol][maxDirectionsPerPol];


  // WaveformInfo maxWaveform;  // do we want this for all ? 

  EventFlags flags;

  SourceHypothesis sun;
  SourceHypothesis wais;
  SourceHypothesis ldb;

  // Adu5Pat pat;
  // RawAnitaHeader header;

  AnitaEventSummary();
  AnitaEventSummary(const RawAnitaHeader* header);
  AnitaEventSummary(const RawAnitaHeader* header, UsefulAdu5Pat* pat);  
  void setTriggerInfomation(const RawAnitaHeader* header);
  void setSourceInformation(UsefulAdu5Pat* pat);  
  void zeroInternals();
  virtual ~AnitaEventSummary() { ; } 

  
  private: 

    ClassDef(AnitaEventSummary, 9); 
}; 





#endif 
