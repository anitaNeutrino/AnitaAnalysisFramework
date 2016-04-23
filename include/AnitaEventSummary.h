#ifndef _ANITA_EVENT_SUMMARY_H_
#define _ANITA_EVENT_SUMMARY_H_

#include "TObject.h" 
#include "Adu5Pat.h" 
#include "RawAnitaHeader.h" 
#include <iostream>


class AnitaEventSummary 
{
public: 
  static const Int_t maxDirectionsPerPol = 3; 

  class PointingHypothesis 
  {
  public: 
    Double_t phi;  // peak phi, degrees
    Double_t theta; // peak theta, degrees
    Double_t value; // peak value
    Double_t snr; // snr of peak
    Double_t hwAngle; // angle with respect to triggering phi sector
    Bool_t triggered; 
    Bool_t masked; 

    Double_t latitude;// on continent, or -9999 if doesn't intersect
    Double_t longitude;// on continent, or -9999 if doesn't intersect
    Double_t altitude;// on continent, or -9999 if doesn't intersect
    
    Double_t sigma_theta;  
    Double_t sigma_phi; 
    Double_t rho; 

    virtual ~PointingHypothesis() {;}
    ClassDef(PointingHypothesis,4); 
  }; 

  class WaveformInfo
  {
  public: 
    Double_t snr; 
    Double_t peakHilbert; 
    Double_t peakVal; 
    Double_t bandwidth;
    Int_t numAntennasInCoherent;     
    virtual ~WaveformInfo() {; } 

    ClassDef(WaveformInfo, 1); 
  }; 

  class EventFlags
  {
  public: 
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


  Int_t run; 
  UInt_t eventNumber; 
  
  Int_t nPeaks[AnitaPol::kNotAPol]; //Number of peaks actually found; this might be less than maxDirectionsPerPol 
  PointingHypothesis peak[AnitaPol::kNotAPol][maxDirectionsPerPol]; 
  WaveformInfo coherent[AnitaPol::kNotAPol][maxDirectionsPerPol];
  WaveformInfo deconvolved[AnitaPol::kNotAPol][maxDirectionsPerPol];


  // WaveformInfo maxWaveform;  // do we want this for all ? 

  EventFlags flags; 
  // Adu5Pat gps;
  // RawAnitaHeader header;

  AnitaEventSummary();
  AnitaEventSummary(const RawAnitaHeader* header);//, const Adu5Pat* pat);
  void setTriggerInfomation(const RawAnitaHeader* header);
  void zeroInternals();
  virtual ~AnitaEventSummary() { ; } 

  
private: 

  ClassDef(AnitaEventSummary, 5); 
}; 





#endif 
