#ifndef _ANITA_EVENT_SUMMARY_H_
#define _ANITA_EVENT_SUMMARY_H_

#include "TObject.h" 
#include "Adu5Pat.h" 
#include "RawAnitaHeader.h" 



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

    Double_t latitude;
    Double_t longitude;
    Double_t altitude;
    
    /* covariance matrix?*/ 
    //    Double_t var_theta;  
    //    Double_t var_phi; 
    //    Double_t covar; 

    virtual ~PointingHypothesis() {;}
    ClassDef(PointingHypothesis,2); 
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
    Int_t isSoftwareTrigger; 
    Int_t isPayloadBlast; 
    Int_t nadirFlag; 
    Int_t strongCWFlag;
    Int_t isHPolTrigger;
    Int_t isVPolTrigger;

    CalPulser pulser; 
    Bool_t isVarner; 
    Bool_t isVarner2; 
    virtual ~EventFlags() {;}

    ClassDef(EventFlags,2); 
  }; 


  Int_t run; 
  UInt_t eventNumber; 
  
  PointingHypothesis peak[AnitaPol::kNotAPol][maxDirectionsPerPol]; 
  WaveformInfo coherent[AnitaPol::kNotAPol][maxDirectionsPerPol];
  WaveformInfo deconvolved[AnitaPol::kNotAPol][maxDirectionsPerPol];
  // WaveformInfo maxWaveform;  // do we want this for all ? 

  EventFlags flags; 
  // Adu5Pat gps;
  // RawAnitaHeader header;


  AnitaEventSummary();
  AnitaEventSummary(const RawAnitaHeader* header);//, const Adu5Pat* pat);
  void tagTriggerAsHPolOrVPol(const RawAnitaHeader* header);
  void zeroInternals();
  virtual ~AnitaEventSummary() { ; } 

  
private: 

  ClassDef(AnitaEventSummary, 3); 
}; 





#endif 
