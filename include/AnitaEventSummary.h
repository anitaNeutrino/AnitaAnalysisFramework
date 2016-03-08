#ifndef _ANITA_EVENT_SUMMARY_H_
#define _ANITA_EVENT_SUMMARY_H_

#include "TObject.h" 
#include "Adu5Pat.h" 
#include "RawAnitaHeader.h" 



class AnitaEventSummary 
{
  public: 
    static const UInt_t MaxPointingHypotheses = 5; 

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

        /* covariance matrix?*/ 
    //    Double_t var_theta;  
    //    Double_t var_phi; 
    //    Double_t covar; 

        virtual ~PointingHypothesis() {;}
      ClassDef(PointingHypothesis,1); 
    }; 

  class WaveformInfo
  {
    public: 
      Double_t snr; 
      Double_t peakHilbert; 
      Double_t peakVal; 
      Double_t bandwidth; 
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

        CalPulser pulser; 
        Bool_t isVarner; 
        Bool_t isVarner2; 
        virtual ~EventFlags() {;}

      ClassDef(EventFlags,1); 
    }; 


    UInt_t runNumber; 
    UInt_t eventNumber; 

    PointingHypothesis hypotheses[MaxPointingHypotheses]; 
    UInt_t nantennas_in_coherent; 
    WaveformInfo coherent[MaxPointingHypotheses]; 
    WaveformInfo deconvolved[MaxPointingHypotheses]; 
    WaveformInfo maxWaveform;  // do we want this for all ? 

    EventFlags flags; 
    Adu5Pat gps; 
    RawAnitaHeader header; 

    virtual ~AnitaEventSummary() { ; } 

  private: 

    ClassDef(AnitaEventSummary, 1); 
}; 





#endif 
