#ifndef _ANITA_EVENT_SUMMARY_H_
#define _ANITA_EVENT_SUMMARY_H_

#include "TObject.h" 
#include "Adu5Pat.h" 
#include "RawAnitaHeader.h" 



class AnitaEventSummary 
{
  public: 
    const UInt_t MaxPointingHypotheses = 5; 

    class PointingHypothesis 
    {
      public: 
        Double_t phi;  // peak phi
        Double_t theta; // peak theta
        Double_t value; // peak value

        /* covariance matrix ?*/ 
    //    Double_t var_theta;  
    //    Double_t var_phi; 
    //    Double_t covar; 

      ClassDef(PointingHypothesis,1); 
    }; 

   class WaveformInfo
  {
    public: 
      Double_t snr; 
      Double_t peakOfHilbert; 
      Double_t peakVal; 

      ClassDef(WaveformInfo, 1); 
  }; 

   class EventFlags
    {
      public: 
        enum CalPulser
        {
          NONE, 
          WAIS, 
          MCMURDO, 
          TD 
        }; 

        Int_t isGood; 
        Int_t isRF; 
        Int_t isSoftwareTrigger; 
        Int_t isPayloadBlast; 
        Int_t nadirFlag; 

        CalPulser pulser; 
        Bool_t isVarner; 
        Bool_t isVarner2; 

      ClassDef(EventFlags,1); 
    }; 


    UInt_t runNumber; 
    UInt_t eventNumber; 

    PointingHypothesis hypotheses[MaxPointingHypotheses]; 
    WaveformInfo coherent[MaxPointingHypotheses]; 
    WaveformInfo maxWaveform;  // do we want this for all ? 

    EventFlags flags; 
    Adu5Pat gps; 
    RawAnitaHeader header; 

  private: 

    ClassDef(AnitaEventSummary, 1); 
}; 





#endif 
