#ifndef _ANITA_EVENT_SUMMARY_H_
#define _ANITA_EVENT_SUMMARY_H_

#include "TObject.h" 
#include "AnitaConventions.h" 
class UsefulAdu5Pat; 
class RawAnitaHeader; 

/** Common analysis output format 
 *  
 *  Needless to say, there's no guarantee that everything will be filled, so be wary if something is 0 (it may not have been filled).   
 *
 * */ 


class AnitaEventSummary 
{
public: 

  /** The maximum number of hypotheses storable per polarization */ 
  static const Int_t maxDirectionsPerPol = 5; 

  /** A pointing hypothesis stores the results of interferometric pointing */ 
  class PointingHypothesis 
  {
  public: 
    PointingHypothesis() { ; }
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
    Double_t chisq; /// chisq/ndof of peak finding process, if available (otherwise zero)


    Double_t theta_adjustment_needed; /// If an event barely missed the ground, it is useful to see the coordinates at which it would hit if theta adjustment by a small amount. This is the calculated small amount that leads to it hitting the ground. 
    Double_t phi_separation; //angular separation from higher value peak in same event. 1000 if highest value event (i.e. first hypothesis) 

    Bool_t triggered; /// was this in a triggered phi sector? 
    Bool_t masked; /// was this in a masked phi sector? 

    ClassDefNV(PointingHypothesis,9); 
  }; 

  /** Stores information about a waveform (coherent or deconvolve) */ 
  class WaveformInfo
  {

  public: 
    WaveformInfo() {; } 
    Double_t snr; ///Signal to Noise of waveform 
    Double_t peakHilbert; /// peak of hilbert envelope
    Double_t peakVal;  /// peak value
    Double_t xPolPeakVal;  // Peak of xpol trace
    Double_t xPolPeakHilbert;  // Peak of xpol hilbert Envelope
    Double_t I,Q,U,V;  // Stokes Parameters
    Double_t bandwidth;  /// bandwidth of power spectrum 

    //Shape parameters, computed using hilbert envelope 
    Double_t riseTime_10_90;  /// Rise time of hilbert env from 10% to 90% of peak
    Double_t riseTime_10_50;  /// Rise time of hilbert env from 10% to 50% of peak
    Double_t fallTime_90_10;  /// Fall time of hilbert env from 90% to 10% of peak
    Double_t fallTime_50_10;  /// Fall time of hilbert env from 50% to 10% of peak
    Double_t width_50_50;   /// Width from first envelope crossing of 50 percent of peak to last 
    Double_t width_10_10;  /// Width from first envelope crossing of 10 percent of peak to last 
    Int_t numAntennasInCoherent; // number of antennas used to make this 

    ClassDefNV(WaveformInfo, 3); 
  }; 


  /** Stores various event flags */
  class EventFlags
  {
  public: 
    EventFlags() {; }
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

    ClassDefNV(EventFlags,3); 
  };

  /** A Source Hypothesis tells us about different potential sources of signals (e.g. calibration pulser) */ 
  class SourceHypothesis
  {
    public:
      SourceHypothesis() { reset(); }
      Double_t theta;
      Double_t phi;
      Double_t distance;

      void reset() { theta = -999; phi = -999; distance = -999; } 
      

      ClassDefNV(SourceHypothesis,1);
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


  AnitaEventSummary();
  AnitaEventSummary(const RawAnitaHeader* header);
  AnitaEventSummary(const RawAnitaHeader* header, UsefulAdu5Pat* pat);  
  void setTriggerInfomation(const RawAnitaHeader* header);
  void setSourceInformation(UsefulAdu5Pat* pat);  
  void zeroInternals();

  
  private: 

    ClassDefNV(AnitaEventSummary, 9); 
}; 





#endif 
