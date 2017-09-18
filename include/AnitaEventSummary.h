#ifndef _ANITA_EVENT_SUMMARY_H_
#define _ANITA_EVENT_SUMMARY_H_

#include "TObject.h" 
#include "AnitaConventions.h"
#include <iostream>

class Adu5Pat;
class UsefulAdu5Pat;
class RawAnitaHeader; 
class TruthAnitaEvent; 

/** 
 * @class AnitaEventSummary
 * @brief Common analysis format between UCorrelator and Acclaim
 * 
 * Two independent analyses will fill most of the variables in these trees.
 * Needless to say, there's no guarantee that everything will be filled,
 * so be wary if something is 0 (it may not have been filled).
 * Also, member variables that have been filled may have slightly different
 * definitions between analyses. You have been warned.
 *
 * This class and its subclasses have utility functions that can be used inside TTree::Draw().
 * For example, with a TTree (called sumTree) of AnitaEventSummaries (called sum) doing
 *
 * sumTree->Draw("sum.highestPeak().dPhiWais()")
 *
 * Should produce a histogram of the reconstructed phi-angle from the WAIS divide cal pulser.
 */

class AnitaEventSummary : public TObject
{
public: 


  //------------------------------------------------------------------------------------
  /*************************************************************************************
   * Static members (for array sizes inside class)
   *************************************************************************************/
  static const Int_t maxDirectionsPerPol = 5; /// The maximum number of hypotheses storable per polarization */ 
  static const Int_t peaksPerSpectrum = 3; /// The maximum number of frequency peaks per waveform spectrum



  //------------------------------------------------------------------------------------
  /*************************************************************************************
   * Sub-classes
   *************************************************************************************/


  /** 
   * @class PointingHypothesis
   * Stores some summary information about the peak of an interferometric map.
   */
  class SourceHypothesis; // Forward class declarations for PointingHypothesis utils
  class PayloadLocation; // Forward class declarations for PointingHypothesis utils

  class PointingHypothesis 
  {
  public: 
    PointingHypothesis() : fContainer(NULL) { ; }
    Double32_t phi;  /// peak phi, degrees
    Double32_t theta; /// peak theta, degrees
    Double32_t value; /// peak value
    Double32_t snr; /// snr of peak
    Double32_t mapRMS; /// rms of interferometric map
    Double32_t mapHistoryVal; /// value of average of the peak location over the past 60 min-bias events
    Double32_t hwAngle; /// angle with respect to triggering phi sector
    Double32_t hwAngleXPol; /// angle with respect to triggering phi sector in opposite polarisation
    Double32_t latitude;/// on continent, or -9999 if doesn't intersect
    Double32_t longitude;/// on continent, or -9999 if doesn't intersect
    Double32_t altitude;/// on continent, or -9999 if doesn't intersect
    Double32_t distanceToSource; /// on continent, or -9999 if doesn't intersect
    Double32_t sigma_theta;  /// error on theta
    Double32_t sigma_phi;  /// error on phi
    Double32_t rho;  /// correlation coefficient between theta and phi
    Double32_t chisq; /// chisq/ndof of peak finding process, if available (otherwise zero)
    Double32_t theta_adjustment_needed; /// If an event barely missed the ground, it is useful to see the coordinates at which it would hit if theta adjustment by a small amount. This is the calculated small amount that leads to it hitting the ground. 
    Double32_t phi_separation; /// angular separation from higher value peak in same event. 1000 if highest value event (i.e. first hypothesis) 
    Double32_t dphi_rough;  /// phi - phi rough
    Double32_t dtheta_rough; /// theta - theta rough
    Bool_t triggered; /// was this in a triggered phi sector? 
    Bool_t triggered_xpol; /// was this in a triggered xpol phi sector?
    Bool_t masked; /// was this in a masked phi sector?
    Bool_t masked_xpol; /// was this in a masked phi xpol sector?

    // Resolution utility functions
    double dPhi(double phi) const;
    double dTheta(double theta, bool different_sign_conventions = false) const;
    double bearing() const;
    double dPhiWais() const;
    double dThetaWais() const;
    double dPhiSun() const;
    double dThetaSun() const;
    double dPhiLDB() const;
    double dThetaLDB() const;
    double dPhiMC() const;
    double dThetaMC() const;
    double minAbsHwAngle() const;
    Bool_t absHwAngleLessThanAbsHwAngleXPol() const;

   private:
    //----------------------------------------------------------------------------------------------------
    // WARNING! This private nonsense is fragile, and a bit hacky.
    // The //! comment after the fContainer member means it does not persist in ROOT.
    // Please do not edit that comment.
    //----------------------------------------------------------------------------------------------------
    mutable AnitaEventSummary* fContainer; //! WARNING! Does not persist! Get access to AnitaEventSummary that contains this PointingHypothesis
    const AnitaEventSummary* getContainer(const char* funcName) const; /// Wraps getting fContainer with a warning if NULL.
    double dPhiSource(const SourceHypothesis& source) const;   // Won't work inside TTree::Draw due to limitations in TTreeFormula, so are private, use e.g. dPhiWais() instead.
    double dThetaSource(const SourceHypothesis& source) const; // Won't work inside TTree::Draw due to limitations in TTreeFormula, so are private, use e.g. dThetaWais() instead.
    void printEvent() const;
    friend class AnitaEventSummary;

    ClassDefNV(PointingHypothesis,14);
  }; 

  /** 
   * @class WaveformInfo
   * @brief Stores information about a coherently summed waveform (filtered/unfiltered/deconvolved)
   * The coherent summing of the waveform corresponds to a direction stored in a PointingHypothesis
   */
  class WaveformInfo
  {

   public: 
    WaveformInfo() {; } 
    Double32_t snr; /// Signal to Noise of waveform 
    Double32_t peakHilbert; /// peak of hilbert envelope
    Double32_t peakVal;  /// peak value
    Double32_t xPolPeakVal;  /// Peak of xpol trace
    Double32_t xPolPeakHilbert;  /// Peak of xpol hilbert Envelope

    Double32_t I,Q,U,V;  // Stokes Parameters


    Double32_t totalPower;  /// Total power in waveform
    Double32_t totalPowerXpol;  /// Total power in xPol waveform

    //spectrum info 
    Double32_t bandwidth[peaksPerSpectrum];  /// bandwidth of each peak (implementation defined, may not be comparable between analyses) 
    Double32_t peakFrequency[peaksPerSpectrum]; /// peak frequency of power spectrum 
    Double32_t peakPower[peaksPerSpectrum]; //power within +/- bandwidth of each peak 
    Double32_t spectrumSlope;  ///  Slope of line fit to spectrum (in log-space, so this is spectral-index) 
    Double32_t spectrumIntercept; /// Intercept of line fit to spectrum (in log-space) 

    //Shape parameters, computed using hilbert envelope 
    // This should probably taken out into its own class 
    Double32_t riseTime_10_90;  /// Rise time of hilbert env from 10% to 90% of peak
    Double32_t riseTime_10_50;  /// Rise time of hilbert env from 10% to 50% of peak
    Double32_t fallTime_90_10;  /// Fall time of hilbert env from 90% to 10% of peak
    Double32_t fallTime_50_10;  /// Fall time of hilbert env from 50% to 10% of peak
    Double32_t width_50_50;   /// Width from first envelope crossing of 50 percent of peak to last 
    Double32_t width_10_10;  /// Width from first envelope crossing of 10 percent of peak to last 
    Double32_t power_10_10;  /// Power enclosed within 10_10 width
    Double32_t power_50_50;  /// Power enclosed within 50_50 width
    Double32_t peakTime;  // Time that peak hilbert env occurs
    Double32_t peakMoments[5];  // moments about Peak  (1st - 5th moments) 

    Double32_t impulsivityMeasure; /// A number that has something to do with how impulsive it is
    Double32_t narrowestWidths[5]; /// Narrowest width containing {10%, 20%, 30%, 40%, 50%} of the total power

    Int_t numAntennasInCoherent; /// Number of antennas used to make this

    Double32_t localMaxToMin; /// Largest value of local max to neighbouring local min (see Acclaim::RootTools::getLocalMaxToMin)
    Double32_t localMaxToMinTime; /// Time between local maxima and minima +ve means max is before min, -ve means min is before max
    Double32_t globalMaxToMin; /// Difference between maximum and minimum voltage
    Double32_t globalMaxToMinTime; /// Time between maximum and minimum volts, +ve means max is before min, -ve means min is before max 


    //some utilities for polarization info
    double linearPolFrac() const;
    double linearPolAngle() const;
    double circPolFrac() const;
    double totalPolFrac() const;

    double standardizedPeakMoment(int i) const;

    inline double skewness(){return standardizedPeakMoment(3);}
    inline double kurtosis(){return standardizedPeakMoment(4);}

    ClassDefNV(WaveformInfo, 11);

    // private:
    //  mutable AnitaEventSummary* fContainer; //! Disgusting hack
  }; 

  /** 
   * @class ChannelInfo
   * @Stores brief information of a channel's waveform
   */
  class ChannelInfo
  {

   public:
    /// Correct indices are set in the AnitaEventSummary constructor
    ChannelInfo() : pol(AnitaPol::kNotAPol), ant(-1) {; }

    Double32_t rms;
    Double32_t avgPower;
    Double32_t snr; /// Signal to Noise of waveform 
    Double32_t peakHilbert; /// peak of hilbert envelope

    Double_t getPhi() const;
    inline Int_t getAnt() const {return ant;}                // could add some errors on -1 here...
    inline AnitaPol::AnitaPol_t getPol() const {return pol;} // could add some errors on 2 here...

   private:
    friend class AnitaEventSummary;
    // WARNING! THE COMMENTS TRAILING THESE AFFECT WHAT HAPPENS IN ROOT!
    // The //! comment initializer means ROOT does not store these variables when writing to files.
    // These variables are filled in the AnitaEventSummary constructor, so each ChannelInfo knows
    // its location in the ChannelInfo channels[2][48] array in AnitaEventSummary
    AnitaPol::AnitaPol_t     pol; //! DOES NOT PERSIST IN ROOT! the polarization
    Int_t                    ant; //! DOES NOT PERSIST IN ROOT! the antenna

    ClassDefNV(ChannelInfo, 4);
  }; 


  /** 
   * @class EventFlags
   * Stores simple integer numbers based on quality cuts, calibration pulser timing, trigger type, etc.
   */
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
    Int_t isStepFunction;

    CalPulser pulser;
    Bool_t isVarner;
    Bool_t isVarner2;

    /** These are used to cut out payload blasts and stuf like that. 
     *  The first element is the total, and then the next are by ring 
     *  So to get the top ring, do 1 + AnitaRing::kTopRing, etc. 
     */
    Double32_t meanPower[1+AnitaRing::kNotARing];
    Double32_t medianPower[1+AnitaRing::kNotARing];
    Double32_t meanPowerFiltered[1+AnitaRing::kNotARing];
    Double32_t medianPowerFiltered[1+AnitaRing::kNotARing];

    Double32_t maxBottomToTopRatio[AnitaPol::kNotAPol]; 
    int maxBottomToTopRatioSector[AnitaPol::kNotAPol]; 
    Double32_t minBottomToTopRatio[AnitaPol::kNotAPol]; 
    int minBottomToTopRatioSector[AnitaPol::kNotAPol]; 


    Int_t nSectorsWhereBottomExceedsTop; 

    /** The fraction of nearby events that are payload blasts */ 
    Double_t blastFraction; 

    ClassDefNV(EventFlags,10);
  };





  /** 
   * @class SourceHypothesis
   * For known sources, such as a calibration pulser
   */
  class SourceHypothesis
  {
   public:
    SourceHypothesis() { reset(); }
    Double32_t theta;
    Double32_t phi;
    Double32_t distance;

    Double32_t mapValue[NUM_POLS];  ///what the instantaneous map value is at this source hypothesis
    Double32_t mapHistoryVal[NUM_POLS]; /// a history of the interferometric map value for the source location

    void reset(); /// sets all the values to nonsense.  Sorry, mapHistoryVal means this is in source now 
      

    ClassDefNV(SourceHypothesis,3);
  };




  /** 
   * @class MCTruth
   * Summary information from the Monte Carlo, true neutrino direction, weight, and energy.
   */
  class MCTruth : public SourceHypothesis
  {
   public:
    MCTruth() { reset(); }
    WaveformInfo wf[AnitaPol::kNotAPol]; 
    double weight;
    double energy;
    void reset(); 

    ClassDefNV(MCTruth,5);
  }; 



  
  /** 
   * @class PayloadLocation
   * Contains the most important 4 numbers in the GPS data.
   */
  class PayloadLocation
  {
   public:
    PayloadLocation() { reset(); }
    PayloadLocation(const Adu5Pat* pat); /// Slightly more useful constructor

    Float_t latitude;
    Float_t longitude;
    Float_t altitude;
    Float_t heading;

    Float_t prevHeading; //useful for determining rotation rate

    void reset() { latitude = -999; longitude = -999; altitude = -999; heading = -999; prevHeading = -999;};
    void update(const Adu5Pat* pat); /// Copy the data from the pat into the object

    ClassDefNV(PayloadLocation,2);
  };  






  //------------------------------------------------------------------------------------
  /*************************************************************************************
   * Public member variables
   *************************************************************************************/
  Int_t run; /// Run
  UInt_t eventNumber; /// Event number
  UInt_t realTime; /// Time of the event
  Int_t nPeaks[AnitaPol::kNotAPol]; /// Number of peaks stored in this AnitaEventSummary (might be less than maxDirectionsPerPol)
  PointingHypothesis peak[AnitaPol::kNotAPol][maxDirectionsPerPol]; /// Summaries of the event peak directions (indices of all WaveformInfo member arrays match peak index)
  WaveformInfo coherent[AnitaPol::kNotAPol][maxDirectionsPerPol]; /// Summaries of the (unfiltered) coherently summed waveforms, array index correponds to entry in peak[][] 
  WaveformInfo deconvolved[AnitaPol::kNotAPol][maxDirectionsPerPol]; /// Summaries of the (unfiltered) de-dispersed coherently summed waveforms, array index correponds to entry in peak[][] 
  WaveformInfo coherent_filtered[AnitaPol::kNotAPol][maxDirectionsPerPol]; /// Summaries of the filtered, coherently summed waveforms, array index correponds to entry in peak[][] 
  WaveformInfo deconvolved_filtered[AnitaPol::kNotAPol][maxDirectionsPerPol]; /// Summaries of the filtered, de-dispersed, coherently summed waveforms, array index correponds to entry in peak[][] 
  ChannelInfo channels[AnitaPol::kNotAPol][NUM_SEAVEYS]; /// Summaries of each channel's waveform.
  EventFlags flags; /// Flags corresponding the event quality, trigger type, calibration pulser timing, etc.
  SourceHypothesis sun; /// Contains location of sun in map coordinates at time of event
  SourceHypothesis wais; /// Contains location of WAIS divide cal pulser in map coordinates at time of event
  SourceHypothesis ldb; /// Contains location of LDB cal pulser in map coordinates at time of event
  MCTruth mc; /// Contains summary information about MC truth, if real data then this filled with constant, unphysical values.
  PayloadLocation anitaLocation; /// Reduced GPS data


  //------------------------------------------------------------------------------------
  /*************************************************************************************
   * Public member functions
   *************************************************************************************/
  // See source file for doxygen comments
  AnitaEventSummary();
  AnitaEventSummary(const RawAnitaHeader* header);
  AnitaEventSummary(const RawAnitaHeader* header, UsefulAdu5Pat* pat, const TruthAnitaEvent * truth = 0);
  void setTriggerInfomation(const RawAnitaHeader* header);
  void setSourceInformation(UsefulAdu5Pat* pat, const TruthAnitaEvent * truth = 0);  
  void zeroInternals();


  // Utilities to find interesting entries in AnitaEventSummary
  AnitaPol::AnitaPol_t highestPol() const;
  Int_t highestPolAsInt() const;
  Int_t highestPeakInd() const;
  const PointingHypothesis& highestPeak() const;
  const WaveformInfo& highestCoherent() const;
  const WaveformInfo& highestDeconvolved() const;
  const WaveformInfo& highestCoherentFiltered() const;
  const WaveformInfo& highestDeconvolvedFiltered() const;

  inline double weight(){return mc.weight > 0 ? mc.weight : 1;} /// Return the weight of the event, always returns 1 for data, the weight from MCTruth otherwise
  AnitaPol::AnitaPol_t mcPol() const;
  int mcPolAsInt() const;
  int mcPeakInd() const;
  const PointingHypothesis& mcPeak() const;
  const WaveformInfo& mcCoherent() const;
  const WaveformInfo& mcDeconvolved() const;
  const WaveformInfo& mcCoherentFiltered() const;
  const WaveformInfo& mcDeconvolvedFiltered() const;

  //------------------------------------------------------------------------------------
 private:

  /*************************************************************************************
   * Private member variables/functions
   *************************************************************************************/ 
  // WARNING! THE COMMENTS TRAILING THESE AFFECT WHAT HAPPENS IN ROOT!
  // The //! comment initializer means ROOT does not store these variables when writing to files.
  // We want to keep it this way, since they are only used to cache the results of the utility
  // funtions to find the highest peak and peak nearest the monte carlo truth info.
  mutable Int_t                fHighestPeakIndex; //! DOES NOT PERSIST IN ROOT! Internal index to cache result of finding highest peak
  mutable AnitaPol::AnitaPol_t fHighestPol;       //! DOES NOT PERSIST IN ROOT! Internal index to cache result of finding highest peak
  mutable Int_t                fMCPeakIndex;      //! DOES NOT PERSIST IN ROOT! Internal index to cache result of finding peak nearest MC
  mutable AnitaPol::AnitaPol_t fMCPol;            //! DOES NOT PERSIST IN ROOT! Internal index to cache result of finding peak nearest MC
  mutable UInt_t               fLastEventNumber;  //! DOES NOT PERSIST IN ROOT! To check for stale caching variables
  /** 
   * Workhorse function to find the highest peak
   * Caches the result in the mutable, non-ROOT-persistent members fHighestPol and fHighestPeakIndex
   */
  void findHighestPeak() const;
  void findMC() const;
  void resetNonPersistent() const;



  ClassDefNV(AnitaEventSummary, 29);
};

#endif 
