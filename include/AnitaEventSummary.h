#ifndef _ANITA_EVENT_SUMMARY_H_
#define _ANITA_EVENT_SUMMARY_H_

#include "TObject.h"
#include "AnitaConventions.h"
#include "Adu5Pat.h"
#include <iostream>

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
  static const Int_t numFracPowerWindows = 5;

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
    Double_t phi;  /// peak phi, degrees
    Double_t theta; /// peak theta, degrees
    Double_t value; /// peak value
    Double_t snr; /// snr of peak
    Double_t mapRMS; /// rms of interferometric map
    Double_t mapHistoryVal; /// value of average of the peak location over the past 60 min-bias events
    Double_t hwAngle; /// angle with respect to triggering phi sector
    Double_t hwAngleXPol; /// angle with respect to triggering phi sector in opposite polarisation
    Double_t latitude;/// on continent, or -9999 if doesn't intersect
    Double_t longitude;/// on continent, or -9999 if doesn't intersect
    Double_t altitude;/// on continent, or -9999 if doesn't intersect
    Double_t distanceToSource; /// on continent, or -9999 if doesn't intersect
    Double_t sigma_theta;  /// error on theta
    Double_t sigma_phi;  /// error on phi
    Double_t rho;  /// correlation coefficient between theta and phi
    Double_t chisq; /// chisq/ndof of peak finding process, if available (otherwise zero)
    Double_t theta_adjustment_needed; /// If an event barely missed the ground, it is useful to see the coordinates at which it would hit if theta adjustment by a small amount. This is the calculated small amount that leads to it hitting the ground.
    Double_t phi_separation; /// angular separation from higher value peak in same event. 1000 if highest value event (i.e. first hypothesis)
    Double_t dphi_rough;  /// phi - phi rough
    Double_t dtheta_rough; /// theta - theta rough
    Bool_t triggered; /// was this in a triggered phi sector?
    Bool_t triggered_xpol; /// was this in a triggered xpol phi sector?
    Bool_t masked; /// was this in a masked phi sector?
    Bool_t masked_xpol; /// was this in a masked phi xpol sector?
    Double_t antennaPeakAverage; /// the average of channel peaks in this direction 

    // Most basic resolution utility functions in payload coordinates relative to ADU5-aft-fore line
    Double_t dPhi(Double_t phi) const;
    Double_t dTheta(Double_t theta, bool different_sign_conventions = false) const;

    // Peak direction relative to north (N->E->S->W)
    Double_t bearing() const;
    Double_t dPhiNorth() const;

    // Find angle from stored source hypotheses
    Double_t dPhiWais() const;
    Double_t dThetaWais() const;
    Double_t dPhiSun() const;
    Double_t dThetaSun() const;
    Double_t dPhiLDB() const;
    Double_t dThetaLDB() const;
    Double_t dPhiMC() const;
    Double_t dThetaMC() const;
    Double_t dPhiTagged() const;    /// See AnitaEventSummary::sourceFromTag()
    Double_t dThetaTagged() const;  /// See AnitaEventSummary::sourceFromTag()

    // Are you within this theta/phi of a stored hypothesis?
    Bool_t closeToMC(Double_t deltaPhiDeg, Double_t deltaThetaDeg) const;
    Bool_t closeToWais(Double_t deltaPhiDeg, Double_t deltaThetaDeg) const;
    Bool_t closeToLDB(Double_t deltaPhiDeg, Double_t deltaThetaDeg) const;
    Bool_t closeToSun(Double_t deltaPhiDeg, Double_t deltaThetaDeg) const;
    Bool_t closeToTagged(Double_t deltaPhiDeg, Double_t deltaThetaDeg) const;

    // peak direction relative to trigger information
    Double_t minAbsHwAngle() const;
    Bool_t absHwAngleLessThanAbsHwAngleXPol() const;

   private:
    //----------------------------------------------------------------------------------------------------
    // WARNING! This private nonsense is fragile, and a bit hacky.
    // The //! comment after the fContainer member means it does not persist in ROOT.
    // Please do not edit that comment.
    //----------------------------------------------------------------------------------------------------
    mutable AnitaEventSummary* fContainer; //! WARNING! Does not persist! Get access to AnitaEventSummary that contains this PointingHypothesis
    const AnitaEventSummary* getContainer(const char* funcName) const; /// Wraps getting fContainer with a warning if NULL.
    Double_t dPhiSource(const SourceHypothesis& source) const;   // Won't work inside TTree::Draw due to limitations in TTreeFormula, so are private, use e.g. dPhiWais() instead.
    Double_t dThetaSource(const SourceHypothesis& source) const; // Won't work inside TTree::Draw due to limitations in TTreeFormula, so are private, use e.g. dThetaWais() instead.
    void printEvent() const;
    friend class AnitaEventSummary;

    ClassDefNV(PointingHypothesis,15);
  };

  /**
   * @class WaveformInfo
   * @brief Stores information about a coherently summed waveform (filtered/unfiltered/deconvolved)
   * The coherent summing of the waveform corresponds to a direction stored in a PointingHypothesis
   * Units are assumed to be in mV, ns,and GHz. 
   */
  class WaveformInfo
  {

   public:
    WaveformInfo() : fContainer(NULL), fLastEventNumberCache(0), nwMeanCache(-1),
                     nwGradCache(-1), nwInterceptCache(-1), nwChisquareCache(-1) {; }
    Double_t snr; //[0,100,16]  /// Signal to Noise of waveform
    Double_t peakHilbert;//[0,4096,21]  /// peak of hilbert envelope
    Double_t peakVal;  //[0,4096,21] /// peak value
    Double_t xPolPeakVal; //[0,4096,21]  /// Peak of xpol trace
    Double_t xPolPeakHilbert; //[0,4096,21]  /// Peak of xpol hilbert Envelope

    Double_t I,Q,U,V;  /// Integral Stokes Parameters (over the entire waveform) 
    Double_t max_dI,max_dQ,max_dU,max_dV; /// instantanteous stokes parameters (computed near max_dI).
    Int_t NPointsMaxStokes; /// The number of points used in the above estimates 

    Double_t totalPower;  /// Total power in waveform
    Double_t totalPowerXpol;  /// Total power in xPol waveform

    //spectrum info
    Double_t bandwidth[peaksPerSpectrum]; //[0,2,16] /// bandwidth of each peak (implementation defined, may not be comparable between analyses)
    Double_t peakFrequency[peaksPerSpectrum]; //[0,2,16] /// peak frequency of power spectrum
    Double_t peakPower[peaksPerSpectrum]; //power within +/- bandwidth of each peak
    Double_t spectrumSlope; // [-100,100,16] ///  Slope of line fit to spectrum (in log-space, so this is spectral-index)
    Double_t spectrumIntercept; // [-200,200,16] /// Intercept of line fit to spectrum (in log-space)

    //Shape parameters, computed using hilbert envelope
    // This should probably taken out into its own class
    Double_t riseTime_10_90; //[0,128,16]  /// Rise time of hilbert env from 10% to 90% of peak
    Double_t riseTime_10_50; //[0,128,16] /// Rise time of hilbert env from 10% to 50% of peak
    Double_t fallTime_90_10; //[0,128,16]/// Fall time of hilbert env from 90% to 10% of peak
    Double_t fallTime_50_10; //[0,128,16] /// Fall time of hilbert env from 50% to 10% of peak
    Double_t width_50_50;  //[0,128,16] /// Width from first envelope crossing of 50 percent of peak to last
    Double_t width_10_10;  //[0,128,16]/// Width from first envelope crossing of 10 percent of peak to last
    Double_t power_10_10;  /// Power enclosed within 10_10 width
    Double_t power_50_50;  /// Power enclosed within 50_50 width
    Double_t peakTime;  //[-128,384,18] /// Time that peak hilbert env occurs
    Double_t peakMoments[5];  /// moments about Peak  (1st - 5th moments)

    Double_t impulsivityMeasure; //[-1,1, 16]  /// A number that has something to do with how impulsive it is
    Double_t fracPowerWindowBegins[numFracPowerWindows]; //[0,128,16] /// Narrowest width containing {10%, 20%, 30%, 40%, 50%} of the total power
    Double_t fracPowerWindowEnds[numFracPowerWindows]; //[0,128,16] /// Narrowest width containing {10%, 20%, 30%, 40%, 50%} of the total power

    Int_t numAntennasInCoherent; /// Number of antennas used to make this

    Double_t localMaxToMin; //[0,4096,21] /// Largest value of local max to neighbouring local min (see Acclaim::RootTools::getLocalMaxToMin)
    Double_t localMaxToMinTime; //[0,100,16] /// Time between local maxima and minima +ve means max is before min, -ve means min is before max
    Double_t globalMaxToMin; //[0,4096,21] /// Difference between maximum and minimum voltage
    Double_t globalMaxToMinTime; //[0,128,16] /// Time between maximum and minimum volts, +ve means max is before min, -ve means min is before max


    //some utilities for polarization info
    Double_t linearPolFrac() const;
    Double_t linearPolAngle() const;
    Double_t circPolFrac() const;
    Double_t totalPolFrac() const;

    Double_t standardizedPeakMoment(Int_t i) const;
    inline Double_t skewness(){return standardizedPeakMoment(3);}
    inline Double_t kurtosis(){return standardizedPeakMoment(4);}

    Double_t fracPowerWindowMean() const;
    Double_t fracPowerWindowGradient() const;
    Double_t fracPowerWindowIntercept() const;
    Double_t fracPowerWindowChisquare() const;

    ClassDefNV(WaveformInfo, 16);

   private:
    friend class AnitaEventSummary;
    void cacheQuantitiesDerivedFromNarrowestWidths() const;
    mutable AnitaEventSummary* fContainer; //! Disgusting hack, do not persist!
    mutable UInt_t fLastEventNumberCache; //! Caching variables, do not persist!
    mutable Double_t nwMeanCache;           //! Caching variables, do not persist!
    mutable Double_t nwGradCache;           //! Caching variables, do not persist!
    mutable Double_t nwInterceptCache;      //! Caching variables, do not persist!
    mutable Double_t nwChisquareCache;      //! Caching variables, do not persist!
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

    Double_t rms; //[0,1024,20]
    Double_t avgPower;
    Double_t snr; //[0,100,16]/// Signal to Noise of waveform
    Double_t peakHilbert;//[0,4096,21]  /// peak of hilbert envelope

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

    ClassDefNV(ChannelInfo, 5);
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
      WAIS,  // is actually Hpol wais in both A3 and A4
      LDB,
      SIPLE,
      TD,
      WAIS_V // the Vpol wais in A4
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
    Double_t meanPower[1+AnitaRing::kNotARing];
    Double_t medianPower[1+AnitaRing::kNotARing];
    Double_t meanPowerFiltered[1+AnitaRing::kNotARing];
    Double_t medianPowerFiltered[1+AnitaRing::kNotARing];

    Double_t maxBottomToTopRatio[AnitaPol::kNotAPol];
    Int_t maxBottomToTopRatioSector[AnitaPol::kNotAPol];
    Double_t minBottomToTopRatio[AnitaPol::kNotAPol];
    Int_t minBottomToTopRatioSector[AnitaPol::kNotAPol];


    Int_t nSectorsWhereBottomExceedsTop;

    /** The fraction of nearby events that are payload blasts */
    Double_t blastFraction;

    ClassDefNV(EventFlags,11);
  };





  /**
   * @class SourceHypothesis
   * For known sources, such as a calibration pulser
   */
  class SourceHypothesis
  {
   public:
    SourceHypothesis() { reset(); }
    Double_t theta;
    Double_t phi;
    Double_t distance;

    Double_t mapValue[NUM_POLS];  ///what the instantaneous map value is at this source hypothesis
    Double_t mapHistoryVal[NUM_POLS]; /// a history of the interferometric map value for the source location

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
    Double_t weight;
    Double_t energy;
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

    /**
     * Convert to an Adu5Pat, (mostly to then instantiate a Adu5Pat)
     * @return an Adu5Pat only with location information
     */
    Adu5Pat pat () const {
      Adu5Pat pat;
      pat.longitude = longitude;
      pat.latitude = latitude;
      pat.altitude = altitude;
      pat.heading = heading;
      return pat;
    }

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
  Bool_t update() const {resetNonPersistent(); return true;}

  // Utilities to find interesting entries in AnitaEventSummary
  AnitaPol::AnitaPol_t highestPol() const;
  Int_t highestPolAsInt() const;
  Int_t highestPeakInd() const;
  const PointingHypothesis& highestPeak() const;
  const WaveformInfo& highestCoherent() const;
  const WaveformInfo& highestDeconvolved() const;
  const WaveformInfo& highestCoherentFiltered() const;
  const WaveformInfo& highestDeconvolvedFiltered() const;
	
  static void setThresholdForMostImpulsive(double threshold); /// value between 0 and 1, will find the brightest peak that is within threshold as impulsive as the most impulsive peak
  AnitaPol::AnitaPol_t mostImpulsivePol(int whichMetric=0) const;
  Int_t mostImpulsivePolAsInt(int whichMetric=0) const;
  Int_t mostImpulsiveInd(int whichMetric=0) const;
  const PointingHypothesis& mostImpulsivePeak(int whichMetric=0) const;
  const WaveformInfo& mostImpulsiveCoherent(int whichMetric=0) const;
  const WaveformInfo& mostImpulsiveDeconvolved(int whichMetric=0) const;
  const WaveformInfo& mostImpulsiveCoherentFiltered(int whichMetric=0) const;
  const WaveformInfo& mostImpulsiveDeconvolvedFiltered(int whichMetric=0) const;

  inline Double_t weight(){return mc.weight > 0 ? mc.weight : 1;} /// Return the weight of the event, always returns 1 for data, the weight from MCTruth otherwise
  AnitaPol::AnitaPol_t trainingPol() const;
  Int_t trainingPolAsInt() const;
  Int_t trainingPeakInd() const;
  const PointingHypothesis& trainingPeak() const;
  const WaveformInfo& trainingCoherent() const;
  const WaveformInfo& trainingDeconvolved() const;
  const WaveformInfo& trainingCoherentFiltered() const;
  const WaveformInfo& trainingDeconvolvedFiltered() const;

  //------------------------------------------------------------------------------------
 private:

  /*************************************************************************************
   * Private member variables/functions
   *************************************************************************************/
  // WARNING! THE COMMENTS TRAILING THESE AFFECT WHAT HAPPENS IN ROOT!
  // The //! comment initializer means ROOT does not store these variables when writing to files.
  // We want to keep it this way, since they are only used to cache the results of the utility
  // funtions to find the highest peak and peak nearest the monte carlo truth info.
  mutable Int_t                fHighestPeakIndex;  //! DOES NOT PERSIST IN ROOT! Internal index to cache result of finding highest peak
  mutable AnitaPol::AnitaPol_t fHighestPol;        //! DOES NOT PERSIST IN ROOT! Internal index to cache result of finding highest peak
  mutable Int_t                fTrainingPeakIndex; //! DOES NOT PERSIST IN ROOT! Internal index to cache result of finding peak training peak (pulser/mc peak)
  mutable AnitaPol::AnitaPol_t fTrainingPol;       //! DOES NOT PERSIST IN ROOT! Internal index to cache result of finding peak training peak (pulser/mc peak)
  mutable Int_t                fMostImpulsiveIndex;//! DOES NOT PERSIST IN ROOT! Internal index to cache result of finding most impulsive waveform
  mutable AnitaPol::AnitaPol_t fMostImpulsivePol;  //! DOES NOT PERSIST IN ROOT! Internal index to cache result of finding most impulsive waveform
  mutable UInt_t               fLastEventNumber;   //! DOES NOT PERSIST IN ROOT! To check for stale caching variables


  void findHighestPeak() const;
  void findTrainingPeak() const;
  void findMostImpulsive(int whichMetric) const;
  void resetNonPersistent() const;
  const SourceHypothesis* sourceFromTag() const;

  ClassDefNV(AnitaEventSummary, 36);
};

#endif
