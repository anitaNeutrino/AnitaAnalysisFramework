#ifndef _FILTERED_ANITA_EVENT_H_
#define _FILTERED_ANITA_EVENT_H_
#include "TObject.h"
#include "AnitaConventions.h"
#include "UsefulAdu5Pat.h"
#include "AnalysisWaveform.h"

class UsefulAnitaEvent;
class Adu5Pat;
class TGraph;
class RawAnitaHeader;
class TCanvas;

class FilterStrategy;


/**
  \brief
  This class is intended to store all the necessary data about an ANITA event for filtering and analysis.
  It stores the raw and filtered waveforms as well as auxilliary information.
*/

class FilteredAnitaEvent
{

  friend class FilterOperation;
  friend class FilterStrategy;

  public:


   /** Create a FilteredAnitaEvent from a useful event, a filter strategy, GPS info and the event header. If keep_all_stages is true, then  all stages of filtering are kept. Note that filtering occurs during construction, not later!*/
   FilteredAnitaEvent(const UsefulAnitaEvent * event, FilterStrategy * strategy, const Adu5Pat * pat, const RawAnitaHeader * header, bool keep_all_stages = false);
   /**
      Created a FilteredAnitaEvent from another FilteredAnitaEvent.
      !THIS IS NOT A COPY CONSTRUCTOR!
      The output of fEv will be used as an input to create this filtered event.
      Useful for doing things like deconvolving the output of another filter.
      If keep_all_stages is true, then  all stages of filtering are kept. Note that filtering occurs during construction, not later!*/
   FilteredAnitaEvent(const FilteredAnitaEvent* fEv, FilterStrategy * strategy, bool keep_all_stages = false);
  

   /** Empty FilteredAnitaEvent. Not particularly useful */
   FilteredAnitaEvent();

   /** Destructor */
   virtual ~FilteredAnitaEvent();

   /** Accessor for raw waveform based on index. Mostly useful if you want to iterate over everything */
   const AnalysisWaveform * getRawGraph(UInt_t i) const { return rawGraphs[i]; }

   /** Accessor for raw waveform based on antenna number and polarization. */
   const AnalysisWaveform * getRawGraph(UInt_t ant, AnitaPol::AnitaPol_t pol) const { return rawGraphsByAntPol[pol][ant]; }

   /** Accessor for raw waveform based on phi, ring, and pol */
   const AnalysisWaveform * getRawGraph(UInt_t phi, AnitaRing::AnitaRing_t ring, AnitaPol::AnitaPol_t pol) const;

   /** Accessor for filtered waveform based on index. Mostly useful if you want to iterate over everything */
   const AnalysisWaveform * getFilteredGraph(UInt_t i) const { return filteredGraphs[i]; }

   /** Accessor for filtered waveform based on antenna number and polarization. */
   const AnalysisWaveform * getFilteredGraph(UInt_t ant, AnitaPol::AnitaPol_t pol) const { return filteredGraphsByAntPol[pol][ant]; }

   /** Accessor for raw waveform based on phi, ring, and pol */
   const AnalysisWaveform * getFilteredGraph(UInt_t phi, AnitaRing::AnitaRing_t ring, AnitaPol::AnitaPol_t pol) const;


   /* If the FilteredAnitaEvent was told to keep_all_stages, you can get the waveform at each stage. Stage 0 is after the first filter operation. */
   const AnalysisWaveform * getFilteredGraphAtStage(UInt_t ant, AnitaPol::AnitaPol_t pol, UInt_t stage) const;

   /** Return the UsefulAnitaEvent this event was constructed with. */
   const UsefulAnitaEvent* getUsefulAnitaEvent() const { return useful; }
   /** Return the Useful GPS info */
   const UsefulAdu5Pat * getGPS() const { return &pat; }

   /** Return the header */
   const RawAnitaHeader * getHeader() const { return header; }

   /** Return the strategy */
   const FilterStrategy* getStrategy() const {return strategy;}

   void plotSummary(TCanvas * chpol = 0, TCanvas * cvpol = 0) const;

   int checkSaturation(ULong64_t *save_hsat  =0, ULong64_t* save_vsat = 0, double threshold=1500) const; 
   
	 int checkStepFunction(Int_t lab = 1, AnitaRing::AnitaRing_t ring = AnitaRing::kMiddleRing, Int_t phiSector = 8, AnitaPol::AnitaPol_t pol = AnitaPol::kVertical) const; 

   /** Various calculations. Don't necessarily have to be in this class. */
   void getAverageSpectrum (TGraph * target, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol ) const;
   void getMedianSpectrum  (TGraph * target, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol , double pctile = 0.5) const;
   double getAveragePower(AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, AnitaRing::AnitaRing_t ring = AnitaRing::kNotARing, bool filtered = false) const;
   double getMedianPower(AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol,  AnitaRing::AnitaRing_t ring = AnitaRing::kNotARing, bool filtered = false) const;
   void getMinMaxRatio(AnitaPol::AnitaPol_t pol, double * max_ratio, double * min_ratio, int* max_sector, int* min_sector, AnitaRing::AnitaRing_t ring1 = AnitaRing::kBottomRing, AnitaRing::AnitaRing_t ring2 = AnitaRing::kTopRing, int nth = 0, int * n_greater = 0) const;

  int getAnitaVersion() const { return anitaVersion; }

  
  protected:
   AnalysisWaveform *rawGraphs[NUM_SEAVEYS*AnitaPol::kNotAPol];
   AnalysisWaveform *rawGraphsByAntPol[AnitaPol::kNotAPol][NUM_SEAVEYS];
   AnalysisWaveform *filteredGraphs[NUM_SEAVEYS*AnitaPol::kNotAPol];
   AnalysisWaveform *filteredGraphsByAntPol[AnitaPol::kNotAPol][NUM_SEAVEYS];


   const UsefulAnitaEvent * useful;
   const FilterStrategy * strategy;
   UsefulAdu5Pat pat;
   const RawAnitaHeader * header;

   int anitaVersion;

   bool keep_all_stages;
   std::vector<AnalysisWaveform *> all_stages[AnitaPol::kNotAPol][NUM_SEAVEYS] ;
   void saveStage(int nreserve);
};





#endif
