#ifndef _FILTERED_ANITA_EVENT_H_
#define _FILTERED_ANITA_EVENT_H_
#include "TObject.h" 
#include "AnitaConventions.h" 

class UsefulAnitaEvent; 
class Adu5Pat; 
class TGraph; 
class RawAnitaHeader; 
class AnalysisWaveform; 

class FilterStrategy; 


class FilteredAnitaEvent 
{

  friend class FilterOperation; 

  public: 
   FilteredAnitaEvent(); 
   FilteredAnitaEvent(const UsefulAnitaEvent * event, FilterStrategy * strategy, const Adu5Pat * pat, const RawAnitaHeader * header); 
   virtual ~FilteredAnitaEvent(); 
   const AnalysisWaveform * getFilteredGraph(UInt_t i) const { return filteredGraphs[i]; }
   const AnalysisWaveform * getFilteredGraph(UInt_t ant, AnitaPol::AnitaPol_t pol) const { return filteredGraphsByAntPol[ant][pol]; }
   const UsefulAnitaEvent* getUsefulAnitaEvent() { return useful; } 
   const Adu5Pat * getGPS() { return pat; } 
   const RawAnitaHeader * getHeader() { return header; } 
  private: 
   AnalysisWaveform *filteredGraphs[NUM_DIGITZED_CHANNELS]; 
   AnalysisWaveform *filteredGraphsByAntPol[2][NUM_DIGITZED_CHANNELS/2]; 

   const UsefulAnitaEvent * useful; 
   const FilterStrategy * strategy; 
   const Adu5Pat * pat; 
   const RawAnitaHeader * header; 

};





#endif

