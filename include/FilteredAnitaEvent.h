#ifndef _FILTERED_ANITA_EVENT_H_
#define _FILTERED_ANITA_EVENT_H_
#include "TObject.h" 
#include "AnitaConventions.h" 
#include "UsefulAdu5Pat.h"

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
   const AnalysisWaveform * getRawGraph(UInt_t i) const { return rawGraphs[i]; }
   const AnalysisWaveform * getRawGraph(UInt_t ant, AnitaPol::AnitaPol_t pol) const { return rawGraphsByAntPol[pol][ant]; }
   const AnalysisWaveform * getFilteredGraph(UInt_t i) const { return filteredGraphs[i]; }
   const AnalysisWaveform * getFilteredGraph(UInt_t ant, AnitaPol::AnitaPol_t pol) const { return filteredGraphsByAntPol[pol][ant]; }
   const UsefulAnitaEvent* getUsefulAnitaEvent() const { return useful; } 
   const UsefulAdu5Pat * getGPS() const { return &pat; } 
   const RawAnitaHeader * getHeader() const { return header; } 
  private: 
   AnalysisWaveform *rawGraphs[NUM_SEAVEYS*AnitaPol::kNotAPol]; 
   AnalysisWaveform *rawGraphsByAntPol[AnitaPol::kNotAPol][NUM_SEAVEYS]; 
   AnalysisWaveform *filteredGraphs[NUM_SEAVEYS*AnitaPol::kNotAPol]; 
   AnalysisWaveform *filteredGraphsByAntPol[AnitaPol::kNotAPol][NUM_SEAVEYS]; 

   const UsefulAnitaEvent * useful; 
   const FilterStrategy * strategy; 
   UsefulAdu5Pat pat; 
   const RawAnitaHeader * header; 

};





#endif

