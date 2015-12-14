#ifndef _FILTERED_ANITA_EVENT_H_
#define _FILTERED_ANITA_EVENT_H_

class UsefulAnitaEvent; 
class Adu5Pat; 
class TGraph; 
class RawAnitaHeader; 

#include "TObject.h" 
class FilterStrategy; 

class FilteredAnitaEvent 
{

  friend class FilterOperation; 

  public: 
   FilteredAnitaEvent(const UsefulAnitaEvent * event, FilterStrategy * strategy, const Adu5Pat * pat, const RawAnitaHeader * header); 
   const TGraph * getFilteredGraph(UInt_t i) const; 
   const UsefulAnitaEvent* getUsefulAnitaEvent() { return useful; } 
   const Adu5Pat * getGPS() { return pat; } 
   const RawAnitaHeader * getHeader() { return header; } 
  private: 
   std::vector<TGraph*> filteredGraphs; 
   const UsefulAnitaEvent * useful; 
   const FilterStrategy * strategy; 
   const Adu5Pat * pat; 
   const RawAnitaHeader * header; 

   ClassDef(FilteredAnitaEvent,1); 
};





#endif

