#ifndef _FILTERED_ANITA_EVENT_H_
#define _FILTERED_ANITA_EVENT_H_

class UsefulAnitaEvent; 
class TGraph; 

#include "TObject.h" 
#include <vector>
#include <list>

class FilterOperation
{
  friend class FilteredAnitaEvent; 
  public: 
    virtual const char * description () const  = 0; 
    virtual void process(FilteredAnitaEvent * event) const; 
}; 



class FilterStrategy
{
  public: 
    FilterStrategy() {; } 
    virtual ~FilterStrategy() {; } 
    void addOperation(const FilterOperation* f); 

  private: 
   std::list<const FilterOperation*> operations;
}; 


class FilteredAnitaEvent 
{
  public: 
    /*Smart version */ 
   FilteredAnitaEvent(const UsefulAnitaEvent * event, const FilterStrategy * strategy, const Adu5Pat * pat, const RawAnitaHeader * header ); 
   const TGraph * getGraph(UInt_t i); 
   
  private: 
   std::vector<TGraph*> filteredGraphs; 
   const UsefulAnitaEvent * useful; 
   const FilterStrategy * strategy; 
   const Adu5Pat * pat; 
   const RawAnitaHeader * header; 
};


#endif

