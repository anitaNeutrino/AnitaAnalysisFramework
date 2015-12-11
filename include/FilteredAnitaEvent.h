#ifndef _FILTERED_ANITA_EVENT_H_
#define _FILTERED_ANITA_EVENT_H_

class UsefulAnitaEvent; 
class Adu5Pat; 
class TGraph; 
class RawAnitaHeader; 

#include "TObject.h" 
#include <vector>
#include <list>

class FilterStrategy; 
class FilteredAnitaEvent 
{
  public: 
   FilteredAnitaEvent(const UsefulAnitaEvent * event, const FilterStrategy * strategy, const Adu5Pat * pat, const RawAnitaHeader * header ); 
   const TGraph * getFilteredGraph(UInt_t i); 
   const UsefulAnitaEvent* getUsefulAnitaEvent() { return useful; } 
   const Adu5Pat * getGPS() { return pat; } 
   const RawAnitaHeader * getHeader() { return header; } 
  private: 
   std::vector<TGraph*> filteredGraphs; 
   const UsefulAnitaEvent * useful; 
   const FilterStrategy * strategy; 
   const Adu5Pat * pat; 
   const RawAnitaHeader * header; 
};


class FilterOperation
{
  public: 
    // short name for operation, will be used for output tree name, if there is one 
    virtual const char * tag () const  = 0; 

    // human readable description, should provide sufficient information to understand what was done 
    virtual const char * description () const  = 0; 

    virtual void process(FilteredAnitaEvent * event, double * tree_vars = 0) const = 0; 


    /* If you ant to output stuff to trees, define the output here */ 
    virtual UInt_t nTreeVars() const  { return 0; } 
    virtual const char *  treeVarName(UInt_t i) const  { (void) i; return ""; } 

}; 



class FilterStrategy
{
  friend class FilteredAnitaEvent; 

  public: 
    FilterStrategy() {; } 
    virtual ~FilterStrategy() {; } 
    void addOperation(const FilterOperation* f) { operations.push_back(f); }

  private: 
   std::list<const FilterOperation*> operations;
}; 



#endif

