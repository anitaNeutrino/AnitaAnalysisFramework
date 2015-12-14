#ifndef _FILTER_STRATEGY_H
#define _FILTER_STRATEGY_H

class TFile; 
class TTree; 
class FilterOperation; 
class FilteredAnitaEvent; 

#include <vector> 
#include <set> 
#include <string> 

class FilterStrategy
{
  friend class FilteredAnitaEvent; 

  public: 
    /** Create a new empty strategy. If a pointer to a TFile is given, then trees may be written to that file if the operations define any output values */ 
    FilterStrategy(TFile * outfile = 0); 
    virtual ~FilterStrategy() {done(); } 
    void addOperation(FilterOperation* f); 
    void process(FilteredAnitaEvent * event); 
    void done(); 

  private: 
     bool started; 
     std::vector<FilterOperation*> operations;
     TFile * f; 
     std::vector<TTree *> trees; 
     std::vector<std::vector<double> > outputStore; 
     std::multiset<std::string> used_ids; 
}; 

#endif 

