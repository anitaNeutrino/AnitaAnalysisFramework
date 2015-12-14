#ifndef _FILTER_OPERATION_H
#define _FILTER_OPERATION_H

class FilteredAnitaEvent; 
class TGraph;
#include <cstddef>

class FilterOperation
{
  public: 
    /** short name for operation, will be used for output tree name, if there is one  */ 
    virtual const char * tag () const  = 0; 

    /** human readable description, should provide sufficient information to understand what was done  */ 
    virtual const char * description () const  = 0; 

    /**  operate on the FilteredAnitaEvent */ 
    virtual void process(FilteredAnitaEvent * event)= 0; 


    /* If you want to output stuff to trees, define the output here */ 
    virtual unsigned nOutputs() const  { return 0; } 
    virtual const char *  outputName(unsigned i) const  { (void) i; return ""; } 
    virtual void fillOutputs(double * vars) const{ (void) vars; return; } 

  protected: 
    TGraph * getGraph(FilteredAnitaEvent *ev, int i); 
    size_t nGraphs(FilteredAnitaEvent *ev); 

}; 


class ConditionalFilterOperation : public FilterOperation
{

  public: 
    ConditionalFilterOperation(bool (*condition)(FilteredAnitaEvent * ev, int trace))
      : fn(condition) { ; } 

    virtual void process(FilteredAnitaEvent * event); 

  protected:
    virtual TGraph * processIf(TGraph * g) = 0; /* return 0 if inplace */ 
    bool (*fn)(FilteredAnitaEvent *, int);
};


/** For filter operations that do the same thing to each Graph */ 
class UniformFilterOperation : public FilterOperation
{

  public: 
    virtual void process(FilteredAnitaEvent * event) ; 

  protected: 
    virtual TGraph * processOne(TGraph * g) = 0; /* return 0 if inplace */ 
}; 





#endif 
