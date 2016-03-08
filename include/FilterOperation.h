#ifndef _FILTER_OPERATION_H
#define _FILTER_OPERATION_H

class FilteredAnitaEvent; 
class AnalysisWaveform; 
#include <cstddef>

class FilterOperation
{
  public: 
    /** short name for operation, will be used for output tree name, if there is one  */ 
    virtual const char * tag () const = 0; 

    /** human readable description, should provide sufficient information to understand what was done  */ 
    virtual const char * description () const = 0; 
 
    /**  operate on the FilteredAnitaEvent */ 
    virtual void process(FilteredAnitaEvent * event)= 0; 


    /* If you want to output stuff to trees, define the output here */ 
    virtual unsigned nOutputs() const  { return 0; } 
    virtual const char *  outputName(unsigned i) const  { (void) i; return ""; } 
    virtual void fillOutputs(double * vars) const{ (void) vars; return; } 
    virtual ~FilterOperation(); 

  protected: 
    AnalysisWaveform * getWf(FilteredAnitaEvent *ev, int i); 
    size_t nGraphs(FilteredAnitaEvent *ev); 
}; 



/** A ConditionalFilterOperation only applies the passed FilterOperation 
 *  if the condition is true. The tag and description are combinations of the 
 *  passed operation and the condition tag and description provided. 
 */ 
class ConditionalFilterOperation : public FilterOperation
{

  public: 
    ConditionalFilterOperation(FilterOperation * operation, 
                               bool (*condition)(FilteredAnitaEvent * ev, int trace), 
                               const char * condition_tag, const char * condition_description) ; 
    
    

    virtual const char * tag() const { return condition_tag; } 
    virtual const char * description () const { return condition_desc; } 

    virtual ~ConditionalFilterOperation(); 
    virtual void process(FilteredAnitaEvent * event); 

  protected:
    bool (*fn)(FilteredAnitaEvent *, int);
    char * condition_tag; 
    char * condition_desc; 
    FilterOperation * fo; 
};


/** For filter operations that do the same thing to each Graph */ 
class UniformFilterOperation : public FilterOperation
{

  public: 
    virtual void process(FilteredAnitaEvent * event) ; 
    virtual ~UniformFilterOperation() {; } 
    virtual void processOne(AnalysisWaveform * g) = 0; 

  protected: 
}; 





#endif 
