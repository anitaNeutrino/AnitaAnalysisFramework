#ifndef _FILTER_OPERATION_H
#define _FILTER_OPERATION_H

class FilteredAnitaEvent; 
class AnalysisWaveform; 
#include "AnitaConventions.h" 

/** A FilteredOperation does things to the waveforms inside a FilteredAnitaEvent 
 *
 *  See anitaGround/cozzyd/UCorrelator for some example implementing filters. 
 *
 **/ 

class FilterOperation
{
  public: 
    /** short name for operation, will be used for output tree name, if there is one  */ 
    virtual const char * tag () const = 0; 

    /** human readable description, should provide sufficient information to understand what was done  */ 
    virtual const char * description () const = 0; 
 
    /**  operate on the FilteredAnitaEvent */ 
    virtual void process(FilteredAnitaEvent * event)= 0; 


    // If you want to output stuff to trees, define the output here 
    
    /** The number of output variables (doubles or double arrays) defined by this operation */
    virtual unsigned nOutputs() const  { return 0; } 

    /** The name of the ith output variable */ 
    virtual const char *  outputName(unsigned i) const  { (void) i; return ""; } 

    /** The length of the ith output variable  (it's a double array of this size)*/ 
    virtual unsigned outputLength(unsigned i) const { (void) i; return 0; } 

    /** Fill the ith output */ 
    virtual void fillOutput(unsigned i, double * v) const{ (void) v; (void) i;  return; } 

    /** Destructor */ 
    virtual ~FilterOperation(); 

  protected: 
    /** Accessor for waveform */ 
    AnalysisWaveform * getWf(FilteredAnitaEvent *ev, int i); 

    /** Accessor for waveform */ 
    AnalysisWaveform * getWf(FilteredAnitaEvent *ev, int ant, AnitaPol::AnitaPol_t pol); 
}; 


/** For filter operations that do the same thing to each waveform */ 
class UniformFilterOperation : public FilterOperation
{

  public: 
    /** Processes an event, calling processOne on each waveform */ 
    virtual void process(FilteredAnitaEvent * event) ; 
    virtual ~UniformFilterOperation() {; } 

    /** You only need to implement this method */ 
    virtual void processOne(AnalysisWaveform * g) = 0; 

}; 



/** A ConditionalFilterOperation only applies the passed FilterOperation 
 *  if the condition is true. The tag and description are combinations of the 
 *  passed operation and the condition tag and description provided. 
 */ 
class ConditionalFilterOperation : public FilterOperation
{

  public: 
    /** Convert a UniformFilterOperation to a conditional operation. The UniformFilterOperation will only be applied
     * to the event if the passed condition function returns true */ 

    ConditionalFilterOperation(UniformFilterOperation * operation, 
                               bool (*condition)(FilteredAnitaEvent * ev, int ant, AnitaPol::AnitaPol_t pol), 
                               const char * condition_tag, const char * condition_description, bool own = false) ; 
    
    

    virtual const char * tag() const { return condition_tag; } 
    virtual const char * description () const { return condition_desc; } 

    virtual ~ConditionalFilterOperation(); 
    virtual void process(FilteredAnitaEvent * event); 

  protected:
    bool (*fn)(FilteredAnitaEvent *, int, AnitaPol::AnitaPol_t);
    char * condition_tag; 
    char * condition_desc; 
    UniformFilterOperation * fo; 
    bool own; 
};






#endif 
