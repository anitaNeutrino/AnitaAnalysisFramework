#ifndef _FILTER_STRATEGY_H
#define _FILTER_STRATEGY_H



class TFile; 
class TTree; 
class FilterOperation; 
class FilteredAnitaEvent; 

#include <vector> 
#include <set> 
#include <string> 

/** 
 *
 * \brief A filter strategy defines the sets of filters that are used and provides some introspection abilities. 
 * At its most basic level, the strategy will serially apply a set of filters 
 *
 */

class FilterStrategy
{
  friend class FilteredAnitaEvent; 
  friend class NoiseMonitor;

  public: 
    /** Create a new empty strategy. If a pointer to a TFile is given, then trees may be written to that file if the operations define any output values */ 
    FilterStrategy(TFile * outfile = 0);

    /** Allow attachment of outfile after construct if haven't called this->process() yet */
    void attachFile(TFile* outfile);

    /** Destructor. Will write to output file if necessary (and deleted any owned operations) */ 
    ~FilterStrategy(); 

    /** Adds an operation to the strategy. This may only be done before any events are processed 
     * @param f  The FilterOperation to add to this strategy
     * @param enable_output Whether or not to enable the Filter Output tree for this operation
     * @param take_ownership  whether or not to delete this operation on destruction of the strategy
     *
     * */ 
    void addOperation(FilterOperation* f, bool enable_output = false, bool take_ownership = false); 

    /** Process an event using this strategy */
    void process(FilteredAnitaEvent * event); 

    /** Retrieve the ith filter operation */ 
    const FilterOperation * getOperation(size_t i) const { return operations[i]; } 

    /** Count the number of operations */ 
    size_t nOperations() const { return operations.size(); } 


    /** Output a string describing the strategy. User responsible for freeing.  */ 
    // char *  describe() const; 

    /** finish */ 
    void done(); 

  private: 
     bool started; 
     std::vector<FilterOperation*> operations;
     std::vector<bool> enable_outputs; 
     std::vector<bool> owns; 
     TFile * f; 
     std::vector<TTree *> trees; 
     std::vector<std::vector<std::vector<double> > >outputStore; 
     std::multiset<std::string> used_ids; 
}; 

#endif 

