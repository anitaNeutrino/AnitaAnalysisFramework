#include "FilteredAnitaEvent.h" 
#include "AnitaConventions.h" 
#include "TGraph.h" 



FilteredAnitaEvent:: FilteredAnitaEvent(const UsefulAnitaEvent * event, const FilterStrategy * strategy, const Adu5Pat * pat, const RawAnitaHeader * header ) 
  : filteredGraphs(NUM_DIGITZED_CHANNELS), 
    useful(event), 
    strategy(strategy), 
    pat(pat), 
    header(header) 
{


  // Initialize the filtered graphs with the raw graphs from Raw Anita Event 
  for (UIint i = 0; i < NUM_DIGITZED_CHANNELS; i++) 
  {
    // UsefulAnitaEvent is not (yet?) const correct 
    filteredGraphs[i] = ((UsefulAnitaEvent*) useful)->getGraph(i); 
  }

  // Loop through the operations and apply them sequentially 
  for (std::list<const FilterOperation *>::const_iterator it = strategy->operations.begin(); it != strategy->operations.end(); it++) 
  {

    (*it)->process(this); 
  }
}








