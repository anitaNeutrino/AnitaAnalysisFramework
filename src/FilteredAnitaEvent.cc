#include "FilteredAnitaEvent.h" 
#include "UsefulAnitaEvent.h"
#include "AnitaConventions.h" 
#include "FilterStrategy.h"


FilteredAnitaEvent:: FilteredAnitaEvent(const UsefulAnitaEvent * event, FilterStrategy * strategy, const Adu5Pat * pat, const RawAnitaHeader * header ) 
  : filteredGraphs(NUM_DIGITZED_CHANNELS), 
    useful(event), 
    strategy(strategy), 
    pat(pat), 
    header(header) 
{

  // Initialize the filtered graphs with the raw graphs from Raw Anita Event 
  for (unsigned i = 0; i < NUM_DIGITZED_CHANNELS; i++) 
  {
    // UsefulAnitaEvent is not (yet?) const correct 
    filteredGraphs[i] = ((UsefulAnitaEvent*) useful)->getGraph((int) i); 
  }

  //tell the strategy to process this
  strategy->process(this); 
}








