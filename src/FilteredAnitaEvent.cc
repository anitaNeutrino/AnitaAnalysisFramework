#include "FilteredAnitaEvent.h" 
#include "UsefulAnitaEvent.h"
#include "AnitaConventions.h" 
#include "FilterStrategy.h"
#include "AnalysisWaveform.h"



//ClassImp(FilteredAnitaEvent); 

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
    filteredGraphs[i] = new AnalysisWaveform (useful->fNumPoints[i], useful->fTimes[i], useful->fVolts[i]); 
  }

  //tell the strategy to process this
  strategy->process(this); 
}



FilteredAnitaEvent::~FilteredAnitaEvent() 
{
  for (unsigned i = 0; i < NUM_DIGITZED_CHANNELS; i++ )
  {
    delete filteredGraphs[i]; 
  }
}







