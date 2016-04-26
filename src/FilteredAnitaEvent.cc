#include "FilteredAnitaEvent.h" 
#include "UsefulAnitaEvent.h"
#include "FilterStrategy.h"
#include "AnalysisWaveform.h"
#include "AnitaGeomTool.h"



//ClassImp(FilteredAnitaEvent); 

FilteredAnitaEvent::FilteredAnitaEvent() 
{
  for (unsigned i = 0; i < 2 * NUM_SEAVEYS; i++ )
    {
    filteredGraphs[i] = NULL; 
    rawGraphs[i] = NULL; 
  }
}



FilteredAnitaEvent:: FilteredAnitaEvent(const UsefulAnitaEvent * event, FilterStrategy * strategy, const Adu5Pat * pat, const RawAnitaHeader * header ) 
  : useful(event), 
    strategy(strategy), 
    pat(pat), 
    header(header) 

{

  // Initialize the filtered graphs with the raw graphs from Raw Anita Event 
  
  int k = 0; 
  for (int pol = AnitaPol::kHorizontal; pol <= AnitaPol::kVertical; pol++)
  {
    for (unsigned ant = 0; ant < NUM_SEAVEYS; ant++) 
    {
      int i = AnitaGeomTool::Instance()->getChanIndexFromAntPol(ant, (AnitaPol::AnitaPol_t) pol); 
      filteredGraphs[k] = new AnalysisWaveform (useful->fNumPoints[i], useful->fTimes[i], useful->fVolts[i]); 
      rawGraphs[k] = new AnalysisWaveform (useful->fNumPoints[i], useful->fTimes[i], useful->fVolts[i]); 
      filteredGraphs[k]->forceEvenSize(260); // do this for correlations 
      rawGraphs[k]->forceEvenSize(260); // do this for correlations 
      filteredGraphsByAntPol[pol][ant] = filteredGraphs[k];  
      rawGraphsByAntPol[pol][ant] = rawGraphs[k];  
      k++; 
    }
  }

  //tell the strategy to process this
  strategy->process(this); 
}



FilteredAnitaEvent::~FilteredAnitaEvent() 
{
  for (unsigned i = 0; i < 2 * NUM_SEAVEYS; i++ )
  {
    delete filteredGraphs[i]; 
    delete rawGraphs[i]; 
  }
}







