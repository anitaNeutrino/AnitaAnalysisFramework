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

FilteredAnitaEvent:: FilteredAnitaEvent(const UsefulAnitaEvent * event, FilterStrategy * strategy, const Adu5Pat * rawpat, const RawAnitaHeader * header, bool save_stages ) 
  : useful(event), 
    strategy(strategy), 
    pat(rawpat), 
    header(header),
    keep_all_stages(save_stages)

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
//      rawGraphs[k]->forceEvenSize(260); // do this for correlations 
      filteredGraphsByAntPol[pol][ant] = filteredGraphs[k];  
      rawGraphsByAntPol[pol][ant] = rawGraphs[k];  
      k++; 
    }
  }

  //tell the strategy to process this
  strategy->process(this); 
}



const AnalysisWaveform * FilteredAnitaEvent::getFilteredGraphAtStage(UInt_t ant, AnitaPol::AnitaPol_t pol, UInt_t stage) const
{
  if (!keep_all_stages)
  {
    fprintf(stderr,"You didn't ask to save all the stages!\n"); 
    return 0; 
  }

  if (stage >= all_stages[pol][ant].size()) return filteredGraphsByAntPol[pol][ant]; 
  else return all_stages[pol][ant][stage]; 
}


void FilteredAnitaEvent::saveStage(int nreserve)
{

  for (int pol = 0; pol < 2; pol++)
  {
    for (int ant = 0; ant < NUM_SEAVEYS; ant++)
    {
      all_stages[pol][ant].reserve(nreserve); 
      all_stages[pol][ant].push_back(new AnalysisWaveform(*filteredGraphsByAntPol[pol][ant])); 
    }
  }
}

FilteredAnitaEvent::~FilteredAnitaEvent() 
{
  for (unsigned pol = 0; pol < 2; pol++)
  {
    for (unsigned ant = 0; ant <  NUM_SEAVEYS; ant++ )
    {
      delete filteredGraphsByAntPol[pol][ant]; 
      delete rawGraphsByAntPol[pol][ant];  

      for (unsigned j = 0; j < all_stages[pol][ant].size() ; j++)
      {
        delete all_stages[pol][ant][j]; 
      }
    }
  }
}







