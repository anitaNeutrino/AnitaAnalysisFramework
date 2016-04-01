#include "FilterOperation.h" 
#include "TGraph.h" 
#include "FilteredAnitaEvent.h" 
#include <string.h> 


FilterOperation::~FilterOperation()
{
}

ConditionalFilterOperation::~ConditionalFilterOperation()
{
 free (condition_tag); 
 free (condition_desc); 
 if (own) delete fo; 
}

ConditionalFilterOperation::ConditionalFilterOperation(UniformFilterOperation * operation, 
                                                       bool (*condition)(FilteredAnitaEvent * ev, int ant, AnitaPol::AnitaPol_t), 
                                                       const char * condition_tag_suffix, const char * condition_description_suffix, bool should_own_operation) 
                                                       : fn(condition), fo(operation), own(should_own_operation)
{
  asprintf(&condition_tag, "%s_%s", fo->tag(), condition_tag_suffix); 
  asprintf(&condition_desc, "%s (if %s) ", fo->description(), condition_description_suffix); 
}




AnalysisWaveform* FilterOperation::getWf(FilteredAnitaEvent *ev, int i) 
{ 
  return ev->filteredGraphs[i]; 
}

AnalysisWaveform* FilterOperation::getWf(FilteredAnitaEvent *ev, int ant, AnitaPol::AnitaPol_t pol) 
{ 
  return ev->filteredGraphsByAntPol[pol][ant]; 
}


void ConditionalFilterOperation::process(FilteredAnitaEvent * ev) 
{
  for (int pol = AnitaPol::kHorizontal; pol <= AnitaPol::kVertical; pol++)
  {
    for (int ant = 0; ant <NUM_SEAVEYS; ant++) 
    {
      if (fn(ev,ant, (AnitaPol::AnitaPol_t)pol))
      {
        fo->processOne(getWf(ev,ant, (AnitaPol::AnitaPol_t)pol)); 
      }
    }
  }
}

void UniformFilterOperation::process(FilteredAnitaEvent * ev) 
{
  for (size_t i = 0; i < NUM_SEAVEYS * 2; i++) 
  {
    processOne(getWf(ev,i)); 
  }
}
