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
}

ConditionalFilterOperation::ConditionalFilterOperation(FilterOperation * operation, 
                                                       bool (*condition)(FilteredAnitaEvent * ev, int trace), 
                                                       const char * condition_tag_suffix, const char * condition_description_suffix) 
                                                       : fn(condition), fo(operation)
{
  asprintf(&condition_tag, "%s_%s", fo->tag(), condition_tag_suffix); 
  asprintf(&condition_desc, "%s (if %s) ", fo->description(), condition_description_suffix); 
}



size_t FilterOperation::nGraphs(FilteredAnitaEvent *ev) 
{ 
  return ev->filteredGraphs.size(); 
}

AnalysisWaveform* FilterOperation::getWf(FilteredAnitaEvent *ev, int i) 
{ 
  return ev->filteredGraphs[i]; 
}



void ConditionalFilterOperation::process(FilteredAnitaEvent * ev) 
{
  for (size_t i = 0; i < nGraphs(ev); i++) 
  {
    if (fn(ev,(int)i))
    {
      fo->process(ev); 
    }
  }
}

void UniformFilterOperation::process(FilteredAnitaEvent * ev) 
{
  for (size_t i = 0; i < nGraphs(ev); i++) 
  {
    processOne(getWf(ev,i)); 
  }
}
