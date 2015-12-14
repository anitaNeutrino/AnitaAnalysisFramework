#include "FilterOperation.h" 
#include "TGraph.h" 
#include "FilteredAnitaEvent.h" 
#include <string.h> 


size_t FilterOperation::nGraphs(FilteredAnitaEvent *ev) 
{ 
  return ev->filteredGraphs.size(); 
}

TGraph* FilterOperation::getGraph(FilteredAnitaEvent *ev, int i) 
{ 
  return ev->filteredGraphs[i]; 
}


void ConditionalFilterOperation::process(FilteredAnitaEvent * ev) 
{
  for (size_t i = 0; i < nGraphs(ev); i++) 
  {
    if (!fn(ev,(int)i)) continue; 
    TGraph * g = getGraph(ev,i); 
    TGraph * gnew = processIf(g); 
    if (gnew) 
    {
      g->Set(gnew->GetN()); 
      memcpy(g->GetX(), gnew->GetX(), gnew->GetN() * sizeof(double)); 
      memcpy(g->GetY(), gnew->GetY(), gnew->GetN() * sizeof(double)); 
      delete gnew; 
    }
  }
}

void UniformFilterOperation::process(FilteredAnitaEvent * ev) 
{
  for (size_t i = 0; i < nGraphs(ev); i++) 
  {
    TGraph * g = getGraph(ev,i); 
    TGraph * gnew = processOne(g); 
    if (gnew) 
    {
      g->Set(gnew->GetN()); 
      memcpy(g->GetX(), gnew->GetX(), gnew->GetN() * sizeof(double)); 
      memcpy(g->GetY(), gnew->GetY(), gnew->GetN() * sizeof(double)); 
      delete gnew; 
    }
  }
}
