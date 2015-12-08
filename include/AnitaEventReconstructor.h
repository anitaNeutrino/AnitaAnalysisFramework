#ifndef _ANITA_EVENT_RECONSTRUCTOR_H
#define _ANITA_EVENT_RECONSTRUCTOR_H

class FilteredAnitaEvent; 
class AnitaEventSummary; 
class AnitaEventReconstructor 
{
  public:
   virtual void process(const FilteredAnitaEvent * ev, AnitaEventSummary * summary) const = 0; 
}; 

#endif 
