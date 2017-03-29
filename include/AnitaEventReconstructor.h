#ifndef _ANITA_EVENT_RECONSTRUCTOR_H
#define _ANITA_EVENT_RECONSTRUCTOR_H

class FilteredAnitaEvent; 
class AnitaEventSummary;
class UsefulAdu5Pat;
class AnitaEventReconstructor 
{
  public:
    virtual void process(const FilteredAnitaEvent * ev, UsefulAdu5Pat* usefulPat, AnitaEventSummary * summary) const = 0; 
    virtual ~AnitaEventReconstructor() { ; }
}; 

#endif 










