#include "AnitaEventSummary.h"

ClassImp(AnitaEventSummary)



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * Takes care of copying the header and GPS info into the event summary
 */
AnitaEventSummary::AnitaEventSummary(const RawAnitaHeader* head,
				     const Adu5Pat* pat){

  this->header = (*head);
  this->gps = (*pat);  
  
}
