#include "AnitaTemplateResults.h"

//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaTemplateResults::AnitaTemplateResults() {
  zeroInternals();
}



void AnitaTemplateResults::zeroInternals(){

  for(Int_t dir=0; dir < maxDirectionsPerPol; dir++)
    {
      memset(&coherentV[dir],0,sizeof(SingleTemplateResult)); 
      memset(&coherentH[dir],0,sizeof(SingleTemplateResult)); 
      
      memset(&deconvolvedV[dir],0,sizeof(SingleTemplateResult)); 
      memset(&deconvolvedH[dir],0,sizeof(SingleTemplateResult)); 
    }
  
  return;
}

