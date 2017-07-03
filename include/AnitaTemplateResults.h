#ifndef ANITA_TEMPLATE_RESULTS
#define ANITA_TEMPLATE_RESULTS

#include "TObject.h"
#include "AnitaEventSummary.h"

class AnitaTemplateResults
{
 public:

  AnitaTemplateResults();

  /** The maximum number of hypotheses storable per polarization */ 
  static const Int_t maxDirectionsPerPol = AnitaEventSummary::maxDirectionsPerPol; 

  /* Number of points on the cone */
  static const int numCRTemplates = 10;

  /*The template correlatin values for a single coherent waveform comparison*/
  class SingleTemplateResult {
  public:
    SingleTemplateResult() {; }
    //impulse response
    Double_t templateImp;
    
    //one for the WAIS template too
    Double_t templateWais;
    
    //and for the bigger multi-coherence-angle one
    Double_t templateCRay[numCRTemplates];
    
    ClassDefNV(SingleTemplateResult,1);
  };

  
  SingleTemplateResult coherentV[maxDirectionsPerPol];
  SingleTemplateResult coherentH[maxDirectionsPerPol];
  
  SingleTemplateResult deconvolvedV[maxDirectionsPerPol];
  SingleTemplateResult deconvolvedH[maxDirectionsPerPol];


  void zeroInternals();
    
 private:
  ClassDefNV(AnitaTemplateResults, 1); 
};

#endif
