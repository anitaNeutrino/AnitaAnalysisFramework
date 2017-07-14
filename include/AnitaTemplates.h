#ifndef ANITA_TEMPLATES
#define ANITA_TEMPLATES

#include "TObject.h"
#include "AnitaEventSummary.h"
#include "FFTtools.h"
#include "TFile.h"
#include "AnalysisWaveform.h"

/* * * * * * * * *
   
   
   Be careful using templates, we might not know what we're looking for


   Anyway this is where all the template stuff will live


 * * * * * * * * * */




/**

   _AnitaTemplateManager_

    Class to import and hold onto the different templates
    
    Saveable to a ROOT file I think if you want to do that sort of thing!
*/

class AnitaTemplateManager : public TObject {

 public:
  /** Default Constructor **/
  AnitaTemplateManager(const int inLength=2048);

  /** Default Destructor **/
  ~AnitaTemplateManager();
  
  /** in case you want to reset for some reason **/
  void zeroInternals();

  /** template length **/
  const int length;
  
  /** Template Storage **/
  //impulse response
  FFTWComplex *theImpTemplateFFT;
  TGraph *theImpTemplate;
  //ANITA3 averaged wais pulses
  FFTWComplex *theWaisTemplateFFT;
  TGraph *theWaisTemplate;
  //ZHAires simulated CRs at various angles
  static const int numCRTemplates = 10; //this is the place it is defined for everything
  FFTWComplex *theCRTemplateFFTs[numCRTemplates];
  TGraph *theCRTemplates[numCRTemplates];
  
  /** Filling the stored templates **/
  void loadTemplates();
  
  /** Check to see if the templates are loaded in bad english */
  bool isTmpltsLoaded() { return kTmpltsLoaded; };
  
  
 private:
  /* check if you loaded stuff already */
  bool kTmpltsLoaded;
  
  
  /* hidden functions for filling the templates one at a time */
  void getImpulseResponseTemplate();
  void getWaisTemplate();
  void getCRTemplates();
  
  
  ClassDefNV(AnitaTemplateManager, 1);
  
};


/**

   _AnitaTemplateSummary_

    Class to store results of any template correlations

*/

class AnitaTemplateSummary
{
 public:

  /** Default Constructor **/
  AnitaTemplateSummary();

  /** Default Destructor **/
  virtual ~AnitaTemplateSummary();

  /** The maximum number of hypotheses storable per polarization */ 
  static const Int_t maxDirectionsPerPol = AnitaEventSummary::maxDirectionsPerPol; 

  /* Number of points on the cone */
  static const int numCRTemplates = AnitaTemplateManager::numCRTemplates;

  /*The template correlatin values for a single coherent waveform comparison*/
  class SingleTemplateResult 
  {
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
  ClassDefNV(AnitaTemplateSummary, 1); 
};


/**

   _AnitaTemplateAnalyzer_
   
   Another class for actually doing the template analysis.  In hindsight it doesn't need to be it's own class
   
   Uses AnitaTemplateManager and coherently summed waveforms from your favorite analysis package
   to fill up an AnitaTemplateSummary
   
*/
class AnitaTemplateAnalyzer
{
 public:
  
  /** Constructor **/
  AnitaTemplateAnalyzer();

  /** Do the analysis! **/
  void doTemplateAnalysis(const AnalysisWaveform *waveform, AnitaTemplateManager *manager, int poli, AnitaTemplateSummary *summary);

  
};



//-----------------
//random utilities
// vvvvvvvvvvvv


/**

   Stuff I use that should maybe be in FFTtools?
   
*/

namespace FFTtools{
  TGraph *normalizeWaveform(TGraph* inGraph);
  double *getCorrelationFromFFT(int length,const FFTWComplex *theFFT1, const FFTWComplex *theFFT2);
}







/**

   windowing functions!

*/

namespace WindowingTools {

  TGraph *windowWave(TGraph*, int&, const int, const int, const int, const int);
  
  
  TGraph *windowCut(TGraph *inGraph,int length);
  
  TGraph *windowDispersed(TGraph *inGraph, int &peakHilbertLoc);
  TGraph *windowDispersed(TGraph *inGraph);
  
  TGraph *windowEField(TGraph *inGraph, int &peakHilbertLoc);
  TGraph *windowEField(TGraph *inGraph);
  

  TGraph *windowWave(TGraph *inGraph, int &peakHilbertLoc,
		     const int upRampLoc = 50,   const int downRampLoc = 600,
		     const int upRampLen = 100, const int downRampLen = 400);

}






#endif
