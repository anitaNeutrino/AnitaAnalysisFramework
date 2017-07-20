#ifndef ANITA_TEMPLATES
#define ANITA_TEMPLATES

#include "TObject.h"
#include "AnitaEventSummary.h"
#include "FFTtools.h"
#include "TFile.h"
#include "AnalysisWaveform.h"
#include "SystemResponse.h"


/* * * * * * * * *
   
   
   Be careful using templates, we might not know what we're looking for


   Anyway this is where all the template stuff will live


 * * * * * * * * * */


/**

   _AnitaTemplateMachine_

    Class to import and hold onto the different templates
    
    Saveable to a ROOT file I think if you want to do that sort of thing!
*/

//they rely on each other so you gotta define it here
class AnitaTemplateSummary;

class AnitaTemplateMachine : public TObject {

 public:
  /** Default Constructor **/
  AnitaTemplateMachine(const int inLength=2048);

  /** Default Destructor **/
  ~AnitaTemplateMachine();
  
  /** in case you want to reset for some reason **/
  void zeroInternals();

  /** template waveform defaults **/
  const int length;
  const double dT = 0.1; //in ns
  
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

  //deconvolved versions for all those
  //impulse response
  FFTWComplex *theImpTemplateFFT_deconv;
  TGraph *theImpTemplate_deconv;
  //ANITA3 averaged wais pulses
  FFTWComplex *theWaisTemplateFFT_deconv;
  TGraph *theWaisTemplate_deconv;
  //ZHAires simulated CRs at various angles
  FFTWComplex *theCRTemplateFFTs_deconv[numCRTemplates];
  TGraph *theCRTemplates_deconv[numCRTemplates];
  

  /** Filling the stored templates **/
  void loadTemplates();
  void deconvolveTemplates(AnitaResponse::DeconvolutionMethod *deconv);
  

  /** Flags */
  //Check to see if the templates are loaded in bad english
  bool isTmpltsLoaded() { return kTmpltsLoaded; };
  //Check to see if you filled the deconvolved templates too
  bool isTmpltsDeconv() { return kTmpltsDeconv; };

  /** use the templates, and a supplied AnalysisWaveform, to fill up the summary */
  void doTemplateAnalysis(const AnalysisWaveform *waveform, int poli, int dir, AnitaTemplateSummary *summary);
  void doDeconvolvedTemplateAnalysis(const AnalysisWaveform *waveform, const AnitaResponse::DeconvolutionMethod *deconv,
				     int poli, int dir, AnitaTemplateSummary *summary);


  /** Write templates to file **/
  void writeTemplatesToFile(TFile *outFile);
  
 private:
  /* flags to see what you still might need to do */
  bool kTmpltsLoaded;
  bool kTmpltsDeconv;
  
  /* hidden functions for filling the templates one at a time */
  void getImpulseResponseTemplate();
  void getWaisTemplate();
  void getCRTemplates();
  
  
  ClassDefNV(AnitaTemplateMachine, 1);
  
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
  static const int numCRTemplates = AnitaTemplateMachine::numCRTemplates;

  /*The template correlatin values for a single coherent waveform comparison*/
  class SingleTemplateResult 
  {
  public:
    SingleTemplateResult() {; }
    //impulse response
    Double_t impulse;
    Double_t impulse_loc;
    
    //one for the WAIS template too
    Double_t wais;
    Double_t wais_loc;
    
    //and for the bigger multi-coherence-angle one
    Double_t cRay[numCRTemplates];
    Double_t cRay_loc[numCRTemplates];
    
    ClassDefNV(SingleTemplateResult,2);
  };

  
  SingleTemplateResult coherent[NUM_POLS][maxDirectionsPerPol];
  
  SingleTemplateResult deconvolved[NUM_POLS][maxDirectionsPerPol];


  void zeroInternals();
    
 private:
  ClassDefNV(AnitaTemplateSummary, 2); 
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