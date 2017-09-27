#include "AnitaTemplates.h"

//---------------------------------------------------------------------------------------------------------
/** AnitaTemplateSummary
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaTemplateSummary::AnitaTemplateSummary() {
  zeroInternals();
}

AnitaTemplateSummary::~AnitaTemplateSummary() {
  //don't need to do anything I don't think
}


void AnitaTemplateSummary::zeroInternals(){

  for (Int_t poli=0; poli<NUM_POLS; poli++) {
    for (Int_t dir=0; dir < maxDirectionsPerPol; dir++)
      {
	memset(&coherent[poli][dir],0,sizeof(SingleTemplateResult)); 
	
	memset(&deconvolved[poli][dir],0,sizeof(SingleTemplateResult)); 
      }
  }

  return;
}




//---------------------------------------------------------------------------------------------------------
/** AnitaTemplateMachine
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaTemplateMachine::AnitaTemplateMachine(int inLength)
    : length(inLength), dT(0.1)
{

  kTmpltsLoaded = false;
  kTmpltsDeconv = false;

  //initialize the TGraphs to NULL at least, the FFTs are allocated in the "get" stuff.
  theImpTemplate = NULL;
  theWaisTemplate = NULL;
  for (int i=0; i < numCRTemplates; i++) {
    theCRTemplates[i] = NULL;
  }

  zeroInternals();
}

AnitaTemplateMachine::~AnitaTemplateMachine() {
  zeroInternals();
}


void AnitaTemplateMachine::zeroInternals() {
  /* in case you want to delete everything */

  if (theImpTemplate != NULL) { delete theImpTemplate; theImpTemplate = NULL; }
  if (theWaisTemplate != NULL) { delete theWaisTemplate; theWaisTemplate = NULL; }
  for (int i=0; i < numCRTemplates; i++) {
    delete theCRTemplates[i];
    theCRTemplates[i] = NULL;
  }

  /* if it is loaded, free the FFT stuff too */
  if (kTmpltsLoaded) {
    free(theImpTemplateFFT);
    free(theWaisTemplateFFT);
    for (int i=0; i < numCRTemplates; i++) {
      free(theCRTemplates[i]);
    }
  }
    
  
  kTmpltsLoaded = false;
  return;
}

void AnitaTemplateMachine::getImpulseResponseTemplate() {
  
  //and get the "averaged" impulse response as the template"
  char* installDir = getenv("ANITA_UTIL_INSTALL_DIR");
  std::stringstream name;
  name.str("");
  name << installDir << "/share/AnitaAnalysisFramework/responses/SingleBRotter/all.imp";
  TGraph *grTemplateRaw = new TGraph(name.str().c_str());
  //waveforms are normally a little over 1024 so lets pad to 2048 (defined in length above)
  TGraph *grTemplatePadded = FFTtools::padWaveToLength(grTemplateRaw,length);
  delete grTemplateRaw;


  //Make it start at zero
  double xStart = grTemplatePadded->GetX()[0];
  for (int pt=0; pt<grTemplatePadded->GetN(); pt++) grTemplatePadded->GetX()[pt] -= xStart;
    

  //then cut it back down with a window function (or not)
  //  int peakHilb = -1;
  //  TGraph *grTemplateCut = windowDispersed(grTemplatePadded,peakHilb);
  //  delete grTemplatePadded;
  //and finally normalize it (last step!)


  TGraph *grTemplate = FFTtools::normalizeWaveform(grTemplatePadded);
  delete grTemplatePadded;

  //give it a name
  grTemplate->SetName("templateImp");

  theImpTemplate = grTemplate;

  //and get the FFT of it as well, since we don't want to do this every single event
  theImpTemplateFFT = FFTtools::doFFT(grTemplate->GetN(),grTemplate->GetY());
  
  return;
}



void AnitaTemplateMachine::getCRTemplates() {

  char* installDir = getenv("ANITA_UTIL_INSTALL_DIR");
  std::stringstream name;
  name.str("");
  name << installDir << "/share/AnitaAnalysisFramework/templates/crTmpltsA3.root";
  TFile *inFile = TFile::Open(name.str().c_str());
  
  for (int i=0; i<numCRTemplates; i++) {
    //want to get graphs 13 through 24 (like in makeTemplate.C)
    int wave = i+13; //peak seems to be at around the 13th one, then by 23 it is basically zero
    name.str("");
    name << "disp" << wave;
    TGraph *grTemplateRaw = (TGraph*)inFile->Get(name.str().c_str());
    
    //waveforms are super long so we can just cut it to the window dimentions
    int peakHilb = -1;
    TGraph *grTemplateCut = WindowingTools::windowCut(grTemplateRaw,length);
    delete grTemplateRaw;

    //Make them start at zero
    double xStart = grTemplateCut->GetX()[0];
    for (int pt=0; pt<grTemplateCut->GetN(); pt++) grTemplateCut->GetX()[pt] -= xStart;
    
    
    //and finally normalize it (last step!)
    TGraph *grTemplate = FFTtools::normalizeWaveform(grTemplateCut);
    delete grTemplateCut;

    //give it a name
    grTemplate->SetName(name.str().c_str());
    
    //    std::cout << "CR Template " << i << " Length: " << grTemplate->GetN() << std::endl;
    
    //and get the FFT of it as well, since we don't want to do this every single event
    FFTWComplex *theTemplateFFT=FFTtools::doFFT(grTemplate->GetN(),grTemplate->GetY());
    theCRTemplates[i] = grTemplate;
    theCRTemplateFFTs[i] = theTemplateFFT;
  }
  
  inFile->Close();

  return;
}


    
void AnitaTemplateMachine::getWaisTemplate() {
  
  //and get the "averaged" impulse response as the template"
  char* installDir = getenv("ANITA_UTIL_INSTALL_DIR");
  std::stringstream name;
  name.str("");
  name << installDir << "/share/AnitaAnalysisFramework/templates/waisTemplateA3.root";
  TFile *inFile = TFile::Open(name.str().c_str());
  TGraph *grTemplateRaw = (TGraph*)inFile->Get("wais01TH");
  //the wais waveform is like N=2832, but most of it is dumb, so cut off the beginning
  //actually just window it!
  int peakHilb = -1;
  TGraph *grTemplateCut = WindowingTools::windowCut(grTemplateRaw,length);
  delete grTemplateRaw;


  //Make it start at zero
  double xStart = grTemplateCut->GetX()[0];
  for (int pt=0; pt<grTemplateCut->GetN(); pt++) grTemplateCut->GetX()[pt] -= xStart;
    

  //and then normalize it
  TGraph *grTemplate = FFTtools::normalizeWaveform(grTemplateCut);
  delete grTemplateCut;

  //give it a name
  grTemplate->SetName("templateWais");

  //  std::cout << "Wais Template Length: " << grTemplate->GetN() << std::endl;

  //and get the FFT of it as well, since we don't want to do this every single event
  theWaisTemplateFFT = FFTtools::doFFT(grTemplate->GetN(),grTemplate->GetY());
  theWaisTemplate = grTemplate;

  return;
}

void AnitaTemplateMachine::loadTemplates() {

  if (kTmpltsLoaded) zeroInternals();


  std::cout << "Loading templates, length=" << length << std::endl;
  getImpulseResponseTemplate();
  getWaisTemplate();
  getCRTemplates();

  kTmpltsLoaded = true;

  return;
}
 


void AnitaTemplateMachine::deconvolveTemplates(AnitaResponse::DeconvolutionMethod *deconv) {

  std::cout << "Deconvolving templates" << std::endl;
  
  std::stringstream name;
  
  if (!isTmpltsLoaded()) {
    std::cout << "Error in AnitaTemplateAnalyzer::deconvolveTemplates(): AnitaTemplateMachine says it has no templates!" << std::endl;
    std::cout << "      Try doing AnitaTemplateMachine::loadTemplates() at some point maybe? " << std::endl;
    return;
  }
  
  double dF = 1./(dT*length);

  int lengthFFT = (length/2 + 1);

  /** Impulse Response */
  std::cout << "0" << std::endl;
  theImpTemplateFFT_deconv = (FFTWComplex*)malloc(lengthFFT*sizeof(FFTWComplex));
  for (int i=0; i<lengthFFT; i++) {
    theImpTemplateFFT_deconv[i] = theImpTemplateFFT[i];
  }				    
  std::cout << "1" << std::endl;
  deconv->deconvolve(lengthFFT,dF,theImpTemplateFFT_deconv,theImpTemplateFFT);
  double* theImpTemplate_deconv_y = FFTtools::doInvFFT(length,theImpTemplateFFT_deconv);
  std::cout << "2" << std::endl;
  theImpTemplate_deconv = new TGraph(length,theImpTemplate->GetX(),theImpTemplate_deconv_y);
  theImpTemplate_deconv->SetName("deconvImp");

  /** Wais */
  theWaisTemplateFFT_deconv = (FFTWComplex*)malloc(lengthFFT*sizeof(FFTWComplex));
  for (int i=0; i<(length/2 + 1); i++) {
    theWaisTemplateFFT_deconv[i] = theWaisTemplateFFT[i];
  }				    
  deconv->deconvolve(lengthFFT,dF,theWaisTemplateFFT_deconv,theImpTemplateFFT);
  double* theWaisTemplate_deconv_y = FFTtools::doInvFFT(length,theWaisTemplateFFT_deconv);
  theWaisTemplate_deconv = new TGraph(length,theWaisTemplate->GetX(),theWaisTemplate_deconv_y);
  theWaisTemplate_deconv->SetName("deconvWais");


  /**CR templates */
  for (int cr = 0; cr<numCRTemplates; cr++ ) {
    theCRTemplateFFTs_deconv[cr] = (FFTWComplex*)malloc(lengthFFT*sizeof(FFTWComplex));
    for (int i=0; i<(length/2 + 1); i++) {
      theCRTemplateFFTs_deconv[cr][i] = theCRTemplateFFTs[cr][i];
    }				    
    deconv->deconvolve(lengthFFT,dF,theCRTemplateFFTs_deconv[cr],theImpTemplateFFT);
    double* theCRTemplate_deconv_y = FFTtools::doInvFFT(length,theCRTemplateFFTs_deconv[cr]);
    theCRTemplates_deconv[cr] = new TGraph(length,theCRTemplates[cr]->GetX(),theCRTemplate_deconv_y);
    name.str("");
    name << "deconvCR" << cr;
    theCRTemplates_deconv[cr]->SetName(name.str().c_str());
  }
  

  kTmpltsDeconv = true;

  return;

}




void AnitaTemplateMachine::doTemplateAnalysis(const AnalysisWaveform *waveform, int poli, int dir, AnitaTemplateSummary *templateSummary) {
  /* copied out of templateSearch.cc */


  if (!isTmpltsLoaded()) {
    std::cout << "Error in AnitaTemplateAnalyzer::doTemplateAnalysis(): AnitaTemplateMachine says it has no templates!" << std::endl;
    std::cout << "      Try doing AnitaTemplateMachine::loadTemplates() at some point maybe? " << std::endl;
    return;
  }


  //I actually want to do SOME filtering though... so sine subtract a single one?
  //      sineSub->processOne(coherentAnalysis,data->header(),);
  const TGraphAligned *coherentAligned = waveform->even();
  TGraph *coherent = new TGraph(coherentAligned->GetN(),coherentAligned->GetX(),coherentAligned->GetY());
  //make sure it is the same sampling rate as the templates and default
  TGraph *coherentResamp = FFTtools::getInterpolatedGraph(coherent,dT);
  delete coherent;
  //make sure it is the same length as the template
  TGraph *coherentPad = FFTtools::padWaveToLength(coherentResamp,length);
  delete coherentResamp;

  //normalize coherently summed waveform
  TGraph *normCoherent = FFTtools::normalizeWaveform(coherentPad);
  delete coherentPad;
  FFTWComplex *coherentFFT=FFTtools::doFFT(length,normCoherent->GetY());
  delete normCoherent; 


  //set up variables for calculating
  double maxValue,minValue,value,value_loc;
  bool value_pol;
  double *dCorr;
    
  //Impulse Response
  dCorr = FFTtools::getCorrelationFromFFT(length,theImpTemplateFFT,coherentFFT);
  maxValue = TMath::MaxElement(length,dCorr);
  minValue = TMath::Abs(TMath::MinElement(length,dCorr));
  if (TMath::Max(maxValue,minValue) == maxValue) {
    value = maxValue;
    value_loc = TMath::LocMax(length,dCorr);
    value_pol = true;
  }
  else {
    value = minValue;
    value_loc = TMath::LocMin(length,dCorr);
    value_pol = false;
  }
  templateSummary->coherent[poli][dir].impulse  = value;
  templateSummary->coherent[poli][dir].impulse_loc = value_loc;
  templateSummary->coherent[poli][dir].impulse_pol = value_pol;

  delete[] dCorr;

  //Wais
  dCorr = FFTtools::getCorrelationFromFFT(length,theWaisTemplateFFT,coherentFFT);
  maxValue = TMath::MaxElement(length,dCorr);
  minValue = TMath::Abs(TMath::MinElement(length,dCorr));
  if (TMath::Max(maxValue,minValue) == maxValue) {
    value = maxValue;
    value_loc = TMath::LocMax(length,dCorr);
    value_pol = true;
  }
  else {
    value = minValue;
    value_loc = TMath::LocMin(length,dCorr);
    value_pol = false;
  }
  templateSummary->coherent[poli][dir].wais  = value;
  templateSummary->coherent[poli][dir].wais_loc = value_loc;
  templateSummary->coherent[poli][dir].wais_pol = value_pol;

  delete[] dCorr;


  //Cosmic Ray Templates
  for (int i=0; i<numCRTemplates; i++) {
    dCorr = FFTtools::getCorrelationFromFFT(length,theCRTemplateFFTs[i],coherentFFT);
    maxValue = TMath::MaxElement(length,dCorr);
    minValue = TMath::Abs(TMath::MinElement(length,dCorr));
    if (TMath::Max(maxValue,minValue) == maxValue) {
      value = maxValue;
      value_loc = TMath::LocMax(length,dCorr);
      value_pol = true;
    }
    else {
      value = minValue;
      value_loc = TMath::LocMin(length,dCorr);
      value_pol = false;
    }
    templateSummary->coherent[poli][dir].cRay[i] = value;
    templateSummary->coherent[poli][dir].cRay_loc[i] = value_loc;
    templateSummary->coherent[poli][dir].cRay_pol[i] = value_pol;

    delete[] dCorr;
  }

  
  delete[] coherentFFT;

  return;
  
}

void AnitaTemplateMachine::doDeconvolvedTemplateAnalysis(const AnalysisWaveform *waveform, 
							 const AnitaResponse::DeconvolutionMethod *deconv,
							 int poli, int dir, AnitaTemplateSummary *templateSummary) {
  /* copied out of templateSearch.cc */


  if (!isTmpltsLoaded()) {
    std::cout << "Error in AnitaTemplateAnalyzer::doTemplateAnalysis(): AnitaTemplateMachine says it has no templates!" << std::endl;
    std::cout << "      Try doing AnitaTemplateMachine::loadTemplates() at some point maybe? " << std::endl;
    return;
  }

  //I actually want to do SOME filtering though... so sine subtract a single one?
  //      sineSub->processOne(deconvolvedAnalysis,data->header(),);
  const TGraphAligned *deconvolvedAligned = waveform->even();
  if (!deconvolvedAligned || deconvolvedAligned->GetN() == 0) {
    return; //apparently sometimes the deconvolved waveform doesn't get filled?
  }
  TGraph *deconvolved = new TGraph(deconvolvedAligned->GetN(),deconvolvedAligned->GetX(),deconvolvedAligned->GetY());
  //make sure it is the same length as the template
  TGraph *deconvolvedPad = FFTtools::padWaveToLength(deconvolved,length);
  delete deconvolved;

  //normalize deconvolvedly summed waveform
  TGraph *normDeconvolved = FFTtools::normalizeWaveform(deconvolvedPad);
  delete deconvolvedPad;
  FFTWComplex *deconvolvedFFT=FFTtools::doFFT(length,normDeconvolved->GetY());
  delete normDeconvolved; 



  //set up variables for calculating
  double maxValue,minValue,value,value_loc;
  double *dCorr;
    
  //Impulse Response
  dCorr = FFTtools::getCorrelationFromFFT(length,theImpTemplateFFT_deconv,deconvolvedFFT);
  maxValue = TMath::MaxElement(length,dCorr);
  minValue = TMath::Abs(TMath::MinElement(length,dCorr));
  if (TMath::Max(maxValue,minValue) == maxValue) {
    value = maxValue;
    value_loc = TMath::LocMax(length,dCorr);
  }
  else {
    value = minValue;
    value_loc = TMath::LocMin(length,dCorr);
  }
  templateSummary->deconvolved[poli][dir].impulse  = value;
  templateSummary->deconvolved[poli][dir].impulse_loc = value_loc;

  delete[] dCorr;

  //Wais
  dCorr = FFTtools::getCorrelationFromFFT(length,theWaisTemplateFFT_deconv,deconvolvedFFT);
  maxValue = TMath::MaxElement(length,dCorr);
  minValue = TMath::Abs(TMath::MinElement(length,dCorr));
  if (TMath::Max(maxValue,minValue) == maxValue) {
    value = maxValue;
    value_loc = TMath::LocMax(length,dCorr);
  }
  else {
    value = minValue;
    value_loc = TMath::LocMin(length,dCorr);
  }
  templateSummary->deconvolved[poli][dir].wais  = value;
  templateSummary->deconvolved[poli][dir].wais_loc = value_loc;

  delete[] dCorr;

  //Cosmic Ray Templates
  for (int i=0; i<numCRTemplates; i++) {
    dCorr = FFTtools::getCorrelationFromFFT(length,theCRTemplateFFTs_deconv[i],deconvolvedFFT);
    maxValue = TMath::MaxElement(length,dCorr);
    minValue = TMath::Abs(TMath::MinElement(length,dCorr));
    if (TMath::Max(maxValue,minValue) == maxValue) {
      value = maxValue;
      value_loc = TMath::LocMax(length,dCorr);
    }
    else {
      value = minValue;
      value_loc = TMath::LocMin(length,dCorr);
    }
    templateSummary->deconvolved[poli][dir].cRay[i]  = value;
    templateSummary->deconvolved[poli][dir].cRay[i] = value_loc;

    delete[] dCorr;
  }
  
  delete[] deconvolvedFFT;
  
  return;
  
}



void AnitaTemplateMachine::writeTemplatesToFile(TFile *outFile) {

  if (!isTmpltsLoaded() || !isTmpltsDeconv()) {
    std::cout << "Error in AnitaTemplateAnalyzer::writeTemplatesToFile(): AnitaTemplateMachine says it is missing templates!" << std::endl;
    std::cout << "      Try doing AnitaTemplateMachine::loadTemplates() and deconvolveTempaltes() at some point maybe? " << std::endl;
    return;
  }
  
  outFile->cd();
  theImpTemplate->Write();
  theImpTemplate_deconv->Write();
  theWaisTemplate->Write();
  theWaisTemplate_deconv->Write();
  for (int i=0; i<numCRTemplates; i++) {
    theCRTemplates[i]->Write();
    theCRTemplates_deconv[i]->Write();
  }


  return;
}




//---------------------------------------------------------------------------------------------------------
/**

   stuff I want to add to FFTtools

*/
TGraph *FFTtools::normalizeWaveform(TGraph *inGraph) {
  
  TGraph *outGraph = (TGraph*)inGraph->Clone();
  
  //normalize it ( as seen in macros/testTemplate.C )
  double waveSum = 0;
  for (int pt=0; pt<outGraph->GetN(); pt++) waveSum += pow(outGraph->GetY()[pt],2);
  for (int pt=0; pt<outGraph->GetN(); pt++) outGraph->GetY()[pt] /= TMath::Sqrt(waveSum / (outGraph->GetN()/4));

  return outGraph;

}


double *FFTtools::getCorrelationFromFFT(int length,const FFTWComplex *theFFT1, const FFTWComplex *theFFT2) 
{


    int newLength=(length/2)+1;
//     std::cout << "newLength " << newLength << std::endl;
    FFTWComplex *tempStep = new FFTWComplex [newLength];
    int no2=length>>1;
    for(int i=0;i<newLength;i++) {
	double reFFT1=theFFT1[i].re;
	double imFFT1=theFFT1[i].im;
	double reFFT2=theFFT2[i].re;
	double imFFT2=theFFT2[i].im;

	//Real part of output 
	tempStep[i].re=(reFFT1*reFFT2+imFFT1*imFFT2)/double(no2/2);
	//Imaginary part of output 
	tempStep[i].im=(imFFT1*reFFT2-reFFT1*imFFT2)/double(no2/2);
    }
//    std::cout << "finished messing around" << std::endl;
    double *theOutput=FFTtools::doInvFFT(length,tempStep);
//    std::cout << "got inverse" << std::endl;
    delete [] tempStep;
    return theOutput;

}



//---------------------------------------------------------------------------------------------------------------------
/**

   Windowing stuff I wrote

 */


TGraph *WindowingTools::windowCut(TGraph *inGraph,int length) {
    /* cut whatever to this length with a nice little tail */
    int tempPeak = -1;
    const int beforeTime = 500;
    const int wings = 10;
    int afterTime = length - beforeTime - 2*wings;
    
    return windowWave(inGraph,tempPeak,beforeTime,afterTime,wings,wings);
  }
  
TGraph *WindowingTools::windowDispersed(TGraph *inGraph, int &peakHilbertLoc) {  
    return windowWave(inGraph,peakHilbertLoc,
		      500,750,20,20);
  }
  
  
TGraph *WindowingTools::windowDispersed(TGraph *inGraph) {
    //overload!
    int tempPeak = -1;
    return windowDispersed(inGraph,tempPeak);
  }


TGraph *WindowingTools::windowEField(TGraph *inGraph, int &peakHilbertLoc) {
  return windowWave(inGraph,peakHilbertLoc,
		    30,30,2,2);
}


TGraph *WindowingTools::windowEField(TGraph *inGraph) {
  //overload!
  int tempPeak = -1;
  return windowEField(inGraph,tempPeak);
}



TGraph *WindowingTools::windowWave(TGraph *inGraph, int &peakHilbertLoc,
				   const int upRampLoc, const int downRampLoc,
				   const int upRampLen, const int downRampLen) {
  //defaults are for the impulse response

  bool debug =false;

  /*

    The noise after the waveform part is useless.  I need to window it to increase the correlation value

    Find the peak of the hilbert envelope, then go 5ns before it (50pts) and then do a hamming maybe like 60 after it?

   */


  //the following are window config params, in POINTS (not nanoseconds)
  // downRampLoc - how far after peak hilbert env to start tail hamming
  // upRampLoc - how far before peak hilbert to start hamming (well close)
  
  // upRampLen: how "long" the hamming should be (half period)
  // downRampLen - how "long" the hamming should be at the end (well close)

  
  //If I don't tell it where to window, it should figure it out
  if (peakHilbertLoc == -1) {
    TGraph *hilbert = FFTtools::getHilbertEnvelope(inGraph);
    peakHilbertLoc = TMath::LocMax(hilbert->GetN(),hilbert->GetY());
    delete hilbert; //no memory leaks!
  }

  int startLoc = peakHilbertLoc - upRampLoc;
  int stopLoc  = peakHilbertLoc + downRampLoc;

  if (stopLoc+downRampLen > inGraph->GetN()) {
    if (debug) std::cout << "****";
    int overrun = (stopLoc+downRampLen) - inGraph->GetN() + 1;
    startLoc -= overrun;
    stopLoc -= overrun;
  }

  if (debug) std::cout << "inGraph->GetN()=" << inGraph->GetN() << " startLoc=" << startLoc << " stopLoc=" << stopLoc;



  TGraph *outGraph = new TGraph();
  for (int pt=0; pt<inGraph->GetN(); pt++) {
    if (pt <= (startLoc-upRampLen)) {
      //     outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],0); //not zero out, CUT
      if (debug) std::cout << pt << " - 1" << std::endl;
    }
    else if (pt > (startLoc-upRampLen) && pt <= startLoc ) {
      int ptMod = pt - (startLoc-upRampLen);
      double modValue = 0.5-(TMath::Cos(ptMod * ( TMath::Pi()/upRampLen ))/2.);
      double value =  modValue * inGraph->GetY()[pt];
      outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],value);
      if (debug) std::cout << pt << " - 2 - " << ptMod << " - " << modValue << std::endl;
    }
    else if (pt > startLoc && pt <= stopLoc) {
      outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],inGraph->GetY()[pt]);
      if (debug) std::cout << pt << " - 3" << std::endl;
    }
    else if (pt > stopLoc && pt <= (stopLoc+downRampLen)) {
      double ptMod = pt - stopLoc;
      double modValue = (1+TMath::Cos(ptMod*( TMath::Pi()/downRampLen ))/2.) - 0.5;
      double value = modValue * inGraph->GetY()[pt];
      outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],value);
      if (debug) std::cout << pt << " - 4 - " << ptMod << " - " << modValue << std::endl;
    }
    else if (pt > stopLoc+downRampLen) {
      //      outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],0); //not zero out, CUT
      if (debug) std::cout << pt << " - 5" << std::endl;
    }
  }

  if (debug) std::cout << " outGraph->GetN()=" << outGraph->GetN() << std::endl;
  return outGraph;

}






