
#include "GeomFilter.h"
#include "AnalysisWaveform.h" 
#include "FFTWComplex.h" 
#include "FFTtools.h"


using namespace std;

void GeometricFilter::process(FilteredAnitaEvent* event) {
  //printf("GeomFilter::process \n");
  double meanFreqVert = 0;
  double meanFreqHoriz = 0;
  int nadirFlag = 0;
  groupFlag = 0;
  currentEvent = event;
  getNotchandBandwidth(nFreq, dbCut, nAntsToUse, nadirFlag, meanFreqVert, meanFreqHoriz);
  for (int a=0; a<notchFreqs.size(); ++a) {
    //printf("antenna %i, %li notches \n", a, notchFreqs[a].size());
    for (int k=0; k<notchFreqs[a].size(); ++k) {
      //printf("  notch %i: freq=%f band=%f phase=%f \n", k, notchFreqs[a][k], notchBands[a][k], notchPhase[a][k]);
    }
  }
//  int printFlag = 0;
  std::vector<double> uniquefreqs;
  std::vector<double> uniqueBands;
  std::vector<double>  uniquePhase;
  double bandWidth = 26;
  int didIFilter = 0;
  int didIFilterAboveSatellite = 0;
  for(int j=0;j<notchFreqs.size();j++){//loop through every antenna
    uniquefreqs = notchFreqs[j];//set vectors needed for functions
    uniqueBands = notchBands[j];
    uniquePhase = notchPhase[j];
    for(int k=0;k<(int)uniquefreqs.size();k++){
      
      if (uniquefreqs[k]!=-1 && printFlag==1) cout<<"WZ going to cut: "<<uniquefreqs[k]<<" MHz in Vert for ant == "<<j<<endl;
      if (uniquefreqs[k]!=-1) didIFilter++;
      if (uniquefreqs[k]!=-1 && uniquefreqs[k]>286) didIFilterAboveSatellite++;
      bandWidth = uniqueBands[k];	
      
      //printf("bandwidth is %f \n", bandWidth);      
      if(uniquefreqs[k]>0) applyAdaptiveFilter_singleAnt(uniquefreqs[k], bandWidth,0,j);//apply filter in Vert
      
    }//k
    
    if(uniquefreqs.size()>0) GeomMethod(j,0,uniquefreqs, uniqueBands,uniquePhase);//apply phase filter 
  }//j

  for(int j=0;j<notchFreqs_Horiz.size();j++){//loop through horizontal
    uniquefreqs = notchFreqs_Horiz[j];
    uniqueBands = notchBands_Horiz[j];
    uniquePhase = notchPhase_Horiz[j];
    for(int k=0;k<(int)uniquefreqs.size();k++){
      
      if (uniquefreqs[k]!=-1 && printFlag==1) cout<<"WZ going to cut: "<<uniquefreqs[k]<<" MHz in Horiz for ant == "<<j<<endl;
      
      bandWidth = uniqueBands[k];
      
      if(uniquefreqs[k]>0) applyAdaptiveFilter_singleAnt(uniquefreqs[k],bandWidth,1,j);//apply filter in Horiz
      
    }//k
    
    

    if(uniquefreqs.size()>0) GeomMethod(j,1,uniquefreqs, uniqueBands,uniquePhase);//apply phase filter in Horiz
  }//j


}

void GeometricFilter::processOne(AnalysisWaveform* wf, const RawAnitaHeader * header, int ant, int pol) {
  printf("GeomFilter::processOne -- not implemented \n");
}

// TODO deprecate nadirFlag?
void GeometricFilter::getNotchandBandwidth(int nfreq,double dbCut,int nAntennasToUse, int nadirFlag, float meanFreqVert, float meanFreqHoriz)  {
  //int printFlag = 0;
//  printFlag=1;
  //printf("GeomFilter::getNotchAndBandwidth \n");   
  //printf("nfreq=%i, dbCut=%f, nAnt=%i, nadirFlag=%i \n", nfreq, dbCut, nAntennasToUse, nadirFlag);     
  
  std::vector<int> whichAntennasToUse (nAntennasToUse,0);
  std::vector< std::vector<double> >antennaFreq(NUM_SEAVEYS, std::vector<double>(50));
  std::vector< std::vector<double> >antennaFreqHoriz(NUM_SEAVEYS, std::vector<double>(50));
  std::vector< std::vector<double> >antennaBandwidth(NUM_SEAVEYS, std::vector<double>(50));
  std::vector< std::vector<double> >antennaBandwidthHoriz(NUM_SEAVEYS, std::vector<double>(50));
  std::vector< std::vector<double> >PeakMag(NUM_SEAVEYS, std::vector<double>(50));
  std::vector< std::vector<double> >PeakMag_backup(NUM_SEAVEYS, std::vector<double>(50));
  std::vector< std::vector<double> >PeakMagHoriz(NUM_SEAVEYS, std::vector<double>(50));
  
  
  std::vector<double> uniquefreqs;
  std::vector<double> uniquefreqs1;
  std::vector<double> uniquefreqsHoriz;
  std::vector<double> uniquefreqsHoriz1;
  
  std::vector<double> uniquebandwidth;
  std::vector<double> uniquebandwidthHoriz;
  
  std::vector<double> uniquePhase;
  std::vector<double> uniquePhaseHoriz;
  
  std::vector<double> uniquePhase_bandwidth;
  std::vector<double> uniquePhaseHoriz_bandwidth;
  
  
  double frequencies[nfreq];
  double frequenciesHoriz[nfreq];
  double bandwidth[nfreq];
  double bandwidthHoriz[nfreq];
  double magPeak[nfreq];
  double magPeakHoriz[nfreq];
  //double peakPhi=0.;
  int antused=0;
  int n=0;
  
  peakPhi = 0;
  if(groupFlag==0) {
    getGroupsofAntennas(nAntennasToUse,nadirFlag);//get unique groups of antennas of size nAntennasToUSe
    groupFlag=1;
  }
  //printf("   %i antenna groups \n", antenna_groups_size);
  for(int antenna_groups=0;antenna_groups<antenna_groups_size;antenna_groups++){//go through all groupings
    
    peakPhi = (double) unique_phis[antenna_groups];
    if (printFlag==1) cout<<"antenna_groups is "<<antenna_groups<<" and peak phi is "<<peakPhi<<"\n";
    for(int k=0;k<nAntennasToUse;k++){
      whichAntennasToUse[k]=antenna_group_holder[antenna_groups][k];
      //cout<<"using antenna "<<whichAntennasToUse[k]<<"\n";
      
    }
    
    //find frequencies to cut, put them into frequencies
    adaptiveFilterPartialPayload(0, dbCut, nfreq,frequencies,bandwidth,magPeak, nAntennasToUse,whichAntennasToUse, meanFreqVert);//Find CW freqs and Bandwidth for group of antennas
    adaptiveFilterPartialPayload(1, dbCut, nfreq,frequenciesHoriz,bandwidthHoriz,magPeakHoriz, nAntennasToUse,whichAntennasToUse, meanFreqHoriz);//IN Horiz
    
    //put freqs cut into each antenna used
    for(int j=0;j<nAntennasToUse;j++){
      antused = whichAntennasToUse[j];
      n=0;
      for(int j1=0;j1<50;j1++){//Bookkeeping. Keep track of what frequencies are being cut for each individual antenna
	      if(antennaFreq[antused][j1]<=0 && n<nfreq){ 
	        antennaFreq[antused][j1]=frequencies[n];
	        antennaBandwidth[antused][j1]=bandwidth[n];
	        PeakMag[antused][j1]=magPeak[n];
	        antennaFreqHoriz[antused][j1]=frequenciesHoriz[n];
	        antennaBandwidthHoriz[antused][j1]=bandwidthHoriz[n];
	        PeakMagHoriz[antused][j1]=magPeakHoriz[n];
	        n++;
	      }//antennaFreq==0	 
      }//j1
    }//j
      
    
    
  }//antenna_groups
  
  notchFreqs.clear(); 
  notchBands.clear(); 
  notchPhase.clear(); 
  notchPhaseBands.clear(); 
  notchBands_Horiz.clear(); 
  notchFreqs_Horiz.clear(); 
  notchPhase_Horiz.clear(); 
  notchPhaseBands_Horiz.clear(); 
  for(int i =0;i<NUM_SEAVEYS;i++){//WE have all possible freqs to cut. Let go through each antenna and make sense of them
    uniquefreqs.clear();//zero stuff
     
    uniquefreqsHoriz.clear();
      
    uniquebandwidth.clear();
    uniquebandwidthHoriz.clear();
    uniquePhase.clear();
    uniquePhaseHoriz.clear();
    uniquePhase_bandwidth.clear();
    uniquePhaseHoriz_bandwidth.clear();
 
     
      //Get rid of duplicated freqs. Also combine overlapping notches to form one big notch

    if (i==44) {
      for (int k=0; k<antennaFreqHoriz[i].size(); ++k) {
        //printf("  antenna %i  frequency %f  \n", i, antennaFreqHoriz[i][k]);
      }
    }
     
    getFrequenciestoCut(i, antennaFreq,antennaBandwidth,PeakMag,uniquefreqs,uniquebandwidth,nfreq,uniquePhase,uniquePhase_bandwidth);//fill uniquefreq with freq to cut
    getFrequenciestoCut(i, antennaFreqHoriz,antennaBandwidthHoriz,PeakMagHoriz,uniquefreqsHoriz,uniquebandwidthHoriz,nfreq,uniquePhaseHoriz,uniquePhaseHoriz_bandwidth);

    //printf("frequencies to cut: antenna %i  hPol \n", i);
    for (int k=0; k<uniquefreqs.size(); ++k) {
      //printf("frequency %f \n", uniquefreqs[k]);
    }
    notchFreqs.push_back(uniquefreqs);//Set up vectors that will be used for actually applying cuts
    notchBands.push_back(uniquebandwidth);
    notchPhase.push_back(uniquePhase);
    notchPhaseBands.push_back(uniquePhase_bandwidth);
    notchFreqs_Horiz.push_back(uniquefreqsHoriz);
    notchBands_Horiz.push_back(uniquebandwidthHoriz);
    notchPhase_Horiz.push_back(uniquePhaseHoriz);
    notchPhaseBands_Horiz.push_back(uniquePhaseHoriz_bandwidth);
      
  }//i
  
  //printf("leaving GeomFilter::getNotchAndBandwidth \n");   
  
}

void GeometricFilter::getGroupsofAntennas(int nAntennasToUse, int nadirFlag){//gets unique groups of Antennas
  //int antenna_group_tmp [nAntennasToUse];//for use in getClosestN function
  std::vector<int> antenna_group (nAntennasToUse,0);//used for sorting and comparison
  int old_group_switch=0;
 
  for(int runoverphi=0;runoverphi<360;runoverphi++){//get all possible antenna groups
   
    old_group_switch=0;
    getClosestNAntennas(nAntennasToUse, runoverphi, antenna_group, nadirFlag);  
    //getClosestNAntennas(nAntennasToUse, runoverphi, antenna_group_tmp, nadirFlag);  
    //for(int i=0;i<nAntennasToUse;i++){
    //  antenna_group[i]=antenna_group_tmp[i];
    //}
    sort(antenna_group.begin(),antenna_group.end());
    
    for(int i0=0;i0<(int)antenna_group_holder.size();i0++){
      if(antenna_group==antenna_group_holder[i0]){
	      old_group_switch=1;
	      break;
      }//check same  
    }//i0 
    
    
    if(old_group_switch==0){
      antenna_group_holder.push_back(antenna_group);
      unique_phis.push_back(runoverphi);
     
    }
   
  }
  
  antenna_groups_size=(int) antenna_group_holder.size();

}

double getMaximum(int n, double *array, int &index){
  double max;
  max=array[0];
  index=0;
  for (int i=1;i<n;i++){
    if (array[i]>max){ 
      max=array[i];
      index=i;
    }
  }
  return max;
}


void GeometricFilter::getClosestNAntennas(int nantennasToUse, double peakPhi, vector<int>& whichAntennasToUse, int nadirFlag)
{

  int nantennasTotal=NUM_SEAVEYS;
  //if (nadirFlag==1) nantennasTotal=NUM_ANTS_WITH_NADIRS;
  //if (nadirFlag==0) nantennasTotal=NUM_ANTS_NO_NADIRS;
  
  double deltaPhiArray[nantennasToUse];
  double deltaPhiThisOne;  
  int index;
  double phi_ant[nantennasTotal];

  
  for (int i=0;i<nantennasToUse;i++){
    deltaPhiArray[i]=360;//set array full of large values
  }
  ULong64_t saturated[2] = {0,0};
  //const UsefulAnitaEvent* usefulEvent = currentEvent->getUsefulAnitaEvent();
  currentEvent->checkSaturation( &saturated[AnitaPol::kHorizontal], 
                                 &saturated[AnitaPol::kVertical], 
                                 1500);   
  // TODO populate saturated channels array by iterating through bit variables
//  printf("saturation hpol: %li   vpol: %li\n", saturated[0], saturated[1]);
  // horrible temporary hack code
  //int saturatedChannels[2*NUM_SEAVEYS] = {0};  
  int polToggle = 0;  // TODO maye iterate on this?
  // end of horrible temporary hack code
  for (int ant=0;ant<nantennasTotal;ant++){
    //cout<<"ant is "<<ant<<"\n";
    //if (ant!=1 && saturatedChannels[ant+polToggle*NUM_SEAVEYS]==0){
    if (saturatedChannels[ant+polToggle*NUM_SEAVEYS]==0){       //TODO really check saturation
      //phi_ant[ant]=fUPGeomTool->getAntPhiPositionRelToAftFore(ant);
      phi_ant[ant]=AnitaGeomTool::Instance()->getAntPhiPositionRelToAftFore(ant);
      deltaPhiThisOne=phi_ant[ant]*180.0/M_PI-peakPhi;
     
      if (deltaPhiThisOne>180.) deltaPhiThisOne-=360.;
      if (deltaPhiThisOne<-180.) deltaPhiThisOne+=360.;
      if (deltaPhiThisOne<0) deltaPhiThisOne*=-1.;
      
      if (deltaPhiThisOne<getMaximum(nantennasToUse, deltaPhiArray, index)){//find maximum of array/replace if new value is lower
	      deltaPhiArray[index]=deltaPhiThisOne;
	      whichAntennasToUse[index]=ant;
      }
    }//ant!=1
  }
}

///////////////////////////////////////
void GeometricFilter::adaptiveFilterPartialPayload(int pol, double dBCut, int nfreq,double *frequencies,double *bandwidth,double *magPeak,int nantennasToUse,vector<int>& whichAntennasToUse, float &mean_freq)//identifies CW frequencies and size of notch for group of antennas
{
  //printf("in GeometricFilter::adaptiveFilterPartialPayload \n");
  // use class variable AnalysisWaveform array "noise" instead
  //if (readBaselineFlag!=1){ //get Baselines for comparison purposes. Only done one time
  //  readBaselineFFTs_Brian();
  //  if (printFlag==1) cout<<"Read in Baseline For Adaptive Filter"<<endl;
  //}  
  //int printFlag = 0;
  const int nPoints = 129;//124
  double magFFT[2000];
  double frequencyArray[2000];
  double deltaF=-1.0;
  double deltaT;
  int length;
  int newLength;
  double *X, *Y;
  
  //int phi_sector =(int) round(peakPhi/7.5);
  
  //int whichAntennasToUse[nantennasToUse];
  int ant;
  for (int i=0;i<2000;i++){
    magFFT[i]=0;
   
    frequencyArray[i]=-1;
  }
 

  double baseline[nPoints]={0};
  
  //GetBaselineperPhi(pol,baseline,nantennasToUse, whichAntennasToUse);//take n antennas and get the avg baseline for those
  
  if(printFlag==1) cout<<"  peakPhi="<<peakPhi<<"\n";
 
  // for (int ant=0;ant<NUM_ANTS_WITH_NADIRS;ant++){  
  for (int ctr=0;ctr<nantennasToUse;ctr++){
    ant=whichAntennasToUse[ctr];
    AnalysisWaveform* rawAWF = getWf(currentEvent, ant, (AnitaPol::AnitaPol_t)pol);
    AnalysisWaveform* thisAWF = new AnalysisWaveform(*rawAWF);
    thisAWF->forceEvenSize(240);
    const TGraphAligned* thisPower = thisAWF->power();
    if (ctr==0) {
      deltaF = thisAWF->deltaF()*1000;
      //deltaT = thisAWF->deltaT();
      newLength = thisPower->GetN();
    }
    //if ((pol!=0 || ant!=1) && saturatedChannels[ant+pol*NUM_SEAVEYS]==0){//get rid of 2V and ant 1 in vert
    // TODO check saturation for real here - just apply the max criterion and forget the saturation[] array
    if (saturatedChannels[ant+pol*NUM_SEAVEYS]==0){// 
      // should not need to do the FFT; 
      //   magFFT wants dB
      /*
      if (pol==0) Y = grEv[ant]->GetY();//set X and Y values for FFT
      else{ Y = grEvHoriz[ant]->GetY(); }
      if (pol==0) X = grEv[ant]->GetX();
      else{ X = grEvHoriz[ant]->GetX(); }      
      deltaT=X[1]-X[0];
      if (pol==0) length=grEv[ant]->GetN();
      else{length=grEvHoriz[ant]->GetN();}      
      newLength=(length/2)+1;
      deltaF=1/(deltaT*length); //GHz
      deltaF*=1e3; //MHz
      FFTWComplex *theFFT=FFTtools::doFFT(length,Y);//do FFT
      */
      // populate magFFT with the frequency spectrum from the FilteredAnitaEvent::AnalysisWaveform (getWf?)
      
      //printf("antenna %i  power graph size is %i, dF=%f  \n", ant, thisPower->GetN(), thisAWF->deltaF());
      //for(int i=0;i<newLength;i++) {
      for(int i=0; i<newLength; i++) {
	      if (ctr==0){
	        if (i==0) frequencyArray[i]=0;
	        if (i>0) frequencyArray[i]=frequencyArray[i-1] + deltaF;//make a freq Array.
	      }
	      //printf("frequency %i: %f \n", i, frequencyArray[i]);
	      if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
	        //magFFT[i]+=theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im;//adding squares of all antennas together
	        //printf("adding in power %i: %f \n", i, thisPower->GetY()[i]);
	        magFFT[i]+= thisPower->GetY()[i];
	      }
	      else {
	        magFFT[i]=-1000;
	       
	      }
      }   
      //delete [] theFFT;
      
    }//pol!=0
    delete thisAWF;
  }//ctr==0

  for(int i=0;i<newLength;i++) {
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      if (pol==0) magFFT[i]=10*log10(sqrt(magFFT[i]/double(nantennasToUse))/10.);//correct RMS
      if (pol==1) magFFT[i]=10*log10(sqrt(magFFT[i]/double(nantennasToUse))/10.);//correct RMS      
    }    
    else {
      magFFT[i]=-1000;
     
    }
  }
 
  //for (int i=0; i<newLength; ++i) {printf("  %i: freq=%f  \n", i, frequencyArray[i]);}
 
  //get baseline average so we can bump FFT around
  double meanBaseline=0;
  double mean=0;
  double firstHalfMeanBaseline=0;
  double secondHalfMeanBaseline=0;
  double firstHalfMean=0;
  double secondHalfMean=0;
  int navg=0;
  int nfirsthalfavg=0;
  int nsecondhalfavg=0;

  // double *bX, *bY;
  double *bX_set,*bY_set;
  
  int n;
  int n_set;
  double old_mean=-1000.;

  vector<int> index_skip (newLength,1);
  // TODO add up the baslines for all antennasToUse
  
  int faSize = noiseSamples[0][0]->GetN();
  double xVals[faSize];
  double yValsH[faSize];
  double yValsV[faSize];
  for (int k=0; k<faSize; ++k) {
    xVals[k] = noiseSamples[0][0]->GetX()[k];
    yValsH[k] = 0;
    yValsV[k] = 0;
  }
    for (int k = 0; k < whichAntennasToUse.size();k++){
    	int a = whichAntennasToUse[k];
    	// check size of noise sample arrays and warn
    	if (noiseSamples[0][a]->GetN() != faSize) {printf(" warning - baseline input array size mismatch (hpol)\n");}
    	if (noiseSamples[1][a]->GetN() != faSize) {printf(" warning - baseline input array size mismatch (vpol)\n");}
    	for (int k=0; k<faSize; ++k) {
      	// hacky finagle sammy: antenna 5 hpol baseline is wonky, just use #6  
      	yValsV[k] += noiseSamples[1][a]->GetY()[k];
      	if (a==4) {
        	yValsH[k] += noiseSamples[0][5]->GetY()[k];
      	} else {
        	yValsH[k] += noiseSamples[0][a]->GetY()[k];
      }
    }
  }
  const TGraphAligned* grHorizCoherentBaseline = new TGraphAligned(faSize, xVals, yValsH);
  const TGraphAligned* grVertCoherentBaseline = new TGraphAligned(faSize, xVals, yValsV);
  // for testing just use antenna 1 baseline for all antennas
  //const TGraphAligned* grHorizCoherentBaseline = noiseSamples[0][0]; 
  //const TGraphAligned* grVertCoherentBaseline = noiseSamples[1][0];
  
  //printf("noise: \n"); for (int k=0; k<grHorizCoherentBaseline->GetN(); ++k) {printf("  %f GHz     %f \n", grHorizCoherentBaseline->GetX()[k], grHorizCoherentBaseline->GetY()[k]);}
  
  if (pol==0){ //vertical
    bY_set=grVertCoherentBaseline->GetY();//should be set to partial payload
    bX_set=grVertCoherentBaseline->GetX();//should be set to partial payload
    n_set=grVertCoherentBaseline->GetN();//should be set to partial payload
  }
  if (pol==1){ //horizontal
    bY_set=grHorizCoherentBaseline->GetY();
    bX_set=grHorizCoherentBaseline->GetX();
    n_set=grHorizCoherentBaseline->GetN();
  }
  double bX[n_set];
  double bY[n_set];

  vector< double> baseline_normal(n_set,0);
  vector< double> baseline_tilt(n_set,0);
  vector <double> baseline_shift(n_set,0);
  vector <double> baseline_shift1(n_set,0);
  int n_max=10;
  
  for(int n0=0;n0<n_max;n0++){//do iterative process to find correct placement of baseline
  
    if( fabs (mean-old_mean) <= 0.01){//if change in baseline is <.01, stop process
      n0 = n_max-1;
    }
    meanBaseline=0.;
    mean=0.;
    firstHalfMeanBaseline=0.;
    secondHalfMeanBaseline=0.;
    firstHalfMean=0.;
    secondHalfMean=0.;
    navg=0;
    nfirsthalfavg=0;
    nsecondhalfavg=0;
     
    for(int i0=0;i0<n_set;i0++){
      baseline_normal[i0]=0.;
      baseline_tilt[i0]=0.;
      baseline_shift[i0]=0.;
      bY[i0] = bY_set[i0];
      bX[i0] = bX_set[i0]*1000;
      n = n_set;
    }

    //for (int i=0; i<n; ++i) {
    //  printf("i=%i, bY=%f \n", i, bY[i]);
    //}
    for (int i=0;i<n;i++){
      baseline_normal[i]=bY[i];
      if (bX[i]>=200 && bX[i]<1200){
        meanBaseline+=bY[i];
        navg++;
      }
      if (bX[i]>=200 && bX[i]<700){ 
        firstHalfMeanBaseline+=bY[i];
        nfirsthalfavg++;
      }
      if (bX[i]>=700 && bX[i]<1200){ 
        secondHalfMeanBaseline+=bY[i];
        nsecondhalfavg++;
      }	
    }
    //printf("    navg=%i \n", navg);
    meanBaseline=meanBaseline/double(navg);
    firstHalfMeanBaseline=firstHalfMeanBaseline/double(nfirsthalfavg);
    secondHalfMeanBaseline=secondHalfMeanBaseline/double(nsecondhalfavg);
    
    navg=0;
    nfirsthalfavg=0;
    nsecondhalfavg=0;
    
   

    //get average of graph in question

    for (int i=0;i<newLength;i++){
      if (frequencyArray[i]>=200 && frequencyArray[i]<1200 && index_skip[i]==1){ 
         mean+=magFFT[i];
        
        navg++;
      }
      if (frequencyArray[i]>=200 && frequencyArray[i]<700 && index_skip[i]==1){
        firstHalfMean+=magFFT[i];
       
        nfirsthalfavg++;
      }
      if (frequencyArray[i]>=700 && frequencyArray[i]<1200 && index_skip[i]==1){
        secondHalfMean+=magFFT[i];
       
        nsecondhalfavg++;
      }

    } 
    mean=mean/double(navg);
    firstHalfMean=firstHalfMean/double(nfirsthalfavg);
    secondHalfMean=secondHalfMean/double(nsecondhalfavg);
    // cout<<"mean is "<<mean<<" old mean is "<<old_mean<<"\n";
    old_mean = mean;
    //printf("     mean=%f, meanBaseline = %f \n", mean, meanBaseline);
    //now bump the average to the baseline average and apply a tilt correction to baseline
    double deltaMean=mean-meanBaseline;
    double deltaMeanFirst=firstHalfMean-firstHalfMeanBaseline-deltaMean;
    double deltaMeanSecond=secondHalfMean-secondHalfMeanBaseline-deltaMean;
    double slope=(deltaMeanFirst-deltaMeanSecond)/500.;
    //printf("   firstHalfMeal=%f, secondHalfMean=%f, firstHalfMeanBaseline=%f, secondHalfMeanBaseline=%f \n", 
    //        firstHalfMean, secondHalfMean, firstHalfMeanBaseline, secondHalfMeanBaseline);
    //printf("   deltaMean=%f, deltaMeanFirst=%f, deltaMeanSecond=%f, slope=%f \n", deltaMean, deltaMeanFirst, deltaMeanSecond, slope);
   
    for(int ctr0=0;ctr0<n;ctr0++){
      baseline_shift[ctr0]=baseline_normal[ctr0]+deltaMean;
      //baseline_shift1[ctr0]=0.;
    }
    for (int i=0; i<n; ++i) {
      //printf("i=%i, bY=%f \n", i, bY[i]);
    }

    //printf("   slope=%f \n", slope);
    for (int ctr=0;ctr<n;ctr++){
      if (bX[ctr]>=200 && bX[ctr]<1200){
        bY[ctr]=bY[ctr]-slope*(bX[ctr]-700.);
        baseline_tilt[ctr]=bY[ctr];        
      }
      else{
        baseline_tilt[ctr]=-1000;
      }
    }

    for (int i=0; i<n; ++i) {
      //printf("i=%i, bY=%f \n", i, bY[i]);
    }
   
   
    for(int i=0;i<n;i++){
      bY[i]= bY[i]+deltaMean;
      //printf("    i=%i, bY=%f dM=%f \n", i, bY[i], deltaMean);
    }
   
    
    //find mean frequency
    mean_freq=0;
    double cum_power=0;
    double avg_power=0;
    for (int i=0;i<newLength;i++){
      if (frequencyArray[i]>200 && frequencyArray[i]<1200){
        avg_power+=pow(10,magFFT[i]/10);
      }
    }

    for (int i=0;i<newLength;i++){
      if (frequencyArray[i]>200 && frequencyArray[i]<1200){
        if (cum_power<avg_power/2.) cum_power+=pow(10,magFFT[i]/10);
        else{
        	mean_freq=frequencyArray[i];
	        break;
        }
      }
    }
    //if (printFlag==1) cout<<"Mean frequency before filtering: "<<mean_freq<<endl;
    if (printFlag==1) printf("Mean frequency before filtering: %f \n ", mean_freq);
    
    //now see if any peaks are ndB above the baseline.
    double deltaMag[newLength];

    int j;
    for (int i=0;i<newLength;i++){
      if (frequencyArray[i]>210 && frequencyArray[i]<1190){
        for (j=0;j<n;j++){
	        if (bX[j]>frequencyArray[i]) break;
        }
        deltaMag[i]=magFFT[i]-bY[j];
        //printf("      i=%i, freq=%f, mag=%f, baseline=%f, delta=%f \n", i, bX[j], magFFT[i], bY[j], deltaMag[i]);
      }
      else deltaMag[i]=-1000;
    }


    int j_index;
    int k_index;
    for(int i=0;i<nfreq;i++){
      int index;
      double maxDelta=getMaximum(newLength,deltaMag,index);
//      printf("  index=%i  freq=%f, maxDelta=%f \n", index, frequencyArray[index], maxDelta);
      if(maxDelta>dBCut){
//        printf("  choosing index=%i  freq=%f,  maxDelta=%f \n", index, frequencyArray[index], maxDelta);
        for(int j=index;j>0;j--){
	        if(deltaMag[j]<dBCut-1){
	          j_index=j;
	         
	          break;
	        }
        }
       
        for(int k=index;k<newLength;k++){
	        if(deltaMag[k]<dBCut-1){
	          k_index=k;
	         
	          break;
	        }
        }
        
        index = (int)ceil((k_index+j_index)/2);
        frequencies[i]=frequencyArray[index];
        bandwidth[i]=frequencyArray[k_index] - frequencyArray[j_index];
        
        magPeak[i]=maxDelta;
        if(bandwidth[i]<deltaF){
	        bandwidth[i]=deltaF;
        }

        for(int reset=j_index;reset<=k_index;reset++){
	  index_skip[reset]=0;
        }

        // cout<<"index is "<<index<<" freq to cut is "<<frequencyArray[index]<<" bandwidth is "<<bandwidth[i]<<"\n";
        for(int clear=j_index;clear<=k_index;clear++){
	  deltaMag[clear]=-1000;
        }
      }
      else{
        magPeak[i]=-1;
        frequencies[i]=-1;
      }
     }//nfreq
   }//n=0
  //if(printFlag==1) cout<<"about to delete grVert/Horiz Coherent \n";
  // no longer have to delete since we're pointing back to class variables
  delete grVertCoherentBaseline;
  delete grHorizCoherentBaseline;
  grVertCoherentBaseline=0;
  grHorizCoherentBaseline=0;

   if(printFlag==1) cout<<" leaving GeometricFilter::adaptiveFilterPartialPayload \n";
 
}

//////////-----------------------------------------------------------------------------------------------------------------------------------------

//////////-----------------------------------------------------------------------------------------------------------------------------------------

void GeometricFilter::getFrequenciestoCut(int antenna,vector< vector<double> > &antennaFreq,vector< vector<double> > &bandwidth, vector< vector<double> > &PeakMag, vector<double> &uniquefreqs,vector<double> &uniquebandwidth, int nfreq, vector<double> &uniquePhase, vector<double> &uniquePhase_bandwidth) {//Get rid of Duplicates and combines overlapping notches
  //printf("in GeometricFilter::getFrequenciesToCut");
  //printf(" antenna %i \n", antenna);
  double max=0.;
  int maxbin=0;
  double bandWidth=0.;
  vector< vector<double> > freq_mag(50, vector<double>(3));
  int new_freq_flag=0;
  int non_zero_flag=0;
  uniquefreqs.clear();
  uniquebandwidth.clear();
  //get freqs to cut

  for(int n=0;n<50;n++){
    //printf("antenna %i  index %i   PeakMag=%f \n", antenna, n, PeakMag[antenna][n]);
    max=0.;
    maxbin=0;
    for(int j=0;j<50;j++){
      if(PeakMag[antenna][j]>max){
	      max=PeakMag[antenna][j];
	      maxbin =j;
        //printf("    new peak maximum found n=%i, j=%i: %f \n", n, maxbin, max);
      }
    }
    if(max>0){
     
      freq_mag[n][0]=antennaFreq[antenna][maxbin];
      freq_mag[n][1]=max;
      freq_mag[n][2]=bandwidth[antenna][maxbin];
    }
    PeakMag[antenna][maxbin]=-1000;
  }//n
 
  for(int k0=0;k0<50;k0++){
    new_freq_flag=0;
    //printf("   n=%i, freq=%f, maxVal=%f, bandwidth=%f \n", k0, freq_mag[k0][0], freq_mag[k0][1], freq_mag[k0][2]);
    if(non_zero_flag==0){
      if(freq_mag[k0][0]> 1 && freq_mag[k0][0]<1210){
        //printf("  pushing unique (1) frequency %f \n", freq_mag[k0][0]);
	      uniquefreqs.push_back(freq_mag[k0][0]); 
	      non_zero_flag=1;
      }
    }
    for(int k1=0;k1<(int)uniquefreqs.size();k1++){
      if((int)freq_mag[k0][0]==(int)uniquefreqs[k1]){
	      new_freq_flag=0;
	      break;
      }
      else{
	      new_freq_flag=1;
      }
    }//k1
    
    if(new_freq_flag==1 && (int)uniquefreqs.size()<nfreq  && freq_mag[k0][0]> 190 && freq_mag[k0][0]<1210){
      //printf("  pushing unique (2) frequency %f \n", freq_mag[k0][0]);
      uniquefreqs.push_back(freq_mag[k0][0]); 
    }
  }//k0
  vector<double> uniqueMag;
  double mag;
  //vector<double> uniquePhase;
  //vector<double> uniquePhase_bandwidth;
  
  if((int)uniquefreqs.size() >0){
    sort(uniquefreqs.begin(),uniquefreqs.end());//sort all the freqs for this antenna

    for(int k0=0;k0<(int)uniquefreqs.size();k0++){
      bandWidth=0.;
      for(int k1=0;k1<50;k1++){
	      if(uniquefreqs[k0]==freq_mag[k1][0]){
	        if(freq_mag[k1][2]>bandWidth){
	          bandWidth=freq_mag[k1][2];
	          mag = freq_mag[k1][1];
	        }
	      }
      }
      if(bandWidth >1) {	
	      uniquebandwidth.push_back(bandWidth);
	      uniqueMag.push_back(mag);
      }
    }
    
    uniquePhase.push_back(uniquefreqs[0]);
    uniquePhase_bandwidth.push_back(uniquebandwidth[0]);
    int kappa=0;
    for(int trial0=0;trial0<(int)uniquefreqs.size();trial0++){
      
      if(uniquefreqs[trial0] < uniquePhase[kappa]-15 || uniquefreqs[trial0]>uniquePhase[kappa]+15){
	if(uniquefreqs[trial0] >200 && uniquefreqs[trial0]<1200){
	  uniquePhase.push_back(uniquefreqs[trial0]);
	  uniquePhase_bandwidth.push_back(uniquebandwidth[trial0]);
	  kappa++;
	}
      }
    }

    
    double middle_freq_value;
    
    double new_bandwidth;
    
    double lower1,lower2;
    double higher1,higher2;
  
    for(int tester=0;tester<(int)uniquefreqs.size()-1;tester++){
      lower1 = uniquefreqs[tester]-uniquebandwidth[tester];
      higher1 = uniquefreqs[tester]+uniquebandwidth[tester];
      
      lower2 = uniquefreqs[tester+1]-uniquebandwidth[tester+1];
      higher2 = uniquefreqs[tester+1]+uniquebandwidth[tester+1];
     
      if(lower1<=lower2-1 && higher1>=higher2+1){//2nd freq to cut lies inside first's bandwidth
	      //	cout<<"first scenario! \n\n";
	      uniquefreqs.erase (uniquefreqs.begin()+tester+1);	    
	      uniquebandwidth.erase (uniquebandwidth.begin()+tester+1);	    
	      tester=tester-1;
      }
      else if(lower2<=lower1-1 && higher2>=higher1+1){//1st freq to cut lies inside 2nds bandwidth
	      //	cout<<"second scenario! \n\n";
	      uniquefreqs.erase (uniquefreqs.begin()+tester);	    
	      uniquebandwidth.erase (uniquebandwidth.begin()+tester);	    
	      tester=tester-1;
      }      
      else if(higher1 >= lower2-1){//overlap. combine.
	      //cout<<"third scenario! \n\n";
	      middle_freq_value = (lower1 + higher2)/2;
	      new_bandwidth = middle_freq_value - lower1;
	      //cout<<"middle_freq_value is "<<middle_freq_value<<" new_bandwidth is "<<new_bandwidth<<"\n";
	      uniquefreqs[tester]=middle_freq_value;
	      uniquebandwidth[tester]=new_bandwidth;
	      uniquefreqs.erase (uniquefreqs.begin()+tester+1);
	      uniquebandwidth.erase (uniquebandwidth.begin()+tester+1);
	      tester=tester-1;
      }
    }
  }//uniquefreqs.size()>0

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeometricFilter::applyAdaptiveFilter_singleAnt(double centerFrequency, double bandWidth, int polFlag,int ant) //function to call coorrect filter and change graphs
{
  
  // TODO instead of grEv as event graph, use unfiltered, evenly-sampled waveform from FilteredAnitaEvent (getWf(currentEvent, ant, pol))
  AnalysisWaveform* thisAWF = getWf(currentEvent, ant, (AnitaPol::AnitaPol_t)polFlag);
  if (polFlag==0){
    if (centerFrequency!=-1){
      //TGraph *gr1 = new TGraph(grEv[ant]->GetN(),grEv[ant]->GetX(),grEv[ant]->GetY());
      //delete grEv[ant];
      const TGraphAligned* thisGrA = thisAWF->even();
      TGraph *gr1 = new TGraph(thisGrA->GetN(),thisGrA->GetX(),thisGrA->GetY());
      // feed it to the filter monster
      TGraph *grNotch=interpolatedFilter(gr1,centerFrequency-bandWidth,centerFrequency+bandWidth);
      //grEv[ant]=new TGraph(grNotch->GetN(),grNotch->GetX(),grNotch->GetY());
      thisAWF->updateEven(grNotch);
      delete gr1;
      delete grNotch;
           
    }//centerfreq
  }//polFlag
  else{
    if (centerFrequency!=-1){
      //TGraph *gr1Horiz=new TGraph(grEvHoriz[ant]->GetN(),grEvHoriz[ant]->GetX(),grEvHoriz[ant]->GetY());
      //delete grEvHoriz[ant];
      const TGraphAligned* thisGrA = thisAWF->even();
      TGraph *gr1Horiz = new TGraph(thisGrA->GetN(),thisGrA->GetX(),thisGrA->GetY());
      // feed it to the filter monster
      TGraph *grNotchHoriz=interpolatedFilter(gr1Horiz,centerFrequency-bandWidth,centerFrequency+bandWidth);
      //grEvHoriz[ant]=new TGraph(grNotchHoriz->GetN(),grNotchHoriz->GetX(),grNotchHoriz->GetY());
      thisAWF->updateEven(grNotchHoriz);
      delete gr1Horiz;
      delete grNotchHoriz;
            
    }//centerfreq
  }//else
  
}//AdaptiveFilter_singleAnt
////////-----------------------------------------------------------------------------------------------------------

////////-----------------------------------------------------------------------------------------------------------
TGraph * GeometricFilter::interpolatedFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq) //interpolated filter. Do interpolation across notch
{
  if(maxFreq>1200){
    maxFreq=1200;
  }
  //TRandom3 random; // for generating random numbers
  //random.SetSeed(fHeadPtr->eventNumber);
  //cout<<"minFreq is "<<minFreq<< " max freq is "<<maxFreq<<"\n";
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();//256
  double frequencyArray[2000];
  double magFFT[2000];

  for(int j0=0;j0<2000;j0++){
    frequencyArray[j0]=0;
    magFFT[j0]=0;
  }
 
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  
  int newLength=(length/2)+1;
  //    double fMax = 1/(2*deltaT);  // In Hz
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  double tempF=0.;
  double startFreq=0.;
  double endFreq=0.;
  int i_start=0;
  int i_end=0;
  for(int i=0;i<newLength;i++){
    //cout<<"minFreq is "<<minFreq<<" maxFreq is "<<maxFreq<<" tempF is "<<tempF<<"\n";
    if(tempF >= (minFreq-deltaF) && tempF <=(minFreq) ){
      //  cout<<"setting starts \n";
      startFreq = tempF;
      i_start = i;
    }
    if(tempF >= (maxFreq) && tempF <= (maxFreq+deltaF)){
      //  cout<<"Setting ends \n";
      endFreq = tempF;
      i_end=i;
    }
    tempF+=deltaF;  
  }
  tempF=0;
  double slope;
  double intercept;
  
 double mag_start = sqrt(pow(theFFT[i_start].re,2) + pow(theFFT[i_start].im,2));
  double mag_end = sqrt(pow(theFFT[i_end].re,2) + pow(theFFT[i_end].im,2));
  slope = (mag_end - mag_start)/(endFreq - startFreq);
  intercept = mag_start - (slope*startFreq);
  
  double mag=0.;
  double phi_old=0.;
 
  double x=0.;
  double y=0.;
 

  for(int i=0;i<newLength;i++) {
   
    if(tempF>startFreq && tempF<endFreq ) {
      phi_old = atan2(theFFT[i].im,theFFT[i].re);//old phase
      
      mag = slope*tempF + intercept;
     
      x = mag*cos(phi_old);
      y = mag*sin(phi_old);
      
      theFFT[i].re = x;
      theFFT[i].im = y;
    }//in inside freq range 
   
    tempF+=deltaF;
  }
  
  double *filteredVals = FFTtools::doInvFFT(length,theFFT);
  TGraph *grFilteredNotch=new TGraph(length,oldX,filteredVals);  //oldX

  // }
 
  delete [] theFFT;
 
  delete [] filteredVals;
 
  return grFilteredNotch; 

}//interpolated notch

/////-------------------------------------------------------------------------------------------------------------------

/////-------------------------------------------------------------------------------------------------------------------

void GeometricFilter::GeomMethod(int ant,int pol,vector<double> Freq,vector<double> bandWidth,vector<double> cutFreqs){ //geom method for phase filtering
  double minFreq;// = Freq - bandWidth;
  double maxFreq;// = Freq + bandWidth;
  double magFFT[2000]={0};
  double magFFT_dB[2000]={0};
  double frequencyArray[2000]={0};
  double deltaT, deltaF;
  int length;
  int newLength;
  double phase_single[2000];
  double phase_single_unshifted[2000];
  double phase_new[2000];
  
  double *times;
  double *volts;
 
   
  double delta; 
  //printf("getting waveform event %i antenna %i pol %i \n", currentEvent->eventNumber,  ant, pol);
  AnalysisWaveform* thisAWF = getWf(currentEvent, ant, (AnitaPol::AnitaPol_t) pol);
  //if(pol==0){
  //  length = grEv[ant]->GetN();
  //}
  //else if(pol==1){
  //  length = grEvHoriz[ant]->GetN();
  //}
  length = thisAWF->even()->GetN();

  //if(pol==0){
  //  times = grEv[ant]->GetX();
  //  volts = grEv[ant]->GetY();   
  //}
  //else if(pol==1){
  //  times = grEvHoriz[ant]->GetX();
  //  volts = grEvHoriz[ant]->GetY();
  //}
  times = thisAWF->even()->GetX();
  volts = thisAWF->even()->GetY();
 
  double average_phase;
  double average_x;
  double average_y;
  int num_avg=0;

  //Get Fourier Transform

  deltaT=times[1]-times[0];
   
  newLength=(length/2)+1;
  deltaF=1/(deltaT*length); //GHz
  deltaF*=1e3; //MHz
  //for (int k=0; k<length; ++k) {printf(" %i: %f \n", k, volts[k]);}  
  FFTWComplex *theFFT=FFTtools::doFFT(length,volts);  // TODO get from the AnalysisWaveform
  
  //cout<<"Freq is "<<Freq[0]<<"\n";
  for(int i=0;i<newLength;i++) {
    if (i==0) frequencyArray[i]=0;
    if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;//make a freq Array. 
    
    
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      
      magFFT[i] = sqrt(pow(theFFT[i].re,2)+pow(theFFT[i].im,2));
      phase_single[i]=atan2(theFFT[i].im,theFFT[i].re);
      phase_single[i]=phase_single[i];//+2*TMath::Pi()*k_unwrap;
      
      phase_single_unshifted[i]=atan2(theFFT[i].im,theFFT[i].re);
      phase_new[i]=phase_single[i];
      
    }//>200 <1200
      
    
    else{
      phase_single[i]=0.;
      phase_single_unshifted[i]=0.;
      magFFT[i]=-1000;
    }
    
  }//i  
  
  for(int i=0;i<newLength;i++) {
    if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
      magFFT_dB[i]=10*log10(sqrt(magFFT[i]/double(1))/10.);//correct RMS
    }    
    else {
      magFFT_dB[i]=-1000;
      
    }
    
  }
  
    double mean_mag;

    double val1;
    double val2;
    double val_max=TMath::Pi()/2.;
    double val_min=TMath::Pi()/2.;
    
    vector<double> lower_bounds;
    vector<double> upper_bounds;

   
    int one_flag=0;
    int two_flag=0;
    
    for(int j=0;j<(int)Freq.size();j++){
      lower_bounds.clear();
      upper_bounds.clear();
      minFreq = Freq[j]-(bandWidth[j]);//start of notch region
      maxFreq = Freq[j]+(bandWidth[j]);//end of notch region
     
      //set vectors that contain where notched regions are
      for(int k=0;k<(int)cutFreqs.size();k++){
	if(cutFreqs[k] > minFreq && cutFreqs[k]<maxFreq){
	  lower_bounds.push_back(cutFreqs[k]-deltaF);//one bin to left of step (needed to make sure we dont include step in average)
	  upper_bounds.push_back(cutFreqs[k]+deltaF);//one bin to right of step (neede to make sure we dont include step in average)
	}
      }//k

      
      lower_bounds.push_back(1200);
      upper_bounds.push_back(1200);
      int k_lower=0;
      int pass_flag=0;
     
      //start process
      for(int i=0;i<newLength;i++){//lowest freq in range to start of notch

	val_max=TMath::Pi()/2.;
	val_min=TMath::Pi()/2.;

	if(frequencyArray[i]>=minFreq && frequencyArray[i]<=maxFreq){//if inside notched region
	  one_flag=0;
	  two_flag=0;
	  pass_flag=0;
	 
	  if(frequencyArray[i] > upper_bounds[k_lower]+1){//increment to next region
	      k_lower++;
	  }
	 
	  if(frequencyArray[i] < lower_bounds[k_lower]) {//can do all freqs up to bin before step
	    pass_flag=1;
	  }
	
	  if(pass_flag==1){//calculate average for method. Want to use up to 7 bins for average, but cannot include where jump occurs so samples is not constant
	    
	    average_x = theFFT[i].re+theFFT[i-1].re + theFFT[i+1].re;
	    average_y = theFFT[i].im+theFFT[i-1].im + theFFT[i+1].im;
	    mean_mag = magFFT[i-1]+magFFT[i]+magFFT[i+1];
	    one_flag=1;
	    num_avg=3;
	  }

	  pass_flag=0;//being careful about flags?
	 
	  if(one_flag==1){
	    pass_flag=1;
	    for(int k=0;k<(int)cutFreqs.size();k++){//can we expand average out one more bin on either side?
	      if(frequencyArray[i-2] > cutFreqs[k]-1 && frequencyArray[i-2] < cutFreqs[k]+1 ){
		pass_flag=0;
	      }
	      if(frequencyArray[i+2] > cutFreqs[k]-1 && frequencyArray[i+2] < cutFreqs[k]+1 ){
		pass_flag=0;
	      }
	    }//k
	    if(pass_flag==1){
	      //yes we can
	      average_x = average_x + theFFT[i-2].re + theFFT[i+2].re;
	      average_y = average_y + theFFT[i-2].im + theFFT[i+2].im;
	      mean_mag = mean_mag + magFFT[i-2] + magFFT[i+2];
	      num_avg+=2;
	      two_flag=1;
	    }
	  }
	  pass_flag=0;
	  
	  if(two_flag==1){
	    pass_flag=1;
	    for(int k=0;k<(int)cutFreqs.size();k++){//if we expanded to 5 bins for average, can we go to 7 bins?
	      if(frequencyArray[i-3] > cutFreqs[k]-1 && frequencyArray[i-3] < cutFreqs[k]+1){
		pass_flag=0;
	      }
	      if(frequencyArray[i+3] > cutFreqs[k]-1 && frequencyArray[i+3] < cutFreqs[k]+1){
		pass_flag=0;
	      }
	    }//k
	    if(pass_flag==1){
	      //yes
	      average_x = average_x + theFFT[i-3].re + theFFT[i+3].re;
	      average_y = average_y + theFFT[i-3].im + theFFT[i+3].im;
	      mean_mag = mean_mag + magFFT[i-3] + magFFT[i+3];
	      num_avg+=2;
	    }
	  }

	  //get average mag and phase
	  average_x = average_x/num_avg;
	  average_y = average_y/num_avg;
	  mean_mag = mean_mag/num_avg;
	  
	  average_phase = atan2(average_y,average_x);
	 
	  delta = abs(average_phase - phase_single[i]);
	  //do solution of phasor equation

	  //DO I NEED BOTH EQUATIONS WITH MY SIMPLIFICATION?
	  //printf("  i=%i, averagePhase=%f, phase_single=%f, delta=%f \n", i, average_phase, phase_single[i], delta);
	  val1 = solveGamma_plus(average_phase, phase_single[i], delta);
	  val2 = solveGamma_minus(average_phase, phase_single[i], delta);
	  if(val1 != val1 && val2 != val2){//something went wrong!
	    //cout<<"Problem! val1,val2 = "<<val1<<" " <<val2<<" mean, phase_single, delta are "<<average_phase<<" "<<phase_single[i]<<" "<<average_x<<" "<<average_y<<"\n";
	  }
	  /////////SINCE CHANGE OF FORMULA, DONT KNOW IF WE NEED THIS ANYMORE.
	   if(cos(val1-mean_mag) > cos(val2-mean_mag)){
	     val_max += val1;
	     val_min += val2;
	   }
	   else{
	     val_max+=val2;
	     val_min+=val1;
	   }
	  
	   //DOES VAL_MIN==VAL_MAX?
	   if(val1 == val1 && val2 == val2 && one_flag==1){//being careful about NaNs
	     if(magFFT[i]>=mean_mag){
	       phase_new[i]=val_max;
	       
	     }
	     if(magFFT[i]<mean_mag){
	       phase_new[i]=val_min;
	      
	     }
	   }
	   
	}//minFreq
	

      }//i
     
    }//j=Freq.size

    //change the Fourier components
    double x;
    double y;
    for(int i=0;i<newLength;i++){
      if (frequencyArray[i]>=200 && frequencyArray[i]<=1200){
	x = magFFT[i]*cos(phase_new[i]);
	y = magFFT[i]*sin(phase_new[i]);

	theFFT[i].re = x;
	theFFT[i].im = y;
	
      }
    }
    //FFT back
    double *filteredVals = FFTtools::doInvFFT(length,theFFT);
    TGraph *grTime = new TGraph(length,times,volts);
    double *times1 = grTime->GetX();
    for(int i=0;i<length;i++){
      volts[i] = filteredVals[i];
    }
    //if(pol==0){
    //  delete grEv[ant];
    //  grEv[ant] = new TGraph(length,times1,filteredVals);
    //}
    //if(pol==1){
    //  delete grEvHoriz[ant];
    //  grEvHoriz[ant] = new TGraph(length,times1,filteredVals);
    //}
    thisAWF->updateEven(grTime);
    delete [] theFFT;
    delete [] filteredVals;
    delete grTime;
   
  
}//geommethod

/////////------------------------------------------------------------------------------------------------------------------------

double GeometricFilter::solveGamma_plus(double theta, double psi, double delta) {
  double gamma;
  double sqrt_val;
  double sin_delta = sin(delta);
  double cos_delta = cos(delta);
  double arg = psi-theta;
  sqrt_val = 1-pow(cos_delta/cos(arg),2);
   sqrt_val = sqrt(sqrt_val);
  gamma = sin(2*arg)*(1+sqrt_val)/(2*sin_delta);
  gamma = acos(gamma);
  gamma = theta + gamma;
 
  return gamma; 
}

//----------------------------------------------------------------------------------------------------------------------------------

double GeometricFilter::solveGamma_minus(double theta, double psi, double delta) {
  double gamma;
  double sqrt_val;
  double sin_delta = sin(delta);
  double cos_delta = cos(delta);
  double arg = psi-theta;
  sqrt_val = 1-pow(cos_delta/cos(arg),2);
  sqrt_val = sqrt(sqrt_val);
  gamma = sin(2*arg)*(1-sqrt_val)/(2*sin_delta);
  gamma = acos(gamma);

  gamma = theta + gamma;

  return gamma; 
}




