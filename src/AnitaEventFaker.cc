#include "AnitaEventFaker.h" 
#include "UsefulAnitaEvent.h" 
#include "SystemResponse.h" 
#include "AnitaGeomTool.h" 
#include "FFTtools.h" 
#include "TRandom3.h" 
#include <fftw3.h> 


const TF2 & AnitaEventFaker::getDefaultHPolAntennaGain() 
{
  static TF2 f ("default_hpol_antenna_gain", "exp(-x*x/(2*(20*20))) * exp(-y*y/(2*(20*20)))", -90,90,-90,90); 
  return f;
}


const TF2 & AnitaEventFaker::getDefaultVPolAntennaGain() 
{
  static TF2 f ("default_vpol_antenna_gain", "exp(-x*x/(2*(20*20))) * exp(-y*y/(2*(20*20)))", -90,90,-90,90); 
  return f; 
}

const TF2 & AnitaEventFaker::getDefaultOffAxisDelay() 
{
  const double c2 = 1.45676e-8; 
  const double c1 = 5.01452e-6; 
  static TF2 f("default_off_axis_delay", "[0] * (x*x+y*y) + [1] * (x*x*x*x + y*y*y*y)", -90,90,-90,90); 
  f.SetParameter(0,c1); 
  f.SetParameter(1,c2); 
  return f; 
}

const TF1 & AnitaEventFaker::getDefaultMagResponse() 
{
  static TF1 f("default_mag_response", "1",0,2); 
  return f; 
}

AnitaEventFaker::AnitaEventFaker(const char * responseDir, const TF1 & mag_response , double delay , double signal_dt) 
      : manager( responseDir, ceil((1./2.6)/ signal_dt)), 
        offAxisGain_hpol(getDefaultHPolAntennaGain()), 
        offAxisGain_vpol(getDefaultVPolAntennaGain()), 
        offAxisDelay(getDefaultOffAxisDelay())
    {
      setSignalFromMagnitudeResponse(mag_response, delay, signal_dt); 
    }

static void normalizeSignal(AnalysisWaveform * sig) 
{
  double scale = sig->even()->getSumV2(); 
  TGraphAligned * g = sig->updateEven(); 
  for (int i = 0; i < g->GetN(); i++)
  {
    g->GetY()[i] /= scale; 
  }

}

void AnitaEventFaker::setSignal(const AnalysisWaveform & sig) 
{
  prototype = sig; 
  normalizeSignal(&prototype); 
  prototype.padEven(0); 

  for (int ant = 0; ant < NUM_SEAVEYS; ant++)
  {
    for (int ipol = 0; ipol < 2; ipol++) 
    {
      signal[ant][ipol] = prototype; 
      manager.response(ipol,ant)->convolveInPlace(&signal[ant][ipol]); 
    }
  }
}


template <class Evalable>  
static AnalysisWaveform * makeSignalFromMagnitudeResponse( const Evalable  & mag_response, double delay, double signal_dt, int N ) 
{
  //find the minimum df in the magnitude response
  double sig_df = 1 / (N*signal_dt); 
  std::vector <double> G(N/2+1); 

  for (int i = 0; i < N/2+1; i++) 
  {
    G[i]=mag_response.Eval(sig_df * i); 
  }

  FFTWComplex * F = FFTtools::makeMinimumPhase(N/2+1, &G[0]); 

  AnalysisWaveform * answer =  new AnalysisWaveform (N, F, sig_df, delay); 

  delete [] F; 

  return answer; 
}


void AnitaEventFaker::setSignalFromMagnitudeResponse(const TF1 & mag_response, double delay , double signal_dt, int npoints ) 
{
  AnalysisWaveform * tmp = makeSignalFromMagnitudeResponse( mag_response, delay, signal_dt, npoints); 
  setSignal(*tmp); 
  delete tmp; 

}

void AnitaEventFaker::setSignalFromMagnitudeResponse(const TGraph& mag_response, double delay , double signal_dt, int npoints ) 
{
  AnalysisWaveform * tmp = makeSignalFromMagnitudeResponse( mag_response, delay, signal_dt, npoints); 
  setSignal(*tmp); 
  delete tmp; 
}



UsefulAnitaEvent * AnitaEventFaker::makePureNoiseEvent(double rms, UsefulAnitaEvent * victim) const
{
//  fprintf(stderr,"WARNING: makePureNoiseEvent() probably doesn't work yet\n"); 
  UsefulAnitaEvent * u = victim; 

  //if we don't have a model event, generate a fake time base for this one 
  if (!victim)
  {
    u = new UsefulAnitaEvent; 
    for (int i = 0; i < NUM_DIGITZED_CHANNELS; i++) 
    {
      double t = 0; 
      u->fNumPoints[i] = NUM_SAMP; 
      for (int j = 0; j < NUM_SAMP; j++) 
      {
        u->fTimes[i][j] = t; 
        double dt = 0; 
        while (dt <= 0) dt = gRandom->Gaus(1/2.6,0.09); //TODO make more general  
        t+=dt; 
      }
    }
    memset(u->fVolts,0,sizeof(u->fVolts)); 
  }

  /* Next, we will generate thermal noise according to the instrument response, upsample it, then evaluate it at the times */ 

  for (int ch = 0; ch < NUM_DIGITZED_CHANNELS; ch++) 
  {
    int chan, surf; 
    AnitaGeomTool::getSurfChanFromChanIndex(ch, surf, chan); 
    if (chan == 8) continue; //this is the clock, ignore it; 
    int ant; 
    AnitaPol::AnitaPol_t pol ; 

    AnitaGeomTool::getAntPolFromSurfChan(surf,chan, ant, pol); 

    //empty waveform, oversampled and a bit longer
    AnalysisWaveform wf(360, 1/2.6, -20); //TODO, these should be configurable! 

    int N = wf.Neven(); 
    TGraphAligned * g = wf.updateEven(); 
    for (int i = 0; i < N; i++) g->GetY()[i] = gRandom->Gaus(0, rms); 

    double df = wf.deltaF(); 
    FFTWComplex * G = wf.updateFreq(); 
    
    for (int i = 0; i < N/2+1; i++) 
    {
      G[i] *= manager.response(pol, ant)->getResponse(i * df); 
    }

    //upsample a bunch
    wf.padFreq(4); 

    //evaluate at the times
    wf.evalEven(u->fNumPoints[ch], u->fTimes[ch], u->fVolts[ch]); 
  }

  return u; 
}



void AnitaEventFaker::addSignal(UsefulAnitaEvent * event, double theta, double phi, double A, 
                                std::complex<double> jones_H, std::complex<double> jones_V) const
{
  

  double th_rad = theta * TMath::DegToRad(); 
  double phi_rad = phi * TMath::DegToRad(); 


  for (int ich = 0; ich < NUM_DIGITZED_CHANNELS; ich++)
  {
    int chan, surf; 
    AnitaGeomTool::getSurfChanFromChanIndex(ich, surf, chan); 
    if (chan == 8) continue; //this is the clock, ignore it; 
    int ant; 
    AnitaPol::AnitaPol_t pol ; 
    AnitaGeomTool::getAntPolFromSurfChan(surf,chan, ant, pol); 
    double sig_t0 = signal[ant][pol].even()->GetX()[0]; 
    double sig_t1 = signal[ant][pol].even()->GetX()[signal[ant][pol].Neven()-1]; 

    //we need to figure out the gain and delay for this channel 
    
    double R = AnitaGeomTool::Instance()->getAntR(ant, pol); 
    double z = AnitaGeomTool::Instance()->getAntZ(ant, pol); 
    double phi0_rad =  AnitaGeomTool::Instance()->getAntPhiPositionRelToAftFore(ant,pol); 
    double phi0=  phi0_rad  * TMath::RadToDeg(); 
    double Greal = A * (pol == AnitaPol::kHorizontal ?  offAxisGain_hpol :  offAxisGain_vpol).Eval(phi-phi0, theta-10); //TODO: don't hardcode 10 
    std::complex<double> G = Greal * (pol == AnitaPol::kHorizontal ? jones_H : jones_V); 

    double ts= (z * tan(th_rad) - R * cos(phi_rad-phi0_rad)) *  1e9 * cos(th_rad) / C_LIGHT; 
    ts+= offAxisDelay.Eval( phi-phi0, theta-10); 

    for (int i = 0; i < event->fNumPoints[ich]; i++) 
    {

      if (event->fTimes[ich][i] - ts< sig_t0) continue; 
      if (event->fTimes[ich][i] - ts> sig_t1) continue; 

      std::complex<double> val = G * std::complex<double> ( signal[ant][pol].evalEven(event->fTimes[ich][i] - ts),
                                                            signal[ant][pol].hilbertTransform()->evalEven(event->fTimes[ich][i] - ts)); 
      event->fVolts[ich][i] += std::real(val); 
    }
  }
}





