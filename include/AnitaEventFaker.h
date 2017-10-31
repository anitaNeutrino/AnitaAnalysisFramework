#ifndef ANITA_EVENT_FAKER_HH
#define ANITA_EVENT_FAKER_HH

/* The event faker adds an artifical signal 
 * on top of an existing event. The signal is the minimum phase
 * waveform with the associated bandwidth; 
 *
 * The intended use is to study pointing resolution with thermal noise events 
 * but can also be used for stokes parameter estimation, etc. 
 *
 * 
 **/ 

class UsefulAnitaEvent;

#include <complex> 
#include "ResponseManager.h" 
#include "TF1.h" 
#include "TF2.h" 
#include "AnalysisWaveform.h" 

class AnitaEventFaker 
{
  public: 
    AnitaEventFaker(const char * responseDir, const TF1 & mag_response = getDefaultMagResponse(), 
                    double delay = 30, double signal_dt = 5e-2); 

    /** 
     * These generates the signal based on the minimum phase signal
     * corresponding to the specified magnitude response. The signal will be convolved with the impulse response. 
     *
     * @param mag_response a function giving the magnitude response (in relative units, as the power will be scaled away )
     * @param delay  the delay to apply to the signal (should be about the trigger delay) 
     * @param signal_dt  the level to interpoalte the signal too. it will then be linearly interpolated on top of the event, so this should be significantly smaller than the sampling period of the events 
     * */ 
    void setSignalFromMagnitudeResponse(const TF1 & mag_response, double delay = 30, double signal_dt = 5e-2, int npoints = 2000); 
    void setSignalFromMagnitudeResponse(const TGraph& mag_response, double delay = 30, double signal_dt = 5e-2, int npoints = 2000); 


    /** This explicitly sets the signal
     * Note that the power will be rescaled to one. 
     * */ 
    void setSignal(const AnalysisWaveform& sig)  ; 




    /** Add a signal to this event coming from that direction
     * @param event the event to add the signal to 
     * @param theta_deg the degrees theta (with positive coming from below) for the signal , in payload coordinates
     * @param phi_deg   the degrees phi of the signal , in payload coordinates
     * @param A         a scale factor for the signal
     * @param rel_Q     double lin_pol_angle 
     *
     */ 
    void addSignal(UsefulAnitaEvent * event, double theta_deg, double phi_deg, 
                   double A = 10,  std::complex<double> jones_H = std::complex<double>( 1,0),  std::complex<double> jones_V = std::complex<double>(0,0) ) const; 

    /** This is just a helper function to make a pure noise event, where each channel is populated with random noise 
     *
     * @param rms  The RMS of the noise, in pre-response units 
     * @param victim If non-zero, this event will be filled with noise instead of constructing a new event. 
     *
     * */  


    UsefulAnitaEvent * makePureNoiseEvent(double rms = 0.1, UsefulAnitaEvent * victim = 0) const; 



    /** Return the current realization of the signal. This is the prottype convolved with each response */ 
    const AnalysisWaveform & getSignal(int ant, int pol) const { return signal[ant][pol]; } 

    const AnalysisWaveform & getPrototype() const { return prototype; } 

    const TF2 & getHPolAntennaGainParameterization() const { return offAxisGain_hpol; } 
    const TF2 & getVPolAntennaGainParameterization() const { return offAxisGain_vpol; } 
    const TF2 & getOffAxisDelay() const { return offAxisDelay; } 

    void setHPolAntennaGainParameterization(const TF2 & f) { offAxisGain_hpol = f; } 
    void setVPolAntennaGainParameterization(const TF2 & f) { offAxisGain_vpol = f; } 
    void setOffAxisDelay(const TF2 & f) { offAxisDelay = f; } 

    /** Retrieves the default magnitude response used. */
    static const TF1 & getDefaultMagResponse(); 

    static const TF2 & getDefaultHPolAntennaGain(); 
    static const TF2 & getDefaultVPolAntennaGain(); 
    static const TF2 & getDefaultOffAxisDelay(); 

  private: 

    AnitaResponse::ResponseManager manager; 
    AnalysisWaveform prototype; 
    AnalysisWaveform signal[NUM_SEAVEYS][2]; 
    TF2 offAxisGain_hpol; 
    TF2 offAxisGain_vpol; 
    TF2 offAxisDelay; 
    


};


#endif
