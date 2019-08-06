#ifndef FREQ_DOMAIN_FUNCTION
#define FREQ_DOMAIN_FUNCTION

#include "AnalysisWaveform.h" 
#include "SystemResponse.h" 
#include <complex> 
#include "TF1.h" 



//note that there are two secret parameters (amplitude and time offset) before the rest of the paramters!!!

class FreqDomainFunction 
{


  public: 

    FreqDomainFunction (  std::complex<double> (*fn) (double f, const double * pars), 
        int npars, const AnitaResponse::AbstractResponse * r = 0,
        double dt = 0.1, double t0 = 0, double t1 = 102.4, bool shift = true);
    double eval(double x, const double * p = 0); 
    std::complex<double> evalFreq(double f, const double * p = 0); 
    double evalPower(double f, const double * p = 0); 
    double evalPhase(double f, const double * p = 0); 
    


    // 0 = no, 1 = use real part and hilbert transform of real part, 2 = use imag part and hilbert transform of imag part
    // 3 = average of 1 and 2 
    void forceCausal(int force) { causal_param = force; wf_init = false; }
    void dedisperseResponse(bool d) { dedisperse_response = d; } 
    void setDebug(bool d) { debug = d; } 
    void setParameters(const double *p); 
    void setResponse(const AnitaResponse::AbstractResponse * r) { response = r; wf_init = false; }

    double operator()(const double *x, const double *p = 0) { return eval(*x,p);  }

    TF1 * makeTF1(const char * name) { TF1 * f = new TF1(name, this, t0,t1, NPar()); f->SetNpx(wf.Neven()); f->SetParameters(parameters()); return f; } 
    const double * parameters() const { return &pars[0]; } 
    unsigned int NPar() const { return pars.size(); } 

    const AnalysisWaveform * waveform() const { return &wf; }

    //how? 0 = use real, 1 = use imag, 2 = average the two 
    static std::complex<double> * makeCausal(int N, const std::complex<double> * in, int how, std::complex<double> * out = 0); 

  private:
    bool do_shift; 
    bool wf_init; 
    bool debug; 
    double t0;
    double t1;
    std::complex<double> (*the_fn) (double f, const double * pars); 
    std::vector<double> pars; 
    const AnitaResponse::AbstractResponse * response;
    AnalysisWaveform wf; 
    int causal_param; 
    bool dedisperse_response; 

};





#endif



