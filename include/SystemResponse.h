#ifndef _UCORRELATOR_SYSTEM_RESPONSE_H 
#define _UCORRELATOR_SYSTEM_RESPONSE_H

  /* oh boy, down the OO rabbit hole we go! 
   *
   *  this file defines methods for deconvolution and system responses
   *
   **/

#include <vector> 
#include <complex>
#include "TH2.h" 
#include "TMutex.h" 
#include "FFTWComplex.h" 
#include <map>

class AnalysisWaveform; 
class TGraph; 
class TF1; 

namespace AnitaResponse{  

 
  class DeconvolutionMethod
  {
  
    public: 
      virtual void deconvolve(AnalysisWaveform *wf, const AnalysisWaveform * response) const = 0; 

      //compatibility interface
      __attribute__((deprecated)) 
      virtual void deconvolve(int Nf, double df, FFTWComplex *Y, const FFTWComplex * response) const; 

      virtual ~DeconvolutionMethod()  { ; } 
  }; 


  class NaiveDeconvolution : public DeconvolutionMethod
  {

    public: 

      virtual void deconvolve(AnalysisWaveform *wf, const AnalysisWaveform *rwf) const ; 

      virtual ~NaiveDeconvolution()  { ; } 

  }; 



  class CLEANDeconvolution : public DeconvolutionMethod 
  {

    public: 
      // restoring function accepts a delta_t as parameter. Will be scaled by restoring component 
      CLEANDeconvolution(const TF1 * restoring_beam, double gain = 0.05, double threshold = 0.1, bool convolve_residuals = true) 
        : r(restoring_beam), gain(gain), threshold(threshold), convolve_residuals(convolve_residuals) {; } 


      virtual void deconvolve(AnalysisWaveform *wf, const AnalysisWaveform *rwf) const  
      {
        deconvolveSavingIntermediate(wf,rwf,0,0,0); 
      }


      virtual void deconvolveSavingIntermediate(AnalysisWaveform *wf, 
                              const AnalysisWaveform * response, 
                              std::vector<double> * components = 0,
                              std::vector<AnalysisWaveform*> *xcorr = 0,
                              std::vector<AnalysisWaveform*> *ys = 0
                              ) const ; 



      static TGraph * makeComponentGraph(std::vector<double> * components) { return 0 ; }





    private: 
      const TF1 * r; 
      double gain; 
      double threshold; 
      bool convolve_residuals;

  }; 

  /*
  class RichardsonLucyDeconvolution : public DeconvolutionMethod
  {
    public: 
      RichardsonLucyDeconvolution(int maxiter = 100, double thresh = 0.01) 
        : maxiter(maxiter), thresh(thresh) {; }
      virtual void deconvolve(size_t N, double df, FFTWComplex * Y, 
                              const FFTWComplex * response) const ; 

    private: 
      int maxiter; 
      double thresh; 

  }; 
  */

  class MinimumPhaseDeconvolution
  {
    public: 
      virtual void deconvolve(AnalysisWaveform *wf, const AnalysisWaveform *rwf) const ; 
      virtual ~MinimumPhaseDeconvolution() { ; } 


  }; 

  /** Wiener Deconvolution where the SNR is given as a TGraph or a TF1 
   *
   * */ 

  class WienerDeconvolution : public DeconvolutionMethod
  {

    public: 
      /** Construct a Wiener Deconvolution where the SNR is given as a TGraph. 
       * This TGraph may be modified externally, but must exist for the lifetime of this object.
       * as the deconvolution does not own it. 
       *
       * Optionally, a pointer to a scale value may be given, which may also be modified externally, and must also exist for the lifetime of this object. */
      WienerDeconvolution (const TGraph * g_snr, const double * scale = 0); 

      /** Construct a Wiener Deconvolution ``automatically'' using the system response, where noise_level is the effective digitizer noise level relative to the response*/ 
      WienerDeconvolution(double noise_level = 1); 

      /** Construct a WienerDeconvolution where the SNR is given as a TF1. This function may be modified externally but must exist for the lifetime of this object. */
      WienerDeconvolution (const TF1 * f); 


      virtual void deconvolve(AnalysisWaveform *wf, const AnalysisWaveform *rwf) const ; 

    protected: 
      virtual double snr(double f, double R2, int N) const; 

    private: 
      const TGraph * snr_graph;
      const double * scale; 
      const TF1 * snr_function; 
      double min,max;
      double noise_level; 
  }; 

  class BandLimitedDeconvolution : public DeconvolutionMethod
  {
    public: 

      BandLimitedDeconvolution(double minfreq, double maxfreq, int edgeorder = 0) 
        : min_freq(minfreq),max_freq(maxfreq), edge_order(edgeorder) 
      {
      }

      virtual void deconvolve(AnalysisWaveform *wf, const AnalysisWaveform *rwf) const ; 

      virtual ~BandLimitedDeconvolution() { ; } 


    private: 
      double min_freq, max_freq;
      int edge_order;
  }; 


  class AllPassDeconvolution : public DeconvolutionMethod
  {
    public: 
      AllPassDeconvolution() { ; } 
      virtual void deconvolve(AnalysisWaveform *wf, const AnalysisWaveform *rwf) const ; 

      virtual ~AllPassDeconvolution() { ; } 

  }; 

  class ImpulseResponseXCorr : public DeconvolutionMethod
  {

    public: 
      ImpulseResponseXCorr() { ; } 
      virtual void deconvolve(AnalysisWaveform *wf, const AnalysisWaveform *rwf) const ; 
      virtual ~ImpulseResponseXCorr() { ; } 

  }; 

  extern DeconvolutionMethod & kDefaultDeconvolution; 

  class AbstractResponse
  {

   public: 
      virtual FFTWComplex getResponse(double f, double angle = 0) const = 0; 
      virtual FFTWComplex * getResponseArray(int N, const double  * f, double angle = 0) const; 
      virtual FFTWComplex * getResponseArray(int N, double df, double angle = 0) const ; 
      virtual double getMagnitude(double f, double angle= 0) const;  
      virtual double getPhase(double f, double angle = 0) const; 

      virtual AnalysisWaveform * impulseResponse(double dt = 1./2.6, int N = 256) const; 
      virtual AnalysisWaveform * convolve(const AnalysisWaveform * wf, double angle = 0) const; 
      virtual AnalysisWaveform * deconvolve(const AnalysisWaveform * wf, const DeconvolutionMethod * method = &kDefaultDeconvolution, double angle = 0) const; 
      virtual void convolveInPlace(AnalysisWaveform * wf, double angle = 0) const; 
      virtual void deconvolveInPlace(AnalysisWaveform * wf, const DeconvolutionMethod * method = &kDefaultDeconvolution, double angle = 0) const; 

      virtual ~AbstractResponse() { ; } 

  };

  /** This class is a bit over-engineered right now since it supports
   *  responses at different angles */ 
  class Response : public AbstractResponse
  {
    public: 
      Response(int NFreq, double df); 
      Response(const TGraph * time_domain, int npad); 
      Response(int Nfreq, double df, int nangles, const double * angles, const FFTWComplex ** responses);  
      Response(int Nfreq, double df, const FFTWComplex * response);  

      void addResponseAtAngle(double angle, const FFTWComplex * response); 

      virtual FFTWComplex getResponse(double f, double angle = 0) const; 
       
      
      const TH2 * getReal() const { return &real; } 
      const TH2 * getImag() const { return &imag; } 

      virtual ~Response() { ; } 

    protected: 
      mutable TMutex lock; 
      int Nfreq; 
      double df; 
      int nangles; 
      std::map<double, FFTWComplex *> responses; 
      mutable TH2D real; 
      mutable TH2D imag; 
      mutable bool dirty; 
      void recompute() const; 
  }; 


  class CompositeResponse : public AbstractResponse
  {
    public:  
      void addResponse(const AbstractResponse * response) { responses.push_back(response); } 
      virtual FFTWComplex getResponse(double f, double angle = 0) const; 
      virtual ~CompositeResponse() { ; } 
   
    private: 
      std::vector<const AbstractResponse * > responses; 
  }; 


}



#endif 

