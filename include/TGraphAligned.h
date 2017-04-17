#ifndef _TGRAPH_ALIGNED_H
#define _TGRAPH_ALIGNED_H

#include "TGraph.h" 

/** 
 *  Slight modification of TGraph where the arrays are guaranteed to be aligned properly. The wisdom
 *  of this class is questionable; in principle it might help the compiler autovectorize, but that seems like 
 *  a fool's errand. For manual vectorization, it can make loads faster in certain cases. 
 *
 */ 




class TGraphAligned : public TGraph {


  /** This is the alignment */ 
#define TGRAPH_ALIGNED_ALIGNMENT 32 
  typedef Double_t * aligned_double_v __attribute__((aligned (TGRAPH_ALIGNED_ALIGNMENT))); 
  //Have to reimplement all interesting constructors, unfortunately, since we
  //unfortunately cannot call the TGraph constructor and have it call our
  //implementation of CtorAllocate 

  public: 
    /** Empty TGraphAligned */ 
    TGraphAligned(); 

    /** Zero'd TGraphAligned of size n*/ 
    TGraphAligned(Int_t n); 
    TGraphAligned & operator=(const TGraphAligned&);
    /** Create from arrays, which are copied */ 
    TGraphAligned(Int_t n, const Double_t * x, const Double_t * y); 

    /** Create a TGraphAligned from a TGraph. Use this in case I haven't yet implemented the TGraph constructor you want; it's just an extra copy. 
     *  I think this should fire also for TGraphAligned, but not sure */ 
    TGraphAligned(const TGraph  &); 


    double * getMoments(int i, double origin, double * moments = 0) const; 

    void dBize(double mindB=-200); 
    void undBize(); 


    /* just copies the numbers from g*/ 
    void adopt(const TGraphAligned * g); 


    void getMeanAndRMS(Double_t * mean, Double_t *  rms, Int_t istart = 0, Int_t iend = -1) const; 
    Double_t getSumV2(Int_t istart = 0, Int_t iend = -1) const; 

    Double_t peakVal(Int_t * location = 0, Int_t istart= 0, Int_t iend = -1, bool abs = false) const; 
    Double_t pk2pk(Int_t nthMax = 0, Int_t nthMin = 0,Int_t * location_max = 0, Int_t * location_min = 0, Int_t istart= 0, Int_t iend = -1) const; 


    aligned_double_v GetX() const { return fX; } 
    aligned_double_v GetY() const { return fY; } 

    /** Destructor */ 
    virtual ~TGraphAligned();

  protected: 

    virtual Double_t **AllocateAlignedArrays(Int_t Narrays, Int_t arraySize); 
    inline Double_t **Allocate(Int_t newsize) { return AllocateAlignedArrays(2, newsize); }
    virtual Bool_t CtorAllocate (void); 

    //needs to be reimplement to avoid delete[] on memalign-allocated pointer, which Valgrind will perhaps unjustly complain about
    virtual void CopyAndRelease(Double_t ** newarrays, Int_t ibegin, Int_t iend, Int_t obegin); 
    
    ClassDef(TGraphAligned,1); 

}; 


#endif 




